#!/usr/bin/env python3
"""Measure one real BOMEX loss-gradient configuration in a fresh CPU process.

Run from the repository root, for example:
    python3 clubb_jax/auto_llm_generated_tests/integration_tests/benchmark_gradient_checkpointing.py -output_dir output/checkpoints
    python3 clubb_jax/auto_llm_generated_tests/integration_tests/benchmark_gradient_checkpointing.py -timesteps 120 -gradient_mode forward -output_dir output/checkpoints

main() initializes the normal persistent clubb_loss_runner, then evaluates all
four metrics for wp2/cloud_frac and case-loss derivatives with respect to C8, beta and
Lscale_tke_multiplier. The first evaluation includes tracing and compilation or
loading from a configured compilation cache; later calls reuse the same session
and inputs. Outputs are JSON metadata and raw NPZ arrays
named by target, mode, timestep count and column count. Compare those arrays
across configurations to check values/gradients; finiteness alone is not a
correctness assertion. This is a manual performance workflow, not a pytest.
memory_sample() reads the Linux resident/swap counters used by the sampler.

The managed JAX launcher owns runtime setup and CPU compiler flags. BOMEX uses
its native vertical grid and dt=60 s, with one average over the second half of
the trajectory. Timesteps must be positive and even. Reverse mode checkpoints
each numerical timestep; forward mode uses ordinary advancement. Linux process
RSS is sampled every 50 ms during each call; JSON records sampled peaks, process-lifetime high-water RSS and swapped
memory. Peaks include compilation, executable caches and allocator retention,
not just autodiff intermediates. Pullback array leaves are also counted as
logical saved-residual bytes, deduplicated by Python array identity; this is a
storage diagnostic, not a measurement of physical buffer allocations.
The temporary Python containers used to count those arrays also contribute to
the measured process peak.
Use a separate invocation per configuration;
multiple configurations in one process would contaminate memory comparisons.
"""
from pathlib import Path
import argparse
import json
import hashlib
import os
import platform
import resource
import sys
import subprocess
import threading
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))


def memory_sample():
    """Return this Linux process's resident and swapped memory in bytes."""
    status = dict(line.split(':', 1) for line in Path('/proc/self/status').read_text().splitlines())
    return {name: int(status[name].split()[0]) * 1024 for name in ('VmRSS', 'VmSwap')}


def main():
    """Run one configuration, assert finite data and retain measured/raw results."""
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument('-help', '-h', action='help')
    parser.add_argument('-timesteps', type=int, default=20)
    parser.add_argument('-gradient_target', choices=('loss', 'metrics'), default='loss')
    parser.add_argument('-gradient_mode', choices=('forward', 'reverse'), default='reverse')
    parser.add_argument('-batch_size', type=int, default=1)
    parser.add_argument('-repetitions', type=int, default=2)
    parser.add_argument('-output_dir', type=Path, required=True)
    args = parser.parse_args()
    if args.timesteps <= 0 or args.timesteps % 2 or args.batch_size < 1 or args.repetitions < 1:
        parser.error('Use a positive even timestep count and positive batch_size/repetitions')
    from clubb_jax.run_jax import ensure_environment
    ensure_environment()
    # Configure the CPU compiler policy before the first JAX backend lookup.
    from clubb_jax import src
    import jax
    import jaxlib
    import numpy as np
    from tuner.case_defaults import read_case_defaults
    from utilities.clubb_loss_runner import clubb_loss_runner
    if '--combiner-disabled=true' not in os.environ.get('XLA_FLAGS', ''):
        raise RuntimeError('The documented CPU compiler workaround must be configured before JAX initialization')
    if jax.default_backend() != 'cpu':
        raise RuntimeError('This RSS benchmark requires the CPU backend')

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    label = f'{args.gradient_target}_{args.gradient_mode}_n{args.timesteps}_c{args.batch_size}'
    settings = read_case_defaults('bomex', {'time_average_range': [args.timesteps * 30, args.timesteps * 60],
                                          'num_time_windows': 1})
    record = {'gradient_target': args.gradient_target, 'case': 'bomex', 'timesteps': args.timesteps, 'checkpointed': args.gradient_mode == 'reverse',
              'gradient_mode': args.gradient_mode, 'columns': args.batch_size,
              'fields': ['wp2', 'cloud_frac'], 'parameters': ['C8', 'beta', 'Lscale_tke_multiplier'],
              'jax': jax.__version__, 'jaxlib': jaxlib.__version__, 'python': platform.python_version(),
              'platform': platform.platform(), 'xla_flags': os.environ.get('XLA_FLAGS', ''),
              'settings': settings, 'runs': [], 'status': 'running',
              'load_start': os.getloadavg(),
              'compilation_cache_enabled': str(jax.config.jax_enable_compilation_cache),
              'compilation_cache_dir': str(jax.config.jax_compilation_cache_dir),
              'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
              'source_sha256': {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                  for name in ('clubb_jax/src/advance_clubb_to_end.py', 'clubb_jax/src/clubb_loss_driver.py',
                               'utilities/clubb_loss_runner.py', 'clubb_jax/src/CLUBB_core/stats_netcdf.py')}}
    path = out / f'{label}.json'
    with clubb_loss_runner('jax', 'bomex', record['fields'], out / f'{label}_session',
                           batch_size=args.batch_size, case_defaults=settings) as runner:
        state = runner.api.active_request.state
        record['dt_seconds'] = state['dt_main']
        record['levels'] = {'zt': state['nzt'], 'zm': state['nzm']}
        record['initialized_array_bytes'] = sum(leaf.size * leaf.dtype.itemsize
            for name, value in state.items() if name != '_initial_state'
            for leaf in jax.tree.leaves(value) if isinstance(leaf, jax.Array))
        params = runner.default_params.copy()
        for name, values in {'C8': np.linspace(.8, 1.2, args.batch_size) if args.batch_size > 1 else 1.,
                             'beta': 1.5, 'Lscale_tke_multiplier': 2.125}.items():
            params[:, runner.param_names.index(name)] = values
        record['initial_memory'] = memory_sample()
        path.write_text(json.dumps(record, indent=2) + '\n')
        original_vjp = jax.vjp
        residual_measurements = []

        def measured_vjp(*primals, **keywords):
            """Count arrays saved by the normal loss pullback without retaining them."""
            result = original_vjp(*primals, **keywords)
            if getattr(primals[0], '__module__', '') != 'clubb_jax.src.clubb_loss_driver':
                return result
            arrays = [leaf for leaf in jax.tree.leaves(result[1]) if isinstance(leaf, jax.Array)]
            unique = {id(leaf): leaf for leaf in arrays}
            residual_measurements.append({
                'array_leaves': len(arrays), 'unique_array_objects': len(unique),
                'logical_bytes': sum(leaf.size * leaf.dtype.itemsize for leaf in arrays),
                'unique_object_bytes': sum(leaf.size * leaf.dtype.itemsize for leaf in unique.values()),
            })
            return result

        for repetition in range(args.repetitions):
            peak = memory_sample()
            stop = threading.Event()

            def sample_memory():
                """Sample RSS while the synchronized loss calculation is active."""
                while not stop.wait(.05):
                    for key, value in memory_sample().items():
                        peak[key] = max(peak[key], value)

            sampler = threading.Thread(target=sample_memory, daemon=True)
            sampler.start()
            print('START', label, repetition, flush=True)
            started = time.perf_counter()
            cpu_started = time.process_time()
            try:
                jax.vjp = measured_vjp
                result = runner.evaluate(params, grad=True, gradient_parameters=record['parameters'],
                    gradient_mode=args.gradient_mode, gradient_target=args.gradient_target, profiles=True)
            finally:
                jax.vjp = original_vjp
                elapsed = time.perf_counter() - started
                cpu_elapsed = time.process_time() - cpu_started
                stop.set()
                sampler.join()
                for key, value in memory_sample().items():
                    peak[key] = max(peak[key], value)
            arrays = {f'{group}_{name}': value for group in ('metrics', 'gradient')
                      for name, value in result[group].items()}
            arrays['loss'] = result['loss']
            arrays.update({f'profile_{name}': profile['model'] for name, profile in result['profiles'].items()})
            from utilities.loss_metrics import loss_metric_values_valid
            for name, value in arrays.items():
                if name.startswith(('metrics_', 'gradient_')):
                    assert loss_metric_values_valid(name.split('_', 1)[1], value).all(), name
                else:
                    assert np.isfinite(value).all(), name
            np.savez_compressed(out / f'{label}_{repetition}.npz', **arrays)
            record['runs'].append({'phase': 'cold' if repetition == 0 else 'warm',
                'wall_seconds': elapsed, 'cpu_seconds': cpu_elapsed, 'sampled_peak_bytes': peak,
                'lifetime_peak_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
                'end_memory': memory_sample(), 'timestep_specializations': len(runner.api.active_request.checkpoint_cache),
                'saved_residuals': residual_measurements[-1] if residual_measurements else None})
            path.write_text(json.dumps(record, indent=2) + '\n')
            print('DONE', label, repetition, json.dumps(record['runs'][-1]), flush=True)
        record['load_end'] = os.getloadavg()
        record['status'] = 'complete'
        path.write_text(json.dumps(record, indent=2) + '\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
