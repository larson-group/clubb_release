#!/usr/bin/env python3
"""Measure mixing-length forward calls and parameter gradients on CPU or GPU.

Run from the repository root; main() prepares the managed JAX environment:

    python3 clubb_jax/tests/benchmark_mixing_length.py
    python3 clubb_jax/tests/benchmark_mixing_length.py -methods 3 \
        -columns 1 8 32 -levels 88 198 -samples 20 -output_file timings.json

Inputs and defaults:
  Time methods 1/2/3, 1/4/32 columns, and 88/198 thermodynamic levels, using
  30 samples per measurement. Columns have a 40 m grid, thlm gradient
  0.002 K/m, rt=0.008 kg/kg, TKE=0.5 m^2/s^2, and entrainment=0.001 1/m.
  PDF launch fraction is 0.35 with controlled component departures. The
  stratified profiles are artificial; they are not saved model trajectories.
  Debug is -1 and statistics are disabled. The launcher environment selects
  the device; e.g. CLUBB_JAX_ACCELERATOR=cuda13 with CUDA_VISIBLE_DEVICES=0.

Outputs and interpretation:
  Compile full calc_Lscale forward and value-and-parameter-gradient kernels
  separately. The scalar gradient loss is mean combined mixing length; its
  timing includes both loss evaluation and reverse differentiation. After
  three warmups, synchronize every call and report median milliseconds plus
  compiled temporary-memory bytes. Emit one JSON row per measurement and,
  with -output_file, a combined JSON file including version/device/layout.
  The file's parent directory must exist. Kernel layout constants are read
  from mixing_length.py; there are no size override options.

  main() is the script entry point. Nonfinite results fail an assertion;
  finite results alone establish neither derivative accuracy nor agreement
  with Fortran. These timings exclude other model physics and file output.
"""
import argparse
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def main():
    """Parse benchmark shapes/methods, prepare JAX, time each compiled kernel, and
    print/write JSON results.
    """
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument('-help', '-h', action='help')
    parser.add_argument('-columns', type=int, nargs='+', default=[1, 4, 32],
                        help='Column counts (default: 1 4 32).')
    parser.add_argument('-levels', type=int, nargs='+', default=[88, 198],
                        help='Thermodynamic level counts (default: 88 198).')
    parser.add_argument('-methods', type=int, nargs='+', choices=(1, 2, 3), default=[1, 2, 3],
                        help='Mixing-length methods to time (default: all three).')
    parser.add_argument('-samples', type=int, default=30,
                        help='Synchronized timed calls per measurement (default: 30).')
    parser.add_argument('-output_file', type=Path, help='Optional JSON results file.')
    args = parser.parse_args()
    if min(*args.columns, args.samples) < 1 or min(args.levels) < 3:
        parser.error('columns/samples must be positive and levels must be at least 3')
    from clubb_jax.run_jax import ensure_environment
    ensure_environment()
    import jax
    import jax.numpy as jnp
    import numpy as np
    from clubb_jax.src.CLUBB_core.constants_clubb import ep1
    from clubb_jax.src.CLUBB_core.err_info import ErrInfo
    from clubb_jax.src.CLUBB_core.error_code import set_debug_level
    from clubb_jax.src.CLUBB_core.grid_class import setup_grid
    from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
    from clubb_jax.src.CLUBB_core import mixing_length
    from clubb_jax.src.CLUBB_core.mixing_length import calc_Lscale
    from clubb_jax.src.CLUBB_core.parameters_tunable import init_clubb_params, PNAME_IDX
    from clubb_jax.src.CLUBB_core.pdf_params import init_pdf_params

    set_debug_level(-1)
    backend = 'CPU' if jax.default_backend() == 'cpu' else 'GPU'
    layout = {name: getattr(mixing_length, f'COHERENT_{backend}_{name.upper()}')
              for name in ('launch_group_size', 'group_min_columns', 'offset_block_size')}
    rows = []
    for ncol in args.columns:
        for nlev in args.levels:
            gr = setup_grid(ngrdcol=ncol, deltaz=40., zm_init=0., zm_top=40.*nlev)
            nzt, nzm = gr.nzt, gr.nzm
            thlm = 298. + .002*gr.zt
            rtm = jnp.full((ncol, nzt), .008)
            exner = jnp.exp(-3e-5*gr.zt)
            cp = jnp.asarray(init_clubb_params(ncol, ''))
            cp = cp.at[:, PNAME_IDX['Lscale_pdf_launch_frac']].set(.35)
            x = dict(thlm=thlm, rtm=rtm, thvm=thlm*(1.+ep1*rtm),
                     em=jnp.full((ncol, nzm), .5), params=cp,
                     pressure=1e5*exner**3.5, exner=exner)
            stats = JaxStats.empty(l_sample=False, names=(), ncol=ncol, max_nlev=nzm)
            for method in args.methods:
                def calculate(x):
                    """Return full calc_Lscale outputs for dynamic state/parameters at
                    the selected fixed shape and method.
                    """
                    pdf = init_pdf_params(nzt, ncol).replace(
                        thl_1=x['thlm']+.7, thl_2=x['thlm']-.4,
                        rt_1=x['rtm']+.001, rt_2=x['rtm']-.0005)
                    return calc_Lscale(nzm, nzt, ncol, gr, False,
                        jnp.full((ncol,), 1e5), jnp.full((ncol,), 1e5),
                        x['pressure'], x['exner'], x['rtm'], x['thlm'], x['thvm'],
                        jnp.full((ncol, nzm), .1), jnp.full((ncol, nzm), 1e-6),
                        jnp.zeros((ncol, nzm)), pdf, x['em'], x['thlm'], 20.,
                        jnp.full((ncol,), .02), jnp.full((ncol,), .01),
                        jnp.full((ncol, nzm), 1e-5), jnp.zeros((ncol, nzt)),
                        .1, 300., jnp.zeros((ncol,)), x['params'], 3, False,
                        method, False, False, False, False,
                        jnp.full((ncol, nzm), 1e-5), stats, ErrInfo.initialized(ncol))

                def loss(params, inputs):
                    """Return mean combined mixing length after replacing the dynamic
                    parameter matrix.
                    """
                    result = calculate(dict(inputs, params=params))
                    return jnp.mean(result[12])

                for operation, function, inputs in [
                    ('forward', jax.jit(calculate), (x,)),
                    ('reverse_params', jax.jit(jax.value_and_grad(loss)), (cp, x)),
                ]:
                    start = time.perf_counter()
                    executable = function.lower(*inputs).compile()
                    compile_seconds = time.perf_counter()-start
                    for _ in range(3):
                        jax.block_until_ready(executable(*inputs))
                    samples = []
                    for _ in range(args.samples):
                        start = time.perf_counter()
                        result = executable(*inputs)
                        jax.block_until_ready(result)
                        samples.append(1000.*(time.perf_counter()-start))
                    assert all(np.isfinite(a).all() for a in jax.tree.leaves(result))
                    memory = executable.memory_analysis()
                    row = dict(method=method, columns=ncol, levels=nzt,
                        operation=operation, compile_seconds=compile_seconds,
                        median_ms=float(np.median(samples)),
                        temporary_bytes=memory.temp_size_in_bytes)
                    if method == 3:
                        row['coherent_layout'] = layout
                    rows.append(row)
                    print(json.dumps(row), flush=True)
            jax.clear_caches()
    if args.output_file:
        args.output_file.write_text(json.dumps(dict(
            jax_version=jax.__version__, devices=[str(d) for d in jax.devices()],
            device_kinds=[d.device_kind for d in jax.devices()],
            coherent_layout=layout,
            profile='thlm gradient .002 K/m, rt=.008, dz=40 m, TKE=.5, mu=.001',
            rows=rows), indent=2)+'\n')


if __name__ == '__main__':
    main()
