#!/usr/bin/env python3
"""Check real JAX gradients through a complete bounded Adam tuning job.

Run from the repository root:
  python tuner/auto_llm_generated_tests/integration_tests/run_adam_gradient_test.py -output_dir output/tests/adam_gradient
  python tuner/auto_llm_generated_tests/integration_tests/run_adam_gradient_test.py -output_dir output/tests/adam_forward -gradient_mode forward
  python tuner/auto_llm_generated_tests/integration_tests/run_adam_gradient_test.py -output_dir output/tests/adam_chains -gradient_mode reverse -batch_size 2

The CPU JAX launcher owns environment setup. The test runs BOMEX and ATEX for
10 steps each (0..600 s), two fields and two averaging windows, using one worker,
one chain per batch column, and two Adam updates. -batch_size defaults to one;
larger batches check the automatic chain count and independent column updates.
C11 starts at 0.41 to avoid the model's
C11==C11b shortcut at the default 0.4, whose executed branch has zero C11
derivative. It requires three scored centers per chain, finite case
and total gradients, a nonzero update, and a retained best score no worse than
the first chain's starting score. It verifies each chain's first update against its transported gradient and checks every
initial total gradient against fresh centered-difference loss runs at the same parameters.
The finite-difference step is 1e-5 in each selected physical coordinate;
rtol=3e-3/atol=1e-6 allows numerical differentiation through a short trajectory.
This establishes the bounded workflow, not long-run convergence or GPU support.
The shared Python unit suite covers linked coordinates, all cross-case policies,
ties, fixed ranges and continuation. Jenkins ownership: the JAX gradient Adam
stage of clubb_tuner runs forward mode with two chains; it is separate from the generated pytest wrappers.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def main():
    """Run the production CLI and independent derivative probes, retaining evidence."""
    parser = argparse.ArgumentParser(description=__doc__, add_help=False, allow_abbrev=False)
    parser.add_argument('-help', action='help')
    parser.add_argument('-output_dir', type=Path, required=True)
    parser.add_argument('-gradient_mode', choices=('auto', 'forward', 'reverse'), default='auto')
    parser.add_argument('-batch_size', type=int, default=1)
    args = parser.parse_args()
    if args.batch_size < 1:
        parser.error('-batch_size must be positive')
    from clubb_jax.run_jax import ensure_environment
    ensure_environment()
    from clubb_jax import src  # Configure the CPU compiler before loading JAX.
    import numpy as np
    from tuner.sample_history import sample_history_paths
    from tuner.taylor_metrics import aggregate_losses
    from utilities.clubb_loss_runner import clubb_loss_runner

    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    job = out / 'job'
    command = [sys.executable, str(ROOT / 'run_scripts/run_tuner_job.py'), '-jax=cpu',
               '-cases', 'bomex:0:600:300', 'atex:0:600:300', '-fields', 'wp2', 'thlm',
               '-param_ranges', 'C8:0.2:0.8', 'C11:0.1:1.0',
               '-strategy', 'adam_gradient:max_updates=2:learning_rate=0.01',
               '-batch_size', str(args.batch_size), '-workers', '1', '-gradient_mode', args.gradient_mode,
               '-seed', '7', '-override', 'C11=0.41', '-run_top', 'never', '-output_job_dir', str(job)]
    with (out / 'tuner.log').open('w') as log:
        subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    result = json.loads((job / 'results.json').read_text())
    request = result['request']
    assert result['state'] == 'finished', result
    chains = args.batch_size
    assert request['strategy']['options']['chain_count'] == chains
    assert result['samples_evaluated'] == 3 * chains
    chunks = [dict(np.load(path)) for path in sample_history_paths(job)]
    sample_ids = np.concatenate([chunk['sample_id'] for chunk in chunks])
    order = np.argsort(sample_ids)
    parameters = np.concatenate([chunk['all_params'] for chunk in chunks])[order]
    total_gradients = np.concatenate([chunk['total_gradient'] for chunk in chunks])[order]
    case_gradients = np.concatenate([chunk['gradient'] for chunk in chunks])[order]
    case_losses = np.concatenate([chunk['case_loss'] for chunk in chunks])[order]
    gradient_names = chunks[0]['gradient_parameters'].tolist()
    cases = chunks[0]['case_names'].tolist()
    param_names = chunks[0]['param_names'].tolist()
    assert total_gradients.shape == (3 * chains, 2)
    assert case_gradients.shape == (3 * chains, 2, 2)
    assert np.isfinite(total_gradients).all() and np.isfinite(case_gradients).all()
    assert np.any(total_gradients[0] != 0.)
    spans = np.array([spec['max'] - spec['min'] for spec in request['parameter_ranges']])
    indices = [param_names.index(name) for name in gradient_names]
    minimum = np.array([spec['min'] for spec in request['parameter_ranges']])
    x = (parameters[:chains, indices] - minimum) / spans
    gx = total_gradients[:chains] * spans
    expected = minimum + spans * np.clip(x - .01 * gx / (np.abs(gx) + 1.e-8), 0., 1.)
    np.testing.assert_allclose(parameters[chains:2 * chains, indices], expected, rtol=1.e-12, atol=1.e-12)
    scores = [aggregate_losses(row.tolist(), [request['case_weights'].get(case, 1.) for case in cases],
                               request['aggregation_mode'], request['aggregation_weights'])['loss'] for row in case_losses]
    starting_loss = float(scores[0])
    best_loss = float(result['best_results'][0]['total_loss'])
    np.testing.assert_allclose(best_loss, min(scores), rtol=1.e-12, atol=1.e-12)
    assert best_loss <= starting_loss, (best_loss, starting_loss)

    # Independent finite differences use fresh loss sessions and the same request
    # configuration, without invoking the optimizer or its gradient-combination helper.
    probes = []
    step = 1.e-5
    for case in cases:
        with clubb_loss_runner('jax', case, request['selected_fields'], out / f'probe_{case}',
                               batch_size=args.batch_size, config=request['config'], override=request['override'],
                               case_defaults=request['case_defaults'][case], loss_weights=request['loss_weights'],
                               field_weights=request['field_weights']) as runner:
            assert list(runner.param_names) == param_names
            values = []
            for index in indices:
                pair = []
                for sign in (-1, 1):
                    point = parameters[:chains].copy()
                    point[:, index] += sign * step
                    pair.append(runner.evaluate(point)['loss'].copy())
                values.append(pair)
            probes.append(values)
    probes = np.asarray(probes)
    assert probes.shape == (len(cases), 2, 2, chains)
    finite_difference = np.empty((chains, 2))
    weights = [request['case_weights'].get(case, 1.) for case in cases]
    for chain in range(chains):
        for index in range(2):
            minus = aggregate_losses(probes[:, index, 0, chain].tolist(), weights, request['aggregation_mode'], request['aggregation_weights'])['loss']
            plus = aggregate_losses(probes[:, index, 1, chain].tolist(), weights, request['aggregation_mode'], request['aggregation_weights'])['loss']
            finite_difference[chain, index] = (plus - minus) / (2 * step)
    assert finite_difference.shape == total_gradients[:chains].shape == (chains, 2)
    np.testing.assert_allclose(total_gradients[:chains], finite_difference, rtol=3.e-3, atol=1.e-6)
    summary = {'status': 'passed', 'command': command, 'gradient_mode': args.gradient_mode,
               'cases': cases, 'fields': request['selected_fields'], 'windows_per_case': 2,
               'scores': scores, 'starting_loss': starting_loss, 'best_loss': best_loss,
               'gradient_parameters': gradient_names,
               'chain_count': chains, 'samples_evaluated': result['samples_evaluated'],
               'initial_total_gradient': total_gradients[0].tolist(),
               'finite_difference': finite_difference[0].tolist(),
               'initial_total_gradients': total_gradients[:chains].tolist(),
               'finite_differences': finite_difference.tolist(),
               'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()}
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
