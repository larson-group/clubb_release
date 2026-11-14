"""Measure registered algorithms on the five Fortran G_unit objectives.

From the repository root, run:
  bash human_reviewed_tests/run_pytests.sh -unit -k algorithm_objectives -q
  bash human_reviewed_tests/run_pytests.sh -unit -k algorithm_objectives -q --tuner-trials=100 --tuner-seed=7

The terminal summary always shows global-minimum hits/runs, global-hit and
improvement percentages, initial/best losses, median gap and mean evaluation counts;
-s is unnecessary. JUnit properties retain per-trial seeds, points and losses.
The first clubb_tuner unit stage and shared clubb_python unit stage own these
checks. Trial i uses seed --tuner-seed+i; the default runs ten trials.

All methods must produce finite bounded proposals, unique sample IDs and finite
losses, make progress, and finish within their declared evaluation budgets.
Optimizers additionally require at least one trial's best loss strictly below
the loss at that function's predefined starting point; samplers do not require
improvement. This distinction uses only the catalog's category. Global-hit rules
are retained from the retired Fortran suite (see algorithm_objectives.py)
and are diagnostic only: 1/3
means one of three complete searches met that function's coordinate/loss rule,
not one iteration or one function out of three. Formula/minimum and independent
analytic-gradient checks plus finite bounded proposals, unique IDs, progress,
budget and exhaustion checks remain required. Reports are saved even if a
later assertion fails; exceptions during a trial are not suppressed.

Every registered method receives empty options and uses its defaults. Common
execution settings are batch size 8, one worker and one synthetic case; gradients
are evaluated/supplied only when the catalog requires them. Actual evaluation
counts differ by method and are reported without equalizing their budgets.
Small trial counts give coarse percentage resolution. No CLUBB case/JAX runtime
is needed.
"""

from collections import deque
import json

import numpy as np
import pytest

from tuner.algorithm_registry import ALGORITHMS, get_algorithm, prepare_strategy
from tuner.auto_llm_generated_tests.pytests.algorithm_objectives import BENCHMARKS, eggholder
from tuner.tuning_strategy import build_tuning_strategy


def _run_objective(algorithm, objective, seed):
    """Drive the production fill/tell interface with an arbitrary scalar function."""
    names = [f'x{i}' for i in range(len(objective.initial))]
    ranges = [{'name': name, 'min': objective.lower, 'max': objective.upper} for name in names]
    requires_gradients = get_algorithm(algorithm).requires_gradients
    request = {'strategy': {'name': algorithm, 'options': {}},
               'backend': 'jax' if requires_gradients else 'fortran',
               'batch_size': 8, 'max_workers': 1, 'cases': ['synthetic'], 'parameter_ranges': ranges}
    prepare_strategy(request)
    strategy = build_tuning_strategy(request['strategy'], param_names=names,
                                     default_params_row=objective.initial, parameter_ranges=ranges, seed=seed)
    budget = request['total_samples']
    best, count, seen, point = np.inf, 0, set(), None
    while not strategy.is_exhausted():
        queue = deque()
        strategy.fill(queue, 16)
        assert queue or strategy.is_exhausted(), 'Strategy stalled without pending feedback'
        feedback = []
        for sample in queue:
            vector = np.asarray(sample['param_row'])
            assert vector.shape == (len(names),)
            assert np.isfinite(vector).all()
            assert ((vector >= objective.lower) & (vector <= objective.upper)).all()
            assert sample['sample_id'] not in seen
            seen.add(sample['sample_id'])
            result = dict(sample)
            if requires_gradients:
                loss, gradient = objective.evaluate(vector, gradient=True)
                assert gradient.shape == vector.shape and np.isfinite(gradient).all()
                result.update(total_gradient=gradient, gradient_parameters=names)
            else:
                loss = objective.evaluate(vector)
            assert np.isfinite(loss)
            if loss < best:
                best, point = loss, vector.tolist()
            result['total_loss'] = loss
            feedback.append(result)
            count += 1
            assert count <= budget, 'Strategy exceeded its declared evaluation budget'
        strategy.tell(list(reversed(feedback)))
    assert count <= strategy.estimated_sample_count()
    assert 0 < count <= budget
    return {'seed': seed, 'best_loss': best, 'best_point': point, 'evaluations': count,
            'improved': bool(best < objective.evaluate(np.array(objective.initial))),
            'gap': best-objective.minimum,
            'global_minimum_hit': bool(objective.hit_test(point, best))}


@pytest.mark.parametrize('objective', BENCHMARKS, ids=lambda objective: objective.name)
def test_objective_minimum_and_gradient_have_independent_references(objective):
    """Check copied formulas at known minima and analytic gradients by finite differences."""
    value, gradient = objective.evaluate(np.array(objective.minimizer), gradient=True)
    assert value == pytest.approx(objective.minimum, abs=1.e-7)
    assert objective.evaluate(np.array(objective.minimizer)) == value
    # A common fractional position avoids special minima, boundaries and Egg cusps.
    point = np.linspace(.31, .57, len(objective.initial)) * (objective.upper-objective.lower) + objective.lower
    value, gradient = objective.evaluate(point, gradient=True)
    assert objective.evaluate(point) == value
    step = 1.e-6
    finite_difference = np.array([
        (objective.evaluate(point+step*direction) - objective.evaluate(point-step*direction))/(2.*step)
        for direction in np.eye(len(point))
    ])
    assert gradient.shape == finite_difference.shape == point.shape
    np.testing.assert_allclose(gradient, finite_difference, rtol=2.e-5, atol=2.e-5)


def test_eggholder_value_only_evaluation_does_not_require_a_cusp_gradient():
    """Keep a finite scalar objective usable where its derivative is undefined."""
    assert eggholder(np.array([0., -47.])) == 0.
    with pytest.raises(ValueError, match='no finite gradient'):
        eggholder(np.array([0., -47.]), gradient=True)


@pytest.mark.parametrize('name,inside,outside,loss_inside,loss_outside', [
    ('goldstein_price', [.01, -1.], [np.nextafter(.01, np.inf), -1.], 3., 3.),
    ('rastrigin', [0., 0., 0., 0., np.nextafter(.1, 0.)], [0., 0., 0., 0., .1], 0., 0.),
    ('himmelblau', [3., 2.], [3., 2.], np.nextafter(.01, 0.), .01),
    ('eggholder', [512., 405.229], [512., 405.23], -959., -959.),
    ('schaffer', [0., 0.], [0., 0.], 1.e-5, np.nextafter(1.e-5, np.inf)),
])
def test_g_unit_hit_boundaries_match_strict_and_inclusive_rules(name, inside, outside, loss_inside, loss_outside):
    """Distinguish the Fortran hit rules at their coordinate or loss boundaries."""
    objective = next(item for item in BENCHMARKS if item.name == name)
    assert objective.hit_test(inside, loss_inside)
    assert not objective.hit_test(outside, loss_outside)


@pytest.mark.parametrize('algorithm', tuple(ALGORITHMS))
@pytest.mark.parametrize('objective', BENCHMARKS, ids=lambda objective: objective.name)
def test_algorithms_search_fortran_objectives(algorithm, objective, record_property, tuner_benchmark):
    """Check every method's interface; require improvement only for optimizer entries."""
    report = {'algorithm': algorithm, 'category': get_algorithm(algorithm).category,
              'objective': objective.name,
              'initial_loss': objective.evaluate(np.array(objective.initial)),
              'hit_condition': objective.hit_description,
              'requested_trials': tuner_benchmark['trials'], 'trials': [], 'status': 'ERROR'}
    tuner_benchmark['results'].append(report)
    try:
        for seed in range(tuner_benchmark['seed'], tuner_benchmark['seed'] + tuner_benchmark['trials']):
            report['trials'].append(_run_objective(algorithm, objective, seed))
        trials = report['trials']
        report['hits'] = sum(trial['global_minimum_hit'] for trial in trials)
        report['global_hit_percent'] = 100. * report['hits'] / len(trials)
        report['improvement_percent'] = 100. * sum(trial['improved'] for trial in trials) / len(trials)
        report['status'] = 'FAIL'
        if report['category'] == 'optimizer':
            assert any(trial['improved'] for trial in trials), (
                f'{algorithm} / {objective.name}: no trial lowered the starting loss '
                f'{report["initial_loss"]:.6g}; best loss '
                f'{min(trial["best_loss"] for trial in trials):.6g}'
            )
        assert all(trial['gap'] >= -1.e-7 for trial in trials), trials
        report['status'] = 'PASS'
    finally:
        record_property('algorithm_benchmark', json.dumps(report))
