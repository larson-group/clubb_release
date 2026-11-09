"""Independent gradient and optimizer checks without model runs.

These provisional tests belong to the shared Python pytest suite (-unit).
Finite differences exercise the public cross-case score; quadratic objectives
exercise linked coordinates, bounded Adam feedback and serialized continuation.
The real JAX tuner workflow is separate under integration_tests.
"""

from collections import deque
import json
import pickle

import numpy as np
import pytest

from tuner.algorithms.adam_gradient import AdamGradientStrategy
from tuner.algorithm_registry import (
    ALGORITHMS, _read_catalog, apply_requirements, get_algorithm,
    list_algorithms, normalize_options, parse_strategy, prepare_strategy,
)
from tuner.gradients import combine_case_gradients, range_gradient
from tuner.taylor_metrics import aggregate_losses
from tuner.tuning_strategy import build_tuning_strategy, normalize_strategy_config


@pytest.mark.parametrize("mode", ["mean_max", "mean_worst_quantile", "quantile_weighted"])
def test_combined_gradient_matches_public_score_finite_difference(mode):
    """Differentiate a weighted seven-case objective away from ranking boundaries."""
    intercept = np.array([1., 3., 5., 8., 12., 16., 21.])
    jacobian = np.array([[.2, -.7], [1., 2.], [-3., .6], [.4, 1.], [2., -2.], [8., 1.], [-1., 5.]])
    case_weights = [1., 0., 2., 4., 3., .5, 2.]
    theta = np.array([.1, -.05])

    def score(point):
        """Use the production score with fixed synthetic case derivatives."""
        return aggregate_losses((intercept + jacobian @ point).tolist(), case_weights, mode)

    result = score(theta)
    actual = combine_case_gradients(result["loss_coefficients"], jacobian)
    h = 1.e-6
    expected = np.array([(score(theta + h * direction)["loss"] - score(theta - h * direction)["loss"]) / (2 * h)
                         for direction in np.eye(2)])
    assert actual.shape == expected.shape == (2,)
    np.testing.assert_allclose(actual, expected, rtol=2.e-8, atol=2.e-9)
    assert result["loss_coefficients"][1] == 0.
    assert sum(result["loss_coefficients"]) == pytest.approx(1.)


def test_ties_and_excluded_case_gradients_have_explicit_behavior():
    """Max ties share influence; stable ranking and zero-weight NaNs stay deterministic."""
    result = aggregate_losses([2., 2., 99.], [1., 3., 0.], "mean_max")
    np.testing.assert_allclose(result["loss_coefficients"], [.35, .65, 0.])
    np.testing.assert_allclose(combine_case_gradients(result["loss_coefficients"], [[2.], [4.], [np.nan]]), [3.3])
    ranked = aggregate_losses([2., 2., 2., 2.], aggregation_mode="quantile_weighted")
    np.testing.assert_allclose(ranked["loss_coefficients"], [.1, .4, .4, .1])
    with pytest.raises(ValueError, match="nonfinite"):
        combine_case_gradients([1., 0.], [[np.nan], [2.]])
    with pytest.raises(ValueError, match="shape"):
        combine_case_gradients([1.], [[1.], [2.]])


def test_linked_and_scaled_gradient_matches_coordinate_finite_difference():
    """Linked physical partials add, range scaling multiplies, and fixed coordinates vanish."""
    ranges = [{"name": "linked", "targets": ["A", "B"], "min": 2., "max": 6.},
              {"name": "C", "min": 1., "max": 1.}]
    x = .3
    p = 2. + 4. * x
    # f(A,B)=A^2+3B^2, evaluated along A=B=p(x).
    actual = range_gradient([6. * p, 999., 2. * p], ["B", "C", "A"], ranges)
    h = 1.e-6
    expected = (4. * (2. + 4. * (x + h))**2 - 4. * (2. + 4. * (x - h))**2) / (2 * h)
    np.testing.assert_allclose(actual, [expected, 0.], rtol=1.e-9)
    with pytest.raises(ValueError, match="Missing"):
        range_gradient([1.], ["A"], ranges)


def _strategy(**options):
    """Create a bounded one-parameter quadratic search with configurable chains."""
    return AdamGradientStrategy(param_names=["A"], default_params_row=[1.5],
                                parameter_ranges=[{"name": "A", "min": 0., "max": 2.}],
                                max_updates=options.pop("max_updates", 12),
                                learning_rate=.05, **options)


def _feedback(sample):
    """Differentiate the quadratic independently in physical parameter units."""
    a = sample["param_row"][0]
    return {"sample_id": sample["sample_id"], "total_loss": (a - .4)**2,
            "total_gradient": [2. * (a - .4)], "gradient_parameters": ["A"]}


def test_adam_uses_gradient_moves_and_scores_the_final_center():
    """The first normalized Adam step is exact, reduces a quadratic and stays bounded."""
    strategy = _strategy(max_updates=4)
    samples = []
    while not strategy.is_exhausted():
        queue = deque()
        strategy.fill(queue, 1)
        assert len(queue) == 1
        sample = queue.popleft()
        samples.append(sample)
        strategy.tell([_feedback(sample)])
    assert len(samples) == strategy.estimated_sample_count() == 5
    np.testing.assert_allclose(samples[1]["param_row"], [1.4], atol=1.e-8)
    assert _feedback(samples[-1])["total_loss"] < _feedback(samples[0])["total_loss"]
    assert all(0. <= sample["param_row"][0] <= 2. for sample in samples)


def test_adam_partial_out_of_order_feedback_survives_checkpoint():
    """An outstanding chain and advanced chains resume with identical proposals."""
    strategy = _strategy(chain_count=3, seed=12)
    queue = deque()
    strategy.fill(queue, 3)
    initial = list(queue)
    strategy.tell([_feedback(initial[2]), _feedback(initial[0])])
    resumed = pickle.loads(pickle.dumps(strategy))
    strategy.tell([_feedback(initial[1])])
    resumed.tell([_feedback(initial[1])])
    left, right = deque(), deque()
    strategy.fill(left, 3)
    resumed.fill(right, 3)
    assert list(left) == list(right)
    assert len(left) == 3
    assert strategy.pending_batch_limit(2) == 2


def test_adam_latin_hypercube_starts_cover_active_ranges_and_preserve_links():
    """Additional chains cover each stratum, reproduce their seed and keep fixed values."""
    options = dict(param_names=['A', 'B', 'C', 'D'], default_params_row=[3., -2., 3., 9.],
                   parameter_ranges=[{'name': 'linked', 'targets': ['A', 'B'], 'min': -1., 'max': 1.},
                                     {'name': 'C', 'min': 2., 'max': 4.},
                                     {'name': 'D', 'min': 7., 'max': 7.}],
                   chain_count=9, seed=17)
    strategy = AdamGradientStrategy(**options)
    starts = np.stack([chain.x for chain in strategy.chains])
    assert starts.shape == (9, 3)
    np.testing.assert_array_equal(starts[0], [1., .5, 0.])
    np.testing.assert_array_equal(starts[:, 2], np.zeros(9))
    for coordinate in range(2):
        np.testing.assert_array_equal(np.sort(np.floor(starts[1:, coordinate] * 8)), np.arange(8))
    repeated = AdamGradientStrategy(**options)
    np.testing.assert_array_equal(starts, np.stack([chain.x for chain in repeated.chains]))
    changed = AdamGradientStrategy(**{**options, 'seed': 18})
    assert not np.array_equal(starts[1:], np.stack([chain.x for chain in changed.chains])[1:])
    queue = deque()
    strategy.fill(queue, 9)
    assert len(queue) == 9
    for sample in queue:
        a, b, c, d = sample['param_row']
        assert a == b and -1. <= a <= 1.
        assert 2. <= c <= 4. and d == 7.


@pytest.mark.parametrize('workers,batch_size', [(1, 1), (1, 4), (3, 2)])
@pytest.mark.parametrize('override', [None, 1, 5])
def test_adam_chain_count_comes_from_capacity_even_for_saved_explicit_counts(workers, batch_size, override):
    """Canonical requests recompute chain counts and budgets from execution capacity."""
    options = {'max_updates': 2}
    if override is not None:
        options['chain_count'] = override
    request = {'strategy': {'name': 'adam_gradient', 'options': options}, 'backend': 'jax',
               'cases': ['bomex', 'atex'], 'batch_size': batch_size, 'max_workers': workers,
               'parameter_ranges': [{'name': 'A', 'min': 0., 'max': 2.}]}
    prepare_strategy(request)
    expected = workers * batch_size
    assert request['strategy']['options']['chain_count'] == expected
    assert request['total_samples'] == 3 * expected
    built = build_tuning_strategy(request['strategy'], param_names=['A'], default_params_row=[1.],
                                  parameter_ranges=request['parameter_ranges'], seed=12)
    assert len(built.chains) == expected
    assert built.estimated_sample_count() == request['total_samples']
    prepare_strategy(request)
    assert request['strategy']['options']['chain_count'] == expected
    request['parameter_ranges'][0]['max'] = 0.
    prepare_strategy(request)
    assert request['total_samples'] == 1


def test_adam_fails_on_bad_gradients_and_evaluates_fixed_ranges_once():
    """Reject an unusable update; all-fixed requests require only a scored center."""
    strategy = _strategy()
    queue = deque()
    strategy.fill(queue, 1)
    feedback = _feedback(queue[0])
    feedback["total_gradient"] = [np.nan]
    with pytest.raises(ValueError, match="Nonfinite"):
        strategy.tell([feedback])
    fixed = AdamGradientStrategy(param_names=["A"], default_params_row=[2.],
                                 parameter_ranges=[{"name": "A", "min": 1., "max": 1.}], chain_count=4)
    queue = deque()
    fixed.fill(queue, 5)
    assert len(queue) == fixed.estimated_sample_count() == 1
    fixed.tell([_feedback(queue[0])])
    assert fixed.is_exhausted()


def test_registry_drives_cli_factory_requirements_and_option_validation():
    """A registered AD strategy is constructible and cannot silently run value-only."""
    strategy = parse_strategy("adam_gradient:max_updates=2:learning_rate=0.03")
    request = {"strategy": strategy, "backend": "jax", "cases": ["bomex"],
               "batch_size": 1, "max_workers": 1, "parameter_ranges": [{"name": "A", "min": 0., "max": 2.}]}
    prepare_strategy(request)
    assert request["grad"] is True
    assert request["total_samples"] == 3
    built = build_tuning_strategy(request["strategy"], param_names=["A"], default_params_row=[1.],
                                  parameter_ranges=request["parameter_ranges"])
    assert isinstance(built, AdamGradientStrategy)
    with pytest.raises(ValueError, match="jax"):
        apply_requirements({"strategy": strategy, "backend": "fortran"})
    for raw in ({"beta1": 1.}, {"epsilon": 0.}, {"max_updates": 1.5}, {"learning_rate": float("inf")}, {"unknown": 1}):
        with pytest.raises(ValueError):
            normalize_options("adam_gradient", raw)
    assert "adam_gradient" in json.dumps(list_algorithms())


def test_new_registration_supports_types_without_frontend_switches(tmp_path, monkeypatch):
    """A new Python file and JSON entry supply parsing, budgets and construction."""
    module_name = "synthetic_tuner_algorithm"
    (tmp_path / (module_name + ".py")).write_text(
        "from tuner.algorithms.random import RandomUniformStrategy\n"
        "class SyntheticStrategy(RandomUniformStrategy):\n"
        "    def __init__(self, count=2, enabled=True, variant='a', **kwargs):\n"
        "        super().__init__(max_samples=count, **kwargs)\n"
        "        self.enabled, self.variant = enabled, variant\n"
        "    @staticmethod\n"
        "    def sample_count(options, ranges):\n"
        "        return options['count']\n"
    )
    catalog_path = tmp_path / "algorithms.json"
    catalog_path.write_text(json.dumps({"synthetic": {
        "implementation": f"{module_name}:SyntheticStrategy",
        "category": "sampler",
        "options": {"count": {"type": "int", "default": 2, "minimum": 1},
                    "enabled": {"type": "bool", "default": True},
                    "variant": {"type": "string", "default": "a", "choices": ["a", "b"]}},
        "cli_order": ["count", "enabled", "variant"]}}))
    monkeypatch.syspath_prepend(str(tmp_path))
    monkeypatch.setitem(ALGORITHMS, "synthetic", _read_catalog(catalog_path)["synthetic"])
    get_algorithm.cache_clear()
    try:
        config = parse_strategy("synthetic:3:false:b")
        assert config["options"] == {"count": 3, "enabled": False, "variant": "b"}
        request = {"strategy": config, "max_workers": 1, "batch_size": 1,
                   "cases": ["case"], "parameter_ranges": []}
        prepare_strategy(request)
        assert request["total_samples"] == 3
        strategy = build_tuning_strategy(config, param_names=[], default_params_row=[], parameter_ranges=[])
        assert strategy.enabled is False and strategy.variant == "b"
        queue = deque()
        strategy.fill(queue, 10)
        assert len(queue) == 3 and strategy.is_exhausted()
    finally:
        get_algorithm.cache_clear()


def test_scheduler_combines_out_of_order_case_results_before_feedback_and_history(tmp_path):
    """Two worker results produce the exact total gradient consumed by Adam and saved to NPZ."""
    from tuner.sample_history import sample_history_paths
    from tuner.taylor_metrics import LOSS_METRIC_NAMES
    from tuner.tuning_scheduler import TuningScheduler

    request = {"backend": "jax", "cases": ["bomex", "atex"], "selected_fields": ["wp2"],
               "batch_size": 1, "max_workers": 1, "aggregation_mode": "mean_max",
               "strategy": {"name": "adam_gradient", "options": {"max_updates": 1}},
               "parameter_ranges": [{"name": "C8", "min": .2, "max": .8},
                                    {"name": "C11", "min": .1, "max": 1.}]}
    scheduler = TuningScheduler(request=request, job_dir=tmp_path, control_path=tmp_path / 'control.json',
                                status_path=tmp_path / 'status.json', results_path=tmp_path / 'results.json')
    scheduler.param_names = ['C8', 'C11']
    scheduler.default_params_row = [.5, .6]
    scheduler.strategy = build_tuning_strategy(request['strategy'], param_names=scheduler.param_names,
                                              default_params_row=scheduler.default_params_row,
                                              parameter_ranges=request['parameter_ranges'])
    scheduler._initialize_sample_history_writer()
    scheduler._fill_pending_samples()
    scheduler._pack_pending_batches()
    batch_id = next(iter(scheduler.batches))
    for case, loss, gradient in [('atex', 3., [7., 11.]), ('bomex', 1., [2., 5.])]:
        scheduler._handle_result({
            'batch_id': batch_id, 'case_name': case, 'field_names': ['wp2'],
            'loss_metric_names': list(LOSS_METRIC_NAMES),
            'loss_metrics_by_metric_window_field_and_column': {name: [[[loss]]] for name in LOSS_METRIC_NAMES},
            'case_loss': [loss], 'field_loss': [[loss]], 'gradient_target': 'loss',
            'gradient_parameters': ['C8', 'C11'], 'gradient': {'loss': [gradient]},
        })
        assert not scheduler.error_message
    completed = scheduler._drain_completed_samples()
    assert len(completed) == 1
    expected = np.array([.3 * 2. + .7 * 7., .3 * 5. + .7 * 11.])
    np.testing.assert_allclose(completed[0]['total_gradient'], expected)
    assert scheduler.strategy.chains[0].updates == 1
    assert scheduler.strategy.chains[0].pending is None
    scheduler._close_sample_history()
    with np.load(sample_history_paths(tmp_path)[0]) as data:
        assert data['gradient'].shape == (1, 2, 2)
        assert data['total_gradient'].shape == (1, 2)
        np.testing.assert_allclose(data['total_gradient'][0], expected)
        np.testing.assert_array_equal(data['gradient'][0], [[2., 5.], [7., 11.]])


def test_retained_results_keep_the_total_gradient(tmp_path):
    """Full results preserve the objective derivative and its physical coordinate order."""
    from tuner.status import write_results
    entry = {"sample_id": 0, "total_loss": 2., "selected_params": {"A": .5}, "all_params": {"A": .5},
             "gradient": {"bomex": {"loss": [3.]}}, "gradient_parameters": ["A"], "total_gradient": [3.]}
    path = tmp_path / 'results.json'
    write_results(path, state='finished', job_dir=tmp_path, request={}, samples_evaluated=1,
                  best_results=[entry], started_at='now', updated_at='now')
    retained = json.loads(path.read_text())['best_results'][0]
    assert retained['total_gradient'] == [3.]
    assert retained['gradient_parameters'] == ['A']


def test_grid_budget_does_not_materialize_the_requested_grid(monkeypatch):
    """A trillion-point preview reports its size without allocating candidate arrays."""
    from tuner.algorithms.resolve import ResolveGridStrategy

    def reject_allocation(*args):
        """Fail if a preview tries to construct points instead of counting them."""
        raise AssertionError('budget preview materialized the grid')

    monkeypatch.setattr(ResolveGridStrategy, '_values_for_range', reject_allocation)
    request = {'strategy': {'name': 'resolve', 'options': {'spacing': 1.e-12}},
               'parameter_ranges': [{'min': 0., 'max': 1.}]}
    prepare_strategy(request)
    assert request['total_samples'] == 1_000_000_000_001
