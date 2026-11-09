"""Check runtime sampling, skipped diagnostics and compilation reuse."""

import importlib

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.matrix_solver_wrapper import tridiag_solve
from clubb_jax.src.CLUBB_core.model_flags import iiPDF_ADG1, iiPDF_new, iiPDF_new_hybrid, tridiag_lu
from clubb_jax.src.CLUBB_core import ErrInfo


def _stats(l_sample=False, names=('profile', 'budget')):
    """Create two-column banks for pointwise and budget accumulation checks."""
    return JaxStats.empty(l_sample=l_sample, names=names, ncol=2, max_nlev=3)


def _assert_equal(actual, expected):
    """Compare each leaf without allowing an unintended shape broadcast."""
    assert jax.tree.structure(actual) == jax.tree.structure(expected)
    for left, right in zip(jax.tree.leaves(actual), jax.tree.leaves(expected)):
        assert np.shape(left) == np.shape(right)
        np.testing.assert_array_equal(left, right)


def test_sampling_toggle_reuses_jit_and_preserves_all_off_sample_state():
    """Alternate sampling with one executable; check sums, counts and budgets."""
    @jax.jit
    def step(stats, value):
        """Collect a profile and one closed budget using the public interface."""
        stats = stats.update('profile', value)
        stats = stats.begin_budget('budget', value)
        stats = stats.update_budget('budget', value * 2.)
        return stats.finalize_budget('budget', value * 4.)

    stats = _stats()
    value = jnp.arange(6., dtype=jnp.float64).reshape(2, 3)
    for l_sample in (False, True, False, True, False):
        previous = stats.begin_timestep(l_sample=l_sample, l_reset_accumulators=False)
        stats = step(previous, value)
        if not l_sample:
            _assert_equal(stats, previous)
    assert step._cache_size() == 1
    np.testing.assert_array_equal(stats.buffers[0][0], 2. * value)
    np.testing.assert_array_equal(stats.buffers[0][1], 10. * value)
    np.testing.assert_array_equal(stats.nsamples[0], np.full((2, 2, 3), 2))
    np.testing.assert_array_equal(stats.l_budget, [False, True])
    np.testing.assert_array_equal(stats.l_in_budget, [False, False])
    np.testing.assert_array_equal(stats.issues, [0, 0])
    # A different registry is intentionally a different static layout.
    step(_stats(names=('profile',)), value)
    assert step._cache_size() == 2


@pytest.mark.parametrize('method', ['update', 'begin_budget', 'update_budget', 'finalize_budget'])
def test_off_sample_budget_calls_do_not_record_errors_or_modify_banks(method):
    """Exercise both legal and illegal budget state while runtime sampling is off."""
    value = jnp.ones((2, 3))
    stats = _stats(True).begin_budget('budget', value)
    stats = stats._new(l_sample=False)
    # update and begin are illegal while a budget is open; the other calls
    # would alter its bank or close it. None should have any off-sample effect.
    call = jax.jit(lambda state: getattr(state, method)('budget', value))
    _assert_equal(call(stats), stats)
    stats = _stats(False)
    _assert_equal(call(stats), stats)


def test_runtime_diagnostic_branch_skips_invalid_values_and_their_derivatives():
    """An off-sample negative square root must not contaminate either AD mode."""
    def score(x, l_sample):
        """Collect one nonlinear diagnostic and return its bank sum."""
        stats = _stats(l_sample)

        def collect(stats):
            """Evaluate the diagnostic inside the sampling branch."""
            return stats.update('profile', jnp.full((2, 3), jnp.sqrt(x)))

        stats = stats.sampled(collect)
        return jnp.sum(stats.buffers[0][0])

    reverse = jax.jit(jax.value_and_grad(score))
    for x, flag, expected, derivative in ((-1., False, 0., 0.), (4., True, 12., 1.5)):
        result, gradient = reverse(x, flag)
        forward_result, tangent = jax.jit(
            lambda x, flag: jax.jvp(lambda y: score(y, flag), (x,), (1.,))
        )(x, flag)
        np.testing.assert_allclose([result, forward_result], [expected, expected])
        np.testing.assert_allclose([gradient, tangent], [derivative, derivative])
    assert reverse._cache_size() == 1


def test_unrequested_diagnostics_are_eliminated_without_evaluating_payload():
    """Static registry checks still skip an absent variable with a runtime flag."""
    @jax.jit
    def step(stats):
        """Only construct a deliberately invalid payload for a requested variable."""
        def collect(stats):
            """The absent variable is rejected during tracing."""
            if stats.var_on_stats_list('absent'):
                raise AssertionError('Unrequested diagnostic was traced')
            # Unknown updates also reject payloads without inspecting their shape.
            return stats.update('absent', object())
        return stats.sampled(collect)

    for flag in (False, True):
        stats = _stats(flag)
        _assert_equal(step(stats), stats)
    assert step._cache_size() == 1


def test_runtime_condition_estimate_preserves_solve_and_skips_off_sample_estimate(monkeypatch):
    """Compare optional tridiagonal diagnostics with an independent NumPy oracle."""
    solver = importlib.import_module('clubb_jax.src.CLUBB_core.matrix_solver_wrapper')
    monkeypatch.setattr(solver, 'clubb_at_least_debug_level', lambda level: False)
    estimates = []
    original_estimate = solver._tridiag_rcond_multiple_lhs

    def estimate(ndim, lhs):
        """Record execution of the expensive estimate, including cached calls."""
        jax.debug.callback(lambda: estimates.append(True), ordered=True)
        return original_estimate(ndim, lhs)

    monkeypatch.setattr(solver, '_tridiag_rcond_multiple_lhs', estimate)
    matrix = np.array([[4., -1., 0.], [-1., 4., -1.], [0., -1., 3.]])
    bands = jnp.array([[-1., -1., 0.], [4., 4., 3.], [0., -1., -1.]])
    rhs = jnp.array([[1., 2., 3.]])
    error = ErrInfo.initialized(1)

    @jax.jit
    def solve(flag):
        """Keep matrix shape fixed while changing only diagnostic sampling."""
        return tridiag_solve('sampling-test', tridiag_lu, 1, 3, bands[:, None, :], rhs, error,
                             l_use_rcond=flag)

    for flag in (False, True, False):
        _, solution, rcond = solve(flag)
        assert solution.shape == (1, 3)
        np.testing.assert_allclose(solution[0], np.linalg.solve(matrix, np.asarray(rhs[0])))
        np.testing.assert_allclose(rcond, [1. / np.linalg.cond(matrix, 1) if flag else 0.])
    assert solve._cache_size() == 1
    jax.effects_barrier()
    assert estimates == [True]


def test_checkpoint_cache_reuses_sampling_toggle(monkeypatch):
    """A replayable timestep keeps one cache entry across on/off sampling."""
    driver = importlib.import_module('clubb_jax.src.advance_clubb_to_end')

    def step(state, itime, time_current, l_sample, l_rad_itime):
        """Advance a small differentiable state while recording sampled values."""
        state['l_sample'] = l_sample
        state['x'] = state['x'] * 2.
        state['_jax_stats'] = state['_jax_stats'].update('profile', state['x'])

    monkeypatch.setattr(driver, '_advance_timestep', step)
    state = dict(x=jnp.ones((2, 3)), _jax_stats=_stats(), l_sample=False)
    cache = {}
    for itime, flag in enumerate((False, True, False, True), 1):
        state['_jax_stats'] = state['_jax_stats'].begin_timestep(
            l_sample=flag, l_reset_accumulators=False)
        driver._checkpoint_timestep(state, itime, float(itime - 1), flag, False, cache)
    assert len(cache) == 1
    np.testing.assert_array_equal(state['_jax_stats'].buffers[0][0], np.full((2, 3), 20.))
    np.testing.assert_array_equal(state['_jax_stats'].nsamples[0][0], np.full((2, 3), 2))


@pytest.mark.parametrize('pdf_type,l_upwind,l_godunov,l_explicit', [
    (iiPDF_ADG1, False, False, False), (iiPDF_ADG1, True, False, False),
    (iiPDF_ADG1, True, True, False),
    (iiPDF_new, False, False, False), (iiPDF_new, True, False, False),
    (iiPDF_new_hybrid, False, False, False), (iiPDF_new_hybrid, True, False, False),
    (iiPDF_ADG1, False, False, True), (iiPDF_ADG1, True, False, True),
])
def test_sampling_preserves_physical_advection_terms(monkeypatch, pdf_type, l_upwind, l_godunov, l_explicit):
    """Sampling changes diagnostic banks, never the scalar or primary LHS/RHS.

    Exercise the ADG1, new-PDF, hybrid, Godunov and explicit branches of
    advance_xp2_xpyp_module.F90:calc_xp2_xpyp_ta_terms. Nonzero scalar fields
    cover callback temporaries shared by physics and sampled diagnostics.
    """
    module = importlib.import_module('clubb_jax.src.CLUBB_core.advance_xp2_xpyp_module')
    from clubb_jax.src.CLUBB_core.grid_class import setup_grid
    from clubb_jax.src.CLUBB_core.pdf_params import init_pdf_implicit_coefs_terms_api

    monkeypatch.setattr(module, 'l_explicit_turbulent_adv_xpyp', l_explicit)
    grid = setup_grid(ngrdcol=2, deltaz=100., zm_init=0., zm_top=500., grid_type=1)
    nc, nt, nm = grid.ngrdcol, grid.nzt, grid.nzm
    args = dict(nzm=nm, nzt=nt, ngrdcol=nc, sclr_dim=1, gr=grid)
    for name in ('wprtp', 'wpthlp', 'rtp2', 'thlp2', 'rtpthlp', 'upwp', 'vpwp',
                 'up2', 'vp2', 'wp2', 'wp3_on_wp2'):
        args[name] = jnp.full((nc, nm), .2)
    for name in ('wprtp2', 'wpthlp2', 'wprtpthlp', 'wp2_zt', 'wp3_on_wp2_zt'):
        args[name] = jnp.full((nc, nt), .2)
    for name in ('wpsclrp', 'sclrp2', 'sclrprtp', 'sclrpthlp'):
        args[name] = jnp.full((nc, nm, 1), .2)
    for name in ('wpsclrp2', 'wpsclrprtp', 'wpsclrpthlp'):
        args[name] = jnp.full((nc, nt, 1), .2)
    args.update(rho_ds_zt=jnp.ones((nc, nt)), rho_ds_zm=jnp.ones((nc, nm)),
                invrs_rho_ds_zm=jnp.ones((nc, nm)), sigma_sqd_w=jnp.full((nc, nm), .1),
                beta=jnp.full((nc,), .8), l_scalar_calc=True, iiPDF_type=pdf_type,
                l_upwind_xpyp_ta=l_upwind, l_godunov_upwind_xpyp_ta=l_godunov,
                pdf_implicit_coefs_terms=jax.tree.map(
                    lambda field: jnp.full_like(field, .15),
                    init_pdf_implicit_coefs_terms_api(nt, nc, 1)))
    names = ('coef_wprtp2_implicit', 'term_wprtp2_explicit',
             'coef_wpthlp2_implicit', 'term_wpthlp2_explicit',
             'coef_wprtpthlp_implicit', 'term_wprtpthlp_explicit')
    stats = JaxStats.empty(l_sample=False, names=names, ncol=nc, max_nlev=nt)
    calculate = jax.jit(lambda stats: module.calc_xp2_xpyp_ta_terms(**args, stats=stats))
    off = calculate(stats)
    on = calculate(stats.begin_timestep(l_sample=True, l_reset_accumulators=False))
    assert calculate._cache_size() == 1
    _assert_equal(off[:-1], on[:-1])
    _assert_equal(off[-1], stats)
    for field in on[:-1]:
        assert np.isfinite(field).all()
    np.testing.assert_array_equal(on[-1].nsamples[0], np.ones((len(names), nc, nt)))
