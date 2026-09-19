"""Differentiate the actual forcing/core/radiation sequence with tau-based Lscale."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from utilities.create_case_namelist import create_case_namelist_file
from clubb_jax.src.clubb_case_initalization import init_clubb_case, clean_up_clubb
from clubb_jax.src.advance_clubb_to_end import advance_clubb_to_end
from clubb_jax.src.CLUBB_core.parameter_indices import iC1, iC1b, iC_invrs_tau_shear

INPUTS = ('thlm', 'rtm', 'um', 'vm', 'wp2', 'up2', 'vp2', 'wp3',
          'rtp2', 'thlp2', 'rtpthlp', 'wprtp', 'wpthlp', 'upwp', 'vpwp', 'clubb_params')
OUTPUT_SCALES = {'thlm': 1e-4, 'rtm': 1e4, 'um': 1., 'vm': 1., 'wp2': 1.,
                 'rtp2': 1e8, 'thlp2': 1., 'rcm': 1e6, 'radht': 1e6}


@pytest.fixture(params=['bomex', 'atex'])
def case_state(request, tmp_path, monkeypatch):
    from clubb_jax.src.CLUBB_core import error_code
    # Initialization sets a process-wide debug level. Restore it after this
    # test, and clear compiled traces which captured the previous setting.
    monkeypatch.setattr(error_code, '_debug_level', error_code._debug_level)
    jax.clear_caches()
    # Tau-based Lscale avoids the data-dependent parcel while loops, whose
    # floating-point carry does not support reverse-mode differentiation.
    # stats='none' disables host statistics collection and file I/O. Stats
    # accumulations use stop_gradient, so the loss uses live model fields.
    # debug=-1 excludes even level-0 host error checks: Python cannot branch
    # on traced err_info arrays. We inspect returned errors outside the trace.
    # max_iters=1 limits this regression to a single full driver timestep.
    path = create_case_namelist_file(
        request.param, tmp_path, stats='none', debug='-1', max_iters=1,
        override='l_diag_Lscale_from_tau=.true.',
    )
    state = init_clubb_case(str(path))
    assert state['flags'].l_diag_Lscale_from_tau
    assert not state['l_stats']
    assert not error_code.clubb_at_least_debug_level(0)
    yield state
    clean_up_clubb(state)
    jax.clear_caches()


def test_whole_driver_step_forward_reverse_and_finite_difference(case_state):
    state = case_state
    inputs = {key: jnp.asarray(state[key]) for key in INPUTS}
    original = {key: np.asarray(value).copy() for key, value in inputs.items()}

    def step(x):
        # The ordinary driver updates its dictionary. Copy it locally to keep
        # repeated JAX evaluations independent; JAX arrays are immutable.
        result = dict(state, **x)
        # Suppress progress printing while tracing. max_steps is a static
        # bound, so differentiation follows one complete configured timestep.
        advance_clubb_to_end(result, l_stdout=False, max_steps=1)
        return result

    def loss(x):
        result = step(x)
        total = sum(scale * jnp.mean(jnp.linspace(1., 2., result[key].shape[-1]) * result[key]**2)
                    for key, scale in OUTPUT_SCALES.items())
        return total, result['err_info']

    value_and_grad = jax.jit(jax.value_and_grad(loss, has_aux=True))
    (value, error), gradient = value_and_grad(inputs)
    assert not error.is_fatal()
    assert jnp.isfinite(value)
    for key, grad in gradient.items():
        assert np.isfinite(grad).all(), key
        np.testing.assert_array_equal(state[key], original[key])

    # Neutral initial layers sit exactly at the clipped buoyancy-root kink.
    # Check transpose/finite-difference agreement slightly away from those
    # boundaries; the unperturbed case above must still have finite gradients.
    inputs['thlm'] = inputs['thlm'] + 1e-3 * jnp.sin(jnp.arange(inputs['thlm'].shape[-1]))
    (value, error), gradient = value_and_grad(inputs)
    assert not error.is_fatal()
    for key, grad in gradient.items():
        assert np.isfinite(grad).all(), key

    # All prognostic fields and parameters participate in this transpose check,
    # including fields initialized at exact clipping boundaries.
    direction = {k: jnp.cos(jnp.arange(v.size, dtype=v.dtype)).reshape(v.shape) *
                 (1e-5 if k in ('rtm', 'wprtp') else 1e-8 if k == 'rtp2' else .01)
                 for k, v in inputs.items()}
    scalar_loss = lambda x: loss(x)[0]
    _, tangent = jax.jit(lambda x, v: jax.jvp(scalar_loss, (x,), (v,)))(inputs, direction)
    reverse = sum(jnp.vdot(gradient[k], direction[k]) for k in inputs)
    np.testing.assert_allclose(tangent, reverse, rtol=1e-9, atol=1e-10)

    # Finite differences only perturb smooth inputs, not zero variances or
    # parameters at bounds, where central differences use another convention.
    smooth_direction = {k: jnp.zeros_like(v) for k, v in inputs.items()}
    smooth_direction['thlm'] = direction['thlm']
    smooth_direction['um'] = direction['um']
    # C1 == C1b selects a constant-coefficient branch; move them together.
    smooth_direction['clubb_params'] = smooth_direction['clubb_params'].at[:, iC1].set(.1).at[:, iC1b].set(.1)
    smooth_direction['clubb_params'] = smooth_direction['clubb_params'].at[:, iC_invrs_tau_shear].set(.01)
    exact = sum(jnp.vdot(gradient[k], smooth_direction[k]) for k in inputs)
    evaluate = jax.jit(scalar_loss)
    h = 1e-3
    plus = {k: v + h * smooth_direction[k] for k, v in inputs.items()}
    minus = {k: v - h * smooth_direction[k] for k, v in inputs.items()}
    finite_difference = (evaluate(plus) - evaluate(minus)) / (2 * h)
    np.testing.assert_allclose(exact, finite_difference, rtol=2e-4, atol=1e-7)

    # Compiling the ordinary driver preserves its eager forward outputs.
    def select(x):
        result = step(x)
        return {k: result[k] for k in set(INPUTS) | set(OUTPUT_SCALES)}
    compiled = jax.jit(select)(inputs)
    eager = select(inputs)
    for key in compiled:
        np.testing.assert_allclose(compiled[key], eager[key], rtol=1e-12, atol=1e-14, err_msg=key)


@pytest.mark.parametrize('debug_level', [-1, 0, 1, 2])
def test_driver_error_checks_follow_debug_level(monkeypatch, debug_level):
    import importlib
    from clubb_jax.src.CLUBB_core import error_code
    from clubb_jax.src.CLUBB_core.err_info import ErrInfo
    driver = importlib.import_module('clubb_jax.src.advance_clubb_to_end')
    monkeypatch.setattr(error_code, '_debug_level', debug_level)
    state = dict(dt_main=60., dt_rad=120., time_initial=0., ifinal=3, time_final=180.,
                 l_stats=False, l_calc_thlp2_rad=False, ngrdcol=1, nzm=2, nzt=1,
                 err_info=ErrInfo.initialized(1))
    for key in ('thlm', 'rtm', 'rcm', 'exner', 'thv_ds_zt', 'thlm_forcing', 'radht'):
        state[key] = jnp.ones((1, 1))
    monkeypatch.setattr(driver, 'calculate_thvm', lambda **kwargs: kwargs['thlm'])
    monkeypatch.setattr(driver, '_prescribe_forcings', lambda state, time: None)
    monkeypatch.setattr(driver, '_advance_clubb_core', lambda state: None)
    radiation_calls = []
    def fatal_radiation(state, time_current, l_rad_itime):
        radiation_calls.append(l_rad_itime)
        state['err_info'] = state['err_info'].set_fatal()
    monkeypatch.setattr(driver, '_advance_radiation', fatal_radiation)
    if debug_level >= 0:
        with pytest.raises(RuntimeError, match='Fatal error in radiation'):
            driver.advance_clubb_to_end(state, l_stdout=False)
        assert radiation_calls == [True]
    else:
        driver.advance_clubb_to_end(state, l_stdout=False)
        assert state['err_info'].is_fatal()  # Flags remain available to callers.
        assert radiation_calls == [True, True, False]
