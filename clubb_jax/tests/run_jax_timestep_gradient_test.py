#!/usr/bin/env python3
"""Differentiate one real BOMEX/ATEX driver step, including a finite-difference check."""
from pathlib import Path
import argparse
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument('-help', '-h', action='help')
    parser.add_argument('-cases', nargs='+', choices=('bomex', 'atex'), default=['bomex', 'atex'])
    args = parser.parse_args()
    from clubb_jax.run_jax import ensure_environment
    ensure_environment()

import jax
import jax.numpy as jnp
import numpy as np
from utilities.create_case_namelist import create_case_namelist_file
from clubb_jax.src.clubb_case_initalization import init_clubb_case, clean_up_clubb
from clubb_jax.src.advance_clubb_to_end import advance_clubb_to_end
from clubb_jax.src.CLUBB_core.parameter_indices import iC1, iC1b, iC_invrs_tau_shear

INPUTS = ('thlm', 'rtm', 'um', 'vm', 'wp2', 'up2', 'vp2', 'wp3',
          'rtp2', 'thlp2', 'rtpthlp', 'wprtp', 'wpthlp', 'upwp', 'vpwp', 'clubb_params')
OUTPUT_SCALES = {'thlm': 1e-4, 'rtm': 1e4, 'um': 1., 'vm': 1., 'wp2': 1.,
                 'rtp2': 1e8, 'thlp2': 1., 'rcm': 1e6, 'radht': 1e6}


def check_driver_step(state):
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



def main(cases):
    from clubb_jax.src.CLUBB_core import error_code
    for case in cases:
        previous_debug = error_code._debug_level
        state = None
        jax.clear_caches()
        try:
            with tempfile.TemporaryDirectory(prefix='clubb-gradient-') as directory:
                # Tau Lscale avoids reverse-mode-incompatible parcel loops; host I/O stays outside the trace.
                path = create_case_namelist_file(
                    case, Path(directory), stats='none', debug='-1', max_iters=1,
                    override='l_diag_Lscale_from_tau=.true.',
                )
                state = init_clubb_case(str(path))
                assert state['flags'].l_diag_Lscale_from_tau
                assert not state['l_stats']
                assert not error_code.clubb_at_least_debug_level(0)
                check_driver_step(state)
                print(f'{case}: full driver gradient/JVP/finite difference/eager checks passed', flush=True)
        finally:
            if state is not None:
                clean_up_clubb(state)
            error_code._debug_level = previous_debug
            jax.clear_caches()
    return 0


if __name__ == '__main__':
    raise SystemExit(main(args.cases))
