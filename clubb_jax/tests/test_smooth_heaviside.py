#!/usr/bin/env python3
"""test_smooth_heaviside.py — validate the JAX smooth_heaviside_peskin port (advance_helper_module.F90).

Peskin smooth Heaviside (Lin & Lee 2005): H=0 for x<-r, H=1 for x>r, else 0.5(1 + x/r + sin(pi x/r)/pi).
Used in the mixing-length and clip_skewness limiters.

ADVERSARIAL FINDING (iter 62): the JAX uses full-precision pi (jnp.pi); the Fortran constants_clubb.F90 uses
TRUNCATED literals pi = 3.141592654 and invrs_pi = 0.31830988618 (~10-11 digits). So the JAX is *more accurate*
and differs from the f2py oracle by ~6.5e-11 — this is the documented "removed accuracy-lowering contrivances"
choice (DESIGN.md tiered-accuracy standard), NOT a transcription error. This test proves exactly that:
  1. The f2py oracle is reproduced bit-exactly (1e-16) by the TRUNCATED-constant closed form.
  2. The JAX matches the FULL-pi closed form bit-exactly, and its only difference from f2py is the
     constant-precision improvement (~6.5e-11).
  3. Endpoints / monotonicity / finite jax.grad.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for p in (_ROOT, _ROOT + "/clubb_python_api"):
    if p not in sys.path:
        sys.path.append(p)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.advance_helper_module import smooth_heaviside_peskin

_PI_TRUNC = 3.141592654          # constants_clubb.F90 pi
_INVRS_PI_TRUNC = 0.31830988618  # constants_clubb.F90 invrs_pi


def _trunc_ref(x, r):
    xr = x / r
    return np.where(x < -r, 0.0, np.where(x > r, 1.0, 0.5 * (1.0 + xr + _INVRS_PI_TRUNC * np.sin(_PI_TRUNC * xr))))


def test_f2py_is_truncated_constants():
    """The f2py oracle == the truncated-constant closed form (1e-16); the full-pi JAX differs only by the
    documented constant-precision improvement (~1e-10)."""
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py smooth_heaviside oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(4)
    worst_trunc = 0.0
    worst_jax = 0.0
    for smth_range in (0.6, 1.0, 0.25):
        x = rng.uniform(-3 * smth_range, 3 * smth_range, (2, 20))
        ref = np.asarray(clubb_f2py.f2py_smooth_heaviside_peskin(np.asfortranarray(x), smth_range))
        worst_trunc = max(worst_trunc, np.max(np.abs(_trunc_ref(x, smth_range) - ref)))
        worst_jax = max(worst_jax, np.max(np.abs(np.asarray(smooth_heaviside_peskin(x, smth_range)) - ref)))
    assert worst_trunc < 1e-15, f"truncated-const ref should reproduce f2py exactly: {worst_trunc:.2e}"
    assert worst_jax < 1e-9, f"JAX (full pi) vs f2py: {worst_jax:.2e}"
    assert worst_jax > 1e-12, "expected a small nonzero diff from the constant-precision improvement"
    print(f"  f2py == truncated-const closed form ({worst_trunc:.1e}); full-pi JAX more accurate "
          f"(differs {worst_jax:.1e})  PASS")


def test_full_pi_closed_form_and_endpoints():
    r = 0.6
    x = np.linspace(-1.5 * r, 1.5 * r, 41)
    got = np.asarray(smooth_heaviside_peskin(x, r))
    xr = x / r
    ref = np.where(x < -r, 0.0, np.where(x > r, 1.0, 0.5 * (1.0 + xr + np.sin(np.pi * xr) / np.pi)))
    assert np.max(np.abs(got - ref)) < 1e-14, "JAX should match the full-pi closed form"
    assert abs(float(np.asarray(smooth_heaviside_peskin(np.array([[0.0]]), r))[0, 0]) - 0.5) < 1e-14
    assert np.all(np.diff(got) >= -1e-14), "not monotone nondecreasing"
    assert got[0] == 0.0 and got[-1] == 1.0, "tails not saturated"
    print("  JAX matches the full-pi closed form + endpoints (H(0)=0.5, saturated, monotone)  PASS")


def test_differentiable():
    r = 0.6
    g = np.asarray(jax.grad(lambda x: jnp.sum(smooth_heaviside_peskin(x, r) ** 2))(
        jnp.asarray(np.linspace(-0.5, 0.5, 11))))
    assert np.isfinite(g).all(), "non-finite grad through smooth_heaviside_peskin"
    print(f"  jax.grad through smooth_heaviside_peskin: finite ({g.size} entries)  PASS")


def main():
    print("test_smooth_heaviside:")
    for t in (test_f2py_is_truncated_constants, test_full_pi_closed_form_and_endpoints, test_differentiable):
        t()
    print("All smooth_heaviside checks PASSED")


if __name__ == "__main__":
    main()
