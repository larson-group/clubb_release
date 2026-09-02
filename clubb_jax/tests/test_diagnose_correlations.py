#!/usr/bin/env python3
"""test_diagnose_correlations.py — validate the JAX diagnose_correlations_module port.

  1. f2py bit-shadow: `clubb_f2py.f2py_diagnose_correlations(iipdf_w, corr_pre, l_calc_w_corr)` on the same
     matrix — bit-to-bit (no grid needed). SKIPs cleanly if clubb_f2py is unbuilt.
  2. Analytic checks for the PDF helpers `calc_mean`/`calc_varnce`/`calc_w_corr` (closed-form values + the
     ±0.99 correlation clip).
  3. A `jax.grad` smoke check through the matrix diagnosis.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
# _ROOT first so `import clubb_jax` resolves to this checkout's package.
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for p in (_ROOT, _ROOT + "/clubb_python_api"):
    if p not in sys.path:
        sys.path.append(p)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.diagnose_correlations_module import (
    diagnose_correlations, calc_mean, calc_varnce, calc_w_corr, corr_array_assertion_checks)


def _corr_matrix(n, seed):
    """A symmetric matrix with unit diagonal and off-diagonals in (-0.9, 0.9) — a valid 'prescribed' input."""
    rng = np.random.default_rng(seed)
    a = rng.uniform(-0.9, 0.9, size=(n, n))
    a = 0.5 * (a + a.T)
    np.fill_diagonal(a, 1.0)
    return np.asfortranarray(a)


def test_helpers():
    # calc_mean(a,x1,x2) = a*x1 + (1-a)*x2
    assert abs(float(calc_mean(0.3, 5.0, 2.0)) - (0.3 * 5 + 0.7 * 2)) < 1e-14
    # calc_varnce = a*((x1-xm)^2+x1p2) + (1-a)*((x2-xm)^2+x2p2)
    a, x1, x2, xm, x1p2, x2p2 = 0.4, 3.0, -1.0, 0.6, 0.5, 0.2
    ref = a * ((x1 - xm) ** 2 + x1p2) + (1 - a) * ((x2 - xm) ** 2 + x2p2)
    assert abs(float(calc_varnce(a, x1, x2, xm, x1p2, x2p2)) - ref) < 1e-13
    # calc_w_corr = clip(wpxp/(max(sx,xt)*max(sw,wt)), ±0.99); the clip must fire for a large cov
    assert abs(float(calc_w_corr(0.5, 1.0, 1.0, 1e-2, 1e-2)) - 0.5) < 1e-14
    assert abs(float(calc_w_corr(100.0, 1.0, 1.0, 1e-2, 1e-2)) - 0.99) < 1e-14   # clipped
    assert abs(float(calc_w_corr(-100.0, 1.0, 1.0, 1e-2, 1e-2)) + 0.99) < 1e-14  # clipped
    print("  helpers calc_mean/calc_varnce/calc_w_corr (incl. ±0.99 clip): closed-form  PASS")


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py diagnose_correlations oracle: SKIP ({type(e).__name__})")
        return
    worst = 0.0
    for n, iipdf_w, seed in ((6, 3, 11), (8, 1, 22), (5, 5, 33), (7, 4, 44)):
        corr_pre = _corr_matrix(n, seed)
        ref = np.asarray(clubb_f2py.f2py_diagnose_correlations(iipdf_w, corr_pre.copy(), 0))
        got = np.asarray(diagnose_correlations(n, iipdf_w, corr_pre, l_calc_w_corr=False))
        rel = np.max(np.abs(got - ref) / (np.abs(ref) + 1e-30))
        worst = max(worst, rel)
        assert rel < 1e-12, f"f2py mismatch (n={n}, iipdf_w={iipdf_w}): rel {rel:.2e}"
    print(f"  f2py diagnose_correlations: bit-match over 4 (n, iipdf_w) configs, worst rel {worst:.2e}  PASS")


def test_differentiable():
    n, iipdf_w = 6, 3
    corr = jnp.asarray(_corr_matrix(n, 7))
    g = jax.grad(lambda C: jnp.sum(diagnose_correlations(n, iipdf_w, C) ** 2))(corr)
    assert np.isfinite(np.asarray(g)).all(), "non-finite grad through diagnose_correlations"
    print(f"  jax.grad through diagnose_correlations: finite ({g.size} entries)  PASS")


def test_corr_array_assertion_checks():
    c = _corr_matrix(5, 9)                       # valid: symmetric, unit diag, off-diag in (-0.9,0.9)
    assert corr_array_assertion_checks(c) is True, "valid correlation matrix rejected"
    bad = c.copy(); bad[0, 1] = bad[1, 0] = 1.5  # off-diagonal out of [-0.99, 0.99]
    assert corr_array_assertion_checks(bad) is False, "out-of-range off-diagonal accepted"
    nd = c.copy(); nd[2, 2] = 0.9                # diagonal != 1
    assert corr_array_assertion_checks(nd) is False, "non-unit-diagonal accepted"
    print("  corr_array_assertion_checks: valid PASS / out-of-range + non-unit-diag FAIL  PASS")


def main():
    print("test_diagnose_correlations:")
    for t in (test_helpers, test_f2py_oracle, test_differentiable, test_corr_array_assertion_checks):
        t()
    print("All diagnose_correlations checks PASSED")


if __name__ == "__main__":
    main()
