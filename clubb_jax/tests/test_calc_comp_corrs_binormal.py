#!/usr/bin/env python3
"""test_calc_comp_corrs_binormal.py — validate the JAX calc_comp_corrs_binormal + smooth_corr_quotient ports.

calc_comp_corrs_binormal (pdf_utilities.F90:766) diagnoses the shared PDF-component correlation of two binormal
variables x and y from their overall covariance <x'y'>, bounded by smooth_corr_quotient. Oracles:
  1. f2py bit-shadow: clubb_f2py.f2py_calc_comp_corrs_binormal on the same fields — both outputs match.
     SKIPs cleanly if clubb_f2py is unbuilt.
  2. Round-trip (oracle-free): build <x'y'> from a chosen corr via the forward covariance formula, recover corr.
  3. smooth_corr_quotient keeps |corr| <= max_mag_correlation even for a huge covariance.
  4. A finite jax.grad.
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

from clubb_jax.src.CLUBB_core.pdf_utilities import calc_comp_corrs_binormal, smooth_corr_quotient

NG, NZ = 2, 6


def _fields(seed):
    rng = np.random.default_rng(seed)
    mu_x_1 = rng.uniform(-2, 2, (NG, NZ)); mu_x_2 = rng.uniform(-2, 2, (NG, NZ))
    mu_y_1 = rng.uniform(-2, 2, (NG, NZ)); mu_y_2 = rng.uniform(-2, 2, (NG, NZ))
    sx1 = rng.uniform(0.1, 2.0, (NG, NZ)); sx2 = rng.uniform(0.1, 2.0, (NG, NZ))
    sy1 = rng.uniform(0.1, 2.0, (NG, NZ)); sy2 = rng.uniform(0.1, 2.0, (NG, NZ))
    a = rng.uniform(0.2, 0.8, (NG, NZ))
    xm = a * mu_x_1 + (1 - a) * mu_x_2
    ym = a * mu_y_1 + (1 - a) * mu_y_2
    xpyp = rng.uniform(-1, 1, (NG, NZ))
    return dict(xpyp=xpyp, xm=xm, ym=ym, mu_x_1=mu_x_1, mu_x_2=mu_x_2, mu_y_1=mu_y_1, mu_y_2=mu_y_2,
                sigma_x_1_sqd=sx1, sigma_x_2_sqd=sx2, sigma_y_1_sqd=sy1, sigma_y_2_sqd=sy2, mixt_frac=a)


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calc_comp_corrs_binormal oracle: SKIP ({type(e).__name__})")
        return
    worst = 0.0
    for seed in (11, 22, 33):
        f = _fields(seed)
        r1, r2 = clubb_f2py.f2py_calc_comp_corrs_binormal(
            f['xpyp'], f['xm'], f['ym'], f['mu_x_1'], f['mu_x_2'], f['mu_y_1'], f['mu_y_2'],
            f['sigma_x_1_sqd'], f['sigma_x_2_sqd'], f['sigma_y_1_sqd'], f['sigma_y_2_sqd'], f['mixt_frac'])
        g1, g2 = calc_comp_corrs_binormal(**f)
        d = max(np.max(np.abs(np.asarray(g1) - np.asarray(r1))),
                np.max(np.abs(np.asarray(g2) - np.asarray(r2))))
        worst = max(worst, d)
        assert d < 1e-12, f"f2py mismatch seed {seed}: {d:.2e}"
    print(f"  f2py calc_comp_corrs_binormal: bit-match (both comps) over 3 seeds, worst {worst:.2e}  PASS")


def test_round_trip():
    rng = np.random.default_rng(5)
    worst = 0.0
    for _ in range(100):
        a = rng.uniform(0.3, 0.7)
        mu_x_1, mu_x_2, mu_y_1, mu_y_2 = rng.uniform(-2, 2, 4)
        sx1, sx2, sy1, sy2 = rng.uniform(0.3, 1.5, 4)   # std devs
        xm = a * mu_x_1 + (1 - a) * mu_x_2
        ym = a * mu_y_1 + (1 - a) * mu_y_2
        corr = rng.uniform(-0.9, 0.9)
        # Forward covariance assembly.
        xpyp = (a * ((mu_x_1 - xm) * (mu_y_1 - ym) + corr * sx1 * sy1)
                + (1 - a) * ((mu_x_2 - xm) * (mu_y_2 - ym) + corr * sx2 * sy2))
        g1, g2 = calc_comp_corrs_binormal(
            np.array([[xpyp]]), np.array([[xm]]), np.array([[ym]]),
            np.array([[mu_x_1]]), np.array([[mu_x_2]]), np.array([[mu_y_1]]), np.array([[mu_y_2]]),
            np.array([[sx1 ** 2]]), np.array([[sx2 ** 2]]), np.array([[sy1 ** 2]]), np.array([[sy2 ** 2]]),
            np.array([[a]]))
        v1 = float(np.asarray(g1).item()); v2 = float(np.asarray(g2).item())
        worst = max(worst, abs(v1 - corr), abs(v1 - v2))
    assert worst < 1e-9, f"round-trip recovery {worst:.2e}"
    print(f"  round-trip: recover corr from assembled <x'y'>, worst {worst:.2e}  PASS")


def test_bound():
    # Huge covariance with tiny variances -> quotient would blow up; smoothing must cap |corr| <= 0.99.
    q = float(np.asarray(smooth_corr_quotient(np.array([[1.0e6]]), np.array([[1.0e-3]]), 1.0e-10)).item())
    assert abs(q) <= 0.99 + 1e-9, f"smooth_corr_quotient exceeded max_mag_correlation: {q}"
    qn = float(np.asarray(smooth_corr_quotient(np.array([[-1.0e6]]), np.array([[1.0e-3]]), 1.0e-10)).item())
    assert abs(qn) <= 0.99 + 1e-9
    print(f"  smooth_corr_quotient bound: |corr| <= 0.99 for a huge covariance ({q:.4f})  PASS")


def test_smooth_corr_quotient_f2py():
    """smooth_corr_quotient (pdf_utilities — the smoothly |corr|<=max_mag clamped numerator/denominator) was only
    checked by the bound property above; validate it directly against the f2py oracle on general num/denom. SKIPs
    if clubb_f2py is unbuilt. (iter 420)"""
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py smooth_corr_quotient oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(4)
    worst = 0.0
    for _ in range(40):
        ng, nz = 2, 8
        num = rng.uniform(-1e3, 1e3, (ng, nz)); den = rng.uniform(1e-4, 1e2, (ng, nz)); thr = 1e-10
        j = np.asarray(smooth_corr_quotient(jnp.asarray(num), jnp.asarray(den), thr))
        f = np.asarray(clubb_f2py.f2py_smooth_corr_quotient(num, den, thr))
        worst = max(worst, float(np.max(np.abs(j - f))))
    assert worst < 1e-12, f"smooth_corr_quotient f2py mismatch {worst:.2e}"
    print(f"  f2py smooth_corr_quotient (40 cases, incl. clamped): bit-match, worst {worst:.2e}  PASS")


def test_differentiable():
    f = _fields(7)
    def loss(xpyp):
        g1, g2 = calc_comp_corrs_binormal(xpyp, f['xm'], f['ym'], f['mu_x_1'], f['mu_x_2'], f['mu_y_1'],
                                          f['mu_y_2'], f['sigma_x_1_sqd'], f['sigma_x_2_sqd'],
                                          f['sigma_y_1_sqd'], f['sigma_y_2_sqd'], f['mixt_frac'])
        return jnp.sum(g1 ** 2 + g2 ** 2)
    g = np.asarray(jax.grad(loss)(jnp.asarray(f['xpyp'])))
    assert np.isfinite(g).all(), "non-finite grad through calc_comp_corrs_binormal"
    print(f"  jax.grad through calc_comp_corrs_binormal: finite ({g.size} entries)  PASS")


def main():
    print("test_calc_comp_corrs_binormal:")
    for t in (test_f2py_oracle, test_round_trip, test_bound, test_smooth_corr_quotient_f2py, test_differentiable):
        t()
    print("All calc_comp_corrs_binormal checks PASSED")


if __name__ == "__main__":
    main()
