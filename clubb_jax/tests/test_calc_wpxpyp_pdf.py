#!/usr/bin/env python3
"""test_calc_wpxpyp_pdf.py — validate the trinormal triple cross-moment <w'x'y'> (pdf_closure decomposition).

`calc_wpxpyp_pdf` (pdf_closure_module.py ↔ pdf_closure_module.F90:calc_wpxpyp_pdf) is the analytic <w'x'y'> of the
two-component TRINORMAL PDF of (w, x, y) the ADG1 closure uses:
    <w'x'y'> = Σ_i weight_i [ dw·dx·dy + corr_xy·σ_x·σ_y·dw + corr_wy·σ_w·σ_y·dx + corr_wx·σ_w·σ_x·dy ]_i
It feeds calc_pdf_higher_order_moments_jax → pdf_closure_driver and was validated only end-to-end. This pins it (1)
against the closed-form transcription, and (2) an INDEPENDENT Monte-Carlo — sample (w,x,y) from the actual 2-component
trivariate-normal mixture (per-component 3×3 covariance via Cholesky, with the three pairwise correlations) and compare
the empirical triple central moment — validating the three correlation cross-terms, not just transcribing. + finite
grad. Companion to iters 559/560 (calc_wp4_pdf, the binormal moments). (iter 561)
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.pdf_closure_module import calc_wpxpyp_pdf

# per-component (w,x,y) means, std-devs, and pairwise corrs (cwx,cwy,cxy) — chosen positive-definite, O(1)-scaled.
_C1 = dict(mw=1.0, mx=0.7, my=-0.4, sw=0.55, sx=0.63, sy=0.71, cwx=0.3, cwy=0.2, cxy=0.25)
_C2 = dict(mw=-0.5, mx=-0.6, my=0.8, sw=0.71, sx=0.59, sy=0.49, cwx=-0.2, cwy=0.3, cxy=-0.15)
_A = 0.42


def _args():
    wm = _A * _C1['mw'] + (1 - _A) * _C2['mw']
    xm = _A * _C1['mx'] + (1 - _A) * _C2['mx']
    ym = _A * _C1['my'] + (1 - _A) * _C2['my']
    return (wm, xm, ym, _C1['mw'], _C2['mw'], _C1['mx'], _C2['mx'], _C1['my'], _C2['my'],
            _C1['sw'] ** 2, _C2['sw'] ** 2, _C1['sx'] ** 2, _C2['sx'] ** 2, _C1['sy'] ** 2, _C2['sy'] ** 2,
            _C1['cwx'], _C2['cwx'], _C1['cwy'], _C2['cwy'], _C1['cxy'], _C2['cxy'], _A)


def test_closed_form():
    a = _args()
    got = float(calc_wpxpyp_pdf(*a))
    wm, xm, ym = a[0], a[1], a[2]
    ref = 0.0
    for w, C in ((_A, _C1), (1 - _A, _C2)):
        dw, dx, dy = C['mw'] - wm, C['mx'] - xm, C['my'] - ym
        ref += w * (dw * dx * dy + C['cxy'] * C['sx'] * C['sy'] * dw
                    + C['cwy'] * C['sw'] * C['sy'] * dx + C['cwx'] * C['sw'] * C['sx'] * dy)
    assert abs(got - ref) / (abs(ref) + 1e-300) < 1e-12, f"calc_wpxpyp_pdf vs closed form: {got} vs {ref}"
    print(f"  calc_wpxpyp_pdf closed form (rel {abs(got-ref)/(abs(ref)+1e-300):.1e})  PASS")


def _chol(C):
    sw, sx, sy = C['sw'], C['sx'], C['sy']
    cov = np.array([[sw * sw, C['cwx'] * sw * sx, C['cwy'] * sw * sy],
                    [C['cwx'] * sw * sx, sx * sx, C['cxy'] * sx * sy],
                    [C['cwy'] * sw * sy, C['cxy'] * sx * sy, sy * sy]])
    return np.linalg.cholesky(cov)


def test_monte_carlo():
    rng = np.random.default_rng(561)
    N = 12_000_000
    a = _args(); got = float(calc_wpxpyp_pdf(*a)); wm, xm, ym = a[0], a[1], a[2]
    L1, L2 = _chol(_C1), _chol(_C2)
    z = rng.standard_normal((3, N))
    s1 = np.array([[_C1['mw']], [_C1['mx']], [_C1['my']]]) + L1 @ z
    s2 = np.array([[_C2['mw']], [_C2['mx']], [_C2['my']]]) + L2 @ z
    pick = rng.random(N) < _A
    w = np.where(pick, s1[0], s2[0]); x = np.where(pick, s1[1], s2[1]); y = np.where(pick, s1[2], s2[2])
    mc = np.mean((w - wm) * (x - xm) * (y - ym))
    rel = abs(got - mc) / (abs(mc) + 1e-300)
    assert rel < 8e-3, f"calc_wpxpyp_pdf vs MC: analytic {got} vs MC {mc} rel {rel:.2e}"
    print(f"  calc_wpxpyp_pdf vs {N//1_000_000}M-sample trinormal mixture <w'x'y'>: rel {rel:.1e} (<8e-3)  PASS")


def test_grad_finite():
    a = list(_args())
    g = jax.grad(lambda c: calc_wpxpyp_pdf(*a[:15], c, *a[16:]) ** 2)(a[15])   # wrt corr_w_x_1
    assert np.isfinite(float(g)), "non-finite grad wrt corr_w_x_1"
    print("  jax.grad(calc_wpxpyp_pdf) wrt corr_w_x_1 finite  PASS")


def main():
    print("test_calc_wpxpyp_pdf:")
    test_closed_form()
    test_monte_carlo()
    test_grad_finite()
    print("All calc_wpxpyp_pdf checks PASSED")


if __name__ == "__main__":
    main()
