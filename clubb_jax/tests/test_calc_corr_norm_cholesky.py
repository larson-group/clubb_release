#!/usr/bin/env python3
"""test_calc_corr_norm_cholesky.py — validate the JAX calc_corr_norm_and_cholesky_factor orchestration.

calc_corr_norm_and_cholesky_factor (setup_clubb_pdf_params.F90:1070) adjusts the prescribed in-cloud and
below-cloud correlation matrices (ADG zeroing, Ncn in-cloud override, eta–hm product estimate), Cholesky-
factorizes each once, then assigns per grid column/level by rc. Oracle:
  1. Reconstruction: at every (i,k) the assigned Cholesky factor L satisfies L Lᵀ == the assigned corr array
     (Cholesky_factor itself is bit-validated vs f2py in test_cholesky_factor).
  2. Structure: corr arrays symmetric + unit diagonal; Cholesky lower-triangular.
  3. The matrix adjustments + rc-selection vs a literal NumPy reference.
  4. A finite jax.grad.
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
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.setup_clubb_pdf_params import calc_corr_norm_and_cholesky_factor
from clubb_jax.src.CLUBB_core.constants_clubb import rc_tol

CHI, ETA, W, NCN = 0, 1, 2, 3
PDF_DIM = 6
HM = [4, 5]
NG, NZT = 3, 4


def _pd_corr(seed):
    """A symmetric positive-definite correlation matrix (unit diagonal)."""
    rng = np.random.default_rng(seed)
    B = rng.standard_normal((PDF_DIM, PDF_DIM))
    cov = B @ B.T + PDF_DIM * np.eye(PDF_DIM)
    d = np.sqrt(np.diag(cov))
    return cov / np.outer(d, d)


def _ref_adjust(cc, cb, iiPDF_type, l_adg):
    """Literal NumPy reproduction of the in-place matrix adjustments."""
    cc = cc.copy(); cb = cb.copy()
    if l_adg and iiPDF_type in (1, 2, 7):
        for (r, c) in ((W, CHI), (W, ETA)):
            cc[r, c] = 0.0; cb[r, c] = 0.0
    cb[NCN, CHI] = cc[NCN, CHI]
    cb[NCN, ETA] = cc[NCN, ETA]
    for j in HM:
        cc[j, ETA] = cc[ETA, CHI] * cc[j, CHI]
        cb[j, ETA] = cb[ETA, CHI] * cb[j, CHI]
    def _symm(M):
        L = np.tril(M)
        return L + L.T - np.diag(np.diag(M))
    return _symm(cc), _symm(cb)


def _inputs(seed):
    rng = np.random.default_rng(seed)
    cc = _pd_corr(seed); cb = _pd_corr(seed + 100)
    rc_1 = rng.choice([0.0, 2e-3], (NG, NZT))
    rc_2 = rng.choice([0.0, 2e-3], (NG, NZT))
    return cc, cb, rc_1, rc_2


def test_reconstruction_and_structure():
    cc, cb, rc_1, rc_2 = _inputs(1)
    corr_1, corr_2, chol_1, chol_2 = calc_corr_norm_and_cholesky_factor(
        cc, cb, rc_1, rc_2, 6, CHI, ETA, W, NCN, True)
    for corr, chol in ((np.asarray(corr_1), np.asarray(chol_1)), (np.asarray(corr_2), np.asarray(chol_2))):
        assert np.max(np.abs(corr - np.swapaxes(corr, -1, -2))) < 1e-14, "corr not symmetric"
        d = corr[:, :, np.arange(PDF_DIM), np.arange(PDF_DIM)]
        assert np.max(np.abs(d - 1.0)) < 1e-14, "corr diag != 1"
        iu = np.triu_indices(PDF_DIM, 1)
        assert np.max(np.abs(chol[:, :, iu[0], iu[1]])) < 1e-14, "Cholesky not lower-triangular"
        recon = np.einsum('ikab,ikcb->ikac', chol, chol)   # L Lᵀ
        assert np.max(np.abs(recon - corr)) < 1e-10, "L Lᵀ != assigned corr"
    print("  reconstruction L Lᵀ = corr + symmetric/unit-diag/lower-tri structure  PASS")


def test_adjustments_and_selection():
    cc, cb, rc_1, rc_2 = _inputs(2)
    iiPDF_type, l_adg = 1, True   # ADG1 -> w-chi/w-eta zeroed
    corr_1, corr_2, _, _ = calc_corr_norm_and_cholesky_factor(
        cc, cb, rc_1, rc_2, iiPDF_type, CHI, ETA, W, NCN, l_adg)
    cc_s, cb_s = _ref_adjust(cc, cb, iiPDF_type, l_adg)
    # Selection: rc>rc_tol -> cloud, else below.
    sel1 = (rc_1 > rc_tol)[:, :, None, None]
    ref_1 = np.where(sel1, cc_s, cb_s)
    sel2 = (rc_2 > rc_tol)[:, :, None, None]
    ref_2 = np.where(sel2, cc_s, cb_s)
    assert np.max(np.abs(np.asarray(corr_1) - ref_1)) < 1e-14, "corr_1 selection/adjust mismatch"
    assert np.max(np.abs(np.asarray(corr_2) - ref_2)) < 1e-14, "corr_2 selection/adjust mismatch"
    # ADG zeroing actually applied.
    assert np.all(np.asarray(corr_1)[:, :, W, CHI] == 0.0) and np.all(np.asarray(corr_1)[:, :, W, ETA] == 0.0)
    print("  matrix adjustments (ADG/Ncn/eta-hm) + rc-selection vs literal NumPy  PASS")


def test_differentiable():
    cc, cb, rc_1, rc_2 = _inputs(3)
    def loss(x):
        c1, c2, l1, l2 = calc_corr_norm_and_cholesky_factor(x, cb, rc_1, rc_2, 6, CHI, ETA, W, NCN, True)
        return jnp.sum(c1 ** 2) + jnp.sum(l1 ** 2) + jnp.sum(c2 ** 2) + jnp.sum(l2 ** 2)
    g = np.asarray(jax.grad(loss)(jnp.asarray(cc)))
    assert np.isfinite(g).all(), "non-finite grad through calc_corr_norm_and_cholesky_factor"
    print(f"  jax.grad through calc_corr_norm_and_cholesky_factor: finite ({g.size} entries)  PASS")


def main():
    print("test_calc_corr_norm_cholesky:")
    for t in (test_reconstruction_and_structure, test_adjustments_and_selection, test_differentiable):
        t()
    print("All calc_corr_norm_and_cholesky_factor checks PASSED")


if __name__ == "__main__":
    main()
