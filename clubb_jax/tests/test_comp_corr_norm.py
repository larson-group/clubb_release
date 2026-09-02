#!/usr/bin/env python3
"""test_comp_corr_norm.py — validate the JAX comp_corr_norm orchestration.

comp_corr_norm (setup_clubb_pdf_params.F90:1273) assembles the lower-triangular, then symmetrized, normal-space
PDF correlation arrays from the component_corr_* routines + calc_corr_w_hm_n. Oracle here is a combination of:
  1. Structural invariants the Fortran guarantees: the output is symmetric with a unit diagonal.
  2. Spot-checks of every distinct assembly rule against the building-block functions called directly
     (chi_eta, w_x, x_hm cloud-twice for Ncn, eta_hm product, w_hm_n_ip from calc_corr_w_hm_n).
  3. A finite jax.grad through the assembled array.
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

from clubb_jax.src.CLUBB_core.setup_clubb_pdf_params import (
    comp_corr_norm, component_corr_chi_eta, component_corr_w_x,
    component_corr_x_hm_n_ip, component_corr_eta_hm_n_ip,
    component_corr_w_hm_n_ip, calc_corr_w_hm_n)

# KK PDF layout: chi=0, eta=1, w=2, Ncn=3, rr=4, Nr=5
CHI, ETA, W, NCN = 0, 1, 2, 3
PDF_DIM, HM_DIM = 6, 2
PDF2HM = np.array([-1, -1, -1, -1, 0, 1])      # rr->hm 0, Nr->hm 1
HM_TOL = np.array([1.0e-10, 1.0e-3])
NCN_TOL = 1.0e-8
NG, NZT = 2, 5


def _inputs(seed):
    rng = np.random.default_rng(seed)
    mu_x_1 = rng.uniform(0.5, 3.0, (NG, NZT, PDF_DIM))
    mu_x_2 = rng.uniform(0.5, 3.0, (NG, NZT, PDF_DIM))
    mu_x_1[:, :, W] = rng.uniform(-1, 1, (NG, NZT))   # w mean can be negative
    mu_x_2[:, :, W] = rng.uniform(-1, 1, (NG, NZT))
    sigma_x_1 = rng.uniform(0.1, 1.0, (NG, NZT, PDF_DIM))
    sigma_x_2 = rng.uniform(0.1, 1.0, (NG, NZT, PDF_DIM))
    sigma_x_1_n = rng.uniform(0.1, 1.0, (NG, NZT, PDF_DIM))
    sigma_x_2_n = rng.uniform(0.1, 1.0, (NG, NZT, PDF_DIM))
    wm_zt = rng.uniform(-1, 1, (NG, NZT))
    rc_1 = rng.choice([0.0, 2e-3], (NG, NZT))         # straddle rc_tol
    rc_2 = rng.choice([0.0, 2e-3], (NG, NZT))
    mixt_frac = rng.uniform(0.3, 0.7, (NG, NZT))
    pf1 = rng.uniform(0.4, 1.0, (NG, NZT))
    pf2 = rng.uniform(0.4, 1.0, (NG, NZT))
    wpNcnp = rng.uniform(-1, 1, (NG, NZT))
    wphmp = rng.uniform(-1e-3, 1e-3, (NG, NZT, HM_DIM))
    # Prescribed correlation arrays (lower triangular meaningful); symmetric-ish random in (-0.8,0.8).
    cc = rng.uniform(-0.8, 0.8, (PDF_DIM, PDF_DIM))
    cb = rng.uniform(-0.8, 0.8, (PDF_DIM, PDF_DIM))
    return dict(mu_x_1=mu_x_1, mu_x_2=mu_x_2, sigma_x_1=sigma_x_1, sigma_x_2=sigma_x_2,
                sigma_x_1_n=sigma_x_1_n, sigma_x_2_n=sigma_x_2_n, wm_zt=wm_zt, rc_1=rc_1, rc_2=rc_2,
                mixt_frac=mixt_frac, precip_frac_1=pf1, precip_frac_2=pf2, wpNcnp_zt=wpNcnp,
                wphydrometp_zt=wphmp, corr_array_n_cloud=cc, corr_array_n_below=cb,
                iiPDF_chi=CHI, iiPDF_eta=ETA, iiPDF_w=W, iiPDF_Ncn=NCN, pdf2hydromet=PDF2HM,
                hydromet_tol=HM_TOL, Ncn_tol_val=NCN_TOL)


def _call(inp, iiPDF_type, l_calc_w_corr, l_fix):
    return comp_corr_norm(iiPDF_type=iiPDF_type, l_calc_w_corr=l_calc_w_corr,
                          l_fix_w_chi_eta_correlations=l_fix, **inp)


def test_structure():
    inp = _inputs(1)
    A1, A2 = _call(inp, iiPDF_type=6, l_calc_w_corr=True, l_fix=True)
    for A in (np.asarray(A1), np.asarray(A2)):
        assert A.shape == (NG, NZT, PDF_DIM, PDF_DIM)
        assert np.max(np.abs(A - np.swapaxes(A, -1, -2))) < 1e-14, "not symmetric"
        diag = A[:, :, np.arange(PDF_DIM), np.arange(PDF_DIM)]
        assert np.max(np.abs(diag - 1.0)) < 1e-14, "diagonal not 1"
        assert np.max(np.abs(A)) <= 1.0 + 1e-12, "correlation magnitude > 1"
    print("  structure: symmetric, unit diagonal, |corr|<=1  PASS")


def test_assembly_rules():
    inp = _inputs(2)
    A1, A2 = _call(inp, iiPDF_type=6, l_calc_w_corr=True, l_fix=True)   # non-ADG type -> w_x cloud/below
    A1 = np.asarray(A1); A2 = np.asarray(A2)
    cc, cb = inp['corr_array_n_cloud'], inp['corr_array_n_below']
    rc_1, rc_2 = inp['rc_1'], inp['rc_2']

    # (eta, chi): component_corr_chi_eta with l_limit=True
    e1, e2 = component_corr_chi_eta(rc_1, rc_2, cc[ETA, CHI], cb[ETA, CHI], True)
    assert np.max(np.abs(A1[:, :, ETA, CHI] - np.asarray(e1))) < 1e-14
    assert np.max(np.abs(A2[:, :, ETA, CHI] - np.asarray(e2))) < 1e-14

    # (w, chi): component_corr_w_x (non-ADG -> cloud/below)
    w1, w2 = component_corr_w_x(rc_1, rc_2, cc[W, CHI], cb[W, CHI], 6, True)
    assert np.max(np.abs(A1[:, :, W, CHI] - np.asarray(w1))) < 1e-14
    assert np.max(np.abs(A2[:, :, W, CHI] - np.asarray(w2))) < 1e-14

    # (Ncn, chi): cloud value used twice -> constant cc[Ncn,chi] everywhere
    assert np.max(np.abs(A1[:, :, NCN, CHI] - cc[NCN, CHI])) < 1e-14
    assert np.max(np.abs(A2[:, :, NCN, CHI] - cc[NCN, CHI])) < 1e-14

    # (hm, eta): product of (eta,chi) and (hm,chi)
    for j in (4, 5):
        assert np.max(np.abs(A1[:, :, j, ETA] - A1[:, :, ETA, CHI] * A1[:, :, j, CHI])) < 1e-14
        assert np.max(np.abs(A2[:, :, j, ETA] - A2[:, :, ETA, CHI] * A2[:, :, j, CHI])) < 1e-14

    # (Ncn, w): component_corr_w_hm_n_ip of calc_corr_w_hm_n diagnosis (l_calc_w_corr=True -> passthrough)
    cwn1, cwn2 = calc_corr_w_hm_n(
        inp['wm_zt'], inp['wpNcnp_zt'], inp['mu_x_1'][:, :, W], inp['mu_x_2'][:, :, W],
        inp['mu_x_1'][:, :, NCN], inp['mu_x_2'][:, :, NCN], inp['sigma_x_1'][:, :, W], inp['sigma_x_2'][:, :, W],
        inp['sigma_x_1'][:, :, NCN], inp['sigma_x_2'][:, :, NCN], inp['sigma_x_1_n'][:, :, NCN],
        inp['sigma_x_2_n'][:, :, NCN], inp['mixt_frac'], np.ones((NG, NZT)), np.ones((NG, NZT)), NCN_TOL)
    p1, p2 = component_corr_w_hm_n_ip(cwn1, rc_1, cwn2, rc_2, cc[NCN, W], cb[NCN, W], True)
    assert np.max(np.abs(A1[:, :, NCN, W] - np.asarray(p1))) < 1e-13, "(Ncn,w) mismatch"
    assert np.max(np.abs(A2[:, :, NCN, W] - np.asarray(p2))) < 1e-13
    print("  assembly rules: chi_eta / w_x / Ncn-chi-twice / eta_hm product / w_hm from calc_corr_w_hm_n  PASS")


def test_adg_zeros_w_corr():
    inp = _inputs(3)
    A1, A2 = _call(inp, iiPDF_type=1, l_calc_w_corr=True, l_fix=True)   # ADG1 -> corr(w,chi)=corr(w,eta)=0
    A1 = np.asarray(A1); A2 = np.asarray(A2)
    assert np.all(A1[:, :, W, CHI] == 0.0) and np.all(A1[:, :, W, ETA] == 0.0), "ADG1 should zero w-chi/w-eta"
    assert np.all(A2[:, :, W, CHI] == 0.0) and np.all(A2[:, :, W, ETA] == 0.0)
    print("  ADG1 PDF type: corr(w,chi)=corr(w,eta)=0  PASS")


def test_no_w_corr_uses_prescribed():
    inp = _inputs(4)
    # l_calc_w_corr=False -> w-hm correlations come from cloud/below prescribed values, not the diagnosis.
    A1, _ = _call(inp, iiPDF_type=6, l_calc_w_corr=False, l_fix=True)
    A1 = np.asarray(A1)
    cc, cb, rc_1 = inp['corr_array_n_cloud'], inp['corr_array_n_below'], inp['rc_1']
    expect = np.where(rc_1 > 1e-6, cc[NCN, W], cb[NCN, W])
    assert np.max(np.abs(A1[:, :, NCN, W] - expect)) < 1e-14, "no-w-corr (Ncn,w) should be cloud/below"
    print("  l_calc_w_corr=False: w-correlations use prescribed cloud/below values  PASS")


def test_differentiable():
    inp = _inputs(5)
    def loss(wphmp):
        i2 = dict(inp); i2['wphydrometp_zt'] = wphmp
        A1, A2 = comp_corr_norm(iiPDF_type=6, l_calc_w_corr=True, l_fix_w_chi_eta_correlations=True, **i2)
        return jnp.sum(A1 ** 2) + jnp.sum(A2 ** 2)
    g = np.asarray(jax.grad(loss)(jnp.asarray(inp['wphydrometp_zt'])))
    assert np.isfinite(g).all(), "non-finite grad through comp_corr_norm"
    print(f"  jax.grad through comp_corr_norm: finite ({g.size} entries)  PASS")


def main():
    print("test_comp_corr_norm:")
    for t in (test_structure, test_assembly_rules, test_adg_zeros_w_corr,
              test_no_w_corr_uses_prescribed, test_differentiable):
        t()
    print("All comp_corr_norm checks PASSED")


if __name__ == "__main__":
    main()
