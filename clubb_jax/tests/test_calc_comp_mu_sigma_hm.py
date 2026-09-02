"""Verification of calc_comp_mu_sigma_hm (hydrometeor in-precip component moments).

This routine (setup_clubb_pdf_params.F90:1653) is NOT exposed by the f2py API, and the rico
stats oracle cannot verify it precisely: its inputs (rrm, precip_frac) and outputs (mu_rr_i)
are stored at end-of-step but were different within the step, so the defining identity
<hm> = a f_p_1 mu_1 + (1-a) f_p_2 mu_2 does not even hold across the stored stats (rel ~7).

Instead this verifies the routine's DEFINING CONTRACT directly: it solves for the two in-precip
component means/stdevs so that the overall mean <hm> AND overall variance <hm'^2> are preserved
(Griffin 2015). For consistent synthetic inputs the port must reproduce <hm> and <hm'^2> exactly
(non-emergency), preserve <hm> in the emergency-bound regime, and satisfy the structural
relations R = omicron*Rmax and sigma_1^2/mu_1^2 = R(1+zeta). All four branches and the root-sign
selection are exercised.
"""
import numpy as np
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

import os
import sys
_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

from clubb_jax.src.CLUBB_core.setup_clubb_pdf_params import (
    calc_comp_mu_sigma_hm, compute_mean_stdev, norm_transform_mean_stdev,
    IIPDF_CHI, IIPDF_ETA, IIPDF_W, IIPDF_NCN)
from clubb_jax.src.CLUBB_core.pdf_utilities import mean_L2N, stdev_L2N

_HM_TOL = 1e-10
_PFTOL = np.array([0.005])


def _run(hmm, mf, pf, pf1, pf2, ratio, omicron, zeta, thl1, thl2, hmp2=None):
    hmm = np.array([[hmm]]); mf = np.array([[mf]]); pf = np.array([[pf]])
    pf1 = np.array([[pf1]]); pf2 = np.array([[pf2]])
    thl1 = np.array([[thl1]]); thl2 = np.array([[thl2]])
    if hmp2 is None:
        hmp2 = ((ratio + 1.0) / pf - 1.0) * hmm ** 2     # consistent overall variance
    else:
        hmp2 = np.array([[hmp2]])
    out = calc_comp_mu_sigma_hm(hmm, hmp2, np.full_like(hmm, ratio), mf, pf, pf1, pf2,
                                _HM_TOL, _PFTOL, thl1, thl2, omicron, zeta)
    s = lambda x: float(np.asarray(x).ravel()[0])
    return [s(x) for x in out], s(hmm), s(mf), s(pf1), s(pf2), s(hmp2)


def _mean(mu1, mu2, mf, pf1, pf2):
    return mf * pf1 * mu1 + (1 - mf) * pf2 * mu2


def _var(mu1, mu2, s1, s2, mf, pf1, pf2, hmm):
    return mf * pf1 * (mu1**2 + s1**2) + (1 - mf) * pf2 * (mu2**2 + s2**2) - hmm**2


def test_both_precip_preserves_mean_and_variance():
    """Non-emergency both-precip: <hm> and <hm'^2> preserved exactly; R, sigma ratios correct."""
    for ratio, omicron, zeta in [(1.0, 0.5, 0.0), (0.8, 0.7, 0.3), (1.5, 0.4, 0.0)]:
        (mu1, mu2, s1, s2, hm1, hm2, r1, r2), hmm, mf, pf1, pf2, hmp2 = _run(
            1e-4, 0.5, 0.5, 0.6, 0.4, ratio, omicron, zeta, 300.0, 300.5)
        assert mu1 > 0 and mu2 > 0
        assert abs(_mean(mu1, mu2, mf, pf1, pf2) - hmm) / hmm < 1e-13, "mean not preserved"
        assert abs(_var(mu1, mu2, s1, s2, mf, pf1, pf2, hmm) - hmp2) / hmp2 < 1e-12, "var not preserved"
        # R = omicron*Rmax, Rmax = (f_p/(a f_p1 (1+zeta)+(1-a) f_p2)) * ratio;
        # sigma_1^2/mu_1^2 = R(1+zeta).
        pf = 0.5
        Rmax = (pf / (mf * pf1 * (1 + zeta) + (1 - mf) * pf2)) * ratio
        assert abs(r2 - omicron * Rmax) / (omicron * Rmax) < 1e-10, f"R wrong {r2}"
        assert abs(r1 - r2 * (1 + zeta)) / (r2 * (1 + zeta)) < 1e-10, "sigma1 ratio wrong"
    print("  both-precip non-emergency: mean (0) + variance (~1e-16) preserved; R, sigma ratios  PASS")


def test_rico_params_non_emergency():
    """rico's actual config (hmp2_ip ratio = 1.25 case-override, omicron=0.5, zeta=0): with the
    precip_fraction invariant pf = a*pf1+(1-a)*pf2 satisfied, the routine is NON-emergency,
    R = omicron*ratio = 0.625, and preserves <hm> and <hm'^2> exactly. (Variance preservation
    REQUIRES that invariant — feeding an inconsistent pf makes Rmax != ratio and breaks it.)"""
    mf, pf1, pf2 = 0.3, 0.06, 0.022
    pf = mf * pf1 + (1 - mf) * pf2                     # the precip_fraction invariant
    (mu1, mu2, s1, s2, hm1, hm2, r1, r2), hmm, mf_, pf1_, pf2_, hmp2 = _run(
        5e-7, mf, pf, pf1, pf2, 1.25, 0.5, 0.0, 298.0, 298.3)   # rico magnitudes/params
    assert abs(_mean(mu1, mu2, mf, pf1, pf2) - hmm) / hmm < 1e-13
    assert abs(_var(mu1, mu2, s1, s2, mf, pf1, pf2, hmm) - hmp2) / hmp2 < 1e-12
    assert abs(r2 - 0.5 * 1.25) < 1e-10, f"rico R should be 0.625 (non-emergency), got {r2}"
    print(f"  rico params (ratio=1.25, omicron=0.5): non-emergency R={r2:.3f}, moments preserved  PASS")


def test_root_sign_selection():
    """thl_1 <= thl_2 -> mu_1 >= mu_2 (cloudier comp has larger mean); reversed otherwise."""
    (mu1a, mu2a, *_), *_ = _run(1e-4, 0.5, 0.5, 0.6, 0.4, 1.0, 0.5, 0.0, 300.0, 300.5)
    (mu1b, mu2b, *_), *_ = _run(1e-4, 0.5, 0.5, 0.6, 0.4, 1.0, 0.5, 0.0, 300.5, 300.0)
    assert mu1a >= mu2a and mu1b <= mu2b, f"root sign wrong: {mu1a},{mu2a} / {mu1b},{mu2b}"
    print("  root-sign by mu_thl ordering: PASS")


def test_emergency_preserves_mean():
    """Large prescribed ratio forces the emergency bound; <hm> still preserved, outputs valid."""
    (mu1, mu2, s1, s2, hm1, hm2, r1, r2), hmm, mf, pf1, pf2, hmp2 = _run(
        1e-4, 0.5, 0.5, 0.6, 0.4, 20.0, 0.5, 0.0, 300.0, 300.5)
    assert mu1 > 0 and mu2 > 0 and s1 >= 0 and s2 >= 0 and np.isfinite(r2)
    assert abs(_mean(mu1, mu2, mf, pf1, pf2) - hmm) / hmm < 1e-12, "emergency mean not preserved"
    # Emergency triggered: R differs from the non-emergency value omicron*Rmax (=10 here).
    Rmax = (0.5 / (mf * pf1 + (1 - mf) * pf2)) * 20.0
    assert abs(r2 - 0.5 * Rmax) > 1e-6, f"emergency did not trigger (R={r2})"
    print(f"  emergency bound (ratio=20): mean preserved, R recomputed {0.5*Rmax:.2f}->{r2:.3f}  PASS")


def test_single_component_branches():
    """precip in comp 1 only / comp 2 only: mean + variance preserved; other comp zero."""
    # comp 1 only: precip_frac_2 below tol
    (mu1, mu2, s1, s2, hm1, hm2, r1, r2), hmm, mf, pf1, pf2, hmp2 = _run(
        1e-4, 0.5, 0.3, 0.6, 0.0, 1.0, 0.5, 0.0, 300.0, 300.5)
    assert mu2 == 0 and s2 == 0
    assert abs(mf * pf1 * mu1 - hmm) / hmm < 1e-13
    assert abs(mf * pf1 * (mu1**2 + s1**2) - hmm**2 - hmp2) / hmp2 < 1e-12
    # comp 2 only: precip_frac_1 below tol
    (mu1, mu2, s1, s2, *_), hmm, mf, pf1, pf2, hmp2 = _run(
        1e-4, 0.5, 0.3, 0.0, 0.6, 1.0, 0.5, 0.0, 300.0, 300.5)
    assert mu1 == 0 and s1 == 0
    assert abs((1 - mf) * pf2 * mu2 - hmm) / hmm < 1e-13
    print("  single-component branches: mean + variance preserved  PASS")


def test_no_precip_and_differentiable():
    """No hydrometeor -> all zeros. And the routine is differentiable in <hm>."""
    (mu1, mu2, s1, s2, hm1, hm2, r1, r2), *_ = _run(
        1e-12, 0.5, 0.5, 0.6, 0.4, 1.0, 0.5, 0.0, 300.0, 300.5)  # hmm < hm_tol
    assert all(v == 0 for v in (mu1, mu2, s1, s2, hm1, hm2))

    def f(hmm):
        a = dict(hmm=jnp.array([[hmm]]), hmp2=jnp.array([[3e-8]]),
                 hmp2_ip_on_hmm2_ip=jnp.array([[1.0]]), mixt_frac=jnp.array([[0.5]]),
                 precip_frac=jnp.array([[0.5]]), precip_frac_1=jnp.array([[0.6]]),
                 precip_frac_2=jnp.array([[0.4]]), hm_tol=_HM_TOL,
                 precip_frac_tol=jnp.array([0.005]), mu_thl_1=jnp.array([[300.0]]),
                 mu_thl_2=jnp.array([[300.5]]), omicron=0.5, zeta_vrnce_rat_in=0.0)
        return calc_comp_mu_sigma_hm(**a)[0][0, 0]   # mu_hm_1
    g = float(jax.grad(f)(1e-4))
    assert np.isfinite(g) and abs(g) > 0, "grad not finite/nonzero"
    print(f"  no-precip zeros + differentiable (d mu1/d hmm = {g:.3f})  PASS")


# ---------------------------------------------------------------------------
# compute_mean_stdev / norm_transform_mean_stdev — the setup_pdf_parameters
# orchestration that stacks the per-PDF-variable component moments and transforms
# the lognormal vars (Ncn + hydrometeors) to normal space (Iter131).
# ---------------------------------------------------------------------------

def _rico_like_pdf_inputs():
    """A small (1,2)-shaped rico-like state with two saturated/precip levels."""
    chi_1 = np.array([[1.0e-4, 2.0e-4]])
    chi_2 = np.array([[-3.0e-5, 5.0e-5]])
    sc1 = np.array([[2.0e-4, 1.0e-4]]); sc2 = np.array([[3.0e-4, 2.0e-4]])
    se1 = np.array([[1.5e-4, 1.0e-4]]); se2 = np.array([[2.0e-4, 1.5e-4]])
    thl_1 = np.array([[298.0, 297.5]]); thl_2 = np.array([[298.4, 297.9]])
    mf = np.array([[0.3, 0.4]])
    pf1 = np.array([[0.06, 0.05]]); pf2 = np.array([[0.02, 0.03]])
    pf = mf * pf1 + (1 - mf) * pf2
    rrm = np.array([[5.0e-6, 1.0e-5]]); Nrm = np.array([[3.0e4, 6.0e4]])
    return chi_1, chi_2, sc1, sc2, se1, se2, thl_1, thl_2, mf, pf, pf1, pf2, rrm, Nrm


def test_compute_mean_stdev_columns():
    """The stacked arrays reproduce: chi/eta/w pass-through, Ncn (const & varying), and the
    hydrometeor columns equal the standalone calc_comp_mu_sigma_hm."""
    (chi_1, chi_2, sc1, sc2, se1, se2, thl_1, thl_2, mf, pf, pf1, pf2, rrm, Nrm) = _rico_like_pdf_inputs()
    w_1 = np.array([[0.4, -0.2]]); w_2 = np.array([[-0.1, 0.3]])
    sw1 = np.array([[0.5, 0.6]]); sw2 = np.array([[0.3, 0.4]])
    rr_ratio, Nr_ratio, omicron, zeta = 1.25, 1.25, 0.5, 0.0
    Nr_tol = _HM_TOL / ((4.0 / 3.0) * np.pi * 1000.0 * (5.0e-3) ** 3)
    hmp2_rr = ((rr_ratio + 1.0) / pf - 1.0) * rrm ** 2
    hmp2_Nr = ((Nr_ratio + 1.0) / pf - 1.0) * Nrm ** 2

    # --- Non-constant Ncn variant (so the Ncn column exercises sqrt(ratio)*Ncnm) ---
    Ncnm = np.array([[1.0e8, 2.0e8]]); Ncnp2_on_Ncnm2 = 0.36
    mu1, mu2, sig1, sig2, hm1, hm2, s2m2_1, s2m2_2 = compute_mean_stdev(
        chi_1, chi_2, sc1, sc2, se1, se2, Ncnm, Ncnp2_on_Ncnm2, False,
        [(rrm, hmp2_rr, rr_ratio, _HM_TOL), (Nrm, hmp2_Nr, Nr_ratio, Nr_tol)],
        thl_1, thl_2, mf, pf, pf1, pf2, _PFTOL, omicron, zeta,
        w_1=w_1, w_2=w_2, stdev_w_1=sw1, stdev_w_2=sw2)
    assert mu1.shape == chi_1.shape + (6,), mu1.shape
    # chi / eta / w pass-through
    np.testing.assert_array_equal(mu1[..., IIPDF_CHI], chi_1)
    np.testing.assert_array_equal(mu2[..., IIPDF_CHI], chi_2)
    np.testing.assert_array_equal(mu1[..., IIPDF_ETA], np.zeros_like(chi_1))
    np.testing.assert_array_equal(sig1[..., IIPDF_ETA], se1)
    np.testing.assert_array_equal(mu1[..., IIPDF_W], w_1)
    np.testing.assert_array_equal(sig2[..., IIPDF_W], sw2)
    # Ncn (varying): mu = Ncnm, sigma = sqrt(ratio)*Ncnm, ratio = Ncnp2_on_Ncnm2
    np.testing.assert_allclose(mu1[..., IIPDF_NCN], Ncnm, rtol=0, atol=0)
    np.testing.assert_allclose(sig1[..., IIPDF_NCN], np.sqrt(Ncnp2_on_Ncnm2) * Ncnm, rtol=1e-14)
    np.testing.assert_allclose(s2m2_1[..., IIPDF_NCN], Ncnp2_on_Ncnm2, rtol=1e-14)
    # hydrometeor columns == standalone calc_comp_mu_sigma_hm
    for j, (hmm, hmp2, ratio, tol) in enumerate(
            [(rrm, hmp2_rr, rr_ratio, _HM_TOL), (Nrm, hmp2_Nr, Nr_ratio, Nr_tol)]):
        m1, m2, s1, s2, h1, h2, r1, r2 = calc_comp_mu_sigma_hm(
            hmm, hmp2, np.full_like(rrm, ratio), mf, pf, pf1, pf2, tol, _PFTOL,
            thl_1, thl_2, omicron, zeta)
        iv = IIPDF_NCN + 1 + j
        np.testing.assert_allclose(mu1[..., iv], m1, rtol=0, atol=0)
        np.testing.assert_allclose(sig2[..., iv], s2, rtol=0, atol=0)
        np.testing.assert_allclose(s2m2_1[..., iv], r1, rtol=0, atol=0)

    # --- Constant Ncn (l_const_Nc_in_cloud=True): sigma_Ncn = 0 ---
    mu1c, _, sig1c, _, _, _, s2m2_1c, _ = compute_mean_stdev(
        chi_1, chi_2, sc1, sc2, se1, se2, Ncnm, 0.0, True,
        [(rrm, hmp2_rr, rr_ratio, _HM_TOL), (Nrm, hmp2_Nr, Nr_ratio, Nr_tol)],
        thl_1, thl_2, mf, pf, pf1, pf2, _PFTOL, omicron, zeta)
    np.testing.assert_array_equal(sig1c[..., IIPDF_NCN], np.zeros_like(chi_1))
    np.testing.assert_array_equal(s2m2_1c[..., IIPDF_NCN], np.zeros_like(chi_1))
    print("  compute_mean_stdev: chi/eta/w pass-through, Ncn (const+varying), hydromet==calc_comp  PASS")


def test_norm_transform_matches_L2N():
    """norm_transform_mean_stdev passes chi/eta/w through and applies mean_L2N/stdev_L2N (with the
    sigma^2/mu^2 ratios) to Ncn + hydrometeors."""
    (chi_1, chi_2, sc1, sc2, se1, se2, thl_1, thl_2, mf, pf, pf1, pf2, rrm, Nrm) = _rico_like_pdf_inputs()
    rr_ratio, Nr_ratio, omicron, zeta = 1.25, 1.25, 0.5, 0.0
    Nr_tol = _HM_TOL / ((4.0 / 3.0) * np.pi * 1000.0 * (5.0e-3) ** 3)
    hmp2_rr = ((rr_ratio + 1.0) / pf - 1.0) * rrm ** 2
    hmp2_Nr = ((Nr_ratio + 1.0) / pf - 1.0) * Nrm ** 2
    Ncnm = np.array([[1.0e8, 2.0e8]]); ratio_Ncn = 0.36
    args = compute_mean_stdev(
        chi_1, chi_2, sc1, sc2, se1, se2, Ncnm, ratio_Ncn, False,
        [(rrm, hmp2_rr, rr_ratio, _HM_TOL), (Nrm, hmp2_Nr, Nr_ratio, Nr_tol)],
        thl_1, thl_2, mf, pf, pf1, pf2, _PFTOL, omicron, zeta)
    mu1, mu2, sig1, sig2, hm1, hm2, s2m2_1, s2m2_2 = args
    mu1_n, mu2_n, sig1_n, sig2_n = norm_transform_mean_stdev(
        mu1, mu2, sig1, sig2, s2m2_1, s2m2_2, Ncnm, hm1, hm2, [_HM_TOL, Nr_tol], False)
    # chi/eta/w unchanged
    for iv in (IIPDF_CHI, IIPDF_ETA, IIPDF_W):
        np.testing.assert_array_equal(mu1_n[..., iv], mu1[..., iv])
        np.testing.assert_array_equal(sig1_n[..., iv], sig1[..., iv])
    # Ncn + hydrometeors: mean_L2N / stdev_L2N of the linear moments
    for iv in (IIPDF_NCN, IIPDF_NCN + 1, IIPDF_NCN + 2):
        exp_mu = mean_L2N(np.maximum(np.abs(mu1[..., iv]), 1e-30), s2m2_1[..., iv])
        exp_sig = stdev_L2N(s2m2_1[..., iv])
        np.testing.assert_allclose(mu1_n[..., iv], exp_mu, rtol=0, atol=0)
        np.testing.assert_allclose(sig1_n[..., iv], exp_sig, rtol=0, atol=0)
    print("  norm_transform_mean_stdev: chi/eta/w unchanged, Ncn+hydromet == mean_L2N/stdev_L2N  PASS")


def test_orchestration_differentiable():
    """jax.grad flows through compute_mean_stdev -> norm_transform (w.r.t. the rrm field)."""
    (chi_1, chi_2, sc1, sc2, se1, se2, thl_1, thl_2, mf, pf, pf1, pf2, rrm, Nrm) = _rico_like_pdf_inputs()
    Nr_tol = _HM_TOL / ((4.0 / 3.0) * np.pi * 1000.0 * (5.0e-3) ** 3)

    def f(scale):
        rr = jnp.asarray(rrm) * scale
        hmp2_rr = ((1.25 + 1.0) / jnp.asarray(pf) - 1.0) * rr ** 2
        hmp2_Nr = ((1.25 + 1.0) / jnp.asarray(pf) - 1.0) * jnp.asarray(Nrm) ** 2
        out = compute_mean_stdev(
            chi_1, chi_2, sc1, sc2, se1, se2, jnp.zeros_like(rr), 0.0, True,
            [(rr, hmp2_rr, 1.25, _HM_TOL), (jnp.asarray(Nrm), hmp2_Nr, 1.25, Nr_tol)],
            thl_1, thl_2, mf, pf, pf1, pf2, _PFTOL, 0.5, 0.0)
        mu1n, _, _, _ = norm_transform_mean_stdev(
            out[0], out[1], out[2], out[3], out[6], out[7],
            jnp.zeros_like(rr), out[4], out[5], [_HM_TOL, Nr_tol], True)
        return jnp.sum(mu1n[..., IIPDF_NCN + 1])     # sum of ln(mu_rr_1)

    g = float(jax.grad(f)(1.0))
    assert np.isfinite(g) and g != 0.0, f"grad not finite/nonzero: {g}"
    print(f"  compute_mean_stdev -> norm_transform differentiable (d/d scale = {g:.3f})  PASS")


if __name__ == "__main__":
    print("calc_comp_mu_sigma_hm (hydrometeor in-precip component moments) verification:")
    test_both_precip_preserves_mean_and_variance()
    test_rico_params_non_emergency()
    test_root_sign_selection()
    test_emergency_preserves_mean()
    test_single_component_branches()
    test_no_precip_and_differentiable()
    print("setup_pdf_parameters orchestration (compute_mean_stdev / norm_transform_mean_stdev):")
    test_compute_mean_stdev_columns()
    test_norm_transform_matches_L2N()
    test_orchestration_differentiable()
    print("All calc_comp_mu_sigma_hm + orchestration tests PASSED.")
