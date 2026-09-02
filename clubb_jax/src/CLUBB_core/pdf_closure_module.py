"""JAX port of pdf_closure_module.F90 - PDF moment-integral closures.

Mirrors clubb_release/src/CLUBB_core/pdf_closure_module.F90. The closure
orchestration lives here as `pdf_closure_driver`; advance_clubb_core.py calls it
rather than inlining it. Alongside the driver this module holds the standalone
analytic moment-integral routines that close the higher-order moments from the
two-component PDF parameters, plus the cloudy-updraft diagnostic:

  calc_wp2xp_pdf / calc_wpxp2_pdf / calc_wp2xp2_pdf / calc_wp4_pdf / calc_wpxpyp_pdf
      — <w'^2 x'>, <w'x'^2>, <w'^2 x'^2>, <w'^4>, <w'x'y'> from the binormal/trinormal
        PDF integrals used directly by the live pdf_closure body.
  calc_w_up_in_cloud — mean cloudy updraft/downdraft vertical velocity (aerosol activation)

Porting deviations:
- `pdf_closure` is not yet a full parity port. The hydrometeor/liquid-ice
  loading branch from Fortran is intentionally omitted, and hydrometeor moment
  inputs are deleted at the top of the routine.
- Fortran stats-gated shortcuts are not duplicated here; the JAX routine
  computes moment fields unconditionally and stats updates are handled outside
  this low-level function.
- Fortran intent(out)/intent(inout) arguments are returned through dictionaries
  or functional derived-type replacements.
- Unsupported `iiPDF_type` values raise `NotImplementedError` rather than
  writing to stderr/error-stopping.
- The Fortran host-model API warning is not directly applicable to this JAX
  internal function, but the comment is preserved below.
- Several helper routines are vectorized over caller-provided array shapes
  instead of using explicit Fortran `nz/ngrdcol` loops.

All pure-jnp -> differentiable where the called PDF path supports it.

References:
  src/CLUBB_core/pdf_closure_module.F90, calc_{wp2xp,wpxp2,wp2xp2,wp4,wpxpyp}_pdf,
  calc_w_up_in_cloud.
"""

from functools import partial

import jax
import jax.numpy as jnp
import jax.scipy.special as jsp

from clubb_jax.src.CLUBB_core.constants_clubb import (
    eps, ep, ep1, ep2, Lv, Rd, Cp, chi_tol, eta_tol, max_num_stdevs, max_mag_correlation,
    min_max_smth_mag, T_freeze_K, sqrt_2, sqrt_2pi, rt_tol, thl_tol, w_tol_sqd,
    w_tol, zero_threshold, rc_tol, p0, kappa,
)
from clubb_jax.src.CLUBB_core.model_flags import (
    iiPDF_ADG1, iiPDF_ADG2, iiPDF_3D_Luhar, iiPDF_new, iiPDF_TSDADG,
    iiPDF_LY93, iiPDF_new_hybrid, l_gamma_Skw,
)
from clubb_jax.src.CLUBB_core.parameter_indices import (
    ibeta, islope_coef_spread_DG_means_w, ipdf_component_stdev_factor_w,
)
from clubb_jax.src.CLUBB_core.saturation import sat_mixrat_ice, sat_mixrat_liq
from clubb_jax.src.CLUBB_core.clip_explicit import clip_rcm
from clubb_jax.src.CLUBB_core.adg1_adg2_3d_luhar_pdf import (
    calc_comp_corrs_binormal, ADG1_pdf_driver, ADG1_ADG2_responder_params,
    ADG2_pdf_driver, Luhar_3D_pdf_driver,
)
from clubb_jax.src.CLUBB_core.new_pdf_main import new_pdf_driver
from clubb_jax.src.CLUBB_core.new_hybrid_pdf_main import new_hybrid_pdf_driver
from clubb_jax.src.CLUBB_core.new_tsdadg_pdf import tsdadg_pdf_driver
from clubb_jax.src.CLUBB_core.LY93_pdf import LY93_driver
from clubb_jax.src.CLUBB_core.Skx_module import Skx_func, compute_gamma_Skw
from clubb_jax.src.CLUBB_core.grid_class import zt2zm, zm2zt
from clubb_jax.src.CLUBB_core.pdf_utilities import calc_corr_chi_x, calc_corr_eta_x
from clubb_jax.src.CLUBB_core.pdf_params import (
    init_pdf_implicit_coefs_terms_api,
)


def _safe_sqrt(value):
    return jnp.sqrt(jnp.maximum(value, 0.0))


def _unsupported_pdf_type_message(iiPDF_type):
    names = {
        iiPDF_ADG1: "iiPDF_ADG1",
        iiPDF_ADG2: "iiPDF_ADG2",
        iiPDF_3D_Luhar: "iiPDF_3D_Luhar",
        iiPDF_new: "iiPDF_new",
        iiPDF_TSDADG: "iiPDF_TSDADG",
        iiPDF_LY93: "iiPDF_LY93",
        iiPDF_new_hybrid: "iiPDF_new_hybrid",
    }
    return f"Unsupported JAX pdf_closure_driver PDF type: {names.get(iiPDF_type, iiPDF_type)}."


def pdf_closure(
    nz, ngrdcol, sclr_dim, sclr_tol, gr,
    hydromet_dim, p_in_Pa, exner, thv_ds,
    wm, wp2, wp3,
    Skw, Skthl_in, Skrt_in, Sku_in, Skv_in,
    rtm, rtp2, wprtp,
    thlm, thlp2, wpthlp,
    um, up2, upwp,
    vm, vp2, vpwp,
    rtpthlp,
    sclrm, wpsclrp, sclrp2,
    sclrprtp, sclrpthlp, Sksclr_in,
    wphydrometp, wp2hmp,
    rtphmp, thlphmp,
    clubb_params, mixt_frac_max_mag,
    saturation_formula,
    stats,
    iiPDF_type,
    l_mix_rat_hm,
    sigma_sqd_w,
    pdf_params, pdf_implicit_coefs_terms,
    err_info,
):
    """Subroutine that computes pdf parameters analytically.

    Based of the original formulation, but with some tweaks
    to remove some of the less realistic assumptions and
    improve transport terms.

    Corrected version that should remove inconsistency

    References:
      The shape of CLUBB's PDF is given by the expression in
      https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:clubb_pdf

      Eqn. 29, 30, 31, 32 & 33  on p. 3547 of
      ``A PDF-Based Model for Boundary Layer Clouds. Part I:
      Method and Model Description'' Golaz, et al. (2002)
      JAS, Vol. 59, pp. 3540--3551.

    !#######################################################################
    !#######################################################################
    ! If you change the argument list of pdf_closure you also have to
    ! change the calls to this function in the host models CAM, WRF, SAM
    ! and GFDL.
    !#######################################################################
    !#######################################################################
    """
    # TODO: port the hydrometeor/liquid-ice-loading branch from the Fortran routine
    # before claiming full pdf_closure parity for hydrometeor prognostic cases.
    del hydromet_dim, wphydrometp, wp2hmp, rtphmp, thlphmp
    # The JAX stats recorder is handled in pdf_closure_driver; this routine still
    # computes moment fields unconditionally instead of using Fortran's stats-gated
    # shortcuts.
    del l_mix_rat_hm, stats, wp3

    if iiPDF_type not in (
        iiPDF_ADG1,
        iiPDF_ADG2,
        iiPDF_3D_Luhar,
        iiPDF_new,
        iiPDF_TSDADG,
        iiPDF_LY93,
        iiPDF_new_hybrid,
    ):
        raise NotImplementedError(_unsupported_pdf_type_message(iiPDF_type))

    sqrt_wp2 = jnp.sqrt(wp2)
    half = jnp.full_like(wp2, 0.5)
    sclr1 = None
    sclr2 = None
    varnce_sclr1 = None
    varnce_sclr2 = None
    if iiPDF_type == iiPDF_ADG1:
        # Calculate the mixture fraction for the multivariate PDF, as well as both
        # PDF component means and both PDF component variances for each of w, rt,
        # theta-l, and passive scalar variables.
        adg1 = ADG1_pdf_driver(
            wm=jnp.asarray(wm),
            rtm=jnp.asarray(rtm),
            thlm=jnp.asarray(thlm),
            um=jnp.asarray(um),
            vm=jnp.asarray(vm),
            wp2=wp2,
            rtp2=rtp2,
            thlp2=thlp2,
            up2=up2,
            vp2=vp2,
            Skw=Skw,
            wprtp=wprtp,
            wpthlp=wpthlp,
            upwp=upwp,
            vpwp=vpwp,
            sqrt_wp2=sqrt_wp2,
            sigma_sqd_w=sigma_sqd_w,
            beta=jnp.asarray(clubb_params)[:, ibeta],
            mixt_frac_max_mag=mixt_frac_max_mag,
        )
    elif iiPDF_type == iiPDF_ADG2:
        # Calculate the mixture fraction for the multivariate PDF, as well as both
        # PDF component means and both PDF component variances for each of w, rt,
        # theta-l, and passive scalar variables.
        (
            w_1, w_2, rt_1, rt_2, thl_1, thl_2,
            varnce_w_1, varnce_w_2, varnce_rt_1, varnce_rt_2,
            varnce_thl_1, varnce_thl_2, mixt_frac,
            alpha_rt, alpha_thl, sigma_sqd_w,
            sclr1, sclr2, varnce_sclr1, varnce_sclr2, _alpha_sclr,
        ) = ADG2_pdf_driver(
            nz, ngrdcol, sclr_dim, sclr_tol,
            wm, rtm, thlm, wp2, rtp2, thlp2,
            Skw, wprtp, wpthlp, sqrt_wp2,
            jnp.asarray(clubb_params)[:, ibeta],
            sclrm, sclrp2, wpsclrp, sclr_dim > 0,
        )
        adg1 = {
            "w_1": w_1,
            "w_2": w_2,
            "varnce_w_1": varnce_w_1,
            "varnce_w_2": varnce_w_2,
            "mixt_frac": mixt_frac,
            "rt_1": rt_1,
            "rt_2": rt_2,
            "varnce_rt_1": varnce_rt_1,
            "varnce_rt_2": varnce_rt_2,
            "alpha_rt": alpha_rt,
            "thl_1": thl_1,
            "thl_2": thl_2,
            "varnce_thl_1": varnce_thl_1,
            "varnce_thl_2": varnce_thl_2,
            "alpha_thl": alpha_thl,
            "u_1": um,
            "u_2": um,
            "varnce_u_1": up2,
            "varnce_u_2": up2,
            "alpha_u": half,
            "v_1": vm,
            "v_2": vm,
            "varnce_v_1": vp2,
            "varnce_v_2": vp2,
            "alpha_v": half,
        }
    elif iiPDF_type == iiPDF_3D_Luhar:
        # Calculate the mixture fraction for the multivariate PDF, as well as both
        # PDF component means and both PDF component variances for each of w, rt,
        # theta-l, and passive scalar variables.
        (
            w_1, w_2, rt_1, rt_2, thl_1, thl_2,
            varnce_w_1, varnce_w_2, varnce_rt_1, varnce_rt_2,
            varnce_thl_1, varnce_thl_2, mixt_frac,
        ) = Luhar_3D_pdf_driver(
            wm, rtm, thlm, wp2, rtp2, thlp2,
            Skw, Skrt_in, Skthl_in, wprtp, wpthlp,
        )
        adg1 = {
            "w_1": w_1,
            "w_2": w_2,
            "varnce_w_1": varnce_w_1,
            "varnce_w_2": varnce_w_2,
            "mixt_frac": mixt_frac,
            "rt_1": rt_1,
            "rt_2": rt_2,
            "varnce_rt_1": varnce_rt_1,
            "varnce_rt_2": varnce_rt_2,
            "alpha_rt": half,
            "thl_1": thl_1,
            "thl_2": thl_2,
            "varnce_thl_1": varnce_thl_1,
            "varnce_thl_2": varnce_thl_2,
            "alpha_thl": half,
            "u_1": um,
            "u_2": um,
            "varnce_u_1": up2,
            "varnce_u_2": up2,
            "alpha_u": half,
            "v_1": vm,
            "v_2": vm,
            "varnce_v_1": vp2,
            "varnce_v_2": vp2,
            "alpha_v": half,
        }
    elif iiPDF_type == iiPDF_new:
        # Calculate the mixture fraction for the multivariate PDF, as well as both
        # PDF component means and both PDF component variances for each of w, rt,
        # theta-l, and passive scalar variables.
        new_pdf = new_pdf_driver(
            wm, rtm, thlm, wp2, rtp2, thlp2,
            Skw, wprtp, wpthlp, rtpthlp,
            clubb_params, Skrt_in, Skthl_in,
            pdf_implicit_coefs_terms,
        )
        pdf_implicit_coefs_terms = new_pdf["pdf_implicit_coefs_terms"]
        adg1 = {
            "w_1": new_pdf["mu_w_1"],
            "w_2": new_pdf["mu_w_2"],
            "varnce_w_1": new_pdf["sigma_w_1_sqd"],
            "varnce_w_2": new_pdf["sigma_w_2_sqd"],
            "mixt_frac": new_pdf["mixt_frac"],
            "rt_1": new_pdf["mu_rt_1"],
            "rt_2": new_pdf["mu_rt_2"],
            "varnce_rt_1": new_pdf["sigma_rt_1_sqd"],
            "varnce_rt_2": new_pdf["sigma_rt_2_sqd"],
            "alpha_rt": half,
            "thl_1": new_pdf["mu_thl_1"],
            "thl_2": new_pdf["mu_thl_2"],
            "varnce_thl_1": new_pdf["sigma_thl_1_sqd"],
            "varnce_thl_2": new_pdf["sigma_thl_2_sqd"],
            "alpha_thl": half,
            "u_1": um,
            "u_2": um,
            "varnce_u_1": up2,
            "varnce_u_2": up2,
            "alpha_u": half,
            "v_1": vm,
            "v_2": vm,
            "varnce_v_1": vp2,
            "varnce_v_2": vp2,
            "alpha_v": half,
        }
    elif iiPDF_type == iiPDF_TSDADG:
        # Calculate the mixture fraction for the multivariate PDF, as well as both
        # PDF component means and both PDF component variances for each of w, rt,
        # theta-l, and passive scalar variables.
        (
            w_1, w_2, rt_1, rt_2, thl_1, thl_2,
            varnce_w_1, varnce_w_2, varnce_rt_1, varnce_rt_2,
            varnce_thl_1, varnce_thl_2, mixt_frac,
        ) = tsdadg_pdf_driver(
            wm, rtm, thlm, wp2, rtp2, thlp2,
            Skw, Skrt_in, Skthl_in, wprtp, wpthlp,
        )
        adg1 = {
            "w_1": w_1,
            "w_2": w_2,
            "varnce_w_1": varnce_w_1,
            "varnce_w_2": varnce_w_2,
            "mixt_frac": mixt_frac,
            "rt_1": rt_1,
            "rt_2": rt_2,
            "varnce_rt_1": varnce_rt_1,
            "varnce_rt_2": varnce_rt_2,
            "alpha_rt": half,
            "thl_1": thl_1,
            "thl_2": thl_2,
            "varnce_thl_1": varnce_thl_1,
            "varnce_thl_2": varnce_thl_2,
            "alpha_thl": half,
            "u_1": um,
            "u_2": um,
            "varnce_u_1": up2,
            "varnce_u_2": up2,
            "alpha_u": half,
            "v_1": vm,
            "v_2": vm,
            "varnce_v_1": vp2,
            "varnce_v_2": vp2,
            "alpha_v": half,
        }
    elif iiPDF_type == iiPDF_LY93:
        # Calculate the mixture fraction for the multivariate PDF, as well as both
        # PDF component means and both PDF component variances for each of w, rt,
        # theta-l, and passive scalar variables.
        (
            w_1, w_2, rt_1, rt_2, thl_1, thl_2,
            varnce_w_1, varnce_w_2, varnce_rt_1, varnce_rt_2,
            varnce_thl_1, varnce_thl_2, mixt_frac,
        ) = LY93_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw, Skrt_in, Skthl_in)
        adg1 = {
            "w_1": w_1,
            "w_2": w_2,
            "varnce_w_1": varnce_w_1,
            "varnce_w_2": varnce_w_2,
            "mixt_frac": mixt_frac,
            "rt_1": rt_1,
            "rt_2": rt_2,
            "varnce_rt_1": varnce_rt_1,
            "varnce_rt_2": varnce_rt_2,
            "alpha_rt": half,
            "thl_1": thl_1,
            "thl_2": thl_2,
            "varnce_thl_1": varnce_thl_1,
            "varnce_thl_2": varnce_thl_2,
            "alpha_thl": half,
            "u_1": um,
            "u_2": um,
            "varnce_u_1": up2,
            "varnce_u_2": up2,
            "alpha_u": half,
            "v_1": vm,
            "v_2": vm,
            "varnce_v_1": vp2,
            "varnce_v_2": vp2,
            "alpha_v": half,
        }
    elif iiPDF_type == iiPDF_new_hybrid:
        # Calculate the mixture fraction for the multivariate PDF, as well as both
        # PDF component means and both PDF component variances for each of w, rt,
        # theta-l, and passive scalar variables.
        gamma_Skw_fnc = compute_gamma_Skw(
            nz, ngrdcol, Skw, clubb_params, l_gamma_Skw,
        )
        new_hybrid = new_hybrid_pdf_driver(
            wm, rtm, thlm, um, vm,
            wp2, rtp2, thlp2, up2, vp2,
            Skw, wprtp, wpthlp, upwp, vpwp,
            gamma_Skw_fnc,
            jnp.asarray(clubb_params)[:, islope_coef_spread_DG_means_w],
            jnp.asarray(clubb_params)[:, ipdf_component_stdev_factor_w],
            Skrt_in, Skthl_in, Sku_in, Skv_in,
            pdf_implicit_coefs_terms,
            sclrm=sclrm if sclr_dim > 0 else None,
            sclrp2=sclrp2 if sclr_dim > 0 else None,
            wpsclrp=wpsclrp if sclr_dim > 0 else None,
            Sksclr=Sksclr_in if sclr_dim > 0 else None,
        )
        pdf_implicit_coefs_terms = new_hybrid["pdf_implicit_coefs_terms"]
        adg1 = {
            "w_1": new_hybrid["mu_w_1"],
            "w_2": new_hybrid["mu_w_2"],
            "varnce_w_1": new_hybrid["sigma_w_1_sqd"],
            "varnce_w_2": new_hybrid["sigma_w_2_sqd"],
            "mixt_frac": new_hybrid["mixt_frac"],
            "rt_1": new_hybrid["mu_rt_1"],
            "rt_2": new_hybrid["mu_rt_2"],
            "varnce_rt_1": new_hybrid["sigma_rt_1_sqd"],
            "varnce_rt_2": new_hybrid["sigma_rt_2_sqd"],
            "alpha_rt": half,
            "thl_1": new_hybrid["mu_thl_1"],
            "thl_2": new_hybrid["mu_thl_2"],
            "varnce_thl_1": new_hybrid["sigma_thl_1_sqd"],
            "varnce_thl_2": new_hybrid["sigma_thl_2_sqd"],
            "alpha_thl": half,
            "u_1": new_hybrid["mu_u_1"],
            "u_2": new_hybrid["mu_u_2"],
            "varnce_u_1": new_hybrid["sigma_u_1_sqd"],
            "varnce_u_2": new_hybrid["sigma_u_2_sqd"],
            "alpha_u": half,
            "v_1": new_hybrid["mu_v_1"],
            "v_2": new_hybrid["mu_v_2"],
            "varnce_v_1": new_hybrid["sigma_v_1_sqd"],
            "varnce_v_2": new_hybrid["sigma_v_2_sqd"],
            "alpha_v": half,
        }

    if sclr_dim > 0:
        if iiPDF_type == iiPDF_ADG1:
            sclr1, sclr2, varnce_sclr1, varnce_sclr2, _alpha_sclr = jax.vmap(
                lambda sclrm_s, sclrp2_s, wpsclrp_s: ADG1_ADG2_responder_params(
                    sclrm_s,
                    sclrp2_s,
                    wp2,
                    sqrt_wp2,
                    wpsclrp_s,
                    adg1["w_1_n"],
                    adg1["w_2_n"],
                    adg1["mixt_frac"],
                    sigma_sqd_w,
                    jnp.asarray(clubb_params)[:, ibeta],
                ),
                in_axes=(2, 2, 2),
                out_axes=-1,
            )(sclrm, sclrp2, wpsclrp)
        elif iiPDF_type == iiPDF_new_hybrid:
            sclr1 = new_hybrid["mu_sclr_1"]
            sclr2 = new_hybrid["mu_sclr_2"]
            varnce_sclr1 = new_hybrid["sigma_sclr_1_sqd"]
            varnce_sclr2 = new_hybrid["sigma_sclr_2_sqd"]
        elif sclr1 is None:
            sclr1 = sclrm
            sclr2 = sclrm
            varnce_sclr1 = sclrp2
            varnce_sclr2 = sclrp2
    else:
        sclr1 = jnp.zeros((ngrdcol, nz, 0), dtype=jnp.float64)
        sclr2 = jnp.zeros((ngrdcol, nz, 0), dtype=jnp.float64)
        varnce_sclr1 = jnp.zeros((ngrdcol, nz, 0), dtype=jnp.float64)
        varnce_sclr2 = jnp.zeros((ngrdcol, nz, 0), dtype=jnp.float64)

    # Calculate the PDF component correlations of rt and thl.
    corr_rt_thl_1, corr_rt_thl_2 = calc_comp_corrs_binormal(
        rtpthlp, rtm, thlm,
        adg1["rt_1"], adg1["rt_2"],
        adg1["thl_1"], adg1["thl_2"],
        adg1["varnce_rt_1"], adg1["varnce_rt_2"],
        adg1["varnce_thl_1"], adg1["varnce_thl_2"],
        adg1["mixt_frac"],
    )
    if iiPDF_type in (iiPDF_ADG1, iiPDF_ADG2, iiPDF_new_hybrid):
        # These PDF types define corr_w_rt_1, corr_w_rt_2, corr_w_thl_1, and
        # corr_w_thl_2 to all have a value of 0, so skip the calculation.
        # The values of corr_u_w_1, corr_u_w_2, corr_v_w_1, and corr_v_w_2 are
        # all defined to be 0, as well.
        corr_w_rt_1 = jnp.zeros_like(adg1["mixt_frac"])
        corr_w_rt_2 = jnp.zeros_like(adg1["mixt_frac"])
        corr_w_thl_1 = jnp.zeros_like(adg1["mixt_frac"])
        corr_w_thl_2 = jnp.zeros_like(adg1["mixt_frac"])
    else:
        # Calculate the PDF component correlations of w and rt.
        corr_w_rt_1, corr_w_rt_2 = calc_comp_corrs_binormal(
            wprtp, wm, rtm,
            adg1["w_1"], adg1["w_2"],
            adg1["rt_1"], adg1["rt_2"],
            adg1["varnce_w_1"], adg1["varnce_w_2"],
            adg1["varnce_rt_1"], adg1["varnce_rt_2"],
            adg1["mixt_frac"],
        )
        # Calculate the PDF component correlations of w and thl.
        corr_w_thl_1, corr_w_thl_2 = calc_comp_corrs_binormal(
            wpthlp, wm, thlm,
            adg1["w_1"], adg1["w_2"],
            adg1["thl_1"], adg1["thl_2"],
            adg1["varnce_w_1"], adg1["varnce_w_2"],
            adg1["varnce_thl_1"], adg1["varnce_thl_2"],
            adg1["mixt_frac"],
        )
    if sclr_dim > 0:
        # Calculate the PDF component correlations of a passive scalar and thl.
        corr_sclr_thl_1, corr_sclr_thl_2 = jax.vmap(
            lambda sclrpthlp_s, sclrm_s, sclr1_s, sclr2_s, varnce_sclr1_s, varnce_sclr2_s: (
                calc_comp_corrs_binormal(
                    sclrpthlp_s,
                    sclrm_s,
                    thlm,
                    sclr1_s,
                    sclr2_s,
                    adg1["thl_1"],
                    adg1["thl_2"],
                    varnce_sclr1_s,
                    varnce_sclr2_s,
                    adg1["varnce_thl_1"],
                    adg1["varnce_thl_2"],
                    adg1["mixt_frac"],
                )
            ),
            in_axes=(2, 2, 2, 2, 2, 2),
            out_axes=-1,
        )(sclrpthlp, sclrm, sclr1, sclr2, varnce_sclr1, varnce_sclr2)
        # Calculate the PDF component correlations of a passive scalar and rt.
        corr_sclr_rt_1, corr_sclr_rt_2 = jax.vmap(
            lambda sclrprtp_s, sclrm_s, sclr1_s, sclr2_s, varnce_sclr1_s, varnce_sclr2_s: (
                calc_comp_corrs_binormal(
                    sclrprtp_s,
                    sclrm_s,
                    rtm,
                    sclr1_s,
                    sclr2_s,
                    adg1["rt_1"],
                    adg1["rt_2"],
                    varnce_sclr1_s,
                    varnce_sclr2_s,
                    adg1["varnce_rt_1"],
                    adg1["varnce_rt_2"],
                    adg1["mixt_frac"],
                )
            ),
            in_axes=(2, 2, 2, 2, 2, 2),
            out_axes=-1,
        )(sclrprtp, sclrm, sclr1, sclr2, varnce_sclr1, varnce_sclr2)
        if iiPDF_type in (iiPDF_ADG1, iiPDF_ADG2, iiPDF_new_hybrid):
            # These PDF types define all PDF component correlations involving w
            # to have a value of 0, so skip the calculation.
            corr_w_sclr_1 = jnp.zeros_like(sclr1)
            corr_w_sclr_2 = jnp.zeros_like(sclr1)
        else:
            # Calculate the PDF component correlations of w and a passive scalar.
            corr_w_sclr_1, corr_w_sclr_2 = jax.vmap(
                lambda wpsclrp_s, sclrm_s, sclr1_s, sclr2_s, varnce_sclr1_s, varnce_sclr2_s: (
                    calc_comp_corrs_binormal(
                        wpsclrp_s,
                        wm,
                        sclrm_s,
                        adg1["w_1"],
                        adg1["w_2"],
                        sclr1_s,
                        sclr2_s,
                        adg1["varnce_w_1"],
                        adg1["varnce_w_2"],
                        varnce_sclr1_s,
                        varnce_sclr2_s,
                        adg1["mixt_frac"],
                    )
                ),
                in_axes=(2, 2, 2, 2, 2, 2),
                out_axes=-1,
            )(wpsclrp, sclrm, sclr1, sclr2, varnce_sclr1, varnce_sclr2)
    else:
        corr_sclr_thl_1 = jnp.zeros((ngrdcol, nz, 0), dtype=jnp.float64)
        corr_sclr_thl_2 = jnp.zeros((ngrdcol, nz, 0), dtype=jnp.float64)
        corr_sclr_rt_1 = jnp.zeros((ngrdcol, nz, 0), dtype=jnp.float64)
        corr_sclr_rt_2 = jnp.zeros((ngrdcol, nz, 0), dtype=jnp.float64)
        corr_w_sclr_1 = jnp.zeros((ngrdcol, nz, 0), dtype=jnp.float64)
        corr_w_sclr_2 = jnp.zeros((ngrdcol, nz, 0), dtype=jnp.float64)

    mf = adg1["mixt_frac"]
    zero_pdf = jnp.zeros_like(mf)

    # Compute higher order moments that include theta_v.

    # First compute some preliminary quantities.
    # "1" denotes first Gaussian; "2" denotes 2nd Gaussian
    # liq water temp (Sommeria & Deardorff 1977 (SD), eqn. 3)
    tl1 = adg1["thl_1"] * exner
    tl2 = adg1["thl_2"] * exner
    rsatl_1 = sat_mixrat_liq(p_in_Pa, tl1, saturation_formula)
    rsatl_2 = sat_mixrat_liq(p_in_Pa, tl2, saturation_formula)
    chi_1, crt_1, cthl_1, stdev_chi_1, stdev_eta_1, covar_chi_eta_1, corr_chi_eta_1 = (
        transform_pdf_chi_eta_component(
            tl1, rsatl_1, adg1["rt_1"], exner,
            adg1["varnce_thl_1"], adg1["varnce_rt_1"],
            corr_rt_thl_1,
        )
    )
    # Calculate cloud fraction component for pdf 1
    cloud_frac_1, rc_1 = calc_liquid_cloud_frac_component(chi_1, stdev_chi_1)
    # Calc ice_supersat_frac
    ice_supersat_frac_1 = calc_ice_cloud_frac_component(
        chi_1, stdev_chi_1, crt_1, rsatl_1, tl1,
        cloud_frac_1, p_in_Pa,
    )
    chi_2, crt_2, cthl_2, stdev_chi_2, stdev_eta_2, covar_chi_eta_2, corr_chi_eta_2 = (
        transform_pdf_chi_eta_component(
            tl2, rsatl_2, adg1["rt_2"], exner,
            adg1["varnce_thl_2"], adg1["varnce_rt_2"],
            corr_rt_thl_2,
        )
    )
    # Calculate cloud fraction component for pdf 2
    cloud_frac_2, rc_2 = calc_liquid_cloud_frac_component(chi_2, stdev_chi_2)
    # Calc ice_supersat_frac
    ice_supersat_frac_2 = calc_ice_cloud_frac_component(
        chi_2, stdev_chi_2, crt_2, rsatl_2, tl2,
        cloud_frac_2, p_in_Pa,
    )
    ice_supersat_frac = (
        mf * ice_supersat_frac_1
        + (1.0 - mf) * ice_supersat_frac_2
    )
    # Compute cloud fraction and mean cloud water mixing ratio.
    # Reference:
    # https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:anl_int_cloud_terms
    cloud_frac = mf * cloud_frac_1 + (1.0 - mf) * cloud_frac_2
    rcm = jnp.maximum(zero_threshold, mf * rc_1 + (1.0 - mf) * rc_2)
    if iiPDF_type in (iiPDF_ADG1, iiPDF_ADG2, iiPDF_new_hybrid):
        # corr_w_rt and corr_w_thl are zero for these pdf types so
        # corr_w_chi and corr_w_eta are zero as well
        corr_w_chi_1 = zero_pdf
        corr_w_chi_2 = zero_pdf
        corr_w_eta_1 = zero_pdf
        corr_w_eta_2 = zero_pdf
    else:
        # Correlation of w and chi for each component.
        corr_w_chi_1 = calc_corr_chi_x(
            crt_1, cthl_1,
            _safe_sqrt(adg1["varnce_rt_1"]),
            _safe_sqrt(adg1["varnce_thl_1"]),
            stdev_chi_1,
            corr_w_rt_1,
            corr_w_thl_1,
        )
        corr_w_chi_2 = calc_corr_chi_x(
            crt_2, cthl_2,
            _safe_sqrt(adg1["varnce_rt_2"]),
            _safe_sqrt(adg1["varnce_thl_2"]),
            stdev_chi_2,
            corr_w_rt_2,
            corr_w_thl_2,
        )
        # Correlation of w and eta for each component.
        corr_w_eta_1 = calc_corr_eta_x(
            crt_1, cthl_1,
            _safe_sqrt(adg1["varnce_rt_1"]),
            _safe_sqrt(adg1["varnce_thl_1"]),
            stdev_eta_1,
            corr_w_rt_1,
            corr_w_thl_1,
        )
        corr_w_eta_2 = calc_corr_eta_x(
            crt_2, cthl_2,
            _safe_sqrt(adg1["varnce_rt_2"]),
            _safe_sqrt(adg1["varnce_thl_2"]),
            stdev_eta_2,
            corr_w_rt_2,
            corr_w_thl_2,
        )

    # Compute moments that depend on theta_v
    #
    # The moments that depend on th_v' are calculated based on an approximated
    # and linearized form of the theta_v equation:
    #
    # theta_v = theta_l + { (R_v/R_d) - 1 } * thv_ds * r_t
    #                   + [ {L_v/(C_p*exner)} - (R_v/R_d) * thv_ds ] * r_c;
    #
    # and therefore:
    #
    # th_v' = th_l' + { (R_v/R_d) - 1 } * thv_ds * r_t'
    #               + [ {L_v/(C_p*exner)} - (R_v/R_d) * thv_ds ] * r_c';
    #
    # where thv_ds is used as a reference value to approximate theta_l.
    #
    # Reference:
    # https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:anl_int_buoy_terms
    #
    # Calculate the contributions to <w'rc'>, <w'^2 rc'>, <rt'rc'>, and
    # <thl'rc'> from the 1st PDF component.
    (
        wprcp_1, wp2rcp_1, rtprcp_1, thlprcp_1,
        uprcp_1, vprcp_1,
    ) = calc_xprcp_component(
        wm, rtm, thlm, um, vm, rcm,
        adg1["w_1"], adg1["rt_1"], adg1["thl_1"], adg1["u_1"], adg1["v_1"],
        adg1["varnce_w_1"], chi_1, stdev_chi_1, stdev_eta_1,
        corr_w_chi_1, corr_chi_eta_1, crt_1, cthl_1, rc_1, cloud_frac_1,
        iiPDF_type,
    )
    # Calculate the contributions to <w'rc'>, <w'^2 rc'>, <rt'rc'>, and
    # <thl'rc'> from the 2nd PDF component.
    (
        wprcp_2, wp2rcp_2, rtprcp_2, thlprcp_2,
        uprcp_2, vprcp_2,
    ) = calc_xprcp_component(
        wm, rtm, thlm, um, vm, rcm,
        adg1["w_2"], adg1["rt_2"], adg1["thl_2"], adg1["u_2"], adg1["v_2"],
        adg1["varnce_w_2"], chi_2, stdev_chi_2, stdev_eta_2,
        corr_w_chi_2, corr_chi_eta_2, crt_2, cthl_2, rc_2, cloud_frac_2,
        iiPDF_type,
    )
    # Calculate <w'rc'>, <w'^2 rc'>, <rt'rc'>, and <thl'rc'>.
    wprcp = mf * wprcp_1 + (1.0 - mf) * wprcp_2
    wp2rcp = mf * wp2rcp_1 + (1.0 - mf) * wp2rcp_2
    rtprcp = mf * rtprcp_1 + (1.0 - mf) * rtprcp_2
    thlprcp = mf * thlprcp_1 + (1.0 - mf) * thlprcp_2
    uprcp = mf * uprcp_1 + (1.0 - mf) * uprcp_2
    vprcp = mf * vprcp_1 + (1.0 - mf) * vprcp_2

    w1, w2 = adg1["w_1"], adg1["w_2"]
    vw1, vw2 = adg1["varnce_w_1"], adg1["varnce_w_2"]
    # Compute higher order moments (these are interactive)
    wp2rtp = calc_wp2xp_pdf(
        wm, rtm, w1, w2,
        adg1["rt_1"], adg1["rt_2"],
        vw1, vw2, adg1["varnce_rt_1"], adg1["varnce_rt_2"],
        corr_w_rt_1, corr_w_rt_2, mf,
    )
    wp2thlp = calc_wp2xp_pdf(
        wm, thlm, w1, w2,
        adg1["thl_1"], adg1["thl_2"],
        vw1, vw2, adg1["varnce_thl_1"], adg1["varnce_thl_2"],
        corr_w_thl_1, corr_w_thl_2, mf,
    )
    wp2up = calc_wp2xp_pdf(
        wm, um, w1, w2,
        adg1["u_1"], adg1["u_2"],
        vw1, vw2, adg1["varnce_u_1"], adg1["varnce_u_2"],
        zero_pdf, zero_pdf, mf,
    )
    # Compute higher order moments (these may be interactive)
    wpup2 = calc_wpxp2_pdf(
        wm, um, w1, w2,
        adg1["u_1"], adg1["u_2"],
        vw1, vw2, adg1["varnce_u_1"], adg1["varnce_u_2"],
        zero_pdf, zero_pdf, mf,
    )
    wpvp2 = calc_wpxp2_pdf(
        wm, vm, w1, w2,
        adg1["v_1"], adg1["v_2"],
        vw1, vw2, adg1["varnce_v_1"], adg1["varnce_v_2"],
        zero_pdf, zero_pdf, mf,
    )
    wp2up2 = calc_wp2xp2_pdf(
        wm, um, w1, w2,
        adg1["u_1"], adg1["u_2"],
        vw1, vw2, adg1["varnce_u_1"], adg1["varnce_u_2"],
        zero_pdf, zero_pdf, mf,
    )
    wp2vp2 = calc_wp2xp2_pdf(
        wm, vm, w1, w2,
        adg1["v_1"], adg1["v_2"],
        vw1, vw2, adg1["varnce_v_1"], adg1["varnce_v_2"],
        zero_pdf, zero_pdf, mf,
    )
    wp4 = calc_wp4_pdf(wm, w1, w2, vw1, vw2, mf)
    wprtp2 = calc_wpxp2_pdf(
        wm, rtm, w1, w2,
        adg1["rt_1"], adg1["rt_2"],
        vw1, vw2, adg1["varnce_rt_1"], adg1["varnce_rt_2"],
        corr_w_rt_1, corr_w_rt_2, mf,
    )
    wpthlp2 = calc_wpxp2_pdf(
        wm, thlm, w1, w2,
        adg1["thl_1"], adg1["thl_2"],
        vw1, vw2, adg1["varnce_thl_1"], adg1["varnce_thl_2"],
        corr_w_thl_1, corr_w_thl_2, mf,
    )
    wprtpthlp = calc_wpxpyp_pdf(
        wm, rtm, thlm, w1, w2,
        adg1["rt_1"], adg1["rt_2"], adg1["thl_1"], adg1["thl_2"],
        vw1, vw2, adg1["varnce_rt_1"], adg1["varnce_rt_2"],
        adg1["varnce_thl_1"], adg1["varnce_thl_2"],
        corr_w_rt_1, corr_w_rt_2, corr_w_thl_1, corr_w_thl_2,
        corr_rt_thl_1, corr_rt_thl_2, mf,
    )

    # Calculate rc_coef, which is the coefficient on <x'rc'> in the <x'thv'> equation.
    rc_coef = Lv / (exner * Cp) - ep2 * thv_ds
    # Calculate <w'thv'>, <w'^2 thv'>, <rt'thv'>, and <thl'thv'>.
    wpthvp = wpthlp + ep1 * thv_ds * wprtp + rc_coef * wprcp
    wp2thvp = wp2thlp + ep1 * thv_ds * wp2rtp + rc_coef * wp2rcp
    rtpthvp = rtpthlp + ep1 * thv_ds * rtp2 + rc_coef * rtprcp
    thlpthvp = thlp2 + ep1 * thv_ds * rtpthlp + rc_coef * thlprcp
    rcp2 = jnp.maximum(
        zero_threshold,
        mf * (chi_1 * rc_1 + cloud_frac_1 * stdev_chi_1 ** 2)
        + (1.0 - mf) * (
            chi_2 * rc_2 + cloud_frac_2 * stdev_chi_2 ** 2
        )
        - rcm ** 2,
    )

    if sclr_dim > 0:
        # Scalar Addition to higher order moments
        wp2sclrp = calc_wp2xp_pdf(
            wm[:, :, None],
            sclrm,
            w1[:, :, None],
            w2[:, :, None],
            sclr1,
            sclr2,
            vw1[:, :, None],
            vw2[:, :, None],
            varnce_sclr1,
            varnce_sclr2,
            corr_w_sclr_1,
            corr_w_sclr_2,
            mf[:, :, None],
        )
        wpsclrp2 = calc_wpxp2_pdf(
            wm[:, :, None],
            sclrm,
            w1[:, :, None],
            w2[:, :, None],
            sclr1,
            sclr2,
            vw1[:, :, None],
            vw2[:, :, None],
            varnce_sclr1,
            varnce_sclr2,
            corr_w_sclr_1,
            corr_w_sclr_2,
            mf[:, :, None],
        )
        wpsclrprtp = calc_wpxpyp_pdf(
            wm[:, :, None],
            sclrm,
            rtm[:, :, None],
            w1[:, :, None],
            w2[:, :, None],
            sclr1,
            sclr2,
            adg1["rt_1"][:, :, None],
            adg1["rt_2"][:, :, None],
            vw1[:, :, None],
            vw2[:, :, None],
            varnce_sclr1,
            varnce_sclr2,
            adg1["varnce_rt_1"][:, :, None],
            adg1["varnce_rt_2"][:, :, None],
            corr_w_sclr_1,
            corr_w_sclr_2,
            zero_pdf[:, :, None],
            zero_pdf[:, :, None],
            corr_sclr_rt_1,
            corr_sclr_rt_2,
            mf[:, :, None],
        )
        wpsclrpthlp = calc_wpxpyp_pdf(
            wm[:, :, None],
            sclrm,
            thlm[:, :, None],
            w1[:, :, None],
            w2[:, :, None],
            sclr1,
            sclr2,
            adg1["thl_1"][:, :, None],
            adg1["thl_2"][:, :, None],
            vw1[:, :, None],
            vw2[:, :, None],
            varnce_sclr1,
            varnce_sclr2,
            adg1["varnce_thl_1"][:, :, None],
            adg1["varnce_thl_2"][:, :, None],
            corr_w_sclr_1,
            corr_w_sclr_2,
            zero_pdf[:, :, None],
            zero_pdf[:, :, None],
            corr_sclr_thl_1,
            corr_sclr_thl_2,
            mf[:, :, None],
        )
        sclrprcp = (
            mf[:, :, None] * ((sclr1 - sclrm) * rc_1[:, :, None])
            + (1.0 - mf[:, :, None]) * ((sclr2 - sclrm) * rc_2[:, :, None])
            + mf[:, :, None] * corr_sclr_rt_1 * crt_1[:, :, None]
            * _safe_sqrt(varnce_sclr1 * adg1["varnce_rt_1"][:, :, None])
            * cloud_frac_1[:, :, None]
            + (1.0 - mf[:, :, None]) * corr_sclr_rt_2 * crt_2[:, :, None]
            * _safe_sqrt(varnce_sclr2 * adg1["varnce_rt_2"][:, :, None])
            * cloud_frac_2[:, :, None]
            - mf[:, :, None] * corr_sclr_thl_1 * cthl_1[:, :, None]
            * _safe_sqrt(varnce_sclr1 * adg1["varnce_thl_1"][:, :, None])
            * cloud_frac_1[:, :, None]
            - (1.0 - mf[:, :, None]) * corr_sclr_thl_2 * cthl_2[:, :, None]
            * _safe_sqrt(varnce_sclr2 * adg1["varnce_thl_2"][:, :, None])
            * cloud_frac_2[:, :, None]
        )
        sclrpthvp = (
            sclrpthlp
            + ep1 * thv_ds[:, :, None] * sclrprtp
            + rc_coef[:, :, None] * sclrprcp
        )
    else:
        wp2sclrp = jnp.zeros((ngrdcol, nz, 0), dtype=jnp.float64)
        wpsclrp2 = jnp.zeros((ngrdcol, nz, 0), dtype=jnp.float64)
        wpsclrprtp = jnp.zeros((ngrdcol, nz, 0), dtype=jnp.float64)
        wpsclrpthlp = jnp.zeros((ngrdcol, nz, 0), dtype=jnp.float64)
        sclrprcp = jnp.zeros((ngrdcol, nz, 0), dtype=jnp.float64)
        sclrpthvp = jnp.zeros((ngrdcol, nz, 0), dtype=jnp.float64)

    if iiPDF_type in (iiPDF_ADG1, iiPDF_ADG2, iiPDF_new_hybrid):
        # Calculate mean cloudy updraft speed.
        w_up_in_cloud, w_down_in_cloud, cloudy_updraft_frac, cloudy_downdraft_frac = (
            calc_w_up_in_cloud(
                adg1["mixt_frac"],
                cloud_frac_1,
                cloud_frac_2,
                adg1["w_1"],
                adg1["w_2"],
                adg1["varnce_w_1"],
                adg1["varnce_w_2"],
            )
        )
    else:
        w_up_in_cloud = jnp.zeros_like(mf)
        w_down_in_cloud = jnp.zeros_like(mf)
        cloudy_updraft_frac = jnp.zeros_like(mf)
        cloudy_downdraft_frac = jnp.zeros_like(mf)

    pdf_params = pdf_params.replace(
        w_1=adg1["w_1"],
        w_2=adg1["w_2"],
        varnce_w_1=adg1["varnce_w_1"],
        varnce_w_2=adg1["varnce_w_2"],
        rt_1=adg1["rt_1"],
        rt_2=adg1["rt_2"],
        varnce_rt_1=adg1["varnce_rt_1"],
        varnce_rt_2=adg1["varnce_rt_2"],
        thl_1=adg1["thl_1"],
        thl_2=adg1["thl_2"],
        varnce_thl_1=adg1["varnce_thl_1"],
        varnce_thl_2=adg1["varnce_thl_2"],
        corr_w_rt_1=corr_w_rt_1,
        corr_w_rt_2=corr_w_rt_2,
        corr_w_thl_1=corr_w_thl_1,
        corr_w_thl_2=corr_w_thl_2,
        corr_rt_thl_1=corr_rt_thl_1,
        corr_rt_thl_2=corr_rt_thl_2,
        alpha_thl=adg1["alpha_thl"],
        alpha_rt=adg1["alpha_rt"],
        crt_1=crt_1,
        crt_2=crt_2,
        cthl_1=cthl_1,
        cthl_2=cthl_2,
        chi_1=chi_1,
        chi_2=chi_2,
        stdev_chi_1=stdev_chi_1,
        stdev_chi_2=stdev_chi_2,
        stdev_eta_1=stdev_eta_1,
        stdev_eta_2=stdev_eta_2,
        covar_chi_eta_1=covar_chi_eta_1,
        covar_chi_eta_2=covar_chi_eta_2,
        corr_w_chi_1=corr_w_chi_1,
        corr_w_chi_2=corr_w_chi_2,
        corr_w_eta_1=corr_w_eta_1,
        corr_w_eta_2=corr_w_eta_2,
        corr_chi_eta_1=corr_chi_eta_1,
        corr_chi_eta_2=corr_chi_eta_2,
        rsatl_1=rsatl_1,
        rsatl_2=rsatl_2,
        rc_1=rc_1,
        rc_2=rc_2,
        cloud_frac_1=cloud_frac_1,
        cloud_frac_2=cloud_frac_2,
        mixt_frac=mf,
        ice_supersat_frac_1=ice_supersat_frac_1,
        ice_supersat_frac_2=ice_supersat_frac_2,
    )

    return (
        pdf_params, pdf_implicit_coefs_terms, err_info,
        wpup2, wpvp2,
        wp2up2, wp2vp2, wp4,
        wprtp2, wp2rtp,
        wpthlp2, wp2thlp, wprtpthlp,
        cloud_frac, ice_supersat_frac,
        rcm, wpthvp, wp2thvp, wp2up,
        rtpthvp, thlpthvp, wprcp, wp2rcp,
        rtprcp, thlprcp, rcp2, uprcp, vprcp,
        w_up_in_cloud, w_down_in_cloud,
        cloudy_updraft_frac, cloudy_downdraft_frac,
        wpsclrprtp, wpsclrp2, sclrpthvp,
        wpsclrpthlp, sclrprcp, wp2sclrp,
        rc_coef,
    )


def transform_pdf_chi_eta_component(tl, rsatl, rt, exner_in,
                                        varnce_thl, varnce_rt, corr_rt_thl):
    """pdf_closure_module.F90:transform_pdf_chi_eta_component.

    s from Lewellen and Yoh 1993 (LY) eqn. 1

    For each normal distribution in the sum of two normal distributions,
    s' = crt * rt'  +  cthl * thl';
    therefore, x's' = crt * x'rt'  +  cthl * x'thl'.
    Larson et al. May, 2001.

    Calculate covariance, correlation, and standard deviation of
    chi and eta for each component
    Include subplume correlation of qt, thl

    Returns the Fortran out-arg order with JAX arrays:
        (chi, crt, cthl, stdev_chi, stdev_eta, covar_chi_eta, corr_chi_eta)
    """
    # SD's beta (eqn. 8)
    beta       = ep * Lv**2 / (Rd * Cp * tl**2)
    invrs      = 1.0 / (1.0 + beta * rsatl)
    # s from Lewellen and Yoh 1993 (LY) eqn. 1
    chi        = (rt - rsatl) * invrs
    # For each normal distribution in the sum of two normal distributions,
    # s' = crt * rt'  +  cthl * thl';
    # therefore, x's' = crt * x'rt'  +  cthl * x'thl'.
    # Larson et al. May, 2001.
    crt        = invrs
    cthl       = ((1.0 + beta * rt) * invrs**2
                  * (Cp / Lv) * beta * rsatl * exner_in)
    # Calculate covariance, correlation, and standard deviation of
    # chi and eta for each component
    # Include subplume correlation of qt, thl
    vrnc_rt_t  = crt**2 * varnce_rt
    vrnc_thl_t = cthl**2 * varnce_thl
    corr_t     = (2.0 * corr_rt_thl * crt * cthl
                  * jnp.sqrt(varnce_rt * varnce_thl))
    vrnc_chi   = vrnc_rt_t - corr_t + vrnc_thl_t
    vrnc_eta   = vrnc_rt_t + corr_t + vrnc_thl_t
    stdev_chi  = _safe_sqrt(vrnc_chi)
    stdev_eta  = _safe_sqrt(vrnc_eta)
    covar_chi_eta = vrnc_rt_t - vrnc_thl_t
    # smooth_corr_quotient (pdf_utilities.F90:1360)
    _denom_thresh = chi_tol * eta_tol
    _smth = min(min_max_smth_mag, _denom_thresh)
    denom = stdev_chi * stdev_eta
    tmp_d = 0.5 * (
        (jnp.abs(covar_chi_eta) / max_mag_correlation + denom)
        + jnp.sqrt(
            (jnp.abs(covar_chi_eta) / max_mag_correlation - denom) ** 2
            + _smth ** 2
        )
    )
    tmp_d = 0.5 * (
        (tmp_d + _denom_thresh)
        + jnp.sqrt((tmp_d - _denom_thresh) ** 2 + _smth ** 2)
    )
    corr_chi_eta = covar_chi_eta / tmp_d
    return chi, crt, cthl, stdev_chi, stdev_eta, covar_chi_eta, corr_chi_eta


def calc_wp4_pdf(wm, w_1, w_2, varnce_w_1, varnce_w_2, mixt_frac):
    """Description:
    Calculates <w'^4> by integrating over the PDF of w.  The integral is:

    <w'^4> = INT(-inf:inf) ( w - <w> )^4 P(w) dw;

    where <w> is the overall mean of w and P(w) is a two-component normal
    distribution of w.  The integrated equation is:

    <w'^4> = mixt_frac * ( 3 * sigma_w_1^4
                           + 6 * ( mu_w_1 - <w> )^2 * sigma_w_1^2
                           + ( mu_w_1 - <w> )^4 )
             + ( 1 - mixt_frac ) * ( 3 * sigma_w_2^4
                                     + 6 * ( mu_w_2 - <w> )^2 * sigma_w_2^2
                                     + ( mu_w_2 - <w> )^4 );

    where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    of w in the 2nd PDF component, sigma_w_1 is the standard deviation of w in
    the 1st PDF component, sigma_w_2 is the standard deviation of w in the 2nd
    PDF component, and mixt_frac is the mixture fraction, which is the weight
    of the 1st PDF component.
    """
    wm = jnp.asarray(wm); a = jnp.asarray(mixt_frac)
    d1 = jnp.asarray(w_1) - wm; d2 = jnp.asarray(w_2) - wm
    v1 = jnp.asarray(varnce_w_1); v2 = jnp.asarray(varnce_w_2)
    # Calculate <w'^4> by integrating over the PDF.
    return (a * (3.0 * v1 ** 2 + 6.0 * d1 ** 2 * v1 + d1 ** 4)
            + (1.0 - a) * (3.0 * v2 ** 2 + 6.0 * d2 ** 2 * v2 + d2 ** 4))


def calc_wp2xp2_pdf(wm, xm, w_1, w_2, x_1, x_2,
                         varnce_w_1, varnce_w_2,
                         varnce_x_1, varnce_x_2,
                         corr_w_x_1, corr_w_x_2,
                         mixt_frac):
    """Description:
    Calculates <w'^2x'^2> by integrating over the PDF of w and x.  The
    integral
    is:

    <w'^2x'^2>
    = INT(-inf:inf) INT(-inf:inf) ( w - <w> )^2 ( x - <x> )^2 P(w,x) dx dw;

    where <w> is the overall mean of w, <x> is the overall mean of x, and
    P(w,x) is a two-component bivariate normal distribution of w and x.  The
    integrated equation is:

    <w'^2x'^2>
      = mixt_frac
         * ( ( mu_w_1 - <w> )**2 * ( ( mu_x_1 - <x> )**2 + sigma_x_1^2 )
         + four * corr_w_x_1 * sigma_w_1 * sigma_x_1 * ( mu_x_1 - <x> ) * (
         mu_w_1 - <w> )
         + ( ( mu_x_1 - <x> )**2 + ( 1 + 2*corr_w_x_1**2 ) * sigma_x_1^2 ) *
         sigma_w_1^2 )
       + ( one - mixt_frac )
         * ( ( mu_w_2 - <w> )**2 * ( ( mu_x_2 - <x> )**2 + sigma_x_2^2 )
         + four * corr_w_x_2 * sigma_w_2 * sigma_x_2 * ( mu_x_2 - <x> ) * (
         mu_w_2 - <w> )
         + ( ( mu_x_2 - <x> )**2 + ( 1 + 2*corr_w_x_2**2 ) * sigma_x_2^2 ) *
         sigma_w_2^2 )

    where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
    component, mu_x_2 is the mean of x in the 2nd PDF component, sigma_w_1 is
    the standard deviation of w in the 1st PDF component, sigma_w_2 is the
    standard deviation of w in the 2nd PDF component, sigma_x_1 is the
    standard deviation of x in the 1st PDF component, sigma_x_2 is the
    standard deviation of x in the 2nd PDF component, corr_w_x_1 is the
    correlation of w and x in the 1st PDF component, corr_w_x_2 is the
    correlation of w and x in the 2nd PDF component, and mixt_frac is the
    mixture fraction, which is the weight of the 1st PDF component.
    """
    one_minus_mf = 1.0 - mixt_frac
    dw_1 = w_1 - wm
    dw_2 = w_2 - wm
    dx_1 = x_1 - xm
    dx_2 = x_2 - xm

    term1 = (dw_1 ** 2 * (dx_1 ** 2 + varnce_x_1)
             + 4.0 * corr_w_x_1 * _safe_sqrt(varnce_w_1 * varnce_x_1) * dx_1 * dw_1
             + (dx_1 ** 2 + (1.0 + 2.0 * corr_w_x_1 ** 2) * varnce_x_1) * varnce_w_1)
    term2 = (dw_2 ** 2 * (dx_2 ** 2 + varnce_x_2)
             + 4.0 * corr_w_x_2 * _safe_sqrt(varnce_w_2 * varnce_x_2) * dx_2 * dw_2
             + (dx_2 ** 2 + (1.0 + 2.0 * corr_w_x_2 ** 2) * varnce_x_2) * varnce_w_2)

    # Calculate <w'x'^2> by integrating over the PDF.
    return mixt_frac * term1 + one_minus_mf * term2


def calc_wp2xp_pdf(wm, xm, w_1, w_2, x_1, x_2, varnce_w_1, varnce_w_2,
                   varnce_x_1, varnce_x_2, corr_w_x_1, corr_w_x_2, mixt_frac):
    """Description:
    Calculates <w'^2 x'> by integrating over the PDF of w and x.  The integral
    is:

    <w'^2 x'>
    = INT(-inf:inf) INT(-inf:inf) ( w - <w> )^2 ( x - <x> ) P(w,x) dx dw;

    where <w> is the overall mean of w, <x> is the overall mean of x, and
    P(w,x) is a two-component bivariate normal distribution of w and x.  The
    integrated equation is:

    <w'^2 x'>
    = mixt_frac * ( ( mu_x_1 - <x> ) * ( ( mu_w_1 - <w> )^2 + sigma_w_1^2 )
                    + 2 * corr_w_x_1 * sigma_w_1 * sigma_x_1
                      * ( mu_w_1 - <w> ) )
      + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> )
                              * ( ( mu_w_2 - <w> )^2 + sigma_w_2^2 )
                              + 2 * corr_w_x_2 * sigma_w_2 * sigma_x_2
                                * ( mu_w_2 - <w> ) );

    where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
    component, mu_x_2 is the mean of x in the 2nd PDF component, sigma_w_1 is
    the standard deviation of w in the 1st PDF component, sigma_w_2 is the
    standard deviation of w in the 2nd PDF component, sigma_x_1 is the
    standard deviation of x in the 1st PDF component, sigma_x_2 is the
    standard deviation of x in the 2nd PDF component, corr_w_x_1 is the
    correlation of w and x in the 1st PDF component, corr_w_x_2 is the
    correlation of w and x in the 2nd PDF component, and mixt_frac is the
    mixture fraction, which is the weight of the 1st PDF component.
    """
    wm = jnp.asarray(wm); xm = jnp.asarray(xm); a = jnp.asarray(mixt_frac)
    dw1 = jnp.asarray(w_1) - wm; dw2 = jnp.asarray(w_2) - wm
    dx1 = jnp.asarray(x_1) - xm; dx2 = jnp.asarray(x_2) - xm
    vw1 = jnp.asarray(varnce_w_1); vw2 = jnp.asarray(varnce_w_2)
    vx1 = jnp.asarray(varnce_x_1); vx2 = jnp.asarray(varnce_x_2)
    c1 = jnp.asarray(corr_w_x_1); c2 = jnp.asarray(corr_w_x_2)
    # Calculate <w'^2 x'> by integrating over the PDF.
    return (a * ((dw1 ** 2 + vw1) * dx1 + 2.0 * c1 * _safe_sqrt(vw1 * vx1) * dw1)
            + (1.0 - a) * ((dw2 ** 2 + vw2) * dx2 + 2.0 * c2 * _safe_sqrt(vw2 * vx2) * dw2))


def calc_wpxp2_pdf(wm, xm, w_1, w_2, x_1, x_2, varnce_w_1, varnce_w_2,
                   varnce_x_1, varnce_x_2, corr_w_x_1, corr_w_x_2, mixt_frac):
    """Description:
    Calculates <w'x'^2> by integrating over the PDF of w and x.  The integral
    is:

    <w'x'^2>
    = INT(-inf:inf) INT(-inf:inf) ( w - <w> ) ( x - <x> )^2 P(w,x) dx dw;

    where <w> is the overall mean of w, <x> is the overall mean of x, and
    P(w,x) is a two-component bivariate normal distribution of w and x.  The
    integrated equation is:

    <w'x'^2>
    = mixt_frac * ( ( mu_w_1 - <w> ) * ( ( mu_x_1 - <x> )^2 + sigma_x_1^2 )
                    + 2 * corr_w_x_1 * sigma_w_1 * sigma_x_1
                      * ( mu_x_1 - <x> ) )
      + ( 1 - mixt_frac ) * ( ( mu_w_2 - <w> )
                              * ( ( mu_x_2 - <x> )^2 + sigma_x_2^2 )
                              + 2 * corr_w_x_2 * sigma_w_2 * sigma_x_2
                                * ( mu_x_2 - <x> ) );

    where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
    component, mu_x_2 is the mean of x in the 2nd PDF component, sigma_w_1 is
    the standard deviation of w in the 1st PDF component, sigma_w_2 is the
    standard deviation of w in the 2nd PDF component, sigma_x_1 is the
    standard deviation of x in the 1st PDF component, sigma_x_2 is the
    standard deviation of x in the 2nd PDF component, corr_w_x_1 is the
    correlation of w and x in the 1st PDF component, corr_w_x_2 is the
    correlation of w and x in the 2nd PDF component, and mixt_frac is the
    mixture fraction, which is the weight of the 1st PDF component.
    """
    wm = jnp.asarray(wm); xm = jnp.asarray(xm); a = jnp.asarray(mixt_frac)
    dw1 = jnp.asarray(w_1) - wm; dw2 = jnp.asarray(w_2) - wm
    dx1 = jnp.asarray(x_1) - xm; dx2 = jnp.asarray(x_2) - xm
    vw1 = jnp.asarray(varnce_w_1); vw2 = jnp.asarray(varnce_w_2)
    vx1 = jnp.asarray(varnce_x_1); vx2 = jnp.asarray(varnce_x_2)
    c1 = jnp.asarray(corr_w_x_1); c2 = jnp.asarray(corr_w_x_2)
    # Calculate <w'x'^2> by integrating over the PDF.
    return (a * (dw1 * (dx1 ** 2 + vx1) + 2.0 * c1 * _safe_sqrt(vw1 * vx1) * dx1)
            + (1.0 - a) * (dw2 * (dx2 ** 2 + vx2) + 2.0 * c2 * _safe_sqrt(vw2 * vx2) * dx2))


def calc_wpxpyp_pdf(wm, xm, ym, w_1, w_2, x_1, x_2, y_1, y_2,
                    varnce_w_1, varnce_w_2, varnce_x_1, varnce_x_2, varnce_y_1, varnce_y_2,
                    corr_w_x_1, corr_w_x_2, corr_w_y_1, corr_w_y_2, corr_x_y_1, corr_x_y_2, mixt_frac):
    """Description:
    Calculates <w'x'y'> by integrating over the PDF of w, x, and y.  The
    integral is:

    <w'x'y'>
    = INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
      ( w - <w> ) ( x - <x> ) ( y - <y> ) P(w,x,y) dy dx dw;

    where <w> is the overall mean of w, <x> is the overall mean of x, <y> is
    the overall mean of y, and P(w,x,y) is a two-component trivariate normal
    distribution of w, x, and y.  The integrated equation is:

    <w'x'y'>
    = mixt_frac
      * ( ( mu_w_1 - <w> ) * ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
          + corr_x_y_1 * sigma_x_1 * sigma_y_1 * ( mu_w_1 - <w> )
          + corr_w_y_1 * sigma_w_1 * sigma_y_1 * ( mu_x_1 - <x> )
          + corr_w_x_1 * sigma_w_1 * sigma_x_1 * ( mu_y_1 - <y> ) )
      + ( 1 - mixt_frac )
        * ( ( mu_w_2 - <w> ) * ( mu_x_2 - <x> ) * ( mu_y_2 - <y> )
            + corr_x_y_2 * sigma_x_2 * sigma_y_2 * ( mu_w_2 - <w> )
            + corr_w_y_2 * sigma_w_2 * sigma_y_2 * ( mu_x_2 - <x> )
            + corr_w_x_2 * sigma_w_2 * sigma_x_2 * ( mu_y_2 - <y> ) );

    where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
    component, mu_x_2 is the mean of x in the 2nd PDF component, mu_y_1 is the
    mean of y in the 1st PDF component, mu_y_2 is the mean of y in the 2nd PDF
    component, sigma_w_1 is the standard deviation of w in the 1st PDF
    component, sigma_w_2 is the standard deviation of w in the 2nd PDF
    component, sigma_x_1 is the standard deviation of x in the 1st PDF
    component, sigma_x_2 is the standard deviation of x in the 2nd PDF
    component, sigma_y_1 is the standard deviation of y in the 1st PDF
    component, sigma_y_2 is the standard deviation of y in the 2nd PDF
    component, corr_w_x_1 is the correlation of w and x in the 1st PDF
    component, corr_w_x_2 is the correlation of w and x in the 2nd PDF
    component, corr_w_y_1 is the correlation of w and y in the 1st PDF
    component, corr_w_y_2 is the correlation of w and y in the 2nd PDF
    component, corr_x_y_1 is the correlation of x and y in the 1st PDF
    component, corr_x_y_2 is the correlation of x and y in the 2nd PDF
    component, and mixt_frac is the mixture fraction, which is the weight of
    the 1st PDF component.
    """
    wm = jnp.asarray(wm); xm = jnp.asarray(xm); ym = jnp.asarray(ym); a = jnp.asarray(mixt_frac)
    dw1 = jnp.asarray(w_1) - wm; dw2 = jnp.asarray(w_2) - wm
    dx1 = jnp.asarray(x_1) - xm; dx2 = jnp.asarray(x_2) - xm
    dy1 = jnp.asarray(y_1) - ym; dy2 = jnp.asarray(y_2) - ym
    vw1 = jnp.asarray(varnce_w_1); vw2 = jnp.asarray(varnce_w_2)
    vx1 = jnp.asarray(varnce_x_1); vx2 = jnp.asarray(varnce_x_2)
    vy1 = jnp.asarray(varnce_y_1); vy2 = jnp.asarray(varnce_y_2)
    cwx1 = jnp.asarray(corr_w_x_1); cwx2 = jnp.asarray(corr_w_x_2)
    cwy1 = jnp.asarray(corr_w_y_1); cwy2 = jnp.asarray(corr_w_y_2)
    cxy1 = jnp.asarray(corr_x_y_1); cxy2 = jnp.asarray(corr_x_y_2)
    # Calculate <w'x'y'> by integrating over the PDF.
    comp1 = (dw1 * dx1 * dy1 + cxy1 * _safe_sqrt(vx1 * vy1) * dw1
             + cwy1 * _safe_sqrt(vw1 * vy1) * dx1 + cwx1 * _safe_sqrt(vw1 * vx1) * dy1)
    comp2 = (dw2 * dx2 * dy2 + cxy2 * _safe_sqrt(vx2 * vy2) * dw2
             + cwy2 * _safe_sqrt(vw2 * vy2) * dx2 + cwx2 * _safe_sqrt(vw2 * vx2) * dy2)
    return a * comp1 + (1.0 - a) * comp2


def calc_liquid_cloud_frac_component(mean_chi, stdev_chi_in):
    """Description:
    Calculates the PDF component cloud water mixing ratio, rc_i, and cloud
    fraction, cloud_frac_i, for the ith PDF component.

    The equation for cloud water mixing ratio, rc, at any point is:

    rc = chi * H(chi);

    and the equation for cloud fraction at a point, fc, is:

    fc = H(chi);

    where where extended liquid water mixing ratio, chi, is equal to cloud
    water mixing ratio, rc, when positive.  When the atmosphere is saturated
    at this point, cloud water is found, and rc = chi, while fc = 1.
    Otherwise, clear air is found at this point, and rc = fc = 0.

    The mean of rc and fc is calculated by integrating over the PDF, such
    that:

    <rc> = INT(-inf:inf) chi * H(chi) * P(chi) dchi; and

    cloud_frac = <fc> = INT(-inf:inf) H(chi) * P(chi) dchi.

    This can be rewritten as:

    <rc> = INT(0:inf) chi * P(chi) dchi; and

    cloud_frac = <fc> = INT(0:inf) P(chi) dchi;

    and further rewritten as:

    <rc> = SUM(i=1,N) mixt_frac_i INT(0:inf) chi * P_i(chi) dchi; and

    cloud_frac = SUM(i=1,N) mixt_frac_i INT(0:inf) P_i(chi) dchi;

    where N is the number of PDF components.  The equation for mean rc in the
    ith PDF component is:

    rc_i = INT(0:inf) chi * P_i(chi) dchi;

    and the equation for cloud fraction in the ith PDF component is:

    cloud_frac_i = INT(0:inf) P_i(chi) dchi.

    The component values are related to the overall values by:

    <rc> = SUM(i=1,N) mixt_frac_i * rc_i; and

    cloud_frac = SUM(i=1,N) mixt_frac_i * cloud_frac_i.
    """
    is_clear = (
        ((jnp.abs(mean_chi) <= eps) & (stdev_chi_in <= chi_tol))
        | (mean_chi < -max_num_stdevs * stdev_chi_in)
    )
    is_full  = mean_chi > max_num_stdevs * stdev_chi_in
    safe_s   = jnp.maximum(stdev_chi_in, 1.0e-100)
    zeta     = mean_chi / safe_s
    cf_mid   = 0.5 * (1.0 + jsp.erf(zeta / sqrt_2))
    rc_mid   = (mean_chi * cf_mid
                + stdev_chi_in * jnp.exp(-0.5 * zeta**2) / sqrt_2pi)
    # The mean of chi is at saturation and does not vary in the ith PDF component
    # The mean of chi is multiple standard deviations above the saturation point.
    # Thus, all cloud in the ith PDF component.
    # The mean of chi is within max_num_stdevs of the saturation point.
    # Thus, layer is partly cloudy, requires calculation.
    cf = jnp.where(is_clear, 0.0, jnp.where(is_full, 1.0, cf_mid))
    rc = jnp.where(is_clear, 0.0, jnp.where(is_full, mean_chi, rc_mid))
    return cf, rc


def calc_ice_cloud_frac_component(mean_chi, stdev_chi_in, crt, rsatl, tl,
                                      cf_liq, p_in_Pa):
    """Description:
    A version of the cloud fraction calculation modified to work
    for layers that are potentially below freezing. If there are
    no below freezing levels, the ice_supersat_frac calculation is
    the same as cloud_frac.

    For the below freezing levels, the saturation point will be
    non-zero, thus we need to calculate chi_at_ice_sat.

    The description of the equations are located in the description
    of calc_liquid_cloud_frac_component.
    """
    # Calculate the saturation mixing ratio of ice
    rsat_ice = sat_mixrat_ice(p_in_Pa, tl)
    # Temperature is freezing, we must compute chi_at_ice_sat and
    # calculate the new cloud_frac_component
    delta    = mean_chi - crt * (rsat_ice - rsatl)   # chi - chi_at_ice_sat
    is_clear = (((jnp.abs(delta) <= eps) & (stdev_chi_in <= chi_tol))
                | (delta < -max_num_stdevs * stdev_chi_in))
    is_full  = delta > max_num_stdevs * stdev_chi_in
    safe_s   = jnp.maximum(stdev_chi_in, 1.0e-100)
    zeta     = delta / safe_s
    ssf_mid  = 0.5 * (1.0 + jsp.erf(zeta / sqrt_2))
    # The mean of chi is at saturation and does not vary in the ith PDF component
    # The mean of chi is multiple standard deviations above the saturation point.
    # Thus, all cloud in the ith PDF component.
    # The mean of chi is within max_num_stdevs of the saturation point.
    # Thus, layer is partly cloudy, requires calculation.
    ssf      = jnp.where(is_clear, 0.0, jnp.where(is_full, 1.0, ssf_mid))
    # If a grid boxes is above freezing, then the calculation is the
    # same as the cloud_frac calculation
    return jnp.where(tl > T_freeze_K, cf_liq, ssf)


def calc_xprcp_component(wm, rtm, thlm, um, vm, rcm,
                             w_i, rt_i, thl_i, u_i, v_i,
                             varnce_w_i, chi_i, stdev_chi_i, stdev_eta_i,
                             corr_w_chi_i, corr_chi_eta_i, crt_i, cthl_i,
                             rc_i, cloud_frac_i, iiPDF_type):
    """Description:
    Calculates the contribution to <w'rc'>, <w'^2 rc'>, <rt'rc'>, and
    <thl'rc'> from the ith PDF component.

    Use equations for PDF component cloud fraction cloud water mixing ratio
    -----------------------------------------------------------------------

    The equation for cloud fraction in the ith PDF component, fc_i, is:

    fc_i = 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) ).

    In the special case that sigma_chi_i = 0, the equation becomes:

    fc_i = | 1; when mu_chi_i > 0;
           | 0; when mu_chi_i <= 0.

    The equation for mean cloud water mixing ratio in the ith PDF component,
    rc_i, is:

    rc_i
    = mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
      + 1/sqrt(2*pi) * sigma_chi_i
        * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) }
    = mu_chi_i * fc_i
      + 1/sqrt(2*pi) * sigma_chi_i
        * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) }.

    The above equations can be substituted into the equations for
    wprcp_contrib_comp_i, wp2rcp_contrib_comp_i, rtprcp_contrib_comp_i, and
    thlprcp_contrib_comp_i.  The new equations are:

    wprcp_contrib_comp_i
    = ( mu_w_i - <w> ) * ( rc_i - <rc> )
      + corr_w_chi_i * sigma_w_i * sigma_chi_i * fc_i;

    wp2rcp_contrib_comp_i
    = ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * ( rc_i - <rc> )
      + 2 * ( mu_w_i - <w> ) * corr_w_chi_i * sigma_w_i * sigma_chi_i * fc_i
      + 1/sqrt(2*pi) * corr_w_chi_i^2 * sigma_w_i^2 * sigma_chi_i
        * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) };

    rtprcp_contrib_comp_i
    = ( mu_rt_i - <rt> ) * ( rc_i - <rc> )
      + ( corr_chi_eta_i * sigma_eta_i + sigma_chi_i ) / ( 2 * crt_i )
        * sigma_chi_i * fc_i; and

    thlprcp_contrib_comp_i
    = ( mu_thl_i - <thl> ) * ( rc_i - <rc> )
      + ( corr_chi_eta_i * sigma_eta_i - sigma_chi_i ) / ( 2 * cthl_i )
        * sigma_chi_i * fc_i.

    JAX port note: this routine returns
    (wprcp, wp2rcp, rtprcp, thlprcp, uprcp, vprcp). The full derivation in
    Fortran is much longer; the active code below corresponds to the compact
    equations above plus the non-ADG corr_w_chi correction.
    """
    drc = rc_i - rcm
    wprcp  = (w_i - wm) * drc
    wp2rcp = ((w_i - wm) ** 2 + varnce_w_i) * drc
    rtprcp = ((rt_i - rtm) * drc
              + (corr_chi_eta_i * stdev_eta_i + stdev_chi_i)
                / (2.0 * crt_i) * stdev_chi_i * cloud_frac_i)
    # Guard against cthl=0 (rsatl=0 limit); cloud_frac=0 masks the result there.
    cthl_safe = jnp.where(cthl_i == 0.0, 1.0, cthl_i)
    thlprcp = ((thl_i - thlm) * drc
               + (corr_chi_eta_i * stdev_eta_i - stdev_chi_i)
                 / (2.0 * cthl_safe) * stdev_chi_i * cloud_frac_i)
    uprcp = (u_i - um) * drc
    vprcp = (v_i - vm) * drc
    # If iiPDF_type isn't iiPDF_ADG1, iiPDF_ADG2, or iiPDF_new_hybrid, so
    # corr_w_chi_i /= 0 (and perhaps corr_u_w_i /= 0).
    if iiPDF_type not in (iiPDF_ADG1, iiPDF_ADG2, iiPDF_new_hybrid):
        # Chi varies significantly in the ith PDF component (stdev_chi > chi_tol)
        # and there is some cloud (0 < cloud_frac <= 1)
        active = (stdev_chi_i > chi_tol) & (cloud_frac_i > 0.0)
        stdev_chi_safe = jnp.where(active, stdev_chi_i, 1.0)
        wprcp = wprcp + jnp.where(
            active,
            corr_w_chi_i * _safe_sqrt(varnce_w_i) * stdev_chi_i * cloud_frac_i,
            0.0,
        )
        wp2rcp = wp2rcp + jnp.where(
            active,
            2.0 * (w_i - wm) * corr_w_chi_i
            * _safe_sqrt(varnce_w_i) * stdev_chi_i * cloud_frac_i
            + corr_w_chi_i ** 2 * varnce_w_i * stdev_chi_i
            * jnp.exp(-chi_i ** 2 / (2.0 * stdev_chi_safe ** 2)) / sqrt_2pi,
            0.0,
        )
    return wprcp, wp2rcp, rtprcp, thlprcp, uprcp, vprcp


def calc_w_up_in_cloud(mixt_frac, cloud_frac_1, cloud_frac_2,
                       w_1, w_2, varnce_w_1, varnce_w_2):
    """Description:
    Subroutine that computes the mean cloudy updraft (and also calculates
    the mean cloudy downdraft).

    In order to activate aerosol, we'd like to feed the activation scheme
    a vertical velocity that's representative of cloudy updrafts. For skewed
    layers, like cumulus layers, this might be an improvement over the square
    root of wp2 that's currently used. At the same time, it would be simpler
    and less expensive than feeding SILHS samples into the aerosol code
    (see larson-group/e3sm#19 and larson-group/e3sm#26).

    The formulas are only valid for certain PDFs in CLUBB (ADG1, ADG2,
    new hybrid), hence we omit calculation if another PDF type is used.

    References: https://www.overleaf.com/project/614a136d47846639af22ae34

    JAX port note: the caller only invokes this for supported PDF types, and
    this function returns
    (w_up_in_cloud, w_down_in_cloud, cloudy_updraft_frac, cloudy_downdraft_frac).
    """
    def _component(w, varnce):
        w = jnp.asarray(w, dtype=jnp.float64)
        stdev = jnp.sqrt(jnp.asarray(varnce, dtype=jnp.float64))
        all_up = w > max_num_stdevs * stdev
        all_down = w < -max_num_stdevs * stdev
        ratio = w / (sqrt_2 * jnp.maximum(eps, stdev))
        erf_r = jsp.erf(ratio)
        exp_neg = jnp.exp(-ratio ** 2)
        w_up_mid = 0.5 * w * (1.0 + erf_r) + (stdev / sqrt_2pi) * exp_neg
        uf_mid = 0.5 * (1.0 + erf_r)
        w_down_mid = 0.5 * w * (1.0 - erf_r) - (stdev / sqrt_2pi) * exp_neg
        # The mean of w in the PDF component is more than
        # max_num_stdevs standard deviations above 0.
        # The entire PDF component is found in an updraft (w > 0).
        w_up = jnp.where(all_up, w, jnp.where(all_down, 0.0, w_up_mid))
        uf = jnp.where(all_up, 1.0, jnp.where(all_down, 0.0, uf_mid))
        # The mean of w in the PDF component is more than
        # max_num_stdevs standard deviations below 0.
        # The entire PDF component is found in a downdraft (w < 0).
        # Otherwise, the PDF component contains both updraft and downdraft.
        w_down = jnp.where(all_up, 0.0, jnp.where(all_down, w, w_down_mid))
        df = 1.0 - uf   # holds in all three branches (Fortran: 1, 0, 1-uf_mid)
        return w_up, uf, w_down, df

    a = jnp.asarray(mixt_frac, dtype=jnp.float64)
    cf1 = jnp.asarray(cloud_frac_1, dtype=jnp.float64); cf2 = jnp.asarray(cloud_frac_2, dtype=jnp.float64)
    w_up_1, uf_1, w_down_1, df_1 = _component(w_1, varnce_w_1)
    w_up_2, uf_2, w_down_2, df_2 = _component(w_2, varnce_w_2)

    # Calculate the total cloudy updraft fraction.
    cloudy_updraft_frac = a * cf1 * uf_1 + (1.0 - a) * cf2 * uf_2
    # Calculate the total cloudy downdraft fraction.
    cloudy_downdraft_frac = a * cf1 * df_1 + (1.0 - a) * cf2 * df_2
    # Calculate the mean vertical velocity found in a cloudy updraft.
    w_up_in_cloud = ((a * cf1 * w_up_1 + (1.0 - a) * cf2 * w_up_2)
                     / jnp.maximum(eps, cloudy_updraft_frac))
    # Calculate the mean vertical velocity found in a cloudy downdraft.
    w_down_in_cloud = ((a * cf1 * w_down_1 + (1.0 - a) * cf2 * w_down_2)
                       / jnp.maximum(eps, cloudy_downdraft_frac))
    return w_up_in_cloud, w_down_in_cloud, cloudy_updraft_frac, cloudy_downdraft_frac


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "nzt",
        "ngrdcol",
        "hydromet_dim",
        "sclr_dim",
        "l_samp_stats_in_pdf_call",
        "iiPDF_type",
        "saturation_formula",
        "l_rtm_nudge",
        "l_trapezoidal_rule_zt",
        "l_trapezoidal_rule_zm",
        "l_call_pdf_closure_twice",
        "l_use_cloud_cover",
        "l_rcm_supersat_adj",
    ),
)


def pdf_closure_driver(
    gr, nzm, nzt, ngrdcol, dt, hydromet_dim, sclr_dim,
    sclr_tol, wprtp, thlm, wpthlp, rtp2, rtp3,
    thlp2, thlp3, rtpthlp, wp2,
    wp3, wm_zm, wm_zt,
    um, up2, upwp, up3,
    vm, vp2, vpwp, vp3,
    p_in_Pa, exner,
    thv_ds_zm, thv_ds_zt, rtm_ref,
    wphydrometp,
    wp2hmp, rtphmp_zt, thlphmp_zt,
    sclrm, wpsclrp, sclrp2,
    sclrprtp, sclrpthlp, sclrp3,
    p_sfc, l_samp_stats_in_pdf_call, stats,
    mixt_frac_max_mag, ts_nudge,
    rtm_min, rtm_nudge_max_altitude,
    clubb_params,
    iiPDF_type,
    saturation_formula,
    l_rtm_nudge,
    l_trapezoidal_rule_zt,
    l_trapezoidal_rule_zm,
    l_call_pdf_closure_twice,
    l_use_cloud_cover,
    l_rcm_supersat_adj,
    l_mix_rat_hm,
    rtm, sigma_sqd_w,
    pdf_implicit_coefs_terms,
    pdf_params,
    pdf_params_zm,
    err_info,
):
    """JAX port of pdf_closure_module.F90:pdf_closure_driver."""
    if iiPDF_type not in (
        iiPDF_ADG1,
        iiPDF_ADG2,
        iiPDF_3D_Luhar,
        iiPDF_new,
        iiPDF_TSDADG,
        iiPDF_LY93,
        iiPDF_new_hybrid,
    ):
        raise NotImplementedError(_unsupported_pdf_type_message(iiPDF_type))

    wp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, wp2, w_tol_sqd)
    rtp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, rtp2, rt_tol ** 2)
    thlp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, thlp2, thl_tol ** 2)
    wprtp_zt = zm2zt(nzm, nzt, ngrdcol, gr, wprtp)
    wpthlp_zt = zm2zt(nzm, nzt, ngrdcol, gr, wpthlp)
    rtpthlp_zt = zm2zt(nzm, nzt, ngrdcol, gr, rtpthlp)
    up2_zt = zm2zt(nzm, nzt, ngrdcol, gr, up2, w_tol_sqd)
    vp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, vp2, w_tol_sqd)
    upwp_zt = zm2zt(nzm, nzt, ngrdcol, gr, upwp)
    vpwp_zt = zm2zt(nzm, nzt, ngrdcol, gr, vpwp)
    sigma_sqd_w_zt = zm2zt(nzm, nzt, ngrdcol, gr, sigma_sqd_w, zero_threshold)

    wp3_zm = zt2zm(nzm, nzt, ngrdcol, gr, wp3)
    rtp3_zm = zt2zm(nzm, nzt, ngrdcol, gr, rtp3)
    thlp3_zm = zt2zm(nzm, nzt, ngrdcol, gr, thlp3)

    Skw_zt = Skx_func(nzt, ngrdcol, wp2_zt, wp3, w_tol, clubb_params)
    Skw_zm = Skx_func(nzm, ngrdcol, wp2, wp3_zm, w_tol, clubb_params)
    Skrt_zt = Skx_func(nzt, ngrdcol, rtp2_zt, rtp3, rt_tol, clubb_params)
    Skrt_zm = Skx_func(nzm, ngrdcol, rtp2, rtp3_zm, rt_tol, clubb_params)
    Skthl_zt = Skx_func(nzt, ngrdcol, thlp2_zt, thlp3, thl_tol, clubb_params)
    Skthl_zm = Skx_func(nzm, ngrdcol, thlp2, thlp3_zm, thl_tol, clubb_params)
    Sku_zt = Skx_func(nzt, ngrdcol, up2_zt, up3, w_tol, clubb_params)
    Skv_zt = Skx_func(nzt, ngrdcol, vp2_zt, vp3, w_tol, clubb_params)

    if sclr_dim > 0:
        wpsclrp_zt = jnp.moveaxis(
            jax.vmap(lambda field: zm2zt(nzm, nzt, ngrdcol, gr, field), in_axes=2)(
                wpsclrp,
            ),
            0,
            -1,
        )
        sclrprtp_zt = jnp.moveaxis(
            jax.vmap(lambda field: zm2zt(nzm, nzt, ngrdcol, gr, field), in_axes=2)(
                sclrprtp,
            ),
            0,
            -1,
        )
        sclrpthlp_zt = jnp.moveaxis(
            jax.vmap(lambda field: zm2zt(nzm, nzt, ngrdcol, gr, field), in_axes=2)(
                sclrpthlp,
            ),
            0,
            -1,
        )
        sclrp2_zt = jnp.moveaxis(
            jax.vmap(
                lambda field, tol: zm2zt(nzm, nzt, ngrdcol, gr, field, tol ** 2),
                in_axes=(2, 0),
            )(sclrp2, sclr_tol),
            0,
            -1,
        )
        Sksclr_zt = jnp.moveaxis(
            jax.vmap(
                lambda sclrp2_s, sclrp3_s, tol: Skx_func(
                    nzt, ngrdcol, sclrp2_s, sclrp3_s, tol, clubb_params,
                ),
                in_axes=(2, 2, 0),
            )(sclrp2_zt, sclrp3, sclr_tol),
            0,
            -1,
        )
    else:
        wpsclrp_zt = jnp.zeros((ngrdcol, nzt, 0), dtype=jnp.float64)
        sclrprtp_zt = jnp.zeros((ngrdcol, nzt, 0), dtype=jnp.float64)
        sclrpthlp_zt = jnp.zeros((ngrdcol, nzt, 0), dtype=jnp.float64)
        sclrp2_zt = jnp.zeros((ngrdcol, nzt, 0), dtype=jnp.float64)
        Sksclr_zt = jnp.zeros((ngrdcol, nzt, 0), dtype=jnp.float64)

    # Interpolate hydrometeor mixed moments to momentum levels.
    if hydromet_dim > 0:
        wphydrometp_zt = jnp.moveaxis(
            jax.vmap(lambda field: zm2zt(nzm, nzt, ngrdcol, gr, field), in_axes=2)(
                wphydrometp,
            ),
            0,
            -1,
        )
    else:
        wphydrometp_zt = jnp.zeros((ngrdcol, nzt, 0), dtype=jnp.float64)

    #----------------------------------------------------------------
    # Call closure scheme
    #----------------------------------------------------------------
    (
        pdf_params, pdf_implicit_coefs_terms, err_info,
        wpup2, wpvp2,
        wp2up2_zt, wp2vp2_zt, wp4_zt,
        wprtp2, wp2rtp,
        wpthlp2, wp2thlp, wprtpthlp,
        cloud_frac, ice_supersat_frac,
        rcm, wpthvp_zt, wp2thvp, wp2up,
        rtpthvp_zt, thlpthvp_zt, wprcp_zt, wp2rcp,
        rtprcp_zt, thlprcp_zt, rcp2_zt, uprcp_zt, vprcp_zt,
        w_up_in_cloud, w_down_in_cloud,
        cloudy_updraft_frac, cloudy_downdraft_frac,
        wpsclrprtp, wpsclrp2, sclrpthvp_zt,
        wpsclrpthlp, sclrprcp_zt, wp2sclrp,
        rc_coef,
    ) = pdf_closure(
        nzt, ngrdcol, sclr_dim, sclr_tol, gr,
        hydromet_dim, p_in_Pa, exner, thv_ds_zt,
        wm_zt, wp2_zt, wp3,
        Skw_zt, Skthl_zt, Skrt_zt, Sku_zt, Skv_zt,
        rtm, rtp2_zt, wprtp_zt,
        thlm, thlp2_zt, wpthlp_zt,
        um, up2_zt, upwp_zt,
        vm, vp2_zt, vpwp_zt,
        rtpthlp_zt,
        sclrm, wpsclrp_zt, sclrp2_zt,
        sclrprtp_zt, sclrpthlp_zt, Sksclr_zt,
        wphydrometp_zt, wp2hmp,
        rtphmp_zt, thlphmp_zt,
        clubb_params, mixt_frac_max_mag,
        saturation_formula,
        stats,
        iiPDF_type,
        l_mix_rat_hm,
        sigma_sqd_w_zt,
        pdf_params, pdf_implicit_coefs_terms,
        err_info,
    )
    del wphydrometp_zt

    zero_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    zero_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)

    wp4 = zt2zm(nzm, nzt, ngrdcol, gr, wp4_zt, zero_threshold)
    wp4 = wp4.at[:, gr.k_lb_zm].set(0.0).at[:, gr.k_ub_zm].set(0.0)

    wpthvp = zt2zm(nzm, nzt, ngrdcol, gr, wpthvp_zt).at[:, gr.k_ub_zm].set(0.0)
    thlpthvp = zt2zm(nzm, nzt, ngrdcol, gr, thlpthvp_zt).at[:, gr.k_ub_zm].set(0.0)
    rtpthvp = zt2zm(nzm, nzt, ngrdcol, gr, rtpthvp_zt).at[:, gr.k_ub_zm].set(0.0)
    wprcp = zt2zm(nzm, nzt, ngrdcol, gr, wprcp_zt).at[:, gr.k_ub_zm].set(0.0)
    rc_coef_zm = zt2zm(nzm, nzt, ngrdcol, gr, rc_coef).at[:, gr.k_ub_zm].set(0.0)
    rtprcp = zt2zm(nzm, nzt, ngrdcol, gr, rtprcp_zt).at[:, gr.k_ub_zm].set(0.0)
    thlprcp = zt2zm(nzm, nzt, ngrdcol, gr, thlprcp_zt).at[:, gr.k_ub_zm].set(0.0)
    uprcp = zt2zm(nzm, nzt, ngrdcol, gr, uprcp_zt).at[:, gr.k_ub_zm].set(0.0)
    vprcp = zt2zm(nzm, nzt, ngrdcol, gr, vprcp_zt).at[:, gr.k_ub_zm].set(0.0)
    wp2up2 = zt2zm(nzm, nzt, ngrdcol, gr, wp2up2_zt).at[:, gr.k_ub_zm].set(0.0)
    wp2vp2 = zt2zm(nzm, nzt, ngrdcol, gr, wp2vp2_zt).at[:, gr.k_ub_zm].set(0.0)

    rcp2 = zt2zm(nzm, nzt, ngrdcol, gr, rcp2_zt, zero_threshold)
    rcp2 = rcp2.at[:, gr.k_ub_zm].set(0.0)

    if sclr_dim > 0:
        sclrpthvp = jnp.moveaxis(
            jax.vmap(
                lambda field: zt2zm(nzm, nzt, ngrdcol, gr, field)
                .at[:, gr.k_ub_zm].set(0.0),
                in_axes=2,
            )(sclrpthvp_zt),
            0,
            -1,
        )
        sclrprcp = jnp.moveaxis(
            jax.vmap(
                lambda field: zt2zm(nzm, nzt, ngrdcol, gr, field)
                .at[:, gr.k_ub_zm].set(0.0),
                in_axes=2,
            )(sclrprcp_zt),
            0,
            -1,
        )
    else:
        sclrpthvp = jnp.zeros_like(wpsclrp)
        sclrprcp = jnp.zeros_like(wpsclrp)

    cloud_frac_zm = zt2zm(nzm, nzt, ngrdcol, gr, cloud_frac)
    cloud_frac_zm = cloud_frac_zm.at[:, gr.k_ub_zm].set(0.0)
    ice_supersat_frac_zm = zero_zm
    rtm_zm = zero_zm
    thlm_zm = zero_zm
    rcm_zm = zero_zm
    wprtp2_zm = zero_zm
    wp2rtp_zm = zero_zm
    wpthlp2_zm = zero_zm
    wp2thlp_zm = zero_zm
    wprtpthlp_zm = zero_zm
    wp2thvp_zm = zero_zm
    wp2up_zm = zero_zm
    wp2rcp_zm = zero_zm
    wpsclrprtp_zm = jnp.zeros((ngrdcol, nzm, sclr_dim), dtype=jnp.float64)
    wpsclrp2_zm = jnp.zeros((ngrdcol, nzm, sclr_dim), dtype=jnp.float64)
    wpsclrpthlp_zm = jnp.zeros((ngrdcol, nzm, sclr_dim), dtype=jnp.float64)
    wp2sclrp_zm = jnp.zeros((ngrdcol, nzm, sclr_dim), dtype=jnp.float64)

    if l_call_pdf_closure_twice:
        if l_rtm_nudge:
            rtm = jnp.where(
                (rtm < rtm_min) & (gr.zt < rtm_nudge_max_altitude),
                rtm + (rtm_ref - rtm) * (dt / ts_nudge),
                rtm,
            )

        (
            pdf_params_zm, err_info,
            rtm_zm, thlm_zm,
            wp2up2, wp2vp2, wp4,
            wprtp2_zm, wp2rtp_zm,
            wpthlp2_zm, wp2thlp_zm, wprtpthlp_zm,
            cloud_frac_zm, ice_supersat_frac_zm,
            rcm_zm, wpthvp, wp2thvp_zm, wp2up_zm,
            rtpthvp, thlpthvp, wprcp, wp2rcp_zm,
            rtprcp, thlprcp, rcp2, uprcp, vprcp,
            wpsclrprtp_zm, wpsclrp2_zm, sclrpthvp,
            wpsclrpthlp_zm, sclrprcp, wp2sclrp_zm,
            rc_coef_zm,
        ) = pdf_closure_driver_zm(
            gr, nzm, nzt, ngrdcol, hydromet_dim,
            sclr_dim, sclr_tol, p_sfc,
            mixt_frac_max_mag, clubb_params,
            iiPDF_type, saturation_formula,
            l_mix_rat_hm, p_in_Pa, thv_ds_zm,
            wm_zm, wp2, wp3_zm,
            Skw_zm, Skthl_zm, Skrt_zm,
            rtm, rtp2, wprtp,
            thlm, thlp2, wpthlp,
            um, up2, upwp, up3,
            vm, vp2, vpwp, vp3,
            rtpthlp, sclrm, wpsclrp, sclrp2,
            sclrprtp, sclrpthlp, sclrp3,
            wphydrometp,
            wp2hmp, rtphmp_zt, thlphmp_zt,
            stats, sigma_sqd_w, pdf_params_zm,
            err_info,
        )

        if not l_trapezoidal_rule_zt:
            cloud_frac_zm = zt2zm(nzm, nzt, ngrdcol, gr, cloud_frac)
            cloud_frac_zm = cloud_frac_zm.at[:, gr.k_ub_zm].set(0.0)

    if l_trapezoidal_rule_zt:
        (
            wprtp2, wpthlp2, wprtpthlp, cloud_frac, ice_supersat_frac,
            rcm, wp2thvp, wp2up, wpsclrprtp, wpsclrp2, wpsclrpthlp,
            wprtp2_zm, wpthlp2_zm, wprtpthlp_zm, cloud_frac_zm,
            ice_supersat_frac_zm, rcm_zm, wp2thvp_zm, wp2up_zm,
            wpsclrprtp_zm, wpsclrp2_zm, wpsclrpthlp_zm,
        ) = trapezoidal_rule_zt(
            nzm, nzt, ngrdcol, sclr_dim, gr,
            l_call_pdf_closure_twice, stats,
            wprtp2, wpthlp2, wprtpthlp, cloud_frac, ice_supersat_frac,
            rcm, wp2thvp, wp2up, wpsclrprtp, wpsclrp2, wpsclrpthlp,
            wprtp2_zm, wpthlp2_zm, wprtpthlp_zm, cloud_frac_zm,
            ice_supersat_frac_zm, rcm_zm, wp2thvp_zm, wp2up_zm,
            wpsclrprtp_zm, wpsclrp2_zm, wpsclrpthlp_zm,
        )
        del wprtp2_zm, wp2rtp_zm, wpthlp2_zm, wp2thlp_zm
        del wprtpthlp_zm, wp2thvp_zm, wp2up_zm, wp2rcp_zm
        del wpsclrprtp_zm, wpsclrp2_zm, wpsclrpthlp_zm, wp2sclrp_zm

    if l_trapezoidal_rule_zm:
        wpthvp, thlpthvp, rtpthvp = trapezoidal_rule_zm(
            nzm, nzt, ngrdcol, gr,
            wpthvp_zt, thlpthvp_zt, rtpthvp_zt,
            wpthvp, thlpthvp, rtpthvp,
        )

    rcm = clip_rcm(nzt, ngrdcol, rtm, "rtm < rcm after pdf_closure", rcm)
    cloud_cover, rcm_in_layer = compute_cloud_cover(
        gr, nzt, ngrdcol, pdf_params, cloud_frac, rcm,
    )
    if l_use_cloud_cover:
        cloud_frac = cloud_cover
        rcm = rcm_in_layer
    cloud_frac = jnp.minimum(1.0, cloud_frac)
    ice_supersat_frac = jnp.minimum(1.0, ice_supersat_frac)

    if l_rtm_nudge and not l_call_pdf_closure_twice:
        rtm = jnp.where(
            (rtm < rtm_min) & (gr.zt < rtm_nudge_max_altitude),
            rtm + (rtm_ref - rtm) * (dt / ts_nudge),
            rtm,
        )

    if l_rcm_supersat_adj:
        from clubb_jax.src.CLUBB_core.T_in_K_module import thlm2T_in_K

        T_in_K = thlm2T_in_K(thlm, exner, rcm)
        rsat = sat_mixrat_liq(p_in_Pa, T_in_K, saturation_formula)
        rel_humidity = (rtm - rcm) / rsat
        rcm_supersat_adj = jnp.where(
            rel_humidity > 1.0,
            (rtm - rcm) - rsat,
            zero_zt,
        )
        rcm = rcm + rcm_supersat_adj
    else:
        rcm_supersat_adj = zero_zt

    Skw_velocity = (
        (1.0 / (1.0 - sigma_sqd_w))
        * (wp3_zm / jnp.maximum(wp2, w_tol_sqd))
    )
    skw_velocity = Skw_velocity

    if l_samp_stats_in_pdf_call:
        for name, value in (
            ("Skw_zt", Skw_zt),
            ("Skw_zm", Skw_zm),
            ("Skthl_zt", Skthl_zt),
            ("Skthl_zm", Skthl_zm),
            ("Skrt_zt", Skrt_zt),
            ("Skrt_zm", Skrt_zm),
            ("uprcp", uprcp),
            ("vprcp", vprcp),
            ("rc_coef", rc_coef),
            ("wp2rcp", wp2rcp),
            ("rtprcp", rtprcp),
            ("rcp2", rcp2),
            ("Skw_velocity", Skw_velocity),
            ("cloud_frac", cloud_frac),
            ("ice_supersat_frac", ice_supersat_frac),
            ("rcm_in_layer", rcm_in_layer),
            ("cloud_cover", cloud_cover),
            ("sigma_sqd_w", sigma_sqd_w),
            ("wpthvp", wpthvp),
            ("wp2thvp", wp2thvp),
            ("wp2up", wp2up),
            ("rtpthvp", rtpthvp),
            ("thlpthvp", thlpthvp),
            ("wprcp", wprcp),
            ("rc_coef_zm", rc_coef_zm),
            ("thlprcp", thlprcp),
            ("wpup2", wpup2),
            ("wpvp2", wpvp2),
            ("wp2up2", wp2up2),
            ("wp2vp2", wp2vp2),
            ("wp4", wp4),
            ("wp2rtp", wp2rtp),
            ("wprtp2", wprtp2),
            ("wp2thlp", wp2thlp),
            ("wpthlp2", wpthlp2),
            ("wprtpthlp", wprtpthlp),
            ("w_up_in_cloud", w_up_in_cloud),
            ("w_down_in_cloud", w_down_in_cloud),
            ("cld_updr_frac", cloudy_updraft_frac),
            ("cld_downdr_frac", cloudy_downdraft_frac),
            ("cloud_frac_zm", cloud_frac_zm),
            ("ice_supersat_frac_zm", ice_supersat_frac_zm),
            ("rtm_zm", rtm_zm),
            ("thlm_zm", thlm_zm),
            ("rcm_zm", rcm_zm),
            ("rcm_supersat_adj", rcm_supersat_adj),
            ("mixt_frac", pdf_params.mixt_frac),
            ("w_1", pdf_params.w_1),
            ("w_2", pdf_params.w_2),
            ("varnce_w_1", pdf_params.varnce_w_1),
            ("varnce_w_2", pdf_params.varnce_w_2),
            ("thl_1", pdf_params.thl_1),
            ("thl_2", pdf_params.thl_2),
            ("varnce_thl_1", pdf_params.varnce_thl_1),
            ("varnce_thl_2", pdf_params.varnce_thl_2),
            ("rt_1", pdf_params.rt_1),
            ("rt_2", pdf_params.rt_2),
            ("varnce_rt_1", pdf_params.varnce_rt_1),
            ("varnce_rt_2", pdf_params.varnce_rt_2),
            ("rc_1", pdf_params.rc_1),
            ("rc_2", pdf_params.rc_2),
            ("rsatl_1", pdf_params.rsatl_1),
            ("rsatl_2", pdf_params.rsatl_2),
            ("cloud_frac_1", pdf_params.cloud_frac_1),
            ("cloud_frac_2", pdf_params.cloud_frac_2),
            ("chi_1", pdf_params.chi_1),
            ("chi_2", pdf_params.chi_2),
            ("stdev_chi_1", pdf_params.stdev_chi_1),
            ("stdev_chi_2", pdf_params.stdev_chi_2),
            ("stdev_eta_1", pdf_params.stdev_eta_1),
            ("stdev_eta_2", pdf_params.stdev_eta_2),
            ("covar_chi_eta_1", pdf_params.covar_chi_eta_1),
            ("covar_chi_eta_2", pdf_params.covar_chi_eta_2),
            ("corr_w_chi_1", pdf_params.corr_w_chi_1),
            ("corr_w_chi_2", pdf_params.corr_w_chi_2),
            ("corr_w_eta_1", pdf_params.corr_w_eta_1),
            ("corr_w_eta_2", pdf_params.corr_w_eta_2),
            ("corr_chi_eta_1", pdf_params.corr_chi_eta_1),
            ("corr_chi_eta_2", pdf_params.corr_chi_eta_2),
            ("corr_w_rt_1", pdf_params.corr_w_rt_1),
            ("corr_w_rt_2", pdf_params.corr_w_rt_2),
            ("corr_w_thl_1", pdf_params.corr_w_thl_1),
            ("corr_w_thl_2", pdf_params.corr_w_thl_2),
            ("corr_rt_thl_1", pdf_params.corr_rt_thl_1),
            ("corr_rt_thl_2", pdf_params.corr_rt_thl_2),
            ("crt_1", pdf_params.crt_1),
            ("crt_2", pdf_params.crt_2),
            ("cthl_1", pdf_params.cthl_1),
            ("cthl_2", pdf_params.cthl_2),
        ):
            stats = stats.update(name, value)

        if stats.var_on_stats_list("chi") or stats.var_on_stats_list("chip2"):
            chi = (
                pdf_params.mixt_frac * pdf_params.chi_1
                + (1.0 - pdf_params.mixt_frac) * pdf_params.chi_2
            )
            stats = stats.update("chi", chi)

        if stats.var_on_stats_list("chip2"):
            chip2 = (
                pdf_params.mixt_frac
                * ((pdf_params.chi_1 - chi) ** 2 + pdf_params.stdev_chi_1 ** 2)
                + (1.0 - pdf_params.mixt_frac)
                * ((pdf_params.chi_2 - chi) ** 2 + pdf_params.stdev_chi_2 ** 2)
            )
            stats = stats.update("chip2", chip2)

        if l_call_pdf_closure_twice:
            for name, value in (
                ("w_1_zm", pdf_params_zm.w_1),
                ("w_2_zm", pdf_params_zm.w_2),
                ("varnce_w_1_zm", pdf_params_zm.varnce_w_1),
                ("varnce_w_2_zm", pdf_params_zm.varnce_w_2),
                ("mixt_frac_zm", pdf_params_zm.mixt_frac),
            ):
                stats = stats.update(name, value)

    if l_samp_stats_in_pdf_call and sclr_dim > 0:
        for sclr in range(sclr_dim):
            sclr_idx = sclr + 1
            stats = stats.update(f"sclr{sclr_idx}pthvp", sclrpthvp[:, :, sclr])
            stats = stats.update(f"sclr{sclr_idx}prcp", sclrprcp[:, :, sclr])
            stats = stats.update(f"wp2sclr{sclr_idx}p", wp2sclrp[:, :, sclr])
            stats = stats.update(f"wpsclr{sclr_idx}p2", wpsclrp2[:, :, sclr])
            stats = stats.update(f"wpsclr{sclr_idx}prtp", wpsclrprtp[:, :, sclr])
            stats = stats.update(f"wpsclr{sclr_idx}pthlp", wpsclrpthlp[:, :, sclr])

    return (
        rtm, pdf_implicit_coefs_terms, pdf_params, pdf_params_zm, err_info,
        rcm, cloud_frac, ice_supersat_frac, wprcp, sigma_sqd_w, wpthvp, wp2thvp,
        wp2up, rtpthvp, thlpthvp, rc_coef, rcm_in_layer, cloud_cover,
        rcp2_zt, thlprcp, rc_coef_zm, sclrpthvp, wpup2, wpvp2, wp2up2,
        wp2vp2, wp4, wp2rtp, wprtp2, wp2thlp, wpthlp2, wprtpthlp, wp2rcp,
        rtprcp, rcp2, uprcp, vprcp, w_up_in_cloud, w_down_in_cloud,
        cloudy_updraft_frac, cloudy_downdraft_frac, skw_velocity,
        cloud_frac_zm, ice_supersat_frac_zm, rtm_zm, thlm_zm, rcm_zm,
        rcm_supersat_adj, wp2sclrp, wpsclrp2, sclrprcp,
        wpsclrprtp, wpsclrpthlp, stats,
    )


def pdf_closure_driver_zm(
    gr, nzm, nzt, ngrdcol, hydromet_dim,
    sclr_dim, sclr_tol, p_sfc,
    mixt_frac_max_mag, clubb_params,
    iiPDF_type, saturation_formula,
    l_mix_rat_hm, p_in_Pa, thv_ds_zm,
    wm_zm, wp2, wp3_zm,
    Skw_zm, Skthl_zm, Skrt_zm,
    rtm, rtp2, wprtp,
    thlm, thlp2, wpthlp,
    um, up2, upwp, up3,
    vm, vp2, vpwp, vp3,
    rtpthlp, sclrm, wpsclrp, sclrp2,
    sclrprtp, sclrpthlp, sclrp3,
    wphydrometp,
    wp2hmp, rtphmp_zt, thlphmp_zt,
    stats, sigma_sqd_w, pdf_params_zm,
    err_info,
):
    """Port of pdf_closure_module.F90:pdf_closure_driver_zm."""
    p_in_Pa_zm = zt2zm(nzm, nzt, ngrdcol, gr, p_in_Pa)
    p_in_Pa_zm = p_in_Pa_zm.at[:, gr.k_lb_zm].set(p_sfc)
    p_in_Pa_zm = p_in_Pa_zm.at[:, gr.k_ub_zm].set(
        jnp.maximum(p_in_Pa_zm[:, gr.k_ub_zm], 0.5 * p_in_Pa[:, gr.k_ub_zt])
    )
    exner_zm = (p_in_Pa_zm / p0) ** kappa

    rtm_zm = zt2zm(nzm, nzt, ngrdcol, gr, rtm, rt_tol)
    thlm_zm = zt2zm(nzm, nzt, ngrdcol, gr, thlm, thl_tol)

    if hydromet_dim > 0:
        rtphmp = jnp.moveaxis(
            jax.vmap(lambda field: zt2zm(nzm, nzt, ngrdcol, gr, field), in_axes=2)(
                rtphmp_zt,
            ),
            0,
            -1,
        )
        thlphmp = jnp.moveaxis(
            jax.vmap(lambda field: zt2zm(nzm, nzt, ngrdcol, gr, field), in_axes=2)(
                thlphmp_zt,
            ),
            0,
            -1,
        )
        wp2hmp_zm = jnp.moveaxis(
            jax.vmap(lambda field: zt2zm(nzm, nzt, ngrdcol, gr, field), in_axes=2)(
                wp2hmp,
            ),
            0,
            -1,
        )
    else:
        rtphmp = jnp.zeros((ngrdcol, nzm, 0), dtype=jnp.float64)
        thlphmp = jnp.zeros((ngrdcol, nzm, 0), dtype=jnp.float64)
        wp2hmp_zm = jnp.zeros((ngrdcol, nzm, 0), dtype=jnp.float64)

    um_zm = zt2zm(nzm, nzt, ngrdcol, gr, um)
    vm_zm = zt2zm(nzm, nzt, ngrdcol, gr, vm)

    up3_zm = zt2zm(nzm, nzt, ngrdcol, gr, up3)
    vp3_zm = zt2zm(nzm, nzt, ngrdcol, gr, vp3)

    Sku_zm = Skx_func(nzm, ngrdcol, up2, up3_zm, w_tol, clubb_params)
    Skv_zm = Skx_func(nzm, ngrdcol, vp2, vp3_zm, w_tol, clubb_params)

    if sclr_dim > 0:
        sclrm_zm = jnp.moveaxis(
            jax.vmap(
                lambda field, tol: zt2zm(nzm, nzt, ngrdcol, gr, field, tol),
                in_axes=(2, 0),
            )(sclrm, sclr_tol),
            0,
            -1,
        )
        sclrp3_zm = jnp.moveaxis(
            jax.vmap(lambda field: zt2zm(nzm, nzt, ngrdcol, gr, field), in_axes=2)(
                sclrp3,
            ),
            0,
            -1,
        )
        Sksclr_zm = jnp.moveaxis(
            jax.vmap(
                lambda sclrp2_s, sclrp3_s, tol: Skx_func(
                    nzm, ngrdcol, sclrp2_s, sclrp3_s, tol, clubb_params,
                ),
                in_axes=(2, 2, 0),
            )(sclrp2, sclrp3_zm, sclr_tol),
            0,
            -1,
        )
    else:
        sclrm_zm = jnp.zeros((ngrdcol, nzm, 0), dtype=jnp.float64)
        sclrp3_zm = jnp.zeros((ngrdcol, nzm, 0), dtype=jnp.float64)
        Sksclr_zm = jnp.zeros((ngrdcol, nzm, 0), dtype=jnp.float64)

    if iiPDF_type in (iiPDF_new, iiPDF_new_hybrid):
        pdf_implicit_coefs_terms_zm = init_pdf_implicit_coefs_terms_api(
            nzm, ngrdcol, sclr_dim,
        )
    else:
        pdf_implicit_coefs_terms_zm = init_pdf_implicit_coefs_terms_api(
            nzm, ngrdcol, sclr_dim,
        )

    (
        pdf_params_zm, pdf_implicit_coefs_terms_zm, err_info,
        wpup2_zm, wpvp2_zm,
        wp2up2, wp2vp2, wp4,
        wprtp2_zm, wp2rtp_zm,
        wpthlp2_zm, wp2thlp_zm, wprtpthlp_zm,
        cloud_frac_zm, ice_supersat_frac_zm,
        rcm_zm, wpthvp, wp2thvp_zm, wp2up_zm,
        rtpthvp, thlpthvp, wprcp, wp2rcp_zm,
        rtprcp, thlprcp, rcp2, uprcp, vprcp,
        w_up_in_cloud_zm, w_down_in_cloud_zm,
        cloudy_updraft_frac_zm, cloudy_downdraft_frac_zm,
        wpsclrprtp_zm, wpsclrp2_zm, sclrpthvp,
        wpsclrpthlp_zm, sclrprcp, wp2sclrp_zm,
        rc_coef_zm,
    ) = pdf_closure(
        nzm, ngrdcol, sclr_dim, sclr_tol, gr,
        hydromet_dim, p_in_Pa_zm, exner_zm, thv_ds_zm,
        wm_zm, wp2, wp3_zm,
        Skw_zm, Skthl_zm, Skrt_zm, Sku_zm, Skv_zm,
        rtm_zm, rtp2, wprtp,
        thlm_zm, thlp2, wpthlp,
        um_zm, up2, upwp,
        vm_zm, vp2, vpwp,
        rtpthlp,
        sclrm_zm, wpsclrp, sclrp2,
        sclrprtp, sclrpthlp, Sksclr_zm,
        wphydrometp, wp2hmp_zm,
        rtphmp, thlphmp,
        clubb_params, mixt_frac_max_mag,
        saturation_formula,
        stats,
        iiPDF_type,
        l_mix_rat_hm,
        sigma_sqd_w,
        pdf_params_zm, pdf_implicit_coefs_terms_zm,
        err_info,
    )
    del pdf_implicit_coefs_terms_zm, wpup2_zm, wpvp2_zm
    del w_up_in_cloud_zm, w_down_in_cloud_zm
    del cloudy_updraft_frac_zm, cloudy_downdraft_frac_zm

    return (
        pdf_params_zm, err_info,
        rtm_zm, thlm_zm,
        wp2up2, wp2vp2, wp4,
        wprtp2_zm, wp2rtp_zm,
        wpthlp2_zm, wp2thlp_zm, wprtpthlp_zm,
        cloud_frac_zm, ice_supersat_frac_zm,
        rcm_zm, wpthvp, wp2thvp_zm, wp2up_zm,
        rtpthvp, thlpthvp, wprcp, wp2rcp_zm,
        rtprcp, thlprcp, rcp2, uprcp, vprcp,
        wpsclrprtp_zm, wpsclrp2_zm, sclrpthvp,
        wpsclrpthlp_zm, sclrprcp, wp2sclrp_zm,
        rc_coef_zm,
    )


def trapezoidal_rule_zt(
    nzm, nzt, ngrdcol, sclr_dim, gr, l_call_pdf_closure_twice, stats,
    wprtp2, wpthlp2, wprtpthlp, cloud_frac, ice_supersat_frac,
    rcm, wp2thvp, wp2up, wpsclrprtp, wpsclrp2, wpsclrpthlp,
    wprtp2_zm, wpthlp2_zm, wprtpthlp_zm, cloud_frac_zm,
    ice_supersat_frac_zm, rcm_zm, wp2thvp_zm, wp2up_zm,
    wpsclrprtp_zm, wpsclrp2_zm, wpsclrpthlp_zm,
):
    """Description:
    This subroutine takes the output variables on the thermo.
    grid and either: interpolates them to the momentum grid, or uses the
    values output from the second call to pdf_closure on momentum levels if
    l_call_pdf_closure_twice is true.  It then calls the function
    trapezoid_zt to recompute the variables on the thermo. grid.

    ldgrant June 2009

    Note:
    The argument variables in the last 5 lines of the subroutine
    (wprtp2_zm through pdf_params_zm) are declared intent(inout) because
    if l_call_pdf_closure_twice is true, these variables will already have
    values from pdf_closure on momentum levels and will not be altered in
    this subroutine.  However, if l_call_pdf_closure_twice is false, these
    variables will not have values yet and will be interpolated to
    momentum levels in this subroutine.
    References:
    None
    """
    if not l_call_pdf_closure_twice:
        # If l_call_pdf_closure_twice is true, the _zm variables already have
        # values from the second call to pdf_closure in advance_clubb_core.
        # If it is false, the variables are interpolated to the _zm levels.

        # Interpolate thermodynamic variables to the momentum grid.
        wprtp2_zm = zt2zm(nzm, nzt, ngrdcol, gr, wprtp2)
        wpthlp2_zm = zt2zm(nzm, nzt, ngrdcol, gr, wpthlp2)
        wprtpthlp_zm = zt2zm(nzm, nzt, ngrdcol, gr, wprtpthlp)
        cloud_frac_zm = zt2zm(nzm, nzt, ngrdcol, gr, cloud_frac)
        ice_supersat_frac_zm = zt2zm(nzm, nzt, ngrdcol, gr, ice_supersat_frac)
        rcm_zm = zt2zm(nzm, nzt, ngrdcol, gr, rcm)
        wp2thvp_zm = zt2zm(nzm, nzt, ngrdcol, gr, wp2thvp)
        wp2up_zm = zt2zm(nzm, nzt, ngrdcol, gr, wp2up)

        # Since top momentum level is higher than top thermo. level,
        # set variables at top momentum level to 0.
        wprtp2_zm = wprtp2_zm.at[:, gr.k_ub_zm].set(0.0)
        wpthlp2_zm = wpthlp2_zm.at[:, gr.k_ub_zm].set(0.0)
        wprtpthlp_zm = wprtpthlp_zm.at[:, gr.k_ub_zm].set(0.0)
        cloud_frac_zm = cloud_frac_zm.at[:, gr.k_ub_zm].set(0.0)
        ice_supersat_frac_zm = ice_supersat_frac_zm.at[:, gr.k_ub_zm].set(0.0)
        rcm_zm = rcm_zm.at[:, gr.k_ub_zm].set(0.0)
        wp2thvp_zm = wp2thvp_zm.at[:, gr.k_ub_zm].set(0.0)
        wp2up_zm = wp2up_zm.at[:, gr.k_ub_zm].set(0.0)

        if sclr_dim > 0:
            wpsclrprtp_zm = jnp.moveaxis(
                jax.vmap(
                    lambda field: zt2zm(nzm, nzt, ngrdcol, gr, field)
                    .at[:, gr.k_ub_zm].set(0.0),
                    in_axes=2,
                )(wpsclrprtp),
                0,
                -1,
            )
            wpsclrp2_zm = jnp.moveaxis(
                jax.vmap(
                    lambda field: zt2zm(nzm, nzt, ngrdcol, gr, field)
                    .at[:, gr.k_ub_zm].set(0.0),
                    in_axes=2,
                )(wpsclrp2),
                0,
                -1,
            )
            wpsclrpthlp_zm = jnp.moveaxis(
                jax.vmap(
                    lambda field: zt2zm(nzm, nzt, ngrdcol, gr, field)
                    .at[:, gr.k_ub_zm].set(0.0),
                    in_axes=2,
                )(wpsclrpthlp),
                0,
                -1,
            )

    if stats.names:
        # Use the trapezoidal rule to recompute the variables on the stats_zt level
        if stats.var_on_stats_list("wprtp2"):
            wprtp2 = calc_trapezoid_zt(nzm, nzt, ngrdcol, gr, wprtp2_zm, wprtp2)
        if stats.var_on_stats_list("wpthlp2"):
            wpthlp2 = calc_trapezoid_zt(nzm, nzt, ngrdcol, gr, wpthlp2_zm, wpthlp2)
        if stats.var_on_stats_list("wprtpthlp"):
            wprtpthlp = calc_trapezoid_zt(
                nzm, nzt, ngrdcol, gr, wprtpthlp_zm, wprtpthlp,
            )
        if sclr_dim > 0:
            wpsclrprtp = jnp.moveaxis(
                jax.vmap(
                    lambda variable_zm, variable_zt: calc_trapezoid_zt(
                        nzm, nzt, ngrdcol, gr, variable_zm, variable_zt,
                    ),
                    in_axes=(2, 2),
                )(wpsclrprtp_zm, wpsclrprtp),
                0,
                -1,
            )
            wpsclrpthlp = jnp.moveaxis(
                jax.vmap(
                    lambda variable_zm, variable_zt: calc_trapezoid_zt(
                        nzm, nzt, ngrdcol, gr, variable_zm, variable_zt,
                    ),
                    in_axes=(2, 2),
                )(wpsclrpthlp_zm, wpsclrpthlp),
                0,
                -1,
            )
            wpsclrp2 = jnp.moveaxis(
                jax.vmap(
                    lambda variable_zm, variable_zt: calc_trapezoid_zt(
                        nzm, nzt, ngrdcol, gr, variable_zm, variable_zt,
                    ),
                    in_axes=(2, 2),
                )(wpsclrp2_zm, wpsclrp2),
                0,
                -1,
            )

    cloud_frac = calc_trapezoid_zt(nzm, nzt, ngrdcol, gr, cloud_frac_zm, cloud_frac)
    ice_supersat_frac = calc_trapezoid_zt(
        nzm, nzt, ngrdcol, gr, ice_supersat_frac_zm, ice_supersat_frac,
    )
    rcm = calc_trapezoid_zt(nzm, nzt, ngrdcol, gr, rcm_zm, rcm)
    wp2thvp = calc_trapezoid_zt(nzm, nzt, ngrdcol, gr, wp2thvp_zm, wp2thvp)
    wp2up = calc_trapezoid_zt(nzm, nzt, ngrdcol, gr, wp2up_zm, wp2up)

    # End of trapezoidal rule
    return (
        wprtp2, wpthlp2, wprtpthlp, cloud_frac, ice_supersat_frac,
        rcm, wp2thvp, wp2up, wpsclrprtp, wpsclrp2, wpsclrpthlp,
        wprtp2_zm, wpthlp2_zm, wprtpthlp_zm, cloud_frac_zm,
        ice_supersat_frac_zm, rcm_zm, wp2thvp_zm, wp2up_zm,
        wpsclrprtp_zm, wpsclrp2_zm, wpsclrpthlp_zm,
    )


def trapezoidal_rule_zm(
    nzm, nzt, ngrdcol, gr, wpthvp_zt, thlpthvp_zt, rtpthvp_zt,
    wpthvp, thlpthvp, rtpthvp,
):
    """Description:
    This subroutine recomputes three variables on the
    momentum grid from pdf_closure -- wpthvp, thlpthvp, and
    rtpthvp -- by calling the function trapezoid_zm.  Only these three
    variables are used in this subroutine because they are the only
    pdf_closure momentum variables used elsewhere in CLUBB.

    The _zt variables are output from the first call to pdf_closure.
    The _zm variables are output from the second call to pdf_closure
    on the momentum levels.
    This is done before the call to this subroutine.

    ldgrant Feb. 2010

    References:
    None
    """
    # Use the trapezoidal rule to recompute the variables on the zm level
    wpthvp = calc_trapezoid_zm(nzm, nzt, ngrdcol, gr, wpthvp_zt, wpthvp)
    thlpthvp = calc_trapezoid_zm(nzm, nzt, ngrdcol, gr, thlpthvp_zt, thlpthvp)
    rtpthvp = calc_trapezoid_zm(nzm, nzt, ngrdcol, gr, rtpthvp_zt, rtpthvp)
    return wpthvp, thlpthvp, rtpthvp


def calc_trapezoid_zt(nzm, nzt, ngrdcol, gr, variable_zm, variable_zt):
    """Description:
    Function which uses the trapezoidal rule from calculus
    to recompute the values for the variables on the thermo. grid which
    are output from the first call to pdf_closure in module clubb_core.

    ldgrant June 2009
    """
    del nzm, nzt, ngrdcol

    def body(level, variable_zt):
        k = gr.k_lb_zt + gr.grid_dir_indx * level
        if gr.grid_dir_indx > 0:
            # Ascending grid
            k_zmp1 = k + 1
            k_zm = k
        else:
            # Descending grid
            k_zmp1 = k
            k_zm = k + 1

        # Trapezoidal rule from calculus
        variable_zt = variable_zt.at[:, k].set(
            0.5
            * (variable_zm[:, k_zmp1] + variable_zt[:, k])
            * (gr.zm[:, k_zmp1] - gr.zt[:, k])
            * gr.grid_dir * gr.invrs_dzt[:, k]
            + 0.5
            * (variable_zt[:, k] + variable_zm[:, k_zm])
            * (gr.zt[:, k] - gr.zm[:, k_zm])
            * gr.grid_dir * gr.invrs_dzt[:, k]
        )
        return variable_zt

    nlev = abs(gr.k_ub_zt - gr.k_lb_zt) + 1
    variable_zt = jax.lax.fori_loop(0, nlev, body, variable_zt)
    return variable_zt


def calc_trapezoid_zm(nzm, nzt, ngrdcol, gr, variable_zt, variable_zm):
    """Description:
    Function which uses the trapezoidal rule from calculus
    to recompute the values for the important variables on the momentum
    grid which are output from pdf_closure in module clubb_core.
    These momentum variables only include wpthvp, thlpthvp, and rtpthvp.

    ldgrant Feb. 2010
    """
    del nzm, nzt, ngrdcol
    start = gr.k_lb_zm + gr.grid_dir_indx
    stop = gr.k_ub_zm - gr.grid_dir_indx

    # Boundary conditions: trapezoidal rule not valid at top zm level, nzmax.
    # Trapezoidal rule also not used at zm level 1.
    def body(level, variable_zm):
        k = start + gr.grid_dir_indx * level
        if gr.grid_dir_indx > 0:
            # Ascending grid
            k_zt = k
            k_ztm1 = k - 1
        else:
            # Descending grid
            k_zt = k - 1
            k_ztm1 = k

        # Trapezoidal rule from calculus
        variable_zm = variable_zm.at[:, k].set(
            0.5
            * (variable_zt[:, k_zt] + variable_zm[:, k])
            * (gr.zt[:, k_zt] - gr.zm[:, k])
            * gr.grid_dir * gr.invrs_dzm[:, k]
            + 0.5
            * (variable_zm[:, k] + variable_zt[:, k_ztm1])
            * (gr.zm[:, k] - gr.zt[:, k_ztm1])
            * gr.grid_dir * gr.invrs_dzm[:, k]
        )
        return variable_zm

    nlev = max((stop - start) * gr.grid_dir_indx + 1, 0)
    variable_zm = jax.lax.fori_loop(0, nlev, body, variable_zm)
    return variable_zm


def compute_cloud_cover(gr, nzt, ngrdcol, pdf_params, cloud_frac, rcm):
    """Description:
    Subroutine to compute cloud cover (the amount of sky
    covered by cloud) and rcm in layer (liquid water mixing ratio in
    the portion of the grid box filled by cloud).

    References:
    Definition of 's' comes from:
    ``The Gaussian Cloud Model Relations'' G. L. Mellor (1977)
    JAS, Vol. 34, pp. 356--358.

    Notes:
    Added July 2009

    JAX port note: the Fortran routine has an unreachable fatal-error branch
    for invalid condition coverage. This JAX form expresses the same reachable
    cases with masks and leaves unchanged values outside them.
    """
    chi_mean = (
        pdf_params.mixt_frac * pdf_params.chi_1
        + (1.0 - pdf_params.mixt_frac) * pdf_params.chi_2
    )
    cloud_cover = cloud_frac
    rcm_in_layer = rcm

    nlev = abs(gr.k_ub_zt - gr.k_lb_zt)

    def body(level, state):
        cloud_cover, rcm_in_layer = state
        k = gr.k_lb_zt + gr.grid_dir_indx * level
        if gr.grid_dir_indx > 0:
            km1 = jnp.maximum(k - 1, 0)
            kp1 = jnp.minimum(k + 1, nzt - 1)
            k_zmp1 = k + 1
            k_zm = k
        else:
            km1 = jnp.minimum(k + 1, nzt - 1)
            kp1 = jnp.maximum(k - 1, 0)
            k_zmp1 = k
            k_zm = k + 1

        no_cloud = rcm[:, k] < rc_tol
        filled = (rcm[:, kp1] >= rc_tol) & (rcm[:, km1] >= rc_tol)

        # First let the cloud fill the entire grid box, then overwrite
        # vert_cloud_frac_upper(k) and/or vert_cloud_frac_lower(k)
        # for a cloud top, cloud base, or one-point cloud.
        upper = jnp.full((ngrdcol,), 0.5, dtype=jnp.float64)
        lower = jnp.full((ngrdcol,), 0.5, dtype=jnp.float64)

        # Cloud top
        top = rcm[:, kp1] < rc_tol
        top_scale = (
            (0.5 / (gr.grid_dir * gr.invrs_dzm[:, k_zmp1]))
            / (gr.zm[:, k_zmp1] - gr.zt[:, k])
        )
        top_frac = top_scale * (
            rcm[:, k] / (rcm[:, k] + jnp.abs(chi_mean[:, kp1]))
        )
        top_frac = jnp.minimum(0.5, top_frac)
        # Make the transition in cloudiness more gradual than using
        # the above min statement alone.
        top_frac = top_frac + (rcm[:, kp1] / rc_tol) * (0.5 - top_frac)
        upper = jnp.where(top, top_frac, upper)

        # Cloud base
        base = rcm[:, km1] < rc_tol
        base_scale = (
            (0.5 / (gr.grid_dir * gr.invrs_dzm[:, k_zm]))
            / (gr.zt[:, k] - gr.zm[:, k_zm])
        )
        base_frac = base_scale * (
            rcm[:, k] / (rcm[:, k] + jnp.abs(chi_mean[:, km1]))
        )
        base_frac = jnp.minimum(0.5, base_frac)
        # Make the transition in cloudiness more gradual than using
        # the above min statement alone.
        base_frac = base_frac + (rcm[:, km1] / rc_tol) * (0.5 - base_frac)
        lower = jnp.where(base, base_frac, lower)

        vert_cloud_frac = upper + lower
        vert_cloud_frac = jnp.maximum(
            cloud_frac[:, k],
            jnp.minimum(1.0, vert_cloud_frac),
        )
        maybe_cover = cloud_frac[:, k] / vert_cloud_frac
        maybe_rcm = rcm[:, k] / vert_cloud_frac

        use_layer = (~no_cloud) & (~filled)
        cloud_cover = cloud_cover.at[:, k].set(
            jnp.where(use_layer, maybe_cover, cloud_cover[:, k])
        )
        rcm_in_layer = rcm_in_layer.at[:, k].set(
            jnp.where(use_layer, maybe_rcm, rcm_in_layer[:, k])
        )
        return cloud_cover, rcm_in_layer

    cloud_cover, rcm_in_layer = jax.lax.fori_loop(
        0, nlev, body, (cloud_cover, rcm_in_layer),
    )
    cloud_cover = cloud_cover.at[:, gr.k_ub_zt].set(cloud_frac[:, gr.k_ub_zt])
    rcm_in_layer = rcm_in_layer.at[:, gr.k_ub_zt].set(rcm[:, gr.k_ub_zt])
    return cloud_cover, rcm_in_layer


__all__ = [
    "calc_wp2xp_pdf", "calc_wpxp2_pdf", "calc_wp2xp2_pdf",
    "calc_wp4_pdf", "calc_wpxpyp_pdf",
    "calc_w_up_in_cloud",
    "transform_pdf_chi_eta_component",
    "calc_liquid_cloud_frac_component",
    "calc_ice_cloud_frac_component",
    "calc_xprcp_component",
    "pdf_closure",
    "trapezoidal_rule_zt",
    "trapezoidal_rule_zm",
    "calc_trapezoid_zt",
    "calc_trapezoid_zm",
    "compute_cloud_cover",
    "pdf_closure_driver_zm",
    "pdf_closure_driver",
]
