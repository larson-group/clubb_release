"""JAX-side entry point for ``src/CLUBB_core/advance_xm_wpxp_module.F90``.

Description:
  Contains the CLUBB advance_xm_wpxp_module scheme. Advance the mean and flux
  terms by one timestep.

References:
  https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:wpxp_eqns

  Eqn. 16 & 17 on p. 3546 of
  ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    Method and Model Description'' Golaz, et al. (2002)
    JAS, Vol. 59, pp. 3540--3551.

See Also
  ``Equations for CLUBB'' Section 5:
    /Implicit solutions for the means and fluxes/

Porting deviations:
- This JAX entry point preserves the Fortran argument and return contract while
  using JAX-side stats, limiter, positive-definite, skewness, and matrix-solver
  helpers inside the jitted region.
- Sponge damping blocks are unsupported here because clubb_case_initalization
  rejects sponge-enabled Python/JAX driver cases before this routine is called.
- error_prints_xm_wpxp is intentionally reduced until full JAX diagnostic state
  is available.
- _solve_one_xm_wpxp preserves the repeated solve structure of the Fortran
  multiple-LHS path while avoiding a larger target-only abstraction.
"""

from functools import partial

import jax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.Skx_module import Skx_func
from clubb_jax.src.CLUBB_core.advance_helper_module import calc_wp3_on_wp2
from clubb_jax.src.CLUBB_core.clip_explicit import (
    clip_covar,
    clip_upwp,
    clip_vpwp,
    clip_wprtp,
    clip_wpsclrp,
    clip_wpthlp,
    upwp_cl_max,
    vpwp_cl_max,
    wprtp_cl_max,
    wpthlp_cl_max,
)
from clubb_jax.src.CLUBB_core.diffusion import diffusion_zt_lhs, diffusion_zm_lhs
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from clubb_jax.src.CLUBB_core.fill_holes import fill_holes_vertical
from clubb_jax.src.CLUBB_core.mean_adv import term_ma_zt_lhs, term_ma_zm_lhs
from clubb_jax.src.CLUBB_core.matrix_solver_wrapper import band_solve
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.mono_flux_limiter import (
    calc_turb_adv_range,
    monotonic_turbulent_flux_limit,
)
from clubb_jax.src.CLUBB_core.pos_definite_module import pos_definite_adj
from clubb_jax.src.CLUBB_core.turbulent_adv_pdf import (
    xpyp_term_ta_pdf_lhs,
    xpyp_term_ta_pdf_lhs_godunov,
    xpyp_term_ta_pdf_rhs,
)
from clubb_jax.src.CLUBB_core.constants_clubb import (
    eps,
    ep1,
    gamma_over_implicit_ts,
    grav,
    one,
    one_half,
    rt_tol,
    rt_tol_mfl,
    thl_tol,
    thl_tol_mfl,
    w_tol,
    w_tol_sqd,
    zero,
    zero_threshold,
)
from clubb_jax.src.CLUBB_core.parameter_indices import (
    ic_K6,
    iC6rt,
    iC6rt_Lscale0,
    iC6rtb,
    iC6rtc,
    iC6thl,
    iC6thl_Lscale0,
    iC6thlb,
    iC6thlc,
    iC7,
    iC7_Lscale0,
    iC7b,
    iC7c,
    iC_uu_shr,
    ialtitude_threshold,
    iwpxp_L_thresh,
)
from clubb_jax.src.CLUBB_core.model_flags import (
    iiPDF_ADG1,
    iiPDF_new,
    iiPDF_new_hybrid,
    l_clip_turb_adv,
    l_explicit_turbulent_adv_wpxp,
    l_force_descending_solves,
    l_pos_def,
    penta_bicgstab,
)
from clubb_jax.src.CLUBB_core.grid_class import (
    ddzt,
    zm2zt,
    zm2zt2zm,
    zt2zm,
    zt2zm2zt,
)
from clubb_jax.src.CLUBB_core import (
    ErrInfo,
    Grid,
    NuVertResDep,
    implicit_coefs_terms,
)

# Parameter Constants
nsub = 2  # Number of subdiagonals in the LHS matrix
nsup = 2  # Number of superdiagonals in the LHS matrix
xm_wpxp_thlm = 1  # Named constant for thlm and wpthlp solving
xm_wpxp_rtm = 2  # Named constant for rtm and wprtp solving
xm_wpxp_scalar = 3  # Named constant for sclrm and wpsclrp solving
xm_wpxp_um = 4  # Named constant for optional um and upwp solving
xm_wpxp_vm = 5  # Named constant for optional vm and vpwp solving

ndiags2 = 2
ndiags3 = 3
ndiags5 = 5


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "nzt",
        "ngrdcol",
        "sclr_dim",
        "iipdf_type",
        "penta_solve_method",
        "tridiag_solve_method",
        "fill_holes_type",
        "l_predict_upwp_vpwp",
        "l_ho_nontrad_coriolis",
        "l_ho_trad_coriolis",
        "l_diffuse_rtm_and_thlm",
        "l_stability_correct_kh_n2_zm",
        "l_godunov_upwind_wpxp_ta",
        "l_upwind_xm_ma",
        "l_uv_nudge",
        "l_diag_lscale_from_tau",
        "l_use_c7_richardson",
        "l_lmm_stepping",
        "l_enable_relaxed_clipping",
        "l_linearize_pbl_winds",
        "l_mono_flux_lim_thlm",
        "l_mono_flux_lim_rtm",
        "l_mono_flux_lim_um",
        "l_mono_flux_lim_vm",
        "l_mono_flux_lim_spikefix",
        "wprtp_cl_num",
        "wpthlp_cl_num",
        "upwp_cl_num",
        "vpwp_cl_num",
        "l_implemented",
    ),
)
def advance_xm_wpxp(
    nzm: int, nzt: int, ngrdcol: int, sclr_dim: int, sclr_tol, gr: Grid, dt: float,
    sigma_sqd_w, wm_zm, wm_zt, wp2, lscale_zm,
    wp3, kh_zt, kh_zm,
    stability_correction,
    invrs_tau_c6_zm, tau_max_zm, wp2rtp, rtpthvp,
    rtm_forcing, wprtp_forcing, rtm_ref, wp2thlp,
    thlpthvp, thlm_forcing, wpthlp_forcing, thlm_ref,
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, thv_ds_zm, rtp2, thlp2,
    w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, mixt_frac_zm,
    l_implemented: bool, em, wp2sclrp, sclrpthvp, sclrm_forcing, sclrp2, cx_fnc_richardson,
    pdf_implicit_coefs_terms: implicit_coefs_terms,
    um_forcing, vm_forcing, ug, vg, wpthvp,
    fcor, fcor_y, um_ref, vm_ref, up2, vp2, uprcp, vprcp, rc_coef_zm,
    clubb_params, nu_vert_res_dep: NuVertResDep, ts_nudge: float,
    iipdf_type: int, penta_solve_method: int, tridiag_solve_method: int, fill_holes_type: int,
    l_predict_upwp_vpwp: bool, l_ho_nontrad_coriolis: bool, l_ho_trad_coriolis: bool,
    l_diffuse_rtm_and_thlm: bool, l_stability_correct_kh_n2_zm: bool,
    l_godunov_upwind_wpxp_ta: bool, l_upwind_xm_ma: bool, l_uv_nudge: bool,
    l_diag_lscale_from_tau: bool, l_use_c7_richardson: bool,
    l_lmm_stepping: bool, l_enable_relaxed_clipping: bool, l_linearize_pbl_winds: bool,
    l_mono_flux_lim_thlm: bool, l_mono_flux_lim_rtm: bool, l_mono_flux_lim_um: bool,
    l_mono_flux_lim_vm: bool, l_mono_flux_lim_spikefix: bool,
    wprtp_cl_num: int, wpthlp_cl_num: int, upwp_cl_num: int, vpwp_cl_num: int,
    stats: JaxStats,
    rtm, wprtp, thlm, wpthlp, sclrm, wpsclrp, um, upwp, vm, vpwp,
    um_pert, vm_pert, upwp_pert, vpwp_pert,
    err_info: ErrInfo,
):
    """Advance the mean and flux terms by one timestep."""
    # The Python API owns scalar storage as NumPy arrays.  Keep the Fortran
    # inout-array behavior while normalizing those arrays at the JAX boundary
    # so scalar cases can use functional indexed updates in eager or jitted
    # execution alike.
    sclrm = jnp.asarray(sclrm)
    wpsclrp = jnp.asarray(wpsclrp)

    l_iter = True
    l_sample = stats.l_sample
    l_perturbed_wind = l_predict_upwp_vpwp and l_linearize_pbl_winds

    wp3_zm = zt2zm(nzm, nzt, ngrdcol, gr, wp3)
    Skw_zm = Skx_func(nzm, ngrdcol, wp2, wp3_zm, w_tol, clubb_params)
    wp3_on_wp2, wp3_on_wp2_zt = calc_wp3_on_wp2(
        nzm, nzt, ngrdcol, gr, wp2, wp3,
    )

    if clubb_at_least_debug_level(0):
        if l_mono_flux_lim_rtm and not l_mono_flux_lim_spikefix:
            err_info = err_info.set_fatal()
            return (
                wprtp_cl_num, wpthlp_cl_num, upwp_cl_num, vpwp_cl_num,
                rtm, wprtp, thlm, wpthlp, sclrm, wpsclrp, um, upwp, vm, vpwp,
                um_pert, vm_pert, upwp_pert, vpwp_pert, err_info, stats,
            )

    l_scalar_calc = sclr_dim > 0
    if iipdf_type == iiPDF_new and not l_explicit_turbulent_adv_wpxp:
        nrhs = 1
    else:
        nrhs = 2 + sclr_dim
        if l_predict_upwp_vpwp:
            nrhs += 2
            if l_perturbed_wind:
                nrhs += 2

    rtm_old = rtm
    thlm_old = thlm
    wprtp_old = wprtp
    wpthlp_old = wpthlp
    sclrm_old = sclrm
    wpsclrp_old = wpsclrp
    um_old = um
    vm_old = vm
    upwp_old = upwp
    vpwp_old = vpwp

    if not l_diag_lscale_from_tau:
        C6rt = clubb_params[:, iC6rt]
        C6rtb = clubb_params[:, iC6rtb]
        C6rtc = clubb_params[:, iC6rtc]
        C6rtc_safe = jnp.where(jnp.abs(C6rtc) > zero, C6rtc, one)
        C6rt_varying = jnp.abs(C6rt - C6rtb) > jnp.abs(C6rt + C6rtb) * eps / 2.0
        C6rt_computed = (
            C6rtb[:, None]
            + (C6rt - C6rtb)[:, None]
            * jnp.exp(-one_half * (jnp.asarray(Skw_zm) / C6rtc_safe[:, None]) ** 2)
        )
        C6rt_Skw_fnc = jnp.where(C6rt_varying[:, None], C6rt_computed, C6rtb[:, None])

        C6thl = clubb_params[:, iC6thl]
        C6thlb = clubb_params[:, iC6thlb]
        C6thlc = clubb_params[:, iC6thlc]
        C6thlc_safe = jnp.where(jnp.abs(C6thlc) > zero, C6thlc, one)
        C6thl_varying = jnp.abs(C6thl - C6thlb) > jnp.abs(C6thl + C6thlb) * eps / 2.0
        C6thl_computed = (
            C6thlb[:, None]
            + (C6thl - C6thlb)[:, None]
            * jnp.exp(-one_half * (jnp.asarray(Skw_zm) / C6thlc_safe[:, None]) ** 2)
        )
        C6thl_Skw_fnc = jnp.where(C6thl_varying[:, None], C6thl_computed, C6thlb[:, None])
        C6rt_Skw_fnc = damp_coefficient(
            nzm, ngrdcol, gr, clubb_params[:, iC6rt], C6rt_Skw_fnc,
            clubb_params[:, iC6rt_Lscale0],
            clubb_params[:, ialtitude_threshold],
            clubb_params[:, iwpxp_L_thresh], lscale_zm,
        )
        C6thl_Skw_fnc = damp_coefficient(
            nzm, ngrdcol, gr, clubb_params[:, iC6thl], C6thl_Skw_fnc,
            clubb_params[:, iC6thl_Lscale0],
            clubb_params[:, ialtitude_threshold],
            clubb_params[:, iwpxp_L_thresh], lscale_zm,
        )
    else:
        C6rt_Skw_fnc = jnp.broadcast_to(clubb_params[:, iC6rt, None], (ngrdcol, nzm)).copy()
        C6thl_Skw_fnc = jnp.broadcast_to(clubb_params[:, iC6thl, None], (ngrdcol, nzm)).copy()

    if l_use_c7_richardson:
        C7_Skw_fnc = jnp.asarray(cx_fnc_richardson, dtype=jnp.float64)
    else:
        C7 = clubb_params[:, iC7]
        C7b = clubb_params[:, iC7b]
        C7c = clubb_params[:, iC7c]
        C7c_safe = jnp.where(jnp.abs(C7c) > zero, C7c, one)
        C7_varying = jnp.abs(C7 - C7b) > jnp.abs(C7 + C7b) * eps / 2.0
        C7_computed = (
            C7b[:, None]
            + (C7 - C7b)[:, None]
            * jnp.exp(-one_half * (jnp.asarray(Skw_zm) / C7c_safe[:, None]) ** 2)
        )
        C7_Skw_fnc = jnp.where(C7_varying[:, None], C7_computed, C7b[:, None])
        C7_Skw_fnc = damp_coefficient(
            nzm, ngrdcol, gr, clubb_params[:, iC7], C7_Skw_fnc,
            clubb_params[:, iC7_Lscale0],
            clubb_params[:, ialtitude_threshold],
            clubb_params[:, iwpxp_L_thresh], lscale_zm,
        )

    if l_sample:
        stats = stats.update("C7_Skw_fnc", C7_Skw_fnc)
        stats = stats.update("C6rt_Skw_fnc", C6rt_Skw_fnc)
        stats = stats.update("C6thl_Skw_fnc", C6thl_Skw_fnc)

    if clubb_at_least_debug_level(0):
        c7_bad = jnp.any(
            (jnp.asarray(C7_Skw_fnc) > one) | (jnp.asarray(C7_Skw_fnc) < zero),
            axis=1,
        )
        err_info = err_info.set_fatal(mask=c7_bad)

    Kw6 = clubb_params[:, ic_K6, None] * kh_zt
    low_lev_effect, high_lev_effect, stats = calc_turb_adv_range(
        nzm, nzt, ngrdcol, gr, dt,
        w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, mixt_frac_zm,
        stats,
    )

    lhs_pr1_wprtp, lhs_pr1_wpthlp, lhs_pr1_wpsclrp = wpxp_term_pr1_lhs(
        nzm, ngrdcol, gr, C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc,
        invrs_tau_c6_zm, l_scalar_calc,
    )

    C6_term = C6rt_Skw_fnc * invrs_tau_c6_zm
    if l_sample:
        stats = stats.update("C6_term", C6_term)

    (
        lhs_ta_wprtp, lhs_ta_wpthlp, lhs_ta_wpup,
        lhs_ta_wpvp, lhs_ta_wpsclrp,
        rhs_ta_wprtp, rhs_ta_wpthlp, rhs_ta_wpup,
        rhs_ta_wpvp, rhs_ta_wpsclrp, stats,
    ) = calc_xm_wpxp_ta_terms(
        nzm, nzt, ngrdcol, sclr_dim, gr, wp2rtp,
        wp2thlp, wp2sclrp,
        rho_ds_zt, invrs_rho_ds_zm, rho_ds_zm,
        sigma_sqd_w, wp3_on_wp2_zt,
        pdf_implicit_coefs_terms,
        iipdf_type,
        l_explicit_turbulent_adv_wpxp, l_predict_upwp_vpwp,
        l_scalar_calc,
        l_godunov_upwind_wpxp_ta,
        stats,
    )

    (
        lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm,
        lhs_tp, lhs_ta_xm, lhs_ac_pr2,
    ) = calc_xm_wpxp_lhs_terms(
        nzm, nzt, ngrdcol, gr, wm_zm, wm_zt, wp2,
        kh_zm, stability_correction, Kw6, C7_Skw_fnc,
        invrs_rho_ds_zt, invrs_rho_ds_zm, rho_ds_zt,
        rho_ds_zm, l_implemented, nu_vert_res_dep,
        l_diffuse_rtm_and_thlm,
        l_stability_correct_kh_n2_zm,
        l_upwind_xm_ma,
    )

    if iipdf_type == iiPDF_new and not l_explicit_turbulent_adv_wpxp:
        (
            wprtp_cl_num, wpthlp_cl_num,
            rtm, wprtp, thlm, wpthlp, sclrm, wpsclrp, err_info, stats,
        ) = solve_xm_wpxp_with_multiple_lhs(
            nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, dt,
            l_iter, nrhs, wm_zt, wp2,
            rtpthvp, rtm_forcing, wprtp_forcing, thlpthvp,
            thlm_forcing, wpthlp_forcing, rho_ds_zm,
            rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,
            thv_ds_zm, rtp2, thlp2, l_implemented,
            sclrpthvp, sclrm_forcing, sclrp2,
            low_lev_effect, high_lev_effect, C7_Skw_fnc,
            lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm,
            lhs_ta_wprtp, lhs_ta_wpthlp, lhs_ta_wpsclrp,
            rhs_ta_wprtp, rhs_ta_wpthlp, rhs_ta_wpsclrp,
            lhs_tp, lhs_ta_xm, lhs_ac_pr2, lhs_pr1_wprtp,
            lhs_pr1_wpthlp, lhs_pr1_wpsclrp,
            penta_solve_method,
            tridiag_solve_method,
            fill_holes_type,
            l_diffuse_rtm_and_thlm,
            l_upwind_xm_ma,
            l_enable_relaxed_clipping,
            l_mono_flux_lim_thlm,
            l_mono_flux_lim_rtm,
            l_mono_flux_lim_um,
            l_mono_flux_lim_vm,
            l_mono_flux_lim_spikefix,
            int(wprtp_cl_num), int(wpthlp_cl_num),
            stats,
            rtm, wprtp, thlm, wpthlp, sclrm, wpsclrp, err_info,
        )
    else:
        (
            wprtp_cl_num, wpthlp_cl_num, upwp_cl_num, vpwp_cl_num,
            rtm, wprtp, thlm, wpthlp, sclrm, wpsclrp, um, upwp, vm, vpwp,
            um_pert, vm_pert, upwp_pert, vpwp_pert, err_info, stats,
        ) = solve_xm_wpxp_with_single_lhs(
            nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, dt,
            l_iter, nrhs,
            wm_zt, wp2, invrs_tau_c6_zm, tau_max_zm,
            rtpthvp, rtm_forcing, wprtp_forcing, thlpthvp,
            thlm_forcing, wpthlp_forcing, rho_ds_zm,
            rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,
            thv_ds_zm, rtp2, thlp2, l_implemented,
            sclrpthvp, sclrm_forcing, sclrp2, um_forcing,
            vm_forcing, ug, vg, uprcp, vprcp, rc_coef_zm, fcor,
            fcor_y, up2, vp2,
            low_lev_effect, high_lev_effect,
            C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc,
            lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm,
            lhs_ta_wprtp,
            rhs_ta_wprtp, rhs_ta_wpthlp, rhs_ta_wpup,
            rhs_ta_wpvp, rhs_ta_wpsclrp,
            lhs_tp, lhs_ta_xm, lhs_ac_pr2, lhs_pr1_wprtp,
            lhs_pr1_wpthlp, lhs_pr1_wpsclrp,
            clubb_params[:, iC_uu_shr],
            penta_solve_method,
            tridiag_solve_method,
            fill_holes_type,
            l_predict_upwp_vpwp,
            l_ho_nontrad_coriolis,
            l_ho_trad_coriolis,
            l_diffuse_rtm_and_thlm,
            l_upwind_xm_ma,
            l_enable_relaxed_clipping,
            l_perturbed_wind,
            l_mono_flux_lim_thlm,
            l_mono_flux_lim_rtm,
            l_mono_flux_lim_um,
            l_mono_flux_lim_vm,
            l_mono_flux_lim_spikefix,
            int(wprtp_cl_num), int(wpthlp_cl_num), int(upwp_cl_num), int(vpwp_cl_num),
            stats,
            rtm, wprtp, thlm, wpthlp,
            sclrm, wpsclrp, um, upwp, vm, vpwp,
            um_pert, vm_pert, upwp_pert, vpwp_pert, err_info,
        )

    if l_lmm_stepping:
        thlm = one_half * (thlm_old + thlm)
        rtm = one_half * (rtm_old + rtm)
        wpthlp = one_half * (wpthlp_old + wpthlp)
        wprtp = one_half * (wprtp_old + wprtp)

        for sclr in range(sclr_dim):
            sclrm = sclrm.at[:, :, sclr].set(
                one_half * (sclrm_old[:, :, sclr] + sclrm[:, :, sclr])
            )
            wpsclrp = wpsclrp.at[:, :, sclr].set(
                one_half * (wpsclrp_old[:, :, sclr] + wpsclrp[:, :, sclr])
            )

        if l_predict_upwp_vpwp:
            um = one_half * (um_old + um)
            vm = one_half * (vm_old + vm)
            upwp = one_half * (upwp_old + upwp)
            vpwp = one_half * (vpwp_old + vpwp)

    # Fortran applies rtm/thlm/uv sponge damping here. The Python/JAX driver
    # rejects sponge-enabled cases in clubb_case_initalization, so these blocks
    # remain unsupported until sponge-layer state is owned by the JAX path.

    if l_predict_upwp_vpwp:
        if l_uv_nudge:
            if l_sample:
                stats = stats.begin_budget("um_ndg", um / dt)
                stats = stats.begin_budget("vm_ndg", vm / dt)
            um = um - ((um - um_ref) * (dt / ts_nudge))
            vm = vm - ((vm - vm_ref) * (dt / ts_nudge))
            if l_sample:
                stats = stats.finalize_budget("um_ndg", um / dt)
                stats = stats.finalize_budget("vm_ndg", vm / dt)

        if l_sample:
            stats = stats.update("um_ref", um_ref)
            stats = stats.update("vm_ref", vm_ref)

    return (
        wprtp_cl_num, wpthlp_cl_num, upwp_cl_num, vpwp_cl_num,
        rtm, wprtp, thlm, wpthlp, sclrm, wpsclrp, um, upwp, vm, vpwp,
        um_pert, vm_pert, upwp_pert, vpwp_pert, err_info, stats,
    )


def xm_wpxp_lhs(
    nzm, nzt, ngrdcol, l_iter, dt,
    l_implemented, lhs_diff_zm, lhs_diff_zt,
    lhs_ma_zm, lhs_ma_zt, lhs_ta_wpxp, lhs_ta_xm,
    lhs_tp, lhs_pr1, lhs_ac_pr2,
    l_diffuse_rtm_and_thlm,
):
    """Compute LHS band diagonal matrix for xm and w'x'."""
    invrs_dt = one / dt
    lhs = jnp.zeros((nsup + nsub + 1, ngrdcol, 2 * nzm - 1), dtype=jnp.float64)

    # Lower (upper) boundary for w'x' on ascending (descending) grid.
    lhs = lhs.at[2, :, 0].set(one)

    # Combine xm and w'x' terms into LHS
    k_xm = slice(1, 2 * nzt + 1, 2)
    lhs = lhs.at[1, :, k_xm].set(jnp.asarray(lhs_ta_xm[0, :, :]))
    lhs = lhs.at[2, :, k_xm].set(invrs_dt)
    lhs = lhs.at[3, :, k_xm].set(jnp.asarray(lhs_ta_xm[1, :, :]))

    k_wpxp = slice(2, 2 * nzm - 2, 2)
    lhs = lhs.at[0, :, k_wpxp].set(
        jnp.asarray(lhs_ma_zm[0, :, 1:nzm - 1])
        + jnp.asarray(lhs_diff_zm[0, :, 1:nzm - 1])
        + gamma_over_implicit_ts * jnp.asarray(lhs_ta_wpxp[0, :, 1:nzm - 1])
    )
    lhs = lhs.at[1, :, k_wpxp].set(jnp.asarray(lhs_tp[0, :, 1:nzm - 1]))
    lhs = lhs.at[2, :, k_wpxp].set(
        jnp.asarray(lhs_ma_zm[1, :, 1:nzm - 1])
        + jnp.asarray(lhs_diff_zm[1, :, 1:nzm - 1])
        + jnp.asarray(lhs_ac_pr2[:, 1:nzm - 1])
        + gamma_over_implicit_ts
        * (jnp.asarray(lhs_ta_wpxp[1, :, 1:nzm - 1]) + jnp.asarray(lhs_pr1[:, 1:nzm - 1]))
    )
    lhs = lhs.at[3, :, k_wpxp].set(jnp.asarray(lhs_tp[1, :, 1:nzm - 1]))
    lhs = lhs.at[4, :, k_wpxp].set(
        jnp.asarray(lhs_ma_zm[2, :, 1:nzm - 1])
        + jnp.asarray(lhs_diff_zm[2, :, 1:nzm - 1])
        + gamma_over_implicit_ts * jnp.asarray(lhs_ta_wpxp[2, :, 1:nzm - 1])
    )

    # Upper (lower) boundary for w'x' on ascending (descending) grid.
    lhs = lhs.at[2, :, 2 * nzm - 2].set(one)

    # LHS time tendency
    if l_iter:
        lhs = lhs.at[2, :, k_wpxp].add(invrs_dt)

    # Calculate diffusion terms for all thermodynamic grid level
    if l_diffuse_rtm_and_thlm:
        lhs = lhs.at[0, :, k_xm].add(jnp.asarray(lhs_diff_zt[0, :, :]))
        lhs = lhs.at[2, :, k_xm].add(jnp.asarray(lhs_diff_zt[1, :, :]))
        lhs = lhs.at[4, :, k_xm].add(jnp.asarray(lhs_diff_zt[2, :, :]))

    # Calculate mean advection terms for all momentum grid level
    if not l_implemented:
        lhs = lhs.at[0, :, k_xm].add(jnp.asarray(lhs_ma_zt[0, :, :]))
        lhs = lhs.at[2, :, k_xm].add(jnp.asarray(lhs_ma_zt[1, :, :]))
        lhs = lhs.at[4, :, k_xm].add(jnp.asarray(lhs_ma_zt[2, :, :]))

    return lhs


def calc_xm_wpxp_lhs_terms(
    nzm, nzt, ngrdcol, gr, wm_zm, wm_zt, wp2,
    Kh_zm, stability_correction, Kw6, C7_Skw_fnc,
    invrs_rho_ds_zt, invrs_rho_ds_zm, rho_ds_zt,
    rho_ds_zm, l_implemented, nu_vert_res_dep,
    l_diffuse_rtm_and_thlm,
    l_stability_correct_Kh_N2_zm,
    l_upwind_xm_ma,
):
    """Calculate various xm and w'x' terms reused by multiple LHS matrices."""
    constant_nu = 0.1
    Kw6_zm = jnp.maximum(
        zt2zm(nzm, nzt, ngrdcol, gr, Kw6),
        zero_threshold,
    )

    # Calculate turbulent advection terms of xm for all grid levels
    lhs_ta_xm = xm_term_ta_lhs(
        nzm, nzt, ngrdcol, gr,
        rho_ds_zm, invrs_rho_ds_zt,
    )

    # Calculate turbulent production terms of w'x' for all grid level
    lhs_tp = wpxp_term_tp_lhs(nzm, ngrdcol, gr, wp2)

    # Calculate accumulation of w'x' and w'x' pressure term 2 of w'x' for all grid level
    # https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:wpxp_pr
    lhs_ac_pr2 = wpxp_terms_ac_pr2_lhs(
        nzm, nzt, ngrdcol, gr, C7_Skw_fnc,
        wm_zt, gr.invrs_dzm,
    )

    # Calculate diffusion terms for all momentum grid level
    lhs_diff_zm = diffusion_zm_lhs(
        nzm, nzt, ngrdcol, gr, Kw6, Kw6_zm, nu_vert_res_dep.nu6,
        invrs_rho_ds_zm, rho_ds_zt,
    )

    # Calculate mean advection terms for all momentum grid level
    lhs_ma_zm = term_ma_zm_lhs(
        nzm, nzt, ngrdcol, wm_zm,
        gr.invrs_dzm, gr.weights_zm2zt,
    )

    lhs_diff_zt = jnp.zeros((ndiags3, ngrdcol, nzt), dtype=jnp.float64)
    lhs_ma_zt = jnp.zeros((ndiags3, ngrdcol, nzt), dtype=jnp.float64)

    # Calculate diffusion terms for all thermodynamic grid level
    if l_diffuse_rtm_and_thlm:
        if l_stability_correct_Kh_N2_zm:
            Kh_N2_zm = jnp.asarray(Kh_zm) / jnp.asarray(stability_correction)
        else:
            Kh_N2_zm = jnp.asarray(Kh_zm)

        K_zm = Kh_N2_zm + constant_nu
        K_zt = jnp.maximum(
            zm2zt(nzm, nzt, ngrdcol, gr, K_zm),
            zero_threshold,
        )
        zeros_array = jnp.zeros(ngrdcol, dtype=jnp.float64)

        lhs_diff_zt = diffusion_zt_lhs(
            nzm, nzt, ngrdcol, gr, K_zm, K_zt, zeros_array,
            invrs_rho_ds_zt, rho_ds_zm,
        )

    # Calculate mean advection terms for all thermodynamic grid level
    if not l_implemented:
        lhs_ma_zt = term_ma_zt_lhs(
            nzm, nzt, ngrdcol, wm_zt, gr.weights_zt2zm,
            gr.invrs_dzt, gr.invrs_dzm,
            l_upwind_xm_ma, gr.grid_dir,
        )

    return lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm, lhs_tp, lhs_ta_xm, lhs_ac_pr2


def xm_wpxp_rhs(
    nzm, nzt, ngrdcol, gr, solve_type, l_iter, dt,
    xm, wpxp, xm_forcing, wpxp_forcing, C7_Skw_fnc,
    xpthvp, rhs_ta, thv_ds_zm,
    lhs_pr1, lhs_ta_wpxp,
    stats,
):
    """Compute RHS vector for xm and w'x'."""
    l_sample = stats.l_sample
    xm = jnp.asarray(xm, dtype=jnp.float64)
    wpxp = jnp.asarray(wpxp, dtype=jnp.float64)
    xm_forcing = jnp.asarray(xm_forcing, dtype=jnp.float64)
    wpxp_forcing = jnp.asarray(wpxp_forcing, dtype=jnp.float64)
    rhs_ta = jnp.asarray(rhs_ta, dtype=jnp.float64)
    lhs_pr1 = jnp.asarray(lhs_pr1, dtype=jnp.float64)
    lhs_ta_wpxp = jnp.asarray(lhs_ta_wpxp, dtype=jnp.float64)

    invrs_dt = one / dt
    rhs_bp_pr3 = wpxp_terms_bp_pr3_rhs(
        nzm, ngrdcol, gr, C7_Skw_fnc, thv_ds_zm, xpthvp,
    )

    rhs = jnp.zeros((ngrdcol, 2 * nzm - 1), dtype=jnp.float64)

    # Index of the momentum boundary levels
    if gr.grid_dir_indx > 0:
        # Ascending Grid
        rhs_lb_idx_zm = 0
        rhs_ub_idx_zm = 2 * nzm - 2
    else:
        # Descending Grid
        rhs_lb_idx_zm = 2 * nzm - 2
        rhs_ub_idx_zm = 0

    # Set lower boundary for w'x'
    rhs = rhs.at[:, rhs_lb_idx_zm].set(wpxp[:, gr.k_lb_zm])

    # RHS time tendency and forcings for xm
    rhs = rhs.at[:, 1:2 * nzt + 1:2].set(
        xm * invrs_dt + xm_forcing
    )

    k_int = slice(1, nzm - 1)
    k_wpxp = slice(2, 2 * nzm - 2, 2)
    k_plus = slice(1 + gr.grid_dir_indx, nzm - 1 + gr.grid_dir_indx)
    k_minus = slice(1 - gr.grid_dir_indx, nzm - 1 - gr.grid_dir_indx)
    rhs_wpxp_int = (
        rhs_bp_pr3[:, k_int]
        + wpxp_forcing[:, k_int]
        + rhs_ta[:, k_int]
        + (one - gamma_over_implicit_ts)
        * (
            -lhs_ta_wpxp[1 - gr.grid_dir_indx, :, k_int] * wpxp[:, k_plus]
            - lhs_ta_wpxp[1, :, k_int] * wpxp[:, k_int]
            - lhs_ta_wpxp[1 + gr.grid_dir_indx, :, k_int] * wpxp[:, k_minus]
            - lhs_pr1[:, k_int] * wpxp[:, k_int]
        )
    )
    if l_iter:
        rhs_wpxp_int = rhs_wpxp_int + wpxp[:, k_int] * invrs_dt
    rhs = rhs.at[:, k_wpxp].set(rhs_wpxp_int)

    # Upper boundary for w'x'
    rhs = rhs.at[:, rhs_ub_idx_zm].set(zero)

    if l_sample:
        if solve_type == xm_wpxp_rtm:
            name_xm_f = "rtm_forcing"
            name_wpxp_bp = "wprtp_bp"
            name_wpxp_pr3 = "wprtp_pr3"
            name_wpxp_f = "wprtp_forcing"
            name_wpxp_ta = "wprtp_ta"
            name_wpxp_pr1 = "wprtp_pr1"
        elif solve_type == xm_wpxp_thlm:
            name_xm_f = "thlm_forcing"
            name_wpxp_bp = "wpthlp_bp"
            name_wpxp_pr3 = "wpthlp_pr3"
            name_wpxp_f = "wpthlp_forcing"
            name_wpxp_ta = "wpthlp_ta"
            name_wpxp_pr1 = "wpthlp_pr1"
        elif solve_type == xm_wpxp_um:
            name_xm_f = ""
            name_wpxp_bp = "upwp_bp"
            name_wpxp_pr3 = "upwp_pr3"
            name_wpxp_f = ""
            name_wpxp_ta = "upwp_ta"
            name_wpxp_pr1 = "upwp_pr1"
        elif solve_type == xm_wpxp_vm:
            name_xm_f = ""
            name_wpxp_bp = "vpwp_bp"
            name_wpxp_pr3 = "vpwp_pr3"
            name_wpxp_f = ""
            name_wpxp_ta = "vpwp_ta"
            name_wpxp_pr1 = "vpwp_pr1"
        else:
            name_xm_f = ""
            name_wpxp_bp = ""
            name_wpxp_pr3 = ""
            name_wpxp_f = ""
            name_wpxp_ta = ""
            name_wpxp_pr1 = ""

        C7_Skw_fnc_zeros = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        C7_Skw_fnc_plus_one = jnp.asarray(C7_Skw_fnc) + one

        # Statistics: explicit contributions for wpxp.
        rhs_bp = wpxp_terms_bp_pr3_rhs(
            nzm, ngrdcol, gr, C7_Skw_fnc_zeros,
            thv_ds_zm, xpthvp,
        )
        rhs_pr3 = wpxp_terms_bp_pr3_rhs(
            nzm, ngrdcol, gr, C7_Skw_fnc_plus_one,
            thv_ds_zm, xpthvp,
        )

        # Keep the bp/pr3 split consistent with the legacy stats decomposition.
        # rhs_bp uses C7_Skw_fnc=0 and rhs_pr3 uses C7_Skw_fnc+1.
        if len(name_wpxp_bp.strip()) > 0:
            stats = stats.update(name_wpxp_bp, rhs_bp)
        if len(name_wpxp_pr3.strip()) > 0:
            stats = stats.update(name_wpxp_pr3, rhs_pr3)

        if len(name_wpxp_f.strip()) > 0:
            wpxp_forcing_stats = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
            wpxp_forcing_stats = wpxp_forcing_stats.at[:, 1:nzm - 1].set(
                wpxp_forcing[:, 1:nzm - 1]
            )
            # w'x' forcing term is completely explicit; call stat_update_var_pt.
            stats = stats.update(name_wpxp_f, wpxp_forcing_stats)

        if len(name_wpxp_ta.strip()) > 0:
            # <w'x'> term ta has both implicit and explicit components.
            stats = stats.begin_budget(name_wpxp_ta, rhs_ta)
            ta_over = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
            ta_over = ta_over.at[:, 1:nzm - 1].set(
                (one - gamma_over_implicit_ts)
                * (
                    -lhs_ta_wpxp[1 - gr.grid_dir_indx, :, k_int] * wpxp[:, k_plus]
                    - lhs_ta_wpxp[1, :, k_int] * wpxp[:, k_int]
                    - lhs_ta_wpxp[1 + gr.grid_dir_indx, :, k_int] * wpxp[:, k_minus]
                )
            )
            stats = stats.update_budget(name_wpxp_ta, ta_over)

        if len(name_wpxp_pr1.strip()) > 0:
            pr1_over = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
            pr1_over = pr1_over.at[:, 1:nzm - 1].set(
                (one - gamma_over_implicit_ts) * lhs_pr1[:, k_int] * wpxp[:, k_int]
            )
            stats = stats.begin_budget(name_wpxp_pr1, pr1_over)

        # Statistics: explicit contributions for xm
        #             (including microphysics/radiation).
        if len(name_xm_f.strip()) > 0:
            stats = stats.update(name_xm_f, xm_forcing)

    return rhs, stats


def calc_xm_wpxp_ta_terms(
    nzm, nzt, ngrdcol, sclr_dim, gr, wp2rtp,
    wp2thlp, wp2sclrp,
    rho_ds_zt, invrs_rho_ds_zm, rho_ds_zm,
    sigma_sqd_w, wp3_on_wp2_zt,
    pdf_implicit_coefs_terms,
    iiPDF_type,
    l_explicit_turbulent_adv_wpxp, l_predict_upwp_vpwp,
    l_scalar_calc,
    l_godunov_upwind_wpxp_ta,
    stats,
):
    """Calculate the turbulent advection terms for LHS and RHS matrices."""
    l_sample = stats.l_sample
    l_dummy_false = False

    coef_wp2rtp_implicit = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    term_wp2rtp_explicit = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    coef_wp2rtp_implicit_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    term_wp2rtp_explicit_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)

    coef_wp2thlp_implicit = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    term_wp2thlp_explicit = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    coef_wp2thlp_implicit_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    term_wp2thlp_explicit_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)

    term_wp2sclrp_explicit = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    term_wp2sclrp_explicit_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)

    sgn_t_vel_wprtp = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    sgn_t_vel_wpthlp = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    sgn_t_vel_wpsclrp = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)

    lhs_ta_wprtp = jnp.zeros((ndiags3, ngrdcol, nzm), dtype=jnp.float64)
    lhs_ta_wpthlp = jnp.zeros((ndiags3, ngrdcol, nzm), dtype=jnp.float64)
    lhs_ta_wpup = jnp.zeros((ndiags3, ngrdcol, nzm), dtype=jnp.float64)
    lhs_ta_wpvp = jnp.zeros((ndiags3, ngrdcol, nzm), dtype=jnp.float64)
    lhs_ta_wpsclrp = jnp.zeros((ndiags3, ngrdcol, nzm, max(sclr_dim, 1)), dtype=jnp.float64)

    rhs_ta_wprtp = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    rhs_ta_wpthlp = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    rhs_ta_wpup = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    rhs_ta_wpvp = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    rhs_ta_wpsclrp = jnp.zeros((ngrdcol, nzm, max(sclr_dim, 1)), dtype=jnp.float64)

    # Set up the implicit coefficients and explicit terms for turbulent
    # advection of <w'rt'>, <w'thl'>, and <w'sclr'>.
    if l_explicit_turbulent_adv_wpxp:

        # The turbulent advection of <w'x'> is handled explicitly
        term_wp2rtp_explicit = jnp.asarray(wp2rtp)
        term_wp2thlp_explicit = jnp.asarray(wp2thlp)

        # Calculate the RHS turbulent advection term for <w'r_t'>
        rhs_ta_wprtp = xpyp_term_ta_pdf_rhs(
            nzm, nzt, ngrdcol, gr, term_wp2rtp_explicit,
            rho_ds_zt, rho_ds_zm,
            invrs_rho_ds_zm,
            l_dummy_false,
            sgn_t_vel_wprtp,
            term_wp2rtp_explicit_zm,
        )

        # Calculate the RHS turbulent advection term for <w'thl'>
        rhs_ta_wpthlp = xpyp_term_ta_pdf_rhs(
            nzm, nzt, ngrdcol, gr, term_wp2thlp_explicit,
            rho_ds_zt, rho_ds_zm,
            invrs_rho_ds_zm,
            l_dummy_false,
            sgn_t_vel_wpthlp,
            term_wp2thlp_explicit_zm,
        )

        for sclr in range(sclr_dim):
            term_wp2sclrp_explicit = jnp.asarray(wp2sclrp)[:, :, sclr]
            rhs_ta_wpsclrp = rhs_ta_wpsclrp.at[:, :, sclr].set(xpyp_term_ta_pdf_rhs(
                nzm, nzt, ngrdcol, gr, term_wp2sclrp_explicit,
                rho_ds_zt, rho_ds_zm,
                invrs_rho_ds_zm,
                l_dummy_false,
                sgn_t_vel_wpsclrp,
                term_wp2sclrp_explicit_zm,
            ))

    else:

        # The turbulent advection of <w'x'> is handled implicitly or
        # semi-implicitly.
        if iiPDF_type == iiPDF_ADG1:

            # The ADG1 PDF is used.
            a1_coef = one / (one - jnp.asarray(sigma_sqd_w))
            a1_coef_zt = jnp.maximum(
                zm2zt(nzm, nzt, ngrdcol, gr, a1_coef),
                zero_threshold,
            )
            coef_wp2rtp_implicit = a1_coef_zt * jnp.asarray(wp3_on_wp2_zt)
            coef_wp2thlp_implicit = coef_wp2rtp_implicit

            if not l_godunov_upwind_wpxp_ta:
                # Calculate the LHS turbulent advection term for <w'r_t'>
                lhs_ta_wprtp = xpyp_term_ta_pdf_lhs(
                    nzm, nzt, ngrdcol, gr, coef_wp2rtp_implicit,
                    rho_ds_zt, rho_ds_zm,
                    invrs_rho_ds_zm,
                    l_dummy_false,
                    sgn_t_vel_wprtp,
                    coef_wp2rtp_implicit_zm,
                )
            else:
                # Godunov-like method for the vertical discretization of ta term
                coef_wp2rtp_implicit = a1_coef_zt * jnp.asarray(wp3_on_wp2_zt)
                coef_wp2thlp_implicit = coef_wp2rtp_implicit
                lhs_ta_wprtp = xpyp_term_ta_pdf_lhs_godunov(
                    nzm, nzt, ngrdcol, gr,
                    coef_wp2rtp_implicit,
                    invrs_rho_ds_zm, rho_ds_zm,
                )

            # For ADG1, the LHS turbulent advection terms for
            # <w'r_t'>, <w'thl'>, <w'sclr'> are all equal
            lhs_ta_wpthlp = jnp.asarray(lhs_ta_wprtp)

            if l_scalar_calc:
                for sclr in range(sclr_dim):
                    lhs_ta_wpsclrp = lhs_ta_wpsclrp.at[:, :, :, sclr].set(lhs_ta_wprtp)

            # The <w'r_t'>, <w'thl'>, <w'sclr'> turbulent advection terms are entirely implicit.
            if l_predict_upwp_vpwp:

                # Predict <u> and <u'w'>, as well as <v> and <v'w'>.
                # These terms are equal to the <w'r_t'> terms as well in this case
                lhs_ta_wpup = jnp.asarray(lhs_ta_wprtp)
                lhs_ta_wpvp = jnp.asarray(lhs_ta_wprtp)

        elif iiPDF_type == iiPDF_new:

            # The new PDF is used.
            coef_wp2rtp_implicit = jnp.asarray(pdf_implicit_coefs_terms.coef_wp2rtp_implicit)
            coef_wp2thlp_implicit = jnp.asarray(pdf_implicit_coefs_terms.coef_wp2thlp_implicit)
            term_wp2rtp_explicit = jnp.asarray(pdf_implicit_coefs_terms.term_wp2rtp_explicit)
            term_wp2thlp_explicit = jnp.asarray(pdf_implicit_coefs_terms.term_wp2thlp_explicit)

            # Calculate the LHS turbulent advection term for <w'rt'>
            lhs_ta_wprtp = xpyp_term_ta_pdf_lhs(
                nzm, nzt, ngrdcol, gr, coef_wp2rtp_implicit,
                rho_ds_zt, rho_ds_zm,
                invrs_rho_ds_zm,
                l_dummy_false,
                sgn_t_vel_wprtp,
                coef_wp2rtp_implicit_zm,
            )

            # Calculate the RHS turbulent advection term for <w'rt'>
            rhs_ta_wprtp = xpyp_term_ta_pdf_rhs(
                nzm, nzt, ngrdcol, gr, term_wp2rtp_explicit,
                rho_ds_zt, rho_ds_zm,
                invrs_rho_ds_zm,
                l_dummy_false,
                sgn_t_vel_wprtp,
                term_wp2rtp_explicit_zm,
            )

            # Calculate the LHS turbulent advection term for <w'thl'>
            lhs_ta_wpthlp = xpyp_term_ta_pdf_lhs(
                nzm, nzt, ngrdcol, gr, coef_wp2thlp_implicit,
                rho_ds_zt, rho_ds_zm,
                invrs_rho_ds_zm,
                l_dummy_false,
                sgn_t_vel_wpthlp,
                coef_wp2thlp_implicit_zm,
            )

            # Calculate the RHS turbulent advection term for <w'thl'>
            rhs_ta_wpthlp = xpyp_term_ta_pdf_rhs(
                nzm, nzt, ngrdcol, gr, term_wp2thlp_explicit,
                rho_ds_zt, rho_ds_zm,
                invrs_rho_ds_zm,
                l_dummy_false,
                sgn_t_vel_wpthlp,
                term_wp2thlp_explicit_zm,
            )

        elif iiPDF_type == iiPDF_new_hybrid:

            # The new hybrid PDF is used.
            coef_wp2rtp_implicit = jnp.asarray(pdf_implicit_coefs_terms.coef_wp2rtp_implicit)
            coef_wp2thlp_implicit = coef_wp2rtp_implicit

            # Calculate the LHS turbulent advection term for <w'rt'>
            lhs_ta_wprtp = xpyp_term_ta_pdf_lhs(
                nzm, nzt, ngrdcol, gr, coef_wp2rtp_implicit,
                rho_ds_zt, rho_ds_zm,
                invrs_rho_ds_zm,
                l_dummy_false,
                sgn_t_vel_wprtp,
                coef_wp2rtp_implicit_zm,
            )

            # For the new hybrid PDF, the LHS turbulent advection terms for
            # <w'r_t'>, <w'thl'>, and <w'sclr'> are all the same.
            lhs_ta_wpthlp = jnp.asarray(lhs_ta_wprtp)

            if l_scalar_calc:
                for sclr in range(sclr_dim):
                    lhs_ta_wpsclrp = lhs_ta_wpsclrp.at[:, :, :, sclr].set(lhs_ta_wprtp)

            if l_predict_upwp_vpwp:

                # Predict <u> and <u'w'>, as well as <v> and <v'w'>.
                # These terms are equal to the <w'r_t'> terms as well in this case
                lhs_ta_wpup = jnp.asarray(lhs_ta_wprtp)
                lhs_ta_wpvp = jnp.asarray(lhs_ta_wprtp)

    if l_sample:
        stats = stats.update("coef_wp2rtp_implicit", coef_wp2rtp_implicit)
        stats = stats.update("term_wp2rtp_explicit", term_wp2rtp_explicit)
        stats = stats.update("coef_wp2thlp_implicit", coef_wp2thlp_implicit)
        stats = stats.update("term_wp2thlp_explicit", term_wp2thlp_explicit)

    return (
        lhs_ta_wprtp, lhs_ta_wpthlp, lhs_ta_wpup,
        lhs_ta_wpvp, lhs_ta_wpsclrp,
        rhs_ta_wprtp, rhs_ta_wpthlp, rhs_ta_wpup,
        rhs_ta_wpvp, rhs_ta_wpsclrp, stats,
    )


def solve_xm_wpxp_with_single_lhs(
    nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, dt,
    l_iter, nrhs,
    wm_zt, wp2, invrs_tau_C6_zm, tau_max_zm,
    rtpthvp, rtm_forcing, wprtp_forcing, thlpthvp,
    thlm_forcing, wpthlp_forcing, rho_ds_zm,
    rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,
    thv_ds_zm, rtp2, thlp2, l_implemented,
    sclrpthvp, sclrm_forcing, sclrp2, um_forcing,
    vm_forcing, ug, vg, uprcp, vprcp, rc_coef_zm, fcor,
    fcor_y, up2, vp2,
    low_lev_effect, high_lev_effect,
    C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc,
    lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm,
    lhs_ta_wpxp,
    rhs_ta_wprtp, rhs_ta_wpthlp, rhs_ta_wpup,
    rhs_ta_wpvp, rhs_ta_wpsclrp,
    lhs_tp, lhs_ta_xm, lhs_ac_pr2, lhs_pr1_wprtp,
    lhs_pr1_wpthlp, lhs_pr1_wpsclrp,
    C_uu_shr,
    penta_solve_method,
    tridiag_solve_method,
    fill_holes_type,
    l_predict_upwp_vpwp,
    l_ho_nontrad_coriolis,
    l_ho_trad_coriolis,
    l_diffuse_rtm_and_thlm,
    l_upwind_xm_ma,
    l_enable_relaxed_clipping,
    l_perturbed_wind,
    l_mono_flux_lim_thlm,
    l_mono_flux_lim_rtm,
    l_mono_flux_lim_um,
    l_mono_flux_lim_vm,
    l_mono_flux_lim_spikefix,
    wprtp_cl_num,
    wpthlp_cl_num, upwp_cl_num, vpwp_cl_num,
    stats,
    rtm, wprtp, thlm, wpthlp,
    sclrm, wpsclrp, um, upwp, vm, vpwp,
    um_pert, vm_pert, upwp_pert, vpwp_pert, err_info,
):
    """Solve all xm_wpxp fields with one shared LHS matrix."""
    l_sample = stats.l_sample
    lhs = xm_wpxp_lhs(
        nzm, nzt, ngrdcol, l_iter, dt,
        l_implemented, lhs_diff_zm, lhs_diff_zt,
        lhs_ma_zm, lhs_ma_zt, lhs_ta_wpxp, lhs_ta_xm,
        lhs_tp, lhs_pr1_wprtp, lhs_ac_pr2,
        l_diffuse_rtm_and_thlm,
    )

    rhs = jnp.zeros((ngrdcol, 2 * nzm - 1, nrhs), dtype=jnp.float64)
    rhs_field, stats = xm_wpxp_rhs(
        nzm, nzt, ngrdcol, gr, xm_wpxp_rtm, l_iter, dt,
        rtm, wprtp, rtm_forcing, wprtp_forcing, C7_Skw_fnc,
        rtpthvp, rhs_ta_wprtp, thv_ds_zm, lhs_pr1_wprtp,
        lhs_ta_wpxp, stats,
    )
    rhs = rhs.at[:, :, 0].set(rhs_field)
    rhs_field, stats = xm_wpxp_rhs(
        nzm, nzt, ngrdcol, gr, xm_wpxp_thlm, l_iter, dt,
        thlm, wpthlp, thlm_forcing, wpthlp_forcing, C7_Skw_fnc,
        thlpthvp, rhs_ta_wpthlp, thv_ds_zm, lhs_pr1_wpthlp,
        lhs_ta_wpxp, stats,
    )
    rhs = rhs.at[:, :, 1].set(rhs_field)

    wpsclrp_forcing = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    for sclr in range(sclr_dim):
        rhs_field, stats = xm_wpxp_rhs(
            nzm, nzt, ngrdcol, gr, xm_wpxp_scalar, l_iter, dt,
            sclrm[:, :, sclr], wpsclrp[:, :, sclr], sclrm_forcing[:, :, sclr],
            wpsclrp_forcing, C7_Skw_fnc,
            sclrpthvp[:, :, sclr], rhs_ta_wpsclrp[:, :, sclr], thv_ds_zm,
            lhs_pr1_wpsclrp, lhs_ta_wpxp, stats,
        )
        rhs = rhs.at[:, :, 2 + sclr].set(rhs_field)

    um_tndcy = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    vm_tndcy = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)

    if l_predict_upwp_vpwp:
        fcor_col = jnp.asarray(fcor)[:, None]
        fcor_y_col = jnp.asarray(fcor_y)[:, None]

        if not l_implemented:
            um_tndcy = jnp.asarray(um_forcing) - fcor_col * (jnp.asarray(vg) - jnp.asarray(vm))
            vm_tndcy = jnp.asarray(vm_forcing) + fcor_col * (jnp.asarray(ug) - jnp.asarray(um))

            if l_sample:
                stats = stats.update("um_gf", -fcor_col * jnp.asarray(vg))
                stats = stats.update("vm_gf", fcor_col * jnp.asarray(ug))
                stats = stats.update("um_cf", fcor_col * jnp.asarray(vm))
                stats = stats.update("vm_cf", -fcor_col * jnp.asarray(um))
                stats = stats.update("um_f", um_forcing)
                stats = stats.update("vm_f", vm_forcing)

        ddzt_um = ddzt(nzm, nzt, ngrdcol, gr, um)
        ddzt_vm = ddzt(nzm, nzt, ngrdcol, gr, vm)

        C_uu_shr_col = jnp.asarray(C_uu_shr)[:, None]
        upwp_forcing = C_uu_shr_col * jnp.asarray(wp2) * ddzt_um
        vpwp_forcing = C_uu_shr_col * jnp.asarray(wp2) * ddzt_vm

        if l_ho_trad_coriolis:
            upwp_forcing = upwp_forcing + fcor_col * jnp.asarray(vpwp)
            vpwp_forcing = vpwp_forcing - fcor_col * jnp.asarray(upwp)

        if l_ho_nontrad_coriolis:
            upwp_forcing = upwp_forcing + fcor_y_col * (jnp.asarray(up2) - jnp.asarray(wp2))

        upwp_forcing_pert = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        vpwp_forcing_pert = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        if l_perturbed_wind:
            ddzt_um_pert = ddzt(nzm, nzt, ngrdcol, gr, um_pert)
            ddzt_vm_pert = ddzt(nzm, nzt, ngrdcol, gr, vm_pert)
            upwp_forcing_pert = C_uu_shr_col * jnp.asarray(wp2) * ddzt_um_pert
            vpwp_forcing_pert = C_uu_shr_col * jnp.asarray(wp2) * ddzt_vm_pert

        if l_sample:
            stats = stats.update("upwp_pr4", C_uu_shr_col * jnp.asarray(wp2) * ddzt_um)
            stats = stats.update("vpwp_pr4", C_uu_shr_col * jnp.asarray(wp2) * ddzt_vm)
            if l_ho_trad_coriolis:
                stats = stats.update("upwp_tct", fcor_col * jnp.asarray(vpwp))
                stats = stats.update("vpwp_tct", -fcor_col * jnp.asarray(upwp))
            if l_ho_nontrad_coriolis:
                stats = stats.update("upwp_nct", fcor_y_col * (jnp.asarray(up2) - jnp.asarray(wp2)))

        tau_C6_zm = jnp.minimum(one / jnp.asarray(invrs_tau_C6_zm), jnp.asarray(tau_max_zm))
        um_smth = zt2zm2zt(nzm, nzt, ngrdcol, gr, um, zt_min=-jnp.inf)
        vm_smth = zt2zm2zt(nzm, nzt, ngrdcol, gr, vm, zt_min=-jnp.inf)

        upthlp = diagnose_upxp(
            nzm, nzt, ngrdcol, gr, upwp, thlm, wpthlp, um_smth,
            C6thl_Skw_fnc, tau_C6_zm, C7_Skw_fnc,
        )
        uprtp = diagnose_upxp(
            nzm, nzt, ngrdcol, gr, upwp, rtm, wprtp, um_smth,
            C6rt_Skw_fnc, tau_C6_zm, C7_Skw_fnc,
        )
        vpthlp = diagnose_upxp(
            nzm, nzt, ngrdcol, gr, vpwp, thlm, wpthlp, vm_smth,
            C6thl_Skw_fnc, tau_C6_zm, C7_Skw_fnc,
        )
        vprtp = diagnose_upxp(
            nzm, nzt, ngrdcol, gr, vpwp, rtm, wprtp, vm_smth,
            C6rt_Skw_fnc, tau_C6_zm, C7_Skw_fnc,
        )

        if l_perturbed_wind:
            upthlp_pert = diagnose_upxp(
                nzm, nzt, ngrdcol, gr, upwp_pert, thlm, wpthlp, um_pert,
                C6thl_Skw_fnc, tau_C6_zm, C7_Skw_fnc,
            )
            uprtp_pert = diagnose_upxp(
                nzm, nzt, ngrdcol, gr, upwp_pert, rtm, wprtp, um_pert,
                C6rt_Skw_fnc, tau_C6_zm, C7_Skw_fnc,
            )
            vpthlp_pert = diagnose_upxp(
                nzm, nzt, ngrdcol, gr, vpwp_pert, thlm, wpthlp, vm_pert,
                C6thl_Skw_fnc, tau_C6_zm, C7_Skw_fnc,
            )
            vprtp_pert = diagnose_upxp(
                nzm, nzt, ngrdcol, gr, vpwp_pert, rtm, wprtp, vm_pert,
                C6rt_Skw_fnc, tau_C6_zm, C7_Skw_fnc,
            )
        else:
            upthlp_pert = uprtp_pert = vpthlp_pert = vprtp_pert = None

        upthvp_tmp = upthlp + ep1 * jnp.asarray(thv_ds_zm) * uprtp + jnp.asarray(rc_coef_zm) * jnp.asarray(uprcp)
        vpthvp_tmp = vpthlp + ep1 * jnp.asarray(thv_ds_zm) * vprtp + jnp.asarray(rc_coef_zm) * jnp.asarray(vprcp)
        upthvp = zm2zt2zm(nzm, nzt, ngrdcol, gr, upthvp_tmp, zm_min=-jnp.inf)
        vpthvp = zm2zt2zm(nzm, nzt, ngrdcol, gr, vpthvp_tmp, zm_min=-jnp.inf)

        if l_perturbed_wind:
            upthvp_pert = (
                upthlp_pert
                + ep1 * jnp.asarray(thv_ds_zm) * uprtp_pert
                + jnp.asarray(rc_coef_zm) * jnp.asarray(uprcp)
            )
            vpthvp_pert = (
                vpthlp_pert
                + ep1 * jnp.asarray(thv_ds_zm) * vprtp_pert
                + jnp.asarray(rc_coef_zm) * jnp.asarray(vprcp)
            )

        if l_sample:
            stats = stats.update("upthlp", upthlp)
            stats = stats.update("uprtp", uprtp)
            stats = stats.update("vpthlp", vpthlp)
            stats = stats.update("vprtp", vprtp)
            stats = stats.update("upthvp", upthvp)
            stats = stats.update("vpthvp", vpthvp)

        rhs_field, stats = xm_wpxp_rhs(
            nzm, nzt, ngrdcol, gr, xm_wpxp_um, l_iter, dt,
            um, upwp, um_tndcy, upwp_forcing, C7_Skw_fnc,
            upthvp, rhs_ta_wpup, thv_ds_zm,
            lhs_pr1_wprtp, lhs_ta_wpxp, stats,
        )
        rhs = rhs.at[:, :, 2 + sclr_dim].set(rhs_field)
        rhs_field, stats = xm_wpxp_rhs(
            nzm, nzt, ngrdcol, gr, xm_wpxp_vm, l_iter, dt,
            vm, vpwp, vm_tndcy, vpwp_forcing, C7_Skw_fnc,
            vpthvp, rhs_ta_wpvp, thv_ds_zm,
            lhs_pr1_wprtp, lhs_ta_wpxp, stats,
        )
        rhs = rhs.at[:, :, 3 + sclr_dim].set(rhs_field)

        if l_perturbed_wind:
            rhs_field, stats = xm_wpxp_rhs(
                nzm, nzt, ngrdcol, gr, xm_wpxp_um, l_iter, dt,
                um_pert, upwp_pert, um_tndcy, upwp_forcing_pert, C7_Skw_fnc,
                upthvp_pert, rhs_ta_wpup, thv_ds_zm,
                lhs_pr1_wprtp, lhs_ta_wpxp, stats,
            )
            rhs = rhs.at[:, :, 4 + sclr_dim].set(rhs_field)
            rhs_field, stats = xm_wpxp_rhs(
                nzm, nzt, ngrdcol, gr, xm_wpxp_vm, l_iter, dt,
                vm_pert, vpwp_pert, vm_tndcy, vpwp_forcing_pert, C7_Skw_fnc,
                vpthvp_pert, rhs_ta_wpvp, thv_ds_zm,
                lhs_pr1_wprtp, lhs_ta_wpxp, stats,
            )
            rhs = rhs.at[:, :, 5 + sclr_dim].set(rhs_field)

    field_pairs = [(rtm, wprtp), (thlm, wpthlp)]
    for sclr in range(sclr_dim):
        field_pairs.append((sclrm[:, :, sclr], wpsclrp[:, :, sclr]))
    if l_predict_upwp_vpwp:
        field_pairs.extend([(um, upwp), (vm, vpwp)])
        if l_perturbed_wind:
            field_pairs.extend([(um_pert, upwp_pert), (vm_pert, vpwp_pert)])

    old_solution = jnp.zeros((ngrdcol, 2 * nzm - 1, nrhs), dtype=jnp.float64)
    if penta_solve_method == penta_bicgstab:
        for irhs, (xm_field, wpxp_field) in enumerate(field_pairs):
            old_solution = old_solution.at[:, 1:2 * (nzm - 1) + 1:2, irhs].set(
                jnp.asarray(xm_field)
            )
            old_solution = old_solution.at[:, 0:2 * nzm - 1:2, irhs].set(
                jnp.asarray(wpxp_field)
            )

    use_rcond = bool(
        l_sample
        and (
            stats.var_on_stats_list("thlm_matrix_condt_num")
            or stats.var_on_stats_list("rtm_matrix_condt_num")
        )
    )
    solution, rcond, err_info = xm_wpxp_solve(
        nzm, ngrdcol, nrhs, gr, penta_solve_method, l_implemented,
        old_solution, lhs, rhs, err_info, use_rcond=use_rcond,
    )

    rtm, wprtp, err_info, wprtp_cl_num, stats = xm_wpxp_clipping_and_stats(
        nzm, nzt, ngrdcol, gr, xm_wpxp_rtm, dt, wp2, rtp2, wm_zt,
        rtm_forcing, rho_ds_zm, rho_ds_zt,
        invrs_rho_ds_zm, invrs_rho_ds_zt,
        rt_tol ** 2, rt_tol, rcond,
        low_lev_effect, high_lev_effect,
        lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,
        lhs_diff_zm, C7_Skw_fnc,
        lhs_tp, lhs_ta_xm, lhs_pr1_wprtp,
        l_implemented, solution[:, :, 0],
        tridiag_solve_method, fill_holes_type,
        l_upwind_xm_ma, l_enable_relaxed_clipping,
        l_mono_flux_lim_thlm, l_mono_flux_lim_rtm,
        l_mono_flux_lim_um, l_mono_flux_lim_vm,
        l_mono_flux_lim_spikefix,
        stats,
        rtm, rt_tol_mfl, wprtp, err_info,
        wpxp_cl_num=wprtp_cl_num,
    )

    thlm, wpthlp, err_info, wpthlp_cl_num, stats = xm_wpxp_clipping_and_stats(
        nzm, nzt, ngrdcol, gr, xm_wpxp_thlm, dt, wp2, thlp2, wm_zt,
        thlm_forcing, rho_ds_zm, rho_ds_zt,
        invrs_rho_ds_zm, invrs_rho_ds_zt,
        thl_tol ** 2, thl_tol, rcond,
        low_lev_effect, high_lev_effect,
        lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,
        lhs_diff_zm, C7_Skw_fnc,
        lhs_tp, lhs_ta_xm, lhs_pr1_wprtp,
        l_implemented, solution[:, :, 1],
        tridiag_solve_method, fill_holes_type,
        l_upwind_xm_ma, l_enable_relaxed_clipping,
        l_mono_flux_lim_thlm, l_mono_flux_lim_rtm,
        l_mono_flux_lim_um, l_mono_flux_lim_vm,
        l_mono_flux_lim_spikefix,
        stats,
        thlm, thl_tol_mfl, wpthlp, err_info,
        wpxp_cl_num=wpthlp_cl_num,
    )

    for sclr in range(sclr_dim):
        sclrm_s, wpsclrp_s, err_info, _wpxp_cl_num, stats = xm_wpxp_clipping_and_stats(
            nzm, nzt, ngrdcol, gr, xm_wpxp_scalar, dt, wp2, sclrp2[:, :, sclr], wm_zt,
            sclrm_forcing[:, :, sclr], rho_ds_zm, rho_ds_zt,
            invrs_rho_ds_zm, invrs_rho_ds_zt,
            sclr_tol[sclr] ** 2, sclr_tol[sclr], rcond,
            low_lev_effect, high_lev_effect,
            lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,
            lhs_diff_zm, C7_Skw_fnc,
            lhs_tp, lhs_ta_xm, lhs_pr1_wprtp,
            l_implemented, solution[:, :, 2 + sclr],
            tridiag_solve_method, fill_holes_type,
            l_upwind_xm_ma, l_enable_relaxed_clipping,
            l_mono_flux_lim_thlm, l_mono_flux_lim_rtm,
            l_mono_flux_lim_um, l_mono_flux_lim_vm,
            l_mono_flux_lim_spikefix,
            stats,
            sclrm[:, :, sclr], sclr_tol[sclr], wpsclrp[:, :, sclr], err_info,
        )
        sclrm = sclrm.at[:, :, sclr].set(sclrm_s)
        wpsclrp = wpsclrp.at[:, :, sclr].set(wpsclrp_s)

    if l_predict_upwp_vpwp:
        um, upwp, err_info, upwp_cl_num, stats = xm_wpxp_clipping_and_stats(
            nzm, nzt, ngrdcol, gr, xm_wpxp_um, dt, wp2, up2, wm_zt,
            um_tndcy, rho_ds_zm, rho_ds_zt,
            invrs_rho_ds_zm, invrs_rho_ds_zt,
            w_tol_sqd, w_tol, rcond,
            low_lev_effect, high_lev_effect,
            lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,
            lhs_diff_zm, C7_Skw_fnc,
            lhs_tp, lhs_ta_xm, lhs_pr1_wprtp,
            l_implemented, solution[:, :, 2 + sclr_dim],
            tridiag_solve_method, fill_holes_type,
            l_upwind_xm_ma, l_enable_relaxed_clipping,
            l_mono_flux_lim_thlm, l_mono_flux_lim_rtm,
            l_mono_flux_lim_um, l_mono_flux_lim_vm,
            l_mono_flux_lim_spikefix,
            stats,
            um, w_tol, upwp, err_info,
            wpxp_cl_num=upwp_cl_num,
        )

        vm, vpwp, err_info, vpwp_cl_num, stats = xm_wpxp_clipping_and_stats(
            nzm, nzt, ngrdcol, gr, xm_wpxp_vm, dt, wp2, vp2, wm_zt,
            vm_tndcy, rho_ds_zm, rho_ds_zt,
            invrs_rho_ds_zm, invrs_rho_ds_zt,
            w_tol_sqd, w_tol, rcond,
            low_lev_effect, high_lev_effect,
            lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,
            lhs_diff_zm, C7_Skw_fnc,
            lhs_tp, lhs_ta_xm, lhs_pr1_wprtp,
            l_implemented, solution[:, :, 3 + sclr_dim],
            tridiag_solve_method, fill_holes_type,
            l_upwind_xm_ma, l_enable_relaxed_clipping,
            l_mono_flux_lim_thlm, l_mono_flux_lim_rtm,
            l_mono_flux_lim_um, l_mono_flux_lim_vm,
            l_mono_flux_lim_spikefix,
            stats,
            vm, w_tol, vpwp, err_info,
            wpxp_cl_num=vpwp_cl_num,
        )

        if l_perturbed_wind:
            um_pert, upwp_pert, err_info, _wpxp_cl_num, stats = xm_wpxp_clipping_and_stats(
                nzm, nzt, ngrdcol, gr, xm_wpxp_um, dt, wp2, up2, wm_zt,
                um_tndcy, rho_ds_zm, rho_ds_zt,
                invrs_rho_ds_zm, invrs_rho_ds_zt,
                w_tol_sqd, w_tol, rcond,
                low_lev_effect, high_lev_effect,
                lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,
                lhs_diff_zm, C7_Skw_fnc,
                lhs_tp, lhs_ta_xm, lhs_pr1_wprtp,
                l_implemented, solution[:, :, 4 + sclr_dim],
                tridiag_solve_method, fill_holes_type,
                l_upwind_xm_ma, l_enable_relaxed_clipping,
                l_mono_flux_lim_thlm, l_mono_flux_lim_rtm,
                l_mono_flux_lim_um, l_mono_flux_lim_vm,
                l_mono_flux_lim_spikefix,
                stats,
                um_pert, w_tol, upwp_pert, err_info,
            )
            vm_pert, vpwp_pert, err_info, _wpxp_cl_num, stats = xm_wpxp_clipping_and_stats(
                nzm, nzt, ngrdcol, gr, xm_wpxp_vm, dt, wp2, vp2, wm_zt,
                vm_tndcy, rho_ds_zm, rho_ds_zt,
                invrs_rho_ds_zm, invrs_rho_ds_zt,
                w_tol_sqd, w_tol, rcond,
                low_lev_effect, high_lev_effect,
                lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,
                lhs_diff_zm, C7_Skw_fnc,
                lhs_tp, lhs_ta_xm, lhs_pr1_wprtp,
                l_implemented, solution[:, :, 5 + sclr_dim],
                tridiag_solve_method, fill_holes_type,
                l_upwind_xm_ma, l_enable_relaxed_clipping,
                l_mono_flux_lim_thlm, l_mono_flux_lim_rtm,
                l_mono_flux_lim_um, l_mono_flux_lim_vm,
                l_mono_flux_lim_spikefix,
                stats,
                vm_pert, w_tol, vpwp_pert, err_info,
            )

    return (
        wprtp_cl_num, wpthlp_cl_num, upwp_cl_num, vpwp_cl_num,
        rtm, wprtp, thlm, wpthlp, sclrm, wpsclrp, um, upwp, vm, vpwp,
        um_pert, vm_pert, upwp_pert, vpwp_pert, err_info, stats,
    )


def solve_xm_wpxp_with_multiple_lhs(
    nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, dt,
    l_iter, nrhs, wm_zt, wp2,
    rtpthvp, rtm_forcing, wprtp_forcing, thlpthvp,
    thlm_forcing, wpthlp_forcing, rho_ds_zm,
    rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,
    thv_ds_zm, rtp2, thlp2, l_implemented,
    sclrpthvp, sclrm_forcing, sclrp2,
    low_lev_effect, high_lev_effect, C7_Skw_fnc,
    lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm,
    lhs_ta_wprtp, lhs_ta_wpthlp, lhs_ta_wpsclrp,
    rhs_ta_wprtp, rhs_ta_wpthlp, rhs_ta_wpsclrp,
    lhs_tp, lhs_ta_xm, lhs_ac_pr2, lhs_pr1_wprtp,
    lhs_pr1_wpthlp, lhs_pr1_wpsclrp,
    penta_solve_method,
    tridiag_solve_method,
    fill_holes_type,
    l_diffuse_rtm_and_thlm,
    l_upwind_xm_ma,
    l_enable_relaxed_clipping,
    l_mono_flux_lim_thlm,
    l_mono_flux_lim_rtm,
    l_mono_flux_lim_um,
    l_mono_flux_lim_vm,
    l_mono_flux_lim_spikefix,
    wprtp_cl_num, wpthlp_cl_num,
    stats,
    rtm, wprtp, thlm, wpthlp, sclrm, wpsclrp, err_info,
):
    """Solve xm_wpxp when LHS matrices differ by field."""
    del nrhs
    l_sample = stats.l_sample
    rtm, wprtp, err_info, wprtp_cl_num, stats = _solve_one_xm_wpxp(
        nzm, nzt, ngrdcol, gr, xm_wpxp_rtm, dt, l_iter, wm_zt, wp2, rtp2,
        rtm_forcing, wprtp_forcing, C7_Skw_fnc, rtpthvp, rhs_ta_wprtp,
        thv_ds_zm, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,
        rt_tol ** 2, rt_tol, rt_tol_mfl, low_lev_effect, high_lev_effect,
        lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm, lhs_ta_wprtp,
        lhs_tp, lhs_ta_xm, lhs_ac_pr2, lhs_pr1_wprtp, l_implemented,
        penta_solve_method, tridiag_solve_method, fill_holes_type,
        l_diffuse_rtm_and_thlm, l_upwind_xm_ma,
        l_enable_relaxed_clipping, l_mono_flux_lim_thlm, l_mono_flux_lim_rtm,
        l_mono_flux_lim_um, l_mono_flux_lim_vm, l_mono_flux_lim_spikefix,
        stats, rtm, wprtp, err_info, wpxp_cl_num=wprtp_cl_num,
    )

    thlm, wpthlp, err_info, wpthlp_cl_num, stats = _solve_one_xm_wpxp(
        nzm, nzt, ngrdcol, gr, xm_wpxp_thlm, dt, l_iter, wm_zt, wp2, thlp2,
        thlm_forcing, wpthlp_forcing, C7_Skw_fnc, thlpthvp, rhs_ta_wpthlp,
        thv_ds_zm, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,
        thl_tol ** 2, thl_tol, thl_tol_mfl, low_lev_effect, high_lev_effect,
        lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm, lhs_ta_wpthlp,
        lhs_tp, lhs_ta_xm, lhs_ac_pr2, lhs_pr1_wpthlp, l_implemented,
        penta_solve_method, tridiag_solve_method, fill_holes_type,
        l_diffuse_rtm_and_thlm, l_upwind_xm_ma,
        l_enable_relaxed_clipping, l_mono_flux_lim_thlm, l_mono_flux_lim_rtm,
        l_mono_flux_lim_um, l_mono_flux_lim_vm, l_mono_flux_lim_spikefix,
        stats, thlm, wpthlp, err_info, wpxp_cl_num=wpthlp_cl_num,
    )

    wpsclrp_forcing = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    for sclr in range(sclr_dim):
        sclrm_s, wpsclrp_s, err_info, _wpxp_cl_num, stats = _solve_one_xm_wpxp(
            nzm, nzt, ngrdcol, gr, xm_wpxp_scalar, dt, l_iter, wm_zt, wp2, sclrp2[:, :, sclr],
            sclrm_forcing[:, :, sclr], wpsclrp_forcing, C7_Skw_fnc,
            sclrpthvp[:, :, sclr], rhs_ta_wpsclrp[:, :, sclr],
            thv_ds_zm, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,
            sclr_tol[sclr] ** 2, sclr_tol[sclr], sclr_tol[sclr],
            low_lev_effect, high_lev_effect,
            lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm, lhs_ta_wpsclrp[:, :, :, sclr],
            lhs_tp, lhs_ta_xm, lhs_ac_pr2, lhs_pr1_wpsclrp, l_implemented,
            penta_solve_method, tridiag_solve_method, fill_holes_type,
            l_diffuse_rtm_and_thlm, l_upwind_xm_ma,
            l_enable_relaxed_clipping, l_mono_flux_lim_thlm, l_mono_flux_lim_rtm,
            l_mono_flux_lim_um, l_mono_flux_lim_vm, l_mono_flux_lim_spikefix,
            stats, sclrm[:, :, sclr], wpsclrp[:, :, sclr], err_info,
        )
        sclrm = sclrm.at[:, :, sclr].set(sclrm_s)
        wpsclrp = wpsclrp.at[:, :, sclr].set(wpsclrp_s)

    return wprtp_cl_num, wpthlp_cl_num, rtm, wprtp, thlm, wpthlp, sclrm, wpsclrp, err_info, stats


def xm_wpxp_solve(
    nzm, ngrdcol, nrhs, gr,
    penta_solve_method, l_implemented,
    old_solution, lhs, rhs, err_info,
    use_rcond=False,
):
    """Solve for xm / w'x' using the band diagonal solver."""
    ndim = 2 * nzm - 1
    lhs_work = jnp.asarray(lhs, dtype=jnp.float64)
    rhs_work = jnp.asarray(rhs, dtype=jnp.float64)
    old_work = jnp.asarray(old_solution, dtype=jnp.float64)

    if l_force_descending_solves and gr.grid_dir_indx > 0:
        lhs_work = lhs_work[::-1, :, ::-1]
        rhs_work = rhs_work[:, ::-1, :]
        if penta_solve_method == penta_bicgstab:
            old_work = old_work[:, ::-1, :]

    err_info, solution, rcond = band_solve(
        "xm_wpxp", penta_solve_method, ngrdcol, nsup,
        nsub, ndim, nrhs, l_implemented,
        lhs_work, rhs_work, err_info,
        old_soln=old_work, use_rcond=use_rcond,
    )

    if l_force_descending_solves and gr.grid_dir_indx > 0:
        solution = jnp.asarray(solution)[:, ::-1, :]

    return solution, jnp.asarray(rcond, dtype=jnp.float64), err_info


def xm_wpxp_clipping_and_stats(
    nzm, nzt, ngrdcol, gr, solve_type, dt, wp2, xp2, wm_zt,
    xm_forcing, rho_ds_zm, rho_ds_zt,
    invrs_rho_ds_zm, invrs_rho_ds_zt,
    xp2_threshold, xm_threshold, rcond,
    low_lev_effect, high_lev_effect,
    lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,
    lhs_diff_zm, C7_Skw_fnc,
    lhs_tp, lhs_ta_xm, lhs_pr1,
    l_implemented, solution,
    tridiag_solve_method,
    fill_holes_type,
    l_upwind_xm_ma,
    l_enable_relaxed_clipping,
    l_mono_flux_lim_thlm,
    l_mono_flux_lim_rtm,
    l_mono_flux_lim_um,
    l_mono_flux_lim_vm,
    l_mono_flux_lim_spikefix,
    stats,
    xm, xm_tol, wpxp, err_info,
    wpxp_cl_num=None,
):
    """Clip solved xm/w'x' fields and finalize implicit stats."""
    l_sample = stats.l_sample
    if solve_type == xm_wpxp_rtm:
        name_xm_ta = "rtm_ta"
        name_xm_ma = "rtm_ma"
        name_xm_pd = "rtm_pd"
        name_xm_cl = "rtm_cl"
        name_xm_matrix_condt_num = "rtm_matrix_condt_num"
        name_wpxp_ma = "wprtp_ma"
        name_wpxp_ta = "wprtp_ta"
        name_wpxp_tp = "wprtp_tp"
        name_wpxp_ac = "wprtp_ac"
        name_wpxp_pr1 = "wprtp_pr1"
        name_wpxp_pr2 = "wprtp_pr2"
        name_wpxp_dp1 = "wprtp_dp1"
        name_wpxp_pd = "wprtp_pd"
    elif solve_type == xm_wpxp_thlm:
        name_xm_ta = "thlm_ta"
        name_xm_ma = "thlm_ma"
        name_xm_pd = ""
        name_xm_cl = "thlm_cl"
        name_xm_matrix_condt_num = "thlm_matrix_condt_num"
        name_wpxp_ma = "wpthlp_ma"
        name_wpxp_ta = "wpthlp_ta"
        name_wpxp_tp = "wpthlp_tp"
        name_wpxp_ac = "wpthlp_ac"
        name_wpxp_pr1 = "wpthlp_pr1"
        name_wpxp_pr2 = "wpthlp_pr2"
        name_wpxp_dp1 = "wpthlp_dp1"
        name_wpxp_pd = ""
    elif solve_type == xm_wpxp_um:
        name_xm_ta = "um_ta"
        name_xm_ma = "um_ma"
        name_xm_pd = ""
        name_xm_cl = ""
        name_xm_matrix_condt_num = ""
        name_wpxp_ma = "upwp_ma"
        name_wpxp_ta = "upwp_ta"
        name_wpxp_tp = "upwp_tp"
        name_wpxp_ac = "upwp_ac"
        name_wpxp_pr1 = "upwp_pr1"
        name_wpxp_pr2 = "upwp_pr2"
        name_wpxp_dp1 = "upwp_dp1"
        name_wpxp_pd = ""
    elif solve_type == xm_wpxp_vm:
        name_xm_ta = "vm_ta"
        name_xm_ma = "vm_ma"
        name_xm_pd = ""
        name_xm_cl = ""
        name_xm_matrix_condt_num = ""
        name_wpxp_ma = "vpwp_ma"
        name_wpxp_ta = "vpwp_ta"
        name_wpxp_tp = "vpwp_tp"
        name_wpxp_ac = "vpwp_ac"
        name_wpxp_pr1 = "vpwp_pr1"
        name_wpxp_pr2 = "vpwp_pr2"
        name_wpxp_dp1 = "vpwp_dp1"
        name_wpxp_pd = ""
    else:
        name_xm_ta = ""
        name_xm_ma = ""
        name_xm_pd = ""
        name_xm_cl = ""
        name_xm_matrix_condt_num = ""
        name_wpxp_ma = ""
        name_wpxp_ta = ""
        name_wpxp_tp = ""
        name_wpxp_ac = ""
        name_wpxp_pr1 = ""
        name_wpxp_pr2 = ""
        name_wpxp_dp1 = ""
        name_wpxp_pd = ""

    xm_old = jnp.asarray(xm, dtype=jnp.float64)

    sol = jnp.asarray(solution, dtype=jnp.float64)
    xm = sol[:, 1:2 * nzt + 1:2]
    wpxp = sol[:, 0:2 * nzm - 1:2]

    if l_sample:
        C7_Skw_fnc_zeros = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        C7_Skw_fnc_plus_one = jnp.asarray(C7_Skw_fnc) + one
        wpxp_ac = wpxp_terms_ac_pr2_lhs(
            nzm, nzt, ngrdcol, gr, C7_Skw_fnc_zeros, wm_zt, gr.invrs_dzm,
        )
        wpxp_pr2 = wpxp_terms_ac_pr2_lhs(
            nzm, nzt, ngrdcol, gr, C7_Skw_fnc_plus_one, wm_zt, gr.invrs_dzm,
        )

        if (
            name_xm_matrix_condt_num.strip()
            and stats.var_on_stats_list(name_xm_matrix_condt_num)
        ):
            stats = stats.update(name_xm_matrix_condt_num, one / jnp.asarray(rcond))

        if not l_implemented:
            lhs_ma_zt = jnp.asarray(lhs_ma_zt, dtype=jnp.float64)
            k_zt = jnp.arange(nzt)
            km1_zt = jnp.clip(k_zt - gr.grid_dir_indx, 0, nzt - 1)
            kp1_zt = jnp.clip(k_zt + gr.grid_dir_indx, 0, nzt - 1)
            tmp_zt_stats = (
                -lhs_ma_zt[1 + gr.grid_dir_indx, :, :] * jnp.take(xm, km1_zt, axis=1)
                - lhs_ma_zt[1, :, :] * xm
                - lhs_ma_zt[1 - gr.grid_dir_indx, :, :] * jnp.take(xm, kp1_zt, axis=1)
            )
            if name_xm_ma.strip():
                stats = stats.update(name_xm_ma, tmp_zt_stats)

        lhs_ta_xm = jnp.asarray(lhs_ta_xm, dtype=jnp.float64)
        tmp_zt_stats = (
            -lhs_ta_xm[1, :, :] * wpxp[:, :nzt]
            - lhs_ta_xm[0, :, :] * wpxp[:, 1:nzt + 1]
        )
        if name_xm_ta.strip():
            stats = stats.update(name_xm_ta, tmp_zt_stats)

        k_int = slice(1, nzm - 1)
        k_minus = slice(1 - gr.grid_dir_indx, nzm - 1 - gr.grid_dir_indx)
        k_plus = slice(1 + gr.grid_dir_indx, nzm - 1 + gr.grid_dir_indx)

        tmp_zm_stats = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        lhs_ma_zm = jnp.asarray(lhs_ma_zm, dtype=jnp.float64)
        tmp_zm_stats = tmp_zm_stats.at[:, 1:nzm - 1].set(
            -lhs_ma_zm[1 + gr.grid_dir_indx, :, k_int] * wpxp[:, k_minus]
            - lhs_ma_zm[1, :, k_int] * wpxp[:, k_int]
            - lhs_ma_zm[1 - gr.grid_dir_indx, :, k_int] * wpxp[:, k_plus]
        )
        if name_wpxp_ma.strip():
            stats = stats.update(name_wpxp_ma, tmp_zm_stats)

        tmp_zm_stats = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        lhs_ta_wpxp = jnp.asarray(lhs_ta_wpxp, dtype=jnp.float64)
        tmp_zm_stats = tmp_zm_stats.at[:, 1:nzm - 1].set(
            -gamma_over_implicit_ts * lhs_ta_wpxp[1 + gr.grid_dir_indx, :, k_int] * wpxp[:, k_minus]
            - gamma_over_implicit_ts * lhs_ta_wpxp[1, :, k_int] * wpxp[:, k_int]
            - gamma_over_implicit_ts * lhs_ta_wpxp[1 - gr.grid_dir_indx, :, k_int] * wpxp[:, k_plus]
        )
        if name_wpxp_ta.strip():
            stats = stats.finalize_budget(name_wpxp_ta, tmp_zm_stats)

        tmp_zm_stats = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        lhs_tp = jnp.asarray(lhs_tp, dtype=jnp.float64)
        tmp_zm_stats = tmp_zm_stats.at[:, 1:nzm - 1].set(
            -lhs_tp[1, :, k_int] * xm[:, :nzt - 1]
            - lhs_tp[0, :, k_int] * xm[:, 1:nzt]
        )
        if name_wpxp_tp.strip():
            stats = stats.update(name_wpxp_tp, tmp_zm_stats)

        tmp_zm_stats = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        wpxp_ac = jnp.asarray(wpxp_ac)
        tmp_zm_stats = tmp_zm_stats.at[:, 1:nzm - 1].set(
            -wpxp_ac[:, k_int] * wpxp[:, k_int]
        )
        if name_wpxp_ac.strip():
            stats = stats.update(name_wpxp_ac, tmp_zm_stats)

        tmp_zm_stats = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        lhs_pr1 = jnp.asarray(lhs_pr1, dtype=jnp.float64)
        tmp_zm_stats = tmp_zm_stats.at[:, 1:nzm - 1].set(
            -gamma_over_implicit_ts * lhs_pr1[:, k_int] * wpxp[:, k_int]
        )
        if name_wpxp_pr1.strip():
            stats = stats.finalize_budget(name_wpxp_pr1, tmp_zm_stats)

        tmp_zm_stats = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        wpxp_pr2 = jnp.asarray(wpxp_pr2)
        tmp_zm_stats = tmp_zm_stats.at[:, 1:nzm - 1].set(
            -wpxp_pr2[:, k_int] * wpxp[:, k_int]
        )
        if name_wpxp_pr2.strip():
            stats = stats.update(name_wpxp_pr2, tmp_zm_stats)

        tmp_zm_stats = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        lhs_diff_zm = jnp.asarray(lhs_diff_zm, dtype=jnp.float64)
        tmp_zm_stats = tmp_zm_stats.at[:, 1:nzm - 1].set(
            -lhs_diff_zm[1 + gr.grid_dir_indx, :, k_int] * wpxp[:, k_minus]
            - lhs_diff_zm[1, :, k_int] * wpxp[:, k_int]
            - lhs_diff_zm[1 - gr.grid_dir_indx, :, k_int] * wpxp[:, k_plus]
        )
        if name_wpxp_dp1.strip():
            stats = stats.update(name_wpxp_dp1, tmp_zm_stats)

    if (
        (l_mono_flux_lim_thlm and solve_type == xm_wpxp_thlm)
        or (l_mono_flux_lim_rtm and solve_type == xm_wpxp_rtm)
        or (l_mono_flux_lim_um and solve_type == xm_wpxp_um)
        or (l_mono_flux_lim_vm and solve_type == xm_wpxp_vm)
    ):
        xm, wpxp, err_info, stats = monotonic_turbulent_flux_limit(
            nzm, nzt, ngrdcol, gr, solve_type, dt,
            xm_old, xp2, wm_zt, xm_forcing, rho_ds_zm, rho_ds_zt,
            invrs_rho_ds_zm, invrs_rho_ds_zt, xp2_threshold, xm_tol,
            l_implemented, low_lev_effect, high_lev_effect, tridiag_solve_method,
            l_upwind_xm_ma, l_mono_flux_lim_spikefix, stats, xm, wpxp, err_info,
        )

    if solve_type == xm_wpxp_rtm and l_pos_def:
        xm_pd = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        wpxp_pd = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        l_do_pos_def = jnp.any(xm < zero) & ~jnp.any(xm_old < zero)

        def pos_definite_branch(fields):
            xm_in, wpxp_in = fields
            return pos_definite_adj(
                nzm, nzt, ngrdcol, gr, dt, xm_in, wpxp_in, xm_old,
            )

        def no_pos_definite_branch(fields):
            xm_in, wpxp_in = fields
            return xm_in, wpxp_in, xm_pd, wpxp_pd

        xm, wpxp, xm_pd, wpxp_pd = jax.lax.cond(
            l_do_pos_def,
            pos_definite_branch,
            no_pos_definite_branch,
            (xm, wpxp),
        )

        if l_sample:
            if name_wpxp_pd.strip():
                stats = stats.update(name_wpxp_pd, wpxp_pd)
            if name_xm_pd.strip():
                stats = stats.update(name_xm_pd, xm_pd)

    if l_sample and name_xm_cl.strip():
        stats = stats.begin_budget(name_xm_cl, xm / dt)

    if fill_holes_type != 0 and solve_type not in (xm_wpxp_um, xm_wpxp_vm):
        xm = fill_holes_vertical(
            nzt, ngrdcol, xm_threshold,
            gr.k_lb_zt, gr.k_ub_zt,
            gr.dzt, rho_ds_zt, gr.grid_dir_indx,
            fill_holes_type, xm,
        )

    if l_sample and name_xm_cl.strip():
        stats = stats.finalize_budget(name_xm_cl, xm / dt)

    if l_enable_relaxed_clipping:
        if solve_type == xm_wpxp_thlm:
            xp2_relaxed = jnp.maximum(0.01, jnp.asarray(xp2))
        else:
            xp2_relaxed = jnp.maximum(1.0e-7, jnp.asarray(xp2))
    else:
        xp2_relaxed = jnp.asarray(xp2, dtype=jnp.float64)

    if solve_type == xm_wpxp_rtm:
        solve_type_cl = clip_wprtp
        wpxp_cl_max = wprtp_cl_max
        clip_name = "wprtp_cl"
    elif solve_type == xm_wpxp_thlm:
        solve_type_cl = clip_wpthlp
        wpxp_cl_max = wpthlp_cl_max
        clip_name = "wpthlp_cl"
    elif solve_type == xm_wpxp_um:
        solve_type_cl = clip_upwp
        wpxp_cl_max = upwp_cl_max
        clip_name = "upwp_cl"
    elif solve_type == xm_wpxp_vm:
        solve_type_cl = clip_vpwp
        wpxp_cl_max = vpwp_cl_max
        clip_name = "vpwp_cl"
    else:
        solve_type_cl = clip_wpsclrp
        wpxp_cl_max = 0
        clip_name = ""

    if wpxp_cl_num is not None and l_sample and clip_name.strip():
        if wpxp_cl_num == 0:
            stats = stats.begin_budget(clip_name, wpxp / dt)
        else:
            stats = stats.update_budget(clip_name, -wpxp / dt)

    if solve_type not in (xm_wpxp_um, xm_wpxp_vm):
        if wpxp_cl_num is not None:
            wpxp_cl_num += 1
        wpxp, wpxp_chnge = clip_covar(
            nzm, ngrdcol, solve_type_cl, wp2, xp2_relaxed, wpxp,
        )
        if wpxp_cl_num is not None and l_sample and clip_name.strip():
            if wpxp_cl_num == wpxp_cl_max:
                stats = stats.finalize_budget(clip_name, wpxp / dt)
            else:
                stats = stats.update_budget(clip_name, wpxp / dt)
    else:
        if wpxp_cl_num is not None:
            wpxp_cl_num += 1
        wpxp, wpxp_chnge = clip_covar(
            nzm, ngrdcol, solve_type_cl, wp2, xp2, wpxp,
        )
        if wpxp_cl_num is not None and l_sample and clip_name.strip():
            if wpxp_cl_num == wpxp_cl_max:
                stats = stats.finalize_budget(clip_name, wpxp / dt)
            else:
                stats = stats.update_budget(clip_name, wpxp / dt)

    if l_clip_turb_adv:
        xm, stats = xm_correction_wpxp_cl(
            nzm, nzt, ngrdcol, solve_type, dt, wpxp_chnge,
            gr.invrs_dzt, stats, xm,
        )

    return xm, wpxp, err_info, wpxp_cl_num, stats


def xm_term_ta_lhs(nzm, nzt, ngrdcol, gr, rho_ds_zm, invrs_rho_ds_zt):
    """Turbulent advection of xm:  implicit portion of the code.

    The d(xm)/dt equation contains a turbulent advection term:

    - (1/rho_ds) * d( rho_ds * w'x' )/dz.

    This term is solved for completely implicitly, such that:

    - (1/rho_ds) * d( rho_ds * w'x'(t+1) )/dz.

    Note:  When the term is brought over to the left-hand side, the sign
           is reversed and the leading "-" in front of the term is changed
           to a "+".

    invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) )
    """
    lhs_ta_xm = jnp.zeros((ndiags2, ngrdcol, nzt), dtype=jnp.float64)
    # Momentum superdiagonal [ x wpxp(k+1,<t+1>) ]
    lhs_ta_xm = lhs_ta_xm.at[0, :, :].set(
        jnp.asarray(invrs_rho_ds_zt) * jnp.asarray(gr.invrs_dzt) * jnp.asarray(rho_ds_zm[:, 1:nzm])
    )
    # Momentum subdiagonal [ x wpxp(k,<t+1>) ]
    lhs_ta_xm = lhs_ta_xm.at[1, :, :].set(
        -jnp.asarray(invrs_rho_ds_zt) * jnp.asarray(gr.invrs_dzt) * jnp.asarray(rho_ds_zm[:, :nzt])
    )
    return lhs_ta_xm


def wpxp_term_tp_lhs(nzm, ngrdcol, gr, wp2):
    """Turbulent production of w'x':  implicit portion of the code.

    The d(w'x')/dt equation contains a turbulent production term:

    - w'^2 d(xm)/dz.

    This term is solved for completely implicitly, such that:

    - w'^2 d(xm(t+1))/dz.

    invrs_dzm(k) = 1 / ( zt(k) - zt(k-1) )
    """
    lhs_tp = jnp.zeros((ndiags2, ngrdcol, nzm), dtype=jnp.float64)
    # Thermodynamic superdiagonal [ x xm(k,<t+1>) ]
    lhs_tp = lhs_tp.at[0, :, 1:nzm - 1].set(
        jnp.asarray(wp2[:, 1:nzm - 1]) * jnp.asarray(gr.invrs_dzm[:, 1:nzm - 1])
    )
    # Thermodynamic subdiagonal [ x xm(k-1,<t+1>) ]
    lhs_tp = lhs_tp.at[1, :, 1:nzm - 1].set(
        -jnp.asarray(wp2[:, 1:nzm - 1]) * jnp.asarray(gr.invrs_dzm[:, 1:nzm - 1])
    )
    return lhs_tp


def wpxp_terms_ac_pr2_lhs(nzm, nzt, ngrdcol, gr, C7_Skw_fnc, wm_zt, invrs_dzm):
    """Accumulation of w'x' and w'x' pressure term 2:  implicit portion of the
    code.

    The d(w'x')/dt equation contains an accumulation term:

    - w'x' dw/dz;

    and pressure term 2:

    - C_7 ( -w'x' dw/dz ).

    Note:  When the term is brought over to the left-hand side, the sign
           is reversed and the leading "-" in front of the term is changed
           to a "+".
    """
    lhs_ac_pr2 = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    # Momentum main diagonal: [ x wpxp(k,<t+1>) ]
    lhs_ac_pr2 = lhs_ac_pr2.at[:, 1:nzm - 1].set(
        (one - jnp.asarray(C7_Skw_fnc[:, 1:nzm - 1]))
        * jnp.asarray(invrs_dzm[:, 1:nzm - 1])
        * (jnp.asarray(wm_zt[:, 1:nzt]) - jnp.asarray(wm_zt[:, :nzt - 1]))
    )
    return lhs_ac_pr2


def wpxp_term_pr1_lhs(
    nzm, ngrdcol, gr, C6rt_Skw_fnc,
    C6thl_Skw_fnc, C7_Skw_fnc,
    invrs_tau_C6_zm, l_scalar_calc,
):
    """Pressure term 1 for w'x':  implicit portion of the code.

    The d(w'x')/dt equation contains pressure term 1:

    - ( C_6 / tau_m ) w'x'.

    This term is solved for completely implicitly, such that:

    - ( C_6 / tau_m ) w'x'(t+1).

    Note:  When the term is brought over to the left-hand side, the sign
           is reversed and the leading "-" in front of the term is changed
           to a "+".
    """
    lhs_pr1_wprtp = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    lhs_pr1_wpthlp = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    lhs_pr1_wpsclrp = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)

    # Momentum main diagonals: [ x wpxp(k,<t+1>) ]
    lhs_pr1_wprtp = lhs_pr1_wprtp.at[:, 1:nzm - 1].set(
        jnp.asarray(C6rt_Skw_fnc[:, 1:nzm - 1]) * jnp.asarray(invrs_tau_C6_zm[:, 1:nzm - 1])
    )
    # Momentum main diagonals: [ x wpxp(k,<t+1>) ]
    lhs_pr1_wpthlp = lhs_pr1_wpthlp.at[:, 1:nzm - 1].set(
        jnp.asarray(C6thl_Skw_fnc[:, 1:nzm - 1]) * jnp.asarray(invrs_tau_C6_zm[:, 1:nzm - 1])
    )

    if l_scalar_calc:
        # Momentum main diagonals: [ x wpxp(k,<t+1>) ]
        lhs_pr1_wpsclrp = lhs_pr1_wpsclrp.at[:, 1:nzm - 1].set(
            jnp.asarray(C7_Skw_fnc[:, 1:nzm - 1]) * jnp.asarray(invrs_tau_C6_zm[:, 1:nzm - 1])
        )

    return lhs_pr1_wprtp, lhs_pr1_wpthlp, lhs_pr1_wpsclrp


def wpxp_terms_bp_pr3_rhs(nzm, ngrdcol, gr, C7_Skw_fnc, thv_ds_zm, xpthvp):
    """Buoyancy production of w'x' and w'x' pressure term 3:  explicit portion of
    the code.

    The d(w'x')/dt equation contains a buoyancy production term:

    + (g/thv_ds) x'th_v';

    and pressure term 3:

    - C_7 ( (g/thv_ds) x'th_v' ).
    """
    rhs_bp_pr3 = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    # Calculate term at all interior grid levels.
    rhs_bp_pr3 = rhs_bp_pr3.at[:, 1:nzm - 1].set(
        (grav / jnp.asarray(thv_ds_zm[:, 1:nzm - 1]))
        * (one - jnp.asarray(C7_Skw_fnc[:, 1:nzm - 1]))
        * jnp.asarray(xpthvp[:, 1:nzm - 1])
    )
    return rhs_bp_pr3


def xm_correction_wpxp_cl(
    nzm, nzt, ngrdcol, solve_type, dt,
    wpxp_chnge, invrs_dzt,
    stats,
    xm,
):
    """Corrects the value of xm if w'x' needed to be clipped, for xm is partially
    based on the derivative of w'x' with respect to altitude.

    The time-tendency equation for xm is:

    d(xm)/dt = -w d(xm)/dz - d(w'x')/dz + d(xm)/dt|_ls;

    where d(xm)/dt|_ls is the rate of change of xm over time due to radiation,
    microphysics, and/or any other large-scale forcing(s).

    Note:  The results of this xm adjustment are highly dependent on the
           order in which clipping and positive-definite adjustments are applied.
    """
    l_sample = stats.l_sample
    wpxp_chnge = jnp.asarray(wpxp_chnge)
    xm = jnp.asarray(xm, dtype=jnp.float64)
    l_clipping_needed = jnp.any(jnp.abs(wpxp_chnge) > eps, axis=1)

    if solve_type == xm_wpxp_rtm:
        name_xm_tacl = "rtm_tacl"
    elif solve_type == xm_wpxp_thlm:
        name_xm_tacl = "thlm_tacl"
    else:
        name_xm_tacl = ""

    invrs_dzt = jnp.asarray(invrs_dzt)
    # The adjustment to xm due to turbulent advection term clipping
    xm_tndcy_wpxp_cl = -invrs_dzt * (wpxp_chnge[:, 1:] - wpxp_chnge[:, :-1])
    xm_tndcy_wpxp_cl = jnp.where(l_clipping_needed[:, None], xm_tndcy_wpxp_cl, zero)
    xm = jnp.where(l_clipping_needed[:, None], xm + xm_tndcy_wpxp_cl * dt, xm)

    if l_sample and name_xm_tacl.strip():
        stats = stats.update(name_xm_tacl, xm_tndcy_wpxp_cl)

    return xm, stats


def damp_coefficient(
    nzm, ngrdcol, gr, coefficient, Cx_Skw_fnc,
    max_coeff_value, altitude_threshold,
    threshold, Lscale_zm,
):
    """Damps a given coefficient linearly based on the value of Lscale.

    For additional information see CLUBB ticket #431.
    """
    coefficient = jnp.asarray(coefficient)
    max_coeff_value = jnp.asarray(max_coeff_value)
    altitude_threshold = jnp.asarray(altitude_threshold)
    threshold = jnp.asarray(threshold)
    Cx_Skw_fnc = jnp.asarray(Cx_Skw_fnc)
    Lscale_zm = jnp.asarray(Lscale_zm)

    mask = (Lscale_zm < threshold[:, None]) & (jnp.asarray(gr.zm) > altitude_threshold[:, None])
    threshold_safe = jnp.where(jnp.abs(threshold) > zero, threshold, one)
    damped = (
        max_coeff_value[:, None]
        + ((coefficient - max_coeff_value) / threshold_safe)[:, None] * Lscale_zm
    )
    return jnp.where(mask, damped, Cx_Skw_fnc)


def diagnose_upxp(
    nzm, nzt, ngrdcol, gr, ypwp, xm, wpxp, ym,
    C6x_Skw_fnc, tau_C6_zm, C7_Skw_fnc,
):
    """Diagnose turbulent horizontal flux of a conserved scalar.

    References:
      None
    """
    ddzt_xm = ddzt(nzm, nzt, ngrdcol, gr, xm)
    ddzt_ym = ddzt(nzm, nzt, ngrdcol, gr, ym)
    ypxp = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    # The value of ypxp is irrelevant to the calculations at the upper and
    # lower boundaries
    ypxp_interior = (
        (jnp.asarray(tau_C6_zm)[:, 1:nzm - 1] / jnp.asarray(C6x_Skw_fnc)[:, 1:nzm - 1])
        * (
            -jnp.asarray(ypwp)[:, 1:nzm - 1] * ddzt_xm[:, 1:nzm - 1]
            - (one - jnp.asarray(C7_Skw_fnc)[:, 1:nzm - 1])
            * (jnp.asarray(wpxp)[:, 1:nzm - 1] * ddzt_ym[:, 1:nzm - 1])
        )
    )
    ypxp = ypxp.at[:, 1:nzm - 1].set(ypxp_interior)
    ypxp = ypxp.at[:, gr.k_lb_zm].set(zero)
    ypxp = ypxp.at[:, gr.k_ub_zm].set(zero)
    return ypxp


def error_prints_xm_wpxp(*args, **kwargs):
    """Prints values of model fields when fatal errors (LU decomp.) occur.

    Porting deviation:
      The JAX implementation prints a concise placeholder because the full
      Fortran diagnostic dump relies on host-side I/O and extensive mutable
      argument state that is intentionally not carried through the jitted path.
    """
    del args, kwargs
    print("Error in advance_xm_wpxp")


# JAX adaptation: this target-only helper compresses repeated multiple-LHS solve
# blocks. Keep it small and limited to the repeated Fortran call pattern; expand
# it back into explicit call sites if it starts hiding field-specific behavior.
def _solve_one_xm_wpxp(
    nzm, nzt, ngrdcol, gr, solve_type, dt, l_iter, wm_zt, wp2, xp2,
    xm_forcing, wpxp_forcing, C7_Skw_fnc, xpthvp, rhs_ta,
    thv_ds_zm, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,
    xp2_threshold, xm_threshold, xm_tol, low_lev_effect, high_lev_effect,
    lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,
    lhs_tp, lhs_ta_xm, lhs_ac_pr2, lhs_pr1, l_implemented,
    penta_solve_method, tridiag_solve_method, fill_holes_type,
    l_diffuse_rtm_and_thlm, l_upwind_xm_ma,
    l_enable_relaxed_clipping, l_mono_flux_lim_thlm, l_mono_flux_lim_rtm,
    l_mono_flux_lim_um, l_mono_flux_lim_vm, l_mono_flux_lim_spikefix,
    stats, xm, wpxp, err_info, wpxp_cl_num=None,
):
    l_sample = stats.l_sample
    lhs = xm_wpxp_lhs(
        nzm, nzt, ngrdcol, l_iter, dt,
        l_implemented, lhs_diff_zm, lhs_diff_zt,
        lhs_ma_zm, lhs_ma_zt, lhs_ta_wpxp, lhs_ta_xm,
        lhs_tp, lhs_pr1, lhs_ac_pr2,
        l_diffuse_rtm_and_thlm,
    )
    rhs_2d, stats = xm_wpxp_rhs(
        nzm, nzt, ngrdcol, gr, solve_type, l_iter, dt,
        xm, wpxp, xm_forcing, wpxp_forcing, C7_Skw_fnc,
        xpthvp, rhs_ta, thv_ds_zm, lhs_pr1, lhs_ta_wpxp,
        stats,
    )
    rhs = rhs_2d[:, :, None]
    old_solution = jnp.zeros((ngrdcol, 2 * nzm - 1, 1), dtype=jnp.float64)
    if penta_solve_method == penta_bicgstab:
        old_solution = old_solution.at[:, 1:2 * (nzm - 1) + 1:2, 0].set(
            jnp.asarray(xm)
        )
        old_solution = old_solution.at[:, 0:2 * nzm - 1:2, 0].set(
            jnp.asarray(wpxp)
        )

    if solve_type == xm_wpxp_rtm:
        cond_name = "rtm_matrix_condt_num"
    elif solve_type == xm_wpxp_thlm:
        cond_name = "thlm_matrix_condt_num"
    else:
        cond_name = ""

    solution, rcond, err_info = xm_wpxp_solve(
        nzm, ngrdcol, 1, gr, penta_solve_method, l_implemented,
        old_solution, lhs, rhs, err_info,
        use_rcond=bool(
            l_sample
            and cond_name.strip()
            and stats.var_on_stats_list(cond_name)
        ),
    )

    return xm_wpxp_clipping_and_stats(
        nzm, nzt, ngrdcol, gr, solve_type, dt, wp2, xp2, wm_zt,
        xm_forcing, rho_ds_zm, rho_ds_zt,
        invrs_rho_ds_zm, invrs_rho_ds_zt,
        xp2_threshold, xm_threshold, rcond,
        low_lev_effect, high_lev_effect,
        lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,
        lhs_diff_zm, C7_Skw_fnc,
        lhs_tp, lhs_ta_xm, lhs_pr1,
        l_implemented, solution[:, :, 0],
        tridiag_solve_method, fill_holes_type,
        l_upwind_xm_ma, l_enable_relaxed_clipping,
        l_mono_flux_lim_thlm, l_mono_flux_lim_rtm,
        l_mono_flux_lim_um, l_mono_flux_lim_vm,
        l_mono_flux_lim_spikefix,
        stats,
        xm, xm_tol, wpxp, err_info,
        wpxp_cl_num=wpxp_cl_num,
    )
