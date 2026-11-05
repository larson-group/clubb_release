"""JAX-side entry point for ``src/CLUBB_core/advance_xp2_xpyp_module.F90``.

Description:
  Contains the subroutine advance_xp2_xpyp and ancillary functions. Prognose
  scalar variances, scalar covariances, and horizontal turbulence components.

References:
  https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:xpyp_eqns
  https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:up2_vp2_eqns

  Eqn. 13, 14, 15 on p. 3545 of
  ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    Method and Model Description'' Golaz, et al. (2002)
    JAS, Vol. 59, pp. 3540--3551.

See also:
  ``Equations for CLUBB'', Section 4:
  /Steady-state solution for the variances/

Porting deviations:
- Sponge damping is intentionally omitted because the JAX case initialization
  rejects sponge-enabled cases before this routine is called.
- xp2_xpyp_single_lhs_valid guards the single-LHS branch. The unsupported
  multiple-LHS path prints a runtime diagnostic and marks err_info fatal.
  It should be revisited when that branch is expressed as JAX control flow.
- Scalar IDs remain Fortran-style 1-based values because they are logical
  category names rather than Python array indices.
"""

from functools import partial

import jax
import sys

import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.advance_helper_module import calc_wp3_on_wp2
from clubb_jax.src.CLUBB_core.clip_explicit import clip_covar, clip_variance
from clubb_jax.src.CLUBB_core.diffusion import diffusion_zm_lhs
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from clubb_jax.src.CLUBB_core.fill_holes import fill_holes_vertical
from clubb_jax.src.CLUBB_core.grid_class import zm2zt, zt2zm
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.mean_adv import term_ma_zm_lhs
from clubb_jax.src.CLUBB_core.matrix_solver_wrapper import tridiag_solve
from clubb_jax.src.CLUBB_core.turbulent_adv_pdf import (
    xpyp_term_ta_pdf_lhs,
    xpyp_term_ta_pdf_lhs_godunov,
    xpyp_term_ta_pdf_rhs,
    xpyp_term_ta_pdf_rhs_godunov,
)
from clubb_jax.src.CLUBB_core.constants_clubb import (
    Cp,
    Lv,
    cloud_frac_min,
    eps,
    gamma_over_implicit_ts,
    grav,
    max_mag_correlation_flux,
    one,
    one_half,
    one_third,
    rt_tol,
    thl_tol,
    two,
    two_thirds,
    w_tol_sqd,
    zero,
    zero_threshold,
)
from clubb_jax.src.CLUBB_core.parameter_indices import (
    ibeta,
    ic_K2,
    ic_K9,
    iC2rt,
    iC2rtthl,
    iC2thl,
    iC4,
    iC14,
    iC_uu_buoy,
    iC_uu_shr,
    irtp2_clip_coef,
)
from clubb_jax.src.CLUBB_core.model_flags import (
    iiPDF_ADG1,
    iiPDF_new,
    iiPDF_new_hybrid,
    l_explicit_turbulent_adv_xpyp,
    l_force_descending_solves,
)
from clubb_jax.src.CLUBB_core import (
    ErrInfo,
    Grid,
    NuVertResDep,
    SclrIdx,
    implicit_coefs_terms,
    pdf_parameter,
)


# Private named constants to avoid string comparisons
xp2_xpyp_rtp2 = 1  # Named constant for rtp2 solves
xp2_xpyp_thlp2 = 2  # Named constant for thlp2 solves
xp2_xpyp_rtpthlp = 3  # Named constant for rtpthlp solves
xp2_xpyp_up2_vp2 = 4  # Named constant for up2_vp2 solves
xp2_xpyp_up2 = 5  # Named constant for up2 solves
xp2_xpyp_vp2 = 6  # Named constant for vp2 solves
xp2_xpyp_scalars = 7  # Named constant for scalar solves
xp2_xpyp_sclrp2 = 8  # Named constant for sclrp2 solves
xp2_xpyp_sclrprtp = 9  # Named constant for sclrprtp solves
xp2_xpyp_sclrpthlp = 10  # Named constant for sclrpthlp solves
xp2_xpyp_single_lhs = 11  # Named constant for single lhs solve

clip_sclrp2 = 14
clip_sclrprtp = 15
clip_sclrpthlp = 16

ndiags3 = 3
min_cloud_frac_mult = 0.10


@partial(
    jax.jit,
    static_argnames=("iiPDF_type", "l_explicit_turbulent_adv_xpyp"),
)
def xp2_xpyp_single_lhs_valid(clubb_params, iiPDF_type, l_explicit_turbulent_adv_xpyp):
    """Return whether the single-LHS solve is valid for the supplied parameters.

    TODO(JAX port): fold this target-only helper back into the Fortran-order
    control flow when both solve branches are JAX-native.
    """
    static_allowed = (
        bool(l_explicit_turbulent_adv_xpyp)
        or (not bool(l_explicit_turbulent_adv_xpyp) and int(iiPDF_type) == iiPDF_ADG1)
    )
    c2rt = clubb_params[:, iC2rt]
    c2thl = clubb_params[:, iC2thl]
    c2rtthl = clubb_params[:, iC2rtthl]
    c2_equivalent = jnp.all(
        (
            jnp.abs(c2rt - c2thl)
            <= jnp.abs(c2rt + c2thl) * one_half * eps
        )
        & (
            jnp.abs(c2rt - c2rtthl)
            <= jnp.abs(c2rt + c2rtthl) * one_half * eps
        )
    )
    l_single_lhs_valid = static_allowed & c2_equivalent
    return l_single_lhs_valid


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "nzt",
        "ngrdcol",
        "sclr_dim",
        "iiPDF_type",
        "tridiag_solve_method",
        "fill_holes_type",
        "l_ho_nontrad_coriolis",
        "l_min_xp2_from_corr_wx",
        "l_C2_cloud_frac",
        "l_upwind_xpyp_ta",
        "l_godunov_upwind_xpyp_ta",
        "l_lmm_stepping",
        "l_implemented",
    ),
)
def advance_xp2_xpyp(
    nzm: int, nzt: int, ngrdcol: int, sclr_dim: int, sclr_tol, gr: Grid, sclr_idx: SclrIdx,
    invrs_tau_xp2_zm, invrs_tau_C4_zm,
    invrs_tau_C14_zm, wm_zm,
    rtm, wprtp, thlm, wpthlp, wpthvp, um, vm,
    wp2, wp3, upwp, vpwp,
    sigma_sqd_w, wprtp2, wpthlp2,
    wprtpthlp, Kh_zt, rtp2_forcing,
    thlp2_forcing, rtpthlp_forcing,
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,
    thv_ds_zm, cloud_frac,
    pdf_implicit_coefs_terms: implicit_coefs_terms,
    dt: float, fcor_y,
    sclrm, wpsclrp,
    wpsclrp2, wpsclrprtp, wpsclrpthlp,
    lhs_splat_wp2,
    clubb_params, nu_vert_res_dep: NuVertResDep,
    iiPDF_type: int,
    tridiag_solve_method: int,
    fill_holes_type: int,
    l_ho_nontrad_coriolis: bool,
    l_min_xp2_from_corr_wx: bool,
    l_C2_cloud_frac: bool,
    l_upwind_xpyp_ta: bool,
    l_godunov_upwind_xpyp_ta: bool,
    l_lmm_stepping: bool,
    l_implemented: bool,
    stats: JaxStats,
    rtp2, thlp2, rtpthlp, up2, vp2,
    sclrp2, sclrprtp, sclrpthlp, err_info: ErrInfo,
):
    """Prognose scalar variances, scalar covariances, and horizontal turbulence components.

    References:
    https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:xpyp_eqns

    https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:up2_vp2_eqns

      Eqn. 13, 14, 15  on p. 3545 of
      ``A PDF-Based Model for Boundary Layer Clouds. Part I:
        Method and Model Description'' Golaz, et al. (2002)
      JAS, Vol. 59, pp. 3540--3551.

    Returns:
        (rtp2, thlp2, rtpthlp, up2, vp2, sclrp2, sclrprtp, sclrpthlp, err_info,
          stats)

    Fortran declarations (advance_xp2_xpyp in src/CLUBB_core/advance_xp2_xpyp_module.F90):

        !------------------------------ Input Variables ------------------------------
        integer, intent(in) :: &
          nzt,      & ! Number of thermodynamic vertical levels
          nzm,      & ! Number of momentum vertical levels
          ngrdcol,  & ! Number of grid columns
          sclr_dim    ! Number of passive scalars

        real( kind = core_rknd ), intent(in), dimension(sclr_dim) :: &
          sclr_tol          ! Threshold(s) on the passive scalars  [units vary]

        type (grid), intent(in) :: &
          gr

        type (sclr_idx_type), intent(in) :: &
          sclr_idx

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) :: &
          invrs_tau_xp2_zm, & ! Inverse time-scale for xp2 on momentum levels [1/s]
          invrs_tau_C4_zm,  & ! Inverse time-scale for C4 terms on mom. levels [1/s]
          invrs_tau_C14_zm, & ! Inverse time-scale for C14 terms on mom. levels [1/s]
          wm_zm,            & ! w-wind component on momentum levels   [m/s]
          wprtp,            & ! <w'r_t'> (momentum levels)            [(m/s)(kg/kg)]
          wpthlp,           & ! <w'th_l'> (momentum levels)           [(m K)/s]
          wpthvp,           & ! <w'th_v'> (momentum levels)           [(m K)/s]
          wp2,              & ! <w'^2> (momentum levels)              [m^2/s^2]
          upwp,             & ! <u'w'> (momentum levels)              [m^2/s^2]
          vpwp,             & ! <v'w'> (momentum levels)              [m^2/s^2]
          sigma_sqd_w,      & ! sigma_sqd_w (momentum levels)         [-]
          rtp2_forcing,     & ! <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
          thlp2_forcing,    & ! <th_l'^2> forcing (momentum levels)   [K^2/s]
          rtpthlp_forcing,  & ! <r_t'th_l'> forcing (momentum levels) [(kg/kg)K/s]
          rho_ds_zm,        & ! Dry, static density on momentum levs. [kg/m^3]
          invrs_rho_ds_zm,  & ! Inv. dry, static density @ mom. levs. [m^3/kg]
          thv_ds_zm           ! Dry, base-state theta_v on mom. levs. [K]

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) :: &
          rtm,              & ! Total water mixing ratio (t-levs)     [kg/kg]
          thlm,             & ! Liquid potential temp. (t-levs)       [K]
          um,               & ! u wind (thermodynamic levels)         [m/s]
          vm,               & ! v wind (thermodynamic levels)         [m/s]
          wp3,              & ! <w'^3> (thermodynamic levels)         [m^3/s^3]
          wprtp2,           & ! <w'r_t'^2> (thermodynamic levels)     [m/s (kg/kg)^2]
          wpthlp2,          & ! <w'th_l'^2> (thermodynamic levels)    [m/s K^2]
          wprtpthlp,        & ! <w'r_t'th_l'> (thermodynamic levels)  [m/s (kg/kg) K]
          Kh_zt,            & ! Eddy diffusivity on thermo. levels    [m^2/s]
          rho_ds_zt,        & ! Dry, static density on thermo. levels [kg/m^3]
          cloud_frac          ! Cloud fraction (thermodynamic levels) [-]

        type(implicit_coefs_terms), intent(in) :: &
          pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

        real( kind = core_rknd ), intent(in) :: &
          dt             ! Model timestep                                [s]

        real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
          fcor_y             ! Nontraditional Coriolis parameter         [s^-1]
                             ! Meridional planetary vorticity. Proportional to cos(latitude)

        ! Passive scalar input
        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt,sclr_dim) :: &
          sclrm,       & ! Mean value; pass. scalar (t-levs.) [{sclr units}]
          wpsclrp2,    & ! <w'sclr'^2> (thermodynamic levels) [m/s{sclr units}^2]
          wpsclrprtp,  & ! <w'sclr'r_t'> (thermo. levels)     [m/s{sclr units)kg/kg]
          wpsclrpthlp    ! <w'sclr'th_l'> (thermo. levels)    [m/s{sclr units}K]

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm,sclr_dim) :: &
          wpsclrp        ! <w'sclr'> (momentum levels)        [m/s{sclr units}]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          lhs_splat_wp2  ! LHS coefficient of wp2 splatting term  [1/s]

        real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
          clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

        type(nu_vertical_res_dep), intent(in) :: &
          nu_vert_res_dep    ! Vertical resolution dependent nu values

        integer, intent(in) :: &
          iiPDF_type,           & ! Selected option for the two-component normal (double
                                  ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                                  ! w, chi, and eta) portion of CLUBB's multivariate,
                                  ! two-component PDF.
          tridiag_solve_method, & ! Specifier for method to solve tridiagonal systems
          fill_holes_type         ! Specifier for which hole filling method to use

        logical, intent(in) :: &
          l_ho_nontrad_coriolis,  & ! Flag to implement the nontraditional Coriolis terms in the
                                    ! prognostic equations of <w'w'>, <u'w'>, and <u'u'>.
          l_min_xp2_from_corr_wx, & ! Flag to base the threshold minimum value of xp2 (rtp2 and thlp2)
                                    ! on keeping the overall correlation of w and x within the limits
                                    ! of -max_mag_correlation_flux to max_mag_correlation_flux.
          l_C2_cloud_frac,        & ! Flag to use cloud fraction to adjust the value of the turbulent
                                    ! dissipation coefficient, C2.
          l_upwind_xpyp_ta,       & ! This flag determines whether we want to use an upwind
                                    ! differencing approximation rather than a centered differencing
                                    ! for turbulent or mean advection terms. It affects rtp2, thlp2,
                                    ! up2, vp2, sclrp2, rtpthlp, sclrprtp, & sclrpthlp.
          l_godunov_upwind_xpyp_ta,  & ! This flag determines whether we want to use a Godunov-like
                                       ! upwind differencing approximation rather than a
                                       ! centered differencing for turbulent or mean advection terms.
                                       ! It affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp,
                                       ! & sclrpthlp.
          l_lmm_stepping, &            ! Apply Linear Multistep Method (LMM) Stepping
          l_implemented                ! True if CLUBB is being implemented and run in a host model

        !------------------- InOut Variables -------------------
        type(stats_type), intent(inout) :: &
          stats

        !------------------------------ Input/Output Variables ------------------------------
        ! An attribute of (inout) is also needed to import the value of the variances
        ! at the surface.  Brian.  12/18/05.
        real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzm) ::  &
          rtp2,    & ! <r_t'^2>                      [(kg/kg)^2]
          thlp2,   & ! <th_l'^2>                     [K^2]
          rtpthlp, & ! <r_t'th_l'>                   [(kg K)/kg]
          up2,     & ! <u'^2>                        [m^2/s^2]
          vp2        ! <v'^2>                        [m^2/s^2]

        ! Passive scalar output
        real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzm,sclr_dim) ::  &
          sclrp2, sclrprtp, sclrpthlp

        type(err_info_type), intent(inout) :: &
          err_info      ! err_info struct containing err_code and err_header
    """
    l_sample = stats.l_sample
    l_scalar_calc = sclr_dim > 0
    initial_rtp2, initial_thlp2, initial_rtpthlp = rtp2, thlp2, rtpthlp
    initial_up2, initial_vp2 = up2, vp2
    initial_sclrp2, initial_sclrprtp, initial_sclrpthlp = (
        sclrp2, sclrprtp, sclrpthlp,
    )
    wp3_on_wp2, wp3_on_wp2_zt = calc_wp3_on_wp2(
        nzm, nzt, ngrdcol, gr, wp2, wp3,
    )
    wp2_zt = jnp.maximum(zm2zt(nzm, nzt, ngrdcol, gr, wp2), w_tol_sqd)

    if l_C2_cloud_frac:
        cloud_frac_zm = zt2zm(nzm, nzt, ngrdcol, gr, cloud_frac)
        cloud_frac_mult = jnp.maximum(min_cloud_frac_mult, cloud_frac_zm)
        C2rt_1d = jnp.where(
            cloud_frac_zm >= cloud_frac_min,
            clubb_params[:, iC2rt][:, None] * cloud_frac_mult,
            clubb_params[:, iC2rt][:, None],
        )
        C2thl_1d = jnp.where(
            cloud_frac_zm >= cloud_frac_min,
            clubb_params[:, iC2thl][:, None] * cloud_frac_mult,
            clubb_params[:, iC2thl][:, None],
        )
        C2rtthl_1d = jnp.where(
            cloud_frac_zm >= cloud_frac_min,
            clubb_params[:, iC2rtthl][:, None] * cloud_frac_mult,
            clubb_params[:, iC2rtthl][:, None],
        )
    else:
        C2rt_1d = jnp.broadcast_to(clubb_params[:, iC2rt][:, None], (ngrdcol, nzm))
        C2thl_1d = jnp.broadcast_to(clubb_params[:, iC2thl][:, None], (ngrdcol, nzm))
        C2rtthl_1d = jnp.broadcast_to(clubb_params[:, iC2rtthl][:, None], (ngrdcol, nzm))

    C2sclr_1d = jnp.broadcast_to(clubb_params[:, iC2rt][:, None], (ngrdcol, nzm))
    C4_1d = two_thirds * jnp.broadcast_to(clubb_params[:, iC4][:, None], (ngrdcol, nzm))
    C14_1d = one_third * jnp.broadcast_to(clubb_params[:, iC14][:, None], (ngrdcol, nzm))

    Kw2 = clubb_params[:, ic_K2][:, None] * Kh_zt
    Kw9 = clubb_params[:, ic_K9][:, None] * Kh_zt
    Kw2_zm = jnp.maximum(zt2zm(nzm, nzt, ngrdcol, gr, Kw2), zero)
    Kw9_zm = jnp.maximum(zt2zm(nzm, nzt, ngrdcol, gr, Kw9), zero)

    (
        lhs_ta_wprtp2, lhs_ta_wpthlp2, lhs_ta_wprtpthlp,
        lhs_ta_wpup2, lhs_ta_wpvp2, lhs_ta_wpsclrp2,
        lhs_ta_wprtpsclrp, lhs_ta_wpthlpsclrp,
        rhs_ta_wprtp2, rhs_ta_wpthlp2, rhs_ta_wprtpthlp,
        rhs_ta_wpup2, rhs_ta_wpvp2, rhs_ta_wpsclrp2,
        rhs_ta_wprtpsclrp, rhs_ta_wpthlpsclrp,
        stats,
    ) = calc_xp2_xpyp_ta_terms(
        nzm, nzt, ngrdcol, sclr_dim, gr,
        wprtp, wprtp2, wpthlp, wpthlp2, wprtpthlp,
        rtp2, thlp2, rtpthlp, upwp, vpwp, up2, vp2, wp2,
        wp2_zt, wpsclrp, wpsclrp2, wpsclrprtp, wpsclrpthlp,
        sclrp2, sclrprtp, sclrpthlp,
        rho_ds_zt, invrs_rho_ds_zm, rho_ds_zm,
        wp3_on_wp2, wp3_on_wp2_zt, sigma_sqd_w,
        pdf_implicit_coefs_terms, l_scalar_calc,
        clubb_params[:, ibeta], iiPDF_type, l_upwind_xpyp_ta,
        l_godunov_upwind_xpyp_ta, stats,
    )

    lhs_diff = diffusion_zm_lhs(
        nzm, nzt, ngrdcol, gr, Kw2, Kw2_zm, nu_vert_res_dep.nu2,
        invrs_rho_ds_zm, rho_ds_zt,
    )
    lhs_ma = term_ma_zm_lhs(
        nzm, nzt, ngrdcol, wm_zm, gr.invrs_dzm, gr.weights_zm2zt,
    )

    l_single_lhs_valid = xp2_xpyp_single_lhs_valid(
        clubb_params, iiPDF_type, l_explicit_turbulent_adv_xpyp,
    )
    l_multiple_lhs_required = jnp.logical_not(l_single_lhs_valid)
    err_info = err_info.set_fatal(l_multiple_lhs_required)
    jax.lax.cond(
        l_multiple_lhs_required,
        lambda _: jax.debug.callback(
            lambda: print(
                "advance_xp2_xpyp requires the multiple-LHS solve, "
                "which is not supported in JAX.\n"
                "A single-LHS solve requires:\n"
                "  (l_explicit_turbulent_adv_xpyp or iiPDF_type == iiPDF_ADG1)\n"
                "  and, in every column, C2thl and C2rtthl must each match C2rt\n"
                "  within the relative epsilon tolerance.",
                file=sys.stderr,
            ),
            ordered=True,
        ),
        lambda _: None,
        operand=None,
    )

    # TODO(JAX port): future lax.cond branch point:
    #   true branch: solve_xp2_xpyp_with_single_lhs
    #   false branch: solve_xp2_xpyp_with_multiple_lhs
    # Until the solve branches are free of Python stats/debug side effects, the
    # multiple-LHS branch is treated as an error and the incoming state is kept.
    rtp2, thlp2, rtpthlp, sclrp2, sclrprtp, sclrpthlp, err_info, stats = (
        solve_xp2_xpyp_with_single_lhs(
            nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, sclr_idx,
            C2rt_1d, invrs_tau_xp2_zm, rtm, thlm, wprtp,
            wpthlp, rtp2_forcing, thlp2_forcing,
            rtpthlp_forcing, sclrm, wpsclrp,
            lhs_ta_wprtp2, lhs_ma, lhs_diff,
            rhs_ta_wprtp2, rhs_ta_wpthlp2,
            rhs_ta_wprtpthlp, rhs_ta_wpsclrp2,
            rhs_ta_wprtpsclrp, rhs_ta_wpthlpsclrp,
            dt, l_scalar_calc, l_lmm_stepping, l_implemented,
            tridiag_solve_method, stats,
            rtp2, thlp2, rtpthlp, sclrp2, sclrprtp, sclrpthlp,
            err_info,
        )
    )

    lhs_diff_uv = diffusion_zm_lhs(
        nzm, nzt, ngrdcol, gr, Kw9, Kw9_zm, nu_vert_res_dep.nu9,
        invrs_rho_ds_zm, rho_ds_zt,
    )

    if l_lmm_stepping:
        up2_old = up2
        vp2_old = vp2

    lhs_dp1_C4 = term_dp1_lhs(nzm, ngrdcol, gr, C4_1d, invrs_tau_C4_zm)
    lhs_dp1_C14 = term_dp1_lhs(nzm, ngrdcol, gr, C14_1d, invrs_tau_C14_zm)
    lhs_dp1 = (lhs_dp1_C4 + lhs_dp1_C14) * gamma_over_implicit_ts

    if iiPDF_type == iiPDF_new_hybrid:
        lhs = xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, lhs_ta_wpup2, lhs_ma, lhs_diff_uv, lhs_dp1)
        rhs, stats = xp2_xpyp_uv_rhs(
            nzm, nzt, ngrdcol, gr, xp2_xpyp_up2, dt,
            wp2, wpthvp, invrs_tau_C4_zm, invrs_tau_C14_zm,
            um, vm, upwp, vpwp, up2, vp2,
            thv_ds_zm, clubb_params[:, iC4], clubb_params[:, iC_uu_shr],
            clubb_params[:, iC_uu_buoy], clubb_params[:, iC14], lhs_splat_wp2,
            lhs_ta_wpup2, rhs_ta_wpup2, lhs_dp1_C4, lhs_dp1_C14, stats,
        )
        if l_ho_nontrad_coriolis:
            rhs = rhs - two * fcor_y[:, None] * upwp
            if l_sample:
                stats = stats.update("up2_nct", -two * fcor_y[:, None] * upwp)
        up2_solution, err_info, stats = xp2_xpyp_solve(
            nzm, ngrdcol, xp2_xpyp_up2_vp2, l_implemented, gr,
            tridiag_solve_method, stats, rhs, lhs, err_info,
        )
        up2 = up2_solution
        if l_lmm_stepping:
            up2 = one_half * (up2_old + up2)
        if l_sample:
            stats = stats_finalize_xp2_xpyp_terms(
                nzm, ngrdcol, xp2_xpyp_up2, gr, up2,
                gamma_over_implicit_ts * lhs_dp1_C14,
                gamma_over_implicit_ts * lhs_dp1_C4,
                lhs_diff_uv, lhs_ta_wpup2, lhs_ma, stats,
            )

        lhs = xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, lhs_ta_wpvp2, lhs_ma, lhs_diff_uv, lhs_dp1)
        rhs, stats = xp2_xpyp_uv_rhs(
            nzm, nzt, ngrdcol, gr, xp2_xpyp_vp2, dt,
            wp2, wpthvp, invrs_tau_C4_zm, invrs_tau_C14_zm,
            vm, um, vpwp, upwp, vp2, up2,
            thv_ds_zm, clubb_params[:, iC4], clubb_params[:, iC_uu_shr],
            clubb_params[:, iC_uu_buoy], clubb_params[:, iC14], lhs_splat_wp2,
            lhs_ta_wpvp2, rhs_ta_wpvp2, lhs_dp1_C4, lhs_dp1_C14, stats,
        )
        vp2_solution, err_info, stats = xp2_xpyp_solve(
            nzm, ngrdcol, xp2_xpyp_up2_vp2, l_implemented, gr,
            tridiag_solve_method, stats, rhs, lhs, err_info,
        )
        vp2 = vp2_solution
        if l_lmm_stepping:
            vp2 = one_half * (vp2_old + vp2)
        if l_sample:
            stats = stats_finalize_xp2_xpyp_terms(
                nzm, ngrdcol, xp2_xpyp_vp2, gr, vp2,
                gamma_over_implicit_ts * lhs_dp1_C14,
                gamma_over_implicit_ts * lhs_dp1_C4,
                lhs_diff_uv, lhs_ta_wpvp2, lhs_ma, stats,
            )
    else:
        lhs = xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, lhs_ta_wpup2, lhs_ma, lhs_diff_uv, lhs_dp1)
        rhs_up2, stats = xp2_xpyp_uv_rhs(
            nzm, nzt, ngrdcol, gr, xp2_xpyp_up2, dt,
            wp2, wpthvp, invrs_tau_C4_zm, invrs_tau_C14_zm,
            um, vm, upwp, vpwp, up2, vp2,
            thv_ds_zm, clubb_params[:, iC4], clubb_params[:, iC_uu_shr],
            clubb_params[:, iC_uu_buoy], clubb_params[:, iC14], lhs_splat_wp2,
            lhs_ta_wpup2, rhs_ta_wpup2, lhs_dp1_C4, lhs_dp1_C14, stats,
        )
        if l_ho_nontrad_coriolis:
            rhs_up2 = rhs_up2 - two * fcor_y[:, None] * upwp
            if l_sample:
                stats = stats.update("up2_nct", -two * fcor_y[:, None] * upwp)
        rhs_vp2, stats = xp2_xpyp_uv_rhs(
            nzm, nzt, ngrdcol, gr, xp2_xpyp_vp2, dt,
            wp2, wpthvp, invrs_tau_C4_zm, invrs_tau_C14_zm,
            vm, um, vpwp, upwp, vp2, up2,
            thv_ds_zm, clubb_params[:, iC4], clubb_params[:, iC_uu_shr],
            clubb_params[:, iC_uu_buoy], clubb_params[:, iC14], lhs_splat_wp2,
            lhs_ta_wpup2, rhs_ta_wpvp2, lhs_dp1_C4, lhs_dp1_C14, stats,
        )
        uv_rhs = jnp.stack((rhs_up2, rhs_vp2), axis=2)
        uv_solution, err_info, stats = xp2_xpyp_solve(
            nzm, ngrdcol, xp2_xpyp_up2_vp2, l_implemented, gr,
            tridiag_solve_method, stats, uv_rhs, lhs, err_info,
        )
        up2 = uv_solution[:, :, 0]
        vp2 = uv_solution[:, :, 1]
        if l_lmm_stepping:
            up2 = one_half * (up2_old + up2)
            vp2 = one_half * (vp2_old + vp2)
        if l_sample:
            stats = stats_finalize_xp2_xpyp_terms(
                nzm, ngrdcol, xp2_xpyp_up2, gr, up2,
                gamma_over_implicit_ts * lhs_dp1_C14,
                gamma_over_implicit_ts * lhs_dp1_C4,
                lhs_diff_uv, lhs_ta_wpup2, lhs_ma, stats,
            )
            stats = stats_finalize_xp2_xpyp_terms(
                nzm, ngrdcol, xp2_xpyp_vp2, gr, vp2,
                gamma_over_implicit_ts * lhs_dp1_C14,
                gamma_over_implicit_ts * lhs_dp1_C4,
                lhs_diff_uv, lhs_ta_wpup2, lhs_ma, stats,
            )

    if fill_holes_type != 0:
        rtp2, stats = pos_definite_variances(
            nzm, ngrdcol, gr, xp2_xpyp_rtp2, fill_holes_type,
            dt, rt_tol ** 2, rho_ds_zm, stats, rtp2,
        )
        thlp2, stats = pos_definite_variances(
            nzm, ngrdcol, gr, xp2_xpyp_thlp2, fill_holes_type,
            dt, thl_tol ** 2, rho_ds_zm, stats, thlp2,
        )
        up2, stats = pos_definite_variances(
            nzm, ngrdcol, gr, xp2_xpyp_up2, fill_holes_type,
            dt, w_tol_sqd, rho_ds_zm, stats, up2,
        )
        vp2, stats = pos_definite_variances(
            nzm, ngrdcol, gr, xp2_xpyp_vp2, fill_holes_type,
            dt, w_tol_sqd, rho_ds_zm, stats, vp2,
        )

    if l_min_xp2_from_corr_wx:
        threshold_array = jnp.maximum(
            rt_tol ** 2,
            wprtp ** 2 / (wp2 * max_mag_correlation_flux ** 2),
        )
    else:
        threshold_array = jnp.full((ngrdcol, nzm), rt_tol ** 2)
    rtp2, stats = clip_variance(
        nzm, ngrdcol, gr, xp2_xpyp_rtp2, dt, threshold_array, stats, rtp2,
    )

    if l_sample:
        stats = stats.begin_budget("rtp2_cl", rtp2 / dt)
    rtm_zm = jnp.maximum(zt2zm(nzm, nzt, ngrdcol, gr, rtm), zero_threshold)
    threshold_array = jnp.maximum(rt_tol ** 2, clubb_params[:, irtp2_clip_coef][:, None] * rtm_zm ** 2)
    rtp2 = jnp.minimum(rtp2, threshold_array)
    if l_sample:
        stats = stats.finalize_budget("rtp2_cl", rtp2 / dt, l_count_sample=False)

    if l_min_xp2_from_corr_wx:
        threshold_array = jnp.maximum(
            thl_tol ** 2,
            wpthlp ** 2 / (wp2 * max_mag_correlation_flux ** 2),
        )
    else:
        threshold_array = jnp.full((ngrdcol, nzm), thl_tol ** 2)
    thlp2, stats = clip_variance(
        nzm, ngrdcol, gr, xp2_xpyp_thlp2, dt, threshold_array, stats, thlp2,
    )

    threshold_array = jnp.full((ngrdcol, nzm), w_tol_sqd)
    up2, stats = clip_variance(
        nzm, ngrdcol, gr, xp2_xpyp_up2, dt, threshold_array, stats, up2,
    )
    if l_sample:
        stats = stats.begin_budget("up2_cl", up2 / dt)
    up2 = jnp.minimum(up2, 1000.0)
    if l_sample:
        stats = stats.finalize_budget("up2_cl", up2 / dt, l_count_sample=False)

    vp2, stats = clip_variance(
        nzm, ngrdcol, gr, xp2_xpyp_vp2, dt, threshold_array, stats, vp2,
    )
    if l_sample:
        stats = stats.begin_budget("vp2_cl", vp2 / dt)
    vp2 = jnp.minimum(vp2, 1000.0)
    if l_sample:
        stats = stats.finalize_budget("vp2_cl", vp2 / dt, l_count_sample=False)

    if l_sample:
        stats = stats.begin_budget("rtpthlp_cl", rtpthlp / dt)
    rtpthlp, _ = clip_covar(nzm, ngrdcol, xp2_xpyp_rtpthlp, rtp2, thlp2, rtpthlp)
    if l_sample:
        stats = stats.finalize_budget("rtpthlp_cl", rtpthlp / dt)

    if l_scalar_calc:
        if fill_holes_type != 0:
            for sclr in range(sclr_dim):
                sclrp2_s, stats = pos_definite_variances(
                    nzm, ngrdcol, gr, xp2_xpyp_sclrp2, fill_holes_type,
                    dt, sclr_tol[sclr] ** 2, rho_ds_zm, stats, sclrp2[:, :, sclr],
                )
                sclrp2 = sclrp2.at[:, :, sclr].set(sclrp2_s)
                if int(sclr_idx.iisclr_rt) == sclr + 1:
                    sclrprtp_s, stats = pos_definite_variances(
                        nzm, ngrdcol, gr, xp2_xpyp_sclrprtp, fill_holes_type,
                        dt, sclr_tol[sclr] ** 2, rho_ds_zm, stats, sclrprtp[:, :, sclr],
                    )
                    sclrprtp = sclrprtp.at[:, :, sclr].set(sclrprtp_s)
                if int(sclr_idx.iisclr_thl) == sclr + 1:
                    sclrpthlp_s, stats = pos_definite_variances(
                        nzm, ngrdcol, gr, xp2_xpyp_sclrpthlp, fill_holes_type,
                        dt, sclr_tol[sclr] ** 2, rho_ds_zm, stats, sclrpthlp[:, :, sclr],
                    )
                    sclrpthlp = sclrpthlp.at[:, :, sclr].set(sclrpthlp_s)

        for sclr in range(sclr_dim):
            threshold_array = jnp.full(
                (ngrdcol, nzm), sclr_tol[sclr] ** 2,
            )
            sclrp2_s, stats = clip_variance(
                nzm, ngrdcol, gr, clip_sclrp2, dt,
                threshold_array, stats, sclrp2[:, :, sclr],
            )
            sclrp2 = sclrp2.at[:, :, sclr].set(sclrp2_s)

        for sclr in range(sclr_dim):
            if int(sclr_idx.iisclr_rt) == sclr + 1:
                threshold_array = jnp.full(
                    (ngrdcol, nzm), sclr_tol[sclr] * rt_tol,
                )
                sclrprtp_s, stats = clip_variance(
                    nzm, ngrdcol, gr, clip_sclrprtp, dt,
                    threshold_array, stats, sclrprtp[:, :, sclr],
                )
                sclrprtp = sclrprtp.at[:, :, sclr].set(sclrprtp_s)
            else:
                sclrprtp_s, _ = clip_covar(
                    nzm, ngrdcol, clip_sclrprtp, sclrp2[:, :, sclr], rtp2, sclrprtp[:, :, sclr],
                )
                sclrprtp = sclrprtp.at[:, :, sclr].set(sclrprtp_s)

        for sclr in range(sclr_dim):
            if int(sclr_idx.iisclr_thl) == sclr + 1:
                threshold_array = jnp.full(
                    (ngrdcol, nzm), sclr_tol[sclr] * thl_tol,
                )
                sclrpthlp_s, stats = clip_variance(
                    nzm, ngrdcol, gr, clip_sclrpthlp, dt,
                    threshold_array, stats, sclrpthlp[:, :, sclr],
                )
                sclrpthlp = sclrpthlp.at[:, :, sclr].set(sclrpthlp_s)
            else:
                sclrpthlp_s, _ = clip_covar(
                    nzm, ngrdcol, clip_sclrpthlp, sclrp2[:, :, sclr], thlp2, sclrpthlp[:, :, sclr],
                )
                sclrpthlp = sclrpthlp.at[:, :, sclr].set(sclrpthlp_s)

    if clubb_at_least_debug_level(0):
        # Host I/O adaptation: full source dump using existing immutable snapshots.
        def print_xp_error(info, values, current, previous):
            import numpy as np
            with np.printoptions(threshold=np.inf):
                print(info.global_error_header_host(), file=sys.stderr)
                print("Error in advance_xp2_xpyp", file=sys.stderr)
                print("Intent(in)", file=sys.stderr)
                print('invrs_tau_xp2_zm = ', values['invrs_tau_xp2_zm'], file=sys.stderr)
                print('invrs_tau_C4_zm = ', values['invrs_tau_C4_zm'], file=sys.stderr)
                print('invrs_tau_C14_zm = ', values['invrs_tau_C14_zm'], file=sys.stderr)
                print('wm_zm = ', values['wm_zm'], file=sys.stderr)
                print('rtm = ', values['rtm'], file=sys.stderr)
                print('wprtp = ', values['wprtp'], file=sys.stderr)
                print('thlm = ', values['thlm'], file=sys.stderr)
                print('wpthlp = ', values['wpthlp'], file=sys.stderr)
                print('wpthvp = ', values['wpthvp'], file=sys.stderr)
                print('um = ', values['um'], file=sys.stderr)
                print('vm = ', values['vm'], file=sys.stderr)
                print('wp2 = ', values['wp2'], file=sys.stderr)
                print('wp3 = ', values['wp3'], file=sys.stderr)
                print('upwp = ', values['upwp'], file=sys.stderr)
                print('vpwp = ', values['vpwp'], file=sys.stderr)
                print('sigma_sqd_w = ', values['sigma_sqd_w'], file=sys.stderr)
                print('Kh_zt = ', values['Kh_zt'], file=sys.stderr)
                print('rtp2_forcing = ', values['rtp2_forcing'], file=sys.stderr)
                print('thlp2_forcing = ', values['thlp2_forcing'], file=sys.stderr)
                print('rtpthlp_forcing = ', values['rtpthlp_forcing'], file=sys.stderr)
                print('rho_ds_zm = ', values['rho_ds_zm'], file=sys.stderr)
                print('rho_ds_zt = ', values['rho_ds_zt'], file=sys.stderr)
                print('invrs_rho_ds_zm = ', values['invrs_rho_ds_zm'], file=sys.stderr)
                print('thv_ds_zm = ', values['thv_ds_zm'], file=sys.stderr)
                print('wp2_zt = ', values['wp2_zt'], file=sys.stderr)
                for sclr in range(sclr_dim):
                    print("sclrm = ", sclr+1, values["sclrm"][:, :, sclr], file=sys.stderr)
                    print("wpsclrp = ", sclr+1, values["wpsclrp"][:, :, sclr], file=sys.stderr)
                print("Intent(In/Out)", file=sys.stderr)
                for name in ("rtp2", "thlp2", "rtpthlp", "up2", "vp2"):
                    if l_lmm_stepping:
                        print(name + " (pre-solve) = ", previous[name], file=sys.stderr)
                    print(name + " = ", current[name], file=sys.stderr)
                for sclr in range(sclr_dim):
                    for name, label in (("sclrp2", "sclrp2"), ("sclrprtp", "sclrprtp"), ("sclrpthlp", "sclrthlp")):
                        if l_lmm_stepping:
                            print(label + " (pre-solve) = ", sclr+1, previous[name][:, :, sclr], file=sys.stderr)
                        print(label + " = ", sclr+1, current[name][:, :, sclr], file=sys.stderr)
        jax.lax.cond(
            err_info.any_fatal(),
            lambda _: jax.debug.callback(
                print_xp_error, err_info,
                {
                    'invrs_tau_xp2_zm': invrs_tau_xp2_zm,
                    'invrs_tau_C4_zm': invrs_tau_C4_zm,
                    'invrs_tau_C14_zm': invrs_tau_C14_zm,
                    'wm_zm': wm_zm,
                    'rtm': rtm,
                    'wprtp': wprtp,
                    'thlm': thlm,
                    'wpthlp': wpthlp,
                    'wpthvp': wpthvp,
                    'um': um,
                    'vm': vm,
                    'wp2': wp2,
                    'wp3': wp3,
                    'upwp': upwp,
                    'vpwp': vpwp,
                    'sigma_sqd_w': sigma_sqd_w,
                    'Kh_zt': Kh_zt,
                    'rtp2_forcing': rtp2_forcing,
                    'thlp2_forcing': thlp2_forcing,
                    'rtpthlp_forcing': rtpthlp_forcing,
                    'rho_ds_zm': rho_ds_zm,
                    'rho_ds_zt': rho_ds_zt,
                    'invrs_rho_ds_zm': invrs_rho_ds_zm,
                    'thv_ds_zm': thv_ds_zm,
                    'wp2_zt': wp2_zt,
                    'sclrm': sclrm,
                    'wpsclrp': wpsclrp,
                },
                {
                    'rtp2': rtp2,
                    'thlp2': thlp2,
                    'rtpthlp': rtpthlp,
                    'up2': up2,
                    'vp2': vp2,
                    'sclrp2': sclrp2,
                    'sclrprtp': sclrprtp,
                    'sclrpthlp': sclrpthlp,
                },
                {
                    'rtp2': initial_rtp2,
                    'thlp2': initial_thlp2,
                    'rtpthlp': initial_rtpthlp,
                    'up2': initial_up2,
                    'vp2': initial_vp2,
                    'sclrp2': initial_sclrp2,
                    'sclrprtp': initial_sclrprtp,
                    'sclrpthlp': initial_sclrpthlp,
                }, ordered=True,
            ),
            lambda _: None, operand=None,
        )

    l_return_initial_state = l_multiple_lhs_required
    rtp2 = jnp.where(l_return_initial_state, initial_rtp2, rtp2)
    thlp2 = jnp.where(l_return_initial_state, initial_thlp2, thlp2)
    rtpthlp = jnp.where(l_return_initial_state, initial_rtpthlp, rtpthlp)
    up2 = jnp.where(l_return_initial_state, initial_up2, up2)
    vp2 = jnp.where(l_return_initial_state, initial_vp2, vp2)
    sclrp2 = jnp.where(l_return_initial_state, initial_sclrp2, sclrp2)
    sclrprtp = jnp.where(l_return_initial_state, initial_sclrprtp, sclrprtp)
    sclrpthlp = jnp.where(l_return_initial_state, initial_sclrpthlp, sclrpthlp)

    return (
        rtp2, thlp2, rtpthlp, up2, vp2,
        sclrp2, sclrprtp, sclrpthlp, err_info, stats,
    )


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "nzt",
        "ngrdcol",
        "sclr_dim",
        "l_scalar_calc",
        "l_lmm_stepping",
        "l_implemented",
        "tridiag_solve_method",
    ),
)
def solve_xp2_xpyp_with_single_lhs(
    nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, sclr_idx,
    C2x, invrs_tau_xp2_zm, rtm, thlm, wprtp,
    wpthlp, rtp2_forcing, thlp2_forcing,
    rtpthlp_forcing, sclrm, wpsclrp,
    lhs_ta, lhs_ma, lhs_diff,
    rhs_ta_wprtp2, rhs_ta_wpthlp2,
    rhs_ta_wprtpthlp, rhs_ta_wpsclrp2,
    rhs_ta_wprtpsclrp, rhs_ta_wpthlpsclrp,
    dt, l_scalar_calc, l_lmm_stepping, l_implemented,
    tridiag_solve_method, stats,
    rtp2, thlp2, rtpthlp, sclrp2, sclrprtp, sclrpthlp, err_info,
):
    """This subroutine generates a single lhs matrix and multiple rhs matricies, then
    solves the system. This should only be used when the lhs matrices for
    <w'r_t'>, <w'th_l'>, <w'rt'thl'> ,<w'sclr'^2>, <w'sclr'rt'>,  and <w'sclr'thl'>
    are identical. Otherwise multiple lhs matrices and multiple solves are required.

    The ADG1 PDF must be in use for this subroutine to have the
    potential to be called.

    Returns:
        (rtp2, thlp2, rtpthlp, sclrp2, sclrprtp, sclrpthlp, err_info, stats)

    Fortran declarations (solve_xp2_xpyp_with_single_lhs in src/CLUBB_core/advance_xp2_xpyp_module.F90):

        ! -------- Input Variables --------
        integer, intent(in) :: &
          nzm,      & ! Number of momentum vertical levels
          nzt,      & ! Number of thermodynamic vertical levels
          ngrdcol,  & ! Number of grid columns
          sclr_dim    ! Number of passive scalars

        real( kind = core_rknd ), intent(in), dimension(sclr_dim) :: &
          sclr_tol          ! Threshold(s) on the passive scalars  [units vary]

        type (grid), intent(in) :: &
          gr

        type (sclr_idx_type), intent(in) :: &
          sclr_idx

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) :: &
          C2x,              &
          invrs_tau_xp2_zm, & ! Inverse time-scale for xp2 on momentum levels [1/s]
          wprtp,            & ! <w'r_t'> (momentum levels)            [(m/s)(kg/kg)]
          wpthlp,           & ! <w'th_l'> (momentum levels)           [(m K)/s]
          rtp2_forcing,     & ! <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
          thlp2_forcing,    & ! <th_l'^2> forcing (momentum levels)   [K^2/s]
          rtpthlp_forcing     ! <r_t'th_l'> forcing (momentum levels) [(kg/kg)K/s]

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) :: &
          rtm,              & ! Total water mixing ratio (t-levs)     [kg/kg]
          thlm                ! Liquid potential temp. (t-levs)       [K]

        logical, intent(in) :: &
          l_scalar_calc, &
          l_lmm_stepping, &
          l_implemented  ! True if CLUBB is being implemented and run in a host model

        real( kind = core_rknd ), intent(in) :: &
          dt             ! Model timestep                                [s]

        ! Passive scalar input
        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt,sclr_dim) :: &
          sclrm          ! Mean value; pass. scalar (t-levs.) [{sclr units}]

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm,sclr_dim) :: &
          wpsclrp        ! <w'sclr'> (momentum levels)        [m/s{sclr units}]

        real( kind = core_rknd ), intent(in), dimension(ndiags3,ngrdcol,nzm) :: &
          lhs_ta,       &
          lhs_diff,     & ! Diffusion contributions to lhs, dissipation term 2
          lhs_ma          ! Mean advection contributions to lhs

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) :: &
          rhs_ta_wprtp2,    & ! For <w'rt'^2>
          rhs_ta_wpthlp2,   & ! For <w'thl'^2>
          rhs_ta_wprtpthlp    ! For <w'rt'thl'>

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm,sclr_dim) :: &
          rhs_ta_wpsclrp2,      & ! For <w'sclr'^2>
          rhs_ta_wprtpsclrp,    & ! For <w'sclr'rt'><w'sclr'^2><w'sclr'thl'>
          rhs_ta_wpthlpsclrp      ! For <w'sclr'thl'>

        integer, intent(in) :: &
          tridiag_solve_method  ! Specifier for method to solve tridiagonal systems

        type(stats_type), intent(inout) :: &
              stats

        ! -------- In/Out Variables --------

        ! An attribute of (inout) is also needed to import the value of the variances
        ! at the surface.  Brian.  12/18/05.
        real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzm) :: &
          rtp2,    & ! <r_t'^2>                      [(kg/kg)^2]
          thlp2,   & ! <th_l'^2>                     [K^2]
          rtpthlp    ! <r_t'th_l'>                   [(kg K)/kg]

        real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzm,sclr_dim) :: &
          sclrp2, sclrprtp, sclrpthlp

        type(err_info_type), intent(inout) :: &
          err_info      ! err_info struct containing err_code and err_header
    """
    l_sample = stats.l_sample
    lhs_dp1 = term_dp1_lhs(nzm, ngrdcol, gr, C2x, invrs_tau_xp2_zm)
    lhs_dp1 = lhs_dp1 * gamma_over_implicit_ts
    lhs = xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, lhs_ta, lhs_ma, lhs_diff, lhs_dp1)
    rhs = jnp.zeros((ngrdcol, nzm, 3 + 3 * sclr_dim))

    rhs_s, stats = xp2_xpyp_rhs(
        nzm, nzt, ngrdcol, gr, xp2_xpyp_rtp2, dt,
        wprtp, wprtp, rtm, rtm, rtp2, rtp2_forcing,
        C2x, invrs_tau_xp2_zm, rt_tol ** 2, lhs_ta, rhs_ta_wprtp2, stats,
    )
    rhs = rhs.at[:, :, 0].set(rhs_s)
    rhs_s, stats = xp2_xpyp_rhs(
        nzm, nzt, ngrdcol, gr, xp2_xpyp_thlp2, dt,
        wpthlp, wpthlp, thlm, thlm, thlp2, thlp2_forcing,
        C2x, invrs_tau_xp2_zm, thl_tol ** 2, lhs_ta, rhs_ta_wpthlp2, stats,
    )
    rhs = rhs.at[:, :, 1].set(rhs_s)
    rhs_s, stats = xp2_xpyp_rhs(
        nzm, nzt, ngrdcol, gr, xp2_xpyp_rtpthlp, dt,
        wprtp, wpthlp, rtm, thlm, rtpthlp, rtpthlp_forcing,
        C2x, invrs_tau_xp2_zm, zero_threshold, lhs_ta, rhs_ta_wprtpthlp, stats,
    )
    rhs = rhs.at[:, :, 2].set(rhs_s)

    if l_scalar_calc:
        for sclr in range(sclr_dim):
            sclrp2_forcing = jnp.zeros((ngrdcol, nzm))
            rhs_s, stats = xp2_xpyp_rhs(
                nzm, nzt, ngrdcol, gr, xp2_xpyp_sclrp2, dt,
                wpsclrp[:, :, sclr], wpsclrp[:, :, sclr],
                sclrm[:, :, sclr], sclrm[:, :, sclr],
                sclrp2[:, :, sclr], sclrp2_forcing,
                C2x, invrs_tau_xp2_zm, sclr_tol[sclr] ** 2,
                lhs_ta, rhs_ta_wpsclrp2[:, :, sclr], stats,
            )
            rhs = rhs.at[:, :, 3 + sclr].set(rhs_s)

            if int(sclr_idx.iisclr_rt) == sclr + 1:
                sclrprtp_forcing = rtp2_forcing
                threshold = rt_tol ** 2
            else:
                sclrprtp_forcing = jnp.zeros((ngrdcol, nzm))
                threshold = zero_threshold
            rhs_s, stats = xp2_xpyp_rhs(
                nzm, nzt, ngrdcol, gr, xp2_xpyp_sclrprtp, dt,
                wpsclrp[:, :, sclr], wprtp,
                sclrm[:, :, sclr], rtm, sclrprtp[:, :, sclr],
                sclrprtp_forcing, C2x, invrs_tau_xp2_zm, threshold,
                lhs_ta, rhs_ta_wprtpsclrp[:, :, sclr], stats,
            )
            rhs = rhs.at[:, :, 3 + sclr + sclr_dim].set(rhs_s)

            if int(sclr_idx.iisclr_thl) == sclr + 1:
                sclrpthlp_forcing = thlp2_forcing
                threshold = thl_tol ** 2
            else:
                sclrpthlp_forcing = jnp.zeros((ngrdcol, nzm))
                threshold = zero_threshold
            rhs_s, stats = xp2_xpyp_rhs(
                nzm, nzt, ngrdcol, gr, xp2_xpyp_sclrpthlp, dt,
                wpsclrp[:, :, sclr], wpthlp,
                sclrm[:, :, sclr], thlm, sclrpthlp[:, :, sclr],
                sclrpthlp_forcing, C2x, invrs_tau_xp2_zm, threshold,
                lhs_ta, rhs_ta_wpthlpsclrp[:, :, sclr], stats,
            )
            rhs = rhs.at[:, :, 3 + sclr + 2 * sclr_dim].set(rhs_s)

    solution, err_info, stats = xp2_xpyp_solve(
        nzm, ngrdcol, xp2_xpyp_single_lhs, l_implemented, gr,
        tridiag_solve_method, stats, rhs, lhs, err_info,
    )
    if clubb_at_least_debug_level(0):
        jax.lax.cond(
            err_info.any_fatal(),
            lambda _: jax.debug.callback(
                lambda info: print(
                    info.global_error_header_host(),
                    'calling xp2_xpyp_solve in solve_xp2_xpyp_with_single_lhs',
                    sep="\n",
                    file=sys.stderr,
                ),
                err_info, ordered=True,
            ),
            lambda _: None, operand=None,
        )


    if l_lmm_stepping:
        rtp2 = one_half * (rtp2 + solution[:, :, 0])
        thlp2 = one_half * (thlp2 + solution[:, :, 1])
        rtpthlp = one_half * (rtpthlp + solution[:, :, 2])
        if sclr_dim > 0:
            for sclr in range(sclr_dim):
                sclrp2 = sclrp2.at[:, :, sclr].set(
                    one_half * (sclrp2[:, :, sclr] + solution[:, :, 3 + sclr])
                )
                sclrprtp = sclrprtp.at[:, :, sclr].set(
                    one_half * (
                        sclrprtp[:, :, sclr] + solution[:, :, 3 + sclr_dim + sclr]
                    )
                )
                sclrpthlp = sclrpthlp.at[:, :, sclr].set(
                    one_half * (
                        sclrpthlp[:, :, sclr] + solution[:, :, 3 + 2 * sclr_dim + sclr]
                    )
                )
    else:
        rtp2 = solution[:, :, 0]
        thlp2 = solution[:, :, 1]
        rtpthlp = solution[:, :, 2]
        if sclr_dim > 0:
            for sclr in range(sclr_dim):
                sclrp2 = sclrp2.at[:, :, sclr].set(solution[:, :, 3 + sclr])
                sclrprtp = sclrprtp.at[:, :, sclr].set(solution[:, :, 3 + sclr_dim + sclr])
                sclrpthlp = sclrpthlp.at[:, :, sclr].set(solution[:, :, 3 + 2 * sclr_dim + sclr])

    if l_sample:
        stats = stats_finalize_xp2_xpyp_terms(
            nzm, ngrdcol, xp2_xpyp_rtp2, gr, rtp2,
            lhs_dp1, lhs_dp1, lhs_diff, lhs_ta, lhs_ma, stats,
        )
        stats = stats_finalize_xp2_xpyp_terms(
            nzm, ngrdcol, xp2_xpyp_thlp2, gr, thlp2,
            lhs_dp1, lhs_dp1, lhs_diff, lhs_ta, lhs_ma, stats,
        )
        stats = stats_finalize_xp2_xpyp_terms(
            nzm, ngrdcol, xp2_xpyp_rtpthlp, gr, rtpthlp,
            lhs_dp1, lhs_dp1, lhs_diff, lhs_ta, lhs_ma, stats,
        )

    return rtp2, thlp2, rtpthlp, sclrp2, sclrprtp, sclrpthlp, err_info, stats


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "nzt",
        "ngrdcol",
        "sclr_dim",
        "iiPDF_type",
        "l_scalar_calc",
        "l_lmm_stepping",
        "l_implemented",
        "tridiag_solve_method",
    ),
)
def solve_xp2_xpyp_with_multiple_lhs(
    nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, sclr_idx,
    C2rt_1d, C2thl_1d, C2rtthl_1d, C2sclr_1d,
    invrs_tau_xp2_zm, rtm, thlm, wprtp, wpthlp,
    rtp2_forcing, thlp2_forcing, rtpthlp_forcing,
    sclrm, wpsclrp,
    lhs_ta_wprtp2, lhs_ta_wpthlp2,
    lhs_ta_wprtpthlp, lhs_ta_wpsclrp2,
    lhs_ta_wprtpsclrp, lhs_ta_wpthlpsclrp,
    lhs_ma, lhs_diff,
    rhs_ta_wprtp2, rhs_ta_wpthlp2,
    rhs_ta_wprtpthlp, rhs_ta_wpsclrp2,
    rhs_ta_wprtpsclrp, rhs_ta_wpthlpsclrp,
    dt, iiPDF_type, l_scalar_calc, l_lmm_stepping,
    l_implemented, tridiag_solve_method, stats,
    rtp2, thlp2, rtpthlp, sclrp2, sclrprtp, sclrpthlp, err_info,
):
    """This subroutine generates different lhs and rhs matrices to solve for.

    Returns:
        (rtp2, thlp2, rtpthlp, sclrp2, sclrprtp, sclrpthlp, err_info, stats)

    Fortran declarations (solve_xp2_xpyp_with_multiple_lhs in src/CLUBB_core/advance_xp2_xpyp_module.F90):

        !------------------------ Input Variables ------------------------
        integer, intent(in) :: &
          nzm,      & ! Number of momentum vertical levels
          nzt,      & ! Number of thermodynamic vertical levels
          ngrdcol,  & ! Number of grid columns
          sclr_dim    ! Number of passive scalars

        real( kind = core_rknd ), intent(in), dimension(sclr_dim) :: &
          sclr_tol          ! Threshold(s) on the passive scalars  [units vary]

        type (grid), intent(in) :: &
          gr

        type (sclr_idx_type), intent(in) :: &
          sclr_idx

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) ::  &
          C2rt_1d, C2thl_1d, C2rtthl_1d, C2sclr_1d, &
          invrs_tau_xp2_zm, & ! Inverse time-scale for xp2 on momentum levels [1/s]
          wprtp,           & ! <w'r_t'> (momentum levels)            [(m/s)(kg/kg)]
          wpthlp,          & ! <w'th_l'> (momentum levels)           [(m K)/s]
          rtp2_forcing,    & ! <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
          thlp2_forcing,   & ! <th_l'^2> forcing (momentum levels)   [K^2/s]
          rtpthlp_forcing    ! <r_t'th_l'> forcing (momentum levels) [(kg/kg)K/s]

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) ::  &
          rtm,             & ! Total water mixing ratio (t-levs)     [kg/kg]
          thlm               ! Liquid potential temp. (t-levs)       [K]

        integer, intent(in) :: &
          iiPDF_type    ! Selected option for the two-component normal (double
                        ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                        ! w, chi, and eta) portion of CLUBB's multivariate,
                        ! two-component PDF.

        logical, intent(in) :: &
          l_scalar_calc, &
          l_lmm_stepping, &
          l_implemented  ! True if CLUBB is being implemented and run in a host model

        real( kind = core_rknd ), intent(in) :: &
          dt             ! Model timestep                                [s]

        ! Passive scalar input
        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt,sclr_dim) ::  &
          sclrm          ! Mean value; pass. scalar (t-levs.) [{sclr units}]

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm,sclr_dim) ::  &
          wpsclrp        ! <w'sclr'> (momentum levels)        [m/s{sclr units}]

        real( kind = core_rknd ), intent(in), dimension(ndiags3,ngrdcol,nzm) :: &
          lhs_ta_wprtp2,    & ! Turbulent advection term for <w'rt'^2>
          lhs_ta_wpthlp2,   & ! Turbulent advection term for <w'thl'^2>
          lhs_ta_wprtpthlp, & ! Turbulent advection term for <w'rtp'thl'>
          lhs_diff,         & ! Diffusion contributions to lhs, dissipation term 2
          lhs_ma              ! Mean advection contributions to lhs

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) :: &
          rhs_ta_wprtp2,    & ! For <w'rt'^2>
          rhs_ta_wpthlp2,   & ! For <w'thl'^2>
          rhs_ta_wprtpthlp    ! For <w'rt'thl'>

        real( kind = core_rknd ), dimension(ndiags3,ngrdcol,nzm,sclr_dim), intent(in) :: &
          lhs_ta_wpsclrp2,    & ! For <w'sclr'^2>
          lhs_ta_wprtpsclrp,  & ! For <w'rt'sclr'>
          lhs_ta_wpthlpsclrp    ! For <w'thl'sclr'>

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm,sclr_dim) :: &
          rhs_ta_wpsclrp2,    & ! For <w'sclr'^2>
          rhs_ta_wprtpsclrp,  & ! For <w'sclr'rt'>
          rhs_ta_wpthlpsclrp    ! For <w'sclr'thl'>

        integer, intent(in) :: &
          tridiag_solve_method  ! Specifier for method to solve tridiagonal systems

        type(stats_type), intent(inout) :: &
          stats

        !------------------------ In/Out Variables ------------------------

        ! Input/Output variables
        ! An attribute of (inout) is also needed to import the value of the variances
        ! at the surface.  Brian.  12/18/05.
        real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzm) :: &
          rtp2,    & ! <r_t'^2>                      [(kg/kg)^2]
          thlp2,   & ! <th_l'^2>                     [K^2]
          rtpthlp    ! <r_t'th_l'>                   [(kg K)/kg]

        real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzm,sclr_dim) :: &
          sclrp2, sclrprtp, sclrpthlp

        type(err_info_type), intent(inout) :: &
          err_info      ! err_info struct containing err_code and err_header
    """
    l_sample = stats.l_sample
    lhs_dp1 = term_dp1_lhs(nzm, ngrdcol, gr, C2rt_1d, invrs_tau_xp2_zm)
    lhs_dp1 = lhs_dp1 * gamma_over_implicit_ts
    lhs = xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, lhs_ta_wprtp2, lhs_ma, lhs_diff, lhs_dp1)
    rhs, stats = xp2_xpyp_rhs(
        nzm, nzt, ngrdcol, gr, xp2_xpyp_rtp2, dt,
        wprtp, wprtp, rtm, rtm, rtp2, rtp2_forcing,
        C2rt_1d, invrs_tau_xp2_zm, rt_tol ** 2, lhs_ta_wprtp2, rhs_ta_wprtp2, stats,
    )
    rtp2_solution, err_info, stats = xp2_xpyp_solve(
        nzm, ngrdcol, xp2_xpyp_rtp2, l_implemented, gr,
        tridiag_solve_method, stats, rhs, lhs, err_info,
    )
    if clubb_at_least_debug_level(0):
        jax.lax.cond(
            err_info.any_fatal(),
            lambda _: jax.debug.callback(
                lambda info: print(
                    info.global_error_header_host(),
                    'calling xp2_xpyp_solve for rtp2 in solve_xp2_xpyp_with_multiple_lhs',
                    sep="\n",
                    file=sys.stderr,
                ),
                err_info, ordered=True,
            ),
            lambda _: None, operand=None,
        )

    rtp2 = one_half * (rtp2 + rtp2_solution) if l_lmm_stepping else rtp2_solution
    if l_sample:
        stats = stats_finalize_xp2_xpyp_terms(
            nzm, ngrdcol, xp2_xpyp_rtp2, gr, rtp2,
            lhs_dp1, lhs_dp1, lhs_diff, lhs_ta_wprtp2, lhs_ma, stats,
        )

    lhs_dp1 = term_dp1_lhs(nzm, ngrdcol, gr, C2thl_1d, invrs_tau_xp2_zm)
    lhs_dp1 = lhs_dp1 * gamma_over_implicit_ts
    lhs = xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, lhs_ta_wpthlp2, lhs_ma, lhs_diff, lhs_dp1)
    rhs, stats = xp2_xpyp_rhs(
        nzm, nzt, ngrdcol, gr, xp2_xpyp_thlp2, dt,
        wpthlp, wpthlp, thlm, thlm, thlp2, thlp2_forcing,
        C2thl_1d, invrs_tau_xp2_zm, thl_tol ** 2, lhs_ta_wpthlp2, rhs_ta_wpthlp2, stats,
    )
    thlp2_solution, err_info, stats = xp2_xpyp_solve(
        nzm, ngrdcol, xp2_xpyp_thlp2, l_implemented, gr,
        tridiag_solve_method, stats, rhs, lhs, err_info,
    )
    if clubb_at_least_debug_level(0):
        jax.lax.cond(
            err_info.any_fatal(),
            lambda _: jax.debug.callback(
                lambda info: print(
                    info.global_error_header_host(),
                    'calling xp2_xpyp_solve for thlp2 in solve_xp2_xpyp_with_multiple_lhs',
                    sep="\n",
                    file=sys.stderr,
                ),
                err_info, ordered=True,
            ),
            lambda _: None, operand=None,
        )

    thlp2 = one_half * (thlp2 + thlp2_solution) if l_lmm_stepping else thlp2_solution
    if l_sample:
        stats = stats_finalize_xp2_xpyp_terms(
            nzm, ngrdcol, xp2_xpyp_thlp2, gr, thlp2,
            lhs_dp1, lhs_dp1, lhs_diff, lhs_ta_wpthlp2, lhs_ma, stats,
        )

    lhs_dp1 = term_dp1_lhs(nzm, ngrdcol, gr, C2rtthl_1d, invrs_tau_xp2_zm)
    lhs_dp1 = lhs_dp1 * gamma_over_implicit_ts
    lhs = xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, lhs_ta_wprtpthlp, lhs_ma, lhs_diff, lhs_dp1)
    rhs, stats = xp2_xpyp_rhs(
        nzm, nzt, ngrdcol, gr, xp2_xpyp_rtpthlp, dt,
        wprtp, wpthlp, rtm, thlm, rtpthlp, rtpthlp_forcing,
        C2rtthl_1d, invrs_tau_xp2_zm, zero_threshold, lhs_ta_wprtpthlp,
        rhs_ta_wprtpthlp, stats,
    )
    rtpthlp_solution, err_info, stats = xp2_xpyp_solve(
        nzm, ngrdcol, xp2_xpyp_rtpthlp, l_implemented, gr,
        tridiag_solve_method, stats, rhs, lhs, err_info,
    )
    if clubb_at_least_debug_level(0):
        jax.lax.cond(
            err_info.any_fatal(),
            lambda _: jax.debug.callback(
                lambda info: print(
                    info.global_error_header_host(),
                    'calling xp2_xpyp_solve for rtpthlp in solve_xp2_xpyp_with_multiple_lhs',
                    sep="\n",
                    file=sys.stderr,
                ),
                err_info, ordered=True,
            ),
            lambda _: None, operand=None,
        )

    rtpthlp = one_half * (rtpthlp + rtpthlp_solution) if l_lmm_stepping else rtpthlp_solution
    if l_sample:
        stats = stats_finalize_xp2_xpyp_terms(
            nzm, ngrdcol, xp2_xpyp_rtpthlp, gr, rtpthlp,
            lhs_dp1, lhs_dp1, lhs_diff, lhs_ta_wprtpthlp, lhs_ma, stats,
        )

    if l_scalar_calc:
        lhs_dp1 = term_dp1_lhs(nzm, ngrdcol, gr, C2sclr_1d, invrs_tau_xp2_zm)
        lhs_dp1 = lhs_dp1 * gamma_over_implicit_ts

        if iiPDF_type != iiPDF_ADG1:
            for sclr in range(sclr_dim):
                sclrp2_forcing = jnp.zeros((ngrdcol, nzm))
                lhs = xp2_xpyp_lhs(
                    nzm, ngrdcol, dt, gr, lhs_ta_wpsclrp2[:, :, :, sclr],
                    lhs_ma, lhs_diff, lhs_dp1,
                )
                rhs, stats = xp2_xpyp_rhs(
                    nzm, nzt, ngrdcol, gr, xp2_xpyp_sclrp2, dt,
                    wpsclrp[:, :, sclr], wpsclrp[:, :, sclr],
                    sclrm[:, :, sclr], sclrm[:, :, sclr],
                    sclrp2[:, :, sclr], sclrp2_forcing,
                    C2sclr_1d, invrs_tau_xp2_zm, sclr_tol[sclr] ** 2,
                    lhs_ta_wpsclrp2[:, :, :, sclr], rhs_ta_wpsclrp2[:, :, sclr],
                    stats,
                )
                solution, err_info, stats = xp2_xpyp_solve(
                    nzm, ngrdcol, xp2_xpyp_scalars, l_implemented, gr,
                    tridiag_solve_method, stats, rhs, lhs, err_info,
                )
                if clubb_at_least_debug_level(0):
                    jax.lax.cond(
                        err_info.any_fatal(),
                        lambda _: jax.debug.callback(
                            lambda info, index=sclr: print(
                                info.global_error_header_host(),
                                'calling xp2_xpyp_solve for sclrp2 with sclr index ' + str(index + 1) + " in solve_xp2_xpyp_with_multiple_lhs",
                                sep="\n",
                                file=sys.stderr,
                            ),
                            err_info, ordered=True,
                        ),
                        lambda _: None, operand=None,
                    )

                sclrp2 = sclrp2.at[:, :, sclr].set(
                    one_half * (sclrp2[:, :, sclr] + solution)
                    if l_lmm_stepping else solution
                )

                if int(sclr_idx.iisclr_rt) == sclr + 1:
                    sclrprtp_forcing = rtp2_forcing
                    threshold = rt_tol ** 2
                else:
                    sclrprtp_forcing = jnp.zeros((ngrdcol, nzm))
                    threshold = zero_threshold
                lhs = xp2_xpyp_lhs(
                    nzm, ngrdcol, dt, gr, lhs_ta_wprtpsclrp[:, :, :, sclr],
                    lhs_ma, lhs_diff, lhs_dp1,
                )
                rhs, stats = xp2_xpyp_rhs(
                    nzm, nzt, ngrdcol, gr, xp2_xpyp_sclrprtp, dt,
                    wpsclrp[:, :, sclr], wprtp,
                    sclrm[:, :, sclr], rtm, sclrprtp[:, :, sclr],
                    sclrprtp_forcing, C2sclr_1d, invrs_tau_xp2_zm, threshold,
                    lhs_ta_wprtpsclrp[:, :, :, sclr], rhs_ta_wprtpsclrp[:, :, sclr],
                    stats,
                )
                solution, err_info, stats = xp2_xpyp_solve(
                    nzm, ngrdcol, xp2_xpyp_scalars, l_implemented, gr,
                    tridiag_solve_method, stats, rhs, lhs, err_info,
                )
                if clubb_at_least_debug_level(0):
                    jax.lax.cond(
                        err_info.any_fatal(),
                        lambda _: jax.debug.callback(
                            lambda info, index=sclr: print(
                                info.global_error_header_host(),
                                'calling xp2_xpyp_solve for sclrprtp with sclr index ' + str(index + 1) + " in solve_xp2_xpyp_with_multiple_lhs",
                                sep="\n",
                                file=sys.stderr,
                            ),
                            err_info, ordered=True,
                        ),
                        lambda _: None, operand=None,
                    )

                sclrprtp = sclrprtp.at[:, :, sclr].set(
                    one_half * (sclrprtp[:, :, sclr] + solution)
                    if l_lmm_stepping else solution
                )

                if int(sclr_idx.iisclr_thl) == sclr + 1:
                    sclrpthlp_forcing = thlp2_forcing
                    threshold = thl_tol ** 2
                else:
                    sclrpthlp_forcing = jnp.zeros((ngrdcol, nzm))
                    threshold = zero_threshold
                lhs = xp2_xpyp_lhs(
                    nzm, ngrdcol, dt, gr, lhs_ta_wpthlpsclrp[:, :, :, sclr],
                    lhs_ma, lhs_diff, lhs_dp1,
                )
                rhs, stats = xp2_xpyp_rhs(
                    nzm, nzt, ngrdcol, gr, xp2_xpyp_sclrpthlp, dt,
                    wpsclrp[:, :, sclr], wpthlp,
                    sclrm[:, :, sclr], thlm, sclrpthlp[:, :, sclr],
                    sclrpthlp_forcing, C2sclr_1d, invrs_tau_xp2_zm, threshold,
                    lhs_ta_wpthlpsclrp[:, :, :, sclr], rhs_ta_wpthlpsclrp[:, :, sclr],
                    stats,
                )
                solution, err_info, stats = xp2_xpyp_solve(
                    nzm, ngrdcol, xp2_xpyp_scalars, l_implemented, gr,
                    tridiag_solve_method, stats, rhs, lhs, err_info,
                )
                if clubb_at_least_debug_level(0):
                    jax.lax.cond(
                        err_info.any_fatal(),
                        lambda _: jax.debug.callback(
                            lambda info, index=sclr: print(
                                info.global_error_header_host(),
                                'calling xp2_xpyp_solve for sclrpthlp with sclr index ' + str(index + 1) + " in solve_xp2_xpyp_with_multiple_lhs",
                                sep="\n",
                                file=sys.stderr,
                            ),
                            err_info, ordered=True,
                        ),
                        lambda _: None, operand=None,
                    )

                sclrpthlp = sclrpthlp.at[:, :, sclr].set(
                    one_half * (sclrpthlp[:, :, sclr] + solution)
                    if l_lmm_stepping else solution
                )
        else:
            lhs = xp2_xpyp_lhs(
                nzm, ngrdcol, dt, gr, lhs_ta_wpsclrp2[:, :, :, 0],
                lhs_ma, lhs_diff, lhs_dp1,
            )
            sclr_rhs = jnp.zeros((ngrdcol, nzm, 3 * sclr_dim))
            for sclr in range(sclr_dim):
                sclrp2_forcing = jnp.zeros((ngrdcol, nzm))
                rhs_s, stats = xp2_xpyp_rhs(
                    nzm, nzt, ngrdcol, gr, xp2_xpyp_sclrp2, dt,
                    wpsclrp[:, :, sclr], wpsclrp[:, :, sclr],
                    sclrm[:, :, sclr], sclrm[:, :, sclr],
                    sclrp2[:, :, sclr], sclrp2_forcing,
                    C2sclr_1d, invrs_tau_xp2_zm, sclr_tol[sclr] ** 2,
                    lhs_ta_wpsclrp2[:, :, :, 0], rhs_ta_wpsclrp2[:, :, sclr],
                    stats,
                )
                sclr_rhs = sclr_rhs.at[:, :, sclr].set(rhs_s)
                if int(sclr_idx.iisclr_rt) == sclr + 1:
                    sclrprtp_forcing = rtp2_forcing
                    threshold = rt_tol ** 2
                else:
                    sclrprtp_forcing = jnp.zeros((ngrdcol, nzm))
                    threshold = zero_threshold
                rhs_s, stats = xp2_xpyp_rhs(
                    nzm, nzt, ngrdcol, gr, xp2_xpyp_sclrprtp, dt,
                    wpsclrp[:, :, sclr], wprtp,
                    sclrm[:, :, sclr], rtm, sclrprtp[:, :, sclr],
                    sclrprtp_forcing, C2sclr_1d, invrs_tau_xp2_zm, threshold,
                    lhs_ta_wpsclrp2[:, :, :, 0], rhs_ta_wprtpsclrp[:, :, sclr],
                    stats,
                )
                sclr_rhs = sclr_rhs.at[:, :, sclr + sclr_dim].set(rhs_s)
                if int(sclr_idx.iisclr_thl) == sclr + 1:
                    sclrpthlp_forcing = thlp2_forcing
                    threshold = thl_tol ** 2
                else:
                    sclrpthlp_forcing = jnp.zeros((ngrdcol, nzm))
                    threshold = zero_threshold
                rhs_s, stats = xp2_xpyp_rhs(
                    nzm, nzt, ngrdcol, gr, xp2_xpyp_sclrpthlp, dt,
                    wpsclrp[:, :, sclr], wpthlp,
                    sclrm[:, :, sclr], thlm, sclrpthlp[:, :, sclr],
                    sclrpthlp_forcing, C2sclr_1d, invrs_tau_xp2_zm, threshold,
                    lhs_ta_wpsclrp2[:, :, :, 0], rhs_ta_wpthlpsclrp[:, :, sclr],
                    stats,
                )
                sclr_rhs = sclr_rhs.at[:, :, sclr + 2 * sclr_dim].set(rhs_s)
            sclr_solution, err_info, stats = xp2_xpyp_solve(
                nzm, ngrdcol, xp2_xpyp_scalars, l_implemented, gr,
                tridiag_solve_method, stats, sclr_rhs, lhs, err_info,
            )
            if clubb_at_least_debug_level(0):
                jax.lax.cond(
                    err_info.any_fatal(),
                    lambda _: jax.debug.callback(
                        lambda info, index=sclr: print(
                            info.global_error_header_host(),
                            'calling xp2_xpyp_solve for sclr with sclr index ' + str(index + 1) + " in solve_xp2_xpyp_with_multiple_lhs",
                            sep="\n",
                            file=sys.stderr,
                        ),
                        err_info, ordered=True,
                    ),
                    lambda _: None, operand=None,
                )

            for sclr in range(sclr_dim):
                if l_lmm_stepping:
                    sclrp2 = sclrp2.at[:, :, sclr].set(
                        one_half * (sclrp2[:, :, sclr] + sclr_solution[:, :, sclr])
                    )
                    sclrprtp = sclrprtp.at[:, :, sclr].set(
                        one_half * (
                            sclrprtp[:, :, sclr] + sclr_solution[:, :, sclr_dim + sclr]
                        )
                    )
                    sclrpthlp = sclrpthlp.at[:, :, sclr].set(
                        one_half * (
                            sclrpthlp[:, :, sclr] + sclr_solution[:, :, 2 * sclr_dim + sclr]
                        )
                    )
                else:
                    sclrp2 = sclrp2.at[:, :, sclr].set(sclr_solution[:, :, sclr])
                    sclrprtp = sclrprtp.at[:, :, sclr].set(sclr_solution[:, :, sclr_dim + sclr])
                    sclrpthlp = sclrpthlp.at[:, :, sclr].set(sclr_solution[:, :, 2 * sclr_dim + sclr])

    return rtp2, thlp2, rtpthlp, sclrp2, sclrprtp, sclrpthlp, err_info, stats


@partial(jax.jit, static_argnames=("nzm", "ngrdcol"))
def xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, lhs_ta, lhs_ma, lhs_diff, lhs_dp1):
    """Compute LHS tridiagonal matrix for a variance or covariance term

    Notes:
        For LHS turbulent advection terms:
            An "over-implicit" weighted time step is applied to this term.
            The weight of the implicit portion of this term is controlled by
            the factor gamma_over_implicit_ts (abbreviated "gamma" in the
            expression below).  A factor is added to the right-hand side of
            the equation in order to balance a weight that is not equal to 1,
            such that:
                 -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
            where X is the variable that is being solved for in a predictive
            equation (<x'^2> or <x'y'> in this case), y(t) is the linearized
            portion of the term that gets treated implicitly, and RHS is the
            portion of the term that is always treated explicitly.  A weight
            of greater than 1 can be applied to make the term more
            numerically stable.

    Returns:
        lhs

    Fortran declarations (xp2_xpyp_lhs in src/CLUBB_core/advance_xp2_xpyp_module.F90):

        !------------------- Input Variables -------------------
        integer, intent(in) :: &
          nzm, &
          ngrdcol

        type( grid ), intent(in) :: &
          gr

        real( kind = core_rknd ), dimension(ndiags3,ngrdcol,nzm), intent(in) :: &
         lhs_ta     ! Turbulent advection contributions to lhs

        real( kind = core_rknd ), dimension(ndiags3,ngrdcol,nzm), intent(in) :: &
          lhs_diff, & ! Diffusion contributions to lhs, dissipation term 2
          lhs_ma      ! Mean advection contributions to lhs

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          lhs_dp1   ! LHS dissipation term 1

        !------------------- Output Variables -------------------
        real( kind = core_rknd ), dimension(ndiags3,ngrdcol,nzm), intent(out) :: &
          lhs         ! Implicit contributions to the term

        real( kind = core_rknd ), intent(in) :: &
          dt
    """
    lhs = jnp.zeros((ndiags3, ngrdcol, nzm))
    interior = slice(1, nzm - 1)
    lhs = lhs.at[0, :, interior].set(
        lhs_diff[0, :, interior]
        + lhs_ma[0, :, interior]
        + lhs_ta[0, :, interior] * gamma_over_implicit_ts
    )
    lhs = lhs.at[1, :, interior].set(
        lhs_diff[1, :, interior]
        + lhs_ma[1, :, interior]
        + lhs_ta[1, :, interior] * gamma_over_implicit_ts
        + lhs_dp1[:, interior]
        + one / dt
    )
    lhs = lhs.at[2, :, interior].set(
        lhs_diff[2, :, interior]
        + lhs_ma[2, :, interior]
        + lhs_ta[2, :, interior] * gamma_over_implicit_ts
    )

    lhs = lhs.at[0, :, gr.k_lb_zm].set(zero)
    lhs = lhs.at[1, :, gr.k_lb_zm].set(one)
    lhs = lhs.at[2, :, gr.k_lb_zm].set(zero)
    lhs = lhs.at[0, :, gr.k_ub_zm].set(zero)
    lhs = lhs.at[1, :, gr.k_ub_zm].set(one)
    lhs = lhs.at[2, :, gr.k_ub_zm].set(zero)
    return lhs


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "ngrdcol",
        "solve_type",
        "l_implemented",
        "tridiag_solve_method",
    ),
)
def xp2_xpyp_solve(
    nzm, ngrdcol, solve_type, l_implemented, gr, tridiag_solve_method,
    stats, rhs, lhs, err_info,
):
    """Solve a tridiagonal system

    Returns:
        (solution, err_info, stats)

    Fortran declarations (xp2_xpyp_solve in src/CLUBB_core/advance_xp2_xpyp_module.F90):

        ! ---------------------- Input variables ----------------------
        integer, intent(in) :: &
          nzm, &
          ngrdcol

        integer, intent(in) :: &
          nrhs  ! Number of right hand side vectors

        logical, intent(in) :: &
          l_implemented  ! True if CLUBB is being implemented and run in a host model

        integer, intent(in) :: &
          solve_type ! Variable(s) description

        type( grid ), intent(in) :: &
          gr

        integer, intent(in) :: &
          tridiag_solve_method  ! Specifier for method to solve tridiagonal systems

        type(stats_type), intent(inout) :: &
          stats

        !------------------- Input/Output variables ------------------
        ! ---------------------- Input/Ouput variables ----------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzm,nrhs), intent(inout) :: &
          rhs  ! Explicit contributions to x variance/covariance term [units vary]

        real( kind = core_rknd ), dimension(ndiags3,ngrdcol,nzm), intent(inout) :: &
          lhs  ! Implicit contributions to x variance/covariance term [units vary]

        type(err_info_type), intent(inout) :: &
          err_info      ! err_info struct containing err_code and err_header

        ! ---------------------- Output Variables ----------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzm,nrhs), intent(out) :: &
          xapxbp ! Computed value of the variable(s) at <t+1> [units vary]
    """
    l_sample = stats.l_sample
    l_single_lhs_solve = False
    l_need_rcond = False

    if solve_type == xp2_xpyp_rtp2:
        solve_type_str = "rtp2"
        l_need_rcond = bool(l_sample and stats.var_on_stats_list("rtp2_matrix_condt_num"))
    elif solve_type == xp2_xpyp_thlp2:
        solve_type_str = "thlp2"
        l_need_rcond = bool(l_sample and stats.var_on_stats_list("thlp2_matrix_condt_num"))
    elif solve_type == xp2_xpyp_rtpthlp:
        solve_type_str = "rtpthlp"
        l_need_rcond = bool(l_sample and stats.var_on_stats_list("rtpthlp_matrix_condt_num"))
    elif solve_type == xp2_xpyp_up2_vp2:
        solve_type_str = "up2_vp2"
        l_need_rcond = bool(l_sample and stats.var_on_stats_list("up2_vp2_matrix_condt_num"))
    elif solve_type == xp2_xpyp_single_lhs:
        l_single_lhs_solve = True
        solve_type_str = "xp2_xpyp_single_lhs"
        l_need_rcond = (
            bool(l_sample and stats.var_on_stats_list("rtp2_matrix_condt_num"))
            or bool(l_sample and stats.var_on_stats_list("thlp2_matrix_condt_num"))
            or bool(l_sample and stats.var_on_stats_list("rtpthlp_matrix_condt_num"))
        )
    else:
        solve_type_str = "scalar"

    lhs_solve = lhs
    rhs_solve = rhs
    if l_force_descending_solves and gr.grid_dir_indx > 0:
        lhs_solve = lhs_solve[::-1, :, ::-1]
        if rhs_solve.ndim == 2:
            rhs_solve = rhs_solve[:, ::-1]
        else:
            rhs_solve = rhs_solve[:, ::-1, :]

    err_info, solution, rcond = tridiag_solve(
        solve_type_str, tridiag_solve_method, ngrdcol, nzm,
        lhs_solve, rhs_solve, err_info,
        l_use_rcond=l_need_rcond,
        l_implemented=l_implemented,
    )
    if clubb_at_least_debug_level(0):
        jax.lax.cond(
            err_info.any_fatal(),
            lambda _: jax.debug.callback(
                lambda info: print(
                    info.global_error_header_host(),
                    'calling tridiag_solve in xp2_xpyp_solve',
                    sep="\n",
                    file=sys.stderr,
                ),
                err_info, ordered=True,
            ),
            lambda _: None, operand=None,
        )


    if l_need_rcond and l_sample:
        cond = one / rcond

        def record_rcond(stats_in):
            if l_single_lhs_solve:
                stats_out = stats_in.update("rtp2_matrix_condt_num", cond)
                stats_out = stats_out.update("thlp2_matrix_condt_num", cond)
                return stats_out.update("rtpthlp_matrix_condt_num", cond)
            if solve_type == xp2_xpyp_rtp2:
                return stats_in.update("rtp2_matrix_condt_num", cond)
            if solve_type == xp2_xpyp_thlp2:
                return stats_in.update("thlp2_matrix_condt_num", cond)
            if solve_type == xp2_xpyp_rtpthlp:
                return stats_in.update("rtpthlp_matrix_condt_num", cond)
            if solve_type == xp2_xpyp_up2_vp2:
                return stats_in.update("up2_vp2_matrix_condt_num", cond)
            return stats_in

        stats = jax.lax.cond(
            err_info.has_fatal(),
            lambda stats_in: stats_in,
            record_rcond,
            stats,
        )

    if l_force_descending_solves and gr.grid_dir_indx > 0:
        if solution.ndim == 2:
            solution = solution[:, ::-1]
        else:
            solution = solution[:, ::-1, :]

    return solution, err_info, stats


@partial(jax.jit, static_argnames=("nzm", "ngrdcol", "solve_type"))
def stats_finalize_xp2_xpyp_terms(
    nzm, ngrdcol, solve_type, gr, xapxbp,
    lhs_dp1, lhs_dp1_pr1, lhs_diff, lhs_ta, lhs_ma, stats,
):
    """
    Fortran declarations (stats_finalize_xp2_xpyp_terms in src/CLUBB_core/advance_xp2_xpyp_module.F90):

        !----------------------- Input variables -----------------------
        integer, intent(in) :: &
          nzm, &
          ngrdcol, &
          solve_type

        type(grid), intent(in) :: &
          gr

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          xapxbp, &
          lhs_dp1, &
          lhs_dp1_pr1

        real( kind = core_rknd ), dimension(ndiags3,ngrdcol,nzm), intent(in) :: &
          lhs_diff, &
          lhs_ta, &
          lhs_ma

        !----------------------- InOut variables -----------------------
        type(stats_type), intent(inout) :: &
          stats
    """
    l_sample = stats.l_sample
    if not l_sample:
        return stats

    if solve_type == xp2_xpyp_rtp2:
        name_dp1, name_dp2, name_ta, name_ma, name_pr1 = "rtp2_dp1", "rtp2_dp2", "rtp2_ta", "rtp2_ma", ""
    elif solve_type == xp2_xpyp_thlp2:
        name_dp1, name_dp2, name_ta, name_ma, name_pr1 = "thlp2_dp1", "thlp2_dp2", "thlp2_ta", "thlp2_ma", ""
    elif solve_type == xp2_xpyp_rtpthlp:
        name_dp1, name_dp2, name_ta, name_ma, name_pr1 = "rtpthlp_dp1", "rtpthlp_dp2", "rtpthlp_ta", "rtpthlp_ma", ""
    elif solve_type == xp2_xpyp_up2:
        name_dp1, name_dp2, name_ta, name_ma, name_pr1 = "up2_dp1", "up2_dp2", "up2_ta", "up2_ma", "up2_pr1"
    elif solve_type == xp2_xpyp_vp2:
        name_dp1, name_dp2, name_ta, name_ma, name_pr1 = "vp2_dp1", "vp2_dp2", "vp2_ta", "vp2_ma", "vp2_pr1"
    else:
        return stats

    grid_dir = int(gr.grid_dir_indx)
    interior = slice(1, nzm - 1)
    km1 = slice(1 - grid_dir, nzm - 1 - grid_dir)
    kp1 = slice(1 + grid_dir, nzm - 1 + grid_dir)

    stats_tmp = jnp.zeros((ngrdcol, nzm))
    stats_tmp = stats_tmp.at[:, interior].set(
        -lhs_dp1[:, interior] * xapxbp[:, interior]
    )
    stats = stats.finalize_budget(name_dp1, stats_tmp)

    if solve_type in (xp2_xpyp_up2, xp2_xpyp_vp2):
        stats_tmp = jnp.zeros((ngrdcol, nzm))
        stats_tmp = stats_tmp.at[:, interior].set(
            -lhs_dp1_pr1[:, interior] * xapxbp[:, interior]
        )
        stats = stats.finalize_budget(name_pr1, stats_tmp)

    stats_tmp = jnp.zeros((ngrdcol, nzm))
    stats_tmp = stats_tmp.at[:, interior].set(
        (-gamma_over_implicit_ts * lhs_ta[1 + grid_dir, :, interior]) * xapxbp[:, km1]
        + (-gamma_over_implicit_ts * lhs_ta[1, :, interior]) * xapxbp[:, interior]
        + (-gamma_over_implicit_ts * lhs_ta[1 - grid_dir, :, interior]) * xapxbp[:, kp1]
    )
    stats = stats.finalize_budget(name_ta, stats_tmp)

    stats_tmp = jnp.zeros((ngrdcol, nzm))
    stats_tmp = stats_tmp.at[:, interior].set(
        (-lhs_diff[1 + grid_dir, :, interior]) * xapxbp[:, km1]
        + (-lhs_diff[1, :, interior]) * xapxbp[:, interior]
        + (-lhs_diff[1 - grid_dir, :, interior]) * xapxbp[:, kp1]
    )
    stats = stats.update(name_dp2, stats_tmp)

    stats_tmp = jnp.zeros((ngrdcol, nzm))
    stats_tmp = stats_tmp.at[:, interior].set(
        (-lhs_ma[1 + grid_dir, :, interior]) * xapxbp[:, km1]
        + (-lhs_ma[1, :, interior]) * xapxbp[:, interior]
        + (-lhs_ma[1 - grid_dir, :, interior]) * xapxbp[:, kp1]
    )
    stats = stats.update(name_ma, stats_tmp)

    return stats


@partial(jax.jit, static_argnames=("nzm", "nzt", "ngrdcol", "solve_type"))
def xp2_xpyp_uv_rhs(
    nzm, nzt, ngrdcol, gr, solve_type, dt,
    wp2, wpthvp,
    invrs_tau_C4_zm, invrs_tau_C14_zm,
    xam, xbm, wpxap, wpxbp, xap2, xbp2,
    thv_ds_zm, C4, C_uu_shr, C_uu_buoy, C14, lhs_splat_wp2,
    lhs_ta, rhs_ta,
    lhs_dp1_C4, lhs_dp1_C14, stats,
):
    """Explicit contributions to u'^2 or v'^2

    Notes:
       For LHS turbulent advection (ta) term:
           An "over-implicit" weighted time step is applied to this term.
           The weight of the implicit portion of this term is controlled by
           the factor gamma_over_implicit_ts (abbreviated "gamma" in the
           expression below).  A factor is added to the right-hand side of
           the equation in order to balance a weight that is not equal to 1,
           such that:
                -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
           where X is the variable that is being solved for in a predictive
           equation (x'^2 or x'y' in this case), y(t) is the linearized
           portion of the term that gets treated implicitly, and RHS is the
           portion of the term that is always treated explicitly.  A weight
           of greater than 1 can be applied to make the term more
           numerically stable.

    Returns:
        (rhs, stats)

    Fortran declarations (xp2_xpyp_uv_rhs in src/CLUBB_core/advance_xp2_xpyp_module.F90):

        !------------------- Input Variables -------------------
        integer, intent(in) :: &
          nzm, &
          nzt, &
          ngrdcol

        type (grid), intent(in) :: gr

        integer, intent(in) :: solve_type

        ! For "over-implicit" weighted time step.
        ! This vector holds output from the LHS (implicit) portion of a term at a
        ! given vertical level.  This output is weighted and applied to the RHS.
        ! This is used if the implicit portion of the term is "over-implicit", which
        ! means that the LHS contribution is given extra weight (>1) in order to
        ! increase numerical stability.  A weighted factor must then be applied to
        ! the RHS in order to balance the weight.
        real( kind = core_rknd ), dimension(ndiags3,ngrdcol,nzm), intent(in) :: &
          lhs_ta     ! LHS turbulent advection term

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          rhs_ta,           & ! RHS turbulent advection terms
          wp2,              & ! w'^2 (momentum levels)                      [m^2/s^2]
          wpthvp,           & ! w'th_v' (momentum levels)                     [K m/s]
          invrs_tau_C4_zm,  & ! Inverse time-scale for C4 terms on mom. levels  [1/s]
          invrs_tau_C14_zm, & ! Inverse time-scale for C14 terms on mom. levels [1/s]
          wpxap,            & ! w'x_a' (momentum levels)                    [m^2/s^2]
          wpxbp,            & ! w'x_b' (momentum levels)                    [m^2/s^2]
          xap2,             & ! x_a'^2 (momentum levels)                    [m^2/s^2]
          xbp2,             & ! x_b'^2 (momentum levels)                    [m^2/s^2]
          thv_ds_zm,        & ! Dry, base-state theta_v on momentum levs.         [K]
          lhs_splat_wp2       ! LHS coefficient of wp2 splatting term           [1/s]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          xam,              & ! x_am (thermodynamic levels)                     [m/s]
          xbm                 ! x_bm (thermodynamic levels)                     [m/s]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          lhs_dp1_C4,  & ! LHS dissipation term 1, for up2 vp2 using C14
          lhs_dp1_C14    ! LHS dissipation term 1, for up2 vp2 using C4

        real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
          C4,        & ! Model parameter C_4                         [-]
          C_uu_shr,  & ! Model parameter C_uu_shr                    [-]
          C_uu_buoy, & ! Model parameter C_uu_buoy                   [-]
          C14          ! Model parameter C_{14}                      [-]

        type(stats_type), intent(inout) :: &
          stats

        !------------------- Output Variables -------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) :: &
          rhs    ! Explicit contributions to x variance/covariance terms

        real( kind = core_rknd ), intent(in) :: &
          dt
    """
    l_sample = stats.l_sample
    if solve_type == xp2_xpyp_vp2:
        name_ta, name_tp, name_dp1 = "vp2_ta", "vp2_tp", "vp2_dp1"
        name_pr1, name_pr2, name_splat = "vp2_pr1", "vp2_pr2", "vp2_splat"
    elif solve_type == xp2_xpyp_up2:
        name_ta, name_tp, name_dp1 = "up2_ta", "up2_tp", "up2_dp1"
        name_pr1, name_pr2, name_splat = "up2_pr1", "up2_pr2", "up2_splat"
    else:
        name_ta = name_tp = name_dp1 = name_pr1 = name_pr2 = name_splat = ""

    grid_dir = int(gr.grid_dir_indx)
    interior = slice(1, nzm - 1)
    km1 = slice(1 - grid_dir, nzm - 1 - grid_dir)
    kp1 = slice(1 + grid_dir, nzm - 1 + grid_dir)

    rhs = jnp.zeros((ngrdcol, nzm))

    rhs_pr1 = term_pr1(nzm, ngrdcol, C4, C14, xbp2, wp2, invrs_tau_C4_zm, invrs_tau_C14_zm)
    rhs_pr2 = term_pr2(
        nzm, nzt, ngrdcol, gr, C_uu_shr, C_uu_buoy,
        thv_ds_zm, wpthvp, wpxap, wpxbp, xam, xbm,
    )
    rhs_term_tp = term_tp_rhs(nzm, nzt, ngrdcol, xam, xam, wpxap, wpxap, gr.invrs_dzm)

    rhs = rhs.at[:, interior].set(
        rhs_ta[:, interior]
        + one_half * lhs_splat_wp2[:, interior] * wp2[:, interior]
        + (
            (one - gamma_over_implicit_ts)
            * (
                -lhs_ta[1 - grid_dir, :, interior] * xap2[:, kp1]
                - lhs_ta[1, :, interior] * xap2[:, interior]
                - lhs_ta[1 + grid_dir, :, interior] * xap2[:, km1]
            )
            + (one - C_uu_shr)[:, None] * rhs_term_tp[:, interior]
            + rhs_pr1[:, interior]
            + (one - gamma_over_implicit_ts)
            * (-lhs_dp1_C4[:, interior] - lhs_dp1_C14[:, interior])
            * xap2[:, interior]
            + rhs_pr2[:, interior]
            + one / dt * xap2[:, interior]
        )
    )

    if l_sample:
        stats_pr1 = term_pr1(
            nzm, ngrdcol, C4, jnp.zeros((ngrdcol,)), xbp2, wp2,
            invrs_tau_C4_zm, invrs_tau_C14_zm,
        )
        stats_pr2 = term_pr1(
            nzm, ngrdcol, jnp.zeros((ngrdcol,)), C14, xbp2, wp2,
            invrs_tau_C4_zm, invrs_tau_C14_zm,
        )

        if bool(l_sample and stats.var_on_stats_list(name_ta)):
            stats_zm_begin = jnp.zeros((ngrdcol, nzm))
            stats_zm_mod = jnp.zeros((ngrdcol, nzm))
            stats_zm_begin = stats_zm_begin.at[:, interior].set(-rhs_ta[:, interior])
            stats_zm_mod = stats_zm_mod.at[:, interior].set(
                (one - gamma_over_implicit_ts)
                * (
                    -lhs_ta[1 - grid_dir, :, interior] * xap2[:, kp1]
                    - lhs_ta[1, :, interior] * xap2[:, interior]
                    - lhs_ta[1 + grid_dir, :, interior] * xap2[:, km1]
                )
            )
            stats = stats.begin_budget(name_ta, stats_zm_begin)
            stats = stats.update_budget(name_ta, stats_zm_mod)

        if bool(l_sample and stats.var_on_stats_list(name_pr1)):
            stats_zm_begin = jnp.zeros((ngrdcol, nzm))
            stats_zm_mod = jnp.zeros((ngrdcol, nzm))
            stats_zm_begin = stats_zm_begin.at[:, interior].set(-stats_pr1[:, interior])
            stats_zm_mod = stats_zm_mod.at[:, interior].set(
                (one - gamma_over_implicit_ts)
                * (-lhs_dp1_C4[:, interior] * xap2[:, interior])
            )
            stats = stats.begin_budget(name_pr1, stats_zm_begin)
            stats = stats.update_budget(name_pr1, stats_zm_mod)

        if bool(l_sample and stats.var_on_stats_list(name_dp1)):
            stats_zm_begin = jnp.zeros((ngrdcol, nzm))
            stats_zm_mod = jnp.zeros((ngrdcol, nzm))
            stats_zm_begin = stats_zm_begin.at[:, interior].set(-stats_pr2[:, interior])
            stats_zm_mod = stats_zm_mod.at[:, interior].set(
                (one - gamma_over_implicit_ts)
                * (-lhs_dp1_C14[:, interior] * xap2[:, interior])
            )
            stats = stats.begin_budget(name_dp1, stats_zm_begin)
            stats = stats.update_budget(name_dp1, stats_zm_mod)

        stats_tmp = jnp.zeros((ngrdcol, nzm))
        stats_tmp = stats_tmp.at[:, 1:nzm - 1].set(rhs_pr2[:, 1:nzm - 1])
        stats = stats.update(name_pr2, stats_tmp)

        stats_tmp = jnp.zeros((ngrdcol, nzm))
        stats_tmp = stats_tmp.at[:, 1:nzm - 1].set(
            (one - C_uu_shr)[:, None] * rhs_term_tp[:, 1:nzm - 1]
        )
        stats = stats.update(name_tp, stats_tmp)

        stats_tmp = jnp.zeros((ngrdcol, nzm))
        stats_tmp = stats_tmp.at[:, 1:nzm - 1].set(
            one_half * lhs_splat_wp2[:, 1:nzm - 1] * wp2[:, 1:nzm - 1]
        )
        stats = stats.update(name_splat, stats_tmp)

    rhs = rhs.at[:, gr.k_lb_zm].set(xap2[:, gr.k_lb_zm])
    rhs = rhs.at[:, gr.k_ub_zm].set(w_tol_sqd)
    return rhs, stats


@partial(jax.jit, static_argnames=("nzm", "nzt", "ngrdcol", "solve_type"))
def xp2_xpyp_rhs(
    nzm, nzt, ngrdcol, gr, solve_type, dt,
    wpxap, wpxbp,
    xam, xbm, xapxbp, xpyp_forcing,
    Cn, invrs_tau_zm, threshold,
    lhs_ta, rhs_ta, stats,
):
    """Explicit contributions to r_t'^2, th_l'^2, r_t'th_l', sclr'r_t',
    sclr'th_l', or sclr'^2.

    Note:
      For LHS turbulent advection term:
         An "over-implicit" weighted time step is applied to this term.
         The weight of the implicit portion of this term is controlled by
         the factor gamma_over_implicit_ts (abbreviated "gamma" in the
         expression below).  A factor is added to the right-hand side of
         the equation in order to balance a weight that is not equal to 1,
         such that:
              -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
         where X is the variable that is being solved for in a predictive
         equation (x'^2 or x'y' in this case), y(t) is the linearized
         portion of the term that gets treated implicitly, and RHS is the
         portion of the term that is always treated explicitly.  A weight
         of greater than 1 can be applied to make the term more
         numerically stable.

    Returns:
        (rhs, stats)

    Fortran declarations (xp2_xpyp_rhs in src/CLUBB_core/advance_xp2_xpyp_module.F90):

        !------------------- Input Variables -------------------
        integer, intent(in) :: &
          nzm, &
          nzt, &
          ngrdcol

        type (grid), intent(in) :: gr

        integer, intent(in) :: solve_type

        real( kind = core_rknd ), intent(in) :: &
          dt                 ! Model timestep                              [s]

        ! For "over-implicit" weighted time step.
        ! This vector holds output from the LHS (implicit) portion of a term at a
        ! given vertical level.  This output is weighted and applied to the RHS.
        ! This is used if the implicit portion of the term is "over-implicit", which
        ! means that the LHS contribution is given extra weight (>1) in order to
        ! increase numerical stability.  A weighted factor must then be applied to
        ! the RHS in order to balance the weight.
        real( kind = core_rknd ), dimension(ndiags3,ngrdcol,nzm), intent(in) :: &
          lhs_ta    ! LHS turbulent advection (ta) term

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          rhs_ta,                  & ! RHS turbulent advection (ta) term
          wpxap,                   & ! w'x_a' (momentum levels)     [m/s{x_a units}]
          wpxbp,                   & ! w'x_b' (momentum levels)     [m/s{x_b units}]
          xapxbp,                  & ! x_a'x_b' (m-levs)          [{x_a un}{x_b un}]
          xpyp_forcing,            & ! <x'y'> forcing (m-levs)      [{x un}{x un}/s]
          invrs_tau_zm,            & ! Time-scale tau on momentum levels       [1/s]
          Cn                         ! Coefficient C_n                           [-]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          xam,                     & ! x_am (thermodynamic levels)     [{x_a units}]
          xbm                        ! x_bm (thermodynamic levels)     [{x_b units}]

        real( kind = core_rknd ), intent(in) :: &
          threshold    ! Smallest allowable mag. value for x_a'x_b'  [{x_am units}
                       !                                              *{x_bm units}]

        type(stats_type), intent(inout) :: &
          stats

        !------------------- Output Variables -------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) :: &
          rhs     ! Explicit contributions to x variance/covariance terms
    """
    l_sample = stats.l_sample
    if solve_type == xp2_xpyp_rtp2:
        name_ta, name_tp, name_tp1, name_tp2, name_dp1, name_f = (
            "rtp2_ta", "rtp2_tp", "", "", "rtp2_dp1", "rtp2_forcing"
        )
    elif solve_type == xp2_xpyp_thlp2:
        name_ta, name_tp, name_tp1, name_tp2, name_dp1, name_f = (
            "thlp2_ta", "thlp2_tp", "", "", "thlp2_dp1", "thlp2_forcing"
        )
    elif solve_type == xp2_xpyp_rtpthlp:
        name_ta, name_tp, name_tp1, name_tp2, name_dp1, name_f = (
            "rtpthlp_ta", "", "rtpthlp_tp1", "rtpthlp_tp2", "rtpthlp_dp1", "rtpthlp_forcing"
        )
    else:
        name_ta = name_tp = name_tp1 = name_tp2 = name_dp1 = name_f = ""

    grid_dir = int(gr.grid_dir_indx)
    interior = slice(1, nzm - 1)
    km1 = slice(1 - grid_dir, nzm - 1 - grid_dir)
    kp1 = slice(1 + grid_dir, nzm - 1 + grid_dir)

    rhs_term_tp = term_tp_rhs(nzm, nzt, ngrdcol, xam, xbm, wpxbp, wpxap, gr.invrs_dzm)
    rhs_term_dp1 = term_dp1_rhs(nzm, ngrdcol, Cn, invrs_tau_zm, threshold)
    lhs_term_dp1 = term_dp1_lhs(nzm, ngrdcol, gr, Cn, invrs_tau_zm)
    rhs = jnp.zeros((ngrdcol, nzm))

    rhs = rhs.at[:, interior].set(
        rhs_ta[:, interior]
        + (one - gamma_over_implicit_ts)
        * (
            -lhs_ta[1 - grid_dir, :, interior] * xapxbp[:, kp1]
            - lhs_ta[1, :, interior] * xapxbp[:, interior]
            - lhs_ta[1 + grid_dir, :, interior] * xapxbp[:, km1]
        )
        + rhs_term_tp[:, interior]
        + rhs_term_dp1[:, interior]
        + (one - gamma_over_implicit_ts)
        * (-lhs_term_dp1[:, interior] * xapxbp[:, interior])
        + xpyp_forcing[:, interior]
        + one / dt * xapxbp[:, interior]
    )

    if l_sample:
        xm_zeros = jnp.zeros((ngrdcol, nzt))
        wpxp_zeros = jnp.zeros((ngrdcol, nzm))
        stats_tp1 = term_tp_rhs(nzm, nzt, ngrdcol, xm_zeros, xbm, wpxp_zeros, wpxap, gr.invrs_dzm)
        stats_tp2 = term_tp_rhs(nzm, nzt, ngrdcol, xam, xm_zeros, wpxbp, wpxp_zeros, gr.invrs_dzm)

        l_ta_stats = bool(l_sample and stats.var_on_stats_list(name_ta))
        l_tp_stats = bool(l_sample and stats.var_on_stats_list(name_tp))
        l_tp1_stats = bool(l_sample and stats.var_on_stats_list(name_tp1))
        l_tp2_stats = bool(l_sample and stats.var_on_stats_list(name_tp2))
        l_dp1_stats = bool(l_sample and stats.var_on_stats_list(name_dp1))
        l_f_stats = bool(l_sample and stats.var_on_stats_list(name_f))

        if l_ta_stats:
            stats_ta_begin = jnp.zeros((ngrdcol, nzm))
            stats_ta_mod = jnp.zeros((ngrdcol, nzm))
            stats_ta_begin = stats_ta_begin.at[:, interior].set(-rhs_ta[:, interior])
            stats_ta_mod = stats_ta_mod.at[:, interior].set(
                (one - gamma_over_implicit_ts)
                * (
                    -lhs_ta[1 - grid_dir, :, interior] * xapxbp[:, kp1]
                    - lhs_ta[1, :, interior] * xapxbp[:, interior]
                    - lhs_ta[1 + grid_dir, :, interior] * xapxbp[:, km1]
                )
            )
            stats = stats.begin_budget(name_ta, stats_ta_begin)
            stats = stats.update_budget(name_ta, stats_ta_mod)

        if l_dp1_stats:
            stats_dp1_begin = jnp.zeros((ngrdcol, nzm))
            stats_dp1_mod = jnp.zeros((ngrdcol, nzm))
            stats_dp1_begin = stats_dp1_begin.at[:, interior].set(
                -rhs_term_dp1[:, interior]
            )
            stats_dp1_mod = stats_dp1_mod.at[:, interior].set(
                (one - gamma_over_implicit_ts)
                * (-lhs_term_dp1[:, interior] * xapxbp[:, interior])
            )
            stats = stats.begin_budget(name_dp1, stats_dp1_begin)
            stats = stats.update_budget(name_dp1, stats_dp1_mod)

        if l_tp_stats:
            stats_tmp = jnp.zeros((ngrdcol, nzm))
            stats_tmp = stats_tmp.at[:, 1:nzm - 1].set(rhs_term_tp[:, 1:nzm - 1])
            stats = stats.update(name_tp, stats_tmp)
        if l_tp1_stats:
            stats_tmp = jnp.zeros((ngrdcol, nzm))
            stats_tmp = stats_tmp.at[:, 1:nzm - 1].set(stats_tp1[:, 1:nzm - 1])
            stats = stats.update(name_tp1, stats_tmp)
        if l_tp2_stats:
            stats_tmp = jnp.zeros((ngrdcol, nzm))
            stats_tmp = stats_tmp.at[:, 1:nzm - 1].set(stats_tp2[:, 1:nzm - 1])
            stats = stats.update(name_tp2, stats_tmp)
        if l_f_stats:
            stats_tmp = jnp.zeros((ngrdcol, nzm))
            stats_tmp = stats_tmp.at[:, 1:nzm - 1].set(xpyp_forcing[:, 1:nzm - 1])
            stats = stats.update(name_f, stats_tmp)

    rhs = rhs.at[:, gr.k_lb_zm].set(xapxbp[:, gr.k_lb_zm])
    rhs = rhs.at[:, gr.k_ub_zm].set(threshold)
    return rhs, stats


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "nzt",
        "ngrdcol",
        "sclr_dim",
        "l_scalar_calc",
        "iiPDF_type",
        "l_upwind_xpyp_ta",
        "l_godunov_upwind_xpyp_ta",
    ),
)
def calc_xp2_xpyp_ta_terms(
    nzm, nzt, ngrdcol, sclr_dim, gr,
    wprtp, wprtp2, wpthlp, wpthlp2, wprtpthlp,
    rtp2, thlp2, rtpthlp, upwp, vpwp, up2, vp2, wp2,
    wp2_zt, wpsclrp, wpsclrp2, wpsclrprtp, wpsclrpthlp,
    sclrp2, sclrprtp, sclrpthlp,
    rho_ds_zt, invrs_rho_ds_zm, rho_ds_zm,
    wp3_on_wp2, wp3_on_wp2_zt, sigma_sqd_w,
    pdf_implicit_coefs_terms, l_scalar_calc,
    beta, iiPDF_type, l_upwind_xpyp_ta,
    l_godunov_upwind_xpyp_ta, stats,
):
    """This procedure calculates all the turbulent advection terms needed by the
    various LHS and RHS matrices. In general, first the implicit coeficients are
    calculated and used to calculate the LHS turbulent advection terms, then the
    explicit terms are calculated and used to calculate the RHS turbulent advection
    terms.

    Returns:
        (lhs_ta_wprtp2, lhs_ta_wpthlp2, lhs_ta_wprtpthlp, lhs_ta_wpup2,
          lhs_ta_wpvp2, lhs_ta_wpsclrp2, lhs_ta_wprtpsclrp, lhs_ta_wpthlpsclrp,
          rhs_ta_wprtp2, rhs_ta_wpthlp2, rhs_ta_wprtpthlp, rhs_ta_wpup2,
          rhs_ta_wpvp2, rhs_ta_wpsclrp2, rhs_ta_wprtpsclrp, rhs_ta_wpthlpsclrp,
          stats)

    Fortran declarations (calc_xp2_xpyp_ta_terms in src/CLUBB_core/advance_xp2_xpyp_module.F90):

        !------------------- Input Variables -------------------
        integer, intent(in) :: &
          nzm,      & ! Number of momentum vertical levels
          nzt,      & ! Number of thermodynamic vertical levels
          ngrdcol,  & ! Number of grid columns
          sclr_dim    ! Number of passive scalars

        type (grid), intent(in) :: &
          gr

        type(implicit_coefs_terms), intent(in) :: &
          pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

        real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim), intent(in) :: &
          wpsclrp         ! <w'sclr'> (momentum levels)        [m/s{sclr units}]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim), intent(in) :: &
          wpsclrp2,     & ! <w'sclr'^2> (thermodynamic levels) [m/s{sclr units}^2]
          wpsclrprtp,   & ! <w'sclr'r_t'> (thermo. levels)     [m/s{sclr units)kg/kg]
          wpsclrpthlp     ! <w'sclr'th_l'> (thermo. levels)    [m/s{sclr units}K]

        ! Passive scalar output
        real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim), intent(in) :: &
          sclrp2,       &
          sclrprtp,     &
          sclrpthlp

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          wp2,              & ! <w'^2> (momentum levels)              [m^2/s^2]
          wprtp,            & ! <w'r_t'> (momentum levels)            [(m/s)(kg/kg)]
          wpthlp,           & ! <w'th_l'> (momentum levels)           [(m K)/s]
          rtp2,             & ! <r_t'^2>                              [(kg/kg)^2]
          thlp2,            & ! <th_l'^2>                             [K^2]
          rtpthlp,          & ! <r_t'th_l'>                           [(kg K)/kg]
          upwp,             & ! <u'w'> (momentum levels)              [m^2/s^2]
          vpwp,             & ! <v'w'> (momentum levels)              [m^2/s^2]
          up2,              & ! <u'^2> (momentum levels)              [m^s/s^2]
          vp2,              & ! <v'^2> (momentum levels)              [m^2/s^2]
          invrs_rho_ds_zm,  & ! Inv. dry, static density @ mom. levs. [m^3/kg]
          rho_ds_zm,        & ! Dry, static density on momentum levs. [kg/m^3]
          wp3_on_wp2,       & ! Smoothed version of <w'^3>/<w'^2> zm  [m/s]
          sigma_sqd_w         ! sigma_sqd_w (momentum levels)         [-]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          wp2_zt,           & ! <w'^2> interpolated to thermo. levels [m^2/s^2]
          wprtp2,           & ! <w'r_t'^2> (thermodynamic levels)     [m/s (kg/kg)^2]
          wpthlp2,          & ! <w'th_l'^2> (thermodynamic levels)    [m/s K^2]
          wprtpthlp,        & ! <w'r_t'th_l'> (thermodynamic levels)  [m/s (kg/kg) K]
          rho_ds_zt,        & ! Dry, static density on thermo. levels [kg/m^3]
          wp3_on_wp2_zt       ! Smoothed version of <w'^3>/<w'^2> zt  [m/s]

        logical, intent(in) :: &
          l_scalar_calc

        real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
          beta    ! CLUBB tunable parameter beta

        integer, intent(in) :: &
          iiPDF_type    ! Selected option for the two-component normal (double
                        ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                        ! w, chi, and eta) portion of CLUBB's multivariate,
                        ! two-component PDF.

        logical, intent(in) :: &
          l_upwind_xpyp_ta, & ! This flag determines whether we want to use an upwind differencing
                              ! approximation rather than a centered differencing for turbulent or
                              ! mean advection terms. It affects rtp2, thlp2, up2, vp2, sclrp2,
                              ! rtpthlp, sclrprtp, & sclrpthlp.
          l_godunov_upwind_xpyp_ta ! This flag determines whether we want to use a Godunov-like upwind
                                   ! approximation rather than a centered differencing for  turbulent
                                   ! or mean advection terms. It affects rtp2, thlp2, up2, vp2, sclrp2,
                                   ! rtpthlp, sclrprtp, & sclrpthlp.

        type(stats_type), intent(inout) :: &
              stats

        !------------------- Output Variables -------------------

        ! Implicit (LHS) turbulent advection terms
        real( kind = core_rknd ), dimension(ndiags3,ngrdcol,nzm), intent(out) :: &
          lhs_ta_wprtp2,    & ! For <w'rt'^2>
          lhs_ta_wpthlp2,   & ! For <w'thl'^2>
          lhs_ta_wprtpthlp, & ! For <w'rt'thl'>
          lhs_ta_wpup2,     & ! For <w'u'^2>
          lhs_ta_wpvp2        ! For <w'v'^2>

        ! Explicit (RHS) turbulent advection terms
        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) :: &
          rhs_ta_wprtp2,    & ! For <w'rt'^2>
          rhs_ta_wpthlp2,   & ! For <w'thl'^2>
          rhs_ta_wprtpthlp, & ! For <w'rt'thl'>
          rhs_ta_wpup2,     & ! For <w'u'^2>
          rhs_ta_wpvp2        ! For <w'v'^2>

        ! Implicit (LHS) turbulent advection terms for scalars
        real( kind = core_rknd ), dimension(ndiags3,ngrdcol,nzm,sclr_dim), intent(out) :: &
          lhs_ta_wpsclrp2,    & ! For <w'sclr'^2>
          lhs_ta_wprtpsclrp,  & ! For <w'rt'sclr'>
          lhs_ta_wpthlpsclrp    ! For <w'thl'sclr'>

        ! Explicit (RHS) turbulent advection terms for scalars
        real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim), intent(out) :: &
          rhs_ta_wpsclrp2,    & ! For <w'sclr'^2>
          rhs_ta_wprtpsclrp,  & ! For <w'rt'sclr'>
          rhs_ta_wpthlpsclrp    ! For <w'thl'sclr'>
    """
    l_sample = stats.l_sample
    lhs_ta_wprtp2 = jnp.zeros((ndiags3, ngrdcol, nzm))
    lhs_ta_wpthlp2 = jnp.zeros((ndiags3, ngrdcol, nzm))
    lhs_ta_wprtpthlp = jnp.zeros((ndiags3, ngrdcol, nzm))
    lhs_ta_wpup2 = jnp.zeros((ndiags3, ngrdcol, nzm))
    lhs_ta_wpvp2 = jnp.zeros((ndiags3, ngrdcol, nzm))
    rhs_ta_wprtp2 = jnp.zeros((ngrdcol, nzm))
    rhs_ta_wpthlp2 = jnp.zeros((ngrdcol, nzm))
    rhs_ta_wprtpthlp = jnp.zeros((ngrdcol, nzm))
    rhs_ta_wpup2 = jnp.zeros((ngrdcol, nzm))
    rhs_ta_wpvp2 = jnp.zeros((ngrdcol, nzm))

    lhs_ta_wpsclrp2 = jnp.zeros((ndiags3, ngrdcol, nzm, sclr_dim))
    lhs_ta_wprtpsclrp = jnp.zeros((ndiags3, ngrdcol, nzm, sclr_dim))
    lhs_ta_wpthlpsclrp = jnp.zeros((ndiags3, ngrdcol, nzm, sclr_dim))
    rhs_ta_wpsclrp2 = jnp.zeros((ngrdcol, nzm, sclr_dim))
    rhs_ta_wprtpsclrp = jnp.zeros((ngrdcol, nzm, sclr_dim))
    rhs_ta_wpthlpsclrp = jnp.zeros((ngrdcol, nzm, sclr_dim))

    coef_wprtp2_implicit = jnp.zeros((ngrdcol, nzt))
    term_wprtp2_explicit = jnp.zeros((ngrdcol, nzt))
    coef_wpthlp2_implicit = jnp.zeros((ngrdcol, nzt))
    term_wpthlp2_explicit = jnp.zeros((ngrdcol, nzt))
    coef_wprtpthlp_implicit = jnp.zeros((ngrdcol, nzt))
    term_wprtpthlp_explicit = jnp.zeros((ngrdcol, nzt))
    coef_wpup2_implicit = jnp.zeros((ngrdcol, nzt))
    term_wpup2_explicit = jnp.zeros((ngrdcol, nzt))
    coef_wpvp2_implicit = jnp.zeros((ngrdcol, nzt))
    term_wpvp2_explicit = jnp.zeros((ngrdcol, nzt))

    coef_wprtp2_implicit_zm = jnp.zeros((ngrdcol, nzm))
    term_wprtp2_explicit_zm = jnp.zeros((ngrdcol, nzm))
    coef_wpthlp2_implicit_zm = jnp.zeros((ngrdcol, nzm))
    term_wpthlp2_explicit_zm = jnp.zeros((ngrdcol, nzm))
    coef_wprtpthlp_implicit_zm = jnp.zeros((ngrdcol, nzm))
    term_wprtpthlp_explicit_zm = jnp.zeros((ngrdcol, nzm))
    coef_wpup2_implicit_zm = jnp.zeros((ngrdcol, nzm))
    term_wpup2_explicit_zm = jnp.zeros((ngrdcol, nzm))
    coef_wpvp2_implicit_zm = jnp.zeros((ngrdcol, nzm))
    term_wpvp2_explicit_zm = jnp.zeros((ngrdcol, nzm))

    sgn_t_vel_rtp2 = jnp.zeros((ngrdcol, nzm))
    sgn_t_vel_thlp2 = jnp.zeros((ngrdcol, nzm))
    sgn_t_vel_rtpthlp = jnp.zeros((ngrdcol, nzm))
    sgn_t_vel_up2 = jnp.zeros((ngrdcol, nzm))
    sgn_t_vel_vp2 = jnp.zeros((ngrdcol, nzm))
    sgn_t_vel_rtp2_zt = jnp.zeros((ngrdcol, nzt))
    sgn_t_vel_thlp2_zt = jnp.zeros((ngrdcol, nzt))
    sgn_t_vel_rtpthlp_zt = jnp.zeros((ngrdcol, nzt))

    a1_coef = one / (one - sigma_sqd_w)
    a1_coef_zt = jnp.maximum(zm2zt(nzm, nzt, ngrdcol, gr, a1_coef), zero_threshold)
    wp_coef = (one - one_third * beta[:, None]) * a1_coef ** 2 * wp3_on_wp2 / wp2
    wp_coef_zt = (
        (one - one_third * beta[:, None])
        * a1_coef_zt ** 2 * wp3_on_wp2_zt / wp2_zt
    )

    if l_explicit_turbulent_adv_xpyp:
        if l_sample:
            coef_wprtp2_implicit = jnp.zeros_like(coef_wprtp2_implicit)
            coef_wpthlp2_implicit = jnp.zeros_like(coef_wpthlp2_implicit)
            coef_wprtpthlp_implicit = jnp.zeros_like(coef_wprtpthlp_implicit)

        if (not l_upwind_xpyp_ta) or l_sample:
            term_wprtp2_explicit = wprtp2
            term_wpthlp2_explicit = wpthlp2
            term_wprtpthlp_explicit = wprtpthlp

        if l_upwind_xpyp_ta:
            term_wprtp2_explicit_zm = zt2zm(nzm, nzt, ngrdcol, gr, wprtp2)
            sgn_t_vel_rtp2 = jnp.where(term_wprtp2_explicit_zm * rtp2 < zero, -one, one)
        rhs_ta_wprtp2 = xpyp_term_ta_pdf_rhs(
            nzm, nzt, ngrdcol, gr, term_wprtp2_explicit,
            rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
            l_upwind_xpyp_ta, sgn_t_vel_rtp2, term_wprtp2_explicit_zm,
        )

        if l_upwind_xpyp_ta:
            term_wpthlp2_explicit_zm = zt2zm(nzm, nzt, ngrdcol, gr, wpthlp2)
            sgn_t_vel_thlp2 = jnp.where(term_wpthlp2_explicit_zm * thlp2 < zero, -one, one)
        rhs_ta_wpthlp2 = xpyp_term_ta_pdf_rhs(
            nzm, nzt, ngrdcol, gr, term_wpthlp2_explicit,
            rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
            l_upwind_xpyp_ta, sgn_t_vel_thlp2, term_wpthlp2_explicit_zm,
        )

        if l_upwind_xpyp_ta:
            term_wprtpthlp_explicit_zm = zt2zm(nzm, nzt, ngrdcol, gr, wprtpthlp)
            sgn_t_vel_rtpthlp = jnp.where(term_wprtpthlp_explicit_zm * rtpthlp < zero, -one, one)
        rhs_ta_wprtpthlp = xpyp_term_ta_pdf_rhs(
            nzm, nzt, ngrdcol, gr, term_wprtpthlp_explicit,
            rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
            l_upwind_xpyp_ta, sgn_t_vel_rtpthlp, term_wprtpthlp_explicit_zm,
        )

        if l_scalar_calc:
            for sclr in range(sclr_dim):
                term_s = jnp.zeros((ngrdcol, nzt))
                term_s_zm = jnp.zeros((ngrdcol, nzm))
                sgn = jnp.zeros((ngrdcol, nzm))
                if l_upwind_xpyp_ta:
                    term_s_zm = zt2zm(nzm, nzt, ngrdcol, gr, wpsclrp2[:, :, sclr])
                    sgn = jnp.where(term_s_zm * sclrp2[:, :, sclr] < zero, -one, one)
                else:
                    term_s = wpsclrp2[:, :, sclr]
                rhs_s = xpyp_term_ta_pdf_rhs(
                    nzm, nzt, ngrdcol, gr, term_s,
                    rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                    l_upwind_xpyp_ta, sgn, term_s_zm,
                )
                rhs_ta_wpsclrp2 = rhs_ta_wpsclrp2.at[:, :, sclr].set(rhs_s)

                term_s = jnp.zeros((ngrdcol, nzt))
                term_s_zm = jnp.zeros((ngrdcol, nzm))
                if l_upwind_xpyp_ta:
                    term_s_zm = zt2zm(nzm, nzt, ngrdcol, gr, wpsclrprtp[:, :, sclr])
                    sgn = jnp.where(term_s_zm * sclrprtp[:, :, sclr] < zero, -one, one)
                else:
                    term_s = wpsclrprtp[:, :, sclr]
                rhs_s = xpyp_term_ta_pdf_rhs(
                    nzm, nzt, ngrdcol, gr, term_s,
                    rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                    l_upwind_xpyp_ta, sgn, term_s_zm,
                )
                rhs_ta_wprtpsclrp = rhs_ta_wprtpsclrp.at[:, :, sclr].set(rhs_s)

                term_s = jnp.zeros((ngrdcol, nzt))
                term_s_zm = jnp.zeros((ngrdcol, nzm))
                if l_upwind_xpyp_ta:
                    term_s_zm = zt2zm(nzm, nzt, ngrdcol, gr, wpsclrpthlp[:, :, sclr])
                    sgn = jnp.where(term_s_zm * sclrpthlp[:, :, sclr] < zero, -one, one)
                else:
                    term_s = wpsclrpthlp[:, :, sclr]
                rhs_s = xpyp_term_ta_pdf_rhs(
                    nzm, nzt, ngrdcol, gr, term_s,
                    rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                    l_upwind_xpyp_ta, sgn, term_s_zm,
                )
                rhs_ta_wpthlpsclrp = rhs_ta_wpthlpsclrp.at[:, :, sclr].set(rhs_s)
    else:
        if iiPDF_type == iiPDF_ADG1:
            if (not l_upwind_xpyp_ta) or l_sample:
                coef_wprtp2_implicit = one_third * beta[:, None] * a1_coef_zt * wp3_on_wp2_zt
                coef_wpthlp2_implicit = coef_wprtp2_implicit
                coef_wprtpthlp_implicit = coef_wprtp2_implicit
            if l_upwind_xpyp_ta:
                coef_wprtp2_implicit_zm = one_third * beta[:, None] * a1_coef * wp3_on_wp2
                sgn_t_vel_rtp2 = jnp.where(wp3_on_wp2 < zero, -one, one)

            if not l_godunov_upwind_xpyp_ta:
                lhs_ta_wprtp2 = xpyp_term_ta_pdf_lhs(
                    nzm, nzt, ngrdcol, gr, coef_wprtp2_implicit,
                    rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                    l_upwind_xpyp_ta, sgn_t_vel_rtp2, coef_wprtp2_implicit_zm,
                )
            else:
                coef_wprtp2_implicit = one_third * beta[:, None] * a1_coef_zt * wp3_on_wp2_zt
                coef_wpthlp2_implicit = coef_wprtp2_implicit
                coef_wprtpthlp_implicit = coef_wprtp2_implicit
                lhs_ta_wprtp2 = xpyp_term_ta_pdf_lhs_godunov(
                    nzm, nzt, ngrdcol, gr, coef_wprtp2_implicit,
                    invrs_rho_ds_zm, rho_ds_zm,
                )
            lhs_ta_wpthlp2 = lhs_ta_wprtp2
            lhs_ta_wprtpthlp = lhs_ta_wprtp2
            if l_scalar_calc:
                lhs_ta_wpsclrp2 = jnp.broadcast_to(
                    lhs_ta_wprtp2[:, :, :, None],
                    (ndiags3, ngrdcol, nzm, sclr_dim),
                )
                lhs_ta_wprtpsclrp = lhs_ta_wpsclrp2
                lhs_ta_wpthlpsclrp = lhs_ta_wpsclrp2

            if (not l_upwind_xpyp_ta) or l_sample:
                wprtp_zt = zm2zt(nzm, nzt, ngrdcol, gr, wprtp)
                wpthlp_zt = zm2zt(nzm, nzt, ngrdcol, gr, wpthlp)
                term_wprtp2_explicit = wp_coef_zt * wprtp_zt ** 2
                term_wpthlp2_explicit = wp_coef_zt * wpthlp_zt ** 2
                term_wprtpthlp_explicit = wp_coef_zt * wprtp_zt * wpthlp_zt
            else:
                wprtp_zt = jnp.zeros((ngrdcol, nzt))
                wpthlp_zt = jnp.zeros((ngrdcol, nzt))

            if l_upwind_xpyp_ta:
                term_wprtp2_explicit_zm = wp_coef * wprtp ** 2
                sgn_t_vel_rtp2 = jnp.where(wp3_on_wp2 < zero, -one, one)
            if not l_godunov_upwind_xpyp_ta:
                rhs_ta_wprtp2 = xpyp_term_ta_pdf_rhs(
                    nzm, nzt, ngrdcol, gr, term_wprtp2_explicit,
                    rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                    l_upwind_xpyp_ta, sgn_t_vel_rtp2, term_wprtp2_explicit_zm,
                )
            else:
                term_wprtp2_explicit_zm = wprtp ** 2
                sgn_t_vel_rtp2_zt = jnp.where(wp_coef_zt < zero, -one, one)
                rhs_ta_wprtp2 = xpyp_term_ta_pdf_rhs_godunov(
                    nzm, nzt, ngrdcol, gr, term_wprtp2_explicit_zm,
                    invrs_rho_ds_zm, sgn_t_vel_rtp2_zt, rho_ds_zm,
                )

            if l_upwind_xpyp_ta:
                term_wpthlp2_explicit_zm = wp_coef * wpthlp ** 2
                sgn_t_vel_thlp2 = jnp.where(wp3_on_wp2 < zero, -one, one)
            if not l_godunov_upwind_xpyp_ta:
                rhs_ta_wpthlp2 = xpyp_term_ta_pdf_rhs(
                    nzm, nzt, ngrdcol, gr, term_wpthlp2_explicit,
                    rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                    l_upwind_xpyp_ta, sgn_t_vel_thlp2, term_wpthlp2_explicit_zm,
                )
            else:
                term_wpthlp2_explicit_zm = wpthlp ** 2
                sgn_t_vel_thlp2_zt = jnp.where(wp_coef_zt < zero, -one, one)
                rhs_ta_wpthlp2 = xpyp_term_ta_pdf_rhs_godunov(
                    nzm, nzt, ngrdcol, gr, term_wpthlp2_explicit_zm,
                    invrs_rho_ds_zm, sgn_t_vel_thlp2_zt, rho_ds_zm,
                )

            if l_upwind_xpyp_ta:
                term_wprtpthlp_explicit_zm = wp_coef * wprtp * wpthlp
                sgn_t_vel_rtpthlp = jnp.where(wp3_on_wp2 < zero, -one, one)
            if not l_godunov_upwind_xpyp_ta:
                rhs_ta_wprtpthlp = xpyp_term_ta_pdf_rhs(
                    nzm, nzt, ngrdcol, gr, term_wprtpthlp_explicit,
                    rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                    l_upwind_xpyp_ta, sgn_t_vel_rtpthlp, term_wprtpthlp_explicit_zm,
                )
            else:
                term_wprtpthlp_explicit_zm = wprtp * wpthlp
                sgn_t_vel_rtpthlp_zt = jnp.where(wp_coef_zt < zero, -one, one)
                rhs_ta_wprtpthlp = xpyp_term_ta_pdf_rhs_godunov(
                    nzm, nzt, ngrdcol, gr, term_wprtpthlp_explicit_zm,
                    invrs_rho_ds_zm, sgn_t_vel_rtpthlp_zt, rho_ds_zm,
                )

            if l_scalar_calc:
                if not l_upwind_xpyp_ta:
                    wpsclrp_zt = jnp.zeros((ngrdcol, nzt, sclr_dim))
                    for sclr in range(sclr_dim):
                        wpsclrp_s = zm2zt(
                            nzm, nzt, ngrdcol, gr, wpsclrp[:, :, sclr],
                        )
                        wpsclrp_zt = wpsclrp_zt.at[:, :, sclr].set(wpsclrp_s)
                else:
                    wpsclrp_zt = jnp.zeros((ngrdcol, nzt, sclr_dim))

                for sclr in range(sclr_dim):
                    if l_upwind_xpyp_ta:
                        term_s_zm = wp_coef * wpsclrp[:, :, sclr] ** 2
                        sgn = jnp.where(wp3_on_wp2 < zero, -one, one)
                        term_s = jnp.zeros((ngrdcol, nzt))
                    else:
                        term_s = wp_coef_zt * wpsclrp_zt[:, :, sclr] ** 2
                        term_s_zm = jnp.zeros((ngrdcol, nzm))
                        sgn = jnp.zeros((ngrdcol, nzm))
                    if not l_godunov_upwind_xpyp_ta:
                        rhs_s = xpyp_term_ta_pdf_rhs(
                            nzm, nzt, ngrdcol, gr, term_s,
                            rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                            l_upwind_xpyp_ta, sgn, term_s_zm,
                        )
                        rhs_ta_wpsclrp2 = rhs_ta_wpsclrp2.at[:, :, sclr].set(rhs_s)
                    else:
                        term_s_zm = wpsclrp[:, :, sclr] ** 2
                        sgn_zt = jnp.where(wp_coef_zt < zero, -one, one)
                        rhs_s = xpyp_term_ta_pdf_rhs_godunov(
                            nzm, nzt, ngrdcol, gr, term_s_zm,
                            invrs_rho_ds_zm, sgn_zt, rho_ds_zm,
                        )
                        rhs_ta_wpsclrp2 = rhs_ta_wpsclrp2.at[:, :, sclr].set(rhs_s)

                    if l_upwind_xpyp_ta:
                        term_s_zm = wp_coef * wpsclrp[:, :, sclr] * wprtp
                        sgn = jnp.where(wp3_on_wp2 < zero, -one, one)
                        term_s = jnp.zeros((ngrdcol, nzt))
                    else:
                        term_s = wp_coef_zt * wpsclrp_zt[:, :, sclr] * wprtp_zt
                        term_s_zm = jnp.zeros((ngrdcol, nzm))
                        sgn = jnp.zeros((ngrdcol, nzm))
                    if not l_godunov_upwind_xpyp_ta:
                        rhs_s = xpyp_term_ta_pdf_rhs(
                            nzm, nzt, ngrdcol, gr, term_s,
                            rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                            l_upwind_xpyp_ta, sgn, term_s_zm,
                        )
                        rhs_ta_wprtpsclrp = rhs_ta_wprtpsclrp.at[:, :, sclr].set(rhs_s)
                    else:
                        term_s_zm = wpsclrp[:, :, sclr] * wprtp
                        sgn_zt = jnp.where(wp_coef_zt < zero, -one, one)
                        rhs_s = xpyp_term_ta_pdf_rhs_godunov(
                            nzm, nzt, ngrdcol, gr, term_s_zm,
                            invrs_rho_ds_zm, sgn_zt, rho_ds_zm,
                        )
                        rhs_ta_wprtpsclrp = rhs_ta_wprtpsclrp.at[:, :, sclr].set(rhs_s)

                    if l_upwind_xpyp_ta:
                        term_s_zm = wp_coef * wpsclrp[:, :, sclr] * wpthlp
                        sgn = jnp.where(wp3_on_wp2 < zero, -one, one)
                        term_s = jnp.zeros((ngrdcol, nzt))
                    else:
                        term_s = wp_coef_zt * wpsclrp_zt[:, :, sclr] * wpthlp_zt
                        term_s_zm = jnp.zeros((ngrdcol, nzm))
                        sgn = jnp.zeros((ngrdcol, nzm))
                    if not l_godunov_upwind_xpyp_ta:
                        rhs_s = xpyp_term_ta_pdf_rhs(
                            nzm, nzt, ngrdcol, gr, term_s,
                            rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                            l_upwind_xpyp_ta, sgn, term_s_zm,
                        )
                        rhs_ta_wpthlpsclrp = rhs_ta_wpthlpsclrp.at[:, :, sclr].set(rhs_s)
                    else:
                        term_s_zm = wpsclrp[:, :, sclr] * wpthlp
                        sgn_zt = jnp.where(wp_coef_zt < zero, -one, one)
                        rhs_s = xpyp_term_ta_pdf_rhs_godunov(
                            nzm, nzt, ngrdcol, gr, term_s_zm,
                            invrs_rho_ds_zm, sgn_zt, rho_ds_zm,
                        )
                        rhs_ta_wpthlpsclrp = rhs_ta_wpthlpsclrp.at[:, :, sclr].set(rhs_s)
        elif iiPDF_type == iiPDF_new:
            if (not l_upwind_xpyp_ta) or l_sample:
                coef_wprtp2_implicit = pdf_implicit_coefs_terms.coef_wprtp2_implicit
                coef_wpthlp2_implicit = pdf_implicit_coefs_terms.coef_wpthlp2_implicit
                coef_wprtpthlp_implicit = pdf_implicit_coefs_terms.coef_wprtpthlp_implicit
            if l_upwind_xpyp_ta:
                coef_wprtp2_implicit_zm = zt2zm(
                    nzm, nzt, ngrdcol, gr, pdf_implicit_coefs_terms.coef_wprtp2_implicit,
                )
                sgn_t_vel_rtp2 = jnp.where(coef_wprtp2_implicit_zm * rtp2 * rtp2 < zero, -one, one)
            lhs_ta_wprtp2 = xpyp_term_ta_pdf_lhs(
                nzm, nzt, ngrdcol, gr, coef_wprtp2_implicit,
                rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                l_upwind_xpyp_ta, sgn_t_vel_rtp2, coef_wprtp2_implicit_zm,
            )
            if l_upwind_xpyp_ta:
                coef_wpthlp2_implicit_zm = zt2zm(
                    nzm, nzt, ngrdcol, gr, pdf_implicit_coefs_terms.coef_wpthlp2_implicit,
                )
                sgn_t_vel_thlp2 = jnp.where(coef_wpthlp2_implicit_zm * thlp2 * thlp2 < zero, -one, one)
            lhs_ta_wpthlp2 = xpyp_term_ta_pdf_lhs(
                nzm, nzt, ngrdcol, gr, coef_wpthlp2_implicit,
                rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                l_upwind_xpyp_ta, sgn_t_vel_thlp2, coef_wpthlp2_implicit_zm,
            )
            if l_upwind_xpyp_ta:
                coef_wprtpthlp_implicit_zm = zt2zm(
                    nzm, nzt, ngrdcol, gr, pdf_implicit_coefs_terms.coef_wprtpthlp_implicit,
                )
                term_wprtpthlp_explicit_zm = zt2zm(
                    nzm, nzt, ngrdcol, gr, pdf_implicit_coefs_terms.term_wprtpthlp_explicit,
                )
                sgn_t_vel_rtpthlp = jnp.where(
                    (
                        (coef_wprtpthlp_implicit_zm * rtpthlp + term_wprtpthlp_explicit_zm)
                        * rtpthlp
                    ) < zero,
                    -one,
                    one,
                )
            lhs_ta_wprtpthlp = xpyp_term_ta_pdf_lhs(
                nzm, nzt, ngrdcol, gr, coef_wprtpthlp_implicit,
                rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                l_upwind_xpyp_ta, sgn_t_vel_rtpthlp, coef_wprtpthlp_implicit_zm,
            )
            if (not l_upwind_xpyp_ta) or l_sample:
                term_wprtp2_explicit = jnp.zeros_like(term_wprtp2_explicit)
                term_wpthlp2_explicit = jnp.zeros_like(term_wpthlp2_explicit)
                term_wprtpthlp_explicit = pdf_implicit_coefs_terms.term_wprtpthlp_explicit

            rhs_ta_wprtpthlp = xpyp_term_ta_pdf_rhs(
                nzm, nzt, ngrdcol, gr, term_wprtpthlp_explicit,
                rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                l_upwind_xpyp_ta, sgn_t_vel_rtpthlp, term_wprtpthlp_explicit_zm,
            )
            rhs_ta_wprtp2 = jnp.zeros_like(rhs_ta_wprtp2)
            rhs_ta_wpthlp2 = jnp.zeros_like(rhs_ta_wpthlp2)
        elif iiPDF_type == iiPDF_new_hybrid:
            if (not l_upwind_xpyp_ta) or l_sample:
                coef_wprtp2_implicit = pdf_implicit_coefs_terms.coef_wprtp2_implicit
                term_wprtp2_explicit = pdf_implicit_coefs_terms.term_wprtp2_explicit

            if l_upwind_xpyp_ta:
                coef_wprtp2_implicit_zm = zt2zm(
                    nzm, nzt, ngrdcol, gr, pdf_implicit_coefs_terms.coef_wprtp2_implicit,
                )
                term_wprtp2_explicit_zm = zt2zm(
                    nzm, nzt, ngrdcol, gr, pdf_implicit_coefs_terms.term_wprtp2_explicit,
                )
                sgn_t_vel_rtp2 = jnp.where(
                    (coef_wprtp2_implicit_zm * rtp2 + term_wprtp2_explicit_zm) * rtp2 < zero,
                    -one,
                    one,
                )
            lhs_ta_wprtp2 = xpyp_term_ta_pdf_lhs(
                nzm, nzt, ngrdcol, gr, coef_wprtp2_implicit,
                rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                l_upwind_xpyp_ta, sgn_t_vel_rtp2, coef_wprtp2_implicit_zm,
            )
            rhs_ta_wprtp2 = xpyp_term_ta_pdf_rhs(
                nzm, nzt, ngrdcol, gr, term_wprtp2_explicit,
                rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                l_upwind_xpyp_ta, sgn_t_vel_rtp2, term_wprtp2_explicit_zm,
            )

            if (not l_upwind_xpyp_ta) or l_sample:
                coef_wpthlp2_implicit = pdf_implicit_coefs_terms.coef_wpthlp2_implicit
                term_wpthlp2_explicit = pdf_implicit_coefs_terms.term_wpthlp2_explicit

            if l_upwind_xpyp_ta:
                coef_wpthlp2_implicit_zm = zt2zm(
                    nzm, nzt, ngrdcol, gr, pdf_implicit_coefs_terms.coef_wpthlp2_implicit,
                )
                term_wpthlp2_explicit_zm = zt2zm(
                    nzm, nzt, ngrdcol, gr, pdf_implicit_coefs_terms.term_wpthlp2_explicit,
                )
                sgn_t_vel_thlp2 = jnp.where(
                    (coef_wpthlp2_implicit_zm * thlp2 + term_wpthlp2_explicit_zm)
                    * thlp2 < zero,
                    -one,
                    one,
                )
            lhs_ta_wpthlp2 = xpyp_term_ta_pdf_lhs(
                nzm, nzt, ngrdcol, gr, coef_wpthlp2_implicit,
                rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                l_upwind_xpyp_ta, sgn_t_vel_thlp2, coef_wpthlp2_implicit_zm,
            )
            rhs_ta_wpthlp2 = xpyp_term_ta_pdf_rhs(
                nzm, nzt, ngrdcol, gr, term_wpthlp2_explicit,
                rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                l_upwind_xpyp_ta, sgn_t_vel_thlp2, term_wpthlp2_explicit_zm,
            )

            if (not l_upwind_xpyp_ta) or l_sample:
                coef_wprtpthlp_implicit = pdf_implicit_coefs_terms.coef_wprtpthlp_implicit
                term_wprtpthlp_explicit = pdf_implicit_coefs_terms.term_wprtpthlp_explicit

            if l_upwind_xpyp_ta:
                coef_wprtpthlp_implicit_zm = zt2zm(
                    nzm, nzt, ngrdcol, gr, pdf_implicit_coefs_terms.coef_wprtpthlp_implicit,
                )
                term_wprtpthlp_explicit_zm = zt2zm(
                    nzm, nzt, ngrdcol, gr, pdf_implicit_coefs_terms.term_wprtpthlp_explicit,
                )
                sgn_t_vel_rtpthlp = jnp.where(
                    (
                        (coef_wprtpthlp_implicit_zm * rtpthlp + term_wprtpthlp_explicit_zm)
                        * rtpthlp
                    ) < zero,
                    -one,
                    one,
                )
            lhs_ta_wprtpthlp = xpyp_term_ta_pdf_lhs(
                nzm, nzt, ngrdcol, gr, coef_wprtpthlp_implicit,
                rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                l_upwind_xpyp_ta, sgn_t_vel_rtpthlp, coef_wprtpthlp_implicit_zm,
            )
            rhs_ta_wprtpthlp = xpyp_term_ta_pdf_rhs(
                nzm, nzt, ngrdcol, gr, term_wprtpthlp_explicit,
                rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                l_upwind_xpyp_ta, sgn_t_vel_rtpthlp, term_wprtpthlp_explicit_zm,
            )

            if (not l_upwind_xpyp_ta) or l_sample:
                coef_wpup2_implicit = pdf_implicit_coefs_terms.coef_wpup2_implicit
                term_wpup2_explicit = pdf_implicit_coefs_terms.term_wpup2_explicit

            if l_upwind_xpyp_ta:
                coef_wpup2_implicit_zm = zt2zm(
                    nzm, nzt, ngrdcol, gr, pdf_implicit_coefs_terms.coef_wpup2_implicit,
                )
                term_wpup2_explicit_zm = zt2zm(
                    nzm, nzt, ngrdcol, gr, pdf_implicit_coefs_terms.term_wpup2_explicit,
                )
                sgn_t_vel_up2 = jnp.where(
                    (coef_wpup2_implicit_zm * up2 + term_wpup2_explicit_zm) * up2 < zero,
                    -one,
                    one,
                )
            lhs_ta_wpup2 = xpyp_term_ta_pdf_lhs(
                nzm, nzt, ngrdcol, gr, coef_wpup2_implicit,
                rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                l_upwind_xpyp_ta, sgn_t_vel_up2, coef_wpup2_implicit_zm,
            )
            rhs_ta_wpup2 = xpyp_term_ta_pdf_rhs(
                nzm, nzt, ngrdcol, gr, term_wpup2_explicit,
                rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                l_upwind_xpyp_ta, sgn_t_vel_up2, term_wpup2_explicit_zm,
            )

            if (not l_upwind_xpyp_ta) or l_sample:
                coef_wpvp2_implicit = pdf_implicit_coefs_terms.coef_wpvp2_implicit
                term_wpvp2_explicit = pdf_implicit_coefs_terms.term_wpvp2_explicit

            if l_upwind_xpyp_ta:
                coef_wpvp2_implicit_zm = zt2zm(
                    nzm, nzt, ngrdcol, gr, pdf_implicit_coefs_terms.coef_wpvp2_implicit,
                )
                term_wpvp2_explicit_zm = zt2zm(
                    nzm, nzt, ngrdcol, gr, pdf_implicit_coefs_terms.term_wpvp2_explicit,
                )
                sgn_t_vel_vp2 = jnp.where(
                    (coef_wpvp2_implicit_zm * vp2 + term_wpvp2_explicit_zm) * vp2 < zero,
                    -one,
                    one,
                )
            lhs_ta_wpvp2 = xpyp_term_ta_pdf_lhs(
                nzm, nzt, ngrdcol, gr, coef_wpvp2_implicit,
                rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                l_upwind_xpyp_ta, sgn_t_vel_vp2, coef_wpvp2_implicit_zm,
            )
            rhs_ta_wpvp2 = xpyp_term_ta_pdf_rhs(
                nzm, nzt, ngrdcol, gr, term_wpvp2_explicit,
                rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                l_upwind_xpyp_ta, sgn_t_vel_vp2, term_wpvp2_explicit_zm,
            )

            if l_scalar_calc:
                for sclr in range(sclr_dim):
                    coef_s = jnp.zeros((ngrdcol, nzt))
                    term_s = jnp.zeros((ngrdcol, nzt))
                    coef_s_zm = jnp.zeros((ngrdcol, nzm))
                    term_s_zm = jnp.zeros((ngrdcol, nzm))
                    sgn = jnp.zeros((ngrdcol, nzm))
                    if (not l_upwind_xpyp_ta) or l_sample:
                        coef_s = pdf_implicit_coefs_terms.coef_wpsclrp2_implicit[:, :, sclr]
                        term_s = pdf_implicit_coefs_terms.term_wpsclrp2_explicit[:, :, sclr]

                    if l_upwind_xpyp_ta:
                        coef_s_zm = zt2zm(
                            nzm, nzt, ngrdcol, gr,
                            pdf_implicit_coefs_terms.coef_wpsclrp2_implicit[:, :, sclr],
                        )
                        term_s_zm = zt2zm(
                            nzm, nzt, ngrdcol, gr,
                            pdf_implicit_coefs_terms.term_wpsclrp2_explicit[:, :, sclr],
                        )
                        sgn = jnp.where(
                            (coef_s_zm * sclrp2[:, :, sclr] + term_s_zm)
                            * sclrp2[:, :, sclr] < zero,
                            -one,
                            one,
                        )
                    lhs_s = xpyp_term_ta_pdf_lhs(
                        nzm, nzt, ngrdcol, gr, coef_s,
                        rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                        l_upwind_xpyp_ta, sgn, coef_s_zm,
                    )
                    rhs_s = xpyp_term_ta_pdf_rhs(
                        nzm, nzt, ngrdcol, gr, term_s,
                        rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                        l_upwind_xpyp_ta, sgn, term_s_zm,
                    )
                    lhs_ta_wpsclrp2 = lhs_ta_wpsclrp2.at[:, :, :, sclr].set(lhs_s)
                    rhs_ta_wpsclrp2 = rhs_ta_wpsclrp2.at[:, :, sclr].set(rhs_s)

                    coef_s = jnp.zeros((ngrdcol, nzt))
                    term_s = jnp.zeros((ngrdcol, nzt))
                    coef_s_zm = jnp.zeros((ngrdcol, nzm))
                    term_s_zm = jnp.zeros((ngrdcol, nzm))
                    sgn = jnp.zeros((ngrdcol, nzm))
                    if (not l_upwind_xpyp_ta) or l_sample:
                        coef_s = pdf_implicit_coefs_terms.coef_wprtpsclrp_implicit[:, :, sclr]
                        term_s = pdf_implicit_coefs_terms.term_wprtpsclrp_explicit[:, :, sclr]

                    if l_upwind_xpyp_ta:
                        coef_s_zm = zt2zm(
                            nzm, nzt, ngrdcol, gr,
                            pdf_implicit_coefs_terms.coef_wprtpsclrp_implicit[:, :, sclr],
                        )
                        term_s_zm = zt2zm(
                            nzm, nzt, ngrdcol, gr,
                            pdf_implicit_coefs_terms.term_wprtpsclrp_explicit[:, :, sclr],
                        )
                        sgn = jnp.where(
                            (coef_s_zm * sclrprtp[:, :, sclr] + term_s_zm)
                            * sclrprtp[:, :, sclr] < zero,
                            -one,
                            one,
                        )
                    lhs_s = xpyp_term_ta_pdf_lhs(
                        nzm, nzt, ngrdcol, gr, coef_s,
                        rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                        l_upwind_xpyp_ta, sgn, coef_s_zm,
                    )
                    rhs_s = xpyp_term_ta_pdf_rhs(
                        nzm, nzt, ngrdcol, gr, term_s,
                        rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                        l_upwind_xpyp_ta, sgn, term_s_zm,
                    )
                    lhs_ta_wprtpsclrp = lhs_ta_wprtpsclrp.at[:, :, :, sclr].set(lhs_s)
                    rhs_ta_wprtpsclrp = rhs_ta_wprtpsclrp.at[:, :, sclr].set(rhs_s)

                    coef_s = jnp.zeros((ngrdcol, nzt))
                    term_s = jnp.zeros((ngrdcol, nzt))
                    coef_s_zm = jnp.zeros((ngrdcol, nzm))
                    term_s_zm = jnp.zeros((ngrdcol, nzm))
                    sgn = jnp.zeros((ngrdcol, nzm))
                    if (not l_upwind_xpyp_ta) or l_sample:
                        coef_s = pdf_implicit_coefs_terms.coef_wpthlpsclrp_implicit[:, :, sclr]
                        term_s = pdf_implicit_coefs_terms.term_wpthlpsclrp_explicit[:, :, sclr]

                    if l_upwind_xpyp_ta:
                        coef_s_zm = zt2zm(
                            nzm, nzt, ngrdcol, gr,
                            pdf_implicit_coefs_terms.coef_wpthlpsclrp_implicit[:, :, sclr],
                        )
                        term_s_zm = zt2zm(
                            nzm, nzt, ngrdcol, gr,
                            pdf_implicit_coefs_terms.term_wpthlpsclrp_explicit[:, :, sclr],
                        )
                        sgn = jnp.where(
                            (coef_s_zm * sclrpthlp[:, :, sclr] + term_s_zm)
                            * sclrpthlp[:, :, sclr] < zero,
                            -one,
                            one,
                        )
                    lhs_s = xpyp_term_ta_pdf_lhs(
                        nzm, nzt, ngrdcol, gr, coef_s,
                        rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                        l_upwind_xpyp_ta, sgn, coef_s_zm,
                    )
                    rhs_s = xpyp_term_ta_pdf_rhs(
                        nzm, nzt, ngrdcol, gr, term_s,
                        rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
                        l_upwind_xpyp_ta, sgn, term_s_zm,
                    )
                    lhs_ta_wpthlpsclrp = lhs_ta_wpthlpsclrp.at[:, :, :, sclr].set(lhs_s)
                    rhs_ta_wpthlpsclrp = rhs_ta_wpthlpsclrp.at[:, :, sclr].set(rhs_s)

    if iiPDF_type != iiPDF_new_hybrid:
        if l_upwind_xpyp_ta:
            coef_wpup2_implicit_zm = one_third * beta[:, None] * a1_coef * wp3_on_wp2
            coef_wpvp2_implicit_zm = coef_wpup2_implicit_zm
            term_wpup2_explicit_zm = wp_coef * upwp ** 2
            term_wpvp2_explicit_zm = wp_coef * vpwp ** 2
            sgn_t_vel_up2 = jnp.where(wp3_on_wp2 < zero, -one, one)
            sgn_t_vel_vp2 = sgn_t_vel_up2
        else:
            upwp_zt = zm2zt(nzm, nzt, ngrdcol, gr, upwp)
            vpwp_zt = zm2zt(nzm, nzt, ngrdcol, gr, vpwp)
            coef_wpup2_implicit = one_third * beta[:, None] * a1_coef_zt * wp3_on_wp2_zt
            coef_wpvp2_implicit = coef_wpup2_implicit
            term_wpup2_explicit = wp_coef_zt * upwp_zt ** 2
            term_wpvp2_explicit = wp_coef_zt * vpwp_zt ** 2

        lhs_ta_wpup2 = xpyp_term_ta_pdf_lhs(
            nzm, nzt, ngrdcol, gr, coef_wpup2_implicit,
            rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
            l_upwind_xpyp_ta, sgn_t_vel_up2, coef_wpup2_implicit_zm,
        )
        lhs_ta_wpvp2 = lhs_ta_wpup2
        rhs_ta_wpup2 = xpyp_term_ta_pdf_rhs(
            nzm, nzt, ngrdcol, gr, term_wpup2_explicit,
            rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
            l_upwind_xpyp_ta, sgn_t_vel_up2, term_wpup2_explicit_zm,
        )
        rhs_ta_wpvp2 = xpyp_term_ta_pdf_rhs(
            nzm, nzt, ngrdcol, gr, term_wpvp2_explicit,
            rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,
            l_upwind_xpyp_ta, sgn_t_vel_vp2, term_wpvp2_explicit_zm,
        )

    if l_sample:
        stats = stats.update("coef_wprtp2_implicit", coef_wprtp2_implicit)
        stats = stats.update("term_wprtp2_explicit", term_wprtp2_explicit)
        stats = stats.update("coef_wpthlp2_implicit", coef_wpthlp2_implicit)
        stats = stats.update("term_wpthlp2_explicit", term_wpthlp2_explicit)
        stats = stats.update("coef_wprtpthlp_implicit", coef_wprtpthlp_implicit)
        stats = stats.update("term_wprtpthlp_explicit", term_wprtpthlp_explicit)

    return (
        lhs_ta_wprtp2, lhs_ta_wpthlp2, lhs_ta_wprtpthlp,
        lhs_ta_wpup2, lhs_ta_wpvp2, lhs_ta_wpsclrp2,
        lhs_ta_wprtpsclrp, lhs_ta_wpthlpsclrp,
        rhs_ta_wprtp2, rhs_ta_wpthlp2, rhs_ta_wprtpthlp,
        rhs_ta_wpup2, rhs_ta_wpvp2, rhs_ta_wpsclrp2,
        rhs_ta_wprtpsclrp, rhs_ta_wpthlpsclrp,
        stats,
    )


@partial(jax.jit, static_argnames=("nzm", "nzt", "ngrdcol"))
def term_tp_rhs(nzm, nzt, ngrdcol, xam, xbm, wpxbp, wpxap, invrs_dzm):
    """Turbulent production of x_a'x_b':  explicit portion of the code.

    The d(x_a'x_b')/dt equation contains a turbulent production term:

    - w'x_b' d(x_am)/dz - w'x_a' d(x_bm)/dz.

    This term is solved for completely explicitly and is discretized as
    follows:

    The values of w'x_a' and w'x_b' are found on the momentum levels, whereas
    the values of x_am and x_bm are found on the thermodynamic levels.  The
    derivatives of both x_am and x_bm are taken over the intermediate
    (central) momentum level.  All of the remaining mathematical operations
    take place at the central momentum level, yielding the desired result.

    ---------xam--------------xbm---------------------------- t(k)

    ===wpxap======d(xam)/dz=========d(xbm)/dz===wpxbp======== m(k)

    ---------xam--------------xbm---------------------------- t(k-1)

    The vertical indices t(k), m(k), and t(k-1) correspond with altitudes
    zt(k), zm(k), and zt(k-1), respectively.  The letter "t" is used for
    thermodynamic levels and the letter "m" is used for momentum levels.

    invrs_dzm(k) = 1 / ( zt(k) - zt(k-1) )

    Returns:
        rhs

    Fortran declarations (term_tp_rhs in src/CLUBB_core/advance_xp2_xpyp_module.F90):

        !------------------------ Input variables ------------------------
        integer, intent(in) :: &
          nzm, &
          nzt, &
          ngrdcol

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          xam,       & ! x_am(k)                     [{x_am units}]
          xbm          ! x_bm(k)                     [{x_bm units}]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          wpxbp,     & ! w'x_b'(k)                   [m/s {x_bm units}]
          wpxap,     & ! w'x_a'(k)                   [m/s {x_am units}]
          invrs_dzm    ! Inverse of grid spacing (k) [1/m]

        !------------------------ Return Variable ------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) :: &
          rhs
    """
    del nzt
    rhs = jnp.zeros((ngrdcol, nzm))
    interior = slice(1, nzm - 1)
    rhs = rhs.at[:, interior].set(
        -wpxbp[:, interior]
        * invrs_dzm[:, interior]
        * (xam[:, 1:] - xam[:, :-1])
        - wpxap[:, interior]
        * invrs_dzm[:, interior]
        * (xbm[:, 1:] - xbm[:, :-1])
    )
    return rhs


@partial(jax.jit, static_argnames=("nzm", "ngrdcol"))
def term_dp1_lhs(nzm, ngrdcol, gr, Cn, invrs_tau_zm):
    """Dissipation term 1 for x_a'x_b':  implicit portion of the code.

    The d(x_a'x_b')/dt equation contains dissipation term 1:

    - ( C_n / tau_zm ) x_a'x_b'.

    For cases where x_a'x_b' is a variance (in other words, where x_a and x_b
    are the same variable), the term is damped to a certain positive
    threshold, such that:

    - ( C_n / tau_zm ) * ( x_a'x_b' - threshold ).

    However, if x_a'x_b' is u'^2 or v'^2, damping to a minimum threshold value
    is part of pressure term 1 and is handled as part of function 'term_pr1'.
    Thus, for u'^2 and v'^2, function 'term_dp1_lhs' is called, but function
    'term_dp1_rhs' is not called, as function 'term_pr1' is called instead.

    For cases where x_a'x_b' is a covariance (in other words, where x_a and
    x_b are different variables), threshold is set to 0, and the expression
    reverts to the form found in the first equation.

    This term is broken into implicit and explicit portions.  The equations
    for u'^2, v'^2, and any covariances only include the implicit portion.
    The implicit portion of this term is:

    - ( C_n / tau_zm ) x_a'x_b'(t+1).

    Note:  When the implicit term is brought over to the left-hand side,
           the sign is reversed and the leading "-" in front of the term
           is changed to a "+".

    The timestep index (t+1) means that the value of x_a'x_b' being used is
    from the next timestep, which is being advanced to in solving the
    d(x_a'x_b')/dt equation.

    The values of x_a'x_b' are found on momentum levels.  The values of
    time-scale tau_zm are also found on momentum levels.

    Note:  For equations that use pressure term 1 (such as the equations for
           u'^2 and v'^2), C_n = ( 2*C_4 + C_14 ) / 3; which combines the
           implicit contributions for dissipation term 1 and pressure term 1
           into one expression.  Otherwise, C_n = C_2.

    Returns:
        lhs

    Fortran declarations (term_dp1_lhs in src/CLUBB_core/advance_xp2_xpyp_module.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nzm, &
          ngrdcol

        type (grid), intent(in) :: gr

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          Cn,          & ! Coefficient C_n                       [-]
          invrs_tau_zm   ! Inverse time-scale tau at momentum levels (k) [1/s]

        !--------------------------- Return Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
          lhs
    """
    lhs = jnp.zeros((ngrdcol, nzm))
    interior = slice(1, nzm - 1)
    # Momentum main diagonal: [ x xapxbp(k,<t+1>) ]
    lhs = lhs.at[:, interior].set(Cn[:, interior] * invrs_tau_zm[:, interior])
    # Set lower boundary to 0
    lhs = lhs.at[:, gr.k_lb_zm].set(zero)
    # Set upper boundary to 0
    lhs = lhs.at[:, gr.k_ub_zm].set(zero)
    return lhs


@partial(jax.jit, static_argnames=("nzm", "ngrdcol"))
def term_dp1_rhs(nzm, ngrdcol, Cn, invrs_tau_zm, threshold):
    """Dissipation term 1 for x_a'x_b':  explicit portion of the code.

    The d(x_a'x_b')/dt equation contains dissipation term 1:

    - ( C_n / tau_zm ) x_a'x_b'.

    For cases where x_a'x_b' is a variance (in other words, where x_a and x_b
    are the same variable), the term is damped to a certain positive
    threshold, such that:

    - ( C_n / tau_zm ) * ( x_a'x_b' - threshold ).

    However, if x_a'x_b' is u'^2 or v'^2, damping to a minimum threshold value
    is part of pressure term 1 and is handled as part of function 'term_pr1'.
    Thus, for u'^2 and v'^2, function 'term_dp1_lhs' is called, but function
    'term_dp1_rhs' is not called, as function 'term_pr1' is called instead.

    For cases where x_a'x_b' is a covariance (in other words, where x_a and
    x_b are different variables), threshold is set to 0, and the expression
    reverts to the form found in the first equation.

    This term is broken into implicit and explicit portions.  The equations
    for u'^2, v'^2, and any covariances only include the implicit portion.
    The explicit portion of this term is:

    + ( C_n / tau_zm ) * threshold.

    The values of time-scale tau_zm and the threshold are found on the
    momentum levels.

    Note:  The equations that use pressure term 1 (such as the equations for
           u'^2 and v'^2) do not call this function.  Thus, within this
           function, C_n = C_2.

    Returns:
        rhs

    Fortran declarations (term_dp1_rhs in src/CLUBB_core/advance_xp2_xpyp_module.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nzm, &
          ngrdcol

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          Cn,             & ! Coefficient C_n                               [-]
          invrs_tau_zm      ! Time-scale tau at momentum levels (k)         [1/s]

        real( kind = core_rknd ) :: &
          threshold         ! Minimum allowable magnitude value of x_a'x_b' [units vary]

        !--------------------------- Return Variable ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) :: &
          rhs
    """
    del ngrdcol
    rhs = Cn[:, :nzm] * invrs_tau_zm[:, :nzm] * threshold
    return rhs


@partial(jax.jit, static_argnames=("nzm", "ngrdcol"))
def term_pr1(nzm, ngrdcol, C4, C14, xbp2, wp2, invrs_tau_C4_zm, invrs_tau_C14_zm):
    """Pressure term 1 for x_a'x_b':  explicit portion of the code.

    Note:  Pressure term 1 is only used when x_a'x_b' is either u'^2 or v'^2.
           For the following description, pressure term 2 for u'^2 is used as
           the example.  Pressure term 2 for v'^2 is the same as pressure
           term 2 for u'^2, except that the v'^2 and u'^2 variables are
           switched.

    The d(u'^2)/dt equation contains dissipation term 1:

    - ( C_4 / tau_zm ) * ( u'^2 - (2/3)*em );

    where em = (1/2) * ( u'^2 + v'^2 + w'^2 );

    and with the substitution applied, dissipation term 1 becomes:

    - ( C_4 / tau_zm ) * ( u'^2 - (1/3) * ( u'^2 + v'^2 + w'^2 ) ).

    The d(u'^2)/dt equation also contains pressure term 1:

    - (2/3) * epsilon;

    where epsilon = C_14 * ( em / tau_zm ).

    Additionally, since pressure term 1 is a damping term, em is damped only
    to it's minimum threshold value, em_min, where:

    em_min = (1/2) * ( u'^2|_min + v'^2|_min + w'^2|_min )
         = (1/2) * ( w_tol^2 + w_tol^2 + w_tol^2 )
         = (3/2) * w_tol^2.

    With the damping threshold applied, epsilon becomes:

    epsilon = C_14 * ( ( em - em_min ) / tau_zm );

    and with all substitutions applied, pressure term 1 becomes:

    - (2/3) * ( C_14 / tau_zm )
            * [ (1/2) * ( u'^2 + v'^2 + w'^2 ) - (3/2) * w_tol^2 ].

    Dissipation term 1 and pressure term 1 are combined and simplify to:

    - [ ( 2*C_4 + C_14 ) / ( 3 * tau_zm ) ] * u'^2
       + [ ( C_4 - C_14 ) / ( 3 * tau_zm ) ] * ( v'^2 + w'^2 )
       + ( C_14 / tau_zm ) * w_tol^2.

    The combined term has both implicit and explicit components.
    The implicit component is:

    - [ ( 2*C_4 + C_14 ) / ( 3 * tau_zm ) ] * u'^2(t+1).

    Note:  When the implicit term is brought over to the left-hand side,
           the sign is reversed and the leading "-" in front of the term
           is changed to a "+".

    Timestep index (t) stands for the index of the current timestep, while
    timestep index (t+1) stands for the index of the next timestep, which is
    being advanced to in solving the d(x_a'x_b')/dt equation.

    The implicit component of the combined dp1 and pr1 term is solved in
    function "term_dp1_lhs" above, where "( 2*C_4 + C_14 ) / 3" is sent in
    as "C_n".

    The explicit component of the combined dp1 and pr1 term is:

    + [ ( C_4 - C_14 ) / ( 3 * tau_zm ) ] * ( v'^2(t) + w'^2(t) )
    + ( C_14 / tau_zm ) * w_tol^2;

    and is discretized as follows:

    The values for v'^2 and w'^2, as well as for tau_zm, are found on the
    momentum levels.  The mathematical operations all take place on the
    momentum levels, yielding the desired result.

    Returns:
        rhs

    Fortran declarations (term_pr1 in src/CLUBB_core/advance_xp2_xpyp_module.F90):

        !------------------------- Input Variables -------------------------
        integer, intent(in) :: &
          nzm, &
          ngrdcol

        real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
          C4,              & ! Model parameter C_4                                 [-]
          C14                ! Model parameter C_14                                [-]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          xbp2,            & ! v'^2(k) (if solving for u'^2) or vice versa         [m^2/s^2]
          wp2,             & ! w'^2(k)                                             [m^2/s^2]
          invrs_tau_C4_zm, & ! Time-scale tau for C4 terms at momentum levels (k)  [1/s]
          invrs_tau_C14_zm   ! Time-scale tau for C14 terms at momentum levels (k) [1/s]

        ! Return Variable
        !------------------------- Return Variables -------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) :: &
          rhs
    """
    rhs = jnp.zeros((ngrdcol, nzm))
    interior = slice(1, nzm - 1)
    rhs = rhs.at[:, interior].set(
        one_third
        * C4[:, None]
        * (xbp2[:, interior] + wp2[:, interior])
        * invrs_tau_C4_zm[:, interior]
        - one_third
        * C14[:, None]
        * (xbp2[:, interior] + wp2[:, interior])
        * invrs_tau_C14_zm[:, interior]
        + C14[:, None] * invrs_tau_C14_zm[:, interior] * w_tol_sqd
    )
    return rhs


@partial(jax.jit, static_argnames=("nzm", "nzt", "ngrdcol"))
def term_pr2(
    nzm, nzt, ngrdcol, gr, C_uu_shr, C_uu_buoy,
    thv_ds_zm, wpthvp, upwp, vpwp, um, vm,
):
    """Pressure term 2 for x_a'x_b':  explicit portion of the code.

    Note:  Pressure term 2 is only used when x_a'x_b' is either u'^2 or v'^2.
           For the following description, pressure term 2 for u'^2 is used as
           the example.  Pressure term 2 for v'^2 is the exact same as
           pressure term 2 for u'^2.

    The d(u'^2)/dt equation contains pressure term 2:

    + (2/3) C_5 [ (g/thv_ds) w'th_v' - u'w' du/dz - v'w' dv/dz ].

    Note that below we have broken up C5 into C_uu_shr for shear terms and
    C_uu_buoy for buoyancy terms.

    This term is solved for completely explicitly and is discretized as
    follows:

    The values of w'th_v', u'w', and v'w' are found on the momentum levels,
    whereas the values of um and vm are found on the thermodynamic levels.
    Additionally, the values of thv_ds_zm are found on the momentum levels.
    The derivatives of both um and vm are taken over the intermediate
    (central) momentum level.  All the remaining mathematical operations take
    place at the central momentum level, yielding the desired result.

    -----um--------------vm---------------------------------------- t(k)

    =upwp====d(um)/dz========d(vm)/dz==vpwp===thv_ds_zm==wpthvp==== m(k)

    -----um--------------vm---------------------------------------- t(k-1)

    The vertical indices t(k), m(k), and t(k-1) correspond with altitudes
    zt(k), zm(k), and zt(k-1), respectively.  The letter "t" is used for
    thermodynamic levels and the letter "m" is used for momentum levels.

    invrs_dzm(k) = 1 / ( zt(k) - zt(k-1) )

    Returns:
        rhs_pr2

    Fortran declarations (term_pr2 in src/CLUBB_core/advance_xp2_xpyp_module.F90):

        ! ------------ Input Variables ------------
        integer, intent(in) :: &
          nzm, &
          nzt, &
          ngrdcol

        type (grid), intent(in) :: gr

        real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
          C_uu_shr,  & ! Model parameter C_uu_shr                       [-]
          C_uu_buoy    ! Model parameter C_uu_buoy                      [-]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          thv_ds_zm, & ! Dry, base-state theta_v at momentum level (k)  [K]
          wpthvp,    & ! w'th_v'(k)                                     [m/K/s]
          upwp,      & ! u'w'(k)                                        [m^2/s^2]
          vpwp         ! v'w'(k)                                        [m^2/s^2]

        ! Note: Entire arrays of um and vm are now required rather than um and vm
        ! only at levels k and k+1.  The entire array is necessary when a vertical
        ! average calculation of d(um)/dz and d(vm)/dz is used. --ldgrant March 2010
        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          um,  & ! mean zonal wind       [m/s]
          vm     ! mean meridional wind  [m/s]

        ! ------------ Output Variable ------------
        real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
          rhs_pr2
    """
    del nzt
    rhs_pr2 = jnp.zeros((ngrdcol, nzm))
    interior = slice(1, nzm - 1)
    rhs_pr2_interior = two_thirds * (
        C_uu_buoy[:, None]
        * (grav / thv_ds_zm[:, interior])
        * wpthvp[:, interior]
        + C_uu_shr[:, None]
        * (
            -upwp[:, interior]
            * gr.invrs_dzm[:, interior]
            * (um[:, 1:] - um[:, :-1])
            - vpwp[:, interior]
            * gr.invrs_dzm[:, interior]
            * (vm[:, 1:] - vm[:, :-1])
        )
    )
    rhs_pr2 = rhs_pr2.at[:, interior].set(
        jnp.maximum(rhs_pr2_interior, zero_threshold)
    )
    # Added by dschanen for ticket #36
    # We have found that when shear generation is zero this term will only be
    # offset by hole-filling (up2_pd/vp2_pd) and reduces turbulence
    # unrealistically at lower altitudes to make up the difference.
    return rhs_pr2


@partial(
    jax.jit,
    static_argnames=("nzm", "ngrdcol", "solve_type", "fill_holes_type"),
)
def pos_definite_variances(
    nzm, ngrdcol, gr, solve_type, fill_holes_type, dt, tolerance, rho_ds_zm,
    stats, xp2_np1,
):
    """Use the hole filling code to make a variance term positive definite

    Returns:
        (xp2_np1, stats)

    Fortran declarations (pos_definite_variances in src/CLUBB_core/advance_xp2_xpyp_module.F90):

        !------------------- Input variables ------------------
        integer, intent(in) :: &
          nzm, &
          ngrdcol

        type (grid), intent(in) :: &
          gr

        integer, intent(in) :: &
          solve_type, &
          fill_holes_type

        real( kind = core_rknd ), intent(in) :: &
          dt        ! Model timestep              [s]

        real( kind = core_rknd ), intent(in) :: &
          tolerance ! Threshold for xp2_np1       [units vary]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          rho_ds_zm    ! Dry, static density on momentum levels         [kg/m^3]

        type(stats_type), intent(inout) :: &
          stats

        real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzm) ::  &
          xp2_np1   ! Variance for <n+1>          [units vary]
    """
    l_sample = stats.l_sample
    if solve_type == xp2_xpyp_rtp2:
        name_pd = "rtp2_pd"
    elif solve_type == xp2_xpyp_thlp2:
        name_pd = "thlp2_pd"
    elif solve_type == xp2_xpyp_up2:
        name_pd = "up2_pd"
    elif solve_type == xp2_xpyp_vp2:
        name_pd = "vp2_pd"
    else:
        name_pd = ""

    if l_sample and name_pd:
        # Store previous value for effect of the positive definite scheme
        stats = stats.begin_budget(name_pd, xp2_np1 / dt)

    if fill_holes_type != 0:
        # Call the hole-filling scheme.
        # The first pass-through should draw from only two levels on either side
        # of the hole.
        # upper_hf_level = nz-1 since we are filling the zm levels
        xp2_np1 = fill_holes_vertical(
            nzm, ngrdcol, tolerance,
            gr.k_lb_zm + gr.grid_dir_indx,
            gr.k_ub_zm - gr.grid_dir_indx,
            gr.dzm, rho_ds_zm, gr.grid_dir_indx,
            fill_holes_type, xp2_np1,
        )

    if l_sample and name_pd:
        stats = stats.finalize_budget(name_pd, xp2_np1 / dt)

    return xp2_np1, stats


def update_xp2_mc(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, dt: float, cloud_frac, rcm, rvm, thlm,
    wm, exner, rrm_evap, pdf_params: pdf_parameter,
    rtp2_mc, thlp2_mc, wprtp_mc, wpthlp_mc, rtpthlp_mc,
):
    """This subroutine is for use when l_morr_xp2_mc = .true.
    The effects of rain evaporation on rtp2 and thlp2 are included by
    assuming rain falls through the moist (cold) portion of the pdf.
    This is accomplished by defining a precip_fraction and assuming a double
    delta shaped pdf, such that the evaporation makes the moist component
    moister and the colder component colder. Calculations are done using
    variables on the zt grid, and the outputs are on the zm grid --storer

    Returns:
        (rtp2_mc, thlp2_mc, wprtp_mc, wpthlp_mc, rtpthlp_mc)

    Fortran declarations (update_xp2_mc in src/CLUBB_core/advance_xp2_xpyp_module.F90):

        ! -------------------- Input Variables --------------------
        integer, intent(in) :: &
          nzm,      & ! Momentum points in the Vertical        [-]
          nzt,      & ! Thermodynamic points in the Vertical        [-]
          ngrdcol

        type(grid), intent(in) :: gr

        real( kind = core_rknd ), intent(in) :: dt ! Model timestep        [s]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          cloud_frac, &       !Cloud fraction                        [-]
          rcm, &              !Cloud water mixing ratio              [kg/kg]
          rvm, &              !Vapor water mixing ratio              [kg/kg]
          thlm, &             !Liquid potential temperature          [K]
          wm, &               !Mean vertical velocity                [m/s]
          exner, &            !Exner function                        [-]
          rrm_evap         !Evaporation of rain                   [kg/kg/s]
                              !It is expected that this variable is negative, as
                              !that is the convention in Morrison microphysics

        type(pdf_parameter), intent(in) :: &
          pdf_params ! PDF parameters

        ! -------------------- Input/Output Variables --------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) :: &
          rtp2_mc, &    !Tendency of <rt'^2> due to evaporation   [(kg/kg)^2/s]
          thlp2_mc, &   !Tendency of <thl'^2> due to evaporation  [K^2/s]
          wprtp_mc, &   !Tendency of <w'rt'> due to evaporation   [m*(kg/kg)/s^2]
          wpthlp_mc, &  !Tendency of <w'thl'> due to evaporation  [m*K/s^2]
          rtpthlp_mc    !Tendency of <rt'thl'> due to evaporation [K*(kg/kg)/s]
    """
    cloud_frac = cloud_frac
    rcm = rcm
    rvm = rvm
    thlm = thlm
    wm = wm
    exner = exner
    rrm_evap = rrm_evap

    precip_frac_double_delta = jnp.zeros((ngrdcol, nzt))
    precip_frac_double_delta = precip_frac_double_delta.at[:, gr.k_ub_zt].set(zero)
    step = -int(gr.grid_dir_indx)
    start = int(gr.k_ub_zt) - int(gr.grid_dir_indx)
    stop = int(gr.k_lb_zt)
    stop_exclusive = stop + (1 if step > 0 else -1)
    ks = jnp.arange(start, stop_exclusive, step, dtype=jnp.int32)

    def precip_frac_step(next_col, k):
        col = jnp.where(
            cloud_frac[:, k] > cloud_frac_min,
            cloud_frac[:, k],
            next_col,
        )
        return col, col

    _, precip_cols = jax.lax.scan(
        precip_frac_step,
        jnp.zeros((ngrdcol,)),
        ks,
    )
    precip_frac_double_delta = precip_frac_double_delta.at[:, ks].set(
        jnp.swapaxes(precip_cols, 0, 1)
    )

    pf_const = jnp.where(
        precip_frac_double_delta > cloud_frac_min,
        (one - precip_frac_double_delta)
        / jnp.where(precip_frac_double_delta > cloud_frac_min, precip_frac_double_delta, one),
        zero,
    )

    temp_rtp2 = (
        pdf_params.mixt_frac * (
            (pdf_params.rt_1 - (rcm + rvm)) ** 2 + pdf_params.varnce_rt_1
        )
        + (one - pdf_params.mixt_frac) * (
            (pdf_params.rt_2 - (rcm + rvm)) ** 2 + pdf_params.varnce_rt_2
        )
    )
    rtp2_mc_zt = (
        rrm_evap ** 2 * pf_const * dt
        + two * jnp.abs(rrm_evap) * jnp.sqrt(temp_rtp2 * pf_const)
    )
    rtp2_mc = zt2zm(nzm, nzt, ngrdcol, gr, rtp2_mc_zt)

    temp_thlp2 = (
        pdf_params.mixt_frac * (
            (pdf_params.thl_1 - thlm) ** 2 + pdf_params.varnce_thl_1
        )
        + (one - pdf_params.mixt_frac) * (
            (pdf_params.thl_2 - thlm) ** 2 + pdf_params.varnce_thl_2
        )
    )
    thlp_fac = Lv / (Cp * exner)
    thlp2_mc_zt = (
        (rrm_evap * thlp_fac) ** 2 * pf_const * dt
        + two * jnp.abs(rrm_evap) * thlp_fac * jnp.sqrt(temp_thlp2 * pf_const)
    )
    thlp2_mc = zt2zm(nzm, nzt, ngrdcol, gr, thlp2_mc_zt)

    temp_wp2 = (
        pdf_params.mixt_frac * (
            (pdf_params.w_1 - wm) ** 2 + pdf_params.varnce_w_1
        )
        + (one - pdf_params.mixt_frac) * (
            (pdf_params.w_2 - wm) ** 2 + pdf_params.varnce_w_2
        )
    )
    sqrt_pf = jnp.sqrt(pf_const)
    wprtp_mc_zt = jnp.abs(rrm_evap) * sqrt_pf * jnp.sqrt(temp_wp2)
    wpthlp_mc_zt = -thlp_fac * jnp.abs(rrm_evap) * sqrt_pf * jnp.sqrt(temp_wp2)
    rtpthlp_mc_zt = (
        -jnp.abs(rrm_evap) * sqrt_pf * (thlp_fac * jnp.sqrt(temp_rtp2) + jnp.sqrt(temp_thlp2))
        - thlp_fac * pf_const * rrm_evap ** 2 * dt
    )

    wprtp_mc = zt2zm(nzm, nzt, ngrdcol, gr, wprtp_mc_zt)
    wpthlp_mc = zt2zm(nzm, nzt, ngrdcol, gr, wpthlp_mc_zt)
    rtpthlp_mc = zt2zm(nzm, nzt, ngrdcol, gr, rtpthlp_mc_zt)

    return rtp2_mc, thlp2_mc, wprtp_mc, wpthlp_mc, rtpthlp_mc
