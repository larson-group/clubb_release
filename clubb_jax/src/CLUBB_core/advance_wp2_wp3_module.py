"""JAX-side entry point for ``src/CLUBB_core/advance_wp2_wp3_module.F90``.

Description:
  Advance w'^2 and w'^3 one timestep.

References:
  https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:wp2_wp3_eqns

  Eqn. 12 & 18 on p. 3545--3546 of
  ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    Method and Model Description'' Golaz, et al. (2002)
    JAS, Vol. 59, pp. 3540--3551.

See also
  ``Equations for CLUBB'', Section 6:
  /Implict solution for the vertical velocity moments/

Porting deviations:
- Sponge damping blocks are unsupported here because clubb_case_initalization
  rejects sponge-enabled Python/JAX driver cases before this routine is called.
- The detailed Fortran error-print diagnostics and diagnostic-only early
  returns are reduced until full JAX diagnostic state is available; fatal
  conditions still mark err_info.
- Stats are threaded explicitly with JaxStats because the Fortran routine uses
  global stats side effects that are not JAX-compatible state.
- Intent(out) and intent(inout) variables are returned explicitly instead of
  being modified through Fortran argument side effects.
"""

from functools import partial

import jax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.Skx_module import Skx_func
from clubb_jax.src.CLUBB_core.advance_helper_module import calc_wp3_on_wp2
from clubb_jax.src.CLUBB_core.clip_explicit import clip_skewness, clip_variance
from clubb_jax.src.CLUBB_core.diffusion import diffusion_zt_lhs, diffusion_zm_lhs
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from clubb_jax.src.CLUBB_core.fill_holes import (
    fill_holes_vertical,
    fill_holes_wp2_from_horz_tke,
)
from clubb_jax.src.CLUBB_core.grid_class import (
    T_ABOVE,
    T_BELOW,
    ddzt,
    zm2zt,
    zm2zt2zm,
    zt2zm,
)
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.mean_adv import term_ma_zt_lhs, term_ma_zm_lhs
from clubb_jax.src.CLUBB_core.matrix_solver_wrapper import band_solve
from clubb_jax.src.CLUBB_core.constants_clubb import (
    eps,
    five,
    four,
    gamma_over_implicit_ts,
    grav,
    max_mag_correlation_flux,
    one,
    one_half,
    one_third,
    three,
    three_halves,
    two,
    two_thirds,
    w_tol,
    w_tol_sqd,
    wp2_max,
    zero,
    zero_threshold,
)
from clubb_jax.src.CLUBB_core.parameter_indices import (
    iC1,
    iC1b,
    iC1c,
    iC4,
    iC8,
    iC8b,
    iC11,
    iC11b,
    iC11c,
    iC12,
    iC_uu_buoy,
    iC_uu_shr,
    iC_wp2_pr_dfsn,
    iC_wp3_pr_dfsn,
    iC_wp3_pr_tp,
    iC_wp3_pr_turb,
    iSkw_max_mag,
    ia3_coef_min,
    ic_K1,
    ic_K8,
)
from clubb_jax.src.CLUBB_core.model_flags import (
    iiPDF_ADG1,
    iiPDF_new,
    iiPDF_new_hybrid,
    l_explicit_turbulent_adv_wp3,
    l_force_descending_solves,
    penta_bicgstab,
)
from clubb_jax.src.CLUBB_core import (
    ErrInfo,
    Grid,
    NuVertResDep,
    implicit_coefs_terms,
)


# Private named constants to avoid string comparisons
clip_wp2 = 12  # Named constant for wp2 clipping.
# NOTE: This must be the same as the clip_wp2 declared in clip_explicit!


# Set logical to true for Crank-Nicholson diffusion scheme
# or to false for completely implicit diffusion scheme.
# Note:  Although Crank-Nicholson diffusion has usually been used for wp2
#        and wp3 in the past, we found that using completely implicit
#        diffusion stabilized the deep convective cases more while having
#        almost no effect on the boundary layer cases.  Brian; 1/4/2008.
l_crank_nich_diff = False

ndiags2 = 2
ndiags3 = 3
ndiags5 = 5


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "nzt",
        "ngrdcol",
        "iiPDF_type",
        "penta_solve_method",
        "fill_holes_type",
        "l_min_wp2_from_corr_wx",
        "l_upwind_xm_ma",
        "l_standard_term_ta",
        "l_partial_upwind_wp3",
        "l_damp_wp2_using_em",
        "l_use_C11_Richardson",
        "l_damp_wp3_Skw_squared",
        "l_lmm_stepping",
        "l_use_tke_in_wp3_pr_turb_term",
        "l_use_tke_in_wp2_wp3_K_dfsn",
        "l_use_wp3_lim_with_smth_Heaviside",
        "l_wp2_fill_holes_tke",
        "l_ho_nontrad_coriolis",
        "l_implemented",
    ),
)
def advance_wp2_wp3(
    nzm: int, nzt: int, ngrdcol: int, gr: Grid, dt,
    sfc_elevation, fcor_y, sigma_sqd_w, wm_zm,
    wm_zt,
    wpup2, wpvp2, wp2up2, wp2vp2, wp4,
    wpthvp, wp2thvp, wp2up, um, vm, upwp, vpwp,
    em, Kh_zm, Kh_zt, invrs_tau_C4_zm,
    invrs_tau_wp3_zt, invrs_tau_C1_zm,
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,
    invrs_rho_ds_zt, thv_ds_zm,
    thv_ds_zt, mixt_frac, Cx_fnc_Richardson,
    lhs_splat_wp2, lhs_splat_wp3,
    pdf_implicit_coefs_terms: implicit_coefs_terms,
    wprtp, wpthlp, rtp2, thlp2,
    clubb_params, nu_vert_res_dep: NuVertResDep,
    iiPDF_type: int,
    penta_solve_method: int,
    fill_holes_type: int,
    l_min_wp2_from_corr_wx: bool,
    l_upwind_xm_ma: bool,
    l_standard_term_ta: bool,
    l_partial_upwind_wp3: bool,
    l_damp_wp2_using_em: bool,
    l_use_C11_Richardson: bool,
    l_damp_wp3_Skw_squared: bool,
    l_lmm_stepping: bool,
    l_use_tke_in_wp3_pr_turb_term: bool,
    l_use_tke_in_wp2_wp3_K_dfsn: bool,
    l_use_wp3_lim_with_smth_Heaviside: bool,
    l_wp2_fill_holes_tke: bool,
    l_ho_nontrad_coriolis: bool,
    l_implemented: bool,
    stats: JaxStats,
    up2, vp2, wp2, wp3, err_info: ErrInfo,
):
    """Advance w'^2 and w'^3 one timestep."""
    del mixt_frac

    wp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, wp2, w_tol_sqd)
    wp3_zm = zt2zm(nzm, nzt, ngrdcol, gr, wp3)
    Skw_zt = Skx_func(nzt, ngrdcol, wp2_zt, wp3, w_tol, clubb_params)
    Skw_zm = Skx_func(nzm, ngrdcol, wp2, wp3_zm, w_tol, clubb_params)
    wp3_on_wp2, wp3_on_wp2_zt = calc_wp3_on_wp2(
        nzm, nzt, ngrdcol, gr, wp2, wp3,
    )
    del wp3_on_wp2_zt

    # Compute the a3 coefficient (formula 25 in `Equations for CLUBB')
    a3_coef = -two * (one - sigma_sqd_w) ** 2 + three
    a3_coef = jnp.maximum(a3_coef, clubb_params[:, ia3_coef_min][:, None])
    a3_coef_zt = zm2zt(nzm, nzt, ngrdcol, gr, a3_coef)

    if stats.l_sample:
        stats = stats.update("a3_coef_zt", a3_coef_zt)
        stats = stats.update("a3_coef", a3_coef)

    if l_crank_nich_diff and l_use_tke_in_wp2_wp3_K_dfsn:
        # General error -> set all entries to clubb_fatal_error
        err_info = err_info.set_fatal()
        return up2, vp2, wp2, wp3, err_info, stats

    # Vince Larson added code to make C11 function of Skw. 13 Mar 2005
    # If this code is used, C11 is no longer relevant, i.e. constants
    #    are hardwired.
    if l_use_C11_Richardson:
        C11_Skw_fnc = zm2zt(
            nzm, nzt, ngrdcol, gr, Cx_fnc_Richardson, zero_threshold,
        )
    else:
        C11 = clubb_params[:, iC11]
        C11b = clubb_params[:, iC11b]
        C11c = clubb_params[:, iC11c]
        C11c_safe = jnp.where(jnp.abs(C11c) > zero, C11c, one)
        C11_varying = jnp.abs(C11 - C11b) > jnp.abs(C11 + C11b) * eps / two
        # Calculate C_{1} and C_{11} as functions of skewness of w.
        # The if..then here is only for computational efficiency -dschanen 2 Sept 08
        C11_Skw_fnc = jnp.where(
            C11_varying[:, None],
            C11b[:, None]
            + (C11 - C11b)[:, None]
            * jnp.exp(-one_half * (Skw_zt / C11c_safe[:, None]) ** 2),
            C11b[:, None],
        )

    C1 = clubb_params[:, iC1]
    C1b = clubb_params[:, iC1b]
    C1c = clubb_params[:, iC1c]
    C1c_safe = jnp.where(jnp.abs(C1c) > zero, C1c, one)
    C1_varying = jnp.abs(C1 - C1b) > jnp.abs(C1 + C1b) * eps / two
    # The if..then here is only for computational efficiency -dschanen 2 Sept 08
    C1_Skw_fnc = jnp.where(
        C1_varying[:, None],
        C1b[:, None]
        + (C1 - C1b)[:, None]
        * jnp.exp(-one_half * (Skw_zm / C1c_safe[:, None]) ** 2),
        C1b[:, None],
    )

    if l_damp_wp2_using_em:
        # Insert 1/3 here to account for the fact that in the dissipation term,
        #   (2/3)*em = (2/3)*(1/2)*(wp2+up2+vp2).  Then we can insert wp2, up2,
        #   and vp2 directly into the dissipation subroutines without prefixing them by (1/3).
        C1_Skw_fnc = one_third * C1_Skw_fnc

    # Set C16_fnc based on Richardson_num
    C16_fnc = Cx_fnc_Richardson[:, :nzt]

    if clubb_at_least_debug_level(0):
        # Assertion check for C11_Skw_fnc
        C11_bad = jnp.any((C11_Skw_fnc > one) | (C11_Skw_fnc < zero), axis=1)
        # Assertion check for C11_Skw_fnc
        C16_bad = jnp.any((C16_fnc > one) | (C16_fnc < zero), axis=1)
        # Error in grid column i -> set ith entry to clubb_fatal_error
        err_info = err_info.set_fatal(mask=C11_bad | C16_bad)

    if stats.l_sample:
        stats = stats.update("C11_Skw_fnc", C11_Skw_fnc)
        stats = stats.update("C1_Skw_fnc", C1_Skw_fnc)

    # Define the Coefficent of Eddy Diffusivity for the wp2 and wp3.
    # Kw1 is used for wp2, which is located on momentum levels.
    # Kw1 is located on thermodynamic levels.
    # Kw1 = c_K1 * Kh_zt ! c_K1 variable is commented out
    Kw1 = clubb_params[:, ic_K1][:, None] * Kh_zt
    # Kw8 is used for wp3, which is located on thermodynamic levels.
    # Kw8 is located on momentum levels.
    # Note: Kw8 is usually defined to be 1/2 of Kh_zm.
    # Kw8 = c_K8 * Kh_zm ! c_K8 variable is commented out
    Kw8 = clubb_params[:, ic_K8][:, None] * Kh_zm

    coef_wp4_implicit = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    a1_coef = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    a1_coef_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)

    if not l_explicit_turbulent_adv_wp3:
        if iiPDF_type == iiPDF_new or iiPDF_type == iiPDF_new_hybrid:
            # Unpack coef_wp4_implicit from pdf_implicit_coefs_terms.
            # Since PDF parameters and the resulting implicit coefficients and
            # explicit terms are calculated on thermodynamic levels, the <w'^4>
            # implicit coefficient needs to be unpacked as coef_wp4_implicit_zt.
            coef_wp4_implicit_zt = pdf_implicit_coefs_terms.coef_wp4_implicit
            # The values of <w'^4> are located on momentum levels.  Interpolate
            # coef_wp4_implicit_zt to momentum levels as coef_wp4_implicit.  The
            # discretization diagram is found in the description section of
            # function wp3_term_ta_new_pdf_lhs below.  These values are always
            # positive.
            coef_wp4_implicit = zt2zm(
                nzm, nzt, ngrdcol, gr, coef_wp4_implicit_zt, zero_threshold,
            )
            # Set the value of coef_wp4_implicit to 0 at the lower boundary and at
            # the upper boundary.  This sets the value of <w'^4> to 0 at the lower
            # and upper boundaries.
            coef_wp4_implicit = coef_wp4_implicit.at[:, gr.k_lb_zm].set(zero)
            coef_wp4_implicit = coef_wp4_implicit.at[:, gr.k_ub_zm].set(zero)
            if stats.l_sample:
                stats = stats.update("coef_wp4_implicit", coef_wp4_implicit)
        elif iiPDF_type == iiPDF_ADG1:
            # Define a_1 and a_3 (both are located on momentum levels).
            # They are variables that are both functions of sigma_sqd_w (where
            # sigma_sqd_w is located on momentum levels).
            a1_coef = one / (one - sigma_sqd_w)
            # Interpolate a_1 from momentum levels to thermodynamic levels.  This
            # will be used for the w'^3 turbulent advection (ta) term.
            a1_coef_zt = zm2zt(
                nzm, nzt, ngrdcol, gr, a1_coef, zero_threshold,
            )

    # Not using pressure term, set to 0
    rhs_pr3_wp3 = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)

    # Initiaize some terms to zero
    wp3_term_ta_lhs_result = jnp.zeros((ndiags5, ngrdcol, nzt), dtype=jnp.float64)
    wp3_pr3_lhs = jnp.zeros((ndiags5, ngrdcol, nzt), dtype=jnp.float64)

    Kw1_zm = zt2zm(nzm, nzt, ngrdcol, gr, Kw1, zero)
    Kw8_zt = zm2zt(nzm, nzt, ngrdcol, gr, Kw8, zero)

    # Experimental bouyancy term
    # Experimental term from CLUBB TRAC ticket #411
    # Compute the vertical derivative of the u and v winds
    if not l_use_tke_in_wp3_pr_turb_term:
        dum_dz = ddzt(nzm, nzt, ngrdcol, gr, um)
        dvm_dz = ddzt(nzm, nzt, ngrdcol, gr, vm)
    else:
        dum_dz = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        dvm_dz = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)

    em_smth = zm2zt2zm(nzm, nzt, ngrdcol, gr, em)
    wp2_smth = zm2zt2zm(nzm, nzt, ngrdcol, gr, wp2)

    rhs_pr_turb_wp3 = wp3_term_pr_turb_rhs(
        nzm, nzt, ngrdcol, gr, clubb_params[:, iC_wp3_pr_turb],
        Kh_zt, wpthvp,
        dum_dz, dvm_dz,
        upwp, vpwp,
        thv_ds_zt,
        rho_ds_zm, invrs_rho_ds_zt,
        em_smth, wp2_smth,
        l_use_tke_in_wp3_pr_turb_term,
    )

    # Calculate term
    rhs_pr_dfsn_wp3 = wp3_term_pr_dfsn_rhs(
        nzm, nzt, ngrdcol, gr,
        clubb_params[:, iC_wp3_pr_dfsn],
        rho_ds_zm, invrs_rho_ds_zt,
        wp2up2, wp2vp2, wp4,
        up2, vp2, wp2,
    )

    rhs_pr_dfsn_wp2 = wp2_term_pr_dfsn_rhs(
        nzm, nzt, ngrdcol, gr, clubb_params[:, iC_wp2_pr_dfsn],
        rho_ds_zt, invrs_rho_ds_zm,
        wpup2, wpvp2, wp3,
    )

    # This part handles the wp2 equation terms.
    lhs_diff_zm = diffusion_zm_lhs(
        nzm, nzt, ngrdcol, gr, Kw1, Kw1_zm, nu_vert_res_dep.nu1,
        invrs_rho_ds_zm, rho_ds_zt,
    )

    # This part handles the wp3 equation terms.
    lhs_diff_zt = diffusion_zt_lhs(
        nzm, nzt, ngrdcol, gr, Kw8, Kw8_zt, nu_vert_res_dep.nu8,
        invrs_rho_ds_zt, rho_ds_zm,
    )

    lhs_diff_zm_crank = jnp.zeros_like(lhs_diff_zm)
    lhs_diff_zt_crank = jnp.zeros_like(lhs_diff_zt)
    if l_crank_nich_diff:
        # Calculate RHS eddy diffusion terms for w'2 and w'3
        lhs_diff_zm_crank = lhs_diff_zm_crank.at[:, :, 1:-1].set(
            one_half * lhs_diff_zm[:, :, 1:-1]
        )
        lhs_diff_zt_crank = lhs_diff_zt_crank.at[:, :, 1:-1].set(
            one_half
            * lhs_diff_zt[:, :, 1:-1]
            * clubb_params[:, iC12][None, :, None]
        )

    # Calculate "over-implicit" pressure terms for w'2 and w'3
    rhs_pr1_wp2 = wp2_term_pr1_rhs(
        nzm, ngrdcol, gr, clubb_params[:, iC4],
        up2, vp2, invrs_tau_C4_zm,
    )

    # Note:  An "over-implicit" weighted time step is applied to the  term.
    #        A weighting factor of greater than 1 may be used to make the
    #        term more numerically stable (see note below for w'^3 RHS
    #        turbulent advection (ta) term).
    lhs_pr1_wp2 = wp2_term_pr1_lhs(
        nzm, ngrdcol, gr,
        clubb_params[:, iC4], invrs_tau_C4_zm,
    )
    # Calculate turbulent production terms of w'^3
    C_wp3_pr_tp = jnp.ones((ngrdcol,), dtype=jnp.float64)
    lhs_adv_tp_wp3 = wp3_term_tp_lhs(
        nzm, nzt, ngrdcol, gr, C_wp3_pr_tp,
        wp2, rho_ds_zm, invrs_rho_ds_zt,
    )

    # Calculate pressure damping of turbulent production of w'^3
    C_wp3_pr_tp = -clubb_params[:, iC_wp3_pr_tp]
    lhs_pr_tp_wp3 = wp3_term_tp_lhs(
        nzm, nzt, ngrdcol, gr, C_wp3_pr_tp,
        wp2, rho_ds_zm, invrs_rho_ds_zt,
    )

    # Sum contributions to turbulent production from standard term & damping
    lhs_tp_wp3 = lhs_adv_tp_wp3 + lhs_pr_tp_wp3

    # Calculate pressure terms 1 for w'^3
    lhs_pr1_wp3 = wp3_term_pr1_lhs(
        nzt, ngrdcol, gr,
        clubb_params[:, iC8], clubb_params[:, iC8b],
        invrs_tau_wp3_zt, Skw_zt,
        l_damp_wp3_Skw_squared,
    )

    # Calculate dissipation terms 1 for w'^2
    lhs_dp1_wp2 = wp2_term_dp1_lhs(
        nzm, ngrdcol, gr,
        C1_Skw_fnc, invrs_tau_C1_zm,
    )

    # Calculate buoyancy production of w'^2 and w'^2 pressure term 2
    rhs_bp_pr2_wp2 = wp2_terms_bp_pr2_rhs(
        nzm, ngrdcol, gr,
        clubb_params[:, iC_uu_buoy],
        thv_ds_zm, wpthvp,
    )

    # Calculate pressure terms 3 for w'^2
    rhs_pr3_wp2 = wp2_term_pr3_rhs(
        nzm, nzt, ngrdcol, gr,
        clubb_params[:, iC_uu_shr],
        clubb_params[:, iC_uu_buoy],
        thv_ds_zm, wpthvp, upwp,
        um, vpwp, vm,
    )

    # Calculate dissipation terms 1 for w'^2
    rhs_dp1_wp2 = wp2_term_dp1_rhs(
        nzm, ngrdcol, gr, C1_Skw_fnc,
        invrs_tau_C1_zm, w_tol_sqd, up2, vp2,
        l_damp_wp2_using_em,
    )

    # Calculate buoyancy production of w'^3 and w'^3 pressure term 2
    rhs_bp1_pr2_wp3 = wp3_terms_bp1_pr2_rhs(
        nzt, ngrdcol, gr, C11_Skw_fnc,
        thv_ds_zt, wp2thvp,
    )

    # Calculate pressure terms 1 for w'^3
    rhs_pr1_wp3 = wp3_term_pr1_rhs(
        nzt, ngrdcol, gr,
        clubb_params[:, iC8], clubb_params[:, iC8b],
        invrs_tau_wp3_zt, Skw_zt, wp3,
        l_damp_wp3_Skw_squared,
    )

    lhs_ta_wp3 = jnp.zeros((ndiags2, ngrdcol, nzt), dtype=jnp.float64)
    rhs_ta_wp3 = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    if l_explicit_turbulent_adv_wp3:
        # The turbulent advection term is being solved explicitly.
        # The w'^3 turbulent advection term is being solved explicitly.
        rhs_ta_wp3 = wp3_term_ta_explicit_rhs(
            nzm, nzt, ngrdcol, gr,
            wp4, rho_ds_zm, invrs_rho_ds_zt,
        )
    else:
        # The turbulent advection term is being solved implicitly.
        if iiPDF_type == iiPDF_ADG1:
            # The ADG1 PDF is used.
            wp3_term_ta_lhs_result = wp3_term_ta_ADG1_lhs(
                nzm, nzt, ngrdcol, gr,
                wp2, a1_coef, a1_coef_zt,
                a3_coef, a3_coef_zt,
                wp3_on_wp2, rho_ds_zm,
                rho_ds_zt, invrs_rho_ds_zt,
                l_standard_term_ta,
                l_partial_upwind_wp3,
            )
        elif iiPDF_type == iiPDF_new or iiPDF_type == iiPDF_new_hybrid:
            # The new PDF or the new hybrid PDF is used.
            lhs_ta_wp3 = wp3_term_ta_new_pdf_lhs(
                nzm, nzt, ngrdcol, gr, coef_wp4_implicit,
                wp2, rho_ds_zm, invrs_rho_ds_zt,
            )
            wp3_term_ta_lhs_result = wp3_term_ta_lhs_result.at[1, :, :].set(
                lhs_ta_wp3[0, :, :]
            )
            wp3_term_ta_lhs_result = wp3_term_ta_lhs_result.at[3, :, :].set(
                lhs_ta_wp3[1, :, :]
            )

    rhs, stats = wp23_rhs(
        nzm, nzt, ngrdcol, gr, dt, fcor_y,
        wp3_term_ta_lhs_result,
        lhs_diff_zm, lhs_diff_zt, lhs_diff_zm_crank, lhs_diff_zt_crank,
        lhs_tp_wp3, lhs_adv_tp_wp3, lhs_pr_tp_wp3,
        lhs_ta_wp3, lhs_dp1_wp2, rhs_dp1_wp2, lhs_pr1_wp2,
        rhs_pr1_wp2, lhs_pr1_wp3, rhs_pr1_wp3, rhs_bp_pr2_wp2,
        rhs_pr_dfsn_wp2, rhs_bp1_pr2_wp3, rhs_pr3_wp2, rhs_pr3_wp3,
        rhs_ta_wp3, rhs_pr_turb_wp3, rhs_pr_dfsn_wp3,
        wp2, wp3, wpup2, wpvp2,
        wpthvp, wp2thvp, wp2up, up2, vp2, upwp,
        C11_Skw_fnc, thv_ds_zm, thv_ds_zt,
        lhs_splat_wp2, lhs_splat_wp3,
        clubb_params,
        iiPDF_type,
        l_use_tke_in_wp2_wp3_K_dfsn,
        l_ho_nontrad_coriolis,
        stats,
    )

    lhs_ma_zm = term_ma_zm_lhs(
        nzm, nzt, ngrdcol, wm_zm,
        gr.invrs_dzm, gr.weights_zm2zt,
    )

    lhs_ma_zt = term_ma_zt_lhs(
        nzm, nzt, ngrdcol, wm_zt, gr.weights_zt2zm,
        gr.invrs_dzt, gr.invrs_dzm,
        l_upwind_xm_ma, gr.grid_dir,
    )

    lhs_diff_zt = lhs_diff_zt * clubb_params[:, iC12][None, :, None]

    if l_crank_nich_diff:
        lhs_diff_zm = lhs_diff_zm.at[:, :, 1:-1].set(
            one_half * lhs_diff_zm[:, :, 1:-1]
        )
        lhs_diff_zt = lhs_diff_zt.at[:, :, 1:-1].set(
            one_half * lhs_diff_zt[:, :, 1:-1]
        )

    lhs_ta_wp2 = wp2_term_ta_lhs(
        nzm, nzt, ngrdcol, gr,
        rho_ds_zt, invrs_rho_ds_zm,
    )

    lhs_ac_pr2_wp2 = wp2_terms_ac_pr2_lhs(
        nzm, nzt, ngrdcol, gr,
        clubb_params[:, iC_uu_shr], wm_zt,
    )

    lhs_ac_pr2_wp3 = wp3_terms_ac_pr2_lhs(
        nzm, nzt, ngrdcol, gr, C11_Skw_fnc, wm_zm,
    )

    lhs = wp23_lhs(
        nzm, nzt, ngrdcol, dt,
        wp3_term_ta_lhs_result,
        lhs_diff_zm, lhs_diff_zt, lhs_ma_zm,
        lhs_ma_zt, lhs_ta_wp2,
        lhs_tp_wp3,
        lhs_ac_pr2_wp2, lhs_ac_pr2_wp3, lhs_dp1_wp2,
        lhs_pr1_wp3, lhs_pr1_wp2, lhs_splat_wp2, lhs_splat_wp3,
    )

    wp2_old = wp2
    wp3_old = wp3

    up2, vp2, wp2, wp3, wp3_zm, wp2_zt, err_info, stats = wp23_solve(
        nzm, nzt, ngrdcol, gr, dt, lhs, rhs,
        lhs_ma_zm, lhs_dp1_wp2, lhs_diff_zm,
        lhs_ta_wp2, lhs_pr1_wp2, lhs_pr1_wp3,
        lhs_diff_zt, lhs_adv_tp_wp3, lhs_pr_tp_wp3,
        wp3_pr3_lhs, lhs_ma_zt,
        wp3_term_ta_lhs_result,
        wm_zm, wm_zt,
        sfc_elevation, C11_Skw_fnc,
        rho_ds_zm,
        upwp, vpwp, wprtp, wpthlp, rtp2, thlp2,
        clubb_params,
        penta_solve_method,
        fill_holes_type,
        l_min_wp2_from_corr_wx,
        l_use_tke_in_wp2_wp3_K_dfsn,
        l_use_wp3_lim_with_smth_Heaviside,
        l_wp2_fill_holes_tke,
        l_implemented,
        stats,
        up2, vp2, wp2, wp3, wp3_zm, wp2_zt, err_info,
    )

    if l_lmm_stepping:
        wp2 = one_half * (wp2_old + wp2)
        wp3 = one_half * (wp3_old + wp3)

    if stats.l_sample:
        stats = stats.update("wp2_zt", wp2_zt)
        stats = stats.update("wp3_zm", wp3_zm)

    return up2, vp2, wp2, wp3, err_info, stats


def wp23_solve(
    nzm, nzt, ngrdcol, gr, dt, lhs, rhs,
    lhs_ma_zm, lhs_dp1_wp2, lhs_diff_zm,
    lhs_ta_wp2, lhs_pr1_wp2, lhs_pr1_wp3,
    lhs_diff_zt, lhs_adv_tp_wp3, lhs_pr_tp_wp3,
    wp3_pr3_lhs, lhs_ma_zt,
    wp3_term_ta_lhs_result,
    wm_zm, wm_zt,
    sfc_elevation, C11_Skw_fnc,
    rho_ds_zm,
    upwp, vpwp, wprtp, wpthlp, rtp2, thlp2,
    clubb_params,
    penta_solve_method,
    fill_holes_type,
    l_min_wp2_from_corr_wx,
    l_use_tke_in_wp2_wp3_K_dfsn,
    l_use_wp3_lim_with_smth_Heaviside,
    l_wp2_fill_holes_tke,
    l_implemented,
    stats,
    up2, vp2, wp2, wp3, wp3_zm, wp2_zt, err_info,
):
    """Decompose, and back substitute the matrix for wp2/wp3.

    References:
    _Equations for CLUBB_ section 6.3
    """
    # Save the value of rhs, which will be overwritten with the solution as
    # part of the solving routine.
    rhs_save = rhs

    old_solut = jnp.zeros((ngrdcol, 2 * nzm - 1), dtype=jnp.float64)
    if penta_solve_method == penta_bicgstab:
        old_solut = old_solut.at[:, 0::2].set(wp2)
        old_solut = old_solut.at[:, 1::2].set(wp3)

    lhs_solve = lhs
    rhs_solve = rhs
    old_solut_solve = old_solut
    if l_force_descending_solves and gr.grid_dir_indx > 0:
        # Matrix solves are bit-different between ascending and descending.
        # This ensures matrices are solved in the same (descending) order,
        # which is useful for ensuring BFBness between grid modes
        # We need to flip in both the vertical dimensions and the bands for the lhs
        lhs_solve = lhs_solve[::-1, :, ::-1]
        rhs_solve = rhs_solve[:, ::-1]
        old_solut_solve = old_solut_solve[:, ::-1]

    l_need_rcond = bool(
        stats.l_sample and stats.var_on_stats_list("wp23_matrix_condt_num")
    )

    # Solve the system with LAPACK
    # Solve the system and return condition number
    # Note: When using lapack this can change the answer slightly
    err_info, solut, rcond = band_solve(
        "wp2_wp3", penta_solve_method,
        ngrdcol, 2, 2, 2 * nzm - 1, 1,
        l_implemented,
        lhs_solve, rhs_solve, err_info,
        old_soln=old_solut_solve,
        use_rcond=l_need_rcond,
    )

    if l_need_rcond and stats.l_sample:
        # Est. of the condition number of the w'^2/w^3 LHS matrix
        stats = stats.update("wp23_matrix_condt_num", one / rcond)

    if l_force_descending_solves and gr.grid_dir_indx > 0:
        # Flip the back to the ascending direction if we forced the solve
        # to be in descending mode
        solut = solut[:, ::-1]

    # Copy result into output arrays and clip
    wp2 = solut[:, 0::2]
    wp3 = solut[:, 1::2]

    if stats.l_sample:
        C_uu_shr_zeros = jnp.zeros((ngrdcol,), dtype=jnp.float64)
        C_uu_shr_plus_one = clubb_params[:, iC_uu_shr] + one
        C11_Skw_fnc_zeros = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        C11_Skw_fnc_plus_one = C11_Skw_fnc + one

        lhs_wp2_ac_term = wp2_terms_ac_pr2_lhs(
            nzm, nzt, ngrdcol, gr, C_uu_shr_zeros, wm_zt,
        )
        lhs_wp2_pr2_term = wp2_terms_ac_pr2_lhs(
            nzm, nzt, ngrdcol, gr, C_uu_shr_plus_one, wm_zt,
        )
        lhs_wp3_ac_term = wp3_terms_ac_pr2_lhs(
            nzm, nzt, ngrdcol, gr, C11_Skw_fnc_zeros, wm_zm,
        )
        lhs_wp3_pr2_term = wp3_terms_ac_pr2_lhs(
            nzm, nzt, ngrdcol, gr, C11_Skw_fnc_plus_one, wm_zm,
        )

        # Finalize implicit contributions for wp2
        # Finalize the implicit pieces of wp2/wp3 budgets after the solve.
        # The companion explicit pieces were opened earlier with begin_budget.
        stats_tmp_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(
            -gamma_over_implicit_ts * lhs_dp1_wp2[:, 1:-1] * wp2[:, 1:-1]
        )
        # w'^2 term dp1 has both implicit and explicit components.
        stats = stats.finalize_budget("wp2_dp1", stats_tmp_zm)

        stats_tmp_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        if gr.grid_dir_indx > 0:
            stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(
                -lhs_diff_zm[2, :, 1:-1] * wp2[:, :-2]
                -lhs_diff_zm[1, :, 1:-1] * wp2[:, 1:-1]
                -lhs_diff_zm[0, :, 1:-1] * wp2[:, 2:]
            )
        else:
            stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(
                -lhs_diff_zm[0, :, 1:-1] * wp2[:, 2:]
                -lhs_diff_zm[1, :, 1:-1] * wp2[:, 1:-1]
                -lhs_diff_zm[2, :, 1:-1] * wp2[:, :-2]
            )
        if l_crank_nich_diff or l_use_tke_in_wp2_wp3_K_dfsn:
            # dp2 has explicit and implicit pieces for CN/tke diffusion.
            # w'^2 term dp2 has both implicit and explicit components.
            stats = stats.finalize_budget("wp2_dp2", stats_tmp_zm)
        else:
            # Otherwise dp2 is fully implicit.
            # w'^2 term dp2 is completely implicit; call stat_update_var_pt.
            stats = stats.update("wp2_dp2", stats_tmp_zm)

        stats_tmp_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(
            -lhs_ta_wp2[1, :, 1:-1] * wp3[:, :-1]
            -lhs_ta_wp2[0, :, 1:-1] * wp3[:, 1:]
        )
        # w'^2 term ta is completely implicit; call stat_update_var_pt.
        stats = stats.update("wp2_ta", stats_tmp_zm)

        stats_tmp_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        if gr.grid_dir_indx > 0:
            stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(
                -lhs_ma_zm[2, :, 1:-1] * wp2[:, :-2]
                -lhs_ma_zm[1, :, 1:-1] * wp2[:, 1:-1]
                -lhs_ma_zm[0, :, 1:-1] * wp2[:, 2:]
            )
        else:
            stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(
                -lhs_ma_zm[0, :, 1:-1] * wp2[:, 2:]
                -lhs_ma_zm[1, :, 1:-1] * wp2[:, 1:-1]
                -lhs_ma_zm[2, :, 1:-1] * wp2[:, :-2]
            )
        # w'^2 term ma is completely implicit; call stat_update_var_pt.
        stats = stats.update("wp2_ma", stats_tmp_zm)

        stats_tmp_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(
            -lhs_wp2_ac_term[:, 1:-1] * wp2[:, 1:-1]
        )
        # w'^2 term ac is completely implicit; call stat_update_var_pt.
        stats = stats.update("wp2_ac", stats_tmp_zm)

        stats_tmp_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(
            -gamma_over_implicit_ts * lhs_pr1_wp2[:, 1:-1] * wp2[:, 1:-1]
        )
        # w'^2 term pr1 has both implicit and explicit components.
        stats = stats.finalize_budget("wp2_pr1", stats_tmp_zm)

        stats_tmp_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(
            -lhs_wp2_pr2_term[:, 1:-1] * wp2[:, 1:-1]
        )
        # w'^2 term pr2 has both implicit and explicit components.
        stats = stats.finalize_budget("wp2_pr2", stats_tmp_zm)

        # Finalize implicit contributions for wp3
        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(
            -gamma_over_implicit_ts * lhs_pr1_wp3[:, 1:-1] * wp3[:, 1:-1]
        )
        # w'^3 term pr1 has both implicit and explicit components.
        stats = stats.finalize_budget("wp3_pr1", stats_tmp_zt)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        if gr.grid_dir_indx > 0:
            stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(
                -lhs_diff_zt[2, :, 1:-1] * wp3[:, :-2]
                -lhs_diff_zt[1, :, 1:-1] * wp3[:, 1:-1]
                -lhs_diff_zt[0, :, 1:-1] * wp3[:, 2:]
            )
        else:
            stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(
                -lhs_diff_zt[0, :, 1:-1] * wp3[:, 2:]
                -lhs_diff_zt[1, :, 1:-1] * wp3[:, 1:-1]
                -lhs_diff_zt[2, :, 1:-1] * wp3[:, :-2]
            )
        if l_crank_nich_diff or l_use_tke_in_wp2_wp3_K_dfsn:
            # dp1 has explicit and implicit pieces for CN/tke diffusion.
            # w'^3 term dp1 has both implicit and explicit components.
            stats = stats.finalize_budget("wp3_dp1", stats_tmp_zt)
        else:
            # Otherwise dp1 is fully implicit.
            # w'^3 term dp1 is completely implicit; call stat_update_var_pt.
            stats = stats.update("wp3_dp1", stats_tmp_zt)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(
            -gamma_over_implicit_ts * wp3_term_ta_lhs_result[4, :, 1:-1] * wp3[:, :-2]
            -gamma_over_implicit_ts * wp3_term_ta_lhs_result[3, :, 1:-1] * wp2[:, 1:-2]
            -gamma_over_implicit_ts * wp3_term_ta_lhs_result[2, :, 1:-1] * wp3[:, 1:-1]
            -gamma_over_implicit_ts * wp3_term_ta_lhs_result[1, :, 1:-1] * wp2[:, 2:-1]
            -gamma_over_implicit_ts * wp3_term_ta_lhs_result[0, :, 1:-1] * wp3[:, 2:]
        )
        # w'^3 term ta has both implicit and explicit components.
        stats = stats.finalize_budget("wp3_ta", stats_tmp_zt)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(
            -gamma_over_implicit_ts * lhs_adv_tp_wp3[1, :, 1:-1] * wp2[:, 1:-2]
            -gamma_over_implicit_ts * lhs_adv_tp_wp3[0, :, 1:-1] * wp2[:, 2:-1]
        )
        # w'^3 term tp has both implicit and explicit components.
        stats = stats.finalize_budget("wp3_tp", stats_tmp_zt)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(
            -gamma_over_implicit_ts * lhs_pr_tp_wp3[1, :, 1:-1] * wp2[:, 1:-2]
            -gamma_over_implicit_ts * lhs_pr_tp_wp3[0, :, 1:-1] * wp2[:, 2:-1]
        )
        # w'^3 term pr_tp same as above tp term but opposite sign.
        stats = stats.finalize_budget("wp3_pr_tp", stats_tmp_zt)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(
            -wp3_pr3_lhs[4, :, 1:-1] * wp3[:, :-2]
            -wp3_pr3_lhs[3, :, 1:-1] * wp2[:, 1:-2]
            -wp3_pr3_lhs[2, :, 1:-1] * wp3[:, 1:-1]
            -wp3_pr3_lhs[1, :, 1:-1] * wp2[:, 2:-1]
            -wp3_pr3_lhs[0, :, 1:-1] * wp3[:, 2:]
        )
        # w'^3 pressure term 3 (pr3) has both implicit and explicit components.
        stats = stats.finalize_budget("wp3_pr3", stats_tmp_zt)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        if gr.grid_dir_indx > 0:
            stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(
                -lhs_ma_zt[2, :, 1:-1] * wp3[:, :-2]
                -lhs_ma_zt[1, :, 1:-1] * wp3[:, 1:-1]
                -lhs_ma_zt[0, :, 1:-1] * wp3[:, 2:]
            )
        else:
            stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(
                -lhs_ma_zt[0, :, 1:-1] * wp3[:, 2:]
                -lhs_ma_zt[1, :, 1:-1] * wp3[:, 1:-1]
                -lhs_ma_zt[2, :, 1:-1] * wp3[:, :-2]
            )
        # w'^3 term ma is completely implicit; call stat_update_var_pt.
        stats = stats.update("wp3_ma", stats_tmp_zt)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(
            -lhs_wp3_ac_term[:, 1:-1] * wp3[:, 1:-1]
        )
        # w'^3 term ac is completely implicit; call stat_update_var_pt.
        stats = stats.update("wp3_ac", stats_tmp_zt)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(
            -lhs_wp3_pr2_term[:, 1:-1] * wp3[:, 1:-1]
        )
        # w'^3 term pr2 has both implicit and explicit components.
        stats = stats.finalize_budget("wp3_pr2", stats_tmp_zt)

    if stats.l_sample:
        # Store previous value for effect of the positive definite scheme
        # Bracket the positive-definite/hole-filling adjustments so the pd
        # budget terms report only the limiter-induced change.
        stats = stats.begin_budget("wp2_pd", wp2 / dt)
        # up2_pd and vp2_pd are also modified in a call to pos_definite_variances
        # in advance_xp2_xpyp, so we need to use stat_modify here instead of stat_begin_update
        stats = stats.begin_budget("up2_pd", up2 / dt)
        stats = stats.begin_budget("vp2_pd", vp2 / dt)

    if fill_holes_type != 0:
        # Use a simple hole filling algorithm
        # upper_hf_level = nzm-1 since we are filling the zm levels
        wp2 = fill_holes_vertical(
            nzm, ngrdcol, w_tol_sqd,
            gr.k_lb_zm + gr.grid_dir_indx,
            gr.k_ub_zm - gr.grid_dir_indx,
            gr.dzm, rho_ds_zm, gr.grid_dir_indx,
            fill_holes_type,
            wp2,
        )

        if l_wp2_fill_holes_tke:
            wp2, up2, vp2 = fill_holes_wp2_from_horz_tke(
                nzm, ngrdcol, w_tol_sqd, 1, nzm - 2,
                wp2, up2, vp2,
            )

    if stats.l_sample:
        # Store updated value for effect of the positive definite scheme
        stats = stats.finalize_budget("wp2_pd", wp2 / dt)
        stats = stats.finalize_budget("up2_pd", up2 / dt, l_count_sample=False)
        stats = stats.finalize_budget("vp2_pd", vp2 / dt, l_count_sample=False)

    # Clip <w'^2> at a minimum and maximum threshold.
    # We attempt to clip extreme values of wp2 to prevent a crash of the
    # type found on the Climate Process Team ticket #49.  Chris Golaz found that
    # instability caused by large wp2 in CLUBB led unrealistic results in AM3.
    # -dschanen 11 Apr 2011

    # The value of <w'^2> is not allowed to become smaller than the threshold
    # value of w_tol^2.  Additionally, that threshold value may be boosted at
    # any grid level in order to keep the overall correlation of w and rt or
    # the overall correlation of w and theta-l between the values of
    # -max_mag_correlation_flux and max_mag_correlation_flux by boosting <w'^2>
    # rather than by limiting the magnitude of <w'rt'> or <w'thl'>.
    if l_min_wp2_from_corr_wx:
        # The overall correlation of w and rt is:
        #
        # corr_w_rt = wprtp / ( sqrt( wp2 ) * sqrt( rtp2 ) );
        #
        # and the overall correlation of w and thl is:
        #
        # corr_w_thl = wpthlp / ( sqrt( wp2 ) * sqrt( thlp2 ) ).
        #
        # Squaring both sides, the equations becomes:
        #
        # corr_w_rt^2 = wprtp^2 / ( wp2 * rtp2 ); and
        #
        # corr_w_thl^2 = wpthlp^2 / ( wp2 * thlp2 ).
        #
        # Using max_mag_correlation_flux for the correlation and then solving for
        # the minimum of wp2, the equation becomes:
        #
        # wp2|_min = max( wprtp^2 / ( rtp2 * max_mag_correlation_flux^2 ),
        #                 wpthlp^2 / ( thlp2 * max_mag_correlation_flux^2 ) ).
        #
        # And similarly for upwp and vpwp
        corr_max_sqd = max_mag_correlation_flux ** 2
        wp2_min_array = jnp.maximum(
            w_tol_sqd,
            wprtp ** 2 / (rtp2 * corr_max_sqd),
        )
        wp2_min_array = jnp.maximum(
            wp2_min_array,
            wpthlp ** 2 / (thlp2 * corr_max_sqd),
        )
        wp2_min_array = jnp.maximum(
            wp2_min_array,
            upwp ** 2 / (up2 * corr_max_sqd),
        )
        wp2_min_array = jnp.maximum(
            wp2_min_array,
            vpwp ** 2 / (vp2 * corr_max_sqd),
        )
        wp2_min_array = jnp.minimum(one, wp2_min_array)
    else:
        # Consider only the minimum tolerance threshold value for wp2.
        wp2_min_array = jnp.full((ngrdcol, nzm), w_tol_sqd, dtype=jnp.float64)

    wp2, stats = clip_variance(
        nzm, ngrdcol, gr, clip_wp2, dt, wp2_min_array,
        stats,
        wp2,
        wp2_max,
    )

    # Interpolate w'^2 from momentum levels to thermodynamic levels.
    # This is used for the clipping of w'^3 according to the value
    # of Sk_w now that w'^2 and w'^3 have been advanced one timestep.
    wp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, wp2, w_tol_sqd)

    # Clip w'^3 by limiting skewness.
    wp3, stats = clip_skewness(
        nzt, ngrdcol, gr, dt, sfc_elevation,
        clubb_params[:, iSkw_max_mag], wp2_zt,
        l_use_wp3_lim_with_smth_Heaviside,
        stats,
        wp3,
    )

    # Compute wp3_zm for output purposes
    wp3_zm = zt2zm(nzm, nzt, ngrdcol, gr, wp3)

    del rhs_save
    return up2, vp2, wp2, wp3, wp3_zm, wp2_zt, err_info, stats


def wp23_lhs(
    nzm, nzt, ngrdcol, dt,
    wp3_term_ta_lhs_result,
    lhs_diff_zm, lhs_diff_zt, lhs_ma_zm,
    lhs_ma_zt, lhs_ta_wp2,
    lhs_tp_wp3,
    lhs_ac_pr2_wp2, lhs_ac_pr2_wp3, lhs_dp1_wp2,
    lhs_pr1_wp3, lhs_pr1_wp2, lhs_splat_wp2, lhs_splat_wp3,
):
    """Compute LHS band diagonal matrix for w'^2 and w'^3.

    This subroutine computes the implicit portion
    of the w'^2 and w'^3 equations.

    Boundary conditions

      Both wp2 and wp3 used fixed-point boundary conditions.
      Therefore, anything set in the above loop at both the upper
      and lower boundaries would be overwritten here.  However, the
      above loop does not extend to the boundary levels.  An array
      with a value of 1 at the main diagonal on the left-hand side
      and with values of 0 at all other diagonals on the left-hand
      side will preserve the right-hand side value at that level.

         wp2(1)  wp3(1)  ...  wp3(nzt)  wp2(nzm)
        [  0.0     0.0          0.0       0.0  ]
        [  0.0     0.0          0.0       0.0  ]
        [  1.0     1.0   ...    1.0       1.0  ]
        [  0.0     0.0          0.0       0.0  ]
        [  0.0     0.0          0.0       0.0  ]

     WARNING: This subroutine has been optimized. Significant changes could
              noticeably  impact computational efficiency. See clubb:ticket:834
    """
    del nzt
    # Calculate invrs_dt
    invrs_dt = one / dt
    lhs = jnp.zeros((ndiags5, ngrdcol, 2 * nzm - 1), dtype=jnp.float64)

    # Lower (upper) boundary for w'2 on ascending (descending) grid.
    lhs = lhs.at[2, :, 0].set(one)
    # Lower (upper) boundary for w'3 on ascending (descending) grid.
    lhs = lhs.at[2, :, 1].set(one)

    # Combine terms to calculate non-boundary lhs values
    # ------ w'2 ------
    # LHS mean advection (ma) and diffusion (diff) terms
    lhs = lhs.at[0, :, 2:-1:2].set(lhs_ma_zm[0, :, 1:-1] + lhs_diff_zm[0, :, 1:-1])
    # LHS turbulent advection (ta) term.
    lhs = lhs.at[1, :, 2:-1:2].set(lhs_ta_wp2[0, :, 1:-1])
    # LHS mean advection (ma) and diffusion (diff) terms
    # LHS accumulation (ac) term and pressure term 2 (pr2).
    # LHS dissipation term 1 (dp1).
    # Note:  An "over-implicit" weighted time step is applied to this term.
    #        A weighting factor of greater than 1 may be used to make the term
    #        more numerically stable (see note below for w'^3 LHS turbulent
    #        advection (ta) term).
    # LHS time tendency.
    lhs = lhs.at[2, :, 2:-1:2].set(
        lhs_ma_zm[1, :, 1:-1] + lhs_diff_zm[1, :, 1:-1]
        + lhs_ac_pr2_wp2[:, 1:-1]
        + gamma_over_implicit_ts * lhs_dp1_wp2[:, 1:-1]
        + invrs_dt
    )
    # LHS turbulent advection (ta) term.
    lhs = lhs.at[3, :, 2:-1:2].set(lhs_ta_wp2[1, :, 1:-1])
    # LHS mean advection (ma) and diffusion (diff) terms
    lhs = lhs.at[4, :, 2:-1:2].set(lhs_ma_zm[2, :, 1:-1] + lhs_diff_zm[2, :, 1:-1])

    # ------ w'3 ------
    # LHS mean advection (ma) and diffusion (diff) terms
    lhs = lhs.at[0, :, 3:-2:2].set(lhs_ma_zt[0, :, 1:-1] + lhs_diff_zt[0, :, 1:-1])
    # LHS turbulent production (tp) term.
    # Note:  An "over-implicit" weighted time step is applied to this term.
    lhs = lhs.at[1, :, 3:-2:2].set(gamma_over_implicit_ts * lhs_tp_wp3[0, :, 1:-1])
    # LHS mean advection (ma) and diffusion (diff) terms
    # LHS accumulation (ac) term and pressure term 2 (pr2).
    # LHS pressure term 1 (pr1).
    # Note:  An "over-implicit" weighted time step is applied to this term.
    # Add implicit splatting
    # LHS time tendency.
    lhs = lhs.at[2, :, 3:-2:2].set(
        lhs_ma_zt[1, :, 1:-1] + lhs_diff_zt[1, :, 1:-1]
        + lhs_ac_pr2_wp3[:, 1:-1]
        + gamma_over_implicit_ts * lhs_pr1_wp3[:, 1:-1]
        + lhs_splat_wp3[:, 1:-1]
        + invrs_dt
    )
    # LHS turbulent production (tp) term.
    # Note:  An "over-implicit" weighted time step is applied to this term.
    lhs = lhs.at[3, :, 3:-2:2].set(gamma_over_implicit_ts * lhs_tp_wp3[1, :, 1:-1])
    # LHS mean advection (ma) and diffusion (diff) terms
    lhs = lhs.at[4, :, 3:-2:2].set(lhs_ma_zt[2, :, 1:-1] + lhs_diff_zt[2, :, 1:-1])

    # Upper (lower) boundary for w'3 on ascending (descending) grid.
    lhs = lhs.at[2, :, 2 * nzm - 3].set(one)
    # Upper (lower) boundary for w'2 on ascending (desceding) grid.
    lhs = lhs.at[2, :, 2 * nzm - 2].set(one)

    # LHS pressure term 1 (pr1) for wp2
    # Note:  An "over-implicit" weighted time step is applied to this term.
    #        A weighting factor of greater than 1 may be used to make the term
    #        more numerically stable (see note below for w'^3 LHS turbulent
    #        advection (ta) term).
    # Reference:
    # https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:wp2_pr
    # Add terms to lhs
    lhs = lhs.at[2, :, 2:-1:2].add(
        gamma_over_implicit_ts * lhs_pr1_wp2[:, 1:-1]
    )

    # Add implicit splatting to wp2
    lhs = lhs.at[2, :, 2:-1:2].add(lhs_splat_wp2[:, 1:-1])

    # LHS turbulent advection (ta) term for wp3
    if not l_explicit_turbulent_adv_wp3:
        # Note:  An "over-implicit" weighted time step is applied to this term.
        #        The weight of the implicit portion of this term is controlled
        #        by the factor gamma_over_implicit_ts (abbreviated "gamma" in
        #        the expression below).  A factor is added to the right-hand
        #        side of the equation in order to balance a weight that is not
        #        equal to 1, such that:
        #             -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
        #        where X is the variable that is being solved for in a
        #        predictive equation (w'^3 in this case), y(t) is the
        #        linearized portion of the term that gets treated implicitly,
        #        and RHS is the portion of the term that is always treated
        #        explicitly (in the case of the w'^3 turbulent advection term,
        #        RHS = 0).  A weight of greater than 1 can be applied to make
        #        the term more numerically stable.
        # Add terms to lhs
        lhs = lhs.at[:, :, 3:-2:2].add(
            gamma_over_implicit_ts * wp3_term_ta_lhs_result[:, :, 1:-1]
        )

    return lhs


def wp23_rhs(
    nzm, nzt, ngrdcol, gr, dt, fcor_y,
    wp3_term_ta_lhs_result,
    lhs_diff_zm, lhs_diff_zt, lhs_diff_zm_crank, lhs_diff_zt_crank,
    lhs_tp_wp3, lhs_adv_tp_wp3, lhs_pr_tp_wp3,
    lhs_ta_wp3, lhs_dp1_wp2, rhs_dp1_wp2, lhs_pr1_wp2,
    rhs_pr1_wp2, lhs_pr1_wp3, rhs_pr1_wp3, rhs_bp_pr2_wp2,
    rhs_pr_dfsn_wp2, rhs_bp1_pr2_wp3, rhs_pr3_wp2, rhs_pr3_wp3,
    rhs_ta_wp3, rhs_pr_turb_wp3, rhs_pr_dfsn_wp3,
    wp2, wp3, wpup2, wpvp2,
    wpthvp, wp2thvp, wp2up, up2, vp2, upwp,
    C11_Skw_fnc, thv_ds_zm, thv_ds_zt,
    lhs_splat_wp2, lhs_splat_wp3,
    clubb_params,
    iiPDF_type,
    l_use_tke_in_wp2_wp3_K_dfsn,
    l_ho_nontrad_coriolis,
    stats,
):
    """Compute RHS vector for w'^2 and w'^3.

    This subroutine computes the explicit portion of
    the w'^2 and w'^3 equations.

    Notes:
         For LHS turbulent advection (ta) term.
            An "over-implicit" weighted time step is applied to this term.
            The weight of the implicit portion of this term is controlled
            by the factor gamma_over_implicit_ts (abbreviated "gamma" in
            the expression below).  A factor is added to the right-hand
            side of the equation in order to balance a weight that is not
            equal to 1, such that:
                 -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
            where X is the variable that is being solved for in a
            predictive equation (w'^3 in this case), y(t) is the
            linearized portion of the term that gets treated implicitly,
            and RHS is the portion of the term that is always treated
            explicitly (in the case of the w'^3 turbulent advection term,
            RHS = 0).  A weight of greater than 1 can be applied to make
            the term more numerically stable.

     WARNING: This subroutine has been optimized. Significant changes could
              noticeably  impact computational efficiency. See clubb:ticket:834
    """
    # Calculate invers_dt
    invrs_dt = one / dt
    # Initialize to zero
    rhs = jnp.zeros((ngrdcol, 2 * nzm - 1), dtype=jnp.float64)

    # Experimental bouyancy term
    # Experimental term from CLUBB TRAC ticket #411
    rhs = rhs.at[:, 3:-2:2].set(rhs_pr_turb_wp3[:, 1:-1] + rhs_pr_dfsn_wp3[:, 1:-1])
    rhs = rhs.at[:, 2:-1:2].add(rhs_pr_dfsn_wp2[:, 1:-1])

    # These lines are for the diffusional term with a Crank-Nicholson
    # time step.  They are not used for completely implicit diffusion.
    if l_crank_nich_diff:
        # Add diffusion terms
        if gr.grid_dir_indx > 0:
            rhs = rhs.at[:, 2:-1:2].add(
                -lhs_diff_zm_crank[2, :, 1:-1] * wp2[:, :-2]
                -lhs_diff_zm_crank[1, :, 1:-1] * wp2[:, 1:-1]
                -lhs_diff_zm_crank[0, :, 1:-1] * wp2[:, 2:]
            )
            rhs = rhs.at[:, 3:-2:2].add(
                -lhs_diff_zt_crank[2, :, 1:-1] * wp3[:, :-2]
                -lhs_diff_zt_crank[1, :, 1:-1] * wp3[:, 1:-1]
                -lhs_diff_zt_crank[0, :, 1:-1] * wp3[:, 2:]
            )
        else:
            rhs = rhs.at[:, 2:-1:2].add(
                -lhs_diff_zm_crank[0, :, 1:-1] * wp2[:, 2:]
                -lhs_diff_zm_crank[1, :, 1:-1] * wp2[:, 1:-1]
                -lhs_diff_zm_crank[2, :, 1:-1] * wp2[:, :-2]
            )
            rhs = rhs.at[:, 3:-2:2].add(
                -lhs_diff_zt_crank[0, :, 1:-1] * wp3[:, 2:]
                -lhs_diff_zt_crank[1, :, 1:-1] * wp3[:, 1:-1]
                -lhs_diff_zt_crank[2, :, 1:-1] * wp3[:, :-2]
            )

    # This code block adds terms to the right-hand side so that TKE is being
    # used in eddy diffusion instead of just wp2 or wp3.  For example, in the
    # wp2 equation, if this flag is false, the eddy diffusion term would
    # normally be completely implicit (hence no right-hand side contribution),
    # and equal to +d/dz((K+nu)d/dz(wp2)).  With this flag set to true, the eddy
    # diffusion term will be +d/dz((K+nu)d/dz(up2+vp2+wp2)), but the up2 and vp2
    # parts are added on here as if they were right-hand side terms. For the wp3
    # equation, with this flag false, the eddy diffusion term is
    # +d/dz((K+nu)d/dz(wp3)), but with this flag true, it will be
    # +d/dz((K+nu)d/dz(wpup2+wpvp2+wp3)).
    if l_use_tke_in_wp2_wp3_K_dfsn:
        if gr.grid_dir_indx > 0:
            rhs = rhs.at[:, 2:-1:2].add(
                -lhs_diff_zm[2, :, 1:-1] * (up2[:, :-2] + vp2[:, :-2])
                -lhs_diff_zm[1, :, 1:-1] * (up2[:, 1:-1] + vp2[:, 1:-1])
                -lhs_diff_zm[0, :, 1:-1] * (up2[:, 2:] + vp2[:, 2:])
            )
            rhs = rhs.at[:, 3:-2:2].add(
                -lhs_diff_zt[2, :, 1:-1] * (wpup2[:, :-2] + wpvp2[:, :-2])
                -lhs_diff_zt[1, :, 1:-1] * (wpup2[:, 1:-1] + wpvp2[:, 1:-1])
                -lhs_diff_zt[2, :, 1:-1] * (wpup2[:, 2:] + wpvp2[:, 2:])
            )
        else:
            rhs = rhs.at[:, 2:-1:2].add(
                -lhs_diff_zm[0, :, 1:-1] * (up2[:, 2:] + vp2[:, 2:])
                -lhs_diff_zm[1, :, 1:-1] * (up2[:, 1:-1] + vp2[:, 1:-1])
                -lhs_diff_zm[2, :, 1:-1] * (up2[:, :-2] + vp2[:, :-2])
            )
            rhs = rhs.at[:, 3:-2:2].add(
                -lhs_diff_zt[0, :, 1:-1] * (wpup2[:, 2:] + wpvp2[:, 2:])
                -lhs_diff_zt[1, :, 1:-1] * (wpup2[:, 1:-1] + wpvp2[:, 1:-1])
                -lhs_diff_zt[0, :, 1:-1] * (wpup2[:, :-2] + wpvp2[:, :-2])
            )

    # Add pressure terms and splat terms
    rhs = rhs.at[:, 2:-1:2].add(rhs_pr1_wp2[:, 1:-1])
    rhs = rhs.at[:, 2:-1:2].add(
        (one - gamma_over_implicit_ts)
        * (-lhs_pr1_wp2[:, 1:-1] * wp2[:, 1:-1])
    )

    if l_ho_nontrad_coriolis:
        # Add the nontraditional Coriolis term for wp2
        # Hing Ong, 19 July 2025
        rhs = rhs.at[:, 2:-1:2].add(two * fcor_y[:, None] * upwp[:, 1:-1])
        # Add the nontraditional Coriolis term for wp3
        # Hing Ong, 1 Septempber 2025
        rhs = rhs.at[:, 3:-2:2].add(three * fcor_y[:, None] * wp2up[:, 1:-1])

    # Combine terms
    # ------ Combine terms for 3rd moment of vertical velocity, <w'^3> ------ !
    # RHS time tendency.
    rhs = rhs.at[:, 3:-2:2].add(invrs_dt * wp3[:, 1:-1])
    # RHS contribution from "over-implicit" turbulent production (tp) term.
    rhs = rhs.at[:, 3:-2:2].add(
        (one - gamma_over_implicit_ts)
        * (
            -lhs_tp_wp3[0, :, 1:-1] * wp2[:, 2:-1]
            -lhs_tp_wp3[1, :, 1:-1] * wp2[:, 1:-2]
        )
    )
    # RHS buoyancy production (bp) term and pressure term 2 (pr2).
    rhs = rhs.at[:, 3:-2:2].add(rhs_bp1_pr2_wp3[:, 1:-1])
    # RHS pressure term 1
    rhs = rhs.at[:, 3:-2:2].add(rhs_pr1_wp3[:, 1:-1])
    # RHS "over implicit" pressure term 1 (pr1).
    rhs = rhs.at[:, 3:-2:2].add(
        (one - gamma_over_implicit_ts)
        * (-lhs_pr1_wp3[:, 1:-1] * wp3[:, 1:-1])
    )

    # ------ Combine terms for 2nd moment of vertical velocity, <w'^2> ------ !
    # RHS time tendency.
    rhs = rhs.at[:, 2:-1:2].add(invrs_dt * wp2[:, 1:-1])
    # RHS buoyancy production (bp) term and pressure term 2 (pr2).
    rhs = rhs.at[:, 2:-1:2].add(rhs_bp_pr2_wp2[:, 1:-1])
    # RHS pressure term 3 (pr3).
    rhs = rhs.at[:, 2:-1:2].add(rhs_pr3_wp2[:, 1:-1])
    # RHS dissipation term 1 (dp1).
    rhs = rhs.at[:, 2:-1:2].add(rhs_dp1_wp2[:, 1:-1])
    # RHS "over implicit" pressure term 1 (pr1).
    rhs = rhs.at[:, 2:-1:2].add(
        (one - gamma_over_implicit_ts)
        * (-lhs_dp1_wp2[:, 1:-1] * wp2[:, 1:-1])
    )

    if l_explicit_turbulent_adv_wp3:
        # The turbulent advection term is being solved explicitly.
        # Add RHS turbulent advection (ta) terms
        rhs = rhs.at[:, 3:-2:2].add(rhs_ta_wp3[:, 1:-1])
    else:
        # The turbulent advection term is being solved implicitly. See note above
        if iiPDF_type == iiPDF_ADG1:
            # The ADG1 PDF is used.
            # Brian -- come up with better method for testing ascending vs descending
            rhs = rhs.at[:, 3:-2:2].add(
                (one - gamma_over_implicit_ts)
                * (
                    -wp3_term_ta_lhs_result[0, :, 1:-1] * wp3[:, 2:]
                    -wp3_term_ta_lhs_result[1, :, 1:-1] * wp2[:, 2:-1]
                    -wp3_term_ta_lhs_result[2, :, 1:-1] * wp3[:, 1:-1]
                    -wp3_term_ta_lhs_result[3, :, 1:-1] * wp2[:, 1:-2]
                    -wp3_term_ta_lhs_result[4, :, 1:-1] * wp3[:, :-2]
                )
            )
        elif iiPDF_type == iiPDF_new or iiPDF_type == iiPDF_new_hybrid:
            # The new PDF or the new hybrid PDF is used.
            # Add terms
            rhs = rhs.at[:, 3:-2:2].add(
                (one - gamma_over_implicit_ts)
                * (
                    -lhs_ta_wp3[0, :, 1:-1] * wp2[:, 2:-1]
                    -lhs_ta_wp3[1, :, 1:-1] * wp2[:, 1:-2]
                )
            )

    # --------- Boundary Conditions ---------

    # Both wp2 and wp3 used fixed-point boundary conditions.
    # Therefore, anything set in the above loop at both the upper
    # and lower boundaries would be overwritten here.  However, the
    # above loop does not extend to the boundary levels.  An array
    # with a value of 1 at the main diagonal on the left-hand side
    # and with values of 0 at all other diagonals on the left-hand
    # side will preserve the right-hand side value at that level.

    # The value of w'^2 at the lower boundary will remain the same.
    # When the lower boundary is at the surface, the surface value of
    # w'^2 is set in subroutine calc_surface_varnce (surface_varnce_module.F).

    # The value of w'^3 at the lower boundary will be 0.

    # The value of w'^2 at the upper boundary will be set to the threshold
    # minimum value of w_tol_sqd.

    # The value of w'^3 at the upper boundary will be set to 0.
    if gr.grid_dir_indx > 0:
        # Ascending Grid
        rhs = rhs.at[:, 0].set(wp2[:, gr.k_lb_zm])
        rhs = rhs.at[:, 1].set(zero)
        rhs = rhs.at[:, 2 * nzt - 1].set(zero)
        rhs = rhs.at[:, 2 * nzm - 2].set(w_tol_sqd)
    else:
        # Descending Grid
        rhs = rhs.at[:, 2 * nzm - 2].set(wp2[:, gr.k_lb_zm])
        rhs = rhs.at[:, 2 * nzt - 1].set(zero)
        rhs = rhs.at[:, 1].set(zero)
        rhs = rhs.at[:, 0].set(w_tol_sqd)

    # --------- Statistics output ---------
    if stats.l_sample:
        C_uu_buoy_zeros = jnp.zeros((ngrdcol,), dtype=jnp.float64)
        C_uu_buoy_plus_one = clubb_params[:, iC_uu_buoy] + one
        C11_Skw_fnc_zeros = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        C11_Skw_fnc_plus_one = C11_Skw_fnc + one

        rhs_bp_wp2 = wp2_terms_bp_pr2_rhs(
            nzm, ngrdcol, gr, C_uu_buoy_zeros, thv_ds_zm, wpthvp,
        )
        rhs_pr2_wp2 = wp2_terms_bp_pr2_rhs(
            nzm, ngrdcol, gr, C_uu_buoy_plus_one, thv_ds_zm, wpthvp,
        )
        rhs_bp1_wp3 = wp3_terms_bp1_pr2_rhs(
            nzt, ngrdcol, gr, C11_Skw_fnc_zeros, thv_ds_zt, wp2thvp,
        )
        rhs_pr2_wp3 = wp3_terms_bp1_pr2_rhs(
            nzt, ngrdcol, gr, C11_Skw_fnc_plus_one, thv_ds_zt, wp2thvp,
        )

        stats_tmp_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(rhs_bp_wp2[:, 1:-1])
        stats = stats.update("wp2_bp", stats_tmp_zm)

        if l_ho_nontrad_coriolis:
            stats_tmp_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
            stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(two * fcor_y[:, None] * upwp[:, 1:-1])
            stats = stats.update("wp2_nct", stats_tmp_zm)

        stats_tmp_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(rhs_pr_dfsn_wp2[:, 1:-1])
        stats = stats.update("wp2_pr_dfsn", stats_tmp_zm)

        stats_tmp_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(-lhs_splat_wp2[:, 1:-1] * wp2[:, 1:-1])
        stats = stats.update("wp2_splat", stats_tmp_zm)

        stats_tmp_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(rhs_pr3_wp2[:, 1:-1])
        stats = stats.update("wp2_pr3", stats_tmp_zm)

        stats_tmp_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(-rhs_pr2_wp2[:, 1:-1])
        stats = stats.begin_budget("wp2_pr2", stats_tmp_zm)

        stats_tmp_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(-rhs_dp1_wp2[:, 1:-1])
        stats = stats.begin_budget("wp2_dp1", stats_tmp_zm)

        stats_tmp_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(
            (one - gamma_over_implicit_ts)
            * (-lhs_dp1_wp2[:, 1:-1] * wp2[:, 1:-1])
        )
        stats = stats.update_budget("wp2_dp1", stats_tmp_zm)

        if l_crank_nich_diff or l_use_tke_in_wp2_wp3_K_dfsn:
            stats_tmp_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
            if l_crank_nich_diff:
                if gr.grid_dir_indx > 0:
                    stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].add(
                        lhs_diff_zm_crank[2, :, 1:-1] * wp2[:, :-2]
                        + lhs_diff_zm_crank[1, :, 1:-1] * wp2[:, 1:-1]
                        + lhs_diff_zm_crank[0, :, 1:-1] * wp2[:, 2:]
                    )
                else:
                    stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].add(
                        lhs_diff_zm_crank[0, :, 1:-1] * wp2[:, 2:]
                        + lhs_diff_zm_crank[1, :, 1:-1] * wp2[:, 1:-1]
                        + lhs_diff_zm_crank[2, :, 1:-1] * wp2[:, :-2]
                    )
            if l_use_tke_in_wp2_wp3_K_dfsn:
                if gr.grid_dir_indx > 0:
                    stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].add(
                        lhs_diff_zm[2, :, 1:-1] * (up2[:, :-2] + vp2[:, :-2])
                        + lhs_diff_zm[1, :, 1:-1] * (up2[:, 1:-1] + vp2[:, 1:-1])
                        + lhs_diff_zm[0, :, 1:-1] * (up2[:, 2:] + vp2[:, 2:])
                    )
                else:
                    stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].add(
                        lhs_diff_zm[0, :, 1:-1] * (up2[:, 2:] + vp2[:, 2:])
                        + lhs_diff_zm[1, :, 1:-1] * (up2[:, 1:-1] + vp2[:, 1:-1])
                        + lhs_diff_zm[2, :, 1:-1] * (up2[:, :-2] + vp2[:, :-2])
                    )
            stats = stats.begin_budget("wp2_dp2", stats_tmp_zm)

        stats_tmp_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(-rhs_pr1_wp2[:, 1:-1])
        stats = stats.begin_budget("wp2_pr1", stats_tmp_zm)
        stats_tmp_zm = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        stats_tmp_zm = stats_tmp_zm.at[:, 1:-1].set(
            (one - gamma_over_implicit_ts)
            * (-lhs_pr1_wp2[:, 1:-1] * wp2[:, 1:-1])
        )
        stats = stats.update_budget("wp2_pr1", stats_tmp_zm)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(rhs_bp1_wp3[:, 1:-1])
        stats = stats.update("wp3_bp1", stats_tmp_zt)

        if l_ho_nontrad_coriolis:
            stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
            stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(three * fcor_y[:, None] * wp2up[:, 1:-1])
            stats = stats.update("wp3_nct", stats_tmp_zt)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(rhs_pr_turb_wp3[:, 1:-1])
        stats = stats.update("wp3_pr_turb", stats_tmp_zt)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(rhs_pr_dfsn_wp3[:, 1:-1])
        stats = stats.update("wp3_pr_dfsn", stats_tmp_zt)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(-lhs_splat_wp3[:, 1:-1] * wp3[:, 1:-1])
        stats = stats.update("wp3_splat", stats_tmp_zt)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(-rhs_pr2_wp3[:, 1:-1])
        stats = stats.begin_budget("wp3_pr2", stats_tmp_zt)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(-rhs_pr1_wp3[:, 1:-1])
        stats = stats.begin_budget("wp3_pr1", stats_tmp_zt)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(
            (one - gamma_over_implicit_ts)
            * (-lhs_pr1_wp3[:, 1:-1] * wp3[:, 1:-1])
        )
        stats = stats.update_budget("wp3_pr1", stats_tmp_zt)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        if l_explicit_turbulent_adv_wp3:
            stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(-rhs_ta_wp3[:, 1:-1])
        elif iiPDF_type == iiPDF_ADG1:
            stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(
                -(one - gamma_over_implicit_ts)
                * (
                    -wp3_term_ta_lhs_result[0, :, 1:-1] * wp3[:, 2:]
                    -wp3_term_ta_lhs_result[1, :, 1:-1] * wp2[:, 2:-1]
                    -wp3_term_ta_lhs_result[2, :, 1:-1] * wp3[:, 1:-1]
                    -wp3_term_ta_lhs_result[3, :, 1:-1] * wp2[:, 1:-2]
                    -wp3_term_ta_lhs_result[4, :, 1:-1] * wp3[:, :-2]
                )
            )
        elif iiPDF_type == iiPDF_new or iiPDF_type == iiPDF_new_hybrid:
            stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(
                -(one - gamma_over_implicit_ts)
                * (
                    -lhs_ta_wp3[0, :, 1:-1] * wp2[:, 2:-1]
                    -lhs_ta_wp3[1, :, 1:-1] * wp2[:, 1:-2]
                )
            )
        stats = stats.begin_budget("wp3_ta", stats_tmp_zt)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(
            -(one - gamma_over_implicit_ts)
            * (
                -lhs_adv_tp_wp3[0, :, 1:-1] * wp2[:, 2:-1]
                -lhs_adv_tp_wp3[1, :, 1:-1] * wp2[:, 1:-2]
            )
        )
        stats = stats.begin_budget("wp3_tp", stats_tmp_zt)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(
            -(one - gamma_over_implicit_ts)
            * (
                -lhs_pr_tp_wp3[0, :, 1:-1] * wp2[:, 2:-1]
                -lhs_pr_tp_wp3[1, :, 1:-1] * wp2[:, 1:-2]
            )
        )
        stats = stats.begin_budget("wp3_pr_tp", stats_tmp_zt)

        stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].set(rhs_pr3_wp3[:, 1:-1])
        stats = stats.begin_budget("wp3_pr3", stats_tmp_zt)

        if l_crank_nich_diff or l_use_tke_in_wp2_wp3_K_dfsn:
            stats_tmp_zt = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
            if l_crank_nich_diff:
                if gr.grid_dir_indx > 0:
                    stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].add(
                        lhs_diff_zt[2, :, 1:-1] * wp3[:, :-2]
                        + lhs_diff_zt[1, :, 1:-1] * wp3[:, 1:-1]
                        + lhs_diff_zt[0, :, 1:-1] * wp3[:, 2:]
                    )
                else:
                    stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].add(
                        lhs_diff_zt[0, :, 1:-1] * wp3[:, 2:]
                        + lhs_diff_zt[1, :, 1:-1] * wp3[:, 1:-1]
                        + lhs_diff_zt[2, :, 1:-1] * wp3[:, :-2]
                    )
            if l_use_tke_in_wp2_wp3_K_dfsn:
                if gr.grid_dir_indx > 0:
                    stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].add(
                        lhs_diff_zt[2, :, 1:-1] * (wpup2[:, :-2] + wpvp2[:, :-2])
                        + lhs_diff_zt[1, :, 1:-1] * (wpup2[:, 1:-1] + wpvp2[:, 1:-1])
                        + lhs_diff_zt[0, :, 1:-1] * (wpup2[:, 2:] + wpvp2[:, 2:])
                    )
                else:
                    stats_tmp_zt = stats_tmp_zt.at[:, 1:-1].add(
                        lhs_diff_zt[0, :, 1:-1] * (wpup2[:, 2:] + wpvp2[:, 2:])
                        + lhs_diff_zt[1, :, 1:-1] * (wpup2[:, 1:-1] + wpvp2[:, 1:-1])
                        + lhs_diff_zt[2, :, 1:-1] * (wpup2[:, :-2] + wpvp2[:, :-2])
                    )
            stats = stats.begin_budget("wp3_dp1", stats_tmp_zt)

    return rhs, stats


def wp2_term_ta_lhs(nzm, nzt, ngrdcol, gr, rho_ds_zt, invrs_rho_ds_zm):
    """Turbulent advection term for w'^2:  implicit portion of the code.

    The d(w'^2)/dt equation contains a turbulent advection term:

    - (1/rho_ds) * d( rho_ds * w'^3 )/dz.

    The term is solved for completely implicitly, such that:

    - (1/rho_ds) * d( rho_ds * w'^3(t+1) )/dz.

    Note:  When the term is brought over to the left-hand side, the sign
           is reversed and the leading "-" in front of the term is changed
           to a "+".

    invrs_dzm(k) = 1 / ( zt(k) - zt(k-1) )
    """
    del nzt
    lhs_ta_wp2 = jnp.zeros((ndiags2, ngrdcol, nzm), dtype=jnp.float64)
    fac = invrs_rho_ds_zm[:, 1:-1] * gr.invrs_dzm[:, 1:-1]
    # Thermodynamic superdiagonal: [ x wp3(k,<t+1>) ]
    lhs_ta_wp2 = lhs_ta_wp2.at[0, :, 1:-1].set(fac * rho_ds_zt[:, 1:])
    # Thermodynamic subdiagonal: [ x wp3(k-1,<t+1>) ]
    lhs_ta_wp2 = lhs_ta_wp2.at[1, :, 1:-1].set(-fac * rho_ds_zt[:, :-1])
    return lhs_ta_wp2


def wp2_terms_ac_pr2_lhs(nzm, nzt, ngrdcol, gr, C_uu_shr, wm_zt):
    """Accumulation of w'^2 and w'^2 pressure term 2:  implicit portion of the
    code.

    The d(w'^2)/dt equation contains an accumulation term:

    - 2 w'^2 dw/dz;

    and pressure term 2:

    - C_5 ( -2 w'^2 dw/dz + 2 (g/th_0) w'th_v' ).

    Note 1:  When the term is brought over to the left-hand side, the sign
             is reversed and the leading "-" in front of the "2" is changed
             to a "+".
    Note 2:  We have broken C5 up into C_uu_shr for this term
             and C_uu_buoy for the buoyancy term.
    """
    del nzt
    lhs_ac_pr2_wp2 = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    # Momentum main diagonal: [ x wp2(k,<t+1>) ]
    lhs_ac_pr2_wp2 = lhs_ac_pr2_wp2.at[:, 1:-1].set(
        (one - C_uu_shr[:, None]) * two * gr.invrs_dzm[:, 1:-1]
        * (wm_zt[:, 1:] - wm_zt[:, :-1])
    )
    return lhs_ac_pr2_wp2


def wp2_term_dp1_lhs(nzm, ngrdcol, gr, C1_Skw_fnc, invrs_tau1m):
    """Dissipation term 1 for w'^2:  implicit portion of the code.

    The d(w'^2)/dt equation contains dissipation term 1:

    - ( C_1 / tau_1m ) w'^2.

    Since w'^2 has a minimum threshold, the term should be damped only to that
    threshold.  The term becomes:

    - ( C_1 / tau_1m ) * ( w'^2 - threshold ).

    Note:  When the implicit term is brought over to the left-hand side, the
           sign is reversed and the leading "-" in front of the term is
           changed to a "+".
    """
    del gr
    lhs_dp1_wp2 = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    # Momentum main diagonal: [ x wp2(k,<t+1>) ]
    lhs_dp1_wp2 = lhs_dp1_wp2.at[:, 1:-1].set(
        C1_Skw_fnc[:, 1:-1] * invrs_tau1m[:, 1:-1]
    )
    return lhs_dp1_wp2


def wp2_term_pr1_lhs(nzm, ngrdcol, gr, C4, invrs_tau_C4_zm):
    """Pressure term 1 for w'^2:  implicit portion of the code.

    The d(w'^2)/dt equation contains pressure term 1:

    - ( C_4 / tau_1m ) * ( w'^2 - (2/3)*em ),

    where em = (1/2) * ( w'^2 + u'^2 + v'^2 ).

    Pressure term 1 has both implicit and explicit components.  The implicit
    portion is:

    - ( C_4 / tau_1m ) * (2/3) * w'^2(t+1);
    """
    del gr
    lhs_pr1_wp2 = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    # Momentum main diagonal: [ x wp2(k,<t+1>) ]
    lhs_pr1_wp2 = lhs_pr1_wp2.at[:, 1:-1].set(
        (two * C4[:, None] * invrs_tau_C4_zm[:, 1:-1]) / three
    )
    return lhs_pr1_wp2


def wp2_terms_bp_pr2_rhs(nzm, ngrdcol, gr, C_uu_buoy, thv_ds_zm, wpthvp):
    """Buoyancy production of w'^2 and w'^2 pressure term 2:  explicit portion of
    the code.

    The d(w'^2)/dt equation contains a buoyancy production term:

    + 2 (g/thv_ds) w'th_v';

    and pressure term 2:

    - C_5 ( -2 w'^2 dw/dz + 2 (g/thv_ds) w'th_v' ).

    Note:  We have broken C5 up into C_uu_shr for the accumulation term
           and C_uu_buoy for the buoyancy term.
    """
    del gr
    rhs_bp_pr2_wp2 = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    # Calculate term at all interior grid levels.
    rhs_bp_pr2_wp2 = rhs_bp_pr2_wp2.at[:, 1:-1].set(
        (one - C_uu_buoy[:, None]) * two
        * (grav / thv_ds_zm[:, 1:-1]) * wpthvp[:, 1:-1]
    )
    return rhs_bp_pr2_wp2


def wp2_term_dp1_rhs(
    nzm, ngrdcol, gr, C1_Skw_fnc,
    invrs_tau1m, threshold, up2, vp2,
    l_damp_wp2_using_em,
):
    """Dissipation term 1 for w'^2:  explicit portion of the code.

    When l_damp_wp2_using_em == .false., then
    Dissipation term 1 for w'^2:  explicit portion of the code.

    if l_damp_wp2_using_em == .true., then
    we damp wp2 using a more standard turbulence closure, -(2/3)*em/tau
    This only works if C1=C14 and l_stability_correct_tau_zm =.false.
    A factor of (1/3) is absorbed into C1.
    The threshold is implicitly set to 0.
    """
    del gr
    rhs_dp1_wp2 = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    # Calculate term at all interior grid levels.
    if l_damp_wp2_using_em:
        rhs_dp1_wp2 = rhs_dp1_wp2.at[:, 1:-1].set(
            -(C1_Skw_fnc[:, 1:-1] * invrs_tau1m[:, 1:-1])
            * (up2[:, 1:-1] + vp2[:, 1:-1])
        )
    else:
        rhs_dp1_wp2 = rhs_dp1_wp2.at[:, 1:-1].set(
            (C1_Skw_fnc[:, 1:-1] * invrs_tau1m[:, 1:-1]) * threshold
        )
    return rhs_dp1_wp2


def wp2_term_pr3_rhs(
    nzm, nzt, ngrdcol, gr,
    C_uu_shr,
    C_uu_buoy,
    thv_ds_zm, wpthvp, upwp,
    um, vpwp, vm,
):
    """Pressure term 3 for w'^2:  explicit portion of the code.

    The d(w'^2)/dt equation contains pressure term 3:

    + (2/3) C_5 [ (g/thv_ds) w'th_v' - u'w' du/dz - v'w' dv/dz ].

    Note that below we have broken up C5 into C_uu_shr for shear terms and
    C_uu_buoy for buoyancy terms.
    """
    del nzt
    rhs_pr3_wp2 = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    # Michael Falk, 2 August 2007
    # Use the following code for standard mixing, with c_k=0.548:
    rhs_pr3_wp2 = rhs_pr3_wp2.at[:, 1:-1].set(
        two_thirds
        * (
            C_uu_buoy[:, None]
            * (grav / thv_ds_zm[:, 1:-1]) * wpthvp[:, 1:-1]
            + C_uu_shr[:, None]
            * (
                -upwp[:, 1:-1] * gr.invrs_dzm[:, 1:-1] * (um[:, 1:] - um[:, :-1])
                -vpwp[:, 1:-1] * gr.invrs_dzm[:, 1:-1] * (vm[:, 1:] - vm[:, :-1])
            )
        )
    )
    # Added by dschanen for ticket #36
    # We have found that when shear generation is zero this term will only be
    # offset by hole-filling (wp2_pd) and reduces turbulence
    # unrealistically at lower altitudes to make up the difference.
    rhs_pr3_wp2 = rhs_pr3_wp2.at[:, 1:-1].set(
        jnp.maximum(rhs_pr3_wp2[:, 1:-1], zero_threshold)
    )
    return rhs_pr3_wp2


def wp2_term_pr1_rhs(nzm, ngrdcol, gr, C4, up2, vp2, invrs_tau_C4_zm):
    """Pressure term 1 for w'^2:  explicit portion of the code.

    Pressure term 1 has both implicit and explicit components.
    The explicit portion is:

    + ( C_4 / tau_1m ) * (1/3) * ( u'^2 + v'^2 );

    and is computed in this function.
    """
    del gr
    rhs_pr1_wp2 = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    # Calculate term at all interior grid levels.
    rhs_pr1_wp2 = rhs_pr1_wp2.at[:, 1:-1].set(
        (C4[:, None] * (up2[:, 1:-1] + vp2[:, 1:-1])
         * invrs_tau_C4_zm[:, 1:-1]) / three
    )
    return rhs_pr1_wp2


def wp2_term_pr_dfsn_rhs(
    nzm, nzt, ngrdcol, gr, C_wp2_pr_dfsn,
    rho_ds_zt, invrs_rho_ds_zm,
    wpup2, wpvp2, wp3,
):
    """Pressure-diffusion RHS for w'^2.

    This term is intended to represent the "diffusion" part of the wp2
    pressure correlation.  The total pressure diffusion term,

      -1 / rho * ( d( <u_k'p'> )/dx_i + d( <u_i'p'> )/dx_k )

    becomes

      -2 / rho * d( <w'p'> )/dz

    for the w'^2 equation.

    References:
      Lumley 1978, p. 170.  See eq. 6.47 and accompanying discussion.
    """
    del nzt
    wpuip2 = wpup2 + wpvp2 + wp3
    rhs_pr_dfsn_wp2 = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    rhs_pr_dfsn_wp2 = rhs_pr_dfsn_wp2.at[:, 1:-1].set(
        C_wp2_pr_dfsn[:, None] * invrs_rho_ds_zm[:, 1:-1]
        * gr.invrs_dzm[:, 1:-1]
        * (
            rho_ds_zt[:, 1:] * wpuip2[:, 1:]
            - rho_ds_zt[:, :-1] * wpuip2[:, :-1]
        )
    )
    if gr.grid_dir_indx > 0:
        # Set lower boundary condition
        rhs_pr_dfsn_wp2 = rhs_pr_dfsn_wp2.at[:, gr.k_lb_zm].set(
            rhs_pr_dfsn_wp2[:, gr.k_lb_zm + gr.grid_dir_indx]
        )
    else:
        # Set lower boundary condition
        rhs_pr_dfsn_wp2 = rhs_pr_dfsn_wp2.at[:, gr.k_lb_zm].set(
            rhs_pr_dfsn_wp2[:, gr.k_lb_zm + gr.grid_dir_indx]
        )
    # Set upper boundary to 0
    rhs_pr_dfsn_wp2 = rhs_pr_dfsn_wp2.at[:, gr.k_ub_zm].set(zero)
    return rhs_pr_dfsn_wp2


def wp3_term_ta_new_pdf_lhs(
    nzm, nzt, ngrdcol, gr, coef_wp4_implicit,
    wp2, rho_ds_zm, invrs_rho_ds_zt,
):
    """Turbulent advection of <w'^3>:  implicit portion of the code.

    This implicit discretization is specifically for the new PDF.

    The d<w'^3>/dt equation contains a turbulent advection term:

    - (1/rho_ds) * d( rho_ds * <w'^4> )/dz.

    A substitution, which is specific to the new PDF, is made in order to
    close the turbulent advection term, such that:

    <w'^4> = coef_wp4_implicit * <w'^2>^2.

    invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) )
    """
    del nzm
    lhs_ta_wp3 = jnp.zeros((ndiags2, ngrdcol, nzt), dtype=jnp.float64)
    # Momentum superdiagonal: [ x wp2(k+1,<t+1>) ]
    lhs_ta_wp3 = lhs_ta_wp3.at[0, :, 1:-1].set(
        invrs_rho_ds_zt[:, 1:-1] * gr.invrs_dzt[:, 1:-1]
        * rho_ds_zm[:, 2:nzt] * coef_wp4_implicit[:, 2:nzt] * wp2[:, 2:nzt]
    )
    # Momentum subdiagonal: [ x wp2(k,<t+1>) ]
    lhs_ta_wp3 = lhs_ta_wp3.at[1, :, 1:-1].set(
        -invrs_rho_ds_zt[:, 1:-1] * gr.invrs_dzt[:, 1:-1]
        * rho_ds_zm[:, 1:nzt - 1] * coef_wp4_implicit[:, 1:nzt - 1] * wp2[:, 1:nzt - 1]
    )
    return lhs_ta_wp3


def wp3_term_ta_ADG1_lhs(
    nzm, nzt, ngrdcol, gr,
    wp2, a1_coef, a1_coef_zt,
    a3_coef, a3_coef_zt,
    wp3_on_wp2, rho_ds_zm,
    rho_ds_zt, invrs_rho_ds_zt,
    l_standard_term_ta,
    l_partial_upwind_wp3,
):
    """Turbulent advection of w'^3:  implicit portion of the code.

    This implicit discretization is specifically for the ADG1 PDF.

    The d(w'^3)/dt equation contains a turbulent advection term:

    - (1/rho_ds) * d( rho_ds * w'^4 )/dz.

    A substitution, which is specific to ADG1, is made in order to close the
    turbulent advection term, such that:

    w'^4 = a_3 * (w'^2)^2  +  a_1 * ( (w'^3)^2 / w'^2 );

    where both a_1 and a_3 are variables that are functions of sigma_sqd_w,
    such that:

    a_1 = 1 / (1 - sigma_sqd_w); and

    a_3 = 3*(sigma_sqd_w)^2 + 6*(1 - sigma_sqd_w)*sigma_sqd_w
          + (1 - sigma_sqd_w)^2.

    https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:wp4_diagnosis
    """
    del nzm
    lhs_ta_wp3 = jnp.zeros((ndiags5, ngrdcol, nzt), dtype=jnp.float64)
    inv = invrs_rho_ds_zt[:, 1:-1]
    idzt = gr.invrs_dzt[:, 1:-1]

    if l_standard_term_ta:
        # The turbulent advection term is discretized normally, in accordance
        # with the model equations found in the documentation and the description
        # listed above.
        if not l_partial_upwind_wp3:
            # All portions of the wp3 turbulent advection term for ADG1 use
            # centered discretization in accordance with description and diagram
            # shown above.
            # Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
            lhs_ta_wp3 = lhs_ta_wp3.at[0, :, 1:-1].set(
                inv * idzt
                * rho_ds_zm[:, 2:nzt] * a1_coef[:, 2:nzt] * wp3_on_wp2[:, 2:nzt]
                * gr.weights_zt2zm[:, 2:nzt, T_ABOVE]
            )
            # Momentum superdiagonal: [ x wp2(k+1,<t+1>) ]
            lhs_ta_wp3 = lhs_ta_wp3.at[1, :, 1:-1].set(
                inv * idzt
                * rho_ds_zm[:, 2:nzt] * a3_coef[:, 2:nzt] * wp2[:, 2:nzt]
            )
            # Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
            lhs_ta_wp3 = lhs_ta_wp3.at[2, :, 1:-1].set(
                inv * idzt
                * (
                    rho_ds_zm[:, 2:nzt] * a1_coef[:, 2:nzt] * wp3_on_wp2[:, 2:nzt]
                    * gr.weights_zt2zm[:, 2:nzt, T_BELOW]
                    - rho_ds_zm[:, 1:nzt - 1] * a1_coef[:, 1:nzt - 1]
                    * wp3_on_wp2[:, 1:nzt - 1]
                    * gr.weights_zt2zm[:, 1:nzt - 1, T_ABOVE]
                )
            )
            # Momentum subdiagonal: [ x wp2(k,<t+1>) ]
            lhs_ta_wp3 = lhs_ta_wp3.at[3, :, 1:-1].set(
                -inv * idzt
                * rho_ds_zm[:, 1:nzt - 1] * a3_coef[:, 1:nzt - 1]
                * wp2[:, 1:nzt - 1]
            )
            # Thermodynamic subdiagonal: [ x wp3(k-1,<t+1>) ]
            lhs_ta_wp3 = lhs_ta_wp3.at[4, :, 1:-1].set(
                -inv * idzt
                * rho_ds_zm[:, 1:nzt - 1] * a1_coef[:, 1:nzt - 1]
                * wp3_on_wp2[:, 1:nzt - 1]
                * gr.weights_zt2zm[:, 1:nzt - 1, T_BELOW]
            )
        else:
            # Partial upwinding of the wp3 turbulent advection term, where the
            # portion of the wp3 turbulent advection term that is linearized in
            # terms of wp2<t+1> is still handled using centered discretization,
            # but the portion of the term that is linearized in terms of wp3<t+1>
            # is handled using an "upwind" discretization that also takes into
            # "winds" that converge or diverge around the central thermodynamic
            # grid level.  Provided by Chris Vogl and Shixuan Zhang.
            grid_dir = gr.grid_dir
            # Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
            lhs_ta_wp3 = lhs_ta_wp3.at[0, :, 1:-1].set(
                inv * idzt * rho_ds_zt[:, 2:nzt] * grid_dir
                * jnp.minimum(
                    grid_dir * a1_coef[:, 2:nzt] * wp3_on_wp2[:, 2:nzt],
                    zero,
                )
            )
            # Momentum superdiagonal: [ x wp2(k+1,<t+1>) ]
            lhs_ta_wp3 = lhs_ta_wp3.at[1, :, 1:-1].set(
                inv * idzt
                * rho_ds_zm[:, 2:nzt] * a3_coef[:, 2:nzt] * wp2[:, 2:nzt]
            )
            # Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
            lhs_ta_wp3 = lhs_ta_wp3.at[2, :, 1:-1].set(
                inv * idzt * rho_ds_zt[:, 1:-1] * grid_dir
                * (
                    jnp.maximum(
                        grid_dir * a1_coef[:, 2:nzt] * wp3_on_wp2[:, 2:nzt],
                        zero,
                    )
                    - jnp.minimum(
                        grid_dir * a1_coef[:, 1:nzt - 1]
                        * wp3_on_wp2[:, 1:nzt - 1],
                        zero,
                    )
                )
            )
            # Momentum subdiagonal: [ x wp2(k,<t+1>) ]
            lhs_ta_wp3 = lhs_ta_wp3.at[3, :, 1:-1].set(
                -inv * idzt
                * rho_ds_zm[:, 1:nzt - 1] * a3_coef[:, 1:nzt - 1]
                * wp2[:, 1:nzt - 1]
            )
            # Thermodynamic subdiagonal: [ x wp3(k-1,<t+1>) ]
            lhs_ta_wp3 = lhs_ta_wp3.at[4, :, 1:-1].set(
                -inv * idzt * rho_ds_zt[:, :nzt - 2] * grid_dir
                * jnp.maximum(
                    grid_dir * a1_coef[:, 1:nzt - 1]
                    * wp3_on_wp2[:, 1:nzt - 1],
                    zero,
                )
            )
    else:
        # Alternate discretization for the turbulent advection term, which
        # contains the term:
        #  - (1/rho_ds) * d [ rho_ds * a_1 * (w'^3)^2 / w'^2 ] / dz.  In order
        # to help stabilize w'^3, a_1 has been pulled outside of the derivative.
        # On the left-hand side of the equation, this effects the thermodynamic
        # superdiagonal (kp1_tdiag), the thermodynamic main diagonal (k_tdiag),
        # and the thermodynamic subdiagonal (km1_tdiag).
        #
        # Additionally, the discretization of the turbulent advection term, which
        # contains the term:
        #  - (1/rho_ds) * d [ rho_ds * a_3 * (w'^2)^2 ] / dz, has been altered to
        # pull a_3 outside of the derivative.  This was done in order to help
        # stabilize w'^3.  On the left-hand side of the equation, this effects
        # the momentum superdiagonal (k_mdiag) and the momentum subdiagonal
        # (km1_mdiag).
        # Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
        lhs_ta_wp3 = lhs_ta_wp3.at[0, :, 1:-1].set(
            inv * a1_coef_zt[:, 1:-1] * idzt
            * rho_ds_zm[:, 2:nzt] * wp3_on_wp2[:, 2:nzt]
            * gr.weights_zt2zm[:, 2:nzt, T_ABOVE]
        )
        # Momentum superdiagonal: [ x wp2(k+1,<t+1>) ]
        lhs_ta_wp3 = lhs_ta_wp3.at[1, :, 1:-1].set(
            inv * a3_coef_zt[:, 1:-1] * idzt
            * rho_ds_zm[:, 2:nzt] * wp2[:, 2:nzt]
        )
        # Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
        lhs_ta_wp3 = lhs_ta_wp3.at[2, :, 1:-1].set(
            inv * a1_coef_zt[:, 1:-1] * idzt
            * (
                rho_ds_zm[:, 2:nzt] * wp3_on_wp2[:, 2:nzt]
                * gr.weights_zt2zm[:, 2:nzt, T_BELOW]
                - rho_ds_zm[:, 1:nzt - 1] * wp3_on_wp2[:, 1:nzt - 1]
                * gr.weights_zt2zm[:, 1:nzt - 1, T_ABOVE]
            )
        )
        # Momentum subdiagonal: [ x wp2(k,<t+1>) ]
        lhs_ta_wp3 = lhs_ta_wp3.at[3, :, 1:-1].set(
            -inv * a3_coef_zt[:, 1:-1] * idzt
            * rho_ds_zm[:, 1:nzt - 1] * wp2[:, 1:nzt - 1]
        )
        # Thermodynamic subdiagonal: [ x wp3(k-1,<t+1>) ]
        lhs_ta_wp3 = lhs_ta_wp3.at[4, :, 1:-1].set(
            -inv * a1_coef_zt[:, 1:-1] * idzt
            * rho_ds_zm[:, 1:nzt - 1] * wp3_on_wp2[:, 1:nzt - 1]
            * gr.weights_zt2zm[:, 1:nzt - 1, T_BELOW]
        )

    return lhs_ta_wp3


def wp3_term_tp_lhs(nzm, nzt, ngrdcol, gr, coef_wp3_tp, wp2, rho_ds_zm, invrs_rho_ds_zt):
    """Turbulent production of w'^3:  implicit portion of the code.

    The d(w'^3)/dt equation contains a turbulent production term:

    + 3 w'^2 dw'^2/dz.

    The pressure scrambling terms damp this term by the coefficient
    coef_wp3_tp.
    """
    del nzm
    lhs_tp_wp3 = jnp.zeros((ndiags2, ngrdcol, nzt), dtype=jnp.float64)
    # Momentum superdiagonal: [ x wp2(k+1,<t+1>) ]
    lhs_tp_wp3 = lhs_tp_wp3.at[0, :, 1:-1].set(
        -coef_wp3_tp[:, None] * three * invrs_rho_ds_zt[:, 1:-1]
        * gr.invrs_dzt[:, 1:-1]
        * rho_ds_zm[:, 2:nzt] * wp2[:, 2:nzt]
        + coef_wp3_tp[:, None] * three_halves * gr.invrs_dzt[:, 1:-1]
        * wp2[:, 2:nzt]
    )
    # Momentum subdiagonal: [ x wp2(k,<t+1>) ]
    lhs_tp_wp3 = lhs_tp_wp3.at[1, :, 1:-1].set(
        coef_wp3_tp[:, None] * three * invrs_rho_ds_zt[:, 1:-1]
        * gr.invrs_dzt[:, 1:-1]
        * rho_ds_zm[:, 1:nzt - 1] * wp2[:, 1:nzt - 1]
        - coef_wp3_tp[:, None] * three_halves * gr.invrs_dzt[:, 1:-1]
        * wp2[:, 1:nzt - 1]
    )
    return lhs_tp_wp3


def wp3_terms_ac_pr2_lhs(nzm, nzt, ngrdcol, gr, C11_Skw_fnc, wm_zm):
    """Accumulation of w'^3 and w'^3 pressure term 2:  implicit portion of the
    code.

    The w'^3 accumulation term is completely implicit, while w'^3 pressure
    term 2 has both implicit and explicit components.
    """
    del nzm
    lhs_ac_pr2_wp3 = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    # Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
    lhs_ac_pr2_wp3 = lhs_ac_pr2_wp3.at[:, 1:-1].set(
        (one - C11_Skw_fnc[:, 1:-1])
        * three * gr.invrs_dzt[:, 1:-1]
        * (wm_zm[:, 2:nzt] - wm_zm[:, 1:nzt - 1])
    )
    return lhs_ac_pr2_wp3


def wp3_term_pr1_lhs(
    nzt, ngrdcol, gr,
    C8, C8b,
    invrs_tau_wp3_zt, Skw_zt,
    l_damp_wp3_Skw_squared,
):
    """Pressure term 1 for w'^3:  implicit portion of the code.

    Pressure term 1 is the term:

    - (C_8/tau_w3t) * ( C_8b * Sk_wt^2 + 1 ) * w'^3;

    where Sk_wt = w'^3 / (w'^2)^(3/2).

    A Taylor Series expansion (truncated after the first derivative term) of
    L(w'^3) around w'^3 = w'^3(t) is used to linearize pressure term 1.
    """
    del gr
    lhs_pr1_wp3 = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    # Calculate term at all interior grid levels.
    if l_damp_wp3_Skw_squared:
        lhs_pr1_wp3 = lhs_pr1_wp3.at[:, 1:-1].set(
            (C8[:, None] * invrs_tau_wp3_zt[:, 1:-1])
            * (three * C8b[:, None] * Skw_zt[:, 1:-1] ** 2 + one)
        )
    else:
        lhs_pr1_wp3 = lhs_pr1_wp3.at[:, 1:-1].set(
            (C8[:, None] * invrs_tau_wp3_zt[:, 1:-1])
            * (five * C8b[:, None] * Skw_zt[:, 1:-1] ** 4 + one)
        )
    return lhs_pr1_wp3


def wp3_term_ta_explicit_rhs(nzm, nzt, ngrdcol, gr, wp4, rho_ds_zm, invrs_rho_ds_zt):
    """Turbulent advection of <w'^3>:  explicit portion of the code.

    This explicit discretization works generally for any PDF.

    The d<w'^3>/dt equation contains a turbulent advection term:

    - (1/rho_ds) * d( rho_ds * <w'^4> )/dz.
    """
    del nzm
    rhs_ta_wp3 = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    # Calculate term at all interior grid levels.
    rhs_ta_wp3 = rhs_ta_wp3.at[:, 1:-1].set(
        -invrs_rho_ds_zt[:, 1:-1] * gr.invrs_dzt[:, 1:-1]
        * (
            rho_ds_zm[:, 2:nzt] * wp4[:, 2:nzt]
            - rho_ds_zm[:, 1:nzt - 1] * wp4[:, 1:nzt - 1]
        )
    )
    return rhs_ta_wp3


def wp3_terms_bp1_pr2_rhs(nzt, ngrdcol, gr, C11_Skw_fnc, thv_ds_zt, wp2thvp):
    """Buoyancy production of w'^3 and w'^3 pressure term 2:  explicit portion of
    the code.

    The d(w'^3)/dt equation contains a buoyancy production term:

    + 3 (g/thv_ds) w'^2th_v';

    and pressure term 2:

    - C_11 ( -3 w'^3 dw/dz + 3 (g/thv_ds) w'^2th_v' ).
    """
    del gr
    rhs_bp1_pr2_wp3 = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    # Calculate term at all interior grid levels.
    rhs_bp1_pr2_wp3 = rhs_bp1_pr2_wp3.at[:, 1:-1].set(
        (one - C11_Skw_fnc[:, 1:-1])
        * three * (grav / thv_ds_zt[:, 1:-1]) * wp2thvp[:, 1:-1]
    )
    return rhs_bp1_pr2_wp3


def wp3_term_pr_turb_rhs(
    nzm, nzt, ngrdcol, gr, C_wp3_pr_turb,
    Kh_zt, wpthvp, dum_dz, dvm_dz,
    upwp, vpwp,
    thv_ds_zt,
    rho_ds_zm, invrs_rho_ds_zt,
    em, wp2,
    l_use_tke_in_wp3_pr_turb_term,
):
    """Pressure-turbulence correlation RHS for w'^3.

    Experimental term from CLUBB TRAC ticket #411. The derivative here is of
    the form used to match LES data.

    This does not appear in Andre et al. 1976 or Bougeault et al. 1981, but
    is based on experiments in matching LES data.

    References:
      None
    """
    del nzm
    rhs_pr_turb_wp3 = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    if not l_use_tke_in_wp3_pr_turb_term:
        rhs_pr_turb_wp3 = rhs_pr_turb_wp3.at[:, 1:-1].set(
            -C_wp3_pr_turb[:, None] * Kh_zt[:, 1:-1] * gr.invrs_dzt[:, 1:-1]
            * (
                grav / thv_ds_zt[:, 1:-1]
                * (wpthvp[:, 2:nzt] - wpthvp[:, 1:nzt - 1])
                - (
                    upwp[:, 2:nzt] * dum_dz[:, 2:nzt]
                    - upwp[:, 1:nzt - 1] * dum_dz[:, 1:nzt - 1]
                )
                - (
                    vpwp[:, 2:nzt] * dvm_dz[:, 2:nzt]
                    - vpwp[:, 1:nzt - 1] * dvm_dz[:, 1:nzt - 1]
                )
            )
        )
    else:
        rhs_pr_turb_wp3 = rhs_pr_turb_wp3.at[:, 1:-1].set(
            -C_wp3_pr_turb[:, None] * invrs_rho_ds_zt[:, 1:-1]
            * gr.invrs_dzt[:, 1:-1]
            * (
                rho_ds_zm[:, 2:nzt] * wp2[:, 2:nzt] * em[:, 2:nzt]
                - rho_ds_zm[:, 1:nzt - 1] * wp2[:, 1:nzt - 1] * em[:, 1:nzt - 1]
            )
        )
    return rhs_pr_turb_wp3


def wp3_term_pr_dfsn_rhs(
    nzm, nzt, ngrdcol, gr, C_wp3_pr_dfsn,
    rho_ds_zm, invrs_rho_ds_zt,
    wp2up2, wp2vp2, wp4,
    up2, vp2, wp2,
):
    """Pressure-diffusion RHS for w'^3.

    This term is intended to represent the "diffusion" part of the total wp3
    pressure correlation.  The total wp3 pressure term, -3w'^2/rho*dp'/dz, can be
    split into

      -3w'^2/rho*dp'/dz = + 3p'/rho*d(w'^2)/dz - 3/rho*d(w'^2p')/dz

    using the product rule.  The second term on the RHS we consider to be the
    diffusion part, calculated by this subroutine.

    References:
      Lumley 1978, p. 170.  See eq. 6.47 and accompanying discussion.
    """
    del nzm
    wp2uip2 = wp2up2 + wp2vp2 + wp4
    wp2_uip2 = wp2 * up2 + wp2 * vp2 + wp2 * wp2
    rhs_pr_dfsn_wp3 = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    rhs_pr_dfsn_wp3 = rhs_pr_dfsn_wp3.at[:, 1:-1].set(
        C_wp3_pr_dfsn[:, None] * invrs_rho_ds_zt[:, 1:-1]
        * gr.invrs_dzt[:, 1:-1]
        * (
            rho_ds_zm[:, 2:nzt] * (wp2uip2[:, 2:nzt] - wp2_uip2[:, 2:nzt])
            - rho_ds_zm[:, 1:nzt - 1]
            * (wp2uip2[:, 1:nzt - 1] - wp2_uip2[:, 1:nzt - 1])
        )
    )
    return rhs_pr_dfsn_wp3


def wp3_term_pr1_rhs(
    nzt, ngrdcol, gr,
    C8, C8b,
    invrs_tau_wp3_zt, Skw_zt, wp3,
    l_damp_wp3_Skw_squared,
):
    """Pressure term 1 for w'^3:  explicit portion of the code.

    Pressure term 1 is the term:

    - (C_8/tau_w3t) * ( C_8b * Sk_wt^2 + 1 ) * w'^3;

    where Sk_wt = w'^3 / (w'^2)^(3/2).

    The explicit portion is:

    + (C_8/tau_w3t) * ( 2 * C_8b * Sk_wt^2 + 1 ) * w'^3(t).
    """
    del gr
    rhs_pr1_wp3 = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    # Calculate term at all interior grid levels.
    if l_damp_wp3_Skw_squared:
        rhs_pr1_wp3 = rhs_pr1_wp3.at[:, 1:-1].set(
            (C8[:, None] * invrs_tau_wp3_zt[:, 1:-1])
            * (two * C8b[:, None] * Skw_zt[:, 1:-1] ** 2) * wp3[:, 1:-1]
        )
    else:
        rhs_pr1_wp3 = rhs_pr1_wp3.at[:, 1:-1].set(
            (C8[:, None] * invrs_tau_wp3_zt[:, 1:-1])
            * (four * C8b[:, None] * Skw_zt[:, 1:-1] ** 4) * wp3[:, 1:-1]
        )
    return rhs_pr1_wp3


__all__ = [
    "advance_wp2_wp3",
    "wp23_solve",
    "wp23_lhs",
    "wp23_rhs",
    "wp2_term_ta_lhs",
    "wp2_terms_ac_pr2_lhs",
    "wp2_term_dp1_lhs",
    "wp2_term_pr1_lhs",
    "wp2_terms_bp_pr2_rhs",
    "wp2_term_dp1_rhs",
    "wp2_term_pr3_rhs",
    "wp2_term_pr1_rhs",
    "wp2_term_pr_dfsn_rhs",
    "wp3_term_ta_new_pdf_lhs",
    "wp3_term_ta_ADG1_lhs",
    "wp3_term_tp_lhs",
    "wp3_terms_ac_pr2_lhs",
    "wp3_term_pr1_lhs",
    "wp3_term_ta_explicit_rhs",
    "wp3_terms_bp1_pr2_rhs",
    "wp3_term_pr_turb_rhs",
    "wp3_term_pr_dfsn_rhs",
    "wp3_term_pr1_rhs",
]
