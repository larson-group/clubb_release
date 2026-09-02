"""JAX port of ``src/CLUBB_core/advance_xp3_module.F90``.

Description:
  Predicts the value of <x'^3> for <rt'^3>, <thl'^3>, and <sclr'^3>.

References:

Porting deviations:
- Stats are threaded explicitly with JaxStats because the Fortran routine uses
  global stats side effects that are not JAX-compatible state.
- Intent(inout) moments are returned explicitly rather than mutated by Fortran
  argument side effects.
"""

from functools import partial

import jax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.Skx_module import Skx_func, xp3_LG_2005_ansatz
from clubb_jax.src.CLUBB_core.constants_clubb import (
    grav,
    one,
    one_half,
    rt_tol,
    thl_tol,
    three,
    w_tol,
    w_tol_sqd,
    zero,
    zero_threshold,
)
from clubb_jax.src.CLUBB_core.model_flags import iiPDF_ADG1
from clubb_jax.src.CLUBB_core.parameter_indices import ixp3_coef_base, ixp3_coef_slope
from clubb_jax.src.CLUBB_core.grid_class import ddzm, zm2zt, zt2zm
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core import Grid


xp3_rtp3 = 1  # Named constant for solving rtp3
xp3_thlp3 = 2  # Named constant for solving thlp3
xp3_sclrp3 = 3  # Named constant for solving sclrp3


@partial(
    jax.jit,
    static_argnames=("nzm", "nzt", "ngrdcol", "sclr_dim", "l_lmm_stepping"),
)
def advance_xp3(
    nzm: int, nzt: int, ngrdcol: int, sclr_dim: int, sclr_tol, gr: Grid, dt,
    rtm, thlm, rtp2, thlp2, wprtp,
    wpthlp, wprtp2, wpthlp2, rho_ds_zm,
    invrs_rho_ds_zt, invrs_tau_zt, tau_max_zt,
    sclrm, sclrp2, wpsclrp, wpsclrp2,
    wp2, wp3, upwp, vpwp, up2, vp2,
    thvm, clubb_params,
    l_lmm_stepping: bool,
    stats: JaxStats,
    rtp3, thlp3, sclrp3, up3, vp3,
):
    """Advance <rt'^3>, <thl'^3>, and <sclr'^3> one model timestep using a
    simplified form of the <x'^3> predictive equation.  The simplified <x'^3>
    equation can either be advanced from its previous value or calculated
    using a steady-state approximation.
    """

    wp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, wp2, w_tol_sqd)

    # Advance <rt'^3> one model timestep or calculate <rt'^3> using a
    # steady-state approximation.
    rtp3, stats = advance_xp3_simplified(
        nzm, nzt, ngrdcol, gr, xp3_rtp3, dt,
        rtm, rtp2, wprtp,
        wprtp2, rho_ds_zm,
        invrs_rho_ds_zt,
        invrs_tau_zt, tau_max_zt,
        rt_tol, l_lmm_stepping,
        stats,
        rtp3,
    )

    # Advance <thl'^3> one model timestep or calculate <thl'^3> using a
    # steady-state approximation.
    thlp3, stats = advance_xp3_simplified(
        nzm, nzt, ngrdcol, gr, xp3_thlp3, dt,
        thlm, thlp2, wpthlp,
        wpthlp2, rho_ds_zm,
        invrs_rho_ds_zt,
        invrs_tau_zt, tau_max_zt,
        thl_tol, l_lmm_stepping,
        stats,
        thlp3,
    )

    # Advance <sclr'^3> one model timestep or calculate <sclr'^3> using a
    # steady-state approximation.
    for sclr in range(sclr_dim):
        sclrp3_sclr, stats = advance_xp3_simplified(
            nzm, nzt, ngrdcol, gr, xp3_sclrp3, dt,
            sclrm[:, :, sclr], sclrp2[:, :, sclr], wpsclrp[:, :, sclr],
            wpsclrp2[:, :, sclr], rho_ds_zm,
            invrs_rho_ds_zt,
            invrs_tau_zt, tau_max_zt,
            sclr_tol[sclr], l_lmm_stepping,
            stats,
            sclrp3[:, :, sclr],
        )
        sclrp3 = sclrp3.at[:, :, sclr].set(sclrp3_sclr)

    # Use a modified form of the Larson and Golaz (2005) ansatz for the
    # ADG1 PDF to calculate <u'^3> and <v'^3> for another type of PDF.
    Skw_zt = Skx_func(nzt, ngrdcol, wp2_zt, wp3, w_tol, clubb_params)

    wpthlp_zt = zm2zt(nzm, nzt, ngrdcol, gr, wpthlp)
    wprtp_zt = zm2zt(nzm, nzt, ngrdcol, gr, wprtp)
    # Positive def. quantity
    thlp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, thlp2, thl_tol ** 2)
    # Positive def. quantity
    rtp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, rtp2, rt_tol ** 2)

    upwp_zt = zm2zt(nzm, nzt, ngrdcol, gr, upwp)
    vpwp_zt = zm2zt(nzm, nzt, ngrdcol, gr, vpwp)
    # Positive def. quantity
    up2_zt = zm2zt(nzm, nzt, ngrdcol, gr, up2, w_tol_sqd)
    # Positive def. quantity
    vp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, vp2, w_tol_sqd)

    thvm_zm = zt2zm(nzm, nzt, ngrdcol, gr, thvm, zero_threshold)
    ddzm_thvm_zm = ddzm(nzm, nzt, ngrdcol, gr, thvm_zm)
    brunt_vaisala_freq_sqd_zt = jnp.maximum(
        (grav / thvm) * ddzm_thvm_zm,
        zero,
    )

    # The xp3_coef_fnc provides some extra tunability to the simple xp3
    # equation.
    # When xp3_coef_fnc goes to 0, the value of Skx goes to the smallest
    # magnitude permitted by the function.  When xp3_coef_fnc goes to 1, the
    # magnitude of Skx becomes huge.
    xp3_coef_fnc = (
        clubb_params[:, ixp3_coef_base][:, None]
        + (one - clubb_params[:, ixp3_coef_slope])[:, None]
        * (
            one
            - jnp.exp(
                brunt_vaisala_freq_sqd_zt
                / clubb_params[:, ixp3_coef_slope][:, None]
            )
        )
    )

    up3 = xp3_LG_2005_ansatz(
        nzt, ngrdcol, Skw_zt, upwp_zt, wp2_zt,
        up2_zt, xp3_coef_fnc,
        clubb_params, w_tol,
    )

    vp3 = xp3_LG_2005_ansatz(
        nzt, ngrdcol, Skw_zt, vpwp_zt, wp2_zt,
        vp2_zt, xp3_coef_fnc,
        clubb_params, w_tol,
    )

    if stats.l_sample:
        stats = stats.update("thlp2_zt", thlp2_zt)
        stats = stats.update("wpthlp_zt", wpthlp_zt)
        stats = stats.update("wprtp_zt", wprtp_zt)
        stats = stats.update("rtp2_zt", rtp2_zt)
        stats = stats.update("up2_zt", up2_zt)
        stats = stats.update("vp2_zt", vp2_zt)
        stats = stats.update("upwp_zt", upwp_zt)
        stats = stats.update("vpwp_zt", vpwp_zt)

    return rtp3, thlp3, sclrp3, up3, vp3, stats


@partial(
    jax.jit,
    static_argnames=("nzm", "nzt", "ngrdcol", "sclr_dim", "iiPDF_type"),
)
def diagnose_xp3(
    nzm: int, nzt: int, ngrdcol: int, sclr_dim: int, sclr_tol, gr: Grid,
    iiPDF_type: int, clubb_params,
    wp2, wp3, thvm,
    wprtp, wpthlp, rtp2, thlp2, upwp, vpwp, up2, vp2,
    sigma_sqd_w, wpsclrp, sclrp2,
    stats: JaxStats,
    rtp3, thlp3, up3, vp3,
    sclrp3,
):
    """Diagnose third-order scalar and horizontal-velocity moments from the
    current second-order moments and PDF width information.

    References:
      Larson and Golaz (2005) for the diagnostic ansatz.
    """

    wp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, wp2, w_tol_sqd)

    # The ADG1 PDF must use this option.
    Skw_zt = Skx_func(nzt, ngrdcol, wp2_zt, wp3, w_tol, clubb_params)

    wpthlp_zt = zm2zt(nzm, nzt, ngrdcol, gr, wpthlp)
    wprtp_zt = zm2zt(nzm, nzt, ngrdcol, gr, wprtp)
    # Positive def. quantity
    thlp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, thlp2, thl_tol ** 2)
    # Positive def. quantity
    rtp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, rtp2, rt_tol ** 2)

    upwp_zt = zm2zt(nzm, nzt, ngrdcol, gr, upwp)
    vpwp_zt = zm2zt(nzm, nzt, ngrdcol, gr, vpwp)
    # Positive def. quantity
    up2_zt = zm2zt(nzm, nzt, ngrdcol, gr, up2, w_tol_sqd)
    # Positive def. quantity
    vp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, vp2, w_tol_sqd)

    if iiPDF_type == iiPDF_ADG1:

        # Use the Larson and Golaz (2005) ansatz for the ADG1 PDF to
        # calculate <rt'^3>, <thl'^3>, <u'^3>, <v'^3>, and <sclr'^3>.
        sigma_sqd_w_zt = zm2zt(
            nzm, nzt, ngrdcol, gr, sigma_sqd_w,
            zero_threshold,
        )

        thlp3 = xp3_LG_2005_ansatz(
            nzt, ngrdcol, Skw_zt, wpthlp_zt, wp2_zt,
            thlp2_zt, sigma_sqd_w_zt,
            clubb_params, thl_tol,
        )

        rtp3 = xp3_LG_2005_ansatz(
            nzt, ngrdcol, Skw_zt, wprtp_zt, wp2_zt,
            rtp2_zt, sigma_sqd_w_zt,
            clubb_params, rt_tol,
        )

        up3 = xp3_LG_2005_ansatz(
            nzt, ngrdcol, Skw_zt, upwp_zt, wp2_zt,
            up2_zt, sigma_sqd_w_zt,
            clubb_params, w_tol,
        )

        vp3 = xp3_LG_2005_ansatz(
            nzt, ngrdcol, Skw_zt, vpwp_zt, wp2_zt,
            vp2_zt, sigma_sqd_w_zt,
            clubb_params, w_tol,
        )

        for sclr in range(sclr_dim):

            wpsclrp_zt = zm2zt(
                nzm, nzt, ngrdcol, gr, wpsclrp[:, :, sclr],
                sclr_tol[sclr] ** 2,
            )
            sclrp2_zt = zm2zt(
                nzm, nzt, ngrdcol, gr, sclrp2[:, :, sclr],
                sclr_tol[sclr] ** 2,
            )

            sclrp3_sclr = xp3_LG_2005_ansatz(
                nzt, ngrdcol, Skw_zt, wpsclrp_zt, wp2_zt,
                sclrp2_zt, sigma_sqd_w_zt,
                clubb_params, sclr_tol[sclr],
            )
            sclrp3 = sclrp3.at[:, :, sclr].set(sclrp3_sclr)

    else:

        # Use a modified form of the Larson and Golaz (2005) ansatz for the
        # ADG1 PDF to calculate <u'^3> and <v'^3> for another type of PDF.
        thvm_zm = zt2zm(nzm, nzt, ngrdcol, gr, thvm, zero_threshold)
        ddzm_thvm_zm = ddzm(nzm, nzt, ngrdcol, gr, thvm_zm)
        brunt_vaisala_freq_sqd_zt = jnp.maximum(
            (grav / thvm) * ddzm_thvm_zm,
            zero,
        )

        # Initialize sigma_sqd_w_zt to zero so we don't break output
        sigma_sqd_w_zt = jnp.zeros((ngrdcol, nzt), dtype=wp2_zt.dtype)

        # The xp3_coef_fnc is used in place of sigma_sqd_w_zt when the
        # ADG1 PDF is not being used.  The xp3_coef_fnc provides some extra
        # tunability to the simple xp3 equation.
        # When xp3_coef_fnc goes to 0, the value of Skx goes to the smallest
        # magnitude permitted by the function.  When xp3_coef_fnc goes to 1,
        # the magnitude of Skx becomes huge.
        # The value of Skx becomes large near cloud top, where there is a
        # higher degree of static stability.  The exp{ } portion of the
        # xp3_coef_fnc allows the xp3_coef_fnc to become larger in regions
        # of high static stability, producing larger magnitude values of Skx.
        xp3_coef_fnc = (
            clubb_params[:, ixp3_coef_base][:, None]
            + (one - clubb_params[:, ixp3_coef_slope])[:, None]
            * (
                one
                - jnp.exp(
                    brunt_vaisala_freq_sqd_zt
                    / clubb_params[:, ixp3_coef_slope][:, None]
                )
            )
        )

        thlp3 = xp3_LG_2005_ansatz(
            nzt, ngrdcol, Skw_zt, wpthlp_zt, wp2_zt,
            thlp2_zt, xp3_coef_fnc,
            clubb_params, thl_tol,
        )

        rtp3 = xp3_LG_2005_ansatz(
            nzt, ngrdcol, Skw_zt, wprtp_zt, wp2_zt,
            rtp2_zt, xp3_coef_fnc,
            clubb_params, rt_tol,
        )

        up3 = xp3_LG_2005_ansatz(
            nzt, ngrdcol, Skw_zt, upwp_zt, wp2_zt,
            up2_zt, xp3_coef_fnc,
            clubb_params, w_tol,
        )

        vp3 = xp3_LG_2005_ansatz(
            nzt, ngrdcol, Skw_zt, vpwp_zt, wp2_zt,
            vp2_zt, xp3_coef_fnc,
            clubb_params, w_tol,
        )

        for sclr in range(sclr_dim):

            wpsclrp_zt = zm2zt(
                nzm, nzt, ngrdcol, gr, wpsclrp[:, :, sclr],
            )
            sclrp2_zt = zm2zt(
                nzm, nzt, ngrdcol, gr, sclrp2[:, :, sclr],
                sclr_tol[sclr] ** 2,
            )

            sclrp3_sclr = xp3_LG_2005_ansatz(
                nzt, ngrdcol, Skw_zt, wpsclrp_zt, wp2_zt,
                sclrp2_zt, xp3_coef_fnc,
                clubb_params, sclr_tol[sclr],
            )
            sclrp3 = sclrp3.at[:, :, sclr].set(sclrp3_sclr)

    if stats.l_sample:
        stats = stats.update("thlp2_zt", thlp2_zt)
        stats = stats.update("wpthlp_zt", wpthlp_zt)
        stats = stats.update("wprtp_zt", wprtp_zt)
        stats = stats.update("rtp2_zt", rtp2_zt)
        stats = stats.update("up2_zt", up2_zt)
        stats = stats.update("vp2_zt", vp2_zt)
        stats = stats.update("upwp_zt", upwp_zt)
        stats = stats.update("vpwp_zt", vpwp_zt)

    return rtp3, thlp3, up3, vp3, sclrp3, stats


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "nzt",
        "ngrdcol",
        "solve_type",
        "l_lmm_stepping",
    ),
)
def advance_xp3_simplified(
    nzm: int, nzt: int, ngrdcol: int, gr: Grid, solve_type: int, dt,
    xm, xp2, wpxp,
    wpxp2, rho_ds_zm,
    invrs_rho_ds_zt,
    invrs_tau_zt, tau_max_zt,
    x_tol, l_lmm_stepping: bool,
    stats: JaxStats,
    xp3,
):
    """Predicts the value of <x'^3> using a simplified form of the <x'^3>
    predictive equation.

    The full predictive equation for <x'^3>, where <x'^3> can be <rt'^3>,
    <thl'^3>, or <sclr'^3>, contains time tendency, mean advection,
    turbulent advection, accumulation, turbulent production, turbulent
    dissipation, diffusion, and microphysics/other forcing.

    The microphysics and turbulent advection terms are both found by
    integration over the subgrid PDF.  This requires new integrated terms.
    The turbulent advection term may need to be made semi-implicit in order
    to aid model stability.  This may be difficult to do for <x'^3>.

    This leaves the following equation:

    d<x'^3>/dt = - 3 * <w'x'^2> * d<x>/dz
                 + 3 * <x'^2> / rho_ds * d( rho_ds * <w'x'> )/dz
                 - ( C_xp3_dissipation / tau ) * <x'^3>.

    When the flag l_predict_xp3 is enabled, the predictive version of <x'^3>
    is used.  When the flag is turned off, the steady-state approximation is
    used.
    """

    # Coefficient in the <x'^3> turbulent dissipation term    [-]
    C_xp3_dissipation = 1.0

    # Flag to either predict <x'^3> or use steady-state approximation.
    l_predict_xp3 = False

    if solve_type == xp3_rtp3:
        # Budget stats for rtp3
        name_bt = "rtp3_bt"
        name_tp = "rtp3_tp"
        name_ac = "rtp3_ac"
        name_dp = "rtp3_dp"
    elif solve_type == xp3_thlp3:
        # Budget stats for thlp3
        name_bt = "thlp3_bt"
        name_tp = "thlp3_tp"
        name_ac = "thlp3_ac"
        name_dp = "thlp3_dp"
    else:
        # Budgets aren't setup for the passive scalars
        name_bt = ""
        name_tp = ""
        name_ac = ""
        name_dp = ""

    if l_predict_xp3:
        if stats.l_sample:
            stats = stats.begin_budget(name_bt, xp3 / dt)

    # Initialize variables
    term_tp = jnp.zeros((ngrdcol, nzt), dtype=xp3.dtype)
    term_ac = jnp.zeros((ngrdcol, nzt), dtype=xp3.dtype)

    # Interpolate <x> to momentum levels.
    xm_zm = zt2zm(nzm, nzt, ngrdcol, gr, xm, zero_threshold)

    # Interpolate <x'^2> to thermodynamic levels.
    xp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, xp2, x_tol ** 2)  # Positive definite quantity

    # Fortran uses kp1 = min(k+1,nzt), which maps the interior Python slice to
    # the adjacent upper momentum levels at indices 1:nzt.

    # Calculate the <x'^3> turbulent production (tp) term.
    term_tp_interior = term_tp_rhs(
        xp2_zt[:, :nzt - 1], wpxp[:, 1:nzt], wpxp[:, :nzt - 1],
        rho_ds_zm[:, 1:nzt], rho_ds_zm[:, :nzt - 1],
        invrs_rho_ds_zt[:, :nzt - 1],
        gr.invrs_dzt[:, :nzt - 1],
    )

    # Calculate the <x'^3> accumulation (ac) term.
    term_ac_interior = term_ac_rhs(
        xm_zm[:, 1:nzt], xm_zm[:, :nzt - 1], wpxp2[:, :nzt - 1],
        gr.invrs_dzt[:, :nzt - 1],
    )

    term_tp = term_tp.at[:, :nzt - 1].set(term_tp_interior)
    term_ac = term_ac.at[:, :nzt - 1].set(term_ac_interior)

    if l_predict_xp3:

        xp3_interior = (
            ((xp3[:, :nzt - 1] / dt) + term_tp_interior + term_ac_interior)
            / (
                (one / dt)
                + (
                    C_xp3_dissipation
                    * invrs_tau_zt[:, :nzt - 1]
                )
            )
        )

        if l_lmm_stepping:
            xp3_interior = one_half * (xp3[:, :nzt - 1] + xp3_interior)

    else:

        # Calculate <x'^3> using the steady-state approximation.
        xp3_interior = (
            jnp.minimum(
                one / invrs_tau_zt[:, :nzt - 1],
                tau_max_zt[:, :nzt - 1],
            )
            * one / C_xp3_dissipation
            * (term_tp_interior + term_ac_interior)
        )

    xp3 = xp3.at[:, :nzt - 1].set(xp3_interior)

    # Set Upper Boundary Condition
    xp3 = xp3.at[:, nzt - 1].set(zero)

    if stats.l_sample:
        stats = stats.update(name_tp, term_tp)
        stats = stats.update(name_ac, term_ac)
        stats = stats.update(
            name_dp,
            -(C_xp3_dissipation * invrs_tau_zt) * xp3,
        )
        if l_predict_xp3:
            stats = stats.finalize_budget(name_bt, xp3 / dt)

    return xp3, stats


@jax.jit
def term_tp_rhs(
    xp2_zt, wpxpp1, wpxp,
    rho_ds_zmp1, rho_ds_zm,
    invrs_rho_ds_zt,
    invrs_dzt,
):
    """Turbulent production of <x'^3>:  explicit portion of the code.

    The d<x'^3>/dt equation contains a turbulent production term:

    + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz.

    The <x'^3> turbulent production term is completely explicit and is
    discretized as follows:

    The values of <x'^3> are found on the thermodynamic levels, while the
    values of <w'x'> and <x'^2> are found on the momentum levels.

    invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) )
    """

    # The <x'^3> turbulent production term.
    term_tp = (
        + three * xp2_zt * invrs_rho_ds_zt
        * invrs_dzt * (rho_ds_zmp1 * wpxpp1 - rho_ds_zm * wpxp)
    )

    return term_tp


@jax.jit
def term_ac_rhs(
    xm_zmp1, xm_zm, wpxp2,
    invrs_dzt,
):
    """Accumulation of <x'^3>:  explicit portion of the code.

    The d<x'^3>/dt equation contains an accumulation term:

    - 3 * <w'x'^2> * d<x>/dz.

    The <x'^3> accumulation term is completely explicit and is discretized as
    follows:

    The values of <x'^3>, <x>, and <w'x'^2> are found on thermodynamic levels.

    invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) )
    """

    # The <x'^3> accumulation term.
    term_ac = -three * wpxp2 * invrs_dzt * (xm_zmp1 - xm_zm)

    return term_ac
