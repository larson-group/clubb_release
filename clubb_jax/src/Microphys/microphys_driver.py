"""Scheme tendency dispatch from microphys_driver.F90.

JAX adaptation: columns are batched, Fortran inout/out values are returned.
"""

from functools import partial

import jax
import jax.numpy as jnp
from clubb_jax.src.Microphys import parameters_microphys as parameters
from clubb_jax.src.Microphys.KK_microphys_module import (
    KK_local_microphys_driver,
    KK_upscaled_microphys,
)
from clubb_jax.src.Microphys.cloud_sed_module import cloud_drop_sed
from clubb_jax.src.CLUBB_core.grid_class import zt2zm, zm2zt, zm2zt2zm
from clubb_jax.src.CLUBB_core.Skx_module import Skx_func
from clubb_jax.src.CLUBB_core.constants_clubb import w_tol, w_tol_sqd


# -----------------------------------------------------------------------------
@partial(
    jax.jit,
    static_argnames=(
        'ngrdcol', 'time_current', 'pdf_dim', 'hydromet_dim',
        'runtype', 'saturation_formula',
        'l_lh_importance_sampling', 'l_lh_instant_var_covar_src',
    ),
)
def calc_microphys_scheme_tendcies(
    gr, ngrdcol, dt, time_current,            # In
    pdf_dim, hydromet_dim, runtype,           # In
    thlm, p_in_Pa, exner, rho,                # In
    rho_zm, rtm,                              # In
    rcm, cloud_frac, wm_zt, wm_zm, wp2, wp3,  # In
    clubb_params,                             # In
    hydromet, Nc_in_cloud,                    # In
    hm_metadata,                              # In
    pdf_params, hydromet_pdf_params,          # In
    precip_fracs,                             # In
    X_nl_all_levs, X_mixt_comp_all_levs,      # In
    lh_sample_point_weights,                  # In
    mu_x_1_n, mu_x_2_n,                       # In
    sigma_x_1_n, sigma_x_2_n,                 # In
    corr_array_1_n, corr_array_2_n,           # In
    lh_rt_clipped, lh_thl_clipped,            # In
    lh_rc_clipped, lh_rv_clipped,             # In
    lh_Nc_clipped,                            # In
    l_lh_importance_sampling,                 # In
    l_lh_instant_var_covar_src,               # In
    saturation_formula,                       # In
    stats,                                    # InOut
    Nccnm,                                    # InOut
):
    """Call a microphysics scheme and output microphysics tendencies for the predictive variables.

    Return the source state/tendency outputs followed by a per-column mask for
    diagnostic ERROR STOPs, which the host transfers into ErrInfo.

    Arguments:
        gr: Grid coordinates, interpolation weights and vertical metrics.
        ngrdcol: Number of grid columns.
        dt: Timestep [s]
        time_current: Current time [s]
        pdf_dim: Number of variables in the multivariate PDF.
        hydromet_dim: Number of precipitating hydrometeor fields.
        runtype: Name of the run, for case specific effects.
        thlm: Liquid potential temp. [K]
        p_in_Pa: Pressure [Pa]
        exner: Exner function [-]
        rho: Density on thermodynamic levels [kg/m^3]
        rho_zm: Density on momentum levels [kg/m^3]
        rtm: Total water mixing ratio [kg/kg]
        rcm: Liquid water mixing ratio [kg/kg]
        cloud_frac: Cloud fraction [-]
        wm_zt: w wind component on thermodynamic levels [m/s]
        wm_zm: w wind component on momentum levels [m/s]
        wp2: w'^2 on the momentum grid [m^2/s^2]
        wp3: w'^3 on the thermo. grid [m^3/s^3]
        clubb_params: Column-dependent tunable CLUBB parameters.
        hydromet: Hydrometeor mean, < h_m > (thermodynamic levels) [units]
        Nc_in_cloud: Mean (in-cloud) cloud droplet concentration [num/kg]
        hm_metadata: Hydrometeor/PDF names, zero-based species indices and tolerances.
        pdf_params: PDF parameters
        hydromet_pdf_params: PDF parameters
        precip_fracs: Precipitation fractions [-]
        X_nl_all_levs: Normally and lognormally distributed hydrometeors and other variables
        X_mixt_comp_all_levs: Which mixture component the sample is in
        lh_sample_point_weights: Weights for cloud weighted sampling
        mu_x_1_n: Mean array (normal space): PDF vars. (comp. 1) [un. vary]
        mu_x_2_n: Mean array (normal space): PDF vars. (comp. 2) [un. vary]
        sigma_x_1_n: Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
        sigma_x_2_n: Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]
        corr_array_1_n: Corr. array (normal space) of PDF vars. (comp. 1) [-]
        corr_array_2_n: Corr. array (normal space) of PDF vars. (comp. 2) [-]
        lh_rt_clipped: rt generated from silhs sample points
        lh_thl_clipped: thl generated from silhs sample points
        lh_rc_clipped: rc generated from silhs sample points
        lh_rv_clipped: rv generated from silhs sample points
        lh_Nc_clipped: Nc generated from silhs sample points
        l_lh_importance_sampling: Do importance sampling (SILHS) [-]
        l_lh_instant_var_covar_src: Produce instantaneous var/covar tendencies [-]
        saturation_formula: Choice of liquid/ice saturation formula.
        stats: Immutable statistics state; return its updated value.
        Nccnm: Cloud condensation nuclei concentration (COAMPS) [num/kg]
    """

    # We must initialize intent(out) variables. If they do not otherwise get
    # assigned, they must still have defined values.
    # Description:
    # Call a microphysics scheme and output microphysics tendencies for the
    # predictive variables.
    # References:
    # H. Morrison, J. A. Curry, and V. I. Khvorostyanov, 2005: A new double-
    # moment microphysics scheme for application in cloud and
    # climate models. Part 1: Description. J. Atmos. Sci., 62, 1665-1677.
    #
    # Khairoutdinov, M. and Kogan, Y.: A new cloud physics parameterization in a
    # large-eddy simulation model of marine stratocumulus, Mon. Wea. Rev., 128,
    # 229-243, 2000.
    # -----------------------------------------------------------------------

    unit_sample_weight_2d = jnp.ones_like(rcm)
    # TODO(port-mirror): propagate source ERROR STOP status as an extra return
    # for the host's ErrInfo until JAX supports these stops inside kernels.
    l_error = jnp.zeros(ngrdcol, dtype=bool)
    hydromet_mc = jnp.zeros_like(hydromet)
    Ncm_mc = rcm_mc = rvm_mc = thlm_mc = jnp.zeros_like(rcm)
    hydromet_vel_zt = hydromet_vel_covar_zt_impc = hydromet_vel_covar_zt_expc = jnp.zeros_like(
        hydromet
    )
    wprtp_mc = wpthlp_mc = rtp2_mc = thlp2_mc = rtpthlp_mc = jnp.zeros((gr.ngrdcol, gr.nzm))

    # Calculate Skw_zm for use in advance_microphys.
    wp3_zm = zt2zm(gr.nzm, gr.nzt, ngrdcol, gr, wp3)
    Skw_zm = Skx_func(
        gr.nzm, ngrdcol, wp2, wp3_zm,  # In
        w_tol, clubb_params,           # In
    )

    # Smooth by interpolating to thermodynamic levels and back.
    Skw_zm_smooth = zm2zt2zm(gr.nzm, gr.nzt, ngrdcol, gr, Skw_zm)
    wp2_zt = zm2zt(gr.nzm, gr.nzt, ngrdcol, gr, wp2, w_tol_sqd)

    # Return if there is delay between model start and microphysics start.
    if time_current < parameters.microphys_start_time:
        return (
            stats,
            Nccnm,
            hydromet_mc,
            Ncm_mc,
            rcm_mc,
            rvm_mc,
            thlm_mc,
            hydromet_vel_zt,
            hydromet_vel_covar_zt_impc,
            hydromet_vel_covar_zt_expc,
            wprtp_mc,
            wpthlp_mc,
            rtp2_mc,
            thlp2_mc,
            rtpthlp_mc,
            Skw_zm_smooth,
            l_error,
        )
    if not runtype:
        raise ValueError("Runtype is null, which should not happen")
    # Calculate the updated mean cloud droplet concentration from Nc_in_cloud
    # and the updated cloud fraction.
    Ncm_microphys = Nc_in_cloud * cloud_frac

    # Determine 's' from Mellor (1977).
    chi = pdf_params.mixt_frac * pdf_params.chi_1 + (1.0 - pdf_params.mixt_frac) * pdf_params.chi_2

    # Compute standard deviation of vertical velocity in the grid column.
    wtmp = jnp.sqrt(wp2_zt)

    # Morrison delta_zt(k) is zt(k+1)-zt(k), which is CLUBB dzm(k+1).
    delta_zt = gr.dzm[:, 1:]

    # COAMPS and GFDL activation retain their existing initialization gates.
    # The Morrison/KK source sampling call blocks are identical; their common
    # block is kept visible here, with scheme-specific statistics below.
    if parameters.lh_microphys_type != parameters.lh_microphys_disabled:
        from clubb_jax.src.Microphys.lh_microphys_driver_module import lh_microphys_driver
        from clubb_jax.src.CLUBB_core.stats_clubb_utilities import stats_accumulate_lh_tend

        if parameters.microphys_scheme == "morrison":
            from clubb_jax.src.Microphys.morrison_microphys_module import morrison_microphys_driver

            microphys_sub = morrison_microphys_driver
        else:
            microphys_sub = KK_local_microphys_driver
        (
            stats,
            hydromet_mc,
            hydromet_vel_zt,
            Ncm_mc,
            rcm_mc,
            rvm_mc,
            thlm_mc,
            rtp2_mc,
            thlp2_mc,
            wprtp_mc,
            wpthlp_mc,
            rtpthlp_mc,
            lh_AKm,
            AKm,
            AKstd,
            AKstd_cld,
            lh_rcm_avg,
            AKm_rcm,
            AKm_rcc,
            l_error,
        ) = lh_microphys_driver(
            gr, ngrdcol, dt, gr.nzt, gr.nzm, parameters.lh_num_samples,  # In
            pdf_dim, hydromet_dim, hm_metadata,                          # In
            X_nl_all_levs, lh_sample_point_weights,                      # In
            pdf_params, precip_fracs, p_in_Pa, exner, rho,               # In
            rcm, delta_zt, cloud_frac,                                   # In
            hydromet, X_mixt_comp_all_levs,                              # In
            lh_rt_clipped, lh_thl_clipped,                               # In
            lh_rc_clipped, lh_rv_clipped,                                # In
            lh_Nc_clipped,                                               # In
            l_lh_importance_sampling,                                    # In
            l_lh_instant_var_covar_src,                                  # In
            saturation_formula,                                          # In
            stats,                                                       # InOut
            microphys_sub,                                               # In
        )
        stats = stats_accumulate_lh_tend(
            gr, ngrdcol, hydromet_dim, hm_metadata,  # In
            hydromet_mc, Ncm_mc,                     # In
            thlm_mc, rvm_mc, rcm_mc,                 # In
            lh_AKm, AKm, AKstd, AKstd_cld,           # In
            lh_rcm_avg, AKm_rcm, AKm_rcc,            # In
            stats,                                   # InOut
        )
        if parameters.microphys_scheme == "khairoutdinov_kogan":
            stats = stats.update("lh_Vrr", hydromet_vel_zt[..., hm_metadata.iirr])
            stats = stats.update("lh_VNr", hydromet_vel_zt[..., hm_metadata.iiNr])
    if (
        parameters.microphys_scheme == "morrison"
        and parameters.lh_microphys_type != parameters.lh_microphys_interactive
    ):
        from clubb_jax.src.Microphys.morrison_microphys_module import morrison_microphys_driver

        # Source l_morr_xp2_mc is gated at initialization. Its ordinary path
        # clears sampled variance tendencies before computing mean tendencies.
        rtp2_mc = thlp2_mc = wprtp_mc = wpthlp_mc = rtpthlp_mc = jnp.zeros((ngrdcol, gr.nzm))
        (
            stats,
            hydromet_mc,
            hydromet_vel_zt,
            Ncm_mc,
            rcm_mc,
            rvm_mc,
            thlm_mc,
            rrm_auto_diag,
            rrm_accr_diag,
            rrm_evap_diag,
            Nrm_auto_diag,
            Nrm_evap_diag,
        ) = morrison_microphys_driver(
            gr, ngrdcol, dt, gr.nzt,                                 # In
            hydromet_dim, hm_metadata,                               # In
            False, thlm, wm_zt, p_in_Pa,                             # In
            exner, rho, cloud_frac, wtmp,                            # In
            delta_zt, rcm, Ncm_microphys, chi, rtm - rcm, hydromet,  # In
            saturation_formula,                                      # In
            unit_sample_weight_2d,                                   # In
            stats,                                                   # InOut
        )
    elif (
        parameters.microphys_scheme == "khairoutdinov_kogan"
        and parameters.lh_microphys_type != parameters.lh_microphys_interactive
    ):
        if parameters.l_local_kk:
            # The source local-KK call replaces mean tendencies only. If
            # l_var_covar_src is enabled, its sampled variance sources remain.
            (
                stats,
                hydromet_mc,
                hydromet_vel_zt,
                Ncm_mc,
                rcm_mc,
                rvm_mc,
                thlm_mc,
                rrm_auto_diag,
                rrm_accr_diag,
                rrm_evap_diag,
                Nrm_auto_diag,
                Nrm_evap_diag,
            ) = KK_local_microphys_driver(
                gr, ngrdcol, dt, gr.nzt,                    # In
                hydromet_dim, hm_metadata,                  # In
                False,                                      # In
                thlm, wm_zt, p_in_Pa, exner, rho,           # In
                cloud_frac, wtmp, delta_zt, rcm,            # In
                Ncm_microphys, chi, rtm - rcm, hydromet,    # In
                saturation_formula, unit_sample_weight_2d,  # In
                stats,                                      # InOut
            )
        else:
            (
                stats,
                hydromet_mc,
                hydromet_vel_zt,
                rcm_mc,
                rvm_mc,
                thlm_mc,
                hydromet_vel_covar_zt_impc,
                hydromet_vel_covar_zt_expc,
                wprtp_mc,
                wpthlp_mc,
                rtp2_mc,
                thlp2_mc,
                rtpthlp_mc,
            ) = KK_upscaled_microphys(
                gr, ngrdcol, dt, gr.nzt, gr.nzm,     # In
                pdf_dim, hydromet_dim, hm_metadata,  # In
                wm_zt, rtm, thlm, p_in_Pa,           # In
                exner, rho, rcm,                     # In
                pdf_params, hydromet_pdf_params,     # In
                precip_fracs,                        # In
                hydromet,                            # In
                mu_x_1_n, mu_x_2_n,                  # In
                sigma_x_1_n, sigma_x_2_n,            # In
                corr_array_1_n, corr_array_2_n,      # In
                saturation_formula,                  # In
                stats,                               # InOut
            )
            if parameters.l_silhs_KK_convergence_adj_mean:
                hydromet_vel_covar_zt_impc = jnp.zeros_like(hydromet)
                hydromet_vel_covar_zt_expc = jnp.zeros_like(hydromet)
    elif parameters.microphys_scheme not in ("none", "morrison", "khairoutdinov_kogan"):
        raise ValueError(f"Unsupported microphysics scheme: {parameters.microphys_scheme}")
    # Source sedimentation statistics occur after either interactive or
    # ordinary microphysics, using whichever tendencies feed the model.
    if parameters.microphys_scheme in ("morrison", "khairoutdinov_kogan"):
        stats = stats.update(
            "Vrr", zt2zm(gr.nzm, gr.nzt, ngrdcol, gr, hydromet_vel_zt[..., hm_metadata.iirr])
        )
        if parameters.microphys_scheme == "khairoutdinov_kogan":
            stats = stats.update(
                "VNr", zt2zm(gr.nzm, gr.nzt, ngrdcol, gr, hydromet_vel_zt[..., hm_metadata.iiNr])
            )
    stats = stats.update("Nccnm", Nccnm)
    if parameters.l_gfdl_activation:
        raise NotImplementedError("GFDL activation core is disabled")
    # Cloud water sedimentation.
    if parameters.l_cloud_sed:
        stats, rcm_mc, thlm_mc = cloud_drop_sed(
            gr, gr.ngrdcol, rcm, Ncm_microphys,      # In
            rho_zm, rho, exner, parameters.sigma_g,  # In
            stats, rcm_mc, thlm_mc,                  # InOut
        )
    return (
        stats,
        Nccnm,
        hydromet_mc,
        Ncm_mc,
        rcm_mc,
        rvm_mc,
        thlm_mc,
        hydromet_vel_zt,
        hydromet_vel_covar_zt_impc,
        hydromet_vel_covar_zt_expc,
        wprtp_mc,
        wpthlp_mc,
        rtp2_mc,
        thlp2_mc,
        rtpthlp_mc,
        Skw_zm_smooth,
        l_error,
    )
