"""Sampled microphysics driver from lh_microphys_driver_module.F90.
"""

import jax.numpy as jnp
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from clubb_jax.src.SILHS.est_kessler_microphys_module import est_kessler_microphys_api
from clubb_jax.src.Microphys.estimate_scm_microphys_module import est_silhs_tndcy


# -----------------------------------------------------------------------------
def lh_microphys_driver(
    gr, ngrdcol, dt, nzt, nzm, num_samples,         # In
    pdf_dim, hydromet_dim, hm_metadata,             # In
    X_nl_all_levs, lh_sample_point_weights,         # In
    pdf_params, precip_fracs, p_in_Pa, exner, rho,  # In
    rcm, delta_zt, cloud_frac,                      # In
    hydromet, X_mixt_comp_all_levs,                 # In
    lh_rt_clipped, lh_thl_clipped,                  # In
    lh_rc_clipped, lh_rv_clipped,                   # In
    lh_Nc_clipped,                                  # In
    l_lh_importance_sampling,                       # In
    l_lh_instant_var_covar_src,                     # In
    saturation_formula,                             # In
    stats,                                          # InOut
    microphys_sub,                                  # In
):
    """Estimate microphysics changes given subcolumns of thlm, rtm,
    et cetera from the subcolumn generator.

    Return source outputs followed by the Kessler diagnostic failure mask.

    Arguments:
        gr: Grid metadata and coordinate arrays
        ngrdcol: Number of model columns
        dt: Model timestep [s]
        nzt: Number of thermodynamic vertical model levels
        nzm: Number of momentum vertical model levels
        num_samples: Number of calls to microphysics per timestep (normally=2)
        pdf_dim: Number of variables to sample
        hydromet_dim: Number of precipitating hydrometeor fields.
        hm_metadata: Hydrometeor/PDF variable index metadata
        X_nl_all_levs: Sample that is transformed ultimately to normal-lognormal
        lh_sample_point_weights: Weight given the individual sample points
        pdf_params: PDF parameters [units vary]
        precip_fracs: Precipitation fractions [-]
        p_in_Pa: Pressure [Pa]
        exner: Exner function [-]
        rho: Density on thermo. grid [kg/m^3]
        rcm: Liquid water mixing ratio [kg/kg]
        delta_zt: Change in meters with height [m]
        cloud_frac: Cloud fraction [-]
        hydromet: Hydrometeor species [units vary]
        X_mixt_comp_all_levs: Which mixture component we're in
        lh_rt_clipped: rt generated from silhs sample points
        lh_thl_clipped: thl generated from silhs sample points
        lh_rc_clipped: rc generated from silhs sample points
        lh_rv_clipped: rv generated from silhs sample points
        lh_Nc_clipped: Nc generated from silhs sample points
        l_lh_importance_sampling: Do importance sampling (SILHS) [-]
        l_lh_instant_var_covar_src: Produce instantaneous var/covar tendencies [-]
        saturation_formula: Integer that stores the saturation formula to be used
        stats: JAX statistics state; updated state is returned
        microphys_sub: Static microphysics procedure (local KK or Morrison)
    """
    # ---- Begin Code ----
    # TODO(port-mirror): compiled JAX cannot perform source ERROR STOPs;
    # append per-column diagnostic status to the returned state for the host.
    l_error = jnp.zeros(ngrdcol, dtype=bool)

    # Perform LH and analytic Kessler calculations as a diagnostic test of SILHS.
    if clubb_at_least_debug_level(2):
        (
            lh_AKm, AKm, AKstd, AKstd_cld, AKm_rcm, AKm_rcc, lh_rcm_avg, l_error,
        ) = est_kessler_microphys_api(
            nzt, num_samples, pdf_dim, ngrdcol,             # In
            X_nl_all_levs, pdf_params, rcm, cloud_frac,     # In
            X_mixt_comp_all_levs, lh_sample_point_weights,  # In
            l_lh_importance_sampling,                       # In
        )
    else:
        lh_AKm = AKm = AKstd = AKstd_cld = lh_rcm_avg = AKm_rcm = AKm_rcc = jnp.zeros_like(rcm)

    # Call the Latin-hypercube microphysics driver for microphys_sub.
    (
        stats,
        lh_hydromet_mc,
        lh_hydromet_vel,
        lh_Ncm_mc,
        lh_rvm_mc,
        lh_rcm_mc,
        lh_thlm_mc,
        lh_rtp2_mc,
        lh_thlp2_mc,
        lh_wprtp_mc,
        lh_wpthlp_mc,
        lh_rtpthlp_mc,
    ) = est_silhs_tndcy(
        gr, ngrdcol, dt, nzt, nzm, num_samples,                        # In
        pdf_dim, hydromet_dim, hm_metadata,                            # In
        X_nl_all_levs, X_mixt_comp_all_levs, lh_sample_point_weights,  # In
        pdf_params, precip_fracs, p_in_Pa, exner, rho,                 # In
        delta_zt, hydromet, rcm,                                       # In
        lh_rt_clipped, lh_thl_clipped,                                 # In
        lh_rc_clipped, lh_rv_clipped,                                  # In
        lh_Nc_clipped,                                                 # In
        l_lh_instant_var_covar_src,                                    # In
        saturation_formula,                                            # In
        stats,                                                         # InOut
        microphys_sub,                                                 # In
    )
    return (
        stats,
        lh_hydromet_mc,
        lh_hydromet_vel,
        lh_Ncm_mc,
        lh_rcm_mc,
        lh_rvm_mc,
        lh_thlm_mc,
        lh_rtp2_mc,
        lh_thlp2_mc,
        lh_wprtp_mc,
        lh_wpthlp_mc,
        lh_rtpthlp_mc,
        lh_AKm,
        AKm,
        AKstd,
        AKstd_cld,
        lh_rcm_avg,
        AKm_rcm,
        AKm_rcc,
        l_error,
    )
