"""Microphysical variance/covariance sources from lh_microphys_var_covar_module.F90."""

from clubb_jax.src.SILHS.math_utilities import (
    compute_sample_mean,
    compute_sample_variance,
    compute_sample_covariance,
)


# -----------------------------------------------------------------------------
def lh_microphys_var_covar_driver_api(
    nzt, num_samples, ngrdcol, dt, lh_sample_point_weights,  # In
    pdf_params, lh_rt_all, lh_thl_all, lh_w_all,             # In
    lh_rcm_mc_all, lh_rvm_mc_all, lh_thlm_mc_all,            # In
    l_lh_instant_var_covar_src,                              # In
):
    """Computes the effect of microphysics on gridbox variances and covariances
    for all model columns, including single-column runs with ngrdcol=1.
    More description:
    The equations for the (co)variance microphysical tendencies, when
    integrated forward in time explicitly, are:
    rtp2_mc    = 2*covar(rt,rt_mc) + dt*var(rt_mc)
    thlp2_mc   = 2*covar(thl,thl_mc) + dt*var(thl_mc)
    wprtp_mc   = covar(w,rt_mc)
    wpthlp_mc  = covar(w,thl_mc)
    rtpthlp_mc = covar(thl,rt_mc) + covar(rt,thl_mc) + dt*covar(rt_mc,thl_mc)
    This code can optionally take the limit of these equations at an
    infinitesimally small time step, such that the terms involving
    dt drop out. This configuration agrees with the KK upscaled analytic
    solution. (See clubb:ticket:753 for more discussion on this.)

    Arguments:
        nzt: Number of vertical levels
        num_samples: Number of SILHS sample points
        ngrdcol: Number of model columns
        dt: Model time step [s]
        lh_sample_point_weights: Weight of SILHS sample points
        pdf_params: The PDF parameters
        lh_rt_all: SILHS samples of total water [kg/kg]
        lh_thl_all: SILHS samples of potential temperature [K]
        lh_w_all: SILHS samples of vertical velocity [m/s]
        lh_rcm_mc_all: SILHS microphys. tendency of rcm [kg/kg/s]
        lh_rvm_mc_all: SILHS microphys. tendency of rvm [kg/kg/s]
        lh_thlm_mc_all: SILHS microphys. tendency of thlm [K/s]
        l_lh_instant_var_covar_src: Produce instantaneous var/covar tendencies [-]
    """
    # ---- Begin Code ----

    lh_rt_mc_all = lh_rcm_mc_all + lh_rvm_mc_all

    # Calculate means, variances, and covariances needed for the tendency terms.
    mean_rt = (
        pdf_params.mixt_frac * pdf_params.rt_1 + (1.0 - pdf_params.mixt_frac) * pdf_params.rt_2
    )
    mean_thl = (
        pdf_params.mixt_frac * pdf_params.thl_1 + (1.0 - pdf_params.mixt_frac) * pdf_params.thl_2
    )
    mean_w = pdf_params.mixt_frac * pdf_params.w_1 + (1.0 - pdf_params.mixt_frac) * pdf_params.w_2
    mean_rt_mc = compute_sample_mean(
        nzt, num_samples, ngrdcol,              # In
        lh_sample_point_weights, lh_rt_mc_all,  # In
    )
    covar_rt_rt_mc = compute_sample_covariance(
        nzt, num_samples, ngrdcol,                    # In
        lh_sample_point_weights, lh_rt_all, mean_rt,  # In
        lh_rt_mc_all, mean_rt_mc,                     # In
    )
    mean_thl_mc = compute_sample_mean(
        nzt, num_samples, ngrdcol,                # In
        lh_sample_point_weights, lh_thlm_mc_all,  # In
    )
    covar_thl_thl_mc = compute_sample_covariance(
        nzt, num_samples, ngrdcol,                      # In
        lh_sample_point_weights, lh_thl_all, mean_thl,  # In
        lh_thlm_mc_all, mean_thl_mc,                    # In
    )
    covar_w_rt_mc = compute_sample_covariance(
        nzt, num_samples, ngrdcol,                  # In
        lh_sample_point_weights, lh_w_all, mean_w,  # In
        lh_rt_mc_all, mean_rt_mc,                   # In
    )
    covar_w_thl_mc = compute_sample_covariance(
        nzt, num_samples, ngrdcol,                  # In
        lh_sample_point_weights, lh_w_all, mean_w,  # In
        lh_thlm_mc_all, mean_thl_mc,                # In
    )
    covar_thl_rt_mc = compute_sample_covariance(
        nzt, num_samples, ngrdcol,                      # In
        lh_sample_point_weights, lh_thl_all, mean_thl,  # In
        lh_rt_mc_all, mean_rt_mc,                       # In
    )
    covar_rt_thl_mc = compute_sample_covariance(
        nzt, num_samples, ngrdcol,                    # In
        lh_sample_point_weights, lh_rt_all, mean_rt,  # In
        lh_thlm_mc_all, mean_thl_mc,                  # In
    )

    # Compute the microphysical variance and covariance tendencies.
    lh_rtp2_mc_zt = 2.0 * covar_rt_rt_mc
    lh_thlp2_mc_zt = 2.0 * covar_thl_thl_mc
    lh_wprtp_mc_zt = covar_w_rt_mc
    lh_wpthlp_mc_zt = covar_w_thl_mc
    lh_rtpthlp_mc_zt = covar_thl_rt_mc + covar_rt_thl_mc
    if not l_lh_instant_var_covar_src:
        # Variances and covariances for timestep-dependent terms.
        # These terms arise when rtm and thlm are integrated forward explicitly.
        # KK upscaled omits them; including them prevents convergence with that
        # analytic solution (see clubb:ticket:753).
        var_rt_mc = compute_sample_variance(
            nzt, num_samples, ngrdcol,                          # In
            lh_rt_mc_all, lh_sample_point_weights, mean_rt_mc,  # In
        )
        var_thl_mc = compute_sample_variance(
            nzt, num_samples, ngrdcol,                             # In
            lh_thlm_mc_all, lh_sample_point_weights, mean_thl_mc,  # In
        )
        covar_rt_mc_thl_mc = compute_sample_covariance(
            nzt, num_samples, ngrdcol,                          # In
            lh_sample_point_weights, lh_rt_mc_all, mean_rt_mc,  # In
            lh_thlm_mc_all, mean_thl_mc,                        # In
        )

        # Add timestep-dependent terms.
        lh_rtp2_mc_zt += dt * var_rt_mc
        lh_thlp2_mc_zt += dt * var_thl_mc
        lh_rtpthlp_mc_zt += dt * covar_rt_mc_thl_mc
    return (
        lh_rtp2_mc_zt,
        lh_thlp2_mc_zt,
        lh_wprtp_mc_zt,
        lh_wpthlp_mc_zt,
        lh_rtpthlp_mc_zt,
    )
