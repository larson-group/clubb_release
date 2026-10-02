"""Disabled interfaces from lh_microphys_driver_module.F90.

Initialization rejects this feature. These source-signature dummies must never
silently simulate an enabled feature; replace them when its driver is ported.
"""


def lh_microphys_driver(gr, ngrdcol, dt, nzt, nzm,
    num_samples, pdf_dim, hydromet_dim, hm_metadata, X_nl_all_levs,
    lh_sample_point_weights, pdf_params, precip_fracs, p_in_Pa, exner,
    rho, rcm, delta_zt, cloud_frac, hydromet,
    X_mixt_comp_all_levs, lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, lh_rv_clipped,
    lh_Nc_clipped, l_lh_importance_sampling, l_lh_instant_var_covar_src, saturation_formula, stats,
    lh_hydromet_mc, lh_hydromet_vel, lh_Ncm_mc, lh_rcm_mc, lh_rvm_mc,
    lh_thlm_mc, lh_rtp2_mc, lh_thlp2_mc, lh_wprtp_mc, lh_wpthlp_mc,
    lh_rtpthlp_mc, lh_AKm, AKm, AKstd, AKstd_cld,
    lh_rcm_avg, AKm_rcm, AKm_rcc, microphys_sub):
    # Description:
    #   Computes an estimate of the change due to microphysics given a set of
    #   subcolumns of thlm, rtm, et cetera from the subcolumn generator.
    #
    # References:
    #   None
    #---------------------------------------------------------------------------
    #
    raise NotImplementedError("SILHS is disabled in the JAX standalone")
