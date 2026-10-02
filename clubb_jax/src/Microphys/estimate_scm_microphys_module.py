"""Disabled interfaces from estimate_scm_microphys_module.F90.

Initialization rejects this feature. These source-signature dummies must never
silently simulate an enabled feature; replace them when its driver is ported.
"""


def est_silhs_tndcy(gr, ngrdcol, dt, nzt, nzm,
    num_samples, pdf_dim, hydromet_dim, hm_metadata, X_nl_all_levs,
    X_mixt_comp_all_levs, lh_sample_point_weights, pdf_params, precip_fracs, p_in_Pa,
    exner, rho, dzq, hydromet, rcm,
    lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, lh_rv_clipped, lh_Nc_clipped,
    l_lh_instant_var_covar_src, saturation_formula, stats, lh_hydromet_mc, lh_hydromet_vel,
    lh_Ncm_mc, lh_rvm_mc, lh_rcm_mc, lh_thlm_mc, lh_rtp2_mc,
    lh_thlp2_mc, lh_wprtp_mc, lh_wpthlp_mc, lh_rtpthlp_mc, microphys_sub):
    # Description:
    #   Estimate the tendency of a microphysics scheme via latin hypercube sampling
    #
    # References:
    #   None
    #-------------------------------------------------------------------------------
    #
    raise NotImplementedError("SILHS sampled microphysics is disabled in the JAX standalone")


def adjust_KK_src_means(dt, nzt, ngrdcol, exner, rcm,
    rrm, Nrm, hydromet, hydromet_dim, iiri,
    rrm_auto, rrm_accr, rrm_evap, Nrm_auto, Nrm_evap,
    lh_Vrr, lh_VNr, rrm_mc, Nrm_mc, rvm_mc,
    rcm_mc, thlm_mc, rrm_src_adj, Nrm_src_adj, rrm_evap_adj,
    Nrm_evap_adj):
    # Description:
    #   Adjusts the means of microphysics terms for KK microphysics by calling the
    #   KK microphysics adjustment subroutine for every model column.
    #
    # References:
    #   clubb:ticket:558
    #-----------------------------------------------------------------------------
    #
    raise NotImplementedError("SILHS sampled microphysics is disabled in the JAX standalone")
