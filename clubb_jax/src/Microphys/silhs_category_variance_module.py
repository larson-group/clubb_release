"""Disabled interfaces from silhs_category_variance_module.F90.

Initialization rejects this feature. These source-signature dummies must never
silently simulate an enabled feature; replace them when its driver is ported.
"""


def silhs_category_variance_driver(ngrdcol, nzt, num_samples, pdf_dim, hydromet_dim,
    hm_metadata, X_nl_all_levs, X_mixt_comp_all_levs, lh_hydromet_mc_all, lh_sample_point_weights,
    pdf_params, precip_fracs, stats):
    # Description:
    #   Compute the variance of a microphysics variable in each importance
    #   category for all model columns.
    #---------------------------------------------------------------------------
    #
    raise NotImplementedError("SILHS category diagnostics is disabled in the JAX standalone")


def silhs_sample_category_variance(ngrdcol, nzt, num_samples, pdf_dim, X_nl_all_levs,
    X_mixt_comp_all_levs, samples_all, lh_sample_point_weights, pdf_params, precip_fracs,
    hm_metadata, root_weight_mean_sq_cat):
    # Description:
    #   Compute category-conditioned root-mean-square values in every column.
    #---------------------------------------------------------------------------
    #
    raise NotImplementedError("SILHS category diagnostics is disabled in the JAX standalone")
