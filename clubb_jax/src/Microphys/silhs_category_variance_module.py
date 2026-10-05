"""Category diagnostics from silhs_category_variance_module.F90.
"""

import jax
import jax.numpy as jnp
from clubb_jax.src.SILHS.silhs_importance_sample_module import (
    define_importance_categories,
    compute_category_real_probs,
    determine_sample_categories,
)
from clubb_jax.src.CLUBB_core.advance_helper_module import sqrt_clipped


# -----------------------------------------------------------------------------
def silhs_category_variance_driver(
    ngrdcol, nzt, num_samples, pdf_dim, hydromet_dim, hm_metadata,  # In
    X_nl_all_levs, X_mixt_comp_all_levs, lh_hydromet_mc_all,        # In
    lh_sample_point_weights, pdf_params, precip_fracs,              # In
    stats,                                                          # InOut
):
    """Compute the variance of a microphysics variable in each importance category for all model
    columns.

    Arguments:
        ngrdcol: Number of model columns
        nzt: Number of height levels
        num_samples: Number of SILHS sample points
        pdf_dim: Number of variates in X_nl
        hydromet_dim: Number of elements of hydromet array
        hm_metadata: Hydrometeor/PDF variable index metadata
        X_nl_all_levs: SILHS samples at all height levels
        X_mixt_comp_all_levs: Mixture component (1 or 2) of each sample point
        lh_hydromet_mc_all: Tendencies of hydrometeors at all sample points
        lh_sample_point_weights: Weight of SILHS sample points
        pdf_params: The PDF parameters
        precip_fracs: Precipitation fractions [-]
        stats: JAX statistics state; updated state is returned
    """
    # Sample rrm_mc from the hydrometeor tendencies generated at each SILHS point.
    samples_all = lh_hydromet_mc_all[..., hm_metadata.iirr]
    root_weight_mean_sq_cat = silhs_sample_category_variance(
        ngrdcol, nzt, num_samples, pdf_dim, X_nl_all_levs,  # In
        X_mixt_comp_all_levs, samples_all,                  # In
        lh_sample_point_weights, pdf_params, precip_fracs,  # In
        hm_metadata,                                        # In
    )
    if stats.l_sample:
        for icat in range(8):
            stats = stats.update(f"silhs_var_cat_{icat+1}", root_weight_mean_sq_cat[..., icat])
    return stats


# -----------------------------------------------------------------------------
def silhs_sample_category_variance(
    ngrdcol, nzt, num_samples, pdf_dim, X_nl_all_levs,  # In
    X_mixt_comp_all_levs, samples_all,                  # In
    lh_sample_point_weights, pdf_params, precip_fracs,  # In
    hm_metadata,                                        # In
):
    """Compute category-conditioned root-mean-square values in every column.

    Arguments:
        ngrdcol: Number of model columns
        nzt: Number of height levels
        num_samples: Number of SILHS sample points
        pdf_dim: Number of variates in X_nl
        X_nl_all_levs: SILHS samples at all height levels
        X_mixt_comp_all_levs: Mixture component (1 or 2) of each sample point
        samples_all: Sample points of variable to compute variance of
        lh_sample_point_weights: Weight of SILHS sample points
        pdf_params: The PDF parameters
        precip_fracs: Precipitation fractions [-]
        hm_metadata: Hydrometeor/PDF variable index metadata
    """
    # Keep classifications consistent with define_importance_categories.
    # The source classification/PDF-probability loops are independent and batched.
    importance_categories = define_importance_categories()
    int_sample_category = determine_sample_categories(
        num_samples, pdf_dim, hm_metadata,            # In
        X_nl_all_levs,                                # In
        X_mixt_comp_all_levs, importance_categories,  # In
    )
    category_real_probs = jax.vmap(
        jax.vmap(
            lambda cf1, cf2, m, p1, p2: compute_category_real_probs(
                importance_categories,  # In
                cf1, cf2,               # In
                m,                      # In
                p1, p2,                 # In
            )
        )
    )(
        pdf_params.cloud_frac_1,
        pdf_params.cloud_frac_2,
        pdf_params.mixt_frac,
        precip_fracs.precip_frac_1,
        precip_fracs.precip_frac_2,
    )

    # Accumulate the weighted squared tendency in each category and average
    # over all samples, then condition on the category's real PDF probability.
    root_weight_mean_sq_cat = (
        jnp.sum(
            lh_sample_point_weights[..., None]
            * samples_all[..., None] ** 2
            * (int_sample_category[..., None] == jnp.arange(8)),
            axis=1,
        )
        / num_samples
    )

    # Categories with no PDF probability retain the source missing value -999.
    root_weight_mean_sq_cat = jnp.where(
        category_real_probs > 0.0,
        sqrt_clipped(
            root_weight_mean_sq_cat / jnp.where(category_real_probs > 0.0, category_real_probs, 1.0)
        ),
        -999.0,
    )
    return root_weight_mean_sq_cat
