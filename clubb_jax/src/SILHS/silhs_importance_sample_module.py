"""Importance sampling from silhs_importance_sample_module.F90.

Category indices are zero based. Native JAX keys replace implicit generator
state; fixed category/cluster metadata permits vectorized sample operations.
"""

from typing import NamedTuple
import jax
import jax.numpy as jnp
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from clubb_jax.src.SILHS import parameters_silhs
from clubb_jax.src.SILHS.generate_uniform_sample_module import (
    rand_permute,
    choose_permuted_random,
)

num_importance_categories = 8


class importance_category_type(NamedTuple):
    """Cloud, precipitation and mixture-component membership for each category."""

    l_in_cloud: object
    l_in_precip: object
    l_in_component_1: object


# -----------------------------------------------------------------------------
def importance_sampling_driver(
    num_samples,                                                # In
    cloud_frac_1, cloud_frac_2,                                 # In
    mixt_frac,                                                  # In
    precip_frac_1, precip_frac_2,                               # In
    cluster_allocation_strategy, l_lh_clustered_sampling,       # In
    l_lh_limit_weights, l_lh_var_frac, l_lh_normalize_weights,  # In
    X_u_chi_one_lev, X_u_dp1_one_lev, X_u_dp2_one_lev,          # InOut
    key,                                                        # In
):
    """Apply importance sampling to a single vertical level.

    Source references: CLUBB ticket 736;
    https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:importance_sampling.

    Arguments:
        num_samples: Number of SILHS sample points
        cloud_frac_1: Cloud fraction in PDF component 1 [-]
        cloud_frac_2: Cloud fraction in PDF component 2 [-]
        mixt_frac: Weight of PDF component 1 [-]
        precip_frac_1: Precipitation fraction in PDF component 1 [-]
        precip_frac_2: Precipitation fraction in PDF component 2 [-]
        cluster_allocation_strategy: Strategy for distributing sample points
        l_lh_clustered_sampling: Use prescribed probability sampling with clusters (SILHS)
        l_lh_limit_weights: Ensure weights stay under a given value
        l_lh_var_frac: Prescribe variance fractions
        l_lh_normalize_weights: Normalize weights to sum to num_samples
        X_u_chi_one_lev: Uniform chi variates scaled by importance sampling
        X_u_dp1_one_lev: Uniform component selectors scaled by importance sampling
        X_u_dp2_one_lev: Uniform precipitation selectors scaled by importance sampling
        key: Native JAX random key for these draws
    """
    importance_categories = define_importance_categories()
    l_error = jnp.array(False)
    category_real_probs = compute_category_real_probs(
        importance_categories,         # In
        cloud_frac_1, cloud_frac_2,    # In
        mixt_frac,                     # In
        precip_frac_1, precip_frac_2,  # In
    )

    # Select the cluster allocation strategy from the configured integer option.
    if l_lh_clustered_sampling:
        if cluster_allocation_strategy == parameters_silhs.eight_cluster_allocation_opt:
            category_prescribed_probs = eight_cluster_allocation(
                importance_categories, category_real_probs, l_lh_var_frac
            )
        elif cluster_allocation_strategy == parameters_silhs.four_cluster_allocation_opt:
            category_prescribed_probs = four_cluster_no_precip(
                importance_categories, category_real_probs, l_lh_var_frac
            )
        elif cluster_allocation_strategy == parameters_silhs.two_cluster_cp_nocp_opt:
            category_prescribed_probs = two_cluster_cp_nocp(
                importance_categories, category_real_probs, l_lh_var_frac
            )
        else:
            raise ValueError("Unsupported allocation strategy in importance_sampling_driver")
        if l_lh_limit_weights:
            category_prescribed_probs, l_error = limit_category_weights(
                category_real_probs,       # In
                category_prescribed_probs,  # InOut
            )
    else:
        # Place approximately half the sample points in cloud.
        category_prescribed_probs = cloud_importance_sampling(
            importance_categories, category_real_probs,  # In
            cloud_frac_1, cloud_frac_2,                  # In
            mixt_frac,                                   # In
        )

    # Compute the weight of each category, then draw stratified category selectors.
    category_sample_weights = compute_category_sample_weights(
        category_real_probs, category_prescribed_probs
    )
    rand_vect = generate_strat_uniform_variate(num_samples, key)
    int_sample_category = pick_sample_categories(num_samples, category_prescribed_probs, rand_vect)

    # A traced kernel cannot ERROR STOP. Preserve invalid-selection status while
    # using the source's first-category fallback to keep array indexing valid.
    invalid_category = jnp.any(int_sample_category < 0)
    int_sample_category = jnp.maximum(int_sample_category, 0)

    # Scale and translate each point to reside in its selected category.
    # JAX applies the source sample loop to the whole sample vector.
    category = importance_category_type(*(x[int_sample_category] for x in importance_categories))
    X_u_chi_one_lev, X_u_dp1_one_lev, X_u_dp2_one_lev = scale_sample_to_category(
        category, cloud_frac_1, cloud_frac_2,               # In
        mixt_frac,                                          # In
        precip_frac_1, precip_frac_2,                       # In
        X_u_chi_one_lev, X_u_dp1_one_lev, X_u_dp2_one_lev,  # InOut
    )

    # Pick the weight for each sample point, and normalize if enabled.
    lh_sample_point_weights = category_sample_weights[int_sample_category]
    if l_lh_normalize_weights:
        weights_sum = jnp.sum(lh_sample_point_weights)
        lh_sample_point_weights = lh_sample_point_weights * num_samples / weights_sum
    # Assertion status is returned rather than raising from traced array data.
    l_error |= invalid_category
    if clubb_at_least_debug_level(2):
        l_error |= importance_sampling_assertions(
            num_samples, importance_categories, category_real_probs,                         # In
            category_prescribed_probs, category_sample_weights, X_u_chi_one_lev,             # In
            X_u_dp1_one_lev, X_u_dp2_one_lev, lh_sample_point_weights, int_sample_category,  # In
            cloud_frac_1, cloud_frac_2,                                                      # In
            mixt_frac,                                                                       # In
            precip_frac_1, precip_frac_2,                                                    # In
            l_lh_normalize_weights,                                                          # In
        )
    return (
        X_u_chi_one_lev,
        X_u_dp1_one_lev,
        X_u_dp2_one_lev,
        lh_sample_point_weights,
        l_error,
    )


# -----------------------------------------------------------------------------
def define_importance_categories():
    """Creates a vector of size num_importance_categories that defines the eight
    importance sampling categories.
    """
    # Source order: cloudy/clear precipitating points, then cloudy/clear
    # nonprecipitating points. Each pair contains components 1 and 2.
    return importance_category_type(
        jnp.array([1, 1, 0, 0, 1, 1, 0, 0], bool),
        jnp.array([1, 1, 1, 1, 0, 0, 0, 0], bool),
        jnp.array([1, 0, 1, 0, 1, 0, 1, 0], bool),
    )


# -----------------------------------------------------------------------------
def compute_category_real_probs(
    importance_categories,         # In
    cloud_frac_1, cloud_frac_2,    # In
    mixt_frac,                     # In
    precip_frac_1, precip_frac_2,  # In
):
    """Computes the real PDF probability associated with each importance sampling
    category.
    For example, if a category is in cloud, out of precipitation, and in mixture
    component two, then without importance sampling, the probability that a
    point will appear in that category is:
    P(cloud,noprecip,comp2) = cloud_frac_2 * (1-precip_frac_2) * (1-mixt_frac)

    Arguments:
        importance_categories: A list of importance categories
        cloud_frac_1: Cloud fraction in PDF component 1 [-]
        cloud_frac_2: Cloud fraction in PDF component 2 [-]
        mixt_frac: Weight of PDF component 1 [-]
        precip_frac_1: Precipitation fraction in PDF component 1 [-]
        precip_frac_2: Precipitation fraction in PDF component 2 [-]
    """
    # Determine the component and its cloud and precipitation fractions.
    cloud_frac_i = jnp.where(importance_categories.l_in_component_1, cloud_frac_1, cloud_frac_2)
    precip_frac_i = jnp.where(importance_categories.l_in_component_1, precip_frac_1, precip_frac_2)
    component_factor = jnp.where(importance_categories.l_in_component_1, mixt_frac, 1.0 - mixt_frac)

    # Determine the cloud and precipitation factors for this category.
    cloud_factor = jnp.where(importance_categories.l_in_cloud, cloud_frac_i, 1.0 - cloud_frac_i)
    precip_factor = jnp.where(importance_categories.l_in_precip, precip_frac_i, 1.0 - precip_frac_i)

    # Compute the real PDF probability of the category.
    return component_factor * cloud_factor * precip_factor


# -----------------------------------------------------------------------------
def compute_category_sample_weights(category_real_probs, category_prescribed_probs):
    """Compute the sample point weights for a sample point in each category based
    on the PDF probability and the modified probability from importance
    sampling

    Arguments:
        category_real_probs: The actual PDF probability of each category
        category_prescribed_probs: The modified probability of each category due to importance
            sampling
    """
    # A category with no probability of being sampled has an irrelevant weight.
    absent = jnp.abs(category_prescribed_probs) < jnp.finfo(jnp.float64).eps
    return jnp.where(
        absent,
        1.0,
        category_real_probs / jnp.where(absent, 1.0, category_prescribed_probs),
    )


# -----------------------------------------------------------------------------
def limit_category_weights(
    category_real_probs,       # In
    category_prescribed_probs,  # InOut
):
    """Modifies the category prescribed probabilities such that an
    invariant maximum weight is achieved. The weight for each category
    j is equal to category_real_probs(j) / category_prescribed_probs(j).
    In this subroutine, category_prescribed_probs(j) is increased in
    each necessary category such that no category has a weight that
    exceeds the maximum.

    Return adjusted probabilities and source fatal status for the sampler.

    Arguments:
        category_real_probs: The PDF probability for each importance category
        category_prescribed_probs: Prescribed probability for each category; these will be
            modified to respect the maximum weight.
    """
    # Source constants: bound weights in categories with appreciable PDF mass.
    max_weight = 2.0
    real_prob_thresh_for_transfer = 1.0e-8
    min_presc_probs = jnp.where(
        category_real_probs >= real_prob_thresh_for_transfer,
        category_real_probs / max_weight,
        category_real_probs,
    )
    min_presc_prob_diff = category_prescribed_probs - min_presc_probs
    total_diff_under = jnp.sum(jnp.maximum(-min_presc_prob_diff, 0.0))

    # Sum the positive differences, then transfer mass to deficient categories.
    total_diff_over = jnp.sum(jnp.maximum(min_presc_prob_diff, 0.0))
    # The prescribed probabilities cannot be adjusted to achieve the maximum
    # weight. Return the source ERROR STOP status through the sampler's ErrInfo.
    l_error = total_diff_under > total_diff_over
    category_prescribed_probs = jnp.where(
        min_presc_prob_diff < 0.0,
        min_presc_probs,
        category_prescribed_probs
        - jnp.maximum(min_presc_prob_diff, 0.0)
        / jnp.where(total_diff_over > 0.0, total_diff_over, 1.0)
        * total_diff_under,
    )

    # As an assertion check, make sure this subroutine actually performed its
    # task successfully. Preserve the source's debug-level-2 guard.
    if clubb_at_least_debug_level(2):
        weight = category_real_probs / jnp.where(
            category_prescribed_probs > 0.0, category_prescribed_probs, 1.0
        )
        l_error |= jnp.any(
            (category_prescribed_probs > 0.0)
            & (category_real_probs >= real_prob_thresh_for_transfer)
            & (weight > max_weight)
        )
    return category_prescribed_probs, l_error


# -----------------------------------------------------------------------------
def pick_sample_categories(num_samples, category_prescribed_probs, rand_vect):
    """Picks a category for each sample point, based on the given probabilities,
    such that the distribution of categories of the sample points
    approximates the probabilities that are given.

    Arguments:
        num_samples: Number of sample points to be picked
        category_prescribed_probs: Prescribed probability for each category
        rand_vect: A sample of num_samples values from the uniform distribution in the range
            (0,1). This will be used to pick the categories.
    """
    # Each entry is the sum of the probabilities of all preceding categories.
    category_cumulative_probs = jnp.concatenate(
        (jnp.zeros(1), jnp.cumsum(category_prescribed_probs[:-1]))
    )

    # Pick categories from those intervals; include rand_vect == 1 in the last
    # category, as in the source fix for CLUBB ticket 805.
    selected = jnp.searchsorted(category_cumulative_probs, rand_vect, side="right") - 1
    valid = (selected >= 0) & (selected < 8) & (rand_vect >= 0.0) & (rand_vect <= 1.0)

    # JAX adaptation: return a sentinel for the caller's error-status path;
    # debug < 0 keeps the source first-category fallback.
    if clubb_at_least_debug_level(0):

        def report_invalid(_):
            jax.debug.print("Invalid rand_vect number in pick_sample_categories: {r}", r=rand_vect)

        jax.lax.cond(jnp.any(~valid), report_invalid, lambda _: None, operand=None)
        return jnp.where(valid, selected, -1)
    return jnp.where(valid, selected, 0)


# -----------------------------------------------------------------------------
def scale_sample_to_category(
    category, cloud_frac_1, cloud_frac_2,  # In
    mixt_frac,                             # In
    precip_frac_1, precip_frac_2,          # In
    X_u_chi, X_u_dp1, X_u_dp2,             # InOut
):
    """Scale and transpose a sample point to reside in the specified category

    Arguments:
        category: Scale the sample point to reside in this category
        cloud_frac_1: Cloud fraction in PDF component 1 [-]
        cloud_frac_2: Cloud fraction in PDF component 2 [-]
        mixt_frac: Weight of PDF component 1 [-]
        precip_frac_1: Precipitation fraction in PDF component 1 [-]
        precip_frac_2: Precipitation fraction in PDF component 2 [-]
        X_u_chi: Uniform sample of extended cloud-water mixing ratio
        X_u_dp1: Uniform sample of the d+1 variate
        X_u_dp2: Uniform sample of the d+2 variate
    """
    # Scale dp1 into (0, mixt_frac) for component 1, or (mixt_frac, 1)
    # for component 2. Select the corresponding cloud/precipitation fractions.
    X_u_dp1 = jnp.where(
        category.l_in_component_1,
        X_u_dp1 * mixt_frac,
        X_u_dp1 * (1.0 - mixt_frac) + mixt_frac,
    )
    cloud_frac_i = jnp.where(category.l_in_component_1, cloud_frac_1, cloud_frac_2)
    precip_frac_i = jnp.where(category.l_in_component_1, precip_frac_1, precip_frac_2)

    # Scale dp2 into (0, precip_frac_i) in precipitation, or
    # (precip_frac_i, 1) outside precipitation.
    X_u_dp2 = jnp.where(
        category.l_in_precip,
        X_u_dp2 * precip_frac_i,
        X_u_dp2 * (1.0 - precip_frac_i) + precip_frac_i,
    )

    # Scale chi into (1 - cloud_frac_i, 1) in cloud, or
    # (0, 1 - cloud_frac_i) in clear air.
    X_u_chi = jnp.where(
        category.l_in_cloud,
        X_u_chi * cloud_frac_i + (1.0 - cloud_frac_i),
        X_u_chi * (1.0 - cloud_frac_i),
    )
    return X_u_chi, X_u_dp1, X_u_dp2


# -----------------------------------------------------------------------------
def two_cluster_cp_nocp(importance_categories, category_real_probs, l_lh_var_frac):
    """Clusters importance categories into two clusters: categories that contain
    either cloud or precipitation (or both), and clusters that contain
    neither.

    Source reference: CLUBB ticket 740.

    Arguments:
        importance_categories: A list of importance categories
        category_real_probs: The real probability for each category
        l_lh_var_frac: Prescribe variance fractions
    """
    # Six cloud-or-precipitation categories; two clear, nonprecipitating ones.
    # Fixed membership replaces the source loop; -1 pads inactive cluster slots.
    return compute_clust_category_probs(
        category_real_probs, 2, 6,                                                   # In
        jnp.array([6, 2]), jnp.array([[0, 1, 2, 3, 4, 5], [6, 7, -1, -1, -1, -1]]),  # In
        jnp.array([1.0, 0.0]),                                                       # In
        l_lh_var_frac,                                                               # In
    )


# -----------------------------------------------------------------------------
def eight_cluster_allocation(importance_categories, category_real_probs, l_lh_var_frac):
    """Clusters importance categories such that each of the eight importance
    categories has its own cluster. Effectively, there are no clusters.

    Source reference: CLUBB ticket 752.

    Arguments:
        importance_categories: A list of importance categories
        category_real_probs: The real probability for each category
        l_lh_var_frac: Prescribe variance fractions
    """
    from dataclasses import fields

    # Dataclass field order follows the source eight-category probability mapping.
    probs = jnp.array(
        tuple(
            getattr(parameters_silhs.eight_cluster_presc_probs, f.name)
            for f in fields(parameters_silhs.eight_cluster_presc_probs_type)
        )
    )
    return compute_clust_category_probs(
        category_real_probs, 8, 1,                              # In
        jnp.ones(8, jnp.int32), jnp.arange(8)[:, None], probs,  # In
        l_lh_var_frac,                                          # In
    )


# -----------------------------------------------------------------------------
def four_cluster_no_precip(importance_categories, category_real_probs, l_lh_var_frac):
    """Clusters categories into four clusters for the four combinations of
    cloud/no cloud and comp 1/comp 2. Precip fraction is effectively
    ignored.

    Source reference: CLUBB ticket 752.

    Arguments:
        importance_categories: A list of importance categories
        category_real_probs: The real probability for each category
        l_lh_var_frac: Prescribe variance fractions
    """
    # Pair precipitating and nonprecipitating categories within each cloud/
    # component cluster, preserving the prescribed 0.3, 0.3, 0.2, 0.2 allocation.
    return compute_clust_category_probs(
        category_real_probs, 4, 2,                                    # In
        jnp.full(4, 2), jnp.array([[0, 4], [1, 5], [2, 6], [3, 7]]),  # In
        jnp.array([0.3, 0.3, 0.2, 0.2]),                              # In
        l_lh_var_frac,                                                # In
    )


# -----------------------------------------------------------------------------
def compute_clust_category_probs(
    category_real_probs, num_clusters, max_num_categories_in_cluster,  # In
    num_categories_in_cluster, cluster_categories, cluster_fractions,  # In
    l_lh_var_frac,                                                     # In
):
    """Calls clust_cat_probs_frm_presc_prb or clust_cat_probs_frm_var_fracs!

    Arguments:
        category_real_probs: The real probability for each category
        num_clusters: The number of clusters to sample from
        max_num_categories_in_cluster: The max number of categories in each cluster
        num_categories_in_cluster: The number of categories in each cluster
        cluster_categories: An integer matrix containing indices corresponding to the members
            of the clusters
        cluster_fractions: Prescribed fraction of some sort for each cluster
        l_lh_var_frac: Prescribe variance fractions
    """
    if l_lh_var_frac:
        return clust_cat_probs_frm_var_fracs(
            category_real_probs, num_clusters, max_num_categories_in_cluster,  # In
            num_categories_in_cluster, cluster_categories, cluster_fractions,  # In
        )
    return clust_cat_probs_frm_presc_prb(
        category_real_probs, num_clusters, max_num_categories_in_cluster,  # In
        num_categories_in_cluster, cluster_categories, cluster_fractions,  # In
    )


# -----------------------------------------------------------------------------
def clust_cat_probs_frm_var_fracs(
    category_real_probs, num_clusters, max_num_categories_in_cluster,           # In
    num_categories_in_cluster, cluster_categories, cluster_variance_fractions,  # In
):
    """This is a generalized algorithm that takes as input a set of "clusters"
    of the importance categories and a variance fraction for each
    cluster, and computes the prescribed probabilities for each category.

    Arguments:
        category_real_probs: The real probability for each category
        num_clusters: The number of clusters to sample from
        max_num_categories_in_cluster: The max number of categories in each cluster
        num_categories_in_cluster: The number of categories in each cluster
        cluster_categories: An integer matrix containing indices corresponding to the members
            of the clusters
        cluster_variance_fractions: Prescribed variance fraction for each cluster
    """
    # Compute the total PDF probability for each cluster. Ignore padded slots.
    valid = jnp.arange(max_num_categories_in_cluster)[None, :] < num_categories_in_cluster[:, None]
    cluster_real_probs = jnp.sum(
        jnp.where(valid, category_real_probs[cluster_categories], 0.0), axis=1
    )

    # Compute the sum of p_j * f_j across clusters.
    pdf_prob_var_frac_prod_sum = jnp.sum(cluster_real_probs * cluster_variance_fractions)
    eps = jnp.finfo(jnp.float64).eps

    # Prescribe probabilities proportional to PDF probability times variance
    # fraction. If no variance has nonzero PDF mass, fall back to ordinary sampling.
    cluster_prescribed_probs = jnp.where(
        jnp.abs(pdf_prob_var_frac_prod_sum) < eps,
        cluster_real_probs,
        cluster_real_probs
        * cluster_variance_fractions
        / jnp.where(jnp.abs(pdf_prob_var_frac_prod_sum) < eps, 1.0, pdf_prob_var_frac_prod_sum),
    )

    # Split each cluster probability among its categories in PDF proportions.
    category_prescribed_probs = jnp.zeros(8)
    for icluster in range(num_clusters):
        mask = jnp.any(
            (jnp.arange(8)[:, None] == cluster_categories[icluster]) & valid[icluster],
            axis=1,
        )
        category_prescribed_probs += jnp.where(
            mask & (jnp.abs(cluster_real_probs[icluster]) >= eps),
            category_real_probs
            / jnp.where(
                jnp.abs(cluster_real_probs[icluster]) < eps,
                1.0,
                cluster_real_probs[icluster],
            )
            * cluster_prescribed_probs[icluster],
            0.0,
        )
    return category_prescribed_probs


# -----------------------------------------------------------------------------
def clust_cat_probs_frm_presc_prb(
    category_real_probs, num_clusters, max_num_categories_in_cluster,         # In
    num_categories_in_cluster, cluster_categories, cluster_prescribed_probs,  # In
):
    """This is a generalized algorithm that takes as input a set of "clusters"
    of the importance categories and a prescribed probability for each
    cluster, and computes the prescribed probabilities for each category
    such that the sum of the prescribed probabilities of every category
    within a cluster is equal to the prescribed probability of the cluster.

    Arguments:
        category_real_probs: The real probability for each category
        num_clusters: The number of clusters to sample from
        max_num_categories_in_cluster: The max number of categories in each cluster
        num_categories_in_cluster: The number of categories in each cluster
        cluster_categories: An integer matrix containing indices corresponding to the members
            of the clusters
        cluster_prescribed_probs: Prescribed probability sum for each cluster
    """
    # Compute total PDF probability for each cluster, excluding padded slots.
    valid = jnp.arange(max_num_categories_in_cluster)[None, :] < num_categories_in_cluster[:, None]
    cluster_real_probs = jnp.sum(
        jnp.where(valid, category_real_probs[cluster_categories], 0.0), axis=1
    )

    # Do not importance sample clusters with extremely small PDF probability.
    l_cluster_presc_prob_modified = cluster_real_probs < parameters_silhs.importance_prob_thresh
    cluster_prescribed_probs_mod = jnp.where(
        l_cluster_presc_prob_modified, cluster_real_probs, cluster_prescribed_probs
    )

    # Distribute extra prescribed probability from thresholding to other clusters.
    # The probability difference may be negative.
    nonzero_real_clust_sum = jnp.sum(
        jnp.where(l_cluster_presc_prob_modified, 0.0, cluster_real_probs)
    )
    presc_prob_difference = jnp.sum(
        jnp.where(
            l_cluster_presc_prob_modified,
            cluster_prescribed_probs - cluster_prescribed_probs_mod,
            0.0,
        )
    )
    cluster_prescribed_probs_mod += jnp.where(
        l_cluster_presc_prob_modified,
        0.0,
        presc_prob_difference
        * cluster_real_probs
        / jnp.where(nonzero_real_clust_sum > 0.0, nonzero_real_clust_sum, 1.0),
    )

    # Split prescribed cluster probabilities into category probabilities.
    category_prescribed_probs = jnp.zeros(8)
    for icluster in range(num_clusters):
        mask = jnp.any(
            (jnp.arange(8)[:, None] == cluster_categories[icluster]) & valid[icluster],
            axis=1,
        )
        category_prescribed_probs += jnp.where(
            mask,
            jnp.where(
                l_cluster_presc_prob_modified[icluster],
                category_real_probs,
                category_real_probs
                / jnp.where(
                    cluster_real_probs[icluster] > 0.0,
                    cluster_real_probs[icluster],
                    1.0,
                )
                * cluster_prescribed_probs_mod[icluster],
            ),
            0.0,
        )
    return category_prescribed_probs


# -----------------------------------------------------------------------------
def cloud_importance_sampling(
    importance_categories, category_real_probs,  # In
    cloud_frac_1, cloud_frac_2,                  # In
    mixt_frac,                                   # In
):
    """Applies cloud weighted sampling such that approximately half of all
    sample points land in cloud and half land out of cloud!

    Arguments:
        importance_categories: A list of importance categories
        category_real_probs: The actual PDF probability for each category
        cloud_frac_1: Cloud fraction in PDF component 1 [-]
        cloud_frac_2: Cloud fraction in PDF component 2 [-]
        mixt_frac: Weight of PDF component 1 [-]
    """
    cloud_frac = mixt_frac * cloud_frac_1 + (1.0 - mixt_frac) * cloud_frac_2
    scale = jnp.where(importance_categories.l_in_cloud, 2.0 * cloud_frac, 2.0 * (1.0 - cloud_frac))
    return jnp.where(
        (cloud_frac >= 0.001) & (cloud_frac < 0.5),
        category_real_probs / jnp.where(scale > 0.0, scale, 1.0),
        category_real_probs,
    )


# -----------------------------------------------------------------------------
def importance_sampling_assertions(
    num_samples, importance_categories, category_real_probs,                         # In
    category_prescribed_probs, category_sample_weights, X_u_chi_one_lev,             # In
    X_u_dp1_one_lev, X_u_dp2_one_lev, lh_sample_point_weights, int_sample_category,  # In
    cloud_frac_1, cloud_frac_2,                                                      # In
    mixt_frac,                                                                       # In
    precip_frac_1, precip_frac_2,                                                    # In
    l_lh_normalize_weights,                                                          # In
):
    """Various assertion checks for importance sampling are performed here.

    Arguments:
        num_samples: Number of SILHS sample points
        importance_categories: The defined importance categories
        category_real_probs: The real PDF probabilities for each category
        category_prescribed_probs: Prescribed probability for each category
        category_sample_weights: Sample weight for each category
        X_u_chi_one_lev: Samples of chi in uniform space
        X_u_dp1_one_lev: Samples of the dp1 variate
        X_u_dp2_one_lev: Samples of the dp2 variate
        lh_sample_point_weights: Weights of samples
        int_sample_category: An integer for each sample corresponding to the category picked
            for the sample
        cloud_frac_1: Cloud fraction in PDF component 1 [-]
        cloud_frac_2: Cloud fraction in PDF component 2 [-]
        mixt_frac: Weight of PDF component 1 [-]
        precip_frac_1: Precipitation fraction in PDF component 1 [-]
        precip_frac_2: Precipitation fraction in PDF component 2 [-]
        l_lh_normalize_weights: Normalize weights to sum to num_samples
    """
    # Assert that real, prescribed and weighted prescribed probabilities sum to 1.
    tolerance = num_importance_categories * jnp.finfo(jnp.float64).eps
    l_error = (
        (jnp.abs(jnp.sum(category_real_probs) - 1.0) > tolerance)
        | (jnp.abs(jnp.sum(category_prescribed_probs) - 1.0) > tolerance)
        | (jnp.abs(jnp.sum(category_sample_weights * category_prescribed_probs) - 1.0) > tolerance)
    )

    # If enabled, assert that the sample weights average to 1.
    if l_lh_normalize_weights:
        l_error |= (
            jnp.abs(jnp.mean(lh_sample_point_weights) - 1.0)
            > num_samples * jnp.finfo(jnp.float64).eps
        )

    # Verify mixture component, cloud and precipitation intervals for each sample.
    category = importance_category_type(*(x[int_sample_category] for x in importance_categories))
    cloud_frac_i = jnp.where(category.l_in_component_1, cloud_frac_1, cloud_frac_2)
    precip_frac_i = jnp.where(category.l_in_component_1, precip_frac_1, precip_frac_2)
    l_error |= jnp.any(
        jnp.where(
            category.l_in_component_1,
            X_u_dp1_one_lev > mixt_frac,
            X_u_dp1_one_lev < mixt_frac,
        )
    )
    l_error |= jnp.any(
        jnp.where(
            category.l_in_cloud,
            X_u_chi_one_lev < 1.0 - cloud_frac_i,
            X_u_chi_one_lev > 1.0 - cloud_frac_i,
        )
    )
    l_error |= jnp.any(
        jnp.where(
            category.l_in_precip,
            X_u_dp2_one_lev > precip_frac_i,
            X_u_dp2_one_lev < precip_frac_i,
        )
    )
    return l_error


# -----------------------------------------------------------------------------
def cloud_weighted_sampling_driver(
    num_samples, p_matrix_chi, p_matrix_dp1,  # In
    cloud_frac_1, cloud_frac_2, mixt_frac,    # In
    X_u_chi, X_u_dp1,                         # InOut
    key,                                      # In
):
    """Performs importance sampling such that half of sample points are in cloud

    Arguments:
        num_samples: Number of SILHS sample points
        p_matrix_chi: Permutation of 0..num_samples-1 for chi
        p_matrix_dp1: Elements from p_matrix for dp1 element
        cloud_frac_1: Cloud fraction in PDF component 1
        cloud_frac_2: Cloud fraction in PDF component 2
        mixt_frac: Weight of first gaussian component
        X_u_chi: Samples of chi in uniform space
        X_u_dp1: Samples of the dp1 variate for determining mixture component
        key: Native JAX random key for these draws
    """
    cloud_frac = mixt_frac * cloud_frac_1 + (1.0 - mixt_frac) * cloud_frac_2
    half = num_samples // 2

    # The lower half of dp1 permutation entries selects cloudy mixture draws;
    # the upper half selects clear draws. Each half remains stratified.
    cloud_indices = jnp.nonzero(p_matrix_dp1 < half, size=half)[0]
    clear_indices = jnp.nonzero(p_matrix_dp1 >= half, size=half)[0]
    cloud_key, clear_key, chi_key = jax.random.split(key, 3)
    mixt_rand_cloud = choose_permuted_random(half, p_matrix_dp1[cloud_indices], cloud_key)
    mixt_rand_clear = choose_permuted_random(half, p_matrix_dp1[clear_indices] - half, clear_key)

    # The chi permutation assigns half the points to cloudy air and half to clear.
    l_cloudy_sample = p_matrix_chi >= half
    cloud_rank = jnp.cumsum(l_cloudy_sample) - 1
    clear_rank = jnp.cumsum(~l_cloudy_sample) - 1
    mixt_rand_element = jnp.where(
        l_cloudy_sample,
        mixt_rand_cloud[jnp.maximum(cloud_rank, 0)],
        mixt_rand_clear[jnp.maximum(clear_rank, 0)],
    )
    X_u_dp1_scaled, X_u_chi_scaled = choose_X_u_scaled(
        l_cloudy_sample,               # In
        p_matrix_chi, num_samples,     # In
        cloud_frac_1, cloud_frac_2,    # In
        mixt_frac, mixt_rand_element,  # In
        chi_key,                       # In
    )

    # Apply cloud/clear weights only in the source cloud-fraction sampling range.
    enabled = (cloud_frac >= 0.001) & (cloud_frac < 0.5)
    lh_sample_point_weights = jnp.where(
        enabled,
        jnp.where(l_cloudy_sample, 2.0 * cloud_frac, 2.0 - 2.0 * cloud_frac),
        1.0,
    )
    return (
        jnp.where(enabled, X_u_chi_scaled, X_u_chi),
        jnp.where(enabled, X_u_dp1_scaled, X_u_dp1),
        lh_sample_point_weights,
    )


# -----------------------------------------------------------------------------
def generate_strat_uniform_variate(num_samples, key):
    """Generates a stratified uniform sample for a single variable

    Arguments:
        num_samples: Number of SILHS sample points
        key: Native JAX random key for these draws
    """
    # Randomly permute the strata, then draw one real variate within each box.
    perm_key, draw_key = jax.random.split(key)
    pvect = rand_permute(num_samples, perm_key)
    return choose_permuted_random(num_samples, pvect, draw_key)


# -----------------------------------------------------------------------------
def choose_X_u_scaled(
    l_cloudy_sample,                # In
    p_matrix_element, num_samples,  # In
    cloud_frac_1, cloud_frac_2,     # In
    mixt_frac, mixt_rand_element,   # In
    key,                            # In
):
    """Find a clear or cloudy point for sampling.

    Arguments:
        l_cloudy_sample: Whether this is a cloudy or clear-air sample point
        p_matrix_element: Integer from 0..num_samples-1 for this sample
        num_samples: Total number of calls to the microphysics
        cloud_frac_1: Cloud fraction associated with mixture component 1 [-]
        cloud_frac_2: Cloud fraction associated with mixture component 2 [-]
        mixt_frac: Mixture fraction [-]
        mixt_rand_element: Random number (0,1) for determining mixture component
        key: Native JAX random key for these draws
    """
    # Determine the conditional mixture fraction given cloud/clear membership.
    cld_comp1_frac = mixt_frac * cloud_frac_1
    cld_comp2_frac = (1.0 - mixt_frac) * cloud_frac_2
    nocld_comp1_frac = mixt_frac - cld_comp1_frac
    nocld_comp2_frac = (1.0 - mixt_frac) - cld_comp2_frac
    denom = jnp.where(
        l_cloudy_sample,
        cld_comp1_frac + cld_comp2_frac,
        nocld_comp1_frac + nocld_comp2_frac,
    )
    conditional_mixt_frac = jnp.where(
        l_cloudy_sample, cld_comp1_frac, nocld_comp1_frac
    ) / jnp.where(denom > 0.0, denom, 1.0)

    # Determine the mixture component from the conditional mixture fraction.
    first = mixt_rand_element < conditional_mixt_frac

    # Rescale the selector to (0, 1), then map dp1 into the chosen component.
    mixt_rand_element_scaled = jnp.where(
        first,
        mixt_rand_element / jnp.where(conditional_mixt_frac > 0.0, conditional_mixt_frac, 1.0),
        (mixt_rand_element - conditional_mixt_frac)
        / jnp.where(conditional_mixt_frac < 1.0, 1.0 - conditional_mixt_frac, 1.0),
    )
    X_u_dp1_element = jnp.where(
        first,
        mixt_rand_element_scaled * mixt_frac,
        mixt_rand_element_scaled * (1.0 - mixt_frac) + mixt_frac,
    )

    # Range the permutation element to its half, draw a stratified chi variate,
    # and scale/translate it into cloudy or clear air in the chosen component.
    chi_rand_element = choose_permuted_random(
        num_samples // 2, p_matrix_element % (num_samples // 2), key
    )
    cloud_frac_i = jnp.where(first, cloud_frac_1, cloud_frac_2)
    X_u_chi_element = jnp.where(
        l_cloudy_sample,
        cloud_frac_i * chi_rand_element + (1.0 - cloud_frac_i),
        (1.0 - cloud_frac_i) * chi_rand_element,
    )
    return X_u_dp1_element, X_u_chi_element


# -----------------------------------------------------------------------------
def determine_sample_categories(
    num_samples, pdf_dim, hm_metadata,           # In
    X_nl_one_lev,                                # In
    X_mixt_comp_one_lev, importance_categories,  # In
):
    """Determines the importance category of each sample.

    Arguments:
        num_samples: Number of SILHS sample points
        pdf_dim: Number of variates in X_nl
        hm_metadata: Hydrometeor/PDF variable index metadata
        X_nl_one_lev: SILHS sample vector at one height level
        X_mixt_comp_one_lev: Mixture component label (1 or 2) of each sample point
    """
    if hm_metadata.iiPDF_rr < 0:
        raise ValueError("iiPDF_rr must be sampled for the category sampler")
    # Decode cloud, precipitation and component flags into the source category
    # order. Category indices are zero based; component labels remain 1 and 2.
    l_in_cloud = X_nl_one_lev[..., hm_metadata.iiPDF_chi] >= 0.0
    l_in_precip = X_nl_one_lev[..., hm_metadata.iiPDF_rr] > 0.0
    l_in_component_1 = X_mixt_comp_one_lev == 1
    return (
        (~l_in_precip).astype(jnp.int32) * 4
        + (~l_in_cloud).astype(jnp.int32) * 2
        + (~l_in_component_1).astype(jnp.int32)
    )
