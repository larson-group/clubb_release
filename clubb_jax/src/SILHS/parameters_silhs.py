"""Parameters and flags from parameters_silhs.F90."""

from dataclasses import dataclass, fields

# Cluster allocation strategies.
# All eight categories, effectively no clustering.
eight_cluster_allocation_opt = 1
# Four clusters for cloud/no cloud and component 1/2; precipitation is ignored.
four_cluster_allocation_opt = 2
# Two clusters: categories with either cloud or precip, and those with neither.
two_cluster_cp_nocp_opt = 3


@dataclass(frozen=True)
class eight_cluster_presc_probs_type:
    """Sample-point allocation for clustered sampling (l_lh_clustered_sampling)."""

    cloud_precip_comp1: float = 0.15
    cloud_precip_comp2: float = 0.15
    nocloud_precip_comp1: float = 0.15
    nocloud_precip_comp2: float = 0.15
    cloud_noprecip_comp1: float = 0.15
    cloud_noprecip_comp2: float = 0.15
    nocloud_noprecip_comp1: float = 0.05
    nocloud_noprecip_comp2: float = 0.05


@dataclass(frozen=True)
class silhs_config_flags_type:
    """Flags for the SILHS sampling code; immutable JAX static metadata."""

    cluster_allocation_strategy: int = two_cluster_cp_nocp_opt
    l_lh_importance_sampling: bool = True    # Do importance sampling.
    l_Lscale_vert_avg: bool = False          # Vertically average Lscale.
    l_lh_straight_mc: bool = False           # No LH or importance sampling.
    l_lh_clustered_sampling: bool = True     # Prescribed sampling with clusters.
    l_rcm_in_cloud_k_lh_start: bool = True   # Start at maximum within-cloud rcm.
    l_random_k_lh_start: bool = False        # Random start between the rcm maxima.
    l_max_overlap_in_cloud: bool = True      # Maximum vertical overlap in cloud.
    l_lh_instant_var_covar_src: bool = True  # Instantaneous var/covar tendencies.
    l_lh_limit_weights: bool = True          # Bound sample-point weights.
    l_lh_var_frac: bool = False              # Prescribe variance fractions.
    l_lh_normalize_weights: bool = True      # Weights sum to num_samples.
    l_lh_deterministic_test: bool = False    # Ordered strata/repeating draws for repeatable testing.


# Prescribed probabilities for l_lh_clustered_sampling = True.
eight_cluster_presc_probs = eight_cluster_presc_probs_type()
importance_prob_thresh = 1.0e-8  # Minimum category PDF probability for sampling.
vert_decorr_coef = 0.03          # Empirical vertical decorrelation constant [-].

# Uniform samples stay in [3.e-8, 1 - 3.e-8]: the inverse-CDF algorithm only
# has single-precision accuracy, even though it accepts double-precision input.
single_prec_thresh = 3.0e-8


# -----------------------------------------------------------------------------
def set_default_silhs_config_flags_api():
    """Set all SILHS flags to a default setting.

    Return the source output flags in their declaration order. Dataclass field
    defaults mirror the assignments in set_default_silhs_config_flags_api.
    """
    return tuple(
        getattr(silhs_config_flags_type(), f.name) for f in fields(silhs_config_flags_type)
    )


# -----------------------------------------------------------------------------
def initialize_silhs_config_flags_type_api(
    cluster_allocation_strategy,  # In
    l_lh_importance_sampling,     # In
    l_Lscale_vert_avg,            # In
    l_lh_straight_mc,             # In
    l_lh_clustered_sampling,      # In
    l_rcm_in_cloud_k_lh_start,    # In
    l_random_k_lh_start,          # In
    l_max_overlap_in_cloud,       # In
    l_lh_instant_var_covar_src,   # In
    l_lh_limit_weights,           # In
    l_lh_var_frac,                # In
    l_lh_normalize_weights,       # In
    l_lh_deterministic_test,      # In
):
    """Initialize the silhs_config_flags_type.

    Arguments:
        cluster_allocation_strategy: Eight-, four- or two-cluster allocation option
        l_lh_importance_sampling: Limit noise by performing importance sampling
        l_Lscale_vert_avg: Deprecated source flag for vertically averaged Lscale
        l_lh_straight_mc: Use true Monte Carlo sampling with no Latin hypercube sampling and
            no importance sampling
        l_lh_clustered_sampling: Use the "new" SILHS importance sampling scheme with
            prescribed probabilities
        l_rcm_in_cloud_k_lh_start: Determine k_lh_start based on maximum within-cloud rcm
        l_random_k_lh_start: Place k_lh_start at a random grid level between maximum rcm and
            maximum rcm_in_cloud
        l_max_overlap_in_cloud: Assume maximum vertical overlap when grid-box rcm exceeds
            cloud threshold
        l_lh_instant_var_covar_src: Produces "instantaneous" variance-covariance microphysical
            source terms, ignoring discretization effects
        l_lh_limit_weights: Limit SILHS sample point weights for stability
        l_lh_var_frac: Prescribe variance fractions
        l_lh_normalize_weights: Scale sample point weights to sum to num_samples (the "ratio
            estimate")
        l_lh_deterministic_test: Use deterministic sampling inputs for repeatable testing
    """
    return silhs_config_flags_type(
        cluster_allocation_strategy,
        l_lh_importance_sampling,
        l_Lscale_vert_avg,
        l_lh_straight_mc,
        l_lh_clustered_sampling,
        l_rcm_in_cloud_k_lh_start,
        l_random_k_lh_start,
        l_max_overlap_in_cloud,
        l_lh_instant_var_covar_src,
        l_lh_limit_weights,
        l_lh_var_frac,
        l_lh_normalize_weights,
        l_lh_deterministic_test,
    )


# -----------------------------------------------------------------------------
def print_silhs_config_flags_api(iunit, silhs_config_flags):
    """Prints the silhs_config_flags.

    Arguments:
        iunit: The file to write to
        silhs_config_flags: Derived type holding all configurable SILHS flags
    """
    for f in fields(silhs_config_flags):
        print(f"{f.name} = {getattr(silhs_config_flags, f.name)}", file=iunit)
