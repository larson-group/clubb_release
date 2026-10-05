"""SILHS driver from latin_hypercube_driver_module.F90.

JAX adaptations: random keys/permutation arrays are explicit, k_lh_start is
zero based, output/inout values are returns, and vertical recurrences use scan.
Fortran CUDA/MT95 generator branches are replaced by native JAX randomness.
"""

import jax
import jax.numpy as jnp
from clubb_jax.src.CLUBB_core.constants_clubb import (
    rc_tol,
    rt_tol,
    cloud_frac_min,
    chi_tol,
    eps,
)
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from clubb_jax.src.CLUBB_core.index_mapping import hydromet2pdf_idx
from clubb_jax.src.SILHS.parameters_silhs import single_prec_thresh
from clubb_jax.src.SILHS.generate_uniform_sample_module import (
    rand_uniform_real,
    generate_uniform_lh_sample,
)
from clubb_jax.src.SILHS.silhs_importance_sample_module import (
    importance_sampling_driver,
    define_importance_categories,
)
from clubb_jax.src.SILHS.transform_to_pdf_module import (
    transform_uniform_samples_to_pdf,
    chi_eta_2_rtthl,
    ltqnorm,
)
from clubb_jax.src.SILHS.math_utilities import (
    compute_sample_mean,
    compute_sample_variance,
    rand_integer_in_range,
)

l_output_2D_lognormal_dist = False
l_output_2D_uniform_dist = False
l_lh_old_cloud_weighted = False


# -----------------------------------------------------------------------------
def generate_silhs_sample(
    iter, pdf_dim, num_samples, sequence_length, nzt, ngrdcol,  # In
    l_calc_weights_all_levs_itime,                              # In
    gr, pdf_params, delta_zm, Lscale,                           # In
    lh_seed, hm_metadata,                                       # In
    mu1, mu2, sigma1, sigma2,                                   # In
    corr_cholesky_mtx_1, corr_cholesky_mtx_2,                   # In
    precip_fracs, silhs_config_flags,                           # In
    vert_decorr_coef,                                           # In
    err_info,                                                   # InOut
    stats,                                                      # InOut
    sampling_state,                                             # InOut
):
    """Generate sample points of moisture, temperature, et cetera for the purpose
    of computing tendencies with a microphysics or radiation scheme.

    Source reference: https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:overview_silhs.

    Arguments:
        iter: Model iteration (time step) number
        pdf_dim: Number of variables to sample
        num_samples: Number of samples per variable
        sequence_length: nt_repeat/num_samples; number of timesteps before sequence repeats
        nzt: Number of thermodynamic vertical model levels
        ngrdcol: Number of grid columns
        l_calc_weights_all_levs_itime: Use independent samples/weights at every level
            when true; otherwise vertically correlate the starting-level sample
        gr: Grid variable type
        pdf_params: PDF parameters [units vary]
        delta_zm: Difference in momentum altitudes [m]
        Lscale: Turbulent mixing length [m]
        lh_seed: Random number generator seed
        hm_metadata: Hydrometeor/PDF variable index metadata
        mu1: Means of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>) [units vary]
        mu2: Means of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>) [units vary]
        sigma1: Stdevs of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>) [units
            vary]
        sigma2: Stdevs of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>) [units
            vary]
        corr_cholesky_mtx_1: Correlations Cholesky matrix (1st comp.) [-]
        corr_cholesky_mtx_2: Correlations Cholesky matrix (2nd comp.) [-]
        precip_fracs: Precipitation fractions [-]
        silhs_config_flags: Flags for the SILHS sampling code [-]
        vert_decorr_coef: Empirically defined de-correlation constant [-]
        err_info: err_info struct containing err_code and err_header
        stats: JAX statistics state; updated state is returned
        sampling_state: Case-owned permutation and prior iteration; updated state is returned
    """
    if silhs_config_flags.l_lh_importance_sampling and sequence_length != 1:
        raise ValueError("Importance sampling requires sequence length one")
    if silhs_config_flags.l_Lscale_vert_avg:
        raise ValueError("l_Lscale_vert_avg is deprecated in Fortran")

    # Extra uniform variates select PDF component (dp1) and precipitation (dp2).
    d_uniform_extra = 2

    # Compute the PDF cloud-water mean used within SILHS.
    rcm_pdf = (
        pdf_params.mixt_frac * pdf_params.rc_1 + (1.0 - pdf_params.mixt_frac) * pdf_params.rc_2
    )

    # Compute the starting vertical level for Latin-hypercube sampling.
    k_lh_start = compute_k_lh_start(
        gr, nzt, ngrdcol, rcm_pdf, pdf_params,         # In
        silhs_config_flags.l_rcm_in_cloud_k_lh_start,  # In
        silhs_config_flags.l_random_k_lh_start,        # In
        lh_seed,                                       # In
    )

    # Row-wise multiply each lower triangular correlation matrix by the standard
    # deviations to obtain the covariance Cholesky factors.
    Sigma_Cholesky1 = jnp.transpose(corr_cholesky_mtx_1 * sigma1[..., None], (3, 0, 1, 2))
    Sigma_Cholesky2 = jnp.transpose(corr_cholesky_mtx_2 * sigma2[..., None], (3, 0, 1, 2))

    # Vertical correlation for arbitrary overlap from Lscale and level spacing.
    X_vert_corr = jnp.exp(-vert_decorr_coef * (gr.grid_dir * delta_zm / Lscale))
    if silhs_config_flags.l_max_overlap_in_cloud:
        X_vert_corr = jnp.where(rcm_pdf > rc_tol, 1.0, X_vert_corr)
    if clubb_at_least_debug_level(1):
        err_info = err_info.set_fatal(
            jnp.any(
                (X_vert_corr > 1.0) | (X_vert_corr < 0.0) | ~jnp.isfinite(X_vert_corr),
                axis=1,
            )
        )

    # Generate random draws, then the uniform sample and its importance weights.
    rand_pool = generate_random_pool(
        nzt, ngrdcol, pdf_dim, num_samples, d_uniform_extra,  # In
        lh_seed, gr,                                          # In
        silhs_config_flags.l_lh_deterministic_test,              # In
    )
    key = jax.random.fold_in(jax.random.PRNGKey(lh_seed), 1)
    X_u_all_levs, lh_sample_point_weights, sampling_state, l_error = generate_all_uniform_samples(
        iter, pdf_dim, d_uniform_extra, num_samples, sequence_length,  # In
        nzt, ngrdcol, k_lh_start, X_vert_corr, rand_pool,              # In
        hm_metadata.iiPDF_chi,                                         # In
        pdf_params.cloud_frac_1,                                       # In
        pdf_params.cloud_frac_2,                                       # In
        pdf_params.mixt_frac, precip_fracs,                            # In
        silhs_config_flags.cluster_allocation_strategy,                # In
        silhs_config_flags.l_lh_importance_sampling,                   # In
        silhs_config_flags.l_lh_straight_mc,                           # In
        silhs_config_flags.l_lh_clustered_sampling,                    # In
        silhs_config_flags.l_lh_limit_weights,                         # In
        silhs_config_flags.l_lh_var_frac,                              # In
        silhs_config_flags.l_lh_normalize_weights,                     # In
        silhs_config_flags.l_lh_deterministic_test,                    # In
        l_calc_weights_all_levs_itime,                                 # In
        key,                                                           # In
        sampling_state,                                                # InOut
    )

    # Determine each sample's mixture component and precipitation membership.
    # Component labels remain 1 and 2, even though array indices are zero based.
    first = X_u_all_levs[..., pdf_dim] < pdf_params.mixt_frac[:, None, :]
    X_mixt_comp_all_levs = jnp.where(first, 1, 2).astype(jnp.int32)
    cloud_frac = jnp.where(
        first, pdf_params.cloud_frac_1[:, None, :], pdf_params.cloud_frac_2[:, None, :]
    )
    precip_frac = jnp.where(
        first,
        precip_fracs.precip_frac_1[:, None, :],
        precip_fracs.precip_frac_2[:, None, :],
    )
    l_in_precip = X_u_all_levs[..., pdf_dim + 1] < precip_frac

    # Transform the uniform sample to the desired normal-lognormal PDF.
    X_nl_all_levs = transform_uniform_samples_to_pdf(
        nzt, ngrdcol, num_samples, pdf_dim, d_uniform_extra,  # In
        hm_metadata,                                          # In
        Sigma_Cholesky1, Sigma_Cholesky2,                     # In
        mu1, mu2, X_mixt_comp_all_levs,                       # In
        X_u_all_levs, cloud_frac,                             # In
        l_in_precip,                                          # In
    )

    # Accumulate diagnostics that require uniform-space sample information.
    if stats is not None and stats.l_sample:
        stats = stats_accumulate_uniform_lh(
            nzt, num_samples, ngrdcol, l_in_precip,                                      # In
            X_mixt_comp_all_levs, X_u_all_levs[..., hm_metadata.iiPDF_chi], pdf_params,  # In
            lh_sample_point_weights, k_lh_start,                                         # In
            stats,                                                                       # InOut
        )
    # Source 2D output flags are compile-time false. Host output initialization
    # remains in latin_hypercube_2D_output_api; I/O never occurs in this kernel.
    if clubb_at_least_debug_level(2):
        l_error |= jnp.any((X_u_all_levs <= 0.0) | (X_u_all_levs >= 1.0))
        l_error |= jnp.any(
            assert_consistent_cloud_frac(
                pdf_params.chi_1, pdf_params.chi_2,                # In
                pdf_params.cloud_frac_1, pdf_params.cloud_frac_2,  # In
                pdf_params.stdev_chi_1, pdf_params.stdev_chi_2,    # In
            )
        )
        l_error |= assert_correct_cloud_normal(
            num_samples,                                # In
            X_u_all_levs[..., hm_metadata.iiPDF_chi],   # In
            X_nl_all_levs[..., hm_metadata.iiPDF_chi],  # In
            X_mixt_comp_all_levs,                       # In
            pdf_params.cloud_frac_1[:, None, :],        # In
            pdf_params.cloud_frac_2[:, None, :],        # In
        )
    # Source unconditional ERROR STOP checks also report failure at negative
    # debug levels. Assertion checks above retain their source debug guards;
    # the host decides whether to stop or retain failed-column status.
    err_info = err_info.set_fatal(l_error)
    return (
        err_info,
        X_nl_all_levs,
        X_mixt_comp_all_levs,
        lh_sample_point_weights,
        stats,
        sampling_state,
    )


# -----------------------------------------------------------------------------
def generate_random_pool(
    nzt, ngrdcol, pdf_dim, num_samples, d_uniform_extra,  # In
    lh_seed, gr,                                          # In
    l_lh_deterministic_test,                              # In
):
    """Populate the (column, sample, level, variate) random pool.

    Native JAX randomness replaces the source CPU/MT95 and CUDA branches.
    The explicit timestep seed preserves reproducibility for restart callers.

    Arguments:
        nzt: Number of vertical levels
        ngrdcol: Number of grid columns
        pdf_dim: Variates
        num_samples: Number of samples
        d_uniform_extra: Uniform variates included in uniform sample but not in
            normal/lognormal sample
        lh_seed: Timestep-dependent native random seed
        gr: Grid variable type
        l_lh_deterministic_test: Use a repeating overlap pool for repeatable testing
    """
    if l_lh_deterministic_test:
        # Repeatable testing only; numerically acceptable sampling is not guaranteed.
        # Cycle through (0.125, 0.625, 0.375, 0.875), offset by the sum of
        # zero-based column, sample, level and variate indices modulo four.
        # For ngrdcol=4, the first three samples at the first level and variate
        # are (one row per sample, columns 1:4):
        #   0.125  0.625  0.375  0.875
        #   0.625  0.375  0.875  0.125
        #   0.375  0.875  0.125  0.625
        test_uniform_draws = jnp.array((0.125, 0.625, 0.375, 0.875))
        i = jnp.arange(ngrdcol, dtype=jnp.int32)[:, None, None, None]
        sample = jnp.arange(num_samples, dtype=jnp.int32)[None, :, None, None]
        k = jnp.arange(nzt, dtype=jnp.int32)[None, None, :, None]
        p = jnp.arange(pdf_dim + d_uniform_extra, dtype=jnp.int32)[None, None, None, :]
        return test_uniform_draws[(i + sample + k + p) % 4]
    return rand_uniform_real(
        jax.random.PRNGKey(lh_seed),
        (ngrdcol, num_samples, nzt, pdf_dim + d_uniform_extra),
    )


# -----------------------------------------------------------------------------
def generate_all_uniform_samples(
    iter, pdf_dim, d_uniform_extra, num_samples, sequence_length,  # In
    nzt, ngrdcol, k_lh_start, X_vert_corr, rand_pool,              # In
    iiPDF_chi,                                                     # In
    cloud_frac_1,                                                  # In
    cloud_frac_2,                                                  # In
    mixt_frac, precip_fracs,                                       # In
    cluster_allocation_strategy,                                   # In
    l_lh_importance_sampling,                                      # In
    l_lh_straight_mc,                                              # In
    l_lh_clustered_sampling,                                       # In
    l_lh_limit_weights,                                            # In
    l_lh_var_frac,                                                 # In
    l_lh_normalize_weights,                                        # In
    l_lh_deterministic_test,                                       # In
    l_calc_weights_all_levs_itime,                                 # In
    key,                                                           # In
    sampling_state,                                                # InOut
):
    """Generates uniform samples for all vertical levels, samples, and variates.
    Apply Latin-hypercube and importance sampling where configured.

    Reference: V. E. Larson and D. P. Schanen (2013), The Subgrid Importance
    Latin Hypercube Sampler (SILHS): a multivariate subcolumn generator.

    Arguments:
        iter: Model iteration number
        pdf_dim: Number of variates in CLUBB's PDF
        d_uniform_extra: Uniform variates included in uniform sample but not in
            normal/lognormal sample
        num_samples: Number of SILHS sample points
        sequence_length: Number of timesteps before new sample points are picked
        k_lh_start: Zero-based starting vertical level in each column
        X_vert_corr: Vertical correlation between adjacent levels [-]
        rand_pool: Array of randomly generated numbers
        cloud_frac_1: The PDF parameters at k_lh_start
        cloud_frac_2: The PDF parameters at k_lh_start
        mixt_frac: Weight of PDF component 1 [-]
        precip_fracs: Precipitation fractions [-]
        cluster_allocation_strategy: Strategy for distributing sample points
        l_lh_importance_sampling: Do importance sampling (SILHS)
        l_lh_straight_mc: Do not apply LH or importance sampling at all (SILHS)
        l_lh_clustered_sampling: Use prescribed probability sampling with clusters (SILHS)
        l_lh_limit_weights: Ensure weights stay under a given value
        l_lh_var_frac: Prescribe variance fractions
        l_lh_normalize_weights: Normalize weights to sum to num_samples
        l_lh_deterministic_test: Ordered strata and repeating draws for repeatable testing
        l_calc_weights_all_levs_itime: Compute independent sampling/weights at every level
        key: Native JAX random key for these draws
        sampling_state: Case-owned permutation and prior iteration; updated state is returned
    """
    # Sanity check precedes both source sampling branches.
    if l_lh_old_cloud_weighted and num_samples % 2:
        raise ValueError("Old cloud-weighted sampling requires an even sample count")

    # Straight Monte Carlo starts with independent draws and equal weights.
    if l_lh_straight_mc:
        X_u_all_levs = jnp.clip(rand_pool, single_prec_thresh, 1.0 - single_prec_thresh)
        if not l_calc_weights_all_levs_itime:
            # Generate uniform sample at other grid levels by vertically
            # correlating the starting sample, also in the straight-MC branch.
            X_u_all_levs = compute_arb_overlap(
                nzt, ngrdcol, num_samples, pdf_dim, d_uniform_extra,  # In
                k_lh_start, X_vert_corr, rand_pool,                   # In
                X_u_all_levs,                                         # InOut
            )
        return (
            X_u_all_levs,
            jnp.ones((ngrdcol, num_samples, nzt)),
            sampling_state,
            jnp.array(False),
        )
    # This scan replaces the source column loop (or level/column loops for the
    # all-level-weight option). It carries the source's shared permutation array.
    count = ngrdcol * nzt if l_calc_weights_all_levs_itime else ngrdcol

    def sample_column(sampling_state, index):
        i = index % ngrdcol
        k = index // ngrdcol if l_calc_weights_all_levs_itime else k_lh_start[i]
        column_key = jax.random.fold_in(key, index)
        draw_key, importance_key = jax.random.split(column_key)
        X_u, sampling_state = generate_uniform_lh_sample(
            iter, num_samples, sequence_length, pdf_dim + d_uniform_extra,  # In
            l_lh_deterministic_test,                                        # In
            draw_key,                                                       # In
            sampling_state,                                                 # InOut
        )

        # Without importance sampling, every point has weight one.
        weights = jnp.ones(num_samples)
        l_error = jnp.array(False)
        if l_lh_importance_sampling:
            if l_lh_old_cloud_weighted:
                from clubb_jax.src.SILHS.silhs_importance_sample_module import (
                    cloud_weighted_sampling_driver,
                )

                chi, comp, weights = cloud_weighted_sampling_driver(
                    num_samples, sampling_state.one_height_time_matrix[:, iiPDF_chi],  # In
                    sampling_state.one_height_time_matrix[:, pdf_dim],                 # In
                    cloud_frac_1[i, k], cloud_frac_2[i, k], mixt_frac[i, k],           # In
                    X_u[:, iiPDF_chi], X_u[:, pdf_dim],                                # InOut
                    importance_key,                                                    # In
                )
                X_u = X_u.at[:, iiPDF_chi].set(chi).at[:, pdf_dim].set(comp)
            else:
                chi, comp, prec, weights, l_error = importance_sampling_driver(
                    num_samples,                                                         # In
                    cloud_frac_1[i, k], cloud_frac_2[i, k],                              # In
                    mixt_frac[i, k],                                                     # In
                    precip_fracs.precip_frac_1[i, k], precip_fracs.precip_frac_2[i, k],  # In
                    cluster_allocation_strategy, l_lh_clustered_sampling,                # In
                    l_lh_limit_weights, l_lh_var_frac, l_lh_normalize_weights,           # In
                    X_u[:, iiPDF_chi], X_u[:, pdf_dim], X_u[:, pdf_dim + 1],             # InOut
                    importance_key,                                                      # In
                )
                X_u = (
                    X_u.at[:, iiPDF_chi]
                    .set(chi)
                    .at[:, pdf_dim]
                    .set(comp)
                    .at[:, pdf_dim + 1]
                    .set(prec)
                )

        # Clip uniform sample points to the expected open interval.
        return sampling_state, (
            jnp.clip(X_u, single_prec_thresh, 1.0 - single_prec_thresh),
            weights,
            l_error,
        )

    sampling_state, (samples, weights, errors) = jax.lax.scan(
        sample_column, sampling_state, jnp.arange(count)
    )

    if l_calc_weights_all_levs_itime:
        # Independent Latin-hypercube/importance sampling at every level.
        X_u_all_levs = samples.reshape(
            nzt, ngrdcol, num_samples, pdf_dim + d_uniform_extra
        ).transpose(1, 2, 0, 3)
        lh_sample_point_weights = weights.reshape(nzt, ngrdcol, num_samples).transpose(1, 2, 0)
    else:
        X_u_all_levs = (
            jnp.zeros_like(rand_pool).at[jnp.arange(ngrdcol), :, k_lh_start, :].set(samples)
        )

        # Generate the other levels by vertically correlating the starting sample.
        # See https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:vert_corr.
        X_u_all_levs = compute_arb_overlap(
            nzt, ngrdcol, num_samples, pdf_dim, d_uniform_extra,  # In
            k_lh_start, X_vert_corr, rand_pool,                   # In
            X_u_all_levs,                                         # InOut
        )
        lh_sample_point_weights = jnp.broadcast_to(weights[:, :, None], (ngrdcol, num_samples, nzt))
    return X_u_all_levs, lh_sample_point_weights, sampling_state, jnp.any(errors)


# -----------------------------------------------------------------------------
def compute_k_lh_start(
    gr, nzt, ngrdcol, rcm_pdf, pdf_params,  # In
    l_rcm_in_cloud_k_lh_start,              # In
    l_random_k_lh_start,                    # In
    lh_seed,                                # In
):
    """Determines the starting SILHS sample level

    Arguments:
        gr: Grid variable type
        nzt: Number of vertical levels
        ngrdcol: Number of grid columns
        rcm_pdf: Liquid water mixing ratio [kg/kg]
        pdf_params: PDF parameters [units vary]
        l_rcm_in_cloud_k_lh_start: Determine k_lh_start based on maximum within-cloud rcm
        l_random_k_lh_start: k_lh_start found randomly between max rcm and rcm_in_cloud
        lh_seed: Random number generator seed
    """
    # Source midpoint fallback, selecting the same physical level on either grid.
    default = nzt // 2 - 1 if gr.grid_dir_indx > 0 else nzt - nzt // 2
    cloud_frac_pdf = (
        pdf_params.mixt_frac * pdf_params.cloud_frac_1
        + (1.0 - pdf_params.mixt_frac) * pdf_params.cloud_frac_2
    )

    # Locate the maxima of within-cloud and grid-box cloud-water mixing ratio.
    rcm_in_cloud = rcm_pdf / jnp.maximum(cloud_frac_pdf, cloud_frac_min)
    k_lh_start_rcm_in_cloud = jnp.where(
        jnp.max(rcm_in_cloud, axis=1) > 0.0, jnp.argmax(rcm_in_cloud, axis=1), default
    )
    k_lh_start_rcm = jnp.where(jnp.max(rcm_pdf, axis=1) > 0.0, jnp.argmax(rcm_pdf, axis=1), default)
    if l_random_k_lh_start:
        k_lh_start = rand_integer_in_range(
            jnp.minimum(k_lh_start_rcm, k_lh_start_rcm_in_cloud),
            jnp.maximum(k_lh_start_rcm, k_lh_start_rcm_in_cloud),

            # Keep the native height draw independent of the random pool
            # (seed key) and Latin-hypercube draws (fold-in tag 1).
            jax.random.fold_in(jax.random.PRNGKey(lh_seed), 2),
        )

        # Reflect the randomized index so opposite grid directions select the
        # same physical level, as in the source comparison convention.
        if gr.grid_dir_indx < 0:
            k_lh_start = k_lh_start_rcm_in_cloud + k_lh_start_rcm - k_lh_start
        return k_lh_start
    return k_lh_start_rcm_in_cloud if l_rcm_in_cloud_k_lh_start else k_lh_start_rcm


# -----------------------------------------------------------------------------
def clip_transform_silhs_output(
    nzt, ngrdcol, num_samples,           # In
    pdf_dim, hydromet_dim, hm_metadata,  # In
    X_mixt_comp_all_levs,                # In
    X_nl_all_levs,                       # InOut
    pdf_params, l_use_Ncn_to_Nc,         # In
):
    """Derives from the SILHS sampling structure X_nl_all_levs the variables
    rt, thl, rc, rv, and Nc, for all sample points and height levels.

    Arguments:
        nzt: Number of vertical levels
        ngrdcol: Number of grid columns
        num_samples: Number of SILHS sample points
        pdf_dim: Number of variates in X_nl
        hydromet_dim: Number of hydrometeor species
        hm_metadata: Hydrometeor/PDF variable index metadata
        X_mixt_comp_all_levs: Which component this sample is in (1 or 2)
        X_nl_all_levs: SILHS sample points [units vary]
        pdf_params: The PDF parameters
        l_use_Ncn_to_Nc: Whether to call Ncn_to_Nc (.true.) or not (.false.); Ncn_to_Nc might
            cause problems with the MG microphysics since the changes made here (Nc-tendency)
            are not fed into the microphysics
    """
    # Calculate rt and thl from chi and eta (source reference: CLUBB ticket 751).
    lh_rt_clipped, lh_thl_clipped = chi_eta_2_rtthl(
        nzt, ngrdcol, num_samples,                                                             # In
        pdf_params.rt_1, pdf_params.thl_1,                                                     # In
        pdf_params.rt_2, pdf_params.thl_2,                                                     # In
        pdf_params.crt_1, pdf_params.cthl_1,                                                   # In
        pdf_params.crt_2, pdf_params.cthl_2,                                                   # In
        pdf_params.chi_1, pdf_params.chi_2,                                                    # In
        X_nl_all_levs[..., hm_metadata.iiPDF_chi], X_nl_all_levs[..., hm_metadata.iiPDF_eta],  # In
        X_mixt_comp_all_levs,                                                                  # In
    )

    # Clip rt, then rc = chi * H(chi); retain at least rt_tol vapor.
    lh_rt_clipped = jnp.maximum(lh_rt_clipped, rt_tol)
    chi = X_nl_all_levs[..., hm_metadata.iiPDF_chi]
    lh_rc_clipped = jnp.minimum(jnp.maximum(chi, 0.0), lh_rt_clipped - rt_tol)

    # Vapor is the residual total water after cloud liquid.
    lh_rv_clipped = lh_rt_clipped - lh_rc_clipped
    lh_Nc_clipped = X_nl_all_levs[..., hm_metadata.iiPDF_Ncn]

    # Nc = Ncn * H(chi) when the cloud-nuclei-to-droplet conversion is enabled.
    if l_use_Ncn_to_Nc:
        lh_Nc_clipped = jnp.where(chi > 0.0, lh_Nc_clipped, 0.0)
    # Fortran l_clip_hydromet_samples is a compile-time false constant.
    return (
        X_nl_all_levs,
        lh_rt_clipped,
        lh_thl_clipped,
        lh_rc_clipped,
        lh_rv_clipped,
        lh_Nc_clipped,
    )


# -----------------------------------------------------------------------------
def assert_consistent_cloud_frac(
    chi_1, chi_2,                # In
    cloud_frac_1, cloud_frac_2,  # In
    stdev_chi_1, stdev_chi_2,    # In
):
    """Performs an assertion check that cloud_frac_i is consistent with chi_i and
    stdev_chi_i in pdf_params for each PDF component.

    Arguments:
        chi_1: Mean chi in PDF component 1 [kg/kg]
        chi_2: Mean chi in PDF component 2 [kg/kg]
        cloud_frac_1: Cloud fraction in PDF component 1 [-]
        cloud_frac_2: Cloud fraction in PDF component 2 [-]
        stdev_chi_1: Standard deviation of chi in PDF component 1 [kg/kg]
        stdev_chi_2: Standard deviation of chi in PDF component 2 [kg/kg]
    """
    return assert_consistent_cf_component(
        chi_1, stdev_chi_1, cloud_frac_1
    ) | assert_consistent_cf_component(chi_2, stdev_chi_2, cloud_frac_2)


# -----------------------------------------------------------------------------
def assert_consistent_cf_component(mu_chi_i, sigma_chi_i, cloud_frac_i):
    """Performs an assertion check that cloud_frac_i is consistent with chi_i and
    stdev_chi_i for a PDF component.
    The SILHS sample generation process relies on precisely a cloud_frac
    amount of mass in the cloudy portion of the PDF of chi, that is, where
    chi > 0. In other words, the probability that chi > 0 should be exactly
    cloud_frac.
    Stated even more mathematically, CDF_chi(0) = 1 - cloud_frac, where
    CDF_chi is the cumulative distribution function of chi. This can be
    expressed as invCDF_chi(1 - cloud_frac) = zero.
    This subroutine uses ltqnorm, which is apparently a fancy name for the
    inverse cumulative distribution function of the standard normal
    distribution.

    Arguments:
        mu_chi_i: Mean of chi in a PDF component
        sigma_chi_i: Standard deviation of chi in a PDF component
        cloud_frac_i: Cloud fraction in a PDF component
    """
    # Match the special zero-cloud condition in calc_cloud_frac_component.
    # This check must change if that condition changes in PDF closure.
    omitted = (jnp.abs(mu_chi_i) <= eps) & (sigma_chi_i <= chi_tol)

    # Check either end of the cloud-fraction tolerance box only where ltqnorm
    # accepts the argument; clipped masked inputs keep JAX branches finite.
    left = cloud_frac_i - 5.0e-6
    right = cloud_frac_i + 5.0e-6
    left_chi = ltqnorm(1.0 - jnp.clip(left, 1.0e-5, 1.0 - 1.0e-5)) * sigma_chi_i + mu_chi_i
    right_chi = ltqnorm(1.0 - jnp.clip(right, 1.0e-5, 1.0 - 1.0e-5)) * sigma_chi_i + mu_chi_i
    return ~omitted & (
        ((left >= 1.0e-5) & (left_chi <= 0.0)) | ((right <= 1.0 - 1.0e-5) & (right_chi > 0.0))
    )


# -----------------------------------------------------------------------------
def assert_correct_cloud_normal(
    num_samples,   # In
    X_u_chi,       # In
    X_nl_chi,      # In
    X_mixt_comp,   # In
    cloud_frac_1,  # In
    cloud_frac_2,  # In
):
    """Asserts that all SILHS sample points that are in cloud in uniform space
    are in cloud in normal space, and that all SILHS sample points that are
    in clear air in uniform space are in clear air in normal space.

    Arguments:
        num_samples: Number of SILHS sample points
        X_u_chi: Samples of chi in uniform space
        X_nl_chi: Samples of chi in normal space
        X_mixt_comp: PDF component of each sample
        cloud_frac_1: Cloud fraction in PDF component 1
        cloud_frac_2: Cloud fraction in PDF component 2
    """
    cloud_frac_i = jnp.where(X_mixt_comp == 1, cloud_frac_1, cloud_frac_2)
    clear = X_u_chi < 1.0 - cloud_frac_i
    return jnp.any(
        jnp.where(
            clear,
            X_nl_chi > 1000.0 * single_prec_thresh,
            X_nl_chi <= -1000.0 * single_prec_thresh,
        )
        | (X_u_chi >= 1.0)
    )


# -----------------------------------------------------------------------------
def latin_hypercube_2D_output_api(
    nzt, zt, pdf_dim, num_samples, hm_metadata,  # In
    stats, err_info,                             # InOut
):
    """Create/open optional SILHS 2D output through the host statistics writer.

    Arguments:
        hm_metadata: Hydrometeor/PDF variable index metadata
        stats: JAX statistics state; updated state is returned
        err_info: err_info struct containing err_code and err_header
    """
    # Host-only API; both source output flags are false in the normal build.
    if stats is None:
        return stats, err_info
    names = [""] * pdf_dim
    for name in (
        "chi",
        "eta",
        "w",
        "rr",
        "ri",
        "rs",
        "rg",
        "Nr",
        "Ncn",
        "Ni",
        "Ns",
        "Ng",
    ):
        index = getattr(hm_metadata, "iiPDF_" + name)
        if index >= 0:
            names[index] = name
    stats.stats_lh_samples_init(
        num_samples,
        nzt,
        names if l_output_2D_lognormal_dist else (),
        (
            names + ["dp1", "dp2", "X_mixt_comp", "lh_sample_point_weights"]
            if l_output_2D_uniform_dist
            else ()
        ),
        zt,
    )
    return stats, err_info


# -----------------------------------------------------------------------------
def compute_arb_overlap(
    nzt, ngrdcol, num_samples, pdf_dim, d_uniform_extra,  # In
    k_lh_start, vert_corr, rand_pool,                     # In
    X_u_all_levs,                                         # InOut
):
    """Re-computes X_u (uniform sample) using an arbitrary correlation specified
    by X_vert_corr (which can vary with height).
    This is an improved algorithm that doesn't require us to convert from a
    unifrom distribution to a Gaussian distribution and back again.

    Arguments:
        nzt: Vertical levels
        ngrdcol: Columns
        num_samples: Number of SILHS sample points
        pdf_dim: Number of variates in CLUBB's PDF
        d_uniform_extra: Uniform variates included in uniform sample but not in
            normal/lognormal sample
        k_lh_start: Zero-based starting vertical level in each column
        vert_corr: Vertical correlation between k points in range [0,1] [-]
        rand_pool: Array of randomly generated numbers
        X_u_all_levs: Uniform variates at all levels [-]. Starting-level values must
            already be populated; the overlap recurrence fills the other levels.
    """
    # Recompute upward and downward from the starting level for all variates.
    # Each scan carries the previous level's sample, as in the source loops.
    start = X_u_all_levs[jnp.arange(ngrdcol), :, k_lh_start, :]

    def step(unbounded_point, k):
        active = k > k_lh_start
        half_width = 1.0 - vert_corr[:, k]
        candidate = (
            unbounded_point
            - half_width[:, None, None]
            + 2.0 * half_width[:, None, None] * rand_pool[:, :, k, :]
        )

        # Fold points outside [single_prec_thresh, 1 - single_prec_thresh]
        # back into the valid interval.
        candidate = jnp.where(
            candidate > 1.0 - single_prec_thresh,
            2.0 - candidate - 2.0 * single_prec_thresh,
            jnp.where(
                candidate < single_prec_thresh,
                -candidate + 2.0 * single_prec_thresh,
                candidate,
            ),
        )
        unbounded_point = jnp.where(active[:, None, None], candidate, unbounded_point)
        return unbounded_point, unbounded_point

    _, upper = jax.lax.scan(step, start, jnp.arange(nzt))

    def step_down(unbounded_point, k):
        active = k < k_lh_start
        half_width = 1.0 - vert_corr[:, k]
        candidate = (
            unbounded_point
            - half_width[:, None, None]
            + 2.0 * half_width[:, None, None] * rand_pool[:, :, k, :]
        )
        candidate = jnp.where(
            candidate > 1.0 - single_prec_thresh,
            2.0 - candidate - 2.0 * single_prec_thresh,
            jnp.where(
                candidate < single_prec_thresh,
                -candidate + 2.0 * single_prec_thresh,
                candidate,
            ),
        )
        unbounded_point = jnp.where(active[:, None, None], candidate, unbounded_point)
        return unbounded_point, unbounded_point

    # Downward recurrence uses the same reflection and random displacement.
    _, lower = jax.lax.scan(step_down, start, jnp.arange(nzt - 1, -1, -1))
    upper = upper.transpose(1, 2, 0, 3)
    lower = lower[::-1].transpose(1, 2, 0, 3)
    return jnp.where(
        jnp.arange(nzt)[None, None, :, None] >= k_lh_start[:, None, None, None],
        upper,
        lower,
    )


# -----------------------------------------------------------------------------
def stats_accumulate_lh_api(
    gr, nzt, ngrdcol, num_samples, pdf_dim, rho_ds_zt,  # In
    hydromet_dim, hm_metadata,                          # In
    lh_sample_point_weights, X_nl_all_levs,             # In
    lh_rt_clipped, lh_thl_clipped,                      # In
    lh_rc_clipped, lh_rv_clipped,                       # In
    lh_Nc_clipped,                                      # In
    stats,                                              # InOut
):
    """Clip subcolumns from latin hypercube and create stats for diagnostic
    purposes.

    Arguments:
        gr: Grid metadata and coordinate arrays
        nzt: Number of vertical model levels
        ngrdcol: Number of model columns
        num_samples: Number of calls to microphysics per timestep (normally=2)
        pdf_dim: Number of variables to sample
        rho_ds_zt: Dry, static density (thermo. levs.) [kg/m^3]
        hydromet_dim: Number of hydrometeor species
        hm_metadata: Hydrometeor/PDF variable index metadata
        X_nl_all_levs: Sample that is transformed ultimately to normal-lognormal
        lh_rt_clipped: rt generated from silhs sample points
        lh_thl_clipped: thl generated from silhs sample points
        lh_rc_clipped: rc generated from silhs sample points
        lh_rv_clipped: rv generated from silhs sample points
        lh_Nc_clipped: Nc generated from silhs sample points
        stats: JAX statistics state; updated state is returned
    """
    if stats is None or not stats.l_sample:
        return stats
    # Weighted sample means; weights are one when importance sampling is off.
    lh_rcm = compute_sample_mean(
        nzt, num_samples, ngrdcol,               # In
        lh_sample_point_weights, lh_rc_clipped,  # In
    )
    stats = stats.update("lh_rcm", lh_rcm)
    stats = stats.update("lh_lwp", jnp.sum(rho_ds_zt * lh_rcm * gr.dzt, axis=1))
    weights_sum = jnp.sum(lh_sample_point_weights, axis=(1, 2))
    stats = stats.update("lh_sample_weights_sum", weights_sum)
    stats = stats.update("lh_sample_weights_avg", weights_sum / (num_samples * nzt))
    lh_thlm = compute_sample_mean(
        nzt, num_samples, ngrdcol,                # In
        lh_sample_point_weights, lh_thl_clipped,  # In
    )
    stats = stats.update("lh_thlm", lh_thlm)
    lh_rvm = compute_sample_mean(
        nzt, num_samples, ngrdcol,               # In
        lh_sample_point_weights, lh_rv_clipped,  # In
    )
    stats = stats.update("lh_rvm", lh_rvm)
    stats = stats.update("lh_vwp", jnp.sum(rho_ds_zt * lh_rvm * gr.dzt, axis=1))
    w = X_nl_all_levs[..., hm_metadata.iiPDF_w]
    lh_wm = compute_sample_mean(nzt, num_samples, ngrdcol, lh_sample_point_weights, w)
    stats = stats.update("lh_wm", lh_wm)
    lh_hydromet = jnp.zeros((ngrdcol, nzt, hydromet_dim))
    hydromet_all_points, Ncn_all_points = copy_X_nl_into_hydromet_all_pts(
        nzt, pdf_dim, num_samples, ngrdcol, X_nl_all_levs,  # In
        hydromet_dim, hm_metadata, lh_hydromet,             # In
    )
    lh_hydromet = compute_sample_mean(
        nzt, num_samples, ngrdcol,                                # In
        lh_sample_point_weights[..., None], hydromet_all_points,  # In
    )
    lh_Ncnm = compute_sample_mean(
        nzt, num_samples, ngrdcol,                # In
        lh_sample_point_weights, Ncn_all_points,  # In
    )
    stats = stats.update("lh_Ncnm", lh_Ncnm)
    lh_Ncm = compute_sample_mean(
        nzt, num_samples, ngrdcol,               # In
        lh_sample_point_weights, lh_Nc_clipped,  # In
    )
    stats = stats.update("lh_Ncm", lh_Ncm)

    # Weighted and unweighted Latin-hypercube estimates of cloud fraction.
    cloud = (X_nl_all_levs[..., hm_metadata.iiPDF_chi] > 0.0).astype(jnp.float64)
    stats = stats.update(
        "lh_cloud_frac",
        compute_sample_mean(nzt, num_samples, ngrdcol, lh_sample_point_weights, cloud),
    )
    stats = stats.update("lh_cloud_frac_unweighted", jnp.mean(cloud, axis=1))

    # Estimate chi, its variance, and eta.
    chi = X_nl_all_levs[..., hm_metadata.iiPDF_chi]
    lh_chi = compute_sample_mean(nzt, num_samples, ngrdcol, lh_sample_point_weights, chi)
    stats = stats.update("lh_chi", lh_chi)
    stats = stats.update(
        "lh_chip2",
        compute_sample_variance(
            nzt, num_samples, ngrdcol,             # In
            chi, lh_sample_point_weights, lh_chi,  # In
        ),
    )
    stats = stats.update(
        "lh_eta",
        compute_sample_mean(
            nzt, num_samples, ngrdcol,                                           # In
            lh_sample_point_weights, X_nl_all_levs[..., hm_metadata.iiPDF_eta],  # In
        ),
    )

    # Variances of velocity, cloud/total water and liquid potential temperature.
    for name, samples, mean in (
        ("lh_wp2_zt", w, lh_wm),
        ("lh_rcp2_zt", lh_rc_clipped, lh_rcm),
        ("lh_rtp2_zt", lh_rt_clipped, lh_rvm + lh_rcm),
        ("lh_thlp2_zt", lh_thl_clipped, lh_thlm),
    ):
        stats = stats.update(
            name,
            compute_sample_variance(
                nzt, num_samples, ngrdcol,               # In
                samples, lh_sample_point_weights, mean,  # In
            ),
        )
    if hm_metadata.iirr >= 0:
        stats = stats.update(
            "lh_rrp2_zt",
            compute_sample_variance(
                nzt, num_samples, ngrdcol,                                            # In
                hydromet_all_points[..., hm_metadata.iirr], lh_sample_point_weights,  # In
                lh_hydromet[..., hm_metadata.iirr],                                   # In
            ),
        )
    if hm_metadata.iiPDF_Ncn >= 0:
        stats = stats.update(
            "lh_Ncnp2_zt",
            compute_sample_variance(
                nzt, num_samples, ngrdcol,                         # In
                Ncn_all_points, lh_sample_point_weights, lh_Ncnm,  # In
            ),
        )
    stats = stats.update(
        "lh_Ncp2_zt",
        compute_sample_variance(
            nzt, num_samples, ngrdcol,                       # In
            lh_Nc_clipped, lh_sample_point_weights, lh_Ncm,  # In
        ),
    )
    if hm_metadata.iiNr >= 0:
        stats = stats.update(
            "lh_Nrp2_zt",
            compute_sample_variance(
                nzt, num_samples, ngrdcol,                                            # In
                hydromet_all_points[..., hm_metadata.iiNr], lh_sample_point_weights,  # In
                lh_hydromet[..., hm_metadata.iiNr],                                   # In
            ),
        )

    # Diagnostic averages of the sample points fed to microphysics.
    for name, index in (
        ("rrm", hm_metadata.iirr),
        ("rim", hm_metadata.iiri),
        ("rsm", hm_metadata.iirs),
        ("rgm", hm_metadata.iirg),
        ("Nrm", hm_metadata.iiNr),
        ("Nim", hm_metadata.iiNi),
        ("Nsm", hm_metadata.iiNs),
        ("Ngm", hm_metadata.iiNg),
    ):
        if index >= 0:
            stats = stats.update("lh_" + name, lh_hydromet[..., index])
    return stats


# -----------------------------------------------------------------------------
def stats_accumulate_uniform_lh(
    nzt, num_samples, ngrdcol, l_in_precip_all_levs,     # In
    X_mixt_comp_all_levs, X_u_chi_all_levs, pdf_params,  # In
    lh_sample_point_weights, k_lh_start,                 # In
    stats,                                               # InOut
):
    """Samples statistics that cannot be deduced from the normal-lognormal
    SILHS sample (X_nl_all_levs)

    Arguments:
        nzt: Number of vertical levels
        num_samples: Number of SILHS sample points
        ngrdcol: Number of grid columns
        l_in_precip_all_levs: Boolean variables indicating whether a sample is in
            precipitation at a given height level
        X_mixt_comp_all_levs: Integers indicating which mixture component a sample is in at a
            given height level
        X_u_chi_all_levs: Uniform value of chi
        pdf_params: The PDF parameters
        lh_sample_point_weights: The weight of each sample
        k_lh_start: Zero-based starting level for preferential cloud sampling
        stats: JAX statistics state; updated state is returned
    """
    if not stats.l_sample:
        return stats
    int_in_precip = l_in_precip_all_levs.astype(jnp.float64)
    int_mixt_comp = (X_mixt_comp_all_levs == 1).astype(jnp.float64)
    stats = stats.update(
        "lh_precip_frac",
        compute_sample_mean(
            nzt, num_samples, ngrdcol,               # In
            lh_sample_point_weights, int_in_precip,  # In
        ),
    )
    stats = stats.update("lh_precip_frac_unweighted", jnp.mean(int_in_precip, axis=1))
    stats = stats.update(
        "lh_mixt_frac",
        compute_sample_mean(
            nzt, num_samples, ngrdcol,               # In
            lh_sample_point_weights, int_mixt_comp,  # In
        ),
    )
    stats = stats.update("lh_mixt_frac_unweighted", jnp.mean(int_mixt_comp, axis=1))

    # Diagnostics retain the source one-based level, represented as a real.
    stats = stats.update("k_lh_start", (k_lh_start + 1).astype(jnp.float64))
    cloud_frac_i = jnp.where(
        X_mixt_comp_all_levs == 1,
        pdf_params.cloud_frac_1[:, None, :],
        pdf_params.cloud_frac_2[:, None, :],
    )
    l_in_cloud = X_u_chi_all_levs > 1.0 - cloud_frac_i
    category = (
        (~l_in_precip_all_levs).astype(jnp.int32) * 4
        + (~l_in_cloud).astype(jnp.int32) * 2
        + (X_mixt_comp_all_levs != 1).astype(jnp.int32)
    )
    for icategory in range(8):
        # Microphysics is not run at the lower boundary level.
        frac = jnp.mean(category == icategory, axis=1).at[:, 0].set(0.0)
        stats = stats.update(f"lh_samp_frac_{icategory+1}", frac)
    return stats


# -----------------------------------------------------------------------------
def copy_X_nl_into_hydromet_all_pts(
    nzt, pdf_dim, num_samples, ngrdcol, X_nl_all_levs,  # In
    hydromet_dim, hm_metadata, hydromet,                # In
):
    """Copy the points from the latin hypercube sample to an array with just the
    hydrometeors, for every model column.

    Arguments:
        hydromet_dim: Number of hydrometeor species
        hm_metadata: Hydrometeor/PDF variable index metadata
        hydromet: Hydrometeor species [units vary]
    """
    # Use mean fields for unsampled species; overwrite with PDF sample values
    # for sampled hydrometeor mixing ratios and number concentrations.
    hydromet_all_points = jnp.broadcast_to(
        hydromet[:, None], (ngrdcol, num_samples, nzt, hydromet_dim)
    )
    for ivar in range(hydromet_dim):
        pdf_hydromet_idx = hydromet2pdf_idx(ivar, hm_metadata)
        if pdf_hydromet_idx >= 0:
            hydromet_all_points = hydromet_all_points.at[..., ivar].set(
                X_nl_all_levs[..., pdf_hydromet_idx]
            )
    Ncn_all_points = (
        X_nl_all_levs[..., hm_metadata.iiPDF_Ncn]
        if hm_metadata.iiPDF_Ncn >= 0
        else jnp.zeros((ngrdcol, num_samples, nzt))
    )
    return hydromet_all_points, Ncn_all_points
