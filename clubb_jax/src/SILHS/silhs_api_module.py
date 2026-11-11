"""Public SILHS API from silhs_api_module.F90.

Sampling storage is an explicit required JAX argument and return; output-only
Fortran arguments are returned. Random generators are native JAX operations.
"""

from clubb_jax.src.SILHS.parameters_silhs import (
    silhs_config_flags_type,
    set_default_silhs_config_flags_api,
    initialize_silhs_config_flags_type_api,
    print_silhs_config_flags_api,
)
from clubb_jax.src.SILHS.latin_hypercube_driver_module import (
    generate_silhs_sample,
    clip_transform_silhs_output,
    latin_hypercube_2D_output_api,
    stats_accumulate_lh_api,
)
from clubb_jax.src.SILHS.est_kessler_microphys_module import est_kessler_microphys_api
from clubb_jax.src.SILHS.lh_microphys_var_covar_module import (
    lh_microphys_var_covar_driver_api,
)


# -----------------------------------------------------------------------------
def generate_silhs_sample_api(
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
    """Generate SILHS samples for all model columns, including ngrdcol=1.

    The Fortran API owns an OpenACC data region. JAX owns device placement;
    sampling_state explicitly carries the source threadprivate permutation.
    Return err_info, normal/lognormal samples, component IDs, sample weights,
    stats, and updated sampling_state, in that order.

    Arguments:
        iter: Model iteration number
        pdf_dim: Number of variables to sample
        num_samples: Number of samples per variable
        sequence_length: nt_repeat/num_samples; number of timesteps before sequence repeats.
        nzt: Number of vertical model levels
        ngrdcol: Number of grid columns
        l_calc_weights_all_levs_itime: Use independent samples/weights at every level
            when true; otherwise vertically correlate the starting-level sample
        gr: Grid variable type
        pdf_params: PDF parameters [units vary]
        delta_zm: Difference in moment. altitudes [m]
        Lscale: Turbulent Mixing Length [m]
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
        silhs_config_flags: Static configuration flags for SILHS sampling
        vert_decorr_coef: Empirically defined de-correlation constant [-]
        err_info: err_info struct containing err_code and err_header
        stats: JAX statistics state; updated state is returned
        sampling_state: Case-owned permutation and prior iteration; updated state is returned
    """
    return generate_silhs_sample(
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
    )


# -----------------------------------------------------------------------------
def clip_transform_silhs_output_api(
    nzt, ngrdcol, num_samples,           # In
    pdf_dim, hydromet_dim, hm_metadata,  # In
    X_mixt_comp_all_levs,                # In
    X_nl_all_levs,                       # InOut
    pdf_params, l_use_Ncn_to_Nc,         # In
):
    """Compute clipped SILHS fields for all columns, including ngrdcol=1.

    JAX adaptation: return the source inout X_nl_all_levs and output rt, thl,
    rc, rv, Nc arrays instead of mutating arrays inside an OpenACC data region.

    Arguments:
        nzt: Number of vertical levels
        ngrdcol: Number of grid columns
        num_samples: Number of SILHS sample points
        pdf_dim: Number of variates in X_nl_one_lev
        hydromet_dim: Number of hydrometeor species
        hm_metadata: Hydrometeor/PDF variable index metadata
        X_mixt_comp_all_levs: Which component this sample is in (1 or 2)
        X_nl_all_levs: SILHS sample points [units vary]
        pdf_params: The PDF parameters
    """
    return clip_transform_silhs_output(
        nzt, ngrdcol, num_samples,           # In
        pdf_dim, hydromet_dim, hm_metadata,  # In
        X_mixt_comp_all_levs,                # In
        X_nl_all_levs,                       # InOut
        pdf_params, l_use_Ncn_to_Nc,         # In
    )
