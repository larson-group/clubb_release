"""Hydrometeor PDF preparation from pdf_hydromet_microphys_wrapper.F90.

JAX adaptation: outputs and stats are functional returns, columns are batched. Sampling storage
is an explicit JAX input/return in place of threadprivate arrays.
"""

from functools import partial

import jax
import jax.numpy as jnp
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from clubb_jax.src.CLUBB_core.setup_clubb_pdf_params import setup_pdf_parameters_api
from clubb_jax.src.Microphys.mixed_moment_PDF_integrals import hydrometeor_mixed_moments
from clubb_jax.src.Microphys import parameters_microphys
from clubb_jax.src.CLUBB_core.hydromet_pdf_parameter_module import init_hydromet_pdf_params


# -----------------------------------------------------------------------------
@partial(
    jax.jit,
    static_argnames=(
        'ngrdcol', 'pdf_dim', 'hydromet_dim',
        'clubb_config_flags', 'silhs_config_flags', 'l_rad_itime',
    ),
)
def pdf_hydromet_microphys_prep(
    gr, ngrdcol, pdf_dim, hydromet_dim,               # In
    itime, vert_decorr_coef,                          # In
    Nc_in_cloud, cloud_frac, ice_supersat_frac,       # In
    rho_ds_zt, Lscale, Kh_zm, hydromet, wphydrometp,  # In
    corr_array_n_cloud, corr_array_n_below,           # In
    hm_metadata, pdf_params, clubb_params,            # In
    clubb_config_flags, silhs_config_flags,           # In
    l_rad_itime,                                      # In
    stats,                                            # InOut
    err_info,                                         # InOut
    precip_fracs,                                     # InOut
    sampling_state,                                   # InOut
):
    """Prepare the hydrometeor PDF, mixed moments and optional SILHS subcolumns.

    sampling_state carries the case-owned permutation between timesteps. Return updated
    statistics/error state, PDF fields, samples/clipped fields and sampling state. Fatal boundaries
    follow the source early returns.

    Arguments:
        gr: Grid metadata and coordinate arrays
        ngrdcol: Number of grid columns
        pdf_dim: Number of variables in the correlation array
        hydromet_dim: Number of hydrometeor species
        itime: Current model timestep index, used by the sampling sequence.
        vert_decorr_coef: Empirically defined de-correlation constant [-]
        Nc_in_cloud: Mean (in-cloud) cloud droplet conc. [num/kg]
        cloud_frac: Cloud fraction [-]
        ice_supersat_frac: Ice supersaturation fraction [-]
        rho_ds_zt: Dry, base-state density on thermo. levs. [kg/m^3]
        Lscale: Turbulent Mixing Length [m]
        Kh_zm: Eddy diffusivity coef. on momentum levels [m^2/s]
        hydromet: Mean of hydrometeor, hm (overall) (t-levs.) [units]
        wphydrometp: Covariance < w'h_m' > (momentum levels) [(m/s)units]
        corr_array_n_cloud: Prescribed normal space corr. array in cloud [-]
        corr_array_n_below: Prescribed normal space corr. array below cl. [-]
        hm_metadata: Hydrometeor/PDF variable index metadata
        pdf_params: PDF parameters [units vary]
        clubb_params: Array of CLUBB's tunable parameters [units vary]
        clubb_config_flags: Derived type holding all configurable CLUBB flags
        silhs_config_flags: Static configuration flags for SILHS sampling
        l_rad_itime: Source radiation-timestep flag, retained by the microphysics interface.
        stats: JAX statistics state; updated state is returned
        err_info: err_info struct containing err_code and err_header
        precip_fracs: Precipitation fractions [-]
        sampling_state: Case-owned permutation and prior iteration; updated state is returned
    """
    # Setup the PDF parameters.
    hydromet_pdf_params = init_hydromet_pdf_params()

    # Pure outputs are defined even on the source's inactive/early-return paths.
    hydrometp2 = jnp.zeros((ngrdcol, gr.nzm, hydromet_dim))
    mu_x_1_n = mu_x_2_n = sigma_x_1_n = sigma_x_2_n = jnp.zeros((ngrdcol, gr.nzt, pdf_dim))
    corr_array_1_n = corr_array_2_n = corr_cholesky_mtx_1 = corr_cholesky_mtx_2 = jnp.zeros(
        (ngrdcol, gr.nzt, pdf_dim, pdf_dim)
    )
    rtphmp_zt = thlphmp_zt = wp2hmp = jnp.zeros_like(hydromet)
    if parameters_microphys.microphys_scheme != "none":
        # The source standalone passes its (ngrdcol,nparams) array to this
        # wrapper's rank-one (nparams) dummy. Preserve Fortran sequence association:
        # the PDF uses the first nparams elements in column-major storage order.
        # A coordinated source/API change is needed before using per-column values.
        clubb_params = jnp.broadcast_to(
            clubb_params.T.reshape(-1)[: clubb_params.shape[1]], clubb_params.shape
        )
        (
            err_info,
            hydrometp2,
            mu_x_1_n,
            mu_x_2_n,
            sigma_x_1_n,
            sigma_x_2_n,
            corr_array_1_n,
            corr_array_2_n,
            corr_cholesky_mtx_1,
            corr_cholesky_mtx_2,
            precip_fracs,
            hydromet_pdf_params,
            stats,
        ) = setup_pdf_parameters_api(
            gr, gr.nzm, gr.nzt, ngrdcol, pdf_dim,             # In
            hydromet_dim,                                     # In
            Nc_in_cloud, cloud_frac, Kh_zm,                   # In
            ice_supersat_frac, hydromet, wphydrometp,         # In
            corr_array_n_cloud, corr_array_n_below,           # In
            hm_metadata,                                      # In
            pdf_params,                                       # In
            clubb_params,                                     # In
            clubb_config_flags.iiPDF_type,                    # In
            clubb_config_flags.l_use_precip_frac,             # In
            clubb_config_flags.l_diagnose_correlations,       # In
            clubb_config_flags.l_calc_w_corr,                 # In
            clubb_config_flags.l_const_Nc_in_cloud,           # In
            clubb_config_flags.l_fix_w_chi_eta_correlations,  # In
            err_info,                                         # InOut
            precip_fracs,                                     # InOut
            hydromet_pdf_params,                              # Out
            stats,                                            # InOut
        )

        # JAX adaptation: the source returns on fatal PDF setup before mixed
        # moments and their statistics. Undefined pure outputs are returned as zero.
        rtphmp_zt, thlphmp_zt, wp2hmp, stats = jax.lax.cond(
            err_info.any_fatal() & clubb_at_least_debug_level(0),
            lambda _: (
                jnp.zeros_like(hydromet),
                jnp.zeros_like(hydromet),
                jnp.zeros_like(hydromet),
                stats,
            ),
            lambda _: hydrometeor_mixed_moments(
                gr, gr.ngrdcol, gr.nzt, pdf_dim, hydromet_dim,  # In
                hydromet, hm_metadata,                          # In
                mu_x_1_n, mu_x_2_n,                             # In
                sigma_x_1_n, sigma_x_2_n,                       # In
                corr_array_1_n, corr_array_2_n,                 # In
                pdf_params, hydromet_pdf_params,                # In
                precip_fracs,                                   # In
                stats,                                          # InOut
            ),
            operand=None,
        )
    # Compute subcolumns if enabled, following source call order.
    if parameters_microphys.lh_microphys_type != parameters_microphys.lh_microphys_disabled:
        from clubb_jax.src.SILHS.silhs_api_module import (
            generate_silhs_sample_api,
            clip_transform_silhs_output_api,
        )
        from clubb_jax.src.SILHS.latin_hypercube_driver_module import stats_accumulate_lh_api

        def compute_subcolumns(_):
            # Timestep-dependent seed: source convention supports repeatable restarts.
            lh_seed_custom = jnp.asarray(parameters_microphys.lh_seed * itime, dtype=jnp.uint32)
            sample_err, X_nl, X_comp, weights, sample_stats, next_state = generate_silhs_sample_api(
                itime, pdf_dim, parameters_microphys.lh_num_samples,       # In
                parameters_microphys.lh_sequence_length, gr.nzt, ngrdcol,  # In
                False,                                                     # In
                gr, pdf_params, gr.dzt, Lscale,                            # In
                lh_seed_custom, hm_metadata,                               # In
                mu_x_1_n, mu_x_2_n, sigma_x_1_n, sigma_x_2_n,              # In
                corr_cholesky_mtx_1, corr_cholesky_mtx_2,                  # In
                precip_fracs, silhs_config_flags,                          # In
                vert_decorr_coef,                                          # In
                err_info,                                                  # InOut
                stats,                                                     # InOut
                sampling_state,                                            # InOut
            )

            # Preserve the source fatal return between generation and clipping.
            # Outputs undefined on that path are explicit zeros in JAX.
            def clip_subcolumns(_):
                clipped, rt, thl, rc, rv, Nc = clip_transform_silhs_output_api(
                    gr.nzt, ngrdcol, parameters_microphys.lh_num_samples,  # In
                    pdf_dim, hydromet_dim, hm_metadata,                    # In
                    X_comp,                                                # In
                    X_nl,                                                  # InOut
                    pdf_params, True,                                      # In
                )
                clipped_stats = stats_accumulate_lh_api(
                    gr, gr.nzt, ngrdcol, parameters_microphys.lh_num_samples, pdf_dim,  # In
                    rho_ds_zt,                                                          # In
                    hydromet_dim, hm_metadata,                                          # In
                    weights, clipped,                                                   # In
                    rt, thl,                                                            # In
                    rc, rv,                                                             # In
                    Nc,                                                                 # In
                    sample_stats,                                                       # InOut
                )
                return clipped_stats, clipped, rt, thl, rc, rv, Nc

            def failed_sampling(_):
                zero = jnp.zeros_like(weights)
                return sample_stats, X_nl, zero, zero, zero, zero, zero

            sample_stats, X_nl, rt, thl, rc, rv, Nc = jax.lax.cond(
                sample_err.any_fatal() & clubb_at_least_debug_level(0),
                failed_sampling,
                clip_subcolumns,
                operand=None,
            )
            return sample_stats, sample_err, X_nl, X_comp, weights, rt, thl, rc, rv, Nc, next_state

        def failed_setup(_):
            zero = jnp.zeros((ngrdcol, parameters_microphys.lh_num_samples, gr.nzt))
            return (
                stats,
                err_info,
                jnp.zeros((*zero.shape, pdf_dim)),
                jnp.zeros(zero.shape, dtype=jnp.int32),
                zero,
                zero,
                zero,
                zero,
                zero,
                zero,
                sampling_state,
            )

        (
            stats,
            err_info,
            X_nl_all_levs,
            X_mixt_comp_all_levs,
            lh_sample_point_weights,
            lh_rt_clipped,
            lh_thl_clipped,
            lh_rc_clipped,
            lh_rv_clipped,
            lh_Nc_clipped,
            sampling_state,
        ) = jax.lax.cond(
            err_info.any_fatal() & clubb_at_least_debug_level(0),
            failed_setup,
            compute_subcolumns,
            operand=None,
        )
    else:
        X_nl_all_levs = jnp.zeros((ngrdcol, 0, gr.nzt, pdf_dim))
        X_mixt_comp_all_levs = jnp.zeros((ngrdcol, 0, gr.nzt), dtype=jnp.int32)
        lh_sample_point_weights = jnp.zeros((ngrdcol, 0, gr.nzt))
        lh_rt_clipped = lh_thl_clipped = lh_rc_clipped = lh_rv_clipped = lh_Nc_clipped = jnp.zeros(
            (ngrdcol, 0, gr.nzt)
        )
    return (
        stats,
        err_info,
        hydrometp2,
        mu_x_1_n,
        mu_x_2_n,
        sigma_x_1_n,
        sigma_x_2_n,
        corr_array_1_n,
        corr_array_2_n,
        corr_cholesky_mtx_1,
        corr_cholesky_mtx_2,
        precip_fracs,
        rtphmp_zt,
        thlphmp_zt,
        wp2hmp,
        X_nl_all_levs,
        X_mixt_comp_all_levs,
        lh_sample_point_weights,
        lh_rt_clipped,
        lh_thl_clipped,
        lh_rc_clipped,
        lh_rv_clipped,
        lh_Nc_clipped,
        hydromet_pdf_params,
        sampling_state,
    )
