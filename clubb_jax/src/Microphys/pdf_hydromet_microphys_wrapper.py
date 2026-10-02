"""Hydrometeor PDF preparation from pdf_hydromet_microphys_wrapper.F90.

JAX adaptation: outputs and stats are functional returns, columns are batched.
Disabled sampling outputs have zero extent; SILHS cannot be enabled here.
"""
import jax
import jax.numpy as jnp
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from clubb_jax.src.CLUBB_core.setup_clubb_pdf_params import setup_pdf_parameters_api
from clubb_jax.src.Microphys.mixed_moment_PDF_integrals import hydrometeor_mixed_moments
from clubb_jax.src.Microphys import parameters_microphys
from clubb_jax.src.CLUBB_core.hydromet_pdf_parameter_module import init_hydromet_pdf_params


def pdf_hydromet_microphys_prep(gr, ngrdcol, pdf_dim, hydromet_dim,
                               itime, vert_decorr_coef,
                               Nc_in_cloud, cloud_frac, ice_supersat_frac,
                               rho_ds_zt, Lscale, Kh_zm, hydromet, wphydrometp,
                               corr_array_n_cloud, corr_array_n_below,
                               hm_metadata, pdf_params, clubb_params,
                               clubb_config_flags, silhs_config_flags,
                               l_rad_itime, stats, err_info, precip_fracs):
    # Setup the PDF parameters.
    hydromet_pdf_params = init_hydromet_pdf_params()
    # Pure outputs are defined even on the source's inactive/early-return paths.
    hydrometp2 = jnp.zeros((ngrdcol, gr.nzm, hydromet_dim))
    mu_x_1_n = mu_x_2_n = sigma_x_1_n = sigma_x_2_n = jnp.zeros((ngrdcol, gr.nzt, pdf_dim))
    corr_array_1_n = corr_array_2_n = corr_cholesky_mtx_1 = corr_cholesky_mtx_2 = jnp.zeros((ngrdcol, gr.nzt, pdf_dim, pdf_dim))
    rtphmp_zt = thlphmp_zt = wp2hmp = jnp.zeros_like(hydromet)
    if parameters_microphys.microphys_scheme != 'none':
        # The source standalone passes its (ngrdcol,nparams) array to this
        # wrapper's rank-one (nparams) dummy. Preserve Fortran sequence association:
        # the PDF uses the first nparams elements in column-major storage order.
        # A coordinated source/API change is needed before using per-column values.
        clubb_params = jnp.broadcast_to(
            clubb_params.T.reshape(-1)[:clubb_params.shape[1]], clubb_params.shape)
        (err_info, hydrometp2, mu_x_1_n, mu_x_2_n, sigma_x_1_n, sigma_x_2_n,
         corr_array_1_n, corr_array_2_n, corr_cholesky_mtx_1, corr_cholesky_mtx_2,
         precip_fracs, hydromet_pdf_params, stats) = setup_pdf_parameters_api(
            gr, gr.nzm, gr.nzt, ngrdcol, pdf_dim, hydromet_dim,
            Nc_in_cloud, cloud_frac, Kh_zm, ice_supersat_frac, hydromet, wphydrometp,
            corr_array_n_cloud, corr_array_n_below, hm_metadata, pdf_params, clubb_params,
            clubb_config_flags.iiPDF_type, clubb_config_flags.l_use_precip_frac,
            clubb_config_flags.l_diagnose_correlations, clubb_config_flags.l_calc_w_corr,
            clubb_config_flags.l_const_Nc_in_cloud, clubb_config_flags.l_fix_w_chi_eta_correlations,
            err_info, precip_fracs, hydromet_pdf_params, stats)
        # JAX adaptation: the source returns on fatal PDF setup before mixed
        # moments and their statistics. Undefined pure outputs are returned as zero.
        rtphmp_zt, thlphmp_zt, wp2hmp, stats = jax.lax.cond(
            err_info.any_fatal() & clubb_at_least_debug_level(0),
            lambda _: (jnp.zeros_like(hydromet), jnp.zeros_like(hydromet),
                       jnp.zeros_like(hydromet), stats),
            lambda _: hydrometeor_mixed_moments(gr, gr.ngrdcol, gr.nzt, pdf_dim, hydromet_dim,
                hydromet, hm_metadata, mu_x_1_n, mu_x_2_n, sigma_x_1_n,
                sigma_x_2_n, corr_array_1_n, corr_array_2_n, pdf_params, hydromet_pdf_params,
                precip_fracs, stats), operand=None)
    # Compute subcolumns if enabled.
    if parameters_microphys.lh_microphys_type != parameters_microphys.lh_microphys_disabled:
        raise NotImplementedError('SILHS sampling is disabled in the JAX standalone')
    X_nl_all_levs = jnp.zeros((ngrdcol, 0, gr.nzt, pdf_dim))
    X_mixt_comp_all_levs = jnp.zeros((ngrdcol, 0, gr.nzt), dtype=jnp.int32)
    lh_sample_point_weights = jnp.zeros((ngrdcol, 0, gr.nzt))
    lh_rt_clipped = lh_thl_clipped = lh_rc_clipped = lh_rv_clipped = lh_Nc_clipped = jnp.zeros((ngrdcol, 0, gr.nzt))
    return (stats, err_info, hydrometp2, mu_x_1_n, mu_x_2_n,
            sigma_x_1_n, sigma_x_2_n, corr_array_1_n, corr_array_2_n,
            corr_cholesky_mtx_1, corr_cholesky_mtx_2, precip_fracs,
            rtphmp_zt, thlphmp_zt, wp2hmp, X_nl_all_levs, X_mixt_comp_all_levs,
            lh_sample_point_weights, lh_rt_clipped, lh_thl_clipped,
            lh_rc_clipped, lh_rv_clipped, lh_Nc_clipped, hydromet_pdf_params)
