"""User-facing CLUBB API via F2PY.

Functions in this module are what users call directly. Push/pull of
Fortran-derived-type data is handled transparently by
clubb_python.derived_types converters.
"""
from clubb_python.derived_types import (
    Grid,
    ConfigFlags,
    CONFIG_FLAG_NAMES,
    SclrIdx,
    NuVertResDep,
    ErrInfo,
    HmMetadata,
    pdf_parameter,
    implicit_coefs_terms,
)
from clubb_python.CLUBB_core.advance_clubb_core import advance_clubb_core
from clubb_python.CLUBB_core.advance_helper import (
    calculate_thlp2_rad,
    calc_ri_zm,
    vertical_avg,
    vertical_integral,
    pvertinterp,
    smooth_heaviside_peskin,
    smooth_max,
    smooth_min,
    calc_xpwp,
    lscale_width_vert_avg,
    wp2_term_splat_lhs,
    wp3_term_splat_lhs,
    calc_brunt_vaisala_freq_sqd,
    compute_cx_fnc_richardson,
    calc_stability_correction,
)
from clubb_python.CLUBB_core.advance_windm_edsclrm import advance_windm_edsclrm
from clubb_python.CLUBB_core.advance_wp2_wp3 import advance_wp2_wp3
from clubb_python.CLUBB_core.advance_xm_wpxp import advance_xm_wpxp
from clubb_python.CLUBB_core.advance_xp2_xpyp import advance_xp2_xpyp, update_xp2_mc
from clubb_python.CLUBB_core.advance_xp3 import advance_xp3
from clubb_python.CLUBB_core.calendar import (
    julian2gregorian_date,
    compute_current_date,
    gregorian2julian_day,
    leap_year,
)
from clubb_python.CLUBB_core.calc_roots import cubic_solve, quadratic_solve
from clubb_python.CLUBB_core.calc_pressure import init_pressure, calculate_thvm, hydrostatic
from clubb_python.CLUBB_core.clip_explicit import (
    clip_rcm,
    clip_covar,
    clip_variance,
    clip_covars_denom,
    clip_skewness,
)
from clubb_python.CLUBB_core.diffusion import diffusion_zt_lhs, diffusion_zm_lhs
from clubb_python.CLUBB_core.fill_holes import (
    fill_holes_vertical,
    fill_holes_wp2_from_horz_tke,
)
from clubb_python.prescribe_forcings import prescribe_forcings
from clubb_python.CLUBB_core.grid_class import (
    setup_grid,
    zm2zt2zm,
    zt2zm2zt,
    zt2zm,
    zm2zt,
    ddzm,
    ddzt,
    cleanup_grid,
    setup_grid_heights,
    read_grid_heights,
)
from clubb_python.CLUBB_core.model_flags import (
    get_default_config_flags,
    init_config_flags,
)
from clubb_python.CLUBB_core.array_index import (
    set_sclr_idx,
    get_sclr_idx,
)
from clubb_python.CLUBB_core.err_info_type_module import (
    init_err_info,
    get_err_code,
    cleanup_err_info,
    set_err_info_values,
)
from clubb_python.CLUBB_core.parameters_tunable import (
    init_clubb_params,
    get_param_names,
    calc_derrived_params,
    check_parameters,
)
from clubb_python.CLUBB_core.pdf_parameter_module import (
    init_pdf_params,
    init_pdf_params_zm,
    init_pdf_implicit,
    init_pdf_implicit_coefs_terms,
    zero_pdf_params,
    zero_pdf_implicit_coefs_terms,
)
from clubb_python.CLUBB_core.error_code import set_debug_level, reset_err_code, initialize_error_headers
from clubb_python.CLUBB_core.file_functions import file_read_1d
from clubb_python.CLUBB_core.interpolation import lin_interpolate_two_points, mono_cubic_interp, zlinterp_fnc
from clubb_python.CLUBB_core.index_mapping import (
    hydromet2pdf_idx,
    mvr_hm_max,
    nx2rx_hm_idx,
    pdf2hydromet_idx,
    rx2nx_hm_idx,
)
from clubb_python.CLUBB_core.lapack_interfaces import (
    lapack_gbsv,
    lapack_gbsvx,
    lapack_gtsv,
    lapack_gtsvx,
    lapack_laqsy,
    lapack_poequ,
    lapack_potrf,
    lapack_syev,
)
from clubb_python.CLUBB_core.lapack_wrap import (
    lapack_band_solve,
    lapack_band_solvex,
    lapack_tridiag_solve,
    lapack_tridiag_solvex,
)
from clubb_python.CLUBB_core.mean_adv import term_ma_zt_lhs, term_ma_zm_lhs
from clubb_python.CLUBB_core.matrix_operations import cholesky_factor, mirror_lower_triangular_matrix
from clubb_python.CLUBB_core.matrix_solver_wrapper import band_solve, tridiag_solve
from clubb_python.CLUBB_core.mixing_length import set_lscale_max, diagnose_lscale_from_tau, calc_lscale_directly
from clubb_python.CLUBB_core.adg1_adg2_3d_luhar_pdf import (
    calc_luhar_params,
    close_luhar_pdf,
    adg1_w_closure,
    adg2_pdf_driver,
    adg1_pdf_driver,
    luhar_3d_pdf_driver,
)
from clubb_python.CLUBB_core.ly93_pdf import calc_params_ly93, ly93_driver
from clubb_python.CLUBB_core.mono_flux_limiter import calc_turb_adv_range, monotonic_turbulent_flux_limit
from clubb_python.CLUBB_core.nc_ncn_eqns import nc_in_cloud_to_ncnm
from clubb_python.CLUBB_core.precipitation_fraction import precip_fraction
from clubb_python.CLUBB_core.numerical_check import (
    check_clubb_settings,
    calculate_spurious_source,
    sfc_varnce_check,
    parameterization_check,
    pdf_closure_check,
    length_check,
)
from clubb_python.CLUBB_core.new_hybrid_pdf import (
    calculate_w_params,
    calculate_responder_params,
    calculate_coef_wp4_implicit,
    calc_coef_wp2xp_implicit,
    calc_coefs_wpxp2_semiimpl,
)
from clubb_python.CLUBB_core.new_hybrid_pdf_main import new_hybrid_pdf_driver
from clubb_python.CLUBB_core.new_pdf import (
    calc_setter_var_params,
    calc_responder_params,
    calc_limits_f_x_responder,
    calc_coef_wp4_implicit,
    calc_coef_wpxp2_implicit,
    calc_coefs_wp2xp_semiimpl,
    calc_coefs_wpxpyp_semiimpl,
)
from clubb_python.CLUBB_core.new_pdf_main import new_pdf_driver
from clubb_python.CLUBB_core.new_tsdadg_pdf import (
    calc_l_x_skx_fnc,
    calc_setter_parameters_tsdadg,
    tsdadg_pdf_driver,
)
from clubb_python.CLUBB_core.corr_varnce_module import assert_corr_symmetric
from clubb_python.CLUBB_core.diagnose_correlations_module import diagnose_correlations, calc_cholesky_corr_mtx_approx
from clubb_python.CLUBB_core.pdf_utilities import (
    smooth_corr_quotient,
    calc_comp_corrs_binormal,
    mean_l2n,
    stdev_l2n,
    corr_nn2nl,
    corr_nn2ll,
    compute_mean_binormal,
    compute_variance_binormal,
    calc_corr_chi_x,
    calc_corr_eta_x,
)
from clubb_python.CLUBB_core.pos_definite_module import pos_definite_adj
from clubb_python.CLUBB_core.penta_bicgstab_solver import penta_bicgstab_solve
from clubb_python.CLUBB_core.penta_lu_solver import penta_lu_solve
from clubb_python.CLUBB_core.remapping_module import remap_vals_to_target
from clubb_python.CLUBB_core.skx_module import skx_func, lg_2005_ansatz, xp3_lg_2005_ansatz
from clubb_python.CLUBB_core.sigma_sqd_w_module import compute_sigma_sqd_w
from clubb_python.CLUBB_core.sponge_layer_damping import sponge_damp_xm, sponge_damp_xp2, sponge_damp_xp3
from clubb_python.CLUBB_core.sfc_varnce_module import calc_sfc_varnce
from clubb_python.CLUBB_core.stats_clubb_utilities import stats_accumulate
from clubb_python.CLUBB_core.stats_netcdf import (
    init_stats,
    finalize_stats,
    get_stats_config,
    get_stats_var_meta,
    get_stats_var_data,
    set_stats_var_data,
    stats_begin_timestep,
    stats_end_timestep,
    stats_update,
    stats_begin_budget,
    stats_update_budget,
    stats_finalize_budget,
    stats_update_scalar,
    stats_update_1d,
    stats_update_2d,
    stats_begin_budget_scalar,
    stats_begin_budget_1d,
    stats_begin_budget_2d,
    stats_update_budget_scalar,
    stats_update_budget_1d,
    stats_update_budget_2d,
    stats_finalize_budget_scalar,
    stats_finalize_budget_1d,
    stats_finalize_budget_2d,
    var_on_stats_list,
    stats_update_grid,
    stats_lh_samples_init,
    stats_lh_samples_write_lognormal,
    stats_lh_samples_write_uniform,
)
from clubb_python.CLUBB_core.t_in_k_module import thlm2t_in_k
from clubb_python.CLUBB_core.tridiag_lu_solver import tridiag_lu_solve
from clubb_python.CLUBB_core.turbulent_adv_pdf import (
    xpyp_term_ta_pdf_lhs,
    xpyp_term_ta_pdf_lhs_godunov,
    xpyp_term_ta_pdf_rhs,
    xpyp_term_ta_pdf_rhs_godunov,
)
from clubb_python.CLUBB_core.saturation import sat_mixrat_liq, sat_mixrat_ice, rcm_sat_adj
from clubb_python.CLUBB_core.pdf_closure import (
    pdf_closure_driver,
    calc_wp4_pdf,
    calc_wp2xp_pdf,
    calc_wpxp2_pdf,
    calc_wpxpyp_pdf,
    calc_w_up_in_cloud,
)
from clubb_python.radiation import cos_solar_zen, sunray_sw, set_simplified_radiation_params
