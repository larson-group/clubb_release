"""Microphysics lifecycle from microphys_init_cleanup.F90.

Host adaptation: ``namelist_file`` may be the driver's already parsed namelist
mapping. Fortran module parameters retain module ownership. File-unit/report
arguments are retained, while the standalone driver owns output-file lifetime.
"""
import sys
from pathlib import Path
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from dataclasses import fields
import importlib
import jax
from clubb_jax.src.Microphys import parameters_microphys as parameters
from clubb_jax.src.Microphys.KK_microphys import parameters_KK
from clubb_jax.src.Microphys.Morrison_microphys import module_mp_graupel
from clubb_jax.src.Input_fields.namelist import read_namelist
from clubb_jax.src.Input_fields.corr_varnce_input_reader import read_correlation_matrix
from clubb_jax.src.CLUBB_core.corr_varnce_module import (
    init_pdf_hydromet_arrays_api, setup_corr_varnce_array_api,
    hmp2_ip_on_hmm2_ip_slope_type, hmp2_ip_on_hmm2_ip_intrcpt_type,
)


def init_microphys(iunit, runtype, namelist_file, case_info_file,
                   host_dx, host_dy, clubb_params, l_diagnose_correlations,
                   l_const_Nc_in_cloud, l_fix_w_chi_eta_correlations):
    # Set default values, then read in the namelist.
    # Reset module state between independently initialized cases.
    # Namelist values are module-owned constants captured at trace time. A new
    # case must not reuse executables traced with the previous case's physics.
    # Description:
    # Set indices to the various hydrometeor species and define hydromet_dim for
    # the purposes of allocating memory.
    # References:
    # None
    #-----------------------------------------------------------------------
    jax.clear_caches()
    importlib.reload(parameters)
    importlib.reload(parameters_KK)
    # Reset the source module defaults for independently initialized JAX cases.
    module_mp_graupel.NNUCCD_REDUCE_COEF = 1.0
    module_mp_graupel.NNUCCC_REDUCE_COEF = 1.0

    #--------------------------------------------------------------------------
    # Parameters for NNUCCD & NNUCCC coefficients on clex9_oct14 case
    #--------------------------------------------------------------------------
    if runtype.strip() == 'clex9_oct14':
        module_mp_graupel.NNUCCD_REDUCE_COEF = .01  # Reduce NNUCCD by factor of 100 for clex9_oct14
        module_mp_graupel.NNUCCC_REDUCE_COEF = .01  # Reduce NNUCCC by factor of 100 for clex9_oct14

    cfg = namelist_file if isinstance(namelist_file, dict) else read_namelist(str(namelist_file))
    cfg = {name.lower(): value for name, value in cfg.items()}
    # Only module-owned entries of /microphysics_setting/ are writable.
    # Other namelists in the driver's merged mapping must not alter enum constants.
    for name in (
        'microphys_scheme', 'l_cloud_sed', 'sigma_g', 'l_ice_microphys',
        'l_graupel', 'l_hail', 'l_var_covar_src', 'l_upwind_diff_sed',
        'l_seifert_beheng', 'l_predict_Nc', 'specify_aerosol', 'l_subgrid_w',
        'l_arctic_nucl', 'l_cloud_edge_activation', 'l_fix_pgam',
        'l_in_cloud_Nc_diff', 'lh_microphys_type', 'l_local_kk',
        'lh_num_samples', 'lh_sequence_length', 'lh_seed',
        'l_silhs_KK_convergence_adj_mean', 'microphys_start_time', 'Nc0_in_cloud',
    ):
        if name.lower() in cfg:
            setattr(parameters, name, cfg[name.lower()])
    # l_gfdl_activation belongs to /gfdl_activation_setting/.
    parameters.l_gfdl_activation = bool(cfg.get('l_gfdl_activation', False))
    # TODO(port-mirror): source namelist/correlation-matrix reports to
    # case_info_file are not yet emitted by the host output layer. Restore them
    # there when that layer accepts microphysics configuration reports.
    parameters.lh_microphys_type = {'disabled': 3, 'interactive': 1, 'non-interactive': 2}.get(
        parameters.lh_microphys_type, parameters.lh_microphys_type)
    parameters_KK.C_evap = float(cfg.get('c_evap', parameters_KK.C_evap))
    parameters_KK.r_0 = float(cfg.get('r_0', parameters_KK.r_0))
    from clubb_jax.src.Microphys.KK_microphys import parabolic_cylinder
    parabolic_cylinder.l_high_accuracy_parab_cyl_fnc = bool(cfg.get('l_high_accuracy_parab_cyl_fnc', False))
    if cfg.get('l_morr_xp2_mc', False):
        raise ValueError('l_morr_xp2_mc is not implemented in the JAX Morrison interface')
    for name in ('l_hail', 'l_seifert_beheng', 'l_arctic_nucl', 'l_cloud_edge_activation', 'l_fix_pgam'):
        if cfg.get(name, False):
            raise ValueError(f'{name} is not implemented in the JAX Morrison configuration')
    if parameters.lh_microphys_type != parameters.lh_microphys_disabled:
        raise ValueError('SILHS microphysics is not supported')
    if parameters.l_gfdl_activation:
        raise ValueError('GFDL activation core is not supported')
    if parameters.microphys_scheme not in ('none', 'khairoutdinov_kogan', 'morrison'):
        raise ValueError(f'Unsupported microphys_scheme: {parameters.microphys_scheme}')

    # Set indices to the various hydrometeor species and define hydromet_dim.
    iirr = iiNr = iiri = iiNi = iirs = iiNs = iirg = iiNg = -1
    hydromet_dim = 0
    if parameters.microphys_scheme == 'morrison':
        # GRAUPEL_INIT constructs a pure parameter mapping when the core is
        # traced, rather than mutating module arrays here. It reads the flags
        # set above; clear_caches prevents stale case constants. Aerosol ccn/aer
        # namelist constants remain unavailable; the guard below rejects their
        # active use. TODO: map these constants when activation is implemented.
        if parameters.specify_aerosol not in ('morrison_no_aerosol', 'morrison_power_law', 'morrison_lognormal'):
            raise ValueError('Unknown Morrison aerosol mode')
        if parameters.l_predict_Nc and parameters.specify_aerosol != 'morrison_no_aerosol':
            raise ValueError('Predicted Morrison Nc with aerosol activation is not implemented in the JAX configuration')
        iirr, iiNr, hydromet_dim = 0, 1, 2
        if parameters.l_ice_microphys:
            iiri, iiNi, iirs, iiNs, hydromet_dim = 2, 3, 4, 5, 6
            if parameters.l_graupel:
                iirg, iiNg, hydromet_dim = 6, 7, 8
        if parameters.l_cloud_sed:
            raise ValueError('Morrison includes cloud sedimentation; l_cloud_sed must be false')
        parameters.l_hydromet_sed = (False,) * hydromet_dim
    elif parameters.microphys_scheme == 'khairoutdinov_kogan':
        if parameters.l_predict_Nc:
            raise ValueError('Khairoutdinov-Kogan does not support l_predict_Nc')
        iirr, iiNr, hydromet_dim = 0, 1, 2
        parameters.l_hydromet_sed = (True, True)
    else:
        parameters.l_predict_Nc = False
        parameters.l_hydromet_sed = ()

    # Initialize hydrometeor metadata and prescribed in-precipitation variances.
    # f90nml gives nested mappings; the fallback parser uses percent-separated keys.
    slope = {}
    intercept = {}
    for field in fields(hmp2_ip_on_hmm2_ip_slope_type):
        name = field.name
        group = cfg.get('hmp2_ip_on_hmm2_ip_slope', {})
        if name.lower() in group:
            slope[name] = group[name.lower()]
        if f'hmp2_ip_on_hmm2_ip_slope%{name.lower()}' in cfg:
            slope[name] = cfg[f'hmp2_ip_on_hmm2_ip_slope%{name.lower()}']
        group = cfg.get('hmp2_ip_on_hmm2_ip_intrcpt', {})
        if name.lower() in group:
            intercept[name] = group[name.lower()]
        if f'hmp2_ip_on_hmm2_ip_intrcpt%{name.lower()}' in cfg:
            intercept[name] = cfg[f'hmp2_ip_on_hmm2_ip_intrcpt%{name.lower()}']
    hm_metadata, pdf_dim = init_pdf_hydromet_arrays_api(
        host_dx, host_dy, hydromet_dim, iirr, iiNr, iiri, iiNi, iirs, iiNs, iirg, iiNg,
        float(cfg.get('ncnp2_on_ncnm2', 1.0)),
        hmp2_ip_on_hmm2_ip_slope_type(**slope), hmp2_ip_on_hmm2_ip_intrcpt_type(**intercept),
    )
    corr_input_path = Path(__file__).resolve().parents[3] / 'input/case_setups'
    corr_file_path_cloud = corr_input_path / f'{runtype}_corr_array_cloud.in'
    corr_file_path_below = corr_input_path / f'{runtype}_corr_array_below.in'
    if corr_file_path_cloud.exists() and corr_file_path_below.exists():
        corr_array_n_cloud_in = read_correlation_matrix(iunit, corr_file_path_cloud, pdf_dim, hm_metadata, None)
        corr_array_n_below_in = read_correlation_matrix(iunit, corr_file_path_below, pdf_dim, hm_metadata, None)
        corr_array_n_cloud, corr_array_n_below = setup_corr_varnce_array_api(
            pdf_dim, hm_metadata, l_fix_w_chi_eta_correlations,
            corr_array_n_cloud_in, corr_array_n_below_in)
    else:
        if clubb_at_least_debug_level(1):
            print(f'Warning: missing correlation input file(s): {corr_file_path_cloud} '
                  f'and/or {corr_file_path_below}', file=sys.stderr)
            print('The default correlation arrays will be used.', file=sys.stderr)
        corr_array_n_cloud, corr_array_n_below = setup_corr_varnce_array_api(
            pdf_dim, hm_metadata, l_fix_w_chi_eta_correlations)
    # TODO(port-mirror): SILHS configuration is a disabled sentinel until its
    # driver is ported; no sampling routine may consume this value.
    silhs_config_flags = None
    vert_decorr_coef_out = float(cfg.get('vert_decorr_coef', 0.0))
    return (hydromet_dim, pdf_dim, hm_metadata, silhs_config_flags,
            vert_decorr_coef_out, corr_array_n_cloud, corr_array_n_below)


def cleanup_microphys():
    # Description:
    # De-allocate arrays used by the microphysics
    # References:
    # None
    #-----------------------------------------------------------------------
    parameters.l_hydromet_sed = ()
