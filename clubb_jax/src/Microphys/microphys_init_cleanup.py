"""Microphysics lifecycle from microphys_init_cleanup.F90.

Host adaptation: ``namelist_file`` may be the driver's already parsed namelist mapping. Fortran
module parameters retain module ownership. File-unit/report arguments are retained, while the
standalone driver owns output-file lifetime.
"""

import sys
from pathlib import Path
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from dataclasses import fields
import importlib
import jax
import jax.numpy as jnp
import numpy as np
from clubb_jax.src.SILHS import parameters_silhs
from clubb_jax.src.Microphys import parameters_microphys as parameters
from clubb_jax.src.Microphys.KK_microphys import parameters_KK
from clubb_jax.src.Microphys.Morrison_microphys import module_mp_graupel
from clubb_jax.src.Input_fields.namelist import read_namelist
from clubb_jax.src.Input_fields.corr_varnce_input_reader import read_correlation_matrix
from clubb_jax.src.CLUBB_core.corr_varnce_module import (
    init_pdf_hydromet_arrays_api,
    setup_corr_varnce_array_api,
    hmp2_ip_on_hmm2_ip_slope_type,
    hmp2_ip_on_hmm2_ip_intrcpt_type,
)


# -----------------------------------------------------------------------------
def init_microphys(
    iunit, runtype, namelist_file, case_info_file,  # In
    host_dx, host_dy,                               # In
    clubb_params,                                   # In
    l_diagnose_correlations,                        # In
    l_const_Nc_in_cloud,                            # InOut
    l_fix_w_chi_eta_correlations,                   # InOut
):
    """Initialize scheme metadata, correlation arrays and SILHS configuration.

    Host namelist input may be an already parsed mapping or a file path. Configuration remains
    module owned, as in Fortran, and is captured at trace time. Return dimensions, metadata, flags
    and correlation arrays.

    Arguments:
        iunit: Fortran unit-number input retained for interface correspondence; Python opens
            the path.
        runtype: Benchmark case name; selects prescribed correlation input files.
        namelist_file: File name
        case_info_file: Existing simulation info file (plain text); None omits file output
        host_dx: Host horizontal grid spacing in x [m].
        host_dy: Host horizontal grid spacing in y [m].
        clubb_params: Column-dependent tunable CLUBB parameters.
        l_diagnose_correlations: Diagnose correlations instead of using fixed ones
        l_const_Nc_in_cloud: Use a constant cloud droplet conc. within cloud (K&K)
        l_fix_w_chi_eta_correlations: Use a fixed correlation for s and t Mellor(chi/eta)
    """
    # Set default values, then read in the namelist.
    # Reset module state between independently initialized cases.
    # Namelist values are module-owned constants captured at trace time. A new
    # case must not reuse executables traced with the previous case's physics.
    # Description:
    # Set indices to the various hydrometeor species and define hydromet_dim for
    # the purposes of allocating memory.
    # References:
    # None
    # -----------------------------------------------------------------------
    jax.clear_caches()
    importlib.reload(parameters)
    importlib.reload(parameters_KK)
    importlib.reload(parameters_silhs)

    # Reset the source module defaults for independently initialized JAX cases.
    module_mp_graupel.NNUCCD_REDUCE_COEF = 1.0
    module_mp_graupel.NNUCCC_REDUCE_COEF = 1.0

    # --------------------------------------------------------------------------
    # Parameters for NNUCCD & NNUCCC coefficients on clex9_oct14 case
    # --------------------------------------------------------------------------
    if runtype.strip() == "clex9_oct14":
        module_mp_graupel.NNUCCD_REDUCE_COEF = (
            0.01  # Reduce NNUCCD by factor of 100 for clex9_oct14
        )
        module_mp_graupel.NNUCCC_REDUCE_COEF = (
            0.01  # Reduce NNUCCC by factor of 100 for clex9_oct14
        )

    cfg = namelist_file if isinstance(namelist_file, dict) else read_namelist(str(namelist_file))
    cfg = {name.lower(): value for name, value in cfg.items()}

    # Only module-owned entries of /microphysics_setting/ are writable.
    # Other namelists in the driver's merged mapping must not alter enum constants.
    for name in (
        "microphys_scheme",
        "l_cloud_sed",
        "sigma_g",
        "l_ice_microphys",
        "l_graupel",
        "l_hail",
        "l_var_covar_src",
        "l_upwind_diff_sed",
        "l_seifert_beheng",
        "l_predict_Nc",
        "specify_aerosol",
        "l_subgrid_w",
        "l_arctic_nucl",
        "l_cloud_edge_activation",
        "l_fix_pgam",
        "l_in_cloud_Nc_diff",
        "lh_microphys_type",
        "l_local_kk",
        "lh_num_samples",
        "lh_sequence_length",
        "lh_seed",
        "l_silhs_KK_convergence_adj_mean",
        "microphys_start_time",
        "Nc0_in_cloud",
    ):
        if name.lower() in cfg:
            setattr(parameters, name, cfg[name.lower()])
    # l_gfdl_activation belongs to /gfdl_activation_setting/.
    parameters.l_gfdl_activation = bool(cfg.get("l_gfdl_activation", False))

    parameters.lh_microphys_type = {"disabled": 3, "interactive": 1, "non-interactive": 2}.get(
        parameters.lh_microphys_type, parameters.lh_microphys_type
    )
    parameters_KK.C_evap = float(cfg.get("c_evap", parameters_KK.C_evap))
    parameters_KK.r_0 = float(cfg.get("r_0", parameters_KK.r_0))
    from clubb_jax.src.Microphys.KK_microphys import parabolic_cylinder

    parabolic_cylinder.l_high_accuracy_parab_cyl_fnc = bool(
        cfg.get("l_high_accuracy_parab_cyl_fnc", False)
    )
    if cfg.get("l_morr_xp2_mc", False):
        raise ValueError("l_morr_xp2_mc is not implemented in the JAX Morrison interface")
    for name in (
        "l_hail",
        "l_seifert_beheng",
        "l_arctic_nucl",
        "l_cloud_edge_activation",
        "l_fix_pgam",
    ):
        if cfg.get(name, False):
            raise ValueError(f"{name} is not implemented in the JAX Morrison configuration")
    if parameters.lh_microphys_type not in (1, 2, 3):
        raise ValueError("Unknown lh_microphys_type")

    # Validate sampled-microphysics configuration before allocation/timestepping.
    if parameters.lh_microphys_type != parameters.lh_microphys_disabled:
        if parameters.microphys_scheme not in ("khairoutdinov_kogan", "morrison"):
            raise ValueError("SILHS requires KK or Morrison microphysics")
        if (
            parameters.microphys_scheme == "khairoutdinov_kogan"
            and not parameters.l_local_kk
            and parameters.lh_microphys_type == parameters.lh_microphys_interactive
        ):
            raise ValueError("Interactive SILHS requires l_local_kk = true for KK microphysics")
        if parameters.lh_num_samples < 1 or parameters.lh_sequence_length < 1:
            raise ValueError("SILHS sample count and sequence length must be positive")
    if (
        parameters.l_silhs_KK_convergence_adj_mean
        and parameters.microphys_scheme != "khairoutdinov_kogan"
    ):
        raise ValueError("l_silhs_KK_convergence_adj_mean requires KK microphysics")
    if parameters.l_gfdl_activation:
        raise ValueError("GFDL activation core is not supported")
    if parameters.microphys_scheme not in ("none", "khairoutdinov_kogan", "morrison"):
        raise ValueError(f"Unsupported microphys_scheme: {parameters.microphys_scheme}")

    # Initialize hydrometeor metadata and prescribed in-precipitation variances.
    # f90nml gives nested mappings; the fallback parser uses percent-separated keys.
    slope = {}
    intercept = {}
    for field in fields(hmp2_ip_on_hmm2_ip_slope_type):
        name = field.name
        group = cfg.get("hmp2_ip_on_hmm2_ip_slope", {})
        if name.lower() in group:
            slope[name] = group[name.lower()]
        if f"hmp2_ip_on_hmm2_ip_slope%{name.lower()}" in cfg:
            slope[name] = cfg[f"hmp2_ip_on_hmm2_ip_slope%{name.lower()}"]
        group = cfg.get("hmp2_ip_on_hmm2_ip_intrcpt", {})
        if name.lower() in group:
            intercept[name] = group[name.lower()]
        if f"hmp2_ip_on_hmm2_ip_intrcpt%{name.lower()}" in cfg:
            intercept[name] = cfg[f"hmp2_ip_on_hmm2_ip_intrcpt%{name.lower()}"]

    # Initialize source-owned SILHS configuration before the timestep path.
    silhs_config_flags = parameters_silhs.silhs_config_flags_type(
        **{
            f.name: cfg[f.name.lower()]
            for f in fields(parameters_silhs.silhs_config_flags_type)
            if f.name.lower() in cfg
        }
    )
    parameters_silhs.importance_prob_thresh = float(cfg.get("importance_prob_thresh", 1.0e-8))
    parameters_silhs.vert_decorr_coef = float(cfg.get("vert_decorr_coef", 0.03))

    # f90nml retains derived types as mappings; the fallback reader keeps '%'
    # names flat. Both represent the source namelist derived-type fields.
    prescribed_probs = cfg.get("eight_cluster_presc_probs", {})
    parameters_silhs.eight_cluster_presc_probs = parameters_silhs.eight_cluster_presc_probs_type(
        **{
            f.name: cfg.get(
                "eight_cluster_presc_probs%" + f.name,
                prescribed_probs.get(
                    f.name, getattr(parameters_silhs.eight_cluster_presc_probs, f.name)
                ),
            )
            for f in fields(parameters_silhs.eight_cluster_presc_probs_type)
        }
    )

    # Write the source microphysics namelist report to the screen and append
    # it to the standalone-owned setup file. Python streams replace file units
    # and write_text overloads; field order follows the Fortran call block.
    if clubb_at_least_debug_level(1):
        report = [
            "--------------------------------------------------",
            "&microphysics_setting",
            "--------------------------------------------------",
        ]
        for name in (
            "microphys_scheme", "l_cloud_sed", "sigma_g", "l_graupel",
            "l_hail", "l_seifert_beheng", "l_predict_Nc",
        ):
            report.append(f"{name} = {getattr(parameters, name)}")
        report.append(f"l_const_Nc_in_cloud = {l_const_Nc_in_cloud}")
        for name in (
            "specify_aerosol", "l_subgrid_w", "l_arctic_nucl",
            "l_cloud_edge_activation", "l_fix_pgam", "l_in_cloud_Nc_diff",
            "l_var_covar_src", "l_upwind_diff_sed",
        ):
            report.append(f"{name} = {getattr(parameters, name)}")
        lh_microphys_type = {
            parameters.lh_microphys_interactive: "interactive",
            parameters.lh_microphys_non_interactive: "non-interactive",
            parameters.lh_microphys_disabled: "disabled",
        }[parameters.lh_microphys_type]
        report.append(f"lh_microphys_type = {lh_microphys_type}")
        for name in ("lh_num_samples", "lh_sequence_length", "lh_seed"):
            report.append(f"{name} = {getattr(parameters, name)}")
        report.extend([
            f"l_fix_w_chi_eta_correlations = {l_fix_w_chi_eta_correlations}",
            f"l_silhs_KK_convergence_adj_mean = {parameters.l_silhs_KK_convergence_adj_mean}",
            f"importance_prob_thresh = {parameters_silhs.importance_prob_thresh}",
            f"host_dx = {host_dx}",
            f"host_dy = {host_dy}",
        ])
        for name, values in (
            ("hmp2_ip_on_hmm2_ip_slope", hmp2_ip_on_hmm2_ip_slope_type(**slope)),
            ("hmp2_ip_on_hmm2_ip_intrcpt", hmp2_ip_on_hmm2_ip_intrcpt_type(**intercept)),
        ):
            # Retain the source report's repeated Ni entry.
            for species in ("rr", "ri", "rs", "rg", "Nr", "Ni", "Ni", "Ng"):
                report.append(f"{name}%{species} = {getattr(values, species)}")
        report.extend([
            f"Ncnp2_on_Ncnm2 = {float(cfg.get('ncnp2_on_ncnm2', 1.0))}",
            f"C_evap = {parameters_KK.C_evap}",
            f"r_0 = {parameters_KK.r_0}",
            f"microphys_start_time = {parameters.microphys_start_time}",
            f"Nc0_in_cloud = {parameters.Nc0_in_cloud}",
        ])
        # Morrison's source aerosol namelist values have single precision.
        # Their active use retains the initialization support gates below.
        for name, default in (
            ("ccnconst", 120.0), ("ccnexpnt", 0.4),
            ("aer_rm1", 0.011e-6), ("aer_rm2", 0.06e-6),
            ("aer_n1", 125.0e6), ("aer_n2", 65.0e6),
            ("aer_sig1", 1.2), ("aer_sig2", 1.7), ("pgam_fixed", 5.0),
        ):
            report.append(f"{name} = {float(np.float32(cfg.get(name, default)))}")
        from clubb_jax.src.CLUBB_core.precipitation_fraction import precip_frac_calc_type

        report.extend([
            f"precip_frac_calc_type = {precip_frac_calc_type}",
            "--------------------------------------------------",
            "&SILHS_setting",
            "--------------------------------------------------",
        ])
        for line in report:
            print(line)
        if case_info_file is not None:
            # The standalone creates the file; lifecycle initialization only
            # appends its configuration, as in the source status='old' open.
            with Path(case_info_file).open("r+") as report_file:
                report_file.seek(0, 2)
                for line in report:
                    print(line, file=report_file)
                parameters_silhs.print_silhs_config_flags_api(
                    report_file, silhs_config_flags,  # In
                )

    # Set indices to the various hydrometeor species and define hydromet_dim.
    iirr = iiNr = iiri = iiNi = iirs = iiNs = iirg = iiNg = -1
    hydromet_dim = 0
    if parameters.microphys_scheme == "morrison":
        # GRAUPEL_INIT constructs a pure parameter mapping when the core is
        # traced, rather than mutating module arrays here. It reads the flags
        # set above; clear_caches prevents stale case constants. Aerosol ccn/aer
        # namelist constants remain unavailable; the guard below rejects their
        # active use. TODO: map these constants when activation is implemented.
        if parameters.specify_aerosol not in (
            "morrison_no_aerosol",
            "morrison_power_law",
            "morrison_lognormal",
        ):
            raise ValueError("Unknown Morrison aerosol mode")
        if parameters.l_predict_Nc and parameters.specify_aerosol != "morrison_no_aerosol":
            raise ValueError(
                "Predicted Morrison Nc with aerosol activation is not implemented "
                "in the JAX configuration"
            )
        iirr, iiNr, hydromet_dim = 0, 1, 2
        if parameters.l_ice_microphys:
            iiri, iiNi, iirs, iiNs, hydromet_dim = 2, 3, 4, 5, 6
            if parameters.l_graupel:
                iirg, iiNg, hydromet_dim = 6, 7, 8
        if parameters.l_cloud_sed:
            raise ValueError("Morrison includes cloud sedimentation; l_cloud_sed must be false")
        parameters.l_hydromet_sed = (False,) * hydromet_dim
    elif parameters.microphys_scheme == "khairoutdinov_kogan":
        if parameters.l_predict_Nc:
            raise ValueError("Khairoutdinov-Kogan does not support l_predict_Nc")
        iirr, iiNr, hydromet_dim = 0, 1, 2
        parameters.l_hydromet_sed = (True, True)
    else:
        parameters.l_predict_Nc = False
        parameters.l_hydromet_sed = ()

    hm_metadata, pdf_dim = init_pdf_hydromet_arrays_api(
        host_dx, host_dy, hydromet_dim,                # In
        iirr, iiNr, iiri, iiNi,                        # In
        iirs, iiNs, iirg, iiNg,                        # In
        float(cfg.get("ncnp2_on_ncnm2", 1.0)),         # In
        hmp2_ip_on_hmm2_ip_slope_type(**slope),        # In
        hmp2_ip_on_hmm2_ip_intrcpt_type(**intercept),  # In
    )
    corr_input_path = Path(__file__).resolve().parents[3] / "input/case_setups"
    corr_file_path_cloud = corr_input_path / f"{runtype}_corr_array_cloud.in"
    corr_file_path_below = corr_input_path / f"{runtype}_corr_array_below.in"
    if corr_file_path_cloud.exists() and corr_file_path_below.exists():
        corr_array_n_cloud_in = read_correlation_matrix(
            iunit, corr_file_path_cloud,  # In
            pdf_dim, hm_metadata,         # In
            None,                         # InOut
        )
        corr_array_n_below_in = read_correlation_matrix(
            iunit, corr_file_path_below,  # In
            pdf_dim, hm_metadata,         # In
            None,                         # InOut
        )
        corr_array_n_cloud, corr_array_n_below = setup_corr_varnce_array_api(
            pdf_dim, hm_metadata,          # In
            l_fix_w_chi_eta_correlations,  # In
            corr_array_n_cloud_in,         # In
            corr_array_n_below_in,         # In
        )
    else:
        if clubb_at_least_debug_level(1):
            print(
                f"Warning: missing correlation input file(s): {corr_file_path_cloud} "
                f"and/or {corr_file_path_below}",
                file=sys.stderr,
            )
            print("The default correlation arrays will be used.", file=sys.stderr)
        corr_array_n_cloud, corr_array_n_below = setup_corr_varnce_array_api(
            pdf_dim, hm_metadata, l_fix_w_chi_eta_correlations
        )

    # These approximate physical-space correlations are a guide: diagnosed
    # variances can change them during timestepping. The source prints them
    # only for fixed variance ratios (zeta_vrnce_rat = 0) and active microphysics.
    from clubb_jax.src.CLUBB_core.parameter_indices import iomicron, izeta_vrnce_rat

    if (
        clubb_at_least_debug_level(1)
        and parameters.microphys_scheme != "none"
        and abs(float(clubb_params[0, izeta_vrnce_rat])) < np.finfo(float).eps
    ):
        from clubb_jax.src.CLUBB_core.index_mapping import pdf2hydromet_idx
        from clubb_jax.src.CLUBB_core.pdf_utilities import stdev_L2N
        from clubb_jax.src.CLUBB_core.setup_clubb_pdf_params import denorm_transform_corr
        from clubb_jax.src.CLUBB_core.matrix_operations import mirror_lower_triangular_matrix

        # The standalone source uses one column and one level for this report.
        sigma2_on_mu2_ip_cloud = jnp.zeros((1, 1, pdf_dim))
        if not l_const_Nc_in_cloud:
            sigma2_on_mu2_ip_cloud = sigma2_on_mu2_ip_cloud.at[..., hm_metadata.iiPDF_Ncn].set(
                hm_metadata.Ncnp2_on_Ncnm2
            )
        for ivar in range(hm_metadata.iiPDF_Ncn + 1, pdf_dim):
            sigma2_on_mu2_ip_cloud = sigma2_on_mu2_ip_cloud.at[..., ivar].set(
                clubb_params[0, iomicron]
                * hm_metadata.hmp2_ip_on_hmm2_ip[pdf2hydromet_idx(ivar, hm_metadata)]
            )
        sigma2_on_mu2_ip_below = sigma2_on_mu2_ip_cloud
        sigma_x_n_cloud = stdev_L2N(sigma2_on_mu2_ip_cloud)
        sigma_x_n_below = sigma_x_n_cloud
        corr_array_cloud, corr_array_below = denorm_transform_corr(
            sigma_x_n_cloud, sigma_x_n_below,                 # In
            sigma2_on_mu2_ip_cloud, sigma2_on_mu2_ip_below,  # In
            corr_array_n_cloud[None, None, :, :],             # In
            corr_array_n_below[None, None, :, :],             # In
            hm_metadata.iiPDF_chi, hm_metadata.iiPDF_eta,    # In
            hm_metadata.iiPDF_w, hm_metadata.iiPDF_Ncn,      # In
        )
        corr_array_cloud = mirror_lower_triangular_matrix(corr_array_cloud[0, 0])
        corr_array_below = mirror_lower_triangular_matrix(corr_array_below[0, 0])

        report = []
        for label, matrix in (("in cloud", corr_array_cloud), ("below cloud", corr_array_below)):
            report.append(f"Correlation array (approximate); {label}:")
            for row in np.asarray(matrix):
                report.append("".join(f"{value:7.3f}" for value in row))
        for line in report:
            print(line)
        if case_info_file is not None:
            with Path(case_info_file).open("r+") as report_file:
                report_file.seek(0, 2)
                for line in report:
                    print(line, file=report_file)

    if parameters.lh_microphys_type != parameters.lh_microphys_disabled:
        # These deterministic inputs do not replace independent importance/start-level draws.
        if silhs_config_flags.l_lh_deterministic_test and (
            silhs_config_flags.l_lh_importance_sampling or silhs_config_flags.l_random_k_lh_start
        ):
            raise ValueError(
                "Deterministic SILHS testing requires importance sampling and random starts disabled"
            )
        if silhs_config_flags.l_lh_importance_sampling and parameters.lh_sequence_length != 1:
            raise ValueError("SILHS importance sampling requires lh_sequence_length = 1")
        if silhs_config_flags.cluster_allocation_strategy not in (1, 2, 3):
            raise ValueError("Unsupported SILHS cluster allocation strategy")
        if silhs_config_flags.l_Lscale_vert_avg:
            raise ValueError("l_Lscale_vert_avg is deprecated in Fortran")
    vert_decorr_coef_out = parameters_silhs.vert_decorr_coef
    return (
        hydromet_dim,
        pdf_dim,
        hm_metadata,
        silhs_config_flags,
        vert_decorr_coef_out,
        corr_array_n_cloud,
        corr_array_n_below,
    )


# -----------------------------------------------------------------------------
def cleanup_microphys():
    """Release module-owned sedimentation flags; case arrays follow JAX lifetime.
    """
    # Description:
    # De-allocate arrays used by the microphysics
    # References:
    # None
    # -----------------------------------------------------------------------
    parameters.l_hydromet_sed = ()
