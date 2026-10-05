"""Structural checks for the Fortran-to-JAX forcing mirrors."""

from __future__ import annotations
from utilities.output_paths import REPO_ROOT as _REPO_ROOT

import ast


BENCHMARK_CASES = _REPO_ROOT / "clubb_jax" / "src" / "Benchmark_cases"
INPUT_FIELDS = _REPO_ROOT / "clubb_jax" / "src" / "Input_fields"
CLUBB_CORE = _REPO_ROOT / "clubb_jax" / "src" / "CLUBB_core"
ADVANCE_CLUBB_TO_END = _REPO_ROOT / "clubb_jax" / "src" / "advance_clubb_to_end.py"


# Source order for the routines supported by the standalone JAX driver. The
# dycore-only forcing routine is intentionally excluded as documented at the
# top of time_dependent_input.py.
SUPPORTED_SOURCE_ROUTINE_ORDER = {
    "arm.py": ["arm_sfclyr"],
    "arm_0003.py": ["arm_0003_sfclyr"],
    "arm_3year.py": ["arm_3year_sfclyr"],
    "arm_97.py": ["arm_97_sfclyr"],
    "astex_a209.py": ["astex_a209_tndcy", "astex_a209_sfclyr"],
    "atex.py": ["calc_forcings", "atex_tndcy", "atex_sfclyr"],
    "atex_long.py": ["calc_forcings", "atex_long_tndcy", "atex_long_sfclyr"],
    "bomex.py": ["bomex_tndcy", "bomex_sfclyr"],
    "clex9_oct14.py": ["clex9_oct14_read_t_dependent"],
    "cloud_feedback.py": ["cloud_feedback_sfclyr"],
    "cobra.py": ["cobra_sfclyr"],
    "diag_ustar_module.py": ["diag_ustar"],
    "dycoms2_rf01.py": ["dycoms2_rf01_tndcy", "dycoms2_rf01_sfclyr"],
    "dycoms2_rf02.py": ["dycoms2_rf02_tndcy", "dycoms2_rf02_sfclyr"],
    "ekman.py": ["ekman_sfclyr"],
    "fire.py": ["fire_sfclyr"],
    "gabls2.py": ["gabls2_tndcy", "gabls2_sfclyr"],
    "gabls3.py": ["gabls3_sfclyr"],
    "gabls3_night.py": ["gabls3_night_sfclyr", "psi_h", "gm1", "gh1", "fm1", "fh1", "landflx"],
    "jun25.py": ["jun25_altocu_read_t_dependent"],
    "lba.py": ["lba_tndcy", "lba_sfclyr"],
    "mpace_a.py": ["mpace_a_tndcy", "mpace_a_sfclyr", "mpace_a_init"],
    "mpace_b.py": ["mpace_b_tndcy", "mpace_b_sfclyr"],
    "neutral_case.py": ["neutral_case_sfclyr"],
    "nov11.py": ["nov11_altocu_rtm_adjust", "nov11_altocu_read_t_dependent"],
    "prescribe_forcings.py": ["prescribe_forcings", "read_surface_var_for_bc"],
    "rico.py": ["rico_tndcy", "rico_sfclyr"],
    "sfc_flux.py": [
        "compute_momentum_flux",
        "compute_ubar",
        "compute_ht_mostr_flux",
        "compute_wpthlp_sfc",
        "compute_wprtp_sfc",
        "set_sclr_sfc_rtm_thlm",
        "convert_sens_ht_to_km_s",
        "convert_latent_ht_to_m_s",
    ],
    "spec_hum_to_mixing_ratio.py": [
        "flux_spec_hum_to_mixing_ratio",
        "force_spec_hum_to_mixing_ratio",
    ],
    "time_dependent_input.py": [
        "initialize_t_dependent_input",
        "finalize_t_dependent_input",
        "initialize_t_dependent_sfc",
        "initialize_t_dependent_forcings",
        "finalize_t_dependent_forcings",
        "finalize_t_dependent_sfc",
        "read_to_grid",
        "apply_time_dependent_forcings_from_array",
        "apply_time_dependent_forcings",
        "time_select",
    ],
    "twp_ice.py": ["twp_ice_sfclyr"],
    "wangara.py": ["wangara_tndcy", "wangara_sfclyr"],
}

PRESCRIBE_FORCINGS_ARGUMENTS = [
    "gr",
    "nzm",
    "nzt",
    "ngrdcol",
    "sclr_dim",
    "edsclr_dim",
    "sclr_idx",
    "runtype",
    "sfctype",
    "time_current",
    "time_initial",
    "dt",
    "um",
    "vm",
    "thlm",
    "p_in_Pa",
    "exner",
    "rho",
    "rho_zm",
    "thvm",
    "veg_T_in_K",
    "l_modify_bc_for_cnvg_test",
    "saturation_formula",
    "stats",
    "rtm",
    "wm_zm",
    "wm_zt",
    "ug",
    "vg",
    "um_ref",
    "vm_ref",
    "thlm_forcing",
    "rtm_forcing",
    "um_forcing",
    "vm_forcing",
    "wprtp_forcing",
    "wpthlp_forcing",
    "rtp2_forcing",
    "thlp2_forcing",
    "rtpthlp_forcing",
    "wpsclrp",
    "sclrm_forcing",
    "edsclrm_forcing",
    "wpthlp_sfc",
    "wprtp_sfc",
    "upwp_sfc",
    "vpwp_sfc",
    "T_sfc",
    "p_sfc",
    "sens_ht",
    "latent_ht",
    "wpsclrp_sfc",
    "wpedsclrp_sfc",
    "err_info",
]

PRESCRIBE_FORCINGS_RESULTS = [
    "stats",
    "rtm",
    "wm_zm",
    "wm_zt",
    "ug",
    "vg",
    "um_ref",
    "vm_ref",
    "thlm_forcing",
    "rtm_forcing",
    "um_forcing",
    "vm_forcing",
    "wprtp_forcing",
    "wpthlp_forcing",
    "rtp2_forcing",
    "thlp2_forcing",
    "rtpthlp_forcing",
    "wpsclrp",
    "sclrm_forcing",
    "edsclrm_forcing",
    "wpthlp_sfc",
    "wprtp_sfc",
    "upwp_sfc",
    "vpwp_sfc",
    "T_sfc",
    "p_sfc",
    "sens_ht",
    "latent_ht",
    "wpsclrp_sfc",
    "wpedsclrp_sfc",
    "err_info",
]


HOST_BOUNDARY_FILES = {
    "mpace_a.py",  # mpace_a_init reads the source data files.
    "time_dependent_input.py",  # Source-equivalent input initialization and parsing.
}
DEVICE_PHYSICS_FILES = tuple(
    filename
    for filename in SUPPORTED_SOURCE_ROUTINE_ORDER
    if filename not in HOST_BOUNDARY_FILES
)


def test_device_physics_does_not_import_numpy_or_raw_f2py():
    for filename in DEVICE_PHYSICS_FILES:
        source = (BENCHMARK_CASES / filename).read_text()
        tree = ast.parse(source)
        imported_modules = {
            alias.name
            for node in tree.body
            if isinstance(node, ast.Import)
            for alias in node.names
        }
        imported_from = {
            node.module or ""
            for node in tree.body
            if isinstance(node, ast.ImportFrom)
        }
        assert "numpy" not in imported_modules, filename
        assert not any("clubb_f2py" in name for name in imported_modules | imported_from), filename
