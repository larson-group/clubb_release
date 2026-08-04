"""Human-maintained CLUBB settings semantics and no-build validation.

This is the single review surface for settings that are exposed by the Run
tab, Tune, CLI namelists, and the local agent interface.  It deliberately
does *not* scan source files or generate a second coverage product: when a
setting is added or its ownership changes, update the short audit lists and
the explicit rules below in the same review as the Fortran change.

Fortran remains the runtime oracle.  ``parameters_tunable.F90`` owns the
parameter order and hard rejection limits, while ``numerical_check.F90`` and
the scientific owner modules decide whether an actual run may proceed.  This
module is the dependency-free mirror used to make tooling responsive before a
Python API build exists.  The compiled test in
``tests/run_clubb_settings_validation_test.py`` checks the common rejecting
cases and all parameter bounds after ``./compile.py -python``.

The audit is intentionally split into two readable parts:

* ``PARAMETER_PRIMARY_OWNERS`` and ``CONFIG_FLAG_PRIMARY_OWNERS`` record the
  one-time review of every public tunable/configurable setting.  An owner is a
  primary implementation location, not a claim that no other routine reads
  the value.
* the constants headed ``MODE/CONSTRAINT RULES`` record the settings that
  change whether another setting has an effect, must match another value, or
  is incompatible.  These constants drive ``resolve_clubb_settings``; the
  comments are therefore the maintained human-readable rule list.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math
import sys
from typing import Any

# Keep this list in exactly the order of ``params_list`` in
# ``src/CLUBB_core/parameters_tunable.F90``.  The order is part of the F2PY
# and tuner interface, so this intentionally remains a plain, inspectable
# tuple rather than being discovered dynamically.
PARAMETER_NAMES = tuple(
    """
    C1 C1b C1c C2rt C2thl C2rtthl C4 C_uu_shr C_uu_buoy C6rt C6rtb C6rtc
    C6thl C6thlb C6thlc C7 C7b C7c C8 C8b C10 C11 C11b C11c C12 C13 C14
    C_wp2_pr_dfsn C_wp3_pr_tp C_wp3_pr_turb C_wp3_pr_dfsn C_wp2_splat
    C6rt_Lscale0 C6thl_Lscale0 C7_Lscale0 wpxp_L_thresh c_K c_K1 nu1 c_K2
    nu2 c_K6 nu6 c_K8 nu8 c_K9 nu9 nu10 c_K_hm c_K_hmb K_hm_min_coef nu_hm
    slope_coef_spread_DG_means_w pdf_component_stdev_factor_w
    coef_spread_DG_means_rt coef_spread_DG_means_thl gamma_coef gamma_coefb
    gamma_coefc mu beta lmin_coef omicron zeta_vrnce_rat
    upsilon_precip_frac_rat lambda0_stability_coef mult_coef taumin taumax
    Lscale_mu_coef Lscale_pert_coef alpha_corr Skw_denom_coef c_K10 c_K10h
    thlp2_rad_coef thlp2_rad_cloud_frac_thresh up2_sfc_coef Skw_max_mag
    C_invrs_tau_bkgnd C_invrs_tau_sfc C_invrs_tau_shear C_invrs_tau_N2
    C_invrs_tau_N2_wp2 C_invrs_tau_N2_xp2 C_invrs_tau_N2_wpxp
    C_invrs_tau_N2_clear_wp3 C_invrs_tau_wpxp_Ri C_invrs_tau_wpxp_N2_thresh
    xp3_coef_base xp3_coef_slope altitude_threshold rtp2_clip_coef Cx_min
    Cx_max Richardson_num_min Richardson_num_max a3_coef_min a_const bv_efold
    wpxp_Ri_exp z_displace
    """.split()
)


# ---------------------------------------------------------------------------
# HAND-MAINTAINED SETTINGS AUDIT
# ---------------------------------------------------------------------------
# Every public parameter appears in exactly one primary-owner row.  These are
# review notes, not an executable source scanner.  Values may be passed
# through several routines; the source listed is the first scientific owner to
# inspect when changing a setting's semantics.
PARAMETER_PRIMARY_OWNERS: dict[str, tuple[str, tuple[str, ...]]] = {
    "retained unused parameters": (
        "src/CLUBB_core/parameters_tunable.F90",
        ("C10", "C13"),
    ),
    "wp2/wp3 closure and pressure": (
        "src/CLUBB_core/advance_wp2_wp3_module.F90",
        (
            "C1", "C1b", "C1c", "C_uu_shr", "C_uu_buoy", "C8", "C8b",
            "C11", "C11b", "C11c", "C12", "C14",
            "C_wp2_pr_dfsn", "C_wp3_pr_tp",
            "C_wp3_pr_turb", "C_wp3_pr_dfsn", "C_wp2_splat", "a3_coef_min",
            "Skw_max_mag",
        ),
    ),
    "wpxp closure and scalar skewness": (
        "src/CLUBB_core/advance_xm_wpxp_module.F90",
        (
            "C2rt", "C2thl", "C2rtthl", "C6rt", "C6rtb", "C6rtc",
            "C6thl", "C6thlb", "C6thlc", "C7", "C7b", "C7c",
            "C6rt_Lscale0", "C6thl_Lscale0", "C7_Lscale0",
            "wpxp_L_thresh", "altitude_threshold", "alpha_corr",
            "Skw_denom_coef", "Cx_min", "Cx_max", "Richardson_num_min",
            "Richardson_num_max", "wpxp_Ri_exp",
        ),
    ),
    "xp2/xpyp and scalar transport": (
        "src/CLUBB_core/advance_xp2_xpyp_module.F90",
        ("C4", "beta", "rtp2_clip_coef"),
    ),
    "diffusion and wind moments": (
        "src/CLUBB_core/advance_*_module.F90",
        (
            "c_K", "c_K1", "nu1", "c_K2", "nu2",
            "c_K6", "nu6", "c_K8", "nu8", "c_K9", "nu9", "nu10",
            "c_K_hm", "c_K_hmb", "K_hm_min_coef", "nu_hm", "c_K10",
            "c_K10h", "up2_sfc_coef",
        ),
    ),
    "mixing length and inverse tau": (
        "src/CLUBB_core/mixing_length.F90",
        (
            "mu", "lmin_coef", "lambda0_stability_coef", "mult_coef",
            "taumin", "taumax", "Lscale_mu_coef", "Lscale_pert_coef",
            "C_invrs_tau_bkgnd", "C_invrs_tau_sfc", "C_invrs_tau_shear",
            "C_invrs_tau_N2", "C_invrs_tau_N2_wp2", "C_invrs_tau_N2_xp2",
            "C_invrs_tau_N2_wpxp", "C_invrs_tau_N2_clear_wp3",
            "C_invrs_tau_wpxp_Ri", "C_invrs_tau_wpxp_N2_thresh", "z_displace",
            "bv_efold",
        ),
    ),
    "PDF and cloud diagnostics": (
        "src/CLUBB_core/pdf_closure_module.F90",
        (
            "omicron", "zeta_vrnce_rat", "upsilon_precip_frac_rat",
            "gamma_coef", "gamma_coefb", "gamma_coefc", "a_const",
        ),
    ),
    "radiative / xp3 auxiliary closures": (
        "src/CLUBB_core/advance_helper_module.F90",
        (
            "thlp2_rad_coef", "thlp2_rad_cloud_frac_thresh", "xp3_coef_base",
            "xp3_coef_slope",
        ),
    ),
    "PDF type 4 (new PDF)": (
        "src/CLUBB_core/new_pdf_main.F90",
        ("coef_spread_DG_means_rt", "coef_spread_DG_means_thl"),
    ),
    "PDF type 7 (hybrid PDF)": (
        "src/CLUBB_core/new_hybrid_pdf_main.F90",
        ("slope_coef_spread_DG_means_w", "pdf_component_stdev_factor_w"),
    ),
}

# Every ConfigFlags member in ``model_flags.F90`` is listed once.  These
# concise groups deliberately record the subsystem rather than duplicating
# model_flags' excellent per-flag prose comments.
CONFIG_FLAG_PRIMARY_OWNERS: dict[str, tuple[str, tuple[str, ...]]] = {
    "PDF, cloud, and correlation choices": (
        "src/CLUBB_core/{pdf_closure_module,setup_clubb_pdf_params}.F90",
        (
            "iiPDF_type", "ipdf_call_placement", "saturation_formula",
            "l_use_precip_frac", "l_use_cloud_cover", "l_diagnose_correlations",
            "l_calc_w_corr", "l_fix_w_chi_eta_correlations", "l_C2_cloud_frac",
            "l_call_pdf_closure_twice", "l_vert_avg_closure",
            "l_trapezoidal_rule_zt", "l_trapezoidal_rule_zm",
        ),
    ),
    "moment advancement, closures, and clipping": (
        "src/CLUBB_core/advance_*_module.F90",
        (
            "l_predict_upwp_vpwp", "l_min_wp2_from_corr_wx",
            "l_min_xp2_from_corr_wx", "l_enable_relaxed_clipping",
            "l_diffuse_rtm_and_thlm", "l_do_expldiff_rtm_thlm",
            "l_upwind_xpyp_ta", "l_upwind_xm_ma", "l_standard_term_ta",
            "l_partial_upwind_wp3", "l_godunov_upwind_wpxp_ta",
            "l_godunov_upwind_xpyp_ta", "l_damp_wp2_using_em",
            "l_damp_wp3_Skw_squared", "l_use_tke_in_wp3_pr_turb_term",
            "l_use_tke_in_wp2_wp3_K_dfsn", "l_use_wp3_lim_with_smth_Heaviside",
            "l_modify_limiters_for_cnvg_test", "l_smooth_Heaviside_tau_wpxp",
            "l_wp2_fill_holes_tke", "l_linearize_pbl_winds",
            "l_mono_flux_lim_thlm", "l_mono_flux_lim_rtm", "l_mono_flux_lim_um",
            "l_mono_flux_lim_vm", "l_mono_flux_lim_spikefix",
        ),
    ),
    "mixing length and stability": (
        "src/CLUBB_core/mixing_length.F90",
        (
            "l_stability_correct_Kh_N2_zm", "l_stability_correct_tau_zm",
            "l_diag_Lscale_from_tau", "l_Lscale_plume_centered",
            "l_Lscale_vert_avg", "l_use_C7_Richardson", "l_use_C11_Richardson",
            "l_use_shear_Richardson", "l_brunt_vaisala_freq_moist",
            "l_use_thvm_in_bv_freq", "l_prescribed_avg_deltaz",
            "l_vary_convect_depth",
        ),
    ),
    "dynamics, grid, and host coupling": (
        "src/CLUBB_core/{clubb_driver,grid_*}.F90",
        (
            "l_ho_nontrad_coriolis", "l_ho_trad_coriolis", "l_uv_nudge",
            "l_rtm_nudge", "l_host_applies_sfc_fluxes", "l_add_dycore_grid",
            "l_lmm_stepping", "l_e3sm_config", "grid_remap_method",
            "grid_adapt_in_time_method", "fill_holes_type", "penta_solve_method",
            "tridiag_solve_method",
        ),
    ),
    "microphysics and SILHS": (
        "src/Microphys/microphys_init_cleanup.F90",
        (
            "l_calc_thlp2_rad", "l_rcm_supersat_adj", "l_const_Nc_in_cloud",
            "l_lh_clustered_sampling", "l_lh_importance_sampling",
            "l_lh_instant_var_covar_src", "l_lh_limit_weights",
            "l_lh_normalize_weights", "l_lh_straight_mc", "l_lh_var_frac",
            "l_max_overlap_in_cloud", "l_random_k_lh_start",
            "l_rcm_in_cloud_k_lh_start", "cluster_allocation_strategy",
        ),
    ),
}

# Keep this in the order printed by the configurable-flags namelist parser.
# It is deliberately explicit, so a reviewer must classify a new ConfigFlags
# member instead of relying on an automatic source scan to guess its meaning.
CONFIG_FLAG_NAMES = tuple(
    """
    cluster_allocation_strategy fill_holes_type grid_adapt_in_time_method
    grid_remap_method iiPDF_type ipdf_call_placement l_add_dycore_grid
    l_brunt_vaisala_freq_moist l_C2_cloud_frac l_call_pdf_closure_twice
    l_damp_wp2_using_em l_damp_wp3_Skw_squared l_diag_Lscale_from_tau
    l_diffuse_rtm_and_thlm l_do_expldiff_rtm_thlm l_e3sm_config
    l_enable_relaxed_clipping l_godunov_upwind_wpxp_ta
    l_godunov_upwind_xpyp_ta l_ho_nontrad_coriolis l_ho_trad_coriolis
    l_host_applies_sfc_fluxes l_lh_clustered_sampling l_lh_importance_sampling
    l_lh_instant_var_covar_src l_lh_limit_weights l_lh_normalize_weights
    l_lh_straight_mc l_lh_var_frac l_linearize_pbl_winds l_lmm_stepping
    l_Lscale_plume_centered l_Lscale_vert_avg l_max_overlap_in_cloud
    l_min_wp2_from_corr_wx l_min_xp2_from_corr_wx l_modify_limiters_for_cnvg_test
    l_mono_flux_lim_rtm l_mono_flux_lim_spikefix l_mono_flux_lim_thlm
    l_mono_flux_lim_um l_mono_flux_lim_vm l_partial_upwind_wp3
    l_predict_upwp_vpwp l_prescribed_avg_deltaz l_random_k_lh_start
    l_rcm_in_cloud_k_lh_start l_rcm_supersat_adj l_smooth_Heaviside_tau_wpxp
    l_stability_correct_Kh_N2_zm l_stability_correct_tau_zm l_standard_term_ta
    l_trapezoidal_rule_zm l_trapezoidal_rule_zt l_upwind_xm_ma
    l_upwind_xpyp_ta l_use_C11_Richardson l_use_C7_Richardson
    l_use_cloud_cover l_use_precip_frac l_use_shear_Richardson
    l_use_thvm_in_bv_freq l_use_tke_in_wp2_wp3_K_dfsn
    l_use_tke_in_wp3_pr_turb_term l_use_wp3_lim_with_smth_Heaviside
    l_vary_convect_depth l_vert_avg_closure l_wp2_fill_holes_tke
    penta_solve_method saturation_formula tridiag_solve_method
    """.split()
)

# These are owned outside ConfigFlags but are part of a complete settings
# decision.  Keep them visible because several rules below require the extra
# runtime/build context.
EXTERNAL_SETTING_OWNERS: dict[str, str] = {
    "l_input_fields": "src/clubb_driver.F90",
    "l_implemented": "src/clubb_driver.F90",
    "l_explicit_turbulent_adv_wpxp": "src/CLUBB_core/model_flags.F90",
    "l_avg_Lscale": "src/CLUBB_core/mixing_length.F90",
    "l_predict_Nc": "src/Microphys/microphys_init_cleanup.F90",
    "microphys_scheme": "src/Microphys/microphys_init_cleanup.F90",
    "rad_scheme": "src/Radiation/radiation_module.F90",
    "l_silhs_rad": "src/Radiation/radiation_module.F90",
    "l_soil_veg": "src/Radiation/radiation_module.F90",
}


def audit_gaps() -> dict[str, tuple[str, ...]]:
    """Return audit-list omissions/duplicates; used by a small unit test.

    This is intentionally *not* a source scanner.  It simply prevents edits
    to this maintained list from quietly dropping or duplicating a public
    parameter or configurable flag.
    """
    parameter_seen = [name for _, names in PARAMETER_PRIMARY_OWNERS.values() for name in names]
    flag_seen = [name for _, names in CONFIG_FLAG_PRIMARY_OWNERS.values() for name in names]
    return {
        "missing_parameters": tuple(sorted(set(PARAMETER_NAMES) - set(parameter_seen))),
        "duplicate_parameters": tuple(sorted(name for name in set(parameter_seen) if parameter_seen.count(name) > 1)),
        "missing_flags": tuple(sorted(set(CONFIG_FLAG_NAMES) - set(flag_seen))),
        "duplicate_flags": tuple(sorted(name for name in set(flag_seen) if flag_seen.count(name) > 1)),
    }


# ---------------------------------------------------------------------------
# MODE / CONSTRAINT RULES
# ---------------------------------------------------------------------------
# All rule IDs below are cited in ``parameter_activity`` and validation
# messages.  This is the human-readable list to update when Fortran semantics
# change.  "inactive" means the supplied value has no effect; it is not an
# error.  Requirements and conflicts remain errors in ``validate_clubb_settings``.

# Retained in the namelist for compatibility, but no model tendency reads them.
UNUSED_PARAMETERS = ("C10", "C13")

# Each group is physically equal in all valid configurations.  The UI and
# tuner represent a group as one coordinate and expand it before execution.
LINKED_PARAMETER_GROUPS = (
    ("C6rt", "C6thl"),
    ("C6rtb", "C6thlb"),
    ("C6rtc", "C6thlc"),
    ("C6rt_Lscale0", "C6thl_Lscale0"),
)

# Parameters unique to one PDF implementation.  They are inactive for every
# other iiPDF_type, not merely discouraged.
PDF_PARAMETER_OWNERS = {
    4: ("coef_spread_DG_means_rt", "coef_spread_DG_means_thl"),
    7: ("slope_coef_spread_DG_means_w", "pdf_component_stdev_factor_w"),
}

# l_use_C7_Richardson replaces the entire tunable C7 function, including its
# direct-Lscale damping length.  l_use_C11_Richardson analogously replaces the
# C11 function used by wp3 closure.
RICHARDSON_REPLACED_PARAMETERS = {
    "l_use_C7_Richardson": ("C7", "C7b", "C7c", "C7_Lscale0"),
    "l_use_C11_Richardson": ("C11", "C11b", "C11c"),
}

# Direct Lscale and inverse-tau Lscale are mutually exclusive implementations.
# The tau coefficients are read only by diagnose_Lscale_from_tau.  Conversely,
# these direct-Lscale perturbation controls and C6 Lscale damping lengths are
# bypassed when tau diagnoses Lscale.  C7's altitude/threshold damping remains
# separate and is therefore deliberately not listed here.
INVERSE_TAU_PARAMETERS = (
    "C_invrs_tau_bkgnd", "C_invrs_tau_sfc", "C_invrs_tau_shear",
    "C_invrs_tau_N2", "C_invrs_tau_N2_wp2", "C_invrs_tau_N2_xp2",
    "C_invrs_tau_N2_wpxp", "C_invrs_tau_N2_clear_wp3",
    "C_invrs_tau_wpxp_Ri", "C_invrs_tau_wpxp_N2_thresh",
)
DIRECT_LSCALE_ONLY_PARAMETERS = (
    "Lscale_mu_coef", "Lscale_pert_coef", "C6rt_Lscale0", "C6thl_Lscale0",
)

# These options only enter the indicated scientific branch.
BRANCH_ONLY_PARAMETERS = {
    "l_use_precip_frac": ("omicron", "zeta_vrnce_rat"),
    "l_calc_thlp2_rad": ("thlp2_rad_coef", "thlp2_rad_cloud_frac_thresh"),
}

# Flag-only consequences discovered in the same audit.  These are kept as
# notes rather than forcibly rewritten: some are compile/build-context
# dependent, and CLUBB's runtime owner remains responsible for any fatal
# decision.  A UI may show an entry as ``no effect`` using this information.
FLAG_EFFECT_RULES = (
    (
        "l_partial_upwind_wp3",
        "has no effect unless iiPDF_type=1 (ADG1) and l_standard_term_ta=true",
        "src/CLUBB_core/advance_wp2_wp3_module.F90",
    ),
    (
        "l_Lscale_plume_centered",
        "requires vertical Lscale averaging; the owner rejects it when l_avg_Lscale=false",
        "src/CLUBB_core/mixing_length.F90",
    ),
    (
        "l_predict_upwp_vpwp",
        "requires implicit wpxp turbulent advection and PDF 1, 7, or 9",
        "src/CLUBB_core/numerical_check.F90",
    ),
    (
        "l_min_xp2_from_corr_wx / l_enable_relaxed_clipping",
        "must have opposite values; both false is retained as a Fortran warning",
        "src/CLUBB_core/numerical_check.F90",
    ),
    (
        "l_damp_wp2_using_em",
        "requires C1=C14 and l_stability_correct_tau_zm=false",
        "src/CLUBB_core/numerical_check.F90",
    ),
)

# Each entry is [minimum, maximum], matching Fortran's 2 x nparams table.
# ``tiny(1._core_rknd)`` and ``1-epsilon(1._core_rknd)`` are encoded here with
# Python's IEEE double equivalents; the parity test verifies this assumption.
_TINY = sys.float_info.min
_ONE_MINUS_EPS = 1.0 - sys.float_info.epsilon
PARAMETER_HARD_BOUNDS: tuple[tuple[float | None, float | None], ...] = tuple(
    (0.0, None) for _ in PARAMETER_NAMES
)


def _with_bounds(
    bounds: tuple[tuple[float | None, float | None], ...],
    names: str,
    low: float | None,
    high: float | None,
) -> tuple[tuple[float | None, float | None], ...]:
    updated = list(bounds)
    for name in names.split():
        updated[PARAMETER_NAMES.index(name)] = (low, high)
    return tuple(updated)


PARAMETER_HARD_BOUNDS = _with_bounds(PARAMETER_HARD_BOUNDS, "C_uu_shr C_uu_buoy C7 C7b C11 C11b upsilon_precip_frac_rat", 0.0, 1.0)
PARAMETER_HARD_BOUNDS = _with_bounds(PARAMETER_HARD_BOUNDS, "slope_coef_spread_DG_means_w pdf_component_stdev_factor_w", _TINY, None)
PARAMETER_HARD_BOUNDS = _with_bounds(PARAMETER_HARD_BOUNDS, "coef_spread_DG_means_rt coef_spread_DG_means_thl", 0.0, _ONE_MINUS_EPS)
PARAMETER_HARD_BOUNDS = _with_bounds(PARAMETER_HARD_BOUNDS, "beta", 0.0, 3.0)
PARAMETER_HARD_BOUNDS = _with_bounds(PARAMETER_HARD_BOUNDS, "omicron", _TINY, 1.0)
PARAMETER_HARD_BOUNDS = _with_bounds(PARAMETER_HARD_BOUNDS, "zeta_vrnce_rat", -1.0 + sys.float_info.epsilon, None)
PARAMETER_HARD_BOUNDS = _with_bounds(PARAMETER_HARD_BOUNDS, "a3_coef_min", 1.0, 3.0)


def parameter_hard_bounds() -> dict[str, dict[str, float | None]]:
    """Return UI-friendly bounds keyed by canonical parameter name."""
    return {
        name: {"min": low, "max": high}
        for name, (low, high) in zip(PARAMETER_NAMES, PARAMETER_HARD_BOUNDS, strict=True)
    }


def parse_finite_number(value: Any) -> float | None:
    """Parse one Fortran-style scalar, returning ``None`` for non-numbers."""
    try:
        parsed = float(str(value).strip().replace("D", "E").replace("d", "e"))
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def format_setting_value(value: Any) -> str:
    """Return a stable namelist/display spelling for one scalar setting."""
    if isinstance(value, bool):
        return ".true." if value else ".false."
    return str(value).strip()


def settings_values_equal(left: Any, right: Any, *, boolean: bool = False) -> bool:
    """Compare settings once, including Fortran numeric spelling and tolerance."""
    if boolean:
        return _as_bool(left) == _as_bool(right)
    left_number = parse_finite_number(left)
    right_number = parse_finite_number(right)
    if left_number is not None and right_number is not None:
        return _fortran_equal(left_number, right_number)
    return format_setting_value(left).lower() == format_setting_value(right).lower()


def default_range_for_value(
    default: float,
    hard_bounds: Mapping[str, float | None] | None,
) -> tuple[float, float]:
    """Return the historical finite Tune interval clipped to hard bounds."""
    default = float(default)
    relative_min, relative_max = sorted((default / 4.0, default * 4.0))
    lower = (hard_bounds or {}).get("min")
    upper = (hard_bounds or {}).get("max")
    low = float(lower) if lower is not None else relative_min
    high = float(upper) if upper is not None else (4.0 if default == 0.0 else relative_max)
    if low > high:
        raise ValueError(f"No feasible default Tune interval: lower {low:g} exceeds upper {high:g}")
    return low, high


def default_tunable_ranges(
    defaults: Mapping[str, Any],
    hard_bounds: Mapping[str, Mapping[str, float | None]] | None = None,
) -> dict[str, dict[str, float]]:
    """Derive one default/hard-bound-aware Tune interval per numeric setting."""
    bounds = hard_bounds or parameter_hard_bounds()
    ranges: dict[str, dict[str, float]] = {}
    for raw_name, raw_value in dict(defaults or {}).items():
        name = canonical_parameter_name(raw_name)
        value = parse_finite_number(raw_value)
        if value is None:
            continue
        low, high = default_range_for_value(value, bounds.get(name))
        ranges[name] = {"default": value, "min": low, "max": high}
    return ranges


def build_settings_schema(
    flag_defaults: Mapping[str, Any],
    parameter_defaults: Mapping[str, Mapping[str, Any]],
    parameter_metadata: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Build a JSON-safe, immutable-by-convention settings interpretation schema.

    Tabs may cache this per selected tunable config.  It deliberately contains
    scientific semantics and defaults but no Dash components or callback IDs.
    """
    normalized_flags = _canonicalize_setting_names(flag_defaults, _FLAG_NAME_ALIASES)
    normalized_parameters = {
        str(file_name): _canonicalize_setting_names(values, _PARAMETER_NAME_ALIASES)
        for file_name, values in dict(parameter_defaults or {}).items()
    }
    metadata = [
        {"file": str(item.get("file") or ""), "name": str(item.get("name") or "")}
        for item in parameter_metadata or []
        if item.get("file") and item.get("name")
    ]
    tunable_defaults = normalized_parameters.get("tunable", {})
    return {
        "flag_defaults": normalized_flags,
        "parameter_defaults": normalized_parameters,
        "parameter_metadata": metadata,
        "hard_bounds": parameter_hard_bounds(),
        "tunable_default_ranges": default_tunable_ranges(tunable_defaults),
    }


def setting_key(file_name: str, name: str) -> str:
    """Return the stable UI-neutral key for a scalar configuration setting."""
    return f"{str(file_name)}:{str(name)}"


def values_by_setting_key(component_ids: Sequence[Mapping[str, Any]], values: Sequence[Any]) -> dict[str, Any]:
    """Map UI-neutral ``file``/``name`` identifiers to values without DOM order."""
    return {
        setting_key(component_id.get("file"), component_id.get("name")): value
        for component_id, value in zip(component_ids or [], values or [])
        if isinstance(component_id, Mapping) and component_id.get("file") and component_id.get("name")
    }


def values_by_name(
    component_ids: Sequence[Mapping[str, Any]],
    values: Sequence[Any],
    *,
    name_key: str = "name",
) -> dict[str, Any]:
    """Map component-like identifiers to values without assuming DOM order."""
    return {
        str(component_id.get(name_key)): value
        for component_id, value in zip(component_ids or [], values or [])
        if isinstance(component_id, Mapping) and component_id.get(name_key)
    }


@dataclass(frozen=True)
class ValidationIssue:
    """One stable, display-friendly validation finding."""

    code: str
    message: str
    severity: str = "error"
    setting: str | None = None
    column: int | None = None


@dataclass(frozen=True)
class SettingsResolution:
    """Declarative consequences of a partial effective CLUBB configuration.

    ``forced_parameters`` contains only required equal-value followers.
    ``forced_flags`` remains a compatibility field and is intentionally empty:
    flag conflicts are reported rather than rewritten.  ``allowed_flag_values``
    represents a restriction without an arbitrary choice.  Parameter state is
    deliberately conservative: only parameters known to be exclusive to a PDF
    family are disabled.
    """

    effective_flags: dict[str, Any]
    forced_flags: dict[str, Any]
    forced_parameters: dict[str, Any]
    allowed_flag_values: dict[str, list[Any]]
    conflicting_setting_keys: list[str]
    flag_states: dict[str, dict[str, str]]
    flag_relationships: list[dict[str, Any]]
    parameter_states: dict[str, dict[str, str]]
    coupled_parameters: list[dict[str, Any]]
    issues: tuple[ValidationIssue, ...]

    def as_dict(self) -> dict[str, Any]:
        """Return JSON-safe structured data for Dash, Tune, and MCP."""
        return {
            "effective_flags": dict(self.effective_flags),
            "forced_flags": dict(self.forced_flags),
            "forced_parameters": dict(self.forced_parameters),
            "allowed_flag_values": {name: list(values) for name, values in self.allowed_flag_values.items()},
            "conflicting_setting_keys": list(self.conflicting_setting_keys),
            "flag_states": {name: dict(state) for name, state in self.flag_states.items()},
            "flag_relationships": [dict(item) for item in self.flag_relationships],
            "parameter_states": {name: dict(state) for name, state in self.parameter_states.items()},
            "coupled_parameters": [dict(item) for item in self.coupled_parameters],
            "issues": [issue.__dict__ for issue in self.issues],
        }


# Keep boolean relationships declarative so the Run and Tune interfaces can
# present the same constraints without re-encoding them in layout callbacks.
# ``required_values`` scales to rules in which one switch constrains several
# companions; the current UI preserves one checkbox per member.
_FLAG_RELATIONSHIPS = (
    {
        "members": ("l_min_xp2_from_corr_wx", "l_enable_relaxed_clipping"),
        "relation": "opposite",
        "label": "exclusive",
        "description": "Enabling either xp2 clipping mode clears the other.",
    },
    {
        "members": ("l_damp_wp2_using_em", "l_stability_correct_tau_zm"),
        "relation": "requires_value",
        "required_values": {"l_damp_wp2_using_em": True, "l_stability_correct_tau_zm": False},
        "switch_on_true": {
            "l_damp_wp2_using_em": {"l_stability_correct_tau_zm": False},
            "l_stability_correct_tau_zm": {"l_damp_wp2_using_em": False},
        },
        "label": "linked",
        "description": "EM wp2 damping requires stability correction off.",
    },
)

_FLAG_NAME_ALIASES = {
    "iipdf_type": "iiPDF_type",
    "l_predict_upwp_vpwp": "l_predict_upwp_vpwp",
    "l_damp_wp2_using_em": "l_damp_wp2_using_em",
    "l_stability_correct_tau_zm": "l_stability_correct_tau_zm",
    "l_min_xp2_from_corr_wx": "l_min_xp2_from_corr_wx",
    "l_enable_relaxed_clipping": "l_enable_relaxed_clipping",
    "l_use_c7_richardson": "l_use_C7_Richardson",
    "l_use_c11_richardson": "l_use_C11_Richardson",
    "l_diag_lscale_from_tau": "l_diag_Lscale_from_tau",
    "l_use_precip_frac": "l_use_precip_frac",
    "l_calc_thlp2_rad": "l_calc_thlp2_rad",
}
_PARAMETER_NAME_ALIASES = {name.lower(): name for name in PARAMETER_NAMES}


def canonical_flag_name(name: str) -> str:
    """Return the canonical spelling of a known case-insensitive flag name."""
    text = str(name).strip()
    return _FLAG_NAME_ALIASES.get(text.lower(), text)


def canonical_parameter_name(name: str) -> str:
    """Return the canonical spelling of a known case-insensitive parameter name."""
    text = str(name).strip()
    return _PARAMETER_NAME_ALIASES.get(text.lower(), text)


def _mark_inactive(
    states: dict[str, tuple[str, str, str | None]],
    names: Sequence[str],
    reason: str,
    rule_id: str,
) -> None:
    """Apply an inactive state without overriding an explicit unused state."""
    for name in names:
        if states[name][0] != "unused":
            states[name] = ("inactive-mode", reason, rule_id)


def parameter_activity(
    pdf_type: int | None = None,
    flags: Mapping[str, Any] | None = None,
) -> dict[str, tuple[str, str, str | None]]:
    """Resolve parameter activity from the maintained rule list.

    ``pdf_type`` is retained as a convenience for existing callers.  Supplying
    ``flags`` enables the finer branch rules; missing flags intentionally do
    not infer that an optional branch is disabled.  That makes a partial edit
    conservative while a resolved Run/Tune configuration is precise.
    """
    normalized = _canonicalize_setting_names(flags, _FLAG_NAME_ALIASES)
    if pdf_type is None and "iiPDF_type" in normalized:
        try:
            pdf_type = int(normalized["iiPDF_type"])
        except (TypeError, ValueError):
            pdf_type = None

    states = {
        name: ("active", "Used by its primary CLUBB owner.", None)
        for name in PARAMETER_NAMES
    }
    for name in UNUSED_PARAMETERS:
        states[name] = (
            "unused",
            "Retained in the tunable namelist but currently not read by a model tendency or PDF routine.",
            "PARAM-UNUSED",
        )

    if pdf_type is not None:
        for owner, names in PDF_PARAMETER_OWNERS.items():
            if pdf_type != owner:
                _mark_inactive(
                    states,
                    names,
                    f"Only used by PDF type {owner}; current PDF type is {pdf_type}.",
                    f"PDF{owner}-OWNERS",
                )
            else:
                for name in names:
                    if states[name][0] != "unused":
                        states[name] = ("active", f"Used by PDF type {owner}.", f"PDF{owner}-OWNERS")

    for flag_name, names in RICHARDSON_REPLACED_PARAMETERS.items():
        if _as_bool(normalized.get(flag_name, False)):
            _mark_inactive(
                states,
                names,
                f"{flag_name}=true replaces this tunable closure function with its Richardson formulation.",
                "RICHARDSON-CLOSURE",
            )

    if _as_bool(normalized.get("l_diag_Lscale_from_tau", False)):
        _mark_inactive(
            states,
            DIRECT_LSCALE_ONLY_PARAMETERS,
            "l_diag_Lscale_from_tau=true uses inverse-tau Lscale instead of this direct-Lscale control.",
            "TAU-LSCALE-MODE",
        )
    elif "l_diag_Lscale_from_tau" in normalized:
        _mark_inactive(
            states,
            INVERSE_TAU_PARAMETERS,
            "l_diag_Lscale_from_tau=false bypasses inverse-tau Lscale coefficients.",
            "DIRECT-LSCALE-MODE",
        )

    for flag_name, names in BRANCH_ONLY_PARAMETERS.items():
        if flag_name in normalized and not _as_bool(normalized[flag_name]):
            _mark_inactive(
                states,
                names,
                f"{flag_name}=false bypasses this parameterized branch.",
                "BRANCH-ONLY",
            )
    return states


def is_independently_tunable(state: Mapping[str, Any] | None) -> bool:
    """Return whether a semantic state may become an independent Tune axis."""
    return str((state or {}).get("state") or "active") == "active"


def _canonicalize_setting_names(values: Mapping[str, Any] | None, aliases: Mapping[str, str]) -> dict[str, Any]:
    """Accept Fortran's case-insensitive identifiers, return canonical names."""
    normalized: dict[str, Any] = {}
    for raw_name, value in dict(values or {}).items():
        name = str(raw_name).strip()
        normalized[aliases.get(name.lower(), name)] = value
    return normalized


def linked_parameter_groups(flags: Mapping[str, Any] | None = None) -> tuple[tuple[str, ...], ...]:
    """Return the equal-value parameter groups required by these flags.

    This is intentionally a small declarative view of the same rules used by
    :func:`resolve_clubb_settings`.  Consumers such as Dash and the tuner can
    present one logical control for a group without re-encoding a scientific
    compatibility rule in their own UI code.
    """

    normalized_flags = _canonicalize_setting_names(flags, _FLAG_NAME_ALIASES)
    groups: list[tuple[str, ...]] = [tuple(pair) for pair in LINKED_PARAMETER_GROUPS]
    if _as_bool(normalized_flags.get("l_damp_wp2_using_em", False)):
        groups.append(("C1", "C14"))
    return tuple(groups)


def linked_parameter_members(group_key: str) -> tuple[str, ...]:
    """Decode a serialized linked-group key produced by presentation code."""
    return tuple(part for part in str(group_key or "").split("=") if part)


def apply_linked_parameter_values(
    values: Mapping[str, Any] | None,
    linked_values: Mapping[str, Any] | None,
    *,
    file_name: str = "tunable",
) -> dict[str, Any]:
    """Expand logical equality controls into physical UI-neutral setting keys.

    ``values`` uses :func:`setting_key` keys.  ``linked_values`` maps a
    serialized group such as ``"C6rt=C6thl"`` to its single logical value.
    Keeping this transformation here means launch callers do not need their
    own interpretation of CLUBB's equality declarations.
    """
    expanded = dict(values or {})
    for group_key, value in dict(linked_values or {}).items():
        for name in linked_parameter_members(group_key):
            expanded[setting_key(file_name, canonical_parameter_name(name))] = value
    return expanded


def flag_relationship_groups() -> tuple[dict[str, Any], ...]:
    """Return UI-safe declarations for boolean compatibility groups."""
    return tuple(
        {
            **relation,
            "members": list(relation["members"]),
            "required_values": dict(relation.get("required_values") or {}),
            "switch_on_true": {
                str(name): dict(values)
                for name, values in dict(relation.get("switch_on_true") or {}).items()
            },
        }
        for relation in _FLAG_RELATIONSHIPS
    )


def resolve_clubb_settings(
    flags: Mapping[str, Any] | None = None,
    parameters: Mapping[str, Any] | None = None,
    *,
    auto_correct: bool = False,
) -> SettingsResolution:
    """Resolve forced settings, restrictions, and known parameter activity.

    Callers may supply a partial set.  Flag incompatibilities are reported,
    never silently rewritten.  The Run UI applies only declared checkbox-pair
    switch rules itself; arbitrary callers receive the same explicit conflict.
    ``auto_correct`` remains accepted for API compatibility but does not alter
    flags.  Feed effective values to ``validate_clubb_settings`` for the
    normal bound and complete-consistency checks.
    """
    supplied_flags = _canonicalize_setting_names(flags, _FLAG_NAME_ALIASES)
    supplied_parameters = _canonicalize_setting_names(parameters, _PARAMETER_NAME_ALIASES)
    effective_flags = dict(supplied_flags)
    forced_flags: dict[str, Any] = {}
    forced_parameters: dict[str, Any] = {}
    allowed: dict[str, list[Any]] = {}
    issues: list[ValidationIssue] = []

    pdf_type: int | None
    try:
        pdf_type = int(effective_flags.get("iiPDF_type")) if "iiPDF_type" in effective_flags else None
    except (TypeError, ValueError):
        pdf_type = None
        issues.append(ValidationIssue("invalid_integer_flag", "iiPDF_type must be an integer.", setting="iiPDF_type"))

    # numerical_check.F90 permits implicit wind-flux prediction only for
    # ADG1 (1) and the new hybrid PDF (7).
    if _as_bool(effective_flags.get("l_predict_upwp_vpwp", False)) and pdf_type is None:
        allowed["iiPDF_type"] = [1, 7]

    coupled: list[dict[str, Any]] = []
    for members in linked_parameter_groups(effective_flags):
        left, right = members
        supplied = [name for name in (left, right) if name in supplied_parameters]
        if len(supplied) == 1:
            master = supplied[0]
            follower = right if master == left else left
            forced_parameters[follower] = supplied_parameters[master]
            coupled.append({"members": [left, right], "master": master, "follower": follower, "relation": "equal"})
        elif len(supplied) == 2:
            coupled.append({"members": [left, right], "master": None, "follower": None, "relation": "equal"})
        else:
            coupled.append({"members": [left, right], "master": None, "follower": None, "relation": "equal"})

    flag_states: dict[str, dict[str, str]] = {}

    parameter_states = {
        name: {"state": state, "reason": reason, "rule_id": rule_id}
        for name, (state, reason, rule_id) in parameter_activity(pdf_type, effective_flags).items()
    }
    for name, value in forced_parameters.items():
        parameter_states[name] = {"state": "linked", "reason": f"Must equal its linked partner ({value!r}).", "rule_id": "PARAM-C6-EQUALITY"}

    validation = validate_clubb_settings(supplied_parameters, effective_flags)
    issues.extend(validation)
    conflict_keys: set[str] = set()
    for issue in validation:
        if issue.severity != "error":
            continue
        if issue.code == "predict_wind_pdf_type":
            conflict_keys.update({"flags:iiPDF_type", "flag:l_predict_upwp_vpwp"})
        elif issue.code == "em_wp2_damping_inconsistent":
            conflict_keys.update({"flag:l_damp_wp2_using_em", "flag:l_stability_correct_tau_zm"})
        elif issue.code == "incompatible_xp2_clipping":
            conflict_keys.update({"flag:l_min_xp2_from_corr_wx", "flag:l_enable_relaxed_clipping"})
        elif issue.setting:
            conflict_keys.add(f"flag:{issue.setting}")
    return SettingsResolution(
        effective_flags=effective_flags,
        forced_flags=forced_flags,
        forced_parameters=forced_parameters,
        allowed_flag_values=allowed,
        conflicting_setting_keys=sorted(conflict_keys),
        flag_states=flag_states,
        flag_relationships=list(flag_relationship_groups()),
        parameter_states=parameter_states,
        coupled_parameters=coupled,
        issues=tuple(issues),
    )


def linked_flag_updates(flag_name: str, enabled: bool) -> dict[str, bool]:
    """Return only declared checkbox-pair changes for one user action.

    This is intentionally separate from validation: arbitrary incompatible
    flags remain selected and are reported as errors, while explicit UI groups
    may synchronize their members' checked state.
    """
    if not enabled:
        return {}
    name = canonical_flag_name(flag_name)
    for relationship in flag_relationship_groups():
        members = [str(member) for member in relationship.get("members") or []]
        if name not in members:
            continue
        if relationship.get("relation") == "opposite":
            return {member: False for member in members if member != name}
        return {
            str(target): bool(value)
            for target, value in dict((relationship.get("switch_on_true") or {}).get(name) or {}).items()
        }
    return {}


def _control_state(
    *,
    changed: bool,
    invalid: bool,
    inactive: bool = False,
    linked: bool = False,
    reason: str = "",
) -> dict[str, Any]:
    """Return one semantic state with explicit precedence for renderers."""
    if invalid:
        state = "invalid"
    elif inactive:
        state = "inactive"
    elif linked:
        state = "linked-changed" if changed else "linked"
    else:
        state = "changed" if changed else "default"
    return {"state": state, "changed": bool(changed), "reason": reason}


def evaluate_settings(
    schema: Mapping[str, Any],
    *,
    flag_values: Mapping[str, Any] | None = None,
    parameter_values: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Evaluate all settings semantics from one config schema and value maps.

    The result is intentionally presentation-ready but UI-neutral: callers
    receive normalized values, overrides, relationships, validation findings,
    and one semantic state per physical setting.  Dash should render this data
    rather than comparing defaults or interpreting CLUBB relationships itself.
    """
    schema_flags = _canonicalize_setting_names(schema.get("flag_defaults"), _FLAG_NAME_ALIASES)
    schema_parameters = {
        str(file_name): _canonicalize_setting_names(values, _PARAMETER_NAME_ALIASES)
        for file_name, values in dict(schema.get("parameter_defaults") or {}).items()
    }
    metadata = [item for item in schema.get("parameter_metadata") or [] if isinstance(item, Mapping)]
    supplied_flags = _canonicalize_setting_names(flag_values, _FLAG_NAME_ALIASES)
    supplied_values = dict(parameter_values or {})

    current_flags = dict(schema_flags)
    current_flags.update(supplied_flags)
    current_parameters = {file_name: dict(values) for file_name, values in schema_parameters.items()}
    for item in metadata:
        file_name = str(item.get("file") or "")
        name = canonical_parameter_name(str(item.get("name") or ""))
        key = setting_key(file_name, name)
        if key in supplied_values:
            current_parameters.setdefault(file_name, {})[name] = supplied_values[key]
    current_flags.update(current_parameters.get("flags", {}))

    all_tunable = dict(current_parameters.get("tunable", {}))
    resolution = resolve_clubb_settings(current_flags, all_tunable)
    resolution_data = resolution.as_dict()
    invalid_keys = set(resolution.conflicting_setting_keys)
    for issue in resolution.issues:
        if issue.severity != "error" or not issue.setting:
            continue
        parameter_name = canonical_parameter_name(issue.setting)
        flag_name = canonical_flag_name(issue.setting)
        if flag_name in schema_flags:
            invalid_keys.add(setting_key("flag", flag_name))
        for item in metadata:
            file_name = str(item.get("file") or "")
            name = canonical_parameter_name(str(item.get("name") or ""))
            if name == parameter_name:
                invalid_keys.add(setting_key(file_name, name))
    control_states: dict[str, dict[str, Any]] = {}
    overrides = {"flags": {}, "tunable": {}, "silhs": {}}

    for name, default in schema_flags.items():
        current = current_flags.get(name, default)
        changed = not settings_values_equal(current, default, boolean=True)
        key = setting_key("flag", name)
        control_states[key] = _control_state(
            changed=changed,
            invalid=key in invalid_keys,
            reason=next((issue.message for issue in resolution.issues if issue.setting == name), ""),
        )
        if changed:
            overrides["flags"][name] = format_setting_value(_as_bool(current))

    active_groups = [list(group) for group in linked_parameter_groups(current_flags)]
    linked_members = {member for group in active_groups for member in group}
    linked_group_states: dict[str, dict[str, Any]] = {}
    for item in metadata:
        file_name = str(item.get("file") or "")
        name = canonical_parameter_name(str(item.get("name") or ""))
        if not file_name or not name:
            continue
        key = setting_key(file_name, name)
        default = schema_parameters.get(file_name, {}).get(name, "")
        current = current_parameters.get(file_name, {}).get(name, default)
        changed = not settings_values_equal(current, default)
        parameter_state = resolution.parameter_states.get(name, {}) if file_name == "tunable" else {}
        inactive = not is_independently_tunable(parameter_state)
        invalid = key in invalid_keys
        control_states[key] = _control_state(
            changed=changed,
            invalid=invalid,
            inactive=inactive,
            linked=file_name == "tunable" and name in linked_members,
            reason=str(parameter_state.get("reason") or ""),
        )
        if changed:
            overrides.setdefault(file_name, {})[name] = format_setting_value(current)

    for group in active_groups:
        group_key = "=".join(group)
        member_keys = [setting_key("tunable", name) for name in group]
        member_states = [control_states.get(key, {}) for key in member_keys]
        linked_group_states[group_key] = _control_state(
            changed=any(bool(item.get("changed")) for item in member_states),
            invalid=any(item.get("state") == "invalid" for item in member_states),
            linked=True,
            reason="These parameters are constrained to the same value by CLUBB.",
        )

    return {
        **resolution_data,
        "normalized_flags": current_flags,
        "normalized_parameters": current_parameters,
        "overrides": overrides,
        "control_states": control_states,
        "linked_parameter_groups": active_groups,
        "linked_group_states": linked_group_states,
        "changed_setting_keys": sorted(key for key, state in control_states.items() if state.get("changed")),
        "tunable_default_ranges": dict(schema.get("tunable_default_ranges") or {}),
    }


def _as_bool(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"t", ".t.", "true", ".true.", "1", "yes"}
    return bool(value)


def _finite_values(value: Any) -> list[float]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        value = [value]
    values: list[float] = []
    for item in value:
        number = float(str(item).replace("D", "E").replace("d", "e"))
        if not math.isfinite(number):
            raise ValueError("must be finite")
        values.append(number)
    return values


def _fortran_equal(left: float, right: float) -> bool:
    """Use the relative equality test written in parameters_tunable.F90."""
    return abs(left - right) <= abs(left + right) * 0.5 * sys.float_info.epsilon


def validate_clubb_settings(
    parameters: Mapping[str, Any],
    flags: Mapping[str, Any] | None = None,
    *,
    lmin: float = 1.0,
    l_implemented: bool = False,
    l_input_fields: bool = False,
    explicit_turbulent_adv_wpxp: bool = False,
    advance_orders: Mapping[str, int] | None = None,
    sponge_flags: Mapping[str, Any] | None = None,
    is_gpu: bool = False,
    l_avg_Lscale: bool | None = None,
    rad_scheme: str | None = None,
    microphys_scheme: str | None = None,
) -> list[ValidationIssue]:
    """Validate a complete effective CLUBB settings set without F2PY.

    Missing parameters are ignored so callers may validate a partial edit.  A
    production caller should pass its resolved namelist values.  Missing flags
    use the relevant Fortran defaults (false, PDF type 1, and option 1).
    """
    issues: list[ValidationIssue] = []
    values: dict[str, list[float]] = {}
    for name, raw in parameters.items():
        if name not in PARAMETER_NAMES:
            issues.append(ValidationIssue("unknown_parameter", f"Unknown tunable parameter {name!r}.", setting=name))
            continue
        try:
            values[name] = _finite_values(raw)
        except (TypeError, ValueError) as exc:
            issues.append(ValidationIssue("invalid_parameter_value", f"{name} {exc}.", setting=name))

    try:
        parsed_lmin = float(lmin)
        if not math.isfinite(parsed_lmin):
            raise ValueError("must be finite")
        if parsed_lmin < 1.0:
            issues.append(ValidationIssue("lmin_too_small", "lmin must be at least 1 m.", setting="lmin"))
    except (TypeError, ValueError):
        issues.append(ValidationIssue("invalid_lmin", "lmin must be numeric.", setting="lmin"))

    for name, series in values.items():
        low, high = PARAMETER_HARD_BOUNDS[PARAMETER_NAMES.index(name)]
        for column, value in enumerate(series):
            if low is not None and value < low:
                issues.append(ValidationIssue("parameter_below_min", f"{name}={value:g} is below the allowed minimum {low:g}.", setting=name, column=column))
            if high is not None and value > high:
                issues.append(ValidationIssue("parameter_above_max", f"{name}={value:g} exceeds the allowed maximum {high:g}.", setting=name, column=column))

    for left, right in LINKED_PARAMETER_GROUPS:
        if left not in values or right not in values:
            continue
        left_values, right_values = values[left], values[right]
        count = max(len(left_values), len(right_values))
        if len(left_values) not in {1, count} or len(right_values) not in {1, count}:
            issues.append(ValidationIssue("parameter_column_shape", f"{left} and {right} must have matching scalar or column counts.", setting=left))
            continue
        for column in range(count):
            a, b = left_values[min(column, len(left_values) - 1)], right_values[min(column, len(right_values) - 1)]
            if not _fortran_equal(a, b):
                issues.append(ValidationIssue("required_equal_parameters", f"{left} and {right} must be equal.", setting=left, column=column))

    flag_values = dict(flags or {})
    def flag(name: str, default: Any = False) -> Any:
        return flag_values.get(name, default)

    def integer_flag(name: str, default: int) -> int | None:
        try:
            return int(flag(name, default))
        except (TypeError, ValueError):
            issues.append(ValidationIssue("invalid_integer_flag", f"{name} must be an integer.", setting=name))
            return None

    pdf_type = integer_flag("iiPDF_type", 1)
    if pdf_type is None:
        pass
    elif pdf_type < 1 or pdf_type > 7:
        issues.append(ValidationIssue("invalid_pdf_type", "iiPDF_type must be between 1 and 7.", setting="iiPDF_type"))
    elif pdf_type in {2, 3, 4, 5, 6} and not l_input_fields:
        issues.append(ValidationIssue("pdf_requires_input_fields", f"iiPDF_type={pdf_type} requires l_input_fields=true.", setting="iiPDF_type"))
    saturation_formula = integer_flag("saturation_formula", 1)
    if saturation_formula is not None and not 1 <= saturation_formula <= 4:
        issues.append(ValidationIssue("invalid_saturation_formula", "saturation_formula must be between 1 and 4.", setting="saturation_formula"))
    call_placement = integer_flag("ipdf_call_placement", 1)
    if call_placement is not None and not 1 <= call_placement <= 2:
        issues.append(ValidationIssue("invalid_pdf_call_placement", "ipdf_call_placement must be 1 or 2.", setting="ipdf_call_placement"))
    if is_gpu and l_input_fields:
        issues.append(ValidationIssue("gpu_input_fields", "l_input_fields=true is not usable in a GPU build.", setting="l_input_fields"))

    if _as_bool(flag("l_damp_wp2_using_em")):
        if _as_bool(flag("l_stability_correct_tau_zm")) or (
            "C1" in values and "C14" in values and any(
                not _fortran_equal(values["C1"][min(i, len(values["C1"]) - 1)], values["C14"][min(i, len(values["C14"]) - 1)])
                for i in range(max(len(values["C1"]), len(values["C14"])))
            )
        ):
            issues.append(ValidationIssue("em_wp2_damping_inconsistent", "l_damp_wp2_using_em requires C1=C14 and l_stability_correct_tau_zm=false.", setting="l_damp_wp2_using_em"))
    if _as_bool(flag("l_predict_upwp_vpwp")):
        if explicit_turbulent_adv_wpxp:
            issues.append(ValidationIssue("predict_wind_requires_implicit_ta", "l_predict_upwp_vpwp is incompatible with explicit wpxp turbulent advection.", setting="l_predict_upwp_vpwp"))
        if pdf_type is not None and pdf_type not in {1, 7}:
            issues.append(ValidationIssue("predict_wind_pdf_type", "l_predict_upwp_vpwp requires PDF type 1 or 7.", setting="l_predict_upwp_vpwp"))
    if _as_bool(flag("l_min_xp2_from_corr_wx")) and _as_bool(flag("l_enable_relaxed_clipping")):
        issues.append(ValidationIssue("incompatible_xp2_clipping", "l_min_xp2_from_corr_wx and l_enable_relaxed_clipping cannot both be true.", setting="l_min_xp2_from_corr_wx"))
    elif not _as_bool(flag("l_min_xp2_from_corr_wx")) and not _as_bool(flag("l_enable_relaxed_clipping")):
        issues.append(ValidationIssue("xp2_clipping_warning", "Both xp2 clipping choices are false; Fortran currently warns but does not reject this.", severity="warning", setting="l_min_xp2_from_corr_wx"))

    if _as_bool(flag("l_Lscale_plume_centered")) and l_avg_Lscale is False:
        issues.append(ValidationIssue("lscale_plume_centered_requires_avg", "l_Lscale_plume_centered requires l_avg_Lscale=true.", setting="l_Lscale_plume_centered"))

    radiation_context_known = rad_scheme is not None or "rad_scheme" in flag_values
    active_rad_scheme = str(rad_scheme if rad_scheme is not None else flag("rad_scheme", "none")).strip().lower()
    if radiation_context_known and _as_bool(flag("l_calc_thlp2_rad")):
        if _as_bool(flag("l_silhs_rad")):
            issues.append(ValidationIssue("thlp2_rad_silhs_rad_conflict", "l_calc_thlp2_rad and l_silhs_rad cannot both be true.", setting="l_calc_thlp2_rad"))
        if active_rad_scheme == "none":
            issues.append(ValidationIssue("thlp2_rad_requires_radiation", "l_calc_thlp2_rad requires rad_scheme other than 'none'.", setting="l_calc_thlp2_rad"))
    if radiation_context_known and _as_bool(flag("l_soil_veg")) and active_rad_scheme == "none":
        issues.append(ValidationIssue("soil_veg_requires_radiation", "l_soil_veg requires rad_scheme other than 'none'.", setting="l_soil_veg"))

    # These two established microphysics branches have opposite Nc ownership.
    # The caller supplies the scheme because it is not a ConfigFlags member.
    normalized_microphysics = str(microphys_scheme or "").strip().lower()
    if normalized_microphysics == "coamps" and not _as_bool(flag("l_predict_Nc")):
        issues.append(ValidationIssue("coamps_requires_predict_nc", "COAMPS microphysics requires l_predict_Nc=true.", setting="l_predict_Nc"))
    if normalized_microphysics in {"kk", "khairoutdinov_kogan"} and _as_bool(flag("l_predict_Nc")):
        issues.append(ValidationIssue("kk_requires_diagnosed_nc", "KK microphysics requires l_predict_Nc=false.", setting="l_predict_Nc"))

    orders = {"order_xm_wpxp": 1, "order_wp2_wp3": 2, "order_xp2_xpyp": 3, "order_windm": 4}
    orders.update(advance_orders or {})
    active_orders: dict[int, str] = {}
    for name, raw_order in orders.items():
        try:
            order = int(raw_order)
        except (TypeError, ValueError):
            issues.append(ValidationIssue("invalid_advance_order", f"{name} must be -1 or between 1 and 4.", setting=name))
            continue
        if order != -1 and order not in {1, 2, 3, 4}:
            issues.append(ValidationIssue("invalid_advance_order", f"{name} must be -1 or between 1 and 4.", setting=name))
        elif order != -1 and order in active_orders:
            issues.append(ValidationIssue("duplicate_advance_order", f"{name} and {active_orders[order]} use order {order}; active advance orders must be unique.", setting=name))
        elif order != -1:
            active_orders[order] = name

    if _as_bool(flag("l_diag_Lscale_from_tau")):
        for name in ("C1", "C1b", "C2rt", "C2thl", "C2rtthl", "C6rt", "C6rtb", "C6thl", "C6thlb", "C14"):
            if name in values and any(value != 1.0 for value in values[name]):
                issues.append(ValidationIssue("tau_lscale_parameter_warning", f"{name} should be 1 when l_diag_Lscale_from_tau is true; Fortran currently warns only.", severity="warning", setting=name))
    if l_implemented:
        for name in ("l_rtm_nudge", "l_uv_nudge"):
            if _as_bool(flag(name)):
                issues.append(ValidationIssue("host_incompatible_flag", f"{name} must be false when l_implemented is true.", setting=name))
        for name, enabled in (sponge_flags or {}).items():
            if _as_bool(enabled):
                issues.append(ValidationIssue("host_incompatible_sponge", f"{name} sponge damping must be false when l_implemented is true.", setting=name))
    return issues


def validation_errors(*args: Any, **kwargs: Any) -> list[ValidationIssue]:
    """Convenience wrapper returning only rejecting findings."""
    return [issue for issue in validate_clubb_settings(*args, **kwargs) if issue.severity == "error"]
