"""CLUBB case initialization and cleanup utilities."""
import math
from pathlib import Path

import numpy as np
import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.config_flags import ConfigFlags
from clubb_jax.src.CLUBB_core.constants_clubb import eps
from clubb_jax.src.CLUBB_core.grid_class import setup_grid as py_setup_grid
from clubb_jax.src.CLUBB_core.sclr_idx import SclrIdx
from clubb_jax.src.CLUBB_core.pdf_params import (
    init_pdf_implicit_coefs_terms_api,
    init_pdf_params as init_pdf_params_py,
)
from clubb_jax.src.CLUBB_core.err_info import ErrInfo
from clubb_jax.src.CLUBB_core.calc_pressure import calculate_thvm
from clubb_jax.src.CLUBB_core.error_code import (
    clubb_at_least_debug_level,
    set_debug_level as set_jax_debug_level,
)
from clubb_jax.src.CLUBB_core.grid_class import zt2zm
from clubb_jax.src.CLUBB_core.model_flags import get_default_config_flags
from clubb_jax.src.CLUBB_core.numerical_check import check_clubb_settings
from clubb_jax.src.CLUBB_core.parameters_tunable import (
    calc_derived_params,
    get_param_names,
    init_clubb_params,
)
from clubb_jax.src.CLUBB_core.saturation import rcm_sat_adj, sat_mixrat_liq
from clubb_jax.src.Input_fields.hydrostatic_module import hydrostatic

# I/O
from clubb_jax.src.Input_fields.grid_file import read_grid_file
from clubb_jax.src.CLUBB_core.stats_netcdf import StatsWriter
from clubb_jax.src.Input_fields.namelist import read_namelist
from clubb_jax.src.Input_fields.sounding import (
    read_sounding,
    interpolate_sounding,
    read_scalar_sounding,
    interpolate_scalar_sounding,
)
from clubb_jax.src.Input_fields.surface import read_surface
from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.mpace_a import mpace_a_init
from clubb_jax.src.Radiation.soil_vegetation import initialize_soil_veg
from clubb_jax.src.Radiation.parameters_radiation import initialize_radiation_parameters


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


# ── Physical constants (from constants_clubb.F90, standalone block) ──────
Cp = 1004.67
Lv = 2.5e6
Rd = 287.04
Rv = 461.5
ep = Rd / Rv                # 0.621993...  (must match Fortran's Rd/Rv exactly)
ep1 = (1.0 - ep) / ep      # ~0.608
ep2 = 1.0 / ep              # ~1.608
kappa = Rd / Cp
grav = 9.81
p0 = 1.0e5
omega_planet = 7.292e-5
radians_per_deg = math.pi / 180.0
rt_tol = 1.0e-8
thl_tol = 1.0e-2
w_tol = 2.0e-2
em_min = 1.5 * w_tol**2
cloud_frac_min = 0.005
Nc0_in_cloud = 100.0e6      # [num/m^3]

_CLOUD_FEEDBACK_CASES = {
    "cloud_feedback_s6",
    "cloud_feedback_s6_p2k",
    "cloud_feedback_s11",
    "cloud_feedback_s11_p2k",
    "cloud_feedback_s12",
    "cloud_feedback_s12_p2k",
}


def _initialize_em_profile(runtype: str, gr, um: np.ndarray):
    """Mirror initialize_clubb() case-based em setup from Fortran."""
    runtype = str(runtype).strip()
    zm = gr.zm
    ngrdcol, nzm = zm.shape
    em = np.full((ngrdcol, nzm), em_min, dtype=np.float64)
    um_out = um.copy()

    def _set_cloud_top_profile(cloud_top: float, em_max_val: float):
        e = np.where(zm < cloud_top, em_max_val, em_min)
        if nzm > 1:
            e[:, 0] = e[:, 1]
        e[:, -1] = em_min
        return e

    em_min_cases = {"bomex", "ekman", "atex_long", "arm"}
    em_one_topmin_cases = {
        "generic", "arm_97", "twp_ice", "arm_0003", "arm_3year",
        "dycoms2_rf02", "gabls3",
    } | _CLOUD_FEEDBACK_CASES
    em_point1_topmin_cases = {"lba", "cobra"}
    fixed_cloud_top_cases = {
        "astex_a209": (700.0, 1.0),
        "fire": (700.0, 4.5),
        "dycoms2_rf01": (800.0, 1.1),
        "mpace_b": (1300.0, 1.0),
        "rico": (1500.0, 1.0),
    }
    offset_cloud_top_cases = {
        "nov11_altocu": 2800.0,
        "clex9_nov02": 2200.0,
        "clex9_oct14": 3500.0,
    }

    if runtype in em_one_topmin_cases:
        em[:, :] = 1.0
        em[:, -1] = em_min
    elif runtype in em_min_cases:
        em[:, :] = em_min
    elif runtype == "atex":
        um_out = np.maximum(um_out, -8.0)
        em[:, :] = em_min
    elif runtype in fixed_cloud_top_cases:
        cloud_top, em_max = fixed_cloud_top_cases[runtype]
        em = _set_cloud_top_profile(cloud_top, em_max)
    elif runtype in offset_cloud_top_cases:
        em = _set_cloud_top_profile(offset_cloud_top_cases[runtype] + float(zm[0, 0]), 0.01)
    elif runtype == "jun25_altocu":
        em[:, :] = 0.01
        if nzm > 1:
            em[:, 0] = em[:, 1]
        em[:, -1] = em_min
    elif runtype in em_point1_topmin_cases:
        em[:, :] = 0.1
        em[:, -1] = em_min
    elif runtype == "gabls2":
        cloud_top = 800.0
        em = np.where(zm < cloud_top, 0.5 * (1.0 - (zm / cloud_top)), em_min)
        if nzm > 1:
            em[:, 0] = em[:, 1]
        em[:, -1] = em_min
    elif runtype == "gabls3_night":
        em[:, :] = 1.0
    elif runtype == "coriolis_test":
        depth = (zm[:, -1] - zm[:, 0])[:, None]
        em = np.sin(np.pi * zm / depth) * (w_tol**2) * 6.0

    return em, um_out


def _initialize_turbulence_state(runtype: str, gr, dt_main: float,
                                 fcor_y: np.ndarray, um: np.ndarray):
    """Mirror initialize_clubb() em/wp2/up2/vp2/upwp initialization."""
    em, um_adj = _initialize_em_profile(runtype, gr, um)

    wp2 = (2.0 / 3.0) * em
    up2 = (2.0 / 3.0) * em
    vp2 = (2.0 / 3.0) * em
    upwp = np.zeros_like(em)

    if str(runtype).strip() == "coriolis_test":
        w_tol_sqd = w_tol**2
        wp2 = (1.0 / 3.0) * em + w_tol_sqd
        up2 = (3.0 / 3.0) * em + w_tol_sqd
        vp2 = (2.0 / 3.0) * em + w_tol_sqd
        em = em + 1.5 * w_tol_sqd
        upwp = 0.5 * dt_main * fcor_y[:, None] * (up2 - wp2)

    return em, wp2, up2, vp2, upwp, um_adj


def run_clubb(namelist_path: str, l_stdout: bool = True):
    """Run CLUBB standalone for a case described by a namelist file.

    Args:
        namelist_path: path to *_model.in file
        l_stdout: print timestep info to stdout
    """
    from clubb_jax.src.advance_clubb_to_end import advance_clubb_to_end

    state = init_clubb_case(namelist_path)
    try:
        num_batches = state['total_param_sets'] // state['ngrdcol']
        for batch_num in range(1, num_batches + 1):
            # Reset fields and select the corresponding parameter/output slice.
            set_case_initial_conditions(state, batch_num=batch_num)
            advance_clubb_to_end(state, l_stdout=l_stdout)
    finally:
        clean_up_clubb(state)
    return state


def _resolve_stats_registry_path(namelist_path: str, cfg: dict) -> Path:
    """Resolve the stats registry file path.

    Priority:
      1) explicit namelist key `stats_registry`
      2) the runfile itself, if it contains `&clubb_stats_nl`
      3) repository default `input/stats/standard_stats.in`
    """
    configured = str(cfg.get('stats_registry', '')).strip()
    if configured:
        p = Path(configured)
        if not p.is_absolute():
            p = Path(namelist_path).resolve().parent / p
        return p.resolve()

    runfile = Path(namelist_path).resolve()
    if '&clubb_stats_nl' in runfile.read_text().lower():
        return runfile

    return _repo_root() / "input" / "stats" / "standard_stats.in"


def _resolve_case_input_path(namelist_dir: Path, runtype: str, suffix: str) -> Path:
    """Resolve case input files for either model.in or aggregated CASE.in runs."""
    candidate = namelist_dir / f"{runtype}{suffix}"
    if candidate.exists():
        return candidate

    fallback = _repo_root() / "input" / "case_setups" / f"{runtype}{suffix}"
    if fallback.exists():
        return fallback

    raise FileNotFoundError(
        f"Required case input file not found: {candidate} (also checked {fallback})"
    )


def _clean_namelist_path(path_value) -> str:
    """Normalize namelist path strings (strip quotes and whitespace)."""
    return str(path_value).strip().strip("'\"")


def _validate_scalar_column_names(names, idx_rt: int, idx_thl: int, idx_co2: int, label: str):
    """Validate scalar column order against namelist scalar-index mapping."""
    for col_idx, name in enumerate(names, start=1):
        if name == 'CO2[ppmv]' and idx_co2 > 0 and col_idx != idx_co2:
            raise ValueError(f"{label}: iisclr/iiedsclr_CO2 index does not match column order.")
        if name == 'rt[kg/kg]' and idx_rt > 0 and col_idx != idx_rt:
            raise ValueError(f"{label}: iisclr/iiedsclr_rt index does not match column order.")
        if name in {'thm[K]', 'thlm[K]', 'T[K]'} and idx_thl > 0 and col_idx != idx_thl:
            raise ValueError(f"{label}: iisclr/iiedsclr_thl index does not match column order.")


def _resolve_grid_file_path(namelist_dir: Path, grid_path_value) -> Path:
    """Resolve a grid filename from namelist conventions to an existing path."""
    raw = _clean_namelist_path(grid_path_value)
    if not raw:
        raise ValueError("Grid filename is empty.")

    p = Path(raw)
    repo_root = _repo_root()
    candidates = []
    if p.is_absolute():
        candidates.append(p)
    else:
        candidates.append(Path.cwd() / p)
        candidates.append(namelist_dir / p)
        candidates.append(repo_root / p)
        if raw.startswith("../input/"):
            candidates.append(repo_root / raw[3:])

    seen = set()
    for candidate in candidates:
        key = str(candidate)
        if key in seen:
            continue
        seen.add(key)
        if candidate.exists():
            return candidate.resolve()

    raise FileNotFoundError(
        f"Grid file not found: {raw}. Checked: "
        + ", ".join(str(c) for c in candidates)
    )


# =========================================================================
# Feature gate
# =========================================================================

def _check_unsupported_features(cfg: dict, flags, microphys_scheme: str,
                                rad_scheme: str, l_calc_thlp2_rad: bool):
    """Check for namelist settings that the JAX driver does not support.

    Raises ValueError with a clear message listing all unsupported features
    that are enabled, so the user can fix them all at once.
    """
    errors = []

    # --- Microphysics ---
    if microphys_scheme not in {"none", "khairoutdinov_kogan", "morrison"}:
        errors.append(
            f"microphys_scheme = '{microphys_scheme}' is not supported "
            "(supported: none, khairoutdinov_kogan, morrison)."
        )

    for name in ('l_gfdl_activation',):
        if bool(cfg.get(name, False)):
            errors.append(f"{name} is not yet supported by the microphysics interface")

    # --- Radiation ---
    supported_rad = {"none", "simplified", "simplified_bomex", "lba"}
    if rad_scheme not in supported_rad:
        errors.append(
            f"rad_scheme = '{rad_scheme}' is not supported "
            f"(supported: {', '.join(sorted(supported_rad))})."
        )

    if l_calc_thlp2_rad and rad_scheme == "none":
        errors.append(
            "l_calc_thlp2_rad = true is incompatible with rad_scheme = 'none'."
        )

    # --- Sponge damping ---
    _SPONGE_FIELDS = ["thlm", "rtm", "uv", "wp2", "wp3", "up2_vp2"]
    sponge_enabled = [
        f for f in _SPONGE_FIELDS
        if bool(cfg.get(f'{f}_sponge_damp_settings%l_sponge_damping', False))
    ]
    if sponge_enabled:
        names = ", ".join(sponge_enabled)
        errors.append(
            f"Sponge damping is enabled for [{names}] but is not supported "
            "(sponge_damp routines are not called from the Python driver)."
        )

    # --- SILHS / Latin Hypercube sampling ---
    lh_type = str(cfg.get('lh_microphys_type', 'disabled')).strip().lower()
    if lh_type not in ('disabled', 'interactive', 'non-interactive'):
        errors.append(f"Unknown lh_microphys_type = '{lh_type}'")
    if bool(cfg.get('l_silhs_rad', False)):
        errors.append("l_silhs_rad = true is not supported (SILHS radiation is not ported).")

    # --- Input fields (time-dependent forcing from files) ---
    if bool(cfg.get('l_input_fields', False)):
        errors.append("l_input_fields = true is not supported.")

    # --- Generalized grid test ---
    if bool(cfg.get('l_test_grid_generalization', False)):
        errors.append("l_test_grid_generalization = true is not supported.")

    # --- Adaptive gridding ---
    # grid_adapt_in_time_method > 0 means some form of adaptation is active.
    grid_adapt = int(cfg.get('grid_adapt_in_time_method', 0))
    if grid_adapt > 0:
        errors.append(
            f"grid_adapt_in_time_method = {grid_adapt} is not supported "
            "(only 0 / no adaptation is implemented)."
        )

    if errors:
        msg = "Python driver does not support the following enabled features:\n"
        msg += "\n".join(f"  - {e}" for e in errors)
        raise ValueError(msg)


# =========================================================================
# Initialization
# =========================================================================

def init_clubb_case(namelist_path: str) -> dict:
    """Initialize a CLUBB case from a namelist file.

    Returns a dict containing all model state arrays and config.
    """
    cfg = read_namelist(namelist_path)
    namelist_dir = Path(namelist_path).resolve().parent

    # Unpack key config values
    # Resolve total parameter count and runtime batch size before allocation.
    total_param_sets = int(cfg['ngrdcol'])
    batch_size = int(cfg.get('batch_size', -1))
    if batch_size == -1:
        batch_size = total_param_sets
    if total_param_sets < 1:
        raise ValueError('ngrdcol in &multicol_def must be >= 1')
    if batch_size < 1 or batch_size > total_param_sets:
        raise ValueError('batch_size in &multicol_def must lie between 1 and ngrdcol')
    if total_param_sets % batch_size != 0:
        raise ValueError('ngrdcol in &multicol_def must be evenly divisible by batch_size')
    ngrdcol = batch_size
    if bool(cfg.get('l_restart', False)) and total_param_sets > batch_size:
        # TODO: select the active batch's saved columns in the restart reader.
        # It currently restores the whole saved column dimension in one read.
        raise ValueError('JAX restart does not yet support runtime batching')
    nzmax = cfg['nzmax']
    grid_type = cfg['grid_type']
    dt_main = cfg['dt_main']
    dt_rad = cfg['dt_rad']
    runtype = cfg['runtype']
    sclr_dim = cfg['sclr_dim']
    edsclr_dim = cfg['edsclr_dim']

    # ── 1. Initialize error handling ────────────────────────────────────
    set_jax_debug_level(cfg['debug_level'])
    err_info = ErrInfo.initialized(ngrdcol=ngrdcol)

    # ── 2. Get config flags ─────────────────────────────────────────────
    flags = get_default_config_flags()
    # Override from namelist (configurable_clubb_flags_nl)
    flag_overrides = {}
    for name in ConfigFlags._fields:
        if name.lower() in cfg:
            flag_overrides[name] = cfg[name.lower()]
    if flag_overrides:
        d = flags._asdict()
        d.update(flag_overrides)
        flags = ConfigFlags(**d)

    saturation_formula = flags.saturation_formula
    microphys_scheme = str(cfg.get('microphys_scheme', 'none')).strip().strip("'\"").lower()
    l_cloud_sed = bool(cfg.get('l_cloud_sed', False))
    sigma_g = float(cfg.get('sigma_g', 1.5))
    rad_scheme = str(cfg.get('rad_scheme', 'none')).strip().strip("'\"").lower()
    l_calc_thlp2_rad = bool(cfg.get('l_calc_thlp2_rad', flags.l_calc_thlp2_rad))

    _check_unsupported_features(cfg, flags, microphys_scheme, rad_scheme, l_calc_thlp2_rad)

    # ── 3. Read sounding ────────────────────────────────────────────────
    snd_path = _resolve_case_input_path(namelist_dir, runtype, "_sounding.in")
    snd = read_sounding(str(snd_path))
    temperature_type = snd['temperature_type']
    subs_type = snd['subs_type']

    # ── 4. Set up grid ──────────────────────────────────────────────────
    deltaz = np.full(ngrdcol, cfg['deltaz_nl'])
    zm_init = np.full(ngrdcol, cfg['zm_init_nl'])
    zm_top = np.full(ngrdcol, cfg['zm_top_nl'])
    sfc_elevation = np.full(ngrdcol, cfg['sfc_elevation_nl'])

    zt_grid_fname = _clean_namelist_path(cfg.get('zt_grid_fname', ''))
    zm_grid_fname = _clean_namelist_path(cfg.get('zm_grid_fname', ''))

    # For grid_type 1 (even spacing), heights are computed by setup_grid.
    # For stretched grids, these arrays are loaded from *.grd files.
    momentum_heights = None
    thermodynamic_heights = None

    if grid_type == 1:
        if zt_grid_fname or zm_grid_fname:
            raise ValueError(
                "grid_type=1 requires both zt_grid_fname and zm_grid_fname to be empty."
            )
    elif grid_type == 2:
        if zm_grid_fname:
            raise ValueError("grid_type=2 requires zm_grid_fname to be empty.")
        if not zt_grid_fname:
            raise ValueError("grid_type=2 requires zt_grid_fname.")

        zt_grid_path = _resolve_grid_file_path(namelist_dir, zt_grid_fname)
        zt_levels = read_grid_file(str(zt_grid_path))
        expected = nzmax - 1
        if zt_levels.size != expected:
            raise ValueError(
                f"zt grid file {zt_grid_path} has {zt_levels.size} levels; "
                f"expected nzmax-1={expected}."
            )
        thermodynamic_heights = np.tile(zt_levels[None, :], (ngrdcol, 1))
    elif grid_type == 3:
        if zt_grid_fname:
            raise ValueError("grid_type=3 requires zt_grid_fname to be empty.")
        if not zm_grid_fname:
            raise ValueError("grid_type=3 requires zm_grid_fname.")

        zm_grid_path = _resolve_grid_file_path(namelist_dir, zm_grid_fname)
        zm_levels = read_grid_file(str(zm_grid_path))
        expected = nzmax
        if zm_levels.size != expected:
            raise ValueError(
                f"zm grid file {zm_grid_path} has {zm_levels.size} levels; "
                f"expected nzmax={expected}."
            )
        momentum_heights = np.tile(zm_levels[None, :], (ngrdcol, 1))
    else:
        raise ValueError(f"Unsupported grid_type: {grid_type}")

    gr = py_setup_grid(
        ngrdcol=ngrdcol,
        deltaz=deltaz,
        zm_init=zm_init,
        zm_top=zm_top,
        l_ascending_grid=True,
        grid_type=grid_type,
        momentum_heights=momentum_heights,
        thermodynamic_heights=thermodynamic_heights,
    )

    # Use grid dimensions from Python grid construction.
    nzm = gr.nzm
    nzt = gr.nzt

    print(f"nzm = {nzm} -- nzt = {nzt}")

    # ── 5. Interpolate sounding onto grid ───────────────────────────────
    # Use first column's zt for interpolation
    zt_1d = gr.zt[0, :]  # (nzt,)
    use_cubic_ic = bool(cfg.get('l_modify_ic_with_cubic_int', False))
    snd_interp = interpolate_sounding(snd, zt_1d, use_cubic=use_cubic_ic)

    # Build 2D arrays (ngrdcol, nzt) by broadcasting 1D profile
    thlm = np.tile(snd_interp['theta'], (ngrdcol, 1))
    rtm = np.tile(snd_interp['rt'], (ngrdcol, 1))
    um = np.tile(snd_interp['u'], (ngrdcol, 1))
    vm = np.tile(snd_interp['v'], (ngrdcol, 1))
    ug = np.tile(snd_interp['ug'], (ngrdcol, 1))
    vg = np.tile(snd_interp['vg'], (ngrdcol, 1))
    wm_zt = np.tile(snd_interp['w'], (ngrdcol, 1))
    p_in_Pa = np.zeros((ngrdcol, nzt))  # will be computed

    # ── 6. Initialize pressure / thermodynamic variables ────────────────
    p_sfc = np.full(ngrdcol, cfg['p_sfc_nl'])
    T0 = cfg['t0']
    fcor_nl = cfg['fcor_nl']
    lat_vals = cfg['lat_vals']

    fcor = np.full(ngrdcol, fcor_nl)
    fcor_y = np.full(ngrdcol, 2.0 * omega_planet * math.cos(lat_vals * radians_per_deg))

    # Compute initial thvm (approximation: thvm = thlm * (1 + ep1 * rv))
    # where rv = rtm / (1 + rtm)
    thvm = thlm * (1.0 + ep1 * (rtm / (1.0 + rtm)))

    # Hydrostatic pressure
    result = hydrostatic(thvm=thvm, p_sfc=p_sfc, gr=gr)
    p_in_Pa, p_in_Pa_zm, exner, exner_zm, rho, rho_zm = result
    p_in_Pa = np.array(p_in_Pa, copy=True)
    p_in_Pa_zm = np.array(p_in_Pa_zm, copy=True)
    exner = np.array(exner, copy=True)
    exner_zm = np.array(exner_zm, copy=True)
    rho = np.array(rho, copy=True)
    rho_zm = np.array(rho_zm, copy=True)

    # Convert temperature type
    if temperature_type in ('thm[K]', 'T[K]'):
        # theta sounding — need to compute rcm and convert to thlm
        thm = thlm.copy()
        # rcm = max(rtm - rsat(p, T), 0)
        T_in_K = thm * exner
        rsat = np.array(sat_mixrat_liq(
            p_in_Pa=p_in_Pa,
            T_in_K=T_in_K,
            saturation_formula=saturation_formula,
        ), copy=True)
        rcm = np.maximum(rtm - rsat, 0.0)
        thlm = thm - Lv / (Cp * exner) * rcm
    elif temperature_type == 'thlm[K]':
        # Already liquid potential temperature
        rcm = np.array(rcm_sat_adj(thlm, rtm, p_in_Pa, exner, saturation_formula), copy=True)
        thm = thlm + Lv / (Cp * exner) * rcm
    else:
        raise ValueError(f"Unknown temperature_type: {temperature_type}")

    # Recompute thvm and hydrostatic with corrected thlm
    # NOTE: Fortran passes thm (not thlm) as the 5th arg to calculate_thvm
    thvm = np.array(calculate_thvm(
        nzt=nzt, ngrdcol=ngrdcol, thlm=thlm, rtm=rtm, rcm=rcm, exner=exner,
        thv_ds_zt=thm * (1.0 + ep2 * (rtm - rcm))**kappa,
    ), copy=True)
    result = hydrostatic(thvm=thvm, p_sfc=p_sfc, gr=gr)
    p_in_Pa, p_in_Pa_zm, exner, exner_zm, rho, rho_zm = result
    p_in_Pa = np.array(p_in_Pa, copy=True)
    p_in_Pa_zm = np.array(p_in_Pa_zm, copy=True)
    exner = np.array(exner, copy=True)
    exner_zm = np.array(exner_zm, copy=True)
    rho = np.array(rho, copy=True)
    rho_zm = np.array(rho_zm, copy=True)

    # Compute dry static density (anelastic base state)
    # NOTE: thm was already computed before the 2nd hydrostatic call
    # (from sounding for thm[K] case, or from thlm+Lv/(Cp*exner)*rcm for
    # thlm[K] case). Do NOT recompute it here with the updated exner —
    # the Fortran uses the original thm throughout.
    rv = rtm - rcm  # water vapor mixing ratio
    p_dry = p_in_Pa / (1.0 + ep2 * rv)
    exner_dry = (p_dry / p0)**kappa
    th_dry = thm * (1.0 + ep2 * rv)**kappa
    rho_dry = p_dry / (Rd * th_dry * exner_dry)

    rho_ds_zt = rho_dry.copy()
    thv_ds_zt = th_dry.copy()
    invrs_rho_ds_zt = 1.0 / rho_ds_zt

    # Momentum level versions via zt2zm interpolation
    rv_zm = np.array(zt2zm(nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, gr=gr, azt=rv), copy=True)
    rv_zm = np.maximum(rv_zm, 0.0)
    thm_zm = np.array(zt2zm(nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, gr=gr, azt=thm), copy=True)

    # rtm_sfc: linearly interpolate sounding rt to the zm surface level,
    # matching Fortran read_sounding which interpolates to gr%zm(1).
    zm_sfc = gr.zm[0, 0]  # surface momentum level (z=0 typically)
    z_snd = snd['z']
    rt_snd = snd['rt']
    valid_rt = rt_snd > -998.0
    if np.sum(valid_rt) >= 2 and zm_sfc >= z_snd[valid_rt][0]:
        rtm_sfc = float(np.interp(zm_sfc, z_snd[valid_rt], rt_snd[valid_rt]))
    else:
        rtm_sfc = float(rtm[0, 0])  # fallback: use lowest zt level
    pd_sfc = p_sfc / (1.0 + ep2 * rtm_sfc)

    p_dry_zm = p_in_Pa_zm / (1.0 + ep2 * rv_zm)
    p_dry_zm[:, 0] = pd_sfc
    exner_dry_zm = (p_dry_zm / p0)**kappa
    th_dry_zm = thm_zm * (1.0 + ep2 * rv_zm)**kappa
    rho_dry_zm = p_dry_zm / (Rd * th_dry_zm * exner_dry_zm)

    rho_ds_zm = rho_dry_zm.copy()
    thv_ds_zm = th_dry_zm.copy()
    invrs_rho_ds_zm = 1.0 / rho_ds_zm

    # ── 7. Subsidence / vertical wind ───────────────────────────────────
    if subs_type == 'omega[Pa/s]':
        wm_zt = -wm_zt / (grav * rho)
        wm_zt[:, -1] = 0.0

    wm_zm = np.array(zt2zm(nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, gr=gr, azt=wm_zt), copy=True)
    wm_zm[:, 0] = 0.0
    wm_zm[:, -1] = 0.0

    # ── 8. Initialize PDF and tunable parameters ────────────────────────
    clubb_params_all = init_clubb_params(total_param_sets, filename=namelist_path)
    clubb_params = clubb_params_all[:ngrdcol]
    pdf_params = init_pdf_params_py(nzt, ngrdcol)
    pdf_params_zm = init_pdf_params_py(nzm, ngrdcol)   # NB: Fortran uses nzm for pdf_params_zm
    pdf_implicit_coefs_terms = init_pdf_implicit_coefs_terms_api(nzt, ngrdcol, sclr_dim)

    # Scalar indices (mirror initialize_clubb defaults/namelist overrides).
    iisclr_rt = int(cfg.get('iisclr_rt', -1))
    iisclr_thl = int(cfg.get('iisclr_thl', -1))
    iisclr_co2 = int(cfg.get('iisclr_co2', -1))
    iiedsclr_rt = int(cfg.get('iiedsclr_rt', -1))
    iiedsclr_thl = int(cfg.get('iiedsclr_thl', -1))
    iiedsclr_co2 = int(cfg.get('iiedsclr_co2', -1))
    sclr_idx = SclrIdx(
        iisclr_rt=iisclr_rt,
        iisclr_thl=iisclr_thl,
        iisclr_CO2=iisclr_co2,
        iiedsclr_rt=iiedsclr_rt,
        iiedsclr_thl=iiedsclr_thl,
        iiedsclr_CO2=iiedsclr_co2,
    )
    nu_vert_res_dep, lmin, mixt_frac_max_mag = calc_derived_params(
        gr, ngrdcol, grid_type, deltaz,               # In
        clubb_params, flags.l_prescribed_avg_deltaz,  # In
    )
    if float(lmin) < 1.0:
        raise ValueError('lmin is < 1.0')

    err_info = check_clubb_settings(
        ngrdcol=ngrdcol,
        params=clubb_params,
        config_flags=flags,
        err_info=err_info,
        l_implemented=False,
        l_input_fields=False,
    )
    # Reset error codes set by check_clubb_settings warnings
    # (they are non-fatal but would cause advance_clubb_core to bail out)
    err_info = err_info.reset_code()

    # ── 9. Initialize TKE / variances (case-specific, Fortran-like) ────
    em, wp2, up2, vp2, upwp, um = _initialize_turbulence_state(
        runtype=runtype,
        gr=gr,
        dt_main=dt_main,
        fcor_y=fcor_y,
        um=um,
    )

    # ── 10. Initialize remaining prognostic arrays to zero ──────────────

    # Reference profiles (matches initialize_clubb logic in Fortran driver)
    uv_sponge_enabled = bool(cfg.get('uv_sponge_damp_settings%l_sponge_damping', False))
    if flags.l_uv_nudge or uv_sponge_enabled:
        um_ref = um.copy()
        vm_ref = vm.copy()
    else:
        um_ref = np.zeros((ngrdcol, nzt))
        vm_ref = np.zeros((ngrdcol, nzt))
    thlm_ref = np.zeros((ngrdcol, nzt))
    rtm_ref = np.zeros((ngrdcol, nzt))

    # Cloud properties
    nc0_in_cloud = float(cfg.get('nc0_in_cloud', Nc0_in_cloud))
    Nc_in_cloud = nc0_in_cloud / rho
    cloud_frac = np.zeros((ngrdcol, nzt))
    Ncm = np.where(rcm > 0, Nc_in_cloud, Nc_in_cloud * cloud_frac_min)

    # Resolve standalone output names once, before microphysics initialization
    # appends its source configuration report to the case-info file.
    repo_root = _repo_root()
    stats_prefix = str(cfg.get('fname_prefix', '')).strip() or runtype
    output_dir_raw = str(cfg.get("output_dir", "")).strip().strip("'\"")
    if output_dir_raw:
        output_dir_path = Path(output_dir_raw)
        if not output_dir_path.is_absolute():
            output_dir_path = (namelist_dir / output_dir_path).resolve()
    else:
        output_dir_path = repo_root / "output"
    case_info_file = output_dir_path / f"{stats_prefix}_setup.txt"
    if cfg['debug_level'] >= 1:
        case_info_file.parent.mkdir(parents=True, exist_ok=True)
        case_info_file.write_text("")

    from clubb_jax.src.Microphys.microphys_init_cleanup import init_microphys
    (
        hydromet_dim, pdf_dim, hm_metadata, silhs_config_flags, vert_decorr_coef,
        corr_array_n_cloud, corr_array_n_below,
    ) = init_microphys(
        0, runtype, cfg, case_info_file,      # In
        1.0e6, 1.0e6,                        # In
        clubb_params,                        # In
        flags.l_diagnose_correlations,       # In
        flags.l_const_Nc_in_cloud,           # InOut
        flags.l_fix_w_chi_eta_correlations,  # InOut
    )
    from clubb_jax.src.Microphys import parameters_microphys
    from clubb_jax.src.SILHS.latin_hypercube_arrays import LatinHypercubeArrays
    num_samples = (
        parameters_microphys.lh_num_samples
        if parameters_microphys.lh_microphys_type != parameters_microphys.lh_microphys_disabled
        else 0
    )
    # JAX adaptation of threadprivate SILHS storage: fixed-shape case-owned
    # permutation arrays are carried through compiled sampling calls.
    sampling_state = LatinHypercubeArrays(
        jnp.zeros(
            (
                parameters_microphys.lh_num_samples * parameters_microphys.lh_sequence_length,
                pdf_dim + 2,
            ),
            dtype=jnp.int32,
        ),
        jnp.array(0, dtype=jnp.int32),
    )
    # Keep a padded trailing extent for zero-species arrays because JAX kernels
    # cannot index a physically empty axis. The logical *_dim values remain authoritative.
    hm_dim_transport = max(hydromet_dim, 1)
    l_mix_rat_hm = hm_metadata.l_mix_rat_hm if hydromet_dim else np.zeros((hm_dim_transport,), dtype=bool)
    wphydrometp = np.zeros((ngrdcol, nzm, hm_dim_transport))
    wp2hmp = np.zeros((ngrdcol, nzt, hm_dim_transport))
    rtphmp_zt = np.zeros((ngrdcol, nzt, hm_dim_transport))
    thlphmp_zt = np.zeros((ngrdcol, nzt, hm_dim_transport))

    sc_dim_transport = max(sclr_dim, 1)
    edsc_dim_transport = max(edsclr_dim, 1)
    sclr_tol = np.array(cfg.get('sclr_tol_nl', [])[:sclr_dim], dtype=np.float64)
    if len(sclr_tol) < sclr_dim:
        sclr_tol = np.pad(sclr_tol, (0, sclr_dim - len(sclr_tol)), constant_values=1e-8)
    sclrm = np.zeros((ngrdcol, nzt, sc_dim_transport))
    sclrp2 = np.zeros((ngrdcol, nzm, sc_dim_transport))
    if sclr_dim > 0:
        sclrp2[:, :, :sclr_dim] = sclr_tol[:sclr_dim].reshape(1, 1, sclr_dim) ** 2
    sclrp3 = np.zeros((ngrdcol, nzt, sc_dim_transport))
    sclrprtp = np.zeros((ngrdcol, nzm, sc_dim_transport))
    sclrpthlp = np.zeros((ngrdcol, nzm, sc_dim_transport))
    sclrpthvp = np.zeros((ngrdcol, nzm, sc_dim_transport))
    wpsclrp = np.zeros((ngrdcol, nzm, sc_dim_transport))
    sclrm_forcing = jnp.zeros((ngrdcol, nzt, sc_dim_transport))
    wpsclrp_sfc = jnp.zeros((ngrdcol, sc_dim_transport))

    edsclrm = np.zeros((ngrdcol, nzt, edsc_dim_transport))
    edsclrm_forcing = jnp.zeros((ngrdcol, nzt, edsc_dim_transport))
    wpedsclrp_sfc = jnp.zeros((ngrdcol, edsc_dim_transport))

    # Initialize scalar means from dedicated scalar sounding files.
    if sclr_dim > 0:
        sclr_path = _resolve_case_input_path(namelist_dir, runtype, "_sclr_sounding.in")
        sclr_raw = read_scalar_sounding(str(sclr_path), sclr_dim)
        _validate_scalar_column_names(
            sclr_raw['names'], iisclr_rt, iisclr_thl, iisclr_co2, label='sclr_sounding'
        )
        sclr_zt = interpolate_scalar_sounding(
            snd['z'], sclr_raw['data'], zt_1d, use_cubic=use_cubic_ic
        )
        sclrm[:, :, :sclr_dim] = np.tile(sclr_zt[None, :, :], (ngrdcol, 1, 1))

    if edsclr_dim > 0:
        edsclr_path = _resolve_case_input_path(namelist_dir, runtype, "_edsclr_sounding.in")
        edsclr_raw = read_scalar_sounding(str(edsclr_path), edsclr_dim)
        _validate_scalar_column_names(
            edsclr_raw['names'], iiedsclr_rt, iiedsclr_thl, iiedsclr_co2, label='edsclr_sounding'
        )
        edsclr_zt = interpolate_scalar_sounding(
            snd['z'], edsclr_raw['data'], zt_1d, use_cubic=use_cubic_ic
        )
        edsclrm[:, :, :edsclr_dim] = np.tile(edsclr_zt[None, :, :], (ngrdcol, 1, 1))

    # ── 11. Read surface file ───────────────────────────────────────────
    sfc_path = None
    try:
        sfc_path = _resolve_case_input_path(namelist_dir, runtype, "_sfc.in")
    except FileNotFoundError:
        sfc_path = None
    sfc_data = None
    if sfc_path is not None and sfc_path.exists():
        sfc_data = read_surface(str(sfc_path))

    # ── 12. Time controls ───────────────────────────────────────────────
    time_initial = cfg['time_initial']
    time_final = cfg['time_final']
    ifinal = int(math.floor((time_final - time_initial) / dt_main))
    iinit = 1
    if bool(cfg.get('l_restart', False)):
        time_restart = float(cfg.get('time_restart', 0.0))
        if abs(math.fmod(time_restart - time_initial, dt_main)) > eps:
            raise ValueError("(time_restart-time_initial) is not a multiple of dt_main")
        if not time_initial <= time_restart < time_final:
            raise ValueError("time_restart must lie in [time_initial, time_final)")
        # The value is increased by 1 to synchronize with restart data.
        iinit = math.floor((time_restart - time_initial) / dt_main) + 1
        restart_path_case = Path(_clean_namelist_path(cfg['restart_path_case']))
        if not restart_path_case.is_absolute():
            # Native restart paths are relative to the repository root, while
            # generated output_dir paths are relative to the aggregate namelist.
            restart_path_case = _repo_root() / restart_path_case
    stats_nsamp = int(round(cfg['stats_tsamp'] / dt_main))
    stats_nout = int(round(cfg['stats_tout'] / dt_main))

    # ── 13. Initialize stats ────────────────────────────────────────────
    l_stats = bool(cfg['l_stats'])
    stats_registry_path = _resolve_stats_registry_path(namelist_path, cfg)
    stats_output_path = output_dir_path / f"{stats_prefix}_stats.nc"
    # An empty source stats_output_filename keeps statistics in memory.
    if 'stats_output_filename' in cfg:
        stats_filename = str(cfg['stats_output_filename']).strip()
        stats_output_path = output_dir_path / stats_filename if stats_filename else None
    if l_stats and stats_output_path is not None and bool(cfg.get('l_restart', False)):
        from clubb_jax.src.Input_fields import input_fields

        # Opening a writer would truncate a restart reference in the same path.
        for restart_stats_path in input_fields.set_filenames(restart_path_case):
            if (
                stats_output_path.resolve() == restart_stats_path.resolve()
                or (
                    stats_output_path.exists() and restart_stats_path.exists()
                    and stats_output_path.samefile(restart_stats_path)
                )
            ):
                raise ValueError("Restart reference and output statistics must use different paths")

    stats_writer = None
    if l_stats:
        if not stats_registry_path.exists():
            raise FileNotFoundError(f"Stats registry file not found: {stats_registry_path}")
        if stats_output_path is not None:
            stats_output_path.parent.mkdir(parents=True, exist_ok=True)
        stats_writer = StatsWriter(
            registry_path=str(stats_registry_path),
            output_path=str(stats_output_path) if stats_output_path is not None else '',
            nzt=nzt,
            nzm=nzm,
            ngrdcol=ngrdcol,
            zt=np.asarray(gr.zt[0, :]),
            zm=np.asarray(gr.zm[0, :]),
            stats_tsamp=float(cfg['stats_tsamp']),
            stats_tout=float(cfg['stats_tout']),
            dt_main=float(dt_main),
            day=int(cfg['day']),
            month=int(cfg['month']),
            year=int(cfg['year']),
            time_initial=float(time_initial),
            stats_tstart=float(cfg.get('stats_tstart', time_initial)),
            stats_tend=float(cfg.get('stats_tend', time_final)),
            ncol_total=total_param_sets,
            clubb_params_vals=np.asarray(clubb_params_all),
            param_names=get_param_names(),
            sclr_dim=sclr_dim,
            edsclr_dim=edsclr_dim,
            hydromet_list=hm_metadata.hydromet_list,
        )
        if not stats_writer.enabled:
            raise RuntimeError("stats_init completed but stats are not enabled")

    if (
        total_param_sets > batch_size and stats_writer is not None
        and stats_output_path is not None
        and parameters_microphys.lh_microphys_type != parameters_microphys.lh_microphys_disabled
    ):
        # Source batch-mode NetCDF output does not support SILHS sample output.
        stats_writer.finalize()
        raise ValueError('Batch-mode stats NetCDF output does not yet support SILHS sample output')

    if parameters_microphys.lh_microphys_type != parameters_microphys.lh_microphys_disabled:
        from clubb_jax.src.SILHS.silhs_api_module import latin_hypercube_2D_output_api

        # Setup 2D output of all subcolumns (if enabled). The source API also
        # defines SILHS coordinates when both sample-output flags are false.
        stats_writer, err_info = latin_hypercube_2D_output_api(
            nzt, gr.zt[0, :], pdf_dim, num_samples, hm_metadata,  # In
            stats_writer, err_info,                             # InOut
        )
        if err_info.is_fatal():
            raise RuntimeError("Fatal error calling latin_hypercube_2D_output_api in init_clubb_case")

    # ── 14. Zero PDF params ─────────────────────────────────────────────
    pdf_params = init_pdf_params_py(nzt, ngrdcol)
    pdf_params_zm = init_pdf_params_py(nzm, ngrdcol)

    # ── 15. Clear any accumulated error codes from init ──────────────
    err_info = err_info.reset_code()

    # The Fortran driver initializes these fields for every case before the
    # timestep loop; prescribe_forcings therefore receives veg_T_in_K even when
    # the selected case does not use the interactive soil/vegetation scheme.
    deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K = initialize_soil_veg(ngrdcol)
    radiation_parameters = initialize_radiation_parameters(cfg, rad_scheme, namelist_dir)

    # ── Build state dict ────────────────────────────────────────────────
    state = dict(
        # Config
        cfg=cfg, flags=flags, gr=gr, namelist_dir=str(namelist_dir),
        runtype=runtype, ngrdcol=ngrdcol, nzt=nzt, nzm=nzm,
        dt_main=dt_main, dt_rad=dt_rad,
        time_initial=time_initial, time_final=time_final,
        iinit=iinit, ifinal=ifinal, l_stats=l_stats, stats_writer=stats_writer,
        stats_nsamp=stats_nsamp, stats_nout=stats_nout,
        stats_registry_path=str(stats_registry_path),
        stats_output_path=str(stats_output_path) if stats_output_path is not None else '',
        total_param_sets=total_param_sets, clubb_params_all=clubb_params_all,
        saturation_formula=saturation_formula,
        sfctype=int(cfg['sfctype']),
        microphys_scheme=microphys_scheme,
        l_cloud_sed=l_cloud_sed,
        sigma_g=sigma_g,
        nc0_in_cloud=nc0_in_cloud,
        rad_scheme=rad_scheme,
        radiation_parameters=radiation_parameters,
        day=int(cfg['day']), month=int(cfg['month']), year=int(cfg['year']),
        lat_vals=float(cfg['lat_vals']), lon_vals=float(cfg['lon_vals']),
        l_calc_thlp2_rad=l_calc_thlp2_rad,
        hydromet_dim=hydromet_dim, sclr_dim=sclr_dim, edsclr_dim=edsclr_dim,
        T0=T0, lmin=lmin, mixt_frac_max_mag=mixt_frac_max_mag,
        ts_nudge=cfg['ts_nudge'],
        rtm_min=cfg['rtm_min'],
        rtm_nudge_max_altitude=cfg['rtm_nudge_max_altitude'],
        l_t_dependent=bool(cfg.get('l_t_dependent', False)),
        l_ignore_forcings=bool(cfg.get('l_ignore_forcings', False)),
        l_input_xpwp_sfc=bool(cfg.get('l_input_xpwp_sfc', False)),
        iisclr_rt=iisclr_rt, iisclr_thl=iisclr_thl, iisclr_co2=iisclr_co2,
        iiedsclr_rt=iiedsclr_rt, iiedsclr_thl=iiedsclr_thl, iiedsclr_co2=iiedsclr_co2,
        sclr_idx=sclr_idx,
        nu_vert_res_dep=nu_vert_res_dep,
        pdf_params=pdf_params,
        pdf_params_zm=pdf_params_zm,
        pdf_implicit_coefs_terms=pdf_implicit_coefs_terms,
        err_info=err_info,
        l_modify_bc_for_cnvg_test=bool(cfg.get('l_modify_bc_for_cnvg_test', False)),
        sfc_data=sfc_data,
        # 1D arrays
        fcor=fcor, fcor_y=fcor_y, sfc_elevation=sfc_elevation,
        p_sfc=p_sfc,
        deep_soil_T_in_K=deep_soil_T_in_K,
        sfc_soil_T_in_K=sfc_soil_T_in_K,
        veg_T_in_K=veg_T_in_K,
        host_dx=np.full(ngrdcol, 1.0e6),
        host_dy=np.full(ngrdcol, 1.0e6),
        upwp_sfc_pert=np.zeros((ngrdcol,)),
        vpwp_sfc_pert=np.zeros((ngrdcol,)),
        sclr_tol=sclr_tol,
        l_mix_rat_hm=l_mix_rat_hm,
        clubb_params=clubb_params,
        # Prognostic zt (ngrdcol, nzt)
        um=um, vm=vm, thlm=thlm, rtm=rtm,
        up3=np.zeros((ngrdcol, nzt)), vp3=np.zeros((ngrdcol, nzt)),
        rtp3=np.zeros((ngrdcol, nzt)), thlp3=np.zeros((ngrdcol, nzt)),
        wp3=np.zeros((ngrdcol, nzt)),
        p_in_Pa=p_in_Pa, exner=exner, rcm=rcm,
        cloud_frac=cloud_frac,
        wp2thvp=np.zeros((ngrdcol, nzt)), wp2up=np.zeros((ngrdcol, nzt)),
        wp2rtp=np.zeros((ngrdcol, nzt)), wp2thlp=np.zeros((ngrdcol, nzt)),
        wpup2=np.zeros((ngrdcol, nzt)), wpvp2=np.zeros((ngrdcol, nzt)),
        ice_supersat_frac=np.zeros((ngrdcol, nzt)),
        um_pert=np.zeros((ngrdcol, nzt)), vm_pert=np.zeros((ngrdcol, nzt)),
        # Prognostic zm (ngrdcol, nzm)
        upwp=upwp, vpwp=np.zeros((ngrdcol, nzm)),
        up2=up2, vp2=vp2,
        wprtp=np.zeros((ngrdcol, nzm)), wpthlp=np.zeros((ngrdcol, nzm)),
        rtp2=np.full((ngrdcol, nzm), rt_tol**2),
        thlp2=np.full((ngrdcol, nzm), thl_tol**2),
        rtpthlp=np.zeros((ngrdcol, nzm)),
        wp2=wp2,
        wpthvp=np.zeros((ngrdcol, nzm)), rtpthvp=np.zeros((ngrdcol, nzm)),
        thlpthvp=np.zeros((ngrdcol, nzm)),
        uprcp=np.zeros((ngrdcol, nzm)), vprcp=np.zeros((ngrdcol, nzm)),
        rc_coef_zm=np.zeros((ngrdcol, nzm)),
        wp4=np.zeros((ngrdcol, nzm)),
        wp2up2=np.zeros((ngrdcol, nzm)), wp2vp2=np.zeros((ngrdcol, nzm)),
        upwp_pert=np.zeros((ngrdcol, nzm)), vpwp_pert=np.zeros((ngrdcol, nzm)),
        # Forcing arrays
        thlm_forcing=np.zeros((ngrdcol, nzt)), rtm_forcing=np.zeros((ngrdcol, nzt)),
        um_forcing=np.zeros((ngrdcol, nzt)), vm_forcing=np.zeros((ngrdcol, nzt)),
        wprtp_forcing=np.zeros((ngrdcol, nzm)), wpthlp_forcing=np.zeros((ngrdcol, nzm)),
        rtp2_forcing=np.zeros((ngrdcol, nzm)), thlp2_forcing=np.zeros((ngrdcol, nzm)),
        rtpthlp_forcing=np.zeros((ngrdcol, nzm)),
        # Meteorological profiles
        wm_zt=wm_zt, wm_zm=wm_zm,
        rho=rho, rho_zm=rho_zm,
        Ncm=Ncm, Nc_in_cloud=Nc_in_cloud,
        rho_ds_zt=rho_ds_zt, rho_ds_zm=rho_ds_zm,
        invrs_rho_ds_zt=invrs_rho_ds_zt, invrs_rho_ds_zm=invrs_rho_ds_zm,
        thv_ds_zt=thv_ds_zt, thv_ds_zm=thv_ds_zm,
        thvm=thvm,
        radht=np.zeros((ngrdcol, nzt)),
        # The source restart reader restores radht only. Fluxes and separate
        # heating caches retain these zeros until the next radiation update.
        radht_SW=np.zeros((ngrdcol, nzt)),
        radht_LW=np.zeros((ngrdcol, nzt)),
        Frad=np.zeros((ngrdcol, nzm)),
        Frad_SW=np.zeros((ngrdcol, nzm)),
        Frad_LW=np.zeros((ngrdcol, nzm)),
        Frad_SW_up=np.zeros((ngrdcol, nzm)),
        Frad_LW_up=np.zeros((ngrdcol, nzm)),
        Frad_SW_down=np.zeros((ngrdcol, nzm)),
        Frad_LW_down=np.zeros((ngrdcol, nzm)),
        rcm_mc=np.zeros((ngrdcol, nzt)),
        thlm_mc=np.zeros((ngrdcol, nzt)),
        rfrzm=np.zeros((ngrdcol, nzt)),
        # Reference profiles
        um_ref=um_ref, vm_ref=vm_ref,
        thlm_ref=thlm_ref, rtm_ref=rtm_ref,
        ug=ug, vg=vg,
        # Hydromet
        wphydrometp=wphydrometp,
        wp2hmp=wp2hmp, rtphmp_zt=rtphmp_zt, thlphmp_zt=thlphmp_zt,
        hydromet=np.zeros((ngrdcol, nzt, hydromet_dim)),
        hm_metadata=hm_metadata, pdf_dim=pdf_dim,
        silhs_config_flags=silhs_config_flags, vert_decorr_coef=vert_decorr_coef,
        sampling_state=sampling_state,
        corr_array_n_cloud=corr_array_n_cloud, corr_array_n_below=corr_array_n_below,
        hydrometp2=jnp.zeros((ngrdcol, nzm, hydromet_dim)),
        K_hm=jnp.zeros((ngrdcol, nzm, hydromet_dim)),
        hydromet_vel_zt=jnp.zeros((ngrdcol, nzt, hydromet_dim)),
        Nccnm=jnp.zeros((ngrdcol, nzt)),
        rvm_mc=jnp.zeros((ngrdcol, nzt)),
        wprtp_mc=jnp.zeros((ngrdcol, nzm)), wpthlp_mc=jnp.zeros((ngrdcol, nzm)),
        rtp2_mc=jnp.zeros((ngrdcol, nzm)), thlp2_mc=jnp.zeros((ngrdcol, nzm)),
        rtpthlp_mc=jnp.zeros((ngrdcol, nzm)),
        X_nl_all_levs=np.zeros((ngrdcol, num_samples, nzt, pdf_dim)),
        X_mixt_comp_all_levs=np.zeros((ngrdcol, num_samples, nzt), dtype=np.int32),
        lh_rt_clipped=np.zeros((ngrdcol, num_samples, nzt)),
        lh_thl_clipped=np.zeros((ngrdcol, num_samples, nzt)),
        lh_rc_clipped=np.zeros((ngrdcol, num_samples, nzt)),
        lh_rv_clipped=np.zeros((ngrdcol, num_samples, nzt)),
        lh_Nc_clipped=np.zeros((ngrdcol, num_samples, nzt)),
        lh_sample_point_weights=np.zeros((ngrdcol, num_samples, nzt)),
        # Scalars
        sclrm=sclrm, sclrp2=sclrp2, sclrp3=sclrp3,
        sclrprtp=sclrprtp, sclrpthlp=sclrpthlp, sclrpthvp=sclrpthvp,
        wpsclrp=wpsclrp,
        sclrm_forcing=sclrm_forcing,
        wpsclrp_sfc=wpsclrp_sfc,
        edsclrm=edsclrm, edsclrm_forcing=edsclrm_forcing,
        wpedsclrp_sfc=wpedsclrp_sfc,
        # Surface fluxes (will be set by forcings)
        wpthlp_sfc=np.zeros((ngrdcol,)),
        wprtp_sfc=np.zeros((ngrdcol,)),
        upwp_sfc=np.zeros((ngrdcol,)),
        vpwp_sfc=np.zeros((ngrdcol,)),
        T_sfc=np.full(ngrdcol, float(cfg.get('t_sfc_nl', 288.0))),
        sens_ht=float(cfg.get('sens_ht', 0.0)),
        latent_ht=float(cfg.get('latent_ht', 0.0)),
        # Output / diagnostic
        thlprcp=np.zeros((ngrdcol, nzm)),
    )

    from clubb_jax.src.CLUBB_core.hydromet_pdf_parameter_module import init_precip_fracs
    state['precip_fracs'] = init_precip_fracs(nzt, ngrdcol)

    # Initialize Time Dependent Input
    time_dependent_input.l_t_dependent = state['l_t_dependent']
    time_dependent_input.l_ignore_forcings = state['l_ignore_forcings']
    time_dependent_input.l_input_xpwp_sfc = state['l_input_xpwp_sfc']
    if state['l_t_dependent']:
        time_dependent_input.initialize_t_dependent_input(
            0,
            runtype,
            nzt,
            np.asarray(gr.zt[0, :]),
            np.asarray(p_in_Pa)[0, :],
            int(cfg.get('grid_adapt_in_time_method', 0)),
        )

    if runtype == 'mpace_a':
        mpace_a_init(0, _repo_root() / 'input' / 'case_setups' / 'mpace_a_forcings')

    if bool(cfg.get('l_restart', False)):
        # initialize_clubb includes reference-profile setup as well as sounding
        # reads. Execute it first; the restart overwrites the initial sounding.
        # The source reader omits scalar state and some frozen hydrometeors;
        # see RESTART_PORT_NOTES.md for exact-continuation limits.
        from clubb_jax.src.Input_fields import input_fields

        input_fields.clubb_day = state['day']
        input_fields.clubb_month = state['month']
        input_fields.clubb_year = state['year']
        input_fields.l_soil_veg = radiation_parameters.l_soil_veg
        for name in ('em', 'tau_zm', 'Kh_zm', 'sigma_sqd_w'):
            state[name] = jnp.zeros((ngrdcol, nzm))
        for name in ('tau_zt', 'Kh_zt', 'sigma_sqd_w_zt'):
            state[name] = jnp.zeros((ngrdcol, nzt))

        try:
            (
                state["um"], state["upwp"], state["vm"], state["vpwp"], state["up2"], state["vp2"],
                state["rtm"],
                state["wprtp"], state["thlm"], state["wpthlp"], state["rtp2"], state["rtp3"],
                state["thlp2"], state["thlp3"], state["rtpthlp"], state["wp2"], state["wp3"],
                state["p_in_Pa"], state["exner"], state["rcm"], state["cloud_frac"],
                state["wpthvp"], state["wp2thvp"], state["wp2up"], state["rtpthvp"],
                state["thlpthvp"],
                state["wp2rtp"], state["wp2thlp"], state["uprcp"], state["vprcp"],
                state["rc_coef_zm"], state["wp4"], state["wpup2"], state["wpvp2"], state["wp2up2"],
                state["wp2vp2"], state["ice_supersat_frac"],
                state["wm_zt"], state["rho"], state["rho_zm"], state["rho_ds_zm"],
                state["rho_ds_zt"], state["thv_ds_zm"], state["thv_ds_zt"],
                state["thlm_forcing"], state["rtm_forcing"], state["wprtp_forcing"],
                state["wpthlp_forcing"], state["rtp2_forcing"],
                state["thlp2_forcing"], state["rtpthlp_forcing"],
                state["hydromet"], state["hydrometp2"], state["wphydrometp"],
                state["Ncm"], state["Nccnm"], state["thvm"], state["em"], state["tau_zm"],
                state["tau_zt"],
                state["Kh_zt"], state["Kh_zm"], state["ug"], state["vg"],
                state["thlprcp"],
                state["sigma_sqd_w"], state["sigma_sqd_w_zt"], state["radht"],
                state["deep_soil_T_in_K"], state["sfc_soil_T_in_K"], state["veg_T_in_K"],
                state["pdf_params"], state["pdf_params_zm"],
                state["rcm_mc"], state["rvm_mc"], state["thlm_mc"],
                state["wprtp_mc"], state["wpthlp_mc"], state["rtp2_mc"],
                state["thlp2_mc"], state["rtpthlp_mc"],
                state["wpthlp_sfc"], state["wprtp_sfc"], state["upwp_sfc"], state["vpwp_sfc"],
            ) = restart_clubb(
                gr, hydromet_dim, hm_metadata,                                             # In
                restart_path_case, time_restart,                                           # In
                state["um"], state["upwp"], state["vm"], state["vpwp"], state["up2"], state["vp2"],       # InOut
                state["rtm"],                                                                             # InOut
                state["wprtp"], state["thlm"], state["wpthlp"], state["rtp2"], state["rtp3"],             # InOut
                state["thlp2"], state["thlp3"], state["rtpthlp"], state["wp2"], state["wp3"],             # InOut
                state["p_in_Pa"], state["exner"], state["rcm"], state["cloud_frac"],                      # InOut
                state["wpthvp"], state["wp2thvp"], state["wp2up"], state["rtpthvp"],                      # InOut
                state["thlpthvp"],                                                                        # InOut
                state["wp2rtp"], state["wp2thlp"], state["uprcp"], state["vprcp"],                        # InOut
                state["rc_coef_zm"], state["wp4"], state["wpup2"], state["wpvp2"], state["wp2up2"],       # InOut
                state["wp2vp2"], state["ice_supersat_frac"],                                              # InOut
                state["wm_zt"], state["rho"], state["rho_zm"], state["rho_ds_zm"],                        # InOut
                state["rho_ds_zt"], state["thv_ds_zm"], state["thv_ds_zt"],                               # InOut
                state["thlm_forcing"], state["rtm_forcing"], state["wprtp_forcing"],                      # InOut
                state["wpthlp_forcing"], state["rtp2_forcing"],                                           # InOut
                state["thlp2_forcing"], state["rtpthlp_forcing"],                                         # InOut
                state["hydromet"], state["hydrometp2"], state["wphydrometp"],                             # InOut
                state["Ncm"], state["Nccnm"], state["thvm"], state["em"], state["tau_zm"],                # InOut
                state["tau_zt"],                                                                          # InOut
                state["Kh_zt"], state["Kh_zm"], state["ug"], state["vg"],                                 # InOut
                state["thlprcp"],                                                                         # InOut
                state["sigma_sqd_w"], state["sigma_sqd_w_zt"], state["radht"],                            # InOut
                state["deep_soil_T_in_K"], state["sfc_soil_T_in_K"], state["veg_T_in_K"],                 # InOut
                state["pdf_params"], state["pdf_params_zm"],                                              # InOut
                state["rcm_mc"], state["rvm_mc"], state["thlm_mc"],                                       # Out
                state["wprtp_mc"], state["wpthlp_mc"], state["rtp2_mc"],                                  # Out
                state["thlp2_mc"], state["rtpthlp_mc"],                                                   # Out
                state["wpthlp_sfc"], state["wprtp_sfc"], state["upwp_sfc"], state["vpwp_sfc"],            # Out
            )
        except Exception:
            # Close the host output handle when restart reads fail.
            if stats_writer is not None:
                stats_writer.finalize()
            raise

        # Calculate reciprocals from the dry densities read from the input file.
        state['invrs_rho_ds_zm'] = 1.0 / state['rho_ds_zm']
        state['invrs_rho_ds_zt'] = 1.0 / state['rho_ds_zt']

        if parameters_microphys.lh_microphys_type != parameters_microphys.lh_microphys_disabled:
            # JAX adaptation: reconstruct the case-owned permutation using the
            # seed of its last reshuffle. No model steps or physics are replayed.
            # The source reseeds each timestep but stores this permutation only
            # in memory; native JAX keys let initialization recover it exactly.
            from clubb_jax.src.SILHS.generate_uniform_sample_module import generate_uniform_lh_sample

            sequence_length = parameters_microphys.lh_sequence_length
            last_iter = iinit - 1
            if sequence_length > 1 and last_iter > 0:
                reshuffle_iter = ((last_iter - 1) // sequence_length) * sequence_length + 1
                lh_seed_custom = jnp.asarray(parameters_microphys.lh_seed * reshuffle_iter, dtype=jnp.uint32)
                key = jax.random.fold_in(jax.random.PRNGKey(lh_seed_custom), 1)
                key = jax.random.fold_in(key, ngrdcol - 1)
                draw_key, _ = jax.random.split(key)
                _, sampling_state = generate_uniform_lh_sample(
                    reshuffle_iter, parameters_microphys.lh_num_samples, sequence_length, pdf_dim + 2,  # In
                    silhs_config_flags.l_lh_deterministic_test,                            # In
                    draw_key,                                                              # In
                    state['sampling_state'],                                                              # InOut
                )
                state['sampling_state'] = LatinHypercubeArrays(
                    sampling_state.one_height_time_matrix,
                    jnp.asarray(last_iter, dtype=jnp.int32),
                )

    print(f"Initialized {runtype} case: nzm={nzm}, nzt={nzt}, ngrdcol={ngrdcol}")
    print(f"  dt_main={dt_main}s, time={time_initial}s to {time_final}s, {ifinal} steps")

    # Adaptation: immutable JAX leaves replace the source's per-field initial
    # allocations. Capture restored restart/SILHS state before any advance.
    state = jax.tree_util.tree_map(
        lambda value: jnp.asarray(value) if isinstance(value, np.ndarray) else value,
        state,
    )
    state['_initial_state'] = dict(state)
    return state


def set_case_initial_conditions(
    state: dict, clubb_params_in=None, batch_num=None,
):
    """Reset advanced fields and statistics for a rerun or runtime batch.

    Mirrors ``set_case_initial_conditions`` in ``src/clubb_driver.F90``.
    ``state`` owns the source driver globals and its InOut error container.
    The immutable initialization snapshot replaces repeated field assignments;
    radiation/PDF/microphysics caches and restart sampling state reset with it.

    Args:
        state: initialized driver fields and resources [InOut].
        clubb_params_in: optional replacement (ngrdcol, nparams) matrix [In].
        batch_num: optional one-based logical runtime batch selector [In].
    """
    initial_state = state['_initial_state']
    stats_writer = initial_state['stats_writer']
    if stats_writer is not None:
        stats_writer.reset()

    if batch_num is not None:
        num_batches = state['total_param_sets'] // state['ngrdcol']
        if batch_num < 1 or batch_num > num_batches:
            raise ValueError(
                'set_case_initial_conditions batch_num is outside the valid runtime batch range'
            )
        if batch_num > 1 and stats_writer is not None:
            stats_writer.start_next_batch()

    # Preserve newly supplied parameters or select an internally stored batch.
    # Without either optional input, the source retains the current matrix.
    if clubb_params_in is not None:
        if clubb_params_in.shape != initial_state['clubb_params'].shape:
            raise ValueError(
                'clubb_params_in must match the initialized (ngrdcol,nparams) shape'
            )
        clubb_params = clubb_params_in
    elif batch_num is not None:
        batch_start = (batch_num - 1) * state['ngrdcol']
        batch_end = batch_start + state['ngrdcol']
        clubb_params = initial_state['clubb_params_all'][batch_start:batch_end, :]
    else:
        clubb_params = state['clubb_params']

    # Restore case fields, zeroed diagnostics and the absolute restart clock.
    # Resource objects stay host-owned; statistics reset through their API above.
    state.clear()
    state.update(initial_state)
    state['_initial_state'] = initial_state
    state['clubb_params'] = jnp.asarray(
        clubb_params, dtype=initial_state['clubb_params'].dtype,
    )

    # Recompute parameter-dependent derived state for the next run.
    state['nu_vert_res_dep'], state['lmin'], state['mixt_frac_max_mag'] = calc_derived_params(
        state['gr'], state['ngrdcol'], state['cfg']['grid_type'],                    # In
        jnp.full((state['ngrdcol'],), state['cfg']['deltaz_nl']),                     # In
        state['clubb_params'], state['flags'].l_prescribed_avg_deltaz,               # In
    )
    if float(state['lmin']) < 1.0:
        raise ValueError('lmin is < 1.0')

    # Re-run the parameter sanity checks when debugging, as at initial setup.
    if clubb_at_least_debug_level(1):
        # TODO: the existing settings API is keyword-only; retain its contract
        # until its signature mirrors the source's explicit argument groups.
        state['err_info'] = check_clubb_settings(
            ngrdcol=state['ngrdcol'], params=state['clubb_params'],                  # In
            l_implemented=False, l_input_fields=False, config_flags=state['flags'], # In
            err_info=state['err_info'],                                            # InOut
        )
        if state['err_info'].is_fatal():
            raise RuntimeError(
                'Fatal error calling check_clubb_settings in set_case_initial_conditions'
            )
    return state['err_info']


def clean_up_clubb(state: dict):
    """Clean up stats state."""
    if state['l_stats'] and state.get('stats_writer') is not None:
        state['stats_writer'].finalize()
    if state['l_t_dependent']:
        time_dependent_input.finalize_t_dependent_input()
    from clubb_jax.src.Microphys.microphys_init_cleanup import cleanup_microphys
    cleanup_microphys()
    print("CLUBB cleanup complete.")


# -----------------------------------------------------------------------------
def restart_clubb(
    gr, hydromet_dim, hm_metadata,                                                         # In
    restart_path_case, time_restart,                                                       # In
    um, upwp, vm, vpwp, up2, vp2, rtm,                                                     # InOut
    wprtp, thlm, wpthlp, rtp2, rtp3,                                                       # InOut
    thlp2, thlp3, rtpthlp, wp2, wp3,                                                       # InOut
    p_in_Pa, exner, rcm, cloud_frac,                                                       # InOut
    wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp,                                             # InOut
    wp2rtp, wp2thlp, uprcp, vprcp,                                                         # InOut
    rc_coef_zm, wp4, wpup2, wpvp2, wp2up2,                                                 # InOut
    wp2vp2, ice_supersat_frac,                                                             # InOut
    wm_zt, rho, rho_zm, rho_ds_zm,                                                         # InOut
    rho_ds_zt, thv_ds_zm, thv_ds_zt,                                                       # InOut
    thlm_forcing, rtm_forcing, wprtp_forcing,                                              # InOut
    wpthlp_forcing, rtp2_forcing,                                                          # InOut
    thlp2_forcing, rtpthlp_forcing,                                                        # InOut
    hydromet, hydrometp2, wphydrometp,                                                     # InOut
    Ncm, Nccnm, thvm, em, tau_zm, tau_zt,                                                  # InOut
    Kh_zt, Kh_zm, ug, vg,                                                                  # InOut
    thlprcp,                                                                               # InOut
    sigma_sqd_w, sigma_sqd_w_zt, radht,                                                    # InOut
    deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K,                                         # InOut
    pdf_params, pdf_params_zm,                                                             # InOut
    rcm_mc, rvm_mc, thlm_mc,                                                               # Out
    wprtp_mc, wpthlp_mc, rtp2_mc,                                                          # Out
    thlp2_mc, rtpthlp_mc,                                                                  # Out
    wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,                                             # Out
):
    """Initialize CLUBB to a designated point in the submitted netCDF file.

    Source: clubb_driver.F90:restart_clubb. Inputs/returns retain source order;
    JAX profile arrays include the column axis, and PDF fields are immutable.
    NetCDF reads run entirely during host initialization.

    Arguments (source order; profiles include all columns):
        gr: CLUBB thermodynamic/momentum grid [m].
        hydromet_dim: Number of hydrometeor species [-].
        hm_metadata: Hydrometeor/PDF species and variable indices [-].
        restart_path_case: Path to netCDF data for restart
        time_restart: Time of model restart [s].
        um: eastward grid-mean wind component (thermo. levs.)  [m/s]
        upwp: u'w' (momentum levels)                         [m^2/s^2]
        vm: northward grid-mean wind component (thermo. levs.) [m/s]
        vpwp: v'w' (momentum levels)                         [m^2/s^2]
        up2: u'^2 (momentum levels)                         [m^2/s^2]
        vp2: v'^2 (momentum levels)                         [m^2/s^2]
        rtm: total water mixing ratio, r_t (thermo. levels) [kg/kg]
        wprtp: w' r_t' (momentum levels)                      [kg/kg m/s]
        thlm: liq. water pot. temp., th_l (thermo. levels)   [K]
        wpthlp: w'th_l' (momentum levels)                      [(m/s) K]
        rtp2: r_t'^2 (momentum levels)                       [(kg/kg)^2]
        rtp3: r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
        thlp2: th_l'^2 (momentum levels)                      [K^2]
        thlp3: th_l'^3 (thermodynamic levels)                 [K^3]
        rtpthlp: r_t'th_l' (momentum levels)                    [(kg/kg) K]
        wp2: w'^2 (momentum levels)                         [m^2/s^2]
        wp3: w'^3 (thermodynamic levels)                    [m^3/s^3]
        p_in_Pa: Air pressure (thermodynamic levels)            [Pa]
        exner: Exner function (thermodynamic levels)          [-]
        rcm: cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
        cloud_frac: cloud fraction (thermodynamic levels)          [-]
        wpthvp: < w' th_v' > (momentum levels)                 [kg/kg K]
        wp2thvp: < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]
        wp2up: < w'^2 u' > (thermodynamic levels)             [m^3/s^3]
        rtpthvp: < r_t' th_v' > (momentum levels)               [kg/kg K]
        thlpthvp: < th_l' th_v' > (momentum levels)              [K^2]
        wp2rtp: w'^2 rt' (thermodynamic levels)      [m^2/s^2 kg/kg]
        wp2thlp: w'^2 thl' (thermodynamic levels)     [m^2/s^2 K]
        uprcp: < u' r_c' > (momentum levels)        [(m/s)(kg/kg)]
        vprcp: < v' r_c' > (momentum levels)        [(m/s)(kg/kg)]
        rc_coef_zm: Coef of X'r_c' in Eq. (34) (m-levs.) [K/(kg/kg)]
        wp4: w'^4 (momentum levels)               [m^4/s^4]
        wpup2: w'u'^2 (thermodynamic levels)        [m^3/s^3]
        wpvp2: w'v'^2 (thermodynamic levels)        [m^3/s^3]
        wp2up2: w'^2 u'^2 (momentum levels)          [m^4/s^4]
        wp2vp2: w'^2 v'^2 (momentum levels)          [m^4/s^4]
        ice_supersat_frac: ice cloud fraction (thermo. levels)  [-]
        wm_zt: vertical mean wind component on thermo. levels  [m/s]
        rho: Air density on thermodynamic levels             [kg/m^3]
        rho_zm: Air density on momentum levels               [kg/m^3]
        rho_ds_zm: Dry, static density on momentum levels       [kg/m^3]
        rho_ds_zt: Dry, static density on thermo. levels           [kg/m^3]
        thv_ds_zm: Dry, base-state theta_v on momentum levels   [K]
        thv_ds_zt: Dry, base-state theta_v on thermo levels        [K]
        thlm_forcing: liquid potential temp. forcing (thermo. levels) [K/s]
        rtm_forcing: total water forcing (thermo. levels)      [(kg/kg)/s]
        wprtp_forcing: total water turbulent flux forcing (m-levs) [m*K/s^2]
        wpthlp_forcing: liq pot temp turb flux forcing (m-levs)[m(kg/kg)/s^2]
        rtp2_forcing: total water variance forcing (m-levs)   [(kg/kg)^2/s]
        thlp2_forcing: liq pot temp variance forcing (m-levs)  [K^2/s]
        rtpthlp_forcing: <r_t'th_l'> covariance forcing (m-levs) [K*(kg/kg)/s]
        hydromet: Array of hydrometeors                [hm units]
        hydrometp2: Variance of a hydrometeor (m-levs.)  [<hm units>^2]
        wphydrometp: Covariance of w and a hydrometeor    [(m/s) <hm units>]
        Ncm: Mean cloud droplet concentration, <N_c> (t-levs.)    [num/kg]
        Nccnm: Cloud condensation nuclei concentration (COAMPS/MG)  [num/kg]
        thvm: Virtual potential temperature                        [K]
        em: Turbulent Kinetic Energy (TKE)                       [m^2/s^2]
        tau_zm: Eddy dissipation time scale on momentum levels       [s]
        tau_zt: Eddy dissipation time scale on thermodynamic levels  [s]
        Kh_zt: Eddy diffusivity coefficient on thermodynamic levels [m^2/s]
        Kh_zm: Eddy diffusivity coefficient on momentum levels      [m^2/s]
        ug: u geostrophic wind                                   [m/s]
        vg: v geostrophic wind                                   [m/s]
        thlprcp: Inout
        sigma_sqd_w: PDF width parameter (momentum levels)                [-]
        sigma_sqd_w_zt: PDF width parameter interpolated to t-levs.          [-]
        radht: SW + LW heating rate                                 [K/s]
        deep_soil_T_in_K: Deep soil temperature [K].
        sfc_soil_T_in_K: Surface soil temperature [K].
        veg_T_in_K: Vegetation temperature [K].
        pdf_params: PDF parameters (thermodynamic levels)    [units vary]
        pdf_params_zm: PDF parameters on momentum levels        [units vary]
        rcm_mc: Tendency of liquid water due to microphysics      [kg/kg/s]
        rvm_mc: Tendency of vapor water due to microphysics       [kg/kg/s]
        thlm_mc: Tendency of liquid pot. temp. due to microphysics [K/s]
        wprtp_mc: Microphysics tendency for <w'rt'>   [m*(kg/kg)/s^2]
        wpthlp_mc: Microphysics tendency for <w'thl'>  [m*K/s^2]
        rtp2_mc: Microphysics tendency for <rt'^2>   [(kg/kg)^2/s]
        thlp2_mc: Microphysics tendency for <thl'^2>  [K^2/s]
        rtpthlp_mc: Microphysics tendency for <rt'thl'> [K*(kg/kg)/s]
        wpthlp_sfc: w'theta_l' surface flux   [(m K)/s]
        wprtp_sfc: w'rt' surface flux        [(m kg)/(kg s)]
        upwp_sfc: u'w' at surface           [m^2/s^2]
        vpwp_sfc: v'w' at surface           [m^2/s^2]
    """
    from clubb_jax.src.Input_fields import input_fields
    from clubb_jax.src.Microphys import parameters_microphys

    # Inform inputfields module.
    input_fields.l_input_um = True
    input_fields.l_input_vm = True
    input_fields.l_input_rtm = True
    input_fields.l_input_thlm = True
    input_fields.l_input_wp2 = True
    input_fields.l_input_ug = True
    input_fields.l_input_vg = True
    input_fields.l_input_rcm = True
    input_fields.l_input_wm_zt = True
    input_fields.l_input_exner = True
    input_fields.l_input_em = True
    input_fields.l_input_p = True
    input_fields.l_input_rho = True
    input_fields.l_input_rho_zm = True
    input_fields.l_input_rho_ds_zm = True
    input_fields.l_input_rho_ds_zt = True
    input_fields.l_input_thv_ds_zm = True
    input_fields.l_input_thv_ds_zt = True
    input_fields.l_input_Lscale = True
    input_fields.l_input_Lscale_up = True
    input_fields.l_input_Lscale_down = True
    input_fields.l_input_Kh_zt = True
    input_fields.l_input_Kh_zm = True
    input_fields.l_input_tau_zm = True
    input_fields.l_input_tau_zt = True
    input_fields.l_input_thvm = True
    input_fields.l_input_wpthvp = True
    input_fields.l_input_wp2thvp = True
    input_fields.l_input_wp2up = True
    input_fields.l_input_rtpthvp = True
    input_fields.l_input_thlpthvp = True
    input_fields.l_input_wp2rtp = True
    input_fields.l_input_wp2thlp = True
    input_fields.l_input_uprcp = True
    input_fields.l_input_vprcp = True
    input_fields.l_input_rc_coef_zm = True
    input_fields.l_input_wp4 = True
    input_fields.l_input_wpup2 = True
    input_fields.l_input_wpvp2 = True
    input_fields.l_input_wp2up2 = True
    input_fields.l_input_wp2vp2 = True
    input_fields.l_input_iss_frac = True
    input_fields.l_input_w_1 = True
    input_fields.l_input_w_2 = True
    input_fields.l_input_varnce_w_1 = True
    input_fields.l_input_varnce_w_2 = True
    input_fields.l_input_rt_1 = True
    input_fields.l_input_rt_2 = True
    input_fields.l_input_varnce_rt_1 = True
    input_fields.l_input_varnce_rt_2 = True
    input_fields.l_input_thl_1 = True
    input_fields.l_input_thl_2 = True
    input_fields.l_input_varnce_thl_1 = True
    input_fields.l_input_varnce_thl_2 = True
    input_fields.l_input_mixt_frac = True
    input_fields.l_input_chi_1 = True
    input_fields.l_input_chi_2 = True
    input_fields.l_input_stdev_chi_1 = True
    input_fields.l_input_stdev_chi_2 = True
    input_fields.l_input_rc_1 = True
    input_fields.l_input_rc_2 = True
    input_fields.l_input_w_1_zm = True
    input_fields.l_input_w_2_zm = True
    input_fields.l_input_varnce_w_1_zm = True
    input_fields.l_input_varnce_w_2_zm = True
    input_fields.l_input_mixt_frac_zm = True
    input_fields.l_input_radht = True

    microphys_scheme = parameters_microphys.microphys_scheme
    if microphys_scheme == "coamps":
        input_fields.l_input_rrm = True
        input_fields.l_input_rsm = True
        input_fields.l_input_rim = True
        input_fields.l_input_rgm = True
        input_fields.l_input_Nccnm = True
        input_fields.l_input_Ncm = True
        input_fields.l_input_Nrm = True
        input_fields.l_input_Nim = True

    elif microphys_scheme == "morrison":
        input_fields.l_input_rrm = True
        input_fields.l_input_Nrm = True
        if parameters_microphys.l_ice_microphys:
            input_fields.l_input_rsm = True
            input_fields.l_input_rim = True
            input_fields.l_input_Nim = True
            if parameters_microphys.l_graupel:
                input_fields.l_input_rgm = True
            else:
                input_fields.l_input_rgm = False
        else:
            input_fields.l_input_rsm = False
            input_fields.l_input_rim = False
            input_fields.l_input_Nim = False
            input_fields.l_input_rgm = False
        input_fields.l_input_Nccnm = False
        if parameters_microphys.l_predict_Nc:
            input_fields.l_input_Ncm = True
        else:
            input_fields.l_input_Ncm = False

    elif microphys_scheme == "khairoutdinov_kogan":
        input_fields.l_input_rrm = True
        input_fields.l_input_rsm = False
        input_fields.l_input_rim = False
        input_fields.l_input_rgm = False
        input_fields.l_input_Nccnm = False
        input_fields.l_input_Ncm = False
        input_fields.l_input_Nrm = True
        input_fields.l_input_Nim = False

    else:
        input_fields.l_input_rrm = False
        input_fields.l_input_rsm = False
        input_fields.l_input_rim = False
        input_fields.l_input_rgm = False
        input_fields.l_input_Nccnm = False
        input_fields.l_input_Ncm = False
        input_fields.l_input_Nrm = False
        input_fields.l_input_Nim = False

    # Source soil/vegetation flag supplied from the case-owned radiation config.
    input_fields.l_input_veg_T_in_K = input_fields.l_soil_veg
    input_fields.l_input_deep_soil_T_in_K = input_fields.l_soil_veg
    input_fields.l_input_sfc_soil_T_in_K = input_fields.l_soil_veg

    input_fields.l_input_wprtp = True
    input_fields.l_input_wpthlp = True
    input_fields.l_input_wp3 = True
    input_fields.l_input_rtp2 = True
    input_fields.l_input_rtp3 = True
    input_fields.l_input_thlp2 = True
    input_fields.l_input_thlp3 = True
    input_fields.l_input_rtpthlp = True
    input_fields.l_input_upwp = True
    input_fields.l_input_vpwp = True
    input_fields.l_input_thlm_forcing = True
    input_fields.l_input_rtm_forcing = True
    input_fields.l_input_up2 = True
    input_fields.l_input_vp2 = True
    input_fields.l_input_sigma_sqd_w = True
    input_fields.l_input_cloud_frac = True
    input_fields.l_input_sigma_sqd_w_zt = True
    input_fields.l_input_wprtp_forcing = True
    input_fields.l_input_wpthlp_forcing = True
    input_fields.l_input_rtp2_forcing = True
    input_fields.l_input_thlp2_forcing = True
    input_fields.l_input_rtpthlp_forcing = True
    input_fields.l_input_thlprcp = True
    input_fields.l_input_rcm_mc = True
    input_fields.l_input_rvm_mc = True
    input_fields.l_input_thlm_mc = True
    input_fields.l_input_wprtp_mc = True
    input_fields.l_input_wpthlp_mc = True
    input_fields.l_input_rtp2_mc = True
    input_fields.l_input_thlp2_mc = True
    input_fields.l_input_rtpthlp_mc = True

    stat_files = input_fields.set_filenames(restart_path_case)
    # Determine the nearest timestep in the netCDF file to the restart time.
    timestep = input_fields.compute_timestep(stat_files[0], True, time_restart)

    # Read data from stats files.
    (
        um, upwp, vm, vpwp,
        up2, vp2, rtm,
        wprtp, thlm, wpthlp,
        rtp2, rtp3,
        thlp2, thlp3, rtpthlp,
        wp2, wp3,
        p_in_Pa, exner, rcm, cloud_frac,
        wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp,
        wp2rtp, wp2thlp, uprcp, vprcp,
        rc_coef_zm, wp4, wpup2,
        wpvp2, wp2up2,
        wp2vp2, ice_supersat_frac,
        wm_zt, rho, rho_zm, rho_ds_zm,
        rho_ds_zt, thv_ds_zm, thv_ds_zt,
        thlm_forcing, rtm_forcing, wprtp_forcing,
        wpthlp_forcing, rtp2_forcing,
        thlp2_forcing, rtpthlp_forcing,
        hydromet, hydrometp2, wphydrometp,
        Ncm, Nccnm, thvm, em,
        tau_zm, tau_zt,
        Kh_zt, Kh_zm, ug, vg,
        thlprcp,
        sigma_sqd_w, sigma_sqd_w_zt, radht,
        deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K,
        pdf_params, pdf_params_zm
    ) = input_fields.stat_fields_reader(
        gr, timestep, hydromet_dim, hm_metadata,
        microphys_scheme, parameters_microphys.l_predict_Nc,
        um, upwp, vm, vpwp,
        up2, vp2, rtm,
        wprtp, thlm, wpthlp,
        rtp2, rtp3,
        thlp2, thlp3, rtpthlp,
        wp2, wp3,
        p_in_Pa, exner, rcm, cloud_frac,
        wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp,
        wp2rtp, wp2thlp, uprcp, vprcp,
        rc_coef_zm, wp4, wpup2,
        wpvp2, wp2up2,
        wp2vp2, ice_supersat_frac,
        wm_zt, rho, rho_zm, rho_ds_zm,
        rho_ds_zt, thv_ds_zm, thv_ds_zt,
        thlm_forcing, rtm_forcing, wprtp_forcing,
        wpthlp_forcing, rtp2_forcing,
        thlp2_forcing, rtpthlp_forcing,
        hydromet, hydrometp2, wphydrometp,
        Ncm, Nccnm, thvm, em,
        tau_zm, tau_zt,
        Kh_zt, Kh_zm, ug, vg,
        thlprcp,
        sigma_sqd_w, sigma_sqd_w_zt, radht,
        deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K,
        pdf_params, pdf_params_zm
    )

    rcm_mc, l_read_error = input_fields.get_clubb_variable_interpolated(
        input_fields.l_input_rcm_mc, stat_files[0], "rcm_mc", gr.nzt, timestep,            # In
        gr.zt[0, :],                                                                       # In
        rcm_mc,                                                                            # InOut
    )
    if l_read_error:
        raise ValueError("Failed to read rcm_mc for CLUBB restart")

    rvm_mc, l_read_error = input_fields.get_clubb_variable_interpolated(
        input_fields.l_input_rvm_mc, stat_files[0], "rvm_mc", gr.nzt, timestep,            # In
        gr.zt[0, :],                                                                       # In
        rvm_mc,                                                                            # InOut
    )
    if l_read_error:
        raise ValueError("Failed to read rvm_mc for CLUBB restart")

    thlm_mc, l_read_error = input_fields.get_clubb_variable_interpolated(
        input_fields.l_input_thlm_mc, stat_files[0], "thlm_mc", gr.nzt, timestep,          # In
        gr.zt[0, :],                                                                       # In
        thlm_mc,                                                                           # InOut
    )
    if l_read_error:
        raise ValueError("Failed to read thlm_mc for CLUBB restart")

    wprtp_mc, l_read_error = input_fields.get_clubb_variable_interpolated(
        input_fields.l_input_wprtp_mc, stat_files[1], "wprtp_mc", gr.nzm, timestep,        # In
        gr.zm[0, :],                                                                       # In
        wprtp_mc,                                                                          # InOut
    )
    if l_read_error:
        raise ValueError("Failed to read wprtp_mc for CLUBB restart")

    wpthlp_mc, l_read_error = input_fields.get_clubb_variable_interpolated(
        input_fields.l_input_wpthlp_mc, stat_files[1], "wpthlp_mc", gr.nzm, timestep,      # In
        gr.zm[0, :],                                                                       # In
        wpthlp_mc,                                                                         # InOut
    )
    if l_read_error:
        raise ValueError("Failed to read wpthlp_mc for CLUBB restart")

    rtp2_mc, l_read_error = input_fields.get_clubb_variable_interpolated(
        input_fields.l_input_rtp2_mc, stat_files[1], "rtp2_mc", gr.nzm, timestep,          # In
        gr.zm[0, :],                                                                       # In
        rtp2_mc,                                                                           # InOut
    )
    if l_read_error:
        raise ValueError("Failed to read rtp2_mc for CLUBB restart")

    thlp2_mc, l_read_error = input_fields.get_clubb_variable_interpolated(
        input_fields.l_input_thlp2_mc, stat_files[1], "thlp2_mc", gr.nzm, timestep,        # In
        gr.zm[0, :],                                                                       # In
        thlp2_mc,                                                                          # InOut
    )
    if l_read_error:
        raise ValueError("Failed to read thlp2_mc for CLUBB restart")

    rtpthlp_mc, l_read_error = input_fields.get_clubb_variable_interpolated(
        input_fields.l_input_rtpthlp_mc, stat_files[1], "rtpthlp_mc", gr.nzm, timestep,    # In
        gr.zm[0, :],                                                                       # In
        rtpthlp_mc,                                                                        # InOut
    )
    if l_read_error:
        raise ValueError("Failed to read rtpthlp_mc for CLUBB restart")

    wpthlp_sfc = wpthlp[:, 0]
    wprtp_sfc = wprtp[:, 0]
    upwp_sfc = upwp[:, 0]
    vpwp_sfc = vpwp[:, 0]

    return (
        um, upwp, vm, vpwp, up2, vp2, rtm,
        wprtp, thlm, wpthlp, rtp2, rtp3,
        thlp2, thlp3, rtpthlp, wp2, wp3,
        p_in_Pa, exner, rcm, cloud_frac,
        wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp,
        wp2rtp, wp2thlp, uprcp, vprcp,
        rc_coef_zm, wp4, wpup2, wpvp2, wp2up2,
        wp2vp2, ice_supersat_frac,
        wm_zt, rho, rho_zm, rho_ds_zm,
        rho_ds_zt, thv_ds_zm, thv_ds_zt,
        thlm_forcing, rtm_forcing, wprtp_forcing,
        wpthlp_forcing, rtp2_forcing,
        thlp2_forcing, rtpthlp_forcing,
        hydromet, hydrometp2, wphydrometp,
        Ncm, Nccnm, thvm, em, tau_zm, tau_zt,
        Kh_zt, Kh_zm, ug, vg,
        thlprcp,
        sigma_sqd_w, sigma_sqd_w_zt, radht,
        deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K,
        pdf_params, pdf_params_zm,
        rcm_mc, rvm_mc, thlm_mc,
        wprtp_mc, wpthlp_mc, rtp2_mc,
        thlp2_mc, rtpthlp_mc,
        wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc
    )
