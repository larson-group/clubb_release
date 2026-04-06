"""CLUBB case initialization and cleanup utilities."""
import math
from pathlib import Path

import numpy as np

# Python API
from clubb_python import clubb_api
from clubb_python.derived_types.config_flags import ConfigFlags
from clubb_python.derived_types.grid_class import setup_grid as py_setup_grid
from clubb_python.derived_types.sclr_idx import SclrIdx
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.pdf_params import (
    init_pdf_implicit_coefs_terms_api,
    init_pdf_params as init_pdf_params_py,
)

# I/O
from clubb_python_driver.io.grid_file import read_grid_file
from clubb_python_driver.io.namelist import read_namelist
from clubb_python_driver.io.sounding import (
    read_sounding,
    interpolate_sounding,
    read_scalar_sounding,
    interpolate_scalar_sounding,
)
from clubb_python_driver.io.surface import read_surface

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
                                 fcor_y: np.ndarray, l_tke_aniso: bool,
                                 um: np.ndarray):
    """Mirror initialize_clubb() em/wp2/up2/vp2/upwp initialization."""
    em, um_adj = _initialize_em_profile(runtype, gr, um)

    if l_tke_aniso:
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
    else:
        wp2 = (2.0 / 3.0) * em
        up2 = np.zeros_like(em)
        vp2 = np.zeros_like(em)
        upwp = np.zeros_like(em)

    return em, wp2, up2, vp2, upwp, um_adj


def run_clubb(namelist_path: str, l_stdout: bool = True):
    """Run CLUBB standalone for a case described by a namelist file.

    Args:
        namelist_path: path to *_model.in file
        l_stdout: print timestep info to stdout
    """
    from clubb_python_driver.advance_clubb_to_end import advance_clubb_to_end

    state = init_clubb_case(namelist_path)
    advance_clubb_to_end(state, l_stdout=l_stdout)
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

    repo_root = Path(__file__).resolve().parents[1]
    return repo_root / "input" / "stats" / "standard_stats.in"


def _resolve_case_input_path(namelist_dir: Path, runtype: str, suffix: str) -> Path:
    """Resolve case input files for either model.in or aggregated CASE.in runs."""
    candidate = namelist_dir / f"{runtype}{suffix}"
    if candidate.exists():
        return candidate

    repo_root = Path(__file__).resolve().parents[1]
    fallback = repo_root / "input" / "case_setups" / f"{runtype}{suffix}"
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
    repo_root = Path(__file__).resolve().parents[1]
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
    """Check for namelist settings that the Python driver does not support.

    Raises ValueError with a clear message listing all unsupported features
    that are enabled, so the user can fix them all at once.
    """
    errors = []

    # --- Microphysics ---
    if microphys_scheme != "none":
        errors.append(
            f"microphys_scheme = '{microphys_scheme}' is not supported "
            "(only 'none' is implemented)."
        )

    # --- Cloud water sedimentation ---
    # cloud_drop_sed is called from the Fortran driver loop but not from
    # the Python driver; rcm_mc / thlm_mc stay zero when this is enabled.
    if bool(cfg.get('l_cloud_sed', False)):
        errors.append(
            "l_cloud_sed = true is not supported "
            "(cloud_drop_sed is not called from the Python driver)."
        )

    # --- Radiation ---
    supported_rad = {"none", "simplified", "simplified_bomex"}
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
    # Sponge damping is applied in the Fortran driver loop after
    # advance_clubb_core; the Python driver does not call these routines.
    _SPONGE_FIELDS = [
        "thlm", "rtm", "uv", "wp2", "wp3", "up2_vp2",
    ]
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
    if lh_type != "disabled":
        errors.append(
            f"lh_microphys_type = '{lh_type}' is not supported "
            "(SILHS sampling is not implemented in the Python driver)."
        )
    if bool(cfg.get('l_silhs_rad', False)):
        errors.append("l_silhs_rad = true is not supported (SILHS is not available).")

    # --- Soil / vegetation ---
    if bool(cfg.get('l_soil_veg', False)):
        errors.append("l_soil_veg = true is not supported.")

    # --- Restarts ---
    if bool(cfg.get('l_restart', False)):
        errors.append("l_restart = true is not supported (no GrADS restart I/O).")

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
    ngrdcol = cfg['ngrdcol']
    nzmax = cfg['nzmax']
    grid_type = cfg['grid_type']
    dt_main = cfg['dt_main']
    dt_rad = cfg['dt_rad']
    runtype = cfg['runtype']
    sclr_dim = cfg['sclr_dim']
    edsclr_dim = cfg['edsclr_dim']

    # ── 1. Initialize error handling ────────────────────────────────────
    clubb_api.init_err_info(ngrdcol)
    clubb_api.set_debug_level(cfg['debug_level'])

    # ── 2. Get config flags ─────────────────────────────────────────────
    flags = clubb_api.get_default_config_flags()
    # Override from namelist (configurable_clubb_flags_nl)
    flag_overrides = {}
    for name in ConfigFlags._fields:
        if name.lower() in cfg:
            flag_overrides[name] = cfg[name.lower()]
    if flag_overrides:
        d = flags._asdict()
        d.update(flag_overrides)
        flags = ConfigFlags(**d)
    clubb_api.init_config_flags(flags)

    saturation_formula = flags.saturation_formula
    microphys_scheme = str(cfg.get('microphys_scheme', 'none')).strip().strip("'\"").lower()
    l_cloud_sed = bool(cfg.get('l_cloud_sed', False))
    sigma_g = float(cfg.get('sigma_g', 1.5))
    rad_scheme = str(cfg.get('rad_scheme', 'none')).strip().strip("'\"").lower()
    l_calc_thlp2_rad = bool(cfg.get('l_calc_thlp2_rad', flags.l_calc_thlp2_rad))

    _check_unsupported_features(cfg, flags, microphys_scheme, rad_scheme, l_calc_thlp2_rad)

    # Initialize Fortran radiation-module state used by sunray_sw.
    clubb_api.set_simplified_radiation_params(
        f0=float(cfg.get('f0', 0.0)),
        f1=float(cfg.get('f1', 0.0)),
        kappa=float(cfg.get('kappa', 0.0)),
        eff_drop_radius=float(cfg.get('eff_drop_radius', 1.0e-5)),
        alvdr=float(cfg.get('alvdr', 0.1)),
        gc=float(cfg.get('gc', 0.85)),
        omega=float(cfg.get('omega', 0.992)),
        l_rad_above_cloud=bool(cfg.get('l_rad_above_cloud', False)),
        l_sw_radiation=bool(cfg.get('l_sw_radiation', False)),
        l_fix_cos_solar_zen=bool(cfg.get('l_fix_cos_solar_zen', False)),
        fs_values=cfg.get('fs_values', [0.0]),
        cos_solar_zen_values=cfg.get('cos_solar_zen_values', [-999.0]),
        cos_solar_zen_times=cfg.get('cos_solar_zen_times', [-999.0]),
    )

    # ── 3. Read sounding ────────────────────────────────────────────────
    snd_path = _resolve_case_input_path(namelist_dir, runtype, "_sounding.in")
    snd = read_sounding(str(snd_path))
    theta_type = snd['theta_type']
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
    result = clubb_api.hydrostatic(
        gr=gr, ngrdcol=ngrdcol, nzt=nzt, nzm=nzm, thvm=thvm, p_sfc=p_sfc
    )
    p_in_Pa, p_in_Pa_zm, exner, exner_zm, rho, rho_zm = result

    # Convert temperature type
    if theta_type in ('thm[K]', 'T[K]'):
        # theta sounding — need to compute rcm and convert to thlm
        thm = thlm.copy()
        # rcm = max(rtm - rsat(p, T), 0)
        T_in_K = thm * exner
        rsat = clubb_api.sat_mixrat_liq(
            gr=gr, nz=nzt, ngrdcol=ngrdcol, p_in_Pa=p_in_Pa, T_in_K=T_in_K,
            saturation_formula=saturation_formula
        )
        rcm = np.maximum(rtm - rsat, 0.0)
        thlm = thm - Lv / (Cp * exner) * rcm
    elif theta_type == 'thlm[K]':
        # Already liquid potential temperature
        rcm = np.zeros_like(thlm)
        for k in range(nzt):
            for i in range(ngrdcol):
                rcm[i, k] = clubb_api.rcm_sat_adj(
                    thlm[i, k], rtm[i, k], p_in_Pa[i, k], exner[i, k],
                    saturation_formula)
        thm = thlm + Lv / (Cp * exner) * rcm
    else:
        raise ValueError(f"Unknown theta_type: {theta_type}")

    # Recompute thvm and hydrostatic with corrected thlm
    # NOTE: Fortran passes thm (not thlm) as the 5th arg to calculate_thvm
    thvm = clubb_api.calculate_thvm(
        nzt=nzt, ngrdcol=ngrdcol, thlm=thlm, rtm=rtm, rcm=rcm, exner=exner,
        thv_ds_zt=thm * (1.0 + ep2 * (rtm - rcm))**kappa,
    )
    result = clubb_api.hydrostatic(
        gr=gr, ngrdcol=ngrdcol, nzt=nzt, nzm=nzm, thvm=thvm, p_sfc=p_sfc
    )
    p_in_Pa, p_in_Pa_zm, exner, exner_zm, rho, rho_zm = result

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
    rv_zm = clubb_api.zt2zm(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azt=rv)
    rv_zm = np.maximum(rv_zm, 0.0)
    thm_zm = clubb_api.zt2zm(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azt=thm)

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

    wm_zm = clubb_api.zt2zm(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azt=wm_zt)
    wm_zm[:, 0] = 0.0
    wm_zm[:, -1] = 0.0

    # ── 8. Initialize PDF and tunable parameters ────────────────────────
    clubb_params = clubb_api.init_clubb_params(ngrdcol, iunit=10, filename=namelist_path)
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
    clubb_api.set_sclr_idx(
        iisclr_rt, iisclr_thl, iisclr_co2,
        iiedsclr_rt, iiedsclr_thl, iiedsclr_co2,
    )

    nu_vert_res_dep, lmin, mixt_frac_max_mag = clubb_api.calc_derrived_params(
        gr=gr,
        ngrdcol=ngrdcol,
        grid_type=grid_type,
        deltaz=deltaz,
        clubb_params=clubb_params,
        nu_vert_res_dep=None,
        l_prescribed_avg_deltaz=False,
    )
    err_info = ErrInfo(ngrdcol=ngrdcol)

    err_info = clubb_api.check_clubb_settings(
        ngrdcol=ngrdcol,
        params=clubb_params,
        config_flags=flags,
        err_info=err_info,
        l_implemented=False,
        l_input_fields=False,
    )
    err_info = clubb_api.check_parameters(
        ngrdcol=ngrdcol,
        clubb_params=clubb_params,
        lmin=lmin,
        err_info=err_info,
    )
    # Reset error codes set by check_clubb_settings warnings
    # (they are non-fatal but would cause advance_clubb_core to bail out)
    clubb_api.reset_err_code()

    # ── 9. Initialize TKE / variances (case-specific, Fortran-like) ────
    em, wp2, up2, vp2, upwp, um = _initialize_turbulence_state(
        runtype=runtype,
        gr=gr,
        dt_main=dt_main,
        fcor_y=fcor_y,
        l_tke_aniso=bool(flags.l_tke_aniso),
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

    # Transport scalar/hydromet arrays across F2PY with a padded trailing extent.
    # The logical *_dim values remain authoritative and are used inside Fortran.
    hydromet_dim = 0
    hm_dim_transport = max(hydromet_dim, 1)
    l_mix_rat_hm = np.zeros((hm_dim_transport,), dtype=bool)
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
    sclrm_forcing = np.zeros((ngrdcol, nzt, sc_dim_transport))
    wpsclrp_sfc = np.zeros((ngrdcol, sc_dim_transport))

    edsclrm = np.zeros((ngrdcol, nzt, edsc_dim_transport))
    edsclrm_forcing = np.zeros((ngrdcol, nzt, edsc_dim_transport))
    wpedsclrp_sfc = np.zeros((ngrdcol, edsc_dim_transport))

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
    stats_nsamp = int(round(cfg['stats_tsamp'] / dt_main))
    stats_nout = int(round(cfg['stats_tout'] / dt_main))

    # ── 13. Initialize stats ────────────────────────────────────────────
    l_stats = bool(cfg['l_stats'])
    repo_root = Path(__file__).resolve().parents[1]
    stats_registry_path = _resolve_stats_registry_path(namelist_path, cfg)
    stats_prefix = str(cfg.get('fname_prefix', '')).strip() or runtype
    output_dir_raw = str(cfg.get("output_dir", "")).strip().strip("'\"")
    if output_dir_raw:
        output_dir_path = Path(output_dir_raw)
        if not output_dir_path.is_absolute():
            output_dir_path = (namelist_dir / output_dir_path).resolve()
    else:
        output_dir_path = repo_root / "output"
    stats_output_path = output_dir_path / f"{stats_prefix}_stats.nc"

    if l_stats:
        if not stats_registry_path.exists():
            raise FileNotFoundError(f"Stats registry file not found: {stats_registry_path}")
        stats_output_path.parent.mkdir(parents=True, exist_ok=True)
        err_info = clubb_api.init_stats(
            registry_path=str(stats_registry_path),
            output_path=str(stats_output_path),
            ncol=ngrdcol,
            stats_tsamp=float(cfg['stats_tsamp']),
            stats_tout=float(cfg['stats_tout']),
            dt_main=float(dt_main),
            day_in=int(cfg['day']),
            month_in=int(cfg['month']),
            year_in=int(cfg['year']),
            time_initial=float(time_initial),
            nzt=nzt,
            zt=gr.zt[0, :],
            nzm=nzm,
            zm=gr.zm[0, :],
            clubb_params=clubb_params,
            param_names=clubb_api.get_param_names(),
            err_info=err_info,
            sclr_dim=sclr_dim,
            edsclr_dim=edsclr_dim,
        )
        err = clubb_api.get_err_code(ngrdcol)
        if np.any(err != 0):
            raise RuntimeError(f"stats_init failed with err_code={err}")
        stats_enabled = bool(clubb_api.get_stats_config()[0])
        if not stats_enabled:
            raise RuntimeError("stats_init completed but stats are not enabled")

    # ── 14. Zero PDF params ─────────────────────────────────────────────
    pdf_params = init_pdf_params_py(nzt, ngrdcol)
    pdf_params_zm = init_pdf_params_py(nzm, ngrdcol)

    # ── 15. Clear any accumulated error codes from init ──────────────
    clubb_api.reset_err_code()

    # ── Build state dict ────────────────────────────────────────────────
    state = dict(
        # Config
        cfg=cfg, flags=flags, gr=gr, namelist_dir=str(namelist_dir),
        runtype=runtype, ngrdcol=ngrdcol, nzt=nzt, nzm=nzm,
        dt_main=dt_main, dt_rad=dt_rad,
        time_initial=time_initial, time_final=time_final,
        ifinal=ifinal, l_stats=l_stats,
        stats_nsamp=stats_nsamp, stats_nout=stats_nout,
        stats_registry_path=str(stats_registry_path),
        stats_output_path=str(stats_output_path),
        saturation_formula=saturation_formula,
        sfctype=int(cfg['sfctype']),
        microphys_scheme=microphys_scheme,
        l_cloud_sed=l_cloud_sed,
        sigma_g=sigma_g,
        nc0_in_cloud=nc0_in_cloud,
        rad_scheme=rad_scheme,
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
        thvm=thvm, radht=np.zeros((ngrdcol, nzt)),
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

    print(f"Initialized {runtype} case: nzm={nzm}, nzt={nzt}, ngrdcol={ngrdcol}")
    print(f"  dt_main={dt_main}s, time={time_initial}s to {time_final}s, {ifinal} steps")

    return state

def clean_up_clubb(state: dict):
    """Clean up Fortran state."""
    if state['l_stats']:
        state['err_info'] = clubb_api.finalize_stats(err_info=state['err_info'])
    clubb_api.cleanup_grid(gr=state['gr'])
    clubb_api.cleanup_err_info(err_info=state['err_info'])
    print("CLUBB cleanup complete.")
