"""Read CLUBB namelist files (clubb.in / *_model.in).

Prefer ``f90nml`` when it is available. Fall back to a simple local parser
that is sufficient for CLUBB's generated aggregate namelists.
"""
from __future__ import annotations

import re

try:
    import f90nml  # type: ignore
except ModuleNotFoundError:
    f90nml = None


# ── CLUBB default values (mirrors clubb_driver.F90 init_clubb_case) ──────

_DEFAULTS = dict(
    # model_setting
    runtype="bomex",
    nzmax=75,
    grid_type=1,
    deltaz_nl=40.0,
    zm_init_nl=0.0,
    zm_top_nl=3500.0,
    zt_grid_fname="",
    zm_grid_fname="",
    day=22, month=6, year=1969,
    lat_vals=15.0, lon_vals=-56.5,
    sfc_elevation_nl=0.0,
    time_initial=0.0,
    time_final=21600.0,
    dt_main=60.0,
    dt_rad=600.0,
    sfctype=0,
    t_sfc_nl=288.0,
    p_sfc_nl=1000.0e2,
    sens_ht=0.0,
    latent_ht=0.0,
    fcor_nl=1.0e-4,
    t0=300.0,
    ts_nudge=86400.0,
    rtm_min=1.0e-300,  # epsilon(rtm_min) in Fortran
    rtm_nudge_max_altitude=10000.0,
    forcings_file_path="",
    l_t_dependent=False,
    l_input_xpwp_sfc=False,
    l_ignore_forcings=False,
    l_modify_ic_with_cubic_int=False,
    l_modify_bc_for_cnvg_test=False,
    l_restart=False,
    restart_path_case="none",
    time_restart=0.0,
    l_input_fields=False,
    debug_level=2,
    sclr_dim=0,
    edsclr_dim=0,
    iisclr_thl=-1,
    iisclr_rt=-1,
    iisclr_co2=-1,
    iiedsclr_thl=-1,
    iiedsclr_rt=-1,
    iiedsclr_co2=-1,
    sclr_tol_nl=[],
    l_rtm_nudge=False,
    l_diagnose_correlations=False,
    l_calc_w_corr=False,
    rad_scheme="none",
    # multicol_def
    ngrdcol=1,
    # stats_setting
    l_stats=False,
    fname_prefix="",
    stats_tsamp=60.0,
    stats_tout=60.0,
    stats_fmt="netcdf",
    l_allow_small_stats_tout=False,
    output_dir="../output/",
)


def read_namelist(path: str) -> dict:
    """Read a CLUBB namelist file and return a flat dict of settings.

    All keys are lowercased. Defaults are applied for missing keys.
    Sponge damping settings are returned as nested dicts under their
    original namelist variable names.
    """
    nml = _read_namelist_groups(path)

    # Merge all namelist groups into one flat dict
    cfg = dict(_DEFAULTS)
    for group_name, group_data in nml.items():
        for key, value in group_data.items():
            cfg[key.lower()] = value

    return cfg


def _read_namelist_groups(path: str) -> dict:
    """Read namelist groups with f90nml when available, else a simple parser."""
    if f90nml is not None:
        return f90nml.read(path)
    return _parse_namelist_text(path)


def _parse_namelist_text(path: str) -> dict:
    """Parse the subset of Fortran namelist syntax used by CLUBB input files."""
    groups: dict[str, dict] = {}
    current_group: dict | None = None

    with open(path, encoding="utf-8") as handle:
        for raw_line in handle:
            line = re.sub(r"!.*", "", raw_line).strip()
            if not line:
                continue
            if line.startswith("&"):
                group_name = line[1:].strip().lower()
                current_group = {}
                groups[group_name] = current_group
                continue
            if line == "/":
                current_group = None
                continue
            if current_group is None or "=" not in line:
                continue

            key, raw_value = line.split("=", 1)
            current_group[key.strip().lower()] = _parse_value(raw_value.strip())

    return groups


def _parse_value(raw_value: str):
    """Parse a scalar or comma-separated list from a namelist assignment."""
    value = raw_value.rstrip(",").strip()
    parts = _split_top_level_commas(value)
    if len(parts) > 1:
        return [_parse_scalar(part) for part in parts]
    return _parse_scalar(value)


def _split_top_level_commas(text: str) -> list[str]:
    """Split a value list on commas while respecting quoted strings."""
    parts: list[str] = []
    buf: list[str] = []
    quote: str | None = None

    for char in text:
        if quote is None and char in {"'", '"'}:
            quote = char
        elif quote == char:
            quote = None

        if char == "," and quote is None:
            token = "".join(buf).strip()
            if token:
                parts.append(token)
            buf = []
            continue
        buf.append(char)

    token = "".join(buf).strip()
    if token:
        parts.append(token)
    return parts


def _parse_scalar(text: str):
    """Convert a Fortran-style literal into a Python scalar."""
    token = text.strip()
    lower = token.lower()

    if lower in {".true.", "true"}:
        return True
    if lower in {".false.", "false"}:
        return False
    if len(token) >= 2 and token[0] == token[-1] and token[0] in {"'", '"'}:
        return token[1:-1]

    numeric = token.replace("d", "e").replace("D", "E")
    try:
        if re.fullmatch(r"[+-]?\d+", numeric):
            return int(numeric)
        return float(numeric)
    except ValueError:
        return token
