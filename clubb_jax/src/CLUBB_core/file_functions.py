"""Read CLUBB input files: JAX port of ``src/CLUBB_core/file_functions.F90``.

Description:
  Reads numbers from plain data files and settings from Fortran namelist files.

Usage:
  file_read_1d(...), file_read_2d(...) -- read whitespace-separated values from
    a data file into a 1-D or 2-D array (used by the M-PACE A forcing
    initialization).
  read_namelist(path) -- read a CLUBB case namelist into one flat dict of
    settings, with the driver defaults for settings the file leaves out.
  _read_namelist_groups(path) -- read a namelist file into a dict per namelist
    group, for callers that need one group, such as the tunable parameters.

Understanding:
  The Fortran reads namelists with the language's own ``read(nml=...)``
  statement in each module. Python has no such statement, so the namelist
  reader lives here with the other file-reading routines. It uses ``f90nml``
  when it is installed, and otherwise a small parser for the namelist syntax
  that CLUBB input files use.
"""

from __future__ import annotations

import re
import sys

import numpy as np
try:
    import f90nml  # type: ignore
except ModuleNotFoundError:
    f90nml = None


# ===============================================================================
def file_read_1d(file_unit, path_and_filename, num_datapts, entries_per_line):
    """This subroutine reads in values from a data file with a number of
    rows and a declared number of columns (entries_per_line) of data.
    It reads in the data in the form of:
    1 ==> (row 1, column 1); 2 ==> (row 1, column 2); etc.

    Example:  a diagram of a data file with 18 total data points
              (DP1 to DP18), with 4 data points per row.

             i = 1     i = 2     i = 3     i = 4
          ---------------------------------------
    k = 1 |   DP1       DP2       DP3       DP4
          |
    k = 2 |   DP5       DP6       DP7       DP8
          |
    k = 3 |   DP9       DP10      DP11      DP12
          |
    k = 4 |   DP13      DP14      DP15      DP16
          |
    k = 5 |   DP17      DP18

    See Michael Falk's comments below for more information.

    Returns:
        variable

    Fortran declarations (file_read_1d in src/CLUBB_core/file_functions.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
         file_unit,          & ! Unit number of file being read.
         num_datapts,        & ! Total number of data points being read in.
         entries_per_line   ! Number of data points
        ! on one line of the file being read.

        character(*), intent(in) :: &
         path_and_filename  ! Path to file and filename of file being read.

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(num_datapts), intent(out) :: &
         variable           ! Data values output into variable
    """
    # ---- Begin Code ----
    # Python path-based reads do not use the Fortran unit number. NumPy handles
    # both the full and partial final lines as whitespace-separated input.
    del file_unit, entries_per_line

    # Open data file.
    try:
        variable = np.fromfile(path_and_filename, sep=" ", count=num_datapts)
    except OSError:
        print("CLUBB encountered an error trying to open " + str(path_and_filename), file=sys.stderr)
        raise
    if variable.size != num_datapts:
        raise ValueError(
            f"Expected {num_datapts} values in {path_and_filename}, "
            f"found {variable.size}"
        )
    return variable


# ===============================================================================
def file_read_2d(
    device,
    file_path,
    file_dimension1,
    file_dimension2,
    file_per_line,
):
    """Michael Falk wrote this routine to read data files in a particular format for mpace_a.
    The 2d mpace_a files list the (file_dimension2) values on a given vertical level, then
    moves to the next level to list its values.
    Each line has a specific number of values, until the last line on a level, which
    is short-- it has the last few values and then a line break.  The next line, beginning
    the next level, is full-sized again.  24 September 2007

    Fortran declarations (file_read_2d in src/CLUBB_core/file_functions.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          device, &
          file_dimension1, &
          file_dimension2, &
          file_per_line

        character(*), intent(in) :: &
          file_path

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(file_dimension1,file_dimension2), intent(out) :: &
          variable
    """
    # ---- Begin Code ----
    # Python path-based reads do not use Fortran units or per-line bookkeeping.
    del device, file_per_line
    num_datapts = file_dimension1 * file_dimension2
    variable = np.fromfile(file_path, sep=" ", count=num_datapts)
    if variable.size != num_datapts:
        raise ValueError(
            f"Expected {num_datapts} values in {file_path}, found {variable.size}"
        )
    return variable.reshape(file_dimension1, file_dimension2)


# ===============================================================================
# ── CLUBB default values (mirrors clubb_driver.F90 init_clubb_case) ──────

_DEFAULTS = dict(
    # model_setting
    runtype="bomex",
    sounding_file="",
    sclr_sounding_file="",
    edsclr_sounding_file="",
    ozone_sounding_file="",
    forcing_file="",
    sfc_file="",
    corr_file_cloud="",
    corr_file_below="",
    iop_file="",
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
            parts.append(token)
            buf = []
            continue
        buf.append(char)

    token = "".join(buf).strip()
    parts.append(token)
    return parts


def _parse_scalar(text: str):
    """Convert a Fortran-style literal into a Python scalar."""
    token = text.strip()
    lower = token.lower()

    # A null namelist entry leaves the corresponding initialized value unchanged.
    if not token:
        return None

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
