"""Host-side replacement for the Fortran stats NetCDF output infrastructure.

Together with ``jax_stats.JaxStats``, this mirrors ``stats_netcdf.F90``:
  - begin_timestep(itime)  → set sampling/output cadence flags
  - JaxStats updates       → accumulate values and budgets in compiled code
  - end_timestep(t, stats) → transfer, validate, average, and write one record
  - finalize()             → close the NetCDF file

Variable shapes in the Python/JAX codebase are (ngrdcol, nz) for profile
variables and (ngrdcol,) for surface variables.  The NetCDF dimensions in the
output file are (time, nz, col) for profiles and (time, col) for surface
variables, matching the Fortran-produced reference files.
"""

from __future__ import annotations

import itertools
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import numpy as np

from clubb_jax.src.CLUBB_core.stats_metadata import GRID_ORDER, StatsLayout

try:
    import netCDF4 as nc
    _HAS_NC = True
except ImportError:
    _HAS_NC = False


# Grid type identifiers (matching Fortran GRID_* enum)
GRID_ZT = "zt"
GRID_ZM = "zm"
GRID_SFC = "sfc"
GRID_LH_ZT = "lh_zt"
GRID_LH_SFC = "lh_sfc"
GRID_RAD_ZT = "rad_zt"
GRID_RAD_ZM = "rad_zm"

_GRID_IDS = {
    GRID_ZT: 1,
    GRID_ZM: 2,
    GRID_RAD_ZT: 3,
    GRID_RAD_ZM: 4,
    GRID_SFC: 5,
    GRID_LH_ZT: 6,
    GRID_LH_SFC: 7,
}

_PROFILE_GRIDS = {GRID_ZT, GRID_ZM, GRID_LH_ZT, GRID_RAD_ZT, GRID_RAD_ZM}
_VALID_GRIDS = _PROFILE_GRIDS | {GRID_SFC, GRID_LH_SFC}


@dataclass(frozen=True)
class StatsDefinition:
    """Static metadata for one expanded stats-registry entry."""

    name: str
    grid: str
    units: str
    long_name: str


def _hydromet_codes(hydromet_list: Iterable[str]) -> tuple[str, ...]:
    return tuple(str(name).strip()[:2] for name in hydromet_list if str(name).strip())


def _expanded_names(
    name: str,
    *,
    sclr_dim: int,
    edsclr_dim: int,
    hydromet_list: Iterable[str],
) -> tuple[str, ...]:
    """Mirror ``stats_expand_registry`` from ``stats_netcdf.F90``."""
    hm = _hydromet_codes(hydromet_list)

    paired_hm = {
        "hm_i": ("", ""),
        "mu_hm_i": ("mu_", ""),
        "sigma_hm_i": ("sigma_", ""),
        "corr_w_hm_i": ("corr_w_", ""),
        "corr_chi_hm_i": ("corr_chi_", ""),
        "corr_eta_hm_i": ("corr_eta_", ""),
        "corr_Ncn_hm_i": ("corr_Ncn_", ""),
        "mu_hm_i_n": ("mu_", "_n"),
        "sigma_hm_i_n": ("sigma_", "_n"),
        "corr_w_hm_i_n": ("corr_w_", "_n"),
        "corr_chi_hm_i_n": ("corr_chi_", "_n"),
        "corr_eta_hm_i_n": ("corr_eta_", "_n"),
        "corr_Ncn_hm_i_n": ("corr_Ncn_", "_n"),
    }
    if name in paired_hm:
        prefix, suffix = paired_hm[name]
        return tuple(
            f"{prefix}{code}_{component}{suffix}"
            for code in hm
            for component in (1, 2)
        )

    paired_ncn = {
        "mu_Ncn_i": "mu_Ncn",
        "sigma_Ncn_i": "sigma_Ncn",
        "corr_w_Ncn_i": "corr_w_Ncn",
        "corr_chi_Ncn_i": "corr_chi_Ncn",
        "corr_eta_Ncn_i": "corr_eta_Ncn",
        "mu_Ncn_i_n": "mu_Ncn",
        "sigma_Ncn_i_n": "sigma_Ncn",
        "corr_w_Ncn_i_n": "corr_w_Ncn",
        "corr_chi_Ncn_i_n": "corr_chi_Ncn",
        "corr_eta_Ncn_i_n": "corr_eta_Ncn",
    }
    if name in paired_ncn:
        suffix = "_n" if name.endswith("_n") else ""
        return tuple(f"{paired_ncn[name]}_{component}{suffix}" for component in (1, 2))

    if name in {"corr_hmx_hmy_i", "corr_hmx_hmy_i_n"}:
        suffix = "_n" if name.endswith("_n") else ""
        return tuple(
            f"corr_{left}_{right}_{component}{suffix}"
            for left, right in itertools.combinations(hm, 2)
            for component in (1, 2)
        )

    hydromet_single = {
        "wp2hmp": lambda code: f"wp2{code}p",
        "K_hm": lambda code: f"K_hm_{code}",
        "hydrometp2": lambda code: f"{code}p2",
        "wphydrometp": lambda code: f"wp{code}p",
        "hmp2_zt": lambda code: f"{code}p2_zt",
        "rtphmp": lambda code: f"rtp{code}p",
        "thlphmp": lambda code: f"thlp{code}p",
    }
    if name in hydromet_single:
        return tuple(hydromet_single[name](code) for code in hm)
    if name == "hmxphmyp":
        return tuple(f"{left}p{right}p" for left, right in itertools.combinations(hm, 2))

    scalar_templates = {
        "sclrm": "sclr{i}m",
        "sclrm_f": "sclr{i}m_f",
        "sclrprtp": "sclr{i}prtp",
        "sclrp2": "sclr{i}p2",
        "sclrpthvp": "sclr{i}pthvp",
        "sclrpthlp": "sclr{i}pthlp",
        "sclrprcp": "sclr{i}prcp",
        "wpsclrp": "wpsclr{i}p",
        "wpsclrp2": "wpsclr{i}p2",
        "wp2sclrp": "wp2sclr{i}p",
        "wpsclrprtp": "wpsclr{i}prtp",
        "wpsclrpthlp": "wpsclr{i}pthlp",
    }
    if name in scalar_templates:
        return tuple(scalar_templates[name].format(i=i) for i in range(1, sclr_dim + 1))

    ed_scalar_templates = {
        "edsclrm": "edsclr{i}m",
        "edsclrm_f": "edsclr{i}m_f",
        "wpedsclrp": "wpedsclr{i}p",
    }
    if name in ed_scalar_templates:
        return tuple(ed_scalar_templates[name].format(i=i) for i in range(1, edsclr_dim + 1))

    if name == "silhs_variance_category":
        return tuple(f"silhs_var_cat_{i}" for i in range(1, 9))
    if name == "lh_samp_frac_category":
        return tuple(f"lh_samp_frac_{i}" for i in range(1, 9))

    return (name,)


def _parse_registry(
    registry_path: str,
    sclr_dim: int = 0,
    edsclr_dim: int = 0,
    hydromet_list: Iterable[str] = (),
) -> Dict[str, Tuple[str, str, str]]:
    """Parse and expand a CLUBB stats registry while preserving source order."""
    definitions: list[StatsDefinition] = []
    path = Path(registry_path)
    text = path.read_text()
    entry_start = re.compile(r"^\s*entry\s*\(", re.IGNORECASE)
    entry_pattern = re.compile(
        r"^\s*entry\s*\(\s*\d+\s*\)\s*=\s*"
        r"(?P<quote>['\"])(?P<value>.*)(?P=quote)\s*,?\s*(?:!.*)?$",
        re.IGNORECASE,
    )
    for line_number, line in enumerate(text.splitlines(), start=1):
        if entry_start.match(line) is None:
            continue
        match = entry_pattern.match(line)
        if match is None:
            raise ValueError(
                f"Malformed stats registry entry at {path}:{line_number}: {line.strip()}"
            )

        fields = tuple(field.strip() for field in match.group("value").split("|"))
        if len(fields) == 4:
            name, grid, units, long_name = fields
        elif len(fields) == 3:
            name, grid, long_name = fields
            units = ""
        else:
            raise ValueError(
                f"Invalid stats registry entry at {path}:{line_number}; expected "
                "name | grid | units | long_name."
            )
        if not name or not grid:
            raise ValueError(
                f"Invalid stats registry entry at {path}:{line_number}; "
                "name and grid must be non-empty."
            )
        grid = {
            "nzt": GRID_ZT,
            "nzm": GRID_ZM,
            "lh_zm": GRID_ZM,
        }.get(grid, grid)
        if grid not in _VALID_GRIDS:
            raise ValueError(
                f"Invalid stats grid {grid!r} for variable {name!r} "
                f"at {path}:{line_number}."
            )
        for expanded_name in _expanded_names(
            name,
            sclr_dim=sclr_dim,
            edsclr_dim=edsclr_dim,
            hydromet_list=hydromet_list,
        ):
            definitions.append(StatsDefinition(expanded_name, grid, units, long_name))

    registry: Dict[str, Tuple[str, str, str]] = {}
    for definition in definitions:
        if definition.name not in registry:
            registry[definition.name] = (
                definition.grid,
                definition.units,
                definition.long_name,
            )
    return registry


def _format_time_units(day: int, month: int, year: int, time_initial: float) -> str:
    """Produce the UDUNITS reference-time string matching Fortran format_date."""
    del time_initial
    return f"seconds since {year:04d}-{month:02d}-{day:02d} 00:00:00.0"


class StatsWriter:
    """Manage CLUBB stats metadata, cadence, and host NetCDF output.

    JAX accumulation is owned exclusively by ``JaxStats``; this object handles
    the remaining host-side responsibilities of ``stats_netcdf.F90``.
    """

    def __init__(
        self,
        registry_path: str,
        output_path: str,
        nzt: int,
        nzm: int,
        ngrdcol: int,
        zt: np.ndarray,
        zm: np.ndarray,
        stats_tsamp: float,
        stats_tout: float,
        dt_main: float,
        day: int,
        month: int,
        year: int,
        time_initial: float,
        clubb_params_vals: Optional[np.ndarray] = None,
        param_names: Optional[list] = None,
        sclr_dim: int = 0,
        edsclr_dim: int = 0,
        hydromet_list: Iterable[str] = (),
        rad_zt: Optional[np.ndarray] = None,
        rad_zm: Optional[np.ndarray] = None,
        stats_tstart: Optional[float] = None,
        stats_tend: Optional[float] = None,
        ncol_total: Optional[int] = None,
        output_zt: Optional[np.ndarray] = None,
        output_zm: Optional[np.ndarray] = None,
        grid_remap_method: Optional[int] = None,
    ):
        self.nzt = nzt
        self.nzm = nzm
        self.ngrdcol = ngrdcol
        self.ncol_total = ngrdcol if ncol_total is None else int(ncol_total)
        if self.ncol_total < self.ngrdcol:
            raise ValueError("ncol_total must be greater than or equal to ngrdcol.")
        self.dt_main = dt_main
        self.time_initial = float(time_initial)
        self.stats_tstart = self.time_initial if stats_tstart is None else float(stats_tstart)
        self.stats_tend = None if stats_tend is None else float(stats_tend)

        if dt_main <= 0.0:
            raise ValueError("dt_main must be positive.")
        if stats_tsamp <= 0.0 or stats_tout <= 0.0:
            raise ValueError("stats_tsamp and stats_tout must be positive.")
        if self.stats_tstart < self.time_initial:
            raise ValueError("stats_tstart must be greater than or equal to time_initial.")
        if self.stats_tend is not None and self.stats_tend <= self.stats_tstart:
            raise ValueError("stats_tend must be greater than stats_tstart.")
        for label, interval in (
            ("stats_tsamp", stats_tsamp),
            ("stats_tout", stats_tout),
            ("stats_tstart", self.stats_tstart - self.time_initial),
        ):
            ratio = float(interval) / dt_main
            if not np.isclose(ratio, round(ratio), rtol=0.0, atol=1.0e-10):
                raise ValueError(f"{label} must align to a dt_main boundary.")
        if self.stats_tend is not None:
            ratio = (self.stats_tend - self.time_initial) / dt_main
            if not np.isclose(ratio, round(ratio), rtol=0.0, atol=1.0e-10):
                raise ValueError("stats_tend must align to a dt_main boundary.")

        # Sampling cadence (matching Fortran stats_init_api)
        self.stats_nsamp = max(1, round(stats_tsamp / dt_main))
        self.stats_nout = max(1, round(stats_tout / dt_main))
        if self.stats_nout % self.stats_nsamp != 0:
            raise ValueError("stats_tout must be an integer multiple of stats_tsamp.")
        self.samples_per_write = max(1, self.stats_nout // self.stats_nsamp)

        if (rad_zt is None) != (rad_zm is None):
            raise ValueError("rad_zt and rad_zm must be supplied together.")
        if (clubb_params_vals is None) != (param_names is None):
            raise ValueError("clubb_params_vals and param_names must be supplied together.")
        if clubb_params_vals is not None:
            params = np.asarray(clubb_params_vals)
            if params.ndim != 2 or params.shape[0] != self.ncol_total:
                raise ValueError("clubb_params_vals must have ncol_total rows.")
            if len(param_names) != params.shape[1]:
                raise ValueError("param_names must match the parameter-array width.")

        source_zt = np.asarray(zt).reshape(-1)[:nzt]
        source_zm = np.asarray(zm).reshape(-1)[:nzm]
        if source_zt.size != nzt or source_zm.size != nzm or nzt != nzm - 1:
            raise ValueError("Stats source grids must satisfy nzt = nzm - 1.")

        # Parse the variable registry
        self.registry = _parse_registry(
            registry_path,
            sclr_dim,
            edsclr_dim,
            hydromet_list,
        )
        if rad_zt is None:
            self.registry = {
                name: definition
                for name, definition in self.registry.items()
                if definition[0] not in {GRID_RAD_ZT, GRID_RAD_ZM}
            }
        self._grid_levels = {
            GRID_ZT: source_zt,
            GRID_ZM: source_zm,
            GRID_LH_ZT: source_zt,
            GRID_RAD_ZT: (
                np.asarray(rad_zt).reshape(-1) if rad_zt is not None else source_zt
            ),
            GRID_RAD_ZM: (
                np.asarray(rad_zm).reshape(-1) if rad_zm is not None else source_zm
            ),
        }
        remap_args = (output_zt, output_zm, grid_remap_method)
        if any(value is not None for value in remap_args) and not all(
            value is not None for value in remap_args
        ):
            raise ValueError("output_zt, output_zm, and grid_remap_method must be supplied together.")
        if output_zt is not None:
            output_zt = np.asarray(output_zt).reshape(-1)
            output_zm = np.asarray(output_zm).reshape(-1)
            if output_zt.size != output_zm.size - 1:
                raise ValueError("Stats output grids must satisfy nzt = nzm - 1.")
        self.grid_remap_method = (
            None if grid_remap_method is None or not output_path else int(grid_remap_method)
        )
        self._output_grid_levels = {
            GRID_ZT: (
                self._grid_levels[GRID_ZT]
                if output_zt is None
                else output_zt
            ),
            GRID_ZM: (
                self._grid_levels[GRID_ZM]
                if output_zm is None
                else output_zm
            ),
        }
        self._remap_state = None

        # Sampling flags
        self.l_sample = False
        self.l_last_sample = False
        self.l_reset = False
        self.enabled = True
        self._time_index = 0  # Number of NetCDF records already written.
        self.active_batch_offset = 0

        # NetCDF file setup
        self._ncid: Optional[object] = None
        self._varids: Dict[str, object] = {}
        self._time_varid: Optional[object] = None
        self._time_bnds_varid: Optional[object] = None
        self._lh_initialized = False
        self._lh_num_samples = 0
        self._lh_nzt = 0
        self._lh_nl_names: tuple[str, ...] = ()
        self._lh_u_names: tuple[str, ...] = ()
        self._lh_nl_varids: tuple[object, ...] = ()
        self._lh_u_varids: tuple[object, ...] = ()

        if output_path and not _HAS_NC:
            raise RuntimeError("netCDF4 is required when stats output_path is non-empty.")
        if _HAS_NC and output_path:
            self._open_netcdf(
                output_path, zt, zm, day, month, year, time_initial,
                clubb_params_vals, param_names,
            )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def begin_timestep(self, itime: int) -> Tuple[bool, bool]:
        """Begin a model step using the existing zero-based Python interface."""
        self.l_reset = False
        if not self.enabled:
            self.l_sample = False
            self.l_last_sample = False
            return False, False

        step_end_time = self.time_initial + (int(itime) + 1) * self.dt_main
        in_window = step_end_time > self.stats_tstart
        if self.stats_tend is not None:
            in_window = in_window and step_end_time <= self.stats_tend
        if not in_window:
            self.l_sample = False
            self.l_last_sample = False
            return False, False

        step_offset = round((step_end_time - self.stats_tstart) / self.dt_main)
        self.l_sample = step_offset % self.stats_nsamp == 0
        self.l_last_sample = step_offset % self.stats_nout == 0

        if step_offset == 1 or (step_offset - 1) % self.stats_nout == 0:
            self.l_reset = True

        return self.l_sample, self.l_last_sample

    def get_jax_layout(self) -> StatsLayout:
        """Return the public immutable metadata needed for JAX bank allocation."""
        names = tuple(self.registry)
        grids = tuple(self.registry[name][0] for name in names)
        grid_nlev = tuple(
            1 if grid in {GRID_SFC, GRID_LH_SFC} else int(self._grid_levels[grid].size)
            for grid in GRID_ORDER
        )
        return StatsLayout(
            names=names,
            grids=grids,
            grid_nlev=grid_nlev,
            samples_per_write=self.samples_per_write,
        )

    def end_timestep(self, t_out: float, jax_stats):
        """Average one output window and write a NetCDF record.

        Accumulation is owned exclusively by the immutable JAX stats banks.
        """
        if not self.l_last_sample:
            return jax_stats
        if not self.enabled:
            return jax_stats

        buffers, nsamples, issues, l_in_budget = jax_stats.to_host_banks()
        jax_stats.report_issues(issues, l_in_budget)
        values = tuple(
            (
                name,
                buffers[jax_stats.name_to_slot[name][0]][jax_stats.name_to_slot[name][1]],
                nsamples[jax_stats.name_to_slot[name][0]][jax_stats.name_to_slot[name][1]],
            )
            for name in self.registry
        )

        for name, _, cnt in values:
            if np.any(cnt > jax_stats.samples_per_write):
                warnings.warn(f"stats oversampling warning for {name}", RuntimeWarning)

        self._time_index += 1
        if self._ncid is not None and self._time_varid is not None:
            self._time_varid[self._time_index - 1] = t_out
            step_offset = round((float(t_out) - self.stats_tstart) / self.dt_main)
            window_start_step = ((step_offset - 1) // self.stats_nout) * self.stats_nout
            window_start_time = self.stats_tstart + window_start_step * self.dt_main
            self._time_bnds_varid[self._time_index - 1, :] = (
                window_start_time,
                t_out,
            )

        if self._ncid is not None:
            for name, buf, cnt in values:
                if np.all(cnt <= 0):
                    continue
                grid = self.registry[name][0]
                avg = np.where(cnt > 0, buf / np.maximum(cnt, 1), 0.0)
                avg = self._remap_values(grid, avg)
                varid = self._varids.get(name)
                if varid is None:
                    continue
                t = self._time_index - 1
                col_start = self.active_batch_offset
                col_end = col_start + self.ngrdcol
                if grid in _PROFILE_GRIDS:
                    varid[t, :, col_start:col_end] = avg.T
                else:
                    varid[t, col_start:col_end] = avg[:, 0]

        # Periodically flush dirty HDF5 chunks to disk so the cache does not grow
        # unbounded over a long run. Without this, a per-step-stats run of a case with
        # a large variable set (e.g. Morrison: ~500 vars × many levels) accumulates
        # hundreds of records' worth of un-flushed chunks and OOMs the process mid-write
        # (→ a truncated/corrupt NetCDF). Every-20-records bounds the dirty cache while
        # keeping the sync I/O negligible next to the physics. (Iter322 — mpace_a OOM.)
        if self._ncid is not None and self._time_index % 20 == 0:
            self._ncid.sync()

        # Reset flags
        self.l_sample = False
        self.l_last_sample = False
        return jax_stats.clear_output_diagnostics()

    def reset(self) -> None:
        """Mirror ``stats_reset_api`` for another run or output batch."""
        if self.active_batch_offset + self.ngrdcol >= self.ncol_total:
            self.active_batch_offset = 0
        self.l_sample = False
        self.l_last_sample = False
        self.l_reset = False
        self._time_index = 0

    def start_next_batch(self) -> None:
        """Advance the active NetCDF column slice by one runtime batch."""
        new_offset = self.active_batch_offset + self.ngrdcol
        if new_offset + self.ngrdcol > self.ncol_total:
            raise ValueError("New stats batch exceeds ncol_total.")
        self.active_batch_offset = new_offset

    def get_source_grid(self) -> tuple[np.ndarray, np.ndarray]:
        """Return copies of the runtime zt and zm grids."""
        if self._remap_state is not None:
            zt, zm = self._remap_state[:2]
            return zt.copy(), zm.copy()
        zt = np.broadcast_to(
            self._grid_levels[GRID_ZT],
            (self.ngrdcol, self._grid_levels[GRID_ZT].size),
        )
        zm = np.broadcast_to(
            self._grid_levels[GRID_ZM],
            (self.ngrdcol, self._grid_levels[GRID_ZM].size),
        )
        return zt.copy(), zm.copy()

    def update_grid(
        self,
        zt_src: np.ndarray,
        zm_src: np.ndarray,
        rho_vals: np.ndarray,
        rho_levels: np.ndarray,
        p_sfc: np.ndarray,
    ) -> None:
        """Update source-grid and density data used by output remapping."""
        if self.grid_remap_method is None:
            return
        zt_src = np.asarray(zt_src, dtype=np.float64)
        zm_src = np.asarray(zm_src, dtype=np.float64)
        if zt_src.ndim == 1:
            zt_src = np.broadcast_to(zt_src, (self.ngrdcol, zt_src.size)).copy()
        if zm_src.ndim == 1:
            zm_src = np.broadcast_to(zm_src, (self.ngrdcol, zm_src.size)).copy()
        rho_vals = np.asarray(rho_vals, dtype=np.float64)
        rho_levels = np.asarray(rho_levels, dtype=np.float64)
        p_sfc = np.asarray(p_sfc, dtype=np.float64)
        if (
            zt_src.ndim != 2
            or zm_src.ndim != 2
            or rho_vals.ndim != 2
            or rho_levels.ndim != 2
            or zt_src.shape[0] != self.ngrdcol
            or zm_src.shape[0] != self.ngrdcol
            or rho_vals.shape[0] != self.ngrdcol
            or rho_levels.shape != rho_vals.shape
            or p_sfc.shape != (self.ngrdcol,)
            or zt_src.shape[1] != zm_src.shape[1] - 1
        ):
            raise ValueError("Invalid per-column grid or density shape for stats remapping.")
        self._grid_levels[GRID_ZT] = zt_src[0].copy()
        self._grid_levels[GRID_ZM] = zm_src[0].copy()
        self._remap_state = (
            zt_src,
            zm_src,
            rho_vals,
            rho_levels,
            p_sfc,
        )

    def _remap_values(self, grid: str, values: np.ndarray) -> np.ndarray:
        if self.grid_remap_method is None or grid not in {GRID_ZT, GRID_ZM}:
            return values
        if self._remap_state is None:
            raise RuntimeError("stats_update_grid must be called before remapped output is written.")

        from clubb_jax.src.CLUBB_core.remapping_module import remap_vals_to_target

        zt_src, zm_src, rho_vals, rho_levels, p_sfc = self._remap_state
        target_zt = np.broadcast_to(
            self._output_grid_levels[GRID_ZT],
            (self.ngrdcol, self._output_grid_levels[GRID_ZT].size),
        )
        target_zm = np.broadcast_to(
            self._output_grid_levels[GRID_ZM],
            (self.ngrdcol, self._output_grid_levels[GRID_ZM].size),
        )
        if grid == GRID_ZT:
            source_edges = zm_src
            target_edges = target_zm
        else:
            source_edges = np.concatenate((zm_src[:, :1], zt_src, zm_src[:, -1:]), axis=1)
            target_edges = np.concatenate((target_zm[:, :1], target_zt, target_zm[:, -1:]), axis=1)
        return np.asarray(
            remap_vals_to_target(
                source_edges,
                target_edges,
                rho_vals,
                rho_levels,
                values,
                p_sfc,
                self.grid_remap_method,
                iv=1,
            )
        )

    def get_stats_config(self) -> tuple[int, ...]:
        """Return the scalar configuration tuple exposed by the former API."""
        return (
            int(self.enabled),
            self.ngrdcol,
            len(self.registry),
            self.stats_nsamp,
            self.stats_nout,
            self.samples_per_write,
            self._time_index,
            int(self.l_sample),
            int(self.l_last_sample),
        )

    def get_stats_var_meta(self, ivar: int) -> tuple[str, str, str, str, int, int]:
        """Return metadata for one variable using zero-based Python indexing."""
        try:
            name = tuple(self.registry)[int(ivar)]
        except IndexError as exc:
            raise IndexError(f"Stats variable index {ivar} is out of range.") from exc
        grid, units, long_name = self.registry[name]
        nz = 1 if grid in {GRID_SFC, GRID_LH_SFC} else self._grid_levels[grid].size
        return name, grid, units, long_name, _GRID_IDS[grid], nz

    def stats_lh_samples_init(
        self,
        num_samples: int,
        nzt: int,
        nl_var_names: Iterable[str],
        u_var_names: Iterable[str],
        zt_vals: np.ndarray,
    ) -> None:
        """Define optional SILHS sample fields in the active stats file."""
        if not self.enabled or self._ncid is None or num_samples <= 0 or nzt <= 0:
            return
        num_samples = int(num_samples)
        nzt = int(nzt)
        nl_names = tuple(f"lh_nl_{str(name).strip()}" for name in nl_var_names)
        u_names = tuple(f"lh_u_{str(name).strip()}" for name in u_var_names)
        zt_vals = np.asarray(zt_vals, dtype=np.float64)
        if zt_vals.shape != (nzt,):
            raise ValueError(f"SILHS zt values must have shape ({nzt},).")
        if self._lh_initialized:
            if self._lh_num_samples != num_samples or self._lh_nzt != nzt:
                raise ValueError("SILHS stats initialization does not match the existing definition.")
            return

        ds = self._ncid
        if "lh_sample_number" in ds.dimensions:
            if len(ds.dimensions["lh_sample_number"]) != num_samples:
                raise ValueError("Existing lh_sample_number dimension has the wrong size.")
        else:
            ds.createDimension("lh_sample_number", num_samples)
        if "lh_zt" in ds.dimensions:
            if len(ds.dimensions["lh_zt"]) != nzt:
                raise ValueError("Existing lh_zt dimension has the wrong size.")
            lh_zt = ds.variables["lh_zt"]
        else:
            ds.createDimension("lh_zt", nzt)
            lh_zt = ds.createVariable("lh_zt", "f8", ("lh_zt",))
            lh_zt.units = "m"
        lh_zt[:] = zt_vals

        dims = ("time", "lh_zt", "lh_sample_number", "col")
        self._lh_nl_varids = tuple(
            ds.createVariable(name, "f8", dims, fill_value=0.0) for name in nl_names
        )
        self._lh_u_varids = tuple(
            ds.createVariable(name, "f8", dims, fill_value=0.0) for name in u_names
        )
        self._lh_num_samples = num_samples
        self._lh_nzt = nzt
        self._lh_nl_names = nl_names
        self._lh_u_names = u_names
        self._lh_initialized = True

    def stats_lh_samples_write_lognormal(self, samples: np.ndarray) -> None:
        """Write SILHS lognormal samples shaped ``(col, sample, zt, var)``."""
        if not self.enabled or self._ncid is None or not self._lh_initialized or not self.l_last_sample:
            return
        samples = np.asarray(samples, dtype=np.float64)
        expected = (
            self.ngrdcol,
            self._lh_num_samples,
            self._lh_nzt,
            len(self._lh_nl_varids),
        )
        if samples.shape != expected:
            raise ValueError(f"SILHS lognormal samples must have shape {expected}.")
        t = self._time_index
        col = slice(self.active_batch_offset, self.active_batch_offset + self.ngrdcol)
        for ivar, varid in enumerate(self._lh_nl_varids):
            varid[t, :, :, col] = samples[:, :, :, ivar].transpose(2, 1, 0)

    def stats_lh_samples_write_uniform(
        self,
        uniform_vals: np.ndarray,
        mixture_comp: np.ndarray,
        sample_weights: np.ndarray,
    ) -> None:
        """Write SILHS uniform, mixture-component, and sample-weight fields."""
        if not self.enabled or self._ncid is None or not self._lh_initialized or not self.l_last_sample:
            return
        uniform_vals = np.asarray(uniform_vals, dtype=np.float64)
        mixture_comp = np.asarray(mixture_comp)
        sample_weights = np.asarray(sample_weights, dtype=np.float64)
        base_shape = (self.ngrdcol, self._lh_num_samples, self._lh_nzt)
        if (
            uniform_vals.ndim != 4
            or uniform_vals.shape[:3] != base_shape
            or mixture_comp.shape != base_shape
            or sample_weights.shape != base_shape
        ):
            raise ValueError("SILHS uniform sample arrays have inconsistent shapes.")
        if uniform_vals.shape[3] + 2 != len(self._lh_u_varids):
            raise ValueError("SILHS uniform variable count does not match its definition.")
        t = self._time_index
        col = slice(self.active_batch_offset, self.active_batch_offset + self.ngrdcol)
        for ivar in range(uniform_vals.shape[3]):
            self._lh_u_varids[ivar][t, :, :, col] = uniform_vals[:, :, :, ivar].transpose(2, 1, 0)
        self._lh_u_varids[-2][t, :, :, col] = mixture_comp.transpose(2, 1, 0)
        self._lh_u_varids[-1][t, :, :, col] = sample_weights.transpose(2, 1, 0)

    def finalize(self) -> None:
        """Flush and close the NetCDF file."""
        if self._ncid is not None:
            try:
                self._ncid.close()
            except Exception:
                pass
            self._ncid = None
        self.enabled = False

    # ------------------------------------------------------------------
    # NetCDF helpers
    # ------------------------------------------------------------------

    def _open_netcdf(
        self,
        output_path: str,
        zt: np.ndarray,
        zm: np.ndarray,
        day: int,
        month: int,
        year: int,
        time_initial: float,
        clubb_params_vals: Optional[np.ndarray],
        param_names: Optional[list],
    ) -> None:
        """Create and define the output NetCDF file (mirrors stats_open_netcdf)."""
        ds = nc.Dataset(output_path, "w", format="NETCDF4_CLASSIC")
        self._ncid = ds

        # Dimensions
        ds.createDimension("time", None)  # unlimited
        ds.createDimension("bnds", 2)
        ds.createDimension("col", self.ncol_total)
        ds.createDimension("sfc", 1)
        ds.createDimension("zt", self._output_grid_levels[GRID_ZT].size)
        ds.createDimension("zm", self._output_grid_levels[GRID_ZM].size)

        if any(grid == GRID_LH_ZT for grid, _, _ in self.registry.values()):
            ds.createDimension("lh_zt", self._grid_levels[GRID_LH_ZT].size)
        has_radiation_stats = any(
            grid in {GRID_RAD_ZT, GRID_RAD_ZM}
            for grid, _, _ in self.registry.values()
        )
        if has_radiation_stats:
            ds.createDimension("rad_zt", self._grid_levels[GRID_RAD_ZT].size)
            ds.createDimension("rad_zm", self._grid_levels[GRID_RAD_ZM].size)

        # Coordinate variables
        tv = ds.createVariable("time", "f8", ("time",))
        tv.units = _format_time_units(day, month, year, time_initial)
        tv.bounds = "time_bnds"
        self._time_varid = tv

        tbv = ds.createVariable("time_bnds", "f8", ("time", "bnds"))
        tbv.units = tv.units
        tbv.interval_semantics = "(start, end]"
        self._time_bnds_varid = tbv

        ds.model_time_initial = self.time_initial
        ds.stats_tstart = self.stats_tstart
        ds.stats_time_window_mode = "open_ended" if self.stats_tend is None else "explicit"
        if self.stats_tend is not None:
            ds.stats_tend = self.stats_tend
        ds.stats_tsamp = self.stats_nsamp * self.dt_main
        ds.stats_tout = self.stats_nout * self.dt_main

        cv = ds.createVariable("col", "f8", ("col",))
        cv.units = "index"
        cv[:] = np.arange(1, self.ncol_total + 1, dtype=np.float64)

        zmv = ds.createVariable("zm", "f8", ("zm",))
        zmv.units = "m"
        zmv[:] = self._output_grid_levels[GRID_ZM]

        ztv = ds.createVariable("zt", "f8", ("zt",))
        ztv.units = "m"
        ztv[:] = self._output_grid_levels[GRID_ZT]

        if "lh_zt" in ds.dimensions:
            lhv = ds.createVariable("lh_zt", "f8", ("lh_zt",))
            lhv.units = "m"
            lhv[:] = self._grid_levels[GRID_LH_ZT]

        if has_radiation_stats:
            rad_zt_var = ds.createVariable("rad_zt", "f8", ("rad_zt",))
            rad_zm_var = ds.createVariable("rad_zm", "f8", ("rad_zm",))
            rad_zt_var.units = "m"
            rad_zm_var.units = "m"
            rad_zt_var[:] = self._grid_levels[GRID_RAD_ZT]
            rad_zm_var[:] = self._grid_levels[GRID_RAD_ZM]

        # clubb_params metadata
        if clubb_params_vals is not None and param_names is not None:
            nparam = len(param_names)
            ds.createDimension("param", nparam)
            pstrlen = 32
            ds.createDimension("param_strlen", pstrlen)

            pv = ds.createVariable("param", "f8", ("param",))
            pv.units = "index"
            pv[:] = np.arange(1, nparam + 1, dtype=np.float64)

            pnv = ds.createVariable("param_name", "S1", ("param", "param_strlen"))
            for i, pn in enumerate(param_names):
                padded = pn.ljust(pstrlen)[:pstrlen]
                pnv[i, :] = np.array(list(padded), dtype="S1")

            cpv = ds.createVariable("clubb_params", "f8", ("param", "col"))
            cpv.long_name = "clubb_params"
            vals = np.asarray(clubb_params_vals, dtype=np.float64)
            # vals shape: (ngrdcol, nparam) → transpose to (nparam, ngrdcol)
            cpv[:, :vals.shape[0]] = vals.T

        # All registered stats variables
        grid_to_dim = {
            GRID_ZT: ("time", "zt", "col"),
            GRID_ZM: ("time", "zm", "col"),
            GRID_LH_ZT: ("time", "lh_zt", "col"),
            GRID_SFC: ("time", "col"),
            GRID_LH_SFC: ("time", "col"),
            GRID_RAD_ZT: ("time", "rad_zt", "col"),
            GRID_RAD_ZM: ("time", "rad_zm", "col"),
        }
        for name, (grid, units, long_name) in self.registry.items():
            dims = grid_to_dim.get(grid)
            if dims is None:
                continue
            v = ds.createVariable(name, "f8", dims, fill_value=0.0)
            v.long_name = long_name
            v.units = units
            self._varids[name] = v
