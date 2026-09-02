"""JAX-native CLUBB statistics accumulation.

The public methods intentionally mirror ``stats_netcdf.F90`` call sites. Names
and grid slots are static trace-time metadata; dynamic pytree leaves contain
only compact grid-grouped accumulation banks and budget diagnostics. NetCDF I/O
is handled by the host ``StatsWriter`` at output boundaries.
"""

from __future__ import annotations

from functools import lru_cache
from types import MappingProxyType
import warnings

import jax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
from clubb_jax.src.CLUBB_core.stats_metadata import GRID_ORDER, StatsLayout

configure_jax_precision()
import jax.numpy as jnp
import numpy as np


ISSUE_UPDATE_AFTER_BUDGET = 1
ISSUE_BEGIN_TWICE = 2
ISSUE_UPDATE_WITHOUT_BEGIN = 4
ISSUE_FINALIZE_WITHOUT_BEGIN = 8

_GRID_TO_ID = {grid: index for index, grid in enumerate(GRID_ORDER)}


def _replace(items: tuple, index: int, value) -> tuple:
    return items[:index] + (value,) + items[index + 1:]


@lru_cache(maxsize=None)
def _slot_metadata(names: tuple[str, ...], grids: tuple[str, ...]):
    """Build immutable trace-time lookup maps once per registry layout."""
    next_slot = [0] * len(GRID_ORDER)
    name_to_id: dict[str, int] = {}
    name_to_slot: dict[str, tuple[int, int]] = {}
    for var_id, (name, grid) in enumerate(zip(names, grids)):
        bank_id = _GRID_TO_ID[grid]
        name_to_id[name] = var_id
        name_to_slot[name] = (bank_id, next_slot[bank_id])
        next_slot[bank_id] += 1
    return MappingProxyType(name_to_id), MappingProxyType(name_to_slot)


@jax.tree_util.register_pytree_node_class
class JaxStats:
    """Functional, JIT-compatible implementation of CLUBB stats updates."""

    def __init__(
        self,
        *,
        l_sample: bool,
        names: tuple[str, ...],
        grids: tuple[str, ...],
        ncol: int,
        grid_nlev: tuple[int, ...],
        samples_per_write: int,
        buffers: tuple,
        nsamples: tuple,
        l_budget,
        l_in_budget,
        issues,
    ):
        self.l_sample = bool(l_sample)
        self.names = tuple(names)
        self.grids = tuple(grids)
        self.ncol = int(ncol)
        self.grid_nlev = tuple(int(value) for value in grid_nlev)
        self.samples_per_write = int(samples_per_write)
        self.buffers = tuple(buffers)
        self.nsamples = tuple(nsamples)
        self.l_budget = l_budget
        self.l_in_budget = l_in_budget
        self.issues = issues

        self.name_to_id, self.name_to_slot = _slot_metadata(self.names, self.grids)

    @classmethod
    def empty(
        cls,
        *,
        l_sample: bool,
        names: tuple[str, ...],
        ncol: int,
        max_nlev: int,
        grids: tuple[str, ...] | None = None,
        grid_nlev: tuple[int, ...] | None = None,
        samples_per_write: int = 1,
    ):
        """Create empty grid-grouped accumulation banks."""
        names = tuple(names)
        if grids is None:
            grids = tuple("zt" for _ in names)
        if len(grids) != len(names):
            raise ValueError("grids must contain one entry per stats variable.")
        if any(grid not in GRID_ORDER for grid in grids):
            raise ValueError("Unsupported grid in JaxStats metadata.")

        if grid_nlev is None:
            grid_nlev = (
                int(max_nlev),
                int(max_nlev),
                1,
                int(max_nlev),
                1,
                int(max_nlev),
                int(max_nlev),
            )
        grid_counts = tuple(sum(grid == key for grid in grids) for key in GRID_ORDER)
        buffers = tuple(
            jnp.zeros((count, int(ncol), int(nlev)), dtype=jnp.float64)
            for count, nlev in zip(grid_counts, grid_nlev)
        )
        nsamples = tuple(
            jnp.zeros((count, int(ncol), int(nlev)), dtype=jnp.int32)
            for count, nlev in zip(grid_counts, grid_nlev)
        )
        return cls(
            l_sample=l_sample,
            names=names,
            grids=tuple(grids),
            ncol=ncol,
            grid_nlev=tuple(grid_nlev),
            samples_per_write=samples_per_write,
            buffers=buffers,
            nsamples=nsamples,
            l_budget=jnp.zeros((len(names),), dtype=bool),
            l_in_budget=jnp.zeros((len(names),), dtype=bool),
            issues=jnp.zeros((len(names),), dtype=jnp.int32),
        )

    @classmethod
    def from_layout(cls, layout: StatsLayout, *, ncol: int):
        """Build fixed-shape banks from public immutable host metadata."""
        return cls.empty(
            l_sample=False,
            names=layout.names,
            grids=layout.grids,
            ncol=int(ncol),
            max_nlev=max(layout.grid_nlev, default=1),
            grid_nlev=layout.grid_nlev,
            samples_per_write=layout.samples_per_write,
        )

    def tree_flatten(self):
        children = (
            self.buffers,
            self.nsamples,
            self.l_budget,
            self.l_in_budget,
            self.issues,
        )
        aux_data = (
            self.l_sample,
            self.names,
            self.grids,
            self.ncol,
            self.grid_nlev,
            self.samples_per_write,
        )
        return children, aux_data

    @classmethod
    def tree_unflatten(cls, aux_data, children):
        l_sample, names, grids, ncol, grid_nlev, samples_per_write = aux_data
        buffers, nsamples, l_budget, l_in_budget, issues = children
        return cls(
            l_sample=l_sample,
            names=names,
            grids=grids,
            ncol=ncol,
            grid_nlev=grid_nlev,
            samples_per_write=samples_per_write,
            buffers=buffers,
            nsamples=nsamples,
            l_budget=l_budget,
            l_in_budget=l_in_budget,
            issues=issues,
        )

    def _new(
        self,
        *,
        l_sample=None,
        buffers=None,
        nsamples=None,
        l_budget=None,
        l_in_budget=None,
        issues=None,
    ):
        return type(self)(
            l_sample=self.l_sample if l_sample is None else bool(l_sample),
            names=self.names,
            grids=self.grids,
            ncol=self.ncol,
            grid_nlev=self.grid_nlev,
            samples_per_write=self.samples_per_write,
            buffers=self.buffers if buffers is None else buffers,
            nsamples=self.nsamples if nsamples is None else nsamples,
            l_budget=self.l_budget if l_budget is None else l_budget,
            l_in_budget=self.l_in_budget if l_in_budget is None else l_in_budget,
            issues=self.issues if issues is None else issues,
        )

    def begin_timestep(self, *, l_sample: bool, reset_accumulators: bool):
        """Set static sampling state and optionally clear the output-window banks."""
        if not reset_accumulators:
            return self._new(l_sample=l_sample)
        buffers = tuple(jnp.zeros_like(value) for value in self.buffers)
        nsamples = tuple(jnp.zeros_like(value) for value in self.nsamples)
        return self._new(
            l_sample=l_sample,
            buffers=buffers,
            nsamples=nsamples,
            l_in_budget=jnp.zeros_like(self.l_in_budget),
            issues=jnp.zeros_like(self.issues),
        )

    def var_on_stats_list(self, name: str) -> bool:
        return bool(name) and name.strip() in self.name_to_id

    def _record_issue(self, issues, var_id: int, condition, bit: int):
        old_value = issues[var_id]
        new_value = jnp.bitwise_or(old_value, jnp.asarray(bit, dtype=jnp.int32))
        return issues.at[var_id].set(jnp.where(condition, new_value, old_value))

    def _apply_payload(
        self,
        *,
        name: str,
        values,
        sign: int,
        count_sample: bool,
        allowed,
        icol: int | None,
        level: int | None,
        buffers: tuple,
        nsamples: tuple,
        issues,
    ):
        var_id = self.name_to_id[name]
        bank_id, slot = self.name_to_slot[name]
        bank = buffers[bank_id]
        count_bank = nsamples[bank_id]
        nlev = self.grid_nlev[bank_id]
        arr = jax.lax.stop_gradient(jnp.asarray(values, dtype=jnp.float64))
        if arr.ndim > 2:
            raise ValueError(f"Stats values must have rank <= 2; got rank {arr.ndim}.")

        if arr.ndim == 0:
            if icol is None:
                raise ValueError("Scalar stats operations require icol.")
            row = int(icol)
            column = 0 if level is None else int(level)
            if row < 0 or row >= self.ncol or column < 0 or column >= nlev:
                raise ValueError(f"Stats selector is out of bounds for {name!r}.")
            candidate = bank.at[slot, row, column].add(sign * arr)
            count_candidate = count_bank.at[slot, row, column].add(int(count_sample))
        elif arr.ndim == 1:
            if level is not None:
                raise ValueError("Rank-1 stats operations do not accept level.")
            if icol is None:
                if nlev != 1 or arr.shape != (self.ncol,):
                    raise ValueError(
                        f"Rank-1 all-column stats update for {name!r} requires "
                        f"surface shape ({self.ncol},); got {arr.shape}."
                    )
                candidate = bank.at[slot, :, 0].add(sign * arr)
                count_candidate = count_bank.at[slot, :, 0].add(int(count_sample))
            else:
                row = int(icol)
                if row < 0 or row >= self.ncol:
                    raise ValueError(f"Stats column selector is out of bounds for {name!r}.")
                if arr.shape != (nlev,):
                    raise ValueError(
                        f"Rank-1 column stats update for {name!r} requires shape "
                        f"({nlev},); got {arr.shape}."
                    )
                candidate = bank.at[slot, row, :].add(sign * arr)
                count_candidate = count_bank.at[slot, row, :].add(int(count_sample))
        else:
            if icol is not None or level is not None:
                raise ValueError("Rank-2 stats operations do not accept icol or level.")
            if arr.shape != (self.ncol, nlev):
                raise ValueError(
                    f"Rank-2 stats update for {name!r} requires "
                    f"({self.ncol}, {nlev}); got {arr.shape}."
                )
            candidate = bank.at[slot, :, :].add(sign * arr)
            count_candidate = count_bank.at[slot, :, :].add(int(count_sample))

        bank = jnp.where(allowed, candidate, bank)
        count_bank = jnp.where(allowed, count_candidate, count_bank)
        return (
            _replace(buffers, bank_id, bank),
            _replace(nsamples, bank_id, count_bank),
            issues,
        )

    def update(self, name: str, values, *, icol: int | None = None, level: int | None = None):
        """Mirror ``stats_update`` with direct JAX-bank accumulation."""
        clean_name = name.strip()
        if not self.l_sample or not clean_name or clean_name not in self.name_to_id:
            return self
        var_id = self.name_to_id[clean_name]
        allowed = jnp.logical_not(self.l_budget[var_id])
        issues = self._record_issue(
            self.issues, var_id, jnp.logical_not(allowed), ISSUE_UPDATE_AFTER_BUDGET
        )
        buffers, nsamples, issues = self._apply_payload(
            name=clean_name,
            values=values,
            sign=1,
            count_sample=True,
            allowed=allowed,
            icol=icol,
            level=level,
            buffers=self.buffers,
            nsamples=self.nsamples,
            issues=issues,
        )
        return self._new(buffers=buffers, nsamples=nsamples, issues=issues)

    def begin_budget(self, name: str, values, *, icol: int | None = None):
        """Mirror ``stats_begin_budget``."""
        clean_name = name.strip()
        if not self.l_sample or not clean_name or clean_name not in self.name_to_id:
            return self
        var_id = self.name_to_id[clean_name]
        allowed = jnp.logical_not(self.l_in_budget[var_id])
        issues = self._record_issue(
            self.issues, var_id, jnp.logical_not(allowed), ISSUE_BEGIN_TWICE
        )
        buffers, nsamples, issues = self._apply_payload(
            name=clean_name,
            values=values,
            sign=-1,
            count_sample=False,
            allowed=allowed,
            icol=icol,
            level=None,
            buffers=self.buffers,
            nsamples=self.nsamples,
            issues=issues,
        )
        l_budget = self.l_budget.at[var_id].set(jnp.where(allowed, True, self.l_budget[var_id]))
        l_in_budget = self.l_in_budget.at[var_id].set(
            jnp.where(allowed, True, self.l_in_budget[var_id])
        )
        return self._new(
            buffers=buffers,
            nsamples=nsamples,
            l_budget=l_budget,
            l_in_budget=l_in_budget,
            issues=issues,
        )

    def update_budget(
        self,
        name: str,
        values,
        *,
        icol: int | None = None,
        level: int | None = None,
    ):
        """Mirror ``stats_update_budget``."""
        clean_name = name.strip()
        if not self.l_sample or not clean_name or clean_name not in self.name_to_id:
            return self
        var_id = self.name_to_id[clean_name]
        allowed = self.l_in_budget[var_id]
        issues = self._record_issue(
            self.issues, var_id, jnp.logical_not(allowed), ISSUE_UPDATE_WITHOUT_BEGIN
        )
        buffers, nsamples, issues = self._apply_payload(
            name=clean_name,
            values=values,
            sign=1,
            count_sample=False,
            allowed=allowed,
            icol=icol,
            level=level,
            buffers=self.buffers,
            nsamples=self.nsamples,
            issues=issues,
        )
        return self._new(buffers=buffers, nsamples=nsamples, issues=issues)

    def finalize_budget(
        self,
        name: str,
        values,
        *,
        icol: int | None = None,
        l_count_sample: bool = True,
    ):
        """Mirror ``stats_finalize_budget``."""
        clean_name = name.strip()
        if not self.l_sample or not clean_name or clean_name not in self.name_to_id:
            return self
        var_id = self.name_to_id[clean_name]
        allowed = self.l_in_budget[var_id]
        issues = self._record_issue(
            self.issues, var_id, jnp.logical_not(allowed), ISSUE_FINALIZE_WITHOUT_BEGIN
        )
        buffers, nsamples, issues = self._apply_payload(
            name=clean_name,
            values=values,
            sign=1,
            count_sample=bool(l_count_sample),
            allowed=allowed,
            icol=icol,
            level=None,
            buffers=self.buffers,
            nsamples=self.nsamples,
            issues=issues,
        )
        l_in_budget = self.l_in_budget.at[var_id].set(
            jnp.where(allowed, False, self.l_in_budget[var_id])
        )
        return self._new(
            buffers=buffers,
            nsamples=nsamples,
            l_in_budget=l_in_budget,
            issues=issues,
        )

    def to_host_banks(self):
        """Transfer all output banks and compact diagnostics in one device operation."""
        buffers, nsamples, issues, l_in_budget = jax.device_get(
            (self.buffers, self.nsamples, self.issues, self.l_in_budget)
        )
        return (
            tuple(np.asarray(value) for value in buffers),
            tuple(np.asarray(value) for value in nsamples),
            np.asarray(issues),
            np.asarray(l_in_budget),
        )

    def report_issues(self, issues: np.ndarray, l_in_budget: np.ndarray) -> None:
        """Emit Fortran-equivalent diagnostics after leaving compiled code."""
        messages = (
            (ISSUE_UPDATE_AFTER_BUDGET, "stats_update called for budget variable"),
            (ISSUE_BEGIN_TWICE, "stats budget begin twice for"),
            (ISSUE_UPDATE_WITHOUT_BEGIN, "stats budget update without begin for"),
            (ISSUE_FINALIZE_WITHOUT_BEGIN, "stats budget finalize without begin for"),
        )
        for var_id, name in enumerate(self.names):
            for bit, message in messages:
                if int(issues[var_id]) & bit:
                    warnings.warn(f"{message} {name}", RuntimeWarning)
            if bool(l_in_budget[var_id]):
                warnings.warn(f"stats budget left mid-update at end_timestep for {name}", RuntimeWarning)

    def clear_output_diagnostics(self):
        return self._new(
            l_in_budget=jnp.zeros_like(self.l_in_budget),
            issues=jnp.zeros_like(self.issues),
        )

    def reset(self):
        """Mirror a full stats reset while retaining static registry metadata."""
        return self._new(
            l_sample=False,
            buffers=tuple(jnp.zeros_like(value) for value in self.buffers),
            nsamples=tuple(jnp.zeros_like(value) for value in self.nsamples),
            l_budget=jnp.zeros_like(self.l_budget),
            l_in_budget=jnp.zeros_like(self.l_in_budget),
            issues=jnp.zeros_like(self.issues),
        )
