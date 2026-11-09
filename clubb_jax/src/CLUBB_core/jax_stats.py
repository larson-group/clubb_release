"""JAX-native CLUBB statistics accumulation.

The public methods intentionally mirror ``stats_netcdf.F90`` call sites. Names
and grid slots are static trace-time metadata; dynamic pytree leaves contain
the sampling flag, compact grid-grouped banks and budget diagnostics. NetCDF I/O
is handled by the host ``StatsWriter`` at output boundaries.

Observed Fortran differences retained in this diagnostic pass:
- Source budget messages are reported from accumulated issue bits at the host output boundary; exact occurrence counts and source-call timing are not retained.
"""

from __future__ import annotations

from functools import lru_cache
from types import MappingProxyType
import warnings

import jax

from clubb_jax.src.CLUBB_core.clubb_precision import core_rknd
from clubb_jax.src.CLUBB_core.stats_metadata import GRID_ORDER, StatsLayout

import jax.numpy as jnp
import numpy as np


ISSUE_UPDATE_AFTER_BUDGET = 1
ISSUE_BEGIN_TWICE = 2
ISSUE_UPDATE_WITHOUT_BEGIN = 4
ISSUE_FINALIZE_WITHOUT_BEGIN = 8

_GRID_TO_ID = {grid: index for index, grid in enumerate(GRID_ORDER)}


def _replace(items: tuple, index: int, value) -> tuple:
    replaced = items[:index] + (value,) + items[index + 1:]
    return replaced


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
        self.l_sample = l_sample
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
            jnp.zeros((count, int(ncol), int(nlev)))
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
            self.l_sample,
            self.buffers,
            self.nsamples,
            self.l_budget,
            self.l_in_budget,
            self.issues,
        )
        aux_data = (
            self.names,
            self.grids,
            self.ncol,
            self.grid_nlev,
            self.samples_per_write,
        )
        return children, aux_data

    @classmethod
    def tree_unflatten(cls, aux_data, children):
        names, grids, ncol, grid_nlev, samples_per_write = aux_data
        l_sample, buffers, nsamples, l_budget, l_in_budget, issues = children
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
            l_sample=self.l_sample if l_sample is None else l_sample,
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

    def begin_timestep(self, *, l_sample: bool, l_reset_accumulators: bool):
        """Set runtime sampling state and optionally clear the output-window banks."""
        if not l_reset_accumulators:
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
        var_on_stats_list = bool(name) and name.strip() in self.name_to_id
        return var_on_stats_list

    def sampled(self, accumulate):
        """Apply a stats-only callback on sampled timesteps, returning new stats.

        The registry remains static. Within the callback sampling is known to
        be enabled, so nested stats updates need no additional runtime guards.
        Host calls with a Python flag keep their ordinary conditional behavior;
        traced calls use one scalar JAX conditional, including the calculations
        of diagnostic values. Both paths preserve the caller's sampling leaf.
        """
        if not self.names:
            return self
        if isinstance(self.l_sample, bool):
            return accumulate(self) if self.l_sample else self

        def accumulate_sample(stats):
            """Trace the sampled branch with nested sampling checks resolved."""
            updated = accumulate(stats._new(l_sample=True))
            return updated._new(l_sample=stats.l_sample)

        return self.sampled_value(accumulate_sample, self)

    def sampled_value(self, compute, previous, *, l_required: bool = False):
        """Compute a sampled or required field, otherwise return its prior value.

        ``previous`` supplies the shape and dtype of the unchanged branch. Use
        this for diagnostic arrays needed by a later stats block; stats-only
        updates should use ``sampled`` to resolve nested sampling checks. Set the
        static ``l_required`` flag when physics needs the same field even on
        unsampled steps.
        """
        if l_required:
            return compute(previous)
        if isinstance(self.l_sample, bool):
            return compute(previous) if self.l_sample else previous
        return jax.lax.cond(self.l_sample, compute, lambda value: value, previous)

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
        l_count_sample: bool,
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
        arr = jnp.asarray(values, dtype=core_rknd)
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
            count_candidate = count_bank.at[slot, row, column].add(int(l_count_sample))
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
                count_candidate = count_bank.at[slot, :, 0].add(int(l_count_sample))
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
                count_candidate = count_bank.at[slot, row, :].add(int(l_count_sample))
        else:
            if icol is not None or level is not None:
                raise ValueError("Rank-2 stats operations do not accept icol or level.")
            if arr.shape != (self.ncol, nlev):
                raise ValueError(
                    f"Rank-2 stats update for {name!r} requires "
                    f"({self.ncol}, {nlev}); got {arr.shape}."
                )
            candidate = bank.at[slot, :, :].add(sign * arr)
            count_candidate = count_bank.at[slot, :, :].add(int(l_count_sample))

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
        if not clean_name or clean_name not in self.name_to_id:
            return self

        def accumulate(stats):
            """Accumulate the requested variable on a sampled timestep."""
            var_id = stats.name_to_id[clean_name]
            allowed = jnp.logical_not(stats.l_budget[var_id])
            issues = stats._record_issue(
                stats.issues, var_id, jnp.logical_not(allowed), ISSUE_UPDATE_AFTER_BUDGET
            )
            buffers, nsamples, issues = stats._apply_payload(
                name=clean_name,
                values=values,
                sign=1,
                l_count_sample=True,
                allowed=allowed,
                icol=icol,
                level=level,
                buffers=stats.buffers,
                nsamples=stats.nsamples,
                issues=issues,
            )
            return stats._new(buffers=buffers, nsamples=nsamples, issues=issues)

        return self.sampled(accumulate)

    def average_subtimesteps(self, previous, num_subtimesteps):
        """Finish source stats_update(sub_timestep_average=.true.) calls.

        JAX adaptation: average only the sample-loop increments, preserving
        earlier timestep/window statistics and counting one aggregate update.
        """
        buffers = tuple(
            old + (new - old) / num_subtimesteps for old, new in zip(previous.buffers, self.buffers)
        )
        nsamples = tuple(
            old + (new - old) // num_subtimesteps for old, new in zip(previous.nsamples, self.nsamples)
        )
        return self._new(buffers=buffers, nsamples=nsamples)

    def begin_budget(self, name: str, values, *, icol: int | None = None):
        """Mirror ``stats_begin_budget``."""
        clean_name = name.strip()
        if not clean_name or clean_name not in self.name_to_id:
            return self

        def accumulate(stats):
            """Accumulate the requested variable on a sampled timestep."""
            var_id = stats.name_to_id[clean_name]
            allowed = jnp.logical_not(stats.l_in_budget[var_id])
            issues = stats._record_issue(
                stats.issues, var_id, jnp.logical_not(allowed), ISSUE_BEGIN_TWICE
            )
            buffers, nsamples, issues = stats._apply_payload(
                name=clean_name,
                values=values,
                sign=-1,
                l_count_sample=False,
                allowed=allowed,
                icol=icol,
                level=None,
                buffers=stats.buffers,
                nsamples=stats.nsamples,
                issues=issues,
            )
            l_budget = stats.l_budget.at[var_id].set(jnp.where(allowed, True, stats.l_budget[var_id]))
            l_in_budget = stats.l_in_budget.at[var_id].set(
                jnp.where(allowed, True, stats.l_in_budget[var_id])
            )
            return stats._new(
                buffers=buffers,
                nsamples=nsamples,
                l_budget=l_budget,
                l_in_budget=l_in_budget,
                issues=issues,
            )

        return self.sampled(accumulate)

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
        if not clean_name or clean_name not in self.name_to_id:
            return self

        def accumulate(stats):
            """Accumulate the requested variable on a sampled timestep."""
            var_id = stats.name_to_id[clean_name]
            allowed = stats.l_in_budget[var_id]
            issues = stats._record_issue(
                stats.issues, var_id, jnp.logical_not(allowed), ISSUE_UPDATE_WITHOUT_BEGIN
            )
            buffers, nsamples, issues = stats._apply_payload(
                name=clean_name,
                values=values,
                sign=1,
                l_count_sample=False,
                allowed=allowed,
                icol=icol,
                level=level,
                buffers=stats.buffers,
                nsamples=stats.nsamples,
                issues=issues,
            )
            return stats._new(buffers=buffers, nsamples=nsamples, issues=issues)

        return self.sampled(accumulate)

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
        if not clean_name or clean_name not in self.name_to_id:
            return self

        def accumulate(stats):
            """Accumulate the requested variable on a sampled timestep."""
            var_id = stats.name_to_id[clean_name]
            allowed = stats.l_in_budget[var_id]
            issues = stats._record_issue(
                stats.issues, var_id, jnp.logical_not(allowed), ISSUE_FINALIZE_WITHOUT_BEGIN
            )
            buffers, nsamples, issues = stats._apply_payload(
                name=clean_name,
                values=values,
                sign=1,
                l_count_sample=bool(l_count_sample),
                allowed=allowed,
                icol=icol,
                level=None,
                buffers=stats.buffers,
                nsamples=stats.nsamples,
                issues=issues,
            )
            l_in_budget = stats.l_in_budget.at[var_id].set(
                jnp.where(allowed, False, stats.l_in_budget[var_id])
            )
            return stats._new(
                buffers=buffers,
                nsamples=nsamples,
                l_in_budget=l_in_budget,
                issues=issues,
            )

        return self.sampled(accumulate)

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
                warnings.warn(f"stats warning: budget left mid-update at end_timestep for {name}", RuntimeWarning)

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
