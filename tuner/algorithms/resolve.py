"""Cover a bounded search region with an evenly spaced Cartesian grid.

ResolveGridStrategy evaluates every combination of one-dimensional grid values.
It requires loss values for ranking but never adapts its proposals to them.
For bounds l_j <= u_j and requested maximum spacing h > 0, coordinate j has
n_j = ceil((u_j-l_j)/h) intervals and n_j+1 points. Nonfixed coordinates use
p[j,k] = l_j + (u_j-l_j)*k/n_j, k=0,...,n_j. Fixed ranges contribute one point.
The actual spacing is therefore at most h, and both endpoints are included.
For [0,1] and h=0.3, the points are 0, 0.25, 0.5, 0.75, 1.

Usage
-----
Construct with names, defaults, logical ranges and optional spacing. Omitted
spacing is half the smallest nonzero range width, in physical parameter units.
For ranges [0,2] and [10,16], h=1 gives 3*7=21 grid points. Fixed ranges do
not constrain the spacing; an entirely fixed search evaluates one point.
prepare_options derives this default for requests and previews; construction
also supports omitted spacing directly. fill(queue, capacity)
advances the Cartesian iterator with the last coordinate changing fastest.
Inherited tell ignores scores. is_exhausted reports when iteration ends;
estimated_sample_count and sample_count give the product of coordinate counts.
The seed argument is accepted by the shared constructor contract but unused.
_values_for_range and _point_count implement the coordinate grid/count policy.

Only coordinate arrays and the Cartesian iterator are retained; the complete
multidimensional list is not materialized. Work still grows as the product of
coordinate counts, so a fine grid becomes impractical in many dimensions.
BaseTuningStrategy expands linked coordinates and preserves untuned defaults.
Current checkpoint serialization uses the iterator's native pickle support;
Python-version restrictions on that iterator can limit continuation.
"""

from __future__ import annotations
from collections import deque
from typing import Iterable
import itertools
import math

from tuner.algorithms.base import BaseTuningStrategy


class ResolveGridStrategy(BaseTuningStrategy):
    """Deterministic full-space sampler using evenly resolved parameter ranges."""

    def __init__(
        self,
        *,
        param_names: Iterable[str],
        default_params_row: Iterable[float],
        parameter_ranges: Iterable[dict],
        spacing: float | None = None,
        seed: int | None = None,
    ):
        """Construct coordinate grids and a lazy product iterator using physical-unit spacing.
        """
        super().__init__(
            param_names=param_names,
            default_params_row=default_params_row,
            parameter_ranges=parameter_ranges,
        )
        self.spacing = self._resolve_spacing(spacing, self.parameter_ranges)
        self.grid_names = [spec["name"] for spec in self.parameter_ranges]
        self.grid_values = [
            self._values_for_range(float(spec["min"]), float(spec["max"]), self.spacing)
            for spec in self.parameter_ranges
        ]
        self.total_samples = math.prod(len(values) for values in self.grid_values)
        self._grid_iter = itertools.product(*self.grid_values)
        self._exhausted = self.total_samples == 0

    @staticmethod
    def _resolve_spacing(spacing, ranges) -> float:
        """Use explicit spacing or half the smallest nonfixed physical range width.

        Fixed-only and empty searches have no intervals to resolve: unit spacing
        leaves their single Cartesian point unchanged. Reject invalid explicit
        spacing and bounds before counting or allocating grid values.
        """
        spans = [float(spec["max"]) - float(spec["min"]) for spec in ranges]
        if any(not math.isfinite(span) or span < 0.0 for span in spans):
            raise ValueError("resolve ranges require finite min <= max")
        if spacing is None:
            active = [span for span in spans if span > 0.0]
            spacing = min(active) / 2.0 if active else 1.0
        spacing = float(spacing)
        if not math.isfinite(spacing) or spacing <= 0.0:
            raise ValueError("spacing must be finite and > 0")
        return spacing

    @staticmethod
    def prepare_options(options, request):
        """Save explicit or range-derived spacing before budgeting and construction."""
        options["spacing"] = ResolveGridStrategy._resolve_spacing(
            options.get("spacing"), request.get("parameter_ranges", []))

    @staticmethod
    def _values_for_range(min_value: float, max_value: float, spacing: float) -> list[float]:
        """Return one coordinate grid, including both bounds and a single value when fixed.

        The interval count rounds upward, so the resulting uniform spacing never
        exceeds the requested positive maximum.
        """
        point_count = ResolveGridStrategy._point_count(min_value, max_value, spacing)
        if point_count == 1:
            return [float(min_value)]
        span = float(max_value) - float(min_value)
        intervals = point_count - 1
        values = [
            float(min_value + span * idx / intervals)
            for idx in range(intervals)
        ]
        values.append(float(max_value))
        return values

    @staticmethod
    def _point_count(min_value: float, max_value: float, spacing: float) -> int:
        """Count a coordinate grid without allocating values; reject inverted bounds.
        """
        span = float(max_value) - float(min_value)
        if span < 0.0:
            raise ValueError("resolve ranges require min <= max")
        return 1 if span == 0.0 else max(1, int(math.ceil(span / spacing))) + 1

    def fill(self, pending_samples: deque, capacity: int) -> None:
        """Emit consecutive Cartesian points up to queue capacity.

        StopIteration marks exhaustion. A partially consumed product is kept between
        calls, so different batch capacities preserve the same proposal order.
        """
        capacity = int(capacity)
        while len(pending_samples) < capacity and not self._exhausted:
            try:
                combo = next(self._grid_iter)
            except StopIteration:
                self._exhausted = True
                return
            selected_params = {
                name: float(value)
                for name, value in zip(self.grid_names, combo)
            }
            pending_samples.append(self._sample_from_selected_params(selected_params))

    def is_exhausted(self) -> bool:
        """Return whether the Cartesian iterator has reached its end.
        """
        return self._exhausted

    def estimated_sample_count(self) -> int | None:
        """Return the number of candidate evaluations, excluding separately scored baselines.
        """
        return self.total_samples

    @staticmethod
    def sample_count(options, ranges):
        """Multiply coordinate counts without allocating any grid or product points.
        """
        ranges = list(ranges)
        spacing = ResolveGridStrategy._resolve_spacing(options.get("spacing"), ranges)
        return math.prod(ResolveGridStrategy._point_count(
            float(spec["min"]), float(spec["max"]), spacing) for spec in ranges)
