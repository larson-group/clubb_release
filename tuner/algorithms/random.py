"""Explore bounded parameter ranges with independent uniform samples.

RandomUniformStrategy is a value-only reference sampler: it proposes points
without learning from their losses. For logical coordinate j, each candidate
uses p_j = lower_j + U_j * (upper_j-lower_j), with independent U_j uniform in
[0,1]. The standard-library generator can include an endpoint through floating
rounding. Fixed ranges always return their fixed value. Independent draws can
cluster; unlike Latin hypercube sampling they do not guarantee interval coverage.

Usage
-----
Construct with parameter names, defaults, ranges, optional seed and max_samples
(default 200). fill(queue, capacity) proposes additional rows; inherited tell
ignores feedback. is_exhausted and estimated_sample_count report the finite
budget; max_samples=None permits unlimited direct/programmatic sampling.
sample_count reports the budget to the registry before construction.

BaseTuningStrategy assigns monotonically increasing IDs, expands linked ranges
(for example C6rt=C6thl) and retains defaults for untuned physical parameters.
No starting default point is injected; baseline scoring is a scheduler task.
The RNG and emitted count survive current-format pickling, so batching or
continuation does not restart the sequence. No gradients, model objects or
objective callback are stored. The scheduler can evaluate any scalar objective
and send its results to tell without changing this sampler.
"""

from __future__ import annotations
from collections import deque
from typing import Iterable
import random

from tuner.algorithms.base import BaseTuningStrategy


class RandomUniformStrategy(BaseTuningStrategy):
    """Uniform random sampler over selected tunable-parameter ranges."""

    def __init__(
        self,
        *,
        param_names: Iterable[str],
        default_params_row: Iterable[float],
        parameter_ranges: Iterable[dict],
        max_samples: int | None = 200,
        seed: int | None = None,
    ):
        """Start a reproducible uniform sequence and optional evaluation limit.

        Defaults and ranges use physical parameter units. The seed affects only
        proposals; changing worker or batch counts does not change the RNG sequence.
        """
        super().__init__(
            param_names=param_names,
            default_params_row=default_params_row,
            parameter_ranges=parameter_ranges,
        )
        self.max_samples = None if max_samples is None else int(max_samples)
        self.rng = random.Random(seed)

    def fill(self, pending_samples: deque, capacity: int) -> None:
        """Append bounded full parameter rows until capacity or the sample budget is reached.

        The queue may already contain proposals. Each new point consumes one RNG draw
        per logical range and one unique sample ID; it needs no completed feedback.
        """
        capacity = int(capacity)
        while len(pending_samples) < capacity and not self.is_exhausted():
            selected_params = {}
            for spec in self.parameter_ranges:
                name = spec["name"]
                selected_params[name] = self.rng.uniform(float(spec["min"]), float(spec["max"]))
            pending_samples.append(self._sample_from_selected_params(selected_params))

    def is_exhausted(self) -> bool:
        """Return whether the finite proposal budget is used; unlimited sampling stays active.
        """
        return self.max_samples is not None and self.next_sample_id >= self.max_samples

    def estimated_sample_count(self) -> int | None:
        """Return the number of candidate evaluations, excluding separately scored baselines.
        """
        return self.max_samples

    @staticmethod
    def sample_count(options, _ranges):
        """Read the pre-construction budget from normalized catalog options.
        """
        return options.get("max_samples")
