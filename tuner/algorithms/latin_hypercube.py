"""Cover parameter ranges with a fixed Latin hypercube ensemble.

LatinHypercubeStrategy proposes ``num_points`` candidates without using loss
feedback. Unlike independent uniform sampling, every parameter visits each of
N equal-width intervals exactly once, where N is the total ensemble size.
For coordinate j and ensemble member i, the normalized sample is

    x[i,j] = (permutation_j[i] + uniform[i,j]) / N,
    p[i,j] = lower[j] + (upper[j] - lower[j]) * x[i,j].

The independent permutations randomize which intervals are paired across
coordinates; the uniform offsets randomize positions within those intervals.
This stratifies each individual coordinate, not every multidimensional cell.
For example, ten points over C8=[0.2,0.8] put one C8 value in each interval of
width 0.06. The shared seed controls permutations and offsets.

Usage
-----
The catalog registers ``latin_hypercube:num_points=200`` for CLI and Dash.
Construct the strategy with parameter names, defaults, ranges, seed and point
count. ``fill(queue, capacity)`` emits the next ensemble members; inherited
``tell(completed)`` ignores their scores. ``is_exhausted`` and
``estimated_sample_count`` describe the finite ensemble. ``sample_count``
reports its size before construction. There is no annealing or refinement
stage: ranking and retaining the best losses belong to the tuner.

The entire normalized design is drawn once, so batch size and feedback order
cannot restart its stratification. It uses O(N*d) memory for d logical ranges.
Fixed ranges stay constant; linked physical parameters share their logical
coordinate. Untuned parameters retain defaults through BaseTuningStrategy.
No default point is injected into the ensemble because that would disturb the
Latin hypercube property. The tuner scores its default baseline separately.
The stored design and cursor support ordinary current-format checkpoints.

Reference: McKay, Beckman and Conover (1979), Technometrics 21(2), 239–245.
"""

from collections import deque

import numpy as np

from tuner.algorithms.base import BaseTuningStrategy, latin_hypercube_points


class LatinHypercubeStrategy(BaseTuningStrategy):
    """Propose a seeded stratified ensemble with one sample in every marginal stratum."""

    def __init__(self, *, param_names, default_params_row, parameter_ranges,
                 num_points=200, seed=None):
        """Draw one complete normalized design and start its emission cursor at zero."""
        super().__init__(param_names=param_names, default_params_row=default_params_row,
                         parameter_ranges=parameter_ranges)
        self.num_points = int(num_points)
        if self.num_points < 1:
            raise ValueError("num_points must be >= 1")
        self.points = latin_hypercube_points(np.random.default_rng(seed), self.num_points,
                                             len(self.parameter_ranges))

    def fill(self, pending_samples: deque, capacity: int) -> None:
        """Emit consecutive design rows up to queue capacity, expanding linked targets."""
        while len(pending_samples) < int(capacity) and not self.is_exhausted():
            point = self.points[self.next_sample_id]
            selected = {
                spec["name"]: spec["min"] + (spec["max"] - spec["min"]) * coordinate
                for spec, coordinate in zip(self.parameter_ranges, point)
            }
            pending_samples.append(self._sample_from_selected_params(selected))

    def is_exhausted(self) -> bool:
        """Return whether all requested ensemble members have been proposed."""
        return self.next_sample_id >= self.num_points

    def estimated_sample_count(self) -> int:
        """Return the total number of ensemble evaluations, including fixed ranges."""
        return self.num_points

    @staticmethod
    def sample_count(options, _ranges):
        """Return the catalog-validated ensemble size before constructing the design."""
        return options["num_points"]
