"""Define the objective-independent proposal/feedback contract for tuning algorithms.

An algorithm chooses parameter vectors; a caller evaluates its own scalar
objective. No CLUBB case, loss function or JAX runtime is embedded in this
interface. For example, the same strategy can receive a BOMEX/ATEX combined
model loss or the value of f(x,y)=x*x+y*y from a test harness.

Usage and sample layout
-----------------------
Subclass BaseTuningStrategy and implement fill, tell when adaptive, and
is_exhausted. Construction receives param_names (the full physical parameter
order), default_params_row (one value per name), and parameter_ranges (logical
search coordinates with name, min, max and optional targets). A range such as
{name: 'C6rt', targets: ['C6rt','C6thl'], min: 0, max: 4} drives both physical
parameters together. _sample_from_selected_params expands logical values into
those physical targets, retains every untuned default and assigns a unique
sample_id. Returned param_row follows param_names; selected_params lists the
physical targets changed by the search; all_params names the full row.

fill(pending_samples, capacity) appends proposals until that queue reaches its
maximum length. The caller drains/evaluates proposals and returns scored sample
dictionaries to tell(completed_samples). Every result has sample_id and
scalar total_loss. A gradient method additionally requires total_gradient in
gradient_parameters physical-name order. Store adaptive context by sample ID:
feedback can arrive in partial groups or a different order from proposal order.
The caller combines case losses/gradients before feedback; algorithms should
not reconstruct the model objective from profile diagnostics.

Completion and scheduling
-------------------------
is_exhausted must remain false while an adaptive method still needs feedback
or an updated final point to be scored. estimated_sample_count reports the
candidate budget after construction; pending_batch_limit optionally limits
useful queued batches. prepare_options derives scheduling settings from a
request before construction; sample_count supplies a budget for previews.
Their inherited defaults make no scheduling changes and report no budget.
The whole strategy state must support pickle continuation in the current format.

latin_hypercube_points is the shared NumPy design generator for LHS and Adam
starts. Each coordinate independently permutes N equal strata and draws one
uniform offset per stratum: x[i,j]=(permutation_j[i]+U[i,j])/N. It returns an
(N,d) normalized array in [0,1), requires N>0 and supports d=0 empty rows.
NumPy loads only when the helper is called; metadata-only catalog discovery
therefore needs no scientific environment.
"""

from __future__ import annotations
from collections import deque
from typing import Iterable


def latin_hypercube_points(rng, count: int, dimensions: int):
    """Draw a normalized (count, dimensions) Latin hypercube using the supplied NumPy RNG.

    Independent coordinate permutations cover each equal-width stratum once.
    Uniform offsets randomize within strata. A positive count is required; zero
    dimensions returns count empty rows. The caller scales values to physical units.
    """
    import numpy as np

    strata = np.empty((count, dimensions), dtype=int)
    for dimension in range(dimensions):
        strata[:, dimension] = rng.permutation(count)
    return (strata + rng.random((count, dimensions))) / count


class BaseTuningStrategy:
    """Common helpers for sample-generating tuning strategies."""

    def __init__(
        self,
        *,
        param_names: Iterable[str],
        default_params_row: Iterable[float],
        parameter_ranges: Iterable[dict],
    ):
        """Validate full defaults and range target names, then start proposal IDs at zero.

        This base owns names, linking and row expansion. Catalog/request validation
        owns numerical option bounds; subclasses own their optimizer state.
        """
        self.param_names = list(param_names)
        self.default_params_row = [float(value) for value in default_params_row]
        self.param_index = {name: idx for idx, name in enumerate(self.param_names)}
        self.parameter_ranges = []
        for raw_spec in parameter_ranges:
            spec = dict(raw_spec)
            name = str(spec.get("name", "")).strip()
            raw_targets = spec.get("targets", [name])
            if isinstance(raw_targets, str):
                raw_targets = [raw_targets]
            targets = [str(target).strip() for target in raw_targets or [] if str(target).strip()]
            if not name or not targets:
                raise ValueError("Each parameter range requires a name and at least one target")
            spec["name"] = name
            spec["targets"] = targets
            self.parameter_ranges.append(spec)
        self.next_sample_id = 0

        if len(self.default_params_row) != len(self.param_names):
            raise ValueError("default_params_row length must match param_names length")

        missing = sorted(
            {
                target
                for spec in self.parameter_ranges
                for target in spec["targets"]
                if target not in self.param_index
            }
        )
        if missing:
            raise ValueError("Unknown tunable parameter(s): " + ", ".join(sorted(missing)))

    def _sample_from_selected_params(self, selected_params: dict[str, float]) -> dict:
        """Expand logical range values into linked physical targets and a complete parameter row.

        Input keys name logical ranges, not necessarily individual physical targets.
        Assign an incrementing ID and preserve untuned defaults; do not evaluate loss.
        """
        expanded = {}
        for spec in self.parameter_ranges:
            name = spec["name"]
            if name not in selected_params:
                continue
            value = float(selected_params[name])
            for target in spec["targets"]:
                expanded[target] = value
        row = list(self.default_params_row)
        for name, value in expanded.items():
            row[self.param_index[name]] = float(value)

        sample_id = self.next_sample_id
        self.next_sample_id += 1
        return {
            "sample_id": sample_id,
            "param_row": row,
            # Persist physical names so status, result reruns, and namelist
            # rendering reproduce every equality-constrained assignment.
            "selected_params": {name: float(value) for name, value in expanded.items()},
            "all_params": {name: float(row[idx]) for idx, name in enumerate(self.param_names)},
        }

    def fill(self, pending_samples: deque, capacity: int) -> None:
        """Append proposals until the pending queue reaches capacity; subclasses must implement this.
        """
        raise NotImplementedError

    def tell(self, _completed_samples: list[dict]) -> None:
        """Consume completed scalar-objective results; the default ignores feedback.

        Adaptive subclasses use sample_id and total_loss, optionally total_gradient
        and gradient_parameters, to update only the corresponding pending state.
        """
        return None

    def is_exhausted(self) -> bool:
        """Report whether proposal generation and required feedback are complete; default is unlimited.
        """
        return False

    def estimated_sample_count(self) -> int | None:
        """Return a known candidate budget after construction, or None for an unknown budget.
        """
        return None

    def pending_batch_limit(self, _batch_size: int) -> int | None:
        """Optionally cap in-flight batches when extra work cannot advance adaptive chains.
        """
        return None

    @staticmethod
    def prepare_options(options: dict, request: dict) -> None:
        """Derive execution-dependent options in place from the request; the default changes nothing.
        """
        return None

    @staticmethod
    def sample_count(options: dict, ranges: list[dict]) -> int | None:
        """Return a pre-construction candidate budget from normalized options and ranges, or None.
        """
        return None
