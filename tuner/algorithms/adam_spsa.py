"""Minimize a scalar objective with Adam using paired loss probes instead of derivatives.

AdamSPSAStrategy implements simultaneous perturbation stochastic approximation
(SPSA): each gradient estimate perturbs every active normalized coordinate at
once. It needs only total_loss, so it works with Fortran or JAX value evaluation
and with synthetic objectives. No model gradient is requested or consumed.

Coordinates, probes and update equations
---------------------------------------
Physical values satisfy p_j = lower_j + span_j*x_j with x_j in [0,1]. For each
chain and update, draw K=spsa_pairs independent directions Delta_k whose active
entries are equally likely -1 or +1 (fixed entries are zero). At perturbation
c, evaluate F at x_k^- = clip(x-c*Delta_k) and x_k^+ = clip(x+c*Delta_k).
Average the usable paired estimates:

    ghat = mean_k [(F(x_k^+) - F(x_k^-))/(2*c) * Delta_k].

Here F is the already combined scalar loss. The denominator remains 2*c even
when bounds clip a probe; estimates near bounds can therefore be biased.
A pair with a nonfinite loss is omitted. If no pairs are usable, ghat is zero;
the update still advances and existing momentum can still move the center.

With t starting at one, update m=0.9*m+0.1*ghat and
v=0.999*v+0.001*ghat**2, divide each by its bias correction (1-beta**t), then
x = clip(x-learning_rate*mhat/(sqrt(vhat)+1e-8), 0, 1). Fixed coordinates remain
zero in normalized storage. Bounds projection does not clear moments.
learning_rate and perturbation both use fractions of each parameter's range;
linked physical targets share one perturbed logical coordinate.

Scheduling and starts
---------------------
Each chain follows initial center -> K probe pairs -> update -> ... -> final
center. Intermediate updated centers are not separately scored; their nearby
probe losses and the initial/final centers remain candidates for ranking.
T updates cost 2 + 2*K*T evaluations per chain, compared with T+1 for gradient
Adam. The initial center is scored even though SPSA needs only the probe pairs.
The first chain uses clipped defaults; additional starts use the shared Latin
hypercube sampler. Separate child random streams control starts and directions.

prepare_options requires batch_size divisible by 2*K. A batch holds
chains_per_batch=batch_size/(2*K), and concurrent_batches=ceil(workers/cases)
cohorts produce chain_count=chains_per_batch*concurrent_batches. A cohort waits
for every member to be ready before enqueueing its complete probe block.
This scheduling differs from gradient Adam's workers*batch_size chain count.
Sample-ID context keeps asynchronous probe results attached to their chain,
pair and sign. All-fixed ranges collapse to one evaluation and one chain.

Usage
-----
Construct with names, defaults, ranges, seed, user options and the derived
scheduling values from prepare_options. fill queues centers/probes; tell
records paired total_loss values and advances ready chains. is_exhausted,
estimated_sample_count and pending_batch_limit expose progress/budgets; static
sample_count gives a pre-construction estimate. Private helpers implement
starts, proposal expansion, loss validation and moment updates. _AdamChain
stores each phase, moments and unfinished probe block for pickle continuation.
The scheduler retains the best scored candidate independently of the chain;
individual steps need not improve and no global-minimum guarantee is implied.

References: Spall (1992), IEEE Transactions on Automatic Control 37(3), 332–341
(SPSA); Kingma and Ba (2014), arXiv:1412.6980 (Adam).
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
import math
from typing import Iterable

import numpy as np

from tuner.algorithms.base import BaseTuningStrategy, latin_hypercube_points


@dataclass
class _AdamChain:
    chain_id: int
    x: np.ndarray
    m: np.ndarray
    v: np.ndarray
    phase: str = "initial_ready"
    update_count: int = 0
    pending_sample_ids: set[int] = field(default_factory=set)
    probe_losses: dict[tuple[int, int], float] = field(default_factory=dict)
    directions: list[np.ndarray] = field(default_factory=list)


class AdamSPSAStrategy(BaseTuningStrategy):
    """Run independent projected Adam chains from batched SPSA gradients."""

    beta1 = 0.9
    beta2 = 0.999
    epsilon = 1.0e-8

    def __init__(
        self,
        *,
        param_names: Iterable[str],
        default_params_row: Iterable[float],
        parameter_ranges: Iterable[dict],
        max_updates: int,
        learning_rate: float,
        perturbation: float,
        spsa_pairs: int,
        chains_per_batch: int,
        concurrent_batches: int,
        chain_count: int,
        seed: int | None = None,
    ):
        """Initialize projected Adam chains, independent random streams and complete probe cohorts.

        Ranges are logical search coordinates; defaults and full proposal rows use
        physical parameter units. Derived cohort sizes must come from prepare_options.
        """
        super().__init__(
            param_names=param_names,
            default_params_row=default_params_row,
            parameter_ranges=parameter_ranges,
        )
        self.max_updates = int(max_updates)
        self.learning_rate = float(learning_rate)
        self.perturbation = float(perturbation)
        self.spsa_pairs = int(spsa_pairs)
        self.chains_per_batch = int(chains_per_batch)
        self.concurrent_batches = int(concurrent_batches)
        self.chain_count = int(chain_count)
        self.dimension = len(self.parameter_ranges)
        self.minimum = np.asarray([float(spec["min"]) for spec in self.parameter_ranges])
        self.maximum = np.asarray([float(spec["max"]) for spec in self.parameter_ranges])
        self.span = self.maximum - self.minimum
        self.active = self.span > 0.0
        if not self.active.any():
            self.chain_count = 1
            self.chains_per_batch = 1
            self.concurrent_batches = 1

        seed_sequence = np.random.SeedSequence(seed)
        lhs_seed, direction_seed = seed_sequence.spawn(2)
        self.rng = np.random.default_rng(direction_seed)
        starts = self._starting_points(lhs_seed)
        self.chains = [
            _AdamChain(
                chain_id=chain_id,
                x=starts[chain_id],
                m=np.zeros(self.dimension, dtype=float),
                v=np.zeros(self.dimension, dtype=float),
            )
            for chain_id in range(self.chain_count)
        ]
        self.cohorts = [
            self.chains[index:index + self.chains_per_batch]
            for index in range(0, len(self.chains), self.chains_per_batch)
        ]
        self._pending_context: dict[int, dict] = {}

    def fill(self, pending_samples: deque, capacity: int) -> None:
        """Queue initial centers, complete ready probe cohorts, then final centers.

        A probe block is emitted only when queue capacity can hold its full batch
        width. Outstanding context prevents duplicate proposals while feedback waits.
        """
        capacity = int(capacity)
        for chain in self.chains:
            if chain.phase == "initial_ready" and len(pending_samples) < capacity:
                self._append_center(pending_samples, chain, "initial")

        batch_width = 2 * self.spsa_pairs * self.chains_per_batch
        for cohort in self.cohorts:
            if not all(chain.phase == "probe_ready" for chain in cohort):
                continue
            if len(pending_samples) + batch_width > capacity:
                return
            for chain in cohort:
                self._append_probes(pending_samples, chain)

        for chain in self.chains:
            if chain.phase == "final_ready" and len(pending_samples) < capacity:
                self._append_center(pending_samples, chain, "final")

    def tell(self, completed_samples: list[dict]) -> None:
        """Attach losses to their chain/pair/sign and update only after its whole probe block arrives.

        Initial/final center feedback changes phase without supplying an SPSA
        gradient. Sorting IDs makes processing deterministic within a delivered group.
        """
        for completed in sorted(completed_samples, key=lambda item: int(item["sample_id"])):
            context = self._pending_context.pop(int(completed["sample_id"]), None)
            if context is None:
                continue
            chain = self.chains[int(context["chain_id"])]
            chain.pending_sample_ids.discard(int(completed["sample_id"]))
            kind = context["kind"]
            if kind == "initial":
                chain.phase = "probe_ready" if self.active.any() else "done"
            elif kind == "final":
                chain.phase = "done"
            else:
                probe_key = (int(context["pair"]), int(context["sign"]))
                chain.probe_losses[probe_key] = self._loss(completed)
                if not chain.pending_sample_ids:
                    self._update_chain(chain)

    def is_exhausted(self) -> bool:
        """Return whether every chain has completed its final scored center.
        """
        return all(chain.phase == "done" for chain in self.chains)

    def estimated_sample_count(self) -> int | None:
        """Count initial/final centers and 2*spsa_pairs probes per update; fixed ranges cost one.
        """
        if not self.active.any():
            return 1
        return self.chain_count * (2 + 2 * self.spsa_pairs * self.max_updates)

    def pending_batch_limit(self, _batch_size: int) -> int | None:
        """Return the derived number of independently active probe cohorts.
        """
        return self.concurrent_batches

    def _starting_points(self, seed_sequence: np.random.SeedSequence) -> list[np.ndarray]:
        """Build clipped-default and Latin hypercube centers using only the initialization stream.
        """
        default = np.zeros(self.dimension, dtype=float)
        for index, spec in enumerate(self.parameter_ranges):
            if self.active[index]:
                value = self.default_params_row[self.param_index[spec["name"]]]
                default[index] = np.clip((value - self.minimum[index]) / self.span[index], 0.0, 1.0)
        starts = [default]
        if self.chain_count > 1:
            active_count = int(np.count_nonzero(self.active))
            sample_count = self.chain_count - 1
            rng = np.random.default_rng(seed_sequence)
            samples = latin_hypercube_points(rng, sample_count, active_count)
            for sample in samples:
                point = np.zeros(self.dimension, dtype=float)
                point[self.active] = sample
                starts.append(point)
        return starts

    def _append_center(self, pending_samples: deque, chain: _AdamChain, kind: str) -> None:
        """Emit a physical initial/final center and register its phase and pending sample ID.
        """
        sample = self._sample_from_x(chain.x)
        sample_id = int(sample["sample_id"])
        self._pending_context[sample_id] = {"kind": kind, "chain_id": chain.chain_id}
        chain.pending_sample_ids = {sample_id}
        chain.phase = f"{kind}_waiting"
        pending_samples.append(sample)

    def _append_probes(self, pending_samples: deque, chain: _AdamChain) -> None:
        """Draw independent sign directions and emit both clipped probes for each pair.

        Reset per-update probe state, attach sample IDs to pair/sign context, and
        wait for all returned losses before changing the chain center.
        """
        chain.directions = []
        chain.probe_losses = {}
        chain.pending_sample_ids = set()
        for pair in range(self.spsa_pairs):
            direction = np.zeros(self.dimension, dtype=float)
            direction[self.active] = self.rng.choice(
                (-1.0, 1.0), size=int(np.count_nonzero(self.active))
            )
            chain.directions.append(direction)
            for sign in (-1, 1):
                point = np.clip(chain.x + sign * self.perturbation * direction, 0.0, 1.0)
                sample = self._sample_from_x(point)
                sample_id = int(sample["sample_id"])
                self._pending_context[sample_id] = {
                    "kind": "probe",
                    "chain_id": chain.chain_id,
                    "pair": pair,
                    "sign": sign,
                }
                chain.pending_sample_ids.add(sample_id)
                pending_samples.append(sample)
        chain.phase = "probe_waiting"

    def _update_chain(self, chain: _AdamChain) -> None:
        """Average finite paired differences and advance bias-corrected Adam moments.

        The module equations define the normalized estimate, clipped-probe bias and
        zero-valid-pair behavior. Finish by selecting final_ready or probe_ready.
        """
        gradient = np.zeros(self.dimension, dtype=float)
        valid_pairs = 0
        for pair, direction in enumerate(chain.directions):
            loss_minus = chain.probe_losses.get((pair, -1), math.inf)
            loss_plus = chain.probe_losses.get((pair, 1), math.inf)
            if not math.isfinite(loss_minus) or not math.isfinite(loss_plus):
                continue
            gradient += (loss_plus - loss_minus) / (2.0 * self.perturbation) * direction
            valid_pairs += 1
        if valid_pairs:
            gradient /= valid_pairs

        chain.update_count += 1
        chain.m = self.beta1 * chain.m + (1.0 - self.beta1) * gradient
        chain.v = self.beta2 * chain.v + (1.0 - self.beta2) * gradient**2
        m_hat = chain.m / (1.0 - self.beta1**chain.update_count)
        v_hat = chain.v / (1.0 - self.beta2**chain.update_count)
        chain.x = np.clip(
            chain.x - self.learning_rate * m_hat / (np.sqrt(v_hat) + self.epsilon),
            0.0,
            1.0,
        )
        chain.x[~self.active] = 0.0
        chain.phase = "final_ready" if chain.update_count >= self.max_updates else "probe_ready"

    def _sample_from_x(self, point: np.ndarray) -> dict:
        """Expand normalized logical coordinates into full physical rows with linked targets.
        """
        physical = self.minimum + self.span * point
        return self._sample_from_selected_params(
            {
                spec["name"]: float(physical[index])
                for index, spec in enumerate(self.parameter_ranges)
            }
        )

    @staticmethod
    def _loss(completed: dict) -> float:
        """Read a finite scalar total_loss; unusable values become infinity and invalidate that pair.
        """
        try:
            loss = float(completed["total_loss"])
        except (KeyError, TypeError, ValueError):
            return math.inf
        return loss if math.isfinite(loss) else math.inf

    @staticmethod
    def prepare_options(options, request):
        """Derive probe cohorts from batch columns, workers and case count; reject partial-pair widths.
        """
        columns_per_chain = 2 * options["spsa_pairs"]
        if request["batch_size"] % columns_per_chain:
            raise ValueError(
                "Adam requires batch_size to be divisible by 2 * spsa_pairs "
                f"({request['batch_size']} is not divisible by {columns_per_chain})"
            )
        chains_per_batch = request["batch_size"] // columns_per_chain
        concurrent_batches = math.ceil(request["max_workers"] / len(request["cases"]))
        options.update(chains_per_batch=chains_per_batch, concurrent_batches=concurrent_batches,
                       chain_count=chains_per_batch * concurrent_batches)

    @staticmethod
    def sample_count(options, _ranges):
        """Count 2+2*spsa_pairs*max_updates candidates per derived chain before construction.
        """
        return options["chain_count"] * (2 + 2 * options["spsa_pairs"] * options["max_updates"])
