"""Search a bounded scalar objective with adaptive simulated-annealing chains.

EnhancedSimulatedAnnealingStrategy proposes points and accepts total_loss through
the same fill/tell interface as the other tuners. It follows the movement,
acceptance and stage adaptation of enhanced_simann::esa_driver. This Python
implementation manages independent asynchronous chains; the Fortran source is
the algorithm reference, not a claim of bitwise parity with its RNG sequence.

Coordinates and acceptance
--------------------------
Chain coordinates are physical logical-range values, not normalized Adam
coordinates. Initial locations form a seeded Latin hypercube across the ranges.
Each coordinate starts with step size 0.35*(upper-lower). A proposal changes a
selected subset of coordinates by a random signed fraction of their step and
clips each moved value to its bounds. A coordinate with no previous movement
attempt is forced into the selection; an empty selection becomes all coordinates.
The normal selection probability is 0.2 + 0.6*improvement_ratio[j].

Given old loss E and proposal E_try, delta=E_try-E. A non-increasing proposal
is accepted. An uphill proposal is accepted with probability exp(-delta/T),
where T is the chain temperature. This permits temporary increases to escape
local minima. xopt/nrgy_opt retain the best accepted point separately from the
current xstart/old_nrgy. The scheduler ranks all scored candidates separately.
Temperature therefore depends on the scale and offset of the objective;
there is no promise that a setting suited to one loss works for another.

Stage adaptation and stopping
------------------------------
A stage ends once total improvements reach 12*d or total attempts reach 50*d,
where d is the number of logical ranges. If accumulated stage loss is positive
and finite, the cooling factor is clamp(min_loss*attempts/sum_loss, 0.3, 0.7);
otherwise it is 0.7. The temperature is multiplied by that factor. For each
coordinate, ratio=improvements/attempts (zero when unattempted), and
step *= 1 + stp_adjst_factor*(ratio-stp_adjst_shift). With default shift=0.7
and factor=1, a ratio of 0.5 contracts a coordinate's step by 20%, while a
ratio of zero contracts it by 70%. This improves local resolution as successful
moves become less frequent; it can also limit exploration of other minima.
Stage counters clear and
search restarts from that chain's best point. f_tol is the relative improvement
threshold used to reset the count of stages without significant improvement.
The stopping check uses T<=max_final_temp and more than three such stages, or
the per-chain max_iters evaluation cap. The initial evaluation counts as one.

Usage
-----
Construct with names, defaults, physical ranges, seed and catalog options.
fill queues one ready proposal per chain; tell pairs returned total_loss with
its sample-ID context and accepts/rejects it. Feedback may be partial/reordered.
is_exhausted reports all chains done; estimated_sample_count and sample_count
report the maximum evaluation budget. prepare_options defaults chains to
workers*batch_size; an explicit SimAnn chain_count is supported.
pending_batch_limit bounds useful in-flight batches. Private helpers implement
starts, coordinate selection/movement, acceptance, stage cooling and summaries.
_SimAnnChain retains each temperature, counters and best/current positions for
current-format checkpoints. Unlike Adam, this method does not reserve a clipped
default start; the scheduler evaluates configured baselines separately.

Reference: historical CLUBB esa_driver and its cited simulated-annealing method:
https://github.com/larson-group/clubb/blob/73f0f6ad1b8c1b1e31ad780aae2d2c853a0674e5/src/Numerical_recipes/enhanced_simann.F90
See that revision for exact
acceptance and adaptation policy.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import math
import random
from typing import Iterable

from tuner.algorithms.base import BaseTuningStrategy


@dataclass
class _SimAnnChain:
    """One independent ``esa_driver`` state."""

    chain_id: int
    xinit: list[float]
    xopt: list[float]
    xstart: list[float]
    stp_cur: list[float]
    improved: list[float]
    attempted: list[float]
    improvement_ratio: list[float]
    temp: float
    iter: int = 1
    stages_w_no_improve: int = 0
    tot_nrgy: float = 0.0
    einit: float | None = None
    nrgy_opt: float | None = None
    old_nrgy: float | None = None
    min_nrgy: float | None = None
    initialized: bool = False
    exhausted: bool = False
    pending_sample_id: int | None = None


class EnhancedSimulatedAnnealingStrategy(BaseTuningStrategy):
    """Adaptive sampler modeled on ``enhanced_simann::esa_driver``."""

    no_improve_max = 3
    use_max = 0.8
    use_min = 0.2
    min_cool = 0.7
    max_cool = 0.3
    n1 = 12
    n2 = 50

    def __init__(
        self,
        *,
        param_names: Iterable[str],
        default_params_row: Iterable[float],
        parameter_ranges: Iterable[dict],
        initial_temp: float = 1.0,
        max_final_temp: float = 1.0e-12,
        max_iters: int = 200,
        chain_count: int = 1,
        stp_adjst_shift: float = 0.7,
        stp_adjst_factor: float = 1.0,
        f_tol: float = 1.0e-4,
        seed: int | None = None,
    ):
        """Initialize independent physical-coordinate chains and annealing controls.

        Each chain starts at a Latin hypercube point with a 35%-of-range proposal
        step. The RNG is shared by this strategy and preserved by checkpoints.
        """
        super().__init__(
            param_names=param_names,
            default_params_row=default_params_row,
            parameter_ranges=parameter_ranges,
        )
        self.rng = random.Random(seed)
        self.max_iters = int(max_iters)
        self.chain_count = max(1, int(chain_count))
        self.initial_temp = float(initial_temp)
        self.max_final_temp = float(max_final_temp)
        self.stp_adjst_shift = float(stp_adjst_shift)
        self.stp_adjst_factor = float(stp_adjst_factor)
        self.f_tol = float(f_tol)

        self.x_names = [spec["name"] for spec in self.parameter_ranges]
        self.xmin = [float(spec["min"]) for spec in self.parameter_ranges]
        self.xmax = [float(spec["max"]) for spec in self.parameter_ranges]
        self.vars = len(self.x_names)
        self.stp_ini = [
            (self.xmax[idx] - self.xmin[idx]) * 0.35
            for idx in range(self.vars)
        ]

        self.chains = [
            self._new_chain(chain_id, xinit)
            for chain_id, xinit in enumerate(self._latin_hypercube_points(self.chain_count))
        ]
        self._pending_context: dict[int, dict] = {}
        self._next_chain_index = 0

        self.xopt: list[float] | None = None
        self.nrgy_opt: float | None = None
        self.iter = 1
        self.temp = self.initial_temp
        self._sync_public_state()

    def fill(self, pending_samples: deque, capacity: int) -> None:
        """Queue one proposal for each ready chain until capacity or exhaustion.

        First score a chain's initial point, then draw movement proposals. Store
        chain ID, proposal kind and selected coordinates under each sample ID.
        """
        capacity = int(capacity)
        while len(pending_samples) < capacity and not self.is_exhausted():
            chain = self._next_available_chain()
            if chain is None:
                return

            if not chain.initialized:
                sample = self._sample_from_x(chain.xinit, chain.chain_id)
                self._pending_context[sample["sample_id"]] = {
                    "kind": "initial",
                    "chain_id": chain.chain_id,
                    "xtry": list(chain.xinit),
                }
            else:
                xpartition = self._select_partition(chain)
                xtry = self._exec_movement(xpartition, chain.xstart, chain.stp_cur)
                sample = self._sample_from_x(xtry, chain.chain_id)
                self._pending_context[sample["sample_id"]] = {
                    "kind": "move",
                    "chain_id": chain.chain_id,
                    "xpartition": xpartition,
                    "xtry": xtry,
                }
            chain.pending_sample_id = sample["sample_id"]
            pending_samples.append(sample)

    def tell(self, completed_samples: list[dict]) -> None:
        """Accept returned losses into the matching chain and refresh public best-state summaries.

        Only complete feedback for the pending sample advances that chain; other
        chains may remain outstanding. Missing/nonfinite losses are treated as infinity.
        """
        for completed in sorted(completed_samples, key=lambda item: int(item["sample_id"])):
            context = self._pending_context.pop(int(completed["sample_id"]), None)
            if context is None:
                continue
            chain = self.chains[int(context["chain_id"])]
            chain.pending_sample_id = None
            new_nrgy = self._loss_from_completed(completed)
            if context["kind"] == "initial":
                self._initialize_energy(chain, new_nrgy)
            else:
                self._accept_or_reject_movement(chain, context, new_nrgy)
            self._sync_public_state()

    def is_exhausted(self) -> bool:
        """Return whether every annealing chain has reached its stop condition.
        """
        return all(chain.exhausted for chain in self.chains)

    def estimated_sample_count(self) -> int | None:
        """Return the maximum candidate count; temperature stopping may finish sooner.
        """
        return self.chain_count * self.max_iters

    def pending_batch_limit(self, batch_size: int) -> int | None:
        """Return enough batches to hold one outstanding candidate per annealing chain.
        """
        return max(1, math.ceil(self.chain_count / max(1, int(batch_size))))

    def _new_chain(self, chain_id: int, xinit: list[float]) -> _SimAnnChain:
        """Create unscored initial/best/current positions and empty movement counters.
        """
        return _SimAnnChain(
            chain_id=chain_id,
            xinit=list(xinit),
            xopt=list(xinit),
            xstart=list(xinit),
            stp_cur=list(self.stp_ini),
            improved=[0.0] * self.vars,
            attempted=[0.0] * self.vars,
            improvement_ratio=[1.0] * self.vars,
            temp=self.initial_temp,
            exhausted=self.vars == 0 or self.max_iters < 1,
        )

    def _latin_hypercube_points(self, count: int) -> list[list[float]]:
        """Draw one start per chain with one sample per equal interval in each physical range.

        Each coordinate uses shuffled stratified fractions. This standard-library
        RNG sequence is retained independently of the NumPy helper used by Adam/LHS.
        """
        if self.vars == 0:
            return [[] for _ in range(count)]

        points = [[0.0] * self.vars for _ in range(count)]
        for dim in range(self.vars):
            span = self.xmax[dim] - self.xmin[dim]
            fractions = [
                (idx + self.rng.random()) / count
                for idx in range(count)
            ]
            self.rng.shuffle(fractions)
            for chain_id, fraction in enumerate(fractions):
                points[chain_id][dim] = self.xmin[dim] + span * fraction
        return points

    def _next_available_chain(self) -> _SimAnnChain | None:
        """Select a nonexhausted chain without outstanding work, rotating fairly through chains.
        """
        if not self.chains:
            return None
        for offset in range(len(self.chains)):
            idx = (self._next_chain_index + offset) % len(self.chains)
            chain = self.chains[idx]
            if not chain.exhausted and chain.pending_sample_id is None:
                self._next_chain_index = (idx + 1) % len(self.chains)
                return chain
        return None

    def _sample_from_x(self, x_values: Iterable[float], chain_id: int) -> dict:
        """Expand a chain's physical logical coordinates into a full linked-parameter sample.
        """
        sample = self._sample_from_selected_params(
            {
                name: float(value)
                for name, value in zip(self.x_names, x_values)
            }
        )
        sample["chain_id"] = int(chain_id)
        return sample

    def _initialize_energy(self, chain: _SimAnnChain, einit: float) -> None:
        """Record the first scored loss as initial/current/best and enforce the evaluation cap.
        """
        chain.einit = einit
        chain.nrgy_opt = einit
        chain.old_nrgy = einit
        chain.min_nrgy = einit
        chain.initialized = True
        if chain.iter >= self.max_iters:
            chain.exhausted = True

    def _select_partition(self, chain: _SimAnnChain) -> list[bool]:
        """Choose movement coordinates using improvement ratios, forcing unattempted ones.
        """
        xpartition = []
        for k in range(self.vars):
            rand = self.rng.random()
            threshold = self.use_min + (self.use_max - self.use_min) * chain.improvement_ratio[k]
            xpartition.append(rand < threshold)

        xpartition = [
            True if abs(attempted) < 1.0e-15 else selected
            for selected, attempted in zip(xpartition, chain.attempted)
        ]
        if not any(xpartition):
            xpartition = [True] * self.vars
        return xpartition

    def _exec_movement(self, xpartition: list[bool], x: list[float], step: list[float]) -> list[float]:
        """Move selected coordinates by random signed fractions of their steps and clip to bounds.
        """
        xtry = list(x)
        for k, selected in enumerate(xpartition):
            if selected:
                xrand = [self.rng.random(), self.rng.random(), self.rng.random()]
                if xrand[0] >= 0.5:
                    xtmp = xtry[k] + xrand[1] * step[k]
                else:
                    xtmp = xtry[k] - xrand[2] * step[k]

                if xtmp > self.xmax[k]:
                    xtmp = self.xmax[k]
                elif xtmp < self.xmin[k]:
                    xtmp = self.xmin[k]
                xtry[k] = xtmp
        return xtry

    def _accept_or_reject_movement(self, chain: _SimAnnChain, context: dict, new_nrgy: float) -> None:
        """Process a scored move, then cool/adapt a completed stage if needed.

        Downhill moves always replace the current point; uphill moves use the
        Boltzmann acceptance probability. Track per-coordinate attempts/improvements,
        update best accepted state, enforce max_iters, then apply the stage equations
        described above. The next stage restarts from the best accepted point.
        """
        if chain.old_nrgy is None or chain.nrgy_opt is None or chain.min_nrgy is None:
            raise RuntimeError("Enhanced simulated annealing received a movement before initial energy")

        xpartition = list(context["xpartition"])
        xtry = list(context["xtry"])
        for k, selected in enumerate(xpartition):
            if selected:
                chain.attempted[k] += 1.0

        delta_nrgy = new_nrgy - chain.old_nrgy
        chain.tot_nrgy += new_nrgy
        chain.iter += 1

        if delta_nrgy <= 0.0:
            for k, selected in enumerate(xpartition):
                if selected:
                    chain.improved[k] += 1.0
            chain.old_nrgy = new_nrgy
            chain.xstart = list(xtry)
            if new_nrgy < chain.min_nrgy:
                chain.min_nrgy = new_nrgy
            if new_nrgy < chain.nrgy_opt:
                if chain.nrgy_opt - new_nrgy > self.f_tol * chain.nrgy_opt:
                    chain.stages_w_no_improve = 0
                chain.xopt = list(xtry)
                chain.nrgy_opt = new_nrgy
        else:
            accept_probability = 0.0 if chain.temp <= 0.0 else math.exp(-delta_nrgy / chain.temp)
            if self.rng.random() <= accept_probability:
                chain.old_nrgy = new_nrgy
                chain.xstart = list(xtry)

        if chain.iter >= self.max_iters:
            chain.exhausted = True
            return

        if sum(chain.improved) < self.n1 * self.vars and sum(chain.attempted) < self.n2 * self.vars:
            return

        if chain.tot_nrgy > 0.0 and math.isfinite(chain.tot_nrgy):
            decay = chain.min_nrgy * sum(chain.attempted) / chain.tot_nrgy
            if not math.isfinite(decay):
                decay = self.min_cool
            chain.temp = chain.temp * max(min(decay, self.min_cool), self.max_cool)
        else:
            chain.temp = chain.temp * self.min_cool

        chain.improvement_ratio = [
            improved / attempted if attempted != 0.0 else 0.0
            for improved, attempted in zip(chain.improved, chain.attempted)
        ]

        chain.stp_cur = [
            step * (self.stp_adjst_factor * (ratio - self.stp_adjst_shift) + 1.0)
            for step, ratio in zip(chain.stp_cur, chain.improvement_ratio)
        ]

        if chain.temp <= self.max_final_temp and chain.stages_w_no_improve > self.no_improve_max:
            chain.exhausted = True
            return

        chain.stages_w_no_improve += 1
        chain.improved = [0.0] * self.vars
        chain.attempted = [0.0] * self.vars
        chain.tot_nrgy = 0.0
        chain.xstart = list(chain.xopt)
        chain.min_nrgy = chain.nrgy_opt

    def _sync_public_state(self) -> None:
        """Expose the best initialized chain plus aggregate iteration/temperature summaries.
        """
        initialized = [chain for chain in self.chains if chain.nrgy_opt is not None]
        if initialized:
            best_chain = min(initialized, key=lambda chain: float(chain.nrgy_opt))
            self.xopt = list(best_chain.xopt)
            self.nrgy_opt = float(best_chain.nrgy_opt)
        self.iter = max((chain.iter for chain in self.chains), default=1)
        self.temp = min((chain.temp for chain in self.chains), default=self.initial_temp)

    @staticmethod
    def _loss_from_completed(completed: dict) -> float:
        """Read total_loss as a finite scalar, treating unusable feedback as infinity.
        """
        try:
            value = float(completed["total_loss"])
        except (KeyError, TypeError, ValueError):
            return math.inf
        return value if math.isfinite(value) else math.inf

    @staticmethod
    def prepare_options(options, request):
        """Default the chain count to worker*batch capacity, preserving explicit SimAnn counts.
        """
        options.setdefault("chain_count", max(1, request["max_workers"] * request["batch_size"]))

    @staticmethod
    def sample_count(options, _ranges):
        """Return the maximum chain_count*max_iters candidate budget.
        """
        return options["max_iters"] * options["chain_count"]
