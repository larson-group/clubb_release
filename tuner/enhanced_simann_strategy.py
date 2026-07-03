"""Enhanced simulated annealing strategy for the CLUBB tuner."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import math
import random
from typing import Iterable

from tuner.tuning_strategy import BaseTuningStrategy


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
        max_iters: int = 2000,
        chain_count: int = 1,
        stp_adjst_shift: float = 0.5,
        stp_adjst_factor: float = 1.0,
        f_tol: float = 1.0e-4,
        seed: int | None = None,
    ):
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
        return all(chain.exhausted for chain in self.chains)

    def estimated_sample_count(self) -> int | None:
        return self.chain_count * self.max_iters

    def _new_chain(self, chain_id: int, xinit: list[float]) -> _SimAnnChain:
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
        sample = self._sample_from_selected_params(
            {
                name: float(value)
                for name, value in zip(self.x_names, x_values)
            }
        )
        sample["chain_id"] = int(chain_id)
        return sample

    def _initialize_energy(self, chain: _SimAnnChain, einit: float) -> None:
        chain.einit = einit
        chain.nrgy_opt = einit
        chain.old_nrgy = einit
        chain.min_nrgy = einit
        chain.initialized = True
        if chain.iter >= self.max_iters:
            chain.exhausted = True

    def _select_partition(self, chain: _SimAnnChain) -> list[bool]:
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
        initialized = [chain for chain in self.chains if chain.nrgy_opt is not None]
        if initialized:
            best_chain = min(initialized, key=lambda chain: float(chain.nrgy_opt))
            self.xopt = list(best_chain.xopt)
            self.nrgy_opt = float(best_chain.nrgy_opt)
        self.iter = max((chain.iter for chain in self.chains), default=1)
        self.temp = min((chain.temp for chain in self.chains), default=self.initial_temp)

    @staticmethod
    def _loss_from_completed(completed: dict) -> float:
        try:
            value = float(completed["total_loss"])
        except (KeyError, TypeError, ValueError):
            return math.inf
        return value if math.isfinite(value) else math.inf
