"""Parallel multi-start Adam using black-box SPSA gradients."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
import math
from typing import Iterable

import numpy as np

from tuner.tuning_strategy import BaseTuningStrategy


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
        return all(chain.phase == "done" for chain in self.chains)

    def estimated_sample_count(self) -> int | None:
        if not self.active.any():
            return 1
        return self.chain_count * (2 + 2 * self.spsa_pairs * self.max_updates)

    def pending_batch_limit(self, _batch_size: int) -> int | None:
        return self.concurrent_batches

    def _starting_points(self, seed_sequence: np.random.SeedSequence) -> list[np.ndarray]:
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
            strata = np.column_stack([
                rng.permutation(sample_count)
                for _ in range(active_count)
            ])
            samples = (strata + rng.random((sample_count, active_count))) / sample_count
            for sample in samples:
                point = np.zeros(self.dimension, dtype=float)
                point[self.active] = sample
                starts.append(point)
        return starts

    def _append_center(self, pending_samples: deque, chain: _AdamChain, kind: str) -> None:
        sample = self._sample_from_x(chain.x)
        sample_id = int(sample["sample_id"])
        self._pending_context[sample_id] = {"kind": kind, "chain_id": chain.chain_id}
        chain.pending_sample_ids = {sample_id}
        chain.phase = f"{kind}_waiting"
        pending_samples.append(sample)

    def _append_probes(self, pending_samples: deque, chain: _AdamChain) -> None:
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
        physical = self.minimum + self.span * point
        return self._sample_from_selected_params(
            {
                spec["name"]: float(physical[index])
                for index, spec in enumerate(self.parameter_ranges)
            }
        )

    @staticmethod
    def _loss(completed: dict) -> float:
        try:
            loss = float(completed["total_loss"])
        except (KeyError, TypeError, ValueError):
            return math.inf
        return loss if math.isfinite(loss) else math.inf
