"""Parameter proposal strategies for the CLUBB tuner."""

from __future__ import annotations

from collections import deque
import itertools
import math
import random
from typing import Iterable


VALID_STRATEGY_NAMES = {"random", "resolve"}


def _float_from_raw(value, label: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be numeric") from exc


def _int_from_raw(value, label: str) -> int:
    try:
        numeric_value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be an integer") from exc
    if int(numeric_value) != numeric_value:
        raise ValueError(f"{label} must be an integer")
    return int(numeric_value)


def normalize_strategy_config(request: dict) -> dict:
    """Return a canonical strategy config while accepting legacy max_samples."""
    raw_strategy = request.get("strategy")
    if raw_strategy is None:
        raw_strategy = {
            "name": "random",
            "options": {},
        }
        if request.get("max_samples") is not None:
            raw_strategy["options"]["max_samples"] = request.get("max_samples")
    elif isinstance(raw_strategy, str):
        raw_strategy = {"name": raw_strategy, "options": {}}
    elif not isinstance(raw_strategy, dict):
        raise ValueError("strategy must be an object")

    name = str(raw_strategy.get("name", "")).strip().lower()
    if not name:
        raise ValueError("strategy.name is required")
    if name not in VALID_STRATEGY_NAMES:
        raise ValueError(f"Unknown tuning strategy: {name}")

    raw_options = raw_strategy.get("options", {})
    if raw_options is None:
        raw_options = {}
    if not isinstance(raw_options, dict):
        raise ValueError("strategy.options must be an object")

    options = {}
    if name == "random":
        unknown = set(raw_options) - {"max_samples"}
        if unknown:
            raise ValueError("Unknown random strategy option(s): " + ", ".join(sorted(unknown)))
        if raw_options.get("max_samples") is not None:
            max_samples = _int_from_raw(raw_options.get("max_samples"), "random max_samples")
            if max_samples < 1:
                raise ValueError("random max_samples must be >= 1")
            options["max_samples"] = max_samples

    if name == "resolve":
        unknown = set(raw_options) - {"spacing"}
        if unknown:
            raise ValueError("Unknown resolve strategy option(s): " + ", ".join(sorted(unknown)))
        if raw_options.get("spacing") is None:
            raise ValueError("resolve spacing is required")
        spacing = _float_from_raw(raw_options.get("spacing"), "resolve spacing")
        if spacing <= 0.0:
            raise ValueError("resolve spacing must be > 0")
        options["spacing"] = spacing

    return {"name": name, "options": options}


class BaseTuningStrategy:
    """Common helpers for sample-generating tuning strategies."""

    def __init__(
        self,
        *,
        param_names: Iterable[str],
        default_params_row: Iterable[float],
        parameter_ranges: Iterable[dict],
    ):
        self.param_names = list(param_names)
        self.default_params_row = [float(value) for value in default_params_row]
        self.param_index = {name: idx for idx, name in enumerate(self.param_names)}
        self.parameter_ranges = [dict(spec) for spec in parameter_ranges]
        self.next_sample_id = 0

        if len(self.default_params_row) != len(self.param_names):
            raise ValueError("default_params_row length must match param_names length")

        missing = [spec["name"] for spec in self.parameter_ranges if spec["name"] not in self.param_index]
        if missing:
            raise ValueError("Unknown tunable parameter(s): " + ", ".join(sorted(missing)))

    def _sample_from_selected_params(self, selected_params: dict[str, float]) -> dict:
        row = list(self.default_params_row)
        for name, value in selected_params.items():
            row[self.param_index[name]] = float(value)

        sample_id = self.next_sample_id
        self.next_sample_id += 1
        return {
            "sample_id": sample_id,
            "param_row": row,
            "selected_params": {name: float(value) for name, value in selected_params.items()},
            "all_params": {name: float(row[idx]) for idx, name in enumerate(self.param_names)},
        }

    def fill(self, pending_samples: deque, capacity: int) -> None:
        """Append sample dictionaries until pending_samples reaches capacity."""
        raise NotImplementedError

    def tell(self, _completed_samples: list[dict]) -> None:
        """Accept completed samples; v1 strategies are non-adaptive."""
        return None

    def is_exhausted(self) -> bool:
        """Return whether the strategy has no more samples to propose."""
        return False

    def estimated_sample_count(self) -> int | None:
        """Return a finite sample count estimate when the strategy has one."""
        return None


class RandomUniformStrategy(BaseTuningStrategy):
    """Uniform random sampler over selected tunable-parameter ranges."""

    def __init__(
        self,
        *,
        param_names: Iterable[str],
        default_params_row: Iterable[float],
        parameter_ranges: Iterable[dict],
        max_samples: int | None = None,
        seed: int | None = None,
    ):
        super().__init__(
            param_names=param_names,
            default_params_row=default_params_row,
            parameter_ranges=parameter_ranges,
        )
        self.max_samples = None if max_samples is None else int(max_samples)
        self.rng = random.Random(seed)

    def fill(self, pending_samples: deque, capacity: int) -> None:
        capacity = int(capacity)
        while len(pending_samples) < capacity and not self.is_exhausted():
            selected_params = {}
            for spec in self.parameter_ranges:
                name = spec["name"]
                selected_params[name] = self.rng.uniform(float(spec["min"]), float(spec["max"]))
            pending_samples.append(self._sample_from_selected_params(selected_params))

    def is_exhausted(self) -> bool:
        return self.max_samples is not None and self.next_sample_id >= self.max_samples

    def estimated_sample_count(self) -> int | None:
        return self.max_samples


class ResolveGridStrategy(BaseTuningStrategy):
    """Deterministic full-space sampler using evenly resolved parameter ranges."""

    def __init__(
        self,
        *,
        param_names: Iterable[str],
        default_params_row: Iterable[float],
        parameter_ranges: Iterable[dict],
        spacing: float,
    ):
        super().__init__(
            param_names=param_names,
            default_params_row=default_params_row,
            parameter_ranges=parameter_ranges,
        )
        self.spacing = float(spacing)
        if self.spacing <= 0.0:
            raise ValueError("spacing must be > 0")
        self.grid_names = [spec["name"] for spec in self.parameter_ranges]
        self.grid_values = [
            self._values_for_range(float(spec["min"]), float(spec["max"]), self.spacing)
            for spec in self.parameter_ranges
        ]
        self.total_samples = math.prod(len(values) for values in self.grid_values)
        self._grid_iter = itertools.product(*self.grid_values)
        self._exhausted = self.total_samples == 0

    @staticmethod
    def _values_for_range(min_value: float, max_value: float, spacing: float) -> list[float]:
        span = float(max_value) - float(min_value)
        if span < 0.0:
            raise ValueError("resolve ranges require min <= max")
        if span == 0.0:
            return [float(min_value)]

        intervals = max(1, int(math.ceil(span / spacing)))
        values = [
            float(min_value + span * idx / intervals)
            for idx in range(intervals)
        ]
        values.append(float(max_value))
        return values

    def fill(self, pending_samples: deque, capacity: int) -> None:
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
        return self._exhausted

    def estimated_sample_count(self) -> int | None:
        return self.total_samples


def build_tuning_strategy(
    strategy_config: dict,
    *,
    param_names: Iterable[str],
    default_params_row: Iterable[float],
    parameter_ranges: Iterable[dict],
    seed: int | None = None,
):
    """Instantiate the requested sample-generation strategy."""
    strategy_config = normalize_strategy_config({"strategy": strategy_config})
    name = strategy_config["name"]
    options = strategy_config.get("options", {})

    if name == "random":
        return RandomUniformStrategy(
            param_names=param_names,
            default_params_row=default_params_row,
            parameter_ranges=parameter_ranges,
            max_samples=options.get("max_samples"),
            seed=seed,
        )

    if name == "resolve":
        return ResolveGridStrategy(
            param_names=param_names,
            default_params_row=default_params_row,
            parameter_ranges=parameter_ranges,
            spacing=options["spacing"],
        )

    raise ValueError(f"Unknown tuning strategy: {name}")
