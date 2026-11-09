"""Normalize strategy requests and construct their registered implementations.

normalize_strategy_config validates a strategy object and its options through
the JSON catalog. build_tuning_strategy constructs the selected class with the
canonical parameter ranges, defaults and seed. Strategy
implementations and their shared interface live in tuner.algorithms.
"""

from __future__ import annotations
from typing import Iterable

from tuner.algorithm_registry import ALGORITHMS, get_algorithm, normalize_options


def normalize_strategy_config(request: dict) -> dict:
    """Validate the strategy object and return its catalog-normalized options."""
    raw_strategy = request.get("strategy")
    if not isinstance(raw_strategy, dict):
        raise ValueError("strategy must be an object")

    name = str(raw_strategy.get("name", "")).strip().lower()
    if not name:
        raise ValueError("strategy.name is required")
    if name not in ALGORITHMS:
        raise ValueError(f"Unknown tuning strategy: {name}")

    raw_options = raw_strategy.get("options", {})
    if not isinstance(raw_options, dict):
        raise ValueError("strategy.options must be an object")

    return {"name": name, "options": normalize_options(name, raw_options)}



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
    spec = get_algorithm(strategy_config["name"])
    return spec.strategy_class(
        param_names=param_names,
        default_params_row=default_params_row,
        parameter_ranges=parameter_ranges,
        **strategy_config["options"],
        seed=seed,
    )
