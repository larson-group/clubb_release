"""Describe and construct the tuning algorithms available to every frontend.

Add a class under tuner/algorithms and an entry in algorithms/algorithms.json.
Option and Algorithm are parsed views of that catalog, not declarations that
implementations must repeat. get_algorithm/list_algorithms read metadata without
loading implementations; normalize_options validates user inputs, and
prepare_strategy asks the selected class for execution-dependent settings.
parse_strategy accepts either positional colon options or named colon options,
for example adam_gradient:max_updates=10:learning_rate=0.01.

The catalog and validation use only the standard library; implementations load
when a strategy is prepared or constructed.
Strategies still communicate with the scheduler through fill/tell and must keep
pickleable state for graceful-stop continuation. No model code runs here.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
from importlib import import_module
import math
import json
from pathlib import Path


CATALOG_PATH = Path(__file__).with_name("algorithms") / "algorithms.json"


def _read_catalog(path: Path) -> dict:
    """Read the algorithm-name-to-metadata mapping from a JSON catalog."""
    catalog = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(catalog, dict):
        raise ValueError("Algorithm catalog must be a JSON object")
    return catalog


ALGORITHMS = _read_catalog(CATALOG_PATH)


@dataclass(frozen=True)
class Option:
    """Represent and validate one catalog option for CLI and frontend use."""

    name: str
    kind: str
    label: str
    description: str
    default: object = None
    required: bool = False
    minimum: float | None = None
    maximum: float | None = None
    exclusive_minimum: bool = False
    exclusive_maximum: bool = False
    choices: tuple = ()
    ui_default: object = None

    def parse(self, raw, algorithm: str):
        """Convert one input to its declared type and reject invalid values."""
        label = f"{algorithm} {self.name}"
        if self.kind in {"int", "float"}:
            try:
                value = float(raw)
            except (ValueError, TypeError) as exc:
                raise ValueError(f"{label} must be numeric") from exc
            if isinstance(raw, bool) or not math.isfinite(value):
                raise ValueError(f"{label} must be finite and numeric")
            if self.kind == "int":
                if int(value) != value:
                    raise ValueError(f"{label} must be an integer")
                value = int(value)
            if self.minimum is not None and (
                value < self.minimum or (self.exclusive_minimum and value == self.minimum)
            ):
                operator = ">" if self.exclusive_minimum else ">="
                raise ValueError(f"{label} must be {operator} {self.minimum:g}")
            if self.maximum is not None and (value > self.maximum or (self.exclusive_maximum and value == self.maximum)):
                operator = "<" if self.exclusive_maximum else "<="
                raise ValueError(f"{label} must be {operator} {self.maximum:g}")
        elif self.kind == "bool":
            if isinstance(raw, bool):
                value = raw
            elif isinstance(raw, str) and raw.lower() in {"true", "false"}:
                value = raw.lower() == "true"
            else:
                raise ValueError(f"{label} must be true or false")
        elif self.kind == "string":
            if not isinstance(raw, str):
                raise ValueError(f"{label} must be a string")
            value = raw
        else:
            raise ValueError(f"Unsupported option type: {self.kind}")
        if self.choices and value not in self.choices:
            raise ValueError(f"{label} must be one of {self.choices}")
        return value


@dataclass(frozen=True)
class Algorithm:
    """Expose parsed catalog metadata and lazily load its implementation class."""

    name: str
    label: str
    description: str
    implementation: str
    category: str
    options: tuple[Option, ...]
    cli_order: tuple[str, ...] = ()
    requires_gradients: bool = False
    internal_options: tuple[str, ...] = ()

    @property
    def strategy_class(self):
        """Import the catalog's class only when constructing or preparing work."""
        from tuner.algorithms.base import BaseTuningStrategy

        module, attribute = self.implementation.split(":")
        strategy = getattr(import_module(module), attribute)
        if not isinstance(strategy, type) or not issubclass(strategy, BaseTuningStrategy):
            raise ValueError(f"{self.implementation} must subclass BaseTuningStrategy")
        return strategy


@lru_cache(maxsize=None)
def get_algorithm(name: str) -> Algorithm:
    """Parse one catalog entry without importing strategy or model runtimes."""
    if name not in ALGORITHMS:
        raise ValueError(f"Unknown tuning strategy: {name}")
    metadata = dict(ALGORITHMS[name])
    if metadata.get("category") not in ("optimizer", "sampler"):
        raise ValueError(f"{name} category must be optimizer or sampler")
    options = []
    for option_name, raw in metadata.pop("options", {}).items():
        fields = dict(raw)
        kind = fields.pop("type")
        if kind not in {"int", "float", "string", "bool"}:
            raise ValueError(f"Unsupported option type for {name}.{option_name}: {kind}")
        fields["choices"] = tuple(fields.get("choices", ()))
        option = Option(option_name, kind,
                        fields.pop("label", option_name.replace("_", " ").capitalize()),
                        fields.pop("description", ""), **fields)
        for default in (option.default, option.ui_default):
            if default is not None:
                option.parse(default, name)
        options.append(option)
    metadata["cli_order"] = tuple(metadata.get("cli_order", ()))
    metadata["internal_options"] = tuple(metadata.get("internal_options", ()))
    metadata.setdefault("label", name.replace("_", " ").capitalize())
    metadata.setdefault("description", "")
    if not set(metadata["cli_order"]) <= {option.name for option in options}:
        raise ValueError(f"{name} cli_order contains undeclared options")
    if not isinstance(metadata.get("requires_gradients", False), bool):
        raise ValueError(f"{name} requires_gradients must be a boolean")
    return Algorithm(name=name, options=tuple(options), **metadata)


def list_algorithms() -> list[dict]:
    """Return JSON-compatible descriptions for CLI, Dash and other clients."""
    return [
        {"name": spec.name, "label": spec.label, "description": spec.description,
         "category": spec.category,
         "requires_gradients": spec.requires_gradients,
         "backends": ["jax"] if spec.requires_gradients else ["fortran", "jax"],
         "options": [asdict(option) for option in spec.options],
         "cli_order": list(spec.cli_order)}
        for spec in (get_algorithm(name) for name in ALGORITHMS)
    ]


def normalize_options(name: str, raw: dict, *, defaults: bool = True) -> dict:
    """Validate declared options, preserving optional execution-derived inputs."""
    spec = get_algorithm(name)
    if not isinstance(raw, dict):
        raise ValueError("strategy.options must be an object")
    unknown = set(raw) - {option.name for option in spec.options} - set(spec.internal_options)
    if unknown:
        raise ValueError(f"Unknown {name} strategy option(s): " + ", ".join(sorted(unknown)))
    result = {}
    for option in spec.options:
        value = raw.get(option.name)
        if value is None and defaults:
            value = option.default
        if value is None:
            if option.required:
                raise ValueError(f"{name} {option.name} is required")
            continue
        result[option.name] = option.parse(value, name)
    for key in spec.internal_options:
        if raw.get(key) is not None:
            result[key] = Option(key, "int", key, "Derived scheduling value", minimum=1).parse(raw[key], name)
    return result


def parse_strategy(text: str) -> dict:
    """Parse positional or name=value colon options."""
    name, *parts = [part.strip() for part in str(text).split(":")]
    name = name.lower()
    spec = get_algorithm(name)
    if any("=" in part for part in parts):
        raw = {}
        for part in parts:
            if "=" not in part:
                raise ValueError("Use either positional or name=value strategy options")
            key, value = part.split("=", 1)
            if key in raw:
                raise ValueError(f"Duplicate strategy option: {key}")
            raw[key] = value
    else:
        if len(parts) > len(spec.cli_order):
            raise ValueError(f"{name} accepts positional options: {', '.join(spec.cli_order)}")
        raw = {key: value for key, value in zip(spec.cli_order, parts) if value}
    if set(raw) & set(spec.internal_options):
        raise ValueError(f"{name} derived scheduling options cannot be set through -strategy")
    return {"name": name, "options": normalize_options(name, raw, defaults=False)}


def apply_requirements(request: dict) -> None:
    """Enable mandatory gradients and reject incompatible backends before launch."""
    if not isinstance(request.get("grad", False), bool):
        raise ValueError("grad must be a boolean")
    strategy = request.get("strategy")
    if strategy is None:
        name = "random"
    elif isinstance(strategy, dict):
        name = strategy.get("name", "random")
    else:
        raise ValueError("strategy must be an object")
    if get_algorithm(str(name).strip().lower()).requires_gradients:
        if str(request.get("backend", "fortran")).strip().lower() != "jax":
            raise ValueError(f"{name} requires the jax backend")
        request["grad"] = True


def prepare_strategy(request: dict) -> dict:
    """Normalize strategy options, apply requirements and derive scheduling inputs."""
    from tuner.tuning_strategy import normalize_strategy_config

    strategy = normalize_strategy_config(request)
    spec = get_algorithm(strategy["name"])
    apply_requirements(request)
    strategy_class = spec.strategy_class
    strategy_class.prepare_options(strategy["options"], request)
    request["strategy"] = strategy
    request["total_samples"] = strategy_class.sample_count(
        strategy["options"], request.get("parameter_ranges", []))
    return strategy
