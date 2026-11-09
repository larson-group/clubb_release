"""Render the tuner's registered algorithms and their typed option controls.

build_algorithm_controls builds every panel in registry order, preserving edited
values when switching methods. option_ids/options_from_values/option_values map
Dash's wildcard values to canonical strategy options and back for saved requests.
algorithm_summary previews requirements and evaluation budgets using the same
registry that validates CLI and worker requests. Scientific behavior stays in
tuner; this module only adapts its declarations to Dash components.
"""

import json

from dash import dcc, html

from tuner.algorithm_registry import ALGORITHMS, apply_requirements, get_algorithm, prepare_strategy


def option_ids():
    """Return option component IDs in the layout's registry/declaration order."""
    return [{"type": "tune-algorithm-option", "algorithm": name, "name": option.name}
            for name in ALGORITHMS for option in get_algorithm(name).options]


def option_values(name=None, options=None):
    """Hydrate one saved strategy while supplying defaults for other panels."""
    result = []
    for algorithm_name in ALGORITHMS:
        for option in get_algorithm(algorithm_name).options:
            default = option.ui_default if option.ui_default is not None else option.default
            result.append((options or {}).get(option.name, default) if name == algorithm_name else default)
    return result


def options_from_values(name, values):
    """Read the selected panel without leaking options from other algorithms."""
    return {identifier["name"]: value for identifier, value in zip(option_ids(), values or [])
            if identifier["algorithm"] == name and value is not None}


def build_algorithm_controls(selected, options=None):
    """Build typed controls from declarations; no algorithm-specific IDs are needed."""
    values = iter(option_values(selected, options))
    panels = []
    for name in ALGORITHMS:
        spec = get_algorithm(name)
        controls = []
        for option in spec.options:
            identifier = {"type": "tune-algorithm-option", "algorithm": name, "name": option.name}
            value = next(values)
            if option.choices or option.kind == "bool":
                choices = option.choices or (True, False)
                control = dcc.Dropdown(id=identifier, options=[{"label": str(item), "value": item} for item in choices],
                                       value=value, clearable=not option.required,
                                       className="tune-algorithm-choice")
            else:
                kwargs = {"type": "number" if option.kind in {"int", "float"} else "text"}
                if option.kind in {"int", "float"}:
                    kwargs.update(min=option.minimum, max=option.maximum,
                                  step=1 if option.kind == "int" else "any")
                control = dcc.Input(id=identifier, value=value, debounce=True,
                                    className="tune-algorithm-input", **kwargs)
            controls.append(html.Div([html.Label(option.label, htmlFor=json.dumps(identifier, sort_keys=True, separators=(",", ":"))), control],
                                     title=option.description, className="tune-algorithm-option"))
        panels.append(html.Div([
            html.Div(spec.description, className="tune-algorithm-description"),
            html.Div(controls, className="tune-algorithm-fields"),
        ], id={"type": "tune-algorithm-panel", "name": name},
                               style={"display": "block" if selected == name else "none", "marginTop": "10px"}))
    return panels


def algorithm_backend_error(name, backend):
    """Explain an incompatible runtime using the registry's requirement validation."""
    if name not in ALGORITHMS:
        return ""
    try:
        apply_requirements({"strategy": {"name": name}, "backend": backend})
    except ValueError:
        return f"{get_algorithm(name).label} requires JAX gradients. Select JAX in the runner chooser."
    return ""


def algorithm_summary(name, values, batch_size, max_workers, cases, ranges, *, backend):
    """Preview the registered budget and requirements, reporting invalid options."""
    backend_error = algorithm_backend_error(name, backend)
    if backend_error:
        return backend_error
    try:
        spec = get_algorithm(name)
        request = {"strategy": {"name": name, "options": options_from_values(name, values)},
                   "backend": backend, "batch_size": int(batch_size), "max_workers": int(max_workers),
                   "cases": list(cases or []), "parameter_ranges": ranges}
        if request["batch_size"] < 1 or request["max_workers"] < 1 or not request["cases"]:
            return "Enter valid scheduling and case settings."
        prepare_strategy(request)
        count = request["total_samples"]
        budget = f"{count} candidate evaluations" if count is not None else "No fixed evaluation limit"
        return budget + (" · requires JAX gradients" if spec.requires_gradients else "")
    except (ValueError, TypeError, KeyError) as exc:
        return str(exc)


def options_ready(name, values, batch_size, max_workers, cases, *, backend):
    """Validate options and execution layout without evaluating any model."""
    try:
        request = {"strategy": {"name": name, "options": options_from_values(name, values)},
                   "backend": backend, "batch_size": int(batch_size), "max_workers": int(max_workers),
                   "cases": list(cases or []), "parameter_ranges": []}
        if min(request["batch_size"], request["max_workers"]) < 1 or not request["cases"]:
            return False
        prepare_strategy(request)
        return True
    except (ValueError, TypeError, KeyError):
        return False
