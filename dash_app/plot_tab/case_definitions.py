"""Load pyplotgen case defaults needed by the Dash plots tab."""

import ast
import os


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CASE_DEFINITIONS_PATH = os.path.join(
    REPO_ROOT, "postprocessing", "pyplotgen", "config", "Case_definitions.py"
)
UNKNOWN = object()


def _legacy_literal_value(node):
    """Support Python versions where string/number literals predate ast.Constant."""
    legacy_str = getattr(ast, "Str", None)
    if legacy_str is not None and isinstance(node, legacy_str):
        return node.s

    legacy_num = getattr(ast, "Num", None)
    if legacy_num is not None and isinstance(node, legacy_num):
        return node.n

    return UNKNOWN


def _eval_node(node, env, case_defs):
    """Evaluate the small AST subset used by pyplotgen case definitions."""
    if isinstance(node, ast.Constant):
        return node.value
    legacy_value = _legacy_literal_value(node)
    if legacy_value is not UNKNOWN:
        return legacy_value
    if isinstance(node, ast.Name):
        if node.id in env:
            return env[node.id]
        if node.id in case_defs:
            return case_defs[node.id]
        return node.id
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        value = _eval_node(node.operand, env, case_defs)
        return -value if isinstance(value, (int, float)) else UNKNOWN
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        left = _eval_node(node.left, env, case_defs)
        right = _eval_node(node.right, env, case_defs)
        if left is UNKNOWN or right is UNKNOWN:
            return UNKNOWN
        if isinstance(left, str) or isinstance(right, str):
            return f"{left}{right}"
        if isinstance(left, (int, float)) and isinstance(right, (int, float)):
            return left + right
        return UNKNOWN
    if isinstance(node, ast.Dict):
        result = {}
        for key_node, val_node in zip(node.keys, node.values):
            key = _eval_node(key_node, env, case_defs)
            val = _eval_node(val_node, env, case_defs)
            if key is UNKNOWN or val is UNKNOWN:
                continue
            result[key] = val
        return result
    if isinstance(node, (ast.List, ast.Tuple)):
        items = []
        for elt in node.elts:
            value = _eval_node(elt, env, case_defs)
            if value is not UNKNOWN:
                items.append(value)
        return items
    return UNKNOWN


def _extract_case_names(node):
    """Extract a flat list of case variable names from ALL_CASES-style assignments."""
    if isinstance(node, ast.Name):
        return [node.id]
    if isinstance(node, (ast.List, ast.Tuple)):
        names = []
        for elt in node.elts:
            if isinstance(elt, ast.Name):
                names.append(elt.id)
        return names
    return []


def _parse_cases(case_path):
    """Parse pyplotgen Case_definitions.py and return the ordered raw case dicts."""
    with open(case_path, "r", encoding="utf-8") as handle:
        tree = ast.parse(handle.read(), filename=case_path)

    env = {}
    case_defs = {}
    case_order = []
    all_cases_node = None
    cases_to_plot_node = None

    for node in tree.body:
        if not isinstance(node, ast.Assign) or not node.targets:
            continue
        target = node.targets[0]
        if not isinstance(target, ast.Name):
            continue
        name = target.id
        if name == "ALL_CASES":
            all_cases_node = node.value
            continue
        if name == "CASES_TO_PLOT":
            cases_to_plot_node = node.value
            continue
        value = _eval_node(node.value, env, case_defs)
        if value is UNKNOWN:
            continue
        env[name] = value
        if isinstance(value, dict) and "name" in value:
            case_defs[name] = value
            case_order.append(name)

    case_names = []
    if cases_to_plot_node is not None:
        case_names = _extract_case_names(cases_to_plot_node)
        if case_names == ["ALL_CASES"] and all_cases_node is not None:
            case_names = _extract_case_names(all_cases_node)
    elif all_cases_node is not None:
        case_names = _extract_case_names(all_cases_node)
    if not case_names:
        case_names = case_order

    return [case_defs[name] for name in case_names if name in case_defs]


def load_case_definitions():
    """Return the pyplotgen case defaults needed by the Dash plots tab."""
    raw_cases = _parse_cases(CASE_DEFINITIONS_PATH)
    return [
        {
            "name": case.get("name", ""),
            "start_time": case.get("start_time"),
            "end_time": case.get("end_time"),
            "height_min_value": case.get("height_min_value"),
            "height_max_value": case.get("height_max_value"),
            "sam_benchmark_file": case.get("sam_benchmark_file"),
            "coamps_benchmark_file": case.get("coamps_benchmark_file"),
            "wrf_benchmark_file": case.get("wrf_benchmark_file"),
            "var_groups": list(case.get("var_groups") or []),
        }
        for case in raw_cases
    ]
