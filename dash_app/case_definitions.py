import ast
import glob
import os

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CASE_DEFINITIONS_PATH = os.path.join(
    REPO_ROOT, "postprocessing", "pyplotgen", "config", "Case_definitions.py"
)

LEGACY_BENCHMARK_ROOT = "/home/pub/les_and_clubb_benchmark_runs/"
UNKNOWN = object()


def _find_first_match(root, subdir, filename):
    pattern = os.path.join(root, subdir, "**", filename)
    matches = glob.glob(pattern, recursive=True)
    return matches[0] if matches else None


def _normalize_path(path, repo_root):
    if not path:
        return None
    path = os.path.expanduser(path)
    if os.path.exists(path):
        return path

    if path.startswith(LEGACY_BENCHMARK_ROOT):
        if os.path.exists(path):
            return path
        return None

    if path.startswith(os.sep):
        candidate = os.path.join(repo_root, path.lstrip(os.sep))
        if os.path.exists(candidate):
            return candidate

    basename = os.path.basename(path)
    for subdir in ["output", "les_runs", "sam_benchmark_runs"]:
        candidate = os.path.join(repo_root, subdir, basename)
        if os.path.exists(candidate):
            return candidate

    for subdir in ["sam_benchmark_runs", "les_runs"]:
        candidate = _find_first_match(repo_root, subdir, basename)
        if candidate:
            return candidate

    return None


def _dict_pick(entry, keys):
    if not isinstance(entry, dict):
        return []
    return [entry.get(key) for key in keys if entry.get(key)]


def _choose_file(repo_root, candidates):
    for candidate in candidates:
        normalized = _normalize_path(candidate, repo_root)
        if normalized and os.path.exists(normalized):
            return normalized
    return None


def _choose_files(repo_root, candidates):
    results = []
    seen = set()
    for candidate in candidates:
        normalized = _normalize_path(candidate, repo_root)
        if normalized and os.path.exists(normalized) and normalized not in seen:
            results.append(normalized)
            seen.add(normalized)
    return results


def _eval_node(node, env, case_defs):
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.Str):
        return node.s
    if isinstance(node, ast.Num):
        return node.n
    if isinstance(node, ast.Name):
        if node.id in env:
            return env[node.id]
        if node.id in case_defs:
            return case_defs[node.id]
        return UNKNOWN
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        val = _eval_node(node.operand, env, case_defs)
        if isinstance(val, (int, float)):
            return -val
        return UNKNOWN
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
            if key is UNKNOWN:
                continue
            val = _eval_node(val_node, env, case_defs)
            if val is UNKNOWN:
                continue
            result[key] = val
        return result
    if isinstance(node, (ast.List, ast.Tuple)):
        items = []
        for elt in node.elts:
            val = _eval_node(elt, env, case_defs)
            if val is UNKNOWN:
                continue
            items.append(val)
        return items
    return UNKNOWN


def _extract_case_names(node):
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
    with open(case_path, "r", encoding="utf-8") as handle:
        tree = ast.parse(handle.read(), filename=case_path)

    env = {}
    case_defs = {}
    case_order = []
    all_cases_node = None
    cases_to_plot_node = None

    for node in tree.body:
        if isinstance(node, ast.Assign):
            if not node.targets:
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

    cases = []
    for name in case_names:
        case = case_defs.get(name)
        if case:
            cases.append(case)
    return cases


def load_case_definitions(repo_root=None):
    if repo_root is None:
        repo_root = REPO_ROOT

    raw_cases = _parse_cases(CASE_DEFINITIONS_PATH)
    cases = []
    for case in raw_cases:
        name = case.get("name", "")

        les_candidates = []
        les_candidates += _dict_pick(case.get("coamps_benchmark_file"), ["sm", "sw"])
        les_candidates += _dict_pick(case.get("sam_benchmark_file"), ["sam_benchmark"])
        les_candidates += _dict_pick(case.get("wrf_benchmark_file"), ["lasso_benchmark"])
        les_candidates += _dict_pick(case.get("sam_file"), ["sam"])
        les_candidates += _dict_pick(case.get("cam_file"), ["cam"])
        les_candidates += _dict_pick(case.get("e3sm_file"), ["e3sm"])

        clubb_candidates = []
        clubb_candidates += _dict_pick(case.get("clubb_file"), ["zt", "zm", "sfc", "subcolumns"])
        clubb_candidates += _dict_pick(case.get("wrf_file"), ["zt", "zm", "sfc", "subcolumns"])

        les_file = _choose_file(repo_root, les_candidates)
        clubb_files = _choose_files(repo_root, clubb_candidates)
        clubb_file = clubb_files[0] if clubb_files else None
        available = bool(les_file and clubb_files)

        cases.append(
            {
                "name": name,
                "start_time": case.get("start_time"),
                "end_time": case.get("end_time"),
                "les_file": les_file,
                "clubb_file": clubb_file,
                "clubb_files": clubb_files,
                "available": available,
            }
        )

    return cases
