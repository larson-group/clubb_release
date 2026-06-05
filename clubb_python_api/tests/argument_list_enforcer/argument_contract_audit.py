#!/usr/bin/env python3
"""Check public Python API and internal F2PY argument contracts."""

from __future__ import annotations

import argparse
import ast
from collections import Counter
import re
import sys
from dataclasses import dataclass
from pathlib import Path


def _strip_fortran_comments(line: str) -> str:
    if line.lstrip().startswith("#"):
        return ""

    in_single = False
    in_double = False
    out = []
    for ch in line:
        if ch == "'" and not in_double:
            in_single = not in_single
        elif ch == '"' and not in_single:
            in_double = not in_double
        elif ch == "!" and not in_single and not in_double:
            break
        out.append(ch)
    return "".join(out)


def merge_fortran_lines(text: str) -> str:
    """Merge continuation lines and remove comments."""
    lines = [_strip_fortran_comments(ln.rstrip()) for ln in text.splitlines()]
    merged: list[str] = []
    cur = ""
    for raw in lines:
        line = raw.rstrip()
        if not line.strip():
            if cur:
                merged.append(cur)
                cur = ""
            continue
        if not cur:
            cur = line
            continue

        if cur.rstrip().endswith("&") or line.lstrip().startswith("&"):
            left = cur.rstrip()
            if left.endswith("&"):
                left = left[:-1].rstrip()
            right = line.lstrip()
            if right.startswith("&"):
                right = right[1:].lstrip()
            cur = f"{left} {right}".strip()
        else:
            merged.append(cur)
            cur = line

    if cur:
        merged.append(cur)

    return "\n".join(merged)


def split_args(arg_text: str) -> list[str]:
    """Split comma-separated arguments while respecting nested parentheses."""
    args: list[str] = []
    depth = 0
    cur: list[str] = []
    for ch in arg_text:
        if ch == "(":
            depth += 1
            cur.append(ch)
        elif ch == ")":
            depth = max(0, depth - 1)
            cur.append(ch)
        elif ch == "," and depth == 0:
            token = "".join(cur).strip()
            if token:
                args.append(token)
            cur = []
        else:
            cur.append(ch)
    tail = "".join(cur).strip()
    if tail:
        args.append(tail)
    return args


def normalize_name(name: str) -> str:
    return re.sub(r"\s+", "", name).strip("&").lower()


@dataclass
class PyfRoutine:
    name: str
    args: list[str]
    out_args: set[str]
    return_args: set[str]
    intents: dict[str, str]


@dataclass
class FortranRoutine:
    name: str
    kind: str
    args: list[str]
    intents: dict[str, str]
    optional_args: set[str]
    file: Path


@dataclass
class F2pyWrapper:
    file: Path
    symbol: str
    source_routine: str
    signature_args: list[str]
    call_args: list[str]


@dataclass
class PythonCall:
    file: Path
    function_name: str
    symbol: str
    signature_args: list[str]
    call_positional_args: list[str]
    call_keyword_args: list[str]


@dataclass
class OrderIssue:
    python_api_fn: str
    pyf_symbol: str
    issue: str
    python_order: list[str]
    pyf_order: list[str]


INTERNAL_ONLY_ARGS = {
    "sclr_dim_transport",
    "edsclr_dim_transport",
    "hydromet_dim_transport",
}

MISSING_SOURCE_MAPPING_EXCEPTIONS = {
    "f2py_lapack_gtsv",
    "f2py_lapack_gtsvx",
    "f2py_lapack_gbsv",
    "f2py_lapack_gbsvx",
    "f2py_lapack_potrf",
    "f2py_lapack_poequ",
    "f2py_lapack_laqsy",
    "f2py_lapack_syev",
    "get_err_code",
    "f2py_set_stats_var_data",
    "f2py_get_param_names",
    "f2py_get_stats_config",
    "f2py_get_stats_var_meta",
    "f2py_get_stats_var_data",
    "f2py_set_simplified_radiation_params",
}

ARG_DIFF_EXCEPTIONS = {
    "stats",
}

DERIVED_SOURCE_OUTPUT_EXCEPTIONS = {
    # f2py_pdf_closure_driver exposes values derived from source outputs after
    # pdf_closure_driver returns.
    "f2py_pdf_closure_driver": {
        "rc_coef",
        "wp2rcp",
        "rtprcp",
        "rcp2",
        "skw_velocity",
        "cloud_frac_zm",
        "ice_supersat_frac_zm",
        "rtm_zm",
        "thlm_zm",
        "rcm_zm",
        "rcm_supersat_adj",
        "sclrprcp",
    },
    # f2py_compute_gamma_skw calls compute_gamma_skw on both vertical grids.
    "f2py_compute_gamma_skw": {
        "gamma_skw_fnc_zt",
    },
}


def _is_internal_only_arg(arg: str) -> bool:
    return arg in INTERNAL_ONLY_ARGS or arg.endswith("_dim_transport")


def _is_missing_source_mapping_exception(symbol: str) -> bool:
    return symbol in MISSING_SOURCE_MAPPING_EXCEPTIONS


def _is_ignored_symbol(symbol: str) -> bool:
    return symbol.startswith("f2py_lapack_")


def _should_ignore_missing_source_mapping(pc: PythonCall, pyf_sub: PyfRoutine) -> bool:
    if _is_missing_source_mapping_exception(pc.symbol):
        return True
    return not pc.signature_args and not pyf_sub.args


def _filter_arg_diff_exceptions(args: list[str]) -> list[str]:
    return [arg for arg in args if arg not in ARG_DIFF_EXCEPTIONS]


def _is_derived_source_output_exception(symbol: str, arg: str) -> bool:
    return arg in DERIVED_SOURCE_OUTPUT_EXCEPTIONS.get(symbol, set())


def _classify_public_api_signature_issue(actual: list[str], expected: list[str]) -> str | None:
    filtered_actual = _filter_arg_diff_exceptions(actual)
    filtered_expected = _filter_arg_diff_exceptions(expected)

    if filtered_actual == filtered_expected:
        return None

    if set(filtered_actual) == set(filtered_expected):
        return "PUBLIC_API_SIGNATURE_ORDER_MISMATCH"

    return "PUBLIC_API_SIGNATURE_ARG_MISMATCH"


def _python_api_fn_label(repo_root: Path, pc: PythonCall) -> str:
    return f"{pc.file.relative_to(repo_root)}::{pc.function_name}"


def _find_generic_dispatch_wrappers(repo_root: Path, py_calls: list[PythonCall]) -> set[str]:
    by_function: dict[str, set[str]] = {}
    for pc in py_calls:
        label = _python_api_fn_label(repo_root, pc)
        by_function.setdefault(label, set()).add(pc.symbol)
    return {
        label for label, symbols in by_function.items()
        if len(symbols) > 1
    }


def _extract_routine_blocks(merged: str) -> list[tuple[str, str, str]]:
    blocks: list[tuple[str, str, str]] = []
    starts = list(
        re.finditer(
            r"(?im)^\s*(subroutine|(?:[a-z0-9_(),=\s]+)?function)\s+([a-z0-9_]+)\s*\(",
            merged,
        )
    )
    for i, mt in enumerate(starts):
        start = mt.start()
        end = starts[i + 1].start() if i + 1 < len(starts) else len(merged)
        kind = "function" if "function" in mt.group(1).lower() else "subroutine"
        blocks.append((kind, mt.group(2).lower(), merged[start:end]))
    return blocks


def _extract_subroutine_blocks(merged: str) -> list[tuple[str, str]]:
    return [(name, block) for kind, name, block in _extract_routine_blocks(merged) if kind == "subroutine"]


def _declared_arg_names(declaration_text: str) -> list[str]:
    names: list[str] = []
    for raw in split_args(declaration_text):
        token = re.sub(r"\([^)]*\)", "", raw.strip()).strip()
        token = token.split("=", 1)[0].strip()
        if re.match(r"^[a-z_][a-z0-9_]*$", token, flags=re.I):
            names.append(normalize_name(token))
    return names


def parse_fortran_arg_intents(block: str) -> dict[str, str]:
    intents: dict[str, str] = {}
    for mt in re.finditer(
        r"(?im)^\s*[^\n]*\bintent\s*\(\s*(inout|in\s*,\s*out|out|in)\s*\)[^\n]*::\s*(.+)$",
        block,
    ):
        intent = re.sub(r"\s+", "", mt.group(1).lower())
        if intent == "in,out":
            intent = "inout"
        for name in _declared_arg_names(mt.group(2)):
            intents.setdefault(name, intent)
    return intents


def parse_fortran_optional_args(block: str) -> set[str]:
    optional_args: set[str] = set()
    for mt in re.finditer(
        r"(?im)^\s*[^\n]*\boptional\b[^\n]*::\s*(.+)$",
        block,
    ):
        for name in _declared_arg_names(mt.group(1)):
            optional_args.add(name)
    return optional_args


def public_fortran_args(routine: FortranRoutine, pyf_out_args: set[str]) -> list[str]:
    return [
        arg for arg in routine.args
        if routine.intents.get(arg) != "out"
        and not (arg not in routine.intents and arg in pyf_out_args)
    ]


def expected_public_args_for_python(
    routine: FortranRoutine,
    pyf_out_args: set[str],
    python_signature_args: list[str],
) -> list[str]:
    python_arg_set = set(python_signature_args)
    expected: list[str] = []
    for arg in public_fortran_args(routine, pyf_out_args):
        if arg in routine.optional_args and arg not in python_arg_set:
            continue
        expected.append(arg)
    return expected


def parse_pyf_routines(repo_root: Path) -> dict[str, PyfRoutine]:
    out: dict[str, PyfRoutine] = {}
    pyf = repo_root / "clubb_python_api" / "clubb_f2py.pyf"
    merged = merge_fortran_lines(pyf.read_text(errors="ignore"))
    for name, block in _extract_subroutine_blocks(merged):
        sig_match = re.search(
            rf"(?im)^\s*subroutine\s+{re.escape(name)}\s*\(([^)]*)\)",
            block,
        )
        args: list[str] = []
        if sig_match:
            args = [normalize_name(a) for a in split_args(sig_match.group(1)) if a.strip()]

        intents = parse_fortran_arg_intents(block)
        out_args = {arg for arg, intent in intents.items() if intent == "out"}
        return_args = {
            arg for arg, intent in intents.items()
            if intent in {"out", "inout"}
        }

        out[name] = PyfRoutine(
            name=name,
            args=args,
            out_args=out_args,
            return_args=return_args,
            intents=intents,
        )

    return out


def parse_fortran_routines(repo_root: Path) -> dict[str, FortranRoutine]:
    """Parse current Fortran source routine signatures in src."""
    out: dict[str, FortranRoutine] = {}
    src_dir = repo_root / "src"
    for path in sorted(src_dir.rglob("*.F90")):
        if "G_unit_test" in path.parts or "G_unit_tests" in path.parts:
            continue
        merged = merge_fortran_lines(path.read_text(errors="ignore"))
        for kind, name, block in _extract_routine_blocks(merged):
            sig_match = re.search(
                rf"(?im)^\s*(?:subroutine|(?:[a-z0-9_(),=\s]+)?function)\s+{re.escape(name)}\s*\(([^)]*)\)",
                block,
            )
            if not sig_match:
                continue
            args = [normalize_name(a) for a in split_args(sig_match.group(1)) if a.strip()]
            intents = parse_fortran_arg_intents(block)
            optional_args = parse_fortran_optional_args(block)
            out[name] = FortranRoutine(
                name=name,
                kind=kind,
                args=args,
                intents=intents,
                optional_args=optional_args,
                file=path,
            )
    return out


def parse_fortran_interfaces(repo_root: Path) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    src_dir = repo_root / "src"
    for path in sorted(src_dir.rglob("*.F90")):
        if "G_unit_test" in path.parts or "G_unit_tests" in path.parts:
            continue
        merged = merge_fortran_lines(path.read_text(errors="ignore"))
        for mt in re.finditer(
            r"(?ims)^\s*interface\s+([a-z0-9_]+)\s*(.*?)^\s*end\s+interface(?:\s+\1)?\s*$",
            merged,
        ):
            interface_name = normalize_name(mt.group(1))
            body = mt.group(2)
            members: list[str] = []
            for proc_mt in re.finditer(r"(?im)^\s*module\s+procedure\s+(.+)$", body):
                members.extend(
                    normalize_name(name) for name in split_args(proc_mt.group(1)) if name.strip()
                )
            if members:
                out[interface_name] = members
    return out


def _extract_fortran_invocations(
    block: str, fortran: dict[str, FortranRoutine], interfaces: dict[str, list[str]], symbol: str
) -> list[tuple[str, list[str]]]:
    call_candidates: list[tuple[str, list[str]]] = []

    for call_match in re.finditer(r"(?im)\bcall\s+([a-z0-9_]+)\b\s*(\()?", block):
        candidate = normalize_name(call_match.group(1))
        if candidate == symbol:
            continue
        if call_match.group(2):
            arg_start = call_match.end()
            depth = 1
            idx = arg_start
            while idx < len(block) and depth > 0:
                if block[idx] == "(":
                    depth += 1
                elif block[idx] == ")":
                    depth -= 1
                idx += 1
            call_args = [
                _normalize_fortran_call_arg(arg)
                for arg in split_args(block[arg_start:idx - 1])
                if arg.strip()
            ]
        else:
            call_args = []
        if candidate in fortran or candidate in interfaces:
            call_candidates.append((candidate, call_args))

    for fn_match in re.finditer(
        r"(?im)^\s*[a-z0-9_%]+(?:\s*\([^=\n]*\))?\s*=\s*([a-z0-9_]+)\s*\(",
        block,
    ):
        candidate = normalize_name(fn_match.group(1))
        if candidate == symbol or (candidate not in fortran and candidate not in interfaces):
            continue
        arg_start = fn_match.end()
        depth = 1
        idx = arg_start
        while idx < len(block) and depth > 0:
            if block[idx] == "(":
                depth += 1
            elif block[idx] == ")":
                depth -= 1
            idx += 1
        call_args = [
            _normalize_fortran_call_arg(arg)
            for arg in split_args(block[arg_start:idx - 1])
            if arg.strip()
        ]
        call_candidates.append((candidate, call_args))

    return call_candidates


def _resolve_interface_member(
    interface_name: str,
    call_args: list[str],
    fortran: dict[str, FortranRoutine],
    interfaces: dict[str, list[str]],
) -> str:
    members = interfaces.get(interface_name, [])
    if not members:
        return ""

    nargs = len(call_args)
    scored: list[tuple[int, int, int, str]] = []
    for member in members:
        routine = fortran.get(member)
        if routine is None:
            continue
        required = len([arg for arg in routine.args if arg not in routine.optional_args])
        total = len(routine.args)
        compatible = required <= nargs <= total
        score = (
            int(compatible),
            min(nargs, total),
            required,
            member,
        )
        scored.append(score)

    if not scored:
        return ""

    scored.sort(reverse=True)
    return scored[0][3]


def _normalize_fortran_call_arg(arg: str) -> str:
    token = arg.strip()
    if "=" in token:
        lhs, _rhs = token.split("=", 1)
        return normalize_name(lhs)
    return normalize_name(re.sub(r"\([^)]*\)", "", token).strip())


def parse_f2py_wrappers(
    repo_root: Path,
    fortran: dict[str, FortranRoutine],
    interfaces: dict[str, list[str]],
) -> dict[str, F2pyWrapper]:
    """Map F2PY wrapper subroutines to the Fortran routine they call."""
    out: dict[str, F2pyWrapper] = {}
    wrappers_dir = repo_root / "clubb_python_api" / "f2py_fortran_wrappers"
    for path in sorted(wrappers_dir.glob("*wrapper.F90")):
        merged = merge_fortran_lines(path.read_text(errors="ignore"))
        for symbol, block in _extract_subroutine_blocks(merged):
            sig_match = re.search(
                rf"(?im)^\s*subroutine\s+{re.escape(symbol)}\s*\(([^)]*)\)",
                block,
            )
            if not sig_match:
                continue
            signature_args = [
                normalize_name(a) for a in split_args(sig_match.group(1)) if a.strip()
            ]

            call_candidates = _extract_fortran_invocations(block, fortran, interfaces, symbol)
            if call_candidates:
                source_base = symbol.removeprefix("f2py_")
                resolved_candidates: list[tuple[str, list[str], str]] = []
                for candidate, call_args in call_candidates:
                    resolved_candidate = candidate
                    if candidate not in fortran:
                        resolved_candidate = _resolve_interface_member(
                            candidate, call_args, fortran, interfaces
                        )
                    if resolved_candidate:
                        resolved_candidates.append((resolved_candidate, call_args, candidate))
                if not resolved_candidates:
                    continue
                source_routine, call_args, source_name = max(
                    resolved_candidates,
                    key=lambda item: (
                        item[2] == source_base,
                        item[0] == source_base,
                        source_base.startswith(f"{item[2]}_"),
                        source_base.startswith(f"{item[0]}_"),
                        len(item[1]),
                        len(item[0]),
                    ),
                )
            else:
                source_routine = ""
                call_args = []
            if not source_routine:
                continue

            out[symbol] = F2pyWrapper(
                file=path,
                symbol=symbol,
                source_routine=source_routine,
                signature_args=signature_args,
                call_args=call_args,
            )
    return out


def infer_source_routine(symbol: str, fortran: dict[str, FortranRoutine]) -> str:
    """Infer source routine for direct f2py bindings without an adapter wrapper."""
    base = symbol.removeprefix("f2py_")
    suffixes = (
        "_same_grid",
        "_scalar_array",
        "_array_scalar",
        "_single_rhs_multiple_lhs",
        "_multiple_rhs_lhs",
        "_single_rhs_lhs",
        "_multiple_lhs",
        "_scalar",
        "_array",
        "_2d",
        "_1d",
    )
    aliases = [base]
    changed = True
    while changed:
        changed = False
        current = aliases[-1]
        for suffix in suffixes:
            if current.endswith(suffix):
                aliases.append(current[: -len(suffix)])
                changed = True
                break

    for alias in aliases:
        if alias in fortran:
            return alias
        if f"{alias}_api" in fortran:
            return f"{alias}_api"

    candidates = [
        name for name in fortran
        for alias in aliases
        if alias.startswith(f"{name}_") or name.startswith(f"{alias}_")
    ]
    if not candidates:
        return ""

    candidates.sort(key=len, reverse=True)
    if len(candidates) > 1 and len(candidates[0]) == len(candidates[1]):
        return ""
    return candidates[0]


def _extract_call_arg_ident(expr: ast.AST) -> str | None:
    if isinstance(expr, ast.Name):
        return normalize_name(expr.id)

    if isinstance(expr, ast.UnaryOp):
        return _extract_call_arg_ident(expr.operand)

    if isinstance(expr, ast.BinOp):
        left = _extract_call_arg_ident(expr.left)
        right = _extract_call_arg_ident(expr.right)
        if left is not None and isinstance(expr.right, ast.Constant):
            return left
        if right is not None and isinstance(expr.left, ast.Constant):
            return right
        return None

    if isinstance(expr, ast.Call):
        if expr.args:
            return _extract_call_arg_ident(expr.args[0])
        return None

    return None


def parse_python_calls(repo_root: Path) -> list[PythonCall]:
    out: list[PythonCall] = []
    py_roots = [
        repo_root / "clubb_python_api" / "clubb_python" / "CLUBB_core",
        repo_root / "clubb_python_api" / "clubb_python",
    ]
    seen_files: set[Path] = set()
    files: list[Path] = []
    for root in py_roots:
        if not root.exists():
            continue
        for path in sorted(root.glob("*.py")):
            if path.name == "__init__.py" or path in seen_files:
                continue
            seen_files.add(path)
            files.append(path)

    for path in files:
        text = path.read_text(errors="ignore")
        try:
            tree = ast.parse(text)
        except SyntaxError:
            continue

        for node in tree.body:
            if not isinstance(node, ast.FunctionDef) or node.name.startswith("_"):
                continue
            signature_args = [
                normalize_name(arg.arg)
                for arg in [*node.args.args, *node.args.kwonlyargs]
            ]

            for call in ast.walk(node):
                if not isinstance(call, ast.Call):
                    continue
                if not isinstance(call.func, ast.Attribute):
                    continue
                if not isinstance(call.func.value, ast.Name) or call.func.value.id != "clubb_f2py":
                    continue

                out.append(
                    PythonCall(
                        file=path,
                        function_name=node.name,
                        symbol=call.func.attr.lower(),
                        signature_args=signature_args,
                        call_positional_args=[
                            ident
                            for arg in call.args
                            if (ident := _extract_call_arg_ident(arg)) is not None
                        ],
                        call_keyword_args=[
                            normalize_name(kw.arg)
                            for kw in call.keywords
                            if kw.arg is not None
                        ],
                    )
                )

    return out


def _shared_order(lhs: list[str], rhs: list[str]) -> tuple[list[str], list[str]]:
    rhs_set = set(rhs)
    lhs_shared = [x for x in lhs if x in rhs_set]

    lhs_set = set(lhs)
    rhs_shared = [x for x in rhs if x in lhs_set]
    return lhs_shared, rhs_shared


def find_order_issues(repo_root: Path) -> list[OrderIssue]:
    return find_contract_issues(repo_root)


def _ordered_extras(actual: list[str], expected: list[str]) -> list[str]:
    expected_set = set(expected)
    return [arg for arg in actual if arg not in expected_set]


def _ordered_missing(actual: list[str], expected: list[str]) -> list[str]:
    actual_set = set(actual)
    return [arg for arg in expected if arg not in actual_set]


def _source_return_args(source: FortranRoutine) -> list[str]:
    return [
        arg for arg in source.args
        if source.intents.get(arg) in {"out", "inout"}
    ]


def _print_issue_details(issue: OrderIssue) -> None:
    filtered_actual = _filter_arg_diff_exceptions(issue.python_order)
    filtered_expected = _filter_arg_diff_exceptions(issue.pyf_order)

    if issue.issue.startswith("PUBLIC_API_SIGNATURE_ORDER_MISMATCH"):
        print(f"  ORDER_ACTUAL:   {filtered_actual}")
        print(f"  ORDER_EXPECTED: {filtered_expected}")
        return

    extra = _ordered_extras(filtered_actual, filtered_expected)
    missing = _ordered_missing(filtered_actual, filtered_expected)

    if extra:
        print(f"  EXTRA:   {extra}")
    if missing:
        print(f"  MISSING: {missing}")
    if not extra and not missing:
        print(f"  ORDER_ACTUAL:   {filtered_actual}")
        print(f"  ORDER_EXPECTED: {filtered_expected}")


def find_contract_issues(repo_root: Path) -> list[OrderIssue]:
    pyf = parse_pyf_routines(repo_root)
    py_calls = parse_python_calls(repo_root)
    fortran = parse_fortran_routines(repo_root)
    interfaces = parse_fortran_interfaces(repo_root)
    f2py_wrappers = parse_f2py_wrappers(repo_root, fortran, interfaces)
    generic_dispatch_wrappers = _find_generic_dispatch_wrappers(repo_root, py_calls)
    issues: list[OrderIssue] = []

    for pc in py_calls:
        if _is_ignored_symbol(pc.symbol):
            continue
        python_api_fn = _python_api_fn_label(repo_root, pc)
        pyf_sub = pyf.get(pc.symbol)
        if pyf_sub is None:
            issues.append(
                OrderIssue(
                    python_api_fn=python_api_fn,
                    pyf_symbol=pc.symbol,
                    issue="MISSING_PYF_SYMBOL",
                    python_order=pc.call_positional_args,
                    pyf_order=[],
                )
            )
            continue

        wrapper = f2py_wrappers.get(pc.symbol)
        if wrapper is None:
            source_routine = infer_source_routine(pc.symbol, fortran)
            if not source_routine:
                if not _should_ignore_missing_source_mapping(pc, pyf_sub):
                    issues.append(
                        OrderIssue(
                            python_api_fn=python_api_fn,
                            pyf_symbol=pc.symbol,
                            issue="MISSING_FORTRAN_SOURCE_MAPPING",
                            python_order=pc.signature_args,
                            pyf_order=pyf_sub.args,
                        )
                    )
                continue
            wrapper = F2pyWrapper(
                file=Path("<direct-pyf>"),
                symbol=pc.symbol,
                source_routine=source_routine,
                signature_args=pyf_sub.args,
                call_args=pyf_sub.args,
            )

        source = fortran.get(wrapper.source_routine)
        if source is None:
            issues.append(
                OrderIssue(
                    python_api_fn=python_api_fn,
                    pyf_symbol=pc.symbol,
                    issue=f"MISSING_FORTRAN_SOURCE::{wrapper.source_routine}",
                    python_order=pc.signature_args,
                    pyf_order=pyf_sub.args,
                )
            )
            continue

        if source.kind == "subroutine":
            source_return_args = _source_return_args(source)
            source_return_arg_set = set(source_return_args)
            wrapper_call_arg_set = set(wrapper.call_args)
            extra_return_args = [
                arg for arg in pyf_sub.args
                if arg in pyf_sub.return_args
                and arg not in source_return_arg_set
                and arg not in wrapper_call_arg_set
                and not _is_internal_only_arg(arg)
                and not _is_derived_source_output_exception(pc.symbol, arg)
            ]
            if extra_return_args:
                issues.append(
                    OrderIssue(
                        python_api_fn=python_api_fn,
                        pyf_symbol=pc.symbol,
                        issue=(
                            "F2PY_EXTRA_SOURCE_OUTPUT_ARG::"
                            f"{source.file.relative_to(repo_root)}::{source.name}"
                        ),
                        python_order=extra_return_args,
                        pyf_order=source_return_args,
                    )
                )

        leaked_internal_args = [arg for arg in pc.signature_args if _is_internal_only_arg(arg)]
        if leaked_internal_args:
            issues.append(
                OrderIssue(
                    python_api_fn=python_api_fn,
                    pyf_symbol=pc.symbol,
                    issue="PUBLIC_API_FORBIDDEN_INTERNAL_ARG",
                    python_order=leaked_internal_args,
                    pyf_order=[],
                )
            )

        expected_public_args = expected_public_args_for_python(
            source, pyf_sub.out_args, pc.signature_args
        )
        mismatch_kind = _classify_public_api_signature_issue(
            pc.signature_args, expected_public_args
        )
        if mismatch_kind is not None:
            if (
                mismatch_kind == "PUBLIC_API_SIGNATURE_ARG_MISMATCH"
                and python_api_fn in generic_dispatch_wrappers
            ):
                mismatch_kind = "GENERIC_DISPATCH_WRAPPER"
            issues.append(
                OrderIssue(
                    python_api_fn=python_api_fn,
                    pyf_symbol=pc.symbol,
                    issue=(
                        f"{mismatch_kind}::"
                        f"{source.file.relative_to(repo_root)}::{source.name}"
                    ),
                    python_order=pc.signature_args,
                    pyf_order=expected_public_args,
                )
            )

        if pyf_sub.args != wrapper.signature_args:
            wrapper_file = (
                str(wrapper.file)
                if wrapper.file.is_absolute()
                else str(wrapper.file.relative_to(repo_root))
                if wrapper.file != Path("<direct-pyf>")
                else "<direct-pyf>"
            )
            issues.append(
                OrderIssue(
                    python_api_fn=python_api_fn,
                    pyf_symbol=pc.symbol,
                    issue=f"PYF_WRAPPER_SIGNATURE_MISMATCH::{wrapper_file}",
                    python_order=pyf_sub.args,
                    pyf_order=wrapper.signature_args,
                )
            )

        pyf_positional_order = [
            arg for arg in pyf_sub.args
            if (
                arg not in pyf_sub.out_args
                and arg not in pc.call_keyword_args
                and not _is_internal_only_arg(arg)
            )
        ]

        signature_forward_args = [
            arg for arg in pc.signature_args
            if arg in pyf_positional_order and arg not in pc.call_keyword_args
        ]
        expected_forward_args = [arg for arg in pyf_positional_order if arg in pc.signature_args]
        if signature_forward_args != expected_forward_args:
            issues.append(
                OrderIssue(
                    python_api_fn=python_api_fn,
                    pyf_symbol=pc.symbol,
                    issue="PUBLIC_API_TO_PYF_ORDER_MISMATCH",
                    python_order=signature_forward_args,
                    pyf_order=expected_forward_args,
                )
            )

        call_shared, pyf_shared = _shared_order(pc.call_positional_args, pyf_positional_order)
        if call_shared != pyf_shared:
            issues.append(
                OrderIssue(
                    python_api_fn=python_api_fn,
                    pyf_symbol=pc.symbol,
                    issue="CALL_ORDER",
                    python_order=call_shared,
                    pyf_order=pyf_shared,
                )
            )
    deduped_issues: list[OrderIssue] = []
    seen_issue_keys: set[tuple[str, str, str, tuple[str, ...], tuple[str, ...]]] = set()
    for issue in issues:
        key = (
            issue.python_api_fn,
            issue.pyf_symbol,
            issue.issue,
            tuple(issue.python_order),
            tuple(issue.pyf_order),
        )
        if key in seen_issue_keys:
            continue
        seen_issue_keys.add(key)
        deduped_issues.append(issue)
    return deduped_issues


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Check Python/F2PY wrapper argument contracts against Fortran."
    )
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args()

    root = args.repo_root.resolve()
    issues = find_order_issues(root)
    if not issues:
        print("No argument contract issues found.")
        return 0

    for issue in issues:
        issue_type = issue.issue.split("::", 1)[0]
        prefix = "WARNING" if issue_type == "GENERIC_DISPATCH_WRAPPER" else "ERROR"
        print(f"{prefix} {issue.issue}: {issue.python_api_fn} -> {issue.pyf_symbol}")
        _print_issue_details(issue)
    print("\nSummary:")
    for issue_type, count in sorted(Counter(issue.issue.split("::", 1)[0] for issue in issues).items()):
        print(f"  {issue_type}: {count}")
    blocking_issues = [
        issue for issue in issues
        if issue.issue.split("::", 1)[0] != "GENERIC_DISPATCH_WRAPPER"
    ]
    return 1 if blocking_issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
