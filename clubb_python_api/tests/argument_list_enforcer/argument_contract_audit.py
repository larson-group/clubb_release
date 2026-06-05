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


@dataclass
class FortranRoutine:
    name: str
    args: list[str]
    intents: dict[str, str]
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


def _is_internal_only_arg(arg: str) -> bool:
    return arg in INTERNAL_ONLY_ARGS or arg.endswith("_dim_transport")


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
            intents[name] = intent
    return intents


def public_fortran_args(routine: FortranRoutine, pyf_out_args: set[str]) -> list[str]:
    return [
        arg for arg in routine.args
        if routine.intents.get(arg) != "out"
        and not (arg not in routine.intents and arg in pyf_out_args)
    ]


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

        out_args: set[str] = set()
        for mt in re.finditer(r"(?im)^\s*[^\n]*\bintent\s*\(\s*out\s*\)[^\n]*::\s*(.+)$", block):
            for raw in split_args(mt.group(1)):
                token = re.sub(r"\([^)]*\)", "", raw.strip()).strip()
                if re.match(r"^[a-z_][a-z0-9_]*$", token, flags=re.I):
                    out_args.add(normalize_name(token))

        out[name] = PyfRoutine(name=name, args=args, out_args=out_args)

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
            out[name] = FortranRoutine(name=name, args=args, intents=intents, file=path)
    return out


def parse_f2py_wrappers(repo_root: Path, fortran: dict[str, FortranRoutine]) -> dict[str, F2pyWrapper]:
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

            call_candidates: list[tuple[str, list[str]]] = []
            for call_match in re.finditer(r"(?im)\bcall\s+([a-z0-9_]+)\b\s*(\()?", block):
                candidate = normalize_name(call_match.group(1))
                if candidate in fortran and candidate != symbol:
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
                            normalize_name(re.sub(r"\([^)]*\)", "", arg).strip())
                            for arg in split_args(block[arg_start:idx - 1])
                            if arg.strip()
                        ]
                    else:
                        call_args = []
                    call_candidates.append((candidate, call_args))
            if call_candidates:
                source_base = symbol.removeprefix("f2py_")
                source_routine, call_args = max(
                    call_candidates,
                    key=lambda item: (
                        item[0] == source_base,
                        source_base.startswith(f"{item[0]}_"),
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
            signature_args = [normalize_name(arg.arg) for arg in node.args.args]

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


def find_contract_issues(repo_root: Path) -> list[OrderIssue]:
    pyf = parse_pyf_routines(repo_root)
    py_calls = parse_python_calls(repo_root)
    fortran = parse_fortran_routines(repo_root)
    f2py_wrappers = parse_f2py_wrappers(repo_root, fortran)
    issues: list[OrderIssue] = []

    for pc in py_calls:
        python_api_fn = f"{pc.file.relative_to(repo_root)}::{pc.function_name}"
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

        expected_public_args = public_fortran_args(source, pyf_sub.out_args)
        if pc.signature_args != expected_public_args:
            issues.append(
                OrderIssue(
                    python_api_fn=python_api_fn,
                    pyf_symbol=pc.symbol,
                    issue=(
                        "PUBLIC_API_SIGNATURE_MISMATCH::"
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
    return issues


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
        print(f"{issue.issue}: {issue.python_api_fn} -> {issue.pyf_symbol}")
        print(f"  actual:   {issue.python_order}")
        print(f"  expected: {issue.pyf_order}")
    print("\nSummary:")
    for issue_type, count in sorted(Counter(issue.issue.split("::", 1)[0] for issue in issues).items()):
        print(f"  {issue_type}: {count}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
