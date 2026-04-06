#!/usr/bin/env python3
"""Check that Python wrapper calls preserve F2PY argument order."""

from __future__ import annotations

import argparse
import ast
import re
import sys
from dataclasses import dataclass
from pathlib import Path


def _strip_fortran_comments(line: str) -> str:
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
    return re.sub(r"\s+", "", name).lower()


@dataclass
class PyfRoutine:
    name: str
    args: list[str]
    out_args: set[str]


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


def _extract_subroutine_blocks(merged: str) -> list[tuple[str, str]]:
    blocks: list[tuple[str, str]] = []
    starts = list(re.finditer(r"(?im)^\s*subroutine\s+([a-z0-9_]+)\s*\(([^)]*)\)", merged))
    for i, mt in enumerate(starts):
        start = mt.start()
        end = starts[i + 1].start() if i + 1 < len(starts) else len(merged)
        blocks.append((mt.group(1).lower(), merged[start:end]))
    return blocks


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
    pyf = parse_pyf_routines(repo_root)
    py_calls = parse_python_calls(repo_root)
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

        pyf_positional_order = [
            arg for arg in pyf_sub.args
            if arg not in pyf_sub.out_args and arg not in pc.call_keyword_args
        ]

        sig_shared, pyf_sig_shared = _shared_order(pc.signature_args, pyf_positional_order)
        if sig_shared != pyf_sig_shared:
            issues.append(
                OrderIssue(
                    python_api_fn=python_api_fn,
                    pyf_symbol=pc.symbol,
                    issue="PYTHON_API_ORDER",
                    python_order=sig_shared,
                    pyf_order=pyf_sig_shared,
                )
            )
            continue

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
            continue
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description="Check Python wrapper argument order against F2PY.")
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args()

    root = args.repo_root.resolve()
    issues = find_order_issues(root)
    if not issues:
        print("No argument-order issues found.")
        return 0

    for issue in issues:
        print(f"{issue.issue}: {issue.python_api_fn} -> {issue.pyf_symbol}")
        print(f"  python: {issue.python_order}")
        print(f"  pyf:    {issue.pyf_order}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
