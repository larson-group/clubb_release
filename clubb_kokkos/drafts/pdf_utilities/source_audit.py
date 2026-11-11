"""Audit Native/draft contracts and arithmetic trees without executing physics.

This source comparison does not compile or test the draft. It compares exact
argument order, assignments and branch predicates after syntax-only mappings,
including expansion of each Native integer square into a parenthesized product.
"""
from pathlib import Path
import ast
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
NATIVE_PATH = ROOT / "src/CLUBB_core/pdf_utilities.F90"
TARGET_PATH = HERE / "pdf_utilities.hpp"
NATIVE = NATIVE_PATH.read_text()
TARGET = TARGET_PATH.read_text()
OUTLINE = json.loads((HERE / "source_outline.json").read_text())


def f_statements(text):
    """Join Native continuations, omitting comments and compiler directives."""
    result = []
    pending = ""
    for raw in text.splitlines():
        code = raw.split("!", 1)[0].strip()
        if not code:
            continue
        pending += " " + code.strip("& ")
        if not code.endswith("&"):
            result.append(pending.strip())
            pending = ""
    if pending:
        raise ValueError("unfinished Native statement")
    return result


def cpp_code(text):
    """Remove line comments before scanning target definitions/expressions."""
    return re.sub(r"//[^\n]*", "", text)


def definitions(text):
    """Collect draft definitions without confusing the forward declaration."""
    text = cpp_code(text)
    result = []
    pattern = re.compile(r"(?:KOKKOS_INLINE_FUNCTION|inline)\s+(real_t|double|void|Real2D)\s+(\w+)\s*\(([^{};]*?)\)\s*\{", re.S)
    for match in pattern.finditer(text):
        depth = 1
        end = match.end()
        while depth:
            if text[end] == "{":
                depth += 1
            elif text[end] == "}":
                depth -= 1
            end += 1
        args = [arg.strip().split()[-1] for arg in match.group(3).split(",")]
        result.append({"name": match.group(2), "args": args, "return_type": match.group(1), "body": text[match.end():end - 1]})
    return result


class NativeSquares(ast.NodeTransformer):
    """Translate only integer **2; preserve every surrounding arithmetic group."""
    def visit_BinOp(self, node):
        self.generic_visit(node)
        if isinstance(node.op, ast.Pow) and isinstance(node.right, ast.Constant) and node.right.value == 2:
            return ast.BinOp(left=node.left, op=ast.Mult(), right=node.left)
        return node

    def visit_Call(self, node):
        self.generic_visit(node)
        if isinstance(node.func, ast.Name) and node.func.id == "tiny":
            return ast.Name(id="native_tiny", ctx=ast.Load())
        return node


def tree(expr):
    """Normalize syntax/intrinsic spelling, retaining arithmetic association."""
    expr = expr.replace("Kokkos::", "")
    expr = re.sub(r"std::numeric_limits<real_t>::min\(\)", "native_tiny", expr)
    expr = re.sub(r"tiny\([^()]*\)", "native_tiny", expr)
    expr = re.sub(r"static_cast<double>\(([^()]*)\)", r"as_double(\1)", expr)
    expr = re.sub(r"real\(\s*max_mag_correlation,\s*kind\s*=\s*dp\s*\)", "as_double(max_mag_correlation)", expr, flags=re.I)
    expr = expr.replace(".and.", " and ").replace(".or.", " or ").replace("&&", " and ").replace("||", " or ")
    node = ast.parse(" ".join(expr.split()), mode="eval")
    node = NativeSquares().visit(node)
    return ast.dump(node, include_attributes=False)


def assignments(statements):
    """Read scientific assignments while ignoring declarations/loop indices."""
    result = []
    for statement in statements:
        match = re.match(r"^(\w+(?:\([^()]*\))?)\s*=(?!=)\s*(.*)$", statement.strip(), re.S)
        if match:
            result.append((tree(match.group(1)), tree(match.group(2))))
    return result


def cpp_assignments(body):
    """Find statements beginning with the Native scalar/array destination."""
    return assignments(match.group(0).strip().rstrip(";") for match in re.finditer(r"(?m)^\s*\w+(?:\([^()]*\))?\s*=(?!=)[^;]*;", body))


def f_conditions(statements):
    result = []
    for statement in statements:
        match = re.match(r"(?:if|elseif)\s*\((.*)\)\s*then$", statement, re.I)
        if match:
            result.append(tree(match.group(1)))
    return result


def cpp_conditions(body):
    result = []
    for match in re.finditer(r"\bif\s*\(", body):
        depth = 1
        end = match.end()
        while depth:
            if body[end] == "(":
                depth += 1
            elif body[end] == ")":
                depth -= 1
            end += 1
        result.append(tree(body[match.end():end - 1]))
    return result


entries = definitions(TARGET)
expected_names = [routine["name"] for routine in OUTLINE["routines"]]
if [entry["name"] for entry in entries[:19]] != expected_names:
    raise ValueError("Native definition order differs")
if [entry["name"] for entry in entries[19:]] != ["compute_variance_binormal", "calc_corr_chi_x", "calc_corr_eta_x"]:
    raise ValueError("unexpected View adaptation inventory")
reviews = []
for expected, entry in zip(OUTLINE["routines"], entries):
    name = expected["name"]
    if entry["args"] != expected["arguments"]:
        raise ValueError((name, "arguments", entry["args"], expected["arguments"]))
    start = re.search(r"^  (?:elemental function|subroutine) " + name + r"\b", NATIVE, re.M).start()
    end = NATIVE.index("end " + expected["kind"] + " " + name, start)
    source_statements = f_statements(NATIVE[start:end])
    native_assigns = assignments(source_statements)
    target_assigns = cpp_assignments(entry["body"])
    if native_assigns != target_assigns:
        raise ValueError((name, "assignment expression trees", native_assigns, target_assigns))
    if f_conditions(source_statements) != cpp_conditions(entry["body"]):
        raise ValueError((name, "branch predicates"))
    if name.endswith("_dp") and entry["return_type"] != "double":
        raise ValueError((name, "double result"))
    reviews.append({"name": name, "arguments": len(entry["args"]), "assignments": len(native_assigns), "branch_predicates": len(f_conditions(source_statements)), "signature_assignment_branch_trees_match": True})

helper_native = (ROOT / "src/CLUBB_core/advance_helper_module.F90").read_text()
start = helper_native.index("  function smooth_max_scalars(")
end = helper_native.index("  end function smooth_max_scalars", start)
helper_entries = definitions((HERE / "advance_helper_module.hpp").read_text())
if assignments(f_statements(helper_native[start:end])) != cpp_assignments(helper_entries[0]["body"]):
    raise ValueError("smooth_max scalar arithmetic differs")
if helper_entries[0]["args"] != ["input_var1", "input_var2", "smth_coef"]:
    raise ValueError("smooth_max scalar contract differs")

report = {
    "status": "source-audited draft; no compilation, execution, integration or BFB validation",
    "source": str(NATIVE_PATH.relative_to(ROOT)),
    "source_sha256": hashlib.sha256(NATIVE_PATH.read_bytes()).hexdigest(),
    "target_sha256": hashlib.sha256(TARGET_PATH.read_bytes()).hexdigest(),
    "native_routines": 19,
    "native_order_preserved": True,
    "routine_reviews": reviews,
    "device_scalar_definitions": 15,
    "array_definitions": 4,
    "needed_rank2_elemental_adaptations": [entry["name"] for entry in entries[19:]],
    "partial_helper": {"source": "src/CLUBB_core/advance_helper_module.F90:1723–1761", "name": "smooth_max_scalars", "argument_and_expression_trees_match": True, "generic_dispatch": "smooth_max scalar overload only"},
    "manual_review": [
        "Read all 1423 Native lines; JAX used for scope context only",
        "All 19 signatures, required arguments, return precision and Native relative order",
        "All assigned scientific expressions, branch predicates and complete grouping, including every integer square",
        "Tiny guard retained only where Native has it; no JAX sqrt clipping, log1p/expm1 rewrites or new denominator floors",
        "Ordered upper/lower correlation checks and negative hydrometeor variance check preserved",
        "Native comments/equations/references and argument meanings retained; stage spacing formatted separately",
        "Three conformable rank2 execution-space overloads return call-owned Real2D, permit strided Views and call scalar formulas",
        "No stats, UDT mutation, saved Native state, new runtime framework or scalar Fortran crossings",
        "Two sequential smooth_max calls and final quotient preserved in per-cell order",
        "No new explicit FMA or compiler-specific arithmetic fix; BFB evidence pending root integration"
    ],
    "material_exclusions": ["No shared production/CMake/generator edits", "No tests or builds", "Only scalar smooth_max subset of advance_helper_module included; root owns dependency dedup/promotion"],
    "outstanding": ["Root integration", "Compile and maintained short model comparisons", "Diagnose/record any demonstrated compiler-specific arithmetic differences"]
}
(HERE / "source_audit.json").write_text(json.dumps(report, indent=2) + "\n")
print("Source audit:19 contracts/orders, assignment trees and branch trees match; helper matches. Numerical validation pending.")
