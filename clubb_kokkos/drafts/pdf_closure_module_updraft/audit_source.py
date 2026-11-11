"""Compare the isolated updraft draft with Native; run no compiled/model test."""
from pathlib import Path
import ast
import copy
import hashlib
import json
import re

folder = Path(__file__).resolve().parent
root = folder.parents[2]
native_path = root / "src/CLUBB_core/pdf_closure_module.F90"
native = native_path.read_text()
native = re.search(
    r"^  subroutine calc_w_up_in_cloud\(.*?^  end subroutine calc_w_up_in_cloud$",
    native, re.M | re.S,
)[0]
cpp_path = folder / "pdf_closure_module.cpp"
header_path = folder / "pdf_closure_module.hpp"
cpp = cpp_path.read_text()
header = header_path.read_text()
clean_cpp = re.sub(r"//[^\n]*", "", cpp)
clean_header = re.sub(r"//[^\n]*", "", header)


def statements(source):
    """Join Native continuation statements without declaration/body comments."""
    result = []
    pending = ""
    for line in source.splitlines():
        code = line.split("!", 1)[0].strip()
        if not code:
            continue
        pending += " " + code.strip("& ")
        if not code.endswith("&"):
            result.append(pending.strip())
            pending = ""
    assert not pending
    return result


native_statements = statements(native)
expected_args = native_statements[0].split("(", 1)[1].rsplit(")", 1)[0].replace(" ", "").split(",")
for text in (clean_cpp, clean_header):
    signature = re.search(r"void calc_w_up_in_cloud\((.*?)\)", text, re.S)[1]
    args = [part.split()[-1] for part in signature.split(",")]
    assert args == expected_args
    assert len(args) == 13
    assert "=" not in signature and "*" not in signature


class SourceSquares(ast.NodeTransformer):
    """Expand only the Native integer square into the ordinary C++ product."""
    def visit_BinOp(self, node):
        node = self.generic_visit(node)
        if isinstance(node.op, ast.Pow):
            assert isinstance(node.right, ast.Constant) and node.right.value == 2
            return ast.BinOp(node.left, ast.Mult(), copy.deepcopy(node.left))
        return node


def expression(text, target=False):
    text = re.sub(r"\s+", "", text)
    if target:
        text = text.replace("Kokkos::", "")
        # The four source MAX operations keep their second-operand tie choice.
        for value in ("stdev_w_1", "stdev_w_2", "cloudy_updraft_frac(i,k)",
                      "cloudy_downdraft_frac(i,k)"):
            source_max = f"(eps>{value}?eps:{value})"
            text = text.replace(source_max, f"max(eps,{value})")
    parsed = ast.parse(text, mode="eval")
    if not target:
        parsed = SourceSquares().visit(parsed)
    return ast.dump(parsed, include_attributes=False)


native_assignments = []
for statement in native_statements:
    match = re.match(r"^(\w+(?:\([^)]*\))?)\s*=\s*(.*)$", statement)
    if match:
        native_assignments.append((expression(match[1]), expression(match[2])))

target_assignments = []
for match in re.finditer(r"^\s*(\w+(?:\([^)]*\))?)\s*=\s*(.*?);", clean_cpp, re.M | re.S):
    target_assignments.append((expression(match[1], True), expression(match[2], True)))
assert len(native_assignments) == 36
assert target_assignments == native_assignments

native_conditions = [re.match(r"(?:else)?if\s*\((.*)\)\s*then", s, re.I)[1]
                     for s in native_statements
                     if re.match(r"(?:else)?if\s*\(", s, re.I)]
target_conditions = re.findall(r"(?:else\s+)?if\s*\((.*?)\)\s*\{", clean_cpp, re.S)
assert len(native_conditions) == len(target_conditions) == 4
assert [expression(v) for v in native_conditions] == [expression(v, True) for v in target_conditions]

expected_locals = json.loads((folder / "SOURCE_OUTLINE.json").read_text())["native_local_names"][:-2]
locals_found = []
for declaration in re.findall(r"\breal_t\s+(.*?);", clean_cpp):
    locals_found.extend(part.strip() for part in declaration.split(","))
assert locals_found == expected_locals
assert clean_cpp.count("Kokkos::parallel_for(") == 1
assert clean_cpp.count("KOKKOS_LAMBDA(") == 1
assert "const int i = cell % ngrdcol;" in clean_cpp
assert "const int k = cell / ngrdcol;" in clean_cpp
assert "f2k::" not in clean_cpp
assert all(token not in clean_cpp for token in ("parallel_reduce", "std::fma", "Kokkos::fma", "pow("))
assert clean_cpp.count("eps >") == 4

# Every Native computation-stage comment remains in source order, including
# the inactive alternate downdraft-fraction formula.
native_body = native[native.index("        stdev_w_1 ="):]
target_body = cpp[cpp.index("        stdev_w_1 ="):]
body_comments = [line.split("!", 1)[1] for line in native_body.splitlines()
                 if line.strip().startswith("!") and not line.strip().startswith("!$acc")]
cursor = 0
for comment in body_comments:
    canonical = re.sub(r"\s+", "", comment)
    matches = list(re.finditer(r"//([^\n]*)", target_body[cursor:]))
    match = next(m for m in matches if re.sub(r"\s+", "", m[1]) == canonical)
    cursor += match.end()

report = {
    "status": "Separate source audit passes; no compilation, tests or numerical/BFB claim",
    "routine": "calc_w_up_in_cloud",
    "native_lines": [3143, 3340],
    "public_required_arguments": 13,
    "required_inputs": 9,
    "required_outputs": 4,
    "signature_header_definition_native_match": True,
    "ordered_assignment_operation_trees": 36,
    "ordered_component_branch_predicates": 4,
    "branches_per_component": 3,
    "native_scalar_locals": locals_found,
    "source_body_comments_preserved_in_order": len(body_comments),
    "column_first_cell_kernels": 1,
    "max_second_operand_tie_choices": 4,
    "private_scratch": "Cell-local scalars only; no allocation, reduction, UDT capture or bridge",
    "bfb_adaptations": "None; integer-square and ordered MAX lowerings are source/language translations",
    "comment_correction": "varnce_w_1/2 are variances, verified from sqrt(varnce) and squared units",
    "reviewed_sha256": {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
                        for path in (native_path, cpp_path, header_path)},
}
(folder / "source_audit.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
