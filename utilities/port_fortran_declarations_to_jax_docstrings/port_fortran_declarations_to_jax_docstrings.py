#!/usr/bin/env python3
"""Copy Fortran argument declarations into the docstrings of JAX ports.

Each JAX function in clubb_jax/src/CLUBB_core that ports a Fortran routine
documents that routine's Input, Input/Output and Output declaration blocks
in its docstring. This script finds the Fortran counterpart of every JAX
function, extracts those blocks, and checks or rewrites the docstring
section that holds them. It also lists each routine's Description and
References comments, from which -write_description writes the docstring
description together with recorded edits that adapt Fortran-only details.

Usage:
  port_fortran_declarations_to_jax_docstrings.py [-check]
      Report JAX functions whose Fortran declaration section is missing or
      differs from the Fortran source. Exits with 1 when any differ.
  port_fortran_declarations_to_jax_docstrings.py -write_declaration
      Insert or replace the declaration section of each docstring. Nothing
      else in the Python files changes; every rewritten file is parsed and
      compared with the original with docstrings removed before it is saved.
  port_fortran_declarations_to_jax_docstrings.py -write_description [-edits PATH]
      Rewrite the docstring description from the Fortran Description and
      References comments, a Returns block read from the code, and the
      per-function edits in custom_ports_of_fortran_descriptions.jsonc, then write
      the section.
  port_fortran_declarations_to_jax_docstrings.py -manifest [-output_file PATH]
      Write a JSON list with, for each matched JAX function, the Fortran
      routine, its Description and References text, the declaration
      section, the current docstring and any matching problems.
  -modules NAME[,NAME] limits every mode to those JAX module stems, e.g.
  -modules new_pdf,advance_windm_edsclrm_module.

Understanding:
  Matching.  A JAX module is paired with the Fortran file of the same stem.
  A function matches the routine of the same name, ignoring case, or a
  generic interface of that name. For an interface, the specific routine
  whose dummy arguments share the most names with the JAX arguments is
  used; ties go to the routine with fewer arguments that the JAX function
  lacks, then to the higher array rank. SPECIFIC_ROUTINE names the routine
  for interfaces whose specifics still tie. Names found in no paired file
  fall back to a unique match among all Fortran files. Private helpers
  (leading underscore) and nested functions are ignored.

  Description.  The comment block after a "Description:" label, or the
  first comment block after the routine header when there is no label.
  References come from a "References:" label. Both have "!" removed.

  Declaration section.  The declarations between `implicit none` and the
  first executable statement are classified by intent. Every argument
  declaration is copied with the comment lines above it, including the
  section dividers; local variable declarations, local-variable labels and
  compiler directives are left out. The Fortran text keeps its relative
  indentation and is placed at the end of the docstring after the line

      Fortran declarations (<routine> in <path>):

  -write_declaration replaces everything from that line to the end of the docstring,
  keeps the text above it, and puts the marker on the line after the
  opening quotes when there is no other text. Routines without arguments
  get no section.
  A docstring that contains a backslash is written as a raw string.

  Routines need dashed section dividers (see the Fortran formatting guide)
  for the copied block to be labelled; routines without any are flagged.
"""

import argparse
import ast
import inspect
import io
import json
import re
import sys
import tokenize
from pathlib import Path

if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from input.scripts.namelist import load_json

CLUBB_ROOT = Path(__file__).resolve().parents[2]
FORTRAN_DIR = CLUBB_ROOT / "src" / "CLUBB_core"
JAX_DIR = CLUBB_ROOT / "clubb_jax" / "src" / "CLUBB_core"
DESCRIPTION_EDITS = Path(__file__).resolve().parent / "custom_ports_of_fortran_descriptions.jsonc"

# Which specific routine to document for a generic interface that the
# argument names cannot resolve.
#
# A Fortran generic interface (e.g. smooth_min) groups several specific
# routines; the compiler picks one at each call from the types and ranks of
# the actual arguments. For the interfaces below, the specifics that take
# nz and ngrdcol have the same argument names and differ only in which
# inputs are arrays: smooth_min_arrays takes two (ngrdcol, nz) arrays,
# smooth_min_array_scalar an array and a scalar, and smooth_min_sclr_idx a
# scalar and an array. The single JAX function covers all of them through
# broadcasting, so no one specific is "the" counterpart. Document the most
# general one, which takes arrays for every input.
SPECIFIC_ROUTINE = {
    "smooth_min": "smooth_min_arrays",
    "smooth_max": "smooth_max_arrays",
}

SECTION_MARKER = "Fortran declarations ("
SECTION_LINE = re.compile(r"^\s*Fortran declarations \(")

# ---------------------------------------------------------------------------
# Fortran parsing
# ---------------------------------------------------------------------------

ROUTINE_START = re.compile(
    r"^\s*(?:(?:pure|elemental|recursive|impure|module)\s+)*"
    r"(?:(?:real|integer|logical|character|type|complex|double\s+precision)"
    r"\s*(?:\([^)]*\))?\s+)?"
    r"(subroutine|function)\s+(\w+)", re.I)
ROUTINE_END = re.compile(r"^\s*end\s*(subroutine|function)\b", re.I)
RESULT_CLAUSE = re.compile(r"result\s*\(\s*(\w+)\s*\)", re.I)
IMPLICIT_NONE = re.compile(r"^\s*implicit\s+none", re.I)
DECLARATION = re.compile(
    r"^\s*(type\s*\(|class\s*\(|(real|integer|logical|character|complex|"
    r"double\s+precision|dimension|save|data|intrinsic|external|procedure|"
    r"parameter)\b)", re.I)
INTENT = re.compile(r"intent\s*\(\s*(in\s*out|inout|in|out)\s*\)", re.I)
DIVIDER = re.compile(r"^\s*!\s*-{3,}\s*([A-Za-z/() ]+?)\s*-{3,}\s*$")
INTERFACE_START = re.compile(r"^\s*interface\s+(\w+)\s*$", re.I)
INTERFACE_END = re.compile(r"^\s*end\s+interface", re.I)
MODULE_PROCEDURE = re.compile(r"^\s*module\s+procedure\s+(.*)$", re.I)
LOCAL_LABEL = re.compile(r"^\s*!.*\blocal\b", re.I)
DIMENSION = re.compile(r"dimension\s*\(([^)]*)\)", re.I)


def code_part(line):
    """Return a Fortran line without its trailing comment."""
    return line.split("!", 1)[0]


def is_continued(line):
    return code_part(line).rstrip().endswith("&")


def is_comment(line):
    return line.strip().startswith("!")


def is_directive(line):
    """True for OpenACC/OpenMP sentinels and preprocessor lines."""
    stripped = line.strip().lower()
    return stripped.startswith("!$") or stripped.startswith("#")


def parse_fortran_file(path):
    """Return (lines, routines, interfaces) for one Fortran source file.

    Each routine is a dict with its name, kind, line range, dummy argument
    names, result name and the analysed declaration statements. Interfaces
    map a lower-case generic name to its lower-case specific routine names.
    """
    lines = path.read_text().splitlines()
    routines, interfaces, open_routines = [], {}, []
    generic = None

    for index, line in enumerate(lines):

        if generic is None:
            match = INTERFACE_START.match(line)
            if match:
                generic = match.group(1).lower()
                interfaces[generic] = []
                continue
        else:
            match = MODULE_PROCEDURE.match(line)
            if match:
                interfaces[generic] += [
                    name.strip().lower()
                    for name in code_part(match.group(1)).replace("&", "").split(",")
                    if name.strip()]
            if INTERFACE_END.match(line):
                generic = None
            continue

        if is_comment(line):
            continue

        if ROUTINE_END.match(line):
            if open_routines:
                routine = open_routines.pop()
                routine["end"] = index
                routines.append(routine)
            continue

        match = ROUTINE_START.match(line)
        if match and not re.match(r"^\s*end\b", line, re.I):
            header_end, header = index, code_part(line)
            while is_continued(lines[header_end]):
                header_end += 1
                while is_directive(lines[header_end]) or is_comment(lines[header_end]):
                    header_end += 1
                header += code_part(lines[header_end])
            kind, name = match.group(1).lower(), match.group(2)
            arguments = re.search(r"\(([^)]*)\)", header[match.end():])
            result = RESULT_CLAUSE.search(header)
            open_routines.append({
                "name": name,
                "kind": kind,
                "start": index,
                "header_end": header_end,
                "dummies": [argument.strip().lower()
                            for argument in arguments.group(1).replace("&", " ").split(",")
                            if argument.strip()] if arguments else [],
                "result": result.group(1) if result else (name if kind == "function" else None),
            })

    for routine in routines:
        analyse_declarations(lines, routine)

    return lines, routines, interfaces


def analyse_declarations(lines, routine):
    """Classify the declaration statements of one routine.

    Adds "implicit" (line of implicit none or None), "exec_start" (first
    executable line) and "statements": dicts with the first and last line
    and the intent class (in, inout, out, opt_* , result, nointent for an
    argument without intent, local, or mixed).
    """
    routine["statements"] = []
    routine["implicit"] = next(
        (k for k in range(routine["start"], routine["end"]) if IMPLICIT_NONE.match(lines[k])),
        None)
    if routine["implicit"] is None:
        routine["exec_start"] = routine["header_end"] + 1
        return

    result = (routine["result"] or "").lower()
    k = routine["implicit"] + 1
    while k < routine["end"]:
        line = lines[k]
        if not line.strip() or is_comment(line) or is_directive(line):
            k += 1
            continue
        if not DECLARATION.match(line):
            break

        first, text = k, [line]
        while is_continued(lines[k]):
            k += 1
            while is_comment(lines[k]) or is_directive(lines[k]) or not lines[k].strip():
                k += 1
            text.append(lines[k])

        joined = "\n".join(code_part(part) for part in text)
        attributes, _, entities = joined.partition("::")
        names = declared_names(entities)

        intent = INTENT.search(attributes)
        if intent:
            intent_class = intent.group(1).lower().replace(" ", "")
            if "optional" in attributes.lower():
                intent_class = "opt_" + intent_class
        elif result and result in names:
            intent_class = "result"
        elif names and all(name in routine["dummies"] for name in names):
            intent_class = "nointent"
        elif any(name in routine["dummies"] for name in names):
            intent_class = "mixed"
        else:
            intent_class = "local"

        rank = max((dims.count(",") + 1 for dims in DIMENSION.findall(attributes)), default=0)
        routine["statements"].append(
            {"first": first, "last": k, "intent": intent_class, "rank": rank})
        k += 1

    routine["exec_start"] = k


def declared_names(entities):
    """Return the lower-case names declared after "::" in a statement.

    Array bounds and initializers are removed first, so the comma in
    "alfa(i1:i2,km)" does not make "km" look like a declared name.
    """
    text = entities.replace("\n", " ").replace("&", " ")
    depth, kept = 0, []
    for char in text:
        if char in "([":
            depth += 1
        elif char in ")]":
            depth -= 1
        elif depth == 0:
            kept.append(char)
    return [re.match(r"\s*(\w+)", entity).group(1).lower()
            for entity in "".join(kept).split(",") if re.match(r"\s*\w", entity)]


def declaration_section(lines, routine):
    """Return the argument declarations of a routine as text lines.

    Each argument statement is copied with the comment lines between the
    previous statement and itself. Local statements, local-variable labels
    and directives are dropped, runs of blank lines are collapsed, and the
    common indentation is removed.
    """
    section, previous_end = [], routine["implicit"]
    for statement in routine["statements"]:
        if statement["intent"] != "local":
            for k in range(previous_end + 1, statement["last"] + 1):
                line = lines[k].rstrip()
                if k < statement["first"] and (is_directive(line) or LOCAL_LABEL.match(line)):
                    continue
                if is_comment(line) and line.strip().lower().startswith("!$"):
                    continue
                section.append(line)
        previous_end = statement["last"]

    collapsed = []
    for line in section:
        if not line.strip() and (not collapsed or not collapsed[-1].strip()):
            continue
        collapsed.append(line)
    while collapsed and not collapsed[-1].strip():
        collapsed.pop()

    indent = min((len(line) - len(line.lstrip()) for line in collapsed if line.strip()), default=0)
    return [line[indent:] for line in collapsed]


def indentation(line):
    return len(line) - len(line.lstrip())


def dedent_comment(body):
    """Remove the indentation of comment text taken from after "!".

    Lines lose their common indentation. When the first line is still
    indented, because some line follows "!" directly (e.g. "!k=1 |" rows of
    an ASCII matrix), each paragraph loses up to that much more, but no
    more than its own smallest indentation, so diagrams stay aligned.
    """
    common = min((indentation(line) for line in body if line.strip()), default=0)
    body = [line[common:] for line in body]
    first = indentation(next((line for line in body if line.strip()), ""))
    if not first:
        return body

    result, paragraph = [], []
    for line in body + [""]:
        if line.strip():
            paragraph.append(line)
            continue
        if paragraph:
            shift = min(first, min(indentation(l) for l in paragraph))
            result += [l[shift:] for l in paragraph]
            paragraph = []
        result.append(line)
    return result[:-1]


def ends_comment_block(lines, k):
    """True if line k is a divider or a dash rule that closes a comment block.

    A rule written "! ---" directly below a line of text of about the same
    length underlines a heading inside the description and does not close it.
    """
    if DIVIDER.match(lines[k]):
        return True
    # Closing rules start at most one space after "!"; indented dash lines
    # are part of ASCII diagrams.
    if not re.match(r"^\s*! ?-{5,}\s*$", lines[k]):
        return False
    above = lines[k - 1].strip().lstrip("!").strip() if k > 0 else ""
    rule = lines[k].strip().lstrip("!").strip()
    underline = lines[k].strip().startswith("! ") and is_comment(lines[k - 1]) and \
        above and set(above) != {"-"} and abs(len(rule) - len(above)) <= 2
    return not underline


def comment_block(lines, start, keyword):
    """Return the comment text that follows a "<keyword>:" line, or None.

    The block ends at a divider or closing dash rule, the other heading
    ("Description:" or "Reference(s):") or the first line that is not a
    comment.
    """
    for k in range(start, len(lines)):
        if not is_comment(lines[k]):
            if lines[k].strip() and not lines[k].strip().lower().startswith("use"):
                return None
            continue
        text = lines[k].strip().lstrip("!")
        if not re.match(rf"\s*{keyword}\s*:", text, re.I):
            continue

        body = []
        inline = re.sub(rf"^\s*{keyword}\s*:", "", text, flags=re.I).strip()
        if inline:
            body.append(inline)
        for j in range(k + 1, len(lines)):
            line = lines[j]
            if not is_comment(line) or ends_comment_block(lines, j):
                break
            content = line.strip()[1:]
            if re.match(r"\s*(description|references?)\s*:", content, re.I):
                break
            body.append(content.rstrip())
        while body and not body[-1].strip():
            body.pop()
        while body and not body[0].strip():
            body.pop(0)
        return dedent_comment(body)
    return None


def unlabelled_description(lines, start):
    """Return the first comment block after a routine header, or None.

    Used when a routine describes itself without a "Description:" label.
    The block ends at a "References:" line, a divider or a non-comment line.
    """
    k = start
    while k < len(lines) and not lines[k].strip():
        k += 1
    body = []
    while k < len(lines) and (is_comment(lines[k]) or not lines[k].strip()):
        if ends_comment_block(lines, k) or re.match(r"^\s*!\s*references?\s*:", lines[k], re.I):
            break
        body.append(lines[k].strip()[1:].rstrip() if lines[k].strip() else "")
        k += 1
    while body and not body[-1].strip():
        body.pop()
    if not body:
        return None
    return dedent_comment(body)


def routine_comments(lines, routine):
    """Return (description, references, labelled) for a routine.

    labelled is False when the description had no "Description:" label and
    was taken from the first comment block after the header.
    """
    start = routine["header_end"] + 1
    stop = routine["implicit"] if routine["implicit"] is not None else routine["end"]
    window = lines[:stop]
    description = comment_block(window, start, "Description")
    labelled = description is not None
    if not labelled:
        description = unlabelled_description(window, start)
    return description, comment_block(window, start, "References?"), labelled


# ---------------------------------------------------------------------------
# Matching JAX functions with Fortran routines
# ---------------------------------------------------------------------------

class FortranIndex:
    """All routines and generic interfaces in one Fortran directory."""

    def __init__(self, fortran_dir):
        self.files = {}
        for path in sorted(fortran_dir.glob("*.F90")):
            lines, routines, interfaces = parse_fortran_file(path)
            self.files[path.stem] = {
                "path": path,
                "lines": lines,
                "routines": {routine["name"].lower(): routine for routine in routines},
                "interfaces": interfaces,
            }

    def candidates(self, stem, name):
        """Return [(file stem, routine)] matching a name within one file."""
        entry = self.files.get(stem)
        if entry is None:
            return []
        if name in entry["routines"]:
            return [(stem, entry["routines"][name])]
        return [(stem, entry["routines"][specific])
                for specific in entry["interfaces"].get(name, [])
                if specific in entry["routines"]]

    def find(self, stem, name, jax_arguments):
        """Return (file stem, routine, problems) for a JAX function name."""
        problems = []
        found = self.candidates(stem, name)
        if not found:
            matches = [other for other in self.files if self.candidates(other, name)]
            if len(matches) == 1:
                found = self.candidates(matches[0], name)
                if stem in self.files:
                    problems.append(f"matched in {matches[0]}.F90, not the paired file")
            elif len(matches) > 1:
                return None, None, [f"name found in several files: {', '.join(matches)}"]
        if not found:
            return None, None, []

        if len(found) > 1:
            def score(item):
                routine = item[1]
                shared = len(set(jax_arguments) & set(routine["dummies"]))
                extra = len(set(routine["dummies"]) - set(jax_arguments))
                rank = max((s["rank"] for s in routine["statements"]), default=0)
                return (shared, -extra, rank)
            found.sort(key=score, reverse=True)
            chosen = SPECIFIC_ROUTINE.get(name)
            if chosen:
                found.sort(key=lambda item: item[1]["name"].lower() != chosen)
            elif score(found[0]) == score(found[1]):
                problems.append("generic interface: no single best specific routine")
            problems.append(
                f"generic interface {name}: used {found[0][1]['name']} out of "
                f"{', '.join(item[1]['name'] for item in found)}")
        return found[0][0], found[0][1], problems


def jax_functions(path):
    """Return the module-level public functions of a JAX module."""
    tree = ast.parse(path.read_text())
    return [node for node in tree.body
            if isinstance(node, ast.FunctionDef) and not node.name.startswith("_")]


def own_returns(node):
    """Return the return statements of a function, excluding nested ones."""
    found, pending = [], list(node.body)
    while pending:
        child = pending.pop()
        if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda, ast.ClassDef)):
            continue
        if isinstance(child, ast.Return):
            found.append(child)
        pending.extend(ast.iter_child_nodes(child))
    return found


def returned_values(node, routine):
    """Describe what a JAX function returns, or None when it varies.

    Returns a list of names for "return a, b", ["dict with keys a, b"] for a
    dict literal, and the Fortran result name for a function that returns a
    single expression.
    """
    shapes = set()
    for statement in own_returns(node):
        value = statement.value
        if isinstance(value, ast.Name):
            shapes.add((value.id,))
        elif isinstance(value, ast.Tuple) and all(isinstance(e, ast.Name) for e in value.elts):
            shapes.add(tuple(e.id for e in value.elts))
        elif isinstance(value, ast.Dict) and value.keys and \
                all(isinstance(k, ast.Constant) and isinstance(k.value, str) for k in value.keys):
            shapes.add(("dict with keys " + ", ".join(k.value for k in value.keys),))
        elif value is not None and not isinstance(value, ast.Tuple) and routine["kind"] == "function":
            shapes.add((routine["result"],))
        else:
            return None
    return list(shapes.pop()) if len(shapes) == 1 else None


def function_arguments(node):
    args = node.args
    return [arg.arg.lower() for arg in args.posonlyargs + args.args + args.kwonlyargs]


# ---------------------------------------------------------------------------
# Docstring editing
# ---------------------------------------------------------------------------

def docstring_node(node):
    if node.body and isinstance(node.body[0], ast.Expr) and \
       isinstance(node.body[0].value, ast.Constant) and isinstance(node.body[0].value.value, str):
        return node.body[0]
    return None


def display_path(path):
    """Return a path relative to the CLUBB checkout when it lies inside it."""
    try:
        return str(path.relative_to(CLUBB_ROOT))
    except ValueError:
        return str(path)


def section_text(stem, routine, declarations, fortran_dir):
    path = display_path(fortran_dir / f"{stem}.F90")
    return [f"{SECTION_MARKER}{routine['name']} in {path}):", ""] + \
           ["    " + line if line else "" for line in declarations]


def current_section(docstring):
    """Return the declaration-section lines of a docstring, or None."""
    lines = inspect.cleandoc(docstring).splitlines()
    for k, line in enumerate(lines):
        if SECTION_LINE.match(line):
            return [l.rstrip() for l in lines[k:]]
    return None


def header_end_line(source_lines, node):
    """Return the 0-based line holding the colon that ends a def header."""
    start = node.lineno - 1
    text = "\n".join(source_lines[start:node.body[0].lineno]) + "\n"
    depth = 0
    for token in tokenize.generate_tokens(io.StringIO(text).readline):
        if token.type == tokenize.OP:
            if token.string in "([{":
                depth += 1
            elif token.string in ")]}":
                depth -= 1
            elif token.string == ":" and depth == 0:
                return start + token.start[0] - 1
    raise ValueError(f"cannot find the end of the header of {node.name}")


def rendered_docstring(prose_lines, section, indent):
    """Return the docstring source lines for prose followed by the section."""
    body = list(prose_lines)
    while body and not body[-1].strip():
        body.pop()
    # Without prose the marker starts on its own line so that cleandoc keeps
    # the indentation of the declarations relative to it.
    body.append("")
    body += section
    prefix = "r" if any("\\" in line for line in body) else ""
    lines = [f'{indent}{prefix}"""{body[0]}'] + \
            [f"{indent}{line}" if line else "" for line in body[1:]]
    lines.append(f'{indent}"""')
    return lines


def prose_of(docstring):
    """Return the docstring lines before the declaration section, dedented."""
    lines = inspect.cleandoc(docstring).splitlines()
    for k, line in enumerate(lines):
        if SECTION_LINE.match(line):
            return [l.rstrip() for l in lines[:k]]
    return [l.rstrip() for l in lines]


def strip_docstrings(tree):
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef, ast.Module)):
            if docstring_node(node) is not None:
                node.body = node.body[1:]
    return ast.dump(tree)


def rewrite_module(path, sections, prose=None):
    """Write the declaration sections, and optionally the prose, of one module.

    sections maps function names to section lines; prose maps function
    names to new prose lines (the current prose is kept for the others).
    Returns the names whose docstring changed.
    """
    source = path.read_text()
    source_lines = source.splitlines()
    changed, edits = [], []

    for node in jax_functions(path):
        if node.name not in sections:
            continue
        doc = docstring_node(node)
        indent = " " * node.body[0].col_offset
        text = (prose or {}).get(node.name)
        if text is None:
            text = prose_of(doc.value.value) if doc else []
        new_lines = rendered_docstring(text, sections[node.name], indent)
        if doc:
            first, last = doc.lineno - 1, doc.end_lineno
        else:
            first = last = header_end_line(source_lines, node) + 1
        if source_lines[first:last] != new_lines:
            edits.append((first, last, new_lines))
            changed.append(node.name)

    for first, last, new_lines in sorted(edits, reverse=True):
        source_lines[first:last] = new_lines

    new_source = "\n".join(source_lines) + ("\n" if source.endswith("\n") else "")
    if strip_docstrings(ast.parse(source)) != strip_docstrings(ast.parse(new_source)):
        raise RuntimeError(f"{path}: rewriting docstrings would change code; not written")
    if changed:
        path.write_text(new_source)
    return changed


# ---------------------------------------------------------------------------
# Docstring prose from the Fortran comments
# ---------------------------------------------------------------------------

NO_REFERENCES = {"", "none", "none.", "(none)"}


EDIT_KEYS = {"drop", "replace", "append"}


def edited_text(lines, edit, name):
    """Apply the "drop", "replace" and "append" edits to prose lines."""
    for needle in edit.get("drop", []):
        if not any(needle in line for line in lines):
            raise ValueError(f"{name}: drop {needle!r} matches no line")
        lines = [line for line in lines if needle not in line]
    text = "\n".join(lines)
    for old, new in edit.get("replace", []):
        if old not in text:
            raise ValueError(f"{name}: replace {old!r} matches nothing")
        text = text.replace(old, new)

    result = []
    for line in text.split("\n"):
        if not line.strip() and (not result or not result[-1].strip()):
            continue
        result.append(line.rstrip())
    while result and not result[-1].strip():
        result.pop()
    if edit.get("append"):
        result += [""] + edit["append"].split("\n")
    return result


def fortran_prose(record, current, edit):
    """Return (prose lines, origin) for one function's docstring.

    The prose is the full Fortran description followed by its references,
    adjusted by the "drop", "replace" and "append" edits. A function without
    a Fortran description keeps its current prose.
    """
    unknown = set(edit) - EDIT_KEYS
    if unknown:
        raise ValueError(f"{record['function']}: unknown edit keys {', '.join(sorted(unknown))}")
    description = record["description"]
    if not description:
        return current, "kept: no Fortran description"
    lines, origin = list(description), "fortran"

    references = record["references"] or []
    if not all(line.strip().lower() in NO_REFERENCES for line in references):
        lines += ["", "References:"] + references
    if any(key in edit for key in ("drop", "replace", "append")):
        origin += " with edits"
    return edited_text(lines, edit, record["function"]) + returns_text(record["returns"]), origin


def returns_text(values):
    """Return a "Returns:" paragraph listing the returned values."""
    if not values:
        return []
    text = values[0] if len(values) == 1 else "(" + ", ".join(values) + ")"
    lines, line = [], "   "
    for word in text.split(" "):
        if len(line) + len(word) + 1 > 76:
            lines.append(line)
            line = "    " + " " * (text.startswith("(") or text.startswith("dict"))
        line += " " + word
    return ["", "Returns:"] + lines + [line]


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def collect(index, jax_dir, fortran_dir, modules):
    """Return one record per JAX function with a Fortran counterpart."""
    records = []
    for path in sorted(jax_dir.glob("*.py")):
        if path.stem == "__init__" or (modules and path.stem not in modules):
            continue
        for node in jax_functions(path):
            stem, routine, problems = index.find(path.stem, node.name.lower(),
                                                 function_arguments(node))
            if routine is None:
                if problems:
                    records.append({"jax_file": display_path(path),
                                    "function": node.name, "line": node.lineno,
                                    "problems": problems})
                continue
            lines = index.files[stem]["lines"]
            if routine["implicit"] is None:
                problems.append("no implicit none; declarations not found")
            elif not any(DIVIDER.match(lines[k])
                         for k in range(routine["implicit"], routine["exec_start"])):
                problems.append("no section dividers in the Fortran declarations")
            if any(s["intent"] == "mixed" for s in routine["statements"]):
                problems.append("a declaration mixes arguments and locals")
            description, references, labelled = routine_comments(lines, routine)
            if description is None:
                problems.append("no description comment")
            elif not labelled:
                problems.append("description has no Description: label")
            declarations = declaration_section(lines, routine) if routine["implicit"] is not None else []
            if routine["implicit"] is not None and not declarations:
                problems.append("no argument declarations; no section written")
            doc = docstring_node(node)
            records.append({
                "jax_file": display_path(path),
                "function": node.name,
                "line": node.lineno,
                "fortran_file": display_path(fortran_dir / f"{stem}.F90"),
                "routine": routine["name"],
                "fortran_line": routine["start"] + 1,
                "description": description,
                "references": references,
                "section": section_text(stem, routine, declarations, fortran_dir) if declarations else None,
                "docstring": doc.value.value if doc else None,
                "returns": returned_values(node, routine),
                "problems": problems,
            })
    return records


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="port_fortran_declarations_to_jax_docstrings.py",
        description="Check or write the Fortran declaration section of JAX docstrings.")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("-check", action="store_true",
                      help="Report missing or outdated sections (default)")
    mode.add_argument("-write_declaration", action="store_true",
                      help="Insert or replace the declaration sections in the JAX files")
    mode.add_argument("-manifest", action="store_true",
                      help="Write matched functions, Fortran comments and sections as JSON")
    mode.add_argument("-write_description", action="store_true",
                      help="Rewrite the docstring descriptions from the Fortran comments and "
                           "-edits, then the declaration sections")
    parser.add_argument("-edits", type=Path, default=DESCRIPTION_EDITS,
                        help="JSONC edits for -write_description, {module: {function: edit}} "
                             f"(default: {DESCRIPTION_EDITS.relative_to(CLUBB_ROOT)})")
    parser.add_argument("-output_file", type=Path,
                        help="Manifest destination (default: standard output)")
    parser.add_argument("-modules", default="",
                        help="Comma-separated JAX module stems to process (default: all)")
    parser.add_argument("-fortran_dir", type=Path, default=FORTRAN_DIR,
                        help=f"Fortran source directory (default: {FORTRAN_DIR.relative_to(CLUBB_ROOT)})")
    parser.add_argument("-jax_dir", type=Path, default=JAX_DIR,
                        help=f"JAX source directory (default: {JAX_DIR.relative_to(CLUBB_ROOT)})")
    args = parser.parse_args(argv)
    if args.output_file and not args.manifest:
        parser.error("-output_file is only used with -manifest")

    fortran_dir, jax_dir = args.fortran_dir.resolve(), args.jax_dir.resolve()
    modules = {stem for stem in args.modules.split(",") if stem}
    records = collect(FortranIndex(fortran_dir), jax_dir, fortran_dir, modules)
    matched = [record for record in records if "routine" in record]

    if args.manifest:
        text = json.dumps(records, indent=1) + "\n"
        if args.output_file:
            args.output_file.write_text(text)
            print(f"Wrote {len(matched)} matched functions to {args.output_file}")
        else:
            sys.stdout.write(text)
        return 0

    for record in records:
        for problem in record["problems"]:
            print(f"note: {record['jax_file']}:{record['line']} {record['function']}: {problem}")

    if args.write_description:
        edits = load_json(args.edits)
        by_file = {}
        for record in matched:
            if record["section"] is not None:
                by_file.setdefault(record["jax_file"], []).append(record)
        known = {(Path(f).stem, r["function"]) for f, rs in by_file.items() for r in rs}
        unknown = [f"{m}.{f}" for m, fs in edits.items() for f in fs
                   if (m, f) not in known and (not modules or m in modules)]
        if unknown:
            parser.error(f"edits name unmatched functions: {', '.join(sorted(unknown))}")
        total = 0
        for jax_file, file_records in sorted(by_file.items()):
            module_edits = edits.get(Path(jax_file).stem, {})
            sections, prose = {}, {}
            for record in file_records:
                current = prose_of(record["docstring"]) if record["docstring"] else []
                try:
                    text, origin = fortran_prose(record, current,
                                                 module_edits.get(record["function"], {}))
                except ValueError as error:
                    parser.error(f"{jax_file}: {error}; update {args.edits}")
                if origin.startswith("kept:"):
                    print(f"note: {jax_file}:{record['line']} {record['function']}: {origin[6:]}")
                sections[record["function"]], prose[record["function"]] = record["section"], text
            changed = rewrite_module(CLUBB_ROOT / jax_file, sections, prose)
            total += len(changed)
            for name in changed:
                print(f"updated: {jax_file} {name}")
        print(f"{total} docstrings updated, {len(matched)} matched functions")
        return 0

    if args.write_declaration:
        by_file = {}
        for record in matched:
            if record["section"] is None:
                continue
            by_file.setdefault(record["jax_file"], {})[record["function"]] = record["section"]
        total = 0
        for jax_file, sections in sorted(by_file.items()):
            changed = rewrite_module(CLUBB_ROOT / jax_file, sections)
            total += len(changed)
            for name in changed:
                print(f"updated: {jax_file} {name}")
        print(f"{total} docstrings updated, {len(matched)} matched functions")
        return 0

    outdated = 0
    matched = [record for record in matched if record["section"] is not None]
    for record in matched:
        current = current_section(record["docstring"] or "")
        if current != record["section"]:
            outdated += 1
            state = "missing" if current is None else "differs from Fortran"
            print(f"{record['jax_file']}:{record['line']} {record['function']}: "
                  f"declaration section {state} ({record['fortran_file']} {record['routine']})")
    print(f"{outdated} of {len(matched)} matched functions need updating")
    return 1 if outdated else 0


if __name__ == "__main__":
    sys.exit(main())
