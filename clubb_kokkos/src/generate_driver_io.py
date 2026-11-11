#!/usr/bin/env python3
"""Generate permanent C++ bindings for the standalone namelist groups.

Run this file to regenerate driver_io.hpp after a Native namelist changes.
driver_groups reads the declarations and group membership in clubb_driver,
parameters_tunable and inputfields. render_header is also used by the
microphysics binding generator. Runtime parsing belongs to namelist::Group;
defaults, initialization order and error decisions belong to the callers.

Shared generation utilities module_declarations/type_declarations resolve
Native kinds and record fields; cpp_argument retains the existing C++ formal
types. Binding and Namelist hold this source metadata without runtime state.

The declaration reader is shared with the existing source generators. This
generator emits no Fortran wrappers, field buffers or mutable model storage.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
import re
import sys
import textwrap

SOURCE = Path(__file__).resolve().parent
ROOT = SOURCE.parents[1]
sys.path.insert(0, str(SOURCE / "f2k"))
from fortran_procedures import Argument, declaration, source_files, statements
from generate_derived_types import split_top, write_generated


@dataclass
class Binding:
    """One Native namelist member/component and its existing C++ argument."""

    name: str
    expression: str
    argument: Argument
    extent: str = ""


@dataclass
class Namelist:
    """A namelist reader's current argument contract and Native binding order."""

    name: str
    source: str
    arguments: list[Argument]
    bindings: list[Binding]
    filename: str
    unused: tuple[str, ...] = ()
    reader_name: str = ""


def module_declarations() -> dict[tuple[str, str], Argument]:
    """Resolve declarations only; USE aliases are handled by each source owner."""
    result = {}
    for _, text in source_files():
        match = re.search(r"(?im)^\s*module\s+(?!procedure)(\w+)", text)
        if not match:
            continue
        body = re.split(r"(?im)^\s*contains\s*$", text, maxsplit=1)[0]
        for statement in statements(body):
            try:
                fields = declaration(statement)
            except ValueError:
                continue
            for name, argument in fields.items():
                result.setdefault((match[1].lower(), name), argument)
    return result


def type_declarations(path: Path, name: str) -> dict[str, Argument]:
    """Read a small namelist record in its Native declaration order."""
    body = re.search(
        r"(?mis)^\s*type\s+" + name + r"\b(.*?)^\s*end\s*type",
        path.read_text(),
    )[1]
    result = {}
    for statement in statements(body):
        result.update(declaration(statement))
    return result


def cpp_argument(argument: Argument, array_type: str | None) -> str:
    """Retain the current generated C++ argument type and ownership contract."""
    if array_type:
        return f"{array_type} {argument.name}"
    bases = {
        "integer": "int",
        "logical": "bool",
        "real": "clubb_kokkos::real_t",
        "double": "double",
        "float": "float",
        "character": "std::string",
        "sponge_damp_settings": "clubb_kokkos::sponge_damp_settings",
    }
    value = bases[argument.kind]
    if argument.intent != "in":
        return f"{value}& {argument.name}"
    if argument.kind == "character":
        return f"const {value}& {argument.name}"
    return f"{value} {argument.name}"


def render_header(description: list[str], groups: list[Namelist], includes: list[str]) -> str:
    """Emit only typed bindings, with Native names and fixed character widths."""
    lines = description + ["#pragma once", '#include "namelist.hpp"']
    lines += [f'#include "{name}"' for name in includes]
    lines += [
        "#include <stdexcept>",
        "#include <string>",
        "#include <utility>",
        "",
        "namespace clubb_kokkos {",
        "",
    ]
    for group in groups:
        templates = {
            a.name: f"Array{i}" for i, a in enumerate(a for a in group.arguments if a.dims)
        }
        if templates:
            template = "template <" + ", ".join("class " + t for t in templates.values()) + ">"
            lines += textwrap.wrap(template, width=96, subsequent_indent="          ")
        lines += [f"inline void {group.reader_name or 'read_' + group.name}("]
        for index, argument in enumerate(group.arguments):
            last = index + 1 == len(group.arguments)
            intent = {"in": "In", "out": "Out", "inout": "InOut"}[argument.intent]
            lines += [
                "    "
                + cpp_argument(argument, templates.get(argument.name))
                + (" // " if last else ", // ")
                + intent
            ]
        lines += [
            ") {",
            "  // Description:",
            f"  // Read &{group.name} from initialized caller values; omitted entries keep them.",
            f"  // Native declaration: {group.source}.",
        ]
        if group.unused:
            lines += [
                "  // The C++ reader owns a call-local file stream; Native unit arguments remain in the API."
            ]
            lines += [f"  (void){name};" for name in group.unused]
        lines += [f'  namelist::Group group("{group.name}");']
        for binding in group.bindings:
            expression = binding.expression
            if binding.extent:
                expression = f"Kokkos::subview({expression}, std::make_pair(0, {binding.extent}))"
            width = (
                ", " + binding.argument.character_length
                if binding.argument.kind == "character"
                else ""
            )
            line = f'  group.bind("{binding.name}", {expression}{width});'
            if len(line) <= 100:
                lines.append(line)
            else:
                lines += [f'  group.bind("{binding.name}",', f"             {expression}{width});"]
        if any(a.name == "read_status" for a in group.arguments):
            lines += [f"  read_status = group.read_file({group.filename});"]
        else:
            lines += [
                f"  if (group.read_file({group.filename}) != 0)",
                "    throw std::runtime_error(group.error());",
            ]
        lines += ["}", ""]
    return "\n".join(lines + ["} // namespace clubb_kokkos", ""])


def driver_groups() -> list[Namelist]:
    """Outline the six existing driver bindings from their Native declarations."""
    source = "src/clubb_driver.F90"
    code = statements((ROOT / source).read_text())
    start = next(
        i for i, s in enumerate(code) if s.lower().startswith("subroutine init_clubb_case(")
    )
    end = next(
        i
        for i in range(start, len(code))
        if code[i].lower().startswith("end subroutine init_clubb_case")
    )
    declarations = {}
    for statement in code[:end]:
        declarations.update(declaration(statement))
    modules = module_declarations()
    for module, names in {
        "time_dependent_input": "l_t_dependent l_input_xpwp_sfc l_ignore_forcings",
        "soil_vegetation": "l_soil_veg",
    }.items():
        for name in names.split():
            declarations[name] = modules[module, name]
    sponge = type_declarations(
        ROOT / "src/CLUBB_core/sponge_layer_damping.F90", "sponge_damp_settings"
    )
    result = []
    for name in ["multicol_def", "model_setting", "stats_setting", "configurable_clubb_flags_nl"]:
        statement = next(
            s for s in code[start:end] if re.match(r"namelist\s*/" + name + "/", s, re.I)
        )
        members = split_top(statement.split("/", 2)[2])
        arguments = []
        records = []
        bindings = []
        for member in members:
            if "sponge_damp_settings" in member.lower():
                records.append(Argument(member, "sponge_damp_settings", (), "inout", False, False))
                bindings += [
                    Binding(member + "%" + field.name, member + "." + field.name, field)
                    for field in sponge.values()
                ]
            else:
                argument = replace(declarations[member.lower()], intent="inout")
                arguments.append(argument)
                bindings.append(
                    Binding(
                        member,
                        argument.name,
                        argument,
                        "sclr_max" if member.lower() == "sclr_tol_nl" else "",
                    )
                )
        arguments = (
            [Argument("runfile", "character", (), "in", False, False, "*")] + arguments + records
        )
        arguments += [Argument("read_status", "integer", (), "out", False, False)]
        result.append(Namelist(name, source, arguments, bindings, "runfile"))

    source = "src/CLUBB_core/parameters_tunable.F90"
    text = re.search(
        r"(?ims)^  subroutine init_clubb_params_api\b.*?^  end subroutine init_clubb_params_api\b",
        (ROOT / source).read_text(),
    )[0]
    declarations = {}
    for statement in statements(text):
        declarations.update(declaration(statement))
    statement = next(
        s for s in statements(text) if re.match(r"namelist\s*/clubb_params_nl/", s, re.I)
    )
    members = split_top(statement.split("/", 2)[2])
    arguments = [
        Argument("runfile", "character", (), "in", False, False, "*"),
        replace(declarations["ngrdcol"], intent="in"),
        replace(declarations["iunit"], intent="in"),
    ]
    parameters = [replace(declarations[name.lower()], intent="inout") for name in members]
    arguments += parameters + [Argument("read_status", "integer", (), "out", False, False)]
    result.append(
        Namelist(
            "clubb_params_nl",
            source,
            arguments,
            [Binding(name, a.name, a, "ngrdcol") for name, a in zip(members, parameters)],
            "runfile",
            ("iunit",),
        )
    )

    source = "src/Input_fields/input_fields.F90"
    full = (ROOT / source).read_text()
    body = re.search(
        r"(?ims)^  subroutine inputfields_init\b.*?^  end subroutine inputfields_init", full
    )[0]
    code = statements(body)
    declarations = {}
    for statement in statements(full.split("  contains")[0]) + code:
        declarations.update(declaration(statement))
    statement = next(s for s in code if s.lower().startswith("namelist /setfields/"))
    members = split_top(statement.split("/", 2)[2])
    parameters = [
        replace(declarations[name.lower()], name=name, intent="inout") for name in members
    ]
    arguments = [
        Argument("iunit", "integer", (), "in", False, False),
        Argument("namelist_filename", "character", (), "in", False, False, "*"),
    ] + parameters
    result.append(
        Namelist(
            "setfields",
            source,
            arguments,
            [Binding(name, a.name, a) for name, a in zip(members, parameters)],
            "namelist_filename",
            ("iunit",),
            "read_setfields_nl",
        )
    )
    return result


def main() -> None:
    """Regenerate the permanent six-group header without changing other owners."""
    groups = driver_groups()
    text = render_header(
        [
            "// Read standalone settings into caller-owned C++ values.",
            "//",
            "// Public readers: read_multicol_def, read_model_setting, read_stats_setting,",
            "// read_configurable_clubb_flags_nl, read_clubb_params_nl and read_setfields_nl.",
            "// Call after supplying the Native defaults. Bindings retain Native member",
            "// names, character widths and explicit array bounds; namelist::Group owns",
            "// parsing and preserves unassigned/null entries. No case state is retained.",
            "// Generated by generate_driver_io.py from the Native declarations.",
        ],
        groups,
        ["CLUBB_core/sponge_layer_damping.hpp", "Input_fields/sounding.hpp"],
    )
    write_generated(SOURCE / "driver_io.hpp", text)
    print("Permanent driver namelist bindings: 6")


if __name__ == "__main__":
    main()
