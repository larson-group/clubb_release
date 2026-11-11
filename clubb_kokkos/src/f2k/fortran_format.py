"""Keep the three generated f2k modules readable with shared descriptor comments.

This handles the declaration forms emitted by the bridge generators, rather
than arbitrary Fortran. Native signatures and conversion statements stay in
their generators; this pass groups declarations and separates code stages.
"""

from __future__ import annotations

import re
import textwrap

DESCRIPTIONS = {
    "input": "Caller-owned input field descriptors [-]",
    "output": "Returned descriptors; C++ releases their C-owned buffers [-]",
    "field": "Descriptor of one packed field [-]",
    "rank": "Number of array dimensions [-]",
    "extent": "Array extents in native dimension order [-]",
    "nbytes": "Output buffer size [bytes]",
    "char_len": "Fixed character width; zero for numeric fields [characters]",
    "ip": "Borrowed view of packed integer data [-]",
    "bp": "Borrowed view of packed logical bytes [-]",
    "rp": "Borrowed view of packed real data; units follow the native field",
    "dp": "Borrowed view of packed double-precision data; units follow the native field",
    "cp": "Borrowed view of packed character bytes [-]",
    "character_buffer": "Borrowed view of packed character bytes [-]",
    "n": "Number of packed elements [-]",
    "j": "Packed array or field index [-]",
    "k": "Character or secondary array index [-]",
    "h": "Hydrometeor species index [-]",
    "ncol": "Number of columns [-]",
    "nzt": "Number of vertical levels [-]",
    "cells": "Number of column-level cells [-]",
    "offset": "Offset of one species component in the packed array [-]",
    "extents": "Output extents in native dimension order [-]",
    "field_shape": "Expected field extents in native dimension order [-]",
    "unit": "Call-local namelist file unit [-]",
    "name": "Explicit-length C setting name [-]",
    "text": "Native setting name reconstructed from C bytes [-]",
    "private_host_lat": "Private error-context column latitudes [degrees north]",
    "private_host_lon": "Private error-context column longitudes [degrees east]",
    "private_chunk_idx": "Private error-context chunk index [-]",
    "private_mpi_rank": "Private error-context MPI rank [-]",
}
SCRATCH = {
    "ip",
    "bp",
    "rp",
    "dp",
    "cp",
    "character_buffer",
    "n",
    "j",
    "k",
    "extents",
    "field_shape",
}


def code_part(line: str) -> str:
    """Strip a Fortran comment while preserving quoted exclamation marks."""
    return re.split(r"!", re.sub(r"'[^']*'|\"[^\"]*\"", "''", line), maxsplit=1)[0]


def entities(text: str) -> list[str]:
    """Split emitted declaration entities without splitting array dimensions."""
    return re.split(r",(?![^()]*\))", text)


def spaced_body(lines: list[str]) -> list[str]:
    """Separate control flow and complete continued statements as visual units."""
    result = []
    continued = []
    for line in lines:
        continued.append(line.rstrip())
        if code_part(line).rstrip().endswith("&"):
            continue

        first = continued[0].strip()
        control = bool(
            re.match(
                r"(?:if\b|else\b|end\s+(?:if|do|select)\b|do\b|select\b|case\b)",
                first,
                re.I,
            )
        )
        separate = control or len(continued) > 1 or first.startswith("!")
        if separate and result and result[-1]:
            result.append("")
        result.extend(continued)
        if control or len(continued) > 1:
            result.append("")
        continued = []
    if continued:
        raise ValueError("unfinished generated Fortran continuation")
    return result


def format_routine(lines: list[str]) -> list[str]:
    """Document emitted descriptors/temporaries and remove unused scratch locals."""
    start = next((i for i, line in enumerate(lines) if "::" in line), None)
    if start is None:
        return lines
    end = start
    while end < len(lines) - 1:
        line = lines[end].strip()
        if (
            "::" in line
            or line.startswith("namelist ")
            or not line
            or line.startswith("!")
        ):
            end += 1
            continue
        break
    # Namelist and USE continuations are single list entries before joining.
    declarations, namelists = [], []
    # Comments introducing the first converted field belong to the body.
    while end > start and (
        not lines[end - 1].strip() or lines[end - 1].lstrip().startswith("!")
    ):
        end -= 1
    body = lines[end:-1]
    used = "\n".join(code_part(line) for line in body)
    for line in lines[start:end]:
        if line.strip().startswith("namelist "):
            namelists.append(line)
            continue
        if "::" not in line:
            continue
        spec, values = line.split("::", 1)
        for value in entities(values):
            value, _, comment = value.partition("!")
            value = value.strip()
            name = re.match(r"\w+", value).group()
            if name in SCRATCH and not re.search(rf"\b{name}\b", used, re.I):
                continue
            intent = re.search(r"intent\((inout|in|out)\)", spec, re.I)
            group = intent.group(1).lower() if intent else "local"
            description = comment.strip() or DESCRIPTIONS.get(name)
            if not description:
                description = (
                    "Call-local native argument; units and shape follow the source contract"
                    if name.startswith("v_") or name == "value"
                    else "Call-local bridge value [-]"
                )
            formatted = spec.rstrip() + " :: " + value
            formatted = formatted.replace("real(core_rknd)", "real( kind = core_rknd )")
            if len(formatted) + len(description) + 3 <= 110:
                formatted += " ! " + description
            else:
                formatted = (
                    "\n".join(
                        textwrap.wrap(
                            description,
                            width=104,
                            initial_indent="    ! ",
                            subsequent_indent="    ! ",
                        )
                    )
                    + "\n"
                    + formatted
                )
            declarations.append((group, name, formatted))
    header = lines[:start]
    if not any(line.strip() == "implicit none" for line in header):
        header += ["", "    implicit none"]
    result = []
    for line in header:
        if line.lstrip().startswith("use ") and result and result[-1]:
            result.append("")
        result.append(line)
    result.append("")
    for group, label in [
        ("in", "Input Variables"),
        ("inout", "Input/Output Variables"),
        ("out", "Output Variables"),
        ("local", "Local Variables"),
    ]:
        formals = re.search(r"subroutine\s+\w+\((.*?)\)", lines[0], re.I)
        order = {name.strip(): i for i, name in enumerate(formals.group(1).split(","))}
        members = [
            line
            for _, name, line in sorted(
                (item for item in declarations if item[0] == group),
                key=lambda item: order.get(item[1], len(order)),
            )
        ]
        if members:
            result += [
                f"    !--------------------------- {label} ---------------------------",
                *members,
                "",
            ]
    result += namelists
    if namelists:
        result.append("")
    result += [
        "    !--------------------------- Begin Code ---------------------------",
        "",
    ]
    result += spaced_body(body)
    return result + ["", lines[-1], ""]


def format_generated_fortran(text: str) -> str:
    """Format complete generated subroutines without changing conversion bodies."""
    # Preserve entries with embedded USE/namelist continuations while detecting
    # the simple declaration region emitted by the generators.
    lines = text.splitlines()
    joined = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if re.match(r"\s*(?:use\b|namelist\s*/)", line, re.I):
            while code_part(lines[i]).rstrip().endswith("&"):
                i += 1
                line += "\n" + lines[i]
        joined.append(line)
        i += 1
    result, current = [], []
    for line in joined:
        if re.match(r"\s*subroutine\b", line, re.I):
            current = [line]
        elif current:
            current.append(line)
            if re.match(r"\s*end subroutine\b", line, re.I):
                result.extend(format_routine(current))
                current = []
        else:
            if line.lstrip().startswith("use ") and result and result[-1]:
                result.append("")
            result.append(line)
    if current:
        raise ValueError("unfinished generated Fortran subroutine")
    return re.sub(r"\n{3,}", "\n\n", "\n".join(result)).rstrip() + "\n"
