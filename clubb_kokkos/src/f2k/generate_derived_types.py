#!/usr/bin/env python3
"""Generate temporary C++ mirrors and fieldwise Fortran converters.

Fortran declarations are the inventory. Generation fails on a declaration it
does not understand so newly added components cannot silently disappear.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from fortran_format import format_generated_fortran

ROOT = Path(__file__).resolve().parents[3]
KOKKOS = ROOT / "clubb_kokkos" / "src"
SOURCES = {
    "grid": "CLUBB_core/grid_class",
    "sclr_idx_type": "CLUBB_core/array_index",
    "clubb_config_flags_type": "CLUBB_core/model_flags",
    "nu_vertical_res_dep": "CLUBB_core/parameters_tunable",
    "pdf_parameter": "CLUBB_core/pdf_parameter_module",
    "implicit_coefs_terms": "CLUBB_core/pdf_parameter_module",
    "err_info_type": "CLUBB_core/err_info_type_module",
    "hm_metadata_type": "CLUBB_core/corr_varnce_module",
    "silhs_config_flags_type": "SILHS/parameters_silhs",
    "hydromet_pdf_parameter": "CLUBB_core/hydromet_pdf_parameter_module",
    "precipitation_fractions": "CLUBB_core/hydromet_pdf_parameter_module",
}
EXPECTED_COUNTS = [16, 6, 66, 7, 49, 30, 7, 26, 13, 18, 5]

# C++ storage slots become Native ordinals only at this temporary bridge.
ZERO_BASED_SCALAR_INDICES = {
    "grid": frozenset({"k_lb_zm", "k_lb_zt", "k_ub_zm", "k_ub_zt"}),
    "hm_metadata_type": frozenset({
        "iiNg",
        "iiNi",
        "iiNr",
        "iiNs",
        "iiPDF_Ncn",
        "iiPDF_Ng",
        "iiPDF_Ni",
        "iiPDF_Nr",
        "iiPDF_Ns",
        "iiPDF_chi",
        "iiPDF_eta",
        "iiPDF_rg",
        "iiPDF_ri",
        "iiPDF_rr",
        "iiPDF_rs",
        "iiPDF_w",
        "iirg",
        "iiri",
        "iirr",
        "iirs",
    }),
    "sclr_idx_type": frozenset({
        "iiedsclr_CO2",
        "iiedsclr_rt",
        "iiedsclr_thl",
        "iisclr_CO2",
        "iisclr_rt",
        "iisclr_thl",
    }),
}


def write_generated(path: Path, text: str) -> None:
    """Preserve modification times when regeneration would not change a file."""
    if path.suffix == ".F90":
        text = format_generated_fortran(text)
    if not path.exists() or path.read_text() != text:
        path.write_text(text)


@dataclass
class Field:
    name: str
    kind: str
    rank: int
    allocatable: bool
    dims: tuple[str, ...]
    width: int
    default: str | None
    private: bool


def split_top(value: str) -> list[str]:
    parts, depth, start = [], 0, 0
    for pos, char in enumerate(value):
        depth += (char == "(") - (char == ")")
        if char == "," and depth == 0:
            parts.append(value[start:pos].strip())
            start = pos + 1
    parts.append(value[start:].strip())
    return parts


def parse_type(name: str, source: str) -> list[Field]:
    text = (ROOT / "src" / f"{source}.F90").read_text()
    match = re.search(
        rf"(?im)^\s*type\s+{name}\s*$(.*?)^\s*end\s*type", text, re.S | re.M | re.I
    )
    if match is None:
        raise ValueError(f"missing Fortran type {name}")
    statements, pending = [], ""
    for raw in match.group(1).splitlines():
        code = raw.split("!", 1)[0].strip()
        if not code:
            continue
        continued = code.endswith("&")
        pending += " " + code.strip("& ")
        if not continued:
            statements.append(pending.strip())
            pending = ""
    if pending:
        raise ValueError(f"unfinished declaration in {name}: {pending}")
    fields = []
    for stmt in statements:
        if "::" not in stmt:
            raise ValueError(f"unsupported declaration in {name}: {stmt}")
        decl, vars_text = stmt.split("::", 1)
        attrs = split_top(decl)
        type_spec = attrs[0].lower().replace(" ", "")
        if type_spec == "integer":
            kind, width = "int", 0
        elif type_spec == "logical":
            kind, width = "logical", 0
        elif type_spec.startswith("real(") and "core_rknd" in type_spec:
            kind, width = "real", 0
        elif type_spec.startswith("character(len="):
            kind, width = "char", int(re.search(r"len=(\d+)", type_spec).group(1))
        else:
            raise ValueError(f"unsupported component type in {name}: {decl}")
        dimension = next(
            (
                a
                for a in attrs[1:]
                if a.lower().replace(" ", "").startswith("dimension(")
            ),
            None,
        )
        dims = (
            tuple(split_top(dimension[dimension.index("(") + 1 : -1]))
            if dimension
            else ()
        )
        is_alloc = any(a.lower() == "allocatable" for a in attrs[1:])
        is_private = any(a.lower() == "private" for a in attrs[1:])
        for var in split_top(vars_text):
            m = re.fullmatch(r"([A-Za-z_]\w*)(?:\s*=\s*(.*))?", var)
            if not m:
                raise ValueError(f"unsupported component in {name}: {var}")
            fields.append(
                Field(m[1], kind, len(dims), is_alloc, dims, width, m[2], is_private)
            )
    return fields


def cpp_type(field: Field) -> str:
    base = {"int": "int", "logical": "std::uint8_t", "real": "real_t", "char": "char"}[
        field.kind
    ]
    if field.rank:
        rank = field.rank + (field.kind == "char")
        return f"Kokkos::View<{base}{'*' * rank}, Kokkos::LayoutLeft>"
    if field.kind == "char":
        return "std::string"
    return base


def cpp_default(field: Field, type_name: str) -> str:
    if field.rank:
        return ""
    if type_name in {"hm_metadata_type", "sclr_idx_type"} and field.name in ZERO_BASED_SCALAR_INDICES[type_name]:
        # C++ previously zero-initialized undeclared Native integer defaults.
        # Preserve that value and declared negative sentinels bijectively.
        native_value = 0 if field.default is None else int(field.default)
        return f" = {native_value - 1}"
    if field.default is None:
        return "{}"
    value = field.default.lower().replace("_core_rknd", "")
    if value in (".true.", ".false."):
        value = "1" if value == ".true." else "0"
    return f" = {value}"


def hydromet_dimension() -> str:
    """Read the fixed hydrometeor extent used by mirrors and collection storage."""
    source = ROOT / "src/CLUBB_core/hydromet_pdf_parameter_module.F90"
    match = re.search(r"max_hydromet_dim\s*=\s*(\d+)", source.read_text(), re.I)
    if match is None:
        raise ValueError("missing max_hydromet_dim")
    return match.group(1)


def hydromet_collection(fields: list[Field]) -> list[str]:
    """Emit fieldwise device storage for the driver's column-by-level values."""
    lines = [
        "// Device storage for Fortran hydromet_pdf_params(column, level).",
        "// Each component is one View so adjacent columns stay contiguous.",
        "struct hydromet_pdf_parameter_collection {",
        f"  static constexpr int max_hydromet_dim = {hydromet_dimension()};",
        "",
        "  int ngrdcol{};",
        "  int nzt{};",
    ]
    for field in fields:
        if field.kind != "real" or field.allocatable:
            raise ValueError(f"unsupported hydromet collection field: {field.name}")
        rank = field.rank + 2
        lines.append(
            f"  Kokkos::View<real_t{'*' * rank}, Kokkos::LayoutLeft> {field.name};"
        )
    lines += [
        "",
        "  void allocate(int columns, int levels) {",
        "    // Allocate one View per Fortran component with columns contiguous.",
        "    ngrdcol = columns;",
        "    nzt = levels;",
    ]
    for field in fields:
        dims = ", ".join(["columns", "levels", *field.dims])
        lines += [
            f"    {field.name} = decltype({field.name})(",
            f'        "hydromet_{field.name}", {dims});',
        ]
    lines += [
        "  }",
        "",
        "  void zero() {",
        "    // Reset every component before the next logical runtime batch.",
    ]
    for field in fields:
        lines.append(f"    Kokkos::deep_copy({field.name}, real_t{{}});")
    lines += ["  }", "};", ""]
    return lines


def emit_headers(types: dict[str, list[Field]]) -> None:
    by_source: dict[str, list[str]] = {}
    for name, source in SOURCES.items():
        by_source.setdefault(source, []).append(name)
    for source, names in by_source.items():
        target = KOKKOS / f"{source}.hpp"
        lines = [
            "#pragma once",
            "",
            "#include <Kokkos_Core.hpp>",
            "#include <cstdint>",
            "#include <string>",
            '#include "clubb_real.hpp"',
            "",
            "namespace clubb_kokkos {",
            "",
        ]
        generated_marker = {
            "CLUBB_core/grid_class": "GRID TYPE",
            "CLUBB_core/model_flags": "MODEL FLAGS TYPE",
            "CLUBB_core/parameters_tunable": "TUNABLE PARAMETERS TYPE",
            "CLUBB_core/corr_varnce_module": "HM METADATA TYPE",
            "CLUBB_core/err_info_type_module": "ERR INFO TYPE",
            "CLUBB_core/hydromet_pdf_parameter_module": "HYDROMET PDF TYPES",
            "CLUBB_core/pdf_parameter_module": "PDF PARAMETER TYPES",
            "SILHS/parameters_silhs": "SILHS CONFIG FLAGS TYPE",
        }.get(source)
        if generated_marker:
            lines.append(f"// BEGIN GENERATED {generated_marker}")
        if source == "CLUBB_core/model_flags":
            native = (ROOT / "src" / f"{source}.F90").read_text()
            lines += ["// Native compile-time options used by the ported core."]
            for name in (
                "order_xm_wpxp",
                "order_xp2_xpyp",
                "order_wp2_wp3",
                "order_windm",
                "iiPDF_ADG1",
                "iiPDF_ADG2",
                "iiPDF_3D_Luhar",
                "iiPDF_new",
                "iiPDF_TSDADG",
                "iiPDF_LY93",
                "iiPDF_new_hybrid",
                "saturation_bolton",
                "saturation_gfdl",
                "saturation_flatau",
                "saturation_lookup",
                "I_sat_sphum",
                "ipdf_pre_advance_fields",
                "ipdf_post_advance_fields",
                "ipdf_pre_post_advance_fields",
                "l_advance_xp3",
                "l_morr_xp2_mc",
                "l_evaporate_cold_rcm",
                "l_silhs_rad",
                "l_use_boussinesq",
                "l_test_grid_generalization",
                "l_explicit_turbulent_adv_wp3",
                "l_upwind_Kh_dp_term",
                "global_fill",
                "sliding_window",
                "l_pos_def",
                "l_clip_turb_adv",
                "l_explicit_turbulent_adv_wpxp",
                "l_force_descending_solves",
                "l_explicit_turbulent_adv_xpyp",
                "l_gamma_Skw",
                "l_cubic_interp",
                "l_quintic_poly_interp",
                "Lscale_direct",
                "Lscale_tau",
                "Lscale_coherence",
                "cons_ullrich_remap",
                "ppm_remap",
                "no_grid_adaptation",
                "lapack",
                "penta_lu",
                "tridiag_lu",
                "penta_bicgstab",
            ):
                match = re.search(
                    rf"(?im)^\s*{name}\s*=\s*(\d+|\.false\.|\.true\.)(?=\s*(?:,|!|$))",
                    native,
                )
                if match is None:
                    raise ValueError(f"unsupported native constant: {name}")
                value = match.group(1).lower()
                logical = value.startswith(".")
                if logical:
                    value = "true" if value == ".true." else "false"
                lines.append(
                    f"inline constexpr {'bool' if logical else 'int'} {name} = {value};"
                )
            lines.append("")
        for name in names:
            lines += [
                f"// Mirrors src/{source}.F90; arrays use Fortran-order LayoutLeft.",
                f"struct {name} {{",
            ]
            for field in types[name]:
                lines.append(f"  {cpp_type(field)} {field.name}{cpp_default(field, name)};")
            lines += ["};", ""]
            if name == "hydromet_pdf_parameter":
                lines += hydromet_collection(types[name])
        if generated_marker:
            begin_marker = f"// BEGIN GENERATED {generated_marker}"
            end_marker = f"// END GENERATED {generated_marker}"
            lines.append(end_marker)
            if target.exists():
                # Permanent APIs surround this generated type declaration.
                # Regeneration owns only the field inventory between the markers.
                existing = target.read_text()
                if existing.count(begin_marker) != 1 or existing.count(end_marker) != 1:
                    raise ValueError(f"{target.name} needs one generated-type marker pair")
                begin = existing.index(begin_marker)
                end = existing.index(end_marker, begin) + len(end_marker)
                generated = "\n".join(lines)
                generated = generated[generated.index(begin_marker) :]
                write_generated(target, existing[:begin] + generated + existing[end:])
                continue
        lines += ["}  // namespace clubb_kokkos", ""]
        target.parent.mkdir(parents=True, exist_ok=True)
        write_generated(target, "\n".join(lines))


def emit_real_header() -> None:
    write_generated(
        KOKKOS / "clubb_real.hpp",
        """#pragma once

namespace clubb_kokkos {
#if CLUBB_REAL_TYPE == 4
using real_t = float;
#elif CLUBB_REAL_TYPE == 8
using real_t = double;
#else
#error Unsupported CLUBB_REAL_TYPE
#endif
}  // namespace clubb_kokkos
""",
    )


def native(field: Field) -> str:
    return f"private_{field.name}" if field.private else f"value%{field.name}"


def extent_expr(field: Field, target: str) -> str:
    dims = [f"size({target}, {i}, kind=c_int64_t)" for i in range(1, field.rank + 1)]
    return "[" + ", ".join((dims + ["1_c_int64_t"] * (4 - len(dims)))) + "]"


def field_shape_assignments(target: str, rank: int) -> list[str]:
    # NVFORTRAN 26.9 mispacks SHAPE(..., kind=c_int64_t) for derived-type components.
    return [
        f"      field_shape({i}) = size({target}, {i}, kind=c_int64_t)"
        for i in range(1, rank + 1)
    ]


def fortran_input(field: Field, index: int) -> list[str]:
    target = native(field)
    kind_ptr = {"int": "ip", "logical": "bp", "real": "rp", "char": "cp"}[field.kind]
    lines = [
        f"    ! {field.name}",
        f"    if (input({index})%rank /= {field.rank}) error stop 'f2k rank: {field.name}'",
    ]
    if field.rank:
        lines += [f"    if (input({index})%allocated /= 0) then"]
        dims = [f"int(input({index})%extent({i}))" for i in range(1, field.rank + 1)]
        if field.allocatable or field.private:
            lines += [f"      allocate({target}({', '.join(dims)}))"]
        else:
            lines += [
                "      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]"
            ]
            lines += field_shape_assignments(target, field.rank)
            lines += [
                f"      if (any(input({index})%extent /= field_shape)) &",
                f"        error stop 'f2k shape: {field.name}'",
            ]
        count = f"int(product(input({index})%extent(1:{field.rank})))"
        if field.kind == "char":
            lines += [
                f"      if (input({index})%char_len /= {field.width}) error stop 'f2k width: {field.name}'",
                f"      n = {count}",
                f"      if (n > 0) then",
                f"        if (.not. c_associated(input({index})%data)) error stop 'f2k data: {field.name}'",
                f"        call c_f_pointer(input({index})%data, cp, [n * {field.width}])",
                f"        do j = 1, n",
                f"          do k = 1, {field.width}",
                f"            {target}(j)(k:k) = cp((j-1)*{field.width}+k)",
                "          end do",
                "        end do",
                "      end if",
            ]
        else:
            lines += [
                f"      n = {count}",
                "      if (n > 0) then",
                f"        if (.not. c_associated(input({index})%data)) error stop 'f2k data: {field.name}'",
                f"        call c_f_pointer(input({index})%data, {kind_ptr}, [n])",
            ]
            expr = f"{kind_ptr}"
            if field.kind == "logical":
                expr = "bp /= 0_c_int8_t"
            lines += [
                f"        {target} = reshape({expr}, shape({target}))",
                "      end if",
            ]
        if not field.allocatable and not field.private:
            lines += [
                "    else",
                f"      error stop 'f2k missing fixed array: {field.name}'",
            ]
        lines += ["    end if"]
    else:
        lines += [
            f"    if (.not. c_associated(input({index})%data)) error stop 'f2k scalar: {field.name}'"
        ]
        if field.kind == "char":
            lines += [
                f"    if (input({index})%char_len /= {field.width}) error stop 'f2k width: {field.name}'",
                f"    call c_f_pointer(input({index})%data, cp, [{field.width}])",
                f"    do k = 1, {field.width}",
                f"      {target}(k:k) = cp(k)",
                "    end do",
            ]
        else:
            lines += [f"    call c_f_pointer(input({index})%data, {kind_ptr}, [1])"]
            expr = f"{kind_ptr}(1)"
            if field.kind == "logical":
                expr = "bp(1) /= 0_c_int8_t"
            lines += [f"    {target} = {expr}"]
    return lines


def fortran_output(field: Field, index: int) -> list[str]:
    target = native(field)
    kind_ptr = {"int": "ip", "logical": "bp", "real": "rp", "char": "cp"}[field.kind]
    bytes_expr = {
        "int": "c_sizeof(0_c_int)",
        "logical": "c_sizeof(0_c_int8_t)",
        "real": "c_sizeof(0.0_core_rknd)",
        "char": "1_c_size_t",
    }[field.kind]
    lines = [f"    ! {field.name}"]
    if field.rank:
        if field.allocatable or field.private:
            lines += [f"    if (allocated({target})) then"]
        else:
            lines += ["    block"]
        lines += [
            f"      n = size({target})",
            "      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]",
        ]
        lines += field_shape_assignments(target, field.rank)
        lines += [
            f"      call allocate_field({field.rank}, field_shape, & ! In",
            f"                          int(n, c_size_t) * {field.width if field.kind == 'char' else 1} * {bytes_expr}, {field.width}, & ! In",
            f"                          output({index})) ! Out",
        ]
        if field.kind == "char":
            lines += [
                f"      call c_f_pointer(output({index})%data, cp, [n * {field.width}])",
                "      do j = 1, n",
                f"        do k = 1, {field.width}",
                f"          cp((j-1)*{field.width}+k) = {target}(j)(k:k)",
                "        end do",
                "      end do",
            ]
        else:
            lines += [f"      call c_f_pointer(output({index})%data, {kind_ptr}, [n])"]
            expr = f"reshape({target}, [n])"
            if field.kind == "logical":
                expr = f"merge(1_c_int8_t, 0_c_int8_t, reshape({target}, [n]))"
            lines += [f"      {kind_ptr} = {expr}"]
        if field.allocatable or field.private:
            lines += [
                "    else",
                f"      call empty_field({field.rank}, {field.width}, & ! In",
                f"                       output({index})) ! Out",
                "    end if",
            ]
        else:
            lines += ["    end block"]
    else:
        size_expr = f"{field.width}" if field.kind == "char" else "1"
        lines += [
            "    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In",
            f"                        int({size_expr}, c_size_t) * {bytes_expr}, {field.width}, & ! In",
            f"                        output({index})) ! Out",
            f"    call c_f_pointer(output({index})%data, {kind_ptr}, [{size_expr}])",
        ]
        if field.kind == "char":
            lines += [
                f"    do k = 1, {field.width}",
                f"      cp(k) = {target}(k:k)",
                "    end do",
            ]
        else:
            expr = (
                target
                if field.kind != "logical"
                else f"merge(1_c_int8_t, 0_c_int8_t, {target})"
            )
            lines += [f"    {kind_ptr}(1) = {expr}"]
    return lines


def collection_fortran(fields: list[Field]) -> list[str]:
    """Map one native array of fixed-size UDTs to column-first field Views."""
    name = "hydromet_pdf_parameter_collection"
    width = hydromet_dimension()
    lines = []
    for direction in ("to", "from"):
        if direction == "to":
            signature = f"  subroutine to_native_{name}(input, value)"
            args = [
                "    type(wire_field), intent(in) :: input(*)",
                "    type(hydromet_pdf_parameter), allocatable, intent(out) :: value(:,:)",
            ]
        else:
            signature = f"  subroutine from_native_{name}(value, output)"
            args = [
                "    type(hydromet_pdf_parameter), intent(in) :: value(:,:)",
                "    type(wire_field), intent(out) :: output(*)",
            ]
        lines += [
            signature,
            "    ! Description: Transfer the driver collection by component without saved state.",
            *args,
            "    real(core_rknd), pointer :: rp(:)",
            "    integer :: ncol, nzt, cells, n, h, j, offset",
            "    integer(c_int64_t) :: field_shape(4)",
        ]
        if direction == "to":
            lines += [
                "    ncol = int(input(1)%extent(1))",
                "    nzt = int(input(1)%extent(2))",
                "    allocate(value(ncol,nzt))",
            ]
        else:
            lines += ["    ncol = size(value,1)", "    nzt = size(value,2)"]
        lines += ["    cells = ncol*nzt"]
        for index, field in enumerate(fields, 1):
            dims = ["int(ncol,c_int64_t)", "int(nzt,c_int64_t)"] + [
                f"{width}_c_int64_t"
            ] * field.rank
            dims += ["1_c_int64_t"] * (4 - len(dims))
            rank = field.rank + 2
            lines += [
                f"    ! {field.name}",
                f"    field_shape = [{', '.join(dims)}]",
                f"    n = cells*{int(width) ** field.rank}",
            ]
            if direction == "to":
                lines += [
                    f"    if (input({index})%rank /= {rank} .or. input({index})%allocated == 0) &",
                    f"      error stop 'f2k collection allocation: {field.name}'",
                    f"    if (any(input({index})%extent /= field_shape)) error stop 'f2k collection shape'",
                    "    if (n > 0) then",
                    f"      call c_f_pointer(input({index})%data,rp,[n])",
                ]
            else:
                lines += [
                    f"    call allocate_field({rank}, field_shape, & ! In",
                    "                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In",
                    f"                        output({index})) ! Out",
                    "    if (n > 0) then",
                    f"      call c_f_pointer(output({index})%data,rp,[n])",
                ]
            if field.rank == 2:
                lines += [
                    f"      do j = 1, {width}",
                    f"        do h = 1, {width}",
                    f"          offset = ((j-1)*{width}+h-1)*cells",
                ]
                member = f"value(:,:)%{field.name}(h,j)"
                indent = "          "
            elif field.rank == 1:
                lines += [f"      do h = 1, {width}", "        offset = (h-1)*cells"]
                member = f"value(:,:)%{field.name}(h)"
                indent = "        "
            else:
                lines += ["      offset = 0"]
                member = f"value(:,:)%{field.name}"
                indent = "      "
            if direction == "to":
                lines.append(
                    f"{indent}{member} = reshape(rp(offset+1:offset+cells),[ncol,nzt])"
                )
            else:
                lines.append(
                    f"{indent}rp(offset+1:offset+cells) = reshape({member},[cells])"
                )
            if field.rank == 2:
                lines += ["        end do", "      end do"]
            elif field.rank == 1:
                lines += ["      end do"]
            lines += ["    end if"]
        lines += [f"  end subroutine {direction}_native_{name}", ""]
    lines += [
        f"  subroutine f2k_roundtrip_{name}_c(input,output) bind(C)",
        "    ! Description: Check both collection converters with a local native array.",
        "    type(wire_field), intent(in) :: input(*)",
        "    type(wire_field), intent(out) :: output(*)",
        "    type(hydromet_pdf_parameter), allocatable :: value(:,:)",
        f"    call to_native_{name}(input,value)",
        f"    call from_native_{name}(value,output)",
        f"  end subroutine f2k_roundtrip_{name}_c",
        "",
    ]
    return lines


def collection_cpp(fields: list[Field]) -> list[str]:
    """Emit native collection packing using the same array transfer implementation."""
    name = "hydromet_pdf_parameter_collection"
    lines = [
        f'extern "C" void f2k_roundtrip_{name}_c(const WireField*, WireField*);',
        "",
        f"Transfer to_fortran(const clubb_kokkos::{name}& value) {{",
        f"  Transfer transfer({len(fields)});",
    ]
    for index, field in enumerate(fields):
        shape = (
            [-1, -1] + [int(hydromet_dimension())] * field.rank + [1] * (2 - field.rank)
        )
        lines.append(
            f"  transfer.fields[{index}] = pack_view(value.{field.name}, transfer, "
            f'{{{", ".join(map(str,shape))}}}, "{field.name}");'
        )
    lines += [
        "  return transfer;",
        "}",
        "",
        f"Transfer make_output(const clubb_kokkos::{name}&) {{",
        f"  return Transfer({len(fields)}, true);",
        "}",
        "",
        f"void to_kokkos(clubb_kokkos::{name}& value, const Transfer& output) {{",
        f'  if (output.fields.size() != {len(fields)}) throw std::runtime_error("f2k collection field count");',
    ]
    for index, field in enumerate(fields):
        shape = (
            [-1, -1] + [int(hydromet_dimension())] * field.rank + [1] * (2 - field.rank)
        )
        lines.append(
            f"  unpack_view(value.{field.name}, output.fields[{index}], "
            f'{{{", ".join(map(str,shape))}}}, "{field.name}");'
        )
    lines += [
        "  value.ngrdcol = static_cast<int>(value.hm_1.extent(0));",
        "  value.nzt = static_cast<int>(value.hm_1.extent(1));",
        "}",
        "",
        f"void roundtrip(clubb_kokkos::{name}& value) {{",
        "  auto input = to_fortran(value);",
        "  auto output = make_output(value);",
        f"  f2k_roundtrip_{name}_c(input.data(), output.data());",
        "  to_kokkos(value, output);",
        "}",
        "",
    ]
    return lines


def emit_fortran(types: dict[str, list[Field]]) -> None:
    lines = [
        "! Temporary fieldwise conversions for the C++/Fortran call boundary.",
        "! to_native_* makes a local Fortran value from C field descriptors;",
        "! from_native_* returns every field in C-owned buffers. The BIND(C)",
        "! roundtrip procedures exercise both directions without saved state.",
        "module f2k_derived_type_bridge",
        "  use, intrinsic :: iso_c_binding, only: &",
        "    c_ptr, c_null_ptr, c_int, c_int64_t, c_int8_t, c_size_t, c_char, &",
        "    c_f_pointer, c_associated, c_sizeof",
        "  use clubb_precision, only: core_rknd",
    ]
    for source in dict.fromkeys(SOURCES.values()):
        module = Path(source).name
        names = [name for name, src in SOURCES.items() if src == source]
        if module == "err_info_type_module":
            names += ["f2k_get_err_info_private_api", "f2k_set_err_info_private_api"]
        lines += [f"  use {module}, only: &", "    " + ", &\n    ".join(names)]
    lines += [
        "  implicit none",
        "  private",
        "  public :: wire_field, allocate_field, empty_field",
        "  public :: "
        + ", &\n    ".join(f"to_native_{name}, from_native_{name}" for name in types),
        "  public :: to_native_hydromet_pdf_parameter_collection, from_native_hydromet_pdf_parameter_collection",
        "  type, bind(C) :: wire_field",
        "    type(c_ptr) :: data",
        "    integer(c_int64_t) :: extent(4)",
        "    integer(c_int) :: rank, allocated, char_len",
        "  end type wire_field",
        "  interface",
        "    function c_malloc(nbytes) bind(C, name='malloc') result(ptr)",
        "      import c_size_t, c_ptr",
        "      integer(c_size_t), value :: nbytes",
        "      type(c_ptr) :: ptr",
        "    end function c_malloc",
        "  end interface",
        "contains",
        "",
        "  subroutine allocate_field(rank, extent, nbytes, char_len, field)",
        "    ! Description: Allocate C-owned output bytes for one converted field.",
        "    type(wire_field), intent(out) :: field",
        "    integer, intent(in) :: rank, char_len",
        "    integer(c_int64_t), intent(in) :: extent(4)",
        "    integer(c_size_t), intent(in) :: nbytes",
        "    field%data = c_malloc(max(1_c_size_t, nbytes))",
        "    if (.not. c_associated(field%data)) error stop 'f2k malloc failed'",
        "    field%rank = rank",
        "    field%extent = extent",
        "    field%allocated = 1",
        "    field%char_len = char_len",
        "  end subroutine allocate_field",
        "",
        "  subroutine empty_field(rank, char_len, field)",
        "    ! Description: Mark an unallocated Fortran component in the C descriptor.",
        "    type(wire_field), intent(out) :: field",
        "    integer, intent(in) :: rank, char_len",
        "    field%data = c_null_ptr",
        "    field%extent = 0_c_int64_t",
        "    field%rank = rank",
        "    field%allocated = 0",
        "    field%char_len = char_len",
        "  end subroutine empty_field",
        "",
    ]
    for name, fields in types.items():
        lines += [
            f"  subroutine to_native_{name}(input, value)",
            f"    ! Description: Copy each C field into a local {name} value.",
            "    type(wire_field), intent(in) :: input(*)",
            f"    type({name}), intent(out) :: value",
            "    integer(c_int), pointer :: ip(:)",
            "    integer(c_int8_t), pointer :: bp(:)",
            "    real(core_rknd), pointer :: rp(:)",
            "    character(c_char), pointer :: cp(:)",
            "    integer :: n, j, k",
            "    integer(c_int64_t) :: field_shape(4)",
        ]
        if name == "err_info_type":
            lines += [
                "    real(core_rknd), allocatable :: private_host_lat(:), private_host_lon(:)",
                "    integer :: private_chunk_idx, private_mpi_rank",
            ]
        for index, field in enumerate(fields, 1):
            lines += fortran_input(field, index)
        if name == "err_info_type":
            lines += [
                "    call f2k_set_err_info_private_api(private_host_lat, private_host_lon, & ! In",
                "                                      private_chunk_idx, private_mpi_rank, & ! In",
                "                                      value)                                  ! InOut",
            ]
        lines += [f"  end subroutine to_native_{name}", ""]
        lines += [
            f"  subroutine from_native_{name}(value, output)",
            f"    ! Description: Copy each {name} field into C-owned output bytes.",
            f"    type({name}), intent(in) :: value",
            "    type(wire_field), intent(out) :: output(*)",
            "    integer(c_int), pointer :: ip(:)",
            "    integer(c_int8_t), pointer :: bp(:)",
            "    real(core_rknd), pointer :: rp(:)",
            "    character(c_char), pointer :: cp(:)",
            "    integer :: n, j, k",
            "    integer(c_int64_t) :: field_shape(4)",
        ]
        if name == "err_info_type":
            lines += [
                "    real(core_rknd), allocatable :: private_host_lat(:), private_host_lon(:)",
                "    integer :: private_chunk_idx, private_mpi_rank",
                "    call f2k_get_err_info_private_api(value, &                             ! In",
                "                                      private_host_lat, private_host_lon, & ! Out",
                "                                      private_chunk_idx, private_mpi_rank) ! Out",
            ]
        for index, field in enumerate(fields, 1):
            lines += fortran_output(field, index)
        lines += [f"  end subroutine from_native_{name}", ""]
        lines += [
            f"  subroutine f2k_roundtrip_{name}_c(input, output) bind(C)",
            "    ! Description: Exercise both converters with one local Fortran value.",
            "    type(wire_field), intent(in) :: input(*)",
            "    type(wire_field), intent(out) :: output(*)",
            f"    type({name}) :: value",
            f"    call to_native_{name}(input, value)",
            f"    call from_native_{name}(value, output)",
            f"  end subroutine f2k_roundtrip_{name}_c",
            "",
        ]
    lines += collection_fortran(types["hydromet_pdf_parameter"])
    lines += ["end module f2k_derived_type_bridge", ""]
    write_generated(KOKKOS / "f2k" / "derived_type_bridge.F90", "\n".join(lines))


def cpp_shape(field: Field) -> str:
    if not field.rank:
        return "{1, 1, 1, 1}"
    if field.allocatable:
        return "{" + ", ".join(["-1"] * field.rank + ["1"] * (4 - field.rank)) + "}"
    dims = [hydromet_dimension() if d == "max_hydromet_dim" else d for d in field.dims]
    return "{" + ", ".join(dims + ["1"] * (4 - field.rank)) + "}"


def emit_cpp(types: dict[str, list[Field]]) -> None:
    path = KOKKOS / "f2k" / "derived_type_bridge.cpp"
    prefix = path.read_text().split("// Generated per-type field lists follow.", 1)[0]
    lines = ["// Generated per-type field lists follow.", ""]
    for name, fields in types.items():
        lines += [
            f'extern "C" void f2k_roundtrip_{name}_c(const WireField*, WireField*);',
            "",
            f"Transfer to_fortran(const clubb_kokkos::{name}& value) {{",
            f"  Transfer transfer({len(fields)});",
        ]
        for index, field in enumerate(fields):
            member = f"value.{field.name}"
            if field.kind == "char" and field.rank:
                call = (
                    f'pack_char_view({member}, {field.width}, transfer, "{field.name}")'
                )
            elif field.kind == "char":
                call = f"pack_string({member}, {field.width}, transfer)"
            elif field.rank:
                call = (
                    f'pack_view({member}, transfer, {cpp_shape(field)}, "{field.name}")'
                )
            else:
                scalar = f"{member} + 1" if field.name in ZERO_BASED_SCALAR_INDICES.get(name, ()) else member
                call = f"pack_scalar({scalar}, transfer)"
            lines.append(f"  transfer.fields[{index}] = {call};")
        lines += [
            "  return transfer;",
            "}",
            "",
            f"Transfer make_output(const clubb_kokkos::{name}&) {{",
            f"  return Transfer({len(fields)}, true);",
            "}",
            "",
            f"void to_kokkos(clubb_kokkos::{name}& value, const Transfer& output) {{",
            f'  if (output.fields.size() != {len(fields)}) throw std::runtime_error("f2k field count");',
        ]
        for index, field in enumerate(fields):
            member = f"value.{field.name}"
            arg = f"output.fields[{index}]"
            if field.kind == "char" and field.rank:
                call = (
                    f'unpack_char_view({member}, {arg}, {field.width}, "{field.name}")'
                )
            elif field.kind == "char":
                call = f'unpack_string({member}, {arg}, {field.width}, "{field.name}")'
            elif field.rank:
                call = (
                    f'unpack_view({member}, {arg}, {cpp_shape(field)}, "{field.name}")'
                )
            else:
                call = f'unpack_scalar({member}, {arg}, "{field.name}")'
            lines.append(f"  {call};")
            if not field.rank and field.name in ZERO_BASED_SCALAR_INDICES.get(name, ()):
                lines.append(f"  {member} -= 1;")
        lines += [
            "}",
            "",
            f"void roundtrip(clubb_kokkos::{name}& value) {{",
            "  auto input = to_fortran(value);",
            "  auto output = make_output(value);",
            f"  f2k_roundtrip_{name}_c(input.data(), output.data());",
            "  to_kokkos(value, output);",
            "}",
            "",
        ]
    lines += collection_cpp(types["hydromet_pdf_parameter"])
    lines += ["}  // namespace f2k", ""]
    write_generated(path, prefix + "\n".join(lines))


def emit_test(types: dict[str, list[Field]]) -> None:
    path = (
        ROOT
        / "clubb_kokkos"
        / "auto_llm_generated_tests"
        / "derived_type_roundtrip.cpp"
    )
    prefix = path.read_text().split("// Generated fieldwise checks follow.", 1)[0]
    lines = ["// Generated fieldwise checks follow."]
    for name, fields in types.items():
        if name in {"hm_metadata_type", "sclr_idx_type"}:
            lines += ["    {", f"      clubb_kokkos::{name} defaults;", "      auto packed = f2k::to_fortran(defaults);"]
            for index, field in enumerate(fields):
                if field.name not in ZERO_BASED_SCALAR_INDICES[name]:
                    continue
                native_value = 0 if field.default is None else int(field.default)
                lines += [
                    f'      same_scalar({native_value - 1}, defaults.{field.name}, "{name}%{field.name} default slot");',
                    "      {",
                    "        int native_value;",
                    f"        std::memcpy(&native_value, packed.fields[{index}].data, sizeof(native_value));",
                    f'        same_scalar({native_value}, native_value, "{name}%{field.name} default wire ordinal");',
                    "      }",
                ]
            lines += ["    }", "    for (int native_index : {-1, 0, 1, 2, 19}) {", f"      clubb_kokkos::{name} value;"]
            for field in fields:
                if field.name in ZERO_BASED_SCALAR_INDICES[name]:
                    lines.append(f"      value.{field.name} = native_index - 1;")
            lines += ["      auto packed = f2k::to_fortran(value);"]
            for index, field in enumerate(fields):
                if field.name in ZERO_BASED_SCALAR_INDICES[name]:
                    lines += ["      {", "        int wire_index;",
                        f"        std::memcpy(&wire_index, packed.fields[{index}].data, sizeof(wire_index));",
                        f'        same_scalar(native_index, wire_index, "{name}%{field.name} wire ordinal");', "      }"]
            lines += ["      f2k::roundtrip(value);"]
            for field in fields:
                if field.name in ZERO_BASED_SCALAR_INDICES[name]:
                    lines.append(f'      same_scalar(native_index - 1, value.{field.name}, "{name}%{field.name} sentinel roundtrip");')
            lines += ["    }"]
        lines += ["    {", f"      clubb_kokkos::{name} value;"]
        for index, field in enumerate(fields, 1):
            member = f"value.{field.name}"
            if field.rank:
                if field.kind == "char":
                    dims = [str(field.width), "2"]
                elif field.allocatable:
                    dims = [str(i + 2) for i in range(field.rank)]
                else:
                    dims = [
                        hydromet_dimension() if d == "max_hydromet_dim" else d
                        for d in field.dims
                    ]
                lines += [
                    f'      {member} = decltype({member})("{field.name}", {", ".join(dims)});',
                    f"      fill({member}, {index});",
                ]
            elif field.kind == "char":
                lines += [
                    f"      {member} = std::string({field.width}, '{chr(65 + index % 26)}');"
                ]
            elif field.kind == "logical":
                lines += [f"      {member} = {index % 2};"]
            elif field.kind == "real":
                lines += [
                    f"      {member} = static_cast<clubb_kokkos::real_t>({index} + 0.125);"
                ]
            else:
                lines += [f"      {member} = {100 + index};"]
        if name == "grid":
            lines += [
                "      value.grid_dir = static_cast<clubb_kokkos::real_t>(-0.0);",
                "      auto host_zm = Kokkos::create_mirror_view(value.zm);",
                "      host_zm(0, 0) = static_cast<clubb_kokkos::real_t>(-0.0);",
                "      Kokkos::deep_copy(value.zm, host_zm);",
            ]
        lines += ["      auto before = value;", "      f2k::roundtrip(value);"]
        for field in fields:
            compare = "same_view" if field.rank else "same_scalar"
            lines += [
                f'      {compare}(before.{field.name}, value.{field.name}, "{name}%{field.name}");'
            ]
        first_alloc = next((f for f in fields if f.allocatable), None)
        if first_alloc is not None:
            lines += [
                f"      value.{first_alloc.name} = {{}};",
                "      f2k::roundtrip(value);",
                f"      if (value.{first_alloc.name}.is_allocated())",
                f'        throw std::runtime_error("{name}%{first_alloc.name} allocation");',
            ]
            if first_alloc.kind != "char":
                dims = ["0"] + [str(i + 2) for i in range(1, first_alloc.rank)]
                lines += [
                    f'      value.{first_alloc.name} = decltype(value.{first_alloc.name})("empty", {", ".join(dims)});',
                    "      f2k::roundtrip(value);",
                    f"      if (!value.{first_alloc.name}.is_allocated() || value.{first_alloc.name}.extent(0) != 0)",
                    f'        throw std::runtime_error("{name}%{first_alloc.name} empty shape");',
                ]
        if name == "err_info_type":
            lines += [
                "      value.err_header = {};",
                "      f2k::roundtrip(value);",
                "      if (value.err_header.is_allocated())",
                '        throw std::runtime_error("err_info_type%err_header allocation");',
                '      value.err_header = decltype(value.err_header)("empty_headers", 200, 0);',
                "      f2k::roundtrip(value);",
                "      if (!value.err_header.is_allocated() || value.err_header.extent(1) != 0)",
                '        throw std::runtime_error("err_info_type%err_header empty shape");',
            ]
        lines += ["    }"]
    lines += ["  }", "  Kokkos::finalize();", "  return 0;", "}", ""]
    write_generated(path, prefix + "\n".join(lines))


def main() -> None:
    types = {name: parse_type(name, source) for name, source in SOURCES.items()}
    for (name, fields), expected in zip(types.items(), EXPECTED_COUNTS, strict=True):
        if len(fields) != expected:
            raise ValueError(f"{name}: expected {expected} fields, got {len(fields)}")
    emit_real_header()
    emit_headers(types)
    emit_fortran(types)
    emit_cpp(types)
    emit_test(types)


if __name__ == "__main__":
    main()
