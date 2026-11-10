"""Check the PDF and implicit-coefficient type fields against canonical Fortran declarations."""
from utilities.output_paths import REPO_ROOT as _REPO_ROOT

import os
_ROOT = str(_REPO_ROOT)

from clubb_jax.src.CLUBB_core.pdf_parameter_module import pdf_parameter, implicit_coefs_terms


_PDF_MODULE_F90 = os.path.join(_ROOT, "src", "CLUBB_core", "pdf_parameter_module.F90")


def _fortran_type_fields(type_name):
    """Parse the field names (in order) of a Fortran derived `type <type_name> … end type` from pdf_parameter_module.F90.
    Returns None if the source is absent or the block isn't found."""
    import re
    if not os.path.exists(_PDF_MODULE_F90):
        return None
    lines = open(_PDF_MODULE_F90).read().splitlines()
    try:
        i0 = next(i for i, l in enumerate(lines) if re.match(rf"\s*type\s+{type_name}\s*$", l))
        i1 = next(i for i, l in enumerate(lines) if re.match(rf"\s*end\s+type\s+{type_name}", l))
    except StopIteration:
        return None
    fields, in_decl = [], False
    for ln in lines[i0:i1]:
        code = ln.split("!", 1)[0]
        if "::" in code:                               # declaration header — fields may trail on this line
            in_decl = True
            code = code.split("::", 1)[1]
        elif not in_decl:
            continue
        for tok in code.replace("&", "").split(","):
            t = tok.strip()
            if re.match(r"^[A-Za-z]\w*$", t):
                fields.append(t)
    return fields


def _check_type_mirror(type_name, jax_tuple, min_fields):
    """Assert the JAX NamedTuple `_fields` == the Fortran `type_name` fields, field-for-field, in order."""
    fort = _fortran_type_fields(type_name)
    if fort is None:
        raise AssertionError(f"  {type_name} field mirror vs Fortran: required input missing (source/type-block absent)")
    jax = list(jax_tuple._fields)
    assert len(fort) >= min_fields, f"parsed only {len(fort)} {type_name} fields from the F90 — extraction broke"
    assert fort == jax, (f"JAX {type_name} fields diverge from the Fortran type (set/order):\n"
                         f"  in Fortran not JAX: {sorted(set(fort) - set(jax))}\n"
                         f"  in JAX not Fortran: {sorted(set(jax) - set(fort))}\n"
                         f"  (first order diff at index "
                         f"{next((i for i,(a,b) in enumerate(zip(fort,jax)) if a!=b), 'n/a')})")
    print(f"  {type_name} fields mirror the Fortran type exactly ({len(fort)} fields, same order)  PASS")


def test_pdf_parameter_fields_mirror_fortran_type():
    """The JAX `pdf_parameter` NamedTuple has the SAME fields, in the SAME order, as the Fortran `pdf_parameter` type
    (pdf_parameter_module.F90) — the PDF closure reads and writes these fields by name. Source-grounded. (iter 475)"""
    _check_type_mirror("pdf_parameter", pdf_parameter, min_fields=40)


def test_implicit_coefs_terms_fields_mirror_fortran_type():
    """The JAX `implicit_coefs_terms` NamedTuple mirrors the Fortran `implicit_coefs_terms` type (pdf_parameter_module.F90)
    — the implicit wp2/wpxp/xp2 closure coefs+terms the advance routines consume; same field-set/order requirement.
    (iter 476)"""
    _check_type_mirror("implicit_coefs_terms", implicit_coefs_terms, min_fields=20)
