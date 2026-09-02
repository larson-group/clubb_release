"""Round-trip validation for the kept Fortran-mirror routines pack_pdf_params_api / unpack_pdf_params_api.

These mirror `pdf_parameter_module.F90:{pack,unpack}_pdf_params_api` (single-column pack/unpack of the 47-slot
`pdf_parameter` ⇄ a flat (nz,47) array). They are retained for the name-mirror but otherwise unused in the
JAX (which carries pdf_params as a namedtuple), so without this they had no coverage. A pack/unpack pair must be
mutual inverses: `pack(unpack(arr)) == arr` for an arbitrary (nz,47) array. (Iter 384 — added after iter 383
removed the dead bare pack_pdf_params/unpack_pdf_params + pack_implicit_coefs_* that had no Fortran equivalent.)

Pure-Python (numpy + the derived-type module); never SKIPs.
"""
import numpy as np

import os
import sys
_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

from clubb_jax.src.CLUBB_core.pdf_params import (
    pack_pdf_params_api, unpack_pdf_params_api, init_pdf_params, pdf_parameter, implicit_coefs_terms)


def test_pack_unpack_api_roundtrip():
    nz = 12
    template = init_pdf_params(nz, 1)              # a valid 1-column pdf_parameter to fill
    # distinct value in every (level, slot) so a mis-indexed pack/unpack would not round-trip
    arr = np.arange(nz * 47, dtype=np.float64).reshape(nz, 47)
    params = unpack_pdf_params_api(arr, nz, template)        # array -> pdf_parameter
    arr2 = np.asarray(pack_pdf_params_api(params, nz))       # pdf_parameter -> array
    # pack returns (nz, 47) for the single column
    assert arr2.shape == arr.shape, f"pack shape {arr2.shape} != {arr.shape}"
    max_abs = float(np.max(np.abs(arr2 - arr)))
    assert max_abs == 0.0, f"pack(unpack(arr)) != arr — max abs diff {max_abs} (mis-indexed slot?)"
    print(f"  pack/unpack_pdf_params_api round-trip exact over {nz}×47 slots (max diff 0.0)  PASS")


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
        print(f"  {type_name} field mirror vs Fortran: SKIP (source/type-block absent)")
        return
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
    (pdf_parameter_module.F90) — the whole PDF closure reads/writes these by name + pack/unpack indexes them by position,
    so a missing/extra/renamed/reordered field would mis-populate the closure. Source-grounded. (iter 475)"""
    _check_type_mirror("pdf_parameter", pdf_parameter, min_fields=40)


def test_implicit_coefs_terms_fields_mirror_fortran_type():
    """The JAX `implicit_coefs_terms` NamedTuple mirrors the Fortran `implicit_coefs_terms` type (pdf_parameter_module.F90)
    — the implicit wp2/wpxp/xp2 closure coefs+terms the advance routines consume; same field-set/order requirement.
    (iter 476)"""
    _check_type_mirror("implicit_coefs_terms", implicit_coefs_terms, min_fields=20)


if __name__ == "__main__":
    print("test_pdf_params_pack_roundtrip:")
    test_pack_unpack_api_roundtrip()
    test_pdf_parameter_fields_mirror_fortran_type()
    test_implicit_coefs_terms_fields_mirror_fortran_type()
    print("Done.")
