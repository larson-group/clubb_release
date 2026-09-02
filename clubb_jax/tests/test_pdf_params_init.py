#!/usr/bin/env python3
"""test_pdf_params_init.py — validate the PDF-container alloc+zero routines (pdf_parameter_module.F90).

The closure carries its state in the `pdf_parameter` (49 fields) and `implicit_coefs_terms` (30 fields) NamedTuples.
test_pdf_params_pack_roundtrip pins their FIELD LAYOUT vs the Fortran `type` definitions; this pins the alloc+zero
INIT logic of the Fortran-mirrored helpers `init_pdf_params` / `zero_pdf_params_api` /
`init_pdf_implicit_coefs_terms_api` / `zero_pdf_implicit_coefs_terms_api` — every array field zeroed at shape
(ngrdcol, nz), with the passive-scalar fields None when sclr_dim==0 and (ngrdcol, nz, sclr_dim) when > 0 (mirroring
the Fortran `if (sclr_dim>0)` alloc). Pure-Python, never SKIPs. (iter 522)
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np

from clubb_jax.src.CLUBB_core.pdf_params import (
    init_pdf_params, zero_pdf_params_api,
    init_pdf_implicit_coefs_terms_api, zero_pdf_implicit_coefs_terms_api,
)

_NG, _NZ = 2, 11
_DIMS = ("ngrdcol", "nz", "sclr_dim")  # integer metadata fields, not arrays


def _check_all_zero_2d(obj, expect_shape, allow_none_for=()):
    """Every array field is an all-zero (ngrdcol, nz) array; named fields may be None (sclr_dim==0 path)."""
    for f in obj._fields:
        if f in _DIMS:
            continue
        v = getattr(obj, f)
        if v is None:
            assert f in allow_none_for, f"field {f} unexpectedly None"
            continue
        a = np.asarray(v)
        assert a.shape == expect_shape, f"field {f} shape {a.shape} != {expect_shape}"
        assert float(np.max(np.abs(a))) == 0.0, f"field {f} not zero-initialised"


def test_pdf_parameter_init():
    p = init_pdf_params(_NZ, _NG)
    assert p.ngrdcol == _NG and p.nz == _NZ
    _check_all_zero_2d(p, (_NG, _NZ))
    # zero_pdf_params_api re-derives an all-zero container from the dims
    z = zero_pdf_params_api(p)
    assert z.ngrdcol == _NG and z.nz == _NZ
    _check_all_zero_2d(z, (_NG, _NZ))
    print(f"  init_pdf_params / zero_pdf_params_api: all {len(p._fields) - 2} fields zero @ (ngrdcol, nz)  PASS")


def test_implicit_coefs_terms_init():
    # sclr_dim == 0 → the 8 passive-scalar fields are None, the rest zero @ (ngrdcol, nz)
    sclr_fields = ("coef_wp2sclrp_implicit", "term_wp2sclrp_explicit", "coef_wpsclrp2_implicit",
                   "term_wpsclrp2_explicit", "coef_wprtpsclrp_implicit", "term_wprtpsclrp_explicit",
                   "coef_wpthlpsclrp_implicit", "term_wpthlpsclrp_explicit")
    c0 = init_pdf_implicit_coefs_terms_api(_NZ, _NG, 0)
    assert c0.sclr_dim == 0
    _check_all_zero_2d(c0, (_NG, _NZ), allow_none_for=sclr_fields)
    for f in sclr_fields:
        assert getattr(c0, f) is None, f"{f} must be None when sclr_dim==0"

    # sclr_dim == 2 → the 8 passive-scalar fields are zero @ (ngrdcol, nz, 2)
    c2 = init_pdf_implicit_coefs_terms_api(_NZ, _NG, 2)
    assert c2.sclr_dim == 2
    for f in sclr_fields:
        a = np.asarray(getattr(c2, f))
        assert a.shape == (_NG, _NZ, 2) and float(np.max(np.abs(a))) == 0.0, f"{f} shape/zero (sclr_dim=2)"

    # zero_pdf_implicit_coefs_terms_api re-inits from the existing dims (idempotent on shape)
    z2 = zero_pdf_implicit_coefs_terms_api(c2)
    assert z2.sclr_dim == 2 and np.asarray(z2.coef_wp4_implicit).shape == (_NG, _NZ)
    print("  init/zero_pdf_implicit_coefs_terms_api: zeroed @ (ngrdcol,nz); sclr fields None@0 / zero@>0  PASS")


def main():
    print("test_pdf_params_init:")
    test_pdf_parameter_init()
    test_implicit_coefs_terms_init()
    print("All pdf_params init checks PASSED")


if __name__ == "__main__":
    main()
