"""Validate the JAX index_mapping port against the index_mapping.F90 logic.

These are pure integer/constant index maps over hm_metadata (no FP oracle needed):
the Fortran semantics are exact equality tests on the metadata indices, so this
checks the JAX reproduces them — including the absent-species (-1) guard, which the
Fortran gets for free from positive 1-based indexing but the 0-based JAX must guard
explicitly so a -1 query never matches a -1 (absent) metadata field.
"""
from utilities.output_paths import REPO_ROOT as _REPO_ROOT
import os
_ROOT = str(_REPO_ROOT)

from clubb_jax.src.CLUBB_core.corr_varnce_module import HmMetadata
from clubb_jax.src.CLUBB_core.index_mapping import (
    pdf2hydromet_idx, hydromet2pdf_idx, rx2Nx_hm_idx, Nx2rx_hm_idx, mvr_hm_max)
from clubb_jax.src.CLUBB_core.constants_clubb import (
    mvr_rain_max, mvr_ice_max, mvr_snow_max, mvr_graupel_max)


def test_warm_rain_only_pdf():
    """KK rain-only PDF: hydromet array [rr=0, Nr=1]; PDF array chi,eta,w,Ncn,rr=4,Nr=5."""
    hm = HmMetadata(hydromet_dim=2, iirr=0, iiNr=1, iiPDF_rr=4, iiPDF_Nr=5)
    assert pdf2hydromet_idx(4, hm) == 0 and pdf2hydromet_idx(5, hm) == 1
    assert pdf2hydromet_idx(2, hm) == -1            # w is not a hydromet
    assert hydromet2pdf_idx(0, hm) == 4 and hydromet2pdf_idx(1, hm) == 5
    assert rx2Nx_hm_idx(0, hm) == 1 and Nx2rx_hm_idx(1, hm) == 0
    assert mvr_hm_max(0, hm) == mvr_rain_max and mvr_hm_max(1, hm) == mvr_rain_max
    print("  warm rain-only PDF round-trips + mvr  PASS")


def test_absent_species_guard():
    """A -1 (absent) query must never spuriously match an absent (-1) metadata field."""
    hm = HmMetadata(hydromet_dim=2, iirr=0, iiNr=1, iiPDF_rr=4, iiPDF_Nr=5)
    assert pdf2hydromet_idx(-1, hm) == -1
    assert hydromet2pdf_idx(-1, hm) == -1
    assert rx2Nx_hm_idx(-1, hm) == -1 and Nx2rx_hm_idx(-1, hm) == -1
    assert mvr_hm_max(-1, hm) == 0.0
    print("  absent-species (-1) guard  PASS")


def test_full_species():
    """All four species present (rain/ice/snow/graupel) — every mvr + round-trip branch."""
    hm = HmMetadata(hydromet_dim=8, iirr=0, iiNr=1, iiri=2, iiNi=3,
                    iirs=4, iiNs=5, iirg=6, iiNg=7)
    hm.iiPDF_rr, hm.iiPDF_Nr = 4, 5
    hm.iiPDF_ri, hm.iiPDF_Ni = 6, 7
    hm.iiPDF_rs, hm.iiPDF_Ns = 8, 9
    hm.iiPDF_rg, hm.iiPDF_Ng = 10, 11
    assert mvr_hm_max(0, hm) == mvr_rain_max and mvr_hm_max(2, hm) == mvr_ice_max
    assert mvr_hm_max(4, hm) == mvr_snow_max and mvr_hm_max(6, hm) == mvr_graupel_max
    assert mvr_hm_max(3, hm) == mvr_ice_max   # concentration index maps to same species
    assert hydromet2pdf_idx(6, hm) == 10 and pdf2hydromet_idx(10, hm) == 6
    assert rx2Nx_hm_idx(6, hm) == 7 and Nx2rx_hm_idx(7, hm) == 6
    print("  full four-species mvr + round-trips  PASS")


def test_sclr_idx_fields_mirror_fortran_type():
    """The JAX `SclrIdx` NamedTuple mirrors the Fortran `sclr_idx_type` (array_index.F90) field-for-field, in order — the
    scalar/eddy-scalar index struct (iisclr_rt/thl/CO2, iiedsclr_*). Scalars are sclr_dim=0 in the validated cases, but a
    drifted index would mis-route scalar transport if ever enabled, so the type-mirror is worth pinning. Source-grounded.
    SKIPs if array_index.F90 is absent. (iter 478)"""
    import re
    from clubb_jax.src.CLUBB_core.sclr_idx import SclrIdx
    f90 = os.path.join(_ROOT, "src", "CLUBB_core", "array_index.F90")
    if not os.path.exists(f90):
        raise AssertionError("  SclrIdx field mirror vs Fortran: required input missing (array_index.F90 absent)")
    lines = open(f90).read().splitlines()
    try:
        i0 = next(i for i, l in enumerate(lines) if re.match(r"\s*type\s+sclr_idx_type\s*$", l))
        i1 = next(i for i, l in enumerate(lines) if re.match(r"\s*end\s+type\s+sclr_idx_type", l))
    except StopIteration:
        raise AssertionError("  SclrIdx field mirror: required input missing (type block not found)")
    fort, in_decl = [], False
    for ln in lines[i0:i1]:
        code = ln.split("!", 1)[0]
        if "::" in code:
            in_decl = True
            code = code.split("::", 1)[1]
        elif not in_decl:
            continue
        for tok in code.replace("&", "").split(","):
            t = tok.strip()
            if re.match(r"^[A-Za-z]\w*$", t):
                fort.append(t)
    jax = list(SclrIdx._fields)
    assert fort and fort == jax, (f"SclrIdx fields diverge from sclr_idx_type:\n  Fortran {fort}\n  JAX     {jax}")
    print(f"  SclrIdx fields mirror sclr_idx_type exactly ({len(fort)} fields, same order)  PASS")
