"""Test wrappers for routines from CLUBB_core/index_mapping.F90."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.hm_metadata import HmMetadata


def _make_hm_metadata():
    """Create a compact hydrometeor metadata object with distinct indices."""
    return HmMetadata(
        iirr=1,
        iiri=2,
        iirs=3,
        iirg=4,
        iiNr=5,
        iiNi=6,
        iiNs=7,
        iiNg=8,
        iiPDF_rr=11,
        iiPDF_ri=12,
        iiPDF_rs=13,
        iiPDF_rg=14,
        iiPDF_Nr=15,
        iiPDF_Ni=16,
        iiPDF_Ns=17,
        iiPDF_Ng=18,
        hydromet_tol=np.zeros(8, dtype=np.float64),
    )


def test_pdf_and_hydrometeor_index_maps_are_consistent():
    """pdf2hydromet_idx and hydromet2pdf_idx should invert each defined hydrometeor."""
    hm = _make_hm_metadata()

    pairs = [
        (hm.iiPDF_rr, hm.iirr),
        (hm.iiPDF_ri, hm.iiri),
        (hm.iiPDF_rs, hm.iirs),
        (hm.iiPDF_rg, hm.iirg),
        (hm.iiPDF_Nr, hm.iiNr),
        (hm.iiPDF_Ni, hm.iiNi),
        (hm.iiPDF_Ns, hm.iiNs),
        (hm.iiPDF_Ng, hm.iiNg),
    ]

    for pdf_idx, hydromet_idx in pairs:
        assert clubb_api.pdf2hydromet_idx(pdf_idx, hm) == hydromet_idx
        assert clubb_api.hydromet2pdf_idx(hydromet_idx, hm) == pdf_idx


def test_rx_nx_index_maps_are_consistent():
    """rx2nx_hm_idx and nx2rx_hm_idx should map paired species correctly."""
    hm = _make_hm_metadata()

    assert clubb_api.rx2nx_hm_idx(hm.iirr, hm) == hm.iiNr
    assert clubb_api.rx2nx_hm_idx(hm.iiri, hm) == hm.iiNi
    assert clubb_api.rx2nx_hm_idx(hm.iirs, hm) == hm.iiNs
    assert clubb_api.rx2nx_hm_idx(hm.iirg, hm) == hm.iiNg

    assert clubb_api.nx2rx_hm_idx(hm.iiNr, hm) == hm.iirr
    assert clubb_api.nx2rx_hm_idx(hm.iiNi, hm) == hm.iiri
    assert clubb_api.nx2rx_hm_idx(hm.iiNs, hm) == hm.iirs
    assert clubb_api.nx2rx_hm_idx(hm.iiNg, hm) == hm.iirg


def test_mvr_hm_max_matches_fortran_constants():
    """mvr_hm_max should return the species-specific CLUBB constants."""
    hm = _make_hm_metadata()

    np.testing.assert_allclose(clubb_api.mvr_hm_max(hm.iirr, hm), 5.0e-3)
    np.testing.assert_allclose(clubb_api.mvr_hm_max(hm.iiNi, hm), 1.3e-4)
    np.testing.assert_allclose(clubb_api.mvr_hm_max(hm.iirs, hm), 1.0e-2)
    np.testing.assert_allclose(clubb_api.mvr_hm_max(hm.iiNg, hm), 2.0e-2)
    np.testing.assert_allclose(clubb_api.mvr_hm_max(999, hm), 0.0)
