"""User-facing wrappers for routines from CLUBB_core/index_mapping.F90."""

from clubb_python.derived_types.hm_metadata import HmMetadata


MVR_RAIN_MAX = 5.0e-3
MVR_ICE_MAX = 1.3e-4
MVR_SNOW_MAX = 1.0e-2
MVR_GRAUPEL_MAX = 2.0e-2


def pdf2hydromet_idx(pdf_idx: int, hm_metadata: HmMetadata) -> int:
    """Map a PDF-array hydrometeor index into the hydrometeor-array index."""
    if pdf_idx == hm_metadata.iiPDF_rr:
        return hm_metadata.iirr
    if pdf_idx == hm_metadata.iiPDF_Nr:
        return hm_metadata.iiNr
    if pdf_idx == hm_metadata.iiPDF_rs:
        return hm_metadata.iirs
    if pdf_idx == hm_metadata.iiPDF_Ns:
        return hm_metadata.iiNs
    if pdf_idx == hm_metadata.iiPDF_rg:
        return hm_metadata.iirg
    if pdf_idx == hm_metadata.iiPDF_Ng:
        return hm_metadata.iiNg
    if pdf_idx == hm_metadata.iiPDF_ri:
        return hm_metadata.iiri
    if pdf_idx == hm_metadata.iiPDF_Ni:
        return hm_metadata.iiNi
    return 0


def hydromet2pdf_idx(hydromet_idx: int, hm_metadata: HmMetadata) -> int:
    """Map a hydrometeor-array index into the matching PDF-array hydrometeor index."""
    if hydromet_idx == hm_metadata.iirr:
        return hm_metadata.iiPDF_rr
    if hydromet_idx == hm_metadata.iiNr:
        return hm_metadata.iiPDF_Nr
    if hydromet_idx == hm_metadata.iiri:
        return hm_metadata.iiPDF_ri
    if hydromet_idx == hm_metadata.iiNi:
        return hm_metadata.iiPDF_Ni
    if hydromet_idx == hm_metadata.iirs:
        return hm_metadata.iiPDF_rs
    if hydromet_idx == hm_metadata.iiNs:
        return hm_metadata.iiPDF_Ns
    if hydromet_idx == hm_metadata.iirg:
        return hm_metadata.iiPDF_rg
    if hydromet_idx == hm_metadata.iiNg:
        return hm_metadata.iiPDF_Ng
    return 0


def rx2nx_hm_idx(rx_idx: int, hm_metadata: HmMetadata) -> int:
    """Map a hydrometeor mixing-ratio index to the paired concentration index."""
    if rx_idx == hm_metadata.iirr:
        return hm_metadata.iiNr
    if rx_idx == hm_metadata.iiri:
        return hm_metadata.iiNi
    if rx_idx == hm_metadata.iirs:
        return hm_metadata.iiNs
    if rx_idx == hm_metadata.iirg:
        return hm_metadata.iiNg
    return 0


def nx2rx_hm_idx(nx_idx: int, hm_metadata: HmMetadata) -> int:
    """Map a hydrometeor concentration index to the paired mixing-ratio index."""
    if nx_idx == hm_metadata.iiNr:
        return hm_metadata.iirr
    if nx_idx == hm_metadata.iiNi:
        return hm_metadata.iiri
    if nx_idx == hm_metadata.iiNs:
        return hm_metadata.iirs
    if nx_idx == hm_metadata.iiNg:
        return hm_metadata.iirg
    return 0


def mvr_hm_max(hydromet_idx: int, hm_metadata: HmMetadata) -> float:
    """Return the maximum allowable mean-volume radius for a hydrometeor species."""
    if hydromet_idx == hm_metadata.iirr or hydromet_idx == hm_metadata.iiNr:
        return MVR_RAIN_MAX
    if hydromet_idx == hm_metadata.iiri or hydromet_idx == hm_metadata.iiNi:
        return MVR_ICE_MAX
    if hydromet_idx == hm_metadata.iirs or hydromet_idx == hm_metadata.iiNs:
        return MVR_SNOW_MAX
    if hydromet_idx == hm_metadata.iirg or hydromet_idx == hm_metadata.iiNg:
        return MVR_GRAUPEL_MAX
    return 0.0
