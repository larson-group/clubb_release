"""Verification of corr_varnce_module — the prescribed PDF-variable correlation arrays."""
import numpy as np


from clubb_jax.src.CLUBB_core.corr_varnce_module import set_corr_arrays_to_default, def_corr_idx, get_corr_var_index, HmMetadata, II_CHI, II_ETA, II_W, II_NCN, II_RR, II_NR

# KK PDF layout [chi, eta, w, Ncn, rr, Nr] -> default-table columns, derived via def_corr_idx
# (mirrors the Fortran set_corr_arrays_to_default <- def_corr_idx chain).
_KK_LAYOUT = HmMetadata(hydromet_dim=2, iiPDF_rr=II_RR, iiPDF_Nr=II_NR)
KK_PDF_TO_DEF = tuple(def_corr_idx(i, _KK_LAYOUT) for i in range(6))

# The 6×6 prescribed lower-triangular arrays, hand-extracted from the Fortran default tables.
# Row index = pdf var j, col index = pdf var i; entry [j,i] (j>i) = corr(j,i).
#                            chi    eta    w     Ncn    rr     Nr
_EXPECTED_CLOUD = np.array([[1.0,  0.0,  0.0,  0.0,  0.0,  0.0],   # chi
                            [-.6,  1.0,  0.0,  0.0,  0.0,  0.0],   # eta
                            [.09,  .027, 1.0,  0.0,  0.0,  0.0],   # w
                            [.09,  .027, .34,  1.0,  0.0,  0.0],   # Ncn
                            [.788, .114, .315, 0.0,  1.0,  0.0],   # rr
                            [.675, .115, .270, 0.0,  .821, 1.0]])  # Nr
_EXPECTED_BELOW = _EXPECTED_CLOUD.copy()
_EXPECTED_BELOW[II_ETA, II_CHI] = 0.3   # only difference: chi-eta -0.6 -> 0.3


def test_set_corr_arrays_to_default():
    """The built 6×6 cloud/below prescribed arrays match the hand-extracted Fortran block."""
    cloud, below = set_corr_arrays_to_default(6, KK_PDF_TO_DEF)
    np.testing.assert_array_equal(cloud, _EXPECTED_CLOUD)
    np.testing.assert_array_equal(below, _EXPECTED_BELOW)
    # cloud and below differ ONLY in chi-eta (no rate-relevant entry).
    diff = np.abs(cloud - below)
    assert np.isclose(diff[II_ETA, II_CHI], 0.9)
    assert np.count_nonzero(diff) == 1   # exactly one entry differs
    print("  set_corr_arrays_to_default: 6×6 cloud/below == Fortran block; differ only in chi-eta  PASS")


def test_def_corr_idx():
    """def_corr_idx maps each PDF variable to its default-table column; -1 for no match."""
    md = HmMetadata(hydromet_dim=2, iiPDF_rr=II_RR, iiPDF_Nr=II_NR)
    assert tuple(def_corr_idx(i, md) for i in range(6)) == (II_CHI, II_ETA, II_W, II_NCN, II_RR, II_NR)
    assert def_corr_idx(II_RR, md) == II_RR and def_corr_idx(II_NR, md) == II_NR
    assert def_corr_idx(99, md) == -1     # no PDF variable at that index
    # A metadata without rain hydrometeors (warm 4-var PDF) maps only chi/eta/w/Ncn.
    md4 = HmMetadata(hydromet_dim=0)
    assert tuple(def_corr_idx(i, md4) for i in range(4)) == (II_CHI, II_ETA, II_W, II_NCN)
    assert def_corr_idx(II_RR, md4) == -1   # rr absent (iiPDF_rr defaults to -1)
    print("  def_corr_idx: KK layout -> (chi,eta,w,Ncn,rr,Nr) cols; absent vars -> -1  PASS")


def test_get_corr_var_index():
    """get_corr_var_index maps a PDF-variable NAME to its iiPDF index; -1 for unknown/absent."""
    md = HmMetadata(hydromet_dim=2, iiPDF_rr=II_RR, iiPDF_Nr=II_NR)
    assert get_corr_var_index("chi", md) == II_CHI and get_corr_var_index("eta", md) == II_ETA
    assert get_corr_var_index("w", md) == II_W and get_corr_var_index("Ncn", md) == II_NCN
    assert get_corr_var_index("rr", md) == II_RR and get_corr_var_index("Nr", md) == II_NR
    assert get_corr_var_index("ri", md) == -1     # frozen species absent in the warm KK PDF
    assert get_corr_var_index("bogus", md) == -1  # unknown name
    print("  get_corr_var_index: name -> iiPDF index (chi..Nr); absent/unknown -> -1  PASS")
