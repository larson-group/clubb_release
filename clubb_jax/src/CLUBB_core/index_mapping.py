"""JAX/Python port of `src/CLUBB_core/index_mapping.F90`.

Description:
Functions to map back and forth between the PDF arrays and the hydrometeor
arrays.

The “iiPDF” indices are used to index all PDF variates, including all
hydrometeor variates.
The “ii” indices are used to index hydrometeor arrays.
The “ii” variates are a subset of the “iiPDF” variates.
Conversions between the two sets of indices are done by the
functions pdf2hydromet_idx and hydromet2pdf_idx below.

iiPDF indices:

Included indices:
iiPDF_chi, iiPDF_eta, iiPDF_w, iiPDF_Ncn, iiPDF_rr, & all other hydrometeors

Number of indices:  pdf_dim

Examples of arrays dimensioned by pdf_dim:
mu_x_1_n, corr_array_n_cloud, . . .

Declared as module variables in module array_index

Initialized in subroutine setup_pdf_indices

ii indices:

Included indices:
iirr, iiNr, iiri, iiNi, iirs, iiNs, iirg, iiNg

Number of indices:  hydromet_dim

Examples of arrays dimensioned by hydromet_dim:
hydromet, wphydrometp, . . .

Declared as module variables in module array_index.

Initialized in subroutine init_microphys

References:
  None

Porting deviations:
- JAX/Python indices are 0-based, with -1 denoting an absent hydrometeor/PDF
  variable. The Fortran functions initialize missing results to 0 because 0 is
  never a valid 1-based index.
"""

from clubb_jax.src.CLUBB_core.constants_clubb import (
    mvr_graupel_max,
    mvr_ice_max,
    mvr_rain_max,
    mvr_snow_max,
)


def _hm(hm_metadata, name):
    return getattr(hm_metadata, name, -1)


def pdf2hydromet_idx(pdf_idx, hm_metadata):
    """Returns the hydromet-array index corresponding to ``pdf_idx``.

    Description:
      Returns the position of a specific precipitating hydrometeor corresponding
      to the PDF index (pdf_idx) in the precipitating hydrometeor array
      (hydromet_idx).

    References:
    """
    # Initialize hydromet_idx
    pairs = (
        # Index for rain water mixing ratio, rr.
        ("iiPDF_rr", "iirr"),
        # Index for rain drop concentration, Nr.
        ("iiPDF_Nr", "iiNr"),
        # Index for snow mixing ratio, rs.
        ("iiPDF_rs", "iirs"),
        # Index for snow flake concentration, Ns.
        ("iiPDF_Ns", "iiNs"),
        # Index for graupel mixing ratio, rg.
        ("iiPDF_rg", "iirg"),
        # Index for graupel concentration, Ng.
        ("iiPDF_Ng", "iiNg"),
        # Index for ice mixing ratio, ri.
        ("iiPDF_ri", "iiri"),
        # Index for ice concentration, Ni.
        ("iiPDF_Ni", "iiNi"),
    )
    for pdf_name, hm_name in pairs:
        idx = _hm(hm_metadata, pdf_name)
        if idx >= 0 and pdf_idx == idx:
            return _hm(hm_metadata, hm_name)
    return -1


def hydromet2pdf_idx(hydromet_idx, hm_metadata):
    """Returns the PDF-array index corresponding to ``hydromet_idx``.

    Description:
      Returns the position of a specific precipitating hydrometeor corresponding
      to the precipitating hydrometeor index (hydromet_idx) in the PDF array
      (pdf_idx).

    References:
    """
    # Initialize pdf_idx.
    pairs = (
        # Index for rain water mixing ratio, rr.
        ("iirr", "iiPDF_rr"),
        # Index for rain drop concentration, Nr.
        ("iiNr", "iiPDF_Nr"),
        # Index for ice mixing ratio, ri.
        ("iiri", "iiPDF_ri"),
        # Index for ice concentration, Ni.
        ("iiNi", "iiPDF_Ni"),
        # Index for snow mixing ratio, rs.
        ("iirs", "iiPDF_rs"),
        # Index for snow flake concentration, Ns.
        ("iiNs", "iiPDF_Ns"),
        # Index for graupel mixing ratio, rg.
        ("iirg", "iiPDF_rg"),
        # Index for graupel concentration, Ng.
        ("iiNg", "iiPDF_Ng"),
    )
    for hm_name, pdf_name in pairs:
        idx = _hm(hm_metadata, hm_name)
        if idx >= 0 and hydromet_idx == idx:
            return _hm(hm_metadata, pdf_name)
    return -1


def rx2Nx_hm_idx(rx_idx, hm_metadata):
    """Returns the concentration hydromet index paired with ``rx_idx``.

    Description:
      Returns the position in the hydrometeor array of the specific
      precipitating hydrometeor concentration (Nx_idx) corresponding to the
      precipitating hydrometeor mixing ratio (rx_idx) of the same species of
      precipitating hydrometeor (rain, ice, snow, or graupel).

    References:
    """
    # Initialize Nx_idx.
    pairs = (
        # Index for rain drop concentration, Nr.
        ("iirr", "iiNr"),
        # Index for ice crystal concentration, Ni.
        ("iiri", "iiNi"),
        # Index for snow flake concentration, Ns.
        ("iirs", "iiNs"),
        # Index for graupel concentration, Ng.
        ("iirg", "iiNg"),
    )
    for r_name, n_name in pairs:
        idx = _hm(hm_metadata, r_name)
        if idx >= 0 and rx_idx == idx:
            return _hm(hm_metadata, n_name)
    return -1


def Nx2rx_hm_idx(Nx_idx, hm_metadata):
    """Returns the mixing-ratio hydromet index paired with ``Nx_idx``.

    Description:
      Returns the position in the hydrometeor array of the specific
      precipitating hydrometeor mixing ratio (rx_idx) corresponding to the
      precipitating hydrometeor concentration (Nx_idx) of the same species of
      precipitating hydrometeor (rain, ice, snow, or graupel).

    References:
    """
    # Initialize rx_idx.
    pairs = (
        # Index for rain water mixing ratio, rr.
        ("iiNr", "iirr"),
        # Index for ice mixing ratio, ri.
        ("iiNi", "iiri"),
        # Index for snow mixing ratio, rs.
        ("iiNs", "iirs"),
        # Index for graupel mixing ratio, rg.
        ("iiNg", "iirg"),
    )
    for n_name, r_name in pairs:
        idx = _hm(hm_metadata, n_name)
        if idx >= 0 and Nx_idx == idx:
            return _hm(hm_metadata, r_name)
    return -1


def mvr_hm_max(hydromet_idx, hm_metadata):
    """Return the maximum allowable mean volume radius for ``hydromet_idx``.

    Description:
      Returns the maximum allowable mean volume radius of a specific
      precipitating hydrometeor type (rain, ice, snow, or graupel) corresponding
      to the precipitating hydrometeor index, whether that index is for the
      mixing ratio or concentration associated with that hydrometeor type.

    References:
    """
    # Initialize mvr_hydromet_max.
    species = (
        # Maximum allowable mean volume radius for rain drops.
        (("iirr", "iiNr"), mvr_rain_max),
        # Maximum allowable mean volume radius for ice crystals.
        (("iiri", "iiNi"), mvr_ice_max),
        # Maximum allowable mean volume radius for snow flakes.
        (("iirs", "iiNs"), mvr_snow_max),
        # Maximum allowable mean volume radius for graupel.
        (("iirg", "iiNg"), mvr_graupel_max),
    )
    for (r_name, n_name), mvr_max in species:
        r_idx = _hm(hm_metadata, r_name)
        n_idx = _hm(hm_metadata, n_name)
        if (r_idx >= 0 and hydromet_idx == r_idx) or (
            n_idx >= 0 and hydromet_idx == n_idx
        ):
            return mvr_max
    return 0.0


__all__ = [
    "pdf2hydromet_idx",
    "hydromet2pdf_idx",
    "rx2Nx_hm_idx",
    "Nx2rx_hm_idx",
    "mvr_hm_max",
]
