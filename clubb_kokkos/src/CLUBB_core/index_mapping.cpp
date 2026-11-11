// Map Native PDF and hydrometeor indices and species radius limits.
//
// pdf2hydromet_idx/hydromet2pdf_idx map the two array index sets;
// rx2Nx_hm_idx/Nx2rx_hm_idx map species mixing ratios and concentrations.
// mvr_hm_max returns a maximum mean volume radius [m]. Indices use zero-based
// slots and unmatched result -1; no active-index guard is added.
// These host metadata lookups reuse the generated hm_metadata_type owner.
// Description:
// Functions to map back and forth between the PDF arrays and the hydrometeor
// arrays.

// The “iiPDF” indices are used to index all PDF variates, including all
// hydrometeor variates.
// The “ii” indices are used to index hydrometeor arrays.
// The “ii” variates are a subset of the “iiPDF” variates.
// Conversions between the two sets of indices are done by the
// functions pdf2hydromet_idx and hydromet2pdf_idx below.
//
// ------------------------------------------------------------------------
//
// iiPDF indices:
//
// Included indices:
// iiPDF_chi, iiPDF_eta, iiPDF_w, iiPDF_Ncn, iiPDF_rr, & all other hydrometeors
//
// Number of indices:  pdf_dim
//
// Examples of arrays dimensioned by pdf_dim:
// mu_x_1_n, corr_array_n_cloud, . . .
//
// Stored in hm_metadata_type in corr_varnce_module.
//
// Initialized in init_pdf_hydromet_arrays_api.
//
// ----------------------------------------------------------------------
//
// ii indices:
//
// Included indices:
// iirr, iiNr, iiri, iiNi, iirs, iiNs, iirg, iiNg
//
// Number of indices:  hydromet_dim
//
// Examples of arrays dimensioned by hydromet_dim:
// hydromet, wphydrometp, . . .
//
// Stored in hm_metadata_type in corr_varnce_module.
//
// Initialized in init_pdf_hydromet_arrays_api.
//
// -----------------------------------------------------------------------
//
// References:
//   None
//-------------------------------------------------------------------------

#include "index_mapping.hpp"
#include "constants_clubb.hpp"

namespace clubb_kokkos {

//=============================================================================
int pdf2hydromet_idx(int pdf_idx, const hm_metadata_type& hm_metadata) {
  // Description:
  // Returns the position of a specific precipitating hydrometeor corresponding
  // to the PDF index (pdf_idx) in the precipitating hydrometeor array
  // (hydromet_idx).

  // References:
  //-----------------------------------------------------------------------

  int hydromet_idx; // Zero-based slot, -1 when unmatched.

  // Initialize hydromet_idx
  hydromet_idx = -1;

  if (pdf_idx == hm_metadata.iiPDF_rr) {

    // Index for rain water mixing ratio, rr.
    hydromet_idx = hm_metadata.iirr;

  } else if (pdf_idx == hm_metadata.iiPDF_Nr) {

    // Index for rain drop concentration, Nr.
    hydromet_idx = hm_metadata.iiNr;

  } else if (pdf_idx == hm_metadata.iiPDF_rs) {

    // Index for snow mixing ratio, rs.
    hydromet_idx = hm_metadata.iirs;

  } else if (pdf_idx == hm_metadata.iiPDF_Ns) {

    // Index for snow flake concentration, Ns.
    hydromet_idx = hm_metadata.iiNs;

  } else if (pdf_idx == hm_metadata.iiPDF_rg) {

    // Index for graupel mixing ratio, rg.
    hydromet_idx = hm_metadata.iirg;

  } else if (pdf_idx == hm_metadata.iiPDF_Ng) {

    // Index for graupel concentration, Ng.
    hydromet_idx = hm_metadata.iiNg;

  } else if (pdf_idx == hm_metadata.iiPDF_ri) {

    // Index for ice mixing ratio, ri.
    hydromet_idx = hm_metadata.iiri;

  } else if (pdf_idx == hm_metadata.iiPDF_Ni) {

    // Index for ice concentration, Ni.
    hydromet_idx = hm_metadata.iiNi;
  }

  return hydromet_idx;
}

//=============================================================================
int hydromet2pdf_idx(int hydromet_idx, const hm_metadata_type& hm_metadata) {
  // Description:
  // Returns the position of a specific precipitating hydrometeor corresponding
  // to the precipitating hydrometeor index (hydromet_idx) in the PDF array
  // (pdf_idx).

  // References:
  //-----------------------------------------------------------------------

  int pdf_idx; // Zero-based slot, -1 when unmatched.

  // Initialize pdf_idx.
  pdf_idx = -1;

  if (hydromet_idx == hm_metadata.iirr) {

    // Index for rain water mixing ratio, rr.
    pdf_idx = hm_metadata.iiPDF_rr;

  } else if (hydromet_idx == hm_metadata.iiNr) {

    // Index for rain drop concentration, Nr.
    pdf_idx = hm_metadata.iiPDF_Nr;

  } else if (hydromet_idx == hm_metadata.iiri) {

    // Index for ice mixing ratio, ri.
    pdf_idx = hm_metadata.iiPDF_ri;

  } else if (hydromet_idx == hm_metadata.iiNi) {

    // Index for ice concentration, Ni.
    pdf_idx = hm_metadata.iiPDF_Ni;

  } else if (hydromet_idx == hm_metadata.iirs) {

    // Index for snow mixing ratio, rs.
    pdf_idx = hm_metadata.iiPDF_rs;

  } else if (hydromet_idx == hm_metadata.iiNs) {

    // Index for snow flake concentration, Ns.
    pdf_idx = hm_metadata.iiPDF_Ns;

  } else if (hydromet_idx == hm_metadata.iirg) {

    // Index for graupel mixing ratio, rg.
    pdf_idx = hm_metadata.iiPDF_rg;

  } else if (hydromet_idx == hm_metadata.iiNg) {

    // Index for graupel concentration, Ng.
    pdf_idx = hm_metadata.iiPDF_Ng;
  }

  return pdf_idx;
}

//=============================================================================
int rx2Nx_hm_idx(int rx_idx, const hm_metadata_type& hm_metadata) {
  // Description:
  // Returns the position in the hydrometeor array of the specific
  // precipitating hydrometeor concentration (Nx_idx) corresponding to the
  // precipitating hydrometeor mixing ratio (rx_idx) of the same species of
  // precipitating hydrometeor (rain, ice, snow, or graupel).

  // References:
  //-----------------------------------------------------------------------

  int Nx_idx; // Zero-based slot, -1 when unmatched.

  // Initialize Nx_idx.
  Nx_idx = -1;

  if (rx_idx == hm_metadata.iirr) {

    // Index for rain drop concentration, Nr.
    Nx_idx = hm_metadata.iiNr;

  } else if (rx_idx == hm_metadata.iiri) {

    // Index for ice crystal concentration, Ni.
    Nx_idx = hm_metadata.iiNi;

  } else if (rx_idx == hm_metadata.iirs) {

    // Index for snow flake concentration, Ns.
    Nx_idx = hm_metadata.iiNs;

  } else if (rx_idx == hm_metadata.iirg) {

    // Index for graupel concentration, Ng.
    Nx_idx = hm_metadata.iiNg;
  }

  return Nx_idx;
}

//=============================================================================
int Nx2rx_hm_idx(int Nx_idx, const hm_metadata_type& hm_metadata) {
  // Description:
  // Returns the position in the hydrometeor array of the specific
  // precipitating hydrometeor mixing ratio (rx_idx) corresponding to the
  // precipitating hydrometeor concentration (Nx_idx) of the same species of
  // precipitating hydrometeor (rain, ice, snow, or graupel).

  // References:
  //-----------------------------------------------------------------------

  int rx_idx; // Zero-based slot, -1 when unmatched.

  // Initialize rx_idx.
  rx_idx = -1;

  if (Nx_idx == hm_metadata.iiNr) {

    // Index for rain water mixing ratio, rr.
    rx_idx = hm_metadata.iirr;

  } else if (Nx_idx == hm_metadata.iiNi) {

    // Index for ice mixing ratio, ri.
    rx_idx = hm_metadata.iiri;

  } else if (Nx_idx == hm_metadata.iiNs) {

    // Index for snow mixing ratio, rs.
    rx_idx = hm_metadata.iirs;

  } else if (Nx_idx == hm_metadata.iiNg) {

    // Index for graupel mixing ratio, rg.
    rx_idx = hm_metadata.iirg;
  }

  return rx_idx;
}

//=============================================================================
real_t mvr_hm_max(int hydromet_idx, const hm_metadata_type& hm_metadata) {
  // Description:
  // Returns the maximum allowable mean volume radius of a specific
  // precipitating hydrometeor type (rain, ice, snow, or graupel) corresponding
  // to the precipitating hydrometeor index, whether that index is for the
  // mixing ratio or concentration associated with that hydrometeor type.

  // References:
  //-----------------------------------------------------------------------

  real_t mvr_hydromet_max; // Native result: maximum mean volume radius [m].

  // Initialize mvr_hydromet_max.
  mvr_hydromet_max = zero;

  if (hydromet_idx == hm_metadata.iirr || hydromet_idx == hm_metadata.iiNr) {

    // Maximum allowable mean volume radius for rain drops.
    mvr_hydromet_max = mvr_rain_max;

  } else if (hydromet_idx == hm_metadata.iiri || hydromet_idx == hm_metadata.iiNi) {

    // Maximum allowable mean volume radius for ice crystals.
    mvr_hydromet_max = mvr_ice_max;

  } else if (hydromet_idx == hm_metadata.iirs || hydromet_idx == hm_metadata.iiNs) {

    // Maximum allowable mean volume radius for snow flakes.
    mvr_hydromet_max = mvr_snow_max;

  } else if (hydromet_idx == hm_metadata.iirg || hydromet_idx == hm_metadata.iiNg) {

    // Maximum allowable mean volume radius for graupel.
    mvr_hydromet_max = mvr_graupel_max;
  }

  return mvr_hydromet_max;
}

} // namespace clubb_kokkos
