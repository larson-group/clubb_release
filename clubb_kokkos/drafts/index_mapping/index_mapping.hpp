// Map Native PDF and hydrometeor indices and species radius limits.
//
// pdf2hydromet_idx/hydromet2pdf_idx map the two array index sets;
// rx2Nx_hm_idx/Nx2rx_hm_idx map species mixing ratios and concentrations.
// mvr_hm_max returns a maximum mean volume radius [m]. Indices retain Native
// one-based values and unmatched result0; no active-index guard is added.
// These host metadata lookups reuse the generated hm_metadata_type owner.
#pragma once
#include "corr_varnce_module.hpp"

namespace clubb_kokkos {

int pdf2hydromet_idx(int pdf_idx, const hm_metadata_type& hm_metadata);

int hydromet2pdf_idx(int hydromet_idx, const hm_metadata_type& hm_metadata);

int rx2Nx_hm_idx(int rx_idx, const hm_metadata_type& hm_metadata);

int Nx2rx_hm_idx(int Nx_idx, const hm_metadata_type& hm_metadata);

real_t mvr_hm_max(int hydromet_idx, const hm_metadata_type& hm_metadata);

} // namespace clubb_kokkos
