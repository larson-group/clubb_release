// Declaration-only public Native metadata calls; no bridge implementation here.
// Root must add these to the existing generated f2k owner before promotion and
// remove this staging header. Both calls execute on the host, never in a kernel.
#pragma once
#include "corr_varnce_module.hpp"
#include "clubb_real.hpp"
namespace f2k {
int Nx2rx_hm_idx(int Nx_idx, const clubb_kokkos::hm_metadata_type &hm_metadata);
clubb_kokkos::real_t mvr_hm_max(int hydromet_idx,
                                const clubb_kokkos::hm_metadata_type &hm_metadata);
} // namespace f2k
