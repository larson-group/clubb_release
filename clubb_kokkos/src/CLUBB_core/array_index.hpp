#pragma once

#include <Kokkos_Core.hpp>
#include <cstdint>
#include <string>
#include "clubb_real.hpp"

namespace clubb_kokkos {

// Mirrors src/CLUBB_core/array_index.F90; arrays use Fortran-order LayoutLeft.
struct sclr_idx_type {
  int iisclr_rt    = -1;
  int iisclr_thl   = -1;
  int iisclr_CO2   = -1;
  int iiedsclr_rt  = -1;
  int iiedsclr_thl = -1;
  int iiedsclr_CO2 = -1;
};

}  // namespace clubb_kokkos
