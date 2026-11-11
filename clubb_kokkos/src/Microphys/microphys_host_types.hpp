// Host C++ engine staging types. Model fields remain caller-owned Kokkos Views.
#pragma once
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/T_in_K_module.hpp"
#include "CLUBB_core/saturation.hpp"
#include "CLUBB_core/model_flags.hpp"
#include "CLUBB_core/Skx_module.hpp"
#include "Microphys/parameters_microphys.hpp"
#include "Microphys/KK_microphys/parameters_KK.hpp"
#include "Microphys/Morrison_microphys/module_mp_graupel_state.hpp"
#include "Microphys/microphys_diagnostics.hpp"
#include "CLUBB_core/stats_netcdf.hpp"
#include <cmath>
namespace clubb_kokkos {
struct KK_microphys_adj_terms_type {
  real_t rrm_src_adj  = zero; // Total rain source adjustment [(kg/kg)/s].
  real_t rrm_evap_adj = zero; // Total rain evaporation adjustment [(kg/kg)/s].
  real_t Nrm_src_adj  = zero; // Total number source adjustment [(num/kg)/s].
  real_t Nrm_evap_adj = zero; // Total number evaporation adjustment [(num/kg)/s].
};
using HostReal1D  = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using HostFloat1D = Kokkos::View<float*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using HostInt1D   = Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using HostBool1D  = Kokkos::View<std::uint8_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using HostReal2D  = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using HostFloat2D = Kokkos::View<float**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using HostInt2D   = Kokkos::View<int**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using HostBool2D  = Kokkos::View<std::uint8_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using HostReal3D  = Kokkos::View<real_t***, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using HostFloat3D = Kokkos::View<float***, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using HostInt3D   = Kokkos::View<int***, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using HostBool3D  = Kokkos::View<std::uint8_t***, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using HostReal4D  = Kokkos::View<real_t****, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using HostFloat4D = Kokkos::View<float****, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using HostInt4D   = Kokkos::View<int****, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using HostBool4D  = Kokkos::View<std::uint8_t****, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using HostAdj2D =
    Kokkos::View<KK_microphys_adj_terms_type**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
namespace microphys_host_detail {
template <class View> bool any_true(View values) {
  // Evaluate Native ANY in first-index traversal order, honoring View strides.
  if constexpr (View::rank == 1) {
    for (std::size_t i = 0; i < values.extent(0); ++i)
      if (values(i))
        return true;
  } else {
    static_assert(View::rank == 2);
    for (std::size_t k = 0; k < values.extent(1); ++k)
      for (std::size_t i = 0; i < values.extent(0); ++i)
        if (values(i, k))
          return true;
  }
  return false;
}
template <class View> HostReal1D maxval_levels(View values) {
  // Evaluate Native MAXVAL(dim=2), preserving vertical comparison order.
  HostReal1D result("maxval_levels", values.extent(0));
  for (std::size_t i = 0; i < values.extent(0); ++i) {
    result(i) = values(i, 0);
    for (std::size_t k = 1; k < values.extent(1); ++k)
      if (values(i, k) > result(i))
        result(i) = values(i, k);
  }
  return result;
}
} // namespace microphys_host_detail
} // namespace clubb_kokkos
