// Compile-only interface edge: constant, strided HostSpace inputs must return
// an owned LayoutLeft HostSpace result. This translation unit is never run.
#include "remapping_module.hpp"
#include <type_traits>

namespace clubb_kokkos {
using ConstHost2D = Kokkos::View<const real_t**, Kokkos::LayoutStride, Kokkos::HostSpace>;
using ConstHost1D = Kokkos::View<const real_t*, Kokkos::LayoutStride, Kokkos::HostSpace>;
using OwnedHost2D = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;

OwnedHost2D compile_const_strided_interface(
    int ngrdcol, const grid& gr_source, const grid& gr_target,
    int source_values_idx, ConstHost2D source_values, int target_values_idx,
    int total_idx_rho_lin_spline, ConstHost2D rho_lin_spline_vals,
    ConstHost2D rho_lin_spline_levels, int iv, ConstHost1D p_sfc,
    int grid_remap_method, bool l_zt_variable) {
  auto result = remap_vals_to_target(
      ngrdcol, gr_source, gr_target, source_values_idx, source_values,
      target_values_idx, total_idx_rho_lin_spline, rho_lin_spline_vals,
      rho_lin_spline_levels, iv, p_sfc, grid_remap_method, l_zt_variable);
  static_assert(std::is_same_v<decltype(result), OwnedHost2D>);
  static_assert(decltype(result)::memory_traits::is_unmanaged == 0);
  return result;
}
} // namespace clubb_kokkos
