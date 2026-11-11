// CLUBB status codes and host diagnostic-level accessors.
// clubb_at_least_debug_level_api tests the current thread's threshold;
// set_clubb_debug_level_api replaces it without clamping. The initial level is 0.
// any_equal reduces a status View where it resides and returns one host boolean.
#pragma once

#include <Kokkos_Core.hpp>

// Model-wide status codes from src/CLUBB_core/error_code.F90.
namespace clubb_kokkos {

inline constexpr int clubb_no_error                 = 0;
inline constexpr int clubb_generalized_grd_test_err = 50;
inline constexpr int clubb_fatal_error              = 99;

bool clubb_at_least_debug_level_api(int level); // In
void set_clubb_debug_level_api(int level);     // In


template <class View, class T> bool any_equal(const View& value, T target) {
  // Preserve native ANY semantics without staging the complete status array.
  // Integer OR is independent of execution order, including an empty View.
  using Execution = typename View::execution_space;
  int found       = 0;
  Kokkos::parallel_reduce(
      "error_code::any_equal", Kokkos::RangePolicy<Execution>(0, value.extent(0)),
      KOKKOS_LAMBDA(const std::size_t i, int& result) { result |= value(i) == target; },
      Kokkos::BOr<int>(found));
  return found != 0;
}

} // namespace clubb_kokkos
