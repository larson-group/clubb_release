// Host formatting and comparison reductions used by Native microphysics diagnostics.
#pragma once
#include "CLUBB_core/array_types.hpp"
#include <string>
#include <sstream>
#include <limits>
#include <iomanip>
namespace clubb_kokkos {
namespace microphys_detail {
template <class View> std::string dump_values(View values) {
  // Format Native arrays in first-index-fastest order at a host diagnostic boundary.
  auto host = to_host_space(values);
  std::ostringstream stream;
  stream << std::setprecision(std::numeric_limits<real_t>::max_digits10);
  for (std::size_t j = 0; j < host.size(); ++j)
    stream << host.data()[j] << ' ';
  return stream.str();
}
template <class First, class Second> bool any_less(First first, Second second) {
  // Reduce Native ANY comparisons after staging logical indices on the backend.
  int count = 0;
  if constexpr (Kokkos::is_view<Second>::value) {
    Kokkos::parallel_reduce(
        "microphys::any_less", first.size(),
        KOKKOS_LAMBDA(int p, int& c) {
          const int i = p % first.extent(0), k = p / first.extent(0);
          if (first(i, k) < second(i, k))
            ++c;
        },
        count);
  } else {
    Kokkos::parallel_reduce(
        "microphys::any_less", first.size(),
        KOKKOS_LAMBDA(int p, int& c) {
          const int i = p % first.extent(0), k = p / first.extent(0);
          if (first(i, k) < second)
            ++c;
        },
        count);
  }
  return count != 0;
}
} // namespace microphys_detail
namespace microphys_detail {
inline std::string trim(std::string text) {
  // Native TRIM removes trailing blanks from fixed-width diagnostic names.
  const auto last = text.find_last_not_of(' ');
  text.resize(last == std::string::npos ? 0 : last + 1);
  return text;
}
} // namespace microphys_detail

} // namespace clubb_kokkos
