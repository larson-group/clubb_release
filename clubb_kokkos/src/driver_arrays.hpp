#pragma once

#include "clubb_driver.hpp"
#include "text_writer.hpp"
#include "CLUBB_core/array_types.hpp"
#include <algorithm>
#include <array>
#include <cstdio>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <type_traits>
#include <vector>
#include <string_view>

namespace clubb_kokkos {

inline void append_command_output(const char* command, std::ostream& output) {
  // Run a fixed metadata command. Model paths stay in C++ streams rather than
  // being interpolated into shell commands, preserving filenames with spaces.
  FILE* pipe = ::popen(command, "r");
  if (!pipe)
    throw std::runtime_error("Cannot read CLUBB source metadata");
  std::array<char, 4096> buffer{};
  while (std::fgets(buffer.data(), buffer.size(), pipe))
    output << buffer.data();
  ::pclose(pipe);
}

// These adapters express Fortran value assignment, sequence association and
// character arrays without changing the scientific routine argument order.
inline DriverArray1 host_values(std::initializer_list<real_t> values) {
  // Construct the small array value used by native scalar/array argument association.
  DriverArray1 result("driver values", values.size());
  std::copy(values.begin(), values.end(), result.data());
  return result;
}

template <class View> auto as_columns(const View& source) {
  // View a native initial profile with one leading column, without copying its data.
  static_assert(View::rank == 1 || View::rank == 2);
  if constexpr (View::rank == 1)
    return DriverSlice2(source.data(),
                        Kokkos::LayoutStride(1, 1, source.extent(0), source.stride(0)));
  else
    return DriverSlice3(source.data(),
                        Kokkos::LayoutStride(1, 1, source.extent(0), source.stride(0),
                                             source.extent(1), source.stride(1)));
}

template <class View> std::string text_at(const View& value, int index) {
  // Read one width-first character entry, trimming its native trailing blanks.
  const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, value);
  std::string text(host.extent(0), ' ');
  for (std::size_t i = 0; i < text.size(); ++i)
    text[i] = host(i, index);
  const auto last = text.find_last_not_of(' ');
  text.resize(last == std::string::npos ? 0 : last + 1);
  return text;
}

template <std::size_t N>
std::string text_at(const std::array<std::string_view, N>& value, int index) {
  // Read the same immutable fixed-width metadata from its permanent C++ owner.
  const auto entry = value[index];
  const auto last  = entry.find_last_not_of(' ');
  return std::string(entry.substr(0, last == std::string_view::npos ? 0 : last + 1));
}

inline std::string text_at(const std::vector<std::string>& value, int index) {
  // Input filenames use the same trailing-blank semantics as Native CHARACTER.
  const auto last = value[index].find_last_not_of(' ');
  return value[index].substr(0, last == std::string::npos ? 0 : last + 1);
}

template <std::size_t N>
std::vector<std::string> text_values(const std::array<std::string_view, N>& value) {
  // Keep the declared width/padding when writing parameter-name metadata.
  return std::vector<std::string>(value.begin(), value.end());
}

template <class View> std::vector<std::string> text_values(const View& value) {
  // Fixed native widths also define NetCDF metadata dimensions. Keep padding
  // here; text_at trims only strings used as filenames or console text.
  const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, value);
  std::vector<std::string> strings(host.extent(1), std::string(host.extent(0), ' '));
  for (std::size_t j = 0; j < strings.size(); ++j)
    for (std::size_t i = 0; i < host.extent(0); ++i)
      strings[j][i] = host(i, j);
  return strings;
}

template <class View> auto host_element(const View& value, int index) {
  // Read one host-side diagnostic value from a possibly device-resident View.
  return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, value)(index);
}

template <class View, class T> void set_host_element(const View& value, int index, T element) {
  // Update a diagnostic element while retaining the other columns.
  auto host   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, value);
  host(index) = element;
  Kokkos::deep_copy(value, host);
}

inline grid clone_grid(const grid& source) {
  // Reproduce Fortran derived-type value assignment; View assignment alone aliases arrays.
  grid result = source;
  result.zm   = decltype(result.zm)("grid copy", source.zm.extent(0), source.zm.extent(1));
  Kokkos::deep_copy(result.zm, source.zm);
  result.zt = decltype(result.zt)("grid copy", source.zt.extent(0), source.zt.extent(1));
  Kokkos::deep_copy(result.zt, source.zt);
  result.invrs_dzm = decltype(result.invrs_dzm)("grid copy", source.invrs_dzm.extent(0),
                                                source.invrs_dzm.extent(1));
  Kokkos::deep_copy(result.invrs_dzm, source.invrs_dzm);
  result.invrs_dzt = decltype(result.invrs_dzt)("grid copy", source.invrs_dzt.extent(0),
                                                source.invrs_dzt.extent(1));
  Kokkos::deep_copy(result.invrs_dzt, source.invrs_dzt);
  result.dzm = decltype(result.dzm)("grid copy", source.dzm.extent(0), source.dzm.extent(1));
  Kokkos::deep_copy(result.dzm, source.dzm);
  result.dzt = decltype(result.dzt)("grid copy", source.dzt.extent(0), source.dzt.extent(1));
  Kokkos::deep_copy(result.dzt, source.dzt);
  result.weights_zm2zt = decltype(result.weights_zm2zt)("grid copy", source.weights_zm2zt.extent(0),
                                                        source.weights_zm2zt.extent(1),
                                                        source.weights_zm2zt.extent(2));
  Kokkos::deep_copy(result.weights_zm2zt, source.weights_zm2zt);
  result.weights_zt2zm = decltype(result.weights_zt2zm)("grid copy", source.weights_zt2zm.extent(0),
                                                        source.weights_zt2zm.extent(1),
                                                        source.weights_zt2zm.extent(2));
  Kokkos::deep_copy(result.weights_zt2zm, source.weights_zt2zm);
  return result;
}
} // namespace clubb_kokkos
