// Shared setup output for the C++ standalone and physics configuration readers.
// write_text mirrors labeled values to stdout and an open case log; write_date
// records the run start. Streams are caller-owned and retain no global file state.
#pragma once
#include "CLUBB_core/array_types.hpp"
#include <array>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <type_traits>

namespace clubb_kokkos {
template <class T, std::size_t N>
std::ostream& operator<<(std::ostream& stream, const std::array<T, N>& values) {
  // Print the fixed Native configuration lists in declaration order.
  for (const auto& value : values)
    stream << value << ' ';
  return stream;
}

template <class Data, class... Properties>
std::ostream& operator<<(std::ostream& stream, const Kokkos::View<Data, Properties...>& value) {
  // Pack logical values before printing a padded or strided column window.
  const auto host = to_host_space(value);
  for (std::size_t i = 0; i < host.size(); ++i)
    stream << host.data()[i] << ' ';
  return stream;
}

inline void write_text(const std::string& text, std::ofstream& file) {
  // Write the same setup message to stdout and the open case log.
  std::cout << text << '\n';
  if (file.is_open())
    file << text << '\n';
}

template <class T> void write_text(const std::string& text, const T& value, std::ofstream& file) {
  // Write one labeled setup value to stdout and the open case log.
  if constexpr (std::is_same_v<T, std::uint8_t>) {
    // Logical View elements use bytes; stream them as booleans, not characters.
    write_text(text, static_cast<bool>(value), file);
  } else {
    std::cout << text << value << '\n';
    if (file.is_open())
      file << text << value << '\n';
  }
}

inline void write_date(std::ofstream& file) {
  // Record the wall-clock start time in the native setup-output stage.
  const auto now = std::time(nullptr);
  std::ostringstream text;
  text << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S");
  write_text(text.str(), file);
}

} // namespace clubb_kokkos
