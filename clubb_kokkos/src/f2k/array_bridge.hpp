#pragma once

#include "f2k/derived_type_bridge.hpp"

#include <algorithm>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace f2k {

inline WireField pack_text(const std::string& value, Transfer& transfer, bool copy_values = true,
                           std::size_t minimum_width = 1) {
  // Native output strings receive space to hold ordinary driver metadata.
  const auto width = copy_values ? std::max(minimum_width, value.size())
                                 : std::max<std::size_t>(1024, value.size());
  if (width > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::length_error("f2k text length");
  WireField field{};
  field.allocated = 1;
  field.char_len  = static_cast<int>(width);
  for (auto& extent : field.extent)
    extent = 1;
  auto& bytes = transfer.buffers.emplace_back(width, ' ');
  if (copy_values)
    std::memcpy(bytes.data(), value.data(), value.size());
  field.data = bytes.data();
  return field;
}

inline void unpack_text(std::string& value, const WireField& field) {
  // Discard Fortran's trailing blank padding while preserving leading bytes.
  if (field.rank != 0 || field.char_len < 0 || !field.data)
    throw std::runtime_error("f2k text metadata");
  value.assign(static_cast<const char*>(field.data), field.char_len);
  const auto last = value.find_last_not_of(' ');
  value.resize(last == std::string::npos ? 0 : last + 1);
}

template <class View> void unpack_char_array(View& value, const WireField& field) {
  // Restore native strings into a (character, string-index) View, including padding.
  if (field.rank != 1 || field.char_len < 0)
    throw std::runtime_error("f2k character array metadata");
  value             = View("f2k characters", field.char_len, field.extent[0]);
  auto host         = Kokkos::create_mirror_view(value);
  const auto* bytes = static_cast<const char*>(field.data);
  for (std::int64_t j = 0; j < field.extent[0]; ++j)
    for (int i = 0; i < field.char_len; ++i)
      host(i, j) = bytes[j * field.char_len + i];
  Kokkos::deep_copy(value, host);
}

inline std::size_t field_element_count(const WireField& field) {
  // Reject malformed or oversized extents before allocating a transfer buffer.
  if (field.rank < 0 || field.rank > 4)
    throw std::runtime_error("f2k array rank");
  std::size_t count = 1;
  for (int d = 0; d < field.rank; ++d) {
    if (field.extent[d] < 0 || (count && static_cast<std::uint64_t>(field.extent[d]) >
                                             std::numeric_limits<std::size_t>::max() / count)) {
      throw std::length_error("f2k array extent");
    }
    count *= static_cast<std::size_t>(field.extent[d]);
  }
  return count;
}

template <class T> WireField pack_value(const T& value, Transfer& transfer) {
  // Scalars are copied bitwise into call-owned host storage.
  WireField field{};
  field.allocated = 1;
  for (auto& extent : field.extent)
    extent = 1;
  auto& bytes = transfer.buffers.emplace_back(sizeof(T));
  std::memcpy(bytes.data(), &value, sizeof(T));
  field.data = bytes.data();
  return field;
}

template <class T> void unpack_value(T& value, const WireField& field) {
  // The native wrapper returns an interoperable scalar of the same type.
  if (field.rank != 0 || !field.allocated || !field.data) {
    throw std::runtime_error("f2k missing scalar");
  }
  std::memcpy(&value, field.data, sizeof(T));
}

template <class Array, class View> Array transfer_array(const View& view) {
  // Allocate a contiguous transfer value with the source's logical dimensions.
  // Host/device pack and unpack share this shape conversion; neither retains it.
  static_assert(View::rank >= 1 && View::rank <= 4);
  if constexpr (View::rank == 1)
    return Array(view.label(), view.extent(0));
  else if constexpr (View::rank == 2)
    return Array(view.label(), view.extent(0), view.extent(1));
  else if constexpr (View::rank == 3)
    return Array(view.label(), view.extent(0), view.extent(1), view.extent(2));
  else
    return Array(view.label(), view.extent(0), view.extent(1), view.extent(2), view.extent(3));
}

template <class View>
WireField pack_array(const View& view, Transfer& transfer, bool copy_values = true) {
  // Pack true View indices rather than copying possible allocation padding.
  constexpr int rank = View::rank;
  static_assert(rank >= 1 && rank <= 4);
  WireField field{};
  field.rank      = rank;
  field.allocated = view.is_allocated() ? 1 : 0;
  for (auto& extent : field.extent)
    extent = 1;
  for (int d = 0; d < rank; ++d)
    field.extent[d] = view.extent(d);
  if (!field.allocated || !copy_values)
    return field;
  using T         = typename View::non_const_value_type;
  const auto host = [&] {
    if constexpr (Kokkos::SpaceAccessibility<Kokkos::DefaultHostExecutionSpace,
                                             typename View::memory_space>::accessible) {
      return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, view);
    } else {
      // A cross-space deep_copy cannot read a noncontiguous grid-column
      // subview. Pack on its execution space before transferring to the host.
      using Packed = Kokkos::View<typename View::non_const_data_type, Kokkos::LayoutLeft,
                                  typename View::device_type>;
      using Host =
          Kokkos::View<typename View::non_const_data_type, Kokkos::LayoutLeft, Kokkos::HostSpace>;
      auto result = transfer_array<Host>(view);
      if constexpr (std::is_same_v<typename View::array_layout, Kokkos::LayoutLeft>) {
        if (view.span_is_contiguous()) {
          Kokkos::deep_copy(result, view);
          return result;
        }
      }
      auto packed = transfer_array<Packed>(view);
      Kokkos::deep_copy(packed, view);
      Kokkos::deep_copy(result, packed);
      return result;
    }
  }();
  const auto count = field_element_count(field);
  if (count > std::numeric_limits<std::size_t>::max() / sizeof(T)) {
    throw std::length_error("f2k array bytes");
  }
  auto& bytes     = transfer.buffers.emplace_back(count * sizeof(T));
  std::size_t p   = 0;
  const auto push = [&](const T& value) {
    std::memcpy(bytes.data() + p++ * sizeof(T), &value, sizeof(T));
  };
  for (std::int64_t l = 0; l < (rank >= 4 ? field.extent[3] : 1); ++l)
    for (std::int64_t k = 0; k < (rank >= 3 ? field.extent[2] : 1); ++k)
      for (std::int64_t j = 0; j < (rank >= 2 ? field.extent[1] : 1); ++j)
        for (std::int64_t i = 0; i < field.extent[0]; ++i) {
          if constexpr (rank == 1)
            push(host(i));
          else if constexpr (rank == 2)
            push(host(i, j));
          else if constexpr (rank == 3)
            push(host(i, j, k));
          else
            push(host(i, j, k, l));
        }
  field.data = bytes.empty() ? nullptr : bytes.data();
  return field;
}

template <class View> void unpack_into(View view, const WireField& field, const char* name) {
  // Copy native values into an existing View, including a nonowning column View.
  constexpr int rank = View::rank;
  static_assert(rank >= 1 && rank <= 4);
  if (field.rank != rank || field.char_len != 0)
    throw std::runtime_error("f2k array metadata");
  if (!field.allocated)
    throw std::runtime_error(std::string("f2k missing output: ") + name);
  for (int d = 0; d < rank; ++d)
    if (view.extent(d) != static_cast<std::size_t>(field.extent[d]))
      throw std::runtime_error(std::string("f2k output shape: ") + name);
  if (field_element_count(field) && !field.data)
    throw std::runtime_error("f2k missing array data");
  using Packed = Kokkos::View<typename View::non_const_data_type, Kokkos::LayoutLeft,
                              typename View::device_type>;
  using Host =
      Kokkos::View<typename View::non_const_data_type, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  auto host = [&] {
    if constexpr (Kokkos::SpaceAccessibility<Kokkos::DefaultHostExecutionSpace,
                                             typename View::memory_space>::accessible)
      return Kokkos::create_mirror_view(view);
    else
      return transfer_array<Host>(view);
  }();
  using T         = typename View::non_const_value_type;
  std::size_t p   = 0;
  const auto pull = [&](T& value) {
    std::memcpy(&value, static_cast<const unsigned char*>(field.data) + p++ * sizeof(T), sizeof(T));
  };
  for (std::int64_t l = 0; l < (rank >= 4 ? field.extent[3] : 1); ++l)
    for (std::int64_t k = 0; k < (rank >= 3 ? field.extent[2] : 1); ++k)
      for (std::int64_t j = 0; j < (rank >= 2 ? field.extent[1] : 1); ++j)
        for (std::int64_t i = 0; i < field.extent[0]; ++i) {
          if constexpr (rank == 1)
            pull(host(i));
          else if constexpr (rank == 2)
            pull(host(i, j));
          else if constexpr (rank == 3)
            pull(host(i, j, k));
          else
            pull(host(i, j, k, l));
        }
  if constexpr (Kokkos::SpaceAccessibility<Kokkos::DefaultHostExecutionSpace,
                                           typename View::memory_space>::accessible) {
    Kokkos::deep_copy(view, host);
  } else {
    if constexpr (std::is_same_v<typename View::array_layout, Kokkos::LayoutLeft>) {
      if (view.span_is_contiguous()) {
        Kokkos::deep_copy(view, host);
        return;
      }
    }
    // Restore logical elements on the device so strided destinations retain
    // their allocation and untouched neighboring columns.
    auto packed = transfer_array<Packed>(view);
    Kokkos::deep_copy(packed, host);
    Kokkos::deep_copy(view, packed);
  }
}

template <class View> void unpack_array(View& view, const WireField& field, const char* name) {
  // Allocatable native outputs replace the caller-owned allocation.
  constexpr int rank = View::rank;
  if (!field.allocated) {
    view = View{};
    return;
  }
  if (field.rank != rank || field.char_len != 0)
    throw std::runtime_error("f2k array metadata");
  (void)field_element_count(field);
  if constexpr (rank == 1)
    view = View(name, field.extent[0]);
  else if constexpr (rank == 2)
    view = View(name, field.extent[0], field.extent[1]);
  else if constexpr (rank == 3)
    view = View(name, field.extent[0], field.extent[1], field.extent[2]);
  else
    view = View(name, field.extent[0], field.extent[1], field.extent[2], field.extent[3]);
  unpack_into(view, field, name);
}

} // namespace f2k
