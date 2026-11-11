// Column-first arrays for ported CLUBB kernels.
//
// Real/Bool Views use the default Kokkos execution space. The API entry adapters
// below stage host callers once per call, matching native API data-region
// ownership. Callers already using these device Views share their allocations.
#pragma once

#include "clubb_real.hpp"
#include "scratch_workspace.hpp"
#include <Kokkos_Core.hpp>
#include <cstdint>
#include <type_traits>
#include <utility>

namespace clubb_kokkos {

using Real1D = Kokkos::View<real_t*, Kokkos::LayoutLeft>;
using Real2D = Kokkos::View<real_t**, Kokkos::LayoutLeft>;
using Real3D = Kokkos::View<real_t***, Kokkos::LayoutLeft>;
using Real4D = Kokkos::View<real_t****, Kokkos::LayoutLeft>;
using Bool1D = Kokkos::View<std::uint8_t*, Kokkos::LayoutLeft>;

namespace array_detail {
template <class Result, class Source> Result uninitialized_like(const Source& source) {
  // Allocate an owning snapshot/packing field whose complete logical data is copied next.
  const auto allocation = [&] {
    using Value = typename Result::non_const_value_type;
    if constexpr (std::is_arithmetic_v<Value> || std::is_enum_v<Value>)
      return Kokkos::view_alloc(Kokkos::WithoutInitializing, source.label());
    else
      return Kokkos::view_alloc(source.label());
  }();
  if constexpr (Source::rank == 1)
    return Result(allocation, source.extent(0));
  else if constexpr (Source::rank == 2)
    return Result(allocation, source.extent(0), source.extent(1));
  else if constexpr (Source::rank == 3)
    return Result(allocation, source.extent(0), source.extent(1), source.extent(2));
  else {
    static_assert(Source::rank == 4, "Unsupported core snapshot rank");
    return Result(allocation, source.extent(0), source.extent(1), source.extent(2),
                  source.extent(3));
  }
}

template <class Result, class Source> Result wrap_like(const Source& source) {
  // Describe an already-verified dense mapping without owning or copying its storage.
  if constexpr (Source::rank == 1)
    return Result(source.data(), source.extent(0));
  else if constexpr (Source::rank == 2)
    return Result(source.data(), source.extent(0), source.extent(1));
  else if constexpr (Source::rank == 3)
    return Result(source.data(), source.extent(0), source.extent(1), source.extent(2));
  else {
    static_assert(Source::rank == 4, "Unsupported core snapshot rank");
    return Result(source.data(), source.extent(0), source.extent(1), source.extent(2),
                  source.extent(3));
  }
}

template <class Result, class Source>
Result workspace_like(ScratchWorkspace& workspace, const Source& source, bool zero_initialize) {
  // Match logical extents without opening a child frame that would release the result.
  static_assert(Source::rank >= 1 && Source::rank <= 4);
  static_assert(std::is_arithmetic_v<typename Result::non_const_value_type>,
                "Transient workspace staging supports arithmetic Views only");
  return [&]<std::size_t... Axis>(std::index_sequence<Axis...>) {
    if (zero_initialize)
      return workspace.template view<Result>(source.label(), source.extent(Axis)...);
    return workspace.template view<Result>(Kokkos::WithoutInitializing, source.label(),
                                           source.extent(Axis)...);
  }(std::make_index_sequence<Source::rank>{});
}

template <class View> bool dense_left_mapping(const View& source) {
  // Contiguity alone is insufficient: a dense LayoutRight field may have different strides.
  using Dense =
      Kokkos::View<typename View::const_data_type, Kokkos::LayoutLeft, typename View::memory_space>;
  auto dense = wrap_like<Dense>(source);
  if (source.span() != dense.span())
    return false;
  for (int axis = 0; axis < View::rank; ++axis)
    if (source.stride(axis) != dense.stride(axis))
      return false;
  return true;
}

template <class Destination, class Source>
void copy_transient(ScratchWorkspace& workspace, const Destination& destination,
                    const Source& source) {
  // Copy directly when possible; pack each genuinely strided cross-space side.
  // Transfer an unchanged scalar type first, then cast in destination memory,
  // matching the existing mirror-based conversion order when value kinds differ.
  static_assert(Destination::rank == Source::rank);
  static_assert(std::is_arithmetic_v<typename Destination::non_const_value_type> &&
                std::is_arithmetic_v<typename Source::non_const_value_type>);
  static_assert(
      std::is_same_v<typename Destination::data_type, typename Destination::non_const_data_type>);
  if (!source.size() ||
      static_cast<const void*>(destination.data()) == static_cast<const void*>(source.data()))
    return;
  if constexpr (Kokkos::SpaceAccessibility<typename Destination::execution_space,
                                           typename Source::memory_space>::accessible) {
    Kokkos::deep_copy(destination, source);
  } else {
    using DenseSource = Kokkos::View<typename Source::const_data_type, Kokkos::LayoutLeft,
                                     typename Source::memory_space>;
    DenseSource dense_source;
    if (dense_left_mapping(source)) {
      dense_source = wrap_like<DenseSource>(source);
    } else {
      using PackedSource = Kokkos::View<typename Source::non_const_data_type, Kokkos::LayoutLeft,
                                        typename Source::memory_space>;
      auto packed_source = workspace_like<PackedSource>(workspace, source, false);
      Kokkos::deep_copy(packed_source, source);
      dense_source = packed_source;
    }
    if constexpr (std::is_same_v<typename Destination::non_const_value_type,
                                 typename Source::non_const_value_type>) {
      if (dense_left_mapping(destination)) {
        using DenseDestination = Kokkos::View<typename Destination::data_type, Kokkos::LayoutLeft,
                                              typename Destination::memory_space>;
        auto dense_destination = wrap_like<DenseDestination>(destination);
        Kokkos::deep_copy(dense_destination, dense_source);
        return;
      }
    }
    using PackedDestination = Kokkos::View<typename Source::non_const_data_type, Kokkos::LayoutLeft,
                                           typename Destination::memory_space>;
    auto packed_destination = workspace_like<PackedDestination>(workspace, source, false);
    Kokkos::deep_copy(packed_destination, dense_source);
    Kokkos::deep_copy(destination, packed_destination);
  }
}

} // namespace array_detail

template <class View> auto to_host_space(const View& source) {
  // Return an independent, owning snapshot, even for an already-host source.
  // Dense fields need only the final host copy. Pack genuinely strided device
  // windows before crossing memory spaces; their mappings cannot be copied raw.
  using Host =
      Kokkos::View<typename View::non_const_data_type, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  Host host = array_detail::uninitialized_like<Host>(source);
  if (source.size() == 0)
    return host;
  if constexpr (Kokkos::SpaceAccessibility<typename Host::execution_space,
                                           typename View::memory_space>::accessible) {
    Kokkos::deep_copy(host, source);
    return host;
  }
  bool same_mapping = source.span() == host.span();
  for (int axis = 0; axis < View::rank; ++axis)
    same_mapping = same_mapping && source.stride(axis) == host.stride(axis);
  if (same_mapping) {
    using Dense =
        Kokkos::View<typename View::data_type, Kokkos::LayoutLeft, typename View::memory_space>;
    auto dense = array_detail::wrap_like<Dense>(source);
    Kokkos::deep_copy(host, dense);
    return host;
  }
  using Packed  = Kokkos::View<typename View::non_const_data_type, Kokkos::LayoutLeft,
                               typename View::memory_space>;
  Packed packed = array_detail::uninitialized_like<Packed>(source);
  Kokkos::deep_copy(packed, source);
  Kokkos::deep_copy(host, packed);
  return host;
}

template <class View> auto to_execution_space(const View& source, bool copy_values = true) {
  // Alias compatible device storage; otherwise stage at the host-facing API.
  using Result = Kokkos::View<typename View::non_const_data_type, Kokkos::LayoutLeft>;
  if constexpr (std::is_same_v<typename View::memory_space, typename Result::memory_space> &&
                std::is_same_v<typename View::data_type, typename View::non_const_data_type>) {
    if constexpr (std::is_same_v<typename View::array_layout, Kokkos::LayoutLeft>)
      return Result(source);
    else if constexpr (std::is_same_v<typename View::array_layout, Kokkos::LayoutStride> ||
                       View::rank == 1) {
      // A dense subview can share its allocation too; retain Kokkos' ownership tracker.
      if (array_detail::dense_left_mapping(source))
        return Result(source);
    }
  }
  Result result;
  if (copy_values)
    result = array_detail::uninitialized_like<Result>(source);
  else if constexpr (View::rank == 1)
    result = Result(source.label(), source.extent(0));
  else if constexpr (View::rank == 2)
    result = Result(source.label(), source.extent(0), source.extent(1));
  else if constexpr (View::rank == 3)
    result = Result(source.label(), source.extent(0), source.extent(1), source.extent(2));
  else if constexpr (View::rank == 4)
    result = Result(source.label(), source.extent(0), source.extent(1), source.extent(2),
                    source.extent(3));
  else
    static_assert(View::rank <= 4, "Unsupported core API rank");
  if (copy_values) {
    if constexpr (Kokkos::SpaceAccessibility<typename Result::execution_space,
                                             typename View::memory_space>::accessible) {
      Kokkos::deep_copy(result, source);
    } else {
      // Pack strided host windows before crossing memory spaces. A direct
      // host/device deep_copy requires matching contiguous layouts.
      auto host = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, result);
      Kokkos::deep_copy(host, source);
      Kokkos::deep_copy(result, host);
    }
  }
  return result;
}

template <class View, class ExecutionView>
void from_execution_space(const View& destination, const ExecutionView& source) {
  // Restore native InOut/Out values while retaining the caller's allocation.
  if (source.size() && destination.data() != source.data()) {
    if constexpr (Kokkos::SpaceAccessibility<typename View::execution_space,
                                             typename ExecutionView::memory_space>::accessible) {
      Kokkos::deep_copy(destination, source);
    } else {
      auto host = Kokkos::create_mirror_view_and_copy(typename View::memory_space{}, source);
      Kokkos::deep_copy(destination, host);
    }
  }
}

// Match create_mirror_view's no-copy host alias, or reserve an uninitialized
// host mirror in the caller frame when the source needs separate host storage.
template <class View> auto create_host_mirror(ScratchWorkspace& workspace, const View& source) {
  // These mirrors are temporary LayoutLeft fields, not independent snapshots.
  // Aliasing a CPU source preserves partly written arrays without reading their
  // uninitialized cells into a new buffer. The caller fills separate GPU mirrors.
  static_assert(std::is_same_v<typename View::array_layout, Kokkos::LayoutLeft>);
  using Host = decltype(Kokkos::create_mirror_view(Kokkos::WithoutInitializing, source));
  if constexpr (std::is_same_v<typename View::memory_space, typename Host::memory_space> &&
                std::is_same_v<typename View::data_type, typename Host::data_type>)
    return Host(source);
  else
    return array_detail::workspace_like<Host>(workspace, source, false);
}

// Explicit workspace overloads below return caller-frame temporaries. The
// original owning helpers above remain available for snapshots and UDT fields.
template <class View> auto to_host_space(ScratchWorkspace& workspace, const View& source) {
  // Copy an independent temporary snapshot into the active caller's host frame.
  using Host =
      Kokkos::View<typename View::non_const_data_type, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  auto host = array_detail::workspace_like<Host>(workspace, source, false);
  array_detail::copy_transient(workspace, host, source);
  return host;
}

template <class View>
auto to_execution_space(ScratchWorkspace& workspace, const View& source, bool copy_values = true) {
  // Keep compatible alias semantics; incompatible Out fields still start at zero.
  using Result = Kokkos::View<typename View::non_const_data_type, Kokkos::LayoutLeft>;
  if constexpr (std::is_same_v<typename View::memory_space, typename Result::memory_space> &&
                std::is_same_v<typename View::data_type, typename View::non_const_data_type>) {
    if constexpr (std::is_same_v<typename View::array_layout, Kokkos::LayoutLeft>)
      return Result(source);
    else if constexpr (std::is_same_v<typename View::array_layout, Kokkos::LayoutStride> ||
                       View::rank == 1) {
      // A dense subview can share its allocation too; retain Kokkos' ownership tracker.
      if (array_detail::dense_left_mapping(source))
        return Result(source);
    }
  }
  auto result = array_detail::workspace_like<Result>(workspace, source, !copy_values);
  if (copy_values)
    array_detail::copy_transient(workspace, result, source);
  return result;
}

template <class View>
auto to_execution_space(ScratchWorkspace& workspace, decltype(Kokkos::WithoutInitializing),
                        const View& source) {
  // Reserve a fully overwritten temporary Out field without zeroing or copying.
  // Compatible mutable execution storage still aliases its existing allocation.
  using Result = Kokkos::View<typename View::non_const_data_type, Kokkos::LayoutLeft>;
  if constexpr (std::is_same_v<typename View::memory_space, typename Result::memory_space> &&
                std::is_same_v<typename View::data_type, typename View::non_const_data_type>) {
    if constexpr (std::is_same_v<typename View::array_layout, Kokkos::LayoutLeft>)
      return Result(source);
    else if constexpr (std::is_same_v<typename View::array_layout, Kokkos::LayoutStride> ||
                       View::rank == 1) {
      // A dense subview can share its allocation too; retain Kokkos' ownership tracker.
      if (array_detail::dense_left_mapping(source))
        return Result(source);
    }
  }
  return array_detail::workspace_like<Result>(workspace, source, false);
}

template <class View, class ExecutionView>
void from_execution_space(ScratchWorkspace& workspace, const View& destination,
                          const ExecutionView& source) {
  // Restore data into the existing destination; only cross-space packing is pooled.
  array_detail::copy_transient(workspace, destination, source);
}

} // namespace clubb_kokkos
