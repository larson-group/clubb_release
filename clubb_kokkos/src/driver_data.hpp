// Keep a host field and its execution-space allocation for a driver case.
//
// bind attaches an existing mutable HostSpace View. sync_device and sync_host
// copy only when the other side is newer; modify_host and modify_device record
// the side written by the caller. host_view and execution_view return handles
// without copying. release detaches both allocations for case cleanup.
//
// Before a partial or InOut write, synchronize the side being written. Marking
// a side modified makes it authoritative, so a complete Out overwrite may
// discard the other side. Independently modifying both sides loses data; this
// class does not merge changes. Callers must finish outstanding kernels before
// rebinding or releasing a field. Compatible CPU Views share the host allocation.
#pragma once

#include <Kokkos_Core.hpp>
#include <stdexcept>
#include <type_traits>

namespace clubb_kokkos {

template <class HostView> class HostDeviceField {
public:
  using execution_view_type =
      Kokkos::View<typename HostView::non_const_data_type, Kokkos::LayoutLeft>;

  HostDeviceField()                                  = default;
  HostDeviceField(const HostDeviceField&)            = delete;
  HostDeviceField& operator=(const HostDeviceField&) = delete;
  HostDeviceField(HostDeviceField&&)                 = default;
  HostDeviceField& operator=(HostDeviceField&&)      = default;

  void bind(const HostView& source) {
    // A new host window becomes authoritative. Reuse execution storage when
    // extents match, including when a later parameter batch changes the window.
    if (bound_ && same_mapping(host_, source))
      return;
    const bool same_shape = bound_ && same_extents(host_, source);
    host_                 = source;
    bound_                = true;
    if constexpr (aliases_host) {
      execution_ = execution_view_type(source);
      current_   = Current::both;
    } else {
      if (!same_shape)
        execution_ = allocate<execution_view_type>(source);
      if constexpr (!execution_can_access_host) {
        pack_ = true;
        if constexpr (std::is_same_v<typename HostView::array_layout, Kokkos::LayoutLeft>)
          pack_ = !same_mapping_without_pointer(source, execution_);
        if (pack_ && (!same_shape || packed_.data() == nullptr))
          packed_ = Kokkos::create_mirror_view(execution_);
      }
      current_ = Current::host;
    }
  }

  void sync_device() {
    // Transfer host changes once, retaining both allocations and any packing
    // buffer. Blocking deep_copy completes data movement before subsequent use.
    require_bound();
    if (current_ != Current::host)
      return;
    if constexpr (!aliases_host) {
      if constexpr (execution_can_access_host) {
        Kokkos::deep_copy(execution_, host_);
      } else {
        if (pack_) {
          Kokkos::deep_copy(packed_, host_);
          Kokkos::deep_copy(execution_, packed_);
        } else {
          Kokkos::deep_copy(execution_, host_);
        }
      }
    }
    current_ = Current::both;
  }

  void sync_host() {
    // Restore execution-space changes before host reads or partial host writes.
    require_bound();
    if (current_ != Current::device)
      return;
    if constexpr (!aliases_host) {
      if constexpr (execution_can_access_host) {
        Kokkos::deep_copy(host_, execution_);
      } else {
        if (pack_) {
          Kokkos::deep_copy(packed_, execution_);
          Kokkos::deep_copy(host_, packed_);
        } else {
          Kokkos::deep_copy(host_, execution_);
        }
      }
    }
    current_ = Current::both;
  }

  void modify_host() {
    // Record a host write; sync_host must precede partial writes to stale data.
    require_bound();
    current_ = Current::host;
  }

  void modify_device() {
    // Record an execution-space write, without fencing or transferring data.
    require_bound();
    current_ = Current::device;
  }

  const HostView& host_view() const {
    // Return the bound host handle without synchronizing it.
    require_bound();
    return host_;
  }

  const execution_view_type& execution_view() const {
    // Return the persistent execution handle without synchronizing it.
    require_bound();
    return execution_;
  }

  void release() {
    // Detach storage after the case has finished all kernels using the field.
    host_      = HostView{};
    execution_ = execution_view_type{};
    packed_    = packed_view_type{};
    current_   = Current::both;
    bound_     = false;
    pack_      = false;
  }

private:
  static_assert(HostView::rank >= 1 && HostView::rank <= 4);
  static_assert(std::is_same_v<typename HostView::memory_space, Kokkos::HostSpace>);
  static_assert(
      std::is_same_v<typename HostView::data_type, typename HostView::non_const_data_type>);
  using packed_view_type =
      Kokkos::View<typename HostView::non_const_data_type, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  static constexpr bool aliases_host =
      std::is_same_v<typename execution_view_type::memory_space, Kokkos::HostSpace> &&
      std::is_same_v<typename HostView::array_layout, Kokkos::LayoutLeft>;
  static constexpr bool execution_can_access_host =
      Kokkos::SpaceAccessibility<typename execution_view_type::execution_space,
                                 Kokkos::HostSpace>::accessible;
  enum class Current { both, host, device };

  template <class OtherView>
  static bool same_extents(const HostView& source, const OtherView& destination) {
    // Compare dimensions without assuming contiguous windows or nonzero sizes.
    for (int i = 0; i < HostView::rank; ++i)
      if (source.extent(i) != destination.extent(i))
        return false;
    return true;
  }

  template <class OtherView>
  static bool same_mapping_without_pointer(const HostView& source, const OtherView& destination) {
    // Cross-space raw copies require the same layout, including any padding.
    if (!same_extents(source, destination) || source.span() != destination.span())
      return false;
    for (int i = 0; i < HostView::rank; ++i)
      if (source.stride(i) != destination.stride(i))
        return false;
    return true;
  }

  static bool same_mapping(const HostView& source, const HostView& destination) {
    // Pointer and strides distinguish batches and same-size differently strided windows.
    return source.data() == destination.data() && same_mapping_without_pointer(source, destination);
  }

  template <class View> static View allocate(const HostView& source) {
    // Preserve Kokkos initialization semantics when creating a new case field.
    if constexpr (HostView::rank == 1)
      return View("driver case/" + source.label(), source.extent(0));
    else if constexpr (HostView::rank == 2)
      return View("driver case/" + source.label(), source.extent(0), source.extent(1));
    else if constexpr (HostView::rank == 3)
      return View("driver case/" + source.label(), source.extent(0), source.extent(1),
                  source.extent(2));
    else
      return View("driver case/" + source.label(), source.extent(0), source.extent(1),
                  source.extent(2), source.extent(3));
  }

  void require_bound() const {
    // Fail at the ownership boundary rather than returning an unbound handle.
    if (!bound_)
      throw std::logic_error("HostDeviceField is not bound to a host View");
  }

  HostView host_;
  execution_view_type execution_;
  packed_view_type packed_;
  Current current_{Current::both};
  bool bound_{};
  bool pack_{};
};

} // namespace clubb_kokkos
