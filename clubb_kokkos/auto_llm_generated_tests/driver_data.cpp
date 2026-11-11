// Check persistent field transfer, aliasing, strided windows and rebinding.
//
// Run after linking Kokkos; the same checks use the default CPU or GPU execution
// space. Values are integers representable in double so the test isolates data
// movement from floating-point arithmetic differences.
#include "driver_data.hpp"

#include <atomic>
#include <iostream>
#include <stdexcept>

namespace {
std::atomic<std::uint64_t> allocations{};
std::atomic<std::uint64_t> copies{};

void count_allocation(Kokkos::Tools::SpaceHandle, const char*, const void*, std::uint64_t) {
  // Observe allocations through Kokkos Tools without changing field ownership.
  ++allocations;
}

void count_copy(Kokkos::Tools::SpaceHandle, const char*, const void*, Kokkos::Tools::SpaceHandle,
                const char*, const void*, std::uint64_t) {
  // Count transfer calls so a redundant synchronized update is detectable.
  ++copies;
}

void require(bool condition, const char* message) {
  // Report the failed ownership or transfer invariant to the invoking test.
  if (!condition)
    throw std::runtime_error(message);
}

template <class HostView> void check_rank(const HostView& host) {
  // Check all supported ranks, including an empty field, with both transfer directions.
  clubb_kokkos::HostDeviceField<HostView> field;
  Kokkos::deep_copy(host, 3.0);
  field.bind(host);
  field.sync_device();
  const auto device = field.execution_view();
  const auto saved  = device.data();
  auto check        = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, device);
  for (std::size_t i = 0; i < host.extent(0); ++i) {
    if constexpr (HostView::rank == 1) {
      require(check(i) == 3.0, "rank1 upload");
    } else {
      for (std::size_t j = 0; j < host.extent(1); ++j) {
        if constexpr (HostView::rank == 2) {
          require(check(i, j) == 3.0, "rank2 upload");
        } else {
          for (std::size_t k = 0; k < host.extent(2); ++k) {
            if constexpr (HostView::rank == 3) {
              require(check(i, j, k) == 3.0, "rank3 upload");
            } else {
              for (std::size_t l = 0; l < host.extent(3); ++l)
                require(check(i, j, k, l) == 3.0, "rank4 upload");
            }
          }
        }
      }
    }
  }
  Kokkos::deep_copy(device, 7.0);
  field.modify_device();
  field.bind(host);
  field.sync_host();
  for (std::size_t i = 0; i < host.extent(0); ++i) {
    if constexpr (HostView::rank == 1) {
      require(host(i) == 7.0, "rank1 download");
    } else {
      for (std::size_t j = 0; j < host.extent(1); ++j) {
        if constexpr (HostView::rank == 2) {
          require(host(i, j) == 7.0, "rank2 download");
        } else {
          for (std::size_t k = 0; k < host.extent(2); ++k) {
            if constexpr (HostView::rank == 3) {
              require(host(i, j, k) == 7.0, "rank3 download");
            } else {
              for (std::size_t l = 0; l < host.extent(3); ++l)
                require(host(i, j, k, l) == 7.0, "rank4 download");
            }
          }
        }
      }
    }
  }
  const auto before_allocations = allocations.load();
  const auto before_copies      = copies.load();
  field.bind(host);
  field.sync_device();
  field.sync_device();
  field.sync_host();
  require(allocations == before_allocations, "synchronized field allocated storage");
  require(copies == before_copies, "synchronized field transferred data");
  require(field.execution_view().data() == saved, "sync or identical bind reallocated");
  if constexpr (std::is_same_v<typename decltype(device)::memory_space, Kokkos::HostSpace>)
    require(device.data() == host.data(), "CPU dense field must alias host");
  field.release();
  bool rejected = false;
  try {
    field.execution_view();
  } catch (const std::logic_error&) {
    rejected = true;
  }
  require(rejected, "released field accepted access");
}

void check_strides() {
  // Update only the selected columns and preserve adjacent source storage.
  using HostArray = Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  using HostSlice = Kokkos::View<double**, Kokkos::LayoutStride, Kokkos::HostSpace>;
  HostArray storage("strided", 8, 4);
  for (int j = 0; j < 4; ++j)
    for (int i = 0; i < 8; ++i)
      storage(i, j) = i + 10 * j;
  HostSlice first = Kokkos::subview(storage, std::pair(1, 4), Kokkos::ALL);
  clubb_kokkos::HostDeviceField<HostSlice> field;
  field.bind(first);
  field.sync_device();
  auto device = field.execution_view();
  auto saved  = device.data();
  Kokkos::parallel_for(
      "field_update", Kokkos::RangePolicy<>(0, 12),
      KOKKOS_LAMBDA(int n) { device(n % 3, n / 3) = 2 * device(n % 3, n / 3) + 1; });
  field.modify_device();
  field.sync_host();
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 8; ++i) {
      const double expected = i >= 1 && i < 4 ? 2 * (i + 10 * j) + 1 : i + 10 * j;
      require(storage(i, j) == expected, "strided download or neighbor corruption");
    }
  }
  first(0, 0) = 91;
  field.modify_host();
  field.sync_device();
  auto check = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, device);
  require(check(0, 0) == 91, "strided host modification");

  HostSlice second        = Kokkos::subview(storage, std::pair(4, 7), Kokkos::ALL);
  auto before_allocations = allocations.load();
  field.bind(second);
  field.sync_device();
  require(allocations == before_allocations, "same-size batch allocated storage");
  require(field.execution_view().data() == saved, "same-size batch reallocated");
  check = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, device);
  require(check(0, 0) == 4 && check(2, 3) == 36, "new batch did not upload");

  HostSlice changed_stride(storage.data(), Kokkos::LayoutStride(3, 1, 4, 3));
  before_allocations = allocations.load();
  field.bind(changed_stride);
  field.sync_device();
  require(allocations == before_allocations, "stride-only rebind allocated storage");
  require(field.execution_view().data() == saved, "stride-only rebind reallocated");
  check = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, device);
  for (int j = 0; j < 4; ++j)
    for (int i = 0; i < 3; ++i)
      require(check(i, j) == changed_stride(i, j), "stride-only rebind did not upload");

  HostSlice resized = Kokkos::subview(storage, std::pair(0, 2), Kokkos::ALL);
  field.bind(resized);
  field.sync_device();
  require(field.execution_view().extent(0) == 2 && field.execution_view().extent(1) == 4,
          "shape change not reflected");
  field.release();
}
} // namespace

int main(int argc, char** argv) {
  // Run ownership edges in one Kokkos lifetime and report exceptions as failures.
  Kokkos::initialize(argc, argv);
  Kokkos::Tools::Experimental::set_allocate_data_callback(count_allocation);
  Kokkos::Tools::Experimental::set_begin_deep_copy_callback(count_copy);
  int status = 0;
  try {
    check_rank(Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::HostSpace>("rank1", 5));
    check_rank(Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::HostSpace>("rank2", 3, 4));
    check_rank(Kokkos::View<double***, Kokkos::LayoutLeft, Kokkos::HostSpace>("rank3", 2, 3, 4));
    check_rank(
        Kokkos::View<double****, Kokkos::LayoutLeft, Kokkos::HostSpace>("rank4", 2, 3, 4, 2));
    check_rank(Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::HostSpace>("empty", 0, 3));
    check_strides();
    std::cout << "HostDeviceField ownership checks pass\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    status = 1;
  }
  Kokkos::Tools::Experimental::set_allocate_data_callback(nullptr);
  Kokkos::Tools::Experimental::set_begin_deep_copy_callback(nullptr);
  Kokkos::finalize();
  return status;
}
