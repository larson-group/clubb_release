// Check nested scratch reservations, initialization, growth and return lifetimes.
//
// Run with Serial or CUDA Kokkos. Integer-valued doubles isolate ownership,
// initialization and ordered reuse from floating-point differences. CTest checks
// warm allocation counts, nested/returned lifetimes, invalid extents and cleanup.
#include "scratch_workspace.hpp"

#include <atomic>
#include <cstring>
#include <iostream>
#include <stdexcept>

namespace {
using Real1D     = Kokkos::View<double*, Kokkos::LayoutLeft>;
using Real2D     = Kokkos::View<double**, Kokkos::LayoutLeft>;
using Real3D     = Kokkos::View<double***, Kokkos::LayoutLeft>;
using Real4D     = Kokkos::View<double****, Kokkos::LayoutLeft>;
using HostReal2D = Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
std::atomic<std::uint64_t> allocations{}, releases{}, copies{};

void count_allocation(Kokkos::Tools::SpaceHandle, const char* label, const void*, std::uint64_t) {
  // Ignore observation mirrors; count only chunks owned by this workspace.
  if (label && std::strncmp(label, "scratch/", 8) == 0)
    ++allocations;
}
void count_release(Kokkos::Tools::SpaceHandle, const char* label, const void*, std::uint64_t) {
  // Check that clear returns retained chunks, excluding observation mirrors.
  if (label && std::strncmp(label, "scratch/", 8) == 0)
    ++releases;
}
void count_copy(Kokkos::Tools::SpaceHandle, const char*, const void*, Kokkos::Tools::SpaceHandle,
                const char*, const void*, std::uint64_t) {
  // Detect initialization/copy work on otherwise empty or released scopes.
  ++copies;
}
void require(bool condition, const char* message) {
  // Identify a broken reservation invariant in CTest output.
  if (!condition)
    throw std::runtime_error(message);
}

template <class View> void require_filled(const View& view, double value) {
  // Copy a contiguous arithmetic View only at a test observation boundary.
  auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, view);
  for (std::size_t i = 0; i < view.size(); ++i)
    require(host.data()[i] == value, "unexpected scratch field value");
}

Real2D returned_result(clubb_kokkos::ScratchWorkspace& scratch, int count) {
  // Reserve the result in the caller's scope before opening helper-local storage.
  auto result = scratch.view<Real2D>("returned_result::result", count, 2);
  auto local  = scratch.scope();
  auto temporary =
      local.view<Real2D>(Kokkos::WithoutInitializing, "returned_result::temporary", count, 2);
  Kokkos::parallel_for(
      "returned_result::fill", count * 2, KOKKOS_LAMBDA(int i) {
        temporary(i % count, i / count) = 9;
        result(i % count, i / count)    = temporary(i % count, i / count) + 4;
      });
  return result;
}

void check_nested_and_returned(clubb_kokkos::ScratchWorkspace& scratch) {
  // Grow beneath a live parent field and retain a helper result after its local scope.
  auto parent         = scratch.scope();
  auto field          = parent.view<Real2D>("parent", 5, 3);
  const auto* pointer = field.data();
  Kokkos::deep_copy(field, 7.0);
  {
    auto child = scratch.scope();
    auto large = child.view<Real3D>("child large", 16, 17, 18);
    require_filled(large, 0.0);
    Kokkos::deep_copy(large, 23.0);
    require(field.data() == pointer, "growth relocated a live parent field");
    bool rejected = false;
    try {
      parent.view<Real1D>("parent while child active", 1);
    } catch (const std::logic_error&) {
      rejected = true;
    }
    require(rejected, "parent allocated while child scope was active");
  }
  require_filled(field, 7.0);
  auto result = returned_result(scratch, 5);
  {
    auto sibling   = scratch.scope();
    auto temporary = sibling.view<Real4D>("sibling", 2, 3, 4, 5);
    require_filled(temporary, 0.0);
    Kokkos::deep_copy(temporary, 31.0);
  }
  require_filled(field, 7.0);
  require_filled(result, 13.0);
}

void check_reuse_and_initialization(clubb_kokkos::ScratchWorkspace& scratch) {
  // Queue writes, release metadata only, then initialize the same reservation.
  // No test fence is inserted between the two default-queue scopes.
  double* saved = nullptr;
  {
    auto local = scratch.scope();
    auto view  = local.view<Real1D>(Kokkos::WithoutInitializing, "reuse", 33);
    saved      = view.data();
    require(reinterpret_cast<std::uintptr_t>(saved) % 128 == 0, "scratch is misaligned");
    Kokkos::parallel_for("reuse::write", 33, KOKKOS_LAMBDA(int i) { view(i) = 19.; });
  }
  {
    auto local = scratch.scope();
    auto view  = local.view<Real1D>(Kokkos::WithoutInitializing, "reuse retained", 33);
    require(view.data() == saved, "same-shaped call did not reuse reservation");
    require_filled(view, 19.0);
    Kokkos::parallel_for("reuse::queued_write", 33, KOKKOS_LAMBDA(int i) { view(i) = 29.; });
  }
  {
    auto local = scratch.scope();
    auto view  = local.view<Real1D>("reuse zeroed", 33);
    require(view.data() == saved, "initialization changed reservation");
    Kokkos::parallel_for("reuse::consume_zero", 33, KOKKOS_LAMBDA(int i) { view(i) += 1.; });
    require_filled(view, 1.0);
  }
  {
    auto local = scratch.scope();
    auto host  = local.view<HostReal2D>("host", 3, 7);
    require_filled(host, 0.0);
    host(2, 6) = 41.;
  }
  {
    auto local = scratch.scope();
    auto host  = local.view<HostReal2D>(Kokkos::WithoutInitializing, "host retained", 3, 7);
    require(host(2, 6) == 41., "host scratch was not retained without initialization");
  }
}

void check_edges(clubb_kokkos::ScratchWorkspace& scratch) {
  // Empty fields reserve/initialize nothing; programmer errors fail at the boundary.
  bool rejected = false;
  try {
    scratch.view<Real1D>("no caller scope", 1);
  } catch (const std::logic_error&) {
    rejected = true;
  }
  require(rejected, "output allocation accepted without caller scope");
  auto local                    = scratch.scope();
  const auto before_allocations = allocations.load(), before_copies = copies.load();
  auto empty = local.view<Real3D>("empty", 3, 0, 2);
  require(empty.extent(0) == 3 && empty.extent(1) == 0 && empty.extent(2) == 2,
          "empty field shape changed");
  require(empty.data() == nullptr && empty.size() == 0, "empty field reserved bytes");
  require(allocations == before_allocations && copies == before_copies,
          "empty field allocated or initialized");
  rejected = false;
  try {
    local.view<Real2D>("negative", -1, 2);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  require(rejected, "negative extent was accepted");
  rejected = false;
  try {
    local.view<Real2D>("overflow", std::numeric_limits<std::size_t>::max(), 2);
  } catch (const std::length_error&) {
    rejected = true;
  }
  require(rejected, "overflowing extent product was accepted");
  rejected = false;
  try {
    scratch.clear();
  } catch (const std::logic_error&) {
    rejected = true;
  }
  require(rejected, "clear accepted an active scope");
  using HostInt = Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  auto integers = local.view<HostInt>("integers", 17);
  require_filled(integers, 0.0);
}

} // namespace

int main(int argc, char** argv) {
  // Check one case's warm reuse, then release chunks before Kokkos finalization.
  Kokkos::ScopeGuard kokkos(argc, argv);
  Kokkos::Tools::Experimental::set_allocate_data_callback(count_allocation);
  Kokkos::Tools::Experimental::set_deallocate_data_callback(count_release);
  Kokkos::Tools::Experimental::set_begin_deep_copy_callback(count_copy);
  int status = 0;
  try {
    clubb_kokkos::ScratchWorkspace scratch(256, 256);
    check_edges(scratch);
    check_nested_and_returned(scratch);
    check_reuse_and_initialization(scratch);
    const auto after_warmup = allocations.load();
    for (int repeat = 0; repeat < 3; ++repeat) {
      check_edges(scratch);
      check_nested_and_returned(scratch);
      check_reuse_and_initialization(scratch);
    }
    require(allocations == after_warmup, "warm routine path allocated new Kokkos chunks");
    Kokkos::fence("scratch test final completion");
    const auto before_clear = releases.load();
    scratch.clear();
    require(releases > before_clear, "clear did not release retained storage");
    require(scratch.reserved_device_bytes() == 0 && scratch.reserved_host_bytes() == 0,
            "clear retained capacity");
    std::cout << "ScratchWorkspace ownership checks pass\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    status = 1;
  }
  Kokkos::Tools::Experimental::set_allocate_data_callback(nullptr);
  Kokkos::Tools::Experimental::set_deallocate_data_callback(nullptr);
  Kokkos::Tools::Experimental::set_begin_deep_copy_callback(nullptr);
  return status;
}
