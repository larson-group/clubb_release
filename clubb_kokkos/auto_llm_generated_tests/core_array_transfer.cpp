// Check native constant bits and the host/device boundary used by the core API.
// Strided column windows update only their selected columns. Rank-three/four
// fields and zero-scalar fields retain their values and logical dimensions.
#include "scratch_workspace_owner.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/fill_holes.hpp"
#include "CLUBB_core/model_flags.hpp"
#include "CLUBB_core/parameter_indices.hpp"
#include "f2k/array_bridge.hpp"

#include <array>
#include <cstring>
#include <iostream>
#include <stdexcept>

extern "C" void f2k_core_constants_c(clubb_kokkos::real_t*, int*, int*);

int main(int argc, char* argv[]) {
  // Use the compiled native constants as the oracle, then mutate device copies
  // of full/strided host fields and verify both returned values and guard columns.
  try {
    Kokkos::ScopeGuard kokkos(argc, argv);
    clubb_kokkos::test::ScratchWorkspaceOwner scratch_owner;
    auto& workspace = scratch_owner.workspace;
    using namespace clubb_kokkos;
    std::array<real_t, 6> native_values;
    std::array<int, 6> native_indices;
    std::array<int, 10> native_flags;
    f2k_core_constants_c(native_values.data(), native_indices.data(), native_flags.data());
    const std::array<real_t, 6> values{zero, two, w_tol, w_tol_sqd, em_min, zero_threshold};
    if (std::memcmp(values.data(), native_values.data(), sizeof(values)))
      throw std::runtime_error("core constant bits differ from Fortran");
    const std::array<int, 6> indices{iC_wp2_splat + 1,  ic_K + 1,     ilambda0_stability_coef + 1,
                                     iup2_sfc_coef + 1, ia_const + 1, ibv_efold + 1};
    if (indices != native_indices)
      throw std::runtime_error("native parameter index mapping");
    const std::array<int, 10> flags{order_xm_wpxp,
                                    order_xp2_xpyp,
                                    order_wp2_wp3,
                                    order_windm,
                                    iiPDF_ADG1,
                                    ipdf_pre_advance_fields,
                                    ipdf_post_advance_fields,
                                    ipdf_pre_post_advance_fields,
                                    l_advance_xp3,
                                    Lscale_tau};
    if (flags != native_flags)
      throw std::runtime_error("native model constants");

    Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace> host("host", 5, 4);
    for (int k = 0; k < 4; ++k)
      for (int i = 0; i < 5; ++i)
        host(i, k) = real_t(i + 10 * k);
    auto columns = Kokkos::subview(host, std::make_pair(1, 4), Kokkos::ALL);
    auto device  = to_execution_space(columns);
    Kokkos::parallel_for(
        "core_array_transfer::strided_columns", Kokkos::RangePolicy<>(0, 12),
        KOKKOS_LAMBDA(const int cell) {
          device(cell % 3, cell / 3) = real_t(2) * device(cell % 3, cell / 3) + real_t(1);
        });
    from_execution_space(columns, device);
    for (int k = 0; k < 4; ++k)
      for (int i = 0; i < 5; ++i) {
        const real_t before   = real_t(i + 10 * k);
        const real_t expected = (i == 0 || i == 4) ? before : real_t(2) * before + real_t(1);
        if (host(i, k) != expected)
          throw std::runtime_error("strided transfer overwrote values/guards");
      }

    Real2D existing("existing", 3, 4);
    auto alias = to_execution_space(existing);
    if (alias.data() != existing.data())
      throw std::runtime_error("device input was copied");

    // Dense LayoutStride windows share storage and retain its allocation lifetime.
    auto dense_window       = Kokkos::subview(existing, Kokkos::ALL, std::make_pair(1, 3));
    auto dense_alias        = to_execution_space(dense_window);
    auto scoped_dense_alias = to_execution_space(workspace, dense_window);
    auto out_dense_alias = to_execution_space(workspace, Kokkos::WithoutInitializing, dense_window);
    if (dense_alias.data() != dense_window.data() ||
        scoped_dense_alias.data() != dense_window.data() ||
        out_dense_alias.data() != dense_window.data())
      throw std::runtime_error("dense execution subview was copied");
    auto retained_alias = [] {
      Real2D owner("retained_dense_owner", 3, 4);
      Kokkos::deep_copy(owner, real_t(17));
      Kokkos::View<real_t**, Kokkos::LayoutStride> strided(owner);
      return to_execution_space(strided);
    }();
    if (!retained_alias.use_count() || to_host_space(retained_alias)(2, 3) != real_t(17))
      throw std::runtime_error("dense execution alias lost allocation ownership");

    Kokkos::deep_copy(existing, real_t(5));
    auto snapshot = to_host_space(existing);
    Kokkos::deep_copy(existing, real_t(9));
    if (snapshot(1, 2) != real_t(5) || snapshot.data() == existing.data())
      throw std::runtime_error("device snapshot aliases the original field");
    // Mirrors retain CPU aliases; device mirrors reserve separate host scratch.
    auto host_alias_mirror = create_host_mirror(workspace, host);
    if (host_alias_mirror.data() != host.data())
      throw std::runtime_error("workspace host mirror lost the CPU alias");
    auto execution_mirror = create_host_mirror(workspace, existing);
    if (execution_mirror.extent(0) != existing.extent(0) ||
        execution_mirror.extent(1) != existing.extent(1))
      throw std::runtime_error("workspace mirror changed logical shape");
    if constexpr (std::is_same_v<typename Real2D::memory_space,
                                 typename decltype(execution_mirror)::memory_space>) {
      if (execution_mirror.data() != existing.data())
        throw std::runtime_error("execution-space host mirror lost its alias");
    } else if (execution_mirror.data() == existing.data())
      throw std::runtime_error("device mirror aliases inaccessible device storage");
    Kokkos::View<const real_t**, Kokkos::LayoutLeft> const_execution(existing);
    auto const_mirror = create_host_mirror(workspace, const_execution);
    if (const_mirror.data() == existing.data())
      throw std::runtime_error("mutable mirror aliases a const input");
    auto host_snapshot  = to_host_space(host);
    host_snapshot(1, 2) = real_t(-99);
    if (host(1, 2) == real_t(-99) || host_snapshot.data() == host.data())
      throw std::runtime_error("host snapshot aliases the original field");

    Kokkos::View<real_t***, Kokkos::LayoutLeft, Kokkos::HostSpace> scalars("scalars", 3, 4, 2);
    for (int s = 0; s < 2; ++s)
      for (int k = 0; k < 4; ++k)
        for (int i = 0; i < 3; ++i)
          scalars(i, k, s) = real_t(i + 10 * k + 100 * s);
    auto scalar_device = to_execution_space(scalars);
    Kokkos::parallel_for(
        "core_array_transfer::scalars", Kokkos::RangePolicy<>(0, 24),
        KOKKOS_LAMBDA(const int cell) {
          scalar_device(cell % 3, (cell / 3) % 4, cell / 12) += real_t(1);
        });
    from_execution_space(scalars, scalar_device);
    for (int s = 0; s < 2; ++s)
      for (int k = 0; k < 4; ++k)
        for (int i = 0; i < 3; ++i)
          if (scalars(i, k, s) != real_t(i + 10 * k + 100 * s + 1))
            throw std::runtime_error("rank-three transfer");

    // Hydrometeor PDF correlation matrices add two species axes. Check a
    // strided column window through execution staging and diagnostic packing.
    Kokkos::View<real_t****, Kokkos::LayoutLeft, Kokkos::HostSpace> correlations("correlations", 5,
                                                                                 4, 2, 2);
    for (int n = 0; n < 2; ++n)
      for (int m = 0; m < 2; ++m)
        for (int k = 0; k < 4; ++k)
          for (int i = 0; i < 5; ++i)
            correlations(i, k, m, n) = real_t(i + 10 * k + 100 * m + 1000 * n);
    auto correlation_columns =
        Kokkos::subview(correlations, std::make_pair(1, 4), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
    auto correlation_execution = to_execution_space(correlation_columns);
    Kokkos::parallel_for(
        "core_array_transfer::correlations", Kokkos::RangePolicy<>(0, 48),
        KOKKOS_LAMBDA(const int cell) {
          correlation_execution(cell % 3, (cell / 3) % 4, (cell / 12) % 2, cell / 24) += real_t(7);
        });
    from_execution_space(correlation_columns, correlation_execution);
    auto correlation_window     = Kokkos::subview(correlation_execution, std::make_pair(0, 2),
                                                  Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
    const auto correlation_host = to_host_space(correlation_window);
    for (int n = 0; n < 2; ++n)
      for (int m = 0; m < 2; ++m)
        for (int k = 0; k < 4; ++k) {
          for (int i = 0; i < 5; ++i) {
            const real_t before   = real_t(i + 10 * k + 100 * m + 1000 * n);
            const real_t expected = (i == 0 || i == 4) ? before : before + real_t(7);
            if (correlations(i, k, m, n) != expected)
              throw std::runtime_error("rank-four transfer changed values/guard columns");
          }
          for (int i = 0; i < 2; ++i)
            if (correlation_host(i, k, m, n) != real_t(i + 1 + 10 * k + 100 * m + 1000 * n + 7))
              throw std::runtime_error("rank-four diagnostic window has wrong indices");
        }

    // A one-column window keeps its four-column parent's leading dimension.
    // Detect and fill a high-level hole without reading or changing guard columns.
    Real2D fill_field("padded_fill_field", 4, 12);
    Real2D fill_spacing("padded_fill_spacing", 4, 12);
    Real2D fill_density("padded_fill_density", 4, 12);
    auto fill_host = Kokkos::create_mirror_view(fill_field);
    for (int k = 0; k < 12; ++k)
      for (int i = 0; i < 4; ++i)
        fill_host(i, k) = i == 2 ? zero : real_t(10 * i + k);
    fill_host(2, 8) = real_t(3);
    fill_host(2, 9) = real_t(-1);
    Kokkos::deep_copy(fill_field, fill_host);
    Kokkos::deep_copy(fill_spacing, one);
    Kokkos::deep_copy(fill_density, one);
    auto fill_column = Kokkos::subview(fill_field, std::make_pair(2, 3), Kokkos::ALL);
    auto fill_dz     = Kokkos::subview(fill_spacing, std::make_pair(2, 3), Kokkos::ALL);
    auto fill_rho    = Kokkos::subview(fill_density, std::make_pair(2, 3), Kokkos::ALL);
    fill_holes_vertical_api(workspace, 12, 1, zero, 0, 11, fill_dz, fill_rho, 1, global_fill,
                            fill_column);
    auto filled_host = to_host_space(fill_field);
    for (int k = 0; k < 12; ++k)
      for (int i = 0; i < 4; ++i) {
        const real_t expected = i == 2 ? (k == 8 ? real_t(2) : zero) : real_t(10 * i + k);
        if (filled_host(i, k) != expected)
          throw std::runtime_error("padded column hole filling changed values/guard columns");
      }

    Kokkos::View<real_t***, Kokkos::LayoutLeft, Kokkos::HostSpace> empty("empty", 3, 4, 0);
    auto empty_device = to_execution_space(empty);
    from_execution_space(empty, empty_device);
    if (empty_device.extent(0) != 3 || empty_device.extent(1) != 4 || empty_device.extent(2) != 0)
      throw std::runtime_error("zero-scalar shape was changed");

    // SILHS receives the first grid column as a strided device profile. Pack
    // that profile and return native-style bytes into a different column,
    // checking both directions and the untouched surrounding columns.
    auto grid_device = to_execution_space(host);
    auto profile     = Kokkos::subview(grid_device, 1, Kokkos::ALL);
    f2k::Transfer transfer(1);
    auto field = f2k::pack_array(profile, transfer);
    for (int k = 0; k < 4; ++k) {
      real_t value{};
      auto* bytes = static_cast<unsigned char*>(field.data) + k * sizeof(real_t);
      std::memcpy(&value, bytes, sizeof(value));
      if (value != host(1, k))
        throw std::runtime_error("strided device profile packed the wrong column");
      value += real_t(1000);
      std::memcpy(bytes, &value, sizeof(value));
    }
    auto destination = Kokkos::subview(grid_device, 3, Kokkos::ALL);
    f2k::unpack_into(destination, field, "grid profile");
    const auto returned = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, grid_device);
    for (int k = 0; k < 4; ++k)
      for (int i = 0; i < 5; ++i)
        if (returned(i, k) != (i == 3 ? host(1, k) + real_t(1000) : host(i, k)))
          throw std::runtime_error("strided device output changed guard columns");
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
