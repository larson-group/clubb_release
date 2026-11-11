// Compare every interpolation routine with the compiled current Fortran source.
// Scalar stencils cover endpoints/interior, signed/flat slopes and extrapolation.
// Profiles cover exact/tolerance matches, clamping, pressure negation, one-point
// sources and native zero-tail behavior. Both host and GPU execution are exact.
#include "CLUBB_core/interpolation.hpp"

#include <array>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

using clubb_kokkos::real_t;
extern "C" {
void f2k_interpolation_scalar_c(int, const real_t*, const int*, real_t*);
void f2k_interpolation_profile_c(int, int, const real_t*, const real_t*, const real_t*,
                                 const real_t*, const real_t*, real_t*, int*, real_t*);
void f2k_interpolation_equal_error_c();
}

template <class Result, class Reference>
void compare(const char* name, Result result, Reference reference) {
  // Check shape and every bit after transferring a real host/device routine result.
  // Pack a strided result before crossing from device to host memory.
  Kokkos::View<real_t*, Kokkos::LayoutLeft, typename Result::memory_space> packed("comparison",
                                                                                  result.extent(0));
  Kokkos::deep_copy(packed, result);
  auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, packed);
  if (host.extent(0) != reference.extent(0))
    throw std::runtime_error("interpolation shape mismatch");
  for (std::size_t i = 0; i < host.extent(0); ++i)
    if (std::memcmp(&host(i), &reference(i), sizeof(real_t))) {
      std::cerr << name << " point=" << i << std::hexfloat << " C++=" << host(i)
                << " native=" << reference(i) << '\n';
      throw std::runtime_error("interpolation differs from native bits");
    }
}

int main(int argc, char* argv[]) {
  // Native ERROR STOP is checked in a separate invocation before CUDA starts.
  if (argc > 1 && std::string(argv[1]) == "-native_equal") {
    f2k_interpolation_equal_error_c();
    return 0;
  }
  try {
    Kokkos::ScopeGuard kokkos(argc, argv);
    using namespace clubb_kokkos;
    using Host1D    = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
    using Host2D    = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
    constexpr int n = 384;
    Host2D data("data", 9, n), reference("reference", 3, n), host_result("host_result", 3, n);
    Kokkos::View<int**, Kokkos::LayoutLeft, Kokkos::HostSpace> levels("levels", 4, n);
    for (int j = 0; j < n; ++j) {
      const int branch   = j % 4;
      levels(0, j)       = branch == 0 ? 0 : branch == 3 ? 2 : -1;
      levels(1, j)       = 0;
      levels(2, j)       = 1;
      levels(3, j)       = branch == 1 ? 1 : 2;
      const real_t shift = real_t(j / 4) * real_t(0.31);
      data(1, j)         = shift;
      data(2, j)         = shift + real_t(0.7);
      data(3, j)         = shift + real_t(2.3);
      data(4, j)         = shift + real_t(4.1);
      if (branch == 0)
        data(1, j) = data(2, j);
      if (branch == 1)
        data(4, j) = data(3, j);
      data(0, j) = data(2, j) + (data(3, j) - data(2, j)) * (real_t((j / 4) % 11) / real_t(7));
      for (int q = 0; q < 4; ++q)
        data(5 + q, j) = real_t(0.137) * real_t(((j / 4 + q * q * 7) % 19) - 9);
      if (j % 16 < 4)
        for (int q = 0; q < 4; ++q)
          data(5 + q, j) = (q % 2) ? real_t(-0.0) : real_t(0.0);
    }
    f2k_interpolation_scalar_c(n, data.data(), levels.data(), reference.data());
    auto device_data = to_execution_space(data);
    Kokkos::View<int**, Kokkos::LayoutLeft> device_levels("device_levels", 4, n);
    Kokkos::deep_copy(device_levels, levels);
    Real2D device_result("device_result", 3, n);
    Kokkos::parallel_for(
        "interpolation::scalars", Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(const int j) {
          device_result(0, j) =
              lin_interpolate_two_points(device_data(0, j), device_data(3, j), device_data(2, j),
                                         device_data(7, j), device_data(6, j));
          device_result(1, j) = linear_interp_factor((device_data(0, j) - device_data(2, j)) /
                                                         (device_data(3, j) - device_data(2, j)),
                                                     device_data(7, j), device_data(6, j));
          device_result(2, j) = mono_cubic_interp(
              device_data(0, j), device_levels(0, j), device_levels(1, j), device_levels(2, j),
              device_levels(3, j), device_data(1, j), device_data(2, j), device_data(3, j),
              device_data(4, j), device_data(5, j), device_data(6, j), device_data(7, j),
              device_data(8, j));
        });
    for (int j = 0; j < n; ++j) {
      host_result(0, j) =
          lin_interpolate_two_points(data(0, j), data(3, j), data(2, j), data(7, j), data(6, j));
      host_result(1, j) = linear_interp_factor(
          (data(0, j) - data(2, j)) / (data(3, j) - data(2, j)), data(7, j), data(6, j));
      host_result(2, j) = mono_cubic_interp(
          data(0, j), levels(0, j), levels(1, j), levels(2, j), levels(3, j), data(1, j),
          data(2, j), data(3, j), data(4, j), data(5, j), data(6, j), data(7, j), data(8, j));
    }
    for (int q = 0; q < 3; ++q) {
      compare("GPU scalar", Kokkos::subview(device_result, q, Kokkos::ALL),
              Kokkos::subview(reference, q, Kokkos::ALL));
      compare("host scalar", Kokkos::subview(host_result, q, Kokkos::ALL),
              Kokkos::subview(reference, q, Kokkos::ALL));
    }

    for (int nsrc : {1, 2, 3, 6}) {
      constexpr int nout = 12;
      Host1D grid("grid", nsrc), values("values", nsrc), target("target", nout);
      Host1D pressure("pressure", nsrc), pressure_out("pressure_out", nout),
          on_grid("on_grid", nsrc);
      Host2D outputs("outputs", nout, 3);
      Kokkos::View<int*, Kokkos::HostSpace> search("search", nout);
      for (int k = 0; k < nsrc; ++k) {
        grid(k)     = nsrc == 3 ? real_t(k) * real_t(0.4e-6) : real_t(k * k) * real_t(0.17);
        values(k)   = real_t(0.731) * real_t(k * k - 3 * k + 1);
        pressure(k) = -grid(k);
      }
      std::array<real_t, nout> targets{real_t(-1),        real_t(0),    real_t(0.03), real_t(0.17),
                                       real_t(0.1700003), real_t(0.68), real_t(0.91), real_t(1.53),
                                       real_t(2.72),      real_t(4.25), real_t(9),    real_t(0.1)};
      if (nsrc == 3)
        targets = {real_t(-1e-6),  real_t(0),      real_t(0.2e-6), real_t(0.4e-6),
                   real_t(0.7e-6), real_t(0.8e-6), real_t(1.1e-6), real_t(1.6e-6),
                   real_t(2e-6),   real_t(3e-6),   real_t(9),      real_t(0.1e-6)};
      for (int k = 0; k < nout; ++k) {
        target(k)       = targets[k];
        pressure_out(k) = -target(k);
      }
      f2k_interpolation_profile_c(nout, nsrc, target.data(), grid.data(), values.data(),
                                  pressure_out.data(), pressure.data(), outputs.data(),
                                  search.data(), on_grid.data());
      auto grid_d = to_execution_space(grid), values_d = to_execution_space(values),
           target_d       = to_execution_space(target);
      auto pressure_d     = to_execution_space(pressure),
           pressure_out_d = to_execution_space(pressure_out);
      compare("GPU height", zlinterp_fnc(nout, nsrc, target_d, grid_d, values_d),
              Kokkos::subview(outputs, Kokkos::ALL, 0));
      compare("GPU pressure", plinterp_fnc(nout, nsrc, pressure_out_d, pressure_d, values_d),
              Kokkos::subview(outputs, Kokkos::ALL, 1));
      compare("GPU clamped", lin_interp_between_grids(nout, nsrc, target_d, grid_d, values_d),
              Kokkos::subview(outputs, Kokkos::ALL, 2));
      compare("host height", zlinterp_fnc(nout, nsrc, target, grid, values),
              Kokkos::subview(outputs, Kokkos::ALL, 0));
      compare("host pressure", plinterp_fnc(nout, nsrc, pressure_out, pressure, values),
              Kokkos::subview(outputs, Kokkos::ALL, 1));
      compare("host clamped", lin_interp_between_grids(nout, nsrc, target, grid, values),
              Kokkos::subview(outputs, Kokkos::ALL, 2));

      // The driver passes column profiles from a column-first two-dimensional
      // grid. Exercise those strides in both memory spaces, including guards.
      Host2D source_columns("source_columns", 4, nsrc), value_columns("value_columns", 4, nsrc),
          target_columns("target_columns", 4, nout);
      Kokkos::deep_copy(source_columns, real_t(-999));
      Kokkos::deep_copy(value_columns, real_t(-999));
      Kokkos::deep_copy(target_columns, real_t(-999));
      Kokkos::deep_copy(Kokkos::subview(source_columns, 1, Kokkos::ALL), grid);
      Kokkos::deep_copy(Kokkos::subview(value_columns, 1, Kokkos::ALL), values);
      Kokkos::deep_copy(Kokkos::subview(target_columns, 1, Kokkos::ALL), target);
      auto source_columns_d = to_execution_space(source_columns),
           value_columns_d  = to_execution_space(value_columns),
           target_columns_d = to_execution_space(target_columns);
      auto source_stride    = Kokkos::subview(source_columns, 1, Kokkos::ALL),
           value_stride     = Kokkos::subview(value_columns, 1, Kokkos::ALL),
           target_stride    = Kokkos::subview(target_columns, 1, Kokkos::ALL);
      auto source_stride_d  = Kokkos::subview(source_columns_d, 1, Kokkos::ALL),
           value_stride_d   = Kokkos::subview(value_columns_d, 1, Kokkos::ALL),
           target_stride_d  = Kokkos::subview(target_columns_d, 1, Kokkos::ALL);
      compare("host strided",
              lin_interp_between_grids(nout, nsrc, target_stride, source_stride, value_stride),
              Kokkos::subview(outputs, Kokkos::ALL, 2));
      compare(
          "GPU strided",
          lin_interp_between_grids(nout, nsrc, target_stride_d, source_stride_d, value_stride_d),
          Kokkos::subview(outputs, Kokkos::ALL, 2));
      Kokkos::View<int*> search_d("search_d", nout);
      Kokkos::parallel_for(
          "interpolation::search", Kokkos::RangePolicy<>(0, nout),
          KOKKOS_LAMBDA(const int k) { search_d(k) = binary_search(nsrc, grid_d, target_d(k)); });
      auto search_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, search_d);
      for (int k = 0; k < nout; ++k) {
        const int expected = search(k) < 0 ? -1 : search(k) - 1;
        if (search_h(k) != expected || binary_search(nsrc, grid, target(k)) != expected)
          throw std::runtime_error("upper-bracket index mismatch");
      }
      if (nsrc > 1)
        for (int k = 0; k < nsrc; ++k) {
          real_t got = real_t(-999);
          lin_interpolate_on_grid_api(nsrc, source_stride, value_stride, grid(k), got);
          if (std::memcmp(&got, &on_grid(k), sizeof(real_t)))
            throw std::runtime_error("on-grid scalar mismatch");
          lin_interpolate_on_grid_api(nsrc, source_stride_d, value_stride_d, grid(k), got);
          if (std::memcmp(&got, &on_grid(k), sizeof(real_t)))
            throw std::runtime_error("GPU on-grid scalar mismatch");
        }
      const auto source_return = to_host_space(source_columns_d),
                 values_return = to_host_space(value_columns_d),
                 target_return = to_host_space(target_columns_d);
      if (std::memcmp(source_return.data(), source_columns.data(),
                      source_columns.size() * sizeof(real_t)) ||
          std::memcmp(values_return.data(), value_columns.data(),
                      value_columns.size() * sizeof(real_t)) ||
          std::memcmp(target_return.data(), target_columns.data(),
                      target_columns.size() * sizeof(real_t)))
        throw std::runtime_error("interpolation changed input/guard columns");
    }
    bool caught = false;
    try {
      (void)lin_interpolate_two_points(real_t(0), real_t(1), real_t(1), real_t(2), real_t(3));
    } catch (const std::runtime_error&) {
      caught = true;
    }
    if (!caught)
      throw std::runtime_error("coincident host heights were not rejected");
    Host1D x("x", 2), t("t", 2);
    x(0)   = real_t(1);
    x(1)   = real_t(0);
    caught = false;
    try {
      real_t value;
      lin_interpolate_on_grid_api(2, x, t, real_t(0.5), value);
    } catch (const std::runtime_error&) {
      caught = true;
    }
    if (!caught)
      throw std::runtime_error("unsorted scalar grid was not rejected");
    caught = false;
    try {
      (void)lin_interp_between_grids(2, 2, t, x, t);
    } catch (const std::runtime_error&) {
      caught = true;
    }
    if (!caught)
      throw std::runtime_error("unsorted profile grid was not rejected");
    x(0)   = real_t(0);
    x(1)   = real_t(1);
    caught = false;
    try {
      real_t value;
      lin_interpolate_on_grid_api(2, x, t, real_t(2), value);
    } catch (const std::runtime_error&) {
      caught = true;
    }
    if (!caught)
      throw std::runtime_error("out-of-range scalar grid was not rejected");
    std::cout << "All interpolation routines match native host/GPU bits and diagnostic guards\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
