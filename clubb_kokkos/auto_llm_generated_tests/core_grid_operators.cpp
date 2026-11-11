// Compare every staggered-grid operator overload with current compiled Fortran.
// Native setup supplies the same stretched grids to both implementations. Tests
// cover column-varying heights/fields, both directions, optional floors, boundary
// extensions, smoothing order, scalar/profile/columns and preserved input bytes.
#include "CLUBB_core/grid_class.hpp"
#include "f2k/derived_type_bridge.hpp"

#include <cstring>
#include <iostream>
#include <stdexcept>

using clubb_kokkos::real_t;
extern "C" {
void f2k_grid_operators_c(int, int, int, const real_t*, const real_t*, const real_t*, const real_t*,
                          real_t, real_t*, real_t*, real_t*, real_t*, f2k::WireField*);
void f2k_grid_flip_c(int, const real_t*, real_t*);
}

template <class Result, class Reference>
void compare(const char* name, Result result, Reference reference) {
  // Read actual GPU results and compare shape and every floating-point bit.
  auto host = clubb_kokkos::to_host_space(result);
  for (int d = 0; d < Result::rank; ++d)
    if (host.extent(d) != reference.extent(d))
      throw std::runtime_error("grid operator shape mismatch");
  const int columns = host.extent(0);
  const int levels  = Result::rank == 1 ? 1 : host.extent(1);
  for (int k = 0; k < levels; ++k)
    for (int i = 0; i < columns; ++i) {
      real_t got, expected;
      if constexpr (Result::rank == 1) {
        got      = host(i);
        expected = reference(i);
      } else {
        got      = host(i, k);
        expected = reference(i, k);
      }
      if (std::memcmp(&got, &expected, sizeof(real_t))) {
        std::cerr << name << " i=" << i << " k=" << k << std::hexfloat << " C++=" << got
                  << " native=" << expected << '\n';
        throw std::runtime_error("grid operator differs from native bits");
      }
    }
}

int main(int argc, char* argv[]) {
  // Exercise actual native grid construction and the independent GPU operators.
  try {
    Kokkos::ScopeGuard kokkos(argc, argv);
    using namespace clubb_kokkos;
    ScratchWorkspace workspace;
    using Host1D       = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
    using Host2D       = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
    using Host3D       = Kokkos::View<real_t***, Kokkos::LayoutLeft, Kokkos::HostSpace>;
    int configurations = 0;
    for (int nzm : {3, 4, 7, 33})
      for (int nc : {1, 5, 129})
        for (int ascending : {0, 1})
          for (real_t minimum : {real_t(-0.0), real_t(0.237)}) {
            auto workspace_scope = workspace.scope();
            const int nzt = nzm - 1;
            Host2D momentum("momentum", nc, nzm), thermo("thermo", nc, nzt), azm("azm", nc, nzm),
                azt("azt", nc, nzt);
            for (int k = 0; k < nzm; ++k)
              for (int i = 0; i < nc; ++i) {
                const int height_level = ascending ? k : nzm - 1 - k;
                momentum(i, k)         = real_t(3 + i) + real_t(0.71) * real_t(height_level) +
                                         real_t(0.137) * real_t(height_level * height_level) +
                                         real_t(0.013) * real_t(i * height_level);
                azm(i, k)              = real_t(0.317) * real_t((i * 5 + k * k * 3) % 23 - 11);
                if (i % 8 == 0 && k % 3 == 0)
                  azm(i, k) = real_t(-0.0);
              }
            for (int k = 0; k < nzt; ++k)
              for (int i = 0; i < nc; ++i) {
                thermo(i, k) = momentum(i, k) + (momentum(i, k + 1) - momentum(i, k)) *
                                                    (real_t(0.27) + real_t(i % 5) * real_t(0.08));
                azt(i, k)    = real_t(0.183) * real_t((i * 7 + k * k * 5) % 19 - 9);
                if (i % 7 == 0 && k % 3 == 0)
                  azt(i, k) = real_t(-0.0);
              }
            Host3D output_m("output_m", nc, nzm, 5), output_t("output_t", nc, nzt, 5);
            Host2D profile_m("profile_m", nzm, 5), profile_t("profile_t", nzt, 5);
            f2k::Transfer wire(16, true);
            f2k_grid_operators_c(nzm, nc, ascending, momentum.data(), thermo.data(), azm.data(),
                                 azt.data(), minimum, output_m.data(), output_t.data(),
                                 profile_m.data(), profile_t.data(), wire.data());
            grid gr;
            f2k::to_kokkos(gr, wire);
            // Check both language directions independently: symmetric roundtrip
            // equality alone would also pass if both endpoint adapters were absent.
            if (gr.k_lb_zm != (ascending ? 0 : nzm - 1) ||
                gr.k_lb_zt != (ascending ? 0 : nzt - 1) ||
                gr.k_ub_zm != (ascending ? nzm - 1 : 0) || gr.k_ub_zt != (ascending ? nzt - 1 : 0))
              throw std::runtime_error("Native grid endpoints were not converted to C++ slots");
            const auto packed_grid = f2k::to_fortran(gr);
            for (int endpoint = 10; endpoint < 14; ++endpoint) {
              int native_value, packed_value;
              std::memcpy(&native_value, wire.fields[endpoint].data, sizeof(int));
              std::memcpy(&packed_value, packed_grid.fields[endpoint].data, sizeof(int));
              if (packed_value != native_value)
                throw std::runtime_error("C++ grid endpoint was not converted to Native ordinal");
            }
            auto m = to_execution_space(azm), t = to_execution_space(azt);
            compare("zt2zm", zt2zm_api(nzm, nzt, nc, gr, t),
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 0));
            compare("zt2zm floor", zt2zm_api(nzm, nzt, nc, gr, t, minimum),
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 1));
            compare("zm2zt2zm", zm2zt2zm(nzm, nzt, nc, gr, m),
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 2));
            compare("zm2zt2zm floor", zm2zt2zm(nzm, nzt, nc, gr, m, minimum),
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 3));
            compare("ddzt", ddzt(nzm, nzt, nc, gr, t),
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 4));
            compare("zm2zt", zm2zt_api(nzm, nzt, nc, gr, m),
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 0));
            compare("zm2zt floor", zm2zt_api(nzm, nzt, nc, gr, m, minimum),
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 1));
            compare("zt2zm2zt", zt2zm2zt(nzm, nzt, nc, gr, t),
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 2));
            compare("zt2zm2zt floor", zt2zm2zt(nzm, nzt, nc, gr, t, minimum),
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 3));
            compare("ddzm", ddzm(nzm, nzt, nc, gr, m),
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 4));
            // Explicit outputs must match the same Native oracle, including
            // boundary floors, descending grids and signed zeros.
            Real2D filled_m("filled_m", nc, nzm), filled_t("filled_t", nc, nzt);
            zt2zm_api(nzm, nzt, nc, gr, t, filled_m, minimum);
            compare("zt2zm existing output", filled_m,
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 1));
            zm2zt_api(nzm, nzt, nc, gr, m, filled_t, minimum);
            compare("zm2zt existing output", filled_t,
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 1));
            ddzt(nzm, nzt, nc, gr, t, filled_m);
            compare("ddzt existing output", filled_m,
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 4));
            ddzm(nzm, nzt, nc, gr, m, filled_t);
            compare("ddzm existing output", filled_t,
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 4));
            zt2zm_api(workspace, nzm, nzt, nc, gr, t, filled_m, minimum);
            compare("scoped direct zt2zm_api", filled_m,
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 1));
            auto scoped_zt2zm_api = zt2zm_api(workspace, nzm, nzt, nc, gr, t, minimum);
            // Later reservations must not reuse storage belonging to the returned field.
            auto later_zt2zm_api = workspace_scope.view<Real2D>("later reservation", nc, nzm);
            Kokkos::deep_copy(later_zt2zm_api, real_t(57));
            compare("scoped returned zt2zm_api", scoped_zt2zm_api,
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 1));
            zm2zt_api(workspace, nzm, nzt, nc, gr, m, filled_t, minimum);
            compare("scoped direct zm2zt_api", filled_t,
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 1));
            auto scoped_zm2zt_api = zm2zt_api(workspace, nzm, nzt, nc, gr, m, minimum);
            // Later reservations must not reuse storage belonging to the returned field.
            auto later_zm2zt_api = workspace_scope.view<Real2D>("later reservation", nc, nzt);
            Kokkos::deep_copy(later_zm2zt_api, real_t(57));
            compare("scoped returned zm2zt_api", scoped_zm2zt_api,
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 1));
            ddzt(workspace, nzm, nzt, nc, gr, t, filled_m);
            compare("scoped direct ddzt", filled_m,
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 4));
            auto scoped_ddzt = ddzt(workspace, nzm, nzt, nc, gr, t);
            // Later reservations must not reuse storage belonging to the returned field.
            auto later_ddzt = workspace_scope.view<Real2D>("later reservation", nc, nzm);
            Kokkos::deep_copy(later_ddzt, real_t(57));
            compare("scoped returned ddzt", scoped_ddzt,
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 4));
            ddzm(workspace, nzm, nzt, nc, gr, m, filled_t);
            compare("scoped direct ddzm", filled_t,
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 4));
            auto scoped_ddzm = ddzm(workspace, nzm, nzt, nc, gr, m);
            // Later reservations must not reuse storage belonging to the returned field.
            auto later_ddzm = workspace_scope.view<Real2D>("later reservation", nc, nzt);
            Kokkos::deep_copy(later_ddzm, real_t(57));
            compare("scoped returned ddzm", scoped_ddzm,
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 4));
            zt2zm2zt(workspace, nzm, nzt, nc, gr, t, filled_t, minimum);
            compare("scoped direct zt2zm2zt", filled_t,
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 3));
            auto scoped_zt2zm2zt = zt2zm2zt(workspace, nzm, nzt, nc, gr, t, minimum);
            // Later reservations must not reuse storage belonging to the returned field.
            auto later_zt2zm2zt = workspace_scope.view<Real2D>("later reservation", nc, nzt);
            Kokkos::deep_copy(later_zt2zm2zt, real_t(57));
            compare("scoped returned zt2zm2zt", scoped_zt2zm2zt,
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 3));
            zm2zt2zm(workspace, nzm, nzt, nc, gr, m, filled_m, minimum);
            compare("scoped direct zm2zt2zm", filled_m,
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 3));
            auto scoped_zm2zt2zm = zm2zt2zm(workspace, nzm, nzt, nc, gr, m, minimum);
            // Later reservations must not reuse storage belonging to the returned field.
            auto later_zm2zt2zm = workspace_scope.view<Real2D>("later reservation", nc, nzm);
            Kokkos::deep_copy(later_zm2zt2zm, real_t(57));
            compare("scoped returned zm2zt2zm", scoped_zm2zt2zm,
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 3));
            Host2D host_filled_m("scoped host output m", nc, nzm),
                host_filled_t("scoped host output t", nc, nzt);
            zt2zm_api(workspace, nzm, nzt, nc, gr, azt, host_filled_m, minimum);
            compare("scoped host output zt2zm_api", host_filled_m,
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 1));
            compare("scoped host returned zt2zm_api",
                    zt2zm_api(workspace, nzm, nzt, nc, gr, azt, minimum),
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 1));
            zm2zt_api(workspace, nzm, nzt, nc, gr, azm, host_filled_t, minimum);
            compare("scoped host output zm2zt_api", host_filled_t,
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 1));
            compare("scoped host returned zm2zt_api",
                    zm2zt_api(workspace, nzm, nzt, nc, gr, azm, minimum),
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 1));
            ddzt(workspace, nzm, nzt, nc, gr, azt, host_filled_m);
            compare("scoped host output ddzt", host_filled_m,
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 4));
            compare("scoped host returned ddzt", ddzt(workspace, nzm, nzt, nc, gr, azt),
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 4));
            ddzm(workspace, nzm, nzt, nc, gr, azm, host_filled_t);
            compare("scoped host output ddzm", host_filled_t,
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 4));
            compare("scoped host returned ddzm", ddzm(workspace, nzm, nzt, nc, gr, azm),
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 4));
            zt2zm2zt(workspace, nzm, nzt, nc, gr, azt, host_filled_t, minimum);
            compare("scoped host output zt2zm2zt", host_filled_t,
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 3));
            compare("scoped host returned zt2zm2zt",
                    zt2zm2zt(workspace, nzm, nzt, nc, gr, azt, minimum),
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 3));
            zm2zt2zm(workspace, nzm, nzt, nc, gr, azm, host_filled_m, minimum);
            compare("scoped host output zm2zt2zm", host_filled_m,
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 3));
            compare("scoped host returned zm2zt2zm",
                    zm2zt2zm(workspace, nzm, nzt, nc, gr, azm, minimum),
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 3));
            // Both smoothers permit an explicit output that aliases the original input.
            Kokkos::deep_copy(filled_m, m);
            zm2zt2zm(workspace, nzm, nzt, nc, gr, filled_m, filled_m, minimum);
            compare("scoped in-place momentum smoother", filled_m,
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 3));
            Kokkos::deep_copy(filled_t, t);
            zt2zm2zt(workspace, nzm, nzt, nc, gr, filled_t, filled_t, minimum);
            compare("scoped in-place thermo smoother", filled_t,
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 3));
            // Generic adapters stage genuinely strided host fields and preserve adjacent values.
            Host3D strided_source("scoped strided source", 2, nc, nzm);
            Host3D strided_result("scoped strided result", 2, nc, nzm);
            Kokkos::deep_copy(strided_result, real_t(93));
            auto source_window = Kokkos::subview(strided_source, 0, Kokkos::ALL, Kokkos::ALL);
            auto result_window = Kokkos::subview(strided_result, 0, Kokkos::ALL, Kokkos::ALL);
            Kokkos::deep_copy(source_window, azm);
            zm2zt2zm(workspace, nzm, nzt, nc, gr, source_window, result_window, minimum);
            compare("scoped strided host smoother", result_window,
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 3));
            compare("scoped strided host returned derivative",
                    ddzm(workspace, nzm, nzt, nc, gr, source_window),
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 4));
            for (int k = 0; k < nzm; ++k)
              for (int i = 0; i < nc; ++i)
                if (strided_result(1, i, k) != real_t(93))
                  throw std::runtime_error("scoped grid adapter overwrote adjacent storage");
            // Overlapping input/output must retain the original snapshot rather
            // than reading values overwritten by neighboring work items.
            Real2D overlap("overlap", nc, nzm);
            auto input_window = Kokkos::subview(overlap, Kokkos::ALL, Kokkos::make_pair(0, nzt));
            Kokkos::deep_copy(input_window, t);
            zt2zm_api(nzm, nzt, nc, gr, input_window, overlap);
            compare("zt2zm overlapping output", overlap,
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 0));
            Kokkos::deep_copy(input_window, t);
            ddzt(nzm, nzt, nc, gr, input_window, overlap);
            compare("ddzt overlapping output", overlap,
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 4));
            Kokkos::deep_copy(input_window, t);
            zt2zm_api(workspace, nzm, nzt, nc, gr, input_window, overlap);
            compare("scoped overlapping interpolation", overlap,
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 0));
            Kokkos::deep_copy(input_window, t);
            ddzt(workspace, nzm, nzt, nc, gr, input_window, overlap);
            compare("scoped overlapping derivative", overlap,
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 4));
            compare("host zt2zm", zt2zm_api(nzm, nzt, nc, gr, azt),
                    Kokkos::subview(output_m, Kokkos::ALL, Kokkos::ALL, 0));
            compare("host zm2zt", zm2zt_api(nzm, nzt, nc, gr, azm),
                    Kokkos::subview(output_t, Kokkos::ALL, Kokkos::ALL, 0));
            Real1D mp("mp", nzm), tp("tp", nzt);
            Kokkos::deep_copy(mp, Kokkos::subview(m, 0, Kokkos::ALL));
            Kokkos::deep_copy(tp, Kokkos::subview(t, 0, Kokkos::ALL));
            compare("zt2zm profile", zt2zm_api(gr, tp), Kokkos::subview(profile_m, Kokkos::ALL, 0));
            compare("zt2zm profile floor", zt2zm_api(gr, tp, minimum),
                    Kokkos::subview(profile_m, Kokkos::ALL, 1));
            compare("zm2zt profile", zm2zt_api(gr, mp), Kokkos::subview(profile_t, Kokkos::ALL, 0));
            compare("zm2zt profile floor", zm2zt_api(gr, mp, minimum),
                    Kokkos::subview(profile_t, Kokkos::ALL, 1));
            compare("ddzt profile", ddzt(gr, tp), Kokkos::subview(profile_m, Kokkos::ALL, 4));
            compare("ddzm profile", ddzm(gr, mp), Kokkos::subview(profile_t, Kokkos::ALL, 4));
            for (int k = 0; k < nzm; ++k) {
              real_t got = zt2zm_api(gr, tp, k);
              if (std::memcmp(&got, &profile_m(k, 2), sizeof(real_t)))
                throw std::runtime_error("momentum scalar mismatch");
              got = zt2zm_api(gr, tp, k, minimum);
              if (std::memcmp(&got, &profile_m(k, 3), sizeof(real_t)))
                throw std::runtime_error("momentum scalar floor mismatch");
            }
            for (int k = 0; k < nzt; ++k) {
              real_t got = zm2zt_api(gr, mp, k);
              if (std::memcmp(&got, &profile_t(k, 2), sizeof(real_t)))
                throw std::runtime_error("thermodynamic scalar mismatch");
              got = zm2zt_api(gr, mp, k, minimum);
              if (std::memcmp(&got, &profile_t(k, 3), sizeof(real_t)))
                throw std::runtime_error("thermodynamic scalar floor mismatch");
            }
            compare("unchanged momentum", m, azm);
            compare("unchanged thermo", t, azt);
            ++configurations;
          }
    for (int n : {0, 1, 7, 33}) {
      Host1D x("flip_input", n), reference("flip_reference", n);
      for (int k = 0; k < n; ++k)
        x(k) = k % 3 ? real_t(k) * real_t(0.137) : real_t(-0.0);
      f2k_grid_flip_c(n, x.data(), reference.data());
      compare("flip", flip(to_execution_space(x), n), reference);
    }
    Kokkos::fence();
    std::cout << configurations
              << " grid configurations match native bits for all operator interfaces\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
