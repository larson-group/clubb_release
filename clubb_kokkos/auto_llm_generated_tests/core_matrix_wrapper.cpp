// Compare all five matrix-wrapper solve surfaces with compiled native Fortran.
// LAPACK methods and optional conditions must fail before solving. NaN checks
// preserve debug/host gates and mark all columns fatal, including healthy ones.
// JAX-excluded methods and unsupported band widths must fail without fallback.
#include "scratch_workspace_owner.hpp"
#include "CLUBB_core/matrix_solver_wrapper.hpp"
#include "CLUBB_core/error_code.hpp"
#include "CLUBB_core/model_flags.hpp"

#include <array>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using clubb_kokkos::real_t;
extern "C" void f2k_matrix_reference_c(int, int, int, int, int, int, int, int, real_t*, real_t*,
                                       real_t*, real_t*, int*);

void solve_cpp(clubb_kokkos::ScratchWorkspace& workspace, int variant, int method, bool implemented,
               clubb_kokkos::Real3D lhs, clubb_kokkos::Real3D rhs,
               clubb_kokkos::err_info_type& err_info, clubb_kokkos::Real3D soln,
               clubb_kokkos::Real1D* rcond) {
  // Select one real public C++ overload; singleton Views reproduce the native
  // one-matrix/profile argument association without changing any array values.
  auto workspace_scope = workspace.scope();
  using namespace clubb_kokkos;
  const int nc = rhs.extent(0), ndim = rhs.extent(1), nrhs = rhs.extent(2);
  const int nsup = lhs.extent(1) / 2;
  switch (variant) {
  case 1:
    band_solve(workspace, "matrix test", method, nc, nsup, nsup, ndim, implemented, lhs,
               Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, 0), err_info,
               Kokkos::subview(soln, Kokkos::ALL, Kokkos::ALL, 0), nullptr, rcond);
    break;
  case 2:
    band_solve(workspace, "matrix test", method, nc, nsup, nsup, ndim, nrhs, implemented, lhs, rhs,
               err_info, soln, nullptr, rcond);
    break;
  case 3:
    if (nc != 1)
      throw std::runtime_error("single-matrix fixture has multiple columns");
    tridiag_solve(workspace, "matrix test", method, ndim, Real2D(lhs.data(), 3, ndim),
                  Real1D(rhs.data(), ndim), err_info, Real1D(soln.data(), ndim), rcond);
    break;
  case 4:
    tridiag_solve(workspace, "matrix test", method, nc, ndim, implemented, lhs,
                  Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, 0), err_info,
                  Kokkos::subview(soln, Kokkos::ALL, Kokkos::ALL, 0), rcond);
    break;
  case 5:
    tridiag_solve(workspace, "matrix test", method, nc, ndim, nrhs, implemented, lhs, rhs, err_info,
                  soln, rcond);
    break;
  default:
    throw std::runtime_error("invalid wrapper fixture variant");
  }
}

int main(int argc, char* argv[]) {
  // Compare all defined finite outputs and exact statuses against source code,
  // separately checking scope rejections that intentionally differ from native.
  try {
    Kokkos::ScopeGuard kokkos(argc, argv);
    clubb_kokkos::test::ScratchWorkspaceOwner scratch_owner;
    auto& workspace = scratch_owner.workspace;
    using namespace clubb_kokkos;
    using Host3D = Kokkos::View<real_t***, Kokkos::LayoutLeft, Kokkos::HostSpace>;
    using Host1D = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
    const std::array<std::array<int, 2>, 3> gates{{{0, 0}, {0, 1}, {1, 0}}};
    for (int variant = 1; variant <= 5; ++variant)
      for (int ndim : {4, 9})
        for (int inject_nan : {0, 1})
          for (auto gate : gates)
            for (int method : {2, 7}) {
              // One unknown-method case per surface verifies native diagnostics
              // without inspecting undefined soln.
              if (method == 7 && (ndim != 4 || inject_nan || gate != gates[0]))
                continue;
              const int nc     = variant == 3 ? 1 : 4;
              const int nrhs   = (variant == 2 || variant == 5) ? 3 : 1;
              const int nbands = variant <= 2 ? 5 : 3, diagonal = nbands / 2;
              Host3D lhs("lhs", nc, nbands, ndim), native_lhs("native_lhs", nbands, nc, ndim);
              Host3D rhs("rhs", nc, ndim, nrhs), native_rhs("native_rhs", nc, ndim, nrhs);
              Host3D native_soln("native_soln", nc, ndim, nrhs);
              Host1D native_rcond("native_rcond", nc);
              std::vector<int> native_codes(nc);
              for (int k = 0; k < ndim; ++k)
                for (int i = 0; i < nc; ++i) {
                  for (int b = 0; b < nbands; ++b) {
                    const int col = k - (b - diagonal);
                    lhs(i, b, k)  = (b == diagonal) ? real_t(4) + real_t(0.01) * real_t(i + k)
                                    : (col >= 0 && col < ndim)
                                        ? real_t(-0.013) * real_t(1 + (i + b + k) % 7)
                                        : real_t(9876.5);
                    native_lhs(b, i, k) = lhs(i, b, k);
                  }
                  for (int n = 0; n < nrhs; ++n)
                    native_rhs(i, k, n) = rhs(i, k, n) = real_t(0.173) * real_t(i + k + n + 1);
                }
              if (inject_nan)
                native_rhs(nc - 1, 1, nrhs - 1) = rhs(nc - 1, 1, nrhs - 1) =
                    std::numeric_limits<real_t>::quiet_NaN();
              auto lhs_d = to_execution_space(lhs), rhs_d = to_execution_space(rhs);
              Real3D soln("soln", nc, ndim, nrhs);
              Real1D rcond("rcond", nc);
              Kokkos::deep_copy(rcond, real_t(1234));
              err_info_type err_info;
              err_info.err_code = decltype(err_info.err_code)("codes", nc);
              auto codes_seed   = Kokkos::create_mirror_view(err_info.err_code);
              for (int i = 0; i < nc; ++i)
                codes_seed(i) = i == 0 ? 50 : 0;
              Kokkos::deep_copy(err_info.err_code, codes_seed);
              err_info.err_header_global = "matrix wrapper test";
              // The reference sets its own Native threshold; set the independent
              // C++ owner to the same gate before comparing the two solvers.
              set_clubb_debug_level_api(gate[0]);
              f2k_matrix_reference_c(variant, method, ndim, nrhs, nc, 0, gate[0], gate[1],
                                     native_lhs.data(), native_rhs.data(), native_soln.data(),
                                     native_rcond.data(), native_codes.data());
              solve_cpp(workspace, variant, method, gate[1], lhs_d, rhs_d, err_info, soln, nullptr);
              auto soln_h  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, soln);
              auto lhs_h   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, lhs_d);
              auto rhs_h   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, rhs_d);
              auto rcond_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, rcond);
              auto codes_h =
                  Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, err_info.err_code);
              for (int i = 0; i < nc; ++i) {
                if (codes_h(i) != native_codes[i])
                  throw std::runtime_error("matrix-wrapper debug/host/status mismatch");
                if (std::memcmp(&rcond_h(i), &native_rcond(i), sizeof(real_t)))
                  throw std::runtime_error("matrix-wrapper rcond mismatch");
              }
              for (int k = 0; k < ndim; ++k)
                for (int i = 0; i < nc; ++i)
                  for (int b = 0; b < nbands; ++b)
                    if (std::memcmp(&lhs_h(i, b, k), &native_lhs(b, i, k), sizeof(real_t)))
                      throw std::runtime_error(
                          "matrix-wrapper changed actual lhs during diagnostics");
              if (std::memcmp(rhs_h.data(), native_rhs.data(), rhs.size() * sizeof(real_t)))
                throw std::runtime_error("matrix-wrapper changed actual rhs during diagnostics");
              if (method == 2)
                for (int n = 0; n < nrhs; ++n)
                  for (int k = 0; k < ndim; ++k)
                    for (int i = 0; i < nc; ++i) {
                      const real_t got = soln_h(i, k, n), expected = native_soln(i, k, n);
                      if (std::isnan(got) != std::isnan(expected) ||
                          (!std::isnan(expected) && std::memcmp(&got, &expected, sizeof(real_t))))
                        throw std::runtime_error("matrix-wrapper solution mismatch");
                    }
            }

    for (int debug_level : {0, 2}) {
      set_clubb_debug_level_api(debug_level);
      for (int variant = 1; variant <= 5; ++variant)
        for (int request = 0; request < 3; ++request) {
          // Reject both explicit LAPACK and LU-plus-condition requests, with
          // and without debug checks, without reaching any LU factorization.
          const int nc       = variant == 3 ? 1 : 3, ndim = 4;
          const int nrhs     = (variant == 2 || variant == 5) ? 3 : 1;
          const bool use_rcond = request != 0;
          const int method   = request == 1 ? 2 : lapack;
          Real3D lhs("unsupported_lhs", nc, variant <= 2 ? 5 : 3, ndim);
          Real3D rhs("rhs", nc, ndim, nrhs), soln("soln", nc, ndim, nrhs);
          Real1D rcond("rcond", nc);
          Kokkos::deep_copy(lhs, real_t(42));
          Kokkos::deep_copy(rhs, real_t(17));
          Kokkos::deep_copy(soln, real_t(9876.5));
          Kokkos::deep_copy(rcond, real_t(1234));
          err_info_type err_info;
          err_info.err_code = decltype(err_info.err_code)("codes", nc);
          Kokkos::deep_copy(err_info.err_code, 50);
          solve_cpp(workspace, variant, method, false, lhs, rhs, err_info, soln,
                    use_rcond ? &rcond : nullptr);
          const auto codes_h = to_host_space(err_info.err_code);
          const auto soln_h  = to_host_space(soln);
          const auto lhs_h   = to_host_space(lhs);
          const auto rhs_h   = to_host_space(rhs);
          const auto rcond_h = to_host_space(rcond);
          for (int i = 0; i < nc; ++i) {
            if (codes_h(i) != clubb_fatal_error)
              throw std::runtime_error("LAPACK request was not fatal for every column");
            if (rcond_h(i) != (use_rcond ? real_t(0) : real_t(1234)))
              throw std::runtime_error("unsupported condition estimate was not cleared");
          }
          for (std::size_t cell = 0; cell < soln_h.size(); ++cell)
            if (soln_h.data()[cell] != real_t(0) || rhs_h.data()[cell] != real_t(17))
              throw std::runtime_error("LAPACK request solved or changed rhs");
          for (std::size_t cell = 0; cell < lhs_h.size(); ++cell)
            if (lhs_h.data()[cell] != real_t(42))
              throw std::runtime_error("LAPACK request changed lhs");
        }
    }
    set_clubb_debug_level_api(0);
    for (int variant : {1, 2}) {
      const int nc = 3, ndim = 4, nrhs = variant == 2 ? 3 : 1;
      Real3D lhs("unsupported_lhs", nc, 5, ndim);
      Real3D rhs("rhs", nc, ndim, nrhs), soln("soln", nc, ndim, nrhs);
      Kokkos::deep_copy(soln, real_t(9876.5));
      err_info_type err_info;
      err_info.err_code = decltype(err_info.err_code)("codes", nc);
      solve_cpp(workspace, variant, penta_bicgstab, false, lhs, rhs, err_info, soln, nullptr);
      auto codes_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, err_info.err_code);
      auto soln_h  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, soln);
      for (int i = 0; i < nc; ++i)
        if (codes_h(i) != clubb_fatal_error)
          throw std::runtime_error("excluded solver was not rejected");
      for (std::size_t cell = 0; cell < soln_h.size(); ++cell)
        if (soln_h.data()[cell] != real_t(0))
          throw std::runtime_error("excluded solver returned a solution");
    }
    for (int variant : {1, 2}) {
      Real3D lhs("wrong_bands", 2, 3, 4), rhs("rhs", 2, 4, variant == 2 ? 3 : 1);
      Real3D soln("soln", 2, 4, variant == 2 ? 3 : 1);
      err_info_type err_info;
      err_info.err_code = decltype(err_info.err_code)("codes", 2);
      solve_cpp(workspace, variant, penta_lu, false, lhs, rhs, err_info, soln, nullptr);
      auto codes_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, err_info.err_code);
      if (codes_h(0) != clubb_fatal_error || codes_h(1) != clubb_fatal_error)
        throw std::runtime_error("non-pentadiagonal band solve was not rejected");
    }
    std::cout
        << "Supported LU interfaces match native; LAPACK requests fail before solving\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
