// Supported matrix solve dispatch; see matrix_solver_wrapper.hpp for interfaces.
#include "scratch_workspace.hpp"
#include "matrix_solver_wrapper.hpp"
#include "penta_lu_solver.hpp"
#include "tridiag_lu_solver.hpp"
#include "model_flags.hpp"
#include "error_code.hpp"

#include <cmath>
#include <iostream>

namespace clubb_kokkos {

static void check_nan_soln_1d(const std::string& solve_name, int solve_method, int ndim,
                              Real1D soln,              // In
                              err_info_type& err_info); // InOut
static void check_nan_soln_2d(const std::string& solve_name, int solve_method, int ngrdcol,
                              int ndim, Real2D soln,    // In
                              err_info_type& err_info); // InOut
static void check_nan_soln_3d(const std::string& solve_name, int solve_method, int ngrdcol,
                              int ndim, int nrhs, Real3D soln, // In
                              err_info_type& err_info);        // InOut

void band_solve_single_rhs_multiple_lhs(ScratchWorkspace& workspace, const std::string& solve_name,
                                        int penta_solve_method, int ngrdcol, int nsup, int nsub,
                                        int ndim,
                                        bool l_implemented,                              // In
                                        Real3D lhs, Real2D rhs, err_info_type& err_info, // InOut
                                        Real2D soln,                                     // Out
                                        [[maybe_unused]] const Real2D* old_soln, // Optional In
                                        Real1D* rcond) {                         // Optional Out
  // Solve a pentadiagonal system with one RHS for each column.
  // A LAPACK method or condition estimate request returns a fatal error.

  // ----------------------- Begin Code -----------------------
  // LAPACK solvers and condition diagnostics are unsupported in Kokkos.
  auto workspace_scope = workspace.scope();
  if (rcond || penta_solve_method == lapack) {
    std::cerr << err_info.err_header_global << "\n" << solve_name
              << "\nLAPACK is not supported in Kokkos.\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    Kokkos::deep_copy(soln, real_t(0));
    if (rcond)
      Kokkos::deep_copy(*rcond, real_t(0));
    return;
  }

  // JAX-supported band solves are pentadiagonal; reject other widths before
  // calling the five-band LU routine instead of reinterpreting its storage.
  if (nsup != 2 || nsub != 2) {
    std::cerr << "WARNING: matrix_solver_wrapper.band_solve_single_rhs_multiple_lhs only supports "
                 "pentadiagonal systems; got nsup="
              << nsup << ", nsub=" << nsub << ".\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    Kokkos::deep_copy(soln, real_t(0));
    return;
  }

  if (penta_solve_method == penta_lu) {
    // Solve the system with a penta-diagonal specific LU decomp.
    penta_lu_solve(workspace, ndim, ngrdcol, // In
                   lhs,                      // InOut
                   rhs,                      // In
                   soln);                    // Out
  } else if (penta_solve_method == penta_bicgstab) {
    // The native iterative path uses old_soln one RHS at a time because
    // RHS magnitudes vary greatly. JAX leaves that solver unsupported.
    std::cerr << "WARNING: matrix_solver_wrapper.band_solve_single_rhs_multiple_lhs does not "
                 "support the penta_bicgstab solver path.\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    Kokkos::deep_copy(soln, real_t(0));
    return;
  } else {
    // The solve method should match one of the above.
    if (clubb_at_least_debug_level_api(0)) {
      std::cerr << err_info.err_header_global << "\nError in band_solve_single_rhs_multiple_lhs: \n"
                << "  no case for penta_solve_method = " << penta_solve_method << '\n';
      // General error -> set all entries to clubb_fatal_error.
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    }
    return;
  }

  if (clubb_at_least_debug_level_api(1) || l_implemented) {
    check_nan_soln_2d(solve_name, penta_solve_method, ngrdcol, ndim, soln, // In
                      err_info);                                           // InOut
  }
}

void band_solve_multiple_rhs_lhs(ScratchWorkspace& workspace, const std::string& solve_name,
                                 int penta_solve_method, int ngrdcol, int nsup, int nsub, int ndim,
                                 int nrhs, bool l_implemented,                    // In
                                 Real3D lhs, Real3D rhs, err_info_type& err_info, // InOut
                                 Real3D soln,                                     // Out
                                 [[maybe_unused]] const Real3D* old_soln,         // Optional In
                                 Real1D* rcond) {                                 // Optional Out
  // Solve a pentadiagonal system with multiple RHSs for each column.
  // A LAPACK method or condition estimate request returns a fatal error.

  // ----------------------- Begin Code -----------------------
  // LAPACK solvers and condition diagnostics are unsupported in Kokkos.
  auto workspace_scope = workspace.scope();
  if (rcond || penta_solve_method == lapack) {
    std::cerr << err_info.err_header_global << "\n" << solve_name
              << "\nLAPACK is not supported in Kokkos.\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    Kokkos::deep_copy(soln, real_t(0));
    if (rcond)
      Kokkos::deep_copy(*rcond, real_t(0));
    return;
  }

  // JAX-supported band solves are pentadiagonal; reject other widths before
  // calling the five-band LU routine instead of reinterpreting its storage.
  if (nsup != 2 || nsub != 2) {
    std::cerr << "WARNING: matrix_solver_wrapper.band_solve_multiple_rhs_lhs only supports "
                 "pentadiagonal systems; got nsup="
              << nsup << ", nsub=" << nsub << ".\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    Kokkos::deep_copy(soln, real_t(0));
    return;
  }

  if (penta_solve_method == penta_lu) {
    // Solve the system with a penta-diagonal specific LU decomp.
    penta_lu_solve(workspace, ndim, nrhs, ngrdcol, // In
                   lhs,                            // InOut
                   rhs,                            // In
                   soln);                          // Out
  } else if (penta_solve_method == penta_bicgstab) {
    // The native iterative path uses old_soln one RHS at a time because
    // RHS magnitudes vary greatly. JAX leaves that solver unsupported.
    std::cerr << "WARNING: matrix_solver_wrapper.band_solve_multiple_rhs_lhs does not support the "
                 "penta_bicgstab solver path.\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    Kokkos::deep_copy(soln, real_t(0));
    return;
  } else {
    // The solve method should match one of the above.
    if (clubb_at_least_debug_level_api(0)) {
      std::cerr << err_info.err_header_global << "\nError in band_solve_multiple_rhs_lhs: \n"
                << "  no case for penta_solve_method = " << penta_solve_method << '\n';
      // General error -> set all entries to clubb_fatal_error.
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    }
    return;
  }

  if (clubb_at_least_debug_level_api(1) || l_implemented) {
    check_nan_soln_3d(solve_name, penta_solve_method, ngrdcol, ndim, nrhs, soln, // In
                      err_info);                                                 // InOut
  }
}

void tridiag_solve_single_rhs_lhs(ScratchWorkspace& workspace, const std::string& solve_name,
                                  int tridiag_solve_method,
                                  int ndim,                                        // In
                                  Real2D lhs, Real1D rhs, err_info_type& err_info, // InOut
                                  Real1D soln,                                     // Out
                                  Real1D* rcond) {                                 // Optional Out
  // Solve a tridiagonal system with one RHS and one matrix.
  // A LAPACK method or condition estimate request returns a fatal error.

  // ----------------------- Begin Code -----------------------
  // LAPACK solvers and condition diagnostics are unsupported in Kokkos.
  auto workspace_scope = workspace.scope();
  if (rcond || tridiag_solve_method == lapack) {
    std::cerr << err_info.err_header_global << "\n" << solve_name
              << "\nLAPACK is not supported in Kokkos.\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    Kokkos::deep_copy(soln, real_t(0));
    if (rcond)
      Kokkos::deep_copy(*rcond, real_t(0));
    return;
  }

  if (tridiag_solve_method == tridiag_lu) {
    // Solve the system with a tri-diagonal specific LU decomp.
    tridiag_lu_solve(workspace, ndim, // In
                     lhs,             // InOut
                     rhs,             // In
                     soln);           // Out
  } else {
    // The solve method should match one of the above.
    if (clubb_at_least_debug_level_api(0)) {
      std::cerr << err_info.err_header_global << "\nError in tridiag_solve_single_rhs_lhs: \n"
                << "  no case for tridiag_solve_method = " << tridiag_solve_method << '\n';
      // General error -> set all entries to clubb_fatal_error.
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    }
    return;
  }

  if (clubb_at_least_debug_level_api(1)) {
    check_nan_soln_1d(solve_name, tridiag_solve_method, ndim, soln, // In
                      err_info);                                    // InOut
  }
}

void tridiag_solve_single_rhs_multiple_lhs(ScratchWorkspace& workspace,
                                           const std::string& solve_name, int tridiag_solve_method,
                                           int ngrdcol, int ndim, bool l_implemented,       // In
                                           Real3D lhs, Real2D rhs, err_info_type& err_info, // InOut
                                           Real2D soln,                                     // Out
                                           Real1D* rcond) { // Optional Out
  // Solve a tridiagonal system with one RHS for each column.
  // A LAPACK method or condition estimate request returns a fatal error.

  // ----------------------- Begin Code -----------------------
  // LAPACK solvers and condition diagnostics are unsupported in Kokkos.
  auto workspace_scope = workspace.scope();
  if (rcond || tridiag_solve_method == lapack) {
    std::cerr << err_info.err_header_global << "\n" << solve_name
              << "\nLAPACK is not supported in Kokkos.\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    Kokkos::deep_copy(soln, real_t(0));
    if (rcond)
      Kokkos::deep_copy(*rcond, real_t(0));
    return;
  }

  if (tridiag_solve_method == tridiag_lu) {
    // Solve the system with a tri-diagonal specific LU decomp.
    tridiag_lu_solve(workspace, ndim, ngrdcol, // In
                     lhs,                      // InOut
                     rhs,                      // In
                     soln);                    // Out
  } else {
    // The solve method should match one of the above.
    if (clubb_at_least_debug_level_api(0)) {
      std::cerr << err_info.err_header_global
                << "\nError in tridiag_solve_single_rhs_multiple_lhs: \n"
                << "  no case for tridiag_solve_method = " << tridiag_solve_method << '\n';
      // General error -> set all entries to clubb_fatal_error.
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    }
    return;
  }

  if (clubb_at_least_debug_level_api(1) || l_implemented) {
    check_nan_soln_2d(solve_name, tridiag_solve_method, ngrdcol, ndim, soln, // In
                      err_info);                                             // InOut
  }
}

void tridiag_solve_multiple_rhs_lhs(ScratchWorkspace& workspace, const std::string& solve_name,
                                    int tridiag_solve_method, int ngrdcol, int ndim, int nrhs,
                                    bool l_implemented,                              // In
                                    Real3D lhs, Real3D rhs, err_info_type& err_info, // InOut
                                    Real3D soln,                                     // Out
                                    Real1D* rcond) {                                 // Optional Out
  // Solve a tridiagonal system with multiple RHSs for each column.
  // A LAPACK method or condition estimate request returns a fatal error.

  // ----------------------- Begin Code -----------------------
  // LAPACK solvers and condition diagnostics are unsupported in Kokkos.
  auto workspace_scope = workspace.scope();
  if (rcond || tridiag_solve_method == lapack) {
    std::cerr << err_info.err_header_global << "\n" << solve_name
              << "\nLAPACK is not supported in Kokkos.\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    Kokkos::deep_copy(soln, real_t(0));
    if (rcond)
      Kokkos::deep_copy(*rcond, real_t(0));
    return;
  }

  if (tridiag_solve_method == tridiag_lu) {
    // Solve the system with a tri-diagonal specific LU decomp.
    tridiag_lu_solve(workspace, ndim, nrhs, ngrdcol, // In
                     lhs,                            // InOut
                     rhs,                            // In
                     soln);                          // Out
  } else {
    // The solve method should match one of the above.
    if (clubb_at_least_debug_level_api(0)) {
      std::cerr << err_info.err_header_global << "\nError in tridiag_solve_multiple_rhs_lhs: \n"
                << "  no case for tridiag_solve_method = " << tridiag_solve_method << '\n';
      // General error -> set all entries to clubb_fatal_error.
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    }
    return;
  }

  if (clubb_at_least_debug_level_api(1) || l_implemented) {
    check_nan_soln_3d(solve_name, tridiag_solve_method, ngrdcol, ndim, nrhs, soln, // In
                      err_info);                                                   // InOut
  }
}

static void check_nan_soln_1d(const std::string& solve_name, int solve_method, int ndim,
                              Real1D soln,               // In
                              err_info_type& err_info) { // InOut
  // Report every NaN at its native one-based coordinates, then mark every
  // column fatal if any is present. A valid solution leaves prior status intact.
  bool l_has_nan = false;
  int nan_count  = 0; // Integer OR/count reduction; no floating-point reordering.

  // Detect on the device; copy values for host printing only when a NaN is found.
  Kokkos::parallel_reduce(
      "check_nan_soln_1d", Kokkos::RangePolicy<>(0, ndim),
      KOKKOS_LAMBDA(const int cell, int& count) {
        if (Kokkos::isnan(soln(cell)))
          ++count;
      },
      nan_count);
  l_has_nan = nan_count != 0;

  if (l_has_nan) {
    const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, soln);
    for (int k = 0; k < ndim; ++k)
      if (std::isnan(host(k)))
        std::cerr << "ERROR: NaN detected: soln(" << k + 1 << ") = " << host(k) << '\n';
    std::cerr << err_info.err_header_global << "\nsolve_method =" << solve_method << '\n'
              << solve_name << " NaNs in solution\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
  }
}

static void check_nan_soln_2d(const std::string& solve_name, int solve_method, int ngrdcol,
                              int ndim, Real2D soln,     // In
                              err_info_type& err_info) { // InOut
  // Report every NaN at its native one-based coordinates, then mark every
  // column fatal if any is present. A valid solution leaves prior status intact.
  bool l_has_nan = false;
  int nan_count  = 0; // Integer OR/count reduction; no floating-point reordering.

  // Detect on the device; copy values for host printing only when a NaN is found.
  Kokkos::parallel_reduce(
      "check_nan_soln_2d", Kokkos::RangePolicy<>(0, ngrdcol * ndim),
      KOKKOS_LAMBDA(const int cell, int& count) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        if (Kokkos::isnan(soln(i, k)))
          ++count;
      },
      nan_count);
  l_has_nan = nan_count != 0;

  if (l_has_nan) {
    const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, soln);
    for (int k = 0; k < ndim; ++k)
      for (int i = 0; i < ngrdcol; ++i)
        if (std::isnan(host(i, k)))
          std::cerr << "ERROR: NaN detected: soln(" << i + 1 << "," << k + 1 << ") = " << host(i, k)
                    << '\n';
    std::cerr << err_info.err_header_global << "\nsolve_method =" << solve_method << '\n'
              << solve_name << " NaNs in solution\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
  }
}

static void check_nan_soln_3d(const std::string& solve_name, int solve_method, int ngrdcol,
                              int ndim, int nrhs, Real3D soln, // In
                              err_info_type& err_info) {       // InOut
  // Report every NaN at its native one-based coordinates, then mark every
  // column fatal if any is present. A valid solution leaves prior status intact.
  bool l_has_nan = false;
  int nan_count  = 0; // Integer OR/count reduction; no floating-point reordering.

  // Detect on the device; copy values for host printing only when a NaN is found.
  Kokkos::parallel_reduce(
      "check_nan_soln_3d", Kokkos::RangePolicy<>(0, ngrdcol * ndim * nrhs),
      KOKKOS_LAMBDA(const int cell, int& count) {
        const int i = cell % ngrdcol;
        const int k = (cell / ngrdcol) % ndim;
        const int j = cell / (ngrdcol * ndim);
        if (Kokkos::isnan(soln(i, k, j)))
          ++count;
      },
      nan_count);
  l_has_nan = nan_count != 0;

  if (l_has_nan) {
    const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, soln);
    for (int j = 0; j < nrhs; ++j)
      for (int k = 0; k < ndim; ++k)
        for (int i = 0; i < ngrdcol; ++i)
          if (std::isnan(host(i, k, j)))
            std::cerr << "ERROR: NaN detected: soln(" << i + 1 << "," << k + 1 << "," << j + 1
                      << ") = " << host(i, k, j) << '\n';
    std::cerr << err_info.err_header_global << "\nsolve_method =" << solve_method << '\n'
              << solve_name << " NaNs in solution\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
  }
}

} // namespace clubb_kokkos
