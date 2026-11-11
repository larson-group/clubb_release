// Reject LAPACK solves; see lapack_wrap.hpp for the retained interfaces.
#include "lapack_wrap.hpp"
#include "error_code.hpp"

#include <iostream>

namespace clubb_kokkos {

void lapack_tridiag_solvex(const std::string& solve_type, [[maybe_unused]] int ndim,
                           [[maybe_unused]] int nrhs, [[maybe_unused]] int ngrdcol, // In
                           [[maybe_unused]] Real3D lhs, [[maybe_unused]] Real3D rhs,
                           err_info_type& err_info, // InOut
                           Real3D soln, Real1D rcond) { // Out
  // Reject the unsupported tridiagonal LAPACK solve without changing lhs/rhs.
  // Mark every column fatal and clear the solution and condition estimate.
  std::cerr << err_info.err_header_global << "\nlapack_tridiag_solvex: " << solve_type
            << "\nLAPACK is not supported in Kokkos.\n";
  Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
  Kokkos::deep_copy(soln, real_t(0));
  Kokkos::deep_copy(rcond, real_t(0));
}

void lapack_band_solvex(const std::string& solve_type, [[maybe_unused]] int nsup,
                        [[maybe_unused]] int nsub, [[maybe_unused]] int ndim,
                        [[maybe_unused]] int nrhs, [[maybe_unused]] int ngrdcol, // In
                        [[maybe_unused]] Real3D lhs, [[maybe_unused]] Real3D rhs,
                        err_info_type& err_info, // InOut
                        Real3D soln, Real1D rcond) { // Out
  // Reject the unsupported band LAPACK solve without changing lhs/rhs.
  // Mark every column fatal and clear the solution and condition estimate.
  std::cerr << err_info.err_header_global << "\nlapack_band_solvex: " << solve_type
            << "\nLAPACK is not supported in Kokkos.\n";
  Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
  Kokkos::deep_copy(soln, real_t(0));
  Kokkos::deep_copy(rcond, real_t(0));
}

} // namespace clubb_kokkos
