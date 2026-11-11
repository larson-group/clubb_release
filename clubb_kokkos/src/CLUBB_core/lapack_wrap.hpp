// Reject unsupported LAPACK solves through the retained solver interfaces.
//
// lapack_tridiag_solvex and lapack_band_solvex set every err_info column to
// clubb_fatal_error and print that LAPACK is not supported in Kokkos. They
// preserve lhs/rhs and clear the solution and reciprocal condition estimate.
// C++ matrices use (column, band, row), and rhs/soln use (column, row, RHS).
// Optional condition diagnostics are unsupported along with LAPACK solvers.
#pragma once

#include "array_types.hpp"
#include "err_info_type_module.hpp"
#include <string>

namespace clubb_kokkos {

void lapack_tridiag_solvex(const std::string& solve_type, int ndim, int nrhs, int ngrdcol, // In
                           Real3D lhs, Real3D rhs, err_info_type& err_info,                // InOut
                           Real3D soln, Real1D rcond);                                     // Out

void lapack_band_solvex(const std::string& solve_type, int nsup, int nsub, // In
                        int ndim, int nrhs, int ngrdcol,                   // In
                        Real3D lhs, Real3D rhs, err_info_type& err_info,   // InOut
                        Real3D soln, Real1D rcond);                        // Out

} // namespace clubb_kokkos
