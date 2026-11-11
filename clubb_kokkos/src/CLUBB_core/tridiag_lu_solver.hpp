// Solve tridiagonal band systems using the native unpivoted LU algorithm.
//
// tridiag_lu_solve selects the native single/multiple-LHS and RHS overloads.
// Each column has an independent matrix. Multicolumn lhs uses (column, band,
// row), with native signed bands mapped by band_cpp = band_native + 1.
// This swaps the first two native dimensions so adjacent GPU threads access
// contiguous columns. rhs/soln use (column, row[, RHS]); a single matrix uses
// (band, row) and rank-one rhs/soln. lhs is declared InOut as in Fortran but
// these routines leave its values unchanged. Inputs must have nonzero LU pivots.
// ndim must be at least 2; no pivoting or threshold repair is applied.
#pragma once

#include "scratch_workspace.hpp"
#include "array_types.hpp"

namespace clubb_kokkos {

void tridiag_lu_solve_single_rhs_lhs(ScratchWorkspace& workspace,          // InOut
                                     int ndim,                             // In
                                     Real2D lhs,                           // InOut
                                     Real1D rhs,                           // In
                                     Real1D soln);                         // Out
void tridiag_lu_solve_single_rhs_multiple_lhs(ScratchWorkspace& workspace, // InOut
                                              int ndim, int ngrdcol,       // In
                                              Real3D lhs,                  // InOut
                                              Real2D rhs,                  // In
                                              Real2D soln);                // Out
void tridiag_lu_solve_multiple_rhs_lhs(ScratchWorkspace& workspace,        // InOut
                                       int ndim, int nrhs, int ngrdcol,    // In
                                       Real3D lhs,                         // InOut
                                       Real3D rhs,                         // In
                                       Real3D soln);                       // Out

// C++ overloads provide the native generic interface without changing its arguments.
inline void tridiag_lu_solve(ScratchWorkspace& workspace, // InOut
                             int ndim,                    // In
                             Real2D lhs,                  // InOut
                             Real1D rhs,                  // In
                             Real1D soln) {               // Out
  // Select the native rank-specific solve.
  auto workspace_scope = workspace.scope();
  tridiag_lu_solve_single_rhs_lhs(workspace, ndim, lhs, rhs, soln);
}

inline void tridiag_lu_solve(ScratchWorkspace& workspace, // InOut
                             int ndim, int ngrdcol,       // In
                             Real3D lhs,                  // InOut
                             Real2D rhs,                  // In
                             Real2D soln) {               // Out
  // Select the native rank-specific solve.
  auto workspace_scope = workspace.scope();
  tridiag_lu_solve_single_rhs_multiple_lhs(workspace, ndim, ngrdcol, lhs, rhs, soln);
}

inline void tridiag_lu_solve(ScratchWorkspace& workspace,     // InOut
                             int ndim, int nrhs, int ngrdcol, // In
                             Real3D lhs,                      // InOut
                             Real3D rhs,                      // In
                             Real3D soln) {                   // Out
  // Select the native rank-specific solve.
  auto workspace_scope = workspace.scope();
  tridiag_lu_solve_multiple_rhs_lhs(workspace, ndim, nrhs, ngrdcol, lhs, rhs, soln);
}

} // namespace clubb_kokkos
