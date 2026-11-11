// Solve pentadiagonal band systems using the native unpivoted LU algorithm.
//
// penta_lu_solve selects the native single/multiple-RHS overloads.
// Each column has an independent matrix. Multicolumn lhs uses (column, band,
// row), with native signed bands mapped by band_cpp = band_native + 2.
// This swaps the first two native dimensions so adjacent GPU threads access
// contiguous columns. rhs/soln use (column, row[, RHS]).
// lhs is declared InOut as in Fortran but
// these routines leave its values unchanged. Inputs must have nonzero LU pivots.
// ndim must be at least 4; no pivoting or threshold repair is applied.
#pragma once

#include "scratch_workspace.hpp"
#include "array_types.hpp"

namespace clubb_kokkos {

void penta_lu_solve_single_rhs_multiple_lhs(ScratchWorkspace& workspace, // InOut
                                            int ndim, int ngrdcol,       // In
                                            Real3D lhs,                  // InOut
                                            Real2D rhs,                  // In
                                            Real2D soln);                // Out
void penta_lu_solve_multiple_rhs_lhs(ScratchWorkspace& workspace,        // InOut
                                     int ndim, int nrhs, int ngrdcol,    // In
                                     Real3D lhs,                         // InOut
                                     Real3D rhs,                         // In
                                     Real3D soln);                       // Out

// C++ overloads provide the native generic interface without changing its arguments.
inline void penta_lu_solve(ScratchWorkspace& workspace, // InOut
                           int ndim, int ngrdcol,       // In
                           Real3D lhs,                  // InOut
                           Real2D rhs,                  // In
                           Real2D soln) {               // Out
  // Select the native rank-specific solve.
  auto workspace_scope = workspace.scope();
  penta_lu_solve_single_rhs_multiple_lhs(workspace, ndim, ngrdcol, lhs, rhs, soln);
}

inline void penta_lu_solve(ScratchWorkspace& workspace,     // InOut
                           int ndim, int nrhs, int ngrdcol, // In
                           Real3D lhs,                      // InOut
                           Real3D rhs,                      // In
                           Real3D soln) {                   // Out
  // Select the native rank-specific solve.
  auto workspace_scope = workspace.scope();
  penta_lu_solve_multiple_rhs_lhs(workspace, ndim, nrhs, ngrdcol, lhs, rhs, soln);
}

} // namespace clubb_kokkos
