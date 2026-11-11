// Dispatch supported band solves and preserve native error diagnostics.
//
// band_solve handles column-wise pentadiagonal systems; tridiag_solve also
// handles the single-matrix form. C++ lhs uses (column,band,row), with bands
// shifted +2/+1 from native signed bands; a single lhs is (band,row). RHS/soln
// use (column,row[,RHS]) or a rank-one single-matrix profile. LAPACK methods
// and optional rcond requests return fatal err_info before any solve runs.
// BiCGStab methods and non-pentadiagonal band choices are also rejected.
// Native debug/host gates and all-column NaN fatal status
// are preserved. old_soln remains the native optional input for the excluded
// BiCGStab method. All array arguments are execution-space Views.
#pragma once

#include "scratch_workspace.hpp"
#include "array_types.hpp"
#include "err_info_type_module.hpp"
#include <string>

namespace clubb_kokkos {

void band_solve_single_rhs_multiple_lhs(ScratchWorkspace& workspace, const std::string& solve_name,
                                        int penta_solve_method, int ngrdcol, int nsup, int nsub,
                                        int ndim,
                                        bool l_implemented,                              // In
                                        Real3D lhs, Real2D rhs, err_info_type& err_info, // InOut
                                        Real2D soln,                                     // Out
                                        const Real2D* old_soln = nullptr,  // Optional In
                                        Real1D* rcond          = nullptr); // Optional Out

void band_solve_multiple_rhs_lhs(ScratchWorkspace& workspace, const std::string& solve_name,
                                 int penta_solve_method, int ngrdcol, int nsup, int nsub, int ndim,
                                 int nrhs, bool l_implemented,                    // In
                                 Real3D lhs, Real3D rhs, err_info_type& err_info, // InOut
                                 Real3D soln,                                     // Out
                                 const Real3D* old_soln = nullptr,                // Optional In
                                 Real1D* rcond          = nullptr);               // Optional Out

void tridiag_solve_single_rhs_lhs(ScratchWorkspace& workspace, const std::string& solve_name,
                                  int tridiag_solve_method,
                                  int ndim,                                        // In
                                  Real2D lhs, Real1D rhs, err_info_type& err_info, // InOut
                                  Real1D soln,                                     // Out
                                  Real1D* rcond = nullptr);                        // Optional Out

void tridiag_solve_single_rhs_multiple_lhs(ScratchWorkspace& workspace,
                                           const std::string& solve_name, int tridiag_solve_method,
                                           int ngrdcol, int ndim, bool l_implemented,       // In
                                           Real3D lhs, Real2D rhs, err_info_type& err_info, // InOut
                                           Real2D soln,                                     // Out
                                           Real1D* rcond = nullptr); // Optional Out

void tridiag_solve_multiple_rhs_lhs(ScratchWorkspace& workspace, const std::string& solve_name,
                                    int tridiag_solve_method, int ngrdcol, int ndim, int nrhs,
                                    bool l_implemented,                              // In
                                    Real3D lhs, Real3D rhs, err_info_type& err_info, // InOut
                                    Real3D soln,                                     // Out
                                    Real1D* rcond = nullptr);                        // Optional Out

// C++ overloads implement the native generic interfaces.
inline void band_solve(ScratchWorkspace& workspace, const std::string& solve_name,
                       int penta_solve_method, int ngrdcol, int nsup, int nsub, int ndim,
                       bool l_implemented, Real3D lhs, Real2D rhs, err_info_type& err_info,
                       Real2D soln, const Real2D* old_soln = nullptr, Real1D* rcond = nullptr) {
  // Dispatch the native generic by matrix/RHS rank.
  auto workspace_scope = workspace.scope();
  band_solve_single_rhs_multiple_lhs(workspace, solve_name, penta_solve_method, // In
                                     ngrdcol, nsup, nsub, ndim, l_implemented,  // In
                                     lhs, rhs, err_info,                        // InOut
                                     soln,                                      // Out
                                     old_soln,                                  // Optional In
                                     rcond);                                    // Optional Out
}

inline void band_solve(ScratchWorkspace& workspace, const std::string& solve_name,
                       int penta_solve_method, int ngrdcol, int nsup, int nsub, int ndim, int nrhs,
                       bool l_implemented, Real3D lhs, Real3D rhs, err_info_type& err_info,
                       Real3D soln, const Real3D* old_soln = nullptr, Real1D* rcond = nullptr) {
  // Dispatch the native generic by matrix/RHS rank.
  auto workspace_scope = workspace.scope();
  band_solve_multiple_rhs_lhs(workspace, solve_name, penta_solve_method,      // In
                              ngrdcol, nsup, nsub, ndim, nrhs, l_implemented, // In
                              lhs, rhs, err_info,                             // InOut
                              soln,                                           // Out
                              old_soln,                                       // Optional In
                              rcond);                                         // Optional Out
}

inline void tridiag_solve(ScratchWorkspace& workspace, const std::string& solve_name,
                          int tridiag_solve_method, int ndim, Real2D lhs, Real1D rhs,
                          err_info_type& err_info, Real1D soln, Real1D* rcond = nullptr) {
  // Dispatch the native generic by matrix/RHS rank.
  auto workspace_scope = workspace.scope();
  tridiag_solve_single_rhs_lhs(workspace, solve_name, tridiag_solve_method, // In
                               ndim,                                        // In
                               lhs, rhs, err_info,                          // InOut
                               soln,                                        // Out
                               rcond);                                      // Optional Out
}

inline void tridiag_solve(ScratchWorkspace& workspace, const std::string& solve_name,
                          int tridiag_solve_method, int ngrdcol, int ndim, bool l_implemented,
                          Real3D lhs, Real2D rhs, err_info_type& err_info, Real2D soln,
                          Real1D* rcond = nullptr) {
  // Dispatch the native generic by matrix/RHS rank.
  auto workspace_scope = workspace.scope();
  tridiag_solve_single_rhs_multiple_lhs(workspace, solve_name, tridiag_solve_method, // In
                                        ngrdcol, ndim, l_implemented,                // In
                                        lhs, rhs, err_info,                          // InOut
                                        soln,                                        // Out
                                        rcond);                                      // Optional Out
}

inline void tridiag_solve(ScratchWorkspace& workspace, const std::string& solve_name,
                          int tridiag_solve_method, int ngrdcol, int ndim, int nrhs,
                          bool l_implemented, Real3D lhs, Real3D rhs, err_info_type& err_info,
                          Real3D soln, Real1D* rcond = nullptr) {
  // Dispatch the native generic by matrix/RHS rank.
  auto workspace_scope = workspace.scope();
  tridiag_solve_multiple_rhs_lhs(workspace, solve_name, tridiag_solve_method, // In
                                 ngrdcol, ndim, nrhs, l_implemented,          // In
                                 lhs, rhs, err_info,                          // InOut
                                 soln,                                        // Out
                                 rcond);                                      // Optional Out
}

} // namespace clubb_kokkos
