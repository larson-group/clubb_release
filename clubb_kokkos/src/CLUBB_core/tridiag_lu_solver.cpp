// Native tridiagonal LU decompositions and substitutions.
// See the header for C++ array layout. The derivation below uses native indices.
// Description:
//   These routines solve lhs*soln=rhs using LU decomp.
//
//   LHS is stored in band diagonal form.
//     lhs = | lhs(0,1)  lhs(-1,1)      0          0          0          0      0
//           | lhs(1,2)  lhs( 0,2)  lhs(-1,2)      0          0          0      0
//           |     0     lhs( 1,3)  lhs( 0,3)  lhs(-1,3)      0          0      0
//           |     0         0      lhs( 1,4)  lhs( 0,4)  lhs(-1,4)      0      0
//           |     0         0          0      lhs( 1,5)  lhs( 0,5)  lhs(-1,5)  0 ...
//           |    ...
//
//    U is stored in band diagonal form
//     U = |   1    upper(1)  0       0       0       0       0
//         |   0      1     upper(2)  0       0       0       0
//         |   0      0       1     upper(3)  0       0       0
//         |   0      0       0       1     upper(4)  0       0
//         |   0      0       0       0       1     upper(5)  0  ...
//         |  ...
//
//   L is also stored in band diagonal form, but the lowest most band is equivalent to the
//   lowermost band of LHS, thus we don't need to store it
//     L = | l_diag(1)       0           0            0          0        0
//         |  lower(2)      l_diag(2)    0            0          0        0
//         |    0           lower(3)   l_diag(3)      0          0        0
//         |    0            0         lower(4)    l_diag(4)     0        0
//         |    0            0           0         lower(5)    l_diag(5)   0   ...
//         |  ...
//
//
//   To perform the LU decomposition, we go element by element.
//   First we start by noting that we want lhs=LU, so the first step of calculating
//   L*U, by multiplying the first row of L by the columns of U, gives us
//
//     l_diag(1)*1        = lhs( 0,1)  =>  l_diag(1)   = lhs( 0,1)
//     l_diag(1)*upper(1) = lhs(-1,1)  =>  upper(1)    = lhs(-1,1) / l_diag(1)
//
//   Now that we're passed the k=1 step, each following step uses all the bands,
//   allowing us to write the general step
//
//     lower(k)*1 = lhs(1,k)                         =>  lower(k)    = lhs( 1,k)
//     lower(k)*upper(k-1)+l_diag(k)*1  = lhs( 0,k)  =>  l_diag(k) = lhs( 0,k)
//                                                                   - lower(k)*upper(k-1)
//     l_diag(k)*upper(k)               = lhs(-1,k)  =>  upper(k)  = lhs(-1,k) / l_diag(k)
//
//   This general step is done for k from 2 to ndim-1 (do k = 2, ndim-), and the last
//   step is tweaked similarly to the first one, where we disclude the upper band
//   since it becomes no longer relevant. Note from this general step that the lower band
//   is always equivalent to first subdiagonal band of lhs, thus we do not need to
//   calculate or store lower. Also note that we only ever need l_diag so that we can divide
//   by it, so instead we compute lower_diag_invrs to reduce divide operations.
//
//   After L and U are computed, normally we do forward substitution using L,
//   then backward substitution using U to find the solution. This is replicated
//   for every right hand side we want to solve for.
//
//
// References:
//   none
//------------------------------------------------------------------------

#include "scratch_workspace.hpp"
#include "tridiag_lu_solver.hpp"

namespace clubb_kokkos {

void tridiag_lu_solve_single_rhs_lhs(ScratchWorkspace& workspace, // InOut
                                     int ndim,                    // In
                                     Real2D lhs,                  // InOut
                                     Real1D rhs,                  // In
                                     Real1D soln) {               // Out
  // Written for single RHS and single LHS. One GPU work item performs the
  // ordered decomposition and both substitutions without interlevel parallelism.

  // ----------------------- Local Variables -----------------------
  auto workspace_scope    = workspace.scope();
  Real1D upper            = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                         "tridiag_lu_solve_single_rhs_lhs::upper",
                                                         ndim); // First U band [-]
  Real1D lower_diag_invrs = workspace_scope.view<Real1D>(
      Kokkos::WithoutInitializing, "tridiag_lu_solve_single_rhs_lhs::lower_diag_invrs",
      ndim); // Inverse L diagonal [inverse lhs units]

  // ----------------------- Begin Code -----------------------
  Kokkos::parallel_for(
      "tridiag_lu_solve_single_rhs_lhs::solve", Kokkos::RangePolicy<>(0, 1),
      KOKKOS_LAMBDA(const int) {
        lower_diag_invrs(0) = real_t(1.0) / lhs(1, 0);
        upper(0)            = lower_diag_invrs(0) * lhs(0, 0);

        for (int k = 1; k < ndim - 1; ++k) {
          lower_diag_invrs(k) = real_t(1.0) / (lhs(1, k) - lhs(2, k) * upper(k - 1));
          upper(k)            = lower_diag_invrs(k) * lhs(0, k);
        }

        lower_diag_invrs(ndim - 1) =
            real_t(1.0) / (lhs(1, ndim - 1) - lhs(2, ndim - 1) * upper(ndim - 2));

        soln(0) = lower_diag_invrs(0) * rhs(0);

        for (int k = 1; k < ndim; ++k) {
          soln(k) = lower_diag_invrs(k) * (rhs(k) - lhs(2, k) * soln(k - 1));
        }

        for (int k = ndim - 2; k >= 0; --k) {
          soln(k) = soln(k) - upper(k) * soln(k + 1);
        }
      });
}

void tridiag_lu_solve_single_rhs_multiple_lhs(ScratchWorkspace& workspace, // InOut
                                              int ndim, int ngrdcol,       // In
                                              Real3D lhs,                  // InOut
                                              Real2D rhs,                  // In
                                              Real2D soln) {               // Out
  // Written for single RHS and multiple LHS. One work item per column performs
  // the ordered decomposition and substitutions. Independent native column
  // sweeps, including boundary rows, share one launch.

  // ----------------------- Local Variables -----------------------
  auto workspace_scope = workspace.scope();
  Real2D upper = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                              "tridiag_lu_solve_single_rhs_multiple_lhs::upper",
                                              ngrdcol, ndim); // First U band [-]
  Real2D lower_diag_invrs = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "tridiag_lu_solve_single_rhs_multiple_lhs::lower_diag_invrs",
      ngrdcol, ndim); // Inverse L diagonal [inverse lhs units]

  // ----------------------- Begin Code -----------------------
  Kokkos::parallel_for(
      "tridiag_lu_solve_single_rhs_multiple_lhs::solve", Kokkos::RangePolicy<>(0, ngrdcol),
      KOKKOS_LAMBDA(const int i) {
        lower_diag_invrs(i, 0) = real_t(1.0) / lhs(i, 1, 0);
        upper(i, 0)            = lower_diag_invrs(i, 0) * lhs(i, 0, 0);

        for (int k = 1; k < ndim - 1; ++k) {
          lower_diag_invrs(i, k) = real_t(1.0) / (lhs(i, 1, k) - lhs(i, 2, k) * upper(i, k - 1));
          upper(i, k)            = lower_diag_invrs(i, k) * lhs(i, 0, k);
        }

        lower_diag_invrs(i, ndim - 1) =
            real_t(1.0) / (lhs(i, 1, ndim - 1) - lhs(i, 2, ndim - 1) * upper(i, ndim - 2));

        soln(i, 0) = lower_diag_invrs(i, 0) * rhs(i, 0);

        for (int k = 1; k < ndim; ++k) {
          soln(i, k) = lower_diag_invrs(i, k) * (rhs(i, k) - lhs(i, 2, k) * soln(i, k - 1));
        }

        for (int k = ndim - 2; k >= 0; --k) {
          soln(i, k) = soln(i, k) - upper(i, k) * soln(i, k + 1);
        }
      });
}

void tridiag_lu_solve_multiple_rhs_lhs(ScratchWorkspace& workspace,     // InOut
                                       int ndim, int nrhs, int ngrdcol, // In
                                       Real3D lhs,                      // InOut
                                       Real3D rhs,                      // In
                                       Real3D soln) {                   // Out
  // Written for multiple RHS and multiple LHS. Decompose once per column,
  // then solve each RHS in its own work item. Vertical sweeps retain native
  // arithmetic order; independent boundary/interior stages share each launch.

  // ----------------------- Local Variables -----------------------
  auto workspace_scope = workspace.scope();
  Real2D upper = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                              "tridiag_lu_solve_multiple_rhs_lhs::upper", ngrdcol,
                                              ndim); // First U band [-]
  Real2D lower_diag_invrs = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "tridiag_lu_solve_multiple_rhs_lhs::lower_diag_invrs", ngrdcol,
      ndim); // Inverse L diagonal [inverse lhs units]

  // ----------------------- Begin Code -----------------------
  Kokkos::parallel_for(
      "tridiag_lu_solve_multiple_rhs_lhs::decompose", Kokkos::RangePolicy<>(0, ngrdcol),
      KOKKOS_LAMBDA(const int i) {
        lower_diag_invrs(i, 0) = real_t(1.0) / lhs(i, 1, 0);
        upper(i, 0)            = lower_diag_invrs(i, 0) * lhs(i, 0, 0);

        for (int k = 1; k < ndim - 1; ++k) {
          lower_diag_invrs(i, k) = real_t(1.0) / (lhs(i, 1, k) - lhs(i, 2, k) * upper(i, k - 1));
          upper(i, k)            = lower_diag_invrs(i, k) * lhs(i, 0, k);
        }

        lower_diag_invrs(i, ndim - 1) =
            real_t(1.0) / (lhs(i, 1, ndim - 1) - lhs(i, 2, ndim - 1) * upper(i, ndim - 2));
      });

  Kokkos::parallel_for(
      "tridiag_lu_solve_multiple_rhs_lhs::solve", Kokkos::RangePolicy<>(0, ngrdcol * nrhs),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int j = cell / ngrdcol;

        soln(i, 0, j) = lower_diag_invrs(i, 0) * rhs(i, 0, j);

        for (int k = 1; k < ndim; ++k) {
          soln(i, k, j) =
              lower_diag_invrs(i, k) * (rhs(i, k, j) - lhs(i, 2, k) * soln(i, k - 1, j));
        }

        for (int k = ndim - 2; k >= 0; --k) {
          soln(i, k, j) = soln(i, k, j) - upper(i, k) * soln(i, k + 1, j);
        }
      });
}

} // namespace clubb_kokkos
