// Native pentadiagonal LU decompositions and substitutions.
// See the header for C++ array layout. The derivation below uses native indices.
// Description:
//   These routines solve lhs*soln=rhs using LU decomp.
//
//   LHS is stored in band diagonal form.
//     lhs = | lhs(0,1)  lhs(-1,1)  lhs(-2,1)     0           0          0          0
//           | lhs(1,2)  lhs( 0,2)  lhs(-1,2)  lhs(-2,2)      0          0          0
//           | lhs(2,3)  lhs( 1,3)  lhs( 0,3)  lhs(-1,3)  lhs(-2,3)      0          0
//           |     0     lhs( 2,4)  lhs( 1,4)  lhs( 0,4)  lhs(-1,4)  lhs(-2,4)      0
//           |     0         0      lhs( 2,5)  lhs( 1,5)  lhs( 0,5)  lhs(-1,5)  lhs(-2,5) ...
//           |    ...
//
//    U is stored in band diagonal form
//     U = |   1    upper_1(1)  upper_2(1)      0           0           0           0
//         |   0        1       upper_1(2)  upper_2(2)      0           0           0
//         |   0        0           1       upper_1(3)  upper_2(3)      0           0
//         |   0        0           0           1       upper_1(4)  upper_2(4)      0
//         |   0        0           0           0           1       upper_1(5)  upper_2(5)  ...
//         |  ...
//
//   L is also stored in band diagonal form, but the lowest most band is equivalent to the
//   lowermost band of LHS, thus we don't need to store it
//     L = | l_diag(1)        0            0           0             0       0
//         | lower_1(2)    l_diag(2)       0           0             0       0
//         |   l_2(3)      lower_1(3)   l_diag(3)      0             0       0
//         |     0           l_2(4)     lower_1(4)  l_diag(4)        0       0
//         |     0             0        l_2(5)      lower_1(5)    l_diag(5)  0   ...
//         |  ...
//
//
//   To perform the LU decomposition, we go element by element.
//   First we start by noting that we want lhs=LU, so the first step of calculating
//   L*U, by multiplying the first row of L by the columns of U, gives us
//
//     l_diag(1)*1          = lhs( 0,1)  =>  l_diag(1)  = lhs( 0,1)
//     l_diag(1)*upper_1(1) = lhs(-1,1)  =>  upper_1(1) = lhs(-1,1) / l_diag(1)
//     l_diag(1)*upper_2(1) = lhs(-2,1)  =>  upper_2(1) = lhs(-2,1) / l_diag(1)
//
//   Multiplying the second row of L by U now we get
//
//     lower_1(2)*1                               = lhs(1,2)   =>  lower_1(2)  = lhs(1,2)
//     lower_1(2)*upper_1(1)+l_diag(2)*1          = lhs(0,2)   =>  l_diag(2)
//                                                                 = lhs(0,2)
//                                                                   - lower_1(2)*upper_1(1)
//     lower_1(2)*upper_2(1)+l_diag(2)*upper_1(2) = lhs(-1,2)  =>  upper_1(2)
//                                                                 = ( lhs(-1,2)
//                                                                   - lower_1(2)*upper_2(1) )
//                                                                   / l_diag(2)
//     l_diag(2)*upper_2(2)                       = lhs(-2,2)  =>  upper_2(2)  =
//                                                                             lhs(-2,2)
//                                                                             / l_diag(2)
//
//   Now that we're passed the k=1 and k=2 steps, each following step uses all the bands,
//   allowing us to write the general step
//
//     l_2(k)*1                                    = lhs(2,k)   =>  l_2(k)     = lhs(2,k)
//     l_2(k)*upper_1(k-2)+lower_1(k)*1            = lhs(1,k)   =>  lower_1(k)
//                                                                  = lhs(1,k)
//                                                                    - l_2(k)*upper_1(k-2)
//     l_2(k)*upper_2(k-2)+lower_1(k)*upper_1(k-1) = lhs( 0,k)  =>  l_diag(k)  = lhs(0,k)
//                    +l_diag(k)*1                                               - l_2(k)
//                                                                               *upper_2(k-2)
//                                                                               + lower_1(k)
//                                                                               *upper_1(k-1)
//
//     lower_1(k)*upper_2(k-1)
//     + l_diag(k)*upper_1(k) = lhs(-1,k)  =>  upper_1(k) = ( lhs(-1,k) - lower_1(k)*upper_2(k-1) )
//                                                          / l_diag(k)
//     l_diag(k)*upper_2(k)   = lhs(-2,k)  =>  upper_2(k) = lhs(-2,k) / l_diag(k)
//
//
//   This general step is done for k from 3 to ndim-2 (do k = 3, ndim-2), and the last two
//   steps are tweaked similarly to the first two, where we disclude one then two bands
//   since they become no longer relevant. Note from this general step that the l_2 band
//   is always equivalent to second subdiagonal band of lhs, thus we do not need to
//   calculate or store l_2. Also note that we only ever need l_diag so that we can divide
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
#include "penta_lu_solver.hpp"

namespace clubb_kokkos {

void penta_lu_solve_single_rhs_multiple_lhs(ScratchWorkspace& workspace, // InOut
                                            int ndim, int ngrdcol,       // In
                                            Real3D lhs,                  // InOut
                                            Real2D rhs,                  // In
                                            Real2D soln) {               // Out
  // Written for single RHS and multiple LHS. One work item per column performs
  // the ordered decomposition and substitutions. Independent native column
  // sweeps, including boundary rows, share one launch.

  // ----------------------- Local Variables -----------------------
  auto workspace_scope = workspace.scope();
  Real2D upper_1 = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                "penta_lu_solve_single_rhs_multiple_lhs::upper_1",
                                                ngrdcol, ndim); // First U band [-]
  Real2D upper_2 = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                "penta_lu_solve_single_rhs_multiple_lhs::upper_2",
                                                ngrdcol, ndim); // Second U band [-]
  Real2D lower_diag_invrs = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "penta_lu_solve_single_rhs_multiple_lhs::lower_diag_invrs",
      ngrdcol, ndim); // Inverse L diagonal [inverse lhs units]
  Real2D lower_1 = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                "penta_lu_solve_single_rhs_multiple_lhs::lower_1",
                                                ngrdcol, ndim); // First L band [lhs units]
  Real2D lower_2 = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                "penta_lu_solve_single_rhs_multiple_lhs::lower_2",
                                                ngrdcol, ndim); // Second L band [lhs units]

  // ----------------------- Begin Code -----------------------
  Kokkos::parallel_for(
      "penta_lu_solve_single_rhs_multiple_lhs::solve", Kokkos::RangePolicy<>(0, ngrdcol),
      KOKKOS_LAMBDA(const int i) {
        lower_diag_invrs(i, 0) = real_t(1.0) / lhs(i, 2, 0);
        upper_1(i, 0)          = lower_diag_invrs(i, 0) * lhs(i, 1, 0);
        upper_2(i, 0)          = lower_diag_invrs(i, 0) * lhs(i, 0, 0);

        lower_1(i, 1)          = lhs(i, 3, 1);
        lower_diag_invrs(i, 1) = real_t(1.0) / (lhs(i, 2, 1) - lower_1(i, 1) * upper_1(i, 0));
        upper_1(i, 1) = lower_diag_invrs(i, 1) * (lhs(i, 1, 1) - lower_1(i, 1) * upper_2(i, 0));
        upper_2(i, 1) = lower_diag_invrs(i, 1) * lhs(i, 0, 1);

        for (int k = 2; k < ndim - 2; ++k) {
          lower_2(i, k)          = lhs(i, 4, k);
          lower_1(i, k)          = lhs(i, 3, k) - lower_2(i, k) * upper_1(i, k - 2);
          lower_diag_invrs(i, k) = real_t(1.0) / (lhs(i, 2, k) - lower_2(i, k) * upper_2(i, k - 2) -
                                                  lower_1(i, k) * upper_1(i, k - 1));
          upper_1(i, k) =
              lower_diag_invrs(i, k) * (lhs(i, 1, k) - lower_1(i, k) * upper_2(i, k - 1));
          upper_2(i, k) = lower_diag_invrs(i, k) * lhs(i, 0, k);
        }

        lower_2(i, ndim - 2) = lhs(i, 4, ndim - 2);
        lower_1(i, ndim - 2) = lhs(i, 3, ndim - 2) - lower_2(i, ndim - 2) * upper_1(i, ndim - 4);
        lower_diag_invrs(i, ndim - 2) =
            real_t(1.0) / (lhs(i, 2, ndim - 2) - lower_2(i, ndim - 2) * upper_2(i, ndim - 4) -
                           lower_1(i, ndim - 2) * upper_1(i, ndim - 3));
        upper_1(i, ndim - 2) = lower_diag_invrs(i, ndim - 2) *
                               (lhs(i, 1, ndim - 2) - lower_1(i, ndim - 2) * upper_2(i, ndim - 3));

        lower_2(i, ndim - 1) = lhs(i, 4, ndim - 1);
        lower_1(i, ndim - 1) = lhs(i, 3, ndim - 1) - lower_2(i, ndim - 1) * upper_1(i, ndim - 3);
        lower_diag_invrs(i, ndim - 1) =
            real_t(1.0) / (lhs(i, 2, ndim - 1) - lower_2(i, ndim - 1) * upper_2(i, ndim - 3) -
                           lower_1(i, ndim - 1) * upper_1(i, ndim - 2));

        soln(i, 0) = lower_diag_invrs(i, 0) * rhs(i, 0);
        soln(i, 1) = lower_diag_invrs(i, 1) * (rhs(i, 1) - lower_1(i, 1) * soln(i, 0));

        for (int k = 2; k < ndim; ++k) {
          soln(i, k) = lower_diag_invrs(i, k) * (rhs(i, k) - lower_2(i, k) * soln(i, k - 2) -
                                                 lower_1(i, k) * soln(i, k - 1));
        }

        soln(i, ndim - 2) = soln(i, ndim - 2) - upper_1(i, ndim - 2) * soln(i, ndim - 1);

        for (int k = ndim - 3; k >= 0; --k) {
          soln(i, k) = soln(i, k) - upper_1(i, k) * soln(i, k + 1) - upper_2(i, k) * soln(i, k + 2);
        }
      });
}

void penta_lu_solve_multiple_rhs_lhs(ScratchWorkspace& workspace,     // InOut
                                     int ndim, int nrhs, int ngrdcol, // In
                                     Real3D lhs,                      // InOut
                                     Real3D rhs,                      // In
                                     Real3D soln) {                   // Out
  // Written for multiple RHS and multiple LHS. Decompose once per column,
  // then solve each RHS in its own work item. Vertical sweeps retain native
  // arithmetic order; independent boundary/interior stages share each launch.

  // ----------------------- Local Variables -----------------------
  auto workspace_scope = workspace.scope();
  Real2D upper_1 = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                "penta_lu_solve_multiple_rhs_lhs::upper_1", ngrdcol,
                                                ndim); // First U band [-]
  Real2D upper_2 = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                "penta_lu_solve_multiple_rhs_lhs::upper_2", ngrdcol,
                                                ndim); // Second U band [-]
  Real2D lower_diag_invrs = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "penta_lu_solve_multiple_rhs_lhs::lower_diag_invrs", ngrdcol,
      ndim); // Inverse L diagonal [inverse lhs units]
  Real2D lower_1 = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                "penta_lu_solve_multiple_rhs_lhs::lower_1", ngrdcol,
                                                ndim); // First L band [lhs units]
  Real2D lower_2 = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                "penta_lu_solve_multiple_rhs_lhs::lower_2", ngrdcol,
                                                ndim); // Second L band [lhs units]

  // ----------------------- Begin Code -----------------------
  Kokkos::parallel_for(
      "penta_lu_solve_multiple_rhs_lhs::decompose", Kokkos::RangePolicy<>(0, ngrdcol),
      KOKKOS_LAMBDA(const int i) {
        lower_diag_invrs(i, 0) = real_t(1.0) / lhs(i, 2, 0);
        upper_1(i, 0)          = lower_diag_invrs(i, 0) * lhs(i, 1, 0);
        upper_2(i, 0)          = lower_diag_invrs(i, 0) * lhs(i, 0, 0);
        lower_1(i, 1)          = lhs(i, 3, 1);
        lower_diag_invrs(i, 1) = real_t(1.0) / (lhs(i, 2, 1) - lower_1(i, 1) * upper_1(i, 0));
        upper_1(i, 1) = lower_diag_invrs(i, 1) * (lhs(i, 1, 1) - lower_1(i, 1) * upper_2(i, 0));
        upper_2(i, 1) = lower_diag_invrs(i, 1) * lhs(i, 0, 1);

        for (int k = 2; k < ndim - 2; ++k) {
          lower_2(i, k)          = lhs(i, 4, k);
          lower_1(i, k)          = lhs(i, 3, k) - lower_2(i, k) * upper_1(i, k - 2);
          lower_diag_invrs(i, k) = real_t(1.0) / (lhs(i, 2, k) - lower_2(i, k) * upper_2(i, k - 2) -
                                                  lower_1(i, k) * upper_1(i, k - 1));
          upper_1(i, k) =
              lower_diag_invrs(i, k) * (lhs(i, 1, k) - lower_1(i, k) * upper_2(i, k - 1));
          upper_2(i, k) = lower_diag_invrs(i, k) * lhs(i, 0, k);
        }

        lower_2(i, ndim - 2) = lhs(i, 4, ndim - 2);
        lower_1(i, ndim - 2) = lhs(i, 3, ndim - 2) - lower_2(i, ndim - 2) * upper_1(i, ndim - 4);
        lower_diag_invrs(i, ndim - 2) =
            real_t(1.0) / (lhs(i, 2, ndim - 2) - lower_2(i, ndim - 2) * upper_2(i, ndim - 4) -
                           lower_1(i, ndim - 2) * upper_1(i, ndim - 3));
        upper_1(i, ndim - 2) = lower_diag_invrs(i, ndim - 2) *
                               (lhs(i, 1, ndim - 2) - lower_1(i, ndim - 2) * upper_2(i, ndim - 3));

        lower_2(i, ndim - 1) = lhs(i, 4, ndim - 1);
        lower_1(i, ndim - 1) = lhs(i, 3, ndim - 1) - lower_2(i, ndim - 1) * upper_1(i, ndim - 3);
        lower_diag_invrs(i, ndim - 1) =
            real_t(1.0) / (lhs(i, 2, ndim - 1) - lower_2(i, ndim - 1) * upper_2(i, ndim - 3) -
                           lower_1(i, ndim - 1) * upper_1(i, ndim - 2));
      });

  Kokkos::parallel_for(
      "penta_lu_solve_multiple_rhs_lhs::solve", Kokkos::RangePolicy<>(0, ngrdcol * nrhs),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int j = cell / ngrdcol;

        soln(i, 0, j) = lower_diag_invrs(i, 0) * rhs(i, 0, j);
        soln(i, 1, j) = lower_diag_invrs(i, 1) * (rhs(i, 1, j) - lower_1(i, 1) * soln(i, 0, j));

        for (int k = 2; k < ndim; ++k) {
          soln(i, k, j) =
              lower_diag_invrs(i, k) * (rhs(i, k, j) - lower_2(i, k) * soln(i, k - 2, j) -
                                        lower_1(i, k) * soln(i, k - 1, j));
        }

        soln(i, ndim - 2, j) = soln(i, ndim - 2, j) - upper_1(i, ndim - 2) * soln(i, ndim - 1, j);

        for (int k = ndim - 3; k >= 0; --k) {
          soln(i, k, j) =
              soln(i, k, j) - upper_1(i, k) * soln(i, k + 1, j) - upper_2(i, k) * soln(i, k + 2, j);
        }
      });
}

} // namespace clubb_kokkos
