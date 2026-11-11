// Compute Native left-hand-side operator coefficients.
//
// Physical inputs and scratch use column-first LayoutLeft Views. The complete
// Out operator keeps Native (band,column,level) order, with three bands ordered
// super/main/sub. Keep any solver-layout transpose at the later solver boundary.
// Description:
// Module mean_adv computes the mean advection terms for all of the
// time-tendency (prognostic) equations in the CLUBB parameterization.  All of
// the mean advection terms are solved for completely implicitly, and therefore
// become part of the left-hand side of their respective equations.
//
// Function term_ma_zt_lhs handles the mean advection terms for the variables
// located at thermodynamic grid levels.  These variables are:  rtm, thlm, wp3,
// all hydrometeor species, and sclrm.
//
// Function term_ma_zm_lhs handles the mean advection terms for the variables
// located at momentum grid levels.  The variables are:  wprtp, wpthlp, wp2,
// rtp2, thlp2, rtpthlp, up2, vp2, wpsclrp, sclrprtp, sclrpthlp, and sclrp2.
#include "mean_adv.hpp"
#include "constants_clubb.hpp"
#include <algorithm>
namespace clubb_kokkos {
namespace {
constexpr int ndiags3 = 3;
} // namespace
//=============================================================================
// wm_zt: wm_zt                        [m/s]
// weights_zt2zm: Thermodynamic-to-momentum interpolation weights (column,momentum
// level,upper/lower category) [-]; Native categories1:2
// invrs_dzt: Inverse of grid spacing      [1/m]
// invrs_dzm: Inverse of grid spacing      [1/m]
// l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent or mean advection terms. It
// affects rtm, thlm, sclrm, um and vm.
// grid_dir: Grid direction multiplier (+1 = ascending; -1 = descending)
// lhs_ma: Mean advection contributions to lhs    [1/s]
void term_ma_zt_lhs(int nzm, int nzt, int ngrdcol, Real2D wm_zt, Real3D weights_zt2zm, // In
                    Real2D invrs_dzt, Real2D invrs_dzm,                                // In
                    bool l_upwind_xm_ma, real_t grid_dir,                              // In
                    Real3D lhs_ma)                                                     // Out
{
  // Description:
  // Mean advection of var_zt:  implicit portion of the code.
  //
  // The variable "var_zt" stands for a variable that is located at
  // thermodynamic grid levels.
  //
  // The d(var_zt)/dt equation contains a mean advection term:
  //
  // - w * d(var_zt)/dz.
  //
  // This term is solved for completely implicitly, such that:
  //
  // - w * d( var_zt(t+1) )/dz.
  //
  // Note:  When the term is brought over to the left-hand side, the sign
  //        is reversed and the leading "-" in front of the term is changed to
  //        a "+".
  //
  // The timestep index (t+1) means that the value of var_zt being used is from
  // the next timestep, which is being advanced to in solving the d(var_zt)/dt
  // equation.
  //
  // This term is discretized as follows:
  //
  // The values of var_zt are found on the thermodynamic levels, as are the
  // values of wm_zt (mean vertical velocity on thermodynamic levels).  The
  // variable var_zt is interpolated to the intermediate momentum levels.  The
  // derivative of the interpolated values is taken over the central
  // thermodynamic level.  The derivative is multiplied by wm_zt at the central
  // thermodynamic level to get the desired result.
  //
  // -----var_zt---------------------------------------------- t(k+1)
  //
  // =================var_zt(interp)========================== m(k+1)
  //
  // -----var_zt---------------------d(var_zt)/dz-----wm_zt--- t(k)
  //
  // =================var_zt(interp)========================== m(k)
  //
  // -----var_zt---------------------------------------------- t(k-1)
  //
  // The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
  // with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
  // The letter "t" is used for thermodynamic levels and the letter "m" is used
  // for momentum levels.
  //
  // invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) )
  //
  //
  // Special discretization for upper and lower boundary levels:
  //
  // Zero derivative method: the derivative d(var_zt)/dz is set to 0 over the
  // model top and the model bottom.
  //
  // This method corresponds with the "zero-flux" boundary condition option
  // for eddy diffusion, where d(var_zt)/dz is set to 0 across the upper
  // boundary.
  //
  // In order to discretize the upper boundary condition, consider a new level
  // outside the model (thermodynamic level gr%nzt+1) just above the upper
  // boundary level (thermodynamic level gr%nzt).  The value of var_zt at the
  // level just outside the model is defined to be the same as the value of
  // var_zt at thermodynamic level gr%nzt.  Therefore, the value of
  // d(var_zt)/dz between the level just outside the model and the uppermost
  // thermodynamic level is 0, staying consistent with the zero-flux boundary
  // condition option for the eddy diffusion portion of the code.  Therefore,
  // the value of var_zt at momentum level gr%nzm, which is the upper boundary
  // of the model, would be the same as the value of var_zt at the uppermost
  // thermodynamic level.
  //
  // The values of var_zt are found on the thermodynamic levels, as are the
  // values of wm_zt (mean vertical velocity on the thermodynamic levels).  The
  // variable var_zt is interpolated to momentum level gr%nzm-1, based on
  // the values of var_zt at thermodynamic levels gr%nzt and gr%nzt-1.  The
  // value of var_zt at momentum level gr%nzm is set equal to the value of
  // var_zt at thermodynamic level gr%nzt, as described above.  The derivative
  // of the set and interpolated values, d(var_zt)/dz, is taken over the
  // central thermodynamic level.  The derivative is multiplied by wm_zt at the
  // central thermodynamic level to get the desired result.
  //
  // For the following diagram, k = gr%nzm, which is the uppermost level of
  // the model:
  //
  // --[var_zt(k+1) = var_zt(k)]----(level outside model)----- t(k+1)
  //
  // ==[var_zt(top) = var_zt(k)]===[d(var_zt)/dz|_(top) = 0]== m(k+1) Boundary
  //
  // -----var_zt(k)------------------d(var_zt)/dz-----wm_zt--- t(k)
  //
  // =================var_zt(interp)========================== m(k)
  //
  // -----var_zt(k-1)----------------------------------------- t(k-1)
  //
  // where (top) stands for the grid index of momentum level k = gr%nzm, which
  // is the upper boundary of the model.
  //
  // An analogous discretization occurs at the lower boundary.
  // References:
  //   None
  //-----------------------------------------------------------------------
  // Constant parameters
  // -------------------------- Input Variables --------------------------
  // -------------------------- Return Variable --------------------------
  // -------------------------- Local Variables --------------------------
  constexpr int kp1_tdiag = 1; // Thermodynamic superdiagonal index.
  constexpr int k_tdiag   = 2; // Thermodynamic main diagonal index.
  constexpr int km1_tdiag = 3; // Thermodynamic subdiagonal index.
  constexpr int t_above   = 1; // Index for upper thermodynamic level grid weight.
  constexpr int t_below   = 2; // Index for lower thermodynamic level grid weight.
  int i;
  int k; // Vertical level index

  //-------------------------- Begin Code --------------------------

  // Use centered differencing
  if (!l_upwind_xm_ma) {

    // Most of the interior model; normal conditions.
    Kokkos::parallel_for(
        "term_ma_zt_lhs::loop1", std::max(0, ((nzt - 1 - (2)) / (1) + 1)) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 2 + remaining % (((nzt - 1 - (2)) / (1) + 1));
          remaining /= (((nzt - 1 - (2)) / (1) + 1));

          // Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
          lhs_ma(kp1_tdiag - 1, i - 1, k - 1) = +wm_zt(i - 1, k - 1) * invrs_dzt(i - 1, k - 1) *
                                                weights_zt2zm(i - 1, k + 1 - 1, t_above - 1);

          // Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs_ma(k_tdiag - 1, i - 1, k - 1) = +wm_zt(i - 1, k - 1) * invrs_dzt(i - 1, k - 1) *
                                              (weights_zt2zm(i - 1, k + 1 - 1, t_below - 1) -
                                               weights_zt2zm(i - 1, k - 1, t_above - 1));

          // Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
          lhs_ma(km1_tdiag - 1, i - 1, k - 1) = -wm_zt(i - 1, k - 1) * invrs_dzt(i - 1, k - 1) *
                                                weights_zt2zm(i - 1, k - 1, t_below - 1);
        });

    // Upper Boundary for an Ascending Grid
    // or Lower Boundary for a Descending Grid.

    // Special discretization for zero derivative method, where the
    // derivative d(var_zt)/dz over the model boundary is set to 0.
    Kokkos::parallel_for(
        "term_ma_zt_lhs::loop2", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);

          // Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
          lhs_ma(kp1_tdiag - 1, i - 1, nzt - 1) = zero;

          // Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs_ma(k_tdiag - 1, i - 1, nzt - 1) =
              +wm_zt(i - 1, nzt - 1) * invrs_dzt(i - 1, nzt - 1) *
              (one - weights_zt2zm(i - 1, nzm - 1 - 1, t_above - 1));

          // Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
          lhs_ma(km1_tdiag - 1, i - 1, nzt - 1) = -wm_zt(i - 1, nzt - 1) *
                                                  invrs_dzt(i - 1, nzt - 1) *
                                                  weights_zt2zm(i - 1, nzm - 1 - 1, t_below - 1);
        });

    // Lower Boundary for an Ascending Grid
    // or Upper Boundary for a Descending Grid.

    // Special discretization for zero derivative method, where the
    // derivative d(var_zt)/dz over the model boundary is set to 0.
    Kokkos::parallel_for(
        "term_ma_zt_lhs::loop3", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);

          // Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
          lhs_ma(kp1_tdiag - 1, i - 1, 1 - 1) = +wm_zt(i - 1, 1 - 1) * invrs_dzt(i - 1, 1 - 1) *
                                                weights_zt2zm(i - 1, 2 - 1, t_above - 1);

          // Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs_ma(k_tdiag - 1, i - 1, 1 - 1) = -wm_zt(i - 1, 1 - 1) * invrs_dzt(i - 1, 1 - 1) *
                                              (one - weights_zt2zm(i - 1, 2 - 1, t_below - 1));

          // Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
          lhs_ma(km1_tdiag - 1, i - 1, 1 - 1) = zero;
        });

  } else {
    // l_upwind_xm_ma == .true.; use "upwind" differencing

    // Most of the interior model; normal conditions.
    Kokkos::parallel_for(
        "term_ma_zt_lhs::loop4", std::max(0, ((nzt - 1 - (2)) / (1) + 1)) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 2 + remaining % (((nzt - 1 - (2)) / (1) + 1));
          remaining /= (((nzt - 1 - (2)) / (1) + 1));
          if (grid_dir * wm_zt(i - 1, k - 1) >= zero) {

            // Mean wind is in an upward direction and an ascending grid is used
            // or mean wind is in a downward direction and a descending grid
            // is used.

            // Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
            lhs_ma(kp1_tdiag - 1, i - 1, k - 1) = zero;

            // Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
            lhs_ma(k_tdiag - 1, i - 1, k - 1) = +wm_zt(i - 1, k - 1) * invrs_dzm(i - 1, k - 1);

            // Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
            lhs_ma(km1_tdiag - 1, i - 1, k - 1) = -wm_zt(i - 1, k - 1) * invrs_dzm(i - 1, k - 1);

          } else {
            // grid_dir * wm_zt < 0

            // Mean wind is in a downward direction and an ascending grid is
            // used or mean wind is in an upward direction and a descending grid
            // is used.

            // Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
            lhs_ma(kp1_tdiag - 1, i - 1, k - 1) =
                +wm_zt(i - 1, k - 1) * invrs_dzm(i - 1, k + 1 - 1);

            // Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
            lhs_ma(k_tdiag - 1, i - 1, k - 1) = -wm_zt(i - 1, k - 1) * invrs_dzm(i - 1, k + 1 - 1);

            // Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
            lhs_ma(km1_tdiag - 1, i - 1, k - 1) = zero;

            // wm_zt > 0
          }
        });

    // Upper Boundary for an Ascending Grid
    // or Lower Boundary for a Descending Grid.
    Kokkos::parallel_for(
        "term_ma_zt_lhs::loop5", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          if (grid_dir * wm_zt(i - 1, nzt - 1) >= zero) {

            // Mean wind is in an upward direction and an ascending grid is used or
            // mean wind is in a downward direction and a descending grid is used.

            // Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
            lhs_ma(kp1_tdiag - 1, i - 1, nzt - 1) = zero;

            // Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
            lhs_ma(k_tdiag - 1, i - 1, nzt - 1) =
                +wm_zt(i - 1, nzt - 1) * invrs_dzm(i - 1, nzm - 1 - 1);

            // Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
            lhs_ma(km1_tdiag - 1, i - 1, nzt - 1) =
                -wm_zt(i - 1, nzt - 1) * invrs_dzm(i - 1, nzm - 1 - 1);

          } else {
            // grid_dir * wm_zt < 0

            // Mean wind is in a downward direction and an ascending grid is used
            // or mean wind is in an upward direction and a descending grid is used

            // Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
            lhs_ma(kp1_tdiag - 1, i - 1, nzt - 1) = zero;

            // Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
            lhs_ma(k_tdiag - 1, i - 1, nzt - 1) = zero;

            // Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
            lhs_ma(km1_tdiag - 1, i - 1, nzt - 1) = zero;

            // grid_dir * wm_zt > 0
          }
        });

    // Lower Boundary for an Ascending Grid
    // or Upper Boundary for a Descending Grid.
    Kokkos::parallel_for(
        "term_ma_zt_lhs::loop6", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          if (grid_dir * wm_zt(i - 1, 1 - 1) >= zero) {

            // Mean wind is in an upward direction and an ascending grid is used or
            // mean wind is in a downward direction and a descending grid is used.

            // Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
            lhs_ma(kp1_tdiag - 1, i - 1, 1 - 1) = zero;

            // Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
            lhs_ma(k_tdiag - 1, i - 1, 1 - 1) = zero;

            // Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
            lhs_ma(km1_tdiag - 1, i - 1, 1 - 1) = zero;

          } else {
            // grid_dir * wm_zt < 0

            // Mean wind is in a downward direction and an ascending grid is used
            // or mean wind is in an upward direction and a descending grid is used

            // Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
            lhs_ma(kp1_tdiag - 1, i - 1, 1 - 1) = +wm_zt(i - 1, 1 - 1) * invrs_dzm(i - 1, 2 - 1);

            // Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
            lhs_ma(k_tdiag - 1, i - 1, 1 - 1) = -wm_zt(i - 1, 1 - 1) * invrs_dzm(i - 1, 2 - 1);

            // Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
            lhs_ma(km1_tdiag - 1, i - 1, 1 - 1) = zero;

            // grid_dir * wm_zt > 0
          }
        });

    // l_upwind_xm_ma
  }

  return;
}

//=============================================================================
// wm_zm: wm_zm                        [m/s]
// invrs_dzm: Inverse of grid spacing      [1/m]
// weights_zm2zt: Momentum-to-thermodynamic interpolation weights (column,thermodynamic
// level,upper/lower category) [-]; Native categories1:2
// lhs_ma: Mean advection contributions to lhs  [1/s]
void term_ma_zm_lhs(int nzm, int nzt, int ngrdcol, Real2D wm_zm, // In
                    Real2D invrs_dzm, Real3D weights_zm2zt,      // In
                    Real3D lhs_ma)                               // Out
{
  // Description:
  // Mean advection of var_zm:  implicit portion of the code.
  //
  // The variable "var_zm" stands for a variable that is located at momentum
  // grid levels.
  //
  // The d(var_zm)/dt equation contains a mean advection term:
  //
  // - w * d(var_zm)/dz.
  //
  // This term is solved for completely implicitly, such that:
  //
  // - w * d( var_zm(t+1) )/dz.
  //
  // Note:  When the term is brought over to the left-hand side, the sign
  //        is reversed and the leading "-" in front of the term is changed to
  //        a "+".
  //
  // The timestep index (t+1) means that the value of var_zm being used is from
  // the next timestep, which is being advanced to in solving the d(var_zm)/dt
  // equation.
  //
  // This term is discretized as follows:
  //
  // The values of var_zm are found on the momentum levels, as are the values
  // of wm_zm (mean vertical velocity on momentum levels).  The variable var_zm
  // is interpolated to the intermediate thermodynamic levels.  The derivative
  // of the interpolated values is taken over the central momentum level.  The
  // derivative is multiplied by wm_zm at the central momentum level to get the
  // desired result.
  //
  // =====var_zm============================================== m(k+1)
  //
  // -----------------var_zm(interp)-------------------------- t(k)
  //
  // =====var_zm=====================d(var_zm)/dz=====wm_zm=== m(k)
  //
  // -----------------var_zm(interp)-------------------------- t(k-1)
  //
  // =====var_zm============================================== m(k-1)
  //
  // The vertical indices m(k+1), t(k), m(k), t(k-1), and m(k-1) correspond
  // with altitudes zm(k+1), zt(k), zm(k), zt(k-1), and zm(k-1), respectively.
  // The letter "t" is used for thermodynamic levels and the letter "m" is used
  // for momentum levels.
  //
  // invrs_dzm(k) = 1 / ( zt(k) - zt(k-1) )
  // References:
  //-----------------------------------------------------------------------
  // Constant parameters
  // -------------------------- Input Variables --------------------------
  // -------------------------- Return Variable --------------------------
  // -------------------------- Local Variables
  constexpr int kp1_mdiag = 1; // Momentum superdiagonal index.
  constexpr int k_mdiag   = 2; // Momentum main diagonal index.
  constexpr int km1_mdiag = 3; // Momentum subdiagonal index.
  constexpr int m_above   = 1; // Index for upper momentum level grid weight.
  constexpr int m_below   = 2; // Index for lower momentum level grid weight.
  int i;
  int k;
  int b; // Vertical level index

  // -------------------------- Begin Code --------------------------

  // Set boundary array to 0
  // This is the Lower Boundary for an Ascending Grid
  // or the Upper Boundary for a Descending Grid.
  Kokkos::parallel_for(
      "term_ma_zm_lhs::loop1", std::max(0, ngrdcol) * std::max(0, ndiags3),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int b   = 1 + remaining % (ndiags3);
        remaining /= (ndiags3);
        const int i = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        lhs_ma(b - 1, i - 1, 1 - 1) = zero;
      });

  // Most of the interior model; normal conditions.
  Kokkos::parallel_for(
      "term_ma_zm_lhs::loop2", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);

        // Momentum superdiagonal: [ x var_zm(k+1,<t+1>) ]
        lhs_ma(kp1_mdiag - 1, i - 1, k - 1) = +wm_zm(i - 1, k - 1) * invrs_dzm(i - 1, k - 1) *
                                              weights_zm2zt(i - 1, k - 1, m_above - 1);

        // Momentum main diagonal: [ x var_zm(k,<t+1>) ]
        lhs_ma(k_mdiag - 1, i - 1, k - 1) = +wm_zm(i - 1, k - 1) * invrs_dzm(i - 1, k - 1) *
                                            (weights_zm2zt(i - 1, k - 1, m_below - 1) -
                                             weights_zm2zt(i - 1, k - 1 - 1, m_above - 1));

        // Momentum subdiagonal: [ x var_zm(k-1,<t+1>) ]
        lhs_ma(km1_mdiag - 1, i - 1, k - 1) = -wm_zm(i - 1, k - 1) * invrs_dzm(i - 1, k - 1) *
                                              weights_zm2zt(i - 1, k - 1 - 1, m_below - 1);
      });

  // Set boundary array to 0
  // This is the Upper Boundary for an Ascending Grid
  // or the Lower Boundary for a Descending Grid.
  Kokkos::parallel_for(
      "term_ma_zm_lhs::loop3", std::max(0, ngrdcol) * std::max(0, ndiags3),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int b   = 1 + remaining % (ndiags3);
        remaining /= (ndiags3);
        const int i = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        lhs_ma(b - 1, i - 1, nzm - 1) = zero;
      });

  return;
}

} // namespace clubb_kokkos
