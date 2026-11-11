// Compute Native left-hand-side operator coefficients.
//
// Physical inputs and scratch use column-first LayoutLeft Views. The complete
// Out operator keeps Native (band,column,level) order, with three bands ordered
// super/main/sub. Keep any solver-layout transpose at the later solver boundary.
// Description:
// Module diffusion computes the eddy diffusion terms for all of the
// time-tendency (prognostic) equations in the CLUBB parameterization.  Most of
// the eddy diffusion terms are solved for completely implicitly, and therefore
// become part of the left-hand side of their respective equations.  However,
// wp2 and wp3 have an option to use a Crank-Nicholson eddy diffusion scheme,
// which has both implicit and explicit components.
//
// Function diffusion_zt_lhs handles the eddy diffusion terms for the variables
// located at thermodynamic grid levels.  These variables are:  wp3 and all
// hydrometeor species.  The variables um and vm also use the Crank-Nicholson
// eddy-diffusion scheme for their turbulent advection term.
//
// Function diffusion_zm_lhs handles the eddy diffusion terms for the variables
// located at momentum grid levels.  The variables are:  wprtp, wpthlp, wp2,
// rtp2, thlp2, rtpthlp, up2, vp2, wpsclrp, sclrprtp, sclrpthlp, and sclrp2.
#include "diffusion.hpp"
#include "constants_clubb.hpp"
#include "model_flags.hpp" // Selected l_upwind_Kh_dp_term owner remains pending.
#include <algorithm>
namespace clubb_kokkos {
//=============================================================================
// K_zm: Coef. of eddy diffusivity at momentum levels   [m^2/s]
// K_zt: Coef. of eddy diffusivity at thermo. levels    [m^2/s]
// nu: Background constant coef. of eddy diffusivity  [m^2/s]
// invrs_rho_ds_ztzxt: Inverse dry statis density on thermo. levels   [m^3/kg]
// rho_ds_zm: Dry statis density on momentum levels          [kg/m^3]
// lhs: LHS coefficient of diffusion term    [1/s]
void diffusion_zt_lhs(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D K_zm, Real2D K_zt, // In
                      Real1D nu,                                                               // In
                      Real2D invrs_rho_ds_ztzxt, Real2D rho_ds_zm,                             // In
                      Real3D lhs) // Out
{
  // Description:
  // Vertical eddy diffusion of var_zt:  implicit portion of the code.
  //
  // The variable "var_zt" stands for a variable that is located at
  // thermodynamic grid levels.
  //
  // The d(var_zt)/dt equation contains an eddy diffusion term:
  //
  // + d [ ( K_zm + nu ) * d(var_zt)/dz ] / dz.
  //
  // This term is usually solved for completely implicitly, such that:
  //
  // + d [ ( K_zm + nu ) * d( var_zt(t+1) )/dz ] / dz.
  //
  // However, when a Crank-Nicholson scheme is used, the eddy diffusion term
  // has both implicit and explicit components, such that:
  //
  // + (1/2) * d [ ( K_zm + nu ) * d( var_zt(t+1) )/dz ] / dz
  //    + (1/2) * d [ ( K_zm + nu ) * d( var_zt(t) )/dz ] / dz;
  //
  // for which the implicit component is:
  //
  // + (1/2) * d [ ( K_zm + nu ) * d( var_zt(t+1) )/dz ] / dz.
  //
  // Note:  When the implicit term is brought over to the left-hand side,
  //        the  sign is reversed and the leading "+" in front of the term
  //        is changed to a "-".
  //
  // Timestep index (t) stands for the index of the current timestep, while
  // timestep index (t+1) stands for the index of the next timestep, which is
  // being advanced to in solving the d(var_zt)/dt equation.
  //
  // The implicit portion of this term is discretized as follows:
  //
  // The values of var_zt are found on the thermodynamic levels, while the
  // values of K_zm are found on the momentum levels.  The derivatives (d/dz)
  // of var_zt are taken over the intermediate momentum levels.  At the
  // intermediate momentum levels, d(var_zt)/dz is multiplied by ( K_zm + nu ).
  // Then, the derivative of the whole mathematical expression is taken over
  // the central thermodynamic level, which yields the desired result.
  //
  // --var_zt------------------------------------------------- t(k+1)
  //
  // ==========d(var_zt)/dz==(K_zm+nu)======================== m(k+1)
  //
  // --var_zt-------------------d[(K_zm+nu)*d(var_zt)/dz]/dz-- t(k)
  //
  // ==========d(var_zt)/dz==(K_zm+nu)======================== m(k)
  //
  // --var_zt------------------------------------------------- t(k-1)
  //
  // The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
  // with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
  // The letter "t" is used for thermodynamic levels and the letter "m" is used
  // for momentum levels.
  //
  // invrs_dzt(k)   = 1 / ( zm(k+1) - zm(k) )
  // invrs_dzm(k+1) = 1 / ( zt(k+1) - zt(k) )
  // invrs_dzm(k)   = 1 / ( zt(k) - zt(k-1) )
  //
  // Note:  This function only computes the general implicit form:
  //        + d [ ( K_zm + nu ) * d( var_zt(t+1) )/dz ] / dz.
  //        For a Crank-Nicholson scheme, the left-hand side result of this
  //        function will have to be multiplied by (1/2).  For a
  //        Crank-Nicholson scheme, the right-hand side (explicit) component
  //        needs to be computed by multiplying the left-hand side results by
  //        (1/2), reversing the sign on each left-hand side element, and then
  //        multiplying each element by the appropriate var_zt(t) value from
  //        the appropriate vertical level.
  //
  //
  // Boundary Conditions:
  //
  // 1) Zero-flux boundary conditions.
  //    This function is set up to use zero-flux boundary conditions at both
  //    the lower boundary level and the upper boundary level.  The flux, F,
  //    is the amount of var_zt flowing normal through the boundary per unit
  //    time per unit surface area.  The derivative of the flux effects the
  //    time-tendency of var_zt, such that:
  //
  //    d(var_zt)/dt = -dF/dz.
  //
  //    For the 2nd-order eddy-diffusion term, +d[(K_zm+nu)*d(var_zt)/dz]/dz,
  //    the flux is:
  //
  //    F = -(K_zm+nu)*d(var_zt)/dz.
  //
  //    In order to have zero-flux boundary conditions, the derivative of
  //    var_zt, d(var_zt)/dz, needs to equal 0 at both the lower boundary and
  //    the upper boundary.
  //
  //    In order to discretize the lower boundary condition, consider a new
  //    level outside the model (thermodynamic level 0) just below the lower
  //    boundary level (thermodynamic level 1).  The value of var_zt at the
  //    level just outside the model is defined to be the same as the value of
  //    var_zt at the lower boundary level.  Therefore, the value of
  //    d(var_zt)/dz between the level just outside the model and the lower
  //    boundary level is 0, satisfying the zero-flux boundary condition.  The
  //    other value for d(var_zt)/dz (between thermodynamic level 2 and
  //    thermodynamic level 1) is taken over the intermediate momentum level
  //    (momentum level 2), where it is multiplied by the factor
  //    ( K_zm(2) + nu ).  Then, the derivative of the whole expression is
  //    taken over the central thermodynamic level.
  //
  //    -var_zt(2)-------------------------------------------- t(2)
  //
  //    ==========d(var_zt)/dz=======(K_zm(2)+nu)============= m(2)
  //
  //    -var_zt(1)---------------d[(K_zm+nu)*d(var_zt)/dz]/dz- t(1)
  //
  //    =========[d(var_zt)/dz = 0]==(K_zm(1)+nu)============= m(1) Boundary
  //
  //    -[var_zt(0) = var_zt(1)]-----(level outside model)---- t(0)
  //
  //    The result is dependent only on values of K_zm found at momentum
  //    level 2 and values of var_zt found at thermodynamic levels 1 and 2.
  //    Thus, it only affects 2 diagonals on the left-hand side matrix.
  //
  //    The same method can be used to discretize the upper boundary by
  //    considering a new level outside the model just above the upper boundary
  //    level.
  //
  // 2) Fixed-point boundary conditions.
  //    Many equations in the model use fixed-point boundary conditions rather
  //    than zero-flux boundary conditions.  This means that the value of
  //    var_zt stays the same over the course of the timestep at the lower
  //    boundary, as well as at the upper boundary.
  //
  //    In order to discretize the boundary conditions for equations requiring
  //    fixed-point boundary conditions, either:
  //    a) in the parent subroutine or function (that calls this function),
  //       loop over all vertical levels from the second-lowest to the
  //       second-highest, ignoring the boundary levels.  Then set the values
  //       at the boundary levels in the parent subroutine; or
  //    b) in the parent subroutine or function, loop over all vertical levels
  //       and then overwrite the results at the boundary levels.
  //
  //    Either way, at the boundary levels, an array with a value of 1 at the
  //    main diagonal on the left-hand side and with values of 0 at all other
  //    diagonals on the left-hand side will preserve the right-hand side value
  //    at that level, thus satisfying the fixed-point boundary conditions.
  //
  //
  // Conservation Properties:
  //
  // When zero-flux boundary conditions are used, this technique of
  // discretizing the eddy diffusion term leads to conservative differencing.
  // When conservative differencing is in place, the column totals for each
  // column in the left-hand side matrix (for the eddy diffusion term) should
  // be equal to 0.  This ensures that the total amount of the quantity var_zt
  // over the entire vertical domain is being conserved, meaning that nothing
  // is lost due to diffusional effects.
  //
  // To see that this conservation law is satisfied, compute the eddy diffusion
  // of var_zt and integrate vertically.  In discretized matrix notation (where
  // "i" stands for the matrix column and "j" stands for the matrix row):
  //
  //  0 = Sum_j Sum_i ( 1/invrs_dzt )_i
  //                     ( invrs_dzt * ((K_zm+nu)*invrs_dzm) )_ij (var_zt)_j.
  //
  // The left-hand side matrix, ( invrs_dzt * ((K_zm+nu)*invrs_dzm) )_ij, is
  // partially written below.  The sum over i in the above equation removes
  // invrs_dzt everywhere from the matrix below.  The sum over j leaves the
  // column totals that are desired.
  //
  // Left-hand side matrix contributions from eddy diffusion term; first four
  // vertical levels:
  //
  //     -------------------------------------------------------------------------->
  // k=1 | +invrs_dzt(k)          -invrs_dzt(k)                          0
  //    |   *(K_zm(k+1)+nu)        *(K_zm(k+1)+nu)
  //    |     *invrs_dzm(k+1)        *invrs_dzm(k+1)
  //    |
  // k=2 | -invrs_dzt(k)          +invrs_dzt(k)              -invrs_dzt(k)
  //    |   *(K_zm(k)+nu)          *[ (K_zm(k+1)+nu)          *(K_zm(k+1)+nu)
  //    |     *invrs_dzm(k)            *invrs_dzm(k+1)          *invrs_dzm(k+1)
  //    |                            +(K_zm(k)+nu)
  //    |                              *invrs_dzm(k) ]
  //    |
  // k=3 |         0              -invrs_dzt(k)              +invrs_dzt(k)
  //    |                          *(K_zm(k)+nu)              *[ (K_zm(k+1)+nu)
  //    |                            *invrs_dzm(k)                *invrs_dzm(k+1)
  //    |                                                       +(K_zm(k)+nu)
  //    |                                                         *invrs_dzm(k) ]
  //    |
  // k=4 |         0                        0                -invrs_dzt(k)
  //    |                                                     *(K_zm(k)+nu)
  //    |                                                       *invrs_dzm(k)
  //   \ /
  //
  // Note:  The superdiagonal term from level 3 and both the main diagonal and
  //        superdiagonal terms from level 4 are not shown on this diagram.
  //
  // Note:  The matrix shown is a tridiagonal matrix.  For a band diagonal
  //        matrix (with 5 diagonals), there would be an extra row between each
  //        of the rows shown and an extra column between each of the columns
  //        shown.  However, for the purposes of the var_zt eddy diffusion
  //        term, those extra row and column values are all 0, and the
  //        conservation properties of the matrix aren't effected.
  //
  // If fixed-point boundary conditions are used, the matrix entries at
  // level 1 (k=1) read:  1   0   0; which means that conservative differencing
  // is not in play.  The total amount of var_zt over the entire vertical
  // domain is not being conserved, as amounts of var_zt may be fluxed out
  // through the upper boundary or lower boundary through the effects of
  // diffusion.
  //
  // Brian Griffin.  April 26, 2008.
  // References:
  // None
  //-----------------------------------------------------------------------
  // Constant parameters
  //------------------------ Input Variables ------------------------
  //------------------------ Return Variable ------------------------
  //------------------------ Local Variable ------------------------
  constexpr int kp1_tdiag = 1; // Thermodynamic superdiagonal index.
  constexpr int k_tdiag   = 2; // Thermodynamic main diagonal index.
  constexpr int km1_tdiag = 3; // Thermodynamic subdiagonal index.
  int i;
  int k; // Vertical level index
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_upwind(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "diffusion_zt_lhs::lhs_upwind"), 3, ngrdcol,
      nzt); // LHS coefficient diffusion term due to upwind  [1/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> K_zm_nu(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "diffusion_zt_lhs::K_zm_nu"), ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rho_K_zm_nu(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "diffusion_zt_lhs::rho_K_zm_nu"), ngrdcol,
      nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> drhoKdz_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "diffusion_zt_lhs::drhoKdz_zt"), ngrdcol,
      nzt);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> ddzm_rho_K_zm_nu(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "diffusion_zt_lhs::ddzm_rho_K_zm_nu"),
      ngrdcol, nzt);
  const auto gr_grid_dir  = gr.grid_dir;
  const auto gr_invrs_dzm = gr.invrs_dzm;
  const auto gr_invrs_dzt = gr.invrs_dzt;

  //------------------------ Begin code ------------------------

  Kokkos::parallel_for(
      "diffusion_zt_lhs::loop1", std::max(0, nzm) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzm);
        remaining /= (nzm);
        K_zm_nu(i - 1, k - 1) = K_zm(i - 1, k - 1) + nu(i - 1);
      });

  if (l_upwind_Kh_dp_term) {

    Kokkos::parallel_for(
        "diffusion_zt_lhs::loop2", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          rho_K_zm_nu(i - 1, k - 1) = rho_ds_zm(i - 1, k - 1) * K_zm_nu(i - 1, k - 1);
        });

    // calculate the dKh_zt/dz
    Kokkos::deep_copy(ddzm_rho_K_zm_nu, ddzm(nzm, nzt, ngrdcol, gr, rho_K_zm_nu));

    Kokkos::parallel_for(
        "diffusion_zt_lhs::loop3", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzt);
          remaining /= (nzt);
          drhoKdz_zt(i - 1, k - 1) =
              -invrs_rho_ds_ztzxt(i - 1, k - 1) * ddzm_rho_K_zm_nu(i - 1, k - 1);
        });

    // extra terms with upwind scheme
    // k = 1 (bottom level); lower boundary level
    Kokkos::parallel_for(
        "diffusion_zt_lhs::loop4", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          lhs_upwind(kp1_tdiag - 1, i - 1, 1 - 1) =
              +gr_grid_dir *
              ((gr_grid_dir * drhoKdz_zt(i - 1, 1 - 1)) < (zero)
                   ? (gr_grid_dir * drhoKdz_zt(i - 1, 1 - 1))
                   : (zero)) *
              gr_invrs_dzm(i - 1, 2 - 1);
          lhs_upwind(k_tdiag - 1, i - 1, 1 - 1) = -gr_grid_dir *
                                                  ((gr_grid_dir * drhoKdz_zt(i - 1, 1 - 1)) < (zero)
                                                       ? (gr_grid_dir * drhoKdz_zt(i - 1, 1 - 1))
                                                       : (zero)) *
                                                  gr_invrs_dzm(i - 1, 2 - 1);
          lhs_upwind(km1_tdiag - 1, i - 1, 1 - 1) = zero;
        });

    // Most of the interior model; normal conditions.
    Kokkos::parallel_for(
        "diffusion_zt_lhs::loop5", std::max(0, ((nzt - 1 - (2)) / (1) + 1)) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 2 + remaining % (((nzt - 1 - (2)) / (1) + 1));
          remaining /= (((nzt - 1 - (2)) / (1) + 1));

          // Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
          lhs_upwind(kp1_tdiag - 1, i - 1, k - 1) =
              +gr_grid_dir *
              ((gr_grid_dir * drhoKdz_zt(i - 1, k - 1)) < (zero)
                   ? (gr_grid_dir * drhoKdz_zt(i - 1, k - 1))
                   : (zero)) *
              gr_invrs_dzm(i - 1, k + 1 - 1);

          // Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs_upwind(k_tdiag - 1, i - 1, k - 1) =
              -gr_grid_dir *
                  ((gr_grid_dir * drhoKdz_zt(i - 1, k - 1)) < (zero)
                       ? (gr_grid_dir * drhoKdz_zt(i - 1, k - 1))
                       : (zero)) *
                  gr_invrs_dzm(i - 1, k + 1 - 1) +
              gr_grid_dir *
                  ((gr_grid_dir * drhoKdz_zt(i - 1, k - 1)) > (zero)
                       ? (gr_grid_dir * drhoKdz_zt(i - 1, k - 1))
                       : (zero)) *
                  gr_invrs_dzm(i - 1, k - 1);

          // Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
          lhs_upwind(km1_tdiag - 1, i - 1, k - 1) =
              -gr_grid_dir *
              ((gr_grid_dir * drhoKdz_zt(i - 1, k - 1)) > (zero)
                   ? (gr_grid_dir * drhoKdz_zt(i - 1, k - 1))
                   : (zero)) *
              gr_invrs_dzm(i - 1, k - 1);
        });

    // k = nzt (top level); upper boundary level.
    // Only relevant if zero-flux boundary conditions are used.
    Kokkos::parallel_for(
        "diffusion_zt_lhs::loop6", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          lhs_upwind(kp1_tdiag - 1, i - 1, nzt - 1) = zero;
          lhs_upwind(k_tdiag - 1, i - 1, nzt - 1) =
              +gr_grid_dir *
              ((gr_grid_dir * drhoKdz_zt(i - 1, nzt - 1)) > (zero)
                   ? (gr_grid_dir * drhoKdz_zt(i - 1, nzt - 1))
                   : (zero)) *
              gr_invrs_dzm(i - 1, nzm - 1 - 1);
          lhs_upwind(km1_tdiag - 1, i - 1, nzt - 1) =
              -gr_grid_dir *
              ((gr_grid_dir * drhoKdz_zt(i - 1, nzt - 1)) > (zero)
                   ? (gr_grid_dir * drhoKdz_zt(i - 1, nzt - 1))
                   : (zero)) *
              gr_invrs_dzm(i - 1, nzm - 1 - 1);
        });
  }

  // k = 1 (bottom level); lower boundary level.
  // Only relevant if zero-flux boundary conditions are used.
  // These k=1 lines currently do not have any effect on model results.
  // These k=1 level of this "lhs" array is not fed into
  // the final LHS matrix that will be used to solve for the next timestep.

  if (!l_upwind_Kh_dp_term) {

    Kokkos::parallel_for(
        "diffusion_zt_lhs::loop7", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);

          // Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
          lhs(kp1_tdiag - 1, i - 1, 1 - 1) = -gr_invrs_dzt(i - 1, 1 - 1) *
                                             invrs_rho_ds_ztzxt(i - 1, 1 - 1) *
                                             (K_zm(i - 1, 2 - 1) + nu(i - 1)) *
                                             rho_ds_zm(i - 1, 2 - 1) * gr_invrs_dzm(i - 1, 2 - 1);

          // Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs(k_tdiag - 1, i - 1, 1 - 1) = +gr_invrs_dzt(i - 1, 1 - 1) *
                                           invrs_rho_ds_ztzxt(i - 1, 1 - 1) *
                                           (K_zm(i - 1, 2 - 1) + nu(i - 1)) *
                                           rho_ds_zm(i - 1, 2 - 1) * gr_invrs_dzm(i - 1, 2 - 1);

          // Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
          lhs(km1_tdiag - 1, i - 1, 1 - 1) = zero;
        });

  } else {

    Kokkos::parallel_for(
        "diffusion_zt_lhs::loop8", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);

          // Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
          lhs(kp1_tdiag - 1, i - 1, 1 - 1) = -gr_invrs_dzt(i - 1, 1 - 1) *
                                                 (K_zt(i - 1, 1 - 1) + nu(i - 1)) *
                                                 gr_invrs_dzm(i - 1, 2 - 1) +
                                             lhs_upwind(kp1_tdiag - 1, i - 1, 1 - 1);

          // Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs(k_tdiag - 1, i - 1, 1 - 1) = +gr_invrs_dzt(i - 1, 1 - 1) *
                                               (K_zt(i - 1, 1 - 1) + nu(i - 1)) *
                                               gr_invrs_dzm(i - 1, 2 - 1) +
                                           lhs_upwind(k_tdiag - 1, i - 1, 1 - 1);

          // Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
          lhs(km1_tdiag - 1, i - 1, 1 - 1) = zero + lhs_upwind(km1_tdiag - 1, i - 1, 1 - 1);
        });
  }

  if (!l_upwind_Kh_dp_term) {

    // Most of the interior model; normal conditions.
    Kokkos::parallel_for(
        "diffusion_zt_lhs::loop9", std::max(0, ((nzt - 1 - (2)) / (1) + 1)) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 2 + remaining % (((nzt - 1 - (2)) / (1) + 1));
          remaining /= (((nzt - 1 - (2)) / (1) + 1));

          // Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
          lhs(kp1_tdiag - 1, i - 1, k - 1) =
              -gr_invrs_dzt(i - 1, k - 1) * invrs_rho_ds_ztzxt(i - 1, k - 1) *
              K_zm_nu(i - 1, k + 1 - 1) * rho_ds_zm(i - 1, k + 1 - 1) *
              gr_invrs_dzm(i - 1, k + 1 - 1);

          // Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs(k_tdiag - 1, i - 1, k - 1) =
              +gr_invrs_dzt(i - 1, k - 1) * invrs_rho_ds_ztzxt(i - 1, k - 1) *
              (K_zm_nu(i - 1, k + 1 - 1) * rho_ds_zm(i - 1, k + 1 - 1) *
                   gr_invrs_dzm(i - 1, k + 1 - 1) +
               K_zm_nu(i - 1, k - 1) * rho_ds_zm(i - 1, k - 1) * gr_invrs_dzm(i - 1, k - 1));

          // Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
          lhs(km1_tdiag - 1, i - 1, k - 1) =
              -gr_invrs_dzt(i - 1, k - 1) * invrs_rho_ds_ztzxt(i - 1, k - 1) *
              K_zm_nu(i - 1, k - 1) * rho_ds_zm(i - 1, k - 1) * gr_invrs_dzm(i - 1, k - 1);
        });

  } else {

    Kokkos::parallel_for(
        "diffusion_zt_lhs::loop10", std::max(0, ((nzt - 1 - (2)) / (1) + 1)) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 2 + remaining % (((nzt - 1 - (2)) / (1) + 1));
          remaining /= (((nzt - 1 - (2)) / (1) + 1));

          // Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
          lhs(kp1_tdiag - 1, i - 1, k - 1) = -gr_invrs_dzt(i - 1, k - 1) *
                                                 (K_zt(i - 1, k - 1) + nu(i - 1)) *
                                                 gr_invrs_dzm(i - 1, k + 1 - 1) +
                                             lhs_upwind(kp1_tdiag - 1, i - 1, k - 1);

          // Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs(k_tdiag - 1, i - 1, k - 1) =
              +gr_invrs_dzt(i - 1, k - 1) * (K_zt(i - 1, k - 1) + nu(i - 1)) *
                  (gr_invrs_dzm(i - 1, k + 1 - 1) + gr_invrs_dzm(i - 1, k - 1)) +
              lhs_upwind(k_tdiag - 1, i - 1, k - 1);

          // Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
          lhs(km1_tdiag - 1, i - 1, k - 1) = -gr_invrs_dzt(i - 1, k - 1) *
                                                 (K_zt(i - 1, k - 1) + nu(i - 1)) *
                                                 gr_invrs_dzm(i - 1, k - 1) +
                                             lhs_upwind(km1_tdiag - 1, i - 1, k - 1);
        });
  }

  // k = nzt (top level); upper boundary level.
  // Only relevant if zero-flux boundary conditions are used.

  if (!l_upwind_Kh_dp_term) {

    Kokkos::parallel_for(
        "diffusion_zt_lhs::loop11", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);

          // Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
          lhs(kp1_tdiag - 1, i - 1, nzt - 1) = zero;

          // Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs(k_tdiag - 1, i - 1, nzt - 1) =
              +gr_invrs_dzt(i - 1, nzt - 1) * invrs_rho_ds_ztzxt(i - 1, nzt - 1) *
              (K_zm(i - 1, nzm - 1 - 1) + nu(i - 1)) * rho_ds_zm(i - 1, nzm - 1 - 1) *
              gr_invrs_dzm(i - 1, nzm - 1 - 1);

          // Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
          lhs(km1_tdiag - 1, i - 1, nzt - 1) =
              -gr_invrs_dzt(i - 1, nzt - 1) * invrs_rho_ds_ztzxt(i - 1, nzt - 1) *
              (K_zm(i - 1, nzm - 1 - 1) + nu(i - 1)) * rho_ds_zm(i - 1, nzm - 1 - 1) *
              gr_invrs_dzm(i - 1, nzm - 1 - 1);
        });

  } else {

    Kokkos::parallel_for(
        "diffusion_zt_lhs::loop12", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);

          // Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
          lhs(kp1_tdiag - 1, i - 1, nzt - 1) = zero + lhs_upwind(kp1_tdiag - 1, i - 1, nzt - 1);

          // Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs(k_tdiag - 1, i - 1, nzt - 1) = +gr_invrs_dzt(i - 1, nzt - 1) *
                                                 (K_zt(i - 1, nzt - 1) + nu(i - 1)) *
                                                 gr_invrs_dzm(i - 1, nzm - 1 - 1) +
                                             lhs_upwind(k_tdiag - 1, i - 1, nzt - 1);

          // Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
          lhs(km1_tdiag - 1, i - 1, nzt - 1) = -gr_invrs_dzt(i - 1, nzt - 1) *
                                                   (K_zt(i - 1, nzt - 1) + nu(i - 1)) *
                                                   gr_invrs_dzm(i - 1, nzm - 1 - 1) +
                                               lhs_upwind(km1_tdiag - 1, i - 1, nzt - 1);
        });
  }

  return;
}

//=============================================================================
// K_zt: Coef. of eddy diffusivity at thermo. levels   [m^2/s]
// K_zm: Coef. of eddy diffusivity at momentum levels [m^2/s]
// nu: Background constant coef. of eddy diffusivity [m^2/s]
// invrs_rho_ds_zm: Inverse dry, static density on momentum levels [m^3/kg]
// rho_ds_zt: Dry, static density on thermodynamic levels [kg/m^3]
// lhs: LHS coefficient of diffusion term    [1/s]
void diffusion_zm_lhs(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D K_zt, Real2D K_zm, // In
                      Real1D nu,                                                               // In
                      Real2D invrs_rho_ds_zm, Real2D rho_ds_zt,                                // In
                      Real3D lhs) // Out
{
  // Description:
  // Vertical eddy diffusion of var_zm:  implicit portion of the code.
  //
  // The variable "var_zm" stands for a variable that is located at momentum
  // grid levels.
  //
  // The d(var_zm)/dt equation contains an eddy diffusion term:
  //
  // + d [ ( K_zt + nu ) * d(var_zm)/dz ] / dz.
  //
  // This term is usually solved for completely implicitly, such that:
  //
  // + d [ ( K_zt + nu ) * d( var_zm(t+1) )/dz ] / dz.
  //
  // However, when a Crank-Nicholson scheme is used, the eddy diffusion term
  // has both implicit and explicit components, such that:
  //
  // + (1/2) * d [ ( K_zt + nu ) * d( var_zm(t+1) )/dz ] / dz
  //    + (1/2) * d [ ( K_zt + nu ) * d( var_zm(t) )/dz ] / dz;
  //
  // for which the implicit component is:
  //
  // + (1/2) * d [ ( K_zt + nu ) * d( var_zm(t+1) )/dz ] / dz.
  //
  // Note:  When the implicit term is brought over to the left-hand side,
  //        the sign is reversed and the leading "+" in front of the term
  //        is changed to a "-".
  //
  // Timestep index (t) stands for the index of the current timestep, while
  // timestep index (t+1) stands for the index of the next timestep, which is
  // being advanced to in solving the d(var_zm)/dt equation.
  //
  // The implicit portion of this term is discretized as follows:
  //
  // The values of var_zm are found on the momentum levels, while the values of
  // K_zt are found on the thermodynamic levels.  The derivatives (d/dz) of
  // var_zm are taken over the intermediate thermodynamic levels.  At the
  // intermediate thermodynamic levels, d(var_zm)/dz is multiplied by
  // ( K_zt + nu ).  Then, the derivative of the whole mathematical expression
  // is taken over the central momentum level, which yields the desired result.
  //
  // ==var_zm================================================= m(k+1)
  //
  // ----------d(var_zm)/dz--(K_zt+nu)------------------------ t(k)
  //
  // ==var_zm===================d[(K_zt+nu)*d(var_zm)/dz]/dz== m(k)
  //
  // ----------d(var_zm)/dz--(K_zt+nu)------------------------ t(k-1)
  //
  // ==var_zm================================================= m(k-1)
  //
  // The vertical indices m(k+1), t(k), m(k), t(k-1), and m(k-1) correspond
  // with altitudes zm(k+1), zt(k), zm(k), zt(k-1), and zm(k-1), respectively.
  // The letter "t" is used for thermodynamic levels and the letter "m" is used
  // for momentum levels.
  //
  // invrs_dzm(k)   = 1 / ( zt(k) - zt(k-1) )
  // invrs_dzt(k)   = 1 / ( zm(k+1) - zm(k) )
  // invrs_dzt(k-1) = 1 / ( zm(k) - zm(k-1) )
  //
  // Note:  This function only computes the general implicit form:
  //        + d [ ( K_zt + nu ) * d( var_zm(t+1) )/dz ] / dz.
  //        For a Crank-Nicholson scheme, the left-hand side result of this
  //        function will have to be multiplied by (1/2).  For a
  //        Crank-Nicholson scheme, the right-hand side (explicit) component
  //        needs to be computed by multiplying the left-hand side results by
  //        (1/2), reversing the sign on each left-hand side element, and then
  //        multiplying each element by the appropriate var_zm(t) value from
  //        the appropriate vertical level.
  //
  //
  // Boundary Conditions:
  //
  // 1) Zero-flux boundary conditions.
  //    This function is set up to use zero-flux boundary conditions at both
  //    the lower boundary level and the upper boundary level.  The flux, F,
  //    is the amount of var_zm flowing normal through the boundary per unit
  //    time per unit surface area.  The derivative of the flux effects the
  //    time-tendency of var_zm, such that:
  //
  //    d(var_zm)/dt = -dF/dz.
  //
  //    For the 2nd-order eddy-diffusion term, +d[(K_zt+nu)*d(var_zm)/dz]/dz,
  //    the flux is:
  //
  //    F = -(K_zt+nu)*d(var_zm)/dz.
  //
  //    In order to have zero-flux boundary conditions, the derivative of
  //    var_zm, d(var_zm)/dz, needs to equal 0 at both the lower boundary and
  //    the upper boundary.
  //
  //    In order to discretize the lower boundary condition, consider a new
  //    level outside the model (momentum level 0) just below the lower
  //    boundary level (momentum level 1).  The value of var_zm at the level
  //    just outside the model is defined to be the same as the value of var_zm
  //    at the lower boundary level.  Therefore, the value of d(var_zm)/dz
  //    between the level just outside the model and the lower boundary level
  //    is 0, satisfying the zero-flux boundary condition.  The other value for
  //    d(var_zm)/dz (between momentum level 2 and momentum level 1) is taken
  //    over the intermediate thermodynamic level (thermodynamic level 1),
  //    where it is multiplied by the factor ( K_zt(1) + nu ).  Then, the
  //    derivative of the whole expression is taken over the central momentum
  //    level.
  //
  //    =var_zm(2)============================================ m(2)
  //
  //    ----------d(var_zm)/dz==(K_zt(1)+nu)------------------ t(1)
  //
  //    =var_zm(1)===============d[(K_zt+nu)*d(var_zm)/dz]/dz= m(1) Boundary
  //
  //    ----------[d(var_zm)/dz = 0]-(level outside model)---- t(0)
  //
  //    =[var_zm(0) = var_zm(1)]=====(level outside model)==== m(0)
  //
  //    The result is dependent only on values of K_zt found at thermodynamic
  //    level 1 and values of var_zm found at momentum levels 1 and 2.  Thus,
  //    it only affects 2 diagonals on the left-hand side matrix.
  //
  //    The same method can be used to discretize the upper boundary by
  //    considering a new level outside the model just above the upper boundary
  //    level.
  //
  // 2) Fixed-point boundary conditions.
  //    Many equations in the model use fixed-point boundary conditions rather
  //    than zero-flux boundary conditions.  This means that the value of
  //    var_zm stays the same over the course of the timestep at the lower
  //    boundary, as well as at the upper boundary.
  //
  //    In order to discretize the boundary conditions for equations requiring
  //    fixed-point boundary conditions, either:
  //    a) in the parent subroutine or function (that calls this function),
  //       loop over all vertical levels from the second-lowest to the
  //       second-highest, ignoring the boundary levels.  Then set the values
  //       at the boundary levels in the parent subroutine; or
  //    b) in the parent subroutine or function, loop over all vertical levels
  //       and then overwrite the results at the boundary levels.
  //
  //    Either way, at the boundary levels, an array with a value of 1 at the
  //    main diagonal on the left-hand side and with values of 0 at all other
  //    diagonals on the left-hand side will preserve the right-hand side value
  //    at that level, thus satisfying the fixed-point boundary conditions.
  //
  //
  // Conservation Properties:
  //
  // When zero-flux boundary conditions are used, this technique of
  // discretizing the eddy diffusion term leads to conservative differencing.
  // When conservative differencing is in place, the column totals for each
  // column in the left-hand side matrix (for the eddy diffusion term) should
  // be equal to 0.  This ensures that the total amount of the quantity var_zm
  // over the entire vertical domain is being conserved, meaning that nothing
  // is lost due to diffusional effects.
  //
  // To see that this conservation law is satisfied, compute the eddy diffusion
  // of var_zm and integrate vertically.  In discretized matrix notation (where
  // "i" stands for the matrix column and "j" stands for the matrix row):
  //
  //  0 = Sum_j Sum_i ( 1/invrs_dzm )_i
  //                     ( invrs_dzm * ((K_zt+nu)*invrs_dzt) )_ij (var_zm)_j.
  //
  // The left-hand side matrix, ( invrs_dzm * ((K_zt+nu)*invrs_dzt) )_ij, is
  // partially written below.  The sum over i in the above equation removes
  // invrs_dzm everywhere from the matrix below.  The sum over j leaves the
  // column totals that are desired.
  //
  // Left-hand side matrix contributions from eddy diffusion term; first four
  // vertical levels:
  //
  //     ---------------------------------------------------------------------->
  // k=1 | +invrs_dzm(k)          -invrs_dzm(k)                      0
  //    |   *(K_zt(k)+nu)          *(K_zt(k)+nu)
  //    |     *invrs_dzt(k)          *invrs_dzt(k)
  //    |
  // k=2 | -invrs_dzm(k)          +invrs_dzm(k)            -invrs_dzm(k)
  //    |   *(K_zt(k-1)+nu)        *[ (K_zt(k)+nu)          *(K_zt(k)+nu)
  //    |     *invrs_dzt(k-1)          *invrs_dzt(k)          *invrs_dzt(k)
  //    |                            +(K_zt(k-1)+nu)
  //    |                              *invrs_dzt(k-1) ]
  //    |
  // k=3 |          0             -invrs_dzm(k)            +invrs_dzm(k)
  //    |                          *(K_zt(k-1)+nu)          *[ (K_zt(k)+nu)
  //    |                            *invrs_dzt(k-1)            *invrs_dzt(k)
  //    |                                                     +(K_zt(k-1)+nu)
  //    |                                                       *invrs_dzt(k-1) ]
  //    |
  // k=4 |          0                       0              -invrs_dzm(k)
  //    |                                                   *(K_zt(k-1)+nu)
  //    |                                                     *invrs_dzt(k-1)
  //   \ /
  //
  // Note:  The superdiagonal term from level 3 and both the main diagonal and
  //        superdiagonal terms from level 4 are not shown on this diagram.
  //
  // Note:  The matrix shown is a tridiagonal matrix.  For a band diagonal
  //        matrix (with 5 diagonals), there would be an extra row between each
  //        of the rows shown and an extra column between each of the columns
  //        shown.  However, for the purposes of the var_zm eddy diffusion
  //        term, those extra row and column values are all 0, and the
  //        conservation properties of the matrix aren't effected.
  //
  // If fixed-point boundary conditions are used, the matrix entries at
  // level 1 (k=1) read:  1   0   0; which means that conservative differencing
  // is not in play.  The total amount of var_zm over the entire vertical
  // domain is not being conserved, as amounts of var_zm may be fluxed out
  // through the upper boundary or lower boundary through the effects of
  // diffusion.
  //
  // Brian Griffin.  April 26, 2008.
  // References:
  // None
  //-----------------------------------------------------------------------
  // Constant parameters
  // Input Variables
  // Return Variable
  // Local Variable
  constexpr int kp1_mdiag = 1; // Momentum superdiagonal index.
  constexpr int k_mdiag   = 2; // Momentum main diagonal index.
  constexpr int km1_mdiag = 3; // Momentum subdiagonal index.
  int i;
  int k; // Vertical level index
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_upwind(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "diffusion_zm_lhs::lhs_upwind"), 3, ngrdcol,
      nzm); // LHS coefficient diffusion term due to upwind  [1/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> drhoKdz_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "diffusion_zm_lhs::drhoKdz_zm"), ngrdcol,
      nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> ddzt_rho_K_zt_nu(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "diffusion_zm_lhs::ddzt_rho_K_zt_nu"),
      ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rho_K_zt_nu(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "diffusion_zm_lhs::rho_K_zt_nu"), ngrdcol,
      nzt);
  const auto gr_grid_dir  = gr.grid_dir;
  const auto gr_invrs_dzt = gr.invrs_dzt;
  const auto gr_invrs_dzm = gr.invrs_dzm;

  //------------Begin code------------------------------

  Kokkos::parallel_for(
      "diffusion_zm_lhs::loop1", std::max(0, nzt) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzt);
        remaining /= (nzt);
        rho_K_zt_nu(i - 1, k - 1) = rho_ds_zt(i - 1, k - 1) * (K_zt(i - 1, k - 1) + nu(i - 1));
      });

  // calculate the dKh_zm/dz
  Kokkos::deep_copy(ddzt_rho_K_zt_nu, ddzt(nzm, nzt, ngrdcol, gr, rho_K_zt_nu));

  Kokkos::parallel_for(
      "diffusion_zm_lhs::loop2", std::max(0, nzm) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzm);
        remaining /= (nzm);
        drhoKdz_zm(i - 1, k - 1) = -invrs_rho_ds_zm(i - 1, k - 1) * ddzt_rho_K_zt_nu(i - 1, k - 1);
      });

  // extra terms with upwind scheme
  // k = 1 (bottom level); lowere boundary level
  Kokkos::parallel_for(
      "diffusion_zm_lhs::loop3", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        lhs_upwind(kp1_mdiag - 1, i - 1, 1 - 1) = +gr_grid_dir *
                                                  ((gr_grid_dir * drhoKdz_zm(i - 1, 1 - 1)) < (zero)
                                                       ? (gr_grid_dir * drhoKdz_zm(i - 1, 1 - 1))
                                                       : (zero)) *
                                                  gr_invrs_dzt(i - 1, 1 - 1);
        lhs_upwind(k_mdiag - 1, i - 1, 1 - 1)   = -gr_grid_dir *
                                                  ((gr_grid_dir * drhoKdz_zm(i - 1, 1 - 1)) < (zero)
                                                       ? (gr_grid_dir * drhoKdz_zm(i - 1, 1 - 1))
                                                       : (zero)) *
                                                  gr_invrs_dzt(i - 1, 1 - 1);
        lhs_upwind(km1_mdiag - 1, i - 1, 1 - 1) = zero;
      });

  // Most of the interior model; normal conditions.
  Kokkos::parallel_for(
      "diffusion_zm_lhs::loop4", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        // Momentum superdiagonal: [ x var_zm(k+1,<t+1>) ]
        lhs_upwind(kp1_mdiag - 1, i - 1, k - 1) = +gr_grid_dir *
                                                  ((gr_grid_dir * drhoKdz_zm(i - 1, k - 1)) < (zero)
                                                       ? (gr_grid_dir * drhoKdz_zm(i - 1, k - 1))
                                                       : (zero)) *
                                                  gr_invrs_dzt(i - 1, k - 1);
        // Momentum main diagonal: [ x var_zm(k,<t+1>) ]
        lhs_upwind(k_mdiag - 1, i - 1, k - 1) =
            -gr_grid_dir *
                ((gr_grid_dir * drhoKdz_zm(i - 1, k - 1)) < (zero)
                     ? (gr_grid_dir * drhoKdz_zm(i - 1, k - 1))
                     : (zero)) *
                gr_invrs_dzt(i - 1, k - 1) +
            gr_grid_dir *
                ((gr_grid_dir * drhoKdz_zm(i - 1, k - 1)) > (zero)
                     ? (gr_grid_dir * drhoKdz_zm(i - 1, k - 1))
                     : (zero)) *
                gr_invrs_dzt(i - 1, k - 1 - 1);
        // Momentum subdiagonal: [ x var_zm(k-1,<t+1>) ]
        lhs_upwind(km1_mdiag - 1, i - 1, k - 1) = -gr_grid_dir *
                                                  ((gr_grid_dir * drhoKdz_zm(i - 1, k - 1)) > (zero)
                                                       ? (gr_grid_dir * drhoKdz_zm(i - 1, k - 1))
                                                       : (zero)) *
                                                  gr_invrs_dzt(i - 1, k - 1 - 1);
      });

  // k = nzm (top level); upper boundary level.
  // Only relevant if zero-flux boundary conditions are used.
  Kokkos::parallel_for(
      "diffusion_zm_lhs::loop5", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        lhs_upwind(kp1_mdiag - 1, i - 1, nzm - 1) = zero;
        lhs_upwind(k_mdiag - 1, i - 1, nzm - 1) =
            +gr_grid_dir *
            ((gr_grid_dir * drhoKdz_zm(i - 1, nzm - 1)) > (zero)
                 ? (gr_grid_dir * drhoKdz_zm(i - 1, nzm - 1))
                 : (zero)) *
            gr_invrs_dzt(i - 1, nzt - 1);
        lhs_upwind(km1_mdiag - 1, i - 1, nzm - 1) =
            -gr_grid_dir *
            ((gr_grid_dir * drhoKdz_zm(i - 1, nzm - 1)) > (zero)
                 ? (gr_grid_dir * drhoKdz_zm(i - 1, nzm - 1))
                 : (zero)) *
            gr_invrs_dzt(i - 1, nzt - 1);
      });

  // k = 1; lower boundary level at surface.
  // Only relevant if zero-flux boundary conditions are used.
  // These k=1 lines currently do not have any effect on model results.
  // These k=1 level of this "lhs" array is not fed into
  // the final LHS matrix that will be used to solve for the next timestep.

  if (!l_upwind_Kh_dp_term) {

    Kokkos::parallel_for(
        "diffusion_zm_lhs::loop6", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);

          // Momentum superdiagonal: [ x var_zm(k+1,<t+1>) ]
          lhs(kp1_mdiag - 1, i - 1, 1 - 1) = -gr_invrs_dzm(i - 1, 1 - 1) *
                                             invrs_rho_ds_zm(i - 1, 1 - 1) *
                                             (K_zt(i - 1, 1 - 1) + nu(i - 1)) *
                                             rho_ds_zt(i - 1, 1 - 1) * gr_invrs_dzt(i - 1, 1 - 1);

          // Momentum main diagonal: [ x var_zm(k,<t+1>) ]
          lhs(k_mdiag - 1, i - 1, 1 - 1) = +gr_invrs_dzm(i - 1, 1 - 1) *
                                           invrs_rho_ds_zm(i - 1, 1 - 1) *
                                           (K_zt(i - 1, 1 - 1) + nu(i - 1)) *
                                           rho_ds_zt(i - 1, 1 - 1) * gr_invrs_dzt(i - 1, 1 - 1);

          // Momentum subdiagonal: [ x var_zm(k-1,<t+1>) ]
          lhs(km1_mdiag - 1, i - 1, 1 - 1) = zero;
        });

  } else {

    Kokkos::parallel_for(
        "diffusion_zm_lhs::loop7", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);

          // Momentum superdiagonal: [ x var_zm(k+1,<t+1>) ]
          lhs(kp1_mdiag - 1, i - 1, 1 - 1) = -gr_invrs_dzm(i - 1, 1 - 1) *
                                                 (K_zm(i - 1, 1 - 1) + nu(i - 1)) *
                                                 gr_invrs_dzt(i - 1, 1 - 1) +
                                             lhs_upwind(kp1_mdiag - 1, i - 1, 1 - 1);

          // Momentum main diagonal: [ x var_zm(k,<t+1>) ]
          lhs(k_mdiag - 1, i - 1, 1 - 1) = +gr_invrs_dzm(i - 1, 1 - 1) *
                                               (K_zm(i - 1, 1 - 1) + nu(i - 1)) *
                                               gr_invrs_dzt(i - 1, 1 - 1) +
                                           lhs_upwind(k_mdiag - 1, i - 1, 1 - 1);

          // Momentum subdiagonal: [ x var_zm(k-1,<t+1>) ]
          lhs(km1_mdiag - 1, i - 1, 1 - 1) = zero + lhs_upwind(km1_mdiag - 1, i - 1, 1 - 1);
        });
  }

  // Most of the interior model; normal conditions.
  if (!l_upwind_Kh_dp_term) {

    Kokkos::parallel_for(
        "diffusion_zm_lhs::loop8", std::max(0, nzm - 2) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 2 + remaining % (nzm - 2);
          remaining /= (nzm - 2);

          // Momentum superdiagonal: [ x var_zm(k+1,<t+1>) ]
          lhs(kp1_mdiag - 1, i - 1, k - 1) = -gr_invrs_dzm(i - 1, k - 1) *
                                             invrs_rho_ds_zm(i - 1, k - 1) *
                                             (K_zt(i - 1, k - 1) + nu(i - 1)) *
                                             rho_ds_zt(i - 1, k - 1) * gr_invrs_dzt(i - 1, k - 1);

          // Momentum main diagonal: [ x var_zm(k,<t+1>) ]
          lhs(k_mdiag - 1, i - 1, k - 1) =
              +gr_invrs_dzm(i - 1, k - 1) * invrs_rho_ds_zm(i - 1, k - 1) *
              ((K_zt(i - 1, k - 1) + nu(i - 1)) * rho_ds_zt(i - 1, k - 1) *
                   gr_invrs_dzt(i - 1, k - 1) +
               (K_zt(i - 1, k - 1 - 1) + nu(i - 1)) * rho_ds_zt(i - 1, k - 1 - 1) *
                   gr_invrs_dzt(i - 1, k - 1 - 1));

          // Momentum subdiagonal: [ x var_zm(k-1,<t+1>) ]
          lhs(km1_mdiag - 1, i - 1, k - 1) =
              -gr_invrs_dzm(i - 1, k - 1) * invrs_rho_ds_zm(i - 1, k - 1) *
              (K_zt(i - 1, k - 1 - 1) + nu(i - 1)) * rho_ds_zt(i - 1, k - 1 - 1) *
              gr_invrs_dzt(i - 1, k - 1 - 1);
        });

  } else {

    Kokkos::parallel_for(
        "diffusion_zm_lhs::loop9", std::max(0, nzm - 2) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 2 + remaining % (nzm - 2);
          remaining /= (nzm - 2);

          // Momentum superdiagonal: [ x var_zm(k+1,<t+1>) ]
          lhs(kp1_mdiag - 1, i - 1, k - 1) = -gr_invrs_dzm(i - 1, k - 1) *
                                                 (K_zm(i - 1, k - 1) + nu(i - 1)) *
                                                 gr_invrs_dzt(i - 1, k - 1) +
                                             lhs_upwind(kp1_mdiag - 1, i - 1, k - 1);

          // Momentum main diagonal: [ x var_zm(k,<t+1>) ]
          lhs(k_mdiag - 1, i - 1, k - 1) =
              +gr_invrs_dzm(i - 1, k - 1) * (K_zm(i - 1, k - 1) + nu(i - 1)) *
                  (gr_invrs_dzt(i - 1, k - 1) + gr_invrs_dzt(i - 1, k - 1 - 1)) +
              lhs_upwind(k_mdiag - 1, i - 1, k - 1);

          // Momentum subdiagonal: [ x var_zm(k-1,<t+1>) ]
          lhs(km1_mdiag - 1, i - 1, k - 1) = -gr_invrs_dzm(i - 1, k - 1) *
                                                 (K_zm(i - 1, k - 1) + nu(i - 1)) *
                                                 gr_invrs_dzt(i - 1, k - 1 - 1) +
                                             lhs_upwind(km1_mdiag - 1, i - 1, k - 1);
        });
  }

  // k = nzm (top level); upper boundary level.
  // Only relevant if zero-flux boundary conditions are used.

  if (!l_upwind_Kh_dp_term) {

    Kokkos::parallel_for(
        "diffusion_zm_lhs::loop10", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          // Momentum superdiagonal: [ x var_zm(k+1,<t+1>) ]
          lhs(kp1_mdiag - 1, i - 1, nzm - 1) = zero;

          // Momentum main diagonal: [ x var_zm(k,<t+1>) ]
          lhs(k_mdiag - 1, i - 1, nzm - 1) =
              +gr_invrs_dzm(i - 1, nzm - 1) * invrs_rho_ds_zm(i - 1, nzm - 1) *
              (K_zt(i - 1, nzt - 1) + nu(i - 1)) * rho_ds_zt(i - 1, nzt - 1) *
              gr_invrs_dzt(i - 1, nzt - 1);

          // Momentum subdiagonal: [ x var_zm(k-1,<t+1>) ]
          lhs(km1_mdiag - 1, i - 1, nzm - 1) =
              -gr_invrs_dzm(i - 1, nzm - 1) * invrs_rho_ds_zm(i - 1, nzm - 1) *
              (K_zt(i - 1, nzt - 1) + nu(i - 1)) * rho_ds_zt(i - 1, nzt - 1) *
              gr_invrs_dzt(i - 1, nzt - 1);
        });

  } else {

    Kokkos::parallel_for(
        "diffusion_zm_lhs::loop11", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);

          // Momentum superdiagonal: [ x var_zm(k+1,<t+1>) ]
          lhs(kp1_mdiag - 1, i - 1, nzm - 1) = zero + lhs_upwind(kp1_mdiag - 1, i - 1, nzm - 1);

          // Momentum main diagonal: [ x var_zm(k,<t+1>) ]
          lhs(k_mdiag - 1, i - 1, nzm - 1) = +gr_invrs_dzm(i - 1, nzm - 1) *
                                                 (K_zm(i - 1, nzm - 1) + nu(i - 1)) *
                                                 gr_invrs_dzt(i - 1, nzt - 1) +
                                             lhs_upwind(k_mdiag - 1, i - 1, nzm - 1);

          // Momentum subdiagonal: [ x var_zm(k-1,<t+1>) ]
          lhs(km1_mdiag - 1, i - 1, nzm - 1) = -gr_invrs_dzm(i - 1, nzm - 1) *
                                                   (K_zm(i - 1, nzm - 1) + nu(i - 1)) *
                                                   gr_invrs_dzt(i - 1, nzt - 1) +
                                               lhs_upwind(km1_mdiag - 1, i - 1, nzm - 1);
        });
  }

  return;
}

} // namespace clubb_kokkos
