// Compute PDF-based turbulent advection for variances, covariances and fluxes.
//
// xpyp_term_ta_pdf_lhs and its Godunov variant fill three-band implicit
// operators; xpyp_term_ta_pdf_rhs and its Godunov variant fill explicit
// tendencies. Physical fields use column-first LayoutLeft Views. LHS output
// retains Native (band,column,momentum level) order, with super/main/sub bands.
// Both storage endpoints are zero; grid direction selects upwind differences.
// Generic entry overloads retain caller-owned Views through the existing array
// owner. All arguments are required.
// Description:
// Calculates the turbulent advection term in the predictive equation for a
// variance or covariance where turbulent advection is calculated by
// integrating over the PDF.  This includes the following predictive fields:
// <w'thl'>, <w'rt'>, <rt'^2>, <thl'^2>, and <rt'thl'>, as well as passive
// scalar fields <w'sclr'>, <sclr'^2>, <sclr'rt'>, and <sclr'thl'>.  CLUBB
// does not produce a PDF for horizontal wind components u and v.  However, the
// horizontal wind variances <u'^2> and <v'^2> still use this code, as well.

#include "turbulent_adv_pdf.hpp"
#include "constants_clubb.hpp"

namespace clubb_kokkos {

//=============================================================================
// coef_wpxpyp_implicit: Coef. of <x'y'> in <w'x'y'>; t-levs  [m/s]
// rho_ds_zt: Dry, static density at t-levels      [kg/m^3]
// rho_ds_zm: Dry, static density at m-levs       [kg/m^3]
// invrs_rho_ds_zm: Inv dry, static density at m-levels  [m^3/kg]
// l_upwind_xpyp_turbulent_adv: Flag to use "upwind" discretization
// sgn_turbulent_vel: Sign of the turbulent velocity      [-]
// coef_wpxpyp_implicit_zm: coef_wpxpyp_implicit interp m-levs  [m/s]
// lhs_ta: LHS coefficient of xpyp turbulent advection  [1/s]
// lhs_ta storage: (band,column,momentum level), super/main/sub bands [1/s].
void xpyp_term_ta_pdf_lhs(int nzm, int nzt, int ngrdcol, const grid& gr,
                          Real2D coef_wpxpyp_implicit,        // In
                          Real2D rho_ds_zt, Real2D rho_ds_zm, // In
                          Real2D invrs_rho_ds_zm,             // In
                          bool l_upwind_xpyp_turbulent_adv,   // In
                          Real2D sgn_turbulent_vel,           // In
                          Real2D coef_wpxpyp_implicit_zm,     // In
                          Real3D lhs_ta)                      // Out
{
  // Description:
  // Turbulent advection of <w'x'>, <x'^2>, and <x'y'>:  implicit portion of
  // the code.
  //
  // 1) <w'x'>
  //
  // The d<w'x'>/dt equation contains a turbulent advection term:
  //
  // - (1/rho_ds) * d( rho_ds * <w'^2 x'> )/dz.
  //
  // The value of <w'^2 x'> is found by integrating over the multivariate PDF
  // of w and x, as detailed in function calc_wp2xp_pdf, which is found in
  // module pdf_closure_module in pdf_closure_module.F90.
  //
  // The equation obtained for <w'^2 x'> is written in terms of PDF parameters.
  // Substitutions that are specific to the type of PDF used are made for the
  // PDF parameters in order to write the <w'x'> turbulent advection term in
  // the following form:
  //
  // <w'^2 x'> = coef_wp2xp_implicit * <w'x'> + term_wp2xp_explicit.
  //
  // For the ADG1 PDF, coef_wp2xp_implicit is a_1 * < w'^3 > / < w'^2 >, where
  // a_1 is 1 / ( 1 - sigma_sqd_w ).  The value of term_wp2xp_explicit is 0, as
  // the <w'x'> turbulent advection term is entirely implicit.
  //
  // For the new PDF, the calculations of both coef_wp2xp_implicit and
  // term_wp2xp_explicit are detailed in function calc_coefs_wp2xp_semiimpl,
  // which is found in module new_pdf in new_pdf.F90.
  //
  // For the new hybrid PDF, the calculation of coef_wp2xp_implicit is
  // detailed in function calc_coef_wp2xp_implicit, which is found in module
  // new_hybrid_pdf in new_hybrid_pdf.F90.  The value of term_wp2xp_explicit
  // is 0, as the <w'x'> turbulent advection term is entirely implicit.
  //
  // For explicit turbulent advection, the value of coef_wp2xp_implicit is 0
  // and the value of term_wp2xp_explicit is <w'^2 x'>, as calculated by
  // retaining the equation for <w'^2 x'> that is written in terms of PDF
  // parameters.  This is a general form that can be used with any type of PDF.
  //
  // The <w'x'> turbulent advection term is rewritten as:
  //
  // - (1/rho_ds)
  //   * d( rho_ds * ( coef_wp2xp_implicit * <w'x'> + term_wp2xp_explicit ) )
  //     /dz.
  //
  // The variable <w'x'> is evaluated at the (t+1) timestep, which allows the
  // <w'x'> turbulent advection term to be expressed semi-implicitly as:
  //
  // - (1/rho_ds) * d( rho_ds * coef_wp2xp_implicit * <w'x'>(t+1) )/dz
  // - (1/rho_ds) * d( rho_ds * term_wp2xp_explicit )/dz.
  //
  // The implicit portion of <w'x'> turbulent advection term is:
  //
  // - (1/rho_ds) * d( rho_ds * coef_wp2xp_implicit * <w'x'>(t+1) )/dz.
  //
  // Note:  When the implicit term is brought over to the left-hand side, the
  //        sign is reversed and the leading "-" in front of the implicit
  //        d[ ] / dz term is changed to a "+".
  //
  // The timestep index (t+1) means that the value of <w'x'> being used is from
  // the next timestep, which is being advanced to in solving the d<w'x'>/dt
  // equation.
  //
  // 2) <x'^2>
  //
  // The d<x'^2>/dt equation contains a turbulent advection term:
  //
  // - (1/rho_ds) * d( rho_ds * <w'x'^2> )/dz;
  //
  // The value of <w'x'^2> is found by integrating over the multivariate PDF of
  // w and x, as detailed in function calc_wpxp2_pdf, which is found in module
  // pdf_closure_module in pdf_closure_module.F90.
  //
  // The equation obtained for <w'x'^2> is written in terms of PDF parameters.
  // Substitutions that are specific to the type of PDF used are made for the
  // PDF parameters in order to write the <x'^2> turbulent advection term in
  // the following form:
  //
  // <w'x'^2> = coef_wpxp2_implicit * <x'^2> + term_wpxp2_explicit.
  //
  // For the ADG1 PDF, the value of coef_wpxp2_implicit is
  // (1/3)*beta * a_1 * < w'^3 > / < w'^2 >.  The value of term_wpxp2_explicit
  // is (1-(1/3)*beta) * a_1^2 * < w'x' >^2 * < w'^3 > / < w'^2 >^2, where
  // a_1 is 1 / ( 1 - sigma_sqd_w ).
  //
  // For the new PDF, the calculation of coef_wpxp2_implicit is detailed in
  // function calc_coef_wpxp2_implicit, which is found in module new_pdf in
  // new_pdf.F90.  The value of term_wpxp2_explicit is 0, as the <x'^2>
  // turbulent advection term is entirely implicit.
  //
  // For the new hybrid PDF, the calculation of both coef_wpxp2_implicit and
  // term_wpxp2_explicit are detailed in subroutine calc_coefs_wpxp2_semiimpl,
  // which is found in module new_hybrid_pdf in new_hybrid_pdf.F90.
  //
  // For explicit turbulent advection, the value of coef_wpxp2_implicit is 0
  // and the value of term_wpxp2_explicit is <w'x'^2>, as calculated by
  // retaining the equation for <w'x'^2> that is written in terms of PDF
  // parameters.  This is a general form that can be used with any type of PDF.
  //
  // The <x'^2> turbulent advection term is rewritten as:
  //
  // - (1/rho_ds)
  //   * d( rho_ds * ( coef_wpxp2_implicit * <x'^2> + term_wpxp2_explicit ) )
  //     /dz;
  //
  // The variable <x'^2> is evaluated at the (t+1) timestep, which allows the
  // <x'^2> turbulent advection term to be expressed semi-implicitly as:
  //
  // - (1/rho_ds) * d( rho_ds * coef_wpxp2_implicit * <x'^2>(t+1) )/dz;
  // - (1/rho_ds) * d( rho_ds * term_wpxp2_explicit )/dz.
  //
  // The implicit portion of <x'^2> turbulent advection term is:
  //
  // - (1/rho_ds) * d( rho_ds * coef_wpxp2_implicit * <x'^2>(t+1) )/dz.
  //
  // Note:  When the implicit term is brought over to the left-hand side, the
  //        sign is reversed and the leading "-" in front of all implicit
  //        d[ ] / dz terms is changed to a "+".
  //
  // The timestep index (t+1) means that the value of <x'^2> being used is from
  // the next timestep, which is being advanced to in solving the d<x'^2>/dt
  // equation.
  //
  // 3) <x'y'>
  //
  // The d<x'y'>/dt equation contains a turbulent advection term:
  //
  // - (1/rho_ds) * d( rho_ds * <w'x'y'> )/dz.
  //
  // The value of <w'x'y'> is found by integrating over the multivariate PDF of
  // w, x, and y, as detailed in function calc_wpxpyp_pdf, which is found in
  // module pdf_closure_module in pdf_closure_module.F90.
  //
  // The equation obtained for <w'x'y'> is written in terms of PDF parameters.
  // Substitutions that are specific to the type of PDF used are made for the
  // PDF parameters in order to write the <x'y'> turbulent advection term in
  // the following form:
  //
  // <w'x'y'> = coef_wpxpyp_implicit * <x'y'> + term_wpxpyp_explicit.
  //
  // For the ADG1 PDF, the value of coef_wpxpyp_implicit is
  // (1/3)*beta * a_1 * < w'^3 > / < w'^2 >.  The value of term_wpxpyp_explicit
  // is (1-(1/3)*beta) * a_1^2 * < w'x' > * < w'y' > * < w'^3 > / < w'^2 >^2,
  // where a_1 is 1 / ( 1 - sigma_sqd_w ).
  //
  // For the new PDF, the calculation of both coef_wpxpyp_implicit and
  // term_wpxpyp_explicit are detailed in function calc_coefs_wpxpyp_semiimpl,
  // which is found in module new_pdf in new_pdf.F90.
  //
  // For the new hybrid PDF, the calculation of both coef_wpxpyp_implicit
  // and term_wpxpyp_explicit are detailed in subroutine
  // calc_coefs_wpxpyp_semiimpl, which is found in module new_hybrid_pdf in
  // new_hybrid_pdf.F90.
  //
  // For explicit turbulent advection, the value of coef_wpxpyp_implicit is 0
  // and the value of term_wpxpyp_explicit is <w'x'y'>, as calculated by
  // retaining the equation for <w'x'y'> that is written in terms of PDF
  // parameters.  This is a general form that can be used with any type of PDF.
  //
  // The <x'y'> turbulent advection term is rewritten as:
  //
  // - (1/rho_ds)
  //   * d( rho_ds * ( coef_wpxpyp_implicit * <x'y'> + term_wpxpyp_explicit ) )
  //     /dz.
  //
  // The variable <x'y'> is evaluated at the (t+1) timestep, which allows the
  // <x'y'> turbulent advection term to be expressed semi-implicitly as:
  //
  // - (1/rho_ds) * d( rho_ds * coef_wpxpyp_implicit * <x'y'>(t+1) )/dz
  // - (1/rho_ds) * d( rho_ds * term_wpxpyp_explicit )/dz.
  //
  // The implicit portion of <x'y'> turbulent advection term is:
  //
  // - (1/rho_ds) * d( rho_ds * coef_wpxpyp_implicit * <x'y'>(t+1) )/dz.
  //
  // Note:  When the implicit term is brought over to the left-hand side, the
  //        sign is reversed and the leading "-" in front of all implicit
  //        d[ ] / dz terms is changed to a "+".
  //
  // The timestep index (t+1) means that the value of <x'y'> being used is from
  // the next timestep, which is being advanced to in solving the d<x'y'>/dt
  // equation.
  //
  // When x and y are the same variable, <x'y'> reduces to <x'^2> and <w'x'y'>
  // reduces to <w'x'^2>.  Likewise, when y is set equal to w, <x'y'> becomes
  // <w'x'> and <w'x'y'> reduces to <w'^2 x'>.  The discretization and the code
  // used in this function will be written generally in terms of <x'y'> and
  // coef_wpxpyp_implicit, but also applies to <x'^2> and coef_wpxp2_implicit,
  // as well as to <w'x'> and coef_wp2xp_implicit.
  //
  // The implicit discretization of this term is as follows:
  //
  // 1) Centered Discretization
  //
  // The values of <x'y'> are found on the momentum levels, while the values of
  // coef_wpxpyp_implicit are found on the thermodynamic levels, which is where
  // they were originally calculated by the PDF.  Additionally, the values of
  // rho_ds_zt are found on the thermodynamic levels, and the values of
  // invrs_rho_ds_zm are found on the momentum levels.  The values of <x'y'>
  // are interpolated to the intermediate thermodynamic levels as <x'y'>|_zt.
  // At the thermodynamic levels, the values of coef_wpxpyp_implicit are
  // multiplied by <x'y'>|_zt, and their product is multiplied by rho_ds_zt.
  // Then, the derivative (d/dz) of that expression is taken over the central
  // momentum level, where it is multiplied by -invrs_rho_ds_zm.  This yields
  // the desired result.
  //
  // =xpyp============================================================== m(k+1)
  //
  // -xpyp_zt(interp)-------coef_wpxpyp_implicit---------rho_ds_zt------ t(k)
  //
  // =xpyp=d(rho_ds_zt*coef_wpxpyp_implicit*xpyp_zt)/dz=invrs_rho_ds_zm= m(k)
  //
  // -xpyp_zt(interp)-------coef_wpxpyp_implicit---------rho_ds_zt------ t(k-1)
  //
  // =xpyp============================================================== m(k-1)
  //
  // The vertical indices m(k+1), t(k), m(k), t(k-1), and m(k-1) correspond
  // with altitudes zm(k+1), zt(k), zm(k), zt(k-1), and zm(k-1), respectively.
  // The letter "t" is used for thermodynamic levels and the letter "m" is used
  // for momentum levels.
  //
  // invrs_dzm(k) = 1 / ( zt(k) - zt(k-1) )
  //
  // 2) "Upwind" Discretization.
  //
  // The values of <x'y'> are found on the momentum levels.  The values of
  // coef_wpxpyp_implicit are originally calculated by the PDF on the
  // thermodynamic levels.  They are interpolated to the intermediate momentum
  // levels as coef_wpxpyp_implicit_zm.  Additionally, the values of rho_ds_zm
  // and the values of invrs_rho_ds_zm are found on the momentum levels.  The
  // sign of the turbulent velocity is found on the central momentum level.  At
  // the momentum levels, the values of coef_wpxpyp_implicit_zm are multiplied
  // by <x'y'>, and their product is multiplied by rho_ds_zm.  Then, the
  // derivative (d/dz) of that expression is taken.  When the sign of the
  // turbulent velocity is positive, the "wind" is coming from below, and the
  // derivative involves the central momentum level and the momentum level
  // immediately below it.  When the sign of the turbulent velocity is
  // negative, the "wind" is coming from above, and the derivative involves the
  // central momentum level and the momenum level immediately above it.  After
  // the derivative is taken, it is multiplied by -invrs_rho_ds_zm at the
  // central momentum level.  This yields the desired result.
  //
  // The turbulent velocity for <x'y'> is <w'x'y'> / <x'y'>, which has units of
  // m/s.  The sign of the turbulent velocity is sgn( <w'x'y'> / <x'y'> ),
  // where:
  //
  // sgn(x) = | 1; when x >= 0
  //          | -1; when x < 0.
  //
  // The sign of the turbulent velocity can also be rewritten as
  // sgn( <w'x'y'> ) / sgn( <x'y'> ).  When a variance (<x'^2>) is being solved
  // for, y = x, and sgn( <x'^2> ) is always 1.  The sign of the turbulent
  // velocity reduces to simply sgn( <w'x'^2> ).
  //
  // ---------coef_wpxpyp_implicit-------------------------------------- t(k+1)
  //
  // =xpyp=======coef_wpxpyp_implicit_zm(interp.)=======rho_ds_zm======= m(k+1)
  //
  // ---------coef_wpxpyp_implicit-------------------------------------- t(k)
  //
  // =xpyp===coef_wpxpyp_implicit_zm(interp.)=rho_ds_zm=invrs_rho_ds_zm= m(k)
  //
  // ---------coef_wpxpyp_implicit-------------------------------------- t(k-1)
  //
  // =xpyp=======coef_wpxpyp_implicit_zm(interp.)=======rho_ds_zm======= m(k-1)
  //
  // ---------coef_wpxpyp_implicit-------------------------------------- t(k-2)
  //
  // The vertical indices t(k+1), m(k+1), t(k), m(k), t(k-1), m(k-1), and
  // t(k-2) correspond with altitudes zt(k+1), zm(k+1), zt(k), zm(k), zt(k-1),
  // zm(k-1), and zt(k-2), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) ); and
  // invrs_dzt(k-1) = 1 / ( zm(k) - zm(k-1) ).

  // References:
  //-----------------------------------------------------------------------

  constexpr int ndiags3   = 3;
  constexpr int kp1_mdiag = 0; // Momentum superdiagonal index.
  constexpr int k_mdiag   = 1; // Momentum main diagonal index.
  constexpr int km1_mdiag = 2; // Momentum subdiagonal index.
  constexpr int m_above   = 0; // Upper momentum-level weight category.
  constexpr int m_below   = 1; // Lower momentum-level weight category.
  // Capture component Views/scalars rather than the owning grid object.
  const auto invrs_dzm = gr.invrs_dzm;
  const auto invrs_dzt = gr.invrs_dzt;
  // Weight slots 0:1 correspond to Native categories 1:2.
  const auto weights_zm2zt = gr.weights_zm2zt;
  const real_t grid_dir    = gr.grid_dir;
  using Policy             = Kokkos::RangePolicy<typename Real2D::execution_space>;
  (void)nzt;

  // Both endpoints are independent of the interior; zero them in one launch.
  Kokkos::parallel_for(
      "xpyp_term_ta_pdf_lhs::0", Policy(0, ngrdcol * ndiags3), KOKKOS_LAMBDA(int cell) {
        const int i           = cell / ndiags3;
        const int b           = cell % ndiags3;
        lhs_ta(b, i, 0)       = zero;
        lhs_ta(b, i, nzm - 1) = zero;
      });
  if (!l_upwind_xpyp_turbulent_adv) {
    // Centered discretization.
    Kokkos::parallel_for(
        "xpyp_term_ta_pdf_lhs::1", Policy(0, (nzm - 2) * ngrdcol), KOKKOS_LAMBDA(int cell) {
          const int k = cell / ngrdcol + 1;
          const int i = cell % ngrdcol;
          // Momentum superdiagonal: [ x xpyp(k+1,<t+1>) ]
          lhs_ta(kp1_mdiag, i, k) = invrs_rho_ds_zm(i, k) * invrs_dzm(i, k) * rho_ds_zt(i, k) *
                                    coef_wpxpyp_implicit(i, k) * weights_zm2zt(i, k, m_above);
          // Momentum main diagonal: [ x xpyp(k,<t+1>) ]
          lhs_ta(k_mdiag, i, k) =
              invrs_rho_ds_zm(i, k) * invrs_dzm(i, k) *
              (rho_ds_zt(i, k) * coef_wpxpyp_implicit(i, k) * weights_zm2zt(i, k, m_below) -
               rho_ds_zt(i, k - 1) * coef_wpxpyp_implicit(i, k - 1) *
                   weights_zm2zt(i, k - 1, m_above));
          // Momentum subdiagonal: [ x xpyp(k-1,<t+1>) ]
          lhs_ta(km1_mdiag, i, k) = -invrs_rho_ds_zm(i, k) * invrs_dzm(i, k) * rho_ds_zt(i, k - 1) *
                                    coef_wpxpyp_implicit(i, k - 1) *
                                    weights_zm2zt(i, k - 1, m_below);
        });
  } else { // l_upwind_xpyp_turbulent_adv
    // "Upwind" discretization
    Kokkos::parallel_for(
        "xpyp_term_ta_pdf_lhs::2", Policy(0, (nzm - 2) * ngrdcol), KOKKOS_LAMBDA(int cell) {
          const int k = cell / ngrdcol + 1;
          const int i = cell % ngrdcol;
          if (grid_dir * sgn_turbulent_vel(i, k) > zero) {
            // The "wind" is blowing upward for an ascending grid or the "wind"
            // is blowing downward for a descending grid.
            // Momentum superdiagonal: [ x xpyp(k+1,<t+1>) ]
            lhs_ta(kp1_mdiag, i, k) = zero;
            // Momentum main diagonal: [ x xpyp(k,<t+1>) ]
            lhs_ta(k_mdiag, i, k) = invrs_rho_ds_zm(i, k) * invrs_dzt(i, k - 1) * rho_ds_zm(i, k) *
                                    coef_wpxpyp_implicit_zm(i, k);
            // Momentum subdiagonal: [ x xpyp(k-1,<t+1>) ]
            lhs_ta(km1_mdiag, i, k) = -invrs_rho_ds_zm(i, k) * invrs_dzt(i, k - 1) *
                                      rho_ds_zm(i, k - 1) * coef_wpxpyp_implicit_zm(i, k - 1);
          } else { // Remaining turbulent directions, including zero.
            // The "wind" is blowing downward for an ascending grid or the "wind"
            // is blowing upward for a descending grid.
            // Momentum superdiagonal: [ x xpyp(k+1,<t+1>) ]
            lhs_ta(kp1_mdiag, i, k) = invrs_rho_ds_zm(i, k) * invrs_dzt(i, k) *
                                      rho_ds_zm(i, k + 1) * coef_wpxpyp_implicit_zm(i, k + 1);
            // Momentum main diagonal: [ x xpyp(k,<t+1>) ]
            lhs_ta(k_mdiag, i, k) = -invrs_rho_ds_zm(i, k) * invrs_dzt(i, k) * rho_ds_zm(i, k) *
                                    coef_wpxpyp_implicit_zm(i, k);
            // Momentum subdiagonal: [ x xpyp(k-1,<t+1>) ]
            lhs_ta(km1_mdiag, i, k) = zero;
          } // sgn_turbulent_vel
        });
  }
}

//=============================================================================
// coef_wpxpyp_implicit: Coef. of <x'y'> in <w'x'y'>; t-lev [m/s]
// invrs_rho_ds_zm: Inv dry, static density @ m-level [m^3/kg]
// rho_ds_zm: Dry, static density at m-lev [kg/m^3]
// lhs_ta storage: (band,column,momentum level), super/main/sub bands [1/s].
void xpyp_term_ta_pdf_lhs_godunov(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                                  Real2D coef_wpxpyp_implicit,                   // In
                                  Real2D invrs_rho_ds_zm, Real2D rho_ds_zm,      // In
                                  Real3D lhs_ta)                                 // Out
{
  // Description:
  //   This subroutine is a revised version of xpyp_term_ta_pdf_lhs. The
  //   revisions are made to use the  Godunov-like upwind scheme for the
  //   vertical discretization of the turbulent advection term. This subroutine
  //   fills one rank-3 (band,column,momentum level) operator array, with
  //   super/main/sub bands and zero values at both storage endpoints.
  //
  // All arguments are required; storage interior levels are 1 through nzm-2.
  //
  // Notes:
  //   This subroutine exists for testing of Godunov-like upwind scheme.
  //   THIS SUBROUTINE DOES NOT HANDLE BOUNDARY CONDITIONS AND SETS THEM TO 0
  //---------------------------------------------------------------------------------------------

  constexpr int ndiags3   = 3;
  constexpr int kp1_mdiag = 0; // Momentum superdiagonal index.
  constexpr int k_mdiag   = 1; // Momentum main diagonal index.
  constexpr int km1_mdiag = 2; // Momentum subdiagonal index.
  // Capture component Views/scalars rather than the owning grid object.
  const auto invrs_dzm  = gr.invrs_dzm;
  const real_t grid_dir = gr.grid_dir;
  using Policy          = Kokkos::RangePolicy<typename Real2D::execution_space>;
  (void)nzt;

  // Both endpoints are independent of the interior; zero them in one launch.
  Kokkos::parallel_for(
      "xpyp_term_ta_pdf_lhs_godunov::0", Policy(0, ngrdcol * ndiags3), KOKKOS_LAMBDA(int cell) {
        const int i           = cell / ndiags3;
        const int b           = cell % ndiags3;
        lhs_ta(b, i, 0)       = real_t(0.0);
        lhs_ta(b, i, nzm - 1) = real_t(0.0);
      });
  // Godunov-like upwind discretization
  Kokkos::parallel_for(
      "xpyp_term_ta_pdf_lhs_godunov::1", Policy(0, (nzm - 2) * ngrdcol), KOKKOS_LAMBDA(int cell) {
        const int k = cell / ngrdcol + 1;
        const int i = cell % ngrdcol;
        // Momentum superdiagonal: [ x xpyp(k+1,<t+1>) ]
        lhs_ta(kp1_mdiag, i, k) = invrs_rho_ds_zm(i, k) * invrs_dzm(i, k) * rho_ds_zm(i, k + 1) *
                                  grid_dir *
                                  Kokkos::min(real_t(0.0), grid_dir * coef_wpxpyp_implicit(i, k));
        // Momentum main diagonal: [ x xpyp(k,<t+1>) ]
        lhs_ta(k_mdiag, i, k) =
            invrs_rho_ds_zm(i, k) * invrs_dzm(i, k) * rho_ds_zm(i, k) * grid_dir *
            (Kokkos::max(real_t(0.0), grid_dir * coef_wpxpyp_implicit(i, k)) -
             Kokkos::min(real_t(0.0), grid_dir * coef_wpxpyp_implicit(i, k - 1)));
        // Momentum subdiagonal: [ x xpyp(k-1,<t+1>) ]
        lhs_ta(km1_mdiag, i, k) =
            -invrs_rho_ds_zm(i, k) * invrs_dzm(i, k) * rho_ds_zm(i, k - 1) * grid_dir *
            Kokkos::max(real_t(0.0), grid_dir * coef_wpxpyp_implicit(i, k - 1));
      });
}

//=============================================================================
// term_wpxpyp_explicit: RHS: <w'x'y'> eq; t-levs      [m/s(x un)(y un)]
// rho_ds_zt: Dry, static density at t-levs          [kg/m^3]
// rho_ds_zm: Dry, static density at m-levs       [kg/m^3]
// invrs_rho_ds_zm: Inv dry, static density at m-levs   [m^3/kg]
// l_upwind_xpyp_turbulent_adv: Flag to use "upwind" discretization
// sgn_turbulent_vel: Sign of the turbulent velocity           [-]
// term_wpxpyp_explicit_zm: term_wpxpyp_expl. zm       [m/s(x un)(y un)]
// rhs_ta: RHS portion of xpyp turbulent advection      [(x un)(y un)/s]
void xpyp_term_ta_pdf_rhs(int nzm, int nzt, int ngrdcol, const grid& gr,
                          Real2D term_wpxpyp_explicit,        // In
                          Real2D rho_ds_zt, Real2D rho_ds_zm, // In
                          Real2D invrs_rho_ds_zm,             // In
                          bool l_upwind_xpyp_turbulent_adv,   // In
                          Real2D sgn_turbulent_vel,           // In
                          Real2D term_wpxpyp_explicit_zm,     // In
                          Real2D rhs_ta)                      // Out
{
  // Description:
  // Turbulent advection of <w'x'>, <x'^2>, and <x'y'>:  explicit portion of
  // the code.
  //
  // 1) <w'x'>
  //
  // The d<w'x'>/dt equation contains a turbulent advection term:
  //
  // - (1/rho_ds) * d( rho_ds * <w'^2 x'> )/dz.
  //
  // The value of <w'^2 x'> is found by integrating over the multivariate PDF
  // of w and x, as detailed in function calc_wp2xp_pdf, which is found in
  // module pdf_closure_module in pdf_closure_module.F90.
  //
  // The equation obtained for <w'^2 x'> is written in terms of PDF parameters.
  // Substitutions that are specific to the type of PDF used are made for the
  // PDF parameters in order to write the <w'x'> turbulent advection term in
  // the following form:
  //
  // <w'^2 x'> = coef_wp2xp_implicit * <w'x'> + term_wp2xp_explicit.
  //
  // For the ADG1 PDF, coef_wp2xp_implicit is a_1 * < w'^3 > / < w'^2 >, where
  // a_1 is 1 / ( 1 - sigma_sqd_w ).  The value of term_wp2xp_explicit is 0, as
  // the <w'x'> turbulent advection term is entirely implicit.
  //
  // For the new PDF, the calculations of both coef_wp2xp_implicit and
  // term_wp2xp_explicit are detailed in function calc_coefs_wp2xp_semiimpl,
  // which is found in module new_pdf in new_pdf.F90.
  //
  // For the new hybrid PDF, the calculation of coef_wp2xp_implicit is
  // detailed in function calc_coef_wp2xp_implicit, which is found in module
  // new_hybrid_pdf in new_hybrid_pdf.F90.  The value of term_wp2xp_explicit
  // is 0, as the <w'x'> turbulent advection term is entirely implicit.
  //
  // For explicit turbulent advection, the value of coef_wp2xp_implicit is 0
  // and the value of term_wp2xp_explicit is <w'^2 x'>, as calculated by
  // retaining the equation for <w'^2 x'> that is written in terms of PDF
  // parameters.  This is a general form that can be used with any type of PDF.
  //
  // The <w'x'> turbulent advection term is rewritten as:
  //
  // - (1/rho_ds)
  //   * d( rho_ds * ( coef_wp2xp_implicit * <w'x'> + term_wp2xp_explicit ) )
  //     /dz.
  //
  // The variable <w'x'> is evaluated at the (t+1) timestep, which allows the
  // <w'x'> turbulent advection term to be expressed semi-implicitly as:
  //
  // - (1/rho_ds) * d( rho_ds * coef_wp2xp_implicit * <w'x'>(t+1) )/dz
  // - (1/rho_ds) * d( rho_ds * term_wp2xp_explicit )/dz.
  //
  // The explicit portion of <w'x'> turbulent advection term is:
  //
  // - (1/rho_ds) * d( rho_ds * term_wp2xp_explicit )/dz.
  //
  // 2) <x'^2>
  //
  // The d<x'^2>/dt equation contains a turbulent advection term:
  //
  // - (1/rho_ds) * d( rho_ds * <w'x'^2> )/dz;
  //
  // The value of <w'x'^2> is found by integrating over the multivariate PDF of
  // w and x, as detailed in function calc_wpxp2_pdf, which is found in module
  // pdf_closure_module in pdf_closure_module.F90.
  //
  // The equation obtained for <w'x'^2> is written in terms of PDF parameters.
  // Substitutions that are specific to the type of PDF used are made for the
  // PDF parameters in order to write the <x'^2> turbulent advection term in
  // the following form:
  //
  // <w'x'^2> = coef_wpxp2_implicit * <x'^2> + term_wpxp2_explicit.
  //
  // For the ADG1 PDF, the value of coef_wpxp2_implicit is
  // (1/3)*beta * a_1 * < w'^3 > / < w'^2 >.  The value of term_wpxp2_explicit
  // is (1-(1/3)*beta) * a_1^2 * < w'x' >^2 * < w'^3 > / < w'^2 >^2, where
  // a_1 is 1 / ( 1 - sigma_sqd_w ).
  //
  // For the new PDF, the calculation of coef_wpxp2_implicit is detailed in
  // function calc_coef_wpxp2_implicit, which is found in module new_pdf in
  // new_pdf.F90.  The value of term_wpxp2_explicit is 0, as the <x'^2>
  // turbulent advection term is entirely implicit.
  //
  // For the new hybrid PDF, the calculation of both coef_wpxp2_implicit and
  // term_wpxp2_explicit are detailed in subroutine calc_coefs_wpxp2_semiimpl,
  // which is found in module new_hybrid_pdf in new_hybrid_pdf.F90.
  //
  // For explicit turbulent advection, the value of coef_wpxp2_implicit is 0
  // and the value of term_wpxp2_explicit is <w'x'^2>, as calculated by
  // retaining the equation for <w'x'^2> that is written in terms of PDF
  // parameters.  This is a general form that can be used with any type of PDF.
  //
  // The <x'^2> turbulent advection term is rewritten as:
  //
  // - (1/rho_ds)
  //   * d( rho_ds * ( coef_wpxp2_implicit * <x'^2> + term_wpxp2_explicit ) )
  //     /dz;
  //
  // The variable <x'^2> is evaluated at the (t+1) timestep, which allows the
  // <x'^2> turbulent advection term to be expressed semi-implicitly as:
  //
  // - (1/rho_ds) * d( rho_ds * coef_wpxp2_implicit * <x'^2>(t+1) )/dz;
  // - (1/rho_ds) * d( rho_ds * term_wpxp2_explicit )/dz.
  //
  // The explicit portion of <x'^2> turbulent advection term is:
  //
  // - (1/rho_ds) * d( rho_ds * term_wpxp2_explicit )/dz.
  //
  // 3) <x'y'>
  //
  // The d<x'y'>/dt equation contains a turbulent advection term:
  //
  // - (1/rho_ds) * d( rho_ds * <w'x'y'> )/dz.
  //
  // The value of <w'x'y'> is found by integrating over the multivariate PDF of
  // w, x, and y, as detailed in function calc_wpxpyp_pdf, which is found in
  // module pdf_closure_module in pdf_closure_module.F90.
  //
  // The equation obtained for <w'x'y'> is written in terms of PDF parameters.
  // Substitutions that are specific to the type of PDF used are made for the
  // PDF parameters in order to write the <x'y'> turbulent advection term in
  // the following form:
  //
  // <w'x'y'> = coef_wpxpyp_implicit * <x'y'> + term_wpxpyp_explicit.
  //
  // For the ADG1 PDF, the value of coef_wpxpyp_implicit is
  // (1/3)*beta * a_1 * < w'^3 > / < w'^2 >.  The value of term_wpxpyp_explicit
  // is (1-(1/3)*beta) * a_1^2 * < w'x' > * < w'y' > * < w'^3 > / < w'^2 >^2,
  // where a_1 is 1 / ( 1 - sigma_sqd_w ).
  //
  // For the new PDF, the calculation of both coef_wpxpyp_implicit and
  // term_wpxpyp_explicit are detailed in function calc_coefs_wpxpyp_semiimpl,
  // which is found in module new_pdf in new_pdf.F90.
  //
  // For the new hybrid PDF, the calculation of both coef_wpxpyp_implicit
  // and term_wpxpyp_explicit are detailed in subroutine
  // calc_coefs_wpxpyp_semiimpl, which is found in module new_hybrid_pdf in
  // new_hybrid_pdf.F90.
  //
  // For explicit turbulent advection, the value of coef_wpxpyp_implicit is 0
  // and the value of term_wpxpyp_explicit is <w'x'y'>, as calculated by
  // retaining the equation for <w'x'y'> that is written in terms of PDF
  // parameters.  This is a general form that can be used with any type of PDF.
  //
  // The <x'y'> turbulent advection term is rewritten as:
  //
  // - (1/rho_ds)
  //   * d( rho_ds * ( coef_wpxpyp_implicit * <x'y'> + term_wpxpyp_explicit ) )
  //     /dz.
  //
  // The variable <x'y'> is evaluated at the (t+1) timestep, which allows the
  // <x'y'> turbulent advection term to be expressed semi-implicitly as:
  //
  // - (1/rho_ds) * d( rho_ds * coef_wpxpyp_implicit * <x'y'>(t+1) )/dz
  // - (1/rho_ds) * d( rho_ds * term_wpxpyp_explicit )/dz.
  //
  // The explicit portion of <x'y'> turbulent advection term is:
  //
  // - (1/rho_ds) * d( rho_ds * term_wpxpyp_explicit )/dz.
  //
  // When x and y are the same variable, <x'y'> reduces to <x'^2> and <w'x'y'>
  // reduces to <w'x'^2>.  Likewise, when y is set equal to w, <x'y'> becomes
  // <w'x'> and <w'x'y'> reduces to <w'^2 x'>.  The discretization and the code
  // used in this function will be written generally in terms of <x'y'> and
  // term_wpxpyp_explicit, but also applies to <x'^2> and term_wpxp2_explicit,
  // as well as to <w'x'> and term_wp2xp_explicit.
  //
  // The explicit discretization of this term is as follows:
  //
  // 1) Centered Discretization
  //
  // The values of <x'y'> are found on the momentum levels, while the values of
  // term_wpxpyp_explicit are found on the thermodynamic levels, which is where
  // they were originally calculated by the PDF.  Additionally, the values of
  // rho_ds_zt are found on the thermodynamic levels, and the values of
  // invrs_rho_ds_zm are found on the momentum levels.  At the thermodynamic
  // levels, the values of term_wpxpyp_explicit are multiplied by rho_ds_zt.
  // Then, the derivative (d/dz) of that expression is taken over the central
  // momentum level, where it is multiplied by -invrs_rho_ds_zm.  This yields
  // the desired result.
  //
  // ---------term_wpxpyp_explicitp1-------rho_ds_ztp1------------------ t(k)
  //
  // ==xpyp==d( rho_ds_zt * term_wpxpyp_explicit )/dz==invrs_rho_ds_zm== m(k)
  //
  // ---------term_wpxpyp_explicit---------rho_ds_zt-------------------- t(k-1)
  //
  // The vertical indices t(k), m(k), and t(k-1) correspond with altitudes
  // zt(k), zm(k), and zt(k-1), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // invrs_dzm(k) = 1 / ( zt(k) - zt(k-1) )
  //
  // 2) "Upwind" Discretization.
  //
  // The values of <x'y'> are found on the momentum levels.  The values of
  // term_wpxpyp_explicit are originally calculated by the PDF on the
  // thermodynamic levels.  They are interpolated to the intermediate momentum
  // levels as term_wpxpyp_explicit_zm.  Additionally, the values of rho_ds_zm
  // and the values of invrs_rho_ds_zm are found on the momentum levels.  The
  // sign of the turbulent velocity is found on the central momentum level.  At
  // the momentum levels, the values of term_wpxpyp_explicit_zm are multiplied
  // by rho_ds_zm.  Then, the derivative (d/dz) of that expression is taken.
  // When the sign of the turbulent velocity is positive, the "wind" is coming
  // from below, and the derivative involves the central momentum level and the
  // momentum level immediately below it.  When the sign of the turbulent
  // velocity is negative, the "wind" is coming from above, and the derivative
  // involves the central momentum level and the momenum level immediately
  // above it.  After the derivative is taken, it is multiplied by
  // -invrs_rho_ds_zm at the central momentum level.  This yields the desired
  // result.
  //
  // The turbulent velocity for <x'y'> is <w'x'y'> / <x'y'>, which has units of
  // m/s.  The sign of the turbulent velocity is sgn( <w'x'y'> / <x'y'> ),
  // where:
  //
  // sgn(x) = | 1; when x >= 0
  //          | -1; when x < 0.
  //
  // The sign of the turbulent velocity can also be rewritten as
  // sgn( <w'x'y'> ) / sgn( <x'y'> ).  When a variance (<x'^2>) is being solved
  // for, y = x, and sgn( <x'^2> ) is always 1.  The sign of the turbulent
  // velocity reduces to simply sgn( <w'x'^2> ).
  //
  // ---------term_wpxpyp_explicit-------------------------------------- t(k+1)
  //
  // ========term_wpxpyp_explicit_zm(interp.)=======rho_ds_zm=========== m(k+1)
  //
  // ---------term_wpxpyp_explicit-------------------------------------- t(k)
  //
  // =xpyp===term_wpxpyp_explicit_zm(interp.)=rho_ds_zm=invrs_rho_ds_zm= m(k)
  //
  // ---------term_wpxpyp_explicit-------------------------------------- t(k-1)
  //
  // ========term_wpxpyp_explicit_zm(interp.)=======rho_ds_zm=========== m(k-1)
  //
  // ---------term_wpxpyp_explicit-------------------------------------- t(k-2)
  //
  // The vertical indices t(k+1), m(k+1), t(k), m(k), t(k-1), m(k-1), and
  // t(k-2) correspond with altitudes zt(k+1), zm(k+1), zt(k), zm(k), zt(k-1),
  // zm(k-1), and zt(k-2), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) ); and
  // invrs_dzt(k-1) = 1 / ( zm(k) - zm(k-1) ).

  // References:
  //-----------------------------------------------------------------------

  // Capture component Views/scalars rather than the owning grid object.
  const auto invrs_dzm  = gr.invrs_dzm;
  const auto invrs_dzt  = gr.invrs_dzt;
  const real_t grid_dir = gr.grid_dir;
  using Policy          = Kokkos::RangePolicy<typename Real2D::execution_space>;
  (void)nzt;

  // Both endpoints are independent of the interior; zero them in one launch.
  Kokkos::parallel_for(
      "xpyp_term_ta_pdf_rhs::0", Policy(0, ngrdcol), KOKKOS_LAMBDA(int cell) {
        const int i        = cell;
        rhs_ta(i, 0)       = zero;
        rhs_ta(i, nzm - 1) = zero;
      });
  if (!l_upwind_xpyp_turbulent_adv) {
    // Centered discretization.
    Kokkos::parallel_for(
        "xpyp_term_ta_pdf_rhs::1", Policy(0, (nzm - 2) * ngrdcol), KOKKOS_LAMBDA(int cell) {
          const int k  = cell / ngrdcol + 1;
          const int i  = cell % ngrdcol;
          rhs_ta(i, k) = -invrs_rho_ds_zm(i, k) * invrs_dzm(i, k) *
                         (rho_ds_zt(i, k) * term_wpxpyp_explicit(i, k) -
                          rho_ds_zt(i, k - 1) * term_wpxpyp_explicit(i, k - 1));
        });
  } else { // l_upwind_xpyp_turbulent_adv
    // "Upwind" discretization
    Kokkos::parallel_for(
        "xpyp_term_ta_pdf_rhs::2", Policy(0, (nzm - 2) * ngrdcol), KOKKOS_LAMBDA(int cell) {
          const int k = cell / ngrdcol + 1;
          const int i = cell % ngrdcol;
          if (grid_dir * sgn_turbulent_vel(i, k) > zero) {
            // The "wind" is blowing upward for an ascending grid or the "wind"
            // is blowing downward for a descending grid.
            rhs_ta(i, k) = -invrs_rho_ds_zm(i, k) * invrs_dzt(i, k - 1) *
                           (rho_ds_zm(i, k) * term_wpxpyp_explicit_zm(i, k) -
                            rho_ds_zm(i, k - 1) * term_wpxpyp_explicit_zm(i, k - 1));
          } else { // Remaining turbulent directions, including zero.
            // The "wind" is blowing downward for an ascending grid or the
            // "wind" is blowing upward for a descending grid.
            rhs_ta(i, k) = -invrs_rho_ds_zm(i, k) * invrs_dzt(i, k) *
                           (rho_ds_zm(i, k + 1) * term_wpxpyp_explicit_zm(i, k + 1) -
                            rho_ds_zm(i, k) * term_wpxpyp_explicit_zm(i, k));
          } // sgn_turbulent_vel
        });
  }
}

//=============================================================================
// term_wpxpyp_explicit_zm: RHS: <w'x'y'> eq; m-lev(k)   [m/s(x un)(y un)]
// invrs_rho_ds_zm: Inv dry, static density at m-lev (k) [m^3/kg]
// sgn_turbulent_vel: Sign of turbulent velocity on thermodynamic levels here [-].
// rho_ds_zm: Dry, static density at m-lev (k) [kg/m^3]
// rhs_ta: Explicit turbulent-advection tendency [(x units)(y units)/s].
void xpyp_term_ta_pdf_rhs_godunov(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                                  Real2D term_wpxpyp_explicit_zm,                // In
                                  Real2D invrs_rho_ds_zm,                        // In
                                  Real2D sgn_turbulent_vel,                      // In
                                  Real2D rho_ds_zm,                              // In
                                  Real2D rhs_ta)                                 // Out
{
  // Description:
  //   This subroutine intends to add godunov upwind difference scheme based
  //   on xpyp_term_ta_pdf_rhs.  The revisions are made to use the Godunov-like
  //   upwind scheme for the vertical discretization.
  //   This subroutine returns an array of values for every grid level.
  //
  // All arguments are required; storage interior levels are 1 through nzm-2.
  //
  // Notes:
  //   This subroutine exists for testing of Godunov-like upwind scheme.
  //   THIS SUBROUTINE DOES NOT HANDLE BOUNDARY CONDITIONS AND SETS THEM TO 0
  //--------------------------------------------------------------------------------------------

  // Capture component Views/scalars rather than the owning grid object.
  const auto invrs_dzm  = gr.invrs_dzm;
  const real_t grid_dir = gr.grid_dir;
  using Policy          = Kokkos::RangePolicy<typename Real2D::execution_space>;
  (void)nzt;

  // Both endpoints are independent of the interior; zero them in one launch.
  Kokkos::parallel_for(
      "xpyp_term_ta_pdf_rhs_godunov::0", Policy(0, ngrdcol), KOKKOS_LAMBDA(int cell) {
        const int i        = cell;
        rhs_ta(i, 0)       = real_t(0.0);
        rhs_ta(i, nzm - 1) = real_t(0.0);
      });
  Kokkos::parallel_for(
      "xpyp_term_ta_pdf_rhs_godunov::1", Policy(0, (nzm - 2) * ngrdcol), KOKKOS_LAMBDA(int cell) {
        const int k  = cell / ngrdcol + 1;
        const int i  = cell % ngrdcol;
        rhs_ta(i, k) = -invrs_rho_ds_zm(i, k) * invrs_dzm(i, k) * grid_dir *
                       (Kokkos::min(real_t(0.0), grid_dir * sgn_turbulent_vel(i, k)) *
                            rho_ds_zm(i, k + 1) * term_wpxpyp_explicit_zm(i, k + 1) +
                        Kokkos::max(real_t(0.0), grid_dir * sgn_turbulent_vel(i, k)) *
                            rho_ds_zm(i, k) * term_wpxpyp_explicit_zm(i, k) -
                        Kokkos::min(real_t(0.0), grid_dir * sgn_turbulent_vel(i, k - 1)) *
                            rho_ds_zm(i, k) * term_wpxpyp_explicit_zm(i, k) -
                        Kokkos::max(real_t(0.0), grid_dir * sgn_turbulent_vel(i, k - 1)) *
                            rho_ds_zm(i, k - 1) * term_wpxpyp_explicit_zm(i, k - 1));
      });
}

} // namespace clubb_kokkos
