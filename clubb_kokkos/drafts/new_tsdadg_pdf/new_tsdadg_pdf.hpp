// Close the trivariate skewness-dependent analytic double Gaussian PDF.
//
// tsdadg_pdf_driver chooses w, total water rt or liquid-water potential temperature
// thl as the setter using the largest absolute skewness, then closes the other
// variables with the same mixture fraction. calc_setter_parameters and
// calc_L_x_Skx_fnc expose the Native scalar setter and mean-spread equations.
// Component means retain input units; component variances have squared units.
//
// The public driver retains its Native nz-profile signature and accepts existing
// rank-one column subviews. Fields retain their column-first owner; local scratch
// stays with this call. Scalar helpers are device-callable; warnings print on the
// host in Native level/field order after an integer-only warning count. There is
// no persistent state or Fortran call. Compilation/numerical evidence is pending.
// Source: src/CLUBB_core/new_tsdadg_pdf.F90, including the full derivations below.
#pragma once

#include "array_types.hpp"
#include "constants_clubb.hpp"
#include <Kokkos_MathematicalFunctions.hpp>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

namespace clubb_kokkos {

// Forward declarations retain the Native definition order and public/private API.
// The named detail scope matches the existing saturation/new-PDF header pattern
// for a Native private device helper, without header-local function identities.
KOKKOS_INLINE_FUNCTION void calc_setter_parameters(real_t xm, real_t xp2, real_t Skx,
                                                   real_t sgn_wpxp,                    // In
                                                   real_t big_L_x_1, real_t big_L_x_2, // In
                                                   real_t& mu_x_1, real_t& mu_x_2,
                                                   real_t& sigma_x_1_sqd,                    // Out
                                                   real_t& sigma_x_2_sqd, real_t& mixt_frac, // Out
                                                   real_t& coef_sigma_x_1_sqd,               // Out
                                                   real_t& coef_sigma_x_2_sqd                // Out
);
KOKKOS_INLINE_FUNCTION void calc_L_x_Skx_fnc(real_t Skx, real_t sgn_wpxp,            // In
                                             real_t small_l_x_1, real_t small_l_x_2, // In
                                             real_t& big_L_x_1, real_t& big_L_x_2    // Out
);
namespace new_tsdadg_pdf_detail {
KOKKOS_INLINE_FUNCTION void calc_respnder_parameters(real_t xm, real_t xp2, real_t Skx,
                                                     real_t sgn_wpxp,                    // In
                                                     real_t mixt_frac, real_t big_L_x_1, // In
                                                     real_t& mu_x_1, real_t& mu_x_2,     // Out
                                                     real_t& sigma_x_1_sqd,
                                                     real_t& sigma_x_2_sqd,      // Out
                                                     real_t& coef_sigma_x_1_sqd, // Out
                                                     real_t& coef_sigma_x_2_sqd  // Out
);
} // namespace new_tsdadg_pdf_detail

// =============================================================================
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23>
void tsdadg_pdf_driver(int nz, Array0 wm, Array1 rtm, Array2 thlm, Array3 wp2, Array4 rtp2,
                       Array5 thlp2,                                                        // In
                       Array6 Skw, Array7 Skrt, Array8 Skthl, Array9 wprtp, Array10 wpthlp, // In
                       Array11 mu_w_1, Array12 mu_w_2,                                      // Out
                       Array13 mu_rt_1, Array14 mu_rt_2,                                    // Out
                       Array15 mu_thl_1, Array16 mu_thl_2,                                  // Out
                       Array17 sigma_w_1_sqd, Array18 sigma_w_2_sqd,                        // Out
                       Array19 sigma_rt_1_sqd, Array20 sigma_rt_2_sqd,                      // Out
                       Array21 sigma_thl_1_sqd, Array22 sigma_thl_2_sqd,                    // Out
                       Array23 mixt_frac                                                    // Out
) {
  // Description:
  // Selects which variable is used to set the mixture fraction for the PDF
  // ("the setter") and which variables are handled after that mixture fraction
  // has been set ("the responders").  Traditionally, w has been used to set
  // the PDF.  However, here, the variable with the greatest magnitude of
  // skewness is used to set the PDF.
  //
  // References:
  //-----------------------------------------------------------------------
  // Profile means/variances use [m/s]/[m^2/s^2] for w, [kg/kg]/[(kg/kg)^2]
  // for rt and [K]/[K^2] for thl; all skewness and spread coefficients are [-].
  static_assert(
      Array0::rank == 1 && Array1::rank == 1 && Array2::rank == 1 && Array3::rank == 1 &&
      Array4::rank == 1 && Array5::rank == 1 && Array6::rank == 1 && Array7::rank == 1 &&
      Array8::rank == 1 && Array9::rank == 1 && Array10::rank == 1 && Array11::rank == 1 &&
      Array12::rank == 1 && Array13::rank == 1 && Array14::rank == 1 && Array15::rank == 1 &&
      Array16::rank == 1 && Array17::rank == 1 && Array18::rank == 1 && Array19::rank == 1 &&
      Array20::rank == 1 && Array21::rank == 1 && Array22::rank == 1 && Array23::rank == 1);
  using Execution = typename Array0::execution_space;
  // The Native per-level arrays below are consumed only at the same level.
  // Retain their names as distinct work-item scalars in the one level kernel.
  const real_t small_l_w_1 = real_t(.75), small_l_w_2 = real_t(.5);
  const real_t small_l_rt_1 = real_t(.75), small_l_rt_2 = real_t(.5);
  const real_t small_l_thl_1 = real_t(.75), small_l_thl_2 = real_t(.5);
  Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Array0::memory_space> negative_variances(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "tsdadg_pdf_driver::negative_variances"), nz,
      6);
  int warning_count = 0;
  Kokkos::parallel_reduce(
      "tsdadg_pdf_driver", Kokkos::RangePolicy<Execution>(0, nz),
      KOKKOS_LAMBDA(const int k, int& count) {
        real_t big_L_w_1, big_L_w_2, big_L_rt_1, big_L_rt_2, big_L_thl_1, big_L_thl_2;
        real_t sgn_wprtp, sgn_wpthlp, sgn_wp2;
        real_t coef_sigma_w_1_sqd, coef_sigma_w_2_sqd;
        real_t coef_sigma_rt_1_sqd, coef_sigma_rt_2_sqd;
        real_t coef_sigma_thl_1_sqd, coef_sigma_thl_2_sqd;
        // Calculate sgn( <w'rt'> ). Zero takes the Native nonnegative branch.
        if (wprtp(k) >= zero) {
          sgn_wprtp = one;
        } else {
          sgn_wprtp = -one;
        }
        // Calculate sgn( <w'thl'> ).
        if (wpthlp(k) >= zero) {
          sgn_wpthlp = one;
        } else {
          sgn_wpthlp = -one;
        }
        // The sign of the variance of w is always positive.
        sgn_wp2 = one;
        calc_L_x_Skx_fnc(Skw(k), sgn_wp2,          // In
                         small_l_w_1, small_l_w_2, // In
                         big_L_w_1, big_L_w_2);    // Out

        calc_L_x_Skx_fnc(Skrt(k), sgn_wprtp,         // In
                         small_l_rt_1, small_l_rt_2, // In
                         big_L_rt_1, big_L_rt_2);    // Out

        calc_L_x_Skx_fnc(Skthl(k), sgn_wpthlp,         // In
                         small_l_thl_1, small_l_thl_2, // In
                         big_L_thl_1, big_L_thl_2);    // Out

        // The variable with the greatest magnitude of skewness will be the setter
        // variable and the other variables will be responder variables.
        if (Kokkos::abs(Skw(k)) >= Kokkos::abs(Skrt(k)) &&
            Kokkos::abs(Skw(k)) >= Kokkos::abs(Skthl(k))) {

          // The variable w has the greatest magnitude of skewness.

          calc_setter_parameters(wm(k), wp2(k),                  // In
                                 Skw(k), sgn_wp2,                // In
                                 big_L_w_1, big_L_w_2,           // In
                                 mu_w_1(k), mu_w_2(k),           // Out
                                 sigma_w_1_sqd(k),               // Out
                                 sigma_w_2_sqd(k), mixt_frac(k), // Out
                                 coef_sigma_w_1_sqd,             // Out
                                 coef_sigma_w_2_sqd);            // Out

          new_tsdadg_pdf_detail::calc_respnder_parameters(rtm(k), rtp2(k),          // In
                                                          Skrt(k), sgn_wprtp,       // In
                                                          mixt_frac(k), big_L_rt_1, // In
                                                          mu_rt_1(k), mu_rt_2(k),   // Out
                                                          sigma_rt_1_sqd(k),        // Out
                                                          sigma_rt_2_sqd(k),        // Out
                                                          coef_sigma_rt_1_sqd,      // Out
                                                          coef_sigma_rt_2_sqd);     // Out

          new_tsdadg_pdf_detail::calc_respnder_parameters(thlm(k), thlp2(k),         // In
                                                          Skthl(k), sgn_wpthlp,      // In
                                                          mixt_frac(k), big_L_thl_1, // In
                                                          mu_thl_1(k), mu_thl_2(k),  // Out
                                                          sigma_thl_1_sqd(k),        // Out
                                                          sigma_thl_2_sqd(k),        // Out
                                                          coef_sigma_thl_1_sqd,      // Out
                                                          coef_sigma_thl_2_sqd);     // Out

        } else if (Kokkos::abs(Skrt(k)) > Kokkos::abs(Skw(k)) &&
                   Kokkos::abs(Skrt(k)) >= Kokkos::abs(Skthl(k))) {

          // The variable rt has the greatest magnitude of skewness.

          calc_setter_parameters(rtm(k), rtp2(k),                 // In
                                 Skrt(k), sgn_wprtp,              // In
                                 big_L_rt_1, big_L_rt_2,          // In
                                 mu_rt_1(k), mu_rt_2(k),          // Out
                                 sigma_rt_1_sqd(k),               // Out
                                 sigma_rt_2_sqd(k), mixt_frac(k), // Out
                                 coef_sigma_rt_1_sqd,             // Out
                                 coef_sigma_rt_2_sqd);            // Out

          new_tsdadg_pdf_detail::calc_respnder_parameters(wm(k), wp2(k),           // In
                                                          Skw(k), sgn_wp2,         // In
                                                          mixt_frac(k), big_L_w_1, // In
                                                          mu_w_1(k), mu_w_2(k),    // Out
                                                          sigma_w_1_sqd(k),        // Out
                                                          sigma_w_2_sqd(k),        // Out
                                                          coef_sigma_w_1_sqd,      // Out
                                                          coef_sigma_w_2_sqd);     // Out

          new_tsdadg_pdf_detail::calc_respnder_parameters(thlm(k), thlp2(k),         // In
                                                          Skthl(k), sgn_wpthlp,      // In
                                                          mixt_frac(k), big_L_thl_1, // In
                                                          mu_thl_1(k), mu_thl_2(k),  // Out
                                                          sigma_thl_1_sqd(k),        // Out
                                                          sigma_thl_2_sqd(k),        // Out
                                                          coef_sigma_thl_1_sqd,      // Out
                                                          coef_sigma_thl_2_sqd);     // Out

        } else {

          // The variable thl has the greatest magnitude of skewness.

          calc_setter_parameters(thlm(k), thlp2(k),                // In
                                 Skthl(k), sgn_wpthlp,             // In
                                 big_L_thl_1, big_L_thl_2,         // In
                                 mu_thl_1(k), mu_thl_2(k),         // Out
                                 sigma_thl_1_sqd(k),               // Out
                                 sigma_thl_2_sqd(k), mixt_frac(k), // Out
                                 coef_sigma_thl_1_sqd,             // Out
                                 coef_sigma_thl_2_sqd);            // Out

          new_tsdadg_pdf_detail::calc_respnder_parameters(wm(k), wp2(k),           // In
                                                          Skw(k), sgn_wp2,         // In
                                                          mixt_frac(k), big_L_w_1, // In
                                                          mu_w_1(k), mu_w_2(k),    // Out
                                                          sigma_w_1_sqd(k),        // Out
                                                          sigma_w_2_sqd(k),        // Out
                                                          coef_sigma_w_1_sqd,      // Out
                                                          coef_sigma_w_2_sqd);     // Out

          new_tsdadg_pdf_detail::calc_respnder_parameters(rtm(k), rtp2(k),          // In
                                                          Skrt(k), sgn_wprtp,       // In
                                                          mixt_frac(k), big_L_rt_1, // In
                                                          mu_rt_1(k), mu_rt_2(k),   // Out
                                                          sigma_rt_1_sqd(k),        // Out
                                                          sigma_rt_2_sqd(k),        // Out
                                                          coef_sigma_rt_1_sqd,      // Out
                                                          coef_sigma_rt_2_sqd);     // Out
        }

        // Retain the Native pre-clipping diagnostic payload and field order.
        negative_variances(k, 0) = zero;
        if (sigma_w_1_sqd(k) < zero) {
          negative_variances(k, 0) = sigma_w_1_sqd(k);
          ++count;
          sigma_w_1_sqd(k) = zero;
          coef_sigma_w_1_sqd = zero;
        }
        // Retain the Native pre-clipping diagnostic payload and field order.
        negative_variances(k, 1) = zero;
        if (sigma_w_2_sqd(k) < zero) {
          negative_variances(k, 1) = sigma_w_2_sqd(k);
          ++count;
          sigma_w_2_sqd(k) = zero;
          coef_sigma_w_2_sqd = zero;
        }
        // Retain the Native pre-clipping diagnostic payload and field order.
        negative_variances(k, 2) = zero;
        if (sigma_rt_1_sqd(k) < zero) {
          negative_variances(k, 2) = sigma_rt_1_sqd(k);
          ++count;
          sigma_rt_1_sqd(k) = zero;
          coef_sigma_rt_1_sqd = zero;
        }
        // Retain the Native pre-clipping diagnostic payload and field order.
        negative_variances(k, 3) = zero;
        if (sigma_rt_2_sqd(k) < zero) {
          negative_variances(k, 3) = sigma_rt_2_sqd(k);
          ++count;
          sigma_rt_2_sqd(k) = zero;
          coef_sigma_rt_2_sqd = zero;
        }
        // Retain the Native pre-clipping diagnostic payload and field order.
        negative_variances(k, 4) = zero;
        if (sigma_thl_1_sqd(k) < zero) {
          negative_variances(k, 4) = sigma_thl_1_sqd(k);
          ++count;
          sigma_thl_1_sqd(k) = zero;
          coef_sigma_thl_1_sqd = zero;
        }
        // Retain the Native pre-clipping diagnostic payload and field order.
        negative_variances(k, 5) = zero;
        if (sigma_thl_2_sqd(k) < zero) {
          negative_variances(k, 5) = sigma_thl_2_sqd(k);
          ++count;
          sigma_thl_2_sqd(k) = zero;
          coef_sigma_thl_2_sqd = zero;
        }
      },
      warning_count);
  // Stage warning payload only when at least one variance is negative. Print
  // levels first, then w1/w2/rt1/rt2/thl1/thl2, using pre-clipping values.
  if (warning_count != 0) {
    const auto warnings_host = to_host_space(negative_variances);
    for (int k = 0; k < nz; ++k) {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      if (warnings_host(k, 0) < zero) {
        diagnostic << "WARNING:  New TSDADG PDF.  The variance of w in the 1st PDF component is "
                      "negative and is being clipped to 0."
                   << '\n';
        diagnostic << "sigma_w_1^2 (before clipping) = " << warnings_host(k, 0) << '\n';
      }
      if (warnings_host(k, 1) < zero) {
        diagnostic << "WARNING:  New TSDADG PDF.  The variance of w in the 2nd PDF component is "
                      "negative and is being clipped to 0."
                   << '\n';
        diagnostic << "sigma_w_2^2 (before clipping) = " << warnings_host(k, 1) << '\n';
      }
      if (warnings_host(k, 2) < zero) {
        diagnostic << "WARNING:  New TSDADG PDF.  The variance of rt in the 1st PDF component is "
                      "negative and is being clipped to 0."
                   << '\n';
        diagnostic << "sigma_rt_1^2 (before clipping) = " << warnings_host(k, 2) << '\n';
      }
      if (warnings_host(k, 3) < zero) {
        diagnostic << "WARNING:  New TSDADG PDF.  The variance of rt in the 2nd PDF component is "
                      "negative and is being clipped to 0."
                   << '\n';
        diagnostic << "sigma_rt_2^2 (before clipping) = " << warnings_host(k, 3) << '\n';
      }
      if (warnings_host(k, 4) < zero) {
        diagnostic << "WARNING:  New TSDADG PDF.  The variance of thl in the 1st PDF component is "
                      "negative and is being clipped to 0."
                   << '\n';
        diagnostic << "sigma_thl_1^2 (before clipping) = " << warnings_host(k, 4) << '\n';
      }
      if (warnings_host(k, 5) < zero) {
        diagnostic << "WARNING:  New TSDADG PDF.  The variance of thl in the 2nd PDF component is "
                      "negative and is being clipped to 0."
                   << '\n';
        diagnostic << "sigma_thl_2^2 (before clipping) = " << warnings_host(k, 5) << '\n';
      }
      std::cerr << diagnostic.str();
    }
  }
}

//=============================================================================
//
// DESCRIPTION OF THE METHOD FOR THE VARIABLE THAT SETS THE MIXTURE FRACTION
// =========================================================================
//
// There are five PDF parameters that need to be calculated for the setting
// variable, which are mu_x_1 (the mean of x in the 1st PDF component), mu_x_2
// (the mean of x in the 2nd PDF component), sigma_x_1 (the standard deviation
// of x in the 1st PDF component), sigma_x_2 (the standard deviation of x in
// the 2nd PDF component), and mixt_frac (the mixture fraction).  In order to
// solve for these five parameters, five equations are needed.  These five
// equations are the equations for <x>, <x'^2>, and <x'^3> as found by
// integrating over the PDF.  Additionally, two more equations, which involve
// tunable parameters L_x_1 and L_x_2, and which are used to help control the
// spread of the PDF component means in the 1st PDF component and the 2nd PDF
// component, respectively, are used in this equation set.  The five equations
// are:
//
// <x> = mixt_frac * mu_x_1 + ( 1 - mixt_frac ) * mu_x_2;
//
// <x'^2> = mixt_frac * ( ( mu_x_1 - <x> )^2 + sigma_x_1^2 )
//          + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> )^2 + sigma_x_2^2 );
//
// <x'^3> = mixt_frac * ( mu_x_1 - <x> )
//                    * ( ( mu_x_1 - <x> )^2 + 3 * sigma_x_1^2 )
//          + ( 1 - mixt_frac ) * ( mu_x_2 - <x> )
//                              * ( ( mu_x_2 - <x> )^2 + 3 * sigma_x_2^2 );
//
// mu_x_1 - <x> = L_x_1
//                * sqrt( ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
//                        / ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
//                * sqrt( <x'^2> ) * sgn( <w'x'> ); and
//
// mu_x_2 - <x> = -L_x_2
//                 * sqrt( ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
//                         / ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
//                 * sqrt( <x'^2> ) * sgn( <w'x'> );
//
// where 0 <= L_x_1 <= 1, 0 <= L_x_2 <= 1, Skx is the skewness of x, such that
// Skx = <x'^3> / <x'^2>^(3/2), and sgn( <w'x'> ) is the sign of <w'x'>, such
// that:
//
// sgn( <w'x'> ) = |  1, when <w'x'> >= 0;
//                 | -1, when <w'x'> < 0.
//
// The resulting equations for the five PDF parameters are:
//
// mu_x_1 = <x> + L_x_1
//                * sqrt( ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
//                        / ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
//                * sqrt( <x'^2> ) * sgn( <w'x'> );
//
// mu_x_2 = <x> - L_x_2
//                * sqrt( ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
//                        / ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
//                * sqrt( <x'^2> ) * sgn( <w'x'> );
//
// mixt_frac = 1 / ( 1 + abs( mu_x_1_nrmlized / mu_x_2_nrmlized ) );
//
// sigma_x_1 = sqrt( ( 1 - mixt_frac * mu_x_1_nrmlized^2
//                     - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
//                     + ( 1 - mixt_frac )
//                       * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
//                           - mu_x_1_nrmlized^2 / 3
//                           + mu_x_2_nrmlized^2 / 3 ) )
//                   * <x'^2> ); and
//
// sigma_x_2 = sqrt( ( 1 - mixt_frac * mu_x_1_nrmlized^2
//                     - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
//                     - mixt_frac
//                       * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
//                           - mu_x_1_nrmlized^2 / 3
//                           + mu_x_2_nrmlized^2 / 3 ) )
//                   * <x'^2> ); where
//
// mu_x_1_nrmlized = ( mu_x_1 - <x> ) / sqrt( <x'^2> ); and
//
// mu_x_2_nrmlized = ( mu_x_2 - <x> ) / sqrt( <x'^2> ).
//
//
// Notes:
//
// This method does NOT work for all values of L_x_1 and L_x_2 (where
// 0 <= L_x_1 <= 1 and 0 <= L_x_2 <= 1).  Only a subregion of this parameter
// space produces valid results.
//
// When both L_x_1 = 0 and L_x_2 = 0, mu_x_1 = mu_x_2 = <x> (which can only
// happen when Skx = 0).  In this scenario, the above equations for mixt_frac,
// sigma_x_1, and sigma_x_2 are all undefined.  In this special case, the
// distribution reduces to a single Gaussian, so the following values are set:
// mixt_frac = 1/2, sigma_x_1 = sqrt( <x'^2> ), and sigma_x_2 = sqrt( <x'^2> ).
//
//
// Tunable parameters:
//
// The parameter L_x_1 controls the 1st PDF component mean while L_x_2 controls
// the 2nd PDF component mean.  The equations involving the tunable parameters
// L_x_1 and L_x_2 (the mu_x_1 and mu_x_2 equations) are based on the values of
// mu_x_1 and mu_x_2 when sigma_x_1 = sigma_x_2 = 0.  In this scenario, the
// equation for mixture fraction reduces to:
//
// mixt_frac = (1/2) * ( 1 +/- Skx / sqrt( 4 + Skx^2 ) ).
//
// The +/- is dependent on the sign of ( mu_x_1 - <x> ) vs. ( mu_x_2 - <x> ).
// This is dependent on sgn( <w'x'> ), and the mixture fraction equation is
// written as:
//
// mixt_frac = (1/2) * ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ).
//
// Meanwhile, the equation for 1 - mixt_frac is:
//
// 1 - mixt_frac = (1/2) * ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ).
//
// When sigma_x_1 = sigma_x_2 = 0, the equations for mu_x_1 and mu_x_2 are:
//
// mu_x_1 = <x> + sqrt( ( 1 - mixt_frac ) / mixt_frac ) * sqrt( <x'^2> )
//                * sgn( <w'x'> ); and
//
// mu_x_2 = <x> - sqrt( mixt_frac / ( 1 - mixt_frac ) ) * sqrt( <x'^2> )
//                * sgn( <w'x'> ).
//
// Substituting the equations for mixt_frac and 1 - mixt_frac into the
// equations for mu_x_1 and mu_x_2 (when sigma_x_1 = sigma_x_2 = 0), the
// equations for mu_x_1 and mu_x_2 become:
//
// mu_x_1 = <x> + sqrt( ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
//                      / ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
//                * sqrt( <x'^2> ) * sgn( <w'x'> ); and
//
// mu_x_2 = <x> - sqrt( ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
//                      / ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
//                * sqrt( <x'^2> ) * sgn( <w'x'> ).
//
// These equations represent the maximum deviation of mu_x_1 and mu_x_2 from
// the overall mean, <x>.  The range of parameters of L_x_i is 0 <= L_x_i <= 1.
// When L_x_1 = L_x_2 = 0, the value of mu_x_1 = mu_x_2 = <x> (and the
// distribution becomes a single Gaussian).  When L_x_i = 1, the value of
// mu_x_i - <x> is at its maximum magnitude.
//
// The values of L_x_1 and L_x_2 are also calculated by skewness functions.
// Those functions are:
//
// L_x_1 = l_x_1 * abs( Skx ) / sqrt( 4 + Skx^2 ); and
// L_x_2 = l_x_2 * abs( Skx ) / sqrt( 4 + Skx^2 );
//
// where both l_x_1 and l_x_2 are tunable parameters.
//
// As previously stated, this method does not work for all combinations of
// L_x_1 and L_x_2, but rather only for a subregion of parameter space.  This
// applies to l_x_1 and l_x_2, as well.  The conditions on l_x_1 and l_x_2 are:
//
// 2/3 < l_x_1 < 1 and 0 < l_x_2 < 1; when Skx * sgn( <w'x'> ) >= 0; and
// 0 < l_x_1 < 1 and 2/3 < l_x_2 < 1; when Skx * sgn( <w'x'> ) < 0.
//
// The condition that l_x_1 > 2/3 prevents a negative PDF component variance
// when Skx = 0.
//
//
// Equations for PDF component standard deviations:
//
// The equations for the PDF component standard deviations can also be written
// as:
//
// sigma_x_1 = sqrt( coef_sigma_x_1_sqd * <x'^2> ); and
//
// sigma_x_2 = sqrt( coef_sigma_x_2_sqd * <x'^2> ); where
//
// coef_sigma_x_1_sqd = 1 - mixt_frac * mu_x_1_nrmlized^2
//                      - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
//                      + ( 1 - mixt_frac )
//                        * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
//                            - mu_x_1_nrmlized^2 / 3
//                            + mu_x_2_nrmlized^2 / 3 ); and
//
// coef_sigma_x_2_sqd = 1 - mixt_frac * mu_x_1_nrmlized^2
//                      - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
//                      - mixt_frac
//                        * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
//                            - mu_x_1_nrmlized^2 / 3
//                            + mu_x_2_nrmlized^2 / 3 ).
//
// The above equations can be substituted into an equation for a variable that
// has been derived by integrating over the PDF.  Many variables like this are
// used in parts of the predictive equation set.  These substitutions allow
// some terms to solved implicitly or semi-implicitly in the predictive
// equations.
//
//
//=============================================================================
KOKKOS_INLINE_FUNCTION void calc_setter_parameters(real_t xm, real_t xp2, real_t Skx,
                                                   real_t sgn_wpxp,                    // In
                                                   real_t big_L_x_1, real_t big_L_x_2, // In
                                                   real_t& mu_x_1, real_t& mu_x_2,
                                                   real_t& sigma_x_1_sqd,                    // Out
                                                   real_t& sigma_x_2_sqd, real_t& mixt_frac, // Out
                                                   real_t& coef_sigma_x_1_sqd,               // Out
                                                   real_t& coef_sigma_x_2_sqd                // Out
) {
  // Description:
  // Calculates the PDF component means, the PDF component variances,
  // and the mixture fraction for the variable that sets the PDF.
  //
  // References:
  //-----------------------------------------------------------------------
  // xm/mu_x_i have units of x; xp2/sigma_x_i_sqd have squared units.
  // Skx, sgn_wpxp, L coefficients, normalized means and variance factors are [-].
  real_t mu_x_1_nrmlized, mu_x_2_nrmlized;                     // Normalized component means [-].
  real_t factor_plus, factor_minus;                            // Dimensionless mean-spread factors.
  real_t sqrt_factor_plus_ov_minus, sqrt_factor_minus_ov_plus; // Square-root ratios [-].
  real_t mu_x_1_nrmlized_thresh; // Minimum-magnitude normalized mean [-].
  // Calculate the factors in the PDF component mean equations.
  factor_plus = one + Skx * sgn_wpxp / Kokkos::sqrt(four + Skx * Skx);

  factor_minus = one - Skx * sgn_wpxp / Kokkos::sqrt(four + Skx * Skx);

  sqrt_factor_plus_ov_minus = Kokkos::sqrt(factor_plus / factor_minus);
  sqrt_factor_minus_ov_plus = Kokkos::sqrt(factor_minus / factor_plus);

  // Calculate the normalized mean of x in the 1st PDF component.
  mu_x_1_nrmlized = big_L_x_1 * sqrt_factor_plus_ov_minus * sgn_wpxp;

  // Calculate the normalized mean of x in the 2nd PDF component.
  mu_x_2_nrmlized = -big_L_x_2 * sqrt_factor_minus_ov_plus * sgn_wpxp;

  // Calculate the mean of x in the 1st PDF component.
  mu_x_1 = xm + mu_x_1_nrmlized * Kokkos::sqrt(xp2);

  // Calculate the mean of x in the 2nd PDF component.
  mu_x_2 = xm + mu_x_2_nrmlized * Kokkos::sqrt(xp2);

  // Calculate the mixture fraction.
  if (Kokkos::abs(mu_x_1_nrmlized) >= eps && Kokkos::abs(mu_x_2_nrmlized) >= eps) {
    mixt_frac = one / (one + Kokkos::abs(mu_x_1_nrmlized / mu_x_2_nrmlized));
  } else if (Kokkos::abs(mu_x_1_nrmlized) >= eps && Kokkos::abs(mu_x_2_nrmlized) < eps) {
    mixt_frac = one / (one + Kokkos::abs(mu_x_1_nrmlized / eps));
  } else if (Kokkos::abs(mu_x_1_nrmlized) < eps && Kokkos::abs(mu_x_2_nrmlized) >= eps) {
    mixt_frac = one / (one + Kokkos::abs(eps / mu_x_2_nrmlized));
  } else { // abs( mu_x_1_nrmlized ) < eps and abs( mu_x_2_nrmlized ) < eps
    mixt_frac = one_half;
  }

  // Use a minimum magnitude value of mu_x_1_nrmlized in the denominator of a
  // term in the PDF component variance equations in order to prevent a
  // divide-by-zero error.
  if (mu_x_1_nrmlized >= zero) {
    mu_x_1_nrmlized_thresh = Kokkos::max(mu_x_1_nrmlized, real_t(eps));
  } else { // mu_x_1_nrmlized < 0
    mu_x_1_nrmlized_thresh = Kokkos::min(mu_x_1_nrmlized, real_t(-eps));
  } // mu_x_1_nrmlized >= 0

  // Calculate the variance of x in the 1st PDF component.
  coef_sigma_x_1_sqd = one - mixt_frac * (mu_x_1_nrmlized * mu_x_1_nrmlized) -
                       (one - mixt_frac) * (mu_x_2_nrmlized * mu_x_2_nrmlized) +
                       (one - mixt_frac) * (Skx / (three * mixt_frac * mu_x_1_nrmlized_thresh) -
                                            mu_x_1_nrmlized * mu_x_1_nrmlized / three +
                                            mu_x_2_nrmlized * mu_x_2_nrmlized / three);

  sigma_x_1_sqd = coef_sigma_x_1_sqd * xp2;

  // Calculate the variance of x in the 2nd PDF component.
  coef_sigma_x_2_sqd = one - mixt_frac * (mu_x_1_nrmlized * mu_x_1_nrmlized) -
                       (one - mixt_frac) * (mu_x_2_nrmlized * mu_x_2_nrmlized) -
                       mixt_frac * (Skx / (three * mixt_frac * mu_x_1_nrmlized_thresh) -
                                    mu_x_1_nrmlized * mu_x_1_nrmlized / three +
                                    mu_x_2_nrmlized * mu_x_2_nrmlized / three);

  sigma_x_2_sqd = coef_sigma_x_2_sqd * xp2;
}

// =============================================================================
KOKKOS_INLINE_FUNCTION void calc_L_x_Skx_fnc(real_t Skx, real_t sgn_wpxp,            // In
                                             real_t small_l_x_1, real_t small_l_x_2, // In
                                             real_t& big_L_x_1, real_t& big_L_x_2    // Out
) {
  // Description:
  // Calculates the values of big_L_x_1 and big_L_x_2 as functions of Skx.
  //
  // References:
  //-----------------------------------------------------------------------
  // Skewness, covariance sign, input l parameters and output L parameters are [-].
  real_t Skx_fnc_factor; // Skewness-dependent spread factor [-].
  // The values of L_x_1 and L_x_2 are calculated by skewness functions.
  // Those functions are:

  // L_x_1 = l_x_1 * abs( Skx ) / sqrt( 4 + Skx^2 ); and
  // L_x_2 = l_x_2 * abs( Skx ) / sqrt( 4 + Skx^2 ).

  // The conditions on l_x_1 and l_x_2 are:

  // 2/3 < l_x_1 < 1 and 0 < l_x_2 < 1; when Skx * sgn( <w'x'> ) >= 0; and
  // 0 < l_x_1 < 1 and 2/3 < l_x_2 < 1; when Skx * sgn( <w'x'> ) < 0.

  // For simplicity, this can also be accomplished by setting 2/3 < l_x_1 < 1
  // and 0 < l_x_2 < 1, and then using the following equations.

  // When Skx * sgn( <w'x'> ) >= 0:
  // L_x_1 = l_x_1 * abs( Skx ) / sqrt( 4 + Skx^2 ); and
  // L_x_2 = l_x_2 * abs( Skx ) / sqrt( 4 + Skx^2 );

  // otherwise, when Skx * sgn( <w'x'> ) < 0, switch l_x_1 and l_x_2:
  // L_x_1 = l_x_2 * abs( Skx ) / sqrt( 4 + Skx^2 ); and
  // L_x_2 = l_x_1 * abs( Skx ) / sqrt( 4 + Skx^2 ).
  Skx_fnc_factor = Kokkos::abs(Skx) / Kokkos::sqrt(four + Skx * Skx);

  if (Skx * sgn_wpxp >= zero) {
    big_L_x_1 = small_l_x_1 * Skx_fnc_factor;
    big_L_x_2 = small_l_x_2 * Skx_fnc_factor;
  } else { // Skx * sgn_wpxp < 0
    big_L_x_1 = small_l_x_2 * Skx_fnc_factor;
    big_L_x_2 = small_l_x_1 * Skx_fnc_factor;
  }
}

//=============================================================================
//
// DESCRIPTION OF THE METHOD FOR EACH RESPONDING VARIABLE
// ======================================================
//
// In order to find equations for the four PDF parameters for each responding
// variable, which are mu_x_1, mu_x_2, sigma_x_1, and sigma_x_2 (where x stands
// for a responding variable here), four equations are needed.  These four
// equations are the equations for <x>, <x'^2>, and <x'^3> as found by
// integrating over the PDF.  Additionally, one more equation, which involves
// tunable parameter L_x_1, and which is used to help control the mean of the
// 1st PDF component, is used in this equation set.  The four equations are:
//
// <x> = mixt_frac * mu_x_1 + ( 1 - mixt_frac ) * mu_x_2;
//
// <x'^2> = mixt_frac * ( ( mu_x_1 - <x> )^2 + sigma_x_1^2 )
//          + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> )^2 + sigma_x_2^2 );
//
// <x'^3> = mixt_frac * ( mu_x_1 - <x> )
//                    * ( ( mu_x_1 - <x> )^2 + 3 * sigma_x_1^2 )
//          + ( 1 - mixt_frac ) * ( mu_x_2 - <x> )
//                              * ( ( mu_x_2 - <x> )^2 + 3 * sigma_x_2^2 ); and
//
// mu_x_1 - <x> = L_x_1
//                * sqrt( ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
//                        / ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
//                * sqrt( <x'^2> ) * sgn( <w'x'> );
//
// where 0 <= L_x_1 <= 1, Skx is the skewness of x, such that
// Skx = <x'^3> / <x'^2>^(3/2), and sgn( <w'x'> ) is given by:
//
// sgn( <w'x'> ) = |  1, when <w'x'> >= 0;
//                 | -1, when <w'x'> < 0.
//
// The resulting equations for the four PDF parameters are:
//
// mu_x_1 = <x> + L_x_1
//                * sqrt( ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
//                        / ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
//                * sqrt( <x'^2> ) * sgn( <w'x'> );
//
// mu_x_2 = <x> - ( mixt_frac / ( 1 - mixt_frac ) ) * ( mu_x_1 - <x> );
//
// sigma_x_1 = sqrt( ( 1 - mixt_frac * mu_x_1_nrmlized^2
//                     - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
//                     + ( 1 - mixt_frac )
//                       * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
//                           - mu_x_1_nrmlized^2 / 3
//                           + mu_x_2_nrmlized^2 / 3 ) )
//                   * <x'^2> ); and
//
// sigma_x_2 = sqrt( ( 1 - mixt_frac * mu_x_1_nrmlized^2
//                     - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
//                     - mixt_frac
//                       * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
//                           - mu_x_1_nrmlized^2 / 3
//                           + mu_x_2_nrmlized^2 / 3 ) )
//                   * <x'^2> ); where
//
// mu_x_1_nrmlized = ( mu_x_1 - <x> ) / sqrt( <x'^2> ); and
//
// mu_x_2_nrmlized = ( mu_x_2 - <x> ) / sqrt( <x'^2> ).
//
//
// Notes:
//
// When L_x_1 = 0, mu_x_1 = mu_x_2 = <x>, sigma_x_1^2 = sigma_x_2^2 = <x'^2>,
// and the distribution reduces to a single Gaussian.
//
//
// Equations for PDF component standard deviations:
//
// The equations for the PDF component standard deviations can also be written
// as:
//
// sigma_x_1 = sqrt( coef_sigma_x_1_sqd * <x'^2> ); and
//
// sigma_x_2 = sqrt( coef_sigma_x_2_sqd * <x'^2> ); where
//
// coef_sigma_x_1_sqd = 1 - mixt_frac * mu_x_1_nrmlized^2
//                      - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
//                      + ( 1 - mixt_frac )
//                        * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
//                            - mu_x_1_nrmlized^2 / 3
//                            + mu_x_2_nrmlized^2 / 3 ); and
//
// coef_sigma_x_2_sqd = 1 - mixt_frac * mu_x_1_nrmlized^2
//                      - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
//                      - mixt_frac
//                        * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
//                            - mu_x_1_nrmlized^2 / 3
//                            + mu_x_2_nrmlized^2 / 3 ).
//
// The above equations can be substituted into an equation for a variable that
// has been derived by integrating over the PDF.  Many variables like this are
// used in parts of the predictive equation set.  These substitutions allow
// some terms to solved implicitly or semi-implicitly in the predictive
// equations.
//
//
//=============================================================================
namespace new_tsdadg_pdf_detail {
KOKKOS_INLINE_FUNCTION void calc_respnder_parameters(real_t xm, real_t xp2, real_t Skx,
                                                     real_t sgn_wpxp,                    // In
                                                     real_t mixt_frac, real_t big_L_x_1, // In
                                                     real_t& mu_x_1, real_t& mu_x_2,     // Out
                                                     real_t& sigma_x_1_sqd,
                                                     real_t& sigma_x_2_sqd,      // Out
                                                     real_t& coef_sigma_x_1_sqd, // Out
                                                     real_t& coef_sigma_x_2_sqd  // Out
) {
  // Description:
  // Calculates component means and variances for a responding variable using
  // the supplied mixture fraction. The mixture is input and is unchanged.
  //
  // References:
  //-----------------------------------------------------------------------
  // xm/mu_x_i have units of x; xp2/sigma_x_i_sqd have squared units.
  // Skx, covariance sign, mixture, L and normalized quantities are dimensionless.
  real_t mu_x_1_nrmlized, mu_x_2_nrmlized; // Normalized component means [-].
  real_t factor_plus, factor_minus;        // Dimensionless mean-spread factors.
  real_t sqrt_factor_plus_ov_minus;        // Square-root factor ratio [-].
  real_t mu_x_1_nrmlized_thresh;           // Minimum-magnitude normalized mean [-].
  // Calculate the factors in the PDF component mean equations.
  factor_plus = one + Skx * sgn_wpxp / Kokkos::sqrt(four + Skx * Skx);

  factor_minus = one - Skx * sgn_wpxp / Kokkos::sqrt(four + Skx * Skx);

  sqrt_factor_plus_ov_minus = Kokkos::sqrt(factor_plus / factor_minus);

  // Calculate the normalized mean of x in the 1st PDF component.
  mu_x_1_nrmlized = big_L_x_1 * sqrt_factor_plus_ov_minus * sgn_wpxp;

  // Calculate the normalized mean of x in the 2nd PDF component.
  mu_x_2_nrmlized = -(mixt_frac / (one - mixt_frac)) * mu_x_1_nrmlized;

  // Calculate the mean of x in the 1st PDF component.
  mu_x_1 = xm + mu_x_1_nrmlized * Kokkos::sqrt(xp2);

  // Calculate the mean of x in the 2nd PDF component.
  mu_x_2 = xm + mu_x_2_nrmlized * Kokkos::sqrt(xp2);

  // Use a minimum magnitude value of mu_x_1_nrmlized in the denominator of a
  // term in the PDF component variance equations in order to prevent a
  // divide-by-zero error.
  if (mu_x_1_nrmlized >= zero) {
    mu_x_1_nrmlized_thresh = Kokkos::max(mu_x_1_nrmlized, real_t(eps));
  } else { // mu_x_1_nrmlized < 0
    mu_x_1_nrmlized_thresh = Kokkos::min(mu_x_1_nrmlized, real_t(-eps));
  } // mu_x_1_nrmlized >= 0

  // Calculate the variance of x in the 1st PDF component.
  coef_sigma_x_1_sqd = one - mixt_frac * (mu_x_1_nrmlized * mu_x_1_nrmlized) -
                       (one - mixt_frac) * (mu_x_2_nrmlized * mu_x_2_nrmlized) +
                       (one - mixt_frac) * (Skx / (three * mixt_frac * mu_x_1_nrmlized_thresh) -
                                            mu_x_1_nrmlized * mu_x_1_nrmlized / three +
                                            mu_x_2_nrmlized * mu_x_2_nrmlized / three);

  sigma_x_1_sqd = coef_sigma_x_1_sqd * xp2;

  // Calculate the variance of x in the 2nd PDF component.
  coef_sigma_x_2_sqd = one - mixt_frac * (mu_x_1_nrmlized * mu_x_1_nrmlized) -
                       (one - mixt_frac) * (mu_x_2_nrmlized * mu_x_2_nrmlized) -
                       mixt_frac * (Skx / (three * mixt_frac * mu_x_1_nrmlized_thresh) -
                                    mu_x_1_nrmlized * mu_x_1_nrmlized / three +
                                    mu_x_2_nrmlized * mu_x_2_nrmlized / three);

  sigma_x_2_sqd = coef_sigma_x_2_sqd * xp2;
}
} // namespace new_tsdadg_pdf_detail

} // namespace clubb_kokkos
