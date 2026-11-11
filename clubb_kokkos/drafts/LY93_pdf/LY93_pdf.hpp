// Close a two-component PDF with the Lewellen and Yoh (1993) method.
//
// LY93_driver finds one mixture fraction from the largest absolute skewness,
// then diagnoses component means and variances for w, total water rt and liquid
// water potential temperature thl. calc_mixt_frac_LY93 exposes the iterative
// mixture calculation; calc_params_LY93 closes one variable with a given mixture.
// Means retain their input units; variances have squared input units.
//
// The Native APIs take one nz-level profile. Generic rank-one Views accept the
// existing column subviews without changing column-first LayoutLeft ownership.
// Scratch and returned mixtures belong to each call and retain input memory
// space. Bisection is serial within each independent level, without a new cap.
// Compilation, numerical agreement and performance remain unverified.
// Reference: Lewellen and Yoh (1993), J. Atmos. Sci. 50(9), 1228–1237.
#pragma once

#include "array_types.hpp"
#include "constants_clubb.hpp"
#include <Kokkos_MathematicalFunctions.hpp>

namespace clubb_kokkos {

// Forward declarations preserve the Native public procedure order.
template <class SkMax>
Kokkos::View<real_t*, Kokkos::LayoutLeft, typename SkMax::memory_space>
calc_mixt_frac_LY93(int nz, SkMax Sk_max);

template <class Xm, class Xp2, class SkX, class MixtFrac, class MuX1, class MuX2, class SigmaX1,
          class SigmaX2>
void calc_params_LY93(int nz, Xm xm, Xp2 xp2, SkX Skx, MixtFrac mixt_frac, // In
                      MuX1 mu_x_1, MuX2 mu_x_2,                            // Out
                      SigmaX1 sigma_x_1_sqd, SigmaX2 sigma_x_2_sqd);       // Out

// =============================================================================
template <class Wm, class Rtm, class Thlm, class Wp2, class Rtp2, class Thlp2, class SkW,
          class SkRt, class SkThl, class MuW1, class MuW2, class MuRt1, class MuRt2, class MuThl1,
          class MuThl2, class SigmaW1, class SigmaW2, class SigmaRt1, class SigmaRt2,
          class SigmaThl1, class SigmaThl2, class MixtFrac>
void LY93_driver(int nz, Wm wm, Rtm rtm, Thlm thlm, Wp2 wp2, Rtp2 rtp2, // In
                 Thlp2 thlp2, SkW Skw, SkRt Skrt, SkThl Skthl,          // In
                 MuW1 mu_w_1, MuW2 mu_w_2,                              // Out
                 MuRt1 mu_rt_1, MuRt2 mu_rt_2,                          // Out
                 MuThl1 mu_thl_1, MuThl2 mu_thl_2,                      // Out
                 SigmaW1 sigma_w_1_sqd, SigmaW2 sigma_w_2_sqd,          // Out
                 SigmaRt1 sigma_rt_1_sqd, SigmaRt2 sigma_rt_2_sqd,      // Out
                 SigmaThl1 sigma_thl_1_sqd, SigmaThl2 sigma_thl_2_sqd,  // Out
                 MixtFrac mixt_frac) {                                  // Out
  // Description:
  // Calculates the mixture fraction and the PDF component means and PDF
  // component variances of w, rt, and theta-l following Lewellen and Yoh.
  //
  // References:
  // Lewellen, W. S. and Yoh, S., 1993. Binormal Model of Ensemble Partial
  // Cloudiness. J. Atmos. Sci., 50, 9, 1228--1237.
  // ---------------------------------------------------------------------------
  // wm/wp2/Skw: mean [m/s], variance [m^2/s^2] and skewness [-] of w.
  // rtm/rtp2/Skrt: mean [kg/kg], variance [(kg/kg)^2] and skewness [-] of rt.
  // thlm/thlp2/Skthl: mean [K], variance [K^2] and skewness [-] of thl.
  // Component means/variances retain these units; mixt_frac is dimensionless.
  static_assert(Wm::rank == 1 && Rtm::rank == 1 && Thlm::rank == 1 && Wp2::rank == 1 &&
                Rtp2::rank == 1 && Thlp2::rank == 1 && SkW::rank == 1 && SkRt::rank == 1 &&
                SkThl::rank == 1);
  using Execution = typename Wm::execution_space;
  Kokkos::View<real_t*, Kokkos::LayoutLeft, typename Wm::memory_space> Sk_max(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "LY93_driver::Sk_max"), nz);

  // Find the maximum of the magnitudes of skewness.
  Kokkos::parallel_for(
      "LY93_driver::Sk_max", Kokkos::RangePolicy<Execution>(0, nz), KOKKOS_LAMBDA(const int k) {
        Sk_max(k) = Kokkos::max(Kokkos::max(Kokkos::abs(Skw(k)), Kokkos::abs(Skrt(k))),
                                Kokkos::abs(Skthl(k)));
      });

  // Calculate mixture fraction.
  Kokkos::deep_copy(mixt_frac, calc_mixt_frac_LY93(nz, Sk_max));

  // Calculate the PDF parameters for w.
  calc_params_LY93(nz, wm, wp2, Skw, mixt_frac,   // In
                   mu_w_1, mu_w_2,                // Out
                   sigma_w_1_sqd, sigma_w_2_sqd); // Out

  // Calculate the PDF parameters for rt.
  calc_params_LY93(nz, rtm, rtp2, Skrt, mixt_frac,  // In
                   mu_rt_1, mu_rt_2,                // Out
                   sigma_rt_1_sqd, sigma_rt_2_sqd); // Out

  // Calculate the PDF parameters for thl.
  calc_params_LY93(nz, thlm, thlp2, Skthl, mixt_frac, // In
                   mu_thl_1, mu_thl_2,                // Out
                   sigma_thl_1_sqd, sigma_thl_2_sqd); // Out
}

// =============================================================================
template <class SkMax>
Kokkos::View<real_t*, Kokkos::LayoutLeft, typename SkMax::memory_space>
calc_mixt_frac_LY93(int nz, SkMax Sk_max) {
  // Description:
  // Calculates mixture fraction iteratively according to Lewellen and Yoh.
  // Sk_max is the maximum absolute skewness [-]; the returned mixture is [-].
  //
  // References:
  // Eq. (21) of Lewellen, W. S. and Yoh, S., 1993. Binormal Model of Ensemble
  // Partial Cloudiness. J. Atmos. Sci., 50, 9, 1228--1237.
  // ---------------------------------------------------------------------------
  static_assert(SkMax::rank == 1);
  using Execution = typename SkMax::execution_space;
  Kokkos::View<real_t*, Kokkos::LayoutLeft, typename SkMax::memory_space> mixt_frac(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_mixt_frac_LY93::mixt_frac"), nz);

  // Tolerance for mixture fraction in solver [-]. Native initializes this local
  // once but never changes it; no mutable saved storage is needed in C++.
  constexpr real_t LY_mixt_frac_tol = real_t(1.0e-4);
  Kokkos::parallel_for(
      "calc_mixt_frac_LY93", Kokkos::RangePolicy<Execution>(0, nz), KOKKOS_LAMBDA(const int k) {
        real_t mixt_frac_low;   // Low mixture bound in the iterative solver [-].
        real_t mixt_frac_high;  // High mixture bound in the iterative solver [-].
        real_t expr_equal_zero; // mixture^6 - Sk_max^2*(1-mixture) [-].
        if (Sk_max(k) > real_t(0.84)) {
          mixt_frac_low = one_half;
          mixt_frac_high = one;
          while (true) { // Solve iteratively for mixture fraction.
            mixt_frac(k) = one_half * (mixt_frac_low + mixt_frac_high);
            expr_equal_zero =
                Kokkos::pow(mixt_frac(k), 6) - (Sk_max(k) * Sk_max(k)) * (one - mixt_frac(k));
            if (Kokkos::abs(expr_equal_zero) < LY_mixt_frac_tol) {
              // Mixture fraction has been solved for within the specified
              // tolerance. Preserve the Native convergence test and no cap.
              break;
            } else {
              if (expr_equal_zero > zero) {
                mixt_frac_high = mixt_frac(k);
              } else { // expr_equal_zero < 0
                mixt_frac_low = mixt_frac(k);
              }
            }
          }
        } else { // Sk_max <= 0.84
          mixt_frac(k) = three_fourths;
        }
      });
  return mixt_frac;
}

// =============================================================================
template <class Xm, class Xp2, class SkX, class MixtFrac, class MuX1, class MuX2, class SigmaX1,
          class SigmaX2>
void calc_params_LY93(int nz, Xm xm, Xp2 xp2, SkX Skx, MixtFrac mixt_frac, // In
                      MuX1 mu_x_1, MuX2 mu_x_2,                            // Out
                      SigmaX1 sigma_x_1_sqd, SigmaX2 sigma_x_2_sqd) {      // Out
  // Description:
  // Calculates the PDF component means and PDF component variances for
  // variable x according to Lewellen and Yoh. xm has the units of x, xp2 has
  // squared units, and Skx/mixt_frac are dimensionless. Outputs retain the
  // mean/variance units of their corresponding inputs.
  //
  // References:
  // Eq. (14), Eq. (15), Eq. (16), Eq. (17), and Eq. (18) of
  // Lewellen, W. S. and Yoh, S., 1993. Binormal Model of Ensemble Partial
  // Cloudiness. J. Atmos. Sci., 50, 9, 1228--1237.
  // ---------------------------------------------------------------------------
  static_assert(Xm::rank == 1 && Xp2::rank == 1 && SkX::rank == 1 && MixtFrac::rank == 1 &&
                MuX1::rank == 1 && MuX2::rank == 1 && SigmaX1::rank == 1 && SigmaX2::rank == 1);
  using Execution = typename Xm::execution_space;
  // All source array operations are independent by level. Keep sgn_Skx and B_x
  // as distinct work-item scalars, preserving formula order in a single kernel.
  Kokkos::parallel_for(
      "calc_params_LY93", Kokkos::RangePolicy<Execution>(0, nz), KOKKOS_LAMBDA(const int k) {
        real_t sgn_Skx; // Sign of Skx [-].
        real_t B_x;     // Spread of the PDF component means [units of x].

        // Find the sign of Skx.
        if (Skx(k) >= zero) {
          sgn_Skx = one;
        } else { // Skx < 0
          sgn_Skx = -one;
        }

        // Calculate B_x, the LY function for the spread of the PDF component means.
        B_x = sgn_Skx * Kokkos::sqrt(xp2(k)) *
              Kokkos::pow(Kokkos::abs(Skx(k)) / (one - mixt_frac(k)), one_third);

        // Calculate the mean of x in the 1st PDF component.
        mu_x_1(k) = xm(k) - B_x * (one - mixt_frac(k));

        // Calculate the mean of x in the 2nd PDF component.
        mu_x_2(k) = xm(k) + B_x * mixt_frac(k);

        // Calculate the variance of x in the 1st PDF component.
        sigma_x_1_sqd(k) = xp2(k) - (B_x * B_x) * (one - mixt_frac(k)) *
                                        (one + mixt_frac(k) + (mixt_frac(k) * mixt_frac(k))) /
                                        (three * mixt_frac(k));

        // Calculate the variance of x in the 2nd PDF component.
        sigma_x_2_sqd(k) =
            xp2(k) + (B_x * B_x) * ((one - mixt_frac(k)) * (one - mixt_frac(k))) / three;
      });
}

} // namespace clubb_kokkos
