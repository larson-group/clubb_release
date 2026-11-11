// Shared scalar utilities used by the Khairoutdinov-Kogan PDF formulas.
//
// factorial preserves Native's floating-point gamma implementation; it does
// not replace it with an integer product or a different standard math library.
#pragma once
#include "Microphys/KK_microphys/Parabolic.hpp"
#include "Microphys/KK_microphys/parameters_KK.hpp"
#include "CLUBB_core/saturation.hpp"
namespace clubb_kokkos {
KOKKOS_INLINE_FUNCTION double factorial(int num) {
  // Calculate factorial through the Native floating-point gamma routine.
  return parabolic::gamma(static_cast<double>((num + 1)));
}
} // namespace clubb_kokkos

namespace clubb_kokkos::kk_engine {
inline double Dv_fnc(double order, double argument, const parameters_KK& kk) {
  // Compute Native Dv using Algorithm850 U/V and its configured tolerance.
  // All supported source call sites have order<=-2 and a>=1.5; arbitrary
  // direct calls outside that proven closure fail at the Parabolic entry.
  constexpr int scaling  = 0;
  constexpr double limit = 1.e308;
  parabolic::ParabolicArray<2> uaxx, vaxx;
  int ierr;
  const double epss = kk.l_high_accuracy_parab_cyl_fnc ? 1.e-15 : 1.e-4;
  double Dv_fnc;
  if (argument <= zero_dp) {
    parabolic::parab(-order - one_half_dp, -argument, scaling, uaxx, vaxx, ierr, epss);
    Dv_fnc = vaxx(0) / ((one_dp / pi_dp) * parabolic::gamma(-order)) -
             Kokkos::sin(pi_dp * (-order - one_half_dp)) * uaxx(0);
  } else {
    parabolic::parab(-order - one_half_dp, argument, scaling, uaxx, vaxx, ierr, epss);
    Dv_fnc = uaxx(0);
  }
  if (ierr != 0)
    Dv_fnc = limit;
  return Dv_fnc;
}
} // namespace clubb_kokkos::kk_engine

namespace clubb_kokkos::kk_engine {
inline real_t G_T_p(real_t T_in_K, real_t p_in_Pa, int saturation_formula,
                    const parameters_KK& kk) {
  // Description:
  // Calculate G(T,p) from the drop radius growth equation.  This is used as a
  // coefficient in the KK evaporation equation.  The equation for G(T,p) is
  // taken from Rogers and Yau (1989).
  // References:
  // Eq. (22) of Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
  // Parameterization in a Large-Eddy Simulation Model of Marine Stratocumulus.
  // Mon. Wea. Rev., 128, 229--243.
  //
  // Eq. (7.17) and Eq. (7.18) of Rogers, R. R. and M. K. Yau, 1989:  A Short
  // Course in Cloud Physics. 3rd edition, Butterworth-Heinemann, 290 pp.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t G_T_p;
  real_t Ka;
  real_t Dv;
  real_t Fk;
  real_t Fd;
  real_t esatv;
  real_t Celsius;

  // Temperature in degrees Celsius.
  Celsius = (T_in_K - T_freeze_K);

  // Coefficient of thermal conductivity of air.
  Ka = ((real_t(5.69) + (real_t(0.017) * Celsius)) * real_t(0.00001));
  Ka = ((real_t(4.1868) * one_hundred) * Ka);

  // Coefficient of diffusion of water vapor in air.
  Dv = ((real_t(0.221) * Kokkos::pow((T_in_K / real_t(273.16)), real_t(1.94))) *
        (real_t(101325.0) / p_in_Pa));
  Dv = (Dv / real_t(10000.0));

  // Calculate saturation mixing ratio and saturation vapor pressure.
  sat_vapor_press_liq(T_in_K, saturation_formula, esatv);

  // The values of F_k and F_d are found in Rogers and Yau (1989);
  // Eq. 7.17 and 7.18.
  Fk = ((((Lv / (Rv * T_in_K)) - one) * (Lv * rho_lw)) / (Ka * T_in_K));
  Fd = (((rho_lw * Rv) * T_in_K) / (Dv * esatv));

  // Calculate G(T,p).
  G_T_p = (one / (Fk + Fd));

  return G_T_p;
}
} // namespace clubb_kokkos::kk_engine
