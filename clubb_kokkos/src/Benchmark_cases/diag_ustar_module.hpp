// Diagnose friction velocity from wind, roughness and surface buoyancy flux.
//
// diag_ustar is a scalar, device-callable similarity calculation. Supply the
// wind speed [m/s] at height z [m], momentum roughness z0 [m], and kinematic
// buoyancy flux bflx [m^2/s^3]; it returns friction velocity [m/s]. Four
// iterations apply for non-negligible flux, with the Native large-stability
// escape retained. The original code cautions that it was not thoroughly tested.
//
// Reference: Businger (1973), Turbulent Transfer in the Atmospheric Surface
// Layer, Workshop on Micrometeorology, pp. 67–100. Bjorn Stevens, March 1999.
#pragma once

#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/constants_clubb.hpp"

namespace clubb_kokkos {

// ----------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION real_t diag_ustar(real_t z, real_t bflx, real_t wnd, real_t z0) {
  // Description:
  // Return ustar using phi_m = 1 + am*zeta for stable conditions and
  // phi_m = (1 - bm*zeta)^(-1/4) for unstable conditions, where zeta = z/lmo
  // and lmo is the Obukhov length. The supplied flux is in kinematic units.
  constexpr real_t am = real_t(4.8);
  constexpr real_t bm = real_t(19.3);

  const real_t lnz  = Kokkos::log(z / z0);
  const real_t klnz = vonk / lnz;
  const real_t c1   = pi / real_t(2.0) - real_t(3.0) * Kokkos::log(real_t(2.0));

  real_t ustar = wnd * klnz;
  if (Kokkos::abs(bflx) > real_t(1.e-6)) {
    for (int iterate = 1; iterate <= 4; ++iterate) {
      const real_t lmo  = -(ustar * ustar * ustar) / (vonk * bflx);
      const real_t zeta = z / lmo;
      if (zeta > real_t(0.0)) {
        if (zeta > real_t(1.e10)) { // Native large-zeta escape, dschanen/UWM.
          ustar = real_t(1.e-10);
          break;
        } else {
          ustar = vonk * wnd / (lnz + am * zeta);
        }
      } else {
        const real_t x    = Kokkos::sqrt(Kokkos::sqrt(real_t(1.0) - bm * zeta));
        const real_t psi1 = real_t(2.0) * Kokkos::log(real_t(1.0) + x) +
                            Kokkos::log(real_t(1.0) + x * x) - real_t(2.0) * Kokkos::atan(x) + c1;
        ustar             = wnd * vonk / (lnz - psi1);
      }
    }
  }

  return ustar;
}

} // namespace clubb_kokkos
