// Compute surface fluxes for Hing Ong's neutral Ekman-layer case.
//
// ekman_sfclyr sets heat/moisture fluxes to zero, prescribes friction velocity
// 0.3 m/s, and computes momentum fluxes from the supplied winds and wind speed
// [m/s]. Output momentum fluxes have units [(m/s)^2]. All fields contain one
// value per column; this routine owns no persistent state.
#pragma once

#include "sfc_flux.hpp"

namespace clubb_kokkos {

//----------------------------------------------------------------------
inline void ekman_sfclyr(int ngrdcol,
                         Real1D um_sfc, Real1D vm_sfc, Real1D ubar,
                         Real1D upwp_sfc, Real1D vpwp_sfc,
                         Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D ustar) {
  // Description:
  // Compute the momentum, heat and moisture fluxes for a neutral Ekman layer.
  // The prescribed friction velocity supplies the momentum-flux calculation.
  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "ekman_sfclyr", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        // No heat and moisture fluxes.
        wpthlp_sfc(i) = real_t(0.0);
        wprtp_sfc(i)  = real_t(0.0);

        // Set ustar.
        ustar(i) = real_t(0.3);
      });

  compute_momentum_flux(ngrdcol, um_sfc, vm_sfc, ubar, ustar,
                         upwp_sfc, vpwp_sfc);
}

} // namespace clubb_kokkos
