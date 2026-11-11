// Compute surface fluxes for the shear-driven neutral-boundary-layer case.
//
// neutral_case_sfclyr prescribes heat flux 0.05 K m/s through 80880 s and zero
// thereafter, with no water flux and friction velocity 0.5 m/s. Supplied winds
// and wind speed [m/s] determine output momentum fluxes [(m/s)^2]. Fields hold
// one value per column; no persistent state is stored.
// Reference: Moeng and Sullivan (1994), following George Bryan's CM1 case.
#pragma once

#include "sfc_flux.hpp"

namespace clubb_kokkos {

//----------------------------------------------------------------------
inline void neutral_case_sfclyr(int ngrdcol, double time,
                                Real1D um_sfc, Real1D vm_sfc, Real1D ubar,
                                Real1D upwp_sfc, Real1D vpwp_sfc,
                                Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D ustar) {
  // Description:
  // Compute surface fluxes for the shear-driven neutral case. Preserve the
  // actual Native switch at 80880 s (the source's 3000-s comment is stale).
  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "neutral_case_sfclyr", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        // Turn off the prescribed heat flux after the Native cutoff time.
        if (time > real_t(80880.0))
          wpthlp_sfc(i) = real_t(0.0);
        else
          wpthlp_sfc(i) = real_t(0.05);

        // No moisture in this case.
        wprtp_sfc(i) = real_t(0.0);

        // Compute ustar from the prescribed constant.
        ustar(i) = real_t(0.5);
      });

  compute_momentum_flux(ngrdcol, um_sfc, vm_sfc, ubar, ustar,
                         upwp_sfc, vpwp_sfc);
}

} // namespace clubb_kokkos
