// Compute bulk surface exchange using the GABLS3 vegetation temperature.
//
// gabls3_sfclyr supplies heat [K m/s], water[(kg/kg) m/s] and friction velocity
// [m/s] from wind [m/s], vegetation temperature [K], near-surface theta_l [K],
// water [kg/kg], lowest height [m] and Exner [-]. It retains the specified tenfold
// water-flux adjustment. Views hold one value per column; no state persists.
// Reference: GABLS3 specification, http://www.knmi.nl/samenw/gabls/
#pragma once

#include "scratch_workspace.hpp"
#include "diag_ustar_module.hpp"
#include "sfc_flux.hpp"

namespace clubb_kokkos {

//-----------------------------------------------------------------------
inline void gabls3_sfclyr(ScratchWorkspace& workspace, int ngrdcol, Real1D ubar, Real1D veg_T_in_K,
                          Real1D thlm_sfc, Real1D rtm_sfc, Real1D lowest_level, Real1D exner_sfc,
                          Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D ustar) {
  // Description:
  // Compute GABLS3 surface heat/water fluxes from vegetation temperature,
  // then diagnose momentum exchange from the resulting buoyancy flux.
  auto workspace_scope = workspace.scope();
  Real1D offset        = workspace_scope.view<Real1D>("gabls3_sfclyr::offset", ngrdcol);
  Real1D C_10          = workspace_scope.view<Real1D>("gabls3_sfclyr::C_10", ngrdcol);
  constexpr real_t z0 = real_t(0.15); // Momentum roughness height [m].
  using Policy        = Kokkos::RangePolicy<typename Real1D::execution_space>;

  // Compute heat and moisture fluxes.
  Kokkos::parallel_for(
      "gabls3_sfclyr::coefficients", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        offset(i) = real_t(9.9e-3);
        C_10(i)   = real_t(0.00195);
      });
  compute_wpthlp_sfc(ngrdcol, C_10, ubar, thlm_sfc, veg_T_in_K, exner_sfc, wpthlp_sfc);
  compute_wprtp_sfc(ngrdcol, C_10, ubar, rtm_sfc, offset, wprtp_sfc);

  // Compute momentum fluxes, preserving the Native water-flux adjustment.
  Kokkos::parallel_for(
      "gabls3_sfclyr::momentum", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        wprtp_sfc(i)                = wprtp_sfc(i) * real_t(10.0);
        const real_t veg_theta_in_K = veg_T_in_K(i) / exner_sfc(i);
        const real_t bflx           = wpthlp_sfc(i) * grav / veg_theta_in_K;
        ustar(i)                    = diag_ustar(lowest_level(i), bflx, ubar(i), z0);
      });
}

} // namespace clubb_kokkos
