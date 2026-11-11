// Compute surface fluxes and convert supplied heat fluxes to kinematic units.
//
// compute_momentum_flux and compute_ubar diagnose momentum flux and wind speed.
// compute_ht_mostr_flux interpolates the supplied time series of sensible/latent
// heat flux [W/m^2]. compute_wpthlp_sfc and compute_wprtp_sfc use bulk transfer
// coefficients to calculate heat [K m/s] and water [(kg/kg) m/s] fluxes.
// set_sclr_sfc_rtm_thlm supplies the corresponding passive/eddy-scalar fluxes.
// convert_sens_ht_to_km_s and convert_latent_ht_to_m_s convert scalar heat fluxes
// [W/m^2] using surface air density [kg/m^3]. Views are column-first; no persistent
// state is stored here. Interpolating scalar host results copies only its input
// time series to host memory.
#pragma once

#include "CLUBB_core/array_index.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/interpolation.hpp"

namespace clubb_kokkos {

//===============================================================================
inline void compute_momentum_flux(int ngrdcol, Real1D um_sfc, Real1D vm_sfc,
                                  Real1D ubar, Real1D ustar,
                                  Real1D upwp_sfc, Real1D vpwp_sfc) {
  // Description:
  // Compute upward turbulent momentum fluxes [(m/s)^2] from winds, wind speed
  // and surface friction velocity [m/s], one value per column.
  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;

  // Compute momentum fluxes.
  Kokkos::parallel_for(
      "compute_momentum_flux", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        upwp_sfc(i) = -um_sfc(i) * (ustar(i) * ustar(i)) / ubar(i);
        vpwp_sfc(i) = -vm_sfc(i) * (ustar(i) * ustar(i)) / ubar(i);
      });
}

//===============================================================================
inline void compute_ubar(int ngrdcol, Real1D um_sfc, Real1D vm_sfc, Real1D ubar) {
  // Description:
  // Compute the surface wind speed [m/s] from the wind at zt(2), subject to
  // the Native 0.25 m/s minimum speed used by bulk surface exchange.
  constexpr real_t ubmin = real_t(0.25); // Minimum speed [m/s].
  using Policy           = Kokkos::RangePolicy<typename Real1D::execution_space>;

  Kokkos::parallel_for(
      "compute_ubar", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        ubar(i) = Kokkos::max(ubmin, Kokkos::sqrt(um_sfc(i) * um_sfc(i) + vm_sfc(i) * vm_sfc(i)));
      });
}

//===============================================================================
inline void
compute_ht_mostr_flux(ScratchWorkspace& workspace, double time_in, int ntimes,
                      Real1D latent_ht_given, Real1D sens_ht_given,
                      Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> time_sfc_given,
                      real_t& heat_flx, real_t& moisture_flx) {
  // Description:
  // Interpolate ARM heat/moisture fluxes [W/m^2] at time_in [s], holding the
  // first or last supplied flux outside the available times. Return scalars.
  auto workspace_scope      = workspace.scope();
  const auto latent_ht_host = to_host_space(workspace, latent_ht_given);
  const auto sens_ht_host   = to_host_space(workspace, sens_ht_given);
  const auto time_sfc_host  = time_sfc_given;

  // Default initialization.
  heat_flx          = real_t(0.0);
  moisture_flx      = real_t(0.0);
  const real_t time = static_cast<real_t>(time_in);

  // Compute heat and moisture fluxes from ARM data in W/m^2.
  if (time <= time_sfc_host(0)) {
    heat_flx     = sens_ht_host(0);
    moisture_flx = latent_ht_host(0);
  } else if (time >= time_sfc_host(ntimes - 1)) {
    heat_flx     = sens_ht_host(ntimes - 1);
    moisture_flx = latent_ht_host(ntimes - 1);
  } else {
    int i1 = 0;
    while (i1 < ntimes - 1) {
      const int i2 = i1 + 1;
      if (time >= time_sfc_host(i1) && time < time_sfc_host(i2)) {
        const real_t time_frac =
            (time - time_sfc_host(i1)) / (time_sfc_host(i2) - time_sfc_host(i1));
        heat_flx     = linear_interp_factor(time_frac, sens_ht_host(i2), sens_ht_host(i1));
        moisture_flx = linear_interp_factor(time_frac, latent_ht_host(i2), latent_ht_host(i1));
        i1 = ntimes;
      }
      i1 = i2;
    }
  }
}

//===============================================================================
inline void compute_wpthlp_sfc(int ngrdcol, Real1D Cd, Real1D ubar, Real1D thlm_sfc,
                               Real1D T_sfc, Real1D exner_sfc, Real1D wpthlp_sfc) {
  // Description:
  // Determine surface heat flux [K m/s] from a transfer coefficient, wind
  // speed [m/s], near-surface theta_l [K], temperature [K] and surface Exner.
  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "compute_wpthlp_sfc", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        wpthlp_sfc(i) = -Cd(i) * ubar(i) * (thlm_sfc(i) - T_sfc(i) / exner_sfc(i));
      });
}

//===============================================================================
inline void compute_wprtp_sfc(int ngrdcol, Real1D Cd, Real1D ubar, Real1D rtm_sfc,
                              Real1D adjustment, Real1D wprtp_sfc) {
  // Description:
  // Determine surface moisture flux [(kg/kg) m/s] using a bulk coefficient,
  // wind speed [m/s], near-surface water ratio and its adjustment [kg/kg].
  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "compute_wprtp_sfc", Policy(0, ngrdcol),
      KOKKOS_LAMBDA(int i) { wprtp_sfc(i) = -Cd(i) * ubar(i) * (rtm_sfc(i) - adjustment(i)); });
}

//===============================================================================
inline void set_sclr_sfc_rtm_thlm(int ngrdcol, int sclr_dim, int edsclr_dim,
                                  const sclr_idx_type& sclr_idx,
                                  Real1D wpthlp_sfc, Real1D wprtp_sfc,
                                  Real2D wpsclrp_sfc, Real2D wpedsclrp_sfc) {
  // Description:
  // Initialize passive/eddy-scalar surface fluxes and assign the theta_l/rt
  // passive copies. Preserve Native's final eddy-scalar zeroing and its order.
  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  if (sclr_dim > 0) {
    Kokkos::parallel_for(
        "set_sclr_sfc_rtm_thlm::zero_passive", Policy(0, ((ngrdcol * sclr_dim))),
        KOKKOS_LAMBDA(int index) {
          const int i          = ((index % ngrdcol));
          const int sclr       = ((index / ngrdcol));
          wpsclrp_sfc(i, sclr) = real_t(0.0);
        });
    Kokkos::parallel_for(
        "set_sclr_sfc_rtm_thlm::copies", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
          // Let passive scalars be equal to rt and theta_l for now.
          if (sclr_idx.iisclr_thl > -1)
            wpsclrp_sfc(i, sclr_idx.iisclr_thl) = wpthlp_sfc(i);
          if (sclr_idx.iisclr_rt > -1)
            wpsclrp_sfc(i, sclr_idx.iisclr_rt) = wprtp_sfc(i);
          if (sclr_idx.iiedsclr_thl > -1)
            wpedsclrp_sfc(i, sclr_idx.iiedsclr_thl) = wpthlp_sfc(i);
          if (sclr_idx.iiedsclr_rt > -1)
            wpedsclrp_sfc(i, sclr_idx.iiedsclr_rt) = wprtp_sfc(i);
        });
  }
  if (edsclr_dim > 0) {
    Kokkos::parallel_for(
        "set_sclr_sfc_rtm_thlm::zero_eddy", Policy(0, ((edsclr_dim * ngrdcol))),
        KOKKOS_LAMBDA(int index) {
          const int i              = ((index % ngrdcol));
          const int edsclr         = ((index / ngrdcol));
          wpedsclrp_sfc(i, edsclr) = real_t(0.0);
        });
  }
}

//==============================================================================
KOKKOS_INLINE_FUNCTION real_t convert_sens_ht_to_km_s(real_t sens_ht, real_t rho_sfc) {
  // Description:
  // Convert sensible heat flux [W/m^2] to the kinematic heat flux [K m/s]
  // using the surface density [kg/m^3] and specific heat Cp [J/(kg K)].
  return sens_ht / (rho_sfc * Cp);
}

//==============================================================================
KOKKOS_INLINE_FUNCTION real_t convert_latent_ht_to_m_s(real_t latent_ht, real_t rho_sfc) {
  // Description:
  // Convert latent heat flux [W/m^2] to water flux [(kg/kg) m/s] using the
  // surface density [kg/m^3] and latent heat of vaporization Lv [J/kg].
  return latent_ht / (rho_sfc * Lv);
}

} // namespace clubb_kokkos
