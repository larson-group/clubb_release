// Supply zero large-scale tendencies and diurnal surface exchange for the LBA case.
//
// lba_tndcy fills theta_l [K/s] and water [(kg/kg)/s] tendencies and
// selected scalar copies on column-first thermodynamic fields. lba_sfclyr
// supplies heat [K m/s], water [(kg/kg) m/s] and friction velocity [m/s].
// The surface modulation uses elapsed time [s] and local density [kg/m^3],
// theta_l [K], height [m] and wind speed [m/s]; no persistent state is stored.
// Reference: Grabowski et al. (2005), GCSS working group 4 case 4 specification.
#pragma once

#include "CLUBB_core/array_index.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "sfc_flux.hpp"
#include "time_dependent_input.hpp"
#include "diag_ustar_module.hpp"

namespace clubb_kokkos {

//----------------------------------------------------------------------
inline void lba_tndcy(int ngrdcol, int sclr_dim, int edsclr_dim, const sclr_idx_type& sclr_idx,
                      const grid& gr,
                      Real2D thlm_forcing, Real2D rtm_forcing,
                      Real3D sclrm_forcing, Real3D edsclrm_forcing) {
  // Description:
  // Set the case's large-scale theta_l and water tendencies and corresponding
  // passive/eddy scalar copies. Other scalar components remain unchanged.
  const int nzt = gr.nzt;
  using Policy  = Kokkos::RangePolicy<typename Real2D::execution_space>;

  Kokkos::parallel_for(
      "lba_tndcy::profiles", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
        const int i = ((index % ngrdcol));
        const int k = ((index / ngrdcol));
        // Large-scale temperature and advective moisture tendencies.
        thlm_forcing(i, k) = real_t(0.0);
        rtm_forcing(i, k)  = real_t(0.0);
      });
  if (sclr_dim > 0) {
    Kokkos::parallel_for(
        "lba_tndcy::passive", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
          const int i = ((index % ngrdcol));
          const int k = ((index / ngrdcol));
          // Test scalars with theta_l and rt if desired.
          if (sclr_idx.iisclr_thl > -1)
            sclrm_forcing(i, k, sclr_idx.iisclr_thl) = thlm_forcing(i, k);
          if (sclr_idx.iisclr_rt > -1)
            sclrm_forcing(i, k, sclr_idx.iisclr_rt) = rtm_forcing(i, k);
        });
  }
  if (edsclr_dim > 0) {
    Kokkos::parallel_for(
        "lba_tndcy::eddy", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
          const int i = ((index % ngrdcol));
          const int k = ((index / ngrdcol));
          if (sclr_idx.iiedsclr_thl > -1)
            edsclrm_forcing(i, k, sclr_idx.iiedsclr_thl) = thlm_forcing(i, k);
          if (sclr_idx.iiedsclr_rt > -1)
            edsclrm_forcing(i, k, sclr_idx.iiedsclr_rt) = rtm_forcing(i, k);
        });
  }
}

//----------------------------------------------------------------------
inline void lba_sfclyr(int ngrdcol, double time_current, double time_initial,
                       Real1D z, Real1D rho_sfc, Real1D thlm_sfc, Real1D ubar,
                       Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D ustar) {
  // Description:
  // Compute LBA surface heat/moisture fluxes using Table A.1 of the case
  // specification, then diagnose friction velocity from the buoyancy flux.
  constexpr real_t z0 = real_t(0.035); // ARM momentum roughness height [m].

  // Compute heat and moisture fluxes from Table A.1.
  const double time = time_current - time_initial;
  const real_t ft   = Kokkos::max(
      real_t(0.0),
      Kokkos::cos(real_t(0.5) * pi *
                  ((real_t(5.25) - static_cast<real_t>(time) / sec_per_hr) / real_t(5.25))));
  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "lba_sfclyr", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        // The coefficients and powers are the specified case constants.
        wpthlp_sfc(i) =
            convert_sens_ht_to_km_s(real_t(270.0) * Kokkos::pow(ft, real_t(1.5)), rho_sfc(i));
        wprtp_sfc(i) =
            convert_latent_ht_to_m_s(real_t(554.0) * Kokkos::pow(ft, real_t(1.3)), rho_sfc(i));
        const real_t bflx = grav / thlm_sfc(i) * wpthlp_sfc(i);

        // Compute ustar.
        ustar(i) = diag_ustar(z(i), bflx, ubar(i), z0);
      });
}

} // namespace clubb_kokkos
