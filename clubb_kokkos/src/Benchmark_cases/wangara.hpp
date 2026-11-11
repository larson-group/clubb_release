// Supply large-scale and surface forcing for Wangara day 33.
//
// wangara_tndcy zeros subsidence [m/s] and heat/water tendencies [K/s,(kg/kg)/s]
// and copies selected scalars. wangara_sfclyr computes diurnal heat [K m/s],
// water[(kg/kg) m/s] and friction velocity [m/s] from UTC time [s]. It accepts
// local Australian EST times only from 07:30 through 17:30. Fields retain column-
// first LayoutLeft storage; no persistent state is stored.
// Reference: Golaz et al. (2002), A PDF-Based Model for Boundary Layer Clouds,
// Part II: Model Results, JAS 59,3552–3571.
#pragma once

#include "CLUBB_core/array_index.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include <iostream>
#include <stdexcept>

namespace clubb_kokkos {

//----------------------------------------------------------------------
inline void wangara_tndcy(int ngrdcol, const grid& gr, int sclr_dim, int edsclr_dim, const sclr_idx_type& sclr_idx,
                          Real2D wm_zt, Real2D wm_zm,
                          Real2D thlm_forcing, Real2D rtm_forcing,
                          Real3D sclrm_forcing, Real3D edsclrm_forcing) {
  // Description:
  // Set zero theta/water tendencies and subsidence for the Wangara case,
  // then copy the requested theta_l/rt passive and eddy-scalar tendencies.
  const int nzt = gr.nzt;
  using Policy  = Kokkos::RangePolicy<typename Real2D::execution_space>;

  // No large-scale subsidence on the momentum grid.
  Kokkos::deep_copy(wm_zm, real_t(0.0));
  Kokkos::parallel_for(
      "wangara_tndcy::profiles", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
        const int i = ((index % ngrdcol));
        const int k = ((index / ngrdcol));
        // No large-scale subsidence, water tendency or cooling.
        wm_zt(i, k)        = real_t(0.0);
        rtm_forcing(i, k)  = real_t(0.0);
        thlm_forcing(i, k) = real_t(0.0);
      });
  if (sclr_dim > 0) {
    Kokkos::parallel_for(
        "wangara_tndcy::passive", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
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
        "wangara_tndcy::eddy", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
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
inline void wangara_sfclyr(int ngrdcol, double time,
                           Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D ustar) {
  // Description:
  // Compute Wangara day 33 heat/moisture fluxes using local Australian EST
  // time, rejecting times outside 07:30–17:30. Return prescribed ustar too.
  const double time_utc   = Kokkos::fmod(time, static_cast<double>(sec_per_day));
  const double est_offset = 36000.0;
  const double time_est   = Kokkos::fmod(time_utc + est_offset, static_cast<double>(sec_per_day));
  if (time_est < 27000.0 || time_est > 63000.0) {
    std::cerr << "wangara_sfclyr: error local time must be between 730 and 1730._core_rknd\n";
    std::cerr << "time_est = " << time_est << '\n';
    throw std::runtime_error("");
  }
  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "wangara_sfclyr", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        // Declare the prescribed value of ustar.
        ustar(i) = real_t(0.13);

        // Compute heat and moisture fluxes from the specified diurnal coefficients.
        wpthlp_sfc(i) =
            real_t(0.18) *
            Kokkos::cos((static_cast<real_t>(time_est) - real_t(45000.0)) / real_t(36000.0) * pi);
        wprtp_sfc(i) = real_t(1.3e-4) * wpthlp_sfc(i);
      });
}

} // namespace clubb_kokkos
