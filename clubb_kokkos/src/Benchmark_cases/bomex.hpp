// Set large-scale tendencies and surface fluxes for the GCSS BOMEX case.
//
// bomex_tndcy prescribes moisture advection on the thermodynamic grid and copies
// theta_l/rt forcing into selected scalar fields. bomex_sfclyr interpolates
// heat/specific-humidity surface fluxes, converts the water flux to mixing-ratio
// units and prescribes ustar=0.28 m/s. Model fields use (column, level[,scalar])
// LayoutLeft Views. Tendencies are [K/s]/[(kg/kg)/s]; surface fluxes are
// [K m/s]/[(kg/kg) m/s]. Only the surface routine needs the appended case-owned
// time-series context; no state or scratch persists in these procedures.
// Reference: GCSS BOMEX specification,
// http://www.knmi.nl/~siebesma/gcss/bomexcomp.init.html
#pragma once

#include "scratch_workspace.hpp"
#include "CLUBB_core/array_index.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "CLUBB_core/interpolation.hpp"
#include "spec_hum_to_mixing_ratio.hpp"
#include "time_dependent_input.hpp"

namespace clubb_kokkos {

//----------------------------------------------------------------------
inline void bomex_tndcy(ScratchWorkspace& workspace, int ngrdcol, int sclr_dim, int edsclr_dim,
                        const sclr_idx_type& sclr_idx, const grid& gr, Real2D rtm,
                        bool l_readiopdata, Real2D thlm_forcing, Real2D rtm_forcing,
                        Real3D sclrm_forcing, Real3D edsclrm_forcing) {
  // Description:
  // Set theta_l and water tendencies for BOMEX, converting the specified
  // specific-humidity advection to mixing-ratio forcing and copying scalars.
  // The Native readiopdata block is commented out; its formal remains unused.
  auto workspace_scope = workspace.scope();
  (void)l_readiopdata;
  const int nzt = gr.nzt;
  const auto zt = gr.zt;
  Real2D qtm_forcing = workspace_scope.view<Real2D>("bomex_tndcy::qtm_forcing", ngrdcol, nzt);
  using Policy = Kokkos::RangePolicy<typename Real2D::execution_space>;

  // Large-scale advective moisture tendency. BOMEX specifies total-water
  // specific humidity rather than mixing ratio.
  Kokkos::parallel_for(
      "bomex_tndcy::profiles", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
        const int i = ((index % ngrdcol));
        const int k = ((index / ngrdcol));
        if (zt(i, k) >= real_t(0.0) && zt(i, k) < real_t(300.0)) {
          qtm_forcing(i, k) = -real_t(1.2e-8);
        } else if (zt(i, k) >= real_t(300.0) && zt(i, k) < real_t(500.0)) {
          qtm_forcing(i, k) = -real_t(1.2e-8) * (real_t(1.0) - (zt(i, k) - real_t(300.0)) /
                                                                   (real_t(500.0) - real_t(300.0)));
        } else {
          qtm_forcing(i, k) = real_t(0.0);
        }

        // Radiative theta-l tendency.
        thlm_forcing(i, k) = real_t(0.0);
      });

  // Convert total-water specific-humidity forcing to mixing-ratio forcing.
  force_spec_hum_to_mixing_ratio(ngrdcol, nzt, rtm, qtm_forcing, rtm_forcing);

  // Test scalars with theta_l and rt if desired.
  if (sclr_dim > 0) {
    Kokkos::parallel_for(
        "bomex_tndcy::passive", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
          const int i = ((index % ngrdcol));
          const int k = ((index / ngrdcol));
          if (sclr_idx.iisclr_thl > -1)
            sclrm_forcing(i, k, sclr_idx.iisclr_thl) = thlm_forcing(i, k);
          if (sclr_idx.iisclr_rt > -1)
            sclrm_forcing(i, k, sclr_idx.iisclr_rt) = rtm_forcing(i, k);
        });
  }
  if (edsclr_dim > 0) {
    Kokkos::parallel_for(
        "bomex_tndcy::eddy", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
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
inline void bomex_sfclyr(ScratchWorkspace& workspace, int ngrdcol, double time, Real1D rtm_sfc,
                         Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D ustar,
                         const time_dependent_input_state& input) {
  // Description:
  // Compute GCSS BOMEX surface heat/water fluxes and friction velocity. The
  // supplied water flux is specific humidity, so convert it using rtm_sfc.
  auto workspace_scope = workspace.scope();
  real_t time_frac;
  int before_time, after_time;

  // Compute heat and moisture fluxes.
  time_select(time, static_cast<int>(input.time_sfc_given.extent(0)),
              input.time_sfc_given, before_time, after_time, time_frac);
  const auto heat = to_host_space(workspace, input.wpthlp_sfc_given);
  const real_t wpthlp_sfc_calc =
      linear_interp_factor(time_frac, heat(after_time), heat(before_time));

  // BOMEX specifies the surface water flux in total-water specific humidity.
  const auto water = to_host_space(workspace, input.wpqtp_sfc_given);
  const real_t wpqtp_sfc_calc =
      linear_interp_factor(time_frac, water(after_time), water(before_time));

  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "bomex_sfclyr", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        // Declare the prescribed value of ustar.
        ustar(i)      = real_t(0.28);
        wpthlp_sfc(i) = wpthlp_sfc_calc;
        // Convert specific-humidity flux to total-water mixing-ratio flux.
        const real_t factor = real_t(1.0) + rtm_sfc(i);
        wprtp_sfc(i)        = factor * factor * wpqtp_sfc_calc;
      });
}

} // namespace clubb_kokkos
