// Supply zero large-scale tendencies, top subsidence and surface exchange for DYCOMS II RF02.
//
// dycoms2_rf02_tndcy fills theta_l [K/s] and water [(kg/kg)/s] tendencies and
// selected scalar copies on column-first thermodynamic fields. dycoms2_rf02_sfclyr
// supplies heat [K m/s], water [(kg/kg) m/s] and friction velocity [m/s].
// The tendency routine also zeroes vertical wind at the model top on both
// grids. Surface series belong to the appended case context; density is 1.21 kg/m^3.
// Reference: Wyant et al. (2007), JGR 112 D24204; Stevens et al. (2003), BAMS 84,579–593.
#pragma once

#include "CLUBB_core/array_index.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "sfc_flux.hpp"
#include "time_dependent_input.hpp"

namespace clubb_kokkos {

//----------------------------------------------------------------------
inline void dycoms2_rf02_tndcy(int ngrdcol, int sclr_dim, int edsclr_dim, const sclr_idx_type& sclr_idx,
                               const grid& gr, Real2D wm_zt, Real2D wm_zm,
                               Real2D thlm_forcing, Real2D rtm_forcing,
                               Real3D sclrm_forcing, Real3D edsclrm_forcing) {
  // Description:
  // Set the case's large-scale theta_l and water tendencies and corresponding
  // passive/eddy scalar copies. Other scalar components remain unchanged.
  const int nzt = gr.nzt;
  using Policy  = Kokkos::RangePolicy<typename Real2D::execution_space>;
  const int nzm = gr.nzm;

  // Zero top subsidence to avoid importing values above the model domain with
  // the one-sided thermodynamic advection derivative; keep both grids consistent.
  Kokkos::parallel_for(
      "dycoms2_rf02_tndcy::top", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        wm_zt(i, ((nzt - 1))) = real_t(0.0);
        wm_zm(i, ((nzm - 1))) = real_t(0.0);
      });

  Kokkos::parallel_for(
      "dycoms2_rf02_tndcy::profiles", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
        const int i = ((index % ngrdcol));
        const int k = ((index / ngrdcol));
        // Large-scale temperature and advective moisture tendencies.
        thlm_forcing(i, k) = real_t(0.0);
        rtm_forcing(i, k)  = real_t(0.0);
      });
  if (sclr_dim > 0) {
    Kokkos::parallel_for(
        "dycoms2_rf02_tndcy::passive", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
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
        "dycoms2_rf02_tndcy::eddy", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
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
inline void dycoms2_rf02_sfclyr(ScratchWorkspace& workspace, int ngrdcol, double time,
                                Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D ustar,
                                const time_dependent_input_state& input) {
  // Description:
  // Compute GCSS DYCOMS II RF02 heat and moisture fluxes using the supplied
  // physical fluxes and prescribed surface air density [kg/m^3].
  auto workspace_scope     = workspace.scope();
  constexpr real_t rho_sfc = real_t(1.21);
  real_t time_frac;
  int before_time, after_time;
  time_select(time, static_cast<int>(input.time_sfc_given.extent(0)),
              input.time_sfc_given, before_time, after_time, time_frac);
  const auto sensible = to_host_space(workspace, input.sens_ht_given);
  const real_t sens_ht =
      linear_interp_factor(time_frac, sensible(after_time), sensible(before_time));
  const auto latent      = to_host_space(workspace, input.latent_ht_given);
  const real_t latent_ht = linear_interp_factor(time_frac, latent(after_time), latent(before_time));

  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "dycoms2_rf02_sfclyr", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        // Declare the prescribed value of ustar.
        ustar(i)      = real_t(0.25);
        wpthlp_sfc(i) = convert_sens_ht_to_km_s(sens_ht, rho_sfc);
        wprtp_sfc(i)  = convert_latent_ht_to_m_s(latent_ht, rho_sfc);
      });
}

} // namespace clubb_kokkos
