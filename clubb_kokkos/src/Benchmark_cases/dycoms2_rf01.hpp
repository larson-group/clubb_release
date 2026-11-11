// Supply zero large-scale tendencies and fixed/bulk surface exchange for DYCOMS II RF01.
//
// dycoms2_rf01_tndcy fills theta_l [K/s] and water [(kg/kg)/s] tendencies and
// selected scalar copies on column-first thermodynamic fields. dycoms2_rf01_sfclyr
// supplies heat [K m/s], water [(kg/kg) m/s] and friction velocity [m/s].
// Surface temperature [K] is also returned; sfctype selects supplied fluxes (0)
// or bulk coefficients (1). Surface series belong to the appended case context.
// Reference: GCSS DYCOMS II RF01 specification, http://www.atmos.ucla.edu/~bstevens/dycoms/rf01/rf01.html.
#pragma once

#include "scratch_workspace.hpp"
#include "CLUBB_core/array_index.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "sfc_flux.hpp"
#include "time_dependent_input.hpp"
#include "CLUBB_core/saturation.hpp"
#include <iostream>

namespace clubb_kokkos {

//----------------------------------------------------------------------
inline void dycoms2_rf01_tndcy(int ngrdcol, const grid& gr, int sclr_dim, int edsclr_dim, const sclr_idx_type& sclr_idx,
                               Real2D thlm_forcing, Real2D rtm_forcing,
                               Real3D sclrm_forcing, Real3D edsclrm_forcing) {
  // Description:
  // Set the case's large-scale theta_l and water tendencies and corresponding
  // passive/eddy scalar copies. Other scalar components remain unchanged.
  const int nzt = gr.nzt;
  using Policy  = Kokkos::RangePolicy<typename Real2D::execution_space>;

  Kokkos::parallel_for(
      "dycoms2_rf01_tndcy::profiles", Policy(0, (ngrdcol * nzt)), KOKKOS_LAMBDA(int index) {
        const int i = (index % ngrdcol);
        const int k = (index / ngrdcol);
        // Large-scale temperature and advective moisture tendencies.
        thlm_forcing(i, k) = real_t(0.0);
        rtm_forcing(i, k)  = real_t(0.0);
      });
  if (sclr_dim > 0) {
    Kokkos::parallel_for(
        "dycoms2_rf01_tndcy::passive", Policy(0, (ngrdcol * nzt)), KOKKOS_LAMBDA(int index) {
          const int i = (index % ngrdcol);
          const int k = (index / ngrdcol);
          // Test scalars with theta_l and rt if desired.
          if (sclr_idx.iisclr_thl > -1)
            sclrm_forcing(i, k, sclr_idx.iisclr_thl) = thlm_forcing(i, k);
          if (sclr_idx.iisclr_rt > -1)
            sclrm_forcing(i, k, sclr_idx.iisclr_rt) = rtm_forcing(i, k);
        });
  }
  if (edsclr_dim > 0) {
    Kokkos::parallel_for(
        "dycoms2_rf01_tndcy::eddy", Policy(0, (ngrdcol * nzt)), KOKKOS_LAMBDA(int index) {
          const int i = (index % ngrdcol);
          const int k = (index / ngrdcol);
          if (sclr_idx.iiedsclr_thl > -1)
            edsclrm_forcing(i, k, sclr_idx.iiedsclr_thl) = thlm_forcing(i, k);
          if (sclr_idx.iiedsclr_rt > -1)
            edsclrm_forcing(i, k, sclr_idx.iiedsclr_rt) = rtm_forcing(i, k);
        });
  }
}

//======================================================================
inline void dycoms2_rf01_sfclyr(ScratchWorkspace& workspace, int ngrdcol, double time, int sfctype,
                                Real1D p_sfc, Real1D exner_sfc, Real1D ubar, Real1D thlm_sfc,
                                Real1D rtm_sfc, Real1D rho_sfc, int saturation_formula,
                                Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D ustar, Real1D T_sfc,
                                const time_dependent_input_state& input) {
  // Description:
  // Compute GCSS DYCOMS II RF01 heat/water surface fluxes using the supplied
  // physical fluxes or the bulk exchange coefficient, according to sfctype.
  auto workspace_scope = workspace.scope();
  real_t time_frac;
  int before_time, after_time;
  time_select(time, static_cast<int>(input.time_sfc_given.extent(0)),
              input.time_sfc_given, before_time, after_time, time_frac);
  const auto temperature = to_host_space(workspace, input.T_sfc_given);
  const real_t T_sfc_interp =
      linear_interp_factor(time_frac, temperature(after_time), temperature(before_time));

  // Compute heat and moisture fluxes.
  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  if (sfctype == 0) {
    const auto sensible = to_host_space(workspace, input.sens_ht_given);
    const real_t sens_ht =
        linear_interp_factor(time_frac, sensible(after_time), sensible(before_time));
    const auto latent = to_host_space(workspace, input.latent_ht_given);
    const real_t latent_ht =
        linear_interp_factor(time_frac, latent(after_time), latent(before_time));
    Kokkos::parallel_for(
        "dycoms2_rf01_sfclyr::fixed", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
          ustar(i)      = real_t(0.25);
          T_sfc(i)      = T_sfc_interp;
          wpthlp_sfc(i) = convert_sens_ht_to_km_s(sens_ht, rho_sfc(i));
          wprtp_sfc(i)  = convert_latent_ht_to_m_s(latent_ht, rho_sfc(i));
        });
  } else if (sfctype == 1) {
    Kokkos::parallel_for(
        "dycoms2_rf01_sfclyr::bulk", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
          ustar(i) = real_t(0.25);
          T_sfc(i) = T_sfc_interp;
          const real_t Cd   = real_t(0.0011);
          const real_t rsat = sat_mixrat_liq_api(p_sfc(i), T_sfc(i), saturation_formula);
          // Consume the column-local coefficients without intermediate fields.
          wpthlp_sfc(i) = -Cd * ubar(i) * (thlm_sfc(i) - T_sfc(i) / exner_sfc(i));
          wprtp_sfc(i)  = -Cd * ubar(i) * (rtm_sfc(i) - rsat);
        });

  } else {
    std::cerr << "Invalid sfctype value = " << sfctype << '\n';
    throw std::runtime_error("");
  }
}

} // namespace clubb_kokkos
