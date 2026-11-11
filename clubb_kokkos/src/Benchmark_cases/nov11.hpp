// Adjust above-cloud water and interpolate heat input for the November 11 case.
//
// nov11_altocu_rtm_adjust multiplies water above 2900 m relative to the bottom
// boundary by 0.89 at the first step at or after one hour. It preserves the timestep
// window when 3600 s is not divisible by dt. nov11_altocu_read_t_dependent returns
// supplied sensible/latent heat [W/m^2] at time [s]. Model rtm uses (column, level)
// [kg/kg]; the appended context owns the surface series. The bottom is not the
// surface, so there is no surface-layer routine. The obsolete Native tendency
// routine is entirely commented out and remains absent here.
// Reference: Larson et al. (2006), JGR 111 D19207, doi:10.1029/2005JD007002.
#pragma once

#include "CLUBB_core/grid_class.hpp"
#include "CLUBB_core/interpolation.hpp"
#include "time_dependent_input.hpp"

namespace clubb_kokkos {

//-----------------------------------------------------------------------
inline void nov11_altocu_rtm_adjust(int ngrdcol, const grid& gr,
                                    double time, double time_initial, real_t dt,
                                    Real2D rtm) {
  // Description:
  // Make the one-time 0.89 adjustment to total water above cloud one hour after
  // initialization. This replaces the same active block from the obsolete
  // nov11_altocu_tndcy, preserving the original timestep window.
  if (time >= time_initial + double(3600.0) &&
      time < time_initial + double(3600.0) + static_cast<double>(dt)) {
    const int nzt = gr.nzt;
    const auto zt = gr.zt;
    const auto zm = gr.zm;
    using Policy  = Kokkos::RangePolicy<typename Real2D::execution_space>;
    Kokkos::parallel_for(
        "nov11_altocu_rtm_adjust", Policy(0, ngrdcol * nzt), KOKKOS_LAMBDA(int index) {
          const int i = index % ngrdcol;
          const int k = index / ngrdcol;
          if (zt(i, k) > real_t(2900.0) + zm(i, 0))
            rtm(i, k) = real_t(0.89) * rtm(i, k);
        });
  }
}

//-----------------------------------------------------------------------
inline void nov11_altocu_read_t_dependent(ScratchWorkspace& workspace, double time, real_t& sens_ht,
                                          real_t& latent_ht,
                                          const time_dependent_input_state& input) {
  // Description:
  // Interpolate the case's sensible/latent heat input at the requested time
  // and return both scalar fluxes in W/m^2. Bracketing errors follow time_select.
  auto workspace_scope = workspace.scope();
  real_t time_frac;
  int before_time, after_time;

  // Interpolate the values read from time_dependent_input.
  time_select(time, static_cast<int>(input.time_sfc_given.extent(0)),
              input.time_sfc_given, before_time, after_time, time_frac);
  const auto sensible = to_host_space(workspace, input.sens_ht_given);
  const auto latent   = to_host_space(workspace, input.latent_ht_given);
  sens_ht   = linear_interp_factor(time_frac, sensible(after_time), sensible(before_time));
  latent_ht = linear_interp_factor(time_frac, latent(after_time), latent(before_time));
}

} // namespace clubb_kokkos
