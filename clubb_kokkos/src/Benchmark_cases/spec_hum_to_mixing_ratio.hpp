// Convert moisture fluxes and tendencies from specific humidity to mixing ratio.
//
// flux_spec_hum_to_mixing_ratio converts one momentum-level flux per column.
// force_spec_hum_to_mixing_ratio converts a (column, thermodynamic level) field.
// Inputs are total-water mixing ratio [kg/kg] and specific-humidity flux
// [(kg/kg) m/s] or tendency [(kg/kg)/s]; outputs have the corresponding mixing-
// ratio units. The linearization neglects higher-order fluctuations, as in
// Native. Fields retain column-first LayoutLeft storage and caller ownership.
#pragma once

#include "CLUBB_core/array_types.hpp"

namespace clubb_kokkos {

//=============================================================================
template <class RtmView, class FluxView, class ResultView>
void flux_spec_hum_to_mixing_ratio(int ngrdcol, const RtmView& rtm_zm,
                                   const FluxView& wpqtp, const ResultView& wprtp) {
  // Description:
  // Convert w'q_t' to w'r_t' using (1 + r_tm)^2 * w'q_t'. The relation is
  // linearized by neglecting comparatively small higher-order terms. It is
  // valid at any momentum level, with a surface flux being the usual use.
  const auto rtm_zm_exec = to_execution_space(rtm_zm);
  const auto wpqtp_exec  = to_execution_space(wpqtp);
  const auto wprtp_exec  = to_execution_space(wprtp, false);
  using Policy           = Kokkos::RangePolicy<typename Real1D::execution_space>;

  // Solve for flux in terms of total water mixing ratio.
  Kokkos::parallel_for(
      "flux_spec_hum_to_mixing_ratio", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        const real_t factor = real_t(1.0) + rtm_zm_exec(i);
        wprtp_exec(i)       = factor * factor * wpqtp_exec(i);
      });
  from_execution_space(wprtp, wprtp_exec);
}

//=============================================================================
template <class RtmView, class ForcingView, class ResultView>
void force_spec_hum_to_mixing_ratio(int ngrdcol, int nzt, const RtmView& rtm,
                                    const ForcingView& qtm_forcing,
                                    const ResultView& rtm_forcing) {
  // Description:
  // Convert d(q_tm)/dt|_f to d(r_tm)/dt|_f on thermodynamic levels using
  // (1 + r_tm)^2 * d(q_tm)/dt|_f, where r_tm is the mean total-water ratio.
  const auto rtm_exec         = to_execution_space(rtm);
  const auto qtm_forcing_exec = to_execution_space(qtm_forcing);
  const auto rtm_forcing_exec = to_execution_space(rtm_forcing, false);
  using Policy                = Kokkos::RangePolicy<typename Real2D::execution_space>;

  // Solve for forcing in terms of total water mixing ratio.
  Kokkos::parallel_for(
      "force_spec_hum_to_mixing_ratio", Policy(0, ngrdcol * nzt), KOKKOS_LAMBDA(int index) {
        const int i            = index % ngrdcol;
        const int k            = index / ngrdcol;
        const real_t factor    = real_t(1.0) + rtm_exec(i, k);
        rtm_forcing_exec(i, k) = factor * factor * qtm_forcing_exec(i, k);
      });
  from_execution_space(rtm_forcing, rtm_forcing_exec);
}

} // namespace clubb_kokkos
