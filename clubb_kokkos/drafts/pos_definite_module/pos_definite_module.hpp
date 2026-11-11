// Limit scalar fluxes while retaining their conservative field update.
//
// pos_definite_adj changes field_np1 and flux_np1 in place and fills field_pd
// and flux_pd budgets. Field arrays are (column,thermodynamic level), flux
// arrays are (column,momentum level), all LayoutLeft. The physical-direction
// limiter loop and storage-endpoint fluxes follow Native grid metadata. Two
// separate column ANY checks select budgets before/after the field update.
// All scratch is owned per call. See INTEGRATION.md for pending validation.
#pragma once
#include "array_types.hpp"
#include "grid_class.hpp"

namespace clubb_kokkos {

// dt: Timestep [s].
// field_np1: Field at n+1 (e.g. rtm [kg/kg]); updated in place.
// flux_np1: Flux applied to the field; replaced with limited flux.
// field_n: Field at n, prior to n+1, with the same units as field_np1.
// field_pd: Budget of the change in the mean field due to the scheme [field units/s].
// flux_pd: Budget of the flux change due to the scheme [flux units/s].
void pos_definite_adj(int nzm, int nzt, int ngrdcol, const grid& gr, real_t dt, // In
                      Real2D field_np1, Real2D flux_np1,                        // InOut
                      Real2D field_n,                                           // In
                      Real2D field_pd, Real2D flux_pd);                         // Out

// Preserve caller-owned sections through the existing shared array owner.
template <class FieldNp1, class FluxNp1, class FieldN, class FieldPd, class FluxPd>
void pos_definite_adj(int nzm, int nzt, int ngrdcol, const grid& gr, real_t dt, // In
                      FieldNp1 field_np1, FluxNp1 flux_np1,                     // InOut
                      FieldN field_n,                                           // In
                      FieldPd field_pd, FluxPd flux_pd)                         // Out
{
  // Stage only this public array boundary; concrete calls retain core storage.
  auto field_np1_exec = to_execution_space(field_np1);
  auto flux_np1_exec = to_execution_space(flux_np1);
  auto field_n_exec = to_execution_space(field_n);
  auto field_pd_exec = to_execution_space(field_pd, false);
  auto flux_pd_exec = to_execution_space(flux_pd, false);
  pos_definite_adj(nzm, nzt, ngrdcol, gr, dt, field_np1_exec, flux_np1_exec, field_n_exec,
                   field_pd_exec, flux_pd_exec);
  from_execution_space(field_np1, field_np1_exec);
  from_execution_space(flux_np1, flux_np1_exec);
  from_execution_space(field_pd, field_pd_exec);
  from_execution_space(flux_pd, flux_pd_exec);
}

} // namespace clubb_kokkos
