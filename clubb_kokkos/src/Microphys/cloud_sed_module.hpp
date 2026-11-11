// Compute cloud-droplet sedimentation and add its water/temperature tendencies.
//
// cloud_drop_sed stages host/strided fields into column-first LayoutLeft Views
// on the Kokkos execution space. It interpolates cloud water/number to momentum
// levels, differentiates the downward droplet flux and updates rcm_mc/thlm_mc
// in place. Statistics accumulate through the case-owned Native backend;
// all scratch is call-local.
// Reference: Ackerman et al. (2009), MWR, doi:10.1175/2008MWR2582.1.
#pragma once
#include "scratch_workspace.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "CLUBB_core/stats_netcdf.hpp"
namespace clubb_kokkos {
void cloud_drop_sed(ScratchWorkspace& workspace,                             // InOut
                    const grid& gr, int ngrdcol, Real2D rcm, Real2D Ncm,     // In
                    Real2D rho_zm, Real2D rho, Real2D exner, real_t sigma_g, // In
                    stats_type& stats, Real2D rcm_mc, Real2D thlm_mc);       // InOut

template <class Field_rcm, class Field_Ncm, class Field_rho_zm, class Field_rho, class Field_exner,
          class Field_rcm_mc, class Field_thlm_mc>
void cloud_drop_sed(ScratchWorkspace& workspace,                                           // InOut
                    const grid& gr, int ngrdcol, Field_rcm rcm, Field_Ncm Ncm,             // In
                    Field_rho_zm rho_zm, Field_rho rho, Field_exner exner, real_t sigma_g, // In
                    stats_type& stats, Field_rcm_mc rcm_mc, Field_thlm_mc thlm_mc) {       // InOut
  // Stage host/strided callers through the shared array owner. Both tendencies
  // are Native InOut fields: preserve their input values and caller allocations.
  auto workspace_scope   = workspace.scope();
  auto rcm_mc_execution  = to_execution_space(workspace, rcm_mc);
  auto thlm_mc_execution = to_execution_space(workspace, thlm_mc);
  cloud_drop_sed(workspace, gr, ngrdcol, to_execution_space(workspace, rcm),
                 to_execution_space(workspace, Ncm),                                        // In
                 to_execution_space(workspace, rho_zm), to_execution_space(workspace, rho), // In
                 to_execution_space(workspace, exner), sigma_g,                             // In
                 stats, rcm_mc_execution, thlm_mc_execution);                               // InOut
  from_execution_space(workspace, rcm_mc, rcm_mc_execution);
  from_execution_space(workspace, thlm_mc, thlm_mc_execution);
}
} // namespace clubb_kokkos
