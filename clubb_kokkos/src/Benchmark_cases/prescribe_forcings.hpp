// Coordinate case-specific large-scale tendencies and surface boundary fluxes.
//
// prescribe_forcings resets timestep forcing, applies time-dependent or legacy
// case tendencies, derives the surface input, dispatches case fluxes and updates
// the C++ statistics context. The Native argument order is kept;
// appended time-dependent and MPACE contexts replace implicit module state.
// Model fields use (column, level[,scalar]) LayoutLeft Views. Heights [m], times [s],
// thermodynamics [K,kg/kg,Pa,kg/m^3] and winds [m/s] determine tendencies and fluxes
// in the corresponding field units per second or multiplied by m/s. Passive-scalar
// slots are zero-based; -1 marks an absent scalar. Temporary Views use the caller's
// workspace scopes; retained storage belongs to the case.
// The private boundary helper preserves the default lowest-level sampling and
// the optional fixed 25 m convergence-test interpolation in Native call order.
#pragma once

#include "scratch_workspace.hpp"
#include "CLUBB_core/array_index.hpp"
#include "CLUBB_core/err_info_type_module.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "mpace_a.hpp"
#include "time_dependent_input.hpp"
#include "CLUBB_core/stats_netcdf.hpp"
#include <string>

namespace clubb_kokkos {

void prescribe_forcings(ScratchWorkspace& workspace,                                 // InOut
                        const grid& gr, int nzm, int nzt, int ngrdcol,               // In
                        int sclr_dim, int edsclr_dim, const sclr_idx_type& sclr_idx, // In
                        bool l_do_expldiff_rtm_thlm,                                 // In
                        const std::string& runtype, int sfctype,                     // In
                        double time_current, double time_initial, real_t dt,         // In
                        bool l_readiopdata,                                          // In
                        Real2D um, Real2D vm, Real2D thlm,                           // In
                        Real2D p_in_Pa, Real2D exner, Real2D rho, Real2D rho_zm, Real2D thvm, // In
                        Real1D veg_T_in_K,                                                    // In
                        bool l_modify_bc_for_cnvg_test, int saturation_formula,               // In
                        stats_type& stats,                                                 // InOut
                        bool l_add_dycore_grid, int grid_remap_method,                     // In
                        int total_idx_rho_lin_spline, Real2D rho_lin_spline_vals,          // In
                        Real2D rho_lin_spline_levels, const grid& gr_dycore,               // In
                        Real2D rtm, Real2D wm_zm, Real2D wm_zt, Real2D ug, Real2D vg,      // InOut
                        Real2D um_ref, Real2D vm_ref,                                      // InOut
                        Real2D thlm_forcing, Real2D rtm_forcing, Real2D um_forcing,        // InOut
                        Real2D vm_forcing, Real2D wprtp_forcing, Real2D wpthlp_forcing,    // InOut
                        Real2D rtp2_forcing, Real2D thlp2_forcing, Real2D rtpthlp_forcing, // InOut
                        Real3D wpsclrp, Real3D sclrm_forcing, Real3D edsclrm_forcing,      // InOut
                        Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D upwp_sfc,
                        Real1D vpwp_sfc,                                                    // InOut
                        Real1D T_sfc, Real1D p_sfc, real_t& sens_ht, real_t& latent_ht,     // InOut
                        Real2D wpsclrp_sfc, Real2D wpedsclrp_sfc, err_info_type& err_info,  // InOut
                        const time_dependent_input_state& input, const mpace_a_state& mpace // In
);

} // namespace clubb_kokkos
