// Store time-dependent forcing and surface data for one standalone case.
//
// initialize_t_dependent_input reads a case's forcing/surface tables, and
// finalize_t_dependent_input releases them. apply_time_dependent_forcings and
// apply_time_dependent_forcings_from_dycore interpolate and apply the profiles;
// the latter remaps an existing dycore grid onto the physics grid. time_select
// identifies adjacent supplied times as zero-based slots and their interpolation fraction.
//
// Native threadprivate module variables belong to time_dependent_input_state,
// owned by the C++ driver. Public surface series remain named as in Native.
// Time metadata resides on the host for scalar control. Profile records retain
// Native (height,time) order; model fields use column-first LayoutLeft.
// Helpers receive this case state explicitly where Native accesses module storage.
#pragma once

#include "scratch_workspace.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/array_index.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "Input_fields/input_reader.hpp"
#include <vector>

namespace clubb_kokkos {

inline constexpr int nforcings = 10; // Number of columns in the input file.

struct time_dependent_input_state {
  Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>
      time_sfc_given;      // Surface time samples [s].
  Real1D latent_ht_given;  // Latent heat flux [W/m^2].
  Real1D sens_ht_given;    // Sensible heat flux [W/m^2].
  Real1D thlm_sfc_given;   // Surface liquid-water potential temperature [K].
  Real1D rtm_sfc_given;    // Surface total-water mixing ratio [kg/kg].
  Real1D CO2_sfc_given;    // Surface carbon-dioxide input, file units.
  Real1D upwp_sfc_given;   // Surface u-momentum flux [(m/s)^2].
  Real1D vpwp_sfc_given;   // Surface v-momentum flux [(m/s)^2].
  Real1D T_sfc_given;      // Surface temperature [K].
  Real1D wpthlp_sfc_given; // Surface heat flux [K m/s].
  Real1D wpqtp_sfc_given;  // Surface specific-humidity flux [(kg/kg) m/s].

  std::vector<two_dim_read_var> t_dependent_forcing_data =
      std::vector<two_dim_read_var>(nforcings); // Data interpolated to the current grid.
  std::vector<two_dim_read_var> t_dependent_forcing_data_f_grid =
      std::vector<two_dim_read_var>(nforcings); // Raw forcing table for grid adaptation.
  Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>
      forcing_times;                            // Immutable host time samples [s].
  one_dim_read_var dimension_var;               // Other table dimension, normally time [s].

  bool l_t_dependent{};             // Read and use time-dependent information.
  bool l_input_xpwp_sfc{};          // Read surface momentum fluxes.
  bool l_ignore_forcings{};         // Ignore the forcing table for this case.
  bool l_sfc_already_initialized{}; // Surface input was initialized previously.
};

struct readiopdata_state;

void initialize_t_dependent_input(ScratchWorkspace& workspace, // InOut
                                  int iunit, const std::string& runtype, bool l_readiopdata,
                                  Real2D rho,                                               // In
                                  int nz, Real1D grid, Real1D p_in_Pa,                      // In
                                  bool l_add_dycore_grid, int grid_adapt_in_time_method,    // In
                                  time_dependent_input_state& input, readiopdata_state& iop // InOut
);
void finalize_t_dependent_input(time_dependent_input_state& input);

void apply_time_dependent_forcings(ScratchWorkspace& workspace, int ngrdcol, int nzm, int nzt,
                                   int sclr_dim, int edsclr_dim, const sclr_idx_type& sclr_idx,
                                   bool l_do_expldiff_rtm_thlm, const grid& gr, double time,
                                   Real2D rtm, Real2D rho, Real2D exner, Real2D thlm_f,
                                   Real2D rtm_f, Real2D um_ref, Real2D vm_ref, Real2D um_f,
                                   Real2D vm_f, Real2D wm_zt, Real2D wm_zm, Real2D ug, Real2D vg,
                                   Real3D sclrm_forcing, Real3D edsclrm_forcing,
                                   const time_dependent_input_state& input);
void apply_time_dependent_forcings_from_dycore(
    ScratchWorkspace& workspace, int ngrdcol, int nzm, int nzt, int sclr_dim, int edsclr_dim,
    const sclr_idx_type& sclr_idx, bool l_do_expldiff_rtm_thlm, const grid& gr,
    const grid& gr_dycore, double time, Real2D rtm, Real2D rho, Real2D exner, int grid_remap_method,
    int total_idx_rho_lin_spline, Real2D rho_lin_spline_vals, Real2D rho_lin_spline_levels,
    Real1D p_sfc, Real2D thlm_f, Real2D rtm_f, Real2D um_ref, Real2D vm_ref, Real2D um_f,
    Real2D vm_f, Real2D wm_zt, Real2D wm_zm, Real2D ug, Real2D vg, Real3D sclrm_forcing,
    Real3D edsclrm_forcing, const time_dependent_input_state& input);

void time_select(
    double time, int nvar,                                                         // In
    Kokkos::View<const real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> time_array, // In
    int& before_time, int& after_time, real_t& time_frac                           // Out
);

// Host or strided driver fields cross this storage boundary once per call.
template <class Field_rho, class Field_grid, class Field_p_in_Pa>
void initialize_t_dependent_input(ScratchWorkspace& workspace, // InOut
                                  int iunit, const std::string& runtype, bool l_readiopdata,
                                  Field_rho rho,                                            // In
                                  int nz, Field_grid grid, Field_p_in_Pa p_in_Pa,           // In
                                  bool l_add_dycore_grid, int grid_adapt_in_time_method,    // In
                                  time_dependent_input_state& input, readiopdata_state& iop // InOut
) {
  // Description:
  // Stage the Host/strided input density, grid and pressure for initialization.
  // The case contexts own the resulting forcing/surface and IOP data.
  auto workspace_scope   = workspace.scope();
  auto rho_execution     = to_execution_space(workspace, rho);
  auto grid_execution    = to_execution_space(workspace, grid);
  auto p_in_Pa_execution = to_execution_space(workspace, p_in_Pa);
  initialize_t_dependent_input(workspace, iunit, runtype, l_readiopdata, rho_execution, nz,
                               grid_execution, p_in_Pa_execution, l_add_dycore_grid,
                               grid_adapt_in_time_method, input, iop);
}

} // namespace clubb_kokkos
