// Store observed profiles and prescribe MPACE-A case forcing.
//
// mpace_a_init reads 38 levels and 139 times of case data before the tendency or
// surface calls. mpace_a_tndcy interpolates six observed profiles in time and
// space to supply temperature/water forcing [K/s,(kg/kg)/s], reference winds [m/s]
// and scalar copies, with zero vertical motion. mpace_a_sfclyr converts supplied
// sensible/latent heat [W/m^2] to heat [K m/s]/water[(kg/kg) m/s] using density [kg/m^3].
// Model fields are column-first; observed tables retain Native(level, time)
// order. The final context replaces threadprivate Native state. Physical table
// reads remain explicit file_read_1d/2d I/O services; all interpolation and
// scientific conversions execute in C++.
// References: ARM MPACE intercomparison; Liou, Wallace and Hobbs, Shettle and Weinman;
// Michael Falk, 2006–2007 and Joshua Fasching, 2007–2008 source notes.
#pragma once

#include "scratch_workspace.hpp"
#include "CLUBB_core/array_index.hpp"
#include "CLUBB_core/grid_class.hpp"
#include <string>

namespace clubb_kokkos {

struct mpace_a_state {
  static constexpr int file_ntimes  = 139;
  static constexpr int file_nlevels = 38;
  static constexpr int per_line     = 5;
  Real1D file_pressure;  // Pressure coordinates retained in the raw file units.
  Real1D file_heights;   // Observed heights [m].
  Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> file_times; // Observed times [s].
  Real2D dTdt_forcing;   // Horizontal temperature tendency [K/hr].
  Real2D dqdt_forcing;   // Horizontal specific-humidity tendency[g/(kg hr)].
  Real2D vertT_forcing;  // Vertical temperature tendency [K/hr].
  Real2D vertq_forcing;  // Vertical specific-humidity tendency[g/(kg hr)].
  Real2D um_obs;         // Observed u wind [m/s].
  Real2D vm_obs;         // Observed v wind [m/s].
  Real1D file_latent_ht; // Surface latent heat [W/m^2].
  Real1D file_sens_ht;   // Surface sensible heat [W/m^2].
};

void mpace_a_tndcy(ScratchWorkspace& workspace, int ngrdcol, int sclr_dim, int edsclr_dim,
                   const sclr_idx_type& sclr_idx, const grid& gr, double time, Real2D p_in_Pa,
                   Real2D wm_zt, Real2D wm_zm, Real2D thlm_forcing, Real2D rtm_forcing,
                   Real2D um_hoc_grid, Real2D vm_hoc_grid, Real3D sclrm_forcing,
                   Real3D edsclrm_forcing, const mpace_a_state& input);
void mpace_a_sfclyr(ScratchWorkspace& workspace, int ngrdcol, double time, Real1D rho_sfc,
                    Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D ustar, const mpace_a_state& input);
void mpace_a_init(int iunit, const std::string& file_path, mpace_a_state& input);

} // namespace clubb_kokkos
