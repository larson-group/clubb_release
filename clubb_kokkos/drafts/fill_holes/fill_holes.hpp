// Fill negative or sub-threshold fields while conserving the available mass or TKE.
//
// fill_holes_vertical_api dispatches global/sliding-window vertical filling;
// fill_holes_wp2_from_horz_tke transfers horizontal TKE at the same level.
// hole_filling_hm_one_lev and fill_holes_hydromet_api redistribute frozen species;
// fill_holes_wv transfers residual holes to vapor/latent-heat tendencies.
// fill_holes_driver_api coordinates those stages and Native budget calls.
// clip_hydromet_conc_mvr enforces particle-radius limits; setup_stats_names
// supplies species budgets/fall speeds. Arrays use LayoutLeft; batched physical
// fields are column first, while single-column APIs retain Native level/species
// axes. Grid/metadata indices remain Native1based. Other vertical methods are
// outside actual JAX-supported scope. Draft only; see INTEGRATION.md.
#pragma once
#include "array_types.hpp"
#include "grid_class.hpp"
#include "corr_varnce_module.hpp"
#include "f2k/stats_bridge.hpp"
#include "constants_clubb.hpp"
#include <string>
#include <limits>
#include <Kokkos_MathematicalFunctions.hpp>
namespace clubb_kokkos {
// nz: Number of grid levels [-]
// ngrdcol: Number of columns [-]
// threshold: A threshold (e.g. w_tol*w_tol) below which field must not fall
// [Units vary; same as field]
// lower_hf_level: Lower grid level of global hole-filling range      []
// upper_hf_level: Upper grid level of global hole-filling range      []
// dz: Spacing between thermodynamic grid levels; centered over momentum grid levels OR Spcaing
// between momentum grid levels; centered over thermodynamic grid levels
// rho_ds: Dry, static density on thermodynamic or momentum levels    [kg/m^3]
// grid_dir_indx: Grid direction index (+1 ascending; -1 descending)
// fill_holes_type: Option for which type of hole filler to use
// field: The field (e.g. wp2) that contains holes [Units same as threshold]
void fill_holes_vertical_api(int nz, int ngrdcol, real_t threshold,       // In
                             int lower_hf_level, int upper_hf_level,      // In
                             Real2D dz, Real2D rho_ds, int grid_dir_indx, // In
                             int fill_holes_type,                         // In
                             Real2D field)                                // InOut
    ;

// nz: Number of grid levels [-]
// ngrdcol: Number of columns [-]
// threshold: A threshold (e.g. w_tol*w_tol) below which field must not fall
// [Units vary; same as field]
// lower_hf_level: Lower grid level of global hole-filling range      []
// upper_hf_level: Upper grid level of global hole-filling range      []
// wp2: The field (here wp2) that contains holes [m^2/s^2]
// up2: Zonal TKE component supplying the hole [m^2/s^2]
// vp2: The fields (here up2, vp2) that will be used to fill holes [m^2/s^2]
void fill_holes_wp2_from_horz_tke(int nz, int ngrdcol, real_t threshold,  // In
                                  int lower_hf_level, int upper_hf_level, // In
                                  Real2D wp2, Real2D up2, Real2D vp2)     // InOut
    ;

// num_hm_fill: number of hydrometeors involved
// hm_one_lev: Same-phase species mixing ratios at one level [kg/kg]
// hm_one_lev_filled: Same-phase mixing ratios after redistribution [kg/kg]
void hole_filling_hm_one_lev(int num_hm_fill, Real1D hm_one_lev, // In
                             Real1D hm_one_lev_filled)           // Out
    ;

// nzt: Number of thermodynamic levels [-]
// hydromet_dim: Number of hydrometeor species [-]
// hydromet: Hydrometeor fields (level,species) [units vary]
// l_frozen_hm: if true, then the hydrometeor is frozen; otherwise liquid
// l_mix_rat_hm: if true, then the quantity is a hydrometeor mixing ratio
// hydromet_filled: Hydrometeor fields after frozen-species redistribution [units vary]
void fill_holes_hydromet_api(int nzt, int hydromet_dim, Real2D hydromet, // In
                             Bool1D l_frozen_hm, Bool1D l_mix_rat_hm,    // In
                             Real2D hydromet_filled)                     // Out
    ;

// nzt: Number of thermodynamic levels [-]
// dt: Timestep         [s]
// exner: Exner function                            [-]
// hydromet_name: Fixed-width species name [-]
// rvm_mc: Water-vapor tendency [kg/kg/s]
// thlm_mc: Liquid-water potential-temperature tendency [K/s]
// hydromet: Hydrometeor array                         [units vary]
void fill_holes_wv(int nzt, real_t dt, Real1D exner, const std::string &hydromet_name, // In
                   Real1D rvm_mc, Real1D thlm_mc, Real1D hydromet)                     // InOut
    ;

// gr: Native grid coordinates and spacings [-,m]
// ngrdcol: Number of model columns
// nzt: Number of vertical levels
// dt: Timestep         [s]
// hydromet_dim: Number of hydrometeor fields
// hm_metadata: Hydrometeor species indices, flags and tolerances [-,units vary]
// l_fill_holes_hm: Enable same-level frozen-species redistribution [-]
// rho_ds_zt: Dry, static density on thermo. levels    [kg/m^3]
// exner: Exner function                           [-]
// stats: Native statistics context, mutated by sampled budget calls [-]
// thlm_mc: Microphysics contributions to liquid potential temp [K/s]
// rvm_mc: Microphysics contributions to vapor water           [kg/kg/s]
// hydromet: Mean of hydrometeor fields    [units vary]
void fill_holes_driver_api(const grid &gr, int ngrdcol, int nzt, real_t dt, int hydromet_dim, // In
                           const hm_metadata_type &hm_metadata, bool l_fill_holes_hm,         // In
                           Real2D rho_ds_zt, Real2D exner, int fill_holes_type,               // In
                           f2k::stats_type &stats, Real2D thlm_mc, Real2D rvm_mc,
                           Real3D hydromet) // InOut
    ;

// nzt: Number of thermodynamic levels [-]
// hydromet_dim: Number of hydrometeor fields
// hm_metadata: Hydrometeor species indices, flags and tolerances [-,units vary]
// hydromet: Mean of hydrometeor fields    [units vary]
// hydromet_clipped: Clipped mean of hydrometeor fields    [units vary]
void clip_hydromet_conc_mvr(int nzt, int hydromet_dim, const hm_metadata_type &hm_metadata, // In
                            Real2D hydromet,                                                // In
                            Real2D hydromet_clipped)                                        // Out
    ;

// ihm: Native1based species index [-]
// hydromet_dim: Number of hydrometeor species [-]
// hydromet_list: Fixed-width species names (character,species) [-]
// max_velocity: Maximum sedimentation velocity [m/s]
void setup_stats_names(int ihm, int hydromet_dim,                               // In
                       Kokkos::View<char **, Kokkos::LayoutLeft> hydromet_list, // In
                       real_t &max_velocity,                                    // InOut
                       std::string *name_bt = nullptr,
                       std::string *name_hf = nullptr,   // Optional Out
                       std::string *name_wvhf = nullptr, // Optional Out
                       std::string *name_cl = nullptr,
                       std::string *name_mc = nullptr) // Optional Out
    ;

// Device form of the same Native one-level API, used by the batched frozen fill.
// It preserves the ordered species sums; host wrappers own warning I/O.
template <class Input, class Output>
KOKKOS_INLINE_FUNCTION void hole_filling_hm_one_lev(int num_hm_fill, Input hm_one_lev, // In
                                                    Output hm_one_lev_filled)          // Out
{
  real_t total_hole = real_t(0.0), total_mass = real_t(0.0);
  int num_neg_hm = 0;
  // Initialize and determine missing mass and available same-phase material.
  for (int i = 1; i <= num_hm_fill; ++i)
    hm_one_lev_filled(i - 1) = real_t(0.0);
  for (int i = 1; i <= num_hm_fill; ++i) {
    if (hm_one_lev(i - 1) < zero) {
      total_hole = total_hole + hm_one_lev(i - 1);
      num_neg_hm = num_neg_hm + 1;
    } else {
      total_mass = total_mass + hm_one_lev(i - 1);
    }
  }
  // No material: leave the field unchanged.
  if (Kokkos::abs(total_mass) < eps) {
    for (int i = 1; i <= num_hm_fill; ++i)
      hm_one_lev_filled(i - 1) = hm_one_lev(i - 1);
    return;
  }
  // Insufficient mass fills holes partly; otherwise remove them and scale donors.
  for (int i = 1; i <= num_hm_fill; ++i) {
    if (Kokkos::abs(total_hole) > total_mass)
      hm_one_lev_filled(i - 1) =
          (hm_one_lev(i - 1) < zero ? hm_one_lev(i - 1) : zero) * (one + total_mass / total_hole);
    else
      hm_one_lev_filled(i - 1) =
          (hm_one_lev(i - 1) > zero ? hm_one_lev(i - 1) : zero) * (one + total_hole / total_mass);
  }
}

// Moment callers pass strided scalar/column sections. Adapt once at the public
// boundary using the existing array owner; the Native axis meanings do not change.
template <class Dz, class Rho, class Field>
void fill_holes_vertical_api(int nz, int ngrdcol, real_t threshold, int lower_hf_level,
                             int upper_hf_level, Dz dz, Rho rho_ds, int grid_dir_indx,
                             int fill_holes_type, Field field) {
  auto dz_exec = to_execution_space(dz);
  auto rho_exec = to_execution_space(rho_ds);
  auto field_exec = to_execution_space(field);
  fill_holes_vertical_api(nz, ngrdcol, threshold, lower_hf_level, upper_hf_level, dz_exec, rho_exec,
                          grid_dir_indx, fill_holes_type, field_exec);
  from_execution_space(field, field_exec);
}
template <class W, class U, class V>
void fill_holes_wp2_from_horz_tke(int nz, int ngrdcol, real_t threshold, int lower_hf_level,
                                  int upper_hf_level, W wp2, U up2, V vp2) {
  auto wp2_exec = to_execution_space(wp2);
  auto up2_exec = to_execution_space(up2);
  auto vp2_exec = to_execution_space(vp2);
  fill_holes_wp2_from_horz_tke(nz, ngrdcol, threshold, lower_hf_level, upper_hf_level, wp2_exec,
                               up2_exec, vp2_exec);
  from_execution_space(wp2, wp2_exec);
  from_execution_space(up2, up2_exec);
  from_execution_space(vp2, vp2_exec);
}

} // namespace clubb_kokkos
