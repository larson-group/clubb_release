// Initialize model profiles from Native-compatible sounding tables.
//
// All original arguments retain their order. The three appended contexts own
// former Native module state for radiation/extended profiles and IOP input.
// Fields use (column, level[, scalar]) Views. Temperature follows temperature_type
// [K] until driver conversion; water is [kg/kg], pressure [Pa] and winds [m/s].
// Subsidence follows subs_type [m/s or Pa/s]; surface fields use [K] and [kg/kg].
#pragma once
#include "scratch_workspace.hpp"
#include "CLUBB_core/array_index.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "Radiation/extended_atmosphere_module.hpp"
#include "readiopdata_module.hpp"
#include <string>
namespace clubb_kokkos {
inline constexpr int sclr_max = 1000; // Native public maximum scalar count.
void read_sounding(ScratchWorkspace& workspace, int ngrdcol, int sclr_dim, int edsclr_dim,
                   const sclr_idx_type& sclr_idx,                                         // In
                   const grid& gr, int iunit, const std::string& runtype, real_t zm_init, // In
                   int saturation_formula,                                                // In
                   bool l_readiopdata, bool l_modify_ic_with_cubic_int,                   // In
                   real_t& p_sfc,                                                         // InOut
                   Real2D thlm, std::string& temperature_type, Real2D rtm, Real2D um, Real2D vm,
                   Real2D ugm, // Out
                   Real2D vgm, // Out
                   std::string& altitude_type, Real2D press, std::string& subs_type,
                   Real2D wm,                                                     // Out
                   Real1D rtm_sfc, Real1D thlm_sfc, Real3D sclrm, Real3D edsclrm, // Out
                   const parameters_radiation_state& parameters_radiation, // Native module In
                   extended_atmosphere_state& extended_atmosphere,         // Native module InOut
                   readiopdata_state& iop);                                // Native module InOut
void read_profile(ScratchWorkspace& workspace, const grid& gr, const std::string& fname, // In
                  Real1D x);                                                             // Out

template <class Thlm, class Rtm, class Um, class Vm, class Ugm, class Vgm, class Pressure, class Wm,
          class RtmSfc, class ThlmSfc, class Sclr, class Edsclr>
void read_sounding(ScratchWorkspace& workspace, int ngrdcol, int sclr_dim, int edsclr_dim,
                   const sclr_idx_type& sclr_idx,                                         // In
                   const grid& gr, int iunit, const std::string& runtype, real_t zm_init, // In
                   int saturation_formula,                                                // In
                   bool l_readiopdata, bool l_modify_ic_with_cubic_int,                   // In
                   real_t& p_sfc,                                                         // InOut
                   Thlm thlm, std::string& temperature_type, Rtm rtm, Um um, Vm vm, Ugm ugm,
                   Vgm vgm,                                                                   // Out
                   std::string& altitude_type, Pressure press, std::string& subs_type, Wm wm, // Out
                   RtmSfc rtm_sfc, ThlmSfc thlm_sfc, Sclr sclrm, Edsclr edsclrm,              // Out
                   const parameters_radiation_state& parameters_radiation, // Native module In
                   extended_atmosphere_state& extended_atmosphere,         // Native module InOut
                   readiopdata_state& iop)                                 // Native module InOut
{
  auto workspace_scope = workspace.scope();
  auto thlm_exec       = to_execution_space(workspace, thlm, false);
  auto rtm_exec        = to_execution_space(workspace, rtm, false);
  auto um_exec         = to_execution_space(workspace, um, false);
  auto vm_exec         = to_execution_space(workspace, vm, false);
  auto ugm_exec        = to_execution_space(workspace, ugm, false);
  auto vgm_exec        = to_execution_space(workspace, vgm, false);
  auto press_exec      = to_execution_space(workspace, press, false);
  auto wm_exec         = to_execution_space(workspace, wm, false);
  auto rtm_sfc_exec    = to_execution_space(workspace, rtm_sfc, false);
  auto thlm_sfc_exec   = to_execution_space(workspace, thlm_sfc, false);
  auto sclrm_exec      = to_execution_space(workspace, sclrm, false);
  auto edsclrm_exec    = to_execution_space(workspace, edsclrm, false);
  read_sounding(workspace, ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr, iunit, runtype, zm_init,
                saturation_formula, l_readiopdata, l_modify_ic_with_cubic_int, p_sfc, thlm_exec,
                temperature_type, rtm_exec, um_exec, vm_exec, ugm_exec, vgm_exec, altitude_type,
                press_exec, subs_type, wm_exec, rtm_sfc_exec, thlm_sfc_exec, sclrm_exec,
                edsclrm_exec, parameters_radiation, extended_atmosphere, iop);
  from_execution_space(workspace, thlm, thlm_exec);
  from_execution_space(workspace, rtm, rtm_exec);
  from_execution_space(workspace, um, um_exec);
  from_execution_space(workspace, vm, vm_exec);
  from_execution_space(workspace, ugm, ugm_exec);
  from_execution_space(workspace, vgm, vgm_exec);
  from_execution_space(workspace, press, press_exec);
  from_execution_space(workspace, wm, wm_exec);
  from_execution_space(workspace, rtm_sfc, rtm_sfc_exec);
  from_execution_space(workspace, thlm_sfc, thlm_sfc_exec);
  from_execution_space(workspace, sclrm, sclrm_exec);
  from_execution_space(workspace, edsclrm, edsclrm_exec);
}
template <class Profile>
void read_profile(ScratchWorkspace& workspace, const grid& gr, const std::string& fname,
                  Profile x) {
  auto workspace_scope = workspace.scope();
  auto x_exec          = to_execution_space(workspace, x, false);
  read_profile(workspace, gr, fname, x_exec);
  from_execution_space(workspace, x, x_exec);
}
} // namespace clubb_kokkos
