// Case-owned IOP sounding/surface/forcing interpretation.
//
// Native default REAL fields stay float even in a real8 CLUBB build. The state
// replaces readiopdata_module saved fields; NetCDF resources remain call-local.
// Public scientific outputs use the permanent core View owner. Optional surface
// outputs use null pointers for Native ABSENT arguments, in original order.
#pragma once
#include "scratch_workspace.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/grid_class.hpp"
#include <string>
namespace clubb_kokkos {

struct readiopdata_state {
  bool masterproc              = true;
  bool scamiop_from_global_cam = false;
  std::string iopfilepath      = "../input/case_data/BOMEX_5day_4scam.nc";
  bool l_dims                  = false;
  bool l_snd                   = false;
  bool l_sfc                   = false;
  bool l_frc                   = false;
  Kokkos::View<float*, Kokkos::LayoutLeft> dplevs; // Native default-real levels.
  int ntime{};
  int nlev{};
  Kokkos::View<int*, Kokkos::LayoutLeft> tsec; // Native integer times[s].
  bool get_add_surface_data{};
  bool have_omega{};
  Kokkos::View<float**, Kokkos::LayoutLeft> dtls; // Native(level,time).
  Kokkos::View<float**, Kokkos::LayoutLeft> dqls; // Native(level,time).
  // The former core helper's logging flag is case-owned too. Every source
  // entry assigns it from its required masterproc argument before use.
  bool masterproc_{};
};

void readiopdata_dims(ScratchWorkspace& workspace, real_t& lat_val, real_t& lon_val, // Out
                      readiopdata_state& iop);               // Native module InOut
void readiopdata_snd(ScratchWorkspace& workspace,            // InOut
                     real_t zm_init, int saturation_formula, // In
                     int& nlevels, Real1D u, Real1D v, Real1D thlm, Real1D rtm, Real1D p, // Out
                     real_t& psfc, Real1D z, Real1D ug, Real1D vg, Real1D subs,           // Out
                     std::string& subs_type,                                              // Out
                     readiopdata_state& iop); // Native module InOut
void readiopdata_sfc(ScratchWorkspace& workspace, Real1D* tsfc, Real1D* ubar,
                     Real1D* shsfc,                           // Optional Out
                     Real1D* lhsfc, Real1D* timesfc,          // Optional Out
                     readiopdata_state& iop);                 // Native module InOut
void readiopdata_frc(ScratchWorkspace& workspace,             // InOut
                     Real1D rtm, const grid& gr,              // In
                     Real1D thlm_forcing, Real1D rtm_forcing, // Out
                     readiopdata_state& iop);                 // Native module InOut

// Host/strided profile boundaries reuse the authoritative array owner.
template <class Profile>
void readiopdata_snd(ScratchWorkspace& workspace,            // InOut
                     real_t zm_init, int saturation_formula, // In
                     int& nlevels, Profile u, Profile v, Profile thlm, Profile rtm,
                     Profile p,                                                     // Out
                     real_t& psfc, Profile z, Profile ug, Profile vg, Profile subs, // Out
                     std::string& subs_type, readiopdata_state& iop) // Out, Native module InOut
{
  auto workspace_scope = workspace.scope();
  auto u_exec          = to_execution_space(workspace, u, false);
  auto v_exec          = to_execution_space(workspace, v, false);
  auto thlm_exec       = to_execution_space(workspace, thlm, false);
  auto rtm_exec        = to_execution_space(workspace, rtm, false);
  auto p_exec          = to_execution_space(workspace, p, false);
  auto z_exec          = to_execution_space(workspace, z, false);
  auto ug_exec         = to_execution_space(workspace, ug, false);
  auto vg_exec         = to_execution_space(workspace, vg, false);
  auto subs_exec       = to_execution_space(workspace, subs, false);
  readiopdata_snd(workspace, zm_init, saturation_formula, nlevels, u_exec, v_exec, thlm_exec,
                  rtm_exec, p_exec, psfc, z_exec, ug_exec, vg_exec, subs_exec, subs_type, iop);
  from_execution_space(workspace, u, u_exec);
  from_execution_space(workspace, v, v_exec);
  from_execution_space(workspace, thlm, thlm_exec);
  from_execution_space(workspace, rtm, rtm_exec);
  from_execution_space(workspace, p, p_exec);
  from_execution_space(workspace, z, z_exec);
  from_execution_space(workspace, ug, ug_exec);
  from_execution_space(workspace, vg, vg_exec);
  from_execution_space(workspace, subs, subs_exec);
}
template <class Rtm, class Forcing>
void readiopdata_frc(ScratchWorkspace& workspace,               // InOut
                     Rtm rtm, const grid& gr,                   // In
                     Forcing thlm_forcing, Forcing rtm_forcing, // Out
                     readiopdata_state& iop)                    // Native module InOut
{
  auto workspace_scope = workspace.scope();
  auto thlm_exec       = to_execution_space(workspace, thlm_forcing, false);
  auto rtm_exec        = to_execution_space(workspace, rtm_forcing, false);
  readiopdata_frc(workspace, to_execution_space(workspace, rtm), gr, thlm_exec, rtm_exec, iop);
  from_execution_space(workspace, thlm_forcing, thlm_exec);
  from_execution_space(workspace, rtm_forcing, rtm_exec);
}
} // namespace clubb_kokkos
