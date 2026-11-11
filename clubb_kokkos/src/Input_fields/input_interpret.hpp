// Interpret sounding coordinates, temperature and subsidence headers.
#pragma once
#include "scratch_workspace.hpp"
#include "input_reader.hpp"
namespace clubb_kokkos {
void read_z_profile(ScratchWorkspace& workspace, int nvar, int nsize,
                    const std::vector<one_dim_read_var>& retVars, real_t p_sfc,
                    real_t zm_init,                                        // In
                    int saturation_formula,                                // In
                    Real1D z, Real1D p_in_Pa, std::string& altitude_type); // Out
void read_theta_profile(ScratchWorkspace& workspace, int nvar, int nsize,
                        const std::vector<one_dim_read_var>& retVars, // In
                        std::string& temperature_type, Real1D theta); // Out
void read_subs_profile(ScratchWorkspace& workspace, int nvar, int nsize,
                       const std::vector<one_dim_read_var>& retVars, // In
                       std::string& subs_type, Real1D subs);         // Out

// Native row/profile actuals cross the existing shared array boundary.
template <class Height, class Pressure>
void read_z_profile(ScratchWorkspace& workspace, int nvar, int nsize,
                    const std::vector<one_dim_read_var>& retVars, real_t p_sfc,
                    real_t zm_init,                                         // In
                    int saturation_formula,                                 // In
                    Height z, Pressure p_in_Pa, std::string& altitude_type) // Out
{
  auto workspace_scope = workspace.scope();
  auto z_exec          = to_execution_space(workspace, z, false);
  auto p_exec          = to_execution_space(workspace, p_in_Pa, false);
  read_z_profile(workspace, nvar, nsize, retVars, p_sfc, zm_init, // In
                 saturation_formula,                              // In
                 z_exec, p_exec, altitude_type);                  // Out
  from_execution_space(workspace, z, z_exec);
  from_execution_space(workspace, p_in_Pa, p_exec);
}
template <class Temperature>
void read_theta_profile(ScratchWorkspace& workspace, int nvar, int nsize,
                        const std::vector<one_dim_read_var>& retVars,     // In
                        std::string& temperature_type, Temperature theta) // Out
{
  auto workspace_scope = workspace.scope();
  auto theta_exec      = to_execution_space(workspace, theta, false);
  read_theta_profile(workspace, nvar, nsize, retVars, temperature_type, theta_exec);
  from_execution_space(workspace, theta, theta_exec);
}
template <class Subsidence>
void read_subs_profile(ScratchWorkspace& workspace, int nvar, int nsize,
                       const std::vector<one_dim_read_var>& retVars, // In
                       std::string& subs_type, Subsidence subs)      // Out
{
  auto workspace_scope = workspace.scope();
  auto subs_exec       = to_execution_space(workspace, subs, false);
  read_subs_profile(workspace, nvar, nsize, retVars, subs_type, subs_exec);
  from_execution_space(workspace, subs, subs_exec);
}
} // namespace clubb_kokkos
