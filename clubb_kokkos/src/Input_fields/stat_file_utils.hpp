// Average and interpolate input NetCDF profiles onto CLUBB vertical grids.
//
// stat_file_average/interval read inclusive Native record ranges and return
// execution-space profiles. The count/vertical-level queries explicitly open and
// close a descriptor; l_netcdf_file checks the .nc suffix. The two grid helpers
// return inclusive zero-based bounds/brackets, with -1 unset and a 1e-6 m exact test.
// Out-of-domain averages are +infinity except an optional lower ghost value.
// Case workspace owns transient buffers; file descriptors own their altitude arrays.
// Owned returns remain independent; direct-output forms fill caller storage.
#pragma once
#include "input_netcdf.hpp"
#include <vector>
namespace clubb_kokkos {
// Execution-space profile mappings accept an existing strided column directly.
using InputProfile      = Kokkos::View<real_t*, Kokkos::LayoutStride>;
using ConstInputProfile = Kokkos::View<const real_t*, Kokkos::LayoutStride>;
using InputIndices      = Kokkos::View<int*, Kokkos::LayoutLeft>;

Real1D stat_file_average(ScratchWorkspace& workspace,             // InOut
                         const std::string& filename, int out_nz, // In
                         int t1, int t2, ConstInputProfile out_heights,
                         const std::string& variable_name,      // In
                         int npower, bool l_spec_bound_cond,    // In
                         bool& l_error,                         // Out
                         const inputfields_state& inputfields); // Native module In
Real1D stat_file_average_interval(ScratchWorkspace& workspace, const std::string& filename, int nz,
                                  const std::vector<int>& t,                 // In
                                  const std::string& variable_name,          // In
                                  ConstInputProfile out_heights, int npower, // In
                                  bool& l_error,                             // Out
                                  const inputfields_state& inputfields);     // Native module In
int stat_file_num_vertical_levels(ScratchWorkspace& workspace, const std::string& varname,
                                  const std::string& filename,           // In
                                  const inputfields_state& inputfields); // Native module In
Real1D stat_file_vertical_levels(ScratchWorkspace& workspace, const std::string& varname,
                                 const std::string& filename, int nz,   // In
                                 const inputfields_state& inputfields); // Native module In
bool l_netcdf_file(const std::string& filename);                        // In
void CLUBB_levels_within_LES_domain(ScratchWorkspace& workspace, const stat_file& fread_var,
                                    ConstInputProfile CLUBB_grid,               // In
                                    int& k_lowest_input, int& k_highest_input); // Out
void LES_grid_to_CLUBB_grid(ScratchWorkspace& workspace, const stat_file& fread_var,
                            ConstInputProfile CLUBB_grid, int k,    // In
                            int& exact_lev_idx, int& lower_lev_idx, // Out
                            int& upper_lev_idx, bool& l_lin_int);   // Out

// Direct-output overloads retain the Native required inputs before the result.
void stat_file_average(ScratchWorkspace& workspace, const std::string& filename, int out_nz, int t1,
                       int t2, ConstInputProfile out_heights, const std::string& variable_name,
                       int npower, bool l_spec_bound_cond, bool& l_error,
                       const inputfields_state& inputfields, InputProfile result);
void stat_file_average_interval(ScratchWorkspace& workspace, const std::string& filename, int nz,
                                const std::vector<int>& t, const std::string& variable_name,
                                ConstInputProfile out_heights, int npower, bool& l_error,
                                const inputfields_state& inputfields, InputProfile result);
void stat_file_vertical_levels(ScratchWorkspace& workspace, const std::string& varname,
                               const std::string& filename, int nz,
                               const inputfields_state& inputfields, InputProfile result);
void LES_grid_to_CLUBB_grid(const stat_file& fread_var, ConstInputProfile CLUBB_grid,
                            int first_level, int last_level, InputIndices exact_lev_idx,
                            InputIndices lower_lev_idx, InputIndices upper_lev_idx,
                            Bool1D l_lin_int);

template <class Heights>
Real1D stat_file_average(ScratchWorkspace& workspace,             // InOut
                         const std::string& filename, int out_nz, // In
                         int t1, int t2, Heights out_heights,
                         const std::string& variable_name,     // In
                         int npower, bool l_spec_bound_cond,   // In
                         bool& l_error,                        // Out
                         const inputfields_state& inputfields) // Native module In
{
  // Stage host heights for the execution-space reader; preserve the owned result.
  auto scope         = workspace.scope();
  const auto heights = [&]() {
    if constexpr (Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                             typename Heights::memory_space>::accessible)
      return ConstInputProfile(out_heights);
    else
      return ConstInputProfile(to_execution_space(workspace, out_heights));
  }();
  return stat_file_average(workspace, filename, out_nz, t1, t2, heights, variable_name, npower,
                           l_spec_bound_cond, l_error, inputfields);
}
template <class Heights>
Real1D stat_file_average_interval(ScratchWorkspace& workspace, const std::string& filename, int nz,
                                  const std::vector<int>& t,            // In
                                  const std::string& variable_name,     // In
                                  Heights out_heights, int npower,      // In
                                  bool& l_error,                        // Out
                                  const inputfields_state& inputfields) // Native module In
{
  // Stage host heights for the execution-space reader; preserve the owned result.
  auto scope         = workspace.scope();
  const auto heights = [&]() {
    if constexpr (Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                             typename Heights::memory_space>::accessible)
      return ConstInputProfile(out_heights);
    else
      return ConstInputProfile(to_execution_space(workspace, out_heights));
  }();
  return stat_file_average_interval(workspace, filename, nz, t, variable_name, heights, npower,
                                    l_error, inputfields);
}
template <class Heights>
void CLUBB_levels_within_LES_domain(ScratchWorkspace& workspace, const stat_file& fread_var,
                                    Heights CLUBB_grid,                        // In
                                    int& k_lowest_input, int& k_highest_input) // Out
{
  // Stage a host grid for the Native ordered bound search.
  auto scope         = workspace.scope();
  const auto heights = [&]() {
    if constexpr (Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                             typename Heights::memory_space>::accessible)
      return ConstInputProfile(CLUBB_grid);
    else
      return ConstInputProfile(to_execution_space(workspace, CLUBB_grid));
  }();
  CLUBB_levels_within_LES_domain(workspace, fread_var, heights, k_lowest_input, k_highest_input);
}
template <class Heights>
void LES_grid_to_CLUBB_grid(ScratchWorkspace& workspace, const stat_file& fread_var,
                            Heights CLUBB_grid, int k,              // In
                            int& exact_lev_idx, int& lower_lev_idx, // Out
                            int& upper_lev_idx, bool& l_lin_int)    // Out
{
  // The scalar boundary remains available; model loops use the batched overload.
  auto scope         = workspace.scope();
  const auto heights = [&]() {
    if constexpr (Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                             typename Heights::memory_space>::accessible)
      return ConstInputProfile(CLUBB_grid);
    else
      return ConstInputProfile(to_execution_space(workspace, CLUBB_grid));
  }();
  LES_grid_to_CLUBB_grid(workspace, fread_var, heights, k, exact_lev_idx, lower_lev_idx,
                         upper_lev_idx, l_lin_int);
}
} // namespace clubb_kokkos
