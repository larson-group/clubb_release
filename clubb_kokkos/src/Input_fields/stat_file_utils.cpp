// Average and interpolate input NetCDF fields while preserving Native order.
#include "stat_file_utils.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/interpolation.hpp"
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace clubb_kokkos {

void stat_file_average(ScratchWorkspace& workspace,             // InOut
                       const std::string& filename, int out_nz, // In
                       int t1, int t2, ConstInputProfile out_heights,
                       const std::string& variable_name,                  // In
                       int npower, bool l_spec_bound_cond, bool& l_error, // In, In, Out
                       const inputfields_state& inputfields,
                       InputProfile result) { // Module In, Out
  // Average a NetCDF variable over inclusive Native records t1:t2. Out-of-domain
  // values are positive infinity, except the optional lowest ghost value.
  // Interpolation, boundary handling and accumulation are fused per level;
  // records retain their increasing source order and exact scalar formulas.
  // The result must not alias the height grid, which remains fixed across records.
  auto scope = workspace.scope();
  stat_file faverage;
  int k_lowest_input, k_highest_input;
  l_error                 = false;
  const int num_timesteps = (t2 - t1) + 1;
  Kokkos::deep_copy(result, real_t(0.0));
  open_netcdf_read(workspace, variable_name, filename, faverage, l_error, inputfields);
  if (l_error)
    return;
  const int file_nz = faverage.z.extent_int(0);
  auto file_variable =
      scope.view<Real1D>(Kokkos::WithoutInitializing, "stat_file_average::file", file_nz);
  bool l_interpolate;
  if (out_nz != file_nz) {
    l_interpolate = true;
  } else {
    const auto file_z = faverage.z;
    int differing     = 0;
    Kokkos::parallel_reduce(
        "stat_file_average::same_grid", out_nz,
        KOKKOS_LAMBDA(int k, int& count) {
          if (Kokkos::abs(file_z(k) - out_heights(k)) >
              Kokkos::abs(file_z(k) + out_heights(k)) * eps / 2)
            count += 1;
        },
        differing);
    l_interpolate = differing != 0;
  }
  // Keep the Native domain errors even when both grids already match.
  CLUBB_levels_within_LES_domain(workspace, faverage, out_heights, k_lowest_input, k_highest_input);
  InputIndices exact_lev_idx, lower_lev_idx, upper_lev_idx;
  Bool1D l_lin_int;
  if (l_interpolate) {
    const int levels = k_highest_input - k_lowest_input + 1;
    exact_lev_idx =
        scope.view<InputIndices>(Kokkos::WithoutInitializing, "stat_file_average::exact", levels);
    lower_lev_idx =
        scope.view<InputIndices>(Kokkos::WithoutInitializing, "stat_file_average::lower", levels);
    upper_lev_idx =
        scope.view<InputIndices>(Kokkos::WithoutInitializing, "stat_file_average::upper", levels);
    l_lin_int =
        scope.view<Bool1D>(Kokkos::WithoutInitializing, "stat_file_average::linear", levels);
    // Pure bracket searches are unnecessary when the grids already match.
    LES_grid_to_CLUBB_grid(faverage, out_heights, k_lowest_input, k_highest_input, exact_lev_idx,
                           lower_lev_idx, upper_lev_idx, l_lin_int);
  }
  const auto file_z = faverage.z;
  // Native TRANSFER(2139095040,1.0_sp) is positive float infinity, not NaN.
  const real_t PosInf = real_t(std::numeric_limits<float>::infinity());
  for (int t = t1; t <= t2; ++t) {
    const bool l_convert_to_MKS = true;
    get_netcdf_var(workspace, faverage, variable_name, t, l_convert_to_MKS, file_variable, l_error);
    if (l_error) {
      std::cerr << "stat_file_average: get_var failed for " << variable_name << " in " << filename
                << " at time=" << t << '\n';
      return;
    }
    Kokkos::parallel_for(
        "stat_file_average::interpolate_accumulate", out_nz, KOKKOS_LAMBDA(int k) {
          real_t interp_variable;
          if (!l_interpolate) {
            interp_variable = file_variable(k);
          } else if (k < k_lowest_input) {
            // Native linear lower extension is disabled; retain only the ghost option.
            if (k == 0 && l_spec_bound_cond)
              interp_variable = file_variable(0);
            else
              interp_variable = PosInf;
          } else if (k > k_highest_input) {
            // Native constant upper extension is disabled too.
            interp_variable = PosInf;
          } else {
            const int index = k - k_lowest_input;
            if (l_lin_int(index))
              interp_variable = lin_interpolate_two_points(
                  out_heights(k), file_z(upper_lev_idx(index)), file_z(lower_lev_idx(index)),
                  file_variable(upper_lev_idx(index)), file_variable(lower_lev_idx(index)));
            else
              interp_variable = file_variable(exact_lev_idx(index));
          }
          // Apply an exponent (for the tuner), after the same interpolation.
          if (npower != 1)
            result(k) = result(k) + Kokkos::pow(interp_variable, npower);
          else
            result(k) = result(k) + interp_variable;
        });
  }
  // Preserve explicit close placement, including Native's earlier error returns.
  close_netcdf_read(faverage);
  Kokkos::parallel_for(
      "stat_file_average::divide", out_nz,
      KOKKOS_LAMBDA(int k) { result(k) = result(k) / real_t(num_timesteps); });
}

Real1D stat_file_average(ScratchWorkspace& workspace,             // InOut
                         const std::string& filename, int out_nz, // In
                         int t1, int t2, ConstInputProfile out_heights,
                         const std::string& variable_name,                  // In
                         int npower, bool l_spec_bound_cond, bool& l_error, // In, In, Out
                         const inputfields_state& inputfields) {            // Native module In
  // Return independent owned storage; temporary file buffers use the case frame.
  Real1D result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "stat_file_average"), out_nz);
  stat_file_average(workspace, filename, out_nz, t1, t2, out_heights, variable_name, npower,
                    l_spec_bound_cond, l_error, inputfields, result);
  return result;
}

void stat_file_average_interval(ScratchWorkspace& workspace, const std::string& filename, int nz,
                                const std::vector<int>& t, // In
                                const std::string& variable_name, ConstInputProfile out_heights,
                                int npower, // In
                                bool& l_error, const inputfields_state& inputfields,
                                InputProfile result) { // Out, Module In, Out
  // Read several Native record intervals and combine them in the source order.
  // t contains valid record pairs, terminated by zero; tdim remains their count.
  auto scope         = workspace.scope();
  constexpr int tmax = std::numeric_limits<int>::max();
  l_error            = false;
  if (t.size() > std::size_t(tmax) || t.size() < 2) {
    std::cerr << "stat_file_average_interval: Invalid time interval\n";
    l_error = true;
    return;
  }
  int tdim = -1;
  for (int i = 0; i < int(t.size()); ++i) {
    if (t[i] == 0)
      break;
    tdim = i + 1;
  }
  bool l_spec_bound_cond = false;
  stat_file_average(workspace, filename, nz, t[0], t[1], out_heights, variable_name, npower,
                    l_spec_bound_cond, l_error, inputfields, result);
  const int first_interval = t[1] - t[0];
  Kokkos::parallel_for(
      "stat_file_average_interval::first", nz,
      KOKKOS_LAMBDA(int k) { result(k) = result(k) * real_t(first_interval); });
  int divisor = t[1] - t[0];
  if (l_error)
    return;
  auto stat_file_temp =
      scope.view<Real1D>(Kokkos::WithoutInitializing, "stat_file_average_interval::temp", nz);
  for (int i = 2; i < tdim; i += 2) {
    l_spec_bound_cond = false;
    stat_file_average(workspace, filename, nz, t[i], t[i + 1], out_heights, variable_name, npower,
                      l_spec_bound_cond, l_error, inputfields, stat_file_temp);
    const int interval = t[i + 1] - t[i];
    Kokkos::parallel_for(
        "stat_file_average_interval::add", nz,
        KOKKOS_LAMBDA(int k) { result(k) = result(k) + stat_file_temp(k) * real_t(interval); });
    divisor = divisor + (t[i + 1] - t[i]);
  }
  Kokkos::parallel_for(
      "stat_file_average_interval::divide", nz,
      KOKKOS_LAMBDA(int k) { result(k) = result(k) / real_t(divisor); });
}

Real1D stat_file_average_interval(ScratchWorkspace& workspace, const std::string& filename, int nz,
                                  const std::vector<int>& t, // In
                                  const std::string& variable_name, ConstInputProfile out_heights,
                                  int npower, // In
                                  bool& l_error,
                                  const inputfields_state& inputfields) { // Out, Native module In
  // Return independent owned storage through the same interval accumulator.
  Real1D result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "stat_file_average_interval"), nz);
  stat_file_average_interval(workspace, filename, nz, t, variable_name, out_heights, npower,
                             l_error, inputfields, result);
  return result;
}

int stat_file_num_vertical_levels(ScratchWorkspace& workspace, const std::string& varname,
                                  const std::string& filename,            // In
                                  const inputfields_state& inputfields) { // Native module In
  // Return the number of input levels, converting the inclusive last slot to a count.
  stat_file file_nz;
  bool l_error;
  open_netcdf_read(workspace, varname, filename, file_nz, l_error, inputfields);
  if (l_error) {
    std::cerr << "Error opening " << filename << '\n';
    throw std::runtime_error("STOP in stat_file_num_vertical_levels");
  }
  const int result = file_nz.iz + 1;
  close_netcdf_read(file_nz);
  return result;
}

void stat_file_vertical_levels(ScratchWorkspace& workspace, const std::string& varname,
                               const std::string& filename, // In
                               int nz, const inputfields_state& inputfields,
                               InputProfile result) { // In, Module In, Out
  // Copy the first nz input heights [m] into the caller's existing profile.
  stat_file file_levels;
  bool l_error;
  open_netcdf_read(workspace, varname, filename, file_levels, l_error, inputfields);
  if (l_error) {
    std::cerr << "Error opening " << filename << '\n';
    throw std::runtime_error("STOP in stat_file_vertical_levels");
  }
  const auto levels = file_levels.z;
  Kokkos::parallel_for(
      "stat_file_vertical_levels::copy", nz, KOKKOS_LAMBDA(int k) { result(k) = levels(k); });
  close_netcdf_read(file_levels);
}

Real1D stat_file_vertical_levels(ScratchWorkspace& workspace, const std::string& varname,
                                 const std::string& filename, // In
                                 int nz,
                                 const inputfields_state& inputfields) { // In, Native module In
  // Return an independently owned copy of the first nz input heights [m].
  Real1D result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "stat_file_vertical_levels"), nz);
  stat_file_vertical_levels(workspace, varname, filename, nz, inputfields, result);
  return result;
}

bool l_netcdf_file(const std::string& filename) {
  // Return true for the Native .nc suffix; stop for an unrecognized file type.
  const std::size_t length = filename.find_last_not_of(' ') + 1;
  if (length >= 3 && filename.substr(length - 3, 3) == ".nc")
    return true;
  std::cerr << "Unrecognized file type: " << filename.substr(0, length) << '\n';
  throw std::runtime_error("STOP in l_netcdf_file");
}

void CLUBB_levels_within_LES_domain(ScratchWorkspace& workspace, const stat_file& fread_var,
                                    ConstInputProfile CLUBB_grid,                // In
                                    int& k_lowest_input, int& k_highest_input) { // Out
  // Find inclusive zero-based CLUBB bounds in the LES domain. Preserve the
  // descending-lowest then ascending-highest search and Native fatal-check order.
  auto scope   = workspace.scope();
  const auto z = fread_var.z;
  const int ia = fread_var.ia, iz = fread_var.iz, kmax = CLUBB_grid.extent_int(0) - 1;
  auto bounds = scope.view<InputIndices>(Kokkos::WithoutInitializing,
                                         "CLUBB_levels_within_LES_domain::bounds", 3);
  Kokkos::parallel_for(
      "CLUBB_levels_within_LES_domain", 1, KOKKOS_LAMBDA(int) {
        bounds(0) = -1;
        bounds(1) = -1;
        bounds(2) = 0;
        int k     = kmax;
        while (true) {
          if (CLUBB_grid(k) < z(ia)) {
            if (k == kmax) {
              bounds(2) = 1; // Bottom of LES domain is above top of CLUBB domain.
              return;
            } else {
              bounds(0) = k + 1;
              break;
            }
          } else if (k == 0) {
            bounds(0) = 0;
            break;
          } else {
            --k;
          }
        }
        k = 0;
        while (true) {
          if (CLUBB_grid(k) > z(iz)) {
            if (k == 0) {
              bounds(2) = 2; // Top of LES domain is below bottom of CLUBB domain.
              return;
            } else {
              bounds(1) = k - 1;
              break;
            }
          } else if (k == kmax) {
            bounds(1) = kmax;
            break;
          } else {
            ++k;
          }
        }
      });
  const auto host = to_host_space(workspace, bounds);
  if (host(2) == 1) {
    std::cerr << "The lowest LES input level is above the top of the CLUBB model domain.\n";
    throw std::runtime_error("Error in CLUBB_levels_within_LES_domain");
  }
  if (host(2) == 2) {
    std::cerr << "The highest LES input level is below the bottom of the CLUBB model domain.\n";
    throw std::runtime_error("Error in CLUBB_levels_within_LES_domain");
  }
  k_lowest_input  = host(0);
  k_highest_input = host(1);
}

void LES_grid_to_CLUBB_grid(const stat_file& fread_var, ConstInputProfile CLUBB_grid,
                            int first_level, int last_level, // In
                            InputIndices exact_lev_idx, InputIndices lower_lev_idx,
                            InputIndices upper_lev_idx, // Out
                            Bool1D l_lin_int) {         // Out
  // Search independent CLUBB levels in one launch. Each work item retains the
  // Native increasing LES scan, exact-match test, unset slots and early exits.
  // Output slot zero corresponds to first_level; all bounds are inclusive.
  const auto z = fread_var.z;
  const int ia = fread_var.ia, iz = fread_var.iz;
  Kokkos::parallel_for(
      "LES_grid_to_CLUBB_grid::levels", last_level - first_level + 1, KOKKOS_LAMBDA(int index) {
        const int k = first_level + index;
        int exact = -1, lower = -1, upper = -1;
        bool linear = false;
        for (int j = ia; j <= iz; ++j) {
          if (Kokkos::abs(z(j) - CLUBB_grid(k)) < real_t(1e-6)) {
            exact  = j;
            linear = false;
          } else if (z(j) < CLUBB_grid(k)) {
            lower = j;
          } else {
            upper = j;
          }
          if (exact >= 0)
            break;
          if (upper == lower + 1) {
            linear = true;
            break;
          }
        }
        exact_lev_idx(index) = exact;
        lower_lev_idx(index) = lower;
        upper_lev_idx(index) = upper;
        l_lin_int(index)     = linear;
      });
}

void LES_grid_to_CLUBB_grid(ScratchWorkspace& workspace, const stat_file& fread_var,
                            ConstInputProfile CLUBB_grid, int k, // In
                            int& exact_lev_idx, int& lower_lev_idx, int& upper_lev_idx,
                            bool& l_lin_int) { // Out
  // Preserve scalar host results for isolated queries through the batched scanner.
  // Model loops pass existing result arrays and avoid this host boundary entirely.
  auto scope  = workspace.scope();
  auto exact  = scope.view<InputIndices>(Kokkos::WithoutInitializing, "LES_grid::exact", 1);
  auto lower  = scope.view<InputIndices>(Kokkos::WithoutInitializing, "LES_grid::lower", 1);
  auto upper  = scope.view<InputIndices>(Kokkos::WithoutInitializing, "LES_grid::upper", 1);
  auto linear = scope.view<Bool1D>(Kokkos::WithoutInitializing, "LES_grid::linear", 1);
  LES_grid_to_CLUBB_grid(fread_var, CLUBB_grid, k, k, exact, lower, upper, linear);
  exact_lev_idx = to_host_space(workspace, exact)(0);
  lower_lev_idx = to_host_space(workspace, lower)(0);
  upper_lev_idx = to_host_space(workspace, upper)(0);
  l_lin_int     = to_host_space(workspace, linear)(0) != 0;
}
} // namespace clubb_kokkos
