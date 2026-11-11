// Collect named model statistics and write one NetCDF file per case.
//
// Initialize with stats_init_api, bracket each model step with
// stats_begin_timestep_api/stats_end_timestep_api, then explicitly finalize.
// stats_update adds registered scalar/surface/profile samples. The three
// budget APIs subtract the initial state, add contributions, and add the final
// state. var_on_stats_list queries registration without changing the cache.
// stats_reset_api/start_next_stats_batch_api reuse a file for successive column
// batches. stats_update_grid refreshes adaptive remapping inputs;
// stats_get_source_grid_api returns copies of the current source heights.
// The three stats_lh_samples_* APIs define/write optional SILHS subcolumns.
//
// Sampling uses the Native (tstart,tend] window and one-based model itime.
// Explicit C++ column/level indices are zero-based. Persistent accumulation
// arrays use column-first LayoutLeft HostSpace storage and Native stat_rknd
// (double), independently of core precision. Device/strided inputs stage only
// after the Native sampling/registration/error gates. An empty output path
// keeps the same registry/buffers without NetCDF writes. Snapshots deep-copy
// values but borrow file identifiers and cannot mutate or close the live file.
// Transient snapshots and file buffers reuse the private case workspace;
// finalization releases it after service frames and blocking transfers complete.
#pragma once

#include "CLUBB_core/array_types.hpp"
#include "scratch_workspace.hpp"
#include "CLUBB_core/err_info_type_module.hpp"
#include "CLUBB_core/grid_class.hpp"

#include <algorithm>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace clubb_kokkos {

// Native clubb_precision: stat_rknd/time_precision=selected_real_kind(p=12).
using stat_t          = double;
using StatsReal1D     = Kokkos::View<stat_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using StatsReal2D     = Kokkos::View<stat_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using StatsHostReal1D = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using StatsHostReal2D = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using StatsHostReal3D = Kokkos::View<real_t***, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using StatsHostReal4D = Kokkos::View<real_t****, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using StatsHostInt1D  = Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using StatsHostInt2D  = Kokkos::View<int**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using StatsHostInt3D  = Kokkos::View<int***, Kokkos::LayoutLeft, Kokkos::HostSpace>;

namespace stats_netcdf_detail {
inline constexpr int NAME_LEN            = 64;
inline constexpr int GRID_LEN            = 16;
inline constexpr int UNITS_LEN           = 32;
inline constexpr int DESC_LEN            = 128;
inline constexpr int REG_LINE_LEN        = 512;
inline constexpr int NML_REG_MAX_ENTRIES = 12000;
inline constexpr int GRID_ZT             = 1;
inline constexpr int GRID_ZM             = 2;
inline constexpr int GRID_RAD_ZT         = 3;
inline constexpr int GRID_RAD_ZM         = 4;
inline constexpr int GRID_SFC            = 5;
inline constexpr int GRID_LH_ZT          = 6;
inline constexpr int GRID_LH_SFC         = 7;
} // namespace stats_netcdf_detail

// One tracked variable; string assignments retain Native fixed-width truncation.
struct stats_var {
  stats_var() = default;
  stats_var(const stats_var& source);
  stats_var& operator=(const stats_var& source);
  stats_var(stats_var&&) noexcept            = default;
  stats_var& operator=(stats_var&&) noexcept = default;
  std::string name{};      // Name of variable [NAME_LEN].
  std::string grid{};      // Grid token [GRID_LEN].
  std::string units{};     // Variable unit string [UNITS_LEN].
  std::string long_name{}; // Description [DESC_LEN].
  int grid_id      = stats_netcdf_detail::GRID_ZT;
  int varid        = -1;
  int nz           = 0;
  int out_nz       = 0;
  bool l_budget    = false;
  bool l_in_budget = false;
  StatsReal2D buffer;      // (column,level) accumulated values.
  StatsHostInt2D nsamples; // Pointwise sample counts.
  StatsReal2D sub_timestep_sum;
  StatsHostInt2D sub_timestep_nsamples;
  bool l_sub_timestep_active     = false;
  bool l_regular_timestep_active = false;
};

struct stats_grid_module {
  bool is_initialized   = false;
  int grid_remap_method = 0;
  clubb_kokkos::grid gr_source;
  clubb_kokkos::grid gr_output;
  StatsHostReal2D rho_lin_spline_vals;   // (column,index) density [kg m-3].
  StatsHostReal2D rho_lin_spline_levels; // (column,index) spline heights [m].
  StatsHostReal1D p_sfc;                 // Surface pressure [Pa].
};

struct stats_lh_2d {
  bool is_initialized = false;
  int sample_dimid    = -1;
  int lh_zt_dimid     = -1;
  int lh_zt_varid     = -1;
  int num_samples     = 0;
  int nzt             = 0;
  int n_nl_vars       = 0;
  int n_u_vars        = 0;
  StatsHostInt1D nl_varids;
  StatsHostInt1D u_varids;
  std::vector<std::string> nl_names;
  std::vector<std::string> u_names;
};

struct stats_varname_lookup_cache {
  int head             = 0; // Zero-based expected-hit cursor.
  int cache_len        = 0;
  int reject_head      = 0; // Zero-based expected-miss cursor.
  int reject_cache_len = 0;
  StatsHostInt1D cache; // Zero-based registry slots; -1 is unset.
  std::vector<std::string> reject_cache;
};

struct stats_type {
  bool enabled                  = false;
  bool l_sample                 = false;
  bool l_last_sample            = false;
  bool l_netcdf_output          = true;
  bool l_different_output_grid  = false;
  bool l_output_rad_files       = false;
  bool l_open_ended_time_window = false;
  int time_index                = 0;
  int nvars                     = 0;
  int ncol_batch                = 0;
  int ncol_total                = 0;
  int active_batch_offset       = 0;
  int stats_nsamp               = 1;
  int stats_nout                = 1;
  int samples_per_write         = 1;
  double dt_main                = 0;
  double time_initial           = 0;
  double tstart                 = 0;
  double tend                   = 0;
  std::vector<stats_var> vars;
  stats_grid_module grid;
  stats_lh_2d lh_2d;
  stats_varname_lookup_cache lookup;
  int ncid               = -1;
  int time_dimid         = -1;
  int bnds_dimid         = -1;
  int time_varid         = -1;
  int time_bnds_varid    = -1;
  int col_dimid          = -1;
  int col_varid          = -1;
  int zt_dimid           = -1;
  int zm_dimid           = -1;
  int rad_zt_dimid       = -1;
  int rad_zm_dimid       = -1;
  int sfc_dimid          = -1;
  int zt_varid           = -1;
  int zm_varid           = -1;
  int rad_zt_varid       = -1;
  int rad_zm_varid       = -1;
  int param_strlen_dimid = -1;
  int param_dimid        = -1;
  int param_varid        = -1;
  int param_name_varid   = -1;
  int clubb_params_varid = -1;

  // Language ownership metadata; vectors cannot distinguish unallocated/size0.
  bool vars_allocated = false;
  // Native's process-SAVE warning is case-owned here: warn once per context,
  // retaining this flag across resets/reinitialization without global storage.
  bool l_warned_overflow = false;
  stats_type()           = default;
  ~stats_type();
  stats_type(const stats_type&)            = delete;
  stats_type& operator=(const stats_type&) = delete;
  void require_mutable() const;
  ScratchWorkspace& workspace();
  friend void copy_stats_snapshot(const stats_type&, stats_type&);

private:
  ScratchWorkspace workspace_; // Case-owned transient buffers; excluded from snapshots.
  bool snapshot_ = false;
};

void copy_stats_snapshot(const stats_type& source, stats_type& destination);

// Optional Native initialization arguments grouped because C++ has no keywords.
struct StatsInitOptions {
  std::optional<StatsHostReal2D> clubb_params;
  std::optional<std::vector<std::string>> param_names;
  std::optional<StatsHostReal1D> rad_zt;
  std::optional<StatsHostReal1D> rad_zm;
  std::optional<std::vector<std::string>> hydromet_list;
  std::optional<int> sclr_dim;
  std::optional<int> edsclr_dim;
  std::optional<StatsHostReal1D> output_zt;
  std::optional<StatsHostReal1D> output_zm;
  std::optional<int> grid_remap_method;
  std::optional<int> ncol_total;
  std::optional<double> stats_tstart;
  std::optional<double> stats_tend;
};

void stats_init_api(const std::string& registry_path, const std::string& output_path,
                    int ncol_batch, real_t stats_tsamp, real_t stats_tout, real_t dt_main,
                    int day_in, int month_in, int year_in, double time_initial, StatsHostReal1D zt,
                    StatsHostReal1D zm, stats_type& stats, err_info_type& err_info,
                    const StatsInitOptions& options = {});
void stats_finalize_api(stats_type& stats, err_info_type& err_info);
void stats_reset_api(stats_type& stats);
void start_next_stats_batch_api(stats_type& stats, err_info_type& err_info);
void stats_get_source_grid_api(const stats_type& stats, StatsHostReal2D& zt, StatsHostReal2D& zm);
void stats_begin_timestep_api(int itime, stats_type& stats);
void stats_end_timestep_api(real_t time_value, stats_type& stats, err_info_type& err_info);
bool var_on_stats_list(const stats_type& stats, const std::string& name);

void stats_update(const std::string& name, real_t values, stats_type& stats,
                  std::optional<int> icol = std::nullopt, std::optional<int> level = std::nullopt,
                  bool sub_timestep_average = false);
void stats_begin_budget(const std::string& name, real_t values, stats_type& stats,
                        std::optional<int> icol = std::nullopt);
void stats_update_budget(const std::string& name, real_t values, stats_type& stats,
                         std::optional<int> icol  = std::nullopt,
                         std::optional<int> level = std::nullopt);
void stats_finalize_budget(const std::string& name, real_t values, stats_type& stats,
                           std::optional<int> icol = std::nullopt, bool l_count_sample = true);
void stats_update_grid(StatsHostReal2D zt_src, StatsHostReal2D zm_src, StatsHostReal2D rho_vals,
                       StatsHostReal2D rho_levels, StatsHostReal1D p_sfc, stats_type& stats);
void stats_lh_samples_init(int num_samples, int nzt, const std::vector<std::string>& nl_var_names,
                           const std::vector<std::string>& u_var_names, StatsHostReal1D zt_vals,
                           stats_type& stats, err_info_type& err_info);
void stats_lh_samples_write_lognormal(StatsHostReal4D samples, stats_type& stats,
                                      err_info_type& err_info);
void stats_lh_samples_write_uniform(StatsHostReal4D uniform_vals, StatsHostInt3D mixture_comp,
                                    StatsHostReal3D sample_weights, stats_type& stats,
                                    err_info_type& err_info);

namespace stats_netcdf_detail {
// Native private routines used by the compiled bodies and generic templates.
void stats_type_initialize(int ncol_batch, const std::vector<stats_var>& defs, int ndefs,
                           StatsHostReal1D zt, StatsHostReal1D zm, stats_type& stats, int& ierr,
                           std::optional<StatsHostReal1D> rad_zt = std::nullopt,
                           std::optional<StatsHostReal1D> rad_zm = std::nullopt);
void stats_open_netcdf(const std::string& output_path, const std::string& time_units,
                       StatsHostReal1D zt, StatsHostReal1D zm, stats_type& stats, int& ierr,
                       std::optional<StatsHostReal2D> clubb_params         = std::nullopt,
                       std::optional<std::vector<std::string>> param_names = std::nullopt,
                       std::optional<StatsHostReal1D> rad_zt               = std::nullopt,
                       std::optional<StatsHostReal1D> rad_zm               = std::nullopt);
void stats_prepare_sub_timestep_average(stats_var& var);
void stats_flush_sub_timestep_averages(stats_type& stats);
void stats_avg_and_write(real_t time_value, stats_type& stats, int& ierr);
void stats_grid_init(ScratchWorkspace& stats_workspace, int ncol_batch, StatsHostReal1D zt_src,
                     StatsHostReal1D zm_src, StatsHostReal1D zt_tgt, StatsHostReal1D zm_tgt,
                     int grid_remap_method, stats_grid_module& grid_ctx, int& ierr);
int stats_get_id(stats_type& stats, const std::string& name);
void stats_lookup_append(stats_type& stats, int id);
void stats_lookup_reject_append(stats_type& stats, const std::string& name);
void format_date(int day_in, int month_in, int year_in, double time_in, std::string& date);
void stats_read_registry_namelist(ScratchWorkspace& stats_workspace, const std::string& path,
                                  std::vector<stats_var>& defs, int& ndefs, int& ierr);
void stats_expand_registry(std::optional<int> sclr_dim, std::optional<int> edsclr_dim,
                           std::vector<stats_var>& defs, int& ndefs,
                           const std::optional<std::vector<std::string>>& hydromet_list,
                           bool& l_warned_overflow);
void add_expanded_def(const stats_var& base_def, const std::string& name,
                      std::vector<stats_var>& out_defs, int& nout, bool& l_warned_overflow);
void parse_registry_line(const std::string& line, stats_var& def, bool& l_valid);
void split_registry_fields(const std::string& line, std::vector<std::string>& fields, int& nfields);
bool is_line_comment(const std::string& line);
} // namespace stats_netcdf_detail

template <class Zt, class Zm>
void stats_init_api(const std::string& registry_path, const std::string& output_path,
                    int ncol_batch, real_t stats_tsamp, real_t stats_tout, real_t dt_main,
                    int day_in, int month_in, int year_in, double time_initial, const Zt& zt,
                    const Zm& zm, stats_type& stats, err_info_type& err_info,
                    const StatsInitOptions& options = {}) {
  // Only required grid vectors need staging; optional metadata is already host-owned.
  stats.require_mutable();
  auto& stats_workspace = stats.workspace();
  auto stats_scope      = stats_workspace.scope();
  stats_init_api(registry_path, output_path, ncol_batch, stats_tsamp, stats_tout, dt_main, day_in,
                 month_in, year_in, time_initial, to_host_space(stats_workspace, zt),
                 to_host_space(stats_workspace, zm), stats, err_info, options);
}

// Native generic entry bodies are templates only to retain each caller's
// memory space/layout until the original update-host boundary is reached.
namespace stats_netcdf_detail {
template <class View>
void stats_update_1d(const std::string& name, const View& values, stats_type& stats,
                     std::optional<int> icol, bool sub_timestep_average) {
  //   Add sampled 1D values to the accumulation buffer for a named variable.
  //   For profile variables with icol present this updates one column profile;
  //   without icol this updates all surface columns.
  stats.require_mutable();
  if (!stats.enabled)
    return;
  const int id = stats_get_id(stats, name);
  if (id < 0)
    return;
  if (!stats.l_sample)
    return;
  auto& var = stats.vars[id];
  if (var.l_budget) {
    std::cerr << "stats error: called stat_update for variable currently in budget: " << var.name
              << '\n';
    return;
  }
  auto& stats_workspace  = stats.workspace();
  auto stats_scope       = stats_workspace.scope();
  const auto values_host = to_host_space(stats_workspace, values); // Native update host(values).
  if (sub_timestep_average) {
    if (var.l_regular_timestep_active) {
      std::cerr << "stats error: mixed regular and sub-timestep updates for " << var.name << '\n';
      return;
    }
    stats_prepare_sub_timestep_average(var);
    var.l_sub_timestep_active = true;
    if (icol) {
      for (int k = 0; k < var.buffer.extent_int(1); ++k) {
        var.sub_timestep_sum(*icol, k) =
            var.sub_timestep_sum(*icol, k) + static_cast<real_t>(values_host(k));
        var.sub_timestep_nsamples(*icol, k) = var.sub_timestep_nsamples(*icol, k) + 1;
      }
    } else {
      for (int c = 0; c < var.buffer.extent_int(0); ++c) {
        var.sub_timestep_sum(c, 0) =
            var.sub_timestep_sum(c, 0) + static_cast<real_t>(values_host(c));
        var.sub_timestep_nsamples(c, 0) = var.sub_timestep_nsamples(c, 0) + 1;
      }
    }
    return;
  }
  if (var.l_sub_timestep_active) {
    std::cerr << "stats error: mixed sub-timestep and regular updates for " << var.name << '\n';
    return;
  }
  var.l_regular_timestep_active = true;
  if (icol) {
    for (int k = 0; k < var.buffer.extent_int(1); ++k) {
      var.buffer(*icol, k)   = var.buffer(*icol, k) + static_cast<real_t>(values_host(k));
      var.nsamples(*icol, k) = var.nsamples(*icol, k) + 1;
    }
  } else {
    for (int c = 0; c < var.buffer.extent_int(0); ++c) {
      var.buffer(c, 0)   = var.buffer(c, 0) + static_cast<real_t>(values_host(c));
      var.nsamples(c, 0) = var.nsamples(c, 0) + 1;
    }
  }
  int max_count = 0;
  for (int k = 0; k < var.nsamples.extent_int(1); ++k)
    for (int c = 0; c < var.nsamples.extent_int(0); ++c)
      max_count = std::max(max_count, var.nsamples(c, k));
  if (max_count > stats.samples_per_write) {
    std::cerr << "stats oversampling warning for " << var.name << ": max_count=" << max_count
              << " expected <=" << stats.samples_per_write << '\n';
  }
}

template <class View>
void stats_update_2d(const std::string& name, const View& values, stats_type& stats,
                     bool sub_timestep_average) {
  //   Add sampled full-field 2D values to the accumulation buffer for a named
  //   profile variable. This path updates all columns and all levels at once.
  stats.require_mutable();
  if (!stats.enabled)
    return;
  const int id = stats_get_id(stats, name);
  if (id < 0)
    return;
  if (!stats.l_sample)
    return;
  auto& var = stats.vars[id];
  if (var.l_budget) {
    std::cerr << "stats error: called stat_update for variable currently in budget: " << var.name
              << '\n';
    return;
  }
  if (values.extent(0) != var.buffer.extent(0) || values.extent(1) != var.buffer.extent(1)) {
    std::cerr << "stats shape mismatch for " << var.name << ": input nz =" << values.extent(0)
              << " expected nz =" << var.buffer.extent(0) << ": input ngrdcol =" << values.extent(1)
              << " expected ngrdcol " << var.buffer.extent(1) << '\n';
    return;
  }
  auto& stats_workspace  = stats.workspace();
  auto stats_scope       = stats_workspace.scope();
  const auto values_host = to_host_space(stats_workspace, values); // Native update host(values).
  if (sub_timestep_average) {
    if (var.l_regular_timestep_active) {
      std::cerr << "stats error: mixed regular and sub-timestep updates for " << var.name << '\n';
      return;
    }
    stats_prepare_sub_timestep_average(var);
    var.l_sub_timestep_active = true;
    for (int k = 0; k < var.buffer.extent_int(1); ++k) {
      for (int c = 0; c < var.buffer.extent_int(0); ++c) {
        var.sub_timestep_sum(c, k) =
            var.sub_timestep_sum(c, k) + static_cast<real_t>(values_host(c, k));
        var.sub_timestep_nsamples(c, k) = var.sub_timestep_nsamples(c, k) + 1;
      }
    }
    return;
  }
  if (var.l_sub_timestep_active) {
    std::cerr << "stats error: mixed sub-timestep and regular updates for " << var.name << '\n';
    return;
  }
  var.l_regular_timestep_active = true;
  for (int k = 0; k < var.buffer.extent_int(1); ++k) {
    for (int c = 0; c < var.buffer.extent_int(0); ++c) {
      var.buffer(c, k)   = var.buffer(c, k) + static_cast<real_t>(values_host(c, k));
      var.nsamples(c, k) = var.nsamples(c, k) + 1;
    }
  }
  int max_count = 0;
  for (int k = 0; k < var.nsamples.extent_int(1); ++k)
    for (int c = 0; c < var.nsamples.extent_int(0); ++c)
      max_count = std::max(max_count, var.nsamples(c, k));
  if (max_count > stats.samples_per_write) {
    std::cerr << "stats oversampling warning for " << var.name << ": max_count=" << max_count
              << " expected <=" << stats.samples_per_write << '\n';
  }
}

template <class View>
void stats_begin_budget_1d(const std::string& name, const View& values, stats_type& stats,
                           std::optional<int> icol) {
  //   Begin budget accumulation for 1D payloads by subtracting the current
  //   vector state from the buffer.
  stats.require_mutable();
  if (!stats.enabled)
    return;
  const int id = stats_get_id(stats, name);
  if (id < 0)
    return;
  if (!stats.l_sample)
    return;
  auto& var = stats.vars[id];
  if (var.l_in_budget) {
    std::cerr << "stats budget begin twice for " << var.name << '\n';
    return;
  }
  auto& stats_workspace  = stats.workspace();
  auto stats_scope       = stats_workspace.scope();
  const auto values_host = to_host_space(stats_workspace, values); // Native update host(values).
  // Subtract the current state from the buffer to begin budget accumulation.
  if (icol) {
    for (int k = 0; k < var.buffer.extent_int(1); ++k) {
      var.buffer(*icol, k) = var.buffer(*icol, k) - static_cast<real_t>(values_host(k));
    }
  } else {
    for (int c = 0; c < var.buffer.extent_int(0); ++c) {
      var.buffer(c, 0) = var.buffer(c, 0) - static_cast<real_t>(values_host(c));
    }
  }
  var.l_budget    = true;
  var.l_in_budget = true;
}

template <class View>
void stats_begin_budget_2d(const std::string& name, const View& values, stats_type& stats) {
  //   Begin budget accumulation for full-field 2D payloads by subtracting the
  //   current state from the full buffer.
  stats.require_mutable();
  if (!stats.enabled)
    return;
  const int id = stats_get_id(stats, name);
  if (id < 0)
    return;
  if (!stats.l_sample)
    return;
  auto& var = stats.vars[id];
  if (var.l_in_budget) {
    std::cerr << "stats budget begin twice for " << var.name << '\n';
    return;
  }
  auto& stats_workspace  = stats.workspace();
  auto stats_scope       = stats_workspace.scope();
  const auto values_host = to_host_space(stats_workspace, values); // Native update host(values).
  // Subtract the current state from the buffer to begin budget accumulation.
  for (int k = 0; k < var.buffer.extent_int(1); ++k) {
    for (int c = 0; c < var.buffer.extent_int(0); ++c) {
      var.buffer(c, k) = var.buffer(c, k) - static_cast<real_t>(values_host(c, k));
    }
  }
  var.l_budget    = true;
  var.l_in_budget = true;
}

template <class View>
void stats_update_budget_1d(const std::string& name, const View& values, stats_type& stats,
                            std::optional<int> icol) {
  //   Add intermediate budget contributions from 1D payloads during an active
  //   budget window.
  stats.require_mutable();
  if (!stats.enabled)
    return;
  if (!stats.l_sample)
    return;
  const int id = stats_get_id(stats, name);
  if (id < 0)
    return;
  auto& var = stats.vars[id];
  if (!var.l_in_budget) {
    std::cerr << "stats budget update without begin for " << var.name << '\n';
    return;
  }
  auto& stats_workspace  = stats.workspace();
  auto stats_scope       = stats_workspace.scope();
  const auto values_host = to_host_space(stats_workspace, values); // Native update host(values).
  if (icol) {
    for (int k = 0; k < var.buffer.extent_int(1); ++k) {
      var.buffer(*icol, k) = var.buffer(*icol, k) + static_cast<real_t>(values_host(k));
    }
  } else {
    for (int c = 0; c < var.buffer.extent_int(0); ++c) {
      var.buffer(c, 0) = var.buffer(c, 0) + static_cast<real_t>(values_host(c));
    }
  }
}

template <class View>
void stats_update_budget_2d(const std::string& name, const View& values, stats_type& stats) {
  //   Add intermediate budget contributions from full-field 2D payloads during
  //   an active budget window.
  stats.require_mutable();
  if (!stats.enabled)
    return;
  if (!stats.l_sample)
    return;
  const int id = stats_get_id(stats, name);
  if (id < 0)
    return;
  auto& var              = stats.vars[id];
  auto& stats_workspace  = stats.workspace();
  auto stats_scope       = stats_workspace.scope();
  const auto values_host = to_host_space(stats_workspace, values); // Native update host(values).
  if (!var.l_in_budget) {
    std::cerr << "stats budget update without begin for " << var.name << '\n';
    return;
  }
  for (int k = 0; k < var.buffer.extent_int(1); ++k) {
    for (int c = 0; c < var.buffer.extent_int(0); ++c) {
      var.buffer(c, k) = var.buffer(c, k) + static_cast<real_t>(values_host(c, k));
    }
  }
}

template <class View>
void stats_finalize_budget_1d(const std::string& name, const View& values, stats_type& stats,
                              std::optional<int> icol, bool l_count_sample) {
  //   Close an active budget window using 1D payloads, and optionally
  //   increment per-element sample counters for averaging.
  stats.require_mutable();
  if (!stats.enabled)
    return;
  if (!stats.l_sample)
    return;
  const int id = stats_get_id(stats, name);
  if (id < 0)
    return;
  auto& var = stats.vars[id];
  // Don't sample if l_count_sample is false; modify without counting a sample.
  const bool do_count_sample = l_count_sample;
  if (!var.l_in_budget) {
    std::cerr << "stats budget finalize without begin for " << var.name << '\n';
    return;
  }
  auto& stats_workspace  = stats.workspace();
  auto stats_scope       = stats_workspace.scope();
  const auto values_host = to_host_space(stats_workspace, values); // Native update host(values).
  if (icol) {
    for (int k = 0; k < var.buffer.extent_int(1); ++k) {
      var.buffer(*icol, k) = var.buffer(*icol, k) + static_cast<real_t>(values_host(k));
      if (do_count_sample)
        var.nsamples(*icol, k) = var.nsamples(*icol, k) + 1;
    }
  } else {
    for (int c = 0; c < var.buffer.extent_int(0); ++c) {
      var.buffer(c, 0) = var.buffer(c, 0) + static_cast<real_t>(values_host(c));
      if (do_count_sample)
        var.nsamples(c, 0) = var.nsamples(c, 0) + 1;
    }
  }
  var.l_budget    = true;
  var.l_in_budget = false;
  if (do_count_sample) {
    int max_count = 0;
    for (int k = 0; k < var.nsamples.extent_int(1); ++k)
      for (int c = 0; c < var.nsamples.extent_int(0); ++c)
        max_count = std::max(max_count, var.nsamples(c, k));
    if (max_count > stats.samples_per_write) {
      std::cerr << "stats oversampling warning for " << var.name << ": max_count=" << max_count
                << " expected <=" << stats.samples_per_write << '\n';
    }
  }
}

template <class View>
void stats_finalize_budget_2d(const std::string& name, const View& values, stats_type& stats,
                              bool l_count_sample) {
  //   Close an active budget window using full-field 2D payloads, and
  //   optionally increment per-element sample counters for averaging.
  stats.require_mutable();
  if (!stats.enabled)
    return;
  if (!stats.l_sample)
    return;
  const int id = stats_get_id(stats, name);
  if (id < 0)
    return;
  auto& var = stats.vars[id];
  // Don't sample if l_count_sample is false; modify without counting a sample.
  const bool do_count_sample = l_count_sample;
  if (!var.l_in_budget) {
    std::cerr << "stats budget finalize without begin for " << var.name << '\n';
    return;
  }
  auto& stats_workspace  = stats.workspace();
  auto stats_scope       = stats_workspace.scope();
  const auto values_host = to_host_space(stats_workspace, values); // Native update host(values).
  for (int k = 0; k < var.buffer.extent_int(1); ++k) {
    for (int c = 0; c < var.buffer.extent_int(0); ++c) {
      var.buffer(c, k) = var.buffer(c, k) + static_cast<real_t>(values_host(c, k));
      if (do_count_sample)
        var.nsamples(c, k) = var.nsamples(c, k) + 1;
    }
  }
  var.l_budget    = true;
  var.l_in_budget = false;
  if (do_count_sample) {
    int max_count = 0;
    for (int k = 0; k < var.nsamples.extent_int(1); ++k)
      for (int c = 0; c < var.nsamples.extent_int(0); ++c)
        max_count = std::max(max_count, var.nsamples(c, k));
    if (max_count > stats.samples_per_write) {
      std::cerr << "stats oversampling warning for " << var.name << ": max_count=" << max_count
                << " expected <=" << stats.samples_per_write << '\n';
    }
  }
}

} // namespace stats_netcdf_detail

template <class View, std::enable_if_t<Kokkos::is_view_v<View> && View::rank == 1, int> = 0>
void stats_update(const std::string& name, const View& values, stats_type& stats,
                  std::optional<int> icol = std::nullopt, bool sub_timestep_average = false) {
  // Dispatch the Native generic without changing its optional argument presence.
  stats_netcdf_detail::stats_update_1d(name, values, stats, icol, sub_timestep_average);
}

template <class View, std::enable_if_t<Kokkos::is_view_v<View> && View::rank == 2, int> = 0>
void stats_update(const std::string& name, const View& values, stats_type& stats,
                  bool sub_timestep_average = false) {
  // Dispatch the Native generic without changing its optional argument presence.
  stats_netcdf_detail::stats_update_2d(name, values, stats, sub_timestep_average);
}

template <class View, std::enable_if_t<Kokkos::is_view_v<View> && View::rank == 1, int> = 0>
void stats_begin_budget(const std::string& name, const View& values, stats_type& stats,
                        std::optional<int> icol = std::nullopt) {
  // Dispatch the Native generic without changing its optional argument presence.
  stats_netcdf_detail::stats_begin_budget_1d(name, values, stats, icol);
}

template <class View, std::enable_if_t<Kokkos::is_view_v<View> && View::rank == 2, int> = 0>
void stats_begin_budget(const std::string& name, const View& values, stats_type& stats) {
  // Dispatch the Native generic without changing its optional argument presence.
  stats_netcdf_detail::stats_begin_budget_2d(name, values, stats);
}

template <class View, std::enable_if_t<Kokkos::is_view_v<View> && View::rank == 1, int> = 0>
void stats_update_budget(const std::string& name, const View& values, stats_type& stats,
                         std::optional<int> icol = std::nullopt) {
  // Dispatch the Native generic without changing its optional argument presence.
  stats_netcdf_detail::stats_update_budget_1d(name, values, stats, icol);
}

template <class View, std::enable_if_t<Kokkos::is_view_v<View> && View::rank == 2, int> = 0>
void stats_update_budget(const std::string& name, const View& values, stats_type& stats) {
  // Dispatch the Native generic without changing its optional argument presence.
  stats_netcdf_detail::stats_update_budget_2d(name, values, stats);
}

template <class View, std::enable_if_t<Kokkos::is_view_v<View> && View::rank == 1, int> = 0>
void stats_finalize_budget(const std::string& name, const View& values, stats_type& stats,
                           std::optional<int> icol = std::nullopt, bool l_count_sample = true) {
  // Dispatch the Native generic without changing its optional argument presence.
  stats_netcdf_detail::stats_finalize_budget_1d(name, values, stats, icol, l_count_sample);
}

template <class View, std::enable_if_t<Kokkos::is_view_v<View> && View::rank == 2, int> = 0>
void stats_finalize_budget(const std::string& name, const View& values, stats_type& stats,
                           bool l_count_sample = true) {
  // Dispatch the Native generic without changing its optional argument presence.
  stats_netcdf_detail::stats_finalize_budget_2d(name, values, stats, l_count_sample);
}

template <class Zt, class Zm, class Rho, class Levels, class Surface>
void stats_update_grid(const Zt& zt_src, const Zm& zm_src, const Rho& rho_vals,
                       const Levels& rho_levels, const Surface& p_sfc, stats_type& stats) {
  // Native early gates precede the device-to-host remapping metadata transfer.
  stats.require_mutable();
  if (!stats.enabled || !stats.l_different_output_grid || !stats.grid.is_initialized)
    return;
  if (zt_src.extent_int(0) != stats.ncol_batch || zm_src.extent_int(0) != stats.ncol_batch ||
      rho_vals.extent_int(0) != stats.ncol_batch || rho_levels.extent_int(0) != stats.ncol_batch ||
      p_sfc.extent_int(0) != stats.ncol_batch)
    return;
  if (zt_src.extent_int(1) != zm_src.extent_int(1) - 1 ||
      rho_levels.extent_int(1) != rho_vals.extent_int(1))
    return;
  auto& stats_workspace = stats.workspace();
  auto stats_scope      = stats_workspace.scope();
  stats_update_grid(to_host_space(stats_workspace, zt_src), to_host_space(stats_workspace, zm_src),
                    to_host_space(stats_workspace, rho_vals),
                    to_host_space(stats_workspace, rho_levels),
                    to_host_space(stats_workspace, p_sfc), stats);
}

template <class Heights>
void stats_lh_samples_init(int num_samples, int nzt, const std::vector<std::string>& nl_var_names,
                           const std::vector<std::string>& u_var_names, const Heights& zt_vals,
                           stats_type& stats, err_info_type& err_info) {
  // Define/write coordinates only for an active NetCDF context.
  stats.require_mutable();
  if (!stats.enabled || !stats.l_netcdf_output || num_samples <= 0 || nzt <= 0 ||
      stats.zt_dimid < 0)
    return;
  auto& stats_workspace = stats.workspace();
  auto stats_scope      = stats_workspace.scope();
  stats_lh_samples_init(num_samples, nzt, nl_var_names, u_var_names,
                        to_host_space(stats_workspace, zt_vals), stats, err_info);
}

template <class Samples>
void stats_lh_samples_write_lognormal(const Samples& samples, stats_type& stats,
                                      err_info_type& err_info) {
  // Native sample-output gates precede staging; values retain col/sample/level/var order.
  stats.require_mutable();
  if (!stats.enabled || !stats.l_netcdf_output || !stats.lh_2d.is_initialized ||
      !stats.l_last_sample || stats.lh_2d.n_nl_vars <= 0)
    return;
  if (samples.extent_int(0) != stats.ncol_batch ||
      samples.extent_int(1) != stats.lh_2d.num_samples ||
      samples.extent_int(2) != stats.lh_2d.nzt || samples.extent_int(3) != stats.lh_2d.n_nl_vars)
    return;
  auto& stats_workspace = stats.workspace();
  auto stats_scope      = stats_workspace.scope();
  stats_lh_samples_write_lognormal(to_host_space(stats_workspace, samples), stats, err_info);
}

template <class Uniform, class Mixture, class Weights>
void stats_lh_samples_write_uniform(const Uniform& uniform_vals, const Mixture& mixture_comp,
                                    const Weights& sample_weights, stats_type& stats,
                                    err_info_type& err_info) {
  // Native uniform-output shape checks also avoid unnecessary device transfers.
  stats.require_mutable();
  if (!stats.enabled || !stats.l_netcdf_output || !stats.lh_2d.is_initialized ||
      !stats.l_last_sample || stats.lh_2d.n_u_vars <= 0)
    return;
  if (uniform_vals.extent_int(0) != stats.ncol_batch ||
      uniform_vals.extent_int(1) != stats.lh_2d.num_samples ||
      uniform_vals.extent_int(2) != stats.lh_2d.nzt ||
      mixture_comp.extent_int(0) != stats.ncol_batch ||
      mixture_comp.extent_int(1) != stats.lh_2d.num_samples ||
      mixture_comp.extent_int(2) != stats.lh_2d.nzt ||
      sample_weights.extent_int(0) != stats.ncol_batch ||
      sample_weights.extent_int(1) != stats.lh_2d.num_samples ||
      sample_weights.extent_int(2) != stats.lh_2d.nzt ||
      stats.lh_2d.n_u_vars != uniform_vals.extent_int(3) + 2)
    return;
  auto& stats_workspace = stats.workspace();
  auto stats_scope      = stats_workspace.scope();
  stats_lh_samples_write_uniform(to_host_space(stats_workspace, uniform_vals),
                                 to_host_space(stats_workspace, mixture_comp),
                                 to_host_space(stats_workspace, sample_weights), stats, err_info);
}

} // namespace clubb_kokkos
