// Collect and save CLUBB statistics; see stats_netcdf.hpp for the public API.
//
// Mirrors src/CLUBB_core/stats_netcdf.F90 in routine order. Accumulation and
// counters remain host-owned; NetCDF-C receives Native-ordered definitions and
// packed slabs. Existing C++ grid/remapping routines handle output-grid changes.
// File ownership stays with the case; read-only snapshots own independent data.
#include "stats_netcdf.hpp"
#include "CLUBB_core/error_code.hpp"
#include "netcdf.hpp"
#include <cmath>
#include <limits>
#include <utility>
#include "CLUBB_core/calendar.hpp"
#include "CLUBB_core/remapping_module.hpp"
#include "namelist.hpp"

namespace clubb_kokkos {
namespace {

template <class View> View stats_value_copy(const View& source) {
  // Native allocatable assignment copies values, whereas a Kokkos View aliases.
  View destination;
  if (!source.is_allocated())
    return destination;
  if constexpr (View::rank == 1)
    destination =
        View(Kokkos::view_alloc(Kokkos::WithoutInitializing, source.label()), source.extent(0));
  else if constexpr (View::rank == 2)
    destination = View(Kokkos::view_alloc(Kokkos::WithoutInitializing, source.label()),
                       source.extent(0), source.extent(1));
  else if constexpr (View::rank == 3)
    destination = View(Kokkos::view_alloc(Kokkos::WithoutInitializing, source.label()),
                       source.extent(0), source.extent(1), source.extent(2));
  Kokkos::deep_copy(destination, source);
  return destination;
}

} // namespace

stats_var::stats_var(const stats_var& source) {
  // Copy Native allocatable values rather than sharing Kokkos allocations.
  *this = source;
}

stats_var& stats_var::operator=(const stats_var& source) {
  // Reproduce Native derived-type assignment, including allocatable buffers.
  if (this == &source)
    return *this;
  name                      = source.name;
  grid                      = source.grid;
  units                     = source.units;
  long_name                 = source.long_name;
  grid_id                   = source.grid_id;
  varid                     = source.varid;
  nz                        = source.nz;
  out_nz                    = source.out_nz;
  l_budget                  = source.l_budget;
  l_in_budget               = source.l_in_budget;
  buffer                    = stats_value_copy(source.buffer);
  nsamples                  = stats_value_copy(source.nsamples);
  sub_timestep_sum          = stats_value_copy(source.sub_timestep_sum);
  sub_timestep_nsamples     = stats_value_copy(source.sub_timestep_nsamples);
  l_sub_timestep_active     = source.l_sub_timestep_active;
  l_regular_timestep_active = source.l_regular_timestep_active;
  return *this;
}

stats_type::~stats_type() {
  // Explicit finalization reports errors; destruction also closes partial owners.
  if (snapshot_)
    return;
  err_info_type err_info;
  stats_finalize_api(*this, err_info);
}

void stats_type::require_mutable() const {
  // Snapshots borrow file identifiers and may only participate in read queries.
  if (snapshot_)
    throw std::logic_error("A stats snapshot cannot update live I/O");
}

ScratchWorkspace& stats_type::workspace() {
  // Supply case-owned temporary storage only to mutable service operations.
  require_mutable();
  return workspace_;
}

void copy_stats_snapshot(const stats_type& source, stats_type& destination) {
  // A Native value snapshot owns copied data and never owns the live file.
  if (&source == &destination)
    throw std::logic_error("Cannot snapshot stats onto itself");
  if (!destination.snapshot_ && (destination.enabled || destination.ncid >= 0))
    throw std::logic_error("Cannot overwrite an active stats owner with a snapshot");
  // Snapshot payloads own values; temporary pools are neither copied nor borrowed.
  destination.workspace_.clear();
  // Mark borrowing before any allocation can throw while copying the values.
  destination.snapshot_                    = true;
  destination.enabled                      = source.enabled;
  destination.l_sample                     = source.l_sample;
  destination.l_last_sample                = source.l_last_sample;
  destination.l_netcdf_output              = source.l_netcdf_output;
  destination.l_different_output_grid      = source.l_different_output_grid;
  destination.l_output_rad_files           = source.l_output_rad_files;
  destination.l_open_ended_time_window     = source.l_open_ended_time_window;
  destination.time_index                   = source.time_index;
  destination.nvars                        = source.nvars;
  destination.ncol_batch                   = source.ncol_batch;
  destination.ncol_total                   = source.ncol_total;
  destination.active_batch_offset          = source.active_batch_offset;
  destination.stats_nsamp                  = source.stats_nsamp;
  destination.stats_nout                   = source.stats_nout;
  destination.samples_per_write            = source.samples_per_write;
  destination.dt_main                      = source.dt_main;
  destination.time_initial                 = source.time_initial;
  destination.tstart                       = source.tstart;
  destination.tend                         = source.tend;
  destination.vars                         = source.vars;
  destination.grid                         = source.grid;
  destination.lh_2d                        = source.lh_2d;
  destination.lookup                       = source.lookup;
  destination.ncid                         = source.ncid;
  destination.time_dimid                   = source.time_dimid;
  destination.bnds_dimid                   = source.bnds_dimid;
  destination.time_varid                   = source.time_varid;
  destination.time_bnds_varid              = source.time_bnds_varid;
  destination.col_dimid                    = source.col_dimid;
  destination.col_varid                    = source.col_varid;
  destination.zt_dimid                     = source.zt_dimid;
  destination.zm_dimid                     = source.zm_dimid;
  destination.rad_zt_dimid                 = source.rad_zt_dimid;
  destination.rad_zm_dimid                 = source.rad_zm_dimid;
  destination.sfc_dimid                    = source.sfc_dimid;
  destination.zt_varid                     = source.zt_varid;
  destination.zm_varid                     = source.zm_varid;
  destination.rad_zt_varid                 = source.rad_zt_varid;
  destination.rad_zm_varid                 = source.rad_zm_varid;
  destination.param_strlen_dimid           = source.param_strlen_dimid;
  destination.param_dimid                  = source.param_dimid;
  destination.param_varid                  = source.param_varid;
  destination.param_name_varid             = source.param_name_varid;
  destination.clubb_params_varid           = source.clubb_params_varid;
  destination.vars_allocated               = source.vars_allocated;
  destination.l_warned_overflow            = source.l_warned_overflow;
  destination.grid.gr_source.zm            = stats_value_copy(source.grid.gr_source.zm);
  destination.grid.gr_source.zt            = stats_value_copy(source.grid.gr_source.zt);
  destination.grid.gr_source.invrs_dzm     = stats_value_copy(source.grid.gr_source.invrs_dzm);
  destination.grid.gr_source.invrs_dzt     = stats_value_copy(source.grid.gr_source.invrs_dzt);
  destination.grid.gr_source.dzm           = stats_value_copy(source.grid.gr_source.dzm);
  destination.grid.gr_source.dzt           = stats_value_copy(source.grid.gr_source.dzt);
  destination.grid.gr_source.weights_zm2zt = stats_value_copy(source.grid.gr_source.weights_zm2zt);
  destination.grid.gr_source.weights_zt2zm = stats_value_copy(source.grid.gr_source.weights_zt2zm);
  destination.grid.gr_output.zm            = stats_value_copy(source.grid.gr_output.zm);
  destination.grid.gr_output.zt            = stats_value_copy(source.grid.gr_output.zt);
  destination.grid.gr_output.invrs_dzm     = stats_value_copy(source.grid.gr_output.invrs_dzm);
  destination.grid.gr_output.invrs_dzt     = stats_value_copy(source.grid.gr_output.invrs_dzt);
  destination.grid.gr_output.dzm           = stats_value_copy(source.grid.gr_output.dzm);
  destination.grid.gr_output.dzt           = stats_value_copy(source.grid.gr_output.dzt);
  destination.grid.gr_output.weights_zm2zt = stats_value_copy(source.grid.gr_output.weights_zm2zt);
  destination.grid.gr_output.weights_zt2zm = stats_value_copy(source.grid.gr_output.weights_zt2zm);
  destination.grid.rho_lin_spline_vals     = stats_value_copy(source.grid.rho_lin_spline_vals);
  destination.grid.rho_lin_spline_levels   = stats_value_copy(source.grid.rho_lin_spline_levels);
  destination.grid.p_sfc                   = stats_value_copy(source.grid.p_sfc);
  destination.lh_2d.nl_varids              = stats_value_copy(source.lh_2d.nl_varids);
  destination.lh_2d.u_varids               = stats_value_copy(source.lh_2d.u_varids);
  destination.lookup.cache                 = stats_value_copy(source.lookup.cache);
}

namespace {
// These call-local character/array adapters replace Fortran intrinsic and
// assumed-shape I/O conventions; they retain no registry or file state.
std::string stats_io_trim(const std::string& value) {
  // Fortran TRIM removes trailing blanks, retaining leading blanks and tabs.
  return value.substr(0, value.find_last_not_of(' ') + 1);
}
std::string stats_io_adjustl(const std::string& value) {
  // Fortran ADJUSTL moves leading blanks to the trailing character padding.
  const auto first = value.find_first_not_of(' ');
  return first == std::string::npos ? std::string{} : value.substr(first);
}
std::string stats_io_character(const std::string& value, std::size_t width) {
  // Fixed-width assignment truncates first; C++ scalar strings omit final padding.
  return stats_io_trim(value.substr(0, width));
}
template <class View>
int stats_io_put_var(ScratchWorkspace& stats_workspace, int ncid, int varid, const View& values,
                     const std::vector<int>& start, const std::vector<int>& count) {
  // The C API needs packed storage even for a padded/strided HostSpace actual.
  // Every caller supplies the Native start/count, avoiding file-shape guessing.
  if constexpr (std::is_same_v<typename View::memory_space, Kokkos::HostSpace>) {
    if (array_detail::dense_left_mapping(values))
      return netcdf::put_var(ncid, varid, values.data(), static_cast<int>(values.size()), start,
                             count);
  }
  auto stats_scope  = stats_workspace.scope();
  const auto packed = to_host_space(stats_workspace, values);
  return netcdf::put_var(ncid, varid, packed.data(), static_cast<int>(packed.size()), start, count);
}
} // namespace

//----------------------------------------------------------------------------
//                   Initialization/Finalize
//----------------------------------------------------------------------------
void stats_init_api(const std::string& registry_path, const std::string& output_path,
                    int ncol_batch, real_t stats_tsamp, real_t stats_tout, real_t dt_main,
                    int day_in, int month_in, int year_in, double time_initial, StatsHostReal1D zt,
                    StatsHostReal1D zm, stats_type& stats, err_info_type& err_info,
                    const StatsInitOptions& options) {
  // Initialize one output context: validate, expand registry, configure grids,
  // allocate runtime buffers, and define the optional NetCDF file in Native order.
  using namespace stats_netcdf_detail;
  stats.require_mutable();
  std::vector<stats_var> defs;
  std::string time_units;
  int ndefs = 0, grid_err = 0, ierr = 0, write_idx = 0;
  double effective_stats_tstart, effective_stats_tend;

  if (stats.enabled) {
    std::cerr << "stats error: stats_init_api called on enabled stats variable\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  ierr                          = 0;
  stats.l_different_output_grid = false;
  if (dt_main <= real_t(0)) {
    std::cerr << "stats init: dt_main must be positive\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  if (stats_tsamp <= real_t(0) || stats_tout <= real_t(0)) {
    std::cerr << "stats init: stats_tsamp and stats_tout must be positive\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  effective_stats_tstart = time_initial;
  if (options.stats_tstart)
    effective_stats_tstart = *options.stats_tstart;
  effective_stats_tend = effective_stats_tstart;
  if (options.stats_tend)
    effective_stats_tend = *options.stats_tend;
  if (effective_stats_tstart < time_initial) {
    std::cerr << "stats init: stats_tstart must be >= time_initial\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  if (options.stats_tend) {
    if (effective_stats_tend <= effective_stats_tstart) {
      std::cerr << "stats init: invalid stats_tstart/stats_tend window\n";
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
  }
  if (std::abs(static_cast<real_t>(effective_stats_tstart - time_initial) / dt_main -
               static_cast<real_t>(std::lround(
                   static_cast<real_t>(effective_stats_tstart - time_initial) / dt_main))) >
      real_t(1.e-8)) {
    std::cerr << "stats init: stats_tstart must align to a dt_main boundary\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  if (options.stats_tend) {
    if (std::abs(static_cast<real_t>(effective_stats_tend - time_initial) / dt_main -
                 static_cast<real_t>(std::lround(
                     static_cast<real_t>(effective_stats_tend - time_initial) / dt_main))) >
        real_t(1.e-8)) {
      std::cerr << "stats init: stats_tend must align to a dt_main boundary\n";
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
  }
  if (std::abs(stats_tsamp / dt_main - static_cast<real_t>(std::lround(stats_tsamp / dt_main))) >
      real_t(1.e-8)) {
    std::cerr << "stats init: stats_tsamp must be an integer multiple of dt_main\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  if (std::abs(stats_tout / dt_main - static_cast<real_t>(std::lround(stats_tout / dt_main))) >
      real_t(1.e-8)) {
    std::cerr << "stats init: stats_tout must be an integer multiple of dt_main\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  stats.stats_nsamp       = std::max(1, static_cast<int>(std::lround(stats_tsamp / dt_main)));
  stats.stats_nout        = std::max(1, static_cast<int>(std::lround(stats_tout / dt_main)));
  stats.samples_per_write = std::max(1, (stats.stats_nout / stats.stats_nsamp));
  stats.ncol_total        = ncol_batch;
  if (options.ncol_total)
    stats.ncol_total = *options.ncol_total;
  stats.active_batch_offset      = 0;
  stats.dt_main                  = static_cast<double>(dt_main);
  stats.time_initial             = time_initial;
  stats.tstart                   = effective_stats_tstart;
  stats.tend                     = effective_stats_tend;
  stats.l_open_ended_time_window = !options.stats_tend.has_value();
  if (stats.ncol_total < ncol_batch) {
    std::cerr << "stats init: ncol_total must be >= runtime ncol_batch\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  // We need stats_nout/stats_nsamp to be an integer.
  if ((stats.stats_nout % stats.stats_nsamp) != 0) {
    std::cerr << "stats init: stats_nout not a multiple of stats_nsamp\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  // Either both rad_zt and rad_zm can be defined, or neither.
  if (options.rad_zt.has_value() != options.rad_zm.has_value()) {
    std::cerr << "stats error: rad_zt and rad_zm must be used together, or not at all\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  } else {
    stats.l_output_rad_files = options.rad_zt.has_value() && options.rad_zm.has_value();
  }
  // Either both clubb_params and param_names can be defined, or neither.
  if (options.clubb_params.has_value() != options.param_names.has_value()) {
    std::cerr << "stats error: clubb_params and param_names must be used together\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  if (options.clubb_params) {
    if (options.clubb_params->extent_int(0) != stats.ncol_total) {
      std::cerr << "stats init: clubb_params metadata must have ncol_total rows\n";
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
  }
  // All three remapping options must be present together.
  if (options.grid_remap_method.has_value() != options.output_zt.has_value() ||
      options.grid_remap_method.has_value() != options.output_zm.has_value()) {
    std::cerr << "stats init: inconsistent initialization for grid remapping\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  // Figure out which variables are defined for output.
  stats_read_registry_namelist(stats.workspace(), registry_path, defs, ndefs, ierr);
  if (ierr != 0) {
    std::cerr << err_info.err_header_global << '\n' << "stats init: error reading registry\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  // Expand shorthand variables, e.g. sclrm -> sclr1m,sclr2m,... .
  stats_expand_registry(options.sclr_dim, options.edsclr_dim, defs, ndefs, options.hydromet_list,
                        stats.l_warned_overflow);
  // Without radiation grids, prune unsupported radiation entries in source order.
  if (!stats.l_output_rad_files) {
    write_idx = 0;
    for (int i = 0; i < ndefs; ++i) {
      if (defs[i].grid_id == GRID_RAD_ZT || defs[i].grid_id == GRID_RAD_ZM)
        continue;
      if (write_idx != i)
        defs[write_idx] = defs[i];
      ++write_idx;
    }
    ndefs = write_idx;
  }
  // Midnight on the case start date is the absolute model-clock NetCDF origin.
  format_date(day_in, month_in, year_in, 0., time_units);
  // Empty output path collects statistics without NetCDF I/O.
  stats.l_netcdf_output = output_path.find_last_not_of(' ') != std::string::npos;
  // Always retain the source grid, including memory-only consumers.
  stats.l_different_output_grid =
      stats.l_netcdf_output && options.grid_remap_method && options.output_zt && options.output_zm;
  if (stats.l_different_output_grid) {
    stats_grid_init(stats.workspace(), ncol_batch, zt, zm, *options.output_zt, *options.output_zm,
                    *options.grid_remap_method, stats.grid, grid_err);
  } else {
    stats_grid_init(stats.workspace(), ncol_batch, zt, zm, zt, zm, 0, stats.grid, grid_err);
  }
  if (grid_err != 0) {
    std::cerr << err_info.err_header_global << '\n'
              << "stats init: failed to initialize grid metadata\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  // Runtime buffers use the source grid; output sizes use the requested output grid.
  if (stats.l_different_output_grid && stats.l_netcdf_output) {
    stats_type_initialize(ncol_batch, defs, ndefs, *options.output_zt, *options.output_zm, stats,
                          ierr, options.rad_zt, options.rad_zm);
  } else {
    stats_type_initialize(ncol_batch, defs, ndefs, zt, zm, stats, ierr, options.rad_zt,
                          options.rad_zm);
  }
  if (ierr != 0) {
    std::cerr << err_info.err_header_global << '\n'
              << "stats init: failed to initialize runtime state\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  if (stats.l_netcdf_output) {
    // Native's OpenMP critical region is unnecessary for the supported
    // standalone path: host statistics calls execute sequentially. Concurrent
    // host callers would need one NetCDF I/O lock before being supported.
    if (stats.l_different_output_grid) {
      stats_open_netcdf(output_path, time_units, *options.output_zt, *options.output_zm, stats,
                        ierr, options.clubb_params, options.param_names, options.rad_zt,
                        options.rad_zm);
    } else {
      stats_open_netcdf(output_path, time_units, zt, zm, stats, ierr, options.clubb_params,
                        options.param_names, options.rad_zt, options.rad_zm);
    }
    if (ierr != 0) {
      std::cerr << err_info.err_header_global << '\n'
                << "stats init: failed to open netcdf output\n";
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
  }
  stats.enabled = true;
}

void stats_finalize_api(stats_type& stats, err_info_type& err_info) {
  // Release all partial/runtime allocations and close the owned file.
  stats.require_mutable();
  for (auto& var : stats.vars) {
    var.buffer                = {};
    var.nsamples              = {};
    var.sub_timestep_sum      = {};
    var.sub_timestep_nsamples = {};
  }
  stats.vars.clear();
  stats.vars_allocated = false;
  if (stats.lh_2d.is_initialized) {
    stats.lh_2d.nl_varids = {};
    stats.lh_2d.u_varids  = {};
    stats.lh_2d.nl_names.clear();
    stats.lh_2d.u_names.clear();
    stats.lh_2d.is_initialized = false;
  }
  if (stats.grid.is_initialized) {
    stats.grid.rho_lin_spline_vals   = {};
    stats.grid.rho_lin_spline_levels = {};
    stats.grid.p_sfc                 = {};
    stats.grid.gr_source.zm          = {};
    stats.grid.gr_source.zt          = {};
    stats.grid.gr_output.zm          = {};
    stats.grid.gr_output.zt          = {};
    stats.grid.is_initialized        = false;
  }
  stats.lookup.cache = {};
  stats.lookup.reject_cache.clear();
  stats.lookup.head              = 0;
  stats.lookup.cache_len         = 0;
  stats.lookup.reject_head       = 0;
  stats.lookup.reject_cache_len  = 0;
  stats.enabled                  = false;
  stats.l_sample                 = false;
  stats.l_last_sample            = false;
  stats.l_netcdf_output          = true;
  stats.l_open_ended_time_window = false;
  stats.nvars                    = 0;
  stats.ncol_batch               = 0;
  stats.ncol_total               = 0;
  stats.active_batch_offset      = 0;
  stats.time_index               = 0;
  // Every service frame has ended and its blocking transfers have completed.
  stats.workspace().clear();
  if (stats.ncid >= 0) {
    const int ret_code = netcdf::close(stats.ncid);
    if (ret_code != netcdf::nc_noerr) {
      std::cerr << err_info.err_header_global << '\n' << "stats finalize: netcdf close failed\n";
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
  }
  stats.ncid = -1;
}

void stats_reset_api(stats_type& stats) {
  // Clear per-run samples while preserving the registry/file and active batch rules.
  stats.require_mutable();
  if (!stats.enabled)
    return;
  if ((stats.active_batch_offset + stats.ncol_batch) >= stats.ncol_total)
    stats.active_batch_offset = 0;
  stats.l_sample                = false;
  stats.l_last_sample           = false;
  stats.time_index              = 0;
  stats.lookup.head             = 0;
  stats.lookup.cache_len        = 0;
  stats.lookup.reject_head      = 0;
  stats.lookup.reject_cache_len = 0;
  for (auto& var : stats.vars) {
    if (var.buffer.is_allocated())
      Kokkos::deep_copy(var.buffer, stat_t(0));
    if (var.nsamples.is_allocated())
      Kokkos::deep_copy(var.nsamples, 0);
    if (var.sub_timestep_sum.is_allocated())
      Kokkos::deep_copy(var.sub_timestep_sum, stat_t(0));
    if (var.sub_timestep_nsamples.is_allocated())
      Kokkos::deep_copy(var.sub_timestep_nsamples, 0);
    var.l_sub_timestep_active     = false;
    var.l_regular_timestep_active = false;
    var.l_budget                  = false;
    var.l_in_budget               = false;
  }
}

void start_next_stats_batch_api(stats_type& stats, err_info_type& err_info) {
  // Advance the output slice by one runtime batch without resizing buffers.
  stats.require_mutable();
  if (!stats.enabled)
    return;
  const int new_offset = (stats.active_batch_offset + stats.ncol_batch);
  if (stats.ncol_batch + new_offset > stats.ncol_total) {
    std::cerr << "stats start next batch: new batch exceeds total stats columns\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  stats.active_batch_offset = new_offset;
}

void stats_get_source_grid_api(const stats_type& stats, StatsHostReal2D& zt, StatsHostReal2D& zm) {
  // Return independent copies of current runtime/source heights, never output heights.
  if (!stats.grid.is_initialized) {
    zt = StatsHostReal2D("stats source zt", 0, 0);
    zm = StatsHostReal2D("stats source zm", 0, 0);
    return;
  }
  zt = to_host_space(stats.grid.gr_source.zt);
  zm = to_host_space(stats.grid.gr_source.zm);
}

namespace stats_netcdf_detail {
void stats_type_initialize(int ncol_batch, const std::vector<stats_var>& defs, int ndefs,
                           StatsHostReal1D zt, StatsHostReal1D zm, stats_type& stats, int& ierr,
                           std::optional<StatsHostReal1D> rad_zt,
                           std::optional<StatsHostReal1D> rad_zm) {
  // Build Native per-variable state, caches and source/output vertical sizes.
  ierr             = 0;
  stats.nvars      = ndefs;
  stats.ncol_batch = ncol_batch;
  stats.vars.clear();
  stats.vars.resize(stats.nvars);
  stats.vars_allocated = true;
  stats.lookup.cache   = StatsHostInt1D("stats lookup hit cache", std::max(16, stats.nvars));
  Kokkos::deep_copy(stats.lookup.cache, -1);
  stats.lookup.reject_cache.assign(std::max(16, stats.nvars), "");
  stats.lookup.head             = 0;
  stats.lookup.cache_len        = 0;
  stats.lookup.reject_head      = 0;
  stats.lookup.reject_cache_len = 0;
  const int nzt = zt.extent_int(0), nzm = zm.extent_int(0);
  int nrad_zt = 0, nrad_zm = 0;
  if (rad_zt && rad_zm) {
    nrad_zt = rad_zt->extent_int(0);
    nrad_zm = rad_zm->extent_int(0);
  }
  // Level counts depend on grid family and whether output remapping is active.
  for (int i = 0; i < ndefs; ++i) {
    auto& var = stats.vars[i];
    var       = defs[i];
    switch (var.grid_id) {
    case GRID_ZT:
      var.out_nz = nzt;
      var.nz     = stats.l_different_output_grid ? stats.grid.gr_source.nzt : var.out_nz;
      break;
    case GRID_ZM:
      var.out_nz = nzm;
      var.nz     = stats.l_different_output_grid ? stats.grid.gr_source.nzm : var.out_nz;
      break;
    case GRID_RAD_ZT:
      var.out_nz = nrad_zt;
      var.nz     = nrad_zt;
      break;
    case GRID_RAD_ZM:
      var.out_nz = nrad_zm;
      var.nz     = nrad_zm;
      break;
    case GRID_SFC:
    case GRID_LH_SFC:
      var.out_nz = 1;
      var.nz     = 1;
      break;
    case GRID_LH_ZT:
      // lh_zt uses source length even when adaptive remapping is active.
      if (stats.l_different_output_grid) {
        var.out_nz = stats.grid.gr_source.nzt;
        var.nz     = stats.grid.gr_source.nzt;
      } else {
        var.out_nz = nzt;
        var.nz     = nzt;
      }
      break;
    default:
      ierr = -1;
      std::cerr << "stats error in grid definition for: " << var.name << '\n';
      return;
    }
    var.buffer   = StatsReal2D("stats buffer " + var.name, stats.ncol_batch, var.nz);
    var.nsamples = StatsHostInt2D("stats sample counts " + var.name, stats.ncol_batch, var.nz);
    Kokkos::deep_copy(var.buffer, stat_t(0));
    Kokkos::deep_copy(var.nsamples, 0);
  }
}
} // namespace stats_netcdf_detail

namespace stats_netcdf_detail {
// -----------------------------------------------------------------------------
void stats_open_netcdf(const std::string& output_path, const std::string& time_units, // In
                       StatsHostReal1D zt, StatsHostReal1D zm,                        // In
                       stats_type& stats, int& ierr,                                  // InOut, Out
                       std::optional<StatsHostReal2D> clubb_params,                   // Optional In
                       std::optional<std::vector<std::string>> param_names,           // Optional In
                       std::optional<StatsHostReal1D> rad_zt,                         // Optional In
                       std::optional<StatsHostReal1D> rad_zm) {                       // Optional In
  // Description:
  //   Create/overwrite the NetCDF file, define the prepared runtime variables,
  //   and write coordinate and parameter metadata arrays in Native order.
  // output_path/time_units: File path and UDUNITS reference time [-].
  // zt/zm/rad_zt/rad_zm: Thermodynamic/momentum coordinate heights [m].
  // stats: Case-owned definitions and returned NetCDF identifiers [-].
  // ierr: NetCDF/library error code, zero on success [-].
  // clubb_params/param_names: Parameter matrix and fixed-width names [units vary].
  auto& stats_workspace = stats.workspace();
  auto stats_scope      = stats_workspace.scope();
  int ret_code, nzt, nzm, nrad_zt, nrad_zm, nparams, lh_nzt, nparam_strlen;

  // --------------------------- Begin Code ---------------------------
  ierr = 0;
  if (!stats.vars_allocated) {
    ierr = -1;
    std::cerr << "stats netcdf error: runtime state not prepared\n";
    return;
  }

  // Define size of zt (nzt) and zm (nzm).
  nzt = zt.extent_int(0);
  nzm = zm.extent_int(0);

  // Define size of rad_zt (nrad_zt) and rad_zm (nrad_zm).
  if (rad_zt && rad_zm) {
    nrad_zt = rad_zt->extent_int(0);
    nrad_zm = rad_zm->extent_int(0);
  } else {
    nrad_zt = 0;
    nrad_zm = 0;
  }

  // Define size of clubb_params (nparams).
  if (clubb_params && param_names)
    nparams = clubb_params->extent_int(1);
  else
    nparams = 0;

  // lh_zt variables keep source-grid heights when output remapping is active.
  if (stats.l_different_output_grid)
    lh_nzt = stats.grid.gr_source.nzt;
  else
    lh_nzt = nzt;

  // NF90_CLOBBER replaces any existing file. Preserve Native status overwrites
  // after the checked create call; later definition statuses are not promoted.
  ret_code = netcdf::create(stats_io_trim(output_path), netcdf::nc_clobber, stats.ncid);
  if (ret_code != netcdf::nc_noerr) {
    std::cerr << "stats netcdf error at nf90_create: " << netcdf::strerror(ret_code) << '\n';
    ierr = ret_code;
    return;
  }

  // Define time and explicit bounds metadata for averaged records.
  ret_code = netcdf::def_dim(stats.ncid, "time", netcdf::nc_unlimited, stats.time_dimid);
  ret_code = netcdf::def_dim(stats.ncid, "bnds", 2, stats.bnds_dimid);
  ret_code =
      netcdf::def_var(stats.ncid, "time", netcdf::nc_double, {stats.time_dimid}, stats.time_varid);
  ret_code = netcdf::put_att(stats.ncid, stats.time_varid, "units", stats_io_trim(time_units));
  ret_code = netcdf::put_att(stats.ncid, stats.time_varid, "bounds", "time_bnds");
  ret_code = netcdf::def_var(stats.ncid, "time_bnds", netcdf::nc_double,
                             {stats.bnds_dimid, stats.time_dimid}, stats.time_bnds_varid);
  ret_code = netcdf::put_att(stats.ncid, stats.time_bnds_varid, "units", stats_io_trim(time_units));
  ret_code =
      netcdf::put_att(stats.ncid, stats.time_bnds_varid, "interval_semantics", "(start, end]");
  ret_code = netcdf::put_att(stats.ncid, netcdf::nc_global, "model_time_initial",
                             real_t(stats.time_initial));
  ret_code = netcdf::put_att(stats.ncid, netcdf::nc_global, "stats_tstart", real_t(stats.tstart));
  if (stats.l_open_ended_time_window) {
    ret_code =
        netcdf::put_att(stats.ncid, netcdf::nc_global, "stats_time_window_mode", "open_ended");
  } else {
    ret_code = netcdf::put_att(stats.ncid, netcdf::nc_global, "stats_time_window_mode", "explicit");
    ret_code = netcdf::put_att(stats.ncid, netcdf::nc_global, "stats_tend", real_t(stats.tend));
  }
  ret_code = netcdf::put_att(stats.ncid, netcdf::nc_global, "stats_tsamp",
                             real_t(stats.stats_nsamp) * real_t(stats.dt_main));
  ret_code = netcdf::put_att(stats.ncid, netcdf::nc_global, "stats_tout",
                             real_t(stats.stats_nout) * real_t(stats.dt_main));

  // Define columns, surface and momentum/thermodynamic coordinates.
  ret_code = netcdf::def_dim(stats.ncid, "col", stats.ncol_total, stats.col_dimid);
  ret_code =
      netcdf::def_var(stats.ncid, "col", netcdf::nc_double, {stats.col_dimid}, stats.col_varid);
  ret_code = netcdf::put_att(stats.ncid, stats.col_varid, "units", "index");
  ret_code = netcdf::def_dim(stats.ncid, "sfc", 1, stats.sfc_dimid);
  ret_code = netcdf::def_dim(stats.ncid, "zm", nzm, stats.zm_dimid);
  ret_code = netcdf::def_var(stats.ncid, "zm", netcdf::nc_double, {stats.zm_dimid}, stats.zm_varid);
  ret_code = netcdf::put_att(stats.ncid, stats.zm_varid, "units", "m");
  ret_code = netcdf::def_dim(stats.ncid, "zt", nzt, stats.zt_dimid);
  ret_code = netcdf::def_var(stats.ncid, "zt", netcdf::nc_double, {stats.zt_dimid}, stats.zt_varid);
  ret_code = netcdf::put_att(stats.ncid, stats.zt_varid, "units", "m");

  bool has_lh_zt = false;
  for (const auto& var : stats.vars)
    has_lh_zt = has_lh_zt || var.grid_id == GRID_LH_ZT;
  if (has_lh_zt) {
    ret_code = netcdf::def_dim(stats.ncid, "lh_zt", lh_nzt, stats.lh_2d.lh_zt_dimid);
    ret_code = netcdf::def_var(stats.ncid, "lh_zt", netcdf::nc_double, {stats.lh_2d.lh_zt_dimid},
                               stats.lh_2d.lh_zt_varid);
    ret_code = netcdf::put_att(stats.ncid, stats.lh_2d.lh_zt_varid, "units", "m");
  }

  // Define radiation coordinates only when both radiation grids are present.
  if (nrad_zt > 0 && nrad_zm > 0) {
    ret_code = netcdf::def_dim(stats.ncid, "rad_zt", nrad_zt, stats.rad_zt_dimid);
    ret_code = netcdf::def_var(stats.ncid, "rad_zt", netcdf::nc_double, {stats.rad_zt_dimid},
                               stats.rad_zt_varid);
    ret_code = netcdf::put_att(stats.ncid, stats.rad_zt_varid, "units", "m");
    ret_code = netcdf::def_dim(stats.ncid, "rad_zm", nrad_zm, stats.rad_zm_dimid);
    ret_code = netcdf::def_var(stats.ncid, "rad_zm", netcdf::nc_double, {stats.rad_zm_dimid},
                               stats.rad_zm_varid);
    ret_code = netcdf::put_att(stats.ncid, stats.rad_zm_varid, "units", "m");
  }

  // Define parameter index, fixed-width names and parameter values.
  if (nparams > 0) {
    // Existing text_values/pack_strings retain one common width, including
    // trailing blanks and names beyond the selected nparams prefix.
    nparam_strlen = 1;
    for (const auto& name : *param_names)
      nparam_strlen = std::max(nparam_strlen, static_cast<int>(name.size()));
    ret_code = netcdf::def_dim(stats.ncid, "param_strlen", nparam_strlen, stats.param_strlen_dimid);
    ret_code = netcdf::def_dim(stats.ncid, "param", nparams, stats.param_dimid);
    ret_code = netcdf::def_var(stats.ncid, "param", netcdf::nc_double, {stats.param_dimid},
                               stats.param_varid);
    ret_code = netcdf::put_att(stats.ncid, stats.param_varid, "units", "index");
    ret_code =
        netcdf::def_var(stats.ncid, "param_name", netcdf::nc_char,
                        {stats.param_strlen_dimid, stats.param_dimid}, stats.param_name_varid);
    ret_code = netcdf::def_var(stats.ncid, "clubb_params", netcdf::nc_double,
                               {stats.col_dimid, stats.param_dimid}, stats.clubb_params_varid);
    ret_code = netcdf::put_att(stats.ncid, stats.clubb_params_varid, "long_name", "clubb_params");
  }

  // Define individual variables in registry order.
  for (int i = 0; i < stats.nvars; ++i) {
    auto& var = stats.vars[i];
    switch (var.grid_id) {
    case GRID_ZT:
      ret_code = netcdf::def_var(stats.ncid, stats_io_trim(var.name), netcdf::nc_double,
                                 {stats.col_dimid, stats.zt_dimid, stats.time_dimid}, var.varid);
      break;
    case GRID_ZM:
      ret_code = netcdf::def_var(stats.ncid, stats_io_trim(var.name), netcdf::nc_double,
                                 {stats.col_dimid, stats.zm_dimid, stats.time_dimid}, var.varid);
      break;
    case GRID_RAD_ZT:
      ret_code =
          netcdf::def_var(stats.ncid, stats_io_trim(var.name), netcdf::nc_double,
                          {stats.col_dimid, stats.rad_zt_dimid, stats.time_dimid}, var.varid);
      break;
    case GRID_RAD_ZM:
      ret_code =
          netcdf::def_var(stats.ncid, stats_io_trim(var.name), netcdf::nc_double,
                          {stats.col_dimid, stats.rad_zm_dimid, stats.time_dimid}, var.varid);
      break;
    case GRID_SFC:
    case GRID_LH_SFC:
      ret_code = netcdf::def_var(stats.ncid, stats_io_trim(var.name), netcdf::nc_double,
                                 {stats.col_dimid, stats.time_dimid}, var.varid);
      break;
    case GRID_LH_ZT:
      ret_code =
          netcdf::def_var(stats.ncid, stats_io_trim(var.name), netcdf::nc_double,
                          {stats.col_dimid, stats.lh_2d.lh_zt_dimid, stats.time_dimid}, var.varid);
      break;
    default:
      ierr = -1;
      std::cerr << "stats error in grid definition for: " << stats_io_trim(var.name) << '\n';
      return;
    }
    ret_code = netcdf::def_var_fill(stats.ncid, var.varid, 0, real_t(0));
    ret_code = netcdf::put_att(stats.ncid, var.varid, "long_name", stats_io_trim(var.long_name));
    ret_code = netcdf::put_att(stats.ncid, var.varid, "units", stats_io_trim(var.units));
  }

  // Leave define mode before writing coordinate and parameter values.
  ret_code                   = netcdf::enddef(stats.ncid);
  StatsHostReal1D col_values = stats_scope.view<StatsHostReal1D>(
      Kokkos::WithoutInitializing, "stats_open_netcdf::col", stats.ncol_total);
  for (int i = 0; i < stats.ncol_total; ++i)
    col_values(i) = real_t(((i + 1)));
  ret_code = stats_io_put_var(stats.workspace(), stats.ncid, stats.col_varid, col_values, {1},
                              {stats.ncol_total});
  ret_code = stats_io_put_var(stats.workspace(), stats.ncid, stats.zm_varid, zm, {1}, {nzm});
  ret_code = stats_io_put_var(stats.workspace(), stats.ncid, stats.zt_varid, zt, {1}, {nzt});

  // Native attempts this coordinate write even without a defined lh_zt ID;
  // its status is deliberately overwritten, like all surrounding writes.
  if (stats.l_different_output_grid) {
    const auto source_zt = Kokkos::subview(stats.grid.gr_source.zt, 0, Kokkos::ALL());
    ret_code = stats_io_put_var(stats.workspace(), stats.ncid, stats.lh_2d.lh_zt_varid, source_zt,
                                {1}, {stats.grid.gr_source.nzt});
  } else {
    ret_code =
        stats_io_put_var(stats.workspace(), stats.ncid, stats.lh_2d.lh_zt_varid, zt, {1}, {nzt});
  }
  if (nrad_zt > 0 && nrad_zm > 0) {
    ret_code = stats_io_put_var(stats.workspace(), stats.ncid, stats.rad_zt_varid, *rad_zt, {1},
                                {nrad_zt});
    ret_code = stats_io_put_var(stats.workspace(), stats.ncid, stats.rad_zm_varid, *rad_zm, {1},
                                {nrad_zm});
  }
  if (nparams > 0) {
    StatsHostReal1D indices = stats_scope.view<StatsHostReal1D>(
        Kokkos::WithoutInitializing, "stats_open_netcdf::param", nparams);
    for (int i = 0; i < nparams; ++i)
      indices(i) = real_t(((i + 1)));
    ret_code =
        stats_io_put_var(stats.workspace(), stats.ncid, stats.param_varid, indices, {1}, {nparams});
    auto names = stats_scope.view<Kokkos::View<char*, Kokkos::LayoutLeft, Kokkos::HostSpace>>(
        Kokkos::WithoutInitializing, "stats_open_netcdf::names",
        static_cast<std::size_t>(nparam_strlen) * nparams);
    Kokkos::deep_copy(names, ' ');
    for (int i = 0; i < nparams; ++i) {
      const auto& name = param_names->at(i); // Native requires this selected prefix to exist.
      std::copy(name.begin(), name.end(),
                names.data() + static_cast<std::size_t>(i) * nparam_strlen);
    }
    ret_code = netcdf::put_var(stats.ncid, stats.param_name_varid, names.data(),
                               static_cast<int>(names.size()), {1, 1}, {nparam_strlen, nparams});
    ret_code = stats_io_put_var(stats.workspace(), stats.ncid, stats.clubb_params_varid,
                                *clubb_params, {1, 1}, {clubb_params->extent_int(0), nparams});
  }
  (void)ret_code;
}
// -----------------------------------------------------------------------------
} // namespace stats_netcdf_detail

//----------------------------------------------------------------------------
//                       Begin/End timestep
//----------------------------------------------------------------------------
void stats_begin_timestep_api(int itime, stats_type& stats) {
  // Set sampling/write flags from the Native one-based model timestep and replay
  // expected lookup order from its first position on each eligible step.
  stats.require_mutable();
  if (!stats.enabled) {
    stats.l_sample      = false;
    stats.l_last_sample = false;
    return;
  }
  stats.l_sample           = false;
  stats.l_last_sample      = false;
  const double sample_time = stats.time_initial + static_cast<double>(itime) * stats.dt_main;
  if (sample_time <= stats.tstart ||
      (!stats.l_open_ended_time_window && sample_time > stats.tend)) {
    stats.lookup.head        = 0;
    stats.lookup.reject_head = 0;
    return;
  }
  const int step_offset      = static_cast<int>(std::lround(
      static_cast<real_t>(sample_time - stats.tstart) / static_cast<real_t>(stats.dt_main)));
  const double expected_time = stats.tstart + static_cast<double>(step_offset) * stats.dt_main;
  if (std::abs(sample_time - expected_time) > std::numeric_limits<double>::epsilon() * 10.) {
    stats.lookup.head        = 0;
    stats.lookup.reject_head = 0;
    return;
  }
  // Start each output window clean. Sampling uses the open/closed (tstart,tend] interval.
  if (stats.stats_nout > 0) {
    if (step_offset == 1 || ((step_offset - 1) % stats.stats_nout) == 0) {
      for (int i = 0; i < stats.nvars; ++i) {
        Kokkos::deep_copy(stats.vars[i].buffer, stat_t(0));
        Kokkos::deep_copy(stats.vars[i].nsamples, 0);
      }
    }
  }
  // A final partial output window is intentionally left unwritten.
  stats.l_sample      = (step_offset % stats.stats_nsamp) == 0;
  stats.l_last_sample = (step_offset % stats.stats_nout) == 0;
  for (int i = 0; i < stats.nvars; ++i) {
    stats.vars[i].l_sub_timestep_active     = false;
    stats.vars[i].l_regular_timestep_active = false;
  }
  stats.lookup.head        = 0;
  stats.lookup.reject_head = 0;
}

void stats_end_timestep_api(real_t time_value, stats_type& stats, err_info_type& err_info) {
  // Flush substep averages, check unfinished budgets, then write at output cadence.
  using namespace stats_netcdf_detail;
  stats.require_mutable();
  int ierr = 0;
  if (!stats.enabled)
    return;
  stats_flush_sub_timestep_averages(stats);
  if (!stats.l_last_sample)
    return;
  for (int i = 0; i < stats.nvars; ++i) {
    if (!stats.vars[i].l_in_budget)
      continue;
    std::cerr << "stats warning: budget left mid-update at end_timestep for " << stats.vars[i].name
              << '\n';
    stats.vars[i].l_in_budget = false;
  }
  stats_avg_and_write(time_value, stats, ierr);
  if (ierr != 0)
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
  stats.l_sample      = false;
  stats.l_last_sample = false;
}

namespace stats_netcdf_detail {
void stats_prepare_sub_timestep_average(stats_var& var) {
  // Allocate and zero pointwise within-timestep sums/counts on first use.
  if (!var.sub_timestep_sum.is_allocated()) {
    var.sub_timestep_sum =
        StatsReal2D("stats substep sum " + var.name, var.buffer.extent(0), var.buffer.extent(1));
    var.sub_timestep_nsamples = StatsHostInt2D("stats substep counts " + var.name,
                                               var.buffer.extent(0), var.buffer.extent(1));
    Kokkos::deep_copy(var.sub_timestep_sum, stat_t(0));
    Kokkos::deep_copy(var.sub_timestep_nsamples, 0);
  }
}

void stats_flush_sub_timestep_averages(stats_type& stats) {
  // Count each contributed within-timestep mean as one sample at that point.
  for (int i = 0; i < stats.nvars; ++i) {
    auto& var = stats.vars[i];
    if (!var.l_sub_timestep_active)
      continue;
    for (int k = 0; k < var.buffer.extent_int(1); ++k) {
      for (int c = 0; c < var.buffer.extent_int(0); ++c) {
        if (var.sub_timestep_nsamples(c, k) > 0) {
          var.buffer(c, k) =
              var.buffer(c, k) +
              var.sub_timestep_sum(c, k) / static_cast<stat_t>(var.sub_timestep_nsamples(c, k));
          var.nsamples(c, k) = var.nsamples(c, k) + 1;
        }
      }
    }
    Kokkos::deep_copy(var.sub_timestep_sum, stat_t(0));
    Kokkos::deep_copy(var.sub_timestep_nsamples, 0);
    var.l_sub_timestep_active = false;
    int max_count             = 0;
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

//----------------------------------------------------------------------------
//                         Sampling/Writing helpers
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//                         Sampling/Writing subroutines
//----------------------------------------------------------------------------
void stats_update(const std::string& name, real_t values, stats_type& stats,
                  std::optional<int> icol, std::optional<int> level, bool sub_timestep_average) {
  // Add one scalar at a surface-column or a profile-column/level point.
  using namespace stats_netcdf_detail;
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
  if (sub_timestep_average) {
    if (var.l_regular_timestep_active) {
      std::cerr << "stats error: mixed regular and sub-timestep updates for " << var.name << '\n';
      return;
    }
    stats_prepare_sub_timestep_average(var);
    var.l_sub_timestep_active = true;
    if (level) {
      // Native requires icol whenever level is present; indices are zero-based here.
      var.sub_timestep_sum(icol.value(), *level) =
          var.sub_timestep_sum(icol.value(), *level) + values;
      var.sub_timestep_nsamples(icol.value(), *level) =
          var.sub_timestep_nsamples(icol.value(), *level) + 1;
    } else if (icol) {
      var.sub_timestep_sum(*icol, 0)      = var.sub_timestep_sum(*icol, 0) + values;
      var.sub_timestep_nsamples(*icol, 0) = var.sub_timestep_nsamples(*icol, 0) + 1;
    }
    return;
  }
  if (var.l_sub_timestep_active) {
    std::cerr << "stats error: mixed sub-timestep and regular updates for " << var.name << '\n';
    return;
  }
  var.l_regular_timestep_active = true;
  if (level) {
    // Updating a specific level and column in a profile field.
    var.buffer(icol.value(), *level)   = var.buffer(icol.value(), *level) + values;
    var.nsamples(icol.value(), *level) = var.nsamples(icol.value(), *level) + 1;
  } else if (icol) {
    // Updating a specific column in a surface field.
    var.buffer(*icol, 0)   = var.buffer(*icol, 0) + values;
    var.nsamples(*icol, 0) = var.nsamples(*icol, 0) + 1;
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

void stats_begin_budget(const std::string& name, real_t values, stats_type& stats,
                        std::optional<int> icol) {
  // Start a scalar budget by subtracting the current surface-column value.
  using namespace stats_netcdf_detail;
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
  // Native's optional icol is required on this scalar branch.
  var.buffer(icol.value(), 0) = var.buffer(icol.value(), 0) - values;
  var.l_budget                = true;
  var.l_in_budget             = true;
}

void stats_update_budget(const std::string& name, real_t values, stats_type& stats,
                         std::optional<int> icol, std::optional<int> level) {
  // Add one intermediate contribution at a surface or profile point.
  using namespace stats_netcdf_detail;
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
  if (level) {
    var.buffer(icol.value(), *level) = var.buffer(icol.value(), *level) + values;
  } else if (icol) {
    var.buffer(*icol, 0) = var.buffer(*icol, 0) + values;
  }
}

void stats_finalize_budget(const std::string& name, real_t values, stats_type& stats,
                           std::optional<int> icol, bool l_count_sample) {
  // Add the scalar closing state and optionally increment its sample count.
  using namespace stats_netcdf_detail;
  stats.require_mutable();
  if (!stats.enabled)
    return;
  if (!stats.l_sample)
    return;
  const int id = stats_get_id(stats, name);
  if (id < 0)
    return;
  auto& var                  = stats.vars[id];
  const bool do_count_sample = l_count_sample;
  if (!var.l_in_budget) {
    std::cerr << "stats budget finalize without begin for " << var.name << '\n';
    return;
  }
  var.buffer(icol.value(), 0) = var.buffer(icol.value(), 0) + values;
  if (do_count_sample)
    var.nsamples(icol.value(), 0) = var.nsamples(icol.value(), 0) + 1;
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

namespace stats_netcdf_detail {
// -----------------------------------------------------------------------------
void stats_avg_and_write(real_t time_value, stats_type& stats, int& ierr) { // In, InOut, Out
  // Description:
  //   Average pointwise sampled buffers, optionally remap zt/zm profiles,
  //   and write one NetCDF record. The averaged buffers remain in the context.
  // time_value: Model time for this record [s].
  // stats: Registry, buffers, sample counters and active file slice [-].
  // ierr: Nonzero on the source's checked variable-write failure [-].
  auto& stats_workspace = stats.workspace();
  auto stats_scope      = stats_workspace.scope();
  int ret_code, grid_id, col_start, step_offset, window_start_step;
  int iv;
  StatsHostReal1D time_buf      = stats_scope.view<StatsHostReal1D>(Kokkos::WithoutInitializing,
                                                                    "stats_avg_and_write::time", 1);
  StatsHostReal2D time_bnds_buf = stats_scope.view<StatsHostReal2D>(
      Kokkos::WithoutInitializing, "stats_avg_and_write::bounds", 2, 1);
  double window_start_time; // Start of this output window [s].

  // --------------------------- Begin Code ---------------------------
  ierr = 0;
  if (!stats.enabled)
    return;
  col_start = (stats.active_batch_offset + 1);

  // Increment the record before writing its time coordinate.
  stats.time_index = (stats.time_index + 1);
  ret_code         = netcdf::nc_noerr;
  if (stats.l_netcdf_output) {
    time_buf(0) = time_value;
    ret_code    = stats_io_put_var(stats.workspace(), stats.ncid, stats.time_varid, time_buf,
                                   {stats.time_index}, {1});
    step_offset = static_cast<int>(
        std::round(real_t(double(time_value) - stats.tstart) / real_t(stats.dt_main)));
    window_start_step   = (stats.stats_nout * ((step_offset - 1) / stats.stats_nout));
    window_start_time   = stats.tstart + double(window_start_step) * stats.dt_main;
    time_bnds_buf(0, 0) = real_t(window_start_time);
    time_bnds_buf(1, 0) = time_value;
    ret_code = stats_io_put_var(stats.workspace(), stats.ncid, stats.time_bnds_varid, time_bnds_buf,
                                {1, stats.time_index}, {2, 1});
  }

  // Registry order determines both the write sequence and failure diagnostics.
  for (int i = 0; i < stats.nvars; ++i) {
    auto& var     = stats.vars[i];
    int max_count = std::numeric_limits<int>::lowest();
    for (int k = 0; k < var.nsamples.extent_int(1); ++k)
      for (int c = 0; c < var.nsamples.extent_int(0); ++c)
        max_count = std::max(max_count, var.nsamples(c, k));
    if (max_count < 1)
      continue;

    // Surface values have one less file dimension than profile values.
    grid_id = var.grid_id;
    if (grid_id == GRID_SFC || grid_id == GRID_LH_SFC) {
      for (int c = 0; c < var.nsamples.extent_int(0); ++c)
        if (var.nsamples(c, 0) > 0)
          var.buffer(c, 0) = var.buffer(c, 0) / stat_t(var.nsamples(c, 0));
      if (stats.l_netcdf_output) {
        const auto values = Kokkos::subview(var.buffer, Kokkos::ALL(), 0);
        ret_code          = stats_io_put_var(stats.workspace(), stats.ncid, var.varid, values,
                                             {col_start, stats.time_index}, {stats.ncol_batch, 1});
      }
    } else {
      for (int k = 0; k < var.nsamples.extent_int(1); ++k)
        for (int c = 0; c < var.nsamples.extent_int(0); ++c)
          if (var.nsamples(c, k) > 0)
            var.buffer(c, k) = var.buffer(c, k) / stat_t(var.nsamples(c, k));

      // Remap source-grid averages only for registered zt/zm profiles.
      if (stats.l_different_output_grid && (grid_id == GRID_ZM || grid_id == GRID_ZT)) {
        iv = 1;
        if (stats.l_netcdf_output) {
          // Native REAL(buffer,core_rknd) occurs before remapping; stat_t may
          // differ from core real_t. The existing API owns device/host staging.
          auto remap_scope              = stats_workspace.scope();
          StatsHostReal2D source_values = remap_scope.view<StatsHostReal2D>(
              Kokkos::WithoutInitializing, "stats_avg_and_write::core_buffer", var.buffer.extent(0),
              var.buffer.extent(1));
          for (int k = 0; k < source_values.extent_int(1); ++k)
            for (int c = 0; c < source_values.extent_int(0); ++c)
              source_values(c, k) = real_t(var.buffer(c, k));
          const auto remapped = remap_vals_to_target(
              stats_workspace, stats.ncol_batch, stats.grid.gr_source, stats.grid.gr_output, // In
              var.nz, source_values, var.out_nz,                                             // In
              stats.grid.rho_lin_spline_vals.extent_int(1),                                  // In
              stats.grid.rho_lin_spline_vals, stats.grid.rho_lin_spline_levels,              // In
              iv, stats.grid.p_sfc, stats.grid.grid_remap_method,                            // In
              grid_id == GRID_ZT);                                                           // In
          ret_code =
              stats_io_put_var(stats.workspace(), stats.ncid, var.varid, remapped,
                               {col_start, 1, stats.time_index}, {stats.ncol_batch, var.out_nz, 1});
        }
      } else if (stats.l_netcdf_output) {
        ret_code =
            stats_io_put_var(stats.workspace(), stats.ncid, var.varid, var.buffer,
                             {col_start, 1, stats.time_index}, {stats.ncol_batch, var.out_nz, 1});
      }
    }

    if (stats.l_netcdf_output && ret_code != netcdf::nc_noerr) {
      std::cerr << "stats write error for " << stats_io_trim(var.name)
                << " grid=" << stats_io_trim(var.grid) << " nz=" << var.nz
                << " out_nz=" << var.out_nz << " ncol_batch=" << stats.ncol_batch
                << " ncol_total=" << stats.ncol_total << " time_index=" << stats.time_index << '\n';
      std::cerr << "stats write error: " << netcdf::strerror(ret_code) << '\n';
      ierr = ret_code;
      return;
    }
  }
}
// -----------------------------------------------------------------------------

// ------------------------- Grid Remapping Helpers -------------------------
void stats_grid_init(ScratchWorkspace& stats_workspace, int ncol_batch, StatsHostReal1D zt_src,
                     StatsHostReal1D zm_src,                         // In
                     StatsHostReal1D zt_tgt, StatsHostReal1D zm_tgt, // In
                     int grid_remap_method, stats_grid_module& grid_ctx,
                     int& ierr) { // In, InOut, Out
  // Description:
  //   Allocate source/output altitude arrays for all columns and store the
  //   remapping selector. Native's body performs no additional size validation.
  // ncol_batch: Runtime column count [-]. zt/zm source/target: Heights [m].
  // grid_remap_method: Remapping selector [-]. grid_ctx: Caller-owned grids [-].
  // ierr: Initialization status, zero in this source body [-].
  auto stats_scope = stats_workspace.scope();
  int nzm_src, nzt_src, nzm_tgt, nzt_tgt;

  ierr                       = 0;
  nzm_src                    = zm_src.extent_int(0);
  nzt_src                    = zt_src.extent_int(0);
  nzm_tgt                    = zm_tgt.extent_int(0);
  nzt_tgt                    = zt_tgt.extent_int(0);
  grid_ctx.grid_remap_method = grid_remap_method;

  grid_ctx.gr_source.zm = {};
  grid_ctx.gr_source.zt = {};
  grid_ctx.gr_output.zm = {};
  grid_ctx.gr_output.zt = {};

  grid_ctx.gr_source.nzm    = nzm_src;
  grid_ctx.gr_source.nzt    = nzt_src;
  grid_ctx.gr_source.zm     = Real2D("stats_grid_init::source_zm", ncol_batch, nzm_src);
  grid_ctx.gr_source.zt     = Real2D("stats_grid_init::source_zt", ncol_batch, nzt_src);
  StatsHostReal2D source_zm = stats_scope.view<StatsHostReal2D>(
      Kokkos::WithoutInitializing, "stats_grid_init::source_zm_host", ncol_batch, nzm_src);
  StatsHostReal2D source_zt = stats_scope.view<StatsHostReal2D>(
      Kokkos::WithoutInitializing, "stats_grid_init::source_zt_host", ncol_batch, nzt_src);
  for (int k = 0; k < nzm_src; ++k)
    for (int c = 0; c < ncol_batch; ++c)
      source_zm(c, k) = zm_src(k);
  Kokkos::deep_copy(grid_ctx.gr_source.zm, source_zm);
  for (int k = 0; k < nzt_src; ++k)
    for (int c = 0; c < ncol_batch; ++c)
      source_zt(c, k) = zt_src(k);
  Kokkos::deep_copy(grid_ctx.gr_source.zt, source_zt);

  grid_ctx.gr_output.nzm    = nzm_tgt;
  grid_ctx.gr_output.nzt    = nzt_tgt;
  grid_ctx.gr_output.zm     = Real2D("stats_grid_init::target_zm", ncol_batch, nzm_tgt);
  grid_ctx.gr_output.zt     = Real2D("stats_grid_init::target_zt", ncol_batch, nzt_tgt);
  StatsHostReal2D target_zm = stats_scope.view<StatsHostReal2D>(
      Kokkos::WithoutInitializing, "stats_grid_init::target_zm_host", ncol_batch, nzm_tgt);
  StatsHostReal2D target_zt = stats_scope.view<StatsHostReal2D>(
      Kokkos::WithoutInitializing, "stats_grid_init::target_zt_host", ncol_batch, nzt_tgt);
  for (int k = 0; k < nzm_tgt; ++k)
    for (int c = 0; c < ncol_batch; ++c)
      target_zm(c, k) = zm_tgt(k);
  Kokkos::deep_copy(grid_ctx.gr_output.zm, target_zm);
  for (int k = 0; k < nzt_tgt; ++k)
    for (int c = 0; c < ncol_batch; ++c)
      target_zt(c, k) = zt_tgt(k);
  Kokkos::deep_copy(grid_ctx.gr_output.zt, target_zt);

  grid_ctx.is_initialized = true;
}
} // namespace stats_netcdf_detail

// -----------------------------------------------------------------------------
void stats_update_grid(StatsHostReal2D zt_src, StatsHostReal2D zm_src,       // In
                       StatsHostReal2D rho_vals, StatsHostReal2D rho_levels, // In
                       StatsHostReal1D p_sfc, stats_type& stats) {           // In, InOut
  // Description:
  //   Refresh the per-column source heights, density spline and surface
  //   pressure used to remap profile output at a write boundary.
  // zt_src/zm_src/rho_levels: Heights [m]. rho_vals: Density [kg m-3].
  // p_sfc: Surface pressure [Pa]. stats: Case-owned remapping state [-].
  stats.require_mutable();
  int nzt, nzm, nrho;
  if (!stats.enabled)
    return;
  if (!stats.l_different_output_grid)
    return;
  if (!stats.grid.is_initialized)
    return;
  if (zt_src.extent_int(0) != stats.ncol_batch)
    return;
  if (zm_src.extent_int(0) != stats.ncol_batch)
    return;
  if (rho_vals.extent_int(0) != stats.ncol_batch)
    return;
  if (rho_levels.extent_int(0) != stats.ncol_batch)
    return;
  if (p_sfc.extent_int(0) != stats.ncol_batch)
    return;

  nzt  = zt_src.extent_int(1);
  nzm  = zm_src.extent_int(1);
  nrho = rho_vals.extent_int(1);
  if (nzt != (nzm - 1))
    return;
  if (rho_levels.extent_int(1) != nrho)
    return;

  if (stats.grid.gr_source.nzt != nzt || stats.grid.gr_source.nzm != nzm) {
    stats.grid.gr_source.zt  = {};
    stats.grid.gr_source.zm  = {};
    stats.grid.gr_source.nzt = nzt;
    stats.grid.gr_source.nzm = nzm;
    stats.grid.gr_source.zt  = Real2D("stats_update_grid::zt", stats.ncol_batch, nzt);
    stats.grid.gr_source.zm  = Real2D("stats_update_grid::zm", stats.ncol_batch, nzm);
  }
  if (!stats.grid.rho_lin_spline_vals.is_allocated() ||
      stats.grid.rho_lin_spline_vals.extent_int(1) != nrho) {
    stats.grid.rho_lin_spline_vals   = {};
    stats.grid.rho_lin_spline_levels = {};
    stats.grid.rho_lin_spline_vals =
        StatsHostReal2D("stats_update_grid::rho", stats.ncol_batch, nrho);
    stats.grid.rho_lin_spline_levels =
        StatsHostReal2D("stats_update_grid::rho_levels", stats.ncol_batch, nrho);
  }
  if (!stats.grid.p_sfc.is_allocated())
    stats.grid.p_sfc = StatsHostReal1D("stats_update_grid::p_sfc", stats.ncol_batch);

  Kokkos::deep_copy(stats.grid.gr_source.zt, zt_src);
  Kokkos::deep_copy(stats.grid.gr_source.zm, zm_src);
  Kokkos::deep_copy(stats.grid.rho_lin_spline_vals, rho_vals);
  Kokkos::deep_copy(stats.grid.rho_lin_spline_levels, rho_levels);
  Kokkos::deep_copy(stats.grid.p_sfc, p_sfc);
}

// ------------------------- SILHS samples handlers -------------------------
void stats_lh_samples_init(int num_samples, int nzt,                     // In
                           const std::vector<std::string>& nl_var_names, // In
                           const std::vector<std::string>& u_var_names,  // In
                           StatsHostReal1D zt_vals,                      // In
                           stats_type& stats, err_info_type& err_info) { // InOut
  // Description:
  //   Define SILHS subcolumn dimensions/variables in the main file and write
  //   their source-grid coordinate. Empty variable-name lists remain present.
  // num_samples/nzt: Sample and level counts [-]. Names: Field identifiers [-].
  // zt_vals: SILHS heights [m]. stats/err_info: Case file and diagnostic state [-].
  using stats_netcdf_detail::NAME_LEN;
  [[maybe_unused]] int ret_code; // Native deliberately overwrites unchecked statuses.
  int nnl, nu;
  std::string var_name;
  stats.require_mutable();

  if (!stats.enabled)
    return;
  if (!stats.l_netcdf_output)
    return;
  if (num_samples < 1 || nzt < 1)
    return;
  if (stats.zt_dimid < 0)
    return;
  nnl = static_cast<int>(nl_var_names.size());
  nu  = static_cast<int>(u_var_names.size());
  if (stats.lh_2d.is_initialized) {
    if (stats.lh_2d.num_samples != num_samples || stats.lh_2d.nzt != nzt) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr << "stats lh init mismatch with existing definition\n";
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    }
    return;
  }

  // Native OMPcritical maps to the supported sequential host-call stream.
  // No concurrent-host I/O guarantee or global file registry is added.
  ret_code = netcdf::redef(stats.ncid);
  ret_code = netcdf::def_dim(stats.ncid, "lh_sample_number", num_samples, stats.lh_2d.sample_dimid);
  if (stats.lh_2d.lh_zt_dimid < 0)
    ret_code = netcdf::def_dim(stats.ncid, "lh_zt", nzt, stats.lh_2d.lh_zt_dimid);
  if (stats.lh_2d.lh_zt_varid < 0) {
    ret_code = netcdf::def_var(stats.ncid, "lh_zt", netcdf::nc_double, {stats.lh_2d.lh_zt_dimid},
                               stats.lh_2d.lh_zt_varid);
    ret_code = netcdf::put_att(stats.ncid, stats.lh_2d.lh_zt_varid, "units", "m");
  }
  stats.lh_2d.num_samples = num_samples;
  stats.lh_2d.nzt         = nzt;
  stats.lh_2d.n_nl_vars   = nnl;
  stats.lh_2d.n_u_vars    = nu;

  if (nnl > 0) {
    stats.lh_2d.nl_varids = StatsHostInt1D("stats_lh_samples_init::nl_varids", nnl);
    stats.lh_2d.nl_names.resize(nnl);
    for (int i = 0; i < nnl; ++i) {
      var_name =
          stats_io_character("lh_nl_" + stats_io_trim(stats_io_adjustl(nl_var_names[i])), NAME_LEN);
      stats.lh_2d.nl_names[i] = var_name;
      ret_code                = netcdf::def_var(
          stats.ncid, stats_io_trim(var_name), netcdf::nc_double,
          {stats.col_dimid, stats.lh_2d.sample_dimid, stats.lh_2d.lh_zt_dimid, stats.time_dimid},
          stats.lh_2d.nl_varids(i));
      ret_code = netcdf::def_var_fill(stats.ncid, stats.lh_2d.nl_varids(i), 0, real_t(0));
    }
  }
  if (nu > 0) {
    stats.lh_2d.u_varids = StatsHostInt1D("stats_lh_samples_init::u_varids", nu);
    stats.lh_2d.u_names.resize(nu);
    for (int i = 0; i < nu; ++i) {
      var_name =
          stats_io_character("lh_u_" + stats_io_trim(stats_io_adjustl(u_var_names[i])), NAME_LEN);
      stats.lh_2d.u_names[i] = var_name;
      ret_code               = netcdf::def_var(
          stats.ncid, stats_io_trim(var_name), netcdf::nc_double,
          {stats.col_dimid, stats.lh_2d.sample_dimid, stats.lh_2d.lh_zt_dimid, stats.time_dimid},
          stats.lh_2d.u_varids(i));
      ret_code = netcdf::def_var_fill(stats.ncid, stats.lh_2d.u_varids(i), 0, real_t(0));
    }
  }
  ret_code = netcdf::enddef(stats.ncid);
  ret_code =
      stats_io_put_var(stats.workspace(), stats.ncid, stats.lh_2d.lh_zt_varid, zt_vals, {1}, {nzt});
  stats.lh_2d.is_initialized = true;
  (void)ret_code;
}

// -----------------------------------------------------------------------------
void stats_lh_samples_write_lognormal(StatsHostReal4D samples,                      // In
                                      stats_type& stats, err_info_type& err_info) { // InOut
  // Description:
  //   Write (column,sample,zt,var) lognormal samples to each stored
  //   (column,sample,zt,time) variable at the next time record.
  // samples: Field-dependent units. stats/err_info: File and diagnostic state [-].
  int t, ret_code;
  stats.require_mutable();
  if (!stats.enabled)
    return;
  if (!stats.l_netcdf_output)
    return;
  if (!stats.lh_2d.is_initialized)
    return;
  if (!stats.l_last_sample)
    return;
  if (stats.lh_2d.n_nl_vars < 1)
    return;
  if (samples.extent_int(0) != stats.ncol_batch)
    return;
  if (samples.extent_int(1) != stats.lh_2d.num_samples)
    return;
  if (samples.extent_int(2) != stats.lh_2d.nzt)
    return;
  if (samples.extent_int(3) != stats.lh_2d.n_nl_vars)
    return;

  t = (stats.time_index + 1);
  for (int v = 0; v < stats.lh_2d.n_nl_vars; ++v) {
    const auto values = Kokkos::subview(samples, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(),
                                        std::make_pair(v, ((v + 1))));
    ret_code = stats_io_put_var(stats.workspace(), stats.ncid, stats.lh_2d.nl_varids(v), values,
                                {1, 1, 1, t},
                                {stats.ncol_batch, stats.lh_2d.num_samples, stats.lh_2d.nzt, 1});
    if (ret_code != netcdf::nc_noerr) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr << "stats lh lognormal write failed for " << stats_io_trim(stats.lh_2d.nl_names[v])
                << '\n';
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
  }
}

// -----------------------------------------------------------------------------
void stats_lh_samples_write_uniform(StatsHostReal4D uniform_vals,                 // In
                                    StatsHostInt3D mixture_comp,                  // In
                                    StatsHostReal3D sample_weights,               // In
                                    stats_type& stats, err_info_type& err_info) { // InOut
  // Description:
  //   Write uniform variates followed by mixture indices and sample weights
  //   in Native variable order, all with (column,sample,zt,time) file layout.
  // uniform_vals/sample_weights: Fractions [0–1]. mixture_comp: Mixture index [-].
  // stats/err_info: Case-owned file and diagnostic state [-].
  int t, ret_code, dp2;
  stats.require_mutable();
  auto& stats_workspace = stats.workspace();
  auto stats_scope      = stats_workspace.scope();
  if (!stats.enabled)
    return;
  if (!stats.l_netcdf_output)
    return;
  if (!stats.lh_2d.is_initialized)
    return;
  if (!stats.l_last_sample)
    return;
  if (stats.lh_2d.n_u_vars < 1)
    return;
  if (uniform_vals.extent_int(0) != stats.ncol_batch)
    return;
  if (uniform_vals.extent_int(1) != stats.lh_2d.num_samples)
    return;
  if (uniform_vals.extent_int(2) != stats.lh_2d.nzt)
    return;
  if (mixture_comp.extent_int(0) != stats.ncol_batch)
    return;
  if (mixture_comp.extent_int(1) != stats.lh_2d.num_samples)
    return;
  if (mixture_comp.extent_int(2) != stats.lh_2d.nzt)
    return;
  if (sample_weights.extent_int(0) != stats.ncol_batch)
    return;
  if (sample_weights.extent_int(1) != stats.lh_2d.num_samples)
    return;
  if (sample_weights.extent_int(2) != stats.lh_2d.nzt)
    return;
  dp2 = uniform_vals.extent_int(3);
  if (stats.lh_2d.n_u_vars != (dp2 + 2))
    return;

  t = (stats.time_index + 1);
  for (int v = 0; v < dp2; ++v) {
    const auto values = Kokkos::subview(uniform_vals, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(),
                                        std::make_pair(v, ((v + 1))));
    ret_code = stats_io_put_var(stats.workspace(), stats.ncid, stats.lh_2d.u_varids(v), values,
                                {1, 1, 1, t},
                                {stats.ncol_batch, stats.lh_2d.num_samples, stats.lh_2d.nzt, 1});
    if (ret_code != netcdf::nc_noerr) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr << "stats lh uniform write failed for " << stats_io_trim(stats.lh_2d.u_names[v])
                << '\n';
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
  }

  {
    auto mixture_scope = stats_workspace.scope();
    // Materialize Native REAL(mixture_comp,core_rknd) and SPREAD(dim4,ncopies1).
    StatsHostReal4D mixture_values = mixture_scope.view<StatsHostReal4D>(
        Kokkos::WithoutInitializing, "stats_lh_samples_write_uniform::mixture", stats.ncol_batch,
        stats.lh_2d.num_samples, stats.lh_2d.nzt, 1);
    for (int k = 0; k < stats.lh_2d.nzt; ++k)
      for (int s = 0; s < stats.lh_2d.num_samples; ++s)
        for (int c = 0; c < stats.ncol_batch; ++c)
          mixture_values(c, s, k, 0) = real_t(mixture_comp(c, s, k));
    ret_code = stats_io_put_var(stats.workspace(), stats.ncid, stats.lh_2d.u_varids(dp2),
                                mixture_values, {1, 1, 1, t},
                                {stats.ncol_batch, stats.lh_2d.num_samples, stats.lh_2d.nzt, 1});
    if (ret_code != netcdf::nc_noerr) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr << "stats lh uniform write failed for " << stats_io_trim(stats.lh_2d.u_names[dp2])
                << '\n';
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
  }
  {
    auto weight_scope             = stats_workspace.scope();
    StatsHostReal4D weight_values = weight_scope.view<StatsHostReal4D>(
        Kokkos::WithoutInitializing, "stats_lh_samples_write_uniform::weights", stats.ncol_batch,
        stats.lh_2d.num_samples, stats.lh_2d.nzt, 1);
    for (int k = 0; k < stats.lh_2d.nzt; ++k)
      for (int s = 0; s < stats.lh_2d.num_samples; ++s)
        for (int c = 0; c < stats.ncol_batch; ++c)
          weight_values(c, s, k, 0) = sample_weights(c, s, k);
    ret_code = stats_io_put_var(stats.workspace(), stats.ncid, stats.lh_2d.u_varids((dp2 + 1)),
                                weight_values, {1, 1, 1, t},
                                {stats.ncol_batch, stats.lh_2d.num_samples, stats.lh_2d.nzt, 1});
    if (ret_code != netcdf::nc_noerr) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr << "stats lh uniform write failed for "
                << stats_io_trim(stats.lh_2d.u_names[(dp2 + 1)]) << '\n';
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
  }
}

//----------------------------------------------------------------------------
//                         Helpers
//----------------------------------------------------------------------------
bool var_on_stats_list(const stats_type& stats, const std::string& name) {
  // Read-only registration query; unlike stats_get_id it does not truncate the needle.
  if (!stats.enabled || name.find_last_not_of(' ') == std::string::npos)
    return false;
  const auto first         = name.find_first_not_of(' ');
  const std::string needle = name.substr(first, name.find_last_not_of(' ') - first + 1);
  for (int i = 0; i < stats.nvars; ++i) {
    const std::string registered =
        stats.vars[i].name.substr(0, stats.vars[i].name.find_last_not_of(' ') + 1);
    if (registered == needle)
      return true;
  }
  return false;
}

namespace stats_netcdf_detail {
int stats_get_id(stats_type& stats, const std::string& name) {
  // Resolve a zero-based registry slot using ordered hit/reject caches.
  stats.require_mutable();
  int id = -1;
  if (!stats.enabled || name.find_last_not_of(' ') == std::string::npos)
    return id;
  const auto first = name.find_first_not_of(' ');
  std::string needle =
      name.substr(first, name.find_last_not_of(' ') - first + 1).substr(0, NAME_LEN);
  needle = needle.substr(0, needle.find_last_not_of(' ') + 1);
  // Fast path: expected lookup at the current zero-based cache head.
  if (stats.lookup.head < stats.lookup.cache_len) {
    const int cached_id = stats.lookup.cache(stats.lookup.head);
    const auto& value   = stats.vars[cached_id].name;
    if (value.substr(0, value.find_last_not_of(' ') + 1) == needle) {
      id = cached_id;
      ++stats.lookup.head;
      return id;
    }
  }
  // Fast reject path: expected missing name at current reject-cache head.
  if (stats.lookup.reject_head < stats.lookup.reject_cache_len) {
    const auto& value = stats.lookup.reject_cache[stats.lookup.reject_head];
    if (value.substr(0, value.find_last_not_of(' ') + 1) == needle) {
      ++stats.lookup.reject_head;
      return id;
    }
  }
  // Fallback reject path: a repeated miss not in the expected order.
  for (int i = 0; i < stats.lookup.reject_cache_len; ++i) {
    const auto& value = stats.lookup.reject_cache[i];
    if (value.substr(0, value.find_last_not_of(' ') + 1) == needle)
      return id;
  }
  // Fallback: linear search, preserving the source cache insertion/cursor rules.
  for (int i = 0; i < stats.nvars; ++i) {
    const auto& value = stats.vars[i].name;
    if (value.substr(0, value.find_last_not_of(' ') + 1) == needle) {
      id = i;
      if (stats.lookup.head >= stats.lookup.cache_len) {
        stats_lookup_append(stats, id);
        ++stats.lookup.head;
      } else if (stats.lookup.cache(stats.lookup.head) == id) {
        ++stats.lookup.head;
      }
      return id;
    }
  }
  stats_lookup_reject_append(stats, needle);
  return id;
}

void stats_lookup_append(stats_type& stats, int id) {
  // Grow and append one zero-based registry slot; unused cache entries are -1.
  if (!stats.lookup.cache.is_allocated()) {
    const int new_cap  = std::max(16, stats.nvars);
    stats.lookup.cache = StatsHostInt1D("stats lookup hit cache", new_cap);
    Kokkos::deep_copy(stats.lookup.cache, -1);
    stats.lookup.cache_len = 0;
  } else if (stats.lookup.cache_len >= stats.lookup.cache.extent_int(0)) {
    const int new_cap = std::max(16, (2 * stats.lookup.cache.extent_int(0)));
    StatsHostInt1D tmp("stats grown lookup hit cache", new_cap);
    Kokkos::deep_copy(tmp, -1);
    for (int i = 0; i < stats.lookup.cache_len; ++i)
      tmp(i) = stats.lookup.cache(i);
    stats.lookup.cache = std::move(tmp);
  }
  ++stats.lookup.cache_len;
  stats.lookup.cache((stats.lookup.cache_len - 1)) = id;
}

void stats_lookup_reject_append(stats_type& stats, const std::string& name) {
  // Grow and append a normalized missing name in source-observed order.
  if (name.find_last_not_of(' ') == std::string::npos)
    return;
  if (stats.lookup.reject_cache.empty()) {
    const int new_cap = std::max(16, stats.nvars);
    stats.lookup.reject_cache.assign(new_cap, "");
    stats.lookup.reject_cache_len = 0;
  } else if (stats.lookup.reject_cache_len >= static_cast<int>(stats.lookup.reject_cache.size())) {
    const int new_cap = std::max(16, 2 * static_cast<int>(stats.lookup.reject_cache.size()));
    std::vector<std::string> tmp(new_cap);
    for (int i = 0; i < stats.lookup.reject_cache_len; ++i)
      tmp[i] = stats.lookup.reject_cache[i];
    stats.lookup.reject_cache = std::move(tmp);
  }
  ++stats.lookup.reject_cache_len;
  const auto first = name.find_first_not_of(' ');
  stats.lookup.reject_cache[(stats.lookup.reject_cache_len - 1)] =
      name.substr(first, name.find_last_not_of(' ') - first + 1).substr(0, NAME_LEN);
}
} // namespace stats_netcdf_detail

namespace stats_netcdf_detail {
// -----------------------------------------------------------------------------
void format_date(int day_in, int month_in, int year_in, double time_in, // In
                 std::string& date) {                                   // Out
  // Description:
  //   Convert a model start date/time into the UDUNITS reference-time string.
  //   Adapted from Chris Golaz's original GrADS version. C++ fixed fields
  //   replace the Native internal CHARACTER writes without normalizing seconds.
  // day/month/year: Gregorian start date [-]. time_in/st_time: Elapsed time [s].
  // date: "seconds since YYYY-MM-DD HH:MM:SS.0", fixed35characters [-].
  int iday, imonth, iyear;
  double st_time;
  compute_current_date_api(day_in, month_in, year_in, time_in, // In
                           iday, imonth, iyear, st_time);      // Out

  // Native Iw.w pads with zeros and prints stars when the sign/digits overflow.
  const auto integer_field = [](int value, std::size_t width) {
    // Format the fixed-width integer substring used by Native internal writes.
    const auto text = std::to_string(value);
    return value < 0 || text.size() > width ? std::string(width, '*')
                                            : std::string(width - text.size(), '0') + text;
  };
  date = "seconds since YYYY-MM-DD HH:MM:00.0";
  date.replace(14, 4, integer_field(iyear, 4));
  date.replace(19, 2, integer_field(imonth, 2));
  date.replace(22, 2, integer_field(iday, 2));
  date.replace(25, 2, integer_field(static_cast<int>(std::floor(st_time / 3600.0)), 2));
  date.replace(
      28, 2,
      integer_field(static_cast<int>((static_cast<int>(std::round(st_time)) % 3600) / 60), 2));
  date.replace(
      31, 2,
      integer_field(
          static_cast<int>(std::round(
              ((double(static_cast<int>(std::round(st_time)) % 3600) / 60.0) -
               double(static_cast<int>((static_cast<int>(std::round(st_time)) % 3600) / 60))) *
              60.0)),
          2));
}

// ----------------------------- Namelist parsing -----------------------------
void stats_read_registry_namelist(ScratchWorkspace& stats_workspace,                     // InOut
                                  const std::string& path,                               // In
                                  std::vector<stats_var>& defs, int& ndefs, int& ierr) { // Out
  // Description:
  //   Read fixed-width entry(i) strings from &clubb_stats_nl, filter empty or
  //   comment entries, and parse the active definitions in their original order.
  // path: Registry file [-]. defs/ndefs: Parsed definitions/count [-].
  // ierr: Portable read/allocation status; invalid registry line returns -1 [-].
  auto stats_scope = stats_workspace.scope();
  using Entry      = Kokkos::View<char**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  Entry entry;
  int ios, count;
  bool l_parse_ok;
  defs.clear(); // Native allocatable INTENT(OUT) deallocates the prior array.
  ierr  = 0;
  ndefs = 0;
  try {
    entry =
        stats_scope.view<Entry>(Kokkos::WithoutInitializing, "stats_read_registry_namelist::entry",
                                REG_LINE_LEN, NML_REG_MAX_ENTRIES);
  } catch (const std::exception&) {
    ierr = 1; // Native ALLOCATE STAT failure; vendor-specific code is not reproduced.
    return;
  }
  Kokkos::deep_copy(entry, ' ');
  namelist::Group group("clubb_stats_nl");
  group.bind("entry", entry);
  ios = group.read_file(path);
  if (ios != 0) {
    ierr = ios;
    return;
  }

  count = 0;
  for (int i = 0; i < NML_REG_MAX_ENTRIES; ++i) {
    const std::string line(entry.data() + static_cast<std::size_t>(i) * REG_LINE_LEN, REG_LINE_LEN);
    if (stats_io_trim(line).empty())
      continue;
    if (is_line_comment(line))
      continue;
    count = (count + 1);
  }
  if (count < 1)
    return;

  defs.resize(count);
  ndefs = 0;
  for (int i = 0; i < NML_REG_MAX_ENTRIES; ++i) {
    const std::string line(entry.data() + static_cast<std::size_t>(i) * REG_LINE_LEN, REG_LINE_LEN);
    if (stats_io_trim(line).empty())
      continue;
    if (is_line_comment(line))
      continue;
    ndefs = (ndefs + 1);
    parse_registry_line(line, defs[(ndefs - 1)], l_parse_ok); // In, Out
    if (!l_parse_ok) {
      ierr = -1;
      std::cerr << "Invalid stats registry entry (expected: name | grid | units | long_name)\n";
      std::cerr << stats_io_trim(line) << '\n';
      return;
    }
  }
}

// -----------------------------------------------------------------------------
void stats_expand_registry(
    std::optional<int> sclr_dim, std::optional<int> edsclr_dim,   // In
    std::vector<stats_var>& defs, int& ndefs,                     // InOut
    const std::optional<std::vector<std::string>>& hydromet_list, // Optional In
    bool& l_warned_overflow) {                                    // InOut
  // Description:
  //   Expand registry shorthand into explicit hydrometeor/scalar definitions,
  //   retaining first-name wins, source iteration order and the fixed capacity.
  // sclr_dim/edsclr_dim: Optional forwarded scalar counts [-]. defs/ndefs: Registry [-].
  // hydromet_list: Names in their Native fixed-width character order [-].
  // l_warned_overflow: Case-owned warning metadata, retained for context lifetime.
  // Native SAVE warns once per process; this explicitly approved adaptation warns
  // once per stats context, without changing truncation, deduplication or output.
  std::vector<stats_var> out_defs;
  stats_var base_def;
  std::string idx_str;
  int max_defs, k, nhm;

  nhm = 0;
  if (hydromet_list)
    nhm = static_cast<int>(hydromet_list->size());
  max_defs = std::max(1, (ndefs + 4096));
  out_defs.resize(max_defs);
  k = 0;
  for (int i = 0; i < ndefs; ++i) {
    base_def = defs[i];
    if (stats_io_trim(base_def.name) == "hm_i") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def, stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_1",
                         out_defs, k, l_warned_overflow);
        add_expanded_def(base_def, stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_2",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "mu_hm_i") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def, "mu_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_1",
                         out_defs, k, l_warned_overflow);
        add_expanded_def(base_def, "mu_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_2",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "mu_Ncn_i") {
      add_expanded_def(base_def, "mu_Ncn_1", out_defs, k, l_warned_overflow);
      add_expanded_def(base_def, "mu_Ncn_2", out_defs, k, l_warned_overflow);
    } else if (stats_io_trim(base_def.name) == "mu_hm_i_n") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def,
                         "mu_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_1_n",
                         out_defs, k, l_warned_overflow);
        add_expanded_def(base_def,
                         "mu_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_2_n",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "mu_Ncn_i_n") {
      add_expanded_def(base_def, "mu_Ncn_1_n", out_defs, k, l_warned_overflow);
      add_expanded_def(base_def, "mu_Ncn_2_n", out_defs, k, l_warned_overflow);
    } else if (stats_io_trim(base_def.name) == "sigma_hm_i") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def,
                         "sigma_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_1",
                         out_defs, k, l_warned_overflow);
        add_expanded_def(base_def,
                         "sigma_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_2",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "sigma_Ncn_i") {
      add_expanded_def(base_def, "sigma_Ncn_1", out_defs, k, l_warned_overflow);
      add_expanded_def(base_def, "sigma_Ncn_2", out_defs, k, l_warned_overflow);
    } else if (stats_io_trim(base_def.name) == "sigma_hm_i_n") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def,
                         "sigma_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_1_n",
                         out_defs, k, l_warned_overflow);
        add_expanded_def(base_def,
                         "sigma_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_2_n",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "sigma_Ncn_i_n") {
      add_expanded_def(base_def, "sigma_Ncn_1_n", out_defs, k, l_warned_overflow);
      add_expanded_def(base_def, "sigma_Ncn_2_n", out_defs, k, l_warned_overflow);
    } else if (stats_io_trim(base_def.name) == "corr_w_hm_i") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def,
                         "corr_w_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_1",
                         out_defs, k, l_warned_overflow);
        add_expanded_def(base_def,
                         "corr_w_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_2",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "corr_w_Ncn_i") {
      add_expanded_def(base_def, "corr_w_Ncn_1", out_defs, k, l_warned_overflow);
      add_expanded_def(base_def, "corr_w_Ncn_2", out_defs, k, l_warned_overflow);
    } else if (stats_io_trim(base_def.name) == "corr_chi_hm_i") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def,
                         "corr_chi_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_1",
                         out_defs, k, l_warned_overflow);
        add_expanded_def(base_def,
                         "corr_chi_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_2",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "corr_chi_Ncn_i") {
      add_expanded_def(base_def, "corr_chi_Ncn_1", out_defs, k, l_warned_overflow);
      add_expanded_def(base_def, "corr_chi_Ncn_2", out_defs, k, l_warned_overflow);
    } else if (stats_io_trim(base_def.name) == "corr_eta_hm_i") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def,
                         "corr_eta_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_1",
                         out_defs, k, l_warned_overflow);
        add_expanded_def(base_def,
                         "corr_eta_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_2",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "corr_eta_Ncn_i") {
      add_expanded_def(base_def, "corr_eta_Ncn_1", out_defs, k, l_warned_overflow);
      add_expanded_def(base_def, "corr_eta_Ncn_2", out_defs, k, l_warned_overflow);
    } else if (stats_io_trim(base_def.name) == "corr_Ncn_hm_i") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def,
                         "corr_Ncn_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_1",
                         out_defs, k, l_warned_overflow);
        add_expanded_def(base_def,
                         "corr_Ncn_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_2",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "corr_hmx_hmy_i") {
      for (int j = 0; j < nhm; ++j) {
        for (int j2 = j + 1; j2 < nhm; ++j2) {
          add_expanded_def(
              base_def,
              "corr_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_" +
                  stats_io_trim(hydromet_list->at(j2).substr(
                      0, std::min(std::size_t(2), stats_io_trim(hydromet_list->at(j2)).size()))) +
                  "_1",
              out_defs, k, l_warned_overflow);
          add_expanded_def(
              base_def,
              "corr_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_" +
                  stats_io_trim(hydromet_list->at(j2).substr(
                      0, std::min(std::size_t(2), stats_io_trim(hydromet_list->at(j2)).size()))) +
                  "_2",
              out_defs, k, l_warned_overflow);
        }
      }
    } else if (stats_io_trim(base_def.name) == "corr_w_hm_i_n") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def,
                         "corr_w_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_1_n",
                         out_defs, k, l_warned_overflow);
        add_expanded_def(base_def,
                         "corr_w_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_2_n",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "corr_w_Ncn_i_n") {
      add_expanded_def(base_def, "corr_w_Ncn_1_n", out_defs, k, l_warned_overflow);
      add_expanded_def(base_def, "corr_w_Ncn_2_n", out_defs, k, l_warned_overflow);
    } else if (stats_io_trim(base_def.name) == "corr_chi_hm_i_n") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def,
                         "corr_chi_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_1_n",
                         out_defs, k, l_warned_overflow);
        add_expanded_def(base_def,
                         "corr_chi_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_2_n",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "corr_chi_Ncn_i_n") {
      add_expanded_def(base_def, "corr_chi_Ncn_1_n", out_defs, k, l_warned_overflow);
      add_expanded_def(base_def, "corr_chi_Ncn_2_n", out_defs, k, l_warned_overflow);
    } else if (stats_io_trim(base_def.name) == "corr_eta_hm_i_n") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def,
                         "corr_eta_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_1_n",
                         out_defs, k, l_warned_overflow);
        add_expanded_def(base_def,
                         "corr_eta_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_2_n",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "corr_eta_Ncn_i_n") {
      add_expanded_def(base_def, "corr_eta_Ncn_1_n", out_defs, k, l_warned_overflow);
      add_expanded_def(base_def, "corr_eta_Ncn_2_n", out_defs, k, l_warned_overflow);
    } else if (stats_io_trim(base_def.name) == "corr_Ncn_hm_i_n") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def,
                         "corr_Ncn_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_1_n",
                         out_defs, k, l_warned_overflow);
        add_expanded_def(base_def,
                         "corr_Ncn_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_2_n",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "corr_hmx_hmy_i_n") {
      for (int j = 0; j < nhm; ++j) {
        for (int j2 = j + 1; j2 < nhm; ++j2) {
          add_expanded_def(
              base_def,
              "corr_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_" +
                  stats_io_trim(hydromet_list->at(j2).substr(
                      0, std::min(std::size_t(2), stats_io_trim(hydromet_list->at(j2)).size()))) +
                  "_1_n",
              out_defs, k, l_warned_overflow);
          add_expanded_def(
              base_def,
              "corr_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "_" +
                  stats_io_trim(hydromet_list->at(j2).substr(
                      0, std::min(std::size_t(2), stats_io_trim(hydromet_list->at(j2)).size()))) +
                  "_2_n",
              out_defs, k, l_warned_overflow);
        }
      }
    } else if (stats_io_trim(base_def.name) == "wp2hmp") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def, "wp2" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "p",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "K_hm") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def, "K_hm_" + stats_io_trim(hydromet_list->at(j).substr(0, 2)),
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "hydrometp2") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def, stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "p2",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "wphydrometp") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def, "wp" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "p",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "hmp2_zt") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def, stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "p2_zt",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "rtphmp") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def, "rtp" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "p",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "thlphmp") {
      for (int j = 0; j < nhm; ++j) {
        add_expanded_def(base_def, "thlp" + stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "p",
                         out_defs, k, l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "hmxphmyp") {
      for (int j = 0; j < nhm; ++j) {
        for (int j2 = j + 1; j2 < nhm; ++j2) {
          add_expanded_def(
              base_def,
              stats_io_trim(hydromet_list->at(j).substr(0, 2)) + "p" +
                  stats_io_trim(hydromet_list->at(j2).substr(
                      0, std::min(std::size_t(2), stats_io_trim(hydromet_list->at(j2)).size()))) +
                  "p",
              out_defs, k, l_warned_overflow);
        }
      }
    } else if (stats_io_trim(base_def.name) == "sclrm") {
      for (int j = 0; j < sclr_dim.value(); ++j) {
        idx_str = std::to_string(((j + 1)));
        add_expanded_def(base_def, "sclr" + stats_io_trim(idx_str) + "m", out_defs, k,
                         l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "sclrm_f") {
      for (int j = 0; j < sclr_dim.value(); ++j) {
        idx_str = std::to_string(((j + 1)));
        add_expanded_def(base_def, "sclr" + stats_io_trim(idx_str) + "m_f", out_defs, k,
                         l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "edsclrm") {
      for (int j = 0; j < edsclr_dim.value(); ++j) {
        idx_str = std::to_string(((j + 1)));
        add_expanded_def(base_def, "edsclr" + stats_io_trim(idx_str) + "m", out_defs, k,
                         l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "edsclrm_f") {
      for (int j = 0; j < edsclr_dim.value(); ++j) {
        idx_str = std::to_string(((j + 1)));
        add_expanded_def(base_def, "edsclr" + stats_io_trim(idx_str) + "m_f", out_defs, k,
                         l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "sclrprtp") {
      for (int j = 0; j < sclr_dim.value(); ++j) {
        idx_str = std::to_string(((j + 1)));
        add_expanded_def(base_def, "sclr" + stats_io_trim(idx_str) + "prtp", out_defs, k,
                         l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "sclrp2") {
      for (int j = 0; j < sclr_dim.value(); ++j) {
        idx_str = std::to_string(((j + 1)));
        add_expanded_def(base_def, "sclr" + stats_io_trim(idx_str) + "p2", out_defs, k,
                         l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "sclrpthvp") {
      for (int j = 0; j < sclr_dim.value(); ++j) {
        idx_str = std::to_string(((j + 1)));
        add_expanded_def(base_def, "sclr" + stats_io_trim(idx_str) + "pthvp", out_defs, k,
                         l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "sclrpthlp") {
      for (int j = 0; j < sclr_dim.value(); ++j) {
        idx_str = std::to_string(((j + 1)));
        add_expanded_def(base_def, "sclr" + stats_io_trim(idx_str) + "pthlp", out_defs, k,
                         l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "sclrprcp") {
      for (int j = 0; j < sclr_dim.value(); ++j) {
        idx_str = std::to_string(((j + 1)));
        add_expanded_def(base_def, "sclr" + stats_io_trim(idx_str) + "prcp", out_defs, k,
                         l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "wpsclrp") {
      for (int j = 0; j < sclr_dim.value(); ++j) {
        idx_str = std::to_string(((j + 1)));
        add_expanded_def(base_def, "wpsclr" + stats_io_trim(idx_str) + "p", out_defs, k,
                         l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "wpsclrp2") {
      for (int j = 0; j < sclr_dim.value(); ++j) {
        idx_str = std::to_string(((j + 1)));
        add_expanded_def(base_def, "wpsclr" + stats_io_trim(idx_str) + "p2", out_defs, k,
                         l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "wp2sclrp") {
      for (int j = 0; j < sclr_dim.value(); ++j) {
        idx_str = std::to_string(((j + 1)));
        add_expanded_def(base_def, "wp2sclr" + stats_io_trim(idx_str) + "p", out_defs, k,
                         l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "wpsclrprtp") {
      for (int j = 0; j < sclr_dim.value(); ++j) {
        idx_str = std::to_string(((j + 1)));
        add_expanded_def(base_def, "wpsclr" + stats_io_trim(idx_str) + "prtp", out_defs, k,
                         l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "wpsclrpthlp") {
      for (int j = 0; j < sclr_dim.value(); ++j) {
        idx_str = std::to_string(((j + 1)));
        add_expanded_def(base_def, "wpsclr" + stats_io_trim(idx_str) + "pthlp", out_defs, k,
                         l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "silhs_variance_category") {
      for (int j = 0; j < 8; ++j) {
        idx_str = std::to_string(((j + 1)));
        add_expanded_def(base_def, "silhs_var_cat_" + stats_io_trim(idx_str), out_defs, k,
                         l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "lh_samp_frac_category") {
      for (int j = 0; j < 8; ++j) {
        idx_str = std::to_string(((j + 1)));
        add_expanded_def(base_def, "lh_samp_frac_" + stats_io_trim(idx_str), out_defs, k,
                         l_warned_overflow);
      }
    } else if (stats_io_trim(base_def.name) == "wpedsclrp") {
      for (int j = 0; j < edsclr_dim.value(); ++j) {
        idx_str = std::to_string(((j + 1)));
        add_expanded_def(base_def, "wpedsclr" + stats_io_trim(idx_str) + "p", out_defs, k,
                         l_warned_overflow);
      }
    } else {
      add_expanded_def(base_def, base_def.name, out_defs, k, l_warned_overflow);
    }
  }

  // Native DEALLOCATE requires allocated defs. An empty/comment-only input
  // returns an unallocated array from the reader; retain that fatal boundary.
  if (defs.empty())
    throw std::runtime_error("stats registry expansion: Native definitions are unallocated");
  defs.clear();
  defs.resize(std::max(1, k));
  for (int i = 0; i < k; ++i)
    defs[i] = out_defs[i];
  ndefs = k;
  out_defs.clear();
}

// -----------------------------------------------------------------------------
void add_expanded_def(const stats_var& base_def, const std::string& name, // In
                      std::vector<stats_var>& out_defs, int& nout,        // InOut
                      bool& l_warned_overflow) {                          // InOut
  // Description:
  //   Append a definition unless the fixed-width name is empty/already present
  //   or the preallocated expansion array is full.
  // base_def/name: Template and expanded variable name [-].
  // out_defs/nout: Destination slots and current active count [-].
  // l_warned_overflow: Caller-owned once-per-context warning flag. Native's
  //   SAVE flag is once-per-process; registry truncation and order are unchanged.
  stats_var def;
  std::string needle;
  needle = stats_io_character(stats_io_trim(name), NAME_LEN);
  if (stats_io_trim(needle).empty())
    return;
  for (int i = 0; i < nout; ++i)
    if (stats_io_trim(out_defs[i].name) == needle)
      return;
  if (nout >= static_cast<int>(out_defs.size())) {
    if (!l_warned_overflow) {
      std::cerr << "stats registry expansion overflow; some template vars may be skipped\n";
      l_warned_overflow = true;
    }
    return;
  }

  def                  = base_def;
  def.name             = needle;
  def.varid            = -1;
  def.nz               = 0;
  nout                 = (nout + 1);
  out_defs[(nout - 1)] = def;
}

// -----------------------------------------------------------------------------
void parse_registry_line(const std::string& line,         // In
                         stats_var& def, bool& l_valid) { // Out
  // Description:
  //   Parse name | grid | units | long_name into a variable definition and
  //   map the grid token. Legacy three-field entries omit the units column.
  // line: Raw registry line [-]. def: Metadata with Native width truncation [-].
  // l_valid: Whether the expected fields and nonempty name/grid were found [-].
  std::vector<std::string> fields(4);
  int nfields;

  // INTENT(OUT) deallocates components and initializes defaulted Native fields.
  def             = stats_var{};
  def.name        = "";
  def.grid        = "";
  def.grid_id     = GRID_ZT;
  def.units       = "";
  def.long_name   = "";
  def.varid       = -1;
  def.nz          = 0;
  def.out_nz      = 0;
  def.l_budget    = false;
  def.l_in_budget = false;
  fields.assign(4, "");
  l_valid = false;
  split_registry_fields(line, fields, nfields); // In, InOut, Out
  if (nfields != 3 && nfields != 4)
    return;
  def.name = stats_io_character(stats_io_trim(stats_io_adjustl(fields[0])), NAME_LEN);
  def.grid = stats_io_character(stats_io_trim(stats_io_adjustl(fields[1])), GRID_LEN);
  if (nfields == 4) {
    def.units     = stats_io_character(stats_io_trim(stats_io_adjustl(fields[2])), UNITS_LEN);
    def.long_name = stats_io_character(stats_io_trim(stats_io_adjustl(fields[3])), DESC_LEN);
  } else {
    // Legacy three-field compatibility: name | grid | long_name.
    def.long_name = stats_io_character(stats_io_trim(stats_io_adjustl(fields[2])), DESC_LEN);
  }
  if (stats_io_trim(def.name).empty() || stats_io_trim(def.grid).empty())
    return;
  const auto token = stats_io_trim(def.grid);
  if (token == "zt" || token == "nzt")
    def.grid_id = GRID_ZT;
  else if (token == "lh_zt")
    def.grid_id = GRID_LH_ZT;
  else if (token == "zm" || token == "lh_zm" || token == "nzm")
    def.grid_id = GRID_ZM;
  else if (token == "rad_zt")
    def.grid_id = GRID_RAD_ZT;
  else if (token == "rad_zm")
    def.grid_id = GRID_RAD_ZM;
  else if (token == "sfc")
    def.grid_id = GRID_SFC;
  else if (token == "lh_sfc")
    def.grid_id = GRID_LH_SFC;
  else
    def.grid_id = GRID_ZT;
  l_valid = true;
}

// -----------------------------------------------------------------------------
void split_registry_fields(const std::string& line,                          // In
                           std::vector<std::string>& fields, int& nfields) { // InOut, Out
  // Description:
  //   Split a pipe-delimited line, trim blanks and retain the caller's field
  //   count/capacity. The sole Native caller uses four256-character fields.
  // line/fields: Raw text and fixed-width parsed fields [-]. nfields: Total count [-].
  int lenline, istart;
  char ch;
  nfields = 0;
  std::fill(fields.begin(), fields.end(), "");
  lenline = static_cast<int>(stats_io_trim(line).size());
  if (lenline < 1)
    return;
  istart = 0;
  for (int i = 0; i < lenline; ++i) {
    ch = line[i];
    if (ch == '|') {
      nfields = (nfields + 1);
      if (nfields <= static_cast<int>(fields.size())) {
        if (i > istart)
          fields[(nfields - 1)] = stats_io_character(
              stats_io_trim(stats_io_adjustl(line.substr(istart, i - istart))), 256);
        else
          fields[(nfields - 1)] = "";
      }
      istart = i + 1;
    }
  }
  nfields = (nfields + 1);
  if (nfields <= static_cast<int>(fields.size())) {
    if (lenline > istart)
      fields[(nfields - 1)] = stats_io_character(
          stats_io_trim(stats_io_adjustl(line.substr(istart, lenline - istart))), 256);
    else
      fields[(nfields - 1)] = "";
  }
}

// -----------------------------------------------------------------------------
bool is_line_comment(const std::string& line) { // In
  // Description:
  //   Recognize blank lines or a first-character #/! marker in registry text.
  //   Leading blanks before a marker are significant, as in the Native body.
  if (stats_io_trim(line).empty())
    return true;
  const char c = line[0];
  return c == '#' || c == '!';
}
} // namespace stats_netcdf_detail

} // namespace clubb_kokkos
