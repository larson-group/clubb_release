#include "CLUBB_core/stats_netcdf.hpp"
#include "CLUBB_core/error_code.hpp"
#include "netcdf.hpp"

#include <Kokkos_Core.hpp>

#include <fstream>
#include <iostream>
#include <stdexcept>

extern "C" int native_verify_stats_output_c();

namespace {

void check_error(const clubb_kokkos::err_info_type& err_info) {
  // Detect initialization, batch-selection, write or finalization errors.
  auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, err_info.err_code);
  for (std::size_t column = 0; column < host.extent(0); ++column) {
    if (host(column) == clubb_kokkos::clubb_fatal_error) {
      throw std::runtime_error("C++ stats reported a fatal error");
    }
  }
}

void check_lookup_and_host_snapshot() {
  // Exercise slot zero, ordered/reordered cache paths, cache growth, and an
  // overlapping HostSpace input whose traversal requires an independent copy.
  using namespace clubb_kokkos;
  using stats_netcdf_detail::stats_get_id;
  stats_type stats;
  if (stats_get_id(stats, "field") != -1 || stats.lookup.head != 0 || stats.lookup.reject_head != 0)
    throw std::runtime_error("disabled stats lookup zero-coordinate contract");
  stats.enabled = true;
  stats.nvars   = 3;
  stats.vars.resize(3);
  stats.vars_allocated = true;
  stats.vars[0].name   = "field";
  stats.vars[1].name   = "surface";
  stats.vars[2].name   = "budget";
  if (stats_get_id(stats, "field") != 0 || stats_get_id(stats, "surface") != 1 ||
      stats_get_id(stats, "field") != 0 || stats.lookup.cache_len != 3)
    throw std::runtime_error("stats registry/cache slots must be zero-based");
  stats.lookup.head = 0;
  if (stats_get_id(stats, "budget") != 2 || stats.lookup.head != 0 ||
      stats_get_id(stats, "field") != 0 || stats.lookup.head != 1)
    throw std::runtime_error("stats reordered/expected hit-cache paths changed");
  if (stats_get_id(stats, "not_registered") != -1 || stats_get_id(stats, "not_registered") != -1 ||
      stats.lookup.reject_head != 1)
    throw std::runtime_error("stats reject cache slot/sentinel contract");
  stats.lookup.head = stats.lookup.cache_len;
  for (int i = 0; i < 32; ++i)
    if (stats_get_id(stats, "field") != 0)
      throw std::runtime_error("stats hit-cache growth changed a slot");
  for (int i = 3; i < stats.lookup.cache_len; ++i)
    if (stats.lookup.cache(i) != 0)
      throw std::runtime_error("stats hit-cache growth changed stored zero slots");
  stats.l_sample          = true;
  stats.samples_per_write = 1;
  auto& variable          = stats.vars[0];
  variable.buffer         = StatsReal2D("stats host alias buffer", 2, 2);
  variable.nsamples       = StatsHostInt2D("stats host alias counts", 2, 2);
  for (int k = 0; k < 2; ++k)
    for (int c = 0; c < 2; ++c)
      variable.buffer(c, k) = stat_t(1 + c + 2 * k);
  using Transpose = Kokkos::View<stat_t**, Kokkos::LayoutStride, Kokkos::HostSpace>;
  Transpose values(variable.buffer.data(), Kokkos::LayoutStride(2, 2, 2, 1));
  stats_update("field", values, stats);
  if (variable.buffer(0, 0) != 2 || variable.buffer(1, 0) != 5 || variable.buffer(0, 1) != 5 ||
      variable.buffer(1, 1) != 8)
    throw std::runtime_error("Host stats inputs must retain independent snapshot values");
  err_info_type error;
  stats_finalize_api(stats, error);
  if (stats.workspace().reserved_device_bytes() != 0 ||
      stats.workspace().reserved_host_bytes() != 0)
    throw std::runtime_error("stats finalization retained temporary workspace chunks");
}

void check_output_grid_remap(clubb_kokkos::err_info_type& err_info) {
  // Write two source cells into one output cell. Unit density gives equal
  // masses, so independent columns must contain exactly (2+6)/2=4 and
  // (4+12)/2=8. This checks source/output dimensions, the remap metadata update,
  // source-grid retrieval, and the saved remapped slab; surface fields stay 9/11.
  using namespace clubb_kokkos;
  std::ofstream registry("stats_netcdf_remap_test.in");
  registry << "&clubb_stats_nl\n"
              " entry(1) = 'field | zt | K | Remapped profile'\n"
              " entry(2) = 'surface | sfc | K | Unremapped surface'\n/\n";
  registry.close();

  StatsHostReal1D zt("remap_source_zt", 2), zm("remap_source_zm", 3);
  zt(0) = 32;
  zt(1) = 96;
  zm(0) = 0;
  zm(1) = 64;
  zm(2) = 128;
  StatsHostReal1D output_zt("remap_output_zt", 1), output_zm("remap_output_zm", 2);
  output_zt(0) = 64;
  output_zm(0) = 0;
  output_zm(1) = 128;
  StatsInitOptions options;
  options.output_zt         = output_zt;
  options.output_zm         = output_zm;
  options.grid_remap_method = 1;
  options.stats_tstart      = 0;
  options.stats_tend        = 1;
  stats_type stats;
  stats_init_api("stats_netcdf_remap_test.in", "stats_netcdf_remap.nc", 2, 1, 1, 1, 1, 1, 2000, 0,
                 zt, zm, stats, err_info, options);
  check_error(err_info);
  if (!stats.l_different_output_grid || stats.vars[0].nz != 2 || stats.vars[0].out_nz != 1)
    throw std::runtime_error("stats remap source/output dimensions");

  StatsHostReal2D source_zt("remap_source_zt_columns", 2, 2);
  StatsHostReal2D source_zm("remap_source_zm_columns", 2, 3);
  StatsHostReal2D rho("remap_rho", 2, 2), rho_levels("remap_rho_levels", 2, 2);
  StatsHostReal1D p_sfc("remap_pressure_reference", 2);
  for (int c = 0; c < 2; ++c) {
    for (int k = 0; k < 2; ++k) {
      source_zt(c, k)  = zt(k);
      rho(c, k)        = 1;
      rho_levels(c, k) = real_t(128 * k);
    }
    for (int k = 0; k < 3; ++k)
      source_zm(c, k) = zm(k);
    // Ullrich weights depend only on pressure differences. A zero reference
    // makes these artificial equal-mass cells yield exact binary weights 1/2,
    // isolating stats routing from cancellation in a large pressure offset.
    p_sfc(c) = 0;
  }
  stats_update_grid(source_zt, source_zm, rho, rho_levels, p_sfc, stats);
  StatsHostReal2D returned_zt, returned_zm;
  stats_get_source_grid_api(stats, returned_zt, returned_zm);
  if (returned_zt.extent_int(1) != 2 || returned_zm.extent_int(1) != 3)
    throw std::runtime_error("stats remap returned the output grid as source");
  for (int c = 0; c < 2; ++c) {
    for (int k = 0; k < 2; ++k)
      if (returned_zt(c, k) != source_zt(c, k))
        throw std::runtime_error("stats remap source thermodynamic heights");
    for (int k = 0; k < 3; ++k)
      if (returned_zm(c, k) != source_zm(c, k))
        throw std::runtime_error("stats remap source momentum heights");
  }
  StatsHostReal2D field("remap_profile", 2, 2);
  StatsHostReal1D surface("remap_surface", 2);
  for (int c = 0; c < 2; ++c) {
    field(c, 0) = real_t(2 * (c + 1));
    field(c, 1) = real_t(6 * (c + 1));
    surface(c)  = real_t(9 + 2 * c);
  }
  stats_begin_timestep_api(1, stats);
  stats_update("field", field, stats);
  stats_update("surface", surface, stats);
  stats_end_timestep_api(real_t(1), stats, err_info);
  stats_finalize_api(stats, err_info);
  check_error(err_info);

  const auto require_netcdf = [](int status) {
    // Make failures in this edge's saved-file reads visible to CTest.
    if (status != netcdf::nc_noerr)
      throw std::runtime_error(netcdf::strerror(status));
  };
  int ncid, dimid, length, varid;
  char dimension_name[32];
  require_netcdf(netcdf::open("stats_netcdf_remap.nc", ncid));
  require_netcdf(netcdf::inq_dimid(ncid, "zt", dimid));
  require_netcdf(netcdf::inquire_dimension(ncid, dimid, dimension_name, 32, length));
  if (length != 1)
    throw std::runtime_error("stats remap saved the source zt dimension");
  require_netcdf(netcdf::inq_dimid(ncid, "zm", dimid));
  require_netcdf(netcdf::inquire_dimension(ncid, dimid, dimension_name, 32, length));
  if (length != 2)
    throw std::runtime_error("stats remap saved the source zm dimension");
  const int start[] = {1, 1, 1}, count[] = {2, 1, 1};
  double saved[2];
  require_netcdf(netcdf::inq_varid(ncid, "field", varid));
  require_netcdf(netcdf::get_var(ncid, varid, 2, 3, start, count, saved));
  if (saved[0] != 4.0 || saved[1] != 8.0)
    throw std::runtime_error("stats remap equal-mass saved profile");
  require_netcdf(netcdf::inq_varid(ncid, "surface", varid));
  require_netcdf(netcdf::get_var(ncid, varid, 2, 2, start, count, saved));
  if (saved[0] != 9.0 || saved[1] != 11.0)
    throw std::runtime_error("stats remap changed the surface field");
  require_netcdf(netcdf::close(ncid));
}

} // namespace

int main(int argc, char* argv[]) {
  // Send GPU fields through two simultaneous C++ contexts. Check averaged
  // NetCDF averages/budget tendencies, batch columns, scalar/profile dispatch
  // and context isolation. Sample output also verifies all four array dimensions.
  // Snapshots permit registry queries but cannot update the owner's live I/O.
  try {
    Kokkos::ScopeGuard kokkos(argc, argv);
    using clubb_kokkos::real_t;
    check_lookup_and_host_snapshot();
    std::ofstream registry("stats_netcdf_test.in");
    registry << "&clubb_stats_nl\n"
                " entry(1) = 'field | zt | K | Test profile'\n"
                " entry(2) = 'surface | sfc | K | Test surface'\n"
                " entry(3) = 'budget | zt | K/s | Test tendency budget'\n"
                " entry(4) = 'substep | sfc | K | Sub-timestep mean'\n/\n";
    registry.close();

    Kokkos::View<real_t*, Kokkos::LayoutLeft> zt("zt", 2), zm("zm", 3);
    auto zt_host = Kokkos::create_mirror_view(zt);
    auto zm_host = Kokkos::create_mirror_view(zm);
    zt_host(0)   = 50;
    zt_host(1)   = 150;
    zm_host(0)   = 0;
    zm_host(1)   = 100;
    zm_host(2)   = 200;
    Kokkos::deep_copy(zt, zt_host);
    Kokkos::deep_copy(zm, zm_host);

    clubb_kokkos::err_info_type err_info;
    err_info.err_code = decltype(err_info.err_code)("stats_error", 2);
    Kokkos::deep_copy(err_info.err_code, clubb_kokkos::clubb_no_error);
    clubb_kokkos::stats_type a, b;
    if (a.enabled || b.enabled) {
      throw std::runtime_error("C++ stats context ownership");
    }

    {
      clubb_kokkos::stats_type snapshot;
      clubb_kokkos::copy_stats_snapshot(a, snapshot);
      if (clubb_kokkos::var_on_stats_list(snapshot, "field")) {
        throw std::runtime_error("disabled stats snapshot has an active registry");
      }
    }

    Kokkos::View<real_t**, Kokkos::LayoutLeft> parameters("parameters", 4, 1);
    Kokkos::deep_copy(parameters, real_t(7));
    clubb_kokkos::StatsInitOptions options;
    options.clubb_params  = clubb_kokkos::to_host_space(parameters);
    options.param_names   = std::vector<std::string>{"test_parameter"};
    options.hydromet_list = std::vector<std::string>{"rrm"};
    options.sclr_dim      = 0;
    options.edsclr_dim    = 0;
    options.ncol_total    = 4;
    options.stats_tstart  = 0;
    options.stats_tend    = 2;
    clubb_kokkos::stats_init_api("stats_netcdf_test.in", "stats_netcdf_a.nc", 2, 1, 2, 1, 1, 1,
                                 2000, 0, zt, zm, a, err_info, options);
    check_error(err_info);

    {
      clubb_kokkos::stats_type snapshot;
      clubb_kokkos::copy_stats_snapshot(a, snapshot);
      if (!snapshot.enabled)
        throw std::runtime_error("stats snapshot lost enabled state");
      if (snapshot.vars[0].buffer.data() == a.vars[0].buffer.data())
        throw std::runtime_error("stats snapshot aliases its owner buffers");
      snapshot.vars[0].buffer(0, 0) = 123;
      if (a.vars[0].buffer(0, 0) != 0)
        throw std::runtime_error("stats snapshot mutated its owner buffers");
      if (!clubb_kokkos::var_on_stats_list(snapshot, "field") ||
          clubb_kokkos::var_on_stats_list(snapshot, "not_registered") ||
          clubb_kokkos::var_on_stats_list(snapshot, "")) {
        throw std::runtime_error("stats snapshot registry query");
      }
      bool rejected = false;
      try {
        clubb_kokkos::stats_reset_api(snapshot);
      } catch (const std::logic_error&) {
        rejected = true;
      }
      if (!rejected)
        throw std::runtime_error("stats snapshot exposed mutable I/O ownership");
      rejected = false;
      try {
        (void)snapshot.workspace();
      } catch (const std::logic_error&) {
        rejected = true;
      }
      if (!rejected)
        throw std::runtime_error("stats snapshot exposed a mutable temporary pool");
    }
    // Later writes and final NetCDF verification also establish that destroying
    // the snapshot did not close the original context's file.
    clubb_kokkos::stats_init_api("stats_netcdf_test.in", "stats_netcdf_b.nc", 2, 1, 2, 1, 1, 1,
                                 2000, 0, zt, zm, b, err_info);
    check_error(err_info);

    // Native defines sample coordinates even when both output name lists are
    // empty. A separate context checks populated names and sample-field writes.
    clubb_kokkos::stats_lh_samples_init(2, 2, {}, {}, zt, a, err_info);
    clubb_kokkos::stats_lh_samples_init(2, 2, {"scalar"}, {"scalar", "component", "weight"}, zt, b,
                                        err_info);
    check_error(err_info);
    clubb_kokkos::Real4D nl_samples("nl_samples", 2, 2, 2, 1);
    clubb_kokkos::Real4D uniform_vals("uniform_vals", 2, 2, 2, 1);
    clubb_kokkos::Real3D sample_weights("sample_weights", 2, 2, 2);
    Kokkos::View<int***, Kokkos::LayoutLeft> mixture_comp("mixture_comp", 2, 2, 2);

    Kokkos::View<real_t**, Kokkos::LayoutLeft> field("field", 2, 2);
    Kokkos::View<real_t**, Kokkos::LayoutLeft> budget_before("budget_before", 2, 2);
    Kokkos::View<real_t**, Kokkos::LayoutLeft> budget_after("budget_after", 2, 2);
    Kokkos::deep_copy(budget_before, real_t(3));
    Kokkos::deep_copy(budget_after, real_t(5));
    Kokkos::View<real_t*, Kokkos::LayoutLeft> surface("surface", 2), profile("profile", 2);
    // An unsampled or unregistered field must return before staging its payload.
    clubb_kokkos::stats_update("field", clubb_kokkos::Real2D{}, a);
    for (int batch = 0; batch < 2; ++batch) {
      clubb_kokkos::stats_reset_api(a);
      if (batch > 0)
        clubb_kokkos::start_next_stats_batch_api(a, err_info);
      check_error(err_info);
      for (int step = 1; step <= 2; ++step) {
        const int base = batch * 100 + 2 * (step - 1);
        Kokkos::parallel_for(
            "stats_field_values", Kokkos::RangePolicy<>(0, 4), KOKKOS_LAMBDA(const int linear) {
              const int column     = linear % 2;
              const int level      = linear / 2;
              field(column, level) = real_t(base + column + 10 * level);
              if (level == 0)
                surface(column) = real_t(base + column);
            });
        clubb_kokkos::stats_begin_timestep_api(step, a);
        if (!a.l_sample || a.l_last_sample != (step == 2)) {
          throw std::runtime_error("C++ stats sampling cadence");
        }
        clubb_kokkos::stats_update("not_registered", clubb_kokkos::Real2D{}, a);
        // Native budgets subtract the initial state and add the final state.
        // Every sampled cell must average to exactly 5 - 3 = 2.
        clubb_kokkos::stats_begin_budget("budget", budget_before, a);
        if (batch == 0) {
          clubb_kokkos::stats_update("field", field, a);
        } else {
          const auto host_field = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, field);
          clubb_kokkos::stats_update("field", host_field, a);
        }
        // Exercise direct execution-space and host-mirror surface Views through
        // the shared rank-one adapter; the final NetCDF check verifies both.
        if (batch == 0) {
          clubb_kokkos::stats_update("surface", surface, a);
        } else {
          const auto host_surface =
              Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, surface);
          clubb_kokkos::stats_update("surface", host_surface, a);
        }
        // The rank-one fourth optional is a column index. An absent index and
        // explicit fifth argument request the within-timestep mean for all columns.
        Kokkos::deep_copy(surface, real_t(2));
        clubb_kokkos::stats_update("substep", surface, a, std::nullopt, true);
        Kokkos::deep_copy(surface, real_t(6));
        clubb_kokkos::stats_update("substep", surface, a, std::nullopt, true);
        clubb_kokkos::stats_finalize_budget("budget", budget_after, a);
        clubb_kokkos::stats_end_timestep_api(real_t(step), a, err_info);
        check_error(err_info);

        if (batch == 0) {
          const int bbase = 200 + 2 * (step - 1);
          Kokkos::parallel_for(
              "stats_profile_values", Kokkos::RangePolicy<>(0, 2),
              KOKKOS_LAMBDA(const int level) { profile(level) = real_t(bbase + 10 * level); });
          clubb_kokkos::stats_begin_timestep_api(step, b);
          clubb_kokkos::stats_update("field", profile, b, 0);
          for (int level = 0; level < 2; ++level) {
            clubb_kokkos::stats_update("field", real_t(bbase + 1 + 10 * level), b, 1, level);
          }
          clubb_kokkos::stats_update("surface", real_t(bbase), b, 0);
          clubb_kokkos::stats_update("surface", real_t(bbase + 1), b, 1);
          Kokkos::parallel_for(
              "stats_sample_values", Kokkos::RangePolicy<>(0, 8), KOKKOS_LAMBDA(const int linear) {
                const int column    = linear % 2;
                const int sample    = (linear / 2) % 2;
                const int level     = linear / 4;
                const real_t offset = real_t(10 * step + column + 2 * sample + 4 * level);
                nl_samples(column, sample, level, 0)   = offset;
                uniform_vals(column, sample, level, 0) = real_t(20) + offset;
                sample_weights(column, sample, level)  = real_t(40) + offset;
                mixture_comp(column, sample, level)    = 1 + column;
              });
          clubb_kokkos::stats_lh_samples_write_lognormal(nl_samples, b, err_info);
          clubb_kokkos::stats_lh_samples_write_uniform(uniform_vals, mixture_comp, sample_weights,
                                                       b, err_info);
          clubb_kokkos::stats_end_timestep_api(real_t(step), b, err_info);
          check_error(err_info);
        }
      }
    }

    clubb_kokkos::stats_finalize_api(a, err_info);
    clubb_kokkos::stats_finalize_api(b, err_info);
    check_error(err_info);
    if (a.enabled || b.enabled || native_verify_stats_output_c() != 0) {
      throw std::runtime_error("C++ stats NetCDF verification");
    }
    check_output_grid_remap(err_info);
    clubb_kokkos::cleanup_err_info_api(err_info);
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
