// Construct staggered grids and apply their operators; see grid_class.hpp for the public APIs.
#include "grid_class.hpp"
#include "file_functions.hpp"
#include <array>
#include "constants_clubb.hpp"
#include "error_code.hpp"
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include "interpolation.hpp"
#include "model_flags.hpp"

#include <stdexcept>

namespace clubb_kokkos {

// Native private band selectors use zero-based C++ indices.
static constexpr int t_above = 0, t_below = 1, m_above = 0, m_below = 1;

static Real2D redirect_interpolated_azm_2D(int nzm, int nzt, int ngrdcol, const grid& gr,
                                           Real2D azt, std::optional<real_t> zm_min);
static Real2D redirect_interpolated_azt_2D(int nzm, int nzt, int ngrdcol, const grid& gr,
                                           Real2D azm, std::optional<real_t> zt_min);
static void linear_interpolated_azm_2D(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt,
                                       Real2D linear_interpolated_azm,
                                       std::optional<real_t> zm_min);
static void cubic_interpolated_azm_2D(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt,
                                      Real2D cubic_interpolated_azm_2D);
static void linear_interpolated_azt_2D(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm,
                                       Real2D linear_interpolated_azt,
                                       std::optional<real_t> zt_min);
static void cubic_interpolated_azt_2D(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm,
                                      Real2D cubic_interpolated_azt_2D);
static void gradzm_2D(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm, Real2D gradzm_2D);
static void gradzt_2D(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt, Real2D gradzt_2D);

static void calc_zt2zm_weights(int nzm, int nzt, int ngrdcol, // In
                               grid& gr);                     // InOut
static void calc_zm2zt_weights(int nzt, int ngrdcol,          // In
                               grid& gr);                     // InOut

static bool overlapping_fields(Real2D source, Real2D output) {
  // Preserve returned-field snapshot semantics when an explicit output overlaps.
  if (source.size() == 0 || output.size() == 0)
    return false;
  const auto source_begin = reinterpret_cast<std::uintptr_t>(source.data());
  const auto output_begin = reinterpret_cast<std::uintptr_t>(output.data());
  if (source_begin <= output_begin)
    return output_begin - source_begin < source.span() * sizeof(real_t);
  return source_begin - output_begin < output.span() * sizeof(real_t);
}

// =============================================================================
void setup_grid(int nzmax, int ngrdcol, Real1D sfc_elevation, bool l_implemented, // In
                bool l_ascending_grid, int grid_type, Real1D deltaz, Real1D zm_init,
                Real1D zm_top,                                         // In
                Real2D momentum_heights, Real2D thermodynamic_heights, // In
                grid& gr, err_info_type& err_info                      // InOut
) {
  // Description:
  //   Grid Constructor
  //
  //   This subroutine sets up the CLUBB vertical grid.
  //
  // References:
  //   ``Equations for CLUBB'',  Sec. 8,  Grid Configuration.
  //-----------------------------------------------------------------------
  // nzmax: Declared momentum-level extent before standalone trimming [#].
  // ngrdcol: Number of grid columns [#].
  // sfc_elevation: Ground elevation above mean sea level [m].
  // l_implemented: True when a host model supplies both sets of level heights [-].
  // l_ascending_grid: True for ascending storage, false for descending storage [-].
  // grid_type: Standalone uniform(1), thermodynamic-input(2), or momentum-input(3) grid [-].
  // deltaz: Standalone uniform vertical spacing by column [m].
  // zm_init: Lower-boundary momentum altitude by column [m].
  // zm_top: Maximum requested momentum altitude by column [m].
  // momentum_heights: Input momentum altitudes, shape (ngrdcol,nzmax) [m].
  // thermodynamic_heights: Input thermodynamic altitudes, shape (ngrdcol,nzmax-1) [m].
  // gr: Caller-owned grid geometry, interpolation weights and physical metadata [units vary].
  // err_info: Caller-owned error status and diagnostic headers [-].
  constexpr int NWARNING = 250; // Warn when the declared grid exceeds this extent [#].
  int begin_height_idx_zm, end_height_idx_zm; // Selected inclusive momentum-input slots [-].
  int begin_height_idx_zt, end_height_idx_zt; // Selected inclusive thermodynamic-input slots [-].
  int k; // Zero-based level index used in first-column searches [-].

  if (nzmax > NWARNING && clubb_at_least_debug_level_api(1)) {
    std::cerr << "Warning:  running with vertical grid which is larger than" << NWARNING
              << "levels.\n";
    std::cerr << "This may take a lot of CPU time and memory.\n";
  }

  gr.nzm = nzmax;
  gr.nzt = nzmax - 1;
  begin_height_idx_zm = 0;
  end_height_idx_zm   = gr.nzm - 1;
  begin_height_idx_zt = 0;
  end_height_idx_zt   = gr.nzt - 1;

  if (!l_implemented) {
    // CLUBB cannot use different vertical extents per column. Like Native,
    // determine the standalone extent only from the first input column.
    if (grid_type == 1) {
      const auto zm_top_host = to_host_space(zm_top);
      const auto zm_init_host = to_host_space(zm_init);
      const auto deltaz_host = to_host_space(deltaz);
      gr.nzm = static_cast<int>(
          std::floor((zm_top_host(0) - zm_init_host(0) + deltaz_host(0)) / deltaz_host(0)));
      gr.nzt = gr.nzm - 1;
    } else if (grid_type == 2) {
      const auto thermodynamic_heights_host = to_host_space(thermodynamic_heights);
      const auto zm_init_host = to_host_space(zm_init);
      const auto zm_top_host = to_host_space(zm_top);

      // Locate the input thermodynamic level associated with the lower boundary.
      if (l_ascending_grid) {
        k = gr.nzt - 1;
        while (thermodynamic_heights_host(0, k) >= zm_init_host(0) && k > 0)
          k = k - 1;
      } else {
        k = 0;
        while (thermodynamic_heights_host(0, k) >= zm_init_host(0) && k < gr.nzt - 1)
          k = k + 1;
      }

      if (thermodynamic_heights_host(0, k) <= zm_init_host(0)) {
        const auto header = to_host_space(err_info.err_header);
        std::cerr.write(header.data(), static_cast<std::streamsize>(header.extent(0)));
        std::cerr << '\n' << "Stretched zt grid cannot fulfill zm_init requirement\n";
        std::cerr << "It is recommended to check the zm grid for all columns of the chunk\n";
        Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
        return;
      } else {
        if (l_ascending_grid)
          begin_height_idx_zt = k;
        else
          end_height_idx_zt = k;
      }

      // Locate the input thermodynamic level associated with the upper boundary.
      if (l_ascending_grid) {
        k = gr.nzt - 1;
        while (thermodynamic_heights_host(0, k) > zm_top_host(0) && k > 0)
          k = k - 1;
      } else {
        k = 0;
        while (thermodynamic_heights_host(0, k) > zm_top_host(0) && k < gr.nzt - 1)
          k = k + 1;
      }

      if (zm_top_host(0) < thermodynamic_heights_host(0, k)) {
        const auto header = to_host_space(err_info.err_header);
        std::cerr.write(header.data(), static_cast<std::streamsize>(header.extent(0)));
        std::cerr << '\n' << "Stretched zt grid cannot fulfill zm_top requirement\n";
        std::cerr << "It is recommended to check the zt grid for all columns of the chunk\n";
        Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
        return;
      } else {
        if (l_ascending_grid)
          end_height_idx_zt = k;
        else
          begin_height_idx_zt = k;
        gr.nzt = end_height_idx_zt - begin_height_idx_zt + 1;
        gr.nzm = gr.nzt + 1;
        // Momentum input is unused here; keep Native's consistent paired bounds.
        begin_height_idx_zm = begin_height_idx_zt;
        end_height_idx_zm = end_height_idx_zt + 1;
      }
    } else if (grid_type == 3) {
      const auto momentum_heights_host = to_host_space(momentum_heights);
      const auto zm_init_host = to_host_space(zm_init);
      const auto zm_top_host = to_host_space(zm_top);

      // Locate the input momentum level associated with the lower boundary.
      if (l_ascending_grid) {
        k = 0;
        while (momentum_heights_host(0, k) < zm_init_host(0) && k < gr.nzm - 1)
          k = k + 1;
      } else {
        k = gr.nzm - 1;
        while (momentum_heights_host(0, k) < zm_init_host(0) && k > 0)
          k = k - 1;
      }

      if (momentum_heights_host(0, k) < zm_init_host(0)) {
        const auto header = to_host_space(err_info.err_header);
        std::cerr.write(header.data(), static_cast<std::streamsize>(header.extent(0)));
        std::cerr << '\n' << "Stretched zm grid cannot fulfill zm_init requirement\n";
        std::cerr << "It is recommended to check the zm grid for all columns of the chunk\n";
        Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
        return;
      } else {
        if (l_ascending_grid)
          begin_height_idx_zm = k;
        else
          end_height_idx_zm = k;
      }

      // Locate the input momentum level associated with the upper boundary.
      if (l_ascending_grid) {
        k = gr.nzm - 1;
        while (momentum_heights_host(0, k) > zm_top_host(0) && k > 0)
          k = k - 1;
      } else {
        k = 0;
        while (momentum_heights_host(0, k) > zm_top_host(0) && k < gr.nzm - 1)
          k = k + 1;
      }

      if (momentum_heights_host(0, k) > zm_top_host(0)) {
        const auto header = to_host_space(err_info.err_header);
        std::cerr.write(header.data(), static_cast<std::streamsize>(header.extent(0)));
        std::cerr << '\n' << "Stretched zm grid cannot fulfill zm_top requirement\n";
        std::cerr << "It is recommended to check the zt grid for all columns of the chunk\n";
        Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
        return;
      } else {
        if (l_ascending_grid)
          end_height_idx_zm = k;
        else
          begin_height_idx_zm = k;
        gr.nzm = end_height_idx_zm - begin_height_idx_zm + 1;
        gr.nzt = gr.nzm - 1;
        // Thermodynamic input is unused here; retain the paired Native bounds.
        begin_height_idx_zt = begin_height_idx_zm;
        end_height_idx_zt = end_height_idx_zm - 1;
      }
    }
  }

  // Native ALLOCATE does not replace an allocated component. Guard each object
  // in source order, preserving earlier successful allocations if a later one fails.
  // WithoutInitializing matches Native ALLOCATE and avoids redundant zero fills.
  try {
    if (gr.zm.is_allocated())
      throw std::logic_error("setup_grid: zm already allocated");
    gr.zm = decltype(gr.zm)(Kokkos::view_alloc(Kokkos::WithoutInitializing, "grid::zm"), ngrdcol,
                            gr.nzm);
    if (gr.zt.is_allocated())
      throw std::logic_error("setup_grid: zt already allocated");
    gr.zt = decltype(gr.zt)(Kokkos::view_alloc(Kokkos::WithoutInitializing, "grid::zt"), ngrdcol,
                            gr.nzt);
    if (gr.dzm.is_allocated())
      throw std::logic_error("setup_grid: dzm already allocated");
    gr.dzm = decltype(gr.dzm)(Kokkos::view_alloc(Kokkos::WithoutInitializing, "grid::dzm"), ngrdcol,
                              gr.nzm);
    if (gr.dzt.is_allocated())
      throw std::logic_error("setup_grid: dzt already allocated");
    gr.dzt = decltype(gr.dzt)(Kokkos::view_alloc(Kokkos::WithoutInitializing, "grid::dzt"), ngrdcol,
                              gr.nzt);
    if (gr.invrs_dzm.is_allocated())
      throw std::logic_error("setup_grid: invrs_dzm already allocated");
    gr.invrs_dzm = decltype(gr.invrs_dzm)(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "grid::invrs_dzm"), ngrdcol, gr.nzm);
    if (gr.invrs_dzt.is_allocated())
      throw std::logic_error("setup_grid: invrs_dzt already allocated");
    gr.invrs_dzt = decltype(gr.invrs_dzt)(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "grid::invrs_dzt"), ngrdcol, gr.nzt);
    if (gr.weights_zm2zt.is_allocated())
      throw std::logic_error("setup_grid: weights_zm2zt already allocated");
    gr.weights_zm2zt = decltype(gr.weights_zm2zt)(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "grid::weights_zm2zt"), ngrdcol, gr.nzt, 2);
    if (gr.weights_zt2zm.is_allocated())
      throw std::logic_error("setup_grid: weights_zt2zm already allocated");
    gr.weights_zt2zm = decltype(gr.weights_zt2zm)(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "grid::weights_zt2zm"), ngrdcol, gr.nzm, 2);
  } catch (const std::exception&) {
    std::cerr << err_info.err_header_global << '\n';
    std::cerr << "In setup_grid: allocation of grid variables failed.\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }

  // Shape adapters materialize only the selected Native level sections.
  setup_grid_heights(gr.nzm, gr.nzt, ngrdcol,         // In
                     l_implemented, l_ascending_grid, // In
                     grid_type, deltaz, zm_init,      // In
                     to_execution_space(Kokkos::subview(
                         momentum_heights, Kokkos::ALL,
                         std::make_pair(begin_height_idx_zm, end_height_idx_zm + 1))),
                     to_execution_space(Kokkos::subview(
                         thermodynamic_heights, Kokkos::ALL,
                         std::make_pair(begin_height_idx_zt, end_height_idx_zt + 1))),
                     gr, err_info // InOut
  );

  if (clubb_at_least_debug_level_api(0)) {
    const auto err_code = err_info.err_code;
    int fatal_count = 0;
    Kokkos::parallel_reduce(
        "setup_grid::height_error", ngrdcol,
        KOKKOS_LAMBDA(const int i, int& count) {
          if (err_code(i) == clubb_fatal_error)
            ++count;
        },
        fatal_count);
    if (fatal_count != 0) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr << "Fatal error in setup_grid calling setup_grid_heights\n";
      return;
    }
  }

  // Native checks storage level1 even for descending storage; do not substitute
  // the physical lower-boundary metadata. Report only the first failing column.
  const auto zm = gr.zm;
  int first_bad_column;
  Kokkos::parallel_reduce(
      "setup_grid::surface_elevation", ngrdcol,
      KOKKOS_LAMBDA(const int i, int& first) {
        if (sfc_elevation(i) > zm(i, 0) && i < first)
          first = i;
      },
      Kokkos::Min<int>(first_bad_column));
  if (first_bad_column < ngrdcol) {
    const auto header = to_host_space(err_info.err_header);
    const auto zm_host = to_host_space(zm);
    const auto sfc_elevation_host = to_host_space(sfc_elevation);
    std::ostringstream diagnostic;
    diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    diagnostic.write(header.data() + first_bad_column * header.extent(0),
                     static_cast<std::streamsize>(header.extent(0)));
    diagnostic << '\n';
    diagnostic << "The altitude of the lowest momentum level, gr%zm(1,1), must be at or above the "
                  "altitude of the surface, sfc_elevation.  The lowest model momentum level cannot "
                  "be below the surface.\n";
    diagnostic << "Altitude of lowest momentum level =" << zm_host(first_bad_column, 0) << '\n';
    diagnostic << "Altitude of the surface =" << sfc_elevation_host(first_bad_column) << '\n';
    std::cerr << diagnostic.str();
    Kokkos::deep_copy(Kokkos::subview(err_info.err_code, first_bad_column), clubb_fatal_error);
    return;
  }
}

// =============================================================================
void cleanup_grid_api(grid& gr) { // InOut
  // Description:
  //   Release the caller-owned grid arrays in Native component order.
  // References: None.
  // gr: Grid whose reference-counted Views are released; metadata is retained [units vary].
  // Native DEALLOCATE STAT can print "Grid deallocation failed.". Kokkos View
  // release has no corresponding STAT result; no fabricated status is introduced.
  gr.zm = {};
  gr.zt = {};
  gr.dzm = {};
  gr.dzt = {};
  gr.invrs_dzm = {};
  gr.invrs_dzt = {};
  gr.weights_zm2zt = {};
  gr.weights_zt2zm = {};
}

// =============================================================================
void setup_grid_heights(int nzm, int nzt, int ngrdcol,                         // In
                        bool l_implemented, bool l_ascending_grid,             // In
                        int grid_type, Real1D deltaz, Real1D zm_init,          // In
                        Real2D momentum_heights, Real2D thermodynamic_heights, // In
                        grid& gr, err_info_type& err_info                      // InOut
) {
  // Description:
  //   Sets the heights and interpolation weights of the column.
  //   This is seperated from setup_grid for those host models that have heights
  //   that vary with time.
  // References:
  //   None
  //------------------------------------------------------------------------------
  // nzm: Number of momentum levels [#].
  // nzt: Number of thermodynamic levels [#].
  // ngrdcol: Number of grid columns [#].
  // l_implemented: Use supplied host momentum and thermodynamic heights [-].
  // l_ascending_grid: Ascending or descending grid-storage convention [-].
  // grid_type: Standalone uniform(1), thermodynamic-input(2), or momentum-input(3) grid [-].
  // deltaz: Uniform vertical spacing by column [m].
  // zm_init: Lower-boundary momentum altitude by column [m].
  // momentum_heights: Input momentum altitudes, shape (ngrdcol,nzm) [m].
  // thermodynamic_heights: Input thermodynamic altitudes, shape (ngrdcol,nzt) [m].
  // gr: Existing allocated grid fields and physical boundary metadata [units vary].
  // err_info: Caller-owned error status and diagnostic headers [-].

  // Native's invalid standalone branch writes no geometry before returning.
  if (!l_implemented && grid_type != 1 && grid_type != 2 && grid_type != 3) {
    std::cerr << err_info.err_header_global << '\n';
    std::cerr << "Invalid grid type: " << grid_type << ".  Valid options are 1, 2, or 3.\n";
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }

  // Keep the owning grid on the host; capture only its component Views.
  const auto zm = gr.zm;
  const auto zt = gr.zt;
  const auto dzm = gr.dzm;
  const auto dzt = gr.dzt;
  const auto invrs_dzm = gr.invrs_dzm;
  const auto invrs_dzt = gr.invrs_dzt;
  // Each column's vertical work respects the Native geometry/spacing order.
  // Columns are independent, so a serial level loop per column also keeps
  // adjacent columns contiguous across device work items.
  Kokkos::parallel_for(
      "setup_grid_heights", ngrdcol, KOKKOS_LAMBDA(const int i) {
        if (!l_implemented) {
          if (grid_type == 1) {
            // Uniform momentum heights, followed by centered thermodynamic heights.
            if (l_ascending_grid) {
              for (int k = 0; k < nzm; ++k)
                zm(i, k) = zm_init(i) + real_t(k) * deltaz(i);
            } else {
              for (int k = 0; k < nzm; ++k)
                zm(i, k) = zm_init(i) + real_t(nzm - k - 1) * deltaz(i);
            }
            for (int k = 0; k < nzt; ++k)
              zt(i, k) = real_t(0.5) * (zm(i, k + 1) + zm(i, k));
          } else if (grid_type == 2) {
            // Stretched thermodynamic input; construct internal momentum midpoints.
            for (int k = 0; k < nzt; ++k)
              zt(i, k) = thermodynamic_heights(i, k);
            for (int k = 1; k < nzm - 1; ++k)
              zm(i, k) = real_t(0.5) * (zt(i, k) + zt(i, k - 1));
            if (l_ascending_grid) {
              zm(i, 0) = zm_init(i);
              zm(i, nzm - 1) = zt(i, nzt - 1) + real_t(0.5) * (zt(i, nzt - 1) - zt(i, nzt - 2));
            } else {
              zm(i, nzm - 1) = zm_init(i);
              zm(i, 0) = zt(i, 0) + real_t(0.5) * (zt(i, 0) - zt(i, 1));
            }
          } else if (grid_type == 3) {
            // Stretched momentum input, followed by thermodynamic midpoints.
            for (int k = 0; k < nzm; ++k)
              zm(i, k) = momentum_heights(i, k);
            for (int k = 0; k < nzt; ++k)
              zt(i, k) = real_t(0.5) * (zm(i, k + 1) + zm(i, k));
          }
        } else {
          // Implemented host supplies both staggered sets of heights.
          for (int k = 0; k < nzm; ++k)
            zm(i, k) = momentum_heights(i, k);
          for (int k = 0; k < nzt; ++k)
            zt(i, k) = thermodynamic_heights(i, k);
        }

        // Thermodynamic spacing centered on momentum levels; descending values
        // remain negative as in Native.
        for (int k = 1; k < nzm - 1; ++k)
          dzm(i, k) = zt(i, k) - zt(i, k - 1);
        if (l_ascending_grid) {
          dzm(i, 0) = real_t(2) * (zt(i, 0) - zm(i, 0));
          dzm(i, nzm - 1) = dzm(i, nzm - 2);
        } else {
          dzm(i, nzm - 1) = real_t(2) * (zm(i, nzm - 1) - zt(i, nzt - 1));
          dzm(i, 0) = dzm(i, 1);
        }

        // Momentum spacing centered on thermodynamic levels.
        for (int k = 0; k < nzt; ++k)
          dzt(i, k) = zm(i, k + 1) - zm(i, k);
        for (int k = 0; k < nzm; ++k)
          invrs_dzm(i, k) = one / dzm(i, k);
        for (int k = 0; k < nzt; ++k)
          invrs_dzt(i, k) = one / dzt(i, k);
      });

  // Preserve the Native call order before assigning direction/boundary metadata.
  calc_zm2zt_weights(nzt, ngrdcol,      // In
                     gr);               // InOut
  calc_zt2zm_weights(nzm, nzt, ngrdcol, // In
                     gr);               // InOut

  if (l_ascending_grid) {
    gr.k_lb_zm       = 0;
    gr.k_lb_zt       = 0;
    gr.k_ub_zm       = nzm - 1;
    gr.k_ub_zt       = nzt - 1;
    gr.grid_dir_indx = 1;
    gr.grid_dir = one;
  } else {
    gr.k_lb_zm       = nzm - 1;
    gr.k_lb_zt       = nzt - 1;
    gr.k_ub_zm       = 0;
    gr.k_ub_zt       = 0;
    gr.grid_dir_indx = -1;
    gr.grid_dir = -one;
  }
}

// =============================================================================
void read_grid_heights(
    int nzmax, int grid_type,                                                           // In
    const std::string& zm_grid_fname,                                                   // In
    const std::string& zt_grid_fname, int file_unit,                                    // In
    Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> momentum_heights,      // Out
    Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> thermodynamic_heights, // Out
    err_info_type& err_info) {                                                          // InOut
  // Description:
  //   Read stretched-grid altitude levels and check the input for all three
  //   grid types. Type2 supplies thermodynamic heights; type3 supplies momentum
  //   heights. The unused output remains zero. Heights are in metres.
  //------------------------------------------------------------------------
  // nzmax: Declared number of momentum levels [#].
  // grid_type: Uniform(1), stretched thermodynamic(2), or stretched momentum(3) [-].
  // zm_grid_fname/zt_grid_fname: Input paths for momentum/thermodynamic altitudes.
  // file_unit: Required Native logical unit; this call owns its file handle [-].
  // momentum_heights: Momentum input altitudes, extent nzmax [m].
  // thermodynamic_heights: Thermodynamic input altitudes, extent nzmax-1 [m].
  // err_info: Caller-owned diagnostic header and per-column error status [-].
  using file_functions_detail::formatted_file;
  using file_functions_detail::read_item;
  int zt_level_count, zm_level_count;
  int input_status;
  real_t generic_input_item = real_t(0);
  int k;

  // Declare the momentum and thermodynamic level altitudes zero until overwritten.
  Kokkos::deep_copy(Kokkos::subview(momentum_heights, std::make_pair(0, nzmax)), real_t(0));
  Kokkos::deep_copy(Kokkos::subview(thermodynamic_heights, std::make_pair(0, nzmax - 1)),
                    real_t(0));

  if (grid_type == 1) {
    // Evenly-spaced grid: heights follow constant spacing and a bottom starting point.
    // Neither input filename may be set. Fortran blank-pads character comparisons.
    if (zm_grid_fname.find_first_not_of(' ') != std::string::npos) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr << "An evenly-spaced grid has been selected.  Please reset zm_grid_fname to ''.\n";
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
    if (zt_grid_fname.find_first_not_of(' ') != std::string::npos) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr << "An evenly-spaced grid has been selected.  Please reset zt_grid_fname to ''.\n";
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
  } else if (grid_type == 2) {
    // Stretched thermodynamic levels determine halfway momentum levels,
    // following the SAM-style grid. No momentum input file may be selected.
    if (zm_grid_fname.find_first_not_of(' ') != std::string::npos) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr << "Thermodynamic level altitudes have been selected "
                   "for use in a stretched (unevenly-spaced) grid. "
                   " Please reset zm_grid_fname to ''.\n";
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
    {
      // Count thermodynamic altitudes through one scalar READ per record.
      formatted_file file(zt_grid_fname);
      zt_level_count = 0;
      while (true) {
        std::array<read_item, 1> item{read_item(generic_input_item)};
        input_status = file.read(item);
        if (input_status < 0)
          break;
        if (input_status > 0) {
          std::cerr << "Error reading thermodynamic level input file.\n";
          Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
          break;
        }
        ++zt_level_count;
      }
    }
    // Native also returns here when a preexisting column already has fatal status.
    const bool fatal_error = any_equal(err_info.err_code, clubb_fatal_error);
    if (fatal_error) {
      std::cerr << err_info.err_header_global << '\n';
      return;
    }
    if (zt_level_count != nzmax - 1) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr
          << "The number of thermodynamic grid altitudes listed in file "
          << zt_grid_fname.substr(0, zt_grid_fname.find_last_not_of(' ') + 1)
          << " does not match the number of CLUBB grid levels specified in the model.in file.\n";
      std::cerr << "Number of thermodynamic grid altitudes listed:  " << zt_level_count << '\n';
      std::cerr << "Number of CLUBB thermo grid levels specified (nzmax-1):  " << nzmax - 1 << '\n';
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
    // Read the thermodynamic altitudes after checking the declared level count.
    file_read_1d(file_unit, zt_grid_fname, nzmax - 1, 1, // In
                 thermodynamic_heights);                 // Out
    for (k = 1; k <= nzmax - 2; ++k) {
      if (thermodynamic_heights(k) <= thermodynamic_heights(k - 1)) {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << err_info.err_header_global << '\n';
        diagnostic << "The declared thermodynamic level grid altitudes are not increasing "
                      "in height as grid level index increases.\n";
        diagnostic << "Grid index:  " << k
                   << ";  Thermodynamic level altitude:  " << thermodynamic_heights(k - 1) << '\n';
        diagnostic << "Grid index:  " << k + 1
                   << ";  Thermodynamic level altitude:  " << thermodynamic_heights(k) << '\n';
        std::cerr << diagnostic.str();
        Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
        return;
      }
    }
  } else if (grid_type == 3) {
    // Stretched momentum levels determine halfway thermodynamic levels,
    // following the WRF-style grid. No thermodynamic input file may be selected.
    if (zt_grid_fname.find_first_not_of(' ') != std::string::npos) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr << "Momentum level altitudes have been selected "
                   "for use in a stretched (unevenly-spaced) grid. "
                   " Please reset zt_grid_fname to ''.\n";
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
    {
      formatted_file file(zm_grid_fname);
      // Find the number of momentum level altitudes listed in the file.
      zm_level_count = 0;
      while (true) {
        std::array<read_item, 1> item{read_item(generic_input_item)};
        input_status = file.read(item);
        if (input_status < 0)
          break;
        if (input_status > 0) {
          std::cerr << err_info.err_header_global << '\n';
          std::cerr << "Error reading momentum level input file.\n";
          Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
          return;
        }
        ++zm_level_count;
      }
    }
    if (zm_level_count != nzmax) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr
          << "The number of momentum grid altitudes listed in file "
          << zm_grid_fname.substr(0, zm_grid_fname.find_last_not_of(' ') + 1)
          << " does not match the number of CLUBB grid levels specified in the model.in file.\n";
      std::cerr << "Number of momentum grid altitudes listed:  " << zm_level_count << '\n';
      std::cerr << "Number of CLUBB grid levels specified:  " << nzmax << '\n';
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
    // Read the momentum altitudes after checking the declared level count.
    file_read_1d(file_unit, zm_grid_fname, nzmax, 1, // In
                 momentum_heights);                  // Out
    for (k = 1; k < nzmax; ++k) {
      if (momentum_heights(k) <= momentum_heights(k - 1)) {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << err_info.err_header_global << '\n';
        diagnostic << "The declared momentum level grid altitudes are not increasing "
                      "in height as grid level index increases.\n";
        diagnostic << "Grid index:  " << k
                   << ";  Momentum level altitude:  " << momentum_heights(k - 1) << '\n';
        diagnostic << "Grid index:  " << k + 1
                   << ";  Momentum level altitude:  " << momentum_heights(k) << '\n';
        std::cerr << diagnostic.str();
        Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
        return;
      }
    }
  }
}

// =============================================================================
static real_t redirect_interpolated_azm_k(const grid& gr, Real1D azt, int k,
                                          std::optional<real_t> zm_min) {
  // Interpolate a thermodynamic profile and return one momentum-level value.
  // Native sequence association exposes the contiguous profile as one column.
  Real2D redirect_interpolated_azm_col =
      redirect_interpolated_azm_2D(gr.nzm, gr.nzt, 1, gr, Real2D(azt.data(), 1, gr.nzt), zm_min);
  return to_host_space(redirect_interpolated_azm_col)(0, k);
}

// =============================================================================
static Real1D redirect_interpolated_azm_1D(const grid& gr, Real1D azt,
                                           std::optional<real_t> zm_min) {
  // Interpolate a full thermodynamic profile to momentum levels.
  Real1D redirect_interpolated_azm_1D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "redirect_interpolated_azm_1D"), gr.nzm);
  // The owning rank-one result outlives both nonowning rank-two aliases.
  zt2zm_api(gr.nzm, gr.nzt, 1, gr, Real2D(azt.data(), 1, gr.nzt),
            Real2D(redirect_interpolated_azm_1D.data(), 1, gr.nzm), zm_min);
  return redirect_interpolated_azm_1D;
}

// =============================================================================
static Real2D redirect_interpolated_azm_2D(int nzm, int nzt, int ngrdcol, const grid& gr,
                                           Real2D azt, std::optional<real_t> zm_min) {
  // Interpolate thermodynamic fields to momentum levels using the native flag.
  Real2D redirect_interpolated_azm_2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "redirect_interpolated_azm_2D"), ngrdcol,
      nzm);
  zt2zm_api(nzm, nzt, ngrdcol, gr, azt, redirect_interpolated_azm_2D, zm_min);
  return redirect_interpolated_azm_2D;
}

// =============================================================================
static real_t redirect_interpolated_azt_k(const grid& gr, Real1D azm, int k,
                                          std::optional<real_t> zt_min) {
  // Interpolate a momentum profile and return one thermodynamic-level value.
  Real2D redirect_interpolated_azt_col =
      redirect_interpolated_azt_2D(gr.nzm, gr.nzt, 1, gr, Real2D(azm.data(), 1, gr.nzm), zt_min);
  return to_host_space(redirect_interpolated_azt_col)(0, k);
}

// =============================================================================
static Real1D redirect_interpolated_azt_1D(const grid& gr, Real1D azm,
                                           std::optional<real_t> zt_min) {
  // Interpolate a full momentum profile to thermodynamic levels.
  Real1D redirect_interpolated_azt_1D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "redirect_interpolated_azt_1D"), gr.nzt);
  zm2zt_api(gr.nzm, gr.nzt, 1, gr, Real2D(azm.data(), 1, gr.nzm),
            Real2D(redirect_interpolated_azt_1D.data(), 1, gr.nzt), zt_min);
  return redirect_interpolated_azt_1D;
}

// =============================================================================
static Real2D redirect_interpolated_azt_2D(int nzm, int nzt, int ngrdcol, const grid& gr,
                                           Real2D azm, std::optional<real_t> zt_min) {
  // Interpolate momentum fields to thermodynamic levels using the native flag.
  Real2D redirect_interpolated_azt_2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "redirect_interpolated_azt_2D"), ngrdcol,
      nzt);
  zm2zt_api(nzm, nzt, ngrdcol, gr, azm, redirect_interpolated_azt_2D, zt_min);
  return redirect_interpolated_azt_2D;
}

// =============================================================================
static void linear_interpolated_azm_2D(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt,
                                       Real2D linear_interpolated_azm,
                                       std::optional<real_t> zm_min) {
  // Interpolate a thermodynamic field to momentum levels on a stretched grid.
  // Boundary extensions need not be positive; the optional minimum clips them.
  // Kernel captures use only the grid's component Views and plain optional
  // scalar values, avoiding a copy of the complete owning grid or std::optional.
  const auto weights_zt2zm = gr.weights_zt2zm;
  const int grid_dir_indx  = gr.grid_dir_indx;
  const bool l_min         = zm_min.has_value();
  const real_t minimum     = zm_min.value_or(real_t(0));

  // Interior interpolation and independent boundary work share one kernel.
  // Each work item uses the original operand order; no neighboring output is read.
  // NVHPC -O0 rounds the above-level product, then contracts the below-level
  // product and sum. Explicit fma retains that choice on CUDA; its implicit
  // contraction can otherwise fuse the other product and change one bit.
  Kokkos::parallel_for(
      "linear_interpolated_azm_2D", Kokkos::RangePolicy<>(0, ngrdcol * nzm),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        if (k > 0 && k < nzm - 1) {
          // Central momentum level between two successive thermodynamic levels.
          linear_interpolated_azm(i, k) = Kokkos::fma(weights_zt2zm(i, k, t_below), azt(i, k - 1),
                                                      weights_zt2zm(i, k, t_above) * azt(i, k));
        } else if (k == 0) {
          // Level1 is the lower ascending boundary or upper descending boundary.
          if (grid_dir_indx == 1)
            linear_interpolated_azm(i, k) = azt(i, 0);
          else
            // Linear extension from the first two thermodynamic levels.
            linear_interpolated_azm(i, k) = Kokkos::fma(weights_zt2zm(i, k, t_below), azt(i, 0),
                                                        weights_zt2zm(i, k, t_above) * azt(i, 1));
        } else {
          // Final momentum level is the upper ascending/lower descending boundary.
          if (grid_dir_indx == 1)
            // Linear extension from the top two thermodynamic levels.
            linear_interpolated_azm(i, k) =
                Kokkos::fma(weights_zt2zm(i, k, t_below), azt(i, nzt - 2),
                            weights_zt2zm(i, k, t_above) * azt(i, nzt - 1));
          else
            linear_interpolated_azm(i, k) = azt(i, nzt - 1);
        }
        // Native MAX selects the minimum on ties, including opposite signed zeros.
        if (l_min)
          linear_interpolated_azm(i, k) =
              linear_interpolated_azm(i, k) > minimum ? linear_interpolated_azm(i, k) : minimum;
      });
}

// =============================================================================
Real2D zt2zm2zt(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt,
                std::optional<real_t> zt_min) {
  // Smooth a thermodynamic field by interpolating to momentum levels and back.
  // Apply the optional minimum only on the second interpolation.

  // Interpolate azt to momentum levels.
  Real2D azt_zm = zt2zm_api(nzm, nzt, ngrdcol, gr, azt);

  // Interpolate back to thermodynamic levels.
  Real2D zt2zm2zt = zm2zt_api(nzm, nzt, ngrdcol, gr, azt_zm, zt_min);
  return zt2zm2zt;
}

// =============================================================================
Real2D zm2zt2zm(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm,
                std::optional<real_t> zm_min) {
  // Smooth a momentum field by interpolating to thermodynamic levels and back.

  // Interpolate azm to thermodynamic levels.
  Real2D azm_zt = zm2zt_api(nzm, nzt, ngrdcol, gr, azm);

  // Interpolate back to momentum levels, applying the optional minimum here.
  Real2D zm2zt2zm = zt2zm_api(nzm, nzt, ngrdcol, gr, azm_zt, zm_min);
  return zm2zt2zm;
}

// =============================================================================
static void cubic_interpolated_azm_2D(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt,
                                      Real2D cubic_interpolated_azm_2D) {
  // Interpolate a thermodynamic field to all momentum levels using Steffen's
  // monotonic cubic spline, implemented by Tak Yamaguchi.
  const auto zm = gr.zm, zt = gr.zt;

  Kokkos::parallel_for(
      "cubic_interpolated_azm_2D", Kokkos::RangePolicy<>(0, ngrdcol * nzm),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        int km1, k00, kp1, kp2;

        // k levels follow Tak's find_indices (dschanen,24 Oct2011), shifted -1.
        if (k == nzm - 2) {
          km1 = nzt - 3;
          kp1 = nzt - 1;
          kp2 = nzt - 1;
          k00 = nzt - 2;
        } else if (k == nzm - 1) { // Linear extension to the final momentum level.
          km1 = nzt - 1;
          kp1 = nzt - 1;
          kp2 = nzt - 1;
          k00 = nzt - 2;
        } else if (k == 0) { // Linear extension to the first momentum level.
          km1 = 1;
          kp1 = 1;
          kp2 = 1;
          k00 = 0;
        } else if (k == 1) {
          kp1 = 1;
          kp2 = 2;
          k00 = 0;
          km1 = 0;
        } else {
          km1 = k - 2;
          kp1 = k;
          kp2 = k + 1;
          k00 = k - 1;
        }

        // Do the actual interpolation using a cubic monotonic spline.
        cubic_interpolated_azm_2D(i, k) =
            mono_cubic_interp(zm(i, k), km1, k00, kp1, kp2, zt(i, km1), zt(i, k00), zt(i, kp1),
                              zt(i, kp2), azt(i, km1), azt(i, k00), azt(i, kp1), azt(i, kp2));
      });
}

// =============================================================================
// =============================================================================
static void calc_zt2zm_weights(int nzm, int nzt, int ngrdcol, // In
                               grid& gr) {                    // InOut
  // Description:
  // Function used to help in an interpolation of a variable (var_zt) located
  // on the thermodynamic grid levels (azt) to the momentum grid levels (azm).
  // This function computes a weighting factor for both the upper thermodynamic
  // level (k) and the lower thermodynamic level (k-1) applied to the central
  // momentum level (k).  For the boundary momentum grid level k=1, a weighting
  // factor for both the thermodynamic level at k=1 and the thermodynamic level
  // at k=2 are calculated based on the use of a linear extension. Likewise,
  // for the boundary momentum grid level k=nzm, a weighting factor for both
  // the thermodynamic level at k=nzt and the thermodynamic level at k=nzt-1
  // are calculated based on the use of a linear extension.  This function
  // outputs the weighting factors at a single momentum grid level (k).  The
  // formulation used is compatible with a stretched (unevenly-spaced) grid.
  // The weights are defined as follows:
  //
  // ---var_zt(k)--------------------------------------------- t(k)
  //                       azt_weight(t_above) = factor
  // ===========var_zt(interp)================================ m(k)
  //                       azt_weight(t_below) = 1 - factor
  // ---var_zt(k-1)------------------------------------------- t(k-1)
  //
  // The vertical indices t(k), m(k), and t(k-1) correspond with altitudes
  // zt(k), zm(k), and zt(k-1), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // NOTE: For an ascending grid, the grid diagram shown above has increasing
  // indices with increasing altitude, such that zt(k) > zt(k-1). For a
  // descending grid, the grid diagram shown above has increasing indices with
  // decreasing altitude, such that zt(k) < zt(k-1).
  //
  // For all momentum levels 1 < k < nzm:
  //
  // The formula for a linear interpolation is given by:
  //
  // var_zt( interp to zm(k) )
  // = [ ( var_zt(k) - var_zt(k-1) ) / ( zt(k) - zt(k-1) ) ]
  //     * ( zm(k) - zt(k-1) ) + var_zt(k-1);
  //
  // which can be rewritten as:
  //
  // var_zt( interp to zm(k) )
  // = [ ( zm(k) - zt(k-1) ) / ( zt(k) - zt(k-1) ) ]
  //     * ( var_zt(k) - var_zt(k-1) ) + var_zt(k-1).
  //
  // Furthermore, the formula can be rewritten as:
  //
  // var_zt( interp to zm(k) )
  // = factor * var_zt(k) + ( 1 - factor ) * var_zt(k-1);
  //
  // where:
  //
  // factor = ( zm(k) - zt(k-1) ) / ( zt(k) - zt(k-1) ).
  //
  // One of the important uses of this function is in situations where the
  // variable to be interpolated is being treated IMPLICITLY in an equation.
  // Usually, the variable to be interpolated is involved in a derivative (such
  // as d(var_zt)/dz in the diagram below).  For the term of the equation
  // containing the derivative, grid weights are needed for two interpolations,
  // rather than just one interpolation.  Thus, four grid weights (labeled
  // A(k), B(k), C(k), and D(k) in the diagram below) are needed.
  //
  // ---var_zt(k+1)------------------------------------------- t(k+1)
  //                                       A(k)
  // ===========var_zt(interp)================================ m(k+1)
  //                                       B(k) = 1 - A(k)
  // ---var_zt(k)-----------d(var_zt)/dz---------------------- t(k)
  //                                       C(k)
  // ===========var_zt(interp)================================ m(k)
  //                                       D(k) = 1 - C(k)
  // ---var_zt(k-1)------------------------------------------- t(k-1)
  //
  // The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
  // with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
  // The letter "t" is used for thermodynamic levels and the letter "m" is used
  // for momentum levels.
  //
  // The grid weights, indexed around the central thermodynamic level (k), are
  // defined as follows:
  //
  // A(k) = ( zm(k+1) - zt(k) ) / ( zt(k+1) - zt(k) );
  //
  // which is the same as "factor" for the interpolation to momentum
  // level (k+1).  In the code, this interpolation is referenced as
  // gr%weights_zt2zm(t_above,mkp1), which can be read as "grid weight in a
  // zt2zm interpolation of the thermodynamic level above momentum level (k+1)
  // (applied to momentum level (k+1))".
  //
  // B(k) = 1 - [ ( zm(k+1) - zt(k) ) / ( zt(k+1) - zt(k) ) ]
  //      = 1 - A(k);
  //
  // which is the same as "1 - factor" for the interpolation to momentum
  // level (k+1).  This equation can also be rewritten as:
  //
  // B(k) = ( zt(k+1) - zm(k+1) ) / ( zt(k+1) - zt(k) )
  //
  // In the code, this interpolation is referenced as
  // gr%weights_zt2zm(t_below,mkp1), which can be read as "grid weight in a
  // zt2zm interpolation of the thermodynamic level below momentum level (k+1)
  // (applied to momentum level (k+1))".
  //
  // C(k) = ( zm(k) - zt(k-1) ) / ( zt(k) - zt(k-1) );
  //
  // which is the same as "factor" for the interpolation to momentum
  // level (k).  In the code, this interpolation is referenced as
  // gr%weights_zt2zm(t_above,mk), which can be read as "grid weight in a
  // zt2zm interpolation of the thermodynamic level above momentum level (k)
  // (applied to momentum level (k))".
  //
  // D(k) = 1 - [ ( zm(k) - zt(k-1) ) / ( zt(k) - zt(k-1) ) ]
  //      = 1 - C(k);
  //
  // which is the same as "1 - factor" for the interpolation to momentum
  // level (k).  This equation can also be rewritten as:
  //
  // D(k) = ( zt(k) - zm(k) ) / ( zt(k) - zt(k-1) )
  //
  // In the code, this interpolation is referenced as
  // gr%weights_zt2zm(t_below,mk), which can be read as "grid weight in a
  // zt2zm interpolation of the thermodynamic level below momentum level (k)
  // (applied to momentum level (k))".
  //
  // Additionally, as long as the central thermodynamic level (k) in the above
  // scenario is not a boundary thermodynamic level (k /= nzt and k /= 1), the
  // four weighting factors have the following relationships:  A(k) = C(k+1)
  // and B(k) = D(k+1).
  //
  //
  // Special condition for grid level, k = 1:
  //
  // The boundary momentum grid level k = 1 is outside the thermodynamic grid
  // level k = 1.  Thus, a linear extension is used at this level.
  //
  // The formula for a linear extension to momentum level k = 1 is given by:
  //
  // var_zt( extend to zm(k) )
  // = [ ( var_zt(k+1) - var_zt(k) ) / ( zt(k+1) - zt(k) ) ]
  //     * ( zm(k) - zt(k) ) + var_zt(k);
  //
  // which can be rewritten as:
  //
  // var_zt( extend to zm(k) )
  // = [ ( zm(k) - zt(k) ) / ( zt(k+1) - zt(k) ) ]
  //     * ( var_zt(k+1) - var_zt(k) ) + var_zt(k).
  //
  // Furthermore, the formula can be rewritten as:
  //
  // var_zt( extend to zm(k) )
  // = factor * var_zt(k+1) + ( 1 - factor ) * var_zt(k);
  //
  // where:
  //
  // factor = ( zm(k) - zt(k) ) / ( zt(k+1) - zt(k) ).
  //
  // Finally, since here k = 1:
  //
  // factor = ( zm(1) - zt(1) ) / ( zt(2) - zt(1) ).
  //
  // Due to the fact that a linear extension is being used, the value of factor
  // will be less than 0.  The weight of thermodynamic level k = 2 on momentum
  // level k = 1 equals the value of factor.  The weight of thermodynamic level
  // k = 1 on momentum level k = 1 equals 1 - factor, which is greater than 1.
  // However, the sum of the weights equals 1.
  //
  //
  // Special condition for grid level, k = nzm:
  //
  // The boundary momentum grid level k = nzm is outside the thermodynamic grid
  // level k = nzt.  Thus, a linear extension is used at this level.
  //
  // The formula for a linear extension to momentum level k = nzm
  // is given by:
  //
  // var_zt( extend to zm(k) )
  // = [ ( var_zt(k-1) - var_zt(k-2) ) / ( zt(k-1) - zt(k-2) ) ]
  //     * ( zm(k) - zt(k-2) ) + var_zt(k-2);
  //
  // which can be rewritten as:
  //
  // var_zt( extend to zm(k) )
  // = [ ( zm(k) - zt(k-2) ) / ( zt(k-1) - zt(k-2) ) ]
  //     * ( var_zt(k-1) - var_zt(k-2) ) + var_zt(k-2).
  //
  // Furthermore, the formula can be rewritten as:
  //
  // var_zt( extend to zm(k) )
  // = factor * var_zt(k-1) + ( 1 - factor ) * var_zt(k-2);
  //
  // where:
  //
  // factor = ( zm(k) - zt(k-2) ) / ( zt(k-1) - zt(k-2) ).
  //
  // Finally, since here k = nzm and nzt = nzm - 1:
  //
  // factor = ( zm(nzm) - zt(nzt-1) ) / ( zt(nzt) - zt(nzt-1) ).
  //
  // Due to the fact that a linear extension is being used, the value of factor
  // will be greater than 1.  The weight of thermodynamic level k = nzt on
  // momentum level k = nzm equals the value of factor.  The weight of
  // thermodynamic level k = nzt-1 on momentum level k = nzm equals
  // 1 - factor, which is less than 0.  However, the sum of the two weights
  // equals 1.
  //
  //
  // Brian Griffin; September 12, 2008.
  // Updated for removal of "ghost" thermodynamic level that was below the
  // surface; July 17, 2023.
  // Updated to generalize the grid for both ascending and descending grid
  // directions; September 24, 2024.
  //
  //-----------------------------------------------------------------------
  // nzm, nzt: Momentum and thermodynamic level counts [#].
  // ngrdcol: Number of columns [#].
  // gr: Heights [m] and interpolation weights [-] in the existing grid.
  const auto zm = gr.zm;
  const auto zt = gr.zt;
  const auto weights_zt2zm = gr.weights_zt2zm;
  // Interior interpolation and the two independent linear extensions share
  // one launch. Native requires at least two thermodynamic levels.
  Kokkos::parallel_for(
      "calc_zt2zm_weights", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        if (k > 0 && k < nzm - 1) {
          // Weight of upper thermodynamic level on momentum level.
          weights_zt2zm(i, k, t_above) = (zm(i, k) - zt(i, k - 1)) / (zt(i, k) - zt(i, k - 1));
          // Weight of lower thermodynamic level on momentum level.
          weights_zt2zm(i, k, t_below) = (zt(i, k) - zm(i, k)) / (zt(i, k) - zt(i, k - 1));
        } else if (k == 0) {
          // Storage level 1 is outside thermodynamic level 1: extend from
          // thermodynamic levels 1 and 2. This is the lower boundary ascending
          // and the upper boundary descending; factor is less than zero.
          weights_zt2zm(i, 0, t_above) = (zm(i, 0) - zt(i, 0)) / (zt(i, 1) - zt(i, 0));
          weights_zt2zm(i, 0, t_below) = (zt(i, 1) - zm(i, 0)) / (zt(i, 1) - zt(i, 0));
        } else {
          // Storage level nzm is outside thermodynamic level nzt: extend from
          // levels nzt and nzt-1. This is the upper boundary ascending and the
          // lower boundary descending; factor is greater than one.
          weights_zt2zm(i, nzm - 1, t_above) =
              (zm(i, nzm - 1) - zt(i, nzt - 2)) / (zt(i, nzt - 1) - zt(i, nzt - 2));
          weights_zt2zm(i, nzm - 1, t_below) =
              (zt(i, nzt - 1) - zm(i, nzm - 1)) / (zt(i, nzt - 1) - zt(i, nzt - 2));
        }
      });
}

static void linear_interpolated_azt_2D(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm,
                                       Real2D linear_interpolated_azt,
                                       std::optional<real_t> zt_min) {
  // Interpolate a momentum field to all thermodynamic levels on a stretched grid.
  const auto weights_zm2zt = gr.weights_zm2zt;
  const bool l_min         = zt_min.has_value();
  const real_t minimum     = zt_min.value_or(real_t(0));

  // Every thermodynamic level lies between two successive momentum levels.
  // Fuse the independent native minimum pass into this same cell's update.
  // Retain NVHPC's rounded above product plus fused below product, as above.
  Kokkos::parallel_for(
      "linear_interpolated_azt_2D", Kokkos::RangePolicy<>(0, ngrdcol * nzt),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        linear_interpolated_azt(i, k) = Kokkos::fma(weights_zm2zt(i, k, m_below), azm(i, k),
                                                    weights_zm2zt(i, k, m_above) * azm(i, k + 1));
        if (l_min)
          linear_interpolated_azt(i, k) =
              linear_interpolated_azt(i, k) > minimum ? linear_interpolated_azt(i, k) : minimum;
      });
}

// =============================================================================
static void cubic_interpolated_azt_2D(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm,
                                      Real2D cubic_interpolated_azt_2D) {
  // Interpolate a momentum field to all thermodynamic levels with Steffen's
  // monotonic cubic spline, implemented by Tak Yamaguchi.
  const auto zm = gr.zm, zt = gr.zt;

  Kokkos::parallel_for(
      "cubic_interpolated_azt_2D", Kokkos::RangePolicy<>(0, ngrdcol * nzt),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        int km1, k00, kp1, kp2;

        // k levels follow Tak's find_indices (dschanen,24 Oct2011), shifted -1.
        if (k == nzt - 1) {
          km1 = nzm - 3;
          kp1 = nzm - 1;
          kp2 = nzm - 1;
          k00 = nzm - 2;
        } else if (k == 0) {
          km1 = 0;
          kp1 = 1;
          kp2 = 2;
          k00 = 0;
        } else {
          km1 = k - 1;
          kp1 = k + 1;
          kp2 = k + 2;
          k00 = k;
        }

        // Do the actual interpolation using a cubic monotonic spline.
        cubic_interpolated_azt_2D(i, k) =
            mono_cubic_interp(zt(i, k), km1, k00, kp1, kp2, zm(i, km1), zm(i, k00), zm(i, kp1),
                              zm(i, kp2), azm(i, km1), azm(i, k00), azm(i, kp1), azm(i, kp2));
      });
}

// =============================================================================
// =============================================================================
static void calc_zm2zt_weights(int nzt, int ngrdcol, // In
                               grid& gr) {           // InOut
  // Description:
  // Function used to help in an interpolation of a variable (var_zm) located
  // on the momentum grid levels (azm) to the thermodynamic grid levels (azt).
  // This function computes a weighting factor for both the upper momentum
  // level (k+1) and the lower momentum level (k) applied to the central
  // thermodynamic level (k).  This function outputs the weighting factors at
  // a single thermodynamic grid level (k).   The formulation used is
  // compatible with a stretched (unevenly-spaced) grid.  The weights are
  // defined as follows:
  //
  // ===var_zm(k+1)=========================================== m(k+1)
  //                       azm_weight(m_above) = factor
  // -----------var_zm(interp)-------------------------------- t(k)
  //                       azm_weight(m_below) = 1 - factor
  // ===var_zm(k)============================================= m(k)
  //
  // The vertical indices m(k+1), t(k), and m(k) correspond with altitudes
  // zm(k+1), zt(k), and zm(k), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // NOTE: For an ascending grid, the grid diagram shown above has increasing
  // indices with increasing altitude, such that zm(k+1) > zm(k). For a
  // descending grid, the grid diagram shown above has increasing indices with
  // decreasing altitude, such that zm(k+1) < zm(k).
  //
  // The formula for a linear interpolation is given by:
  //
  // var_zm( interp to zt(k) )
  // = [ ( var_zm(k+1) - var_zm(k) ) / ( zm(k+1) - zm(k) ) ]
  //     * ( zt(k) - zm(k) ) + var_zm(k);
  //
  // which can be rewritten as:
  //
  // var_zm( interp to zt(k) )
  // = [ ( zt(k) - zm(k) ) / ( zm(k+1) - zm(k) ) ]
  //     * ( var_zm(k+1) - var_zm(k) ) + var_zm(k).
  //
  // Furthermore, the formula can be rewritten as:
  //
  // var_zm( interp to zt(k) )
  // = factor * var_zm(k+1) + ( 1 - factor ) * var_zm(k);
  //
  // where:
  //
  // factor = ( zt(k) - zm(k) ) / ( zm(k+1) - zm(k) ).
  //
  // One of the important uses of this function is in situations where the
  // variable to be interpolated is being treated IMPLICITLY in an equation.
  // Usually, the variable to be interpolated is involved in a derivative (such
  // as d(var_zm)/dz in the diagram below).  For the term of the equation
  // containing the derivative, grid weights are needed for two interpolations,
  // rather than just one interpolation.  Thus, four grid weights (labeled
  // A(k), B(k), C(k), and D(k) in the diagram below) are needed.
  //
  // ===var_zm(k+1)=========================================== m(k+1)
  //                                       A(k)
  // -----------var_zm(interp)-------------------------------- t(k)
  //                                       B(k) = 1 - A(k)
  // ===var_zm(k)===========d(var_zm)/dz====================== m(k)
  //                                       C(k)
  // -----------var_zm(interp)-------------------------------- t(k-1)
  //                                       D(k) = 1 - C(k)
  // ===var_zm(k-1)=========================================== m(k-1)
  //
  // The vertical indices m(k+1), t(k), m(k), t(k-1), and m(k-1) correspond
  // with altitudes zm(k+1), zt(k), zm(k), zt(k-1), and zm(k-1), respectively.
  // The letter "t" is used for thermodynamic levels and the letter "m" is used
  // for momentum levels.
  //
  // The grid weights, indexed around the central momentum level (k), are
  // defined as follows:
  //
  // A(k) = ( zt(k) - zm(k) ) / ( zm(k+1) - zm(k) );
  //
  // which is the same as "factor" for the interpolation to thermodynamic
  // level (k).  In the code, this interpolation is referenced as
  // gr%weights_zm2zt(m_above,tk), which can be read as "grid weight in a
  // zm2zt interpolation of the momentum level above thermodynamic
  // level (k) (applied to thermodynamic level (k))".
  //
  // B(k) = 1 - [ ( zt(k) - zm(k) ) / ( zm(k+1) - zm(k) ) ]
  //      = 1 - A(k);
  //
  // which is the same as "1 - factor" for the interpolation to thermodynamic
  // level (k).  This equation can also be rewritten as:
  //
  // B(k) = ( zm(k+1) - zt(k) ) / ( zm(k+1) - zm(k) ).
  //
  // In the code, this interpolation is referenced as
  // gr%weights_zm2zt(m_below,tk), which can be read as "grid weight in a
  // zm2zt interpolation of the momentum level below thermodynamic
  // level (k) (applied to thermodynamic level (k))".
  //
  // C(k) = ( zt(k-1) - zm(k-1) ) / ( zm(k) - zm(k-1) );
  //
  // which is the same as "factor" for the interpolation to thermodynamic
  // level (k-1).  In the code, this interpolation is referenced as
  // gr%weights_zm2zt(m_above,tkm1), which can be read as "grid weight in a
  // zm2zt interpolation of the momentum level above thermodynamic level (k-1)
  // (applied to thermodynamic level (k-1))".
  //
  // D(k) = 1 - [ ( zt(k-1) - zm(k-1) ) / ( zm(k) - zm(k-1) ) ]
  //      = 1 - C(k);
  //
  // which is the same as "1 - factor" for the interpolation to thermodynamic
  // level (k-1).  The equation can also be rewritten as:
  //
  // D(k) = ( zm(k) - zt(k-1) ) / ( zm(k) - zm(k-1) )
  //
  // In the code, this interpolation is referenced as
  // gr%weights_zm2zt(m_below,tkm1), which can be read as "grid weight in a
  // zm2zt interpolation of the momentum level below thermodynamic level (k-1)
  // (applied to thermodynamic level (k-1))".
  //
  // Additionally, as long as the central momentum level (k) in the above
  // scenario is not the lowermost momentum level or the uppermost momentum
  // level (k /= 1 and k /= gr%nz), the four weighting factors have the
  // following relationships:  A(k) = C(k+1) and B(k) = D(k+1).
  //
  //
  // Brian Griffin; September 12, 2008.
  // Updated for removal of "ghost" thermodynamic level that was below the
  // surface; July 18, 2023.
  // Updated to generalize the grid for both ascending and descending grid
  // directions; September 24, 2024.
  //
  //-----------------------------------------------------------------------
  // nzt: Number of thermodynamic levels [#].
  // ngrdcol: Number of columns [#].
  // gr: Heights [m] and interpolation weights [-] in the existing grid.
  const auto zm = gr.zm;
  const auto zt = gr.zt;
  const auto weights_zm2zt = gr.weights_zm2zt;
  // The thermodynamic levels are always between two momentum levels.
  // Linear interpolation is used.
  Kokkos::parallel_for(
      "calc_zm2zt_weights", ngrdcol * nzt, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        // Weight of upper momentum level on thermodynamic level.
        weights_zm2zt(i, k, m_above) = (zt(i, k) - zm(i, k)) / (zm(i, k + 1) - zm(i, k));
        // Weight of lower momentum level on thermodynamic level.
        weights_zm2zt(i, k, m_below) = (zm(i, k + 1) - zt(i, k)) / (zm(i, k + 1) - zm(i, k));
      });
}

static void gradzm_2D(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm, Real2D gradzm_2D) {
  // Differentiate a momentum-level field across thermodynamic levels.
  // The returned field has input units divided by metres.
  const auto invrs_dzt = gr.invrs_dzt;

  // Loop over all thermodynamic grid levels.
  Kokkos::parallel_for(
      "gradzm_2D", Kokkos::RangePolicy<>(0, ngrdcol * nzt), KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        gradzm_2D(i, k) = (azm(i, k + 1) - azm(i, k)) * invrs_dzt(i, k);
      });
}

// =============================================================================
static Real1D gradzm_1D(const grid& gr, Real1D azm) {
  // Return a momentum profile's vertical derivative on thermodynamic levels.
  Real1D gradzm_1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "gradzm_1D"), gr.nzt);
  gradzm_2D(gr.nzm, gr.nzt, 1, gr, Real2D(azm.data(), 1, gr.nzm),
            Real2D(gradzm_1D.data(), 1, gr.nzt));
  return gradzm_1D;
}

// =============================================================================
static void gradzt_2D(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt, Real2D gradzt_2D) {
  // Differentiate a thermodynamic-level field across momentum levels.
  const auto invrs_dzm = gr.invrs_dzm;

  // Interior momentum levels and boundary extensions share a single kernel.
  // Native boundaries copy the neighboring derivative. Repeating that same
  // expression here avoids a dependency and retains its exact operand order.
  Kokkos::parallel_for(
      "gradzt_2D", Kokkos::RangePolicy<>(0, ngrdcol * nzm), KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        const int interior = k == 0 ? 1 : k == nzm - 1 ? nzm - 2 : k;
        gradzt_2D(i, k)    = (azt(i, interior) - azt(i, interior - 1)) * invrs_dzm(i, interior);
      });
}

// =============================================================================
static Real1D gradzt_1D(const grid& gr, Real1D azt) {
  // Return a thermodynamic profile's derivative on momentum levels, extending
  // each boundary from its adjacent interior derivative like the 2D interface.
  Real1D gradzt_1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "gradzt_1D"), gr.nzm);
  gradzt_2D(gr.nzm, gr.nzt, 1, gr, Real2D(azt.data(), 1, gr.nzt),
            Real2D(gradzt_1D.data(), 1, gr.nzm));
  return gradzt_1D;
}

// =============================================================================
Real1D flip(Real1D x, int xdim) {
  // Reverse a profile for packages such as BUGSrad that store heights in the
  // opposite order. tmp owns the reversed result, preserving the input.
  Real1D tmp("flip::tmp", xdim);
  Kokkos::parallel_for(
      "flip", Kokkos::RangePolicy<>(0, xdim),
      KOKKOS_LAMBDA(const int indx) { tmp(indx) = x(xdim - 1 - indx); });
  return tmp;
}

// Native generic dispatch is explicit in C++; no additional scientific routine.
Real2D zt2zm_api(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt,
                 std::optional<real_t> zm_min) {
  // Select the native column/profile overload.
  return redirect_interpolated_azm_2D(nzm, nzt, ngrdcol, gr, azt, zm_min);
}
void zt2zm_api(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt, // In
               Real2D azm,                                                // Out
               std::optional<real_t> zm_min) {                            // Optional In
  // Fill caller-owned momentum levels with the native interpolation/floor rules.
  if (overlapping_fields(azt, azm)) {
    Kokkos::deep_copy(azm, zt2zm_api(nzm, nzt, ngrdcol, gr, azt, zm_min));
    return;
  }
  if (l_quintic_poly_interp && !l_cubic_interp)
    throw std::runtime_error("Error: Model flag l_quintic_poly_interp should not be true if "
                             "l_cubic_interp is false.");
  if (l_cubic_interp && nzm >= 4)
    cubic_interpolated_azm_2D(nzm, nzt, ngrdcol, gr, azt, azm);
  else
    linear_interpolated_azm_2D(nzm, nzt, ngrdcol, gr, azt, azm, zm_min);
}
Real1D zt2zm_api(const grid& gr, Real1D azt, std::optional<real_t> zm_min) {
  // Select the native profile overload.
  return redirect_interpolated_azm_1D(gr, azt, zm_min);
}
Real1D zt2zm_api(const grid& gr, Real1D azt, real_t zm_min) {
  // A real floor selects the profile overload; an integer selects a level.
  return redirect_interpolated_azm_1D(gr, azt, zm_min);
}
real_t zt2zm_api(const grid& gr, Real1D azt, int k, std::optional<real_t> zm_min) {
  // Select one zero-based level through the native scalar overload.
  return redirect_interpolated_azm_k(gr, azt, k, zm_min);
}
Real2D zm2zt_api(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm,
                 std::optional<real_t> zt_min) {
  // Select the native column/profile overload.
  return redirect_interpolated_azt_2D(nzm, nzt, ngrdcol, gr, azm, zt_min);
}
void zm2zt_api(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm, // In
               Real2D azt,                                                // Out
               std::optional<real_t> zt_min) {                            // Optional In
  // Fill caller-owned thermodynamic levels; cubic paths retain native floor behavior.
  if (overlapping_fields(azm, azt)) {
    Kokkos::deep_copy(azt, zm2zt_api(nzm, nzt, ngrdcol, gr, azm, zt_min));
    return;
  }
  if (l_quintic_poly_interp && !l_cubic_interp)
    throw std::runtime_error("Error: Model flag l_quintic_poly_interp should not be true if "
                             "l_cubic_interp is false.");
  if (l_cubic_interp && nzm >= 2)
    cubic_interpolated_azt_2D(nzm, nzt, ngrdcol, gr, azm, azt);
  else
    linear_interpolated_azt_2D(nzm, nzt, ngrdcol, gr, azm, azt, zt_min);
}
Real1D zm2zt_api(const grid& gr, Real1D azm, std::optional<real_t> zt_min) {
  // Select the native profile overload.
  return redirect_interpolated_azt_1D(gr, azm, zt_min);
}
Real1D zm2zt_api(const grid& gr, Real1D azm, real_t zt_min) {
  // Preserve the native real-floor/integer-level generic distinction.
  return redirect_interpolated_azt_1D(gr, azm, zt_min);
}
real_t zm2zt_api(const grid& gr, Real1D azm, int k, std::optional<real_t> zt_min) {
  // Select one zero-based level through the native scalar overload.
  return redirect_interpolated_azt_k(gr, azm, k, zt_min);
}
Real2D ddzm(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm) {
  // Select the native column/profile derivative overload.
  Real2D derivative(Kokkos::view_alloc(Kokkos::WithoutInitializing, "gradzm_2D"), ngrdcol, nzt);
  ddzm(nzm, nzt, ngrdcol, gr, azm, derivative);
  return derivative;
}
void ddzm(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm, // In
          Real2D derivative) {                                       // Out
  // Differentiate momentum levels directly into existing thermodynamic storage.
  if (overlapping_fields(azm, derivative)) {
    Kokkos::deep_copy(derivative, ddzm(nzm, nzt, ngrdcol, gr, azm));
    return;
  }
  gradzm_2D(nzm, nzt, ngrdcol, gr, azm, derivative);
}
Real1D ddzm(const grid& gr, Real1D azm) {
  // Select the native profile derivative overload.
  return gradzm_1D(gr, azm);
}
Real2D ddzt(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt) {
  // Select the native column/profile derivative overload.
  Real2D derivative(Kokkos::view_alloc(Kokkos::WithoutInitializing, "gradzt_2D"), ngrdcol, nzm);
  ddzt(nzm, nzt, ngrdcol, gr, azt, derivative);
  return derivative;
}
void ddzt(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt, // In
          Real2D derivative) {                                       // Out
  // Differentiate thermodynamic levels directly into existing momentum storage.
  if (overlapping_fields(azt, derivative)) {
    Kokkos::deep_copy(derivative, ddzt(nzm, nzt, ngrdcol, gr, azt));
    return;
  }
  gradzt_2D(nzm, nzt, ngrdcol, gr, azt, derivative);
}
Real1D ddzt(const grid& gr, Real1D azt) {
  // Select the native profile derivative overload.
  return gradzt_1D(gr, azt);
}

// Caller-frame overloads delegate every numerical operation to the existing APIs.
Real2D zt2zm_api(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
                 Real2D azt, std::optional<real_t> zm_min) {
  // Reserve the result in the caller frame; the numerical operator is shared.
  auto azm = workspace.view<Real2D>(Kokkos::WithoutInitializing, "redirect_interpolated_azm_2D",
                                    ngrdcol, nzm);
  zt2zm_api(nzm, nzt, ngrdcol, gr, azt, azm, zm_min);
  return azm;
}

void zt2zm_api(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
               Real2D azt, Real2D azm, std::optional<real_t> zm_min) {
  // Keep source values intact when an explicit output overlaps its input.
  if (overlapping_fields(azt, azm)) {
    auto temporary =
        workspace.view<Real2D>(Kokkos::WithoutInitializing, "zt2zm_api::overlap", ngrdcol, nzm);
    zt2zm_api(nzm, nzt, ngrdcol, gr, azt, temporary, zm_min);
    Kokkos::deep_copy(azm, temporary);
    return;
  }
  zt2zm_api(nzm, nzt, ngrdcol, gr, azt, azm, zm_min);
}
Real2D zm2zt_api(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
                 Real2D azm, std::optional<real_t> zt_min) {
  // Reserve the result in the caller frame; the numerical operator is shared.
  auto azt = workspace.view<Real2D>(Kokkos::WithoutInitializing, "redirect_interpolated_azt_2D",
                                    ngrdcol, nzt);
  zm2zt_api(nzm, nzt, ngrdcol, gr, azm, azt, zt_min);
  return azt;
}

void zm2zt_api(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
               Real2D azm, Real2D azt, std::optional<real_t> zt_min) {
  // Keep source values intact when an explicit output overlaps its input.
  if (overlapping_fields(azm, azt)) {
    auto temporary =
        workspace.view<Real2D>(Kokkos::WithoutInitializing, "zm2zt_api::overlap", ngrdcol, nzt);
    zm2zt_api(nzm, nzt, ngrdcol, gr, azm, temporary, zt_min);
    Kokkos::deep_copy(azt, temporary);
    return;
  }
  zm2zt_api(nzm, nzt, ngrdcol, gr, azm, azt, zt_min);
}
Real2D ddzm(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
            Real2D azm) {
  // Reserve the result in the caller frame; the numerical operator is shared.
  auto derivative = workspace.view<Real2D>(Kokkos::WithoutInitializing, "gradzm_2D", ngrdcol, nzt);
  ddzm(nzm, nzt, ngrdcol, gr, azm, derivative);
  return derivative;
}

void ddzm(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm,
          Real2D derivative) {
  // Keep source values intact when an explicit output overlaps its input.
  if (overlapping_fields(azm, derivative)) {
    auto temporary =
        workspace.view<Real2D>(Kokkos::WithoutInitializing, "ddzm::overlap", ngrdcol, nzt);
    ddzm(nzm, nzt, ngrdcol, gr, azm, temporary);
    Kokkos::deep_copy(derivative, temporary);
    return;
  }
  ddzm(nzm, nzt, ngrdcol, gr, azm, derivative);
}
Real2D ddzt(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
            Real2D azt) {
  // Reserve the result in the caller frame; the numerical operator is shared.
  auto derivative = workspace.view<Real2D>(Kokkos::WithoutInitializing, "gradzt_2D", ngrdcol, nzm);
  ddzt(nzm, nzt, ngrdcol, gr, azt, derivative);
  return derivative;
}

void ddzt(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt,
          Real2D derivative) {
  // Keep source values intact when an explicit output overlaps its input.
  if (overlapping_fields(azt, derivative)) {
    auto temporary =
        workspace.view<Real2D>(Kokkos::WithoutInitializing, "ddzt::overlap", ngrdcol, nzm);
    ddzt(nzm, nzt, ngrdcol, gr, azt, temporary);
    Kokkos::deep_copy(derivative, temporary);
    return;
  }
  ddzt(nzm, nzt, ngrdcol, gr, azt, derivative);
}
Real2D zt2zm2zt(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
                Real2D azt, std::optional<real_t> zt_min) {
  // Keep both interpolation results in the active caller frame.
  auto smoothed =
      workspace.view<Real2D>(Kokkos::WithoutInitializing, "zt2zm2zt::result", ngrdcol, nzt);
  zt2zm2zt(workspace, nzm, nzt, ngrdcol, gr, azt, smoothed, zt_min);
  return smoothed;
}

void zt2zm2zt(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
              Real2D azt, Real2D smoothed, std::optional<real_t> zt_min) {
  // Complete the first interpolation before the second writes any overlapping input.
  // Preserve the Native floor on the second interpolation only.
  auto azt_zm =
      workspace.view<Real2D>(Kokkos::WithoutInitializing, "zt2zm2zt::intermediate", ngrdcol, nzm);
  zt2zm_api(nzm, nzt, ngrdcol, gr, azt, azt_zm);
  zm2zt_api(nzm, nzt, ngrdcol, gr, azt_zm, smoothed, zt_min);
}
Real2D zm2zt2zm(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
                Real2D azm, std::optional<real_t> zm_min) {
  // Keep both interpolation results in the active caller frame.
  auto smoothed =
      workspace.view<Real2D>(Kokkos::WithoutInitializing, "zm2zt2zm::result", ngrdcol, nzm);
  zm2zt2zm(workspace, nzm, nzt, ngrdcol, gr, azm, smoothed, zm_min);
  return smoothed;
}

void zm2zt2zm(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
              Real2D azm, Real2D smoothed, std::optional<real_t> zm_min) {
  // Complete the first interpolation before the second writes any overlapping input.
  // Preserve the Native floor on the second interpolation only.
  auto azm_zt =
      workspace.view<Real2D>(Kokkos::WithoutInitializing, "zm2zt2zm::intermediate", ngrdcol, nzt);
  zm2zt_api(nzm, nzt, ngrdcol, gr, azm, azm_zt);
  zt2zm_api(nzm, nzt, ngrdcol, gr, azm_zt, smoothed, zm_min);
}

} // namespace clubb_kokkos
