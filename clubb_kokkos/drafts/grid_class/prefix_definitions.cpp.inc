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
  int begin_height_idx_zm, end_height_idx_zm; // Selected Native momentum-input bounds [-].
  int begin_height_idx_zt, end_height_idx_zt; // Selected Native thermodynamic-input bounds [-].
  int k; // Native one-based level index used in first-column searches [-].

  if (nzmax > NWARNING && f2k::clubb_at_least_debug_level_api(1)) {
    std::cerr << "Warning:  running with vertical grid which is larger than" << NWARNING
              << "levels.\n";
    std::cerr << "This may take a lot of CPU time and memory.\n";
  }

  gr.nzm = nzmax;
  gr.nzt = nzmax - 1;
  begin_height_idx_zm = 1;
  end_height_idx_zm = gr.nzm;
  begin_height_idx_zt = 1;
  end_height_idx_zt = gr.nzt;

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
        k = gr.nzt;
        while (thermodynamic_heights_host(0, k - 1) >= zm_init_host(0) && k > 1)
          k = k - 1;
      } else {
        k = 1;
        while (thermodynamic_heights_host(0, k - 1) >= zm_init_host(0) && k < gr.nzt)
          k = k + 1;
      }

      if (thermodynamic_heights_host(0, k - 1) <= zm_init_host(0)) {
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
        k = gr.nzt;
        while (thermodynamic_heights_host(0, k - 1) > zm_top_host(0) && k > 1)
          k = k - 1;
      } else {
        k = 1;
        while (thermodynamic_heights_host(0, k - 1) > zm_top_host(0) && k < gr.nzt)
          k = k + 1;
      }

      if (zm_top_host(0) < thermodynamic_heights_host(0, k - 1)) {
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
        k = 1;
        while (momentum_heights_host(0, k - 1) < zm_init_host(0) && k < gr.nzm)
          k = k + 1;
      } else {
        k = gr.nzm;
        while (momentum_heights_host(0, k - 1) < zm_init_host(0) && k > 1)
          k = k - 1;
      }

      if (momentum_heights_host(0, k - 1) < zm_init_host(0)) {
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
        k = gr.nzm;
        while (momentum_heights_host(0, k - 1) > zm_top_host(0) && k > 1)
          k = k - 1;
      } else {
        k = 1;
        while (momentum_heights_host(0, k - 1) > zm_top_host(0) && k < gr.nzm)
          k = k + 1;
      }

      if (momentum_heights_host(0, k - 1) > zm_top_host(0)) {
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
                         std::make_pair(begin_height_idx_zm - 1, end_height_idx_zm))),
                     to_execution_space(Kokkos::subview(
                         thermodynamic_heights, Kokkos::ALL,
                         std::make_pair(begin_height_idx_zt - 1, end_height_idx_zt))),
                     gr, err_info // InOut
  );

  if (f2k::clubb_at_least_debug_level_api(0)) {
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
              for (int k = 1; k <= nzm; ++k)
                zm(i, k - 1) = zm_init(i) + real_t(k - 1) * deltaz(i);
            } else {
              for (int k = 1; k <= nzm; ++k)
                zm(i, k - 1) = zm_init(i) + real_t(nzm - k) * deltaz(i);
            }
            for (int k = 1; k <= nzt; ++k)
              zt(i, k - 1) = real_t(0.5) * (zm(i, k) + zm(i, k - 1));
          } else if (grid_type == 2) {
            // Stretched thermodynamic input; construct internal momentum midpoints.
            for (int k = 1; k <= nzt; ++k)
              zt(i, k - 1) = thermodynamic_heights(i, k - 1);
            for (int k = 2; k <= nzm - 1; ++k)
              zm(i, k - 1) = real_t(0.5) * (zt(i, k - 1) + zt(i, k - 2));
            if (l_ascending_grid) {
              zm(i, 0) = zm_init(i);
              zm(i, nzm - 1) = zt(i, nzt - 1) + real_t(0.5) * (zt(i, nzt - 1) - zt(i, nzt - 2));
            } else {
              zm(i, nzm - 1) = zm_init(i);
              zm(i, 0) = zt(i, 0) + real_t(0.5) * (zt(i, 0) - zt(i, 1));
            }
          } else if (grid_type == 3) {
            // Stretched momentum input, followed by thermodynamic midpoints.
            for (int k = 1; k <= nzm; ++k)
              zm(i, k - 1) = momentum_heights(i, k - 1);
            for (int k = 1; k <= nzt; ++k)
              zt(i, k - 1) = real_t(0.5) * (zm(i, k) + zm(i, k - 1));
          }
        } else {
          // Implemented host supplies both staggered sets of heights.
          for (int k = 1; k <= nzm; ++k)
            zm(i, k - 1) = momentum_heights(i, k - 1);
          for (int k = 1; k <= nzt; ++k)
            zt(i, k - 1) = thermodynamic_heights(i, k - 1);
        }

        // Thermodynamic spacing centered on momentum levels; descending values
        // remain negative as in Native.
        for (int k = 2; k <= nzm - 1; ++k)
          dzm(i, k - 1) = zt(i, k - 1) - zt(i, k - 2);
        if (l_ascending_grid) {
          dzm(i, 0) = real_t(2) * (zt(i, 0) - zm(i, 0));
          dzm(i, nzm - 1) = dzm(i, nzm - 2);
        } else {
          dzm(i, nzm - 1) = real_t(2) * (zm(i, nzm - 1) - zt(i, nzt - 1));
          dzm(i, 0) = dzm(i, 1);
        }

        // Momentum spacing centered on thermodynamic levels.
        for (int k = 1; k <= nzt; ++k)
          dzt(i, k - 1) = zm(i, k) - zm(i, k - 1);
        for (int k = 1; k <= nzm; ++k)
          invrs_dzm(i, k - 1) = one / dzm(i, k - 1);
        for (int k = 1; k <= nzt; ++k)
          invrs_dzt(i, k - 1) = one / dzt(i, k - 1);
      });

  // Preserve the Native call order before assigning direction/boundary metadata.
  calc_zm2zt_weights(nzt, ngrdcol,      // In
                     gr);               // InOut
  calc_zt2zm_weights(nzm, nzt, ngrdcol, // In
                     gr);               // InOut

  if (l_ascending_grid) {
    gr.k_lb_zm = 1;
    gr.k_lb_zt = 1;
    gr.k_ub_zm = nzm;
    gr.k_ub_zt = nzt;
    gr.grid_dir_indx = 1;
    gr.grid_dir = one;
  } else {
    gr.k_lb_zm = nzm;
    gr.k_lb_zt = nzt;
    gr.k_ub_zm = 1;
    gr.k_ub_zt = 1;
    gr.grid_dir_indx = -1;
    gr.grid_dir = -one;
  }
}
