// Merge in clubb_api_module.hpp at the original Native wrapper positions.
// Both APIs keep Native error checks after the core calls.
template <class Momentum, class Thermodynamic>
void setup_grid_heights_api(bool l_implemented, bool l_ascending_grid,    // In
                            int grid_type, real_t deltaz, real_t zm_init, // In
                            Momentum momentum_heights,                    // In
                            grid& gr,                                     // InOut
                            Thermodynamic thermodynamic_heights,          // In
                            err_info_type& err_info_api) {                // InOut
  // Description:
  //   Reset an allocated single-column grid's heights and interpolation weights.
  using Host1D = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  Host1D deltaz_col("setup_grid_heights_api::deltaz_col", 1);
  Host1D zm_init_col("setup_grid_heights_api::zm_init_col", 1);
  Real2D momentum_heights_col("setup_grid_heights_api::momentum_heights_col", 1, gr.nzm);
  Real2D thermodynamic_heights_col("setup_grid_heights_api::thermodynamic_heights_col", 1, gr.nzt);
  deltaz_col(0) = deltaz;
  zm_init_col(0) = zm_init;
  Kokkos::deep_copy(Kokkos::subview(momentum_heights_col, 0, Kokkos::ALL),
                    to_execution_space(momentum_heights));
  Kokkos::deep_copy(Kokkos::subview(thermodynamic_heights_col, 0, Kokkos::ALL),
                    to_execution_space(thermodynamic_heights));
  setup_grid_heights(gr.nzm, gr.nzt, 1,                               // In
                     l_implemented, l_ascending_grid,                 // In
                     grid_type, deltaz_col, zm_init_col,              // In
                     momentum_heights_col, thermodynamic_heights_col, // In
                     gr, err_info_api);                               // InOut
  const auto err_code = to_host_space(err_info_api.err_code);
  for (std::size_t i = 0; i < err_code.extent(0); ++i) {
    if (err_code(i) == clubb_fatal_error) {
      std::cerr << err_info_api.err_header_global << '\n';
      std::cerr << "Error in CLUBB calling setup_grid_heights\n";
      break;
    }
  }
}

template <class Spacing, class Start, class Momentum, class Thermodynamic>
void setup_grid_heights_api(int nzm, int nzt, int ngrdcol,                // In
                            bool l_implemented, bool l_ascending_grid,    // In
                            int grid_type, Spacing deltaz, Start zm_init, // In
                            Momentum momentum_heights,                    // In
                            grid& gr,                                     // InOut
                            Thermodynamic thermodynamic_heights,          // In
                            err_info_type& err_info_api) {                // InOut
  // Description:
  //   Reset allocated column grids' heights and interpolation weights.
  setup_grid_heights(nzm, nzt, ngrdcol,                       // In
                     l_implemented, l_ascending_grid,         // In
                     grid_type, deltaz, zm_init,              // In
                     momentum_heights, thermodynamic_heights, // In
                     gr, err_info_api);                       // InOut
  const auto err_code = to_host_space(err_info_api.err_code);
  for (std::size_t i = 0; i < err_code.extent(0); ++i) {
    if (err_code(i) == clubb_fatal_error) {
      std::cerr << err_info_api.err_header_global << '\n';
      std::cerr << "Error in CLUBB calling setup_grid_heights\n";
      break;
    }
  }
}

template <class Momentum, class Thermodynamic>
void setup_grid_api(int nzmax, real_t sfc_elevation, bool l_implemented,            // In
                    bool l_ascending_grid, int grid_type,                           // In
                    real_t deltaz, real_t zm_init, real_t zm_top,                   // In
                    Momentum momentum_heights, Thermodynamic thermodynamic_heights, // In
                    grid& gr, err_info_type& err_info_api) {                        // InOut
  // Description:
  //   Allocate and initialize the grid for a single-column caller.
  using Host1D = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  Host1D sfc_elevation_col("setup_grid_api::sfc_elevation_col", 1);
  Host1D deltaz_col("setup_grid_api::deltaz_col", 1);
  Host1D zm_init_col("setup_grid_api::zm_init_col", 1);
  Host1D zm_top_col("setup_grid_api::zm_top_col", 1);
  Real2D momentum_heights_col("setup_grid_api::momentum_heights_col", 1, nzmax);
  Real2D thermodynamic_heights_col("setup_grid_api::thermodynamic_heights_col", 1, nzmax - 1);
  sfc_elevation_col(0) = sfc_elevation;
  deltaz_col(0) = deltaz;
  zm_init_col(0) = zm_init;
  zm_top_col(0) = zm_top;
  Kokkos::deep_copy(Kokkos::subview(momentum_heights_col, 0, Kokkos::ALL),
                    to_execution_space(momentum_heights));
  Kokkos::deep_copy(Kokkos::subview(thermodynamic_heights_col, 0, Kokkos::ALL),
                    to_execution_space(thermodynamic_heights));
  setup_grid(nzmax, 1, sfc_elevation_col, l_implemented,                       // In
             l_ascending_grid, grid_type, deltaz_col, zm_init_col, zm_top_col, // In
             momentum_heights_col, thermodynamic_heights_col,                  // In
             gr, err_info_api);                                                // InOut
  const auto err_code = to_host_space(err_info_api.err_code);
  for (std::size_t i = 0; i < err_code.extent(0); ++i) {
    if (err_code(i) == clubb_fatal_error) {
      std::cerr << err_info_api.err_header_global << '\n';
      std::cerr << "Error in CLUBB calling setup_grid\n";
      break;
    }
  }
}

template <class Surface, class Spacing, class Start, class Top, class Momentum, class Thermodynamic>
void setup_grid_api(int nzmax, int ngrdcol, Surface sfc_elevation, bool l_implemented, // In
                    bool l_ascending_grid, int grid_type, Spacing deltaz, Start zm_init,
                    Top zm_top,                                                     // In
                    Momentum momentum_heights, Thermodynamic thermodynamic_heights, // In
                    grid& gr, err_info_type& err_info_api) {                        // InOut
  // Description:
  //   Allocate and initialize the grid for a column chunk.
  setup_grid(nzmax, ngrdcol, sfc_elevation, l_implemented,         // In
             l_ascending_grid, grid_type, deltaz, zm_init, zm_top, // In
             momentum_heights, thermodynamic_heights,              // In
             gr, err_info_api);                                    // InOut
  const auto err_code = to_host_space(err_info_api.err_code);
  for (std::size_t i = 0; i < err_code.extent(0); ++i) {
    if (err_code(i) == clubb_fatal_error) {
      std::cerr << err_info_api.err_header_global << '\n';
      std::cerr << "Error in CLUBB calling setup_grid\n";
      break;
    }
  }
}
