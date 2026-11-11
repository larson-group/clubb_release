// Merge after the authoritative generated grid type, before existing operators.
void setup_grid(int nzmax, int ngrdcol, Real1D sfc_elevation, bool l_implemented, // In
                bool l_ascending_grid, int grid_type, Real1D deltaz, Real1D zm_init,
                Real1D zm_top,                                                 // In
                Real2D momentum_heights, Real2D thermodynamic_heights,         // In
                grid& gr, err_info_type& err_info);                            // InOut
void cleanup_grid_api(grid& gr);                                               // InOut
void setup_grid_heights(int nzm, int nzt, int ngrdcol,                         // In
                        bool l_implemented, bool l_ascending_grid,             // In
                        int grid_type, Real1D deltaz, Real1D zm_init,          // In
                        Real2D momentum_heights, Real2D thermodynamic_heights, // In
                        grid& gr, err_info_type& err_info);                    // InOut

// Shape/memory-space adapters retain the caller-owned grid and error state.
// Compatible execution Views alias; host and strided inputs stage once.
template <class Surface, class Spacing, class Start, class Top, class Momentum, class Thermodynamic>
void setup_grid(int nzmax, int ngrdcol, Surface sfc_elevation, bool l_implemented, // In
                bool l_ascending_grid, int grid_type, Spacing deltaz, Start zm_init,
                Top zm_top,                                                     // In
                Momentum momentum_heights, Thermodynamic thermodynamic_heights, // In
                grid& gr, err_info_type& err_info) {                            // InOut
  setup_grid(nzmax, ngrdcol, to_execution_space(sfc_elevation), l_implemented,  // In
             l_ascending_grid, grid_type, to_execution_space(deltaz),           // In
             to_execution_space(zm_init), to_execution_space(zm_top),           // In
             to_execution_space(momentum_heights), to_execution_space(thermodynamic_heights), // In
             gr, err_info); // InOut
}
template <class Spacing, class Start, class Momentum, class Thermodynamic>
void setup_grid_heights(int nzm, int nzt, int ngrdcol,                                   // In
                        bool l_implemented, bool l_ascending_grid,                       // In
                        int grid_type, Spacing deltaz, Start zm_init,                    // In
                        Momentum momentum_heights, Thermodynamic thermodynamic_heights,  // In
                        grid& gr, err_info_type& err_info) {                             // InOut
  setup_grid_heights(nzm, nzt, ngrdcol,                                                  // In
                     l_implemented, l_ascending_grid,                                    // In
                     grid_type, to_execution_space(deltaz), to_execution_space(zm_init), // In
                     to_execution_space(momentum_heights),                               // In
                     to_execution_space(thermodynamic_heights),                          // In
                     gr, err_info);                                                      // InOut
}
