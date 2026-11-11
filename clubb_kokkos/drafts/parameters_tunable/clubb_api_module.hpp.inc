// Merge at the Native calc_derived_params_api single/multi-column wrapper positions.
template <class Vector>
void calc_derived_params_api(const grid& gr, int grid_type, real_t deltaz,       // In
                             Vector clubb_params,                                // In
                             bool l_prescribed_avg_deltaz,                       // In
                             nu_vertical_res_dep& nu_vert_res_dep, real_t& lmin, // InOut
                             real_t& mixt_frac_max_mag                           // InOut
) {
  // Description:
  //   Adjust background eddy diffusivity for vertical spacing in one column.
  // gr: Grid definition [-]. grid_type: Uniform(1), thermodynamic(2), momentum(3) [-].
  // deltaz: Change per height level [m]. clubb_params: Tunable parameter vector [-].
  // l_prescribed_avg_deltaz: Use deltaz as the average spacing when true [-].
  // nu_vert_res_dep: Resolution-dependent diffusivities [m^2/s].
  // lmin: Minimum mixing length [m]. mixt_frac_max_mag: Mixture-fraction limit [-].
  using Host1D = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  Host1D deltaz_col("calc_derived_params_api::deltaz_col", 1);
  Host1D lmin_col(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_derived_params_api::lmin_col"), 1);
  Host1D mixt_frac_max_mag_col(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                                  "calc_derived_params_api::mixt_frac_max_mag_col"),
                               1);
  Real2D clubb_params_col("calc_derived_params_api::clubb_params_col", 1, nparams);
  deltaz_col(0) = deltaz;
  Kokkos::deep_copy(Kokkos::subview(clubb_params_col, 0, Kokkos::ALL),
                    to_execution_space(clubb_params));
  // Native passes its original rank-1 vector by sequence association. C++ uses
  // the corresponding one-column array, with identical values and ordering.
  calc_derived_params(gr, 1, grid_type, deltaz_col, // In
                      clubb_params_col,             // In
                      l_prescribed_avg_deltaz,      // In
                      nu_vert_res_dep, lmin_col,    // InOut
                      mixt_frac_max_mag_col);       // InOut
  lmin = lmin_col(0);
  mixt_frac_max_mag = mixt_frac_max_mag_col(0);
}

template <class Spacing, class Matrix, class MinimumLength, class MixtureLimit>
void calc_derived_params_api(const grid& gr, int ngrdcol, int grid_type, Spacing deltaz, // In
                             Matrix clubb_params,                                        // In
                             bool l_prescribed_avg_deltaz,                               // In
                             nu_vertical_res_dep& nu_vert_res_dep, MinimumLength lmin,   // InOut
                             MixtureLimit mixt_frac_max_mag                              // InOut
) {
  // Description:
  //   Adjust background eddy diffusivity for vertical spacing in a column chunk.
  // gr: Grid definition [-]. ngrdcol: Number of columns [#].
  // grid_type: Uniform(1), thermodynamic(2), momentum(3) input grid [-].
  // deltaz: Change per height level [m]. clubb_params: Tunable parameter matrix [-].
  // l_prescribed_avg_deltaz: Use deltaz as the average spacing when true [-].
  // nu_vert_res_dep: Resolution-dependent diffusivities [m^2/s].
  // lmin: Minimum mixing length [m]. mixt_frac_max_mag: Mixture-fraction limit [-].
  calc_derived_params(gr, ngrdcol, grid_type, deltaz, // In
                      clubb_params,                   // In
                      l_prescribed_avg_deltaz,        // In
                      nu_vert_res_dep, lmin,          // InOut
                      mixt_frac_max_mag);             // InOut
}
