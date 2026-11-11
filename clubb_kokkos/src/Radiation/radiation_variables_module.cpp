// Case-owned radiation fields; see the matching header for lifecycle usage.
#include "radiation_variables_module.hpp"

namespace clubb_kokkos {

void setup_radiation_variables(int ngrdcol, int nzm, int nzt,                 // In
                               radiation_variables_state& radiation_variables // Native module Out
) {
  // Allocate the five net/directional fluxes and four component diagnostics,
  // then initialize both vertical grids to zero in Native loop order.
  radiation_variables.Frad =
      Real2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::Frad"), ngrdcol, nzm);
  radiation_variables.Frad_SW_up = Real2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::Frad_SW_up"), ngrdcol, nzm);
  radiation_variables.Frad_LW_up = Real2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::Frad_LW_up"), ngrdcol, nzm);
  radiation_variables.Frad_SW_down = Real2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::Frad_SW_down"), ngrdcol, nzm);
  radiation_variables.Frad_LW_down = Real2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::Frad_LW_down"), ngrdcol, nzm);
  radiation_variables.radht_LW =
      Real2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::radht_LW"), ngrdcol, nzt);
  radiation_variables.radht_SW =
      Real2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::radht_SW"), ngrdcol, nzt);
  radiation_variables.Frad_SW =
      Real2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::Frad_SW"), ngrdcol, nzm);
  radiation_variables.Frad_LW =
      Real2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::Frad_LW"), ngrdcol, nzm);
  const auto Frad         = radiation_variables.Frad;
  const auto Frad_SW_up   = radiation_variables.Frad_SW_up;
  const auto Frad_LW_up   = radiation_variables.Frad_LW_up;
  const auto Frad_SW_down = radiation_variables.Frad_SW_down;
  const auto Frad_LW_down = radiation_variables.Frad_LW_down;
  const auto radht_LW     = radiation_variables.radht_LW;
  const auto radht_SW     = radiation_variables.radht_SW;
  const auto Frad_SW      = radiation_variables.Frad_SW;
  const auto Frad_LW      = radiation_variables.Frad_LW;
  Kokkos::parallel_for(
      "setup_radiation_variables::momentum", Kokkos::RangePolicy<>(0, ngrdcol * nzm),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        Frad(i, k)         = real_t(0);
        Frad_SW_up(i, k)   = real_t(0);
        Frad_LW_up(i, k)   = real_t(0);
        Frad_SW_down(i, k) = real_t(0);
        Frad_LW_down(i, k) = real_t(0);
        Frad_SW(i, k)      = real_t(0);
        Frad_LW(i, k)      = real_t(0);
      });
  Kokkos::parallel_for(
      "setup_radiation_variables::thermodynamic", Kokkos::RangePolicy<>(0, ngrdcol * nzt),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        radht_LW(i, k) = real_t(0);
        radht_SW(i, k) = real_t(0);
      });
}

void reset_radiation_variables(radiation_variables_state& radiation_variables) {
  // Reset allocated regular radiation fields; do nothing before setup.
  if (!radiation_variables.Frad.is_allocated())
    return;
  const auto Frad         = radiation_variables.Frad;
  const auto Frad_SW_up   = radiation_variables.Frad_SW_up;
  const auto Frad_LW_up   = radiation_variables.Frad_LW_up;
  const auto Frad_SW_down = radiation_variables.Frad_SW_down;
  const auto Frad_LW_down = radiation_variables.Frad_LW_down;
  const auto radht_LW     = radiation_variables.radht_LW;
  const auto radht_SW     = radiation_variables.radht_SW;
  const auto Frad_SW      = radiation_variables.Frad_SW;
  const auto Frad_LW      = radiation_variables.Frad_LW;
  const int ngrdcol = Frad.extent_int(0), nzm = Frad.extent_int(1);
  const int nzt = radht_LW.extent_int(1);
  Kokkos::parallel_for(
      "reset_radiation_variables::momentum", Kokkos::RangePolicy<>(0, ngrdcol * nzm),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        Frad(i, k)         = real_t(0);
        Frad_SW_up(i, k)   = real_t(0);
        Frad_LW_up(i, k)   = real_t(0);
        Frad_SW_down(i, k) = real_t(0);
        Frad_LW_down(i, k) = real_t(0);
        Frad_SW(i, k)      = real_t(0);
        Frad_LW(i, k)      = real_t(0);
      });
  Kokkos::parallel_for(
      "reset_radiation_variables::thermodynamic", Kokkos::RangePolicy<>(0, ngrdcol * nzt),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        radht_LW(i, k) = real_t(0);
        radht_SW(i, k) = real_t(0);
      });
}

void setup_bugsrad_variables(
    ScratchWorkspace& workspace,                            // InOut
    int ngrdcol, int nzmax, Real1D zt_grid, Real1D zm_grid, // In
    Real1D zm_grid_spacing, Real1D p_in_Pa_zm,              // In
    radiation_variables_state& radiation_variables,         // Native module Out
    const parameters_radiation_state& parameters_radiation, // Native module In
    extended_atmosphere_state& extended_atmosphere          // Native module InOut
) {
  // Determine the atmosphere extension, then allocate/zero its radiation fields.
  // This lifecycle contract does not call the excluded BUGSrad physics engine.
  int extended_atmos_bottom_level_local, extended_atmos_top_level_local;
  int extended_atmos_range_size_local, lin_int_buffer_local;
  determine_extended_atmos_bounds(workspace, nzmax, zt_grid, zm_grid, zm_grid_spacing,
                                  p_in_Pa_zm, // In
                                  extended_atmos_bottom_level_local,
                                  extended_atmos_top_level_local,                        // Out
                                  extended_atmos_range_size_local, lin_int_buffer_local, // Out
                                  parameters_radiation, extended_atmosphere);
  radiation_variables.extended_atmos_bottom_level = extended_atmos_bottom_level_local;
  radiation_variables.extended_atmos_top_level    = extended_atmos_top_level_local;
  radiation_variables.extended_atmos_range_size   = extended_atmos_range_size_local;
  radiation_variables.lin_int_buffer              = lin_int_buffer_local;
  const int rad_zt_dim = (nzmax - 1) + lin_int_buffer_local + extended_atmos_range_size_local;
  const int rad_zm_dim = (nzmax - 1) + lin_int_buffer_local + extended_atmos_range_size_local + 1;
  radiation_variables.T_in_K = RadiationDouble2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::T_in_K"), ngrdcol, rad_zt_dim);
  radiation_variables.rcil = RadiationDouble2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::rcil"), ngrdcol, rad_zt_dim);
  radiation_variables.o3l = RadiationDouble2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::o3l"), ngrdcol, rad_zt_dim);
  radiation_variables.rsm_rad = RadiationDouble2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::rsm_rad"), ngrdcol, rad_zt_dim);
  radiation_variables.rcm_in_cloud_rad = RadiationDouble2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::rcm_in_cloud_rad"), ngrdcol,
      rad_zt_dim);
  radiation_variables.cloud_frac_rad = RadiationDouble2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::cloud_frac_rad"), ngrdcol,
      rad_zt_dim);
  radiation_variables.ice_supersat_frac_rad = RadiationDouble2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::ice_supersat_frac_rad"), ngrdcol,
      rad_zt_dim);
  radiation_variables.p_in_mb = RadiationDouble2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::p_in_mb"), ngrdcol, rad_zt_dim);
  radiation_variables.sp_humidity =
      RadiationDouble2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::sp_humidity"),
                        ngrdcol, rad_zt_dim);
  radiation_variables.radht_SW_rad =
      RadiationDouble2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::radht_SW_rad"),
                        ngrdcol, rad_zt_dim);
  radiation_variables.radht_LW_rad =
      RadiationDouble2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::radht_LW_rad"),
                        ngrdcol, rad_zt_dim);
  radiation_variables.Frad_uLW = RadiationDouble2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::Frad_uLW"), ngrdcol, rad_zm_dim);
  radiation_variables.Frad_dLW = RadiationDouble2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::Frad_dLW"), ngrdcol, rad_zm_dim);
  radiation_variables.Frad_uSW = RadiationDouble2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::Frad_uSW"), ngrdcol, rad_zm_dim);
  radiation_variables.Frad_dSW = RadiationDouble2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::Frad_dSW"), ngrdcol, rad_zm_dim);
  radiation_variables.fdswcl = RadiationDouble2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::fdswcl"), ngrdcol, rad_zm_dim);
  radiation_variables.fuswcl = RadiationDouble2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::fuswcl"), ngrdcol, rad_zm_dim);
  radiation_variables.fdlwcl = RadiationDouble2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::fdlwcl"), ngrdcol, rad_zm_dim);
  radiation_variables.fulwcl = RadiationDouble2D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "radiation::fulwcl"), ngrdcol, rad_zm_dim);
  const auto T_in_K                = radiation_variables.T_in_K;
  const auto rcil                  = radiation_variables.rcil;
  const auto o3l                   = radiation_variables.o3l;
  const auto rsm_rad               = radiation_variables.rsm_rad;
  const auto rcm_in_cloud_rad      = radiation_variables.rcm_in_cloud_rad;
  const auto cloud_frac_rad        = radiation_variables.cloud_frac_rad;
  const auto ice_supersat_frac_rad = radiation_variables.ice_supersat_frac_rad;
  const auto p_in_mb               = radiation_variables.p_in_mb;
  const auto sp_humidity           = radiation_variables.sp_humidity;
  const auto radht_SW_rad          = radiation_variables.radht_SW_rad;
  const auto radht_LW_rad          = radiation_variables.radht_LW_rad;
  const auto Frad_uLW              = radiation_variables.Frad_uLW;
  const auto Frad_dLW              = radiation_variables.Frad_dLW;
  const auto Frad_uSW              = radiation_variables.Frad_uSW;
  const auto Frad_dSW              = radiation_variables.Frad_dSW;
  const auto fdswcl                = radiation_variables.fdswcl;
  const auto fuswcl                = radiation_variables.fuswcl;
  const auto fdlwcl                = radiation_variables.fdlwcl;
  const auto fulwcl                = radiation_variables.fulwcl;
  Kokkos::parallel_for(
      "setup_bugsrad_variables::thermodynamic", Kokkos::RangePolicy<>(0, ngrdcol * rad_zt_dim),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        T_in_K(i, k)                = 0.0;
        rcil(i, k)                  = 0.0;
        o3l(i, k)                   = 0.0;
        rsm_rad(i, k)               = 0.0;
        rcm_in_cloud_rad(i, k)      = 0.0;
        cloud_frac_rad(i, k)        = 0.0;
        ice_supersat_frac_rad(i, k) = 0.0;
        radht_SW_rad(i, k)          = 0.0;
        radht_LW_rad(i, k)          = 0.0;
        p_in_mb(i, k)               = 0.0;
        sp_humidity(i, k)           = 0.0;
      });
  Kokkos::parallel_for(
      "setup_bugsrad_variables::momentum", Kokkos::RangePolicy<>(0, ngrdcol * rad_zm_dim),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        Frad_uLW(i, k) = 0.0;
        Frad_dLW(i, k) = 0.0;
        Frad_uSW(i, k) = 0.0;
        Frad_dSW(i, k) = 0.0;
        fdswcl(i, k)   = 0.0;
        fuswcl(i, k)   = 0.0;
        fdlwcl(i, k)   = 0.0;
        fulwcl(i, k)   = 0.0;
      });
}

void reset_bugsrad_variables(radiation_variables_state& radiation_variables) {
  // Reset allocated extension diagnostics; retain the source's early return.
  if (!radiation_variables.T_in_K.is_allocated())
    return;
  const auto T_in_K                = radiation_variables.T_in_K;
  const auto rcil                  = radiation_variables.rcil;
  const auto o3l                   = radiation_variables.o3l;
  const auto rsm_rad               = radiation_variables.rsm_rad;
  const auto rcm_in_cloud_rad      = radiation_variables.rcm_in_cloud_rad;
  const auto cloud_frac_rad        = radiation_variables.cloud_frac_rad;
  const auto ice_supersat_frac_rad = radiation_variables.ice_supersat_frac_rad;
  const auto p_in_mb               = radiation_variables.p_in_mb;
  const auto sp_humidity           = radiation_variables.sp_humidity;
  const auto radht_SW_rad          = radiation_variables.radht_SW_rad;
  const auto radht_LW_rad          = radiation_variables.radht_LW_rad;
  const auto Frad_uLW              = radiation_variables.Frad_uLW;
  const auto Frad_dLW              = radiation_variables.Frad_dLW;
  const auto Frad_uSW              = radiation_variables.Frad_uSW;
  const auto Frad_dSW              = radiation_variables.Frad_dSW;
  const auto fdswcl                = radiation_variables.fdswcl;
  const auto fuswcl                = radiation_variables.fuswcl;
  const auto fdlwcl                = radiation_variables.fdlwcl;
  const auto fulwcl                = radiation_variables.fulwcl;
  const int ngrdcol = T_in_K.extent_int(0), rad_zt_dim = T_in_K.extent_int(1);
  const int rad_zm_dim = Frad_uLW.extent_int(1);
  Kokkos::parallel_for(
      "reset_bugsrad_variables::thermodynamic", Kokkos::RangePolicy<>(0, ngrdcol * rad_zt_dim),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        T_in_K(i, k)                = 0.0;
        rcil(i, k)                  = 0.0;
        o3l(i, k)                   = 0.0;
        rsm_rad(i, k)               = 0.0;
        rcm_in_cloud_rad(i, k)      = 0.0;
        cloud_frac_rad(i, k)        = 0.0;
        ice_supersat_frac_rad(i, k) = 0.0;
        p_in_mb(i, k)               = 0.0;
        sp_humidity(i, k)           = 0.0;
        radht_SW_rad(i, k)          = 0.0;
        radht_LW_rad(i, k)          = 0.0;
      });
  Kokkos::parallel_for(
      "reset_bugsrad_variables::momentum", Kokkos::RangePolicy<>(0, ngrdcol * rad_zm_dim),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        Frad_uLW(i, k) = 0.0;
        Frad_dLW(i, k) = 0.0;
        Frad_uSW(i, k) = 0.0;
        Frad_dSW(i, k) = 0.0;
        fdswcl(i, k)   = 0.0;
        fuswcl(i, k)   = 0.0;
        fdlwcl(i, k)   = 0.0;
        fulwcl(i, k)   = 0.0;
      });
}

void cleanup_radiation_variables(radiation_variables_state& radiation_variables) {
  // Release regular radiation fields only after the source allocation guard.
  if (radiation_variables.Frad.is_allocated()) {
    radiation_variables.Frad         = Real2D{};
    radiation_variables.Frad_SW_up   = Real2D{};
    radiation_variables.Frad_LW_up   = Real2D{};
    radiation_variables.Frad_SW_down = Real2D{};
    radiation_variables.Frad_LW_down = Real2D{};
    radiation_variables.radht_LW     = Real2D{};
    radiation_variables.radht_SW     = Real2D{};
    radiation_variables.Frad_SW      = Real2D{};
    radiation_variables.Frad_LW      = Real2D{};
  }
}

void cleanup_bugsrad_variables(radiation_variables_state& radiation_variables) {
  // Reset the extension indices, then release allocated BUGS diagnostic arrays.
  radiation_variables.extended_atmos_bottom_level = -1;
  radiation_variables.extended_atmos_top_level    = 0;
  radiation_variables.extended_atmos_range_size   = 0;
  radiation_variables.lin_int_buffer              = 0;
  if (radiation_variables.T_in_K.is_allocated()) {
    radiation_variables.T_in_K                = RadiationDouble2D{};
    radiation_variables.rcil                  = RadiationDouble2D{};
    radiation_variables.o3l                   = RadiationDouble2D{};
    radiation_variables.rsm_rad               = RadiationDouble2D{};
    radiation_variables.rcm_in_cloud_rad      = RadiationDouble2D{};
    radiation_variables.cloud_frac_rad        = RadiationDouble2D{};
    radiation_variables.ice_supersat_frac_rad = RadiationDouble2D{};
    radiation_variables.p_in_mb               = RadiationDouble2D{};
    radiation_variables.sp_humidity           = RadiationDouble2D{};
    radiation_variables.radht_SW_rad          = RadiationDouble2D{};
    radiation_variables.radht_LW_rad          = RadiationDouble2D{};
    radiation_variables.Frad_uLW              = RadiationDouble2D{};
    radiation_variables.Frad_dLW              = RadiationDouble2D{};
    radiation_variables.Frad_uSW              = RadiationDouble2D{};
    radiation_variables.Frad_dSW              = RadiationDouble2D{};
    radiation_variables.fdswcl                = RadiationDouble2D{};
    radiation_variables.fuswcl                = RadiationDouble2D{};
    radiation_variables.fdlwcl                = RadiationDouble2D{};
    radiation_variables.fulwcl                = RadiationDouble2D{};
  }
}

} // namespace clubb_kokkos
