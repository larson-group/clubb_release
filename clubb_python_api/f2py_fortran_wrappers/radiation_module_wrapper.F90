! radiation_module_wrapper.F90 — minimal wrappers for simplified radiation support

subroutine f2py_cos_solar_zen(day, month, year, current_time, lat_vals, lon_vals, cos_zen)

  use clubb_precision, only: core_rknd, time_precision
  use cos_solar_zen_module, only: cos_solar_zen

  implicit none

  integer, intent(in) :: day, month, year
  real(time_precision), intent(in) :: current_time
  real(core_rknd), intent(in) :: lat_vals, lon_vals
  real(core_rknd), intent(out) :: cos_zen

  cos_zen = real(cos_solar_zen(day, month, year, current_time, lat_vals, lon_vals), kind=core_rknd)

end subroutine f2py_cos_solar_zen

subroutine f2py_set_simplified_radiation_params( &
    f0_in, f1_in, kappa_in, eff_drop_radius_in, alvdr_in, gc_in, omega_in, &
    l_rad_above_cloud_in, l_sw_radiation_in, l_fix_cos_solar_zen_in, &
    nparam_in, fs_values_in, cos_solar_zen_values_in, cos_solar_zen_times_in)

  use clubb_precision, only: core_rknd, dp
  use parameters_radiation, only: &
    f0, f1, kappa, eff_drop_radius, alvdr, gc, omega, &
    l_rad_above_cloud, l_sw_radiation, l_fix_cos_solar_zen, nparam, &
    fs_values, cos_solar_zen_values, cos_solar_zen_times

  implicit none

  real(core_rknd), intent(in) :: f0_in, f1_in, kappa_in, eff_drop_radius_in, alvdr_in, gc_in, omega_in
  logical, intent(in) :: l_rad_above_cloud_in, l_sw_radiation_in, l_fix_cos_solar_zen_in
  integer, intent(in) :: nparam_in
  real(core_rknd), dimension(20), intent(in) :: fs_values_in, cos_solar_zen_values_in, cos_solar_zen_times_in
  integer :: ncopy

  f0 = f0_in
  f1 = f1_in
  kappa = kappa_in
  eff_drop_radius = eff_drop_radius_in
  alvdr = real(alvdr_in, kind=dp)
  gc = gc_in
  omega = omega_in

  l_rad_above_cloud = l_rad_above_cloud_in
  l_sw_radiation = l_sw_radiation_in
  l_fix_cos_solar_zen = l_fix_cos_solar_zen_in

  ncopy = max(1, min(20, nparam_in))
  nparam = ncopy

  fs_values = 0.0_core_rknd
  cos_solar_zen_values = -999.0_core_rknd
  cos_solar_zen_times = -999.0_core_rknd

  fs_values(1:ncopy) = fs_values_in(1:ncopy)
  cos_solar_zen_values(1:ncopy) = cos_solar_zen_values_in(1:ncopy)
  cos_solar_zen_times(1:ncopy) = cos_solar_zen_times_in(1:ncopy)

end subroutine f2py_set_simplified_radiation_params


subroutine f2py_sunray_sw(ngrdcol, nzt, fs0, amu0, rho, rcm, frad_sw)

  use clubb_precision, only: core_rknd
  use rad_lwsw_module, only: sunray_sw
  use parameters_radiation, only: eff_drop_radius, omega, alvdr, gc
  use derived_type_storage, only: stored_grid

  implicit none

  integer, intent(in) :: ngrdcol, nzt
  real(core_rknd), intent(in) :: fs0, amu0
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: rho, rcm
  real(core_rknd), dimension(ngrdcol, nzt+1), intent(out) :: frad_sw

  logical, parameter :: l_center = .true.

  call sunray_sw(ngrdcol, nzt, rcm, rho, amu0, &
                 stored_grid%dzt, stored_grid%zm, stored_grid%zt, &
                 eff_drop_radius, real(alvdr, kind=core_rknd), &
                 gc, fs0, omega, l_center, &
                 frad_sw)

end subroutine f2py_sunray_sw
