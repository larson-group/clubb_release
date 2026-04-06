! advance_helper_module_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module advance_helper_module

subroutine f2py_calculate_thlp2_rad(ngrdcol, nzm, nzt, &
    rcm, thlprcp, radht, clubb_params, thlp2_forcing)

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use derived_type_storage, only: stored_grid
  use advance_helper_module, only: calculate_thlp2_rad

  implicit none

  integer, intent(in) :: ngrdcol, nzm, nzt
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: rcm
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: thlprcp
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: radht
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params
  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: thlp2_forcing

  call calculate_thlp2_rad(ngrdcol, nzm, nzt, stored_grid, &
    rcm, thlprcp, radht, clubb_params, thlp2_forcing)

end subroutine f2py_calculate_thlp2_rad

subroutine f2py_calc_ri_zm(nzm, ngrdcol, bv_freq_sqd, shear, lim_bv, lim_shear, &
    ri_zm)

  use clubb_precision, only: core_rknd
  use advance_helper_module, only: calc_Ri_zm

  implicit none

  integer, intent(in) :: nzm, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: bv_freq_sqd, shear
  real(core_rknd), intent(in) :: lim_bv, lim_shear
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: ri_zm

  call calc_Ri_zm(nzm, ngrdcol, bv_freq_sqd, shear, lim_bv, lim_shear, ri_zm)

end subroutine f2py_calc_ri_zm

subroutine f2py_vertical_avg(total_idx, rho_ds, field, dz, avg_out)

  use clubb_precision, only: core_rknd
  use advance_helper_module, only: vertical_avg

  implicit none

  integer, intent(in) :: total_idx
  real(core_rknd), dimension(total_idx), intent(in) :: rho_ds, field, dz
  real(core_rknd), intent(out) :: avg_out

  avg_out = vertical_avg(total_idx, rho_ds, field, dz)

end subroutine f2py_vertical_avg

subroutine f2py_vertical_integral(total_idx, rho_ds, field, dz, integral_out)

  use clubb_precision, only: core_rknd
  use advance_helper_module, only: vertical_integral

  implicit none

  integer, intent(in) :: total_idx
  real(core_rknd), dimension(total_idx), intent(in) :: rho_ds, field, dz
  real(core_rknd), intent(out) :: integral_out

  integral_out = vertical_integral(total_idx, rho_ds, field, dz)

end subroutine f2py_vertical_integral

subroutine f2py_pvertinterp(nzt, ngrdcol, p_mid, p_out, input_var, interp_var)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use advance_helper_module, only: pvertinterp

  implicit none

  integer, intent(in) :: nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: p_mid, input_var
  real(core_rknd), intent(in) :: p_out
  real(core_rknd), dimension(ngrdcol), intent(out) :: interp_var

  call pvertinterp(nzt, ngrdcol, stored_grid, p_mid, p_out, input_var, interp_var)

end subroutine f2py_pvertinterp

subroutine f2py_smooth_heaviside_peskin(nz, ngrdcol, input, smth_range, &
    smth_output)

  use clubb_precision, only: core_rknd
  use advance_helper_module, only: smooth_heaviside_peskin

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: input
  real(core_rknd), intent(in) :: smth_range
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: smth_output

  smth_output = smooth_heaviside_peskin(nz, ngrdcol, input, smth_range)

end subroutine f2py_smooth_heaviside_peskin

subroutine f2py_smooth_max_scalar_array(nz, ngrdcol, input_var1, input_var2, &
    smth_coef, output_var)

  use clubb_precision, only: core_rknd
  use advance_helper_module, only: smooth_max

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), intent(in) :: input_var1, smth_coef
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: input_var2
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: output_var

  output_var = smooth_max(nz, ngrdcol, input_var1, input_var2, smth_coef)

end subroutine f2py_smooth_max_scalar_array

subroutine f2py_smooth_max_array_scalar(nz, ngrdcol, input_var1, input_var2, &
    smth_coef, output_var)

  use clubb_precision, only: core_rknd
  use advance_helper_module, only: smooth_max

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: input_var1
  real(core_rknd), intent(in) :: input_var2, smth_coef
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: output_var

  output_var = smooth_max(nz, ngrdcol, input_var1, input_var2, smth_coef)

end subroutine f2py_smooth_max_array_scalar

subroutine f2py_smooth_min_scalar_array(nz, ngrdcol, input_var1, input_var2, &
    smth_coef, output_var)

  use clubb_precision, only: core_rknd
  use advance_helper_module, only: smooth_min

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), intent(in) :: input_var1, smth_coef
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: input_var2
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: output_var

  output_var = smooth_min(nz, ngrdcol, input_var1, input_var2, smth_coef)

end subroutine f2py_smooth_min_scalar_array

subroutine f2py_smooth_min_array_scalar(nz, ngrdcol, input_var1, input_var2, &
    smth_coef, output_var)

  use clubb_precision, only: core_rknd
  use advance_helper_module, only: smooth_min

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: input_var1
  real(core_rknd), intent(in) :: input_var2, smth_coef
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: output_var

  output_var = smooth_min(nz, ngrdcol, input_var1, input_var2, smth_coef)

end subroutine f2py_smooth_min_array_scalar

subroutine f2py_calc_xpwp_2d(nzm, nzt, ngrdcol, km_zm, xm, xpwp)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use advance_helper_module, only: calc_xpwp

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: km_zm
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: xm
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: xpwp

  call calc_xpwp(nzm, nzt, ngrdcol, stored_grid, km_zm, xm, xpwp)

end subroutine f2py_calc_xpwp_2d

subroutine f2py_lscale_width_vert_avg(nzm, ngrdcol, smth_type, var_profile, &
    lscale_zm, rho_ds_zm, var_below_ground_value, avg_out)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use advance_helper_module, only: Lscale_width_vert_avg

  implicit none

  integer, intent(in) :: nzm, ngrdcol, smth_type
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: var_profile
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: lscale_zm, rho_ds_zm
  real(core_rknd), intent(in) :: var_below_ground_value
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: avg_out

  avg_out = Lscale_width_vert_avg(nzm, ngrdcol, stored_grid, smth_type, &
    var_profile, lscale_zm, rho_ds_zm, var_below_ground_value)

end subroutine f2py_lscale_width_vert_avg

subroutine f2py_wp2_term_splat_lhs(nzm, nzt, ngrdcol, c_wp2_splat, &
    brunt_vaisala_freq_sqd_splat, lhs_splat_wp2)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use advance_helper_module, only: wp2_term_splat_lhs

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol), intent(in) :: c_wp2_splat
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: brunt_vaisala_freq_sqd_splat
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: lhs_splat_wp2

  call wp2_term_splat_lhs(nzm, nzt, ngrdcol, stored_grid, c_wp2_splat, &
    brunt_vaisala_freq_sqd_splat, lhs_splat_wp2)

end subroutine f2py_wp2_term_splat_lhs

subroutine f2py_wp3_term_splat_lhs(nzm, nzt, ngrdcol, c_wp2_splat, &
    brunt_vaisala_freq_sqd_splat, lhs_splat_wp3)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use advance_helper_module, only: wp3_term_splat_lhs

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol), intent(in) :: c_wp2_splat
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: brunt_vaisala_freq_sqd_splat
  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: lhs_splat_wp3

  call wp3_term_splat_lhs(nzm, nzt, ngrdcol, stored_grid, c_wp2_splat, &
    brunt_vaisala_freq_sqd_splat, lhs_splat_wp3)

end subroutine f2py_wp3_term_splat_lhs

subroutine f2py_calc_brunt_vaisala_freq_sqd(nzm, nzt, ngrdcol, &
    thlm, exner, rtm, rcm, p_in_Pa, thvm, ice_supersat_frac, &
    saturation_formula, &
    l_brunt_vaisala_freq_moist, l_use_thvm_in_bv_freq, &
    l_modify_limiters_for_cnvg_test, &
    bv_efold, T0, &
    brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed, &
    brunt_vaisala_freq_sqd_dry, brunt_vaisala_freq_sqd_moist, &
    brunt_vaisala_freq_sqd_smth)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use advance_helper_module, only: calc_brunt_vaisala_freq_sqd

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    thlm, exner, rtm, rcm, p_in_Pa, thvm, ice_supersat_frac
  integer, intent(in) :: saturation_formula
  logical, intent(in) :: l_brunt_vaisala_freq_moist
  logical, intent(in) :: l_use_thvm_in_bv_freq
  logical, intent(in) :: l_modify_limiters_for_cnvg_test
  real(core_rknd), dimension(ngrdcol), intent(in) :: bv_efold
  real(core_rknd), intent(in) :: T0
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: &
    brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed, &
    brunt_vaisala_freq_sqd_dry, brunt_vaisala_freq_sqd_moist, &
    brunt_vaisala_freq_sqd_smth

  call calc_brunt_vaisala_freq_sqd(nzm, nzt, ngrdcol, stored_grid, &
    thlm, exner, rtm, rcm, p_in_Pa, thvm, ice_supersat_frac, &
    saturation_formula, &
    l_brunt_vaisala_freq_moist, &
    l_use_thvm_in_bv_freq, &
    l_modify_limiters_for_cnvg_test, &
    bv_efold, T0, &
    brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed, &
    brunt_vaisala_freq_sqd_dry, brunt_vaisala_freq_sqd_moist, &
    brunt_vaisala_freq_sqd_smth)

end subroutine f2py_calc_brunt_vaisala_freq_sqd

subroutine f2py_compute_cx_fnc_richardson(nzm, nzt, ngrdcol, &
    Lscale_zm, ddzt_umvm_sqd, rho_ds_zm, &
    brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed, &
    clubb_params, &
    l_use_shear_Richardson, l_modify_limiters_for_cnvg_test, &
    Cx_fnc_Richardson)

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use derived_type_storage, only: stored_grid
  use advance_helper_module, only: compute_Cx_fnc_Richardson

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: &
    Lscale_zm, ddzt_umvm_sqd, rho_ds_zm, &
    brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params
  logical, intent(in) :: l_use_shear_Richardson
  logical, intent(in) :: l_modify_limiters_for_cnvg_test
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: Cx_fnc_Richardson

  call compute_Cx_fnc_Richardson(nzm, nzt, ngrdcol, stored_grid, &
    Lscale_zm, ddzt_umvm_sqd, rho_ds_zm, &
    brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed, &
    clubb_params, &
    l_use_shear_Richardson, &
    l_modify_limiters_for_cnvg_test, &
    Cx_fnc_Richardson)

end subroutine f2py_compute_cx_fnc_richardson

subroutine f2py_calc_stability_correction(nzm, nzt, ngrdcol, &
    thlm, Lscale_zm, em, exner, rtm, rcm, &
    p_in_Pa, thvm, ice_supersat_frac, &
    lambda0_stability_coef, bv_efold, T0, &
    saturation_formula, &
    l_brunt_vaisala_freq_moist, l_use_thvm_in_bv_freq, &
    l_modify_limiters_for_cnvg_test, &
    stability_correction)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use advance_helper_module, only: calc_stability_correction

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    thlm, exner, rtm, rcm, p_in_Pa, thvm, ice_supersat_frac
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: Lscale_zm, em
  real(core_rknd), dimension(ngrdcol), intent(in) :: lambda0_stability_coef, bv_efold
  real(core_rknd), intent(in) :: T0
  integer, intent(in) :: saturation_formula
  logical, intent(in) :: l_brunt_vaisala_freq_moist
  logical, intent(in) :: l_use_thvm_in_bv_freq
  logical, intent(in) :: l_modify_limiters_for_cnvg_test
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: stability_correction

  call calc_stability_correction(nzm, nzt, ngrdcol, stored_grid, &
    thlm, Lscale_zm, em, exner, rtm, rcm, &
    p_in_Pa, thvm, ice_supersat_frac, &
    lambda0_stability_coef, bv_efold, T0, &
    saturation_formula, &
    l_brunt_vaisala_freq_moist, &
    l_use_thvm_in_bv_freq, &
    l_modify_limiters_for_cnvg_test, &
    stability_correction)

end subroutine f2py_calc_stability_correction
