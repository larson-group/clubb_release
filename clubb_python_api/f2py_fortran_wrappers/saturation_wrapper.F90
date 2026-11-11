! saturation_wrapper.F90 — wrappers organized by source module

subroutine f2py_sat_mixrat_liq_2d(nz, ngrdcol, p_in_Pa, T_in_K, &
    saturation_formula, rsat)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use saturation, only: sat_mixrat_liq_api

  implicit none

  integer, intent(in) :: nz, ngrdcol, saturation_formula
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: p_in_Pa, T_in_K
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: rsat

  rsat = sat_mixrat_liq_api(nz, ngrdcol, stored_grid, p_in_Pa, T_in_K, &
    saturation_formula)

end subroutine f2py_sat_mixrat_liq_2d

subroutine f2py_rcm_sat_adj(thlm, rtm, p_in_Pa, exner, &
    saturation_formula, rcm_out)

  use clubb_precision, only: core_rknd
  use saturation, only: rcm_sat_adj

  implicit none

  real(core_rknd), intent(in) :: thlm, rtm, p_in_Pa, exner
  integer, intent(in) :: saturation_formula
  real(core_rknd), intent(out) :: rcm_out

  rcm_out = rcm_sat_adj(thlm, rtm, p_in_Pa, exner, saturation_formula)

end subroutine f2py_rcm_sat_adj

subroutine f2py_sat_mixrat_ice_2d(nz, ngrdcol, p_in_Pa, T_in_K, &
    saturation_formula, rsat)

  use clubb_precision, only: core_rknd
  use saturation, only: sat_mixrat_ice

  implicit none

  integer, intent(in) :: nz, ngrdcol, saturation_formula
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: p_in_Pa, T_in_K
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: rsat

  rsat = sat_mixrat_ice(nz, ngrdcol, p_in_Pa, T_in_K, saturation_formula)

end subroutine f2py_sat_mixrat_ice_2d
