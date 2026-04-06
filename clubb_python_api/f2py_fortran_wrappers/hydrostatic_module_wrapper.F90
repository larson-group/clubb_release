! hydrostatic_module_wrapper.F90 — wrappers organized by source module

subroutine f2py_hydrostatic(ngrdcol, nzt, nzm, thvm, p_sfc, &
    p_in_Pa, p_in_Pa_zm, exner, exner_zm, rho, rho_zm)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use hydrostatic_module, only: hydrostatic

  implicit none

  integer, intent(in) :: ngrdcol, nzt, nzm
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: thvm
  real(core_rknd), dimension(ngrdcol), intent(in) :: p_sfc
  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: p_in_Pa, exner, rho
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: p_in_Pa_zm, exner_zm, rho_zm

  call hydrostatic(ngrdcol, stored_grid, thvm, p_sfc, &
    p_in_Pa, p_in_Pa_zm, exner, exner_zm, rho, rho_zm)

end subroutine f2py_hydrostatic

