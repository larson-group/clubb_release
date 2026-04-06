! diffusion_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module diffusion

subroutine f2py_diffusion_zt_lhs(nzm, nzt, ngrdcol, k_zm, k_zt, nu, &
    invrs_rho_ds_zt, rho_ds_zm, lhs)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use diffusion, only: diffusion_zt_lhs

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: k_zm, rho_ds_zm
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: k_zt, invrs_rho_ds_zt
  real(core_rknd), dimension(ngrdcol), intent(in) :: nu
  real(core_rknd), dimension(3, ngrdcol, nzt), intent(out) :: lhs

  call diffusion_zt_lhs(nzm, nzt, ngrdcol, stored_grid, k_zm, k_zt, nu, &
    invrs_rho_ds_zt, rho_ds_zm, lhs)

end subroutine f2py_diffusion_zt_lhs

subroutine f2py_diffusion_zm_lhs(nzm, nzt, ngrdcol, k_zt, k_zm, nu, &
    invrs_rho_ds_zm, rho_ds_zt, lhs)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use diffusion, only: diffusion_zm_lhs

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: k_zt, rho_ds_zt
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: k_zm, invrs_rho_ds_zm
  real(core_rknd), dimension(ngrdcol), intent(in) :: nu
  real(core_rknd), dimension(3, ngrdcol, nzm), intent(out) :: lhs

  call diffusion_zm_lhs(nzm, nzt, ngrdcol, stored_grid, k_zt, k_zm, nu, &
    invrs_rho_ds_zm, rho_ds_zt, lhs)

end subroutine f2py_diffusion_zm_lhs
