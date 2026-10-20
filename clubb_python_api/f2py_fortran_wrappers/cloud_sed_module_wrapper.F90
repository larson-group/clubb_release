! cloud_sed_module_wrapper.F90 — wrappers for module cloud_sed_module (src/Microphys)

subroutine f2py_cloud_drop_sed(nzm, nzt, ngrdcol, rcm, Ncm, rho_zm, rho, exner, &
    sigma_g, rcm_mc, thlm_mc)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid, stored_stats
  use cloud_sed_module, only: cloud_drop_sed

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: rcm, Ncm, rho, exner
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: rho_zm
  real(core_rknd), intent(in) :: sigma_g
  real(core_rknd), dimension(ngrdcol, nzt), intent(inout) :: rcm_mc, thlm_mc

  call cloud_drop_sed(stored_grid, ngrdcol, rcm, Ncm, rho_zm, rho, exner, sigma_g, &
    stored_stats, rcm_mc, thlm_mc)

end subroutine f2py_cloud_drop_sed
