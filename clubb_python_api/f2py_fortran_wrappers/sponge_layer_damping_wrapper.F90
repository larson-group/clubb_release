! sponge_layer_damping_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module sponge_layer_damping

subroutine f2py_sponge_damp_xm(nzm, nzt, dt, zt, zm, xm_ref, xm, &
    tau_sponge_damp, sponge_layer_depth, xm_p)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use sponge_layer_damping, only: sponge_damp_xm, sponge_damp_profile

  implicit none

  integer, intent(in) :: nzm, nzt
  real(core_rknd), intent(in) :: dt
  real(core_rknd), dimension(nzt), intent(in) :: zt, xm_ref, xm, tau_sponge_damp
  real(core_rknd), dimension(nzm), intent(in) :: zm
  real(core_rknd), intent(in) :: sponge_layer_depth
  real(core_rknd), dimension(nzt), intent(out) :: xm_p

  type(sponge_damp_profile) :: profile

  allocate( profile%tau_sponge_damp(nzt) )
  profile%tau_sponge_damp = tau_sponge_damp
  profile%sponge_layer_depth = sponge_layer_depth
  xm_p = sponge_damp_xm(stored_grid, nzm, nzt, dt, zt, zm, xm_ref, xm, profile)
  deallocate( profile%tau_sponge_damp )

end subroutine f2py_sponge_damp_xm

subroutine f2py_sponge_damp_xp2(nzm, dt, zm, xp2, x_tol_sqd, &
    tau_sponge_damp, sponge_layer_depth, xp2_damped)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use sponge_layer_damping, only: sponge_damp_xp2, sponge_damp_profile

  implicit none

  integer, intent(in) :: nzm
  real(core_rknd), intent(in) :: dt
  real(core_rknd), dimension(nzm), intent(in) :: zm, xp2, tau_sponge_damp
  real(core_rknd), intent(in) :: x_tol_sqd, sponge_layer_depth
  real(core_rknd), dimension(nzm), intent(out) :: xp2_damped

  type(sponge_damp_profile) :: profile

  allocate( profile%tau_sponge_damp(nzm) )
  profile%tau_sponge_damp = tau_sponge_damp
  profile%sponge_layer_depth = sponge_layer_depth
  xp2_damped = sponge_damp_xp2(stored_grid, nzm, dt, zm, xp2, x_tol_sqd, profile)
  deallocate( profile%tau_sponge_damp )

end subroutine f2py_sponge_damp_xp2

subroutine f2py_sponge_damp_xp3(nzm, nzt, dt, z, zm, xp3, &
    tau_sponge_damp, sponge_layer_depth, xp3_damped)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use sponge_layer_damping, only: sponge_damp_xp3, sponge_damp_profile

  implicit none

  integer, intent(in) :: nzm, nzt
  real(core_rknd), intent(in) :: dt
  real(core_rknd), dimension(nzt), intent(in) :: z, xp3, tau_sponge_damp
  real(core_rknd), dimension(nzm), intent(in) :: zm
  real(core_rknd), intent(in) :: sponge_layer_depth
  real(core_rknd), dimension(nzt), intent(out) :: xp3_damped

  type(sponge_damp_profile) :: profile

  allocate( profile%tau_sponge_damp(nzt) )
  profile%tau_sponge_damp = tau_sponge_damp
  profile%sponge_layer_depth = sponge_layer_depth
  xp3_damped = sponge_damp_xp3(stored_grid, nzm, nzt, dt, z, zm, xp3, profile)
  deallocate( profile%tau_sponge_damp )

end subroutine f2py_sponge_damp_xp3
