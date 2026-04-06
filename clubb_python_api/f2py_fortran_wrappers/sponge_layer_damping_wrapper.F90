! sponge_layer_damping_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module sponge_layer_damping

subroutine f2py_sponge_damp_xm(nzm, nzt, ngrdcol, dt, zt, zm, xm_ref, xm, &
    tau_sponge_damp, sponge_layer_depth, xm_p)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use sponge_layer_damping, only: sponge_damp_xm, sponge_damp_profile

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), intent(in) :: dt
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: zt, xm_ref, xm, tau_sponge_damp
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: zm
  real(core_rknd), dimension(ngrdcol), intent(in) :: sponge_layer_depth
  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: xm_p

  type(sponge_damp_profile) :: profile
  integer :: i

  allocate( profile%tau_sponge_damp(nzt) )
  do i = 1, ngrdcol
    profile%tau_sponge_damp = tau_sponge_damp(i,:)
    profile%sponge_layer_depth = sponge_layer_depth(i)
    xm_p(i,:) = sponge_damp_xm(stored_grid, nzm, nzt, dt, zt(i,:), zm(i,:), &
      xm_ref(i,:), xm(i,:), profile)
  end do
  deallocate( profile%tau_sponge_damp )

end subroutine f2py_sponge_damp_xm

subroutine f2py_sponge_damp_xp2(nzm, ngrdcol, dt, zm, xp2, x_tol_sqd, &
    tau_sponge_damp, sponge_layer_depth, xp2_damped)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use sponge_layer_damping, only: sponge_damp_xp2, sponge_damp_profile

  implicit none

  integer, intent(in) :: nzm, ngrdcol
  real(core_rknd), intent(in) :: dt
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: zm, xp2, tau_sponge_damp
  real(core_rknd), dimension(ngrdcol), intent(in) :: x_tol_sqd, sponge_layer_depth
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: xp2_damped

  type(sponge_damp_profile) :: profile
  integer :: i

  allocate( profile%tau_sponge_damp(nzm) )
  do i = 1, ngrdcol
    profile%tau_sponge_damp = tau_sponge_damp(i,:)
    profile%sponge_layer_depth = sponge_layer_depth(i)
    xp2_damped(i,:) = sponge_damp_xp2(stored_grid, nzm, dt, zm(i,:), xp2(i,:), &
      x_tol_sqd(i), profile)
  end do
  deallocate( profile%tau_sponge_damp )

end subroutine f2py_sponge_damp_xp2

subroutine f2py_sponge_damp_xp3(nzm, nzt, ngrdcol, dt, z, zm, xp3, &
    tau_sponge_damp, sponge_layer_depth, xp3_damped)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use sponge_layer_damping, only: sponge_damp_xp3, sponge_damp_profile

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), intent(in) :: dt
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: z, xp3, tau_sponge_damp
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: zm
  real(core_rknd), dimension(ngrdcol), intent(in) :: sponge_layer_depth
  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: xp3_damped

  type(sponge_damp_profile) :: profile
  integer :: i

  allocate( profile%tau_sponge_damp(nzt) )
  do i = 1, ngrdcol
    profile%tau_sponge_damp = tau_sponge_damp(i,:)
    profile%sponge_layer_depth = sponge_layer_depth(i)
    xp3_damped(i,:) = sponge_damp_xp3(stored_grid, nzm, nzt, dt, z(i,:), zm(i,:), &
      xp3(i,:), profile)
  end do
  deallocate( profile%tau_sponge_damp )

end subroutine f2py_sponge_damp_xp3

