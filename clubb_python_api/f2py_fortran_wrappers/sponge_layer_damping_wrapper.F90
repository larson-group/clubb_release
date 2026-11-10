! sponge_layer_damping_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module sponge_layer_damping

subroutine f2py_sponge_damp_xm( ngrdcol, nzm, nzt, dt, zt, zm, xm_ref, xm, & ! In
                                tau_sponge_damp, sponge_layer_depth, & ! In
                                xm_p ) ! Out

  ! Description:
  !   Call sponge_damp_xm with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use sponge_layer_damping, only: sponge_damp_xm, sponge_damp_profile

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  integer, intent(in) :: nzm ! Number of momentum levels [-]
  integer, intent(in) :: nzt ! Number of thermodynamic levels [-]
  real(core_rknd), intent(in) :: dt ! Model Timestep  [s]
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: zt ! Thermodynamic grid heights [m]
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: zm ! Momentum grid heights [m]
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: xm_ref ! Reference mean field [units vary]
  ! Mean field before damping [units vary]
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: xm
  ! Column/level damping time scales [s]
  real(core_rknd), dimension(ngrdcol,nzt), intent(in) :: tau_sponge_damp
  real(core_rknd), dimension(ngrdcol), intent(in) :: sponge_layer_depth ! Sponge layer depths [m]

  !--------------------------- Output Variables ---------------------------
  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: xm_p ! Damped value of xm  [units_vary]

  !--------------------------- Local Variables ---------------------------
  type(sponge_damp_profile) :: profile ! Reconstructed column damping profiles

  !--------------------------- Begin Code ---------------------------

  allocate( profile%tau_sponge_damp(ngrdcol,nzt), profile%sponge_layer_depth(ngrdcol) )
  profile%tau_sponge_damp = tau_sponge_damp
  profile%sponge_layer_depth = sponge_layer_depth
  xm_p = sponge_damp_xm(ngrdcol, nzm, nzt, stored_grid, dt, zt, zm, xm_ref, xm, profile)
  deallocate( profile%tau_sponge_damp, profile%sponge_layer_depth )

end subroutine f2py_sponge_damp_xm

subroutine f2py_sponge_damp_xp2( ngrdcol, nzm, dt, zm, xp2, x_tol_sqd, & ! In
                                 tau_sponge_damp, sponge_layer_depth, & ! In
                                 xp2_damped ) ! Out

  ! Description:
  !   Call sponge_damp_xp2 with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use sponge_layer_damping, only: sponge_damp_xp2, sponge_damp_profile

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  integer, intent(in) :: nzm ! Number of momentum levels [-]
  real(core_rknd), intent(in) :: dt ! Model Timestep  [s]
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: zm ! Momentum grid heights [m]
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: xp2 ! Variance before damping [units vary]
  real(core_rknd), intent(in) :: x_tol_sqd ! Minimum variance [units vary]
  ! Column/level damping time scales [s]
  real(core_rknd), dimension(ngrdcol,nzm), intent(in) :: tau_sponge_damp
  real(core_rknd), dimension(ngrdcol), intent(in) :: sponge_layer_depth ! Sponge layer depths [m]

  !--------------------------- Output Variables ---------------------------
  ! Variance of x, <x'^2>, after damping   [units vary]
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: xp2_damped

  !--------------------------- Local Variables ---------------------------
  type(sponge_damp_profile) :: profile ! Reconstructed column damping profiles

  !--------------------------- Begin Code ---------------------------

  allocate( profile%tau_sponge_damp(ngrdcol,nzm), profile%sponge_layer_depth(ngrdcol) )
  profile%tau_sponge_damp = tau_sponge_damp
  profile%sponge_layer_depth = sponge_layer_depth
  xp2_damped = sponge_damp_xp2(ngrdcol, nzm, stored_grid, dt, zm, xp2, x_tol_sqd, profile)
  deallocate( profile%tau_sponge_damp, profile%sponge_layer_depth )

end subroutine f2py_sponge_damp_xp2

subroutine f2py_sponge_damp_xp3( ngrdcol, nzm, nzt, dt, z, zm, xp3, tau_sponge_damp, & ! In
                                 sponge_layer_depth, & ! In
                                 xp3_damped ) ! Out

  ! Description:
  !   Call sponge_damp_xp3 with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use sponge_layer_damping, only: sponge_damp_xp3, sponge_damp_profile

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  integer, intent(in) :: nzm ! Number of momentum levels [-]
  integer, intent(in) :: nzt ! Number of thermodynamic levels [-]
  real(core_rknd), intent(in) :: dt ! Model Timestep  [s]
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: z ! Thermodynamic grid heights [m]
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: zm ! Momentum grid heights [m]
  ! Third moment before damping [units vary]
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: xp3
  ! Column/level damping time scales [s]
  real(core_rknd), dimension(ngrdcol,nzt), intent(in) :: tau_sponge_damp
  real(core_rknd), dimension(ngrdcol), intent(in) :: sponge_layer_depth ! Sponge layer depths [m]

  !--------------------------- Output Variables ---------------------------
  ! <x'^3> after damping   [units vary]
  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: xp3_damped

  !--------------------------- Local Variables ---------------------------
  type(sponge_damp_profile) :: profile ! Reconstructed column damping profiles

  !--------------------------- Begin Code ---------------------------

  allocate( profile%tau_sponge_damp(ngrdcol,nzt), profile%sponge_layer_depth(ngrdcol) )
  profile%tau_sponge_damp = tau_sponge_damp
  profile%sponge_layer_depth = sponge_layer_depth
  xp3_damped = sponge_damp_xp3(ngrdcol, nzm, nzt, stored_grid, dt, z, zm, xp3, profile)
  deallocate( profile%tau_sponge_damp, profile%sponge_layer_depth )

end subroutine f2py_sponge_damp_xp3
