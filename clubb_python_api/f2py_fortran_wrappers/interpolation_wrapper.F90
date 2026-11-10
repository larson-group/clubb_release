! interpolation_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module interpolation

subroutine f2py_lin_interpolate_two_points(height_int, height_high, height_low, &
    var_high, var_low, var_interp)

  use clubb_precision, only: core_rknd
  use interpolation, only: lin_interpolate_two_points

  implicit none

  real(core_rknd), intent(in) :: height_int, height_high, height_low
  real(core_rknd), intent(in) :: var_high, var_low
  real(core_rknd), intent(out) :: var_interp

  var_interp = lin_interpolate_two_points(height_int, height_high, height_low, &
    var_high, var_low)

end subroutine f2py_lin_interpolate_two_points

subroutine f2py_mono_cubic_interp(z_in, km1, k00, kp1, kp2, &
    zm1, z00, zp1, zp2, fm1, f00, fp1, fp2, f_out)

  use clubb_precision, only: core_rknd
  use interpolation, only: mono_cubic_interp

  implicit none

  real(core_rknd), intent(in) :: z_in
  integer, intent(in) :: km1, k00, kp1, kp2
  real(core_rknd), intent(in) :: zm1, z00, zp1, zp2
  real(core_rknd), intent(in) :: fm1, f00, fp1, fp2
  real(core_rknd), intent(out) :: f_out

  f_out = mono_cubic_interp(z_in, km1, k00, kp1, kp2, &
    zm1, z00, zp1, zp2, fm1, f00, fp1, fp2)

end subroutine f2py_mono_cubic_interp

subroutine f2py_zlinterp( ngrdcol, dim_out, dim_src, grid_out, grid_src, var_src, & ! In
                          var_out ) ! Out

  ! Description:
  !   Call zlinterp with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use interpolation, only: zlinterp

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  integer, intent(in) :: dim_out ! Number of target levels [-]
  integer, intent(in) :: dim_src ! Number of source levels [-]
  real(core_rknd), dimension(ngrdcol, dim_out), intent(in) :: grid_out ! [m]
  real(core_rknd), dimension(ngrdcol, dim_src), intent(in) :: grid_src ! [m]
  real(core_rknd), dimension(ngrdcol, dim_src), intent(in) :: var_src ! [units vary]

  !--------------------------- Output Variables ---------------------------
  real(core_rknd), dimension(ngrdcol, dim_out), intent(out) :: var_out ! [units vary]

  !--------------------------- Begin Code ---------------------------

  var_out = zlinterp(ngrdcol, dim_out, dim_src, grid_out, grid_src, var_src)

end subroutine f2py_zlinterp
