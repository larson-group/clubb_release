! Use the complete native interpolation module as the numerical oracle.
subroutine f2k_interpolation_scalar_c(n, data, levels, output) bind(C)
  ! Description: Evaluate both linear formulas and the monotone cubic on supplied
  ! stencils. Relative level indices choose endpoint/interior/extrapolation paths.
  use, intrinsic :: iso_c_binding, only: c_int
  use clubb_precision, only: core_rknd
  use interpolation, only: lin_interpolate_two_points, linear_interp_factor, mono_cubic_interp
  implicit none

  !--------------------------- Input Variables ---------------------------
  integer(c_int), value, intent(in) :: n ! Number of stencils [-]
  real( kind = core_rknd ), intent(in) :: data(9,n) ! Target, four heights and four field values
  integer(c_int), intent(in) :: levels(4,n) ! Native branch-selecting level indices [-]

  !--------------------------- Output Variables ---------------------------
  real( kind = core_rknd ), intent(out) :: output(3,n) ! Linear, factor and cubic values

  !--------------------------- Local Variables ---------------------------
  integer :: j ! Stencil index [-]

  !--------------------------- Begin Code ---------------------------
  do j = 1, n
    output(1,j) = lin_interpolate_two_points(data(1,j), data(4,j), data(3,j), &
                                             data(8,j), data(7,j))
    output(2,j) = linear_interp_factor((data(1,j)-data(3,j))/(data(4,j)-data(3,j)), &
                                        data(8,j), data(7,j))
    output(3,j) = mono_cubic_interp(data(1,j), levels(1,j), levels(2,j), &
                                     levels(3,j), levels(4,j), &
                                     data(2,j), data(3,j), data(4,j), data(5,j), &
                                     data(6,j), data(7,j), data(8,j), data(9,j))
  end do
end subroutine f2k_interpolation_scalar_c

subroutine f2k_interpolation_profile_c(nout, nsrc, grid_out, grid_src, values, &
                                       pressure_out, pressure_src, output, search, on_grid) bind(C)
  ! Description: Compare zero-fill and clamped profiles, pressure-grid negation,
  ! upper-bracket searches and the scalar on-grid interface on valid coordinates.
  use, intrinsic :: iso_c_binding, only: c_int
  use clubb_precision, only: core_rknd
  use error_code, only: set_clubb_debug_level_api
  use interpolation, only: binary_search, plinterp_fnc, zlinterp_fnc, &
                           lin_interp_between_grids, lin_interpolate_on_grid_api
  implicit none

  !--------------------------- Input Variables ---------------------------
  integer(c_int), value, intent(in) :: nout, nsrc ! Target and source profile sizes [-]
  real( kind = core_rknd ), intent(in) :: grid_out(nout), grid_src(nsrc) ! Heights [m]
  real( kind = core_rknd ), intent(in) :: values(nsrc) ! Field values [units vary]
  real( kind = core_rknd ), intent(in) :: pressure_out(nout), pressure_src(nsrc) ! Pressures [Pa]

  !--------------------------- Output Variables ---------------------------
  real( kind = core_rknd ), intent(out) :: output(nout,3) ! Height, pressure and clamped values
  integer(c_int), intent(out) :: search(nout) ! Native one-based upper brackets or -1
  real( kind = core_rknd ), intent(out) :: on_grid(nsrc) ! Exact source-point scalar calls

  !--------------------------- Local Variables ---------------------------
  real( kind = core_rknd ) :: xlist(nsrc), tlist(nsrc) ! Legacy inout inputs, unchanged by routine
  integer :: i ! Profile index [-]

  !--------------------------- Begin Code ---------------------------
  call set_clubb_debug_level_api(0) ! In
  output(:,1) = zlinterp_fnc(nout, nsrc, grid_out, grid_src, values)
  output(:,2) = plinterp_fnc(nout, nsrc, pressure_out, pressure_src, values)
  output(:,3) = lin_interp_between_grids(nout, nsrc, grid_out, grid_src, values)
  do i = 1, nout
    search(i) = binary_search(nsrc, grid_src, grid_out(i))
  end do
  if (nsrc > 1) then
    xlist = grid_src
    tlist = values
    do i = 1, nsrc
      call lin_interpolate_on_grid_api(nsrc, & ! In
                                        xlist, tlist, & ! InOut
                                        grid_src(i), & ! In
                                        on_grid(i)) ! InOut
    end do
  end if
end subroutine f2k_interpolation_profile_c

subroutine f2k_interpolation_equal_error_c() bind(C)
  ! Description: Confirm the canonical host routine rejects coincident heights.
  use clubb_precision, only: core_rknd
  use interpolation, only: lin_interpolate_two_points
  implicit none
  real( kind = core_rknd ) :: value ! Unreachable result [units vary]
  value = lin_interpolate_two_points(0._core_rknd, 1._core_rknd, 1._core_rknd, &
                                      2._core_rknd, 3._core_rknd)
end subroutine f2k_interpolation_equal_error_c
