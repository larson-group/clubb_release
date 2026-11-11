! Compare staggered-grid operators with current compiled native procedures.
subroutine f2k_grid_operators_c(nzm, ngrdcol, ascending, momentum, thermo, azm, azt, minimum, &
                                 output_m, output_t, profile_m, profile_t, &
                                 grid_output) bind(C)
  ! Description: Construct a native host grid with column-varying stretched
  ! heights, then reference every interpolation, smoothing and derivative overload.
  ! Return that same grid through the real fieldwise converter for the GPU calls.
  use, intrinsic :: iso_c_binding, only: c_int
  use clubb_precision, only: core_rknd
  use grid_class, only: grid, setup_grid_heights, zt2zm_api, zm2zt_api, &
                       zt2zm2zt, zm2zt2zm, ddzm, ddzt
  use err_info_type_module, only: err_info_type
  use f2k_derived_type_bridge, only: wire_field, from_native_grid
  implicit none

  !--------------------------- Input Variables ---------------------------
  integer(c_int), value, intent(in) :: nzm, ngrdcol, ascending ! Grid shape/direction [-]
  real( kind = core_rknd ), intent(in) :: momentum(ngrdcol,nzm) ! Momentum heights [m]
  real( kind = core_rknd ), intent(in) :: thermo(ngrdcol,nzm-1) ! Thermodynamic heights [m]
  real( kind = core_rknd ), intent(in) :: azm(ngrdcol,nzm), azt(ngrdcol,nzm-1) ! Field values
  real( kind = core_rknd ), value, intent(in) :: minimum ! Optional field floor

  !--------------------------- Output Variables ---------------------------
  real( kind = core_rknd ), intent(out) :: output_m(ngrdcol,nzm,5) ! Momentum operator results
  real( kind = core_rknd ), intent(out) :: output_t(ngrdcol,nzm-1,5) ! Thermo operator results
  real( kind = core_rknd ), intent(out) :: profile_m(nzm,5), profile_t(nzm-1,5) ! 1D/scalar results
  type(wire_field), intent(out) :: grid_output(16) ! Caller-owned converted native grid

  !--------------------------- Local Variables ---------------------------
  type(grid) :: gr ! Native setup and operator grid
  type(err_info_type) :: err_info ! Setup status
  integer :: k, nzt ! Vertical index and thermodynamic level count [-]
  real( kind = core_rknd ) :: deltaz(ngrdcol), zm_init(ngrdcol) ! Unused host-grid settings [m]

  !--------------------------- Begin Code ---------------------------
  nzt = nzm-1
  gr%nzm = nzm
  gr%nzt = nzt
  allocate(gr%zm(ngrdcol,nzm), gr%zt(ngrdcol,nzt), &
           gr%dzm(ngrdcol,nzm), gr%dzt(ngrdcol,nzt), &
           gr%invrs_dzm(ngrdcol,nzm), gr%invrs_dzt(ngrdcol,nzt), &
           gr%weights_zm2zt(ngrdcol,nzt,2), gr%weights_zt2zm(ngrdcol,nzm,2))
  allocate(err_info%err_code(ngrdcol))
  err_info%err_code = 0
  deltaz = 0._core_rknd
  zm_init = 0._core_rknd
  call setup_grid_heights(nzm, nzt, ngrdcol, .true., ascending /= 0, & ! In
                            1, deltaz, zm_init, momentum, thermo, & ! In
                            gr, err_info) ! InOut
  if (any(err_info%err_code /= 0)) error stop "native grid setup failed"

  output_m(:,:,1) = zt2zm_api(nzm, nzt, ngrdcol, gr, azt)
  output_m(:,:,2) = zt2zm_api(nzm, nzt, ngrdcol, gr, azt, minimum)
  output_m(:,:,3) = zm2zt2zm(nzm, nzt, ngrdcol, gr, azm)
  output_m(:,:,4) = zm2zt2zm(nzm, nzt, ngrdcol, gr, azm, minimum)
  output_m(:,:,5) = ddzt(nzm, nzt, ngrdcol, gr, azt)
  output_t(:,:,1) = zm2zt_api(nzm, nzt, ngrdcol, gr, azm)
  output_t(:,:,2) = zm2zt_api(nzm, nzt, ngrdcol, gr, azm, minimum)
  output_t(:,:,3) = zt2zm2zt(nzm, nzt, ngrdcol, gr, azt)
  output_t(:,:,4) = zt2zm2zt(nzm, nzt, ngrdcol, gr, azt, minimum)
  output_t(:,:,5) = ddzm(nzm, nzt, ngrdcol, gr, azm)
  profile_m(:,1) = zt2zm_api(gr, azt(1,:))
  profile_m(:,2) = zt2zm_api(gr, azt(1,:), minimum)
  profile_m(:,5) = ddzt(gr, azt(1,:))
  profile_t(:,1) = zm2zt_api(gr, azm(1,:))
  profile_t(:,2) = zm2zt_api(gr, azm(1,:), minimum)
  profile_t(:,5) = ddzm(gr, azm(1,:))
  do k = 1, nzm
    profile_m(k,3) = zt2zm_api(gr, azt(1,:), k)
    profile_m(k,4) = zt2zm_api(gr, azt(1,:), k, minimum)
  end do
  do k = 1, nzt
    profile_t(k,3) = zm2zt_api(gr, azm(1,:), k)
    profile_t(k,4) = zm2zt_api(gr, azm(1,:), k, minimum)
  end do
  call from_native_grid(gr, grid_output)
end subroutine f2k_grid_operators_c

subroutine f2k_grid_flip_c(n, x, output) bind(C)
  ! Description: Reference profile reversal, including an empty input profile.
  use, intrinsic :: iso_c_binding, only: c_int
  use clubb_precision, only: core_rknd
  use grid_class, only: flip
  implicit none

  !--------------------------- Input Variables ---------------------------
  integer(c_int), value, intent(in) :: n ! Profile length [-]
  real( kind = core_rknd ), intent(in) :: x(n) ! Input field values

  !--------------------------- Output Variables ---------------------------
  real( kind = core_rknd ), intent(out) :: output(n) ! Reversed field values

  !--------------------------- Begin Code ---------------------------
  output = flip(x, n)
end subroutine f2k_grid_flip_c
