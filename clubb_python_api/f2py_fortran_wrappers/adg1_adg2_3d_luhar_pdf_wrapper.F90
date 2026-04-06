! adg1_adg2_3d_luhar_pdf_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module adg1_adg2_3d_luhar_pdf

subroutine f2py_calc_luhar_params(nz, ngrdcol, skx, wpxp, xp2, x_tol_sqd, &
    mixt_frac, big_m, small_m)

  use clubb_precision, only: core_rknd
  use adg1_adg2_3d_luhar_pdf, only: calc_luhar_params

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: skx, wpxp, xp2
  real(core_rknd), intent(in) :: x_tol_sqd
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mixt_frac, big_m, small_m

  integer :: i

  do i = 1, ngrdcol
    call calc_luhar_params(nz, skx(i,:), wpxp(i,:), xp2(i,:), x_tol_sqd, &
      mixt_frac(i,:), big_m(i,:), small_m(i,:))
  end do

end subroutine f2py_calc_luhar_params

subroutine f2py_close_luhar_pdf(nz, ngrdcol, xm, xp2, mixt_frac, small_m, wpxp, &
    x_tol_sqd, sigma_sqd_x_1, sigma_sqd_x_2, varnce_x_1, varnce_x_2, x_1_n, &
    x_2_n, x_1, x_2)

  use clubb_precision, only: core_rknd
  use adg1_adg2_3d_luhar_pdf, only: close_luhar_pdf

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: xm, xp2, mixt_frac, small_m, wpxp
  real(core_rknd), intent(in) :: x_tol_sqd
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
    sigma_sqd_x_1, sigma_sqd_x_2, varnce_x_1, varnce_x_2, x_1_n, x_2_n, x_1, x_2

  integer :: i

  do i = 1, ngrdcol
    call close_luhar_pdf(nz, xm(i,:), xp2(i,:), mixt_frac(i,:), small_m(i,:), wpxp(i,:), &
      x_tol_sqd, sigma_sqd_x_1(i,:), sigma_sqd_x_2(i,:), varnce_x_1(i,:), &
      varnce_x_2(i,:), x_1_n(i,:), x_2_n(i,:), x_1(i,:), x_2(i,:))
  end do

end subroutine f2py_close_luhar_pdf

subroutine f2py_adg1_w_closure(nz, ngrdcol, wm, wp2, skw, sigma_sqd_w, sqrt_wp2, &
    mixt_frac_max_mag, w_1, w_2, w_1_n, w_2_n, varnce_w_1, varnce_w_2, mixt_frac)

  use clubb_precision, only: core_rknd
  use adg1_adg2_3d_luhar_pdf, only: adg1_w_closure
  use derived_type_storage, only: stored_err_info

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wm, wp2, skw, sigma_sqd_w, sqrt_wp2
  real(core_rknd), intent(in) :: mixt_frac_max_mag
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
    w_1, w_2, w_1_n, w_2_n, varnce_w_1, varnce_w_2, mixt_frac

  call adg1_w_closure(nz, ngrdcol, wm, wp2, skw, sigma_sqd_w, sqrt_wp2, mixt_frac_max_mag, stored_err_info, &
    w_1, w_2, w_1_n, w_2_n, varnce_w_1, varnce_w_2, mixt_frac)

end subroutine f2py_adg1_w_closure

subroutine f2py_adg2_pdf_driver(nz, ngrdcol, sclr_dim, sclr_dim_transport, sclr_tol, wm, rtm, thlm, &
    wp2, rtp2, thlp2, skw, wprtp, wpthlp, sqrt_wp2, beta, sclrm, sclrp2, wpsclrp, &
    l_scalar_calc, w_1, w_2, rt_1, rt_2, thl_1, thl_2, varnce_w_1, varnce_w_2, &
    varnce_rt_1, varnce_rt_2, varnce_thl_1, varnce_thl_2, mixt_frac, alpha_rt, &
    alpha_thl, sigma_sqd_w, sclr_1, sclr_2, varnce_sclr_1, varnce_sclr_2, alpha_sclr)

  use clubb_precision, only: core_rknd
  use adg1_adg2_3d_luhar_pdf, only: adg2_pdf_driver
  use derived_type_storage, only: stored_err_info

  implicit none

  integer, intent(in) :: nz, ngrdcol, sclr_dim, sclr_dim_transport
  logical, intent(in) :: l_scalar_calc
  real(core_rknd), dimension(sclr_dim_transport), intent(in) :: sclr_tol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    wm, rtm, thlm, wp2, rtp2, thlp2, skw, wprtp, wpthlp, sqrt_wp2
  real(core_rknd), dimension(ngrdcol), intent(in) :: beta
  real(core_rknd), dimension(ngrdcol, nz, sclr_dim_transport), intent(in) :: sclrm, sclrp2, wpsclrp
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
    w_1, w_2, rt_1, rt_2, thl_1, thl_2, varnce_w_1, varnce_w_2, &
    varnce_rt_1, varnce_rt_2, varnce_thl_1, varnce_thl_2, mixt_frac, alpha_rt, alpha_thl, sigma_sqd_w
  real(core_rknd), dimension(ngrdcol, nz, sclr_dim_transport), intent(out) :: &
    sclr_1, sclr_2, varnce_sclr_1, varnce_sclr_2, alpha_sclr

  call adg2_pdf_driver(nz, ngrdcol, sclr_dim, sclr_tol(1:sclr_dim), wm, rtm, thlm, wp2, rtp2, thlp2, skw, &
    wprtp, wpthlp, sqrt_wp2, beta, sclrm(:,:,1:sclr_dim), sclrp2(:,:,1:sclr_dim), &
    wpsclrp(:,:,1:sclr_dim), l_scalar_calc, stored_err_info, &
    w_1, w_2, rt_1, rt_2, thl_1, thl_2, varnce_w_1, varnce_w_2, varnce_rt_1, varnce_rt_2, &
    varnce_thl_1, varnce_thl_2, mixt_frac, alpha_rt, alpha_thl, sigma_sqd_w, &
    sclr_1(:,:,1:sclr_dim), sclr_2(:,:,1:sclr_dim), varnce_sclr_1(:,:,1:sclr_dim), &
    varnce_sclr_2(:,:,1:sclr_dim), alpha_sclr(:,:,1:sclr_dim))

end subroutine f2py_adg2_pdf_driver

subroutine f2py_adg1_pdf_driver(nz, ngrdcol, sclr_dim, sclr_dim_transport, sclr_tol, wm, rtm, thlm, um, vm, &
    wp2, rtp2, thlp2, up2, vp2, skw, wprtp, wpthlp, upwp, vpwp, sqrt_wp2, sigma_sqd_w, &
    beta, mixt_frac_max_mag, sclrm, sclrp2, wpsclrp, l_scalar_calc, w_1, w_2, rt_1, rt_2, &
    thl_1, thl_2, u_1, u_2, v_1, v_2, varnce_w_1, varnce_w_2, varnce_rt_1, varnce_rt_2, &
    varnce_thl_1, varnce_thl_2, varnce_u_1, varnce_u_2, varnce_v_1, varnce_v_2, mixt_frac, &
    alpha_rt, alpha_thl, alpha_u, alpha_v, sclr_1, sclr_2, varnce_sclr_1, varnce_sclr_2, alpha_sclr)

  use clubb_precision, only: core_rknd
  use adg1_adg2_3d_luhar_pdf, only: adg1_pdf_driver
  use derived_type_storage, only: stored_err_info

  implicit none

  integer, intent(in) :: nz, ngrdcol, sclr_dim, sclr_dim_transport
  logical, intent(in) :: l_scalar_calc
  real(core_rknd), dimension(sclr_dim_transport), intent(in) :: sclr_tol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    wm, rtm, thlm, um, vm, wp2, rtp2, thlp2, up2, vp2, skw, wprtp, wpthlp, upwp, vpwp, sqrt_wp2, sigma_sqd_w
  real(core_rknd), dimension(ngrdcol), intent(in) :: beta
  real(core_rknd), intent(in) :: mixt_frac_max_mag
  real(core_rknd), dimension(ngrdcol, nz, sclr_dim_transport), intent(in) :: sclrm, sclrp2, wpsclrp
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
    w_1, w_2, rt_1, rt_2, thl_1, thl_2, u_1, u_2, v_1, v_2, &
    varnce_w_1, varnce_w_2, varnce_rt_1, varnce_rt_2, varnce_thl_1, varnce_thl_2, &
    varnce_u_1, varnce_u_2, varnce_v_1, varnce_v_2, mixt_frac, alpha_rt, alpha_thl, alpha_u, alpha_v
  real(core_rknd), dimension(ngrdcol, nz, sclr_dim_transport), intent(out) :: &
    sclr_1, sclr_2, varnce_sclr_1, varnce_sclr_2, alpha_sclr

  call adg1_pdf_driver(nz, ngrdcol, sclr_dim, sclr_tol(1:sclr_dim), wm, rtm, thlm, um, vm, wp2, rtp2, thlp2, &
    up2, vp2, skw, wprtp, wpthlp, upwp, vpwp, sqrt_wp2, sigma_sqd_w, beta, mixt_frac_max_mag, &
    sclrm(:,:,1:sclr_dim), sclrp2(:,:,1:sclr_dim), wpsclrp(:,:,1:sclr_dim), l_scalar_calc, stored_err_info, &
    w_1, w_2, rt_1, rt_2, thl_1, thl_2, u_1, u_2, &
    v_1, v_2, varnce_w_1, varnce_w_2, varnce_rt_1, varnce_rt_2, varnce_thl_1, varnce_thl_2, &
    varnce_u_1, varnce_u_2, varnce_v_1, varnce_v_2, mixt_frac, alpha_rt, alpha_thl, alpha_u, &
    alpha_v, sclr_1(:,:,1:sclr_dim), sclr_2(:,:,1:sclr_dim), varnce_sclr_1(:,:,1:sclr_dim), &
    varnce_sclr_2(:,:,1:sclr_dim), alpha_sclr(:,:,1:sclr_dim))

end subroutine f2py_adg1_pdf_driver

subroutine f2py_luhar_3d_pdf_driver(nz, ngrdcol, wm, rtm, thlm, wp2, rtp2, thlp2, skw, skrt, &
    skthl, wprtp, wpthlp, w_1, w_2, rt_1, rt_2, thl_1, thl_2, varnce_w_1, varnce_w_2, &
    varnce_rt_1, varnce_rt_2, varnce_thl_1, varnce_thl_2, mixt_frac)

  use clubb_precision, only: core_rknd
  use adg1_adg2_3d_luhar_pdf, only: luhar_3d_pdf_driver

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    wm, rtm, thlm, wp2, rtp2, thlp2, skw, skrt, skthl, wprtp, wpthlp
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
    w_1, w_2, rt_1, rt_2, thl_1, thl_2, varnce_w_1, varnce_w_2, &
    varnce_rt_1, varnce_rt_2, varnce_thl_1, varnce_thl_2, mixt_frac

  integer :: i

  do i = 1, ngrdcol
    call luhar_3d_pdf_driver(nz, wm(i,:), rtm(i,:), thlm(i,:), wp2(i,:), rtp2(i,:), thlp2(i,:), &
      skw(i,:), skrt(i,:), skthl(i,:), wprtp(i,:), wpthlp(i,:), &
      w_1(i,:), w_2(i,:), rt_1(i,:), rt_2(i,:), thl_1(i,:), thl_2(i,:), &
      varnce_w_1(i,:), varnce_w_2(i,:), varnce_rt_1(i,:), varnce_rt_2(i,:), &
      varnce_thl_1(i,:), varnce_thl_2(i,:), mixt_frac(i,:))
  end do

end subroutine f2py_luhar_3d_pdf_driver
