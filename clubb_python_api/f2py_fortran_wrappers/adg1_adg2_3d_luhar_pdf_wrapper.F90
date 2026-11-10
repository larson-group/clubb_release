! adg1_adg2_3d_luhar_pdf_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module adg1_adg2_3d_luhar_pdf

subroutine f2py_calc_luhar_params( nz, ngrdcol, skx, wpxp, xp2, x_tol_sqd, & ! In
                                   mixt_frac, big_m, small_m ) ! Out

  ! Description:
  !   Call calc_luhar_params with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use adg1_adg2_3d_luhar_pdf, only: calc_luhar_params

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: nz ! Number of vertical levels [-]
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  ! Skewness of x ( <x'^3> / <x'2>^(3/2) )      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: skx
  ! Covariance of w and x                       [m/s (x units)]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wpxp
  ! Variance of x                               [(x units)^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: xp2
  real(core_rknd), intent(in) :: x_tol_sqd ! Tolerance value of x squared           [(x units)^2]

  !--------------------------- Output Variables ---------------------------
  ! Mixture fraction                       [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mixt_frac
  ! Luhar's M                              [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: big_m
  ! Luhar's m                              [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: small_m

  !--------------------------- Begin Code ---------------------------

  call calc_luhar_params(skx, wpxp, xp2, x_tol_sqd, & ! In
    mixt_frac, big_m, small_m) ! Out

end subroutine f2py_calc_luhar_params

subroutine f2py_close_luhar_pdf( nz, ngrdcol, xm, xp2, mixt_frac, small_m, wpxp, & ! In
                                 x_tol_sqd, & ! In
                                 sigma_sqd_x_1, sigma_sqd_x_2, varnce_x_1, varnce_x_2, & ! Out
                                 x_1_n, x_2_n, x_1, x_2 ) ! Out

  ! Description:
  !   Call close_luhar_pdf with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use adg1_adg2_3d_luhar_pdf, only: close_luhar_pdf

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: nz ! Number of vertical levels [-]
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  ! Mean (overall) of x, <x>                    [(x units)]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: xm
  ! Variance (overall) of x, <x'^2>             [(x units)^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: xp2
  ! Mixture fraction                            [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: mixt_frac
  ! Luhar's small m                             [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: small_m
  ! Covariance of w and x, <w'x'>               [m/s (x units)]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wpxp
  ! Tolerance value of x squared                [(x units)^2]
  real(core_rknd), intent(in) :: x_tol_sqd

  !--------------------------- Output Variables ---------------------------
  ! Normalized width parameter of x (1st PDF component) [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_sqd_x_1
  ! Normalized width parameter of x (1st PDF component) [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_sqd_x_2
  ! Variance of x (1st PDF component)         [(x units)^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: varnce_x_1
  ! Variance of x (2nd PDF component)         [(x units)^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: varnce_x_2
  ! Normalized mean of x (1st PDF component)            [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: x_1_n
  ! Normalized mean of x (2nd PDF component)            [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: x_2_n
  ! Mean of x (1st PDF component)               [(x units)]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: x_1
  ! Mean of x (2nd PDF component)               [(x units)]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: x_2

  !--------------------------- Begin Code ---------------------------

  call close_luhar_pdf( xm, xp2, mixt_frac, small_m, & ! In
                        wpxp, x_tol_sqd, & ! In
                        sigma_sqd_x_1, sigma_sqd_x_2, & ! Out
                        varnce_x_1, varnce_x_2, & ! Out
                        x_1_n, x_2_n, & ! Out
                        x_1, x_2 ) ! Out

end subroutine f2py_close_luhar_pdf

subroutine f2py_adg1_w_closure(nz, ngrdcol, wm, wp2, skw, sigma_sqd_w, sqrt_wp2, &
    mixt_frac_max_mag, w_1, w_2, w_1_n, w_2_n, varnce_w_1, varnce_w_2, mixt_frac)

  use clubb_precision, only: core_rknd
  use adg1_adg2_3d_luhar_pdf, only: adg1_w_closure
  use derived_type_storage, only: stored_err_info

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wm, wp2, skw, sigma_sqd_w, sqrt_wp2
  real(core_rknd), dimension(ngrdcol), intent(in) :: mixt_frac_max_mag
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
  real(core_rknd), dimension(ngrdcol), intent(in) :: mixt_frac_max_mag
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

subroutine f2py_luhar_3d_pdf_driver( nz, ngrdcol, wm, rtm, thlm, wp2, rtp2, thlp2, & ! In
                                     skw, skrt, skthl, wprtp, wpthlp, & ! In
                                     w_1, w_2, rt_1, rt_2, thl_1, thl_2, varnce_w_1, & ! Out
                                     varnce_w_2, varnce_rt_1, varnce_rt_2, & ! Out
                                     varnce_thl_1, varnce_thl_2, mixt_frac ) ! Out

  ! Description:
  !   Call luhar_3d_pdf_driver with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use adg1_adg2_3d_luhar_pdf, only: luhar_3d_pdf_driver

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: nz ! Number of vertical levels [-]
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  ! Mean of w-wind component (vertical velocity)  [m/s]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wm
  ! Mean of total water mixing ratio              [kg/kg]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: rtm
  ! Mean of liquid water potential temperature    [K]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: thlm
  ! Variance of w (overall)                       [m^2/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wp2
  ! Variance of r_t (overall)                     [(kg/kg)^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: rtp2
  ! Variance of th_l (overall)                    [K^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: thlp2
  ! Skewness of w                                 [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: skw
  ! Skewness of rt                                [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: skrt
  ! Skewness of th_l                              [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: skthl
  ! Covariance of w and r_t                       [(kg/kg)(m/s)]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wprtp
  ! Covariance of w and th_l                      [K(m/s)]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wpthlp

  !--------------------------- Output Variables ---------------------------
  ! Mean of w (1st PDF component)                      [m/s]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: w_1
  ! Mean of w (2nd PDF component)                      [m/s]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: w_2
  ! Mean of r_t (1st PDF component)                  [kg/kg]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: rt_1
  ! Mean of r_t (2nd PDF component)                  [kg/kg]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: rt_2
  ! Mean of th_l (1st PDF component)                     [K]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: thl_1
  ! Mean of th_l (2nd PDF component)                     [K]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: thl_2
  ! Variance of w (1st PDF component)              [m^2/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: varnce_w_1
  ! Variance of w (2nd PDF component)              [m^2/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: varnce_w_2
  ! Variance of r_t (1st PDF component)          [kg^2/kg^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: varnce_rt_1
  ! Variance of r_t (2nd PDF component)          [kg^2/kg^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: varnce_rt_2
  ! Variance of th_l (1st PDF component)               [K^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: varnce_thl_1
  ! Variance of th_l (2nd PDF component)               [K^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: varnce_thl_2
  ! Mixture fraction (weight of 1st PDF component)       [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mixt_frac

  !--------------------------- Begin Code ---------------------------

  call luhar_3d_pdf_driver( nz, ngrdcol, & ! In
                            wm, rtm, thlm, wp2, rtp2, thlp2, & ! In
                            skw, skrt, skthl, wprtp, wpthlp, & ! In
                            w_1, w_2, & ! Out
                            rt_1, rt_2, & ! Out
                            thl_1, thl_2, & ! Out
                            varnce_w_1, varnce_w_2, & ! Out
                            varnce_rt_1, varnce_rt_2, & ! Out
                            varnce_thl_1, varnce_thl_2, & ! Out
                            mixt_frac ) ! Out

end subroutine f2py_luhar_3d_pdf_driver
