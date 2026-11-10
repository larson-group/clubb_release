! new_tsdadg_pdf_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module new_tsdadg_pdf

subroutine f2py_calc_l_x_skx_fnc(skx, sgn_wpxp, small_l_x_1, small_l_x_2, &
    big_l_x_1, big_l_x_2)

  use clubb_precision, only: core_rknd
  use new_tsdadg_pdf, only: calc_L_x_Skx_fnc

  implicit none

  real(core_rknd), intent(in) :: skx, sgn_wpxp, small_l_x_1, small_l_x_2
  real(core_rknd), intent(out) :: big_l_x_1, big_l_x_2

  call calc_L_x_Skx_fnc(skx, sgn_wpxp, small_l_x_1, small_l_x_2, big_l_x_1, big_l_x_2)

end subroutine f2py_calc_l_x_skx_fnc

subroutine f2py_calc_setter_parameters_tsdadg(xm, xp2, skx, sgn_wpxp, &
    big_l_x_1, big_l_x_2, mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, &
    mixt_frac, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd)

  use clubb_precision, only: core_rknd
  use new_tsdadg_pdf, only: calc_setter_parameters

  implicit none

  real(core_rknd), intent(in) :: xm, xp2, skx, sgn_wpxp, big_l_x_1, big_l_x_2
  real(core_rknd), intent(out) :: &
    mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, mixt_frac, &
    coef_sigma_x_1_sqd, coef_sigma_x_2_sqd

  call calc_setter_parameters( xm, xp2, skx, sgn_wpxp, & ! In
                               big_l_x_1, big_l_x_2, & ! In
                               mu_x_1, mu_x_2, & ! Out
                               sigma_x_1_sqd, sigma_x_2_sqd, & ! Out
                               mixt_frac, & ! Out
                               coef_sigma_x_1_sqd, & ! Out
                               coef_sigma_x_2_sqd ) ! Out

end subroutine f2py_calc_setter_parameters_tsdadg

subroutine f2py_tsdadg_pdf_driver( nz, ngrdcol, wm, rtm, thlm, wp2, rtp2, thlp2, skw, & ! In
                                   skrt, skthl, wprtp, wpthlp, & ! In
                                   mu_w_1, mu_w_2, mu_rt_1, mu_rt_2, mu_thl_1, & ! Out
                                   mu_thl_2, sigma_w_1_sqd, sigma_w_2_sqd, & ! Out
                                   sigma_rt_1_sqd, sigma_rt_2_sqd, sigma_thl_1_sqd, & ! Out
                                   sigma_thl_2_sqd, mixt_frac ) ! Out

  ! Description:
  !   Call tsdadg_pdf_driver with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use new_tsdadg_pdf, only: tsdadg_pdf_driver

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: nz ! Number of vertical levels [-]
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  ! Mean of w (overall)                [m/s]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wm
  ! Mean of rt (overall)               [kg/kg]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: rtm
  ! Mean of thl (overall)              [K]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: thlm
  ! Variance of w (overall)            [m^2/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wp2
  ! Variance of rt (overall)           [kg^2/kg^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: rtp2
  ! Variance of thl (overall)          [K^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: thlp2
  ! Skewness of w (overall)            [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: skw
  ! Skewness of rt (overall)           [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: skrt
  ! Skewness of thl (overall)          [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: skthl
  ! Covariance of w and rt (overall)   [(m/s)kg/kg]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wprtp
  ! Covariance of w and thl (overall)  [(m/s)K]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wpthlp

  !--------------------------- Output Variables ---------------------------
  ! Mean of w (1st PDF component)        [m/s]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mu_w_1
  ! Mean of w (2nd PDF component)        [m/s]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mu_w_2
  ! Mean of rt (1st PDF component)       [kg/kg]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mu_rt_1
  ! Mean of rt (2nd PDF component)       [kg/kg]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mu_rt_2
  ! Mean of thl (1st PDF component)      [K]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mu_thl_1
  ! Mean of thl (2nd PDF component)      [K]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mu_thl_2
  ! Variance of w (1st PDF component)    [m^2/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_w_1_sqd
  ! Variance of w (2nd PDF component)    [m^2/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_w_2_sqd
  ! Variance of rt (1st PDF component)   [kg^2/kg^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_rt_1_sqd
  ! Variance of rt (2nd PDF component)   [kg^2/kg^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_rt_2_sqd
  ! Variance of thl (1st PDF component)  [K^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_thl_1_sqd
  ! Variance of thl (2nd PDF component)  [K^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_thl_2_sqd
  ! Mixture fraction                     [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mixt_frac

  !--------------------------- Begin Code ---------------------------

  call tsdadg_pdf_driver( nz, ngrdcol, & ! In
                          wm, rtm, thlm, wp2, rtp2, thlp2, & ! In
                          skw, skrt, skthl, wprtp, wpthlp, & ! In
                          mu_w_1, mu_w_2, & ! Out
                          mu_rt_1, mu_rt_2, & ! Out
                          mu_thl_1, mu_thl_2, & ! Out
                          sigma_w_1_sqd, sigma_w_2_sqd, & ! Out
                          sigma_rt_1_sqd, sigma_rt_2_sqd, & ! Out
                          sigma_thl_1_sqd, sigma_thl_2_sqd, & ! Out
                          mixt_frac ) ! Out

end subroutine f2py_tsdadg_pdf_driver
