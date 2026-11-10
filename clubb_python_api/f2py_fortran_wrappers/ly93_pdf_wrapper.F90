! ly93_pdf_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module ly93_pdf

subroutine f2py_calc_params_ly93( nz, ngrdcol, xm, xp2, skx, mixt_frac, & ! In
                                  mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd ) ! Out

  ! Description:
  !   Call calc_params_ly93 with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use LY93_pdf, only: calc_params_LY93

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: nz ! Number of vertical levels [-]
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  ! Mean of x (overall)        [units vary]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: xm
  ! Variance of x (overall)    [(units vary)^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: xp2
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: skx ! Skewness of x (overall)    [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: mixt_frac ! Mixture fraction           [-]

  !--------------------------- Output Variables ---------------------------
  ! Mean of x (1st PDF component)        [units vary]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mu_x_1
  ! Mean of x (2nd PDF component)        [units vary]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mu_x_2
  ! Variance of x (1st PDF component)    [(units vary)^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_x_1_sqd
  ! Variance of x (2nd PDF component)    [(units vary)^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_x_2_sqd

  !--------------------------- Begin Code ---------------------------

  call calc_params_LY93( nz, ngrdcol, xm, xp2, skx, mixt_frac, & ! In
                         mu_x_1, mu_x_2, sigma_x_1_sqd, & ! Out
                         sigma_x_2_sqd ) ! Out

end subroutine f2py_calc_params_ly93

subroutine f2py_ly93_driver( nz, ngrdcol, wm, rtm, thlm, wp2, rtp2, thlp2, skw, skrt, & ! In
                             skthl, & ! In
                             mu_w_1, mu_w_2, mu_rt_1, mu_rt_2, mu_thl_1, mu_thl_2, & ! Out
                             sigma_w_1_sqd, sigma_w_2_sqd, sigma_rt_1_sqd, & ! Out
                             sigma_rt_2_sqd, sigma_thl_1_sqd, sigma_thl_2_sqd, & ! Out
                             mixt_frac ) ! Out

  ! Description:
  !   Call ly93_driver with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use LY93_pdf, only: LY93_driver

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: nz ! Number of vertical levels [-]
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  ! Mean of w (overall)                      [m/s]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wm
  ! Mean of rt (overall)                     [kg/kg]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: rtm
  ! Mean of thl (overall)                    [K]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: thlm
  ! Variance of w (overall)                  [m^2/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wp2
  ! Variance of rt (overall)                 [kg^2/kg^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: rtp2
  ! Variance of thl (overall)                [K^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: thlp2
  ! Skewness of w (overall)                  [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: skw
  ! Skewness of rt (overall)                 [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: skrt
  ! Skewness of thl (overall)                [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: skthl

  !--------------------------- Output Variables ---------------------------
  ! Mean of w (1st PDF component)          [m/s]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mu_w_1
  ! Mean of w (2nd PDF component)          [m/s]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mu_w_2
  ! Mean of rt (1st PDF component)         [kg/kg]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mu_rt_1
  ! Mean of rt (2nd PDF component)         [kg/kg]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mu_rt_2
  ! Mean of thl (1st PDF component)        [K]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mu_thl_1
  ! Mean of thl (2nd PDF component)        [K]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mu_thl_2
  ! Variance of w (1st PDF component)      [m^2/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_w_1_sqd
  ! Variance of w (2nd PDF component)      [m^2/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_w_2_sqd
  ! Variance of rt (1st PDF component)     [m^2/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_rt_1_sqd
  ! Variance of rt (2nd PDF component)     [m^2/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_rt_2_sqd
  ! Variance of thl (1st PDF component)    [m^2/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_thl_1_sqd
  ! Variance of thl (2nd PDF component)    [m^2/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_thl_2_sqd
  ! Mixture fraction                       [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mixt_frac

  !--------------------------- Begin Code ---------------------------

  call LY93_driver( nz, ngrdcol, & ! In
                    wm, rtm, thlm, wp2, rtp2, thlp2, & ! In
                    skw, skrt, skthl, & ! In
                    mu_w_1, mu_w_2, & ! Out
                    mu_rt_1, mu_rt_2, & ! Out
                    mu_thl_1, mu_thl_2, & ! Out
                    sigma_w_1_sqd, sigma_w_2_sqd, & ! Out
                    sigma_rt_1_sqd, sigma_rt_2_sqd, & ! Out
                    sigma_thl_1_sqd, sigma_thl_2_sqd, & ! Out
                    mixt_frac ) ! Out

end subroutine f2py_ly93_driver
