! new_pdf_main_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module new_pdf_main

subroutine f2py_new_pdf_driver(nz, ngrdcol, wm, rtm, thlm, wp2, rtp2, thlp2, skw, &
    wprtp, wpthlp, rtpthlp, clubb_params, skrt, skthl, mu_w_1, mu_w_2, mu_rt_1, &
    mu_rt_2, mu_thl_1, mu_thl_2, sigma_w_1_sqd, sigma_w_2_sqd, sigma_rt_1_sqd, &
    sigma_rt_2_sqd, sigma_thl_1_sqd, sigma_thl_2_sqd, mixt_frac)

  ! Description:
  !   Call new_pdf_driver with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use new_pdf_main, only: new_pdf_driver
  use derived_type_storage, only: stored_pdf_implicit_coefs_terms

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: nz ! Number of vertical levels [-]
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  ! Mean of w (overall)                 [m/s]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wm
  ! Mean of rt (overall)                [kg/kg]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: rtm
  ! Mean of thl (overall)               [K]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: thlm
  ! Variance of w (overall)             [m^2/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wp2
  ! Variance of rt (overall)            [kg^2/kg^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: rtp2
  ! Variance of thl (overall)           [K^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: thlp2
  ! Skewness of w (overall)             [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: skw
  ! Covariance of w and rt (overall)    [(m/s)kg/kg]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wprtp
  ! Covariance of w and thl (overall)   [(m/s)K]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wpthlp
  ! Covariance of rt and thl (overall)  [(kg/kg)K]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: rtpthlp
  ! Array of CLUBB's tunable parameters    [units vary]
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params

  !--------------------------- Input/Output Variables ---------------------------
  real(core_rknd), dimension(ngrdcol, nz), intent(inout) :: skrt ! InOut
  real(core_rknd), dimension(ngrdcol, nz), intent(inout) :: skthl ! InOut

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

  call new_pdf_driver( nz, ngrdcol, & ! In
                       wm, rtm, thlm, wp2, rtp2, thlp2, & ! In
                       skw, wprtp, wpthlp, rtpthlp, & ! In
                       clubb_params, & ! In
                       skrt, skthl, stored_pdf_implicit_coefs_terms, & ! InOut
                       mu_w_1, mu_w_2, & ! Out
                       mu_rt_1, mu_rt_2, & ! Out
                       mu_thl_1, mu_thl_2, & ! Out
                       sigma_w_1_sqd, sigma_w_2_sqd, & ! Out
                       sigma_rt_1_sqd, sigma_rt_2_sqd, & ! Out
                       sigma_thl_1_sqd, sigma_thl_2_sqd, & ! Out
                       mixt_frac ) ! Out

end subroutine f2py_new_pdf_driver

