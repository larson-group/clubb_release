! new_pdf_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module new_pdf

subroutine f2py_calc_setter_var_params( nz, ngrdcol, xm, xp2, skx, sgn_wpxp, f_x, & ! In
                                        zeta_x, & ! In
                                        mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, & ! Out
                                        mixt_frac, coef_sigma_x_1_sqd, & ! Out
                                        coef_sigma_x_2_sqd ) ! Out

  ! Description:
  !   Call calc_setter_var_params with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use new_pdf, only: calc_setter_var_params

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: nz ! Number of vertical levels [-]
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  ! Mean of x (overall)                             [units vary]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: xm
  ! Variance of x (overall)                     [(units vary)^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: xp2
  ! Skewness of x                                            [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: skx
  ! Sign of the covariance of w and x (overall)              [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: sgn_wpxp
  ! Parameter for the spread of the PDF component means of x [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: f_x
  ! Parameter for the PDF component variances of x           [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: zeta_x

  !--------------------------- Output Variables ---------------------------
  ! Mean of x (1st PDF component)                  [units vary]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mu_x_1
  ! Mean of x (2nd PDF component)                  [units vary]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mu_x_2
  ! Standard deviation of x (1st PDF component)    [units vary]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_x_1
  ! Standard deviation of x (2nd PDF component)    [units vary]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_x_2
  ! Mixture fraction                                        [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mixt_frac
  ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: coef_sigma_x_1_sqd
  ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: coef_sigma_x_2_sqd

  !--------------------------- Begin Code ---------------------------

  call calc_setter_var_params( nz, ngrdcol, & ! In
                               xm, xp2, skx, & ! In
                               sgn_wpxp, f_x, zeta_x, & ! In
                               mu_x_1, mu_x_2, & ! Out
                               sigma_x_1, sigma_x_2, & ! Out
                               mixt_frac, & ! Out
                               coef_sigma_x_1_sqd, coef_sigma_x_2_sqd ) ! Out

end subroutine f2py_calc_setter_var_params

subroutine f2py_calc_responder_params( nz, ngrdcol, xm, xp2, skx, sgn_wpxp, f_x, & ! In
                                       mixt_frac, & ! In
                                       mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, & ! Out
                                       coef_sigma_x_1_sqd, coef_sigma_x_2_sqd ) ! Out

  ! Description:
  !   Call calc_responder_params with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use new_pdf, only: calc_responder_params

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: nz ! Number of vertical levels [-]
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  ! Mean of x (overall)                             [units vary]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: xm
  ! Variance of x (overall)                     [(units vary)^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: xp2
  ! Skewness of x                                            [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: skx
  ! Sign of the covariance of w and x (overall)              [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: sgn_wpxp
  ! Parameter for the spread of the PDF component means of x [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: f_x
  ! Mixture fraction                                         [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: mixt_frac

  !--------------------------- Output Variables ---------------------------
  ! Mean of x (1st PDF component)              [units vary]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mu_x_1
  ! Mean of x (2nd PDF component)              [units vary]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: mu_x_2
  ! Variance of x (1st PDF component)      [(units vary)^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_x_1_sqd
  ! Variance of x (2nd PDF component)      [(units vary)^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: sigma_x_2_sqd
  ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: coef_sigma_x_1_sqd
  ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: coef_sigma_x_2_sqd

  !--------------------------- Begin Code ---------------------------

  call calc_responder_params( nz, ngrdcol, & ! In
                              xm, xp2, skx, & ! In
                              sgn_wpxp, f_x, mixt_frac, & ! In
                              mu_x_1, mu_x_2, & ! Out
                              sigma_x_1_sqd, sigma_x_2_sqd, & ! Out
                              coef_sigma_x_1_sqd, coef_sigma_x_2_sqd ) ! Out

end subroutine f2py_calc_responder_params

subroutine f2py_calc_limits_f_x_responder( nz, ngrdcol, mixt_frac, skx, sgn_wpxp, & ! In
                                           max_skx2_pos_skx_sgn_wpxp, & ! In
                                           max_skx2_neg_skx_sgn_wpxp, & ! In
                                           min_f_x, max_f_x ) ! Out

  ! Description:
  !   Call calc_limits_f_x_responder with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use new_pdf, only: calc_limits_f_x_responder

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: nz ! Number of vertical levels [-]
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  ! Mixture fraction                 [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: mixt_frac
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: skx ! Skewness of x                    [-]
  ! Sign of covariance of w and x    [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: sgn_wpxp
  ! Maximum Skx^2 when Skx*sgn(<w'x'>) >= 0 [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: max_skx2_pos_skx_sgn_wpxp
  ! Maximum Skx^2 when Skx*sgn(<w'x'>) < 0  [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: max_skx2_neg_skx_sgn_wpxp

  !--------------------------- Output Variables ---------------------------
  ! Minimum allowable value of F_x    [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: min_f_x
  ! Maximum allowable value of F_x    [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: max_f_x

  !--------------------------- Begin Code ---------------------------

  call calc_limits_f_x_responder(nz, ngrdcol, mixt_frac, skx, sgn_wpxp, & ! In
    max_skx2_pos_skx_sgn_wpxp, max_skx2_neg_skx_sgn_wpxp, & ! In
    min_f_x, max_f_x) ! Out

end subroutine f2py_calc_limits_f_x_responder

subroutine f2py_calc_coef_wp4_implicit( nz, ngrdcol, mixt_frac, f_w, & ! In
                                        coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, & ! In
                                        coef_wp4_implicit ) ! Out

  ! Description:
  !   Call calc_coef_wp4_implicit with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use new_pdf, only: calc_coef_wp4_implicit

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: nz ! Number of vertical levels [-]
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  ! Mixture fraction                               [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: mixt_frac
  ! Parameter: spread of the PDF comp. means of w  [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: f_w
  ! sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: coef_sigma_w_1_sqd
  ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: coef_sigma_w_2_sqd

  !--------------------------- Output Variables ---------------------------
  ! Coef.: <w'^4> = coef_wp4_implicit * <w'^2>^2    [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: coef_wp4_implicit

  !--------------------------- Begin Code ---------------------------

  coef_wp4_implicit = calc_coef_wp4_implicit( nz, ngrdcol, mixt_frac, f_w, & ! In
                                              coef_sigma_w_1_sqd, & ! In
                                              coef_sigma_w_2_sqd ) ! In

end subroutine f2py_calc_coef_wp4_implicit

subroutine f2py_calc_coef_wpxp2_implicit( nz, ngrdcol, wp2, xp2, wpxp, sgn_wpxp, & ! In
                                          mixt_frac, f_w, f_x, coef_sigma_w_1_sqd, & ! In
                                          coef_sigma_w_2_sqd, coef_sigma_x_1_sqd, & ! In
                                          coef_sigma_x_2_sqd, & ! In
                                          coef_wpxp2_implicit ) ! Out

  ! Description:
  !   Call calc_coef_wpxp2_implicit with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use new_pdf, only: calc_coef_wpxp2_implicit

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: nz ! Number of vertical levels [-]
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  ! Variance of w (overall)                  [m^2/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wp2
  ! Variance of x (overall)           [(units vary)^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: xp2
  ! Covariance of w and x           [m/s (units vary)]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wpxp
  ! Sign of the covariance of w and x              [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: sgn_wpxp
  ! Mixture fraction                               [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: mixt_frac
  ! Parameter: spread of the PDF comp. means of w  [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: f_w
  ! Parameter: spread of the PDF comp. means of x  [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: f_x
  ! sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: coef_sigma_w_1_sqd
  ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: coef_sigma_w_2_sqd
  ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: coef_sigma_x_1_sqd
  ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: coef_sigma_x_2_sqd

  !--------------------------- Output Variables ---------------------------
  ! Coef.: <w'x'^2> = coef_wpxp2_implicit * <x'^2> [m/s]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: coef_wpxp2_implicit

  !--------------------------- Begin Code ---------------------------

  coef_wpxp2_implicit = calc_coef_wpxp2_implicit( nz, ngrdcol, wp2, xp2, wpxp, sgn_wpxp, & ! In
                                                  mixt_frac, f_w, f_x, & ! In
                                                  coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, & ! In
                                                  coef_sigma_x_1_sqd, & ! In
                                                  coef_sigma_x_2_sqd ) ! In

end subroutine f2py_calc_coef_wpxp2_implicit

subroutine f2py_calc_coefs_wp2xp_semiimpl( nz, ngrdcol, wp2, xp2, sgn_wpxp, mixt_frac, & ! In
                                           f_w, f_x, coef_sigma_w_1_sqd, & ! In
                                           coef_sigma_w_2_sqd, coef_sigma_x_1_sqd, & ! In
                                           coef_sigma_x_2_sqd, & ! In
                                           coef_wp2xp_implicit, term_wp2xp_explicit ) ! Out

  ! Description:
  !   Call calc_coefs_wp2xp_semiimpl with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use new_pdf, only: calc_coefs_wp2xp_semiimpl

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: nz ! Number of vertical levels [-]
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  ! Variance of w (overall)                  [m^2/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wp2
  ! Variance of x (overall)           [(units vary)^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: xp2
  ! Sign of the covariance of w and x              [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: sgn_wpxp
  ! Mixture fraction                               [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: mixt_frac
  ! Parameter: spread of the PDF comp. means of w  [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: f_w
  ! Parameter: spread of the PDF comp. means of x  [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: f_x
  ! sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: coef_sigma_w_1_sqd
  ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: coef_sigma_w_2_sqd
  ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: coef_sigma_x_1_sqd
  ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: coef_sigma_x_2_sqd

  !--------------------------- Output Variables ---------------------------
  ! Coefficient that is multiplied by <w'x'>    [m/s]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: coef_wp2xp_implicit
  ! Term that is on the RHS    [m^2/s^2 (units vary)]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: term_wp2xp_explicit

  !--------------------------- Begin Code ---------------------------

  call calc_coefs_wp2xp_semiimpl( nz, ngrdcol, & ! In
                                  wp2, xp2, sgn_wpxp, mixt_frac, f_w, f_x, & ! In
                                  coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, & ! In
                                  coef_sigma_x_1_sqd, coef_sigma_x_2_sqd, & ! In
                                  coef_wp2xp_implicit, term_wp2xp_explicit ) ! Out

end subroutine f2py_calc_coefs_wp2xp_semiimpl

subroutine f2py_calc_coefs_wpxpyp_semiimpl( nz, ngrdcol, wp2, xp2, yp2, wpxp, wpyp, & ! In
                                            sgn_wpxp, sgn_wpyp, mixt_frac, f_w, f_x, & ! In
                                            f_y, coef_sigma_w_1_sqd, & ! In
                                            coef_sigma_w_2_sqd, coef_sigma_x_1_sqd, & ! In
                                            coef_sigma_x_2_sqd, coef_sigma_y_1_sqd, & ! In
                                            coef_sigma_y_2_sqd, & ! In
                                            coef_wpxpyp_implicit, & ! Out
                                            term_wpxpyp_explicit ) ! Out

  ! Description:
  !   Call calc_coefs_wpxpyp_semiimpl with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use new_pdf, only: calc_coefs_wpxpyp_semiimpl

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: nz ! Number of vertical levels [-]
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  ! Variance of w (overall)                  [m^2/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wp2
  ! x'^2 (momentum levels)                   [units vary]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: xp2
  ! Variance of y, y'^2 (momentum levels)         [{y units}^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: yp2
  ! Covariance of w and x (overall)     [m/s(x units)]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wpxp
  ! Covariance of w and y (overall)     [m/s(y units)]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wpyp
  ! Sgn(<w'x'>); 1 when <w'x'> >= 0 or -1 when <w'x'> < 0 [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: sgn_wpxp
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: sgn_wpyp ! In
  ! Mixture fraction                               [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: mixt_frac
  ! Parameter: spread of the PDF comp. means of w  [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: f_w
  ! Parameter for the spread of the PDF component means of x [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: f_x
  ! Parameter: spread of the PDF comp. means of y  [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: f_y
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: coef_sigma_w_1_sqd ! Out
  ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>     [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: coef_sigma_w_2_sqd
  ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: coef_sigma_x_1_sqd
  ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: coef_sigma_x_2_sqd
  ! sigma_y_1^2 = coef_sigma_y_1_sqd * <y'^2>      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: coef_sigma_y_1_sqd
  ! sigma_y_2^2 = coef_sigma_y_2_sqd * <y'^2>      [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: coef_sigma_y_2_sqd

  !--------------------------- Output Variables ---------------------------
  ! Coefficient that is multiplied by <x'y'>   [m/s]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: coef_wpxpyp_implicit
  ! Term that is on the RHS  [m/s(x units)(y units)]
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: term_wpxpyp_explicit

  !--------------------------- Begin Code ---------------------------

  call calc_coefs_wpxpyp_semiimpl( nz, ngrdcol, & ! In
                                   wp2, xp2, yp2, wpxp, wpyp, & ! In
                                   sgn_wpxp, sgn_wpyp, & ! In
                                   mixt_frac, f_w, f_x, f_y, & ! In
                                   coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, & ! In
                                   coef_sigma_x_1_sqd, coef_sigma_x_2_sqd, & ! In
                                   coef_sigma_y_1_sqd, coef_sigma_y_2_sqd, & ! In
                                   coef_wpxpyp_implicit, term_wpxpyp_explicit ) ! Out

end subroutine f2py_calc_coefs_wpxpyp_semiimpl
