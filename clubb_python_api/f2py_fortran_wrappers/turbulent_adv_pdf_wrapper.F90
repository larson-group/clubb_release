subroutine f2py_xpyp_term_ta_pdf_lhs(nzm, nzt, ngrdcol, coef_wpxpyp_implicit, rho_ds_zt, &
    rho_ds_zm, invrs_rho_ds_zm, l_upwind_xpyp_turbulent_adv, sgn_turbulent_vel, &
    coef_wpxpyp_implicit_zm, lhs_ta)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use turbulent_adv_pdf, only: xpyp_term_ta_pdf_lhs

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  logical, intent(in) :: l_upwind_xpyp_turbulent_adv
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: coef_wpxpyp_implicit, rho_ds_zt
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: &
    rho_ds_zm, invrs_rho_ds_zm, sgn_turbulent_vel, coef_wpxpyp_implicit_zm
  real(core_rknd), dimension(3, ngrdcol, nzm), intent(out) :: lhs_ta

  call xpyp_term_ta_pdf_lhs(nzm, nzt, ngrdcol, stored_grid, coef_wpxpyp_implicit, rho_ds_zt, &
    rho_ds_zm, invrs_rho_ds_zm, l_upwind_xpyp_turbulent_adv, sgn_turbulent_vel, &
    coef_wpxpyp_implicit_zm, lhs_ta)

end subroutine f2py_xpyp_term_ta_pdf_lhs

subroutine f2py_xpyp_term_ta_pdf_lhs_godunov(nzm, nzt, ngrdcol, coef_wpxpyp_implicit, &
    invrs_rho_ds_zm, rho_ds_zm, lhs_ta)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use turbulent_adv_pdf, only: xpyp_term_ta_pdf_lhs_godunov

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: coef_wpxpyp_implicit
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: invrs_rho_ds_zm, rho_ds_zm
  real(core_rknd), dimension(3, ngrdcol, nzm), intent(out) :: lhs_ta

  call xpyp_term_ta_pdf_lhs_godunov(nzm, nzt, ngrdcol, stored_grid, coef_wpxpyp_implicit, &
    invrs_rho_ds_zm, rho_ds_zm, lhs_ta)

end subroutine f2py_xpyp_term_ta_pdf_lhs_godunov

subroutine f2py_xpyp_term_ta_pdf_rhs(nzm, nzt, ngrdcol, term_wpxpyp_explicit, rho_ds_zt, &
    rho_ds_zm, invrs_rho_ds_zm, l_upwind_xpyp_turbulent_adv, sgn_turbulent_vel, &
    term_wpxpyp_explicit_zm, rhs_ta)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use turbulent_adv_pdf, only: xpyp_term_ta_pdf_rhs

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  logical, intent(in) :: l_upwind_xpyp_turbulent_adv
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: term_wpxpyp_explicit, rho_ds_zt
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: &
    rho_ds_zm, invrs_rho_ds_zm, sgn_turbulent_vel, term_wpxpyp_explicit_zm
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: rhs_ta

  call xpyp_term_ta_pdf_rhs(nzm, nzt, ngrdcol, stored_grid, term_wpxpyp_explicit, rho_ds_zt, &
    rho_ds_zm, invrs_rho_ds_zm, l_upwind_xpyp_turbulent_adv, sgn_turbulent_vel, &
    term_wpxpyp_explicit_zm, rhs_ta)

end subroutine f2py_xpyp_term_ta_pdf_rhs

subroutine f2py_xpyp_term_ta_pdf_rhs_godunov(nzm, nzt, ngrdcol, term_wpxpyp_explicit_zm, &
    invrs_rho_ds_zm, sgn_turbulent_vel, rho_ds_zm, rhs_ta)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use turbulent_adv_pdf, only: xpyp_term_ta_pdf_rhs_godunov

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: term_wpxpyp_explicit_zm, invrs_rho_ds_zm, rho_ds_zm
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: sgn_turbulent_vel
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: rhs_ta

  call xpyp_term_ta_pdf_rhs_godunov(nzm, nzt, ngrdcol, stored_grid, term_wpxpyp_explicit_zm, &
    invrs_rho_ds_zm, sgn_turbulent_vel, rho_ds_zm, rhs_ta)

end subroutine f2py_xpyp_term_ta_pdf_rhs_godunov
