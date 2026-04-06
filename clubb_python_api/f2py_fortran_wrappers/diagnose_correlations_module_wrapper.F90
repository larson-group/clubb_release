subroutine f2py_diagnose_correlations(pdf_dim, iipdf_w, corr_array_pre, l_calc_w_corr, corr_array)

  use clubb_precision, only: core_rknd
  use diagnose_correlations_module, only: diagnose_correlations

  implicit none

  integer, intent(in) :: pdf_dim, iipdf_w
  logical, intent(in) :: l_calc_w_corr
  real(core_rknd), dimension(pdf_dim, pdf_dim), intent(in) :: corr_array_pre
  real(core_rknd), dimension(pdf_dim, pdf_dim), intent(out) :: corr_array

  call diagnose_correlations(pdf_dim, iipdf_w, corr_array_pre, l_calc_w_corr, corr_array)

end subroutine f2py_diagnose_correlations

subroutine f2py_calc_cholesky_corr_mtx_approx(n_variables, iipdf_w, corr_matrix, &
    corr_cholesky_mtx, corr_mtx_approx)

  use clubb_precision, only: core_rknd
  use diagnose_correlations_module, only: calc_cholesky_corr_mtx_approx

  implicit none

  integer, intent(in) :: n_variables, iipdf_w
  real(core_rknd), dimension(n_variables, n_variables), intent(in) :: corr_matrix
  real(core_rknd), dimension(n_variables, n_variables), intent(out) :: &
    corr_cholesky_mtx, corr_mtx_approx

  call calc_cholesky_corr_mtx_approx(n_variables, iipdf_w, corr_matrix, corr_cholesky_mtx, corr_mtx_approx)

end subroutine f2py_calc_cholesky_corr_mtx_approx
