subroutine f2py_assert_corr_symmetric(pdf_dim, corr_array_n)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_err_info
  use corr_varnce_module, only: assert_corr_symmetric

  implicit none

  integer, intent(in) :: pdf_dim
  real(core_rknd), dimension(pdf_dim, pdf_dim), intent(in) :: corr_array_n

  call assert_corr_symmetric(corr_array_n, pdf_dim, stored_err_info)

end subroutine f2py_assert_corr_symmetric
