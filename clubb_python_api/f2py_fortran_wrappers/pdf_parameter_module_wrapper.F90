! pdf_parameter_module_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module pdf_parameter_module

subroutine f2py_zero_pdf_params(which)

  use derived_type_storage, only: stored_pdf_params, stored_pdf_params_zm
  use pdf_parameter_module, only: zero_pdf_params_api

  implicit none

  integer, intent(in) :: which  ! 1 = pdf_params, 2 = pdf_params_zm

  if (which == 2) then
    call zero_pdf_params_api(stored_pdf_params_zm)
  else
    call zero_pdf_params_api(stored_pdf_params)
  end if

end subroutine f2py_zero_pdf_params

subroutine f2py_zero_pdf_implicit_coefs_terms()

  use derived_type_storage, only: stored_pdf_implicit_coefs_terms
  use pdf_parameter_module, only: zero_pdf_implicit_coefs_terms_api

  implicit none

  call zero_pdf_implicit_coefs_terms_api(stored_pdf_implicit_coefs_terms)

end subroutine f2py_zero_pdf_implicit_coefs_terms

subroutine f2py_init_pdf_params(nz, ngrdcol) &
  bind(C, name="f2py_init_pdf_params_")

  use derived_type_storage, only: stored_pdf_params
  use pdf_parameter_module, only: init_pdf_params

  implicit none

  integer, intent(in) :: nz, ngrdcol

  call init_pdf_params(nz, ngrdcol, stored_pdf_params)

end subroutine f2py_init_pdf_params

subroutine f2py_init_pdf_params_zm(nz, ngrdcol) &
  bind(C, name="f2py_init_pdf_params_zm_")

  use derived_type_storage, only: stored_pdf_params_zm
  use pdf_parameter_module, only: init_pdf_params

  implicit none

  integer, intent(in) :: nz, ngrdcol

  call init_pdf_params(nz, ngrdcol, stored_pdf_params_zm)

end subroutine f2py_init_pdf_params_zm

subroutine f2py_init_pdf_implicit(nz, ngrdcol, sclr_dim) &
  bind(C, name="f2py_init_pdf_implicit_")

  use derived_type_storage, only: stored_pdf_implicit_coefs_terms
  use pdf_parameter_module, only: init_pdf_implicit_coefs_terms_api

  implicit none

  integer, intent(in) :: nz, ngrdcol, sclr_dim

  call init_pdf_implicit_coefs_terms_api(nz, ngrdcol, sclr_dim, &
                                         stored_pdf_implicit_coefs_terms)

end subroutine f2py_init_pdf_implicit
