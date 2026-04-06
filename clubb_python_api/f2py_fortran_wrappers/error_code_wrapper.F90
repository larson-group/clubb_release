! error_code_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module error_code

subroutine f2py_set_clubb_debug_level(level)

  use error_code, only: set_clubb_debug_level_api

  implicit none

  integer, intent(in) :: level

  call set_clubb_debug_level_api(level)

end subroutine f2py_set_clubb_debug_level

subroutine f2py_reset_err_code() &
  bind(C, name="f2py_reset_err_code_")

  use error_code, only: clubb_no_error
  use derived_type_storage, only: stored_err_info

  implicit none

  stored_err_info%err_code = clubb_no_error

end subroutine f2py_reset_err_code

subroutine f2py_initialize_error_headers()

  use error_code, only: initialize_error_headers

  implicit none

  call initialize_error_headers

end subroutine f2py_initialize_error_headers
