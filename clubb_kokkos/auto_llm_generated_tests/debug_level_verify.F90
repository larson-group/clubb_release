! Test-only access to Native threshold semantics and retained stats metadata.
subroutine native_debug_set_c( level ) bind(C)
  ! Description: Set only the Native test/backend copy of the threshold.
  use, intrinsic :: iso_c_binding, only: c_int
  use error_code, only: set_clubb_debug_level_api
  implicit none
  integer(c_int), value, intent(in) :: level ! Diagnostic threshold [-]
  call set_clubb_debug_level_api( level ) ! In
end subroutine native_debug_set_c

integer(c_int) function native_debug_check_c( level ) bind(C) result(enabled)
  ! Description: Return the actual Native comparison as a C integer.
  use, intrinsic :: iso_c_binding, only: c_int
  use error_code, only: clubb_at_least_debug_level_api
  implicit none
  integer(c_int), value, intent(in) :: level ! Requested threshold [-]
  enabled = 0
  if ( clubb_at_least_debug_level_api(level) ) enabled = 1
end function native_debug_check_c
