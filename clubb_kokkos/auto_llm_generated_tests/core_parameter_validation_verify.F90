! Test-only entry to Native parameter validation; no production bridge is added.
subroutine native_check_parameters_edge_c( ngrdcol, nparams_in, parameters, & ! In
                                            codes ) bind(C)                  ! InOut
  ! Description:
  !   Validate caller-provided column/parameter data with call-local Native
  !   error state. Return column codes without resetting their initial values.
  use, intrinsic :: iso_c_binding, only: c_int
  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use parameters_tunable, only: check_parameters_api
  use err_info_type_module, only: err_info_type
  implicit none

  !--------------------------- Input Variables ---------------------------
  integer(c_int), value, intent(in) :: &
    ngrdcol, nparams_in ! Column and parameter counts [#]
  real( kind = core_rknd ), intent(in) :: &
    parameters(ngrdcol,nparams_in) ! Caller-owned tunable matrix [-]

  !--------------------------- Input/Output Variables ---------------------------
  integer(c_int), intent(inout) :: &
    codes(ngrdcol) ! Initial and final Native status by column [-]

  !--------------------------- Local Variables ---------------------------
  type(err_info_type) :: err_info ! Call-local diagnostic headers and codes [-]
  integer :: i ! Column index [#]

  !--------------------------- Begin Code ---------------------------
  if ( nparams_in /= nparams ) error stop "parameter validator test shape mismatch"
  allocate( err_info%err_code(ngrdcol), err_info%err_header(ngrdcol) )
  err_info%err_code = codes
  do i = 1, ngrdcol
    write( err_info%err_header(i), '(A,I0)' ) "parameter validation column ", i
  end do
  call check_parameters_api( ngrdcol, parameters, & ! In
                              err_info )             ! InOut
  codes = err_info%err_code
end subroutine native_check_parameters_edge_c
