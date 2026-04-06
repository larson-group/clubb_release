! file_functions_wrapper.F90 — wrappers organized by source module

subroutine f2py_file_read_1d(file_unit, path_and_filename, num_datapts, &
    entries_per_line, variable)

  use clubb_precision, only: core_rknd
  use file_functions, only: file_read_1d

  implicit none

  integer, intent(in) :: file_unit, num_datapts, entries_per_line
  character(len=*), intent(in) :: path_and_filename
  real(core_rknd), dimension(num_datapts), intent(out) :: variable

  call file_read_1d(file_unit, path_and_filename, num_datapts, entries_per_line, variable)

end subroutine f2py_file_read_1d
