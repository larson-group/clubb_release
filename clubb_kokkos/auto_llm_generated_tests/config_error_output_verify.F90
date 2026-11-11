! Test-only access to Native fixed-width headers; no production bridge is added.
subroutine native_error_headers_c( ngrdcol, chunk_idx, mpi_rank, lat, lon, & ! In
                                    global_header, column_headers ) bind(C) ! Out
  ! Description:
  !   Initialize call-local Native error context and return its fixed character
  !   records for edge cases that normal successful SCM runs do not print.
  use, intrinsic :: iso_c_binding, only: c_int, c_char
  use clubb_precision, only: core_rknd
  use err_info_type_module, only: err_info_type, init_err_info_api, cleanup_err_info_api

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer(c_int), value, intent(in) :: &
    ngrdcol, chunk_idx, mpi_rank ! Column count, chunk index and MPI rank [-]
  real( kind = core_rknd ), intent(in) :: &
    lat(ngrdcol), lon(ngrdcol) ! Raw input latitude/longitude [degrees]

  !--------------------------- Output Variables ---------------------------
  character( kind = c_char ), intent(out) :: &
    global_header(214), column_headers(200,ngrdcol) ! Native fixed-width records [-]

  !--------------------------- Local Variables ---------------------------
  type(err_info_type) :: err_info ! Local context released before returning [-]
  integer :: i, c ! Column and character positions [#]

  !--------------------------- Begin Code ---------------------------
  call init_err_info_api( ngrdcol, chunk_idx, mpi_rank, lat, lon, & ! In
                           err_info )                            ! InOut
  do c = 1, 214
    global_header(c) = err_info%err_header_global(c:c)
  end do
  do i = 1, ngrdcol
    do c = 1, 200
      column_headers(c,i) = err_info%err_header(i)(c:c)
    end do
  end do
  call cleanup_err_info_api( err_info ) ! InOut
end subroutine native_error_headers_c
