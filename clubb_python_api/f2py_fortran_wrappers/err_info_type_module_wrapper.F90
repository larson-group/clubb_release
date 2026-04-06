! err_info_type_module_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module err_info_type_module

subroutine f2py_cleanup_err_info()

  use derived_type_storage, only: &
    stored_err_info, &
    stored_err_info_chunk_idx, stored_err_info_mpi_rank, &
    stored_err_info_lat, stored_err_info_lon
  use err_info_type_module, only: cleanup_err_info_api

  implicit none

  call cleanup_err_info_api(stored_err_info)
  if (allocated(stored_err_info_lat)) deallocate(stored_err_info_lat)
  if (allocated(stored_err_info_lon)) deallocate(stored_err_info_lon)
  stored_err_info_chunk_idx = 1
  stored_err_info_mpi_rank = 0

end subroutine f2py_cleanup_err_info

subroutine f2py_set_err_info_values(ngrdcol, chunk_idx, mpi_rank, lat, lon)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: &
    stored_err_info, &
    stored_err_info_chunk_idx, stored_err_info_mpi_rank, &
    stored_err_info_lat, stored_err_info_lon
  use err_info_type_module, only: set_err_info_values_api

  implicit none

  integer, intent(in) :: ngrdcol, chunk_idx, mpi_rank
  real(core_rknd), dimension(ngrdcol), intent(in) :: lat, lon

  call set_err_info_values_api(ngrdcol, stored_err_info, &
    chunk_idx_in=chunk_idx, mpi_rank_in=mpi_rank, &
    lat_in=lat, lon_in=lon)

  if (.not. allocated(stored_err_info_lat) .or. .not. allocated(stored_err_info_lon)) then
    if (allocated(stored_err_info_lat)) deallocate(stored_err_info_lat)
    if (allocated(stored_err_info_lon)) deallocate(stored_err_info_lon)
    allocate(stored_err_info_lat(ngrdcol))
    allocate(stored_err_info_lon(ngrdcol))
  else if (size(stored_err_info_lat) /= ngrdcol .or. size(stored_err_info_lon) /= ngrdcol) then
    deallocate(stored_err_info_lat)
    deallocate(stored_err_info_lon)
    allocate(stored_err_info_lat(ngrdcol))
    allocate(stored_err_info_lon(ngrdcol))
  end if

  stored_err_info_chunk_idx = chunk_idx
  stored_err_info_mpi_rank = mpi_rank
  stored_err_info_lat = lat
  stored_err_info_lon = lon

end subroutine f2py_set_err_info_values

subroutine f2py_get_err_info_values(ngrdcol, chunk_idx, mpi_rank, lat, lon)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: &
    stored_err_info_chunk_idx, stored_err_info_mpi_rank, &
    stored_err_info_lat, stored_err_info_lon

  implicit none

  integer, intent(in) :: ngrdcol
  integer, intent(out) :: chunk_idx, mpi_rank
  real(core_rknd), dimension(ngrdcol), intent(out) :: lat, lon
  integer :: ncopy

  chunk_idx = stored_err_info_chunk_idx
  mpi_rank = stored_err_info_mpi_rank
  lat = 0.0_core_rknd
  lon = 0.0_core_rknd

  if (allocated(stored_err_info_lat) .and. allocated(stored_err_info_lon)) then
    ncopy = min(ngrdcol, size(stored_err_info_lat), size(stored_err_info_lon))
    if (ncopy > 0) then
      lat(1:ncopy) = stored_err_info_lat(1:ncopy)
      lon(1:ncopy) = stored_err_info_lon(1:ncopy)
    end if
  end if

end subroutine f2py_get_err_info_values

subroutine f2py_init_err_info(ngrdcol) &
  bind(C, name="f2py_init_err_info_")

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: &
    stored_err_info, &
    stored_err_info_chunk_idx, stored_err_info_mpi_rank, &
    stored_err_info_lat, stored_err_info_lon
  use err_info_type_module, only: init_default_err_info_api, &
                                  cleanup_err_info_api

  implicit none

  integer, intent(in) :: ngrdcol

  ! Deallocate if previously allocated (safe to call multiple times)
  if (allocated(stored_err_info%err_code)) then
    call cleanup_err_info_api(stored_err_info)
  end if

  call init_default_err_info_api(ngrdcol, stored_err_info)
  if (allocated(stored_err_info_lat)) deallocate(stored_err_info_lat)
  if (allocated(stored_err_info_lon)) deallocate(stored_err_info_lon)
  allocate(stored_err_info_lat(ngrdcol))
  allocate(stored_err_info_lon(ngrdcol))
  stored_err_info_lat = 0.0_core_rknd
  stored_err_info_lon = 0.0_core_rknd
  stored_err_info_chunk_idx = 1
  stored_err_info_mpi_rank = 0

end subroutine f2py_init_err_info
