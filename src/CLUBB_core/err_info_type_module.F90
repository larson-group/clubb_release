!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module err_info_type_module

  ! Description:
  !   Contains derived data type 'err_code'.
  !   Used for passing error info along the call tree.
  !-----------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd

  implicit none

  public :: err_info_type, init_err_info_api, init_default_err_info_api, &
            set_err_info_values_api, cleanup_err_info_api

  private :: update_headers

  ! Derived data types to pass error info along call tree
  ! We don't want the model info variables to be accessible from the outside
  ! since the error headers need to be reset each time one of these four variables changes:
  ! host_lat, host_lon, chunk_idx, and mpi_rank
  ! ->
  ! USE THE SUBROUTINE `SET_ERR_INFO_VALUES_API` TO SET MODEL INFO VALUES FOR ERR_INFO
  ! ERR_CODE CAN BE SET IN THE USUAL WAY BY ASSIGNMENT
  type err_info_type

    ! How many dimensions, if any, do we need?
    ! Could use thread private
    integer, allocatable, public, dimension(:) :: &
      err_code            ! Error code indicating error status of CLUBB

    ! NOTE: CAM uses selected_real_kind(12) for these variables
    ! SAM uses core_rknd in clubb_sgs
    ! WRF uses default reals
    ! So we might lose precision here.
    ! Although it shouldn't matter that much since we're not doing calculations with these values
    ! and even single precision should be good enough for lat/lon values.
    real(kind=core_rknd), allocatable, private, dimension(:) :: &
      host_lat,              & ! Latitude, -90 to 90 degrees
      host_lon                 ! Longitude, -180 to 180 degrees

    integer, private :: &
      chunk_idx,        & ! Index of the chunk of columns that CLUBB is running on
      mpi_rank            ! Process rank if the run uses MPI

    character(len=200), allocatable, public, dimension(:) :: err_header

    character(len=214), public :: err_header_global
  end type err_info_type

  contains

  !-----------------------------------------------------------------------
  subroutine init_default_err_info_api(ngrdcol, err_info)
    ! Description: Allocate and fill an err_info_type variables with default values
    !              Call this if you do not have real values for err_info
    !
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
      unused_var      ! Constant

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      ngrdcol             ! Number of grid columns for this CLUBB run

    !--------------------------- Input/Output Variables ---------------------------
    type (err_info_type), intent(inout) :: &
      err_info

    !--------------------------- Local Variables ---------------------------

    ! Dummy default values for parallelization info
    integer :: &
      mpi_rank_dummy = 1,   &
      chunk_idx_dummy = 1

    ! ---------------------- Begin Code ----------------------

    ! Call to allocate and fill variables with values:
    ! err_code is set to clubb_no_error
    ! MPI rank and chunk index set to 1
    ! All latitude and longitude entries set to -999.
    call init_err_info_api(ngrdcol, mpi_rank_dummy, chunk_idx_dummy, &
                           reshape((/unused_var/), (/ngrdcol/), (/unused_var/)), &
                           reshape((/unused_var/), (/ngrdcol/), (/unused_var/)), &
                           err_info)

  end subroutine init_default_err_info_api

  !=============================================================================
  subroutine init_err_info_api(ngrdcol, chunk_idx_in, mpi_rank_in, lat_in, lon_in, err_info)
    ! Description: Allocate and fill an err_info_type variables with input values
    !              Call this if you do not have real values for err_info
    !
    !-----------------------------------------------------------------------

    use error_code, only : &
      clubb_no_error

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      ngrdcol,          & ! Number of grid columns for this CLUBB run
      chunk_idx_in,     & ! Index of the chunk of columns that CLUBB is running on
      mpi_rank_in         ! MPI rank of the process

    real(kind=core_rknd), dimension(ngrdcol), intent(in) :: &
      lat_in,           & ! Latitude, -90 to 90 degrees
      lon_in              ! Longitude, -180 to 180 degrees

    !--------------------------- Input/Output Variables ---------------------------
    type (err_info_type), intent(inout) :: &
      err_info

    !--------------------------- Local Variables ---------------------------
    integer :: i
    ! ---------------------- Begin Code ----------------------

    ! Allocate arrays
    allocate( err_info%host_lat(1:ngrdcol))
    allocate( err_info%host_lon(1:ngrdcol))
    allocate( err_info%err_code(1:ngrdcol))
    allocate( err_info%err_header(1:ngrdcol))

    ! Populate err_code array
    do i = 1, ngrdcol
      err_info%err_code(i) = clubb_no_error
    end do

    ! Populate rest of info variables and headers
    call set_err_info_values_api(ngrdcol, err_info, chunk_idx_in, mpi_rank_in, lat_in, lon_in)

  end subroutine init_err_info_api

  !=============================================================================

  subroutine set_err_info_values_api(ngrdcol, err_info, chunk_idx_in, mpi_rank_in, lat_in, lon_in)
    ! Description: Set the info variables and rewrite the headers
    !
    ! USE THIS SUBROUTINE TO SET INFO VALUES FOR ERR_INFO
    ! ERR_CODE CAN BE SET IN THE USUAL WAY BY ASSIGNMENT
    !-----------------------------------------------------------------------
    use constants_clubb, only: &
        fstderr  ! Constant(s)

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      ngrdcol             ! Number of grid columns for this CLUBB run

    integer, intent(in), optional :: &
      chunk_idx_in,     & ! Index of the chunk of columns that CLUBB is running on
      mpi_rank_in         ! MPI rank of the process

    real(kind=core_rknd), dimension(ngrdcol), intent(in), optional :: &
      lat_in,           & ! Latitude, -90 to 90 degrees
      lon_in              ! Longitude, -180 to 180 degrees

    !--------------------------- Input/Output Variables ---------------------------
    type (err_info_type), intent(inout) :: &
      err_info

    !--------------------------- Local Variables ---------------------------
    integer :: i

    real(kind=core_rknd) :: &
      tmp

    ! ---------------------- Begin Code ----------------------

    ! Check if err_info variables were allocated before accessing them here
    if ( .not. allocated(err_info%err_header) .or. &
         .not. allocated(err_info%host_lat) .or. &
         .not. allocated(err_info%err_header) ) then
      write(fstderr,*) "Warning! Error assigning values to allocatable arrays in CLUBB's ", &
                       "err_info struct. Arrays were not allocated before use. They will be ", &
                       "allocated now, but it is advisable to call an init_api ", &
                       "subroutines first!"
      ! Assume that all allocatables are not allocated since they only can be allocated together
      ! And since we don't know if all info is passed into this subroutine, we can only call
      ! the default version of init
      call init_default_err_info_api(ngrdcol, err_info)
    end if
    ! Populate parallelization info, if present

    !! MPI
    if ( present(mpi_rank_in) ) then
      err_info%mpi_rank = mpi_rank_in
    end if

    !! Chunk info
    if ( present(chunk_idx_in) ) then
      err_info%chunk_idx = chunk_idx_in
    end if

    ! Populate latitude info from input, if present
    if ( present(lat_in) ) then

      do i = 1, ngrdcol
        ! Adjust latitude value if necessary
        if ( lat_in(i) > 90._core_rknd ) then
            ! If latitude is greater than 90, assume that the range for the input is [0,180]
            ! -> Subtract 90 to get range [-90,90]
            tmp = lat_in(i) - 90._core_rknd
        else
            tmp = lat_in(i)
        end if
        err_info%host_lat(i) = tmp
      end do

    end if

    ! Populate longitude info from input, if present
    if ( present(lon_in) ) then
      do i = 1, ngrdcol
        ! Adjust latitude value if necessary
        if ( lon_in(i) > 180._core_rknd ) then
            ! If longitude is greater than 180, assume that the range for the input is [0,360]
            ! -> Subtract 180 to get range [-180,180]
            tmp = lon_in(i) - 180._core_rknd
        else
            tmp = lon_in(i)
        end if
        err_info%host_lon(i) = tmp
      end do
    end if

    ! Write updated info to headers
    call update_headers(ngrdcol, err_info)

  end subroutine set_err_info_values_api
  !=============================================================================

  subroutine update_headers(ngrdcol, err_info)
    ! Description: Rewrite the headers from the current info values
    !
    !-----------------------------------------------------------------------

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      ngrdcol             ! Number of grid columns for this CLUBB run

    !--------------------------- Input/Output Variables ---------------------------
    type (err_info_type), intent(inout) :: &
      err_info

    !--------------------------- Local Variables ---------------------------
    integer :: i
    ! ---------------------- Begin Code ----------------------

    ! Populate error header character string that we want to append
    ! to the error messages throughout CLUBB

    ! Global (column-independent) header
    write(err_info%err_header_global,'(A38, I3, A15, I3, A, A28, A33, F6.2, A1, F6.2, A, A39, A40)') &
                                     "Fatal error in CLUBB - MPI Process ", err_info%mpi_rank, &
                                     " / Chunk index ", err_info%chunk_idx, &
                                     NEW_LINE('a'), "Column index not available. ", &
                                     "Latitude/Longitude for column 1: ", err_info%host_lat(1), &
                                     "/", err_info%host_lon(1), NEW_LINE('a'), &
                                     "Latitude range from -90 to 90 degrees, ", &
                                     "Longitude range from -180 to 180 degrees"

    ! Populate column-specific headers
    do i = 1, ngrdcol
      write(err_info%err_header(i),'(A38, I3, A15, I3, A, A22, I2, A11, F6.2, A13, F6.2, A, A39, A40)') &
                                   "Fatal error in CLUBB in - MPI Process ", err_info%mpi_rank, &
                                   " / Chunk index ", err_info%chunk_idx, &
                                   NEW_LINE('a'), "Grid column index i = ", i, &
                                   ", Latitude ", err_info%host_lat(i), &
                                   " / Longitude ", err_info%host_lat(i), NEW_LINE('a'), &
                                   "Latitude range from -90 to 90 degrees, ", &
                                   "Longitude range from -180 to 180 degrees"
    end do

  end subroutine update_headers

  !=============================================================================
  subroutine cleanup_err_info_api(err_info)
    ! Description:
    !   De-allocates the memory for err_info
    !
    ! References:
    !   None
    !------------------------------------------------------------------------------
    use constants_clubb, only: &
        fstderr  ! Constant(s)

    implicit none

    type(err_info_type), intent(inout) :: &
      err_info

    ! Local Variable(s)
    integer :: ierr

    ! ----- Begin Code -----

    ! Check if err_info variables were allocated before accessing them here
    if ( allocated(err_info%err_header) ) then
      ! Deallocate memory for column headers
      deallocate( err_info%err_header, stat=ierr )
      if ( ierr /= 0 ) then
        write(fstderr,*) "Deallocation of failed err_info%err_header failed."
      end if
    end if

    if ( allocated(err_info%err_code) ) then
      ! Deallocate memory for column headers
      deallocate( err_info%err_code, stat=ierr )
      if ( ierr /= 0 ) then
        write(fstderr,*) "Deallocation of failed err_info%err_code failed."
      end if
    end if

    if ( allocated(err_info%host_lat) ) then
      ! Deallocate memory for column headers
      deallocate( err_info%host_lat, stat=ierr )
      if ( ierr /= 0 ) then
        write(fstderr,*) "Deallocation of failed err_info%host_lat failed."
      end if
    end if

    if ( allocated(err_info%host_lon) ) then
      ! Deallocate memory for column headers
      deallocate( err_info%host_lon, stat=ierr )
      if ( ierr /= 0 ) then
        write(fstderr,*) "Deallocation of failed err_info%host_lon failed."
      end if
    end if

    return

  end subroutine cleanup_err_info_api

end module err_info_type_module
