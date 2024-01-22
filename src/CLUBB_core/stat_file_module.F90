!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module stat_file_module
 

! Description:
!   Contains two derived types for describing the contents and location of
!   either NetCDF or GrADS files.
!-------------------------------------------------------------------------------
   use clubb_precision, only: & 
       stat_rknd,  & ! Variable
       time_precision, &
       core_rknd
 
   implicit none

   public :: grid_avg_variable, samples_of_variable, stat_file

   ! These are used in a 2D or 3D host model to output multiple columns
   ! Set clubb_i and clubb_j according to the column within the host model;
   ! The indices must not exceed nlon (for i) or nlat (for j).
   !
   !!!!!!!!! Unused Warning !!!!!!!!!!!!
   ! Currently we do not save stats output with lat and lon indexing, and
   ! because of this, we hardcode these as 1. When implementing the multiple 
   ! column stats outputting in the future, we should collapse the lat and 
   ! lon dimension, and use a generic "column" dimension instead.
   integer, parameter, public :: &
    clubb_i = 1, &
    clubb_j = 1

   private ! Default scope

  ! Structures to hold the description of a variable

   type grid_avg_variable
     ! Pointer to the array
     real(kind=stat_rknd), dimension(:,:,:), pointer :: ptr

     character(len = 30) :: name        ! Variable name
     character(len = 100) :: description ! Variable description
     character(len = 25) :: units       ! Variable units

     integer :: indx ! NetCDF module Id for var / GrADS index

     logical :: l_silhs ! If true, we sample this variable once for each SILHS
                        ! sample point per timestep, rather than just once per
                        ! timestep.
   end type grid_avg_variable
!$omp declare mapper (grid_avg_variable::x) map ( &
!$omp  x%ptr &
!$omp , x%name &
!$omp , x%description &
!$omp , x%units &
!$omp , x%indx &
!$omp , x%l_silhs &
!$omp )

   type samples_of_variable
     ! Pointer to the array
     real(kind=stat_rknd), dimension(:,:,:,:), pointer :: ptr

     character(len = 30) :: name        ! Variable name
     character(len = 100) :: description ! Variable description
     character(len = 25) :: units       ! Variable units

     integer :: indx ! NetCDF module Id for var / GrADS index

     logical :: l_silhs ! If true, we sample this variable once for each SILHS
                        ! sample point per timestep, rather than just once per
                        ! timestep.
   end type samples_of_variable
!$omp declare mapper (samples_of_variable::x) map ( &
!$omp  x%ptr &
!$omp , x%name &
!$omp , x%description &
!$omp , x%units &
!$omp , x%indx &
!$omp , x%l_silhs &
!$omp )

  ! Structure to hold the description of a NetCDF output file
  ! This makes the new code as compatible as possible with the
  ! GrADS output code

   type stat_file

     ! File information

     character(len = 200) ::  &
       fname,   & ! File name without suffix
       fdir    ! Path where fname resides

     integer :: iounit  ! This number is used internally by the
                        ! NetCDF module to track the data set, or by
                        ! GrADS to track the actual file unit.
     integer :: &
       nrecord, & ! Number of records written
       ntimes     ! Number of times written

     logical :: &
       l_defined,  &  ! Whether nf90_enddef() has been called
       l_byte_swapped ! Is this a file in the opposite byte ordering?

     ! NetCDF datafile dimensions indices (Samp*Id for SILHS samples)
     integer ::  & 
       SampDimId, LatDimId, LongDimId, AltDimId, TimeDimId, &
       SampVarId, LatVarId, LongVarId, AltVarId, TimeVarId

     ! Grid information

     integer :: ia, iz  ! Vertical extent

     integer :: nlat, nlon ! The number of points in the X and Y

     ! Number of SILHS samples (i.e. subcolumns).  Initialized to zero
     ! to be safe, but will be updated if appropriate
     integer :: nsamp = 0

     real( kind = core_rknd ), dimension(:), allocatable ::  & 
       z ! Height of vertical levels [m]

     ! Time information

     integer :: day, month, year ! Date of starting time

     real( kind = core_rknd ), dimension(:), allocatable :: & 
       lat_vals, & ! Latitude                   [Degrees N]
       lon_vals, & ! Longitude                  [Degrees E]
       samp_idx   ! SILHS subcolumn index

     real( kind = core_rknd ) :: & 
       dtwrite ! Interval between output    [Seconds]

     real( kind = time_precision ) ::  & 
       time    ! Start time                 [Seconds]

     ! Statistical Variables

     integer :: nvar  ! Number of variables for this file

     type (grid_avg_variable), dimension(:), allocatable ::  &
       grid_avg_var ! List and variable description

     type (samples_of_variable), dimension(:), allocatable :: &
       samples_of_var

   end type stat_file
!$omp declare mapper (stat_file::x) map ( &
!$omp  x%fdir &
!$omp , x%iounit &
!$omp , x%ntimes &
!$omp , x%l_byte_swapped &
!$omp , x%sampdimid &
!$omp , x%latdimid &
!$omp , x%longdimid &
!$omp , x%altdimid &
!$omp , x%timedimid &
!$omp , x%sampvarid &
!$omp , x%latvarid &
!$omp , x%longvarid &
!$omp , x%altvarid &
!$omp , x%timevarid &
!$omp , x%ia &
!$omp , x%iz &
!$omp , x%nlat &
!$omp , x%nlon &
!$omp , x%nsamp &
!$omp , x%z &
!$omp , x%day &
!$omp , x%month &
!$omp , x%year &
!$omp , x%samp_idx &
!$omp , x%dtwrite &
!$omp , x%time &
!$omp , x%nvar &
!$omp , x%grid_avg_var &
!$omp , x%samples_of_var &
!$omp )

 end module stat_file_module


