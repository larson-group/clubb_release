!-----------------------------------------------------------------------
! $Id: outputfile.F90,v 1.2 2008-07-28 19:34:43 faschinj Exp $
module outputfile_class
#ifdef STATS

!     Description:
!     Contains two derived types for either NetCDF or GrADS files.
!-----------------------------------------------------------------------
   use stats_precision, only: & 
       stat_rknd,  & ! Variable
       time_precision
 
   implicit none

   public :: variable, outputfile

   private ! Default scope
   
  ! Structure to hold the description of a variable

   type variable
     ! Pointer to the array
     real(kind=stat_rknd), dimension(:), pointer :: ptr 

     character(len = 15) :: name        ! Variable name
     character(len = 50) :: description ! Variable description
     character(len = 20) :: units       ! Variable units

     integer :: Id                      ! NetCDF module Id for var
   end type variable

  ! Structure to hold the description of a NetCDF output file
  ! This makes the new code as compatible as possible with the
  ! GrADS output code

   type outputfile

  ! File information

     character(len = 200) ::  & 
     fname,   & ! File name without suffix
     fdir    ! Path where fname resides

     integer :: iounit  ! This number is used internally by the 
                        ! NetCDF module to track the data set, or by 
                        ! GrADS to track the actual file unit.
                                   
     integer :: nrecord  ! Number of records written
     integer :: ntimes   ! Number of times written
     logical :: ldefined ! Whether nf90_enddef() has been called

  ! NetCDF datafile dimensions indices
     integer ::  & 
     LatDimId, LongDimId, AltDimId, TimeDimId, & 
     LatVarId, LongVarId, AltVarId, TimeVarId


  ! Grid information

     integer :: ia, iz  ! Vertical extent

     real, dimension(:), pointer ::  & 
     z ! Height of vertical levels [m]

  ! Time information

     integer :: day, month, year ! Date of starting time

     real ::  & 
     rlat,    & ! Latitude                   [Degrees N]
     rlon    ! Longitude                  [Degrees E]

     real(kind=time_precision) :: & 
     dtwrite ! Interval between output    [Seconds]

     real(kind=time_precision) ::  & 
     time    ! Start time                 [Seconds]

  ! Statistical Variables

     integer :: nvar  ! Number of variables for this file

     type (variable), dimension(:), pointer ::  & 
     var ! List and variable description

   end type outputfile

#endif
 end module outputfile_class
