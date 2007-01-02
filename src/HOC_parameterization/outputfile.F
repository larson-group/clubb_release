!-----------------------------------------------------------------------
! $Id: outputfile.F,v 1.3 2007-01-02 23:27:08 dschanen Exp $

      module outputfile_class

!     Description:
!     Contains two derived types for either NetCDF or GrADS files.
!-----------------------------------------------------------------------
       
         implicit none

        ! Structure to hold the description of a variable

         type variable
           real, dimension(:), pointer :: ptr ! Pointer to the array

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

           character(len = 200) :: 
     .     fname,  ! File name without suffix
     .     fdir    ! Path where fname resides

           integer :: iounit  ! This number is used internally by the 
                              ! NetCDF module to track the data set, or by 
                              ! GrADS to track the actual file unit.
                                         
           integer :: nrecord  ! Number of records written
           integer :: ntimes   ! Number of times written
           logical :: ldefined ! Whether nf90_enddef() has been called

        ! NetCDF datafile dimensions indices
           integer :: 
     .     LatDimId, LongDimId, AltDimId, TimeDimId,
     .     LatVarId, LongVarId, AltVarId, TimeVarId


        ! Grid information

           integer :: ia, iz  ! Vertical extent

           real, dimension(:), pointer :: 
     .     z ! Height of vertical levels [m]

        ! Time information

           integer :: day, month, year ! Date of starting time

           real :: 
     .     rlat,   ! Latitude                   [Degrees N]
     .     rlon,   ! Longitude                  [Degrees E]
     .     time,   ! Start time                 [Seconds]
     .     dtwrite ! Interval between output    [Seconds]

        ! Statistical Variables

           integer :: nvar  ! Number of variables for this file

           type (variable), dimension(:), pointer :: 
     .     var ! List and variable description

         end type outputfile

       end module outputfile_class
