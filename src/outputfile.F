!-----------------------------------------------------------------------
! $Id: outputfile.F,v 1.2 2006-03-06 17:03:49 dschanen Exp $
!     MODULE outputfile_class

!     contains a structure for either NetCDF or GrADS
!-----------------------------------------------------------------------

      module outputfile_class
       
         implicit none

!        Structure to hold description of a variable

         type variable
           real, pointer :: ptr(:)          ! pointer to the array
           character (len = 15) name        ! variable name
           character (len = 50) description ! variable description
           character (len = 20) units       ! variable units
           integer Id                       ! NetCDF module Id for var
         end type variable

!        Structure to hold the description of a NetCDF output file
!        This makes the new code as compatible as possible with the old

         type outputfile

!          File information

           character (len = 200) fname   ! file name without suffix
           character (len = 200) fdir    ! path where fname resides
           integer iounit  ! this number is used internally by the 
                           ! netcdf module to track the data set, or by 
                           ! grads to track the actual file unit
                                         
           integer nrecord               ! number of records written
           integer ntimes                ! number of times written
           logical ldefined              ! whether nf90_enddef() has been called

!          NetCDF datafile dimensions indices
           integer LatDimId, LongDimId, AltDimId, TimeDimId
           integer LatVarId, LongVarId, AltVarId, TimeVarId


!          Grid information

           integer ia, iz                ! extent
           real, pointer :: z(:)         ! height of vertical levels

!          Time information

           integer day, month, year      ! date of starting time
           real time                     ! start time in seconds
           real dtwrite                  ! interval between outputs in seconds

!           Variables

           integer nvar                  ! number of variables for this file
           type (variable), pointer :: var(:) ! list and variable description

         end type outputfile

       end module outputfile_class
