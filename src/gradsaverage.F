      module grads_common

      contains
!----------------------------------------------------------------------
!     Function GRADS_AVERAGE( )
!     Average a GRADS file variable over the interval t1 to t2
!----------------------------------------------------------------------
      function grads_average( filename, nz, t1, t2, variable_name, 
     .                        error )

      use grid_class
      use inputfile_class

      implicit none

!---- Arguments -------------------------------------------------------
      character(len=*), intent(in) :: filename
      integer, intent(in)          :: nz             ! number of z-levels
      integer, intent(in)          :: t1             ! beginning timestep to look at
      integer, intent(in)          :: t2             ! ending timestep to look at
      character(len=*), intent(in) :: variable_name  ! which variable to read
      logical, intent(out)         :: error          ! error status from reading 

!---- Internal Variables ----------------------------------------------

      integer i                             ! les_array loop variable
      integer t                             ! timestep loop variable

      integer num_timesteps                 ! steps between t1 and t2 

      type (inputgrads) faverage            ! Data file

      real, dimension(nz)  :: grads_temp

!---- Return ----------------------------------------------------------
      real, dimension(nz)  :: grads_average

!-----------------------------------------------------------------------
!  Open grads file
!-----------------------------------------------------------------------
      call open_grads_read(10, filename, faverage)

!-----------------------------------------------------------------------
! initialize variables
!-----------------------------------------------------------------------
      num_timesteps = ( t2 - t1 ) + 1
      grads_temp    = 0.
      grads_average = 0.

!-----------------------------------------------------------------------
! read in floating point values
!-----------------------------------------------------------------------
        do t = t1, t2
          call get_var( faverage, variable_name, t, 
     .                  grads_temp(1:nz), error ) 
          if ( error ) then
             print *, "grads_average: get_var failed for " 
     .         // trim( variable_name ) // " in " // trim( filename )
     .         // " at time=", t
             return 
          endif
          grads_average(1:nz) = grads_average(1:nz) + grads_temp(1:nz) 
        enddo

!-----------------------------------------------------------------------
! close data file
!-----------------------------------------------------------------------
      call close_grads_read( faverage )

!-----------------------------------------------------------------------
! Take average over num_timesteps
!-----------------------------------------------------------------------
      grads_average(1:nz) = grads_average(1:nz) / real(num_timesteps)

      return
      end function grads_average

!-------------------------------------------------------------------------
!   Function GRADS_AVERAGE_INTERVAL( )
!   Take several averages over an interval.
!-------------------------------------------------------------------------
      function grads_average_interval( filename, nz, t, variable_name,
     .                                 error )

      implicit none

!---- Arguments -------------------------------------------------------
      character(len=*), intent(in)         :: filename
      integer, intent(in)                  :: nz             ! number of z-levels
      integer, dimension(:), intent(in)    :: t(:)           ! timesteps to look at
      character(len=*), intent(in)         :: variable_name  ! which variable to read
      logical, intent(out)                 :: error

!---- Return Type  ----------------------------------------------------
      real, dimension(nz)                  :: grads_average_interval

!---- Internal Variables ----------------------------------------------
      integer                              :: i       ! loop variable 
      integer                              :: tdim    ! dimension of t variable
      real, dimension(nz)                  :: grads_temp
      integer                              :: divisor

      do i=1, size( t )
        if ( t( i ) == 0 ) exit
        tdim = i
      enddo

      grads_average_interval = grads_average( filename, nz, t(1), t(2), 
     .                                        variable_name, error ) *
     .                                      ( t(2) - t(1) )
      divisor = t(2) - t(1)

      if ( error ) return

      do i=3, tdim, 2 
        grads_temp = grads_average( filename, nz, t(i), t(i+1), 
     .                              variable_name, error )
        grads_average_interval = grads_average_interval + grads_temp *
     .                           ( t(i+1) - t(i) )
        divisor =  divisor +  ( t(i+1) - t(i) )
      enddo

      grads_average_interval(1:nz) = grads_average_interval(1:nz)
     .                                 / real(divisor)

      return
      end function grads_average_interval

!-------------------------------------------------------------------------
!   Function GRADS_ZLVL( )
!   returns a scalar of the number of z-levels in a file
!-------------------------------------------------------------------------
      integer function grads_zlvl( filename )

      use grid_class
      use inputfile_class

      implicit none

      character(len=*), intent(in) :: filename

      type (inputgrads) fz            ! Data file

      call open_grads_read( 10, filename, fz )

      grads_zlvl = fz%iz

      call close_grads_read( fz )

      return
      end function grads_zlvl

!-------------------------------------------------------------------------
      end module grads_common
