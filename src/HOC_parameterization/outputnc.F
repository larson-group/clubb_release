! $Id: outputnc.F,v 1.3 2006-12-20 03:16:10 dschanen Exp $
!-----------------------------------------------------------------------

       module outputnc

!      Description:
!      Functions and subroutines for writing NetCDF files

!      References:
!      <http://www.unidata.ucar.edu/software/netcdf/docs/>
!-----------------------------------------------------------------------

       implicit none

       public :: open_nc, write_nc, close_nc

       private :: define_nc, write_grid, first_write, format_date

       contains
!-----------------------------------------------------------------------
       subroutine open_nc( ncf, unit, fdir, fname, ia, iz, zgrid, 
     .                     day, month, year, rlat, rlon,
     .                     time, dtwrite, nvar )

!      Description:
!      Defines the structure used to reference the file `ncf'

!      References:
!      None
!-----------------------------------------------------------------------
       use netcdf
       use outputfile_class
       use constants, only: fstderr

       implicit none

       ! Input Variables
       character(len=*), intent(in) :: 
     .  fdir,  ! Directory name of file
     .  fname  ! File name

       integer, intent(in) :: 
     .  unit,             ! Ignored; here for compatibility with GrADS writing
     .  day, month, year, ! Time
     .  ia, iz,           ! First and last grid point?
     .  nvar              ! Number of variables

       real, intent(in) :: 
     .  rlat,  ! Latitude                       [degrees_E]
     .  rlon,  ! Longitude                      [degrees_N]
     .  time,  ! Current time                   [s]
     .  dtwrite! Time between write intervals   [s]

       real, dimension(:), intent(in) :: 
     .  zgrid  ! The model grid                 [m]

       ! Input/output Variables
       type (outputfile), intent(inout) :: ncf

       ! Local Variables
       integer :: stat  ! Error status
       integer :: k     ! Array index

       ! Initialization for NetCDF
       ncf%ldefined = .false.

       ! Define file (compatability with GrADS writing)
       ncf%fdir   = fdir
       ncf%fname  = fname
       ncf%ia     = ia
       ncf%iz     = iz
       ncf%day    = day
       ncf%month  = month
       ncf%year   = year
       ncf%rlat   = rlat
       ncf%rlon   = rlon
 
       ncf%dtwrite = dtwrite
       ncf%nvar    = nvar

       ! From open_grads.
       ! Must be in case of a reversed grid as in COAMPS
       if ( ia <= iz ) then
         do k=1,iz-ia+1
           ncf%z(k) = zgrid(ia+k-1)
         enddo
       else ! Always this for HOC
         do k=1,ia-iz+1
           ncf%z(k) = zgrid(ia-k+1)
         enddo
       endif

       ! Create NetCDF dataset: enter define mode  
       stat = nf90_create( path = trim( fdir )//trim( fname )//'.nc', 
     .                     cmode = NF90_CLOBBER, ! overwrite existing file
     .                     ncid = ncf%iounit )  
       if ( stat /= NF90_NOERR ) then
         write(unit=fstderr,fmt=*) "Error opening file: ", 
     .     trim( fdir )//trim( fname )//'.nc',
     .     trim( nf90_strerror( stat ) )
         stop
       endif

       call define_nc( ncf%iounit, ncf%ia, ncf%iz, ncf%z,
     .                 ncf%day, ncf%month, ncf%year, ncf%time,
     .                 ncf%dtwrite, ncf%nvar, ncf%LatDimId, 
     .                 ncf%LongDimId, ncf%AltDimId, ncf%TimeDimId,
     .                 ncf%LatVarId,  ncf%LongVarId, ncf%AltVarId, 
     .                 ncf%TimeVarId )



       return
       end subroutine open_nc

!-----------------------------------------------------------------------

       subroutine write_nc( ncf )

!      Description:
!      Writes some data to the NetCDF dataset, but doesn't close it.
!-----------------------------------------------------------------------

       use netcdf
       use outputfile_class
       use constants, only: fstderr

       implicit none

       ! Input
       type (outputfile), intent(inout) :: ncf    ! The file

       ! Local Variables
       integer, dimension(:), allocatable :: stat ! Error status
       real(kind=8), dimension(1) :: time         ! Time          [s]

       integer :: i ! Array index

       ncf%ntimes = ncf%ntimes + 1

       if ( .not. ncf%ldefined ) then
         call first_write( ncf ) ! finalize the variable definitions
         call write_grid( ncf )  ! define lat., long., and grid
         ncf%ldefined = .true.
       endif

       allocate( stat( ncf%nvar ) )

       time = int( ncf%ntimes * dble( ncf%dtwrite / 60.0 ) ) !  minutes(rounded)
!      time = dble( ncf%ntimes ) * ncf%dtwrite ! seconds

       stat(1) = nf90_put_var( ncid=ncf%iounit, varid=ncf%TimeVarId, 
     .                         values=time(1), start=(/ncf%ntimes/) ) 
       if ( stat(1) /= NF90_NOERR ) then
         stop "time put() failed"
       end if

       do i = 1, ncf%nvar, 1
!        provide values for the variables
!        stat(i) = nf90_put_var( ncid=ncf%iounit, varid=ncf%var(i)%Id, 
!    .          values=reshape( ncf%var(i)%ptr(ncf%ia:ncf%iz), 
!    .                          (/1, 1, ncf%iz, 1/ ) ) ,
!    .          start=(/1,1,1,ncf%ntimes/) )
       ! Work around for a performance issue on pgf90
         stat(i) 
     .   = nf90_put_var( ncid=ncf%iounit, varid=ncf%var(i)%Id, 
     .                   values=ncf%var(i)%ptr(ncf%ia:ncf%iz), 
     .                   start=(/1,1,1,ncf%ntimes/),
     .                   count=(/1,1,ncf%iz,1/) )
       end do ! i=1..nvar

       if ( any (stat /= NF90_NOERR ) ) then 
         do i=1,ncf%nvar,1
           if( stat(i) /= NF90_NOERR ) then
             write(unit=fstderr,fmt=*) ncf%var(i)%name, 
     .         trim( nf90_strerror( stat(i) ) )
           end if
         end do
         stop "nf90_put_var error"
       end if


       deallocate( stat )

       return 
       end subroutine write_nc
!-----------------------------------------------------------------------
!      SUBROUTINE define_nc()

!      Used internally to create a definition for the NetCDF dataset
!-----------------------------------------------------------------------
       subroutine define_nc( ncid, ia, iz, zgrid, day, month, 
     .                       year, time, dtwrite, nvar,
     .                       LatDimId, LongDimId, AltDimId, TimeDimId,
     .                       LatVarId, LongVarId, AltVarId, TimeVarId )
       use netcdf
       use constants, only: fstderr

       implicit none

       ! Constant parameters
       integer, parameter :: 
     .  nlat  = 1,  ! Number of points in the E/W direction
     .  nlong = 1   ! Number of points in the N/S direction

       ! Input Variables
       integer, intent(in) :: 
     .  day, month, year, ! Time of year
     .  ncid,             ! Number used by NetCDF for ref. the file
     .  ia, iz,           ! First and last grid point in the z?
     .  nvar              ! Number of variables

       real, intent(in) :: 
     .  time, dtwrite  ! Current time & write interval   [s]

       real, intent(in), dimension(:) :: 
     .  zgrid   ! HOC grid      [m]

       ! Output Variables
       integer, intent(out) :: 
     .  LatDimId, LongDimId, AltDimId, TimeDimId  ! NetCDF id's for dimensions

       ! NetCDF id's for data (e.g. longitude) associated with each dimension
       integer, intent(out) :: 
     .  LatVarId, LongVarId, AltVarId, TimeVarId  

       ! Local variables
       integer :: stat
       character(len=35) :: TimeUnits

       ! Define the dimensions for the variables
       stat = 
     . nf90_def_dim( ncid, "longitude", nlong, LongDimId )
       if ( stat /= NF90_NOERR ) then
         write(fstderr,*) "Error defining longitude: ",
     .     trim( nf90_strerror( stat ) )
         stop
       endif

       stat = 
     . nf90_def_dim( ncid, "latitude", nlat, LatDimId )
       if ( stat /= NF90_NOERR ) then
         write(fstderr,*) "Error defining latitude: ",
     .     trim( nf90_strerror( stat ) )
         stop
       endif

       stat = 
     . nf90_def_dim( ncid, "altitude", iz, AltDimId )
       if ( stat /= NF90_NOERR ) then
         write(fstderr,*) "Error defining altitude: ",
     .     trim( nf90_strerror( stat ) )
         stop
       endif

       stat = 
     . nf90_def_dim( ncid, "time", NF90_UNLIMITED, TimeDimId ) 
       if ( stat /= NF90_NOERR ) then
         write(fstderr,*) "Error defining time",
     .     trim( nf90_strerror( stat ) )
         stop
       endif

       ! Define the initial variables for the dimensions
       stat = nf90_def_var( ncid, "longitude", NF90_FLOAT,
     .                      (/LongDimId/), LongVarId )

       stat = nf90_def_var( ncid, "latitude", NF90_FLOAT,
     .                      (/LatDimId/), LatVarId )

       stat = nf90_def_var( ncid, "altitude", NF90_FLOAT,
     . (/AltDimId/), AltVarId )

       ! grads2nc stores time as a double prec. value, so we follow that
       stat = nf90_def_var( ncid, "time", NF90_DOUBLE,
     .                      (/TimeDimId/), TimeVarId )

       ! Assign attribute values

       ! Time attribute
       stat = nf90_put_att( ncid, TimeVarId, "cartesian_axis", "T" )
       if ( stat /= NF90_NOERR ) then
         write(fstderr,*) "Error defining time: ",
     .     trim( nf90_strerror( stat ) )
         stop
       endif

       call format_date( day, month, year, time, TimeUnits )

       stat = nf90_put_att( ncid, TimeVarId, "units", TimeUnits )
       if ( stat /= NF90_NOERR ) then
         write(fstderr,*) "Error defining time: ",
     .     trim( nf90_strerror( stat ) )
         stop
       endif

       stat = nf90_put_att( ncid, TimeVarId, "ipositive", 1 )
       if ( stat /= NF90_NOERR ) then
         write(fstderr,*) "Error defining time: ",
     .     trim( nf90_strerror( stat ) )
         stop
       endif

       stat = nf90_put_att( ncid, TimeVarId, "calendar_type", 
     .        "Gregorian" )
       if ( stat /= NF90_NOERR ) then
         write(fstderr,*) "Error defining time",
     .     trim( nf90_strerror( stat ) )
         stop
       endif

       ! Location

       stat = nf90_put_att( ncid, LongVarId, "cartesian_axis", "X" )

       stat = nf90_put_att( ncid, LongVarId, "units", 
     .         "degrees_E" )

       stat = nf90_put_att( ncid, LongVarId, "ipositive", 
     .         1 )

       stat = nf90_put_att( ncid, LatVarId, "cartesian_axis", 
     .         "Y" )

       stat = nf90_put_att( ncid, LatVarId, "units", 
     .         "degrees_N" )

       stat = nf90_put_att( ncid, LatVarId, "ipositive", 
     .         1 )

       ! Altitude
       stat = nf90_put_att( ncid, AltVarId, "cartesian_axis", 
     .         "Z" )

       stat = nf90_put_att( ncid, AltVarId, "units", 
     .         "meters" )

       stat = nf90_put_att( ncid, AltVarId, "positive", 
     .         "up" )

       stat = nf90_put_att( ncid, AltVarId, "ipositive", 
     .         1 )

       return
       end subroutine define_nc

!-----------------------------------------------------------------------
       subroutine close_nc( ncf )

!      Description:
!      Close a previously opened stats file.

!      Notes:
!      I assume nf90_close() exists so that the NetCDF libraries can do a
!      form of buffered I/O, but I don't know the implementation
!      details. -dschanen
!-----------------------------------------------------------------------

       use outputfile_class
       use netcdf, only: nf90_close, nf90_strerror, NF90_NOERR
       use constants, only: fstderr

       implicit none

       ! Input/Output Variables
       type (outputfile), intent(inout) :: ncf

       ! Local Variables
       integer :: stat

       stat = nf90_close( ncf%iounit )
       if ( stat /= NF90_NOERR ) then
         write(fstderr,*) "Error closing file "// 
     .     trim( ncf%fname )//": ", trim( nf90_strerror( stat ) )
         stop
       endif

       return
       end subroutine close_nc

!-----------------------------------------------------------------------
       subroutine first_write( ncf )

!      Description:
!      Used on the first call to write_nc to finalize definitions
!      for the dataset, including the attributes for variable records
!-----------------------------------------------------------------------

       use netcdf
       use outputfile_class
       use constants, only: fstderr

       implicit none

       ! Input/Output Variables
       type (outputfile), intent(inout) :: ncf

       ! Local Variables
       integer, dimension(:), allocatable :: stat

       integer :: i     ! Array index
       logical :: error ! Error stat

       ! Range for NetCDF variables
       real(kind=4), dimension(2) :: var_range

       ! Dimensions for variables
       integer, dimension(4) :: var_dim

!-----------------------------------------------------------------------
!      Typical valid ranges (IEEE 754)

!      real(kind=4): +/- 3.4028235E+38
!      real(kind=8): +/- 1.797693134862316E+308
!      real(kind=16):+/- 1.189731495357231765085759326628007E+4932

!      We use a 4 byte data model for NetCDF and GrADS to save disk space
!-----------------------------------------------------------------------
       var_range(1) = -huge( var_range(1) ) 
       var_range(2) =  huge( var_range(2) )
!      var_range = (/ -1.e31, 1.e31 /)

!      Explanation:  The NetCDF documentation claims the NF90_UNLIMITED
!      variable should be the first dimension, but def_var is somehow
!      inverted and requires the opposite.  After writing, these 
!      dimensions are all in the opposite order of this in the file.
!      -dschanen

       var_dim(1) = ncf%LongDimId
       var_dim(2) = ncf%LatDimId
       var_dim(3) = ncf%AltDimId
       var_dim(4) = ncf%TimeDimId ! The NF90_UNLIMITED dimension

       allocate( stat( ncf%nvar ) )

       error = .false.

       do i = 1, ncf%nvar, 1
!        stat(i) = nf90_def_var( ncf%iounit, trim( ncf%var(i)%name ),
!    .             NF90_FLOAT, (/ncf%TimeDimId, ncf%AltDimId, 
!    .             ncf%LatDimId, ncf%LongDimId/), ncf%var(i)%Id )
         stat(i) = nf90_def_var( ncf%iounit, trim( ncf%var(i)%name ),
     .             NF90_FLOAT, var_dim(:), ncf%var(i)%Id )
         if ( stat(i) /= NF90_NOERR ) then
           write(fstderr,*) "Error defining variable ", 
     .       ncf%var(i)%name //": ", trim( nf90_strerror( stat(i) ) )
           error = .true.
         endif

         stat(i) = nf90_put_att( ncf%iounit, ncf%var(i)%Id,
     .             "valid_range", var_range(1:2) )
         if ( stat(i) /= NF90_NOERR ) then
           write(fstderr,*) "Error defining valid range",
     .       trim( nf90_strerror( stat(i) ) )
           error = .true.
         endif

         stat(i) = nf90_put_att( ncf%iounit, ncf%var(i)%Id, "title", 
     .             trim( ncf%var(i)%description ) )
         if ( stat(i) /= NF90_NOERR ) then
           write(fstderr,*) "Error in description",
     .       trim( nf90_strerror( stat(i) ) )
           error = .true.
         endif

         stat(i) = nf90_put_att( ncf%iounit, ncf%var(i)%Id, "units", 
     .             trim( ncf%var(i)%units ) )
         if ( stat(i) /= NF90_NOERR ) then
           write(fstderr,*) "Error in units",
     .       trim( nf90_strerror( stat(i) ) )
           error = .true.
         endif
       end do

       if ( error ) stop "Error in definition"

       stat(1) = nf90_enddef( ncf%iounit ) ! end definitions
       if ( stat(1) /= NF90_NOERR ) then
         write(fstderr,*) "Error finalizing definitions",
     .              trim( nf90_strerror( stat(1) ) )
         stop
       endif

       deallocate( stat )

       return
       end subroutine first_write

!-----------------------------------------------------------------------
       subroutine write_grid( ncf )

!      Description:
!      Writes inforation about latitude, longitude and the grid
!-----------------------------------------------------------------------

       use netcdf
       use outputfile_class
       use constants, only: fstderr

       implicit none

!      Constants
!      real, parameter, dimension(1) :: deg_east  = 0.0
!      real, parameter, dimension(1) :: deg_north = 0.0

!      Input
       type (outputfile), intent(inout) :: ncf

       integer :: stat

       stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%AltVarId, 
     .                      values=ncf%z(ncf%ia:ncf%iz) )
       if ( stat /= NF90_NOERR ) then
         write(fstderr,*) "Error entering grid: ", 
     .     trim( nf90_strerror( stat ) )
         stop
       endif

       stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%LongVarId, 
     .                      values=ncf%rlon ) 
       if ( stat /= NF90_NOERR ) then
         write(fstderr,*) "Error entering longitude: ", 
     .     trim( nf90_strerror( stat ) )
         stop
       endif

       stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%LatVarId, 
     .                      values=ncf%rlat ) 
       if ( stat /= NF90_NOERR ) then
         write(fstderr,*) "Error entering latitude: ", 
     .     trim( nf90_strerror( stat ) )
         stop
       endif

       return
       end subroutine write_grid

!-----------------------------------------------------------------------

       subroutine format_date( iday, imonth, iyear, st_time, date )

!      Description:
!      Put the model date in a format that udunits and NetCDF can easily
!      handle.  GrADSnc is dumb and apparently cannot handle time
!      intervals < 1 minute.

!      Notes:
!      Adapted from the original GrADS version written by Chris Golaz.
!      Uses Fortran `internal' files to write the string output.
!-----------------------------------------------------------------------
       implicit none

       ! Input Variables
       integer, intent(in) :: 
     . iday, imonth, iyear  ! Integer for day, month and year.

       real, intent(in) :: st_time ! Start time

       ! Output Variables
       character(len=35), intent(out) :: date

!      date(1:14) = "minutes since "
       date = "minutes since YYYY-MM-DD HH:MM:00.0"
!      date = "seconds since YYYY-MM-DD HH:MM:00.0"
       write(date(15:18),'(i4.4)') iyear
!      write(date(19),'(a1)') '-'
       write(date(20:21),'(i2.2)') imonth
!      write(date(22),'(a1)') '-'
       write(date(23:24),'(i2.2)') iday
!      write(date(25),'(a1)') ' '
       write(date(26:27),'(i2.2)') int(st_time/3600.)
!      write(date(28),'(a1)') ":"
       write(date(29:30),'(i2.2)') int( mod(int(st_time),3600) / 60. )
!      write(date(30:35),'(a5)') ":00.0"

       return
       end subroutine format_date

       end module outputnc
