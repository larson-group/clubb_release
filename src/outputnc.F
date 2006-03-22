!-----------------------------------------------------------------------
! $Id: outputnc.F,v 1.5 2006-03-22 16:56:59 dschanen Exp $
!      MODULE outputnc

!      functions and subroutines for writing NetCDF files
!-----------------------------------------------------------------------
       module outputnc

       implicit none

       public open_nc, write_nc, close_nc
       private define_nc, write_grid, first_write, format_date

       contains
!-----------------------------------------------------------------------
!      SUBROUTINE open_nc()

!      Defines the structure used to reference the file `ncf'
!-----------------------------------------------------------------------
       subroutine open_nc( ncf, unit, fdir, fname, ia, iz, zgrid, 
     .                     day, month, year, time, dtwrite, nvar )
       use netcdf
       use outputfile_class

       implicit none

!      Input
       type (outputfile), intent(inout) :: ncf

       integer, intent(in) :: unit ! (ignored)

       character(len=*), intent(in) :: fdir, fname
       integer, intent(in) :: day, month, year

       integer, intent(in) :: ia, iz
       integer, intent(in) :: nvar

       real, intent(in)    :: zgrid(:)
       real, intent(in)    :: time, dtwrite

!      Internal
       integer stat
       integer k

!      Initialization for NetCDF
       ncf%ldefined = .false.

!      Define file (compatability)

       ncf%fdir   = fdir
       ncf%fname  = fname
       ncf%ia     = ia
       ncf%iz     = iz
       ncf%day    = day
       ncf%month  = month
       ncf%year   = year
 
       ncf%dtwrite = dtwrite
       ncf%nvar    = nvar

!      I'm not entirely sure the purpose of this code -dschanen
       if ( ia <= iz ) then
         do k=1,iz-ia+1
           ncf%z(k) = zgrid(ia+k-1)
         enddo
       else
         do k=1,ia-iz+1
           ncf%z(k) = zgrid(ia-k+1)
         enddo
       endif

!      Create NetCDF dataset: enter define mode  
       stat = nf90_create( path = trim( fdir )//trim( fname )//'.nc', 
     .                     cmode = NF90_CLOBBER, ! overwrite existing file
     .                     ncid = ncf%iounit )  
       if ( stat /= NF90_NOERR ) then
         write(*,*) "Error opening file ", 
     .              trim( fdir )//trim( fname)//'.nc',
     .              trim( nf90_strerror(stat) )
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
!      SUBROUTINE write_nc()

!      Writes some data to the NetCDF dataset, but doesn't close it.
!-----------------------------------------------------------------------

       subroutine write_nc( ncf )

       use netcdf
       use outputfile_class

       implicit none

!      Input
       type (outputfile), intent(inout) :: ncf

!      Internal
       integer, dimension(:), allocatable :: stat
       real(kind=8), dimension(1) :: time

       integer i

       ncf%ntimes = ncf%ntimes + 1

       if ( .not. ncf%ldefined ) then
         call first_write( ncf ) ! finalize the variable definitions
         call write_grid( ncf )  ! define lat., long., and grid
         ncf%ldefined = .true.
       endif

       allocate( stat( ncf%nvar ) )
       time = int( ncf%ntimes * dble(ncf%dtwrite / 60.0) ) !  minutes(rounded)
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
         stat(i) = nf90_put_var( ncid=ncf%iounit, varid=ncf%var(i)%Id, 
     .          values=ncf%var(i)%ptr(ncf%ia:ncf%iz), 
     .          start=(/1,1,1,ncf%ntimes/),
     .          count=(/1,1,ncf%iz,1/) )
       end do
       if ( any (stat /= NF90_NOERR ) ) then 
         do i=1,ncf%nvar,1
           if( stat(i) /= NF90_NOERR ) then
             write(*,*) ncf%var(i)%name, trim( nf90_strerror(stat(i)) )
           end if
         end do
         stop "nf90_put_var error"
       endif


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

       implicit none

!      Constants
       integer, parameter :: nlong = 1
       integer, parameter :: nlat  = 1

!      Input
       integer, intent(in) :: day, month, year

       integer, intent(in) :: ncid
       integer, intent(in) :: ia, iz
       integer, intent(in) :: nvar

       real, intent(in)    :: zgrid(:)
       real, intent(in)    :: time, dtwrite

!      Output
       integer, intent(out) :: LatDimId, LongDimId, AltDimId, TimeDimId
       integer, intent(out) :: LatVarId, LongVarId, AltVarId, TimeVarId

!      Internal
       integer stat
       character(len=35) TimeUnits

       ! Define the dimensions for the variables
       stat = 
     . nf90_def_dim( ncid, "longitude", nlong, LongDimId )
       if ( stat /= NF90_NOERR ) then
         write(*,*) "Error defining longitude",
     .               trim( nf90_strerror(stat) )
         stop
       endif

       stat = 
     . nf90_def_dim( ncid, "latitude", nlat, LatDimId )
       if ( stat /= NF90_NOERR ) then
         write(*,*) "Error defining latitude",
     .               trim( nf90_strerror(stat) )
         stop
       endif

       stat = 
     . nf90_def_dim( ncid, "altitude", iz, AltDimId )
       if ( stat /= NF90_NOERR ) then
         write(*,*) "Error defining altitude",
     .               trim( nf90_strerror(stat) )
         stop
       endif

       stat = 
     . nf90_def_dim( ncid, "time", NF90_UNLIMITED, TimeDimId ) 
       if ( stat /= NF90_NOERR ) then
         write(*,*) "Error defining time",
     .               trim( nf90_strerror(stat) )
         stop
       endif

       ! Define the initial variables for the dimensions
       stat = nf90_def_var( ncid, "longitude", NF90_FLOAT,
     . (/LongDimId/), LongVarId )

       stat = nf90_def_var( ncid, "latitude", NF90_FLOAT,
     . (/LatDimId/), LatVarId )

       stat = nf90_def_var( ncid, "altitude", NF90_FLOAT,
     . (/AltDimId/), AltVarId )

       ! grads2nc stores time as a double prec. value, so we follow that
       stat = nf90_def_var( ncid, "time", NF90_DOUBLE,
     . (/TimeDimId/), TimeVarId )

       ! Assign attribute values

       ! Time
       stat = nf90_put_att( ncid, TimeVarId, "cartesian_axis", 
     .        "T" )
       if ( stat /= NF90_NOERR ) then
         write(*,*) "Error defining time",
     .               trim( nf90_strerror(stat) )
         stop
       endif

       call format_date( day, month, year, time, TimeUnits )

       stat = nf90_put_att( ncid, TimeVarId, "units", 
     .                      TimeUnits )
       if ( stat /= NF90_NOERR ) then
         write(*,*) "Error defining time",
     .               trim( nf90_strerror(stat) )
         stop
       endif

       stat = nf90_put_att( ncid, TimeVarId, "ipositive", 
     .         1 )
       if ( stat /= NF90_NOERR ) then
         write(*,*) "Error defining time",
     .               trim( nf90_strerror(stat) )
         stop
       endif

       stat = nf90_put_att( ncid, TimeVarId, "calendar_type", 
     .        "Gregorian" )
       if ( stat /= NF90_NOERR ) then
         write(*,*) "Error defining time",
     .               trim( nf90_strerror(stat) )
         stop
       endif

       ! Location

       stat = nf90_put_att( ncid, LongVarId, "cartesian_axis", 
     .         "X" )

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
!      SUBROUTINE close_nc()

!      I assume nf90_close exists so that the NetCDF libraries can do a
!      form of buffered I/O, but I don't know the implementation details
!-----------------------------------------------------------------------
       subroutine close_nc( ncf )

       use outputfile_class
       use netcdf

       implicit none

!      Input
       type (outputfile), intent(inout) :: ncf

!      Internal
       integer stat

       stat = nf90_close( ncf%iounit )
       if ( stat /= NF90_NOERR ) then
         write(*,*) "Error closing file:"//trim( ncf%fname ),
     .               trim( nf90_strerror(stat) )
         stop
       endif

       return
       end subroutine close_nc

!-----------------------------------------------------------------------
!      SUBROUTINE first_write()

!      Used on the first call to write_nc to finalize definitions
!      for the dataset, including the attributes for variable records
!-----------------------------------------------------------------------
       subroutine first_write( ncf )

       use netcdf
       use outputfile_class

       implicit none

!      Input
       type (outputfile), intent(inout) :: ncf

!      Internal
       integer, dimension(:), allocatable :: stat
       integer i
       logical error

!      Range for netcdf variables
       real(kind=4), dimension(2) :: var_range

!      Dimensions for variables
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
           write(*,*) "Error defining variable ", 
     .                ncf%var(i)%name,
     .                trim( nf90_strerror(stat(i)) )
           error = .true.
         endif

         stat(i) = nf90_put_att( ncf%iounit, ncf%var(i)%Id,
     .             "valid_range", var_range(1:2) )
         if ( stat(i) /= NF90_NOERR ) then
           write(*,*) "Error defining valid range",
     .                trim( nf90_strerror(stat(i)) )
           error = .true.
         endif

         stat(i) = nf90_put_att( ncf%iounit, ncf%var(i)%Id, "title", 
     .             trim( ncf%var(i)%description ) )
         if ( stat(i) /= NF90_NOERR ) then
           write(*,*) "Error in description",
     .                trim( nf90_strerror(stat(i)) )
           error = .true.
         endif

         stat(i) = nf90_put_att( ncf%iounit, ncf%var(i)%Id, "units", 
     .             trim( ncf%var(i)%units ) )
         if ( stat(i) /= NF90_NOERR ) then
           write(*,*) "Error in units",
     .                trim( nf90_strerror(stat(i)) )
           error = .true.
         endif
       end do

       if ( error ) stop "Error in definition"

       stat(1) = nf90_enddef( ncf%iounit ) ! end definitions
       if ( stat(1) /= NF90_NOERR ) then
         write(*,*) "Error finalizing definitions",
     .              trim( nf90_strerror(stat(1)) )
         stop
       endif

       deallocate( stat )

       return
       end subroutine first_write

!-----------------------------------------------------------------------
!      SUBROUTINE write_grid( )

!      Writes inforation about latitude, longitude and the grid
!-----------------------------------------------------------------------
       subroutine write_grid( ncf )

       use netcdf
       use outputfile_class

       implicit none

!      Constants
       real, parameter, dimension(1) :: deg_east  = 0.0
       real, parameter, dimension(1) :: deg_north = 0.0

!      Input
       type (outputfile), intent(inout) :: ncf

       integer stat

       stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%AltVarId, 
     .                      values=ncf%z(ncf%ia:ncf%iz) )
       if ( stat /= NF90_NOERR ) then
         write(*,*) "Error entering grid", 
     .              trim( nf90_strerror(stat) )
         stop
       endif

       stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%LongVarId, 
     .                      values=deg_north(1:1) ) 
       if ( stat /= NF90_NOERR ) then
         write(*,*) "Error entering long.", 
     .              trim( nf90_strerror(stat) )
         stop
       endif

       stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%LatVarId, 
     .                      values=deg_east(1:1) ) 
       if ( stat /= NF90_NOERR ) then
         write(*,*) "Error entering lat.", 
     .              trim( nf90_strerror(stat) )
         stop
       endif

       return
       end subroutine write_grid

!-----------------------------------------------------------------------
!      SUBROUTINE format_date()

!      Adapted from the original GrADS version written by Chris Golaz
!-----------------------------------------------------------------------
       subroutine format_date( iday, imonth, iyear, st_time, date )
       implicit none

       integer, intent(in) :: iday, imonth, iyear
       real, intent(in)    :: st_time
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
