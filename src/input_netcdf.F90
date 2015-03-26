!-----------------------------------------------------------------------
! $Id$

module input_netcdf
#ifdef NETCDF
  use constants_clubb, only:  & 
    fstdout,    & ! Variable(s) 
    fstderr,    &
    sec_per_day,&
    sec_per_hr, &
    sec_per_min 

  implicit none

  private ! Default Scope

  public :: get_netcdf_var, open_netcdf_read, & 
    close_netcdf_read

  contains

!-----------------------------------------------------------------------
  subroutine open_netcdf_read( test_variable, path, ncf, l_error )

! Description:
!   Open a netCDF data set in read-only mode
!-----------------------------------------------------------------------
    use netcdf, only: &
      nf90_inq_varid, & ! Function(s)
      nf90_get_var, &
      nf90_get_att, &
      nf90_inquire_variable, &
      nf90_open, &
      nf90_strerror, &
      nf90_inquire_dimension

    use netcdf, only: &
      NF90_NOWRITE,      & ! Constant(s) 
      NF90_NOERR,        & 
      NF90_MAX_VAR_DIMS, &
      NF90_FLOAT, &
      NF90_DOUBLE, &
      NF90_MAX_NAME

    use stat_file_module, only: stat_file ! Type

    use clubb_precision, only: &
      time_precision, &
      core_rknd

    use calendar, only: &
      gregorian2julian_date

    use clubb_model_settings, only: &
      day, &
      month, &
      year

    implicit none

    ! Input Variables
    character(len=*), intent(in) :: &
      test_variable ! Variable used to determine the dimensions in the file

    character(len=*), intent(in) ::  & 
      path ! The file name

    ! Output Variables
    type(stat_file), intent(out) :: ncf ! Derived type for the netCDF file

    logical, intent(out) :: l_error

    ! Local Variables
    integer :: i, ierr, itmp, varid, xtype

    integer :: xdim, ydim, ndims
    
    real( kind = time_precision ) :: hours, minutes, seconds, multiplier, delta_d

    integer :: length, netcdf_day, netcdf_month, netcdf_year, &
               clubb_idate, netcdf_idate
    
    real( kind = core_rknd ), dimension(2) :: write_times

    character(len=80) :: time

    character(len=NF90_MAX_NAME) :: zname, time_name, dim_name

    integer, dimension(NF90_MAX_VAR_DIMS) :: dimIds

    integer :: clubb_day, clubb_month, clubb_year

!-----------------------------------------------------------------------

    ! ---- Begin Code ----
    ! Initialize clubb_day, month, and year
    clubb_day = day
    clubb_month = month
    clubb_year = year

    ! Initialize l_error to false
    l_error = .false.

    ierr = nf90_open( path=trim( path ), mode=NF90_NOWRITE, ncid=ncf%iounit )

    if ( ierr /= NF90_NOERR ) then
      write(fstderr,*) "Error opening netCDF file: "//trim( path )
      write(fstderr,*) nf90_strerror( ierr )
      l_error = .true.
      return 
    end if

    ierr = nf90_inq_varid( ncid=ncf%iounit, name=test_variable, varid=varid )

    if ( ierr /= NF90_NOERR ) then
      write(fstderr,*) nf90_strerror( ierr )
      l_error = .true.
      return
    end if

    ierr = nf90_inquire_variable( ncid=ncf%iounit, varid=varid, xtype=xtype, &
                                  ndims=ndims, dimIds=dimIds )

    if ( ierr /= NF90_NOERR .or. ndims /= 4 .or. &
                            .not. (xtype == NF90_DOUBLE .or. xtype == NF90_FLOAT) ) then
            write(fstderr,*) "input_netcdf.open_netcdf_read : The netCDF data doesn't "// &
            "conform to expected precision, shape, or dimensions"
      l_error = .true.
      return

    end if

    xdim = -1
    ydim = -1
    ncf%iz = -1
    ncf%ntimes = -1

    do i = 1, ndims

      ierr = nf90_inquire_dimension( ncid=ncf%iounit, dimId=dimIds(i), name=dim_name, len=itmp )
      if ( ierr /= NF90_NOERR ) then
        write(fstderr,*) nf90_strerror( ierr )
        l_error = .true.
        return
      end if

      select case ( trim( dim_name ) )
      case ( "Z", "z", "altitude", "height" )
        ncf%ia = 1
        ncf%iz = itmp
        ncf%AltDimId = dimIds(i)
        zname = dim_name
      case ( "X", "x", "longitude" )
        xdim = itmp
        ncf%LongDimId = dimIds(i)
      case ( "Y", "y", "latitude" )
        ydim = itmp
        ncf%LatDimId = dimIds(i)
      case ( "T", "t", "time" )
        ncf%ntimes = itmp
        ncf%TimeDimId = dimIds(i)
        time_name = dim_name
      case default
        l_error = .true.
        return
      end select

    end do

    if ( ydim /= 1 .or. xdim /= 1 ) then
      write(fstderr,*) "input_netcdf.open_netcdf_read : The netCDF data doesn't "// &
        "conform to the expected X or Y dimension"
      l_error = .true.
      return
    end if
      

    if ( .not. l_error ) then

      ! Allocate altitudes
      allocate( ncf%z(ncf%iz) )

      ! Determine the altitudes
      ierr = nf90_inq_varid( ncid=ncf%iounit, name=trim( zname ), varid=ncf%AltVarId )
      if ( ierr /= NF90_NOERR ) then
        write(fstderr,*) nf90_strerror( ierr )
        l_error = .true.
        return
      end if

!     ierr = nf90_inquire_variable( ncid=ncf%iounit, varid=ncf%AltVarId, & ! In
!                                   ndims=itmp ) ! Out

!     if ( ierr /= NF90_NOERR ) then
!       write(fstderr,*) nf90_strerror( ierr )
!       l_error = .true.
!       return
!     end if

      ierr = nf90_get_var( ncid=ncf%iounit, varid=ncf%AltVarId, & ! In
                           start=(/1/), count=(/ncf%iz/), & ! In
                           values=ncf%z(:) ) ! Out

      if ( ierr /= NF90_NOERR ) then
        write(fstderr,*) nf90_strerror( ierr )
        l_error = .true.
        return
      end if

      ierr = nf90_inq_varid( ncid=ncf%iounit, name=trim( time_name ), varid=ncf%TimeVarId )
      if ( ierr /= NF90_NOERR ) then
        write(fstderr,*) nf90_strerror( ierr )
        l_error = .true.
        return
      end if

!     ierr = nf90_inquire_variable( ncid=ncf%iounit, varid=ncf%TimeVarId, & ! In
!                                   ndims=itmp ) ! Out

      ierr = nf90_get_att( ncid=ncf%iounit, varid=ncf%TimeVarId, name="units", & ! In
                           values=time ) ! Out

      if ( ierr /= NF90_NOERR ) then
        write(fstderr,*) nf90_strerror( ierr )
        l_error = .true.
        return
      end if

      ! Read starting time from the "units" attribute of the netcdf file
      time = trim( time )
      length = len( trim( time ) )
      if ( length < 21 ) then
        write(fstderr,*) "The NetCDF file does not have a proper time unit &
                         &specification. The ""units"" attribute for the &
                         &time variable must be in the form:"
        write(fstderr,*) "TIMEUNITS since YYYY-MM-DD HH:MM:SS.S"
        l_error = .true.
        return
      end if

      read(time( length-20:length-17), *) netcdf_year
      read(time( length-15:length-14), *) netcdf_month
      read(time( length-12:length-11), *) netcdf_day

      read(time( length-10:length-8 ), *) hours
      read(time( length-6:length-5 ), *) minutes
      read(time( length-3:length - 2 ), *) seconds

      ncf%year = netcdf_year
      ncf%month = netcdf_month
      ncf%day = netcdf_day

      ! Compute delta_d, the difference in days between the netcdf file's start
      ! date and CLUBB's start date (netcdf_date - clubb_date). This difference
      ! should be added to ncf%time to fix an inconsistency that is caused when
      ! the netcdf date differs from CLUBB's date.
      ! July 2013 Eric Raut
      clubb_idate = gregorian2julian_date( clubb_day, clubb_month, clubb_year)
      netcdf_idate = gregorian2julian_date( netcdf_day, netcdf_month, netcdf_year)
      delta_d = real( netcdf_idate - clubb_idate, kind=time_precision )

      ncf%time = (hours * real(sec_per_hr,kind=time_precision)) + &
                 (minutes * real(sec_per_min,kind=time_precision)) &
                + seconds + (delta_d * real(sec_per_day,kind=time_precision))

      ! Get rid of compiler warning
      multiplier = 0.0_time_precision

      ! Figure out what units Time is in so dtwrite can be set correctly
      select case ( time( 1:index ( time, ' ' ) ) )
      case ( "hours" )
        multiplier = real(sec_per_hr, kind=time_precision)

      case ( "minutes" )
        multiplier = real(sec_per_min, kind=time_precision)

      case ( "seconds" )
        multiplier = 1._time_precision

      case default
        multiplier = 0._time_precision  ! Initialized to eliminate g95 compiler warning -meyern
        l_error = .true.
        return

      end select

      ierr = nf90_get_var( ncid=ncf%iounit, varid=ncf%TimeVarId, & ! In
                           start=(/1/), count=(/2/), & ! In
                           values=write_times ) ! Out

      if ( ierr /= NF90_NOERR ) then
        write(fstderr,*) nf90_strerror( ierr )
        l_error = .true.
        return
      end if

      ncf%dtwrite =  (write_times(2) - write_times(1)) * real(multiplier,kind=core_rknd)

    end if

    return
  end subroutine open_netcdf_read

!----------------------------------------------------------------------
  subroutine get_netcdf_var( ncf, varname, itime, l_convert_to_MKS, x, l_error )

! Description:
!   Read in values from a netCDF file.

! Assumptions:
!   We assume a NF90_FLOAT the same as a kind=4 real in Fortran, and
!   That the variables will obey COARDS conventions and have 4 dimensions.
!----------------------------------------------------------------------

    use netcdf, only: &
      nf90_inq_varid, & ! Function(s)
      nf90_get_var, &
      nf90_get_att, &
      nf90_inquire_variable, &
      nf90_strerror

    use netcdf, only: &
      NF90_NOERR, & ! Constant(s)
      NF90_MAX_VAR_DIMS, &
      NF90_DOUBLE, &
      NF90_FLOAT

    use stat_file_module, only: &
      stat_file ! Type(s)
    use constants_clubb, only: &
      g_per_kg, &
      sec_per_day

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    type (stat_file), intent(in) :: &
      ncf ! NetCDF file structure

    character(len=*), intent(in) ::  & 
      varname ! The variable name as it occurs in the control file

    integer, intent(in) :: & 
      itime ! Obtain variable varname at time itime [m]

    logical, intent(in) :: &
      l_convert_to_MKS ! Convert variables to MKS upon reading

    ! Output
    real( kind = core_rknd ), dimension(:), intent(out) ::  & 
      x ! Result variable
     
    logical, intent(out) :: l_error

    ! Local Variables

    integer :: varid, ierr, xtype, ndims

    ! Local Variables
    integer, dimension(NF90_MAX_VAR_DIMS) :: dimIds

    character(len=10) :: &
      units ! The units on the variable

    real( kind = core_rknd ), dimension(:,:,:,:), allocatable :: & 
      x4 ! The variable from the file


    ! ---- Begin Code ----


    ! Initialize l_error to false
    l_error = .false.

    ierr = nf90_inq_varid( ncid=ncf%iounit, name=varname, varid=varid )

    if ( ierr /= NF90_NOERR ) then
      write(fstderr,*) nf90_strerror( ierr )
      l_error = .true.
      return
    end if

    ierr = nf90_inquire_variable( ncid=ncf%iounit, varid=varid, xtype=xtype, &
                                  ndims=ndims, dimIds=dimIds )

    if ( ierr /= NF90_NOERR .or. ndims /= 4 .or. &
                            .not. (xtype == NF90_DOUBLE .or. xtype == NF90_FLOAT) ) then
            write(fstderr,*) "input_netcdf.get_var : The netCDF data doesn't "// &
            "conform to expected precision, shape, or dimensions"
      l_error = .true.
      return

    else

      if ( itime > ncf%ntimes ) then
        write(fstderr,*) "input_netcdf.get_var: The time specified exceeds the netCDF data"
        l_error = .true.
        return
      end if

      allocate( x4(1,1,ncf%iz,1) )

    end if

    ! Obtain a vertical column at time itime
    ierr = nf90_get_var( ncid=ncf%iounit, varid=varid, values=x4, &
      start=(/1,1,1,itime/), count=(/1,1,ncf%iz,1/) )

    if ( ierr /= NF90_NOERR ) then
      write(fstderr,*) nf90_strerror( ierr )
      l_error = .true.
      return
    end if

    x = real( x4(1,1,:,1), kind = core_rknd )

    if ( l_convert_to_MKS ) then

      ierr = nf90_get_att( ncid=ncf%iounit, varid=varid, name='units', values=units )

      if ( ierr /= NF90_NOERR ) then
        write(fstderr,*) nf90_strerror( ierr )
        l_error = .true.
        return

      else
        select case ( trim( units ) )
        case ( "g/kg" ) 
          x = x / g_per_kg

        case ( "W/m2" ) 
           l_error = .true.
           write(fstderr,*) "get_netcdf_var: Unable to convert variables of this type to MKS units"
           return

        case ( "K/day" ) 
          x = x / sec_per_day

        case default
          ! Do nothing

        end select ! units

      end if ! ierr /= NF90_NOERR

    end if ! l_convert_to_MKS

    return
  end subroutine get_netcdf_var


!-----------------------------------------------------------------------
  subroutine close_netcdf_read( ncf )

! Description:
!   Close a previously opened netCDF file
!-----------------------------------------------------------------------
    use netcdf, only: nf90_close, nf90_strerror ! Function(s)

    use netcdf, only: NF90_NOERR ! Constant(s)

    use stat_file_module, only: &
      stat_file ! Type(s)

    implicit none

    ! Input Variables
    type(stat_file), intent(in) :: &
      ncf

    integer :: ierr
!-----------------------------------------------------------------------
    ! Close file
    ierr = nf90_close( ncid=ncf%iounit )

    if ( ierr /= NF90_NOERR ) then
      write(fstderr,*) "Error closing a netCDF file"
      write(fstderr,*) nf90_strerror( ierr )
      stop "Fatal error."
    end if

    return
  end subroutine close_netcdf_read

#endif /* NETCDF */

end module input_netcdf
