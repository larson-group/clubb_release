! $Id$
!-----------------------------------------------------------------------
module output_netcdf
#ifdef NETCDF

! Description:
!   Functions and subroutines for writing NetCDF files

! References:
!   <http://www.unidata.ucar.edu/software/netcdf/docs/>
!-----------------------------------------------------------------------

  implicit none

  public :: open_netcdf, write_netcdf, close_netcdf

  private :: define_netcdf, write_grid, first_write, format_date

  private ! Default scope

  contains
!-----------------------------------------------------------------------
  subroutine open_netcdf( unit, fdir, fname, ia, iz, zgrid,  & 
                          day, month, year, rlat, rlon, & 
                          time, dtwrite, nvar, ncf )

! Description:
!   Defines the structure used to reference the file `ncf'

! References:
!   None
!-----------------------------------------------------------------------
    use netcdf, only: & 
      NF90_CLOBBER, & ! Variable(s)
      NF90_NOERR,   & 
      nf90_create,  & ! Procedure
      nf90_strerror

    use output_file_module, only: & 
      outputfile ! Type

    use stats_precision, only:  & 
      time_precision ! Variable(s)

    use constants, only:  & 
        fstderr ! Variable(s)

    implicit none

    ! Input Variables
    character(len=*), intent(in) ::  & 
     fdir,   & ! Directory name of file
     fname     ! File name

    integer, intent(in) ::  & 
     unit,              & ! Ignored; here for compatibility with GrADS writing
     day, month, year,  & ! Time
     ia, iz,            & ! First and last grid point?
     nvar                 ! Number of variables

    real, intent(in) ::  & 
     rlat,   & ! Latitude                        [degrees_E]
     rlon      ! Longitude                       [degrees_N]

    real(kind=time_precision), intent(in) :: & 
     dtwrite ! Time between write intervals   [s]

    real(kind=time_precision), intent(in) ::  & 
     time   ! Current time                    [s]

    real, dimension(:), intent(in) ::  & 
     zgrid  ! The model grid                  [m]

    ! Input/output Variables
    type (outputfile), intent(inout) :: ncf

    ! Local Variables
    integer :: stat  ! Error status
    integer :: k     ! Array index
   
    ! Make the compiler warning about unused variables go away
    if ( .false. ) k = unit 

    ! Initialization for NetCDF
    ncf%l_defined = .false.

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
    ncf%time   = time

    ncf%dtwrite = dtwrite
    ncf%nvar    = nvar

   ! From open_grads.
    ! This probably for the case of a reversed grid as in COAMPS
    if ( ia <= iz ) then
      do k=1,iz-ia+1
        ncf%z(k) = zgrid(ia+k-1)
      end do
    else ! Always this for CLUBB
      do k=1,ia-iz+1
        ncf%z(k) = zgrid(ia-k+1)
      end do
    end if

    ! Create NetCDF dataset: enter define mode
    stat = nf90_create( path = trim( fdir )//trim( fname )//'.nc',  & 
                        cmode = NF90_CLOBBER,  & ! overwrite existing file
                        ncid = ncf%iounit )
    if ( stat /= NF90_NOERR ) then
      write(unit=fstderr,fmt=*) "Error opening file: ",  & 
        trim( fdir )//trim( fname )//'.nc', & 
        trim( nf90_strerror( stat ) )
      stop
    end if

    call define_netcdf( ncf%iounit, ncf%iz, ncf%day, ncf%month, ncf%year, ncf%time, &
                        ncf%LatDimId, ncf%LongDimId, ncf%AltDimId, ncf%TimeDimId, & 
                        ncf%LatVarId, ncf%LongVarId, ncf%AltVarId, ncf%TimeVarId )



    return
  end subroutine open_netcdf

!-----------------------------------------------------------------------

  subroutine write_netcdf( ncf )

!      Description:
!      Writes some data to the NetCDF dataset, but doesn't close it.
!-----------------------------------------------------------------------

    use netcdf, only: & 
        NF90_NOERR,  & ! Variable(s)
        nf90_put_var,  & ! Procedure
        nf90_strerror

    use output_file_module, only: & 
        outputfile ! Variable
    use constants, only:  & 
        fstderr ! Variable

    implicit none

! Input
    type (outputfile), intent(inout) :: ncf    ! The file

! Local Variables
    integer, dimension(:), allocatable :: stat ! Error status
    real(kind=8), dimension(1) :: time         ! Time          [s]

    integer :: i ! Array index

    ncf%ntimes = ncf%ntimes + 1

    if ( .not. ncf%l_defined ) then
      call first_write( ncf ) ! finalize the variable definitions
      call write_grid( ncf )  ! define lat., long., and grid
      ncf%l_defined = .true.
    endif

    allocate( stat( ncf%nvar ) )

    time = nint( ncf%ntimes * dble( ncf%dtwrite / 60.0 ) ) !  minutes(rounded)
!      time = dble( ncf%ntimes ) * ncf%dtwrite ! seconds

    stat(1) = nf90_put_var( ncid=ncf%iounit, varid=ncf%TimeVarId,  & 
                            values=time(1), start=(/ncf%ntimes/) )
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
      stat(i)  & 
      = nf90_put_var( ncid=ncf%iounit, varid=ncf%var(i)%Id,  & 
                      values=ncf%var(i)%ptr(ncf%ia:ncf%iz),  & 
                      start=(/1,1,1,ncf%ntimes/), & 
                      count=(/1,1,ncf%iz,1/) )
    end do ! i=1..nvar

    if ( any (stat /= NF90_NOERR ) ) then
      do i=1,ncf%nvar,1
        if( stat(i) /= NF90_NOERR ) then
          write(unit=fstderr,fmt=*) ncf%var(i)%name,  & 
            trim( nf90_strerror( stat(i) ) )
        end if
      end do
      stop "nf90_put_var error"
    end if


    deallocate( stat )

    return
  end subroutine write_netcdf

!-----------------------------------------------------------------------
  subroutine define_netcdf( ncid, iz, day, month, year, time, & 
                            LatDimId, LongDimId, AltDimId, TimeDimId, & 
                            LatVarId, LongVarId, AltVarId, TimeVarId )

! Description:
!   Used internally to create a definition for the NetCDF dataset
!
! References:
!   None
!-----------------------------------------------------------------------
    use netcdf, only: & 
      NF90_NOERR,   & ! Constants
      NF90_FLOAT, & 
      NF90_DOUBLE, & 
      NF90_UNLIMITED

    use netcdf, only: & 
      nf90_def_dim,  & ! Functions
      nf90_strerror, & 
      nf90_def_var, & 
      nf90_put_att

    use stats_precision, only:  & 
      time_precision ! Variable(s)

    use constants, only:  & 
      fstderr ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter ::  & 
      nlat  = 1,   & ! Number of points in the N/S direction
      nlong = 1      ! Number of points in the E/W direction

    ! Input Variables
    integer, intent(in) ::  & 
      day, month, year,  & ! Time of year
      ncid,              & ! Number used by NetCDF for ref. the file
      iz                   ! Dimension in z

    real(kind=time_precision), intent(in) ::  & 
      time    ! Current model time [s]

    ! Output Variables
    integer, intent(out) ::  & 
      LatDimId, LongDimId, AltDimId, TimeDimId  ! NetCDF id's for dimensions

    ! NetCDF id's for data (e.g. longitude) associated with each dimension
    integer, intent(out) ::  & 
      LatVarId, LongVarId, AltVarId, TimeVarId

    ! Local variables
    integer :: stat
    character(len=35) :: TimeUnits

    ! ---- Begin Code ----

    ! Define the dimensions for the variables
    stat = nf90_def_dim( ncid, "longitude", nlong, LongDimId )

    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining longitude: ", & 
        trim( nf90_strerror( stat ) )
      stop
    endif

    stat =  nf90_def_dim( ncid, "latitude", nlat, LatDimId )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining latitude: ", & 
        trim( nf90_strerror( stat ) )
      stop
    endif

    stat = nf90_def_dim( ncid, "altitude", iz, AltDimId )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining altitude: ", & 
      trim( nf90_strerror( stat ) )
      stop
    endif

    stat =  nf90_def_dim( ncid, "time", NF90_UNLIMITED, TimeDimId )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time", & 
        trim( nf90_strerror( stat ) )
      stop
    endif

    ! Define the initial variables for the dimensions
    stat = nf90_def_var( ncid, "longitude", NF90_FLOAT, & 
                         (/LongDimId/), LongVarId )

    stat = nf90_def_var( ncid, "latitude", NF90_FLOAT, & 
                         (/LatDimId/), LatVarId )

    stat = nf90_def_var( ncid, "altitude", NF90_FLOAT, & 
                        (/AltDimId/), AltVarId )

    ! grads2nc stores time as a double prec. value, so we follow that
    stat = nf90_def_var( ncid, "time", NF90_DOUBLE, & 
                         (/TimeDimId/), TimeVarId )

    ! Assign attribute values

    ! Time attribute
    stat = nf90_put_att( ncid, TimeVarId, "cartesian_axis", "T" )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time: ", trim( nf90_strerror( stat ) )
      stop
    endif

    call format_date( day, month, year, time, TimeUnits )

    stat = nf90_put_att( ncid, TimeVarId, "units", TimeUnits )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time: ", trim( nf90_strerror( stat ) )
      stop
    endif

    stat = nf90_put_att( ncid, TimeVarId, "ipositive", 1 )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time: ", trim( nf90_strerror( stat ) )
      stop
    endif

    stat = nf90_put_att( ncid, TimeVarId, "calendar_type", "Gregorian" )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time", trim( nf90_strerror( stat ) )
      stop
    endif

    ! Define Location
    ! X & Y coordinates
    stat = nf90_put_att( ncid, LongVarId, "cartesian_axis", "X" )

    stat = nf90_put_att( ncid, LongVarId, "units",  "degrees_E" )

    stat = nf90_put_att( ncid, LongVarId, "ipositive",  1 )

    stat = nf90_put_att( ncid, LatVarId, "cartesian_axis",  "Y" )

    stat = nf90_put_att( ncid, LatVarId, "units", "degrees_N" )

    stat = nf90_put_att( ncid, LatVarId, "ipositive", 1 )

    ! Altitude, Z coordinate
    stat = nf90_put_att( ncid, AltVarId, "cartesian_axis",  "Z" )

    stat = nf90_put_att( ncid, AltVarId, "units", "meters" )

    stat = nf90_put_att( ncid, AltVarId, "positive",  "up" )

    stat = nf90_put_att( ncid, AltVarId, "ipositive", 1 )

    return
  end subroutine define_netcdf

!-----------------------------------------------------------------------
  subroutine close_netcdf( ncf )

!      Description:
!      Close a previously opened stats file.

!      Notes:
!      I assume nf90_close() exists so that the NetCDF libraries can do a
!      form of buffered I/O, but I don't know the implementation
!      details. -dschanen
!-----------------------------------------------------------------------

    use output_file_module, only: & 
        outputfile ! Type
    use netcdf, only: & 
        NF90_NOERR,  & ! Variable
        nf90_close,  & ! Procedure(s)
        nf90_strerror

    use constants, only:  & 
        fstderr  ! Variable

    implicit none

! Input/Output Variables
    type (outputfile), intent(inout) :: ncf

! Local Variables
    integer :: stat

    stat = nf90_close( ncf%iounit )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error closing file "//  & 
        trim( ncf%fname )//": ", trim( nf90_strerror( stat ) )
      stop
    endif

    return
  end subroutine close_netcdf

!-----------------------------------------------------------------------
  subroutine first_write( ncf )

! Description:
!   Used on the first call to write_nc to finalize definitions
!   for the dataset, including the attributes for variable records.
! References:
!   None
!-----------------------------------------------------------------------

    use netcdf, only: & 
      NF90_NOERR,  & ! Constants
      NF90_FLOAT,  & 
      NF90_GLOBAL, &
      nf90_def_var,  & ! Procedure(s)
      nf90_strerror, & 
      nf90_put_att, & 
      nf90_enddef

    use output_file_module, only: &
      outputfile ! Derived type

    use constants, only:  &
      fstderr ! Variable

    use parameters_model, only: &
      T0, &       ! Real variables
      ts_nudge, &
      sclrtol    ! Real array variable

    use parameters_tunable, only: &
      params_list ! Variable names (characters)

    use parameters_tunable, only: &
      get_parameters ! Subroutine

    use parameter_indices, only: &
      nparams ! Integer

    use model_flags, only: &
      l_local_kk, & ! Logicals
      l_pos_def, &
      l_hole_fill, &
      l_clip_semi_implicit, &
      l_3pt_sqd_dfsn, &
      l_standard_term_ta, &
      l_single_C2_Skw, &
      l_gamma_Skw, &
      l_bugsrad, &
      l_uv_nudge, &
      l_tke_aniso

    use parameters_microphys, only: &
      micro_scheme, & ! Variable(s)
      l_cloud_sed

    implicit none

    ! Input/Output Variables
    type (outputfile), intent(inout) :: ncf

    ! Local Variables
    integer, dimension(:), allocatable :: stat

    real, dimension(nparams) :: params ! Tunable parameters

    integer :: i     ! Array index
    logical :: l_error ! Error stat

    character(len=10) :: current_time
    character(len=8)  :: current_date
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

! var_range = (/ -1.e31, 1.e31 /)

! Explanation:  The NetCDF documentation claims the NF90_UNLIMITED
!   variable should be the first dimension, but def_var is somehow
!   inverted and requires the opposite.  After writing, these
!   dimensions are all in the opposite order of this in the file.
!   -dschanen

    var_dim(1) = ncf%LongDimId
    var_dim(2) = ncf%LatDimId
    var_dim(3) = ncf%AltDimId
    var_dim(4) = ncf%TimeDimId ! The NF90_UNLIMITED dimension

    allocate( stat( ncf%nvar ) )

    l_error = .false.

    do i = 1, ncf%nvar, 1
!        stat(i) = nf90_def_var( ncf%iounit, trim( ncf%var(i)%name ),
!    .             NF90_FLOAT, (/ncf%TimeDimId, ncf%AltDimId,
!    .             ncf%LatDimId, ncf%LongDimId/), ncf%var(i)%Id )
      stat(i) = nf90_def_var( ncf%iounit, trim( ncf%var(i)%name ), & 
                NF90_FLOAT, var_dim(:), ncf%var(i)%Id )
      if ( stat(i) /= NF90_NOERR ) then
        write(fstderr,*) "Error defining variable ",  & 
          ncf%var(i)%name //": ", trim( nf90_strerror( stat(i) ) )
        l_error = .true.
      endif

      stat(i) = nf90_put_att( ncf%iounit, ncf%var(i)%Id, & 
                "valid_range", var_range(1:2) )
      if ( stat(i) /= NF90_NOERR ) then
        write(fstderr,*) "Error defining valid range", & 
          trim( nf90_strerror( stat(i) ) )
        l_error = .true.
      endif

      stat(i) = nf90_put_att( ncf%iounit, ncf%var(i)%Id, "title",  & 
                trim( ncf%var(i)%description ) )
      if ( stat(i) /= NF90_NOERR ) then
        write(fstderr,*) "Error in description", & 
          trim( nf90_strerror( stat(i) ) )
        l_error = .true.
      endif

      stat(i) = nf90_put_att( ncf%iounit, ncf%var(i)%Id, "units",  & 
                trim( ncf%var(i)%units ) )
      if ( stat(i) /= NF90_NOERR ) then
        write(fstderr,*) "Error in units", & 
          trim( nf90_strerror( stat(i) ) )
        l_error = .true.
      endif
    end do

    if ( l_error ) stop "Error in definition"

    deallocate( stat )

    allocate( stat(4) )

    ! Define global attributes of the file, for reproducing the results and
    ! determining how a run was configured
    stat(1) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "Conventions", "COARDS" )
    stat(2) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "model", "CLUBB" )

    ! Figure out when the model is producing this file
    call date_and_time( current_date, current_time )

    stat(3) = nf90_put_att( &
                         ncf%iounit, NF90_GLOBAL, "created_on", &
                         current_date(1:4)//'-'//current_date(5:6)//'-'// &
                         current_date(7:8)//' '// &
                         current_time(1:2)//':'//current_time(3:4) )

    stat(4) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "micro_scheme", &
                            trim( micro_scheme ) )

    if ( any( stat /= NF90_NOERR ) ) then
      write(fstderr,*) "Error writing model information"
      do i = 1, size( stat ), 1
        write(fstderr,*) trim( nf90_strerror( stat(i) ) )
      end do
      stop
    end if

    ! Write the model flags to the file
    deallocate( stat )
    allocate( stat(12) ) ! # of model flags

    stat(1) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_local_kk", lchar( l_local_kk ) )
    stat(2) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_pos_def", lchar( l_pos_def ) )
    stat(3) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_hole_fill", lchar( l_hole_fill ) )
    stat(4) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_clip_semi_implicit", &
      lchar( l_clip_semi_implicit ) )
    stat(5) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_3pt_sqd_dfsn", &
      lchar( l_3pt_sqd_dfsn ) )
    stat(6) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_standard_term_ta", &
      lchar( l_standard_term_ta ) )
    stat(7) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_single_C2_Skw", &
      lchar( l_single_C2_Skw ) )
    stat(8) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_gamma_Skw", lchar( l_gamma_Skw ) )
    stat(9) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_bugsrad", lchar( l_bugsrad ) )
    stat(10) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_cloud_sed", lchar( l_cloud_sed ) )
    stat(11) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_uv_nudge", lchar( l_uv_nudge ) )
    stat(12) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_tke_aniso", lchar( l_tke_aniso ) )

    if ( any( stat /= NF90_NOERR ) ) then
      write(fstderr,*) "Error writing model flags"
      do i = 1, size( stat ), 1
        write(fstderr,*) i, trim( nf90_strerror( stat(i) ) )
      end do
      stop
    end if

    ! Write model parameter values to the file
    deallocate( stat )
    allocate( stat(nparams) )

    stat(1) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "T0", T0 )
    stat(2) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "ts_nudge", ts_nudge )
    stat(3) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "sclrtol", sclrtol )

    call get_parameters( params )

    do i = 1, nparams, 1
      stat(i) = nf90_put_att( ncf%iounit, NF90_GLOBAL, params_list(i), params(i) )
    end do

    if ( any( stat /= NF90_NOERR ) ) then
      write(fstderr,*) "Error writing parameters"
      do i = 1, nparams, 1
        write(fstderr,*) i, trim( nf90_strerror( stat(i) ) )
      end do
      stop
    end if

    stat(1) = nf90_enddef( ncf%iounit ) ! end definitions
    if ( stat(1) /= NF90_NOERR ) then
      write(fstderr,*) "Error finalizing definitions", & 
        trim( nf90_strerror( stat(1) ) )
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

    use netcdf, only: & 
        NF90_NOERR,   & ! Variable(s)
        nf90_put_var,  & ! Procedure(s)
        nf90_strerror
    use output_file_module, only: & 
        outputfile ! Type
    use constants, only:  & 
        fstderr ! Variable

    implicit none

!      Constants
!      real, parameter, dimension(1) :: deg_east  = 0.0
!      real, parameter, dimension(1) :: deg_north = 0.0

!      Input
    type (outputfile), intent(inout) :: ncf

    integer :: stat

    stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%AltVarId,  & 
                         values=ncf%z(ncf%ia:ncf%iz) )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error entering grid: ",  & 
        trim( nf90_strerror( stat ) )
      stop
    endif

    stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%LongVarId,  & 
                         values=ncf%rlon )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error entering longitude: ",  & 
        trim( nf90_strerror( stat ) )
      stop
    endif

    stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%LatVarId,  & 
                         values=ncf%rlat )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error entering latitude: ",  & 
        trim( nf90_strerror( stat ) )
      stop
    endif

    return
  end subroutine write_grid

!-----------------------------------------------------------------------

  subroutine format_date & 
             ( day_in, month_in, year_in, time_in, date )

!      Description:
!      Put the model date in a format that udunits and NetCDF can easily
!      handle.  GrADSnc is dumb and apparently cannot handle time
!      intervals < 1 minute.

!      Notes:
!      Adapted from the original GrADS version written by Chris Golaz.
!      Uses Fortran `internal' files to write the string output.
!-----------------------------------------------------------------------

    use calendar, only: compute_current_date

    use stats_precision, only:  & 
        time_precision ! Variable(s)

    implicit none

! Input Variables
    integer, intent(in) ::  & 
      day_in,           & ! Day of Month at Model Start   [dd]
      month_in,         & ! Month of Year at Model Start  [mm]
      year_in             ! Year at Model Start         [yyyy]

    real(kind=time_precision), intent(in) :: time_in ! Start time [s]

! Output Variables
    character(len=35), intent(out) :: date

    integer::  & 
    iday, imonth, iyear  ! Integer for day, month and year.
    real(kind=time_precision) :: st_time ! Start time [s]

    call compute_current_date( day_in, month_in,  & 
                               year_in, & 
                               time_in, & 
                               iday, imonth, & 
                               iyear, & 
                               st_time )

!      date(1:14) = "minutes since "
    date = "minutes since YYYY-MM-DD HH:MM:00.0"
!      date = "seconds since YYYY-MM-DD HH:MM:00.0"
    write(date(15:18),'(i4.4)') iyear
!      write(date(19),'(a1)') '-'
    write(date(20:21),'(i2.2)') imonth
!      write(date(22),'(a1)') '-'
    write(date(23:24),'(i2.2)') iday
!      write(date(25),'(a1)') ' '
    write(date(26:27),'(i2.2)') floor(st_time/3600)
!      write(date(28),'(a1)') ":"
    write(date(29:30),'(i2.2)') int( mod(nint(st_time),3600) / 60 )
!      write(date(30:35),'(a5)') ":00.0"

    return
  end subroutine format_date

!===============================================================================
  character function lchar( l_input )
! Description:
!   Cast a logical to a character data type
! References:
!   None
!-------------------------------------------------------------------------------

    implicit none

    logical, intent(in) :: l_input

    if ( l_input ) then
      lchar = 'T'
    else
      lchar = 'F'
    end if

    return
  end function lchar

#endif /*NETCDF*/
end module output_netcdf
