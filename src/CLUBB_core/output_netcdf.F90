!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module output_netcdf

! Description:
!   Functions and subroutines for writing NetCDF files

! References:
!   <http://www.unidata.ucar.edu/software/netcdf/docs/>
!-------------------------------------------------------------------------------

  implicit none

#ifdef NETCDF

  public :: open_netcdf_for_writing, write_netcdf, &
            write_netcdf_w_diff_output_gr, close_netcdf, &
            format_date, first_write, output_multi_col_fields

  private :: define_netcdf, write_grid, write_netcdf_helper

  ! Constant parameters
  ! This will truncate all timesteps smaller than 1 mn to a minute for 
  ! the purposes of viewing the data in grads
  logical, parameter, private :: &
    l_grads_netcdf_boost_ts = .false. 

  private ! Default scope

  contains

!-------------------------------------------------------------------------------
  subroutine open_netcdf_for_writing( nlat, nlon, fdir, fname, ia, iz, zgrid,  & 
                                      day, month, year, lat_vals, lon_vals, & 
                                      time, dtwrite, nvar, &
                                      ncf, err_code, &
                                      nsamp ) ! optional

! Description:
!   Defines the structure used to reference the file `ncf'

! References:
!   None
!-------------------------------------------------------------------------------
    use netcdf, only: & 
        NF90_CLOBBER, & ! Variable(s)
        NF90_NOERR,   & 
        nf90_create,  & ! Procedure
        nf90_strerror

    use stat_file_module, only: & 
        stat_file ! Type

    use clubb_precision, only:  & 
        time_precision, & ! Variable(s)
        core_rknd

    use constants_clubb, only:  & 
        fstderr ! Variable(s)

    use error_code, only: &
        clubb_fatal_error     ! Constant

    implicit none

    ! Input Variables
    character(len=*), intent(in) ::  & 
      fdir,   & ! Directory name of file
      fname     ! File name

    integer, intent(in) ::  & 
      nlat, nlon,       & ! Number of points in the X and Y
      day, month, year, & ! Time
      ia, iz,           & ! First and last grid point
      nvar                ! Number of variables

    real( kind = core_rknd ), dimension(nlat), intent(in) ::  & 
      lat_vals ! Latitudes   [degrees_E]

    real( kind = core_rknd ), dimension(nlon), intent(in) ::  & 
      lon_vals ! Longitudes  [degrees_N]

    real( kind = core_rknd ), intent(in) :: & 
      dtwrite ! Time between write intervals   [s]

    real( kind = time_precision ), intent(in) ::  & 
     time   ! Current time                    [s]

    real( kind = core_rknd ), dimension(:), intent(in) ::  & 
      zgrid  ! The model grid                  [m]

    ! Input/output Variables
    type (stat_file), intent(inout) :: ncf

    integer, intent(inout) :: &
      err_code      ! Error code catching and relaying any errors occurring in this subroutine

    ! Number of SILHS samples, used only for SILHS sample outputting
    integer, optional, intent(in) :: nsamp

    ! Local Variables
    integer :: stat  ! Error status
    integer :: k     ! Array index

    ! ---- Begin Code ----

    ncf%nvar    = nvar

    ! If there is no data to write, then return
    if ( ncf%nvar == 0 ) then
      return
    end if

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
    ncf%nlat   = nlat
    ncf%nlon   = nlon
    ncf%time   = time

    ncf%dtwrite = dtwrite
    ncf%ntimes = 0

! According to Chris Vogl, netcdf can handle time steps < 1 min.  So this check is unneeded.
!    ! Check to make sure the timestep is appropriate. The GrADS program does not support an
!    ! output timestep less than 1 minute.  Other programs can read netCDF files like this
!    if ( dtwrite < sec_per_min ) then
!      write(fstderr,*) "Warning: GrADS program requires an output timestep of at least &
!                       &one minute, but the requested output timestep &
!                       &(stats_tout) is less than one minute."
!      if ( .not. l_allow_small_stats_tout ) then
!        write(fstderr,*) "To override this warning, set l_allow_small_stats_tout = &
!                         &.true. in the stats_setting namelist in the &
!                         &appropriate *_model.in file."
!        write(fstderr,*) "Fatal error in open_netcdf_for_writing"
!        err_code = clubb_fatal_error
!        return
!      end if
!    end if ! dtwrite < sec_per_min

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

    allocate( ncf%lat_vals(1:nlat), ncf%lon_vals(1:nlon) )

    ncf%lat_vals = lat_vals
    ncf%lon_vals = lon_vals

    ! If nsamp is present, SILHS samples are being handled.  Therefore set
    ! ncf%nsamp and ncf%samp_idx.  samp_idx holds the SILHS sample indices.
    if ( present(nsamp) ) then
      ncf%nsamp = nsamp
      allocate( ncf%samp_idx(1:nsamp) )
      forall( k=1:nsamp )
        ncf%samp_idx(k) = real( k, kind = core_rknd )
      end forall
    endif

    ! Create NetCDF dataset: enter define mode
    stat = nf90_create( path = trim( fdir )//trim( fname )//'.nc',  & 
                        cmode = NF90_CLOBBER,  & ! overwrite existing file
                        ncid = ncf%iounit )

    if ( stat /= NF90_NOERR ) then
      write(unit=fstderr,fmt=*) "Error opening file: ",  & 
        trim( fdir )//trim( fname )//'.nc', & 
        trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    call define_netcdf( ncf%iounit, ncf%nlat, ncf%nlon, ncf%iz, ncf%nsamp, & ! In
                  ncf%day, ncf%month, ncf%year, ncf%time, & ! In
                  ncf%SampDimId, ncf%LatDimId, ncf%LongDimId, ncf%AltDimId, ncf%TimeDimId, &! Out
                  ncf%SampVarId, ncf%LatVarId, ncf%LongVarId, ncf%AltVarId, ncf%TimeVarId, &! Out
                  err_code )                                                                ! Inout

    return
  end subroutine open_netcdf_for_writing

!-------------------------------------------------------------------------------

  subroutine write_netcdf_helper( l_different_output_grid, &
                                  gr_source, gr_target, &
                                  zm_zt, &
                                  total_idx_rho_lin_spline, &
                                  rho_lin_spline_vals, &
                                  rho_lin_spline_levels, &
                                  ncf, err_code )

! Description:
!   Writes some data to the NetCDF dataset, but doesn't close it.
!
! References:
!   None   
!-------------------------------------------------------------------------------

    use netcdf, only: & 
        NF90_NOERR,  & ! Variable(s)
        nf90_put_var,  & ! Procedure
        nf90_strerror

    use stat_file_module, only: & 
        stat_file ! Variable

    use constants_clubb, only:  & 
        fstderr, & ! Variable
        sec_per_min

    use error_code, only: &
        clubb_fatal_error     ! Constant

    use clubb_precision, only: &
        stat_rknd, &
        core_rknd, &
        time_precision ! Constant

    use grid_adaptation_module, only: &
        remap_zm_values, remap_zt_values

    use grid_class, only: grid ! Type

    implicit none

    logical, intent(in) :: &
      l_different_output_grid  ! use different grid to output values to file
    
    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of points in the rho spline

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, &  ! the rho values for constructing the spline for remapping;
                              ! only used if l_different_output_gr is .true.
      rho_lin_spline_levels   ! the levels at which the rho values are given;
                             ! only used if l_different_output_gr is .true.

    type( grid ), intent(in) :: &
      gr_source, & ! the grid where the values are currently given on;
                   ! only used if l_different_output_gr is .true.
      gr_target    ! the grid where the values should be remapped to;
                   ! only used if l_different_output_gr is .true.

    character(len=2), intent(in) :: zm_zt ! whether the file is the zm or zt file

    type (stat_file), intent(inout) :: ncf    ! The file

    integer, intent(inout) :: &
      err_code      ! Error code catching and relaying any errors occurring in this subroutine

    ! Local Variables
    integer, dimension(:), allocatable :: stat ! Error status

    real(kind=8), dimension(1) :: time         ! Time          [s]

    real(kind=stat_rknd), dimension(:,:,:,:), allocatable :: &
      grid_avg_var_diff_gr  ! only needed if l_different_output_grid = .true., to store the
                            ! remapped values on the output grid

    real(kind=stat_rknd), dimension(:,:,:,:,:), allocatable :: &
      samples_of_var_diff_gr  ! only needed if l_different_output_grid = .true., to store the
                              ! remapped values on the output grid

    integer :: i, &     ! Array index
               samp, lon, lat ! loop vars

    ! ---- Begin Code ----
    ! TODO remove if, for removing compiler warnings about accessing uninitialized objects
    if ( l_different_output_grid ) then
      ! TODO is ncf%iz always the count or do I have to take ncf%ia - ncf%iz ?
      allocate( grid_avg_var_diff_gr(1,ncf%nlon,ncf%nlat,ncf%iz) )
      allocate( samples_of_var_diff_gr(1,ncf%nsamp,ncf%nlon,ncf%nlat,ncf%iz) )
    endif

    ! If there is no data to write, then return
    if ( ncf%nvar == 0 ) then
      return
    end if

    ncf%ntimes = ncf%ntimes + 1

    allocate( stat( ncf%nvar ) )
    if ( l_grads_netcdf_boost_ts ) then
      time = real( nint( real(ncf%ntimes, kind=time_precision) &
                            * real(ncf%dtwrite / sec_per_min, time_precision) ), &
                              kind=time_precision ) ! minutes(rounded)
    else
      time = real( ncf%ntimes, kind=time_precision ) &
           * real( ncf%dtwrite, kind=time_precision )  ! seconds
    end if

    stat(1) = nf90_put_var( ncid=ncf%iounit, varid=ncf%TimeVarId,  & 
                            values=time(1), start=(/ncf%ntimes/) )
    if ( stat(1) /= NF90_NOERR ) then
      write(fstderr,*) "time variable nf90_put_var failed"
      err_code = clubb_fatal_error
      return
    end if

    ! If the grid_avg_var are allocated, then print to 4d netcdf.
    ! Otherwise, if the samples_of_var are allocated, print to 5d
    do i = 1, ncf%nvar, 1
      if ( allocated(ncf%grid_avg_var) ) then
        if ( l_different_output_grid ) then
          do lon = 1, ncf%nlon
            do lat = 1, ncf%nlat
              if ( zm_zt == 'zm' ) then
                grid_avg_var_diff_gr(:,lon,lat,:) &
                  = remap_zm_values( 1, gr_source, gr_target, &
                                     total_idx_rho_lin_spline, &
                                     rho_lin_spline_vals, &
                                     rho_lin_spline_levels, &
                                     ncf%grid_avg_var(i)%ptr(lon,lat,ncf%ia:gr_source%nzm ) &
                                    )
              else if ( zm_zt == 'zt' ) then
                grid_avg_var_diff_gr(:,lon,lat,:) &
                  = remap_zt_values( 1, gr_source, gr_target, &
                                     total_idx_rho_lin_spline, &
                                     rho_lin_spline_vals, &
                                     rho_lin_spline_levels, &
                                     ncf%grid_avg_var(i)%ptr(lon,lat,ncf%ia:gr_source%nzt ) &
                                    )
              else
                error stop 'Invalid value for zm_zt in write_netcdf_helper()'
              endif
            end do
          end do
          stat(i)  &
          = nf90_put_var( ncid=ncf%iounit, varid=ncf%grid_avg_var(i)%indx,  &
                          values=grid_avg_var_diff_gr(1,:,:,:),  &
                          start=(/1,1,1,ncf%ntimes/), &
                          count=(/ncf%nlon,ncf%nlat,ncf%iz,1/) )
        else
          stat(i)  &
          = nf90_put_var( ncid=ncf%iounit, varid=ncf%grid_avg_var(i)%indx,  &
                          values=ncf%grid_avg_var(i)%ptr(:,:,ncf%ia:ncf%iz),  &
                          start=(/1,1,1,ncf%ntimes/), &
                          count=(/ncf%nlon,ncf%nlat,ncf%iz,1/) )
        endif
      elseif ( allocated(ncf%samples_of_var) ) then
        if ( l_different_output_grid ) then
          do samp = 1, ncf%nsamp
            do lon = 1, ncf%nlon
              do lat = 1, ncf%nlat
                if ( zm_zt == 'zm' ) then
                  samples_of_var_diff_gr(:,samp,lon,lat,:) &
                    = remap_zm_values( 1, gr_source, gr_target, &
                                       total_idx_rho_lin_spline, &
                                       rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       ncf%samples_of_var(i)%ptr(samp,lon,lat,ncf%ia:gr_source%nzm)&
                                      )
                else if ( zm_zt == 'zt' ) then
                  samples_of_var_diff_gr(:,samp,lon,lat,:) &
                    = remap_zt_values( 1, gr_source, gr_target, &
                                       total_idx_rho_lin_spline, &
                                       rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       ncf%samples_of_var(i)%ptr(samp,lon,lat,ncf%ia:gr_source%nzt) &
                                      )
                else
                  error stop 'Invalid value for zm_zt in write_netcdf_helper()'
                endif
              end do
            end do
          end do
          stat(i)  &
          = nf90_put_var( ncid=ncf%iounit, varid=ncf%samples_of_var(i)%indx,  &
                          values=samples_of_var_diff_gr(1,:,:,:,:),  &
                          start=(/1,1,1,1,ncf%ntimes/), &
                          count=(/ncf%nsamp,ncf%nlon,ncf%nlat,ncf%iz,1/) )
        else
          stat(i)  &
          = nf90_put_var( ncid=ncf%iounit, varid=ncf%samples_of_var(i)%indx,  &
                          values=ncf%samples_of_var(i)%ptr(:,:,:,ncf%ia:ncf%iz),  &
                          start=(/1,1,1,1,ncf%ntimes/), &
                          count=(/ncf%nsamp,ncf%nlon,ncf%nlat,ncf%iz,1/) )
        endif
      endif
    enddo ! i=1..nvar

    if ( any (stat /= NF90_NOERR ) ) then
      do i=1,ncf%nvar,1
        if( stat(i) /= NF90_NOERR ) then
          if ( allocated(ncf%grid_avg_var) ) then
            write(unit=fstderr,fmt=*) ncf%grid_avg_var(i)%name,  &
              trim( nf90_strerror( stat(i) ) )
          elseif ( allocated(ncf%samples_of_var) ) then
            write(unit=fstderr,fmt=*) ncf%samples_of_var(i)%name,  &
              trim( nf90_strerror( stat(i) ) )
          endif
        end if
      end do
      write(fstderr,*) "nf90_put_var error"
      err_code = clubb_fatal_error
      return
    end if

    deallocate( stat )

    ! TODO remove if, for removing compiler warnings about accessing uninitialized objects
    if ( l_different_output_grid ) then
      deallocate( grid_avg_var_diff_gr )
      deallocate( samples_of_var_diff_gr )
    endif

    return
  end subroutine write_netcdf_helper

  !-------------------------------------------------------------------------------

  subroutine write_netcdf( ncf, err_code )

  ! Description:
  !   Writes some data to the NetCDF dataset, but doesn't close it.
  !
  ! References:
  !   None   
  !-------------------------------------------------------------------------------

    use stat_file_module, only: & 
        stat_file ! Variable

    use grid_class, only: grid ! Type

    use clubb_precision, only: core_rknd

    implicit none

    type (stat_file), intent(inout) :: ncf    ! The file

    integer, intent(inout) :: &
      err_code      ! Error code catching and relaying any errors occurring in this subroutine

    ! Local Variables
    logical :: &
      l_different_output_grid ! use different grid to write values to file

    integer :: &
      total_idx_rho_lin_spline_placeholder ! number of points in the rho spline

    real( kind = core_rknd ), dimension(:), allocatable :: &
      rho_lin_spline_vals_placeholder, &  ! the rho values for constructing the spline for
                                          ! remapping; only used if l_different_output_gr is .true.
      rho_lin_spline_levels_placeholder   ! the levels at which the rho values are given;
                                          ! only used if l_different_output_gr is .true.

    type( grid ) :: &
      gr_source_placeholder, & ! the grid where the values are currently given on;
                               ! only used if l_different_output_gr is .true.
      gr_target_placeholder    ! the grid where the values should be remapped to;
                               ! only used if l_different_output_gr is .true.

    character(len=2) :: zm_zt_placeholder ! whether the file is the zm or zt file
                                          ! only used if l_different_output_gr is .true.

    ! ---- Begin Code ----
    l_different_output_grid = .false.
    zm_zt_placeholder = ''
    total_idx_rho_lin_spline_placeholder = 1

    allocate( rho_lin_spline_vals_placeholder(total_idx_rho_lin_spline_placeholder) )
    allocate( rho_lin_spline_levels_placeholder(total_idx_rho_lin_spline_placeholder) )

    call write_netcdf_helper( l_different_output_grid, &
                              gr_source_placeholder, gr_target_placeholder, &
                              zm_zt_placeholder, &
                              total_idx_rho_lin_spline_placeholder, &
                              rho_lin_spline_vals_placeholder, &
                              rho_lin_spline_levels_placeholder, &
                              ncf, err_code )
     
    deallocate( rho_lin_spline_vals_placeholder )
    deallocate( rho_lin_spline_levels_placeholder )

  end subroutine write_netcdf

  !-------------------------------------------------------------------------------

  subroutine write_netcdf_w_diff_output_gr( gr_source, gr_target, &
                                            zm_zt, &
                                            total_idx_rho_lin_spline, &
                                            rho_lin_spline_vals, &
                                            rho_lin_spline_levels, &
                                            ncf, err_code )

  ! Description:
  !   Writes some data to the NetCDF dataset, but doesn't close it.
  !   The data is remapped to the given output grid before it is written to file.
  !
  ! References:
  !   None   
  !-------------------------------------------------------------------------------

    use stat_file_module, only: & 
        stat_file ! Variable

    use grid_class, only: grid ! Type

    use clubb_precision, only: core_rknd

    implicit none

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of points in the rho spline

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, &  ! the rho values for constructing the spline for remapping;
                              ! only used if l_different_output_gr is .true.
      rho_lin_spline_levels   ! the levels at which the rho values are given;
                             ! only used if l_different_output_gr is .true.

    type( grid ), intent(in) :: &
      gr_source, & ! the grid where the values are currently given on;
                   ! only used if l_different_output_gr is .true.
      gr_target    ! the grid where the values should be remapped to;
                   ! only used if l_different_output_gr is .true.

    character(len=2), intent(in) :: zm_zt ! whether the file is the zm or zt file

    type (stat_file), intent(inout) :: ncf    ! The file

    integer, intent(inout) :: &
      err_code      ! Error code catching and relaying any errors occurring in this subroutine

    ! Local variables
    logical :: &
      l_different_output_grid ! use different grid to write values to file

    ! ---- Begin Code ----
    l_different_output_grid = .true.

    call write_netcdf_helper( l_different_output_grid, &
                              gr_source, gr_target, &
                              zm_zt, &
                              total_idx_rho_lin_spline, &
                              rho_lin_spline_vals, &
                              rho_lin_spline_levels, &
                              ncf, err_code )

  end subroutine write_netcdf_w_diff_output_gr

!-------------------------------------------------------------------------------
  subroutine define_netcdf( ncid, nlat, nlon, iz, nsamp, &
                            day, month, year, time, & 
                            SampDimId, LatDimId, LongDimId, AltDimId, TimeDimId, &
                            SampVarId, LatVarId, LongVarId, AltVarId, TimeVarId, &
                            err_code )

! Description:
!   Used internally to create a definition for the NetCDF dataset
!
! References:
!   None
!-------------------------------------------------------------------------------
    use netcdf, only: & 
        NF90_NOERR,   & ! Constants
        NF90_DOUBLE, & 
        NF90_UNLIMITED

    use netcdf, only: & 
        nf90_def_dim,  & ! Functions
        nf90_strerror, & 
        nf90_def_var, & 
        nf90_put_att

    use clubb_precision, only:  & 
        time_precision ! Variable(s)

    use constants_clubb, only:  & 
        fstderr ! Variable(s)

    use error_code, only: &
        clubb_fatal_error     ! Constant

    implicit none

    integer, intent(in) ::  & 
      nlat,   & ! Number of points in the N/S direction
      nlon,   & ! Number of points in the E/W direction
      nsamp     ! Number of SILHS samples

    ! Input Variables
    integer, intent(in) ::  & 
      day, month, year,  & ! Time of year
      ncid,              & ! Number used by NetCDF for ref. the file
      iz                   ! Dimension in z

    real(kind=time_precision), intent(in) ::  & 
      time    ! Current model time [s]

    ! Output Variables
    integer, intent(out) ::  &
    ! NetCDF id's for dimensions, including for SILHS samples if needed
      SampDimId, LatDimId, LongDimId, AltDimId, TimeDimId

    ! NetCDF id's for data (e.g. longitude) associated with each dimension,
    ! including for SILHS samples if needed
    integer, intent(out) ::  & 
      SampVarId, LatVarId, LongVarId, AltVarId, TimeVarId

    ! Input/Output Variables
    integer, intent(inout) :: &
      err_code      ! Error code catching and relaying any errors occurring in this subroutine

    ! Local variables
    integer :: stat
    character(len=35) :: TimeUnits

    ! ---- Begin Code ----

    ! Define the dimensions for the variables.
    ! Start with SILHS samples so this dimension is listed first in the netCDF
    ! file.  Since ncf not present to test allocation, test using nsamp. nsamp
    ! is initialized to zero, will only be nonzero if printing SILHS samples.
    if ( nsamp > 0 ) then
      ! Define SILHS sample dimension
      stat =  nf90_def_dim( ncid, "lh_sample_number", nsamp, SampDimId )
      if ( stat /= NF90_NOERR ) then
        write(fstderr,*) "Error defining lh_sample_number: ", &
          trim( nf90_strerror( stat ) )
        err_code = clubb_fatal_error
        return
      end if
      ! Define SILHS sample number variable
      stat = nf90_def_var( ncid, "lh_sample_number", NF90_DOUBLE, &
                          (/SampDimId/), SampVarId )
      ! Attributes for SILHS sample number variable
      stat = nf90_put_att( ncid, SampVarId, "description", "SILHS sample (i.e. subcolumn) index" )
      stat = nf90_put_att( ncid, SampVarId, "units", "number" )
    endif !if nsamp>0

    stat = nf90_def_dim( ncid, "longitude", nlon, LongDimId )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining longitude: ", & 
        trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    stat =  nf90_def_dim( ncid, "latitude", nlat, LatDimId )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining latitude: ", & 
        trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    stat = nf90_def_dim( ncid, "altitude", iz, AltDimId )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining altitude: ", & 
      trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    stat =  nf90_def_dim( ncid, "time", NF90_UNLIMITED, TimeDimId )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time: ", & 
        trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    ! Define the initial variables for the dimensions
    ! Longitude = deg_E = X
    stat = nf90_def_var( ncid, "longitude", NF90_DOUBLE, & 
                         (/LongDimId/), LongVarId )

    ! Latitude = deg_N = Y
    stat = nf90_def_var( ncid, "latitude", NF90_DOUBLE, & 
                         (/LatDimId/), LatVarId )

    ! Altitude = meters above the surface = Z
    stat = nf90_def_var( ncid, "altitude", NF90_DOUBLE, & 
                        (/AltDimId/), AltVarId )

    ! grads2nc stores time as a double prec. value, so we follow that
    stat = nf90_def_var( ncid, "time", NF90_DOUBLE, & 
                         (/TimeDimId/), TimeVarId )

    ! Assign attribute values

    ! Time attribute
    stat = nf90_put_att( ncid, TimeVarId, "cartesian_axis", "T" )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time: ", trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    call format_date( day, month, year, time, & ! intent(in)
                      TimeUnits ) ! intent(out)

    stat = nf90_put_att( ncid, TimeVarId, "units", TimeUnits )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time: ", trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    stat = nf90_put_att( ncid, TimeVarId, "ipositive", 1 )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time: ", trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    stat = nf90_put_att( ncid, TimeVarId, "calendar_type", "Gregorian" )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time", trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

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

!-------------------------------------------------------------------------------
  subroutine close_netcdf( ncf )

! Description:
!   Close a previously opened stats file.

! Notes:
!   I assume nf90_close() exists so that the NetCDF libraries can do a
!   form of buffered I/O, but I don't know the implementation
!   details. -dschanen
!-------------------------------------------------------------------------------

    use stat_file_module, only: & 
        stat_file ! Type

    use netcdf, only: & 
        NF90_NOERR,  & ! Variable
        nf90_close,  & ! Procedure(s)
        nf90_strerror

    use constants_clubb, only:  & 
        fstderr  ! Variable

    implicit none

    ! Input/Output Variables
    type (stat_file), intent(inout) :: ncf

    ! Local Variables
    integer :: stat

    ! ---- Begin Code ----

    ! If there is no data to write, then return
    if ( ncf%nvar == 0 ) then
      return
    end if

    stat = nf90_close( ncf%iounit )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error closing file "//  & 
        trim( ncf%fname )//": ", trim( nf90_strerror( stat ) )
      error stop "Fatal error"
    end if

    return
  end subroutine close_netcdf

!-------------------------------------------------------------------------------
  subroutine first_write( clubb_params, sclr_dim, sclr_tol, &
                          l_uv_nudge, &
                          l_tke_aniso, &
                          l_standard_term_ta, &
                          ncf, err_code )

! Description:
!   Used on the first call to write_nc to finalize definitions
!   for the dataset, including the attributes for variable records.
! References:
!   None
!-------------------------------------------------------------------------------

    use netcdf, only: & 
        NF90_NOERR,  & ! Constants
        NF90_FLOAT,  &
        NF90_DOUBLE, & 
        NF90_GLOBAL, &
        nf90_def_var,  & ! Procedure(s)
        nf90_strerror, & 
        nf90_put_att, & 
        nf90_enddef

    use stat_file_module, only: &
        stat_file ! Derived type

    use constants_clubb, only:  &
        fstderr ! Variable

    use parameters_model, only: &
        T0, &       ! Real variables
        ts_nudge

    use parameters_tunable, only: &
        params_list ! Variable names (characters)

    use parameter_indices, only: &
        nparams ! Integer

    use model_flags, only: &
        l_pos_def, &
        l_hole_fill, &
        l_gamma_Skw

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_fatal_error     ! Constant

    implicit none

    ! External
    intrinsic :: date_and_time, huge, selected_real_kind, size, any, trim

    ! Enabling l_output_file_run_date allows the date and time that the netCDF
    ! output file is created to be included in the netCDF output file.
    ! Disabling l_output_file_run_date means that this information will not be
    ! included in the netCDF output file.  The advantage of disabling this
    ! output is that it allows for a check for binary differences between two
    ! netCDF output files.
    logical, parameter :: &
      l_output_file_run_date = .false.

    ! Input Variables
    integer, intent(in) :: &
      sclr_dim

    real( kind = core_rknd ), dimension(sclr_dim), intent(in) :: &
      sclr_tol

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    logical, intent(in) :: &
      l_uv_nudge,         & ! For wind speed nudging
      l_tke_aniso,        & ! For anisotropic turbulent kinetic energy, i.e. TKE = 1/2
                            ! (u'^2 + v'^2 + w'^2)
      l_standard_term_ta    ! Use the standard discretization for the turbulent advection terms.
                            ! Setting to .false. means that a_1 and a_3 are pulled outside of the
                            ! derivative in advance_wp2_wp3_module.F90 and in
                            ! advance_xp2_xpyp_module.F90.

    ! Input/Output Variables
    type (stat_file), intent(inout) :: ncf

    integer, intent(inout) :: &
      err_code      ! Error code catching and relaying any errors occurring in this subroutine

    ! Local Variables
    integer, dimension(:), allocatable :: stat

    integer :: netcdf_precision ! Level of precision for netCDF output

    integer :: i     ! Array index

    character(len=10) :: current_time
    character(len=8)  :: current_date
    ! Range for NetCDF variables
    real( kind = core_rknd ), dimension(2) :: var_range

    ! Dimensions for variables
    integer, allocatable, dimension(:) :: var_dim

!-------------------------------------------------------------------------------
!      Typical valid ranges (IEEE 754)

!      real(kind=4): +/- 3.4028235E+38
!      real(kind=8): +/- 1.797693134862316E+308
!      real(kind=16):+/- 1.189731495357231765085759326628007E+4932

!-------------------------------------------------------------------------------

    ! ---- Begin Code ----

    ! If there is no data to write, then return
    if ( ncf%nvar == 0 ) then
      return
    end if

    var_range(1) = -huge( var_range(1) )
    var_range(2) =  huge( var_range(2) )

! var_range = (/ -1.e31, 1.e31 /)

! Explanation:  The NetCDF documentation claims the NF90_UNLIMITED
!   variable should be the first dimension, but def_var is somehow
!   inverted and requires the opposite.  After writing, these
!   dimensions are all in the opposite order of this in the file.
!   -dschanen

    ! If samples_of_var is allocated, print to 5d netcdf, otherwise 4d.
    if ( allocated(ncf%samples_of_var) ) then
      allocate( var_dim(1:5) )
      var_dim(1) = ncf%SampDimId
      i = 1
    else
      allocate( var_dim(1:4) )
      i = 0
    endif

    var_dim(i+1) = ncf%LongDimId ! X
    var_dim(i+2) = ncf%LatDimId  ! Y
    var_dim(i+3) = ncf%AltDimId  ! Z
    var_dim(i+4) = ncf%TimeDimId ! The NF90_UNLIMITED dimension

    allocate( stat( ncf%nvar ) )

    select case (core_rknd)
      case ( selected_real_kind( p=5 ) )
        netcdf_precision = NF90_FLOAT
      case ( selected_real_kind( p=12 ) )
        netcdf_precision = NF90_DOUBLE
      case default
        netcdf_precision = NF90_DOUBLE
    end select

    ! Specify whether "grid_avg_var" or "samples_of_var"
    do i = 1, ncf%nvar, 1
      if ( allocated(ncf%grid_avg_var) ) then
        stat(i) = nf90_def_var( ncf%iounit, trim( ncf%grid_avg_var(i)%name ), &
                    netcdf_precision, var_dim(:), ncf%grid_avg_var(i)%indx )
      elseif ( allocated(ncf%samples_of_var) ) then
        stat(i) = nf90_def_var( ncf%iounit, trim( ncf%samples_of_var(i)%name ), &
                    netcdf_precision, var_dim(:), ncf%samples_of_var(i)%indx )
      endif
      if ( stat(i) /= NF90_NOERR ) then
        write(fstderr,*) "Error defining variable ",  & 
          ncf%grid_avg_var(i)%name //": ", trim( nf90_strerror( stat(i) ) )
        err_code = clubb_fatal_error
        return
      end if

      if ( allocated(ncf%grid_avg_var) ) then
        stat(i) = nf90_put_att( ncf%iounit, ncf%grid_avg_var(i)%indx, &
                    "valid_range", var_range(1:2) )
      elseif ( allocated(ncf%samples_of_var) ) then
        stat(i) = nf90_put_att( ncf%iounit, ncf%samples_of_var(i)%indx, &
                    "valid_range", var_range(1:2) )
      endif
      if ( stat(i) /= NF90_NOERR ) then
        write(fstderr,*) "Error defining valid range", & 
          trim( nf90_strerror( stat(i) ) )
        err_code = clubb_fatal_error
        return
      end if

      if ( allocated(ncf%grid_avg_var) ) then
        stat(i) = nf90_put_att( ncf%iounit, ncf%grid_avg_var(i)%indx, "long_name",  &
                  trim( ncf%grid_avg_var(i)%description ) )
      elseif ( allocated(ncf%samples_of_var) ) then
        stat(i) = nf90_put_att( ncf%iounit, ncf%samples_of_var(i)%indx, "long_name",  &
                  trim( ncf%samples_of_var(i)%description ) )
      endif
      if ( stat(i) /= NF90_NOERR ) then
        write(fstderr,*) "Error in description", & 
          trim( nf90_strerror( stat(i) ) )
        err_code = clubb_fatal_error
        return
      end if

      if ( allocated(ncf%grid_avg_var) ) then
        stat(i) = nf90_put_att( ncf%iounit, ncf%grid_avg_var(i)%indx, "units",  &
                  trim( ncf%grid_avg_var(i)%units ) )
      elseif ( allocated(ncf%samples_of_var) ) then
        stat(i) = nf90_put_att( ncf%iounit, ncf%samples_of_var(i)%indx, "units",  &
                  trim( ncf%samples_of_var(i)%units ) )
      endif
      if ( stat(i) /= NF90_NOERR ) then
        write(fstderr,*) "Error in units", & 
          trim( nf90_strerror( stat(i) ) )
        err_code = clubb_fatal_error
        return
      end if
    end do

    deallocate( stat )

    if ( l_output_file_run_date ) then
      allocate( stat(3) )
    else
      allocate( stat(2) )
    end if

    ! Define global attributes of the file, for reproducing the results and
    ! determining how a run was configured
    stat(1) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "Conventions", "COARDS" )
    stat(2) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "model", "CLUBB" )

    if ( l_output_file_run_date ) then

      ! Enabling l_output_file_run_date allows the date and time that the
      ! netCDF output file is created to be included in the netCDF output file.
      ! Disabling l_output_file_run_date means that this information will not
      ! be included in the netCDF output file.  The advantage of disabling this
      ! output is that it allows for a check for binary differences between two
      ! netCDF output files.

      ! Figure out when the model is producing this file
      call date_and_time( current_date, current_time )

      stat(3) = nf90_put_att(ncf%iounit, NF90_GLOBAL, "created_on", &
                             current_date(1:4)//'-'//current_date(5:6)//'-'// &
                             current_date(7:8)//' '// &
                             current_time(1:2)//':'//current_time(3:4) )

    end if ! l_output_file_run_date

    if ( any( stat /= NF90_NOERR ) ) then
      write(fstderr,*) "Error writing model information"
      do i = 1, size( stat ), 1
        write(fstderr,*) trim( nf90_strerror( stat(i) ) )
      end do
      err_code = clubb_fatal_error
      return
    end if

    ! Write the model flags to the file
    deallocate( stat )
    allocate( stat(6) ) ! # of model flags

    stat(1) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_pos_def", lchar( l_pos_def ) )
    stat(2) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_hole_fill", lchar( l_hole_fill ) )
    stat(3) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_standard_term_ta", &
      lchar( l_standard_term_ta ) )
    stat(4) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_gamma_Skw", lchar( l_gamma_Skw ) )
    stat(5) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_uv_nudge", lchar( l_uv_nudge ) )
    stat(6) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_tke_aniso", lchar( l_tke_aniso ) )

    if ( any( stat /= NF90_NOERR ) ) then
      write(fstderr,*) "Error writing model flags"
      do i = 1, size( stat ), 1
        write(fstderr,*) i, trim( nf90_strerror( stat(i) ) )
      end do
      err_code = clubb_fatal_error
      return
    end if

    ! Write model parameter values to the file
    deallocate( stat )
    allocate( stat(nparams) )

    stat(1) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "T0", T0 )
    stat(2) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "ts_nudge", ts_nudge )
    stat(3) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "sclr_tol", sclr_tol )

    do i = 1, nparams, 1
      stat(i) = nf90_put_att( ncf%iounit, NF90_GLOBAL, params_list(i), clubb_params(i) )
    end do

    if ( any( stat /= NF90_NOERR ) ) then
      write(fstderr,*) "Error writing parameters"
      do i = 1, nparams, 1
        write(fstderr,*) i, trim( nf90_strerror( stat(i) ) )
      end do
      err_code = clubb_fatal_error
      return
    end if

    stat(1) = nf90_enddef( ncf%iounit ) ! end definitions
    if ( stat(1) /= NF90_NOERR ) then
      write(fstderr,*) "Error finalizing definitions", & 
        trim( nf90_strerror( stat(1) ) )
      err_code = clubb_fatal_error
      return
    end if

    deallocate( stat )

    call write_grid( ncf, err_code )  ! define lat., long., and grid intent(inout)
    ncf%l_defined = .true.

    return
  end subroutine first_write

!-------------------------------------------------------------------------------
  subroutine write_grid( ncf, err_code )

! Description:
!   Writes inforation about latitude, longitude and the grid
! References:
!   None
!-------------------------------------------------------------------------------

    use netcdf, only: & 
        NF90_NOERR,   & ! Variable(s)
        nf90_put_var,  & ! Procedure(s)
        nf90_strerror

    use stat_file_module, only: & 
        stat_file ! Type

    use constants_clubb, only:  & 
        fstderr ! Variable

    use error_code, only: &
        clubb_fatal_error   ! Constant

    implicit none

    ! Input Variable(s)
    type (stat_file), intent(inout) :: ncf

    integer, intent(inout) :: &
      err_code      ! Error code catching and relaying any errors occurring in this subroutine

    integer :: stat

    ! ---- Begin Code ----

    stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%AltVarId, &
                         values=ncf%z(ncf%ia:ncf%iz) )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error entering grid: ", &
        trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%LongVarId, &
                         values=ncf%lon_vals )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error entering longitude: ", &
        trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%LatVarId, &
                         values=ncf%lat_vals )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error entering latitude: ", &
        trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    ! Write the SILHS sample indices if samples_of_var allocated
    if ( allocated(ncf%samples_of_var) ) then
      stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%SampVarId, &
                           values=ncf%samp_idx )
      if ( stat /= NF90_NOERR ) then
        write(fstderr,*) "Error entering grid: ", &
          trim( nf90_strerror( stat ) )
        err_code = clubb_fatal_error
        return
      end if
    endif

    return
  end subroutine write_grid

!-------------------------------------------------------------------------------

  subroutine format_date &
             ( day_in, month_in, year_in, time_in, &
               date )

! Description:
!   Put the model date in a format that udunits and NetCDF can easily
!   handle.  GrADSnc is dumb and apparently cannot handle time
!   intervals < 1 minute.

! Notes:
!   Adapted from the original GrADS version written by Chris Golaz.
!   Uses Fortran `internal' files to write the string output.
!-------------------------------------------------------------------------------

    use calendar, only:  &
        compute_current_date ! Procedure(s)

    use clubb_precision, only:  & 
        time_precision ! Variable(s)

    implicit none

    ! External
    intrinsic :: floor, int, mod, nint

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

    call compute_current_date( day_in, month_in,  & ! intent(in)
                               year_in, &  ! intent(in)
                               time_in, & ! intent(in)
                               iday, imonth, & ! intent(out)
                               iyear, &  ! intent(out)
                               st_time ) ! intent(out)

    if ( .not. l_grads_netcdf_boost_ts ) then
      date = "seconds since YYYY-MM-DD HH:MM:00.0"
    else
      date = "minutes since YYYY-MM-DD HH:MM:00.0"
    end if
    write(date(15:18),'(i4.4)') iyear
    write(date(20:21),'(i2.2)') imonth
    write(date(23:24),'(i2.2)') iday
    write(date(26:27),'(i2.2)') floor( st_time / 3600._time_precision )
    write(date(29:30),'(i2.2)') int( mod( nint( st_time ),3600 ) / 60 )
    
    if ( .not. l_grads_netcdf_boost_ts ) then
      write(date(32:33),'(i2.2)') nint(((real(mod( nint( st_time ),3600),kind=time_precision) / &
                     60._time_precision) - (real(int(mod( nint( st_time ),3600 ) / 60 ), & 
                                               kind=time_precision) ) )*60._time_precision)
    end if

    return
  end subroutine format_date

!===============================================================================
  character function lchar( l_input )
! Description:
!   Cast a logical to a character data type.
!
! References:
!   None
!-------------------------------------------------------------------------------

    implicit none

    ! Input Variable
    logical, intent(in) :: l_input

    ! ---- Begin Code ----

    if ( l_input ) then
      lchar = 'T'
    else
      lchar = 'F'
    end if

    return
  end function lchar

  subroutine output_multi_col_fields( nzm, nzt, ngrdcol, sclr_dim, edsclr_dim, &
                                      calls_per_out, l_output_double_prec, l_last_timestep, &
                                      gr, dt, output_file_prefix, &
                                      day, month, year, time_initial, &
                                      um, vm, up3, vp3, rtm, thlm, rtp3, thlp3, wp3, upwp, vpwp, &
                                      up2, vp2, wprtp, wpthlp, rtp2, thlp2, rtpthlp, wp2, &
                                      sclrm, sclrp3, wpsclrp, sclrp2, sclrprtp, sclrpthlp, &
                                      p_in_Pa, exner, rcm, cloud_frac, wp2thvp, wpthvp, rtpthvp, &
                                      thlpthvp, sclrpthvp, wp2rtp, wp2thlp, wpup2, wpvp2, &
                                      ice_supersat_frac, uprcp, vprcp, rc_coef_zm, wp4, wp2up2, &
                                      wp2vp2, um_pert, vm_pert, upwp_pert, vpwp_pert, edsclrm, &
                                      rcm_in_layer, cloud_cover, w_up_in_cloud, w_down_in_cloud, &
                                      cloudy_updraft_frac, cloudy_downdraft_frac, wprcp, &
                                      invrs_tau_zm, Kh_zt, Kh_zm, thlprcp )
    !
    ! Description:
    !   This subroutine outputs netcdf files with multiple columns.
    !   The dimensions of the variables in the files are:
    !     (ngrdcol,nzm,time) for ${output_file_prefix}_multi_col_zm.nc
    !     (ngrdcol,nzt,time) for ${output_file_prefix}_multi_col_zt.nc
    !
    ! NOTE:
    !   Currently most fields are not output, but passed in anyway in case
    !   they need to be output in the future.
    !
    ! References:
    !   Issue #1033 discusses a precessor procedure which duplicated then 
    !   output a _multi_col.nc netcdf file similar to the ones this now produces.
    !   https://github.com/larson-group/clubb/issues/1033
    !----------------------------------------------------------------------------
    
    use netcdf, only: &
      nf90_create, &
      nf90_def_dim, &
      nf90_def_var, &
      nf90_put_att, &
      nf90_put_var, &
      nf90_close, &
      nf90_open, &
      NF90_CLOBBER, &
      NF90_UNLIMITED, &
      NF90_DOUBLE, &
      NF90_FLOAT, &
      NF90_INT, &
      NF90_WRITE, &
      nf90_noerr

    use grid_class, only: &
      grid

    use clubb_precision, only: &
      core_rknd, &
      time_precision, &
      sp

    use constants_clubb, only: &
      one, &
      zero
                                    
    implicit none     

    type( grid ), intent(in) :: &
      gr

    integer, intent(in) :: &
      ngrdcol, &
      nzt, &
      nzm, &
      sclr_dim, &
      edsclr_dim, &
      calls_per_out

    logical, intent(in) :: &
      l_output_double_prec, &
      l_last_timestep

    real( kind = core_rknd ), intent(in) :: &
      dt

    character(len=100) :: &
      output_file_prefix

    integer, intent(in) ::  & 
      day, month, year ! Used to define start date the of simulation

    real(kind = time_precision ), intent(in) :: & 
      time_initial  ! Time of start of simulation     [s]

    ! These are prognostic or are planned to be in the future
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) ::  &
      um,                     & ! eastward grid-mean wind component (thermodynamic levels)   [m/s]
      vm,                     & ! northward grid-mean wind component (thermodynamic levels)   [m/s]
      up3,                    & ! u'^3 (thermodynamic levels)                    [m^3/s^3]
      vp3,                    & ! v'^3 (thermodynamic levels)                    [m^3/s^3]
      rtm,                    & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      thlm,                   & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      rtp3,                   & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
      thlp3,                  & ! th_l'^3 (thermodynamic levels)                 [K^3]
      wp3,                    & ! w'^3 (thermodynamic levels)                    [m^3/s^3]
      p_in_Pa,                & ! Air pressure (thermodynamic levels)       [Pa]
      exner,                  & ! Exner function (thermodynamic levels)     [-]
      rcm,                    & ! cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
      cloud_frac,             & ! cloud fraction (thermodynamic levels)          [-]
      wp2thvp,                & ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]
      wp2rtp,                 & ! w'^2 rt' (thermodynamic levels)      [m^2/s^2 kg/kg]
      wp2thlp,                & ! w'^2 thl' (thermodynamic levels)     [m^2/s^2 K]
      wpup2,                  & ! w'u'^2 (thermodynamic levels)        [m^3/s^3]
      wpvp2,                  & ! w'v'^2 (thermodynamic levels)        [m^3/s^3]
      ice_supersat_frac,      & ! ice cloud fraction (thermo. levels)  [-]
      um_pert,                & ! perturbed <u>       [m/s]
      vm_pert,                & ! perturbed <v>       [m/s]
      rcm_in_layer,           & ! rcm within cloud layer                          [kg/kg]
      cloud_cover,            & ! cloud cover                                     [-]
      w_up_in_cloud,          & ! Average cloudy updraft velocity       [m/s]
      w_down_in_cloud,        & ! Average cloudy downdraft velocity     [m/s]
      cloudy_updraft_frac,    & ! cloudy updraft fraction               [-]
      cloudy_downdraft_frac,  & ! cloudy downdraft fraction             [-]
      Kh_zt                     ! Eddy diffusivity coefficient on thermodynamic levels   [m^2/s]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) ::  &
      upwp,                   & ! u'w' (momentum levels)                         [m^2/s^2]
      vpwp,                   & ! v'w' (momentum levels)                         [m^2/s^2]
      up2,                    & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2,                    & ! v'^2 (momentum levels)                         [m^2/s^2]
      wprtp,                  & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
      wpthlp,                 & ! w'th_l' (momentum levels)                      [(m/s) K]
      rtp2,                   & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      thlp2,                  & ! th_l'^2 (momentum levels)                      [K^2]
      rtpthlp,                & ! r_t'th_l' (momentum levels)                    [(kg/kg) K]
      wp2,                    &  ! w'^2 (momentum levels)                         [m^2/s^2]
      wpthvp,                 & ! < w' th_v' > (momentum levels)                 [kg/kg K]
      rtpthvp,                & ! < r_t' th_v' > (momentum levels)               [kg/kg K]
      thlpthvp,               & ! < th_l' th_v' > (momentum levels)              [K^2]
      uprcp,                  & ! < u' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      vprcp,                  & ! < v' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      rc_coef_zm,             & ! Coef of X'r_c' in Eq. (34) (m-levs.) [K/(kg/kg)]
      wp4,                    & ! w'^4 (momentum levels)               [m^4/s^4]
      wp2up2,                 & ! w'^2 u'^2 (momentum levels)          [m^4/s^4]
      wp2vp2,                 & ! w'^2 v'^2 (momentum levels)          [m^4/s^4]
      upwp_pert,              & ! perturbed <u'w'>    [m^2/s^2]
      vpwp_pert,              & ! perturbed <v'w'>    [m^2/s^2] 
      wprcp,                  & ! w'r_c' (momentum levels)              [(kg/kg) m/s]
      invrs_tau_zm,           & ! One divided by tau on zm levels       [1/s]
      Kh_zm,                  & ! Eddy diffusivity coefficient on momentum levels        [m^2/s]
      thlprcp                   ! thl'rc'              [K kg/kg]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt,sclr_dim) :: &
      sclrm,     & ! Passive scalar mean (thermo. levels) [units vary]
      sclrp3       ! sclr'^3 (thermodynamic levels)       [{units vary}^3]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm,sclr_dim) :: &
      wpsclrp,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
      sclrp2,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
      sclrprtp,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp, & ! sclr'thl' (momentum levels)          [{units vary} K]
      sclrpthvp    ! < sclr' th_v' > (momentum levels)   [units vary]

    ! Eddy passive scalar variable
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt,edsclr_dim) :: &
      edsclrm   ! Eddy passive scalar grid-mean (thermo. levels)   [units vary]

    !------------------------ Local Variables ------------------------

    ! Index of zm variables
    integer, dimension(15), save :: &
      ind_zm

    ! Index of zt variables
    integer, dimension(15), save :: &
      ind_zt

    ! Variables we need to save for the netcdf writing
    integer, save :: &
      ncid_zm, &
      column_id_zm,     vertical_id_zm,     time_id_zm, &
      column_var_id_zm, vertical_var_id_zm, time_var_id_zm, &
      ncid_zt, &
      column_id_zt,     vertical_id_zt,     time_id_zt, &
      column_var_id_zt, vertical_var_id_zt, time_var_id_zt
    
    character(len=100), save :: &
      multicol_nc_file_zm, &
      multicol_nc_file_zt

    character(len=35) :: &
      TimeUnits

    real( kind=time_precision ) :: &
      time

    integer, dimension(3) :: &
      var_dim_zm, var_dim_zt

    integer :: i, k, n, status
    
    real( kind = core_rknd ), dimension(:,:), allocatable, save :: &
      wpthlp_out, &
      wprtp_out, &
      wp2_out, &
      thlp2_out, &
      rtp2_out, &
      rtpthlp_out, &
      upwp_out, &
      vpwp_out, &
      up2_out, &
      vp2_out

    real( kind = core_rknd ), dimension(:,:), allocatable, save :: &
      wp3_out, &
      rcm_out, &
      cloud_frac_out, &
      rtm_out, &
      thlm_out

    integer, dimension(ngrdcol) :: &
      column_list

    integer, save :: &
      n_calls   = 0, &  ! Number of total calls to this procedure
      n_outputs = 1, &  ! Number of times netcdf writing has been done
      n_samples = 0     ! Number of samples since last netcdf write

    integer, dimension(3) :: &
      start

    integer, dimension(2) :: &
      count

    real( kind = core_rknd ) :: &
      sample_weight

    integer, save :: &
        NF90_PREC

    !------------------------ Begin Code ------------------------

    ! If this is the first call, we have to create and define the netcdf files
    if ( n_calls == 0 ) then

      ! Set either double or single precision
      if ( l_output_double_prec ) then
        NF90_PREC = NF90_DOUBLE
      else
        NF90_PREC = NF90_FLOAT
      end if

      multicol_nc_file_zm = trim(output_file_prefix) // "_multi_col_zm.nc"
      multicol_nc_file_zt = trim(output_file_prefix) // "_multi_col_zt.nc"

      ! Get a formatted date
      call format_date( day, month, year, time_initial, & ! intent(in)
                        TimeUnits ) ! intent(out)

      !=================================== zm file definition ===================================

      ! Create the zm netcdf file
      status = nf90_create( path = multicol_nc_file_zm,  & 
                            cmode = NF90_CLOBBER,  & ! overwrite existing file
                            ncid = ncid_zm )

      ! Define the variable dimensions, columns, altitude, and time
      status = nf90_def_dim( ncid_zm, 'columns',   ngrdcol,        column_id_zm )
      status = nf90_def_dim( ncid_zm, 'altitude',  nzm,            vertical_id_zm )
      status = nf90_def_dim( ncid_zm, 'time',      NF90_UNLIMITED, time_id_zm )

      status = nf90_def_var( ncid_zm, "columns",   NF90_INT,      (/column_id_zm/),    column_var_id_zm )
      status = nf90_def_var( ncid_zm, "altitude",  NF90_PREC,  (/vertical_id_zm/),  vertical_var_id_zm )
      status = nf90_def_var( ncid_zm, "time",      NF90_PREC,      (/time_id_zm/),      time_var_id_zm )

      ! Define the attributes for the time variable
      status = nf90_put_att( ncid_zm, time_var_id_zm, "cartesian_axis", "T" )
      status = nf90_put_att( ncid_zm, time_var_id_zm, "units",          TimeUnits )
      status = nf90_put_att( ncid_zm, time_var_id_zm, "ipositive",      1 )
      status = nf90_put_att( ncid_zm, time_var_id_zm, "calendar_type", "Gregorian" )

      ! Define the attributes for the column variable
      status = nf90_put_att( ncid_zm, column_var_id_zm, "cartesian_axis", "X" )
      status = nf90_put_att( ncid_zm, column_var_id_zm, "units",          "degrees_E" )
      status = nf90_put_att( ncid_zm, column_var_id_zm, "ipositive",      1 )

      ! Define the attributes for the altitude variable
      status = nf90_put_att( ncid_zm, vertical_var_id_zm, "cartesian_axis",   "Z" )
      status = nf90_put_att( ncid_zm, vertical_var_id_zm, "units",            "meters" )
      status = nf90_put_att( ncid_zm, vertical_var_id_zm, "positive",         "up" )
      status = nf90_put_att( ncid_zm, vertical_var_id_zm, "ipositive",        1 )

      ! Define the variable dimensions, these are the dimension of the 
      ! variables we want to save
      var_dim_zm(1) = column_id_zm
      var_dim_zm(2) = vertical_id_zm
      var_dim_zm(3) = time_id_zm

      ! Define the zm variables to save
      status = nf90_def_var( ncid_zm, trim("wpthlp"),      NF90_PREC, var_dim_zm(:), ind_zm(1) )
      status = nf90_def_var( ncid_zm, trim("wprtp"),       NF90_PREC, var_dim_zm(:), ind_zm(2) )
      status = nf90_def_var( ncid_zm, trim("wp2"),         NF90_PREC, var_dim_zm(:), ind_zm(3) )
      status = nf90_def_var( ncid_zm, trim("thlp2"),       NF90_PREC, var_dim_zm(:), ind_zm(4) )
      status = nf90_def_var( ncid_zm, trim("rtp2"),        NF90_PREC, var_dim_zm(:), ind_zm(5) )
      status = nf90_def_var( ncid_zm, trim("rtpthlp"),     NF90_PREC, var_dim_zm(:), ind_zm(6) )
      status = nf90_def_var( ncid_zm, trim("upwp"),        NF90_PREC, var_dim_zm(:), ind_zm(7) )
      status = nf90_def_var( ncid_zm, trim("vpwp"),        NF90_PREC, var_dim_zm(:), ind_zm(8) )
      status = nf90_def_var( ncid_zm, trim("up2"),         NF90_PREC, var_dim_zm(:), ind_zm(9) )
      status = nf90_def_var( ncid_zm, trim("vp2"),         NF90_PREC, var_dim_zm(:), ind_zm(10) )

      ! End definition of file
      call nf_enddef(ncid_zm)

      !=================================== zt file definition ===================================
      
      ! Create the zt netcdf file
      status = nf90_create( path = multicol_nc_file_zt,  & 
                            cmode = NF90_CLOBBER,  & ! overwrite existing file
                            ncid = ncid_zt )

      ! Define the variable dimensions, columns, altitude, and time
      status = nf90_def_dim( ncid_zt, 'columns',   ngrdcol,        column_id_zt )
      status = nf90_def_dim( ncid_zt, 'altitude',  nzt,            vertical_id_zt )
      status = nf90_def_dim( ncid_zt, 'time',      NF90_UNLIMITED, time_id_zt )

      status = nf90_def_var( ncid_zt, "columns",   NF90_INT,      (/column_id_zt/),    column_var_id_zt )
      status = nf90_def_var( ncid_zt, "altitude",  NF90_PREC,  (/vertical_id_zt/),  vertical_var_id_zt )
      status = nf90_def_var( ncid_zt, "time",      NF90_PREC,      (/time_id_zt/),      time_var_id_zt )

      ! Define the attributes for the time variable
      status = nf90_put_att( ncid_zt, vertical_var_id_zt, "cartesian_axis",   "Z" )
      status = nf90_put_att( ncid_zt, vertical_var_id_zt, "units",            "meters" )
      status = nf90_put_att( ncid_zt, vertical_var_id_zt, "positive",         "up" )
      status = nf90_put_att( ncid_zt, vertical_var_id_zt, "ipositive",        1 )

      ! Define the attributes for the column variable
      status = nf90_put_att( ncid_zt, column_var_id_zt, "cartesian_axis", "X" )
      status = nf90_put_att( ncid_zt, column_var_id_zt, "units",          "degrees_E" )
      status = nf90_put_att( ncid_zt, column_var_id_zt, "ipositive",      1 )

      ! Define the attributes for the altitude variable
      status = nf90_put_att( ncid_zt, time_var_id_zt, "cartesian_axis", "T" )
      status = nf90_put_att( ncid_zt, time_var_id_zt, "units",          TimeUnits )
      status = nf90_put_att( ncid_zt, time_var_id_zt, "ipositive",      1 )
      status = nf90_put_att( ncid_zt, time_var_id_zt, "calendar_type", "Gregorian" )

      ! Define the variable dimensions, these are the dimension of the 
      ! variables we want to save
      var_dim_zt(1) = column_id_zt
      var_dim_zt(2) = vertical_id_zt
      var_dim_zt(3) = time_id_zt

      ! Define the zt variables to save
      status = nf90_def_var( ncid_zt, trim("wp3"),         NF90_PREC, var_dim_zt(:), ind_zt(1) )
      status = nf90_def_var( ncid_zt, trim("rcm"),         NF90_PREC, var_dim_zt(:), ind_zt(2) )
      status = nf90_def_var( ncid_zt, trim("cloud_frac"),  NF90_PREC, var_dim_zt(:), ind_zt(3) )
      status = nf90_def_var( ncid_zt, trim("rtm"),         NF90_PREC, var_dim_zt(:), ind_zt(4) )
      status = nf90_def_var( ncid_zt, trim("thlm"),        NF90_PREC, var_dim_zt(:), ind_zt(5) )

      ! End definition of file
      call nf_enddef(ncid_zt)

      !=================================== end definitions ===================================

      ! Store the data defining vertical levels
      status = nf90_put_var( ncid_zm, vertical_var_id_zm, gr%zm(1,:) )
      status = nf90_put_var( ncid_zt, vertical_var_id_zt, gr%zt(1,:)  )

      ! Store the data for the columns, just label each column in order 
      do i = 1, ngrdcol
        column_list(i) = i
      end do

      status = nf90_put_var( ncid_zm, column_var_id_zm, column_list, start=(/1/), count=(/ngrdcol/) )
      status = nf90_put_var( ncid_zt, column_var_id_zt, column_list, start=(/1/), count=(/ngrdcol/) )

      !=================================== Allocate Output Arrays ===================================

      allocate( wpthlp_out(ngrdcol,nzm), &
                wprtp_out(ngrdcol,nzm), &
                wp2_out(ngrdcol,nzm), &
                thlp2_out(ngrdcol,nzm), &
                rtp2_out(ngrdcol,nzm), &
                rtpthlp_out(ngrdcol,nzm), &
                upwp_out(ngrdcol,nzm), &
                vpwp_out(ngrdcol,nzm), &
                up2_out(ngrdcol,nzm), &
                vp2_out(ngrdcol,nzm) )

      allocate( wp3_out(ngrdcol,nzt), &
                rcm_out(ngrdcol,nzt), &
                cloud_frac_out(ngrdcol,nzt), &
                rtm_out(ngrdcol,nzt), &
                thlm_out(ngrdcol,nzt) )

      !$acc enter data create( wpthlp_out, wprtp_out, wp2_out, thlp2_out, rtp2_out, &
      !$acc                    rtpthlp_out, upwp_out, vpwp_out, up2_out, vp2_out, wp3_out, &
      !$acc                    rcm_out, cloud_frac_out, rtm_out, thlm_out )

    end if

    ! If there haven't been any samples in this output period yet, zero the output arrays
    if ( n_samples == 0 ) then

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          wpthlp_out(i,k)   = zero
          wprtp_out(i,k)    = zero
          wp2_out(i,k)      = zero
          thlp2_out(i,k)    = zero
          rtp2_out(i,k)     = zero
          rtpthlp_out(i,k)  = zero
          upwp_out(i,k)     = zero
          vpwp_out(i,k)     = zero
          up2_out(i,k)      = zero
          vp2_out(i,k)      = zero
        end do
      end do

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          wp3_out(i,k)        = zero
          rcm_out(i,k)        = zero
          cloud_frac_out(i,k) = zero
          rtm_out(i,k)        = zero
          thlm_out(i,k)       = zero
        end do
      end do

    end if

    ! Increment sample count and call count
    n_samples = n_samples + 1
    n_calls   = n_calls + 1

    !  time = n_calls * dt
    time = real( n_calls, kind=time_precision ) * real( dt, kind=time_precision )  ! seconds

    ! Add the new zm sample
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol
        wpthlp_out(i,k)   = wpthlp_out(i,k)  + wpthlp(i,k)
        wprtp_out(i,k)    = wprtp_out(i,k)   + wprtp(i,k)
        wp2_out(i,k)      = wp2_out(i,k)     + wp2(i,k)
        thlp2_out(i,k)    = thlp2_out(i,k)   + thlp2(i,k)
        rtp2_out(i,k)     = rtp2_out(i,k)    + rtp2(i,k)
        rtpthlp_out(i,k)  = rtpthlp_out(i,k) + rtpthlp(i,k)
        upwp_out(i,k)     = upwp_out(i,k)    + upwp(i,k)
        vpwp_out(i,k)     = vpwp_out(i,k)    + vpwp(i,k)
        up2_out(i,k)      = up2_out(i,k)     + up2(i,k)
        vp2_out(i,k)      = vp2_out(i,k)     + vp2(i,k)
      end do
    end do

    ! Add the new zt sample
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol
        wp3_out(i,k)        = wp3_out(i,k)        + wp3(i,k)
        rcm_out(i,k)        = rcm_out(i,k)        + rcm(i,k)
        cloud_frac_out(i,k) = cloud_frac_out(i,k) + cloud_frac(i,k)
        rtm_out(i,k)        = rtm_out(i,k)        + rtm(i,k)
        thlm_out(i,k)       = thlm_out(i,k)       + thlm(i,k)
      end do
    end do

    ! Perform a netcdf write every "calls_per_out" timsteps, or if this is the last timestep
    if ( mod( n_calls, calls_per_out ) == 0 .or. l_last_timestep ) then

      ! Calculate the sample weight, simply the inverse of the number of samples
      sample_weight = one / real( n_samples, kind = core_rknd )

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          wpthlp_out(i,k)   = sample_weight * wpthlp_out(i,k)  
          wprtp_out(i,k)    = sample_weight * wprtp_out(i,k)   
          wp2_out(i,k)      = sample_weight * wp2_out(i,k)     
          thlp2_out(i,k)    = sample_weight * thlp2_out(i,k)   
          rtp2_out(i,k)     = sample_weight * rtp2_out(i,k)    
          rtpthlp_out(i,k)  = sample_weight * rtpthlp_out(i,k) 
          upwp_out(i,k)     = sample_weight * upwp_out(i,k)    
          vpwp_out(i,k)     = sample_weight * vpwp_out(i,k)    
          up2_out(i,k)      = sample_weight * up2_out(i,k)     
          vp2_out(i,k)      = sample_weight * vp2_out(i,k)     
        end do
      end do

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          wp3_out(i,k)        = sample_weight * wp3_out(i,k)        
          rcm_out(i,k)        = sample_weight * rcm_out(i,k)        
          cloud_frac_out(i,k) = sample_weight * cloud_frac_out(i,k) 
          rtm_out(i,k)        = sample_weight * rtm_out(i,k)        
          thlm_out(i,k)       = sample_weight * thlm_out(i,k)       
        end do
      end do

      !$acc update host( wpthlp_out, wprtp_out, wp2_out, thlp2_out, rtp2_out, rtpthlp_out, &
      !$acc              upwp_out, vpwp_out, up2_out, vp2_out, wp3_out, rcm_out, cloud_frac_out, rtm_out, thlm_out )

      ! Update the time variables
      status = nf90_put_var( ncid_zm, time_var_id_zm, time, (/n_outputs/) )  
      status = nf90_put_var( ncid_zt, time_var_id_zt, time, (/n_outputs/) )

      start = (/ 1, 1, n_outputs /)
      
      if ( l_output_double_prec ) then

        count = (/ ngrdcol, nzm /)

        ! Update the zm variables we want, needs to match the variables we defined 
        status = nf90_put_var( ncid_zm, ind_zm(1),       wpthlp_out, start, count )
        status = nf90_put_var( ncid_zm, ind_zm(2),        wprtp_out, start, count )
        status = nf90_put_var( ncid_zm, ind_zm(3),          wp2_out, start, count )
        status = nf90_put_var( ncid_zm, ind_zm(4),        thlp2_out, start, count )
        status = nf90_put_var( ncid_zm, ind_zm(5),         rtp2_out, start, count )
        status = nf90_put_var( ncid_zm, ind_zm(6),      rtpthlp_out, start, count )
        status = nf90_put_var( ncid_zm, ind_zm(7),         upwp_out, start, count )
        status = nf90_put_var( ncid_zm, ind_zm(8),         vpwp_out, start, count )
        status = nf90_put_var( ncid_zm, ind_zm(9),          up2_out, start, count )
        status = nf90_put_var( ncid_zm, ind_zm(10),         vp2_out, start, count )

        count = (/ ngrdcol, nzt /)

        ! Update the zt variables we want, needs to match the variables we defined 
        status = nf90_put_var( ncid_zt, ind_zt(1),          wp3_out, start, count )
        status = nf90_put_var( ncid_zt, ind_zt(2),          rcm_out, start, count )
        status = nf90_put_var( ncid_zt, ind_zt(3),   cloud_frac_out, start, count )
        status = nf90_put_var( ncid_zt, ind_zt(4),          rtm_out, start, count )
        status = nf90_put_var( ncid_zt, ind_zt(5),         thlm_out, start, count )

      else

        count = (/ ngrdcol, nzm /)

        ! Update the zm variables we want, needs to match the variables we defined 
        status = nf90_put_var( ncid_zm, ind_zm(1),  real(  wpthlp_out, kind=sp), start, count )
        status = nf90_put_var( ncid_zm, ind_zm(2),  real(   wprtp_out, kind=sp), start, count )
        status = nf90_put_var( ncid_zm, ind_zm(3),  real(     wp2_out, kind=sp), start, count )
        status = nf90_put_var( ncid_zm, ind_zm(4),  real(   thlp2_out, kind=sp), start, count )
        status = nf90_put_var( ncid_zm, ind_zm(5),  real(    rtp2_out, kind=sp), start, count )
        status = nf90_put_var( ncid_zm, ind_zm(6),  real( rtpthlp_out, kind=sp), start, count )
        status = nf90_put_var( ncid_zm, ind_zm(7),  real(    upwp_out, kind=sp), start, count )
        status = nf90_put_var( ncid_zm, ind_zm(8),  real(    vpwp_out, kind=sp), start, count )
        status = nf90_put_var( ncid_zm, ind_zm(9),  real(     up2_out, kind=sp), start, count )
        status = nf90_put_var( ncid_zm, ind_zm(10), real(     vp2_out, kind=sp), start, count )

        count = (/ ngrdcol, nzt /)

        ! Update the zt variables we want, needs to match the variables we defined 
        status = nf90_put_var( ncid_zt, ind_zt(1), real(        wp3_out, kind=sp), start, count )
        status = nf90_put_var( ncid_zt, ind_zt(2), real(        rcm_out, kind=sp), start, count )
        status = nf90_put_var( ncid_zt, ind_zt(3), real( cloud_frac_out, kind=sp), start, count )
        status = nf90_put_var( ncid_zt, ind_zt(4), real(        rtm_out, kind=sp), start, count )
        status = nf90_put_var( ncid_zt, ind_zt(5), real(       thlm_out, kind=sp), start, count )
        
      end if
        
      ! Reset the samples counter to 0, and increment the output counter
      n_samples = 0
      n_outputs = n_outputs + 1

    end if

    if ( l_last_timestep ) then

      ! Close netcdf file
      status = nf90_close( ncid = ncid_zm )
      status = nf90_close( ncid = ncid_zt )

      !$acc exit data delete( wpthlp_out, wprtp_out, wp2_out, thlp2_out, rtp2_out, &
      !$acc                   rtpthlp_out, upwp_out, vpwp_out, up2_out, vp2_out, wp3_out, &
      !$acc                   rcm_out, cloud_frac_out, rtm_out, thlm_out )

    end if

  end subroutine output_multi_col_fields
#endif

end module output_netcdf
