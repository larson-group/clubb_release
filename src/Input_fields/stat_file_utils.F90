! $Id$
module stat_file_utils

  implicit none

  public :: stat_file_average, stat_file_average_interval, &
    stat_file_num_vertical_levels, stat_file_vertical_levels, &
    LES_grid_to_CLUBB_grid, CLUBB_levels_within_LES_domain, l_netcdf_file
        
  private ! Default Scope

  contains

!----------------------------------------------------------------------
  function stat_file_average( filename, out_nz, ngrdcol, & ! In
                              t1, t2, out_heights, & ! In
                              variable_name, & ! In
                              npower, l_spec_bound_cond, & ! In
                              l_error ) ! Out
! Description:
!   Read each shared file profile once and average it on every column grid.
!   Average a netCDF file variable over the interval t1 to t2

! References:
!   None
!----------------------------------------------------------------------

    use constants_clubb, only: fstderr, eps ! Variable(s)

    use stat_file_module, only: stat_file ! Type(s)

    use input_netcdf, only: &
      open_netcdf_read, get_netcdf_var, & ! Procedures
      close_netcdf_read

!   use extrapolation, only: &
!     lin_ext_zt_bottom

    use interpolation, only: &
      lin_interpolate_two_points

    use clubb_precision, only: &
      core_rknd, & ! Variable(s)
      sp

    implicit none

    ! External
    intrinsic :: transfer

    ! Input Variables
    character(len=*), intent(in) ::  & 
      filename ! Name of the file

    integer, intent(in) ::  & 
      out_nz,  & ! Number of output vertical levels [-]
      ngrdcol, & ! Number of output columns [-]
      t1,      & ! Beginning timestep to look at
      t2         ! Ending timestep to look at

    real( kind = core_rknd ), dimension(ngrdcol,out_nz), intent(in) :: &
      out_heights ! Heights of the output grid [m]

    character(len=*), intent(in) :: & 
      variable_name ! Name of the variable to read in

    integer, intent(in) :: & 
      npower ! Exponent operator, must be 1 or 2

    ! Special boundary condition for when the lowest model level is a ghost point
    ! comes from a fine vertical resolution run.  E.g.
    !     Stat file levels    Model levels
    !           50m                50m
    !           25m
    !          -25m
    !                             -50m
    ! So we need to set the lowest model level to the -25m value.
    logical, dimension(ngrdcol), intent(in) :: &
      l_spec_bound_cond 

    ! Output Variable
    logical, intent(out) :: & 
      l_error ! error status for this function

    ! Return Variable for function
    real( kind = core_rknd ), dimension(ngrdcol,out_nz) :: stat_file_average

    ! Local Variables
    type (stat_file) :: faverage ! Data file derived type

    real( kind = core_rknd ), allocatable, dimension(:) :: file_variable ! Temporary variable

    real( kind = core_rknd ), dimension(ngrdcol,out_nz) ::  interp_variable ! Temporary variable

    integer :: & 
      t, &  ! Timestep loop index
      k, &  ! Vertical loop index
      i     ! Column loop index

    integer :: &
      file_nz ! Number of vertical levels in the file

    logical, dimension(ngrdcol) :: l_interpolate

    logical, dimension(ngrdcol,out_nz) :: l_lin_int

    integer, dimension(ngrdcol,out_nz) :: &
      exact_lev_idx, &
      lower_lev_idx, &
      upper_lev_idx

    integer :: & 
      num_timesteps ! Steps between t1 and t2

    integer, dimension(ngrdcol) :: &
      k_lowest_input, & ! First output level within the file domain [-]
      k_highest_input   ! Last output level within the file domain [-]

    logical :: &
      l_convert_to_MKS ! convert inputs to MKS units

    real( kind = sp ) :: &
      PosInf = transfer( 2139095040, 1.0_sp ) ! 2139095040 is a magic number
    
!-----------------------------------------------------------------------

    ! Initialize variables
    l_error = .false.

    num_timesteps = ( t2 - t1 ) + 1
    stat_file_average = 0.0_core_rknd

    ! Open netCDF file
    call open_netcdf_read( variable_name, filename, faverage, l_error )

    if ( l_error ) return

    ! Determine variable size
    file_nz = size( faverage%z ) 

    allocate( file_variable(file_nz) ) 

    ! Determine the interpolation needed on each column grid.
    do i = 1, ngrdcol
      if ( out_nz /= file_nz ) then
        l_interpolate(i) = .true.
      else
        l_interpolate(i) = any( abs(faverage%z(:)-out_heights(i,:)) > &
                               abs(faverage%z(:)+out_heights(i,:))*eps/2 )
      end if

    end do

    call CLUBB_levels_within_LES_domain( faverage, out_nz, ngrdcol, out_heights, & ! In
                                         k_lowest_input, & ! Out
                                         k_highest_input ) ! Out

    do k = 1, out_nz
      do i = 1, ngrdcol
        if ( k < k_lowest_input(i) .or. k > k_highest_input(i) ) cycle
        call LES_grid_to_CLUBB_grid( faverage, out_heights(i,k), & ! In
                                     exact_lev_idx(i,k), lower_lev_idx(i,k), & ! Out
                                     upper_lev_idx(i,k), l_lin_int(i,k) ) ! Out
      end do
    end do

    ! Read in variables from netCDF file
    do t = t1, t2

      l_convert_to_MKS = .true.
      call get_netcdf_var( faverage, variable_name, t, l_convert_to_MKS, &
                           file_variable(1:file_nz), l_error )

      if ( l_error ) then
        write(fstderr,*) "stat_file_average: get_var failed for "  & 
          //trim( variable_name )//" in "//trim( filename )//" at time=", t
        return
      end if

      ! Interpolate within each column's domain and retain the legacy boundary policy.
      do k = 1, out_nz
        do i = 1, ngrdcol
          if ( .not. l_interpolate(i) ) then
            interp_variable(i,k) = file_variable(k)
          else if ( k < k_lowest_input(i) ) then
            if ( k == 1 .and. l_spec_bound_cond(i) ) then
              interp_variable(i,k) = file_variable(1)
            else
              interp_variable(i,k) = real(PosInf, kind = core_rknd)
            end if
          else if ( k > k_highest_input(i) ) then
            interp_variable(i,k) = real(PosInf, kind = core_rknd)
          else if ( l_lin_int(i,k) ) then
            interp_variable(i,k) = lin_interpolate_two_points( out_heights(i,k), &
              faverage%z(upper_lev_idx(i,k)), faverage%z(lower_lev_idx(i,k)), &
              file_variable(upper_lev_idx(i,k)), file_variable(lower_lev_idx(i,k)) )
          else
            interp_variable(i,k) = file_variable(exact_lev_idx(i,k))
          end if
        end do
      end do

      ! Apply an exponent (for the tuner), retaining the per-column time sum order.
      if ( npower /= 1 ) then
        stat_file_average = stat_file_average + interp_variable**npower
      else
        stat_file_average = stat_file_average + interp_variable
      end if

    end do ! t = t1, t2

    ! Close the netCDF file
    call close_netcdf_read( faverage )

    ! Take average over num_timesteps
    stat_file_average = stat_file_average / &
               real( num_timesteps, kind = core_rknd )

    return
  end function stat_file_average

!-------------------------------------------------------------------------
  function stat_file_average_interval & 
           ( filename, nz, t, variable_name, & 
             out_heights, npower, l_error )

! Description:
!   Reads in netCDF data from a file and then takes several averages
!   over an interval.

! References:
!   None

! Notes:
!   The variable t is assumed size, which needs to be used with caution.
!-------------------------------------------------------------------------
    use constants_clubb, only: fstderr ! Variable(s)

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Constant Parameters
    integer, parameter ::  & 
      tmax = huge( 1 ) ! Sanity check for huge t array

    ! Input Variables
    character(len=*), intent(in) ::  & 
      filename ! Name of the file

    integer, intent(in) ::  & 
      nz     ! Number of vertical grid levels in the host model

    integer, dimension(:), intent(in) ::  & 
      t ! Timesteps to use for taking an average over an interval

    character(len=*), intent(in) ::  & 
      variable_name  ! Name of the variable to read in

    real( kind = core_rknd ), dimension(nz), intent(in) :: out_heights

    integer, intent(in) ::  & 
      npower ! exponent applied to data retrieved from the file (1 or 2)

    ! Output Variables
    logical, intent(out) ::  & 
      l_error ! status of this function

    ! Return Variables
    real( kind = core_rknd ), dimension(nz) ::  & 
      stat_file_average_interval

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) :: stat_file_temp

    real( kind = core_rknd ), dimension(1,nz) :: averaged_profile ! Shared reference profile

    integer ::  & 
      i,       & ! Loop variable 
      tdim,    & ! Dimension to read over for t variable
      divisor

    logical, dimension(1) :: l_spec_bound_cond ! special boundary condition for stat_file_average

!-------------------------------------------------------------------------
    l_error = .false.

    ! Sanity check
    if ( size( t ) > tmax .or. size( t ) < 2 ) then
      write(unit=fstderr,fmt=*)  & 
        "stat_file_average_interval: Invalid time interval"
      l_error = .true.
      return
    end if

    tdim = -1
    do i = 1, size( t )
      if ( t( i ) == 0 ) exit
      tdim = i
    end do

    l_spec_bound_cond = .false.
    averaged_profile = stat_file_average( filename, nz, 1, & ! In
                                          t(1), t(2), reshape(out_heights, [1,nz]), & ! In
                                          variable_name, & ! In
                                          npower, l_spec_bound_cond, & ! In
                                          l_error ) ! Out
    stat_file_average_interval = averaged_profile(1,:) * real( t(2) - t(1), kind = core_rknd )

    divisor = t(2) - t(1)

    if ( l_error ) return

    do i=3, tdim, 2
      l_spec_bound_cond = .false.
      averaged_profile = stat_file_average( filename, nz, 1, & ! In
                                            t(i), t(i+1), reshape(out_heights, [1,nz]), & ! In
                                            variable_name, & ! In
                                            npower, l_spec_bound_cond, & ! In
                                            l_error ) ! Out
      stat_file_temp = averaged_profile(1,:)
      stat_file_average_interval  & 
        = stat_file_average_interval + stat_file_temp * real( t(i+1) - t(i), kind = core_rknd )
      divisor = divisor + ( t(i+1) - t(i) )
    end do

    stat_file_average_interval(1:nz)  = stat_file_average_interval(1:nz) /&
         real( divisor, kind = core_rknd )

    return
  end function stat_file_average_interval

!-------------------------------------------------------------------------
  integer function stat_file_num_vertical_levels( varname, filename )

!       Description:
!       Returns a integer for the number of vertical levels in a file

!       References:
!       None
!-------------------------------------------------------------------------

    use stat_file_module, only: stat_file ! Type(s)

    use constants_clubb, only: fstderr

    use input_netcdf, only: open_netcdf_read, close_netcdf_read ! Procedure(s)

    implicit none

    ! Input Variables
    character(len=*), intent(in) ::  & 
      filename, & ! File name
      varname     ! Variable name

    ! Local Variables
    type (stat_file) :: file_nz ! Data file

    logical :: l_error

    ! ---- Begin Code ----

    call open_netcdf_read( varname, filename, file_nz, l_error )

    if ( l_error ) then
      write(fstderr,*) "Error opening "// filename
      error stop
    end if

    ! Set return variable
    stat_file_num_vertical_levels = file_nz%iz

    ! Close file
    call close_netcdf_read( file_nz )

    return
  end function stat_file_num_vertical_levels

!-------------------------------------------------------------------------
  function stat_file_vertical_levels( varname, filename, nz )

    use stat_file_module, only: stat_file ! Type(s)

    use constants_clubb, only: fstderr

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use input_netcdf, only: open_netcdf_read, close_netcdf_read ! Procedure(s)

    implicit none

    ! Input Variables
    character(len=*), intent(in) ::  & 
      filename, & ! File name
      varname     ! Variable name

    integer, intent(in) :: &
      nz ! Number of vertical levels

    ! Output Variables
    real( kind = core_rknd ), dimension(nz) :: &
      stat_file_vertical_levels ! A vector containing the heights of each vertical level [m]

    ! Local Variables
    type (stat_file) :: file_levels  ! Data file

    logical :: l_error

    ! ---- Begin Code ----

    call open_netcdf_read( varname, filename, file_levels, l_error )

    if ( l_error ) then
      write(fstderr,*) "Error opening "// filename
      error stop
    end if

    ! Set return variable
    stat_file_vertical_levels(1:nz) = file_levels%z(1:nz)

    ! Close file
    call close_netcdf_read( file_levels )

    return
  end function stat_file_vertical_levels
!-------------------------------------------------------------------------
  logical function l_netcdf_file( filename )

    use constants_clubb, only : fstderr

    implicit none

    ! External
    intrinsic :: len, trim

    ! Input Variables
    character(len=*), intent(in) :: &
      filename

    ! Local Variables
    integer :: length

    ! ---- Begin Code ----

    length = len( trim( filename ) )

    if ( filename(length-2:length) == ".nc" ) then
      l_netcdf_file = .true.
    else
      l_netcdf_file = .false.
      write(fstderr,*) "Unrecognized file type: "//trim( filename )
      error stop
    end if

    return
  end function l_netcdf_file

!===============================================================================
  subroutine CLUBB_levels_within_LES_domain( fread_var, nz, ngrdcol, CLUBB_grid, & ! In
                                             k_lowest_input, & ! Out
                                             k_highest_input ) ! Out

    ! Description:
    !   Find each column's lowest and highest grid levels within the LES domain.
    !   The input column grids are ascending, as in the profile reader.
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use stat_file_module, only: stat_file ! Type

    use constants_clubb, only: fstderr ! Error output unit

    use clubb_precision, only: core_rknd ! Precision

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(stat_file), intent(in) :: fread_var ! LES grid and valid level indices

    integer, intent(in) :: &
      nz, & ! Number of CLUBB vertical levels [-]
      ngrdcol ! Number of columns [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      CLUBB_grid ! Column heights on thermodynamic or momentum levels [m]

    !--------------------------- Output Variables ---------------------------
    integer, dimension(ngrdcol), intent(out) :: &
      k_lowest_input, & ! Lowest CLUBB level within the LES domain [-]
      k_highest_input ! Highest CLUBB level within the LES domain [-]

    !--------------------------- Local Variables ---------------------------
    integer :: i, k ! Column and vertical indices

    !--------------------------- Begin Code ---------------------------

    if ( any( CLUBB_grid(:,nz) < fread_var%z(fread_var%ia) ) ) then
      write(fstderr,*) "The lowest LES input level is above the top ", &
                       "of the CLUBB model domain."
      error stop "Error in CLUBB_levels_within_LES_domain"
    end if

    if ( any( CLUBB_grid(:,1) > fread_var%z(fread_var%iz) ) ) then
      write(fstderr,*) "The highest LES input level is below the bottom ", &
                       "of the CLUBB model domain."
      error stop "Error in CLUBB_levels_within_LES_domain"
    end if

    k_lowest_input = 1
    k_highest_input = nz

    do k = 1, nz
      do i = 1, ngrdcol
        if ( CLUBB_grid(i,k) < fread_var%z(fread_var%ia) ) k_lowest_input(i) = k + 1
        if ( CLUBB_grid(i,k) > fread_var%z(fread_var%iz) ) &
          k_highest_input(i) = min(k_highest_input(i), k - 1)
      end do
    end do

  end subroutine CLUBB_levels_within_LES_domain

!===============================================================================
  pure subroutine LES_grid_to_CLUBB_grid( fread_var, CLUBB_height,  &
                                          exact_lev_idx, lower_lev_idx,  &
                                          upper_lev_idx, l_lin_int )

    ! Description:
    !   Finds the level on the LES grid that is exactly even with the CLUBB
    !   grid level (either thermodynamic or momentum grid) that is input
    !   height. Else, it finds the two LES levels that sandwich the CLUBB
    !   height that is input.

    !-----------------------------------------------------------------------

    use stat_file_module, only:  &
        stat_file  ! Variable type

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables.
    type(stat_file), intent(in) ::  &
      fread_var  ! Information about LES run.

    real( kind = core_rknd ), intent(in) ::  &
      CLUBB_height ! Altitude of one CLUBB grid level
                 ! (either thermodynamic or momentum grid levels)  [m]

    ! Output Variables.
    integer, intent(out) ::  &
      exact_lev_idx, & ! In case of an exact match, index of LES level that is
                       ! exactly even with CLUBB level k.
      lower_lev_idx, & ! In case linear interpolation is needed, index of LES
                       ! level that is immediately below CLUBB level k.
      upper_lev_idx    ! In case linear interpolation is needed, index of LES
                       ! level that is immediately above CLUBB level k.

    logical, intent(out) ::  &
      l_lin_int  ! Flag that is turned on if linear interpolation is needed.

    ! Local Variable.
    integer :: j

    ! Initialize the output quantities.
    exact_lev_idx = 0
    lower_lev_idx = 0
    upper_lev_idx = 0
    l_lin_int     = .false.

    ! Initialize LES vertical grid loop index, j, at the lowest LES grid index,
    ! which is fread_var%ia.

    do j = fread_var%ia, fread_var%iz, 1

      if ( abs( fread_var%z(j) - CLUBB_height ) < 1e-6_core_rknd ) then

        ! There is an LES level altitude at LES level j that is an exact
        ! match to the CLUBB level altitude at CLUBB grid level k.
        exact_lev_idx = j
        l_lin_int = .false.

      else if ( fread_var%z(j) < CLUBB_height ) then

        ! The LES level altitude at LES level j is lower than the CLUBB level
        ! altitude at CLUBB grid level k.
        lower_lev_idx = j

      else  ! fread_var%z(j) > CLUBB_height

        ! The LES level altitude at LES level j is higher than the CLUBB level
        ! altitude at CLUBB grid level k.
        upper_lev_idx = j

      end if

      if ( exact_lev_idx > 0 ) exit  ! An exact answer has been found,
      ! exit the loop.

      if ( upper_lev_idx == lower_lev_idx + 1 ) then

        ! CLUBB level k has been found between two successive LES levels.
        ! Linear interpolation is needed.  An answer has been found, exit
        ! the loop.
        l_lin_int = .true.
        exit

      end if

    end do

    return
  end subroutine LES_grid_to_CLUBB_grid

end module stat_file_utils
