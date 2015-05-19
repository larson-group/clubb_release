! $Id$
module stat_file_utils

  implicit none

  public :: stat_file_average, stat_file_average_interval, &
    stat_file_num_vertical_levels, stat_file_vertical_levels, &
    LES_grid_to_CLUBB_grid, CLUBB_levels_within_LES_domain, l_netcdf_file
        
  private ! Default Scope

  contains

!----------------------------------------------------------------------
  function stat_file_average( filename, out_nz, &
                              t1, t2, out_heights, variable_name,  & 
                              npower, l_spec_bound_cond, l_error )
! Description:
!   Average a GrADS or netCDF file variable over the interval t1 to t2

! References:
!   None
!----------------------------------------------------------------------

    use constants_clubb, only: fstderr ! Variable(s)

    use parameters_model, only: PosInf ! Variable(s)

    use stat_file_module, only: stat_file ! Type(s)

    use input_grads, only: &
      open_grads_read, get_grads_var,  & ! Procedures
      close_grads_read

#ifdef NETCDF
    use input_netcdf, only: &
      open_netcdf_read, get_netcdf_var, & ! Procedures
      close_netcdf_read
#endif

!   use extrapolation, only: &
!     lin_ext_zt_bottom

    use interpolation, only: &
      lin_interpolate_two_points

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: transfer

#ifdef _OPENMP
    integer :: omp_get_thread_num ! Function
#endif

    ! Input Variables
    character(len=*), intent(in) ::  & 
      filename ! Name of the file

    integer, intent(in) ::  & 
      out_nz,  & ! Number of vertial levels in the GrADS file.
      t1,      & ! Beginning timestep to look at
      t2         ! Ending timestep to look at

    real( kind = core_rknd ), dimension(out_nz), intent(in) :: &
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
    logical, intent(in) :: &
      l_spec_bound_cond 

    ! Output Variable
    logical, intent(out) :: & 
      l_error ! error status for this function

    ! Return Variable for function
    real( kind = core_rknd ), dimension(out_nz) :: stat_file_average

    ! Local Variables
    type (stat_file) :: faverage ! Data file derived type

    real( kind = core_rknd ), allocatable, dimension(:) :: file_variable ! Temporary variable

    real( kind = core_rknd ), dimension(out_nz) ::  interp_variable ! Temporary variable

    integer :: & 
      t, &  ! Timestep loop index
      k     ! Vertical loop index

    integer :: &
      file_nz, & ! Number of vertical levels in the file
      iunit      ! File I/O unit

    logical :: l_interpolate, l_grads_file

    logical, dimension(out_nz) :: l_lin_int

    integer, dimension(out_nz) :: &
      exact_lev_idx, &
      lower_lev_idx, &
      upper_lev_idx

    integer :: & 
      num_timesteps,  & ! steps between t1 and t2
      k_lowest_input, &
      k_highest_input

    logical :: &
      l_convert_to_MKS ! convert inputs to MKS units

!-----------------------------------------------------------------------

    ! Initialize variables
    l_error = .false.

    num_timesteps = ( t2 - t1 ) + 1
    stat_file_average = 0.0_core_rknd

    ! Determine file type
    l_grads_file = .not. l_netcdf_file( filename ) 
    
    ! Open GraDS file
    if ( l_grads_file ) then
#ifdef _OPENMP
    iunit = omp_get_thread_num( ) + 10 ! Known magic number
#else
    iunit = 10
#endif
      call open_grads_read( iunit, filename, faverage, l_error )

    else

#ifdef NETCDF
      call open_netcdf_read( variable_name, filename, faverage, l_error )

#else
      write(fstderr,*) "This version of CLUBB was not compiled with netCDF support"
      l_error = .true.
#endif

    end if ! l_grads_file

    if ( l_error ) return

    ! Determine variable size
    file_nz = size( faverage%z ) 

    allocate( file_variable(file_nz) ) 

    ! Do we need to interpolate?
    if ( out_nz /= file_nz ) then
      l_interpolate = .true.
    else if ( any( faverage%z(:) /= out_heights(:) ) ) then
      l_interpolate = .true.
    else
      l_interpolate = .false.
    end if

    call CLUBB_levels_within_LES_domain( faverage, out_heights,  &
                                         k_lowest_input, k_highest_input )

    do  k = k_lowest_input, k_highest_input, 1
      ! CLUBB vertical level k is found at an altitude that is within the
      ! domain of the LES output.
      call LES_grid_to_CLUBB_grid( faverage, out_heights, k,  &
                                   exact_lev_idx(k), lower_lev_idx(k),  &
                                   upper_lev_idx(k), l_lin_int(k) )
    end do

    ! Read in variables from GrADS file
    do t = t1, t2

      if ( l_grads_file ) then
        call get_grads_var( faverage, variable_name, t,  & 
                            file_variable(1:file_nz), l_error )
      else
#ifdef NETCDF
        l_convert_to_MKS = .true.
        call get_netcdf_var( faverage, variable_name, t, l_convert_to_MKS, &
                             file_variable(1:file_nz), l_error )
#endif
      end if


      if ( l_error ) then
        write(fstderr,*) "stat_file_average: get_var failed for "  & 
          //trim( variable_name )//" in "//trim( filename )//" at time=", t
        return
      end if

      ! Interpolate as needed using the inputfields code
      if ( l_interpolate ) then

        do k = k_lowest_input, k_highest_input, 1
          if ( l_lin_int(k) ) then
            interp_variable(k) = lin_interpolate_two_points( out_heights(k), &
              faverage%z(upper_lev_idx(k)), faverage%z(lower_lev_idx(k)), &
              file_variable(upper_lev_idx(k)), file_variable(lower_lev_idx(k)) )
          else
            interp_variable(k) = file_variable(exact_lev_idx(k))
          end if
        end do

        do k = k_lowest_input-1, 1, -1
          ! Do a linear extension on the lower points
!         interp_variable(k) = lin_ext_zt_bottom( interp_variable(k+4), &
!           interp_variable(k+1), out_heights(k+2), out_heights(k+1), out_heights(k) )

          ! Set the ghost point if we have a special boundary condition
          if ( k == 1 .and. l_spec_bound_cond ) then
            interp_variable(k) = file_variable(1)
          else
            ! Set undefined points to NaN
            interp_variable(k) = real(PosInf, kind = core_rknd)
          end if
        end do

        do k = k_highest_input+1, out_nz, 1
          ! Do the the dum-dum thing and set points above output domain 
          ! to be constant with height
!         interp_variable(k) = interp_variable(k-1)

          ! Set undefined points to NaN
          interp_variable(k) = real(PosInf, kind = core_rknd)
        end do

      else
        interp_variable(1:out_nz) = file_variable(1:out_nz)

      end if


      ! Apply an exponent (for the tuner)
      if ( npower /= 1 ) then
        stat_file_average(1:out_nz) = stat_file_average(1:out_nz) &
        + interp_variable(1:out_nz)**npower

      else
        stat_file_average(1:out_nz) = stat_file_average(1:out_nz) &
        + interp_variable(1:out_nz)

      end if

    end do ! t = t1, t2

    if ( l_grads_file ) then
      ! Close the GrADS file
      call close_grads_read( faverage )
    else
#ifdef NETCDF
      ! Close the netCDF file
      call close_netcdf_read( faverage )
#endif
    end if

    ! Take average over num_timesteps
    stat_file_average(1:out_nz) = stat_file_average(1:out_nz) / &
               real( num_timesteps, kind = core_rknd )

    return
  end function stat_file_average

!-------------------------------------------------------------------------
  function stat_file_average_interval & 
           ( filename, nz, t, variable_name, & 
             out_heights, npower, l_error )

! Description:
!   Reads in GrADS data from a file and then takes several averages
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

    integer ::  & 
      i,       & ! Loop variable 
      tdim,    & ! Dimension to read over for t variable
      divisor

    logical :: l_spec_bound_cond ! special boundary condition for stat_file_average

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
    stat_file_average_interval & 
    = stat_file_average( filename, nz, &
                         t(1), t(2), out_heights, variable_name, &
                         npower, l_spec_bound_cond, l_error )  & 
          * real( t(2) - t(1), kind = core_rknd )

    divisor = t(2) - t(1)

    if ( l_error ) return

    do i=3, tdim, 2
      l_spec_bound_cond = .false.
      stat_file_temp = stat_file_average( filename, nz, &
                               t(i), t(i+1), out_heights, variable_name, &
                               npower, l_spec_bound_cond, l_error )
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

!       Notes:
!       Somewhat inefficient, since it reads in the entire .ctl file to
!       determine the number of levels
!-------------------------------------------------------------------------

    use stat_file_module, only: stat_file ! Type(s)

    use input_grads, only: open_grads_read, close_grads_read ! Procedure(s)

    use constants_clubb, only: fstderr

#ifdef NETCDF
    use input_netcdf, only: open_netcdf_read, close_netcdf_read ! Procedure(s)
#endif

    implicit none

    ! External
#ifdef _OPENMP
    integer :: omp_get_thread_num ! Function
#endif

    ! Input Variables
    character(len=*), intent(in) ::  & 
      filename, & ! File name
      varname     ! Variable name

    ! Local Variables
    type (stat_file) :: file_nz ! Data file

    logical :: l_grads_file, l_error

    integer :: iunit ! File unit number

    ! ---- Begin Code ----

    l_grads_file = .not. l_netcdf_file( filename )  

    if ( l_grads_file ) then

#ifdef _OPENMP
    iunit = omp_get_thread_num( ) + 10 ! Known magic number
#else
    iunit = 10
#endif
      ! Read in the control file
      call open_grads_read( iunit, filename, file_nz, l_error )

    else

#ifdef NETCDF
      call open_netcdf_read( varname, filename, file_nz, l_error )
#else
      write(fstderr,*) "This version of CLUBB was not compiled with netCDF support"
      l_error = .true.
#endif
      if ( l_error ) then
        write(fstderr,*) "Error opening "// filename
        stop
      end if

    end if

    ! Set return variable
    stat_file_num_vertical_levels = file_nz%iz

    ! Close file
    if ( l_grads_file ) then
      call close_grads_read( file_nz )
    else
#ifdef NETCDF
      call close_netcdf_read( file_nz )
#endif
    end if

    return
  end function stat_file_num_vertical_levels

!-------------------------------------------------------------------------
  function stat_file_vertical_levels( varname, filename, nz )

    use stat_file_module, only: stat_file ! Type(s)

    use input_grads, only: open_grads_read, close_grads_read ! Procedure(s)

    use constants_clubb, only: fstderr

    use clubb_precision, only: &
      core_rknd ! Variable(s)

#ifdef NETCDF
    use input_netcdf, only: open_netcdf_read, close_netcdf_read ! Procedure(s)
#endif

    implicit none

    ! External
#ifdef _OPENMP
    integer :: omp_get_thread_num ! Function
#endif

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

    logical :: l_grads_file, l_error

    integer :: iunit ! file unit number

    ! ---- Begin Code ----

    l_grads_file = .not. l_netcdf_file( filename )  

    if ( l_grads_file ) then
#ifdef _OPENMP
    iunit = omp_get_thread_num( ) + 10 ! Known magic number
#else
    iunit = 10
#endif
      ! Read in the control file
      call open_grads_read( iunit, filename, file_levels, l_error )

    else

#ifdef NETCDF
      call open_netcdf_read( varname, filename, file_levels, l_error )
#else
      write(fstderr,*) "This version of CLUBB was not compiled with netCDF support"
      l_error = .true.
#endif
      if ( l_error ) then
        write(fstderr,*) "Error opening "// filename
        stop
      end if

    end if

    ! Set return variable
    stat_file_vertical_levels(1:nz) = file_levels%z(1:nz)

    ! Close file
    if ( l_grads_file ) then
      call close_grads_read( file_levels )
    else
#ifdef NETCDF
      call close_netcdf_read( file_levels )
#endif
    end if

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
    else if ( filename(length-3:length) == ".ctl" ) then
      l_netcdf_file = .false.
    else
      l_netcdf_file = .false.
      write(fstderr,*) "Unrecognized file type: "//trim( filename )
      stop
    end if

    return
  end function l_netcdf_file

!===============================================================================
  subroutine CLUBB_levels_within_LES_domain( fread_var, CLUBB_grid,  &
                                             k_lowest_input, k_highest_input )

    ! Description:
    ! This subroutine finds both the lowest and the highest CLUBB grid levels
    ! (for either the thermodynamic grid or the momentum grid) that are with the
    ! domain of the LES grid.

    ! References:
    !   None
    !-----------------------------------------------------------------------


    use stat_file_module, only:  &
        stat_file  ! Variable type

    use constants_clubb, only:  &
        fstderr ! Constant

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables.
    type(stat_file), intent(in) ::  &
      fread_var  ! Information about LES run.

    real( kind = core_rknd ), dimension(:), intent(in) ::  &
      CLUBB_grid ! Altitude of CLUBB grid levels
                 ! (either thermodynamic or momentum grid levels)  [m]

    ! Output Variables
    integer, intent(out) ::  &
      k_lowest_input,  & ! The lowest CLUBB level that's within the LES domain.
      k_highest_input    ! The highest CLUBB level that's within the LES domain.

    ! Local Variable
    integer :: k, kmax  ! Array index

    ! ---- Begin Code ----

    ! Find the lowest CLUBB level that falls within the LES domain.
    k    = size( CLUBB_grid )
    kmax = size( CLUBB_grid )
    do
      if ( CLUBB_grid(k) < fread_var%z(fread_var%ia) ) then

        if ( k == kmax ) then
          ! The bottom of the LES domain is above the top of the CLUBB
          ! domain.
          write(fstderr,*) "The lowest LES input level is above the top ",  &
                           "of the CLUBB model domain."
          stop "Error in CLUBB_levels_within_LES_domain"
        else
          ! Level k is the first CLUBB level below the LES domain.  Thus, the
          ! lowest CLUBB level within the LES domain has the index k + 1.
          k_lowest_input = k + 1
          exit
        endif

      elseif ( k == 1 ) then

        ! The bottom CLUBB level is within the LES domain.
        k_lowest_input = 1
        exit

      else   ! k > 1 and k <= kmax; level not yet found.

        ! Increment one more CLUBB vertical level down.
        k = k - 1

      endif
    enddo

    ! Find the highest CLUBB level that falls within the LES domain.
    k = 1
    do
      if ( CLUBB_grid(k) > fread_var%z(fread_var%iz) ) then

        if ( k == 1 ) then
          ! The top of the LES domain is below the bottom of the CLUBB
          ! domain.
          write(fstderr,*) "The highest LES input level is below the ",  &
                           "bottom of the CLUBB model domain."
          stop "Error in CLUBB_levels_within_LES_domain"
        else
          ! Level k is the first CLUBB level above the LES domain.  Thus, the
          ! highest CLUBB level within the LES domain has the index k - 1.
          k_highest_input = k - 1
          exit
        endif

      elseif ( k == kmax ) then

        ! The top CLUBB level is within the LES domain.
        k_highest_input = kmax 
        exit

      else   ! k < kmax and k >= 1; level not yet found.

        ! Increment one more CLUBB vertical level up.
        k = k + 1

      endif
    enddo

    return
  end subroutine CLUBB_levels_within_LES_domain

!===============================================================================
  pure subroutine LES_grid_to_CLUBB_grid( fread_var, CLUBB_grid, k,  &
                                          exact_lev_idx, lower_lev_idx,  &
                                          upper_lev_idx, l_lin_int )

    ! Description:
    !   Finds the level on the LES grid that is exactly even with the CLUBB
    !   grid level (either thermodynamic or momentum grid) that is input
    !   (level k).  Else, it finds the two LES levels that sandwich the CLUBB
    !   grid level that is input.

    !-----------------------------------------------------------------------

    use stat_file_module, only:  &
        stat_file  ! Variable type

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables.
    type(stat_file), intent(in) ::  &
      fread_var  ! Information about LES run.

    real( kind = core_rknd ), dimension(:), intent(in) ::  &
      CLUBB_grid ! Altitude of CLUBB grid levels
                 ! (either thermodynamic or momentum grid levels)  [m]

    integer, intent(in) ::  &
      k  ! Index of CLUBB vertical level that is being compared to.

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

      if ( abs( fread_var%z(j) - CLUBB_grid(k) ) < 1e-6_core_rknd ) then

        ! There is an LES level altitude at LES level j that is an exact
        ! match to the CLUBB level altitude at CLUBB grid level k.
        exact_lev_idx = j
        l_lin_int = .false.

      else if ( fread_var%z(j) < CLUBB_grid(k) ) then

        ! The LES level altitude at LES level j is lower than the CLUBB level
        ! altitude at CLUBB grid level k.
        lower_lev_idx = j

      else  ! fread_var%z(j) > CLUBB_grid(k)

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
