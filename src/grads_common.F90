! $Id$
module grads_common

  implicit none

  public :: grads_average, grads_average_interval, &
    grads_num_vertical_levels, grads_vertical_levels

  private ! Default Scope

  contains

!----------------------------------------------------------------------
  function grads_average( filename, out_nz, &
                          t1, t2, out_heights, variable_name,  & 
                          npower, l_error )
! Description:
!   Average a GrADS file variable over the interval t1 to t2

! References:
!   None
!----------------------------------------------------------------------

    use constants, only: fstderr ! Variable(s)

    use stat_file_module, only: stat_file ! Type(s)

    use inputfile_class, only: open_grads_read, get_var,  & ! Procedures
                         close_grads_read
    use inputfields, only: &
      CLUBB_levels_within_LES_domain, &
      LES_grid_to_CLUBB_grid, &
      lin_ext_zt_bottom

    use interpolation, only: &
      lin_int

    implicit none

    ! External
    intrinsic :: transfer

    ! Constant paramters
    integer, parameter :: &
      iunit = 10

    ! Input Variables
    character(len=*), intent(in) ::  & 
      filename ! Name of the file

    integer, intent(in) ::  & 
      out_nz,  & ! Number of vertial levels in the GrADS file.
      t1,      & ! Beginning timestep to look at
      t2         ! Ending timestep to look at

    real, dimension(out_nz), intent(in) :: &
      out_heights ! Heights of the output grid [m]

    character(len=*), intent(in) :: & 
      variable_name ! Name of the variable to read in

    integer, intent(in) :: & 
      npower ! Exponent operator, must be 1 or 2

    ! Output Variable
    logical, intent(out) :: & 
      l_error ! error status for this function

    ! Return Variable for function
    real, dimension(out_nz) :: grads_average

    ! Local Variables
    type (stat_file) :: faverage ! Data file derived type

    real, allocatable, dimension(:) :: file_variable ! Temporary variable

    real, dimension(out_nz) ::  interp_variable ! Temporary variable

    integer :: & 
      t, &  ! Timestep loop index
      k     ! Vertical loop index
    integer :: file_nz
    logical :: l_interpolate

    logical, dimension(out_nz) :: l_lin_int

    integer, dimension(out_nz) :: &
      exact_lev_idx, &
      lower_lev_idx, &
      upper_lev_idx

    integer :: & 
      num_timesteps,  & ! steps between t1 and t2
      k_lowest_input, &
      k_highest_input, &
      nanbits

    ! Set nanbits to NaN
    data nanbits /Z"7F800000"/

!-----------------------------------------------------------------------
    l_error = .false.

    ! Initialize variables
    num_timesteps = ( t2 - t1 ) + 1
    grads_average = 0.0

    ! Open grads file
    call open_grads_read( iunit, filename, faverage )

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
      call get_var( faverage, variable_name, t,  & 
                    file_variable(1:file_nz), l_error )

      if ( l_error ) then
        write(fstderr,*) "grads_average: get_var failed for "  & 
          //trim( variable_name )//" in "//trim( filename )//" at time=", t
        return
      end if

      ! Interpolate as needed using the inputfields code
      if ( l_interpolate ) then

        do k = k_lowest_input, k_highest_input, 1
          if ( l_lin_int(k) ) then
            interp_variable(k) = lin_int( out_heights(k), &
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

          ! Set undefined points to NaN
          interp_variable(k) = transfer( nanbits, interp_variable(k) )
        end do

        do k = k_highest_input+1, out_nz, 1
          ! Do the the dum-dum thing and set points above output domain 
          ! to be constant with height
!         interp_variable(k) = interp_variable(k-1)

          ! Set undefined points to NaN
          interp_variable(k) = transfer( nanbits, interp_variable(k) )
        end do

      else
        interp_variable(1:out_nz) = file_variable(1:out_nz)

      end if


      ! Apply an exponent (for the tuner)
      if ( npower /= 1 ) then
        grads_average(1:out_nz) = grads_average(1:out_nz) + interp_variable(1:out_nz)**npower

      else
        grads_average(1:out_nz) = grads_average(1:out_nz) + interp_variable(1:out_nz)

      end if

    end do ! t = t1, t2

    ! Close the GrADS file
    call close_grads_read( faverage )

    ! Take average over num_timesteps
    grads_average(1:out_nz) = grads_average(1:out_nz) / real( num_timesteps )

    return
  end function grads_average

!-------------------------------------------------------------------------
  function grads_average_interval & 
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
    use constants, only: fstderr ! Variable(s)

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

    real, dimension(nz), intent(in) :: out_heights

    integer, intent(in) ::  & 
      npower ! exponent applied to data retrieved from the file (1 or 2)

    ! Output Variables
    logical, intent(out) ::  & 
      l_error ! status of this function

    ! Return Variables
    real, dimension(nz) ::  & 
      grads_average_interval

    ! Local Variables
    real, dimension(nz) :: grads_temp

    integer ::  & 
      i,       & ! Loop variable 
      tdim,    & ! Dimension to read over for t variable
      divisor

!-------------------------------------------------------------------------
    l_error = .false.

    ! Sanity check
    if ( size( t ) > tmax .or. size( t ) < 2 ) then
      write(unit=fstderr,fmt=*)  & 
        "grads_average_interval: Invalid time interval"
      l_error = .true.
      return
    end if

    tdim = -1
    do i = 1, size( t )
      if ( t( i ) == 0 ) exit
      tdim = i
    end do

    grads_average_interval & 
    = grads_average & 
      ( filename, nz, t(1), t(2), out_heights, variable_name, npower, l_error )  & 
      * ( t(2) - t(1) )

    divisor = t(2) - t(1)

    if ( l_error ) return

    do i=3, tdim, 2
      grads_temp = grads_average & 
                   ( filename, nz, t(i), t(i+1), out_heights,  & 
                     variable_name, npower, l_error )
      grads_average_interval  & 
      = grads_average_interval + grads_temp * ( t(i+1) - t(i) )
      divisor = divisor + ( t(i+1) - t(i) )
    end do

    grads_average_interval(1:nz)  = grads_average_interval(1:nz) / real( divisor )

    return
  end function grads_average_interval

!-------------------------------------------------------------------------
  integer function grads_num_vertical_levels( filename )

!       Description:
!       Returns a integer for the number of vertical levels in a file

!       References:
!       None

!       Notes:
!       Somewhat inefficient, since it reads in the entire .ctl file to
!       determine the number of levels
!-------------------------------------------------------------------------

    use stat_file_module, only: stat_file ! Type(s)

    use inputfile_class, only: open_grads_read, close_grads_read ! Procedure(s)

    implicit none

    ! Input Variables
    character(len=*), intent(in) ::  & 
      filename ! File name

    ! Local Variables
    type (stat_file) :: fz            ! Data file

    ! Read in the control file
    call open_grads_read( 10, filename, fz )

    ! Set return variable
    grads_num_vertical_levels = fz%iz

    ! Close file
    call close_grads_read( fz )

    return
  end function grads_num_vertical_levels

!-------------------------------------------------------------------------
  function grads_vertical_levels( filename, nz )

    use stat_file_module, only: stat_file ! Type(s)

    use inputfile_class, only: open_grads_read, close_grads_read ! Procedure(s)

    implicit none

    ! Input Variables
    character(len=*), intent(in) ::  & 
      filename ! File name

    integer, intent(in) :: &
      nz ! Number of vertical levels

    ! Output Variables
    real, dimension(nz) :: grads_vertical_levels

    ! Local Variables
    type (stat_file) :: fz  ! Data file

    ! Read in the control file
    call open_grads_read( 10, filename, fz )

    ! Set return variable
    grads_vertical_levels(1:nz) = fz%z(1:nz)

    ! Close file
    call close_grads_read( fz )

    return
  end function grads_vertical_levels
!-------------------------------------------------------------------------
end module grads_common
