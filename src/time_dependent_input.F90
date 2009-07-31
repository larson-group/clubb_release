!$Id$
module time_dependent_input
!
!  Description: This module is responsible for managing the reading in and
!  storage of time dependent information for a case.
!
!--------------------------------------------------------------------------------------------------

  use input_reader, only: &
    two_dim_read_var, &
    one_dim_read_var

  implicit none

  public :: initialize_t_dependent_input, finalize_t_dependent_input, time_select, &
            apply_time_dependent_forcings

  private :: initialize_t_dependent_forcings, &
             finalize_t_dependent_forcings,   & 
             initialize_t_dependent_surface,  &
             finalize_t_dependent_surface,    &
             read_to_grid

  integer, parameter :: nCols = 10 ! Number of columns in the input file

  real, public, target, allocatable, dimension(:) :: & ! Module variables used to describe 
    time_sfc_given, &                                  ! the surace over time.
    LH_given,       &
    SH_given,       &
    thlm_sfc_given, &
    rtm_sfc_given,  &
    psfc_given


  type(two_dim_read_var), private, dimension(nCols) :: &
    t_dependent_forcing_data ! Data structure that defines the change in input
                             ! files over time

  type(one_dim_read_var), private :: dimension_var ! Data structure that describes other 
                                                   ! dimension of the two_dim_read_var

  logical, public :: l_t_dependent ! Flag used to determine when
  !                                  time dependent information is read in.
  !                                  It is suggested that the flag be checked
  !                                  before using any of the variables stored
  !                                  in the module.


  ! File path constants
  character(len=*), private, parameter :: input_path = "../input/case_setups/"

  character(len=*), private, parameter :: forcings_path = "_forcings.in"

  character(len=*), private, parameter :: surface_path = "_surface.in"

  private

  contains

  !-------------------------------------------------------------------------------------------------
  subroutine initialize_t_dependent_input( iunit, runtype, grid_size, grid, p_in_Pa )
    !
    !  Description: This subroutine reads in time dependent information about a
    !  case that is stored inside the module.
    !
    !-----------------------------------------------------------------------------------------------

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: iunit ! File I/O

    character(len=*), intent(in) :: runtype ! Runtype

    integer, intent(in) :: grid_size ! Size of the model grid

    real, dimension(grid_size), intent(in) :: grid ! Model grid

    real, dimension(grid_size), intent(in) :: p_in_Pa ! Pressure[Pa]

    ! Begin Code

    call initialize_t_dependent_forcings &
                   ( iunit, input_path//trim(runtype)//forcings_path, grid_size, grid, p_in_Pa )

    call initialize_t_dependent_surface &
                   ( iunit, input_path//trim(runtype)//surface_path )

  end subroutine initialize_t_dependent_input

  !------------------------------------------------------------------------------
  subroutine finalize_t_dependent_input()
    !
    ! Description: This subroutine frees memory stored after initilizing the
    ! time dependent data of this module.
    !
    !-----------------------------------------------------------------------------

    implicit none

    ! Begin Code

    call finalize_t_dependent_forcings()
    call finalize_t_dependent_surface()

  end subroutine finalize_t_dependent_input

  !------------------------------------------------------------------------------
  subroutine initialize_t_dependent_surface( iunit, input_file )
    !
    !  Description: This subroutine reads in a file that details time dependent
    !  input values that vary in one dimension.
    !-----------------------------------------------------------------------------

    use input_reader, only: read_one_dim_file, one_dim_read_var, &
                            fill_blanks_one_dim_vars, read_x_profile

    use input_names, only: &
      time_name,     &
      thetal_name,   &
      rt_name,       &
      LH_name,       &
      SH_name,       &
      pressure_name

    implicit none

    ! Constants
    integer, parameter :: nCols = 6

    ! Input Variable(s)
    integer, intent(in) :: iunit ! File I/O unit

    character(len=*), intent(in) :: input_file ! Path to surface.in file

    ! Local Variable(s)

    type(one_dim_read_var), dimension(nCols) :: retVars

    integer dim_size

    ! Begin Code

    call read_one_dim_file( iunit, nCols, input_file, retVars )

    call fill_blanks_one_dim_vars( nCols, retVars )

    dim_size = size( retVars(1)%values )

    allocate( time_sfc_given( 1:dim_size ) )

    time_sfc_given = read_x_profile( nCols, dim_size, time_name, retVars )

    allocate( LH_given( 1:dim_size ) )

    LH_given = read_x_profile( nCols, dim_size, LH_name, retVars )

    allocate( SH_given( 1:dim_size ) )

    SH_given = read_x_profile( nCols, dim_size, SH_name, retVars )

    allocate( thlm_sfc_given( 1:dim_size ) )

    thlm_sfc_given = read_x_profile( nCols, dim_size, thetal_name, retVars )

    allocate( rtm_sfc_given( 1:dim_size ) )

    rtm_sfc_given = read_x_profile( nCols, dim_size, rt_name, retVars )

    allocate( psfc_given( 1:dim_size ) )

    psfc_given = read_x_profile( nCols, dim_size, pressure_name, retVars )

  end subroutine initialize_t_dependent_surface

  !-------------------------------------------------------------------------------------
  subroutine initialize_t_dependent_forcings( iunit, input_file, grid_size, grid, p_in_Pa )
    !
    !  Description: This subroutine reads in a file that details time dependent
    !  input values that vary in two dimensions.
    !
    !-------------------------------------------------------------------------------------

    use input_reader, only: read_two_dim_file, two_dim_read_var, fill_blanks_two_dim_vars

    use input_names, only: &
      z_name, &
      pressure_name

    implicit none


    ! Input Variable(s)
    integer, intent(in) :: iunit ! File I/O

    character(len=*), intent(in) :: input_file ! Path to the input file

    integer, intent(in) :: grid_size  ! Size of Model Grid [-]

    real, intent(in), dimension(grid_size) :: grid ! Altitudes of Grid [m]

    real, intent(in), dimension(grid_size) :: p_in_Pa ! Pressure [Pa]

    ! Local Variables

    integer :: i, n_f_grid_z, n_f_grid_t

    type(two_dim_read_var), dimension(nCols) :: t_dependent_forcing_data_f_grid

    ! Begin Code


    ! Read in the forcing data from the input file
    call read_two_dim_file( iunit, nCols, input_file, &
                            t_dependent_forcing_data_f_grid, dimension_var )

    n_f_grid_z = size( t_dependent_forcing_data_f_grid(1)%values, 1 )

    n_f_grid_t = size( dimension_var%values )

    ! Fill in blanks with linear interpolation. Whole profiles of -999.9 will
    ! remain that way thus marking them blank.
    call fill_blanks_two_dim_vars( nCols, dimension_var, t_dependent_forcing_data_f_grid )

    do i=1, nCols
      allocate( t_dependent_forcing_data(i)%values(1:grid_size,1:n_f_grid_t) )
    end do

    select case( t_dependent_forcing_data_f_grid(1)%name )
    case( z_name )

      t_dependent_forcing_data(1)%name = z_name
      t_dependent_forcing_data(1)%values(:,1:n_f_grid_t) = spread(grid,2,n_f_grid_t)

    case( pressure_name )

      t_dependent_forcing_data(1)%name = pressure_name
      t_dependent_forcing_data(1)%values(:,1:n_f_grid_t) = spread(-p_in_Pa, 2, n_f_grid_t )
      t_dependent_forcing_data_f_grid(1)%values = -t_dependent_forcing_data_f_grid(1)%values

    case default
      stop "Incompatible grid type in first element of t_dependent_forcings."
    end select

    ! Interpolate the time dependent input data to the appropriate grid.

    do i=2, nCols

      t_dependent_forcing_data(i)%name = t_dependent_forcing_data_f_grid(i)%name
      t_dependent_forcing_data(i)%values = read_to_grid( nCols, n_f_grid_z, n_f_grid_t, &
                         grid_size,t_dependent_forcing_data(1)%values(:,1), &
                         t_dependent_forcing_data_f_grid, t_dependent_forcing_data(i)%name )

    end do


  end subroutine initialize_t_dependent_forcings

  !----------------------------------------------------------
  subroutine finalize_t_dependent_forcings()
    !
    !   Description: Clears memory initialized in initialize_t_dependent_forcings.
    !   This should be called at the end of the model
    !----------------------------------------------------------

    implicit none

    integer i

    ! Begin Code

    do i=1, nCols
      deallocate( t_dependent_forcing_data(i)%values )
    end do

  end subroutine finalize_t_dependent_forcings
  !-------------------------------------------------------------------------------------------------
  subroutine finalize_t_dependent_surface( )
    !
    !  Description: Clears memory initialized in initialize_t_dependent_surface.
    !  This should be called at the end of the model.
    !
    !-----------------------------------------------------------------------------------------------

    implicit none

    ! Begin Code

    deallocate( time_sfc_given )
    deallocate( LH_given )
    deallocate( SH_given )
    deallocate( thlm_sfc_given )
    deallocate( rtm_sfc_given )
    deallocate( psfc_given )

  end subroutine finalize_t_dependent_surface

  !------------------------------------------------------------------------------------------------
  function read_to_grid( ntwo_dim_vars, dim_size, other_dim_size, &
                         grid_size, grid, two_dim_vars, target_name ) result(var)
    !
    !  Description: This is a helper function for doing the translation from the
    !  forcing grid to the model grid.
    !
    !----------------------------------------------------------------------------------------------

    use input_reader, only: read_x_table, two_dim_read_var

    use interpolation, only: zlinterp_fnc

    implicit none

    integer, intent(in) :: &
      ntwo_dim_vars, &
      dim_size, &
      other_dim_size, &
      grid_size

    real, dimension(grid_size), intent(in) :: &
      grid

    type(two_dim_read_var), dimension(ntwo_dim_vars), intent(in) :: &
      two_dim_vars

    character(len=*), intent(in) :: &
      target_name

    real, dimension(dim_size, other_dim_size) :: temp_var

    real, dimension(grid_size, other_dim_size) :: var

    integer i

    ! Begin Code

    temp_var = read_x_table( ntwo_dim_vars,  dim_size, other_dim_size, target_name, two_dim_vars )

    do i=1, other_dim_size
      var(:,i) = zlinterp_fnc( grid_size, dim_size, grid, &
                                    two_dim_vars(1)%values(:,i), temp_var(:,i) )
    end do

    return

  end function read_to_grid

  !-------------------------------------------------------------------------------------------------
  subroutine apply_time_dependent_forcings( time, grid_size, rtm, rho, exner,  &
    thlm_f, rtm_f, um_ref, vm_ref, um_f, vm_f, wm_zt, wm_zm,  ug, vg, &
    sclrm_forcing, edsclrm_forcing )
    !
    !  Description: This subroutine converts the time dependent information stored in
    !  memory (time_dependent_forcing_data) into the format used by CLUBB.
    !
    !-----------------------------------------------------------------------------------------------

    use error_code, only: &
      clubb_debug ! Procedure(s)

    use constants, only: &
      grav, & ! Variable(s)
      sec_per_hr, &
      pascal_per_mb, &
      fstderr

    use interpolation, only: &
      factor_interp ! Procedure(s)

    use input_names, only: &
      z_name, & ! Variable(s)
      pressure_name, &
      temperature_f_name, &
      rt_f_name,&
      sp_humidity_f_name, &
      thetal_f_name, &
      theta_f_name, &
      wm_name, &
      omega_name, &
      um_ref_name, &
      vm_ref_name, &
      um_f_name, &
      vm_f_name, &
      ug_name,&
      vg_name, &
      omega_mb_hr_name

    use grid_class, only : zt2zm ! Procedure(s)

    use stats_precision, only: time_precision ! Variable(s)

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use array_index, only: iisclr_rt, iisclr_thl, & ! Variable(s)
                           iiedsclr_rt, iiedsclr_thl

    implicit none

    ! Input Variable(s)

    real(kind=time_precision), intent(in) :: time ! Model Time [s]

    integer, intent(in) :: grid_size ! Size of the model grid

    real, dimension(grid_size), intent(in) :: &
      exner,   & ! Exner Function                             [-]
      rho,     & ! Air Density                                [kg/m^3]
      rtm        ! Total Water Mixing Ratio                   [kg/kg]

    ! Output Variable(s)

    real, dimension(grid_size), intent(out) :: &
      thlm_f, & ! Potential Temperature forcing     [K/s]
      rtm_f,  & ! Total Water Mixing Ration forcing [kg/kg/s]
      um_ref, & ! um reference                      [m/s]
      vm_ref, & ! vm reference                      [m/s]
      um_f,   & ! um tendency                       [m/s/s]
      vm_f,   & ! vm tendency                       [m/s/s]
      wm_zt,  & ! subsidence on zt grid             [m/s]
      wm_zm,  & ! subsidense on zm grid             [m/s]
      ug,     & ! u geostrophic wind                [m/s]
      vg        ! v geostrophic wind                [m/s]

    real, dimension( grid_size, sclr_dim ), intent(out) :: &
      sclrm_forcing ! Scalar forcing [-]

    real, dimension( grid_size, edsclr_dim ), intent(out) :: &
      edsclrm_forcing ! Edscalar forcing [-]

    ! Local Variable(s)
    integer :: i, j, i1, i2

    real, dimension(grid_size) :: temp_array

    real time_frac

    ! Begin Code

    time_frac = -1.0 ! Default initialization

    call time_select( time, size(dimension_var%values), dimension_var%values, i1, i2 )

    ! Determine time interpolation factor
    time_frac =  real( ( time-dimension_var%values(i1) )/ &
                      ( dimension_var%values(i2) - dimension_var%values(i1) ) )

    if( time_frac == -1.0 ) then
      call clubb_debug(1,"times are not sorted in forcing")
    endif

    ! Parse the values in t_dependent_forcing_data for CLUBB compatible forcing
    ! data.
    do i=2, nCols

      temp_array = factor_interp( time_frac, t_dependent_forcing_data(i)%values(:,i2), &
                                         t_dependent_forcing_data(i)%values(:,i1) )

      ! Check to see if temp_array is an actual profile or a dummy profile
      ! If it is a dummy profile we dont want it to apply itself as it may
      ! overwrite legitimate information from another source.
      if( .not. any( temp_array == -999.9 ) ) then
        select case (t_dependent_forcing_data(i)%name)
        case(temperature_f_name, theta_f_name, thetal_f_name)

          select case(t_dependent_forcing_data(i)%name)
          case(temperature_f_name)

            thlm_f = temp_array / exner

          case(theta_f_name)

            thlm_f = temp_array ! I am not sure on the conversion of this

          case(thetal_f_name)

            thlm_f = temp_array

          end select

          if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_f
          if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_f

        case(rt_f_name, sp_humidity_f_name)

          select case(t_dependent_forcing_data(i)%name)
          case(sp_humidity_f_name)

            rtm_f = temp_array * ( 1. + rtm )**2

          case(rt_f_name )

            rtm_f = temp_array

          end select

          if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_f
          if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_f


        case(um_ref_name)

          um_ref = temp_array

        case(vm_ref_name)

          vm_ref = temp_array

        case(um_f_name)

          um_f = temp_array

        case(vm_f_name)

          vm_f = temp_array

        case(wm_name, omega_name, omega_mb_hr_name)

          select case(t_dependent_forcing_data(i)%name)
          case(wm_name)

            wm_zt = temp_array

          case(omega_name)

            do j=2,grid_size
              wm_zt(j) = - temp_array(j) / (grav * rho(j))
            end do
            wm_zt(1) = 0.0

          case(omega_mb_hr_name)

            do j=2,grid_size

              temp_array(j) = temp_array(j) * pascal_per_mb / real( sec_per_hr )

              wm_zt(j) = - temp_array(j) / (grav * rho(j))

            end do

            wm_zt(1) = 0.0

          end select

          wm_zm = zt2zm( wm_zt )

        case(ug_name)

          ug = temp_array
          ug(1) = ug(2)

        case(vg_name)

          vg = temp_array
          vg(1) = vg(2)

        case default

          write(fstderr, *) "Incompatable forcing type: "//t_dependent_forcing_data(i)%name
          stop

        end select

      end if 

    end do ! 1 .. nCols

  end subroutine apply_time_dependent_forcings

  !------------------------------------------------------------------------------------------------
  subroutine time_select( time, nvar, time_array, left_time, right_time )
    !
    !   Description: This subroutine determines which indexes of the given
    !                time_array should be used when interpolating a value
    !                at the specified time.
    !
    !
    !----------------------------------------------------------------------------------------------

    use stats_precision, only: time_precision ! Variable(s)

    implicit none

    ! Input Variable(s)

    integer, intent(in) :: nvar                     ! Number of array elements [-]

    real(kind=time_precision), intent(in) :: time   ! Target time              [s]

    real, dimension(nvar), intent(in) :: time_array ! Array of times           [s]

    ! Output Variable(s)

    integer, intent(out) :: &
      right_time, &  ! Index of a time later than the target time [-]
      left_time      ! Index of time before the target time       [-]


    ! Local Variable(s)

    integer :: k

    ! Begin Code

    if( time <= time_array(1)) then

      left_time = 1
      right_time = 2

    else if ( time >= time_array(nvar) ) then

      left_time = nvar
      right_time = nvar - 1

    else

      do k=1,nvar-1

        if ((time > time_array(k)) .and. &
         (time <= time_array(k+1))) then

          left_time = k
          right_time = k+1

        end if

      end do

    endif

    return

  end subroutine time_select

end module time_dependent_input
