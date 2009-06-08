!$Id$
module time_dependant_input
!
!  Description: This module is responsible for managing the reading in and
!  storage of time dependant information for a case.
!
!--------------------------------------------------------------------------------------------------


  implicit none

  public :: initialize_t_dependant_input, finalize_t_dependant_input, time_select

  private :: initialize_t_dependant_forcings, &
             finalize_t_dependant_forcings, & 
             initialize_t_dependant_surface, &
             finalize_t_dependant_surface, &
             read_to_grid

  ! array(altitude,time)
  real, public, target, allocatable, dimension(:,:) :: &
    thlm_f_given, &
    rtm_f_given, & 
    um_given, &
    vm_given, &
    um_f_given, &
    vm_f_given, &
    wm_given, &
    ug_given, &
    vg_given

  real, public, target, allocatable, dimension(:) :: &
    time_f_given, &   
    time_sfc_given, &
    LH_given, &
    SH_given, &
    thlm_sfc_given, &
    rtm_sfc_given, &
    psfc_given


  logical, public :: l_t_dependant ! Flag used to determine when
  !                                     time dependant information is read in.
  !                                     It is suggested that the flag be checked
  !                                     before using any of the variables stored
  !                                     in the module.


  ! File path constants
  character(len=*), private, parameter :: input_path = "../input/case_setups/"

  character(len=*), private, parameter :: forcings_path = "_forcings.in"

  character(len=*), private, parameter :: surface_path = "_surface.in"

  private

  contains

  !-------------------------------------------------------------------------------------------------
  subroutine initialize_t_dependant_input( iunit, runtype, grid_size, grid )
    !
    !  Description: This subroutine reads in time dependant information about a
    !  case that is stored inside the module.
    !
    !-----------------------------------------------------------------------------------------------

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: iunit ! File I/O

    character(len=*), intent(in) :: runtype ! Runtype

    integer, intent(in) :: grid_size ! Size of the model grid

    real, dimension(grid_size), intent(in) :: grid ! Model grid

    ! Begin Code

    call initialize_t_dependant_forcings &
                   ( iunit, input_path//trim(runtype)//forcings_path, grid_size, grid )

    call initialize_t_dependant_surface &
                   ( iunit, input_path//trim(runtype)//surface_path )

  end subroutine initialize_t_dependant_input

  !------------------------------------------------------------------------------
  subroutine finalize_t_dependant_input()
    !
    ! Description: This subroutine frees memory stored after initilizing the
    ! time dependant data of this module.
    !
    !-----------------------------------------------------------------------------

    implicit none

    ! Begin Code

    call finalize_t_dependant_forcings()
    call finalize_t_dependant_surface()

  end subroutine

  !------------------------------------------------------------------------------
  subroutine initialize_t_dependant_surface( iunit, input_file )
    !
    !  Description: This subroutine reads in a file that details time dependant
    !  input values that vary in one dimension.
    !-----------------------------------------------------------------------------

    use input_reader, only: read_one_dim_file, one_dim_read_var, &
                            fill_blanks_one_dim_vars

    use sounding, only: read_x_profile

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

    time_sfc_given = read_x_profile( nCols, dim_size, 'Time[s]', retVars )

    allocate( LH_given( 1:dim_size ) )

    LH_given = read_x_profile( nCols, dim_size, 'LH[W\m^2]', retVars )

    allocate( SH_given( 1:dim_size ) )

    SH_given = read_x_profile( nCols, dim_size, 'SH[W\m^2]', retVars )

    allocate( thlm_sfc_given( 1:dim_size ) )

    thlm_sfc_given = read_x_profile( nCols, dim_size, 'thlm[K]', retVars )

    allocate( rtm_sfc_given( 1:dim_size ) )

    rtm_sfc_given = read_x_profile( nCols, dim_size, 'rtm[kg\kg]', retVars )

    allocate( psfc_given( 1:dim_size ) )

    psfc_given = read_x_profile( nCols, dim_size, 'Press[Pa]', retVars )

  end subroutine

  !-------------------------------------------------------------------------------------
  subroutine initialize_t_dependant_forcings( iunit, input_file, grid_size, grid )
    !
    !  Description: This subroutine reads in a file that details time dependant
    !  input values that vary in two dimensions.
    !
    !-------------------------------------------------------------------------------------

    use input_reader, only: read_two_dim_file, two_dim_read_var, one_dim_read_var, &
                            fill_blanks_two_dim_vars

    implicit none

    integer, parameter :: nCols = 10 ! Number of columns in the input file

    ! Input Variable(s)
    integer, intent(in) :: iunit ! File I/O

    character(len=*), intent(in) :: input_file ! Path to the input file

    integer, intent(in) :: grid_size ! Size of the model grid

    real, dimension(grid_size), intent(in) :: grid ! Model grid

    ! Local Variables
    type(two_dim_read_var), dimension(nCols) :: retVars

    type(one_dim_read_var) :: dimension_var

    integer dim_size

    integer other_dim_size

    ! Begin Code

    call read_two_dim_file( iunit, nCols, input_file, retVars, dimension_var )

    call fill_blanks_two_dim_vars( nCols, dimension_var, retVars )

    dim_size = size(retVars(1)%values,1)

    other_dim_size = size(dimension_var%values)

    allocate( thlm_f_given( grid_size, other_dim_size ) )

    thlm_f_given = read_to_grid( nCols, dim_size, other_dim_size, &
                                 grid_size, grid, retVars, 'thlm_f[K\s]' )

    allocate( rtm_f_given( grid_size, other_dim_size ) )

    rtm_f_given = read_to_grid( nCols, dim_size, other_dim_size, &
                                grid_size, grid, retVars, 'rtm_f[kg\kg\s]' )


    allocate( um_given( grid_size, other_dim_size ) )

    um_given = read_to_grid( nCols, dim_size, other_dim_size, &
                             grid_size, grid, retVars, 'um[m\s]' )

    allocate( vm_given( grid_size, other_dim_size ) )

    vm_given = read_to_grid( nCols, dim_size, other_dim_size, &
                              grid_size, grid, retVars, 'vm[m\s]' )


    allocate( um_f_given( grid_size, other_dim_size ) )

    um_f_given = read_to_grid( nCols, dim_size, other_dim_size, &
                               grid_size, grid, retVars, 'um_f[m\s^2]' )

    allocate( vm_f_given( grid_size, other_dim_size ) )

    vm_f_given = read_to_grid( nCols, dim_size, other_dim_size, &
                               grid_size, grid, retVars, 'vm_f[m\s^2]' )

    allocate( wm_given( grid_size, other_dim_size ) )

    wm_given = read_to_grid( nCols, dim_size, other_dim_size, &
                             grid_size, grid, retVars, 'wm[m\s]' )

    allocate( ug_given( grid_size, other_dim_size ) )

    ug_given = read_to_grid( nCols, dim_size, other_dim_size, &
                             grid_size, grid, retVars, 'ug[m\s]' )

    allocate( vg_given( grid_size, other_dim_size ) )

    vg_given = read_to_grid( nCols, dim_size, other_dim_size, &
                             grid_size, grid, retVars, 'vg[m\s]' )


    allocate( time_f_given(other_dim_size ) )

    time_f_given = dimension_var%values

    return

  end subroutine initialize_t_dependant_forcings

  !----------------------------------------------------------
  subroutine finalize_t_dependant_forcings()
    !
    !   Description: Clears memory initialized in initialize_t_dependant_forcings.
    !   This should be called at the end of the model
    !----------------------------------------------------------

    implicit none

    ! Begin Code

    deallocate( thlm_f_given )
    deallocate( rtm_f_given )
    deallocate( um_given )
    deallocate( vm_given )
    deallocate( um_f_given )
    deallocate( vm_f_given )
    deallocate( wm_given )
    deallocate( ug_given )
    deallocate( vg_given )


  end subroutine finalize_t_dependant_forcings
  !-------------------------------------------------------------------------------------------------
  subroutine finalize_t_dependant_surface( )
    !
    !  Description: Clears memory initialized in initialize_t_dependant_surface.
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

  end subroutine finalize_t_dependant_surface

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

end module time_dependant_input
