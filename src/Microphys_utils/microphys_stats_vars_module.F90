!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module microphys_stats_vars_module

! Description:
!   This module contains the derived type microphys_stats_vars_type, which is used to
!   feed output variables out of microphysics schemes

  use clubb_precision, only: &
    core_rknd

  implicit none

  private ! Set Default Scope

  public :: microphys_stats_vars_type, microphys_stats_alloc, microphys_put_var, &
            microphys_get_var, microphys_get_index, &
            microphys_stats_accumulate, microphys_stats_cleanup

  type microphys_stats_vars_type

    logical :: &
      l_allocated = .false. ! This is set to true when the structure is allocated

    integer :: &
      nz, &       ! Number of vertical levels
      num_vars, & ! Number of output variables from microphysics
      alloc_size  ! Size of allocated stats_indices and output_values arrays

    ! An array of statistics indices corresponding to the output variables
    integer, dimension(:), pointer :: &
      stats_indices

    ! Values of the output variables (over the vertical domain)
    real( kind = core_rknd ), dimension(:,:), pointer :: &
      output_values

  end type microphys_stats_vars_type

  contains

  !-----------------------------------------------------------------------
  subroutine microphys_stats_alloc( nz, num_vars, microphys_stats_vars )

  ! Description:
  !   Allocates a new microphysics stats type

  ! References:
  !   None
  !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz, &                  ! Number of vertical levels for statistics output
      num_vars               ! Number of statistics variables desired
                             ! (actual number output can be less)

    ! Input/Output Variables
    type(microphys_stats_vars_type), intent(out) :: &
      microphys_stats_vars   ! Unallocated microphys_stats_vars_type object

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    allocate( microphys_stats_vars%stats_indices(num_vars), &
              microphys_stats_vars%output_values(nz,num_vars) )

    microphys_stats_vars%nz = nz
    microphys_stats_vars%alloc_size = num_vars
    microphys_stats_vars%num_vars = 0 ! Since no variables have been put into structure

    microphys_stats_vars%l_allocated = .true.

    return
  end subroutine microphys_stats_alloc
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine microphys_put_var( var_index, value, microphys_stats_vars )

  ! Description:
  !   Places a variable in the microphys_stats_vars structure

  ! References:
  !   None
  !-----------------------------------------------------------------------

    use constants_clubb, only: &
      fstderr           ! Constant

    implicit none

    ! Input/Output Variables
    type(microphys_stats_vars_type), intent(inout) :: &
      microphys_stats_vars   ! The structure to place the variable into

    ! Input Variables
    integer, intent(in) :: &
      var_index              ! Statistics index of the variable

    real( kind = core_rknd ), dimension(microphys_stats_vars%nz), intent(in) :: &
      value

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    if ( microphys_stats_vars%num_vars == microphys_stats_vars%alloc_size ) then
      ! There is no more room in the structure. Do nothing (for now).
      write(fstderr,*) "Invalid allocation size"
      stop "Fatal error in microphys_put_var"
    end if

    if ( var_index > 0 ) then
      microphys_stats_vars%num_vars = microphys_stats_vars%num_vars + 1

      microphys_stats_vars%stats_indices(microphys_stats_vars%num_vars) = var_index
      microphys_stats_vars%output_values(:,microphys_stats_vars%num_vars) = value
    end if

    return
  end subroutine microphys_put_var
  !-----------------------------------------------------------------------

  subroutine microphys_stats_accumulate( microphys_stats_vars, l_stats_samp, grid )

  ! Description:
  !   Samples all variables stored in the microphys_stats_vars structure by
  !   using stat_update_var. The variables are sampled on the grid given as an
  !   argument to the subroutine. Therefore, there are no global variables used!

  ! References:
  !   None
  !-----------------------------------------------------------------------

    use stats_type_utilities, only: &
      stat_update_var ! Procedure

    use stats_type, only: &
      stats           ! Type

    implicit none

    ! Input Variables
    type(microphys_stats_vars_type), intent(in) :: &
      microphys_stats_vars

    logical, intent(in) :: &
      l_stats_samp              ! Are we sampling this timestep?

    ! Input/Output Variables
    type(stats), intent(inout) :: &
      grid                      ! Which grid to sample variables to

    ! Local Variables
    integer :: i ! Loop variable

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    if ( l_stats_samp ) then

      do i=1, microphys_stats_vars%num_vars

        call stat_update_var( microphys_stats_vars%stats_indices(i), &
                              microphys_stats_vars%output_values(:,i), grid )

      end do ! i=1, microphys_stats_vars%num_vars

    end if ! l_stats_samp

  end subroutine microphys_stats_accumulate
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  function microphys_get_index( stats_index, microphys_stats_vars ) result( var_index )

  ! Description:
  !   Get an index to the arrays in the microphys_stats_vars structure
  !   associated with the given stats index.

  ! References:
  !   Eric Raut
  !-----------------------------------------------------------------------
    implicit none

    ! Input Variables
    integer, intent(in) :: &
      stats_index             ! The stats index of the variable

    type(microphys_stats_vars_type), intent(in) :: &
      microphys_stats_vars    ! The statistics structure

    ! Output Variable
    integer :: &
      var_index               ! The index to the structure associated with
                              ! the given stats index

    ! Local Variables
    integer :: ivar
    logical :: l_found

  !-----------------------------------------------------------------------

    ! Initialize variables
    var_index = 0
    l_found = .false.

    !----- Begin Code -----
    do ivar=1, microphys_stats_vars%num_vars
      if ( microphys_stats_vars%stats_indices(ivar) == stats_index) then

        var_index = ivar
        l_found = .true.
        exit

      end if
    end do

    if ( .not. l_found ) then
      stop "Variable not found in microphys_get_index"
    end if

    return
  end function microphys_get_index
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  function microphys_get_var( stats_index, microphys_stats_vars ) result( stats_var )
  
  ! Description:
  !   Gets the specified statistics variable from the input structure.

  ! References:
  !   Eric Raut
  !-----------------------------------------------------------------------
    implicit none

    ! Input Variables
    integer, intent(in) :: &
      stats_index            ! The index of the variable

    type(microphys_stats_vars_type), intent(in) :: &
      microphys_stats_vars   ! The statistics structure

    ! Output Variable
    real( kind = core_rknd ), dimension(microphys_stats_vars%nz) :: &
      stats_var              ! The output statistics variable

    ! Local Variables
    integer :: ivar
    logical :: l_found

  !-----------------------------------------------------------------------

    ! Initialize variables
    stats_var = 0._core_rknd
    l_found = .false.

    !----- Begin Code -----
    do ivar=1, microphys_stats_vars%num_vars
      if ( microphys_stats_vars%stats_indices(ivar) == stats_index) then

        stats_var = microphys_stats_vars%output_values(:,ivar)
        l_found = .true.
        exit

      end if
    end do

    if ( .not. l_found ) then
      stop "Variable not found in microphys_get_var"
    end if

    return
  end function microphys_get_var
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine microphys_stats_cleanup( microphys_stats_vars )

  ! Description:
  !   Deallocates all (dynamic) memory associated with the
  !   microphys_stats_vars object

  ! References:
  !   None
  !-----------------------------------------------------------------------

    implicit none

    ! Input/Output Variables
    type(microphys_stats_vars_type), intent(inout) :: &
      microphys_stats_vars   ! Unallocated microphys_stats_vars_type object

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    deallocate( microphys_stats_vars%stats_indices, &
                microphys_stats_vars%output_values )

    microphys_stats_vars%l_allocated = .false.
    return
  end subroutine microphys_stats_cleanup
  !-----------------------------------------------------------------------

end module microphys_stats_vars_module
