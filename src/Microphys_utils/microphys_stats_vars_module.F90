! $Id$
module microphys_stats_vars_module

! Description:
!   This module contains the derived type microphys_stats_vars_type, which is used to
!   feed output variables out of microphysics schemes

  use clubb_precision, only: &
    core_rknd

  implicit none

  private ! Set Default Scope

  public :: microphys_stats_vars_type, microphys_stats_alloc

  type microphys_stats_vars_type

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
    type(microphys_stats_vars_type), intent(inout) :: &
      microphys_stats_vars   ! Unallocated microphys_stats_vars_type object

    ! Local Variables
    integer, dimension(:), allocatable, target :: &
      stats_indices_alloc

    real( kind = core_rknd ), dimension(:,:), allocatable, target :: &
      output_values_alloc

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    allocate( stats_indices_alloc(num_vars), output_values_alloc(nz,num_vars) )

    microphys_stats_vars%nz = nz
    microphys_stats_vars%alloc_size = num_vars
    microphys_stats_vars%num_vars = 0 ! Since no variables have been put into structure

    microphys_stats_vars%stats_indices => stats_indices_alloc
    microphys_stats_vars%output_values => output_values_alloc

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
      ! There is no more room in the structure. Do nothing.
      return
    end if

    microphys_stats_vars%num_vars = microphys_stats_vars%num_vars + 1

    microphys_stats_vars%stats_indices(microphys_stats_vars%num_vars) = var_index
    microphys_stats_vars%output_values(:,microphys_stats_vars%num_vars) = value

    return
  end subroutine microphys_put_var
  !-----------------------------------------------------------------------

end module microphys_stats_vars_module
