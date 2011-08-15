!-------------------------------------------------------------------------
! $Id$
module advance_helper_module

! Description:
!   This module contains helper methods for the advance_* modules.
!------------------------------------------------------------------------

  implicit none

  public :: set_boundary_conditions

  private ! Set Default Scope

  contains

  !---------------------------------------------------------------------------
  subroutine set_boundary_conditions( diag_index, low_bound, high_bound, lhs, &
                                      diag_index2, low_bound2, high_bound2 )

  ! Description:
  !   Sets the boundary conditions for a LAPACK matrix.
  !
  ! References:
  !   none
  !---------------------------------------------------------------------------

    implicit none

    integer, intent(in) :: &
      diag_index, low_bound, high_bound ! boundary indexes for the first variable

    integer, intent(in), optional :: &
      diag_index2, low_bound2, high_bound2 ! boundary indexes for the second variable

    real, dimension(:,:), intent(inout) :: &
      lhs ! left hand side of the LAPACK matrix equation

    ! --------------------- BEGIN CODE ----------------------

    if( ( present(low_bound2) .or. present(high_bound2) ) .and. &
         ( .not. present(diag_index2) ) ) then

      stop "Boundary index provided without diag_index."

    end if

    ! Set the lower boundaries for the first variable
    lhs(:,low_bound) = 0.0
    lhs(diag_index,low_bound) = 1.0

    ! Set the upper boundaries for the first variable
    lhs(:,high_bound) = 0.0
    lhs(diag_index,high_bound) = 1.0

    ! Set the lower boundaries for the second variable, if it is provided
    if( present(low_bound2) ) then

      lhs(:,low_bound2) = 0.0
      lhs(diag_index2,low_bound2) = 1.0

    end if

    ! Set the upper boundaries for the second variable, if it is provided
    if( present(high_bound2) ) then

      lhs(:,high_bound2) = 0.0
      lhs(diag_index2,high_bound2) = 1.0

    end if

  end subroutine set_boundary_conditions

end module advance_helper_module
