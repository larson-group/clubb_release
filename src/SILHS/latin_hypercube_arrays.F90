!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module latin_hypercube_arrays

  use clubb_precision, only: &
    dp, & ! double precision
    core_rknd

  implicit none

  public :: cleanup_latin_hypercube_arrays

  private

  integer, allocatable, dimension(:,:), public :: & 
    height_time_matrix ! matrix of rand ints

!$omp threadprivate(height_time_matrix)

  contains

  !-----------------------------------------------------------------------------
  subroutine cleanup_latin_hypercube_arrays( )

    ! Description:
    !   De-allocate latin hypercube arrays
    ! References:
    !   None
    !---------------------------------------------------------------------------
    implicit none

    ! External
    intrinsic :: allocated

    ! ---- Begin Code ----

    if ( allocated( height_time_matrix ) ) then
      deallocate( height_time_matrix )
    end if

    return
  end subroutine cleanup_latin_hypercube_arrays

end module latin_hypercube_arrays
