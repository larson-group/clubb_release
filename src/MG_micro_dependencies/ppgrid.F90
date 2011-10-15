! $Id$
module ppgrid
!
! Dummy module for importing variables into morrison-gettelman microphysics
!---------------------------------------------------------------------------------------------------

  implicit none

  private

  public :: pcols, pver, pverp, init_ppgrid

  integer :: &
    pcols,   &  ! number of points in the “i” direction.
    pver,    &  ! number of points in the vertical
    pverp       ! pver + 1
    
  parameter(pcols = 1) ! CLUBB is a single column model.
    
  contains
  
!================================================================================================
  subroutine init_ppgrid( nz )
    !
    !  Description: Initialize ppgrid variables for being imported into MG microphysics.
    !---------------------------------------------------------------------------------
      
    integer, intent(in) :: nz ! Points in the vertical

    pver  = nz - 1
    pverp = nz

  end subroutine init_ppgrid
  
!================================================================================================
  
end module ppgrid
