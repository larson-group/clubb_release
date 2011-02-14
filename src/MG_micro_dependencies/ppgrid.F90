! $Id$
module ppgrid
!
! Dummy module for importing variables into morrison-gettelman microphysics
!---------------------------------------------------------------------------------------------------

  implicit none

  private

  public :: pcols, pver, pverp, init_ppgrid

  integer :: &
    pcols = 1,   &  ! number of points in the “i” direction. CLUBB is a single column model.
    pver,        &  ! number of points in the vertical
    pverp           ! pver + 1
    
  contains
  
!================================================================================================
  subroutine init_ppgrid
    !
    !  Description: Initialize ppgrid variables for being imported into MG microphysics.
    !               We need to do this in a separate subroutine because Fortran doesn't
    !               permit using 'gr' in the variable initialization statements.
    !
    !---------------------------------------------------------------------------------
    use grid_class, only: gr
      
    pver  = gr%nnzp - 1
    pverp = gr%nnzp

  end subroutine init_ppgrid
  
!================================================================================================
  
end module ppgrid
