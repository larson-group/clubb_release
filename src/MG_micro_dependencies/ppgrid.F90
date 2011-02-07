! $Id$
module ppgrid
!
! Dummy module for importing variables into morrison-gettelman microphysics
!---------------------------------------------------------------------------------------------------

  implicit none

  private

  public :: pcols, pver, pverp
  
  use grid_class, only: gr

  integer :: &
    pcols  = 1  ! number of points in the “i” direction. CLUBB is a single column model.
    pver = gr%nnzp - 1  ! number of points in the vertical
    pverp = gr%nnzp   ! pver + 1
