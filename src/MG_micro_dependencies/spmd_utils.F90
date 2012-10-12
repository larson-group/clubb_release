! $Id$
module spmd_utils
!
! Dummy module for importing variables into morrison_gettelman microphysics
!---------------------------------------------------------------------------------------------------

  implicit none

  private

  public :: masterproc

  ! This variable is not used anywhere in MG, it is just imported. Because of this we
  ! are setting it to a dummy value so MG will compile.
  logical :: masterproc = .false.

end module spmd_utils
