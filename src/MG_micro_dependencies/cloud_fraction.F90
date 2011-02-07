! $Id$
module cloud_fraction
!
! Dummy module for importing variables into morrison-gettelman microphysics
!---------------------------------------------------------------------------------------------------

  implicit none

  private

  public :: cldfrc_getparams

  ! This variable is not used anywhere in MG, it is just imported. Because of this we
  ! are setting it to a dummy value so MG will compile.
  integer :: &
        cldfrc_getparams = 0

end module cloud_fraction
