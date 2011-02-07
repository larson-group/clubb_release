! $Id$
module cam_logfile
!
! Dummy module for importing variables into morrison-gettelman microphysics
!---------------------------------------------------------------------------------------------------

  implicit none

  private

  public :: iulog
  
  ! This variable is not used anywhere in MG, it is just imported. Because of this we
  ! are setting it to a dummy value
  integer :: iulog = 0

end module cam_logfile
