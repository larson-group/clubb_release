! $Id$
module cldwat2m_macro
!
! Dummy module for importing variables into morrison-gettelman microphysics
!---------------------------------------------------------------------------------------------------

  implicit none

  private

  public :: rhmini, rhmaxi
  
  ! These variables are not used anywhere in MG, they are just imported. Because of this we
  ! are setting them to dummy values
  integer :: rhmini = 0, rhmaxi = 0

end module cldwat2m_macro
