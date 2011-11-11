! $Id$
module cldwat2m_macro
!
! Dummy module for importing variables into morrison-gettelman microphysics
!---------------------------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  private

  public :: rhmini, rhmaxi
  
  ! These variables are not used anywhere in MG, they are just imported. Because of this we
  ! are setting them to dummy values
  real(r8), parameter :: rhmini = 0.80_r8 ! Minimum rh for ice cloud fraction
  real(r8), parameter :: rhmaxi = 1.1_r8  ! rhi at which ice cloud fraction = 1

end module cldwat2m_macro
