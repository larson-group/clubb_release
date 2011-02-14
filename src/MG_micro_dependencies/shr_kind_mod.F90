! $Id$
module shr_kind_mod
!
! Dummy module for importing variables into morrison-gettelman microphysics
!---------------------------------------------------------------------------------------------------

  use kinds, only: &
    shr_kind_r8 => real_kind ! Determine if data values are single or double precision
    
  implicit none

  private
    
  public :: shr_kind_r8

end module shr_kind_mod
