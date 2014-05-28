! $Id$
module shr_kind_mod
!
! Dummy module for importing variables into morrison_gettelman microphysics
!---------------------------------------------------------------------------------------------------

! use kinds, only: &
!   shr_kind_r8 => real_kind ! Determine if data values are single or double precision
    
  implicit none

  private
    
  integer, parameter, public :: shr_kind_r8 = 8 ! Added to fix a PGI error -dschanen 6 June 2011

end module shr_kind_mod
