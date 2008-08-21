!-----------------------------------------------------------------------
!$Id$

module diag_ustar_mod

implicit none

public :: diag_ustar

private ! Default Scope

contains
! ----------------------------------------------------------------------
!
! DISCLAIMER : this code appears to be correct but has not been
!              very thouroughly tested. If you do notice any
!              anomalous behaviour then please contact Andy and/or
!              Bjorn
!
! Function diag_ustar:  returns value of ustar using the below
! similarity functions and a specified buoyancy flux (bflx) given in
! kinematic units
!
! phi_m (zeta > 0) =  (1 + am * zeta)
! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
!
! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
!
! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface
! Layer, in Workshop on Micormeteorology, pages 67-100.
!
! Code writen March, 1999 by Bjorn Stevens
!
real function diag_ustar( z, bflx, wnd, z0 ) 

use constants, only: grav, vonk, pi ! Variable(s)
!, eps
implicit none

real, parameter      :: am   =  4.8   !   "          "         "
real, parameter      :: bm   = 19.3   !   "          "         "

real, intent (in)    :: z             ! height where u locates
real, intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
real, intent (in)    :: wnd           ! wind speed at z
real, intent (in)    :: z0            ! momentum roughness height

integer :: iterate
real    :: lnz, klnz, c1, x, psi1, zeta, lmo, ustar

lnz   = log( z / z0 )
klnz  = vonk/lnz
c1    = pi / 2.0 - 3.0*log( 2.0 )

ustar =  wnd*klnz
!      if (bflx /= 0.0) then
if (abs(bflx) > 1.e-6) then
!      if (abs(bflx) > 1.e-4) then
  do iterate=1,4
!          lmo   = -bflx * vonk/(ustar**3 + eps)
    lmo   = -ustar**3 / ( vonk * bflx )
    zeta  = z/lmo
    if (zeta > 0.) then
      ustar =  vonk*wnd  /(lnz + am*zeta)
    else
      x     = sqrt( sqrt( 1.0 - bm*zeta ) )
      psi1  = 2.*log( 1.0+x ) + log( 1.0+x*x ) - 2.*atan( x ) + c1
      ustar = wnd*vonk/(lnz - psi1)
    end if
  end do
end if

diag_ustar = ustar

return
end function diag_ustar

end module diag_ustar_mod
