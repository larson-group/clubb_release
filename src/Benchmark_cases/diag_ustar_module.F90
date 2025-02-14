!-----------------------------------------------------------------------
!$Id$

module diag_ustar_module

  use clubb_precision, only: core_rknd ! Variable(s)

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
  real( kind = core_rknd ) function diag_ustar( z, bflx, wnd, z0 ) 

    !$acc routine seq

    use constants_clubb, only: vonk, pi ! Variable(s)

    use clubb_precision, only: core_rknd ! Variable(s)

    !, eps
    implicit none

    real( kind = core_rknd ), parameter      :: am   =  4.8_core_rknd   !   "          "         "
    real( kind = core_rknd ), parameter      :: bm   = 19.3_core_rknd   !   "          "         "

    real( kind = core_rknd ), intent (in)    :: z             ! height where u locates
    real( kind = core_rknd ), intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
    real( kind = core_rknd ), intent (in)    :: wnd           ! wind speed at z
    real( kind = core_rknd ), intent (in)    :: z0            ! momentum roughness height

    integer :: iterate
    real( kind = core_rknd )    :: lnz, klnz, c1, x, psi1, zeta, lmo, ustar

    lnz   = log( z / z0 )
    klnz  = vonk/lnz
    c1    = pi / 2.0_core_rknd - 3.0_core_rknd*log( 2.0_core_rknd )

    ustar =  wnd*klnz
    !      if (bflx /= 0.0_core_rknd) then
    if (abs(bflx) > 1.e-6_core_rknd) then
    !      if (abs(bflx) > 1._core_rknde-4) then
      do iterate=1,4
    !          lmo   = -bflx * vonk/(ustar**3 + eps)
        lmo   = -ustar**3 / ( vonk * bflx )
        zeta  = z/lmo
        if (zeta > 0._core_rknd) then
          if ( zeta > 1.e10_core_rknd ) then ! -dschanen UWM for large zeta
            ustar = 1e-10_core_rknd
            exit
          else
            ustar =  vonk*wnd  /(lnz + am*zeta)
          end if
        else
          x     = sqrt( sqrt( 1.0_core_rknd - bm*zeta ) )
          psi1  = 2._core_rknd*log( 1.0_core_rknd+x ) + &
            log( 1.0_core_rknd+x*x ) - 2._core_rknd*atan( x ) + c1
          ustar = wnd*vonk/(lnz - psi1)
        end if
      end do ! 1..4
    end if

    diag_ustar = ustar

    return
    
  end function diag_ustar

end module diag_ustar_module
