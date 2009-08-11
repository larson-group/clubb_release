!----------------------------------------------------------------------
! $Id$
module gabls3_night

  !       Description:
  !       Contains subroutines for the GABLS3 LES.
  !----------------------------------------------------------------------

  implicit none

  public :: gabls3_night_sfclyr

  private :: landflx, psi_h, gm1, gh1, fm1, fh1

  private

  contains

  !-----------------------------------------------------------------------
  subroutine gabls3_night_sfclyr( time, um_sfc, vm_sfc,  &
                            thlm_sfc, rtm_sfc, lowest_level, & 
                            upwp_sfc, vpwp_sfc, &
                            wpthlp_sfc, wprtp_sfc, ustar )
    !       Description:
    !       This subroutine computes surface fluxes of horizontal momentum,
    !       heat and moisture according to GCSS ATEX specifications
    !
    !       References:
    !
    !----------------------------------------------------------------------

    use constants, only: kappa, grav, Rd, Cp, p0, Lv ! Variable(s)

    use stats_precision, only: time_precision ! Variable(s)

    use surface_flux, only: compute_momentum_flux ! Procedure(s)

    use time_dependent_input, only: l_t_dependent,  & ! Variable(s)
                                    time_sfc_given, &
                                    thlm_sfc_given, &
                                    rtm_sfc_given,  &
                                    time_select       ! Procedure(s)

    use interpolation, only: factor_interp ! Procedure(s)

    implicit none

    ! Constants
    real, parameter ::  & 
      z0 = 0.15 ! Roughness length  [m]

    ! Input variables
    real(kind=time_precision), intent(in) :: time ! Model time [s]

    real, intent(in) ::  & 
      um_sfc,       & ! um at zt(2)                     [m/s]
      vm_sfc,       & ! vm at zt(2)                     [m/s]
      thlm_sfc,     & ! Theta_l at zt(2)                [K]
      rtm_sfc,      & ! rt at zt(2)                     [kg/kg]
      lowest_level    ! gr%zt(2), height of the lowest
                      ! above-ground gridpoint          [m]

    ! Output variables
    real, intent(out) ::  & 
      upwp_sfc,    & ! turbulent upward flux of u-momentum  [m^2/s^2]
      vpwp_sfc,    & ! turbulent upward flux of v-momentum  [m^2/s^2]
      ustar          ! surface friction velocity            [m/s]

    real, intent(out):: &
      wpthlp_sfc,  & ! w'theta_l' surface flux   [(m K)/s]
      wprtp_sfc      ! w'rt' surface flux        [(m kg)/(kg s)]

    ! Local Variables
    real :: &
      ubar, & ! Average surface wind speed [m/s]
      qs,   & ! Vapor at height z0
      ts      ! Potential temp. at height z0

    real :: time_frac ! time interpolation factor

    integer :: i1, i2

    if( l_t_dependent ) then

      call time_select( time, size(time_sfc_given), time_sfc_given, i1, i2 )

      time_frac = real( (time - time_sfc_given(i1)) /  &  ! at the first time, time_frac=0;
            (time_sfc_given(i2) - time_sfc_given(i1)) )   ! at the second time, time_frac=1.


      ts = factor_interp( time_frac, thlm_sfc_given(i2), thlm_sfc_given(i1) )

      qs = factor_interp( time_frac, rtm_sfc_given(i2), rtm_sfc_given(i1) )

      ! Compute heat and moisture fluxes
      call landflx( thlm_sfc, ts, rtm_sfc, qs, um_sfc, vm_sfc, lowest_level, z0, & ! Intent(in)
                    wpthlp_sfc, wprtp_sfc, ubar, ustar )                           ! Intent(out)

      ! Compute momentum fluxes
      call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, & ! Intent(in)
                                  upwp_sfc, vpwp_sfc )           ! Intent(out)

    end if ! l_t_dependent

    return
  end subroutine gabls3_night_sfclyr

  !-----------------------------------------------------------------------------------------------
  real function psi_h( x, xlmo )
    implicit none
    real, intent(in) :: x
    real, intent(in) :: xlmo

    psi_h = ( -5. * x )/xlmo

  end function psi_h

  !-----------------------------------------------------------------------------------------------
  real function gm1( x )

    implicit none

    real, intent(in) :: x

    gm1 = (1.-15.*x)**0.25
  end function gm1

  !-----------------------------------------------------------------------------------------------
  real function gh1( x )

    implicit none

    real, intent(in) :: x

    gh1=sqrt(1.-9.*x)/0.74

  end function gh1

  !-----------------------------------------------------------------------------------------------
  real function fm1( x )
    implicit none

    real, intent(in) :: x

    real :: pii

    pii=acos(-1.)/2.

    fm1 = 2.*alog((1.+x)/2.)+alog((1.+x*x)/2.)-2.*atan(x)+pii

  end function fm1

  !-----------------------------------------------------------------------------------------------
  real function fh1( x )

    implicit none

    real, intent(in) :: x


    fh1 = 2.*alog((1.+0.74*x)/2.)

  end function fh1

  !------------------------------------------------------------------------------------------------
  subroutine landflx( th, ts, qh, qs, uh, vh, h, z0, &
                      shf, lhf, vel, ustar )
    !
    !  Description: landflx.F90 from SAM 6.7.5
    !
    !----------------------------------------------------------------------------------------------

    use constants, only: eps

    implicit none

    ! Input:

    real, intent(in) :: th ! pot. temperature at height h [K]
    real, intent(in) :: ts ! pot. temperature at z0       [K]
    real, intent(in) :: qh ! vapor at height h            [kg/kg]
    real, intent(in) :: qs ! vapor at z0                  [kg/kg]
    real, intent(in) :: uh ! zonal wind speed at height h [m/s]
    real, intent(in) :: vh ! merid wind speed at height h [m/s]
    real, intent(in) :: h  ! height h                     [m]
    real, intent(in) :: z0 ! friction height              [m]

    ! Output:

    real, intent(out) :: shf   ! sensible heat flux (K m/s)
    real, intent(out) :: lhf   ! latent heat flux (m/s)
    real, intent(out) :: vel
    real, intent(out) :: ustar

    real r, zody
    real a, b, c, d
    real xm, xh, xsi, xsi1, xsi2, fm, fh

    real xlmo

    integer iter


    zody=alog(h/z0)

    vel = sqrt(max(0.5,uh**2+vh**2))
    r=9.81/ts*(th*(1+eps*qh)-ts*(1.+eps*qs))*h/vel**2
    iter=0

    if( r < 0. ) then

      xsi=0.
      iter=iter+1
      xm=gm1(xsi)
      xh=gh1(xsi)
      fm=zody-fm1(xm)
      fh=0.74*(zody-fh1(xh))
      xsi1=r/fh*fm**2
      xsi=xsi1

      xsi=-abs(xsi)
      iter=iter+1
      xm=gm1(xsi)
      xh=gh1(xsi)
      fm=zody-fm1(xm)
      fh=0.74*(zody-fh1(xh))
      xsi1=r/fh*fm**2
      xsi=xsi1

      xsi=-abs(xsi)
      iter=iter+1
      xm=gm1(xsi)
      xh=gh1(xsi)
      fm=zody-fm1(xm)
      fh=0.74*(zody-fh1(xh))
      xsi1=r/fh*fm**2
      xsi=xsi1

    else
      a=4.8*4.8*r-1.00*6.35
      b=(2.*r*4.8-1.00)*zody
      c=r*zody**2
      d=sqrt(b*b-4*a*c)
      xsi1=(-b+d)/a/2.
      xsi2=(-b-d)/a/2.
      xsi=amax1(xsi1,xsi2)
      fm=zody+4.8*xsi
      fh=1.00*(zody+7.8*xsi)
!  	a=4.7*4.7*r-0.74*6.35
!	b=(2.*r*4.7-0.74)*zody
!	c=r*zody**2
!	d=sqrt(b*b-4*a*c)
!	xsi1=(-b+d)/a/2.
!	xsi2=(-b-d)/a/2.
!	xsi=amax1(xsi1,xsi2)
!	fm=zody+4.7*xsi
!	fh=0.74*(zody+6.35*xsi)
    end if


    vel = sqrt(uh**2+vh**2)
! Modification for GABLS3_night
! Specification states how to compute these
! Joshua Fasching January 2009
!shf=0.4**2/fm/fh*vel*(ts-th)
!lhf=0.4**2/fm/fh*vel*(qs-qh)

    ustar = 0.4/fm*vel

    if( xsi >= 0. ) then
      xsi = max(1.e-5,xsi)
    else
      xsi = min(-1.e-5,xsi)
    end if
    xlmo = h/xsi
!------ Modification of GABLS3_night
    shf = ( 0.4 * ustar * (ts-th) ) / &
             ( alog(h/0.25) - psi_h(h, xlmo) + psi_h(0.25, xlmo) )
    lhf = ( 0.4 * ustar * (qs-qh) ) / &
             ( alog(h/0.25) - psi_h(h, xlmo) + psi_h(0.25, xlmo) )
!-----------------------------
    return
  end subroutine landflx

end module gabls3_night
