!----------------------------------------------------------------------
! $Id$
module gabls3_night

  ! Description:
  !   Contains subroutines for the GABLS3 LES.
  !
  ! References:
  !   http://www4.ncsu.edu/~sbasu5/GABLS3/
  !----------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  public :: gabls3_night_sfclyr

  private :: landflx

  private

  contains

  !-----------------------------------------------------------------------
  subroutine gabls3_night_sfclyr( ngrdcol, time, um_sfc, vm_sfc,  &
                                  thlm_sfc, rtm_sfc, lowest_level, & 
                                  upwp_sfc, vpwp_sfc, &
                                  wpthlp_sfc, wprtp_sfc, ustar )
    ! Description:
    !   This subroutine computes surface fluxes of horizontal momentum,
    !   heat and moisture according to GCSS ATEX specifications
    !
    ! References:
    !   http://www4.ncsu.edu/~sbasu5/GABLS3/
    !----------------------------------------------------------------------

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    use sfc_flux, only: compute_momentum_flux ! Procedure(s)

    use time_dependent_input, only: l_t_dependent,    & ! Variable(s)
                                    l_input_xpwp_sfc, &
                                    time_sfc_given,   &
                                    thlm_sfc_given,   &
                                    rtm_sfc_given,    &
                                    upwp_sfc_given,   &
                                    vpwp_sfc_given,   &
                                    time_select         ! Procedure(s)

    use interpolation, only: linear_interp_factor ! Procedure(s)

    implicit none

    ! Constants
    real( kind = core_rknd ), parameter ::  & 
      z0 = 0.15_core_rknd ! Roughness length  [m]

    ! Input variables
    integer, intent(in) :: &
      ngrdcol

    real(kind=time_precision), intent(in) :: &
      time ! Model time [s]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  & 
      um_sfc,       & ! um at zt(2)                     [m/s]
      vm_sfc,       & ! vm at zt(2)                     [m/s]
      thlm_sfc,     & ! Theta_l at zt(2)                [K]
      rtm_sfc,      & ! rt at zt(2)                     [kg/kg]
      lowest_level    ! gr%zt(1,2), height of the lowest
                      ! above-ground gridpoint          [m]

    ! Output variables
    real( kind = core_rknd ), dimension(ngrdcol), intent(out) ::  & 
      upwp_sfc,    & ! turbulent upward flux of u-momentum  [m^2/s^2]
      vpwp_sfc,    & ! turbulent upward flux of v-momentum  [m^2/s^2]
      wpthlp_sfc,  & ! w'theta_l' surface flux   [(m K)/s]
      wprtp_sfc,   & ! w'rt' surface flux        [(m kg)/(kg s)]
      ustar          ! surface friction velocity            [m/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      qs,   & ! Vapor at height z0
      ts      ! Potential temp. at height z0

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      ubar   ! Average surface wind speed [m/s]

    real( kind = core_rknd ) :: &
      time_frac, & ! time interpolation factor
      upwp_sfc_interp, &
      vpwp_sfc_interp

    integer :: before_time, after_time, i

    !$acc enter data create( ubar )

    if( l_t_dependent ) then

      ! Use time_select to determine the time indexes before and after 'time', as well as
      ! the time fraction necessary for linear_interp_factor
      call time_select( time, size( time_sfc_given ), time_sfc_given, &
                        before_time, after_time, time_frac )


      ts = linear_interp_factor( time_frac, thlm_sfc_given(after_time), &
                                 thlm_sfc_given(before_time) )

      qs = linear_interp_factor( time_frac, rtm_sfc_given(after_time), rtm_sfc_given(before_time) )
     
      ! Compute heat and moisture fluxes
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        call landflx( thlm_sfc(i), ts, rtm_sfc(i), qs,   & ! Intent(in)
                      um_sfc(i), vm_sfc(i), lowest_level(i), z0,     & ! Intent(in)
                      wpthlp_sfc(i), wprtp_sfc(i), ubar(i), ustar(i) )    ! Intent(out)
      end do

      if ( l_input_xpwp_sfc ) then

        ! Feed in momentum fluxes
        upwp_sfc_interp = linear_interp_factor( time_frac, upwp_sfc_given(after_time), &
                                            upwp_sfc_given(before_time) )  

        vpwp_sfc_interp = linear_interp_factor( time_frac, vpwp_sfc_given(after_time), &
                                            vpwp_sfc_given(before_time) )

        !$acc parallel loop gang vector default(present)
        do i = 1, ngrdcol
          upwp_sfc(i) = upwp_sfc_interp
          vpwp_sfc(i) = vpwp_sfc_interp
        end do

      else
        ! Compute momentum fluxes
        call compute_momentum_flux( ngrdcol, um_sfc, vm_sfc, ubar, ustar, & ! Intent(in)
                                    upwp_sfc, vpwp_sfc )           ! Intent(out)
      end if ! l_input_xpwp_sfc

    end if ! l_t_dependent

    !$acc exit data delete( ubar )

    return

  end subroutine gabls3_night_sfclyr

  !-----------------------------------------------------------------------------------------------
  real( kind = core_rknd ) function psi_h( x, xlmo )

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none
    real( kind = core_rknd ), intent(in) :: x
    real( kind = core_rknd ), intent(in) :: xlmo ! Monin-Obukhov length [m]

    psi_h = ( -5._core_rknd * x )/xlmo

  end function psi_h

  !-----------------------------------------------------------------------------------------------
  real( kind = core_rknd) function gm1( x )

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

    real( kind = core_rknd ), intent(in) :: x

    gm1 = (1._core_rknd-15._core_rknd*x)**0.25_core_rknd
  end function gm1

  !-----------------------------------------------------------------------------------------------
  real( kind = core_rknd ) function gh1( x )

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

    real( kind = core_rknd ), intent(in) :: x

    gh1=sqrt(1._core_rknd-9._core_rknd*x)/0.74_core_rknd

  end function gh1

  !-----------------------------------------------------------------------------------------------
  real( kind = core_rknd ) function fm1( x )

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

    real( kind = core_rknd ), intent(in) :: x

    real( kind = core_rknd ) :: pii

    pii=acos(-1._core_rknd)/2._core_rknd

    fm1 = 2._core_rknd*real(alog((1.+real(x))/2.), kind = core_rknd)+ &
      real(alog((1.+real(x*x))/2.), kind = core_rknd)-2._core_rknd*atan(x)+pii

  end function fm1

  !-----------------------------------------------------------------------------------------------
  real( kind = core_rknd ) function fh1( x )

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

    real( kind = core_rknd ), intent(in) :: x


    fh1 = 2._core_rknd*real(alog((1.+0.74*real(x))/2.0), kind = core_rknd)

  end function fh1

  !------------------------------------------------------------------------------------------------
  subroutine landflx( th, ts, qh, qs, uh, vh, h, z0, &
                      shf, lhf, vel, ustar )
    !$acc routine seq
    
    !
    !  Description: landflx.F90 from SAM 6.7.5
    !
    !----------------------------------------------------------------------------------------------

    use constants_clubb, only: ep1

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

    ! Input:

    real( kind = core_rknd ), intent(in) :: th ! pot. temperature at height h [K]
    real( kind = core_rknd ), intent(in) :: ts ! pot. temperature at z0       [K]
    real( kind = core_rknd ), intent(in) :: qh ! vapor at height h            [kg/kg]
    real( kind = core_rknd ), intent(in) :: qs ! vapor at z0                  [kg/kg]
    real( kind = core_rknd ), intent(in) :: uh ! zonal wind speed at height h [m/s]
    real( kind = core_rknd ), intent(in) :: vh ! merid wind speed at height h [m/s]
    real( kind = core_rknd ), intent(in) :: h  ! height h                     [m]
    real( kind = core_rknd ), intent(in) :: z0 ! friction height              [m]

    ! Output:

    real( kind = core_rknd ), intent(out) :: shf   ! sensible heat flux (K m/s)
    real( kind = core_rknd ), intent(out) :: lhf   ! latent heat flux (m/s)
    real( kind = core_rknd ), intent(out) :: vel
    real( kind = core_rknd ), intent(out) :: ustar

    real( kind = core_rknd ) r, zody
    real( kind = core_rknd ) a, b, c, d
    real( kind = core_rknd ) xm, xh, xsi, xsi1, xsi2, fm, fh

    real( kind = core_rknd ) xlmo ! Monin-Obukhov length [m]

    integer iter


    zody=real(alog(real(h)/real(z0)), kind = core_rknd)

    vel = sqrt(max(0.5_core_rknd,uh**2+vh**2))
    r=9.81_core_rknd/ts*(th*(1._core_rknd+ep1*qh)-ts*(1._core_rknd+ep1*qs))*h/vel**2
    iter=0

    if( r < 0._core_rknd ) then

      xsi=0._core_rknd
      iter=iter+1
      xm=gm1(xsi)
      xh=gh1(xsi)
      fm=zody-fm1(xm)
      fh=0.74_core_rknd*(zody-fh1(xh))
      xsi1=r/fh*fm**2
      xsi=xsi1

      xsi=-abs(xsi)
      iter=iter+1
      xm=gm1(xsi)
      xh=gh1(xsi)
      fm=zody-fm1(xm)
      fh=0.74_core_rknd*(zody-fh1(xh))
      xsi1=r/fh*fm**2
      xsi=xsi1

      xsi=-abs(xsi)
      iter=iter+1
      xm=gm1(xsi)
      xh=gh1(xsi)
      fm=zody-fm1(xm)
      fh=0.74_core_rknd*(zody-fh1(xh))
      xsi1=r/fh*fm**2
      xsi=xsi1

    else
      a=4.8_core_rknd*4.8_core_rknd*r-1.00_core_rknd*6.35_core_rknd
      b=(2._core_rknd*r*4.8_core_rknd-1.00_core_rknd)*zody
      c=r*zody**2
      d=sqrt(b*b-4._core_rknd*a*c)
      xsi1=(-b+d)/a/2._core_rknd
      xsi2=(-b-d)/a/2._core_rknd
      xsi=real(amax1(real(xsi1),real(xsi2)), kind = core_rknd)
      fm=zody+4.8_core_rknd*xsi
      fh=1.00_core_rknd*(zody+7.8_core_rknd*xsi)
!  	a=4._core_rknd7*4._core_rknd7*r-0._core_rknd74*6._core_rknd35
!	b=(2._core_rknd*r*4._core_rknd7-0._core_rknd74)*zody
!	c=r*zody**2
!	d=sqrt(b*b-4*a*c)
!	xsi1=(-b+d)/a/2.
!	xsi2=(-b-d)/a/2.
!	xsi=amax1(xsi1,xsi2)
!	fm=zody+4.7_core_rknd*xsi
!	fh=0.74_core_rknd*(zody+6.35_core_rknd*xsi)
    end if


    vel = sqrt(uh**2+vh**2)
! Modification for GABLS3_night
! Specification states how to compute these
! Joshua Fasching January 2009
!shf=0.4_core_rknd**2/fm/fh*vel*(ts-th)
!lhf=0.4_core_rknd**2/fm/fh*vel*(qs-qh)

    ustar = 0.4_core_rknd/fm*vel

    if( xsi >= 0._core_rknd ) then
      xsi = max(1.e-5_core_rknd,xsi)
    else
      xsi = min(-1.e-5_core_rknd,xsi)
    end if
    xlmo = h/xsi
!------ Modification of GABLS3_night
    shf = ( 0.4_core_rknd * ustar * (ts-th) ) / &
             ( real(alog(real(h)/0.25), kind = core_rknd) - psi_h(h, xlmo) &
                 + psi_h(0.25_core_rknd, xlmo) )
    lhf = ( 0.4_core_rknd * ustar * (qs-qh) ) / &
             ( real(alog(real(h)/0.25), kind = core_rknd) - psi_h(h, xlmo) &
                 + psi_h(0.25_core_rknd, xlmo) )
!-----------------------------
    return
  end subroutine landflx

end module gabls3_night
