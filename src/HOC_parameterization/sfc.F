! $Id: sfc.F,v 1.8 2007-06-06 15:09:36 dschanen Exp $
#define SCLR_RT 2
#define SCLR_THETA 1

        module surface_var

        implicit none

        private ! Default to private

        public ::
     .  sfc_var

        contains

!------------------------------------------------------------------------
        subroutine sfc_var( upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc,
     .                      wp2_sfc, up2_sfc, vp2_sfc,
     .                      thlp2_sfc, rtp2_sfc, rtpthlp_sfc, 
     .                      wpsclrp_sfc, sclrp2_sfc,
     .                      sclrprtp_sfc, sclrpthlp_sfc )

!       Description:
!       This subroutine computes estimate of the surface thermodynamic 
!       second order moments

!       References:
!       None
!------------------------------------------------------------------------S

        use parameters, only: T0
        use constants

        implicit none

        ! External
        intrinsic :: sqrt, present, max

        ! Constant Parameters
        real, parameter :: 
     .  a = 1.8,
     .  z = 1.0,
     .  ufmin = 0.0001,
     .  sclr_var_coef = 0.25 ! This value is made up! - Vince Larson 12 Jul 2005

        ! Input Variables
        real, intent(in) :: 
     .  upwp_sfc,    ! Surface u momentum flux  [m^2/s^2]
     .  vpwp_sfc,    ! Surface v momentum flux  [m^2/s^2]
     .  wpthlp_sfc,  ! Surface thetal flux      [K m/s]
     .  wprtp_sfc    ! Surface moisture flux    [kg/kg m/s]

        ! Input (Optional) 
        real, optional, intent(in), dimension(sclr_dim) :: 
     .  wpsclrp_sfc ! Passive scalar flux       [units m/s]

        ! Output Variables
        real, intent(out) :: 
     .  wp2_sfc,    ! Vertical velocity variance        [m^2/s^2]
     .  up2_sfc,    ! u'^2                              [m^2/s^2]
     .  vp2_sfc,    ! u'^2                              [m^2/s^2]
     .  thlp2_sfc,  ! thetal variance                   [K^2]
     .  rtp2_sfc,   ! rt variance                       [kg/kg]
     .  rtpthlp_sfc ! thetal rt covariance              [kg K]

        ! Output Variables (Optional) 
        real, optional, intent(out), dimension(sclr_dim) :: 
     .  sclrp2_sfc,   ! Passive scalar variance                 [units^2]
     .  sclrprtp_sfc, ! Passive scalar r_t covariance           [units kg/kg]
     .  sclrpthlp_sfc ! Passive scalar theta_l covariance       [units K]

        ! Local Variables 
        real :: ustar2, wstar
        real :: uf

        integer :: i ! Loop index

        ! Compute ustar^2

        ustar2 = sqrt( upwp_sfc * upwp_sfc + vpwp_sfc * vpwp_sfc )

        ! Compute wstar following Andre et al., 1976

        if ( wpthlp_sfc > 0 ) then
          wstar = ( 1.0/T0 * grav * wpthlp_sfc * z ) ** (1./3.)
        else
          wstar = 0.
        end if

        ! Surface friction velocity following Andre et al. 1978

        uf = sqrt( ustar2 + 0.3 * wstar * wstar ) 
        uf = max( ufmin, uf )

        ! Compute estimate for surface second order moments

        wp2_sfc     =  a * uf**2
        up2_sfc     =  a * uf**2  ! From Andre, et al. 1978
        vp2_sfc     =  a * uf**2  ! "  "
        thlp2_sfc   = 0.1 * a * ( wpthlp_sfc / uf )**2
        rtp2_sfc    = 0.4 * a * ( wprtp_sfc / uf )**2
        rtpthlp_sfc = a * ( wpthlp_sfc / uf ) * ( wprtp_sfc / uf )

        ! Optional mixing scheme scalars
        if ( present( wpsclrp_sfc ) .and. present( sclrprtp_sfc ) 
     .       .and. present( sclrpthlp_sfc ) ) then
          do i=1, sclr_dim
            sclrprtp_sfc(i) 
     .      = a * (wprtp_sfc / uf) * (wpsclrp_sfc(i) / uf)
            sclrpthlp_sfc(i) 
     .      = a * (wpthlp_sfc / uf) * (wpsclrp_sfc(i) / uf)
            sclrp2_sfc(i)
     .      = sclr_var_coef * a * ( wpsclrp_sfc(i) / uf )**2

        ! Note: I have no idea why rtp2 and thlp2 have these coefficients
        ! -dschanen 7/6/05
            ! The following cases used the formulas for rt and theta in order
            ! to test the scalar code.  However, they need to be commented out
            ! because they do not apply to such things as CO2.
!            select case( i ) 
!            case( SCLR_RT )
!              sclrp2_sfc(i)   = 0.4 * a * ( wpsclrp_sfc(i) / uf )**2
!              sclrprtp_sfc(i) = 0.4 * a * ( wpsclrp_sfc(i) / uf )**2
!            case( SCLR_THETA )
!              sclrp2_sfc(i)    =  0.1 * a * ( wpsclrp_sfc(i) / uf )**2
!              sclrpthlp_sfc(i) =  0.1 * a * ( wpsclrp_sfc(i) / uf )**2
!            end select
          end do ! 1,...sclr_dim
        end if

        return
        end subroutine sfc_var

        end module surface_var
