! $Id: sfc.F,v 1.3 2006-11-02 01:45:17 dschanen Exp $
#define SCLR_RT 2
#define SCLR_THETA 1

        module surface

        implicit none

        private ! Default to private

        public ::
     .  sfc_momentum_fluxes,
     .  sfc_thermo_fluxes,
     .  sfc_var

        contains

!------------------------------------------------------------------------
        subroutine sfc_momentum_fluxes( u, v, upwp_sfc, vpwp_sfc )

!       Description:
!       This subroutine computes surface momentum fluxes using aerodynamic
!       formulas.

!       References:
!       None
!------------------------------------------------------------------------

        implicit none

        ! External
        intrinsic :: sqrt

        ! Constant parameter
        real, parameter :: 
     .  ustar = 0.3

        ! Input variables
        real, intent(in) :: 
     .  u, ! u wind first level above ground    [m/s]
     .  v  ! v wind first level above ground    [m/s]

        ! Output Variables
        real, intent(out) :: 
     .  upwp_sfc,! sfc u momentum flux (m^2/s^2)
     .  vpwp_sfc ! sfc v momentum flux (m^2/s^2)

        ! Local Variables
        real :: M ! total wind speed above ground

        ! Computes fluxes

        M = sqrt( u*u + v*v )
        upwp_sfc = - ustar*ustar * u / M
        vpwp_sfc = - ustar*ustar * v / M

        return
        end subroutine sfc_momentum_fluxes

!------------------------------------------------------------------------
        subroutine sfc_thermo_fluxes( u, v, Tsfc, psfc, thlair, rtair,
     .                                wpthlp_sfc, wprtp_sfc,
     .                                sclrair, wpsclrp_sfc )
!       Description:
!       This subroutine computes surface fluxes of heat and moisture 
!       using aerodynamic formulas.

!       References:
!       None
!------------------------------------------------------------------------

        use constants

        implicit none
        
        ! External
        intrinsic :: present, sqrt

        real, external :: rsat

        ! Parameter
        real, parameter :: C = 1.3e-3

        ! Input Variables
        real, intent(in) :: 
     .  u,      ! u wind                        [m/s]
     .  v,      ! u wind                        [m/s]
     .  Tsfc,   ! Surface temperature           [K]
     .  psfc,   ! Surface pressure              [Pa]
     .  thlair, ! theta_l at first model layer  [K]
     .  rtair   ! rt at first model layer       [kg/kg]


        ! Output Variables
        real, intent(out) :: 
     .  wpthlp_sfc,! surface thetal flux        [K m/s]
     .  wprtp_sfc  ! surface moisture flux      [kg/kg m/s]

        ! Input Variables (optional)
        real, optional, intent(in), dimension(sclr_dim) :: 
     .  sclrair ! Passive scalars at first model layer  [units vary]

        ! Output Variables (optional) 
        real, optional, intent(out), dimension(sclr_dim) :: 
     .  wpsclrp_sfc ! surface passive scalar flux       [units m/s]

        ! Local Variables
        real :: M  ! Total wind speed above ground

        ! Compute fluxes
        M = sqrt( u*u + v*v )
        wpthlp_sfc = -C * M * ( thlair - Tsfc * (psfc/p0)**kappa )
        wprtp_sfc  = -C * M * ( rtair - rsat( psfc, Tsfc ) )

        ! Compute passive scalar as == to wprtp & wpthlp for now
        if ( present( wpsclrp_sfc ) .and. present( sclrair ) ) then
          wpsclrp_sfc(SCLR_THETA) = -C * M *
     .               ( sclrair(SCLR_THETA) - Tsfc*(psfc/p0)**kappa )
          wpsclrp_sfc(SCLR_RT) = -C * M *
     .               ( sclrair(SCLR_RT) - rsat(psfc, Tsfc) )
        endif

        return
        end subroutine sfc_thermo_fluxes

!------------------------------------------------------------------------
        subroutine sfc_var( upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc,
     .                      wp2_sfc, thlp2_sfc, rtp2_sfc, rtpthlp_sfc, 
     .                      wpsclrp_sfc, sclrp2_sfc,
     .                      sclrprtp_sfc, sclrpthlp_sfc )
!       Description:
!       This subroutine computes estimate of the surface thermodynamic 
!       second order moments

!       References:
!       None
!------------------------------------------------------------------------

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
!        uf = sqrt( ustar2 + 0.3 * wstar**2 )  !test
        uf = max( ufmin, uf )

        ! Compute estimate for surface second order moments

        wp2_sfc     =  a * uf**2
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
            select case( i ) 
            case( SCLR_RT )
              sclrp2_sfc(i)   = 0.4 * a * ( wpsclrp_sfc(i) / uf )**2
              sclrprtp_sfc(i) = 0.4 * a * ( wpsclrp_sfc(i) / uf )**2
            case( SCLR_THETA )
              sclrp2_sfc(i)    =  0.1 * a * ( wpsclrp_sfc(i) / uf )**2
              sclrpthlp_sfc(i) =  0.1 * a * ( wpsclrp_sfc(i) / uf )**2
            end select
          end do ! 1,...sclr_dim
        end if

        return
        end subroutine sfc_var

        end module surface
