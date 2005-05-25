!------------------------------------------------------------------------
! This subroutine computes surface momentum fluxes using aerodynamic
! formulas.
!
        subroutine sfc_momentum_fluxes( u, v, upwp_sfc, vpwp_sfc )
        implicit none

! Parameter

!       This variable occurs in several parts of hoc, varying from .25 to 3.
!       is it supposed to be constant?  -dschanen 5/25/05
        real, parameter :: ustar = 0.3

! Input

        real, intent(IN) :: u         ! u wind (m/s) first level above ground
        real, intent(IN) :: v         ! v wind (m/s) first level above ground

! Output

        real, intent(OUT) :: upwp_sfc ! sfc u momentum flux (m^2/s^2)
        real, intent(OUT) :: vpwp_sfc ! sfc v momentum flux (m^2/s^2)

! Internal

        real M                        ! total wind speed above ground

! Computes fluxes

        M = sqrt( u*u + v*v )
        upwp_sfc = - ustar*ustar * u / M
        vpwp_sfc = - ustar*ustar * v / M

        return
        end subroutine sfc_momentum_fluxes

!------------------------------------------------------------------------
! This subroutine computes surface fluxes of heat and moisture using
! aerodynamic formulas.
!
        subroutine sfc_thermo_fluxes( u, v, Tsfc, psfc, thlair, rtair,
     +                                wpthlp_sfc, wprtp_sfc )

        use constants

        implicit none
        
! External

        real rsat
        external rsat

! Input

        real, intent(IN) :: u           ! u wind (m/s)
        real, intent(IN) :: v           ! u wind (m/s)
        real, intent(IN) :: Tsfc        ! surface temperature (K)
        real, intent(IN) :: psfc        ! surface pressure (Pa)
        real, intent(IN) :: thlair      ! thetal at first model layer (K)
        real, intent(IN) :: rtair       ! rt at first model layer (kg/kg)

! Output

        real, intent(OUT) :: wpthlp_sfc ! surface thetal flux ( K m/s )
        real, intent(OUT) :: wprtp_sfc  ! surface moisture flux ( kg/kg m/s )

! Parameter

        real, parameter :: C = 1.3e-3

! Internal

        real M                          ! total wind speed above ground

! Compute fluxes

        M = sqrt( u*u + v*v )
        wpthlp_sfc = -C * M * ( thlair - Tsfc * (psfc/p0)**kappa )
        wprtp_sfc  = -C * M * ( rtair - rsat(psfc,Tsfc) )

        return
        end subroutine sfc_thermo_fluxes

!------------------------------------------------------------------------
! This subroutine computes estimate of the surface thermodynamic second
! order moments
!
        subroutine sfc_var( upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc,
     +                      wp2_sfc, thlp2_sfc, rtp2_sfc, rtpthlp_sfc )

        use constants

        implicit none

! Input

        real, intent(IN) :: upwp_sfc     ! surface u momentum flux ( m^2/s^2 )
        real, intent(IN) :: vpwp_sfc     ! surface v momentum flux ( m^2/s^2 )
        real, intent(IN) :: wpthlp_sfc   ! surface thetal flux ( K m/s )
        real, intent(IN) :: wprtp_sfc    ! surface moisture flux ( kg/kg m/s )

! Output

        real, intent(OUT) :: wp2_sfc     ! vertical velocity variance
        real, intent(OUT) :: thlp2_sfc   ! thetal variance
        real, intent(OUT) :: rtp2_sfc    ! rt variance
        real, intent(OUT) :: rtpthlp_sfc ! thetal rt covariance
        
! Internal

        real, parameter :: a = 1.8
        real, parameter :: z = 1.
        real, parameter :: ufmin = 0.0001

        real ustar2, wstar
        real uf

! Compute ustar^2

        ustar2 = sqrt( upwp_sfc * upwp_sfc + vpwp_sfc * vpwp_sfc )

! Compute wstar following Andre et al., 1976

        if ( wpthlp_sfc > 0 ) then
          wstar = ( alpha * grav * wpthlp_sfc * z ) ** (1./3.)
        else
          wstar = 0.
        end if

! Surface friction velocity following Andre et al. 1978

        uf = sqrt( ustar2 + 0.3 * wstar * wstar ) 
        uf = max( ufmin, uf )

! Compute estimate for surface second order moments

        wp2_sfc =  a * uf**2
        thlp2_sfc = 0.1 * a * ( wpthlp_sfc / uf )**2
        rtp2_sfc = 0.4 * a * ( wprtp_sfc / uf )**2
        rtpthlp_sfc = a * ( wpthlp_sfc / uf ) * ( wprtp_sfc / uf )

        return
        end subroutine sfc_var
