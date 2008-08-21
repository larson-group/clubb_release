! $Id: sfc.F90,v 1.1 2008-07-22 16:04:28 faschinj Exp $
#define SCLR_RT 2
#define SCLR_THETA 1

        module surface_var

        implicit none

        private ! Default to private

        public ::  & 
        sfc_var

        contains

!------------------------------------------------------------------------
        subroutine sfc_var( upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc, & 
                            ustar_in, wpsclrp_sfc,  & 
                            wp2_sfc, up2_sfc, vp2_sfc, & 
                            thlp2_sfc, rtp2_sfc, rtpthlp_sfc, err_code, & 
                            sclrp2_sfc, & 
                            sclrprtp_sfc,  & 
                            sclrpthlp_sfc )

!       Description:
!       This subroutine computes estimate of the surface thermodynamic 
!       second order moments

!       References:
!       None
!------------------------------------------------------------------------S

        use parameters, only:  & 
            T0 ! Variable(s)
        use constants, only: & 
            grav,  & ! Variable(s)
            fstderr
        use parameters, only: & 
            sclr_dim  ! Variable(s)
        use numerical_check, only: & 
            sfc_var_check ! Procedure
        use error_code, only:  & 
            clubb_var_equals_NaN,  & ! Variable(s)
            clubb_at_debug_level        
        implicit none

        ! External
        intrinsic :: sqrt, max

        ! Constant Parameters
        real, parameter ::  & 
        a = 1.8, & 
        z = 1.0, & 
! Vince Larson increased ufmin to stabilize arm_97.  24 Jul 2007
!     .  ufmin = 0.0001,
        ufmin = 0.01, & 
! End Vince Larson's change.
! Vince Larson changed in order to make correlations between [-1,1].  31 Jan 2008.
!     .  sclr_var_coef = 0.25, ! This value is made up! - Vince Larson 12 Jul 2005
        sclr_var_coef = 0.4,  & ! This value is made up! - Vince Larson 12 Jul 2005
! End Vince Larson's change
! Vince Larson reduced surface spike in scalar variances associated w/ Andre et al. 1978 scheme
        reduce_coef   = 0.2

        ! Input Variables
        real, intent(in) ::  & 
        upwp_sfc,     & ! Surface u momentum flux   [m^2/s^2]
        vpwp_sfc,     & ! Surface v momentum flux   [m^2/s^2]
        wpthlp_sfc,   & ! Surface thetal flux       [K m/s]
        wprtp_sfc,    & ! Surface moisture flux     [kg/kg m/s]
        ustar_in     ! Surface friction velocity [m/s]

        ! Input (Optional) 
        real, intent(in), dimension(sclr_dim) ::  & 
        wpsclrp_sfc ! Passive scalar flux       [units m/s]

        ! Output Variables
        real, intent(out) ::  & 
        wp2_sfc,     & ! Vertical velocity variance        [m^2/s^2]
        up2_sfc,     & ! u'^2                              [m^2/s^2]
        vp2_sfc,     & ! u'^2                              [m^2/s^2]
        thlp2_sfc,   & ! thetal variance                   [K^2]
        rtp2_sfc,    & ! rt variance                       [kg/kg]
        rtpthlp_sfc ! thetal rt covariance              [kg K]

        integer, intent(out) :: & 
        err_code

        ! Output Variables (Optional) 
        real, intent(out), dimension(sclr_dim) ::  & 
        sclrp2_sfc,    & ! Passive scalar variance                 [units^2]
        sclrprtp_sfc,  & ! Passive scalar r_t covariance           [units kg/kg]
        sclrpthlp_sfc ! Passive scalar theta_l covariance       [units K]

        ! Local Variables 
        real :: ustar2, wstar
        real :: uf

        ! Variables for Andre et al., 1978 parameterization.
        real :: ustar
        real :: Lngth
        real :: zeta

        integer :: i ! Loop index

        ! Logical for Andre et al., 1978 parameterization.
        logical, parameter :: landre_1978 = .false.

        IF ( landre_1978 ) THEN

           ! For cases where ustar isn't specified.
           IF ( ustar_in == 0.0 ) THEN
              ustar = MAX( ( upwp_sfc**2 + vpwp_sfc**2 )**(1.0/4.0), & 
                           ufmin )
           ELSE
              ustar = MAX( ustar_in, ufmin )
           ENDIF

           ! Find Monin-Obukhov Length (Andre et al., 1978, p. 1866).
           Lngth = - ( ustar**3 ) /  & 
                     ( 0.35 * (1.0/T0) * grav * wpthlp_sfc )

           ! Find the value of dimensionless height zeta 
           ! (Andre et al., 1978, p. 1866).
           zeta = z / Lngth

           ! Andre et al, 1978, Eq. 29.
           ! Notes:  1) "reduce_coef" is a reduction coefficient intended to 
           !            make the values of rtp2, thlp2, and rtpthlp smaller 
           !            at the surface.
           !         2) With the reduction coefficient having a value of 0.2,
           !            the surface correlations of both w & rt and w & thl
           !            have a value of about 0.845.  These correlations are
           !            greater if zeta < 0.  The correlations have a value
           !            greater than 1 if zeta <= -0.212.
           !         3) The surface correlation of rt & thl is 1.
           ! Brian Griffin; February 2, 2008.
           IF ( zeta < 0.0 ) THEN
              thlp2_sfc   = reduce_coef  & 
                           * ( wpthlp_sfc**2 / ustar**2 ) & 
                           * 4.0 * ( 1.0 - 8.3*zeta )**(-2.0/3.0)
              rtp2_sfc    = reduce_coef  & 
                           * ( wprtp_sfc**2 / ustar**2 ) & 
                           * 4.0 * ( 1.0 - 8.3*zeta )**(-2.0/3.0)
              rtpthlp_sfc = reduce_coef  & 
                           * ( wprtp_sfc*wpthlp_sfc / ustar**2 ) & 
                           * 4.0 * ( 1.0 - 8.3*zeta )**(-2.0/3.0)
              wp2_sfc     = ( ustar**2 ) & 
                           * ( 1.75 + 2.0*(-zeta)**(2.0/3.0) )
           ELSE
              thlp2_sfc   = reduce_coef  & 
                           * 4.0 * ( wpthlp_sfc**2 / ustar**2 )
              rtp2_sfc    = reduce_coef  & 
                           * 4.0 * ( wprtp_sfc**2 / ustar**2 )
              rtpthlp_sfc = reduce_coef  & 
                           * 4.0 * ( wprtp_sfc*wpthlp_sfc / ustar**2 )
              wp2_sfc     = 1.75 * ustar**2
           ENDIF

           ! Calculate wstar following Andre et al., 1978, p. 1866.
           wstar = ( (1.0/T0) * grav * wpthlp_sfc * z )**(1.0/3.0)

           ! Andre et al, 1978, Eq. 29.
           IF ( wpthlp_sfc > 0.0 ) THEN
              !up2_sfc = 4.0 * ustar**2 + 0.3 * wstar**2
              !vp2_sfc = 1.75 * ustar**2 + 0.3 * wstar**2
              up2_sfc = 1.75 * ustar**2 + 0.3 * wstar**2
              vp2_sfc = 1.75 * ustar**2 + 0.3 * wstar**2
           ELSE
              !up2_sfc = 4.0 * ustar**2
              !vp2_sfc = 1.75 * ustar**2
              up2_sfc = 1.75 * ustar**2
              vp2_sfc = 1.75 * ustar**2
           ENDIF

           ! Passive scalars
           IF ( sclr_dim > 0 ) THEN
              DO i = 1, sclr_dim
                 ! Notes:  1) "reduce_coef" is a reduction coefficient intended 
                 !            to make the values of sclrprtp, sclrpthlp, and 
                 !            sclrp2 smaller at the surface.
                 !         2) With the reduction coefficient having a value 
                 !            of 0.2, the surface correlation of w & sclr has 
                 !            a value of about 0.845.  The correlation is
                 !            greater if zeta < 0.  The correlation has a value
                 !            greater than 1 if zeta <= -0.212.
                 !         3) The surface correlations of both rt & sclr and 
                 !            thl & sclr are 1.
                 ! Brian Griffin; February 2, 2008.
                 IF ( zeta < 0.0 ) THEN
                    sclrprtp_sfc(i)  & 
                    = reduce_coef  & 
                     * ( wpsclrp_sfc(i)*wprtp_sfc / ustar**2 ) & 
                     * 4.0 * ( 1.0 - 8.3*zeta )**(-2.0/3.0)
                    sclrpthlp_sfc(i)  & 
                    = reduce_coef  & 
                     * ( wpsclrp_sfc(i)*wpthlp_sfc / ustar**2 ) & 
                     * 4.0 * ( 1.0 - 8.3*zeta )**(-2.0/3.0)
                    sclrp2_sfc(i)  & 
                    = reduce_coef   & 
                     * ( wpsclrp_sfc(i)**2 / ustar**2 ) & 
                     * 4.0 * ( 1.0 - 8.3*zeta )**(-2.0/3.0)
                 ELSE
                    sclrprtp_sfc(i)  & 
                    = reduce_coef  & 
                     * 4.0 * ( wpsclrp_sfc(i)*wprtp_sfc / ustar**2 )
                    sclrpthlp_sfc(i)  & 
                    = reduce_coef  & 
                     * 4.0 * ( wpsclrp_sfc(i)*wpthlp_sfc / ustar**2 )
                    sclrp2_sfc(i)  & 
                    = reduce_coef & 
                     * 4.0 * ( wpsclrp_sfc(i)**2 / ustar**2 )
                 ENDIF
              ENDDO ! 1,...sclr_dim
           ENDIF

        ELSE ! Previous code.

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
! Vince Larson changed to make correlations between [-1,1]  31 Jan 2008
!           thlp2_sfc   = 0.1 * a * ( wpthlp_sfc / uf )**2
!           rtp2_sfc    = 0.4 * a * ( wprtp_sfc / uf )**2
!           rtpthlp_sfc = a * ( wpthlp_sfc / uf ) * ( wprtp_sfc / uf )
           ! Notes:  1) With "a" having a value of 1.8, the surface 
           !            correlations of both w & rt and w & thl have 
           !            a value of about 0.878.
           !         2) The surface correlation of rt & thl is 0.5.
           ! Brian Griffin; February 2, 2008.
           thlp2_sfc   = 0.4 * a * ( wpthlp_sfc / uf )**2
           rtp2_sfc    = 0.4 * a * ( wprtp_sfc / uf )**2
           rtpthlp_sfc = 0.2 * a *  & 
                           ( wpthlp_sfc / uf ) * ( wprtp_sfc / uf )
! End Vince Larson's change.

           ! Passive scalars
           if ( sclr_dim > 0 ) then
              do i=1, sclr_dim
! Vince Larson changed coeffs to make correlations between [-1,1]. 31 Jan 2008
!                 sclrprtp_sfc(i) 
!     .           = a * (wprtp_sfc / uf) * (wpsclrp_sfc(i) / uf)
!                 sclrpthlp_sfc(i) 
!     .           = a * (wpthlp_sfc / uf) * (wpsclrp_sfc(i) / uf)
!                 sclrp2_sfc(i)
!     .           = sclr_var_coef * a * ( wpsclrp_sfc(i) / uf )**2
                 ! Notes:  1) With "a" having a value of 1.8 and 
                 !            "sclr_var_coef" having a value of 0.4, the 
                 !            surface correlation of w & sclr has a value 
                 !            of about 0.878.
                 !         2) With "sclr_var_coef" having a value of 0.4, 
                 !            the surface correlations of both rt & sclr and
                 !            thl & sclr are 0.5.
                 ! Brian Griffin; February 2, 2008.
                 sclrprtp_sfc(i)  & 
                 = 0.2 * a * (wprtp_sfc / uf) * (wpsclrp_sfc(i) / uf)
                 sclrpthlp_sfc(i)  & 
                 = 0.2 * a * (wpthlp_sfc / uf) * (wpsclrp_sfc(i) / uf)
                 sclrp2_sfc(i) & 
                 = sclr_var_coef * a * ( wpsclrp_sfc(i) / uf )**2
! End Vince Larson's change.

           ! Note: I have no idea why rtp2 and thlp2 have these coefficients
           ! -dschanen 7/6/05
              ! The following cases used the formulas for rt and theta in order
              ! to test the scalar code.  However, they need to be commented out
              ! because they do not apply to such things as CO2.
!             select case( i ) 
!             case( SCLR_RT )
!                sclrp2_sfc(i)   = 0.4 * a * ( wpsclrp_sfc(i) / uf )**2
!                sclrprtp_sfc(i) = 0.4 * a * ( wpsclrp_sfc(i) / uf )**2
!             case( SCLR_THETA )
!                sclrp2_sfc(i)    =  0.1 * a * ( wpsclrp_sfc(i) / uf )**2
!                sclrpthlp_sfc(i) =  0.1 * a * ( wpsclrp_sfc(i) / uf )**2
!             end select
              end do ! 1,...sclr_dim
           end if

        END IF
        
        if ( clubb_at_debug_level( 2 ) ) then

          call sfc_var_check( wp2_sfc, up2_sfc, vp2_sfc,  & 
                              thlp2_sfc, rtp2_sfc, rtpthlp_sfc, & 
                              err_code, & 
                              sclrp2_sfc, sclrprtp_sfc, sclrpthlp_sfc )
        
!       Error reporting
!       Joshua Fasching February 2008
          if ( err_code == clubb_var_equals_NaN ) then

            write(fstderr,*) "Error in sfc_var"
            write(fstderr,*) "Intent(in)"
                    
            write(fstderr,*) "upwp_sfc = ", upwp_sfc 
            write(fstderr,*) "vpwp_sfc = ", vpwp_sfc 
            write(fstderr,*) "wpthlp_sfc = ", wpthlp_sfc
            write(fstderr,*) "wprtp_sfc = ", wprtp_sfc
            write(fstderr,*) "ustar_in = ", ustar_in 
                    
            if ( sclr_dim > 0 ) then 
              write(fstderr,*) "wpsclrp_sfc = ", wpsclrp_sfc
            end if
                    
            write(fstderr,*) "Intent(out)"
                    
            write(fstderr,*) "wp2_sfc = ", wp2_sfc  
            write(fstderr,*) "up2_sfc = ", up2_sfc
            write(fstderr,*) "vp2_sfc = ", vp2_sfc
            write(fstderr,*) "thlp2_sfc = ", thlp2_sfc
            write(fstderr,*) "rtp2_sfc = ", rtp2_sfc
            write(fstderr,*) "rtpthlp_sfc = ", rtpthlp_sfc
                
            if ( sclr_dim > 0 ) then 
              write(fstderr,*) "sclrp2_sfc = ", sclrp2_sfc
              write(fstderr,*) "sclrprtp_sfc = ", sclrprtp_sfc
              write(fstderr,*) "sclrpthlp_sfc = ", sclrpthlp_sfc
            end if

          end if ! err_code == clubb_var_equals_NaN

        end if ! clubb_at_debug_level ( 2 )
        
        return

        end subroutine sfc_var

        end module surface_var
