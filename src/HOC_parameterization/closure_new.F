! $Id: closure_new.F,v 1.13 2008-02-12 22:50:08 faschinj Exp $
        module pdf_closure

        implicit none

        public :: pdf_closure_new

        contains 
!------------------------------------------------------------------------
        subroutine pdf_closure_new
     .             ( p, exner, wm, wp2, wp3, Sc,
     .               rtm, rtp2, wprtp,
     .               thlm, thlp2, wpthlp,
     .               rtpthlp, wp4, wprtp2, wp2rtp,
     .               wpthlp2, wp2thlp, wprtpthlp,
     .               cf, rcm, wpthvp, wp2thvp, rtpthvp, 
     .               thlpthvp, wprcp, wp2rcp, rtprcp, thlprcp,
     .               rcp2, pdf_parms, crt1, crt2, cthl1, cthl2,isValid,
     .               sclrm, wpsclrp,                    !optional
     .               wpsclrprtp, wpsclrp2, sclrprtp,    ! "    "
     .               sclrp2, sclrpthlp, sclrpthvp,      ! "    "
     .               wpsclrpthlp, sclrprcp, wp2sclrp )  ! "    "

!       Description:
!       Subroutine that computes pdf parameters analytically.

!       Based of the original formulation, but with some tweaks
!       to remove some of the less realistic assumptions and
!       improve transport terms.

!       Corrected version that should remove inconsistency

!       References:
!       Eqn. 29, 30, 31, 32 & 33  on p. 3547 of 
!       ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!         Method and Model Description'' Golaz, et al. (2002)
!       JAS, Vol. 59, pp. 3540--3551.
!------------------------------------------------------------------------

        use constants, only: 
     .  pi,      ! The ratio of radii to their circumference
     .  Cp,      ! Dry air specific heat at constant p [J/kg/K]
     .  Lv,      ! Latent heat of vaporization         [J/kg]
     .  Rd,      ! Dry air gas constant                [J/kg/K]
     .  Rv,      ! Water vapor gas constant            [J/kg/K]
     .  ep,      ! Rd / Rv;     ep  = 0.622            [-]
     .  ep1,     ! (1.0-ep)/ep; ep1 = 0.61             [-]
     .  ep2,     ! 1.0/ep;      ep2 = 1.61             [-]
     .  wtol,    ! Tolerance for w                     [m/s]
     .  rttol,   ! Tolerance for r_t                   [kg/kg]
     .  thltol,  ! Tolerance for th_l                  [K]
     .  sstol,   ! Tolerance for pdf parameter         [kg/kg]
     .  sclrtol, ! Array of passive scalar tolerances  [units vary]
     .  sclr_dim, ! Number of passive scalar variables
     .  fstderr
        use parameters, only:
     .  beta,    ! Plume widths for th_l and r_t       [-]
     .  T0       ! Reference temperature               [K]

        use anl_erf, only: erf ! The error function
        
        use numerical_check, only: pdf_closure_new_check
        implicit none

        ! External
        real, external :: rsat

        intrinsic :: sqrt, exp, min, max, abs

        ! Input Variables
        real, intent(in) :: 
     .  p,      ! Pressure.                     [Pa] 
     .  exner,  ! Exner function.               [-]
     .  wm,     ! mean w                        [m/s] 
     .  wp2,    ! w'^2                          [m^2/s^2] 
     .  wp3,    ! w'^3                          [m^3/s^3]
     .  Sc,     ! Width of individual w plumes  [-]
     .  rtm,    ! Mean total water              [kg/kg]
     .  rtp2,   ! Total water mixing ratio      [kg/kg]
     .  wprtp,  ! w' r_t'                       [(kg m)(kg s)]
     .  thlm,   ! Mean th_l                     [K]
     .  thlp2,  ! th_l'^2                       [K^2]
     .  wpthlp, ! w' th_l'                      [(m K)/s]
     .  rtpthlp ! r_t' th_l'                    [(K kg)/kg]

        ! Input (Optional mixing scheme variables)
        real, optional,  dimension(sclr_dim), intent(in) :: 
     .  sclrm,     ! Mean passive scalar        [units vary]
     .  wpsclrp,   ! w' sclr'                   [units vary]
     .  sclrp2,    ! sclr'^2                    [units vary]
     .  sclrprtp,  ! sclr' r_t'                 [units vary]
     .  sclrpthlp  ! sclr' th_l'                [units vary]

        ! Output Variables

        real, intent(out) :: 
     .  wp4,            ! w'^4                  [m^4/s^4]
     .  wprtp2,         ! w' r_t'               [(m kg)/(s kg)]
     .  wp2rtp,         ! w'^2 r_t'             [(m^2 kg)/(s^2 kg)]
     .  wpthlp2,        ! w' th_l'^2            [(m K^2)/s]
     .  wp2thlp,        ! w'^2 th_l'            [(m^2 K)/s^2]
     .  cf,             ! Cloud fraction        [%]
     .  rcm,            ! Mean liquid water     [kg/kg]
     .  wpthvp,         ! Buoyancy flux         [(K m)/s] 
     .  wp2thvp,        ! w'^2 th_v'            [(m^2 K)/s^2]
     .  rtpthvp,        ! r_t' th_v'            [(kg K)/kg]
     .  thlpthvp,       ! th_l' th_v'           [K^2]
     .  wprcp,          ! w' r_c'               [(m kg)/(s kg)]
     .  wp2rcp,         ! w'^2 r_c'             [(m^2 kg)/(s^2 kg)]
     .  rtprcp,         ! r_t' r_c'             [(kg^2)/(kg^2)]
     .  thlprcp,        ! th_l' r_c'            [(K kg)/kg]
     .  rcp2,           ! r_c'^2                [(kg^2)/(kg^2)]
     .  wprtpthlp,      ! w' r_t' th_l'         [(m kg K)/(s kg)]
     .  crt1, crt2, 
     .  cthl1, cthl2

        real, intent(out), dimension(26) ::
     .  pdf_parms       ! pdf paramters         [units vary]

        logical, intent(out) ::
     .  isValid         ! Are the outputs usable numbers?

        ! Output (Optional passive scalar variables)

        real, optional, intent(out), dimension(sclr_dim) :: 
     .  sclrpthvp,
     .  sclrprcp,
     .  wpsclrp2,
     .  wpsclrprtp,
     .  wpsclrpthlp,
     .  wp2sclrp

        ! Local Variables 

        real :: 
     .  a,                      ! pdf parameter
     .  w1, w2, sw1, sw2,       ! pdf parameters
     .  thl1, thl2,             ! pdf parameters
     .  sthl1, sthl2,           ! pdf parameters
     .  rt1, rt2,               ! pdf parameters
     .  srt1, srt2,             ! pdf parameters 
     .  w1_n, w2_n, 
     .  thl1_n, thl2_n, 
     .  rt1_n, rt2_n

        ! Passive scalar variables
        real, dimension(sclr_dim) :: 
     .  sclr1, sclr2, 
     .  ssclr1, ssclr2,
     .  alpha_sclr, 
     .  sclr1_n, sclr2_n,
     .  rsclrthl, rsclrrt
 
        logical :: scalar_calc

        ! Quantities needed to predict higher order moments
        real :: 
     .  tl1, tl2, 
     .  beta1, beta2, 
     .  rsl1, rsl2, 
     .  es1, es2,
     .  ss1, ss2,
     .  s1, s2, 
     .  zeta1, zeta2,
     .  rc1, rc2, 
     .  R1, R2

        ! Sub-plume correlation coefficient between rt, thl
        ! varies between -1 < rrtthl < 1

        real :: rrtthl

        ! Thermodynamic quantity

        real :: BD

        ! Skewness of w

        real :: Skw 

        ! weighting factor for interpolation of srt1,2, sthl1,2

        real :: diffweight

        ! alpha coefficients

        real :: alpha_thl, alpha_rt

!       variables for a generalization of Chris Golaz' closure
!       varies width of plumes in theta_l, rt
        real :: width_factor_1, width_factor_2

        integer :: i   ! Index

!------------------------ Code begins ----------------------------------

        ! Check whether the optional scalars are present.
        ! Rather than test this on an equation to variable basis
        ! as in mixing, this uses an all or nothing test.

        if (      present( sclrm )     .and. present( wpsclrp )
     .      .and. present( sclrprtp )  .and. present( sclrpthvp )
     .      .and. present( sclrpthlp ) .and. present( sclrprcp )
     .      .and. present( wpsclrprtp ).and. present( wpsclrp2 ) 
     .      .and. present( sclrp2 )    .and. present( wpsclrpthlp )
     .      .and. present( wp2sclrp ) ) then
          scalar_calc = .true.
        else
          scalar_calc = .false.
        end if

        BD = Lv / (exner*Cp) - ep2 * T0

        ! If there is no velocity, then use single delta fnc. as pdf
        ! Otherwise compute parameters explicitly.

        if ( wp2 <= wtol**2 )  then

          a      = 0.5
          w1     = wm
          w2     = wm
          sw1    = 0.
          sw2    = 0.
          rt1    = rtm
          rt2    = rtm
          srt1   = 0.
          srt2   = 0.
          thl1   = thlm
          thl2   = thlm
          sthl1  = 0.
          sthl2  = 0.
          rrtthl = 0.

        else ! Explicitly 

          Skw = wp3/wp2**1.5

          if ( abs(Skw) <= 1e-5 ) then
            a = 0.5
          else
            a = 0.5
     .          * ( 
     .              1.0 - Skw/sqrt( 4.0*( 1.0 - Sc )**3 + Skw**2 )
     .            )
          end if

          ! Clip a, 1-a, to avoid dividing by zero
          if (a > 0.99) a = 0.99
          if (a < 0.01) a = 0.01

          w1_n = sqrt( ( (1-a)/a )*(1-Sc) )
          w2_n = -sqrt( ( a/(1-a) )*(1-Sc) )
          w1   = wm + sqrt( wp2 )*w1_n
          w2   = wm + sqrt( wp2 )*w2_n

          sw1  = Sc*wp2
          sw2  = Sc*wp2

          ! Vince Larson added a dimensionless factor so that the
          ! width of plumes in theta_l, rt can vary.  
          ! beta is a constant defined in module constants  
          !   Set 0<beta<3.
          ! beta=1.5 recovers Chris Golaz' simplified formula.
          ! 3 Nov 2003

          width_factor_1 = ( 2.0/3.0 )*beta 
     .                   + 2.0*a*( 1.0 - ( 2.0/3.0 )*beta )
          width_factor_2 = 2.0 - width_factor_1

          if ( thlp2 <= thltol**2 ) then
            thl1      = thlm
            thl2      = thlm
            sthl1     = 0.0
            sthl2     = 0.0
            alpha_thl = 0.5
          else
            thl1_n = - (wpthlp/(sqrt( wp2 )*sqrt( thlp2 )))/w2_n
            thl2_n = - (wpthlp/(sqrt( wp2 )*sqrt( thlp2 )))/w1_n

            thl1   = thlm - ( wpthlp/sqrt( wp2 ) )/w2_n
            thl2   = thlm - ( wpthlp/sqrt( wp2 ) )/w1_n

            alpha_thl
     .      = 0.5
     .        * ( 1.0 - wpthlp*wpthlp /((1.0-Sc)*wp2*thlp2) )
            alpha_thl = max( min( alpha_thl, 1.0 ), 0.0 )

        ! Vince Larson multiplied original expressions by width_factor_1,2
        !   to generalize scalar skewnesses.  05 Nov 03
            sthl1 = ( alpha_thl / a * thlp2 ) * width_factor_1
            sthl2 = ( alpha_thl / (1.-a) * thlp2 ) * width_factor_2

          end if ! thlp2 <= thltol**2

          if ( rtp2 <= rttol**2 ) then
            rt1      = rtm
            rt2      = rtm
            srt1     = 0.0
            srt2     = 0.0
            alpha_rt = 0.5
          else
            rt1_n = -( wprtp / ( sqrt( wp2 )*sqrt( rtp2 ) ) ) / w2_n
            rt2_n = -( wprtp / ( sqrt( wp2 )*sqrt( rtp2 ) ) ) / w1_n

            rt1 = rtm - ( wprtp / sqrt( wp2 ) ) / w2_n
            rt2 = rtm - ( wprtp / sqrt( wp2 ) ) / w1_n

            alpha_rt
     .      = 0.5
     .        * ( 1.0 - wprtp*wprtp /((1.0-Sc)*wp2*rtp2) )
            alpha_rt = max( min( alpha_rt, 1.0 ), 0.0 )

        ! Vince Larson multiplied original expressions by width_factor_1,2
        !   to generalize scalar skewnesses.  05 Nov 03
            srt1 = ( alpha_rt / a * rtp2 ) * width_factor_1
            srt2 = ( alpha_rt / (1.-a) * rtp2 ) * width_factor_2

          end if ! rtp2 <= rttol**2 

        ! Compute pdf parameters for passive scalars
          if ( scalar_calc ) then
            do i = 1, sclr_dim
              if ( sclrp2(i) <= sclrtol(i)**2 ) then
                sclr1(i)      = sclrm(i)
                sclr2(i)      = sclrm(i)
                ssclr1(i)     = 0.0
                ssclr2(i)     = 0.0
                alpha_sclr(i) = 0.5
              else
                sclr1_n(i) = - ( wpsclrp(i) / (sqrt( wp2 ) 
     .                           * sqrt( sclrp2(i) )) )/w2_n
                sclr2_n(i) = - ( wpsclrp(i) / (sqrt( wp2 ) 
     .                           * sqrt( sclrp2(i) )) )/w1_n
  
                sclr1(i) = sclrm(i) 
     .                   - ( wpsclrp(i) / sqrt( wp2 ) ) / w2_n
                sclr2(i) = sclrm(i) 
     .                   - ( wpsclrp(i) / sqrt( wp2 ) ) / w1_n

                alpha_sclr(i) 
     .          = 0.5
     .            * ( 1.0 - wpsclrp(i)*wpsclrp(i) 
     .                      / ((1.0-Sc)*wp2*sclrp2(i)) )
                alpha_sclr(i) = max( min( alpha_sclr(i), 1.0 ), 0.0 )

        ! Vince Larson multiplied original expressions by width_factor_1,2
        !  to generalize scalar skewnesses.  05 Nov 03
                ssclr1(i) = ( alpha_sclr(i) / a * sclrp2(i) ) 
     .                        * width_factor_1
                ssclr2(i) = ( alpha_sclr(i) / (1.-a) * sclrp2(i) ) 
     .                        * width_factor_2
              end if ! sclrp2(i) <= sclrtol(i)**2
            end do ! i=1, sclr_dim
          end if ! scalar_calc

        ! We include sub-plume correlation with coeff rrtthl.

          if ( srt1*sthl1 > 0 .and. srt2*sthl2 > 0) then
            rrtthl = ( rtpthlp 
     .             - a * ( rt1-rtm ) * ( thl1-thlm )
     .             - (1-a) * ( rt2-rtm ) * ( thl2-thlm ) )
     .                 / ( a*sqrt( srt1*sthl1 ) 
     .                     + (1-a)*sqrt( srt2*sthl2 ) )
            if ( rrtthl < -1.0 ) then
              rrtthl = -1.0
            end if
            if ( rrtthl > 1.0 ) then
              rrtthl = 1.0
            end if
          else
            rrtthl = 0.0
          end if ! srt1*sthl1 > 0 .and. srt2*sthl2 > 0

        ! Sub-plume correlation, rsclrthl, between passive scalar
        !      and theta_l.
          if ( scalar_calc ) then
            do i=1, sclr_dim
              if ( ssclr1(i)*sthl1 > 0 .and. ssclr2(i)*sthl2 > 0 ) then
                  rsclrthl(i) = ( sclrpthlp(i) 
     .             - a * ( sclr1(i)-sclrm(i) ) * ( thl1-thlm )
     .             - (1-a) * ( sclr2(i)-sclrm(i) ) * ( thl2-thlm ) )
     .                 / ( a*sqrt( ssclr1(i)*sthl1 ) 
     .                         + (1-a)*sqrt( ssclr2(i)*sthl2 ) )
                if ( rsclrthl(i) < -1.0 ) then
                  rsclrthl(i) = -1.0
                end if
                if ( rsclrthl(i) > 1.0 ) then
                  rsclrthl(i) = 1.0
                end if
              else
                rsclrthl(i) = 0.0
              end if

         ! Sub-plume correlation, rsclrrt, between passive scalar
         !   and total water.

              if ( ssclr1(i)*srt1 > 0 .and. ssclr2(i)*srt2 > 0 ) then
                rsclrrt(i) = ( sclrprtp(i) 
     .           - a * ( sclr1(i)-sclrm(i) ) * ( rt1-rtm )
     .           - (1-a) * ( sclr2(i)-sclrm(i) ) * ( rt2-rtm ) )
     .               / ( a*sqrt( ssclr1(i)*srt1 ) 
     .                   + (1-a)*sqrt( ssclr2(i)*srt2 ) )
                if ( rsclrrt(i) < -1.0 ) then
                  rsclrrt(i) = -1.0
                endif
                if ( rsclrrt(i) > 1.0 ) then
                  rsclrrt(i) = 1.0
                endif
              else
                rsclrrt(i) = 0.0
              end if
            end do ! i=1, sclr_dim
          end if ! scalar_calc 

        end if

        ! Compute higher order moments

        wp2rtp  = a * ( (w1-wm)**2+sw1 ) * ( rt1-rtm )
     .          + (1-a) * ( (w2-wm)**2+sw2 ) * ( rt2-rtm )

        wp2thlp = a * ( (w1-wm)**2+sw1 ) * ( thl1-thlm )
     .          + (1-a) * ( (w2-wm)**2+sw2 ) * ( thl2-thlm )

        wp4     = a * ( 3*sw1**2 + 6*((w1-wm)**2)*sw1 + (w1-wm)**4 )
     .          + (1-a) * ( 3*sw2**2 + 6*((w2-wm)**2)*sw2 + (w2-wm)**4 )

        wprtp2  = a * ( w1-wm )*( (rt1-rtm)**2 + srt1 ) 
     .          + (1-a) * ( w2-wm )*( (rt2-rtm)**2 + srt2)

        wpthlp2 = a * ( w1-wm )*( (thl1-thlm)**2 + sthl1 ) 
     .          + (1-a) * ( w2-wm )*( (thl2-thlm)**2+sthl2 )

        wprtpthlp = a * ( w1-wm )*( (rt1-rtm)*(thl1-thlm) 
     .                + rrtthl*sqrt( srt1*sthl1 ) )
     .            + ( 1-a ) * ( w2-wm )*( (rt2-rtm)*(thl2-thlm)
     .                + rrtthl*sqrt( srt2*sthl2 ) )

        ! Scalar Addition to higher order moments
        if ( scalar_calc ) then
          do i=1, sclr_dim

            wp2sclrp(i)  = a * ( (w1-wm)**2+sw1 )*( sclr1(i)-sclrm(i) )
     .          + (1-a) * ( (w2-wm)**2+sw2 ) * ( sclr2(i)-sclrm(i) )

            wpsclrp2(i) = a * ( w1-wm )*
     .                    ( (sclr1(i)-sclrm(i))**2 + ssclr1(i) ) 
     .                  + (1-a) * ( w2-wm )*
     .                    ( (sclr2(i)-sclrm(i))**2 + ssclr2(i) )

            wpsclrprtp(i) = a * ( w1-wm ) * 
     .      ( ( rt1-rtm )*( sclr1(i)-sclrm(i) ) 
     .         + rsclrrt(i)*sqrt( srt1*ssclr1(i) ) )
     .      + ( 1-a )*( w2-wm )*
     .      ( ( rt2-rtm )*( sclr2(i)-sclrm(i) )
     .         + rsclrrt(i)*sqrt( srt2*ssclr2(i) ) )

            wpsclrpthlp(i) = a * ( w1-wm ) * 
     .      ( ( sclr1(i)-sclrm(i) )*( thl1-thlm ) 
     .         + rsclrthl(i)*sqrt( ssclr1(i)*sthl1 ) )
     .      + ( 1-a ) * ( w2-wm ) *
     .      ( ( sclr2(i)-sclrm(i) )*( thl2-thlm )
     .         + rsclrthl(i)*sqrt( ssclr2(i)*sthl2 ) )


          end do ! i=1, sclr_dim
        end if ! scalar_calc

        ! Compute higher order moments that include theta_v.

        ! First compute some preliminary quantities.
        ! "1" denotes first Gaussian; "2" denotes 2nd Gaussian
        ! liq water temp (Sommeria & Deardorff 1977 (SD), eqn. 3)

        tl1  = thl1*exner
        tl2  = thl2*exner

        rsl1 = rsat(p,tl1)
        rsl2 = rsat(p,tl2)

        ! SD's beta (eqn. 8)
        beta1 = ep * ( Lv/(Rd*tl1) ) * ( Lv/(Cp*tl1) )
        beta2 = ep * ( Lv/(Rd*tl2) ) * ( Lv/(Cp*tl2) )

        ! s from Lewellen and Yoh 1993 (LY) eqn. 1
        s1 = ( rt1 - rsl1 ) / ( 1 + beta1 * rsl1 )
        s2 = ( rt2 - rsl2 ) / ( 1 + beta2 * rsl2 )

        ! Coefficient for s'
        crt1  = 1/( 1 + beta1*rsl1)
        crt2  = 1/( 1 + beta2*rsl2)

        cthl1 = ( (1 + beta1 * rt1) / ( 1 + beta1*rsl1)**2 )
     .          * ( Cp/Lv ) * beta1 * rsl1 * exner
        cthl2 = ( ( 1 + beta2 * rt2) / ( 1 + beta2*rsl2 )**2 )
     .          * ( Cp/Lv ) * beta2 * rsl2 * exner

        ! Standard deviation of s
        ! include subplume correlation of qt, thl
        ! Because of round-off error, 
        ! ss1 (and probably ss2) can become negative when rrtthl=1
        ! One could also write this as a squared term
        ! plus a postive correction; this might be a neater format

        ss1 = sqrt( max( 0.0, ( srt1*crt1**2 + sthl1*cthl1**2 
     .          - 2.0*rrtthl*crt1*sqrt( srt1 )*cthl1*sqrt( sthl1 ) ) 
     .                 ) 
     .            )
        ss2 = sqrt( max( 0.0, ( srt2*crt2**2 + sthl2*cthl2**2
     .          - 2.0*rrtthl*crt2*sqrt( srt2 )*cthl2*sqrt( sthl2 ) ) 
     .                 ) 
     .            )

!        ss1 = sqrt( (sqrt(srt1)*crt1 - sqrt(sthl1)*cthl1)**2 
!     .          + (1.-rrtthl)*2.*crt1*sqrt(srt1)*cthl1*sqrt(sthl1)  )
!        ss2 = sqrt( (sqrt(srt2)*crt2 - sqrt(sthl2)*cthl2)**2
!     .          + (1.-rrtthl)*2.*crt2*sqrt(srt2)*cthl2*sqrt(sthl2)  )


        ! We need to introduce a threshold value for the variance of s

        if ( ss1 > sstol ) then
           zeta1 = s1/ss1
           R1  = 0.5*( 1 + erf( zeta1/sqrt( 2.0 ) ) )
           rc1 = s1*R1+ss1*exp( -0.5*zeta1**2 )/( sqrt( 2*pi ) )
        else
           if ( s1 < 0.0 ) then
              R1  = 0.0
              rc1 = 0.0
           else
              R1  = 1.0
              rc1 = s1
           end if
        end if
        if ( ss2 > sstol ) then
           zeta2 = s2/ss2
           R2    = 0.5*( 1 + erf( zeta2/sqrt( 2.0 ) ) )
           rc2   = s2*R2+ss2*exp( -0.5*zeta2**2 )/( sqrt( 2*pi ) )
        else
           if ( s2 < 0.0 ) then
             R2  = 0.0
             rc2 = 0.0
           else
             R2  = 1.0
             rc2 = s2
           end if ! s2 < 0
        end if ! ss2 > sstol

        ! Compute moments that depend on theta_v

        wp2rcp  = a * ((w1-wm)**2 + sw1)*rc1
     .                 + (1-a) * ((w2-wm)**2 + sw2)*rc2
     .                 - wp2 * (a*rc1+(1-a)*rc2)
        wp2thvp = wp2thlp + ep1*T0*wp2rtp + BD*wp2rcp

        wprcp   = a * (w1-wm)*rc1 + (1-a) * (w2-wm)*rc2
        wpthvp  = wpthlp + ep1*T0*wprtp + BD*wprcp

        ! Account for subplume correlation in qt-thl
        thlprcp  = a * ( (thl1-thlm)*rc1 - (cthl1*sthl1)*R1 )
     .            + (1-a) * ( (thl2-thlm)*rc2 - (cthl2*sthl2)*R2 )
     .            + a*rrtthl*crt1*sqrt( srt1 )*sqrt( sthl1 )*R1
     .            + (1-a)*rrtthl*crt2*sqrt( srt2 )*sqrt( sthl2 )*R2
        thlpthvp = thlp2 + ep1*T0*rtpthlp + BD*thlprcp

        ! Account for subplume correlation in qt-thl
        rtprcp   = a * ( (rt1-rtm)*rc1 + (crt1*srt1)*R1 )
     .            + (1-a) * ( (rt2-rtm)*rc2 + (crt2*srt2)*R2 )
     .            - a*rrtthl*cthl1*sqrt( srt1 )*sqrt( sthl1 )*R1
     .            - (1-a)*rrtthl*cthl2*sqrt( srt2 )*sqrt( sthl2 )*R2

        rtpthvp  = rtpthlp + ep1*T0*rtp2 + BD*rtprcp

        ! Account for subplume correlation between scalar, theta_v.
        ! See Eqs. A13, A8 from Larson et al. (2002) ``Small-scale...''
        !  where the ``scalar'' in this paper is w.
        if ( scalar_calc ) then
          do i=1, sclr_dim
            sclrprcp(i) = a * ( ( sclr1(i)-sclrm(i) ) * rc1 )
     .                  + (1-a) * ( ( sclr2(i)-sclrm(i) ) * rc2 )
     .                  + a*rsclrrt(i) * crt1 
     .                    * sqrt( ssclr1(i) ) * sqrt( srt1 ) * R1
     .                  + (1-a) * rsclrrt(i) * crt2 
     .                    * sqrt( ssclr2(i) ) * sqrt( srt2 ) * R2
     .                  - a * rsclrthl(i) * cthl1 
     .                    * sqrt( ssclr1(i) )* sqrt( sthl1 ) * R1
     .                  - (1-a) * rsclrthl(i) * cthl2 
     .                    * sqrt( ssclr2(i) ) * sqrt( sthl2 ) * R2

            sclrpthvp(i) = sclrpthlp(i) + ep1*T0*sclrprtp(i) 
     .                   + BD*sclrprcp(i)
          end do ! i=1, sclr_dim
        end if ! scalar_calc

        ! Compute mean cloud fraction and cloud water

        cf  = a * R1 + (1-a) * R2
        rcm = a * rc1 + (1-a) * rc2

        ! Note: Brian added the following lines to ensure that there
        ! are never any negative liquid water values (or any negative
        ! cloud fraction values, for that matter).  According to
        ! Vince Larson, the analytic formula should not produce any
        ! negative results, but such computer-induced errors such as
        ! round-off error may produce such a value.  This has been
        ! corrected because Brian found a small negative value of
        ! rcm in the first timestep of the FIRE case.

        cf  = max( 0.0, cf )
        rcm = max( 0.0, rcm )

        ! Compute variance of liquid water mixing ratio.
        ! This is not needed for closure.  Statistical Analysis only.
        rcp2 =      a    *  ( s1*rc1 + R1*ss1**2 )
     .         + ( 1-a ) *  ( s2*rc2 + R2*ss2**2 )
     .         -   rcm**2

        rcp2 = max( 0.0, rcp2 )


        ! Save pdf parameters

        pdf_parms(1)  = w1
        pdf_parms(2)  = w2
        pdf_parms(3)  = sw1
        pdf_parms(4)  = sw2
        pdf_parms(5)  = rt1
        pdf_parms(6)  = rt2
        pdf_parms(7)  = srt1
        pdf_parms(8)  = srt2
        pdf_parms(9)  = thl1
        pdf_parms(10) = thl2
        pdf_parms(11) = sthl1
        pdf_parms(12) = sthl2
        pdf_parms(13) = a
        pdf_parms(14) = rc1
        pdf_parms(15) = rc2
        pdf_parms(16) = rsl1
        pdf_parms(17) = rsl2
        pdf_parms(18) = R1
        pdf_parms(19) = R2
        pdf_parms(20) = s1
        pdf_parms(21) = s2
        pdf_parms(22) = ss1
        pdf_parms(23) = ss2
        pdf_parms(24) = rrtthl
        pdf_parms(25) = alpha_thl
        pdf_parms(26) = alpha_rt
        call pdf_closure_new_check( wp4, wprtp2, wp2rtp, wpthlp2,
     .                           wp2thlp, cf, rcm, wpthvp, wp2thvp,
     .                           rtpthvp, thlpthvp, wprcp, wp2rcp,
     .                           rtprcp, thlprcp, rcp2, wprtpthlp,
     .                           crt1, crt2, cthl1, cthl2, pdf_parms,
     .                           isValid 
#ifdef SCALARS
     .                           ,sclrpthvp, sclrprcp, wpsclrp2,
     .                           wpsclrprtp, wpsclrpthlp,
     .                           wp2sclrp 
#endif /*SCALARS*/
     .                           )

!       Error Reporting
!       Joshua Fasching February 2008
                
        if ( .not. isValid ) then
                
           write(fstderr,*) "Error in pdf_closure_new"
           
           write(fstderr,*) "Intent(in)"
           
           write(fstderr,*) "p = ", p 
           write(fstderr,*) "exner = ", exner
           write(fstderr,*) "wm = ", wm 
           write(fstderr,*) "wp2 = ", wp2 
           write(fstderr,*) "wp3 = ", wp3
           write(fstderr,*) "Sc = ", Sc
           write(fstderr,*) "rtm = ", rtm
           write(fstderr,*) "rtp2 = ", rtp2
           write(fstderr,*) "wprtp = ", wprtp
           write(fstderr,*) "thlm = ", thlm
           write(fstderr,*) "thlp2 = ", thlp2
           write(fstderr,*) "wpthlp = ", wpthlp
           write(fstderr,*) "rtpthlp = ", rtpthlp

          if ( present( sclrm ) ) then
             write(fstderr,*) "sclrm = ", sclrm
          endif
          
          if ( present( wpsclrp ) ) then
             write(fstderr,*) "wpsclrp = ", wpsclrp
          endif
          
          if ( present( sclrp2 ) ) then
             write(fstderr,*) "sclrp2 = ", sclrp2
          endif
          
          if ( present( sclrprtp ) ) then
             write(fstderr,*) "sclrprtp = ", sclrprtp
          endif
          
          if ( present( sclrpthlp ) )then
             write(fstderr,*) "sclrpthlp = ", sclrpthlp
          endif

          write(fstderr,*) "Intent(out)"
          
          write(fstderr,*) "wp4 = ", wp4
          write(fstderr,*) "wprtp2 = ", wprtp2
          write(fstderr,*) "wp2rtp = ", wp2rtp
          write(fstderr,*) "wpthlp2 = ", wpthlp2
          write(fstderr,*) "cf = ", cf
          write(fstderr,*) "rcm = ", rcm
          write(fstderr,*) "wpthvp = ", wpthvp
          write(fstderr,*) "wp2thvp = ", wp2thvp
          write(fstderr,*) "rtpthvp = ", rtpthvp
          write(fstderr,*) "thlpthvp = ", thlpthvp
          write(fstderr,*) "wprcp = ", wprcp
          write(fstderr,*) "wp2rcp = ", wp2rcp
          write(fstderr,*) "rtprcp = ", rtprcp
          write(fstderr,*) "thlprcp = ", thlprcp
          write(fstderr,*) "rcp2 = ", rcp2
          write(fstderr,*) "wprtpthlp = ", wprtpthlp
          write(fstderr,*) "crt1 = ", crt1
          write(fstderr,*) "crt2 = ", crt2
          write(fstderr,*) "cthl1 = ", cthl1
          write(fstderr,*) "cthl2 = ", cthl2
          write(fstderr,*) "pdf_parms = ", pdf_parms
          
          if ( present( sclrpthvp ) )then
             write(fstderr,*) "sclrpthvp = ", sclrpthvp
          endif
          
          if ( present( sclrprcp ) )then
             write(fstderr,*) "sclrprcp = ", sclrprcp
          endif
          
          if ( present( wpsclrp2 ) )then
             write(fstderr,*) "wpsclrp2 = ", wpsclrp2
          endif
          
          if ( present( wpsclrprtp ) )then
             write(fstderr,*) "wpsclrprtp = ", wpsclrprtp
          endif
          
          if ( present( wpsclrpthlp ) )then
             write(fstderr,*) "wpsclrpthlp = ", wpsclrpthlp
          endif

          if ( present( wp2sclrp ) )then
             write(fstderr,*) "wp2sclrp = ", wp2sclrp
          endif
          
        endif
        
        return
        
        end subroutine pdf_closure_new

        end module pdf_closure
