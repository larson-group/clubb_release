!------------------------------------------------------------------------
!
! Subroutine that computes pdf parameters analytically.
!
! Based of the original formulation, but with some tweaks
! to remove some of the less realistic assumptions and
! improve transport terms.
!
! Corrected version that should remove inconsistency
!
        subroutine pdf_closure_new(p,exner,
     .                         wm,wp2,wp3,Sc,
     .                         rtm,rtp2,wprtp,
     .                         thlm,thlp2,wpthlp,
     .                         wp4,
     .                         wprtp2,wp2rtp,
     .                         wpthlp2,wp2thlp,
     .                         rtpthlp,
     .                         wprtpthlp,
     .                         cf,rcm,
     .                         wpthvp,wp2thvp,rtpthvp,thlpthvp,
     .                         wprlp,wp2rlp,rtprlp,thlprlp,
     .                         pdf_parms,
     .                         crt1,crt2,cthl1,cthl2)

        use constants

        implicit none

! External

        real erf
        real rsat
        external rsat

! Input

        real, intent(in) :: p, exner
        real, intent(in) :: wm, wp2, wp3, Sc
        real, intent(in) :: rtm, rtp2, wprtp
        real, intent(in) :: thlm, thlp2, wpthlp
        real, intent(in) :: rtpthlp

! Output

        real, intent(out) :: wp4
        real, intent(out) :: wprtp2, wp2rtp, wpthlp2, wp2thlp
        real, intent(out) :: cf, rcm
        real, intent(out) :: wpthvp, wp2thvp, rtpthvp, thlpthvp
        real, intent(out) :: wprlp, wp2rlp, rtprlp, thlprlp
        real, intent(out) :: wprtpthlp
        real, intent(out) :: pdf_parms(26)
        real, intent(out) :: crt1, crt2, cthl1, cthl2

! Internal 

        real :: a
        real :: w1, w2, sw1, sw2
        real :: thl1, thl2, sthl1, sthl2
        real :: rt1, rt2, srt1, srt2
        real :: w1_n, w2_n, thl1_n, thl2_n, rt1_n, rt2_n

!        quantities needed to predict higher order moments

        real :: tl1, tl2, beta1, beta2, rsl1, rsl2, es1, es2
        real :: ss1, ss2, s1, s2, zeta1, zeta2
        real :: rl1, rl2, R1, R2

!        sub-plume correlation coefficient between rt, thl
!        varies between -1 < rrtthl < 1

        real :: rrtthl

!        thermodynamic quantity

        real :: bd

!        skewness of w

        real :: Skw 

!       weighting factor for interpolation of srt1,2, sthl1,2

        real :: diffweight

!       alpha coefficients

        real :: alpha_thl, alpha_rt

!       variables for a generalization of Chris Golaz' closure
!       varies width of plumes in theta_l, rt
        real :: width_factor_1, width_factor_2



! Code begins -------------------------------------------

        BD = Lv / (exner*cp) - ep2 * T0

!        if no velocity, then use single delta fnc. as pdf
!        else compute parameters explicitly

        if ( wp2 <= wtol**2 )  then

          a = 0.5
          w1 = wm
          w2 = wm
          sw1 = 0.
          sw2 = 0.
          rt1 = rtm
          rt2 = rtm
          srt1 = 0.
          srt2 = 0.
          thl1 = thlm
          thl2 = thlm
          sthl1 = 0.
          sthl2 = 0.
          rrtthl = 0.

        else

          Skw = wp3/wp2**1.5

          if ( abs(Skw) <= 1e-5 ) then
            a = 0.5
          else
            a = 0.5
     .          * ( 
     .              1. - Skw/sqrt(4.*(1.-Sc)**3+Skw**2)
     .            )
          end if

!          clip a, 1-a, to avoid dividing by zero
          if (a > 0.99) a = 0.99
          if (a < 0.01) a = 0.01

          w1_n = sqrt(((1-a)/a)*(1-Sc))
          w2_n = -sqrt((a/(1-a))*(1-Sc))
          w1 = wm + sqrt(wp2)*w1_n
          w2 = wm + sqrt(wp2)*w2_n

          sw1 = Sc*wp2
          sw2 = Sc*wp2

!         Vince Larson added a factor so that width of plumes in 
!         theta_l, rt can vary.  
!         beta is a constant defined in constants.f.  Set 0<beta<3.
!         beta=1.5 recovers Chris Golaz' simplified formula.
!         3 Nov 2003

          width_factor_1 = (2./3.)*beta + 2.*a*(1 - (2./3.)*beta)
          width_factor_2 = 2. - width_factor_1

          if ( thlp2 <= thltol**2 ) then
            thl1 = thlm
            thl2 = thlm
            sthl1 = 0.
            sthl2 = 0.
            alpha_thl = 0.5
          else
            thl1_n = - (wpthlp/(sqrt(wp2)*sqrt(thlp2)))/w2_n
            thl2_n = - (wpthlp/(sqrt(wp2)*sqrt(thlp2)))/w1_n

            thl1 = thlm - (wpthlp/sqrt(wp2))/w2_n
            thl2 = thlm - (wpthlp/sqrt(wp2))/w1_n

            alpha_thl
     .      = 0.5
     .        * ( 1.0 - wpthlp*wpthlp /((1.0-Sc)*wp2*thlp2) )
            alpha_thl = max( min( alpha_thl, 1.0 ), 0.0 )

! Vince Larson multiplied original expressions by width_factor_1,2
!  to generalize scalar skewnesses.  05 Nov 03
            sthl1 = ( alpha_thl / a * thlp2 ) * width_factor_1
            sthl2 = ( alpha_thl / (1.-a) * thlp2 ) * width_factor_2

          end if

          if ( rtp2 <= rttol**2) then
            rt1 = rtm
            rt2 = rtm
            srt1 = 0.
            srt2 = 0.
            alpha_rt = 0.5
          else
            rt1_n = - (wprtp/(sqrt(wp2)*sqrt(rtp2)))/w2_n
            rt2_n = - (wprtp/(sqrt(wp2)*sqrt(rtp2)))/w1_n

            rt1 = rtm - (wprtp/sqrt(wp2))/w2_n
            rt2 = rtm - (wprtp/sqrt(wp2))/w1_n

            alpha_rt
     .      = 0.5
     .        * ( 1.0 - wprtp*wprtp /((1.0-Sc)*wp2*rtp2) )
            alpha_rt = max( min( alpha_rt, 1.0 ), 0.0 )

! Vince Larson multiplied original expressions by width_factor_1,2
!    to generalize scalar skewnesses.  05 Nov 03
            srt1 = ( alpha_rt / a * rtp2 ) * width_factor_1
            srt2 = ( alpha_rt / (1.-a) * rtp2 ) * width_factor_2

        end if

!       We include sub-plume correlation with coeff rrtthl.
!       We assume srt1=srt2, sthl1=sthl2, rrtthl is a fixed const.

          if ( srt1*sthl1 > 0 .and. srt2*sthl2 > 0) then
              rrtthl = ( rtpthlp 
     .             - a * ( rt1-rtm ) * ( thl1-thlm )
     .             - (1-a) * ( rt2-rtm ) * ( thl2-thlm ) )
     .                 / ( a*sqrt(srt1*sthl1) + (1-a)*sqrt(srt2*sthl2) )
                if (rrtthl < -1.) then
                        rrtthl = -1.
                endif
                if (rrtthl > 1.) then
                        rrtthl = 1.
                endif
          else
                rrtthl = 0.
          end if

        end if

!       higher order moments

        wp2rtp = a * ( (w1-wm)**2+sw1 ) * ( rt1-rtm )
     .           + (1-a) * ( (w2-wm)**2+sw2 ) * ( rt2-rtm )

        wp2thlp = a * ( (w1-wm)**2+sw1 ) * ( thl1-thlm )
     .            + (1-a) * ( (w2-wm)**2+sw2 ) * ( thl2-thlm )

        wp4 = a * ( 3*sw1**2 + 6*((w1-wm)**2)*sw1 + (w1-wm)**4 )
     .        + (1-a) * ( 3*sw2**2 + 6*((w2-wm)**2)*sw2 + (w2-wm)**4 )

        wprtp2 = a * ( w1-wm )*( (rt1-rtm)**2 + srt1 ) 
     .           + (1-a) * ( w2-wm )*( (rt2-rtm)**2 + srt2)

        wpthlp2 = a * ( w1-wm )*( (thl1-thlm)**2 + sthl1 ) 
     .            + (1-a) * ( w2-wm )*( (thl2-thlm)**2+sthl2 )

        wprtpthlp = a * ( w1-wm )*( (rt1-rtm)*(thl1-thlm) 
     .                + rrtthl*sqrt(srt1*sthl1) )
     .            + ( 1-a ) * ( w2-wm )*( (rt2-rtm)*(thl2-thlm)
     .                + rrtthl*sqrt(srt2*sthl2) )


!       Compute higher order moments that include theta_v.

!       First compute some preliminary quantities.
!       "1" denotes first Gaussian; "2" denotes 2nd Gaussian
!       liq water temp (Sommeria & Deardorff 1977 (SD), eqn. 3)

        tl1 = thl1*exner
        tl2 = thl2*exner

        rsl1 = rsat(p,tl1)
        rsl2 = rsat(p,tl2)

!       SD's beta (eqn. 8)
        beta1 = ep*(Lv/(R*tl1))*(Lv/(cp*tl1))
        beta2 = ep*(Lv/(R*tl2))*(Lv/(cp*tl2))

!       s from Lewellen and Yoh 1993 (LY) eqn. 1
        s1 = (rt1-rsl1)/(1+beta1*rsl1)
        s2 = (rt2-rsl2)/(1+beta2*rsl2)

!       coefficient for s'
        crt1 = 1/(1+beta1*rsl1)
        crt2 = 1/(1+beta2*rsl2)

        cthl1 = ((1+beta1*rt1)/(1+beta1*rsl1)**2)*(cp/Lv)
     .        *beta1*rsl1*exner
        cthl2 = ((1+beta2*rt2)/(1+beta2*rsl2)**2)*(cp/Lv)
     .        *beta2*rsl2*exner

!     standard deviation of s
!     include subplume correlation of qt, thl
!     Because of round-off error, 
!     ss1 (and probably ss2) can become negative when rrtthl=1
!     One could also write this as a squared term
!     plus a postive correction; this might be a neater format

        ss1 = sqrt(max(0.,(srt1*crt1**2 + sthl1*cthl1**2 
     .          - 2.*rrtthl*crt1*sqrt(srt1)*cthl1*sqrt(sthl1)) ) )
        ss2 = sqrt(max(0.,(srt2*crt2**2 + sthl2*cthl2**2
     .          - 2.*rrtthl*crt2*sqrt(srt2)*cthl2*sqrt(sthl2)) ) )

!        ss1 = sqrt( (sqrt(srt1)*crt1 - sqrt(sthl1)*cthl1)**2 
!     .          + (1.-rrtthl)*2.*crt1*sqrt(srt1)*cthl1*sqrt(sthl1)  )
!        ss2 = sqrt( (sqrt(srt2)*crt2 - sqrt(sthl2)*cthl2)**2
!     .          + (1.-rrtthl)*2.*crt2*sqrt(srt2)*cthl2*sqrt(sthl2)  )


!       we need to introduce a threshold value for the variance of s

        if ( ss1 > sstol ) then
           zeta1 = s1/ss1
           R1 = 0.5*(1+erf(zeta1/sqrt(2.0)))
           rl1 =s1*R1+ss1*exp(-0.5*zeta1**2)/(sqrt(2*pi))
        else
           if ( s1 < 0. ) then
              R1 = 0.
              rl1 = 0.
           else
              R1 = 1.
              rl1 = s1
           end if
        end if
        if ( ss2 > sstol ) then
           zeta2 = s2/ss2
           R2 = 0.5*(1+erf(zeta2/sqrt(2.0)))
           rl2 =s2*R2+ss2*exp(-0.5*zeta2**2)/(sqrt(2*pi))
        else
           if ( s2 < 0. ) then
             R2 = 0.
             rl2 = 0.
           else
             R2 = 1.
             rl2 = s2
           end if
        end if

!       now compute moments that depend on theta_v

        wp2rlp = a * ((w1-wm)**2 + sw1)*rl1
     .                 + (1-a) * ((w2-wm)**2 + sw2)*rl2
     .                 - wp2 * (a*rl1+(1-a)*rl2)
        wp2thvp = wp2thlp + ep1*T0*wp2rtp + BD*wp2rlp

        wprlp = a * (w1-wm)*rl1 + (1-a) * (w2-wm)*rl2
        wpthvp = wpthlp + ep1*T0*wprtp + BD*wprlp

!       account for subplume correlation in qt-thl
        thlprlp = a * ( (thl1-thlm)*rl1 - (cthl1*sthl1)*R1 )
     .            + (1-a) * ( (thl2-thlm)*rl2 - (cthl2*sthl2)*R2 )
     .            + a*rrtthl*crt1*sqrt(srt1)*sqrt(sthl1)*R1
     .            + (1-a)*rrtthl*crt2*sqrt(srt2)*sqrt(sthl2)*R2
        thlpthvp = thlp2 + ep1*T0*rtpthlp + BD*thlprlp

!       account for subplume correlation in qt-thl
        rtprlp = a * ( (rt1-rtm)*rl1 + (crt1*srt1)*R1 )
     .                 + (1-a) * ( (rt2-rtm)*rl2 + (crt2*srt2)*R2 )
     .            - a*rrtthl*cthl1*sqrt(srt1)*sqrt(sthl1)*R1
     .            - (1-a)*rrtthl*cthl2*sqrt(srt2)*sqrt(sthl2)*R2

        rtpthvp = rtpthlp + ep1*T0*rtp2 + BD*rtprlp

!       compute mean cloud fraction and cloud water

        cf = a * R1 + (1-a) * R2
        rcm = a * rl1 + (1-a) * rl2

        ! Note: Brian added the following line to ensure that there
        ! are never any negative liquid water values.  According to
        ! Vince Larson, the analytic formula should not produce any
        ! negative results, but such computer-induced errors such as
        ! round-off error may produce such a value.  This has been
        ! corrected because Brian found a small negative value of
        ! rcm in the first timestep of the FIRE case.

        rcm = MAX(0.0, rcm)

!       save pdf parameters

        pdf_parms(1) = w1
        pdf_parms(2) = w2
        pdf_parms(3) = sw1
        pdf_parms(4) = sw2
        pdf_parms(5) = rt1
        pdf_parms(6) = rt2
        pdf_parms(7) = srt1
        pdf_parms(8) = srt2
        pdf_parms(9) = thl1
        pdf_parms(10) = thl2
        pdf_parms(11) = sthl1
        pdf_parms(12) = sthl2
        pdf_parms(13) = a
        pdf_parms(14) = rl1
        pdf_parms(15) = rl2
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

        return
        end subroutine pdf_closure_new
