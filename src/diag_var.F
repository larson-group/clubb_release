! $Id: diag_var.F,v 1.22 2006-04-14 20:08:07 dschanen Exp $
#define SCLR_THETA 1
#define SCLR_RT 2
!-----------------------------------------------------------------------
!  This is a module because optional variables require it to be.
!-----------------------------------------------------------------------
      module diagnose_variances

        implicit none

        public diag_var

        private calc_aa, calc_bb, calc_cc, calc_dd_var, calc_dd_covar

        contains
!-----------------------------------------------------------------------
! subroutine diag_var( )
! Subprogram to diagnose variances by solving steady-state equations
!-----------------------------------------------------------------------
        subroutine diag_var( tau, wmm, rtm, wprtp, thlm, wpthlp, 
     .                       wp2, wp3, Sc,
     .                       rtp2, thlp2, rtpthlp,
     .                       liter, dt, isValid,
     .                       sclrm, wpsclrp, sclrp2,        ! optional
     .                       sclrprtp, sclrpthlp )          ! optional

        use constants
        use grid_class
#ifdef STATS
        use statistics
#endif
       
        implicit none

! Input variables

        real, intent(in) :: tau(1:gr%nnzp)
        real, intent(in) :: wmm(1:gr%nnzp)
        real, intent(in) :: rtm(1:gr%nnzp)
        real, intent(in) :: wprtp(1:gr%nnzp)
        real, intent(in) :: thlm(1:gr%nnzp)
        real, intent(in) :: wpthlp(1:gr%nnzp)
        real, intent(in) :: wp2(1:gr%nnzp)
        real, intent(in) :: wp3(1:gr%nnzp)
        real, intent(in) :: Sc(1:gr%nnzp)

        logical, intent(in) :: liter
        real, intent(in)    :: dt


! Output variables
!  note: if liter=.true., 
!  then these would need to have an attrib. of type (inout)
!     An attribute of (inout) is also needed to import the
!     value of the variances at the surface.  Brian.  12/18/05.

       real, intent(inout) :: rtp2(1:gr%nnzp)
       real, intent(inout) :: thlp2(1:gr%nnzp)
       real, intent(inout) :: rtpthlp(1:gr%nnzp)

! Input mixing scheme variables (optional)

       real, optional, intent(in), 
     . dimension(gr%nnzp, sclrm_dimension) :: sclrm, wpsclrp

! Output mixing scheme variables (optional)

       real, optional, intent(out), 
     . dimension(gr%nnzp, sclrm_dimension):: sclrp2, sclrprtp, sclrpthlp

c Local constants

!       real, parameter :: eps = 1.e-10  ! now in constants.F
c Misc. Variables -D.S.
       logical, intent(inout) :: isValid

c Local arrays allocated dynamically on the stack

       real a1(1:gr%nnzp)
       real aa(1:gr%nnzp), bb(1:gr%nnzp), cc(1:gr%nnzp), dd(1:gr%nnzp)

c Loop indices

       integer i, k
       integer km1, kp1

#ifdef STATS
       if (lstats_samp) then

         if ( irtp2_bt > 0 ) then
           zm%x(:,irtp2_bt) = zm%x(:,irtp2_bt) - rtp2 / dt
         endif

         if ( ithlp2_bt > 0 ) then
           zm%x(:,ithlp2_bt) = zm%x(:,ithlp2_bt) - thlp2 / dt
         endif

         if ( irtpthlp_bt > 0 ) then
           zm%x(:,irtpthlp_bt) = zm%x(:,irtpthlp_bt) - rtpthlp / dt
         endif

       endif
#endif

c Define a1

       a1 = 1.0 / (1.0-Sc)

c Solve for rtp2

c      prepare tridiagonal system

       aa = calc_aa( a1, wp2, wp3, wmm, nu2, beta, eps )
       bb = calc_bb( a1, wp2, wp3, tau, 
     .               C2rt, dt, nu2, beta, eps, liter )
       cc = calc_cc( a1, wp2, wp3, wmm, nu2, beta, eps )
       dd = calc_dd_var( a1, wp2, wp3, wprtp, rtm, 
     .                   rtp2, dt, beta, eps, liter )

#ifdef STATS
       if (lstats_samp) then

         do k=1, gr%nnzp

           km1   = max(k-1,1)
           kp1   = min(k+1,gr%nnzp)

           if ( irtp2_dp1 > 0 ) then
             zmscr01(k) = - C2rt / tau(k)
           endif
         
           if ( irtp2_dp2 * ithlp2_dp2 * irtpthlp_dp2 > 0 ) then
             zmscr02(k) = + nu2 * gr%dzm(k) * gr%dzt(k)
             zmscr03(k) = 
     .       - nu2 * gr%dzm(k) * ( gr%dzt(kp1) + gr%dzt(k) )
             zmscr04(k) = + nu2 * gr%dzm(k) * gr%dzt(kp1)
           endif

           if ( irtp2_ta * ithlp2_ta * irtpthlp_ta > 0 ) then
             zmscr05(k) = 
     .       + (1.0/3.0) * beta
     .         * (1.0/2.0) * gr%dzm(k)
     .            * ( a1(k) + a1(km1) )
     .            * ( wp3(k) / MAX( wp2(k) + wp2(km1), 2.0*eps ) )
             zmscr06(k) = 
     .       - (1.0/3.0) * beta
     .         * (1.0/2.0) * gr%dzm(k)
     .            * (   ( a1(kp1) + a1(k) )
     .                * ( wp3(kp1) / MAX( wp2(kp1) + wp2(k), 2.0*eps ) )
     .               -  ( a1(k) + a1(km1) )
     .                * ( wp3(k) / MAX( wp2(k) + wp2(km1), 2.0*eps ) )
     .              )
             zmscr07(k) = 
     .       - (1.0/3.0) * beta
     .         * (1.0/2.0) * gr%dzm(k)
     .            * ( a1(kp1) + a1(k) )
     .            * ( wp3(kp1) / MAX( wp2(kp1) + wp2(k), 2.0*eps ) )
           endif

           if ( irtp2_ma * ithlp2_ma * irtpthlp_ma > 0 ) then
             zmscr08(k) = + wmm(k) * (1.0/2.0) * gr%dzm(k)
             zmscr09(k) = - wmm(k) * (1.0/2.0) * gr%dzm(k)
           endif

           if ( irtp2_ta > 0 ) then
             zm%x(k,irtp2_ta) 
     .       = zm%x(k,irtp2_ta)
     .       - ( 1.0 - (1.0/3.0) * beta )
     .         * (1.0/4.0) * gr%dzm(k)
     .            * (   ( a1(kp1) + a1(k) )**2.0
     .                  * wp3(kp1) * ( wprtp(kp1) + wprtp(k) )**2.0
     .                  / ( MAX( wp2(kp1) + wp2(k), 2.0*eps ) )**2.0
     .               -  ( a1(k) + a1(km1) )**2.0
     .                  * wp3(k) * ( wprtp(k) + wprtp(km1) )**2.0
     .                  / ( MAX( wp2(k) + wp2(km1), 2.0*eps ) )**2.0
     .              )
           endif

         enddo

       endif
#endif /*STATS*/
       
c      solve tridiagonal system
       call tridag( aa, bb, cc, dd, rtp2, gr%nnzp, isValid )
       if ( .not. isValid ) then
          print *, "diag_var: rtp2 tridag failed"
          return
       endif

#ifdef STATS
       if (lstats_samp) then

          do k=1, gr%nnzp

           km1 = max(k-1,1)
           kp1 = min(k+1,gr%nnzp)

           if ( irtp2_dp1 > 0 ) then
             zm%x(k,irtp2_dp1) = zm%x(k,irtp2_dp1)
     .       + zmscr01(k) * rtp2(k)
             zm%n(k,irtp2_dp1) = zm%n(k,irtp2_dp1) + 1
           endif
 
           if ( irtp2_dp2 > 0 ) then
             zm%x(k,irtp2_dp2) = zm%x(k,irtp2_dp2)
     .       + zmscr02(k) * rtp2(km1)
     .       + zmscr03(k) * rtp2(k)
     .       + zmscr04(k) * rtp2(kp1)
             zm%n(k,irtp2_dp2) = zm%n(k,irtp2_dp2) + 1
           endif
 
           if ( irtp2_ta > 0 ) then
             zm%x(k,irtp2_ta) = zm%x(k,irtp2_ta)
     .       + zmscr05(k) * rtp2(km1)
     .       + zmscr06(k) * rtp2(k)
     .       + zmscr07(k) * rtp2(kp1)
             zm%n(k,irtp2_ta) = zm%n(k,irtp2_ta) + 1
           endif

           if ( irtp2_ma > 0 ) then
             zm%x(k,irtp2_ma) = zm%x(k,irtp2_ma)
     .       + zmscr08(k) * rtp2(km1)
     .       + zmscr09(k) * rtp2(kp1)
             zm%n(k,irtp2_ma) = zm%n(k,irtp2_ma) + 1
           endif
 
           if ( irtp2_tp > 0 ) then
             zm%x(k,irtp2_tp) = zm%x(k,irtp2_tp)
     .        - 2.0 * wprtp(k) * gr%dzm(k) * ( rtm(kp1) - rtm(k) )
             zm%n(k,irtp2_tp) = zm%n(k,irtp2_tp) + 1
           endif

          enddo

       endif
#endif /*STATS*/

#ifdef STATS
       if (lstats_samp .and. irtp2_cl > 0 ) then
         zm%x(:,irtp2_cl) = zm%x(:,irtp2_cl) - rtp2 / dt
       endif
#endif

c      clip small and negative values
       ! The surface variances are set elsewhere.  They do not need to be
       ! clipped at the lowest level.  Brian Griffin.  12/18/05.
       do k = 2, gr%nnzp, 1

         if ( wp2(k) >= wtol*wtol ) then
           rtp2(k) 
     .     = max( rtp2(k), 
     .            wprtp(k)*wprtp(k)/((1.0-Sc(k))*wp2(k)) + eps,
     .            rttol*rttol + eps )
         else
           rtp2(k) = max( rtp2(k), 0.0 )
         endif

       enddo

#ifdef STATS
       if (lstats_samp .and. irtp2_cl > 0 ) then
         zm%x(:,irtp2_cl) = zm%x(:,irtp2_cl) + rtp2 / dt
         zm%n(:,irtp2_cl) = zm%n(:,irtp2_cl) + 1
       endif
#endif

c Solve for thlp2

c      prepare tridiagonal system
c      note: we can reuse the same aa & cc arrays

        bb = calc_bb( a1, wp2, wp3, tau, 
     .                C2thl, dt, nu2, beta, eps, liter )
        dd = calc_dd_var( a1, wp2, wp3, wpthlp, thlm, 
     .                    thlp2, dt, beta, eps, liter )

#ifdef STATS
       if (lstats_samp) then

         do k=1,gr%nnzp

           km1 = max(k-1,1)
           kp1 = min(k+1,gr%nnzp)

           if ( ithlp2_dp1 > 0 ) then
             zmscr01(k) = - C2thl / tau(k)
           endif
         
           if ( ithlp2_ta > 0 ) then
             zm%x(k,ithlp2_ta) 
     .       = zm%x(k,ithlp2_ta)
     .       - ( 1.0 - (1.0/3.0) * beta )
     .         * (1.0/4.0) * gr%dzm(k)
     .            * (   ( a1(kp1) + a1(k) )**2.0
     .                  * wp3(kp1) * ( wpthlp(kp1) + wpthlp(k) )**2.0
     .                  / ( MAX( wp2(kp1) + wp2(k), 2.0*eps ) )**2.0
     .               -  ( a1(k) + a1(km1) )**2.0
     .                  * wp3(k) * ( wpthlp(k) + wpthlp(km1) )**2.0
     .                  / ( MAX( wp2(k) + wp2(km1), 2.0*eps ) )**2.0
     .              )
           endif

         enddo ! k

       endif
#endif
       
c      solve tridiagonal system
       call tridag( aa, bb, cc, dd, thlp2, gr%nnzp, isValid )
       if ( .not. isValid ) then
          print *, "diag_var: thlp2 tridag failed"
          return
       endif

#ifdef STATS
       if (lstats_samp) then

         do k=1, gr%nnzp

           km1 = max(k-1,1)
           kp1 = min(k+1,gr%nnzp)

           if ( ithlp2_dp1 > 0 ) then
             zm%x(k,ithlp2_dp1) = zm%x(k,ithlp2_dp1)
     .       + zmscr01(k) * thlp2(k)
             zm%n(k,ithlp2_dp1) = zm%n(k,ithlp2_dp1) + 1
           endif
 
           if ( ithlp2_dp2 > 0 ) then
             zm%x(k,ithlp2_dp2) = zm%x(k,ithlp2_dp2)
     .       + zmscr02(k) * thlp2(km1)
     .       + zmscr03(k) * thlp2(k)
     .       + zmscr04(k) * thlp2(kp1)
             zm%n(k,ithlp2_dp2) = zm%n(k,ithlp2_dp2) + 1
           endif
 
           if ( ithlp2_ta > 0 ) then
             zm%x(k,ithlp2_ta) = zm%x(k,ithlp2_ta)
     .       + zmscr05(k) * thlp2(km1)
     .       + zmscr06(k) * thlp2(k)
     .       + zmscr07(k) * thlp2(kp1)
             zm%n(k,ithlp2_ta) = zm%n(k,ithlp2_ta) + 1
           endif

           if ( ithlp2_ma > 0 ) then
             zm%x(k,ithlp2_ma) = zm%x(k,ithlp2_ma)
     .       + zmscr08(k) * thlp2(km1)
     .       + zmscr09(k) * thlp2(kp1)
             zm%n(k,ithlp2_ma) = zm%n(k,ithlp2_ma) + 1
           endif
 
           if ( ithlp2_tp > 0 ) then
             zm%x(k,ithlp2_tp) = zm%x(k,ithlp2_tp)
     .        - 2.0 * wpthlp(k) * gr%dzm(k) * ( thlm(kp1) - thlm(k) )
             zm%n(k,ithlp2_tp) = zm%n(k,ithlp2_tp) + 1

           endif

         enddo ! k

       endif
#endif

#ifdef STATS
       if (lstats_samp .and. ithlp2_cl > 0 ) then
         zm%x(:,ithlp2_cl) = zm%x(:,ithlp2_cl) - thlp2 / dt
       endif
#endif


c      clip small and negative values
       ! The surface variances are set elsewhere.  They do not need to be
       ! clipped at the lowest level.  Brian Griffin.  12/18/05.
       do k = 2, gr%nnzp, 1

         if ( wp2(k) >= wtol*wtol ) then
           thlp2(k) 
     .     = max( thlp2(k), 
     .            wpthlp(k)*wpthlp(k)/((1.0-Sc(k))*wp2(k)) + eps,
     .            thltol*thltol + eps )
         else
           thlp2(k) = max( thlp2(k), 0.0 )
         endif

       enddo ! k

#ifdef STATS
       if (lstats_samp .and. ithlp2_cl > 0 ) then
         zm%x(:,ithlp2_cl) = zm%x(:,ithlp2_cl) + thlp2 / dt
         zm%n(:,ithlp2_cl) = zm%n(:,ithlp2_cl) + 1
       endif
#endif

c Solve for sclrp2

c      prepare tridiagonal system
c      note: we can reuse the same aa, cc arrays

      if (     present( sclrm )
     .   .and. present( wpsclrp ) 
     .   .and. present( sclrp2 ) ) then

       do i=1, sclrm_dimension

         if ( i == SCLR_THETA ) then
           bb = calc_bb( a1, wp2, wp3, tau,
     .                   C2thl, dt, nu2, beta, eps, liter )
         elseif ( i == SCLR_RT ) then
           bb = calc_bb( a1, wp2, wp3, tau,
     .                   C2rt, dt, nu2, beta, eps, liter )
         endif

         dd = calc_dd_var( a1, wp2, wp3, wpsclrp(:,i), sclrm(:,i), 
     .                     sclrp2(:,i), dt, beta, eps, liter )

c        solve tridiagonal system
         call tridag( aa, bb, cc, dd, sclrp2(:,i), gr%nnzp, isValid )
         if ( .not. isValid ) then
           print *, "diag_var: sclrp2 tridag failed"
           return
         endif

c        clip small and negative values
         ! The surface variances are set elsewhere.  They do not need to be
         ! clipped at the lowest level.  Brian Griffin.  12/18/05.
         do k = 2, gr%nnzp, 1

           if ( wp2(k) >= wtol*wtol ) then
             sclrp2(k,i) 
     .       = max( sclrp2(k,i), 
     .              wpsclrp(k,i)*wpsclrp(k,i)/((1.0-Sc(k))*wp2(k))+eps,
     .              sclrtol(i)*sclrtol(i) + eps )
           else
             sclrp2(k,i) = max( sclrp2(k,i), 0.0 )
           endif

         enddo ! k

       enddo ! i

      endif ! if present( sclrs )

! %% Debug
!      sclrp2(:,2) = rtp2
! %% Debug

c Solve for rtpthlp

c      prepare tridiagonal system

       bb = calc_bb( a1, wp2, wp3, tau, 
     .               C2rtthl, dt, nu2, beta, eps, liter )

       dd = calc_dd_covar( a1, wp2, wp3, wprtp, wpthlp, 
     .                     rtm, thlm, rtpthlp, 
     .                     dt, beta, eps, liter ) 

#ifdef STATS

       if (lstats_samp) then

          do k=1, gr%nnzp

             km1 = max(k-1,1)
             kp1 = min(k+1,gr%nnzp)

             if ( irtpthlp_dp1 > 0 ) then
               zmscr01(k) = - C2rtthl / tau(k)
             endif
         
             if ( irtpthlp_ta > 0 ) then
             zm%x(k,irtpthlp_ta) 
     .       = zm%x(k,irtpthlp_ta)
     .       - ( 1.0 - (1.0/3.0) * beta )
     .         * (1.0/4.0) * gr%dzm(k)
     .            * (  ( a1(kp1) + a1(k) )**2.0
     .                 * ( wprtp(kp1) + wprtp(k) ) 
     .                 * ( wpthlp(kp1) + wpthlp(k) )
     .                 * wp3(kp1) 
     .                 / ( MAX( wp2(kp1) + wp2(k), 2.0*eps ) )**2.0
     .               - ( a1(k) + a1(km1) )**2.0
     .                 * ( wprtp(k) + wprtp(km1) ) 
     .                 * ( wpthlp(k) + wpthlp(km1) )
     .                 * wp3(k) 
     .                 / ( MAX( wp2(k) + wp2(km1), 2.0*eps ) )**2.0
     .              )
             endif

          enddo ! k

       endif
#endif
       
c      solve tridiagonal system
       call tridag( aa, bb, cc, dd, rtpthlp, gr%nnzp, isValid )
       if ( .not. isValid ) then
          print *, "diag_var: rtpthlp tridag failed"
          return
       endif

#ifdef STATS
       if (lstats_samp) then

         do k=1, gr%nnzp

           km1 = max(k-1,1)
           kp1 = min(k+1,gr%nnzp)

           if ( irtpthlp_dp1 > 0 ) then
             zm%x(k,irtpthlp_dp1) = zm%x(k,irtpthlp_dp1)
     .       + zmscr01(k) * rtpthlp(k)
             zm%n(k,irtpthlp_dp1) = zm%n(k,irtpthlp_dp1) + 1
           endif
 
           if ( irtpthlp_dp2 > 0 ) then
             zm%x(k,irtpthlp_dp2) = zm%x(k,irtpthlp_dp2)
     .       + zmscr02(k) * rtpthlp(km1)
     .       + zmscr03(k) * rtpthlp(k)
     .       + zmscr04(k) * rtpthlp(kp1)
             zm%n(k,irtpthlp_dp2) = zm%n(k,irtpthlp_dp2) + 1
           endif
 
           if ( irtpthlp_ta > 0 ) then
             zm%x(k,irtpthlp_ta) = zm%x(k,irtpthlp_ta)
     .       + zmscr05(k) * rtpthlp(km1)
     .       + zmscr06(k) * rtpthlp(k)
     .       + zmscr07(k) * rtpthlp(kp1)
             zm%n(k,irtpthlp_ta) = zm%n(k,irtpthlp_ta) + 1
           endif

           if ( irtpthlp_ma > 0 ) then
             zm%x(k,irtpthlp_ma) = zm%x(k,irtpthlp_ma)
     .       + zmscr08(k) * rtpthlp(km1)
     .       + zmscr09(k) * rtpthlp(kp1)
             zm%n(k,irtpthlp_ma) = zm%n(k,irtpthlp_ma) + 1
           endif
 
           if ( irtpthlp_tp1 > 0 ) then
             zm%x(k,irtpthlp_tp1) = zm%x(k,irtpthlp_tp1)
     .       - wprtp(k) * gr%dzm(k) * ( thlm(kp1) - thlm(k) )
             zm%n(k,irtpthlp_tp1) = zm%n(k,irtpthlp_tp1) + 1
           endif

           if ( irtpthlp_tp2 > 0 ) then
             zm%x(k,irtpthlp_tp2) = zm%x(k,irtpthlp_tp2)
     .       - wpthlp(k) * gr%dzm(k) * ( rtm(kp1) - rtm(k) )
             zm%n(k,irtpthlp_tp2) = zm%n(k,irtpthlp_tp2) + 1
           endif

         enddo

       endif
#endif /*STATS*/

#ifdef STATS
       if (lstats_samp) then

         if ( irtp2_bt > 0 ) then
           zm%x(:,irtp2_bt) = zm%x(:,irtp2_bt) + rtp2 / dt
           zm%n(:,irtp2_bt) = zm%n(:,irtp2_bt) + 1
         endif

         if ( ithlp2_bt > 0 ) then
           zm%x(:,ithlp2_bt) = zm%x(:,ithlp2_bt) + thlp2 / dt
           zm%n(:,ithlp2_bt) = zm%n(:,ithlp2_bt) + 1
         endif

         if ( irtpthlp_bt > 0 ) then
           zm%x(:,irtpthlp_bt) = zm%x(:,irtpthlp_bt) + rtpthlp / dt
           zm%n(:,irtpthlp_bt) = zm%n(:,irtpthlp_bt) + 1
         endif

       endif
#endif

c Solve for sclrprtp
      if (     present( sclrm )
     .   .and. present( wpsclrp ) 
     .   .and. present( sclrprtp ) ) then

c     prepare tridiagonal system
       do i=1, sclrm_dimension

         if ( i == SCLR_THETA ) then
           bb = calc_bb( a1, wp2, wp3, tau, 
     .                   C2rtthl, dt, nu2, beta, eps, liter )
         else 
           if ( i == SCLR_RT ) then
             bb = calc_bb( a1, wp2, wp3, tau, 
     .                     C2rt, dt, nu2, beta, eps, liter )
           endif
         endif

         dd = calc_dd_covar( a1, wp2, wp3, wprtp, wpsclrp(i,:),
     .                       rtm, sclrm(i,:), sclrprtp(i,:),
     .                       dt, beta, eps, liter )

c       solve tridiagonal system
         call tridag( aa, bb, cc, dd, sclrprtp(:,i), gr%nnzp, isValid )
         if ( .not. isValid ) then
           print *, "diag_var: sclrprtp tridag failed"
           return
         endif 

       enddo ! i=1, slcrm_dimension

c      clip small and negative values
!        We don't want to do this generally, only to test that the scalars
!        work with their 2nd element is rt
!         ! The surface variances are set elsewhere.  They do not need to be
!         ! clipped at the lowest level.  Brian Griffin.  12/18/05.
!         do k = 2, gr%nnzp, 1
!         if ( wp2(k) >= wtol*wtol ) then
!           sclrprtp(k,2) 
!     .     = max( sclrprtp(k,2), 
!     .            wpsclrp(k,2)*wpsclrp(k,2)/((1.0-Sc(k))*wp2(k)) + eps,
!     .            sclrtol(2)*sclrtol(2) + eps )
!         else
!           sclrprtp(k,2) = max( sclrprtp(k,2), 0.0 )
!         endif
!       enddo ! k
!
      endif ! present (sclr variables)

c Solve for sclrpthlp
      if (     present( sclrm )
     .   .and. present( wpsclrp ) 
     .   .and. present( sclrpthlp ) ) then

c     prepare tridiagonal system
       do i=1, sclrm_dimension

         if ( i == SCLR_THETA ) then
           bb = calc_bb( a1, wp2, wp3, tau, 
     .                   C2thl, dt, nu2, beta, eps, liter )
         else 
           if ( i == SCLR_RT ) then
           bb = calc_bb( a1, wp2, wp3, tau, 
     .                   C2rtthl, dt, nu2, beta, eps, liter )
           endif
         endif

         dd = calc_dd_covar( a1, wp2, wp3, wpthlp, wpsclrp(i,:),
     .                       thlm, sclrm(i,:), sclrpthlp(i,:),
     .                       dt, beta, eps, liter )

c       solve tridiagonal system
         call tridag( aa, bb, cc, dd, sclrpthlp(:,i), gr%nnzp, isValid )
         if ( .not. isValid ) then
           print *, "diag_var: sclrpthlp tridag failed"
           return
         endif 

       enddo ! i=1, slcrm_dimension

c      clip small and negative values
!        We don't want to do this generally, only to test that the scalars
!        work with their 1st element is thl
!         ! The surface variances are set elsewhere.  They do not need to be
!         ! clipped at the lowest level.  Brian Griffin.  12/18/05.
!         do k = 2, gr%nnzp, 1
!         if ( wp2(k) >= wtol*wtol ) then
!           sclrpthlp(k,1) 
!     .     = max( sclrpthlp(k,1), 
!     .            wpsclrp(k,1)*wpsclrp(k,1)/((1.0-Sc(k))*wp2(k)) + eps,
!     .            sclrtol(1)*sclrtol(1) + eps )
!         else
!           sclrpthlp(k,1) = max( sclrpthlp(k,1), 0.0 )
!         endif
!       enddo ! k

      endif ! present (sclr variables)

      return
      end subroutine diag_var

!-----------------------------------------------------------------------
!  FUNCTION calc_aa
!  This only needs to be calculated once per diag_var call
!-----------------------------------------------------------------------
        function calc_aa( a1, wp2, wp3, wmm, nu2, beta, eps )
        use grid_class

        implicit none

!       Input
        real, intent(in), dimension(gr%nnzp) :: a1, wp2, wp3, wmm

        real, intent(in) :: beta, nu2, eps

!       Return type
        real, dimension(gr%nnzp) :: calc_aa

!       Internal
        integer k, km1, kp1

        do k = 2, gr%nnzp-1, 1

           km1   = max(k-1,1)
           kp1   = min(k+1, gr%nnzp)

           calc_aa(k) = 
     .     - wmm(k) * (1.0/2.0) * gr%dzm(k)
     .     - (1.0/3.0) * beta
     .       * (1.0/2.0) * gr%dzm(k)
     .          * ( a1(k) + a1(km1) )
     .          * ( wp3(k) / MAX( wp2(k) + wp2(km1), 2.0*eps ) )
     .     - nu2 * gr%dzm(k) * gr%dzt(k)

        enddo

        ! Boundary Conditions
        ! These are set so that the sfc_var value of rtp2, thlp2, or rtpthlp
        ! can be used at the lowest boundary and the values of those
        ! variables can be set to 0 at the top boundary.  Brian.  12/18/05.
        calc_aa(1) = 0.0
        calc_aa(gr%nnzp) = 0.0

        return
        end function calc_aa

!-----------------------------------------------------------------------
!  FUNCTION calc_bb
!  Needs to be re-calculated when C2 is using an rt or thl value
!-----------------------------------------------------------------------
        function calc_bb( a1, wp2, wp3, tau,
     .                    C2x, dt, nu2, beta, eps, liter )
        use grid_class

        implicit none

!       Input
        real, intent(in), dimension(gr%nnzp) :: a1, wp2, wp3, tau

        real, intent(in)    :: C2x, dt, beta, nu2, eps

        logical, intent(in) :: liter

!       Return type
        real, dimension(gr%nnzp) :: calc_bb

!       Internal
        integer k, km1, kp1

        do k = 2, gr%nnzp-1, 1

           km1 = max(k-1,1)
           kp1 = min(k+1, gr%nnzp)

           calc_bb(k) = 
     .     + (1.0/3.0) * beta
     .       * (1.0/2.0) * gr%dzm(k)
     .          * (   ( a1(kp1) + a1(k) )
     .              * ( wp3(kp1) / MAX( wp2(kp1) + wp2(k), 2.0*eps ) )
     .             -  ( a1(k) + a1(km1) )
     .              * ( wp3(k) / MAX( wp2(k) + wp2(km1), 2.0*eps ) )
     .            )
     .     + C2x / tau(k)
     .     + nu2 * gr%dzm(k) * ( gr%dzt(kp1) + gr%dzt(k) )

           if (liter) calc_bb(k) = calc_bb(k) + 1.0/dt

        enddo

        ! Boundary Conditions
        ! These are set so that the sfc_var value of rtp2, thlp2, or rtpthlp
        ! can be used at the lowest boundary and the values of those
        ! variables can be set to 0 at the top boundary.  Brian.  12/18/05.
        calc_bb(1) = 1.0
        if (liter) calc_bb(1) = calc_bb(1) + 1.0/dt
        calc_bb(gr%nnzp) = 1.0
        if (liter) calc_bb(1) = calc_bb(1) + 1.0/dt

        return
        end function calc_bb

!-----------------------------------------------------------------------
!  FUNCTION calc_cc
!  This only needs to be calculated once per diag_var call
!-----------------------------------------------------------------------
        function calc_cc( a1, wp2, wp3, wmm, nu2, beta, eps )
        use grid_class

        implicit none

!       Input
        real, intent(in), dimension(gr%nnzp) :: a1, wp2, wp3, wmm

        real, intent(in) :: beta, nu2, eps

!       Return type
        real, dimension(gr%nnzp) :: calc_cc

!       Internal
        integer k, km1, kp1

        do k = 2, gr%nnzp-1, 1

           km1 = max(k-1,1)
           kp1 = min(k+1, gr%nnzp)

           calc_cc(k) = 
     .     + wmm(k) * (1.0/2.0) * gr%dzm(k)
     .     + (1.0/3.0) * beta
     .       * (1.0/2.0) * gr%dzm(k)
     .          * ( a1(kp1) + a1(k) )
     .          * ( wp3(kp1) / MAX( wp2(kp1) + wp2(k), 2.0*eps ) )
     .     - nu2 * gr%dzm(k) * gr%dzt(kp1)

        enddo

        ! Boundary Conditions
        ! These are set so that the sfc_var value of rtp2, thlp2, or rtpthlp
        ! can be used at the lowest boundary and the values of those
        ! variables can be set to 0 at the top boundary.  Brian.  12/18/05.
        calc_cc(1) = 0.0
        calc_cc(gr%nnzp) = 0.0

        return
        end function calc_cc

!-----------------------------------------------------------------------
!  FUNCTION calc_dd_covar
!  Calculate dd for the purpose of solving the tridiagonal system 
!    of a covariance term ( e.g rtpthlp )
!-----------------------------------------------------------------------
        function calc_dd_covar( a1, wp2, wp3, wpxap, wpxbp, 
     .                          xam, xbm, xapxbp, 
     .                          dt, beta, eps, liter )
        use grid_class

        implicit none

!       Input
        real, intent(in), dimension(gr%nnzp) :: a1, wp2, wp3
        real, intent(in), dimension(gr%nnzp) :: wpxap, wpxbp
        real, intent(in), dimension(gr%nnzp) :: xam, xbm
        real, intent(in), dimension(gr%nnzp) :: xapxbp

        real, intent(in) :: dt, beta, eps

        logical, intent(in) :: liter

!       Return type
        real, dimension(gr%nnzp) :: calc_dd_covar

!       Internal
        integer k, km1, kp1

        do k = 2, gr%nnzp-1, 1

          km1 = max(k-1,1)
          kp1 = min(k+1, gr%nnzp)

          calc_dd_covar(k) = 
     .    - ( 1.0 - (1.0/3.0) * beta )
     .      * (1.0/4.0) * gr%dzm(k)
     .         * (  ( a1(kp1) + a1(k) )**2.0
     .              * ( wpxap(kp1) + wpxap(k) ) 
     .              * ( wpxbp(kp1) + wpxbp(k) )
     .              * wp3(kp1) 
     .              / ( MAX( wp2(kp1) + wp2(k), 2.0*eps ) )**2.0
     .            - ( a1(k) + a1(km1) )**2.0
     .              * ( wpxap(k) + wpxap(km1) ) 
     .              * ( wpxbp(k) + wpxbp(km1) )
     .              * wp3(k) 
     .              / ( MAX( wp2(k) + wp2(km1), 2.0*eps ) )**2.0
     .           )
     .    - wpxbp(k) * gr%dzm(k) * ( xam(kp1) - xam(k) )
     .    - wpxap(k) * gr%dzm(k) * ( xbm(kp1) - xbm(k) )

          if (liter) calc_dd_covar(k) = calc_dd_covar(k) 
     .                                + 1.0/dt*xapxbp(k)

        enddo

        ! Boundary Conditions
        ! These are set so that the sfc_var value of rtp2, thlp2, or rtpthlp
        ! can be used at the lowest boundary and the values of those
        ! variables can be set to 0 at the top boundary.  Brian.  12/18/05.
        calc_dd_covar(1) = xapxbp(1)
        if (liter) calc_dd_covar(1) = calc_dd_covar(1) 
     .                              + 1.0/dt*xapxbp(1)
        calc_dd_covar(gr%nnzp) = 0.0
        if (liter) calc_dd_covar(gr%nnzp) = calc_dd_covar(gr%nnzp) 
     .                                     + 1.0/dt*xapxbp(gr%nnzp)

        return
        end function calc_dd_covar

!-----------------------------------------------------------------------
!  FUNCTION calc_dd_var
!  Calculate dd for the purpose of solving the tridiagonal system 
!    of a variance term (e.g. rtp2)
!-----------------------------------------------------------------------
        function calc_dd_var( a1, wp2, wp3, wpxp, xm,
     .                        xp2, dt, beta, eps, liter )
        use grid_class

        implicit none

!       Input
        real, intent(in), dimension(gr%nnzp) :: a1, wp2, wp3
        real, intent(in), dimension(gr%nnzp) :: wpxp, xm, xp2

        real, intent(in) :: dt, beta, eps

        logical, intent(in) :: liter

!       Return type
        real, dimension(gr%nnzp) :: calc_dd_var

!       Internal
        integer k, km1, kp1

        do k = 2, gr%nnzp-1, 1

           km1 = max(k-1,1)
           kp1 = min(k+1, gr%nnzp)

           calc_dd_var(k) =
     .     - ( 1.0 - (1.0/3.0) * beta )
     .       * (1.0/4.0) * gr%dzm(k)
     .          * (   ( a1(kp1) + a1(k) )**2.0
     .                * wp3(kp1) * ( wpxp(kp1) + wpxp(k) )**2.0
     .                / ( MAX( wp2(kp1) + wp2(k), 2.0*eps ) )**2.0
     .             -  ( a1(k) + a1(km1) )**2.0
     .                * wp3(k) * ( wpxp(k) + wpxp(km1) )**2.0
     .                / ( MAX( wp2(k) + wp2(km1), 2.0*eps ) )**2.0
     .            )
     .     - 2.0 * wpxp(k) * gr%dzm(k) * ( xm(kp1) - xm(k) )

           if (liter) calc_dd_var(k) = calc_dd_var(k) + 1.0/dt*xp2(k)

        enddo

        ! Boundary Conditions
        ! These are set so that the sfc_var value of rtp2, thlp2, or rtpthlp
        ! can be used at the lowest boundary and the values of those
        ! variables can be set to 0 at the top boundary.  Brian.  12/18/05.
        calc_dd_var(1) = xp2(1)
        if (liter) calc_dd_var(1) = calc_dd_var(1) + 1.0/dt*xp2(1)
        calc_dd_var(gr%nnzp) = 0.0
        if (liter) calc_dd_var(gr%nnzp) = calc_dd_var(gr%nnzp) 
     .                                   + 1.0/dt*xp2(gr%nnzp)

        return
        end function calc_dd_var

!-----------------------------------------------------------------------
      end module diagnose_variances
