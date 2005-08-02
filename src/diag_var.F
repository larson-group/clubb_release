! $Id: diag_var.F,v 1.10 2005-08-02 17:32:39 dschanen Exp $
#define SCLR_THETA 1
#define SCLR_RT 2
!-----------------------------------------------------------------------
!  This is a module because optional variables require it to be.
!-----------------------------------------------------------------------
      module diagnose_variances

        implicit none

        public diag_var

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
!  note: if liter=.true., then these would need to have an attrib. of type (inout)

       real, intent(out) :: rtp2(1:gr%nnzp)
       real, intent(out) :: thlp2(1:gr%nnzp)
       real, intent(out) :: rtpthlp(1:gr%nnzp)

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

       real a1(1:gr%nnzp), a2(1:gr%nnzp)
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

c Define a1 and a2

       a1 = 1.0 / (1.0-Sc)
       a2 = 1.0 / ( (1.0-Sc)*(1.0-Sc) )

c Solve for rtp2

c      prepare tridiagonal system

       do k=1,gr%nnzp

         km1   = max(k-1,1)
         kp1   = min(k+1,gr%nnzp)

         aa(k) = - 0.5 * a1(k)
     .                 * (2./3.) * beta
     .                 * gr%dzm(k)
     .                 * wp3(k)/max(wp2(km1)+wp2(k),2.0*eps)
     .           + 0.5 * wmm(k) * gr%dzm(k)
     .           - nu2 * gr%dzm(k) * gr%dzt(k)

         bb(k) = C2rt / tau(k)
     .           + 0.5 * a1(k)
     .                 * (2./3.) * beta
     .                 * gr%dzm(k)
     .                 * ( wp3(kp1)/max(wp2(k)+wp2(kp1),2.0*eps)
     .                     - wp3(k)/max(wp2(km1)+wp2(k),2.0*eps) )
     .           + nu2 * gr%dzm(k)*( gr%dzt(kp1) + gr%dzt(k) )
         if (liter) bb(k) = bb(k) + 1.0/dt

         cc(k) =   0.5 * a1(k)
     .                 * (2./3.) * beta
     .                 * gr%dzm(k)
     .                 * wp3(kp1)/max(wp2(k)+wp2(kp1),2.0*eps)
     .           - 0.5 * wmm(k) * gr%dzm(k)
     .           - nu2 * gr%dzm(k) * gr%dzt(kp1)

         dd(k) = - 2.0 * wprtp(k) * gr%dzm(k) * (rtm(kp1)-rtm(k))
     .           - 0.5 * a2(k)
     .                 * (2./3.) * ( 2 + (1-beta) )
     .                 * gr%dzm(k)
     .                 * (  wp3(kp1)
     .                      * (wprtp(k)+wprtp(kp1))**2
     .                      / max(wp2(k)+wp2(kp1),2.0*eps)**2
     .                    - wp3(k)
     .                      * (wprtp(km1)+wprtp(k))**2
     .                      / max(wp2(km1)+wp2(k),2.0*eps)**2
     .                   )
         if (liter) dd(k) = dd(k) + 1.0/dt*rtp2(k)

#ifdef STATS
         if (lstats_samp) then

           if ( irtp2_dp1 * ithlp2_dp1 * irtpthlp_dp1 > 0 ) then
             zmscr01(k) = - C2rt / tau(k)
           endif
         
           if ( irtp2_dp2 * ithlp2_dp2 * irtpthlp_dp2 > 0 ) then
             zmscr02(k)
     .       = + nu2 * gr%dzm(k) * gr%dzt(k)
             zmscr03(k)
     .       = - nu2 * gr%dzm(k)*( gr%dzt(kp1) + gr%dzt(k) )
             zmscr04(k)
     .       = + nu2 * gr%dzm(k) * gr%dzt(kp1)
           endif

           if ( irtp2_ta * ithlp2_ta * irtpthlp_ta > 0 ) then
             zmscr05(k)
     .       = + 0.5 * a1(k)
     .               * (2./3.) * beta
     .               * gr%dzm(k)
     .               * wp3(k)/max(wp2(km1)+wp2(k),2.0*eps)
             zmscr06(k)
     .       = - 0.5 * a1(k)
     .               * (2./3.) * beta
     .               * gr%dzm(k)
     .               * ( wp3(kp1)/max(wp2(k)+wp2(kp1),2.0*eps)
     .                   - wp3(k)/max(wp2(km1)+wp2(k),2.0*eps) )
             zmscr07(k)
     .       = - 0.5 * a1(k)
     .               * (2./3.) * beta
     .               * gr%dzm(k)
     .               * wp3(kp1)/max(wp2(k)+wp2(kp1),2.0*eps)
           endif

           if ( irtp2_ma * ithlp2_ma * irtpthlp_ma > 0 ) then
             zmscr08(k) = - 0.5 * wmm(k) * gr%dzm(k)
             zmscr09(k) = + 0.5 * wmm(k) * gr%dzm(k)
           endif

           if ( irtp2_ta > 0 ) then
             zm%x(k,irtp2_ta) 
     .       = zm%x(k,irtp2_ta)
     .         - 0.5 * a2(k)
     .               * (2./3.) * ( 2 + (1-beta) )
     .               * gr%dzm(k)
     .               * (  wp3(kp1)
     .                    * (wprtp(k)+wprtp(kp1))**2
     .                    / max(wp2(k)+wp2(kp1),2.0*eps)**2
     .                  - wp3(k)
     .                    * (wprtp(km1)+wprtp(k))**2
     .                    / max(wp2(km1)+wp2(k),2.0*eps)**2
     .                 )
           endif

         endif
#endif

       enddo
       
c      solve tridiagonal system
       call tridag( aa, bb, cc, dd, rtp2, gr%nnzp, isValid )

       if ( .not. isValid ) then
          print *, "diag_var: rtp2 tridag failed"
          return
       endif
#ifdef STATS
       if (lstats_samp) then

          do k=1,gr%nnzp

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
     .        - 2.0 * wprtp(k) * gr%dzm(k) * (rtm(kp1)-rtm(k))
             zm%n(k,irtp2_tp) = zm%n(k,irtp2_tp) + 1
           endif

          enddo

       endif
#endif

#ifdef STATS
       if (lstats_samp .and. irtp2_cl > 0 ) then
         zm%x(:,irtp2_cl) = zm%x(:,irtp2_cl) - rtp2 / dt
       endif
#endif

c      clip small and negative values

       do k=1,gr%nnzp

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

       do k=1,gr%nnzp

         km1 = max(k-1,1)
         kp1 = min(k+1,gr%nnzp)

         bb(k) = C2thl / tau(k)
     .           + 0.5 * a1(k)
     .                 * (2./3.) * beta
     .                 * gr%dzm(k)
     .                 * ( wp3(kp1)/max(wp2(k)+wp2(kp1),2.0*eps)
     .                     - wp3(k)/max(wp2(km1)+wp2(k),2.0*eps) )
     .           + nu2 * gr%dzm(k)*( gr%dzt(kp1) + gr%dzt(k) )
         if (liter) bb(k) = bb(k) + 1.0/dt

         dd(k) = - 2.0 * wpthlp(k) * gr%dzm(k) * (thlm(kp1)-thlm(k))
     .           - 0.5 * a2(k)
     .                 * (2./3.) * ( 2 + (1-beta) )
     .                 * gr%dzm(k)
     .                 * (  wp3(kp1)
     .                      * (wpthlp(k)+wpthlp(kp1))**2
     .                      / max(wp2(k)+wp2(kp1),2.0*eps)**2
     .                    - wp3(k)
     .                      * (wpthlp(km1)+wpthlp(k))**2
     .                      / max(wp2(km1)+wp2(k),2.0*eps)**2
     .                   )
         if (liter) dd(k) = dd(k) + 1.0/dt*thlp2(k)

#ifdef STATS
         if (lstats_samp) then

           if ( ithlp2_ta > 0 ) then
             zm%x(k,ithlp2_ta) 
     .       = zm%x(k,ithlp2_ta)
     .           - 0.5 * a2(k)
     .                 * (2./3.) * ( 2 + (1-beta) )
     .                 * gr%dzm(k)
     .                 * (  wp3(kp1)
     .                      * (wpthlp(k)+wpthlp(kp1))**2
     .                      / max(wp2(k)+wp2(kp1),2.0*eps)**2
     .                    - wp3(k)
     .                      * (wpthlp(km1)+wpthlp(k))**2
     .                      / max(wp2(km1)+wp2(k),2.0*eps)**2
     .                   )
           endif

         endif
#endif

       enddo
       
c      solve tridiagonal system
       call tridag( aa, bb, cc, dd, thlp2, gr%nnzp, isValid )

       if ( .not. isValid ) then
          print *, "diag_var: thlp2 tridag failed"
          return
       endif

#ifdef STATS
         if (lstats_samp) then

          do k=1,gr%nnzp

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
     .        - 2.0 * wpthlp(k) * gr%dzm(k) * (thlm(kp1)-thlm(k))
             zm%n(k,ithlp2_tp) = zm%n(k,ithlp2_tp) + 1

           endif

          enddo

         endif
#endif

#ifdef STATS
       if (lstats_samp .and. ithlp2_cl > 0 ) then
         zm%x(:,ithlp2_cl) = zm%x(:,ithlp2_cl) - thlp2 / dt
       endif
#endif


c      clip small and negative values

       do k=1,gr%nnzp

         if ( wp2(k) >= wtol*wtol ) then
           thlp2(k) 
     .     = max( thlp2(k), 
     .            wpthlp(k)*wpthlp(k)/((1.0-Sc(k))*wp2(k)) + eps,
     .            thltol*thltol + eps )
         else
           thlp2(k) = max( thlp2(k), 0.0 )
         endif

       enddo

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
       do k=1, gr%nnzp
         km1 = max(k-1,1)
         kp1 = min(k+1,gr%nnzp)

         if ( i == SCLR_THETA ) then
           bb(k) = C2thl / tau(k)
     .             + 0.5 * a1(k)
     .                   * (2./3.) * beta
     .                   * gr%dzm(k)
     .                   * ( wp3(kp1)/max(wp2(k)+wp2(kp1),2.0*eps)
     .                       - wp3(k)/max(wp2(km1)+wp2(k),2.0*eps) )
     .             + nu2 * gr%dzm(k)*( gr%dzt(kp1) + gr%dzt(k) )
           if (liter) bb(k) = bb(k) + 1.0/dt
         else 
           if ( i == SCLR_RT ) then
             bb(k) = C2rt / tau(k)
     .               + 0.5 * a1(k)
     .                     * (2./3.) * beta
     .                     * gr%dzm(k)
     .                     * ( wp3(kp1)/max(wp2(k)+wp2(kp1),2.0*eps)
     .                         - wp3(k)/max(wp2(km1)+wp2(k),2.0*eps) )
     .               + nu2 * gr%dzm(k)*( gr%dzt(kp1) + gr%dzt(k) )
             if (liter) bb(k) = bb(k) + 1.0/dt
           endif
         endif

         dd(k) = - 2.0 * wpsclrp(k,i) * gr%dzm(k)
     .                 * ( sclrm(kp1,i) - sclrm(k,i) )
     .           - 0.5 * a2(k)
     .                 * (2./3.) * ( 2 + (1-beta) )
     .                 * gr%dzm(k)
     .                 * (  wp3(kp1)
     .                      * (wpsclrp(k,i)+wpsclrp(kp1,i))**2
     .                      / max(wp2(k)+wp2(kp1),2.0*eps)**2
     .                    - wp3(k)
     .                      * (wpsclrp(km1,i)+wpsclrp(k,i))**2
     .                      / max(wp2(km1)+wp2(k),2.0*eps)**2
     .                   )
         if (liter) dd(k) = dd(k) + 1.0/dt*sclrp2(k,i)

       enddo ! k=1...nz
       
c      solve tridiagonal system
       call tridag( aa, bb, cc, dd, sclrp2(:,i), gr%nnzp, isValid )

       if ( .not. isValid ) then
          print *, "diag_var: sclrp2 tridag failed"
          return
       endif
c      clip small and negative values

       do k=1,gr%nnzp

         if ( wp2(k) >= wtol*wtol ) then
           sclrp2(k,i) 
     .     = max( sclrp2(k,i), 
     .            wpsclrp(k,i)*wpsclrp(k,i)/((1.0-Sc(k))*wp2(k)) + eps,
     .            sclrtol(i)*sclrtol(i) + eps )
         else
           sclrp2(k,i) = max( sclrp2(k,i), 0.0 )
         endif

       enddo ! k
       enddo ! i
      endif
! %% Debug
!      sclrp2(:,2) = rtp2
! %% Debug

c Solve for rtpthlp

c      prepare tridiagonal system

       do k=1,gr%nnzp

         km1 = max(k-1,1)
         kp1 = min(k+1,gr%nnzp)

         bb(k) = C2rtthl / tau(k)
     .           + 0.5 * a1(k)
     .                 * (2./3.) * beta
     .                 * gr%dzm(k)
     .                 * ( wp3(kp1)/max(wp2(k)+wp2(kp1),2.0*eps)
     .                     - wp3(k)/max(wp2(km1)+wp2(k),2.0*eps) )
     .           + nu2 * gr%dzm(k)*( gr%dzt(kp1) + gr%dzt(k) )
         if (liter) bb(k) = bb(k) + 1.0/dt

         dd(k) = - wpthlp(k) * gr%dzm(k) * (rtm(kp1)-rtm(k))
     .           - wprtp(k)  * gr%dzm(k) * (thlm(kp1)-thlm(k))
     .           - 0.5 * a2(k)
     .                 * (2./3.) * ( 2 + (1-beta) )
     .                 * gr%dzm(k)
     .                 * (  wp3(kp1)
     .                      * (wprtp(k)+wprtp(kp1))
     .                      * (wpthlp(k)+wpthlp(kp1))
     .                      / max(wp2(k)+wp2(kp1),2.0*eps)**2
     .                    - wp3(k)
     .                      * (wprtp(km1)+wprtp(k))
     .                      * (wpthlp(km1)+wpthlp(k))
     .                      / max(wp2(km1)+wp2(k),2.0*eps)**2
     .                   )
         if (liter) dd(k) = dd(k) + 1.0/dt*rtpthlp(k)

#ifdef STATS
         if (lstats_samp) then

           if ( irtpthlp_ta > 0 ) then
             zm%x(k,irtpthlp_ta) 
     .       = zm%x(k,irtpthlp_ta)
     .           - 0.5 * a2(k)
     .                 * (2./3.) * ( 2 + (1-beta) )
     .                 * gr%dzm(k)
     .                 * (  wp3(kp1)
     .                      * (wprtp(k)+wprtp(kp1))
     .                      * (wpthlp(k)+wpthlp(kp1))
     .                      / max(wp2(k)+wp2(kp1),2.0*eps)**2
     .                    - wp3(k)
     .                      * (wprtp(km1)+wprtp(k))
     .                      * (wpthlp(km1)+wpthlp(k))
     .                      / max(wp2(km1)+wp2(k),2.0*eps)**2
     .                   )
           endif

         endif
#endif

       enddo
       
c      solve tridiagonal system
       call tridag( aa, bb, cc, dd, rtpthlp, gr%nnzp, isValid )
       if ( .not. isValid ) then
          print *, "diag_var: rtpthlp tridag failed"
          return
       endif

#ifdef STATS
         if (lstats_samp) then

          do k=1,gr%nnzp

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
     .       - wprtp(k)  * gr%dzm(k) * (thlm(kp1)-thlm(k))
             zm%n(k,irtpthlp_tp1) = zm%n(k,irtpthlp_tp1) + 1
           endif

           if ( irtpthlp_tp2 > 0 ) then
             zm%x(k,irtpthlp_tp2) = zm%x(k,irtpthlp_tp2)
     .       - wpthlp(k) * gr%dzm(k) * (rtm(kp1)-rtm(k))
             zm%n(k,irtpthlp_tp2) = zm%n(k,irtpthlp_tp2) + 1
           endif

          enddo

         endif
#endif

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
        do k=1, gr%nnzp
          km1 = max(k-1,1)
          kp1 = min(k+1,gr%nnzp)

         if ( i == SCLR_THETA ) then
           bb(k) = C2rtthl / tau(k)
     .             + 0.5 * a1(k)
     .                   * (2./3.) * beta
     .                   * gr%dzm(k)
     .                   * ( wp3(kp1)/max(wp2(k)+wp2(kp1),2.0*eps)
     .                       - wp3(k)/max(wp2(km1)+wp2(k),2.0*eps) )
     .             + nu2 * gr%dzm(k)*( gr%dzt(kp1) + gr%dzt(k) )
           if (liter) bb(k) = bb(k) + 1.0/dt
         else 
           if ( i == SCLR_RT ) then
             bb(k) = C2rt / tau(k)
     .               + 0.5 * a1(k)
     .                     * (2./3.) * beta
     .                     * gr%dzm(k)
     .                     * ( wp3(kp1)/max(wp2(k)+wp2(kp1),2.0*eps)
     .                         - wp3(k)/max(wp2(km1)+wp2(k),2.0*eps) )
     .               + nu2 * gr%dzm(k)*( gr%dzt(kp1) + gr%dzt(k) )
             if (liter) bb(k) = bb(k) + 1.0/dt
           endif
         endif

          dd(k) = - wprtp(k) * gr%dzm(k) * (sclrm(kp1,i)-sclrm(k,i)) 
     .            - wpsclrp(k,i) * gr%dzm(k) * (rtm(kp1)-rtm(k))
     .            - 0.5 * a2(k)
     .                  * (2./3.) * ( 2 + (1-beta) )
     .                  * gr%dzm(k)
     .                  * (  wp3(kp1)
     .                       * (wpsclrp(k,i)+wpsclrp(kp1,i))
     .                       * (wprtp(k)+wprtp(kp1))
     .                       / max(wp2(k)+wp2(kp1),2.0*eps)**2
     .                     - wp3(k)
     .                       * (wpsclrp(km1,i)+wpsclrp(k,i))
     .                       * (wprtp(km1)+wprtp(k))
     .                       / max(wp2(km1)+wp2(k),2.0*eps)**2
     .                    )
          if (liter) dd(k) = dd(k) + 1.0/dt*sclrprtp(k,i)

        enddo ! k=1, nz
       
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
!       do k=1, gr%nnzp
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
        do k=1, gr%nnzp
          km1 = max(k-1,1)
          kp1 = min(k+1,gr%nnzp)

         if ( i == SCLR_THETA ) then
           bb(k) = C2thl / tau(k)
     .             + 0.5 * a1(k)
     .                   * (2./3.) * beta
     .                   * gr%dzm(k)
     .                   * ( wp3(kp1)/max(wp2(k)+wp2(kp1),2.0*eps)
     .                       - wp3(k)/max(wp2(km1)+wp2(k),2.0*eps) )
     .             + nu2 * gr%dzm(k)*( gr%dzt(kp1) + gr%dzt(k) )
           if (liter) bb(k) = bb(k) + 1.0/dt
         else 
           if ( i == SCLR_RT ) then
             bb(k) = C2rtthl / tau(k)
     .               + 0.5 * a1(k)
     .                     * (2./3.) * beta
     .                     * gr%dzm(k)
     .                     * ( wp3(kp1)/max(wp2(k)+wp2(kp1),2.0*eps)
     .                         - wp3(k)/max(wp2(km1)+wp2(k),2.0*eps) )
     .               + nu2 * gr%dzm(k)*( gr%dzt(kp1) + gr%dzt(k) )
             if (liter) bb(k) = bb(k) + 1.0/dt
           endif
         endif

          dd(k) = - wpthlp(k) * gr%dzm(k) * (sclrm(kp1,i)-sclrm(k,i)) 
     .            - wpsclrp(k,i) * gr%dzm(k) * (thlm(kp1)-thlm(k))
     .            - 0.5 * a2(k)
     .                  * (2./3.) * ( 2 + (1-beta) )
     .                  * gr%dzm(k)
     .                  * (  wp3(kp1)
     .                       * (wpsclrp(k,i)+wpsclrp(kp1,i))
     .                       * (wpthlp(k)+wpthlp(kp1))
     .                       / max(wp2(k)+wp2(kp1),2.0*eps)**2
     .                     - wp3(k)
     .                       * (wpsclrp(km1,i)+wpsclrp(k,i))
     .                       * (wpthlp(km1)+wpthlp(k))
     .                       / max(wp2(km1)+wp2(k),2.0*eps)**2
     .                    )
          if (liter) dd(k) = dd(k) + 1.0/dt*sclrpthlp(k,i)

        enddo ! k=1, nz
       
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
!       do k=1, gr%nnzp
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
      end module diagnose_variances
