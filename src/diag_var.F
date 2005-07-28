! $Id: diag_var.F,v 1.9 2005-07-28 20:54:10 dschanen Exp $
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

       enddo
       
c      solve tridiagonal system
       call tridag( aa, bb, cc, dd, rtp2, gr%nnzp, isValid )

       if ( .not. isValid ) then
          print *, "diag_var: rtp2 tridag failed"
          return
       endif

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

       enddo
       
c      solve tridiagonal system
       call tridag( aa, bb, cc, dd, thlp2, gr%nnzp, isValid )

       if ( .not. isValid ) then
          print *, "diag_var: thlp2 tridag failed"
          return
       endif
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

       enddo
       
c      solve tridiagonal system
       call tridag( aa, bb, cc, dd, rtpthlp, gr%nnzp, isValid )
       if ( .not. isValid ) then
          print *, "diag_var: rtpthlp tridag failed"
          return
       endif

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
