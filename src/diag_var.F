      module diagnose_variances
        implicit none
        contains
!-----------------------------------------------------------------------
! subroutine diag_var( )
! Subprogram to diagnose variances by solving steady-state equations
!-----------------------------------------------------------------------
        subroutine diag_var( tau, wmm, rtm, wprtp, thlm, wpthlp, 
     .                       wp2, wp3, Sc,
     .                       rtp2, thlp2, rtpthlp,
     .                       liter, dt, isValid,
     .                       sclrm, wpsclrp,              ! optional
     .                       sclrp2, sclrpthlp )          ! optional

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
     . dimension(sclrm_dimension, gr%nnzp) :: sclrm, wpsclrp

! Output mixing scheme variables (optional)

       real, optional, intent(out), 
     . dimension(sclrm_dimension, gr%nnzp) :: sclrp2
       real, optional, intent(out), 
     . dimension(sclrm_dimension, gr%nnzp) :: sclrpthlp

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

         bb(k) = C2 / tau(k)
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
c      note: we can reuse the same aa, bb, cc arrays

       do k=1,gr%nnzp

         km1 = max(k-1,1)
         kp1 = min(k+1,gr%nnzp)

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
     .            wpthlp(k)*wpthlp(k)/((1.0-Sc(k))*wp2(k)) + 1.0e-8,
     .            thltol*thltol + 1.0e-8 )
         else
           thlp2(k) = max( thlp2(k), 0.0 )
         endif

       enddo

c Solve for sclrp2(1), the variance of a passive scalar

c      prepare tridiagonal system
c      note: we can reuse the same aa, bb, cc arrays

      if (     present( sclrm )
     .   .and. present( wpsclrp ) 
     .   .and. present( sclrp2 ) ) then

       do k=1,gr%nnzp

         km1 = max(k-1,1)
         kp1 = min(k+1,gr%nnzp)

         dd(k) = - 2.0 * wpsclrp(1,k) * gr%dzm(k) 
     .                 * ( sclrm(1,kp1)-sclrm(1,k) )
     .           - 0.5 * a2(k)
     .                 * (2./3.) * ( 2 + (1-beta) )
     .                 * gr%dzm(k)
     .                 * (  wp3(kp1)
     .                      * (wpsclrp(1,k)+wpsclrp(1,kp1))**2
     .                      / max( wp2(k)+wp2(kp1), 2.0*eps )**2
     .                    - wp3(k)
     .                      * (wpsclrp(1,km1)+wpsclrp(1,k))**2
     .                      / max( wp2(km1)+wp2(k), 2.0*eps )**2
     .                   )
         if (liter) dd(k) = dd(k) + 1.0/dt*sclrp2(1,k)

       enddo
       
c      solve tridiagonal system
       call tridag( aa, bb, cc, dd, sclrp2(1,:), gr%nnzp, isValid )

       if ( .not. isValid ) then
          print *, "diag_var:  sclrp2 tridag failed"
          return
       endif

c      clip small and negative values

       do k=1, gr%nnzp

         if ( wp2(k) >= wtol*wtol ) then
           sclrp2(1,k) = max( sclrp2(1,k), 
     .                        wpsclrp(1,k)*wpsclrp(1,k)/
     .                        ( ( 1.0-Sc(k) )*wp2(k) ) + eps,
     .                        sclrtol*sclrtol + eps )
         else
           sclrp2(1,k) = max( sclrp2(1,k), 0.0 )
         endif
       enddo ! k=1, nz
      endif ! if present( sclrp2 )

c Solve for rtpthlp

c      prepare tridiagonal system

       do k=1,gr%nnzp

         km1 = max(k-1,1)
         kp1 = min(k+1,gr%nnzp)

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

c Solve for sclrpthlp
      if (     present( sclrm )
     .   .and. present( wpsclrp ) 
     .   .and. present( sclrpthlp ) ) then

c     prepare tridiagonal system
      do i=1, sclrm_dimension
        do k=1,gr%nnzp
          km1 = max(k-1,1)
          kp1 = min(k+1,gr%nnzp)

          dd(k) = - wpthlp(k) * gr%dzm(k) * (sclrm(i,kp1)-sclrm(i,k)) 
     .            - wpsclrp(i,k) * gr%dzm(k) * (thlm(kp1)-thlm(k))
     .            - 0.5 * a2(k)
     .                  * (2./3.) * ( 2 + (1-beta) )
     .                  * gr%dzm(k)
     .                  * (  wp3(kp1)
     .                       * (wpsclrp(i,k)+wpsclrp(i,kp1))
     .                       * (wpthlp(k)+wpthlp(kp1))
     .                       / max(wp2(k)+wp2(kp1),2.0*eps)**2
     .                     - wp3(k)
     .                       * (wpsclrp(i,km1)+wpsclrp(i,k))
     .                       * (wpthlp(km1)+wpthlp(k))
     .                       / max(wp2(km1)+wp2(k),2.0*eps)**2
     .                    )
          if (liter) dd(k) = dd(k) + 1.0/dt*sclrpthlp(i,k)

        enddo ! k=1, nz
       
c       solve tridiagonal system
        call tridag( aa, bb, cc, dd, sclrpthlp(i,:), gr%nnzp, isValid )
        if ( .not. isValid ) then
          print *, "diag_var: sclrpthlp tridag failed"
          return
        endif 
      enddo ! i=1, slcrm_dimension

      endif ! present (sclr variables)

      return
      end subroutine diag_var
!-----------------------------------------------------------------------
      end module diagnose_variances
