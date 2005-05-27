!-----------------------------------------------------------------------
! subroutine diag_var( )

! Subprogram to diagnose variances by solving steady-state equations
!-----------------------------------------------------------------------
        subroutine diag_var( tau, wmm, qtm, wpqtp, thlm, wpthlp, 
     .                       wp2, wp3, Sc,
     .                       qtp2, thlp2, qtpthlp,
     .                       liter, dt, isValid )

       use constants
       use grid_class
       
       implicit none

c Input variables

       real, intent(in) :: tau(1:gr%nnzp)
       real, intent(in) :: wmm(1:gr%nnzp)
       real, intent(in) :: qtm(1:gr%nnzp)
       real, intent(in) :: wpqtp(1:gr%nnzp)
       real, intent(in) :: thlm(1:gr%nnzp)
       real, intent(in) :: wpthlp(1:gr%nnzp)
       real, intent(in) :: wp2(1:gr%nnzp)
       real, intent(in) :: wp3(1:gr%nnzp)
       real, intent(in) :: Sc(1:gr%nnzp)

       logical, intent(in) :: liter
       real, intent(in)    :: dt

c Output variables

       real, intent(out) :: qtp2(1:gr%nnzp)
       real, intent(out) :: thlp2(1:gr%nnzp)
       real, intent(out) :: qtpthlp(1:gr%nnzp)

c Local constants

!       real, parameter :: eps = 1.e-10  ! now in constants.F
c Misc. Variables -D.S.
       logical, intent(inout) :: isValid

c Local arrays allocated dynamically on the stack

       real a1(1:gr%nnzp), a2(1:gr%nnzp)
       real aa(1:gr%nnzp), bb(1:gr%nnzp), cc(1:gr%nnzp), dd(1:gr%nnzp)

c Loop indices

       integer k
       integer km1,kp1

c Define a1 and a2

       a1 = 1.0 / (1.0-Sc)
       a2 = 1.0 / ( (1.0-Sc)*(1.0-Sc) )

c Solve for qtp2

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

         dd(k) = - 2.0 * wpqtp(k) * gr%dzm(k) * (qtm(kp1)-qtm(k))
     .           - 0.5 * a2(k)
     .                 * (2./3.) * ( 2 + (1-beta) )
     .                 * gr%dzm(k)
     .                 * (  wp3(kp1)
     .                      * (wpqtp(k)+wpqtp(kp1))**2
     .                      / max(wp2(k)+wp2(kp1),2.0*eps)**2
     .                    - wp3(k)
     .                      * (wpqtp(km1)+wpqtp(k))**2
     .                      / max(wp2(km1)+wp2(k),2.0*eps)**2
     .                   )
         if (liter) dd(k) = dd(k) + 1.0/dt*qtp2(k)

       enddo
       
c      solve tridiagonal system

       call tridag( aa, bb, cc, dd, qtp2, gr%nnzp, isValid )

       if ( .not. isValid ) then
          print *, "diag_var: tridag failed"
          return
       endif

c      clip small and negative values

       do k=1,gr%nnzp

         if ( wp2(k) >= wtol*wtol ) then
           qtp2(k) 
     .     = max( qtp2(k), 
     .            wpqtp(k)*wpqtp(k)/((1.0-Sc(k))*wp2(k)) + eps,
     .            qttol*qttol + eps )
         else
           qtp2(k) = max( qtp2(k), 0.0 )
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
c       print *, "diag_var line 147"
       call tridag( aa, bb, cc, dd, thlp2, gr%nnzp, isValid )

       if ( .not. isValid ) then
          print *, "diag_var: tridag failed"
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

c Solve for qtpthlp

c      prepare tridiagonal system

       do k=1,gr%nnzp

         km1 = max(k-1,1)
         kp1 = min(k+1,gr%nnzp)

         dd(k) = - wpthlp(k) * gr%dzm(k) * (qtm(kp1)-qtm(k))
     .           - wpqtp(k)  * gr%dzm(k) * (thlm(kp1)-thlm(k))
     .           - 0.5 * a2(k)
     .                 * (2./3.) * ( 2 + (1-beta) )
     .                 * gr%dzm(k)
     .                 * (  wp3(kp1)
     .                      * (wpqtp(k)+wpqtp(kp1))
     .                      * (wpthlp(k)+wpthlp(kp1))
     .                      / max(wp2(k)+wp2(kp1),2.0*eps)**2
     .                    - wp3(k)
     .                      * (wpqtp(km1)+wpqtp(k))
     .                      * (wpthlp(km1)+wpthlp(k))
     .                      / max(wp2(km1)+wp2(k),2.0*eps)**2
     .                   )
         if (liter) dd(k) = dd(k) + 1.0/dt*qtpthlp(k)

       enddo
       
c      solve tridiagonal system
       call tridag(aa,bb,cc,dd,qtpthlp,gr%nnzp, isValid)
       if ( .not. isValid ) then
          print *, "diag_var: tridag failed"
          return
       endif

       return
       end subroutine diag_var
!-----------------------------------------------------------------------
