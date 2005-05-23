      subroutine tau(thlm, thvm, rtm, rlm, cf, p, exner, taut)

      use constants
      use grid_class

      implicit none

! External

      real, external :: rsat

! Input

      real, dimension(1:gr%nnzp), intent(in) :: thlm
      real, dimension(1:gr%nnzp), intent(in) :: thvm
      real, dimension(1:gr%nnzp), intent(in) :: rtm
      real, dimension(1:gr%nnzp), intent(in) :: rlm
      real, dimension(1:gr%nnzp), intent(in) :: cf
      real, dimension(1:gr%nnzp), intent(in) :: p
      real, dimension(1:gr%nnzp), intent(in) :: exner

! Output

      real, dimension(1:gr%nnzp), intent(out) :: taut

! Internal

      integer k
      real n2, n2dry, n2moist
      real thm(1:gr%nnzp), tm(1:gr%nnzp), rs(1:gr%nnzp)
      real dthm(1:gr%nnzp), dthvm(1:gr%nnzp),
     .     drs(1:gr%nnzp), drt(1:gr%nnzp)

! Compute dry potential temparature

      thm = thlm + Lv/(cp*exner) * rlm

! Compute temperature

      tm = thm*exner

! Compute saturation ratio

      do k=1,gr%nnzp
       rs(k) = rsat(p(k),tm(k))
      enddo

! Compute vertical derivative

      dthm = ddzm( zt2zm( thm ) )
      dthm(1) = dthm(2)
      dthm(gr%nnzp) = dthm(gr%nnzp-1)

      dthvm = ddzm( zt2zm( thvm ) )
      dthvm(1) = dthvm(2)
      dthvm(gr%nnzp) = dthvm(gr%nnzp-1)

      drs = ddzm( zt2zm( rs ) )
      drs(1) = drs(2)
      drs(gr%nnzp) = drs(gr%nnzp-1)

      drt = ddzm( zt2zm( rtm ) )
      drt(1) = drt(2)
      drt(gr%nnzp) = drt(gr%nnzp-1)

! Compute tau

      do k=2,gr%nnzp-1


        n2dry = g/T0 * dthvm(k)


        n2moist = g*(  ( 1.0 + (Lv*rs(k))/(R*tm(k)) )
     .                 /( 1.0 + (ep1*Lv*Lv*rs(k))/(cp*R*tm(k)*tm(k)) )
     .                 *( dthm(k)/thm(k) + Lv/(cp*tm(k))*drs(k) )
     .                - drt(k)
     .               )

        if ( rlm(k) < 1.e-8 ) then
          n2 = n2dry
        else
          n2 = n2moist
        endif
!        n2 = cf(k)*n2moist + (1.0-cf(k))*n2dry

        n2 = max( n2, 1.e-8 )
        n2moist = max( n2moist, 1.e-8 )
        n2dry = max( n2dry, 1.e-8 )
        taut(k) = 2.0/( 1.0/300. + 1.0*sqrt(n2) )
!        taut(k) = min( taut(k), taumax )

      enddo
      taut(1) = taut(2)
      taut(gr%nnzp) = taut(gr%nnzp-1)

      return
      end subroutine tau
