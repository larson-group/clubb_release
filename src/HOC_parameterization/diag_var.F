! $Id: diag_var.F,v 1.9 2006-12-20 02:31:07 dschanen Exp $
#define SCLR_THETA 1
#define SCLR_RT 2
!-----------------------------------------------------------------------
      module diagnose_variances

!     Description:
!     Contains subroutine diag_var and ancillaries functions.
!-----------------------------------------------------------------------

        implicit none

        public :: diag_var

        private :: diag_var_lhs
        private :: calc_aa, calc_bb, calc_cc, calc_dd_var, calc_dd_covar

        contains
!-----------------------------------------------------------------------
        subroutine diag_var( taum, wmm, rtm, wprtp, 
     .                       thlm, wpthlp, 
     .                       wp2, wp3, Scm, Skwm,
     .                       rtp2, thlp2, rtpthlp,
     .                       liter, dt, isValid,
     .                       sclrm, wpsclrp, sclrp2,        ! optional
     .                       sclrprtp, sclrpthlp )          ! optional
!       Description:
!       Subprogram to diagnose variances by solving steady-state equations

!       References:
!       Eqn. 13, 14, 15  on p. 3545 of 
!       ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!         Method and Model Description'' Golaz, et al. (2002)
!       JAS, Vol. 59, pp. 3540--3551.

!      See also
!      ``Equations for HOC'', Section 4:
!      /Steady-state solutions for the variances/ 

!-----------------------------------------------------------------------

        use constants
        use grid_class
        use lapack_wrap, only: tridag_solve
#ifdef STATS
        use statistics
#endif /*STATS*/
       
        implicit none

        ! Constants
        integer, parameter :: 
     .  jrtp2    = 1, ! RHS index for rtp2
     .  jthlp2   = 2, ! RHS index for thlp2
     .  jrtpthlp = 3  ! RHS index for rtpthlp

        ! Input variables
        real, intent(in), dimension(gr%nnzp) :: 
     .  taum,  ! Tau on moment. grid            [s]
     .  wmm,   ! w wind on m                    [m/s]
     .  rtm,   ! Total water mixing ratio       [kg/kg]
     .  wprtp, ! w' r_t'                        [(m kg)/(s kg)]
     .  thlm,  ! Liquid potential temp.         [K]
     .  wpthlp,! w' th_l'                       [(m K)/s]
     .  wp2,   ! w'^2                           [m^2/s^2]
     .  wp3,   ! w'^3                           [m^3/s^3]
     .  Scm,   ! Sc on moment. grid             [-]
     .  Skwm   ! Skw on moment. grid            [-]

        logical, intent(in) :: liter ! Whether variances are prognostic

        real, intent(in) :: dt ! Timestep       [s]

        ! Input/Output variables
        ! An attribute of (inout) is also needed to import the
        ! value of the variances at the surface.  Brian.  12/18/05.
        real, intent(inout), dimension(gr%nnzp) :: 
     .  rtp2,   ! r_t'^2                        [(kg/kg)^2]
     .  thlp2,  ! th_l'^2                       [K^2]
     .  rtpthlp ! r_t' th_l'                    [(kg K)/kg]

        ! Output variable for singular matrices

        logical, intent(out) :: isValid

        ! Passive scalar input (optional)
        real, optional, intent(in), dimension(gr%nnzp, sclr_dim) :: 
     .  sclrm, wpsclrp

        ! Passive scalar output (optional)
        real, optional, intent(out), dimension(gr%nnzp, sclr_dim) :: 
     .  sclrp2, sclrprtp, sclrpthlp

        ! Local Variables
        real, dimension(gr%nnzp) :: 
     .  a1 ! a1; See eqn. 20 in `Equations for HOC'
        real, dimension(gr%nnzp) :: 
     .  aa, bb, cc ! sub, main, and super diagonals

        real, dimension(gr%nnzp) :: 
     .  C2_Skw_fnc

        real, allocatable, dimension(:,:) :: 
     .  dd,      ! RHS of tridiagonal system
     .  solution ! Solution to tridiagonal systems

        integer :: nrhs

        ! Location of scalars in dd / solution
        integer, dimension(sclr_dim) :: 
     .  jsclrp2, jsclrprtp, jsclrpthlp 

        logical :: scalar_calc

        ! Loop indices
        integer :: i
        integer :: k, km1, kp1

!-----------------------------------------------------------------------

        if ( present( sclrm ) .and. present( wpsclrp ) .and.
     .       present( sclrp2 ) .and. present( sclrprtp ) .and.
     .       present( sclrpthlp ) ) then

          nrhs = 3 * sclr_dim + 3
          scalar_calc = .true.

          do i = 1, sclr_dim
            jsclrp2(i) = 3 + i
            jsclrprtp(i) = 3 + sclr_dim + i
            jsclrpthlp(i) = 3 + 2*sclr_dim + i
          end do

        else
          nrhs = 3
          scalar_calc = .false.

        end if

        ! rtp2, thlp2, rtpthlp + passive scalars
        allocate( dd(gr%nnzp,nrhs) ) 
        allocate( solution(gr%nnzp,nrhs) )

        ! Define a1
        a1(1:gr%nnzp) = 1.0 / ( 1.0 - Scm(1:gr%nnzp) )

#ifndef DISABLE_C2_SKW
        ! Define C2 as a function of Skw -dschanen 14 April 06
        C2_Skw_fnc(1:gr%nnzp) 
     .  = C2b + (C2 - C2b )
     .    *exp( -(1.0/2.0) * (Skwm(1:gr%nnzp)/C2c)**2 )

#else
        C2_Skw_fnc = C2

#endif /* DISABLE */

#ifdef STATS
        if (lstats_samp) then

          if ( irtp2_bt > 0 ) then
            zm%x(:,irtp2_bt) = zm%x(:,irtp2_bt) - rtp2 / dt
          end if

          if ( ithlp2_bt > 0 ) then
            zm%x(:,ithlp2_bt) = zm%x(:,ithlp2_bt) - thlp2 / dt
          end if

          if ( irtpthlp_bt > 0 ) then
            zm%x(:,irtpthlp_bt) = zm%x(:,irtpthlp_bt) - rtpthlp / dt
          end if

        end if
#endif
        ! Brian added this line of code which can be uncommented
        ! in order to bring back the use of C2rt, C2thl, and 
        ! C2rtthl by overwriting the value of C2_Skw_fnc before
        ! the call to the subroutines for each respective variable.
!        C2_Skw_fnc = 1.5000   ! C2rt

        ! Implicit contributions to terms

        call diag_var_lhs( a1, wp2, wp3, taum, wmm, 
     .                     C2_Skw_fnc, nu2, beta, eps, dt, liter,
     .                     aa, bb, cc )

        ! Explicit contributions to rtp2

        dd(1:gr%nnzp,jrtp2) = 
     .  calc_dd_var( a1, wp2, wp3, wprtp, rtm, 
     .               rtp2, dt, beta, eps, liter )

#ifdef STATS
        if (lstats_samp) then

          do k=1, gr%nnzp

           km1 = max(k-1,1)
           kp1 = min(k+1,gr%nnzp)

           if ( irtp2_ta > 0 ) then
             zm%x(k,irtp2_ta) 
     .       = zm%x(k,irtp2_ta)
     .       - ( 1.0 - (1.0/3.0) * beta )
     .         * (1.0/4.0) * gr%dzm(k)
     .            * (   ( a1(kp1) + a1(k) )**2
     .                  * wp3(kp1) * ( wprtp(kp1) + wprtp(k) )**2
     .                  / ( max( wp2(kp1) + wp2(k), 2.0*eps ) )**2
     .               -  ( a1(k) + a1(km1) )**2
     .                  * wp3(k) * ( wprtp(k) + wprtp(km1) )**2
     .                  / ( max( wp2(k) + wp2(km1), 2.0*eps ) )**2
     .              )
           end if

          end do

        end if ! lstats_samp
#endif
        ! Use old C2rt stuff
!       call tridag_solve( "rtp2", gr%nnzp, 1, cc, bb, aa, dd, 
!    .                     rtp2, isValid )


       ! Brian added this line of code which can be uncommented
       ! in order to bring back the use of C2rt, C2thl, and 
       ! C2rtthl by overwriting the value of C2_Skw_fnc before
       ! the call to the subroutines for each respective variable.
!       C2_Skw_fnc = 1.000   ! C2thl
!       bb = calc_bb( a1, wp2, wp3, taum, 
!     .               C2_Skw_fnc, dt, nu2, beta, eps, liter )

        ! Explicit contributions to thlp2

        dd(1:gr%nnzp,jthlp2) 
     .  = calc_dd_var( a1, wp2, wp3, wpthlp, thlm, 
     .                 thlp2, dt, beta, eps, liter )

#ifdef STATS
       if (lstats_samp) then

         do k=1,gr%nnzp

           km1 = max(k-1,1)
           kp1 = min(k+1,gr%nnzp)

           if ( ithlp2_ta > 0 ) then
             zm%x(k,ithlp2_ta) 
     .       = zm%x(k,ithlp2_ta)
     .       - ( 1.0 - (1.0/3.0) * beta )
     .         * (1.0/4.0) * gr%dzm(k)
     .            * (   ( a1(kp1) + a1(k) )**2
     .                  * wp3(kp1) * ( wpthlp(kp1) + wpthlp(k) )**2
     .                  / ( max( wp2(kp1) + wp2(k), 2.0*eps ) )**2
     .               -  ( a1(k) + a1(km1) )**2
     .                  * wp3(k) * ( wpthlp(k) + wpthlp(km1) )**2
     .                  / ( max( wp2(k) + wp2(km1), 2.0*eps ) )**2
     .              )
           end if ! ithlp2_ta

         end do ! k

       end if ! lstats_samp
#endif
        ! Use old C2thl stuff
!       call tridag_solve( "thlp2", gr%nnzp, 1, cc, bb, aa, dd, 
!    .                     thlp2, isValid )

       ! Brian added this line of code which can be uncommented
       ! in order to bring back the use of C2rt, C2thl, and 
       ! C2rtthl by overwriting the value of C2_Skw_fnc before
       ! the call to the subroutines for each respective variable.
!       C2_Skw_fnc = 2.000   ! C2rtthl
!       bb = calc_bb( a1, wp2, wp3, taum, 
!     .               C2_Skw_fnc, dt, nu2, beta, eps, liter )
       
       ! Explicit contributions to rtpthlp
       dd(:,jrtpthlp) 
     . = calc_dd_covar( a1, wp2, wp3, wprtp, wpthlp, 
     .                  rtm, thlm, rtpthlp, 
     .                  dt, beta, eps, liter ) 

#ifdef STATS

       if (lstats_samp) then

          do k=1, gr%nnzp

             km1 = max(k-1,1)
             kp1 = min(k+1,gr%nnzp)

             if ( irtpthlp_dp1 > 0 ) then
               zmscr01(k) = - C2_Skw_fnc(k) / taum(k)
             endif
         
             if ( irtpthlp_ta > 0 ) then
             zm%x(k,irtpthlp_ta) 
     .       = zm%x(k,irtpthlp_ta)
     .       - ( 1.0 - (1.0/3.0) * beta )
     .         * (1.0/4.0) * gr%dzm(k)
     .            * (  ( a1(kp1) + a1(k) )**2
     .                 * ( wprtp(kp1) + wprtp(k) ) 
     .                 * ( wpthlp(kp1) + wpthlp(k) )
     .                 * wp3(kp1) 
     .                 / ( max( wp2(kp1) + wp2(k), 2.0*eps ) )**2
     .               - ( a1(k) + a1(km1) )**2
     .                 * ( wprtp(k) + wprtp(km1) ) 
     .                 * ( wpthlp(k) + wpthlp(km1) )
     .                 * wp3(k) 
     .                 / ( max( wp2(k) + wp2(km1), 2.0*eps ) )**2
     .              )
             end if ! irtpthlp_ta

          end do ! 1..gr%nnzp

       end if ! lstats_samp
#endif
        ! Use old C2rtthl stuff
!       call tridag_solve( "rtpthlp", gr%nnzp, 1, cc, bb, aa, dd, 
!    .                     rtpthlp, isValid )
       
        if ( scalar_calc ) then

          do i=1, sclr_dim

            ! RHS for sclrp2
            dd(:,jsclrp2(i)) = 
     .      calc_dd_var( a1, wp2, wp3, wpsclrp(:,i), sclrm(:,i), 
     .                   sclrp2(:,i), dt, beta, eps, liter )

            ! RHS for sclrprtp
            dd(:,jsclrprtp(i)) = 
     .      calc_dd_covar( a1, wp2, wp3, wprtp, wpsclrp(:,i),
     .                     rtm, sclrm(:,i), sclrprtp(:,i),
     .                     dt, beta, eps, liter )

            ! RHS for sclrpthlp
            dd(:,jsclrpthlp(i)) = 
     .      calc_dd_covar( a1, wp2, wp3, wpthlp, wpsclrp(:,i),
     .                     thlm, sclrm(:,i), sclrpthlp(:,i),
     .                     dt, beta, eps, liter )

          end do ! 1..sclr_dim

        end if ! present (sclr variables)

        call tridag_solve( "diag_var", gr%nnzp, nrhs, cc, bb, aa, dd, 
     .                     solution, isValid )

        rtp2    = solution(1:gr%nnzp,jrtp2)
        thlp2   = solution(1:gr%nnzp,jthlp2)
        rtpthlp = solution(1:gr%nnzp,jrtpthlp)

#ifdef STATS
        ! Compute rtp2 implicit budget terms
        if (lstats_samp) then

          do k=1, gr%nnzp

           km1 = max(k-1,1)
           kp1 = min(k+1,gr%nnzp)

           if ( irtp2_dp1 > 0 ) then
             zm%x(k,irtp2_dp1) = zm%x(k,irtp2_dp1)
     .       + zmscr01(k) * rtp2(k)
             zm%n(k,irtp2_dp1) = zm%n(k,irtp2_dp1) + 1
           end if
 
           if ( irtp2_dp2 > 0 ) then
             zm%x(k,irtp2_dp2) = zm%x(k,irtp2_dp2)
     .       + zmscr02(k) * rtp2(km1)
     .       + zmscr03(k) * rtp2(k)
     .       + zmscr04(k) * rtp2(kp1)
             zm%n(k,irtp2_dp2) = zm%n(k,irtp2_dp2) + 1
           end if
 
           if ( irtp2_ta > 0 ) then
             zm%x(k,irtp2_ta) = zm%x(k,irtp2_ta)
     .       + zmscr05(k) * rtp2(km1)
     .       + zmscr06(k) * rtp2(k)
     .       + zmscr07(k) * rtp2(kp1)
             zm%n(k,irtp2_ta) = zm%n(k,irtp2_ta) + 1
           end if

           if ( irtp2_ma > 0 ) then
             zm%x(k,irtp2_ma) = zm%x(k,irtp2_ma)
     .       + zmscr08(k) * rtp2(km1)
     .       + zmscr09(k) * rtp2(kp1)
             zm%n(k,irtp2_ma) = zm%n(k,irtp2_ma) + 1
           end if
 
           if ( irtp2_tp > 0 ) then
             zm%x(k,irtp2_tp) = zm%x(k,irtp2_tp)
     .        - 2.0 * wprtp(k) * gr%dzm(k) * ( rtm(kp1) - rtm(k) )
             zm%n(k,irtp2_tp) = zm%n(k,irtp2_tp) + 1
           end if

          end do

       end if
#endif /*STATS*/

#ifdef STATS
       if ( lstats_samp ) then
         ! Compute thlp2 implicit budget terms
         do k=1, gr%nnzp

           km1 = max(k-1,1)
           kp1 = min(k+1,gr%nnzp)

           if ( ithlp2_dp1 > 0 ) then
             zm%x(k,ithlp2_dp1) = zm%x(k,ithlp2_dp1)
     .       + zmscr01(k) * thlp2(k)
             zm%n(k,ithlp2_dp1) = zm%n(k,ithlp2_dp1) + 1
           end if
 
           if ( ithlp2_dp2 > 0 ) then
             zm%x(k,ithlp2_dp2) = zm%x(k,ithlp2_dp2)
     .       + zmscr02(k) * thlp2(km1)
     .       + zmscr03(k) * thlp2(k)
     .       + zmscr04(k) * thlp2(kp1)
             zm%n(k,ithlp2_dp2) = zm%n(k,ithlp2_dp2) + 1
           end if
 
           if ( ithlp2_ta > 0 ) then
             zm%x(k,ithlp2_ta) = zm%x(k,ithlp2_ta)
     .       + zmscr05(k) * thlp2(km1)
     .       + zmscr06(k) * thlp2(k)
     .       + zmscr07(k) * thlp2(kp1)
             zm%n(k,ithlp2_ta) = zm%n(k,ithlp2_ta) + 1
           end if

           if ( ithlp2_ma > 0 ) then
             zm%x(k,ithlp2_ma) = zm%x(k,ithlp2_ma)
     .       + zmscr08(k) * thlp2(km1)
     .       + zmscr09(k) * thlp2(kp1)
             zm%n(k,ithlp2_ma) = zm%n(k,ithlp2_ma) + 1
           end if
 
           if ( ithlp2_tp > 0 ) then
             zm%x(k,ithlp2_tp) = zm%x(k,ithlp2_tp)
     .        - 2.0 * wpthlp(k) * gr%dzm(k) * ( thlm(kp1) - thlm(k) )
             zm%n(k,ithlp2_tp) = zm%n(k,ithlp2_tp) + 1

           end if

         end do ! k

       end if ! lstats
#endif /*STATS*/

#ifdef STATS
        ! Compute rtpthlp implicit budget terms
        if (lstats_samp) then

          do k=1, gr%nnzp

           km1 = max(k-1,1)
           kp1 = min(k+1,gr%nnzp)

           if ( irtpthlp_dp1 > 0 ) then
             zm%x(k,irtpthlp_dp1) = zm%x(k,irtpthlp_dp1)
     .       + zmscr01(k) * rtpthlp(k)
             zm%n(k,irtpthlp_dp1) = zm%n(k,irtpthlp_dp1) + 1
           end if
 
           if ( irtpthlp_dp2 > 0 ) then
             zm%x(k,irtpthlp_dp2) = zm%x(k,irtpthlp_dp2)
     .       + zmscr02(k) * rtpthlp(km1)
     .       + zmscr03(k) * rtpthlp(k)
     .       + zmscr04(k) * rtpthlp(kp1)
             zm%n(k,irtpthlp_dp2) = zm%n(k,irtpthlp_dp2) + 1
           end if
 
           if ( irtpthlp_ta > 0 ) then
             zm%x(k,irtpthlp_ta) = zm%x(k,irtpthlp_ta)
     .       + zmscr05(k) * rtpthlp(km1)
     .       + zmscr06(k) * rtpthlp(k)
     .       + zmscr07(k) * rtpthlp(kp1)
             zm%n(k,irtpthlp_ta) = zm%n(k,irtpthlp_ta) + 1
           end if

           if ( irtpthlp_ma > 0 ) then
             zm%x(k,irtpthlp_ma) = zm%x(k,irtpthlp_ma)
     .       + zmscr08(k) * rtpthlp(km1)
     .       + zmscr09(k) * rtpthlp(kp1)
             zm%n(k,irtpthlp_ma) = zm%n(k,irtpthlp_ma) + 1
           end if
 
           if ( irtpthlp_tp1 > 0 ) then
             zm%x(k,irtpthlp_tp1) = zm%x(k,irtpthlp_tp1)
     .       - wprtp(k) * gr%dzm(k) * ( thlm(kp1) - thlm(k) )
             zm%n(k,irtpthlp_tp1) = zm%n(k,irtpthlp_tp1) + 1
           end if

           if ( irtpthlp_tp2 > 0 ) then
             zm%x(k,irtpthlp_tp2) = zm%x(k,irtpthlp_tp2)
     .       - wpthlp(k) * gr%dzm(k) * ( rtm(kp1) - rtm(k) )
             zm%n(k,irtpthlp_tp2) = zm%n(k,irtpthlp_tp2) + 1
           end if

         end do
       end if
#endif /*STATS*/

#ifdef STATS
        if (lstats_samp .and. irtp2_cl > 0 ) then
          zm%x(:,irtp2_cl) = zm%x(:,irtp2_cl) - rtp2 / dt
        end if
        if (lstats_samp .and. ithlp2_cl > 0 ) then
          zm%x(:,ithlp2_cl) = zm%x(:,ithlp2_cl) - thlp2 / dt
        end if
#endif

        ! Clip small and negative values
        ! The surface variances are set elsewhere.  They do not need to be
        ! clipped at the lowest level.  Brian Griffin.  12/18/05.
        do k = 2, gr%nnzp, 1

          if ( wp2(k) >= wtol*wtol ) then
            rtp2(k) 
     .      = max( rtp2(k), 
     .             wprtp(k)*wprtp(k)/((1.0-Scm(k))*wp2(k)) + eps,
     .             rttol*rttol + eps )
            thlp2(k) 
     .      = max( thlp2(k), 
     .             wpthlp(k)*wpthlp(k)/((1.0-Scm(k))*wp2(k)) + eps,
     .             thltol*thltol + eps )
          else
            thlp2(k) = max( thlp2(k), 0.0 )
            rtp2(k)  = max( rtp2(k), 0.0 )
          end if

        end do ! k=2..gr%nnzp

#ifdef STATS
        if (lstats_samp .and. irtp2_cl > 0 ) then
          zm%x(:,irtp2_cl) = zm%x(:,irtp2_cl) + rtp2 / dt
          zm%n(:,irtp2_cl) = zm%n(:,irtp2_cl) + 1
        end if

        if (lstats_samp .and. ithlp2_cl > 0 ) then
          zm%x(:,ithlp2_cl) = zm%x(:,ithlp2_cl) + thlp2 / dt
          zm%n(:,ithlp2_cl) = zm%n(:,ithlp2_cl) + 1
        end if
#endif

#ifdef STATS
        if (lstats_samp) then

          if ( irtp2_bt > 0 ) then
            zm%x(:,irtp2_bt) = zm%x(:,irtp2_bt) + rtp2 / dt
            zm%n(:,irtp2_bt) = zm%n(:,irtp2_bt) + 1
          end if

          if ( ithlp2_bt > 0 ) then
            zm%x(:,ithlp2_bt) = zm%x(:,ithlp2_bt) + thlp2 / dt
            zm%n(:,ithlp2_bt) = zm%n(:,ithlp2_bt) + 1
          end if

          if ( irtpthlp_bt > 0 ) then
            zm%x(:,irtpthlp_bt) = zm%x(:,irtpthlp_bt) + rtpthlp / dt
            zm%n(:,irtpthlp_bt) = zm%n(:,irtpthlp_bt) + 1
          end if

        end if
#endif /*STATS*/

        if ( scalar_calc ) then
          sclrp2    = solution(1:gr%nnzp,jsclrp2)
          sclrpthlp = solution(1:gr%nnzp,jsclrpthlp)
          sclrprtp  = solution(1:gr%nnzp,jsclrprtp)

          ! Clip scalar variance terms
          do i =1, sclr_dim, 1
            do k = 2, gr%nnzp, 1
              if ( wp2(k) >= wtol*wtol ) then
                sclrp2(k,i) 
     .          = max
     .            ( sclrp2(k,i), 
     .              wpsclrp(k,i)*wpsclrp(k,i)/((1.0-Scm(k))*wp2(k))+eps,
     .              sclrtol(i)*sclrtol(i) + eps )
              else
                sclrp2(k,i) = max( sclrp2(k,i), 0.0 )
              end if ! wp2(k) >= wtol**2
            end do ! k
          end do ! i

        end if ! scalar_calc

        deallocate( dd, solution )

        return
        end subroutine diag_var

!-----------------------------------------------------------------------
        subroutine diag_var_lhs( a1, wp2, wp3, taum, wmm, 
     .                           C2x, nu2, beta, eps, dt, liter,
     .                           aa, bb, cc )
!       Description:
!       Compute LHS banded matrix

!       References:
!       None
!-----------------------------------------------------------------------
        use grid_class
#ifdef STATS
        use statistics
#endif

        implicit none

        ! Input Variables
        real, dimension(gr%nnzp), intent(in) ::
     .  a1,    ! See eqn. 20 in `Equations for HOC'     [-]
     .  wp2,   ! w'^2                                   [m^2/s^2]
     .  wp3,   ! w'^3                                   [m^3/s^3]
     .  taum,  ! Tau on moment. grid                    [s]
     .  wmm,   ! w wind on m.                           [m/s]
     .  C2x    ! C2 (possibly a function of Skw)        [-]

        real, intent(in) ::
     .  dt,    ! Timestep length                        [s]
     .  nu2,   ! Previously determined eddy viscosity.  [-]
     .  beta,  ! Constant parameter                     [-]
     .  eps    ! Small value to prevent a divide by zero[-]

        logical, intent(in) ::
     .  liter  ! Whether the variances are prognostic

        ! Output Variables
        real, dimension(gr%nnzp), intent(out) ::
     .  aa,! Sub diagonal
     .  bb,! Main diagonal
     ,  cc ! Super diagonal

        ! Local Variables
        integer :: k, kp1, km1

        ! Prepare tridiagonal system
        !  Where cc, bb & aa contain the LHS super, main, and sub diagonal,
        !  respectively, and dd contains the RHS(s)

        aa = calc_aa( a1, wp2, wp3, wmm, nu2, beta, eps )
        bb = calc_bb( a1, wp2, wp3, taum, 
     .                C2x, dt, nu2, beta, eps, liter )
        cc = calc_cc( a1, wp2, wp3, wmm, nu2, beta, eps )

#ifdef STATS
       if ( lstats_samp ) then

         do k=1, gr%nnzp

           km1   = max(k-1,1)
           kp1   = min(k+1,gr%nnzp)

           if ( irtp2_dp1 * ithlp2_dp1 * irtpthlp_dp1 > 0 ) then
             zmscr01(k) = - C2x(k) / taum(k)
           end if
         
           if ( irtp2_dp2 * ithlp2_dp2 * irtpthlp_dp2 > 0 ) then
             zmscr02(k) = + nu2 * gr%dzm(k) * gr%dzt(k)
             zmscr03(k) = 
     .       - nu2 * gr%dzm(k) * ( gr%dzt(kp1) + gr%dzt(k) )
             zmscr04(k) = + nu2 * gr%dzm(k) * gr%dzt(kp1)
           end if

           if ( irtp2_ta * ithlp2_ta * irtpthlp_ta > 0 ) then
             zmscr05(k) = 
     .       + (1.0/3.0) * beta
     .         * (1.0/2.0) * gr%dzm(k)
     .            * ( a1(k) + a1(km1) )
     .            * ( wp3(k) / max( wp2(k) + wp2(km1), 2.0*eps ) )
             zmscr06(k) = 
     .       - (1.0/3.0) * beta
     .         * (1.0/2.0) * gr%dzm(k)
     .            * (   ( a1(kp1) + a1(k) )
     .                * ( wp3(kp1) / max( wp2(kp1) + wp2(k), 2.0*eps ) )
     .               -  ( a1(k) + a1(km1) )
     .                * ( wp3(k) / max( wp2(k) + wp2(km1), 2.0*eps ) )
     .              )
             zmscr07(k) = 
     .       - (1.0/3.0) * beta
     .         * (1.0/2.0) * gr%dzm(k)
     .            * ( a1(kp1) + a1(k) )
     .            * ( wp3(kp1) / max( wp2(kp1) + wp2(k), 2.0*eps ) )
           end if

           if ( irtp2_ma * ithlp2_ma * irtpthlp_ma > 0 ) then
             zmscr08(k) = + wmm(k) * (1.0/2.0) * gr%dzm(k)
             zmscr09(k) = - wmm(k) * (1.0/2.0) * gr%dzm(k)
           end if

         end do ! k=1,gr%nnzp

       end if
#endif /*STATS*/

        return
        end subroutine diag_var_lhs

!-----------------------------------------------------------------------
        function calc_aa( a1, wp2, wp3, wmm, nu2, beta, eps )

!       Description:
!       Calculate subdiagonal

!       Notes:
!       This only needs to be calculated once per diag_var call
!-----------------------------------------------------------------------
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

           km1 = max(k-1,1)
           kp1 = min(k+1, gr%nnzp)

           calc_aa(k) = 
     .     - wmm(k) * (1.0/2.0) * gr%dzm(k)
     .     - (1.0/3.0) * beta
     .       * (1.0/2.0) * gr%dzm(k)
     .          * ( a1(k) + a1(km1) )
     .          * ( wp3(k) / max( wp2(k) + wp2(km1), 2.0*eps ) )
     .     - nu2 * gr%dzm(k) * gr%dzt(k)

        end do

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
!  Note: new C2_Skw need to be generated only once
!-----------------------------------------------------------------------
        function calc_bb( a1, wp2, wp3, taum,
     .                    C2x, dt, nu2, beta, eps, liter )
        use grid_class

        implicit none

!       Input
        real, intent(in), dimension(gr%nnzp) :: a1, wp2, wp3, taum, C2x

        real, intent(in) :: dt, beta, nu2, eps

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
     .              * ( wp3(kp1) / max( wp2(kp1) + wp2(k), 2.0*eps ) )
     .             -  ( a1(k) + a1(km1) )
     .              * ( wp3(k) / max( wp2(k) + wp2(km1), 2.0*eps ) )
     .            )
     .     + C2x(k) / taum(k)
     .     + nu2 * gr%dzm(k) * ( gr%dzt(kp1) + gr%dzt(k) )

           if (liter) calc_bb(k) = calc_bb(k) + 1.0/dt

        end do

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
        function calc_cc( a1, wp2, wp3, wmm, nu2, beta, eps )

!       Description:
!       Calculate the superdiagonal

!       Notes:
!       This only needs to be calculated once per diag_var call
!-----------------------------------------------------------------------

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
     .          * ( wp3(kp1) / max( wp2(kp1) + wp2(k), 2.0*eps ) )
     .     - nu2 * gr%dzm(k) * gr%dzt(kp1)

        end do

        ! Boundary Conditions
        ! These are set so that the sfc_var value of rtp2, thlp2, or rtpthlp
        ! can be used at the lowest boundary and the values of those
        ! variables can be set to 0 at the top boundary.  Brian.  12/18/05.
        calc_cc(1) = 0.0
        calc_cc(gr%nnzp) = 0.0

        return
        end function calc_cc

!-----------------------------------------------------------------------
        function calc_dd_covar( a1, wp2, wp3, wpxap, wpxbp, 
     .                          xam, xbm, xapxbp, 
     .                          dt, beta, eps, liter )
!       Description:
!       Calculate RHS for the purpose of solving the tridiagonal system 
!       of a covariance term (e.g rtpthlp).
!-----------------------------------------------------------------------

        use grid_class

        implicit none

        ! Input
        real, intent(in), dimension(gr%nnzp) :: a1, wp2, wp3
        real, intent(in), dimension(gr%nnzp) :: wpxap, wpxbp
        real, intent(in), dimension(gr%nnzp) :: xam, xbm
        real, intent(in), dimension(gr%nnzp) :: xapxbp

        real, intent(in) :: dt, beta, eps

        logical, intent(in) :: liter

        ! Return type
        real, dimension(gr%nnzp) :: calc_dd_covar

        ! Local Variables
        integer k, km1, kp1

        do k = 2, gr%nnzp-1, 1

          km1 = max(k-1,1)
          kp1 = min(k+1, gr%nnzp)

          calc_dd_covar(k) = 
     .    - ( 1.0 - (1.0/3.0) * beta )
     .      * (1.0/4.0) * gr%dzm(k)
     .         * (  ( a1(kp1) + a1(k) )**2
     .              * ( wpxap(kp1) + wpxap(k) ) 
     .              * ( wpxbp(kp1) + wpxbp(k) )
     .              * wp3(kp1) 
     .              / ( max( wp2(kp1) + wp2(k), 2.0*eps ) )**2
     .            - ( a1(k) + a1(km1) )**2
     .              * ( wpxap(k) + wpxap(km1) ) 
     .              * ( wpxbp(k) + wpxbp(km1) )
     .              * wp3(k) 
     .              / ( max( wp2(k) + wp2(km1), 2.0*eps ) )**2
     .           )
     .    - wpxbp(k) * gr%dzm(k) * ( xam(kp1) - xam(k) )
     .    - wpxap(k) * gr%dzm(k) * ( xbm(kp1) - xbm(k) )

          if (liter) calc_dd_covar(k) = calc_dd_covar(k) 
     .                                + 1.0/dt*xapxbp(k)

        end do

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
        function calc_dd_var( a1, wp2, wp3, wpxp, xm,
     .                        xp2, dt, beta, eps, liter )

!       Description:
!       Calculate dd for the purpose of solving the tridiagonal system 
!       of a variance term (e.g. rtp2)
!-----------------------------------------------------------------------

        use grid_class

        implicit none

        ! Input Variables
        real, intent(in), dimension(gr%nnzp) :: a1, wp2, wp3
        real, intent(in), dimension(gr%nnzp) :: wpxp, xm, xp2

        real, intent(in) :: dt, beta, eps

        logical, intent(in) :: liter

        ! Return type
        real, dimension(gr%nnzp) :: calc_dd_var

        ! Local Variables
        integer k, km1, kp1

        do k = 2, gr%nnzp-1, 1

           km1 = max(k-1,1)
           kp1 = min(k+1, gr%nnzp)

           calc_dd_var(k) =
     .     - ( 1.0 - (1.0/3.0) * beta )
     .       * (1.0/4.0) * gr%dzm(k)
     .          * (   ( a1(kp1) + a1(k) )**2
     .                * wp3(kp1) * ( wpxp(kp1) + wpxp(k) )**2
     .                / ( max( wp2(kp1) + wp2(k), 2.0*eps ) )**2
     .             -  ( a1(k) + a1(km1) )**2
     .                * wp3(k) * ( wpxp(k) + wpxp(km1) )**2
     .                / ( max( wp2(k) + wp2(km1), 2.0*eps ) )**2
     .            )
     .     - 2.0 * wpxp(k) * gr%dzm(k) * ( xm(kp1) - xm(k) )

           if (liter) calc_dd_var(k) = calc_dd_var(k) + 1.0/dt*xp2(k)

        end do

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
