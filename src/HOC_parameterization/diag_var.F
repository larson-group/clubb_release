! $Id: diag_var.F,v 1.11 2007-01-30 17:59:07 dschanen Exp $
#define SCLR_THETA 1
#define SCLR_RT 2
!-----------------------------------------------------------------------
      module diagnose_variances

!     Description:
!     Contains subroutine diag_var and ancillaries functions.
!-----------------------------------------------------------------------

        implicit none

        public :: diag_var

        private :: diag_var_lhs, diag_var_solve

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

        ! Constants parameters

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
     .  C2sclr_1d, C2rt_1d, C2thl_1d, C2rtthl_1d

        real, dimension(gr%nnzp) :: 
     .  a1 ! a1; See eqn. 20 in `Equations for HOC'

        real, dimension(3,gr%nnzp) :: 
     .  lhs ! Tridiagonal matrix

        real, dimension(gr%nnzp,1) ::
     .  rhs ! RHS vector of Tridiagonal matrix

        real, dimension(gr%nnzp,sclr_dim*3) :: 
     .  sclr_rhs,     ! RHS of scalar tridiagonal system
     .  sclr_solution ! Solution to tridiagonal system

        logical, dimension(4) ::
     .  Valid_arr

        logical :: scalar_calc

        ! Loop indices
        integer :: i
        integer :: k, km1, kp1


!-----------------------------------------------------------------------
#ifdef SINGLE_C2_SKW

        C2rt_1d(1:gr%nnzp) 
     .  = C2b + (C2-C2b) *exp( -0.5 * (Skwm(1:gr%nnzp)/C2c)**2 )

        C2thl_1d = C2rt_1d
        C2rtthl_1d = C2rt_1d

        C2sclr_1d  = C2rt_1d
#else
        C2rt_1d(1:gr%nnzp)    = C2rt
        C2thl_1d(1:gr%nnzp)   = C2thl
        C2rtthl_1d(1:gr%nnzp) = C2rtthl

        C2sclr_1d(1:gr%nnzp)  = C2rt  ! Use rt value for now
#endif

        if ( present( sclrm ) .and. present( wpsclrp ) .and.
     .       present( sclrp2 ) .and. present( sclrprtp ) .and.
     .       present( sclrpthlp ) ) then

          scalar_calc = .true.

        else
          scalar_calc = .false.

        end if

        ! Define a1
        a1(1:gr%nnzp) = 1.0 / ( 1.0 - Scm(1:gr%nnzp) )

#ifdef STATS
        if ( lstats_samp ) then

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
        ! Initialize to true
        Valid_arr(:) = .true.

        ! Implicit contributions to term rtp2

        call diag_var_lhs( a1, wp2, wp3, taum, wmm, 
     .                     C2rt_1d, nu2, beta, eps, dt, liter,
     .                     lhs )

        ! Explicit contributions to rtp2
        call diag_var_rhs( "rtp2", a1, wp2, wp3, 
     .                     rtm, rtm,  wprtp, wprtp, rtp2,
     .                     dt, beta, eps, liter, rhs )
        
        ! Solve the tridiagonal system

        call diag_var_solve( "rtp2", 1, rhs, 
     .                       lhs, rtp2, Valid_arr(1) )

        ! Implicit contributions to term thlp2

        call diag_var_lhs( a1, wp2, wp3, taum, wmm, 
     .                     C2thl_1d, nu2, beta, eps, dt, liter,
     .                     lhs )

        ! Explicit contributions to thlp2

        call diag_var_rhs( "thlp2", a1, wp2, wp3, 
     .                     thlm, thlm, wpthlp, wpthlp, thlp2,
     .                     dt, beta, eps, liter, rhs )

!       rhs(:,1) = calc_dd_covar( a1, wp2, wp3, wpthlp, wpthlp,
!    .                       thlm, thlm, thlp2,
!    .                       dt, beta, eps, liter )

        ! Solve the tridiagonal system

        call diag_var_solve( "thlp2", 1, rhs, 
     .                       lhs, thlp2, Valid_arr(2) )

        ! Implicit contributions to term rtpthlp

        call diag_var_lhs( a1, wp2, wp3, taum, wmm, 
     .                     C2rtthl_1d, nu2, beta, eps, dt, liter,
     .                     lhs )

        ! Explicit contributions to rtpthlp

        call diag_var_rhs( "rtpthlp", a1, wp2, wp3, 
     .                     rtm, thlm, wprtp, wpthlp, rtpthlp,
     .                     dt, beta, eps, liter, rhs )

!       rhs(:,1) = calc_dd_covar( a1, wp2, wp3, wprtp, wpthlp,
!    .                       rtm, thlm, rtpthlp,
!    .                       dt, beta, eps, liter )

        ! Solve the tridiagonal system

        call diag_var_solve( "rtpthlp", 1, rhs, 
     .                       lhs, rtpthlp, Valid_arr(3) )


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
        if ( lstats_samp .and. irtp2_cl > 0 ) then
          zm%x(:,irtp2_cl) = zm%x(:,irtp2_cl) + rtp2 / dt
          zm%n(:,irtp2_cl) = zm%n(:,irtp2_cl) + 1
        end if

        if ( lstats_samp .and. ithlp2_cl > 0 ) then
          zm%x(:,ithlp2_cl) = zm%x(:,ithlp2_cl) + thlp2 / dt
          zm%n(:,ithlp2_cl) = zm%n(:,ithlp2_cl) + 1
        end if
#endif

#ifdef STATS
        if ( lstats_samp ) then

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

        ! Implicit contributions to passive scalars

          call diag_var_lhs( a1, wp2, wp3, taum, wmm, 
     .                       C2sclr_1d, nu2, beta, eps, dt, liter,
     .                       lhs )


          ! Explicit contributions to passive scalars
          do i = 1, sclr_dim, 1
            call diag_var_rhs
     .           ( "sclrp2", a1, wp2, wp3, 
     .             sclrm(:,i), sclrm(:,i), wpsclrp(:,i), wpsclrp(:,i), 
     .             sclrp2(:,i), dt, beta, eps, liter, sclr_rhs(:,i) )

          end do

          do i = 1, sclr_dim, 1
            call diag_var_rhs
     .           ( "sclrprtp", a1, wp2, wp3, 
     .             sclrm(:,i), rtm, wpsclrp(:,i), wprtp, 
     .             sclrprtp(:,i), dt, beta, eps, liter, 
     .             sclr_rhs(:,i+sclr_dim) )
          end do

          do i = 1, sclr_dim, 1
            call diag_var_rhs
     .           ( "sclrpthlp", a1, wp2, wp3, 
     .             sclrm(:,i), thlm, wpsclrp(:,i), wpthlp, 
     .             sclrprtp(:,i), dt, beta, eps, liter, 
     .             sclr_rhs(:,i+2*sclr_dim) )
          end do

          ! Solve the tridiagonal system

          call diag_var_solve
     .         ( "scalars", 3*sclr_dim, sclr_rhs, 
     .           lhs, sclr_solution, Valid_arr(4) )

          sclrp2(:,1:sclr_dim) = sclr_solution(:,1:sclr_dim)

          sclrprtp(:,1:sclr_dim) 
     .    = sclr_solution(:,sclr_dim+1:2*sclr_dim)

          sclrpthlp(:,1:sclr_dim) 
     .    = sclr_solution(:,2*sclr_dim+1:3*sclr_dim)

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

        ! Check for singular matrices
        if ( any( .not. Valid_arr(:) ) ) then
          isValid = .false.
        end if

        return
        end subroutine diag_var

!-----------------------------------------------------------------------
        subroutine diag_var_lhs( a1, wp2, wp3, taum, wmm, 
     .                           C2, nu2, beta, eps, dt, liter,
     .                           lhs )
!       Description:
!       Compute LHS tridiagonal matrix for a variance or coveriance term

!       References:
!       None
!-----------------------------------------------------------------------
        use grid_class
#ifdef STATS
        use statistics
#endif

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  kp1diag = 1,
     .  kdiag   = 2,
     .  km1diag = 3


        ! Input Variables
        real, dimension(gr%nnzp), intent(in) ::
     .  a1,    ! See eqn. 20 in `Equations for HOC'     [-]
     .  wp2,   ! w'^2                                   [m^2/s^2]
     .  wp3,   ! w'^3                                   [m^3/s^3]
     .  taum,  ! Tau on moment. grid                    [s]
     .  wmm,   ! w wind on m.                           [m/s]
     .  C2     ! C2 coefficient                         [-]

        real, intent(in) ::
     .  dt,    ! Timestep length                        [s]
     .  nu2,   ! Previously determined eddy viscosity.  [-]
     .  beta,  ! Constant parameter                     [-]
     .  eps    ! Small value to prevent a divide by zero[-]

        logical, intent(in) ::
     .  liter  ! Whether the variances are prognostic

        ! Output Variables
        real, dimension(3,gr%nnzp), intent(out) ::
     .  lhs ! Implicit contributions to the term

        ! Local Variables
        real, dimension(3) ::
     .  tmp

        ! Array indices
        integer :: k, kp1, km1

        ! Setup LHS of the tridiagonal system
        do k = 2, gr%nnzp-1, 1
          km1 = max( k-1, 1 )
          kp1 = min( k+1, gr%nnzp )

          lhs(kp1diag:km1diag,k)
     .    = term_dp2( nu2, gr%dzm(km1), gr%dzm(k), gr%dzt(k) ) ! dp2

!         lhs(kdiag,k)
!    .    = lhs(kdiag,k) + term_dp1( C2(k), taum(k) )  ! dp1
          lhs(kdiag,k)
     .    = lhs(kdiag,k) + C2(k) / taum(k) ! dp1

          if ( liter ) then
            lhs(kdiag,k) = lhs(kdiag,k) + ( 1.0 / dt ) ! Time tendency
          end if

          lhs((/kp1diag,km1diag/),k)
     .    = lhs((/kp1diag,km1diag/),k)
     .    + term_ma( wmm(k), gr%dzm(k) ) ! ma

          lhs(kp1diag:km1diag,k)
     .    = lhs(kp1diag:km1diag,k)
     .    + term_ta_lhs
     .      ( a1(kp1), a1(k), a1(km1), wp3(kp1), wp3(k),
     .        wp2(kp1), wp2(k), wp2(km1), gr%dzm(k), beta, eps ) ! ta

#ifdef STATS
         if ( lstats_samp ) then

           if ( irtp2_dp1 * ithlp2_dp1 * irtpthlp_dp1 > 0 ) then
!            tmp(1) = term_dp1( C2(k), taum(k) )
             tmp(1) = C2(k) / taum(k)
             zmscr01(k) = - tmp(1)
           end if
         
           if ( irtp2_dp2 * ithlp2_dp2 * irtpthlp_dp2 > 0 ) then
             tmp(1:3)
     .       = term_dp2( nu2, gr%dzm(km1), gr%dzm(k), gr%dzt(k) )
             zmscr02(k) = -tmp(3)
             zmscr03(k) = -tmp(2)
             zmscr04(k) = -tmp(1)
           end if

           if ( irtp2_ta * ithlp2_ta * irtpthlp_ta > 0 ) then
             tmp(1:3)
     .       = term_ta_lhs
     .         ( a1(kp1), a1(k), a1(km1), wp3(kp1), wp3(k),
     .           wp2(kp1), wp2(k), wp2(km1), gr%dzm(k), beta, eps )
             zmscr05(k) = -tmp(3)
             zmscr06(k) = -tmp(2)
             zmscr07(k) = -tmp(1)
           end if

           if ( irtp2_ma * ithlp2_ma * irtpthlp_ma > 0 ) then
             tmp(1:2)
     .       = term_ma( wmm(k), gr%dzm(k) )
             zmscr08(k) = -tmp(2)
             zmscr09(k) = -tmp(1)
           end if

       end if
#endif /*STATS*/

        end do ! k=2..gr%nnzp-1

        ! Boundary Conditions
        ! These are set so that the sfc_var value of rtp2, thlp2, or rtpthlp
        ! can be used at the lowest boundary and the values of those
        ! variables can be set to 0 at the top boundary.  Brian.  12/18/05.
        lhs(:,1) = 0.0
        lhs(:,gr%nnzp) = 0.0

        lhs(kdiag,1) = 1.0
        lhs(kdiag,gr%nnzp) = 1.0

        if ( liter ) then
          lhs(kdiag,1) = 1.0/dt
          lhs(kdiag,1) = 1.0/dt
        end if



        return
        end subroutine diag_var_lhs

!-----------------------------------------------------------------------
        subroutine diag_var_solve( solve_type, nrhs, rhs, 
     .                             lhs, xapxbp, isValid )

!-----------------------------------------------------------------------
        use lapack_wrap, only: tridag_solve
        use grid_class
#ifdef STATS
        use statistics
#endif

        ! Constant parameters
        integer, parameter ::
     .  kp1diag = 1,
     .  kdiag   = 2,
     .  km1diag = 3

        ! Input variables
        integer, intent(in) :: nrhs  ! Number of right hand side vectors

        character(len=*), intent(in) :: solve_type

        ! Input/Ouput variables
        real, dimension(3,gr%nnzp), intent(inout) ::
     .  lhs     ! Implicit contributions to x variance/covariance term

        real, dimension(gr%nnzp,nrhs), intent(inout) ::
     .  rhs   

        ! Output variables
        real, dimension(gr%nnzp,nrhs), intent(inout) :: 
     .  xapxbp

        logical, intent(out) ::
     .  isValid ! Returns as false in the event of a singular matrix

        integer :: k, kp1, km1

#ifdef STATS
        integer ::
     .  ixapxbp_dp1,
     .  ixapxbp_dp2,
     .  ixapxbp_ta,
     .  ixapxbp_ma

        select case( trim( solve_type ) )
        case( "rtp2" )
          ixapxbp_dp1 = irtp2_dp1
          ixapxbp_dp2 = irtp2_dp2
          ixapxbp_ta  = irtp2_ta
          ixapxbp_ma  = irtp2_ma
        case( "thlp2" )
          ixapxbp_dp1 = ithlp2_dp1
          ixapxbp_dp2 = ithlp2_dp2
          ixapxbp_ta  = ithlp2_ta
          ixapxbp_ma  = ithlp2_ma
        case( "rtpthlp" )
          ixapxbp_dp1 = irtpthlp_dp1
          ixapxbp_dp2 = irtpthlp_dp2
          ixapxbp_ta  = irtpthlp_ta
          ixapxbp_ma  = irtpthlp_ma
        case default ! No budgets for passive scalars
          ixapxbp_dp1 = 0
          ixapxbp_dp2 = 0
          ixapxbp_ta  = 0
          ixapxbp_ma  = 0
        end select
#endif
        ! Solve the system
        call tridag_solve
     .       ( solve_type, gr%nnzp, nrhs, lhs(kp1diag,:), lhs(kdiag,:), 
     .         lhs(km1diag,:), rhs(:,1:nrhs), xapxbp(:,1:nrhs), 
     .         isValid )

#ifdef STATS
        ! Compute implicit budget terms
        if ( lstats_samp ) then

          do k=1, gr%nnzp

            km1 = max( k-1, 1 )
            kp1 = min( k+1, gr%nnzp )

            if ( ixapxbp_dp1 > 0 ) then
              zm%x(k,ixapxbp_dp1) = zm%x(k,ixapxbp_dp1)
     .        + zmscr01(k) * xapxbp(k,1)
              zm%n(k,ixapxbp_dp1) = zm%n(k,ixapxbp_dp1) + 1
            end if
 
            if ( ixapxbp_dp2 > 0 ) then
              zm%x(k,ixapxbp_dp2) = zm%x(k,ixapxbp_dp2)
     .        + zmscr02(k) * xapxbp(km1,1)
     .        + zmscr03(k) * xapxbp(k,1)
     .        + zmscr04(k) * xapxbp(kp1,1)
              zm%n(k,ixapxbp_dp2) = zm%n(k,ixapxbp_dp2) + 1
            end if
 
            if ( ixapxbp_ta > 0 ) then
              zm%x(k,ixapxbp_ta) = zm%x(k,ixapxbp_ta)
     .        + zmscr05(k) * xapxbp(km1,1)
     .        + zmscr06(k) * xapxbp(k,1)
     .        + zmscr07(k) * xapxbp(kp1,1)
              zm%n(k,ixapxbp_ta) = zm%n(k,ixapxbp_ta) + 1
            end if

            if ( ixapxbp_ma > 0 ) then
              zm%x(k,ixapxbp_ma) = zm%x(k,ixapxbp_ma)
     .        + zmscr08(k) * xapxbp(km1,1)
     .        + zmscr09(k) * xapxbp(kp1,1)
              zm%n(k,ixapxbp_ma) = zm%n(k,ixapxbp_ma) + 1
            end if

          end do
        end if
#endif /*STATS*/

        return
        end subroutine diag_var_solve

!-----------------------------------------------------------------------
        subroutine diag_var_rhs( solve_type, a1, wp2, wp3, 
     .                           xam, xbm,  wpxap, wpxbp, xapxbp,
     .                           dt, beta, eps, liter, rhs )
!       Description:

!-----------------------------------------------------------------------
        use grid_class
#ifdef STATS
        use statistics
#endif

        implicit none


        ! Input Variables
        character(len=*), intent(in) :: solve_type

        real, dimension(gr%nnzp), intent(in) ::
     .  a1,     ! a1            [-]
     .  wp2,    ! w'^2          [m^2/s^2]
     .  wp3     ! w'^3          [m^3/s^3]

        real, dimension(gr%nnzp), intent(in) ::
     .  wpxap,  ! w'x'(1)       [units vary]
     .  wpxbp,  ! w'x'(2)       [units vary]
     .  xam,    ! x(1)          [units vary]
     .  xbm,    ! x(2)          [units vary]
     .  xapxbp  ! x'(1) x'(2)   [units vary]

        real, intent(in) ::
     .  dt      ! Timestep      [s]

        real, intent(in) ::
     .  beta,  ! Model parameter                [-]
     .  eps    ! Small value to prevent divz    [-]

        logical, intent(in) ::
     .  liter   ! Whether x is prognostic (T/F)

        real, dimension(gr%nnzp,1), intent(out) ::
     .  rhs     ! Explicit contributions to x variance/covariance terms

        ! Local Variables
        integer :: k, kp1, km1 ! Array indices

#ifdef STATS
        integer ::
     .  ixapxbp_ta,
     .  ixapxbp_tp,
     .  ixapxbp_tp1,
     .  ixapxbp_tp2

        select case( trim( solve_type ) )
        case( "rtp2" )
          ixapxbp_ta  = irtp2_ta
          ixapxbp_tp  = irtp2_tp
          ixapxbp_tp1 = 0
          ixapxbp_tp2 = 0
        case( "thlp2" )
          ixapxbp_ta  = ithlp2_ta
          ixapxbp_tp  = ithlp2_tp
          ixapxbp_tp1 = 0
          ixapxbp_tp2 = 0
        case( "rtpthlp" )
          ixapxbp_ta  = irtpthlp_ta
          ixapxbp_tp  = 0
          ixapxbp_tp1 = irtpthlp_tp1
          ixapxbp_tp2 = irtpthlp_tp2
        case default ! No budgets for passive scalars
          ixapxbp_ta  = 0
          ixapxbp_tp  = 0
          ixapxbp_tp1 = 0
          ixapxbp_tp2 = 0
        end select
#endif
        do k = 2, gr%nnzp-1, 1

          km1 = max( k-1, 1 )
          kp1 = min( k+1, gr%nnzp )

          rhs(k,1)
     .    = term_ta_rhs
     .      ( a1(kp1), a1(k), a1(km1), wp3(kp1), wp3(k),
     .        wp2(kp1), wp2(k), wp2(km1), 
     .        wpxbp(kp1), wpxbp(k), wpxbp(km1),
     .        wpxap(kp1), wpxap(k), wpxap(km1),
     .        gr%dzm(k), beta, eps )
     .    + term_tp( xam(kp1), xam(k), xbm(kp1), xbm(k), 
     .               wpxbp(k), wpxap(k), gr%dzm(k) )

          if ( liter ) then 
            rhs(k,1) = rhs(k,1) + 1.0/dt*xapxbp(k)  ! Time tendency
          end if

        end do ! k=2..gr%nnzp-1

        ! Boundary Conditions
        ! These are set so that the sfc_var value of rtp2, thlp2, or rtpthlp
        ! can be used at the lowest boundary and the values of those
        ! variables can be set to 0 at the top boundary.  Brian.  12/18/05.
        if ( liter ) then 
          rhs(1,1) = xapxbp(1) + 1.0/dt*xapxbp(1)
          rhs(gr%nnzp,1) = 1.0/dt*xapxbp(gr%nnzp)
        else
          rhs(1,1) = xapxbp(1)
          rhs(gr%nnzp,1) = 0.0
        end if

#ifdef STATS
        ! Compute implicit budget terms
        if ( lstats_samp ) then

          do k=1, gr%nnzp

            km1 = max( k-1, 1 )
            kp1 = min( k+1, gr%nnzp )
            if ( ixapxbp_ta > 0 ) then
              zm%x(k,ixapxbp_ta) 
     .          = zm%x(k,ixapxbp_ta)
     .          + term_ta_rhs
     .           ( a1(kp1), a1(k), a1(km1), wp3(kp1), wp3(k),
     .             wp2(kp1), wp2(k), wp2(km1), 
     .             wpxbp(kp1), wpxbp(k), wpxbp(km1),
     .             wpxap(kp1), wpxap(k), wpxap(km1),
     .             gr%dzm(k), beta, eps )

            end if 
            ! rtp2/thlp2 case (1 term)
            if ( ixapxbp_tp > 0 ) then
              zm%x(k,ixapxbp_tp) 
     .          = zm%x(k,ixapxbp_tp)
     .          + term_tp( xam(kp1), xam(k), xbm(kp1), xbm(k), 
     .                     wpxbp(k), wpxap(k), gr%dzm(k) )
              zm%n(k,ixapxbp_tp) = zm%n(k,ixapxbp_tp) + 1
            end if

            ! rtpthlp case (2 terms)
            if ( ixapxbp_tp1 > 0 ) then
              zm%x(k,ixapxbp_tp1) = zm%x(k,ixapxbp_tp1)
     .        - wpxap(k) * gr%dzm(k) * ( xbm(kp1) - xbm(k) )
              zm%n(k,ixapxbp_tp1) = zm%n(k,ixapxbp_tp1) + 1
            end if

            if ( ixapxbp_tp2 > 0 ) then
              zm%x(k,ixapxbp_tp2) = zm%x(k,ixapxbp_tp2)
     .        - wpxbp(k) * gr%dzm(k) * ( xam(kp1) - xam(k) )
              zm%n(k,ixapxbp_tp2) = zm%n(k,ixapxbp_tp2) + 1
            end if

          end do
        end if
#endif /*STATS*/

        return
        end subroutine diag_var_rhs
!-----------------------------------------------------------------------
        pure function term_ta_lhs( a1p1, a1, a1m1, wp3p1, wp3,
     .                             wp2p1, wp2, wp2m1, dzm, beta, eps )
     .  result( lhs )

!       Description:
!       Turbulent advection term (implicit contributions)

!       References:
!-----------------------------------------------------------------------
        implicit none

        ! External
        intrinsic :: max

        ! Constant parameters
        integer, parameter ::
     .  kp1diag = 1,
     .  kdiag   = 2,
     .  km1diag = 3

        ! Input variables
        real, intent(in) ::
     .  a1p1, ! a1(k+1)                         [-]
     .  a1,   ! a(k)                            [-]
     .  a1m1, ! a1(k-1)                         [-]
     .  wp3p1,! w'^3(k+1)                       [m^3/s^3]
     .  wp3,  ! w'^3(k)                         [m^3/s^3]
     .  wp2p1,! w'^2(k+1)                       [m^2/s^2]
     .  wp2,  ! w'^2(k)                         [m^2/s^2]
     .  wp2m1,! w'^2(k+1)                       [m^2/s^2]
     .  dzm,  ! Inverse of the grid spacing     [m]
     .  beta, ! Model parameter                 [-]
     .  eps   ! Model parameter                 [-]

        ! Return Variable
        real, dimension(3) :: lhs

        lhs(kp1diag) 
     .  = + (1.0/3.0) * beta * (1.0/2.0) * dzm *
     .      ( a1p1 + a1 ) * ( wp3p1 / max( wp2p1 + wp2, 2.0*eps ) )

        lhs(kdiag) 
     .  = + (1.0/3.0) * beta * (1.0/2.0) * dzm *
     .      (  ( a1p1 + a1 ) * ( wp3p1 / max( wp2p1 + wp2, 2.0*eps ) )
     .        -( a1 + a1m1 ) * ( wp3 / max( wp2 + wp2m1, 2.0*eps ) ) 
     .      )

        lhs(km1diag) 
     .  = - (1.0/3.0) * beta * (1.0/2.0) * dzm *
     .      ( a1 + a1m1 ) * ( wp3 / max( wp2 + wp2m1, 2.0*eps) )


        return
        end function term_ta_lhs

!-----------------------------------------------------------------------
        pure function term_ta_rhs( a1p1, a1, a1m1, wp3p1, wp3,
     .                             wp2p1, wp2, wp2m1,
     .                             wpxbpp1, wpxbp, wpxbpm1, 
     .                             wpxapp1, wpxap, wpxapm1,
     .                             dzm, beta, eps )
     .  result( rhs )

!       Description:
!       Turbulent advection term (explicit contributions)

!       References:
!-----------------------------------------------------------------------
        implicit none

        ! External
        intrinsic :: max

        ! Input variables
        real, intent(in) ::
     .  a1p1,   ! a1(k+1)                       [-]
     .  a1,     ! a(k)                          [-]
     .  a1m1,   ! a1(k-1)                       [-]
     .  wp3p1,  ! w'^3(k+1)                     [m^3/s^3]
     .  wp3,    ! w'^3(k)                       [m^3/s^3]
     .  wp2p1,  ! w'^2(k+1)                     [m^2/s^2]
     .  wp2,    ! w'^2(k)                       [m^2/s^2]
     .  wp2m1,  ! w'^2(k+1)                     [m^2/s^2]
     .  wpxbpp1,! w'x'(k+1)                     [m/s units]
     .  wpxbp,  ! w'x'(k)                       [m/s units]
     .  wpxbpm1,! w'x'(k-1)                     [m/s units]
     .  wpxapp1,! w'x'(k+1)                     [m/s units]
     .  wpxap,  ! w'x'(k)                       [m/s units]
     .  wpxapm1,! w'x'(k-1)                     [m/s units]
     .  dzm,    ! Inverse of the grid spacing   [m]
     .  beta,   ! Model parameter               [-]
     .  eps     ! Model parameter               [-]

        ! Return Variable
        real :: rhs

        rhs
     .  = - ( 1.0 - (1.0/3.0) * beta ) * (1.0/4.0) * dzm
     .      * (  ( a1p1 + a1 )**2 * ( wpxapp1 + wpxap ) 
     .            * ( wpxbpp1 + wpxbp ) * wp3p1 
     .              / ( max( wp2p1 + wp2, 2.0*eps ) )**2
     .          - ( a1 + a1m1 )**2 * ( wpxap + wpxapm1 ) 
     .            * ( wpxbp + wpxbpm1 ) * wp3
     .              / ( max( wp2 + wp2m1, 2.0*eps ) )**2
     .         )

        return
        end function term_ta_rhs

!-----------------------------------------------------------------------
        pure function term_tp( xamp1, xam, xbmp1, xbm, 
     .                         wpxbp, wpxap, dzm )
     .  result( rhs )

!       Description:
!       Turbulent production

!       References:
!-----------------------------------------------------------------------

        implicit none

        ! Input variables
        real, intent(in) ::
     .  xam,  ! x(k)                            [units vary]
     .  xamp1,! x(k+1)                          [units vary]
     .  xbm,  ! x(k)                            [units vary]
     .  xbmp1,! x(k+1)                          [units vary]
     .  wpxbp,! w'x'(k)                         [m/s units]
     .  wpxap,! w'x'(k)                         [m/s units]
     .  dzm   ! Inverse of the grid spacing     [m]

        ! Return Variable
        real :: rhs

        rhs
     .  = - wpxbp * dzm * ( xamp1 - xam )
     .    - wpxap * dzm * ( xbmp1 - xbm )

        return
        end function term_tp

!-----------------------------------------------------------------------
        pure function term_ma( wmm, dzm ) 
     .  result( lhs )

!       Description:
!       Mean advection term

!       References:
!-----------------------------------------------------------------------
        implicit none

        ! Constant parameters
        integer, parameter ::
     .  kp1diag = 1,
     .  km1diag = 2

        real, intent(in) ::
     .  wmm, ! w wind on m                      [m/s]
     .  dzm  ! Inverse of the grid spacing      [m]

        ! Return Variable
        real, dimension(2) :: lhs

        lhs(kp1diag) = + wmm * 0.5 * dzm
        lhs(km1diag) = - wmm * 0.5 * dzm

        return
        end function term_ma

!-----------------------------------------------------------------------
        pure function term_dp2( nu2, dzm, dzt, dztp1 ) 
     .  result( lhs )

!       Description:
!       Dissipation term

!       References:
!-----------------------------------------------------------------------
        implicit none

        ! Constant parameters
        integer, parameter ::
     .  kp1diag = 1,
     .  kdiag   = 2,
     .  km1diag = 3

        real, intent(in) ::
     .  nu2,  ! Previously determined eddy viscosity     [-]
     .  dzt,  ! Inverse of the grid spacing              [m]
     .  dztp1,! Inverse of the grid spacing              [m]
     .  dzm   ! Inverse of the grid spacing              [m]

        ! Return Variable
        real, dimension(3) :: lhs

        lhs(kp1diag) = - nu2 * dzm * dzt
        lhs(kdiag)   = + nu2 * dzm * ( dztp1 + dzt )
        lhs(km1diag) = - nu2 * dzm * dztp1

        return
        end function term_dp2

!-----------------------------------------------------------------------
        pure function term_dp1( C2, taum ) 
     .  result( lhs )

!       Description:
!       Dissipation term

!       References:
!-----------------------------------------------------------------------
        implicit none

        ! Constant parameters
        integer, parameter ::
     .  kdiag = 1

        real, intent(in) ::
     .  C2, ! Coefficient      [-]
     .  taum ! Tau              [s]

        real, dimension(1) :: lhs

        lhs(kdiag) = + C2 / taum

        return
        end function term_dp1

        end module diagnose_variances
