!-----------------------------------------------------------------------
! $Id: diag_var.F,v 1.52 2008-06-18 18:20:49 faschinj Exp $
#define SCLR_THETA 1
#define SCLR_RT 2
!===============================================================================
      module diagnose_variances

!     Description:
!     Contains the subroutine diag_var and ancillary functions.
!-----------------------------------------------------------------------

        implicit none

        public :: diag_var

        private :: diag_var_lhs, 
     .             diag_var_solve, 
     .             diag_var_uv_rhs,
     .             diag_var_rhs,
     .             term_ta_lhs,
     .             term_ta_rhs,
     .             term_tp,
     .             term_dp1,
     .             term_pr1,
     .             term_pr2
       
        private    ! Set default scope

        contains

!===============================================================================
        subroutine diag_var( taum, wmm, rtm, wprtp, 
     .                       thlm, wpthlp, wpthvp, um, vm,
     .                       wp2, wp3, upwp, vpwp, Scm, Skwm, Kht,
     .                       rtp2, thlp2, rtpthlp,
     .                       up2, vp2,
     .                       liter, dt, err_code,
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
!      /Steady-state solution for the variances/ 
!-----------------------------------------------------------------------

        use constants, only:
     .      sclr_dim, ! Variable(s)
     .      wtol,
     .      rttol,
     .      thltol,
     .      sclrtol,
     .      fstderr
        
        use parameters, 
     .    only: C2, C2b, C2c, C2rt, C2thl, C2rtthl, c_K2, nu2, nu9, ! Variable(s)
     .          c_Ksqd, beta, C4, C14, C5, T0
        
        use grid_class, only:
     .      gr, ! Variable(s)
     .      zm2zt ! Procedure(s)
        
        use stats_prec, only: 
     .      time_prec ! Variable(s)
        
        use lapack_wrap, only: 
     .      tridag_solve ! Procedure(s)
        
        use error_code, only: 
     .      clubb_no_error, ! Variable(s)
     .      lapack_error,   ! Procedure(s)
     .      clubb_at_debug_level
        
#ifdef STATS
        use stats_type, only:
     .      stat_begin_update, stat_end_update ! Procedure(s) 
        
        use stats_hoc, only:
     .      irtp2_bt, ! Variable(s)
     .      irtp2_cl,
     .      ithlp2_bt,
     .      ithlp2_cl,
     .      irtpthlp_bt,
     .      ivp2_bt,
     .      ivp2_cl,
     .      iup2_bt,
     .      iup2_cl,
     .      zm,
     .      lstats_samp
#endif /*STATS*/
       
        implicit none

        ! Input variables
        real, intent(in), dimension(gr%nnzp) :: 
     .  taum,  ! Tau on moment. grid            [s]
     .  wmm,   ! w wind on m                    [m/s]
     .  rtm,   ! Total water mixing ratio       [kg/kg]
     .  wprtp, ! w' r_t'                        [(m kg)/(s kg)]
     .  thlm,  ! Liquid potential temp.         [K]
     .  wpthlp,! w' th_l'                       [(m K)/s]
     .  wpthvp,! w' th_v'                       [(m K)/s]
     .  um,    ! u wind                         [m/s]
     .  vm,    ! v wind                         [m/s]
     .  wp2,   ! w'^2                           [m^2/s^2]
     .  wp3,   ! w'^3                           [m^3/s^3]
     .  upwp,  ! u'w'                           [m^2/s^2]
     .  vpwp,  ! u'w'                           [m^2/s^2]
     .  Scm,   ! Sc on moment. grid             [-]
     .  Skwm,  ! Skw on moment. grid            [-]
     .  Kht    ! Eddy diffusivity on t-lev.     [m^2/s]

        logical, intent(in) :: liter ! Whether variances are prognostic

        real(kind=time_prec), intent(in) :: dt ! Timestep       [s]

        ! Input/Output variables
        ! An attribute of (inout) is also needed to import the
        ! value of the variances at the surface.  Brian.  12/18/05.
        real, intent(inout), dimension(gr%nnzp) :: 
     .  rtp2,   ! r_t'^2                        [(kg/kg)^2]
     .  thlp2,  ! th_l'^2                       [K^2]
     .  rtpthlp,! r_t' th_l'                    [(kg K)/kg]
     .  up2,    ! u'^2                          [m^2/s^2]
     .  vp2     ! v'^2                          [m^2/s^2]

        ! Output variable for singular matrices

        integer, intent(out) :: err_code

        ! Passive scalar input (optional)
        real, optional, intent(in), dimension(gr%nnzp, sclr_dim) :: 
     .  sclrm, wpsclrp

        ! Passive scalar output (optional)
        real, optional, intent(out), dimension(gr%nnzp, sclr_dim) :: 
     .  sclrp2, sclrprtp, sclrpthlp

        ! Local Variables
        real, dimension(gr%nnzp) ::
     .  C2sclr_1d, C2rt_1d, C2thl_1d, C2rtthl_1d, C4_C14_1d

        real, dimension(gr%nnzp) ::
     .  a1 ! a_1 (momentum levels); See eqn. 24 in `Equations for HOC' [-]

        real, dimension(gr%nnzp) ::
     .  a1_zt,     ! a_1 interpolated to thermodynamic levels      [-]
     .  wp2_zt,    ! w'^2 interpolated to thermodynamic levels     [m^2/s^2]
     .  wprtp_zt,  ! w'r_t' interpolated to thermodynamic levels   [(kg/kg) m/s]
     .  wpthlp_zt, ! w'th_l' interpolated to thermodyamnic levels  [K m/s]
     .  upwp_zt,   ! u'w' interpolated to thermodynamic levels     [m^2/s^2]
     .  vpwp_zt,   ! v'w' interpolated to thermodynamic levels     [m^2/s^2]
     .  wpsclrp_zt ! w'sclr' interpolated to thermodynamic levels  [m/s {sclrm units}]

        real, dimension(3,gr%nnzp) :: 
     .  lhs ! Tridiagonal matrix

        real, dimension(gr%nnzp,1) ::
     .  rhs ! RHS vector of Tridiagonal matrix

        real, dimension(gr%nnzp,sclr_dim*3) :: 
     .  sclr_rhs,     ! RHS of scalar tridiagonal system
     .  sclr_solution ! Solution to tridiagonal system

        integer, dimension(5+1) ::
     .  Valid_arr

        ! Eddy Diffusion for Variances and Covariances.
        real, dimension(gr%nnzp) :: Kw2   ! [m^2/s]

        ! Variables used for adding (xapxbp)^2: 3-point average
        ! diffusion coefficient.
        real, dimension(gr%nnzp) ::
     .  rtp2_zt,    ! r_t'^2 interpolated to thermodynamic levels    [kg^2/kg^2]
     .  thlp2_zt,   ! th_l'^2 interpolated to thermodynamic levels   [K^2]
     .  rtpthlp_zt, ! r_t'th_l' interpolated to thermodynamic levels [K kg/kg]
     .  rtp2_zt_sqd_3pt,
     .  thlp2_zt_sqd_3pt,
     .  rtpthlp_zt_sqd_3pt,
     .  Kw2_rtp2,
     .  Kw2_thlp2,
     .  Kw2_rtpthlp

        ! wtol_sqd = the square of the minimum threshold on w,
        !     [wtol_sqd] = m^2 s^{-2}.  Vince Larson 11 Mar 2008.
        real :: wtol_sqd

        logical :: scalar_calc

        ! Loop indices
        integer :: i
        integer :: k, km1, kp1

!-----------------------------------------------------------------------
#ifdef SINGLE_C2_SKW
        ! Use a single value of C2 for all equations.
        C2rt_1d(1:gr%nnzp) 
     .  = C2b + (C2-C2b) *exp( -0.5 * (Skwm(1:gr%nnzp)/C2c)**2 )

        C2thl_1d   = C2rt_1d
        C2rtthl_1d = C2rt_1d

        C2sclr_1d  = C2rt_1d
#else
        ! Use 3 different values of C2 for rtp2, thlp2, rtpthlp.
        C2rt_1d(1:gr%nnzp)    = C2rt
        C2thl_1d(1:gr%nnzp)   = C2thl
        C2rtthl_1d(1:gr%nnzp) = C2rtthl

        C2sclr_1d(1:gr%nnzp)  = C2rt  ! Use rt value for now
#endif
        C4_C14_1d(1:gr%nnzp) = 2.0/3.0 * C4 + ( 1.0/3.0 * C14 )

        ! Are we solving for a passive scalar as well?

        if ( present( sclrm ) .and. present( wpsclrp ) .and.
     .       present( sclrp2 ) .and. present( sclrprtp ) .and.
     .       present( sclrpthlp ) ) then

          scalar_calc = .true.

        else
          scalar_calc = .false.

        end if

        ! Define a_1 (located on momentum levels).
        ! It is a variable that is a function of Sc (where Scm is
        ! located on momentum levels).
        a1(1:gr%nnzp) = 1.0 / ( 1.0 - Scm(1:gr%nnzp) )

        ! wtol_sqd = the square of the minimum threshold on w,
        !     [wtol_sqd] = m^2 s^{-2}.  Vince Larson 11 Mar 2008.
        wtol_sqd = wtol * wtol

        ! Interpolate a_1, w'^2, w'r_t', w'th_l', u'w', and v'w' from 
        ! momentum levels to thermodynamic levels.  These will be used 
        ! for the turbulent advection (ta) terms in each equation.
        a1_zt     = max( zm2zt( a1 ), 0.0 )   ! Positive definite quantity
        wp2_zt    = max( zm2zt( wp2 ), 0.0 )   ! Positive definite quantity
        wprtp_zt  = zm2zt( wprtp )
        wpthlp_zt = zm2zt( wpthlp )
        upwp_zt   = zm2zt( upwp )
        vpwp_zt   = zm2zt( vpwp )

        ! Interpolate r_t'^2, th_l'^2, and r_t'th_l' from momentum 
        ! levels to thermodynamic levels.  These will be used for 
        ! extra diffusion based on a three-point average of (var)^2.
        rtp2_zt    = max( zm2zt( rtp2 ), 0.0 )   ! Positive definite quantity
        thlp2_zt   = max( zm2zt( thlp2 ), 0.0 )   ! Positive definite quantity
        rtpthlp_zt = zm2zt( rtpthlp )

#ifdef STATS
        if ( lstats_samp ) then

          call stat_begin_update( irtp2_bt, real(rtp2 / dt), zm )

          call stat_begin_update( ithlp2_bt, real(thlp2 / dt), zm )

          call stat_begin_update( irtpthlp_bt, real(rtpthlp / dt), zm )

          call stat_begin_update( ivp2_bt, real(vp2 / dt), zm )

          call stat_begin_update( iup2_bt, real(up2 / dt), zm )

        end if
#endif

        ! Initialize tridiagonal solutions to valid

        Valid_arr(:) = clubb_no_error

        ! (xapxbp)^2: 3-point average diffusion coefficient.
        do k = 1, gr%nnzp, 1

           km1 = max( k-1, 1 )
           kp1 = min( k+1, gr%nnzp )

           ! Compute the square of rtp2_zt, averaged over 3 points.  26 Jan 2008
           rtp2_zt_sqd_3pt(k) = ( rtp2_zt(km1)**2 + rtp2_zt(k)**2 
     .                           + rtp2_zt(kp1)**2 ) / 3.0
           ! Account for units (kg/kg)**4  Vince Larson 29 Jan 2008
           rtp2_zt_sqd_3pt(k) = 1e12 * rtp2_zt_sqd_3pt(k)  

           ! Compute the square of thlp2_zt, averaged over 3 points.  26 Jan 2008
           thlp2_zt_sqd_3pt(k) = ( thlp2_zt(km1)**2 + thlp2_zt(k)**2 
     .                            + thlp2_zt(kp1)**2 ) / 3.0

           ! Compute the square of rtpthlp_zt, averaged over 3 points.  26 Jan 2008
           rtpthlp_zt_sqd_3pt(k) 
     .          = ( rtpthlp_zt(km1)**2 + rtpthlp_zt(k)**2 
     .             + rtpthlp_zt(kp1)**2 ) / 3.0
           ! Account for units (kg/kg)**2 Vince Larson 29 Jan 2008
           rtpthlp_zt_sqd_3pt(k) = 1e6 * rtpthlp_zt_sqd_3pt(k)

        enddo

        ! Define the Coefficent of Eddy Diffusivity for the 
        ! Variances and Covariances.
        do k = 1, gr%nnzp, 1
           ! Kw2 is used for all variances and covariances in subroutine
           ! diag_var.  The variances and covariances are located on 
           ! momentum levels.  Kw2 is located on thermodynamic levels.
           ! Kw2 = c_K2 * Kht
           Kw2(k) = c_K2 * Kht(k)
           ! Kw2_rtp2 must have units of m^2/s.  Since rtp2_zt_sqd_3pt has 
           ! units of kg^2/kg^2, c_Ksqd is given units of m^2/[ s (kg^2/kg^2) ]
           ! in this case.
           Kw2_rtp2(k)    = Kw2(k) + c_Ksqd * rtp2_zt_sqd_3pt(k) ! Vince Larson increased by c_Ksqd, 29Jan2008
           ! Kw2_thlp2 must have units of m^2/s.  Since thlp2_zt_sqd_3pt has 
           ! units of K^2, c_Ksqd is given units of m^2/[ s K^2 ] in this case.
           Kw2_thlp2(k)   = Kw2(k) + c_Ksqd * thlp2_zt_sqd_3pt(k) ! Vince Larson increased by c_Ksqd, 29Jan2008
           ! Kw2_rtpthlp must have units of m^2/s.  Since rtpthlp_zt_sqd_3pt has
           ! units of K (kg/kg), c_Ksqd is given units of m^2/[ s K (kg/kg) ] 
           ! in this case.
           Kw2_rtpthlp(k) = Kw2(k) + c_Ksqd * rtpthlp_zt_sqd_3pt(k) ! Vince Larson increased by c_Ksqd, 29Jan2008
        enddo


        !!!!!***** r_t'^2 *****!!!!!
        ! Implicit contributions to term rtp2

        ! Implicit contributions to term rtp2
        call diag_var_lhs( dt, liter, a1, a1_zt, wp2_zt, 
     .                     wp3, taum, wmm, Kw2_rtp2, C2rt_1d, 
     .                     nu2, beta, wtol_sqd, lhs )

        ! Explicit contributions to rtp2
        call diag_var_rhs( "rtp2", dt, liter, a1, a1_zt, 
     .                     wp2_zt, wp3, wprtp, wprtp_zt,
     .                     wprtp, wprtp_zt, rtm, rtm, rtp2,
     .                     beta, wtol_sqd, rhs )
        
        ! Solve the tridiagonal system
        call diag_var_solve( "rtp2", 1, rhs, 
     .                       lhs, rtp2, Valid_arr(1) )


        !!!!!***** th_l'^2 *****!!!!!
        ! Implicit contributions to term thlp2

        ! Implicit contributions to term thlp2
        call diag_var_lhs( dt, liter, a1, a1_zt, wp2_zt, 
     .                     wp3, taum, wmm, Kw2_thlp2, C2thl_1d, 
     .                     nu2, beta, wtol_sqd, lhs )

        ! Explicit contributions to thlp2
        call diag_var_rhs( "thlp2", dt, liter, a1, a1_zt, 
     .                     wp2_zt, wp3, wpthlp, wpthlp_zt,
     .                     wpthlp, wpthlp_zt, thlm, thlm, thlp2,
     .                     beta, wtol_sqd, rhs )

        ! Solve the tridiagonal system
        call diag_var_solve( "thlp2", 1, rhs, 
     .                       lhs, thlp2, Valid_arr(2) )


        !!!!!***** r_t'th_l' *****!!!!!
        ! Implicit contributions to term rtpthlp

        call diag_var_lhs( dt, liter, a1, a1_zt, wp2_zt, 
     .                     wp3, taum, wmm, Kw2_rtpthlp, C2rtthl_1d, 
     .                     nu2, beta, wtol_sqd, lhs )

        ! Explicit contributions to rtpthlp
        call diag_var_rhs( "rtpthlp", dt, liter, a1, a1_zt, 
     .                     wp2_zt, wp3, wprtp, wprtp_zt,
     .                     wpthlp, wpthlp_zt, rtm, thlm, rtpthlp,
     .                     beta, wtol_sqd, rhs )

        ! Solve the tridiagonal system
        call diag_var_solve( "rtpthlp", 1, rhs, 
     .                       lhs, rtpthlp, Valid_arr(3) )


        !!!!!***** u'^2 *****!!!!!
        ! Implicit contributions to term up2

        call diag_var_lhs( dt, liter, a1, a1_zt, wp2_zt, 
     .                     wp3, taum, wmm, Kw2, C4_C14_1d, 
     .                     nu9, beta, wtol_sqd, lhs )

        ! Explicit contributions to up2
        call diag_var_uv_rhs( "up2", dt, liter, a1, a1_zt,
     .                        wp2, wp2_zt, wp3, wpthvp, taum, 
     .                        um, vm, upwp, upwp_zt, vpwp,
     .                        vpwp_zt, up2, vp2, C4, C5, C14,
     .                        T0, beta, wtol_sqd, rhs )

        ! Solve the tridiagonal system
        call diag_var_solve( "up2", 1, rhs,
     .                       lhs, up2, Valid_arr(4) )


        !!!!!***** v'^2 *****!!!!!
        ! Implicit contributions to term vp2

        call diag_var_lhs( dt, liter, a1, a1_zt, wp2_zt, 
     .                     wp3, taum, wmm, Kw2, C4_C14_1d, 
     .                     nu9, beta, wtol_sqd, lhs )

        ! Explicit contributions to vp2
        call diag_var_uv_rhs( "vp2", dt, liter, a1, a1_zt,
     .                        wp2, wp2_zt, wp3, wpthvp, taum, 
     .                        vm, um, vpwp, vpwp_zt, upwp,
     .                        upwp_zt, vp2, up2, C4, C5, C14,
     .                        T0, beta, wtol_sqd, rhs )

        ! Solve the tridiagonal system
        call diag_var_solve( "vp2", 1, rhs, 
     .                       lhs, vp2, Valid_arr(5) )

#ifdef STATS
        if( lstats_samp ) then
          call stat_begin_update( irtp2_cl, real( rtp2 / dt), zm )      

          call stat_begin_update( ithlp2_cl, real( thlp2 / dt), zm )

          call stat_begin_update( iup2_cl, real( up2 / dt ), zm )

          call stat_begin_update( ivp2_cl, real( vp2 / dt ), zm )
        end if
#endif

        ! Clip small and negative values
        ! The surface variances are set elsewhere.  They do not need to be
        ! clipped at the lowest level.  Brian Griffin.  12/18/05.
        do k = 2, gr%nnzp, 1

          if ( wp2(k) >= wtol*wtol ) then
            rtp2(k) 
     .      = max( rtp2(k), 
! Vince Larson commented out because may cause unstable feedback w/ wprtp eqn..
!       The correlation of w and rt may now exceed 1, which would be unphysical.  7 Jul 2007 
!     .             wprtp(k)*wprtp(k)/((1.0-Scm(k))*wp2(k)) + eps,
! End Vince Larson's commenting
     .             rttol*rttol )
            thlp2(k) 
     .      = max( thlp2(k), 
! Vince Larson commented out because may cause unstable feedback w/ wpthlp eqn.. 
!    The correlation of w and thl may now exceed 1, which would be unphysical.  7 Jul 2007 
!     .             wpthlp(k)*wpthlp(k)/((1.0-Scm(k))*wp2(k)) + eps,
! End Vince Larson's commenting
     .             thltol*thltol )
          else
            thlp2(k) = max( thlp2(k), 0.0 )
            rtp2(k)  = max( rtp2(k), 0.0 )
          end if

        end do ! k=2..gr%nnzp

        ! Clipping for up2/vp2
        ! This appears to be unneeded for now.  If these do need to be
        ! clipped in the future, this could be uncommented and budget
        ! could be added. -dschanen 24 April 2007
       where ( vp2 < 0.0 ) vp2 = 0.0
       where ( up2 < 0.0 ) up2 = 0.0

#ifdef STATS
        if ( lstats_samp ) then
           call stat_end_update( irtp2_cl, real( rtp2 / dt ), zm )
           
           call stat_end_update( ithlp2_cl, real( thlp2 / dt ), zm )

           call stat_end_update( iup2_cl, real( up2 / dt ), zm )
           
           call stat_end_update( ivp2_cl, real( vp2 / dt ), zm )
        endif
#endif

#ifdef STATS
        if ( lstats_samp ) then
          call stat_end_update( irtp2_bt, real( rtp2 / dt), zm )

          call stat_end_update( ithlp2_bt, real( thlp2 / dt), zm )

          call stat_end_update( irtpthlp_bt, real( rtpthlp / dt), zm )
          
          call stat_end_update( iup2_bt, real( up2 / dt), zm )

          call stat_end_update( ivp2_bt, real( vp2 / dt), zm )
        end if
#endif /*STATS*/

        if ( scalar_calc ) then

          ! Implicit contributions to passive scalars

          !!!!!***** sclr'^2, sclr'r_t', sclr'th_l' *****!!!!!

          call diag_var_lhs( dt, liter, a1, a1_zt, wp2_zt, 
     .                       wp3, taum, wmm, Kw2, C2sclr_1d, 
     .                       nu2, beta, wtol_sqd, lhs )


          ! Explicit contributions to passive scalars


          !!!!!***** sclr'^2 *****!!!!!

          do i = 1, sclr_dim, 1

            ! Interpolate w'sclr' from momentum levels to thermodynamic 
            ! levels.  These will be used for the turbulent advection (ta) 
            ! terms in each equation.
            wpsclrp_zt = zm2zt( wpsclrp(:,i) )

            call diag_var_rhs( "sclrp2", dt, liter, a1, a1_zt, 
     .                         wp2_zt, wp3, wpsclrp(:,i), 
     .                         wpsclrp_zt, wpsclrp(:,i), wpsclrp_zt, 
     .                         sclrm(:,i), sclrm(:,i), sclrp2(:,i),
     .                         beta, wtol_sqd, sclr_rhs(:,i) )

          enddo


          !!!!!***** sclr'r_t' *****!!!!!

          do i = 1, sclr_dim, 1

            ! Interpolate w'sclr' from momentum levels to thermodynamic 
            ! levels.  These will be used for the turbulent advection (ta) 
            ! terms in each equation.
            wpsclrp_zt = zm2zt( wpsclrp(:,i) )

            call diag_var_rhs( "sclrprtp", dt, liter, a1, a1_zt, 
     .                         wp2_zt, wp3, wpsclrp(:,i), 
     .                         wpsclrp_zt, wprtp, wprtp_zt, sclrm(:,i), 
     .                         rtm, sclrprtp(:,i), beta, wtol_sqd, 
     .                         sclr_rhs(:,i+sclr_dim) )

          enddo


          !!!!!***** sclr'th_l' *****!!!!!

          do i = 1, sclr_dim, 1

            ! Interpolate w'sclr' from momentum levels to thermodynamic 
            ! levels.  These will be used for the turbulent advection (ta) 
            ! terms in each equation.
            wpsclrp_zt = zm2zt( wpsclrp(:,i) )

            call diag_var_rhs( "sclrpthlp", dt, liter, a1, a1_zt, 
     .                         wp2_zt, wp3, wpsclrp(:,i), 
     .                         wpsclrp_zt, wpthlp, wpthlp_zt, 
     .                         sclrm(:,i), thlm, sclrpthlp(:,i), beta, 
     .                         wtol_sqd, sclr_rhs(:,i+2*sclr_dim) )

          enddo

          ! Solve the tridiagonal system

          call diag_var_solve
     .         ( "scalars", 3*sclr_dim, sclr_rhs, 
     .           lhs, sclr_solution, Valid_arr(6) )

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
! Vince Larson commented out because it may cause unstable feedback with wpsclrp eqn..
!     The correlation of w and sclr may now exceed 1, which would be unphysical.  7 Jul 2007
!     .              wpsclrp(k,i)*wpsclrp(k,i)/((1.0-Scm(k))*wp2(k))+eps,
! End Vince Larson's commenting
     .              sclrtol(i)*sclrtol(i) )
              else
                sclrp2(k,i) = max( sclrp2(k,i), 0.0 )
              end if ! wp2(k) >= wtol**2
            end do ! k
          end do ! i

        end if ! scalar_calc

        ! Check for singular matrices
        do i = 1, 5+1
                if(Valid_arr(i) > err_code ) err_code = Valid_arr(i)
        enddo        
                
        if ( lapack_error( err_code ) .and. 
     .       clubb_at_debug_level( 1 ) ) then
                
           write(fstderr,*) "Error in diag_var"
           
           write(fstderr,*) "Intent(in)"
           
           write(fstderr,*) "taum = ", taum
           write(fstderr,*) "wmm = ", wmm
           write(fstderr,*) "rtm = ", rtm
           write(fstderr,*) "wprtp = ", wprtp
           write(fstderr,*) "thlm = ", thlm
           write(fstderr,*) "wpthlp = ", wpthlp
           write(fstderr,*) "wpthvp = ", wpthvp
           write(fstderr,*) "um = ", um
           write(fstderr,*) "vm = ", vm
           write(fstderr,*) "wp2 = ", wp2
           write(fstderr,*) "wp3 = ", wp3
           write(fstderr,*) "upwp = ", upwp
           write(fstderr,*) "vpwp = ", vpwp
           write(fstderr,*) "Scm = ", Scm
           write(fstderr,*) "Skwm = ", Skwm
           write(fstderr,*) "Kht = ", Kht
           
           write(fstderr,*) "Intent(In/Out)"
           
           write(fstderr,*) "rtp2 = ", rtp2
           write(fstderr,*) "thlp2 = ", thlp2
           write(fstderr,*) "rtpthlp = ", rtpthlp
           write(fstderr,*) "up2 = ", up2
           write(fstderr,*) "vp2 = ", vp2
           
           if( present( sclrp2 ) ) then 
              write(fstderr,*) "sclrp2 = ", sclrp2
           endif
           
           if( present( sclrprtp ) ) then
             write(fstderr,*) "sclrprtp = ", sclrprtp
           endif
           
           if( present( sclrpthlp ) ) then
             write(fstderr,*) "sclrthlp = ", sclrpthlp
           endif
           
           write(fstderr,*) "Intent(out)"
           
           if( present( sclrm ) ) then
              write(fstderr,*) "sclrm = ", sclrm
           endif
           
           if( present( wpsclrp ) ) then
              write(fstderr,*) "wpsclrp = ", wpsclrp
           endif
           
        endif 
        
        return
        
        end subroutine diag_var

!===============================================================================
        subroutine diag_var_lhs( dt, liter, a1, a1_zt, wp2_zt, 
     .                           wp3, taum, wmm, Kw2, Cn, 
     .                           nu, beta, wtol_sqd, lhs )
        
!       Description:
!       Compute LHS tridiagonal matrix for a variance or coveriance term

!       References:
!       None
!-----------------------------------------------------------------------

        use grid_class, only:
     .      gr ! Variable(s)

        use stats_prec, only: 
     .      time_prec ! Variable(s)

        use diffusion, only: 
     .      diffusion_zm_lhs ! Procedure(s)

        use mean_adv, only: 
     .      term_ma_zm_lhs ! Procedure(s)

#ifdef STATS
        use stats_hoc, only:
     .      zmscr01,
     .      zmscr02,
     .      zmscr03,
     .      zmscr04,
     .      zmscr05,
     .      zmscr06,
     .      zmscr07,
     .      zmscr08,
     .      zmscr09,
     .      zmscr10,
     .      lstats_samp,
     .      irtp2_ma,
     .      irtp2_ta,
     .      irtp2_dp1,
     .      irtp2_dp2,
     .      ithlp2_ma,
     .      ithlp2_ta,
     .      ithlp2_dp1,
     .      ithlp2_dp2,
     .      irtpthlp_ma,
     .      irtpthlp_ta,
     .      irtpthlp_dp1,
     .      irtpthlp_dp2,
     .      iup2_ma,
     .      iup2_ta,
     .      iup2_dp2,
     .      ivp2_ma,
     .      ivp2_ta,
     .      ivp2_dp2

#endif

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  kp1_mdiag = 1,   ! Momentum superdiagonal index.
     .  k_mdiag   = 2,   ! Momentum main diagonal index.
     .  km1_mdiag = 3    ! Momentum subdiagonal index.

        real(kind=time_prec), intent(in) ::
     .  dt     ! Timestep length                             [s]

        logical, intent(in) ::
     .  liter  ! Whether the variances are prognostic

        ! Input Variables
        real, dimension(gr%nnzp), intent(in) ::
     .  a1,    ! Scm-related term a_1 (momentum levels)      [-]
     .  a1_zt, ! a_1 interpolated to thermodynamic levels    [-]
     .  wp2_zt,! w'^2 interpolated to thermodynamic levels   [m^2/s^2]
     .  wp3,   ! w'^3 (thermodynamic levels)                 [m^3/s^3]
     .  taum,  ! Time-scale tau on momentum levels           [s]
     .  wmm,   ! w wind component on momentum levels         [m/s]
     .  Kw2,   ! Coefficient of eddy diffusivity (all vars.) [m^2/s]
     .  Cn     ! Coefficient C_n                             [-]

        real, intent(in) ::
     .  nu,      ! Background constant coef. of eddy diff.   [-]
     .  beta,    ! Constant model parameter beta             [-]
     .  wtol_sqd ! w wind component tolerance squared        [m^2/s^2]

        ! Output Variables
        real, dimension(3,gr%nnzp), intent(out) ::
     .  lhs ! Implicit contributions to the term

        ! Local Variables

        ! Array indices
        integer :: k, kp1, km1

#ifdef STATS
        real, dimension(3) ::
     .  tmp
#endif /*STATS*/


        ! Setup LHS of the tridiagonal system
        do k = 2, gr%nnzp-1, 1

          km1 = max( k-1, 1 )
          kp1 = min( k+1, gr%nnzp )

          ! LHS eddy diffusion term: dissipation term 2 (dp2).
          lhs(kp1_mdiag:km1_mdiag,k)
     .    = diffusion_zm_lhs( Kw2(k), Kw2(kp1), nu,
     .                        gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )

          ! LHS dissipation term 1 (dp1)
          ! (and pressure term 1 (pr1) for u'^2 and v'^2).
          lhs(k_mdiag,k)
     .    = lhs(k_mdiag,k) + term_dp1( Cn(k), taum(k) )

          ! LHS time tendency.
          if ( liter ) then
            lhs(k_mdiag,k) = lhs(k_mdiag,k) + ( 1.0 / dt )
          endif

          ! LHS mean advection (ma) term.
          lhs(kp1_mdiag:km1_mdiag,k)
     .    = lhs(kp1_mdiag:km1_mdiag,k)
     .    + term_ma_zm_lhs( wmm(k), gr%dzm(k), k )

          ! LHS turbulent advection (ta) term.
          lhs(kp1_mdiag:km1_mdiag,k)
     .    = lhs(kp1_mdiag:km1_mdiag,k)
     .    + term_ta_lhs( a1_zt(kp1), a1(k), a1_zt(k), 
     .                   wp3(kp1), wp3(k), wp2_zt(kp1), wp2_zt(k), 
     .                   gr%dzm(k), beta, wtol_sqd, k )

#ifdef STATS
         if ( lstats_samp ) then

         ! Statistics: implicit contributions for rtp2, thlp2, 
         !             rtpthlp, up2, or vp2.

         if ( irtp2_dp1 + ithlp2_dp1 + irtpthlp_dp1  > 0 ) then
           ! Note:  The statistical implicit contribution to term dp1
           !        (as well as to term pr1) for up2 and vp2 is recorded
           !        in diag_var_uv_rhs because up2 and vp2 use a special
           !        dp1/pr1 combined term.
           tmp(1) = term_dp1( Cn(k), taum(k) )
           zmscr01(k) = -tmp(1)
         endif
         
         if ( irtp2_dp2 + ithlp2_dp2 + irtpthlp_dp2 + 
     .        iup2_dp2 + ivp2_dp2 > 0 ) then
           tmp(1:3)
     .     = diffusion_zm_lhs( Kw2(k), Kw2(kp1), nu,
     .                         gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )
           zmscr02(k) = -tmp(3)
           zmscr03(k) = -tmp(2)
           zmscr04(k) = -tmp(1)
         endif

         if ( irtp2_ta + ithlp2_ta + irtpthlp_ta +
     .        iup2_ta + ivp2_ta > 0 ) then
           tmp(1:3)
     .     = term_ta_lhs( a1_zt(kp1), a1(k), a1_zt(k), 
     .                    wp3(kp1), wp3(k), wp2_zt(kp1), wp2_zt(k), 
     .                    gr%dzm(k), beta, wtol_sqd, k )
           zmscr05(k) = -tmp(3)
           zmscr06(k) = -tmp(2)
           zmscr07(k) = -tmp(1)
         endif

         if ( irtp2_ma + ithlp2_ma + irtpthlp_ma +
     .        iup2_ma + ivp2_ma > 0 ) then
           tmp(1:3)
     .     = term_ma_zm_lhs( wmm(k), gr%dzm(k), k )
           zmscr08(k) = -tmp(3)
           zmscr09(k) = -tmp(2)
           zmscr10(k) = -tmp(1)
         endif

         endif ! lstats_samp
#endif /*STATS*/

        enddo ! k=2..gr%nnzp-1


        ! Boundary Conditions
        ! These are set so that the sfc_var value of the variances and 
        ! covariances can be used at the lowest boundary and the values of those
        ! variables can be set to 0 at the top boundary.
        ! Fixed-point boundary conditions are used for both the variances and 
        ! the covariances.
        lhs(:,1) = 0.0
        lhs(:,gr%nnzp) = 0.0

        lhs(k_mdiag,1) = 1.0
        lhs(k_mdiag,gr%nnzp) = 1.0

        ! This boundary condition was changed by dschanen on 24 April 2007
        ! When we run prognostically we want to preserve the surface value.
!       if ( liter ) then
!         lhs(k_mdiag,1) = 1.0/dt
!         lhs(k_mdiag,1) = 1.0/dt
!       end if


        return
        end subroutine diag_var_lhs

!===============================================================================
        subroutine diag_var_solve( solve_type, nrhs, rhs, 
     .                             lhs, xapxbp, err_code )

!-----------------------------------------------------------------------

        use lapack_wrap, only: 
     .      tridag_solve, ! Variable(s)
     .      tridag_solvex
!    .    , band_solve
        
        use grid_class, only:
     .      gr ! Variable(s)
        
#ifdef STATS
        use stats_type, only:
     .      stat_update_var_pt ! Procedure(s)
        
        use stats_hoc, only:
     .      zm, ! Variable(s) 
     .      sfc,
     .      irtp2_cn,
     .      irtp2_dp1,
     .      irtp2_dp2,
     .      irtp2_ta,
     .      irtp2_ma,
     .      ithlp2_cn,
     .      ithlp2_dp1,
     .      ithlp2_dp2,
     .      ithlp2_ta,
     .      ithlp2_ma,
     .      irtpthlp_cn,
     .      irtpthlp_dp1,
     .      irtpthlp_dp2,
     .      irtpthlp_ta,
     .      irtpthlp_ma,
     .      iup2_cn,
     .      iup2_dp1,
     .      iup2_dp2,
     .      iup2_ta,
     .      iup2_ma,
     .      iup2_pr1,
     .      ivp2_cn,
     .      ivp2_dp1,
     .      ivp2_dp2,
     .      ivp2_ta,
     .      ivp2_ma,
     .      ivp2_pr1,
     .      lstats_samp,
     .      zmscr01,
     .      zmscr02,
     .      zmscr03,
     .      zmscr04,
     .      zmscr05,
     .      zmscr06,
     .      zmscr07,
     .      zmscr08,
     .      zmscr09,
     .      zmscr10,
     .      zmscr11

#endif

        ! Constant parameters
        integer, parameter ::
     .  kp1_mdiag = 1,   ! Momentum superdiagonal index.
     .  k_mdiag   = 2,   ! Momentum main diagonal index.
     .  km1_mdiag = 3    ! Momentum subdiagonal index.

        ! Input variables
        integer, intent(in) :: nrhs  ! Number of right hand side vectors

        character(len=*), intent(in) :: 
     .  solve_type ! Variable(s) description

        ! Input/Ouput variables
        real, dimension(3,gr%nnzp), intent(inout) ::
     .  lhs  ! Implicit contributions to x variance/covariance term

        real, dimension(gr%nnzp,nrhs), intent(inout) ::
     .  rhs  ! Explicit contributions to x variance/covariance term

        ! Output variables
        real, dimension(gr%nnzp,nrhs), intent(inout) :: 
     .  xapxbp ! Computed value of the variable at <t+1> [units vary]

        integer, intent(out) ::
     .  err_code ! Returns an error code in the event of a singular matrix

        ! Local variables
        real :: rcond  ! Est. of the reciprocal of the condition #

        ! Array indices
        integer :: k, kp1, km1

#ifdef STATS
        integer ::
     .  ixapxbp_dp1,
     .  ixapxbp_dp2,
     .  ixapxbp_ta,
     .  ixapxbp_ma,
     .  ixapxbp_pr1,
     .  ixapxbp_cn

        select case ( trim( solve_type ) )
        case ( "rtp2" )
          ixapxbp_cn  = irtp2_cn
          ixapxbp_dp1 = irtp2_dp1
          ixapxbp_dp2 = irtp2_dp2
          ixapxbp_ta  = irtp2_ta
          ixapxbp_ma  = irtp2_ma
          ixapxbp_pr1 = 0
        case ( "thlp2" )
          ixapxbp_cn  = ithlp2_cn
          ixapxbp_dp1 = ithlp2_dp1
          ixapxbp_dp2 = ithlp2_dp2
          ixapxbp_ta  = ithlp2_ta
          ixapxbp_ma  = ithlp2_ma
          ixapxbp_pr1 = 0
        case ( "rtpthlp" )
          ixapxbp_cn  = irtpthlp_cn
          ixapxbp_dp1 = irtpthlp_dp1
          ixapxbp_dp2 = irtpthlp_dp2
          ixapxbp_ta  = irtpthlp_ta
          ixapxbp_ma  = irtpthlp_ma
          ixapxbp_pr1 = 0
        case ( "up2" )
          ixapxbp_cn  = iup2_cn
          ixapxbp_dp1 = iup2_dp1
          ixapxbp_dp2 = iup2_dp2
          ixapxbp_ta  = iup2_ta
          ixapxbp_ma  = iup2_ma
          ixapxbp_pr1 = iup2_pr1
        case ( "vp2" )
          ixapxbp_cn  = ivp2_cn
          ixapxbp_dp1 = ivp2_dp1
          ixapxbp_dp2 = ivp2_dp2
          ixapxbp_ta  = ivp2_ta
          ixapxbp_ma  = ivp2_ma
          ixapxbp_pr1 = ivp2_pr1
        case default ! No budgets for passive scalars
          ixapxbp_cn  = 0
          ixapxbp_dp1 = 0
          ixapxbp_dp2 = 0
          ixapxbp_ta  = 0
          ixapxbp_ma  = 0
          ixapxbp_pr1 = 0
        end select
#endif

#ifdef STATS 
       if ( lstats_samp .and. ixapxbp_cn > 0 ) then
          call tridag_solvex
     .         ( solve_type, gr%nnzp, nrhs, lhs(kp1_mdiag,:), 
     .           lhs(k_mdiag,:), lhs(km1_mdiag,:), rhs(:,1:nrhs), 
     .           xapxbp(:,1:nrhs), rcond, err_code )

          ! Est. of the condition number of the variance LHS matrix 
         call stat_update_var_pt( ixapxbp_cn, 1, 1.0 / rcond, sfc )

        else 
          call tridag_solve
     .         ( solve_type, gr%nnzp, nrhs, lhs(kp1_mdiag,:), 
     .           lhs(k_mdiag,:), lhs(km1_mdiag,:), rhs(:,1:nrhs), 
     .           xapxbp(:,1:nrhs), err_code )
        end if
#else
        ! Solve system using tridag_solve. This uses LAPACK sgtsv,
        ! which relies on Gaussian elimination to decompose the matrix.
        call tridag_solve
     .       ( solve_type, gr%nnzp, nrhs, lhs(kp1_mdiag,:), 
     .         lhs(k_mdiag,:), lhs(km1_mdiag,:), rhs(:,1:nrhs), 
     .         xapxbp(:,1:nrhs), err_code )

        ! Alternative: Use LU decomposition instead.
!       call band_solve
!    .       ( solve_type, 1, 1, gr%nnzp, nrhs, 
!    .         lhs, rhs(:,1:nrhs), xapxbp(:,1:nrhs), isValid )
#endif

#ifdef STATS
        ! Compute implicit budget terms
        if ( lstats_samp ) then

          do k=2, gr%nnzp-1

            km1 = max( k-1, 1 )
            kp1 = min( k+1, gr%nnzp )

            call stat_update_var_pt( ixapxbp_dp1, k,
     .          zmscr01(k) * xapxbp(k,1), zm )

            call stat_update_var_pt( ixapxbp_dp2, k,
     .         zmscr02(k) * xapxbp(km1,1)
     .        + zmscr03(k) * xapxbp(k,1)
     .        + zmscr04(k) * xapxbp(kp1,1), zm )
           
            call stat_update_var_pt( ixapxbp_ta, k,
     .          zmscr05(k) * xapxbp(km1,1)
     .        + zmscr06(k) * xapxbp(k,1)
     .        + zmscr07(k) * xapxbp(kp1,1), zm )

            call stat_update_var_pt( ixapxbp_ma, k,
     .          zmscr08(k) * xapxbp(km1,1)
     .        + zmscr09(k) * xapxbp(k,1)
     .        + zmscr10(k) * xapxbp(kp1,1), zm )

            call stat_update_var_pt( ixapxbp_pr1, k,
     .          zmscr11(k) * xapxbp(k,1), zm )

          end do
        end if
#endif /*STATS*/

        return
        end subroutine diag_var_solve

!===============================================================================
        subroutine diag_var_uv_rhs( solve_type, dt, liter, a1, a1_zt,
     .                              wp2, wp2_zt, wp3, wpthvp, taum, 
     .                              xam, xbm, wpxap, wpxap_zt, wpxbp,
     .                              wpxbp_zt, xap2, xbp2, C4, C5, C14,
     .                              T0, beta, wtol_sqd, rhs )

!       Description:
!       Explicit contributions to u'^2 or v'^2
!-----------------------------------------------------------------------

        use grid_class, only:
     .      gr ! Variable(s)
        
        use constants, only: 
     .      grav ! Variable(s)
        
        use stats_prec, only: 
     .      time_prec ! Variable(s)
        
#ifdef STATS
        use stats_type, only:
     .      stat_modify_pt, stat_update_var_pt ! Procedure(s)

        use stats_hoc, only:
     .      ivp2_ta, ! Variable(s)
     .      ivp2_tp,
     .      ivp2_dp1,
     .      ivp2_pr1,
     .      ivp2_pr2,
     .      iup2_ta,
     .      iup2_tp,
     .      iup2_dp1,
     .      iup2_pr1,
     .      iup2_pr2,
     .      zm,
     .      zmscr01,
     .      zmscr11,
     .      lstats_samp
#endif

        implicit none

        ! Input Variables
        character(len=*), intent(in) :: solve_type

        real(kind=time_prec), intent(in) ::
     .  dt        ! Timestep                                    [s]

        logical, intent(in) ::
     .  liter  ! Whether x is prognostic (T/F)

        real, dimension(gr%nnzp), intent(in) ::
     .  a1,       ! Scm-related term a_1 (momentum levels)      [-]
     .  a1_zt,    ! a_1 interpolated to thermodynamic levels    [-]
     .  wp2,      ! w'^2 (momentum levels)                      [m^2/s^2]
     .  wp2_zt,   ! w'^2 interpolated to thermodynamic levels   [m^2/s^2]
     .  wp3,      ! w'^3 (thermodynamic levels)                 [m^3/s^3]
     .  wpthvp,   ! w'th_v' (momentum levels)                   [K m/s]
     .  taum,     ! Time-scale tau on momentum levels           [s]
     .  xam,      ! x_am (thermodynamic levels)                 [m/s]
     .  xbm,      ! x_bm (thermodynamic levels)                 [m/s]
     .  wpxap,    ! w'x_a' (momentum levels)                    [m^2/s^2]
     .  wpxap_zt, ! w'x_a' interpolated to thermodynamic levels [m^2/s^2]
     .  wpxbp,    ! w'x_b' (momentum levels)                    [m^2/s^2]
     .  wpxbp_zt, ! w'x_b' interpolated to thermodynamic levels [m^2/s^2]
     .  xap2,     ! x_a'^2 (momentum levels)                    [m^2/s^2]
     .  xbp2      ! x_b'^2 (momentum levels)                    [m^2/s^2]

        real, intent(in) ::
     .  C4,      ! Model parameter C_4                          [-]
     .  C5,      ! Model parameter C_5                          [-]
     .  C14,     ! Model parameter C_14                         [-]
     .  T0,      ! Reference temperature                        [K]
     .  beta,    ! Model parameter beta                         [-]
     .  wtol_sqd ! w wind component tolerance squared           [m^2/s^2]

        real, dimension(gr%nnzp,1), intent(out) ::
     .  rhs    ! Explicit contributions to x variance/covariance terms

        ! Local Variables

        ! Array indices
        integer :: k, kp1, km1

#ifdef STATS
        real :: tmp
#endif /*STATS*/

#ifdef STATS
        integer ::
     .  ixapxbp_ta,
     .  ixapxbp_tp,
     .  ixapxbp_dp1,
     .  ixapxbp_pr1,
     .  ixapxbp_pr2


        select case ( trim( solve_type ) )
        case ( "vp2" )
          ixapxbp_ta  = ivp2_ta
          ixapxbp_tp  = ivp2_tp
          ixapxbp_dp1 = ivp2_dp1
          ixapxbp_pr1 = ivp2_pr1
          ixapxbp_pr2 = ivp2_pr2
        case ( "up2" )
          ixapxbp_ta  = iup2_ta
          ixapxbp_tp  = iup2_tp
          ixapxbp_dp1 = iup2_dp1
          ixapxbp_pr1 = iup2_pr1
          ixapxbp_pr2 = iup2_pr2
        case default ! No budgets for passive scalars
          ixapxbp_ta  = 0
          ixapxbp_tp  = 0
          ixapxbp_dp1 = 0
          ixapxbp_pr1 = 0
          ixapxbp_pr2 = 0
        end select
#endif


        do k = 2, gr%nnzp-1, 1

          km1 = max( k-1, 1 )
          kp1 = min( k+1, gr%nnzp )

          rhs(k,1)
          ! RHS turbulent advection (ta) term.
     .    = term_ta_rhs( a1_zt(kp1), a1(k), a1_zt(k), 
     .                   wp3(kp1), wp3(k), wp2_zt(kp1), 
     .                   wp2_zt(k), wpxbp_zt(kp1), wpxbp_zt(k), 
     .                   wpxap_zt(kp1), wpxap_zt(k), gr%dzm(k), 
     .                   beta, wtol_sqd )
          ! RHS turbulent production (tp) term.
     .    + (1.0 - C5) 
     .       * term_tp( xam(kp1), xam(k), xam(kp1), xam(k),
     .                  wpxap(k), wpxap(k), gr%dzm(k) )
          ! RHS pressure term 1 (pr1) (and dissipation term 1 (dp1)).
     .    + term_pr1( C4, C14, xbp2(k), wp2(k), taum(k) )
          ! RHS pressure term 2 (pr2).
     .    + term_pr2( C5, grav, T0, wpthvp(k), wpxap(k), wpxbp(k), 
     .                xam(kp1), xam(k), xbm(kp1), xbm(k), gr%dzm(k) )

          ! RHS time tendency.
          if ( liter ) then 
            rhs(k,1) = rhs(k,1) + 1.0/dt * xap2(k)
          endif

#ifdef STATS
          if ( lstats_samp ) then

          ! Statistics: explicit contributions for up2 or vp2.
               
            call stat_modify_pt( ixapxbp_ta, k,
     .             term_ta_rhs( a1_zt(kp1), a1(k), a1_zt(k), 
     .                         wp3(kp1), wp3(k), wp2_zt(kp1), 
     .                         wp2_zt(k), wpxbp_zt(kp1), wpxbp_zt(k), 
     .                         wpxap_zt(kp1), wpxap_zt(k), gr%dzm(k), 
     .                         beta, wtol_sqd ), zm )
            
            if ( ixapxbp_dp1 > 0 ) then
              ! Note:  The function term_pr1 is the explicit component 
              !        of a semi-implicit solution to dp1 and pr1.
              ! Record the statistical contribution of the implicit
              ! component of term dp1 for up2 or vp2.  This will 
              ! overwrite anything set statistically in diag_var_lhs 
              ! for this term.
              tmp = term_dp1( (2.0/3.0)*C4, taum(k) )
              zmscr01(k) = -tmp
              ! Statistical contribution of the explicit component
              ! of term dp1 for up2 or vp2.
            endif

            call stat_modify_pt( ixapxbp_dp1, k, 
     .              term_pr1( C4, 0.0, xbp2(k), wp2(k), taum(k) ), zm )


            if ( ixapxbp_pr1 > 0 ) then
              ! Note:  The function term_pr1 is the explicit component 
              !        of a semi-implicit solution to dp1 and pr1.
              ! Statistical contribution of the implicit component
              ! of term pr1 for up2 or vp2.
              tmp = term_dp1( (1.0/3.0)*C14, taum(k) )
              zmscr11(k) = -tmp
              ! Statistical contribution of the explicit component
              ! of term pr1 for up2 or vp2.
            endif

            call stat_modify_pt( ixapxbp_pr1, k,
     .              term_pr1( 0.0, C14, xbp2(k), wp2(k), taum(k) ), zm )

            call stat_update_var_pt( ixapxbp_pr2, k,
     .          term_pr2
     .          ( C5, grav, T0, wpthvp(k), wpxap(k), wpxbp(k), 
     .            xam(kp1), xam(k), xbm(kp1), xbm(k), gr%dzm(k) ), zm )

            call stat_update_var_pt( ixapxbp_tp, k,
     .            (1.0 - C5) 
     .             * term_tp( xam(kp1), xam(k), xam(kp1), xam(k),
     .                        wpxap(k), wpxap(k), gr%dzm(k) ), zm )

          endif ! lstats_samp
#endif /*STATS*/

        enddo ! k=2..gr%nnzp-1


        ! Boundary Conditions
        ! These are set so that the sfc_var value of up2 or vp2
        ! can be used at the lowest boundary and the values of those
        ! variables can be set to 0 at the top boundary.
        ! Fixed-point boundary conditions are used for the variances.

        ! This boundary condition was changed by dschanen on 24 April 2007
        ! When we run prognostically we want to preserve the surface value.
!       if ( liter ) then 
!         rhs(1,1) = xap2(1) + 1.0/dt*xap2(1)
!         rhs(gr%nnzp,1) = 1.0/dt*xap2(gr%nnzp)
!       else
        rhs(1,1) = xap2(1)
        rhs(gr%nnzp,1) = 0.0
!       end if


        return
        end subroutine diag_var_uv_rhs

!===============================================================================
        subroutine diag_var_rhs( solve_type, dt, liter, a1, a1_zt, 
     .                           wp2_zt, wp3, wpxap, wpxap_zt,
     .                           wpxbp, wpxbp_zt, xam, xbm, xapxbp,
     .                           beta, wtol_sqd, rhs )

!       Description:
!       Explicit contributions to r_t'^2, th_l'^2, r_t'th_l', 
!                                 sclr'r_t', sclr'th_l', or sclr'^2.
!-----------------------------------------------------------------------

        use grid_class, only:
     .      gr ! Variable(s)
        
        use stats_prec, only: 
     .      time_prec ! Variable(s)
        
#ifdef STATS
        use stats_type, only:
     .      stat_modify_pt, stat_update_var_pt ! Procedure(s)
        
        use stats_hoc, only:
     .      irtp2_ta, ! Variable(s)
     .      irtp2_tp,
     .      ithlp2_ta,
     .      ithlp2_tp,
     .      irtpthlp_ta,
     .      irtpthlp_tp1,
     .      irtpthlp_tp2,
     .      zm,
     .      lstats_samp
#endif

        implicit none

        ! Input Variables
        character(len=*), intent(in) :: solve_type

        real(kind=time_prec), intent(in) ::
     .  dt        ! Timestep                                    [s]

        logical, intent(in) ::
     .  liter   ! Whether x is prognostic (T/F)

        real, dimension(gr%nnzp), intent(in) ::
     .  a1,       ! Scm-related term a_1 (momentum levels)      [-]
     .  a1_zt,    ! a_1 interpolated to thermodynamic levels    [-]
     .  wp2_zt,   ! w'^2 interpolated to thermodynamic levels   [m^2/s^2]
     .  wp3,      ! w'^3 (thermodynamic levels)                 [m^3/s^3]
     .  wpxap,    ! w'x_a' (momentum levels)                    [m/s {x_am units}]
     .  wpxap_zt, ! w'x_a' interpolated to thermodynamic levels [m/s {x_am units}]
     .  wpxbp,    ! w'x_b' (momentum levels)                    [m/s {x_bm units}]
     .  wpxbp_zt, ! w'x_b' interpolated to thermodynamic levels [m/s {x_bm units}]
     .  xam,      ! x_am (thermodynamic levels)                 [{x_am units}]
     .  xbm,      ! x_bm (thermodynamic levels)                 [{x_bm units}]
     .  xapxbp    ! x_a'x_b' (momentum levels)                  [{x_am units}
                  !                                                *{x_bm units}]

        real, intent(in) ::
     .  beta,     ! Model parameter beta                        [-]
     .  wtol_sqd  ! w wind component tolerance squared          [m^2/s^2]

        real, dimension(gr%nnzp,1), intent(out) ::
     .  rhs     ! Explicit contributions to x variance/covariance terms

        ! Local Variables

        ! Array indices
        integer :: k, kp1, km1

#ifdef STATS
        integer ::
     .  ixapxbp_ta,
     .  ixapxbp_tp,
     .  ixapxbp_tp1,
     .  ixapxbp_tp2


        select case ( trim( solve_type ) )
        case ( "rtp2" )
          ixapxbp_ta  = irtp2_ta
          ixapxbp_tp  = irtp2_tp
          ixapxbp_tp1 = 0
          ixapxbp_tp2 = 0
        case ( "thlp2" )
          ixapxbp_ta  = ithlp2_ta
          ixapxbp_tp  = ithlp2_tp
          ixapxbp_tp1 = 0
          ixapxbp_tp2 = 0
        case ( "rtpthlp" )
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
          ! RHS turbulent advection (ta) term.
     .    = term_ta_rhs( a1_zt(kp1), a1(k), a1_zt(k), 
     .                   wp3(kp1), wp3(k), wp2_zt(kp1), 
     .                   wp2_zt(k), wpxbp_zt(kp1), wpxbp_zt(k), 
     .                   wpxap_zt(kp1), wpxap_zt(k), gr%dzm(k), 
     .                   beta, wtol_sqd )
          ! RHS turbulent production (tp) term.
     .    + term_tp( xam(kp1), xam(k), xbm(kp1), xbm(k), 
     .               wpxbp(k), wpxap(k), gr%dzm(k) )

          ! RHS time tendency.
          if ( liter ) then 
            rhs(k,1) = rhs(k,1) + 1.0/dt * xapxbp(k)
          endif

#ifdef STATS
          if ( lstats_samp ) then

          ! Statistics: explicit contributions for rtp2, thlp2, or rtpthlp.

          call stat_modify_pt( ixapxbp_ta, k,
     .        term_ta_rhs( a1_zt(kp1), a1(k), a1_zt(k),
     .                     wp3(kp1), wp3(k), wp2_zt(kp1),
     .                     wp2_zt(k), wpxbp_zt(kp1), wpxbp_zt(k),
     .                     wpxap_zt(kp1), wpxap_zt(k), gr%dzm(k),
     .                     beta, wtol_sqd ), zm )

          ! rtp2/thlp2 case (1 turbulent production term)
          call stat_update_var_pt( ixapxbp_tp, k,
     .        term_tp( xam(kp1), xam(k), xbm(kp1), xbm(k),
     .                 wpxbp(k), wpxap(k), gr%dzm(k) ), zm )
          
          ! rtpthlp case (2 turbulent production terms)
          call stat_update_var_pt( ixapxbp_tp1, k,
     .        term_tp( 0.0, 0.0, xbm(kp1), xbm(k),
     .                 0.0, wpxap(k), gr%dzm(k) ), zm )
          
          call stat_update_var_pt( ixapxbp_tp2, k,
     .        term_tp( xam(kp1), xam(k), 0.0, 0.0,
     .                 wpxbp(k), 0.0, gr%dzm(k) ), zm )
          
          endif ! lstats_samp
#endif /*STATS*/

        enddo ! k=2..gr%nnzp-1


        ! Boundary Conditions
        ! These are set so that the sfc_var value of rtp2, thlp2, or rtpthlp
        ! (or sclrp2, sclrprtp, or sclrpthlp) can be used at the lowest boundary
        ! and the values of those variables can be set to 0 at the top boundary.
        ! Fixed-point boundary conditions are used for both the variances and 
        ! the covariances.

        ! This boundary condition was changed by dschanen on 24 April 2007
        ! When we run prognostically we want to preserve the surface value.
!       if ( .false. ) then
!         rhs(1,1) = xapxbp(1) + 1.0/dt*xapxbp(1)
!         rhs(gr%nnzp,1) = 1.0/dt*xapxbp(gr%nnzp)
!       else
        rhs(1,1) = xapxbp(1)
        rhs(gr%nnzp,1) = 0.0
!       end if

        return
        end subroutine diag_var_rhs

!===============================================================================
        pure function term_ta_lhs( a1_ztp1, a1, a1_zt, 
     .                             wp3p1, wp3, wp2_ztp1, wp2_zt, 
     .                             dzm, beta, wtol_sqd, level )
     .  result( lhs )

!       Description:
!       Turbulent advection of x_a'x_b':  implicit portion of the code.
!
!       The d(x_a'x_b')/dt equation contains a turbulent advection term:
!
!       - d(w'x_a'x_b')/dz.
!
!       A substitution is made in order to close the turbulent advection
!       term, such that:
!
!       w'x_a'x_b' = (1/3)*beta * a_1 * ( w'^3 / w'^2 ) * x_a'x_b'
!                       + (1-(1/3)*beta) * (a_1)^2 * ( w'^3 / (w'^2)^2 ) 
!                         * w'x_a' * w'x_b',
!
!       where a_1 is a variable that is a function of Sc.
!       The turbulent advection term is rewritten as:
!
!       - d [ (1/3)*beta * a_1 * ( w'^3 / w'^2 ) * x_a'x_b'
!                + (1-(1/3)*beta) * (a_1)^2 * ( w'^3 / (w'^2)^2 ) 
!                  * w'x_a' * w'x_b' ]
!         / dz,
!
!       which produces an implicit and an explicit portion of this term.
!       The implicit portion of this term is:
!
!       - d [ (1/3)*beta * a_1 * ( w'^3 / w'^2 ) * x_a'x_b'(t+1) ] / dz.
!
!       Since (1/3)*beta is a constant, it can be pulled outside of the
!       derivative.  The implicit portion of this term becomes:
!
!       - (1/3)*beta * d [ a_1 * ( w'^3 / w'^2 ) * x_a'x_b'(t+1) ] / dz.
!
!       Note:  When the term is brought over to the left-hand side, the
!              sign is reversed and the leading "-" in front of the
!              term is changed to a "+".
!
!       The timestep index (t+1) means that the value of x_a'x_b' being 
!       used is from the next timestep, which is being advanced to in
!       solving the d(x_a'x_b')/dt equation.
!
!       The implicit portion of this term is discretized as follows:
!
!       The values of x_a'x_b' are found on the momentum levels, as are 
!       the values of w'^2 and a_1.  The values of w'^3 are found on the
!       thermodynamic levels.  The variables x_a'x_b', w'^2, and a_1 are 
!       each interpolated to the intermediate thermodynamic levels.  The
!       values of the mathematical expression (called F here) within the
!       dF/dz term are computed on the thermodynamic levels.  Then the
!       derivative (d/dz) of the expression (F) is taken over the 
!       central momentum level, yielding the desired result.
!       In this function, the values of F are as follows:
!
!       F = a_1(t) * ( w'^3(t) / w'^2(t) ) * x_a'x_b'(t+1),
!
!       where the timestep index (t) stands for the index of the current
!       timestep.
!
!
!       =a1p1========wp2p1========xapxbpp1======================= m(k+1)
!
!       ---a1(interp)---wp2(interp)---xapxbp(interp)---wp3p1----- t(k+1)
!
!       =a1==========wp2==========xapxbp=================dF/dz=== m(k)
!
!       ---a1(interp)---wp2(interp)---xapxbp(interp)---wp3------- t(k)
!
!       =a1m1========wp2m1========xapxbpm1======================= m(k-1)
!
!       The vertical indices m(k+1), t(k+1), m(k), t(k), and m(k-1)
!       correspond with altitudes zm(k+1), zt(k+1), zm(k), zt(k),
!       and zm(k-1), respectively.  The letter "t" is used for
!       thermodynamic levels and the letter "m" is used for momentum
!       levels.
!
!       dzm(k) = 1 / ( zt(k+1) - zt(k) )

!       References:
!-----------------------------------------------------------------------

        use grid_class, only: ! gr%weights_zm2zt
     .      gr ! Variable(s)

        implicit none

        ! External
        intrinsic :: max

        ! Constant parameters
        integer, parameter ::
     .  kp1_mdiag = 1,   ! Momentum superdiagonal index.
     .  k_mdiag   = 2,   ! Momentum main diagonal index.
     .  km1_mdiag = 3    ! Momentum subdiagonal index.

        integer, parameter ::
     .  m_above = 1,   ! Index for upper momentum level grid weight.
     .  m_below = 2    ! Index for lower momentum level grid weight.

        ! Input Variables
        real, intent(in) ::
     .  a1_ztp1,    ! a_1 interpolated to thermo. level (k+1)  [-]
     .  a1,         ! a_1(k)                                   [-]
     .  a1_zt,      ! a_1 interpolated to thermo. level (k)    [-]
     .  wp3p1,      ! w'^3(k+1)                                [m^3/s^3]
     .  wp3,        ! w'^3(k)                                  [m^3/s^3]
     .  wp2_ztp1,   ! w'^2 interpolated to thermo. level (k+1) [m^2/s^2]
     .  wp2_zt,     ! w'^2 interpolated to thermo. level (k)   [m^2/s^2]
     .  dzm,        ! Inverse of grid spacing                  [1/m]
     .  beta,       ! Model parameter                          [-]
     .  wtol_sqd    ! w wind component tolerance squared       [m^2/s^2]

        integer, intent(in) ::
     .  level ! Central momentum level (on which calculation occurs).

        ! Return Variable
        real, dimension(3) :: lhs

        ! Local Variables
        integer ::
     .  tkp1, ! Thermodynamic level directly above central momentum level.
     .  tk    ! Thermodynamic level directly below central momentum level.

        ! Thermodynamic level (k+1) is between momentum level (k+1)
        ! and momentum level (k).
        tkp1 = level + 1

        ! Thermodynamic level (k) is between momentum level (k)
        ! and momentum level (k-1).
        tk = level

        ! Brian tried a new discretization for the turbulent advection 
        ! term, for which the implicit portion of the term is:
        ! - d [ a_1 * (1/3)*beta * ( w'^3 / w'^2 ) * x_a'x_b' ] / dz.  
        ! In order to help stabilize x_a'x_b', a_1 has been pulled 
        ! outside the derivative.

        ! Momentum superdiagonal: [ x xapxbp(k+1,<t+1>) ]
        lhs(kp1_mdiag) 
!     .  = + (1.0/3.0) * beta * dzm
!     .      * a1_ztp1 * ( wp3p1 / max( wp2_ztp1, wtol_sqd ) ) 
!     .      * gr%weights_zm2zt(m_above,tkp1)
     .  = + (1.0/3.0) * beta * a1 * dzm
     .      * ( wp3p1 / max( wp2_ztp1, wtol_sqd ) ) 
     .      * gr%weights_zm2zt(m_above,tkp1)

        ! Momentum main diagonal: [ x xapxbp(k,<t+1>) ]
        lhs(k_mdiag) 
!     .  = + (1.0/3.0) * beta * dzm
!     .      * (   a1_ztp1 * ( wp3p1 / max( wp2_ztp1, wtol_sqd ) ) 
!     .            * gr%weights_zm2zt(m_below,tkp1)
!     .          - a1_zt * ( wp3 / max( wp2_zt, wtol_sqd ) )
!     .            * gr%weights_zm2zt(m_above,tk)
!     .        )
     .  = + (1.0/3.0) * beta * a1 * dzm
     .      * (   ( wp3p1 / max( wp2_ztp1, wtol_sqd ) ) 
     .            * gr%weights_zm2zt(m_below,tkp1)
     .          - ( wp3 / max( wp2_zt, wtol_sqd ) )
     .            * gr%weights_zm2zt(m_above,tk)
     .        )

        ! Momentum subdiagonal: [ x xapxbp(k-1,<t+1>) ]
        lhs(km1_mdiag) 
!     .  = - (1.0/3.0) * beta * dzm
!     .      * a1_zt * ( wp3 / max( wp2_zt, wtol_sqd ) )
!     .      * gr%weights_zm2zt(m_below,tk)
     .  = - (1.0/3.0) * beta * a1 * dzm
     .      * ( wp3 / max( wp2_zt, wtol_sqd ) )
     .      * gr%weights_zm2zt(m_below,tk)

        ! End of Brian's a1 change.  14 Feb 2008.

        return
        end function term_ta_lhs

!===============================================================================
        pure function term_ta_rhs( a1_ztp1, a1, a1_zt, 
     .                             wp3p1, wp3, wp2_ztp1, 
     .                             wp2_zt, wpxbp_ztp1, wpxbp_zt, 
     .                             wpxap_ztp1, wpxap_zt, dzm, 
     .                             beta, wtol_sqd )
     .  result( rhs )

!       Description:
!       Turbulent advection of x_a'x_b':  explicit portion of the code.
!
!       The d(x_a'x_b')/dt equation contains a turbulent advection term:
!
!       - d(w'x_a'x_b')/dz.
!
!       A substitution is made in order to close the turbulent advection
!       term, such that:
!
!       w'x_a'x_b' = (1/3)*beta * a_1 * ( w'^3 / w'^2 ) * x_a'x_b'
!                       + (1-(1/3)*beta) * (a_1)^2 * ( w'^3 / (w'^2)^2 )
!                         * w'x_a' * w'x_b',
!
!       where a_1 is a variable that is a function of Sc.
!       The turbulent advection term is rewritten as:
!
!       - d [ (1/3)*beta * a_1 * ( w'^3 / w'^2 ) * x_a'x_b'
!                + (1-(1/3)*beta) * (a_1)^2 * ( w'^3 / (w'^2)^2 )
!                  * w'x_a' * w'x_b' ]
!         / dz,
!
!       which produces an implicit and an explicit portion of this term.
!       The explicit portion of this term is:
!
!       - d [ (1-(1/3)*beta) * (a_1)^2 * ( w'^3 / (w'^2)^2 ) 
!             * w'x_a' * w'x_b' ] / dz.
!
!       Since (1-(1/3)*beta) is a constant, it can be pulled outside of 
!       the derivative.  The explicit portion of this term becomes:
!
!       - (1-(1/3)*beta) * d [ (a_1)^2 * ( w'^3 / (w'^2)^2 ) 
!                              * w'x_a' * w'x_b' ] / dz.
!
!       The explicit portion of this term is discretized as follows:
!
!       The values of w'x_a', w'x_b', w'^2, and a_1 are found on the 
!       momentum levels.  The values of w'^3 are found on the 
!       thermodynamic levels.  The variables w'x_a', w'x_b', w'^2, and 
!       a_1 are each interpolated to the intermediate thermodynamic 
!       levels.  The values of the mathematical expression (called F 
!       here) within the dF/dz term are computed on the thermodynamic 
!       levels.  Then the derivative (d/dz) of the expression (F) is 
!       taken over the central momentum level, yielding the desired 
!       result.
!       In this function, the values of F are as follows:
!
!       F = ( a_1(t) )^2 * ( w'^3(t) / ( w'^2(t) )^2 ) 
!           * w'x_a'(t) * w'x_b'(t),
!
!       where the timestep index (t) stands for the index of the current
!       timestep.
!
!
!       =a1p1=======wp2p1=======wpxapp1=======wpxbpp1============ m(k+1)
!
!       -a1(interp)-wp2(interp)-wpxap(interp)-wpxbp(interp)-wp3p1 t(k+1)
!
!       =a1=========wp2=========wpxap=========wpxbp====dF/dz===== m(k)
!
!       -a1(interp)-wp2(interp)-wpxap(interp)-wpxbp(interp)-wp3-- t(k)
!
!       =a1m1=======wp2m1=======wpxapm1=======wpxbpm1============ m(k-1)
!
!       The vertical indices m(k+1), t(k+1), m(k), t(k), and m(k-1)
!       correspond with altitudes zm(k+1), zt(k+1), zm(k), zt(k),
!       and zm(k-1), respectively.  The letter "t" is used for
!       thermodynamic levels and the letter "m" is used for momentum
!       levels.
!
!       dzm(k) = 1 / ( zt(k+1) - zt(k) )

!       References:
!-----------------------------------------------------------------------

        implicit none

        ! External
        intrinsic :: max

        ! Input variables
        real, intent(in) ::
     .  a1_ztp1,      ! a_1 interpolated to thermo. level (k+1)    [-]
     .  a1,           ! a_1(k)                                     [-]
     .  a1_zt,        ! a_1 interpolated to thermo. level (k)      [-]
     .  wp3p1,        ! w'^3(k+1)                                  [m^3/s^3]
     .  wp3,          ! w'^3(k)                                    [m^3/s^3]
     .  wp2_ztp1,     ! w'^2 interpolated to thermo. level (k+1)   [m^2/s^2]
     .  wp2_zt,       ! w'^2 interpolated to thermo. level (k)     [m^2/s^2]
     .  wpxbp_ztp1,   ! w'x_b' interpolated to thermo. level (k+1) [m/s {x_bm units}]
     .  wpxbp_zt,     ! w'x_b' interpolated to thermo. level (k)   [m/s {x_bm units}]
     .  wpxap_ztp1,   ! w'x_a' interpolated to thermo. level (k+1) [m/s {x_am units}]
     .  wpxap_zt,     ! w'x_a' interpolated to thermo. level (k)   [m/s {x_am units}]
     .  dzm,          ! Inverse of grid spacing                    [1/m]
     .  beta,         ! Model parameter                            [-]
     .  wtol_sqd      ! w wind component tolerance squared         [m^2/s^2]

        ! Return Variable
        real :: rhs

        ! Brian tried a new discretization for the turbulent advection 
        ! term, for which the explicit portion of the term is:
        ! - d [ (a_1)^2 * (1-(1/3)*beta) * ( w'^3 / (w'^2)^2 ) 
        !       * w'x_a' * w'x_b' ] / dz.  
        ! In order to help stabilize x_a'x_b', (a_1)^2 has been pulled 
        ! outside the derivative.

        rhs
!     .  = - ( 1.0 - (1.0/3.0) * beta ) * dzm
!     .      * (   a1_ztp1**2 * wpxap_ztp1 * wpxbp_ztp1
!     .            * ( wp3p1 / max( wp2_ztp1, wtol_sqd )**2 )
!     .          - a1_zt**2 * wpxap_zt * wpxbp_zt 
!     .            * ( wp3 / max( wp2_zt, wtol_sqd )**2 )
!     .        )
     .  = - ( 1.0 - (1.0/3.0) * beta ) * a1**2 * dzm
     .      * (   wpxap_ztp1 * wpxbp_ztp1
     .            * ( wp3p1 / max( wp2_ztp1, wtol_sqd )**2 )
     .          - wpxap_zt * wpxbp_zt
     .            * ( wp3 / max( wp2_zt, wtol_sqd )**2 )
     .        )

        ! End of Brian's a1 change.  14 Feb 2008.

        return
        end function term_ta_rhs

!===============================================================================
        pure function term_tp( xamp1, xam, xbmp1, xbm, 
     .                         wpxbp, wpxap, dzm )
     .  result( rhs )

!       Description:
!       Turbulent production of x_a'x_b':  explicit portion of the code.
!
!       The d(x_a'x_b')/dt equation contains a turbulent production term:
!
!       - w'x_b' d(x_am)/dz - w'x_a' d(x_bm)/dz.
!
!       This term is solved for completely explicitly and is
!       discretized as follows:
!
!       The values of w'x_a' and w'x_b' are found on the momentum levels,
!       whereas the values of x_am and x_bm are found on the 
!       thermodynamic levels.  The derivatives of both x_am and x_bm are
!       taken over the intermediate (central) momentum level.  All of 
!       the remaining mathematical operations take place at the central
!       momentum level, yielding the desired result.
!
!       -------xamp1------------xbmp1---------------------------- t(k+1)
!
!       =wpxap======d(xam)/dz=========d(xbm)/dz===wpxbp========== m(k)
!
!       -------xam--------------xbm------------------------------ t(k)
!
!       The vertical indices t(k+1), m(k), and t(k) correspond with
!       altitudes zt(k+1), zm(k), and zt(k), respectively.  The letter
!       "t" is used for thermodynamic levels and the letter "m" is used
!       for momentum levels.
!
!       dzm(k) = 1 / ( zt(k+1) - zt(k) )

!       References:
!-----------------------------------------------------------------------

        implicit none

        ! Input variables
        real, intent(in) ::
     .  xam,     ! x_am(k)                     [{x_am units}]
     .  xamp1,   ! x_am(k+1)                   [{x_am units}]
     .  xbm,     ! x_bm(k)                     [{x_bm units}]
     .  xbmp1,   ! x_bm(k+1)                   [{x_bm units}]
     .  wpxbp,   ! w'x_b'(k)                   [m/s {x_bm units}]
     .  wpxap,   ! w'x_a'(k)                   [m/s {x_am units}]
     .  dzm      ! Inverse of grid spacing (k) [1/m]

        ! Return Variable
        real :: rhs

        rhs
     .  = - wpxbp * dzm * ( xamp1 - xam )
     .    - wpxap * dzm * ( xbmp1 - xbm )

        return
        end function term_tp

!===============================================================================
        pure function term_dp1( Cn, taum ) 
     .  result( lhs )

!       Description:
!       Dissipation term 1 for x_a'x_b':  implicit portion of the code.
!
!       The d(x_a'x_b')/dt equation contains dissipation term 1:
!
!       - ( C_n / tau_m ) x_a'x_b'.
!
!       This term is solved for completely implicitly, such that:
!
!       - ( C_n / tau_m ) x_a'x_b'(t+1).
!
!       Note:  When the term is brought over to the left-hand side, the
!              sign is reversed and the leading "-" in front of the
!              term is changed to a "+".
!
!       The timestep index (t+1) means that the value of x_a'x_b' being 
!       used is from the next timestep, which is being advanced to in 
!       solving the d(x_a'x_b')/dt equation.
!
!       The values of x_a'x_b' are found on momentum levels.  The values 
!       of time-scale tau_m are also found on momentum levels.
!
!       Note: For equations that use pressure term 1 (such as the 
!             equations for u'^2 and v'^2), C_n = ( 2*C_4 + C_14 ) / 3;
!             which combines the implicit contributions for dissipation
!             term 1 and pressure term 1 into one expression.
!             Otherwise, C_n = C_2.

!       References:
!-----------------------------------------------------------------------

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  k_mdiag   = 1   ! Momentum main diagonal index.

        real, intent(in) ::
     .  Cn,    ! Coefficient C_n                       [-]
     .  taum   ! Time-scale tau at momentum levels (k) [s]

        real :: lhs

        ! Momentum main diagonal: [ x xapxbp(k,<t+1>) ]
        lhs 
     .  = + Cn / taum

        return
        end function term_dp1

!===============================================================================
        pure function term_pr1( C4, C14, xbp2, wp2, taum )
     .  result( rhs )

!       Description:
!       Pressure term 1 for x_a'x_b':  explicit portion of the code.
!
!       Note:  Pressure term 1 is only used when x_a'x_b' is either 
!              u'^2 or v'^2.  For the following description, the 
!              equation for u'^2 is used as the example.  The equation
!              for v'^2 is the same as the equation for u'^2, except
!              that the v'^2 and u'^2 variables are switched.
!
!       The d(u'^2)/dt equation contains dissipation term 1:
!
!       - ( C_4 / tau_m ) * ( u'^2 - (2/3)*em );
!
!       and pressure term 1:
!
!       - (2/3) * ( C_14 / tau_m ) * em,
!
!       where em = (1/2) * ( u'^2 + v'^2 + w'^2 ).
!
!       This simplifies to:
!
!       - [ ( 2*C_4 + C_14 ) / ( 3 * tau_m ) ] * u'^2
!          + [ ( C_4 - C_14 ) / ( 3 * tau_m ) ] * ( v'^2 + w'^2 ).
!
!       The combined term has both implicit and explicit components.
!       The implicit component is:
!
!       - [ ( 2*C_4 + C_14 ) / ( 3 * tau_m ) ] * u'^2(t+1).
!
!       Note:  When the implicit term is brought over to the left-hand
!              side, the sign is reversed and the leading "-" in front
!              of the term is changed to a "+".
!
!       Timestep index (t) stands for the index of the current timestep,
!       while timestep index (t+1) stands for the index of the next 
!       timestep, which is being advanced to in solving the 
!       d(x_a'x_b')/dt equation.
!
!       The implicit component of the combined dp1 and pr1 term is 
!       solved in function "term_dp1" above, 
!       where "( 2*C_4 + C_14 ) / 3" is sent in as "C_n".
!
!       The explicit component of the combined dp1 and pr1 term is:
!
!       + [ ( C_4 - C_14 ) / ( 3 * tau_m ) ] * ( v'^2(t) + w'^2(t) );
!
!       and is discretized as follows:
!
!       The values for v'^2 and w'^2, as well as for tau_m, are found on
!       the momentum levels.  The mathematical operations all take place
!       on the momentum levels, yielding the desired result.

!       References:
!-----------------------------------------------------------------------

        implicit none

        real, intent(in) ::
     .  C4,     ! Model parameter C_4                         [-]
     .  C14,    ! Model parameter C_14                        [-]
     .  xbp2,   ! v'^2(k) (if solving for u'^2) or vice versa [m^2/s^2]
     .  wp2,    ! w'^2(k)                                     [m^2/s^2]
     .  taum    ! Time-scale tau at momentum levels (k)       [s]

        real :: rhs

        rhs = + 1.0/3.0 * ( C4 - C14 ) * ( xbp2 + wp2 ) / taum

        return
        end function term_pr1

!===============================================================================
        pure function term_pr2( C5, grav, T0, wpthvp, upwp, vpwp, 
     .                          ump1, um, vmp1, vm, dzm ) 
     .  result( rhs )

!       Description:
!       Pressure term 2 for x_a'x_b':  explicit portion of the code.
!
!       Note:  Pressure term 2 is only used when x_a'x_b' is either
!              u'^2 or v'^2.  For the following description, the
!              equation for u'^2 is used as the example.  The equation
!              for v'^2 is the same as the equation for u'^2.
!
!       The d(u'^2)/dt equation contains pressure term 2:
!
!       + (2/3) C_5 [ (g/th_0) w'th_v' - u'w' du/dz - v'w' dv/dz ].
!
!       This term is solved for completely explicitly and is
!       discretized as follows:
!
!       The values of w'th_v', u'w', and v'w' are found on the momentum
!       levels, whereas the values of um and vm are found on the
!       thermodynamic levels.  The derivatives of both um and vm are
!       taken over the intermediate (central) momentum level.  All the
!       remaining mathematical operations take place at the central
!       momentum level, yielding the desired result.
!
!       ------ump1------------vmp1------------------------------- t(k+1)
!
!       =upwp======d(um)/dz========d(vm)/dz===vpwp=====wpthvp==== m(k)
!
!       ------um--------------vm--------------------------------- t(k)
!
!       The vertical indices t(k+1), m(k), and t(k) correspond with
!       altitudes zt(k+1), zm(k), and zt(k), respectively.  The letter
!       "t" is used for thermodynamic levels and the letter "m" is used
!       for momentum levels.
!
!       dzm(k) = 1 / ( zt(k+1) - zt(k) )

!       References:
!-----------------------------------------------------------------------

        implicit none

        real, intent(in) ::
     .  C5,       ! Model parameter C_5             [-]
     .  grav,     ! Gravitational acceleration      [m/s^2]
     .  T0,       ! Reference temperature           [K]
     .  wpthvp,   ! w'th_v'(k)                      [m/K/s]
     .  upwp,     ! u'w'(k)                         [m^2/s^2]
     .  vpwp,     ! v'w'(k)                         [m^2/s^2]
     .  ump1,     ! um(k+1)                         [m/s]
     .  um,       ! um(k)                           [m/s]
     .  vmp1,     ! vm(k+1)                         [m/s]
     .  vm,       ! vm(k)                           [m/s]
     .  dzm       ! Inverse of the grid spacing (k) [1/m]

        real :: rhs

        ! As applied to w'2
        rhs = + (2.0/3.0) * C5
     .          * ( ( grav / T0 ) * wpthvp
     .              -upwp * dzm * ( ump1 - um )
     .              -vpwp * dzm * ( vmp1 - vm ) )

        return
        end function term_pr2

!===============================================================================

        end module diagnose_variances
