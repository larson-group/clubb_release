!-----------------------------------------------------------------------
! $Id: mixing.F,v 1.59 2008-06-18 16:43:56 dschanen Exp $

        ! Passive scalars
#define SCLR_THETA 1
#define SCLR_RT 2

!===============================================================================
        module mixing

!       Description:
!       Contains the HOC mixing scheme.  Allows for optional variables.

!       References:
!       None
!-----------------------------------------------------------------------

        implicit none

        private ! Default scope

        public  :: timestep_mixing

        private :: mixing_lhs,
     .             mixing_rhs,
     .             mixing_solve,
     .             xm_term_ta_lhs,
     .             wpxp_term_ta_lhs,
     .             wpxp_term_tp_lhs,
     .             wpxp_terms_ac_pr2_lhs,
     .             wpxp_term_pr1_lhs,
     .             wpxp_terms_bp_pr3_rhs

        ! Parameter Constants
        integer, parameter, private ::
     .  nsub = 2,  ! Number of subdiagonals in the LHS matrix
     .  nsup = 2,  ! Number of superdiagonals in the LHS matrix
     .  nrhs = 1   ! Number of RHS vectors


        contains

!===============================================================================
        subroutine timestep_mixing( dt, Scm, wmm, wmt, wp2, wp3,
     .                              Kht, taum, Skwm, rtpthvp, 
     .                              rtm_forcing, thlpthvp, 
     .                              thlm_forcing, rtp2, thlp2,
     .                              rtm, wprtp, thlm, wpthlp,
     .                              implemented, err_code,
     .                              sclrpthvp, sclrm_forcing,
     .                              sclrp2, sclrm, wpsclrp )

!       Description:
!       Advance the mean and flux terms by one timestep.

!       References:
!       Eqn. 16 & 17 on p. 3546 of 
!       ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!         Method and Model Description'' Golaz, et al. (2002)
!       JAS, Vol. 59, pp. 3540--3551.

!       See Also
!       ``Equations for HOC'' Section 5:
!       /Implicit solutions for the means and fluxes/
!-----------------------------------------------------------------------

        use parameters, only: 
     .      C6rt, ! Variable(s)
     .      C6rtb, 
     .      C6rtc, 
     .      C6thl, 
     .      C6thlb, 
     .      C6thlc,
     .      C7, 
     .      C7b, 
     .      C7c, 
     .      c_K6, 
     .      c_Ksqd

        use constants, only: 
     .      sclr_dim, ! Variable(s)
     .      fstderr

        use grid_class, only:
     .      gr, ! Variable(s)
     .      zm2zt ! Procedure(s)

        use stats_prec, only: 
     .      time_prec ! Variable(s)
        use error_code, only: 
     .      lapack_error, ! Procedure(s)
     .      clubb_at_debug_level    
#ifdef STATS
        use stats_type, only:
     .      stat_begin_update, stat_end_update ! Procedure(s)
        
        use stats_hoc, only:
     .      zt, ! Variable(s)
     .      irtm_cl, 
     .      ithlm_cl,
     .      lstats_samp
#endif /*STATS*/
        
        implicit none

        ! External
        intrinsic :: exp, present

        ! Input Variables
        real(kind=time_prec), intent(in) :: 
     .  dt            ! Timestep                                 [s]

        real, intent(in), dimension(gr%nnzp) ::
     .  Scm,          ! Sc on momentum levels                    [-]
     .  wmm,          ! w wind component on momentum levels      [m/s]
     .  wmt,          ! w wind component on thermodynamic levels [m/s]
     .  wp2,          ! w'^2 (momentum levels)                   [m^2/s^2]
     .  wp3,          ! w'^3 (thermodynamic levels)              [m^3/s^3]
     .  Kht,          ! Eddy diffusivity on thermodynamic levels [m^2/s]
     .  taum,         ! Time-scale tau on momentum levels        [s]
     .  Skwm,         ! Skewness of w on momentum levels         [-]
     .  rtpthvp,      ! r_t'th_v' (momentum levels)              [(kg/kg) K]
     .  rtm_forcing,  ! r_t forcing (thermodynamic levels)       [(kg/kg)/s]
     .  thlpthvp,     ! th_l'th_v' (momentum levels)             [K^2]
     .  thlm_forcing, ! th_l forcing (thermodynamic levels)      [K/s]
! Added for clipping by Vince Larson 29 Sep 2007
     .  rtp2,         ! r_t'^2 (momentum levels)                 [(kg/kg)^2]
     .  thlp2         ! th_l'^2 (momentum levels)                [K^2]
! End of Vince Larson's addition.

        logical, intent(in) :: 
     .  implemented   ! Flag for CLUBB being implemented in a larger model.
         
        ! Input/Output Variables
        real, intent(inout), dimension(gr%nnzp) :: 
     .  rtm,         ! r_t  (total water mixing ratio)           [kg/kg]
     .  wprtp,       ! w'r_t'                                    [(kg/kg) m/s]
     .  thlm,        ! th_l (liquid water potential temperature) [K]
     .  wpthlp       ! w'th_l'                                   [K m/s]

        integer, intent(inout) :: err_code ! Model status 

        ! Additional optional variables for passive scalars
        ! Input Variables
        real, intent(in), optional, dimension(gr%nnzp,sclr_dim) :: 
     .  sclrpthvp, sclrm_forcing !                           [Units vary]
     .  ,sclrp2                  ! For clipping Vince Larson [Units vary]
        ! Input/Output Variables
        real, intent(inout), optional, dimension(gr%nnzp,sclr_dim) :: 
     .  sclrm, wpsclrp !                                     [Units vary]

        ! Local variables
        real, dimension(nsup+nsub+1,2*gr%nnzp) ::
     .  lhs  ! Implicit contributions to wpxp/xm (band diag. matrix) (LAPACK)

        real, dimension(2*gr%nnzp,nrhs) ::
     .  rhs  ! Right-hand side of band diag. matrix. (LAPACK)

        ! Constant parameters as a function of Skw.
        real, dimension(gr%nnzp) :: 
     .  C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc

        ! Eddy Diffusion for wpthlp and wprtp.
        real, dimension(gr%nnzp) :: Kw6   ! wpxp eddy diff. [m^2/s]

        ! Variables used for adding (wpxp)^2: 3-point average
        ! diffusion coefficient.
        real, dimension(gr%nnzp) ::
     .  wprtp_zt,
     .  wpthlp_zt,
     .  wprtp_zt_sqd_3pt,
     .  wpthlp_zt_sqd_3pt,
     .  Kw6_rt,
     .  Kw6_thl

        ! Indices
        integer :: i
        integer :: k, km1, kp1

!-----------------------------------------------------------------------

        ! Compute C6 and C7 as a function of Skw
        C6rt_Skw_fnc(:) = C6rtb + (C6rt-C6rtb)
     .                        *EXP( -0.5 * (Skwm(:)/C6rtc)**2 )

        C6thl_Skw_fnc(:) = C6thlb + (C6thl-C6thlb)
     .                          *EXP( -0.5 * (Skwm(:)/C6thlc)**2 )

        C7_Skw_fnc(:) = C7b + (C7-C7b)
     .                    *EXP( -0.5 * (Skwm(:)/C7c)**2 )

!        C6rt_Skw_fnc = C6rt
!        C6thl_Skw_fnc = C6thl
!        C7_Skw_fnc = C7

        ! (wpxp)^2: 3-point average diffusion coefficient.

        ! Interpolate w'x' (w'r_t' and w'th_l') from momentum levels 
        ! to thermodynamic levels.  This is used for extra diffusion
        ! based on a three-point average of (w'x')^2.
        wprtp_zt  = zm2zt( wprtp )
        wpthlp_zt = zm2zt( wpthlp )

        do k = 1, gr%nnzp, 1

           km1 = max( k-1, 1 )
           kp1 = min( k+1, gr%nnzp )

           ! Compute the square of wprtp_zt, averaged over 3 points.  26 Jan 2008
           wprtp_zt_sqd_3pt(k) = ( wprtp_zt(km1)**2 + wprtp_zt(k)**2
     .                            + wprtp_zt(kp1)**2 ) / 3.0
           ! Account for units of mix ratio (kg/kg)**2   Vince Larson 29 Jan 2008
           wprtp_zt_sqd_3pt(k) = 1e6 * wprtp_zt_sqd_3pt(k)

           ! Compute the square of wpthlp_zt, averaged over 3 points.  26 Jan 2008
           wpthlp_zt_sqd_3pt(k) = ( wpthlp_zt(km1)**2 + wpthlp_zt(k)**2
     .                             + wpthlp_zt(kp1)**2 ) / 3.0

        enddo

!       Define the Coefficent of Eddy Diffusivity for the wpthlp and wprtp.
        do k = 1, gr%nnzp, 1
           ! Kw6 is used for wpthlp and wprtp, which are located on 
           ! momentum levels.
           ! Kw6 is located on thermodynamic levels.
           ! Kw6 = c_K6 * Kht
           Kw6(k) = c_K6 * Kht(k)
           ! Kw6_rt must have units of m^2/s.  Since wprtp_zt_sqd_3pt has 
           ! units of m/s (kg/kg), c_Ksqd is given units of m/(kg/kg) 
           ! in this case.
           Kw6_rt(k)  = Kw6(k) + c_Ksqd * wprtp_zt_sqd_3pt(k)  ! Vince Larson increased by c_Ksqd, 29Jan2008
           ! Kw6_thl must have units of m^2/s.  Since wpthlp_zt_sqd_3pt has 
           ! units of m/s K, c_Ksqd is given units of m/K in this case.
           Kw6_thl(k) = Kw6(k) + c_Ksqd * wpthlp_zt_sqd_3pt(k) ! Vince Larson increased by c_Ksqd, 29Jan2008
        enddo

        ! Setup and decompose matrix

        ! Compute the implicit portion of the rtm and w'rt' equations.
        ! Build the left-hand side matrix.
        call mixing_lhs( .true., dt, Scm, wmm, wmt, wp2, wp3,
     .                   Kw6_rt, taum, C7_Skw_fnc, C6rt_Skw_fnc,
     .                   implemented, lhs )

        ! Compute the explicit portion of the rtm and w'rt' equations.
        ! Build the right-hand side vector.
        call mixing_rhs( "rtm", .true., dt, rtm, wprtp,
     .                   rtm_forcing, C7_Skw_fnc, rtpthvp, rhs )

        ! Solve rtm / w'rt'
        call mixing_solve( "rtm", dt, wp2, rtp2, lhs, rhs,
     .                     rtm, wprtp, err_code )

        if ( lapack_error( err_code ) )  then
     
          write(fstderr,'(a)') "rt mixing failed"
          return
        endif
 
        ! Clipping rtm
        ! Joshua Fasching March 2008
        
#ifdef STATS
        ! Computed value before clipping
        if ( lstats_samp ) then
           call stat_begin_update( irtm_cl, real( rtm / dt ), zt )
        end if 
#endif        
        ! The arm_0003 case produces negative rtm near the tropopause.
        !    To avoid this, we clip rtm.  This is not a good solution,
        !    because it renders rtm non-conserved.  We should look into 
        !    a positive definite advection scheme.
        !    Vince Larson.  13 Nov 2007

        ! The clipping of rtm causes a spurious source of moisture,
        !   particularly in SAM-CLUBB mode, and so this code will be
        !   disabled by default for now.
        !   David Schanen 15 Apr 2008
        do k = 1, gr%nnzp, 1 
          if ( rtm(k) < 0.0 ) then
!           rtm(k) = 0.0
            if ( clubb_at_debug_level( 1 ) ) then
            write(fstderr,*) "rtm < 0 in mixing at k= ", k
            endif 
          endif

        enddo
                
#ifdef STATS
         if ( lstats_samp ) then
            call stat_end_update( irtm_cl, real( rtm / dt ) , zt )      
         endif
#endif

        ! Compute the implicit portion of the thlm and w'thl' equations.
        ! Build the left-hand side matrix.
        call mixing_lhs( .true., dt, Scm, wmm, wmt, wp2, wp3,
     .                   Kw6_thl, taum, C7_Skw_fnc, C6thl_Skw_fnc,
     .                   implemented, lhs )

        ! Compute the explicit portion of the thlm and w'thl' equations.
        ! Build the right-hand side vector.
        call mixing_rhs( "thlm", .true., dt, thlm, wpthlp,
     .                   thlm_forcing, C7_Skw_fnc, thlpthvp, rhs )

        ! Solve for thlm / w'thl'
        call mixing_solve( "thlm", dt, wp2, thlp2, lhs, rhs,
     .                     thlm, wpthlp, err_code )

        if ( lapack_error( err_code ) ) then       
          write(fstderr,'(a)') "thetal mixing failed"
          return
        endif
        
        ! Clipping thlm
        ! Joshua Fasching March 2008
        
#ifdef STATS
        ! Computed value before clipping
        if ( lstats_samp ) then
           call stat_begin_update( ithlm_cl, real(thlm / dt ), zt)
         endif 
#endif        

        ! The value of potential temperature cannot fall below 0,
        ! so we clip accordingly
        do k = 1, gr%nnzp, 1 
          if ( thlm(k) < 0.0 ) then
            thlm(k) = 0.0
            write(fstderr,*) "thlm < 0 in mixing at k= ", k 
          endif
        enddo
        
                
#ifdef STATS
         if ( lstats_samp ) then
            call stat_end_update( ithlm_cl, real( thlm/dt ), zt )
         end if
#endif
        ! End change Joshua Fasching March 2008
        
        ! Solve sclrm / wpsclrp
        if ( present( sclrm )     .and. present( wpsclrp ) .and.
     .       present( sclrpthvp ) .and. present( sclrm_forcing ) ) then

          ! Compute the implicit portion of the sclrm and w'sclr' equations.
          ! Build the left-hand side matrix.
          call mixing_lhs( .true., dt, Scm, wmm, wmt, wp2, wp3,
     .                     Kw6, taum, C7_Skw_fnc, C6rt_Skw_fnc,
     .                     implemented, lhs )

          ! Compute the explicit portion of the sclrm and w'sclr' equations.
          ! Build the right-hand side vector.
          call mixing_rhs( "sclrm_rtm", .true., dt, 
     .                     sclrm(:,SCLR_RT), wpsclrp(:,SCLR_RT), 
     .                     sclrm_forcing(:,SCLR_RT), C7_Skw_fnc, 
     .                     sclrpthvp(:,SCLR_RT), rhs )

          ! Solve for sclrm / w'sclr'
          call mixing_solve( "sclrm_rtm", dt, wp2, 
     .                       sclrp2(:,SCLR_RT), lhs, rhs,
     .                       sclrm(:,SCLR_RT), wpsclrp(:,SCLR_RT), 
     .                       err_code )

          if ( lapack_error( err_code ) ) then    
            write(fstderr,'(a)') "Passive scalar 1 mixing failed"
            return
          endif

          ! Compute the implicit portion of the sclrm and w'sclr' equations.
          ! Build the left-hand side matrix.
          call mixing_lhs( .true., dt, Scm, wmm, wmt, wp2, wp3,
     .                     Kw6, taum, C7_Skw_fnc, C6thl_Skw_fnc,
     .                     implemented, lhs )

          ! Compute the explicit portion of the sclrm and w'sclr' equations.
          ! Build the right-hand side vector.
          call mixing_rhs( "sclrm_theta", .true., dt, 
     .                     sclrm(:,SCLR_THETA), wpsclrp(:,SCLR_THETA), 
     .                     sclrm_forcing(:,SCLR_THETA), C7_Skw_fnc, 
     .                     sclrpthvp(:,SCLR_THETA), rhs )

          ! Solve for sclrm / w'sclr'
          call mixing_solve( "sclrm_theta", dt, wp2, 
     .                       sclrp2(:,SCLR_THETA), lhs, rhs,
     .                       sclrm(:,SCLR_THETA), wpsclrp(:,SCLR_THETA), 
     .                       err_code )

          if ( lapack_error( err_code ) ) then  
            write(fstderr,'(a)') "Passive scalar 2 mixing failed"
            return
          endif
        endif ! passive scalars
        
!       Error Report
!       (This code is unreachable)
!       Joshua Fasching Feb 2008 

!        if ( lapack_error( err_code ) ) then
                
!           write(fstderr,*) "Error in timestep_mixing"
           
!           write(fstderr,*) "Intent(in)"
           

!           write(fstderr,*) "dt = ", dt
!           write(fstderr,*) "taum = ", taum
!           write(fstderr,*) "wmm = ", wmm
!           write(fstderr,*) "wmt = ", wmt
!           write(fstderr,*) "wp2 = ", wp2
!           write(fstderr,*) "wp3 = ", wp3
!           write(fstderr,*) "Scm = ", Scm
!           write(fstderr,*) "Skwm = ", Skwm

           !write(fstderr,*) "dt = ", dt
           !write(fstderr,*) "Scm = ", Scm
           !write(fstderr,*) "wmm = ", wmm
           !write(fstderr,*) "wmt = ", wmt
           !write(fstderr,*) "wp2 = ", wp2
           !write(fstderr,*) "wp3 = ", wp3
           !write(fstderr,*) "Kht = ", Kht
           !write(fstderr,*) "taum = ", taum
           !write(fstderr,*) "Skwm = ", Skwm
           !write(fstderr,*) "rtpthvp = ", rtpthvp
           !write(fstderr,*) "rtm_forcing = ", rtm_forcing
           !write(fstderr,*) "thlpthvp = ", thlpthvp
           !write(fstderr,*) "thlm_forcing = ", thlm_forcing
           !write(fstderr,*) "rtp2 = ", rtp2
           !write(fstderr,*) "thlp2 = ", thlp2
           
!           if( present( sclrpthvp ) ) then
!              write(fstderr,*) "sclrpthvp = ", sclrpthvp
!           endif
           
!           if( present( sclrm_forcing ) ) then
!              write(fstderr,*) "sclrm_forcing = ", sclrm_forcing
!           endif
           
!           if( present( sclrp2 ) ) then
!              write(fstderr,*) "sclrp2 = ", sclrp2
!           endif
           
!           write(fstderr,*) "Intent(inout)"
           
!           write(fstderr,*) "rtm = ", rtm 
!           write(fstderr,*) "wprtp = ", wprtp
!           write(fstderr,*) "thlm = ", thlm
!           write(fstderr,*) "wpthlp =", wpthlp
           
!           if( present( sclrm ) ) then
!              write(fstderr,*) "sclrm = ", sclrm
!           endif
           
!           if( present( wpsclrp ) ) then
!              write(fstderr,*) "wpsclrp = ", wpsclrp
!           endif
           
!        end if
 
        return
        
        end subroutine timestep_mixing

!===============================================================================
        subroutine mixing_lhs( liter, dt, Scm, wmm, wmt, wp2, wp3,
     .                         Kw6, taum, C7_Skw_fnc, C6x_Skw_fnc,
     .                         implemented, lhs )

!       Description:
!       Compute LHS band diagonal matrix for xm and w'x'.
!       This subroutine computes the implicit portion of
!       the xm and w'x' equations.

!       References:
!------------------------------------------------------------------------

        use parameters, only: 
     .      nu6 ! Variable(s)
        
        use constants, only: 
     .      wtol ! Variable(s)

        use grid_class, only: 
     .      gr, ! Variable(s)
     .      zm2zt ! Procedure(s)

        use stats_prec, only: 
     .      time_prec ! Variable(s)

        use diffusion, only: 
     .      diffusion_zm_lhs ! Procedure(s)

        use mean_adv, only:
     .      term_ma_zt_lhs, ! Procedure(s)
     .      term_ma_zm_lhs
#ifdef STATS
        use stats_hoc, only:
     .      ztscr01, ! Variable(s)
     .      ztscr02,
     .      ztscr03,
     .      ztscr04,
     .      ztscr05,
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
     .      zmscr11,
     .      zmscr12,
     .      zmscr13,
     .      zmscr14,
     .      lstats_samp,
     .      ithlm_ma,
     .      ithlm_ta,
     .      irtm_ma,
     .      irtm_ta,
     .      iwpthlp_ma,
     .      iwpthlp_ta,
     .      iwpthlp_tp,
     .      iwpthlp_ac,
     .      iwpthlp_pr1,
     .      iwpthlp_pr2,
     .      iwpthlp_dp1,
     .      iwprtp_ma,
     .      iwprtp_ta,
     .      iwprtp_tp,
     .      iwprtp_ac,
     .      iwprtp_pr1,
     .      iwprtp_pr2,
     .      iwprtp_dp1
#endif /*STATS*/

        implicit none

        ! External
        intrinsic :: min, max

        ! Input variables
        logical, intent(in) :: liter

        real(kind=time_prec), intent(in) :: 
     .  dt          ! Timestep                                 [s]

        real, intent(in), dimension(gr%nnzp) :: 
     .  Scm,        ! Sc on momentum levels                    [-]
     .  wmm,        ! w wind component on momentum levels      [m/s]
     .  wmt,        ! w wind component on thermodynamic levels [m/s]
     .  wp2,        ! w'^2 (momentum levels)                   [m^2/s^2]
     .  wp3,        ! w'^3 (thermodynamic levels)              [m^3/s^3]
     .  Kw6,        ! Coefficient of eddy diffusivity for w'x' [m^2/s]
     .  taum,       ! Time-scale tau on momentum levels        [s]
     .  C7_Skw_fnc, ! C_7 parameter with Sk_w applied          [-]
     .  C6x_Skw_fnc ! C_6x parameter with Sk_w applied         [-]

        logical, intent(in) :: 
     .  implemented ! Flag for CLUBB being implemented in a larger model.

        ! Output Variable
        real, intent(out), dimension(nsup+nsub+1,2*gr%nnzp) :: 
     .  lhs ! Implicit contributions to wpxp/xm (band diag. matrix) (LAPACK)

        ! Local Variables
        real, dimension(gr%nnzp) :: 
     .  a1 ! a_1 (momentum levels); See eqn. 24 in `Equations for HOC' [-]

        real, dimension(gr%nnzp) ::
     .  a1_zt, ! a_1 interpolated to thermodynamic levels              [-]
     .  wp2_zt ! w'^2 interpolated to thermodynamic levels             [m^2/s^2]

        ! wtol_sqd = the square of the minimum threshold on w,
        !     [wtol_sqd] = m^2 s^{-2}.  Vince Larson 11 Mar 2008.
        real :: wtol_sqd

        ! Indices
        integer :: k, km1, kp1
        integer :: k_xm, k_wpxp

#ifdef STATS
        real, dimension(3) :: tmp
#endif /*STATS*/


        ! wtol_sqd = the square of the minimum threshold on w,
        !     [wtol_sqd] = m^2 s^{-2}.  Vince Larson 11 Mar 2008.
        wtol_sqd = wtol * wtol

        ! Define a_1 (located on momentum levels).
        ! It is a variable that is a function of Sc (where Scm is
        ! located on momentum levels).
        a1(1:gr%nnzp) = 1.0 / ( 1.0 - Scm(1:gr%nnzp) )

        ! Interpolate a_1 and w'^2 from momentum levels to thermodynamic
        ! levels.  This will be used for the w'x' turbulent advection
        ! (ta) term.
        a1_zt  = max( zm2zt( a1 ), 0.0 )   ! Positive definite quantity
        wp2_zt = max( zm2zt( wp2 ), 0.0 )   ! Positive definite quantity

        ! Initialize the left-hand side matrix to 0.
        lhs = 0.0

        do k = 2, gr%nnzp-1, 1

          ! Define indices

          km1 = max( k-1, 1 )
          kp1 = min( k+1, gr%nnzp )

          k_xm   = 2*k - 1
          k_wpxp = 2*k


          !!!!!***** xm *****!!!!!

          ! xm: Left-hand side (implicit xm portion of the code).
          !
          ! Thermodynamic subdiagonal (lhs index: 3+2)
          !         [ x xm(k-1,<t+1>) ]
          ! Momentum subdiagonal (lhs index: 3+1)
          !         [ x wpxp(k-1,<t+1>) ]
          ! Thermodynamic main diagonal (lhs index: 3)
          !         [ x xm(k,<t+1>) ]
          ! Momentum superdiagonal (lhs index: 3-1)
          !         [ x wpxp(k,<t+1>) ]
          ! Thermodynamic superdiagonal (lhs index: 3-2)
          !         [ x xm(k+1,<t+1>) ]

          ! LHS mean advection (ma) term.
          if ( .not. implemented ) then

            lhs((/3-2,3,3+2/),k_xm)
     .      = lhs((/3-2,3,3+2/),k_xm)
     .      + term_ma_zt_lhs( wmt(k), gr%dzt(k), k )

          else

            lhs((/3-2,3,3+2/),k_xm)
     .      = lhs((/3-2,3,3+2/),k_xm) + 0.0

          endif

          ! LHS turbulent advection (ta) term.
          lhs((/3-1,3+1/),k_xm)
     .    = lhs((/3-1,3+1/),k_xm)
     .    + xm_term_ta_lhs( gr%dzt(k) )

          ! LHS time tendency.
          lhs(3,k_xm)
     .    = lhs(3,k_xm) + 1.0 / dt

#ifdef STATS
          if (lstats_samp) then

          ! Statistics: implicit contributions for rtm or thlm.

          if ( irtm_ma > 0 .or. ithlm_ma > 0 ) then
            if ( .not. implemented ) then
              tmp(1:3) = 
     .        + term_ma_zt_lhs( wmt(k), gr%dzt(k), k )
              ztscr01(k) = - tmp(3)
              ztscr02(k) = - tmp(2)
              ztscr03(k) = - tmp(1)
            else
              ztscr01(k) = 0.0
              ztscr02(k) = 0.0
              ztscr03(k) = 0.0
            endif
          endif

          if ( irtm_ta > 0 .or. ithlm_ta > 0 ) then
            tmp(1:2) =
     .      + xm_term_ta_lhs( gr%dzt(k) )
            ztscr04(k) = - tmp(2)
            ztscr05(k) = - tmp(1)
          endif

          endif
#endif /*STATS*/


          !!!!!***** w'x' *****!!!!!

          ! w'x': Left-hand side (implicit w'x' portion of the code).
          !
          ! Momentum subdiagonal (lhs index: 3+2)
          !         [ x wpxp(k-1,<t+1>) ]
          ! Thermodynamic subdiagonal (lhs index: 3+1)
          !         [ x xm(k,<t+1>) ]
          ! Momentum main diagonal (lhs index: 3)
          !         [ x wpxp(k,<t+1>) ]
          ! Thermodynamic superdiagonal (lhs index: 3-1)
          !         [ x xm(k+1,<t+1>) ]
          ! Momentum superdiagonal (lhs index: 3-2)
          !         [ x wpxp(k+1,<t+1>) ]

          ! LHS mean advection (ma) term.
          lhs((/3-2,3,3+2/),k_wpxp)
     .    = lhs((/3-2,3,3+2/),k_wpxp)
     .    + term_ma_zm_lhs( wmm(k), gr%dzm(k), k )

          ! LHS turbulent advection (ta) term.
          lhs((/3-2,3,3+2/),k_wpxp)
     .    = lhs((/3-2,3,3+2/),k_wpxp)
     .    + wpxp_term_ta_lhs( wp2_zt(kp1), wp2_zt(k), 
     .                        a1_zt(kp1), a1_zt(k),
     .                        wp3(kp1), wp3(k), gr%dzm(k), 
     .                        wtol_sqd, k )

          ! LHS turbulent production (tp) term.
          lhs((/3-1,3+1/),k_wpxp)
     .    = lhs((/3-1,3+1/),k_wpxp)
     .    + wpxp_term_tp_lhs( wp2(k), gr%dzm(k) )

          ! LHS accumulation (ac) term and pressure term 2 (pr2).
          lhs(3,k_wpxp)
     .    = lhs(3,k_wpxp)
     .    + wpxp_terms_ac_pr2_lhs( C7_Skw_fnc(k), 
     .                             wmt(kp1), wmt(k), gr%dzm(k) )

          ! LHS pressure term 1 (pr1).
          lhs(3,k_wpxp)
     .    = lhs(3,k_wpxp)
     .    + wpxp_term_pr1_lhs( C6x_Skw_fnc(k), taum(k) )

          ! LHS eddy diffusion term: dissipation term 1 (dp1).
          lhs((/3-2,3,3+2/),k_wpxp)
     .    = lhs((/3-2,3,3+2/),k_wpxp)
     .    + diffusion_zm_lhs( Kw6(k), Kw6(kp1), nu6,
     .                        gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )

          ! LHS time tendency.
          if (liter) lhs(3,k_wpxp) = lhs(3,k_wpxp) + 1.0 / dt

#ifdef STATS
          if (lstats_samp) then

          ! Statistics: implicit contributions for wprtp or wpthlp.

          if ( iwprtp_ma > 0 .or. iwpthlp_ma > 0 ) then
            tmp(1:3) =
     .      + term_ma_zm_lhs( wmm(k), gr%dzm(k), k )
            zmscr01(k) = - tmp(3)
            zmscr02(k) = - tmp(2)
            zmscr03(k) = - tmp(1)
          endif

          if ( iwprtp_ta > 0 .or. iwpthlp_ta > 0 ) then
            tmp(1:3) =
     .      + wpxp_term_ta_lhs( wp2_zt(kp1), wp2_zt(k), 
     .                          a1_zt(kp1), a1_zt(k),
     .                          wp3(kp1), wp3(k), gr%dzm(k), 
     .                          wtol_sqd, k )
            zmscr04(k) = - tmp(3)
            zmscr05(k) = - tmp(2)
            zmscr06(k) = - tmp(1)
          endif

          if ( iwprtp_tp > 0 .or. iwpthlp_tp > 0 ) then
            tmp(1:2) =
     .      + wpxp_term_tp_lhs( wp2(k), gr%dzm(k) )
            zmscr07(k) = - tmp(2)
            zmscr08(k) = - tmp(1)
          endif

          if ( iwprtp_ac > 0 .or. iwpthlp_ac > 0 ) then
            zmscr09(k) = 
     .      - wpxp_terms_ac_pr2_lhs( 0.0,
     .                               wmt(kp1), wmt(k), gr%dzm(k) )
          endif

          if ( iwprtp_pr1 > 0 .or. iwpthlp_pr1 > 0 ) then
            zmscr10(k) =
     .      - wpxp_term_pr1_lhs( C6x_Skw_fnc(k), taum(k) )
          endif

          if ( iwprtp_pr2 > 0 .or. iwpthlp_pr2 > 0 ) then
            zmscr11(k) =
     .      - wpxp_terms_ac_pr2_lhs( (1.0+C7_Skw_fnc(k)),
     .                               wmt(kp1), wmt(k), gr%dzm(k) )
          endif

          if ( iwprtp_dp1 > 0 .or. iwpthlp_dp1 > 0 ) then
            tmp(1:3) =
     .      + diffusion_zm_lhs( Kw6(k), Kw6(kp1), nu6,
     .                          gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )
            zmscr12(k) = - tmp(3)
            zmscr13(k) = - tmp(2)
            zmscr14(k) = - tmp(1)
          endif

          endif
#endif /*STATS*/

        enddo ! 2..gr%nnzp-1


        ! Boundary conditions

        ! Both the mean (xm) and the turbulent flux (wpxp) use fixed-point 
        ! boundary conditions.  Therefore, anything set in the above loop 
        ! at both the upper and lower boundaries would be overwritten here.
        ! However, the above loop does not extend to the boundary levels.  
        ! An array with a value of 1 at the main diagonal on the left-hand 
        ! side and with values of 0 at all other diagonals on the left-hand
        ! side will preserve the right-hand side value at that level.
        !
        !   xm(1)  wpxp(1) ... xm(nz) wpxp(nz)
        ! [  0.0     0.0         0.0     0.0  ]
        ! [  0.0     0.0         0.0     0.0  ]
        ! [  1.0     1.0   ...   1.0     1.0  ]
        ! [  0.0     0.0         0.0     0.0  ]
        ! [  0.0     0.0         0.0     0.0  ]

        ! Lower boundary
        k      = 1
        k_xm   = 2*k - 1
        k_wpxp = 2*k

        ! xm
        lhs(:,k_xm)     = 0.0
        lhs(3,k_xm)     = 1.0
        ! w'x'
        lhs(:,k_wpxp)   = 0.0
        lhs(3,k_wpxp)   = 1.0

        ! Upper boundary
        k      = gr%nnzp
        k_xm   = 2*k - 1
        k_wpxp = 2*k

        ! xm
        lhs(:,k_xm)     = 0.0
        lhs(3,k_xm)     = 1.0
        ! w'x'
        lhs(:,k_wpxp)   = 0.0
        lhs(3,k_wpxp)   = 1.0


        return
        end subroutine mixing_lhs

!===============================================================================
        subroutine mixing_rhs( solve_type, liter, dt, xm, wpxp,
     .                         xm_forcing, C7_Skw_fnc, xpthvp, rhs )

!       Description:
!       Compute RHS vector for xm and w'x'.
!       This subroutine computes the explicit portion of
!       the xm and w'x' equations.

!       References:
!------------------------------------------------------------------------

        use grid_class, only:
     .      gr ! Variable(s)

        use stats_prec, only: 
     .      time_prec ! Variable(s)

#ifdef STATS
        use stats_type, only:
     .      stat_update_var_pt
        
        use stats_hoc, only:
     .      zt, ! Variable(s)
     .      zm,
     .      irtm_forcing,
     .      ithlm_forcing,
     .      iwprtp_bp,
     .      iwprtp_pr3,
     .      iwpthlp_bp,
     .      iwpthlp_pr3,
     .      lstats_samp

#endif /*STATS*/

        implicit none

        ! Input Variables
        character(len=*), intent(in) ::
     .  solve_type  ! Variables being solved for.

        logical, intent(in) :: liter

        real(kind=time_prec), intent(in) :: 
     .  dt          ! Timestep                            [s]

        real, dimension(gr%nnzp), intent(in) ::
     .  xm,         ! xm (thermodynamic levels)           [{xm units}]
     .  wpxp,       ! w'x' (momentum levels)              [{xm units} m/s]
     .  xm_forcing, ! xm forcings (thermodynamic levels)  [{xm units}/s]
     .  C7_Skw_fnc, ! C_7 parameter with Sk_w applied     [-]
     .  xpthvp      ! x'th_v' (momentum levels)           [{xm units} K]

        ! Output Variable
        real, intent(out), dimension(2*gr%nnzp,nrhs) :: 
     .  rhs  ! Right-hand side of band diag. matrix. (LAPACK)

        ! Local Variables.
        ! Indices
        integer :: k, k_xm, k_wpxp

#ifdef STATS
        integer ::
     .    ixm_f,
     .    iwpxp_bp,
     .    iwpxp_pr3

        select case ( trim( solve_type ) )
        case ( "rtm" )  ! rtm/wprtp budget terms
          ixm_f     = irtm_forcing
          iwpxp_bp  = iwprtp_bp
          iwpxp_pr3 = iwprtp_pr3
        case ( "thlm" ) ! thlm/wpthlp budget terms
          ixm_f     = ithlm_forcing
          iwpxp_bp  = iwpthlp_bp
          iwpxp_pr3 = iwpthlp_pr3
        case default    ! this includes the sclrm case
          ixm_f     = 0
          iwpxp_bp  = 0
          iwpxp_pr3 = 0
        end select
#endif /*STATS*/


        ! Initialize the right-hand side vector to 0.
        rhs = 0.0

        do k = 2, gr%nnzp-1, 1

          ! Define indices

          k_xm   = 2*k - 1
          k_wpxp = 2*k


          !!!!!***** xm *****!!!!!

          ! xm: Right-hand side (explicit xm portion of the code).

          ! RHS time tendency.
          rhs(k_xm,1)
     .    = rhs(k_xm,1) + xm(k) / dt

          ! RHS xm forcings.
          ! Note: xm forcings include the effects of microphysics,
          !       cloud water sedimentation, radiation, and any
          !       imposed forcings on xm.
          rhs(k_xm,1)
     .    = rhs(k_xm,1) + xm_forcing(k)

#ifdef STATS
          if ( lstats_samp ) then

          ! Statistics: explicit contributions for xm 
          !             (including microphysics/radiation).
           call stat_update_var_pt( ixm_f, k, xm_forcing(k), zt )
          endif ! lstats_samp
#endif /*STATS*/


          !!!!!***** w'x' *****!!!!!

          ! w'x': Right-hand side (explicit w'x' portion of the code).

          ! RHS buoyancy production (bp) term and pressure term 3 (pr3).
          rhs(k_wpxp,1)
     .    = rhs(k_wpxp,1)
     .    + wpxp_terms_bp_pr3_rhs( C7_Skw_fnc(k), xpthvp(k) )

          ! RHS time tendency.
          if ( liter ) rhs(k_wpxp,1) = rhs(k_wpxp,1) + wpxp(k) / dt

#ifdef STATS
          if ( lstats_samp ) then

          ! Statistics: explicit contributions for wpxp.
          call stat_update_var_pt( iwpxp_bp, k,
     .        wpxp_terms_bp_pr3_rhs( 0.0, xpthvp(k) ), zm )

          call stat_update_var_pt( iwpxp_pr3, k,
     .        wpxp_terms_bp_pr3_rhs( (1.0+C7_Skw_fnc(k)),xpthvp(k)), zm)

          endif ! lstats_samp
#endif /*STATS*/

        enddo ! k=2..gr%nnzp-1


        ! Boundary conditions.

        ! Both the mean (xm) and the turbulent flux (wpxp) use fixed-point
        ! boundary conditions.  Therefore, anything set in the above loop
        ! at both the upper and lower boundaries would be overwritten here.
        ! However, the above loop does not extend to the boundary levels.
        ! An array with a value of 1 at the main diagonal on the left-hand
        ! side and with values of 0 at all other diagonals on the left-hand
        ! side will preserve the right-hand side value at that level.

        ! Lower boundary
        k      = 1
        k_xm   = 2*k - 1
        k_wpxp = 2*k
        ! The value of xm at the lower boundary will remain the same.
        ! However, the value of xm at the lower boundary gets overwritten
        ! after the matrix is solved for the next timestep, such 
        ! that xm(1) = xm(2).
        rhs(k_xm,1)   = xm(k)
        ! The value of w'x' at the lower boundary will remain the same.
        ! The surface value of w'x' is set elsewhere 
        ! (case-specific information).
        rhs(k_wpxp,1) = wpxp(k)

        ! Upper boundary
        k      = gr%nnzp
        k_xm   = 2*k - 1
        k_wpxp = 2*k
        ! The value of xm at the upper boundary will remain the same.
        rhs(k_xm,1)   = xm(k)
        ! The value of w'x' at the upper boundary will be 0.
        rhs(k_wpxp,1) = 0.0


        end subroutine mixing_rhs

!===============================================================================
        subroutine mixing_solve( solve_type, dt, wp2, xp2, lhs, rhs,
     .                           xm, wpxp, err_code )

!       Description:
!       Solve for xm / w'x' using the band diagonal solver.

!       References:
!------------------------------------------------------------------------

        use grid_class, only:
     .      gr ! Variable(s)

        use stats_prec, only: 
     .      time_prec ! Variable(s)

        use lapack_wrap, only: 
     .      band_solve, ! Procedure(s)
     .      band_solvex

        use pos_definite_mod, only: 
     .      pos_definite_adj ! Procedure(s)

        use error_code, only: 
     .      lapack_error ! Procedure(s)

#ifdef STATS
        use stats_type, only:
     .      stat_begin_update, stat_update_var_pt, ! Procedure(s)
     .      stat_end_update, stat_update_var
        
        use stats_hoc, only:
     .      zt, ! Variable(s)
     .      zm,
     .      sfc,
     .      irtm_bt,
     .      irtm_ta,
     .      irtm_ma,
     .      irtm_cn,
     .      irtm_pd,
     .      iwprtp_bt,
     .      iwprtp_ma,
     .      iwprtp_ta,
     .      iwprtp_tp,
     .      iwprtp_ac,
     .      iwprtp_pr1,
     .      iwprtp_pr2,
     .      iwprtp_dp1,
     .      iwprtp_cl,
     .      iwprtp_pd,
     .      ithlm_bt,
     .      ithlm_ta,
     .      ithlm_ma,
     .      ithlm_cn,
     .      iwpthlp_bt,
     .      iwpthlp_ma,
     .      iwpthlp_ta,
     .      iwpthlp_tp,
     .      iwpthlp_ac,
     .      iwpthlp_pr1,
     .      iwpthlp_pr2,
     .      iwpthlp_dp1,
     .      iwpthlp_cl,
     .      lstats_samp,
     .      ztscr01,
     .      ztscr02,
     .      ztscr03,
     .      ztscr04,
     .      ztscr05,
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
     .      zmscr11,
     .      zmscr12,
     .      zmscr13,
     .      zmscr14
#endif /*STATS*/

        implicit none

        ! Input Variables
        character(len=*), intent(in) :: 
     .  solve_type  ! Variables being solved for.

        real(kind=time_prec), intent(in) :: 
     .  dt          ! Timestep                      [s]

        real, intent(in), dimension(gr%nnzp) :: 
! Vince Larson added for clipping 20070929
     .  wp2,        ! w'^2 (momentum levels)        [(m/s)^2]
     .  xp2         ! x'^2 (momentum levels)        [{xm units}^2]
! End Vince Larson's addition

        ! Input/Output Variables
        real, intent(inout), dimension(nsup+nsub+1,2*gr%nnzp) ::
     .  lhs  ! Implicit contributions to wpxp/xm (band diag. matrix) (LAPACK)

        real, intent(inout), dimension(2*gr%nnzp,nrhs) :: 
     .  rhs  ! Right-hand side of band diag. matrix. (LAPACK)

        real, intent(inout), dimension(gr%nnzp) :: 
     .  xm,  ! Mean term: xm (thermodynamic levels) [units vary] 
     .  wpxp ! Flux term: w'x' (momentum levels)    [{xm units} m/s]

        ! Output Variable
        integer, intent(out) :: err_code

        ! Local Variables
        real, target, dimension(2*gr%nnzp,nrhs) :: 
     .  solution ! Solution to band diagonal system

        real, dimension(gr%nnzp) ::
     .  xm_n ! Old value of xm for positive definite scheme     [units vary] 

        real, dimension(gr%nnzp) ::
     .  wpxp_pd, xm_pd ! Change in xm and wpxp due to the pos. def. scheme

        real :: rcond ! Est. of the reciprocal of the condition #

        ! Indices
        integer :: k, km1, kp1
        integer :: k_xm, k_wpxp

#ifdef STATS
        integer ::
     .    ixm_bt,
     .    ixm_ta,
     .    ixm_ma,
     .    ixm_mc,
     .    ixm_cn,
     .    ixm_pd,
     .    iwpxp_bt,
     .    iwpxp_ma,
     .    iwpxp_ta,
     .    iwpxp_tp,
     .    iwpxp_ac,
     .    iwpxp_pr1,
     .    iwpxp_pr2,
     .    iwpxp_dp1,
     .    iwpxp_cl,
     .    iwpxp_pd

        select case ( trim( solve_type ) )
        case ( "rtm" ) ! rtm/wprtp budget terms
          ixm_bt    = irtm_bt
          ixm_ta    = irtm_ta
          ixm_ma    = irtm_ma
          ixm_cn    = irtm_cn ! Applies to both the mean and flux terms
          ixm_pd    = irtm_pd
          iwpxp_bt  = iwprtp_bt
          iwpxp_ma  = iwprtp_ma
          iwpxp_ta  = iwprtp_ta
          iwpxp_tp  = iwprtp_tp
          iwpxp_ac  = iwprtp_ac
          iwpxp_pr1 = iwprtp_pr1
          iwpxp_pr2 = iwprtp_pr2
          iwpxp_dp1 = iwprtp_dp1
          iwpxp_cl  = iwprtp_cl
          iwpxp_pd  = iwprtp_pd
        case ( "thlm" ) ! thlm/wpthlp budget terms
          ixm_bt    = ithlm_bt
          ixm_ta    = ithlm_ta
          ixm_ma    = ithlm_ma
          ixm_cn    = ithlm_cn ! Applies to both the mean and flux terms
          ixm_pd    = 0
          iwpxp_bt  = iwpthlp_bt
          iwpxp_ma  = iwpthlp_ma
          iwpxp_ta  = iwpthlp_ta
          iwpxp_tp  = iwpthlp_tp
          iwpxp_ac  = iwpthlp_ac
          iwpxp_pr1 = iwpthlp_pr1
          iwpxp_pr2 = iwpthlp_pr2
          iwpxp_dp1 = iwpthlp_dp1
          iwpxp_cl  = iwpthlp_cl
          iwpxp_pd  = 0
        case default  ! this includes the sclrm case
          ixm_bt    = 0
          ixm_ta    = 0
          ixm_ma    = 0
          ixm_cn    = 0
          ixm_pd    = 0
          iwpxp_bt  = 0
          iwpxp_ma  = 0
          iwpxp_ta  = 0
          iwpxp_tp  = 0
          iwpxp_ac  = 0
          iwpxp_pr1 = 0
          iwpxp_pr2 = 0
          iwpxp_dp1 = 0
          iwpxp_cl  = 0
          iwpxp_pd  = 0
        end select
#endif /*STATS*/

#ifdef STATS
        if ( lstats_samp ) then

          ! xm total time tendency ( 1st calculation)
          call stat_begin_update( ixm_bt, real( xm /dt ), zt )

          ! wpxp total time tendency (1st calculation)
          call stat_begin_update( iwpxp_bt, real( wpxp / dt ), zm )

        end if ! lstats_samp
#endif /*STATS*/


#ifdef STATS
        if ( lstats_samp .and. ixm_cn > 0 ) then
          ! Perform LU decomp and solve system (LAPACK with diagnostics)
          call band_solvex( solve_type, nsup, nsub, 2*gr%nnzp, nrhs,
     .                      lhs, rhs, solution, rcond, err_code )

          ! Est. of the condition number of the mean/flux LHS matrix 
          call stat_update_var_pt( ixm_cn, 1, 1.0 / rcond, sfc )

        else
          ! Perform LU decomp and solve system (LAPACK)
          call band_solve( solve_type, nsup, nsub, 2*gr%nnzp, nrhs,
     .                     lhs, rhs, solution, err_code )
        end if
#else
        ! Perform LU decomp and solve system (LAPACK)
        call band_solve( solve_type, nsup, nsub, 2*gr%nnzp, nrhs,
     .                   lhs, rhs, solution, err_code )

#endif /*STATS*/

        ! Return if the solver has failed
        if ( lapack_error( err_code ) ) return


        ! Copy result into output arrays

        do k=1,gr%nnzp

          k_xm   = 2 * k - 1
          k_wpxp = 2 * k
          
          xm_n(k) = xm(k)

          xm(k)   = solution(k_xm,1)
          wpxp(k) = solution(k_wpxp,1)

        end do ! k=1..gr%nnzp

!       Boundary condition on xm

!        xm(1) = 2. * xm(2) - xm(3)
!        xm(gr%nnzp) = 2. * xm(gr%nnzp-1) - xm(gr%nnzp-2)
         xm(1) = xm(2)
!        xm(gr%nnzp) = xm(gr%nnzp-1)


#ifdef STATS
          if ( lstats_samp ) then

           do k = 2, gr%nnzp-1

             km1 = max( k-1, 1 )
             kp1 = min( k+1, gr%nnzp )

             ! Finalize implicit contributions for xm
             call stat_update_var_pt( ixm_ma, k,
     .           ztscr01(k) * xm(km1)
     .         + ztscr02(k) * xm(k)
     .         + ztscr03(k) * xm(kp1), zt )
             
             call stat_update_var_pt( ixm_ta, k,
     .           ztscr04(k) * wpxp(km1)
     .         + ztscr05(k) * wpxp(k), zt )

             ! Finalize implicit contributions for wpxp
             call stat_update_var_pt( iwpxp_ma, k,
     .           zmscr01(k) * wpxp(km1)
     .         + zmscr02(k) * wpxp(k)
     .         + zmscr03(k) * wpxp(kp1), zm )


             call stat_update_var_pt( iwpxp_ta, k,
     .           zmscr04(k) * wpxp(km1)
     .         + zmscr05(k) * wpxp(k)
     .         + zmscr06(k) * wpxp(kp1), zm )
             
             call stat_update_var_pt( iwpxp_tp, k,
     .           zmscr07(k) * xm(k)
     .         + zmscr08(k) * xm(kp1), zm )

             call stat_update_var_pt( iwpxp_ac, k,
     .           zmscr09(k) * wpxp(k), zm )

             call stat_update_var_pt( iwpxp_pr1, k,
     .           zmscr10(k) * wpxp(k), zm )

             call stat_update_var_pt( iwpxp_pr2, k,
     .           zmscr11(k) * wpxp(k), zm )
             
             call stat_update_var_pt( iwpxp_dp1, k,
     .           zmscr12(k) * wpxp(km1)
     .         + zmscr13(k) * wpxp(k)
     .         + zmscr14(k) * wpxp(kp1), zm )
             
           enddo ! 1..gr%nnzp

          endif ! lstats_samp
#endif /*STATS*/

        ! Apply a flux limiting positive definite scheme if the solution 
        ! for the mean field is negative and we're determining total water
        if ( trim( solve_type ) == "rtm" .and. 
     .       any( xm < 0.0 ) ) then

           call pos_definite_adj( dt, "zt", xm, wpxp,
     .                            xm_n, xm_pd, wpxp_pd )

        else
          ! For stats purposes
          xm_pd   = 0.0
          wpxp_pd = 0.0

        end if ! solve_type == "rtm" and rtm <n+1> less than 0

#ifdef STATS
        if ( lstats_samp ) then

          call stat_update_var( iwpxp_pd, wpxp_pd(1:gr%nnzp), zm )

          call stat_update_var( ixm_pd, xm_pd(1:gr%nnzp), zt )

        end if
#endif
       


! Vince Larson clipped turbulent flux (wpxp) to ensure 
!     that -1 < [correlation of w and x] < 1 .  20070929
#ifdef STATS
        if ( lstats_samp ) then
           call stat_begin_update( iwpxp_cl, real( wpxp / dt ), zm )
        end if 
#endif


        ! Clip wpxp to make sure -0.99 < correlation < 0.99
        where ( wpxp >  0.99 * sqrt( wp2 * xp2 ) )
     .          wpxp =  0.99 * sqrt( wp2 * xp2 ) 
        where ( wpxp < -0.99 * sqrt( wp2 * xp2 ) )
     .          wpxp = -0.99 * sqrt( wp2 * xp2 ) 

#ifdef STATS
        if ( lstats_samp ) then
           call stat_end_update( iwpxp_cl, real( wpxp / dt ), zm )
        end if
#endif
! End of Vince Larson's clipping.

#ifdef STATS
        if ( lstats_samp ) then
          ! xm time tendency (2nd calculation)
          call stat_end_update( ixm_bt, real( xm / dt ), zt )

          ! wpxp time tendency (2nd calculation)
          call stat_end_update( iwpxp_bt, real( wpxp / dt ), zm )

        endif
#endif /*STATS*/

        return
        end subroutine mixing_solve

!===============================================================================
        pure function xm_term_ta_lhs( dzt )
     .  result( lhs )

!       Description:
!       Turbulent advection of xm:  implicit portion of the code.
!
!       The d(xm)/dt equation contains a turbulent advection term:
!
!       - d(w'x')/dz.
!
!       This term is solved for completely implicitly, such that:
!
!       - d( w'x'(t+1) )/dz.
!
!       Note:  When the term is brought over to the left-hand side, the
!              sign is reversed and the leading "-" in front of the
!              term is changed to a "+".
!
!       The timestep index (t+1) means that the value of w'x' being used
!       is from the next timestep, which is being advanced to in solving
!       the d(xm)/dt and d(w'x')/dt equations.
!
!       This term is discretized as follows:
!
!       While the values of xm are found on the thermodynamic levels, the
!       values of w'x' are found on the momentum levels.  The
!       derivative of w'x' is taken over the intermediate (central)
!       thermodynamic level, yielding the desired results.
!
!       ===================wpxp================================== m(k)
!
!       -----------------------------d(wpxp)/dz------------------ t(k)
!
!       ===================wpxpm1================================ m(k-1)
!
!       The vertical indices m(k), t(k), and m(k-1) correspond with 
!       altitudes zm(k), zt(k), and zm(k-1), respectively.  The letter 
!       "t" is used for thermodynamic levels and the letter "m" is used 
!       for momentum levels.
!
!       dzt(k) = 1 / ( zm(k) - zm(k-1) )

!       References:
!-----------------------------------------------------------------------

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  k_mdiag   = 1,   ! Momentum superdiagonal index.
     .  km1_mdiag = 2    ! Momentum subdiagonal index.

        ! Input Variables
        real, intent(in) ::
     .  dzt   ! Inverse of grid spacing (k)   [1/m]

        ! Return Variable
        real, dimension(2) :: lhs

        ! Momentum superdiagonal [ x wpxp(k,<t+1>) ]
        lhs(k_mdiag)
     .  = + dzt

        ! Momentum subdiagonal [ x wpxp(k-1,<t+1>) ]
        lhs(km1_mdiag)
     .  = - dzt

        return
        end function xm_term_ta_lhs

!===============================================================================
        pure function wpxp_term_ta_lhs( wp2_ztp1, wp2_zt, 
     .                                  a1_ztp1, a1_zt,
     .                                  wp3p1, wp3, dzm, 
     .                                  wtol_sqd, level )
     .  result( lhs )

!       Description:
!       Turbulent advection of w'x':  implicit portion of the code.
!
!       The d(w'x')/dt equation contains a turbulent advection term:
!
!       - d(w'^2x')/dz.
!
!       A substitution is made in order to close the turbulent advection
!       term, such that:
!
!       w'^2x' = a_1 * ( w'^3 / w'^2 ) * w'x',
!
!       where a_1 is a variable that is a function of Sc.  The turbulent
!       advection term becomes:
!
!       - d [ a_1 * ( w'^3 / w'^2 ) * w'x' ] / dz.
!
!       This term is solved for completely implicitly, such that:
!
!       - d [ a_1 * ( w'^3 / w'^2 ) * w'x'(t+1) ] / dz.
!
!       Note:  When the term is brought over to the left-hand side, the
!              sign is reversed and the leading "-" in front of the
!              term is changed to a "+".
!
!       The timestep index (t+1) means that the value of w'x' being used
!       is from the next timestep, which is being advanced to in solving
!       the d(w'x')/dt equation.
!
!       This term is discretized as follows:
!
!       The values of w'x', w'^2, and a_1 are found on the momentum 
!       levels, while the values of w'^3 are found on the thermodynamic 
!       levels.  Each of the variables w'x', w'^2, and a_1 are 
!       interpolated to the intermediate thermodynamic levels.  The 
!       values of the mathematical expression (called F here) within the
!       dF/dz term are computed on the thermodynamic levels.  Then, the
!       derivative (d/dz) of the expression (F) is taken over the 
!       central momentum level, yielding the desired result.
!       In this function, the values of F are as follows:
!
!       F = a_1(t) * ( w'^3(t) / w'^2(t) ) * w'x'(t+1),
!
!       where the timestep index (t) stands for the index of the current
!       timestep.
!
!
!       =a1p1========wp2p1========wpxpp1========================= m(k+1)
!
!       -----a1(interp)---wp2(interp)---wpxp(interp)---wp3p1----- t(k+1)
!
!       =a1==========wp2==========wpxp===================dF/dz=== m(k)
!
!       -----a1(interp)---wp2(interp)---wpxp(interp)---wp3------- t(k)
!
!       =a1m1========wp2m1========wpxpm1========================= m(k-1)
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

        use grid_class, only:
     .      gr ! Variable

        implicit none

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
     .  wp2_ztp1,   ! w'^2 interpolated to thermodynamic level (k+1) [m^2/s^2]
     .  wp2_zt,     ! w'^2 interpolated to thermodynamic level (k)   [m^2/s^2]
     .  a1_ztp1,    ! a_1 interpolated to thermodynamic level (k+1)  [-]
     .  a1_zt,      ! a_1 interpolated to thermodynamic level (k)    [-]
     .  wp3p1,      ! w'^3(k+1)                                      [m^3/s^3]
     .  wp3,        ! w'^3(k)                                        [m^3/s^3]
     .  dzm,        ! Inverse of grid spacing (k)                    [1/m]
     .  wtol_sqd    ! w wind component tolerance squared             [m^2/s^2]

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

        ! Note:  The w'x' turbulent advection term, which is
        !        - d [ a_1 * ( w'^3 / w'^2 ) * w'x' ] / dz, still keeps 
        !        the a_1 term inside the derivative, unlike the w'^3 
        !        equation (in wp23.F) and the equations in diag_var.F 
        !        for r_t'^2, th_l'^2, r_t'th_l', u'^2, v'^2, sclr'r_t', 
        !        sclr'th_l', and sclr'^2.  Brian.

        ! Momentum superdiagonal: [ x wpxp(k+1,<t+1>) ]
        lhs(kp1_mdiag)
     .  = + dzm
     .      * a1_ztp1 * ( wp3p1 / max( wp2_ztp1, wtol_sqd ) )
     .      * gr%weights_zm2zt(m_above,tkp1)

        ! Momentum main diagonal: [ x wpxp(k,<t+1>) ]
        lhs(k_mdiag)
     .  = + dzm
     .      * (   a1_ztp1 * ( wp3p1 / max( wp2_ztp1, wtol_sqd ) )
     .            * gr%weights_zm2zt(m_below,tkp1)
     .          - a1_zt * ( wp3 / max( wp2_zt, wtol_sqd ) )
     .            * gr%weights_zm2zt(m_above,tk)
     .        )

        ! Momentum subdiagonal: [ x wpxp(k-1,<t+1>) ]
        lhs(km1_mdiag)
     .  = - dzm
     .      * a1_zt * ( wp3 / max( wp2_zt, wtol_sqd ) )
     .      * gr%weights_zm2zt(m_below,tk)

        return
        end function wpxp_term_ta_lhs

!===============================================================================
        pure function wpxp_term_tp_lhs( wp2, dzm )
     .  result( lhs )

!       Description:
!       Turbulent production of w'x':  implicit portion of the code.
!
!       The d(w'x')/dt equation contains a turbulent production term:
!
!       - w'^2 d(xm)/dz.
!
!       This term is solved for completely implicitly, such that:
!
!       - w'^2 * d( xm(t+1) )/dz.
!
!       Note:  When the term is brought over to the left-hand side, the
!              sign is reversed and the leading "-" in front of the
!              term is changed to a "+".
!
!       The timestep index (t+1) means that the value of xm being used
!       is from the next timestep, which is being advanced to in solving
!       the d(w'x')/dt and d(xm)/dt equations.
!
!       This term is discretized as follows:
!
!       The values of xm are found on thermodynamic levels, while the
!       values of w'^2 are found on momentum levels.  The derivative of
!       xm is taken over the intermediate (central) momentum level, where
!       it is multiplied by w'^2, yielding the desired result.
!
!       ---------------------------xmp1-------------------------- t(k+1)
!
!       ==========wp2=====================d(xm)/dz=============== m(k)
!
!       ---------------------------xm---------------------------- t(k)
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

        ! Constant parameters
        integer, parameter ::
     .  kp1_tdiag = 1,   ! Thermodynamic superdiagonal index.
     .  k_tdiag = 2      ! Thermodynamic subdiagonal index.

        ! Input Variables
        real, intent(in) ::
     .  wp2,   ! w'^2(k)                       [m^2/s^2]
     .  dzm    ! Inverse of grid spacing (k)   [1/m]

        ! Return Variable
        real, dimension(2) :: lhs

        ! Thermodynamic superdiagonal [ x xm(k+1,<t+1>) ]
        lhs(kp1_tdiag)
     .  = + wp2 * dzm

        ! Thermodynamic subdiagonal [ x xm(k,<t+1>) ]
        lhs(k_tdiag)
     .  = - wp2 * dzm

        return
        end function wpxp_term_tp_lhs

!===============================================================================
        pure function wpxp_terms_ac_pr2_lhs( C7_Skw_fnc, 
     .                                       wmtp1, wmt, dzm )
     .  result( lhs )

!       Description:
!       Accumulation of w'x' and w'x' pressure term 2:
!       implicit portion of the code.
!
!       The d(w'x')/dt equation contains an accumulation term:
!
!       - w'x' dw/dz;
!
!       and pressure term 2:
!
!       + C_7 w'x' dw/dz.
!
!       Both the w'x' accumulation term and pressure term 2 are 
!       completely implicit.  The accumulation term and pressure term 2
!       are combined and solved together as:
!
!       - ( 1 - C_7 ) * w'x'(t+1) * dw/dz.
!
!       Note:  When the term is brought over to the left-hand side, the
!              sign is reversed and the leading "-" in front of the 
!              term is changed to a "+".
!
!       The timestep index (t+1) means that the value of w'x' being used
!       is from the next timestep, which is being advanced to in solving
!       the d(w'x')/dt equation.
!
!       The terms are discretized as follows:
!
!       The values of w'x' are found on momentum levels, while the values
!       of wmt (mean vertical velocity on thermodynamic levels) are found
!       on thermodynamic levels.  The vertical derivative of wmt is
!       taken over the intermediate (central) momentum level.  It is then
!       multiplied by w'x' (implicitly calculated at timestep (t+1)) and
!       the coefficients to yield the desired results.
!
!       -------wmtp1--------------------------------------------- t(k+1)
!
!       ===============d(wmt)/dz============wpxp================= m(k)
!
!       -------wmt----------------------------------------------- t(k)
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

        ! Constant parameters
        integer, parameter ::
     .  k_mdiag   = 1    ! Momentum main diagonal index.

        ! Input Variables
        real, intent(in) ::
     .  C7_Skw_fnc, ! C_7 parameter with Sk_w applied (k)   [-]
     .  wmtp1,      ! wmt(k+1)                              [m/s]
     .  wmt,        ! wmt(k)                                [m/s]
     .  dzm         ! Inverse of grid spacing (k)           [1/m]

        ! Return Variable
        real :: lhs

        ! Momentum main diagonal: [ x wpxp(k,<t+1>) ]
        lhs
     .  = + ( 1.0 - C7_Skw_fnc ) * dzm * ( wmtp1 - wmt )

        return
        end function wpxp_terms_ac_pr2_lhs

!===============================================================================
        pure function wpxp_term_pr1_lhs( C6x_Skw_fnc, taum )
     .  result( lhs )

!       Description
!       Pressure term 1 for w'x':  implicit portion of the code.
!
!       The d(w'x')/dt equation contains pressure term 1:
!
!       - ( C_6 / tau_m ) w'x'.
!
!       This term is solved for completely implicitly, such that:
!
!       - ( C_6 / tau_m ) w'x'(t+1)
!
!       Note:  When the term is brought over to the left-hand side, the
!              sign is reversed and the leading "-" in front of the
!              term is changed to a "+".
!
!       The timestep index (t+1) means that the value of w'x' being used
!       is from the next timestep, which is being advanced to in solving
!       the d(w'x')/dt equation.
!
!       The values of w'x' are found on the momentum levels.  The values 
!       of the C_6 skewness function and time-scale tau_m are also found 
!       on the momentum levels.

!       References:
!-----------------------------------------------------------------------

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  k_mdiag   = 1    ! Momentum main diagonal index.

        ! Input Variables
        real, intent(in) ::
     .  C6x_Skw_fnc, ! C_6x parameter with Sk_w applied (k)   [-]
     .  taum         ! Time-scale tau at momentum levels (k)  [s]

        ! Return Variable
        real :: lhs

        ! Momentum main diagonal: [ x wpxp(k,<t+1>) ]
        lhs
     .  = + C6x_Skw_fnc / taum

        return
        end function wpxp_term_pr1_lhs

!===============================================================================
        pure function wpxp_terms_bp_pr3_rhs( C7_Skw_fnc, xpthvp )
     .  result( rhs )

!       Description:
!       Buoyancy production of w'x' and w'x' pressure term 3:
!       explicit portion of the code.
!
!       The d(w'x')/dt equation contains a buoyancy production term:
!
!       + (g/th_0) x'th_v';
!
!       and pressure term 3:
!
!       - C_7 (g/th_0) x'th_v'.
!
!       Both the w'x' buoyancy production term and pressure term 3 are
!       completely explicit.  The buoyancy production term and pressure 
!       term 3 are combined and solved together as:
!
!       + ( 1 - C_7 ) * (g/th_0) * x'th_v'.

!       References:
!-----------------------------------------------------------------------

        use constants, only:
            ! Variable(s) 
     .      grav ! Gravitational acceleration [m/s^2]
        use parameters, only:
            ! Variable(s) 
     .      T0  ! Reference temperature      [K]

        implicit none

        ! Input Variables
        real, intent(in) ::
     .  C7_Skw_fnc, ! C_7 parameter with Sk_w applied (k)   [-]
     .  xpthvp      ! x'th_v'(k)                            [K {xm units}]

        ! Return Variable
        real :: rhs

        rhs
     .  = grav/T0 * ( 1.0 - C7_Skw_fnc ) * xpthvp

        return
        end function wpxp_terms_bp_pr3_rhs

!===============================================================================

        end module mixing
