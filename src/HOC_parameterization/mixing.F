!-----------------------------------------------------------------------
! $Id: mixing.F,v 1.46 2008-04-19 15:55:42 griffinb Exp $

        ! Passive scalars
#define SCLR_THETA 1
#define SCLR_RT 2

!-----------------------------------------------------------------------
        module mixing

!       Description:
!       Contains the HOC mixing scheme.  Allows for optional variables.

!       References:
!       None
!-----------------------------------------------------------------------

        implicit none

        private
        ! Parameter Constants
        integer, parameter ::
     .  nsub = 2,  ! Number of subdiagonals in the LHS matrix
     .  nsup = 2,  ! Number of superdiagonals in the LHS matrix
     .  nrhs = 1   ! Number of RHS vectors

        public  :: timestep_mixing

        private :: mixing_solve, 
     .             mixing_lhs,
     .             wpxp_term_ma_lhs,
     .             wpxp_term_ta_lhs,
     .             wpxp_term_tp_lhs,
     .             wpxp_terms_ac_pr2_lhs,
     .             wpxp_term_pr1_lhs

        contains
!-----------------------------------------------------------------------
        subroutine timestep_mixing
     .             ( dt, taum, wmm, wmt, wp2, wp3, Scm, Skwm, Kht,
     .               rtm, wprtp, rtpthvp, rtm_forcing,
     .               thlm, wpthlp, thlpthvp, thlm_forcing,
     .               rtp2, thlp2, ! Added for clipping Vince Larson 20070929
     .               rhot, rhom, exner, 
     .               err_code, implemented,
     .               sclrm, wpsclrp, sclrpthvp, sclrm_forcing,
     .               sclrp2 )     ! Added for clipping Vince Larson 20070929
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
        use parameters
        use constants
        use grid_class
        use stats_prec, only: time_prec
        use error_code, only: CLUBB_singular_matrix,
     .                        CLUBB_bad_lapack_arg
     
#ifdef STATS
        use stats_hoc
#endif /*STATS*/
        
        implicit none

        ! External
        intrinsic :: exp, present

        ! Input Variables
        real(kind=time_prec), intent(in) :: 
     .  dt           ! Timestep                       [s]

        real, intent(in), dimension(gr%nnzp) ::
     .  taum,        ! Tau on moment. grid                [s]
     .  wmm,         ! w wind on m                        [m/s]
     .  wmt,         ! w wind on t                        [m/s]
     .  wp2,         ! w'^2                               [m^2/s^2]
     .  wp3,         ! w'^3                               [m^3/s^3]
     .  Scm,         ! Sc on moment. grid                 [-]
     .  Skwm,        ! Skw on moment. grid                [-]
     .  Kht,         ! Coefficient of Eddy Diffusivity    [m^2/s]
     .  rtpthvp,     ! r_t'th_v'                          [(kg K)/kg]
     .  rtm_forcing, ! rt forcing                         [(kg/kg)/s]
     .  thlpthvp,    ! th_l' th_v'                        [K^2]
     .  thlm_forcing, ! Large scale tendency for thlm      [K/s]
! Added for clipping by Vince Larson 29 Sep 2007
     .  rtp2,        ! Variance of total water mix rat    [(kg/kg)^2]
     .  thlp2       ! Variance of theta_l (liq pot temp) [K^2]
! End of Vince Larson's addition.

        ! Additions for DYCOMS 2 RF02
        real, intent(in), dimension(gr%nnzp) ::
     .  rhot,        ! rho (density) on thermo. grid  [kg/m^3] 
     .  rhom,        ! rho (density) on moment. grid  [kg/m^3] 
     .  exner        ! Exner function                 [-]

        logical, intent(in) :: 
     .  implemented  ! Flag for if HOC is implemented in a larger model or not.
         
        ! Input/Output Variables
        real, intent(inout), dimension(gr%nnzp) :: 
     .  rtm,         ! r_t Total water mixing ratio   [kg/kg]
     .  wprtp,       ! w' r_t'                        [(m kg)/(s kg)]
     .  thlm,        ! th_l Liquid potential temp.    [K]
     .  wpthlp       ! w' th_l'                       [(m K)/s]

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
        real, dimension(2+2+1,2*gr%nnzp) ::
     .  lhs  ! Band diagonal sparse matrix stored in compact form 

        real, dimension(2*gr%nnzp,nrhs) ::
     .  rhs ! RHS vector for banded matrix

        ! Constant parameters as a function of Skw.
        real, dimension(gr%nnzp) :: 
     .  C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc

        ! Misc. storage
!       real, dimension(gr%nnzp) :: tmp1, tmp2

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

        end do

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

!        call mixing_lhs( .true., dt, taum, wmm, wmt, Kw6, 
        call mixing_lhs( .true., dt, taum, wmm, wmt, Kw6_rt, 
     .                   wp2, wp3, Scm, C6rt_Skw_fnc, C7_Skw_fnc,
     .                   lhs, implemented )

        ! Solve rtm / wprtp
        call mixing_solve( "rtm", .true., dt, taum,
     .                     rtm, wprtp, rtpthvp, rtm_forcing,
     .                     wp2, rtp2, ! For clipping Vince Larson 20070929
     .                     rhot, rhom, exner, 
     .                     lhs, rhs, C7_Skw_fnc, err_code )  

        if ( err_code == CLUBB_singular_matrix .or.
     .       err_code == CLUBB_bad_lapack_arg )  then
     
          write(fstderr,'(a)') "rt mixing failed"
          return
        end if
 
        ! Clipping rtm
        ! Joshua Fasching March 2008
        
#ifdef STATS
        ! Computed value before clipping
        if (lstats_samp .and. irtm_cl > 0) then
          zt%x(:,irtm_cl) = zt%x(:,irtm_cl) - rtm / dt
        endif 
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
            write(fstderr,*) "rtm < 0 in mixing at k= ", k 
          end if
        end do
                
#ifdef STATS
        if (lstats_samp .and. irtm_cl > 0) then
          zt%x(:,irtm_cl) = zt%x(:,irtm_cl) + rtm /dt
          zt%n(:,irtm_cl) = zt%n(:,irtm_cl) + 1
        endif
#endif
        
        ! Prof. Larson's modified mixing_lhs 
!        call mixing_lhs( .true., dt, taum, wmm, wmt, Kw6,
        call mixing_lhs( .true., dt, taum, wmm, wmt, Kw6_thl,
     .                   wp2, wp3, Scm, C6thl_Skw_fnc, C7_Skw_fnc,
     .                   lhs, implemented )


        ! Solve for thlm / wpthlp
        call mixing_solve( "thlm", .true., dt, taum,
     .                     thlm, wpthlp, thlpthvp, thlm_forcing,
     .                     wp2, thlp2, ! For clipping Vince Larson 20070929
     .                     rhot, rhom, exner, 
     .                     lhs, rhs, C7_Skw_fnc, err_code )


        if ( err_code == CLUBB_singular_matrix .or.
     .       err_code == CLUBB_bad_lapack_arg ) then       
          write(fstderr,'(a)') "thetal mixing failed"
          return
        end if
        
        ! Clipping thlm
        ! Joshua Fasching March 2008
        
#ifdef STATS
        ! Computed value before clipping
        if (lstats_samp .and. ithlm_cl > 0) then
          zt%x(:,ithlm_cl) = zt%x(:,ithlm_cl) - thlm / dt
        endif 
#endif        
        ! The value of potential temperature cannot fall below 0,
        ! so we clip accordingly
        do k = 1, gr%nnzp, 1 
          if ( thlm(k) < 0.0 ) then
            thlm(k) = 0.0
            write(fstderr,*) "thlm < 0 in mixing at k= ", k 
          end if
        end do
        
                
#ifdef STATS
        if (lstats_samp .and. ithlm_cl > 0) then
          zt%x(:,ithlm_cl) = zt%x(:,ithlm_cl) + thlm /dt
          zt%n(:,ithlm_cl) = zt%n(:,ithlm_cl) + 1
        endif
#endif
        ! End change Joshua Fasching March 2008
        
        ! Solve sclrm / wpsclrp
        if ( present( sclrm )     .and. present( wpsclrp ) .and.
     .       present( sclrpthvp ) .and. present( sclrm_forcing ) ) then

          call mixing_lhs( .true., dt, taum, wmm, wmt, Kw6,
     .                     wp2, wp3, Scm, C6rt_Skw_fnc, C7_Skw_fnc,
     .                     lhs, implemented )

          call mixing_solve( "sclrm_rtm", .true., dt, taum,
     .                       sclrm(:,SCLR_RT), wpsclrp(:,SCLR_RT),
     .                       sclrpthvp(:,SCLR_RT), 
     .                       sclrm_forcing(:,SCLR_RT),
     .                       wp2, sclrp2(:,SCLR_RT), ! For clipping Vince Larson 20070929
     .                       rhot, rhom, exner, 
     .                       lhs, rhs, C7_Skw_fnc, err_code )

          if ( err_code == CLUBB_singular_matrix .or.
     .         err_code == CLUBB_bad_lapack_arg ) then

    
            write(fstderr,'(a)') "Passive scalar 1 mixing failed"
            return
          end if

          call mixing_lhs( .true., dt, taum, wmm, wmt, Kw6,
     .                     wp2, wp3, Scm, C6thl_Skw_fnc, C7_Skw_fnc,
     .                     lhs, implemented )

          call mixing_solve( "sclrm_theta", .true., dt, taum,
     .                       sclrm(:,SCLR_THETA), wpsclrp(:,SCLR_THETA),
     .                       sclrpthvp(:,SCLR_THETA), 
     .                       sclrm_forcing(:,SCLR_THETA),
     .                       wp2, sclrp2(:,SCLR_THETA), ! For clipping Vince Larson 20070929
     .                       rhot, rhom, exner, 
     .                       lhs, rhs, C7_Skw_fnc, err_code )

          if ( err_code == CLUBB_singular_matrix .or.
     .       err_code == CLUBB_bad_lapack_arg ) then  
            write(fstderr,'(a)') "Passive scalar 2 mixing failed"
            return
          end if
        end if ! passive scalars
        
!       Error Report
!       Joshua Fasching Feb 2008 

        if ( err_code == CLUBB_singular_matrix .or.
     .       err_code == CLUBB_bad_lapack_arg ) then
                
           write(fstderr,*) "Error in timestep_mixing"
           
           write(fstderr,*) "Intent(in)"
           
           write(fstderr,*) "dt = ", dt
           write(fstderr,*) "taum = ", taum
           write(fstderr,*) "wmm = ", wmm
           write(fstderr,*) "wmt = ", wmt
           write(fstderr,*) "wp2 = ", wp2
           write(fstderr,*) "wp3 = ", wp3
           write(fstderr,*) "Scm = ", Scm
           write(fstderr,*) "Skwm = ", Skwm
           write(fstderr,*) "Kht = ", Kht
           write(fstderr,*) "rtpthvp = ", rtpthvp
           write(fstderr,*) "rtm_forcing = ", rtm_forcing
           write(fstderr,*) "thlpthvp = ", thlpthvp
           write(fstderr,*) "thlm_forcing = ", thlm_forcing
           write(fstderr,*) "rtp2 = ", rtp2
           write(fstderr,*) "thlp2 = ", thlp2
           write(fstderr,*) "rhot = ", rhot
           write(fstderr,*) "rhom = ", rhom
           write(fstderr,*) "exner = ", exner
           
           if( present( sclrpthvp ) ) then
              write(fstderr,*) "sclrpthvp = ", sclrpthvp
           endif
           
           if( present( sclrm_forcing ) ) then
              write(fstderr,*) "sclrm_forcing = ", sclrm_forcing
           endif
           
           if( present( sclrp2 ) ) then
              write(fstderr,*) "sclrp2 = ", sclrp2
           endif
           
           write(fstderr,*) "Intent(inout)"
           
           write(fstderr,*) "rtm = ", rtm 
           write(fstderr,*) "wprtp = ", wprtp
           write(fstderr,*) "thlm = ", thlm
           write(fstderr,*) "wpthlp =", wpthlp
           
           if( present( sclrm ) ) then
              write(fstderr,*) "sclrm = ", sclrm
           endif
           
           if( present( wpsclrp ) ) then
              write(fstderr,*) "wpsclrp = ", wpsclrp
           endif
           
        end if
 
        return
        
        end subroutine timestep_mixing

!------------------------------------------------------------------------
        subroutine mixing_lhs( liter, dt, taum, wmm, wmt, Kw6, 
     .                         wp2, wp3, Scm, C6x_Skw_fnc, C7_Skw_fnc,
     .                         lhs, implemented )
!       Description:
!       References:
!------------------------------------------------------------------------

        use parameters
        use constants
        use grid_class
        use stats_prec, only: time_prec
        use diffusion, only: diffusion_mlev
#ifdef STATS
        use stats_hoc
#endif /*STATS*/

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  kp1diag = 1,
     .  kdiag   = 2,
     .  km1diag = 3

        ! External
        intrinsic :: min, max

        ! Input variables
        logical, intent(in) :: liter

        real(kind=time_prec), intent(in) :: dt ! Timestep     [s]

        real, intent(in), dimension(gr%nnzp) :: 
     .  taum,        ! Tau on moment. grid            [s]
     .  wmm,         ! w wind on m                    [m/s]
     .  wmt,         ! w wind on t                    [m/s]
     .  Kw6,         ! Coef. of Eddy Diffusivity      [m^2/s]
     .  wp2,         ! w'^2                           [m^2/s^2]
     .  wp3,         ! w'^3                           [m^3/s^3]
     .  Scm          ! Sc on moment. grid             [-]

        logical, intent(in)    :: implemented

        real, intent(in), dimension(gr%nnzp) :: 
     .  C6x_Skw_fnc,
     .  C7_Skw_fnc

        ! Output variables
        real, intent(out), dimension(2+2+1,2*gr%nnzp) :: lhs ! LHS matrix

        ! Local Variables
        real, dimension(gr%nnzp) :: a1

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

        ! Define a1

        a1(1:gr%nnzp) = 1.0 / ( 1.0 - Scm(1:gr%nnzp) )

        ! Set up left hand side matrix

        lhs = 0.0

        do k=2,gr%nnzp-1

          ! Define indices

          km1 = max( k-1, 1 )
          kp1 = min( k+1, gr%nnzp )

          k_xm   = 2*k - 1
          k_wpxp = 2*k

          ! xm: mean vertical advection

          if ( .not. implemented ) then

            ! Thermodynamic subdiagonal [ x xm(k-1,<t+1>) ]
            lhs(3+2,k_xm) = - wmt(k) * (1.0/2.0) * gr%dzt(k)
            ! Thermodynamic superdiagonal [ x xm(k+1,<t+1>) ]
            lhs(3-2,k_xm) = + wmt(k) * (1.0/2.0) * gr%dzt(k)

          else

            ! Thermodynamic subdiagonal [ x xm(k-1,<t+1>) ]
            lhs(3+2,k_xm) = 0.0
            ! Thermodynamic superdiagonal [ x xm(k+1,<t+1>) ]
            lhs(3-2,k_xm) = 0.0

          end if

!         xm: turbulent transport

          ! Momentum subdiagonal [ x wpxp(k-1,<t+1>) ]
          lhs(3+1,k_xm) = - gr%dzt(k)
          ! Momentum superdiagonal [ x wpxp(k,<t+1>) ]
          lhs(3-1,k_xm) = + gr%dzt(k)

!         xm: time tendency

          ! Thermodynamic main diagonal [ x xm(k,<t+1>) ]
          lhs(3,k_xm) = 1.0 / dt

#ifdef STATS
          if (lstats_samp) then

          ! Implicit contributions for rtm or wprtp
          if ( irtm_ma > 0 .or. ithlm_ma > 0 ) then
            ztscr01(k) = + wmt(k) * (1.0/2.0) * gr%dzt(k)
            ztscr02(k) = - wmt(k) * (1.0/2.0) * gr%dzt(k)
          endif ! ma

          if ( irtm_ta > 0 .or. ithlm_ta > 0 ) then
            ztscr03(k) = + gr%dzt(k)
            ztscr04(k) = - gr%dzt(k)
          endif ! ta

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
          lhs((/3-2,3+2/),k_wpxp)
     .    = lhs((/3-2,3+2/),k_wpxp)
     .    + wpxp_term_ma_lhs( wmm(k), gr%dzm(k) )

          ! LHS turbulent advection (ta) term.
          lhs((/3-2,3,3+2/),k_wpxp)
     .    = lhs((/3-2,3,3+2/),k_wpxp)
     .    + wpxp_term_ta_lhs( wp2(kp1), wp2(k), wp2(km1), 
     .                        a1(kp1), a1(k), a1(km1), 
     .                        wp3(kp1), wp3(k), gr%dzm(k), wtol_sqd )

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
     .    + diffusion_mlev( Kw6(k), Kw6(kp1), nu6,
     .                      gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )

          ! LHS time tendency.
          if (liter) lhs(3,k_wpxp) = lhs(3,k_wpxp) + 1.0 / dt

#ifdef STATS
          if (lstats_samp) then

          ! Statistics: implicit contributions for wprtp or wpthlp.

          if ( iwprtp_ma > 0 .or. iwpthlp_ma > 0 ) then
            tmp(1:2) =
     .      + wpxp_term_ma_lhs( wmm(k), gr%dzm(k) )
            zmscr01(k) = - tmp(2)
            zmscr02(k) = - tmp(1)
          endif ! ma

          if ( iwprtp_ta > 0 .or. iwpthlp_ta > 0 ) then
            tmp(1:3) =
     .      + wpxp_term_ta_lhs( wp2(kp1), wp2(k), wp2(km1),
     .                          a1(kp1), a1(k), a1(km1),
     .                          wp3(kp1), wp3(k), gr%dzm(k), wtol_sqd )
            zmscr03(k) = - tmp(3)
            zmscr04(k) = - tmp(2)
            zmscr05(k) = - tmp(1)
          endif ! ta

          if ( iwprtp_tp > 0 .or. iwpthlp_tp > 0 ) then
            tmp(1:2) =
     .      + wpxp_term_tp_lhs( wp2(k), gr%dzm(k) )
            zmscr06(k) = - tmp(2)
            zmscr07(k) = - tmp(1)
          endif

          if ( iwprtp_ac > 0 .or. iwpthlp_ac > 0 ) then
            zmscr08(k) = 
     .      - wpxp_terms_ac_pr2_lhs( 0.0,
     .                               wmt(kp1), wmt(k), gr%dzm(k) )
          endif

          if ( iwprtp_pr1 > 0 .or. iwpthlp_pr1 > 0 ) then
            zmscr09(k) =
     .      - wpxp_term_pr1_lhs( C6x_Skw_fnc(k), taum(k) )
          endif ! pr1

          if ( iwprtp_pr2 > 0 .or. iwpthlp_pr2 > 0 ) then
            zmscr10(k) =
     .      - wpxp_terms_ac_pr2_lhs( (1.0+C7_Skw_fnc(k)),
     .                               wmt(kp1), wmt(k), gr%dzm(k) )
          endif

          if ( iwprtp_dp1 > 0 .or. iwpthlp_dp1 > 0 ) then
            tmp(1:3) =
     .      + diffusion_mlev( Kw6(k), Kw6(kp1), nu6,
     .                        gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )
            zmscr11(k) = - tmp(3)
            zmscr12(k) = - tmp(2)
            zmscr13(k) = - tmp(1)
          endif

          endif
#endif /*STATS*/

        enddo ! 2..gr%nnzp-1


        ! Boundary condition for the mean and flux

        ! Fixed-point boundary conditions are used for the turbulent fluxes.

        ! Lower boundary
        k      = 1
        k_xm   = 2*k - 1
        k_wpxp = 2*k
        lhs(:,k_xm)     = 0.0
        lhs(3,k_xm)     = 1.0
        lhs(:,k_wpxp)   = 0.0
        lhs(3,k_wpxp)   = 1.0

        ! Upper boundary
        k      = gr%nnzp
        k_xm   = 2*k - 1
        k_wpxp = 2*k
        lhs(:,k_xm)     = 0.0
        lhs(3,k_xm)     = 1.0
        lhs(:,k_wpxp)   = 0.0
        lhs(3,k_wpxp)   = 1.0

        return
        end subroutine mixing_lhs

!------------------------------------------------------------------------
        
        subroutine mixing_solve
     .             ( solve_type, liter, dt, taum,
     .               xm, wpxp, xpthvp, xm_forcing,
     .               wp2, xp2, ! For clipping Vince Larson 20070929
     .               rhot, rhom, exner, 
     .               lhs, rhs, C7_Skw_fnc, err_code )


        use parameters
        use constants
        use grid_class
        use stats_prec, only: time_prec
        use lapack_wrap, only: band_solve, band_solvex
        use pos_definite_mod, only: pos_definite_adj
        use error_code, only: CLUBB_no_error
#ifdef STATS
        use stats_hoc
#endif /*STATS*/

        implicit none


        ! Input
        character(len=*), intent(in) :: 
     .  solve_type ! Variables being solved for

        logical, intent(in) :: 
     .  liter ! Whether the mean and flux terms are prognostic

        real(kind=time_prec), intent(in) :: dt      ! Timestep     [s]

        real, intent(in), dimension(gr%nnzp) :: 
     .  taum,        ! Tau on moment. grid                    [s]
     .  rhot,        ! Rho (density) on thermo. grid -Brian   [kg/m^3]
     .  rhom,        ! Rho (density) on moment. grid -Brian   [kg/m^3]
     .  exner,       ! Exner function -Brian                  [-]
     .  xpthvp,      ! x' theta_v'                            [(units)(K)]
     .  xm_forcing,  ! xm forcing -Brian                      [(units)/s] 
! Vince Larson added for clipping 20070929
     .  wp2,         ! Variance of vertical velocity          [(m/s)^2]
     .  xp2,         ! Variance of scalar (rt, thl, or sclr)  Units vary
! End Vince Larson's addition
     .  C7_Skw_fnc   ! C7 parameter (Skw dependent)           [-] 

        ! Input/Output
        real, intent(inout), dimension(gr%nnzp) :: 
     .  xm,  ! Mean term                                  [units vary] 
     .  wpxp ! Flux term                                  [units vary] 

        real, target, intent(inout), dimension(nsup+nsub+1,2*gr%nnzp) ::
     .  lhs  ! Left hand side (LAPACK)

        real, intent(inout), dimension(2*gr%nnzp,nrhs) :: 
     .  rhs  ! Right hand side (LAPACK)

        ! Output
  
        integer, intent(out) :: err_code

        ! Local Variables
        real, dimension(2*gr%nnzp,nrhs) :: 
     .  solution ! Solution to band diagonal system

        real, dimension(3,gr%nnzp) ::
     .  flux_lhs

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
     .    ixm_f,
     .    ixm_mc,
     .    ixm_cn,
     .    ixm_pd,
     .    iwpxp_bt,
     .    iwpxp_ma,
     .    iwpxp_ta,
     .    iwpxp_tp,
     .    iwpxp_ac,
     .    iwpxp_bp,
     .    iwpxp_pr1,
     .    iwpxp_pr2,
     .    iwpxp_pr3,
     .    iwpxp_dp1,
     .    iwpxp_cl,
     .    iwpxp_pd

        select case ( trim( solve_type ) )
        case ( "rtm" ) ! rtm/wprtp budget terms
          ixm_bt    = irtm_bt
          ixm_ta    = irtm_ta
          ixm_ma    = irtm_ma
          ixm_f     = irtm_forcing
          ixm_cn    = irtm_cn ! Applies to both the mean and flux terms
          ixm_pd    = irtm_pd ! Applies to both the mean and flux terms
          iwpxp_bt  = iwprtp_bt
          iwpxp_ma  = iwprtp_ma
          iwpxp_ta  = iwprtp_ta
          iwpxp_tp  = iwprtp_tp
          iwpxp_ac  = iwprtp_ac
          iwpxp_bp  = iwprtp_bp
          iwpxp_pr1 = iwprtp_pr1
          iwpxp_pr2 = iwprtp_pr2
          iwpxp_pr3 = iwprtp_pr3
          iwpxp_dp1 = iwprtp_dp1
          iwpxp_cl  = iwprtp_cl
          iwpxp_pd  = iwprtp_pd
        case ( "thlm" ) ! thlm/wpthlp budget terms
          ixm_bt    = ithlm_bt
          ixm_ta    = ithlm_ta
          ixm_ma    = ithlm_ma
          ixm_f     = ithlm_forcing
          ixm_cn    = ithlm_cn ! Applies to both the mean and flux terms
          ixm_pd    = 0
          iwpxp_bt  = iwpthlp_bt
          iwpxp_ma  = iwpthlp_ma
          iwpxp_ta  = iwpthlp_ta
          iwpxp_tp  = iwpthlp_tp
          iwpxp_ac  = iwpthlp_ac
          iwpxp_bp  = iwpthlp_bp
          iwpxp_pr1 = iwpthlp_pr1
          iwpxp_pr2 = iwpthlp_pr2
          iwpxp_pr3 = iwpthlp_pr3
          iwpxp_dp1 = iwpthlp_dp1
          iwpxp_cl  = iwpthlp_cl
          iwpxp_pd  = 0
        case default  ! this includes the sclrm case
          ixm_bt    = 0
          ixm_ta    = 0
          ixm_ma    = 0
          ixm_f     = 0
          ixm_cn    = 0
          ixm_pd    = 0
          iwpxp_bt  = 0
          iwpxp_ma  = 0
          iwpxp_ta  = 0
          iwpxp_tp  = 0
          iwpxp_ac  = 0
          iwpxp_bp  = 0
          iwpxp_pr1 = 0
          iwpxp_pr2 = 0
          iwpxp_pr3 = 0
          iwpxp_dp1 = 0
          iwpxp_cl  = 0
          iwpxp_pd  = 0
        end select
#endif /*STATS*/

#ifdef STATS
        if ( lstats_samp ) then

          ! xm total time tendency ( 1st calculation)
          if ( ixm_bt > 0 ) then
            zt%x(:,ixm_bt) = zt%x(:,ixm_bt) - xm / dt
          end if

          ! Explicit contributions for xm (including microphysics/radition)
          if ( ixm_f > 0 ) then
            zt%x(:,ixm_f) = zt%x(:,ixm_f) + xm_forcing
            zt%n(:,ixm_f) = zt%n(:,ixm_f) + 1
          end if

          ! wpxp total time tendency (1st calculation)
          if ( iwpxp_bt > 0 ) then
            zm%x(:,iwpxp_bt) = zm%x(:,iwpxp_bt) - wpxp / dt
          end if

          ! Explicit contributions for wpxp
          if ( iwpxp_bp > 0 ) then
            zm%x(:,iwpxp_bp) = zm%x(:,iwpxp_bp)
     .      + grav/T0 * 1.0 * xpthvp(:)
            zm%n(:,iwpxp_bp) = zm%n(:,iwpxp_bp) + 1
          end if

          if ( iwpxp_pr3 > 0 ) then
            zm%x(:,iwpxp_pr3) = zm%x(:,iwpxp_pr3)
     .      - grav/T0 * C7_Skw_fnc * xpthvp(:)
            zm%n(:,iwpxp_pr3) = zm%n(:,iwpxp_pr3) + 1
          end if
        end if ! lstats_samp
#endif /*STATS*/

        ! Setup matrix for positive definite advection scheme
        if ( trim( solve_type ) == "rtm" ) then

          flux_lhs = 0.0

          do k = 2, gr%nnzp-1, 1
            k_xm   = 2*k - 1
            k_wpxp = 2*k

            flux_lhs(1,k+1) = lhs(2,k_xm)
            !flux_lhs(2,k)   = 0. !  (Already 0)
            flux_lhs(3,k-1) = lhs(4,k_xm)

            ! %%
            !lhs(3-2,k_xm) = 0.0
            !lhs(3+2,k_xm) = 0.0
          end do

        end if

!       Set up RHS vector

        rhs = 0.0


!        do k=2,gr%nnzp
        do k=2,gr%nnzp-1

!         Define indices

          k_xm   = 2*k - 1
          k_wpxp = 2*k

!         xm: time tendency and forcing

!         Brian:: rtm also includes effects on rtm of production
!                   and evaporation of rain (rrm).  3/30/05.
!                 thlm also includes effects on thlm of production
!                   and evaporation of rain (rrm).  3/31/05.

          ! Right-hand side.
          rhs(k_xm,1) 
     .    = xm(k) / dt
     .    + xm_forcing(k)

!         wpxp: time tendency and buoyancy

          ! Right-hand side.
          rhs(k_wpxp,1)
     .    = grav/T0 * (1.0-C7_Skw_fnc(k)) * xpthvp(k)

          if ( liter ) rhs(k_wpxp,1) = rhs(k_wpxp,1) + wpxp(k) / dt

        end do ! k=2..gr%nnzp-1

        ! Boundary condition for the flux and mean

        ! Fixed-point boundary conditions are used for the turbulent fluxes.

        k      = 1
        k_xm   = 2*k - 1
        k_wpxp = 2*k
        rhs(k_xm,1)   = xm(k)
        rhs(k_wpxp,1) = wpxp(k)

        k      = gr%nnzp
        k_xm   = 2*k - 1
        k_wpxp = 2*k
        rhs(k_xm,1)   = xm(k)
        rhs(k_wpxp,1) = 0.0

#ifdef STATS
        if ( lstats_samp .and. ixm_cn > 0 ) then
          ! Perform LU decomp and solve system (LAPACK with diagnostics)
          call band_solvex( solve_type, nsup, nsub, 2*gr%nnzp, nrhs,
     .                      lhs, rhs, solution, rcond, err_code )

          ! Est. of the condition number of the mean/flux LHS matrix 
          sfc%x(1,ixm_cn) = sfc%x(1,ixm_cn) + 1.0 / rcond
          sfc%n(1,ixm_cn) = sfc%n(1,ixm_cn) + 1

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
        if ( err_code /= CLUBB_no_error ) return

        ! Apply a flux limiting positive definite scheme if the solution 
        ! for the mean field is negative and we're determining total water
        if ( trim( solve_type ) == "rtm" ) then

          call pos_definite_adj( dt, flux_lhs, 
     .                           solution(:,1), xm, xm_pd, wpxp_pd )
#ifdef STATS
          if ( lstats_samp ) then
            if ( iwpxp_pd > 0 ) then
              zm%x(:,iwpxp_pd) = zm%x(:,iwpxp_pd) + wpxp_pd(1:gr%nnzp)
              zm%n(:,iwpxp_pd) = zm%n(:,iwpxp_pd) + 1
            end if
            if ( ixm_pd > 0 ) then
              zt%x(:,ixm_pd) = zt%x(:,ixm_pd) + xm_pd(1:gr%nnzp)
              zt%n(:,ixm_pd) = zt%n(:,ixm_pd) + 1
            end if
          end if
#endif
        end if
        

        ! Copy result into output arrays

        do k=1,gr%nnzp

          k_xm   = 2 * k - 1
          k_wpxp = 2 * k

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
             if ( ixm_ma > 0 ) then
               zt%x(k,ixm_ma) = zt%x(k,ixm_ma)
     .         + ztscr01(k) * xm(km1)
     .         + ztscr02(k) * xm(kp1)
               zt%n(k,ixm_ma) = zt%n(k,ixm_ma) + 1
             end if ! ma

             if ( ixm_ta > 0 ) then
               zt%x(k,ixm_ta) = zt%x(k,ixm_ta)
     .         + ztscr03(k) * wpxp(km1)
     .         + ztscr04(k) * wpxp(k)
               zt%n(k,ixm_ta) = zt%n(k,ixm_ta) + 1
             end if ! ta

             ! Finalize implicit contributions for wpxp
             if ( iwpxp_ma > 0 ) then
               zm%x(k,iwpxp_ma) = zm%x(k,iwpxp_ma)
     .         + zmscr01(k) * wpxp(km1)
     .         + zmscr02(k) * wpxp(kp1)
               zm%n(k,iwpxp_ma) = zm%n(k,iwpxp_ma) + 1
             end if

             if ( iwpxp_ta > 0 ) then
               zm%x(k,iwpxp_ta) = zm%x(k,iwpxp_ta)
     .         + zmscr03(k) * wpxp(km1)
     .         + zmscr04(k) * wpxp(k)
     .         + zmscr05(k) * wpxp(kp1)
               zm%n(k,iwpxp_ta) = zm%n(k,iwpxp_ta) + 1
             end if

             if ( iwpxp_tp > 0 ) then
               zm%x(k,iwpxp_tp) = zm%x(k,iwpxp_tp)
     .         + zmscr06(k) * xm(k)
     .         + zmscr07(k) * xm(kp1)
               zm%n(k,iwpxp_tp) = zm%n(k,iwpxp_tp) + 1
             end if

             if ( iwpxp_ac > 0 ) then
               zm%x(k,iwpxp_ac) = zm%x(k,iwpxp_ac)
     .         + zmscr08(k) * wpxp(k)
               zm%n(k,iwpxp_ac) = zm%n(k,iwpxp_ac) + 1
             end if

             if ( iwpxp_pr1 > 0 ) then
               zm%x(k,iwpxp_pr1) = zm%x(k,iwpxp_pr1)
     .         + zmscr09(k) * wpxp(k)
               zm%n(k,iwpxp_pr1) = zm%n(k,iwpxp_pr1) + 1
             end if

             if ( iwpxp_pr2 > 0 ) then
               zm%x(k,iwpxp_pr2) = zm%x(k,iwpxp_pr2)
     .         + zmscr10(k) * wpxp(k)
               zm%n(k,iwpxp_pr2) = zm%n(k,iwpxp_pr2) + 1
             end if

             if ( iwpxp_dp1 > 0 ) then
               zm%x(k,iwpxp_dp1) = zm%x(k,iwpxp_dp1)
     .         + zmscr11(k) * wpxp(km1)
     .         + zmscr12(k) * wpxp(k)
     .         + zmscr13(k) * wpxp(kp1)
               zm%n(k,iwpxp_dp1) = zm%n(k,iwpxp_dp1) + 1
             end if

           end do ! 1..gr%nnzp

          end if ! lstats_samp
#endif /*STATS*/

! Vince Larson clipped turbulent flux (wpxp) to ensure 
!     that -1 < [correlation of w and x] < 1 .  20070929
#ifdef STATS
        if ( lstats_samp .and. iwpxp_cl > 0) then
          zm%x(:,iwpxp_cl) = zm%x(:,iwpxp_cl) - wpxp / dt
        end if
#endif


        ! Clip wpxp to make sure -0.99 < correlation < 0.99
        where ( wpxp >  0.99 * sqrt( wp2 * xp2 ) )
     .          wpxp =  0.99 * sqrt( wp2 * xp2 ) 
        where ( wpxp < -0.99 * sqrt( wp2 * xp2 ) )
     .          wpxp = -0.99 * sqrt( wp2 * xp2 ) 

#ifdef STATS
        if ( lstats_samp .and. iwpxp_cl > 0 ) then
          zm%x(:,iwpxp_cl) = zm%x(:,iwpxp_cl) + wpxp / dt
          zm%n(:,iwpxp_cl) = zm%n(:,iwpxp_cl) + 1
        end if
#endif
! End of Vince Larson's clipping.

#ifdef STATS
        if ( lstats_samp ) then
          ! xm time tendency (2nd calculation)
          if ( ixm_bt > 0 ) then
            zt%x(:,ixm_bt) = zt%x(:,ixm_bt) + xm / dt
            zt%n(:,ixm_bt) = zt%n(:,ixm_bt) + 1
          end if

          ! wpxp time tendency (2nd calculation)
          if ( iwpxp_bt > 0 ) then
            zm%x(:,iwpxp_bt) = zm%x(:,iwpxp_bt) + wpxp / dt
            zm%n(:,iwpxp_bt) = zm%n(:,iwpxp_bt) + 1
          end if

        end if
#endif /*STATS*/

        return
        end subroutine mixing_solve

!===============================================================================
        pure function wpxp_term_ma_lhs( wmm, dzm )
     .  result( lhs )

!       Description:
!       Mean advection of w'x': implicit portion of the code.
!
!       The d(w'x')/dt equation contains a mean advection term:
!
!       - w * d(w'x')/dz.
!
!       This term is solved for completely implicitly, such that:
!
!       - w * d( w'x'(t+1) )/dz.
!
!       Note:  When the term is brought over to the left-hand side, the
!              sign is reversed and the leading "-" in front of the
!              term is changed to a "+".
!
!       The timestep index (t+1) means that the value of w'x' being used
!       is from the timestep being advanced to.
!
!       This term is discretized as follows:
!
!       The values of w'x' are found on the momentum levels, as are the
!       values of wmm (mean vertical velocity on momentum levels).
!       The variable w'x' is interpolated to the intermediate
!       thermodynamic levels.  The derivative of the interpolated values
!       is taken over the central momentum level.  The derivative
!       is multiplied by wmm at the central momentum level to get
!       the desired result.
!
!       ==========wpxpp1========================================= m(k+1)
!
!       -------------------wpxp(interp)-------------------------- t(k+1)
!
!       ==========wpxp====================d(wpxp)/dz=======wmm=== m(k)
!
!       -------------------wpxp(interp)-------------------------- t(k)
!
!       ==========wpxpm1========================================= m(k-1)
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

        ! Constant parameters
        integer, parameter ::
     .  kp1_mdiag = 1,   ! Momentum superdiagonal index.
     .  km1_mdiag = 2    ! Momentum subdiagonal index.

        ! Input Variables
        real, intent(in) ::
     .  wmm,  ! wmm(k)                        [m/s]
     .  dzm   ! Inverse of grid spacing (k)   [1/m]

        ! Return Variable
        real, dimension(2) :: lhs

        ! Momentum superdiagonal: [ x wpxp(k+1,<t+1>) ]
        lhs(kp1_mdiag)
     .  = + wmm * (1.0/2.0) * dzm

        ! Momentum subdiagonal: [ x wpxp(k-1,<t+1>) ]
        lhs(km1_mdiag)
     .  = - wmm * (1.0/2.0) * dzm

        return
        end function wpxp_term_ma_lhs

!===============================================================================
        pure function wpxp_term_ta_lhs( wp2p1, wp2, wp2m1, 
     .                                  a1p1, a1, a1m1, 
     .                                  wp3p1, wp3, dzm, wtol_sqd )
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
!       is from the timestep being advanced to.
!
!       This term is discretized as follows:
!
!       The values of w'x', w'^2, and a_1 are found on the momentum 
!       levels, while the values of w'^3 are found on the thermodynamic 
!       levels.  Each of the variables w'x', w'^2, and a_1 are 
!       interpolated to the intermediate thermodynamic levels.  The 
!       values of the mathematical expression (called Y here) within the
!       dY/dz term are computed on the thermodynamic levels.  Then, the
!       derivative (d/dz) of the expression (Y) is taken over the 
!       central momentum level, yielding the desired result.
!       In this function, the values of Y are as follows:
!
!       Y = a_l * ( w'^3 / w'^2 ) * w'x'(t+1).
!
!
!       =a1p1========wp2p1========wpxpp1========================= m(k+1)
!
!       -----a1(interp)---wp2(interp)---wpxp(interp)---wp3p1----- t(k+1)
!
!       =a1==========wp2==========wpxp===================dY/dz=== m(k)
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

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  kp1_mdiag = 1,   ! Momentum superdiagonal index.
     .  k_mdiag   = 2,   ! Momentum main diagonal index.
     .  km1_mdiag = 3    ! Momentum subdiagonal index.

        ! Input Variables
        real, intent(in) ::
     .  wp2p1,   ! w'^2(k+1)                     [m^2/s^2]
     .  wp2,     ! w'^2(k)                       [m^2/s^2]
     .  wp2m1,   ! w'^2(k-1)                     [m^2/s^2]
     .  a1p1,    ! a1(k+1)                       [-]
     .  a1,      ! a1(k)                         [-]
     .  a1m1,    ! a1(k-1)                       [-]
     .  wp3p1,   ! w'^3(k+1)                     [m^3/s^3]
     .  wp3,     ! w'^3(k)                       [m^3/s^3]
     .  dzm,     ! Inverse of grid spacing (k)   [1/m]
     .  wtol_sqd ! Model parameter wtol^2        [m^2/s^2]

        ! Return Variable
        real, dimension(3) :: lhs

        ! Momentum superdiagonal: [ x wpxp(k+1,<t+1>) ]
        lhs(kp1_mdiag)
     .  = + dzm
     .      * (1.0/2.0) * ( a1p1 + a1 )
     .      * ( wp3p1 / MAX( wp2 + wp2p1, 2.0*wtol_sqd ) )

        ! Momentum main diagonal: [ x wpxp(k,<t+1>) ]
        lhs(k_mdiag)
     .  = + dzm
     .      * (   (1.0/2.0) * ( a1p1 + a1 )
     .            * ( wp3p1 / MAX( wp2 + wp2p1, 2.0*wtol_sqd ) )
     .          -  (1.0/2.0) * ( a1 + a1m1 )
     .            * ( wp3 / MAX( wp2m1 + wp2, 2.0*wtol_sqd ) )
     .        )

        ! Momentum subdiagonal: [ x wpxp(k-1,<t+1>) ]
        lhs(km1_mdiag)
     .  = - dzm
     .      * (1.0/2.0) * ( a1 + a1m1 )
     .      * ( wp3 / MAX( wp2m1 + wp2, 2.0*wtol_sqd ) )

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
!       - w'^2 d(x)/dz.
!
!       This term is solved for completely implicitly, such that:
!
!       - w'^2 * d( x(t+1) )/dz.
!
!       Note:  When the term is brought over to the left-hand side, the
!              sign is reversed and the leading "-" in front of the
!              term is changed to a "+".
!
!       The timestep index (t+1) means that the value of x being used
!       is from the timestep being advanced to.
!
!       This term is discretized as follows:
!
!       The values of x are found on thermodynamic levels, while the
!       values of w'^2 are found on momentum levels.  The derivative of
!       x is taken over the intermediate (central) momentum level, where
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
!       - w'x' dw/dz,
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
!       is from the timestep being advanced to.
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
!       Pressure term 1 for w'x': implicit portion of the code.
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
!       is from the timestep being advanced to.
!
!       The values of w'x' are found on momentum levels.  The values of
!       the C_6 skewness function and time-scale taum are also found on
!       momentum levels.

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

        end module mixing
