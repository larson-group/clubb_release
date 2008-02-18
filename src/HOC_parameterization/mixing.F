!-----------------------------------------------------------------------
! $Id: mixing.F,v 1.31 2008-02-18 14:30:17 griffinb Exp $

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

        private :: mixing_solve, mixing_lhs

        contains
!-----------------------------------------------------------------------
        subroutine timestep_mixing
     .             ( dt, taum, wmm, wmt, wp2, wp3, Scm, Skwm, Kht,
     .               rtm, wprtp, rtpthvp, rtm_forcing,
     .               thlm, wpthlp, thlpthvp, thlm_forcing,
     .               rtp2, thlp2, ! Added for clipping Vince Larson 20070929
     .               rhot, rhom, exner, rcm, Ncm, 
     .               hm_thl_tndcy, hm_rt_tndcy,
     .               Fcsed, sed_rcm,
     .               isValid, implemented,
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

        implicit none

        ! External
        intrinsic :: exp, present

        ! Input Variables
        real, intent(in) :: 
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
     .  exner,       ! Exner function                 [-]
     .  rcm,         ! Liquid water mixing ratio      [kg/kg]
     .  Ncm,         ! Cloud droplet number conc.     [num/kg]
     .  hm_rt_tndcy, ! Change in rt due to microphys. [kg/kg/s]
     .  hm_thl_tndcy ! Change in thl due to microphys.[K/s]

        logical, intent(in) :: 
     .  implemented  ! Flag for if HOC is implemented in a larger model or not.
         
        ! Input/Output Variables
        real, intent(inout), dimension(gr%nnzp) :: 
     .  rtm,         ! r_t Total water mixing ratio   [kg/kg]
     .  wprtp,       ! w' r_t'                        [(m kg)/(s kg)]
     .  thlm,        ! th_l Liquid potential temp.    [K]
     .  wpthlp       ! w' th_l'                       [(m K)/s]

        logical, intent(inout) :: isValid ! Model status 

        ! Output Variables
        ! Addition for DYCOMS_2
        real, intent(out), dimension(gr%nnzp) :: 
     .  Fcsed,       ! Cloud water sedimentation flux       [kg/(m^2 s)]
     .  sed_rcm      ! d(rcm)/dt due to cloud sedimentation [kg/(m^2 s)]

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

        wprtp_zt  = zm2zt(wprtp)
        wpthlp_zt = zm2zt(wpthlp)

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

!        call mixing_lhs( .true., dt, taum, wmm, wmt, Kw6, 
        call mixing_lhs( .true., dt, taum, wmm, wmt, Kw6_rt, 
     .                   wp2, wp3, Scm, C6rt_Skw_fnc, C7_Skw_fnc,
     .                   lhs, implemented )

        ! Solve rtm / wprtp
        call mixing_solve( "rtm", .true., dt, taum,
     .                     rtm, wprtp, rtpthvp, rtm_forcing,
     .                     wp2, rtp2, ! For clipping Vince Larson 20070929
     .                     rhot, rhom, exner, rcm, Ncm,  ! Brian
     .                     hm_rt_tndcy,
     .                     Fcsed, sed_rcm,               ! Brian
     .                     lhs, rhs, C7_Skw_fnc, isValid )  

        if ( .not. isValid ) then
          write(fstderr,'(a)') "rt mixing failed"
          return
        end if

        ! The arm_0003 case produces negative rtm near the tropopause.
        !    To avoid this, we clip rtm.  This is not a good solution,
        !    because it renders rtm non-conserved.  We should look into 
        !    a positive definite advection scheme.
        !    Vince Larson.  13 Nov 2007
        rtm = max( rtm, 0.0 ) 
        ! End of Vince Larson's change.

        ! Prof. Larson's modified mixing_lhs 
!        call mixing_lhs( .true., dt, taum, wmm, wmt, Kw6,
        call mixing_lhs( .true., dt, taum, wmm, wmt, Kw6_thl,
     .                   wp2, wp3, Scm, C6thl_Skw_fnc, C7_Skw_fnc,
     .                   lhs, implemented )


        ! Solve for thlm / wpthlp
        call mixing_solve( "thlm", .true., dt, taum,
     .                     thlm, wpthlp, thlpthvp, thlm_forcing,
     .                     wp2, thlp2, ! For clipping Vince Larson 20070929
     .                     rhot, rhom, exner, rcm, Ncm,  ! Brian
     .                     hm_thl_tndcy,
     .                     Fcsed, sed_rcm,               ! Brian
     .                     lhs, rhs, C7_Skw_fnc, isValid )

        if ( .not. isValid ) then
          write(fstderr,'(a)') "thetal mixing failed"
          return
        end if


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
     .                       rhot, rhom, exner, rcm, Ncm,  ! Brian
     .                       hm_rt_tndcy,
     .                       Fcsed, sed_rcm,               ! Brian
     .                       lhs, rhs, C7_Skw_fnc, isValid )
          if ( .not. isValid ) then
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
     .                       rhot, rhom, exner, rcm, Ncm,  ! Brian
     .                       hm_thl_tndcy,
     .                       Fcsed, sed_rcm,               ! Brian
     .                       lhs, rhs, C7_Skw_fnc, isValid )
          if ( .not. isValid ) then
            write(fstderr,'(a)') "Passive scalar 2 mixing failed"
            return
          end if
        end if ! passive scalars
        
!       Error Report
!       Joshua Fasching Feb 2008          
        if ( .not. isValid ) then
                
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
           write(fstderr,*) "rcm = ", rcm
           write(fstderr,*) "Ncm = ", Ncm
           write(fstderr,*) "hm_rt_tndcy", hm_rt_tndcy
           write(fstderr,*) "hm_thl_tndcy", hm_thl_tndcy
           
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
           
           write(fstderr,*) "Intent(out)"
           
           write(fstderr,*) "Fcsed = ", Fcsed
           write(fstderr,*) "sed_rcm", sed_rcm
           
        endif
 
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
        use diffusion, only: diffusion_mlev
#ifdef STATS
        use hoc_stats
#endif /*STATS*/

        implicit none

        ! External
        intrinsic :: min, max

        ! Input variables
        logical, intent(in) :: liter

        real, intent(in) :: dt ! Timestep             [s]

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
        real, intent(out)      :: lhs(2+2+1,2*gr%nnzp) ! LHS

        ! Local Variables
        real, dimension(gr%nnzp) :: a1

        ! Added by Brian for eddy diffusion.
        real, dimension(3,gr%nnzp) :: diff_lhs_imp
        ! Constant parameters
        integer, parameter ::
     .  kp1diag = 1,
     .  kdiag   = 2,
     .  km1diag = 3

        ! Indices
        integer :: k, km1, kp1
        integer :: k_xm, k_wpxp

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

          ! Diffusion
          diff_lhs_imp(kp1diag:km1diag,k)
     .    = diffusion_mlev( Kw6(k), Kw6(kp1), nu6,
     .                      gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )

!         wpxp: mean vertical advection, turbulent transport, time tendency, 
!         damping and diffusion

          ! Momentum subdiagonal [ x wpxp(k-1,<t+1>) ]
          !lhs(k_wpxp,3-2)
          ! Note:  The turbulent advection term, which is
          !        -d[ a1 * (wp3/wp2) * wpxp ] / dz, still
          !        keeps the a1 term inside the derivative,
          !        unlike the wp3 equation and the equations
          !        in diag_var:  rtp2, thlp2, rtpthlp, up2, 
          !        vp2, sclrprtp, sclrpthlp, and sclrp2.  Brian.
          lhs(3+2,k_wpxp)
     .    = - wmm(k) * (1.0/2.0) * gr%dzm(k)
     .      - gr%dzm(k)
     .        * (1.0/2.0) * ( a1(k) + a1(km1) )
     .        * ( wp3(k) / MAX( wp2(km1) + wp2(k), 2.0*eps ) )
     .      + diff_lhs_imp(3,k)

          ! Momentum main diagonal [ x wpxp(k,<t+1>) ]
          lhs(3,k_wpxp)
     .    = + gr%dzm(k)
     .        * (  (1.0/2.0) * ( a1(kp1) + a1(k) )
     .           * ( wp3(kp1) / MAX( wp2(k) + wp2(kp1), 2.0*eps ) )
     .          -  (1.0/2.0) * ( a1(k) + a1(km1) )
     .           * ( wp3(k) / MAX( wp2(km1) + wp2(k), 2.0*eps ) )
     .          )
     .      + (1.0-C7_Skw_fnc(k)) * gr%dzm(k) * ( wmt(kp1) - wmt(k) )
!     .      + C6x / taum(k)
     .      + C6x_Skw_fnc(k) / taum(k)
     .      + diff_lhs_imp(2,k)
          if (liter) lhs(3,k_wpxp) = lhs(3,k_wpxp) + 1.0 / dt

          ! Momentum superdiagonal [ x wpxp(k+1,<t+1>) ]
          !lhs(k_wpxp,3+2)
          lhs(3-2,k_wpxp)
     .    = + wmm(k) * (1.0/2.0) * gr%dzm(k)
     .      + gr%dzm(k)
     .        * (1.0/2.0) * ( a1(kp1) + a1(k) )
     .        * ( wp3(kp1) / MAX( wp2(k) + wp2(kp1), 2.0*eps ) )
     .      + diff_lhs_imp(1,k)

!         wpxp: production

          ! Thermodynamic subdiagonal [ x xm(k,<t+1>) ]
          lhs(3+1,k_wpxp)
     .    = - wp2(k) * gr%dzm(k)

          ! Thermodynamic superdiagonal [ x xm(k+1,<t+1>) ]
          lhs(3-1,k_wpxp)
     .    = + wp2(k) * gr%dzm(k)

#ifdef STATS
          if (lstats_samp) then

            ! Implicit contributions for rtm or wprtp
            if ( irtm_ma > 0 .or. ithlm_ma > 0 ) then
              ztscr01(k) = + wmt(k) * (1.0/2.0) * gr%dzt(k)
              ztscr02(k) = - wmt(k) * (1.0/2.0) * gr%dzt(k)
            end if ! ma

            if ( irtm_ta > 0 .or. ithlm_ta > 0 ) then
              ztscr03(k) = + gr%dzt(k)
              ztscr04(k) = - gr%dzt(k)
            end if ! ta

            ! Implicit contributions for wprtp or wpthlp
            if ( iwprtp_ma > 0 .or. iwpthlp_ma > 0 ) then
              zmscr01(k) = + wmm(k) * (1.0/2.0) * gr%dzm(k)
              zmscr02(k) = - wmm(k) * (1.0/2.0) * gr%dzm(k)
            end if ! ma

            if ( iwprtp_ta > 0 .or. iwpthlp_ta > 0 ) then
            zmscr03(k) =
     .        + gr%dzm(k)
     .          * (1.0/2.0) * ( a1(k) + a1(km1) )
     .          * ( wp3(k) / MAX( wp2(km1) + wp2(k), 2.0*eps ) )
            zmscr04(k) = 
     .        - gr%dzm(k)
     .          * (  (1.0/2.0) * ( a1(kp1) + a1(k) )
     .             * ( wp3(kp1) / MAX( wp2(k) + wp2(kp1), 2.0*eps ) )
     .            -  (1.0/2.0) * ( a1(k) + a1(km1) )
     .             * ( wp3(k) / MAX( wp2(km1) + wp2(k), 2.0*eps ) )
     .            )
            zmscr05(k) = 
     .        - gr%dzm(k)
     .          * (1.0/2.0) * ( a1(kp1) + a1(k) )
     .          * ( wp3(kp1) / MAX( wp2(k) + wp2(kp1), 2.0*eps ) )
            end if ! ta

            if ( iwprtp_tp > 0 .or. iwpthlp_tp > 0 ) then
              zmscr06(k) = + wp2(k) * gr%dzm(k)
              zmscr07(k) = - wp2(k) * gr%dzm(k)
            end if

            if ( iwprtp_ac > 0 .or. iwpthlp_ac > 0 ) then
              zmscr08(k) = - 1.0 * gr%dzm(k) * ( wmt(kp1) - wmt(k) )
            end if

            if ( iwprtp_pr1 > 0 .or. iwpthlp_pr1 > 0 ) then
!               zmscr09(k) = - C6x / taum(k)
                zmscr09(k) = - C6x_Skw_fnc(k) / taum(k)
            end if ! pr1

            if ( iwprtp_pr2 > 0 .or. iwpthlp_pr2 > 0 ) then
              zmscr10(k) = + C7_Skw_fnc(k) 
     .                       * gr%dzm(k) * ( wmt(kp1) - wmt(k) )
            end if

            if ( iwprtp_dp1 > 0 .or. iwpthlp_dp1 > 0 ) then
              zmscr11(k) = - diff_lhs_imp(3,k)
              zmscr12(k) = - diff_lhs_imp(2,k)
              zmscr13(k) = - diff_lhs_imp(1,k)
            end if

          end if
#endif
        end do ! 2..gr%nnzp-1

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
        
        subroutine mixing_solve( solve_type, liter, dt, taum,
     .                           xm, wpxp, xpthvp, xm_forcing,
     .                           wp2, xp2, ! For clipping Vince Larson 20070929
     .                           rhot, rhom, exner, rcm, Ncm,  ! Brian
     .                           hm_tndcy,
     .                           Fcsed, sed_rcm,               ! Brian
     .                           lhs, rhs, C7_Skw_fnc, isValid )

        use parameters
        use constants
        use grid_class
        use lapack_wrap, only: band_solve, band_solvex
        use model_flags, only:
     .    cloud_sed,   ! Cloud sedimentation flag -Brian
     .    kk_rain,     ! Khairoutdinov and Kogan (2000) drizzle scheme -Brian
! Michael Falk, 27 Jul 2007
     .    licedfs,     ! Simplified ice scheme
! eMFc
     .    lcoamps_micro! COAMPS(R) microphysics
#ifdef STATS
        use hoc_stats
#endif /*STATS*/

        implicit none


        ! Input
        character(len=*), intent(in) :: 
     .  solve_type ! Variables being solved for

        logical, intent(in) :: 
     .  liter ! Whether the mean and flux terms are prognostic

        real, intent(in) :: dt      ! Timestep              [s]

        real, intent(in), dimension(gr%nnzp) :: 
     .  taum,        ! Tau on moment. grid                    [s]
     .  rhot,        ! Rho (density) on thermo. grid -Brian   [kg/m^3]
     .  rhom,        ! Rho (density) on moment. grid -Brian   [kg/m^3]
     .  exner,       ! Exner function -Brian                  [-]
     .  rcm,         ! Liquid water mixing ratio -Brian       [kg/kg] 
     .  Ncm,         ! Cloud droplet number conc. -Brian      [num/kg]
     .  hm_tndcy,    ! Change in r_t or theta_l due to micro  [(units)/s]  
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

        real, intent(inout), dimension(nsup+nsub+1,2*gr%nnzp) :: 
     .  lhs  ! Left hand side (LAPACK)
        real, intent(inout), dimension(2*gr%nnzp,nrhs) :: 
     .  rhs  ! Right hand side (LAPACK)

        ! Output
        real, intent(out), dimension(gr%nnzp) :: 
     .  Fcsed,  ! Cloud water sedimentation flux.-Brian       [kg/(m^2 s)]
     .  sed_rcm ! d(rcm)/dt due to cloud sedimentation.-Brian [kg/(m^2 s)]

        logical, intent(out) :: isValid

        ! Local Variables
        real, dimension(2*gr%nnzp,nrhs) :: 
     .  solution ! Solution to band diagonal system

        real :: rcond ! Est. of the reciprocal of the condition #

        ! sigma_g is the Geometric standard deviation of cloud droplets falling
        ! in a Stokes regime.
        real, parameter :: sigma_g = 1.5  ! The DYCOMS2 RF02 intercomparison
                                          ! specified value is 1.5.
!        real, parameter :: sigma_g = 1.2  ! Bjorn Stevens suggests a value of
!                                          ! 1.2 is more appropriate for RF02.

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
     .    iwpxp_cl

        select case ( trim( solve_type ) )
        case ( "rtm" ) ! rtm/wprtp budget terms
          ixm_bt    = irtm_bt
          ixm_ta    = irtm_ta
          ixm_ma    = irtm_ma
          ixm_f     = irtm_forcing
          ixm_mc    = irtm_mc
          ixm_cn    = irtm_cn ! Affects both the mean and flux terms
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
        case ( "thlm" ) ! thlm/wpthlp budget terms
          ixm_bt    = ithlm_bt
          ixm_ta    = ithlm_ta
          ixm_ma    = ithlm_ma
          ixm_f     = ithlm_forcing
          ixm_mc    = ithlm_mc
          ixm_cn    = ithlm_cn ! Affects both the mean and flux terms
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
        case default  ! this includes the sclrm case
          ixm_bt    = 0
          ixm_ta    = 0
          ixm_ma    = 0
          ixm_f     = 0
          ixm_mc    = 0
          ixm_cn    = 0
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
        end select
#endif /*STATS*/

#ifdef STATS
        if ( lstats_samp ) then

          ! xm total time tendency ( 1st calculation)
          if ( ixm_bt > 0 ) then
            zt%x(:,ixm_bt) = zt%x(:,ixm_bt) - xm / dt
          end if

          ! Explicit contributions for xm (not including drizzle)
          if ( ixm_f > 0 ) then
            zt%x(:,ixm_f) = zt%x(:,ixm_f) + xm_forcing
            zt%n(:,ixm_f) = zt%n(:,ixm_f) + 1
          end if
          ! Explicit contributions for xm from microphysics 
          ! (including drizzle )
! Michael Falk, 27 Jul 2007
!          if ( ixm_mc > 0 .and. ( kk_rain .or. lcoamps_micro ) ) then
          if ( ixm_mc > 0 .and. 
     .        ( kk_rain .or. lcoamps_micro .or. licedfs) ) then
! eMFc
            zt%x(:,ixm_mc) = zt%x(:,ixm_mc) + hm_tndcy
            zt%n(:,ixm_mc) = zt%n(:,ixm_mc) + 1
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

        !=====================================================================
        ! NOTE:  ADDITION OF RAIN EFFECTS AND CLOUD WATER SEDIMENTATION ON
        !        RTM AND THLM.
        !=====================================================================
        !
        ! Equations:  rtm = rvm + rcm
        !            thlm = thm - ( Lv / (Cp*exner) ) * rcm
        !
        ! When water condenses, latent heat is given off and theta (thm)
        ! increases by a factor of ( Lv / (Cp*exner) ) * rcm(condensed).
        ! The opposite effect occurs with evaporation.
        !
        !=====================================================================
        !||     Effect     |  rvm  |  rcm  |  rtm  |  rrm  |  thm  |  thlm  ||
        !|===================================================================|
        !|| Sedimentation  |       |       |       |       |       |        ||
        !|| Effects of     | stays | incr. | incr. | stays | stays | decr.  ||
        !|| Cloud Water.   | same  |       |       | same  | same  |        ||
        !|| sed_rcm > 0    |       |       |       |       |       |        ||
        !|===================================================================|
        !|| Evaporation    |       |       |       |       |       |        ||
        !|| of rain to     | incr. | stays | incr. | decr. | decr. | decr.  ||
        !|| water vapor.   |       | same  |       |       |       |        ||
        !|| cond_rrm < 0   |       |       |       |       |       |        ||
        !|-------------------------------------------------------------------|
        !|| Autoconversion |       |       |       |       |       |        ||
        !|| of cloud water | stays | decr. | decr. | incr. | stays | incr.  ||
        !|| to rain water. | same  |       |       |       | same  |        ||
        !|| auto_rrm > 0   |       |       |       |       |       |        ||
        !|-------------------------------------------------------------------|
        !|| Accretion of   |       |       |       |       |       |        ||
        !|| cloud water by | stays | decr. | decr. | incr. | stays | incr.  ||
        !|| rain water.    | same  |       |       |       | same  |        ||
        !|| accr_rrm > 0   |       |       |       |       |       |        ||
        !=====================================================================
        !
        ! Note: In HOC, cond_rrm will always be either negative or zero.
        !
        ! Overall effects of rain and cloud water sedimentation:
        !
        ! (drtm/dt)t  = (drtm/dt)0 
        !                       + sed_rcm - cond_rrm - auto_rrm - accr_rrm
        ! (dthlm/dt)t = (dthlm/dt)0  -  ( Lv / (Cp*exner) ) 
        !                       * ( sed_rcm - cond_rrm - auto_rrm - accr_rrm )
        !
        ! Note by Brian Griffin.
        !=====================================================================

        ! Code addition by Brian for cloud water sedimentation.
        !
        ! Sedimentation flux of cloud droplets should be treated by assuming
        ! a log-normal size distribution of droplets falling in a Stoke's
        ! regime, in which the sedimentation flux (Fcsed) is given by:
        !
        ! Sedimentation Flux = constant
        !                     *[(3/(4*pi*rho_lw*NcV))^(2/3)]
        !                     *[(rho*rc)^(5/3)]
        !                     *EXP[5*((LOG(sigma_g))^2)]
        ! constant = 1.19 x 10^8 (m^-1 s^-1)
        ! sigma_g:  geometric standard deviation
        !
        ! When written for a mass-dependent cloud droplet concentration, Nc:
        !
        ! Sedimentation Flux = constant
        !                     *[(3/(4*pi*rho_lw*Nc*rho))^(2/3)]
        !                     *[(rho*rc)^(5/3)]
        !                     *EXP[5*((LOG(sigma_g))^2)]
        !
        ! According to the above equation, sedimentation flux
        ! is defined positive downwards.  Therefore, 
        !
        ! (drc/dt)Fcsed = (1.0/rho) * d(Fcsed)/dz

        ! Define cloud water sedimentation flux on momentum levels.
        IF ( cloud_sed ) THEN

           DO k = 2, gr%nnzp-1, 1

              IF ( zt2zm(rcm,k) > 0.0 .AND. zt2zm(Ncm,k) > 0.0 ) THEN
                 Fcsed(k) = 1.19E8
     .                  * (   (  3.0 / ( 4.0*pi*rho_lw
     .                                   *zt2zm(Ncm,k)*rhom(k) )  )
     .                     **(2.0/3.0)   )
     .                  * ( ( rhom(k)*zt2zm(rcm,k) )**(5.0/3.0) )
     .                  * EXP( 5.0*( (LOG(sigma_g))**(2) ) )
              ELSE
                 Fcsed(k) = 0.0
              END IF

           END DO ! k=2..gr%nnzp-1

           ! Boundary conditions.
           Fcsed(1)       = 0.0
           Fcsed(gr%nnzp) = 0.0

           ! Find drc/dt due to cloud water sedimentation flux.
           ! This value is defined on thermodynamic levels.
           ! Fcsed units:  kg (liquid) / [ m^2 * s ]
           ! Multiply by Lv for units of W / m^2.
           ! sed_rcm units:  [ kg (liquid) / kg (air) ] / s
           sed_rcm = (1.0/rhot) * ddzm( Fcsed )

        END IF ! cloud_sed


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

          ! + xm_microphysics -- cloud water sedimentation.
          ! Code addition by Brian for cloud water sedimentation.
          IF ( cloud_sed ) THEN

             select case ( trim( solve_type ) )
             case ( "rtm" ) 
                rhs(k_xm,1) = rhs(k_xm,1) + sed_rcm(k)
             case ( "thlm" )
                rhs(k_xm,1) = rhs(k_xm,1) - ( Lv / (Cp*exner(k)) )
     .                               * sed_rcm(k)
             ! The following cases used the formulas for rt and theta in order
             ! to test the scalar code.  However, they need to be commented out
             ! because they do not apply to such things as CO2.  Cloud water
             ! sedimentation would have no effect on such a parameter.
!             case ( "sclrm_rtm" ) 
!                rhs(k_xm) = rhs(k_xm) + sed_rcm(k)
!             case ( "sclrm_thlm" )
!                rhs(k_xm) = rhs(k_xm) - ( Lv / (Cp*exner(k)) )
!     .                               * sed_rcm(k)
             end select

          END IF ! cloud_sed

          ! + xm_microphysics -- rain effects.
          ! Code addition by Brian for rain.
! Michael Falk, 27 Jul 2007
!          if ( kk_rain .or. lcoamps_micro ) then
          if ( kk_rain .or. lcoamps_micro .or. licedfs ) then
! eMFc
            select case(trim( solve_type ))
            case ( "rtm", "thlm" ) 
            rhs(k_xm,1) = rhs(k_xm,1) + hm_tndcy(k)  ! thlm/rtm only; Brian

            ! The following cases used the formulas for rt and theta in order
            ! to test the scalar code.  However, they need to be commented out
            ! because they do not apply to such things as CO2.  Rain water
            ! effects would have no effect on such a parameter.
!            case ( "sclrm_rtm" ) 
!              rhs(k_xm) = rhs(k_xm) 
!     .                 - cond_rrm(k) - auto_rrm(k) - accr_rrm(k)    ! sclrm = rtm
!            case ( "sclrm_thlm" )
!              rhs(k_xm) = rhs(k_xm) - ( Lv / (Cp*exner(k)) ) *
!     .                    ( -cond_rrm(k) - auto_rrm(k) - accr_rrm(k) )  ! sclrm = thlm
             case default
               ! Add nothing
             end select
           end if ! kk_rain .or lcoamps_micro

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
     .                      lhs, rhs, solution, rcond, isValid )

          ! Est. of the condition number of the mean/flux LHS matrix 
          sfc%x(1,ixm_cn) = sfc%x(1,ixm_cn) + 1.0 / rcond
          sfc%n(1,ixm_cn) = sfc%n(1,ixm_cn) + 1

        else
          ! Perform LU decomp and solve system (LAPACK)
          call band_solve( solve_type, nsup, nsub, 2*gr%nnzp, nrhs,
     .                     lhs, rhs, solution, isValid )
        end if
#else
        ! Perform LU decomp and solve system (LAPACK)
        call band_solve( solve_type, nsup, nsub, 2*gr%nnzp, nrhs,
     .                   lhs, rhs, solution, isValid )

#endif /*STATS*/

        if ( .not. isValid ) return

!       Copy result into output arrays

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

           do k = 1, gr%nnzp
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

        end module mixing
