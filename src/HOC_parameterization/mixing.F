!-----------------------------------------------------------------------
! $Id: mixing.F,v 1.10 2006-12-12 20:21:17 dschanen Exp $

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

        public  :: timestep_mixing

        private :: mixing_solve, mixing_lhs

        contains
!-----------------------------------------------------------------------
        subroutine timestep_mixing
     .             ( dt, taum, wmm, wmt, wp2, wp3, Scm, Skwm,
     .               rtm, wprtp, rtpthvp, rtm_forcing,
     .               thlm, wpthlp, thlpthvp, thlm_forcing,
     .               rhot, rhom, exner, rcm, Ncm, 
     .               hm_thl_tndcy, hm_rt_tndcy,
     .               Fcsed, sed_rcm,
     .               isValid, implemented,
     .               sclrm, wpsclrp, sclrpthvp, sclrm_forcing )
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

        use constants
        use grid_class

        implicit none

        ! External
        intrinsic :: exp, present

        ! Input Variables
        real, intent(in) :: 
     .  dt           ! Timestep                       [s]

        real, intent(in), dimension(gr%nnzp) ::
     .  taum,        ! Tau on moment. grid            [s]
     .  wmm,         ! w wind on m                    [m/s]
     .  wmt,         ! w wind on t                    [m/s]
     .  wp2,         ! w'^2                           [m^2/s^2]
     .  wp3,         ! w'^3                           [m^3/s^3]
     .  Scm,         ! Sc on moment. grid             [-]
     .  Skwm,        ! Skw on moment. grid            [-]
     .  rtpthvp,     ! r_t'th_v'                      [(kg K)/kg]
     .  rtm_forcing, ! rt forcing                     [(kg/kg)/s]
     .  thlpthvp,    ! th_l' th_v'                    [K^2]
     .  thlm_forcing ! Large scale tendency for thlm  [K/s]

        ! Additions for DYCOMS_2
        real, intent(in), dimension(gr%nnzp) ::
     .  rhot,        ! rho (density) on thermo. grid  [kg/m^3] 
     .  rhom,        ! rho (density) on moment. grid  [kg/m^3] 
     .  exner,       ! Exner function                 [-]
     .  rcm,         ! Liquid water mixing ratio      [kg/kg]
     .  Ncm,         ! Cloud droplet number conc.     [num/m^3]
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
        ! Input/Output Variables
        real, intent(inout), optional, dimension(gr%nnzp,sclr_dim) :: 
     .  sclrm, wpsclrp !                                     [Units vary]

        ! Local variables
        real, allocatable, dimension(:,:) :: 
     .  lhs  ! Banded matrix in compact form 
!    .  al  ! LU decomposition of banded matrix 

        real, allocatable, dimension(:) :: 
     .  dd ! RHS for banded matrix

!       integer, allocatable, dimension(:) :: 
!    .  indx ! Index of pivots for LU decomposition

        ! Constant parameters with Skw applied
        real, allocatable, dimension(:) :: 
     .  C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc

        ! Misc. storage
        real, dimension(gr%nnzp) :: tmp1, tmp2
        
        ! Indices
        integer :: i
        integer :: k, kp1

!-----------------------------------------------------------------------

        ! Allocate memory
        allocate( C6rt_Skw_fnc(1:gr%nnzp) )
        allocate( C6thl_Skw_fnc(1:gr%nnzp) )
        allocate( C7_Skw_fnc(1:gr%nnzp) )

        !allocate( lhs(1:2*gr%nnzp,1:5) ) ! (NR)
        allocate( lhs(1:5,1:2*gr%nnzp) )  ! (LAPACK)
        allocate( dd(1:2*gr%nnzp) )
        !allocate( al(1:2*gr%nnzp,1:2) )
        !allocate( indx(1:2*gr%nnzp) )

        C6rt_Skw_fnc = C6rtb + (C6rt-C6rtb)
     .                        *EXP( -0.5 * (Skwm/C6rtc)**2 )
        C6thl_Skw_fnc = C6thlb + (C6thl-C6thlb)
     .                          *EXP( -0.5 * (Skwm/C6thlc)**2 )
        C7_Skw_fnc = C7b + (C7-C7b)
     .                    *EXP( -0.5 * (Skwm/C7c)**2 )

!        C6rt_Skw_fnc = C6rt
!        C6thl_Skw_fnc = C6thl
!        C7_Skw_fnc = C7

!       Setup and decompose matrix

!        call mixing_lhs( .true., dt, taum, wmm, wmt, wp2, wp3, Scm,
!       Prof. Larson's modified mixing_lhs 
        call mixing_lhs( "rtm", .true., dt, taum, wmm, wmt, 
     .                   wp2, wp3, Scm,
     .                   C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc,
     .                   lhs, implemented )

        ! Solve rtm / wprtp
        call mixing_solve( "rtm", .true., dt, taum,
     .                     rtm, wprtp, rtpthvp, rtm_forcing,
     .                     rhot, rhom, exner, rcm, Ncm,  ! Brian
     .                     hm_rt_tndcy,
     .                     Fcsed, sed_rcm,               ! Brian
     .                     lhs, dd, C7_Skw_fnc, isValid )  
        if ( .not. isValid ) then
          print *, "rt mixing failed"
          deallocate( lhs, dd, C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc )
          return
        end if

        ! Prof. Larson's modified mixing_lhs 
        call mixing_lhs( "thlm",.true., dt, taum, wmm, wmt, 
     .                   wp2, wp3, Scm,
     .                   C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc,
     .                   lhs, implemented )

        ! Solve for thlm / wpthlp
        call mixing_solve( "thlm", .true., dt, taum,
     .                     thlm, wpthlp, thlpthvp, thlm_forcing,
     .                     rhot, rhom, exner, rcm, Ncm,  ! Brian
     .                     hm_thl_tndcy,
     .                     Fcsed, sed_rcm,               ! Brian
     .                     lhs, dd, C7_Skw_fnc, isValid )
        if ( .not. isValid ) then
          print *, "theta mixing failed"
          deallocate( lhs, dd, C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc )
          return
        end if

        ! Solve sclrm / wpsclrp
        if ( present( sclrm )     .and. present( wpsclrp ) .and.
     .       present( sclrpthvp ) .and. present( sclrm_forcing ) ) then

          ! Prof. Larson's modified mixing_lhs 
          call mixing_lhs( "rtm",.true., dt, taum, wmm, wmt, 
     .                     wp2, wp3, Scm,
     .                     C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc,
     .                     lhs, implemented )

          call mixing_solve( "sclrm_rtm", .true., dt, taum,
     .                       sclrm(:,SCLR_RT), wpsclrp(:,SCLR_RT),
     .                       sclrpthvp(:,SCLR_RT), 
     .                       sclrm_forcing(:,SCLR_RT),
     .                       rhot, rhom, exner, rcm, Ncm,  ! Brian
     .                       hm_rt_tndcy,
     .                       Fcsed, sed_rcm,               ! Brian
     .                       lhs, dd, C7_Skw_fnc, isValid )
          if ( .not. isValid ) then
            print *, "sclr rt mixing failed"
            deallocate( lhs, dd, C6rt_Skw_fnc, 
     .                  C6thl_Skw_fnc, C7_Skw_fnc )
            return
          end if

          ! Prof. Larson's modified mixing_lhs 
          call mixing_lhs( "thlm",.true., dt, taum, wmm, wmt, 
     .                     wp2, wp3, Scm, 
     .                     C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc,
     .                     lhs, implemented )

          call mixing_solve( "sclrm_theta", .true., dt, taum,
     .                       sclrm(:,SCLR_THETA), wpsclrp(:,SCLR_THETA),
     .                       sclrpthvp(:,SCLR_THETA), 
     .                       sclrm_forcing(:,SCLR_THETA),
     .                       rhot, rhom, exner, rcm, Ncm,  ! Brian
     .                       hm_thl_tndcy,
     .                       Fcsed, sed_rcm,               ! Brian
     .                       lhs, dd, C7_Skw_fnc, isValid )
          if ( .not. isValid ) then
            print *, "sclr theta mixing failed"
            deallocate( lhs, dd, C6rt_Skw_fnc, 
     .                  C6thl_Skw_fnc, C7_Skw_fnc )
            return
          end if
        end if ! passive scalars
 
        ! Free memory
        deallocate( lhs, dd, C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc )

        return
        end subroutine timestep_mixing

!------------------------------------------------------------------------
        subroutine mixing_lhs( solve_type, liter, dt, taum, wmm, wmt,
     .                         wp2, wp3, Scm, 
     .                         C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc,
     .                         lhs, implemented )
!       Description:
!       References:
!------------------------------------------------------------------------

        use constants
        use grid_class
#ifdef STATS
        use statistics
#endif /*STATS*/

        implicit none

        ! Input variables

        character(len=*), intent(in) :: solve_type  ! For new mixing scheme

        logical, intent(in) :: liter

        real, intent(in) :: dt ! Timestep             [s]
        real, intent(in), dimension(gr%nnzp) :: 
     .  taum,        ! Tau on moment. grid            [s]
     .  wmm,         ! w wind on m                    [m/s]
     .  wmt,         ! w wind on t                    [m/s]
     .  wp2,         ! w'^2                           [m^2/s^2]
     .  wp3,         ! w'^3                           [m^3/s^3]
     .  Scm          ! Sc on moment. grid             [-]

        ! Output variables

!       real, intent(out)      :: lhs(1:2*gr%nnzp,1:5) !(NR)
!       real, intent(out)      :: dd(1:2*gr%nnzp)
!       real, intent(out)      :: al(1:2*gr%nnzp,1:2)
!       integer, intent(out)   :: indx(1:2*gr%nnzp)
        real, intent(out)      :: lhs(1:5,1:2*gr%nnzp) !(LAPACK)
        logical, intent(in)    :: implemented

        real, intent(in), dimension(gr%nnzp) :: 
     .  C6rt_Skw_fnc,
     .  C6thl_Skw_fnc,
     .  C7_Skw_fnc


        ! Local Variables
        ! real :: d

        real, dimension(gr%nnzp) :: a1

        ! Indices
        integer :: k, km1, kp1
        integer :: k_xm, k_wpxp

!       Define a1

        a1(1:gr%nnzp) = 1.0 / ( 1.0 - Scm(1:gr%nnzp) )

!       Set up left hand side matrix

        lhs = 0.0

!        do k=2,gr%nnzp
        do k=2,gr%nnzp-1

!         Define indices

          km1 = max( k-1, 1 )
          kp1 = min( k+1, gr%nnzp )

          k_xm   = 2*k - 1
          k_wpxp = 2*k

!         xm: mean vertical advection

          if ( .not. implemented ) then

            ! Thermodynamic subdiagonal [ x xm(k-1,<t+1>) ]
            !lhs(k_xm,3-2) = - wmt(k) * (1.0/2.0) * gr%dzt(k)
            lhs(3+2,k_xm) = - wmt(k) * (1.0/2.0) * gr%dzt(k)
            ! Thermodynamic superdiagonal [ x xm(k+1,<t+1>) ]
            !lhs(k_xm,3+2) = + wmt(k) * (1.0/2.0) * gr%dzt(k)
            lhs(3-2,k_xm) = + wmt(k) * (1.0/2.0) * gr%dzt(k)

          else

            ! Thermodynamic subdiagonal [ x xm(k-1,<t+1>) ]
            !lhs(k_xm,3-2) = 0.0
            lhs(3+2,k_xm) = 0.0
            ! Thermodynamic superdiagonal [ x xm(k+1,<t+1>) ]
            !lhs(k_xm,3+2) = 0.0
            lhs(3-2,k_xm) = 0.0

          end if

!         xm: turbulent transport

          ! Momentum subdiagonal [ x wpxp(k-1,<t+1>) ]
          !lhs(k_xm,3-1) = - gr%dzt(k)
          lhs(3+1,k_xm) = - gr%dzt(k)
          ! Momentum superdiagonal [ x wpxp(k,<t+1>) ]
          !lhs(k_xm,3+1) = + gr%dzt(k)
          lhs(3-1,k_xm) = + gr%dzt(k)

!         xm: time tendency

          ! Thermodynamic main diagonal [ x xm(k,<t+1>) ]
          !lhs(k_xm,3) = 1.0 / dt
          lhs(3,k_xm) = 1.0 / dt

!         wpxp: mean vertical advection, turbulent transport, time tendency, 
!         damping and diffusion

          ! Momentum subdiagonal [ x wpxp(k-1,<t+1>) ]
          !lhs(k_wpxp,3-2)
          lhs(3+2,k_wpxp)
     .    = - wmm(k) * (1.0/2.0) * gr%dzm(k)
     .      - gr%dzm(k)
     .        * (1.0/2.0) * ( a1(k) + a1(km1) )
     .        * ( wp3(k) / MAX( wp2(km1) + wp2(k), 2.0*eps ) )
     .      - nu6 * gr%dzm(k) * gr%dzt(k)      

          ! Momentum main diagonal [ x wpxp(k,<t+1>) ]
          select case( trim( solve_type ) )
          case ("rtm")
            !lhs(k_wpxp,3)
            lhs(3,k_wpxp)
     .      = + gr%dzm(k)
     .          * (  (1.0/2.0) * ( a1(kp1) + a1(k) )
     .             * ( wp3(kp1) / MAX( wp2(k) + wp2(kp1), 2.0*eps ) )
     .            -  (1.0/2.0) * ( a1(k) + a1(km1) )
     .             * ( wp3(k) / MAX( wp2(km1) + wp2(k), 2.0*eps ) )
     .            )
     .        + (1.0-C7_Skw_fnc(k)) * gr%dzm(k) * ( wmt(kp1) - wmt(k) )
!     .        + C6rt / taum(k)
     .        + C6rt_Skw_fnc(k) / taum(k)
     .        + nu6 * gr%dzm(k) * ( gr%dzt(kp1) + gr%dzt(k) )
            !if (liter) lhs(k_wpxp,3) = lhs(k_wpxp,3) + 1.0 / dt
            if (liter) lhs(3,k_wpxp) = lhs(3,k_wpxp) + 1.0 / dt
          case ("thlm")
            !lhs(k_wpxp,3)
            lhs(3,k_wpxp)
     .      = + gr%dzm(k)
     .          * (  (1.0/2.0) * ( a1(kp1) + a1(k) )
     .             * ( wp3(kp1) / MAX( wp2(k) + wp2(kp1), 2.0*eps ) )
     .            -  (1.0/2.0) * ( a1(k) + a1(km1) )
     .             * ( wp3(k) / MAX( wp2(km1) + wp2(k), 2.0*eps ) )
     .            )
     .        + (1.0-C7_Skw_fnc(k)) * gr%dzm(k) * ( wmt(kp1) - wmt(k) )
!     .        + C6thl / taum(k)
     .        + C6thl_Skw_fnc(k) / taum(k)
     .        + nu6 * gr%dzm(k) * ( gr%dzt(kp1) + gr%dzt(k) )
            !if (liter) lhs(k_wpxp,3) = lhs(k_wpxp,3) + 1.0 / dt
            if (liter) lhs(3,k_wpxp) = lhs(3,k_wpxp) + 1.0 / dt
          end select

          ! Momentum superdiagonal [ x wpxp(k+1,<t+1>) ]
          !lhs(k_wpxp,3+2)
          lhs(3-2,k_wpxp)
     .    = + wmm(k) * (1.0/2.0) * gr%dzm(k)
     .      + gr%dzm(k)
     .        * (1.0/2.0) * ( a1(kp1) + a1(k) )
     .        * ( wp3(kp1) / MAX( wp2(k) + wp2(kp1), 2.0*eps ) )
     .      - nu6 * gr%dzm(k) * gr%dzt(kp1)

!         wpxp: production

          ! Thermodynamic subdiagonal [ x xm(k,<t+1>) ]
          !lhs(k_wpxp,3-1)
          lhs(3+1,k_wpxp)
     .    = - wp2(k) * gr%dzm(k)

          ! Thermodynamic superdiagonal [ x xm(k+1,<t+1>) ]
          !lhs(k_wpxp,3+1)
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
              select case (trim( solve_type ) )
              case("rtm")
!               zmscr09(k) = - C6rt / taum(k)
                zmscr09(k) = - C6rt_Skw_fnc(k) / taum(k)
              case("thlm")
!                zmscr09(k) = - C6thl / taum(k)
                zmscr09(k) = - C6thl_Skw_fnc(k) / taum(k)
              case default
                print *, "Unknown solve_type:", solve_type
                stop "invalid mixing_lhs call"
              end select
            end if ! pr1

            if ( iwprtp_pr2 > 0 .or. iwpthlp_pr2 > 0 ) then
              zmscr10(k) = + C7_Skw_fnc(k) 
     .                       * gr%dzm(k) * ( wmt(kp1) - wmt(k) )
            end if

            if ( iwprtp_dp1 > 0 .or. iwpthlp_dp1 > 0 ) then
              zmscr11(k) = + nu6 * gr%dzm(k) * gr%dzt(k)
              zmscr12(k) = - nu6 * gr%dzm(k) 
     .                             * ( gr%dzt(kp1) + gr%dzt(k) )
              zmscr13(k) = + nu6 * gr%dzm(k) * gr%dzt(kp1)
            end if

          end if
#endif
        end do ! 2..gr%nnzp-1

!       Boundary condition for the flux

        k      = 1
        k_xm   = 2*k - 1
        k_wpxp = 2*k
        ! (NR)
        !lhs(k_xm,3-1)   = 0.0
        !lhs(k_xm,3)     = 1.0
        !lhs(k_xm,3+1)   = 0.0
        !lhs(k_wpxp,3-2) = 0.0
        !lhs(k_wpxp,3)   = 1.0
        !lhs(k_wpxp,3+2) = 0.0
        lhs(3+1,k_xm)   = 0.0
        lhs(3,k_xm)     = 1.0
        lhs(3-1,k_xm)   = 0.0
        lhs(3+2,k_wpxp) = 0.0
        lhs(3,k_wpxp)   = 1.0
        lhs(3-2,k_wpxp) = 0.0

        k      = gr%nnzp
        k_xm   = 2*k - 1
        k_wpxp = 2*k
        ! (NR)
        !lhs(k_xm,3-1)   = 0.0
        !lhs(k_xm,3)     = 1.0
        !lhs(k_xm,3+1)   = 0.0
        !lhs(k_wpxp,3-2) = 0.0
        !lhs(k_wpxp,3)   = 1.0
        !lhs(k_wpxp,3+2) = 0.0
        lhs(3+1,k_xm)   = 0.0
        lhs(3,k_xm)     = 1.0
        lhs(3-1,k_xm)   = 0.0
        lhs(3+2,k_wpxp) = 0.0
        lhs(3,k_wpxp)   = 1.0
        lhs(3-2,k_wpxp) = 0.0

!       Perform LU decomposition (NR)

!       call bandec( lhs, 2*gr%nnzp, 2, 2, 2*gr%nnzp, 5, al, 2, indx, d, 
!    .               isValid )

        return
        end subroutine mixing_lhs

!------------------------------------------------------------------------
        
        subroutine mixing_solve( solve_type, liter, dt, taum,
     .                           xm, wpxp, xpthvp, xm_forcing,
     .                           rhot, rhom, exner, rcm, Ncm,  ! Brian
     .                           hm_tndcy,
     .                           Fcsed, sed_rcm,               ! Brian
     .                           lhs, dd, C7_Skw_fnc, isValid )

        use constants
        use grid_class
        use lapack_wrap, only: band_solve
        use model_flags, only:
     .    cloud_sed,   ! Cloud sedimentation flag -Brian
     .    kk_rain,     ! Khairoutdinov and Kogan (2000) drizzle scheme -Brian
     .    lcoamps_micro! COAMPS(R) microphysics
#ifdef STATS
        use statistics
#endif /*STATS*/

        implicit none


        ! Input
        character(len=*), intent(in) :: 
     .  solve_type ! Variables being solved for

        logical, intent(in) :: 
     .  liter ! Whether the mean and flux terms are prognostic

        real, intent(in) :: dt      ! Timestep              [s]

        real, intent(in), dimension(gr%nnzp) :: 
     .  taum,        ! Tau on moment. grid                  [s]
     .  rhot,        ! Rho (density) on thermo. grid -Brian [kg/m^3]
     .  rhom,        ! Rho (density) on moment. grid -Brian [kg/m^3]
     .  exner,       ! Exner function -Brian                [-]
     .  rcm,         ! Liquid water mixing ratio -Brian     [kg/kg] 
     .  Ncm,         ! Cloud droplet number conc. -Brian    [num/m^3]
     .  hm_tndcy,    ! Change in r_t or theta_l due to micro[(units)/s]  
     .  xpthvp,      ! x' theta_v'                          [(units)(K)]
     .  xm_forcing,  ! xm forcing -Brian                    [(units)/s] 
     .  C7_Skw_fnc   ! C7 parameter (Skw dependent)         [-] 

        ! Input/Output
        real, intent(inout), dimension(gr%nnzp) :: 
     .  xm,  ! Mean term                                  [units vary] 
     .  wpxp ! Flux term                                  [units vary] 

        real, intent(inout), dimension(5,2*gr%nnzp) :: 
     .  lhs  ! Left hand side (LAPACK)
        real, intent(inout), dimension(2*gr%nnzp) :: 
     .  dd   ! Right hand side (LAPACK)

        ! These were the LU decomp. and pivot index, respectively
        !real, intent(in)    :: al(1:2*gr%nnzp,1:2) !(NR)
        !integer, intent(in) :: indx(1:2*gr%nnzp)   !(NR)

        ! Output
        real, intent(out), dimension(gr%nnzp) :: 
     .  Fcsed,  ! Cloud water sedimentation flux.-Brian       [kg/(m^2 s)]
     .  sed_rcm ! d(rcm)/dt due to cloud sedimentation.-Brian [kg/(m^2 s)]

        logical, intent(out) :: isValid

        ! Local Variables
        real, dimension(2*gr%nnzp) :: 
     .  solution ! Solution to band diagonal system

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
     .    iwpxp_bt,
     .    iwpxp_ma,
     .    iwpxp_ta,
     .    iwpxp_tp,
     .    iwpxp_ac,
     .    iwpxp_bp,
     .    iwpxp_pr1,
     .    iwpxp_pr2,
     .    iwpxp_pr3,
     .    iwpxp_dp1

        select case( trim( solve_type ) )
        case( "rtm" ) ! rtm/wprtp budget terms
          ixm_bt    = irtm_bt
          ixm_ta    = irtm_ta
          ixm_ma    = irtm_ma
          ixm_f     = irtm_forcing
          ixm_mc    = irtm_mc
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
        case( "thlm" ) ! thlm/wpthlp budget terms
          ixm_bt    = ithlm_bt
          ixm_ta    = ithlm_ta
          ixm_ma    = ithlm_ma
          ixm_f     = ithlm_forcing
          ixm_mc    = ithlm_mc
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
        case default  ! this includes the sclrm case
          ixm_bt    = 0
          ixm_ta    = 0
          ixm_ma    = 0
          ixm_f     = 0
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
          if ( ixm_mc > 0 .and. ( kk_rain .or. lcoamps_micro ) ) then
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
        !                     *[(3/(4*pi*rho_lw*Nc))^(2/3)]
        !                     *[(rho*rc)^(5/3)]
        !                     *EXP[5*((LOG(sigma_g))^2)]
        ! constant = 1.19 x 10^8 (m^-1 s^-1)
        ! sigma_g:  geometric standard deviation = 1.5
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
     .                  * (  ( 3.0 / (4.0*pi*rho_lw*zt2zm(Ncm,k)) )
     .                     **(2.0/3.0)  )
     .                  * ( ( rhom(k)*zt2zm(rcm,k) )**(5.0/3.0) )
     .                  * EXP( 5.0*( (LOG(1.5))**(2) ) )
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


!       Set up RHS vector dd

        dd = 0.0


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
          dd(k_xm) 
     .    = xm(k) / dt
     .    + xm_forcing(k)

          ! + xm_microphysics -- cloud water sedimentation.
          ! Code addition by Brian for cloud water sedimentation.
          IF ( cloud_sed ) THEN

             select case(trim( solve_type ))
             case ( "rtm" ) 
                dd(k_xm) = dd(k_xm) + sed_rcm(k)
             case ( "thlm" )
                dd(k_xm) = dd(k_xm) - ( Lv / (Cp*exner(k)) )
     .                               * sed_rcm(k)
             ! The following cases used the formulas for rt and theta in order
             ! to test the scalar code.  However, they need to be commented out
             ! because they do not apply to such things as CO2.  Cloud water
             ! sedimentation would have no effect on such a parameter.
!             case ( "sclrm_rtm" ) 
!                dd(k_xm) = dd(k_xm) + sed_rcm(k)
!             case ( "sclrm_thlm" )
!                dd(k_xm) = dd(k_xm) - ( Lv / (Cp*exner(k)) )
!     .                               * sed_rcm(k)
             end select

          END IF ! cloud_sed

          ! + xm_microphysics -- rain effects.
          ! Code addition by Brian for rain.
          if ( kk_rain .or. lcoamps_micro ) then
            select case(trim( solve_type ))
            case ( "rtm", "thlm" ) 
            dd(k_xm) = dd(k_xm) + hm_tndcy(k)  ! thlm/rtm only; Brian

            ! The following cases used the formulas for rt and theta in order
            ! to test the scalar code.  However, they need to be commented out
            ! because they do not apply to such things as CO2.  Rain water
            ! effects would have no effect on such a parameter.
!            case ( "sclrm_rtm" ) 
!              dd(k_xm) = dd(k_xm) 
!     .                 - cond_rrm(k) - auto_rrm(k) - accr_rrm(k)    ! sclrm = rtm
!            case ( "sclrm_thlm" )
!              dd(k_xm) = dd(k_xm) - ( Lv / (Cp*exner(k)) ) *
!     .                    ( -cond_rrm(k) - auto_rrm(k) - accr_rrm(k) )  ! sclrm = thlm
             case default
               ! Add nothing
             end select
           end if ! kk_rain .or lcoamps_micro

!         wpxp: time tendency and buoyancy

          ! Right-hand side.
          dd(k_wpxp)
     .    = grav/T0 * (1.0-C7_Skw_fnc(k)) * xpthvp(k)

          if ( liter ) dd(k_wpxp) = dd(k_wpxp) + wpxp(k) / dt

        end do ! k=2..gr%nnzp-1

        ! Boundary condition for the flux and mean

        k      = 1
        k_xm   = 2*k - 1
        k_wpxp = 2*k
        dd(k_xm)   = xm(k)
        dd(k_wpxp) = wpxp(k)

        k      = gr%nnzp
        k_xm   = 2*k - 1
        k_wpxp = 2*k
        dd(k_xm)   = xm(k)
        dd(k_wpxp) = 0.0

        ! Solve system (NR)

!       call banbks( lhs, 2 * gr%nnzp, 2, 2, 
!    .               2 * gr%nnzp, 5, al, 2, indx, dd )

        ! Perform LU decomposition and solve system
        ! Note: band_solvex has the same interface and can compute
        ! additional diagnostics at a higher computational cost.
        call band_solve( solve_type, 2, 2, 2*gr%nnzp, 1,
     .                   lhs, dd, solution, isValid )

        if ( .not. isValid ) return

!       Copy result into output arrays

        do k=1,gr%nnzp

          k_xm   = 2 * k - 1
          k_wpxp = 2 * k

          !xm(k)   = dd(k_xm)
          !wpxp(k) = dd(k_wpxp)
          xm(k)   = solution(k_xm)
          wpxp(k) = solution(k_wpxp)

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
