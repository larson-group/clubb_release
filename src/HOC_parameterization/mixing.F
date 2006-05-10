!-----------------------------------------------------------------------
! $Id: mixing.F,v 1.2 2006-05-10 18:46:16 griffinb Exp $
!
! Module containing the hoc mixing scheme.  Allows for optional variables
! and assumed shape arrays.
!-----------------------------------------------------------------------
#define SCLR_THETA 1
#define SCLR_RT 2

      module mixing
      implicit none

      public timestep_mixing
      private mixing_solve
      private mixing_lhs

      contains
!-----------------------------------------------------------------------
! Subroutine timestep_mixing( )
!-----------------------------------------------------------------------
        subroutine timestep_mixing
     .             ( dt, taum, wmm, wmt, wp2, wp3, Sc,
     .               rtm, wprtp, rtpthvp, rtm_forcing,
     .               thlm, wpthlp, thlpthvp, thlm_forcing,
     .               rhot, rhom, exner, rcm, Ncm, 
     .               cond_rrm, auto_rrm, accr_rrm,
     .               Fcsed, sed_rcm,
     .               cloud_sed, kk_rain, isValid, implemented,
     .               sclrm, wpsclrp, sclrpthvp, sclrm_forcing )

        use constants
        use grid_class

        implicit none

!       Input/Output

        real, intent(in)    :: dt
        real, intent(in)    :: taum(1:gr%nnzp)
        real, intent(in)    :: wmm(1:gr%nnzp)
        real, intent(in)    :: wmt(1:gr%nnzp)
        real, intent(in)    :: wp2(1:gr%nnzp)
        real, intent(in)    :: wp3(1:gr%nnzp)
        real, intent(in)    :: Sc(1:gr%nnzp)
        real, intent(inout) :: rtm(1:gr%nnzp)
        real, intent(inout) :: wprtp(1:gr%nnzp)
        real, intent(in)    :: rtpthvp(1:gr%nnzp)
        real, intent(in)    :: rtm_forcing(1:gr%nnzp)
        real, intent(inout) :: thlm(1:gr%nnzp)
        real, intent(inout) :: wpthlp(1:gr%nnzp)
        real, intent(in)    :: thlpthvp(1:gr%nnzp)
        real, intent(in)    :: thlm_forcing(1:gr%nnzp)
!       Additions for DYCOMS_2
        real, intent(in)    :: rhot(1:gr%nnzp)    
        real, intent(in)    :: rhom(1:gr%nnzp)    
        real, intent(in)    :: exner(1:gr%nnzp)    
        real, intent(in)    :: rcm(1:gr%nnzp)    
        real, intent(in)    :: Ncm(1:gr%nnzp)    
        real, intent(in)    :: cond_rrm(1:gr%nnzp)
        real, intent(in)    :: auto_rrm(1:gr%nnzp)
        real, intent(in)    :: accr_rrm(1:gr%nnzp)
        real, intent(out)   :: Fcsed(1:gr%nnzp)
        real, intent(out)   :: sed_rcm(1:gr%nnzp)
        logical, intent(in) :: cloud_sed
        logical, intent(in) :: kk_rain
!       End Additions for DYCOMS_2
!       Flag for if HOC is implemented in a larger model or not.
        logical, intent(inout) :: implemented
!       Model status 
        logical, intent(inout) :: isValid
!       Additions for new mixing scheme
        real, intent(inout), optional, 
     .    dimension(:,:) :: sclrm, wpsclrp 
        real, intent(in), optional, 
     .    dimension(:,:) :: sclrpthvp, sclrm_forcing
!       End Additions for new mixing scheme

!       Internal
        real, allocatable    :: aa(:,:), dd(:), al(:,:)
        integer, allocatable :: indx(:)
        integer i

        real, dimension(gr%nnzp) :: tmp1, tmp2

        integer k, kp1
        real, allocatable :: Skw_arr(:)
        real, allocatable :: C6rt_Skw_fnc(:), C6thl_Skw_fnc(:)
        real, allocatable :: C7_Skw_fnc(:)

        allocate( Skw_arr(1:gr%nnzp) )
        allocate( C6rt_Skw_fnc(1:gr%nnzp) )
        allocate( C6thl_Skw_fnc(1:gr%nnzp) )
        allocate( C7_Skw_fnc(1:gr%nnzp) )

!       Allocate memory
        allocate( aa(1:2*gr%nnzp,1:5) )
        allocate( dd(1:2*gr%nnzp) )
        allocate( al(1:2*gr%nnzp,1:2) )
        allocate( indx(1:2*gr%nnzp) )

        do k=1,gr%nnzp

          kp1 = min(k+1,gr%nnzp)
          Skw_arr(k) = (0.5*(wp3(k)+wp3(kp1))) /
     .                    ( max(wtol**2,wp2(k)) )**1.5
        enddo

        C6rt_Skw_fnc = C6rtb + (C6rt-C6rtb)
     .                        *EXP( -0.5 * (Skw_arr/C6rtc)**2.0 )
        C6thl_Skw_fnc = C6thlb + (C6thl-C6thlb)
     .                          *EXP( -0.5 * (Skw_arr/C6thlc)**2.0 )
        C7_Skw_fnc = C7b + (C7-C7b)
     .                    *EXP( -0.5 * (Skw_arr/C7c)**2.0 )

!        C6rt_Skw_fnc = C6rt
!        C6thl_Skw_fnc = C6thl
!        C7_Skw_fnc = C7

!       Setup and decompose matrix

!        call mixing_lhs( .true., dt, taum, wmm, wmt, wp2, wp3, Sc,
!   Prof. Larson's modified mixing_lhs 
        call mixing_lhs("rtm", .true., dt, taum, wmm, wmt, wp2, wp3, Sc,
     .                   C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc,
     .                   aa, dd, al, indx, implemented, isValid )
        if ( .not. isValid ) then
           print *, "mixing_lhs failed"
           deallocate( aa, dd, al, indx )
           return
        endif
!        print *, "before"
!        print *, "rt=",  sum( sclrm )   - sum ( rtm ),
!     .           "wp=",  sum( wprtp )   - sum ( wpsclrp ),
!     .           "thvp=",sum( rtpthvp ) - sum ( sclrpthvp )

!       Solve rtm / wprtp
        call mixing_solve( "rtm", .true., dt, taum,
     .                     rtm, wprtp, rtpthvp, rtm_forcing,
     .                     rhot, rhom, exner, rcm, Ncm,  ! Brian
     .                     cond_rrm, auto_rrm, accr_rrm, ! Brian
     .                     Fcsed, sed_rcm,               ! Brian
     .                     cloud_sed, kk_rain,           ! Brian
     .                     aa, dd, al, indx, C7_Skw_fnc )  

!   Prof. Larson's modified mixing_lhs 
        call mixing_lhs("thlm",.true., dt, taum, wmm, wmt, wp2, wp3, Sc,
     .                   C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc,
     .                   aa, dd, al, indx, implemented, isValid )
        if ( .not. isValid ) then
           print *, "mixing_lhs failed"
           deallocate( aa, dd, al, indx )
           return
        endif


!       Solve thlm / wpthlp
        call mixing_solve( "thlm", .true., dt, taum,
     .                     thlm, wpthlp, thlpthvp, thlm_forcing,
     .                     rhot, rhom, exner, rcm, Ncm,  ! Brian
     .                     cond_rrm, auto_rrm, accr_rrm, ! Brian
     .                     Fcsed, sed_rcm,               ! Brian
     .                     cloud_sed, kk_rain,           ! Brian
     .                     aa, dd, al, indx, C7_Skw_fnc )

!       Solve sclrm / wpsclrp
        if ( present( sclrm )     .and. present( wpsclrp ) .and.
     .       present( sclrpthvp ) .and. present( sclrm_forcing ) ) then

!         Prof. Larson's modified mixing_lhs 
          call mixing_lhs("rtm",.true., dt, taum, wmm, wmt, wp2, wp3, 
     .                    Sc, C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc,
     .                    aa, dd, al, indx, implemented, isValid )
          if ( .not. isValid ) then
            print *, "mixing_lhs failed"
            deallocate( aa, dd, al, indx )
            return
          endif


          call mixing_solve( "sclrm_rtm", .true., dt, taum,
     .                       sclrm(:,SCLR_RT), wpsclrp(:,SCLR_RT),
     .                       sclrpthvp(:,SCLR_RT), 
     .                       sclrm_forcing(:,SCLR_RT),
     .                       rhot, rhom, exner, rcm, Ncm,  ! Brian
     .                       cond_rrm, auto_rrm, accr_rrm, ! Brian
     .                       Fcsed, sed_rcm,               ! Brian
     .                       cloud_sed, kk_rain,           ! Brian
     .                       aa, dd, al, indx, C7_Skw_fnc )

!         Prof. Larson's modified mixing_lhs 
          call mixing_lhs("thlm",.true., dt, taum, wmm, wmt, wp2, wp3,
     .                    Sc, C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc,
     .                    aa, dd, al, indx, implemented, isValid )
          if ( .not. isValid ) then
            print *, "mixing_lhs failed"
            deallocate( aa, dd, al, indx )
            return
          endif

          call mixing_solve( "sclrm_thlm", .true., dt, taum,
     .                       sclrm(:,SCLR_THETA), wpsclrp(:,SCLR_THETA),
     .                       sclrpthvp(:,SCLR_THETA), 
     .                       sclrm_forcing(:,SCLR_THETA),
     .                       rhot, rhom, exner, rcm, Ncm,  ! Brian
     .                       cond_rrm, auto_rrm, accr_rrm, ! Brian
     .                       Fcsed, sed_rcm,               ! Brian
     .                       cloud_sed, kk_rain,           ! Brian
     .                       aa, dd, al, indx, C7_Skw_fnc )
        endif
 
!       Free memory

        deallocate( aa, dd, al, indx )
        deallocate( Skw_arr )
        deallocate( C6rt_Skw_fnc )
        deallocate( C6thl_Skw_fnc )
        deallocate( C7_Skw_fnc )

        return
        end subroutine timestep_mixing

!------------------------------------------------------------------------

        subroutine mixing_lhs( solve_type, liter, dt, taum, wmm, wmt,
     .                         wp2, wp3, Sc, 
     .                         C6rt_Skw_fnc, C6thl_Skw_fnc,C7_Skw_fnc,
     .                         aa, dd, al, indx, implemented, isValid )

        use constants
        use grid_class
#ifdef STATS
        use statistics
#endif STATS

        implicit none

!       Input variables

! Prof. Larson's modification for split C6/C2
        character(len=*), intent(in) :: solve_type  ! For new mixing scheme
        logical, intent(in) :: liter

        real, intent(in) :: dt
        real, intent(in) :: taum(1:gr%nnzp)
        real, intent(in) :: wmm(1:gr%nnzp), wmt(1:gr%nnzp) 
        real, intent(in) :: wp2(1:gr%nnzp), wp3(1:gr%nnzp)
        real, intent(in) :: Sc(1:gr%nnzp)

!       Output variables

        real, intent(out)      :: aa(1:2*gr%nnzp,1:5)
        real, intent(out)      :: dd(1:2*gr%nnzp)
        real, intent(out)      :: al(1:2*gr%nnzp,1:2)
        integer, intent(out)   :: indx(1:2*gr%nnzp)
        logical, intent(inout) :: implemented
        logical, intent(inout) :: isValid

        real, intent(in) :: C6rt_Skw_fnc(1:gr%nnzp)
        real, intent(in) :: C6thl_Skw_fnc(1:gr%nnzp)
        real, intent(in) :: C7_Skw_fnc(1:gr%nnzp)

!       Local constants

!        real, parameter    :: eps = 1.e-10   ! in module constants

!       Internal

        real d
        integer k, km1, kp1
        integer k_xm, k_wpxp
        real a1(1:gr%nnzp)

!       Define a1

        a1 = 1.0 / (1.0-Sc)

!       Set up matrix aa

        aa = 0.0

!        do k=2,gr%nnzp
        do k=2,gr%nnzp-1

!         Define indices

          km1 = max(k-1,1)
          kp1 = min(k+1,gr%nnzp)

          k_xm   = 2*k - 1
          k_wpxp = 2*k

!         xm: mean vertical advection

          IF ( .NOT. implemented ) THEN

          ! Thermodynamic subdiagonal [ x xm(k-1,<t+1>) ]
          aa(k_xm,3-2) = - wmt(k) * (1.0/2.0) * gr%dzt(k)
          ! Thermodynamic superdiagonal [ x xm(k+1,<t+1>) ]
          aa(k_xm,3+2) = + wmt(k) * (1.0/2.0) * gr%dzt(k)

          ELSE

          ! Thermodynamic subdiagonal [ x xm(k-1,<t+1>) ]
          aa(k_xm,3-2) = 0.0
          ! Thermodynamic superdiagonal [ x xm(k+1,<t+1>) ]
          aa(k_xm,3+2) = 0.0

          ENDIF

!         xm: turbulent transport

          ! Momentum subdiagonal [ x wpxp(k-1,<t+1>) ]
          aa(k_xm,3-1) = - gr%dzt(k)
          ! Momentum superdiagonal [ x wpxp(k,<t+1>) ]
          aa(k_xm,3+1) = + gr%dzt(k)

!         xm: time tendency

          ! Thermodynamic main diagonal [ x xm(k,<t+1>) ]
          aa(k_xm,3) = 1.0 / dt

!         wpxp: mean vertical advection, turbulent transport, time tendency, 
!         damping and diffusion

          ! Momentum subdiagonal [ x wpxp(k-1,<t+1>) ]
          aa(k_wpxp,3-2)
     .    = - wmm(k) * (1.0/2.0) * gr%dzm(k)
     .      - gr%dzm(k)
     .        * (1.0/2.0) * ( a1(k) + a1(km1) )
     .        * ( wp3(k) / MAX( wp2(km1) + wp2(k), 2.0*eps ) )
     .      - nu6 * gr%dzm(k) * gr%dzt(k)      

          ! Momentum main diagonal [ x wpxp(k,<t+1>) ]
          select case ( trim(solve_type) )
          case ("rtm")
            aa(k_wpxp,3)
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
            if (liter) aa(k_wpxp,3) = aa(k_wpxp,3) + 1.0 / dt
          case ("thlm")
            aa(k_wpxp,3)
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
            if (liter) aa(k_wpxp,3) = aa(k_wpxp,3) + 1.0 / dt
          end select

          ! Momentum superdiagonal [ x wpxp(k+1,<t+1>) ]
          aa(k_wpxp,3+2)
     .    = + wmm(k) * (1.0/2.0) * gr%dzm(k)
     .      + gr%dzm(k)
     .        * (1.0/2.0) * ( a1(kp1) + a1(k) )
     .        * ( wp3(kp1) / MAX( wp2(k) + wp2(kp1), 2.0*eps ) )
     .      - nu6 * gr%dzm(k) * gr%dzt(kp1)

!         wpxp: production

          ! Thermodynamic subdiagonal [ x xm(k,<t+1>) ]
          aa(k_wpxp,3-1)
     .    = - wp2(k) * gr%dzm(k)

          ! Thermodynamic superdiagonal [ x xm(k+1,<t+1>) ]
          aa(k_wpxp,3+1)
     .    = + wp2(k) * gr%dzm(k)

#ifdef STATS
          if (lstats_samp) then

!          Implicit contributions for wprtp or wpthlp
!          qt changed to rt to avoid confusion -dschanen 2 Aug 05

           if ( iwprtp_ma > 0 .or. iwpthlp_ma > 0 ) then
             zmscr01(k) = + wmm(k) * (1.0/2.0) * gr%dzm(k)
             zmscr02(k) = - wmm(k) * (1.0/2.0) * gr%dzm(k)
           endif

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
           endif

           if ( iwprtp_tp > 0 .or. iwpthlp_tp > 0 ) then
             zmscr06(k) = + wp2(k) * gr%dzm(k)
             zmscr07(k) = - wp2(k) * gr%dzm(k)
           endif

           if ( iwprtp_ac > 0 .or. iwpthlp_ac > 0 ) then
             zmscr08(k) = - 1.0 * gr%dzm(k) * ( wmt(kp1) - wmt(k) )
           endif

           if ( iwprtp_pr1 > 0 .or. iwpthlp_pr1 > 0 ) then
             select case (trim( solve_type ) )
             case("rtm")
!               zmscr09(k) = - C6rt / taum(k)
               zmscr09(k) = - C6rt_Skw_fnc(k) / taum(k)
             case("thlm")
!               zmscr09(k) = - C6thl / taum(k)
               zmscr09(k) = - C6thl_Skw_fnc(k) / taum(k)
             case default
               print *, "Unknown solve_type:", solve_type
               stop "invalid mixing_lhs call"
             end select
           endif

           if ( iwprtp_pr2 > 0 .or. iwpthlp_pr2 > 0 ) then
             zmscr10(k) = + C7_Skw_fnc(k) 
     .                      * gr%dzm(k) * ( wmt(kp1) - wmt(k) )
           endif

           if ( iwprtp_dp1 > 0 .or. iwpthlp_dp1 > 0 ) then
             zmscr11(k) = + nu6 * gr%dzm(k) * gr%dzt(k)
             zmscr12(k) = - nu6 * gr%dzm(k) 
     .                            * ( gr%dzt(kp1) + gr%dzt(k) )
             zmscr13(k) = + nu6 * gr%dzm(k) * gr%dzt(kp1)
           endif

          endif
#endif
        enddo

!       Boundary condition for the flux

        k      = 1
        k_xm   = 2*k - 1
        k_wpxp = 2*k
        aa(k_xm,3-1)   = 0.0
        aa(k_xm,3)     = 1.0
        aa(k_xm,3+1)   = 0.0
        aa(k_wpxp,3-2) = 0.0
        aa(k_wpxp,3)   = 1.0
        aa(k_wpxp,3+2) = 0.0

        k      = gr%nnzp
        k_xm   = 2*k - 1
        k_wpxp = 2*k
        aa(k_xm,3-1)   = 0.0
        aa(k_xm,3)     = 1.0
        aa(k_xm,3+1)   = 0.0
        aa(k_wpxp,3-2) = 0.0
        aa(k_wpxp,3)   = 1.0
        aa(k_wpxp,3+2) = 0.0

!       Perform LU decomposition

        call bandec( aa, 2*gr%nnzp, 2, 2, 2*gr%nnzp, 5, al, 2, indx, d, 
     .               isValid )

        return
        end subroutine mixing_lhs

!------------------------------------------------------------------------
        
        subroutine mixing_solve( solve_type, liter, dt, taum,
     .                           xm, wpxp, xpthvp, xm_forcing,
     .                           rhot, rhom, exner, rcm, Ncm,  ! Brian
     .                           cond_rrm, auto_rrm, accr_rrm, ! Brian
     .                           Fcsed, sed_rcm,               ! Brian
     .                           cloud_sed, kk_rain,           ! Brian
     .                           aa, dd, al, indx, C7_Skw_fnc )

        use constants
        use grid_class
#ifdef STATS
        use statistics
#endif STATS

        implicit none

!       Input/Output

        character(len=*), intent(in) :: solve_type  ! For new mixing scheme

        logical, intent(in) :: liter
        real, intent(in)    :: dt
        real, intent(in)    :: taum(1:gr%nnzp)
        real, intent(inout) :: xm(1:gr%nnzp), wpxp(1:gr%nnzp)
        real, intent(in)    :: xpthvp(1:gr%nnzp)
        real, intent(in)    :: xm_forcing(1:gr%nnzp)

        real, intent(in)    :: rhot(1:gr%nnzp)      ! Brian
        real, intent(in)    :: rhom(1:gr%nnzp)      ! Brian
        real, intent(in)    :: exner(1:gr%nnzp)     ! Brian
        real, intent(in)    :: rcm(1:gr%nnzp)       ! Brian
        real, intent(in)    :: Ncm(1:gr%nnzp)       ! Brian
        real, intent(in)    :: cond_rrm(1:gr%nnzp)  ! Brian
        real, intent(in)    :: auto_rrm(1:gr%nnzp)  ! Brian
        real, intent(in)    :: accr_rrm(1:gr%nnzp)  ! Brian
        real, intent(out)   :: Fcsed(1:gr%nnzp)     ! Brian
        real, intent(out)   :: sed_rcm(1:gr%nnzp)   ! Brian
        logical, intent(in) :: cloud_sed            ! Brian
        logical, intent(in) :: kk_rain              ! Brian

        real, intent(in)    :: aa(1:2*gr%nnzp,1:5)
        real, intent(inout) :: dd(1:2*gr%nnzp)
        real, intent(in)    :: al(1:2*gr%nnzp,1:2)
        integer, intent(in) :: indx(1:2*gr%nnzp)

        real, intent(in) :: C7_Skw_fnc(1:gr%nnzp)

!       Internal

        integer k, km1, kp1
        integer k_xm, k_wpxp

#ifdef STATS
        integer
     .     iwpxp_bt,
     .     iwpxp_ma,
     .     iwpxp_ta,
     .     iwpxp_tp,
     .     iwpxp_ac,
     .     iwpxp_bp,
     .     iwpxp_pr1,
     .     iwpxp_pr2,
     .     iwpxp_pr3,
     .     iwpxp_dp1
        select case ( trim(solve_type) )
        case ( "rtm" )
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
        case ( "thlm" )
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
#endif STATS

#ifdef STATS
        if (lstats_samp) then

!         Explicit contributions for wpxp

          if ( iwpxp_bt > 0 ) then
            zm%x(:,iwpxp_bt) = zm%x(:,iwpxp_bt) - wpxp / dt
          endif

          if ( iwpxp_bp > 0 ) then
            zm%x(:,iwpxp_bp) = zm%x(:,iwpxp_bp)
     .      + grav/T0 * 1.0 * xpthvp(:)
            zm%n(:,iwpxp_bp) = zm%n(:,iwpxp_bp) + 1
          endif

          if ( iwpxp_pr3 > 0 ) then
            zm%x(:,iwpxp_pr3) = zm%x(:,iwpxp_pr3)
     .      - grav/T0 * C7_Skw_fnc * xpthvp(:)
            zm%n(:,iwpxp_pr3) = zm%n(:,iwpxp_pr3) + 1
          endif
        endif
#endif

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
     .                  * EXP( 5.0*( (LOG(1.5))**(2.0) ) )
              ELSE
                 Fcsed(k) = 0.0
              ENDIF

           ENDDO

           ! Boundary conditions.
           Fcsed(1)       = 0.0
           Fcsed(gr%nnzp) = 0.0

           ! Find drc/dt due to cloud water sedimentation flux.
           ! This value is defined on thermodynamic levels.
           ! Fcsed units:  kg (liquid) / [ m^2 * s ]
           ! Multiply by Lv for units of W / m^2.
           ! sed_rcm units:  [ kg (liquid) / kg (air) ] / s
           sed_rcm = (1.0/rhot) * ddzm( Fcsed )

        ENDIF


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
             case ( "sclrm_rtm" ) 
                dd(k_xm) = dd(k_xm) + sed_rcm(k)
             case ( "sclrm_thlm" )
                dd(k_xm) = dd(k_xm) - ( Lv / (Cp*exner(k)) )
     .                               * sed_rcm(k)
             end select

          ENDIF

          ! + xm_microphysics -- rain effects.
          ! Code addition by Brian for rain.
          IF ( kk_rain ) THEN

             select case(trim( solve_type ))
             case ( "rtm" ) 
                dd(k_xm) = dd(k_xm)
     .                     - cond_rrm(k) - auto_rrm(k) - accr_rrm(k)    ! rtm only; Brian
             case ( "thlm" )
                dd(k_xm) = dd(k_xm) - ( Lv / (Cp*exner(k)) ) *
     .                    ( -cond_rrm(k) - auto_rrm(k) - accr_rrm(k) )  ! thlm only; Brian
             case ( "sclrm_rtm" ) 
                dd(k_xm) = dd(k_xm) 
     .                     - cond_rrm(k) - auto_rrm(k) - accr_rrm(k)    ! sclrm = rtm
             case ( "sclrm_thlm" )
                dd(k_xm) = dd(k_xm) - ( Lv / (Cp*exner(k)) ) *
     .                    ( -cond_rrm(k) - auto_rrm(k) - accr_rrm(k) )  ! sclrm = thlm
             end select

          ENDIF

!         wpxp: time tendency and buoyancy

          ! Right-hand side.
          dd(k_wpxp)
     .    = grav/T0 * (1.0-C7_Skw_fnc(k)) * xpthvp(k)
          if (liter) dd(k_wpxp) = dd(k_wpxp) + wpxp(k) / dt

        enddo

!       Boundary condition for the flux

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

!       Solve system

        call banbks( aa, 2 * gr%nnzp, 2, 2, 
     .               2 * gr%nnzp, 5, al, 2, indx, dd )

!       Copy result into output arrays

        do k=1,gr%nnzp

          k_xm   = 2 * k - 1
          k_wpxp = 2 * k

          xm(k)   = dd(k_xm)
          wpxp(k) = dd(k_wpxp)

        enddo

!       Boundary condition on xm

!        xm(1) = 2. * xm(2) - xm(3)
!        xm(gr%nnzp) = 2. * xm(gr%nnzp-1) - xm(gr%nnzp-2)
        xm(1) = xm(2)
!        xm(gr%nnzp) = xm(gr%nnzp-1)

#ifdef STATS
          if (lstats_samp) then

!          Finalize implicit contributions for wpxp

           do k = 1, gr%nnzp
             km1 = max(k-1,1)
             kp1 = min(k+1,gr%nnzp)

             if ( iwpxp_ma > 0 ) then
               zm%x(k,iwpxp_ma) = zm%x(k,iwpxp_ma)
     .         + zmscr01(k) * wpxp(km1)
     .         + zmscr02(k) * wpxp(kp1)
               zm%n(k,iwpxp_ma) = zm%n(k,iwpxp_ma) + 1
             endif

             if ( iwpxp_ta > 0 ) then
               zm%x(k,iwpxp_ta) = zm%x(k,iwpxp_ta)
     .         + zmscr03(k) * wpxp(km1)
     .         + zmscr04(k) * wpxp(k)
     .         + zmscr05(k) * wpxp(kp1)
               zm%n(k,iwpxp_ta) = zm%n(k,iwpxp_ta) + 1
             endif

             if ( iwpxp_tp > 0 ) then
               zm%x(k,iwpxp_tp) = zm%x(k,iwpxp_tp)
     .         + zmscr06(k) * xm(k)
     .         + zmscr07(k) * xm(kp1)
               zm%n(k,iwpxp_tp) = zm%n(k,iwpxp_tp) + 1
             endif

             if ( iwpxp_ac > 0 ) then
               zm%x(k,iwpxp_ac) = zm%x(k,iwpxp_ac)
     .         + zmscr08(k) * wpxp(k)
               zm%n(k,iwpxp_ac) = zm%n(k,iwpxp_ac) + 1
             endif

             if ( iwpxp_pr1 > 0 ) then
               zm%x(k,iwpxp_pr1) = zm%x(k,iwpxp_pr1)
     .         + zmscr09(k) * wpxp(k)
               zm%n(k,iwpxp_pr1) = zm%n(k,iwpxp_pr1) + 1
             endif

             if ( iwpxp_pr2 > 0 ) then
               zm%x(k,iwpxp_pr2) = zm%x(k,iwpxp_pr2)
     .         + zmscr10(k) * wpxp(k)
               zm%n(k,iwpxp_pr2) = zm%n(k,iwpxp_pr2) + 1
             endif

             if ( iwpxp_dp1 > 0 ) then
               zm%x(k,iwpxp_dp1) = zm%x(k,iwpxp_dp1)
     .         + zmscr11(k) * wpxp(km1)
     .         + zmscr12(k) * wpxp(k)
     .         + zmscr13(k) * wpxp(kp1)
               zm%n(k,iwpxp_dp1) = zm%n(k,iwpxp_dp1) + 1
             endif

           enddo

          endif
#endif

#ifdef STATS
        if (lstats_samp) then

          if ( iwpxp_bt > 0 ) then
            zm%x(:,iwpxp_bt) = zm%x(:,iwpxp_bt) + wpxp / dt
            zm%n(:,iwpxp_bt) = zm%n(:,iwpxp_bt) + 1
          endif

        endif
#endif STATS

        return
        end subroutine mixing_solve
      end module mixing
