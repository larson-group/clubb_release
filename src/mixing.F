!-----------------------------------------------------------------------
! $Id: mixing.F,v 1.6 2005-07-09 17:50:52 griffinb Exp $

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
     .             ( dt, taum, wmm, wmt, wp2, wp3, Sc, exner,
     .               rtm, wprtp, 
     .               rtpthvp, rtm_forcing,
     .               cond_rrm, auto_rrm, accr_rrm,
     .               thlm, wpthlp, thlpthvp, thlm_forcing, isValid,
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
!       Additions for DYCOMS_2
        real, intent(in)    :: exner(1:gr%nnzp)    
        real, intent(in)    :: cond_rrm(1:gr%nnzp)
        real, intent(in)    :: auto_rrm(1:gr%nnzp)
        real, intent(in)    :: accr_rrm(1:gr%nnzp)
!       End Additions for DYCOMS_2
        real, intent(inout) :: rtm(1:gr%nnzp)
        real, intent(inout) :: wprtp(1:gr%nnzp)
        real, intent(in)    :: rtpthvp(1:gr%nnzp)
        real, intent(in)    :: rtm_forcing(1:gr%nnzp)
        real, intent(inout) :: thlm(1:gr%nnzp)
        real, intent(inout) :: wpthlp(1:gr%nnzp)
        real, intent(in)    :: thlpthvp(1:gr%nnzp)
        real, intent(in)    :: thlm_forcing(1:gr%nnzp)
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

!       Allocate memory
        allocate( aa(1:2*gr%nnzp,1:5) )
        allocate( dd(1:2*gr%nnzp) )
        allocate( al(1:2*gr%nnzp,1:2) )
        allocate( indx(1:2*gr%nnzp) )

!       Setup and decompose matrix

        call mixing_lhs( .true., dt, taum, wmm, wmt, wp2, wp3, Sc,
     .                   aa, dd, al, indx, isValid )
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
        call mixing_solve( "rtm", .true., dt, taum, exner,
     .                     rtm, wprtp, rtpthvp, rtm_forcing,
     .                     cond_rrm, auto_rrm, accr_rrm, ! Brian DYCOMS_2
     .                     aa, dd, al, indx )  

!       Solve thlm / wpthlp
        call mixing_solve( "thlm", .true., dt, taum, exner,
     .                     thlm, wpthlp, thlpthvp, thlm_forcing,
     .                     cond_rrm, auto_rrm, accr_rrm, ! Brian DYCOMS_2
     .                     aa, dd, al, indx )  

!       Solve sclrm / wpsclrp
        if ( present( sclrm )     .and. present( wpsclrp ) .and.
     .       present( sclrpthvp ) .and. present( sclrm_forcing ) ) then
          do i=1, sclrm_dimension
            call mixing_solve( "sclrm", .true., dt, taum, exner,
     .                         sclrm(i,:), wpsclrp(i,:),
     .                         sclrpthvp(i,:), sclrm_forcing(i,:),
     .                         cond_rrm, auto_rrm, accr_rrm,
     .                         aa, dd, al, indx )
          enddo
        endif
!        print *, "after"
!        print *, "rt=",  sum( sclrm )   - sum ( rtm ),
!     .           "wp=",  sum( wprtp )   - sum ( wpsclrp ),
!     .           "thvp=",sum( rtpthvp ) - sum ( sclrpthvp )

 
!       Free memory

        deallocate( aa, dd, al, indx )

        return
        end subroutine timestep_mixing

!------------------------------------------------------------------------

        subroutine mixing_lhs( liter, dt, taum, wmm, wmt, wp2, wp3, Sc, 
     .                         aa, dd, al, indx, isValid )

        use constants
        use grid_class
#ifdef STATS
        use statistics
#endif STATS

        implicit none

!       Input variables

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
        logical, intent(inout) :: isValid


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

          aa(k_xm,3-2) = - 0.5 * wmt(k) * gr%dzt(k)
          aa(k_xm,3+2) = + 0.5 * wmt(k) * gr%dzt(k)

!         xm: turbulent transport

          aa(k_xm,3-1) = - gr%dzt(k)
          aa(k_xm,3+1) = + gr%dzt(k)

!         xm: time tendency

          aa(k_xm,3) = 1.0 / dt

!         wpxp: mean vertical advection, turbulent transport, time tendency, 
!         damping and diffusion

          aa(k_wpxp,3-2)
     .    = - 0.5 * wmm(k) * gr%dzm(k)
     .      - a1(k) * gr%dzm(k)
     .        * wp3(k)/max(wp2(km1)+wp2(k),2.0*eps)
     .      - nu6 * gr%dzm(k) * gr%dzt(k)      

          aa(k_wpxp,3)
     .    = + C6 / taum(k)
     .      + a1(k) * gr%dzm(k)
     .        * (  wp3(kp1)/max(wp2(k)+wp2(kp1),2.0*eps)
     .           - wp3(k)  /max(wp2(km1)+wp2(k),2.0*eps) )
     .      + (1.0-C7) * gr%dzm(k) * ( wmt(kp1) - wmt(k) )
     .      + nu6 * gr%dzm(k)*( gr%dzt(kp1) + gr%dzt(k) )
          if (liter) aa(k_wpxp,3) = aa(k_wpxp,3) + 1.0 / dt

          aa(k_wpxp,3+2)
     .    = + 0.5 * wmm(k) * gr%dzm(k)
     .      + a1(k) * gr%dzm(k)
     .        * wp3(kp1)/max(wp2(k)+wp2(kp1),2.0*eps)
     .      - nu6 * gr%dzm(k) * gr%dzt(kp1)

!         wpxp: production

          aa(k_wpxp,3-1)
     .    = - wp2(k) * gr%dzm(k)

          aa(k_wpxp,3+1)
     .    = + wp2(k) * gr%dzm(k)

#ifdef STATS
          if (lstats_samp) then

!          Implicit contributions for wpqtp or wpthlp

           if ( iwpqtp_ma > 0 .or. iwpthlp_ma > 0 ) then
             zmscr01(k) = + 0.5 * wmm(k) * gr%dzm(k)
             zmscr02(k) = - 0.5 * wmm(k) * gr%dzm(k)
           endif

           if ( iwpqtp_ta > 0 .or. iwpthlp_ta > 0 ) then
            zmscr03(k) =
     .      + a1(k) * gr%dzm(k)
     .        * wp3(k)/max(wp2(km1)+wp2(k),2.0*eps)
            zmscr04(k) = 
     .      - a1(k) * gr%dzm(k)
     .        * (  wp3(kp1)/max(wp2(k)+wp2(kp1),2.0*eps)
     .           - wp3(k)  /max(wp2(km1)+wp2(k),2.0*eps) )
            zmscr05(k) = 
     .      - a1(k) * gr%dzm(k)
     .        * wp3(kp1)/max(wp2(k)+wp2(kp1),2.0*eps)
           endif

           if ( iwpqtp_tp > 0 .or. iwpthlp_tp > 0 ) then
             zmscr06(k) = + wp2(k) * gr%dzm(k)
             zmscr07(k) = - wp2(k) * gr%dzm(k)
           endif

           if ( iwpqtp_ac > 0 .or. iwpthlp_ac > 0 ) then
             zmscr08(k) = - 1.0 * gr%dzm(k) * ( wmt(kp1) - wmt(k) )
           endif

           if ( iwpqtp_pr1 > 0 .or. iwpthlp_pr1 > 0 ) then
             zmscr09(k) = - C6 / taum(k)
           endif

           if ( iwpqtp_pr2 > 0 .or. iwpthlp_pr2 > 0 ) then
             zmscr10(k) = + C7 * gr%dzm(k) * ( wmt(kp1) - wmt(k) )
           endif

           if ( iwpqtp_dp1 > 0 .or. iwpthlp_dp1 > 0 ) then
             zmscr11(k) = + nu6 * gr%dzm(k) * gr%dzt(k)
             zmscr12(k) = - nu6 * gr%dzm(k) *( gr%dzt(kp1) + gr%dzt(k) )
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
        
        subroutine mixing_solve( solve_type, liter, dt, taum, exner, ! ex: Brian
     .                           xm, wpxp, xpthvp, xm_forcing,
     .                           cond_rrm, auto_rrm, accr_rrm, ! Brian
     .                           aa, dd, al, indx )

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
        real, intent(in)    :: exner(1:gr%nnzp) ! Brian
        real, intent(inout) :: xm(1:gr%nnzp), wpxp(1:gr%nnzp)
        real, intent(in)    :: xpthvp(1:gr%nnzp)
        real, intent(in)    :: xm_forcing(1:gr%nnzp)

        real, intent(in)    :: cond_rrm(1:gr%nnzp)
        real, intent(in)    :: auto_rrm(1:gr%nnzp)
        real, intent(in)    :: accr_rrm(1:gr%nnzp)

        real, intent(in)    :: aa(1:2*gr%nnzp,1:5)
        real, intent(inout) :: dd(1:2*gr%nnzp)
        real, intent(in)    :: al(1:2*gr%nnzp,1:2)
        integer, intent(in) :: indx(1:2*gr%nnzp)

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
          iwpxp_bt  = iwpqtp_bt
          iwpxp_ma  = iwpqtp_ma
          iwpxp_ta  = iwpqtp_ta
          iwpxp_tp  = iwpqtp_tp
          iwpxp_ac  = iwpqtp_ac
          iwpxp_bp  = iwpqtp_bp
          iwpxp_pr1 = iwpqtp_pr1
          iwpxp_pr2 = iwpqtp_pr2
          iwpxp_pr3 = iwpqtp_pr3
          iwpxp_dp1 = iwpqtp_dp1
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
        case default
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
     .      - grav/T0 * C7 * xpthvp(:)
            zm%n(:,iwpxp_pr3) = zm%n(:,iwpxp_pr3) + 1
          endif
        endif
#endif

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

          dd(k_xm) 
     .    = xm(k) / dt
     .    + xm_forcing(k)

          ! Code addition by Brian for rain.
          select case(trim( solve_type ))
          case ( "rtm" ) 
            dd(k_xm) = dd(k_xm)
     .                 - cond_rrm(k) - auto_rrm(k) - accr_rrm(k)  ! rtm only; Brian
          case ( "thlm" )
            dd(k_xm) = dd(k_xm) - ( Lv / (Cp*exner(k)) ) *
     .                 ( -cond_rrm(k) - auto_rrm(k) - accr_rrm(k) )  ! thlm only; Brian
          case ( "sclrm" ) 
            dd(k_xm) = dd(k_xm)
     .                 - cond_rrm(k) - auto_rrm(k) - accr_rrm(k)  ! sclrm = rtm for now
          end select

!         wpxp: time tendency and buoyancy

          dd(k_wpxp)
     .    = grav/T0 * (1.0-C7) * xpthvp(k)
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
!#ifdef STATS
!          if (lstats_samp) then
!
!!          Finalize implicit contributions for wpxp
!
!           do k=1,gr%nnzp
!             km1 = max(k-1,1)
!             kp1 = min(k+1,gr%nnzp)
!
!             if ( iwpxp_ma > 0 ) then
!               zm%x(k,iwpxp_ma) = zm%x(k,iwpxp_ma)
!     .         + zmscr01(k) * wpxp(km1)
!     .         + zmscr02(k) * wpxp(kp1)
!               zm%n(k,iwpxp_ma) = zm%n(k,iwpxp_ma) + 1
!             endif
!
!             if ( iwpxp_ta > 0 ) then
!               zm%x(k,iwpxp_ta) = zm%x(k,iwpxp_ta)
!     .         + zmscr03(k) * wpxp(km1)
!     .         + zmscr04(k) * wpxp(k)
!     .         + zmscr05(k) * wpxp(kp1)
!               zm%n(k,iwpxp_ta) = zm%n(k,iwpxp_ta) + 1
!             endif
!
!             if ( iwpxp_tp > 0 ) then
!               zm%x(k,iwpxp_tp) = zm%x(k,iwpxp_tp)
!     .         + zmscr06(k) * xm(k)
!     .         + zmscr07(k) * xm(kp1)
!               zm%n(k,iwpxp_tp) = zm%n(k,iwpxp_tp) + 1
!             endif
!
!             if ( iwpxp_ac > 0 ) then
!               zm%x(k,iwpxp_ac) = zm%x(k,iwpxp_ac)
!     .         + zmscr08(k) * wpxp(k)
!               zm%n(k,iwpxp_ac) = zm%n(k,iwpxp_ac) + 1
!             endif
!
!             if ( iwpxp_pr1 > 0 ) then
!               zm%x(k,iwpxp_pr1) = zm%x(k,iwpxp_pr1)
!     .         + zmscr09(k) * wpxp(k)
!               zm%n(k,iwpxp_pr1) = zm%n(k,iwpxp_pr1) + 1
!             endif
!
!             if ( iwpqtp_pr2 > 0 ) then
!               zm%x(k,iwpxp_pr2) = zm%x(k,iwpxp_pr2)
!     .         + zmscr10(k) * wpxp(k)
!               zm%n(k,iwpxp_pr2) = zm%n(k,iwpxp_pr2) + 1
!             endif
!
!             if ( iwpxp_dp1 > 0 ) then
!               zm%x(k,iwpxp_dp1) = zm%x(k,iwpxp_dp1)
!     .         + zmscr11(k) * wpxp(km1)
!     .         + zmscr12(k) * wpxp(k)
!     .         + zmscr13(k) * wpxp(kp1)
!               zm%n(k,iwpxp_dp1) = zm%n(k,iwpxp_dp1) + 1
!             endif
!
!           enddo
!
!          endif
!#endif STATS

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
