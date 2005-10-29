!------------------------------------------------------------------------
! $Id: wp23.F,v 1.9 2005-10-29 21:25:07 griffinb Exp $

! SUBROUTINE timestep_wp23( )
!  Advance wp2 and wp3 one timestep
!------------------------------------------------------------------------
        subroutine timestep_wp23
     .             (dt, Lscale, taum, taut, Sc, a, Khm, wmt, wmm,
     .              um, upwp, vm, vpwp, wpthvp, wp2thvp,
     .              wp2, wp3, isValid)

        use grid_class
        use constants
#ifdef STATS
        use statistics
#endif STATS

        implicit none

!       Input/Output

        real, intent(in)    :: dt
        real, intent(in)    :: Lscale(1:gr%nnzp)
        real, intent(in)    :: taum(1:gr%nnzp)
        real, intent(in)    :: taut(1:gr%nnzp)
        real, intent(in)    :: Sc(1:gr%nnzp)
        real, intent(in)    :: a(1:gr%nnzp)
        real, intent(in)    :: Khm(1:gr%nnzp)
        real, intent(in)    :: wmt(1:gr%nnzp)
        real, intent(in)    :: wmm(1:gr%nnzp)
        real, intent(in)    :: um(1:gr%nnzp)
        real, intent(in)    :: upwp(1:gr%nnzp)
        real, intent(in)    :: vm(1:gr%nnzp)
        real, intent(in)    :: vpwp(1:gr%nnzp)
        real, intent(in)    :: wpthvp(1:gr%nnzp)
        real, intent(in)    :: wp2thvp(1:gr%nnzp)
        real, intent(inout) :: wp2(1:gr%nnzp)
        real, intent(inout) :: wp3(1:gr%nnzp)
!        real, intent(inout) :: em(1:gr%nnzp)

        logical, intent(inout):: isValid

!       Internal

        real, allocatable :: wp2t(:), wp3t(:), tauw3t(:), Lm(:)
      
        integer k
        real aa, Skw

! Internal variables for C11 function, Vince Larson 13 Mar 2005
! Brian added C1 function.
        integer km1, kp1
        real, allocatable :: Skw_arr_m(:), C1_Skw_fnc(:)
        real, allocatable :: Skw_arr_t(:), C11_Skw_fnc(:)
! End Vince Larson's addition.

!       Allocate memory

        allocate( tauw3t(1:gr%nnzp) )
        allocate( wp2t(1:gr%nnzp) )
        allocate( wp3t(1:gr%nnzp) )
        allocate( Lm(1:gr%nnzp) )

! Allocate variables for C11 function, Vince Larson 13 Mar 2005
! Brian added C1 function
        allocate( Skw_arr_m(1:gr%nnzp) )
        allocate( C1_Skw_fnc(1:gr%nnzp) )
        allocate( Skw_arr_t(1:gr%nnzp) )
        allocate( C11_Skw_fnc(1:gr%nnzp) )
! End Vince Larson's addition.

!       Define L on momentum levels

        Lm = zt2zm( Lscale )
        Lm(1) = Lm(2)
        Lm(gr%nnzp) = Lm(gr%nnzp-1)

!       RHSs

        wp2t = - wmm * ddzt( zm2zt( wp2 ) )
     .         + 2.0*(1.0-C5) * grav * alpha * wpthvp
     .         - 2.0*(1.0-C5) * wp2 * ddzt( wmt )
     .         + C5 * 2./3. * grav * alpha * wpthvp
     .         - C5 * 2./3. * upwp * ddzt( um )
     .         - C5 * 2./3. * vpwp * ddzt( vm )

#ifdef STATS
       if (lstats_samp) then

         if ( iwp2_ma > 0 ) then
           zm%x(:,iwp2_ma) = zm%x(:,iwp2_ma) 
     .     - wmm * ddzt( zm2zt( wp2 ) )
           zm%n(:,iwp2_ma) = zm%n(:,iwp2_ma) + 1
         endif

         if ( iwp2_ac > 0 ) then
           zm%x(:,iwp2_ac) = zm%x(:,iwp2_ac) 
     .     - 2.0 * wp2 * ddzt( wmt )
           zm%n(:,iwp2_ac) = zm%n(:,iwp2_ac) + 1
         endif

         if ( iwp2_bp > 0 ) then
           zm%x(:,iwp2_bp) = zm%x(:,iwp2_bp) 
     .     + 2.0* grav * alpha * wpthvp
           zm%n(:,iwp2_bp) = zm%n(:,iwp2_bp) + 1
         endif

         if ( iwp2_pr2 > 0 ) then
           zm%x(:,iwp2_pr2) = zm%x(:,iwp2_pr2) 
     .     - C5 * ( - 2.0 * wp2 * ddzt( wmt ) 
     .              + 2.0 * grav * alpha * wpthvp
     .            )
           zm%n(:,iwp2_pr2) = zm%n(:,iwp2_pr2) + 1
         endif

         if ( iwp2_pr3 > 0 ) then
           zm%x(:,iwp2_pr3) = zm%x(:,iwp2_pr3) 
     .     + C5 * 2./3. * grav * alpha * wpthvp
     .     - C5 * 2./3. * upwp * ddzt( um )
     .     - C5 * 2./3. * vpwp * ddzt( vm )
           zm%n(:,iwp2_pr3) = zm%n(:,iwp2_pr3) + 1
         endif

       endif
#endif STATS

! Vince Larson added code to make C11 function of Skw. 13 Mar 2005
! If this code is used, C11 is no longer relevant, i.e. constants
!    are hardwired.

! Calculate skewness array
        do k=1,gr%nnzp

          km1 = max(k-1,1)
          kp1 = min(k+1,gr%nnzp)

          Skw_arr_t(k) = wp3(k) /
     .                      ( 0.5*max(wtol**2,wp2(km1)+wp2(k)) )**1.5
          Skw_arr_m(k) = 0.5*(wp3(kp1)+wp3(k)) /
     .                      ( MAX(wtol**2.0, wp2(k))**1.5 )

        enddo

        C11_Skw_fnc = C11b + (C11-C11b)
     .                      *EXP( -0.5 * (Skw_arr_t/C11c)**2.0 )

        C1_Skw_fnc = C1b + (C1-C1b)
     .                    *EXP( -0.5 * (Skw_arr_m/C1c)**2.0 )

!        C11_Skw_fnc = C11
!        C1_Skw_fnc = C1

! Old version with C11 constant
!        wp3t = - wmt * ddzm( zt2zm( wp3 ) )
!     .         - 3.0*(1.0-C11) * wp3 * ddzm( wmm )
!     .         + 3.0*(1.0-C11) * grav * alpha * wp2thvp

! New version with C11 dependent on Skw
        wp3t = - wmt * ddzm( zt2zm( wp3 ) )
     .         - 3.0 * (1.0-C11_Skw_fnc) * wp3 * ddzm( wmm )
     .         + 3.0 * (1.0-C11_Skw_fnc) * grav * alpha * wp2thvp
! End Vince Larson's change.

#ifdef STATS
       if (lstats_samp) then

         if ( iwp3_ma > 0 ) then
           zt%x(:,iwp3_ma) = zt%x(:,iwp3_ma) 
     .     - wmt * ddzm( zt2zm( wp3 ) )
           zt%n(:,iwp3_ma) = zt%n(:,iwp3_ma) + 1
         endif

         if ( iwp3_ac > 0 ) then
           zt%x(:,iwp3_ac) = zt%x(:,iwp3_ac) 
     .     - 3.0 * wp3 * ddzm( wmm )
           zt%n(:,iwp3_ac) = zt%n(:,iwp3_ac) + 1
         endif

         if ( iwp3_bp > 0 ) then
           zt%x(:,iwp3_bp) = zt%x(:,iwp3_bp) 
     .     + 3.0 * grav * alpha * wp2thvp
           zt%n(:,iwp3_bp) = zt%n(:,iwp3_bp) + 1
         endif

         if ( iwp3_pr2 > 0 ) then
           zt%x(:,iwp3_pr2) = zt%x(:,iwp3_pr2) 
!     .     - C11 * ( - 3.0 * wp3 * ddzm( wmm )
     .     - C11_Skw_fnc * ( - 3.0 * wp3 * ddzm( wmm )
     .               + 3.0 * grav * alpha * wp2thvp
     .             )
           zt%n(:,iwp3_pr2) = zt%n(:,iwp3_pr2) + 1
         endif

       endif
#endif STATS

!       Define tauw

!        tauw3t = taut 
!     .           / ( 1.
!     .                   + 3.0 * max( 
!     .                     min(1.-(a-0.01)/(0.05-0.01)
!     .                         ,1.)
!     .                         ,0.) 
!     .                   + 3.0 * max( 
!     .                     min(1.-(a-0.99)/(0.95-0.99)
!     .                         ,1.)
!     .                         ,0.) 
!     .              ) 

!        do k=1,gr%nnzp
!
!          Skw = abs( wp3(k)/max(wp2(k),1.e-8)**1.5 )
!          Skw = min( 5.0, Skw )
!          tauw3t(k) = taut(k) / ( 0.005*Skw**4 + 1.0 )
!
!        enddo

        tauw3t = taut 

!       Solve semi-implicitly

        call wp23_solve( dt, Sc, Lscale, Lm, taum, tauw3t, 0.5*Khm,
     .                   wmt, wmm, wpthvp, wp2thvp,
     .                   wp2t, wp3t, C1_Skw_fnc, 
     .                   wp2, wp3, isValid )
        
        if ( .not. isValid ) then
          deallocate( tauw3t, wp2t, wp3t, Lm)
          deallocate( Skw_arr_t, C11_Skw_fnc )
          deallocate( Skw_arr_m, C1_Skw_fnc )
          return
        endif

!       tke is assumed to be 3/2 of wp2

!        em = 1.5 * wp2   ! moved to hoc_closure_timestep

!       Free memory

!        deallocate( tauw3t, wp2t, wp3t )
        deallocate( tauw3t, wp2t, wp3t, Lm )  ! Brian

! De-allocate variables for C11 function, Vince Larson 13 Mar 2005
! Brian added C1 function.
        deallocate( Skw_arr_t, C11_Skw_fnc )
        deallocate( Skw_arr_m, C1_Skw_fnc )
! End Vince Larson's addition.

        return
        end subroutine timestep_wp23

!------------------------------------------------------------------------
!  subroutine wp23_solve( )
!------------------------------------------------------------------------

        subroutine wp23_solve( dt, Sc, Lscale, Lm, tau1m, tauw3t, Kh,
     .                         wmt, wmm, wpthvp, wp2thvp,
     .                         wp2t, wp3t, C1_Skw_fnc,
     .                         wp2, wp3, isValid)

        use grid_class
        use constants
#ifdef STATS
        use statistics
#endif STATS

        implicit none

!       Input

        real, intent(in) :: dt
        real, intent(in) :: Sc(1:gr%nnzp)
        real, intent(in) :: Lscale(1:gr%nnzp)
        real, intent(in) :: Lm(1:gr%nnzp)
        real, intent(in) :: tau1m(1:gr%nnzp)
        real, intent(in) :: tauw3t(1:gr%nnzp)
        real, intent(in) :: Kh(1:gr%nnzp)
        real, intent(in) :: wmt(1:gr%nnzp)
        real, intent(in) :: wmm(1:gr%nnzp)
        real, intent(in) :: wpthvp(1:gr%nnzp)
        real, intent(in) :: wp2thvp(1:gr%nnzp)
        real, intent(in) :: wp2t(1:gr%nnzp)
        real, intent(in) :: wp3t(1:gr%nnzp)
        real, intent(in) :: C1_Skw_fnc(1:gr%nnzp)

!       Input/Output

        real, intent(inout) :: wp2(1:gr%nnzp)
        real, intent(inout) :: wp3(1:gr%nnzp)
        logical, intent(inout) :: isValid

!       Local constants

!        real, parameter    :: eps = 1.e-10     ! in constants module

!       Internal

        real, allocatable    :: aa(:,:), dd(:), al(:,:)
        integer, allocatable :: indx(:)

!       Internal

        real Skw

        real atmp, ctmp
        real d
        integer k, km1, kp1
        integer k_wp2, k_wp3

        real a1(1:gr%nnzp), a3(1:gr%nnzp)

#ifdef STATS
        if (lstats_samp) then

          if ( iwp2_bt > 0 ) then
            zm%x(:,iwp2_bt) = zm%x(:,iwp2_bt) - wp2 / dt
          endif

          if ( iwp3_bt > 0 ) then
            zt%x(:,iwp3_bt) = zt%x(:,iwp3_bt) - wp3 / dt
          endif

        endif
#endif STATS

!       Allocate memory

        allocate( aa(1:2*gr%nnzp,1:5) )
        allocate( dd(1:2*gr%nnzp) )
        allocate( al(1:2*gr%nnzp,1:2) )
        allocate( indx(1:2*gr%nnzp) )

!       Define a1 and a3

        a1 = 1.0 / (1.0-Sc)
        a3 = 3.0*Sc*Sc 
     .       + 6.0*(1.0-Sc)*Sc 
     .       + (1.0-Sc)*(1.0-Sc) 
     .       - 1.5

!       Set up matrix aa and RHS dd

        aa = 0.0
        dd = 0.0

        do k=1,gr%nnzp

!         Define indices

          km1 = max(k-1,1)
          kp1 = min(k+1,gr%nnzp)

          k_wp3 = 2*k - 1
          k_wp2 = 2*k

          Skw = wp3(k)/(0.5*(wp2(km1)+wp2(k)))**1.5
          skw = min( max( Skw, -4.5 ), 4.5 )

!         wp2

          atmp = - 0.5 * nu1 * gr%dzm(k) * gr%dzt(k)
          ctmp = - 0.5 * nu1 * gr%dzm(k) * gr%dzt(kp1)

          aa(k_wp2,3-2)
     .    = atmp
     .      

          aa(k_wp2,3-1)
     .    = - gr%dzm(k)

          aa(k_wp2,3)
     .    = 1.0 / dt
     .      - atmp - ctmp
!     .      + C1/tau1m(k)
     .      + C1_Skw_fnc(k)/tau1m(k)

          aa(k_wp2,3+1)
     .    = + gr%dzm(k)

          aa(k_wp2,3+2)
     .    = ctmp

          dd(k_wp2)
     .    = - atmp * wp2(km1)
     .      + ( atmp + ctmp + 1.0/dt ) * wp2(k)
     .      - ctmp * wp2(kp1)
     .      + wp2t(k)

#ifdef STATS
          if (lstats_samp) then

!          Explicit contributions for wp2

           if ( iwp2_dp2 > 0 ) then
            zm%x(k,iwp2_dp2) = zm%x(k,iwp2_dp2) 
     .      - atmp * wp2(km1)
     .      + ( atmp + ctmp ) * wp2(k)
     .      - ctmp * wp2(kp1)
           endif

!          Implicit contributions for wp2

           if ( iwp2_dp1 > 0 ) then
!             zmscr01(k) = - C1/tau1m(k)
             zmscr01(k) = - C1_Skw_fnc(k)/tau1m(k)
           endif

           if ( iwp2_dp2 > 0 ) then
             zmscr02(k) = - atmp
             zmscr03(k) = atmp + ctmp
             zmscr04(k) = - ctmp
           endif

           if ( iwp2_ta > 0 ) then
             zmscr05(k) = + gr%dzm(k)
             zmscr06(k) = - gr%dzm(k)
           endif

         endif
#endif STATS

!         wp3
!  Added a new constant C12, a factor of atmp and ctmp 
!  Initially, this new constant will be set to 1.0 -dschanen 9/19/05 
          atmp = -0.5 * C12 * (Kh(km1) + nu8) * gr%dzt(k) * gr%dzm(km1)
          ctmp = -0.5 * C12 * (Kh(k)   + nu8) * gr%dzt(k) * gr%dzm(k)

          aa(k_wp3,3-2)
     .    = atmp
     .      - 0.5 * a1(k) * gr%dzt(k) 
     .            * (wp3(km1)+wp3(k))/max(wp2(km1),eps)

          aa(k_wp3,3-1)
     .    = - 2.0 * a3(k) * gr%dzt(k) * wp2(km1)

          aa(k_wp3,3)
     .    = 1.0 / dt
     .      + C8/tauw3t(k)*(5.0*C8b*Skw**4+1.0)
     .      - atmp - ctmp
     .      + 0.5 * a1(k) * gr%dzt(k) 
     .            * ( (wp3(k)+wp3(kp1))/max(wp2(k),eps)
     .               -(wp3(km1)+wp3(k))/max(wp2(km1),eps) )

          aa(k_wp3,3+1)
     .    = + 2.0 * a3(k) * gr%dzt(k) * wp2(k)

          aa(k_wp3,3+2)
     .    = ctmp
     .      + 0.5 * a1(k) * gr%dzt(k) 
     .            * (wp3(k)+wp3(kp1))/max(wp2(k),eps)

          dd(k_wp3)
     .    = - atmp * wp3(km1)
     .      + ( atmp + ctmp + 1.0/dt ) * wp3(k)
     .      - ctmp * wp3(kp1)
     .      - C8/tauw3t(k)
     .        *(-4.0*C8b*Skw**4)*wp3(k)
     .      + wp3t(k)
     .      + a3(k) * gr%dzt(k) * ( wp2(k)**2 - wp2(km1)**2 )
     .      + 0.25 * a1(k) * gr%dzt(k)
     .            * ( (wp3(k)+wp3(kp1))**2/max(wp2(k),eps)
     .               -(wp3(km1)+wp3(k))**2/max(wp2(km1),eps) )

#ifdef STATS
          if (lstats_samp) then

!          Explicit contributions for wp3

           if ( iwp3_ta > 0 ) then
            zt%x(k,iwp3_ta) = zt%x(k,iwp3_ta) 
     .      + (a3(k)+1.5) * gr%dzt(k) * ( wp2(k)**2 - wp2(km1)**2 )
     .      + 0.25 * a1(k) * gr%dzt(k)
     .            * ( (wp3(k)+wp3(kp1))**2/max(wp2(k),eps)
     .               -(wp3(km1)+wp3(k))**2/max(wp2(km1),eps) )
           endif

           if ( iwp3_tp > 0 ) then
            zt%x(k,iwp3_tp) = zt%x(k,iwp3_tp) 
     .      - 1.5 * gr%dzt(k) * ( wp2(k)**2 - wp2(km1)**2 )
           endif

           if ( iwp3_pr1 > 0 ) then
            zt%x(k,iwp3_pr1) = zt%x(k,iwp3_pr1) 
     .      - C8/tauw3t(k)
     .        *(-4.0*C8b*Skw**4)*wp3(k)
           endif

           if ( iwp3_dp1 > 0 ) then
            zt%x(k,iwp3_dp1) = zt%x(k,iwp3_dp1) 
     .      - atmp * wp3(km1)
     .      + ( atmp + ctmp ) * wp3(k)
     .      - ctmp * wp3(kp1)
           endif

!          Implicit contributions for wp3

           if ( iwp3_ta > 0 ) then
             ztscr05(k) =
     .       + 0.5 * a1(k) * gr%dzt(k) 
     .             * (wp3(km1)+wp3(k))/max(wp2(km1),eps)
             ztscr06(k)
     .       = + 2.0 * (a3(k)+1.5) * gr%dzt(k) * wp2(km1)
             ztscr07(k) =
     .       - 0.5 * a1(k) * gr%dzt(k) 
     .             * ( (wp3(k)+wp3(kp1))/max(wp2(k),eps)
     .                -(wp3(km1)+wp3(k))/max(wp2(km1),eps) )
             ztscr08(k)
     .       = - 2.0 * (a3(k)+1.5) * gr%dzt(k) * wp2(k)
             ztscr09(k) =
     .       - 0.5 * a1(k) * gr%dzt(k) 
     .             * (wp3(k)+wp3(kp1))/max(wp2(k),eps)
           endif

           if ( iwp3_tp > 0 ) then
             ztscr10(k)
     .       = - 2.0 * 1.5 * gr%dzt(k) * wp2(km1)
             ztscr11(k)
     .       = + 2.0 * 1.5 * gr%dzt(k) * wp2(k)
           endif

           if ( iwp3_pr1 > 0 ) then
             ztscr01(k) =
     .       - C8/tauw3t(k)*(5.0*C8b*Skw**4+1.0)    
           endif

           if ( iwp3_dp1 > 0 ) then
             ztscr02(k) = - atmp
             ztscr03(k) = atmp + ctmp
             ztscr04(k) = - ctmp
           endif

          endif
#endif STATS
        enddo

!       Boundary conditions

        k = 1
        k_wp3 = 2*k - 1
        k_wp2 = 2*k

        aa(k_wp2,:) = 0.0
        aa(k_wp2,3) = 1.0
        dd(k_wp2)   = wp2(k)
   
        aa(k_wp3,:) = 0.0
        aa(k_wp3,3) = 1.0
        dd(k_wp3)   = 0.0

        k = gr%nnzp
        k_wp3 = 2*k - 1
        k_wp2 = 2*k

        aa(k_wp2,:) = 0.0
        aa(k_wp2,3) = 1.0
        dd(k_wp2)   = 0.0
   
        aa(k_wp3,:) = 0.0
        aa(k_wp3,3) = 1.0
        dd(k_wp3)   = 0.0

!       Perform LU decomposition

        call bandec(aa,2*gr%nnzp,2,2,2*gr%nnzp,5,al,2,indx,d,isValid)
        if ( .not. isValid ) then
            print *, "LU decomp. failed"
            deallocate( aa, dd, al, indx )
            return
        endif

!       Solve system

        call banbks(aa,2*gr%nnzp,2,2,2*gr%nnzp,5,al,2,indx,dd)

!       Copy result into output arrays and clip

        do k=1,gr%nnzp

          km1 = max(k-1,1)
          kp1 = min(k+1,gr%nnzp)

          k_wp3 = 2*k - 1
          k_wp2 = 2*k

          wp2(k) = dd(k_wp2)
          wp3(k) = dd(k_wp3)

        enddo

#ifdef STATS
        if (lstats_samp) then

!         Finalize implicit contributions for wp2

          do k=1,gr%nnzp
            km1 = max(k-1,1)
            kp1 = min(k+1,gr%nnzp)

            if ( iwp2_dp1 > 0 ) then
             zm%x(k,iwp2_dp1) = zm%x(k,iwp2_dp1)
     .       + zmscr01(k) * wp2(k)
             zm%n(k,iwp2_dp1) = zm%n(k,iwp2_dp1) + 1
            endif

            if ( iwp2_dp2 > 0 ) then
             zm%x(k,iwp2_dp2) = zm%x(k,iwp2_dp2)
     .       + zmscr02(k) * wp2(km1)
     .       + zmscr03(k) * wp2(k)
     .       + zmscr04(k) * wp2(kp1)
             zm%n(k,iwp2_dp2) = zm%n(k,iwp2_dp2) + 1
            endif

            if ( iwp2_ta > 0 ) then
             zm%x(k,iwp2_ta) = zm%x(k,iwp2_ta)
     .       + zmscr05(k) * wp3(k)
     .       + zmscr06(k) * wp3(kp1)
             zm%n(k,iwp2_ta) = zm%n(k,iwp2_ta) + 1
            endif

          enddo

!         Finalize implicit contributions for wp3

          do k=1,gr%nnzp
            km1 = max(k-1,1)
            kp1 = min(k+1,gr%nnzp)

            if ( iwp3_pr1 > 0 ) then
             zt%x(k,iwp3_pr1) = zt%x(k,iwp3_pr1)
     .       + ztscr01(k) * wp3(k)
             zt%n(k,iwp3_pr1) = zt%n(k,iwp3_pr1) + 1
            endif

            if ( iwp3_dp1 > 0 ) then
             zt%x(k,iwp3_dp1) = zt%x(k,iwp3_dp1)
     .       + ztscr02(k) * wp3(km1)
     .       + ztscr03(k) * wp3(k)
     .       + ztscr04(k) * wp3(kp1)
             zt%n(k,iwp3_dp1) = zt%n(k,iwp3_dp1) + 1
            endif

            if ( iwp3_ta > 0 ) then
             zt%x(k,iwp3_ta) = zt%x(k,iwp3_ta)
     .       + ztscr05(k) * wp3(km1)
     .       + ztscr06(k) * wp2(km1)
     .       + ztscr07(k) * wp3(k)
     .       + ztscr08(k) * wp2(k)
     .       + ztscr09(k) * wp3(kp1)
             zt%n(k,iwp3_ta) = zt%n(k,iwp3_ta) + 1
            endif

            if ( iwp3_tp > 0 ) then
             zt%x(k,iwp3_tp) = zt%x(k,iwp3_tp)
     .       + ztscr10(k) * wp2(km1)
     .       + ztscr11(k) * wp2(k)
             zt%n(k,iwp3_tp) = zt%n(k,iwp3_tp) + 1
            endif

          enddo

!         Effect of clipping term

          if ( iwp2_cl > 0 ) then
            zm%x(:,iwp2_cl) = zm%x(:,iwp2_cl) - wp2 / dt
          endif

          if ( iwp3_cl > 0 ) then
            zt%x(:,iwp3_cl) = zt%x(:,iwp3_cl) - wp3 / dt
          endif

        endif
#endif STATS

        ! Clip skewness consistently with a
 
        do k=1,gr%nnzp

         km1 = max(k-1,1)

         atmp = max(0.5*(wp2(km1)+wp2(k)),eps)
         ctmp = atmp**1.5
         if ( wp3(k)/ctmp .gt. 4.5 ) then
            wp3(k) = 4.5 * ctmp
            wp2(k) = atmp
         elseif ( wp3(k)/ctmp .lt. -4.5 ) then
            wp3(k) = -4.5 * ctmp
            wp2(k) = atmp
         endif

         wp2(k) = max( 2./3. * emin, wp2(k) ) 

        enddo

#ifdef STATS
        if (lstats_samp) then
          if ( iwp2_cl > 0 ) then
            zm%x(:,iwp2_cl) = zm%x(:,iwp2_cl) + wp2 / dt
            zm%n(:,iwp2_cl) = zm%n(:,iwp2_cl) + 1
          endif
          if ( iwp3_cl > 0 ) then
            zt%x(:,iwp3_cl) = zt%x(:,iwp3_cl) + wp3 / dt
            zt%n(:,iwp3_cl) = zt%n(:,iwp3_cl) + 1
          endif
        endif
#endif STATS

!       Free memory

        deallocate( aa, dd, al, indx )

#ifdef STATS
        if (lstats_samp) then

          if ( iwp2_bt > 0 ) then
            zm%x(:,iwp2_bt) = zm%x(:,iwp2_bt) + wp2 / dt
            zm%n(:,iwp2_bt) = zm%n(:,iwp2_bt) + 1
          endif

          if ( iwp3_bt > 0 ) then
            zt%x(:,iwp3_bt) = zt%x(:,iwp3_bt) + wp3 / dt
            zt%n(:,iwp3_bt) = zt%n(:,iwp3_bt) + 1
          endif

        endif
#endif STATS

        return
        end subroutine wp23_solve

