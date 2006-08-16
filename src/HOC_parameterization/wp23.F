!------------------------------------------------------------------------
! $Id: wp23.F,v 1.7 2006-08-16 23:39:38 dschanen Exp $

        subroutine timestep_wp23
     .             ( dt, Lscale, taum, taut, Scm, Skwm, Skwt,
     .               a, Kht, wmt, wmm,
     .               um, upwp, vm, vpwp, wpthvp, wp2thvp,
     .               wp2, wp3, isValid )
!       Description:
!       Advance wp2 and wp3 one timestep

!       References:
!       See ``Implict solution for the vertical velocity moments''
!       _Equations for HOC_ (6)
!------------------------------------------------------------------------

        use grid_class
        use constants
        use model_flags
#ifdef STATS
        use statistics
#endif /*STATS*/

        implicit none

        ! External
        external :: wp23_solve

        intrinsic :: exp

        ! Input Variables
        real, intent(in) :: dt ! Timestep       [s]

        real, intent(in), dimension(gr%nnzp) :: 
     .  Lscale,! Mixing Length                  [m]
     .  taum,  ! Tau on moment. grid            [s]
     .  taut,  ! Tau on thermo. grid            [s]
     .  Scm,   ! Sc on moment. grid             [-]
     .  Skwm,  ! Skw on moment. grid            [-]
     .  Skwt,  ! Skw on thermo. grid            [-]
     .  a,     ! PDF parm a (13)                [-]
     .  Kht,   ! Kh Eddy diffusivity on thermo. [m^2/s]
     .  wmt,   ! w wind on t                    [m/s]
     .  wmm,   ! w wind on m                    [m/s]
     .  um,    ! u wind                         [m/s]
     .  upwp,  ! u'w'                           [m^2/s^2]
     .  vm,    ! v wind                         [m/s]
     .  vpwp,  ! u'w'                           [m^2/s^2]
     .  wpthvp,! w'th_v' Buoyancy flux          [(K m)/s]
     .  wp2thvp! w'^2 th_v'                     [(m^2 K)/s^2]

        ! Input/Output
        real, dimension(gr%nnzp), intent(inout) :: 
     .  wp2, ! w'^2                             [m^2/s^2]
     .  wp3  ! w'^3                             [m^3/s^3]

        logical, intent(inout) :: isValid

        ! Local Variables
        real, allocatable :: tauw3t(:), Lm(:)
      
        integer :: k, km1, kp1
        real :: aa

        ! Internal variables for C11 function, Vince Larson 13 Mar 2005
        ! Brian added C1 function.
        real, dimension(:), allocatable :: 
     .  C1_Skw_fnc, 
     .  C11_Skw_fnc
        ! End Vince Larson's addition.

!-----------------------------------------------------------------------

        ! Allocate memory
        allocate( tauw3t(1:gr%nnzp) )
        allocate( Lm(1:gr%nnzp) )

        ! Allocate variables for C11 function, Vince Larson 13 Mar 2005
        ! Brian added C1 function
        allocate( C1_Skw_fnc(1:gr%nnzp) )
        allocate( C11_Skw_fnc(1:gr%nnzp) )
        ! End Vince Larson's addition.

        ! Define L on momentum levels
        Lm = zt2zm( Lscale )
        Lm(1) = Lm(2)
        Lm(gr%nnzp) = Lm(gr%nnzp-1)

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

! Vince Larson added code to make C11 function of Skw. 13 Mar 2005
! If this code is used, C11 is no longer relevant, i.e. constants
!    are hardwired.

! Calculate skewness arrays
!       DO k = 1, gr%nnzp, 1

!         km1 = MAX(k-1,1)
!         kp1 = MIN(k+1,gr%nnzp)

!         Skw_arr_t(k) = wp3(k) /
!    .    ( (1.0/2.0) * MAX( wp2(km1)+wp2(k), wtol**2 ) )**(3.0/2.0)
!         Skw_arr_m(k) = (1.0/2.0) * ( wp3(kp1)+wp3(k) ) /
!    .                   MAX( wp2(k), wtol**2 )**(3.0/2.0)

!       ENDDO

        C11_Skw_fnc(1:gr%nnzp) = 
     .  C11b + (C11-C11b)*EXP( -(1.0/2.0) * (Skwt(1:gr%nnzp)/C11c)**2 )

        C1_Skw_fnc(1:gr%nnzp) = 
     .  C1b + (C1-C1b)*EXP( -(1.0/2.0) * (Skwm(1:gr%nnzp)/C1c)**2 )

!        C11_Skw_fnc = C11
!        C1_Skw_fnc = C1


!       Solve semi-implicitly
!       Note: Kw is defined to be 1/2 of Kht
        call wp23_solve( dt, Scm, Lscale, Lm, taum, tauw3t, 0.5*Kht,
     .                   wmt, wmm, wpthvp, wp2thvp,
     .                   um, upwp, vm, vpwp,
     .                   Skwt, C1_Skw_fnc, C11_Skw_fnc, 
     .                   wp2, wp3, isValid )
        
        if ( .not. isValid ) then
          deallocate( tauw3t, Lm)
          deallocate( C11_Skw_fnc )
          deallocate( C1_Skw_fnc )
          return
        endif

        ! Free memory

!        deallocate( tauw3t )
        deallocate( tauw3t, Lm )  ! Brian

        ! De-allocate variables for C11 function, Vince Larson 13 Mar 2005
        ! Brian added C1 function.
        deallocate( C11_Skw_fnc )
        deallocate( C1_Skw_fnc )
        ! End Vince Larson's addition.

        return
        end subroutine timestep_wp23

!------------------------------------------------------------------------
        subroutine wp23_solve( dt, Scm, Lscale, Lm, tau1m, tauw3t, Kw,
     .                         wmt, wmm, wpthvp, wp2thvp,
     .                         um, upwp, vm, vpwp,
     .                         Skwt, C1_Skw_fnc, C11_Skw_fnc,
     .                         wp2, wp3, isValid)
!       Description:
!       Decompose, and back substitute the matrix for wp2/wp3

!       References:
!       _Equations for HOC_ section 6.3 
!------------------------------------------------------------------------

        use grid_class
        use constants
        use model_flags
#ifdef STATS
        use statistics
#endif /*STATS*/

        implicit none

        ! External
        intrinsic :: max, min, sqrt

        ! Input Variables
        real, intent(in) :: 
     .  dt          ! Timestep                               [s]

        real, intent(in), dimension(gr%nnzp) :: 
     .  Scm,        ! Sc on moment. grid                     [-]
     .  Lscale,     ! Mixing Length                          [m]
     .  Lm,         ! Mixing Length interp. to moment.       [m]
     .  tau1m,      ! Tau on moment. grid                    [s]
     .  tauw3t,     ! Tau on thermo. grid                    [s]
     .  Kw,         ! 1/2 Kht Eddy diffusivity               [m^2/s]
     .  wmt,        ! w wind on t                            [m/s]
     .  wmm,        ! w wind on m                            [m/s]
     .  wpthvp,     ! w'th_v' Buoyancy flux                  [(K m)/s]
     .  wp2thvp,    ! w'^2 th_v'                             [(m^2 K)/s^2]
     .  um,         ! u wind                                 [m/s]
     .  upwp,       ! u'w'                                   [m^2/s^2]
     .  vm,         ! v wind                                 [m/s]
     .  vpwp,       ! u'w'                                   [m^2/s^2]
     .  Skwt,       ! Skw in thermo. grid                    [-]
     .  C1_Skw_fnc, ! C1 parameter with Skw applied          [-]
     .  C11_Skw_fnc ! C11 parameter with Skw applied         [-]

        ! Input/Output Variables
        real, dimension(gr%nnzp), intent(inout) :: 
     .  wp2, ! w'^2                                          [m^2/s^2]
     .  wp3  ! w'^3                                          [m^3/s^3]

        logical, intent(inout) :: isValid

        ! Local Variables
        real, dimension(:,:), allocatable :: 
     .  aa, ! Band matrix
     .  al  ! LU decomposition of the band matrix

        real, dimension(:), allocatable :: 
     .  dd ! RHS of band matrix

        integer, allocatable :: 
     .  indx(:) ! Index of pivots in band matrix decomposition process

        real :: d ! diagnostic variable for LU decomposition

        real, dimension(gr%nnzp) :: 
     .  a1, ! a1; See (20) in _Equations for HOC_             [-]
     .  a3  ! a3; See (22) in _Equations for HOC_             [-]

        real :: atmp, ctmp

        integer :: k, km1, kp1
        integer :: k_wp2, k_wp3


#ifdef STATS
        if (lstats_samp) then

          if ( iwp2_bt > 0 ) then
            zm%x(:,iwp2_bt) = zm%x(:,iwp2_bt) - wp2 / dt
          endif

          if ( iwp3_bt > 0 ) then
            zt%x(:,iwp3_bt) = zt%x(:,iwp3_bt) - wp3 / dt
          endif

        endif
#endif /*STATS*/

!       Allocate memory

        allocate( aa(1:2*gr%nnzp,1:5) )
        allocate( dd(1:2*gr%nnzp) )
        allocate( al(1:2*gr%nnzp,1:2) )
        allocate( indx(1:2*gr%nnzp) )

!       Define a1 and a3 (both are located on momentum levels).

        a1 = 1.0 / (1.0-Scm)
        a3 = 3.0*Scm*Scm 
     .       + 6.0*(1.0-Scm)*Scm 
     .       + (1.0-Scm)*(1.0-Scm) 
     .       - (3.0/2.0)

!       Set up matrix aa and RHS dd

        aa = 0.0
        dd = 0.0

        DO k = 1, gr%nnzp, 1

!         Define indices

          km1 = max(k-1,1)
          kp1 = min(k+1,gr%nnzp)

          k_wp3 = 2*k - 1
          k_wp2 = 2*k

          ! Skewness of w is needed on the thermodynamic levels 
          ! for the wp3 equation
          !Skw = wp3(k) / ( (1.0/2.0)*( wp2(k) + wp2(km1) ) )**(3.0/2.0)
          ! Skewness of w limited to a magnitude of no more than 4.5
          !Skw = MIN( MAX( Skw, -4.5 ), 4.5 )


!         wp2

          ! Diffusional term with Crank-Nicholson time step.
          atmp = - nu1 * (1.0/2.0) * gr%dzm(k) * gr%dzt(k)
          ! Completely implicit diffusional term.
!          atmp = - nu1 * gr%dzm(k) * gr%dzt(k)
          ! Diffusional term with Crank-Nicholson time step.
          ctmp = - nu1 * (1.0/2.0) * gr%dzm(k) * gr%dzt(kp1)
          ! Completely implicit diffusional term.
!          ctmp = - nu1 * gr%dzm(k) * gr%dzt(kp1)


          ! Momentum subdiagonal [ x wp2(k-1,<t+1>) ]
          aa(k_wp2,3-2)
     .    = - wmm(k) * (1.0/2.0) * gr%dzm(k)
     .      + atmp

          ! Thermodynamic subdiagonal [ x wp3(k,<t+1>) ]
          aa(k_wp2,3-1)
     .    = - gr%dzm(k)

          ! Momentum main diagonal [ x wp2(k,<t+1>) ]
          aa(k_wp2,3)
     .    = + 1.0 / dt
     .      + ( 1.0 - C5 ) * 2.0 * gr%dzm(k) * ( wmt(kp1) - wmt(k) )
!     .      + C1 / tau1m(k)
     .      + C1_Skw_fnc(k) / tau1m(k)
     .      - atmp - ctmp

          ! Thermodynamic superdiagonal [ x wp3(k+1,<t+1>) ]
          aa(k_wp2,3+1)
     .    = + gr%dzm(k)

          ! Momentum superdiagonal [ x wp2(k+1,<t+1>) ]
          aa(k_wp2,3+2)
     .    = + wmm(k) * (1.0/2.0) * gr%dzm(k)
     .      + ctmp


          ! Right-hand side
          dd(k_wp2)
     .    = + ( 1.0 / dt ) * wp2(k) 
     .      + ( 1.0 - C5 ) * 2.0 * ( grav / T0 ) * wpthvp(k)
     .      + (2.0/3.0) * C5
     .                  * ( ( grav / T0 ) * wpthvp(k)
     .                     - upwp(k) * gr%dzm(k) * ( um(kp1) - um(k) )
     .                     - vpwp(k) * gr%dzm(k) * ( vm(kp1) - vm(k) )
     .                    )
          ! These lines are for the diffusional term with a Crank-Nicholson 
          ! time step.  They need to be commented out for completely implicit
          ! diffusion to be used.
     .      - atmp * wp2(km1)
     .      + ( atmp + ctmp ) * wp2(k)
     .      - ctmp * wp2(kp1)


#ifdef STATS
          if (lstats_samp) then

!          Explicit contributions for wp2

           if ( iwp2_dp2 > 0 ) then
            zm%x(k,iwp2_dp2) = zm%x(k,iwp2_dp2) 
     .      - atmp * wp2(km1)
     .      + ( atmp + ctmp ) * wp2(k)
     .      - ctmp * wp2(kp1)
           endif

           if ( iwp2_bp > 0 ) then
            zm%x(k,iwp2_bp) = zm%x(k,iwp2_bp) 
     .      + 1.0 * 2.0 * grav * (1.0/T0) * wpthvp(k)
            zm%n(k,iwp2_bp) = zm%n(k,iwp2_bp) + 1
           endif

           if ( iwp2_pr2 > 0 ) then
            zm%x(k,iwp2_pr2) = zm%x(k,iwp2_pr2) 
     .      - C5 * 2.0 * grav * (1.0/T0) * wpthvp(k)
           endif

           if ( iwp2_pr3 > 0 ) then
            zm%x(k,iwp2_pr3) = zm%x(k,iwp2_pr3) 
     .      + (2.0/3.0) * C5
     .                  * (  grav * (1.0/T0) * wpthvp(k)
     .                     - upwp(k) * gr%dzm(k) * ( um(kp1) - um(k) )
     .                     - vpwp(k) * gr%dzm(k) * ( vm(kp1) - vm(k) )
     .                    )
            zm%n(k,iwp2_pr3) = zm%n(k,iwp2_pr3) + 1
           endif

!          Implicit contributions for wp2

           if ( iwp2_dp1 > 0 ) then
!             zmscr01(k) = - C1/tau1m(k)
             zmscr01(k) = - C1_Skw_fnc(k) / tau1m(k)
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

           if ( iwp2_ma > 0 ) then
             zmscr07(k) = + wmm(k) * (1.0/2.0) * gr%dzm(k)
             zmscr08(k) = - wmm(k) * (1.0/2.0) * gr%dzm(k)
           endif

           if ( iwp2_ac > 0 ) then
             zmscr09(k) = 
     .       - 1.0 * 2.0 * gr%dzm(k) * ( wmt(kp1) - wmt(k) )
           endif

           if ( iwp2_pr2 > 0 ) then
             zmscr10(k) = 
     .      + C5 * 2.0 * gr%dzm(k) * ( wmt(kp1) - wmt(k) )
           endif

         endif
#endif /*STATS*/


!         wp3
!  Added a new constant C12, a factor of atmp and ctmp 
!  Initially, this new constant will be set to 1.0 -dschanen 9/19/05 
          ! Diffusional term with Crank-Nicholson time step.
          atmp = - C12 * ( Kw(k) + nu8 ) 
     .                 * (1.0/2.0) * gr%dzt(k) * gr%dzm(km1)
          ! Completely implicit diffusional term.
!          atmp = - C12 * ( Kw(k) + nu8 ) * gr%dzt(k) * gr%dzm(km1)
          ! Diffusional term with Crank-Nicholson time step.
          ctmp = - C12 * ( Kw(k) + nu8 ) 
     .                 * (1.0/2.0) * gr%dzt(k) * gr%dzm(k)
          ! Completely implicit diffusional term.
!          ctmp = - C12 * ( Kw(k) + nu8 ) * gr%dzt(k) * gr%dzm(k)


          ! Thermodynamic subdiagonal [ x wp3(k-1,<t+1>) ]
          aa(k_wp3,3-2) = 
     .    - wmt(k) * (1.0/2.0) * gr%dzt(k)
     .    - (1.0/2.0) * gr%dzt(k)
     .      * a1(km1) * ( wp3(k) + wp3(km1) ) / MAX( wp2(km1), eps )
     .    + atmp

          ! Momentum subdiagonal [ x wp2(k-1,<t+1>) ]
          aa(k_wp3,3-1) = 
!     .    - 2.0 * gr%dzt(k) * a3(km1) * wp2(km1)   ! standard deriv.
!     .    - 2.0 * gr%dzt(k) * a3(k) * wp2(km1)     ! Golaz's no a3 deriv.
     .    - 2.0 * gr%dzt(k) * 0.5*(a3(k)+a3(km1)) * wp2(km1)   ! Centered no a3 deriv.

          ! Thermodynamic main diagonal [ ( x wp3(k,<t+1>) ]
          aa(k_wp3,3) = 
     .    + 1.0 / dt
     .    + (1.0/2.0) * gr%dzt(k)
     .      * (  a1(k) * ( wp3(kp1) + wp3(k) ) / MAX( wp2(k), eps )
     .         - a1(km1) * ( wp3(k) + wp3(km1) ) / MAX( wp2(km1), eps )
     .        )
     .    + ( 1.0 - C11_Skw_fnc(k) )
     .              * 3.0 * gr%dzt(k) * ( wmm(k) - wmm(km1) )
     .    + ( C8 / tauw3t(k) ) * ( 5.0 * C8b * Skwt(k)**4 + 1.0 )
     .    - atmp - ctmp

          ! Momentum superdiagonal [ x wp2(k,<t+1>) ]
          aa(k_wp3,3+1) = 
!     .    + 2.0 * gr%dzt(k) * a3(k) * wp2(k)   ! standard deriv.
!     .    + 2.0 * gr%dzt(k) * a3(k) * wp2(k)   ! Golaz's no a3 deriv.
     .    + 2.0 * gr%dzt(k) * 0.5*(a3(k)+a3(km1)) * wp2(k)   ! Centered no a3 deriv.

          ! Thermodynamic superdiagonal [ x wp3(k+1,<t+1>) ]
          aa(k_wp3,3+2) = 
     .    + wmt(k) * (1.0/2.0) * gr%dzt(k)
     .    + (1.0/2.0) * gr%dzt(k)
     .      * a1(k) * ( wp3(kp1) + wp3(k) ) / MAX( wp2(k), eps )
     .    + ctmp


          ! Right-hand side
          dd(k_wp3) = 
     .    + ( 1.0 / dt ) * wp3(k)
     .    + gr%dzt(k)
!     .      * ( ( a3(k) * wp2(k)**2 ) - ( a3(km1) * wp2(km1)**2 ) ) ! standard deriv.
!     .      * ( ( a3(k) * wp2(k)**2 ) - ( a3(k) * wp2(km1)**2 ) )   ! Golaz's no a3 deriv.
     .      * (  ( 0.5*(a3(k)+a3(km1)) * wp2(k)**2 )        ! Centered no a3 deriv.
     .         - ( 0.5*(a3(k)+a3(km1)) * wp2(km1)**2 )  )   ! Centered no a3 deriv.
     .    + (1.0/4.0) * gr%dzt(k)
     .      * (  a1(k) * ( wp3(kp1) + wp3(k) )**2 
     .                   / MAX( wp2(k), eps )
     .         - a1(km1) * ( wp3(k) + wp3(km1) )**2
     .                     / MAX( wp2(km1), eps )
     .        )
     .    + ( 1.0 - C11_Skw_fnc(k) ) * 3.0 * ( grav/T0 ) * wp2thvp(k)
     .    + ( C8 / tauw3t(k) ) * ( 4.0 * C8b * Skwt(k)**4 ) * wp3(k)
          ! These lines are for the diffusional term with a Crank-Nicholson 
          ! time step.  They need to be commented out for completely implicit
          ! diffusion to be used.
     .    - atmp * wp3(km1)
     .    + ( atmp + ctmp ) * wp3(k)
     .    - ctmp * wp3(kp1)


#ifdef STATS
          if (lstats_samp) then

!         Explicit contributions for wp3

          if ( iwp3_ta > 0 ) then
          zt%x(k,iwp3_ta) = zt%x(k,iwp3_ta) 
     .    + gr%dzt(k)
!     .      * (  ( ( a3(k) + (3.0/2.0) ) * wp2(k)**2 )      ! standard deriv.
!     .         - ( ( a3(km1) + (3.0/2.0) ) * wp2(km1)**2 )  ! standard deriv.
!     .        )                                               ! standard deriv.
!     .      * (  ( ( a3(k) + (3.0/2.0) ) * wp2(k)**2 )      ! Golaz's no a3 deriv.
!     .         - ( ( a3(k) + (3.0/2.0) ) * wp2(km1)**2 )    ! Golaz's no a3 deriv.
!     .        )                                               ! Golaz's no a3 deriv.
     .      * (  ( ( 0.5*(a3(k)+a3(km1)) + (3.0/2.0) )        ! Centered no a3 deriv.
     .             * wp2(k)**2 )                            ! Centered no a3 deriv.
     .         - ( ( 0.5*(a3(k)+a3(km1)) + (3.0/2.0) )        ! Centered no a3 deriv.
     .             * wp2(km1)**2 )                          ! Centered no a3 deriv.
     .        )                                               ! Centered no a3 deriv.
     .    + (1.0/4.0) * gr%dzt(k)
     .      * (  a1(k) * ( wp3(kp1) + wp3(k) )**2 
     .                   / MAX( wp2(k), eps )
     .         - a1(km1) * ( wp3(k) + wp3(km1) )**2 
     .                     / MAX( wp2(km1), eps )
     .        )
          endif

          if ( iwp3_tp > 0 ) then
          zt%x(k,iwp3_tp) = zt%x(k,iwp3_tp) 
     .    - gr%dzt(k)
     .      * (  ( (3.0/2.0) * wp2(k)**2 ) 
     .         - ( (3.0/2.0) * wp2(km1)**2 )
     .        )
          endif

          if ( iwp3_bp > 0 ) then
          zt%x(k,iwp3_bp) = zt%x(k,iwp3_bp) 
     .    + 1.0 * 3.0 * ( grav / T0 )  * wp2thvp(k)
          zt%n(k,iwp3_bp) = zt%n(k,iwp3_bp) + 1
          endif

          if ( iwp3_pr2 > 0 ) then
          zt%x(k,iwp3_pr2) = zt%x(k,iwp3_pr2) 
     .    - C11_Skw_fnc(k) * 3.0 * grav * ( grav / T0 ) * wp2thvp(k)
          endif

          if ( iwp3_pr1 > 0 ) then
          zt%x(k,iwp3_pr1) = zt%x(k,iwp3_pr1) 
     .    + ( C8 / tauw3t(k) ) * ( 4.0 * C8b * Skwt(k)**4 ) * wp3(k)
          endif

          ! These lines are for the diffusional term with a Crank-Nicholson 
          ! time step.  They need to be commented out for completely implicit
          ! diffusion to be used.
          if ( iwp3_dp1 > 0 ) then
          zt%x(k,iwp3_dp1) = zt%x(k,iwp3_dp1) 
     .    - atmp * wp3(km1)
     .    + ( atmp + ctmp ) * wp3(k)
     .    - ctmp * wp3(kp1)
          endif

!         Implicit contributions for wp3

          if ( iwp3_ta > 0 ) then
          ztscr05(k) =
     .    + (1.0/2.0) * gr%dzt(k)
     .      * a1(km1) * ( wp3(k) + wp3(km1) ) / MAX( wp2(km1), eps )
          ztscr06(k) =
!     .    + 2.0 * gr%dzt(k) * ( a3(km1) + (3.0/2.0) ) * wp2(km1)  ! standard deriv.
!     .    + 2.0 * gr%dzt(k) * ( a3(k) + (3.0/2.0) ) * wp2(km1)    ! Golaz's no a3 deriv.
     .    + 2.0 * gr%dzt(k) * ( 0.5*(a3(k)+a3(km1)) + (3.0/2.0) ) ! Centered no a3 deriv.
     .      * wp2(km1)                                            ! Centered no a3 deriv.
          ztscr07(k) =
     .    - (1.0/2.0) * gr%dzt(k)
     .      * (  a1(k) * ( wp3(kp1) + wp3(k) ) / MAX( wp2(k), eps )
     .         - a1(km1) * ( wp3(k) + wp3(km1) ) / MAX( wp2(km1), eps )
     .        )
          ztscr08(k) =
!     .    - 2.0 * gr%dzt(k) * ( a3(k) + (3.0/2.0) ) * wp2(k)      ! standard deriv.
!     .    - 2.0 * gr%dzt(k) * ( a3(k) + (3.0/2.0) ) * wp2(k)      ! Golaz's no a3 deriv.
     .    - 2.0 * gr%dzt(k) * ( 0.5*(a3(k)+a3(km1)) + (3.0/2.0) ) ! Centered a3 deriv.
     .      * wp2(k)                                              ! Centered no a3 deriv.
          ztscr09(k) =
     .    - (1.0/2.0) * gr%dzt(k)
     .      * a1(k) * ( wp3(kp1) + wp3(k) ) / MAX( wp2(k), eps )
          endif

          if ( iwp3_tp > 0 ) then
          ztscr10(k) = 
     .    - 2.0 * gr%dzt(k) * (3.0/2.0) * wp2(km1)
          ztscr11(k) = 
     .    + 2.0 * gr%dzt(k) * (3.0/2.0) * wp2(k)
          endif

          if ( iwp3_ma > 0 ) then
          ztscr12(k) = 
     .    + wmt(k) * (1.0/2.0) * gr%dzt(k)
          ztscr13(k) = 
     .    - wmt(k) * (1.0/2.0) * gr%dzt(k)
          endif

          if ( iwp3_ac > 0 ) then
          ztscr14(k) = 
     .    - 1.0 * 3.0 * gr%dzt(k) * ( wmm(k) - wmm(km1) )
          endif

          if ( iwp3_pr2 > 0 ) then
          ztscr15(k) =
     .    + C11_Skw_fnc(k) * 3.0 * gr%dzt(k) * ( wmm(k) - wmm(km1) )
          endif

          if ( iwp3_pr1 > 0 ) then
          ztscr01(k) =
     .    - ( C8 / tauw3t(k) ) * ( 5.0 * C8b * Skwt(k)**4 + 1.0 )
          endif

          if ( iwp3_dp1 > 0 ) then
          ztscr02(k) = - atmp
          ztscr03(k) = atmp + ctmp
          ztscr04(k) = - ctmp
          endif

          endif
#endif /*STATS*/

        ENDDO


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
          print *, "wp23 LU decomp. failed"
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

            if ( iwp2_ma > 0 ) then
             zm%x(k,iwp2_ma) = zm%x(k,iwp2_ma)
     .       + zmscr07(k) * wp2(km1)
     .       + zmscr08(k) * wp2(kp1)
             zm%n(k,iwp2_ma) = zm%n(k,iwp2_ma) + 1
            endif

            if ( iwp2_ac > 0 ) then
             zm%x(k,iwp2_ac) = zm%x(k,iwp2_ac)
     .       + zmscr09(k) * wp2(k)
             zm%n(k,iwp2_ac) = zm%n(k,iwp2_ac) + 1
            endif

            if ( iwp2_pr2 > 0 ) then
             zm%x(k,iwp2_pr2) = zm%x(k,iwp2_pr2)
     .       + zmscr10(k) * wp2(k)
             zm%n(k,iwp2_pr2) = zm%n(k,iwp2_pr2) + 1
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

            if ( iwp3_ma > 0 ) then
             zt%x(k,iwp3_ma) = zt%x(k,iwp3_ma)
     .       + ztscr12(k) * wp3(km1)
     .       + ztscr13(k) * wp3(kp1)
             zt%n(k,iwp3_ma) = zt%n(k,iwp3_ma) + 1
            endif

            if ( iwp3_ac > 0 ) then
             zt%x(k,iwp3_ac) = zt%x(k,iwp3_ac)
     .       + ztscr14(k) * wp3(k)
             zt%n(k,iwp3_ac) = zt%n(k,iwp3_ac) + 1
            endif

            if ( iwp3_pr2 > 0 ) then
             zt%x(k,iwp3_pr2) = zt%x(k,iwp3_pr2) 
     .       + ztscr15(k) * wp3(k)
             zt%n(k,iwp3_pr2) = zt%n(k,iwp3_pr2) + 1
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
#endif /*STATS*/

        ! Clip Skewness.
 
        DO k = 1, gr%nnzp, 1

           km1 = MAX(k-1,1)

           ! Clipping wp3 at the first layer above ground according to 
           ! Andre et al. (1976b & 1978).
           ! According to Andre et al. (1976b & 1978), wp3 should not
           ! exceed [2*(wp2^3)]^(1/2) at any level.  However, this term
           ! should be multiplied by 0.2 close to the surface to include
           ! surface effects.  There already is a wp3 clipping term in
           ! place for all other altitudes, but this term will be 
           ! included for the surface layer only.
           ! Therefore, the lowest level wp3 should not exceed
           ! 0.2 * SQRT(2) * wp2^(3/2).  Brian Griffin.  12/18/05.

           ! NOTE: Clipping not good; must try to get rid of.  Effects only
           !       Wangara case.

           IF ( gr%zt(k) <= 100.0 .AND. bottom_at_sfc ) THEN ! Clip for 100 m. above ground.

              atmp = MAX( ( (wp2(km1)+wp2(k)) / 2.0 ) , eps )    ! wp2
              ctmp = atmp**(3.0/2.0)                             ! wp2^(3/2)
              IF ( wp3(k) >= 0.2 * SQRT(2.0) * ctmp ) THEN
                 wp3(k) = 0.2 * SQRT(2.0) * ctmp
              ENDIF

           ELSE

              ! Clip skewness consistently with a
              atmp = MAX(0.5*(wp2(km1)+wp2(k)),eps)
              ctmp = atmp**1.5
              IF ( wp3(k)/ctmp .gt. 4.5 ) THEN
                 wp3(k) = 4.5 * ctmp
                 wp2(k) = atmp
              ELSEIF ( wp3(k)/ctmp .lt. -4.5 ) THEN
                 wp3(k) = -4.5 * ctmp
                 wp2(k) = atmp
              ENDIF

           ENDIF

           wp2(k) = max( 2./3. * emin, wp2(k) ) 

        ENDDO

        ! wp3 below the surface is set to 0.
        wp3(1) = 0.0

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
#endif /*STATS*/

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
#endif /*STATS*/

        return
        end subroutine wp23_solve

