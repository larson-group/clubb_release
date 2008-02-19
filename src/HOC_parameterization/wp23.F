!------------------------------------------------------------------------
! $Id: wp23.F,v 1.37 2008-02-19 18:29:02 vlarson Exp $

        subroutine timestep_wp23
     .             ( dt, Lscale, taum, taut, Scm, Skwm, Skwt,
     .               a, Kht, Khm, wmt, wmm, up2, vp2,
     .               um, upwp, vm, vpwp, wpthvp, wp2thvp,
     .               wp2, wp3, isValid )
!       Description:
!       Advance wp2 and wp3 one timestep

!       References:
!       Eqn. 12 & 18 on p. 3545--3546 of
!       ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!         Method and Model Description'' Golaz, et al. (2002)
!       JAS, Vol. 59, pp. 3540--3551.

!       See also 
!       ``Equations for HOC'', Section 6:
!       /Implict solution for the vertical velocity moments/
!------------------------------------------------------------------------

        use grid_class
        use parameters
        use constants
        use model_flags
#ifdef STATS
        use hoc_stats
#endif /*STATS*/

        implicit none

        ! External
        external :: wp23_solve

        intrinsic :: exp

        ! Input Variables
        real, intent(in) :: dt ! Timestep               [s]

        real, intent(in), dimension(gr%nnzp) :: 
     .  Lscale,! Mixing Length                          [m]
     .  taum,  ! Tau on moment. grid                    [s]
     .  taut,  ! Tau on thermo. grid                    [s]
     .  Scm,   ! Sc on moment. grid                     [-]
     .  Skwm,  ! Skw on moment. grid                    [-]
     .  Skwt,  ! Skw on thermo. grid                    [-]
     .  a,     ! PDF parameter a (13)                   [-]
     .  Kht,   ! Kh Eddy diffusivity on thermo. grid    [m^2/s]
     .  Khm,   ! Kh Eddy diffusivity on momentum grid   [m^2/s]
     .  wmt,   ! w wind on the thermo. grid             [m/s]
     .  wmm,   ! w wind on the momentum grid            [m/s]
     .  up2,   ! u'^2                                   [m^2/s^2]
     .  vp2,   ! v'^2                                   [m^2/s^2]
     .  um,    ! u wind (thermo. grid)                  [m/s]
     .  upwp,  ! u'w' (moment. grid)                    [m^2/s^2]
     .  vm,    ! v wind (thermo. grid)                  [m/s]
     .  vpwp,  ! u'w' (moment. grid)                    [m^2/s^2]
     .  wpthvp,! w'th_v' Buoyancy flux                  [(K m)/s]
     .  wp2thvp! w'^2 th_v'                             [(m^2 K)/s^2]

        ! Input/Output
        real, dimension(gr%nnzp), intent(inout) :: 
     .  wp2, ! w'^2 (moment. grid)                      [m^2/s^2]
     .  wp3  ! w'^3 (thermo. grid)                      [m^3/s^3]

        logical, intent(inout) :: isValid ! Diagnostic

        ! Local Variables
        real, dimension(gr%nnzp) :: 
     .  tauw3t, ! Currently just taut           [s]
     .  Lm      ! Length scale (moment. grid)   [m]

        ! Eddy Diffusion for wp2 and wp3.
        real, dimension(gr%nnzp) :: Kw1   ! wp2 eddy diff. [m^2/s]
        real, dimension(gr%nnzp) :: Kw8   ! wp3 eddy diff. [m^2/s]

        ! Note:  wp2_zt and wp2_zt_sqd_3pt, and wp3_zm and wp3_zm_sqd_3pt 
        !        are used to help determine the coefficients of eddy 
        !        diffusivity for wp2 and wp3, respectively.
        real, dimension(gr%nnzp) ::
     .  wp2_zt,         ! wp2 interpolated to thermo. levels [m^2/s^2]
     .  wp3_zm,         ! wp3 interpolated to moment. levels [m^3/s^3]
     .  wp2_zt_sqd_3pt, ! wp2**2, avgd over 3 points         [m^4/s^4]
     .  wp3_zm_sqd_3pt  ! wp3**2, avgd over 3 points         [m^6/s^6]

        ! Internal variables for C11 function, Vince Larson 13 Mar 2005
        ! Brian added C1 function.
        real, dimension(gr%nnzp) :: 
     .  C1_Skw_fnc, !   [-]
     .  C11_Skw_fnc !   [-]
        ! End Vince Larson's addition.
      
        integer :: k, km1, kp1  ! Array indices

!-----------------------------------------------------------------------

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
!        end do

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

!       END DO

        C11_Skw_fnc(1:gr%nnzp) = 
     .  C11b + (C11-C11b)*EXP( -(1.0/2.0) * (Skwt(1:gr%nnzp)/C11c)**2 )

        C1_Skw_fnc(1:gr%nnzp) = 
     .  C1b + (C1-C1b)*EXP( -(1.0/2.0) * (Skwm(1:gr%nnzp)/C1c)**2 )

!        C11_Skw_fnc = C11
!        C1_Skw_fnc = C1

        ! Vince Larson added extra diffusion based on wp2.  21 Dec 2007.
        ! Vince Larson added extra diffusion based on wp3.  15 Dec 2007.

        wp2_zt = zm2zt(wp2)
        wp3_zm = zt2zm(wp3)

        do k = 1, gr%nnzp, 1

           km1 = max( k-1, 1 )
           kp1 = min( k+1, gr%nnzp )

           ! Compute the square of wp2_zt, averaged over 3 points.  15 Dec 2007
           wp2_zt_sqd_3pt(k) 
     .          = ( wp2_zt(km1)**2 + wp2_zt(k)**2 + wp2_zt(kp1)**2 ) / 3.0

           ! Compute the square of wp3_zm, averaged over 3 points.  15 Dec 2007
           wp3_zm_sqd_3pt(k) 
     .          = ( wp3_zm(km1)**2 + wp3_zm(k)**2 + wp3_zm(kp1)**2 ) / 3.0

        enddo

        ! End Vince Larson's addition.

        ! Define the Coefficent of Eddy Diffusivity for the wp2 and wp3.
        do k = 1, gr%nnzp, 1
           ! Kw1 is used for wp2, which is located on momentum levels.
           ! Kw1 is located on thermodynamic levels.
           ! Kw1 = c_K1 * Kht
           Kw1(k) = c_K1 * Kht(k)
           ! Vince Larson added extra diffusion based on wp2.  21 Dec 2007.
           ! Kw1 must have units of m^2/s.  Since wp2_zt_sqd_3pt has units 
           ! of m^4/s^4, c_Ksqd is given units of s^3/m^2 in this case.
           Kw1(k) = Kw1(k) + c_Ksqd * wp2_zt_sqd_3pt(k) ! Vince Larson increased by c_Ksqd. 29Jan2008
           ! End Vince Larson's addition.
           ! Kw8 is used for wp3, which is located on thermodynamic levels.
           ! Kw8 is located on momentum levels.
           ! Note: Kw8 is defined to be 1/2 of Khm.
           ! Kw8 = c_K8 * Khm
           Kw8(k) = c_K8 * Khm(k)
           ! Vince Larson added extra diffusion based on wp3.  15 Dec 2007.
           ! Kw8 must have units of m^2/s.  Since wp3_zm_sqd_3pt has units 
           ! of m^6/s^6, c_Ksqd is given units of s^5/m^4 in this case.
           Kw8(k) = Kw8(k) + c_Ksqd * wp3_zm_sqd_3pt(k)  ! Vince Larson increased by c_Ksqd. 29Jan2008
           ! End Vince Larson's addition.
        enddo

!       Solve semi-implicitly
! Vince Larson increased Kw8 in order to run arm_0003 and arm_97 stably
!  14 Nov 2007.
        call wp23_solve( dt, Scm, Lscale, Lm, taum, tauw3t, Kw1, 
     .                   Kw8, wmt, wmm, wpthvp, wp2thvp,
     .                   um, upwp, vm, vpwp, up2, vp2,
     .                   Skwt, C1_Skw_fnc, C11_Skw_fnc, 
     .                   wp2, wp3, isValid )
!        call wp23_solve( dt, Scm, Lscale, Lm, taum, tauw3t, 10*0.5*Kht,
!     .                   wmt, wmm, wpthvp, wp2thvp,
!     .                   um, upwp, vm, vpwp, up2, vp2,
!     .                   Skwt, C1_Skw_fnc, C11_Skw_fnc, 
!     .                   wp2, wp3, isValid )
! End Vince Larson's change

!       Error output
!       Joshua Fasching Feb 2008
        if ( .not. isValid ) then
                
           write(fstderr,*) "Errors in timestep_wp32"

           write(fstderr,*) "Intent(in)"
           
           write(fstderr,*) "Lscale = ", Lscale
           write(fstderr,*) "taum = ", taum
           write(fstderr,*) "taut = ", taut
           write(fstderr,*) "Scm = ", Scm
           write(fstderr,*) "Skwm = ", Skwm
           write(fstderr,*) "Skwt = ", Skwt
           write(fstderr,*) "a = ", a
           write(fstderr,*) "Kht = ", Kht
           write(fstderr,*) "Khm = ", Khm
           write(fstderr,*) "wmt = ", wmt
           write(fstderr,*) "wmm = ", wmm
           write(fstderr,*) "up2 = ", up2
           write(fstderr,*) "vp2 = ", vp2
           write(fstderr,*) "um = ", um
           write(fstderr,*) "upwp = ", upwp
           write(fstderr,*) "vm = ", vm
           write(fstderr,*) "vpwp = ", vpwp
           write(fstderr,*) "wpthvp = ", wpthvp
           write(fstderr,*) "wp2thvp = ", wp2thvp
           
           write(fstderr,*) "Intent(in/out)"
           
           write(fstderr,*) "wp2 = ", wp2
           write(fstderr,*) "wp3 = ", wp3
           
        endif
        
        return
        
        end subroutine timestep_wp23

!------------------------------------------------------------------------
        subroutine wp23_solve( dt, Scm, Lscale, Lm, tau1m, tauw3t, Kw1,
     .                         Kw8, wmt, wmm, wpthvp, wp2thvp,
     .                         um, upwp, vm, vpwp, up2, vp2,
     .                         Skwt, C1_Skw_fnc, C11_Skw_fnc,
     .                         wp2, wp3, isValid)
!       Description:
!       Decompose, and back substitute the matrix for wp2/wp3

!       References:
!       _Equations for HOC_ section 6.3 
!------------------------------------------------------------------------

        use grid_class
        use parameters
        use constants
        use model_flags
        use diffusion
        use lapack_wrap, only: band_solve, band_solvex
#ifdef STATS
        use hoc_stats
#endif /*STATS*/

        implicit none

        ! External
        intrinsic :: max, min, sqrt

        ! Parameter Constants
        integer, parameter ::
     .  nsub = 2,  ! Number of subdiagonals in the LHS matrix
     .  nsup = 2,  ! Number of superdiagonals in the LHS matrix
     .  nrhs = 1   ! Number of RHS vectors

        ! Input Variables
        real, intent(in) :: 
     .  dt          ! Timestep                               [s]

        real, intent(in), dimension(gr%nnzp) :: 
     .  Scm,        ! Sc on moment. grid                     [-]
     .  Lscale,     ! Mixing Length                          [m]
     .  Lm,         ! Mixing Length interp. to moment.       [m]
     .  tau1m,      ! Tau on moment. grid                    [s]
     .  tauw3t,     ! Tau on thermo. grid                    [s]
     .  Kw1,        ! Eddy diffusivity                       [m^2/s]
     .  Kw8,        ! 1/2 Khm Eddy diffusivity               [m^2/s]
     .  wmt,        ! w wind on t                            [m/s]
     .  wmm,        ! w wind on m                            [m/s]
     .  up2,        ! u'^2                                   [m^s/s^2]
     .  vp2,        ! v'^2                                   [m^s/s^2]
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
        real, dimension(nsup+nsub+1,2*gr%nnzp) :: 
     .  lhs ! Implicit contributions to wp2/wp3 (band diag. matrix)

        real, dimension(2*gr%nnzp) :: 
     .  rhs,  ! RHS of band matrix
     .  solut ! Solution to band diagonal system.

        real, dimension(gr%nnzp) :: 
     .  a1, ! a1; See eqn. 20 in `Equations for HOC'             [-]
     .  a3  ! a3; See eqn. 22 in `Equations for HOC'             [-]

        real :: 
     .  atmp, ctmp, ! Temporary variables for the LHS
     .  rcond  ! Est. of the reciprocal of the condition #

        ! Array indices
        integer :: k, km1, kp1, k_wp2, k_wp3

        ! Added by Brian for eddy diffusion.
        real, dimension(3) :: diff_lhs
        ! Set logical to true for Crank-Nicholson diffusion scheme
        ! or to false for completely implicit diffusion scheme.
        ! Note:  Although Crank-Nicholson diffusion has usually been
        !        used for wp2 and wp3 in the past, we found that using
        !        completely implicit diffusion stabilized the deep
        !        convective cases more while having almost no effect on
        !        the boundary layer cases.  Brian; 1/4/2008.
!        logical, parameter :: lcrank_nich_diff = .true.
        logical, parameter :: lcrank_nich_diff = .false.
        ! Constant parameters
        integer, parameter ::
     .  kp1diag = 1,
     .  kdiag   = 2,
     .  km1diag = 3

#ifdef STATS
        if (lstats_samp) then

          if ( iwp2_bt > 0 ) then
            zm%x(:,iwp2_bt) = zm%x(:,iwp2_bt) - wp2 / dt
          end if

          if ( iwp3_bt > 0 ) then
            zt%x(:,iwp3_bt) = zt%x(:,iwp3_bt) - wp3 / dt
          end if

        end if
#endif /*STATS*/

        ! Define a1 and a3 (both are located on momentum levels).
        ! Note: some compilers appear to interpret the pow function with
        ! a positive integer exponent differently than a repeated
        ! multiply. -dschanen 19 March 2007

        a1 = 1.0 / ( 1.0 - Scm )
        a3 = 3.0 * Scm*Scm
     .       + 6.0*(1.0-Scm)*Scm 
     .       + (1.0-Scm)*(1.0-Scm)
     .       - (3.0/2.0)

        ! Set up matrix lhs and vector rhs

        lhs = 0.0
        rhs = 0.0

!        do k = 1, gr%nnzp, 1
        do k = 2, gr%nnzp-1, 1

!         Define indices

          km1 = max( k-1, 1 )
          kp1 = min( k+1, gr%nnzp )

          k_wp3 = 2*k - 1
          k_wp2 = 2*k

          ! Skewness of w is needed on the thermodynamic levels 
          ! for the wp3 equation
          !Skw = wp3(k) / ( (1.0/2.0)*( wp2(k) + wp2(km1) ) )**(3.0/2.0)
          ! Skewness of w limited to a magnitude of no more than 4.5
          !Skw = MIN( MAX( Skw, -4.5 ), 4.5 )

          ! wp2

          ! Diffusional term with Crank-Nicholson time step.
          ! atmp = - nu1 * (1.0/2.0) * gr%dzm(k) * gr%dzt(k)
          ! Completely implicit diffusional term.
          ! atmp = - nu1 * gr%dzm(k) * gr%dzt(k)
          ! Diffusional term with Crank-Nicholson time step.
          ! ctmp = - nu1 * (1.0/2.0) * gr%dzm(k) * gr%dzt(kp1)
          ! Completely implicit diffusional term.
          ! ctmp = - nu1 * gr%dzm(k) * gr%dzt(kp1)

          if ( lcrank_nich_diff ) then
             ! Eddy diffusion for wp2 using a Crank-Nicholson time step. 
             diff_lhs(kp1diag:km1diag)
     .       = (1.0/2.0)
     .       * diffusion_mlev( Kw1(k), Kw1(kp1), nu1,
     .                         gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )
          else
             ! Eddy diffusion for wp2 using a completely implicit time step. 
             diff_lhs(kp1diag:km1diag)
     .       = diffusion_mlev( Kw1(k), Kw1(kp1), nu1,
     .                         gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )
          endif

          ! Momentum subdiagonal [ x wp2(k-1,<t+1>) ]
          lhs(3+2,k_wp2)
     .    = - wmm(k) * (1.0/2.0) * gr%dzm(k)
     .      + diff_lhs(3)

          ! Thermodynamic subdiagonal [ x wp3(k,<t+1>) ]
          lhs(3+1,k_wp2)
     .    = - gr%dzm(k)

          ! Momentum main diagonal [ x wp2(k,<t+1>) ]
          lhs(3,k_wp2)
     .    = + 1.0 / dt
     .      + ( 1.0 - C5 ) * 2.0 * gr%dzm(k) * ( wmt(kp1) - wmt(k) )
     .      + C1_Skw_fnc(k) / tau1m(k)
     .      + diff_lhs(2)
          if ( lKhm_aniso ) then
            ! Add in this term if we're not assuming tke = 1.5 * wp2
            lhs(3,k_wp2) = lhs(3,k_wp2)
     .      + ( 2.0 * C4 ) / ( 3.0 * tau1m(k) )
          end if

          ! Thermodynamic superdiagonal [ x wp3(k+1,<t+1>) ]
          lhs(3-1,k_wp2)
     .    = + gr%dzm(k)

          ! Momentum superdiagonal [ x wp2(k+1,<t+1>) ]
          lhs(3-2,k_wp2)
     .    = + wmm(k) * (1.0/2.0) * gr%dzm(k)
     .      + diff_lhs(1)


          ! Right-hand side
          rhs(k_wp2)
     .    = + ( 1.0 / dt ) * wp2(k) 
     .      + ( 1.0 - C5 ) * 2.0 * ( grav / T0 ) * wpthvp(k)
! Michael Falk, 2 August 2007
! Use the following code for standard mixing, with c_k=0.548:
!     .      + (2.0/3.0) * C5
!     .                  * ( ( grav / T0 ) * wpthvp(k)
!     .                     - upwp(k) * gr%dzm(k) * ( um(kp1) - um(k) )
!     .                     - vpwp(k) * gr%dzm(k) * ( vm(kp1) - vm(k) )
!     .                    )
! Use the following code for alternate mixing, with c_k=0.1 or 0.2
     .      + (2.0/3.0) * C5
     .                  * ( ( grav / T0 ) * wpthvp(k)
     .                  - 0. * upwp(k) * gr%dzm(k) * ( um(kp1) - um(k) )
     .                  - 0. * vpwp(k) * gr%dzm(k) * ( vm(kp1) - vm(k) )
     .                    )
! eMFc

          if ( lcrank_nich_diff ) then
             ! These lines are for the diffusional term with a Crank-Nicholson 
             ! time step.  They are not used for completely implicit diffusion.
             rhs(k_wp2)   =   rhs(k_wp2)
     .                      - diff_lhs(3) * wp2(km1)
     .                      - diff_lhs(2) * wp2(k)
     .                      - diff_lhs(1) * wp2(kp1)
          endif

          if ( lKhm_aniso ) then
            rhs(k_wp2) = rhs(k_wp2)
     .      + ( C4 * ( up2(k) + vp2(k) ) ) / ( 3.0 * tau1m(k) )
          end if


#ifdef STATS
          if ( lstats_samp ) then

!          Explicit contributions for wp2

           if ( iwp2_dp2 > 0 ) then
            if ( lcrank_nich_diff ) then
             zm%x(k,iwp2_dp2) = zm%x(k,iwp2_dp2) 
     .       - diff_lhs(3) * wp2(km1)
     .       - diff_lhs(2) * wp2(k)
     .       - diff_lhs(1) * wp2(kp1)
            endif
           end if

           if ( iwp2_bp > 0 ) then
            zm%x(k,iwp2_bp) = zm%x(k,iwp2_bp) 
     .      + 1.0 * 2.0 * grav * (1.0/T0) * wpthvp(k)
            zm%n(k,iwp2_bp) = zm%n(k,iwp2_bp) + 1
           end if

           if ( iwp2_pr1 > 0 .and. lKhm_aniso ) then
            zm%x(k,iwp2_pr1) = zm%x(k,iwp2_pr1) 
     .      + ( C4 * ( up2(k) + vp2(k) ) ) / ( 3.0 * tau1m(k) )
           end if

           if ( iwp2_pr2 > 0 ) then
            zm%x(k,iwp2_pr2) = zm%x(k,iwp2_pr2) 
     .      - C5 * 2.0 * grav * (1.0/T0) * wpthvp(k)
           end if

           if ( iwp2_pr3 > 0 ) then
            zm%x(k,iwp2_pr3) = zm%x(k,iwp2_pr3) 
     .      + (2.0/3.0) * C5
     .                  * (  grav * (1.0/T0) * wpthvp(k)
! Vince Larson made budget reflect Michael Falk's zeroing of shear terms. 20 Nov 2007
!     .                     - upwp(k) * gr%dzm(k) * ( um(kp1) - um(k) )
!     .                     - vpwp(k) * gr%dzm(k) * ( vm(kp1) - vm(k) )
!     .                    )
     .                - 0. * upwp(k) * gr%dzm(k) * ( um(kp1) - um(k) )
     .                - 0. * vpwp(k) * gr%dzm(k) * ( vm(kp1) - vm(k) )
     .                    )
! End of Vince Larson' change.
            zm%n(k,iwp2_pr3) = zm%n(k,iwp2_pr3) + 1
           end if

           ! Implicit contributions for wp2

           if ( iwp2_dp1 > 0 ) then
             zmscr01(k) = - C1_Skw_fnc(k) / tau1m(k)
           end if

           if ( iwp2_dp2 > 0 ) then
             zmscr02(k) = - diff_lhs(3)
             zmscr03(k) = - diff_lhs(2)
             zmscr04(k) = - diff_lhs(1)
           end if

           if ( iwp2_ta > 0 ) then
             zmscr05(k) = + gr%dzm(k)
             zmscr06(k) = - gr%dzm(k)
           end if

           if ( iwp2_ma > 0 ) then
             zmscr07(k) = + wmm(k) * (1.0/2.0) * gr%dzm(k)
             zmscr08(k) = - wmm(k) * (1.0/2.0) * gr%dzm(k)
           end if

           if ( iwp2_ac > 0 ) then
             zmscr09(k) = 
     .       - 1.0 * 2.0 * gr%dzm(k) * ( wmt(kp1) - wmt(k) )
           end if

           if ( iwp2_pr2 > 0 ) then
             zmscr10(k) = 
     .      + C5 * 2.0 * gr%dzm(k) * ( wmt(kp1) - wmt(k) )
           end if

           if ( iwp2_pr1 > 0 .and. lKhm_aniso ) then
             zmscr11(k) = - ( 2.0 * C4 ) / ( 3.0 * tau1m(k) )
           end if

         end if
#endif /*STATS*/


!         wp3
!  Added a new constant C12, a factor of atmp and ctmp 
!  Initially, this new constant will be set to 1.0 -dschanen 9/19/05 
          ! Diffusional term with Crank-Nicholson time step.
!          atmp = - C12 * ( Kw(k) + nu8 ) 
!     .                 * (1.0/2.0) * gr%dzt(k) * gr%dzm(km1)
          ! Completely implicit diffusional term.
!          atmp = - C12 * ( Kw(k) + nu8 ) * gr%dzt(k) * gr%dzm(km1)
          ! Diffusional term with Crank-Nicholson time step.
!          ctmp = - C12 * ( Kw(k) + nu8 ) 
!     .                 * (1.0/2.0) * gr%dzt(k) * gr%dzm(k)
          ! Completely implicit diffusional term.
!          ctmp = - C12 * ( Kw(k) + nu8 ) * gr%dzt(k) * gr%dzm(k)

          if ( lcrank_nich_diff ) then
             ! Eddy diffusion for wp3 using a Crank-Nicholson time step. 
             diff_lhs(kp1diag:km1diag)
     .       = C12 * (1.0/2.0)
     .       * diffusion_tlev( Kw8(k), Kw8(km1), nu8,
     .                         gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
          else
             ! Eddy diffusion for wp3 using a completely implicit time step. 
             diff_lhs(kp1diag:km1diag)
     .       = C12 
     .       * diffusion_tlev( Kw8(k), Kw8(km1), nu8,
     .                         gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
          endif

          ! Thermodynamic subdiagonal [ x wp3(k-1,<t+1>) ]
          lhs(3+2,k_wp3) = 
     .    - wmt(k) * (1.0/2.0) * gr%dzt(k)
          ! Brian tried a new discretization for the turbulent advection term, 
          ! which contains the term -d[ a1 * wp3^2 / wp2 ] / dz.  In order to 
          ! help stabilize wp3, a1 has been pulled outside of the derivative.
!     .    - (1.0/2.0) * gr%dzt(k)
!     .      * a1(km1) * ( wp3(k) + wp3(km1) ) / MAX( wp2(km1), eps )
     .    - (1.0/2.0) * ( a1(k) + a1(km1) )
     .      * (1.0/2.0) * gr%dzt(k)
     .      * ( wp3(k) + wp3(km1) ) / MAX( wp2(km1), eps )
          ! End of Brian's a1 change.  Feb. 14, 2008.
     .    + diff_lhs(3)

          ! Momentum subdiagonal [ x wp2(k-1,<t+1>) ]
          lhs(3+1,k_wp3) = 
          ! The discretization of the turbulent advection term containing 
          ! the term -d[ (a3 + 3/2) * wp2^2 ] / dz has been altered to pull 
          ! a3 outside of the derivative.
!     .    - 2.0 * gr%dzt(k) * a3(km1) * wp2(km1)   ! standard deriv.
!     .    - 2.0 * gr%dzt(k) * a3(k) * wp2(km1)     ! Golaz's no a3 deriv.
     .    - 2.0 * gr%dzt(k) * 0.5*(a3(k)+a3(km1)) * wp2(km1)   ! Centered no a3 deriv.
          ! End of code that pulls out a3.

          ! Thermodynamic main diagonal [ ( x wp3(k,<t+1>) ]
          lhs(3,k_wp3) = 
     .    + 1.0 / dt
          ! New discretization for the turbulent advection term with a1 
          ! pulled outside of the derivative.  Brian.
!     .    + (1.0/2.0) * gr%dzt(k)
!     .      * (  a1(k) * ( wp3(kp1) + wp3(k) ) / MAX( wp2(k), eps )
!     .         - a1(km1) * ( wp3(k) + wp3(km1) ) / MAX( wp2(km1), eps )
!     .        )
     .    + (1.0/2.0) * ( a1(k) + a1(km1) )
     .      * (1.0/2.0) * gr%dzt(k)
     .      * (  ( wp3(kp1) + wp3(k) ) / MAX( wp2(k), eps )
     .         - ( wp3(k) + wp3(km1) ) / MAX( wp2(km1), eps )
     .        )
          ! End of Brian's a1 change.  Feb. 14, 2008.
     .    + ( 1.0 - C11_Skw_fnc(k) )
     .              * 3.0 * gr%dzt(k) * ( wmm(k) - wmm(km1) )
     .    + ( C8 / tauw3t(k) ) * ( 5.0 * C8b * Skwt(k)**4 + 1.0 )
     .    + diff_lhs(2)

          ! Momentum superdiagonal [ x wp2(k,<t+1>) ]
          lhs(3-1,k_wp3) = 
          ! Turbulent advection term with a3 outside of the derivative.
!     .    + 2.0 * gr%dzt(k) * a3(k) * wp2(k)   ! standard deriv.
!     .    + 2.0 * gr%dzt(k) * a3(k) * wp2(k)   ! Golaz's no a3 deriv.
     .    + 2.0 * gr%dzt(k) * 0.5*(a3(k)+a3(km1)) * wp2(k)   ! Centered no a3 deriv.
          ! End of code that pulls out a3.

          ! Thermodynamic superdiagonal [ x wp3(k+1,<t+1>) ]
          lhs(3-2,k_wp3) = 
     .    + wmt(k) * (1.0/2.0) * gr%dzt(k)
          ! New discretization for the turbulent advection term with a1 
          ! pulled outside of the derivative.  Brian.
!     .    + (1.0/2.0) * gr%dzt(k)
!     .      * a1(k) * ( wp3(kp1) + wp3(k) ) / MAX( wp2(k), eps )
     .    + (1.0/2.0) * ( a1(k) + a1(km1) )
     .      * (1.0/2.0) * gr%dzt(k)
     .      * ( wp3(kp1) + wp3(k) ) / MAX( wp2(k), eps )
          ! End of Brian's a1 change. Feb. 14, 2008.
     .    + diff_lhs(1)


          ! Right-hand side
          rhs(k_wp3) = 
     .    + ( 1.0 / dt ) * wp3(k)
          ! Turbulent advection term with a3 outside of the derivative.
     .    + gr%dzt(k)
!     .      * ( ( a3(k) * wp2(k)**2 ) - ( a3(km1) * wp2(km1)**2 ) ) ! standard deriv.
!     .      * ( ( a3(k) * wp2(k)**2 ) - ( a3(k) * wp2(km1)**2 ) )   ! Golaz's no a3 deriv.
     .      * (  ( 0.5*(a3(k)+a3(km1)) * wp2(k)**2 )        ! Centered no a3 deriv.
     .         - ( 0.5*(a3(k)+a3(km1)) * wp2(km1)**2 )  )   ! Centered no a3 deriv.
          ! End of code to pull out a3.
          ! New discretization for the turbulent advection term with a1 
          ! pulled outside of the derivative.  Brian.
!     .    + (1.0/4.0) * gr%dzt(k)
!     .      * (  a1(k) * ( wp3(kp1) + wp3(k) )**2 
!     .                   / MAX( wp2(k), eps )
!     .         - a1(km1) * ( wp3(k) + wp3(km1) )**2
!     .                     / MAX( wp2(km1), eps )
!     .        )
     .    + (1.0/2.0) * ( a1(k) + a1(km1) )
     .      * (1.0/4.0) * gr%dzt(k)
     .        * (  ( wp3(kp1) + wp3(k) )**2 
     .             / MAX( wp2(k), eps )
     .           - ( wp3(k) + wp3(km1) )**2
     .             / MAX( wp2(km1), eps )
     .          )
          ! End of Brian's a1 change. Feb. 14, 2008.
     .    + ( 1.0 - C11_Skw_fnc(k) ) * 3.0 * ( grav/T0 ) * wp2thvp(k)
     .    + ( C8 / tauw3t(k) ) * ( 4.0 * C8b * Skwt(k)**4 ) * wp3(k)

          if ( lcrank_nich_diff ) then
             ! These lines are for the diffusional term with a Crank-Nicholson 
             ! time step.  They are not used for completely implicit diffusion.
             rhs(k_wp3)   =   rhs(k_wp3)
     .                      - diff_lhs(3) * wp3(km1)
     .                      - diff_lhs(2) * wp3(k)
     .                      - diff_lhs(1) * wp3(kp1)
          endif


#ifdef STATS
          if (lstats_samp) then

!         Explicit contributions for wp3

          if ( iwp3_ta > 0 ) then
          zt%x(k,iwp3_ta) = zt%x(k,iwp3_ta) 
          ! Turbulent advection term with a3 outside of the derivative.
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
          ! End of code to pull out a3.
          ! New discretization for the turbulent advection term with a1 
          ! pulled outside of the derivative.  Brian.
!     .    + (1.0/4.0) * gr%dzt(k)
!     .      * (  a1(k) * ( wp3(kp1) + wp3(k) )**2 
!     .                   / MAX( wp2(k), eps )
!     .         - a1(km1) * ( wp3(k) + wp3(km1) )**2 
!     .                     / MAX( wp2(km1), eps )
!     .        )
     .    + (1.0/2.0) * ( a1(k) + a1(km1) )
     .      * (1.0/4.0) * gr%dzt(k)
     .        * (  ( wp3(kp1) + wp3(k) )**2 
     .             / MAX( wp2(k), eps )
     .           - ( wp3(k) + wp3(km1) )**2
     .             / MAX( wp2(km1), eps )
     .          )
          ! End of Brian's a1 change. Feb. 14, 2008.
          end if

          if ( iwp3_tp > 0 ) then
          zt%x(k,iwp3_tp) = zt%x(k,iwp3_tp) 
     .    - gr%dzt(k)
     .      * (  ( (3.0/2.0) * wp2(k)**2 ) 
     .         - ( (3.0/2.0) * wp2(km1)**2 )
     .        )
          end if

          if ( iwp3_bp > 0 ) then
          zt%x(k,iwp3_bp) = zt%x(k,iwp3_bp) 
     .    + 1.0 * 3.0 * ( grav / T0 )  * wp2thvp(k)
          zt%n(k,iwp3_bp) = zt%n(k,iwp3_bp) + 1
          end if

          if ( iwp3_pr2 > 0 ) then
          zt%x(k,iwp3_pr2) = zt%x(k,iwp3_pr2) 
     .    - C11_Skw_fnc(k) * 3.0 * ( grav / T0 ) * wp2thvp(k)
          end if

          if ( iwp3_pr1 > 0 ) then
          zt%x(k,iwp3_pr1) = zt%x(k,iwp3_pr1) 
     .    + ( C8 / tauw3t(k) ) * ( 4.0 * C8b * Skwt(k)**4 ) * wp3(k)
          end if

          ! These lines are for the diffusional term with a Crank-Nicholson 
          ! time step.  They need to be commented out for completely implicit
          ! diffusion to be used.
          if ( iwp3_dp1 > 0 ) then
           if ( lcrank_nich_diff ) then
            zt%x(k,iwp3_dp1) = zt%x(k,iwp3_dp1) 
     .      - diff_lhs(3) * wp3(km1)
     .      - diff_lhs(2) * wp3(k)
     .      - diff_lhs(1) * wp3(kp1)
           endif
          end if

!         Implicit contributions for wp3

          if ( iwp3_ta > 0 ) then
            ztscr05(k) =
            ! New discretization for the turbulent advection term with a1 
            ! pulled outside of the derivative.  Brian.
!     .      + (1.0/2.0) * gr%dzt(k)
!     .        * a1(km1) * ( wp3(k) + wp3(km1) ) / MAX( wp2(km1), eps )
     .        + (1.0/2.0) * ( a1(k) + a1(km1) )
     .          * (1.0/2.0) * gr%dzt(k)
     .          * ( wp3(k) + wp3(km1) ) / MAX( wp2(km1), eps )
            ! End of Brian's a1 change. Feb. 14, 2008.
            ztscr06(k) =
           ! Turbulent advection term with a3 outside of the derivative.
!     .    + 2.0 * gr%dzt(k) * ( a3(km1) + (3.0/2.0) ) * wp2(km1)  ! standard deriv.
!     .    + 2.0 * gr%dzt(k) * ( a3(k) + (3.0/2.0) ) * wp2(km1)    ! Golaz's no a3 deriv.
     .        + 2.0 * gr%dzt(k) * ( 0.5*(a3(k)+a3(km1)) + (3.0/2.0) ) ! Centered no a3 deriv.
     .          * wp2(km1)                                            ! Centered no a3 deriv.
            ! End of code to pull out a3.
            ztscr07(k) =
            ! New discretization for the turbulent advection term with a1 
            ! pulled outside of the derivative.  Brian.
!     .      - (1.0/2.0) * gr%dzt(k)
!     .        * ( a1(k) * ( wp3(kp1) + wp3(k) ) / MAX( wp2(k), eps )
!     .           -a1(km1) * ( wp3(k) + wp3(km1) ) / MAX( wp2(km1), eps )
!     .          )
     .        - (1.0/2.0) * ( a1(k) + a1(km1) )
     .          * (1.0/2.0) * gr%dzt(k)
     .          * (  ( wp3(kp1) + wp3(k) ) / MAX( wp2(k), eps )
     .             - ( wp3(k) + wp3(km1) ) / MAX( wp2(km1), eps )
     .            )
            ! End of Brian's a1 change. Feb. 14, 2008.
            ztscr08(k) =
            ! Turbulent advection term with a3 outside of the derivative.
!     .    - 2.0 * gr%dzt(k) * ( a3(k) + (3.0/2.0) ) * wp2(k)      ! standard deriv.
!     .    - 2.0 * gr%dzt(k) * ( a3(k) + (3.0/2.0) ) * wp2(k)      ! Golaz's no a3 deriv.
     .        - 2.0 * gr%dzt(k) * ( 0.5*(a3(k)+a3(km1)) + (3.0/2.0) ) ! Centered no a3 deriv.
     .          * wp2(k)                                              ! Centered no a3 deriv.
            ! End of code to pull out a3.
            ztscr09(k) =
            ! New discretization for the turbulent advection term with a1 
            ! pulled outside of the derivative.  Brian.
!     .      - (1.0/2.0) * gr%dzt(k)
!     .        * a1(k) * ( wp3(kp1) + wp3(k) ) / MAX( wp2(k), eps )
     .        - (1.0/2.0) * ( a1(k) + a1(km1) )
     .          * (1.0/2.0) * gr%dzt(k)
     .          * ( wp3(kp1) + wp3(k) ) / MAX( wp2(k), eps )
            ! End of Brian's a1 change. Feb. 14, 2008.
          end if

          if ( iwp3_tp > 0 ) then
            ztscr10(k) = 
     .      - 2.0 * gr%dzt(k) * (3.0/2.0) * wp2(km1)
            ztscr11(k) = 
     .      + 2.0 * gr%dzt(k) * (3.0/2.0) * wp2(k)
          end if

          if ( iwp3_ma > 0 ) then
            ztscr12(k) = 
     .      + wmt(k) * (1.0/2.0) * gr%dzt(k)
            ztscr13(k) = 
     .      - wmt(k) * (1.0/2.0) * gr%dzt(k)
          end if

          if ( iwp3_ac > 0 ) then
            ztscr14(k) = 
     .      - 1.0 * 3.0 * gr%dzt(k) * ( wmm(k) - wmm(km1) )
          end if

          if ( iwp3_pr2 > 0 ) then
            ztscr15(k) =
     .      + C11_Skw_fnc(k) * 3.0 * gr%dzt(k) * ( wmm(k) - wmm(km1) )
          end if

          if ( iwp3_pr1 > 0 ) then
          ztscr01(k) =
     .    - ( C8 / tauw3t(k) ) * ( 5.0 * C8b * Skwt(k)**4 + 1.0 )
          end if

          if ( iwp3_dp1 > 0 ) then
            ztscr02(k) = - diff_lhs(3)
            ztscr03(k) = - diff_lhs(2)
            ztscr04(k) = - diff_lhs(1)
          end if

          end if
#endif /*STATS*/

        end do

        ! Boundary conditions

        ! Both wp2 and wp3 used fixed-point boundary conditions.  Therefore,
        ! anything set in the above loop is overwritten at both the upper
        ! and lower boundaries.

        !   wp3(1)  wp2(1) ... wp3(nz) wp2(nz)
        ! [  0.0     0.0         0.0     0.0  ]
        ! [  0.0     0.0         0.0     0.0  ]
        ! [  1.0     1.0   ...   1.0     1.0  ]
        ! [  0.0     0.0         0.0     0.0  ]
        ! [  0.0     0.0         0.0     0.0  ]

        ! Lower boundary
        k = 1
        k_wp3 = 2*k - 1
        k_wp2 = 2*k

        lhs(:,k_wp2) = 0.0
        lhs(3,k_wp2) = 1.0
        rhs(k_wp2)   = wp2(k)
   
        lhs(:,k_wp3) = 0.0
        lhs(3,k_wp3) = 1.0
        rhs(k_wp3)   = 0.0

        ! Upper boundary
        k = gr%nnzp
        k_wp3 = 2*k - 1
        k_wp2 = 2*k

        lhs(:,k_wp2) = 0.0
        lhs(3,k_wp2) = 1.0
        rhs(k_wp2)   = 0.0
   
        lhs(:,k_wp3) = 0.0
        lhs(3,k_wp3) = 1.0
        rhs(k_wp3)   = 0.0

#ifdef STATS
        if ( lstats_samp .and. iwp23_cn > 0 ) then
          ! Perform LU decomp and solve system (LAPACK with diagnostics)
          call band_solvex( "wp23", nsup, nsub, 2*gr%nnzp, nrhs,
     .                      lhs, rhs, solut, rcond, isValid )

          ! Est. of the condition number of the w'^2/w^3 LHS matrix 
          sfc%x(1,iwp23_cn) = sfc%x(1,iwp23_cn) + 1.0 / rcond
          sfc%n(1,iwp23_cn) = sfc%n(1,iwp23_cn) + 1

        else
          ! Perform LU decomp and solve system (LAPACK)
          call band_solve( "wp23", nsup, nsub, 2*gr%nnzp, nrhs,
     .                     lhs, rhs, solut, isValid )
        end if
#else
        ! Perform LU decomp and solve system (LAPACK)
        call band_solve( "wp23", nsup, nsub, 2*gr%nnzp, nrhs,
     .                   lhs, rhs, solut, isValid )

#endif /*STATS*/

        if ( .not. isValid ) return


!       Copy result into output arrays and clip

        do k=1,gr%nnzp

          km1 = max( k-1, 1)
          kp1 = min( k+1, gr%nnzp )

          k_wp3 = 2*k - 1
          k_wp2 = 2*k

          wp2(k) = solut(k_wp2)
          wp3(k) = solut(k_wp3)

        end do

#ifdef STATS
        if (lstats_samp) then

!         Finalize implicit contributions for wp2

          do k=1,gr%nnzp
            km1 = max( k-1, 1 )
            kp1 = min( k+1, gr%nnzp )

            if ( iwp2_dp1 > 0 ) then
             zm%x(k,iwp2_dp1) = zm%x(k,iwp2_dp1)
     .       + zmscr01(k) * wp2(k)
             zm%n(k,iwp2_dp1) = zm%n(k,iwp2_dp1) + 1
            end if

            if ( iwp2_dp2 > 0 ) then
             zm%x(k,iwp2_dp2) = zm%x(k,iwp2_dp2)
     .       + zmscr02(k) * wp2(km1)
     .       + zmscr03(k) * wp2(k)
     .       + zmscr04(k) * wp2(kp1)
             zm%n(k,iwp2_dp2) = zm%n(k,iwp2_dp2) + 1
            end if

            if ( iwp2_ta > 0 ) then
             zm%x(k,iwp2_ta) = zm%x(k,iwp2_ta)
     .       + zmscr05(k) * wp3(k)
     .       + zmscr06(k) * wp3(kp1)
             zm%n(k,iwp2_ta) = zm%n(k,iwp2_ta) + 1
            end if

            if ( iwp2_ma > 0 ) then
             zm%x(k,iwp2_ma) = zm%x(k,iwp2_ma)
     .       + zmscr07(k) * wp2(km1)
     .       + zmscr08(k) * wp2(kp1)
             zm%n(k,iwp2_ma) = zm%n(k,iwp2_ma) + 1
            end if

            if ( iwp2_ac > 0 ) then
             zm%x(k,iwp2_ac) = zm%x(k,iwp2_ac)
     .       + zmscr09(k) * wp2(k)
             zm%n(k,iwp2_ac) = zm%n(k,iwp2_ac) + 1
            end if

            if ( iwp2_pr1 > 0 .and. lKhm_aniso ) then
              zm%x(k,iwp2_pr1) = zm%x(k,iwp2_pr1)
     .        + zmscr11(k) * wp2(k)
              zm%n(k,iwp2_pr1) = zm%n(k,iwp2_pr1) + 1
            end if

            if ( iwp2_pr2 > 0 ) then
             zm%x(k,iwp2_pr2) = zm%x(k,iwp2_pr2)
     .       + zmscr10(k) * wp2(k)
             zm%n(k,iwp2_pr2) = zm%n(k,iwp2_pr2) + 1
            end if

          end do

          ! Finalize implicit contributions for wp3

          do k=1,gr%nnzp
            km1 = max( k-1, 1 )
            kp1 = min( k+1, gr%nnzp )

            if ( iwp3_pr1 > 0 ) then
             zt%x(k,iwp3_pr1) = zt%x(k,iwp3_pr1)
     .       + ztscr01(k) * wp3(k)
             zt%n(k,iwp3_pr1) = zt%n(k,iwp3_pr1) + 1
            end if

            if ( iwp3_dp1 > 0 ) then
             zt%x(k,iwp3_dp1) = zt%x(k,iwp3_dp1)
     .       + ztscr02(k) * wp3(km1)
     .       + ztscr03(k) * wp3(k)
     .       + ztscr04(k) * wp3(kp1)
             zt%n(k,iwp3_dp1) = zt%n(k,iwp3_dp1) + 1
            end if

            if ( iwp3_ta > 0 ) then
             zt%x(k,iwp3_ta) = zt%x(k,iwp3_ta)
     .       + ztscr05(k) * wp3(km1)
     .       + ztscr06(k) * wp2(km1)
     .       + ztscr07(k) * wp3(k)
     .       + ztscr08(k) * wp2(k)
     .       + ztscr09(k) * wp3(kp1)
             zt%n(k,iwp3_ta) = zt%n(k,iwp3_ta) + 1
            end if

            if ( iwp3_tp > 0 ) then
             zt%x(k,iwp3_tp) = zt%x(k,iwp3_tp)
     .       + ztscr10(k) * wp2(km1)
     .       + ztscr11(k) * wp2(k)
             zt%n(k,iwp3_tp) = zt%n(k,iwp3_tp) + 1
            end if

            if ( iwp3_ma > 0 ) then
             zt%x(k,iwp3_ma) = zt%x(k,iwp3_ma)
     .       + ztscr12(k) * wp3(km1)
     .       + ztscr13(k) * wp3(kp1)
             zt%n(k,iwp3_ma) = zt%n(k,iwp3_ma) + 1
            end if

            if ( iwp3_ac > 0 ) then
             zt%x(k,iwp3_ac) = zt%x(k,iwp3_ac)
     .       + ztscr14(k) * wp3(k)
             zt%n(k,iwp3_ac) = zt%n(k,iwp3_ac) + 1
            end if

            if ( iwp3_pr2 > 0 ) then
             zt%x(k,iwp3_pr2) = zt%x(k,iwp3_pr2) 
     .       + ztscr15(k) * wp3(k)
             zt%n(k,iwp3_pr2) = zt%n(k,iwp3_pr2) + 1
            end if

          end do

          ! Effect of clipping term

          if ( iwp2_cl > 0 ) then
            zm%x(:,iwp2_cl) = zm%x(:,iwp2_cl) - wp2 / dt
          end if

          if ( iwp3_cl > 0 ) then
            zt%x(:,iwp3_cl) = zt%x(:,iwp3_cl) - wp3 / dt
          end if

        end if
#endif /*STATS*/

        ! Clip Skewness.
 
        DO k = 1, gr%nnzp, 1

           km1 = MAX(k-1,1)

!  Vince Larson commented out the Andre et al clipping to see if we 
!     could avoid using it.  26 Jul 2007
!  Brian and Vince undid the change because Wangara case still needs 
!  the Andre et al clipping.  27 Jul 2007.
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

           ! NOTE: Clipping not good; must try to get rid of.  Affects only
           !       Wangara case.

           IF ( gr%zt(k) <= 100.0 ) THEN ! Clip for 100 m. above ground.

              atmp = MAX( ( (wp2(km1)+wp2(k)) / 2.0 ) , eps )    ! wp2
              ctmp = atmp**(3.0/2.0)                             ! wp2^(3/2)
              IF ( wp3(k) >= 0.2 * SQRT(2.0) * ctmp ) THEN
                 wp3(k) = 0.2 * SQRT(2.0) * ctmp
              ! Vince Larson added clipping for negative wp3 
              !     in order to stabilize arm_3year    12 Dec 2007
              ELSE IF ( wp3(k) <= -0.2 * SQRT(2.0) * ctmp ) THEN
                 wp3(k) = -0.2 * SQRT(2.0) * ctmp
              ! End Vince Larson's addition
              END IF

           ELSE

              ! Clip skewness consistently with a
              ! NOTE: atmp is wp2 interpolated to thermodynamic levels
              !       (with the added insurance that it cannot be negative
              !       or zero).  wp2 should not be set equal to atmp
              !       because atmp is the value at a thermodynamic level,
              !       whereas wp2 is located at a momentum level.  The
              !       only effect of that piece of code is to mistakenly
              !       lift the wp2 value in altitude one-half grid box
              !       every time it is called.  The line is unnecessary,
              !       and can be commented out.  Brian Griffin. 6/30/07.
              atmp = MAX( 0.5*(wp2(km1)+wp2(k)), eps )
              ctmp = atmp**1.5
              IF ( wp3(k)/ctmp > 4.5 ) THEN
                 wp3(k) = 4.5 * ctmp
                 !wp2(k) = atmp
              ELSE IF ( wp3(k)/ctmp < -4.5 ) THEN
                 wp3(k) = -4.5 * ctmp
                 !wp2(k) = atmp
              END IF

           END IF
! End of Vince Larson's commenting of Andre et al clipping.

           wp2(k) = max( 2./3. * emin, wp2(k) ) 

        END DO

#ifdef STATS
        if (lstats_samp) then
          if ( iwp2_cl > 0 ) then
            zm%x(:,iwp2_cl) = zm%x(:,iwp2_cl) + wp2 / dt
            zm%n(:,iwp2_cl) = zm%n(:,iwp2_cl) + 1
          end if
          if ( iwp3_cl > 0 ) then
            zt%x(:,iwp3_cl) = zt%x(:,iwp3_cl) + wp3 / dt
            zt%n(:,iwp3_cl) = zt%n(:,iwp3_cl) + 1
          end if
        end if
#endif /*STATS*/

#ifdef STATS
        if (lstats_samp) then

          if ( iwp2_bt > 0 ) then
            zm%x(:,iwp2_bt) = zm%x(:,iwp2_bt) + wp2 / dt
            zm%n(:,iwp2_bt) = zm%n(:,iwp2_bt) + 1
          end if

          if ( iwp3_bt > 0 ) then
            zt%x(:,iwp3_bt) = zt%x(:,iwp3_bt) + wp3 / dt
            zt%n(:,iwp3_bt) = zt%n(:,iwp3_bt) + 1
          end if

        end if
#endif /*STATS*/

        return
        end subroutine wp23_solve

