!------------------------------------------------------------------------
! $Id: wp23.F,v 1.49 2008-04-10 21:48:35 griffinb Exp $
!========================================================================
        module wp23

        implicit none

        private

        public :: timestep_wp23

        private :: wp23_solve,
     .             wp2_term_ma_lhs,
     .             wp2_term_ta_lhs,
     .             wp2_terms_ac_pr2_lhs,
     .             wp2_term_dp1_lhs,
     .             wp2_term_pr1_lhs,
     .             wp2_terms_bp_pr2_rhs,
     .             wp2_term_pr3_rhs,
     .             wp2_term_pr1_rhs,
     .             wp3_term_ma_lhs,
     .             wp3_terms_ta_tp_lhs,
     .             wp3_terms_ac_pr2_lhs,
     .             wp3_term_pr1_lhs,
     .             wp3_terms_ta_tp_rhs,
     .             wp3_terms_bp_pr2_rhs,
     .             wp3_term_pr1_rhs

        contains

!========================================================================
        subroutine timestep_wp23
     .             ( dt, Lscale, taum, taut, Scm, Skwm, Skwt,
     .               a, Kht, Khm, wmt, wmm, up2, vp2,
     .               um, upwp, vm, vpwp, wpthvp, wp2thvp,
     .               wp2, wp3, err_code )
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
        use stats_prec, only: time_prec
        use error_code, only: CLUBB_singular_matrix,
     .                        CLUBB_bad_lapack_arg
#ifdef STATS
        use stats_hoc
#endif /*STATS*/

        implicit none

        intrinsic :: exp

        ! Input Variables
        real(kind=time_prec), intent(in) :: dt ! Timestep               [s]

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

        integer, intent(inout) :: err_code ! Diagnostic

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
     .                   wp2, wp3, err_code )
!        call wp23_solve( dt, Scm, Lscale, Lm, taum, tauw3t, 10*0.5*Kht,
!     .                   wmt, wmm, wpthvp, wp2thvp,
!     .                   um, upwp, vm, vpwp, up2, vp2,
!     .                   Skwt, C1_Skw_fnc, C11_Skw_fnc, 
!     .                   wp2, wp3, isValid )
! End Vince Larson's change

!       Error output
!       Joshua Fasching Feb 2008
        if ( err_code == CLUBB_singular_matrix .or.
     .       err_code == CLUBB_bad_lapack_arg ) then
                
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

!========================================================================
        subroutine wp23_solve( dt, Scm, Lscale, Lm, tau1m, tauw3t, Kw1,
     .                         Kw8, wmt, wmm, wpthvp, wp2thvp,
     .                         um, upwp, vm, vpwp, up2, vp2,
     .                         Skwt, C1_Skw_fnc, C11_Skw_fnc,
     .                         wp2, wp3, err_code)
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
        use stats_prec, only: time_prec
        use lapack_wrap, only: band_solve, band_solvex
        use error_code, only: CLUBB_singular_matrix,
     .                        CLUBB_bad_lapack_arg
#ifdef STATS
        use stats_hoc
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
        real(kind=time_prec), intent(in) :: 
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

        integer, intent(inout) :: err_code

        ! Local Variables
        real, dimension(nsup+nsub+1,2*gr%nnzp) :: 
     .  lhs ! Implicit contributions to wp2/wp3 (band diag. matrix)

        real, dimension(2*gr%nnzp) :: 
     .  rhs,  ! RHS of band matrix
     .  solut ! Solution to band diagonal system.

        real, dimension(gr%nnzp) :: 
     .  a1, ! a1; See eqn. 20 in `Equations for HOC'             [-]
     .  a3  ! a3; See eqn. 22 in `Equations for HOC'             [-]

        real, dimension(5) :: tmp

        real :: 
     .  atmp, ctmp, ! Temporary variables for the LHS
     .  rcond  ! Est. of the reciprocal of the condition #

        ! Array indices
        integer :: k, km1, kp1, k_wp2, k_wp3

        ! Added by Brian for eddy diffusion.
        real, dimension(3) :: rhs_diff
        ! Set logical to true for Crank-Nicholson diffusion scheme
        ! or to false for completely implicit diffusion scheme.
        ! Note:  Although Crank-Nicholson diffusion has usually been
        !        used for wp2 and wp3 in the past, we found that using
        !        completely implicit diffusion stabilized the deep
        !        convective cases more while having almost no effect on
        !        the boundary layer cases.  Brian; 1/4/2008.
!        logical, parameter :: lcrank_nich_diff = .true.
        logical, parameter :: lcrank_nich_diff = .false.

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


          !!!!!***** w'^2 *****!!!!!

          ! w'^2: Left-hand side (implicit w'^2 portion of the code).
          !
          ! Momentum subdiagonal (lhs index: 3+2)
          !         [ x wp2(k-1,<t+1>) ]
          ! Thermodynamic subdiagonal (lhs index: 3+1)
          !         [ x wp3(k,<t+1>) ]
          ! Momentum main diagonal (lhs index: 3)
          !         [ x wp2(k,<t+1>) ]
          ! Thermodynamic superdiagonal (lhs index: 3-1)
          !         [ x wp3(k+1,<t+1>) ]
          ! Momentum superdiagonal (lhs index: 3-2)
          !         [ x wp2(k+1,<t+1>) ]

          ! LHS time tendency.
          lhs(3,k_wp2)
     .    = + 1.0 / dt

          ! LHS mean advection (ma) term.
          lhs((/3-2,3+2/),k_wp2)
     .    = lhs((/3-2,3+2/),k_wp2)
     .    + wp2_term_ma_lhs( wmm(k), gr%dzm(k) )

          ! LHS turbulent advection (ta) term.
          lhs((/3-1,3+1/),k_wp2)
     .    = lhs((/3-1,3+1/),k_wp2)
     .    + wp2_term_ta_lhs( gr%dzm(k) )

          ! LHS accumulation (ac) term and pressure term 2 (pr2).
          lhs(3,k_wp2)
     .    = lhs(3,k_wp2)
     .    + wp2_terms_ac_pr2_lhs( C5, wmt(kp1), wmt(k), gr%dzm(k)  )

          ! LHS dissipation term 1 (dp1).
          lhs(3,k_wp2)
     .    = lhs(3,k_wp2)
     .    + wp2_term_dp1_lhs( C1_Skw_fnc(k), tau1m(k) )

          ! LHS eddy diffusion term: dissipation term 2 (dp2).
          if ( lcrank_nich_diff ) then
             ! Eddy diffusion for wp2 using a Crank-Nicholson time step.
             lhs((/3-2,3,3+2/),k_wp2)
     .       = lhs((/3-2,3,3+2/),k_wp2)
     .       + (1.0/2.0)
     .       * diffusion_mlev( Kw1(k), Kw1(kp1), nu1,
     .                         gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )
          else
             ! Eddy diffusion for wp2 using a completely implicit time step. 
             lhs((/3-2,3,3+2/),k_wp2)
     .       = lhs((/3-2,3,3+2/),k_wp2)
     .       + diffusion_mlev( Kw1(k), Kw1(kp1), nu1,
     .                         gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )
          endif

          ! LHS pressure term 1 (pr1).
          if ( lKhm_aniso ) then
             ! Add in this term if we're not assuming tke = 1.5 * wp2
             lhs(3,k_wp2)
     .       = lhs(3,k_wp2)
     .       + wp2_term_pr1_lhs( C4, tau1m(k) )
          endif


          ! w'^2: Right-hand side (explicit w'^2 portion of the code).

          ! RHS time tendency.
          rhs(k_wp2)
     .    = + ( 1.0 / dt ) * wp2(k) 

          ! RHS buoyancy production (bp) term and pressure term 2 (pr2).
          rhs(k_wp2)
     .    = rhs(k_wp2)
     .    + wp2_terms_bp_pr2_rhs( C5, wpthvp(k) )

          ! RHS pressure term 3 (pr3).
          rhs(k_wp2)
     .    = rhs(k_wp2)
     .    + wp2_term_pr3_rhs( C5, wpthvp(k), upwp(k), um(kp1), um(k),
     .                        vpwp(k), vm(kp1), vm(k), gr%dzm(k) )

          ! RHS eddy diffusion term: dissipation term 2 (dp2).
          if ( lcrank_nich_diff ) then
             ! These lines are for the diffusional term with a Crank-Nicholson 
             ! time step.  They are not used for completely implicit diffusion.
             rhs_diff(1:3)
     .       = (1.0/2.0)
     .       * diffusion_mlev( Kw1(k), Kw1(kp1), nu1,
     .                         gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )
             rhs(k_wp2)   =   rhs(k_wp2)
     .                      - rhs_diff(3) * wp2(km1)
     .                      - rhs_diff(2) * wp2(k)
     .                      - rhs_diff(1) * wp2(kp1)
          endif

          ! RHS pressure term 1 (pr1).
          if ( lKhm_aniso ) then
             rhs(k_wp2)
     .       = rhs(k_wp2)
     .       + wp2_term_pr1_rhs( C4, up2(k), vp2(k), tau1m(k) )
          endif


#ifdef STATS
          if ( lstats_samp ) then

!          Explicit contributions for wp2

           if ( iwp2_dp2 > 0 ) then
            if ( lcrank_nich_diff ) then
             zm%x(k,iwp2_dp2) = zm%x(k,iwp2_dp2) 
     .       - rhs_diff(3) * wp2(km1)
     .       - rhs_diff(2) * wp2(k)
     .       - rhs_diff(1) * wp2(kp1)
            endif
           end if

           if ( iwp2_bp > 0 ) then
            zm%x(k,iwp2_bp) = zm%x(k,iwp2_bp) 
     .      + wp2_terms_bp_pr2_rhs( 0.0, wpthvp(k) )
            zm%n(k,iwp2_bp) = zm%n(k,iwp2_bp) + 1
           end if

           if ( iwp2_pr1 > 0 .and. lKhm_aniso ) then
            zm%x(k,iwp2_pr1) = zm%x(k,iwp2_pr1) 
     .      + wp2_term_pr1_rhs( C4, up2(k), vp2(k), tau1m(k) )
           end if

           if ( iwp2_pr2 > 0 ) then
            zm%x(k,iwp2_pr2) = zm%x(k,iwp2_pr2) 
     .      + wp2_terms_bp_pr2_rhs( (1.0+C5), wpthvp(k) )
           end if

           if ( iwp2_pr3 > 0 ) then
            zm%x(k,iwp2_pr3) = zm%x(k,iwp2_pr3) 
     .      + wp2_term_pr3_rhs( C5, wpthvp(k), upwp(k), um(kp1), um(k),
     .                          vpwp(k), vm(kp1), vm(k), gr%dzm(k) )
            zm%n(k,iwp2_pr3) = zm%n(k,iwp2_pr3) + 1
           end if

           ! Implicit contributions for wp2

           if ( iwp2_dp1 > 0 ) then
             zmscr01(k) =
     .       - wp2_term_dp1_lhs( C1_Skw_fnc(k), tau1m(k) )
           end if

           if ( iwp2_dp2 > 0 ) then
             if ( lcrank_nich_diff ) then
                ! Eddy diffusion for wp2 using a Crank-Nicholson time step.
                tmp(1:3)
     .          = (1.0/2.0)
     .          * diffusion_mlev( Kw1(k), Kw1(kp1), nu1,
     .                            gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )
             else
                ! Eddy diffusion for wp2 using a completely implicit time step. 
                tmp(1:3)
     .          = diffusion_mlev( Kw1(k), Kw1(kp1), nu1,
     .                            gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )
             endif
             zmscr02(k) = - tmp(3)
             zmscr03(k) = - tmp(2)
             zmscr04(k) = - tmp(1)
           end if

           if ( iwp2_ta > 0 ) then
             tmp(1:2) = 
     .       + wp2_term_ta_lhs( gr%dzm(k) )
             zmscr05(k) = - tmp(2)
             zmscr06(k) = - tmp(1)
           end if

           if ( iwp2_ma > 0 ) then
             tmp(1:2) =
     .       + wp2_term_ma_lhs( wmm(k), gr%dzm(k) )
             zmscr07(k) = - tmp(2)
             zmscr08(k) = - tmp(1)
           end if

           if ( iwp2_ac > 0 ) then
             zmscr09(k) = 
     .       - wp2_terms_ac_pr2_lhs( 0.0, wmt(kp1), wmt(k), gr%dzm(k)  )
           end if

           if ( iwp2_pr2 > 0 ) then
             zmscr10(k) = 
     .       - wp2_terms_ac_pr2_lhs( (1.0+C5), wmt(kp1), wmt(k), 
     .                               gr%dzm(k)  )
           end if

           if ( iwp2_pr1 > 0 .and. lKhm_aniso ) then
             zmscr11(k) = - wp2_term_pr1_lhs( C4, tau1m(k) )
           end if

         end if
#endif /*STATS*/


          !!!!!***** w'^3 *****!!!!!

          ! w'^3: Left-hand side (implicit w'^3 portion of the code).
          !
          ! Thermodynamic subdiagonal (lhs index: 3+2) 
          !         [ x wp3(k-1,<t+1>) ]
          ! Momentum subdiagonal (lhs index: 3+1)
          !         [ x wp2(k-1,<t+1>) ]
          ! Thermodynamic main diagonal (lhs index: 3)
          !         [ x wp3(k,<t+1>) ]
          ! Momentum superdiagonal (lhs index: 3-1)
          !         [ x wp2(k,<t+1>) ]
          ! Thermodynamic superdiagonal (lhs index: 3-2) 
          !         [ x wp3(k+1,<t+1>) ]

          ! LHS time tendency.
          lhs(3,k_wp3)
     .    = + 1.0 / dt

          ! LHS mean advection (ma) term.
          lhs((/3-2,3+2/),k_wp3)
     .    = lhs((/3-2,3+2/),k_wp3)
     .    + wp3_term_ma_lhs( wmt(k), gr%dzt(k) )

          ! LHS turbulent advection (ta) and turbulent production (tp) terms.
          lhs(3-2:3+2,k_wp3)
     .    = lhs(3-2:3+2,k_wp3)
     .    + wp3_terms_ta_tp_lhs( wp3(kp1), wp3(k), wp3(km1), 
     .                           wp2(k), wp2(km1), a1(k), a1(km1),
     .                           a3(k), a3(km1), gr%dzt(k), eps, wtol )

          ! LHS accumulation (ac) term and pressure term 2 (pr2).
          lhs(3,k_wp3)
     .    = lhs(3,k_wp3)
     .    + wp3_terms_ac_pr2_lhs( C11_Skw_fnc(k),
     .                            wmm(k), wmm(km1), gr%dzt(k) )

          ! LHS pressure term 1 (pr1).
          lhs(3,k_wp3)
     .    = lhs(3,k_wp3)
     .    + wp3_term_pr1_lhs( C8, C8b, tauw3t(k), Skwt(k) )

          ! LHS eddy diffusion term: dissipation term 1 (dp1).
          !  Added a new constant, C12.
          !  Initially, this new constant will be set to 1.0 -dschanen 9/19/05 
          if ( lcrank_nich_diff ) then
             ! Eddy diffusion for wp3 using a Crank-Nicholson time step. 
             lhs((/3-2,3,3+2/),k_wp3)
     .       = lhs((/3-2,3,3+2/),k_wp3)
     .       + C12 * (1.0/2.0)
     .       * diffusion_tlev( Kw8(k), Kw8(km1), nu8,
     .                         gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
          else
             ! Eddy diffusion for wp3 using a completely implicit time step. 
             lhs((/3-2,3,3+2/),k_wp3)
     .       = lhs((/3-2,3,3+2/),k_wp3)
     .       + C12 
     .       * diffusion_tlev( Kw8(k), Kw8(km1), nu8,
     .                         gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
          endif


          ! w'^3: Right-hand side (explicit w'^3 portion of the code).

          ! RHS time tendency.
          rhs(k_wp3) = 
     .    + ( 1.0 / dt ) * wp3(k)

          ! RHS turbulent advection (ta) and turbulent production (tp) terms.
          rhs(k_wp3)
     .    = rhs(k_wp3)
     .    + wp3_terms_ta_tp_rhs( wp3(kp1), wp3(k), wp3(km1), 
     .                           wp2(k), wp2(km1), a1(k), a1(km1),
     .                           a3(k), a3(km1), gr%dzt(k), eps, wtol )

          ! RHS buoyancy production (bp) term and pressure term 2 (pr2).
          rhs(k_wp3)
     .    = rhs(k_wp3)
     .    + wp3_terms_bp_pr2_rhs( C11_Skw_fnc(k), wp2thvp(k) )

          ! RHS pressure term 1 (pr1).
          rhs(k_wp3)
     .    = rhs(k_wp3)
     .    + wp3_term_pr1_rhs( C8, C8b, tauw3t(k), Skwt(k), wp3(k) )

          ! RHS eddy diffusion term: dissipation term 1 (dp1).
          if ( lcrank_nich_diff ) then
             ! These lines are for the diffusional term with a Crank-Nicholson 
             ! time step.  They are not used for completely implicit diffusion.
             rhs_diff(1:3)
     .       = C12 * (1.0/2.0)
     .       * diffusion_tlev( Kw8(k), Kw8(km1), nu8,
     .                         gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
             rhs(k_wp3)   =   rhs(k_wp3)
     .                      - rhs_diff(3) * wp3(km1)
     .                      - rhs_diff(2) * wp3(k)
     .                      - rhs_diff(1) * wp3(kp1)
          endif


#ifdef STATS
          if (lstats_samp) then

!         Explicit contributions for wp3

          if ( iwp3_ta > 0 ) then
          zt%x(k,iwp3_ta) = zt%x(k,iwp3_ta) 
     .    + wp3_terms_ta_tp_rhs( wp3(kp1), wp3(k), wp3(km1), 
     .                           wp2(k), wp2(km1), a1(k), a1(km1),
     .                           a3(k)+(3.0/2.0), a3(km1)+(3.0/2.0), 
     .                           gr%dzt(k), eps, wtol )
          end if

          if ( iwp3_tp > 0 ) then
          zt%x(k,iwp3_tp) = zt%x(k,iwp3_tp) 
     .    + wp3_terms_ta_tp_rhs( wp3(kp1), wp3(k), wp3(km1), 
     .                           wp2(k), wp2(km1), 0.0, 0.0,
     .                           0.0-(3.0/2.0), 0.0-(3.0/2.0), 
     .                           gr%dzt(k), eps, wtol )
          end if

          if ( iwp3_bp > 0 ) then
          zt%x(k,iwp3_bp) = zt%x(k,iwp3_bp) 
     .    + wp3_terms_bp_pr2_rhs( 0.0, wp2thvp(k) )
          zt%n(k,iwp3_bp) = zt%n(k,iwp3_bp) + 1
          end if

          if ( iwp3_pr2 > 0 ) then
          zt%x(k,iwp3_pr2) = zt%x(k,iwp3_pr2) 
     .    + wp3_terms_bp_pr2_rhs( (1.0+C11_Skw_fnc(k)), wp2thvp(k) )
          end if

          if ( iwp3_pr1 > 0 ) then
          zt%x(k,iwp3_pr1) = zt%x(k,iwp3_pr1) 
     .    + wp3_term_pr1_rhs( C8, C8b, tauw3t(k), Skwt(k), wp3(k) )
          end if

          ! These lines are for the diffusional term with a Crank-Nicholson 
          ! time step.  They need to be commented out for completely implicit
          ! diffusion to be used.
          if ( iwp3_dp1 > 0 ) then
           if ( lcrank_nich_diff ) then
            zt%x(k,iwp3_dp1) = zt%x(k,iwp3_dp1) 
     .      - rhs_diff(3) * wp3(km1)
     .      - rhs_diff(2) * wp3(k)
     .      - rhs_diff(1) * wp3(kp1)
           endif
          end if

!         Implicit contributions for wp3

          if ( iwp3_ta > 0 ) then
            tmp(1:5) = 
     .      wp3_terms_ta_tp_lhs( wp3(kp1), wp3(k), wp3(km1), 
     .                           wp2(k), wp2(km1), a1(k), a1(km1),
     .                           a3(k)+(3.0/2.0), a3(km1)+(3.0/2.0), 
     .                           gr%dzt(k), eps, wtol )
            ztscr05(k) = -tmp(5)
            ztscr06(k) = -tmp(4)
            ztscr07(k) = -tmp(3)
            ztscr08(k) = -tmp(2)
            ztscr09(k) = -tmp(1)
          end if

          if ( iwp3_tp > 0 ) then
            tmp(1:5) = 
     .      wp3_terms_ta_tp_lhs( wp3(kp1), wp3(k), wp3(km1), 
     .                           wp2(k), wp2(km1), 0.0, 0.0,
     .                           0.0-(3.0/2.0), 0.0-(3.0/2.0), 
     .                           gr%dzt(k), eps, wtol )
            ztscr10(k) = -tmp(4)
            ztscr11(k) = -tmp(2)
          end if

          if ( iwp3_ma > 0 ) then
            tmp(1:2) =
     .      wp3_term_ma_lhs( wmt(k), gr%dzt(k) )
            ztscr12(k) = -tmp(2)
            ztscr13(k) = -tmp(1)
          end if

          if ( iwp3_ac > 0 ) then
            ztscr14(k) = 
     .      - wp3_terms_ac_pr2_lhs( 0.0,
     .                              wmm(k), wmm(km1), gr%dzt(k) )
          end if

          if ( iwp3_pr2 > 0 ) then
            ztscr15(k) =
     .      - wp3_terms_ac_pr2_lhs( (1.0+C11_Skw_fnc(k)),
     .                              wmm(k), wmm(km1), gr%dzt(k) )
          end if

          if ( iwp3_pr1 > 0 ) then
          ztscr01(k) =
     .    - wp3_term_pr1_lhs( C8, C8b, tauw3t(k), Skwt(k) )
          end if

          if ( iwp3_dp1 > 0 ) then
            if ( lcrank_nich_diff ) then
               ! Eddy diffusion for wp3 using a Crank-Nicholson time step. 
               tmp(1:3)
     .         = C12 * (1.0/2.0)
     .         * diffusion_tlev( Kw8(k), Kw8(km1), nu8,
     .                           gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
            else
               ! Eddy diffusion for wp3 using a completely implicit time step. 
               tmp(1:3)
     .         = C12 
     .         * diffusion_tlev( Kw8(k), Kw8(km1), nu8,
     .                           gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
            endif
            ztscr02(k) = - tmp(3)
            ztscr03(k) = - tmp(2)
            ztscr04(k) = - tmp(1)
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
     .                      lhs, rhs, solut, rcond, err_code )

          ! Est. of the condition number of the w'^2/w^3 LHS matrix 
          sfc%x(1,iwp23_cn) = sfc%x(1,iwp23_cn) + 1.0 / rcond
          sfc%n(1,iwp23_cn) = sfc%n(1,iwp23_cn) + 1

        else
          ! Perform LU decomp and solve system (LAPACK)
          call band_solve( "wp23", nsup, nsub, 2*gr%nnzp, nrhs,
     .                     lhs, rhs, solut, err_code )
        end if
#else
        ! Perform LU decomp and solve system (LAPACK)
        call band_solve( "wp23", nsup, nsub, 2*gr%nnzp, nrhs,
     .                   lhs, rhs, solut, err_code )

#endif /*STATS*/

        if ( err_code == CLUBB_singular_matrix .or.
     .        err_code == CLUBB_bad_lapack_arg ) return


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

          do k=2,gr%nnzp-1
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

          do k=2, gr%nnzp-1, 1
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

!========================================================================
        pure function wp2_term_ma_lhs( wmm, dzm )
     .  result( lhs )

!       Description:
!       Mean advection of w'^2: implicit portion of the code.
!
!       The values of w'^2 are found on the momentum levels, as are the 
!       values of wmm (mean vertical velocity on momentum levels).  
!       The variable w'^2 is interpolated to the intermediate
!       thermodynamic levels.  The derivative of the interpolated values
!       is taken over the central momentum level.  The derivative
!       is multiplied by wmm at the central momentum level to get
!       the desired result.  The w'^3 term is solved for implicity.
!
!       ============wp2p1======================================== m(k+1)
!
!       --------------------wp2(interp)-------------------------- t(k+1)
!
!       ============wp2===================d(wp2)/dz========wmm=== m(k)
!
!       --------------------wp2(interp)-------------------------- t(k)
!
!       ============wp2m1======================================== m(k-1)
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

        ! Momentum superdiagonal: [ x wp2(k+1,<ts+1>) ]
        lhs(kp1_mdiag)
     .  = + wmm * (1.0/2.0) * dzm

        ! Momentum subdiagonal: [ x wp2(k-1,<ts+1>) ]
        lhs(km1_mdiag)
     .  = - wmm * (1.0/2.0) * dzm

        return
        end function wp2_term_ma_lhs

!========================================================================
        pure function wp2_term_ta_lhs( dzm )
     .  result( lhs )

!       Description:
!       Turbulent advection term for w'^2: implicit portion of the code.
!
!       While the values of w'^2 are found on the momentum levels, the 
!       values of w'^3 are found on the thermodynamic levels.  The 
!       derivative of w'^3 is taken over the intermediate (central) 
!       momentum level, yielding the desired results.  The value of 
!       w'^3 being used is the value at timestep (ts+1), thus implicit
!       discretization is used.
!
!       -------------------wp3p1--------------------------------- t(k+1)
!
!       =============================d(wp3)/dz=================== m(k)
!
!       -------------------wp3----------------------------------- t(k)
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
     .  k_tdiag   = 2    ! Thermodynamic subdiagonal index.

        ! Input Variables
        real, intent(in) ::
     .  dzm     ! Inverse of grid spacing (k)   [1/m]

        ! Return Variable
        real, dimension(2) :: lhs

        ! Thermodynamic superdiagonal: [ x wp3(k+1,<ts+1>) ]
        lhs(kp1_tdiag)
     .  = + dzm

        ! Thermodynamic subdiagonal: [ x wp3(k,<ts+1>) ]
        lhs(k_tdiag)
     .  = - dzm

        return
        end function wp2_term_ta_lhs

!========================================================================
        pure function wp2_terms_ac_pr2_lhs( C5, wmtp1, wmt, dzm )
     .  result( lhs )

!       Description:
!       Accumulation of w'^2 and w'^2 pressure term 2:
!       implicit portion of the code.
!
!       The w'^2 accumulation term is completely implicit, while w'^2
!       pressure term 2 has both implicit and explicit components.  The
!       accumulation term and the implicit portion of pressure term 2 
!       are combined and solved together.
!
!       The values of w'^2 are found on momentum levels, while the values
!       of wmt (mean vertical velocity on thermodynamic levels) are found
!       on thermodynamic levels.  The vertical derivative of wmt is 
!       taken over the intermediate (central) momentum level.  It is then 
!       multiplied by w^'2 (implicitly calculated at timestep (ts+1)) and 
!       the coefficients to yield the desired results.
!
!       -------wmtp1--------------------------------------------- t(k+1)
!
!       ===============d(wmt)/dz============wp2================== m(k)
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
     .  C5,     ! Model parameter C_5           [-]
     .  wmtp1,  ! wmt(k+1)                      [m/s]
     .  wmt,    ! wmt(k)                        [m/s]
     .  dzm     ! Inverse of grid spacing (k)   [1/m]

        ! Return Variable
        real :: lhs

        ! Momentum main diagonal: [ x wp2(k,<ts+1>) ]
        lhs
     .  = + ( 1.0 - C5 ) * 2.0 * dzm * ( wmtp1 - wmt )

        return
        end function wp2_terms_ac_pr2_lhs

!========================================================================
        pure function wp2_term_dp1_lhs( C1_Skw_fnc, tau1m )
     .  result( lhs )

!       Description:
!       Dissipation term 1 for w'^2: implicit portion of the code.
!
!       The values of w'^2 are found on momentum levels.  The values of 
!       the C_1 skewness function and time-scale tau1m are also found on
!       momentum levels.  The w'^2 dissipation term 1 is completely 
!       implicit.

!       References:
!-----------------------------------------------------------------------

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  k_mdiag   = 1    ! Momentum main diagonal index.

        ! Input Variables
        real, intent(in) ::
     .  C1_Skw_fnc, ! C_1 skewness function (k)            [-]
     .  tau1m       ! Time scale at momentum levels        [s]

        ! Return Variable
        real :: lhs

        ! Momentum main diagonal: [ x wp2(k,<ts+1>) ]
        lhs
     .  = + C1_Skw_fnc / tau1m

        return
        end function wp2_term_dp1_lhs

!========================================================================
        pure function wp2_term_pr1_lhs( C4, tau1m )
     .  result( lhs )

!       Description
!       Pressure term 1 for w'^2: implicit portion of the code.
!
!       Pressure term 1 has both implicit and explicit components.  The
!       values of w'^2 are found on momentum levels, as are the values 
!       of tau1m.  For this portion of the code, w'^2 is solved for 
!       implicitly.

!       References:
!-----------------------------------------------------------------------

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  k_mdiag   = 1    ! Momentum main diagonal index.

        ! Input Variables
        real, intent(in) ::
     .  C4,     ! Model parameter C_4                  [-]
     .  tau1m   ! Time scale at momentum levels        [s]

        ! Return Variable
        real :: lhs

        ! Momentum main diagonal: [ x wp2(k,<ts+1>) ]
        lhs
     .  = + ( 2.0 * C4 ) / ( 3.0 * tau1m )

        return
        end function wp2_term_pr1_lhs

!========================================================================
        pure function wp2_terms_bp_pr2_rhs( C5, wpthvp )
     .  result( rhs )

!       Description:
!       Buoyancy production of w'^2 and w'^2 pressure term 2:
!       explicit portion of the code.
!
!       The w'^2 buoyancy production term is completely explicit, while 
!       w'^2 pressure term 2 has both implicit and explicit components.  
!       The buoyancy production term and the explicit portion of pressure
!       term 2 are combined and solved together.

!       References:
!-----------------------------------------------------------------------

        use constants, only: grav ! Gravitational acceleration [m/s^2]
        use parameters, only: T0  ! Reference temperature      [K]

        implicit none

        ! Input Variables
        real, intent(in) ::
     .  C5,      ! Model parameter C_5 [-]
     .  wpthvp   ! w'th_v'(k)          [K m/s]

        ! Return Variable
        real :: rhs

        rhs
     .  = + ( 1.0 - C5 ) * 2.0 * ( grav / T0 ) * wpthvp

        return
        end function wp2_terms_bp_pr2_rhs

!========================================================================
        pure function wp2_term_pr3_rhs( C5, wpthvp, upwp, ump1, um,
     .                                  vpwp, vmp1, vm, dzm )
     .  result( rhs )

!       Description:
!       Pressure term 3 for w'^2: explicit portion of the code.
!
!       The values of w'th_v', u'w', and v'w' are found on the momentum
!       levels, whereas the values of um and vm are found on the 
!       thermodynamic levels.  The derivatives of both um and vm are 
!       taken over the intermediate (central) momentum level.  All the
!       remaining mathematical operations take place at the central 
!       momentum level, yielding the desired result.  The results 
!       contributes the explicit portion of the w'^2 solution.
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

        use constants, only: grav ! Gravitational acceleration [m/s^2]
        use parameters, only: T0  ! Reference temperature      [K]

        implicit none

        ! Input Variables
        real, intent(in) ::
     .  C5,     ! Model parameter C_5           [-]
     .  wpthvp, ! w'th_v'(k)                    [K m/s]
     .  upwp,   ! u'w'(k)                       [m^2/s^2]
     .  ump1,   ! um(k+1)                       [m/s]
     .  um,     ! um(k)                         [m/s]
     .  vpwp,   ! v'w'(k)                       [m^2/s^2]
     .  vmp1,   ! vm(k+1)                       [m/s]
     .  vm,     ! vm(k)                         [m/s]
     .  dzm     ! Inverse of grid spacing (k)   [1/m]

        ! Return Variable
        real :: rhs

        rhs
! Michael Falk, 2 August 2007
! Use the following code for standard mixing, with c_k=0.548:
!     .  = + (2.0/3.0) * C5
!     .                * ( ( grav / T0 ) * wpthvp
!     .                    - upwp * dzm * ( ump1 - um )
!     .                    - vpwp * dzm * ( vmp1 - vm )
!     .                  )
! Use the following code for alternate mixing, with c_k=0.1 or 0.2
     .  = + (2.0/3.0) * C5
     .                * ( ( grav / T0 ) * wpthvp
     .                    - 0. * upwp * dzm * ( ump1 - um )
     .                    - 0. * vpwp * dzm * ( vmp1 - vm )
     .                  )
! eMFc

        return
        end function wp2_term_pr3_rhs

!========================================================================
        pure function wp2_term_pr1_rhs( C4, up2, vp2, tau1m )
     .  result( rhs )

!       Description:
!       Pressure term 1 for w'^2: explicit portion of the code.
!
!       Pressure term 1 has both implicit and explicit components.  The
!       values of u'^2 and v'^2 are found on momentum levels, as are the
!       values of tau1m.  This portion of the code contributes to the 
!       explicit part of the w'^2 solution.

!       References:
!-----------------------------------------------------------------------

        implicit none

        ! Input Variables
        real, intent(in) ::
     .  C4,     ! Model parameter C_4                  [-]
     .  up2,    ! u'^2(k)                              [m^2/s^2]
     .  vp2,    ! v'^2(k)                              [m^2/s^2]
     .  tau1m   ! Time scale at momentum levels        [s]

        ! Return Variable
        real :: rhs

        rhs
     .  = + ( C4 * ( up2 + vp2 ) ) / ( 3.0 * tau1m )

        return
        end function wp2_term_pr1_rhs

!========================================================================
        pure function wp3_term_ma_lhs( wmt, dzt )
     .  result( lhs )

!       Description:
!       Mean advection of w'^3: implicit portion of the code.
!
!       The values of w'^3 are found on the thermodynamic levels, as 
!       are the values of wmt (mean vertical velocity on thermodynamic 
!       levels).  The variable w'^3 is interpolated to the intermediate
!       momentum levels.  The derivative of the interpolated values is
!       taken over the central thermodynamic level.  The derivative 
!       is multiplied by wmt at the central thermodynamic level to get 
!       the desired result.  The w'^3 term is solved for implicity.
!
!       ------------wp3p1---------------------------------------- t(k+1)
!
!       ====================wp3(interp)========================== m(k)
!
!       ------------wp3-------------------d(wp3)/dz--------wmt--- t(k)
!
!       ====================wp3(interp)========================== m(k-1)
!
!       ------------wp3m1---------------------------------------- t(k-1)
!
!       The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1)
!       correspond with altitudes zt(k+1), zm(k), zt(k), zm(k-1),
!       and zt(k-1), respectively.  The letter "t" is used for
!       thermodynamic levels and the letter "m" is used for momentum
!       levels.
!
!       dzt(k) = 1 / ( zm(k) - zm(k-1) )

!       References:
!-----------------------------------------------------------------------

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  kp1_tdiag = 1,   ! Thermodynamic superdiagonal index.
     .  km1_tdiag = 2    ! Thermodynamic subdiagonal index.

        ! Input Variables
        real, intent(in) :: 
     .  wmt,  ! wmt(k)                        [m/s]
     .  dzt   ! Inverse of grid spacing (k)   [1/m]

        ! Return Variable
        real, dimension(2) :: lhs

        ! Thermodynamic superdiagonal: [ x wp3(k+1,<ts+1>) ]
        lhs(kp1_tdiag)
     .  = + wmt * (1.0/2.0) * dzt

        ! Thermodynamic subdiagonal: [ x wp3(k-1,<ts+1>) ]
        lhs(km1_tdiag)
     .  = - wmt * (1.0/2.0) * dzt

        return
        end function wp3_term_ma_lhs

!========================================================================
        pure function wp3_terms_ta_tp_lhs( wp3p1, wp3, wp3m1, 
     .                                     wp2, wp2m1, a1, a1m1,
     .                                     a3, a3m1, dzt, eps, wtol )
     .  result( lhs )

!       Description:
!       Turbulent advection and turbulent production of w'^3:  
!       implicit portion of the code.
!
!       The values of w'^3 are found on the thermodynamic levels, while
!       the values of w'^2, a_1, and a_3 are found on the momentum 
!       levels.  The variable w'^3 is interpolated to the intermediate 
!       momentum levels.  The values of the mathematical expressions
!       (called X and Y here) within the dX/dz and dY/dz terms are 
!       computed on the momentum levels.  Then, the derivatives (d/dz) 
!       of the expressions (X and Y) are taken over the central 
!       thermodynamic level, yielding the desired result.
!       In this function, the values of X and Y are as follows:
!
!       X = a_3(ts) * 2 * w'^2(ts) * w'^2(ts+1); and
!
!       Y = a_1(ts) * ( 2 * w'^3(ts) * w'^3(ts+1) ) / w'^2(ts);
!
!       where (ts) stands for the index of the current timestep and 
!       (ts+1) stands for the index of the timestep being advanced to.
!       The terms with timestep index (ts+1) are being solved for 
!       implicitly.
!
!       --------------------wp3p1-------------------------------- t(k+1)
!
!       =a3====wp2====a1=========wp3(interp)===================== m(k)
!
!       --------------------wp3------------------(d/dz) terms---- t(k)
!
!       =a3m1==wp2m1==a1m1=======wp3(interp)===================== m(k-1)
!
!       --------------------wp3m1-------------------------------- t(k-1)
!
!       The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1)
!       correspond with altitudes zt(k+1), zm(k), zt(k), zm(k-1), 
!       and zt(k-1), respectively.  The letter "t" is used for
!       thermodynamic levels and the letter "m" is used for momentum
!       levels.
!
!       dzt(k) = 1 / ( zm(k) - zm(k-1) )

!       References:
!-----------------------------------------------------------------------

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  kp1_tdiag = 1,   ! Thermodynamic superdiagonal index.
     .  k_mdiag   = 2,   ! Momentum superdiagonal index.
     .  k_tdiag   = 3,   ! Thermodynamic main diagonal index.
     .  km1_mdiag = 4,   ! Momentum subdiagonal index. 
     .  km1_tdiag = 5    ! Thermodynamic subdiagonal index.

        ! Input Variables
        real, intent(in) :: 
     .  wp3p1,   ! w'^3(k+1)                     [m^3/s^3]
     .  wp3,     ! w'^3(k)                       [m^3/s^3]
     .  wp3m1,   ! w'^3(k-1)                     [m^3/s^3]
     .  wp2,     ! w'^2(k)                       [m^2/s^2]
     .  wp2m1,   ! w'^2(k-1)                     [m^2/s^2]
     .  a1,      ! a1(k)                         [-]
     .  a1m1,    ! a1(k-1)                       [-]
     .  a3,      ! a3(k)                         [-]
     .  a3m1,    ! a3(k-1)                       [-]
     .  dzt,     ! Inverse of grid spacing (k)   [1/m]
     .  eps,     ! Model parameter               [-]
     .  wtol     ! Model parameter               [m/s]

        ! Return Variable
        real, dimension(5) :: lhs

        ! Brian tried a new discretization for the turbulent advection 
        ! term, which contains the term -d[ a_1 * (w'^3)^2 / w'^2 ] / dz.
        ! In order to help stabilize w'^3, a_1 has been pulled outside of
        ! the derivative.  On the left-hand side of the equation, this 
        ! effects the thermodynamic superdiagonal (kp1_tdiag), the 
        ! thermodynamic main diagonal (k_tdiag), and the thermodynamic 
        ! subdiagonal (km1_tdiag).

        ! Additionally, the discretization of the turbulent advection 
        ! term containing the term -d[ (a_3 + 3/2) * (w'^2)^2 ] / dz has 
        ! been altered to pull a3 outside of the derivative.  This was 
        ! done in order to help stabilize w'^3.  On the left-hand side of
        ! the equation, this effects the momentum superdiagonal (k_mdiag)
        ! and the momentum subdiagonal (km1_mdiag).

        ! Thermodynamic superdiagonal: [ x wp3(k+1,<ts+1>) ]
        lhs(kp1_tdiag)
!     .  = + (1.0/2.0) * dzt
!     .      * a1 * ( wp3p1 + wp3 ) / MAX( wp2, eps )
     .  = + (1.0/2.0) * ( a1 + a1m1 )
     .      * (1.0/2.0) * dzt
     .      * ( wp3p1 + wp3 ) / MAX( wp2, wtol**2 )

        ! Momentum superdiagonal: [ x wp2(k,<ts+1>) ]
        lhs(k_mdiag)
!     .  = + 2.0 * dzt * a3 * wp2
     .  = + 2.0 * dzt * 0.5*(a3+a3m1) * wp2

        ! Thermodynamic main diagonal: [ x wp3(k,<ts+1>) ]
        lhs(k_tdiag)
!     .  = + (1.0/2.0) * dzt
!     .      * (  a1 * ( wp3p1 + wp3 ) / MAX( wp2, eps )
!     .         - a1m1 * ( wp3 + wp3m1 ) / MAX( wp2m1, eps )
!     .        )
     .  = + (1.0/2.0) * ( a1 + a1m1 )
     .      * (1.0/2.0) * dzt
     .      * (  ( wp3p1 + wp3 ) / MAX( wp2, wtol**2 )
     .         - ( wp3 + wp3m1 ) / MAX( wp2m1, wtol**2 )
     .        )

        ! Momentum subdiagonal: [ x wp2(k-1,<ts+1>) ]
        lhs(km1_mdiag)
!     .  = - 2.0 * dzt * a3m1 * wp2m1
     .  = - 2.0 * dzt * 0.5*(a3+a3m1) * wp2m1

        ! Thermodynamic subdiagonal: [ x wp3(k-1,<ts+1>) ]
        lhs(km1_tdiag)
!     .  = - (1.0/2.0) * dzt
!     .      * a1m1 * ( wp3 + wp3m1 ) / MAX( wp2m1, eps )
     .  = - (1.0/2.0) * ( a1 + a1m1 )
     .      * (1.0/2.0) * dzt
     .      * ( wp3 + wp3m1 ) / MAX( wp2m1, wtol**2 )

        ! End of code that pulls out a3.
        ! End of Brian's a1 change.  Feb. 14, 2008.

        return
        end function wp3_terms_ta_tp_lhs

!========================================================================
        pure function wp3_terms_ac_pr2_lhs( C11_Skw_fnc, 
     .                                      wmm, wmmm1, dzt )
     .  result( lhs )

!       Description:
!       Accumulation of w'^3 and w'^3 pressure term 2:
!       implicit portion of the code.
!
!       The w'^3 accumulation term is completely implicit, while w'^3
!       pressure term 2 has both implicit and explicit components.  The
!       accumulation term and the implicit portion of pressure term 2 
!       are combined and solved together.
!
!       The values of w'^3 are found on thermodynamic levels, while the
!       values of wmm (mean vertical velocity on momentum levels) are 
!       found on momentum levels.  The vertical derivative of wmm is 
!       taken over the intermediate (central) thermodynamic level.  It
!       is then multiplied by w^'3 (implicitly calculated at timestep 
!       (ts+1)) and the coefficients to yield the desired results.
!
!       =======wmm=============================================== m(k)
!
!       ---------------d(wmm)/dz------------wp3------------------ t(k)
!
!       =======wmmm1============================================= m(k-1)
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
     .  k_tdiag   = 1    ! Thermodynamic main diagonal index.

        ! Input Variables
        real, intent(in) ::
     .  C11_Skw_fnc, ! C_11 skewness function (k)    [-]
     .  wmm,         ! wmm(k)                        [m/s]
     .  wmmm1,       ! wmm(k-1)                      [m/s]
     .  dzt          ! Inverse of grid spacing (k)   [1/m]

        ! Return Variable
        real :: lhs

        ! Thermodynamic main diagonal: [ x wp3(k,<ts+1>) ]
        lhs
     .  = + ( 1.0 - C11_Skw_fnc )
     .      * 3.0 * dzt * ( wmm - wmmm1 )

        return
        end function wp3_terms_ac_pr2_lhs

!========================================================================
        pure function wp3_term_pr1_lhs( C8, C8b, tauw3t, Skwt )
     .  result( lhs )

!       Description:
!       Pressure term 1 for w'^3: implicit portion of the code.
!
!       Pressure term 1 is the term: 
!       - (C_8/tau_w3t) * ( C_8b * Sk_wt^4 + 1 ) * w'^3;
!       which has to be linearized.  After it has been linearized, it 
!       is broken down into implicit (LHS) and explicit (RHS) 
!       components.
!
!       The values of w'^3 are found on the thermodynamic levels, as are
!       the values of tau_w3t, and Sk_wt; where 
!       Sk_wt = w^'3 / (w'^2)^(3/2); and w'^2 is interpolated to 
!       thermodynamic levels.  For this portion of the code, w'^3 is
!       solved for implicitly.

!       References:
!-----------------------------------------------------------------------

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  k_tdiag   = 1    ! Thermodynamic main diagonal index.

        ! Input Variables
        real, intent(in) ::
     .  C8,     ! Model parameter C_8                        [-]
     .  C8b,    ! Model parameter C_8b                       [-]
     .  tauw3t, ! Time scale at thermodynamic levels         [s]
     .  Skwt    ! Skewness of w at thermodynamic levels (k)  [-]

        ! Return Variable
        real :: lhs

        ! Thermodynamic main diagonal: [ x wp3(k,<ts+1>) ]
        lhs
     .  = + ( C8 / tauw3t ) * ( 5.0 * C8b * Skwt**4 + 1.0 )

        return
        end function wp3_term_pr1_lhs

!========================================================================
        pure function wp3_terms_ta_tp_rhs( wp3p1, wp3, wp3m1, 
     .                                     wp2, wp2m1, a1, a1m1,
     .                                     a3, a3m1, dzt, eps, wtol )
     .  result( rhs )

!       Description:
!       Turbulent advection and turbulent production of wp3:
!       explicit portion of the code.
!
!       The values of w'^3 are found on the thermodynamic levels, while
!       the values of w'^2, a_1, and a_3 are found on the momentum 
!       levels.  The variable w'^3 is interpolated to the intermediate 
!       momentum levels.  The values of the mathematical expressions
!       (called X and Y here) within the dX/dz and dY/dz terms are 
!       computed on the momentum levels.  Then, the derivatives (d/dz) 
!       of the expressions (X and Y) are taken over the central 
!       thermodynamic level, yielding the desired result.
!       In this function, the values of X and Y are as follows:
!
!       X = a_3(ts) * (w'^2(ts))^2; and
!
!       Y = a_1(ts) * (w'^3(ts))^2 / w'^2(ts);
!
!       where (ts) stands for the index of the current timestep.
!       The terms found in this function are being solved for 
!       explicitly.
!
!       --------------------wp3p1-------------------------------- t(k+1)
!
!       =a3====wp2====a1=========wp3(interp)===================== m(k)
!
!       --------------------wp3------------------(d/dz) terms---- t(k)
!
!       =a3m1==wp2m1==a1m1=======wp3(interp)===================== m(k-1)
!
!       --------------------wp3m1-------------------------------- t(k-1)
!
!       The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1)
!       correspond with altitudes zt(k+1), zm(k), zt(k), zm(k-1), 
!       and zt(k-1), respectively.  The letter "t" is used for
!       thermodynamic levels and the letter "m" is used for momentum
!       levels.
!
!       dzt(k) = 1 / ( zm(k) - zm(k-1) )

!       References:
!-----------------------------------------------------------------------

        implicit none

        ! Input Variables
        real, intent(in) :: 
     .  wp3p1,   ! w'^3(k+1)                     [m^3/s^3]
     .  wp3,     ! w'^3(k)                       [m^3/s^3]
     .  wp3m1,   ! w'^3(k-1)                     [m^3/s^3]
     .  wp2,     ! w'^2(k)                       [m^2/s^2]
     .  wp2m1,   ! w'^2(k-1)                     [m^2/s^2]
     .  a1,      ! a1(k)                         [-]
     .  a1m1,    ! a1(k-1)                       [-]
     .  a3,      ! a3(k)                         [-]
     .  a3m1,    ! a3(k-1)                       [-]
     .  dzt,     ! Inverse of grid spacing (k)   [1/m]
     .  eps,     ! Model parameter               [-]
     .  wtol     ! Model parameter               [m/s]

        ! Return Variable
        real :: rhs

        ! Brian tried a new discretization for the turbulent advection 
        ! term, which contains the term -d[ a_1 * (w'^3)^2 / w'^2 ] / dz.
        ! In order to help stabilize w'^3, a_1 has been pulled outside of
        ! the derivative.  This effects the right-hand side of the 
        ! equation, as well as the  left-hand side.

        ! Additionally, the discretization of the turbulent advection 
        ! term containing the term -d[ (a_3 + 3/2) * (w'^2)^2 ] / dz has 
        ! been altered to pull a3 outside of the derivative.  This was 
        ! done in order to help stabilize w'^3.  This effects the 
        ! right-hand side of the equation, as well as the left-hand side.

        rhs
!     .  = + dzt
!     .      * ( ( a3 * wp2**2 ) - ( a3m1 * wp2m1**2 ) )
!     .    + (1.0/4.0) * dzt
!     .      * (  a1 * ( wp3p1 + wp3 )**2
!     .                   / MAX( wp2, eps )
!     .         - a1m1 * ( wp3 + wp3m1 )**2
!     .                     / MAX( wp2m1, eps )
!     .        )
     .  = + dzt
     .      * (  ( 0.5*(a3+a3m1) * wp2**2 )
     .         - ( 0.5*(a3+a3m1) * wp2m1**2 )  )
     .    + (1.0/2.0) * ( a1 + a1m1 )
     .      * (1.0/4.0) * dzt
     .        * (  ( wp3p1 + wp3 )**2
     .                / MAX( wp2, wtol**2 )
     .           - ( wp3 + wp3m1 )**2
     .                / MAX( wp2m1, wtol**2 )
     .          )


        ! End of code that pulls out a3.
        ! End of Brian's a1 change.  Feb. 14, 2008.

        return
        end function wp3_terms_ta_tp_rhs

!========================================================================
        pure function wp3_terms_bp_pr2_rhs( C11_Skw_fnc, wp2thvp )
     .  result( rhs )

!       Description:
!       Buoyancy production of w'^3 and w'^3 pressure term 2:
!       explicit portion of the code.
!
!       The w'^3 buoyancy production term is completely explicit, while
!       w'^3 pressure term 2 has both implicit and explicit components. 
!       The buoyancy production term and the explicit portion of pressure
!       term 2 are combined and solved together.

!       References:
!-----------------------------------------------------------------------

        use constants, only: grav ! Gravitational acceleration [m/s^2]
        use parameters, only: T0  ! Reference temperature      [K]

        implicit none

        ! Input Variables
        real, intent(in) ::
     .  C11_Skw_fnc, ! C_11 skewness function (k)    [-]
     .  wp2thvp      ! w'^2 th_v' (k)                [K m^2/s^2]

        ! Return Variable
        real :: rhs

        rhs
     .  = + ( 1.0 - C11_Skw_fnc ) * 3.0 * ( grav/T0 ) * wp2thvp

        return
        end function wp3_terms_bp_pr2_rhs

!========================================================================
        pure function wp3_term_pr1_rhs( C8, C8b, tauw3t, Skwt, wp3 )
     .  result( rhs )

!       Description:
!       Pressure term 1 for w'^3: explicit portion of the code.
!
!       Pressure term 1 is the term: 
!       - (C_8/tau_w3t) * ( C_8b * Sk_wt^4 + 1 ) * w'^3;
!       which has to be linearized.  After it has been linearized, it 
!       is broken down into implicit (LHS) and explicit (RHS) 
!       components.
!
!       The values of w'^3 are found on the thermodynamic levels, as are
!       the values of tau_w3t, and Sk_wt; where 
!       Sk_wt = w^'3 / (w'^2)^(3/2); and w'^2 is interpolated to 
!       thermodynamic levels.  For this portion of the code, w'^3 is
!       solved for explicitly.

!       References:
!-----------------------------------------------------------------------

        implicit none

        ! Input Variables
        real, intent(in) ::
     .  C8,     ! Model parameter C_8                        [-]
     .  C8b,    ! Model parameter C_8b                       [-]
     .  tauw3t, ! Time scale at thermodynamic levels         [s]
     .  Skwt,   ! Skewness of w at thermodynamic levels (k)  [-]
     .  wp3     ! w'^3(k)                                    [m^3/s^3]

        ! Return Variable
        real :: rhs

        rhs
     .  = + ( C8 / tauw3t ) * ( 4.0 * C8b * Skwt**4 ) * wp3

        return
        end function wp3_term_pr1_rhs

!========================================================================

        end module wp23
