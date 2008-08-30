!------------------------------------------------------------------------
! $Id$
!===============================================================================
module wp23

  implicit none

  private ! Default Scope

  public :: timestep_wp23

  private :: wp23_solve, & 
             wp23_lhs, & 
             wp23_rhs, & 
             wp2_term_ta_lhs, & 
             wp2_terms_ac_pr2_lhs, & 
             wp2_term_dp1_lhs, & 
             wp2_term_pr1_lhs, & 
             wp2_terms_bp_pr2_rhs, & 
             wp2_term_dp1_rhs, &
             wp2_term_pr3_rhs, & 
             wp2_term_pr1_rhs, & 
             wp3_terms_ta_tp_lhs, & 
             wp3_terms_ac_pr2_lhs, & 
             wp3_term_pr1_lhs, & 
             wp3_terms_ta_tp_rhs, & 
             wp3_terms_bp_pr2_rhs, & 
             wp3_term_pr1_rhs

contains

  !=============================================================================
  subroutine timestep_wp23( dt, sigma_sqd_w, wm_zm, wm_zt, wpthvp, wp2thvp,  & 
                            um, vm, upwp, vpwp, up2, vp2, Kh_zm, Kh_zt, & 
                            tau_zm, tau_zt, Skw_zm, Skw_zt, a, & 
                            wp2_zt, wp2, wp3, err_code )

    ! Description:
    ! Advance w'^2 and w'^3 one timestep.

    ! References:
    ! Eqn. 12 & 18 on p. 3545--3546 of
    ! ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !   Method and Model Description'' Golaz, et al. (2002)
    ! JAS, Vol. 59, pp. 3540--3551.

    ! See also
    ! ``Equations for HOC'', Section 6:
    ! /Implict solution for the vertical velocity moments/
    !------------------------------------------------------------------------

    use grid_class, only:  & 
        gr,     & ! Variable(s)
        zt2zm,  & ! Procedure(s)
        zm2zt

    use parameters, only:  & 
        C11c,  & ! Variable(s)
        C11b,  & 
        C11,  & 
        C1c,  & 
        C1b,  & 
        C1,  & 
        c_K1,  & 
        c_ksqd,  & 
        c_K8

    use constants, only:  & 
        fstderr    ! Variable(s)

    use stats_precision, only:  & 
        time_precision ! Variable(s)

    use error_code, only:  & 
        lapack_error,  & ! Procedure(s)
        clubb_at_least_debug_level

    implicit none

    intrinsic :: exp

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      dt             ! Model timestep                           [s]

    real, intent(in), dimension(gr%nnzp) ::  & 
      sigma_sqd_w, & ! sigma_sqd_w on momentum levels           [-]
      wm_zm,       & ! w wind component on momentum levels      [m/s]
      wm_zt,       & ! w wind component on thermodynamic levels [m/s]
      wpthvp,      & ! w'th_v' (momentum levels)                [K m/s]
      wp2thvp,     & ! w'^2th_v' (thermodynamic levels)         [K m^2/s^2]
      um,          & ! u wind component (thermodynamic levels)  [m/s]
      vm,          & ! v wind component (thermodynamic levels)  [m/s]
      upwp,        & ! u'w' (momentum levels)                   [m^2/s^2]
      vpwp,        & ! v'w' (momentum levels)                   [m^2/s^2]
      up2,         & ! u'^2 (momentum levels)                   [m^2/s^2]
      vp2,         & ! v'^2 (momentum levels)                   [m^2/s^2]
      Kh_zm,       & ! Eddy diffusivity on momentum levels      [m^2/s]
      Kh_zt,       & ! Eddy diffusivity on thermodynamic levels [m^2/s]
      tau_zm,      & ! Time-scale tau on momentum levels        [s]
      tau_zt,      & ! Time-scale tau on thermodynamic levels   [s]
      Skw_zm,      & ! Skewness of w on momentum levels         [-]
      Skw_zt,      & ! Skewness of w on thermodynamic levels    [-]
      wp2_zt,      & ! w'^2 interpolated to thermodyamic levels [m^2/s^2]
      a              ! PDF parameter "a": pdf_parms(:,13)       [-]

    ! Input/Output
    real, dimension(gr%nnzp), intent(inout) ::  & 
      wp2,  & ! w'^2 (momentum levels)                          [m^2/s^2]
      wp3     ! w'^3 (thermodynamic levels)                     [m^3/s^3]

    integer, intent(inout) :: err_code ! Diagnostic

    ! Local Variables
    real, dimension(gr%nnzp) ::  & 
      tauw3t  ! Currently just tau_zt                           [s]

    ! Eddy Diffusion for w'^2 and w'^3.
    real, dimension(gr%nnzp) :: Kw1    ! w'^2 coef. eddy diff.  [m^2/s]
    real, dimension(gr%nnzp) :: Kw8    ! w'^3 coef. eddy diff.  [m^2/s]

    ! Note:  wp2_zt and wp2_zt_sqd_3pt, and wp3_zm and wp3_zm_sqd_3pt
    !        are used to help determine the coefficients of eddy
    !        diffusivity for wp2 and wp3, respectively.
    real, dimension(gr%nnzp) :: & 
      wp3_zm,          & ! w'^3 interpolated to momentum levels     [m^3/s^3]
      wp2_zt_sqd_3pt,  & ! (w'^2)^2; averaged over 3 points         [m^4/s^4]
      wp3_zm_sqd_3pt     ! (w'^3)^2; averaged over 3 points         [m^6/s^6]

    ! Internal variables for C11 function, Vince Larson 13 Mar 2005
    ! Brian added C1 function.
    real, dimension(gr%nnzp) ::  & 
      C1_Skw_fnc,  & ! C_1 parameter with Sk_w applied              [-]
      C11_Skw_fnc    ! C_11 parameter with Sk_w applied             [-]
    ! End Vince Larson's addition.

    integer :: k, km1, kp1  ! Array indices


    !-----------------------------------------------------------------------



!       Define tauw

!        tauw3t = tau_zt
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
!          tauw3t(k) = tau_zt(k) / ( 0.005*Skw**4 + 1.0 )
!
!        end do

    tauw3t = tau_zt

    ! Vince Larson added code to make C11 function of Skw. 13 Mar 2005
    ! If this code is used, C11 is no longer relevant, i.e. constants
    !    are hardwired.

    ! Calculate C_1 and C_11 as functions of skewness of w.

    C11_Skw_fnc(1:gr%nnzp) =  & 
    C11b + (C11-C11b)*EXP( -(1.0/2.0) * (Skw_zt(1:gr%nnzp)/C11c)**2 )

    C1_Skw_fnc(1:gr%nnzp) =  & 
    C1b + (C1-C1b)*EXP( -(1.0/2.0) * (Skw_zm(1:gr%nnzp)/C1c)**2 )

    !C11_Skw_fnc = C11
    !C1_Skw_fnc = C1

    ! Vince Larson added extra diffusion based on wp2.  21 Dec 2007.
    ! Vince Larson added extra diffusion based on wp3.  15 Dec 2007.

    ! Interpolate w'^3 from thermodynamic levels to momentum levels.
    ! This will be used for the w'^3 turbulent advection (ta) and
    ! turbulent production (tp) combined term.
    ! This is also used for extra diffusion based on a three-point
    ! average of (w'^3)^2.
    wp3_zm = zt2zm( wp3 )

    ! Interpolate w'^2 from momentum levels to thermodynamic levels.
    ! This is used for extra diffusion based on a three-point
    ! average of (w'^2)^2.
    !wp2_zt = max( zm2zt( wp2 ), 0.0 )   ! Positive definite quantity

    do k = 1, gr%nnzp, 1

      km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nnzp )

      ! Compute the square of wp2_zt, averaged over 3 points.  15 Dec 2007
      wp2_zt_sqd_3pt(k)  & 
         = ( wp2_zt(km1)**2 + wp2_zt(k)**2 + wp2_zt(kp1)**2 ) / 3.0

      ! Compute the square of wp3_zm, averaged over 3 points.  15 Dec 2007
      wp3_zm_sqd_3pt(k)  & 
         = ( wp3_zm(km1)**2 + wp3_zm(k)**2 + wp3_zm(kp1)**2 ) / 3.0

    enddo

    ! End Vince Larson's addition.

    ! Define the Coefficent of Eddy Diffusivity for the wp2 and wp3.
    do k = 1, gr%nnzp, 1

      ! Kw1 is used for wp2, which is located on momentum levels.
      ! Kw1 is located on thermodynamic levels.
      ! Kw1 = c_K1 * Kh_zt
      Kw1(k) = c_K1 * Kh_zt(k)
      ! Vince Larson added extra diffusion based on wp2.  21 Dec 2007.
      ! Kw1 must have units of m^2/s.  Since wp2_zt_sqd_3pt has units
      ! of m^4/s^4, c_Ksqd is given units of s^3/m^2 in this case.
      Kw1(k) = Kw1(k) + c_Ksqd * wp2_zt_sqd_3pt(k)
      ! End Vince Larson's addition.

      ! Kw8 is used for wp3, which is located on thermodynamic levels.
      ! Kw8 is located on momentum levels.
      ! Note: Kw8 is defined to be 1/2 of Kh_zm.
      ! Kw8 = c_K8 * Kh_zm
      Kw8(k) = c_K8 * Kh_zm(k)
      ! Vince Larson added extra diffusion based on wp3.  15 Dec 2007.
      ! Kw8 must have units of m^2/s.  Since wp3_zm_sqd_3pt has units
      ! of m^6/s^6, c_Ksqd is given units of s^5/m^4 in this case.
      Kw8(k) = Kw8(k) + c_Ksqd * wp3_zm_sqd_3pt(k)
      ! End Vince Larson's addition.

    enddo

    ! Solve semi-implicitly
    call wp23_solve( dt, sigma_sqd_w, wm_zm, wm_zt, wpthvp, wp2thvp, & ! Intent(in)
                     um, vm, upwp, vpwp, up2, vp2, Kw1, &              ! Intent(in)
                     Kw8, Skw_zt, tau_zm, tauw3t, C1_Skw_fnc, &        ! Intent(in)
                     C11_Skw_fnc, wp3_zm, &                            ! Intent(in)
                     wp2, wp3, err_code )                              ! Intent(inout)

!       Error output
!       Joshua Fasching Feb 2008
    if ( lapack_error( err_code ) .and.  & 
         clubb_at_least_debug_level( 1 ) ) then

      write(fstderr,*) "Errors in timestep_wp32"

      write(fstderr,*) "Intent(in)"

      write(fstderr,*) "dt = ", dt
      write(fstderr,*) "sigma_sqd_w = ", sigma_sqd_w
      write(fstderr,*) "wm_zm = ", wm_zm
      write(fstderr,*) "wm_zt = ", wm_zt
      write(fstderr,*) "wpthvp = ", wpthvp
      write(fstderr,*) "wp2thvp = ", wp2thvp
      write(fstderr,*) "um = ", um
      write(fstderr,*) "vm = ", vm
      write(fstderr,*) "upwp = ", upwp
      write(fstderr,*) "vpwp = ", vpwp
      write(fstderr,*) "up2 = ", up2
      write(fstderr,*) "vp2 = ", vp2
      write(fstderr,*) "Kh_zm = ", Kh_zm
      write(fstderr,*) "Kh_zt = ", Kh_zt
      write(fstderr,*) "tau_zm = ", tau_zm
      write(fstderr,*) "tau_zt = ", tau_zt
      write(fstderr,*) "Skw_zm = ", Skw_zm
      write(fstderr,*) "Skw_zt = ", Skw_zt
      write(fstderr,*) "a = ", a
      write(fstderr,*) "wp2zt = ", wp2_zt

      write(fstderr,*) "Intent(in/out)"

      write(fstderr,*) "wp2 = ", wp2
      write(fstderr,*) "wp3 = ", wp3

    endif

    return

  end subroutine timestep_wp23

  !=============================================================================
  subroutine wp23_solve( dt, sigma_sqd_w, wm_zm, wm_zt, wpthvp, wp2thvp, & 
                         um, vm, upwp, vpwp, up2, vp2, Kw1, & 
                         Kw8, Skw_zt, tau1m, tauw3t, C1_Skw_fnc, & 
                         C11_Skw_fnc, wp3_zm, &
                         wp2, wp3, err_code )

    ! Description:
    ! Decompose, and back substitute the matrix for wp2/wp3

    ! References:
    ! _Equations for HOC_ section 6.3
    !------------------------------------------------------------------------

    use grid_class, only:  & 
        gr,         & ! Variable(s) 
        zm2zt,      & ! Procedure(s)
        ddzt       ! Function

    use constants, only: & 
        emin,       & ! Variables(s)
        eps

    use model_flags, only:  & 
        l_tke_aniso,  & ! Variable(s)
        l_hole_fill

    use stats_precision, only:  & 
        time_precision  ! Variable(s)

    use lapack_wrap, only:  & 
        band_solve,  & ! Procedure(s) 
        band_solvex

    use fill_holes, only: & 
        fill_holes_driver

    use error_code, only:  & 
        lapack_error ! Procedure(s)

    use clip_explicit, only: &
        clip_variance, & ! Procedure(s)
        clip_skewness

    use stats_type, only: & 
        stat_begin_update,  & ! Procedure(s)
        stat_update_var_pt, &
        stat_end_update,  &
        stat_end_update_pt

    use stats_variables, only:  & 
        zm,         & ! Variable(s)
        zt, & 
        sfc, & 
        l_stats_samp, & 
        iwp2_bt, & 
        iwp2_ta, & 
        iwp2_ma, & 
        iwp2_pd, & 
        iwp2_ac, & 
        iwp2_dp1, & 
        iwp2_dp2, & 
        iwp2_pr1, & 
        iwp2_pr2, & 
        iwp3_bt, & 
        iwp3_ta, & 
        iwp3_ma, & 
        iwp3_tp, & 
        iwp3_ac, & 
        iwp3_dp1, & 
        iwp3_pr1, & 
        iwp3_pr2, & 
        iwp23_cn, & 
        zmscr01, & 
        zmscr02, & 
        zmscr03, & 
        zmscr04, & 
        zmscr05, & 
        zmscr06, & 
        zmscr07, & 
        zmscr08, & 
        zmscr09, & 
        zmscr10, & 
        zmscr11, & 
        zmscr12, & 
        ztscr01, & 
        ztscr02, & 
        ztscr03, & 
        ztscr04, & 
        ztscr05, & 
        ztscr06, & 
        ztscr07, & 
        ztscr08, & 
        ztscr09, & 
        ztscr10, & 
        ztscr11, & 
        ztscr12, & 
        ztscr13, & 
        ztscr14, & 
        ztscr15, & 
        ztscr16


    implicit none

! External
    intrinsic :: max, min, sqrt

! Parameter Constants
    integer, parameter :: & 
      nsub = 2,   & ! Number of subdiagonals in the LHS matrix
      nsup = 2,   & ! Number of superdiagonals in the LHS matrix
      nrhs = 1      ! Number of RHS vectors

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      dt              ! Timestep                                  [s]

    real, intent(in), dimension(gr%nnzp) ::  & 
      sigma_sqd_w,  & ! sigma_sqd_w on momentum levels            [-]
      wm_zm,        & ! w wind component on momentum levels       [m/s]
      wm_zt,        & ! w wind component on thermodynamic levels  [m/s]
      wpthvp,       & ! w'th_v' (momentum levels)                 [K m/s]
      wp2thvp,      & ! w'^2th_v' (thermodynamic levels)          [K m^2/s^2]
      um,           & ! u wind component (thermodynamic levels)   [m/s]
      vm,           & ! v wind component (thermodynamic levels)   [m/s]
      upwp,         & ! u'w' (momentum levels)                    [m^2/s^2]
      vpwp,         & ! v'w' (momentum levels)                    [m^2/s^2]
      up2,          & ! u'^2 (momentum levels)                    [m^2/s^2]
      vp2,          & ! v'^2 (momentum levels)                    [m^2/s^2]
      Kw1,          & ! Coefficient of eddy diffusivity for w'^2  [m^2/s]
      Kw8,          & ! Coefficient of eddy diffusivity for w'^3  [m^2/s]
      Skw_zt,       & ! Skewness of w on thermodynamic levels     [-]
      tau1m,        & ! Time-scale tau on momentum levels         [s]
      tauw3t,       & ! Time-scale tau on thermodynamic levels    [s]
      C1_Skw_fnc,   & ! C_1 parameter with Sk_w applied           [-]
      C11_Skw_fnc,  & ! C_11 parameter with Sk_w applied          [-]
      wp3_zm          ! w'^3 interpolated to momentum levels      [m^3/s^3]

    ! Input/Output Variables
    real, dimension(gr%nnzp), intent(inout) ::  & 
      wp2,  & ! w'^2 (momentum levels)                            [m^2/s^2]
      wp3     ! w'^3 (thermodynamic levels)                       [m^3/s^3]

    integer, intent(inout) :: err_code ! Have any errors occured?

    ! Local Variables
    real, dimension(nsup+nsub+1,2*gr%nnzp) ::  & 
      lhs ! Implicit contributions to wp2/wp3 (band diag. matrix)

    real, dimension(2*gr%nnzp) ::  & 
      rhs   ! RHS of band matrix

!        real, target, dimension(2*gr%nnzp) ::
    real, dimension(2*gr%nnzp) ::  & 
      solut ! Solution to band diagonal system.

    real, dimension(gr%nnzp) ::  & 
      a1,  & ! a_1 (momentum levels); See eqn. 24 in `Equations for HOC' [-]
      a3     ! a_3 (momentum levels); See eqn. 26 in `Equations for HOC' [-]

    real, dimension(gr%nnzp) ::  & 
      a1_zt,  & ! a_1 interpolated to thermodynamic levels        [-]
      a3_zt,  & ! a_3 interpolated to thermodynamic levels        [-]
      wp2_zt    ! w'^2 interpolated to thermodyamic levels        [m^2/s^2]

!      real, dimension(gr%nnzp) ::  &
!        wp2_n ! w'^2 at the previous timestep           [m^2/s^2]

    real ::  & 
      rcond  ! Est. of the reciprocal of the condition #

    ! Array indices
    integer :: k, km1, kp1, k_wp2, k_wp3

    ! Set logical to true for Crank-Nicholson diffusion scheme
    ! or to false for completely implicit diffusion scheme.
    ! Note:  Although Crank-Nicholson diffusion has usually been used for wp2 and
    !        wp3 in the past, we found that using completely implicit diffusion
    !        stabilized the deep convective cases more while having almost no effect
    !        on the boundary layer cases.  Brian; 1/4/2008.
!    logical, parameter :: l_crank_nich_diff = .true.
    logical, parameter :: l_crank_nich_diff = .false.

    if (l_stats_samp) then
      call stat_begin_update( iwp2_bt, real(wp2 / dt), zm )

      call stat_begin_update( iwp3_bt, real(wp3 / dt), zt )
    endif


    ! Define a_1 and a_3 (both are located on momentum levels).
    ! They are variables that are both functions of sigma_sqd_w (where
    ! sigma_sqd_w is located on momentum levels).
    ! Note: some compilers appear to interpret the pow function with
    ! a positive integer exponent differently than a repeated
    ! multiply. -dschanen 19 March 2007

    a1 = 1.0 / ( 1.0 - sigma_sqd_w )

    a3 = 3.0 * sigma_sqd_w*sigma_sqd_w & 
         + 6.0*(1.0-sigma_sqd_w)*sigma_sqd_w  & 
         + (1.0-sigma_sqd_w)*(1.0-sigma_sqd_w) & 
         - (3.0/2.0)

    ! Interpolate a_1 and a_3 from momentum levels to thermodynamic
    ! levels.  This will be used for the w'^3 turbulent advection
    ! (ta) and turbulent production (tp) combined term.
    a1_zt  = max( zm2zt( a1 ), 0.0 )   ! Positive definite quantity
    a3_zt  = zm2zt( a3 )

    ! Compute the implicit portion of the w'^2 and w'^3 equations.
    ! Build the left-hand side matrix.
    !call wp23_lhs( dt, wp2, wp3_zm, wm_zm, wm_zt, a1, a3, a1_zt,  &
    call wp23_lhs( dt, wp2, wp3_zm, wm_zm, wm_zt, a1_zt,  &
                   a3_zt, Kw1, Kw8, Skw_zt, tau1m, tauw3t,  & 
                   C1_Skw_fnc, C11_Skw_fnc, l_crank_nich_diff,  & 
                   lhs )

    ! Compute the explicit portion of the w'^2 and w'^3 equations.
    ! Build the right-hand side vector.
    !call wp23_rhs( dt, wp2, wp3, wp3_zm, a1_zt, a1, a3, &
    call wp23_rhs( dt, wp2, wp3, wp3_zm, a1_zt, &
                   a3_zt, wpthvp, wp2thvp, um, vm,  & 
                   upwp, vpwp, up2, vp2, Kw1, Kw8,  & 
                   Skw_zt, tau1m, tauw3t, C1_Skw_fnc, &
                   C11_Skw_fnc, l_crank_nich_diff, rhs )

    ! Solve the system of equations for w'^2 and w'^3.
    if ( l_stats_samp .and. iwp23_cn > 0 ) then

      ! Perform LU decomp and solve system (LAPACK with diagnostics)
      call band_solvex( "wp23", nsup, nsub, 2*gr%nnzp, nrhs, & 
                        lhs, rhs, solut, rcond, err_code )

      ! Est. of the condition number of the w'^2/w^3 LHS matrix
      call stat_update_var_pt( iwp23_cn, 1, 1.0 / rcond, sfc )

    else
      ! Perform LU decomp and solve system (LAPACK)
      call band_solve( "wp23", nsup, nsub, 2*gr%nnzp, nrhs, & 
                       lhs, rhs, solut, err_code )
    end if


    if ( lapack_error( err_code ) ) return


    ! Copy result into output arrays and clip

    do k = 1, gr%nnzp

      km1 = max( k-1, 1)
      kp1 = min( k+1, gr%nnzp )

      k_wp3 = 2*k - 1
      k_wp2 = 2*k

      ! wp2_n(k) = wp2(k) ! For the positive definite scheme

      wp2(k) = solut(k_wp2)
      wp3(k) = solut(k_wp3)

    end do

    if (l_stats_samp) then

    ! Finalize implicit contributions for wp2

      do k = 2, gr%nnzp-1

        km1 = max( k-1, 1 )
        kp1 = min( k+1, gr%nnzp )

        call stat_end_update_pt( iwp2_dp1, k, & 
           zmscr01(k) * wp2(k), zm )

        call stat_update_var_pt( iwp2_dp2, k, & 
           zmscr02(k) * wp2(km1) & 
         + zmscr03(k) * wp2(k) & 
         + zmscr04(k) * wp2(kp1), zm )

        call stat_update_var_pt( iwp2_ta, k, & 
           zmscr05(k) * wp3(k) & 
         + zmscr06(k) * wp3(kp1), zm )

        call stat_update_var_pt( iwp2_ma, k, & 
           zmscr07(k) * wp2(km1) & 
         + zmscr08(k) * wp2(k) & 
         + zmscr09(k) * wp2(kp1), zm )

        call stat_update_var_pt( iwp2_ac, k,  & 
           zmscr10(k) * wp2(k), zm )

        if ( l_tke_aniso ) then
          call stat_end_update_pt( iwp2_pr1, k, & 
             zmscr12(k) * wp2(k), zm )
        endif

        call stat_end_update_pt( iwp2_pr2, k, & 
           zmscr11(k) * wp2(k), zm )

      enddo

      ! Finalize implicit contributions for wp3

      do k = 2, gr%nnzp-1, 1

        km1 = max( k-1, 1 )
        kp1 = min( k+1, gr%nnzp )

        call stat_end_update_pt( iwp3_pr1, k, & 
           ztscr01(k) * wp3(k), zt )

        call stat_update_var_pt( iwp3_dp1, k, & 
           ztscr02(k) * wp3(km1) & 
         + ztscr03(k) * wp3(k) & 
         + ztscr04(k) * wp3(kp1), zt )

        call stat_end_update_pt( iwp3_ta, k, & 
           ztscr05(k) * wp3(km1) & 
         + ztscr06(k) * wp2(km1) & 
         + ztscr07(k) * wp3(k) & 
         + ztscr08(k) * wp2(k) & 
         + ztscr09(k) * wp3(kp1), zt )

        call stat_end_update_pt( iwp3_tp, k,  & 
           ztscr10(k) * wp2(km1) & 
         + ztscr11(k) * wp2(k), zt )

        call stat_update_var_pt( iwp3_ma, k, & 
           ztscr12(k) * wp3(km1) & 
         + ztscr13(k) * wp3(k) & 
         + ztscr14(k) * wp3(kp1), zt )

        call stat_update_var_pt( iwp3_ac, k, & 
           ztscr15(k) * wp3(k), zt )

        call stat_end_update_pt( iwp3_pr2, k, & 
           ztscr16(k) * wp3(k), zt )

      enddo
    endif ! l_stats_samp


    if ( l_stats_samp ) then
      ! Store previous value for effect of the positive definite scheme
      call stat_begin_update( iwp2_pd, real( wp2 / dt ), zm )
    endif

    if ( l_hole_fill .and. any( wp2 < 2./3*emin ) ) then

      ! Use a simple hole filling algorithm
      call fill_holes_driver( 2, 2./3.*emin, "zm", wp2 )

    endif ! wp2

    if ( l_stats_samp ) then
      ! Store updated value for effect of the positive definite scheme
      call stat_end_update( iwp2_pd, real( wp2 / dt ), zm )
    endif


    ! Clip w'^2 at a minimum threshold.
    call clip_variance( "wp2", dt, 2./3.*emin, wp2 )

    ! Interpolate w'^2 from momentum levels to thermodynamic levels.
    ! This is used for the clipping of w'^3 according to the value
    ! of Sk_w now that w'^2 and w'^3 have been advanced one timestep.
    wp2_zt = max( zm2zt( wp2 ), 2./3.*emin )   ! Positive definite quantity

    ! Clip w'^3 by limiting skewness.
    call clip_skewness( dt, wp2_zt, wp3 )


    if (l_stats_samp) then
      call stat_end_update( iwp2_bt, real( wp2 / dt ), zm )

      call stat_end_update( iwp3_bt, real( wp3 / dt ), zt )
    endif


    return
  end subroutine wp23_solve

  !=============================================================================
  !subroutine wp23_lhs( dt, wp2, wp3_zm, wm_zm, wm_zt, a1, a3, a1_zt,  &
  subroutine wp23_lhs( dt, wp2, wp3_zm, wm_zm, wm_zt, a1_zt,  & 
                       a3_zt, Kw1, Kw8, Skw_zt, tau1m, tauw3t,  & 
                       C1_Skw_fnc, C11_Skw_fnc, l_crank_nich_diff,  & 
                       lhs )

    ! Description:
    ! Compute LHS band diagonal matrix for w'^2 and w'^3.
    ! This subroutine computes the implicit portion 
    ! of the w'^2 and w'^3 equations.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr ! Variable

    use parameters, only:  & 
        C4,  & ! Variables
        C5,  & 
        C8,  & 
        C8b, & 
        C12, & 
        nu1, & 
        nu8

    use constants, only:  & 
        wtol,  & ! Variables
        eps

    use model_flags, only: & 
        l_tke_aniso ! Variables

    use diffusion, only: & 
        diffusion_zm_lhs,  & ! Procedures
        diffusion_zt_lhs

    use mean_adv, only: & 
        term_ma_zm_lhs,  & ! Procedures
        term_ma_zt_lhs

    use stats_precision, only: time_precision


    use stats_variables, only:       & 
        zmscr01,    & 
        zmscr02,    & 
        zmscr03,    & 
        zmscr04,    & 
        zmscr05,    & 
        zmscr06,    & 
        zmscr07,    & 
        zmscr08,    & 
        zmscr09,    & 
        zmscr11,    & 
        zmscr10,    & 
        zmscr12,    & 
        ztscr01,    & 
        ztscr02,    & 
        ztscr03,    & 
        ztscr04,    & 
        ztscr05,    & 
        ztscr06,    & 
        ztscr07,    & 
        ztscr08,    & 
        ztscr09,    & 
        ztscr10,    & 
        ztscr11,    & 
        ztscr12,    & 
        ztscr13,    & 
        ztscr14,    & 
        ztscr15,    & 
        ztscr16,    & 
        l_stats_samp, & 
        iwp2_dp1, & 
        iwp2_dp2, & 
        iwp2_ta, & 
        iwp2_ma, & 
        iwp2_ac, & 
        iwp2_pr2, & 
        iwp2_pr1, & 
        iwp3_ta, & 
        iwp3_tp, & 
        iwp3_ma, & 
        iwp3_ac, & 
        iwp3_pr2, & 
        iwp3_pr1, & 
        iwp3_dp1


    implicit none

    ! Parameter Constants
    integer, parameter :: & 
      nsub = 2,   & ! Number of subdiagonals in the LHS matrix.
      nsup = 2      ! Number of superdiagonals in the LHS matrix.

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      dt             ! Timestep length                                [s]

    real, dimension(gr%nnzp), intent(in) ::  & 
      wp2,         & ! w'^2 (momentum levels)                         [m^2/s^2]
      wp3_zm,      & ! w'^3 interpolated to momentum levels           [m^3/s^3]
      wm_zm,       & ! w wind component on momentum levels            [m/s]
      wm_zt,       & ! w wind component on thermodynamic levels       [m/s]
!      a1,          & ! sigma_sqd_w-related term a_1 (momentum levels) [-]
!      a3,          & ! sigma_sqd_w-related term a_1 (momentum levels) [-]
      a1_zt,       & ! a_1 interpolated to thermodynamic levels       [-]
      a3_zt,       & ! a_3 interpolated to thermodynamic levels       [-]
      Kw1,         & ! Coefficient of eddy diffusivity for w'^2       [m^2/s]
      Kw8,         & ! Coefficient of eddy diffusivity for w'^3       [m^2/s]
      Skw_zt,      & ! Skewness of w on thermodynamic levels          [-]
      tau1m,       & ! Time-scale tau on momentum levels              [s]
      tauw3t,      & ! Time-scale tau on thermodynamic levels         [s]
      C1_Skw_fnc,  & ! C_1 parameter with Sk_w applied                [-]
      C11_Skw_fnc    ! C_11 parameter with Sk_w applied               [-]

    logical, intent(in) :: & 
      l_crank_nich_diff   ! Turns on/off Crank-Nicholson diffusion.

    ! Output Variable
    real, dimension(nsup+nsub+1,2*gr%nnzp), intent(out) ::  & 
      lhs ! Implicit contributions to wp2/wp3 (band diag. matrix)

    ! Local Variables

    ! Array indices
    integer :: k, km1, kp1, k_wp2, k_wp3

    real, dimension(5) :: tmp


    ! Initialize the left-hand side matrix to 0.
    lhs = 0.0

    do k = 2, gr%nnzp-1, 1

      ! Define indices

      km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nnzp )

      k_wp3 = 2*k - 1
      k_wp2 = 2*k


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
      lhs(3,k_wp2) & 
      = real( + 1.0 / dt )

      ! LHS mean advection (ma) term.
      lhs((/3-2,3,3+2/),k_wp2) & 
      = lhs((/3-2,3,3+2/),k_wp2) & 
      + term_ma_zm_lhs( wm_zm(k), gr%dzm(k), k )

      ! LHS turbulent advection (ta) term.
      lhs((/3-1,3+1/),k_wp2) & 
      = lhs((/3-1,3+1/),k_wp2) & 
      + wp2_term_ta_lhs( gr%dzm(k) )

      ! LHS accumulation (ac) term and pressure term 2 (pr2).
      lhs(3,k_wp2) & 
      = lhs(3,k_wp2) & 
      + wp2_terms_ac_pr2_lhs( C5, wm_zt(kp1), wm_zt(k), gr%dzm(k)  )

      ! LHS dissipation term 1 (dp1).
      lhs(3,k_wp2) & 
      = lhs(3,k_wp2) & 
      + wp2_term_dp1_lhs( C1_Skw_fnc(k), tau1m(k) )

      ! LHS eddy diffusion term: dissipation term 2 (dp2).
      if ( l_crank_nich_diff ) then
        ! Eddy diffusion for wp2 using a Crank-Nicholson time step.
        lhs((/3-2,3,3+2/),k_wp2) & 
        = lhs((/3-2,3,3+2/),k_wp2) & 
        + (1.0/2.0) & 
        * diffusion_zm_lhs( Kw1(k), Kw1(kp1), nu1, & 
                            gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )
      else
        ! Eddy diffusion for wp2 using a completely implicit time step.
        lhs((/3-2,3,3+2/),k_wp2) & 
        = lhs((/3-2,3,3+2/),k_wp2) & 
        + diffusion_zm_lhs( Kw1(k), Kw1(kp1), nu1, & 
                            gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )
      endif

      ! LHS pressure term 1 (pr1).
      if ( l_tke_aniso ) then
        ! Add in this term if we're not assuming tke = 1.5 * wp2
        lhs(3,k_wp2) & 
        = lhs(3,k_wp2) & 
        + wp2_term_pr1_lhs( C4, tau1m(k) )
      endif

      if ( l_stats_samp ) then

        ! Statistics: implicit contributions for wp2.

        if ( iwp2_dp1 > 0 ) then
          zmscr01(k) = & 
          - wp2_term_dp1_lhs( C1_Skw_fnc(k), tau1m(k) )
        endif

        if ( iwp2_dp2 > 0 ) then
          if ( l_crank_nich_diff ) then
            ! Eddy diffusion for wp2 using a Crank-Nicholson time step.
            tmp(1:3) & 
            = (1.0/2.0) & 
            * diffusion_zm_lhs( Kw1(k), Kw1(kp1), nu1, & 
                              gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )
          else
            ! Eddy diffusion for wp2 using a completely implicit time step.
            tmp(1:3) & 
            = diffusion_zm_lhs( Kw1(k), Kw1(kp1), nu1, & 
                                gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )
          endif

          zmscr02(k) = - tmp(3)
          zmscr03(k) = - tmp(2)
          zmscr04(k) = - tmp(1)

        endif

        if ( iwp2_ta > 0 ) then
          tmp(1:2) =  & 
          + wp2_term_ta_lhs( gr%dzm(k) )
          zmscr05(k) = - tmp(2)
          zmscr06(k) = - tmp(1)
        endif

        if ( iwp2_ma > 0 ) then
          tmp(1:3) = & 
          + term_ma_zm_lhs( wm_zm(k), gr%dzm(k), k )
          zmscr07(k) = - tmp(3)
          zmscr08(k) = - tmp(2)
          zmscr09(k) = - tmp(1)
        endif

        if ( iwp2_ac > 0 ) then
          zmscr10(k) =  & 
          - wp2_terms_ac_pr2_lhs( 0.0, wm_zt(kp1), wm_zt(k), gr%dzm(k)  )
        endif

        if ( iwp2_pr2 > 0 ) then
          zmscr11(k) =  & 
          - wp2_terms_ac_pr2_lhs( (1.0+C5), wm_zt(kp1), wm_zt(k),  & 
                                  gr%dzm(k)  )
        endif

        if ( iwp2_pr1 > 0 .and. l_tke_aniso ) then
          zmscr12(k) = - wp2_term_pr1_lhs( C4, tau1m(k) )
        endif

      endif



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
      lhs(3,k_wp3) & 
      = real( + 1.0 / dt )

      ! LHS mean advection (ma) term.
      lhs((/3-2,3,3+2/),k_wp3) & 
      = lhs((/3-2,3,3+2/),k_wp3) & 
      + term_ma_zt_lhs( wm_zt(k), gr%dzt(k), k )

      ! LHS turbulent advection (ta) and turbulent production (tp) terms.
      lhs(3-2:3+2,k_wp3) & 
      = lhs(3-2:3+2,k_wp3) & 
      + wp3_terms_ta_tp_lhs( wp3_zm(k), wp3_zm(km1),  & 
                             wp2(k), wp2(km1),  &
                         !a1(k), a1(km1), a3(k), a3(km1),  &
                             a1_zt(k), a3_zt(k),  & 
                             gr%dzt(k), wtol, k )

      ! LHS accumulation (ac) term and pressure term 2 (pr2).
      lhs(3,k_wp3) & 
      = lhs(3,k_wp3) & 
      + wp3_terms_ac_pr2_lhs( C11_Skw_fnc(k), & 
                              wm_zm(k), wm_zm(km1), gr%dzt(k) )

      ! LHS pressure term 1 (pr1).
      lhs(3,k_wp3) & 
      = lhs(3,k_wp3) & 
      + wp3_term_pr1_lhs( C8, C8b, tauw3t(k), Skw_zt(k) )

      ! LHS eddy diffusion term: dissipation term 1 (dp1).
      !  Added a new constant, C12.
      !  Initially, this new constant will be set to 1.0 -dschanen 9/19/05
      if ( l_crank_nich_diff ) then
        ! Eddy diffusion for wp3 using a Crank-Nicholson time step.
        lhs((/3-2,3,3+2/),k_wp3) & 
        = lhs((/3-2,3,3+2/),k_wp3) & 
        + C12 * (1.0/2.0) & 
        * diffusion_zt_lhs( Kw8(k), Kw8(km1), nu8, & 
                            gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
      else
        ! Eddy diffusion for wp3 using a completely implicit time step.
        lhs((/3-2,3,3+2/),k_wp3) & 
        = lhs((/3-2,3,3+2/),k_wp3) & 
        + C12  & 
        * diffusion_zt_lhs( Kw8(k), Kw8(km1), nu8, & 
                            gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
      endif

      if (l_stats_samp) then

        ! Statistics: implicit contributions for wp3.

        if ( iwp3_ta > 0 ) then
          tmp(1:5) =  & 
          wp3_terms_ta_tp_lhs( wp3_zm(k), wp3_zm(km1),  & 
                               wp2(k), wp2(km1),  &
                           !a1(k), a1(km1),  &
                           !a3(k)+(3.0/2.0), a3(km1)+(3.0/2.0),  &
                               a1_zt(k), a3_zt(k)+(3.0/2.0),  &
                               gr%dzt(k), wtol, k )
          ztscr05(k) = -tmp(5)
          ztscr06(k) = -tmp(4)
          ztscr07(k) = -tmp(3)
          ztscr08(k) = -tmp(2)
          ztscr09(k) = -tmp(1)
        endif

        if ( iwp3_tp > 0 ) then
          tmp(1:5) =  & 
          wp3_terms_ta_tp_lhs( wp3_zm(k), wp3_zm(km1),  & 
                               wp2(k), wp2(km1),  & 
                           !0.0, 0.0,  &
                           !0.0-(3.0/2.0), 0.0-(3.0/2.0),  &
                               0.0, 0.0-(3.0/2.0),  & 
                               gr%dzt(k), wtol, k )
          ztscr10(k) = -tmp(4)
          ztscr11(k) = -tmp(2)
        endif

        if ( iwp3_ma > 0 ) then
          tmp(1:3) = & 
          term_ma_zt_lhs( wm_zt(k), gr%dzt(k), k )
          ztscr12(k) = -tmp(3)
          ztscr13(k) = -tmp(2)
          ztscr14(k) = -tmp(1)
        endif

        if ( iwp3_ac > 0 ) then
          ztscr15(k) =  & 
          - wp3_terms_ac_pr2_lhs( 0.0, & 
                                  wm_zm(k), wm_zm(km1), gr%dzt(k) )
        endif

        if ( iwp3_pr2 > 0 ) then
          ztscr16(k) = & 
          - wp3_terms_ac_pr2_lhs( (1.0+C11_Skw_fnc(k)), & 
                                  wm_zm(k), wm_zm(km1), gr%dzt(k) )
        endif

        if ( iwp3_pr1 > 0 ) then
          ztscr01(k) = & 
          - wp3_term_pr1_lhs( C8, C8b, tauw3t(k), Skw_zt(k) )
        endif

        if ( iwp3_dp1 > 0 ) then
          if ( l_crank_nich_diff ) then
            ! Eddy diffusion for wp3 using a Crank-Nicholson time step.
            tmp(1:3) & 
            = C12 * (1.0/2.0) & 
            * diffusion_zt_lhs( Kw8(k), Kw8(km1), nu8, & 
                                gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
          else
            ! Eddy diffusion for wp3 using a completely implicit time step.
            tmp(1:3) & 
            = C12  & 
            * diffusion_zt_lhs( Kw8(k), Kw8(km1), nu8, & 
                                gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
          endif

          ztscr02(k) = - tmp(3)
          ztscr03(k) = - tmp(2)
          ztscr04(k) = - tmp(1)

        endif

      endif

    enddo ! k = 2, gr%nnzp-1, 1


    ! Boundary conditions

    ! Both wp2 and wp3 used fixed-point boundary conditions.
    ! Therefore, anything set in the above loop at both the upper
    ! and lower boundaries would be overwritten here.  However, the
    ! above loop does not extend to the boundary levels.  An array
    ! with a value of 1 at the main diagonal on the left-hand side
    ! and with values of 0 at all other diagonals on the left-hand
    ! side will preserve the right-hand side value at that level.
    !
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

    ! w'^2
    lhs(:,k_wp2) = 0.0
    lhs(3,k_wp2) = 1.0
    ! w'^3
    lhs(:,k_wp3) = 0.0
    lhs(3,k_wp3) = 1.0

    ! Upper boundary
    k = gr%nnzp
    k_wp3 = 2*k - 1
    k_wp2 = 2*k

    ! w'^2
    lhs(:,k_wp2) = 0.0
    lhs(3,k_wp2) = 1.0
    ! w'^3
    lhs(:,k_wp3) = 0.0
    lhs(3,k_wp3) = 1.0


    return
  end subroutine wp23_lhs

  !=============================================================================
  !subroutine wp23_rhs( dt, wp2, wp3, wp3_zm, a1_zt, a1, a3, &
  subroutine wp23_rhs( dt, wp2, wp3, wp3_zm, a1_zt, &
                       a3_zt, wpthvp, wp2thvp, um, vm,  & 
                       upwp, vpwp, up2, vp2, Kw1, Kw8,  & 
                       Skw_zt, tau1m, tauw3t, C1_Skw_fnc, &
                       C11_Skw_fnc, l_crank_nich_diff, rhs )

    ! Description:
    ! Compute RHS vector for w'^2 and w'^3.
    ! This subroutine computes the explicit portion of 
    ! the w'^2 and w'^3 equations.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr ! Variable

    use parameters, only:  & 
        C4,  & ! Variables
        C5,  & 
        C8,  & 
        C8b,  & 
        C12,  & 
        nu1, & 
        nu8

    use constants, only: & 
        wtol,  & ! Variable(s) 
        emin,  &
        eps

    use model_flags, only:  & 
        l_tke_aniso ! Variable

    use diffusion, only: & 
        diffusion_zm_lhs,  & ! Procedures
        diffusion_zt_lhs

    use stats_precision, only:  & 
        time_precision ! Variable

    use stats_variables, only:  & 
        l_stats_samp, iwp2_dp1, iwp2_dp2, zm, iwp2_bp,   & ! Variable(s)
        iwp2_pr1, iwp2_pr2, iwp2_pr3, iwp3_ta, zt, & 
        iwp3_tp, iwp3_bp, iwp3_pr2, iwp3_pr1, iwp3_dp1

    use stats_type, only: stat_update_var_pt, stat_begin_update_pt ! Procedure(s)


    implicit none

! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      dt             ! Timestep length                                [s]

    real, dimension(gr%nnzp), intent(in) ::  & 
      wp2,         & ! w'^2 (momentum levels)                         [m^2/s^2]
      wp3,         & ! w'^3 (thermodynamic levels)                    [m^3/s^3]
      wp3_zm,      & ! w'^3 interpolated to momentum levels           [m^3/s^3]
!      a1,          & ! sigma_sqd_w-related term a_1 (momentum levels) [-]
!      a3,          & ! sigma_sqd_w-related term a_1 (momentum levels) [-]
      a1_zt,       & ! a_1 interpolated to thermodynamic levels       [-]
      a3_zt,       & ! a_3 interpolated to thermodynamic levels       [-]
      wpthvp,      & ! w'th_v' (momentum levels)                      [K m/s]
      wp2thvp,     & ! w'^2th_v' (thermodynamic levels)               [K m^2/s^2]
      um,          & ! u wind component (thermodynamic levels)        [m/s]
      vm,          & ! v wind component (thermodynamic levels)        [m/s]
      upwp,        & ! u'w' (momentum levels)                         [m^2/s^2]
      vpwp,        & ! v'w' (momentum levels)                         [m^2/s^2]
      up2,         & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2,         & ! v'^2 (momentum levels)                         [m^2/s^2]
      Kw1,         & ! Coefficient of eddy diffusivity for w'^2       [m^2/s]
      Kw8,         & ! Coefficient of eddy diffusivity for w'^3       [m^2/s]
      Skw_zt,      & ! Skewness of w on thermodynamic levels          [-]
      tau1m,       & ! Time-scale tau on momentum levels              [s]
      tauw3t,      & ! Time-scale tau on thermodynamic levels         [s]
      C1_Skw_fnc,  & ! C_1 parameter with Sk_w applied                [-]
      C11_Skw_fnc    ! C_11 parameter with Sk_w applied               [-]

    logical, intent(in) :: & 
      l_crank_nich_diff   ! Turns on/off Crank-Nicholson diffusion.

    ! Output Variable
    real, dimension(2*gr%nnzp), intent(out) :: & 
      rhs   ! RHS of band matrix

    ! Local Variables

    ! Array indices
    integer :: k, km1, kp1, k_wp2, k_wp3

    ! For use in Crank-Nicholson eddy diffusion.
    real, dimension(3) :: rhs_diff


    ! Initialize the right-hand side vector to 0.
    rhs = 0.0

    do k = 2, gr%nnzp-1, 1


      ! Define indices

      km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nnzp )

      k_wp3 = 2*k - 1
      k_wp2 = 2*k


      !!!!!***** w'^2 *****!!!!!

      ! w'^2: Right-hand side (explicit w'^2 portion of the code).

      ! RHS time tendency.
      rhs(k_wp2) & 
      = real( + ( 1.0 / dt ) * wp2(k) )

      ! RHS buoyancy production (bp) term and pressure term 2 (pr2).
      rhs(k_wp2) & 
      = rhs(k_wp2) & 
      + wp2_terms_bp_pr2_rhs( C5, wpthvp(k) )

      ! RHS pressure term 3 (pr3).
      rhs(k_wp2) & 
      = rhs(k_wp2) & 
      + wp2_term_pr3_rhs( C5, wpthvp(k), upwp(k), um(kp1), um(k), & 
                          vpwp(k), vm(kp1), vm(k), gr%dzm(k) )

      ! RHS dissipation term 1 (dp1).
      rhs(k_wp2) &
      = rhs(k_wp2) &
      + wp2_term_dp1_rhs( C1_Skw_fnc(k), tau1m(k), 2./3.*emin )

      ! RHS eddy diffusion term: dissipation term 2 (dp2).
      if ( l_crank_nich_diff ) then
        ! These lines are for the diffusional term with a Crank-Nicholson
        ! time step.  They are not used for completely implicit diffusion.
        rhs_diff(1:3) & 
        = (1.0/2.0) & 
        * diffusion_zm_lhs( Kw1(k), Kw1(kp1), nu1, & 
                            gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )
        rhs(k_wp2)   =   rhs(k_wp2) & 
                       - rhs_diff(3) * wp2(km1) & 
                       - rhs_diff(2) * wp2(k) & 
                       - rhs_diff(1) * wp2(kp1)
      endif

      ! RHS pressure term 1 (pr1).
      if ( l_tke_aniso ) then
        rhs(k_wp2) & 
        = rhs(k_wp2) & 
        + wp2_term_pr1_rhs( C4, up2(k), vp2(k), tau1m(k) )
      endif

      if ( l_stats_samp ) then

        ! Statistics: explicit contributions for wp2.

        if ( l_crank_nich_diff ) then
          call stat_begin_update_pt( iwp2_dp2, k, & 
            rhs_diff(3) * wp2(km1) & 
          + rhs_diff(2) * wp2(k) & 
          + rhs_diff(1) * wp2(kp1), zm )
        endif

        call stat_update_var_pt( iwp2_bp, k, & 
          wp2_terms_bp_pr2_rhs( 0.0, wpthvp(k) ), zm )


        if ( l_tke_aniso ) then
          call stat_begin_update_pt( iwp2_pr1, k, & 
            -wp2_term_pr1_rhs( C4, up2(k), vp2(k), tau1m(k) ), zm )

        endif

        call stat_begin_update_pt( iwp2_pr2, k, & 
          -wp2_terms_bp_pr2_rhs( (1.0+C5), wpthvp(k) ), zm )

        call stat_begin_update_pt( iwp2_dp1, k, &
          -wp2_term_dp1_rhs( C1_Skw_fnc(k), tau1m(k), 2./3.*emin ), zm )

        call stat_update_var_pt( iwp2_pr3, k, & 
          wp2_term_pr3_rhs( C5, wpthvp(k), upwp(k), um(kp1), um(k), & 
                        vpwp(k), vm(kp1), vm(k), gr%dzm(k) ), zm )

      endif



      !!!!!***** w'^3 *****!!!!!

      ! w'^3: Right-hand side (explicit w'^3 portion of the code).

      ! RHS time tendency.
      rhs(k_wp3) = & 
      real( + ( 1.0 / dt ) * wp3(k) )

      ! RHS turbulent advection (ta) and turbulent production (tp) terms.
      rhs(k_wp3) & 
      = rhs(k_wp3) & 
      + wp3_terms_ta_tp_rhs( wp3_zm(k), wp3_zm(km1),  &
                             wp2(k), wp2(km1),  &
                             !a1(k), a1(km1), a3(k), a3(km1),  &
                             a1_zt(k), a3_zt(k),  &
                             gr%dzt(k), wtol )

      ! RHS buoyancy production (bp) term and pressure term 2 (pr2).
      rhs(k_wp3) & 
      = rhs(k_wp3) & 
      + wp3_terms_bp_pr2_rhs( C11_Skw_fnc(k), wp2thvp(k) )

      ! RHS pressure term 1 (pr1).
      rhs(k_wp3) & 
      = rhs(k_wp3) & 
      + wp3_term_pr1_rhs( C8, C8b, tauw3t(k), Skw_zt(k), wp3(k) )

      ! RHS eddy diffusion term: dissipation term 1 (dp1).
      if ( l_crank_nich_diff ) then
        ! These lines are for the diffusional term with a Crank-Nicholson
        ! time step.  They are not used for completely implicit diffusion.
        rhs_diff(1:3) & 
        = C12 * (1.0/2.0) & 
        * diffusion_zt_lhs( Kw8(k), Kw8(km1), nu8, & 
                            gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
        rhs(k_wp3)   =   rhs(k_wp3) & 
                       - rhs_diff(3) * wp3(km1) & 
                       - rhs_diff(2) * wp3(k) & 
                       - rhs_diff(1) * wp3(kp1)
      endif

      if (l_stats_samp) then

        ! Statistics: explicit contributions for wp3.
        call stat_begin_update_pt( iwp3_ta, k, & 
          -wp3_terms_ta_tp_rhs( wp3_zm(k), wp3_zm(km1),  &
                                wp2(k), wp2(km1),  &
                                !a1(k), a1(km1),  &
                                !a3(k)+(3.0/2.0), a3(km1)+(3.0/2.0),  &
                                a1_zt(k), a3_zt(k)+(3.0/2.0),  &  
                                gr%dzt(k), wtol ), zt )

        call stat_begin_update_pt( iwp3_tp, k,  &
          -wp3_terms_ta_tp_rhs( wp3_zm(k), wp3_zm(km1),  &
                                wp2(k), wp2(km1),  &
                                !0.0, 0.0,  & 
                                !0.0-(3.0/2.0), 0.0-(3.0/2.0),  &
                                0.0, 0.0-(3.0/2.0),  & 
                                gr%dzt(k), wtol ),zt )

        call stat_update_var_pt( iwp3_bp, k, & 
          wp3_terms_bp_pr2_rhs( 0.0, wp2thvp(k) ), zt )

        call stat_begin_update_pt( iwp3_pr2, k, & 
          -wp3_terms_bp_pr2_rhs( (1.0+C11_Skw_fnc(k)), wp2thvp(k) ), & 
          zt )

        call stat_begin_update_pt( iwp3_pr1, k, & 
          -wp3_term_pr1_rhs( C8, C8b, tauw3t(k), Skw_zt(k), wp3(k) ), & 
          zt)

        if ( l_crank_nich_diff ) then
          call stat_begin_update_pt( iwp3_dp1, k, & 
            rhs_diff(3) * wp3(km1) & 
      + rhs_diff(2) * wp3(k) & 
      + rhs_diff(1) * wp3(kp1), zt )

        endif

      endif

    enddo


    ! Boundary conditions

    ! Both wp2 and wp3 used fixed-point boundary conditions.
    ! Therefore, anything set in the above loop at both the upper
    ! and lower boundaries would be overwritten here.  However, the
    ! above loop does not extend to the boundary levels.  An array
    ! with a value of 1 at the main diagonal on the left-hand side
    ! and with values of 0 at all other diagonals on the left-hand
    ! side will preserve the right-hand side value at that level.

    ! Lower boundary
    k = 1
    k_wp3 = 2*k - 1
    k_wp2 = 2*k

    ! The value of w'^2 at the lower boundary will remain the same.
    ! When the lower boundary is at the surface, the surface value of
    ! w'^2 is set in subroutine sfc_var (sfc.F).
    rhs(k_wp2)   = wp2(k)
    ! The value of w'^3 at the lower boundary will be 0.
    rhs(k_wp3)   = 0.0

    ! Upper boundary
    k = gr%nnzp
    k_wp3 = 2*k - 1
    k_wp2 = 2*k

! The value of w'^2 at the upper boundary will be 0.
    rhs(k_wp2)   = 0.0
! The value of w'^3 at the upper boundary will be 0.
    rhs(k_wp3)   = 0.0

    return

  end subroutine wp23_rhs

  !=============================================================================
  pure function wp2_term_ta_lhs( dzm ) & 
  result( lhs )

    ! Description:
    ! Turbulent advection term for w'^2:  implicit portion of the code.
    !
    ! The d(w'^2)/dt equation contains a turbulent advection term:
    !
    ! - d(w'^3)/dz.
    !
    ! The term is solved for completely implicitly, such that:
    !
    ! - d( w'^3(t+1) )/dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign 
    !        is reversed and the leading "-" in front of the term is changed 
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of w'^3 being used is from 
    ! the next timestep, which is being advanced to in solving the d(w'^2)/dt 
    ! and d(w'^3)/dt equations.
    !
    ! This term is discretized as follows:
    !
    ! While the values of w'^2 are found on the momentum levels, the values of 
    ! w'^3 are found on the thermodynamic levels.  The derivative of w'^3 is 
    ! taken over the intermediate (central) momentum level, yielding the desired
    ! results.
    !
    ! -------------------wp3p1--------------------------------- t(k+1)
    !
    ! =============================d(wp3)/dz=================== m(k)
    !
    ! -------------------wp3----------------------------------- t(k)
    !
    ! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes 
    ! zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for 
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! dzm(k) = 1 / ( zt(k+1) - zt(k) )

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2       ! Thermodynamic subdiagonal index.

    ! Input Variables
    real, intent(in) :: & 
      dzm     ! Inverse of grid spacing (k)   [1/m]

    ! Return Variable
    real, dimension(2) :: lhs

    ! Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
    lhs(kp1_tdiag) & 
    = + dzm

    ! Thermodynamic subdiagonal: [ x wp3(k,<t+1>) ]
    lhs(k_tdiag) & 
    = - dzm

    return

  end function wp2_term_ta_lhs

  !=============================================================================
  pure function wp2_terms_ac_pr2_lhs( C5, wm_ztp1, wm_zt, dzm ) & 
  result( lhs )

    ! Description:
    ! Accumulation of w'^2 and w'^2 pressure term 2:  implicit portion of the 
    ! code.
    !
    ! The d(w'^2)/dt equation contains an accumulation term:
    !
    ! - 2 w'^2 dw/dz;
    !
    ! and pressure term 2:
    !
    ! - C_5 ( -2 w'^2 dw/dz + 2 (g/th_0) w'th_v' ).
    !
    ! The w'^2 accumulation term is completely implicit, while w'^2 pressure 
    ! term 2 has both implicit and explicit components.  The accumulation term 
    ! and the implicit portion of pressure term 2 are combined and solved 
    ! together as:
    !
    ! + ( 1 - C_5 ) ( -2 w'^2(t+1) dw/dz ).
    !
    ! Note:  When the term is brought over to the left-hand side, the sign 
    !        is reversed and the leading "-" in front of the "2" is changed 
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of w'^2 being used is from 
    ! the next timestep, which is being advanced to in solving the d(w'^2)/dt 
    ! equation.
    !
    ! The terms are discretized as follows:
    !
    ! The values of w'^2 are found on the momentum levels, while the values of 
    ! wm_zt (mean vertical velocity on thermodynamic levels) are found on the 
    ! thermodynamic levels.  The vertical derivative of wm_zt is taken over the 
    ! intermediate (central) momentum level.  It is then multiplied by w'^2 
    ! (implicitly calculated at timestep (t+1)) and the coefficients to yield 
    ! the desired results.
    !
    ! -------wm_ztp1------------------------------------------- t(k+1)
    !
    ! ===============d(wm_zt)/dz============wp2================ m(k)
    !
    ! -------wm_zt--------------------------------------------- t(k)
    !
    ! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes 
    ! zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for 
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! dzm(k) = 1 / ( zt(k+1) - zt(k) )

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) :: & 
      C5,      & ! Model parameter C_5                            [-]
      wm_ztp1, & ! w wind component at thermodynamic levels (k+1) [m/s]
      wm_zt,   & ! w wind component at thermodynamic levels (k)   [m/s]
      dzm        ! Inverse of grid spacing (k)                    [1/m]

    ! Return Variable
    real :: lhs

    ! Momentum main diagonal: [ x wp2(k,<t+1>) ]
    lhs & 
    = + ( 1.0 - C5 ) * 2.0 * dzm * ( wm_ztp1 - wm_zt )

    return

  end function wp2_terms_ac_pr2_lhs

  !=============================================================================
  pure function wp2_term_dp1_lhs( C1_Skw_fnc, tau1m ) & 
  result( lhs )

    ! Description:
    ! Dissipation term 1 for w'^2:  implicit portion of the code.
    !
    ! The d(w'^2)/dt equation contains dissipation term 1:
    !
    ! - ( C_1 / tau_1m ) w'^2.
    !
    ! Since w'^2 has a minimum threshold, the term should be damped only to that
    ! threshold.  The term becomes:
    !
    ! - ( C_1 / tau_1m ) * ( w'^2 - threshold ).
    !
    ! This term is broken into implicit and explicit portions.  The implicit 
    ! portion of this term is:
    !
    ! - ( C_1 / tau_1m ) w'^2(t+1).
    !
    ! Note:  When the implicit term is brought over to the left-hand side, the
    !        sign is reversed and the leading "-" in front of the term is 
    !        changed to a "+".
    !
    ! The timestep index (t+1) means that the value of w'^2 being used is from 
    ! the next timestep, which is being advanced to in solving the d(w'^2)/dt 
    ! equation.
    !
    ! The values of w'^2 are found on the momentum levels.  The values of the 
    ! C_1 skewness function and time-scale tau1m are also found on the momentum 
    ! levels.

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) :: & 
      C1_Skw_fnc,  & ! C_1 parameter with Sk_w applied (k)   [-]
      tau1m          ! Time-scale tau at momentum levels (k) [s]

    ! Return Variable
    real :: lhs

    ! Momentum main diagonal: [ x wp2(k,<t+1>) ]
    lhs & 
    = + C1_Skw_fnc / tau1m

    return
  end function wp2_term_dp1_lhs

  !=============================================================================
  pure function wp2_term_pr1_lhs( C4, tau1m ) & 
  result( lhs )

    ! Description
    ! Pressure term 1 for w'^2:  implicit portion of the code.
    !
    ! The d(w'^2)/dt equation contains pressure term 1:
    !
    ! - ( C_4 / tau_1m ) * ( w'^2 - (2/3)*em ),
    !
    ! where em = (1/2) * ( w'^2 + u'^2 + v'^2 ).
    !
    ! This simplifies to:
    !
    ! - ( C_4 / tau_1m ) * (2/3) * w'^2
    !    + ( C_4 / tau_1m ) * (1/3) * ( u'^2 + v'^2 ).
    !
    ! Pressure term 1 has both implicit and explicit components.  The implicit
    ! portion is:
    !
    ! - ( C_4 / tau_1m ) * (2/3) * w'^2(t+1);
    !
    ! and is computed in this function.
    !
    ! Note:  When the implicit term is brought over to the left-hand side, the
    !        sign is reversed and the leading "-" in front of the term is 
    !        changed to a "+".
    !
    ! The timestep index (t+1) means that the value of w'^2 being used is from 
    ! the next timestep, which is being advanced to in solving the d(w'^2)/dt 
    ! equation.
    !
    ! The values of w'^2 are found on momentum levels, as are the values of tau1m.

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) :: & 
      C4,   & ! Model parameter C_4                   [-]
      tau1m   ! Time-scale tau at momentum levels (k) [s]

    ! Return Variable
    real :: lhs

    ! Momentum main diagonal: [ x wp2(k,<t+1>) ]
    lhs & 
    = + ( 2.0 * C4 ) / ( 3.0 * tau1m )

    return
  end function wp2_term_pr1_lhs

  !=============================================================================
  pure function wp2_terms_bp_pr2_rhs( C5, wpthvp ) & 
  result( rhs )

    ! Description:
    ! Buoyancy production of w'^2 and w'^2 pressure term 2:  explicit portion of
    ! the code.
    !
    ! The d(w'^2)/dt equation contains a buoyancy production term:
    !
    ! + 2 (g/th_0) w'th_v';
    !
    ! and pressure term 2:
    !
    ! - C_5 ( -2 w'^2 dw/dz + 2 (g/th_0) w'th_v' ).
    !
    ! The w'^2 buoyancy production term is completely explicit, while w'^2 
    ! pressure term 2 has both implicit and explicit components.  The buoyancy 
    ! production term and the explicit portion of pressure term 2 are combined 
    ! and solved together as:
    !
    ! + ( 1 - C_5 ) ( 2 (g/th_0) w'th_v' ).

    ! References:
    !-----------------------------------------------------------------------

    use constants, only:  & 
    ! Variable(s)        
        grav ! Gravitational acceleration [m/s^2]

    use parameters, only: & 
    ! Variable(s) 
        T0  ! Reference temperature      [K]

    implicit none

    ! Input Variables
    real, intent(in) :: & 
      C5,    & ! Model parameter C_5 [-]
      wpthvp   ! w'th_v'(k)          [K m/s]

    ! Return Variable
    real :: rhs

    rhs & 
    = + ( 1.0 - C5 ) * 2.0 * ( grav / T0 ) * wpthvp

    return
  end function wp2_terms_bp_pr2_rhs

  !=============================================================================
  pure function wp2_term_dp1_rhs( C1_Skw_fnc, tau1m, threshold ) & 
  result( rhs )

    ! Description:
    ! Dissipation term 1 for w'^2:  explicit portion of the code.
    !
    ! The d(w'^2)/dt equation contains dissipation term 1:
    !
    ! - ( C_1 / tau_1m ) w'^2.
    !
    ! Since w'^2 has a minimum threshold, the term should be damped only to that
    ! threshold.  The term becomes:
    !
    ! - ( C_1 / tau_1m ) * ( w'^2 - threshold ).
    !
    ! This term is broken into implicit and explicit portions.  The explicit 
    ! portion of this term is:
    !
    ! + ( C_1 / tau_1m ) * threshold.
    !
    ! The values of the C_1 skewness function, time-scale tau1m, and the 
    ! threshold are found on the momentum levels.

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) :: & 
      C1_Skw_fnc,  & ! C_1 parameter with Sk_w applied (k)   [-]
      tau1m,       & ! Time-scale tau at momentum levels (k) [s]
      threshold      ! Minimum allowable value of w'^2       [m^2/s^2]

    ! Return Variable
    real :: rhs

    rhs & 
    = + ( C1_Skw_fnc / tau1m ) * threshold

    return
  end function wp2_term_dp1_rhs

  !=============================================================================
  pure function wp2_term_pr3_rhs( C5, wpthvp, upwp, ump1, um, & 
                                  vpwp, vmp1, vm, dzm ) & 
  result( rhs )

    ! Description:
    ! Pressure term 3 for w'^2:  explicit portion of the code.
    !
    ! The d(w'^2)/dt equation contains pressure term 3:
    !
    ! + (2/3) C_5 [ (g/th_0) w'th_v' - u'w' du/dz - v'w' dv/dz ].
    !
    ! This term is solved for completely explicitly and is discretized as 
    ! follows:
    !
    ! The values of w'th_v', u'w', and v'w' are found on the momentum levels,
    ! whereas the values of um and vm are found on the thermodynamic levels.  
    ! The derivatives of both um and vm are taken over the intermediate 
    ! (central) momentum level.  All the remaining mathematical operations take 
    ! place at the central momentum level, yielding the desired result.
    !
    ! ------ump1------------vmp1------------------------------- t(k+1)
    !
    ! =upwp======d(um)/dz========d(vm)/dz===vpwp=====wpthvp==== m(k)
    !
    ! ------um--------------vm--------------------------------- t(k)
    !
    ! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes 
    ! zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for 
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! dzm(k) = 1 / ( zt(k+1) - zt(k) )

    ! References:
    !-----------------------------------------------------------------------

    use constants, only: & 
    ! Variables 
        grav ! Gravitational acceleration [m/s^2]

    use parameters, only: & 
    ! Variables 
        T0  ! Reference temperature      [K]

    implicit none

    ! Input Variables
    real, intent(in) :: & 
      C5,      & ! Model parameter C_5           [-]
      wpthvp,  & ! w'th_v'(k)                    [K m/s]
      upwp,    & ! u'w'(k)                       [m^2/s^2]
      ump1,    & ! um(k+1)                       [m/s]
      um,      & ! um(k)                         [m/s]
      vpwp,    & ! v'w'(k)                       [m^2/s^2]
      vmp1,    & ! vm(k+1)                       [m/s]
      vm,      & ! vm(k)                         [m/s]
      dzm        ! Inverse of grid spacing (k)   [1/m]

    ! Return Variable
    real :: rhs

    rhs & 
    ! Michael Falk, 2 August 2007
    ! Use the following code for standard mixing, with c_k=0.548:
    = + (2.0/3.0) * C5 & 
                  * ( ( grav / T0 ) * wpthvp & 
                      - upwp * dzm * ( ump1 - um ) & 
                      - vpwp * dzm * ( vmp1 - vm ) & 
                    )
     ! Use the following code for alternate mixing, with c_k=0.1 or 0.2
!    = + (2.0/3.0) * C5 &
!                  * ( ( grav / T0 ) * wpthvp &
!                      - 0. * upwp * dzm * ( ump1 - um ) &
!                      - 0. * vpwp * dzm * ( vmp1 - vm ) &
!                    )
!    eMFc

    return
  end function wp2_term_pr3_rhs

  !=============================================================================
  pure function wp2_term_pr1_rhs( C4, up2, vp2, tau1m ) & 
  result( rhs )

    ! Description:
    ! Pressure term 1 for w'^2:  explicit portion of the code.
    !
    ! The d(w'^2)/dt equation contains pressure term 1:
    !
    ! - ( C_4 / tau_1m ) * ( w'^2 - (2/3)*em );
    !
    ! where em = (1/2) * ( w'^2 + u'^2 + v'^2 ).
    !
    ! This simplifies to:
    !
    ! - ( C_4 / tau_1m ) * (2/3) * w'^2
    !    + ( C_4 / tau_1m ) * (1/3) * ( u'^2 + v'^2 ).
    !
    ! Pressure term 1 has both implicit and explicit components.
    ! The explicit portion is:
    !
    ! + ( C_4 / tau_1m ) * (1/3) * ( u'^2 + v'^2 );
    !
    ! and is computed in this function.
    !
    ! The values of u'^2 and v'^2 are found on momentum levels, as are the 
    ! values of tau1m.

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) :: & 
      C4,   & ! Model parameter C_4                   [-]
      up2,  & ! u'^2(k)                               [m^2/s^2]
      vp2,  & ! v'^2(k)                               [m^2/s^2]
      tau1m   ! Time-scale tau at momentum levels (k) [s]

    ! Return Variable
    real :: rhs

    rhs & 
    = + ( C4 * ( up2 + vp2 ) ) / ( 3.0 * tau1m )

    return
  end function wp2_term_pr1_rhs

  !=============================================================================
  pure function wp3_terms_ta_tp_lhs( wp3_zm, wp3_zmm1,  &
                                     wp2, wp2m1,  &
                                     !a1, a1m1, a3, a3m1,  &
                                     a1_zt, a3_zt,  &
                                     dzt, wtol, level )  &
  result( lhs )

    ! Description:
    ! Turbulent advection and turbulent production of w'^3:  implicit portion of
    ! the code.
    !
    ! The d(w'^3)/dt equation contains a turbulent advection term:
    !
    ! - d(w'^4)/dz;
    !
    ! and a turbulent production term:
    !
    ! + 3 w'^2 d(w'^2)/dz.
    !
    ! A substitution is made in order to close the turbulent advection term, 
    ! such that:
    !
    ! w'^4 = (a_3 + 3/2) * (w'^2)^2  +  a_1 * ( (w'^3)^2 / w'^2 );
    !
    ! where a_1 and a_3 are variables that are both functions of sigma_sqd_w.  
    ! The turbulent production term is rewritten as:
    !
    ! + 3 w'^2 d(w'^2)/dz = + (3/2) d( (w'^2)^2 )/dz.
    !
    ! The turbulent advection and turbulent production terms are combined as:
    !
    ! - d [ a_3 * (w'^2)^2 ] / dz
    !    - d [ a_1 * ( (w'^3)^2 / w'^2 ) ] / dz.
    !
    ! The (w'^2)^2 and (w'^3)^2 terms are both linearized, such that:
    !
    ! ( w'^2(t+1) )^2 = - ( w'^2(t) )^2  +  2 * w'^2(t) * w'^2(t+1);
    ! ( w'^3(t+1) )^2 = - ( w'^3(t) )^2  +  2 * w'^3(t) * w'^3(t+1);
    !
    ! which produces implicit and explicit portions of these terms.  The 
    ! implicit portion of these terms is:
    !
    ! - d [ a_3 * 2 * w'^2(t) * w'^2(t+1) ] / dz
    !    - d [ a_1 * ( 2 * w'^3(t) * w'^3(t+1) ) / w'^2(t) ] / dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign is
    !        reversed and the leading "-" in front of both d[ ] / dz terms is
    !        changed to a "+".
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is 
    ! being advanced to in solving the d(w'^3)/dt and d(w'^2)/dt equations.
    !
    ! The implicit portion of these terms is discretized as follows:
    !
    ! The values of w'^3 are found on the thermodynamic levels, while the values
    ! of w'^2, a_1, and a_3 are found on the momentum levels.  The variable w'^3
    ! is interpolated to the intermediate momentum levels.  The values of the
    ! mathematical expressions (called F and G here) within the dF/dz and dG/dz
    ! terms are computed on the momentum levels.  Then, the derivatives (d/dz) 
    ! of the expressions (F and G) are taken over the central thermodynamic 
    ! level, yielding the desired result.  In this function, the values of F 
    ! and G are as follows:
    !
    ! F = a_3(t) * 2 * w'^2(t) * w'^2(t+1); and
    !
    ! G = a_1(t) * ( 2 * w'^3(t) * w'^3(t+1) ) / w'^2(t).
    !
    !
    ! ----------------------wp3p1------------------------------ t(k+1)
    !
    ! ===a3====wp2====a1=========wp3(interp)=================== m(k)
    !
    ! ----------------------wp3----------------dF/dz---dG/dz--- t(k)
    !
    ! ===a3m1==wp2m1==a1m1=======wp3(interp)=================== m(k-1)
    !
    ! ----------------------wp3m1------------------------------ t(k-1)
    !
    ! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond 
    ! with altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1), respectively. 
    ! The letter "t" is used for thermodynamic levels and the letter "m" is 
    ! used for momentum levels.
    !
    ! dzt(k) = 1 / ( zm(k) - zm(k-1) )

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr ! Variable

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_mdiag   = 2,    & ! Momentum superdiagonal index.
      k_tdiag   = 3,    & ! Thermodynamic main diagonal index.
      km1_mdiag = 4,    & ! Momentum subdiagonal index. 
      km1_tdiag = 5       ! Thermodynamic subdiagonal index.

    integer, parameter :: & 
      t_above = 1,    & ! Index for upper thermodynamic level grid weight.
      t_below = 2       ! Index for lower thermodynamic level grid weight.

    ! Input Variables
    real, intent(in) ::  & 
      wp3_zm,      & ! w'^3 interpolated to momentum level (k)     [m^3/s^3]
      wp3_zmm1,    & ! w'^3 interpolated to momentum level (k-1)   [m^3/s^3]
      wp2,         & ! w'^2(k)                                     [m^2/s^2]
      wp2m1,       & ! w'^2(k-1)                                   [m^2/s^2]
!      a1,          & ! a_1(k)                                      [-]
!      a1m1,        & ! a_1(k-1)                                    [-]
!      a3,          & ! a_3(k)                                      [-]
!      a3m1,        & ! a_3(k-1)                                    [-]
      a1_zt,       & ! a_1 interpolated to thermodynamic level (k) [-]
      a3_zt,       & ! a_3 interpolated to thermodynamic level (k) [-]
      dzt,         & ! Inverse of grid spacing (k)                 [1/m]
      wtol           ! w wind component tolerance                  [m/s]

    integer, intent(in) :: & 
      level ! Central thermodynamic level (on which calculation occurs).

    ! Return Variable
    real, dimension(5) :: lhs

    ! Local Variables
    integer :: & 
      mk,    & ! Momentum level directly above central thermodynamic level.
      mkm1     ! Momentum level directly below central thermodynamic level.

    ! Momentum level (k) is between thermodynamic level (k+1)
    ! and thermodynamic level (k).
    mk = level

    ! Momentum level (k-1) is between thermodynamic level (k)
    ! and thermodynamic level (k-1).
    mkm1 = level - 1

    ! Brian tried a new discretization for the turbulent advection term, which
    ! contains the term:  - d [ a_1 * (w'^3)^2 / w'^2 ] / dz.  In order to help
    ! stabilize w'^3, a_1 has been pulled outside of the derivative.  On the
    ! left-hand side of the equation, this effects the thermodynamic 
    ! superdiagonal (kp1_tdiag), the thermodynamic main diagonal (k_tdiag), and 
    ! the thermodynamic subdiagonal (km1_tdiag).

    ! Additionally, the discretization of the turbulent advection term, which
    ! contains the term:  - d [ (a_3 + 3/2) * (w'^2)^2 ] / dz, has been altered 
    ! to pull a_3 outside of the derivative.  This was done in order to help 
    ! stabilize w'^3.  On the left-hand side of the equation, this effects the 
    ! momentum superdiagonal (k_mdiag) and the momentum subdiagonal (km1_mdiag).

    ! Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
    lhs(kp1_tdiag) & 
!    = + dzt * a1 * ( 2.0 * wp3_zm / max(wp2, wtol**2) ) &
!            * gr%weights_zt2zm(t_above,mk)
    = + a1_zt * dzt * ( 2.0 * wp3_zm / max(wp2, wtol**2) ) & 
              * gr%weights_zt2zm(t_above,mk)

    ! Momentum superdiagonal: [ x wp2(k,<t+1>) ]
    lhs(k_mdiag) & 
!    = + dzt * a3 * 2.0 * wp2
    = + a3_zt * dzt * 2.0 * wp2

    ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
    lhs(k_tdiag) & 
!    = + dzt &
!        * (   a1 * ( 2.0 * wp3_zm / max(wp2, wtol**2) ) &
!              * gr%weights_zt2zm(t_below,mk) &
!            - a1m1 * ( 2.0 * wp3_zmm1 / max(wp2m1, wtol**2) ) &
!              * gr%weights_zt2zm(t_above,mkm1) &
!          )
    = + a1_zt * dzt & 
        * (   ( 2.0 * wp3_zm / max(wp2, wtol**2) ) & 
              * gr%weights_zt2zm(t_below,mk) & 
            - ( 2.0 * wp3_zmm1 / max(wp2m1, wtol**2) ) & 
              * gr%weights_zt2zm(t_above,mkm1) & 
          )

    ! Momentum subdiagonal: [ x wp2(k-1,<t+1>) ]
    lhs(km1_mdiag) & 
!    = - dzt * a3m1 * 2.0 * wp2m1
    = - a3_zt * dzt * 2.0 * wp2m1

    ! Thermodynamic subdiagonal: [ x wp3(k-1,<t+1>) ]
    lhs(km1_tdiag) & 
!    = - dzt * a1m1 * ( 2.0 * wp3_zmm1 / max(wp2m1, wtol**2) ) &
!            * gr%weights_zt2zm(t_below,mkm1)
    = - a1_zt * dzt * ( 2.0 * wp3_zmm1 / max(wp2m1, wtol**2) ) & 
              * gr%weights_zt2zm(t_below,mkm1)

    ! End of code that pulls out a3.
    ! End of Brian's a1 change.  Feb. 14, 2008.

    return
  end function wp3_terms_ta_tp_lhs

  !=============================================================================
  pure function wp3_terms_ac_pr2_lhs( C11_Skw_fnc,  & 
                                      wm_zm, wm_zmm1, dzt ) & 
  result( lhs )

    ! Description:
    ! Accumulation of w'^3 and w'^3 pressure term 2:  implicit portion of the 
    ! code.
    !
    ! The d(w'^3)/dt equation contains an accumulation term:
    !
    ! - 3 w'^3 dw/dz;
    !
    ! and pressure term 2:
    !
    ! - C_11 ( -3 w'^3 dw/dz + 3 (g/th_0) w'^2th_v' ).
    !
    ! The w'^3 accumulation term is completely implicit, while w'^3 pressure 
    ! term 2 has both implicit and explicit components.  The accumulation term 
    ! and the implicit portion of pressure term 2 are combined and solved 
    ! together as:
    !
    ! + ( 1 - C_11 ) ( -3 w'^3(t+1) dw/dz ).
    !
    ! Note:  When the term is brought over to the left-hand side, the sign 
    !        is reversed and the leading "-" in front of the "3" is changed 
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of w'^3 being used is from 
    ! the next timestep, which is being advanced to in solving the d(w'^3)/dt 
    ! equation.
    !
    ! The terms are discretized as follows:
    !
    ! The values of w'^3 are found on thermodynamic levels, while the values of
    ! wm_zm (mean vertical velocity on momentum levels) are found on momentum
    ! levels.  The vertical derivative of wm_zm is taken over the intermediate
    ! (central) thermodynamic level.  It is then multiplied by w'^3 (implicitly
    ! calculated at timestep (t+1)) and the coefficients to yield the desired
    ! results.
    !
    ! =======wm_zm============================================= m(k)
    !
    ! ---------------d(wm_zm)/dz------------wp3---------------- t(k)
    !
    ! =======wm_zmm1=========================================== m(k-1)
    !
    ! The vertical indices m(k), t(k), and m(k-1) correspond with altitudes 
    ! zm(k), zt(k), and zm(k-1), respectively.  The letter "t" is used for 
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! dzt(k) = 1 / ( zm(k) - zm(k-1) )

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) :: & 
      C11_Skw_fnc,  & ! C_11 parameter with Sk_w applied (k)      [-]
      wm_zm,        & ! w wind component at momentum levels (k)   [m/s]
      wm_zmm1,      & ! w wind component at momentum levels (k-1) [m/s]
      dzt             ! Inverse of grid spacing (k)               [1/m]

    ! Return Variable
    real :: lhs

    ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
    lhs & 
    = + ( 1.0 - C11_Skw_fnc ) & 
        * 3.0 * dzt * ( wm_zm - wm_zmm1 )

    return
  end function wp3_terms_ac_pr2_lhs

  !=============================================================================
  pure function wp3_term_pr1_lhs( C8, C8b, tauw3t, Skw_zt ) & 
  result( lhs )

    ! Description:
    ! Pressure term 1 for w'^3:  implicit portion of the code.
    !
    ! Pressure term 1 is the term:
    !
    ! - (C_8/tau_w3t) * ( C_8b * Sk_wt^4 + 1 ) * w'^3;
    !
    ! where Sk_wt = w'^3 / (w'^2)^(3/2).
    !
    ! This term needs to be linearized, so function L(w'^3) is defined to be 
    ! equal to this term (pressure term 1), such that:
    !
    ! L(w'^3) = - (C_8/tau_w3t) * ( C_8b * (w'^3)^5 / (w'^2)^6 + w'^3 ).
    !
    ! A Taylor Series expansion (truncated after the first derivative term) of
    ! L(w'^3) around w'^3 = w'^3(t) is used to linearize pressure term 1.
    ! Evaluating L(w'^3) at w'^3(t+1):
    !
    ! L( w'^3(t+1) ) = L( w'^3(t) )
    !                  + ( d L(w'^3) / d w'^3 )|_(w'^3=w'^3(t))
    !                    * ( w'^3(t+1) - w'^3(t) ).
    !
    ! After evaluating the expression above, the term has become linearized.  It
    ! is broken down into implicit (LHS) and explicit (RHS) components.
    ! The implicit portion is:
    !
    ! - (C_8/tau_w3t) * ( 5 * C_8b * Sk_wt^4 + 1 ) * w'^3(t+1).
    !
    ! Note:  When the term is brought over to the left-hand side, the sign 
    !        is reversed and the leading "-" in front of the term is changed 
    !        to a "+".
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is 
    ! being advanced to in solving the d(w'^3)/dt equation.
    !
    ! The values of w'^3 are found on the thermodynamic levels, as are the 
    ! values of tau_w3t and Sk_wt (in Sk_wt, w'^3 is found on thermodynamic 
    ! levels and w'^2 is interpolated to thermodynamic levels).

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) :: & 
      C8,      & ! Model parameter C_8                        [-]
      C8b,     & ! Model parameter C_8b                       [-]
      tauw3t,  & ! Time-scale tau at thermodynamic levels (k) [s]
      Skw_zt     ! Skewness of w at thermodynamic levels (k)  [-]

    ! Return Variable
    real :: lhs

    ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
    lhs & 
    = + ( C8 / tauw3t ) * ( 5.0 * C8b * Skw_zt**4 + 1.0 )

    return
  end function wp3_term_pr1_lhs

  !=============================================================================
  pure function wp3_terms_ta_tp_rhs( wp3_zm, wp3_zmm1,  & 
                                     wp2, wp2m1,  & 
                                     !a1, a1m1, a3, a3m1,  &
                                     a1_zt, a3_zt,  & 
                                     dzt, wtol )  & 
  result( rhs )

    ! Description:
    ! Turbulent advection and turbulent production of wp3:  explicit portion of 
    ! the code.
    !
    ! The d(w'^3)/dt equation contains a turbulent advection term:
    !
    ! - d(w'^4)/dz;
    !
    ! and a turbulent production term:
    !
    ! + 3 w'^2 d(w'^2)/dz.
    !
    ! A substitution is made in order to close the turbulent advection term, 
    ! such that:
    !
    ! w'^4 = (a_3 + 3/2) * (w'^2)^2  +  a_1 * ( (w'^3)^2 / w'^2 );
    !
    ! where a_1 and a_3 are variables that are both functions of sigma_sqd_w.  
    ! The turbulent production term is rewritten as:
    !
    ! + 3 w'^2 d(w'^2)/dz = + (3/2) d( (w'^2)^2 )/dz.
    !
    ! The turbulent advection and turbulent production terms are combined as:
    !
    ! - d [ a_3 * (w'^2)^2 ] / dz
    !    - d [ a_1 * ( (w'^3)^2 / w'^2 ) ] / dz.
    !
    ! The (w'^2)^2 and (w'^3)^2 terms are both linearized, such that:
    !
    ! ( w'^2(t+1) )^2 = - ( w'^2(t) )^2  +  2 * w'^2(t) * w'^2(t+1);
    ! ( w'^3(t+1) )^2 = - ( w'^3(t) )^2  +  2 * w'^3(t) * w'^3(t+1);
    !
    ! which produces implicit and explicit portions of these terms.  The 
    ! explicit portion of these terms is:
    !
    ! + d [ a_3 * ( w'^2(t) )^2 ] / dz
    !    + d [ a_1 * ( w'^3(t) )^2 / w'^2(t) ] / dz.
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is 
    ! being advanced to in solving the d(w'^3)/dt and d(w'^2)/dt equations.
    !
    ! The explicit portion of these terms is discretized as follows:
    !
    ! The values of w'^3 are found on the thermodynamic levels, while the values
    ! of w'^2, a_1, and a_3 are found on the momentum levels.  The variable w'^3
    ! is interpolated to the intermediate momentum levels.  The values of the
    ! mathematical expressions (called F and G here) within the dF/dz and dG/dz
    ! terms are computed on the momentum levels.  Then, the derivatives (d/dz) 
    ! of the expressions (F and G) are taken over the central thermodynamic 
    ! level, yielding the desired result.  In this function, the values of F 
    ! and G are as follows:
    !
    ! F = a_3(t) * ( w'^2(t) )^2; and
    !
    ! G = a_1(t) * ( w'^3(t) )^2 / w'^2(t).
    !
    !
    ! ----------------------wp3p1------------------------------ t(k+1)
    !
    ! ===a3====wp2====a1=========wp3(interp)=================== m(k)
    !
    ! ----------------------wp3----------------dF/dz---dG/dz--- t(k)
    !
    ! ===a3m1==wp2m1==a1m1=======wp3(interp)=================== m(k-1)
    !
    ! ----------------------wp3m1------------------------------ t(k-1)
    !
    ! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond 
    ! with altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is used
    ! for momentum levels.
    !
    ! dzt(k) = 1 / ( zm(k) - zm(k-1) )

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) ::  & 
      wp3_zm,      & ! w'^3 interpolated to momentum level (k)     [m^3/s^3]
      wp3_zmm1,    & ! w'^3 interpolated to momentum level (k-1)   [m^3/s^3]
      wp2,         & ! w'^2(k)                                     [m^2/s^2]
      wp2m1,       & ! w'^2(k-1)                                   [m^2/s^2]
!      a1,          & ! a_1(k)                                      [-]
!      a1m1,        & ! a_1(k-1)                                    [-]
!      a3,          & ! a_3(k)                                      [-]
!      a3m1,        & ! a_3(k-1)                                    [-]
      a1_zt,       & ! a_1 interpolated to thermodynamic level (k) [-]
      a3_zt,       & ! a_3 interpolated to thermodynamic level (k) [-]
      dzt,         & ! Inverse of grid spacing (k)                 [1/m]
      wtol           ! w wind component tolerance                  [m/s]

    ! Return Variable
    real :: rhs

    ! Brian tried a new discretization for the turbulent advection term, which
    ! contains the term:  - d [ a_1 * (w'^3)^2 / w'^2 ] / dz.  In order to help
    ! stabilize w'^3, a_1 has been pulled outside of the derivative.  This 
    ! effects the right-hand side of the equation, as well as the left-hand 
    ! side.

    ! Additionally, the discretization of the turbulent advection term, which
    ! contains the term:  - d [ (a_3 + 3/2) * (w'^2)^2 ] / dz, has been altered 
    ! to pull a_3 outside of the derivative.  This was done in order to help 
    ! stabilize w'^3.  This effects the right-hand side of the equation, as well
    ! as the left-hand side.

    rhs & 
!    = + dzt &
!        * ( a3 * wp2**2 - a3m1 * wp2m1**2 ) &
!      + dzt &
!        * (   a1 * ( wp3_zm**2 / max(wp2, wtol**2) ) &
!            - a1m1 * ( wp3_zmm1**2 / max(wp2m1, wtol**2) ) &
!          )
    = + a3_zt * dzt & 
        * ( wp2**2 - wp2m1**2 ) & 
      + a1_zt * dzt & 
        * (   ( wp3_zm**2 / max(wp2, wtol**2) ) & 
            - ( wp3_zmm1**2 / max(wp2m1, wtol**2) ) & 
          )

    ! End of code that pulls out a3.
    ! End of Brian's a1 change.  Feb. 14, 2008.

    return
  end function wp3_terms_ta_tp_rhs

  !=============================================================================
  pure function wp3_terms_bp_pr2_rhs( C11_Skw_fnc, wp2thvp ) & 
  result( rhs )

    ! Description:
    ! Buoyancy production of w'^3 and w'^3 pressure term 2:  explicit portion of
    ! the code.
    !
    ! The d(w'^3)/dt equation contains a buoyancy production term:
    !
    ! + 3 (g/th_0) w'^2th_v';
    !
    ! and pressure term 2:
    !
    ! - C_11 ( -3 w'^3 dw/dz + 3 (g/th_0) w'^2th_v' ).
    !
    ! The w'^3 buoyancy production term is completely explicit, while w'^3 
    ! pressure term 2 has both implicit and explicit components.  The buoyancy 
    ! production term and the explicit portion of pressure term 2 are combined 
    ! and solved together as:
    !
    ! + ( 1 - C_ll ) ( 3 (g/th_0) w'^2th_v' ).

    ! References:
    !-----------------------------------------------------------------------

    use constants, only: & 
    ! Variable(s) 
        grav ! Gravitational acceleration [m/s^2]
    use parameters, only:  & 
    ! Variable(s)
        T0  ! Reference temperature      [K]

    implicit none

    ! Input Variables
    real, intent(in) :: & 
      C11_Skw_fnc,  & ! C_11 parameter with Sk_w applied (k)   [-]
      wp2thvp         ! w'^2th_v'(k)                           [K m^2/s^2]

    ! Return Variable
    real :: rhs

    rhs & 
    = + ( 1.0 - C11_Skw_fnc ) * 3.0 * ( grav/T0 ) * wp2thvp

    return
  end function wp3_terms_bp_pr2_rhs

  !=============================================================================
  pure function wp3_term_pr1_rhs( C8, C8b, tauw3t, Skw_zt, wp3 ) & 
  result( rhs )

    ! Description:
    ! Pressure term 1 for w'^3:  explicit portion of the code.
    !
    ! Pressure term 1 is the term:
    !
    ! - (C_8/tau_w3t) * ( C_8b * Sk_wt^4 + 1 ) * w'^3;
    !
    ! where Sk_wt = w'^3 / (w'^2)^(3/2).
    !
    ! This term needs to be linearized, so function L(w'^3) is defined to be 
    ! equal to this term (pressure term 1), such that:
    !
    ! L(w'^3) = - (C_8/tau_w3t) * ( C_8b * (w'^3)^5 / (w'^2)^6 + w'^3 ).
    !
    ! A Taylor Series expansion (truncated after the first derivative term) of
    ! L(w'^3) around w'^3 = w'^3(t) is used to linearize pressure term 1.
    ! Evaluating L(w'^3) at w'^3(t+1):
    !
    ! L( w'^3(t+1) ) = L( w'^3(t) )
    !                  + ( d L(w'^3) / d w'^3 )|_(w'^3=w'^3(t))
    !                    * ( w'^3(t+1) - w'^3(t) ).
    !
    ! After evaluating the expression above, the term has become linearized.  It
    ! is broken down into implicit (LHS) and explicit (RHS) components.
    ! The explicit portion is:
    !
    ! + (C_8/tau_w3t) * ( 4 * C_8b * Sk_wt^4 + 1 ) * w'^3(t).
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is 
    ! being advanced to in solving the d(w'^3)/dt equation.
    !
    ! The values of w'^3 are found on the thermodynamic levels, as are the 
    ! values of tau_w3t and Sk_wt (in Sk_wt, w'^3 is found on thermodynamic 
    ! levels and w'^2 is interpolated to thermodynamic levels).

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) :: & 
      C8,      & ! Model parameter C_8                        [-]
      C8b,     & ! Model parameter C_8b                       [-]
      tauw3t,  & ! Time-scale tau at thermodynamic levels (k) [s]
      Skw_zt,  & ! Skewness of w at thermodynamic levels (k)  [-]
      wp3        ! w'^3(k)                                    [m^3/s^3]

    ! Return Variable
    real :: rhs

    rhs & 
    = + ( C8 / tauw3t ) * ( 4.0 * C8b * Skw_zt**4 ) * wp3

    return
  end function wp3_term_pr1_rhs

!===============================================================================

end module wp23
