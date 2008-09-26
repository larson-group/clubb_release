!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module advance_xp2_xpyp_module

! Description:
!   Contains the subroutine advance_xp2_xpyp and ancillary functions.
!-----------------------------------------------------------------------

  implicit none

  public :: advance_xp2_xpyp

  private :: xp2_xpyp_lhs,  & 
             xp2_xpyp_solve,  & 
             xp2_xpyp_uv_rhs, & 
             xp2_xpyp_rhs, & 
             xp2_xpyp_implicit_stats, & 
             term_ta_lhs, & 
             term_ta_rhs, & 
             term_tp, & 
             term_dp1_lhs, & 
             term_dp1_rhs, & 
             term_pr1, & 
             term_pr2

  private    ! Set default scope

  contains

!===============================================================================
  subroutine advance_xp2_xpyp( tau_zm, wm_zm, rtm, wprtp, & 
                               thlm, wpthlp, wpthvp, um, vm, & 
                               wp2, wp2_zt, wp3, upwp, vpwp, &
                               sigma_sqd_w, Skw_zm, Kh_zt, & 
                               l_iter, dt, & 
                               sclrm, wpsclrp, & 
                               rtp2, thlp2, rtpthlp, & 
                               up2, vp2,  & 
                               err_code, & 
                               sclrp2, sclrprtp, sclrpthlp )

! Description:
!   Subprogram to diagnose variances by solving steady-state equations

! References:
!   Eqn. 13, 14, 15  on p. 3545 of
!   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!     Method and Model Description'' Golaz, et al. (2002)
!   JAS, Vol. 59, pp. 3540--3551.

! See also:
!   ``Equations for HOC'', Section 4:
!   /Steady-state solution for the variances/
!-----------------------------------------------------------------------

    use constants, only: & 
      wtol_sqd,  & ! Variable(s)
      rttol, & 
      thltol, & 
      emin, & 
      fstderr, &
      zero_threshold

    use model_flags, only: & 
      l_hole_fill, &    ! logical constants
      l_single_C2_Skw, &
      l_3pt_sqd_dfsn

    use parameters_tunable, only: &
      C2rt,     & ! Variable(s)
      C2thl,    &
      C2rtthl,  &
      c_K2,     &
      nu2,      &
      c_K9,     &
      nu9,      &
      c_Ksqd,   &
      beta,     &
      C4,       &
      C14,      &
      C5,       &
      T0,       &
      sclr_dim, &
      sclrtol,  &
      C2,       &
      C2b,      &
      C2c

    use grid_class, only: & 
      gr,  & ! Variable(s)
      zm2zt ! Procedure(s)

    use stats_precision, only:  & 
      time_precision ! Variable(s)

    use clip_explicit, only: & 
      clip_covariance,  & ! Procedure(s)
      clip_variance

    use error_code, only:  & 
      clubb_no_error,  & ! Variable(s)
      lapack_error,    & ! Procedure(s)
      clubb_at_least_debug_level

    use stats_type, only: & 
      stat_begin_update, & ! Procedure(s)
      stat_end_update

    use stats_variables, only: & 
      irtp2_bt,  & ! Variable(s)
      ithlp2_bt, & 
      irtpthlp_bt, & 
      ivp2_bt, & 
      iup2_bt, & 
      zm, & 
      l_stats_samp

    use array_index, only: &
      iisclr_rt, &
      iisclr_thl

    implicit none

    ! Input variables
    real, intent(in), dimension(gr%nnzp) ::  & 
      tau_zm,      & ! Tau on momentum grid                          [s]
      wm_zm,       & ! w wind on m                                   [m/s]
      rtm,         & ! Total water mixing ratio                      [kg/kg]
      wprtp,       & ! w' r_t'                                       [(m kg)/(s kg)]
      thlm,        & ! Liquid potential temp.                        [K]
      wpthlp,      & ! w' th_l'                                      [(m K)/s]
      wpthvp,      & ! w' th_v'                                      [(m K)/s]
      um,          & ! u wind                                        [m/s]
      vm,          & ! v wind                                        [m/s]
      wp2,         & ! w'^2                                          [m^2/s^2]
      wp3,         & ! w'^3                                          [m^3/s^3]
      upwp,        & ! u'w'                                          [m^2/s^2]
      vpwp,        & ! v'w'                                          [m^2/s^2]
      sigma_sqd_w, & ! sigma_sqd_w on momentum grid                  [-]
      Skw_zm,      & ! Skw on moment. grid                           [-]
      Kh_zt,       & ! Eddy diffusivity on t-lev.                    [m^2/s]
      wp2_zt         ! w'^2 interpolated to thermodynamic levels     [m^2/s^2]

    logical, intent(in) :: l_iter ! Whether variances are prognostic

    real(kind=time_precision), intent(in) :: &
      dt             ! Model timestep                                [s]

    ! Passive scalar input
    real, intent(in), dimension(gr%nnzp, sclr_dim) ::  & 
      sclrm, wpsclrp

    ! Input/Output variables
    ! An attribute of (inout) is also needed to import the value of the variances
    ! at the surface.  Brian.  12/18/05.
    real, intent(inout), dimension(gr%nnzp) ::  & 
      rtp2,    & ! r_t'^2                        [(kg/kg)^2]
      thlp2,   & ! th_l'^2                       [K^2]
      rtpthlp, & ! r_t' th_l'                    [(kg K)/kg]
      up2,     & ! u'^2                          [m^2/s^2]
      vp2        ! v'^2                          [m^2/s^2]

    ! Output variable for singular matrices
    integer, intent(out) :: err_code

    ! Passive scalar output
    real, intent(inout), dimension(gr%nnzp, sclr_dim) ::  & 
      sclrp2, sclrprtp, sclrpthlp

    ! Local Variables
    real, dimension(gr%nnzp) :: & 
      C2sclr_1d, C2rt_1d, C2thl_1d, C2rtthl_1d, C4_C14_1d

    real, dimension(gr%nnzp) :: & 
      a1 ! a_1 (momentum levels); See eqn. 24 in `Equations for HOC' [-]

    real, dimension(gr%nnzp) :: & 
      upwp_zt,    & ! u'w' interpolated to thermodynamic levels     [m^2/s^2]
      vpwp_zt,    & ! v'w' interpolated to thermodynamic levels     [m^2/s^2]
      wpsclrp_zt    ! w'sclr' interpolated to thermodynamic levels  [m/s {sclrm units}]

    real :: & 
      threshold     ! Minimum value for variances                   [units vary]

    real, dimension(3,gr%nnzp) ::  & 
      lhs ! Tridiagonal matrix

    real, dimension(gr%nnzp,1) :: & 
      rhs ! RHS vector of tridiagonal matrix

    real, dimension(gr%nnzp,2) :: & 
      uv_rhs,    &! RHS vectors of tridiagonal system for up2/vp2
      uv_solution ! Solution to the tridiagonal system for up2/vp2

    real, dimension(gr%nnzp,sclr_dim*3) ::  & 
      sclr_rhs,   & ! RHS vectors of tridiagonal system for the passive scalars
      sclr_solution ! Solution to tridiagonal system for the passive scalars

    integer, dimension(5+1) :: & 
      Valid_arr

    ! Eddy Diffusion for Variances and Covariances.
    real, dimension(gr%nnzp) ::  & 
      Kw2,      & ! For rtp2, thlp2, rtpthlp, and passive scalars  [m^2/s]
      Kw9         ! For up2 and vp2                                [m^2/s]

    real, dimension(gr%nnzp) :: & 
      a1_zt,      & ! a_1 interpolated to thermodynamic levels       [-]
      wprtp_zt,   & ! w'r_t' interpolated to thermodynamic levels    [(kg/kg) m/s]
      wpthlp_zt,  & ! w'th_l' interpolated to thermodyamnic levels   [K m/s]
      rtp2_zt,    & ! r_t'^2 interpolated to thermodynamic levels    [kg^2/kg^2]
      thlp2_zt,   & ! th_l'^2 interpolated to thermodynamic levels   [K^2]
      rtpthlp_zt, & ! r_t'th_l' interpolated to thermodynamic levels [K kg/kg]
      rtp2_zt_sqd_3pt, & 
      thlp2_zt_sqd_3pt, & 
      rtpthlp_zt_sqd_3pt, & 
      Kw2_rtp2, & 
      Kw2_thlp2, & 
      Kw2_rtpthlp

    logical :: l_scalar_calc

    ! Loop indices
    integer :: i
    integer :: k, km1, kp1

!---------------------------- Begin Code --------------------------------------


    if ( l_single_C2_Skw ) then
      ! Use a single value of C2 for all equations.
      C2rt_1d(1:gr%nnzp)  & 
      = C2b + (C2-C2b) *exp( -0.5 * (Skw_zm(1:gr%nnzp)/C2c)**2 )

      C2thl_1d   = C2rt_1d
      C2rtthl_1d = C2rt_1d

      C2sclr_1d  = C2rt_1d
    else
      ! Use 3 different values of C2 for rtp2, thlp2, rtpthlp.
      C2rt_1d(1:gr%nnzp)    = C2rt
      C2thl_1d(1:gr%nnzp)   = C2thl
      C2rtthl_1d(1:gr%nnzp) = C2rtthl

      C2sclr_1d(1:gr%nnzp)  = C2rt  ! Use rt value for now
    endif

    C4_C14_1d(1:gr%nnzp) = 2.0/3.0 * C4 + ( 1.0/3.0 * C14 )

    ! Are we solving for passive scalars as well?
    if ( sclr_dim > 0 ) then
      l_scalar_calc = .true.
    else
      l_scalar_calc = .false.
    endif


    ! Define a_1 (located on momentum levels).
    ! It is a variable that is a function of sigma_sqd_w (where sigma_sqd_w is
    ! located on the momentum levels).
    a1(1:gr%nnzp) = 1.0 / ( 1.0 - sigma_sqd_w(1:gr%nnzp) )


    ! Interpolate a_1, w'r_t', w'th_l', u'w', and v'w' from the momentum levels to
    ! the thermodynamic levels.  These will be used for the turbulent advection (ta)
    ! terms in each equation.
    a1_zt     = max( zm2zt( a1 ), zero_threshold )   ! Positive definite quantity
    wprtp_zt  = zm2zt( wprtp )
    wpthlp_zt = zm2zt( wpthlp )
    upwp_zt   = zm2zt( upwp )
    vpwp_zt   = zm2zt( vpwp )


    if ( l_stats_samp ) then

      call stat_begin_update( irtp2_bt, real(rtp2 / dt), &          ! Intent(in)
                              zm )                                  ! Intent(inout)

      call stat_begin_update( ithlp2_bt, real(thlp2 / dt), &        ! Intent(in)
                              zm )                                  ! Intent(inout)

      call stat_begin_update( irtpthlp_bt, real(rtpthlp / dt), &    ! Intent(in)
                              zm )                                  ! Intent(in/out)

      call stat_begin_update( ivp2_bt, real(vp2 / dt), &            ! Intent(in)
                              zm )                                  ! Intent(inout)

      call stat_begin_update( iup2_bt, real(up2 / dt),  &           ! Intent(in)
                              zm )                                  ! Intent(inout)

    endif

    ! Initialize tridiagonal solutions to valid

    Valid_arr(:) = clubb_no_error


    ! Define the Coefficent of Eddy Diffusivity for the variances and covariances.
    do k = 1, gr%nnzp, 1

      ! Kw2 is used for variances and covariances rtp2, thlp2, rtpthlp, and passive
      ! scalars.  The variances and covariances are located on the momentum levels.
      ! Kw2 is located on the thermodynamic levels.
      ! Kw2 = c_K2 * Kh_zt
      Kw2(k) = c_K2 * Kh_zt(k)

      ! Kw9 is used for variances up2 and vp2.  The variances are located on the
      ! momentum levels.  Kw9 is located on the thermodynamic levels.
      ! Kw9 = c_K9 * Kh_zt
      Kw9(k) = c_K9 * Kh_zt(k)

    enddo

    ! (xapxbp)^2: 3-point average diffusion coefficient.
    if ( l_3pt_sqd_dfsn ) then

      ! Interpolate r_t'^2, th_l'^2, and r_t'th_l' from the momentum levels to the
      ! thermodynamic levels.  These will be used for extra diffusion based on a
      ! three-point average of (var)^2.
      rtp2_zt    = max( zm2zt( rtp2 ), rttol**2 )  ! Positive def. quantity
      thlp2_zt   = max( zm2zt( thlp2 ), thltol**2 )  ! Positive def. quantity
      rtpthlp_zt = zm2zt( rtpthlp )

      do k = 1, gr%nnzp, 1

        km1 = max( k-1, 1 )
        kp1 = min( k+1, gr%nnzp )

        ! Compute the square of rtp2_zt, averaged over 3 points.  26 Jan 2008
        rtp2_zt_sqd_3pt(k) = ( rtp2_zt(km1)**2 + rtp2_zt(k)**2  & 
                               + rtp2_zt(kp1)**2 ) / 3.0
        ! Account for units (kg/kg)**4  Vince Larson 29 Jan 2008
        rtp2_zt_sqd_3pt(k) = 1e12 * rtp2_zt_sqd_3pt(k)

        ! Compute the square of thlp2_zt, averaged over 3 points.  26 Jan 2008
        thlp2_zt_sqd_3pt(k) = ( thlp2_zt(km1)**2 + thlp2_zt(k)**2  & 
                                + thlp2_zt(kp1)**2 ) / 3.0

        ! Compute the square of rtpthlp_zt, averaged over 3 points.  26 Jan 2008
        rtpthlp_zt_sqd_3pt(k) = ( rtpthlp_zt(km1)**2 + rtpthlp_zt(k)**2  & 
                                  + rtpthlp_zt(kp1)**2 ) / 3.0
        ! Account for units (kg/kg)**2 Vince Larson 29 Jan 2008
        rtpthlp_zt_sqd_3pt(k) = 1e6 * rtpthlp_zt_sqd_3pt(k)

      enddo

      ! Define Kw2_rtp2, Kw2_thlp2, and Kw2_rtpthlp
      do k = 1, gr%nnzp, 1

        ! Kw2_rtp2 must have units of m^2/s.  Since rtp2_zt_sqd_3pt has units of
        ! kg^2/kg^2, c_Ksqd is given units of m^2/[ s (kg^2/kg^2) ] in this case.
        Kw2_rtp2(k) = Kw2(k) + c_Ksqd * rtp2_zt_sqd_3pt(k)
        ! Vince Larson increased by c_Ksqd, 29Jan2008

        ! Kw2_thlp2 must have units of m^2/s.  Since thlp2_zt_sqd_3pt has units of
        ! K^2, c_Ksqd is given units of m^2/[ s K^2 ] in this case.
        Kw2_thlp2(k) = Kw2(k) + c_Ksqd * thlp2_zt_sqd_3pt(k)
        ! Vince Larson increased by c_Ksqd, 29Jan2008

        ! Kw2_rtpthlp must have units of m^2/s.  Since rtpthlp_zt_sqd_3pt has
        ! units of K (kg/kg), c_Ksqd is given units of m^2/[ s K (kg/kg) ] in this
        ! case.
        Kw2_rtpthlp(k) = Kw2(k) + c_Ksqd * rtpthlp_zt_sqd_3pt(k)
        ! Vince Larson increased by c_Ksqd, 29Jan2008

      enddo

    else  ! Three-point squared diffusion turned off.

      ! Define Kw2_rtp2, Kw2_thlp2, and Kw2_rtpthlp
      do k = 1, gr%nnzp, 1

        Kw2_rtp2(k)    = Kw2(k)
        Kw2_thlp2(k)   = Kw2(k)
        Kw2_rtpthlp(k) = Kw2(k)

      enddo

    endif  ! l_3pt_sqd_dfsn


    !!!!!***** r_t'^2 *****!!!!!

    ! Implicit contributions to term rtp2
    call xp2_xpyp_lhs( dt, l_iter, wp2_zt, wp3,  &             ! Intent(in)
                       a1, a1_zt, tau_zm, wm_zm, Kw2_rtp2, &   ! Intent(in)
                       C2rt_1d, nu2, beta,           &         ! Intent(in)
                       lhs )                                   ! Intent(out)


    call xp2_xpyp_rhs( "rtp2", dt, l_iter, a1, a1_zt, &     ! Intent(in)
                       wp2_zt, wp3, wprtp, wprtp_zt, &      ! Intent(in)
                       wprtp, wprtp_zt, rtm, rtm, rtp2, &   ! Intent(in)
                       C2rt_1d, tau_zm, rttol**2, beta, &   ! Intent(in)
                       rhs )                                ! Intent(out)

    ! Solve the tridiagonal system
    call xp2_xpyp_solve( "rtp2", 1, &                               ! Intent(in)
                         rhs, lhs, rtp2, &                          ! Intent(inout)
                         Valid_arr(1) )                             ! Intent(out)

    if ( l_stats_samp ) then
      call xp2_xpyp_implicit_stats( "rtp2", rtp2 ) ! Intent(in)
    end if

    !!!!!***** th_l'^2 *****!!!!!

    ! Implicit contributions to term thlp2
    call xp2_xpyp_lhs( dt, l_iter, wp2_zt, wp3,  &                  ! Intent(in)
                       a1, a1_zt, tau_zm, wm_zm, Kw2_thlp2,  &      ! Intent(in)
                       C2thl_1d, nu2, beta,           &             ! Intent(in)
                       lhs )                                        ! Intent(out)

    ! Explicit contributions to thlp2
    call xp2_xpyp_rhs( "thlp2", dt, l_iter, a1, a1_zt, &            ! Intent(in)
                       wp2_zt, wp3, wpthlp, wpthlp_zt, &            ! Intent(in)
                       wpthlp, wpthlp_zt, thlm, thlm, thlp2, &      ! Intent(in)
                       C2thl_1d, tau_zm, thltol**2, beta, &         ! Intent(in)
                       rhs )                                        ! Intent(out)

    ! Solve the tridiagonal system
    call xp2_xpyp_solve( "thlp2", 1, &          ! Intent(in)
                         rhs, lhs, thlp2, &     ! Intent(inout)
                         Valid_arr(2) )         ! Intent(out)

    if ( l_stats_samp ) then
      call xp2_xpyp_implicit_stats( "thlp2", thlp2 ) ! Intent(in)
    end if


    !!!!!***** r_t'th_l' *****!!!!!

    ! Implicit contributions to term rtpthlp
    call xp2_xpyp_lhs( dt, l_iter, wp2_zt, wp3,  &                  ! Intent(in)
                       a1, a1_zt, tau_zm, wm_zm, Kw2_rtpthlp,  &    ! Intent(in)
                       C2rtthl_1d, nu2, beta,           &           ! Intent(in)
                       lhs )                                        ! Intent(out)

    ! Explicit contributions to rtpthlp
    call xp2_xpyp_rhs( "rtpthlp", dt, l_iter, a1, a1_zt, &          ! Intent(in)
                       wp2_zt, wp3, wprtp, wprtp_zt, &              ! Intent(in)
                       wpthlp, wpthlp_zt, rtm, thlm, rtpthlp, &     ! Intent(in)
                       C2rtthl_1d, tau_zm, 0.0, beta, &             ! Intent(in)
                       rhs )                                        ! Intent(out)

    ! Solve the tridiagonal system
    call xp2_xpyp_solve( "rtpthlp", 1, &            ! Intent(in)
                         rhs, lhs, rtpthlp, &       ! Intent(inout)
                         Valid_arr(3) )             ! Intent(out)

    if ( l_stats_samp ) then
      call xp2_xpyp_implicit_stats( "rtpthlp", rtpthlp ) ! Intent(in)
    end if


    !!!!!***** u'^2 / v'^2 *****!!!!!

    ! Implicit contributions to term up2/vp2
    call xp2_xpyp_lhs( dt, l_iter, wp2_zt, wp3,  &             ! Intent(in)
                       a1, a1_zt, tau_zm, wm_zm, Kw9,  &       ! Intent(in)
                       C4_C14_1d, nu9, beta,           &       ! Intent(in)
                       lhs )                                   ! Intent(out)

    ! Explicit contributions to up2
    call xp2_xpyp_uv_rhs( "up2", dt, l_iter, a1, a1_zt, &       ! Intent(in)
                          wp2, wp2_zt, wp3, wpthvp, tau_zm,  &  ! Intent(in)
                          um, vm, upwp, upwp_zt, vpwp, &        ! Intent(in)
                          vpwp_zt, up2, vp2, C4, C5, C14, &     ! Intent(in)
                          T0, beta,           &                 ! Intent(in)
                          uv_rhs(:,1) )                         ! Intent(out)

    ! Explicit contributions to vp2
    call xp2_xpyp_uv_rhs( "vp2", dt, l_iter, a1, a1_zt, &       ! Intent(in)
                          wp2, wp2_zt, wp3, wpthvp, tau_zm,  &  ! Intent(in)
                          vm, um, vpwp, vpwp_zt, upwp, &        ! Intent(in)
                          upwp_zt, vp2, up2, C4, C5, C14, &     ! Intent(in)
                          T0, beta,           &                 ! Intent(in)
                          uv_rhs(:,2) )                         ! Intent(out)

    ! Solve the tridiagonal system
    call xp2_xpyp_solve( "up2_vp2", 2,              & ! Intent(in)
                         uv_rhs, lhs,               & ! Intent(inout)
                         uv_solution, Valid_arr(4) )  ! Intent(out)

    up2(1:gr%nnzp) = uv_solution(1:gr%nnzp,1)
    vp2(1:gr%nnzp) = uv_solution(1:gr%nnzp,2)

    if ( l_stats_samp ) then
      call xp2_xpyp_implicit_stats( "up2", up2 ) ! Intent(in)
      call xp2_xpyp_implicit_stats( "vp2", vp2 ) ! Intent(in)
    end if


    ! Apply the positive definite scheme to variances
    if ( l_hole_fill ) then
      call pos_definite_variances( "rtp2", dt, rttol**2, &   ! Intent(in)
                                   rtp2 )                    ! Intent(inout)
      call pos_definite_variances( "thlp2", dt, thltol**2, & ! Intent(in)
                                   thlp2 )                   ! Intent(inout)
      call pos_definite_variances( "up2", dt, 2./3.*emin, &  ! Intent(in)
                                   up2 )                     ! Intent(inout)
      call pos_definite_variances( "vp2", dt, 2./3.*emin, &  ! Intent(in)
                                   vp2 )                     ! Intent(inout)
    endif


    ! Clipping for r_t'^2

    !threshold = 0.0
    !
    !where ( wp2 >= wtol_sqd ) &
    !   threshold = rttol*rttol

    threshold = rttol**2

    call clip_variance( "rtp2", dt, threshold, & ! Intent(in)
                        rtp2 )                   ! Intent(inout)


    ! Clipping for th_l'^2

    !threshold = 0.0
    !
    !where ( wp2 >= wtol_sqd ) &
    !   threshold = thltol*thltol

    threshold = thltol**2

    call clip_variance( "thlp2", dt, threshold, & ! Intent(in)
                        thlp2 )                   ! Intent(inout)


    ! Clipping for u'^2

    !threshold = 0.0

    threshold = 2./3.*emin

    call clip_variance( "up2", dt, threshold, & ! Intent(in)
                        up2 )                   ! Intent(inout)


    ! Clipping for v'^2

    !threshold = 0.0
    threshold = 2./3.*emin

    call clip_variance( "vp2", dt, threshold, & ! Intent(in)
                        vp2 )                   ! Intent(inout)


    ! Clipping for r_t'th_l'
    ! Clipping r_t'th_l' at each vertical level, based on the
    ! correlation of r_t and th_l at each vertical level, such that:
    ! corr_(r_t,th_l) = r_t'th_l' / [ sqrt(r_t'^2) * sqrt(th_l'^2) ];
    ! -1 <= corr_(r_t,th_l) <= 1.
    ! Since r_t'^2, th_l'^2, and r_t'th_l' are all computed in the
    ! same place, clipping for r_t'th_l' only has to be done once.
    call clip_covariance( "rtpthlp", .true.,  &        ! Intent(in)
                          .true., dt, rtp2, thlp2,  &  ! Intent(in)
                          rtpthlp )                    ! Intent(inout)


    if ( l_stats_samp ) then

      call stat_end_update( irtp2_bt, real( rtp2 / dt), & ! Intent(in)
                            zm )                          ! Intent(inout)

      call stat_end_update( ithlp2_bt, real( thlp2 / dt), & ! Intent(in) 
                            zm )                            ! Intent(inout)

      call stat_end_update( irtpthlp_bt, real( rtpthlp / dt), & ! Intent(in)
                            zm )                                ! Intent(inout)

      call stat_end_update( iup2_bt, real( up2 / dt), & ! Intent(in)
                            zm )                        ! Intent(inout)

      call stat_end_update( ivp2_bt, real( vp2 / dt),& ! Intent(in)
                            zm )                       ! Intent(inout)

    endif


    if ( l_scalar_calc ) then

      ! Implicit contributions to passive scalars

      !!!!!***** sclr'^2, sclr'r_t', sclr'th_l' *****!!!!!

      call xp2_xpyp_lhs( dt, l_iter, wp2_zt, wp3,  &        ! Intent(in) 
                         a1, a1_zt, tau_zm, wm_zm, Kw2,  &  ! Intent(in)
                         C2sclr_1d, nu2, beta,           &  ! Intent(in)
                         lhs )                              ! Intent(out)


      ! Explicit contributions to passive scalars

      do i = 1, sclr_dim, 1

        ! Interpolate w'sclr' from momentum levels to thermodynamic
        ! levels.  These will be used for the turbulent advection (ta)
        ! terms in each equation.
        wpsclrp_zt = zm2zt( wpsclrp(:,i) )

        !!!!!***** sclr'^2 *****!!!!!

        call xp2_xpyp_rhs( "sclrp2", dt, l_iter, a1, a1_zt, &       ! Intent(in)
                           wp2_zt, wp3, wpsclrp(:,i),  &            ! Intent(in)
                           wpsclrp_zt, wpsclrp(:,i), wpsclrp_zt,  & ! Intent(in)
                           sclrm(:,i), sclrm(:,i), sclrp2(:,i), &   ! Intent(in)
                           C2sclr_1d, tau_zm, 0.0, beta, &          ! Intent(in)
                           sclr_rhs(:,i) )                          ! Intent(out)


        !!!!!***** sclr'r_t' *****!!!!!
        if ( i == iisclr_rt ) then
          ! In this case we're trying to emulate rt'^2 with sclr'rt', so we 
          ! handle this as we would a variance, even though generally speaking
          ! the scalar is not rt
          threshold = rttol**2
        else
          threshold = 0.0
        end if

        call xp2_xpyp_rhs( "sclrprtp", dt, l_iter, a1, a1_zt, &        ! Intent(in)
                           wp2_zt, wp3, wpsclrp(:,i),  &               ! Intent(in)
                           wpsclrp_zt, wprtp, wprtp_zt, sclrm(:,i),  & ! Intent(in)
                           rtm, sclrprtp(:,i), C2sclr_1d, tau_zm, &    ! Intent(in)     
                           threshold, beta, &                          ! Intent(in)
                           sclr_rhs(:,i+sclr_dim) )                    ! Intent(out)


        !!!!!***** sclr'th_l' *****!!!!!

        if ( i == iisclr_thl ) then
          ! In this case we're trying to emulate thl'^2 with sclr'thl', so we
          ! handle this as we did with sclr_rt, above.
          threshold = thltol**2
        else
          threshold = 0.0
        end if

        call xp2_xpyp_rhs( "sclrpthlp", dt, l_iter, a1, a1_zt, & ! Intent(in)
                           wp2_zt, wp3, wpsclrp(:,i),  &         ! Intent(in)
                           wpsclrp_zt, wpthlp, wpthlp_zt,  &     ! Intent(in)
                           sclrm(:,i), thlm, sclrpthlp(:,i), &   ! Intent(in)
                           C2sclr_1d, tau_zm, threshold, beta, & ! Intent(in)
                           sclr_rhs(:,i+2*sclr_dim) )            ! Intent(out)
      end do ! 1..sclr_dim


      ! Solve the tridiagonal system

      call xp2_xpyp_solve( "scalars", 3*sclr_dim, &         ! Intent(in)
                           sclr_rhs, lhs, sclr_solution, &  ! Intent(inout)
                           Valid_arr(6) )                   ! Intent(out)

      sclrp2(:,1:sclr_dim) = sclr_solution(:,1:sclr_dim)

      sclrprtp(:,1:sclr_dim) = sclr_solution(:,sclr_dim+1:2*sclr_dim)

      sclrpthlp(:,1:sclr_dim) = sclr_solution(:,2*sclr_dim+1:3*sclr_dim)

      ! Apply hole filling algorithm to the scalar variance terms
      if ( l_hole_fill ) then
        do i = 1, sclr_dim, 1
          call pos_definite_variances( "sclrp2", dt, sclrtol(i)**2, & ! Intent(in)
                                       sclrp2(:,i) )               ! Intent(inout)
          if ( i == iisclr_rt ) then 
             ! Here again, we do this kluge here to make sclr'rt' == rt'^2
            call pos_definite_variances( "sclrprtp", dt, sclrtol(i)**2, & ! Intent(in)
                                         sclrprtp(:,i) )               ! Intent(inout)
          end if
          if ( i == iisclr_thl ) then
            ! As with sclr'rt' above, but for sclr'thl'
            call pos_definite_variances( "sclrpthlp", dt, sclrtol(i)**2, & ! Intent(in)
                                         sclrpthlp(:,i) )               ! Intent(inout)
          end if
        enddo
      endif


      ! Clipping for sclr'^2
      do i = 1, sclr_dim, 1

!     threshold = 0.0
!
!     where ( wp2 >= wtol_sqd ) &
!        threshold = sclrtol(i)*sclrtol(i)

        threshold = sclrtol(i)**2

        call clip_variance( "sclrp2", dt, threshold, & ! Intent(in)
                            sclrp2(:,i) )              ! Intent(inout)

      enddo


      ! Clipping for sclr'r_t'
      ! Clipping sclr'r_t' at each vertical level, based on the
      ! correlation of sclr and r_t at each vertical level, such that:
      ! corr_(sclr,r_t) = sclr'r_t' / [ sqrt(sclr'^2) * sqrt(r_t'^2) ];
      ! -1 <= corr_(sclr,r_t) <= 1.
      ! Since sclr'^2, r_t'^2, and sclr'r_t' are all computed in the
      ! same place, clipping for sclr'r_t' only has to be done once.
      do i = 1, sclr_dim, 1

        if  ( i == iisclr_rt ) then
          ! Treat this like a variance if we're emulating rt
          threshold = sclrtol(i) * rttol

          call clip_variance( "sclrprtp", dt, threshold, & ! Intent(in)
                              sclrprtp(:,i) )              ! Intent(inout)
        else
          call clip_covariance( "sclrprtp", .true.,  &               ! Intent(in) 
                                .true., dt, sclrp2(:,i), rtp2(:), &  ! Intent(in)
                                sclrprtp(:,i) )                      ! Intent(inout)
        end if
      enddo


      ! Clipping for sclr'th_l'
      ! Clipping sclr'th_l' at each vertical level, based on the
      ! correlation of sclr and th_l at each vertical level, such that:
      ! corr_(sclr,th_l) = sclr'th_l' / [ sqrt(sclr'^2) * sqrt(th_l'^2) ];
      ! -1 <= corr_(sclr,th_l) <= 1.
      ! Since sclr'^2, th_l'^2, and sclr'th_l' are all computed in the
      ! same place, clipping for sclr'th_l' only has to be done once.
      do i = 1, sclr_dim, 1
        if ( i == iisclr_thl ) then
          ! As above, but for thl
          threshold = sclrtol(i) * thltol
          call clip_variance( "sclrpthlp", dt, threshold, & ! Intent(in)
                              sclrpthlp(:,i) )              ! Intent(inout)
        else

          call clip_covariance( "sclrpthlp", .true.,  &              ! Intent(in) 
                                .true., dt, sclrp2(:,i), thlp2(:), & ! Intent(in) 
                                sclrpthlp(:,i) )                     ! Intent(inout)
        end if
      enddo

    endif ! l_scalar_calc


    ! Check for singular matrices
    do i = 1, 5+1
      if( Valid_arr(i) > err_code ) err_code = Valid_arr(i)
    enddo

    if ( lapack_error( err_code ) .and.  & 
         clubb_at_least_debug_level( 1 ) ) then

      write(fstderr,*) "Error in advance_xp2_xpyp"

      write(fstderr,*) "Intent(in)"

      write(fstderr,*) "tau_zm = ", tau_zm
      write(fstderr,*) "wm_zm = ", wm_zm
      write(fstderr,*) "rtm = ", rtm
      write(fstderr,*) "wprtp = ", wprtp
      write(fstderr,*) "thlm = ", thlm
      write(fstderr,*) "wpthlp = ", wpthlp
      write(fstderr,*) "wpthvp = ", wpthvp
      write(fstderr,*) "um = ", um
      write(fstderr,*) "vm = ", vm
      write(fstderr,*) "wp2 = ", wp2
      write(fstderr,*) "wp3 = ", wp3
      write(fstderr,*) "upwp = ", upwp
      write(fstderr,*) "vpwp = ", vpwp
      write(fstderr,*) "sigma_sqd_w = ", sigma_sqd_w
      write(fstderr,*) "Skw_zm = ", Skw_zm
      write(fstderr,*) "Kh_zt = ", Kh_zt

      do i = 1, sclr_dim
        write(fstderr,*) "sclrm = ", i, sclrm(:,i)
        write(fstderr,*) "wpsclrp = ", i, wpsclrp(:,i)
      enddo

      write(fstderr,*) "Intent(In/Out)"

      write(fstderr,*) "rtp2 = ", rtp2
      write(fstderr,*) "thlp2 = ", thlp2
      write(fstderr,*) "rtpthlp = ", rtpthlp
      write(fstderr,*) "up2 = ", up2
      write(fstderr,*) "vp2 = ", vp2

      do i = 1, sclr_dim
        write(fstderr,*) "sclrp2 = ", i, sclrp2(:,i)
        write(fstderr,*) "sclrprtp = ", i, sclrprtp(:,i)
        write(fstderr,*) "sclrthlp = ", i, sclrpthlp(:,i)
      enddo

    endif

    return
  end subroutine advance_xp2_xpyp

!===============================================================================
  subroutine xp2_xpyp_lhs( dt, l_iter, wp2_zt, wp3,  & 
                           a1, a1_zt, tau_zm, wm_zm, Kw,  &
                           Cn, nu, beta, lhs )

! Description:
! Compute LHS tridiagonal matrix for a variance or coveriance term

! References:
! None
!-----------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use stats_precision, only:  & 
        time_precision ! Variable(s)

    use diffusion, only:  & 
        diffusion_zm_lhs ! Procedure(s)

    use mean_adv, only:  & 
        term_ma_zm_lhs ! Procedure(s)

    use stats_variables, only: & 
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
        l_stats_samp, & 
        irtp2_ma, & 
        irtp2_ta, & 
        irtp2_dp1, & 
        irtp2_dp2, & 
        ithlp2_ma, & 
        ithlp2_ta, & 
        ithlp2_dp1, & 
        ithlp2_dp2, & 
        irtpthlp_ma, & 
        irtpthlp_ta, & 
        irtpthlp_dp1, & 
        irtpthlp_dp2, & 
        iup2_ma, & 
        iup2_ta, & 
        iup2_dp2, & 
        ivp2_ma, & 
        ivp2_ta, & 
        ivp2_dp2


    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_mdiag = 1, & ! Momentum superdiagonal index.
      k_mdiag   = 2, & ! Momentum main diagonal index.
      km1_mdiag = 3    ! Momentum subdiagonal index.

    real(kind=time_precision), intent(in) :: & 
      dt        ! Timestep length                                [s]

    logical, intent(in) :: & 
      l_iter  ! Whether the variances are prognostic

    ! Input Variables
    real, dimension(gr%nnzp), intent(in) :: & 
      wp2_zt,  & ! w'^2 interpolated to thermodynamic levels      [m^2/s^2]
      wp3,     & ! w'^3 (thermodynamic levels)                    [m^3/s^3]
      a1,      & ! sigma_sqd_w-related term a_1 (momentum levels) [-]
      a1_zt,   & ! a_1 interpolated to thermodynamic levels       [-]
      tau_zm,  & ! Time-scale tau on momentum levels              [s]
      wm_zm,   & ! w wind component on momentum levels            [m/s]
      Kw,      & ! Coefficient of eddy diffusivity (all vars.)    [m^2/s]
      Cn         ! Coefficient C_n                                [-]

    real, intent(in) :: & 
      nu,      & ! Background constant coef. of eddy diff.        [-]
      beta       ! Constant model parameter beta                  [-]

    ! Output Variables
    real, dimension(3,gr%nnzp), intent(out) :: & 
      lhs ! Implicit contributions to the term

    ! Local Variables

    ! Array indices
    integer :: k, kp1
    !integer :: km1

    real, dimension(3) :: & 
      tmp


    ! Setup LHS of the tridiagonal system
    do k = 2, gr%nnzp-1, 1

      !km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nnzp )

      ! LHS eddy diffusion term: dissipation term 2 (dp2).
      lhs(kp1_mdiag:km1_mdiag,k) & 
      = diffusion_zm_lhs( Kw(k), Kw(kp1), nu, & 
                          gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )

      ! LHS dissipation term 1 (dp1)
      ! (and pressure term 1 (pr1) for u'^2 and v'^2).
      lhs(k_mdiag,k) & 
      = lhs(k_mdiag,k) + term_dp1_lhs( Cn(k), tau_zm(k) )

      ! LHS time tendency.
      if ( l_iter ) then
        lhs(k_mdiag,k) = real( lhs(k_mdiag,k) + ( 1.0 / dt ) )
      endif

      ! LHS mean advection (ma) term.
      lhs(kp1_mdiag:km1_mdiag,k) & 
      = lhs(kp1_mdiag:km1_mdiag,k) & 
      + term_ma_zm_lhs( wm_zm(k), gr%dzm(k), k )

      ! LHS turbulent advection (ta) term.
      lhs(kp1_mdiag:km1_mdiag,k) & 
      = lhs(kp1_mdiag:km1_mdiag,k) & 
      + term_ta_lhs( wp3(kp1), wp3(k), wp2_zt(kp1), wp2_zt(k),  &
                     a1_zt(kp1), a1(k), a1_zt(k), gr%dzm(k), beta, k )

      if ( l_stats_samp ) then

        ! Statistics: implicit contributions for rtp2, thlp2,
        !             rtpthlp, up2, or vp2.

        if ( irtp2_dp1 + ithlp2_dp1 + irtpthlp_dp1  > 0 ) then
          ! Note:  The statistical implicit contribution to term dp1
          !        (as well as to term pr1) for up2 and vp2 is recorded
          !        in xp2_xpyp_uv_rhs because up2 and vp2 use a special
          !        dp1/pr1 combined term.
          tmp(1) = term_dp1_lhs( Cn(k), tau_zm(k) )
          zmscr01(k) = -tmp(1)
        endif

        if ( irtp2_dp2 + ithlp2_dp2 + irtpthlp_dp2 +  & 
             iup2_dp2 + ivp2_dp2 > 0 ) then
          tmp(1:3) & 
          = diffusion_zm_lhs( Kw(k), Kw(kp1), nu, & 
                              gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )
          zmscr02(k) = -tmp(3)
          zmscr03(k) = -tmp(2)
          zmscr04(k) = -tmp(1)
        endif

        if ( irtp2_ta + ithlp2_ta + irtpthlp_ta + & 
             iup2_ta + ivp2_ta > 0 ) then
          tmp(1:3) & 
          = term_ta_lhs( wp3(kp1), wp3(k), wp2_zt(kp1), wp2_zt(k),  &
                         a1_zt(kp1), a1(k), a1_zt(k), gr%dzm(k), beta, k )
          zmscr05(k) = -tmp(3)
          zmscr06(k) = -tmp(2)
          zmscr07(k) = -tmp(1)
        endif

        if ( irtp2_ma + ithlp2_ma + irtpthlp_ma + & 
             iup2_ma + ivp2_ma > 0 ) then
          tmp(1:3) & 
          = term_ma_zm_lhs( wm_zm(k), gr%dzm(k), k )
          zmscr08(k) = -tmp(3)
          zmscr09(k) = -tmp(2)
          zmscr10(k) = -tmp(1)
        endif

      endif ! l_stats_samp

    enddo ! k=2..gr%nnzp-1


    ! Boundary Conditions
    ! These are set so that the sfc_var value of the variances and covariances can
    ! be used at the lowest boundary and the values of those variables can be set
    ! to 0 at the top boundary.  Fixed-point boundary conditions are used for both
    ! the variances and the covariances.
    lhs(:,1) = 0.0
    lhs(:,gr%nnzp) = 0.0

    lhs(k_mdiag,1) = 1.0
    lhs(k_mdiag,gr%nnzp) = 1.0

    ! This boundary condition was changed by dschanen on 24 April 2007
    ! When we run prognostically we want to preserve the surface value.
    !if ( l_iter ) then
    !  lhs(k_mdiag,1) = 1.0/dt
    !  lhs(k_mdiag,1) = 1.0/dt
    !endif

    return
  end subroutine xp2_xpyp_lhs

!===============================================================================
  subroutine xp2_xpyp_solve( solve_type, nrhs, rhs, lhs, xapxbp, err_code )

! Description:
!   Solve a tridiagonal system
!
! References:
!   None
!-------------------------------------------------------------------------------

    use lapack_wrap, only:  & 
      tridag_solve,  & ! Variable(s)
      tridag_solvex !, &
!    band_solve

    use grid_class, only: & 
      gr ! Variable(s)

    use stats_type, only: & 
      stat_update_var_pt  ! Procedure(s)

    use stats_variables, only: & 
      sfc, &  ! Derived type
      irtp2_matrix_condt_num, & ! Stat index Variables
      ithlp2_matrix_condt_num, & 
      irtpthlp_matrix_condt_num, & 
      iup2_vp2_matrix_condt_num, & 
      l_stats_samp  ! Logical

    implicit none

    ! External
    intrinsic :: trim

    ! Constant parameters
    integer, parameter :: & 
      kp1_mdiag = 1, & ! Momentum superdiagonal index.
      k_mdiag   = 2, & ! Momentum main diagonal index.
      km1_mdiag = 3    ! Momentum subdiagonal index.

    ! Input variables
    integer, intent(in) :: &
      nrhs  ! Number of right hand side vectors

    character(len=*), intent(in) ::  & 
      solve_type ! Variable(s) description

    ! Input/Ouput variables
    real, dimension(gr%nnzp,nrhs), intent(inout) :: & 
      rhs  ! Explicit contributions to x variance/covariance term [units vary]

    real, dimension(3,gr%nnzp), intent(inout) :: & 
      lhs  ! Implicit contributions to x variance/covariance term [units vary]

    real, dimension(gr%nnzp,nrhs), intent(out) ::  & 
      xapxbp ! Computed value of the variable(s) at <t+1> [units vary]

    integer, intent(out) :: & 
      err_code ! Returns an error code in the event of a singular matrix

    ! Local variables
    real :: rcond  ! Est. of the reciprocal of the condition # on the matrix

    integer ::  ixapxbp_matrix_condt_num ! Stat index

    ! --- Begin Code ---

    select case ( trim( solve_type ) )
    !-----------------------------------------------------------------------------
    ! Note that these are diagnostics from inverting the matrix, not a budget
    !-----------------------------------------------------------------------------
    case ( "rtp2" )
      ixapxbp_matrix_condt_num  = irtp2_matrix_condt_num

    case ( "thlp2" )
      ixapxbp_matrix_condt_num  = ithlp2_matrix_condt_num

    case ( "rtpthlp" )
      ixapxbp_matrix_condt_num  = irtpthlp_matrix_condt_num

    case ( "up2_vp2" )
      ixapxbp_matrix_condt_num  = iup2_vp2_matrix_condt_num

    case default
      ! No condition number is setup for the passive scalars
      ixapxbp_matrix_condt_num  = 0

    end select

    if ( l_stats_samp .and. ixapxbp_matrix_condt_num > 0 ) then
      call tridag_solvex & 
           ( solve_type, gr%nnzp, nrhs, &                                          ! Intent(in) 
             lhs(kp1_mdiag,:), lhs(k_mdiag,:), lhs(km1_mdiag,:), rhs(:,1:nrhs),  & ! Intent(inout)
             xapxbp(:,1:nrhs), rcond, err_code )                                   ! Intent(out)

      ! Est. of the condition number of the variance LHS matrix
      call stat_update_var_pt( ixapxbp_matrix_condt_num, 1, 1.0 / rcond, &  ! Intent(in)
                               sfc )                          ! Intent(inout)

    else
      call tridag_solve & 
           ( solve_type, gr%nnzp, nrhs, lhs(kp1_mdiag,:),  &        ! Intent(in)
             lhs(k_mdiag,:), lhs(km1_mdiag,:), rhs(:,1:nrhs),  &    ! Intent(inout)
             xapxbp(:,1:nrhs), err_code )                           ! Intent(out)
    end if

    return
  end subroutine xp2_xpyp_solve

!===============================================================================
  subroutine xp2_xpyp_implicit_stats( solve_type, xapxbp )

! Description:
!   Finalize implicit contributions for r_t'^2, th_l'^2, r_t'th_l',
!   u'^2, and v'^2.
!
! References:
!   None
!-------------------------------------------------------------------------------
    use grid_class, only: &
      gr ! Derived type variable

    use stats_type, only: & 
      stat_end_update_pt, & ! Procedure(s)
      stat_update_var_pt

    use stats_variables, only: & 
      zm,  & ! Variable(s) 
      irtp2_dp1, & 
      irtp2_dp2, & 
      irtp2_ta, & 
      irtp2_ma, & 
      ithlp2_dp1, & 
      ithlp2_dp2, & 
      ithlp2_ta, & 
      ithlp2_ma, & 
      irtpthlp_dp1, & 
      irtpthlp_dp2, & 
      irtpthlp_ta, & 
      irtpthlp_ma, & 
      iup2_dp1, & 
      iup2_dp2, & 
      iup2_ta, & 
      iup2_ma, & 
      iup2_pr1, & 
      ivp2_dp1, & 
      ivp2_dp2, & 
      ivp2_ta, & 
      ivp2_ma, & 
      ivp2_pr1, & 
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
      zmscr11

    implicit none

    ! External
    intrinsic :: max, min, trim

    ! Input variables
    character(len=*), intent(in) ::  & 
      solve_type ! Variable(s) description

    real, dimension(gr%nnzp), intent(in) ::  & 
      xapxbp ! Computed value of the variable at <t+1> [units vary]

    ! Local variables
    integer :: k, kp1, km1 ! Array indices

    ! Budget indices
    integer :: & 
      ixapxbp_dp1, & 
      ixapxbp_dp2, & 
      ixapxbp_ta, & 
      ixapxbp_ma, & 
      ixapxbp_pr1

    ! --- Begin Code ---

    select case ( trim( solve_type ) )
    case ( "rtp2" )
      ixapxbp_dp1 = irtp2_dp1
      ixapxbp_dp2 = irtp2_dp2
      ixapxbp_ta  = irtp2_ta
      ixapxbp_ma  = irtp2_ma
      ixapxbp_pr1 = 0

    case ( "thlp2" )
      ixapxbp_dp1 = ithlp2_dp1
      ixapxbp_dp2 = ithlp2_dp2
      ixapxbp_ta  = ithlp2_ta
      ixapxbp_ma  = ithlp2_ma
      ixapxbp_pr1 = 0

    case ( "rtpthlp" )
      ixapxbp_dp1 = irtpthlp_dp1
      ixapxbp_dp2 = irtpthlp_dp2
      ixapxbp_ta  = irtpthlp_ta
      ixapxbp_ma  = irtpthlp_ma
      ixapxbp_pr1 = 0

    case ( "up2" )
      ixapxbp_dp1 = iup2_dp1
      ixapxbp_dp2 = iup2_dp2
      ixapxbp_ta  = iup2_ta
      ixapxbp_ma  = iup2_ma
      ixapxbp_pr1 = iup2_pr1

    case ( "vp2" )
      ixapxbp_dp1 = ivp2_dp1
      ixapxbp_dp2 = ivp2_dp2
      ixapxbp_ta  = ivp2_ta
      ixapxbp_ma  = ivp2_ma
      ixapxbp_pr1 = ivp2_pr1

    case default ! No budgets are setup for the passive scalars
      ixapxbp_dp1 = 0
      ixapxbp_dp2 = 0
      ixapxbp_ta  = 0
      ixapxbp_ma  = 0
      ixapxbp_pr1 = 0

    end select

    do k = 2, gr%nnzp-1

      km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nnzp )

      ! x'y' term dp1 has both implicit and explicit components;
      ! call stat_end_update_pt.
      call stat_end_update_pt( ixapxbp_dp1, k, &            ! Intent(in)
                                 zmscr01(k) * xapxbp(k), &  ! Intent(in)
                               zm )                         ! Intent(inout)

      ! x'y' term dp2 is completely implicit; call stat_update_var_pt.
      call stat_update_var_pt( ixapxbp_dp2, k, &            ! Intent(in)
                                 zmscr02(k) * xapxbp(km1) & ! Intent(in)
                               + zmscr03(k) * xapxbp(k) & 
                               + zmscr04(k) * xapxbp(kp1), &
                               zm )                         ! Intent(inout)

      ! x'y' term ta has both implicit and explicit components;
      ! call stat_end_update_pt.
      call stat_end_update_pt( ixapxbp_ta, k, &              ! Intent(in)
                                 zmscr05(k) * xapxbp(km1) &  ! Intent(in)
                               + zmscr06(k) * xapxbp(k) &  
                               + zmscr07(k) * xapxbp(kp1), &
                               zm )                          ! Intent(inout)

      ! x'y' term ma is completely implicit; call stat_update_var_pt.
      call stat_update_var_pt( ixapxbp_ma, k, &              ! Intent(in)
                                 zmscr08(k) * xapxbp(km1) &  ! Intent(in)
                               + zmscr09(k) * xapxbp(k) & 
                               + zmscr10(k) * xapxbp(kp1), &
                               zm )                          ! Intent(inout)

      ! x'y' term pr1 has both implicit and explicit components;
      ! call stat_end_update_pt.
      call stat_end_update_pt( ixapxbp_pr1, k, &            ! Intent(in)
                                 zmscr11(k) * xapxbp(k), &  ! Intent(in)
                               zm )                         ! Intent(inout)

    end do ! k=2..gr%nnzp-1

    return
  end subroutine xp2_xpyp_implicit_stats

!===============================================================================
  subroutine xp2_xpyp_uv_rhs( solve_type, dt, l_iter, a1, a1_zt, & 
                              wp2, wp2_zt, wp3, wpthvp, tau_zm,  & 
                              xam, xbm, wpxap, wpxap_zt, wpxbp, & 
                              wpxbp_zt, xap2, xbp2, C4, C5, C14, & 
                              T0, beta, &
                              rhs )

! Description:
! Explicit contributions to u'^2 or v'^2
!-----------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use constants, only:  & 
        grav ! Variable(s)

    use stats_precision, only:  & 
        time_precision ! Variable(s)

    use stats_type, only: & 
        stat_begin_update_pt, & ! Procedure(s)
        stat_update_var_pt

    use stats_variables, only: & 
        ivp2_ta,  & ! Variable(s)
        ivp2_tp, & 
        ivp2_dp1, & 
        ivp2_pr1, & 
        ivp2_pr2, & 
        iup2_ta, & 
        iup2_tp, & 
        iup2_dp1, & 
        iup2_pr1, & 
        iup2_pr2, & 
        zm, & 
        zmscr01, & 
        zmscr11, & 
        l_stats_samp

    implicit none

! Input Variables
    character(len=*), intent(in) :: solve_type

    real(kind=time_precision), intent(in) :: & 
      dt          ! Model timestep                                 [s]

    logical, intent(in) :: & 
      l_iter  ! Whether x is prognostic (T/F)

    real, dimension(gr%nnzp), intent(in) :: & 
      a1,       & ! sigma_sqd_w-related term a_1 (momentum levels) [-]
      a1_zt,    & ! a_1 interpolated to thermodynamic levels       [-]
      wp2,      & ! w'^2 (momentum levels)                         [m^2/s^2]
      wp2_zt,   & ! w'^2 interpolated to thermodynamic levels      [m^2/s^2]
      wp3,      & ! w'^3 (thermodynamic levels)                    [m^3/s^3]
      wpthvp,   & ! w'th_v' (momentum levels)                      [K m/s]
      tau_zm,   & ! Time-scale tau on momentum levels              [s]
      xam,      & ! x_am (thermodynamic levels)                    [m/s]
      xbm,      & ! x_bm (thermodynamic levels)                    [m/s]
      wpxap,    & ! w'x_a' (momentum levels)                       [m^2/s^2]
      wpxap_zt, & ! w'x_a' interpolated to thermodynamic levels    [m^2/s^2]
      wpxbp,    & ! w'x_b' (momentum levels)                       [m^2/s^2]
      wpxbp_zt, & ! w'x_b' interpolated to thermodynamic levels    [m^2/s^2]
      xap2,     & ! x_a'^2 (momentum levels)                       [m^2/s^2]
      xbp2        ! x_b'^2 (momentum levels)                       [m^2/s^2]

    real, intent(in) :: & 
      C4,       & ! Model parameter C_4                            [-]
      C5,       & ! Model parameter C_5                            [-]
      C14,      & ! Model parameter C_{14}                         [-]
      T0,       & ! Reference temperature                          [K]
      beta        ! Model parameter beta                           [-]

    real, dimension(gr%nnzp,1), intent(out) :: & 
      rhs    ! Explicit contributions to x variance/covariance terms

    ! Local Variables

    ! Array indices
    integer :: k, kp1
    !integer :: km1

    real :: tmp

    integer :: & 
      ixapxbp_ta, & 
      ixapxbp_tp, & 
      ixapxbp_dp1, & 
      ixapxbp_pr1, & 
      ixapxbp_pr2

    !----------------------------- Begin Code ----------------------------------

    select case ( trim( solve_type ) )
    case ( "vp2" )
      ixapxbp_ta  = ivp2_ta
      ixapxbp_tp  = ivp2_tp
      ixapxbp_dp1 = ivp2_dp1
      ixapxbp_pr1 = ivp2_pr1
      ixapxbp_pr2 = ivp2_pr2
    case ( "up2" )
      ixapxbp_ta  = iup2_ta
      ixapxbp_tp  = iup2_tp
      ixapxbp_dp1 = iup2_dp1
      ixapxbp_pr1 = iup2_pr1
      ixapxbp_pr2 = iup2_pr2
    case default ! No budgets for passive scalars
      ixapxbp_ta  = 0
      ixapxbp_tp  = 0
      ixapxbp_dp1 = 0
      ixapxbp_pr1 = 0
      ixapxbp_pr2 = 0
    end select


    do k = 2, gr%nnzp-1, 1

      !km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nnzp )

      rhs(k,1) & 
      ! RHS turbulent advection (ta) term.
      = term_ta_rhs( wp3(kp1), wp3(k), wp2_zt(kp1), wp2_zt(k), &
                     a1_zt(kp1), a1(k), a1_zt(k), wpxbp_zt(kp1), wpxbp_zt(k), &
                     wpxap_zt(kp1), wpxap_zt(k), gr%dzm(k), beta ) &
      ! RHS turbulent production (tp) term.
      + (1.0 - C5)  & 
         * term_tp( xam(kp1), xam(k), xam(kp1), xam(k), & 
                    wpxap(k), wpxap(k), gr%dzm(k) ) & 
      ! RHS pressure term 1 (pr1) (and dissipation term 1 (dp1)).
      + term_pr1( C4, C14, xbp2(k), wp2(k), tau_zm(k) ) & 
      ! RHS pressure term 2 (pr2).
      + term_pr2( C5, grav, T0, wpthvp(k), wpxap(k), wpxbp(k), & 
                  xam(kp1), xam(k), xbm(kp1), xbm(k), gr%dzm(k) )

      ! RHS time tendency.
      if ( l_iter ) then
        rhs(k,1) = real( rhs(k,1) + 1.0/dt * xap2(k) )
      endif

      if ( l_stats_samp ) then

        ! Statistics: explicit contributions for up2 or vp2.

        ! x'y' term ta has both implicit and explicit components; call
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
        ! subtracts the value sent in, reverse the sign on term_ta_rhs.
        call stat_begin_update_pt( ixapxbp_ta, k, &                 ! Intent(in) 
        -term_ta_rhs( wp3(kp1), wp3(k), wp2_zt(kp1), wp2_zt(k), &   ! Intent(in)
                      a1_zt(kp1), a1(k), a1_zt(k), wpxbp_zt(kp1), wpxbp_zt(k), &
                      wpxap_zt(kp1), wpxap_zt(k), gr%dzm(k), beta ), &
                                   zm )                             ! Intent(inout)

        if ( ixapxbp_dp1 > 0 ) then
          ! Note:  The function term_pr1 is the explicit component of a
          !        semi-implicit solution to dp1 and pr1.
          ! Record the statistical contribution of the implicit component of
          ! term dp1 for up2 or vp2.  This will overwrite anything set
          ! statistically in xp2_xpyp_lhs for this term.
          ! Note:  To find the contribution of x'y' term dp1, substitute
          !        (2/3)*C_4 for the C_n input to function term_dp1_lhs.
          tmp = term_dp1_lhs( (2.0/3.0)*C4, tau_zm(k) )
          zmscr01(k) = -tmp
          ! Statistical contribution of the explicit component of term dp1 for
          ! up2 or vp2.
          ! x'y' term dp1 has both implicit and explicit components; call
          ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
          ! subtracts the value sent in, reverse the sign on term_pr1.
          ! Note:  To find the contribution of x'y' term dp1, substitute 0 for
          !        the C_14 input to function term_pr1.
          call stat_begin_update_pt( ixapxbp_dp1, k, &              ! Intent(in)
               -term_pr1( C4, 0.0, xbp2(k), wp2(k), tau_zm(k) ), &  ! Intent(in)
                                     zm )                           ! Intent(inout)
        endif

        if ( ixapxbp_pr1 > 0 ) then
          ! Note:  The function term_pr1 is the explicit component of a
          !        semi-implicit solution to dp1 and pr1.
          ! Statistical contribution of the implicit component of term pr1 for
          ! up2 or vp2.
          ! Note:  To find the contribution of x'y' term pr1, substitute
          !        (1/3)*C_14 for the C_n input to function term_dp1_lhs.
          tmp = term_dp1_lhs( (1.0/3.0)*C14, tau_zm(k) )
          zmscr11(k) = -tmp
          ! Statistical contribution of the explicit component of term pr1 for
          ! up2 or vp2.
          ! x'y' term pr1 has both implicit and explicit components; call
          ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
          ! subtracts the value sent in, reverse the sign on term_pr1.
          ! Note:  To find the contribution of x'y' term pr1, substitute 0 for
          !        the C_4 input to function term_pr1.
          call stat_begin_update_pt( ixapxbp_pr1, k, &                 ! Intent(in)  
               -term_pr1( 0.0, C14, xbp2(k), wp2(k), tau_zm(k) ), &    ! Intent(in)
                                     zm )                              ! Intent(inout)
        endif

        ! x'y' term pr2 is completely explicit; call stat_update_var_pt.
        call stat_update_var_pt( ixapxbp_pr2, k, &                          ! Intent(in)
              term_pr2( C5, grav, T0, wpthvp(k), wpxap(k), wpxbp(k),  &     ! Intent(in)
                        xam(kp1), xam(k), xbm(kp1), xbm(k), gr%dzm(k) ), &  
                                 zm )                                       ! Intent(inout)

        ! x'y' term tp is completely explicit; call stat_update_var_pt.
        call stat_update_var_pt( ixapxbp_tp, k, &                  ! Intent(in) 
              (1.0 - C5) &                                         ! Intent(in)
               * term_tp( xam(kp1), xam(k), xam(kp1), xam(k), &
                          wpxap(k), wpxap(k), gr%dzm(k) ), & 
                                 zm )                              ! Intent(inout)

      endif ! l_stats_samp

    enddo ! k=2..gr%nnzp-1


    ! Boundary Conditions
    ! These are set so that the sfc_var value of up2 or vp2 can be used at the
    ! lowest boundary and the values of those variables can be set to 0 at the top
    ! boundary.  Fixed-point boundary conditions are used for the variances.

    ! This boundary condition was changed by dschanen on 24 April 2007
    ! When we run prognostically we want to preserve the surface value.
    !if ( l_iter ) then
    !  rhs(1,1) = xap2(1) + 1.0/dt*xap2(1)
    !  rhs(gr%nnzp,1) = 1.0/dt*xap2(gr%nnzp)
    !else
       rhs(1,1) = xap2(1)
       rhs(gr%nnzp,1) = 0.0
    !end if

    return
  end subroutine xp2_xpyp_uv_rhs

!===============================================================================
  subroutine xp2_xpyp_rhs( solve_type, dt, l_iter, a1, a1_zt, &
                           wp2_zt, wp3, wpxap, wpxap_zt, & 
                           wpxbp, wpxbp_zt, xam, xbm, xapxbp, & 
                           Cn, tau_zm, threshold, beta, & 
                           rhs )

! Description:
!   Explicit contributions to r_t'^2, th_l'^2, r_t'th_l', sclr'r_t', sclr'th_l',
!   or sclr'^2.
!-------------------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use stats_precision, only:  & 
        time_precision ! Variable(s)

    use stats_type, only: & 
        stat_begin_update_pt, & ! Procedure(s)
        stat_update_var_pt

    use stats_variables, only: & 
        irtp2_ta,  & ! Variable(s)
        irtp2_tp,  & 
        irtp2_dp1, &
        ithlp2_ta,  & 
        ithlp2_tp,  &
        ithlp2_dp1, & 
        irtpthlp_ta,  & 
        irtpthlp_tp1, & 
        irtpthlp_tp2, &
        irtpthlp_dp1, & 
        zm, & 
        l_stats_samp

    implicit none

! Input Variables
    character(len=*), intent(in) :: solve_type

    real(kind=time_precision), intent(in) :: & 
      dt          ! Model timestep                                  [s]

    logical, intent(in) :: & 
      l_iter   ! Whether x is prognostic (T/F)

    real, dimension(gr%nnzp), intent(in) :: & 
      a1,       & ! sigma_sqd_w-related term a_1 (momentum levels)  [-]
      a1_zt,    & ! a_1 interpolated to thermodynamic levels        [-]
      wp2_zt,   & ! w'^2 interpolated to thermodynamic levels       [m^2/s^2]
      wp3,      & ! w'^3 (thermodynamic levels)                     [m^3/s^3]
      wpxap,    & ! w'x_a' (momentum levels)                        [m/s {x_am units}]
      wpxap_zt, & ! w'x_a' interpolated to thermodynamic levels     [m/s {x_am units}]
      wpxbp,    & ! w'x_b' (momentum levels)                        [m/s {x_bm units}]
      wpxbp_zt, & ! w'x_b' interpolated to thermodynamic levels     [m/s {x_bm units}]
      xam,      & ! x_am (thermodynamic levels)                     [{x_am units}]
      xbm,      & ! x_bm (thermodynamic levels)                     [{x_bm units}]
      xapxbp,   & ! x_a'x_b' (momentum levels)                      [{x_am units}*{x_bm units}]
      tau_zm,   & ! Time-scale tau on momentum levels               [s]
      Cn          ! Coefficient C_n                                 [-]

    real, intent(in) :: &
      threshold, & ! Smallest allowable magnitude value for x_a'x_b' [{x_am units}
               !                                                    *{x_bm units}] 
      beta         ! Model parameter beta                            [-]

    real, dimension(gr%nnzp,1), intent(out) :: & 
      rhs     ! Explicit contributions to x variance/covariance terms

    ! Local Variables

    ! Array indices
    integer :: k, kp1
    !integer :: km1

    integer :: & 
      ixapxbp_ta, & 
      ixapxbp_tp, & 
      ixapxbp_tp1, & 
      ixapxbp_tp2, &
      ixapxbp_dp1

    !------------------------------ Begin Code ---------------------------------

    select case ( trim( solve_type ) )
    case ( "rtp2" )
      ixapxbp_ta  = irtp2_ta
      ixapxbp_tp  = irtp2_tp
      ixapxbp_tp1 = 0
      ixapxbp_tp2 = 0
      ixapxbp_dp1 = irtp2_dp1
    case ( "thlp2" )
      ixapxbp_ta  = ithlp2_ta
      ixapxbp_tp  = ithlp2_tp
      ixapxbp_tp1 = 0
      ixapxbp_tp2 = 0
      ixapxbp_dp1 = ithlp2_dp1
    case ( "rtpthlp" )
      ixapxbp_ta  = irtpthlp_ta
      ixapxbp_tp  = 0
      ixapxbp_tp1 = irtpthlp_tp1
      ixapxbp_tp2 = irtpthlp_tp2
      ixapxbp_dp1 = irtpthlp_dp1
    case default ! No budgets for passive scalars
      ixapxbp_ta  = 0
      ixapxbp_tp  = 0
      ixapxbp_tp1 = 0
      ixapxbp_tp2 = 0
      ixapxbp_dp1 = 0
    end select


    do k = 2, gr%nnzp-1, 1

      !km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nnzp )

      rhs(k,1) & 
      ! RHS turbulent advection (ta) term.
      = term_ta_rhs( wp3(kp1), wp3(k), wp2_zt(kp1), wp2_zt(k), &
                     a1_zt(kp1), a1(k), a1_zt(k), wpxbp_zt(kp1), wpxbp_zt(k), &
                     wpxap_zt(kp1), wpxap_zt(k), gr%dzm(k), beta ) &
      ! RHS turbulent production (tp) term.
      + term_tp( xam(kp1), xam(k), xbm(kp1), xbm(k), & 
                 wpxbp(k), wpxap(k), gr%dzm(k) )

      ! RHS dissipation term 1 (dp1)
      rhs(k,1) &
      = rhs(k,1) + term_dp1_rhs( Cn(k), tau_zm(k), threshold )

      ! RHS time tendency.
      if ( l_iter ) then
        rhs(k,1) = real( rhs(k,1) + 1.0/dt * xapxbp(k) )
      endif

      if ( l_stats_samp ) then

        ! Statistics: explicit contributions for rtp2, thlp2, or rtpthlp.

        ! x'y' term ta has both implicit and explicit components; call
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
        ! subtracts the value sent in, reverse the sign on term_ta_rhs.
        call stat_begin_update_pt( ixapxbp_ta, k, &                ! Intent(in) 
        -term_ta_rhs( wp3(kp1), wp3(k), wp2_zt(kp1), wp2_zt(k), &  ! Intent(in)
                      a1_zt(kp1), a1(k), a1_zt(k), wpxbp_zt(kp1), wpxbp_zt(k), &
                      wpxap_zt(kp1), wpxap_zt(k), gr%dzm(k), beta ), &
                                   zm )                            ! Intent(inout)

        ! x'y' term dp1 has both implicit and explicit components; call
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
        ! subtracts the value sent in, reverse the sign on term_dp1_rhs.
        call stat_begin_update_pt( ixapxbp_dp1, k, &           ! Intent(in)
             -term_dp1_rhs( Cn(k), tau_zm(k), threshold ), &   ! Intent(in)
                                   zm )                        ! Intent(inout)

        ! rtp2/thlp2 case (1 turbulent production term)
        ! x'y' term tp is completely explicit; call stat_update_var_pt.
        call stat_update_var_pt( ixapxbp_tp, k, &             ! Intent(in)
              term_tp( xam(kp1), xam(k), xbm(kp1), xbm(k), &  ! Intent(in)
                       wpxbp(k), wpxap(k), gr%dzm(k) ), & 
                                 zm )                         ! Intent(inout)

        ! rtpthlp case (2 turbulent production terms)
        ! x'y' term tp1 is completely explicit; call stat_update_var_pt.
        ! Note:  To find the contribution of x'y' term tp1, substitute 0 for all
        !        the xam inputs and the wpxbp input to function term_tp.
        call stat_update_var_pt( ixapxbp_tp1, k, &    ! Intent(in)
              term_tp( 0.0, 0.0, xbm(kp1), xbm(k), &  ! Intent(in)
                       0.0, wpxap(k), gr%dzm(k) ), &
                                 zm )                 ! Intent(inout)

        ! x'y' term tp2 is completely explicit; call stat_update_var_pt.
        ! Note:  To find the contribution of x'y' term tp2, substitute 0 for all
        !        the xbm inputs and the wpxap input to function term_tp.
        call stat_update_var_pt( ixapxbp_tp2, k, &    ! Intent(in)
              term_tp( xam(kp1), xam(k), 0.0, 0.0, &  ! Intent(in)
                       wpxbp(k), 0.0, gr%dzm(k) ), &
                                 zm )                 ! Intent(inout)

      endif ! l_stats_samp

    enddo ! k=2..gr%nnzp-1


    ! Boundary Conditions
    ! These are set so that the sfc_var value of rtp2, thlp2, or rtpthlp (or sclrp2,
    ! sclrprtp, or sclrpthlp) can be used at the lowest boundary and the values of
    ! those variables can be set to 0 at the top boundary.  Fixed-point boundary
    ! conditions are used for both the variances and the covariances.

    ! This boundary condition was changed by dschanen on 24 April 2007
    ! When we run prognostically we want to preserve the surface value.
    !if ( .false. ) then
    !  rhs(1,1) = xapxbp(1) + 1.0/dt*xapxbp(1)
    !  rhs(gr%nnzp,1) = 1.0/dt*xapxbp(gr%nnzp)
    !else
       rhs(1,1) = xapxbp(1)
       rhs(gr%nnzp,1) = 0.0
    !endif

    return
  end subroutine xp2_xpyp_rhs

!===============================================================================
  pure function term_ta_lhs( wp3p1, wp3, wp2_ztp1, wp2_zt, &
                             a1_ztp1, a1, a1_zt, dzm, beta, level ) & 
  result( lhs )

! Description:
! Turbulent advection of x_a'x_b':  implicit portion of the code.
!
! The d(x_a'x_b')/dt equation contains a turbulent advection term:
!
! - d(w'x_a'x_b')/dz.
!
! A substitution is made in order to close the turbulent advection term,
! such that:
!
! w'x_a'x_b' = (1/3)*beta * a_1 * ( w'^3 / w'^2 ) * x_a'x_b'
!                 + (1-(1/3)*beta) * (a_1)^2 * ( w'^3 / (w'^2)^2 )
!                   * w'x_a' * w'x_b';
!
! where a_1 is a variable that is a function of sigma_sqd_w.  The turbulent
! advection term is rewritten as:
!
! - d [ (1/3)*beta * a_1 * ( w'^3 / w'^2 ) * x_a'x_b'
!          + (1-(1/3)*beta) * (a_1)^2 * ( w'^3 / (w'^2)^2 )
!            * w'x_a' * w'x_b' ]
!   / dz;
!
! which produces an implicit and an explicit portion of this term.  The implicit
! portion of this term is:
!
! - d [ (1/3)*beta * a_1 * ( w'^3 / w'^2 ) * x_a'x_b'(t+1) ] / dz.
!
! Since (1/3)*beta is a constant, it can be pulled outside of the derivative.
! The implicit portion of this term becomes:
!
! - (1/3)*beta * d [ a_1 * ( w'^3 / w'^2 ) * x_a'x_b'(t+1) ] / dz.
!
! Note:  When the term is brought over to the left-hand side, the sign is
!        reversed and the leading "-" in front of the term is changed to a "+".
!
! The timestep index (t+1) means that the value of x_a'x_b' being used is from
! the next timestep, which is being advanced to in solving the d(x_a'x_b')/dt
! equation.
!
! The implicit portion of this term is discretized as follows:
!
! The values of x_a'x_b' are found on the momentum levels, as are the values of
! w'^2 and a_1.  The values of w'^3 are found on the thermodynamic levels.  The
! variables x_a'x_b', w'^2, and a_1 are each interpolated to the intermediate
! thermodynamic levels.  The values of the mathematical expression (called F
! here) within the dF/dz term are computed on the thermodynamic levels.  Then
! the derivative (d/dz) of the expression (F) is taken over the central momentum
! level, yielding the desired result.  In this function, the values of F are as
! follows:
!
! F = a_1(t) * ( w'^3(t) / w'^2(t) ) * x_a'x_b'(t+1);
!
! where the timestep index (t) stands for the index of the current timestep.
!
!
! ==a1p1========wp2p1========xapxbpp1====================== m(k+1)
!
! ----a1(interp)---wp2(interp)---xapxbp(interp)---wp3p1---- t(k+1)
!
! ==a1==========wp2==========xapxbp=================dF/dz== m(k)
!
! ----a1(interp)---wp2(interp)---xapxbp(interp)---wp3------ t(k)
!
! ==a1m1========wp2m1========xapxbpm1====================== m(k-1)
!
! The vertical indices m(k+1), t(k+1), m(k), t(k), and m(k-1) correspond with
! altitudes zm(k+1), zt(k+1), zm(k), zt(k), and zm(k-1), respectively.  The
! letter "t" is used for thermodynamic levels and the letter "m" is used for
! momentum levels.
!
! dzm(k) = 1 / ( zt(k+1) - zt(k) )

! References:
!-----------------------------------------------------------------------

    use grid_class, only:  & ! gr%weights_zm2zt
        gr ! Variable(s)

    use constants, only:  &
        wtol_sqd

    use model_flags, only:  &
        l_standard_term_ta

    implicit none

    ! External
    intrinsic :: max

    ! Constant parameters
    integer, parameter :: & 
      kp1_mdiag = 1, & ! Momentum superdiagonal index.
      k_mdiag   = 2, & ! Momentum main diagonal index.
      km1_mdiag = 3    ! Momentum subdiagonal index.

    integer, parameter :: & 
      m_above = 1, & ! Index for upper momentum level grid weight.
      m_below = 2    ! Index for lower momentum level grid weight.

    ! Input Variables
    real, intent(in) :: & 
      wp3p1,    & ! w'^3(k+1)                                      [m^3/s^3]
      wp3,      & ! w'^3(k)                                        [m^3/s^3]
      wp2_ztp1, & ! w'^2 interpolated to thermodynamic level (k+1) [m^2/s^2]
      wp2_zt,   & ! w'^2 interpolated to thermodynamic level (k)   [m^2/s^2]
      a1_ztp1,  & ! a_1 interpolated to thermodynamic level (k+1)  [-]
      a1,       & ! a_1(k)                                         [-]
      a1_zt,    & ! a_1 interpolated to thermodynamic level (k)    [-]
      dzm,      & ! Inverse of grid spacing                        [1/m]
      beta        ! Model parameter                                [-]

    integer, intent(in) :: & 
      level ! Central momentum level (on which calculation occurs).

    ! Return Variable
    real, dimension(3) :: lhs

    ! Local Variables
    integer :: & 
      tkp1, & ! Thermodynamic level directly above central momentum level.
      tk      ! Thermodynamic level directly below central momentum level.


    ! Thermodynamic level (k+1) is between momentum level (k+1)
    ! and momentum level (k).
    tkp1 = level + 1

    ! Thermodynamic level (k) is between momentum level (k)
    ! and momentum level (k-1).
    tk = level

    if ( l_standard_term_ta ) then

      ! The turbulent advection term is discretized normally, in accordance
      ! with the model equations found in the documentation and the description
      ! listed above.

      ! Momentum superdiagonal: [ x xapxbp(k+1,<t+1>) ]
      lhs(kp1_mdiag)  &
      = + (1.0/3.0) * beta * dzm  &
          * a1_ztp1 * ( wp3p1 / max( wp2_ztp1, wtol_sqd ) )  &
          * gr%weights_zm2zt(m_above,tkp1)

      ! Momentum main diagonal: [ x xapxbp(k,<t+1>) ]
      lhs(k_mdiag)  &
      = + (1.0/3.0) * beta * dzm  &
          * (   a1_ztp1 * ( wp3p1 / max( wp2_ztp1, wtol_sqd ) )  &
                * gr%weights_zm2zt(m_below,tkp1)  &
              - a1_zt * ( wp3 / max( wp2_zt, wtol_sqd ) )  &
                * gr%weights_zm2zt(m_above,tk)  &
            )

      ! Momentum subdiagonal: [ x xapxbp(k-1,<t+1>) ]
      lhs(km1_mdiag)  &
      = - (1.0/3.0) * beta * dzm  &
          * a1_zt * ( wp3 / max( wp2_zt, wtol_sqd ) )  &
          * gr%weights_zm2zt(m_below,tk)

    else

      ! Brian tried a new discretization for the turbulent advection term, for
      ! which the implicit portion of the term is:
      ! - d [ a_1 * (1/3)*beta * ( w'^3 / w'^2 ) * x_a'x_b' ] / dz.  In order
      ! to help stabilize x_a'x_b', a_1 has been pulled outside the derivative.

      ! Momentum superdiagonal: [ x xapxbp(k+1,<t+1>) ]
      lhs(kp1_mdiag)  & 
      = + (1.0/3.0) * beta * a1 * dzm & 
          * ( wp3p1 / max( wp2_ztp1, wtol_sqd ) )  & 
          * gr%weights_zm2zt(m_above,tkp1)

      ! Momentum main diagonal: [ x xapxbp(k,<t+1>) ]
      lhs(k_mdiag)  & 
      = + (1.0/3.0) * beta * a1 * dzm & 
          * (   ( wp3p1 / max( wp2_ztp1, wtol_sqd ) )  & 
                * gr%weights_zm2zt(m_below,tkp1) & 
              - ( wp3 / max( wp2_zt, wtol_sqd ) ) & 
                * gr%weights_zm2zt(m_above,tk) & 
            )

      ! Momentum subdiagonal: [ x xapxbp(k-1,<t+1>) ]
      lhs(km1_mdiag)  & 
      = - (1.0/3.0) * beta * a1 * dzm & 
          * ( wp3 / max( wp2_zt, wtol_sqd ) ) & 
          * gr%weights_zm2zt(m_below,tk)

      ! End of Brian's a1 change.  14 Feb 2008.

    endif


    return
  end function term_ta_lhs

!===============================================================================
  pure function term_ta_rhs( wp3p1, wp3, wp2_ztp1, wp2_zt, &
                             a1_ztp1, a1, a1_zt, wpxbp_ztp1, wpxbp_zt, &
                             wpxap_ztp1, wpxap_zt, dzm, beta ) &
  result( rhs )

! Description:
! Turbulent advection of x_a'x_b':  explicit portion of the code.
!
! The d(x_a'x_b')/dt equation contains a turbulent advection term:
!
! - d(w'x_a'x_b')/dz.
!
! A substitution is made in order to close the turbulent advection term,
! such that:
!
! w'x_a'x_b' = (1/3)*beta * a_1 * ( w'^3 / w'^2 ) * x_a'x_b'
!                 + (1-(1/3)*beta) * (a_1)^2 * ( w'^3 / (w'^2)^2 )
!                   * w'x_a' * w'x_b';
!
! where a_1 is a variable that is a function of sigma_sqd_w.  The turbulent
! advection term is rewritten as:
!
! - d [ (1/3)*beta * a_1 * ( w'^3 / w'^2 ) * x_a'x_b'
!          + (1-(1/3)*beta) * (a_1)^2 * ( w'^3 / (w'^2)^2 )
!            * w'x_a' * w'x_b' ]
!   / dz;
!
! which produces an implicit and an explicit portion of this term.  The explicit
! portion of this term is:
!
! - d [ (1-(1/3)*beta) * (a_1)^2 * ( w'^3 / (w'^2)^2 )
!       * w'x_a' * w'x_b' ] / dz.
!
! Since (1-(1/3)*beta) is a constant, it can be pulled outside of the
! derivative.  The explicit portion of this term becomes:
!
! - (1-(1/3)*beta) * d [ (a_1)^2 * ( w'^3 / (w'^2)^2 )
!                        * w'x_a' * w'x_b' ] / dz.
!
! The explicit portion of this term is discretized as follows:
!
! The values of w'x_a', w'x_b', w'^2, and a_1 are found on the momentum levels.
! The values of w'^3 are found on the thermodynamic levels.  The variables
! w'x_a', w'x_b', w'^2, and a_1 are each interpolated to the intermediate
! thermodynamic levels.  The values of the mathematical expression (called F
! here) within the dF/dz term are computed on the thermodynamic levels.  Then
! the derivative (d/dz) of the expression (F) is taken over the central momentum
! level, yielding the desired result.  In this function, the values of F are as
! follows:
!
! F = ( a_1(t) )^2 * ( w'^3(t) / ( w'^2(t) )^2 )
!     * w'x_a'(t) * w'x_b'(t);
!
! where the timestep index (t) stands for the index of the current timestep.
!
!
! =a1p1=======wp2p1=======wpxapp1=======wpxbpp1============ m(k+1)
!
! -a1(interp)-wp2(interp)-wpxap(interp)-wpxbp(interp)-wp3p1 t(k+1)
!
! =a1=========wp2=========wpxap=========wpxbp====dF/dz===== m(k)
!
! -a1(interp)-wp2(interp)-wpxap(interp)-wpxbp(interp)-wp3-- t(k)
!
! =a1m1=======wp2m1=======wpxapm1=======wpxbpm1============ m(k-1)
!
! The vertical indices m(k+1), t(k+1), m(k), t(k), and m(k-1) correspond with
! altitudes zm(k+1), zt(k+1), zm(k), zt(k), and zm(k-1), respectively.  The
! letter "t" is used for thermodynamic levels and the letter "m" is used for
! momentum levels.
!
! dzm(k) = 1 / ( zt(k+1) - zt(k) )

! References:
!-----------------------------------------------------------------------

    use constants, only:  &
        wtol_sqd

    use model_flags, only:  &
        l_standard_term_ta

    implicit none

! External
    intrinsic :: max

! Input variables
    real, intent(in) :: & 
      wp3p1,      & ! w'^3(k+1)                                  [m^3/s^3]
      wp3,        & ! w'^3(k)                                    [m^3/s^3]
      wp2_ztp1,   & ! w'^2 interpolated to thermo. level (k+1)   [m^2/s^2]
      wp2_zt,     & ! w'^2 interpolated to thermo. level (k)     [m^2/s^2]
      a1_ztp1,    & ! a_1 interpolated to thermo. level (k+1)    [-]
      a1,         & ! a_1(k)                                     [-]
      a1_zt,      & ! a_1 interpolated to thermo. level (k)      [-]
      wpxbp_ztp1, & ! w'x_b' interpolated to thermo. level (k+1) [m/s {x_bm units}]
      wpxbp_zt,   & ! w'x_b' interpolated to thermo. level (k)   [m/s {x_bm units}]
      wpxap_ztp1, & ! w'x_a' interpolated to thermo. level (k+1) [m/s {x_am units}]
      wpxap_zt,   & ! w'x_a' interpolated to thermo. level (k)   [m/s {x_am units}]
      dzm,        & ! Inverse of grid spacing                    [1/m]
      beta          ! Model parameter                            [-]

! Return Variable
    real :: rhs


    if ( l_standard_term_ta ) then

      ! The turbulent advection term is discretized normally, in accordance
      ! with the model equations found in the documentation and the description
      ! listed above.

      rhs & 
      = - ( 1.0 - (1.0/3.0) * beta ) * dzm &
          * (   a1_ztp1**2 * wpxap_ztp1 * wpxbp_ztp1 &
                * ( wp3p1 / max( wp2_ztp1, wtol_sqd )**2 ) &
              - a1_zt**2 * wpxap_zt * wpxbp_zt &
                * ( wp3 / max( wp2_zt, wtol_sqd )**2 ) &
            )

    else

      ! Brian tried a new discretization for the turbulent advection term, for
      ! which the explicit portion of the term is:
      ! - d [ (a_1)^2 * (1-(1/3)*beta) * ( w'^3 / (w'^2)^2 )
      !       * w'x_a' * w'x_b' ] / dz.
      ! In order to help stabilize x_a'x_b', (a_1)^2 has been pulled outside
      ! the derivative.

      rhs & 
      = - ( 1.0 - (1.0/3.0) * beta ) * a1**2 * dzm & 
          * (   wpxap_ztp1 * wpxbp_ztp1 & 
                * ( wp3p1 / max( wp2_ztp1, wtol_sqd )**2 ) & 
              - wpxap_zt * wpxbp_zt & 
                * ( wp3 / max( wp2_zt, wtol_sqd )**2 ) & 
            )

      ! End of Brian's a1 change.  14 Feb 2008.

    endif


    return
  end function term_ta_rhs

!===============================================================================
  pure function term_tp( xamp1, xam, xbmp1, xbm,  & 
                         wpxbp, wpxap, dzm ) & 
  result( rhs )

! Description:
! Turbulent production of x_a'x_b':  explicit portion of the code.
!
! The d(x_a'x_b')/dt equation contains a turbulent production term:
!
! - w'x_b' d(x_am)/dz - w'x_a' d(x_bm)/dz.
!
! This term is solved for completely explicitly and is discretized as follows:
!
! The values of w'x_a' and w'x_b' are found on the momentum levels, whereas the
! values of x_am and x_bm are found on the thermodynamic levels.  The
! derivatives of both x_am and x_bm are taken over the intermediate (central)
! momentum level.  All of the remaining mathematical operations take place at
! the central momentum level, yielding the desired result.
!
! ---------xamp1------------xbmp1-------------------------- t(k+1)
!
! ===wpxap======d(xam)/dz=========d(xbm)/dz===wpxbp======== m(k)
!
! ---------xam--------------xbm---------------------------- t(k)
!
! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes zt(k+1),
! zm(k), and zt(k), respectively.  The letter "t" is used for thermodynamic
! levels and the letter "m" is used for momentum levels.
!
! dzm(k) = 1 / ( zt(k+1) - zt(k) )

! References:
!-----------------------------------------------------------------------

    implicit none

! Input variables
    real, intent(in) :: & 
      xam,   & ! x_am(k)                     [{x_am units}]
      xamp1, & ! x_am(k+1)                   [{x_am units}]
      xbm,   & ! x_bm(k)                     [{x_bm units}]
      xbmp1, & ! x_bm(k+1)                   [{x_bm units}]
      wpxbp, & ! w'x_b'(k)                   [m/s {x_bm units}]
      wpxap, & ! w'x_a'(k)                   [m/s {x_am units}]
      dzm      ! Inverse of grid spacing (k) [1/m]

! Return Variable
    real :: rhs

    rhs & 
    = - wpxbp * dzm * ( xamp1 - xam ) & 
      - wpxap * dzm * ( xbmp1 - xbm )

    return
  end function term_tp

!===============================================================================
  pure function term_dp1_lhs( Cn, tau_zm )  & 
  result( lhs )

! Description:
! Dissipation term 1 for x_a'x_b':  implicit portion of the code.
!
! The d(x_a'x_b')/dt equation contains dissipation term 1:
!
! - ( C_n / tau_zm ) x_a'x_b'.
!
! For cases where x_a'x_b' is a variance (in other words, where x_a and x_b are
! the same variable), the term is damped to a certain positive threshold,
! such that:
!
! - ( C_n / tau_zm ) * ( x_a'x_b' - threshold ).
!
! However, if x_a'x_b' is u'^2 or v'^2, the term is simply damped to 0.  The
! expression reverts to the form found in the first equation.  For u'^2 and
! v'^2, function 'term_dp1_lhs' is called, but function 'term_dp1_rhs' is not
! called.
!
! For cases where x_a'x_b' is a covariance (in other words, where x_a and x_b
! are different variables), threshold is set to 0, and the expression reverts to
! the form found in the first equation.
!
! This term is broken into implicit and explicit portions.  The equations for
! u'^2, v'^2, and any covariances only include the implicit portion.  The
! implicit portion of this term is:
!
! - ( C_n / tau_zm ) x_a'x_b'(t+1).
!
! Note:  When the implicit term is brought over to the left-hand side, the
!        sign is reversed and the leading "-" in front of the term is changed
!        to a "+".
!
! The timestep index (t+1) means that the value of x_a'x_b' being used is from
! the next timestep, which is being advanced to in solving the d(x_a'x_b')/dt
! equation.
!
! The values of x_a'x_b' are found on momentum levels.  The values of time-scale
! tau_zm are also found on momentum levels.
!
! Note:  For equations that use pressure term 1 (such as the equations for u'^2
!        and v'^2), C_n = ( 2*C_4 + C_14 ) / 3; which combines the implicit
!        contributions for dissipation term 1 and pressure term 1 into one
!        expression.  Otherwise, C_n = C_2.

! References:
!-----------------------------------------------------------------------

    implicit none

    real, intent(in) :: & 
      Cn,    & ! Coefficient C_n                       [-]
      tau_zm   ! Time-scale tau at momentum levels (k) [s]

    real :: lhs

! Momentum main diagonal: [ x xapxbp(k,<t+1>) ]
    lhs  & 
    = + Cn / tau_zm

    return
  end function term_dp1_lhs

!===============================================================================
  pure function term_dp1_rhs( Cn, tau_zm, threshold ) &
  result( rhs )

! Description:
! Dissipation term 1 for x_a'x_b':  explicit portion of the code.
!
! The d(x_a'x_b')/dt equation contains dissipation term 1:
!
! - ( C_n / tau_zm ) x_a'x_b'.
!
! For cases where x_a'x_b' is a variance (in other words, where x_a and x_b are
! the same variable), the term is damped to a certain positive threshold,
! such that:
!
! - ( C_n / tau_zm ) * ( x_a'x_b' - threshold ).
!
! However, if x_a'x_b' is u'^2 or v'^2, the term is simply damped to 0.  The
! expression reverts to the form found in the first equation.  For u'^2 and
! v'^2, function 'term_dp1_lhs' is called, but function 'term_dp1_rhs' is not
! called.
!
! For cases where x_a'x_b' is a covariance (in other words, where x_a and x_b
! are different variables), threshold is set to 0, and the expression reverts to
! the form found in the first equation.
!
! This term is broken into implicit and explicit portions.  The equations for
! u'^2, v'^2, and any covariances only include the implicit portion.  The
! explicit portion of this term is:
!
! + ( C_n / tau_zm ) * threshold.
!
! The values of time-scale tau_zm and the threshold are found on the momentum
! levels.
!
! Note:  The equations that use pressure term 1 (such as the equations for u'^2
!        and v'^2) do not call this function.  Thus, within this function,
!        C_n = C_2.

! References:
!-----------------------------------------------------------------------

    implicit none

    real, intent(in) :: &
      Cn,       & ! Coefficient C_n                               [-]
      tau_zm,   & ! Time-scale tau at momentum levels (k)         [s]
      threshold   ! Minimum allowable magnitude value of x_a'x_b' [units vary]

    real :: rhs

    rhs  & 
    = + ( Cn / tau_zm ) * threshold

    return
  end function term_dp1_rhs

!===============================================================================
  pure function term_pr1( C4, C14, xbp2, wp2, tau_zm ) & 
  result( rhs )

! Description:
! Pressure term 1 for x_a'x_b':  explicit portion of the code.
!
! Note:  Pressure term 1 is only used when x_a'x_b' is either u'^2 or v'^2.
!        For the following description, the equation for u'^2 is used as the
!        example.  The equation for v'^2 is the same as the equation for u'^2,
!        except that the v'^2 and u'^2 variables are switched.
!
! The d(u'^2)/dt equation contains dissipation term 1:
!
! - ( C_4 / tau_zm ) * ( u'^2 - (2/3)*em );
!
! and pressure term 1:
!
! - (2/3) * ( C_14 / tau_zm ) * em;
!
! where em = (1/2) * ( u'^2 + v'^2 + w'^2 ).
!
! This simplifies to:
!
! - [ ( 2*C_4 + C_14 ) / ( 3 * tau_zm ) ] * u'^2
!    + [ ( C_4 - C_14 ) / ( 3 * tau_zm ) ] * ( v'^2 + w'^2 ).
!
! The combined term has both implicit and explicit components.
! The implicit component is:
!
! - [ ( 2*C_4 + C_14 ) / ( 3 * tau_zm ) ] * u'^2(t+1).
!
! Note:  When the implicit term is brought over to the left-hand side, the
!        sign is reversed and the leading "-" in front of the term is changed
!        to a "+".
!
! Timestep index (t) stands for the index of the current timestep, while
! timestep index (t+1) stands for the index of the next timestep, which is being
! advanced to in solving the d(x_a'x_b')/dt equation.
!
! The implicit component of the combined dp1 and pr1 term is solved in function
! "term_dp1_lhs" above, where "( 2*C_4 + C_14 ) / 3" is sent in as "C_n".
!
! The explicit component of the combined dp1 and pr1 term is:
!
! + [ ( C_4 - C_14 ) / ( 3 * tau_zm ) ] * ( v'^2(t) + w'^2(t) );
!
! and is discretized as follows:
!
! The values for v'^2 and w'^2, as well as for tau_zm, are found on the momentum
! levels.  The mathematical operations all take place on the momentum levels,
! yielding the desired result.

! References:
!-----------------------------------------------------------------------

    implicit none

    real, intent(in) :: & 
      C4,    & ! Model parameter C_4                         [-]
      C14,   & ! Model parameter C_14                        [-]
      xbp2,  & ! v'^2(k) (if solving for u'^2) or vice versa [m^2/s^2]
      wp2,   & ! w'^2(k)                                     [m^2/s^2]
      tau_zm   ! Time-scale tau at momentum levels (k)       [s]

    real :: rhs

    rhs = + 1.0/3.0 * ( C4 - C14 ) * ( xbp2 + wp2 ) / tau_zm

    return
  end function term_pr1

!===============================================================================
  pure function term_pr2( C5, grav, T0, wpthvp, upwp, vpwp,  & 
                          ump1, um, vmp1, vm, dzm )  & 
  result( rhs )

! Description:
! Pressure term 2 for x_a'x_b':  explicit portion of the code.
!
! Note:  Pressure term 2 is only used when x_a'x_b' is either u'^2 or v'^2.
!        For the following description, the equation for u'^2 is used as the
!        example.  The equation for v'^2 is the same as the equation for u'^2.
!
! The d(u'^2)/dt equation contains pressure term 2:
!
! + (2/3) C_5 [ (g/th_0) w'th_v' - u'w' du/dz - v'w' dv/dz ].
!
! This term is solved for completely explicitly and is discretized as follows:
!
! The values of w'th_v', u'w', and v'w' are found on the momentum levels,
! whereas the values of um and vm are found on the thermodynamic levels.  The
! derivatives of both um and vm are taken over the intermediate (central)
! momentum level.  All the remaining mathematical operations take place at the
! central momentum level, yielding the desired result.
!
! -------ump1------------vmp1------------------------------ t(k+1)
!
! ==upwp======d(um)/dz========d(vm)/dz===vpwp=====wpthvp=== m(k)
!
! -------um--------------vm-------------------------------- t(k)
!
! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes zt(k+1),
! zm(k), and zt(k), respectively.  The letter "t" is used for thermodynamic
! levels and the letter "m" is used for momentum levels.
!
! dzm(k) = 1 / ( zt(k+1) - zt(k) )

! References:
!-----------------------------------------------------------------------

    implicit none

    real, intent(in) :: & 
      C5,     & ! Model parameter C_5             [-]
      grav,   & ! Gravitational acceleration      [m/s^2]
      T0,     & ! Reference temperature           [K]
      wpthvp, & ! w'th_v'(k)                      [m/K/s]
      upwp,   & ! u'w'(k)                         [m^2/s^2]
      vpwp,   & ! v'w'(k)                         [m^2/s^2]
      ump1,   & ! um(k+1)                         [m/s]
      um,     & ! um(k)                           [m/s]
      vmp1,   & ! vm(k+1)                         [m/s]
      vm,     & ! vm(k)                           [m/s]
      dzm       ! Inverse of the grid spacing (k) [1/m]

    real :: rhs

! As applied to w'2
    rhs = + (2.0/3.0) * C5 & 
            * ( ( grav / T0 ) * wpthvp & 
                -upwp * dzm * ( ump1 - um ) & 
                -vpwp * dzm * ( vmp1 - vm ) )

    return
  end function term_pr2

!===============================================================================
  subroutine pos_definite_variances( solve_type, dt, tolerance, & 
                                     xp2_np1 )
! Description:
!   Use the hole filling code to make a variance term positive definite
!--------------------------------------------------------------------------------
    use fill_holes, only: fill_holes_driver
    use grid_class, only: gr
    use stats_precision, only: time_precision

    use stats_variables, only:  & 
        zm, l_stats_samp, & 
        irtp2_pd, ithlp2_pd, iup2_pd, ivp2_pd ! variables
    use stats_type, only:  & 
        stat_begin_update, stat_end_update ! subroutines


    implicit none

    ! External
    intrinsic :: any, real, trim

    ! Input variables
    character(len=*), intent(in) :: & 
      solve_type
    real(kind=time_precision), intent(in) :: & 
      dt        ! Model timestep              [s]

    real, intent(in) :: & 
      tolerance ! Threshold for xp2_np1       [units vary]

    ! Input/Output variables
    real, intent(inout), dimension(gr%nnzp) ::  & 
      xp2_np1   ! Variance for <n+1>          [units vary]

    ! Local variables
    integer :: & 
      ixp2_pd

    select case( trim( solve_type ) )
    case ( "rtp2" )
      ixp2_pd = irtp2_pd
    case ( "thlp2" )
      ixp2_pd = ithlp2_pd
    case ( "up2" )
      ixp2_pd = iup2_pd
    case ( "vp2" )
      ixp2_pd = ivp2_pd
    case default
      ixp2_pd = 0 ! This includes the passive scalars
    end select

    if ( l_stats_samp ) then
      ! Store previous value for effect of the positive definite scheme
      call stat_begin_update( ixp2_pd, real( xp2_np1 / dt ), &    ! Intent(in)
                              zm )                                ! Intent(inout)
    endif


    if ( any( xp2_np1 < tolerance ) ) then

      ! Call the hole-filling scheme.
      ! The first pass-through should draw from only two levels on either side
      ! of the hole.
      call fill_holes_driver( 2, tolerance, "zm", & ! Intent(in) 
                              xp2_np1 )             ! Intent(inout)

    endif

    if ( l_stats_samp ) then
      ! Store previous value for effect of the positive definite scheme
      call stat_end_update( ixp2_pd, real( xp2_np1 / dt ), & ! Intent(in)
                            zm )                             ! Intent(inout)
    endif


    return
  end subroutine pos_definite_variances

!===============================================================================

end module advance_xp2_xpyp_module
