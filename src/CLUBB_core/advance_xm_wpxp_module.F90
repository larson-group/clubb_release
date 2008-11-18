!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module advance_xm_wpxp_module

  ! Description:
  ! Contains the CLUBB advance_xm_wpxp_module scheme.

  ! References:
  ! None
  !-----------------------------------------------------------------------

  implicit none

  private ! Default scope

  public  :: advance_xm_wpxp

  private :: xm_wpxp_lhs, & 
             xm_wpxp_rhs, & 
             xm_wpxp_solve, & 
             xm_term_ta_lhs, & 
             wpxp_term_ta_lhs, & 
             wpxp_term_tp_lhs, & 
             wpxp_terms_ac_pr2_lhs, & 
             wpxp_term_pr1_lhs, & 
             wpxp_terms_bp_pr3_rhs, &
             xm_correction_wpxp_cl

  ! Parameter Constants
  integer, parameter, private :: & 
    nsub = 2, & ! Number of subdiagonals in the LHS matrix
    nsup = 2    ! Number of superdiagonals in the LHS matrix


  contains

  !=============================================================================
  subroutine advance_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, wp3, & 
                              Kh_zt, tau_zm, Skw_zm, rtpthvp,  & 
                              rtm_forcing, thlpthvp,  & 
                              thlm_forcing, rtp2, thlp2, wp2_zt, & 
                              l_implemented, & 
                              sclrpthvp, sclrm_forcing, sclrp2,  & 
                              rtm, wprtp, thlm, wpthlp, & 
                              err_code, & 
                              sclrm, wpsclrp )

    ! Description:
    ! Advance the mean and flux terms by one timestep.

    ! References:
    ! Eqn. 16 & 17 on p. 3546 of
    ! ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !   Method and Model Description'' Golaz, et al. (2002)
    !   JAS, Vol. 59, pp. 3540--3551.

    ! See Also
    ! ``Equations for HOC'' Section 5:
    !   /Implicit solutions for the means and fluxes/
    !-----------------------------------------------------------------------

    use parameters_tunable, only:  & 
        C6rt,  & ! Variable(s)
        C6rtb,  & 
        C6rtc,  & 
        C6thl,  & 
        C6thlb,  & 
        C6thlc, & 
        C7,  & 
        C7b,  & 
        C7c,  & 
        c_K6,  & 
        c_Ksqd

    use constants, only:  & 
        fstderr, &  ! Constant
        max_mag_correlation, &
        zero_threshold

    use parameters_model, only: & 
        sclr_dim  ! Variable(s)

    use grid_class, only: & 
        gr,  & ! Variable(s)
        zm2zt ! Procedure(s)

    use model_flags, only: &
        l_3pt_sqd_dfsn,     & ! Variable(s)
        l_clip_semi_implicit

    use stats_precision, only:  & 
        time_precision ! Variable(s)

    use error_code, only:  & 
        lapack_error,  & ! Procedure(s)
        clubb_at_least_debug_level

    use stats_type, only: & 
        stat_begin_update, stat_end_update ! Procedure(s)

    use stats_variables, only: & 
        zt,  & ! Variable(s)
        irtm_cl,  & 
        ithlm_cl, & 
        irtm_matrix_condt_num, & 
        ithlm_matrix_condt_num, & 
        l_stats_samp


    implicit none

    ! External
    intrinsic :: exp, sqrt

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      dt               ! Timestep                                 [s]

    real, intent(in), dimension(gr%nnzp) :: & 
      sigma_sqd_w,   & ! sigma_sqd_w on momentum levels           [-]
      wm_zm,         & ! w wind component on momentum levels      [m/s]
      wm_zt,         & ! w wind component on thermodynamic levels [m/s]
      wp2,           & ! w'^2 (momentum levels)                   [m^2/s^2]
      wp2_zt,        & ! w'^2 (interpolated to thermo. levels)    [m^2/s^2]
      wp3,           & ! w'^3 (thermodynamic levels)              [m^3/s^3]
      Kh_zt,         & ! Eddy diffusivity on thermodynamic levels [m^2/s]
      tau_zm,        & ! Time-scale tau on momentum levels        [s]
      Skw_zm,        & ! Skewness of w on momentum levels         [-]
      rtpthvp,       & ! r_t'th_v' (momentum levels)              [(kg/kg) K]
      rtm_forcing,   & ! r_t forcing (thermodynamic levels)       [(kg/kg)/s]
      thlpthvp,      & ! th_l'th_v' (momentum levels)             [K^2]
      thlm_forcing,  & ! th_l forcing (thermodynamic levels)      [K/s]
      ! Added for clipping by Vince Larson 29 Sep 2007
      rtp2,          & ! r_t'^2 (momentum levels)                 [(kg/kg)^2]
      thlp2            ! th_l'^2 (momentum levels)                [K^2]
    ! End of Vince Larson's addition.

    logical, intent(in) ::  & 
      l_implemented      ! Flag for CLUBB being implemented in a larger model.


    ! Additional variables for passive scalars
    ! Input Variables
    real, intent(in), dimension(gr%nnzp,sclr_dim) ::  & 
      sclrpthvp, sclrm_forcing,  & !                           [Units vary]
      sclrp2                       ! For clipping Vince Larson [Units vary]

    ! Input/Output Variables
    real, intent(inout), dimension(gr%nnzp) ::  & 
      rtm,       & ! r_t  (total water mixing ratio)           [kg/kg]
      wprtp,     & ! w'r_t'                                    [(kg/kg) m/s]
      thlm,      & ! th_l (liquid water potential temperature) [K]
      wpthlp       ! w'th_l'                                   [K m/s]

    integer, intent(inout) :: err_code ! Model status

    ! Input/Output Variables
    real, intent(inout), dimension(gr%nnzp,sclr_dim) ::  & 
      sclrm, wpsclrp !                                     [Units vary]

    ! Local variables
    real, dimension(nsup+nsub+1,2*gr%nnzp) :: & 
      lhs  ! Implicit contributions to wpxp/xm (band diag. matrix) (LAPACK)

    real, dimension(gr%nnzp) ::  & 
      C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc

    ! Eddy Diffusion for wpthlp and wprtp.
    real, dimension(gr%nnzp) :: Kw6   ! wpxp eddy diff. [m^2/s]

    real, dimension(gr%nnzp) ::  & 
      a1 ! a_1 (momentum levels); See eqn. 24 in `Equations for HOC' [-]

    real, dimension(gr%nnzp) ::  & 
      a1_zt     ! a_1 interpolated to thermodynamic levels           [-]

    ! Variables used for adding (wpxp)^2: 3-point average
    ! diffusion coefficient.
    real, dimension(gr%nnzp) :: & 
      wprtp_zt, & 
      wpthlp_zt, & 
      wprtp_zt_sqd_3pt, & 
      wpthlp_zt_sqd_3pt, & 
      Kw6_rt, & 
      Kw6_thl

    ! Variables used for clipping of w'x' due to correlation
    ! of w with x, such that:
    ! corr_(w,x) = w'x' / [ sqrt(w'^2) * sqrt(x'^2) ];
    ! -1 <= corr_(w,x) <= 1.
    real, dimension(gr%nnzp) :: & 
      wpxp_upper_lim, & ! Keeps correlations from becoming greater than 1.
      wpxp_lower_lim    ! Keeps correlations from becoming less than -1.

    real, dimension(gr%nnzp) :: dummy_1d ! Unreferenced array

    real, allocatable, dimension(:,:) :: & 
      rhs,     &! Right-hand sides of band diag. matrix. (LAPACK)
      solution  ! solution vectors of band diag. matrix. (LAPACK)

    ! Constant parameters as a function of Skw.

    integer :: nrhs ! Number of RHS vectors

    real :: rcond

    ! Indices
    integer :: i
    integer :: k, km1, kp1

    !---------------------------------------------------------------------------

    ! ----- Begin Code -----
    if ( l_clip_semi_implicit .and. l_3pt_sqd_dfsn ) then
      nrhs = 1
    else
      nrhs = 2+sclr_dim
    end if

    ! Allocate rhs and solution vector
    allocate( rhs(2*gr%nnzp,nrhs) )
    allocate( solution(2*gr%nnzp,nrhs) )

    ! Compute C6 and C7 as a function of Skw
    ! The if...then is just here to save compute time
    if ( C6rt /= C6rtb ) then
      C6rt_Skw_fnc(1:gr%nnzp) = C6rtb + (C6rt-C6rtb) & 
        *EXP( -0.5 * (Skw_zm(1:gr%nnzp)/C6rtc)**2 )

    else
      C6rt_Skw_fnc(1:gr%nnzp) = C6rtb

    end if

    if ( C6thl /= C6thlb ) then
      C6thl_Skw_fnc(1:gr%nnzp) = C6thlb + (C6thl-C6thlb) & 
        *EXP( -0.5 * (Skw_zm(1:gr%nnzp)/C6thlc)**2 )

    else
      C6thl_Skw_fnc(1:gr%nnzp) = C6thlb

    end if

    if ( C7 /= C7b ) then
      C7_Skw_fnc(1:gr%nnzp) = C7b + (C7-C7b) & 
        *EXP( -0.5 * (Skw_zm(1:gr%nnzp)/C7c)**2 )

    else
      C7_Skw_fnc(1:gr%nnzp) = C7b

    end if

    !        C6rt_Skw_fnc = C6rt
    !        C6thl_Skw_fnc = C6thl
    !        C7_Skw_fnc = C7


    ! Define the Coefficent of Eddy Diffusivity for the wpthlp and wprtp.
    ! Kw6 is used for wpthlp and wprtp, which are located on momentum levels.
    ! Kw6 is located on thermodynamic levels.
    ! Kw6 = c_K6 * Kh_zt

    Kw6(1:gr%nnzp) = c_K6 * Kh_zt(1:gr%nnzp)


    ! (wpxp)^2: 3-point average diffusion coefficient.
    if ( l_3pt_sqd_dfsn ) then

      ! Interpolate w'x' (w'r_t' and w'th_l') from the momentum levels to the
      ! thermodynamic levels.  This is used for extra diffusion based on a
      ! three-point average of (w'x')^2.
      wprtp_zt  = zm2zt( wprtp )
      wpthlp_zt = zm2zt( wpthlp )

      do k = 1, gr%nnzp, 1

        km1 = max( k-1, 1 )
        kp1 = min( k+1, gr%nnzp )

        ! Compute the square of wprtp_zt, averaged over 3 points.  26 Jan 2008
        wprtp_zt_sqd_3pt(k) = ( wprtp_zt(km1)**2 + wprtp_zt(k)**2 & 
                                + wprtp_zt(kp1)**2 ) / 3.0
        ! Account for units of mix ratio (kg/kg)**2   Vince Larson 29 Jan 2008
        wprtp_zt_sqd_3pt(k) = 1e6 * wprtp_zt_sqd_3pt(k)

        ! Compute the square of wpthlp_zt, averaged over 3 points.
        ! 26 Jan 2008
        wpthlp_zt_sqd_3pt(k) = ( wpthlp_zt(km1)**2 + wpthlp_zt(k)**2 & 
                                 + wpthlp_zt(kp1)**2 ) / 3.0

      enddo

      ! Define Kw6_rt and Kw6_thl
      do k = 1, gr%nnzp, 1

        ! Kw6_rt must have units of m^2/s.  Since wprtp_zt_sqd_3pt has units
        ! of m/s (kg/kg), c_Ksqd is given units of m/(kg/kg) in this case.
        ! Vince Larson increased by c_Ksqd, 29Jan2008
        Kw6_rt(k)  = Kw6(k) + c_Ksqd * wprtp_zt_sqd_3pt(k)

        ! Kw6_thl must have units of m^2/s.  Since wpthlp_zt_sqd_3pt has units
        ! of m/s K, c_Ksqd is given units of m/K in this case.
        Kw6_thl(k) = Kw6(k) + c_Ksqd * wpthlp_zt_sqd_3pt(k)
        ! End Vince Larson's change

      enddo

    else  ! Three-point squared diffusion turned off.

      ! Define Kw6_rt and Kw6_thl
      do k = 1, gr%nnzp, 1

        Kw6_rt(k)  = Kw6(k)
        Kw6_thl(k) = Kw6(k)

      enddo

    endif  ! l_3pt_sqd_dfsn


    ! Define a_1 (located on momentum levels).
    ! It is a variable that is a function of sigma_sqd_w (where sigma_sqd_w is
    ! located on momentum levels).
    a1(1:gr%nnzp) = 1.0 / ( 1.0 - sigma_sqd_w(1:gr%nnzp) )

    ! Interpolate a_1 from momentum levels to thermodynamic levels.  This will
    ! be used for the w'x' turbulent advection (ta) term.
    a1_zt  = max( zm2zt( a1 ), zero_threshold )   ! Positive definite quantity

    ! Setup and decompose matrix for each variable.

    if ( l_clip_semi_implicit .and. l_3pt_sqd_dfsn ) then

      ! Compute the upper and lower limits of w'r_t' at every level,
      ! based on the correlation of w and r_t, such that:
      ! corr_(w,r_t) = w'r_t' / [ sqrt(w'^2) * sqrt(r_t'^2) ];
      ! -1 <= corr_(w,r_t) <= 1.
      if ( l_clip_semi_implicit ) then
        wpxp_upper_lim =  max_mag_correlation * sqrt( wp2 * rtp2 )
        wpxp_lower_lim = -wpxp_upper_lim

      end if

      ! Compute the implicit portion of the r_t and w'r_t' equations.
      ! Build the left-hand side matrix.
      call xm_wpxp_lhs( .true., dt, wprtp, a1_zt, wm_zm, wm_zt, wp2, wp2_zt, & ! Intent(in)
                        wp3, Kw6_rt, tau_zm, C7_Skw_fnc, C6rt_Skw_fnc, &       ! Intent(in)
                        wpxp_upper_lim, wpxp_lower_lim, l_implemented, &       ! Intent(in)
                        lhs )                                                  ! Intent(out)

      ! Compute the explicit portion of the r_t and w'r_t' equations.
      ! Build the right-hand side vector.
      call xm_wpxp_rhs( "rtm", .true., dt, rtm, wprtp, rtm_forcing, & ! Intent(in)
                       C7_Skw_fnc, rtpthvp, wp2_zt, a1_zt, wp3, & ! Intent(in)
                       wpxp_upper_lim, wpxp_lower_lim, &   ! Intent(in)
                       rhs(:,1) )                          ! Intent(out)

      ! Solve r_t / w'r_t'
      if ( l_stats_samp .and. irtm_matrix_condt_num > 0 ) then
        call xm_wpxp_solve( nrhs, &     ! Intent(in)
                            lhs, rhs, &  ! Intent(inout)
                            solution, err_code, rcond ) ! Intent(out)
      else
        call xm_wpxp_solve( nrhs, &     ! Intent(in)
                            lhs, rhs, &  ! Intent(inout)
                            solution, err_code ) ! Intent(out)
      end if

      if ( lapack_error( err_code ) )  then
        write(fstderr,'(a)') "rt LU decomp. failed"
        deallocate( rhs, solution )
        return
      end if

      call xm_wpxp_clipping_and_stats &
           ( "rtm", dt, wp2, rtp2, rtm, wprtp, solution(:,1), rcond )
      ! Clipping rtm
      ! Joshua Fasching March 2008


      ! Computed value before clipping
      if ( l_stats_samp ) then
        call stat_begin_update( irtm_cl, real( rtm / dt ), & ! Intent(in)
                                zt )                         ! Intent(inout)
      end if

      ! The arm_0003 case produces negative rtm near the tropopause.
      !    To avoid this, we clip rtm.  This is not a good solution,
      !    because it renders rtm non-conserved.  We should look into
      !    a positive definite advection scheme.
      !    Vince Larson.  13 Nov 2007

      ! The clipping of rtm causes a spurious source of moisture,
      !   particularly in SAM-CLUBB mode, and so this code will be
      !   disabled by default for now.
      !   David Schanen 15 Apr 2008
      do k = 1, gr%nnzp, 1
        if ( rtm(k) < 0.0 ) then
          !    rtm(k) = 0.0
          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) "rtm < 0 in advance_xm_wpxp_module at k= ", k
          endif
        endif

      enddo

      if ( l_stats_samp ) then
        call stat_end_update( irtm_cl, real( rtm / dt ), & ! Intent(in) 
                              zt )                         ! Intent(inout)
      endif


      ! Compute the upper and lower limits of w'th_l' at every level,
      ! based on the correlation of w and th_l, such that:
      ! corr_(w,th_l) = w'th_l' / [ sqrt(w'^2) * sqrt(th_l'^2) ];
      ! -1 <= corr_(w,th_l) <= 1.
      if ( l_clip_semi_implicit ) then
        wpxp_upper_lim =  max_mag_correlation * sqrt( wp2 * thlp2 )
        wpxp_lower_lim = -wpxp_upper_lim

      end if

      ! Compute the implicit portion of the th_l and w'th_l' equations.
      ! Build the left-hand side matrix.
      call xm_wpxp_lhs( .true., dt, wpthlp, a1_zt, wm_zm, wm_zt, wp2, wp2_zt, &      ! Intent(in)
                        wp3, Kw6_thl, tau_zm, C7_Skw_fnc, C6thl_Skw_fnc, & ! Intent(in)
                        wpxp_upper_lim, wpxp_lower_lim, l_implemented, &           ! Intent(in)
                        lhs )                                                      ! Intent(inout)

      ! Compute the explicit portion of the th_l and w'th_l' equations.
      ! Build the right-hand side vector.
      call xm_wpxp_rhs( "thlm", .true., dt, thlm, wpthlp, thlm_forcing, & ! Intent(in)
                        C7_Skw_fnc, thlpthvp, wp2_zt, a1_zt, wp3, &       ! Intent(in)
                        wpxp_upper_lim, wpxp_lower_lim, &                 ! Intent(in)
                        rhs(:,1) )                                        ! Intent(out)

      ! Solve for th_l / w'th_l'
      if ( l_stats_samp .and. ithlm_matrix_condt_num > 0 ) then
        call xm_wpxp_solve( nrhs, &     ! Intent(in)
                            lhs, rhs, &  ! Intent(inout)
                            solution, err_code, rcond ) ! Intent(out)
      else
        call xm_wpxp_solve( nrhs, &     ! Intent(in)
                            lhs, rhs, &  ! Intent(inout)
                            solution, err_code ) ! Intent(out)
      end if

      if ( lapack_error( err_code ) ) then
        write(fstderr,'(a)') "thetal LU decomp. failed"
        deallocate( rhs, solution )
        return
      endif

      call xm_wpxp_clipping_and_stats &
           ( "thlm", dt, wp2, thlp2, thlm, wpthlp, solution(:,1), rcond )
      ! Clipping thlm
      ! Joshua Fasching March 2008


      ! Computed value before clipping
      if ( l_stats_samp ) then
        call stat_begin_update( ithlm_cl, real(thlm / dt ), & ! Intent(in)
                                zt )                          ! Intent(inout)
      endif

      ! The value of potential temperature cannot fall below 0,
      ! so we clip accordingly
      do k = 1, gr%nnzp, 1
        if ( thlm(k) < 0.0 ) then
          thlm(k) = 0.0
          write(fstderr,*) "thlm < 0 in advance_xm_wpxp at k= ", k
        endif
      enddo


      if ( l_stats_samp ) then
        call stat_end_update( ithlm_cl, real( thlm/dt ), & ! Intent(in)
                              zt )                         ! Intent(inout)
      end if
      ! End change Joshua Fasching March 2008

      ! Solve sclrm / wpsclrp
      ! If sclr_dim is 0, then this loop will execute 0 times.
      do i = 1, sclr_dim, 1

        ! Compute the upper and lower limits of w'sclr' at every level,
        ! based on the correlation of w and sclr, such that:
        ! corr_(w,sclr) = w'sclr' / [ sqrt(w'^2) * sqrt(sclr'^2) ];
        ! -1 <= corr_(w,sclr) <= 1.
        if ( l_clip_semi_implicit ) then
          wpxp_upper_lim(:) =  max_mag_correlation * sqrt( wp2(:) * sclrp2(:,i) )
          wpxp_lower_lim(:) = -wpxp_upper_lim(:)

        end if

        ! Compute the implicit portion of the sclr and w'sclr' equations.
        ! Build the left-hand side matrix.
        call xm_wpxp_lhs( .true., dt, wpsclrp(:,i), a1_zt, wm_zm, wm_zt, wp2, & ! Intent(in)
                          wp2_zt, wp3, Kw6, tau_zm, C7_Skw_fnc, C6rt_Skw_fnc, & ! Intent(in)
                          wpxp_upper_lim, wpxp_lower_lim, l_implemented, &      ! Intent(in)
                          lhs )                                                 ! Intent(out)

        ! Compute the explicit portion of the sclrm and w'sclr' equations.
        ! Build the right-hand side vector.
        call xm_wpxp_rhs( "scalars", .true., dt, sclrm(:,i), wpsclrp(:,i),  &   ! Intent(in)
                          sclrm_forcing(:,i), C7_Skw_fnc, sclrpthvp(:,i), &     ! Intent(in)
                          wp2_zt, a1_zt, wp3, wpxp_upper_lim, wpxp_lower_lim, & ! Intent(in)
                          rhs(:,1) )                                               ! Intent(out)

        ! Solve for sclrm / w'sclr'
        call xm_wpxp_solve( nrhs, &     ! Intent(in)
                            lhs, rhs, &  ! Intent(inout)
                            solution, err_code ) ! Intent(out)

        if ( lapack_error( err_code ) ) then
          write(fstderr,'(a)') "Passive scalar ", i, " LU decomp. failed."
          deallocate( rhs, solution )
          return
        end if

        call xm_wpxp_clipping_and_stats &
             ( "scalars", dt, wp2, sclrp2(:,i), sclrm(:,i), wpsclrp(:,i), solution(:,1), rcond )

      end do ! passive scalars

    else ! Simple case, where l_clip_semi_implicit and l_3pt_sqd_dfsn are both false

      ! Create the lhs once
      call xm_wpxp_lhs( .true., dt, dummy_1d, a1_zt, wm_zm, wm_zt, wp2, wp2_zt, &  ! Intent(in)
                        wp3, Kw6, tau_zm, C7_Skw_fnc, C6rt_Skw_fnc, &  ! Intent(in)
                        dummy_1d, dummy_1d, l_implemented, &  ! Intent(in)
                        lhs ) ! Intent(out)

      ! Compute the explicit portion of the r_t and w'r_t' equations.
      ! Build the right-hand side vector.
      call xm_wpxp_rhs( "rtm", .true., dt, rtm, wprtp, rtm_forcing, &   ! Intent(in)
                       C7_Skw_fnc, rtpthvp, wp2_zt, a1_zt, wp3, & ! Intent(in)
                       wpxp_upper_lim, wpxp_lower_lim, &   ! Intent(in)
                       rhs(:,1) )                          ! Intent(out)

      ! Compute the explicit portion of the th_l and w'th_l' equations.
      ! Build the right-hand side vector.
      call xm_wpxp_rhs( "thlm", .true., dt, thlm, wpthlp, thlm_forcing, &   ! Intent(in)
                        C7_Skw_fnc, thlpthvp, wp2_zt, a1_zt, wp3, & ! Intent(in)
                        wpxp_upper_lim, wpxp_lower_lim, &      ! Intent(in)
                        rhs(:,2) )                             ! Intent(out)

      do i = 1, sclr_dim, 1
        call xm_wpxp_rhs( "scalars", .true., dt, sclrm(:,i), wpsclrp(:,i),  & ! Intent(in)
                           sclrm_forcing(:,i), C7_Skw_fnc, sclrpthvp(:,i), &   ! Intent(in)
                           wp2_zt, a1_zt, wp3, wpxp_upper_lim, wpxp_lower_lim, & ! Intent(in)
                           rhs(:,2+i) )                                        ! Intent(out)
      end do

      ! Solve for all fields
      if ( l_stats_samp .and. ithlm_matrix_condt_num + irtm_matrix_condt_num > 0 ) then
        call xm_wpxp_solve( nrhs, &     ! Intent(in)
                            lhs, rhs, &  ! Intent(inout)
                            solution, err_code, rcond ) ! Intent(out)
      else
        call xm_wpxp_solve( nrhs, &     ! Intent(in)
                            lhs, rhs, &  ! Intent(inout)
                            solution, err_code ) ! Intent(out)
      end if

      if ( lapack_error( err_code ) ) then
        write(fstderr,'(a)') "xm_wpxp LU decomp. failed."
        deallocate( rhs, solution )
        return
      end if

      call xm_wpxp_clipping_and_stats &
           ( "rtm", dt, wp2, rtp2, rtm, wprtp, solution(:,1), rcond )

      call xm_wpxp_clipping_and_stats &
           ( "thlm", dt, wp2, thlp2, thlm, wpthlp, solution(:,2), rcond )

      do i = 1, sclr_dim, 1
        call xm_wpxp_clipping_and_stats &
             ( "scalars", dt, wp2, sclrp2(:,i), sclrm(:,i), wpsclrp(:,i), solution(:,2+i), rcond )
      end do ! 1..sclr_dim

    end if ! l_clip_semi_implicit .and. l_3pt_sqd_dfsn

    deallocate( rhs, solution )

!       Error Report
!       (This code is unreachable)
!       Joshua Fasching Feb 2008

!        if ( lapack_error( err_code ) ) then

!           write(fstderr,*) "Error in advance_xm_wpxp"

!           write(fstderr,*) "Intent(in)"


!           write(fstderr,*) "dt = ", dt
!           write(fstderr,*) "tau_zm = ", tau_zm
!           write(fstderr,*) "wm_zm = ", wm_zm
!           write(fstderr,*) "wm_zt = ", wm_zt
!           write(fstderr,*) "wp2 = ", wp2
!           write(fstderr,*) "wp3 = ", wp3
!           write(fstderr,*) "sigma_sqd_w = ", sigma_sqd_w
!           write(fstderr,*) "Skw_zm = ", Skw_zm

    !write(fstderr,*) "dt = ", dt
    !write(fstderr,*) "sigma_sqd_w = ", sigma_sqd_w
    !write(fstderr,*) "wm_zm = ", wm_zm
    !write(fstderr,*) "wm_zt = ", wm_zt
    !write(fstderr,*) "wp2 = ", wp2
    !write(fstderr,*) "wp3 = ", wp3
    !write(fstderr,*) "Kh_zt = ", Kh_zt
    !write(fstderr,*) "tau_zm = ", tau_zm
    !write(fstderr,*) "Skw_zm = ", Skw_zm
    !write(fstderr,*) "rtpthvp = ", rtpthvp
    !write(fstderr,*) "rtm_forcing = ", rtm_forcing
    !write(fstderr,*) "thlpthvp = ", thlpthvp
    !write(fstderr,*) "thlm_forcing = ", thlm_forcing
    !write(fstderr,*) "rtp2 = ", rtp2
    !write(fstderr,*) "thlp2 = ", thlp2

!           if( present( sclrpthvp ) ) then
!              write(fstderr,*) "sclrpthvp = ", sclrpthvp
!           endif

!           if( present( sclrm_forcing ) ) then
!              write(fstderr,*) "sclrm_forcing = ", sclrm_forcing
!           endif

!           if( present( sclrp2 ) ) then
!              write(fstderr,*) "sclrp2 = ", sclrp2
!           endif

!           write(fstderr,*) "Intent(inout)"

!           write(fstderr,*) "rtm = ", rtm
!           write(fstderr,*) "wprtp = ", wprtp
!           write(fstderr,*) "thlm = ", thlm
!           write(fstderr,*) "wpthlp =", wpthlp

!           if( present( sclrm ) ) then
!              write(fstderr,*) "sclrm = ", sclrm
!           endif

!           if( present( wpsclrp ) ) then
!              write(fstderr,*) "wpsclrp = ", wpsclrp
!           endif

!        end if

    return

  end subroutine advance_xm_wpxp

  !=============================================================================
  subroutine xm_wpxp_lhs( l_iter, dt, wpxp, a1_zt, wm_zm, wm_zt, wp2, wp2_zt, &
                          wp3, Kw6, tau_zm, C7_Skw_fnc, C6x_Skw_fnc, & 
                          wpxp_upper_lim, wpxp_lower_lim, l_implemented, &
                          lhs )

    ! Description:
    ! Compute LHS band diagonal matrix for xm and w'x'.
    ! This subroutine computes the implicit portion of
    ! the xm and w'x' equations.

    ! References:
    !------------------------------------------------------------------------

    use parameters_tunable, only:  & 
        nu6 ! Variable(s)

    use grid_class, only:  & 
        gr,  & ! Variable(s)
        zm2zt ! Procedure(s)

    use constants, only: &
        gamma_over_implicit_ts ! Variable(s)

    use model_flags, only: &
        l_clip_semi_implicit ! Variable(s)

    use stats_precision, only:  & 
        time_precision ! Variable(s)

    use diffusion, only:  & 
        diffusion_zm_lhs ! Procedure(s)

    use mean_adv, only: & 
        term_ma_zt_lhs,  & ! Procedure(s)
        term_ma_zm_lhs

    use clip_semi_implicit, only: & 
        clip_semi_imp_lhs ! Procedure(s)

    use stats_variables, only: & 
        ztscr01,  & ! Variable(s)
        ztscr02, & 
        ztscr03, & 
        ztscr04, & 
        ztscr05, & 
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
        zmscr13, & 
        zmscr14, & 
        zmscr15, & 
        l_stats_samp, & 
        ithlm_ma, & 
        ithlm_ta, & 
        irtm_ma, & 
        irtm_ta, & 
        iwpthlp_ma, & 
        iwpthlp_ta, & 
        iwpthlp_tp, & 
        iwpthlp_ac, & 
        iwpthlp_pr1, & 
        iwpthlp_pr2, & 
        iwpthlp_dp1, & 
        iwpthlp_sicl, & 
        iwprtp_ma, & 
        iwprtp_ta, & 
        iwprtp_tp, & 
        iwprtp_ac, & 
        iwprtp_pr1, & 
        iwprtp_pr2, & 
        iwprtp_dp1, & 
        iwprtp_sicl


    implicit none

    ! External
    intrinsic :: min, max

    ! Constant parameters
    ! Left-hand side matrix diagonal identifiers for
    ! momentum-level variable, w'x'.
    integer, parameter ::  &
      m_kp1_mdiag = 1, & ! Momentum superdiagonal index for w'x'.
      m_kp1_tdiag = 2, & ! Thermodynamic superdiagonal index for w'x'.
      m_k_mdiag   = 3, & ! Momentum main diagonal index for w'x'.
      m_k_tdiag   = 4, & ! Thermodynamic subdiagonal index for w'x'.
      m_km1_mdiag = 5    ! Momentum subdiagonal index for w'x'.

    ! Left-hand side matrix diagonal identifiers for
    ! thermodynamic-level variable, xm.
    integer, parameter ::  &
      t_kp1_tdiag = 1, & ! Thermodynamic superdiagonal index for xm.
      t_k_mdiag   = 2, & ! Momentum superdiagonal index for xm.
      t_k_tdiag   = 3, & ! Thermodynamic main diagonal index for xm.
      t_km1_mdiag = 4, & ! Momentum subdiagonal index for xm.
      t_km1_tdiag = 5    ! Thermodynamic subdiagonal index for xm.

    ! Input variables
    logical, intent(in) :: l_iter

    real(kind=time_precision), intent(in) ::  & 
      dt                 ! Timestep                                 [s]

    real, intent(in), dimension(gr%nnzp) :: & 
      wpxp,            & ! w'x' (momentum levels) at timestep (t)   [{xm units} m/s]
      a1_zt,           & ! a_1 interpolated to thermodynamic levels  [-]
      wm_zm,           & ! w wind component on momentum levels      [m/s]
      wm_zt,           & ! w wind component on thermodynamic levels [m/s]
      wp2,             & ! w'^2 (momentum levels)                   [m^2/s^2]
      wp2_zt,          & ! w'^2 interpolated to thermodynamic levels[m^2/s^2]
      wp3,             & ! w'^3 (thermodynamic levels)              [m^3/s^3]
      Kw6,             & ! Coefficient of eddy diffusivity for w'x' [m^2/s]
      tau_zm,          & ! Time-scale tau on momentum levels        [s]
      C7_Skw_fnc,      & ! C_7 parameter with Sk_w applied          [-]
      C6x_Skw_fnc,     & ! C_6x parameter with Sk_w applied         [-]
      wpxp_upper_lim,  & ! Keeps correlations from becoming > 1.    [units vary]
      wpxp_lower_lim     ! Keeps correlations from becoming < -1.   [units vary]

    logical, intent(in) ::  & 
      l_implemented ! Flag for CLUBB being implemented in a larger model.

    ! Output Variable
    real, intent(out), dimension(nsup+nsub+1,2*gr%nnzp) ::  & 
      lhs ! Implicit contributions to wpxp/xm (band diag. matrix) (LAPACK)

    ! Local Variables

    ! Indices
    !integer :: km1
    integer :: k, kp1
    integer :: k_xm, k_wpxp

    real, dimension(3) :: tmp


    ! Initialize the left-hand side matrix to 0.
    lhs = 0.0

    do k = 2, gr%nnzp-1, 1

      ! Define indices

      !km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nnzp )

      k_xm   = 2*k - 1
      k_wpxp = 2*k


      !!!!!***** xm *****!!!!!

      ! xm: Left-hand side (implicit xm portion of the code).
      !
      ! Thermodynamic subdiagonal (lhs index: t_km1_tdiag)
      !         [ x xm(k-1,<t+1>) ]
      ! Momentum subdiagonal (lhs index: t_km1_mdiag)
      !         [ x wpxp(k-1,<t+1>) ]
      ! Thermodynamic main diagonal (lhs index: t_k_tdiag)
      !         [ x xm(k,<t+1>) ]
      ! Momentum superdiagonal (lhs index: t_k_mdiag)
      !         [ x wpxp(k,<t+1>) ]
      ! Thermodynamic superdiagonal (lhs index: t_kp1_tdiag)
      !         [ x xm(k+1,<t+1>) ]

      ! LHS mean advection (ma) term.
      if ( .not. l_implemented ) then

        lhs((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/),k_xm) & 
        = lhs((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/),k_xm) & 
        + term_ma_zt_lhs( wm_zt(k), gr%dzt(k), k )

      else

        lhs((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/),k_xm) & 
        = lhs((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/),k_xm) + 0.0

      endif

      ! LHS turbulent advection (ta) term.
      lhs((/t_k_mdiag,t_km1_mdiag/),k_xm) & 
      = lhs((/t_k_mdiag,t_km1_mdiag/),k_xm) & 
      + xm_term_ta_lhs( gr%dzt(k) )

      ! LHS time tendency.
      lhs(t_k_tdiag,k_xm) & 
      = real( lhs(t_k_tdiag,k_xm) + 1.0 / dt )

      if (l_stats_samp) then

        ! Statistics: implicit contributions for rtm or thlm.

        if ( irtm_ma > 0 .or. ithlm_ma > 0 ) then
          if ( .not. l_implemented ) then
            tmp(1:3) =  & 
            + term_ma_zt_lhs( wm_zt(k), gr%dzt(k), k )
            ztscr01(k) = - tmp(3)
            ztscr02(k) = - tmp(2)
            ztscr03(k) = - tmp(1)
          else
            ztscr01(k) = 0.0
            ztscr02(k) = 0.0
            ztscr03(k) = 0.0
          endif
        endif

        if ( irtm_ta > 0 .or. ithlm_ta > 0 ) then
          tmp(1:2) = & 
          + xm_term_ta_lhs( gr%dzt(k) )
          ztscr04(k) = - tmp(2)
          ztscr05(k) = - tmp(1)
        endif

      endif


      !!!!!***** w'x' *****!!!!!

      ! w'x': Left-hand side (implicit w'x' portion of the code).
      !
      ! Momentum subdiagonal (lhs index: m_km1_mdiag)
      !         [ x wpxp(k-1,<t+1>) ]
      ! Thermodynamic subdiagonal (lhs index: m_k_tdiag)
      !         [ x xm(k,<t+1>) ]
      ! Momentum main diagonal (lhs index: m_k_mdiag)
      !         [ x wpxp(k,<t+1>) ]
      ! Thermodynamic superdiagonal (lhs index: m_kp1_tdiag)
      !         [ x xm(k+1,<t+1>) ]
      ! Momentum superdiagonal (lhs index: m_kp1_mdiag)
      !         [ x wpxp(k+1,<t+1>) ]

      ! LHS mean advection (ma) term.
      lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wpxp) & 
      = lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wpxp) & 
      + term_ma_zm_lhs( wm_zm(k), gr%dzm(k), k )

      ! LHS turbulent advection (ta) term.
      ! Note:  An "over-implicit" weighted time step is applied to this term.
      !        The weight of the implicit portion of this term is controlled
      !        by the factor gamma_over_implicit_ts (abbreviated "gamma" in the
      !        the equation in order to balance a weight that is not equal to 1,
      !        such that:
      !             -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
      !        where X is the variable that is being solved for in a predictive
      !        equation (w'x' in this case), y(t) is the linearized portion of
      !        the term that gets treated implicitly, and RHS is the portion of
      !        the term that is always treated explicitly (in the case of the
      !        w'x' turbulent advection term, RHS = 0).  A weight of greater
      !        than 1 can be applied to make the term more numerically stable.
      lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wpxp)  & 
      = lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wpxp)  &
      + gamma_over_implicit_ts  &
      * wpxp_term_ta_lhs( wp2_zt(kp1), wp2_zt(k),  & 
                          a1_zt(kp1), a1_zt(k),  & 
                          wp3(kp1), wp3(k), gr%dzm(k), k )

      ! LHS turbulent production (tp) term.
      lhs((/m_kp1_tdiag,m_k_tdiag/),k_wpxp) & 
      = lhs((/m_kp1_tdiag,m_k_tdiag/),k_wpxp) & 
      + wpxp_term_tp_lhs( wp2(k), gr%dzm(k) )

      ! LHS accumulation (ac) term and pressure term 2 (pr2).
      lhs(m_k_mdiag,k_wpxp) & 
      = lhs(m_k_mdiag,k_wpxp) & 
      + wpxp_terms_ac_pr2_lhs( C7_Skw_fnc(k),  & 
                               wm_zt(kp1), wm_zt(k), gr%dzm(k) )

      ! LHS pressure term 1 (pr1).
      lhs(m_k_mdiag,k_wpxp) & 
      = lhs(m_k_mdiag,k_wpxp) & 
      + wpxp_term_pr1_lhs( C6x_Skw_fnc(k), tau_zm(k) )

      ! LHS eddy diffusion term: dissipation term 1 (dp1).
      lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wpxp) & 
      = lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wpxp) & 
      + diffusion_zm_lhs( Kw6(k), Kw6(kp1), nu6, & 
                          gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )

      ! LHS time tendency.
      if (l_iter) lhs(m_k_mdiag,k_wpxp)  &
                  = real( lhs(m_k_mdiag,k_wpxp) + 1.0 / dt )

      ! LHS portion of semi-implicit clipping term.
      if ( l_clip_semi_implicit ) then

        lhs(m_k_mdiag,k_wpxp) & 
        = lhs(m_k_mdiag,k_wpxp) & 
        + clip_semi_imp_lhs( dt, wpxp(k),  & 
                             .true., wpxp_upper_lim(k),  & 
                             .true., wpxp_lower_lim(k) )

      endif

      if (l_stats_samp) then

        ! Statistics: implicit contributions for wprtp or wpthlp.

        if ( iwprtp_ma > 0 .or. iwpthlp_ma > 0 ) then
          tmp(1:3) = & 
          + term_ma_zm_lhs( wm_zm(k), gr%dzm(k), k )
          zmscr01(k) = - tmp(3)
          zmscr02(k) = - tmp(2)
          zmscr03(k) = - tmp(1)
        endif

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note above for LHS turbulent
        !        advection (ta) term).
        if ( iwprtp_ta > 0 .or. iwpthlp_ta > 0 ) then
          tmp(1:3)  &
          = gamma_over_implicit_ts  &
          * wpxp_term_ta_lhs( wp2_zt(kp1), wp2_zt(k),  &
                              a1_zt(kp1), a1_zt(k),  &
                              wp3(kp1), wp3(k), gr%dzm(k), k )
          zmscr04(k) = - tmp(3)
          zmscr05(k) = - tmp(2)
          zmscr06(k) = - tmp(1)
        endif

        if ( iwprtp_tp > 0 .or. iwpthlp_tp > 0 ) then
          tmp(1:2) = & 
          + wpxp_term_tp_lhs( wp2(k), gr%dzm(k) )
          zmscr07(k) = - tmp(2)
          zmscr08(k) = - tmp(1)
        endif

        ! Note:  To find the contribution of w'x' term ac, substitute 0 for the
        !        C_7 skewness function input to function wpxp_terms_ac_pr2_lhs.
        if ( iwprtp_ac > 0 .or. iwpthlp_ac > 0 ) then
          zmscr09(k) =  & 
          - wpxp_terms_ac_pr2_lhs( 0.0, & 
                                   wm_zt(kp1), wm_zt(k), gr%dzm(k) )
        endif

        if ( iwprtp_pr1 > 0 .or. iwpthlp_pr1 > 0 ) then
          zmscr10(k) = & 
          - wpxp_term_pr1_lhs( C6x_Skw_fnc(k), tau_zm(k) )
        endif

        ! Note:  To find the contribution of w'x' term pr2, add 1 to the
        !        C_7 skewness function input to function wpxp_terms_ac_pr2_lhs.
        if ( iwprtp_pr2 > 0 .or. iwpthlp_pr2 > 0 ) then
          zmscr11(k) = & 
          - wpxp_terms_ac_pr2_lhs( (1.0+C7_Skw_fnc(k)), & 
                                   wm_zt(kp1), wm_zt(k), gr%dzm(k) )
        endif

        if ( iwprtp_dp1 > 0 .or. iwpthlp_dp1 > 0 ) then
          tmp(1:3) = & 
          + diffusion_zm_lhs( Kw6(k), Kw6(kp1), nu6, & 
                              gr%dzt(kp1), gr%dzt(k), gr%dzm(k), k )
          zmscr12(k) = - tmp(3)
          zmscr13(k) = - tmp(2)
          zmscr14(k) = - tmp(1)
        endif

        if ( l_clip_semi_implicit ) then
          if ( iwprtp_sicl > 0 .or. iwpthlp_sicl > 0 ) then
            zmscr15(k) = & 
            - clip_semi_imp_lhs( dt, wpxp(k),  & 
                                 .true., wpxp_upper_lim(k),  & 
                                 .true., wpxp_lower_lim(k) )
          endif
        endif

      endif


    enddo ! 2..gr%nnzp-1


    ! Boundary conditions

    ! Both the mean (xm) and the turbulent flux (wpxp) use fixed-point
    ! boundary conditions.  Therefore, anything set in the above loop
    ! at both the upper and lower boundaries would be overwritten here.
    ! However, the above loop does not extend to the boundary levels.
    ! An array with a value of 1 at the main diagonal on the left-hand
    ! side and with values of 0 at all other diagonals on the left-hand
    ! side will preserve the right-hand side value at that level.
    !
    !   xm(1)  wpxp(1) ... xm(nz) wpxp(nz)
    ! [  0.0     0.0         0.0     0.0  ]
    ! [  0.0     0.0         0.0     0.0  ]
    ! [  1.0     1.0   ...   1.0     1.0  ]
    ! [  0.0     0.0         0.0     0.0  ]
    ! [  0.0     0.0         0.0     0.0  ]

    ! Lower boundary
    k      = 1
    k_xm   = 2*k - 1
    k_wpxp = 2*k

    ! xm
    lhs(:,k_xm)           = 0.0
    lhs(t_k_tdiag,k_xm)   = 1.0
    ! w'x'
    lhs(:,k_wpxp)         = 0.0
    lhs(m_k_mdiag,k_wpxp) = 1.0

    ! Upper boundary
    k      = gr%nnzp
    k_xm   = 2*k - 1
    k_wpxp = 2*k

    ! xm
    lhs(:,k_xm)           = 0.0
    lhs(t_k_tdiag,k_xm)   = 1.0
    ! w'x'
    lhs(:,k_wpxp)         = 0.0
    lhs(m_k_mdiag,k_wpxp) = 1.0


    return
  end subroutine xm_wpxp_lhs

  !=============================================================================
  subroutine xm_wpxp_rhs( solve_type, l_iter, dt, xm, wpxp, xm_forcing, & 
                          C7_Skw_fnc, xpthvp, wp2_zt, a1_zt, wp3, & 
                          wpxp_upper_lim, wpxp_lower_lim, &
                          rhs )
    ! Description:
    ! Compute RHS vector for xm and w'x'.
    ! This subroutine computes the explicit portion of
    ! the xm and w'x' equations.

    ! References:
    !------------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use constants, only:  &
        gamma_over_implicit_ts ! Variable(s)

    use model_flags, only: &
        l_clip_semi_implicit ! Variable(s)

    use stats_precision, only:  & 
        time_precision ! Variable(s)

    use clip_semi_implicit, only: & 
        clip_semi_imp_rhs ! Procedure(s)

    use stats_type, only: & 
        stat_update_var_pt, & 
        stat_begin_update_pt

    use stats_variables, only: & 
        zt, & ! Variable(s)
        zm, & 
        irtm_forcing, & 
        ithlm_forcing, & 
        iwprtp_bp, & 
        iwprtp_pr3, & 
        iwprtp_sicl, &
        iwprtp_ta, &
        iwpthlp_bp, & 
        iwpthlp_pr3, & 
        iwpthlp_sicl, &
        iwpthlp_ta, &
        l_stats_samp


    implicit none

    ! Input Variables
    character(len=*), intent(in) :: & 
      solve_type  ! Variables being solved for.

    logical, intent(in) :: l_iter

    real(kind=time_precision), intent(in) ::  & 
      dt                 ! Timestep                               [s]

    real, dimension(gr%nnzp), intent(in) :: & 
      xm,              & ! xm (thermodynamic levels)              [{xm units}]
      wpxp,            & ! w'x' (momentum levels)                 [{xm units} m/s]
      xm_forcing,      & ! xm forcings (thermodynamic levels)     [{xm units}/s]
      C7_Skw_fnc,      & ! C_7 parameter with Sk_w applied        [-]
      xpthvp,          & ! x'th_v' (momentum levels)              [{xm units} K]
      wp2_zt,          & ! w'^2 interpolated to thermodynamic levels [m^2/s^2]
      a1_zt,           & ! a_1 interpolated to thermodynamic levels  [-]
      wp3,             & ! w'^3 (thermodynamic levels)               [m^3/s^3]
      wpxp_upper_lim,  & ! Keeps correlations from becoming > 1.  [units vary]
      wpxp_lower_lim     ! Keeps correlations from becoming < -1. [units vary]

    ! Output Variable
    real, intent(out), dimension(2*gr%nnzp) ::  & 
      rhs  ! Right-hand side of band diag. matrix. (LAPACK)

    ! Local Variables.

    ! For "over-implicit" weighted time step.
    ! This vector holds output from the LHS (implicit) portion of a term at a
    ! given vertical level.  This output is weighted and applied to the RHS.
    ! This is used if the implicit portion of the term is "over-implicit", which
    ! means that the LHS contribution is given extra weight (>1) in order to
    ! increase numerical stability.  A weighted factor must then be applied to
    ! the RHS in order to balance the weight.
    real, dimension(3) :: lhs_fnc_output

    ! Indices
    integer :: k, km1, kp1, k_xm, k_wpxp


    integer :: & 
      ixm_f, & 
      iwpxp_bp, & 
      iwpxp_pr3, & 
      iwpxp_sicl, &
      iwpxp_ta

    select case ( trim( solve_type ) )
    case ( "rtm" )  ! rtm/wprtp budget terms
      ixm_f      = irtm_forcing
      iwpxp_bp   = iwprtp_bp
      iwpxp_pr3  = iwprtp_pr3
      iwpxp_sicl = iwprtp_sicl
      iwpxp_ta   = iwprtp_ta
    case ( "thlm" ) ! thlm/wpthlp budget terms
      ixm_f      = ithlm_forcing
      iwpxp_bp   = iwpthlp_bp
      iwpxp_pr3  = iwpthlp_pr3
      iwpxp_sicl = iwpthlp_sicl
      iwpxp_ta   = iwpthlp_ta
    case default    ! this includes the sclrm case
      ixm_f      = 0
      iwpxp_bp   = 0
      iwpxp_pr3  = 0
      iwpxp_sicl = 0
      iwpxp_ta   = 0
    end select


    ! Initialize the right-hand side vector to 0.
    rhs = 0.0

    do k = 2, gr%nnzp-1, 1

      ! Define indices

      km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nnzp )

      k_xm   = 2*k - 1
      k_wpxp = 2*k


      !!!!!***** xm *****!!!!!

      ! xm: Right-hand side (explicit xm portion of the code).

      ! RHS time tendency.
      rhs(k_xm) = real( rhs(k_xm) + xm(k) / dt )

      ! RHS xm forcings.
      ! Note: xm forcings include the effects of microphysics,
      !       cloud water sedimentation, radiation, and any
      !       imposed forcings on xm.
      rhs(k_xm) = rhs(k_xm) + xm_forcing(k)

      if ( l_stats_samp ) then

        ! Statistics: explicit contributions for xm
        !             (including microphysics/radiation).

        ! xm forcings term is completely explicit; call stat_update_var_pt.
        call stat_update_var_pt( ixm_f, k, xm_forcing(k), zt )

      endif ! l_stats_samp


      !!!!!***** w'x' *****!!!!!

      ! w'x': Right-hand side (explicit w'x' portion of the code).

      ! RHS buoyancy production (bp) term and pressure term 3 (pr3).
      rhs(k_wpxp) & 
      = rhs(k_wpxp) & 
      + wpxp_terms_bp_pr3_rhs( C7_Skw_fnc(k), xpthvp(k) )

      ! RHS time tendency.
      if ( l_iter ) rhs(k_wpxp) =  & 
           real( rhs(k_wpxp) + wpxp(k) / dt )

      ! RHS portion of semi-implicit clipping (sicl) term.
      if ( l_clip_semi_implicit ) then

        rhs(k_wpxp) & 
        = rhs(k_wpxp) & 
        + clip_semi_imp_rhs( dt, wpxp(k), & 
                             .true., wpxp_upper_lim(k), & 
                             .true., wpxp_lower_lim(k) )

      endif

      ! RHS contribution from "over-implicit" weighted time step
      ! for LHS turbulent advection (ta) term.
      !
      ! Note:  An "over-implicit" weighted time step is applied to this term.
      !        The weight of the implicit portion of this term is controlled
      !        by the factor gamma_over_implicit_ts (abbreviated "gamma" in the
      !        expression below).  A factor is added to the right-hand side of
      !        the equation in order to balance a weight that is not equal to 1,
      !        such that:
      !             -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
      !        where X is the variable that is being solved for in a predictive
      !        equation (w'x' in this case), y(t) is the linearized portion of
      !        the term that gets treated implicitly, and RHS is the portion of
      !        the term that is always treated explicitly (in the case of the
      !        w'x' turbulent advection term, RHS = 0).  A weight of greater
      !        than 1 can be applied to make the term more numerically stable.
      lhs_fnc_output(1:3)  &
      = wpxp_term_ta_lhs( wp2_zt(kp1), wp2_zt(k),  &
                          a1_zt(kp1), a1_zt(k),  &
                          wp3(kp1), wp3(k), gr%dzm(k), k )
      rhs(k_wpxp)  &
      = rhs(k_wpxp)  &
      + ( 1.0 - gamma_over_implicit_ts )  &
      * ( - lhs_fnc_output(1) * wpxp(kp1)  &
          - lhs_fnc_output(2) * wpxp(k)  &
          - lhs_fnc_output(3) * wpxp(km1) )


      if ( l_stats_samp ) then

        ! Statistics: explicit contributions for wpxp.

        ! w'x' term bp is completely explicit; call stat_update_var_pt.
        ! Note:  To find the contribution of w'x' term bp, substitute 0 for the
        !        C_7 skewness function input to function wpxp_terms_bp_pr3_rhs.
        call stat_update_var_pt( iwpxp_bp, k, & 
            wpxp_terms_bp_pr3_rhs( 0.0, xpthvp(k) ), zm )

        ! w'x' term pr3 is completely explicit; call stat_update_var_pt.
        ! Note:  To find the contribution of w'x' term pr3, add 1 to the
        !        C_7 skewness function input to function wpxp_terms_bp_pr2_rhs.
        call stat_update_var_pt( iwpxp_pr3, k, & 
            wpxp_terms_bp_pr3_rhs( (1.0+C7_Skw_fnc(k)),xpthvp(k)), zm)

        ! w'x' term sicl has both implicit and explicit components; call
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
        ! subtracts the value sent in, reverse the sign on clip_semi_imp_rhs.
        if ( l_clip_semi_implicit ) then
          call stat_begin_update_pt( iwpxp_sicl, k, & 
             -clip_semi_imp_rhs( dt, wpxp(k), & 
                                 .true., wpxp_upper_lim(k), & 
                                 .true., wpxp_lower_lim(k) ), zm )
        endif

        ! w'x' term ta is normally completely implicit.  However, there is a
        ! RHS contribution from the "over-implicit" weighted time step.  A
        ! weighting factor of greater than 1 may be used to make the term more
        ! numerically stable (see note above for RHS contribution from
        ! "over-implicit" weighted time step for LHS turbulent advection (ta)
        ! term).  Therefore, w'x' term ta has both implicit and explicit
        ! components; call stat_begin_update_pt.  Since stat_begin_update_pt
        ! automatically subtracts the value sent in, reverse the sign on the
        ! input value.
        lhs_fnc_output(1:3)  &
        = wpxp_term_ta_lhs( wp2_zt(kp1), wp2_zt(k),  &
                            a1_zt(kp1), a1_zt(k),  &
                            wp3(kp1), wp3(k), gr%dzm(k), k )
        call stat_begin_update_pt( iwpxp_ta, k, &
              - ( 1.0 - gamma_over_implicit_ts )  &
              * ( - lhs_fnc_output(1) * wpxp(kp1)  &
                  - lhs_fnc_output(2) * wpxp(k)  &
                  - lhs_fnc_output(3) * wpxp(km1) ), zm )


      endif ! l_stats_samp

    enddo ! k=2..gr%nnzp-1


    ! Boundary conditions.

    ! Both the mean (xm) and the turbulent flux (wpxp) use fixed-point
    ! boundary conditions.  Therefore, anything set in the above loop
    ! at both the upper and lower boundaries would be overwritten here.
    ! However, the above loop does not extend to the boundary levels.
    ! An array with a value of 1 at the main diagonal on the left-hand
    ! side and with values of 0 at all other diagonals on the left-hand
    ! side will preserve the right-hand side value at that level.

    ! Lower boundary
    k      = 1
    k_xm   = 2*k - 1
    k_wpxp = 2*k
    ! The value of xm at the lower boundary will remain the same.
    ! However, the value of xm at the lower boundary gets overwritten
    ! after the matrix is solved for the next timestep, such
    ! that xm(1) = xm(2).
    rhs(k_xm)   = xm(k)
    ! The value of w'x' at the lower boundary will remain the same.
    ! The surface value of w'x' is set elsewhere
    ! (case-specific information).
    rhs(k_wpxp) = wpxp(k)

    ! Upper boundary
    k      = gr%nnzp
    k_xm   = 2*k - 1
    k_wpxp = 2*k
    ! The value of xm at the upper boundary will remain the same.
    rhs(k_xm)   = xm(k)
    ! The value of w'x' at the upper boundary will be 0.
    rhs(k_wpxp) = 0.0


  end subroutine xm_wpxp_rhs

  !=============================================================================
  subroutine xm_wpxp_solve( nrhs, lhs, rhs, solution, err_code, rcond )

    ! Description:
    !   Solve for xm / w'x' using the band diagonal solver.

    ! References:
    !   None
    !------------------------------------------------------------------------

    use grid_class, only: & 
      gr ! Variable(s)

    use lapack_wrap, only:  & 
      band_solve,  & ! Procedure(s)
      band_solvex

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nrhs ! Number of rhs vectors

    ! Input/Output Variables
    real, intent(inout), dimension(nsup+nsub+1,2*gr%nnzp) :: & 
      lhs  ! Implicit contributions to wpxp/xm (band diag. matrix in LAPACK storage)

    real, intent(inout), dimension(2*gr%nnzp,nrhs) ::  & 
      rhs      ! Right-hand side of band diag. matrix. (LAPACK storage)

    real, intent(out), dimension(2*gr%nnzp,nrhs) ::  & 
      solution ! Solution to band diagonal system (LAPACK storage)

    ! Output Variables
    integer, intent(out) :: err_code

    real, optional, intent(out) :: rcond ! Est. of the reciprocal of the condition #


    if ( present( rcond ) ) then
      ! Perform LU decomp and solve system (LAPACK with diagnostics)
      call band_solvex( "xm_wpxp", nsup, nsub, 2*gr%nnzp, nrhs, & 
                        lhs, rhs, solution, rcond, err_code )


    else
      ! Perform LU decomp and solve system (LAPACK)
      call band_solve( "xm_wpxp", nsup, nsub, 2*gr%nnzp, nrhs, & 
                       lhs, rhs, solution, err_code )
    end if


    return
  end subroutine xm_wpxp_solve

!===============================================================================
  subroutine xm_wpxp_clipping_and_stats &
             ( solve_type, dt, wp2, xp2, xm, wpxp, solution, rcond )
! Description:
!   Clips and computes implicit stats for an artitrary xm and wpxp
!
! References:
!   None
!-------------------------------------------------------------------------------
    use grid_class, only: & 
      gr ! Variable(s)

    use model_flags, only: &
        l_clip_semi_implicit ! Variable(s)

    use stats_precision, only:  & 
        time_precision ! Variable(s)

    use pos_definite_mod, only:  & 
        pos_definite_adj ! Procedure(s)

    use clip_explicit, only: & 
        clip_covariance ! Procedure(s)

    use model_flags, only: & 
        l_pos_def, &  ! Logical for whether to apply the positive definite scheme to rtm
        l_clip_turb_adv ! Logical for whether to clip xm when wpxp is clipped

    use stats_type, only: & 
        stat_begin_update,  & ! Procedure(s)
        stat_update_var_pt, & 
        stat_end_update_pt, & 
        stat_end_update,  & 
        stat_update_var, & 
        stat_modify

    use stats_variables, only: & 
        zt,  & ! Variable(s)
        zm, & 
        sfc, & 
        irtm_bt, & 
        irtm_ta, & 
        irtm_ma, & 
        irtm_matrix_condt_num, & 
        irtm_pd, & 
        iwprtp_bt, & 
        iwprtp_ma, & 
        iwprtp_ta, & 
        iwprtp_tp, & 
        iwprtp_ac, & 
        iwprtp_pr1, & 
        iwprtp_pr2, & 
        iwprtp_dp1, & 
        iwprtp_pd, & 
        iwprtp_sicl, & 
        ithlm_bt, & 
        ithlm_ta, & 
        ithlm_ma, & 
        ithlm_matrix_condt_num, & 
        iwpthlp_bt, & 
        iwpthlp_ma, & 
        iwpthlp_ta, & 
        iwpthlp_tp, & 
        iwpthlp_ac, & 
        iwpthlp_pr1, & 
        iwpthlp_pr2, & 
        iwpthlp_dp1, & 
        iwpthlp_sicl, & 
        l_stats_samp, & 
        ztscr01, & 
        ztscr02, & 
        ztscr03, & 
        ztscr04, & 
        ztscr05, & 
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
        zmscr13, & 
        zmscr14, & 
        zmscr15

    implicit none

    ! Input Variables
    character(len=*), intent(in) ::  & 
      solve_type  ! Variables being solved for.

    real(kind=time_precision), intent(in) ::  & 
      dt  ! Timestep   [s]

    real, intent(in), dimension(gr%nnzp) ::  & 
      wp2,  & ! w'^2 (momentum levels)        [m^2/s^2]
      xp2     ! x'^2 (momentum levels)        [{xm units}^2]

    real, intent(in), dimension(2*gr%nnzp) :: &
      solution ! The <t+1> value of xm and wpxp   [units vary]

    real, intent(in) :: &
      rcond ! Reciprocal of the estimated condition number (from computing A^-1)

    real, intent(inout), dimension(gr%nnzp) :: & 
      xm, & ! The mean x field  [units vary]
      wpxp  ! The flux of x     [units vary m/s]

    ! Local Variables
    character(len=25) :: & 
      solve_type_cl ! solve_type used for clipping statistics.

    real, dimension(gr%nnzp) :: & 
      xm_n ! Old value of xm for positive definite scheme     [units vary]

    real, dimension(gr%nnzp) :: & 
      wpxp_pd, xm_pd ! Change in xm and wpxp due to the pos. def. scheme

    real, dimension(gr%nnzp) :: &
      wpxp_chnge  ! Net change in w'x' due to clipping        [units vary]

    ! Indices
    integer :: k, km1, kp1
    integer :: k_xm, k_wpxp

    integer :: & 
      ixm_bt, & 
      ixm_ta, & 
      ixm_ma, & 
      ixm_matrix_condt_num, & 
      ixm_pd, & 
      iwpxp_bt, & 
      iwpxp_ma, & 
      iwpxp_ta, & 
      iwpxp_tp, & 
      iwpxp_ac, & 
      iwpxp_pr1, & 
      iwpxp_pr2, & 
      iwpxp_dp1, & 
      iwpxp_pd, & 
      iwpxp_sicl

    ! ----- Begin code ------

    select case ( trim( solve_type ) )
    case ( "rtm" ) ! rtm/wprtp budget terms
      ixm_bt     = irtm_bt
      ixm_ta     = irtm_ta
      ixm_ma     = irtm_ma
      ixm_pd     = irtm_pd
      iwpxp_bt   = iwprtp_bt
      iwpxp_ma   = iwprtp_ma
      iwpxp_ta   = iwprtp_ta
      iwpxp_tp   = iwprtp_tp
      iwpxp_ac   = iwprtp_ac
      iwpxp_pr1  = iwprtp_pr1
      iwpxp_pr2  = iwprtp_pr2
      iwpxp_dp1  = iwprtp_dp1
      iwpxp_pd   = iwprtp_pd
      iwpxp_sicl = iwprtp_sicl

      ! This is a diagnostic from inverting the matrix, not a budget
      ixm_matrix_condt_num = irtm_matrix_condt_num
    case ( "thlm" ) ! thlm/wpthlp budget terms
      ixm_bt     = ithlm_bt
      ixm_ta     = ithlm_ta
      ixm_ma     = ithlm_ma
      ixm_pd     = 0
      iwpxp_bt   = iwpthlp_bt
      iwpxp_ma   = iwpthlp_ma
      iwpxp_ta   = iwpthlp_ta
      iwpxp_tp   = iwpthlp_tp
      iwpxp_ac   = iwpthlp_ac
      iwpxp_pr1  = iwpthlp_pr1
      iwpxp_pr2  = iwpthlp_pr2
      iwpxp_dp1  = iwpthlp_dp1
      iwpxp_pd   = 0
      iwpxp_sicl = iwpthlp_sicl

      ! This is a diagnostic from inverting the matrix, not a budget
      ixm_matrix_condt_num = ithlm_matrix_condt_num

    case default  ! this includes the sclrm case
      ixm_bt     = 0
      ixm_ta     = 0
      ixm_ma     = 0
      ixm_pd     = 0
      iwpxp_bt   = 0
      iwpxp_ma   = 0
      iwpxp_ta   = 0
      iwpxp_tp   = 0
      iwpxp_ac   = 0
      iwpxp_pr1  = 0
      iwpxp_pr2  = 0
      iwpxp_dp1  = 0
      iwpxp_pd   = 0
      iwpxp_sicl = 0

      ixm_matrix_condt_num = 0
    end select


    if ( l_stats_samp ) then


      ! xm total time tendency ( 1st calculation)
      call stat_begin_update( ixm_bt, real( xm /dt ), zt )

      ! wpxp is clipped after xp2 is updated in subroutine advance_xp2_xpyp and
      ! after wp2 is updated in subroutine advance_wp2_wp3.  The overall time
      ! tendency must include the effects of those two clippings, as well.
      ! Therefore, the wpxp total time tendency term is just being modified in
      ! advance_xm_wpxp_module.F90, rather than being entirely contained in
      ! advance_xm_wpxp_module.F90.
      !!  wpxp total time tendency (1st calculation)
      !call stat_begin_update( iwpxp_bt, real( wpxp / dt ), zm )

      ! wpxp total time tendency (1st calculation in advance_xm_wpxp_module.F90)
      call stat_modify( iwpxp_bt, real( -wpxp / dt ), zm )
      ! Brian Griffin; July 5, 2008.

    end if ! l_stats_samp

    ! Copy result into output arrays

    do k=1, gr%nnzp, 1

      k_xm   = 2 * k - 1
      k_wpxp = 2 * k

      xm_n(k) = xm(k)

      xm(k)   = solution(k_xm)
      wpxp(k) = solution(k_wpxp)

    end do ! k=1..gr%nnzp

    ! Boundary condition on xm

    !xm(1) = 2. * xm(2) - xm(3)
    !xm(gr%nnzp) = 2. * xm(gr%nnzp-1) - xm(gr%nnzp-2)
    xm(1) = xm(2)
    !xm(gr%nnzp) = xm(gr%nnzp-1)

    if ( l_stats_samp ) then

      if ( ixm_matrix_condt_num > 0 ) then
        ! Est. of the condition number of the mean/flux LHS matrix
        call stat_update_var_pt( ixm_matrix_condt_num, 1, 1.0 / rcond, sfc )
      end if

      do k = 2, gr%nnzp-1

        km1 = max( k-1, 1 )
        kp1 = min( k+1, gr%nnzp )

        ! Finalize implicit contributions for xm

        ! xm term ma is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( ixm_ma, k, & 
            ztscr01(k) * xm(km1) & 
          + ztscr02(k) * xm(k) & 
          + ztscr03(k) * xm(kp1), zt )

        ! xm term ta is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( ixm_ta, k, & 
            ztscr04(k) * wpxp(km1) & 
          + ztscr05(k) * wpxp(k), zt )

        ! Finalize implicit contributions for wpxp

        ! w'x' term ma is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwpxp_ma, k, & 
            zmscr01(k) * wpxp(km1) & 
          + zmscr02(k) * wpxp(k) & 
          + zmscr03(k) * wpxp(kp1), zm )

        ! w'x' term ta is normally completely implicit.  However, due to the
        ! RHS contribution from the "over-implicit" weighted time step,
        ! w'x' term ta has both implicit and explicit components;
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwpxp_ta, k, & 
            zmscr04(k) * wpxp(km1) & 
          + zmscr05(k) * wpxp(k) & 
          + zmscr06(k) * wpxp(kp1), zm )

        ! w'x' term tp is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwpxp_tp, k, & 
            zmscr07(k) * xm(k) & 
          + zmscr08(k) * xm(kp1), zm )

        ! w'x' term ac is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwpxp_ac, k, & 
            zmscr09(k) * wpxp(k), zm )

        ! w'x' term pr1 is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwpxp_pr1, k, & 
            zmscr10(k) * wpxp(k), zm )

        ! w'x' term pr2 is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwpxp_pr2, k, & 
            zmscr11(k) * wpxp(k), zm )

        ! w'x' term dp1 is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwpxp_dp1, k, & 
            zmscr12(k) * wpxp(km1) & 
          + zmscr13(k) * wpxp(k) & 
          + zmscr14(k) * wpxp(kp1), zm )

        ! w'x' term sicl has both implicit and explicit components;
        ! call stat_end_update_pt.
        if ( l_clip_semi_implicit ) then
          call stat_end_update_pt( iwpxp_sicl, k, & 
              zmscr15(k) * wpxp(k), zm )
        endif

      enddo ! 1..gr%nnzp

    endif ! l_stats_samp


    ! Apply a flux limiting positive definite scheme if the solution
    ! for the mean field is negative and we're determining total water
    if ( trim( solve_type ) == "rtm" .and. l_pos_def .and. any( xm < 0.0 ) ) then

      call pos_definite_adj( dt, "zt", xm, wpxp, & 
                             xm_n, xm_pd, wpxp_pd )

    else
      ! For stats purposes
      xm_pd   = 0.0
      wpxp_pd = 0.0

    end if ! solve_type == "rtm" and rtm <n+1> less than 0

    if ( l_stats_samp ) then

      call stat_update_var( iwpxp_pd, wpxp_pd(1:gr%nnzp), zm )

      call stat_update_var( ixm_pd, xm_pd(1:gr%nnzp), zt )

    end if

    ! Use solve_type to find solve_type_cl, which is used
    ! in subroutine clip_covariance.
    select case ( trim( solve_type ) )
    case ( "rtm" )
      solve_type_cl = "wprtp"
    case ( "thlm" )
      solve_type_cl = "wpthlp"
    case default
      solve_type_cl = "wpsclrp"
    end select

    ! Clipping for w'x'
    ! Clipping w'x' at each vertical level, based on the
    ! correlation of w and x at each vertical level, such that:
    ! corr_(w,x) = w'x' / [ sqrt(w'^2) * sqrt(x'^2) ];
    ! -1 <= corr_(w,x) <= 1.
    ! Since w'^2, x'^2, and w'x' are updated in different places
    ! from each other, clipping for w'x' has to be done three times
    ! (three times each for w'r_t', w'th_l', and w'sclr').  This is
    ! the second instance of w'x' clipping.
    call clip_covariance( solve_type_cl, .false.,  & 
                          .false., dt, wp2, xp2,  & 
                          wpxp, wpxp_chnge )

    ! Adjusting xm based on clipping for w'x'.
    if ( any( wpxp_chnge /= 0.0 ) .and. l_clip_turb_adv ) then
      call xm_correction_wpxp_cl( solve_type, dt, wpxp_chnge, gr%dzt, &
                                  xm )
    endif

    if ( l_stats_samp ) then

      ! xm time tendency (2nd calculation)
      call stat_end_update( ixm_bt, real( xm / dt ), zt )

      ! wpxp is clipped after xp2 is updated in subroutine advance_xp2_xpyp and
      ! after wp2 is updated in subroutine advance_wp2_wp3.  The overall time
      ! tendency must include the effects of those two clippings, as well.
      ! Therefore, the wpxp total time tendency term is just being modified in
      ! advance_xm_wpxp_module.F90, rather than being entirely contained in
      ! advance_xm_wpxp_module.F90.
      !! wpxp time tendency (2nd calculation)
      !call stat_end_update( iwpxp_bt, real( wpxp / dt ), zm )

      ! wpxp time tendency (2nd calculation in advance_xm_wpxp_module.F90)
      call stat_modify( iwpxp_bt, real( wpxp / dt ), zm )
      ! Brian Griffin; July 5, 2008.

    endif

    return
  end subroutine xm_wpxp_clipping_and_stats

  !=============================================================================
  pure function xm_term_ta_lhs( dzt ) & 
  result( lhs )

    ! Description:
    ! Turbulent advection of xm:  implicit portion of the code.
    !
    ! The d(xm)/dt equation contains a turbulent advection term:
    !
    ! - d(w'x')/dz.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - d( w'x'(t+1) )/dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign
    !        is reversed and the leading "-" in front of the term is changed
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of w'x' being used is from
    ! the next timestep, which is being advanced to in solving the d(xm)/dt and
    ! d(w'x')/dt equations.
    !
    ! This term is discretized as follows:
    !
    ! While the values of xm are found on the thermodynamic levels, the values
    ! of w'x' are found on the momentum levels.  The derivative of w'x' is taken
    ! over the intermediate (central) thermodynamic level, yielding the desired
    ! results.
    !
    ! ===================wpxp================================== m(k)
    !
    ! -----------------------------d(wpxp)/dz------------------ t(k)
    !
    ! ===================wpxpm1================================ m(k-1)
    !
    ! The vertical indices m(k), t(k), and m(k-1) correspond with altitudes
    ! zm(k), zt(k), and zm(k-1), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! dzt(k) = 1 / ( zm(k) - zm(k-1) )

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      k_mdiag   = 1,    & ! Momentum superdiagonal index.
      km1_mdiag = 2       ! Momentum subdiagonal index.

    ! Input Variables
    real, intent(in) :: & 
      dzt   ! Inverse of grid spacing (k)   [1/m]

    ! Return Variable
    real, dimension(2) :: lhs

    ! Momentum superdiagonal [ x wpxp(k,<t+1>) ]
    lhs(k_mdiag) & 
    = + dzt

    ! Momentum subdiagonal [ x wpxp(k-1,<t+1>) ]
    lhs(km1_mdiag) & 
    = - dzt

    return
  end function xm_term_ta_lhs

  !=============================================================================
  pure function wpxp_term_ta_lhs( wp2_ztp1, wp2_zt,  & 
                                  a1_ztp1, a1_zt, & 
                                  wp3p1, wp3, dzm, level ) & 
  result( lhs )

    ! Description:
    ! Turbulent advection of w'x':  implicit portion of the code.
    !
    ! The d(w'x')/dt equation contains a turbulent advection term:
    !
    ! - d(w'^2x')/dz.
    !
    ! A substitution is made in order to close the turbulent advection term,
    ! such that:
    !
    ! w'^2x' = a_1 * ( w'^3 / w'^2 ) * w'x',
    !
    ! where a_1 is a variable that is a function of sigma_sqd_w.  The turbulent
    ! advection term becomes:
    !
    ! - d [ a_1 * ( w'^3 / w'^2 ) * w'x' ] / dz.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - d [ a_1 * ( w'^3 / w'^2 ) * w'x'(t+1) ] / dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign
    !        is reversed and the leading "-" in front of the term is changed
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of w'x' being used is from
    ! the next timestep, which is being advanced to in solving the d(w'x')/dt
    ! equation.
    !
    ! This term is discretized as follows:
    !
    ! The values of w'x', w'^2, and a_1 are found on the momentum levels, while
    ! the values of w'^3 are found on the thermodynamic levels.  Each of the
    ! variables w'x', w'^2, and a_1 are interpolated to the intermediate
    ! thermodynamic levels.  The values of the mathematical expression (called F
    ! here) within the dF/dz term are computed on the thermodynamic levels.
    ! Then, the derivative (d/dz) of the expression (F) is taken over the
    ! central momentum level, yielding the desired result.  In this function,
    ! the values of F are as follows:
    !
    ! F = a_1(t) * ( w'^3(t) / w'^2(t) ) * w'x'(t+1);
    !
    ! where the timestep index (t) stands for the index of the current timestep.
    !
    !
    ! =a1p1========wp2p1========wpxpp1========================= m(k+1)
    !
    ! -----a1(interp)---wp2(interp)---wpxp(interp)---wp3p1----- t(k+1)
    !
    ! =a1==========wp2==========wpxp===================dF/dz=== m(k)
    !
    ! -----a1(interp)---wp2(interp)---wpxp(interp)---wp3------- t(k)
    !
    ! =a1m1========wp2m1========wpxpm1========================= m(k-1)
    !
    ! The vertical indices m(k+1), t(k+1), m(k), t(k), and m(k-1) correspond
    ! with altitudes zm(k+1), zt(k+1), zm(k), zt(k), and zm(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is used
    ! for momentum levels.
    !
    ! dzm(k) = 1 / ( zt(k+1) - zt(k) )

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable; gr%weights_zm2zt

    use constants, only: &
        wtol_sqd ! Constant; minimum threshold for w'^2 [m^2/s^2]

!    use model_flags, only:  &
!        l_standard_term_ta

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_mdiag = 1,    & ! Momentum superdiagonal index.
      k_mdiag   = 2,    & ! Momentum main diagonal index.
      km1_mdiag = 3       ! Momentum subdiagonal index.

    integer, parameter :: & 
      m_above = 1,    & ! Index for upper momentum level grid weight.
      m_below = 2       ! Index for lower momentum level grid weight.

    ! Input Variables
    real, intent(in) :: & 
      wp2_ztp1,    & ! w'^2 interpolated to thermodynamic level (k+1) [m^2/s^2]
      wp2_zt,      & ! w'^2 interpolated to thermodynamic level (k)   [m^2/s^2]
      a1_ztp1,     & ! a_1 interpolated to thermodynamic level (k+1)  [-]
      a1_zt,       & ! a_1 interpolated to thermodynamic level (k)    [-]
      wp3p1,       & ! w'^3(k+1)                                      [m^3/s^3]
      wp3,         & ! w'^3(k)                                        [m^3/s^3]
      dzm            ! Inverse of grid spacing (k)                    [1/m]

    integer, intent(in) :: & 
      level ! Central momentum level (on which calculation occurs).

    ! Return Variable
    real, dimension(3) :: lhs

    ! Local Variables
    integer :: & 
      tkp1,  & ! Thermodynamic level directly above central momentum level.
      tk       ! Thermodynamic level directly below central momentum level.

    ! Thermodynamic level (k+1) is between momentum level (k+1)
    ! and momentum level (k).
    tkp1 = level + 1

    ! Thermodynamic level (k) is between momentum level (k)
    ! and momentum level (k-1).
    tk = level

    ! Note:  The w'x' turbulent advection term, which is
    !        - d [ a_1 * ( w'^3 / w'^2 ) * w'x' ] / dz, still keeps the a_1 term
    !        inside the derivative, unlike the w'^3 equation (found in
    !        advance_wp2_wp3_mod.F90) and the equations for r_t'^2, th_l'^2,
    !        r_t'th_l', u'^2, v'^2, sclr'r_t', sclr'th_l', and sclr'^2 (found in
    !        advance_xp2_xpyp_module.F90).  Brian.

    ! Always use the standard discretization for the w'x' turbulent advection
    ! term.  Brian.
    !if ( l_standard_term_ta ) then

    ! The turbulent advection term is discretized normally, in accordance
    ! with the model equations found in the documentation and the description
    ! listed above.

    ! Momentum superdiagonal: [ x wpxp(k+1,<t+1>) ]
    lhs(kp1_mdiag) & 
    = + dzm & 
        * a1_ztp1 * ( wp3p1 / max( wp2_ztp1, wtol_sqd ) ) & 
        * gr%weights_zm2zt(m_above,tkp1)

    ! Momentum main diagonal: [ x wpxp(k,<t+1>) ]
    lhs(k_mdiag) & 
    = + dzm & 
        * (   a1_ztp1 * ( wp3p1 / max( wp2_ztp1, wtol_sqd ) ) & 
              * gr%weights_zm2zt(m_below,tkp1) & 
            - a1_zt * ( wp3 / max( wp2_zt, wtol_sqd ) ) & 
              * gr%weights_zm2zt(m_above,tk) & 
          )

    ! Momentum subdiagonal: [ x wpxp(k-1,<t+1>) ]
    lhs(km1_mdiag) & 
    = - dzm & 
        * a1_zt * ( wp3 / max( wp2_zt, wtol_sqd ) ) & 
        * gr%weights_zm2zt(m_below,tk)

    !endif


    return
  end function wpxp_term_ta_lhs

  !=============================================================================
  pure function wpxp_term_tp_lhs( wp2, dzm ) & 
  result( lhs )

    ! Description:
    ! Turbulent production of w'x':  implicit portion of the code.
    !
    ! The d(w'x')/dt equation contains a turbulent production term:
    !
    ! - w'^2 d(xm)/dz.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - w'^2 * d( xm(t+1) )/dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign
    !        is reversed and the leading "-" in front of the term is changed
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of xm being used is from the
    ! next timestep, which is being advanced to in solving the d(w'x')/dt and
    ! d(xm)/dt equations.
    !
    ! This term is discretized as follows:
    !
    ! The values of xm are found on thermodynamic levels, while the values of
    ! w'^2 are found on momentum levels.  The derivative of xm is taken over the
    ! intermediate (central) momentum level, where it is multiplied by w'^2,
    ! yielding the desired result.
    !
    ! ---------------------------xmp1-------------------------- t(k+1)
    !
    ! ==========wp2=====================d(xm)/dz=============== m(k)
    !
    ! ---------------------------xm---------------------------- t(k)
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
      k_tdiag = 2         ! Thermodynamic subdiagonal index.

    ! Input Variables
    real, intent(in) :: & 
      wp2,    & ! w'^2(k)                       [m^2/s^2]
      dzm       ! Inverse of grid spacing (k)   [1/m]

    ! Return Variable
    real, dimension(2) :: lhs

    ! Thermodynamic superdiagonal [ x xm(k+1,<t+1>) ]
    lhs(kp1_tdiag) & 
    = + wp2 * dzm

    ! Thermodynamic subdiagonal [ x xm(k,<t+1>) ]
    lhs(k_tdiag) & 
    = - wp2 * dzm

    return
  end function wpxp_term_tp_lhs

  !=============================================================================
  pure function wpxp_terms_ac_pr2_lhs( C7_Skw_fnc,  & 
                                       wm_ztp1, wm_zt, dzm ) & 
  result( lhs )

    ! Description:
    ! Accumulation of w'x' and w'x' pressure term 2:  implicit portion of the
    ! code.
    !
    ! The d(w'x')/dt equation contains an accumulation term:
    !
    ! - w'x' dw/dz;
    !
    ! and pressure term 2:
    !
    ! + C_7 w'x' dw/dz.
    !
    ! Both the w'x' accumulation term and pressure term 2 are completely
    ! implicit.  The accumulation term and pressure term 2 are combined and
    ! solved together as:
    !
    ! - ( 1 - C_7 ) * w'x'(t+1) * dw/dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign
    !        is reversed and the leading "-" in front of the term is changed
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of w'x' being used is from
    ! the next timestep, which is being advanced to in solving the d(w'x')/dt
    ! equation.
    !
    ! The terms are discretized as follows:
    !
    ! The values of w'x' are found on momentum levels, while the values of wm_zt
    ! (mean vertical velocity on thermodynamic levels) are found on
    ! thermodynamic levels.  The vertical derivative of wm_zt is taken over the
    ! intermediate (central) momentum level.  It is then multiplied by w'x'
    ! (implicitly calculated at timestep (t+1)) and the coefficients to yield
    ! the desired results.
    !
    ! -------wm_ztp1------------------------------------------- t(k+1)
    !
    ! ===============d(wm_zt)/dz============wpxp=============== m(k)
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
      C7_Skw_fnc,  & ! C_7 parameter with Sk_w applied (k)             [-]
      wm_ztp1,     & ! w wind component on thermodynamic level (k+1)   [m/s]
      wm_zt,       & ! w wind component on thermodynamic level (k)     [m/s]
      dzm            ! Inverse of grid spacing (k)                     [1/m]


    ! Return Variable
    real :: lhs

    ! Momentum main diagonal: [ x wpxp(k,<t+1>) ]
    lhs & 
    = + ( 1.0 - C7_Skw_fnc ) * dzm * ( wm_ztp1 - wm_zt )

    return
  end function wpxp_terms_ac_pr2_lhs

  !=============================================================================
  pure function wpxp_term_pr1_lhs( C6x_Skw_fnc, tau_zm ) & 
  result( lhs )

    ! Description
    ! Pressure term 1 for w'x':  implicit portion of the code.
    !
    ! The d(w'x')/dt equation contains pressure term 1:
    !
    ! - ( C_6 / tau_m ) w'x'.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - ( C_6 / tau_m ) w'x'(t+1)
    !
    ! Note:  When the term is brought over to the left-hand side, the sign
    !        is reversed and the leading "-" in front of the term is changed
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of w'x' being used is from
    ! the next timestep, which is being advanced to in solving the d(w'x')/dt
    ! equation.
    !
    ! The values of w'x' are found on the momentum levels.  The values of the
    ! C_6 skewness function and time-scale tau_m are also found on the momentum
    ! levels.

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) :: & 
      C6x_Skw_fnc,  & ! C_6x parameter with Sk_w applied (k)   [-]
      tau_zm          ! Time-scale tau at momentum level (k)   [s]

    ! Return Variable
    real :: lhs

    ! Momentum main diagonal: [ x wpxp(k,<t+1>) ]
    lhs & 
    = + C6x_Skw_fnc / tau_zm

    return
  end function wpxp_term_pr1_lhs

  !=============================================================================
  pure function wpxp_terms_bp_pr3_rhs( C7_Skw_fnc, xpthvp ) & 
  result( rhs )

    ! Description:
    ! Buoyancy production of w'x' and w'x' pressure term 3:  explicit portion of
    ! the code.
    !
    ! The d(w'x')/dt equation contains a buoyancy production term:
    !
    ! + (g/th_0) x'th_v';
    !
    ! and pressure term 3:
    !
    ! - C_7 (g/th_0) x'th_v'.
    !
    ! Both the w'x' buoyancy production term and pressure term 3 are completely
    ! explicit.  The buoyancy production term and pressure term 3 are combined
    ! and solved together as:
    !
    ! + ( 1 - C_7 ) * (g/th_0) * x'th_v'.

    ! References:
    !-----------------------------------------------------------------------

    use constants, only: & 
    ! Variable(s) 
        grav ! Gravitational acceleration [m/s^2]
    use parameters_model, only: & 
    ! Variable(s) 
        T0  ! Reference temperature      [K]

    implicit none

    ! Input Variables
    real, intent(in) :: & 
      C7_Skw_fnc,  & ! C_7 parameter with Sk_w applied (k)   [-]
      xpthvp         ! x'th_v'(k)                            [K {xm units}]

    ! Return Variable
    real :: rhs

    rhs & 
    = grav/T0 * ( 1.0 - C7_Skw_fnc ) * xpthvp

    return
  end function wpxp_terms_bp_pr3_rhs

  !=============================================================================
  subroutine xm_correction_wpxp_cl( solve_type, dt, wpxp_chnge, dzt, &
                                    xm )

    ! Description:
    ! Corrects the value of xm if w'x' needed to be clipped, for xm is partially
    ! based on the derivative of w'x' with respect to altitude.
    !
    ! The time-tendency equation for xm is:
    !
    ! d(xm)/dt = -w d(xm)/dz - d(w'x')/dz + d(xm)/dt|_ls;
    !
    ! where d(xm)/dt|_ls is the rate of change of xm over time due to radiation,
    ! microphysics, and/or any other large-scale forcing(s).
    !
    ! The time-tendency equation for xm is solved in conjunction with the
    ! time-tendency equation for w'x'.  Both equations are solved together in a
    ! semi-implicit manner.  However, after both equations have been solved (and
    ! thus both xm and w'x' have been advanced to the next timestep with
    ! timestep index {t+1}), the value of covariance w'x' may be clipped at any
    ! level in order to prevent the correlation of w and x from becoming greater
    ! than 1 or less than -1.
    !
    ! The correlation between w and x is:
    !
    ! corr_(w,x) = w'x' / [ sqrt(w'^2) * sqrt(x'^2) ].
    !
    ! The correlation must always have a value between -1 and 1, such that:
    !
    ! -1 <= corr_(w,x) <= 1.
    !
    ! Therefore, there is an upper limit on w'x', such that:
    !
    ! w'x' <=  [ sqrt(w'^2) * sqrt(x'^2) ];
    !
    ! and a lower limit on w'x', such that:
    !
    ! w'x' >= -[ sqrt(w'^2) * sqrt(x'^2) ].
    !
    ! The aforementioned time-tendency equation for xm is based on the value of
    ! w'x' without being clipped (w'x'{t+1}_unclipped), such that:
    !
    ! d(xm)/dt = -w d(xm{t+1})/dz - d(w'x'{t+1}_unclipped)/dz + d(xm{t})/dt|_ls;
    !
    ! where the both the mean advection term, -w d(xm{t+1})/dz, and the
    ! turbulent advection term, -d(w'x'{t+1}_unclipped)/dz, are solved
    ! completely implicitly.  The xm forcing term, +d(xm{t})/dt|_ls, is solved
    ! completely explicitly.
    !
    ! However, if w'x' needs to be clipped after being advanced one timestep,
    ! then xm needs to be altered to reflect the fact that w'x' has a different
    ! value than the value used while both were being solved together.  Ideally,
    ! the xm time-tendency equation that should be used is:
    !
    ! d(xm)/dt = -w d(xm{t+1})/dz - d(w'x'{t+1}_clipped)/dz + d(xm{t})/dt|_ls.
    !
    ! However, w'x'{t+1}_clipped isn't known until after the w'x' and xm
    ! equations have been solved together.  However, a proper adjuster can be
    ! applied to xm through the use of the following relationship:
    !
    ! w'x'{t+1}_clipped = w'x'{t+1}_unclipped + w'x'{t+1}_amount_clipped;
    !
    ! at any given vertical level.
    !
    ! When the expression above is substituted into the preceeding xm
    ! time-tendency equation, the resulting equation for xm time-tendency is:
    !
    ! d(xm)/dt = -w d(xm{t+1})/dz - d(w'x'{t+1}_unclipped)/dz
    !               - d(w'x'{t+1}_amount_clipped)/dz + d(xm{t})/dt|_ls.
    !
    ! Thus, the resulting xm time-tendency equation is the same as the original
    ! xm time-tendency equation, but with added adjuster term:
    !
    ! -d(w'x'{t+1}_amount_clipped)/dz.
    !
    ! Since the adjuster term needs to be applied after xm has already been
    ! solved, it needs to be multiplied by the timestep length and added on to
    ! xm{t+1}, such that:
    !
    ! xm{t+1}_after_adjustment =
    !    xm{t+1}_before_adjustment + ( -d(w'x'{t+1}_amount_clipped)/dz ) * dt.
    !
    ! The adjuster term is discretized as follows:
    !
    ! The values of w'x' are located on the momentum levels.  Thus, the values
    ! of w'x'_amount_clipped are also located on the momentum levels.  The
    ! values of xm are located on the thermodynamic levels.  The derivatives
    ! (d/dz) of w'x'_amount_clipped are taken over the intermediate
    ! thermodynamic levels, where they are applied to xm.
    !
    ! =======wpxp_amount_clipped=============================== m(k)
    !
    ! -----------------------------d(wpxp_amount_clipped)/dz--- t(k)
    !
    ! =======wpxpm1_amount_clipped============================= m(k-1)
    !
    ! The vertical indices m(k), t(k), and m(k-1) correspond with altitudes
    ! zm(k), zt(k), and zm(k-1), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! dzt(k) = 1 / ( zm(k) - zm(k-1) )

    ! Note:  The results of this xm adjustment are highly dependent on the
    !        numerical stability and the smoothness of the w'^2 and x'^2 fields.
    !        An unstable "sawtooth" profile for w'^2 and/or x'^2 causes an
    !        unstable "sawtooth" profile for the upper and lower limits on w'x'.
    !        In turn, this causes an unstable "sawtooth" profile for
    !        w'x'_amount_clipped.  Taking the derivative of that such a "noisy"
    !        field and applying the results to xm causes the xm field to become
    !        more "noisy" and unstable.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr  ! Variable(s); gr%nnzp only.

    use stats_precision, only: &
        time_precision

    use stats_type, only: &
        stat_update_var ! Procedure(s)

    use stats_variables, only: &
        l_stats_samp, & ! Variable(s)
        zt, &
        ithlm_tacl, &
        irtm_tacl

    implicit none

    ! Input Variables
    character(len=*), intent(in) :: &
      solve_type    ! Variable that is being solved for.

    real(kind=time_precision), intent(in) :: &
      dt            ! Model timestep                            [s]

    real, dimension(gr%nnzp), intent(in) :: &
      wpxp_chnge, & ! Amount of change in w'x' due to clipping  [m/s {xm units}]
      dzt           ! Inverse of grid spacing                   [1/m]

    ! Input/Output Variable
    real, dimension(gr%nnzp), intent(inout) :: &
      xm            ! xm (thermodynamic levels)                 [{xm units}]

    ! Local Variables
    real, dimension(gr%nnzp) :: &
      xm_tndcy_wpxp_cl ! d(xm)/dt due to clipping of w'x'       [{xm units}/s]

    integer :: k    ! Array index

    integer :: ixm_tacl  ! Statistical index


    select case ( trim( solve_type ) )
    case ( "rtm" )
      ixm_tacl = irtm_tacl
    case ( "thlm" )
      ixm_tacl = ithlm_tacl
    case default
      ixm_tacl = 0
    end select

    ! Adjusting xm based on clipping for w'x'.
    ! Loop over all thermodynamic levels between the second-lowest and the
    ! highest.
    do k = 2, gr%nnzp, 1
      xm_tndcy_wpxp_cl(k) = - dzt(k) * ( wpxp_chnge(k) - wpxp_chnge(k-1) )
      xm(k) = real( xm(k) + xm_tndcy_wpxp_cl(k) * dt )
    enddo

    if ( l_stats_samp ) then
      ! The adjustment to xm due to turbulent advection term clipping
      ! (xm term tacl) is completely explicit; call stat_update_var.
      call stat_update_var( ixm_tacl, xm_tndcy_wpxp_cl, zt )
    endif


    return

  end subroutine xm_correction_wpxp_cl

!===============================================================================

end module advance_xm_wpxp_module
