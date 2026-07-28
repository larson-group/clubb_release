!-------------------------------------------------------------------------------
module pdf_9_module

  ! PDF type 9 isolates a parcel-derived transport Gaussian from the CLUBB
  ! w-rt covariance, then diagnoses an ordinary ADG1 pair from the exact
  ! residual moments.  Its directional parcel ledger is updated only after
  ! the current PDF has been integrated, so the newly diagnosed geometry is
  ! the candidate consumed by the next PDF call.

  implicit none

  private

  ! Keep this experiment isolated from the parcel launch, buoyancy reference,
  ! and stopping rules.  Setting it false restores grid-mean thermodynamic
  ! entrainment exactly.
  logical, parameter :: &
    l_pdf9_weighted_g3_entrainment = .true.

  public :: pdf_9_driver, &
            diagnose_pdf9_mixing_reach, &
            compute_upward_mix, &
            compute_downward_mix

  contains

  !-----------------------------------------------------------------------------
  subroutine pdf_9_driver(                                           &
      nz, ngrdcol, sclr_dim, sclr_tol,                                 & ! In
      wm, rtm, thlm, um, vm,                                           & ! In
      wp2, wp3, rtp2, thlp2, up2, vp2,                                & ! In
      Skw, wprtp, wpthlp, rtpthlp, upwp, vpwp,                        & ! In
      sigma_sqd_w, beta, mixt_frac_max_mag,                            & ! In
      sclrm, sclrp2, wpsclrp, l_scalar_calc, err_info,                 & ! InOut
      w_1, w_2, rt_1, rt_2, thl_1, thl_2,                             & ! Out
      u_1, u_2, v_1, v_2,                                             & ! Out
      varnce_w_1, varnce_w_2, varnce_rt_1, varnce_rt_2,               & ! Out
      varnce_thl_1, varnce_thl_2, varnce_u_1, varnce_u_2,             & ! Out
      varnce_v_1, varnce_v_2, mixt_frac,                              & ! Out
      alpha_rt, alpha_thl, alpha_u, alpha_v,                           & ! Out
      sclr_1, sclr_2, varnce_sclr_1, varnce_sclr_2, alpha_sclr,       & ! Out
      w_3, rt_3, thl_3,                                                & ! Out
      varnce_w_3, varnce_rt_3, varnce_thl_3,                          & ! Out
      corr_w_rt_3, corr_w_thl_3, corr_rt_thl_3, mixt_frac_3,          & ! Out
      rtpthlp_outer,                                                   & ! Out
      varnce_u_3, varnce_v_3, covar_w_u_3, covar_w_v_3,               & ! Out
      pdf9_candidate_valid )                                           ! In

    use adg1_adg2_3d_luhar_pdf, only: &
      ADG1_pdf_driver

    use clubb_precision, only: &
      core_rknd

    use constants_clubb, only: &
      one, &
      one_half, &
      rt_tol, &
      thl_tol, &
      w_tol_sqd, &
      zero, &
      zero_threshold

    use error_code, only: &
      clubb_fatal_error

    use err_info_type_module, only: &
      err_info_type

    implicit none

    integer, intent(in) :: &
      nz,       & ! Number of vertical levels
      ngrdcol,  & ! Number of grid columns
      sclr_dim    ! Number of passive scalars

    real( kind = core_rknd ), dimension(sclr_dim), intent(in) :: &
      sclr_tol

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      wm, rtm, thlm, um, vm,                                     &
      wp2, wp3, rtp2, thlp2, up2, vp2,                           &
      Skw, wprtp, wpthlp, rtpthlp, upwp, vpwp, sigma_sqd_w

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      beta

    real( kind = core_rknd ), intent(in) :: &
      mixt_frac_max_mag

    real( kind = core_rknd ), dimension(ngrdcol,nz,sclr_dim), intent(in) :: &
      sclrm, sclrp2, wpsclrp

    logical, intent(in) :: &
      l_scalar_calc

    type(err_info_type), intent(inout) :: &
      err_info

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      w_1, w_2, rt_1, rt_2, thl_1, thl_2,                       &
      u_1, u_2, v_1, v_2,                                       &
      varnce_w_1, varnce_w_2, varnce_rt_1, varnce_rt_2,         &
      varnce_thl_1, varnce_thl_2, varnce_u_1, varnce_u_2,       &
      varnce_v_1, varnce_v_2, mixt_frac,                        &
      alpha_rt, alpha_thl, alpha_u, alpha_v,                     &
      mixt_frac_3,                                              &
      rtpthlp_outer,                                            &
      varnce_u_3, varnce_v_3, covar_w_u_3, covar_w_v_3

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(inout) :: &
      w_3, rt_3, thl_3,                                         &
      varnce_w_3, varnce_rt_3, varnce_thl_3,                    &
      corr_w_rt_3, corr_w_thl_3, corr_rt_thl_3

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      pdf9_candidate_valid

    real( kind = core_rknd ), dimension(ngrdcol,nz,sclr_dim), intent(out) :: &
      sclr_1, sclr_2, varnce_sclr_1, varnce_sclr_2, alpha_sclr

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      residual_wm, residual_rtm, residual_thlm,                  &
      wp2_outer, rtp2_outer, thlp2_outer,                        &
      up2_outer, vp2_outer, wprtp_outer, wpthlp_outer,           &
      upwp_outer, vpwp_outer, Skw_outer, sqrt_wp2_outer

    real( kind = core_rknd ), dimension(3) :: &
      grid_mean, packet_mean, residual_mean

    real( kind = core_rknd ), dimension(3,3) :: &
      grid_covariance, packet_covariance, residual_covariance

    real( kind = core_rknd ) :: &
      requested_weight, & ! Analytic G3 weight before feasibility limiting [-]
      packet_covar_w_rt, & ! Within-G3 w-rt covariance [(m/s)(kg/kg)]
      residual_weight      ! Weight of the residual ADG1 pair [-]

    real( kind = core_rknd ), parameter :: &
      pdf9_weight_max = 0.35_core_rknd, & ! Conservative G3 fraction cap [-]
      packet_weight_min = 1.0e-10_core_rknd ! Numerically active G3 weight [-]

    logical :: &
      valid_residual

    integer :: &
      i, k, shrink_count

    ! Preserve a finite candidate geometry from the preceding post-closure
    ! parcel diagnosis, then diagnose its analytic occupancy from w'r_t'.
    ! The residual pair is reconstructed exactly before ADG1 is called.
    do k = 1, nz
      do i = 1, ngrdcol
        if ( pdf9_candidate_valid(i,k) > one_half ) then
          varnce_w_3(i,k) = max( varnce_w_3(i,k), zero_threshold )
          varnce_rt_3(i,k) = max( varnce_rt_3(i,k), zero_threshold**2 )
          varnce_thl_3(i,k) = max( varnce_thl_3(i,k), zero_threshold )
        else
          w_3(i,k) = wm(i,k)
          rt_3(i,k) = rtm(i,k)
          thl_3(i,k) = thlm(i,k)

          varnce_w_3(i,k) = max( wp2(i,k), w_tol_sqd )
          varnce_rt_3(i,k) = max( rtp2(i,k), rt_tol**2 )
          varnce_thl_3(i,k) = max( thlp2(i,k), thl_tol**2 )

          corr_w_rt_3(i,k) = zero
          corr_w_thl_3(i,k) = zero
          corr_rt_thl_3(i,k) = zero
        end if

        grid_mean = [ wm(i,k), rtm(i,k), thlm(i,k) ]
        packet_mean = [ w_3(i,k), rt_3(i,k), thl_3(i,k) ]
        grid_covariance(1,:) = [ wp2(i,k), wprtp(i,k), wpthlp(i,k) ]
        grid_covariance(2,:) = [ wprtp(i,k), rtp2(i,k), rtpthlp(i,k) ]
        grid_covariance(3,:) = [ wpthlp(i,k), rtpthlp(i,k), thlp2(i,k) ]

        packet_covariance = zero
        packet_covariance(1,1) = varnce_w_3(i,k)
        packet_covariance(2,2) = varnce_rt_3(i,k)
        packet_covariance(3,3) = varnce_thl_3(i,k)
        packet_covariance(1,2) = corr_w_rt_3(i,k) &
          * sqrt( varnce_w_3(i,k) * varnce_rt_3(i,k) )
        packet_covariance(2,1) = packet_covariance(1,2)
        packet_covariance(1,3) = corr_w_thl_3(i,k) &
          * sqrt( varnce_w_3(i,k) * varnce_thl_3(i,k) )
        packet_covariance(3,1) = packet_covariance(1,3)
        packet_covariance(2,3) = corr_rt_thl_3(i,k) &
          * sqrt( varnce_rt_3(i,k) * varnce_thl_3(i,k) )
        packet_covariance(3,2) = packet_covariance(2,3)

        packet_covar_w_rt = packet_covariance(1,2)
        if ( pdf9_candidate_valid(i,k) > one_half .and. .not. l_scalar_calc ) then
          requested_weight = diagnose_pdf9_weight(                        &
              wprtp(i,k), packet_covar_w_rt,                               & ! In
              ( w_3(i,k) - wm(i,k) ) * ( rt_3(i,k) - rtm(i,k) ),           & ! In
              pdf9_weight_max )                                             ! In
        else
          requested_weight = zero
        end if

        valid_residual = .false.
        do shrink_count = 0, 40
          call remove_packet_moments( requested_weight,                    & ! In
                                      grid_mean, grid_covariance, wp3(i,k), & ! In
                                      packet_mean, packet_covariance,       & ! In
                                      valid_residual, residual_mean,        & ! Out
                                      residual_covariance, Skw_outer(i,k) )   ! Out
          if ( valid_residual ) then
            valid_residual = residual_covariance(1,1) >= w_tol_sqd &
                             .and. residual_covariance(2,2) >= rt_tol**2 &
                             .and. residual_covariance(3,3) >= thl_tol**2
          end if
          if ( valid_residual ) exit
          requested_weight = one_half * requested_weight
          if ( requested_weight <= packet_weight_min ) requested_weight = zero
        end do

        if ( .not. valid_residual ) then
          requested_weight = zero
          call remove_packet_moments( requested_weight,                    & ! In
                                      grid_mean, grid_covariance, wp3(i,k), & ! In
                                      packet_mean, packet_covariance,       & ! In
                                      valid_residual, residual_mean,        & ! Out
                                      residual_covariance, Skw_outer(i,k) )   ! Out
        end if
        if ( .not. valid_residual ) then
          err_info%err_code(i) = clubb_fatal_error
          cycle
        end if

        mixt_frac_3(i,k) = requested_weight
        residual_weight = one - requested_weight
        residual_wm(i,k) = residual_mean(1)
        residual_rtm(i,k) = residual_mean(2)
        residual_thlm(i,k) = residual_mean(3)
        wp2_outer(i,k) = max( residual_covariance(1,1), w_tol_sqd )
        rtp2_outer(i,k) = max( residual_covariance(2,2), rt_tol**2 )
        thlp2_outer(i,k) = max( residual_covariance(3,3), thl_tol**2 )
        wprtp_outer(i,k) = residual_covariance(1,2)
        wpthlp_outer(i,k) = residual_covariance(1,3)
        rtpthlp_outer(i,k) = residual_covariance(2,3)
        if ( requested_weight <= packet_weight_min ) Skw_outer(i,k) = Skw(i,k)

        varnce_u_3(i,k) = max( up2(i,k), zero )
        varnce_v_3(i,k) = max( vp2(i,k), zero )
        if ( wp2(i,k) > tiny(one) ) then
          covar_w_u_3(i,k) = upwp(i,k) * varnce_w_3(i,k) / wp2(i,k)
          covar_w_v_3(i,k) = vpwp(i,k) * varnce_w_3(i,k) / wp2(i,k)
        else
          covar_w_u_3(i,k) = zero
          covar_w_v_3(i,k) = zero
        end if
        covar_w_u_3(i,k) = clip_covariance( covar_w_u_3(i,k), &
                                            varnce_w_3(i,k), varnce_u_3(i,k) )
        covar_w_v_3(i,k) = clip_covariance( covar_w_v_3(i,k), &
                                            varnce_w_3(i,k), varnce_v_3(i,k) )

        up2_outer(i,k) = max( ( up2(i,k) - requested_weight * varnce_u_3(i,k) ) &
                              / max( residual_weight, tiny(one) ), zero )
        vp2_outer(i,k) = max( ( vp2(i,k) - requested_weight * varnce_v_3(i,k) ) &
                              / max( residual_weight, tiny(one) ), zero )
        upwp_outer(i,k) = ( upwp(i,k) - requested_weight * covar_w_u_3(i,k) ) &
                          / max( residual_weight, tiny(one) )
        vpwp_outer(i,k) = ( vpwp(i,k) - requested_weight * covar_w_v_3(i,k) ) &
                          / max( residual_weight, tiny(one) )
        upwp_outer(i,k) = clip_covariance( upwp_outer(i,k), &
                                           wp2_outer(i,k), up2_outer(i,k) )
        vpwp_outer(i,k) = clip_covariance( vpwp_outer(i,k), &
                                           wp2_outer(i,k), vp2_outer(i,k) )
        sqrt_wp2_outer(i,k) = sqrt( wp2_outer(i,k) )
      end do
    end do

    call ADG1_pdf_driver( nz, ngrdcol, sclr_dim, sclr_tol,             & ! In
                          residual_wm, residual_rtm, residual_thlm, um, vm, & ! In
                          wp2_outer, rtp2_outer, thlp2_outer,          & ! In
                          up2_outer, vp2_outer, Skw_outer,             & ! In
                          wprtp_outer, wpthlp_outer,                   & ! In
                          upwp_outer, vpwp_outer, sqrt_wp2_outer,      & ! In
                          sigma_sqd_w, beta, mixt_frac_max_mag,        & ! In
                          sclrm, sclrp2, wpsclrp, l_scalar_calc,       & ! In
                          err_info,                                   & ! InOut
                          w_1, w_2, rt_1, rt_2, thl_1, thl_2,         & ! Out
                          u_1, u_2, v_1, v_2,                         & ! Out
                          varnce_w_1, varnce_w_2,                     & ! Out
                          varnce_rt_1, varnce_rt_2,                   & ! Out
                          varnce_thl_1, varnce_thl_2,                 & ! Out
                          varnce_u_1, varnce_u_2,                     & ! Out
                          varnce_v_1, varnce_v_2, mixt_frac,          & ! Out
                          alpha_rt, alpha_thl, alpha_u, alpha_v,      & ! Out
                          sclr_1, sclr_2, varnce_sclr_1,              & ! Out
                          varnce_sclr_2, alpha_sclr )                   ! Out

    if ( any(err_info%err_code == clubb_fatal_error) ) then
      return
    end if

  end subroutine pdf_9_driver

  !=============================================================================

  function diagnose_pdf9_weight( grid_covar_w_rt, packet_covar_w_rt, &
                                  center_product, weight_max ) result( weight )

    ! Description:
    !   Finds the smaller physical G3 weight for which the residual two-
    !   Gaussian population has zero aggregate w-rt covariance:
    !
    !     C = a c_3 + a D / (1-a),
    !
    !   where C is grid w'r_t', c_3 is G3's internal covariance, and D is
    !   (w_3-wbar)(rt_3-rtbar).  No nonlinear fit or iteration is involved.

    use clubb_precision, only: &
      core_rknd

    use constants_clubb, only: &
      one, &
      zero

    use, intrinsic :: ieee_arithmetic, only: &
      ieee_is_finite

    implicit none

    real( kind = core_rknd ), intent(in) :: &
      grid_covar_w_rt, &
      packet_covar_w_rt, &
      center_product, &
      weight_max

    real( kind = core_rknd ) :: &
      weight

    real( kind = core_rknd ) :: &
      coefficient, discriminant, root_1, root_2, scale

    weight = zero
    if ( .not. ieee_is_finite(grid_covar_w_rt) &
         .or. .not. ieee_is_finite(packet_covar_w_rt) &
         .or. .not. ieee_is_finite(center_product) ) return

    scale = max( abs(grid_covar_w_rt), abs(packet_covar_w_rt), &
                 abs(center_product), tiny(one) )
    if ( abs(grid_covar_w_rt) <= 1.0e-10_core_rknd * scale ) return
    if ( grid_covar_w_rt * ( packet_covar_w_rt + center_product ) <= zero ) return

    coefficient = grid_covar_w_rt + packet_covar_w_rt + center_product
    if ( abs(packet_covar_w_rt) <= 1.0e-10_core_rknd * scale ) then
      if ( abs(coefficient) > tiny(one) ) weight = grid_covar_w_rt / coefficient
    else
      discriminant = coefficient**2 &
                     - 4.0_core_rknd * packet_covar_w_rt * grid_covar_w_rt
      if ( discriminant < zero ) return
      discriminant = sqrt( max( discriminant, zero ) )
      root_1 = ( coefficient - discriminant ) &
               / ( 2.0_core_rknd * packet_covar_w_rt )
      root_2 = ( coefficient + discriminant ) &
               / ( 2.0_core_rknd * packet_covar_w_rt )
      if ( root_1 > zero .and. root_1 < one ) weight = root_1
      if ( root_2 > zero .and. root_2 < one ) then
        if ( weight <= zero .or. root_2 < weight ) weight = root_2
      end if
    end if

    if ( .not. ieee_is_finite(weight) .or. weight <= zero ) then
      weight = zero
    else
      weight = min( weight, weight_max )
    end if

    return

  end function diagnose_pdf9_weight

  !=============================================================================

  subroutine remove_packet_moments( packet_weight, grid_mean, grid_covariance, & ! In
                                    grid_wp3, packet_mean, packet_covariance,   & ! In
                                    valid, residual_mean, residual_covariance, & ! Out
                                    residual_skewness )                          ! Out

    ! Description:
    !   Removes G3 exactly from the grid mean, covariance, and w third moment,
    !   then verifies that the residual trivariate covariance is realizable.

    use clubb_precision, only: &
      core_rknd

    use constants_clubb, only: &
      one, &
      zero

    use, intrinsic :: ieee_arithmetic, only: &
      ieee_is_finite

    implicit none

    real( kind = core_rknd ), intent(in) :: &
      packet_weight, &
      grid_wp3

    real( kind = core_rknd ), dimension(3), intent(in) :: &
      grid_mean, packet_mean

    real( kind = core_rknd ), dimension(3,3), intent(in) :: &
      grid_covariance, packet_covariance

    logical, intent(out) :: &
      valid

    real( kind = core_rknd ), dimension(3), intent(out) :: &
      residual_mean

    real( kind = core_rknd ), dimension(3,3), intent(out) :: &
      residual_covariance

    real( kind = core_rknd ), intent(out) :: &
      residual_skewness

    integer :: &
      j, l

    real( kind = core_rknd ) :: &
      residual_weight, packet_dw, residual_dw, packet_wp3, residual_wp3, &
      scale, determinant

    real( kind = core_rknd ), parameter :: &
      packet_weight_min = 1.0e-10_core_rknd

    real( kind = core_rknd ), dimension(3,3) :: &
      grid_raw_second, packet_raw_second, residual_raw_second

    valid = .false.
    residual_mean = zero
    residual_covariance = zero
    residual_skewness = zero
    if ( .not. ieee_is_finite(packet_weight) &
         .or. .not. ieee_is_finite(grid_wp3) &
         .or. .not. all(ieee_is_finite(grid_mean)) &
         .or. .not. all(ieee_is_finite(grid_covariance)) ) return
    residual_weight = one - packet_weight
    if ( packet_weight < zero .or. residual_weight <= 0.5_core_rknd ) return

    if ( packet_weight <= packet_weight_min ) then
      residual_mean = grid_mean
      residual_covariance = 0.5_core_rknd &
                            * ( grid_covariance + transpose(grid_covariance) )
    else
      if ( .not. all(ieee_is_finite(packet_mean)) &
           .or. .not. all(ieee_is_finite(packet_covariance)) ) return
      do l = 1, 3
        do j = 1, 3
          grid_raw_second(j,l) = grid_covariance(j,l) &
                                 + grid_mean(j) * grid_mean(l)
          packet_raw_second(j,l) = packet_covariance(j,l) &
                                   + packet_mean(j) * packet_mean(l)
          residual_raw_second(j,l) = ( grid_raw_second(j,l) &
                                      - packet_weight * packet_raw_second(j,l) ) &
                                      / residual_weight
        end do
      end do
      residual_mean = ( grid_mean - packet_weight * packet_mean ) / residual_weight
      do l = 1, 3
        do j = 1, 3
          residual_covariance(j,l) = residual_raw_second(j,l) &
                                     - residual_mean(j) * residual_mean(l)
        end do
      end do
      residual_covariance = 0.5_core_rknd &
                            * ( residual_covariance + transpose(residual_covariance) )
    end if

    if ( .not. all(ieee_is_finite(residual_mean)) &
         .or. .not. all(ieee_is_finite(residual_covariance)) ) return
    scale = max( maxval( abs( [ residual_covariance(1,1), &
                                residual_covariance(2,2), &
                                residual_covariance(3,3) ] ) ), tiny(one) )
    determinant = covariance_determinant( residual_covariance )
    if ( minval( [ residual_covariance(1,1), residual_covariance(2,2), &
                   residual_covariance(3,3) ] ) < -1.0e-12_core_rknd * scale ) return
    if ( residual_covariance(1,1) * residual_covariance(2,2) &
         - residual_covariance(1,2)**2 < -1.0e-12_core_rknd * scale**2 ) return
    if ( residual_covariance(1,1) * residual_covariance(3,3) &
         - residual_covariance(1,3)**2 < -1.0e-12_core_rknd * scale**2 ) return
    if ( residual_covariance(2,2) * residual_covariance(3,3) &
         - residual_covariance(2,3)**2 < -1.0e-12_core_rknd * scale**2 ) return
    if ( packet_weight > packet_weight_min &
         .and. determinant < -1.0e-12_core_rknd * scale**3 ) return

    residual_covariance(1,1) = max( residual_covariance(1,1), zero )
    residual_covariance(2,2) = max( residual_covariance(2,2), zero )
    residual_covariance(3,3) = max( residual_covariance(3,3), zero )

    if ( packet_weight > packet_weight_min ) then
      packet_dw = packet_mean(1) - grid_mean(1)
      packet_wp3 = packet_dw**3 &
                   + 3.0_core_rknd * packet_dw * packet_covariance(1,1)
    else
      packet_wp3 = zero
    end if
    residual_dw = residual_mean(1) - grid_mean(1)
    residual_wp3 = ( grid_wp3 - packet_weight * packet_wp3 ) / residual_weight &
                   - residual_dw**3 &
                   - 3.0_core_rknd * residual_dw * residual_covariance(1,1)
    residual_skewness = residual_wp3 &
                        / max( residual_covariance(1,1), tiny(one) )**1.5_core_rknd
    if ( .not. ieee_is_finite(residual_skewness) ) return
    valid = .true.

    return

  end subroutine remove_packet_moments

  !=============================================================================

  function covariance_determinant( covariance ) result( determinant )

    use clubb_precision, only: &
      core_rknd

    implicit none

    real( kind = core_rknd ), dimension(3,3), intent(in) :: &
      covariance

    real( kind = core_rknd ) :: &
      determinant

    determinant = covariance(1,1) &
                  * ( covariance(2,2) * covariance(3,3) - covariance(2,3)**2 ) &
                  - covariance(1,2) &
                    * ( covariance(1,2) * covariance(3,3) &
                        - covariance(1,3) * covariance(2,3) ) &
                  + covariance(1,3) &
                    * ( covariance(1,2) * covariance(2,3) &
                        - covariance(1,3) * covariance(2,2) )

    return

  end function covariance_determinant

  !=============================================================================

  function clip_covariance( covariance, variance_x, variance_y ) result( clipped )

    use clubb_precision, only: &
      core_rknd

    use constants_clubb, only: &
      zero

    implicit none

    real( kind = core_rknd ), intent(in) :: &
      covariance, variance_x, variance_y

    real( kind = core_rknd ) :: &
      clipped, limit

    limit = sqrt( max( variance_x * variance_y, zero ) )
    clipped = min( max( covariance, -limit ), limit )

    return

  end function clip_covariance

  !=============================================================================

  subroutine compute_upward_mix(                                      &
      nzm, nzt, ngrdcol, gr, l_implemented,                           & ! In
      p_in_Pa, exner, entrain_rtm, entrain_thlm,                      & ! In
      thvm, thv_ds,                                                   & ! In
      parcel_rt_init, parcel_thl_init, parcel_tke_init,                & ! In
      mu, lmin, saturation_formula,                                   & ! In
      err_info,                                                       & ! InOut
      pdf9_lscale_up,                                                 & ! Out
      parcel_tke, parcel_w, parcel_rt, parcel_thl,                    & ! Out
      parcel_buoyancy, parcel_status,                                 & ! Out
      crossing_weight, crossing_mean_w, crossing_mean_rt,            & ! Out
      crossing_mean_thl, crossing_var_w, crossing_var_rt,            & ! Out
      crossing_var_thl, crossing_covar_w_rt, crossing_covar_w_thl,   & ! Out
      crossing_covar_rt_thl )                                          ! Out

    ! Description:
    !   Duplicates the upward half of CLUBB's historical moist, nonlocal
    !   parcel-tracking mixing-length calculation.  The only intentional
    !   difference is that PDF 9 supplies the parcel's launch rt, thl, and
    !   kinetic energy explicitly.  This routine is deliberately independent
    !   so that its parcel trajectory can be extended during PDF-9 prototyping.

    use clubb_precision, only: &
      core_rknd

    use constants_clubb, only: &
      Cp, &
      Rd, &
      ep, &
      ep1, &
      ep2, &
      Lv, &
      grav, &
      fstderr, &
      zero_threshold, &
      eps, &
      one_half, &
      one, &
      two, &
      zero

    use error_code, only: &
      clubb_fatal_error

    use err_info_type_module, only: &
      err_info_type

    use grid_class, only: &
      grid

    use saturation, only: &
      sat_mixrat_liq_api

    implicit none

    real( kind = core_rknd ), parameter :: &
      zlmin = 0.1_core_rknd, & ! Absolute minimum upward reach [m]
      Lscale_sfclyr_depth = 500._core_rknd ! Surface-layer floor depth [m]

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      nzm, & ! Number of momentum levels [-]
      nzt, & ! Number of thermodynamic levels [-]
      ngrdcol, & ! Number of grid columns [-]
      saturation_formula ! Saturation formula selector [-]

    type(grid), intent(in) :: &
      gr ! CLUBB vertical grid

    logical, intent(in) :: &
      l_implemented ! True when CLUBB is embedded in a host model [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      p_in_Pa, & ! Pressure [Pa]
      exner, & ! Exner function [-]
      entrain_rtm, & ! Environmental total-water entrainment target [kg/kg]
      entrain_thlm, & ! Environmental theta-l entrainment target [K]
      thvm, & ! Grid-mean virtual potential temperature [K]
      thv_ds, & ! Dry base-state virtual potential temperature [K]
      parcel_rt_init, & ! Parcel total water at launch [kg/kg]
      parcel_thl_init, & ! Parcel liquid-water potential temperature at launch [K]
      parcel_tke_init ! Parcel launch kinetic energy [m^2/s^2]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      mu ! Fractional entrainment rate [1/m]

    real( kind = core_rknd ), intent(in) :: &
      lmin ! Surface-layer minimum mixing length [m]

    !----------------------- Input/Output Variables ------------------------
    type(err_info_type), intent(inout) :: &
      err_info ! CLUBB error state

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      pdf9_lscale_up, & ! PDF-9 upward parcel reach [m]
      crossing_weight, & ! Sum of positive-w crossing weights [m/s]
      crossing_mean_w, & ! Crossing-weighted vertical velocity [m/s]
      crossing_mean_rt, & ! Crossing-weighted total water [kg/kg]
      crossing_mean_thl, & ! Crossing-weighted liquid-water potential temperature [K]
      crossing_var_w, & ! Crossing-weighted w variance [m2/s2]
      crossing_var_rt, & ! Crossing-weighted rt variance [kg2/kg2]
      crossing_var_thl, & ! Crossing-weighted thl variance [K2]
      crossing_covar_w_rt, & ! Crossing-weighted covariance of w and rt [m/s kg/kg]
      crossing_covar_w_thl, & ! Crossing-weighted covariance of w and thl [K m/s]
      crossing_covar_rt_thl ! Crossing-weighted covariance of rt and thl [K kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,nzt), intent(out) :: &
      parcel_tke, & ! Candidate parcel energy at (launch,destination) [m2/s2]
      parcel_w, & ! Positive parcel speed at (launch,destination) [m/s]
      parcel_rt, & ! Parcel total water at (launch,destination) [kg/kg]
      parcel_thl, & ! Parcel liquid-water potential temperature [K]
      parcel_buoyancy, & ! CAPE gradient at destination [m/s2]
      parcel_status ! 1 reached, -1 first failed crossing, 0 unvisited [-]

    !---------------------------- Local Variables --------------------------
    integer :: &
      i, & ! Grid-column index [-]
      j, & ! Parcel destination-level index [-]
      k, & ! Parcel launch-level index [-]
      start_index, & ! First level used by saturation adjustment [-]
      j_zm, & ! Momentum-level index associated with j [-]
      kp1_zm ! Momentum-level index immediately above k [-]

    integer :: &
      step ! Number of traversed layers from a launch level [-]

    real( kind = core_rknd ) :: &
      tke, & ! Remaining parcel kinetic energy [m^2/s^2]
      CAPE_incr, & ! CAPE increment across one layer [m^2/s^2]
      dCAPE_dz_j, & ! CAPE gradient at level j [m/s^2]
      dCAPE_dz_j_minus_1, & ! CAPE gradient at the preceding level [m/s^2]
      lminh, & ! Height-dependent minimum reach [m]
      thl_par_j, & ! Parcel thl at level j [K]
      rt_par_j, & ! Parcel rt at level j [kg/kg]
      rc_par_j, & ! Parcel cloud water at level j [kg/kg]
      thv_par_j, & ! Parcel virtual potential temperature at level j [K]
      tl_par_j, & ! Parcel liquid-water temperature at level j [K]
      rsatl_par_j, & ! Parcel saturation mixing ratio at level j [kg/kg]
      s_par_j, & ! Parcel extended liquid water at level j [kg/kg]
      Lscale_up_max_alt, & ! Highest reachable altitude from lower launches [m]
      Lv2_coef, & ! Precomputed latent-heating coefficient [K^2]
      tl_par_j_sqd, & ! Squared parcel temperature [K^2]
      invrs_dCAPE_diff, & ! Inverse CAPE-gradient difference [s^2/m]
      invrs_Lscale_sfclyr_depth ! Inverse surface-layer depth [1/m]

    real( kind = core_rknd ) :: &
      candidate_tke, & ! Energy after the candidate layer traversal [m2/s2]
      crossing_wgt, & ! Positive-w proxy weight for crossing statistics [m/s]
      delta_w, & ! Difference from destination crossing mean [m/s]
      delta_rt, & ! Difference from destination crossing mean [kg/kg]
      delta_thl ! Difference from destination crossing mean [K]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      grav_on_thvm, & ! Gravity divided by environmental thv [m/s^2/K]
      Lv_coef, & ! Latent-heating coefficient [K]
      thl_par_j_precalc, & ! Recursive thl environment contribution [K]
      rt_par_j_precalc, & ! Recursive rt environment contribution [kg/kg]
      tl_par_1, & ! Parcel temperature after the first layer [K]
      rt_par_1, & ! Parcel rt after the first layer [kg/kg]
      rsatl_par_1, & ! Parcel saturation mixing ratio after the first layer [kg/kg]
      thl_par_1, & ! Parcel thl after the first layer [K]
      dCAPE_dz_1, & ! Initial CAPE gradient [m/s^2]
      s_par_1, & ! Initial parcel extended liquid water [kg/kg]
      rc_par_1, & ! Initial parcel cloud water [kg/kg]
      CAPE_incr_1, & ! Initial CAPE increment [m^2/s^2]
      thv_par_1, & ! Initial parcel virtual potential temperature [K]
      tke_i, & ! Parcel launch kinetic energy [m^2/s^2]
      tke_state, & ! Energy of the last successfully reached level [m2/s2]
      thl_state, & ! thl at the last successfully reached level [K]
      rt_state, & ! rt at the last successfully reached level [kg/kg]
      dCAPE_state ! CAPE gradient at the last successfully reached level [m/s2]

    logical, dimension(ngrdcol,nzt) :: &
      parcel_active ! True while a launched parcel can traverse another layer

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      exp_mu_dzm, & ! Layer entrainment exponential [-]
      invrs_dzm_on_mu, & ! Layer integration coefficient [-]
      entrain_coef ! Linear-environment entrainment coefficient [-]

    !--------------------------- Begin Code ---------------------------

    do i = 1, ngrdcol
      if ( abs(mu(i)) < eps ) then
        err_info%err_code(i) = clubb_fatal_error
        write(fstderr,*) err_info%err_header(i)
        write(fstderr,*) "mu = ", mu(i)
      end if
    end do

    if ( any(err_info%err_code == clubb_fatal_error) ) then
      write(fstderr,*) "Entrainment rate mu cannot be 0"
      write(fstderr,*) "Fatal error in subroutine compute_upward_mix"
      return
    end if

    do k = 1, nzt
      do i = 1, ngrdcol
        tke_i(i,k) = parcel_tke_init(i,k)
        pdf9_lscale_up(i,k) = zlmin
        grav_on_thvm(i,k) = grav / thvm(i,k)
        Lv_coef(i,k) = Lv / ( exner(i,k) * Cp ) - ep2 * thv_ds(i,k)
      end do
    end do

    do i = 1, ngrdcol
      do k = 1, nzm
        exp_mu_dzm(i,k) = exp( -mu(i) * gr%grid_dir * gr%dzm(i,k) )
        invrs_dzm_on_mu(i,k) = ( gr%grid_dir * gr%invrs_dzm(i,k) ) / mu(i)
        entrain_coef(i,k) = ( one - exp_mu_dzm(i,k) ) * invrs_dzm_on_mu(i,k)
      end do
    end do

    Lv2_coef = ep * Lv**2 / ( Rd * Cp )
    invrs_Lscale_sfclyr_depth = one / Lscale_sfclyr_depth

    ! Precompute the affine environmental contribution for every traversed
    ! layer, exactly as in the historical parcel-at-a-time implementation.
    do j = gr%k_lb_zt+gr%grid_dir_indx, gr%k_ub_zt-gr%grid_dir_indx, gr%grid_dir_indx
      do i = 1, ngrdcol
        if ( gr%grid_dir_indx > 0 ) then
          j_zm = j
        else
          j_zm = j+1
        end if

        thl_par_j_precalc(i,j) = entrain_thlm(i,j) &
          - entrain_thlm(i,j-gr%grid_dir_indx) * exp_mu_dzm(i,j_zm) &
          - ( entrain_thlm(i,j) - entrain_thlm(i,j-gr%grid_dir_indx) ) &
            * entrain_coef(i,j_zm)
        rt_par_j_precalc(i,j) = entrain_rtm(i,j) &
          - entrain_rtm(i,j-gr%grid_dir_indx) * exp_mu_dzm(i,j_zm) &
          - ( entrain_rtm(i,j) - entrain_rtm(i,j-gr%grid_dir_indx) ) &
            * entrain_coef(i,j_zm)
      end do
    end do

    ! The first traversal retains the historical linearized launch formula.
    thl_par_1 = zero
    tl_par_1 = zero
    rt_par_1 = zero
    do j = gr%k_lb_zt+gr%grid_dir_indx, gr%k_ub_zt-gr%grid_dir_indx, gr%grid_dir_indx
      do i = 1, ngrdcol
        if ( gr%grid_dir_indx > 0 ) then
          j_zm = j
        else
          j_zm = j+1
        end if

        thl_par_1(i,j) = entrain_thlm(i,j) &
          - ( entrain_thlm(i,j) - parcel_thl_init(i,j-gr%grid_dir_indx) ) &
            * entrain_coef(i,j_zm)
        tl_par_1(i,j) = thl_par_1(i,j) * exner(i,j)
        rt_par_1(i,j) = entrain_rtm(i,j) &
          - ( entrain_rtm(i,j) - parcel_rt_init(i,j-gr%grid_dir_indx) ) &
            * entrain_coef(i,j_zm)
      end do
    end do

    start_index = gr%k_lb_zt+gr%grid_dir_indx
    rsatl_par_1 = sat_mixrat_liq_api( nzt, ngrdcol, gr, p_in_Pa, tl_par_1, &
                                      saturation_formula, start_index )

    do j = gr%k_lb_zt+gr%grid_dir_indx, gr%k_ub_zt, gr%grid_dir_indx
      do i = 1, ngrdcol
        if ( gr%grid_dir_indx > 0 ) then
          j_zm = j
        else
          j_zm = j+1
        end if

        tl_par_j_sqd = tl_par_1(i,j)**2
        s_par_1(i,j) = ( rt_par_1(i,j) - rsatl_par_1(i,j) ) * tl_par_j_sqd &
          / ( tl_par_j_sqd + Lv2_coef * rsatl_par_1(i,j) )
        rc_par_1(i,j) = max( s_par_1(i,j), zero_threshold )
        thv_par_1(i,j) = thl_par_1(i,j) + ep1 * thv_ds(i,j) * rt_par_1(i,j) &
          + Lv_coef(i,j) * rc_par_1(i,j)
        dCAPE_dz_1(i,j) = grav_on_thvm(i,j) * ( thv_par_1(i,j) - thvm(i,j) )
        CAPE_incr_1(i,j) = one_half * dCAPE_dz_1(i,j) &
          * gr%grid_dir * gr%dzm(i,j_zm)
      end do
    end do

    ! Initialize every launch together.  The third index of each history is
    ! the actual destination level, not merely a traversal counter.
    parcel_tke = zero
    parcel_w = zero
    parcel_rt = zero
    parcel_thl = zero
    parcel_buoyancy = zero
    parcel_status = zero
    crossing_weight = zero
    crossing_mean_w = zero
    crossing_mean_rt = zero
    crossing_mean_thl = zero
    crossing_var_w = zero
    crossing_var_rt = zero
    crossing_var_thl = zero
    crossing_covar_w_rt = zero
    crossing_covar_w_thl = zero
    crossing_covar_rt_thl = zero

    do k = 1, nzt
      do i = 1, ngrdcol
        parcel_active(i,k) = .false.
        tke_state(i,k) = tke_i(i,k)
        rt_state(i,k) = parcel_rt_init(i,k)
        thl_state(i,k) = parcel_thl_init(i,k)
        dCAPE_state(i,k) = zero

        parcel_tke(i,k,k) = tke_i(i,k)
        parcel_w(i,k,k) = sqrt( max( two * tke_i(i,k), zero ) )
        parcel_rt(i,k,k) = parcel_rt_init(i,k)
        parcel_thl(i,k,k) = parcel_thl_init(i,k)
        parcel_status(i,k,k) = one
      end do
    end do

    ! Advance all launch parcels by one layer at a time.  Each step is a flat
    ! (column,launch-level) operation, so every parcel writes only to its own
    ! history slab and no trajectory carries scalar temporaries across launches.
    do step = 1, nzt
      do k = gr%k_lb_zt, gr%k_ub_zt-2*gr%grid_dir_indx, gr%grid_dir_indx
        do i = 1, ngrdcol
          j = k + step * gr%grid_dir_indx

          if ( step == 1 ) then
            if ( gr%grid_dir_indx * j >= gr%grid_dir_indx * gr%k_ub_zt ) cycle

            candidate_tke = tke_i(i,k) + CAPE_incr_1(i,j)
            parcel_tke(i,k,j) = candidate_tke
            parcel_w(i,k,j) = sqrt( max( two * candidate_tke, zero ) )
            parcel_rt(i,k,j) = rt_par_1(i,j)
            parcel_thl(i,k,j) = thl_par_1(i,j)
            parcel_buoyancy(i,k,j) = dCAPE_dz_1(i,j)

            if ( candidate_tke > zero ) then
              parcel_status(i,k,j) = one
              parcel_active(i,k) = .true.
              tke_state(i,k) = candidate_tke
              rt_state(i,k) = rt_par_1(i,j)
              thl_state(i,k) = thl_par_1(i,j)
              dCAPE_state(i,k) = dCAPE_dz_1(i,j)
            else
              parcel_status(i,k,j) = -one
              parcel_active(i,k) = .false.
            end if

          else if ( parcel_active(i,k) ) then
            if ( gr%grid_dir_indx * j >= gr%grid_dir_indx * gr%k_ub_zt ) then
              parcel_active(i,k) = .false.
              cycle
            end if

            if ( gr%grid_dir_indx > 0 ) then
              j_zm = j
            else
              j_zm = j+1
            end if

            thl_par_j = thl_par_j_precalc(i,j) &
              + thl_state(i,k) * exp_mu_dzm(i,j_zm)
            rt_par_j = rt_par_j_precalc(i,j) &
              + rt_state(i,k) * exp_mu_dzm(i,j_zm)
            tl_par_j = thl_par_j * exner(i,j)
            rsatl_par_j = sat_mixrat_liq_api( p_in_Pa(i,j), tl_par_j, &
                                              saturation_formula )
            tl_par_j_sqd = tl_par_j**2
            s_par_j = ( rt_par_j - rsatl_par_j ) * tl_par_j_sqd &
              / ( tl_par_j_sqd + Lv2_coef * rsatl_par_j )
            rc_par_j = max( s_par_j, zero_threshold )
            thv_par_j = thl_par_j + ep1 * thv_ds(i,j) * rt_par_j &
              + Lv_coef(i,j) * rc_par_j
            dCAPE_dz_j = grav_on_thvm(i,j) * ( thv_par_j - thvm(i,j) )
            CAPE_incr = one_half * ( dCAPE_dz_j + dCAPE_state(i,k) ) &
              * gr%grid_dir * gr%dzm(i,j_zm)
            candidate_tke = tke_state(i,k) + CAPE_incr

            parcel_tke(i,k,j) = candidate_tke
            parcel_w(i,k,j) = sqrt( max( two * candidate_tke, zero ) )
            parcel_rt(i,k,j) = rt_par_j
            parcel_thl(i,k,j) = thl_par_j
            parcel_buoyancy(i,k,j) = dCAPE_dz_j

            if ( candidate_tke > zero ) then
              parcel_status(i,k,j) = one
              tke_state(i,k) = candidate_tke
              rt_state(i,k) = rt_par_j
              thl_state(i,k) = thl_par_j
              dCAPE_state(i,k) = dCAPE_dz_j
            else
              parcel_status(i,k,j) = -one
              parcel_active(i,k) = .false.
            end if
          end if
        end do
      end do
    end do

    ! Reconstruct the historical stopping distance from the trajectory ledger,
    ! including the first failed crossing needed for sub-layer interpolation.
    do i = 1, ngrdcol
      Lscale_up_max_alt = zero
      do k = gr%k_lb_zt, gr%k_ub_zt-2*gr%grid_dir_indx, gr%grid_dir_indx
        j = k + gr%grid_dir_indx

        if ( parcel_status(i,k,j) > zero ) then
          do while ( gr%grid_dir_indx * (j+gr%grid_dir_indx) &
                     < gr%grid_dir_indx * gr%k_ub_zt )
            if ( parcel_status(i,k,j+gr%grid_dir_indx) <= zero ) exit
            j = j + gr%grid_dir_indx
          end do

          pdf9_lscale_up(i,k) = zlmin + gr%zt(i,j) - gr%zt(i,k)

          if ( gr%grid_dir_indx * (j+gr%grid_dir_indx) &
               < gr%grid_dir_indx * gr%k_ub_zt ) then
            dCAPE_dz_j_minus_1 = parcel_buoyancy(i,k,j)
            tke = parcel_tke(i,k,j)
            j = j + gr%grid_dir_indx
            dCAPE_dz_j = parcel_buoyancy(i,k,j)

            if ( abs( dCAPE_dz_j - dCAPE_dz_j_minus_1 ) * two <= &
                 abs( dCAPE_dz_j + dCAPE_dz_j_minus_1 ) * eps ) then
              pdf9_lscale_up(i,k) = pdf9_lscale_up(i,k) - tke / dCAPE_dz_j
            else
              if ( gr%grid_dir_indx > 0 ) then
                j_zm = j
              else
                j_zm = j+1
              end if

              invrs_dCAPE_diff = one / ( dCAPE_dz_j - dCAPE_dz_j_minus_1 )
              pdf9_lscale_up(i,k) = pdf9_lscale_up(i,k) &
                - dCAPE_dz_j_minus_1 * invrs_dCAPE_diff &
                  * gr%grid_dir * gr%dzm(i,j_zm) &
                - sqrt( dCAPE_dz_j_minus_1**2 &
                        - two * tke * gr%grid_dir * gr%invrs_dzm(i,j_zm) &
                          * ( dCAPE_dz_j - dCAPE_dz_j_minus_1 ) ) &
                  * invrs_dCAPE_diff * gr%grid_dir * gr%dzm(i,j_zm)
            end if
          end if
        else
          if ( gr%grid_dir_indx > 0 ) then
            kp1_zm = k+1
          else
            kp1_zm = k
          end if

          pdf9_lscale_up(i,k) = zlmin &
            - sqrt( -two * tke_i(i,k) * gr%grid_dir * gr%dzm(i,kp1_zm) &
                    * dCAPE_dz_1(i,k+gr%grid_dir_indx) ) &
              / dCAPE_dz_1(i,k+gr%grid_dir_indx)
        end if

        if ( gr%zt(i,k) + pdf9_lscale_up(i,k) < Lscale_up_max_alt ) then
          pdf9_lscale_up(i,k) = Lscale_up_max_alt - gr%zt(i,k)
        else
          Lscale_up_max_alt = pdf9_lscale_up(i,k) + gr%zt(i,k)
        end if
      end do
    end do

    ! Diagnose a provisional transport ensemble at each destination.  Positive
    ! crossing speed is used as a mass-flux proxy weight; these moments are
    ! diagnostics only and do not yet alter PDF component 3.
    do j = 1, nzt
      do i = 1, ngrdcol
        do k = 1, nzt
          if ( parcel_status(i,k,j) > zero ) then
            crossing_wgt = parcel_w(i,k,j)
            crossing_weight(i,j) = crossing_weight(i,j) + crossing_wgt
            crossing_mean_w(i,j) = crossing_mean_w(i,j) &
              + crossing_wgt * parcel_w(i,k,j)
            crossing_mean_rt(i,j) = crossing_mean_rt(i,j) &
              + crossing_wgt * parcel_rt(i,k,j)
            crossing_mean_thl(i,j) = crossing_mean_thl(i,j) &
              + crossing_wgt * parcel_thl(i,k,j)
          end if
        end do

        if ( crossing_weight(i,j) > zero ) then
          crossing_mean_w(i,j) = crossing_mean_w(i,j) / crossing_weight(i,j)
          crossing_mean_rt(i,j) = crossing_mean_rt(i,j) / crossing_weight(i,j)
          crossing_mean_thl(i,j) = crossing_mean_thl(i,j) / crossing_weight(i,j)

          do k = 1, nzt
            if ( parcel_status(i,k,j) > zero ) then
              crossing_wgt = parcel_w(i,k,j)
              delta_w = parcel_w(i,k,j) - crossing_mean_w(i,j)
              delta_rt = parcel_rt(i,k,j) - crossing_mean_rt(i,j)
              delta_thl = parcel_thl(i,k,j) - crossing_mean_thl(i,j)
              crossing_var_w(i,j) = crossing_var_w(i,j) + crossing_wgt * delta_w**2
              crossing_var_rt(i,j) = crossing_var_rt(i,j) + crossing_wgt * delta_rt**2
              crossing_var_thl(i,j) = crossing_var_thl(i,j) + crossing_wgt * delta_thl**2
              crossing_covar_w_rt(i,j) = crossing_covar_w_rt(i,j) &
                + crossing_wgt * delta_w * delta_rt
              crossing_covar_w_thl(i,j) = crossing_covar_w_thl(i,j) &
                + crossing_wgt * delta_w * delta_thl
              crossing_covar_rt_thl(i,j) = crossing_covar_rt_thl(i,j) &
                + crossing_wgt * delta_rt * delta_thl
            end if
          end do

          crossing_var_w(i,j) = crossing_var_w(i,j) / crossing_weight(i,j)
          crossing_var_rt(i,j) = crossing_var_rt(i,j) / crossing_weight(i,j)
          crossing_var_thl(i,j) = crossing_var_thl(i,j) / crossing_weight(i,j)
          crossing_covar_w_rt(i,j) = crossing_covar_w_rt(i,j) / crossing_weight(i,j)
          crossing_covar_w_thl(i,j) = crossing_covar_w_thl(i,j) / crossing_weight(i,j)
          crossing_covar_rt_thl(i,j) = crossing_covar_rt_thl(i,j) / crossing_weight(i,j)
        end if
      end do
    end do

    do k = 1, nzt
      do i = 1, ngrdcol
        if ( l_implemented ) then
          lminh = max( zero_threshold, &
                       Lscale_sfclyr_depth &
                         - ( gr%zt(i,k) - gr%zm(i,gr%k_lb_zm) ) ) &
                    * lmin * invrs_Lscale_sfclyr_depth
        else
          lminh = max( zero_threshold, Lscale_sfclyr_depth - gr%zt(i,k) ) &
                    * lmin * invrs_Lscale_sfclyr_depth
        end if
        pdf9_lscale_up(i,k) = max( lminh, pdf9_lscale_up(i,k) )
      end do
    end do


    return

  end subroutine compute_upward_mix

  !=============================================================================

  subroutine compute_downward_mix(                                    &
      nzm, nzt, ngrdcol, gr, l_implemented,                           & ! In
      p_in_Pa, exner, entrain_rtm, entrain_thlm,                      & ! In
      thvm, thv_ds,                                                   & ! In
      parcel_rt_init, parcel_thl_init, parcel_tke_init,                & ! In
      mu, lmin, saturation_formula,                                   & ! In
      err_info,                                                       & ! InOut
      pdf9_lscale_down,                                               & ! Out
      parcel_tke, parcel_w, parcel_rt, parcel_thl,                    & ! Out
      parcel_buoyancy, parcel_status,                                 & ! Out
      crossing_weight, crossing_mean_w, crossing_mean_rt,            & ! Out
      crossing_mean_thl, crossing_var_w, crossing_var_rt,            & ! Out
      crossing_var_thl, crossing_covar_w_rt, crossing_covar_w_thl,   & ! Out
      crossing_covar_rt_thl )                                          ! Out

    ! Description:
    !   Duplicates the downward half of CLUBB's historical moist, nonlocal
    !   parcel-tracking mixing-length calculation.  The only intentional
    !   difference is that PDF 9 supplies the parcel's launch rt, thl, and
    !   kinetic energy explicitly.  This routine is deliberately independent
    !   so that its parcel trajectory can be extended during PDF-9 prototyping.

    use clubb_precision, only: &
      core_rknd

    use constants_clubb, only: &
      Cp, &
      Rd, &
      ep, &
      ep1, &
      ep2, &
      Lv, &
      grav, &
      fstderr, &
      zero_threshold, &
      eps, &
      one_half, &
      one, &
      two, &
      zero

    use error_code, only: &
      clubb_fatal_error

    use err_info_type_module, only: &
      err_info_type

    use grid_class, only: &
      grid

    use saturation, only: &
      sat_mixrat_liq_api

    implicit none

    real( kind = core_rknd ), parameter :: &
      zlmin = 0.1_core_rknd, & ! Absolute minimum downward reach [m]
      Lscale_sfclyr_depth = 500._core_rknd ! Surface-layer floor depth [m]

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      nzm, & ! Number of momentum levels [-]
      nzt, & ! Number of thermodynamic levels [-]
      ngrdcol, & ! Number of grid columns [-]
      saturation_formula ! Saturation formula selector [-]

    type(grid), intent(in) :: &
      gr ! CLUBB vertical grid

    logical, intent(in) :: &
      l_implemented ! True when CLUBB is embedded in a host model [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      p_in_Pa, & ! Pressure [Pa]
      exner, & ! Exner function [-]
      entrain_rtm, & ! Environmental total-water entrainment target [kg/kg]
      entrain_thlm, & ! Environmental theta-l entrainment target [K]
      thvm, & ! Grid-mean virtual potential temperature [K]
      thv_ds, & ! Dry base-state virtual potential temperature [K]
      parcel_rt_init, & ! Parcel total water at launch [kg/kg]
      parcel_thl_init, & ! Parcel liquid-water potential temperature at launch [K]
      parcel_tke_init ! Parcel launch kinetic energy [m^2/s^2]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      mu ! Fractional entrainment rate [1/m]

    real( kind = core_rknd ), intent(in) :: &
      lmin ! Surface-layer minimum mixing length [m]

    !----------------------- Input/Output Variables ------------------------
    type(err_info_type), intent(inout) :: &
      err_info ! CLUBB error state

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      pdf9_lscale_down, & ! PDF-9 downward parcel reach [m]
      crossing_weight, & ! Sum of downward-speed crossing weights [m/s]
      crossing_mean_w, & ! Crossing-weighted signed vertical velocity [m/s]
      crossing_mean_rt, & ! Crossing-weighted total water [kg/kg]
      crossing_mean_thl, & ! Crossing-weighted liquid-water potential temperature [K]
      crossing_var_w, & ! Crossing-weighted w variance [m2/s2]
      crossing_var_rt, & ! Crossing-weighted rt variance [kg2/kg2]
      crossing_var_thl, & ! Crossing-weighted thl variance [K2]
      crossing_covar_w_rt, & ! Crossing-weighted covariance of w and rt [m/s kg/kg]
      crossing_covar_w_thl, & ! Crossing-weighted covariance of w and thl [K m/s]
      crossing_covar_rt_thl ! Crossing-weighted covariance of rt and thl [K kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,nzt), intent(out) :: &
      parcel_tke, & ! Candidate parcel energy at (launch,destination) [m2/s2]
      parcel_w, & ! Signed negative parcel speed at (launch,destination) [m/s]
      parcel_rt, & ! Parcel total water at (launch,destination) [kg/kg]
      parcel_thl, & ! Parcel liquid-water potential temperature [K]
      parcel_buoyancy, & ! CAPE gradient at destination [m/s2]
      parcel_status ! 1 reached, -1 first failed crossing, 0 unvisited [-]

    !---------------------------- Local Variables --------------------------
    integer :: &
      i, & ! Grid-column index [-]
      j, & ! Parcel destination-level index [-]
      k, & ! Parcel launch-level index [-]
      start_index, & ! First level used by saturation adjustment [-]
      jp1_zm, & ! Momentum-level index associated with j+1 [-]
      k_zm, & ! Momentum-level index associated with k [-]
      step ! Number of traversed layers from a launch level [-]

    real( kind = core_rknd ) :: &
      tke, & ! Remaining parcel kinetic energy [m^2/s^2]
      CAPE_incr, & ! CAPE increment across one layer [m^2/s^2]
      dCAPE_dz_j, & ! CAPE gradient at level j [m/s^2]
      dCAPE_dz_j_plus_1, & ! CAPE gradient at the preceding level [m/s^2]
      lminh, & ! Height-dependent minimum reach [m]
      thl_par_j, & ! Parcel thl at level j [K]
      rt_par_j, & ! Parcel rt at level j [kg/kg]
      rc_par_j, & ! Parcel cloud water at level j [kg/kg]
      thv_par_j, & ! Parcel virtual potential temperature at level j [K]
      tl_par_j, & ! Parcel liquid-water temperature at level j [K]
      rsatl_par_j, & ! Parcel saturation mixing ratio at level j [kg/kg]
      s_par_j, & ! Parcel extended liquid water at level j [kg/kg]
      Lscale_down_min_alt, & ! Lowest reachable altitude from upper launches [m]
      Lv2_coef, & ! Precomputed latent-heating coefficient [K^2]
      tl_par_j_sqd, & ! Squared parcel temperature [K^2]
      invrs_dCAPE_diff, & ! Inverse CAPE-gradient difference [s^2/m]
      invrs_Lscale_sfclyr_depth ! Inverse surface-layer depth [1/m]

    real( kind = core_rknd ) :: &
      candidate_tke, & ! Energy after a candidate downward traversal [m2/s2]
      crossing_wgt, & ! Positive magnitude of downward crossing speed [m/s]
      delta_w, & ! Difference from destination crossing mean [m/s]
      delta_rt, & ! Difference from destination crossing mean [kg/kg]
      delta_thl ! Difference from destination crossing mean [K]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      grav_on_thvm, & ! Gravity divided by environmental thv [m/s^2/K]
      Lv_coef, & ! Latent-heating coefficient [K]
      thl_par_j_precalc, & ! Recursive thl environment contribution [K]
      rt_par_j_precalc, & ! Recursive rt environment contribution [kg/kg]
      tl_par_1, & ! Parcel temperature after the first layer [K]
      rt_par_1, & ! Parcel rt after the first layer [kg/kg]
      rsatl_par_1, & ! Parcel saturation mixing ratio after the first layer [kg/kg]
      thl_par_1, & ! Parcel thl after the first layer [K]
      dCAPE_dz_1, & ! Initial CAPE gradient [m/s^2]
      s_par_1, & ! Initial parcel extended liquid water [kg/kg]
      rc_par_1, & ! Initial parcel cloud water [kg/kg]
      CAPE_incr_1, & ! Initial CAPE increment [m^2/s^2]
      thv_par_1, & ! Initial parcel virtual potential temperature [K]
      tke_i, & ! Parcel launch kinetic energy [m^2/s^2]
      tke_state, & ! Energy of the last successfully reached level [m2/s2]
      thl_state, & ! thl at the last successfully reached level [K]
      rt_state, & ! rt at the last successfully reached level [kg/kg]
      dCAPE_state ! CAPE gradient at the last successfully reached level [m/s2]

    logical, dimension(ngrdcol,nzt) :: &
      parcel_active ! True while a launched parcel can traverse another layer

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      exp_mu_dzm, & ! Layer entrainment exponential [-]
      invrs_dzm_on_mu, & ! Layer integration coefficient [-]
      entrain_coef ! Linear-environment entrainment coefficient [-]

    !--------------------------- Begin Code ---------------------------

    do i = 1, ngrdcol
      if ( abs(mu(i)) < eps ) then
        err_info%err_code(i) = clubb_fatal_error
        write(fstderr,*) err_info%err_header(i)
        write(fstderr,*) "mu = ", mu(i)
      end if
    end do

    if ( any(err_info%err_code == clubb_fatal_error) ) then
      write(fstderr,*) "Entrainment rate mu cannot be 0"
      write(fstderr,*) "Fatal error in subroutine compute_downward_mix"
      return
    end if

    do k = 1, nzt
      do i = 1, ngrdcol
        tke_i(i,k) = parcel_tke_init(i,k)
        pdf9_lscale_down(i,k) = zlmin
        grav_on_thvm(i,k) = grav / thvm(i,k)
        Lv_coef(i,k) = Lv / ( exner(i,k) * Cp ) - ep2 * thv_ds(i,k)
      end do
    end do

    do i = 1, ngrdcol
      do k = 1, nzm
        exp_mu_dzm(i,k) = exp( -mu(i) * gr%grid_dir * gr%dzm(i,k) )
        invrs_dzm_on_mu(i,k) = ( gr%grid_dir * gr%invrs_dzm(i,k) ) / mu(i)
        entrain_coef(i,k) = ( one - exp_mu_dzm(i,k) ) * invrs_dzm_on_mu(i,k)
      end do
    end do

    Lv2_coef = ep * Lv**2 / ( Rd * Cp )
    invrs_Lscale_sfclyr_depth = one / Lscale_sfclyr_depth

    ! This block is copied from the downward half of compute_mixing_length.
    do j = gr%k_lb_zt, gr%k_ub_zt-gr%grid_dir_indx, gr%grid_dir_indx
      do i = 1, ngrdcol
        if ( gr%grid_dir_indx > 0 ) then
          jp1_zm = j+1
        else
          jp1_zm = j
        end if

        thl_par_j_precalc(i,j) = entrain_thlm(i,j) &
          - entrain_thlm(i,j+gr%grid_dir_indx) * exp_mu_dzm(i,jp1_zm) &
          - ( entrain_thlm(i,j) - entrain_thlm(i,j+gr%grid_dir_indx) ) &
            * entrain_coef(i,jp1_zm)
        rt_par_j_precalc(i,j) = entrain_rtm(i,j) &
          - entrain_rtm(i,j+gr%grid_dir_indx) * exp_mu_dzm(i,jp1_zm) &
          - ( entrain_rtm(i,j) - entrain_rtm(i,j+gr%grid_dir_indx) ) &
            * entrain_coef(i,jp1_zm)
      end do
    end do

    thl_par_1 = zero
    tl_par_1 = zero
    rt_par_1 = zero
    do j = gr%k_lb_zt, gr%k_ub_zt-gr%grid_dir_indx, gr%grid_dir_indx
      do i = 1, ngrdcol
        if ( gr%grid_dir_indx > 0 ) then
          jp1_zm = j+1
        else
          jp1_zm = j
        end if

        thl_par_1(i,j) = entrain_thlm(i,j) &
          - ( entrain_thlm(i,j) - parcel_thl_init(i,j+gr%grid_dir_indx) ) &
            * entrain_coef(i,jp1_zm)
        tl_par_1(i,j) = thl_par_1(i,j) * exner(i,j)
        rt_par_1(i,j) = entrain_rtm(i,j) &
          - ( entrain_rtm(i,j) - parcel_rt_init(i,j+gr%grid_dir_indx) ) &
            * entrain_coef(i,jp1_zm)
      end do
    end do

    start_index = gr%k_lb_zt
    rsatl_par_1 = sat_mixrat_liq_api( nzt, ngrdcol, gr, p_in_Pa, tl_par_1, &
                                      saturation_formula, start_index )

    ! Use the same launch-by-destination ledger convention as the upward
    ! calculation.  Downward velocity is stored with its physical negative
    ! sign; its magnitude is used later as the positive crossing weight.
    parcel_tke = zero
    parcel_w = zero
    parcel_rt = zero
    parcel_thl = zero
    parcel_buoyancy = zero
    parcel_status = zero
    crossing_weight = zero
    crossing_mean_w = zero
    crossing_mean_rt = zero
    crossing_mean_thl = zero
    crossing_var_w = zero
    crossing_var_rt = zero
    crossing_var_thl = zero
    crossing_covar_w_rt = zero
    crossing_covar_w_thl = zero
    crossing_covar_rt_thl = zero

    do k = 1, nzt
      do i = 1, ngrdcol
        parcel_tke(i,k,k) = tke_i(i,k)
        parcel_w(i,k,k) = -sqrt( max( two * tke_i(i,k), zero ) )
        parcel_rt(i,k,k) = parcel_rt_init(i,k)
        parcel_thl(i,k,k) = parcel_thl_init(i,k)
        parcel_status(i,k,k) = one
      end do
    end do

    do j = gr%k_lb_zt, gr%k_ub_zt-gr%grid_dir_indx, gr%grid_dir_indx
      do i = 1, ngrdcol
        if ( gr%grid_dir_indx > 0 ) then
          jp1_zm = j+1
        else
          jp1_zm = j
        end if

        tl_par_j_sqd = tl_par_1(i,j)**2
        s_par_1(i,j) = ( rt_par_1(i,j) - rsatl_par_1(i,j) ) * tl_par_j_sqd &
          / ( tl_par_j_sqd + Lv2_coef * rsatl_par_1(i,j) )
        rc_par_1(i,j) = max( s_par_1(i,j), zero_threshold )
        thv_par_1(i,j) = thl_par_1(i,j) + ep1 * thv_ds(i,j) * rt_par_1(i,j) &
          + Lv_coef(i,j) * rc_par_1(i,j)
        dCAPE_dz_1(i,j) = grav_on_thvm(i,j) * ( thv_par_1(i,j) - thvm(i,j) )
        CAPE_incr_1(i,j) = one_half * dCAPE_dz_1(i,j) &
          * gr%grid_dir * gr%dzm(i,jp1_zm)
      end do
    end do

    do k = 1, nzt
      do i = 1, ngrdcol
        parcel_active(i,k) = .false.
        tke_state(i,k) = tke_i(i,k)
        rt_state(i,k) = parcel_rt_init(i,k)
        thl_state(i,k) = parcel_thl_init(i,k)
        dCAPE_state(i,k) = zero
      end do
    end do

    ! Advance every downward launch by the same traversal count before moving
    ! to the next layer.  The inner operation is a flat (column, launch-level)
    ! loop and each parcel writes only to its own destination in the 3-D ledger.
    do step = 1, nzt
      do k = gr%k_ub_zt, gr%k_lb_zt+gr%grid_dir_indx, -gr%grid_dir_indx
        do i = 1, ngrdcol
          j = k - step * gr%grid_dir_indx

          if ( step == 1 ) then
            candidate_tke = tke_i(i,k) - CAPE_incr_1(i,j)
            parcel_tke(i,k,j) = candidate_tke
            parcel_w(i,k,j) = -sqrt( max( two * candidate_tke, zero ) )
            parcel_rt(i,k,j) = rt_par_1(i,j)
            parcel_thl(i,k,j) = thl_par_1(i,j)
            parcel_buoyancy(i,k,j) = dCAPE_dz_1(i,j)

            if ( candidate_tke > zero ) then
              parcel_status(i,k,j) = one
              parcel_active(i,k) = .true.
              tke_state(i,k) = candidate_tke
              rt_state(i,k) = rt_par_1(i,j)
              thl_state(i,k) = thl_par_1(i,j)
              dCAPE_state(i,k) = dCAPE_dz_1(i,j)
            else
              parcel_status(i,k,j) = -one
              parcel_active(i,k) = .false.
            end if

          else if ( parcel_active(i,k) ) then
            if ( gr%grid_dir_indx * j < gr%grid_dir_indx * gr%k_lb_zt ) then
              parcel_active(i,k) = .false.
              cycle
            end if

            if ( gr%grid_dir_indx > 0 ) then
              jp1_zm = j+1
            else
              jp1_zm = j
            end if

            thl_par_j = thl_par_j_precalc(i,j) &
              + thl_state(i,k) * exp_mu_dzm(i,jp1_zm)
            rt_par_j = rt_par_j_precalc(i,j) &
              + rt_state(i,k) * exp_mu_dzm(i,jp1_zm)
            tl_par_j = thl_par_j * exner(i,j)
            rsatl_par_j = sat_mixrat_liq_api( p_in_Pa(i,j), tl_par_j, &
                                              saturation_formula )
            tl_par_j_sqd = tl_par_j**2
            s_par_j = ( rt_par_j - rsatl_par_j ) * tl_par_j_sqd &
              / ( tl_par_j_sqd + Lv2_coef * rsatl_par_j )
            rc_par_j = max( s_par_j, zero_threshold )
            thv_par_j = thl_par_j + ep1 * thv_ds(i,j) * rt_par_j &
              + Lv_coef(i,j) * rc_par_j
            dCAPE_dz_j = grav_on_thvm(i,j) * ( thv_par_j - thvm(i,j) )
            CAPE_incr = one_half * ( dCAPE_dz_j + dCAPE_state(i,k) ) &
              * gr%grid_dir * gr%dzm(i,jp1_zm)
            candidate_tke = tke_state(i,k) - CAPE_incr

            parcel_tke(i,k,j) = candidate_tke
            parcel_w(i,k,j) = -sqrt( max( two * candidate_tke, zero ) )
            parcel_rt(i,k,j) = rt_par_j
            parcel_thl(i,k,j) = thl_par_j
            parcel_buoyancy(i,k,j) = dCAPE_dz_j

            if ( candidate_tke > zero ) then
              parcel_status(i,k,j) = one
              tke_state(i,k) = candidate_tke
              rt_state(i,k) = rt_par_j
              thl_state(i,k) = thl_par_j
              dCAPE_state(i,k) = dCAPE_dz_j
            else
              parcel_status(i,k,j) = -one
              parcel_active(i,k) = .false.
            end if
          end if
        end do
      end do
    end do

    ! Recover the historical stopping distance from the completed ledger,
    ! retaining the first failed layer for its sub-grid CAPE interpolation.
    do i = 1, ngrdcol
      Lscale_down_min_alt = gr%zt(i,gr%k_ub_zt)
      do k = gr%k_ub_zt, gr%k_lb_zt+gr%grid_dir_indx, -gr%grid_dir_indx
        j = k - gr%grid_dir_indx

        if ( parcel_status(i,k,j) > zero ) then
          do while ( gr%grid_dir_indx * (j-gr%grid_dir_indx) &
                     >= gr%grid_dir_indx * gr%k_lb_zt )
            if ( parcel_status(i,k,j-gr%grid_dir_indx) <= zero ) exit
            j = j - gr%grid_dir_indx
          end do

          pdf9_lscale_down(i,k) = zlmin + gr%zt(i,k) - gr%zt(i,j)

          if ( gr%grid_dir_indx * (j-gr%grid_dir_indx) &
               >= gr%grid_dir_indx * gr%k_lb_zt ) then
            dCAPE_dz_j_plus_1 = parcel_buoyancy(i,k,j)
            tke = parcel_tke(i,k,j)
            j = j - gr%grid_dir_indx
            dCAPE_dz_j = parcel_buoyancy(i,k,j)

            if ( abs( dCAPE_dz_j - dCAPE_dz_j_plus_1 ) * two <= &
                 abs( dCAPE_dz_j + dCAPE_dz_j_plus_1 ) * eps ) then
              pdf9_lscale_down(i,k) = pdf9_lscale_down(i,k) &
                + tke / dCAPE_dz_j
            else
              if ( gr%grid_dir_indx > 0 ) then
                jp1_zm = j+1
              else
                jp1_zm = j
              end if
              invrs_dCAPE_diff = one / ( dCAPE_dz_j - dCAPE_dz_j_plus_1 )
              pdf9_lscale_down(i,k) = pdf9_lscale_down(i,k) &
                - dCAPE_dz_j_plus_1 * invrs_dCAPE_diff &
                  * gr%grid_dir * gr%dzm(i,jp1_zm) &
                + sqrt( dCAPE_dz_j_plus_1**2 &
                        + two * tke * gr%grid_dir * gr%invrs_dzm(i,jp1_zm) &
                          * ( dCAPE_dz_j - dCAPE_dz_j_plus_1 ) ) &
                  * invrs_dCAPE_diff * gr%grid_dir * gr%dzm(i,jp1_zm)
            end if
          end if
        else
          if ( gr%grid_dir_indx > 0 ) then
            k_zm = k
          else
            k_zm = k+1
          end if
          pdf9_lscale_down(i,k) = zlmin &
            + sqrt( two * tke_i(i,k) * gr%grid_dir * gr%dzm(i,k_zm) &
                    * dCAPE_dz_1(i,k-gr%grid_dir_indx) ) &
              / dCAPE_dz_1(i,k-gr%grid_dir_indx)
        end if

        if ( gr%zt(i,k) - pdf9_lscale_down(i,k) > Lscale_down_min_alt ) then
          pdf9_lscale_down(i,k) = gr%zt(i,k) - Lscale_down_min_alt
        else
          Lscale_down_min_alt = gr%zt(i,k) - pdf9_lscale_down(i,k)
        end if
      end do
    end do

    ! Aggregate all successfully sinking parcels at each destination.  The
    ! mean w remains negative, while -w provides a positive speed weight.
    do j = 1, nzt
      do i = 1, ngrdcol
        do k = 1, nzt
          if ( parcel_status(i,k,j) > zero ) then
            crossing_wgt = -parcel_w(i,k,j)
            crossing_weight(i,j) = crossing_weight(i,j) + crossing_wgt
            crossing_mean_w(i,j) = crossing_mean_w(i,j) &
              + crossing_wgt * parcel_w(i,k,j)
            crossing_mean_rt(i,j) = crossing_mean_rt(i,j) &
              + crossing_wgt * parcel_rt(i,k,j)
            crossing_mean_thl(i,j) = crossing_mean_thl(i,j) &
              + crossing_wgt * parcel_thl(i,k,j)
          end if
        end do

        if ( crossing_weight(i,j) > zero ) then
          crossing_mean_w(i,j) = crossing_mean_w(i,j) / crossing_weight(i,j)
          crossing_mean_rt(i,j) = crossing_mean_rt(i,j) / crossing_weight(i,j)
          crossing_mean_thl(i,j) = crossing_mean_thl(i,j) / crossing_weight(i,j)

          do k = 1, nzt
            if ( parcel_status(i,k,j) > zero ) then
              crossing_wgt = -parcel_w(i,k,j)
              delta_w = parcel_w(i,k,j) - crossing_mean_w(i,j)
              delta_rt = parcel_rt(i,k,j) - crossing_mean_rt(i,j)
              delta_thl = parcel_thl(i,k,j) - crossing_mean_thl(i,j)
              crossing_var_w(i,j) = crossing_var_w(i,j) + crossing_wgt * delta_w**2
              crossing_var_rt(i,j) = crossing_var_rt(i,j) + crossing_wgt * delta_rt**2
              crossing_var_thl(i,j) = crossing_var_thl(i,j) + crossing_wgt * delta_thl**2
              crossing_covar_w_rt(i,j) = crossing_covar_w_rt(i,j) &
                + crossing_wgt * delta_w * delta_rt
              crossing_covar_w_thl(i,j) = crossing_covar_w_thl(i,j) &
                + crossing_wgt * delta_w * delta_thl
              crossing_covar_rt_thl(i,j) = crossing_covar_rt_thl(i,j) &
                + crossing_wgt * delta_rt * delta_thl
            end if
          end do

          crossing_var_w(i,j) = crossing_var_w(i,j) / crossing_weight(i,j)
          crossing_var_rt(i,j) = crossing_var_rt(i,j) / crossing_weight(i,j)
          crossing_var_thl(i,j) = crossing_var_thl(i,j) / crossing_weight(i,j)
          crossing_covar_w_rt(i,j) = crossing_covar_w_rt(i,j) / crossing_weight(i,j)
          crossing_covar_w_thl(i,j) = crossing_covar_w_thl(i,j) / crossing_weight(i,j)
          crossing_covar_rt_thl(i,j) = crossing_covar_rt_thl(i,j) / crossing_weight(i,j)
        end if
      end do
    end do

    do k = 1, nzt
      do i = 1, ngrdcol
        if ( l_implemented ) then
          lminh = max( zero_threshold, &
                       Lscale_sfclyr_depth &
                         - ( gr%zt(i,k) - gr%zm(i,gr%k_lb_zm) ) ) &
                    * lmin * invrs_Lscale_sfclyr_depth
        else
          lminh = max( zero_threshold, Lscale_sfclyr_depth - gr%zt(i,k) ) &
                    * lmin * invrs_Lscale_sfclyr_depth
        end if
        pdf9_lscale_down(i,k) = max( lminh, pdf9_lscale_down(i,k) )
      end do
    end do


    return

  end subroutine compute_downward_mix

  !-----------------------------------------------------------------------------
  subroutine diagnose_pdf9_mixing_reach(                              &
      nzm, nzt, ngrdcol, gr, l_implemented,                            & ! In
      p_in_Pa, exner, rtm, thlm, thvm, wp2, wprtp, wpthlp,             & ! In
      thv_ds, lmin, clubb_params, saturation_formula,                  & ! In
      stats, err_info,                                                 & ! InOut
      w_3, rt_3, thl_3,                                                & ! InOut
      varnce_w_3, varnce_rt_3, varnce_thl_3,                          & ! InOut
      corr_w_rt_3, corr_w_thl_3, corr_rt_thl_3,                       & ! InOut
      mixt_frac_3,                                                     & ! In
      pdf9_candidate_valid,                                          & ! InOut
      pdf9_lscale_up, pdf9_lscale_down )                               ! Out

    ! Description:
    !   Diagnoses how far locally representative rising and sinking parcels
    !   can travel.  This is intentionally a diagnostic-only stage: it neither
    !   moves PDF component 3 nor changes CLUBB's operational Lscale.
    !
    ! The launch perturbations are conditional means from a linear moment
    ! model.  For the rising parcel,
    !
    !   w'_p     = +sqrt(2/pi) sqrt(w'^2)
    !   r_t'_p   = (w'r_t'/w'^2) w'_p
    !   theta'_p = (w'theta_l'/w'^2) w'_p,
    !
    ! and the signs reverse for the sinking parcel.  Once launched, the two
    ! independent PDF-9 solvers use duplicated copies of the historical
    ! entrainment, saturation, buoyancy, CAPE integration, and stopping rules.

    use clubb_precision, only: &
      core_rknd

    use constants_clubb, only: &
      one, &
      one_half, &
      pi, &
      two, &
      w_tol_sqd

    use error_code, only: &
      clubb_fatal_error

    use err_info_type_module, only: &
      err_info_type

    use grid_class, only: &
      grid, &
      zm2zt_api

    use parameter_indices, only: &
      imu, &
      nparams

    use stats_netcdf, only: &
      stats_type, &
      stats_update

    implicit none

    integer, intent(in) :: &
      nzm, &
      nzt, &
      ngrdcol, &
      saturation_formula

    type(grid), intent(in) :: &
      gr

    logical, intent(in) :: &
      l_implemented

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      p_in_Pa, &
      exner, &
      rtm, &
      thlm, &
      thvm, &
      thv_ds

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      wp2, &
      wprtp, &
      wpthlp

    real( kind = core_rknd ), intent(in) :: &
      lmin

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params

    type(stats_type), intent(inout) :: &
      stats

    type(err_info_type), intent(inout) :: &
      err_info

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      pdf9_lscale_up, &
      pdf9_lscale_down

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
      w_3, &
      rt_3, &
      thl_3, &
      varnce_w_3, &
      varnce_rt_3, &
      varnce_thl_3, &
      corr_w_rt_3, &
      corr_w_thl_3, &
      corr_rt_thl_3, &
      pdf9_candidate_valid

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      mixt_frac_3

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      wp2_zt, &
      wprtp_zt, &
      wpthlp_zt, &
      launch_w, &
      launch_tke, &
      launch_tke_down, &
      launch_rt_up, &
      launch_rt_down, &
      launch_thl_up, &
      launch_thl_down, &
      entrain_rt_up, &
      entrain_rt_down, &
      entrain_thl_up, &
      entrain_thl_down, &
      entrain_g3_weight, &
      crossing_weight, &
      crossing_mean_w, &
      crossing_mean_rt, &
      crossing_mean_thl, &
      crossing_var_w, &
      crossing_var_rt, &
      crossing_var_thl, &
      crossing_covar_w_rt, &
      crossing_covar_w_thl, &
      crossing_covar_rt_thl, &
      candidate_crossing_count_up, &
      candidate_crossing_count_down, &
      candidate_crossing_count_combined, &
      candidate_weighted_support_up, &
      candidate_weighted_support_down, &
      candidate_weighted_support_combined, &
      candidate_donor_distance_up, &
      candidate_donor_distance_down, &
      candidate_donor_distance_combined, &
      candidate_launch_from_g3, &
      candidate_down_launch_from_g3, &
      candidate_branch_prob_up, &
      candidate_branch_prob_down, &
      candidate_valid_up, &
      candidate_valid_down, &
      candidate_w_up, &
      candidate_rt_up, &
      candidate_thl_up, &
      candidate_var_w_up, &
      candidate_var_rt_up, &
      candidate_var_thl_up, &
      candidate_w_down, &
      candidate_w_down_uncapped, &
      candidate_down_cap_fraction, &
      candidate_destination_sigma_w, &
      candidate_rt_down, &
      candidate_thl_down, &
      candidate_var_w_down, &
      candidate_var_rt_down, &
      candidate_var_thl_down, &
      candidate_covar_w_rt_up, &
      candidate_covar_w_thl_up, &
      candidate_covar_rt_thl_up, &
      candidate_covar_w_rt_down, &
      candidate_covar_w_thl_down, &
      candidate_covar_rt_thl_down, &
      candidate_corr_w_rt_up, &
      candidate_corr_w_thl_up, &
      candidate_corr_rt_thl_up, &
      candidate_corr_w_rt_down, &
      candidate_corr_w_thl_down, &
      candidate_corr_rt_thl_down

    real( kind = core_rknd ), dimension(ngrdcol,nzt,nzt) :: &
      parcel_tke, &
      parcel_w, &
      parcel_rt, &
      parcel_thl, &
      parcel_buoyancy, &
      parcel_status

    integer :: &
      i, &
      j, &
      k

    real( kind = core_rknd ) :: &
      candidate_count, &
      count_up, &
      count_down, &
      support_up, &
      support_down, &
      candidate_support, &
      capped_support, &
      donor_distance, &
      provenance_weight, &
      effective_arrival_w, &
      sigma_w, &
      standardized_w, &
      normal_density, &
      probability_up, &
      probability_down, &
      conditional_w, &
      launch_speed_cap, &
      covariance_w_rt, &
      covariance_w_thl, &
      active_weight, &
      g3_rt_up, &
      g3_rt_down, &
      g3_thl_up, &
      g3_thl_down, &
      delta_w, &
      delta_rt, &
      delta_thl, &
      denom


    wp2_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, wp2(:,:), w_tol_sqd )
    wprtp_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, wprtp(:,:) )
    wpthlp_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, wpthlp(:,:) )

    do k = 1, nzt
      do i = 1, ngrdcol
        launch_w(i,k) = sqrt( two / pi ) * sqrt( max( wp2_zt(i,k), w_tol_sqd ) )
        launch_tke(i,k) = one_half * launch_w(i,k)**2
        launch_tke_down(i,k) = launch_tke(i,k)

        launch_rt_up(i,k) = rtm(i,k) &
          + wprtp_zt(i,k) / max( wp2_zt(i,k), w_tol_sqd ) * launch_w(i,k)
        launch_rt_down(i,k) = rtm(i,k) &
          - wprtp_zt(i,k) / max( wp2_zt(i,k), w_tol_sqd ) * launch_w(i,k)

        launch_thl_up(i,k) = thlm(i,k) &
          + wpthlp_zt(i,k) / max( wp2_zt(i,k), w_tol_sqd ) * launch_w(i,k)
        launch_thl_down(i,k) = thlm(i,k) &
          - wpthlp_zt(i,k) / max( wp2_zt(i,k), w_tol_sqd ) * launch_w(i,k)

        active_weight = max( 0._core_rknd, min( mixt_frac_3(i,k), &
                                                1._core_rknd ) )

        g3_rt_up = rt_3(i,k)
        g3_rt_down = rt_3(i,k)
        g3_thl_up = thl_3(i,k)
        g3_thl_down = thl_3(i,k)

        candidate_launch_from_g3(i,k) = 0._core_rknd
        candidate_down_launch_from_g3(i,k) = 0._core_rknd
        ! The ordinary moment-conditioned launch represents one half of a
        ! symmetric launch population.  A retained G3 below replaces these
        ! defaults with the exact positive- and negative-w probability masses.
        candidate_branch_prob_up(i,k) = one_half
        candidate_branch_prob_down(i,k) = one_half
        if ( k /= gr%k_lb_zt .and. pdf9_candidate_valid(i,k) > one_half &
             .and. varnce_w_3(i,k) > w_tol_sqd ) then
          ! Component 3 now represents the pooled population of rising and
          ! sinking transported air.  Launch each direction from the matching
          ! half of that Gaussian rather than from its (possibly near-zero)
          ! unconditional mean.  Linear Gaussian conditioning carries the
          ! corresponding rt and theta-l shifts along with w.
          sigma_w = sqrt( varnce_w_3(i,k) )
          standardized_w = w_3(i,k) / sigma_w
          normal_density = exp( -one_half * standardized_w**2 ) &
            / sqrt( two * pi )
          probability_up = one_half * ( 1._core_rknd &
            + erf( standardized_w / sqrt(two) ) )
          probability_down = 1._core_rknd - probability_up
          candidate_branch_prob_up(i,k) = probability_up
          candidate_branch_prob_down(i,k) = probability_down

          covariance_w_rt = corr_w_rt_3(i,k) &
            * sqrt( max( varnce_w_3(i,k) * varnce_rt_3(i,k), 0._core_rknd ) )
          covariance_w_thl = corr_w_thl_3(i,k) &
            * sqrt( max( varnce_w_3(i,k) * varnce_thl_3(i,k), 0._core_rknd ) )

          ! A recursively diagnosed component can acquire an extreme tail
          ! mean when one directional probability is small.  Do not let that
          ! algebraic tail speed feed back into a still more energetic parcel
          ! on the next pass.  One local grid-mean RMS is a parameter-free
          ! physical scale; ordinary moment-conditioned launches above remain
          ! unchanged at sqrt(2/pi) times that RMS.
          launch_speed_cap = sqrt( max( wp2_zt(i,k), w_tol_sqd ) )

          if ( probability_up > 1.e-12_core_rknd ) then
            conditional_w = w_3(i,k) &
              + sigma_w * normal_density / probability_up
            conditional_w = max( 0._core_rknd, &
              min( conditional_w, launch_speed_cap ) )
            launch_w(i,k) = conditional_w
            launch_tke(i,k) = one_half * conditional_w**2
            launch_rt_up(i,k) = rt_3(i,k) + covariance_w_rt &
              / varnce_w_3(i,k) * ( conditional_w - w_3(i,k) )
            launch_thl_up(i,k) = thl_3(i,k) + covariance_w_thl &
              / varnce_w_3(i,k) * ( conditional_w - w_3(i,k) )
            g3_rt_up = launch_rt_up(i,k)
            g3_thl_up = launch_thl_up(i,k)
            candidate_launch_from_g3(i,k) = 1._core_rknd
          end if

          if ( probability_down > 1.e-12_core_rknd ) then
            conditional_w = w_3(i,k) &
              - sigma_w * normal_density / probability_down
            conditional_w = min( 0._core_rknd, &
              max( conditional_w, -launch_speed_cap ) )
            launch_tke_down(i,k) = one_half * conditional_w**2
            launch_rt_down(i,k) = rt_3(i,k) + covariance_w_rt &
              / varnce_w_3(i,k) * ( conditional_w - w_3(i,k) )
            launch_thl_down(i,k) = thl_3(i,k) + covariance_w_thl &
              / varnce_w_3(i,k) * ( conditional_w - w_3(i,k) )
            g3_rt_down = launch_rt_down(i,k)
            g3_thl_down = launch_thl_down(i,k)
            candidate_down_launch_from_g3(i,k) = 1._core_rknd
          end if
        end if

        if ( l_pdf9_weighted_g3_entrainment ) then
          entrain_g3_weight(i,k) = active_weight
          entrain_rt_up(i,k) = ( one - active_weight ) * rtm(i,k) &
            + active_weight * g3_rt_up
          entrain_rt_down(i,k) = ( one - active_weight ) * rtm(i,k) &
            + active_weight * g3_rt_down
          entrain_thl_up(i,k) = ( one - active_weight ) * thlm(i,k) &
            + active_weight * g3_thl_up
          entrain_thl_down(i,k) = ( one - active_weight ) * thlm(i,k) &
            + active_weight * g3_thl_down
        else
          entrain_g3_weight(i,k) = 0._core_rknd
          entrain_rt_up(i,k) = rtm(i,k)
          entrain_rt_down(i,k) = rtm(i,k)
          entrain_thl_up(i,k) = thlm(i,k)
          entrain_thl_down(i,k) = thlm(i,k)
        end if
      end do
    end do

    if ( stats%l_sample ) then
      call stats_update( "pdf9_entrain_g3_weight", entrain_g3_weight, stats )
      call stats_update( "pdf9_entrain_rt_up", entrain_rt_up, stats )
      call stats_update( "pdf9_entrain_rt_down", entrain_rt_down, stats )
      call stats_update( "pdf9_entrain_thl_up", entrain_thl_up, stats )
      call stats_update( "pdf9_entrain_thl_down", entrain_thl_down, stats )
    end if

    call compute_upward_mix(                                         &
        nzm, nzt, ngrdcol, gr, l_implemented,                        & ! In
        p_in_Pa, exner, entrain_rt_up, entrain_thl_up,               & ! In
        thvm, thv_ds,                                                & ! In
        launch_rt_up, launch_thl_up, launch_tke,                     & ! In
        clubb_params(:,imu), lmin, saturation_formula,               & ! In
        err_info,                                                    & ! InOut
        pdf9_lscale_up,                                               & ! Out
        parcel_tke, parcel_w, parcel_rt, parcel_thl,                  & ! Out
        parcel_buoyancy, parcel_status,                               & ! Out
        crossing_weight, crossing_mean_w, crossing_mean_rt,          & ! Out
        crossing_mean_thl, crossing_var_w, crossing_var_rt,          & ! Out
        crossing_var_thl, crossing_covar_w_rt, crossing_covar_w_thl, & ! Out
        crossing_covar_rt_thl )                                        ! Out

    if ( any(err_info%err_code == clubb_fatal_error) ) then
      return
    end if

    ! Diagnose the upward half of the transported population.  A parcel's
    ! state has already entrained toward the traversed environment.  Its vote
    ! in the destination Gaussian is additionally weighted by the fraction of
    ! launch identity retained over the donor distance,
    !
    !   q_+(k -> j) = P_+(k) exp( -mu |z_j-z_k| ).
    !
    ! Thus nearby crossings carry more statistical authority without imposing
    ! a hard range cutoff, while a rare conditional branch cannot count as a
    ! full member of the transported population.  Raw crossing counts remain
    ! separate diagnostics.
    do j = 1, nzt
      do i = 1, ngrdcol
        candidate_count = 0._core_rknd
        candidate_support = 0._core_rknd
        candidate_donor_distance_up(i,j) = 0._core_rknd
        candidate_w_up(i,j) = 0._core_rknd
        candidate_rt_up(i,j) = 0._core_rknd
        candidate_thl_up(i,j) = 0._core_rknd
        candidate_var_w_up(i,j) = 0._core_rknd
        candidate_var_rt_up(i,j) = 0._core_rknd
        candidate_var_thl_up(i,j) = 0._core_rknd
        candidate_covar_w_rt_up(i,j) = 0._core_rknd
        candidate_covar_w_thl_up(i,j) = 0._core_rknd
        candidate_covar_rt_thl_up(i,j) = 0._core_rknd

        do k = 1, nzt
          if ( k /= j .and. parcel_status(i,k,j) > 0._core_rknd ) then
            donor_distance = abs( gr%zt(i,j) - gr%zt(i,k) )
            provenance_weight = candidate_branch_prob_up(i,k) &
              * exp( -clubb_params(i,imu) * donor_distance )
            candidate_count = candidate_count + 1._core_rknd
            candidate_support = candidate_support + provenance_weight
            candidate_donor_distance_up(i,j) = candidate_donor_distance_up(i,j) &
              + provenance_weight * donor_distance
            candidate_w_up(i,j) = candidate_w_up(i,j) &
              + provenance_weight * parcel_w(i,k,j)
            candidate_rt_up(i,j) = candidate_rt_up(i,j) &
              + provenance_weight * parcel_rt(i,k,j)
            candidate_thl_up(i,j) = candidate_thl_up(i,j) &
              + provenance_weight * parcel_thl(i,k,j)
          end if
        end do

        candidate_crossing_count_up(i,j) = candidate_count
        candidate_weighted_support_up(i,j) = candidate_support
        if ( j /= gr%k_lb_zt .and. candidate_support > 0._core_rknd ) then
          candidate_donor_distance_up(i,j) = &
            candidate_donor_distance_up(i,j) / candidate_support
          candidate_w_up(i,j) = candidate_w_up(i,j) / candidate_support
          candidate_rt_up(i,j) = candidate_rt_up(i,j) / candidate_support
          candidate_thl_up(i,j) = candidate_thl_up(i,j) / candidate_support

          do k = 1, nzt
            if ( k /= j .and. parcel_status(i,k,j) > 0._core_rknd ) then
              donor_distance = abs( gr%zt(i,j) - gr%zt(i,k) )
              provenance_weight = candidate_branch_prob_up(i,k) &
                * exp( -clubb_params(i,imu) * donor_distance )
              delta_w = parcel_w(i,k,j) - candidate_w_up(i,j)
              delta_rt = parcel_rt(i,k,j) - candidate_rt_up(i,j)
              delta_thl = parcel_thl(i,k,j) - candidate_thl_up(i,j)
              candidate_var_w_up(i,j) = candidate_var_w_up(i,j) &
                + provenance_weight * delta_w**2
              candidate_var_rt_up(i,j) = candidate_var_rt_up(i,j) &
                + provenance_weight * delta_rt**2
              candidate_var_thl_up(i,j) = candidate_var_thl_up(i,j) &
                + provenance_weight * delta_thl**2
              candidate_covar_w_rt_up(i,j) = candidate_covar_w_rt_up(i,j) &
                + provenance_weight * delta_w * delta_rt
              candidate_covar_w_thl_up(i,j) = candidate_covar_w_thl_up(i,j) &
                + provenance_weight * delta_w * delta_thl
              candidate_covar_rt_thl_up(i,j) = candidate_covar_rt_thl_up(i,j) &
                + provenance_weight * delta_rt * delta_thl
            end if
          end do

          candidate_var_w_up(i,j) = candidate_var_w_up(i,j) / candidate_support
          candidate_var_rt_up(i,j) = candidate_var_rt_up(i,j) / candidate_support
          candidate_var_thl_up(i,j) = candidate_var_thl_up(i,j) / candidate_support
          candidate_covar_w_rt_up(i,j) = &
            candidate_covar_w_rt_up(i,j) / candidate_support
          candidate_covar_w_thl_up(i,j) = &
            candidate_covar_w_thl_up(i,j) / candidate_support
          candidate_covar_rt_thl_up(i,j) = &
            candidate_covar_rt_thl_up(i,j) / candidate_support

          denom = sqrt( max( candidate_var_w_up(i,j) &
                            * candidate_var_rt_up(i,j), 0._core_rknd ) )
          if ( denom > 0._core_rknd ) then
            candidate_corr_w_rt_up(i,j) = max( -1._core_rknd, &
              min( 1._core_rknd, candidate_covar_w_rt_up(i,j) / denom ) )
          else
            candidate_corr_w_rt_up(i,j) = 0._core_rknd
          end if

          denom = sqrt( max( candidate_var_w_up(i,j) &
                            * candidate_var_thl_up(i,j), 0._core_rknd ) )
          if ( denom > 0._core_rknd ) then
            candidate_corr_w_thl_up(i,j) = max( -1._core_rknd, &
              min( 1._core_rknd, candidate_covar_w_thl_up(i,j) / denom ) )
          else
            candidate_corr_w_thl_up(i,j) = 0._core_rknd
          end if

          denom = sqrt( max( candidate_var_rt_up(i,j) &
                            * candidate_var_thl_up(i,j), 0._core_rknd ) )
          if ( denom > 0._core_rknd ) then
            candidate_corr_rt_thl_up(i,j) = max( -1._core_rknd, &
              min( 1._core_rknd, candidate_covar_rt_thl_up(i,j) / denom ) )
          else
            candidate_corr_rt_thl_up(i,j) = 0._core_rknd
          end if

          candidate_valid_up(i,j) = 1._core_rknd
        else
          candidate_w_up(i,j) = 0._core_rknd
          candidate_rt_up(i,j) = rtm(i,j)
          candidate_thl_up(i,j) = thlm(i,j)
          candidate_corr_w_rt_up(i,j) = 0._core_rknd
          candidate_corr_w_thl_up(i,j) = 0._core_rknd
          candidate_corr_rt_thl_up(i,j) = 0._core_rknd
          candidate_donor_distance_up(i,j) = 0._core_rknd
          candidate_valid_up(i,j) = 0._core_rknd
        end if
      end do
    end do

    ! Save the upward ledger before the shared work arrays are reused by the
    ! downward calculation.  This avoids doubling the already-large 3-D
    ! diagnostic allocation.
    if ( stats%l_sample ) then
      call stats_update( "pdf9_lscale_up", pdf9_lscale_up, stats )
      call stats_update( "pdf9_up_parcel_tke", parcel_tke, stats )
      call stats_update( "pdf9_up_parcel_w", parcel_w, stats )
      call stats_update( "pdf9_up_parcel_rt", parcel_rt, stats )
      call stats_update( "pdf9_up_parcel_thl", parcel_thl, stats )
      call stats_update( "pdf9_up_parcel_buoyancy", parcel_buoyancy, stats )
      call stats_update( "pdf9_up_parcel_status", parcel_status, stats )
      call stats_update( "pdf9_up_crossing_weight", crossing_weight, stats )
      call stats_update( "pdf9_up_crossing_mean_w", crossing_mean_w, stats )
      call stats_update( "pdf9_up_crossing_mean_rt", crossing_mean_rt, stats )
      call stats_update( "pdf9_up_crossing_mean_thl", crossing_mean_thl, stats )
      call stats_update( "pdf9_up_crossing_var_w", crossing_var_w, stats )
      call stats_update( "pdf9_up_crossing_var_rt", crossing_var_rt, stats )
      call stats_update( "pdf9_up_crossing_var_thl", crossing_var_thl, stats )
      call stats_update( "pdf9_up_crossing_covar_w_rt", crossing_covar_w_rt, stats )
      call stats_update( "pdf9_up_crossing_covar_w_thl", crossing_covar_w_thl, stats )
      call stats_update( "pdf9_up_crossing_covar_rt_thl", crossing_covar_rt_thl, stats )
      call stats_update( "pdf9_candidate_valid_up", candidate_valid_up, stats )
      call stats_update( "pdf9_candidate_launch_from_g3_up", &
                         candidate_launch_from_g3, stats )
      call stats_update( "pdf9_candidate_launch_from_g3_down", &
                         candidate_down_launch_from_g3, stats )
      call stats_update( "pdf9_candidate_crossing_count_up", &
                         candidate_crossing_count_up, stats )
      call stats_update( "pdf9_candidate_weighted_support_up", &
                         candidate_weighted_support_up, stats )
      call stats_update( "pdf9_candidate_donor_distance_up", &
                         candidate_donor_distance_up, stats )
      call stats_update( "pdf9_candidate_w_up", candidate_w_up, stats )
      call stats_update( "pdf9_candidate_rt_up", candidate_rt_up, stats )
      call stats_update( "pdf9_candidate_thl_up", candidate_thl_up, stats )
      call stats_update( "pdf9_candidate_var_w_up", candidate_var_w_up, stats )
      call stats_update( "pdf9_candidate_var_rt_up", candidate_var_rt_up, stats )
      call stats_update( "pdf9_candidate_var_thl_up", candidate_var_thl_up, stats )
      call stats_update( "pdf9_candidate_covar_w_rt_up", &
                         candidate_covar_w_rt_up, stats )
      call stats_update( "pdf9_candidate_covar_w_thl_up", &
                         candidate_covar_w_thl_up, stats )
      call stats_update( "pdf9_candidate_covar_rt_thl_up", &
                         candidate_covar_rt_thl_up, stats )
      call stats_update( "pdf9_candidate_corr_w_rt_up", &
                         candidate_corr_w_rt_up, stats )
      call stats_update( "pdf9_candidate_corr_w_thl_up", &
                         candidate_corr_w_thl_up, stats )
      call stats_update( "pdf9_candidate_corr_rt_thl_up", &
                         candidate_corr_rt_thl_up, stats )
      call stats_update( "pdf9_candidate_branch_prob_up", &
                         candidate_branch_prob_up, stats )
      call stats_update( "pdf9_candidate_branch_prob_down", &
                         candidate_branch_prob_down, stats )
    end if

    call compute_downward_mix(                                       &
        nzm, nzt, ngrdcol, gr, l_implemented,                        & ! In
        p_in_Pa, exner, entrain_rt_down, entrain_thl_down,           & ! In
        thvm, thv_ds,                                                & ! In
        launch_rt_down, launch_thl_down, launch_tke_down,            & ! In
        clubb_params(:,imu), lmin, saturation_formula,               & ! In
        err_info,                                                    & ! InOut
        pdf9_lscale_down,                                             & ! Out
        parcel_tke, parcel_w, parcel_rt, parcel_thl,                  & ! Out
        parcel_buoyancy, parcel_status,                               & ! Out
        crossing_weight, crossing_mean_w, crossing_mean_rt,          & ! Out
        crossing_mean_thl, crossing_var_w, crossing_var_rt,          & ! Out
        crossing_var_thl, crossing_covar_w_rt, crossing_covar_w_thl, & ! Out
        crossing_covar_rt_thl )                                        ! Out

    if ( any(err_info%err_code == clubb_fatal_error) ) then
      return
    end if

    ! Diagnose the sinking half with the same branch-mass and retained-identity
    ! weighting used above.  This keeps upward and downward provenance directly
    ! comparable without promoting a rare Gaussian tail to a full parcel vote.
    do j = 1, nzt
      do i = 1, ngrdcol
        candidate_count = 0._core_rknd
        candidate_support = 0._core_rknd
        candidate_donor_distance_down(i,j) = 0._core_rknd
        candidate_w_down(i,j) = 0._core_rknd
        candidate_w_down_uncapped(i,j) = 0._core_rknd
        candidate_down_cap_fraction(i,j) = 0._core_rknd
        candidate_destination_sigma_w(i,j) = &
          sqrt( max( wp2_zt(i,j), w_tol_sqd ) )
        candidate_rt_down(i,j) = 0._core_rknd
        candidate_thl_down(i,j) = 0._core_rknd
        candidate_var_w_down(i,j) = 0._core_rknd
        candidate_var_rt_down(i,j) = 0._core_rknd
        candidate_var_thl_down(i,j) = 0._core_rknd
        candidate_covar_w_rt_down(i,j) = 0._core_rknd
        candidate_covar_w_thl_down(i,j) = 0._core_rknd
        candidate_covar_rt_thl_down(i,j) = 0._core_rknd
        capped_support = 0._core_rknd

        do k = 1, nzt
          if ( k /= j .and. parcel_status(i,k,j) > 0._core_rknd ) then
            donor_distance = abs( gr%zt(i,j) - gr%zt(i,k) )
            provenance_weight = candidate_branch_prob_down(i,k) &
              * exp( -clubb_params(i,imu) * donor_distance )
            candidate_count = candidate_count + 1._core_rknd
            candidate_support = candidate_support + provenance_weight
            candidate_donor_distance_down(i,j) = &
              candidate_donor_distance_down(i,j) &
              + provenance_weight * donor_distance
            ! Retain the raw trajectory for reach/provenance diagnostics, but
            ! limit how much sinking-speed authority one arrival has when the
            ! destination-level G3 moments are pooled.
            effective_arrival_w = max( parcel_w(i,k,j), &
              -candidate_destination_sigma_w(i,j) )
            if ( effective_arrival_w > parcel_w(i,k,j) ) then
              capped_support = capped_support + provenance_weight
            end if
            candidate_w_down(i,j) = candidate_w_down(i,j) &
              + provenance_weight * effective_arrival_w
            candidate_w_down_uncapped(i,j) = &
              candidate_w_down_uncapped(i,j) &
              + provenance_weight * parcel_w(i,k,j)
            candidate_rt_down(i,j) = candidate_rt_down(i,j) &
              + provenance_weight * parcel_rt(i,k,j)
            candidate_thl_down(i,j) = candidate_thl_down(i,j) &
              + provenance_weight * parcel_thl(i,k,j)
          end if
        end do

        candidate_crossing_count_down(i,j) = candidate_count
        candidate_weighted_support_down(i,j) = candidate_support
        if ( candidate_support > 0._core_rknd ) then
          candidate_donor_distance_down(i,j) = &
            candidate_donor_distance_down(i,j) / candidate_support
          candidate_w_down(i,j) = candidate_w_down(i,j) / candidate_support
          candidate_w_down_uncapped(i,j) = &
            candidate_w_down_uncapped(i,j) / candidate_support
          candidate_down_cap_fraction(i,j) = capped_support / candidate_support
          candidate_rt_down(i,j) = candidate_rt_down(i,j) / candidate_support
          candidate_thl_down(i,j) = candidate_thl_down(i,j) / candidate_support

          do k = 1, nzt
            if ( k /= j .and. parcel_status(i,k,j) > 0._core_rknd ) then
              donor_distance = abs( gr%zt(i,j) - gr%zt(i,k) )
              provenance_weight = candidate_branch_prob_down(i,k) &
                * exp( -clubb_params(i,imu) * donor_distance )
              effective_arrival_w = max( parcel_w(i,k,j), &
                -candidate_destination_sigma_w(i,j) )
              delta_w = effective_arrival_w - candidate_w_down(i,j)
              delta_rt = parcel_rt(i,k,j) - candidate_rt_down(i,j)
              delta_thl = parcel_thl(i,k,j) - candidate_thl_down(i,j)
              candidate_var_w_down(i,j) = candidate_var_w_down(i,j) &
                + provenance_weight * delta_w**2
              candidate_var_rt_down(i,j) = candidate_var_rt_down(i,j) &
                + provenance_weight * delta_rt**2
              candidate_var_thl_down(i,j) = candidate_var_thl_down(i,j) &
                + provenance_weight * delta_thl**2
              candidate_covar_w_rt_down(i,j) = candidate_covar_w_rt_down(i,j) &
                + provenance_weight * delta_w * delta_rt
              candidate_covar_w_thl_down(i,j) = candidate_covar_w_thl_down(i,j) &
                + provenance_weight * delta_w * delta_thl
              candidate_covar_rt_thl_down(i,j) = candidate_covar_rt_thl_down(i,j) &
                + provenance_weight * delta_rt * delta_thl
            end if
          end do

          candidate_var_w_down(i,j) = candidate_var_w_down(i,j) / candidate_support
          candidate_var_rt_down(i,j) = candidate_var_rt_down(i,j) / candidate_support
          candidate_var_thl_down(i,j) = candidate_var_thl_down(i,j) / candidate_support
          candidate_covar_w_rt_down(i,j) = &
            candidate_covar_w_rt_down(i,j) / candidate_support
          candidate_covar_w_thl_down(i,j) = &
            candidate_covar_w_thl_down(i,j) / candidate_support
          candidate_covar_rt_thl_down(i,j) = &
            candidate_covar_rt_thl_down(i,j) / candidate_support

          denom = sqrt( max( candidate_var_w_down(i,j) &
                            * candidate_var_rt_down(i,j), 0._core_rknd ) )
          if ( denom > 0._core_rknd ) then
            candidate_corr_w_rt_down(i,j) = max( -1._core_rknd, &
              min( 1._core_rknd, candidate_covar_w_rt_down(i,j) / denom ) )
          else
            candidate_corr_w_rt_down(i,j) = 0._core_rknd
          end if

          denom = sqrt( max( candidate_var_w_down(i,j) &
                            * candidate_var_thl_down(i,j), 0._core_rknd ) )
          if ( denom > 0._core_rknd ) then
            candidate_corr_w_thl_down(i,j) = max( -1._core_rknd, &
              min( 1._core_rknd, candidate_covar_w_thl_down(i,j) / denom ) )
          else
            candidate_corr_w_thl_down(i,j) = 0._core_rknd
          end if

          denom = sqrt( max( candidate_var_rt_down(i,j) &
                            * candidate_var_thl_down(i,j), 0._core_rknd ) )
          if ( denom > 0._core_rknd ) then
            candidate_corr_rt_thl_down(i,j) = max( -1._core_rknd, &
              min( 1._core_rknd, candidate_covar_rt_thl_down(i,j) / denom ) )
          else
            candidate_corr_rt_thl_down(i,j) = 0._core_rknd
          end if

          candidate_valid_down(i,j) = 1._core_rknd
        else
          candidate_w_down(i,j) = 0._core_rknd
          candidate_w_down_uncapped(i,j) = 0._core_rknd
          candidate_down_cap_fraction(i,j) = 0._core_rknd
          candidate_rt_down(i,j) = rtm(i,j)
          candidate_thl_down(i,j) = thlm(i,j)
          candidate_corr_w_rt_down(i,j) = 0._core_rknd
          candidate_corr_w_thl_down(i,j) = 0._core_rknd
          candidate_corr_rt_thl_down(i,j) = 0._core_rknd
          candidate_donor_distance_down(i,j) = 0._core_rknd
          candidate_valid_down(i,j) = 0._core_rknd
        end if
      end do
    end do

    ! Retain the next G3 candidate from the distance-weighted union of all incoming
    ! upward and downward crossings.  The parallel-axis terms below make this
    ! exactly equivalent to concatenating the two raw ledgers and applying the
    ! same provenance weights, without retaining a second 3-D parcel ledger.
    do j = 1, nzt
      do i = 1, ngrdcol
        count_up = candidate_crossing_count_up(i,j)
        count_down = candidate_crossing_count_down(i,j)
        candidate_count = count_up + count_down
        candidate_crossing_count_combined(i,j) = candidate_count
        support_up = candidate_weighted_support_up(i,j)
        support_down = candidate_weighted_support_down(i,j)
        candidate_support = support_up + support_down
        candidate_weighted_support_combined(i,j) = candidate_support

        if ( j /= gr%k_lb_zt .and. candidate_support > 0._core_rknd ) then
          candidate_donor_distance_combined(i,j) = &
            ( support_up * candidate_donor_distance_up(i,j) &
            + support_down * candidate_donor_distance_down(i,j) ) &
            / candidate_support
          w_3(i,j) = ( support_up * candidate_w_up(i,j) &
                     + support_down * candidate_w_down(i,j) ) / candidate_support
          rt_3(i,j) = ( support_up * candidate_rt_up(i,j) &
                      + support_down * candidate_rt_down(i,j) ) / candidate_support
          thl_3(i,j) = ( support_up * candidate_thl_up(i,j) &
                       + support_down * candidate_thl_down(i,j) ) / candidate_support

          delta_w = candidate_w_up(i,j) - w_3(i,j)
          varnce_w_3(i,j) = support_up &
            * ( candidate_var_w_up(i,j) + delta_w**2 )
          delta_w = candidate_w_down(i,j) - w_3(i,j)
          varnce_w_3(i,j) = ( varnce_w_3(i,j) + support_down &
            * ( candidate_var_w_down(i,j) + delta_w**2 ) ) / candidate_support

          delta_rt = candidate_rt_up(i,j) - rt_3(i,j)
          varnce_rt_3(i,j) = support_up &
            * ( candidate_var_rt_up(i,j) + delta_rt**2 )
          delta_rt = candidate_rt_down(i,j) - rt_3(i,j)
          varnce_rt_3(i,j) = ( varnce_rt_3(i,j) + support_down &
            * ( candidate_var_rt_down(i,j) + delta_rt**2 ) ) / candidate_support

          delta_thl = candidate_thl_up(i,j) - thl_3(i,j)
          varnce_thl_3(i,j) = support_up &
            * ( candidate_var_thl_up(i,j) + delta_thl**2 )
          delta_thl = candidate_thl_down(i,j) - thl_3(i,j)
          varnce_thl_3(i,j) = ( varnce_thl_3(i,j) + support_down &
            * ( candidate_var_thl_down(i,j) + delta_thl**2 ) ) / candidate_support

          candidate_covar_w_rt_up(i,j) = support_up &
            * ( candidate_covar_w_rt_up(i,j) &
              + ( candidate_w_up(i,j) - w_3(i,j) ) &
              * ( candidate_rt_up(i,j) - rt_3(i,j) ) )
          candidate_covar_w_rt_up(i,j) = &
            ( candidate_covar_w_rt_up(i,j) + support_down &
              * ( candidate_covar_w_rt_down(i,j) &
                + ( candidate_w_down(i,j) - w_3(i,j) ) &
                * ( candidate_rt_down(i,j) - rt_3(i,j) ) ) ) / candidate_support

          candidate_covar_w_thl_up(i,j) = support_up &
            * ( candidate_covar_w_thl_up(i,j) &
              + ( candidate_w_up(i,j) - w_3(i,j) ) &
              * ( candidate_thl_up(i,j) - thl_3(i,j) ) )
          candidate_covar_w_thl_up(i,j) = &
            ( candidate_covar_w_thl_up(i,j) + support_down &
              * ( candidate_covar_w_thl_down(i,j) &
                + ( candidate_w_down(i,j) - w_3(i,j) ) &
                * ( candidate_thl_down(i,j) - thl_3(i,j) ) ) ) / candidate_support

          candidate_covar_rt_thl_up(i,j) = support_up &
            * ( candidate_covar_rt_thl_up(i,j) &
              + ( candidate_rt_up(i,j) - rt_3(i,j) ) &
              * ( candidate_thl_up(i,j) - thl_3(i,j) ) )
          candidate_covar_rt_thl_up(i,j) = &
            ( candidate_covar_rt_thl_up(i,j) + support_down &
              * ( candidate_covar_rt_thl_down(i,j) &
                + ( candidate_rt_down(i,j) - rt_3(i,j) ) &
                * ( candidate_thl_down(i,j) - thl_3(i,j) ) ) ) / candidate_support

          denom = sqrt( max( varnce_w_3(i,j) * varnce_rt_3(i,j), 0._core_rknd ) )
          if ( denom > 0._core_rknd ) then
            corr_w_rt_3(i,j) = max( -1._core_rknd, &
              min( 1._core_rknd, candidate_covar_w_rt_up(i,j) / denom ) )
          else
            corr_w_rt_3(i,j) = 0._core_rknd
          end if

          denom = sqrt( max( varnce_w_3(i,j) * varnce_thl_3(i,j), 0._core_rknd ) )
          if ( denom > 0._core_rknd ) then
            corr_w_thl_3(i,j) = max( -1._core_rknd, &
              min( 1._core_rknd, candidate_covar_w_thl_up(i,j) / denom ) )
          else
            corr_w_thl_3(i,j) = 0._core_rknd
          end if

          denom = sqrt( max( varnce_rt_3(i,j) * varnce_thl_3(i,j), 0._core_rknd ) )
          if ( denom > 0._core_rknd ) then
            corr_rt_thl_3(i,j) = max( -1._core_rknd, &
              min( 1._core_rknd, candidate_covar_rt_thl_up(i,j) / denom ) )
          else
            corr_rt_thl_3(i,j) = 0._core_rknd
          end if
          pdf9_candidate_valid(i,j) = 1._core_rknd
        else
          w_3(i,j) = 0._core_rknd
          rt_3(i,j) = rtm(i,j)
          thl_3(i,j) = thlm(i,j)
          varnce_w_3(i,j) = 0._core_rknd
          varnce_rt_3(i,j) = 0._core_rknd
          varnce_thl_3(i,j) = 0._core_rknd
          corr_w_rt_3(i,j) = 0._core_rknd
          corr_w_thl_3(i,j) = 0._core_rknd
          corr_rt_thl_3(i,j) = 0._core_rknd
          candidate_donor_distance_combined(i,j) = 0._core_rknd
          pdf9_candidate_valid(i,j) = 0._core_rknd
        end if
      end do
    end do

    if ( stats%l_sample ) then
      call stats_update( "pdf9_lscale_down", pdf9_lscale_down, stats )
      call stats_update( "pdf9_down_parcel_tke", parcel_tke, stats )
      call stats_update( "pdf9_down_parcel_w", parcel_w, stats )
      call stats_update( "pdf9_down_parcel_rt", parcel_rt, stats )
      call stats_update( "pdf9_down_parcel_thl", parcel_thl, stats )
      call stats_update( "pdf9_down_parcel_buoyancy", parcel_buoyancy, stats )
      call stats_update( "pdf9_down_parcel_status", parcel_status, stats )
      call stats_update( "pdf9_down_crossing_weight", crossing_weight, stats )
      call stats_update( "pdf9_down_crossing_mean_w", crossing_mean_w, stats )
      call stats_update( "pdf9_down_crossing_mean_rt", crossing_mean_rt, stats )
      call stats_update( "pdf9_down_crossing_mean_thl", crossing_mean_thl, stats )
      call stats_update( "pdf9_down_crossing_var_w", crossing_var_w, stats )
      call stats_update( "pdf9_down_crossing_var_rt", crossing_var_rt, stats )
      call stats_update( "pdf9_down_crossing_var_thl", crossing_var_thl, stats )
      call stats_update( "pdf9_down_crossing_covar_w_rt", crossing_covar_w_rt, stats )
      call stats_update( "pdf9_down_crossing_covar_w_thl", crossing_covar_w_thl, stats )
      call stats_update( "pdf9_down_crossing_covar_rt_thl", crossing_covar_rt_thl, stats )
      call stats_update( "pdf9_candidate_valid_down", candidate_valid_down, stats )
      call stats_update( "pdf9_candidate_crossing_count_down", &
                         candidate_crossing_count_down, stats )
      call stats_update( "pdf9_candidate_weighted_support_down", &
                         candidate_weighted_support_down, stats )
      call stats_update( "pdf9_candidate_donor_distance_down", &
                         candidate_donor_distance_down, stats )
      call stats_update( "pdf9_candidate_w_down", candidate_w_down, stats )
      call stats_update( "pdf9_candidate_w_down_uncapped", &
                         candidate_w_down_uncapped, stats )
      call stats_update( "pdf9_candidate_down_cap_fraction", &
                         candidate_down_cap_fraction, stats )
      call stats_update( "pdf9_candidate_destination_sigma_w", &
                         candidate_destination_sigma_w, stats )
      call stats_update( "pdf9_candidate_rt_down", candidate_rt_down, stats )
      call stats_update( "pdf9_candidate_thl_down", candidate_thl_down, stats )
      call stats_update( "pdf9_candidate_var_w_down", candidate_var_w_down, stats )
      call stats_update( "pdf9_candidate_var_rt_down", candidate_var_rt_down, stats )
      call stats_update( "pdf9_candidate_var_thl_down", candidate_var_thl_down, stats )
      call stats_update( "pdf9_candidate_covar_w_rt_down", &
                         candidate_covar_w_rt_down, stats )
      call stats_update( "pdf9_candidate_covar_w_thl_down", &
                         candidate_covar_w_thl_down, stats )
      call stats_update( "pdf9_candidate_covar_rt_thl_down", &
                         candidate_covar_rt_thl_down, stats )
      call stats_update( "pdf9_candidate_corr_w_rt_down", &
                         candidate_corr_w_rt_down, stats )
      call stats_update( "pdf9_candidate_corr_w_thl_down", &
                         candidate_corr_w_thl_down, stats )
      call stats_update( "pdf9_candidate_corr_rt_thl_down", &
                         candidate_corr_rt_thl_down, stats )
      call stats_update( "pdf9_candidate_valid_combined", &
                         pdf9_candidate_valid, stats )
      call stats_update( "pdf9_candidate_crossing_count_combined", &
                         candidate_crossing_count_combined, stats )
      call stats_update( "pdf9_candidate_weighted_support_combined", &
                         candidate_weighted_support_combined, stats )
      call stats_update( "pdf9_candidate_donor_distance_combined", &
                         candidate_donor_distance_combined, stats )
      call stats_update( "pdf9_candidate_w_combined", w_3, stats )
      call stats_update( "pdf9_candidate_rt_combined", rt_3, stats )
      call stats_update( "pdf9_candidate_thl_combined", thl_3, stats )
      call stats_update( "pdf9_candidate_var_w_combined", varnce_w_3, stats )
      call stats_update( "pdf9_candidate_var_rt_combined", varnce_rt_3, stats )
      call stats_update( "pdf9_candidate_var_thl_combined", varnce_thl_3, stats )
      call stats_update( "pdf9_candidate_covar_w_rt_combined", &
                         candidate_covar_w_rt_up, stats )
      call stats_update( "pdf9_candidate_covar_w_thl_combined", &
                         candidate_covar_w_thl_up, stats )
      call stats_update( "pdf9_candidate_covar_rt_thl_combined", &
                         candidate_covar_rt_thl_up, stats )
      call stats_update( "pdf9_candidate_corr_w_rt_combined", corr_w_rt_3, stats )
      call stats_update( "pdf9_candidate_corr_w_thl_combined", &
                         corr_w_thl_3, stats )
      call stats_update( "pdf9_candidate_corr_rt_thl_combined", &
                         corr_rt_thl_3, stats )
    end if


    return

  end subroutine diagnose_pdf9_mixing_reach

end module pdf_9_module
