!-----------------------------------------------------------------------
!  $Id$
!===============================================================================
module stats_clubb_utilities

  implicit none

  private ! Set Default Scope
  public :: stats_accumulate, stats_accumulate_hydromet_api, stats_accumulate_lh_tend

contains

    subroutine stats_accumulate( &
                       nzm, nzt, ngrdcol, sclr_dim, edsclr_dim, &
                       gr, dt, &
                     l_implemented, l_host_applies_sfc_fluxes, &
                     l_stability_correct_tau_zm, &
                     clubb_params, &
                     um, vm, upwp, vpwp, up2, vp2, &
                     thlm, rtm, thlm_before, rtm_before, thlm_forcing, rtm_forcing, &
                     wpthlp_sfc, wprtp_sfc, wprtp, wpthlp, &
                     wp2, wp3, rtp2, rtp3, thlp2, thlp3, &
                     rtpthlp, &
                     p_in_Pa, exner, rho, rho_zm, &
                     rho_ds_zm, rho_ds_zt, thv_ds_zm, thv_ds_zt, &
                     wm_zt, wm_zm, rcm, &
                     cloud_frac, &
                     thvm, ug, vg, &
                     ddzt_umvm_sqd, stability_correction, &
                     Kh_zt, &
                     rsat, &
                     Kh_zm, &
                     em, &
                     sclrm, sclrp2, &
                     sclrprtp, sclrpthlp, sclrm_forcing, &
                     wpsclrp, wpedsclrp, edsclrm, &
                     edsclrm_forcing, &
                     saturation_formula, &
                     stats )

    use constants_clubb, only: &
          cloud_frac_min, &  ! Constant
          eps, &
          w_tol

    use T_in_K_module, only: & 
        thlm2T_in_K_api ! Procedure

    use constants_clubb, only: &
        rc_tol    ! Constant(s)

    use grid_class, only: &
        grid, & ! Type
        zt2zm_api

    use parameter_indices, only: &
        nparams

    use Skx_module, only: &
        Skx_func, &
        compute_gamma_Skw

    use stats_netcdf, only: &
        stats_type, &
        stats_update, &
        var_on_stats_list

    use advance_helper_module, only: &
        calc_wp3_on_wp2, &
        vertical_avg, &     ! Procedure(s)
        vertical_integral

    use interpolation, only: &
        lin_interpolate_two_points             ! Procedure

    use saturation, only: &
        sat_mixrat_ice ! Procedure

    use numerical_check, only: &
        calculate_spurious_source

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    integer, intent(in) :: &
      nzm, nzt, ngrdcol, sclr_dim, edsclr_dim

    type(grid), intent(in) :: &
      gr

    real( kind = core_rknd ), intent(in) :: &
      dt

    logical, intent(in) :: &
      l_implemented, &
      l_host_applies_sfc_fluxes, &
      l_stability_correct_tau_zm

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nparams) :: &
      clubb_params

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) :: &
      um, vm, thlm, rtm, thlm_before, rtm_before, thlm_forcing, rtm_forcing, &
      wp3, rtp3, thlp3

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) :: &
      upwp, vpwp, up2, vp2, wprtp, wpthlp, wp2, rtp2, thlp2, rtpthlp

    real( kind = core_rknd ), intent(in), dimension(ngrdcol) :: &
      wpthlp_sfc, &
      wprtp_sfc

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) :: &
      p_in_Pa, exner, rho, rho_ds_zt, thv_ds_zt, wm_zt

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) :: &
      rho_zm, rho_ds_zm, thv_ds_zm, wm_zm

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) :: &
      rcm, cloud_frac

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) :: &
      thvm, ug, vg, Kh_zt, rsat

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) :: &
      ddzt_umvm_sqd, stability_correction, Kh_zm, em

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt,sclr_dim) :: &
      sclrm, sclrm_forcing

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm,sclr_dim) :: &
      sclrp2, sclrprtp, sclrpthlp, wpsclrp

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt,edsclr_dim) :: &
      edsclrm, edsclrm_forcing

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm,edsclr_dim) :: &
      wpedsclrp

    integer, intent(in) :: &
      saturation_formula

    type(stats_type), intent(inout) :: &
      stats

    ! Local Variables

    integer :: sclr, edsclr, k, i
      
      real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
        T_in_K,        &  ! Absolute temperature         [K]
        rsati,         &  ! Saturation w.r.t ice         [kg/kg]
        rcm_in_cloud     ! rcm in cloud                 [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      shear       ! Wind shear production term   [m^2/s^3]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      dzm         ! Signed momentum-grid layer thickness [m]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      dzt         ! Signed thermodynamic-grid layer thickness [m]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      wp3_on_wp2    ! Ratio of wp3 to wp2 on momentum levels [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      wp3_zm,         & ! Third moment of w on momentum levels [m^3/s^3]
      Skw_zm,         & ! Skewness of w on momentum levels [-]
      gamma_Skw_fnc     ! Gamma as a function of w skewness [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      wp3_on_wp2_zt    ! Ratio of wp3 to wp2 on thermodynamic levels [m/s]

    real( kind = core_rknd ), dimension(ngrdcol):: xtmp

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      rtm_integral_before, &
      rtm_integral_after, &
      rtm_integral_forcing, &
      rtm_flux_top, &
      rtm_flux_sfc, &
      rtm_spur_src, &
      thlm_integral_before, &
      thlm_integral_after, &
      thlm_integral_forcing, &
      thlm_flux_top, &
      thlm_flux_sfc, &
      thlm_spur_src

    character(len=50) :: sclr_idx

    ! ---- Begin Code ----

    ! Sample fields

      if ( stats%l_sample ) then

        dzm(:,:) = gr%grid_dir * gr%dzm(:,:)
        dzt(:,:) = gr%grid_dir * gr%dzt(:,:)

        !$acc enter data create( wp3_on_wp2, wp3_on_wp2_zt )

        call calc_wp3_on_wp2( nzm, nzt, ngrdcol, gr, wp2, wp3, &
                              wp3_on_wp2, wp3_on_wp2_zt )

        !$acc update host( wp2, vp2, up2, wprtp, wpthlp, upwp, vpwp, rtp2, thlp2, &
        !$acc              rtpthlp, rtm, thlm, um, vm, wp3, &
        !$acc              um, vm, upwp, vpwp, up2, vp2, &
        !$acc              thlm, rtm, thlm_before, rtm_before, thlm_forcing, rtm_forcing, &
        !$acc              wpthlp_sfc, wprtp_sfc, wprtp, wpthlp, &
        !$acc              wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &
        !$acc              p_in_Pa, exner, rho, rho_zm, &
        !$acc              rho_ds_zm, rho_ds_zt, thv_ds_zm, thv_ds_zt, &
        !$acc              wm_zt, wm_zm, rcm, cloud_frac, &
        !$acc              thvm, ug, vg, ddzt_umvm_sqd, Kh_zt, &
        !$acc              Kh_zm, &
        !$acc              em, wp3_on_wp2, wp3_on_wp2_zt )

        !$acc update host( stability_correction ) if ( l_stability_correct_tau_zm )

        !$acc update host( sclrm, sclrp2, sclrprtp, sclrpthlp, sclrm_forcing, &
        !$acc              wpsclrp, wpedsclrp ) &
        !$acc if ( sclr_dim > 0 )

        !$acc update host( edsclrm, edsclrm_forcing  ) if ( edsclr_dim > 0 )

        if ( var_on_stats_list( stats, "T_in_K" ) .or.       &
            var_on_stats_list( stats, "rsati" ) ) then
        T_in_K = thlm2T_in_K_api( nzt, ngrdcol, thlm, exner, rcm )
        call stats_update( "T_in_K", T_in_K, stats )
      end if

      if ( var_on_stats_list( stats, "rsati" ) ) then
        rsati = sat_mixrat_ice( nzt, ngrdcol, p_in_Pa, T_in_K, saturation_formula )
        call stats_update( "rsati", rsati, stats )
      end if

      if ( var_on_stats_list( stats, "rcm_in_cloud" ) ) then
        where ( cloud_frac > cloud_frac_min )
            rcm_in_cloud = rcm / cloud_frac
        elsewhere
            rcm_in_cloud = rcm
        endwhere
        call stats_update( "rcm_in_cloud", rcm_in_cloud, stats )
      end if

      if ( var_on_stats_list( stats, "shear" ) ) then
        shear(:,1)   = 0.0_core_rknd
        do k = 2, nzm-1, 1
          do i = 1, ngrdcol
            shear(i,k) = - upwp(i,k) * ( um(i,k) - um(i,k-1) ) * gr%invrs_dzm(i,k)  &
                        - vpwp(i,k) * ( vm(i,k) - vm(i,k-1) ) * gr%invrs_dzm(i,k)
          end do
        enddo
        shear(:,nzm) = 0.0_core_rknd
        call stats_update( "shear", shear, stats )
      end if

      ! stats_zt variables
      call stats_update( "thlm", thlm, stats )
      call stats_update( "thvm", thvm, stats )
      call stats_update( "rtm", rtm, stats )
      call stats_update( "rcm", rcm, stats )
      call stats_update( "um", um, stats )
      call stats_update( "vm", vm, stats )
      call stats_update( "wm_zt", wm_zt, stats )
      call stats_update( "ug", ug, stats )
      call stats_update( "vg", vg, stats )
      call stats_update( "p_in_Pa", p_in_Pa, stats )
      call stats_update( "exner", exner, stats )
      call stats_update( "rho_ds_zt", rho_ds_zt, stats )
      call stats_update( "thv_ds_zt", thv_ds_zt, stats )
      call stats_update( "wp3", wp3, stats )
      call stats_update( "Kh_zt", Kh_zt, stats )
      call stats_update( "rho", rho, stats )
      call stats_update( "rsat", rsat, stats )
      call stats_update( "thlp3", thlp3, stats )
      call stats_update( "rtp3", rtp3, stats )
      call stats_update( "wp3_on_wp2_zt", wp3_on_wp2_zt, stats )

      if ( sclr_dim > 0 ) then
        do sclr = 1, sclr_dim
          write( sclr_idx, * ) sclr
          sclr_idx = adjustl( sclr_idx )
          call stats_update( "sclr"//trim(sclr_idx)//"m", sclrm(:,:,sclr), stats )
          call stats_update( "sclr"//trim(sclr_idx)//"m_f", sclrm_forcing(:,:,sclr), stats )
        end do
      end if
      
      if ( edsclr_dim > 0 ) then
        do edsclr = 1, edsclr_dim
          write( sclr_idx, * ) edsclr
          sclr_idx = adjustl( sclr_idx )
          call stats_update( "edsclr"//trim(sclr_idx)//"m", edsclrm(:,:,edsclr), stats )
          call stats_update( "edsclr"//trim(sclr_idx)//"m_f", edsclrm_forcing(:,:,edsclr), stats )
        end do
      end if

      ! stats_zm variables
      call stats_update( "wm_zm", wm_zm, stats )
      call stats_update( "ddzt_umvm_sqd", ddzt_umvm_sqd, stats )
      call stats_update( "wp2", wp2, stats )
      call stats_update( "rtp2", rtp2, stats )
      call stats_update( "thlp2", thlp2, stats )
      call stats_update( "rtpthlp", rtpthlp, stats )
      call stats_update( "wprtp", wprtp, stats )
      call stats_update( "wpthlp", wpthlp, stats )
      if ( l_stability_correct_tau_zm ) then
        call stats_update( "stability_correction", stability_correction, stats )
      end if
      call stats_update( "Kh_zm", Kh_zm, stats )
      call stats_update( "upwp", upwp, stats )
      call stats_update( "vpwp", vpwp, stats )
      call stats_update( "vp2", vp2, stats )
      call stats_update( "up2", up2, stats )
      call stats_update( "rho_zm", rho_zm, stats )
      call stats_update( "rho_ds_zm", rho_ds_zm, stats )
      call stats_update( "thv_ds_zm", thv_ds_zm, stats )
      call stats_update( "em", em, stats )
      call stats_update( "wp3_on_wp2", wp3_on_wp2, stats )
      call stats_update( "wp3_on_wp2_cfl_num", wp3_on_wp2 * dt / dzm, stats )
      if ( var_on_stats_list( stats, "gamma_Skw_fnc" ) ) then
        !$acc data create( Skw_zm, wp3_zm ) copyout( gamma_Skw_fnc )
        wp3_zm(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, wp3(:,:) )
        call Skx_func( nzm, ngrdcol, wp2, wp3_zm, &
                       w_tol, clubb_params, &
                       Skw_zm )
        call compute_gamma_Skw( nzm, ngrdcol, Skw_zm, clubb_params, & ! In
                                gamma_Skw_fnc )                       ! Out
        !$acc end data
        call stats_update( "gamma_Skw_fnc", gamma_Skw_fnc, stats )
      end if
      if ( sclr_dim > 0 ) then
        do sclr = 1, sclr_dim
          write( sclr_idx, * ) sclr
          sclr_idx = adjustl( sclr_idx )
          call stats_update( "sclr"//trim(sclr_idx)//"p2", sclrp2(:,:,sclr), stats )
          call stats_update( "sclr"//trim(sclr_idx)//"prtp", sclrprtp(:,:,sclr), stats )
          call stats_update( "sclr"//trim(sclr_idx)//"pthlp", sclrpthlp(:,:,sclr), stats )
          call stats_update( "wpsclr"//trim(sclr_idx)//"p", wpsclrp(:,:,sclr), stats )
        end do
      end if
      if ( edsclr_dim > 0 ) then
        do edsclr = 1, edsclr_dim
          write( sclr_idx, * ) edsclr
          sclr_idx = adjustl( sclr_idx )
          call stats_update( "wpedsclr"//trim(sclr_idx)//"p", wpedsclrp(:,:,edsclr), stats )
        end do
      end if

      ! stats_sfc variables
      do i = 1, ngrdcol
        call stats_update( "cc", maxval( cloud_frac(i,:) ), stats, i )
      end do

      if ( var_on_stats_list( stats, "z_cloud_base" ) ) then
        do i = 1, ngrdcol
          k = 1
          do while ( rcm(i,k) < rc_tol .and. k < nzt )
            k = k + 1
          enddo
          if ( k == 1 ) then
            xtmp(i) = gr%zt(i,1)
          elseif ( k > 1 .and. k < nzt ) then
            xtmp(i) = lin_interpolate_two_points( rc_tol, rcm(i,k), rcm(i,k-1), &
                                                  gr%zt(i,k), gr%zt(i,k-1) )
          else
            xtmp(i) = -10.0_core_rknd
          end if
        end do
        call stats_update( "z_cloud_base", xtmp, stats )
      end if

      if ( var_on_stats_list( stats, "lwp" ) ) then
        do i = 1, ngrdcol
          xtmp(i) = vertical_integral( nzt, rho_ds_zt(i,:), rcm(i,:), dzt(i,:) )
        end do
        call stats_update( "lwp", xtmp, stats )
      end if

      if ( var_on_stats_list( stats, "vwp" ) ) then
        do i = 1, ngrdcol
          xtmp(i) = vertical_integral( nzt, rho_ds_zt(i,:), &
                                    ( rtm(i,:) - rcm(i,:) ), dzt(i,:) )
        end do
        call stats_update( "vwp", xtmp, stats )
      end if

      do i = 1, ngrdcol
        call stats_update( "thlm_vert_avg", &
              vertical_avg( nzt, rho_ds_zt(i,:), thlm(i,:), dzt(i,:) ), stats, i)
        call stats_update( "rtm_vert_avg", &
              vertical_avg( nzt, rho_ds_zt(i,:), rtm(i,:), dzt(i,:) ), stats, i)
        call stats_update( "um_vert_avg", &
              vertical_avg( nzt, rho_ds_zt(i,:), um(i,:), dzt(i,:) ), stats, i)
        call stats_update( "vm_vert_avg", &
              vertical_avg( nzt, rho_ds_zt(i,:), vm(i,:), dzt(i,:) ), stats, i)
        call stats_update( "wp2_vert_avg", &
              vertical_avg( nzm, rho_ds_zm(i,:), wp2(i,:), dzm(i,:) ), stats, i)
        call stats_update( "up2_vert_avg", &
              vertical_avg( nzm, rho_ds_zm(i,:), up2(i,:), dzm(i,:) ), stats, i)
        call stats_update( "vp2_vert_avg", &
              vertical_avg( nzm, rho_ds_zm(i,:), vp2(i,:), dzm(i,:) ), stats, i)
        call stats_update( "rtp2_vert_avg", &
              vertical_avg( nzm, rho_ds_zm(i,:), rtp2(i,:), dzm(i,:) ), stats, i)
        call stats_update( "thlp2_vert_avg", &
              vertical_avg( nzm, rho_ds_zm(i,:), thlp2(i,:), dzm(i,:) ), stats, i)
      end do

      if ( var_on_stats_list( stats, "tot_vartn_normlzd_rtm" ) ) then
        do i = 1, ngrdcol
          if ( abs(rtm(i,nzt) - rtm(i,1)) < eps ) then
            xtmp(i) = -999_core_rknd
          else
            xtmp(i) = sum(abs(rtm(i,2:nzt) - rtm(i,1:nzt-1)) / abs(rtm(i,nzt) - rtm(i,1)))
          end if
        end do
        call stats_update( "tot_vartn_normlzd_rtm", xtmp, stats )
      end if

      if ( var_on_stats_list( stats, "tot_vartn_normlzd_thlm" ) ) then
        do i = 1, ngrdcol
          if ( abs(thlm(i,nzt) - thlm(i,1)) < eps ) then
            xtmp(i) = -999_core_rknd
          else
            xtmp(i) = sum(abs(thlm(i,2 : nzt) - thlm(i,1 : nzt-1)) / abs(thlm(i,nzt) - thlm(i,1)))
          end if
        end do
        call stats_update( "tot_vartn_normlzd_thlm", xtmp, stats )
      end if

      if ( var_on_stats_list( stats, "tot_vartn_normlzd_wprtp" ) ) then
        do i = 1, ngrdcol
          if ( abs(wprtp(i,nzm) - wprtp(i,1)) < eps ) then
            xtmp(i) = -999_core_rknd
          else
            xtmp(i) = sum(abs(wprtp(i,2:nzm) - wprtp(i,1:nzm-1)) / abs(wprtp(i,nzm) - wprtp(i,1)))
          end if
        end do
        call stats_update( "tot_vartn_normlzd_wprtp", xtmp, stats )
      end if

      ! Spurious source will only be calculated if rtm_ma and thlm_ma are zero.
      ! Therefore, wm must be zero or l_implemented must be true.
      do i = 1, ngrdcol
        if ( l_implemented .or. &
             (all( abs(wm_zt(i,:)) < eps ) .and. all( abs(wm_zm(i,:)) < eps ))) then
          rtm_flux_top(i) = rho_ds_zm(i,gr%k_ub_zm) * wprtp(i,gr%k_ub_zm)

          if ( .not. l_host_applies_sfc_fluxes ) then
            rtm_flux_sfc(i) = rho_ds_zm(i,gr%k_lb_zm) * wprtp_sfc(i)
          else
            rtm_flux_sfc(i) = 0.0_core_rknd
          end if

          ! Get the vertical integral of rtm before this function begins
          ! so that spurious source can be calculated.
          rtm_integral_before(i)  &
          = vertical_integral( nzt, rho_ds_zt(i,:), rtm_before(i,:), dzt(i,:) )

          rtm_integral_after(i)  &
          = vertical_integral( nzt, rho_ds_zt(i,:), rtm(i,:), dzt(i,:) )

          rtm_integral_forcing(i)  &
          = vertical_integral( nzt, rho_ds_zt(i,:), rtm_forcing(i,:), dzt(i,:) )

          ! Calculate the spurious source for rtm.
          rtm_spur_src(i)  &
          = calculate_spurious_source( rtm_integral_after(i), &
                                       rtm_integral_before(i), &
                                       rtm_flux_top(i), rtm_flux_sfc(i), &
                                       rtm_integral_forcing(i), &
                                       dt )

          thlm_flux_top(i) = rho_ds_zm(i,gr%k_ub_zm) * wpthlp(i,gr%k_ub_zm)

          if ( .not. l_host_applies_sfc_fluxes ) then
            thlm_flux_sfc(i) = rho_ds_zm(i,gr%k_lb_zm) * wpthlp_sfc(i)
          else
            thlm_flux_sfc(i) = 0.0_core_rknd
          end if

          ! Get the vertical integral of thlm before this function begins
          ! so that spurious source can be calculated.
          thlm_integral_before(i)  &
          = vertical_integral( nzt, rho_ds_zt(i,:), thlm_before(i,:), dzt(i,:) )

          thlm_integral_after(i)  &
          = vertical_integral( nzt, rho_ds_zt(i,:), thlm(i,:), dzt(i,:) )

          thlm_integral_forcing(i)  &
          = vertical_integral( nzt, rho_ds_zt(i,:), thlm_forcing(i,:), dzt(i,:) )

          ! Calculate the spurious source for thlm.
          thlm_spur_src(i)  &
          = calculate_spurious_source( thlm_integral_after(i), &
                                       thlm_integral_before(i), &
                                       thlm_flux_top(i), thlm_flux_sfc(i), &
                                       thlm_integral_forcing(i), &
                                       dt )
        else ! If l_implemented is false, we don't want spurious source output
          rtm_spur_src(i) = -9999.0_core_rknd
          thlm_spur_src(i) = -9999.0_core_rknd
        end if
      end do
      call stats_update( "rtm_spur_src", rtm_spur_src, stats )
      call stats_update( "thlm_spur_src", thlm_spur_src, stats )

      !$acc exit data delete( wp3_on_wp2, wp3_on_wp2_zt )
    end if

    return
  end subroutine stats_accumulate

!------------------------------------------------------------------------------
  subroutine stats_accumulate_hydromet_api( gr, hydromet_dim, hm_metadata, & ! intent(in)
                                            hydromet, rho_ds_zt,           & ! intent(in)
                                            stats, icol )                    ! intent(inout/in)
! Description:
!   Compute stats related the hydrometeors

! References:
!   None
!------------------------------------------------------------------------------

    use grid_class, only: &
        grid ! Type
        
    use corr_varnce_module, only: &
        hm_metadata_type

    use advance_helper_module, only: &
        vertical_integral ! Procedure(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use stats_netcdf, only: &
        stats_type, &
        stats_update, &
        var_on_stats_list

    implicit none

    type (grid), intent(in) :: &
      gr

    integer, intent(in) :: &
      hydromet_dim

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    type(stats_type), intent(inout) :: &
      stats

    integer, intent(in) :: &
      icol

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(in) :: &
      hydromet ! All hydrometeors except for rcm        [units vary]

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: &
      rho_ds_zt ! Dry, static density (thermo. levs.)      [kg/m^3]

    ! Local Variables
    real(kind=core_rknd) :: xtmp
    
    ! ---- Begin Code ----
    if ( stats%l_sample ) then

      if ( hm_metadata%iirr > 0 ) then
        call stats_update( "rrm", hydromet(:,hm_metadata%iirr), stats, icol )
      end if

      if ( hm_metadata%iirs > 0 ) then
        call stats_update( "rsm", hydromet(:,hm_metadata%iirs), stats, icol )
      end if 

      if ( hm_metadata%iiri > 0 ) then 
        call stats_update( "rim", hydromet(:,hm_metadata%iiri), stats, icol )
      end if

      if ( hm_metadata%iirg > 0 ) then
        call stats_update( "rgm", hydromet(:,hm_metadata%iirg), stats, icol )
      end if

      if ( hm_metadata%iiNi > 0 ) then
        call stats_update( "Nim", hydromet(:,hm_metadata%iiNi), stats, icol )
      end if

      if ( hm_metadata%iiNr > 0 ) then
        call stats_update( "Nrm", hydromet(:,hm_metadata%iiNr), stats, icol )
      end if

      if ( hm_metadata%iiNs > 0 ) then
        call stats_update( "Nsm", hydromet(:,hm_metadata%iiNs), stats, icol )
      end if

      if ( hm_metadata%iiNg > 0 ) then
        call stats_update( "Ngm", hydromet(:,hm_metadata%iiNg), stats, icol )
      end if

      ! Snow Water Path
      if ( var_on_stats_list( stats, "swp" ) .and. hm_metadata%iirs > 0 ) then
        xtmp &
        = vertical_integral &
               ( gr%nzt, rho_ds_zt, &
                 hydromet(:,hm_metadata%iirs), gr%dzt(1,:) )

        call stats_update( "swp", xtmp, stats, icol )
      end if

      ! Ice Water Path
      if ( var_on_stats_list( stats, "iwp" ) .and. hm_metadata%iiri > 0 ) then
        xtmp &
        = vertical_integral &
               ( gr%nzt, rho_ds_zt, &
                 hydromet(:,hm_metadata%iiri), gr%dzt(1,:) )

        call stats_update( "iwp", xtmp, stats, icol )
      end if

      ! Rain Water Path
      if ( var_on_stats_list( stats, "rwp" ) .and. hm_metadata%iirr > 0 ) then
        xtmp &
        = vertical_integral &
               ( gr%nzt, rho_ds_zt, &
                 hydromet(:,hm_metadata%iirr), gr%dzt(1,:) )

        call stats_update( "rwp", xtmp, stats, icol )
      end if
    end if

    return
  end subroutine stats_accumulate_hydromet_api
!------------------------------------------------------------------------------
  subroutine stats_accumulate_lh_tend( gr, hydromet_dim, hm_metadata, &
                                       lh_hydromet_mc, lh_Ncm_mc, &
                                       lh_thlm_mc, lh_rvm_mc, lh_rcm_mc, &
                                       lh_AKm, AKm, AKstd, AKstd_cld, &
                                       lh_rcm_avg, AKm_rcm, AKm_rcc, &
                                       stats, icol )

! Description:
!   Compute stats for the tendency of latin hypercube sample points.

! References:
!   None
!------------------------------------------------------------------------------

    use grid_class, only: &
        grid ! Type

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use corr_varnce_module, only: &
        hm_metadata_type

    use stats_netcdf, only: &
        stats_type, &
        stats_update

    implicit none

    !----------------------- Input Variables -----------------------
    type (grid), intent(in) :: &
      gr

    integer, intent(in) :: &
      hydromet_dim

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(in) :: &
      lh_hydromet_mc ! Tendency of hydrometeors except for rvm/rcm  [units vary]

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: &
      lh_Ncm_mc,  & ! Tendency of cloud droplet concentration  [num/kg/s]
      lh_thlm_mc, & ! Tendency of liquid potential temperature [kg/kg/s]
      lh_rcm_mc,  & ! Tendency of cloud water                  [kg/kg/s]
      lh_rvm_mc     ! Tendency of vapor                        [kg/kg/s]

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: &
      lh_AKm,     & ! Kessler ac estimate                 [kg/kg/s]
      AKm,        & ! Exact Kessler ac                    [kg/kg/s]
      AKstd,      & ! St dev of exact Kessler ac          [kg/kg/s]
      AKstd_cld,  & ! Stdev of exact w/in cloud ac        [kg/kg/s]
      lh_rcm_avg, & ! Monte Carlo rcm estimate            [kg/kg]
      AKm_rcm,    & ! Kessler ac based on rcm             [kg/kg/s]
      AKm_rcc       ! Kessler ac based on rcm/cloud_frac  [kg/kg/s]
      
    type(stats_type), intent(inout) :: &
      stats

    integer, intent(in) :: &
      icol

    !----------------------- Local Variables -----------------------

    !----------------------- Begin Code -----------------------

    if ( stats%l_sample ) then
      call stats_update( "lh_thlm_mc", lh_thlm_mc(:), stats, icol )
      call stats_update( "lh_rcm_mc", lh_rcm_mc(:), stats, icol )
      call stats_update( "lh_rvm_mc", lh_rvm_mc(:), stats, icol )
      call stats_update( "lh_Ncm_mc", lh_Ncm_mc(:), stats, icol )

      if ( hm_metadata%iirr > 0 ) then
        call stats_update( "lh_rrm_mc", lh_hydromet_mc(:,hm_metadata%iirr), stats, icol )
      end if
      if ( hm_metadata%iirs > 0 ) then
        call stats_update( "lh_rsm_mc", lh_hydromet_mc(:,hm_metadata%iirs), stats, icol )
      end if
      if ( hm_metadata%iiri > 0 ) then
        call stats_update( "lh_rim_mc", lh_hydromet_mc(:,hm_metadata%iiri), stats, icol )
      end if
      if ( hm_metadata%iirg > 0 ) then
        call stats_update( "lh_rgm_mc", lh_hydromet_mc(:,hm_metadata%iirg), stats, icol )
      end if
      if ( hm_metadata%iiNi > 0 ) then
        call stats_update( "lh_Nim_mc", lh_hydromet_mc(:,hm_metadata%iiNi), stats, icol )
      end if
      if ( hm_metadata%iiNr > 0 ) then
        call stats_update( "lh_Nrm_mc", lh_hydromet_mc(:,hm_metadata%iiNr), stats, icol )
      end if
      if ( hm_metadata%iiNs > 0 ) then
        call stats_update( "lh_Nsm_mc", lh_hydromet_mc(:,hm_metadata%iiNs), stats, icol )
      end if
      if ( hm_metadata%iiNg > 0 ) then
        call stats_update( "lh_Ngm_mc", lh_hydromet_mc(:,hm_metadata%iiNg), stats, icol )
      end if

      call stats_update( "AKm", AKm(:), stats, icol )
      call stats_update( "lh_AKm", lh_AKm(:), stats, icol )
      call stats_update( "lh_rcm_avg", lh_rcm_avg(:), stats, icol )
      call stats_update( "AKstd", AKstd(:), stats, icol )
      call stats_update( "AKstd_cld", AKstd_cld(:), stats, icol )
      call stats_update( "AKm_rcm", AKm_rcm(:), stats, icol )
      call stats_update( "AKm_rcc", AKm_rcc(:), stats, icol )
    end if

    return
  end subroutine stats_accumulate_lh_tend

end module stats_clubb_utilities
