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
                     invrs_dzm, zt, dzm, dzt, dt, &
                     um, vm, upwp, vpwp, up2, vp2, &
                     thlm, rtm, wprtp, wpthlp, &
                     wp2, wp3, rtp2, rtp3, thlp2, thlp3, &
                     rtpthlp, &
                     wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, &
                     p_in_Pa, exner, rho, rho_zm, &
                     rho_ds_zm, rho_ds_zt, thv_ds_zm, thv_ds_zt, &
                     wm_zt, wm_zm, rcm, wprcp, rc_coef, &
                     rc_coef_zm, &
                     rcm_zm, rtm_zm, thlm_zm, cloud_frac, &
                     ice_supersat_frac, &
                     cloud_frac_zm, ice_supersat_frac_zm, rcm_in_layer, &
                     cloud_cover, rcm_supersat_adj, sigma_sqd_w, &
                     thvm, ug, vg, Lscale, wpthlp2, wp2thlp, &
                     wprtp2, wp2rtp, &
                     Lscale_up, Lscale_down, tau_zt, Kh_zt, wp2rcp, &
                     wprtpthlp, sigma_sqd_w_zt, rsat, wp2_zt, &
                     thlp2_zt, &
                     wpthlp_zt, wprtp_zt, rtp2_zt, rtpthlp_zt, &
                     up2_zt, &
                     vp2_zt, upwp_zt, vpwp_zt, wpup2, wpvp2, &
                     wp2up2, wp2vp2, wp4, &
                     tau_zm, Kh_zm, thlprcp, &
                     rtprcp, rcp2, em, a3_coef, a3_coef_zt, &
                     wp3_zm, wp3_on_wp2, wp3_on_wp2_zt, Skw_velocity, &
                     w_up_in_cloud, w_down_in_cloud, &
                     cloudy_updraft_frac, cloudy_downdraft_frac, &
                     pdf_params, pdf_params_zm, &
                     sclrm, sclrp2, &
                     sclrprtp, sclrpthlp, sclrm_forcing, sclrpthvp, &
                     wpsclrp, sclrprcp, wp2sclrp, wpsclrp2, &
                     wpsclrprtp, &
                     wpsclrpthlp, wpedsclrp, edsclrm, &
                     edsclrm_forcing, &
                     saturation_formula, &
                     l_call_pdf_closure_twice, &
                     stats )

    use constants_clubb, only: &
        cloud_frac_min, &  ! Constant
        eps

    use pdf_utilities, only: &
        compute_variance_binormal    ! Procedure

    use pdf_parameter_module, only: & 
        pdf_parameter ! Type

    use T_in_K_module, only: & 
        thlm2T_in_K_api ! Procedure

    use constants_clubb, only: & 
        rc_tol, fstderr    ! Constant(s)

    use stats_netcdf, only: &
        stats_type, &
        stats_update, &
        var_on_stats_list

    use advance_helper_module, only: &
        vertical_avg, &     ! Procedure(s)
        vertical_integral

    use interpolation, only: & 
        lin_interpolate_two_points             ! Procedure

    use saturation, only: &
        sat_mixrat_ice ! Procedure

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    integer, intent(in) :: &
      nzm, nzt, ngrdcol, sclr_dim, edsclr_dim

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) :: &
      invrs_dzm, dzm

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) :: &
      zt, dzt

    real( kind = core_rknd ), intent(in) :: &
      dt

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) :: &
      um, vm, thlm, rtm, wp3, rtp3, thlp3, wp2thvp, wp2up

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) :: &
      upwp, vpwp, up2, vp2, wprtp, wpthlp, wp2, rtp2, thlp2, rtpthlp, &
      wpthvp, rtpthvp, thlpthvp

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) :: &
      p_in_Pa, exner, rho, rho_ds_zt, thv_ds_zt, wm_zt

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) :: &
      rho_zm, rho_ds_zm, thv_ds_zm, wm_zm

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) :: &
      rcm, rc_coef, cloud_frac, ice_supersat_frac, rcm_in_layer, cloud_cover, rcm_supersat_adj

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) :: &
      rcm_zm, rtm_zm, thlm_zm, wprcp, rc_coef_zm, cloud_frac_zm, ice_supersat_frac_zm

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) :: &
      sigma_sqd_w

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) :: &
      thvm, ug, vg, Lscale, wpthlp2, wp2thlp, wprtp2, wp2rtp, Lscale_up, Lscale_down, tau_zt, &
      Kh_zt, wp2rcp, wprtpthlp, sigma_sqd_w_zt, rsat

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) :: &
      wp2_zt, thlp2_zt, wpthlp_zt, wprtp_zt, rtp2_zt, rtpthlp_zt, up2_zt, vp2_zt, upwp_zt, vpwp_zt, &
      wpup2, wpvp2, a3_coef_zt, wp3_on_wp2_zt, w_up_in_cloud, w_down_in_cloud, cloudy_updraft_frac, &
      cloudy_downdraft_frac

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) :: &
      wp2up2, wp2vp2, wp4, tau_zm, Kh_zm, thlprcp, rtprcp, rcp2, em, a3_coef, wp3_zm, wp3_on_wp2, &
      Skw_velocity

    type(pdf_parameter), intent(in) :: &
      pdf_params, pdf_params_zm

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt,sclr_dim) :: &
      sclrm, sclrm_forcing, wp2sclrp, wpsclrp2, wpsclrprtp, wpsclrpthlp

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm,sclr_dim) :: &
      sclrprcp, sclrp2, sclrprtp, sclrpthlp, sclrpthvp, wpsclrp

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt,edsclr_dim) :: &
      edsclrm, edsclrm_forcing

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm,edsclr_dim) :: &
      wpedsclrp

    integer, intent(in) :: &
      saturation_formula

    logical, intent(in) :: &
      l_call_pdf_closure_twice

    type(stats_type), intent(inout) :: &
      stats

    ! Local Variables

    integer :: sclr, edsclr, k, i
      
    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      T_in_K,        &  ! Absolute temperature         [K]
      rsati,         &  ! Saturation w.r.t ice         [kg/kg]
      chi,           &  ! Mellor's 's'                 [kg/kg]
      chip2,         &  ! Variance of Mellor's 's'     [kg/kg]
      rcm_in_cloud      ! rcm in cloud                 [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      shear       ! Wind shear production term   [m^2/s^3]

    real( kind = core_rknd ), dimension(ngrdcol):: xtmp

    character(len=50) :: sclr_idx

    ! ---- Begin Code ----

    ! Sample fields

    if ( stats%l_sample ) then

      if ( var_on_stats_list( stats, "T_in_K" ) .or.       &
            var_on_stats_list( stats, "rsati" ) ) then
        T_in_K = thlm2T_in_K_api( nzt, ngrdcol, thlm, exner, rcm )
        call stats_update( "T_in_K", T_in_K, stats )
      end if

      if ( var_on_stats_list( stats, "rsati" ) ) then
        rsati = sat_mixrat_ice( nzt, ngrdcol, p_in_Pa, T_in_K, saturation_formula )
        call stats_update( "rsati", rsati, stats )
      end if

      if ( var_on_stats_list( stats, "chi" ) .or.       &
            var_on_stats_list( stats, "chip2" ) ) then
        chi = pdf_params%mixt_frac * pdf_params%chi_1 &
                    + (1.0_core_rknd-pdf_params%mixt_frac) * pdf_params%chi_2
        call stats_update( "chi", chi, stats )
      end if
      
      if ( var_on_stats_list( stats, "chip2" ) ) then
        chip2 = compute_variance_binormal( chi, pdf_params%chi_1, pdf_params%chi_2, &
                                           pdf_params%stdev_chi_1, pdf_params%stdev_chi_2, &
                                           pdf_params%mixt_frac )
        call stats_update( "chip2", chip2, stats )
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
            shear(i,k) = - upwp(i,k) * ( um(i,k) - um(i,k-1) ) * invrs_dzm(i,k)  &
                        - vpwp(i,k) * ( vm(i,k) - vm(i,k-1) ) * invrs_dzm(i,k)
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
      call stats_update( "cloud_frac", cloud_frac, stats )
      call stats_update( "ice_supersat_frac", ice_supersat_frac, stats )
      call stats_update( "rcm_in_layer", rcm_in_layer, stats )
      call stats_update( "cloud_cover", cloud_cover, stats )
      call stats_update( "rcm_supersat_adj", rcm_supersat_adj, stats )
      call stats_update( "p_in_Pa", p_in_Pa, stats )
      call stats_update( "exner", exner, stats )
      call stats_update( "rho_ds_zt", rho_ds_zt, stats )
      call stats_update( "thv_ds_zt", thv_ds_zt, stats )
      call stats_update( "Lscale", Lscale, stats )
      call stats_update( "wpup2", wpup2, stats )
      call stats_update( "wpvp2", wpvp2, stats )
      call stats_update( "wp3", wp3, stats )
      call stats_update( "wpthlp2", wpthlp2, stats )
      call stats_update( "wp2thlp", wp2thlp, stats )
      call stats_update( "wprtp2", wprtp2, stats )
      call stats_update( "wp2rtp", wp2rtp, stats )
      call stats_update( "Lscale_up", Lscale_up, stats )
      call stats_update( "Lscale_down", Lscale_down, stats )
      call stats_update( "tau_zt", tau_zt, stats )
      call stats_update( "Kh_zt", Kh_zt, stats )
      call stats_update( "wp2thvp", wp2thvp, stats )
      call stats_update( "wp2up", wp2up, stats )
      call stats_update( "wp2rcp", wp2rcp, stats )
      call stats_update( "w_up_in_cloud", w_up_in_cloud, stats )
      call stats_update( "w_down_in_cloud", w_down_in_cloud, stats )
      call stats_update( "cld_updr_frac", cloudy_updraft_frac, stats )
      call stats_update( "cld_downdr_frac", cloudy_downdraft_frac, stats )
      call stats_update( "wprtpthlp", wprtpthlp, stats )
      call stats_update( "rc_coef", rc_coef, stats )
      call stats_update( "sigma_sqd_w_zt", sigma_sqd_w_zt, stats )
      call stats_update( "rho", rho, stats )
      call stats_update( "rsat", rsat, stats )
      call stats_update( "mixt_frac", pdf_params%mixt_frac, stats )
      call stats_update( "w_1", pdf_params%w_1, stats )
      call stats_update( "w_2", pdf_params%w_2, stats )
      call stats_update( "varnce_w_1", pdf_params%varnce_w_1, stats )
      call stats_update( "varnce_w_2", pdf_params%varnce_w_2, stats )
      call stats_update( "thl_1", pdf_params%thl_1, stats )
      call stats_update( "thl_2", pdf_params%thl_2, stats )
      call stats_update( "varnce_thl_1", pdf_params%varnce_thl_1, stats )
      call stats_update( "varnce_thl_2", pdf_params%varnce_thl_2, stats )
      call stats_update( "rt_1", pdf_params%rt_1, stats )
      call stats_update( "rt_2", pdf_params%rt_2, stats )
      call stats_update( "varnce_rt_1", pdf_params%varnce_rt_1, stats )
      call stats_update( "varnce_rt_2", pdf_params%varnce_rt_2, stats )
      call stats_update( "rc_1", pdf_params%rc_1, stats )
      call stats_update( "rc_2", pdf_params%rc_2, stats )
      call stats_update( "rsatl_1", pdf_params%rsatl_1, stats )
      call stats_update( "rsatl_2", pdf_params%rsatl_2, stats )
      call stats_update( "cloud_frac_1", pdf_params%cloud_frac_1, stats )
      call stats_update( "cloud_frac_2", pdf_params%cloud_frac_2, stats )
      call stats_update( "chi_1", pdf_params%chi_1, stats )
      call stats_update( "chi_2", pdf_params%chi_2, stats )
      call stats_update( "stdev_chi_1", pdf_params%stdev_chi_1, stats )
      call stats_update( "stdev_chi_2", pdf_params%stdev_chi_2, stats )
      call stats_update( "stdev_eta_1", pdf_params%stdev_eta_1, stats )
      call stats_update( "stdev_eta_2", pdf_params%stdev_eta_2, stats )
      call stats_update( "covar_chi_eta_1", pdf_params%covar_chi_eta_1, stats )
      call stats_update( "covar_chi_eta_2", pdf_params%covar_chi_eta_2, stats )
      call stats_update( "corr_w_chi_1", pdf_params%corr_w_chi_1, stats )
      call stats_update( "corr_w_chi_2", pdf_params%corr_w_chi_2, stats )
      call stats_update( "corr_w_eta_1", pdf_params%corr_w_eta_1, stats )
      call stats_update( "corr_w_eta_2", pdf_params%corr_w_eta_2, stats )
      call stats_update( "corr_chi_eta_1", pdf_params%corr_chi_eta_1, stats )
      call stats_update( "corr_chi_eta_2", pdf_params%corr_chi_eta_2, stats )
      call stats_update( "corr_w_rt_1", pdf_params%corr_w_rt_1, stats )
      call stats_update( "corr_w_rt_2", pdf_params%corr_w_rt_2, stats )
      call stats_update( "corr_w_thl_1", pdf_params%corr_w_thl_1, stats )
      call stats_update( "corr_w_thl_2", pdf_params%corr_w_thl_2, stats )
      call stats_update( "corr_rt_thl_1", pdf_params%corr_rt_thl_1, stats )
      call stats_update( "corr_rt_thl_2", pdf_params%corr_rt_thl_2, stats )
      call stats_update( "crt_1", pdf_params%crt_1, stats )
      call stats_update( "crt_2", pdf_params%crt_2, stats )
      call stats_update( "cthl_1", pdf_params%cthl_1, stats )
      call stats_update( "cthl_2", pdf_params%cthl_2, stats )
      call stats_update( "wp2_zt", wp2_zt, stats )
      call stats_update( "thlp2_zt", thlp2_zt, stats )
      call stats_update( "thlp3", thlp3, stats )
      call stats_update( "wpthlp_zt", wpthlp_zt, stats )
      call stats_update( "wprtp_zt", wprtp_zt, stats )
      call stats_update( "rtp2_zt", rtp2_zt, stats )
      call stats_update( "rtp3", rtp3, stats )
      call stats_update( "rtpthlp_zt", rtpthlp_zt, stats )
      call stats_update( "up2_zt", up2_zt, stats )
      call stats_update( "vp2_zt", vp2_zt, stats )
      call stats_update( "upwp_zt", upwp_zt, stats )
      call stats_update( "vpwp_zt", vpwp_zt, stats )
      call stats_update( "a3_coef_zt", a3_coef_zt, stats )
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
      call stats_update( "wp2", wp2, stats )
      call stats_update( "wp3_zm", wp3_zm, stats )
      call stats_update( "rtp2", rtp2, stats )
      call stats_update( "thlp2", thlp2, stats )
      call stats_update( "rtpthlp", rtpthlp, stats )
      call stats_update( "wprtp", wprtp, stats )
      call stats_update( "wpthlp", wpthlp, stats )
      call stats_update( "wp2up2", wp2up2, stats )
      call stats_update( "wp2vp2", wp2vp2, stats )
      call stats_update( "wp4", wp4, stats )
      call stats_update( "wpthvp", wpthvp, stats )
      call stats_update( "rtpthvp", rtpthvp, stats )
      call stats_update( "thlpthvp", thlpthvp, stats )
      call stats_update( "tau_zm", tau_zm, stats )
      call stats_update( "Kh_zm", Kh_zm, stats )
      call stats_update( "wprcp", wprcp, stats )
      call stats_update( "rc_coef_zm", rc_coef_zm, stats )
      call stats_update( "thlprcp", thlprcp, stats )
      call stats_update( "rtprcp", rtprcp, stats )
      call stats_update( "rcp2", rcp2, stats )
      call stats_update( "upwp", upwp, stats )
      call stats_update( "vpwp", vpwp, stats )
      call stats_update( "vp2", vp2, stats )
      call stats_update( "up2", up2, stats )
      call stats_update( "rho_zm", rho_zm, stats )
      call stats_update( "sigma_sqd_w", sigma_sqd_w, stats )
      call stats_update( "rho_ds_zm", rho_ds_zm, stats )
      call stats_update( "thv_ds_zm", thv_ds_zm, stats )
      call stats_update( "em", em, stats )
      call stats_update( "Skw_velocity", Skw_velocity, stats )
      call stats_update( "a3_coef", a3_coef, stats )
      call stats_update( "wp3_on_wp2", wp3_on_wp2, stats )
      call stats_update( "wp3_on_wp2_cfl_num", wp3_on_wp2 * dt / dzm, stats )
      call stats_update( "cloud_frac_zm", cloud_frac_zm, stats )
      call stats_update( "ice_supersat_frac_zm", ice_supersat_frac_zm, stats )
      call stats_update( "rcm_zm", rcm_zm, stats )
      call stats_update( "rtm_zm", rtm_zm, stats )
      call stats_update( "thlm_zm", thlm_zm, stats )
      if ( l_call_pdf_closure_twice ) then
        call stats_update( "w_1_zm", pdf_params_zm%w_1, stats )
        call stats_update( "w_2_zm", pdf_params_zm%w_2, stats )
        call stats_update( "varnce_w_1_zm", pdf_params_zm%varnce_w_1, stats )
        call stats_update( "varnce_w_2_zm", pdf_params_zm%varnce_w_2, stats )
        call stats_update( "mixt_frac_zm", pdf_params_zm%mixt_frac, stats )
      end if
      if ( sclr_dim > 0 ) then
        do sclr = 1, sclr_dim
          write( sclr_idx, * ) sclr
          sclr_idx = adjustl( sclr_idx )
          call stats_update( "sclr"//trim(sclr_idx)//"p2", sclrp2(:,:,sclr), stats )
          call stats_update( "sclr"//trim(sclr_idx)//"prtp", sclrprtp(:,:,sclr), stats )
          call stats_update( "sclr"//trim(sclr_idx)//"pthvp", sclrpthvp(:,:,sclr), stats )
          call stats_update( "sclr"//trim(sclr_idx)//"pthlp", sclrpthlp(:,:,sclr), stats )
          call stats_update( "sclr"//trim(sclr_idx)//"prcp", sclrprcp(:,:,sclr), stats )
          call stats_update( "wpsclr"//trim(sclr_idx)//"p", wpsclrp(:,:,sclr), stats )
          call stats_update( "wp2sclr"//trim(sclr_idx)//"p", wp2sclrp(:,:,sclr), stats )
          call stats_update( "wpsclr"//trim(sclr_idx)//"p2", wpsclrp2(:,:,sclr), stats )
          call stats_update( "wpsclr"//trim(sclr_idx)//"prtp", wpsclrprtp(:,:,sclr), stats )
          call stats_update( "wpsclr"//trim(sclr_idx)//"pthlp", wpsclrpthlp(:,:,sclr), stats )
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
            xtmp(i) = zt(i,1)
          elseif ( k > 1 .and. k < nzt ) then
            xtmp(i) = lin_interpolate_two_points( rc_tol, rcm(i,k), rcm(i,k-1), zt(i,k), zt(i,k-1) )
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
        stats_update, &
        var_on_stats_list

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
