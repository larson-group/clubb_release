! numerical_check_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module numerical_check

subroutine f2py_check_clubb_settings(ngrdcol, params, &
    l_implemented, l_input_fields)

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use derived_type_storage, only: stored_config_flags, stored_err_info
  use numerical_check, only: check_clubb_settings_api

  implicit none

  integer, intent(in) :: ngrdcol
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: params
  logical, intent(in) :: l_implemented, l_input_fields

  call check_clubb_settings_api(ngrdcol, params, &
    l_implemented, l_input_fields, &
    stored_config_flags, stored_err_info)

end subroutine f2py_check_clubb_settings

subroutine f2py_calculate_spurious_source(integral_after, integral_before, &
    flux_top, flux_sfc, integral_forcing, dt, spurious_source)

  use clubb_precision, only: core_rknd
  use numerical_check, only: calculate_spurious_source

  implicit none

  real(core_rknd), intent(in) :: integral_after, integral_before
  real(core_rknd), intent(in) :: flux_top, flux_sfc, integral_forcing, dt
  real(core_rknd), intent(out) :: spurious_source

  spurious_source = calculate_spurious_source( &
    integral_after, integral_before, flux_top, flux_sfc, integral_forcing, dt)

end subroutine f2py_calculate_spurious_source

subroutine f2py_sfc_varnce_check(sclr_dim, sclr_dim_transport, wp2_sfc, up2_sfc, vp2_sfc, thlp2_sfc, &
    rtp2_sfc, rtpthlp_sfc, sclrp2_sfc, sclrprtp_sfc, sclrpthlp_sfc)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_err_info
  use numerical_check, only: sfc_varnce_check

  implicit none

  integer, intent(in) :: sclr_dim, sclr_dim_transport
  real(core_rknd), intent(in) :: wp2_sfc, up2_sfc, vp2_sfc, thlp2_sfc, rtp2_sfc, rtpthlp_sfc
  real(core_rknd), dimension(sclr_dim_transport), intent(in) :: sclrp2_sfc, sclrprtp_sfc, sclrpthlp_sfc

  call sfc_varnce_check(sclr_dim, wp2_sfc, up2_sfc, vp2_sfc, thlp2_sfc, rtp2_sfc, rtpthlp_sfc, &
    sclrp2_sfc, sclrprtp_sfc, sclrpthlp_sfc, stored_err_info)

end subroutine f2py_sfc_varnce_check

subroutine f2py_parameterization_check( &
    nzm, nzt, sclr_dim, edsclr_dim, sclr_dim_transport, edsclr_dim_transport, &
    thlm_forcing, rtm_forcing, um_forcing, vm_forcing, wm_zm, wm_zt, p_in_Pa, &
    rho_zm, rho, exner, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, &
    thv_ds_zm, thv_ds_zt, wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, &
    um, upwp, vm, vpwp, up2, vp2, rtm, wprtp, thlm, wpthlp, wp2, wp3, &
    rtp2, thlp2, rtpthlp, prefix, wpsclrp_sfc, wpedsclrp_sfc, sclrm, wpsclrp, &
    sclrp2, sclrprtp, sclrpthlp, sclrm_forcing, edsclrm, edsclrm_forcing)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_err_info
  use numerical_check, only: parameterization_check

  implicit none

  integer, intent(in) :: nzm, nzt, sclr_dim, edsclr_dim, sclr_dim_transport, edsclr_dim_transport
  real(core_rknd), dimension(nzt), intent(in) :: &
    thlm_forcing, rtm_forcing, um_forcing, vm_forcing, wm_zt, p_in_Pa, rho, exner, &
    rho_ds_zt, invrs_rho_ds_zt, thv_ds_zt, um, vm, rtm, thlm, wp3
  real(core_rknd), dimension(nzm), intent(in) :: &
    wm_zm, rho_zm, rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, upwp, vpwp, up2, vp2, &
    wprtp, wpthlp, wp2, rtp2, thlp2, rtpthlp
  real(core_rknd), intent(in) :: &
    wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc
  character(len=*), intent(in) :: prefix
  real(core_rknd), dimension(sclr_dim_transport), intent(in) :: wpsclrp_sfc
  real(core_rknd), dimension(edsclr_dim_transport), intent(in) :: wpedsclrp_sfc
  real(core_rknd), dimension(nzt, sclr_dim_transport), intent(in) :: sclrm, sclrm_forcing
  real(core_rknd), dimension(nzm, sclr_dim_transport), intent(in) :: wpsclrp, sclrp2, sclrprtp, sclrpthlp
  real(core_rknd), dimension(nzt, edsclr_dim_transport), intent(in) :: edsclrm, edsclrm_forcing

  call parameterization_check( &
    nzm, nzt, sclr_dim, edsclr_dim, &
    thlm_forcing, rtm_forcing, um_forcing, vm_forcing, wm_zm, wm_zt, p_in_Pa, &
    rho_zm, rho, exner, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, &
    thv_ds_zm, thv_ds_zt, wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, &
    um, upwp, vm, vpwp, up2, vp2, rtm, wprtp, thlm, wpthlp, wp2, wp3, rtp2, thlp2, &
    rtpthlp, prefix, wpsclrp_sfc, wpedsclrp_sfc, sclrm, wpsclrp, sclrp2, sclrprtp, &
    sclrpthlp, sclrm_forcing, edsclrm, edsclrm_forcing, stored_err_info )

end subroutine f2py_parameterization_check


subroutine f2py_pdf_closure_check(nz, sclr_dim, sclr_dim_transport, &
    wp4, wprtp2, wp2rtp, wpthlp2, wp2thlp, cloud_frac, rcm, wpthvp, wp2thvp, wp2up, &
    rtpthvp, thlpthvp, wprcp, wp2rcp, rtprcp, thlprcp, rcp2, wprtpthlp, crt_1, crt_2, &
    cthl_1, cthl_2, sclrpthvp, sclrprcp, wpsclrp2, wpsclrprtp, wpsclrpthlp, wp2sclrp)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_pdf_params, stored_stats, stored_err_info
  use numerical_check, only: pdf_closure_check

  implicit none

  integer, intent(in) :: nz, sclr_dim, sclr_dim_transport
  real(core_rknd), dimension(nz), intent(in) :: &
    wp4, wprtp2, wp2rtp, wpthlp2, wp2thlp, cloud_frac, rcm, wpthvp, wp2thvp, wp2up, &
    rtpthvp, thlpthvp, wprcp, wp2rcp, rtprcp, thlprcp, rcp2, wprtpthlp, crt_1, crt_2, &
    cthl_1, cthl_2
  real(core_rknd), dimension(nz, sclr_dim_transport), intent(in) :: &
    sclrpthvp, sclrprcp, wpsclrp2, wpsclrprtp, wpsclrpthlp, wp2sclrp

  call pdf_closure_check(nz, sclr_dim, wp4, wprtp2, wp2rtp, wpthlp2, wp2thlp, cloud_frac, rcm, &
    wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, wprcp, wp2rcp, rtprcp, thlprcp, rcp2, &
    wprtpthlp, crt_1, crt_2, cthl_1, cthl_2, stored_pdf_params, sclrpthvp, sclrprcp, wpsclrp2, &
    wpsclrprtp, wpsclrpthlp, wp2sclrp, stored_stats, stored_err_info)

end subroutine f2py_pdf_closure_check

subroutine f2py_length_check(nzt, lscale, lscale_up, lscale_down)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_err_info
  use numerical_check, only: length_check

  implicit none

  integer, intent(in) :: nzt
  real(core_rknd), dimension(nzt), intent(in) :: lscale, lscale_up, lscale_down

  call length_check(nzt, lscale, lscale_up, lscale_down, stored_err_info)

end subroutine f2py_length_check
