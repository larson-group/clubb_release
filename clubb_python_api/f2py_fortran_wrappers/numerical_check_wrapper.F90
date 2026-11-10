! numerical_check_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module numerical_check

subroutine f2py_check_clubb_settings(ngrdcol, edsclr_dim, params, &
    l_implemented, l_input_fields)

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use derived_type_storage, only: stored_sclr_idx, stored_config_flags, stored_err_info
  use numerical_check, only: check_clubb_settings_api

  implicit none

  integer, intent(in) :: ngrdcol, edsclr_dim
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: params
  logical, intent(in) :: l_implemented, l_input_fields

  call check_clubb_settings_api(ngrdcol, edsclr_dim, stored_sclr_idx, params, &
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

subroutine f2py_sfc_varnce_check( ngrdcol, sclr_dim, sclr_dim_transport, wp2_sfc, & ! In
                                  up2_sfc, vp2_sfc, thlp2_sfc, rtp2_sfc, rtpthlp_sfc, & ! In
                                  sclrp2_sfc, sclrprtp_sfc, sclrpthlp_sfc ) ! In

  ! Description:
  !   Call sfc_varnce_check with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_err_info
  use numerical_check, only: sfc_varnce_check

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  integer, intent(in) :: sclr_dim ! Number of passive scalars [-]
  integer, intent(in) :: sclr_dim_transport ! Passive-scalar storage extent [-]
  ! Vertical velocity variance        [m^2/s^2]
  real(core_rknd), dimension(ngrdcol), intent(in) :: wp2_sfc
  ! u'^2                              [m^2/s^2]
  real(core_rknd), dimension(ngrdcol), intent(in) :: up2_sfc
  ! u'^2                              [m^2/s^2]
  real(core_rknd), dimension(ngrdcol), intent(in) :: vp2_sfc
  ! thetal variance                   [K^2]
  real(core_rknd), dimension(ngrdcol), intent(in) :: thlp2_sfc
  ! rt variance                       [(kg/kg)^2]
  real(core_rknd), dimension(ngrdcol), intent(in) :: rtp2_sfc
  ! thetal rt covariance              [kg K/kg]
  real(core_rknd), dimension(ngrdcol), intent(in) :: rtpthlp_sfc
  ! Passive scalar variance                 [units^2]
  real(core_rknd), dimension(ngrdcol, sclr_dim_transport), intent(in) :: sclrp2_sfc
  ! Passive scalar r_t covariance           [units kg/kg]
  real(core_rknd), dimension(ngrdcol, sclr_dim_transport), intent(in) :: sclrprtp_sfc
  ! Passive scalar theta_l covariance       [units K]
  real(core_rknd), dimension(ngrdcol, sclr_dim_transport), intent(in) :: sclrpthlp_sfc

  !--------------------------- Begin Code ---------------------------

  call sfc_varnce_check( ngrdcol, sclr_dim, wp2_sfc, & ! In
                         up2_sfc, vp2_sfc, & ! In
                         thlp2_sfc, rtp2_sfc, & ! In
                         rtpthlp_sfc, & ! In
                         sclrp2_sfc, & ! In
                         sclrprtp_sfc, & ! In
                         sclrpthlp_sfc, & ! In
                         stored_err_info ) ! InOut

end subroutine f2py_sfc_varnce_check

subroutine f2py_parameterization_check( &
    nzm, nzt, ngrdcol, sclr_dim, edsclr_dim, sclr_dim_transport, edsclr_dim_transport, &
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

  integer, intent(in) :: nzm, nzt, ngrdcol, sclr_dim, edsclr_dim, sclr_dim_transport, edsclr_dim_transport
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    thlm_forcing, rtm_forcing, um_forcing, vm_forcing, wm_zt, p_in_Pa, rho, exner, &
    rho_ds_zt, invrs_rho_ds_zt, thv_ds_zt, um, vm, rtm, thlm, wp3
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: &
    wm_zm, rho_zm, rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, upwp, vpwp, up2, vp2, &
    wprtp, wpthlp, wp2, rtp2, thlp2, rtpthlp
  real(core_rknd), dimension(ngrdcol), intent(in) :: &
    wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc
  character(len=*), intent(in) :: prefix
  real(core_rknd), dimension(ngrdcol, sclr_dim_transport), intent(in) :: wpsclrp_sfc
  real(core_rknd), dimension(ngrdcol, edsclr_dim_transport), intent(in) :: wpedsclrp_sfc
  real(core_rknd), dimension(ngrdcol, nzt, sclr_dim_transport), intent(in) :: sclrm, sclrm_forcing
  real(core_rknd), dimension(ngrdcol, nzm, sclr_dim_transport), intent(in) :: wpsclrp, sclrp2, sclrprtp, sclrpthlp
  real(core_rknd), dimension(ngrdcol, nzt, edsclr_dim_transport), intent(in) :: edsclrm, edsclrm_forcing

  call parameterization_check( &
    nzm, nzt, ngrdcol, sclr_dim, edsclr_dim, &
    thlm_forcing, rtm_forcing, um_forcing, vm_forcing, wm_zm, wm_zt, p_in_Pa, &
    rho_zm, rho, exner, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, &
    thv_ds_zm, thv_ds_zt, wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, &
    um, upwp, vm, vpwp, up2, vp2, rtm, wprtp, thlm, wpthlp, wp2, wp3, rtp2, thlp2, &
    rtpthlp, prefix, wpsclrp_sfc, wpedsclrp_sfc, sclrm, wpsclrp, sclrp2, sclrprtp, &
    sclrpthlp, sclrm_forcing, edsclrm, edsclrm_forcing, stored_err_info )

end subroutine f2py_parameterization_check

subroutine f2py_pdf_closure_check( nz, ngrdcol, sclr_dim, sclr_dim_transport, wp4, & ! In
                                   wprtp2, wp2rtp, wpthlp2, wp2thlp, cloud_frac, rcm, & ! In
                                   wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, wprcp, & ! In
                                   wp2rcp, rtprcp, thlprcp, rcp2, wprtpthlp, crt_1, & ! In
                                   crt_2, cthl_1, cthl_2, sclrpthvp, sclrprcp, & ! In
                                   wpsclrp2, wpsclrprtp, wpsclrpthlp, wp2sclrp ) ! In

  ! Description:
  !   Call pdf_closure_check with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_pdf_params, stored_stats, stored_err_info
  use numerical_check, only: pdf_closure_check

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: nz ! Number of vertical levels [-]
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  integer, intent(in) :: sclr_dim ! Number of active passive scalars [-]
  integer, intent(in) :: sclr_dim_transport ! Passive-scalar storage extent [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wp4 ! w'^4                  [m^4/s^4]
  ! w' r_t'               [(m kg)/(s kg)]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wprtp2
  ! w'^2 r_t'             [(m^2 kg)/(s^2 kg)]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wp2rtp
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wpthlp2 ! w' th_l'^2            [(m K^2)/s]
  ! w'^2 th_l'            [(m^2 K)/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wp2thlp
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: cloud_frac ! Cloud fraction        [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: rcm ! Mean liquid water     [kg/kg]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wpthvp ! Buoyancy flux         [(K m)/s]
  ! w'^2 th_v'            [(m^2 K)/s^2]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wp2thvp
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wp2up ! w'^2 u'               [m^3/s^3]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: rtpthvp ! r_t' th_v'            [(kg K)/kg]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: thlpthvp ! th_l' th_v'           [K^2]
  ! w' r_c'               [(m kg)/(s kg)]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wprcp
  ! w'^2 r_c'             [(m^2 kg)/(s^2 kg)]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wp2rcp
  ! r_t' r_c'             [(kg^2)/(kg^2)]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: rtprcp
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: thlprcp ! th_l' r_c'            [(K kg)/kg]
  ! r_c'^2                [(kg^2)/(kg^2)]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: rcp2
  ! w' r_t' th_l'         [(m kg K)/(s kg)]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wprtpthlp
  ! r_t coef. in chi/eta eqns. (1st PDF comp.)        [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: crt_1
  ! r_t coef. in chi/eta eqns. (2nd PDF comp.)        [-]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: crt_2
  ! th_l coef.: chi/eta eqns. (1st PDF comp.) [(kg/kg)/K]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: cthl_1
  ! th_l coef.: chi/eta eqns. (2nd PDF comp.) [(kg/kg)/K]
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: cthl_2
  ! < sclr' th_v' > (momentum levels)   [units vary]
  real(core_rknd), dimension(ngrdcol, nz, sclr_dim_transport), intent(in) :: sclrpthvp
  ! < sclr' r_c' > on momentum grid
  real(core_rknd), dimension(ngrdcol, nz, sclr_dim_transport), intent(in) :: sclrprcp
  real(core_rknd), dimension(ngrdcol, nz, sclr_dim_transport), intent(in) :: wpsclrp2 ! w'sclr'^2
  real(core_rknd), dimension(ngrdcol, nz, sclr_dim_transport), intent(in) :: wpsclrprtp ! w'sclr'rt'
  ! w'sclr'thl'
  real(core_rknd), dimension(ngrdcol, nz, sclr_dim_transport), intent(in) :: wpsclrpthlp
  real(core_rknd), dimension(ngrdcol, nz, sclr_dim_transport), intent(in) :: wp2sclrp ! w'^2 sclr'

  !--------------------------- Begin Code ---------------------------

  call pdf_closure_check( nz, ngrdcol, sclr_dim, wp4, wprtp2, wp2rtp, wpthlp2, wp2thlp, & ! In
                          cloud_frac, rcm, wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, & ! In
                          wprcp, & ! In
                          wp2rcp, rtprcp, thlprcp, rcp2, wprtpthlp, & ! In
                          crt_1, crt_2, cthl_1, & ! In
                          cthl_2, & ! In
                          stored_pdf_params, sclrpthvp, sclrprcp, wpsclrp2, wpsclrprtp, & ! In
                          wpsclrpthlp, & ! In
                          wp2sclrp, stored_stats, & ! In
                          stored_err_info ) ! InOut

end subroutine f2py_pdf_closure_check

subroutine f2py_length_check( nzt, ngrdcol, lscale, lscale_up, lscale_down ) ! In

  ! Description:
  !   Call length_check with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_err_info
  use numerical_check, only: length_check

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: nzt ! Number of thermodynamic levels [-]
  integer, intent(in) :: ngrdcol ! Number of grid columns [-]
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: lscale ! Mixing length                 [m]
  ! Upward mixing length          [m]
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: lscale_up
  ! Downward mixing length        [m]
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: lscale_down

  !--------------------------- Begin Code ---------------------------

  call length_check(nzt, ngrdcol, lscale, lscale_up, lscale_down, & ! In
    stored_err_info) ! InOut

end subroutine f2py_length_check
