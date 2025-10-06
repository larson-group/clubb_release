!------------------------------------------------------------------------
! $Id$
!===============================================================================
module numerical_check

  implicit none

!       Made is_nan_2d public so it may be used
!       for finding code that cause NaNs
!       Joshua Fasching November 2007

!       *_check subroutines were added to ensure that the
!       subroutines they are checking perform correctly
!       Joshua Fasching February 2008

!       rad_clipping has been replaced by rad_check as the new
!       subroutine only reports if there are invalid values.
!       Joshua Fasching March 2008

  private ! Default scope

  public :: invalid_model_arrays, is_nan_2d,  &
            rad_check, parameterization_check, &
            sfc_varnce_check, pdf_closure_check, &
            length_check, is_nan_sclr, calculate_spurious_source, &
            check_clubb_settings_api

  private :: check_negative, check_nan


  ! Abstraction of check_nan
  interface check_nan
    module procedure check_nan_sclr, check_nan_2d
  end interface

  ! Abstraction of check_negative
  interface check_negative
    module procedure check_negative_index!, check_negative_total
  end interface


  contains
!---------------------------------------------------------------------------------
  subroutine length_check( nzt, Lscale, Lscale_up, Lscale_down, err_info )
!
!        Description: This subroutine determines if any of the output
!        variables for the length_new subroutine carry values that
!        are NaNs.
!
!        Joshua Fasching February 2008
!---------------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use err_info_type_module, only: &
      err_info_type     ! Type

    implicit none

    integer, intent(in) :: &
      nzt

    ! Constant Parameters
    character(*), parameter :: proc_name = "compute_mixing_length"

    ! Input Variables
    real( kind = core_rknd ), dimension(nzt), intent(in) :: &
      Lscale,     & ! Mixing length                 [m]
      Lscale_up,  & ! Upward mixing length          [m]
      Lscale_down   ! Downward mixing length        [m]

    ! Input/Output Variables
    type(err_info_type), intent(inout) :: &
      err_info      ! err_info struct containing err_code and err_header
!-----------------------------------------------------------------------------

    call check_nan( Lscale, "Lscale", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( Lscale_up, "Lscale_up", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( Lscale_down, "Lscale_down", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)

    return
  end subroutine length_check

!---------------------------------------------------------------------------
  subroutine pdf_closure_check( nz, sclr_dim, &
                                wp4, wprtp2, wp2rtp, wpthlp2, &
                                wp2thlp, cloud_frac, rcm, wpthvp, wp2thvp, &
                                rtpthvp, thlpthvp, wprcp, wp2rcp, & 
                                rtprcp, thlprcp, rcp2, wprtpthlp, & 
                                crt_1, crt_2, cthl_1, cthl_2, pdf_params, &
                                sclrpthvp, sclrprcp, wpsclrp2, & 
                                wpsclrprtp, wpsclrpthlp, wp2sclrp, &
                                stats_metadata, &
                                err_info )

! Description: This subroutine determines if any of the output
!   variables for the pdf_closure subroutine carry values that
!   are NaNs.
!
! Joshua Fasching February 2008
!---------------------------------------------------------------------------

    use pdf_parameter_module, only: &
        pdf_parameter  ! type

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use stats_variables, only: &
        stats_metadata_type

    use err_info_type_module, only: &
      err_info_type     ! Type

    implicit none

    integer, intent(in) :: &
      nz, &
      sclr_dim

    ! Parameter Constants
    character(len=*), parameter :: proc_name = &
      "pdf_closure"

    ! Input Variables
    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      wp4,             & ! w'^4                  [m^4/s^4]
      wprtp2,          & ! w' r_t'               [(m kg)/(s kg)]
      wp2rtp,          & ! w'^2 r_t'             [(m^2 kg)/(s^2 kg)]
      wpthlp2,         & ! w' th_l'^2            [(m K^2)/s]
      wp2thlp,         & ! w'^2 th_l'            [(m^2 K)/s^2]
      cloud_frac,      & ! Cloud fraction        [-]
      rcm,             & ! Mean liquid water     [kg/kg]
      wpthvp,          & ! Buoyancy flux         [(K m)/s] 
      wp2thvp,         & ! w'^2 th_v'            [(m^2 K)/s^2]
      rtpthvp,         & ! r_t' th_v'            [(kg K)/kg]
      thlpthvp,        & ! th_l' th_v'           [K^2]
      wprcp,           & ! w' r_c'               [(m kg)/(s kg)]
      wp2rcp,          & ! w'^2 r_c'             [(m^2 kg)/(s^2 kg)]
      rtprcp,          & ! r_t' r_c'             [(kg^2)/(kg^2)]
      thlprcp,         & ! th_l' r_c'            [(K kg)/kg]
      rcp2,            & ! r_c'^2                [(kg^2)/(kg^2)]
      wprtpthlp,       & ! w' r_t' th_l'         [(m kg K)/(s kg)]
      crt_1, crt_2,  &
      cthl_1, cthl_2

    type(pdf_parameter), intent(in) :: &
      pdf_params        ! PDF parameters          [units vary]

    ! Input (Optional passive scalar variables)
    real( kind = core_rknd ), dimension(nz,sclr_dim), intent(in) :: &
      sclrpthvp, &
      sclrprcp, &
      wpsclrp2, &
      wpsclrprtp, &
      wpsclrpthlp, &
      wp2sclrp

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    ! Input/Output variables
    type(err_info_type), intent(inout) :: &
      err_info      ! err_info struct containing err_code and err_header

    ! Local variables
    integer:: &
      sclr                ! Scalar loop index

!-------------------------------------------------------------------------------

    ! ---- Begin Code ----

    if ( stats_metadata%iwp4 > 0 ) call check_nan( wp4,"wp4", proc_name, & ! intent(in)
                                                   err_info ) ! intent(inout)
    if ( stats_metadata%iwprtp2 > 0 ) call check_nan( wprtp2,"wprtp2", proc_name, & ! intent(in)
                                                      err_info ) ! intent(inout)
    call check_nan( wp2rtp,"wp2rtp", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    if ( stats_metadata%iwpthlp2 > 0 ) call check_nan( wpthlp2,"wpthlp2", proc_name, & ! intent(in)
                                                       err_info ) ! intent(inout)
    call check_nan( wp2thlp,"wp2thlp", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( cloud_frac,"cloud_frac", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( rcm,"rcm", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( wpthvp, "wpthvp", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( wp2thvp, "wp2thvp", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( rtpthvp, "rtpthvp", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( thlpthvp, "thlpthvp", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( wprcp, "wprcp", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( wp2rcp, "wp2rcp", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( rtprcp, "rtprcp", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( thlprcp, "thlprcp", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    if ( stats_metadata%ircp2 >  0 ) call check_nan( rcp2, "rcp2", proc_name, & ! intent(in)
                                                     err_info ) ! intent(inout)
    if ( stats_metadata%iwprtpthlp > 0 ) &
         call check_nan( wprtpthlp, "wprtpthlp", proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_nan( crt_1, "crt_1", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( crt_2, "crt_2", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( cthl_1, "cthl_1", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( cthl_2, "cthl_2", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    ! Check each PDF parameter at the grid level sent in.
    call check_nan( pdf_params%w_1(1,:), "pdf_params%w_1(1,:)", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%w_2(1,:), "pdf_params%w_2(1,:)", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%varnce_w_1(1,:), "pdf_params%varnce_w_1(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%varnce_w_2(1,:), "pdf_params%varnce_w_2(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%rt_1(1,:), "pdf_params%rt_1(1,:)", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%rt_2(1,:), "pdf_params%rt_2(1,:)", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%varnce_rt_1(1,:), "pdf_params%varnce_rt_1(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%varnce_rt_2(1,:), "pdf_params%varnce_rt_2(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%thl_1(1,:), "pdf_params%thl_1(1,:)", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%thl_2(1,:), "pdf_params%thl_2(1,:)", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%varnce_thl_1(1,:), "pdf_params%varnce_thl_1(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%varnce_thl_2(1,:), "pdf_params%varnce_thl_2(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%mixt_frac(1,:), "pdf_params%mixt_frac(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%corr_w_rt_1(1,:), "pdf_params%corr_w_rt_1(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%corr_w_rt_2(1,:), "pdf_params%corr_w_rt_2(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%corr_w_thl_1(1,:), "pdf_params%corr_w_thl_1(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%corr_w_thl_2(1,:), "pdf_params%corr_w_thl_2(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%corr_rt_thl_1(1,:), "pdf_params%corr_rt_thl_1(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%corr_rt_thl_2(1,:), "pdf_params%corr_rt_thl_2(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%rc_1(1,:), "pdf_params%rc_1(1,:)", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%rc_2(1,:), "pdf_params%rc_2(1,:)", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%rsatl_1(1,:), "pdf_params%rsatl_1(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%rsatl_2(1,:), "pdf_params%rsatl_2(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%cloud_frac_1(1,:), "pdf_params%cloud_frac_1(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%cloud_frac_2(1,:), "pdf_params%cloud_frac_2(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%chi_1(1,:), "pdf_params%chi_1(1,:)", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%chi_2(1,:), "pdf_params%chi_2(1,:)", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%stdev_chi_1(1,:), "pdf_params%stdev_chi_1(1,:)", &! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%stdev_chi_2(1,:), "pdf_params%stdev_chi_2(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%stdev_eta_1(1,:), "pdf_params%stdev_eta_1(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%stdev_eta_2(1,:), "pdf_params%stdev_eta_2(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%covar_chi_eta_1(1,:), "pdf_params%covar_chi_eta_1(1,:)",&!intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%covar_chi_eta_2(1,:), "pdf_params%covar_chi_eta_2(1,:)",&!intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%corr_w_chi_1(1,:), "pdf_params%corr_w_chi_1(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%corr_w_chi_2(1,:), "pdf_params%corr_w_chi_2(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%corr_w_eta_1(1,:), "pdf_params%corr_w_eta_1(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%corr_w_eta_2(1,:), "pdf_params%corr_w_eta_2(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%corr_chi_eta_1(1,:), "pdf_params%corr_chi_eta_1(1,:)", & !intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%corr_chi_eta_2(1,:), "pdf_params%corr_chi_eta_2(1,:)", & !intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%alpha_thl(1,:), "pdf_params%alpha_thl(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%alpha_rt(1,:), "pdf_params%alpha_rt(1,:)", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%ice_supersat_frac_1(1,:), & ! intent(in)
                    "pdf_params%ice_supersat_frac_1(1,:)", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( pdf_params%ice_supersat_frac_2(1,:), & ! intent(in)
                    "pdf_params%ice_supersat_frac_2(1,:)", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)

    if ( sclr_dim > 0 ) then
       do sclr = 1, sclr_dim, 1
          call check_nan( sclrpthvp(:,sclr),"sclrpthvp", & ! intent(in)
                          proc_name, & ! intent(in)
                          err_info ) ! intent(inout)
          call check_nan( sclrprcp(:,sclr), "sclrprcp", & ! intent(in)
                          proc_name, & ! intent(in)
                          err_info ) ! intent(inout)
          call check_nan( wpsclrprtp(:,sclr), "wpsclrprtp", & ! intent(in)
                          proc_name, & ! intent(in)
                          err_info ) ! intent(inout)
          call check_nan( wpsclrp2(:,sclr), "wpsclrp2", & ! intent(in)
                          proc_name, & ! intent(in)
                          err_info ) ! intent(inout)
          call check_nan( wpsclrpthlp(:,sclr), "wpsclrtlp", & ! intent(in)
                          proc_name, & ! intent(in)
                          err_info ) ! intent(inout)
          call check_nan( wp2sclrp(:,sclr), "wp2sclrp", & ! intent(in)
                          proc_name, & ! intent(in)
                          err_info ) ! intent(inout)
       enddo ! sclr = 1, sclr_dim, 1
    endif

    return
  end subroutine pdf_closure_check

!-------------------------------------------------------------------------------
  subroutine parameterization_check(                                        & 
               nzm, nzt, sclr_dim, edsclr_dim,                              & ! intent(in)
               thlm_forcing, rtm_forcing, um_forcing,                       & ! intent(in)
               vm_forcing, wm_zm, wm_zt, p_in_Pa,                           & ! intent(in)
               rho_zm, rho, exner, rho_ds_zm,                               & ! intent(in)
               rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,                 & ! intent(in)
               thv_ds_zm, thv_ds_zt, wpthlp_sfc, wprtp_sfc, upwp_sfc,       & ! intent(in)
               vpwp_sfc, p_sfc, um, upwp, vm, vpwp, up2, vp2,               & ! intent(in)
               rtm, wprtp, thlm, wpthlp, wp2, wp3,                          & ! intent(in)
               rtp2, thlp2, rtpthlp,                                        & ! intent(in)
!              rcm,                                                         &
               prefix,                                                      & ! intent(in)
               wpsclrp_sfc, wpedsclrp_sfc, sclrm, wpsclrp,                  & ! intent(in)
               sclrp2,                                                      & ! intent(in)
               sclrprtp, sclrpthlp, sclrm_forcing, edsclrm,                 & ! intent(in)
               edsclrm_forcing,                                             & ! intent(in)
               err_info )                                                     ! intent(inout)

!
! Description:
!   This subroutine determines what input variables may have NaN values.
!   In addition it checks to see if rho_zm, rho, exner, up2, vp2, rtm, thlm,
!   wp2, rtp2, thlp2, or tau_zm have negative values.
!-------------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level_api,   & ! Procedure
        clubb_no_error,               & ! Constants
        clubb_fatal_error

    use constants_clubb, only: &
        fstderr ! Variable

    use err_info_type_module, only: &
      err_info_type     ! Type

    implicit none

    integer, intent(in) :: &
      nzm, &
      nzt, &
      sclr_dim, &
      edsclr_dim

    ! Constant Parameters
    ! Name of the procedure using parameterization_check
    character(len=18), parameter :: &
      proc_name = "advance_clubb_core"

    ! Input variables
    real( kind = core_rknd ), intent(in), dimension(nzt) :: &
      thlm_forcing,    & ! theta_l forcing (thermodynamic levels)    [K/s]
      rtm_forcing,     & ! r_t forcing (thermodynamic levels)        [(kg/kg)/s]
      um_forcing,      & ! u wind forcing (thermodynamic levels)     [m/s/s]
      vm_forcing,      & ! v wind forcing (thermodynamic levels)     [m/s/s]
      wm_zt,           & ! w mean wind component on thermo. levels   [m/s]
      p_in_Pa,         & ! Air pressure (thermodynamic levels)       [Pa]
      rho,             & ! Air density on thermodynamic levels       [kg/m^3]
      exner,           & ! Exner function (thermodynamic levels)     [-]
      rho_ds_zt,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      thv_ds_zt!,      & ! Dry, base-state theta_v on thermo. levs.  [K]
!     rcm                ! Cloud water mixing ratio  [kg/kg] - Unused

    real( kind = core_rknd ), intent(in), dimension(nzm) :: &
      wm_zm,           & ! w mean wind component on momentum levels  [m/s]
      rho_zm,          & ! Air density on momentum levels            [kg/m^3]
      rho_ds_zm,       & ! Dry, static density on momentum levels    [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      thv_ds_zm          ! Dry, base-state theta_v on momentum levs. [K]

    real( kind = core_rknd ), intent(in) :: &
      wpthlp_sfc,   & ! w' theta_l' at surface.   [(m K)/s]
      wprtp_sfc,    & ! w' r_t' at surface.       [(kg m)/( kg s)]
      upwp_sfc,     & ! u'w' at surface.          [m^2/s^2]
      vpwp_sfc,     & ! v'w' at surface.          [m^2/s^2]
      p_sfc           ! pressure at surface       [Pa]

    ! These are prognostic or are planned to be in the future
    real( kind = core_rknd ), intent(in), dimension(nzt) :: &
      um,      & ! u mean wind component (thermodynamic levels)   [m/s]
      vm,      & ! v mean wind component (thermodynamic levels)   [m/s]
      rtm,     & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      thlm,    & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      wp3        ! w'^3 (thermodynamic levels)                    [m^3/s^3]

    real( kind = core_rknd ), intent(in), dimension(nzm) :: &
      upwp,    & ! u'w' (momentum levels)                         [m^2/s^2]
      vpwp,    & ! v'w' (momentum levels)                         [m^2/s^2]
      up2,     & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2,     & ! v'^2 (momentum levels)                         [m^2/s^2]
      wprtp,   & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
      wpthlp,  & ! w' th_l' (momentum levels)                     [(m/s) K]
      rtp2,    & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      thlp2,   & ! th_l'^2 (momentum levels)                      [K^2]
      rtpthlp, & ! r_t' th_l' (momentum levels)                   [(kg/kg) K]
      wp2        ! w'^2 (momentum levels)                         [m^2/s^2]

    character(len=*), intent(in) :: prefix ! Location where subroutine is called

    real( kind = core_rknd ), intent(in), dimension(sclr_dim) :: &
      wpsclrp_sfc    ! Scalar flux at surface [units m/s]

    real( kind = core_rknd ), intent(in), dimension(edsclr_dim) :: &
      wpedsclrp_sfc ! Eddy-Scalar flux at surface      [units m/s]

    real( kind = core_rknd ), intent(in),dimension(nzt,sclr_dim) :: &
      sclrm,         & ! Passive scalar mean      [units vary]
      sclrm_forcing    ! Passive scalar forcing   [units / s]

    real( kind = core_rknd ), intent(in),dimension(nzm,sclr_dim) :: &
      wpsclrp,       & ! w'sclr'                  [units vary]
      sclrp2,        & ! sclr'^2                  [units vary]
      sclrprtp,      & ! sclr'rt'                 [units vary]
      sclrpthlp        ! sclr'thl'                [units vary]

    real( kind = core_rknd ), intent(in),dimension(nzt,edsclr_dim) :: &
      edsclrm,         & ! Eddy passive scalar mean    [units vary]
      edsclrm_forcing    ! Eddy passive scalar forcing [units / s]

    ! Input/Output variables
    type(err_info_type), intent(inout) :: &
      err_info      ! err_info struct containing err_code and err_header

    ! Local Variables
    integer :: sclr, edsclr ! Loop iterator for the scalars
    integer :: k ! Vertical grid level

!-------- Input Nan Check ----------------------------------------------

    call check_nan( thlm_forcing, "thlm_forcing", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( rtm_forcing,"rtm_forcing", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( um_forcing,"um_forcing", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( vm_forcing,"vm_forcing", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)

    call check_nan( wm_zm, "wm_zm", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( wm_zt, "wm_zt", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( p_in_Pa, "p_in_Pa", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( rho_zm, "rho_zm", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( rho, "rho", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( exner, "exner", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( rho_ds_zm, "rho_ds_zm", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( rho_ds_zt, "rho_ds_zt", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( invrs_rho_ds_zm, "invrs_rho_ds_zm", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( invrs_rho_ds_zt, "invrs_rho_ds_zt", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( thv_ds_zm, "thv_ds_zm", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( thv_ds_zt, "thv_ds_zt", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)

    call check_nan( um, "um", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( upwp, "upwp", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( vm, "vm", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( vpwp, "vpwp", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( up2, "up2", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( vp2, "vp2", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( rtm, "rtm", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( wprtp, "wprtp", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( thlm, "thlm", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( wpthlp, "wpthlp", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( wp2, "wp2", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( wp3, "wp3", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( rtp2, "rtp2", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( thlp2, "thlp2", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( rtpthlp, "rtpthlp", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)

    call check_nan( wpthlp_sfc, "wpthlp_sfc", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( wprtp_sfc, "wprtp_sfc", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( upwp_sfc, "upwp_sfc", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( vpwp_sfc, "vpwp_sfc", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( p_sfc, "p_sfc", prefix//proc_name, & ! intent(in)
                    err_info ) ! intent(inout)

    do sclr = 1, sclr_dim

      call check_nan( sclrm_forcing(:,sclr),"sclrm_forcing", & ! intent(in)
                      prefix//proc_name, & ! intent(in)
                      err_info ) ! intent(inout) ) ! intent(in)

      call check_nan( wpsclrp_sfc(sclr),"wpsclrp_sfc", & ! intent(in)
                      prefix//proc_name, & ! intent(in)
                      err_info ) ! intent(inout)

      call check_nan( sclrm(:,sclr),"sclrm", prefix//proc_name, & ! intent(in)
                      err_info ) ! intent(inout)
      call check_nan( wpsclrp(:,sclr),"wpsclrp", prefix//proc_name, & ! intent(in)
                      err_info ) ! intent(inout)
      call check_nan( sclrp2(:,sclr),"sclrp2", prefix//proc_name, & ! intent(in)
                      err_info ) ! intent(inout)
      call check_nan( sclrprtp(:,sclr),"sclrprtp", prefix//proc_name, & ! intent(in)
                      err_info ) ! intent(inout)
      call check_nan( sclrpthlp(:,sclr),"sclrpthlp", prefix//proc_name, & ! intent(in)
                      err_info ) ! intent(inout)

    end do


    do edsclr = 1, edsclr_dim

      call check_nan( edsclrm_forcing(:,edsclr),"edsclrm_forcing", prefix//proc_name, & ! intent(in)
                      err_info ) ! intent(inout)

      call check_nan( wpedsclrp_sfc(edsclr),"wpedsclrp_sfc", & ! intent(in)
                      prefix//proc_name, & ! intent(in)
                      err_info ) ! intent(inout)

      call check_nan( edsclrm(:,edsclr),"edsclrm", prefix//proc_name, & ! intent(in)
                      err_info ) ! intent(inout)

    enddo

!---------------------------------------------------------------------

    if ( clubb_at_least_debug_level_api( 0 ) ) then
        if ( any(err_info%err_code == clubb_fatal_error) ) then
            return
        end if
    end if

    call check_negative( rtm, 1, nzt, "rtm", prefix//proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( p_in_Pa, 1, nzt, "p_in_Pa", prefix//proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( rho, 1, nzt, "rho", prefix//proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( rho_zm, 1, nzm, "rho_zm", prefix//proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( exner, 1, nzt, "exner", prefix//proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( rho_ds_zm, 1, nzm, "rho_ds_zm", prefix//proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( rho_ds_zt, 1, nzt, "rho_ds_zt", prefix//proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( invrs_rho_ds_zm, 1, nzm, "invrs_rho_ds_zm", & ! intent(in)
                         prefix//proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( invrs_rho_ds_zt, 1, nzt, "invrs_rho_ds_zt", & ! intent(in)
                         prefix//proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( thv_ds_zm, 1, nzm, "thv_ds_zm", prefix//proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( thv_ds_zt, 1, nzt, "thv_ds_zt", prefix//proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( up2, 1, nzm, "up2", prefix//proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( vp2, 1, nzm, "vp2", prefix//proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( wp2, 1, nzm, "wp2", prefix//proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( thlm, 1, nzt, "thlm", prefix//proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( rtp2, 1, nzm, "rtp2", prefix//proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( thlp2, 1, nzm, "thlp2", prefix//proc_name, & ! intent(in)
                         err_info ) ! intent(inout)

    if ( any(err_info%err_code == clubb_fatal_error) .and. prefix == "beginning of " ) then
        ! Reset err_code (needs column nr!)
        err_info%err_code = clubb_no_error   ! Negative value generated by host model,
                                             ! hence ignore error
    end if

    ! Check the first levels for temperatures greater than 200K
    do k=1, min( 10, size(thlm) )
        if ( thlm(k) < 190. ) then
            write(fstderr,*) "Liquid water potential temperature (thlm) < 190K ", &
                             "at grid level k = ", k, ": thlm(",k,") = ", thlm(k)
        end if
    end do 

    return
  end subroutine parameterization_check

!-----------------------------------------------------------------------
  subroutine sfc_varnce_check( sclr_dim, wp2_sfc, up2_sfc, vp2_sfc, thlp2_sfc, &
           rtp2_sfc, rtpthlp_sfc, &
           sclrp2_sfc, sclrprtp_sfc, sclrpthlp_sfc, &
           err_info )
!
!       Description:This subroutine determines if any of the output
!       variables for the calc_surface_varnce subroutine carry values that
!       are nans.
!
!       Joshua Fasching February 2008
!
!
!-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use err_info_type_module, only: &
      err_info_type     ! Type

    implicit none

    ! Constant Parameters
    ! Name of the subroutine calling the check
    character(len=*), parameter :: &
      proc_name = "calc_surface_varnce"

    ! Input Variables
    integer, intent(in) :: &
      sclr_dim

    real( kind = core_rknd ),intent(in) :: &
      wp2_sfc,     & ! Vertical velocity variance        [m^2/s^2]
      up2_sfc,     & ! u'^2                              [m^2/s^2]
      vp2_sfc,     & ! u'^2                              [m^2/s^2]
      thlp2_sfc,   & ! thetal variance                   [K^2]
      rtp2_sfc,    & ! rt variance                       [(kg/kg)^2]
      rtpthlp_sfc    ! thetal rt covariance              [kg K/kg]


    real( kind = core_rknd ), dimension(sclr_dim), intent(in) :: &
      sclrp2_sfc,    & ! Passive scalar variance                 [units^2]
      sclrprtp_sfc,  & ! Passive scalar r_t covariance           [units kg/kg]
      sclrpthlp_sfc    ! Passive scalar theta_l covariance       [units K]

    ! Input/Output variables
    type(err_info_type), intent(inout) :: &
      err_info      ! err_info struct containing err_code and err_header
!-----------------------------------------------------------------------

    ! ---- Begin Code ----

    call check_nan( wp2_sfc, "wp2_sfc", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( up2_sfc, "up2_sfc", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( vp2_sfc, "vp2_sfc", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( thlp2_sfc, "thlp2_sfc", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( rtp2_sfc, "rtp2_sfc", proc_name, & ! intent(in)
                    err_info ) ! intent(inout)
    call check_nan( rtpthlp_sfc, "rtpthlp_sfc", & ! intent(in)
                    proc_name, & ! intent(in)
                    err_info ) ! intent(inout)

    if ( sclr_dim > 0 ) then
      call check_nan( sclrp2_sfc, "sclrp2_sfc", & ! intent(in)
                      proc_name, & ! intent(in)
                      err_info ) ! intent(inout)

      call check_nan( sclrprtp_sfc, "sclrprtp_sfc", & ! intent(in)
                      proc_name, & ! intent(in)
                      err_info ) ! intent(inout)

      call check_nan( sclrpthlp_sfc, "sclrpthlp_sfc", & ! intent(in)
                      proc_name, & ! intent(in)
                      err_info ) ! intent(inout)
    end if

    return
  end subroutine sfc_varnce_check

!-----------------------------------------------------------------------
  subroutine rad_check( nzm, nzt, thlm, rcm, rtm, rim, &
                        cloud_frac, p_in_Pa, exner, rho_zm, &
                        err_info )
! Description:
!   Checks radiation input variables. If they are < 0 it reports
!   to the console.
!------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use err_info_type_module, only: &
      err_info_type     ! Type

    implicit none

    integer, intent(in) :: &
      nzm, &
      nzt

    ! Constant Parameters
    character(len=*), parameter :: &
      proc_name = "Before BUGSrad."

    ! Input/Output variables
    real( kind = core_rknd ), dimension(nzt), intent(in) :: &
      thlm,           & ! Liquid Water Potential Temperature   [K/s]
      rcm,            & ! Liquid Water Mixing Ratio            [kg/kg]
      rtm,            & ! Total Water Mixing Ratio             [kg/kg]
      rim,            & ! Ice Water Mixing Ratio               [kg/kg]
      cloud_frac,     & ! Cloud Fraction                       [-]
      p_in_Pa,        & ! Pressure                             [Pa]
      exner             ! Exner Function                       [-]

    real( kind = core_rknd ), dimension(nzm), intent(in) :: &
      rho_zm            ! Air Density                          [kg/m^3]

    ! Input/Output variables
    type(err_info_type), intent(inout) :: &
      err_info      ! err_info struct containing err_code and err_header

    ! Local variables
    real( kind = core_rknd ), dimension(nzt) :: rvm

!-------------------------------------------------------------------------

    rvm = rtm - rcm

    call check_negative( thlm, 1, nzt, "thlm", proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( rcm, 1, nzt, "rcm", proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( rtm, 1, nzt, "rtm", proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( rvm, 1, nzt, "rvm", proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( rim, 1, nzt, "rim", proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( cloud_frac, 1, nzt,"cloud_frac", proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( p_in_Pa, 1, nzt, "p_in_Pa", proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( exner, 1, nzt, "exner", proc_name, & ! intent(in)
                         err_info ) ! intent(inout)
    call check_negative( rho_zm, 1, nzm, "rho_zm", proc_name, & ! intent(in)
                         err_info ) ! intent(inout)

    return

  end subroutine rad_check

!-----------------------------------------------------------------------
  logical function invalid_model_arrays( nzm, nzt, hydromet_dim, hydromet_list, &
                                         sclr_dim, edsclr_dim, &
                                         um, vm, rtm, wprtp, thlm, wpthlp, &
                                         rtp2, thlp2, rtpthlp, wp2, wp3, &
                                         wp2thvp, rtpthvp, thlpthvp, &
                                         hydromet, sclrm, edsclrm )

!       Description:
!       Checks for invalid floating point values in select model arrays.

!       References:
!       None
!------------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr   ! Constant(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none
    
    integer, intent(in) :: &
      nzm, &
      nzt, &
      hydromet_dim, &
      sclr_dim, &
      edsclr_dim

    character(len=10), dimension(hydromet_dim), intent(in) :: &
      hydromet_list

    real( kind = core_rknd ), dimension(nzt), intent(in) :: &
      um,       & ! eastward grid-mean wind comp. (thermo. levs.)  [m/s]
      vm,       & ! northward grid-mean wind comp. (thermo. levs.) [m/s]
      rtm,      & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      thlm,     & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      wp3,      & ! w'^3 (thermodynamic levels)                    [m^3/s^3]
      wp2thvp     ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]

    real( kind = core_rknd ), dimension(nzm), intent(in) :: &
      wprtp,    & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
      wpthlp,   & ! w'th_l' (momentum levels)                      [(m/s) K]
      rtp2,     & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      thlp2,    & ! th_l'^2 (momentum levels)                      [K^2]
      rtpthlp,  & ! r_t'th_l' (momentum levels)                    [(kg/kg) K]
      wp2,      & ! w'^2 (momentum levels)                         [m^2/s^2]
      rtpthvp,  & ! < r_t' th_v' > (momentum levels)               [kg/kg K]
      thlpthvp    ! < th_l' th_v' > (momentum levels)              [K^2]

    real( kind = core_rknd ), dimension(nzt,hydromet_dim), intent(in) :: &
      hydromet    ! Array of hydrometeors                          [units vary]

    real( kind = core_rknd ), dimension(nzt,sclr_dim), intent(in) :: &
      sclrm    ! Passive scalar mean (thermo. levels)              [units vary]

    real( kind = core_rknd ), dimension(nzt,edsclr_dim), intent(in) :: &
      edsclrm   ! Eddy passive scalar grid-mean (thermo. levels)   [units vary]

    ! Local Variables
    integer :: i, sclr, edsclr

    invalid_model_arrays = .false.

    ! Check whether any variable array contains a NaN for
    ! um, vm, thlm, rtm, rtp2, thlp2, wprtp, wpthlp, rtpthlp,
    ! wp2, & wp3.
    if ( is_nan_2d( um ) ) then
      write(fstderr,*) "NaN in um model array"
!         write(fstderr,*) "um= ", um
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( vm ) ) then
      write(fstderr,*) "NaN in vm model array"
!         write(fstderr,*) "vm= ", vm
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( wp2 ) ) then
      write(fstderr,*) "NaN in wp2 model array"
!         write(fstderr,*) "wp2= ", wp2
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( wp3 ) ) then
      write(fstderr,*) "NaN in wp3 model array"
!         write(fstderr,*) "wp3= ", wp3
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( rtm ) ) then
      write(fstderr,*) "NaN in rtm model array"
!         write(fstderr,*) "rtm= ", rtm
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( thlm ) ) then
      write(fstderr,*) "NaN in thlm model array"
!         write(fstderr,*) "thlm= ", thlm
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( rtp2 ) ) then
      write(fstderr,*) "NaN in rtp2 model array"
!         write(fstderr,*) "rtp2= ", rtp2
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( thlp2 ) ) then
      write(fstderr,*) "NaN in thlp2 model array"
!         write(fstderr,*) "thlp2= ", thlp2
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( wprtp ) ) then
      write(fstderr,*) "NaN in wprtp model array"
!         write(fstderr,*) "wprtp= ", wprtp
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( wpthlp ) ) then
      write(fstderr,*) "NaN in wpthlp model array"
!         write(fstderr,*) "wpthlp= ", wpthlp
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( rtpthlp ) ) then
      write(fstderr,*) "NaN in rtpthlp model array"
!         write(fstderr,*) "rtpthlp= ", rtpthlp
      invalid_model_arrays = .true.
!         return
    end if

    if ( hydromet_dim > 0 ) then
      do i = 1, hydromet_dim, 1
        if ( is_nan_2d( hydromet(:,i) ) ) then
          write(fstderr,*) "NaN in a hydrometeor model array "// &
            trim( hydromet_list(i) )
!             write(fstderr,*) "hydromet= ", hydromet
          invalid_model_arrays = .true.
!             return
        end if
      end do
    end if

!       if ( is_nan_2d( wm_zt ) ) then
!         write(fstderr,*) "NaN in wm_zt model array"
!         write(fstderr,*) "wm_zt= ", wm_zt
!         invalid_model_arrays = .true.
!         return
!       end if

    if ( is_nan_2d( wp2thvp ) ) then
      write(fstderr,*) "NaN in wp2thvp model array"
!         write(fstderr,*) "wp2thvp = ", wp2thvp
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( rtpthvp ) ) then
      write(fstderr,*) "NaN in rtpthvp model array"
!         write(fstderr,*) "rtpthvp = ", rtpthvp
      invalid_model_arrays = .true.
    end if

    if ( is_nan_2d( thlpthvp ) ) then
      write(fstderr,*) "NaN in thlpthvp model array"
!         write(fstderr,*) "thlpthvp = ", thlpthvp
      invalid_model_arrays = .true.
    end if

    do sclr = 1, sclr_dim, 1
      if ( is_nan_2d( sclrm(:,sclr) ) ) then
        write(fstderr,*) "NaN in sclrm", sclr, "model array"
!           write(fstderr,'(a6,i2,a1)') "sclrm(", sclr, ")"
!           write(fstderr,*) sclrm(:,sclr)
        invalid_model_arrays = .true.
      end if
    end do

    do edsclr = 1, edsclr_dim, 1
      if ( is_nan_2d( edsclrm(:,edsclr) ) ) then
        write(fstderr,*) "NaN in edsclrm", edsclr, "model array"
!           write(fstderr,'(a8,i2,a1)') "edsclrm(", edsclr, ")"
!           write(fstderr,*) edsclrm(:,edsclr)
        invalid_model_arrays = .true.
      end if
    end do

    return
  end function invalid_model_arrays

!------------------------------------------------------------------------
  pure logical function is_nan_sclr( xarg )

! Description:
!   Checks if a given scalar real is a NaN, +inf or -inf.

! Notes:
!   I was advised by Andy Vaught to use a data statement and the transfer( )
!   intrinsic rather than using a hex number in a parameter for portability.

!   Certain compiler optimizations may cause variables with invalid
!   results to flush to zero.  Avoid these!
!  -dschanen 16 Dec 2010

!------------------------------------------------------------------------

    use, intrinsic :: ieee_arithmetic

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: xarg

    ! ---- Begin Code ---

    if (.not. ieee_is_finite(xarg) .or. ieee_is_nan(xarg)) then
      ! Try ieee_is_finite ieee_is_nan 
      is_nan_sclr = .true.
    else
      is_nan_sclr = .false.
    end if


    return
  end function is_nan_sclr
!------------------------------------------------------------------------

!------------------------------------------------------------------------
  pure logical function is_nan_2d( x2d )

! Description:
!   Checks if a given real vector is a NaN, +inf or -inf.

!------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: any

    ! Input Variables
    real( kind = core_rknd ), dimension(:), intent(in) :: x2d

    ! Local Variables
    integer :: k

    ! ---- Begin Code ----

    is_nan_2d = .false.

    do k = 1, size( x2d )
      if ( is_nan_sclr( x2d(k) ) ) then
        is_nan_2d = .true.
        exit
      end if
    end do

    return

  end function is_nan_2d


!------------------------------------------------------------------------
  subroutine check_negative_index( var, varstart, varend, varname, operation, err_info )
!
! Description:
!   Checks for negative values in the var array and reports
!   the index in which the negative values occur.
!
!-----------------------------------------------------------------------
    use constants_clubb, only: &
        fstderr ! Variable

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_fatal_error              ! Constant

    use err_info_type_module, only: &
      err_info_type     ! Type

    implicit none

    real( kind = core_rknd ), intent(in) :: var(:)

    integer, intent(in) :: varstart, varend

    character(len=*), intent(in):: &
    varname,     & ! Varible being examined
    operation      ! Procedure calling check_zero

    ! Input/Output variables
    type(err_info_type), intent(inout) :: &
      err_info      ! err_info struct containing err_code and err_header

    ! Local Variable
    integer :: k ! Loop iterator

    do k = varstart, varend

      if ( var(k) < 0.0_core_rknd ) then
        write(fstderr,*) err_info%err_header_global
        write(fstderr,*) varname, " < 0 in ", operation, &
                         " at k = ", k
        ! Set to clubb_fatal_error in single col (col index not available!!)
        err_info%err_code = clubb_fatal_error

      end if

    end do

    return

  end subroutine check_negative_index


!------------------------------------------------------------------------
  subroutine check_nan_2d( var, varname, operation, err_info )
!
!  Description:
!    Checks for a NaN in the var array and reports it.
!
!
!------------------------------------------------------------------------
    use constants_clubb, only: &
        fstderr ! Variable(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_fatal_error              ! Constant

    use err_info_type_module, only: &
      err_info_type     ! Type

    implicit none

    ! External
    intrinsic :: present

    ! Input variables
    real( kind = core_rknd ), intent(in), dimension(:) :: var ! Variable being examined

    character(len=*), intent(in):: &
      varname,  & ! Name of variable
      operation   ! Procedure calling check_nan

    ! Input/Output variables
    type(err_info_type), intent(inout) :: &
      err_info      ! err_info struct containing err_code and err_header

    if ( is_nan_2d( var ) ) then
      write(fstderr,*) err_info%err_header_global
      write(fstderr,*) varname, " is NaN in ",operation
      ! Single col error, but col nr is not available
      err_info%err_code = clubb_fatal_error
    end if

    return
  end subroutine check_nan_2d

!-----------------------------------------------------------------------
  subroutine check_nan_sclr( var, varname, operation, err_info )
!
! Description:
!   Checks for a NaN in the scalar var then reports it.
!
!-----------------------------------------------------------------------
    use constants_clubb, only: &
        fstderr ! Variable

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_fatal_error              ! Constant

    use err_info_type_module, only: &
      err_info_type     ! Type

    implicit none

    ! External
    intrinsic :: present

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: var        ! Variable being examined

    character(len=*), intent(in):: &
      varname,    & ! Name of variable being examined
      operation     ! Procedure calling check_nan

    ! Input/Output variables
    type(err_info_type), intent(inout) :: &
      err_info      ! err_info struct containing err_code and err_header

!--------------------------------------------------------------------
    if ( is_nan_sclr( var ) ) then
      write(fstderr,*) err_info%err_header_global
      write(fstderr,*) varname, " is NaN in ",operation
      ! Single col error, but col nr is not available
      err_info%err_code = clubb_fatal_error
    end if

    return

  end subroutine check_nan_sclr
!-------------------------------------------------------------------------

!-----------------------------------------------------------------------
  function calculate_spurious_source( integral_after, integral_before, &
                                           flux_top, flux_sfc, &
                                           integral_forcing, dt ) &
  result( spurious_source )
!
! Description:
!   Checks whether there is conservation within the column and returns any
!   imbalance as spurious_source where spurious_source is defined negative
!   for a spurious sink.
!
!-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      integral_after, &   ! Vertically-integrated quantity after dt time  [units vary]
      integral_before, &  ! Vertically-integrated quantity before dt time [units vary]
      flux_top, &         ! Total flux at the top of the domain           [units vary]
      flux_sfc, &         ! Total flux at the bottom of the domain        [units vary]
      integral_forcing, & ! Vertically-integrated forcing                 [units vary]
      dt                  ! Timestep size                                 [s]

    ! Return Variable
    real( kind = core_rknd ) :: spurious_source ! [units vary]

!--------------------------------------------------------------------

    ! ---- Begin Code ----

    spurious_source = (integral_after - integral_before) / dt &
                        + flux_top - flux_sfc - integral_forcing

    return

  end function calculate_spurious_source
!-------------------------------------------------------------------------

  !=============================================================================
  subroutine check_clubb_settings_api( ngrdcol,             & ! intent(in)
                                       params,              & ! intent(in)
                                       l_implemented,       & ! intent(in)
                                       l_input_fields,      & ! intent(in)
                                       clubb_config_flags,  & ! intent(in)
                                       err_info )             ! intent(inout)

      ! Description:
      !   Subroutine to set up the model for execution.
      !
      ! References:
      !   None
      !---------------------------------------------------------------------

      use parameter_indices, only:  &
          nparams,      & ! Variable(s)
          iC1,          & ! Constant(s)
          iC1b,         &
          iC2rt,        &
          iC2thl,       &
          iC2rtthl,     &
          iC6rt,        &
          iC6rtb,       &
          iC6thl,       &
          iC6thlb,      &
          iC14

      use parameters_tunable, only: &
          nu_vertical_res_dep    ! Type(s)

      use constants_clubb, only:  &
          fstderr, &  ! Variable(s)
          one, &
          eps

      use error_code, only: &
          clubb_at_least_debug_level_api,  & ! Procedures
          initialize_error_headers,    &
          clubb_fatal_error              ! Constant

      use model_flags, only: &
          clubb_config_flags_type, & ! Type
          iiPDF_ADG1,       & ! Variable(s)
          iiPDF_ADG2,       &
          iiPDF_3D_Luhar,   &
          iiPDF_new,        &
          iiPDF_TSDADG,     &
          iiPDF_LY93,       &
          iiPDF_new_hybrid, &
          saturation_bolton,  &
          saturation_gfdl,    &
          saturation_flatau,  &
          saturation_lookup,  &
          ipdf_pre_advance_fields, &     
          ipdf_post_advance_fields, &    
          ipdf_pre_post_advance_fields, &
          l_explicit_turbulent_adv_wpxp, &
          order_xm_wpxp, &
          order_xp2_xpyp, &
          order_wp2_wp3, &
          order_windm
          
#ifdef CLUBB_GPU
      use model_flags, only: &
          lapack
#endif

      use clubb_precision, only: &
          core_rknd ! Variable(s)

      use sponge_layer_damping, only: &
          thlm_sponge_damp_settings,    & ! Variable(s)
          rtm_sponge_damp_settings,     &
          uv_sponge_damp_settings,      &
          wp2_sponge_damp_settings,     &
          wp3_sponge_damp_settings,     &
          up2_vp2_sponge_damp_settings


      use err_info_type_module, only: &
        err_info_type        ! Type

      implicit none

      !---------------------- Input Variables ----------------------

      integer, intent(in) :: &
        ngrdcol

      real( kind = core_rknd ), intent(in), dimension(ngrdcol,nparams) :: &
        params  ! Including C1, nu1, nu2, etc.

      ! Flags
      logical, intent(in) ::  &
        l_implemented, & ! Flag for whether CLUBB is being run in a host model
        l_input_fields   ! Flag for whether LES input fields are being used

      type(clubb_config_flags_type), intent(in) :: &
        clubb_config_flags

      !---------------------- InOut Variables ----------------------
      type(err_info_type), intent(inout) :: &
        err_info        ! err_info struct containing err_code and err_header

      !---------------------- Begin Code ----------------------

      call initialize_error_headers

#ifdef CLUBB_GPU
      if ( clubb_config_flags%penta_solve_method == lapack ) then
        write(fstderr,*) "WARNING: The penta-diagonal lapack solver is not GPU accelerated"
        write(fstderr,*) " Set penta_solve_method = 2, to use an accelerated penta-diagonal solver"
      end if

      if ( clubb_config_flags%tridiag_solve_method == lapack ) then
        write(fstderr,*) "WARNING: The tri-diagonal lapack solver is not GPU accelerated"
        write(fstderr,*) " Set tridiag_solve_method = 2, to use an accelerated tri-diagonal solver"
      end if

      if ( l_input_fields ) then
        error stop "l_input_fields = .true. not usable when running on GPUs"
      end if
#endif

      ! Sanity check
      if ( clubb_at_least_debug_level_api( 0 ) ) then

        if ( clubb_config_flags%l_damp_wp2_using_em &
           .and. ( any( abs(params(:,iC1) - params(:,iC14)) > &
                        abs(params(:,iC1) + params(:,iC14)) /  2 * eps ) &
                   .or. clubb_config_flags%l_stability_correct_tau_zm ) ) then
          write(fstderr, *) err_info%err_header_global
          write(fstderr,*) "l_damp_wp2_using_em = T requires C1=C14 and" &
                            // " l_stability_correct_tau_zm = F"
          write(fstderr,*) "C1 = ", params(:,iC1)
          write(fstderr,*) "C14 = ", params(:,iC14)
          write(fstderr,*) "l_stability_correct_tau_zm = ", &
                           clubb_config_flags%l_stability_correct_tau_zm
          write(fstderr,*) "Fatal error in check_clubb_settings_api"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
          return
        end if

      end if

      ! Sanity check for the saturation formula
      select case ( clubb_config_flags%saturation_formula )
      case ( saturation_bolton )
        ! Using the Bolton 1980 approximations for SVP over vapor/ice

      case ( saturation_gfdl )
        ! Using the Flatau, et al. polynomial approximation for SVP over vapor/ice

      case ( saturation_flatau )   ! h1g, 2010-06-16
        ! Using the GFDL SVP formula (Goff-Gratch)

        ! Add new saturation formulas after this

      case ( saturation_lookup )
        ! Using the lookup table

      case default
        write(fstderr, *) err_info%err_header_global
        write(fstderr,*) "Unknown approx. of saturation vapor pressure: ", &
           clubb_config_flags%saturation_formula
        write(fstderr,*) "Fatal error in check_clubb_settings_api"
        ! General error -> set all entries to clubb_fatal_error
        err_info%err_code = clubb_fatal_error
        return
      end select

      ! Check for the type of two component normal (double Gaussian) PDF being
      ! used for w, rt, and theta-l (or w, chi, and eta).
      if ( clubb_config_flags%iiPDF_type < iiPDF_ADG1 &
           .or. clubb_config_flags%iiPDF_type > iiPDF_new_hybrid ) then
         write(fstderr, *) err_info%err_header_global
         write(fstderr,*) "Unknown type of double Gaussian PDF selected: ", &
                          clubb_config_flags%iiPDF_type
         write(fstderr,*) "iiPDF_type = ", clubb_config_flags%iiPDF_type
         write(fstderr,*) "Fatal error in check_clubb_settings_api"
         ! General error -> set all entries to clubb_fatal_error
         err_info%err_code = clubb_fatal_error
         return
      endif ! iiPDF_type < iiPDF_ADG1 or iiPDF_type > iiPDF_lY93

      ! The ADG2 and 3D Luhar PDFs can only be used as part of input fields.
      if ( clubb_config_flags%iiPDF_type == iiPDF_ADG2 ) then
         if ( .not. l_input_fields ) then
            write(fstderr, *) err_info%err_header_global
            write(fstderr,*) "The ADG2 PDF can only be used with" &
                             // " input fields (l_input_fields = .true.)."
            write(fstderr,*) "iiPDF_type = ", clubb_config_flags%iiPDF_type
            write(fstderr,*) "l_input_fields = ", l_input_fields
            write(fstderr,*) "Fatal error in check_clubb_settings_api"
            ! General error -> set all entries to clubb_fatal_error
            err_info%err_code = clubb_fatal_error
            return
         endif ! .not. l_input_fields
      endif ! iiPDF_type == iiPDF_ADG2

      if ( clubb_config_flags%iiPDF_type == iiPDF_3D_Luhar ) then
         if ( .not. l_input_fields ) then
            write(fstderr, *) err_info%err_header_global
            write(fstderr,*) "The 3D Luhar PDF can only be used with" &
                             // " input fields (l_input_fields = .true.)."
            write(fstderr,*) "iiPDF_type = ", clubb_config_flags%iiPDF_type
            write(fstderr,*) "l_input_fields = ", l_input_fields
            write(fstderr,*) "Fatal error in check_clubb_settings_api"
            ! General error -> set all entries to clubb_fatal_error
            err_info%err_code = clubb_fatal_error
            return
         endif ! .not. l_input_fields
      endif ! iiPDF_type == iiPDF_3D_Luhar

      ! This also currently applies to the new PDF until it has been fully
      ! implemented.
      if ( clubb_config_flags%iiPDF_type == iiPDF_new ) then
         if ( .not. l_input_fields ) then
            write(fstderr, *) err_info%err_header_global
            write(fstderr,*) "The new PDF can only be used with" &
                             // " input fields (l_input_fields = .true.)."
            write(fstderr,*) "iiPDF_type = ", clubb_config_flags%iiPDF_type
            write(fstderr,*) "l_input_fields = ", l_input_fields
            write(fstderr,*) "Fatal error in check_clubb_settings_api"
            ! General error -> set all entries to clubb_fatal_error
            err_info%err_code = clubb_fatal_error
            return
         endif ! .not. l_input_fields
      endif ! iiPDF_type == iiPDF_new

      ! This also currently applies to the TSDADG PDF until it has been fully
      ! implemented.
      if ( clubb_config_flags%iiPDF_type == iiPDF_TSDADG ) then
         if ( .not. l_input_fields ) then
            write(fstderr, *) err_info%err_header_global
            write(fstderr,*) "The new TSDADG PDF can only be used with" &
                             // " input fields (l_input_fields = .true.)."
            write(fstderr,*) "iiPDF_type = ", clubb_config_flags%iiPDF_type
            write(fstderr,*) "l_input_fields = ", l_input_fields
            write(fstderr,*) "Fatal error in check_clubb_settings_api"
            ! General error -> set all entries to clubb_fatal_error
            err_info%err_code = clubb_fatal_error
            return
         endif ! .not. l_input_fields
      endif ! iiPDF_type == iiPDF_TSDADG

      ! This also applies to Lewellen and Yoh (1993).
      if ( clubb_config_flags%iiPDF_type == iiPDF_LY93 ) then
         if ( .not. l_input_fields ) then
            write(fstderr, *) err_info%err_header_global
            write(fstderr,*) "The Lewellen and Yoh PDF can only be used with" &
                             // " input fields (l_input_fields = .true.)."
            write(fstderr,*) "iiPDF_type = ", clubb_config_flags%iiPDF_type
            write(fstderr,*) "l_input_fields = ", l_input_fields
            write(fstderr,*) "Fatal error in check_clubb_settings_api"
            ! General error -> set all entries to clubb_fatal_error
            err_info%err_code = clubb_fatal_error
            return
         endif ! .not. l_input_fields
      endif ! iiPDF_type == iiPDF_LY93

      ! Check the option for the placement of the call to CLUBB's PDF.
      if ( clubb_config_flags%ipdf_call_placement < ipdf_pre_advance_fields &
           .or. clubb_config_flags%ipdf_call_placement > ipdf_pre_post_advance_fields ) then
         write(fstderr, *) err_info%err_header_global
         write(fstderr,*) "Invalid option selected for ipdf_call_placement: ", &
                          clubb_config_flags%ipdf_call_placement
         write(fstderr,*) "Fatal error in check_clubb_settings_api"
         ! General error -> set all entries to clubb_fatal_error
         err_info%err_code = clubb_fatal_error
         return
      endif

      ! The l_predict_upwp_vpwp flag requires that the ADG1 PDF is used
      ! implicitly in subroutine advance_xm_wpxp.
      if ( clubb_config_flags%l_predict_upwp_vpwp ) then

         ! When l_predict_upwp_vpwp is enabled, the
         ! l_explicit_turbulent_adv_wpxp flag must be turned off.
         ! Otherwise, explicit turbulent advection would require PDF parameters
         ! for u and v to be calculated in PDF closure.  These would be needed
         ! to calculate integrated fields such as wp2up, etc.
         if ( l_explicit_turbulent_adv_wpxp ) then
            write(fstderr, *) err_info%err_header_global
            write(fstderr,*) "The l_explicit_turbulent_adv_wpxp option" &
                             // " is not currently set up for use with the" &
                             // " l_predict_upwp_vpwp code."
            write(fstderr,*) "Fatal error in check_clubb_settings_api"
            ! General error -> set all entries to clubb_fatal_error
            err_info%err_code = clubb_fatal_error
            return
         endif ! l_explicit_turbulent_adv_wpxp

         ! When l_predict_upwp_vpwp is enabled, the PDF type must be set to
         ! the ADG1 PDF or the new hybrid PDF.  The other PDFs are not currently
         ! set up to calculate variables needed for implicit or semi-implicit
         ! turbulent advection, such as coef_wp2up_implicit, etc.
         if ( ( clubb_config_flags%iiPDF_type /= iiPDF_ADG1 ) &
              .and. ( clubb_config_flags%iiPDF_type /= iiPDF_new_hybrid ) ) then
            write(fstderr, *) err_info%err_header_global
            write(fstderr,*) "Currently, only the ADG1 PDF and the new hybrid" &
                             // " PDF are set up for use with the" &
                             // " l_predict_upwp_vpwp code."
            write(fstderr,*) "Fatal error in check_clubb_settings_api"
            ! General error -> set all entries to clubb_fatal_error
            err_info%err_code = clubb_fatal_error
            return
         endif ! iiPDF_type /= iiPDF_ADG1

      endif ! l_predict_upwp_vpwp

      ! The flags l_min_xp2_from_corr_wx and l_enable_relaxed_clipping must
      ! have opposite values.
      if ( ( clubb_config_flags%l_min_xp2_from_corr_wx ) &
         .and. ( clubb_config_flags%l_enable_relaxed_clipping ) ) then
         write(fstderr, *) err_info%err_header_global
         write(fstderr,*) "Invalid configuration: l_min_xp2_from_corr_wx = T " &
                          // "and l_enable_relaxed_clipping = T"
         write(fstderr,*) "They must have opposite values"
         write(fstderr,*) "Fatal error in check_clubb_settings_api"
         ! General error -> set all entries to clubb_fatal_error
         err_info%err_code = clubb_fatal_error
         return
      elseif ( ( .not. clubb_config_flags%l_min_xp2_from_corr_wx ) &
               .and. ( .not. clubb_config_flags%l_enable_relaxed_clipping ) ) then
         write(fstderr, *) err_info%err_header_global
         write(fstderr,*) "Invalid configuration: l_min_xp2_from_corr_wx = F " &
                          // "and l_enable_relaxed_clipping = F"
         write(fstderr,*) "They must have opposite values"
         write(fstderr,*) "Fatal error in check_clubb_settings_api"
         ! General error -> set all entries to clubb_fatal_error
         !err_info%err_code = clubb_fatal_error
         !return
      endif

      ! Checking for the code that orders CLUBB's advance_ subroutines
      if ( order_xm_wpxp < 1 .or. order_xm_wpxp > 4 ) then
         write(fstderr, *) err_info%err_header_global
         write(fstderr,*) "The variable order_xm_wpxp must have a value " &
                          // "between 1 and 4"
         write(fstderr,*) "order_xm_wpxp = ", order_xm_wpxp
         write(fstderr,*) "Fatal error in check_clubb_settings_api"
         ! General error -> set all entries to clubb_fatal_error
         err_info%err_code = clubb_fatal_error
         return
      elseif ( order_xm_wpxp == order_wp2_wp3 &
               .or. order_xm_wpxp == order_xp2_xpyp &
               .or. order_xm_wpxp == order_windm ) then
         write(fstderr, *) err_info%err_header_global
         write(fstderr,*) "The variable order_xm_wpxp has the same value " &
                          // "as another order_ variable.  Please give each " &
                          // "order index a unique value."
         write(fstderr,*) "order_xm_wpxp = ", order_xm_wpxp
         write(fstderr,*) "order_wp2_wp3 = ", order_wp2_wp3
         write(fstderr,*) "order_xp2_xpyp = ", order_xp2_xpyp
         write(fstderr,*) "order_windm = ", order_windm
         write(fstderr,*) "Fatal error in check_clubb_settings_api"
         ! General error -> set all entries to clubb_fatal_error
         err_info%err_code = clubb_fatal_error
         return
      endif

      if ( order_wp2_wp3 < 1 .or. order_wp2_wp3 > 4 ) then
         write(fstderr, *) err_info%err_header_global
         write(fstderr,*) "The variable order_wp2_wp3 must have a value " &
                          // "between 1 and 4"
         write(fstderr,*) "order_wp2_wp3 = ", order_wp2_wp3
         write(fstderr,*) "Fatal error in check_clubb_settings_api"
         ! General error -> set all entries to clubb_fatal_error
         err_info%err_code = clubb_fatal_error
         return
      elseif ( order_wp2_wp3 == order_xm_wpxp &
               .or. order_wp2_wp3 == order_xp2_xpyp &
               .or. order_wp2_wp3 == order_windm ) then
         write(fstderr, *) err_info%err_header_global
         write(fstderr,*) "The variable order_wp2_wp3 has the same value " &
                          // "as another order_ variable.  Please give each " &
                          // "order index a unique value."
         write(fstderr,*) "order_wp2_wp3 = ", order_wp2_wp3
         write(fstderr,*) "order_xm_wpxp = ", order_xm_wpxp
         write(fstderr,*) "order_xp2_xpyp = ", order_xp2_xpyp
         write(fstderr,*) "order_windm = ", order_windm
         write(fstderr,*) "Fatal error in check_clubb_settings_api"
         ! General error -> set all entries to clubb_fatal_error
         err_info%err_code = clubb_fatal_error
         return
      endif

      if ( order_xp2_xpyp < 1 .or. order_xp2_xpyp > 4 ) then
         write(fstderr, *) err_info%err_header_global
         write(fstderr,*) "The variable order_xp2_xpyp must have a value " &
                          // "between 1 and 4"
         write(fstderr,*) "order_xp2_xpyp = ", order_xp2_xpyp
         write(fstderr,*) "Fatal error in check_clubb_settings_api"
         ! General error -> set all entries to clubb_fatal_error
         err_info%err_code = clubb_fatal_error
         return
      elseif ( order_xp2_xpyp == order_wp2_wp3 &
               .or. order_xp2_xpyp == order_xm_wpxp &
               .or. order_xp2_xpyp == order_windm ) then
         write(fstderr, *) err_info%err_header_global
         write(fstderr,*) "The variable order_xp2_xpyp has the same value " &
                          // "as another order_ variable.  Please give each " &
                          // "order index a unique value."
         write(fstderr,*) "order_xp2_xpyp = ", order_xp2_xpyp
         write(fstderr,*) "order_wp2_wp3 = ", order_wp2_wp3
         write(fstderr,*) "order_xm_wpxp = ", order_xm_wpxp
         write(fstderr,*) "order_windm = ", order_windm
         write(fstderr,*) "Fatal error in check_clubb_settings_api"
         ! General error -> set all entries to clubb_fatal_error
         err_info%err_code = clubb_fatal_error
         return
      endif

      if ( order_windm < 1 .or. order_windm > 4 ) then
         write(fstderr, *) err_info%err_header_global
         write(fstderr,*) "The variable order_windm must have a value " &
                          // "between 1 and 4"
         write(fstderr,*) "order_windm = ", order_windm
         write(fstderr,*) "Fatal error in check_clubb_settings_api"
         ! General error -> set all entries to clubb_fatal_error
         err_info%err_code = clubb_fatal_error
         return
      elseif ( order_windm == order_wp2_wp3 &
               .or. order_windm == order_xp2_xpyp &
               .or. order_windm == order_xm_wpxp ) then
         write(fstderr, *) err_info%err_header_global
         write(fstderr,*) "The variable order_windm has the same value " &
                          // "as another order_ variable.  Please give each " &
                          // "order index a unique value."
         write(fstderr,*) "order_windm = ", order_windm
         write(fstderr,*) "order_wp2_wp3 = ", order_wp2_wp3
         write(fstderr,*) "order_xp2_xpyp = ", order_xp2_xpyp
         write(fstderr,*) "order_xm_wpxp = ", order_xm_wpxp
         write(fstderr,*) "Fatal error in check_clubb_settings_api"
         ! General error -> set all entries to clubb_fatal_error
         err_info%err_code = clubb_fatal_error
         return
      endif

      ! Checking that when the l_diag_Lscale_from_tau is enabled, the
      ! relevant Cx tunable parameters are all set to a value of 1 (as
      ! you're supposed to tune the C_invrs_tau_ parameters instead).
      if ( clubb_config_flags%l_diag_Lscale_from_tau ) then

         ! Note: someday when we can successfully run with all these parameters
         ! having a value of 1, the "Warning" messages should be removed and the
         ! "Fatal error" messages should be uncommented.

         ! C1 must have a value of 1
         if ( any(params(:,iC1) > one .or. params(:,iC1) < one) ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C1 must have a value of 1."
            write(fstderr,*) "C1 = ", params(:,iC1)
            write(fstderr,*) "Warning in check_clubb_settings_api"
            !write(fstderr,*) "Fatal error in check_clubb_settings_api"
            ! General error -> set all entries to clubb_fatal_error
            !err_info%err_code = clubb_fatal_error
         endif ! C1 check

         ! C1b must have a value of 1
         if ( any(params(:,iC1b) > one .or. params(:,iC1b) < one) ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C1b must have a value of 1."
            write(fstderr,*) "C1b = ", params(:,iC1b)
            write(fstderr,*) "Warning in check_clubb_settings_api"
            !write(fstderr,*) "Fatal error in check_clubb_settings_api"
            ! General error -> set all entries to clubb_fatal_error
            !err_info%err_code = clubb_fatal_error
         endif ! C1b check

         ! C2rt must have a value of 1
         if ( any(params(:,iC2rt) > one .or. params(:,iC2rt) < one) ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C2rt must have a value of 1."
            write(fstderr,*) "C2rt = ", params(:,iC2rt)
            write(fstderr,*) "Warning in check_clubb_settings_api"
            !write(fstderr,*) "Fatal error in check_clubb_settings_api"
            ! General error -> set all entries to clubb_fatal_error
            !err_info%err_code = clubb_fatal_error
         endif ! C2rt check

         ! C2thl must have a value of 1
         if ( any(params(:,iC2thl) > one .or. params(:,iC2thl) < one) ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C2thl must have a value of 1."
            write(fstderr,*) "C2thl = ", params(:,iC2thl)
            write(fstderr,*) "Warning in check_clubb_settings_api"
            !write(fstderr,*) "Fatal error in check_clubb_settings_api"
            ! General error -> set all entries to clubb_fatal_error
            !err_info%err_code = clubb_fatal_error
         endif ! C2thl check

         ! C2rtthl must have a value of 1
         if ( any(params(:,iC2rtthl) > one .or. params(:,iC2rtthl) < one) ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C2rtthl must have a value of 1."
            write(fstderr,*) "C2rtthl = ", params(:,iC2rtthl)
            write(fstderr,*) "Warning in check_clubb_settings_api"
            !write(fstderr,*) "Fatal error in check_clubb_settings_api"
            ! General error -> set all entries to clubb_fatal_error
            !err_info%err_code = clubb_fatal_error
         endif ! C2rtthl check

         ! C6rt must have a value of 1
         if ( any(params(:,iC6rt) > one .or. params(:,iC6rt) < one) ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C6rt must have a value of 1."
            write(fstderr,*) "C6rt = ", params(:,iC6rt)
            write(fstderr,*) "Warning in check_clubb_settings_api"
            !write(fstderr,*) "Fatal error in check_clubb_settings_api"
            ! General error -> set all entries to clubb_fatal_error
            !err_info%err_code = clubb_fatal_error
         endif ! C6rt check

         ! C6rtb must have a value of 1
         if ( any(params(:,iC6rtb) > one .or. params(:,iC6rtb) < one) ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C6rtb must have a value of 1."
            write(fstderr,*) "C6rtb = ", params(:,iC6rtb)
            write(fstderr,*) "Warning in check_clubb_settings_api"
            !write(fstderr,*) "Fatal error in check_clubb_settings_api"
            ! General error -> set all entries to clubb_fatal_error
            !err_info%err_code = clubb_fatal_error
         endif ! C6rtb check

         ! C6thl must have a value of 1
         if ( any(params(:,iC6thl) > one .or. params(:,iC6thl) < one) ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C6thl must have a value of 1."
            write(fstderr,*) "C6thl = ", params(:,iC6thl)
            write(fstderr,*) "Warning in check_clubb_settings_api"
            !write(fstderr,*) "Fatal error in check_clubb_settings_api"
            ! General error -> set all entries to clubb_fatal_error
            !err_info%err_code = clubb_fatal_error
         endif ! C6thl check

         ! C6thlb must have a value of 1
         if ( any(params(:,iC6thlb) > one .or. params(:,iC6thlb) < one) ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C6thlb must have a value of 1."
            write(fstderr,*) "C6thlb = ", params(:,iC6thlb)
            write(fstderr,*) "Warning in check_clubb_settings_api"
            !write(fstderr,*) "Fatal error in check_clubb_settings_api"
            ! General error -> set all entries to clubb_fatal_error
            !err_info%err_code = clubb_fatal_error
         endif ! C6thlb check

         ! C14 must have a value of 1
         if ( any(params(:,iC14) > one .or. params(:,iC14) < one) ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C14 must have a value of 1."
            write(fstderr,*) "C14 = ", params(:,iC14)
            write(fstderr,*) "Warning in check_clubb_settings_api"
            !write(fstderr,*) "Fatal error in check_clubb_settings_api"
            ! General error -> set all entries to clubb_fatal_error
            !err_info%err_code = clubb_fatal_error
         endif ! C14 check

      endif ! l_diag_Lscale_from_tau

      if ( l_implemented ) then

        if ( clubb_config_flags%l_rtm_nudge ) then
          write(fstderr, *) err_info%err_header_global
          write(fstderr,*) "l_rtm_nudge must be set to .false. when " &
                           // "l_implemented = .true."
          write(fstderr,*) "Fatal error in check_clubb_settings_api"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

        if ( clubb_config_flags%l_uv_nudge ) then
          write(fstderr, *) err_info%err_header_global
          write(fstderr,*) "l_rtm_nudge must be set to .false. when " &
                           // "l_implemented = .true."
          write(fstderr,*) "Fatal error in check_clubb_settings_api"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

        if ( thlm_sponge_damp_settings%l_sponge_damping ) then
          write(fstderr, *) err_info%err_header_global
          write(fstderr,*) "thlm_sponge_damp_settings%l_sponge_damping " &
                           // "must be set to .false. when  l_implemented = .true."
          write(fstderr,*) "Fatal error in check_clubb_settings_api"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

        if ( rtm_sponge_damp_settings%l_sponge_damping ) then
          write(fstderr, *) err_info%err_header_global
          write(fstderr,*) "rtm_sponge_damp_settings%l_sponge_damping " &
                           // "must be set to .false. when  l_implemented = .true."
          write(fstderr,*) "Fatal error in check_clubb_settings_api"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

        if ( uv_sponge_damp_settings%l_sponge_damping ) then
          write(fstderr, *) err_info%err_header_global
          write(fstderr,*) "uv_sponge_damp_settings%l_sponge_damping " &
                           // "must be set to .false. when  l_implemented = .true."
          write(fstderr,*) "Fatal error in check_clubb_settings_api"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

        if ( wp2_sponge_damp_settings%l_sponge_damping ) then
          write(fstderr, *) err_info%err_header_global
          write(fstderr,*) "wp2_sponge_damp_settings%l_sponge_damping " &
                           // "must be set to .false. when  l_implemented = .true."
          write(fstderr,*) "Fatal error in check_clubb_settings_api"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

        if ( wp3_sponge_damp_settings%l_sponge_damping ) then
          write(fstderr, *) err_info%err_header_global
          write(fstderr,*) "wp3_sponge_damp_settings%l_sponge_damping " &
                           // "must be set to .false. when  l_implemented = .true."
          write(fstderr,*) "Fatal error in check_clubb_settings_api"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

        if ( up2_vp2_sponge_damp_settings%l_sponge_damping ) then
          write(fstderr, *) err_info%err_header_global
          write(fstderr,*) "up2_vp2_sponge_damp_settings%l_sponge_damping " &
                           // "must be set to .false. when  l_implemented = .true."
          write(fstderr,*) "Fatal error in check_clubb_settings_api"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

      end if

      return

    end subroutine check_clubb_settings_api

end module numerical_check
