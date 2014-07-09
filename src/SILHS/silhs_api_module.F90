!--------------------------------------------------------------------------------------------------
! $Id$
!==================================================================================================
!
!       ######## ########### ###        ###    ###  ########            ###     ######### #########
!     ###    ###    ###     ###        ###    ### ###    ###         ### ###   ###    ###   ###
!    ###           ###     ###        ###    ### ###               ###   ###  ###    ###   ###
!   ##########    ###     ###        ########## ##########       ########### #########    ###
!         ###    ###     ###        ###    ###        ###       ###     ### ###          ###
! ###    ###    ###     ###        ###    ### ###    ###       ###     ### ###          ###
! ######## ########### ########## ###    ###  ########        ###     ### ###       #########
!
! The SILHS API serves as the doorway through which external models can interact with SILHS.
!
!               PLEASE REMEMBER, IF ANY CODE IS CHANGED IN THIS DOCUMENT,
!                   THE CHANGES MUST BE PROPOGATED TO ALL HOST MODELS.
!
module silhs_api_module

  implicit none

  private

  public  &
    LH_subcolumn_generator_api, &
    LH_subcolumn_generator_mod_api, &
    stats_accumulate_LH_api, &
    latin_hypercube_2D_output_api, &
    latin_hypercube_2D_close_api, &
    lh_microphys_driver_api, &
    cleanup_lh_arrays_api, &
    copy_X_nl_into_hm_all_pts_api, &
    copy_X_nl_into_rc_all_pts_api, &
    est_kessler_microphys_api

contains

  !================================================================================================
  ! LH_subcolumn_generator - Generates sample points of moisture, temperature, et cetera.
  !================================================================================================

  subroutine LH_subcolumn_generator_api( &
    iter, d_variables, num_samples, sequence_length, nz, &
    thlm, pdf_params, wm_zt, delta_zm, rcm, Ncnm, rvm, &
    hydromet, sigma2_on_mu2_ip_array_cloud, sigma2_on_mu2_ip_array_below, &
    corr_array_cloud, corr_array_below, Lscale_vert_avg, &
    X_nl_all_levs, X_mixt_comp_all_levs, lh_rt, lh_thl, &
    lh_sample_point_weights )

    use latin_hypercube_driver_module, only : LH_subcolumn_generator

    use parameters_model, only: hydromet_dim ! Variable

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      iter,            & ! Model iteration number
      d_variables,     & ! Number of variables to sample
      num_samples,     & ! Number of samples per variable
      sequence_length, & ! nt_repeat/num_samples; number of timesteps before sequence repeats.
      nz                 ! Number of vertical model levels

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params ! PDF parameters       [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      thlm,      & ! Liquid potential temperature       [K]
      wm_zt,     & ! Mean w                             [m/s]
      delta_zm     ! Difference in moment. altitudes    [m]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rcm,  & ! Liquid water mixing ratio                              [kg/kg]
      Ncnm, & ! Cloud nuclei concentration (simplified); Nc=Ncn*H(s)   [#/kg]
      rvm     ! Vapor water mixing ratio                               [kg/kg]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      sigma2_on_mu2_ip_array_cloud, & ! sigma_x^2/mu_x^2 ip; cloudy levs. [-]
      sigma2_on_mu2_ip_array_below    ! sigma_x^2/mu_x^2 ip; clear levs.  [-]

    real( kind = core_rknd ), dimension(d_variables,d_variables), intent(in) :: &
      corr_array_cloud, & ! Correlation for hydrometeor species [-]
      corr_array_below

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      Lscale_vert_avg ! 3pt vertical average of Lscale  [m]

    ! Output Variables
    real( kind = dp ), intent(out), dimension(nz,num_samples,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, intent(out), dimension(nz,num_samples) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real( kind = core_rknd ), intent(out), dimension(nz,num_samples) :: &
      lh_rt, lh_thl ! Sample of total water and liquid potential temperature [kg/kg],[K]

    real( kind = core_rknd ), intent(out), dimension(num_samples) :: &
      lh_sample_point_weights

    call LH_subcolumn_generator( &
      iter, d_variables, num_samples, sequence_length, nz, &
      thlm, pdf_params, wm_zt, delta_zm, rcm, Ncnm, rvm, &
      hydromet, sigma2_on_mu2_ip_array_cloud, sigma2_on_mu2_ip_array_below, &
      corr_array_cloud, corr_array_below, Lscale_vert_avg, &
      X_nl_all_levs, X_mixt_comp_all_levs, lh_rt, lh_thl, &
      lh_sample_point_weights )

  end subroutine LH_subcolumn_generator_api

  !================================================================================================
  ! LH_subcolumn_generator_mod - Generates sample points of moisture, temperature, et cetera.
  !================================================================================================

  subroutine LH_subcolumn_generator_mod_api( &
    iter, d_variables, num_samples, sequence_length, nz, & ! In
    pdf_params, delta_zm, rcm, Lscale_vert_avg, & ! In
    mu1, mu2, sigma1, sigma2, & ! In
    corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! In
    hydromet_pdf_params, & ! In
    X_nl_all_levs, X_mixt_comp_all_levs, lh_rt, lh_thl, & ! Out
    lh_sample_point_weights ) ! Out

    use latin_hypercube_driver_module, only : LH_subcolumn_generator_mod

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter ! Type

    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      iter,            & ! Model iteration number
      d_variables,     & ! Number of variables to sample
      num_samples,     & ! Number of samples per variable
      sequence_length, & ! nt_repeat/num_samples; number of timesteps before sequence repeats.
      nz                 ! Number of vertical model levels

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params ! PDF parameters       [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      delta_zm, &  ! Difference in moment. altitudes    [m]
      rcm          ! Liquid water mixing ratio          [kg/kg]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      Lscale_vert_avg ! 3pt vertical average of Lscale  [m]

    ! Output Variables
    real( kind = dp ), intent(out), dimension(nz,num_samples,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, intent(out), dimension(nz,num_samples) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real( kind = core_rknd ), intent(out), dimension(nz,num_samples) :: &
      lh_rt, lh_thl ! Sample of total water and liquid potential temperature [kg/kg],[K]

    real( kind = core_rknd ), intent(out), dimension(num_samples) :: &
      lh_sample_point_weights

    ! More Input Variables!
    real( kind = dp ), dimension(d_variables,d_variables,nz), intent(in) :: &
      corr_cholesky_mtx_1, & ! Correlations Cholesky matrix (1st comp.)  [-]
      corr_cholesky_mtx_2    ! Correlations Cholesky matrix (2nd comp.)  [-]

    real( kind = core_rknd ), dimension(d_variables,nz), intent(in) :: &
      mu1,    & ! Means of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>)  [units vary]
      mu2,    & ! Means of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>)  [units vary]
      sigma1, & ! Stdevs of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>) [units vary]
      sigma2    ! Stdevs of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>) [units vary]

    type(hydromet_pdf_parameter), dimension(nz), intent(in) :: &
      hydromet_pdf_params

    call LH_subcolumn_generator_mod( &
      iter, d_variables, num_samples, sequence_length, nz, & ! In
      pdf_params, delta_zm, rcm, Lscale_vert_avg, & ! In
      mu1, mu2, sigma1, sigma2, & ! In
      corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! In
      hydromet_pdf_params, & ! In
      X_nl_all_levs, X_mixt_comp_all_levs, lh_rt, lh_thl, & ! Out
      lh_sample_point_weights ) ! Out

  end subroutine LH_subcolumn_generator_mod_api

  !================================================================================================
  ! stats_accumulate_LH - Clips subcolumns from latin hypercube and creates stats.
  !================================================================================================

  subroutine stats_accumulate_LH_api( &
    nz, num_samples, d_variables, rho_ds_zt, &
    lh_sample_point_weights, X_nl_all_levs, &
    lh_thl, lh_rt, Nc_in_cloud )

    use latin_hypercube_driver_module, only : stats_accumulate_LH

    use clubb_precision, only: &
      core_rknd, & ! Variable(s)
      dp

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables,     & ! Number of variables to sample
      num_samples,   & ! Number of calls to microphysics per timestep (normally=2)
      nz                 ! Number of vertical model levels

    real( kind = core_rknd ), intent(in), dimension(nz) :: &
      rho_ds_zt  ! Dry, static density (thermo. levs.) [kg/m^3]

    real( kind = core_rknd ), intent(in), dimension(num_samples) :: &
      lh_sample_point_weights

    real( kind = dp ), intent(in), dimension(nz,num_samples,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real( kind = core_rknd ), intent(in), dimension(nz,num_samples) :: &
      lh_thl, & ! Sample of liquid potential temperature [K]
      lh_rt     ! Sample of total water mixing ratio     [kg/kg]

    real( kind = core_rknd), dimension(nz), intent(in) :: &
      ! Constant value of N_c within cloud, to be used with l_const_Nc_in_cloud
      Nc_in_cloud

    call stats_accumulate_LH( &
      nz, num_samples, d_variables, rho_ds_zt, &
      lh_sample_point_weights, X_nl_all_levs, &
      lh_thl, lh_rt, Nc_in_cloud )

  end subroutine stats_accumulate_LH_api

  !================================================================================================
  ! est_kessler_microphys - Computes microphysical grid box means of Kessler autoconversion scheme.
  !================================================================================================

  subroutine est_kessler_microphys_api( &
    nz, num_samples, d_variables, &
    X_nl_all_levs, pdf_params, rcm, cloud_frac, &
    X_mixt_comp_all_levs, lh_sample_point_weights, &
    lh_AKm, AKm, AKstd, AKstd_cld, &
    AKm_rcm, AKm_rcc, lh_rcm_avg )

    use est_kessler_microphys_module, only : est_kessler_microphys

    use pdf_parameter_module, only:  &
      pdf_parameter  ! Type

    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd

    implicit none

    ! Input Variables

    integer, intent(in) :: &
      nz, &          ! Number of vertical levels
      num_samples, & ! Number of sample points
      d_variables    ! Number of variates

    real( kind = dp ), dimension(nz,num_samples,d_variables), intent(in) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      cloud_frac    ! Cloud fraction           [-]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rcm          ! Liquid water mixing ratio                [kg/kg]

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params ! PDF parameters       [units vary]

    integer, dimension(nz,num_samples), intent(in) :: &
      X_mixt_comp_all_levs ! Whether we're in mixture component 1 or 2

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      lh_sample_point_weights ! Weight for cloud weighted sampling

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      lh_AKm,    & ! Monte Carlo estimate of Kessler autoconversion [kg/kg/s]
      AKm,       & ! Exact Kessler autoconversion, AKm,             [kg/kg/s]
      AKstd,     & ! Exact standard deviation of gba Kessler        [kg/kg/s]
      AKstd_cld, & ! Exact w/in cloud std of gba Kessler            [kg/kg/s]
      AKm_rcm,   & ! Exact local gba Kessler auto based on rcm      [kg/kg/s]
      AKm_rcc      ! Exact local gba Kessler based on w/in cloud rc [kg/kg/s]

    ! For comparison, estimate kth liquid water using Monte Carlo
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      lh_rcm_avg ! LH estimate of grid box avg liquid water [kg/kg]

    call est_kessler_microphys( &
      nz, num_samples, d_variables, &
      X_nl_all_levs, pdf_params, rcm, cloud_frac, &
      X_mixt_comp_all_levs, lh_sample_point_weights, &
      lh_AKm, AKm, AKstd, AKstd_cld, &
      AKm_rcm, AKm_rcc, lh_rcm_avg )

  end subroutine est_kessler_microphys_api

  !================================================================================================
  ! The functions and subroutines below are only used by CLUBB_standalone
  !================================================================================================
  !================================================================================================
  ! latin_hypercube_2D_output
  !================================================================================================

  subroutine latin_hypercube_2D_output_api( &
    fname_prefix, fdir, stats_tout, nz, &
    zt, time_initial )

    use latin_hypercube_driver_module, only : latin_hypercube_2D_output

    use clubb_precision, only: &
      time_precision, & ! Constant
      core_rknd

    implicit none

    ! Input Variables
    character(len=*), intent(in) :: &
      fname_prefix, & ! Prefix for file name
      fdir            ! Directory for output

    real(kind=time_precision), intent(in) :: &
      stats_tout, & ! Frequency to write to disk        [s]
      time_initial  ! Initial time                      [s]

    integer, intent(in) :: &
      nz ! Number of vertical levels

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      zt ! Altitudes [m]

    call latin_hypercube_2D_output( &
      fname_prefix, fdir, stats_tout, nz, &
      zt, time_initial )

  end subroutine latin_hypercube_2D_output_api

  !================================================================================================
  ! latin_hypercube_2D_close - Closes a 2D sample file.
  !================================================================================================

  subroutine latin_hypercube_2D_close_api

    use latin_hypercube_driver_module, only : latin_hypercube_2D_close

    implicit none

    call latin_hypercube_2D_close

  end subroutine latin_hypercube_2D_close_api

  !================================================================================================
  ! lh_microphys_driver - Computes an estimate of the change due to microphysics.
  !================================================================================================

  subroutine lh_microphys_driver_api( &
    dt, nz, num_samples, d_variables, &
    X_nl_all_levs, lh_rt, lh_thl, lh_sample_point_weights, &
    pdf_params, p_in_Pa, exner, rho, &
    rcm, delta_zt, cloud_frac, &
    hydromet, X_mixt_comp_all_levs, Nc_in_cloud, &
    lh_hydromet_mc, lh_hydromet_vel, lh_Ncm_mc, &
    lh_rcm_mc, lh_rvm_mc, lh_thlm_mc, &
    lh_rtp2_mc, lh_thlp2_mc, lh_wprtp_mc, &
    lh_wpthlp_mc, lh_rtpthlp_mc, &
    microphys_sub )

    use latin_hypercube_driver_module, only : lh_microphys_driver

    use parameters_model, only: hydromet_dim ! Variable

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type


    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd, &
      time_precision

    implicit none

    ! Interface block
#include "microphys_interface.inc"

    ! Input Variables
    real( kind = time_precision ), intent(in) :: &
      dt ! Model timestep       [s]

    integer, intent(in) :: &
      d_variables,     & ! Number of variables to sample
      num_samples,   & ! Number of calls to microphysics per timestep (normally=2)
      nz               ! Number of vertical model levels

    ! Input Variables
    real( kind = dp ), intent(in), dimension(nz,num_samples,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, intent(in), dimension(nz,num_samples) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real( kind = core_rknd ), intent(in), dimension(nz,num_samples) :: &
      lh_rt, lh_thl ! Sample of total water and liquid potential temperature [kg/kg],[K]

    real( kind = core_rknd ), intent(in), dimension(num_samples) :: &
      lh_sample_point_weights ! Weight given the individual sample points

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params ! PDF parameters       [units vary]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      cloud_frac,  & ! Cloud fraction               [-]
      delta_zt,    & ! Change in meters with height [m]
      rcm,         & ! Liquid water mixing ratio    [kg/kg]
      p_in_Pa,     & ! Pressure                     [Pa]
      exner,       & ! Exner function               [-]
      rho            ! Density on thermo. grid      [kg/m^3]

    real( kind = core_rknd), dimension(nz), intent(in) :: &
      ! Constant value of N_c within cloud, to be used with l_const_Nc_in_cloud
      Nc_in_cloud

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      lh_hydromet_mc, & ! LH estimate of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_vel   ! LH estimate of hydrometeor sedimentation velocity [m/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      lh_Ncm_mc,     & ! LH estimate of time tndcy. of cloud droplet conc.     [num/kg/s]
      lh_rcm_mc,     & ! LH estimate of time tndcy. of liq. water mixing ratio [kg/kg/s]
      lh_rvm_mc,     & ! LH estimate of time tndcy. of vapor water mix. ratio  [kg/kg/s]
      lh_thlm_mc,    & ! LH estimate of time tndcy. of liquid potential temp.  [K/s]
      lh_rtp2_mc,    & ! LH microphysics tendency for <rt'^2>                  [(kg/kg)^2/s]
      lh_thlp2_mc,   & ! LH microphysics tendency for <thl'^2>                 [K^2/s]
      lh_wprtp_mc,   & ! LH microphysics tendency for <w'rt'>                  [m*(kg/kg)/s^2]
      lh_wpthlp_mc,  & ! LH microphysics tendency for <w'thl'>                 [m*K/s^2]
      lh_rtpthlp_mc    ! LH microphysics tendency for <rt'thl'>                [K*(kg/kg)/s]

    call lh_microphys_driver( &
      dt, nz, num_samples, d_variables, &
      X_nl_all_levs, lh_rt, lh_thl, lh_sample_point_weights, &
      pdf_params, p_in_Pa, exner, rho, &
      rcm, delta_zt, cloud_frac, &
      hydromet, X_mixt_comp_all_levs, Nc_in_cloud, &
      lh_hydromet_mc, lh_hydromet_vel, lh_Ncm_mc, &
      lh_rcm_mc, lh_rvm_mc, lh_thlm_mc, &
      lh_rtp2_mc, lh_thlp2_mc, lh_wprtp_mc, &
      lh_wpthlp_mc, lh_rtpthlp_mc, &
      microphys_sub )

  end subroutine lh_microphys_driver_api

  !================================================================================================
  ! cleanup_latin_hypercube_arrays - De-allocate latin hypercube arrays.
  !================================================================================================

  subroutine cleanup_lh_arrays_api

    use latin_hypercube_arrays, only : cleanup_latin_hypercube_arrays

    implicit none

    call cleanup_latin_hypercube_arrays

  end subroutine cleanup_lh_arrays_api

  !================================================================================================
  ! copy_X_nl_into_hydromet_all_pts - Copies latin hypercube sample points to a hydrometeor array.
  !================================================================================================

  subroutine copy_X_nl_into_hm_all_pts_api( &
    nz, d_variables, num_samples, &
    X_nl_all_levs, &
    hydromet, &
    hydromet_all_points, &
    Ncn_all_points )

    use estimate_scm_microphys_module, only : copy_X_nl_into_hydromet_all_pts

    use parameters_model, only: &
      hydromet_dim ! Variable

    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd

    implicit none

    integer, intent(in) :: &
      nz,            & ! Number of vertical levels
      d_variables,   & ! Number of variates
      num_samples    ! Number of calls to microphysics

    real( kind = dp ), dimension(nz,num_samples,d_variables), intent(in) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nz,num_samples,hydromet_dim), intent(out) :: &
      hydromet_all_points ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nz,num_samples), intent(out) :: &
      Ncn_all_points    ! Cloud nuclei conc. (simplified); Nc=Ncn*H(s)   [#/kg]

    call copy_X_nl_into_hydromet_all_pts( &
      nz, d_variables, num_samples, &
      X_nl_all_levs, &
      hydromet, &
      hydromet_all_points, &
      Ncn_all_points )

  end subroutine copy_X_nl_into_hm_all_pts_api

  !================================================================================================
  ! copy_X_nl_into_rc_all_pts - Extracts a sample of rc from X_nl_all_levs for each subcolumn.
  !================================================================================================

  subroutine copy_X_nl_into_rc_all_pts_api( &
    nz, d_variables, num_samples, &
    X_nl_all_levs, lh_rc )

    use estimate_scm_microphys_module, only : copy_X_nl_into_rc_all_pts

    use clubb_precision, only: &
      dp, &
      core_rknd

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz, &            ! Number of vertical levels
      d_variables, &   ! Number of lognormal variates
      num_samples      ! Number of SILHS samples per variate

    real( kind = dp ), dimension(nz,num_samples,d_variables), intent(in) :: &
      X_nl_all_levs    ! Normal-lognormal SILHS sample   [units vary]

    ! Output variables
    real( kind = core_rknd ), dimension(nz,num_samples), intent(out) :: &
      lh_rc            ! SILHS samples of rc       [kg/kg]

    call copy_X_nl_into_rc_all_pts( &
      nz, d_variables, num_samples, &
    X_nl_all_levs, lh_rc )

  end subroutine copy_X_nl_into_rc_all_pts_api

end module silhs_api_module
