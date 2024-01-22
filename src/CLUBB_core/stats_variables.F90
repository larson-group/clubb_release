!-------------------------------------------------------------------------------
! $Id$
!-------------------------------------------------------------------------------

! Description:
!   Holds pointers and other variables for statistics to be written to 
!   GrADS files and netCDF files.
!-------------------------------------------------------------------------------
module stats_variables

  use clubb_precision, only:  & 
      core_rknd ! Variable(s)

  implicit none

  private ! Set Default Scope

  !=============================

  public :: stats_metadata_type
  
  type stats_metadata_type

    ! Sampling and output frequencies
    real( kind = core_rknd ) :: &
      stats_tsamp = 0._core_rknd, & ! Sampling interval   [s]
      stats_tout  = 0._core_rknd ! Output interval     [s]

    logical ::  &
      l_stats            = .false., & ! Main flag to turn statistics on/off
      l_output_rad_files = .false., & ! Flag to turn off radiation statistics output
      l_netcdf           = .false., & ! Output to NetCDF format
      l_grads            = .false., & ! Output to GrADS format
      l_silhs_out        = .false., & ! Output SILHS files (stats_lh_zt and stats_lh_sfc)
      l_allow_small_stats_tout = .false. ! Do not stop if output timestep is too low for
                        ! requested format, e.g. l_grads = .true. and
                        ! stats_tout < 60.0

    logical :: & 
      l_stats_samp = .false., & ! Sample flag for current time step
      l_stats_last = .false.    ! Last time step of output period

    character(len=200) ::  & 
      fname_zt     = '', & ! Name of the stats file for thermodynamic grid fields
      fname_lh_zt  = '', & ! Name of the stats file for LH variables on the stats_zt grid
      fname_lh_sfc = '', & ! Name of the stats file for LH variables on the stats_zt grid
      fname_zm     = '', & ! Name of the stats file for momentum grid fields
      fname_rad_zt = '', & ! Name of the stats file for the stats_zt radiation grid fields
      fname_rad_zm = '', & ! Name of the stats file for the stats_zm radiation grid fields
      fname_sfc    = ''    ! Name of the stats file for surface only fields


    !====================================================================
    !                           ZT Indices
    !====================================================================
    integer :: & 
       ithlm = 0, & 
       ithvm = 0, & 
       irtm = 0, & 
       ircm = 0, &
       irvm = 0, & 
       ium = 0, & 
       ivm = 0, & 
       iwm_zt = 0, &
       iwm_zm = 0, &
       ium_ref = 0,&
       ivm_ref = 0, & 
       iug = 0, & 
       ivg = 0, & 
       icloud_frac = 0, &
       iice_supersat_frac = 0, &
       ircm_in_layer = 0, &
       ircm_in_cloud = 0, &
       icloud_cover = 0, &
       ip_in_Pa = 0, & 
       iexner = 0, & 
       irho_ds_zt = 0, &
       ithv_ds_zt = 0, &
       iLscale = 0, & 
       iwp3 = 0, & 
       ithlp3 = 0, &
       irtp3 = 0, &
       iwpthlp2 = 0, & 
       iwp2thlp = 0, & 
       iwprtp2 = 0, & 
       iwp2rtp = 0, &
       iSkw_zt = 0, &
       iSkthl_zt = 0, &
       iSkrt_zt = 0, &
       ircm_supersat_adj = 0

    integer :: & 
       iLscale_up = 0, & 
       iLscale_down = 0, & 
       iLscale_pert_1 = 0, & 
       iLscale_pert_2 = 0, & 
       itau_zt = 0, &
       iinvrs_tau_zt = 0, & 
       iKh_zt = 0, & 
       iwp2thvp = 0, & 
       iwp2rcp = 0, & 
       iw_up_in_cloud = 0, &
       iw_down_in_cloud = 0, &
       icld_updr_frac = 0,   &
       icld_downdr_frac = 0, &
       iwprtpthlp = 0, &
       irc_coef = 0, &
       isigma_sqd_w_zt = 0, & 
       irho = 0

    integer :: &
       iinvrs_tau_zm       = 0, &
       iinvrs_tau_xp2_zm   = 0, &
       iinvrs_tau_wp2_zm   = 0, &
       iinvrs_tau_wpxp_zm  = 0, &
       iinvrs_tau_wp3_zm   = 0, &
       iinvrs_tau_no_N2_zm = 0, &
       iinvrs_tau_bkgnd = 0,  &
       iinvrs_tau_sfc   = 0,  &
       iinvrs_tau_shear = 0,  &
       itau_no_N2_zm = 0,     & 
       itau_xp2_zm   = 0,     &
       itau_wp2_zm   = 0,     &
       itau_wp3_zm   = 0,     &
       itau_wpxp_zm  = 0


    integer, dimension(:), allocatable :: & 
       ihm_1, &
       ihm_2

    integer :: & 
       iprecip_frac = 0, &
       iprecip_frac_1 = 0, &
       iprecip_frac_2 = 0, &
       iNcnm = 0 

    integer, dimension(:), allocatable :: &
       imu_hm_1,         &
       imu_hm_2,         &
       imu_hm_1_n,       &
       imu_hm_2_n,       &
       isigma_hm_1,      &
       isigma_hm_2,      &
       isigma_hm_1_n,    &
       isigma_hm_2_n,    &
       icorr_w_hm_1,     &
       icorr_w_hm_2,     &
       icorr_chi_hm_1,   &
       icorr_chi_hm_2,   &
       icorr_eta_hm_1,   &
       icorr_eta_hm_2,   &
       icorr_Ncn_hm_1,   &
       icorr_Ncn_hm_2,   &
       icorr_w_hm_1_n,   &
       icorr_w_hm_2_n,   &
       icorr_chi_hm_1_n, &
       icorr_chi_hm_2_n, &
       icorr_eta_hm_1_n, &
       icorr_eta_hm_2_n, &
       icorr_Ncn_hm_1_n, &
       icorr_Ncn_hm_2_n

    integer, dimension(:,:), allocatable :: &
       icorr_hmx_hmy_1,   &
       icorr_hmx_hmy_2,   &
       icorr_hmx_hmy_1_n, &
       icorr_hmx_hmy_2_n

    integer :: &
       imu_Ncn_1 = 0,      &
       imu_Ncn_2 = 0,      &
       imu_Ncn_1_n = 0,    &
       imu_Ncn_2_n = 0,    &
       isigma_Ncn_1 = 0,   &
       isigma_Ncn_2 = 0,   &
       isigma_Ncn_1_n = 0, &
       isigma_Ncn_2_n = 0

    integer :: &
       icorr_w_chi_1_ca = 0, &
       icorr_w_chi_2_ca = 0, &
       icorr_w_eta_1_ca = 0, &
       icorr_w_eta_2_ca = 0, &
       icorr_w_Ncn_1 = 0,  &
       icorr_w_Ncn_2 = 0,  &
       icorr_chi_eta_1_ca = 0, &
       icorr_chi_eta_2_ca = 0, &
       icorr_chi_Ncn_1 = 0,  &
       icorr_chi_Ncn_2 = 0,  &
       icorr_eta_Ncn_1 = 0,  &
       icorr_eta_Ncn_2 = 0

    integer :: &
       icorr_w_Ncn_1_n = 0, &
       icorr_w_Ncn_2_n = 0, &
       icorr_chi_Ncn_1_n = 0, &
       icorr_chi_Ncn_2_n = 0, &
       icorr_eta_Ncn_1_n = 0, &
       icorr_eta_Ncn_2_n = 0

    integer, dimension(:), allocatable :: &
      isilhs_variance_category, &
      ilh_samp_frac_category


    integer :: & 
       iNcm = 0,             & ! Brian
       iNccnm = 0,           & 
       iNc_in_cloud = 0,     &
       iNc_activated = 0,    &
       isnowslope = 0,       & ! Adam Smith, 22 April 2008
       ised_rcm = 0,         & ! Brian
       irsat = 0,            & ! Brian
       irsati = 0,           & 
       irrm = 0,          & ! Brian
       im_vol_rad_rain = 0,  & ! Brian
       im_vol_rad_cloud = 0, & ! COAMPS only. dschanen 6 Dec 2006
       iprecip_rate_zt = 0,    & ! Brian
       iAKm = 0,             & ! analytic Kessler.  Vince Larson 22 May 2005 
       ilh_AKm = 0,          & ! LH Kessler.  Vince Larson  22 May 2005
       iradht = 0,           & ! Radiative heating.
       iradht_LW = 0,        & !   "           "   Long-wave component
       iradht_SW = 0,        & !   "           "   Short-wave component
       irel_humidity = 0

    integer :: & 
       iAKstd = 0,     &
       iAKstd_cld = 0, & 
       iAKm_rcm = 0, & 
       iAKm_rcc = 0


    integer :: & 
     irfrzm = 0

    ! Skewness functions on stats_zt grid
    integer :: &
      iC11_Skw_fnc = 0


    integer :: &
      icloud_frac_zm = 0, &
      iice_supersat_frac_zm = 0, &
      ircm_zm = 0, &
      irtm_zm = 0, &
      ithlm_zm = 0


    integer :: &
      iw_1_zm = 0, &
      iw_2_zm = 0, &
      ivarnce_w_1_zm = 0, &
      ivarnce_w_2_zm = 0, &
      imixt_frac_zm = 0


    integer :: &
      ilh_rcm_avg = 0, &
      ik_lh_start = 0


    integer :: & 
       iNrm = 0,       & ! Rain droplet number concentration
       iNim = 0,       & ! Ice number concentration
       iNsm = 0,    & ! Snow number concentration
       iNgm = 0    ! Graupel number concentration

    integer :: & 
       iT_in_K = 0      ! Absolute temperature

    integer :: &
      ieff_rad_cloud = 0, &
      ieff_rad_ice = 0, &
      ieff_rad_snow = 0, &
      ieff_rad_rain = 0, &
      ieff_rad_graupel = 0


    integer :: & 
      irsm = 0, &
      irgm = 0, & 
      irim = 0, & 
      idiam = 0,           & ! Diameter of ice crystal           [m]
      imass_ice_cryst = 0, & ! Mass of a single ice crystal      [kg]
      ircm_icedfs = 0,     & ! Change in liquid water due to ice [kg/kg/s]
      iu_T_cm = 0         ! Fallspeed of ice crystal in cm/s  [cm s^{-1}]



    ! thlm/rtm budget terms
    integer :: & 
      irtm_bt = 0,         & ! rtm total time tendency
      irtm_ma = 0,         & ! rtm mean advect. term
      irtm_ta = 0,         & ! rtm turb. advect. term
      irtm_forcing = 0,    & ! rtm large scale forcing term
      irtm_mc = 0,         & ! rtm change from microphysics
      irtm_sdmp = 0,       & ! rtm change from sponge damping
      irvm_mc = 0,         & ! rvm change from microphysics
      ircm_mc = 0,         & ! rcm change from microphysics
      ircm_sd_mg_morr = 0, & ! rcm sedimentation tendency
      irtm_mfl = 0,        & ! rtm change due to monotonic flux limiter
      irtm_tacl = 0,       & ! rtm correction from turbulent advection (wprtp) clipping
      irtm_cl = 0,         & ! rtm clipping term
      irtm_pd = 0,         & ! thlm postive definite adj term
      ithlm_bt = 0,        & ! thlm total time tendency
      ithlm_ma = 0,        & ! thlm mean advect. term
      ithlm_ta = 0,        & ! thlm turb. advect. term
      ithlm_forcing = 0,   & ! thlm large scale forcing term
      ithlm_sdmp = 0,      & ! thlm change from sponge damping
      ithlm_mc = 0,        & ! thlm change from microphysics
      ithlm_mfl = 0,       & ! thlm change due to monotonic flux limiter
      ithlm_tacl = 0,      & ! thlm correction from turbulent advection (wpthlp) clipping
      ithlm_cl = 0           ! thlm clipping term


    !monatonic flux limiter diagnostic terms
    integer :: &
      ithlm_mfl_min = 0, &
      ithlm_mfl_max = 0, &
      iwpthlp_enter_mfl = 0, &
      iwpthlp_exit_mfl = 0, &
      iwpthlp_mfl_min = 0, &
      iwpthlp_mfl_max = 0, &
      irtm_mfl_min = 0, &
      irtm_mfl_max = 0, &
      iwprtp_enter_mfl = 0, &
      iwprtp_exit_mfl = 0, &
      iwprtp_mfl_min = 0, &
      iwprtp_mfl_max = 0, &
      ithlm_enter_mfl = 0, &
      ithlm_exit_mfl = 0, &
      ithlm_old = 0, &
      ithlm_without_ta = 0, &
      irtm_enter_mfl = 0, &
      irtm_exit_mfl = 0, &
      irtm_old = 0, &
      irtm_without_ta = 0


    integer :: & 
       iwp3_bt  = 0, & 
       iwp3_ma  = 0, & 
       iwp3_ta  = 0, & 
       iwp3_tp  = 0, & 
       iwp3_ac  = 0, & 
       iwp3_bp1 = 0, &
       iwp3_pr_tp = 0, & 
       iwp3_pr_turb = 0, &
       iwp3_pr_dfsn = 0, & 
       iwp3_pr1 = 0, & 
       iwp3_pr2 = 0, & 
       iwp3_pr3 = 0, &
       iwp3_dp1 = 0, &
       iwp3_sdmp = 0, &
       iwp3_cl  = 0, &
       iwp3_splat = 0


    integer :: &
      irtp3_bt  = 0, &
      irtp3_tp  = 0, &
      irtp3_ac  = 0, &
      irtp3_dp  = 0, &
      ithlp3_bt = 0, &
      ithlp3_tp = 0, &
      ithlp3_ac = 0, &
      ithlp3_dp = 0


    ! Rain mixing ratio budgets
    integer :: & 
       irrm_bt = 0, &
       irrm_ma = 0, &
       irrm_ta = 0, &
       irrm_sd = 0, &
       irrm_ts = 0, &
       irrm_sd_morr = 0, &
       irrm_evap = 0, &
       irrm_auto = 0, &
       irrm_accr = 0, &
       irrm_evap_adj = 0, &
       irrm_src_adj = 0, &
       irrm_mc_nonadj = 0, &
       irrm_mc = 0, &
       irrm_hf = 0, &
       irrm_wvhf = 0, &
       irrm_cl = 0


    integer :: &
       iNrm_bt = 0, &
       iNrm_ma = 0, &
       iNrm_ta = 0, &
       iNrm_sd = 0, &
       iNrm_ts = 0, &
       iNrm_evap = 0, &
       iNrm_auto = 0, &
       iNrm_evap_adj = 0, &
       iNrm_src_adj = 0, &
       iNrm_mc = 0, &
       iNrm_cl = 0



    ! Snow/Ice/Graupel mixing ratio budgets
    integer :: & 
       irsm_bt = 0, & 
       irsm_ma = 0, & 
       irsm_sd = 0, & 
       irsm_sd_morr = 0, &
       irsm_ta = 0, & 
       irsm_mc = 0, & 
       irsm_hf = 0, &
       irsm_wvhf = 0, &
       irsm_cl = 0, &
       irsm_sd_morr_int = 0


    integer :: & 
       irgm_bt = 0, & 
       irgm_ma = 0, & 
       irgm_sd = 0, & 
       irgm_sd_morr = 0, &
       irgm_ta = 0, & 
       irgm_mc = 0, & 
       irgm_hf = 0, &
       irgm_wvhf = 0, &
       irgm_cl = 0


    integer :: & 
       irim_bt = 0, & 
       irim_ma = 0, & 
       irim_sd = 0, & 
       irim_sd_mg_morr = 0, &
       irim_ta = 0, & 
       irim_mc = 0, & 
       irim_hf = 0, &
       irim_wvhf = 0, &
       irim_cl = 0


    integer :: &
      iNsm_bt = 0,  &
      iNsm_ma = 0,  &
      iNsm_sd = 0,  &
      iNsm_ta = 0,  &
      iNsm_mc = 0,  &
      iNsm_cl = 0


    integer :: &
      iNgm_bt = 0, &
      iNgm_ma = 0, &
      iNgm_sd = 0, &
      iNgm_ta = 0, &
      iNgm_mc = 0, &
      iNgm_cl = 0


    integer :: &
      iNim_bt = 0, &
      iNim_ma = 0, &
      iNim_sd = 0, &
      iNim_ta = 0, &
      iNim_mc = 0, &
      iNim_cl = 0


    integer :: &
      iNcm_bt = 0, &
      iNcm_ma = 0, &
      iNcm_ta = 0, &
      iNcm_mc = 0, &
      iNcm_cl = 0, &
      iNcm_act = 0


    ! Covariances between w, r_t, theta_l and KK microphysics tendencies.
    ! Additionally, covariances between r_r and N_r and KK rain drop mean
    ! volume radius.  These are all calculated on thermodynamic grid levels.
    integer :: &
      iw_KK_evap_covar_zt = 0,   & ! Covariance of w and KK evaporation tendency.
      irt_KK_evap_covar_zt = 0,  & ! Covariance of r_t and KK evaporation tendency.
      ithl_KK_evap_covar_zt = 0, & ! Covariance of theta_l and KK evap. tendency.
      iw_KK_auto_covar_zt = 0,   & ! Covariance of w and KK autoconversion tendency.
      irt_KK_auto_covar_zt = 0,  & ! Covariance of r_t and KK autoconversion tendency.
      ithl_KK_auto_covar_zt = 0, & ! Covariance of theta_l and KK autoconv. tendency.
      iw_KK_accr_covar_zt = 0,   & ! Covariance of w and KK accretion tendency.
      irt_KK_accr_covar_zt = 0,  & ! Covariance of r_t and KK accretion tendency.
      ithl_KK_accr_covar_zt = 0, & ! Covariance of theta_l and KK accretion tendency.
      irr_KK_mvr_covar_zt = 0,   & ! Covariance of r_r and KK mean volume radius.
      iNr_KK_mvr_covar_zt = 0,   & ! Covariance of N_r and KK mean volume radius.
      iKK_mvr_variance_zt = 0      ! Variance of KK rain drop mean volume radius.


    ! Wind budgets
    integer :: & 
       ivm_bt = 0, & 
       ivm_ma = 0, & 
       ivm_ta = 0, & 
       ivm_gf = 0, & 
       ivm_cf = 0, &
       ivm_f = 0, &
       ivm_sdmp = 0, &
       ivm_ndg = 0, &
       ivm_mfl = 0


    integer :: & 
       ium_bt = 0, & 
       ium_ma = 0, & 
       ium_ta = 0, & 
       ium_gf = 0, & 
       ium_cf = 0, & 
       ium_f = 0, &
       ium_sdmp = 0, &
       ium_ndg = 0, &
       ium_mfl = 0



    ! PDF parameters
    integer :: & 
       imixt_frac = 0, & 
       iw_1 = 0, & 
       iw_2 = 0, & 
       ivarnce_w_1 = 0, & 
       ivarnce_w_2 = 0, & 
       ithl_1 = 0, & 
       ithl_2 = 0, & 
       ivarnce_thl_1 = 0, & 
       ivarnce_thl_2 = 0, & 
       irt_1 = 0, & 
       irt_2 = 0, & 
       ivarnce_rt_1 = 0, & 
       ivarnce_rt_2 = 0, & 
       irc_1 = 0, & 
       irc_2 = 0, & 
       irsatl_1 = 0, & 
       irsatl_2 = 0, & 
       icloud_frac_1 = 0, & 
       icloud_frac_2 = 0

    integer :: & 
       ichi_1 = 0, &
       ichi_2 = 0, &
       istdev_chi_1 = 0, & 
       istdev_chi_2 = 0, &
       ichip2 = 0, &
       istdev_eta_1 = 0, &
       istdev_eta_2 = 0, &
       icovar_chi_eta_1 = 0, &
       icovar_chi_eta_2 = 0, &
       icorr_w_chi_1 = 0, &
       icorr_w_chi_2 = 0, &
       icorr_w_eta_1 = 0, &
       icorr_w_eta_2 = 0, &
       icorr_chi_eta_1 = 0, &
       icorr_chi_eta_2 = 0, &
       icorr_w_rt_1 = 0, &
       icorr_w_rt_2 = 0, &
       icorr_w_thl_1 = 0, &
       icorr_w_thl_2 = 0, &
       icorr_rt_thl_1 = 0, &
       icorr_rt_thl_2 = 0, &
       icrt_1 = 0, &
       icrt_2 = 0, &
       icthl_1 = 0, &
       icthl_2 = 0

    integer :: &
      iF_w = 0, &
      iF_rt = 0, &
      iF_thl = 0, &
      imin_F_w = 0, &
      imax_F_w = 0, &
      imin_F_rt = 0, &
      imax_F_rt = 0, &
      imin_F_thl = 0, &
      imax_F_thl = 0


    integer :: &
      icoef_wprtp2_implicit = 0, &
      iterm_wprtp2_explicit = 0, &
      icoef_wpthlp2_implicit = 0, &
      iterm_wpthlp2_explicit = 0, &
      icoef_wprtpthlp_implicit = 0, &
      iterm_wprtpthlp_explicit = 0, &
      icoef_wp2rtp_implicit = 0, &
      iterm_wp2rtp_explicit = 0, &
      icoef_wp2thlp_implicit = 0, &
      iterm_wp2thlp_explicit = 0


    integer :: & 
      iwp2_zt = 0, & 
      ithlp2_zt = 0, & 
      iwpthlp_zt = 0, & 
      iwprtp_zt = 0, & 
      irtp2_zt = 0, & 
      irtpthlp_zt = 0, &
      iup2_zt = 0, &
      ivp2_zt = 0, &
      iupwp_zt = 0, &
      ivpwp_zt = 0


    integer, dimension(:), allocatable :: &
      iwp2hmp


    integer, dimension(:), allocatable :: &
      ihydrometp2,  &
      iwphydrometp, &
      irtphmp,      &
      ithlphmp


    integer, dimension(:,:), allocatable :: &
      ihmxphmyp


    integer, dimension(:), allocatable :: &
      ihmp2_zt


    integer :: &
      ichi = 0

    integer, allocatable, dimension(:) :: & 
      isclrm,   & ! Passive scalar mean (1)
      isclrm_f    ! Passive scalar forcing (1)

  ! Used to calculate clear-sky radiative fluxes.
    integer :: &
      ifulwcl = 0, ifdlwcl = 0, ifdswcl = 0, ifuswcl = 0


    integer, allocatable, dimension(:) :: & 
      iedsclrm,   & ! Eddy-diff. scalar term (1)
      iedsclrm_f    ! Eddy-diffusivity scalar forcing (1)


    integer :: &
      ilh_thlm_mc = 0,      & ! Latin hypercube estimate of thlm_mc
      ilh_rvm_mc = 0,       & ! Latin hypercube estimate of rvm_mc
      ilh_rcm_mc = 0,       & ! Latin hypercube estimate of rcm_mc
      ilh_Ncm_mc = 0,       & ! Latin hypercube estimate of Ncm_mc
      ilh_rrm_mc = 0,    & ! Latin hypercube estimate of rrm_mc
      ilh_Nrm_mc = 0,       & ! Latin hypercube estimate of Nrm_mc
      ilh_rsm_mc = 0,    & ! Latin hypercube estimate of rsm_mc
      ilh_Nsm_mc = 0,    & ! Latin hypercube estimate of Nsm_mc
      ilh_rgm_mc = 0, & ! Latin hypercube estimate of rgm_mc
      ilh_Ngm_mc = 0, & ! Latin hypercube estimate of Ngm_mc
      ilh_rim_mc = 0,     & ! Latin hypercube estimate of rim_mc
      ilh_Nim_mc = 0          ! Latin hypercube estimate of Nim_mc

    integer :: &
      ilh_rrm_auto = 0, & ! Latin hypercube estimate of autoconversion
      ilh_rrm_accr = 0, & ! Latin hypercube estimate of accretion
      ilh_rrm_evap = 0, & ! Latin hypercube estimate of evaporation
      ilh_Nrm_auto    = 0, & ! Latin hypercube estimate of Nrm autoconversion
      ilh_Nrm_evap    = 0, & ! Latin hypercube estimate of Nrm evaporation
      ilh_m_vol_rad_rain = 0, &
      ilh_rrm_mc_nonadj = 0


    integer :: &
      ilh_rrm_src_adj  = 0, & ! Latin hypercube estimate of source adjustment (KK only!)
      ilh_rrm_evap_adj = 0, & ! Latin hypercube estimate of evap adjustment (KK only!)
      ilh_Nrm_src_adj     = 0, & ! Latin hypercube estimate of Nrm source adjustmet (KK only!)
      ilh_Nrm_evap_adj    = 0    ! Latin hypercube estimate of Nrm evap adjustment (KK only!)

    integer :: &
      ilh_Vrr = 0, & ! Latin hypercube estimate of rrm sedimentation velocity
      ilh_VNr = 0    ! Latin hypercube estimate of Nrm sedimentation velocity

    integer :: &
      ilh_rrm = 0, &
      ilh_Nrm = 0, &
      ilh_rim = 0, &
      ilh_Nim = 0, &
      ilh_rsm = 0, &
      ilh_Nsm = 0, &
      ilh_rgm = 0, &
      ilh_Ngm = 0, &
      ilh_thlm = 0, &
      ilh_rcm = 0, &
      ilh_Ncm = 0, &
      ilh_Ncnm = 0, &
      ilh_rvm = 0, &
      ilh_wm = 0, &
      ilh_cloud_frac = 0, &
      ilh_chi = 0, &
      ilh_eta = 0, &
      ilh_precip_frac = 0, &
      ilh_mixt_frac = 0


    integer :: &
      ilh_cloud_frac_unweighted  = 0,  &
      ilh_precip_frac_unweighted = 0,  &
      ilh_mixt_frac_unweighted   = 0


    integer :: &
      ilh_wp2_zt = 0, &
      ilh_Nrp2_zt = 0, &
      ilh_Ncnp2_zt = 0, &
      ilh_Ncp2_zt = 0, &
      ilh_rcp2_zt = 0, &
      ilh_rtp2_zt = 0, &
      ilh_thlp2_zt = 0, &
      ilh_rrp2_zt = 0, &
      ilh_chip2 = 0 ! Eric Raut

    ! SILHS covariance estimate indicies
    integer :: &
      ilh_rtp2_mc    = 0, &
      ilh_thlp2_mc   = 0, &
      ilh_wprtp_mc   = 0, &
      ilh_wpthlp_mc  = 0, &
      ilh_rtpthlp_mc = 0


    ! Indices for Morrison budgets
    integer :: &
      iPSMLT = 0, &
      iEVPMS = 0, &
      iPRACS = 0, &
      iEVPMG = 0, &
      iPRACG = 0, &
      iPGMLT = 0, &
      iMNUCCC = 0, &
      iPSACWS = 0, &
      iPSACWI = 0, &
      iQMULTS = 0, &
      iQMULTG = 0, &
      iPSACWG = 0, &
      iPGSACW = 0, &
      iPRD = 0, &
      iPRCI = 0, &
      iPRAI = 0, &
      iQMULTR = 0, &
      iQMULTRG = 0, &
      iMNUCCD = 0, &
      iPRACI = 0, &
      iPRACIS = 0, &
      iEPRD = 0, &
      iMNUCCR = 0, &
      iPIACR = 0, &
      iPIACRS = 0, &
      iPGRACS = 0, &
      iPRDS = 0, &
      iEPRDS = 0, &
      iPSACR = 0, &
      iPRDG = 0, &
      iEPRDG = 0


    ! More indices for Morrison budgets!!
    integer :: &
      iNGSTEN = 0, &
      iNRSTEN = 0, &
      iNISTEN = 0, &
      iNSSTEN = 0, &
      iNCSTEN = 0, &
      iNPRC1 = 0,  &
      iNRAGG = 0,  &
      iNPRACG = 0, &
      iNSUBR = 0,  &
      iNSMLTR = 0, &
      iNGMLTR = 0, &
      iNPRACS = 0, &
      iNNUCCR = 0, &
      iNIACR = 0,  &
      iNIACRS = 0, &
      iNGRACS = 0, &
      iNSMLTS = 0, &
      iNSAGG = 0,  &
      iNPRCI = 0, &
      iNSCNG = 0, &
      iNSUBS = 0, &
      iPRC = 0, &
      iPRA = 0, &
      iPRE = 0


    ! More indices for Morrison budgets!!
    integer :: &
      iPCC = 0, & 
      iNNUCCC = 0, & 
      iNPSACWS = 0, &
      iNPRA = 0, & 
      iNPRC = 0, & 
      iNPSACWI = 0, &
      iNPSACWG = 0, &
      iNPRAI = 0, &
      iNMULTS = 0, & 
      iNMULTG = 0, & 
      iNMULTR = 0, & 
      iNMULTRG = 0, & 
      iNNUCCD = 0, & 
      iNSUBI = 0, & 
      iNGMLTG = 0, &
      iNSUBG = 0, &
      iNACT = 0

    integer :: &
      iSIZEFIX_NR = 0, &
      iSIZEFIX_NC = 0, &
      iSIZEFIX_NI = 0, &
      iSIZEFIX_NS = 0, &
      iSIZEFIX_NG = 0, &
      iNEGFIX_NR = 0, &
      iNEGFIX_NC = 0, &
      iNEGFIX_NI = 0, &
      iNEGFIX_NS = 0, &
      iNEGFIX_NG = 0, &
      iNIM_MORR_CL = 0, &
      iQC_INST = 0, &
      iQR_INST = 0, &
      iQI_INST = 0, &
      iQS_INST = 0, &
      iQG_INST = 0, &
      iNC_INST = 0, &
      iNR_INST = 0, &
      iNI_INST = 0, &
      iNS_INST = 0, &
      iNG_INST = 0, &
      iT_in_K_mc = 0, &
      ihl_on_Cp_residual = 0, &
      iqto_residual = 0


    !====================================================================
    !                           ZM Indices
    !====================================================================
    integer :: & 
       iwp2 = 0, & 
       irtp2 = 0, & 
       ithlp2 = 0, & 
       irtpthlp = 0, & 
       iwprtp = 0, & 
       iwpthlp = 0, &
       iwpup2 = 0, &
       iwpvp2 = 0, &
       iwp2up2 = 0, &
       iwp2vp2 = 0, & 
       iwp4 = 0, & 
       iwpthvp = 0, & 
       irtpthvp = 0, & 
       ithlpthvp = 0, & 
       itau_zm = 0, & 
       iKh_zm = 0, & 
       iwprcp = 0, & 
       irc_coef_zm = 0, &
       ithlprcp = 0, & 
       irtprcp = 0, & 
       ircp2 = 0, & 
       iupwp = 0, & 
       ivpwp = 0, &
       iupthlp = 0, &
       iuprtp = 0, &
       ivpthlp = 0, &
       ivprtp = 0, &
       iupthvp = 0, &
       iuprcp = 0, &
       ivpthvp = 0, &
       ivprcp = 0, &
       iSkw_zm = 0, &
       iSkthl_zm = 0, &
       iSkrt_zm = 0

    integer :: &
       irho_zm = 0, & 
       isigma_sqd_w = 0, &
       irho_ds_zm = 0, &
       ithv_ds_zm = 0, &
       iem = 0, & 
       ishear = 0, & ! Brian
       imean_w_up = 0, &
       imean_w_down = 0, &
       iFrad = 0, & 
       iFrad_LW = 0,   & ! Brian
       iFrad_SW = 0,   & ! Brian
       iFrad_LW_up = 0,   & 
       iFrad_SW_up = 0,   & 
       iFrad_LW_down = 0,   & 
       iFrad_SW_down = 0,   & 
       iFprec = 0,          & ! Brian
       iFcsed = 0             ! Brian

     ! Stability correction applied to Kh_N2_zm (diffusion on rtm and thlm)
     integer :: &
       istability_correction = 0 ! schemena


    integer, dimension(:), allocatable :: &
      iK_hm

    ! Skewness Functions on stats_zm grid
    integer :: &
      igamma_Skw_fnc = 0,  &
      iC6rt_Skw_fnc = 0,   &
      iC6thl_Skw_fnc = 0,  &
      iC7_Skw_fnc = 0,     &
      iC1_Skw_fnc = 0,     &
      ibrunt_vaisala_freq_sqd = 0, &
      ibrunt_vaisala_freq_sqd_splat = 0, &
      ibrunt_vaisala_freq_sqd_mixed = 0, &
      ibrunt_vaisala_freq_sqd_moist = 0, &
      ibrunt_vaisala_freq_sqd_dry = 0, &
      ibrunt_vaisala_freq_sqd_smth = 0, &
      ibrunt_freq_pos = 0, &
      ibrunt_freq_out_cloud = 0, &
      iRi_zm = 0,          &
      ishear_sqd = 0,      &
      iC6_term = 0


    integer :: &
      icoef_wp4_implicit = 0


    ! Covariance of w and cloud droplet concentration, < w'N_c' >
    integer :: &
      iwpNcp = 0


    ! Sedimentation velocities
    integer :: & 
      iVNr = 0,    &
      iVrr = 0,    &
      iVNc = 0,    &
      iVrc = 0,    &
      iVNs = 0, &
      iVrs = 0, &
      iVNi = 0,  &
      iVri = 0,  &
      iVrg = 0


    ! Covariance of sedimentation velocity and hydrometeor, <V_xx'x_x'>.
    integer :: &
      iVrrprrp = 0,         &
      iVNrpNrp = 0,         &
      iVrrprrp_expcalc = 0, &
      iVNrpNrp_expcalc = 0


    integer :: & 
       iwp2_bt = 0, & 
       iwp2_ma = 0, & 
       iwp2_ta = 0, & 
       iwp2_ac = 0, & 
       iwp2_bp = 0, & 
       iwp2_pr1 = 0, & 
       iwp2_pr2 = 0, & 
       iwp2_pr3 = 0, &
       iwp2_pr_dfsn = 0, &  
       iwp2_dp1 = 0, & 
       iwp2_dp2 = 0, &
       iwp2_sdmp = 0, &
       iwp2_pd = 0, & 
       iwp2_cl = 0, &
       iwp2_sf = 0, &
       iwp2_splat = 0


    integer :: & 
       iwprtp_bt = 0,      & 
       iwprtp_ma = 0,      & 
       iwprtp_ta = 0,      & 
       iwprtp_tp = 0,      & 
       iwprtp_ac = 0,      & 
       iwprtp_bp = 0,      & 
       iwprtp_pr1 = 0,     & 
       iwprtp_pr2 = 0,     & 
       iwprtp_pr3 = 0,     & 
       iwprtp_dp1 = 0,     &
       iwprtp_mfl = 0,     &
       iwprtp_cl = 0,      & 
       iwprtp_sicl = 0,    & 
       iwprtp_pd = 0,      &
       iwprtp_forcing = 0, &
       iwprtp_mc = 0


    integer :: & 
       iwpthlp_bt = 0,      & 
       iwpthlp_ma = 0,      & 
       iwpthlp_ta = 0,      & 
       iwpthlp_tp = 0,      & 
       iwpthlp_ac = 0,      & 
       iwpthlp_bp = 0,      & 
       iwpthlp_pr1 = 0,     & 
       iwpthlp_pr2 = 0,     & 
       iwpthlp_pr3 = 0,     & 
       iwpthlp_dp1 = 0,     &
       iwpthlp_mfl = 0,     &
       iwpthlp_cl = 0,      & 
       iwpthlp_sicl = 0,    &
       iwpthlp_forcing = 0, &
       iwpthlp_mc = 0


    integer :: & 
       iupwp_bt = 0,  &
       iupwp_ma = 0,  &
       iupwp_ta = 0,  &
       iupwp_tp = 0,  &
       iupwp_ac = 0,  &
       iupwp_bp = 0,  &
       iupwp_pr1 = 0, &
       iupwp_pr2 = 0, &
       iupwp_pr3 = 0, &
       iupwp_pr4 = 0, &
       iupwp_dp1 = 0, &
       iupwp_mfl = 0, &
       iupwp_cl = 0


    integer :: & 
       ivpwp_bt = 0,  &
       ivpwp_ma = 0,  &
       ivpwp_ta = 0,  &
       ivpwp_tp = 0,  &
       ivpwp_ac = 0,  &
       ivpwp_bp = 0,  &
       ivpwp_pr1 = 0, &
       ivpwp_pr2 = 0, &
       ivpwp_pr3 = 0, &
       ivpwp_pr4 = 0, &
       ivpwp_dp1 = 0, &
       ivpwp_mfl = 0, &
       ivpwp_cl = 0


  !    Dr. Golaz's new variance budget terms
  !    qt was changed to rt to avoid confusion

    integer :: & 
       irtp2_bt = 0,      & 
       irtp2_ma = 0,      & 
       irtp2_ta = 0,      & 
       irtp2_tp = 0,      & 
       irtp2_dp1 = 0,     & 
       irtp2_dp2 = 0,     & 
       irtp2_pd = 0,      & 
       irtp2_cl = 0,      &
       irtp2_sf = 0,      &
       irtp2_forcing = 0, &
       irtp2_mc = 0
       

    integer :: & 
       ithlp2_bt = 0,      & 
       ithlp2_ma = 0,      & 
       ithlp2_ta = 0,      & 
       ithlp2_tp = 0,      & 
       ithlp2_dp1 = 0,     & 
       ithlp2_dp2 = 0,     & 
       ithlp2_pd = 0,      & 
       ithlp2_cl = 0,      &
       ithlp2_sf = 0,      &
       ithlp2_forcing = 0, &
       ithlp2_mc = 0


    integer :: & 
      irtpthlp_bt = 0,      & 
      irtpthlp_ma = 0,      & 
      irtpthlp_ta = 0,      & 
      irtpthlp_tp1 = 0,     & 
      irtpthlp_tp2 = 0,     & 
      irtpthlp_dp1 = 0,     & 
      irtpthlp_dp2 = 0,     & 
      irtpthlp_cl = 0,      &
      irtpthlp_sf = 0,      &
      irtpthlp_forcing = 0, &
      irtpthlp_mc = 0


    integer :: & 
      iup2 = 0, & 
      ivp2 = 0


    integer :: & 
      iup2_bt = 0, & 
      iup2_ta = 0, & 
      iup2_tp = 0, & 
      iup2_ma = 0, & 
      iup2_dp1 = 0, & 
      iup2_dp2 = 0, & 
      iup2_pr1 = 0, & 
      iup2_pr2 = 0, & 
      iup2_sdmp = 0, & 
      iup2_pd = 0, & 
      iup2_cl = 0, &
      iup2_sf = 0, &
      iup2_splat = 0, &
      ivp2_bt = 0, & 
      ivp2_ta = 0, & 
      ivp2_tp = 0, & 
      ivp2_ma = 0, & 
      ivp2_dp1 = 0, & 
      ivp2_dp2 = 0, & 
      ivp2_pr1 = 0, & 
      ivp2_pr2 = 0, & 
      ivp2_sdmp = 0, & 
      ivp2_pd = 0, & 
      ivp2_cl = 0, &
      ivp2_sf = 0, &
      ivp2_splat = 0


  !       Passive scalars.  Note that floating point roundoff may make
  !       mathematically equivalent variables different values.
    integer,allocatable, dimension(:) :: & 
      isclrprtp,           & ! sclr'(1)rt'     / rt'^2
      isclrp2,             & ! sclr'(1)^2      / rt'^2
      isclrpthvp,          & ! sclr'(1)th_v'   / rt'th_v' 
      isclrpthlp,          & ! sclr'(1)th_l'   / rt'th_l' 
      isclrprcp,           & ! sclr'(1)rc'     / rt'rc'
      iwpsclrp,            & ! w'slcr'(1)      / w'rt'
      iwp2sclrp,           & ! w'^2 sclr'(1)   / w'^2 rt'
      iwpsclrp2,           & ! w'sclr'(1)^2    / w'rt'^2
      iwpsclrprtp,         & ! w'sclr'(1)rt'   / w'rt'^2
      iwpsclrpthlp           ! w'sclr'(1)th_l' / w'rt'th_l'


    integer, allocatable, dimension(:) :: & 
       iwpedsclrp ! eddy sclr'(1)w'


    !====================================================================
    !                         RAD ZT Indices
    !====================================================================
    integer :: &
      iT_in_K_rad = 0, &
      ircil_rad = 0, &
      io3l_rad = 0, &
      irsm_rad = 0, &
      ircm_in_cloud_rad = 0, &
      icloud_frac_rad = 0, & 
      iice_supersat_frac_rad = 0, &
      iradht_rad = 0, &
      iradht_LW_rad = 0, &
      iradht_SW_rad = 0, &
      ip_in_mb_rad = 0, &
      isp_humidity_rad = 0


    !====================================================================
    !                       RAD ZM Indices
    !====================================================================
    integer :: &
      iFrad_LW_rad = 0, &
      iFrad_SW_rad = 0, &
      iFrad_SW_up_rad = 0, &
      iFrad_LW_up_rad = 0, &
      iFrad_SW_down_rad = 0, &
      iFrad_LW_down_rad = 0


    !====================================================================
    !                           SFC Indices
    !====================================================================
    integer :: & 
      iustar = 0, &
      isoil_heat_flux = 0,&
      iveg_T_in_K = 0,&
      isfc_soil_T_in_K = 0, &
      ideep_soil_T_in_K = 0,& 
      ilh = 0, & 
      ish = 0, & 
      icc = 0, & 
      ilwp = 0, &
      ivwp = 0, &        ! nielsenb
      iiwp = 0, &        ! nielsenb
      iswp = 0, &        ! nielsenb
      irwp = 0, &
      iz_cloud_base = 0, & 
      iz_inversion = 0, & 
      iprecip_rate_sfc = 0,    &    ! Brian
      irain_flux_sfc = 0,   &    ! Brian
      irrm_sfc = 0, & ! Brian
      iwpthlp_sfc = 0, &
      iprecip_frac_tol = 0

    integer :: &
      iwprtp_sfc = 0, &
      iupwp_sfc = 0, &
      ivpwp_sfc = 0, &
      ithlm_vert_avg = 0, &
      irtm_vert_avg = 0, &
      ium_vert_avg = 0, &
      ivm_vert_avg = 0, &
      iwp2_vert_avg = 0, & ! nielsenb
      iup2_vert_avg = 0, &
      ivp2_vert_avg = 0, &
      irtp2_vert_avg = 0, &
      ithlp2_vert_avg = 0, &
      iT_sfc = 0         ! kcwhite

    integer :: &
      itot_vartn_normlzd_rtm = 0, &
      itot_vartn_normlzd_thlm = 0, &
      itot_vartn_normlzd_wprtp = 0

    integer :: & 
      iwp23_matrix_condt_num = 0, & 
      irtm_matrix_condt_num = 0, & 
      ithlm_matrix_condt_num = 0, & 
      irtp2_matrix_condt_num = 0, & 
      ithlp2_matrix_condt_num = 0, & 
      irtpthlp_matrix_condt_num = 0, & 
      iup2_vp2_matrix_condt_num = 0, & 
      iwindm_matrix_condt_num = 0

    integer :: & 
      imorr_snow_rate = 0


    integer :: &
      irtm_spur_src = 0,    &
      ithlm_spur_src = 0


    integer :: &
      iSkw_velocity = 0, & ! Skewness velocity
      iwp3_zm = 0, &
      ithlp3_zm = 0, &
      irtp3_zm = 0, &
      ia3_coef = 0, &
      ia3_coef_zt = 0

    integer :: &
      iwp3_on_wp2 = 0,         & ! w'^3 / w'^2                   [m/s]
      iwp3_on_wp2_zt = 0,      & ! w'^3 / w'^2                   [m/s]
      iwp3_on_wp2_cfl_num = 0    ! ( w'^3 / w'^2 ) * dt / dz     [-]

    integer :: & 
      ilh_morr_snow_rate = 0

    integer :: & 
      ilh_vwp = 0, &
      ilh_lwp = 0, &
      ilh_sample_weights_sum = 0, &
      ilh_sample_weights_avg = 0


    integer :: &
      icloud_frac_refined = 0, &
      ircm_refined = 0

    integer :: &
      irtp2_from_chi = 0

  end type stats_metadata_type
!$omp declare mapper (stats_metadata_type::x) map ( &
!$omp  x%stats_tout &
!$omp , x%l_allow_small_stats_tout &
!$omp , x%l_stats_last &
!$omp , x%fname_sfc &
!$omp , x%ithlm &
!$omp , x%ithvm &
!$omp , x%irtm &
!$omp , x%ircm &
!$omp , x%irvm &
!$omp , x%ium &
!$omp , x%ivm &
!$omp , x%iwm_zt &
!$omp , x%iwm_zm &
!$omp , x%ium_ref &
!$omp , x%ivm_ref &
!$omp , x%iug &
!$omp , x%ivg &
!$omp , x%icloud_frac &
!$omp , x%iice_supersat_frac &
!$omp , x%ircm_in_layer &
!$omp , x%ircm_in_cloud &
!$omp , x%icloud_cover &
!$omp , x%ip_in_pa &
!$omp , x%iexner &
!$omp , x%irho_ds_zt &
!$omp , x%ithv_ds_zt &
!$omp , x%ilscale &
!$omp , x%iwp3 &
!$omp , x%ithlp3 &
!$omp , x%irtp3 &
!$omp , x%iwpthlp2 &
!$omp , x%iwp2thlp &
!$omp , x%iwprtp2 &
!$omp , x%iwp2rtp &
!$omp , x%iskw_zt &
!$omp , x%iskthl_zt &
!$omp , x%iskrt_zt &
!$omp , x%ircm_supersat_adj &
!$omp , x%ilscale_up &
!$omp , x%ilscale_down &
!$omp , x%ilscale_pert_1 &
!$omp , x%ilscale_pert_2 &
!$omp , x%itau_zt &
!$omp , x%iinvrs_tau_zt &
!$omp , x%ikh_zt &
!$omp , x%iwp2thvp &
!$omp , x%iwp2rcp &
!$omp , x%iw_up_in_cloud &
!$omp , x%iw_down_in_cloud &
!$omp , x%icld_updr_frac &
!$omp , x%icld_downdr_frac &
!$omp , x%iwprtpthlp &
!$omp , x%irc_coef &
!$omp , x%isigma_sqd_w_zt &
!$omp , x%irho &
!$omp , x%iinvrs_tau_zm &
!$omp , x%iinvrs_tau_xp2_zm &
!$omp , x%iinvrs_tau_wp2_zm &
!$omp , x%iinvrs_tau_wpxp_zm &
!$omp , x%iinvrs_tau_wp3_zm &
!$omp , x%iinvrs_tau_no_n2_zm &
!$omp , x%iinvrs_tau_bkgnd &
!$omp , x%iinvrs_tau_sfc &
!$omp , x%iinvrs_tau_shear &
!$omp , x%itau_no_n2_zm &
!$omp , x%itau_xp2_zm &
!$omp , x%itau_wp2_zm &
!$omp , x%itau_wp3_zm &
!$omp , x%itau_wpxp_zm &
!$omp , x%ihm_1 &
!$omp , x%ihm_2 &
!$omp , x%iprecip_frac &
!$omp , x%iprecip_frac_1 &
!$omp , x%iprecip_frac_2 &
!$omp , x%incnm &
!$omp , x%imu_hm_1 &
!$omp , x%imu_hm_2 &
!$omp , x%imu_hm_1_n &
!$omp , x%imu_hm_2_n &
!$omp , x%isigma_hm_1 &
!$omp , x%isigma_hm_2 &
!$omp , x%isigma_hm_1_n &
!$omp , x%isigma_hm_2_n &
!$omp , x%icorr_w_hm_1 &
!$omp , x%icorr_w_hm_2 &
!$omp , x%icorr_chi_hm_1 &
!$omp , x%icorr_chi_hm_2 &
!$omp , x%icorr_eta_hm_1 &
!$omp , x%icorr_eta_hm_2 &
!$omp , x%icorr_ncn_hm_1 &
!$omp , x%icorr_ncn_hm_2 &
!$omp , x%icorr_w_hm_1_n &
!$omp , x%icorr_w_hm_2_n &
!$omp , x%icorr_chi_hm_1_n &
!$omp , x%icorr_chi_hm_2_n &
!$omp , x%icorr_eta_hm_1_n &
!$omp , x%icorr_eta_hm_2_n &
!$omp , x%icorr_ncn_hm_1_n &
!$omp , x%icorr_ncn_hm_2_n &
!$omp , x%icorr_hmx_hmy_1 &
!$omp , x%icorr_hmx_hmy_2 &
!$omp , x%icorr_hmx_hmy_1_n &
!$omp , x%icorr_hmx_hmy_2_n &
!$omp , x%imu_ncn_1 &
!$omp , x%imu_ncn_2 &
!$omp , x%imu_ncn_1_n &
!$omp , x%imu_ncn_2_n &
!$omp , x%isigma_ncn_1 &
!$omp , x%isigma_ncn_2 &
!$omp , x%isigma_ncn_1_n &
!$omp , x%isigma_ncn_2_n &
!$omp , x%icorr_w_chi_1_ca &
!$omp , x%icorr_w_chi_2_ca &
!$omp , x%icorr_w_eta_1_ca &
!$omp , x%icorr_w_eta_2_ca &
!$omp , x%icorr_w_ncn_1 &
!$omp , x%icorr_w_ncn_2 &
!$omp , x%icorr_chi_eta_1_ca &
!$omp , x%icorr_chi_eta_2_ca &
!$omp , x%icorr_chi_ncn_1 &
!$omp , x%icorr_chi_ncn_2 &
!$omp , x%icorr_eta_ncn_1 &
!$omp , x%icorr_eta_ncn_2 &
!$omp , x%icorr_w_ncn_1_n &
!$omp , x%icorr_w_ncn_2_n &
!$omp , x%icorr_chi_ncn_1_n &
!$omp , x%icorr_chi_ncn_2_n &
!$omp , x%icorr_eta_ncn_1_n &
!$omp , x%icorr_eta_ncn_2_n &
!$omp , x%isilhs_variance_category &
!$omp , x%ilh_samp_frac_category &
!$omp , x%inccnm &
!$omp , x%inc_in_cloud &
!$omp , x%inc_activated &
!$omp , x%irsati &
!$omp , x%irel_humidity &
!$omp , x%iakstd &
!$omp , x%iakstd_cld &
!$omp , x%iakm_rcm &
!$omp , x%iakm_rcc &
!$omp , x%irfrzm &
!$omp , x%ic11_skw_fnc &
!$omp , x%icloud_frac_zm &
!$omp , x%iice_supersat_frac_zm &
!$omp , x%ircm_zm &
!$omp , x%irtm_zm &
!$omp , x%ithlm_zm &
!$omp , x%iw_1_zm &
!$omp , x%iw_2_zm &
!$omp , x%ivarnce_w_1_zm &
!$omp , x%ivarnce_w_2_zm &
!$omp , x%imixt_frac_zm &
!$omp , x%ilh_rcm_avg &
!$omp , x%ik_lh_start &
!$omp , x%ingm &
!$omp , x%it_in_k &
!$omp , x%ieff_rad_cloud &
!$omp , x%ieff_rad_ice &
!$omp , x%ieff_rad_snow &
!$omp , x%ieff_rad_rain &
!$omp , x%ieff_rad_graupel &
!$omp , x%irsm &
!$omp , x%irgm &
!$omp , x%irim &
!$omp , x%iu_t_cm &
!$omp , x%ithlm_cl &
!$omp , x%ithlm_mfl_min &
!$omp , x%ithlm_mfl_max &
!$omp , x%iwpthlp_enter_mfl &
!$omp , x%iwpthlp_exit_mfl &
!$omp , x%iwpthlp_mfl_min &
!$omp , x%iwpthlp_mfl_max &
!$omp , x%irtm_mfl_min &
!$omp , x%irtm_mfl_max &
!$omp , x%iwprtp_enter_mfl &
!$omp , x%iwprtp_exit_mfl &
!$omp , x%iwprtp_mfl_min &
!$omp , x%iwprtp_mfl_max &
!$omp , x%ithlm_enter_mfl &
!$omp , x%ithlm_exit_mfl &
!$omp , x%ithlm_old &
!$omp , x%ithlm_without_ta &
!$omp , x%irtm_enter_mfl &
!$omp , x%irtm_exit_mfl &
!$omp , x%irtm_old &
!$omp , x%irtm_without_ta &
!$omp , x%iwp3_bt &
!$omp , x%iwp3_ma &
!$omp , x%iwp3_ta &
!$omp , x%iwp3_tp &
!$omp , x%iwp3_ac &
!$omp , x%iwp3_bp1 &
!$omp , x%iwp3_pr_tp &
!$omp , x%iwp3_pr_turb &
!$omp , x%iwp3_pr_dfsn &
!$omp , x%iwp3_pr1 &
!$omp , x%iwp3_pr2 &
!$omp , x%iwp3_pr3 &
!$omp , x%iwp3_dp1 &
!$omp , x%iwp3_sdmp &
!$omp , x%iwp3_cl &
!$omp , x%iwp3_splat &
!$omp , x%irtp3_bt &
!$omp , x%irtp3_tp &
!$omp , x%irtp3_ac &
!$omp , x%irtp3_dp &
!$omp , x%ithlp3_bt &
!$omp , x%ithlp3_tp &
!$omp , x%ithlp3_ac &
!$omp , x%ithlp3_dp &
!$omp , x%irrm_bt &
!$omp , x%irrm_ma &
!$omp , x%irrm_ta &
!$omp , x%irrm_sd &
!$omp , x%irrm_ts &
!$omp , x%irrm_sd_morr &
!$omp , x%irrm_evap &
!$omp , x%irrm_auto &
!$omp , x%irrm_accr &
!$omp , x%irrm_evap_adj &
!$omp , x%irrm_src_adj &
!$omp , x%irrm_mc_nonadj &
!$omp , x%irrm_mc &
!$omp , x%irrm_hf &
!$omp , x%irrm_wvhf &
!$omp , x%irrm_cl &
!$omp , x%inrm_bt &
!$omp , x%inrm_ma &
!$omp , x%inrm_ta &
!$omp , x%inrm_sd &
!$omp , x%inrm_ts &
!$omp , x%inrm_evap &
!$omp , x%inrm_auto &
!$omp , x%inrm_evap_adj &
!$omp , x%inrm_src_adj &
!$omp , x%inrm_mc &
!$omp , x%inrm_cl &
!$omp , x%irsm_bt &
!$omp , x%irsm_ma &
!$omp , x%irsm_sd &
!$omp , x%irsm_sd_morr &
!$omp , x%irsm_ta &
!$omp , x%irsm_mc &
!$omp , x%irsm_hf &
!$omp , x%irsm_wvhf &
!$omp , x%irsm_cl &
!$omp , x%irsm_sd_morr_int &
!$omp , x%irgm_bt &
!$omp , x%irgm_ma &
!$omp , x%irgm_sd &
!$omp , x%irgm_sd_morr &
!$omp , x%irgm_ta &
!$omp , x%irgm_mc &
!$omp , x%irgm_hf &
!$omp , x%irgm_wvhf &
!$omp , x%irgm_cl &
!$omp , x%irim_bt &
!$omp , x%irim_ma &
!$omp , x%irim_sd &
!$omp , x%irim_sd_mg_morr &
!$omp , x%irim_ta &
!$omp , x%irim_mc &
!$omp , x%irim_hf &
!$omp , x%irim_wvhf &
!$omp , x%irim_cl &
!$omp , x%insm_bt &
!$omp , x%insm_ma &
!$omp , x%insm_sd &
!$omp , x%insm_ta &
!$omp , x%insm_mc &
!$omp , x%insm_cl &
!$omp , x%ingm_bt &
!$omp , x%ingm_ma &
!$omp , x%ingm_sd &
!$omp , x%ingm_ta &
!$omp , x%ingm_mc &
!$omp , x%ingm_cl &
!$omp , x%inim_bt &
!$omp , x%inim_ma &
!$omp , x%inim_sd &
!$omp , x%inim_ta &
!$omp , x%inim_mc &
!$omp , x%inim_cl &
!$omp , x%incm_bt &
!$omp , x%incm_ma &
!$omp , x%incm_ta &
!$omp , x%incm_mc &
!$omp , x%incm_cl &
!$omp , x%incm_act &
!$omp , x%ikk_mvr_variance_zt &
!$omp , x%ivm_bt &
!$omp , x%ivm_ma &
!$omp , x%ivm_ta &
!$omp , x%ivm_gf &
!$omp , x%ivm_cf &
!$omp , x%ivm_f &
!$omp , x%ivm_sdmp &
!$omp , x%ivm_ndg &
!$omp , x%ivm_mfl &
!$omp , x%ium_bt &
!$omp , x%ium_ma &
!$omp , x%ium_ta &
!$omp , x%ium_gf &
!$omp , x%ium_cf &
!$omp , x%ium_f &
!$omp , x%ium_sdmp &
!$omp , x%ium_ndg &
!$omp , x%ium_mfl &
!$omp , x%imixt_frac &
!$omp , x%iw_1 &
!$omp , x%iw_2 &
!$omp , x%ivarnce_w_1 &
!$omp , x%ivarnce_w_2 &
!$omp , x%ithl_1 &
!$omp , x%ithl_2 &
!$omp , x%ivarnce_thl_1 &
!$omp , x%ivarnce_thl_2 &
!$omp , x%irt_1 &
!$omp , x%irt_2 &
!$omp , x%ivarnce_rt_1 &
!$omp , x%ivarnce_rt_2 &
!$omp , x%irc_1 &
!$omp , x%irc_2 &
!$omp , x%irsatl_1 &
!$omp , x%irsatl_2 &
!$omp , x%icloud_frac_1 &
!$omp , x%icloud_frac_2 &
!$omp , x%ichi_1 &
!$omp , x%ichi_2 &
!$omp , x%istdev_chi_1 &
!$omp , x%istdev_chi_2 &
!$omp , x%ichip2 &
!$omp , x%istdev_eta_1 &
!$omp , x%istdev_eta_2 &
!$omp , x%icovar_chi_eta_1 &
!$omp , x%icovar_chi_eta_2 &
!$omp , x%icorr_w_chi_1 &
!$omp , x%icorr_w_chi_2 &
!$omp , x%icorr_w_eta_1 &
!$omp , x%icorr_w_eta_2 &
!$omp , x%icorr_chi_eta_1 &
!$omp , x%icorr_chi_eta_2 &
!$omp , x%icorr_w_rt_1 &
!$omp , x%icorr_w_rt_2 &
!$omp , x%icorr_w_thl_1 &
!$omp , x%icorr_w_thl_2 &
!$omp , x%icorr_rt_thl_1 &
!$omp , x%icorr_rt_thl_2 &
!$omp , x%icrt_1 &
!$omp , x%icrt_2 &
!$omp , x%icthl_1 &
!$omp , x%icthl_2 &
!$omp , x%if_w &
!$omp , x%if_rt &
!$omp , x%if_thl &
!$omp , x%imin_f_w &
!$omp , x%imax_f_w &
!$omp , x%imin_f_rt &
!$omp , x%imax_f_rt &
!$omp , x%imin_f_thl &
!$omp , x%imax_f_thl &
!$omp , x%icoef_wprtp2_implicit &
!$omp , x%iterm_wprtp2_explicit &
!$omp , x%icoef_wpthlp2_implicit &
!$omp , x%iterm_wpthlp2_explicit &
!$omp , x%icoef_wprtpthlp_implicit &
!$omp , x%iterm_wprtpthlp_explicit &
!$omp , x%icoef_wp2rtp_implicit &
!$omp , x%iterm_wp2rtp_explicit &
!$omp , x%icoef_wp2thlp_implicit &
!$omp , x%iterm_wp2thlp_explicit &
!$omp , x%iwp2_zt &
!$omp , x%ithlp2_zt &
!$omp , x%iwpthlp_zt &
!$omp , x%iwprtp_zt &
!$omp , x%irtp2_zt &
!$omp , x%irtpthlp_zt &
!$omp , x%iup2_zt &
!$omp , x%ivp2_zt &
!$omp , x%iupwp_zt &
!$omp , x%ivpwp_zt &
!$omp , x%iwp2hmp &
!$omp , x%ihydrometp2 &
!$omp , x%iwphydrometp &
!$omp , x%irtphmp &
!$omp , x%ithlphmp &
!$omp , x%ihmxphmyp &
!$omp , x%ihmp2_zt &
!$omp , x%ichi &
!$omp , x%isclrm_f &
!$omp , x%ifulwcl &
!$omp , x%ifdlwcl &
!$omp , x%ifdswcl &
!$omp , x%ifuswcl &
!$omp , x%iedsclrm_f &
!$omp , x%ilh_nim_mc &
!$omp , x%ilh_m_vol_rad_rain &
!$omp , x%ilh_rrm_mc_nonadj &
!$omp , x%ilh_nrm_evap_adj &
!$omp , x%ilh_vnr &
!$omp , x%ilh_rrm &
!$omp , x%ilh_nrm &
!$omp , x%ilh_rim &
!$omp , x%ilh_nim &
!$omp , x%ilh_rsm &
!$omp , x%ilh_nsm &
!$omp , x%ilh_rgm &
!$omp , x%ilh_ngm &
!$omp , x%ilh_thlm &
!$omp , x%ilh_rcm &
!$omp , x%ilh_ncm &
!$omp , x%ilh_ncnm &
!$omp , x%ilh_rvm &
!$omp , x%ilh_wm &
!$omp , x%ilh_cloud_frac &
!$omp , x%ilh_chi &
!$omp , x%ilh_eta &
!$omp , x%ilh_precip_frac &
!$omp , x%ilh_mixt_frac &
!$omp , x%ilh_cloud_frac_unweighted &
!$omp , x%ilh_precip_frac_unweighted &
!$omp , x%ilh_mixt_frac_unweighted &
!$omp , x%ilh_wp2_zt &
!$omp , x%ilh_nrp2_zt &
!$omp , x%ilh_ncnp2_zt &
!$omp , x%ilh_ncp2_zt &
!$omp , x%ilh_rcp2_zt &
!$omp , x%ilh_rtp2_zt &
!$omp , x%ilh_thlp2_zt &
!$omp , x%ilh_rrp2_zt &
!$omp , x%ilh_chip2 &
!$omp , x%ilh_rtp2_mc &
!$omp , x%ilh_thlp2_mc &
!$omp , x%ilh_wprtp_mc &
!$omp , x%ilh_wpthlp_mc &
!$omp , x%ilh_rtpthlp_mc &
!$omp , x%ipsmlt &
!$omp , x%ievpms &
!$omp , x%ipracs &
!$omp , x%ievpmg &
!$omp , x%ipracg &
!$omp , x%ipgmlt &
!$omp , x%imnuccc &
!$omp , x%ipsacws &
!$omp , x%ipsacwi &
!$omp , x%iqmults &
!$omp , x%iqmultg &
!$omp , x%ipsacwg &
!$omp , x%ipgsacw &
!$omp , x%iprd &
!$omp , x%iprci &
!$omp , x%iprai &
!$omp , x%iqmultr &
!$omp , x%iqmultrg &
!$omp , x%imnuccd &
!$omp , x%ipraci &
!$omp , x%ipracis &
!$omp , x%ieprd &
!$omp , x%imnuccr &
!$omp , x%ipiacr &
!$omp , x%ipiacrs &
!$omp , x%ipgracs &
!$omp , x%iprds &
!$omp , x%ieprds &
!$omp , x%ipsacr &
!$omp , x%iprdg &
!$omp , x%ieprdg &
!$omp , x%ingsten &
!$omp , x%inrsten &
!$omp , x%inisten &
!$omp , x%inssten &
!$omp , x%incsten &
!$omp , x%inprc1 &
!$omp , x%inragg &
!$omp , x%inpracg &
!$omp , x%insubr &
!$omp , x%insmltr &
!$omp , x%ingmltr &
!$omp , x%inpracs &
!$omp , x%innuccr &
!$omp , x%iniacr &
!$omp , x%iniacrs &
!$omp , x%ingracs &
!$omp , x%insmlts &
!$omp , x%insagg &
!$omp , x%inprci &
!$omp , x%inscng &
!$omp , x%insubs &
!$omp , x%iprc &
!$omp , x%ipra &
!$omp , x%ipre &
!$omp , x%ipcc &
!$omp , x%innuccc &
!$omp , x%inpsacws &
!$omp , x%inpra &
!$omp , x%inprc &
!$omp , x%inpsacwi &
!$omp , x%inpsacwg &
!$omp , x%inprai &
!$omp , x%inmults &
!$omp , x%inmultg &
!$omp , x%inmultr &
!$omp , x%inmultrg &
!$omp , x%innuccd &
!$omp , x%insubi &
!$omp , x%ingmltg &
!$omp , x%insubg &
!$omp , x%inact &
!$omp , x%isizefix_nr &
!$omp , x%isizefix_nc &
!$omp , x%isizefix_ni &
!$omp , x%isizefix_ns &
!$omp , x%isizefix_ng &
!$omp , x%inegfix_nr &
!$omp , x%inegfix_nc &
!$omp , x%inegfix_ni &
!$omp , x%inegfix_ns &
!$omp , x%inegfix_ng &
!$omp , x%inim_morr_cl &
!$omp , x%iqc_inst &
!$omp , x%iqr_inst &
!$omp , x%iqi_inst &
!$omp , x%iqs_inst &
!$omp , x%iqg_inst &
!$omp , x%inc_inst &
!$omp , x%inr_inst &
!$omp , x%ini_inst &
!$omp , x%ins_inst &
!$omp , x%ing_inst &
!$omp , x%it_in_k_mc &
!$omp , x%ihl_on_cp_residual &
!$omp , x%iqto_residual &
!$omp , x%iwp2 &
!$omp , x%irtp2 &
!$omp , x%ithlp2 &
!$omp , x%irtpthlp &
!$omp , x%iwprtp &
!$omp , x%iwpthlp &
!$omp , x%iwpup2 &
!$omp , x%iwpvp2 &
!$omp , x%iwp2up2 &
!$omp , x%iwp2vp2 &
!$omp , x%iwp4 &
!$omp , x%iwpthvp &
!$omp , x%irtpthvp &
!$omp , x%ithlpthvp &
!$omp , x%itau_zm &
!$omp , x%ikh_zm &
!$omp , x%iwprcp &
!$omp , x%irc_coef_zm &
!$omp , x%ithlprcp &
!$omp , x%irtprcp &
!$omp , x%ircp2 &
!$omp , x%iupwp &
!$omp , x%ivpwp &
!$omp , x%iupthlp &
!$omp , x%iuprtp &
!$omp , x%ivpthlp &
!$omp , x%ivprtp &
!$omp , x%iupthvp &
!$omp , x%iuprcp &
!$omp , x%ivpthvp &
!$omp , x%ivprcp &
!$omp , x%iskw_zm &
!$omp , x%iskthl_zm &
!$omp , x%iskrt_zm &
!$omp , x%irho_zm &
!$omp , x%isigma_sqd_w &
!$omp , x%irho_ds_zm &
!$omp , x%ithv_ds_zm &
!$omp , x%iem &
!$omp , x%imean_w_up &
!$omp , x%imean_w_down &
!$omp , x%ifrad &
!$omp , x%ifrad_lw_up &
!$omp , x%ifrad_sw_up &
!$omp , x%ifrad_lw_down &
!$omp , x%ifrad_sw_down &
!$omp , x%ifcsed &
!$omp , x%istability_correction &
!$omp , x%ik_hm &
!$omp , x%igamma_skw_fnc &
!$omp , x%ic6rt_skw_fnc &
!$omp , x%ic6thl_skw_fnc &
!$omp , x%ic7_skw_fnc &
!$omp , x%ic1_skw_fnc &
!$omp , x%ibrunt_vaisala_freq_sqd &
!$omp , x%ibrunt_vaisala_freq_sqd_splat &
!$omp , x%ibrunt_vaisala_freq_sqd_mixed &
!$omp , x%ibrunt_vaisala_freq_sqd_moist &
!$omp , x%ibrunt_vaisala_freq_sqd_dry &
!$omp , x%ibrunt_vaisala_freq_sqd_smth &
!$omp , x%ibrunt_freq_pos &
!$omp , x%ibrunt_freq_out_cloud &
!$omp , x%iri_zm &
!$omp , x%ishear_sqd &
!$omp , x%ic6_term &
!$omp , x%icoef_wp4_implicit &
!$omp , x%iwpncp &
!$omp , x%ivnr &
!$omp , x%ivrr &
!$omp , x%ivnc &
!$omp , x%ivrc &
!$omp , x%ivns &
!$omp , x%ivrs &
!$omp , x%ivni &
!$omp , x%ivri &
!$omp , x%ivrg &
!$omp , x%ivrrprrp &
!$omp , x%ivnrpnrp &
!$omp , x%ivrrprrp_expcalc &
!$omp , x%ivnrpnrp_expcalc &
!$omp , x%iwp2_bt &
!$omp , x%iwp2_ma &
!$omp , x%iwp2_ta &
!$omp , x%iwp2_ac &
!$omp , x%iwp2_bp &
!$omp , x%iwp2_pr1 &
!$omp , x%iwp2_pr2 &
!$omp , x%iwp2_pr3 &
!$omp , x%iwp2_pr_dfsn &
!$omp , x%iwp2_dp1 &
!$omp , x%iwp2_dp2 &
!$omp , x%iwp2_sdmp &
!$omp , x%iwp2_pd &
!$omp , x%iwp2_cl &
!$omp , x%iwp2_sf &
!$omp , x%iwp2_splat &
!$omp , x%iwprtp_bt &
!$omp , x%iwprtp_ma &
!$omp , x%iwprtp_ta &
!$omp , x%iwprtp_tp &
!$omp , x%iwprtp_ac &
!$omp , x%iwprtp_bp &
!$omp , x%iwprtp_pr1 &
!$omp , x%iwprtp_pr2 &
!$omp , x%iwprtp_pr3 &
!$omp , x%iwprtp_dp1 &
!$omp , x%iwprtp_mfl &
!$omp , x%iwprtp_cl &
!$omp , x%iwprtp_sicl &
!$omp , x%iwprtp_pd &
!$omp , x%iwprtp_forcing &
!$omp , x%iwprtp_mc &
!$omp , x%iwpthlp_bt &
!$omp , x%iwpthlp_ma &
!$omp , x%iwpthlp_ta &
!$omp , x%iwpthlp_tp &
!$omp , x%iwpthlp_ac &
!$omp , x%iwpthlp_bp &
!$omp , x%iwpthlp_pr1 &
!$omp , x%iwpthlp_pr2 &
!$omp , x%iwpthlp_pr3 &
!$omp , x%iwpthlp_dp1 &
!$omp , x%iwpthlp_mfl &
!$omp , x%iwpthlp_cl &
!$omp , x%iwpthlp_sicl &
!$omp , x%iwpthlp_forcing &
!$omp , x%iwpthlp_mc &
!$omp , x%iupwp_bt &
!$omp , x%iupwp_ma &
!$omp , x%iupwp_ta &
!$omp , x%iupwp_tp &
!$omp , x%iupwp_ac &
!$omp , x%iupwp_bp &
!$omp , x%iupwp_pr1 &
!$omp , x%iupwp_pr2 &
!$omp , x%iupwp_pr3 &
!$omp , x%iupwp_pr4 &
!$omp , x%iupwp_dp1 &
!$omp , x%iupwp_mfl &
!$omp , x%iupwp_cl &
!$omp , x%ivpwp_bt &
!$omp , x%ivpwp_ma &
!$omp , x%ivpwp_ta &
!$omp , x%ivpwp_tp &
!$omp , x%ivpwp_ac &
!$omp , x%ivpwp_bp &
!$omp , x%ivpwp_pr1 &
!$omp , x%ivpwp_pr2 &
!$omp , x%ivpwp_pr3 &
!$omp , x%ivpwp_pr4 &
!$omp , x%ivpwp_dp1 &
!$omp , x%ivpwp_mfl &
!$omp , x%ivpwp_cl &
!$omp , x%irtp2_bt &
!$omp , x%irtp2_ma &
!$omp , x%irtp2_ta &
!$omp , x%irtp2_tp &
!$omp , x%irtp2_dp1 &
!$omp , x%irtp2_dp2 &
!$omp , x%irtp2_pd &
!$omp , x%irtp2_cl &
!$omp , x%irtp2_sf &
!$omp , x%irtp2_forcing &
!$omp , x%irtp2_mc &
!$omp , x%ithlp2_bt &
!$omp , x%ithlp2_ma &
!$omp , x%ithlp2_ta &
!$omp , x%ithlp2_tp &
!$omp , x%ithlp2_dp1 &
!$omp , x%ithlp2_dp2 &
!$omp , x%ithlp2_pd &
!$omp , x%ithlp2_cl &
!$omp , x%ithlp2_sf &
!$omp , x%ithlp2_forcing &
!$omp , x%ithlp2_mc &
!$omp , x%irtpthlp_bt &
!$omp , x%irtpthlp_ma &
!$omp , x%irtpthlp_ta &
!$omp , x%irtpthlp_tp1 &
!$omp , x%irtpthlp_tp2 &
!$omp , x%irtpthlp_dp1 &
!$omp , x%irtpthlp_dp2 &
!$omp , x%irtpthlp_cl &
!$omp , x%irtpthlp_sf &
!$omp , x%irtpthlp_forcing &
!$omp , x%irtpthlp_mc &
!$omp , x%iup2 &
!$omp , x%ivp2 &
!$omp , x%iup2_bt &
!$omp , x%iup2_ta &
!$omp , x%iup2_tp &
!$omp , x%iup2_ma &
!$omp , x%iup2_dp1 &
!$omp , x%iup2_dp2 &
!$omp , x%iup2_pr1 &
!$omp , x%iup2_pr2 &
!$omp , x%iup2_sdmp &
!$omp , x%iup2_pd &
!$omp , x%iup2_cl &
!$omp , x%iup2_sf &
!$omp , x%iup2_splat &
!$omp , x%ivp2_bt &
!$omp , x%ivp2_ta &
!$omp , x%ivp2_tp &
!$omp , x%ivp2_ma &
!$omp , x%ivp2_dp1 &
!$omp , x%ivp2_dp2 &
!$omp , x%ivp2_pr1 &
!$omp , x%ivp2_pr2 &
!$omp , x%ivp2_sdmp &
!$omp , x%ivp2_pd &
!$omp , x%ivp2_cl &
!$omp , x%ivp2_sf &
!$omp , x%ivp2_splat &
!$omp , x%iwpsclrpthlp &
!$omp , x%iwpedsclrp &
!$omp , x%it_in_k_rad &
!$omp , x%ircil_rad &
!$omp , x%io3l_rad &
!$omp , x%irsm_rad &
!$omp , x%ircm_in_cloud_rad &
!$omp , x%icloud_frac_rad &
!$omp , x%iice_supersat_frac_rad &
!$omp , x%iradht_rad &
!$omp , x%iradht_lw_rad &
!$omp , x%iradht_sw_rad &
!$omp , x%ip_in_mb_rad &
!$omp , x%isp_humidity_rad &
!$omp , x%ifrad_lw_rad &
!$omp , x%ifrad_sw_rad &
!$omp , x%ifrad_sw_up_rad &
!$omp , x%ifrad_lw_up_rad &
!$omp , x%ifrad_sw_down_rad &
!$omp , x%ifrad_lw_down_rad &
!$omp , x%iustar &
!$omp , x%isoil_heat_flux &
!$omp , x%iveg_t_in_k &
!$omp , x%isfc_soil_t_in_k &
!$omp , x%ideep_soil_t_in_k &
!$omp , x%ilh &
!$omp , x%ish &
!$omp , x%icc &
!$omp , x%ilwp &
!$omp , x%irwp &
!$omp , x%iz_cloud_base &
!$omp , x%iz_inversion &
!$omp , x%iwpthlp_sfc &
!$omp , x%iprecip_frac_tol &
!$omp , x%iwprtp_sfc &
!$omp , x%iupwp_sfc &
!$omp , x%ivpwp_sfc &
!$omp , x%ithlm_vert_avg &
!$omp , x%irtm_vert_avg &
!$omp , x%ium_vert_avg &
!$omp , x%ivm_vert_avg &
!$omp , x%iup2_vert_avg &
!$omp , x%ivp2_vert_avg &
!$omp , x%irtp2_vert_avg &
!$omp , x%ithlp2_vert_avg &
!$omp , x%it_sfc &
!$omp , x%itot_vartn_normlzd_rtm &
!$omp , x%itot_vartn_normlzd_thlm &
!$omp , x%itot_vartn_normlzd_wprtp &
!$omp , x%iwp23_matrix_condt_num &
!$omp , x%irtm_matrix_condt_num &
!$omp , x%ithlm_matrix_condt_num &
!$omp , x%irtp2_matrix_condt_num &
!$omp , x%ithlp2_matrix_condt_num &
!$omp , x%irtpthlp_matrix_condt_num &
!$omp , x%iup2_vp2_matrix_condt_num &
!$omp , x%iwindm_matrix_condt_num &
!$omp , x%imorr_snow_rate &
!$omp , x%irtm_spur_src &
!$omp , x%ithlm_spur_src &
!$omp , x%iwp3_zm &
!$omp , x%ithlp3_zm &
!$omp , x%irtp3_zm &
!$omp , x%ia3_coef &
!$omp , x%ia3_coef_zt &
!$omp , x%iwp3_on_wp2_cfl_num &
!$omp , x%ilh_morr_snow_rate &
!$omp , x%ilh_vwp &
!$omp , x%ilh_lwp &
!$omp , x%ilh_sample_weights_sum &
!$omp , x%ilh_sample_weights_avg &
!$omp , x%icloud_frac_refined &
!$omp , x%ircm_refined &
!$omp , x%irtp2_from_chi &
!$omp )

end module stats_variables


