!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module stats_zt

  implicit none

  private ! Default Scope

  public :: stats_init_zt

! Constant parameters
  integer, parameter, public :: nvarmax_zt = 750 ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_zt( vars_zt, l_error )

! Description:
!   Initializes array indices for zt

! Note:
!   All code that is within subroutine stats_init_zt, including variable
!   allocation code, is not called if l_stats is false.  This subroutine is
!   called only when l_stats is true.

!-----------------------------------------------------------------------

    use constants_clubb, only:  &
        fstderr ! Constant(s)

    use stats_variables, only: & 
        ithlm,  & ! Variable(s)
        iT_in_K, & 
        ithvm, & 
        irtm, & 
        ircm, &
        irfrzm, &
        irvm, & 
        ium, & 
        ivm, & 
        iwm_zt, & 
        ium_ref, &
        ivm_ref, &
        iug, & 
        ivg, & 
        icloud_frac, &
        iice_supersat_frac, &
        ircm_in_layer, &
        ircm_in_cloud, &
        icloud_cover, &
        ip_in_Pa, & 
        iexner, & 
        irho_ds_zt, &
        ithv_ds_zt, &
        iLscale
 
    use stats_variables, only: & 
        iwp3, & ! Variable(s)
        iwpthlp2, & 
        iwp2thlp, & 
        iwprtp2, & 
        iwp2rtp, & 
        iLscale_up, & 
        iLscale_down, & 
        itau_zt, & 
        iKh_zt, & 
        iwp2thvp, & 
        iwp2rcp, & 
        iwprtpthlp, & 
        isigma_sqd_w_zt

    use stats_variables, only: & 
        ihm1, & ! Variable(s)
        ihm2, &
        iLWP1, &
        iLWP2, &
        iprecip_frac, &
        iprecip_frac_1, &
        iprecip_frac_2, &
        iNcnm

    use stats_variables, only: &
        imu_hm_1,       & ! Variable(s)
        imu_hm_2,       &
        imu_Ncn_1,      &
        imu_Ncn_2,      &
        imu_hm_1_n,     &
        imu_hm_2_n,     &
        imu_Ncn_1_n,    &
        imu_Ncn_2_n,    &
        isigma_hm_1,    &
        isigma_hm_2,    &
        isigma_Ncn_1,   &
        isigma_Ncn_2,   &
        isigma_hm_1_n,  &
        isigma_hm_2_n,  &
        isigma_Ncn_1_n, &
        isigma_Ncn_2_n

    use stats_variables, only: &
        icorr_whm_1,    & ! Variable(s)
        icorr_whm_2,    &
        icorr_wNcn_1,   &
        icorr_wNcn_2,   &
        icorr_shm_1,    &
        icorr_shm_2,    &
        icorr_sNcn_1,   &
        icorr_sNcn_2,   &
        icorr_thm_1,    &
        icorr_thm_2,    &
        icorr_tNcn_1,   &
        icorr_tNcn_2,   &
        icorr_Ncnhm_1,  &
        icorr_Ncnhm_2,  &
        icorr_hmxhmy_1, &
        icorr_hmxhmy_2

    use stats_variables, only: &
        icorr_whm_1_n,    & ! Variable(s)
        icorr_whm_2_n,    &
        icorr_wNcn_1_n,   &
        icorr_wNcn_2_n,   &
        icorr_shm_1_n,    &
        icorr_shm_2_n,    &
        icorr_sNcn_1_n,   &
        icorr_sNcn_2_n,   &
        icorr_thm_1_n,    &
        icorr_thm_2_n,    &
        icorr_tNcn_1_n,   &
        icorr_tNcn_2_n,   &
        icorr_Ncnhm_1_n,  &
        icorr_Ncnhm_2_n,  &
        icorr_hmxhmy_1_n, &
        icorr_hmxhmy_2_n

    use stats_variables, only: & ! janhft 09/25/12
        icorr_sw,   & ! Variable(s)
        icorr_srr,  &
        icorr_sNr,  &
        icorr_sNcn, &
        icorr_rrNr, &
        icorr_wrr,  &
        icorr_wNr,  &
        icorr_wNcn

    use stats_variables, only: & 
        irel_humidity, &
        irho, & 
        iNcm, &
        iNc_in_cloud, &
        iNc_activated, &
        iNccnm, & 
        isnowslope, & 
        ised_rcm, & 
        irsat, & 
        irsati, & 
        irrainm, & 
        iNrm, & 
        iprecip_rate_zt, & 
        iradht, & 
        iradht_LW, & 
        iradht_SW, & 
        idiam, & 
        imass_ice_cryst, & 
        ircm_icedfs, & 
        iu_T_cm, & 
        im_vol_rad_rain, & 
        im_vol_rad_cloud, & 
        irsnowm, & 
        irgraupelm, & 
        iricem

    use stats_variables, only: & 
        ieff_rad_cloud, &
        ieff_rad_ice, &
        ieff_rad_snow, &
        ieff_rad_rain, &
        ieff_rad_graupel

    use stats_variables, only: & 
        irtm_bt, & 
        irtm_ma, & 
        irtm_ta, & 
        irtm_forcing, & 
        irtm_mc, &
        irtm_sdmp, &
        ircm_mc, &
        ircm_sd_mg_morr, &
        irvm_mc, &
        irtm_mfl, &
        irtm_tacl, & 
        irtm_cl, & 
        irtm_pd, & 
        ithlm_bt, & 
        ithlm_ma, & 
        ithlm_ta, & 
        ithlm_forcing, & 
        ithlm_mc, &
        ithlm_sdmp

    use stats_variables, only: &
        ithlm_mfl, &
        ithlm_tacl, &
        ithlm_cl, &
        iwp3_bt, & 
        iwp3_ma, & 
        iwp3_ta, & 
        iwp3_tp, & 
        iwp3_ac, & 
        iwp3_bp1, & 
        iwp3_bp2, & 
        iwp3_pr1, & 
        iwp3_pr2, & 
        iwp3_dp1, &
        iwp3_4hd, & 
        iwp3_cl

    ! Monotonic flux limiter diagnostic variables
    use stats_variables, only: &
        ithlm_mfl_min, &
        ithlm_mfl_max, &
        irtm_mfl_min, &
        irtm_mfl_max, &
        ithlm_enter_mfl, &
        ithlm_exit_mfl, &
        ithlm_old, &
        ithlm_without_ta, &
        irtm_enter_mfl, &
        irtm_exit_mfl, &
        irtm_old, &
        irtm_without_ta

    use stats_variables, only: & 
        irrainm_bt, & 
        irrainm_ma, & 
        irrainm_ta, &
        irrainm_sd, &
        irrainm_ts, &
        irrainm_sd_morr, &
        irrainm_cond, & 
        irrainm_auto, & 
        irrainm_accr, & 
        irrainm_cond_adj, & 
        irrainm_src_adj, & 
        irrainm_mc, & 
        irrainm_hf

    use stats_variables, only: &
        irrainm_wvhf, & 
        irrainm_cl, & 
        iNrm_bt, & 
        iNrm_ma, & 
        iNrm_ta, & 
        iNrm_sd, & 
        iNrm_ts, & 
        iNrm_cond, & 
        iNrm_auto, & 
        iNrm_cond_adj, & 
        iNrm_src_adj, & 
        iNrm_mc, & 
        iNrm_cl

    use stats_variables, only: & 
        irsnowm_bt, & 
        irsnowm_ma, & 
        irsnowm_sd, &
        irsnowm_sd_morr, &
        irsnowm_ta, &
        irsnowm_mc, & 
        irsnowm_hf, & 
        irsnowm_wvhf, & 
        irsnowm_cl, & 
        irgraupelm_bt, & 
        irgraupelm_ma, & 
        irgraupelm_sd, &
        irgraupelm_sd_morr, &
        irgraupelm_ta, & 
        irgraupelm_mc

    use stats_variables, only: &
        irgraupelm_hf, & 
        irgraupelm_wvhf, & 
        irgraupelm_cl, & 
        iricem_bt, & 
        iricem_ma, & 
        iricem_sd, &
        iricem_sd_mg_morr, &
        iricem_ta, & 
        iricem_mc, & 
        iricem_hf, &
        iricem_wvhf, &
        iricem_cl
 
    use stats_variables, only: & 
        ivm_bt, & 
        ivm_ma, & 
        ivm_gf, & 
        ivm_cf, & 
        ivm_ta, &
        ivm_f, & 
        ivm_sdmp, &
        ivm_ndg, &
        ium_bt, & 
        ium_ma, & 
        ium_gf, & 
        ium_cf, & 
        ium_ta, &
        ium_f, &
        ium_sdmp, &
        ium_ndg

    use stats_variables, only: & 
        imixt_frac, & ! Variable(s) 
        iw1, & 
        iw2, & 
        ivarnce_w1, & 
        ivarnce_w2, & 
        ithl1, & 
        ithl2, & 
        ivarnce_thl1, & 
        ivarnce_thl2, & 
        irt1, & 
        irt2, & 
        ivarnce_rt1, & 
        ivarnce_rt2, & 
        irc1, & 
        irc2, & 
        irsl1, & 
        irsl2, & 
        icloud_frac1, & 
        icloud_frac2

    use stats_variables, only: &
        is1, & 
        is2, & 
        istdev_s1, & 
        istdev_s2, &
        isp2,  &
        istdev_t1, &
        istdev_t2, &
        icovar_st_1, &
        icovar_st_2, &
        icorr_st_1, &
        icorr_st_2, &
        irrtthl, &
        icrt1, &
        icrt2, &
        icthl1, &
        icthl2


    use stats_variables, only: & 
        iwp2_zt, & 
        ithlp2_zt, & 
        iwpthlp_zt, & 
        iwprtp_zt, & 
        irtp2_zt, & 
        irtpthlp_zt, &
        iup2_zt, &
        ivp2_zt, &
        iupwp_zt, &
        ivpwp_zt

    use stats_variables, only: & 
        ihmp2_zt
 
    use stats_variables, only: & 
        zt, & 
        isclrm, & 
        isclrm_f, & 
        iedsclrm, & 
        iedsclrm_f

    use stats_variables, only: & 
        iNsnowm, & ! Variable(s)
        iNrm, &
        iNgraupelm, &
        iNim, & 
        iNsnowm_bt, &
        iNsnowm_mc, &
        iNsnowm_ma, &
        iNsnowm_ta, &
        iNsnowm_sd, &
        iNsnowm_cl, &
        iNgraupelm_bt, &
        iNgraupelm_mc, &
        iNgraupelm_ma, &
        iNgraupelm_ta, &
        iNgraupelm_sd, &
        iNgraupelm_cl, &
        iNim_bt, &
        iNim_mc, &
        iNim_ma, &
        iNim_ta, &
        iNim_sd, &
        iNim_cl

    use stats_variables, only: & 
        iNcm_bt, &
        iNcm_mc, &
        iNcm_ma, &
        iNcm_ta, &
        iNcm_cl, &
        iNcm_act

    use stats_variables, only: &
        iw_KK_evap_covar_zt,   &
        irt_KK_evap_covar_zt,  &
        ithl_KK_evap_covar_zt, &
        iw_KK_auto_covar_zt,   &
        irt_KK_auto_covar_zt,  &
        ithl_KK_auto_covar_zt, &
        iw_KK_accr_covar_zt,   &
        irt_KK_accr_covar_zt,  &
        ithl_KK_accr_covar_zt, &
        irr_KK_mvr_covar_zt,   &
        iNr_KK_mvr_covar_zt,   &
        iKK_mvr_variance_zt

    use stats_variables, only: &
        iC11_Skw_fnc, & ! Variable(s)
        is_mellor, &
        iwp3_on_wp2_zt, &
        ia3_coef_zt
      
    use stats_variables, only: &
        iLscale_pert_1, & ! Variable(s)
        iLscale_pert_2

    use stats_variables, only: &
        iPSMLT,  & ! Variable(s)
        iEVPMS,  &
        iPRACS,  &
        iEVPMG,  &
        iPRACG,  &
        iPGMLT,  &
        iMNUCCC, &
        iPSACWS, &
        iPSACWI, &
        iQMULTS, &
        iQMULTG, &
        iPSACWG, &
        iPGSACW, &
        iPRD,    &
        iPRCI,   &
        iPRAI,   &
        iQMULTR, &
        iQMULTRG,&
        iMNUCCD, &
        iPRACI,  &
        iPRACIS, &
        iEPRD,   &
        iMNUCCR, &
        iPIACR,  &
        iPIACRS, &
        iPGRACS, &
        iPRDS,   &
        iEPRDS,  &
        iPSACR,  &
        iPRDG,   &
        iEPRDG

    use stats_variables, only: &
        iNGSTEN, & ! Lots of variable(s)
        iNRSTEN, &
        iNISTEN, &
        iNSSTEN, &
        iNCSTEN, &
        iNPRC1,  &
        iNRAGG,  &
        iNPRACG, &
        iNSUBR,  &
        iNSMLTR, &
        iNGMLTR, &
        iNPRACS, &
        iNNUCCR, &
        iNIACR,  &
        iNIACRS, &
        iNGRACS, &    
        iNSMLTS, &
        iNSAGG,  &
        iNPRCI, &
        iNSCNG, &
        iNSUBS, &
        iPRC, &
        iPRA, &
        iPRE

    use stats_variables, only: &
        iPCC, &
        iNNUCCC, &
        iNPSACWS, &
        iNPRA, &
        iNPRC, &
        iNPSACWI, &
        iNPSACWG, &
        iNPRAI, &
        iNMULTS, &
        iNMULTG, &
        iNMULTR, &
        iNMULTRG, &
        iNNUCCD, &
        iNSUBI, &
        iNGMLTG, &
        iNSUBG, &
        iNACT, &
        iSIZEFIX_NR, &
        iSIZEFIX_NC, &
        iSIZEFIX_NI, &
        iSIZEFIX_NS, &
        iSIZEFIX_NG, &
        iNEGFIX_NR, &
        iNEGFIX_NC, &
        iNEGFIX_NI, &
        iNEGFIX_NS, &
        iNEGFIX_NG, &
        iNIM_MORR_CL, &
        iQC_INST, &
        iQR_INST, &
        iQI_INST, &
        iQS_INST, & 
        iQG_INST, &
        iNC_INST, &
        iNR_INST, &
        iNI_INST, &
        iNS_INST, & 
        iNG_INST, &
        iT_in_K_mc

    use stats_variables, only: &
        iwp2hmp, & ! Variable(s)
        icloud_frac_refined, &
        ircm_refined, &
        ihl_on_Cp_residual, &
        iqto_residual

    use stats_type, only: & 
        stat_assign ! Procedure

    use parameters_model, only: &
        hydromet_dim, & ! Variable(s)
        sclr_dim,     &
        edsclr_dim

    use parameters_microphys, only: &
        hydromet_list  ! Variable(s)

    use array_index, only: &
        l_mix_rat_hm  ! Variable(s)

    implicit none

    ! External
    intrinsic :: trim

    ! Local Constants

    ! Input Variable
    character(len= * ), dimension(nvarmax_zt), intent(in) :: vars_zt

    ! Input / Output Variable        
    logical, intent(inout) :: l_error

    ! Local Varables
    integer :: tot_zt_loops

    integer :: i, j, k

    integer :: hm_idx, hmx_idx, hmy_idx

    character(len=10) :: hm_type, hmx_type, hmy_type

    logical :: l_found

    character(len=50) :: sclr_idx

    ! The default initialization for array indices for zt is zero (see module
    ! stats_variables)

    ! Allocate and initialize hydrometeor statistical variables.
    allocate( ihm1(1:hydromet_dim) )
    allocate( ihm2(1:hydromet_dim) )
    allocate( imu_hm_1(1:hydromet_dim) )
    allocate( imu_hm_2(1:hydromet_dim) )
    allocate( imu_hm_1_n(1:hydromet_dim) )
    allocate( imu_hm_2_n(1:hydromet_dim) )
    allocate( isigma_hm_1(1:hydromet_dim) )
    allocate( isigma_hm_2(1:hydromet_dim) )
    allocate( isigma_hm_1_n(1:hydromet_dim) )
    allocate( isigma_hm_2_n(1:hydromet_dim) )

    allocate( icorr_whm_1(1:hydromet_dim) )
    allocate( icorr_whm_2(1:hydromet_dim) )
    allocate( icorr_shm_1(1:hydromet_dim) )
    allocate( icorr_shm_2(1:hydromet_dim) )
    allocate( icorr_thm_1(1:hydromet_dim) )
    allocate( icorr_thm_2(1:hydromet_dim) )
    allocate( icorr_Ncnhm_1(1:hydromet_dim) )
    allocate( icorr_Ncnhm_2(1:hydromet_dim) )
    allocate( icorr_hmxhmy_1(1:hydromet_dim,1:hydromet_dim) )
    allocate( icorr_hmxhmy_2(1:hydromet_dim,1:hydromet_dim) )

    allocate( icorr_whm_1_n(1:hydromet_dim) )
    allocate( icorr_whm_2_n(1:hydromet_dim) )
    allocate( icorr_shm_1_n(1:hydromet_dim) )
    allocate( icorr_shm_2_n(1:hydromet_dim) )
    allocate( icorr_thm_1_n(1:hydromet_dim) )
    allocate( icorr_thm_2_n(1:hydromet_dim) )
    allocate( icorr_Ncnhm_1_n(1:hydromet_dim) )
    allocate( icorr_Ncnhm_2_n(1:hydromet_dim) )
    allocate( icorr_hmxhmy_1_n(1:hydromet_dim,1:hydromet_dim) )
    allocate( icorr_hmxhmy_2_n(1:hydromet_dim,1:hydromet_dim) )

    allocate( ihmp2_zt(1:hydromet_dim) )

    ihm1(:) = 0
    ihm2(:) = 0
    imu_hm_1(:) = 0
    imu_hm_2(:) = 0
    imu_hm_1_n(:) = 0
    imu_hm_2_n(:) = 0
    isigma_hm_1(:) = 0
    isigma_hm_2(:) = 0
    isigma_hm_1_n(:) = 0
    isigma_hm_2_n(:) = 0

    icorr_whm_1(:) = 0
    icorr_whm_2(:) = 0
    icorr_shm_1(:) = 0
    icorr_shm_2(:) = 0
    icorr_thm_1(:) = 0
    icorr_thm_2(:) = 0
    icorr_Ncnhm_1(:) = 0
    icorr_Ncnhm_2(:) = 0
    icorr_hmxhmy_1(:,:) = 0
    icorr_hmxhmy_2(:,:) = 0

    icorr_whm_1_n(:) = 0
    icorr_whm_2_n(:) = 0
    icorr_shm_1_n(:) = 0
    icorr_shm_2_n(:) = 0
    icorr_thm_1_n(:) = 0
    icorr_thm_2_n(:) = 0
    icorr_Ncnhm_1_n(:) = 0
    icorr_Ncnhm_2_n(:) = 0
    icorr_hmxhmy_1_n(:,:) = 0
    icorr_hmxhmy_2_n(:,:) = 0

    ihmp2_zt(:) = 0

    ! Allocate and then zero out passive scalar arrays
    allocate( isclrm(1:sclr_dim) )
    allocate( isclrm_f(1:sclr_dim) )

    isclrm(:)     = 0
    isclrm_f(:)   = 0

    allocate( iedsclrm(1:edsclr_dim) )
    allocate( iedsclrm_f(1:edsclr_dim) )

    iedsclrm(:)   = 0
    iedsclrm_f(:) = 0

    ! Assign pointers for statistics variables zt using stat_assign

    tot_zt_loops = zt%nn

    if ( any( vars_zt == "hmi" ) ) then
       ! Correct for number of variables found under "hmi".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "hmi" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "mu_hm_i" ) ) then
       ! Correct for number of variables found under "mu_hm_i".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "mu_hm_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "mu_hm_i_n" ) ) then
       ! Correct for number of variables found under "mu_hm_i_n".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "mu_hm_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "sigma_hm_i" ) ) then
       ! Correct for number of variables found under "sigma_hm_i".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "sigma_hm_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "sigma_hm_i_n" ) ) then
       ! Correct for number of variables found under "sigma_hm_i_n".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "sigma_hm_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_whm_i" ) ) then
       ! Correct for number of variables found under "corr_whm_i".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "corr_whm_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_shm_i" ) ) then
       ! Correct for number of variables found under "corr_shm_i".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "corr_shm_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_thm_i" ) ) then
       ! Correct for number of variables found under "corr_thm_i".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "corr_thm_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_Ncnhm_i" ) ) then
       ! Correct for number of variables found under "corr_Ncnhm_i".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "corr_Ncnhm_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_hmxhmy_i" ) ) then
       ! Correct for number of variables found under "corr_hmxhmy_i".
       ! Subtract 2 (1st PDF component and 2nd PDF component) multipled by the
       ! number of correlations between two hydrometeors, which is found by:
       ! (1/2) * hydromet_dim * ( hydromet_dim - 1 ); from the loop size.
       tot_zt_loops = tot_zt_loops - hydromet_dim * ( hydromet_dim - 1 )
       ! Add 1 for "corr_hmxhmy_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_whm_i_n" ) ) then
       ! Correct for number of variables found under "corr_whm_i_n".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "corr_whm_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_shm_i_n" ) ) then
       ! Correct for number of variables found under "corr_shm_i_n".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "corr_shm_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_thm_i_n" ) ) then
       ! Correct for number of variables found under "corr_thm_i_n".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "corr_thm_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_Ncnhm_i_n" ) ) then
       ! Correct for number of variables found under "corr_Ncnhm_i_n".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "corr_Ncnhm_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_hmxhmy_i_n" ) ) then
       ! Correct for number of variables found under "corr_hmxhmy_i_n".
       ! Subtract 2 (1st PDF component and 2nd PDF component) multipled by the
       ! number of normalized correlations between two hydrometeors, which is
       ! found by:  (1/2) * hydromet_dim * ( hydromet_dim - 1 );
       ! from the loop size.
       tot_zt_loops = tot_zt_loops - hydromet_dim * ( hydromet_dim - 1 )
       ! Add 1 for "corr_hmxhmy_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif

    if ( any( vars_zt == "hmp2_zt" ) ) then
       ! Correct for number of variables found under "hmp2_zt".
       ! Subtract 1 from the loop size for each hydrometeor.
       tot_zt_loops = tot_zt_loops - hydromet_dim
       ! Add 1 for "hmp2_zt" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif

    k = 1

    do i = 1, tot_zt_loops

      select case ( trim( vars_zt(i) ) )
      case ('thlm')
        ithlm = k
        call stat_assign( var_index=ithlm, var_name="thlm", &
             var_description="Liquid water potential temperature (theta_l) [K]", var_units="K", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('T_in_K')
        iT_in_K = k
        call stat_assign( var_index=iT_in_K, var_name="T_in_K", &
             var_description="Absolute temperature [K]", var_units="K", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('thvm')
        ithvm = k
        call stat_assign( var_index=ithvm, var_name="thvm", &
             var_description="Virtual potential temperature [K]", var_units="K", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('rtm')
        irtm = k

        call stat_assign( var_index=irtm, var_name="rtm", &
             var_description="Total (vapor+liquid) water mixing ratio [kg/kg]", &
             var_units="kg/kg", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('rcm')
        ircm = k
        call stat_assign( var_index=ircm, var_name="rcm", &
             var_description="Cloud water mixing ratio [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rfrzm')
        irfrzm = k
        call stat_assign( var_index=irfrzm, var_name="rfrzm", &
             var_description="Total ice phase water mixing ratio [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rvm')
        irvm = k
        call stat_assign( var_index=irvm, var_name="rvm", &
             var_description="Vapor water mixing ratio [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1
      case ('rel_humidity')
        irel_humidity = k
        call stat_assign( var_index=irel_humidity, var_name="rel_humidity", &
             var_description="Relative humidity w.r.t. liquid (range [0,1]) [-]", &
             var_units="[-]", l_silhs=.false., grid_kind=zt )
        k = k + 1
      case ('um')
        ium = k
        call stat_assign( var_index=ium, var_name="um", &
             var_description="East-west (u) wind [m/s]", var_units="m/s", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1
      case ('vm')
        ivm = k
        call stat_assign( var_index=ivm, var_name="vm", &
             var_description="North-south (v) wind [m/s]", var_units="m/s", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1
      case ('wm_zt')
        iwm_zt = k
        call stat_assign( var_index=iwm_zt, var_name="wm", &
             var_description="Vertical (w) wind [m/s]", var_units="m/s", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1
      case ('um_ref')
        ium_ref = k
        call stat_assign( var_index=ium_ref, var_name="um_ref", &
             var_description="reference u wind (m/s) [m/s]", var_units="m/s", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1
      case ('vm_ref')
        ivm_ref = k
        call stat_assign( var_index=ivm_ref, var_name="vm_ref", &
             var_description="reference v wind (m/s) [m/s]", var_units="m/s", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1
      case ('ug')
        iug = k
        call stat_assign( var_index=iug, var_name="ug", &
             var_description="u geostrophic wind [m/s]", var_units="m/s", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1
      case ('vg')
        ivg = k
        call stat_assign( var_index=ivg, var_name="vg", &
             var_description="v geostrophic wind [m/s]", var_units="m/s", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1
      case ('cloud_frac')
        icloud_frac = k
        call stat_assign( var_index=icloud_frac, var_name="cloud_frac", &
             var_description="Cloud fraction (between 0 and 1) [-]", var_units="count", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1
      
      case ('ice_supersat_frac')
        iice_supersat_frac = k
        call stat_assign( var_index=iice_supersat_frac, var_name="ice_supersat_frac", &
             var_description="Ice cloud fraction (between 0 and 1) [-]", var_units="count", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rcm_in_layer')
        ircm_in_layer = k
        call stat_assign( var_index=ircm_in_layer, var_name="rcm_in_layer", &
             var_description="rcm in cloud layer [kg/kg]", var_units="kg/kg", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('rcm_in_cloud')
        ircm_in_cloud = k
        call stat_assign( var_index=ircm_in_cloud, var_name="rcm_in_cloud", &
             var_description="in-cloud value of rcm (for microphysics) [kg/kg]", &
             var_units="kg/kg", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('cloud_cover')
        icloud_cover = k
        call stat_assign( var_index=icloud_cover, var_name="cloud_cover", &
             var_description="Cloud cover (between 0 and 1) [-]", var_units="count", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1
      case ('p_in_Pa')
        ip_in_Pa = k
        call stat_assign( var_index=ip_in_Pa, var_name="p_in_Pa", &
             var_description="Pressure [Pa]", var_units="Pa", l_silhs=.false., grid_kind=zt )
        k = k + 1
      case ('exner')
        iexner = k
        call stat_assign( var_index=iexner, var_name="exner", &
             var_description="Exner function = (p/p0)**(rd/cp) [-]", var_units="count", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1
      case ('rho_ds_zt')
        irho_ds_zt = k
        call stat_assign( var_index=irho_ds_zt, var_name="rho_ds_zt", &
             var_description="Dry, static, base-state density [kg/m^3]", var_units="kg m^{-3}", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1
      case ('thv_ds_zt')
        ithv_ds_zt = k
        call stat_assign( var_index=ithv_ds_zt, var_name="thv_ds_zt", &
             var_description="Dry, base-state theta_v [K]", var_units="K", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1
      case ('Lscale')
        iLscale = k
        call stat_assign( var_index=iLscale, var_name="Lscale", &
             var_description="Mixing length [m]", var_units="m", l_silhs=.false., grid_kind=zt )
        k = k + 1
      case ('thlm_forcing')
        ithlm_forcing = k
        call stat_assign( var_index=ithlm_forcing, var_name="thlm_forcing", &
             var_description="thlm budget: thetal forcing (includes thlm_mc and radht) [K s^{-1}]",&
             var_units="K s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1
      case ('thlm_mc')
        ithlm_mc = k
        call stat_assign( var_index=ithlm_mc, var_name="thlm_mc", &
             var_description="Change in thlm due to microphysics (not in budget) [K s^{-1}]", &
             var_units="K s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1
      case ('rtm_forcing')
        irtm_forcing = k
        call stat_assign( var_index=irtm_forcing, var_name="rtm_forcing", &
             var_description="rtm budget: rt forcing (includes rtm_mc) [kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rtm_mc')
        irtm_mc = k
        call stat_assign( var_index=irtm_mc, var_name="rtm_mc", &
             var_description="Change in rt due to microphysics (not in budget) &
             &[kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rvm_mc')
        irvm_mc = k
        call stat_assign( var_index=irvm_mc, var_name="rvm_mc", &
             var_description="Time tendency of vapor mixing ratio due to microphysics [kg/kg/s]", &
             var_units="kg/(kg s)", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rcm_mc')
        ircm_mc = k
        call stat_assign( var_index=ircm_mc, var_name="rcm_mc", &
             var_description="Time tendency of liquid water mixing ratio due microphysics &
             &[kg/kg/s]", &
             var_units="kg/kg/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rcm_sd_mg_morr')
        ircm_sd_mg_morr = k
        call stat_assign( var_index=ircm_sd_mg_morr, var_name="rcm_sd_mg_morr", &
             var_description="rcm sedimentation when using morrision or MG microphysics &
             &(not in budget, included in rcm_mc) [kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('thlm_mfl_min')
        ithlm_mfl_min = k
        call stat_assign( var_index=ithlm_mfl_min, var_name="thlm_mfl_min", &
             var_description="Minimum allowable thlm [K]", var_units="K", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('thlm_mfl_max')
        ithlm_mfl_max = k
        call stat_assign( var_index=ithlm_mfl_max, var_name="thlm_mfl_max", &
             var_description="Maximum allowable thlm [K]", var_units="K", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('thlm_enter_mfl')
        ithlm_enter_mfl = k
        call stat_assign( var_index=ithlm_enter_mfl, var_name="thlm_enter_mfl", &
             var_description="Thlm before flux-limiter [K]", var_units="K", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('thlm_exit_mfl')
        ithlm_exit_mfl = k
        call stat_assign( var_index=ithlm_exit_mfl, var_name="thlm_exit_mfl", &
             var_description="Thlm exiting flux-limiter [K]", var_units="K", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('thlm_old')
        ithlm_old = k
        call stat_assign( var_index=ithlm_old, var_name="thlm_old", &
             var_description="Thlm at previous timestep [K]", var_units="K", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('thlm_without_ta')
        ithlm_without_ta = k
        call stat_assign( var_index=ithlm_without_ta, var_name="thlm_without_ta", &
             var_description="Thlm without turbulent advection contribution [K]", var_units="K", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rtm_mfl_min')
        irtm_mfl_min = k
        call stat_assign( var_index=irtm_mfl_min, var_name="rtm_mfl_min", &
             var_description="Minimum allowable rtm  [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rtm_mfl_max')
        irtm_mfl_max = k
        call stat_assign( var_index=irtm_mfl_max, var_name="rtm_mfl_max", &
             var_description="Maximum allowable rtm  [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rtm_enter_mfl')
        irtm_enter_mfl = k
        call stat_assign( var_index=irtm_enter_mfl, var_name="rtm_enter_mfl", &
             var_description="Rtm before flux-limiter  [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rtm_exit_mfl')
        irtm_exit_mfl = k
        call stat_assign( var_index=irtm_exit_mfl, var_name="rtm_exit_mfl", &
             var_description="Rtm exiting flux-limiter  [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rtm_old')
        irtm_old = k
        call stat_assign( var_index=irtm_old, var_name="rtm_old", &
             var_description="Rtm at previous timestep  [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rtm_without_ta')
        irtm_without_ta = k
        call stat_assign( var_index=irtm_without_ta, var_name="rtm_without_ta", &
             var_description="Rtm without turbulent advection contribution  [kg/kg]", &
             var_units="kg/kg", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('wp3')
        iwp3 = k
        call stat_assign( var_index=iwp3, var_name="wp3", &
             var_description="w third order moment [m^3/s^3]", var_units="m^3/s^3", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('wpthlp2')
        iwpthlp2 = k
        call stat_assign( var_index=iwpthlp2, var_name="wpthlp2", &
             var_description="w'thl'^2 [(m K^2)/s]", var_units="(m K^2)/s", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('wp2thlp')
        iwp2thlp = k
        call stat_assign( var_index=iwp2thlp, var_name="wp2thlp", &
             var_description="w'^2thl' [(m^2 K)/s^2]", var_units="(m^2 K)/s^2", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('wprtp2')
        iwprtp2 = k
        call stat_assign( var_index=iwprtp2, var_name="wprtp2", &
             var_description="w'rt'^2 [(m kg)/(s kg)]", var_units="(m kg)/(s kg)", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('wp2rtp')
        iwp2rtp = k
        call stat_assign( var_index=iwp2rtp, var_name="wp2rtp", &
             var_description="w'^2rt' [(m^2 kg)/(s^2 kg)]", var_units="(m^2 kg)/(s^2 kg)", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Lscale_up')
        iLscale_up = k
        call stat_assign( var_index=iLscale_up, var_name="Lscale_up", &
             var_description="Upward mixing length [m]", var_units="m", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('Lscale_down')
        iLscale_down = k
        call stat_assign( var_index=iLscale_down, var_name="Lscale_down", &
             var_description="Downward mixing length [m]", var_units="m", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('Lscale_pert_1')
        iLscale_pert_1 = k
        call stat_assign( var_index=iLscale_pert_1, var_name="Lscale_pert_1", &
             var_description="Mixing length using a perturbed value of rtm/thlm [m]", &
             var_units="m", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Lscale_pert_2')
        iLscale_pert_2 = k
        call stat_assign( var_index=iLscale_pert_2, var_name="Lscale_pert_2", &
             var_description="Mixing length using a perturbed value of rtm/thlm [m]", &
             var_units="m", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('tau_zt')
        itau_zt = k
        call stat_assign( var_index=itau_zt, var_name="tau_zt", &
             var_description="Dissipation time [s]", var_units="s", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('Kh_zt')
        iKh_zt = k
        call stat_assign( var_index=iKh_zt, var_name="Kh_zt", &
             var_description="Eddy diffusivity [m^2/s]", var_units="m^2/s", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('wp2thvp')
        iwp2thvp = k
        call stat_assign( var_index=iwp2thvp, var_name="wp2thvp", &
             var_description="w'^2thv' [K m^2/s^2]", var_units="K m^2/s^2", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('wp2rcp')
        iwp2rcp = k
        call stat_assign( var_index=iwp2rcp, var_name="wp2rcp", &
             var_description="w'^2rc' [(m^2 kg)/(s^2 kg)]", var_units="(m^2 kg)/(s^2 kg)", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('wprtpthlp')
        iwprtpthlp = k
        call stat_assign( var_index=iwprtpthlp, var_name="wprtpthlp", &
             var_description="w'rt'thl' [(m kg K)/(s kg)]", var_units="(m kg K)/(s kg)", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('sigma_sqd_w_zt')
        isigma_sqd_w_zt = k
        call stat_assign( var_index=isigma_sqd_w_zt, var_name="sigma_sqd_w_zt", &
             var_description="Nondimensionalized w variance of Gaussian component [-]", &
             var_units="-", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rho')
        irho = k
        call stat_assign( var_index=irho, var_name="rho", var_description="Air density [kg/m^3]", &
             var_units="kg m^{-3}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Ncm')           ! Brian
        iNcm = k
        call stat_assign( var_index=iNcm, var_name="Ncm", &
             var_description="Cloud droplet number concentration [num/kg]", var_units="num/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Nc_in_cloud')
        iNc_in_cloud = k

        call stat_assign( var_index=iNc_in_cloud, var_name="Nc_in_cloud", &
             var_description="In cloud droplet concentration [num/kg]", var_units="num/kg", &
             l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('Nc_activated')
        iNc_activated = k

        call stat_assign( var_index=iNc_activated, var_name="Nc_activated", &
             var_description="Droplets activated by GFDL activation [num/kg]", &
             var_units="num/kg", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('Nccnm')
        iNccnm = k
        call stat_assign( var_index=iNccnm, var_name="Nccnm", &
             var_description="Cloud condensation nuclei concentration (COAMPS/MG) [num/kg]", &
             var_units="num/kg", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Nim')           ! Brian
        iNim = k
        call stat_assign( var_index=iNim, var_name="Nim", &
             var_description="Ice crystal number concentration [num/kg]", var_units="num/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('snowslope')     ! Adam Smith, 22 April 2008
        isnowslope = k
        call stat_assign( var_index=isnowslope, var_name="snowslope", &
             var_description="COAMPS microphysics snow slope parameter [1/m]", var_units="1/m", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Nsnowm')        ! Adam Smith, 22 April 2008
        iNsnowm = k
        call stat_assign( var_index=iNsnowm, var_name="Nsnowm", &
             var_description="Snow particle number concentration [num/kg]", var_units="num/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Ngraupelm')
        iNgraupelm = k
        call stat_assign( var_index=iNgraupelm, var_name="Ngraupelm", &
             var_description="Graupel number concentration  [num/kg]", var_units="num/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('sed_rcm')       ! Brian
        ised_rcm = k
        call stat_assign( var_index=ised_rcm, var_name="sed_rcm", &
             var_description="d(rcm)/dt due to cloud sedimentation [kg / (m^2 s)]", &
             var_units="kg / [m^2 s]", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rsat')           ! Brian
        irsat = k
        call stat_assign( var_index=irsat, var_name="rsat", &
             var_description="Saturation mixing ratio over liquid [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rsati')
        irsati = k
        call stat_assign( var_index=irsati, var_name="rsati", &
             var_description="Saturation mixing ratio over ice [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rrainm')           ! Brian
        irrainm = k
        call stat_assign( var_index=irrainm, var_name="rrainm", &
             var_description="Rain water mixing ratio [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rsnowm')
        irsnowm = k
        call stat_assign( var_index=irsnowm, var_name="rsnowm", &
             var_description="Snow water mixing ratio [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('ricem')
        iricem = k
        call stat_assign( var_index=iricem, var_name="ricem", &
             var_description="Pristine ice water mixing ratio [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rgraupelm')
        irgraupelm = k
        call stat_assign( var_index=irgraupelm, var_name="rgraupelm", &
             var_description="Graupel water mixing ratio [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Nrm')           ! Brian
        iNrm = k
        call stat_assign( var_index=iNrm, var_name="Nrm", &
             var_description="Rain drop number concentration [num/kg]", var_units="num/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('m_vol_rad_rain')  ! Brian
        im_vol_rad_rain = k
        call stat_assign( var_index=im_vol_rad_rain, var_name="mvrr", &
             var_description="Rain drop mean volume radius [m]", var_units="m", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('m_vol_rad_cloud')
        im_vol_rad_cloud = k
        call stat_assign( var_index=im_vol_rad_cloud, var_name="m_vol_rad_cloud", &
             var_description="Cloud drop mean volume radius [m]", var_units="m", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('eff_rad_cloud')
        ieff_rad_cloud = k
        call stat_assign( var_index=ieff_rad_cloud, var_name="eff_rad_cloud", &
             var_description="Cloud drop effective volume radius [microns]", var_units="microns", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('eff_rad_ice')
        ieff_rad_ice = k

        call stat_assign( var_index=ieff_rad_ice, var_name="eff_rad_ice", &
             var_description="Ice effective volume radius [microns]", var_units="microns", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('eff_rad_snow')
        ieff_rad_snow = k
        call stat_assign( var_index=ieff_rad_snow, var_name="eff_rad_snow", &
             var_description="Snow effective volume radius [microns]", var_units="microns", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('eff_rad_rain')
        ieff_rad_rain = k
        call stat_assign( var_index=ieff_rad_rain, var_name="eff_rad_rain", &
             var_description="Rain drop effective volume radius [microns]", var_units="microns", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('eff_rad_graupel')
        ieff_rad_graupel = k
        call stat_assign( var_index=ieff_rad_graupel, var_name="eff_rad_graupel", &
             var_description="Graupel effective volume radius [microns]", var_units="microns", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('precip_rate_zt')     ! Brian
        iprecip_rate_zt = k

        call stat_assign( var_index=iprecip_rate_zt, var_name="precip_rate_zt", &
             var_description="Rain rate [mm/day]", var_units="mm/day", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('radht')
        iradht = k

        call stat_assign( var_index=iradht, var_name="radht", &
             var_description="Total (sw+lw) radiative heating rate [K/s]", var_units="K/s", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('radht_LW')
        iradht_LW = k

        call stat_assign( var_index=iradht_LW, var_name="radht_LW", &
             var_description="Long-wave radiative heating rate [K/s]", var_units="K/s", &
             l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('radht_SW')
        iradht_SW = k
        call stat_assign( var_index=iradht_SW, var_name="radht_SW", &
             var_description="Short-wave radiative heating rate [K/s]", var_units="K/s", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('diam')
        idiam = k

        call stat_assign( var_index=idiam, var_name="diam", &
             var_description="Ice crystal diameter [m]", var_units="m", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('mass_ice_cryst')
        imass_ice_cryst = k
        call stat_assign( var_index=imass_ice_cryst, var_name="mass_ice_cryst", &
             var_description="Mass of a single ice crystal [kg]", var_units="kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rcm_icedfs')

        ircm_icedfs = k
        call stat_assign( var_index=ircm_icedfs, var_name="rcm_icedfs", &
             var_description="Change in liquid due to ice [kg/kg/s]", var_units="kg/kg/s", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('u_T_cm')
        iu_T_cm = k
        call stat_assign( var_index=iu_T_cm, var_name="u_T_cm", &
             var_description="Ice crystal fallspeed [cm s^{-1}]", var_units="cm s^{-1}", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rtm_bt')
        irtm_bt = k

        call stat_assign( var_index=irtm_bt, var_name="rtm_bt", &
             var_description="rtm budget: rtm time tendency [kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rtm_ma')
        irtm_ma = k

        call stat_assign( var_index=irtm_ma, var_name="rtm_ma", &
             var_description="rtm budget: rtm vertical mean advection [kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rtm_ta')
        irtm_ta = k

        call stat_assign( var_index=irtm_ta, var_name="rtm_ta", &
             var_description="rtm budget: rtm turbulent advection [kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rtm_mfl')
        irtm_mfl = k

        call stat_assign( var_index=irtm_mfl, var_name="rtm_mfl", &
             var_description="rtm budget: rtm correction due to monotonic flux limiter &
             &[kg kg^{-1} s^{-1}]", var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rtm_tacl')
        irtm_tacl = k

        call stat_assign( var_index=irtm_tacl, var_name="rtm_tacl", &
             var_description="rtm budget: rtm correction due to ta term (wprtp) clipping &
             &[kg kg^{-1} s^{-1}]", var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('rtm_cl')
        irtm_cl = k

        call stat_assign( var_index=irtm_cl, var_name="rtm_cl", &
             var_description="rtm budget: rtm clipping [kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )

        k = k + 1
      case ('rtm_sdmp')
        irtm_sdmp = k

        call stat_assign( var_index=irtm_sdmp, var_name="rtm_sdmp", &
             var_description="rtm budget: rtm correction due to sponge damping &
             &[kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1


      case ('rtm_pd')
        irtm_pd = k

        call stat_assign( var_index=irtm_pd, var_name="rtm_pd", &
             var_description="rtm budget: rtm positive definite adjustment [kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('thlm_bt')
        ithlm_bt = k

        call stat_assign( var_index=ithlm_bt, var_name="thlm_bt", &
             var_description="thlm budget: thlm time tendency [K s^{-1}]", var_units="K s^{-1}", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('thlm_ma')
        ithlm_ma = k

        call stat_assign( var_index=ithlm_ma, var_name="thlm_ma", &
             var_description="thlm budget: thlm vertical mean advection [K s^{-1}]", &
             var_units="K s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('thlm_sdmp')
        ithlm_sdmp = k

        call stat_assign( var_index=ithlm_sdmp, var_name="thlm_sdmp", &
             var_description="thlm budget: thlm correction due to sponge damping [K s^{-1}]", &
             var_units="K s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1


      case ('thlm_ta')
        ithlm_ta = k

        call stat_assign( var_index=ithlm_ta, var_name="thlm_ta", &
             var_description="thlm budget: thlm turbulent advection [K s^{-1}]", &
             var_units="K s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('thlm_mfl')
        ithlm_mfl = k

        call stat_assign( var_index=ithlm_mfl, var_name="thlm_mfl", &
             var_description="thlm budget: thlm correction due to monotonic flux limiter &
             &[K s^{-1}]", &
             var_units="K s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('thlm_tacl')
        ithlm_tacl = k

        call stat_assign( var_index=ithlm_tacl, var_name="thlm_tacl", &
             var_description="thlm budget: thlm correction due to ta term (wpthlp) clipping &
             &[K s^{-1}]", &
             var_units="K s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('thlm_cl')
        ithlm_cl = k

        call stat_assign( var_index=ithlm_cl, var_name="thlm_cl", &
             var_description="thlm budget: thlm_cl [K s^{-1}]", var_units="K s^{-1}", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('wp3_bt')
        iwp3_bt = k

        call stat_assign( var_index=iwp3_bt, var_name="wp3_bt", &
             var_description="wp3 budget: wp3 time tendency [m^{3} s^{-4}]", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('wp3_ma')
        iwp3_ma = k

        call stat_assign( var_index=iwp3_ma, var_name="wp3_ma", &
             var_description="wp3 budget: wp3 vertical mean advection [m^{3} s^{-4}]", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('wp3_ta')
        iwp3_ta = k

        call stat_assign( var_index=iwp3_ta, var_name="wp3_ta", &
             var_description="wp3 budget: wp3 turbulent advection [m^{3} s^{-4}]", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('wp3_tp')
        iwp3_tp = k
        call stat_assign( var_index=iwp3_tp, var_name="wp3_tp", &
             var_description="wp3 budget: wp3 turbulent transport [m^{3} s^{-4}]", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('wp3_ac')
        iwp3_ac = k
        call stat_assign( var_index=iwp3_ac, var_name="wp3_ac", &
             var_description="wp3 budget: wp3 accumulation term [m^{3} s^{-4}]", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('wp3_bp1')
        iwp3_bp1 = k
        call stat_assign( var_index=iwp3_bp1, var_name="wp3_bp1", &
             var_description="wp3 budget: wp3 buoyancy production [m^{3} s^{-4}]", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('wp3_bp2')
        iwp3_bp2 = k
        call stat_assign( var_index=iwp3_bp2, var_name="wp3_bp2", &
             var_description="wp3 budget: wp3 2nd buoyancy production term [m^{3} s^{-4}]", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('wp3_pr1')
        iwp3_pr1 = k
        call stat_assign( var_index=iwp3_pr1, var_name="wp3_pr1", &
             var_description="wp3 budget: wp3 pressure term 1 [m^{3} s^{-4}]", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('wp3_pr2')
        iwp3_pr2 = k
        call stat_assign( var_index=iwp3_pr2, var_name="wp3_pr2", &
             var_description="wp3 budget: wp3 pressure term 2 [m^{3} s^{-4}]", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('wp3_dp1')
        iwp3_dp1 = k
        call stat_assign( var_index=iwp3_dp1, var_name="wp3_dp1", &
             var_description="wp3 budget: wp3 dissipation term 1 [m^{3} s^{-4}]", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('wp3_4hd')
        iwp3_4hd = k
        call stat_assign( var_index=iwp3_4hd, var_name="wp3_4hd", &
             var_description="wp3 budget: wp3 4th-order hyper-diffusion [m^{3} s^{-4}]", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('wp3_cl')
        iwp3_cl = k
        call stat_assign( var_index=iwp3_cl, var_name="wp3_cl", &
             var_description="wp3 budget: wp3 clipping term [m^{3} s^{-4}]", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rrainm_bt')
        irrainm_bt = k
        call stat_assign( var_index=irrainm_bt, var_name="rrainm_bt", &
             var_description="rrainm budget: rrainm time tendency [kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rrainm_ma')
        irrainm_ma = k

        call stat_assign( var_index=irrainm_ma, var_name="rrainm_ma", &
             var_description="rrainm budget: rrainm vertical mean advection [kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rrainm_sd')
        irrainm_sd = k

        call stat_assign( var_index=irrainm_sd, var_name="rrainm_sd", &
             var_description="rrainm budget: rrainm sedimentation [kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rrainm_ts')
        irrainm_ts = k

        call stat_assign( var_index=irrainm_ts, var_name="rrainm_ts", &
             var_description="rrainm budget: rrainm turbulent sedimentation [kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rrainm_sd_morr')
        irrainm_sd_morr = k

        call stat_assign( var_index=irrainm_sd_morr, var_name="rrainm_sd_morr", &
             var_description="rrainm sedimentation when using morrision microphysics &
             &(not in budget, included in rrainm_mc) [kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('rrainm_ta')
        irrainm_ta = k

        call stat_assign( var_index=irrainm_ta, var_name="rrainm_ta", &
             var_description="rrainm budget: rrainm turbulent advection [kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rrainm_cond')
        irrainm_cond = k

        call stat_assign( var_index=irrainm_cond, var_name="rrainm_cond", &
             var_description="rrainm evaporation rate [kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rrainm_auto')
        irrainm_auto = k

        call stat_assign( var_index=irrainm_auto, var_name="rrainm_auto", &
             var_description="rrainm autoconversion rate [kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rrainm_accr')
        irrainm_accr = k
        call stat_assign( var_index=irrainm_accr, var_name="rrainm_accr", &
             var_description="rrainm accretion rate [kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rrainm_cond_adj')
        irrainm_cond_adj = k

        call stat_assign( var_index=irrainm_cond_adj, var_name="rrainm_cond_adj", &
             var_description="rrainm evaporation adjustment due to over-evaporation &
             &[kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rrainm_src_adj')
        irrainm_src_adj = k

        call stat_assign( var_index=irrainm_src_adj, var_name="rrainm_src_adj", &
             var_description="rrainm source term adjustment due to over-depletion &
             &[kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rrainm_hf')
        irrainm_hf = k
        call stat_assign( var_index=irrainm_hf, var_name="rrainm_hf", &
             var_description="rrainm budget: rrainm hole-filling term [kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rrainm_wvhf')
        irrainm_wvhf = k
        call stat_assign( var_index=irrainm_wvhf, var_name="rrainm_wvhf", &
             var_description="rrainm budget: rrainm water vapor hole-filling term &
             &[kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rrainm_cl')
        irrainm_cl = k
        call stat_assign( var_index=irrainm_cl, var_name="rrainm_cl", &
             var_description="rrainm budget: rrainm clipping term [kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rrainm_mc')
        irrainm_mc = k

        call stat_assign( var_index=irrainm_mc, var_name="rrainm_mc", &
             var_description="rrainm budget: Change in rrainm due to microphysics &
             &[kg kg^{-1} s^{-1}]", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Nrm_bt')
        iNrm_bt = k
        call stat_assign( var_index=iNrm_bt, var_name="Nrm_bt", &
             var_description="Nrm budget: Nrm time tendency [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('Nrm_ma')
        iNrm_ma = k

        call stat_assign( var_index=iNrm_ma, var_name="Nrm_ma", &
             var_description="Nrm budget: Nrm vertical mean advection [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Nrm_sd')
        iNrm_sd = k

        call stat_assign( var_index=iNrm_sd, var_name="Nrm_sd", &
             var_description="Nrm budget: Nrm sedimentation [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('Nrm_ts')
        iNrm_ts = k

        call stat_assign( var_index=iNrm_ts, var_name="Nrm_ts", &
             var_description="Nrm budget: Nrm turbulent sedimentation [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Nrm_ta')
        iNrm_ta = k
        call stat_assign( var_index=iNrm_ta, var_name="Nrm_ta", &
             var_description="Nrm budget: Nrm turbulent advection [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('Nrm_cond')
        iNrm_cond = k

        call stat_assign( var_index=iNrm_cond, var_name="Nrm_cond", &
             var_description="Nrm evaporation rate [(num/kg)/s]", var_units="(num/kg)/s", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Nrm_auto')
        iNrm_auto = k

        call stat_assign( var_index=iNrm_auto, var_name="Nrm_auto", &
             var_description="Nrm autoconversion rate [(num/kg)/s]", var_units="(num/kg)/s", &
             l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('Nrm_cond_adj')
        iNrm_cond_adj = k

        call stat_assign( var_index=iNrm_cond_adj, var_name="Nrm_cond_adj", &
             var_description="Nrm evaporation adjustment due to over-evaporation [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Nrm_src_adj')
        iNrm_src_adj = k

        call stat_assign( var_index=iNrm_src_adj, var_name="Nrm_src_adj", &
             var_description="Nrm source term adjustment due to over-depletion [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Nrm_cl')
        iNrm_cl = k
        call stat_assign( var_index=iNrm_cl, var_name="Nrm_cl", &
             var_description="Nrm budget: Nrm clipping term [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Nrm_mc')
        iNrm_mc = k
        call stat_assign( var_index=iNrm_mc, var_name="Nrm_mc", &
             var_description="Nrm budget: Change in Nrm due to microphysics (Not in budget) &
             &[(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rsnowm_bt')
        irsnowm_bt = k
        call stat_assign( var_index=irsnowm_bt, var_name="rsnowm_bt", &
             var_description="rsnowm budget: rsnowm time tendency [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('rsnowm_ma')
        irsnowm_ma = k

        call stat_assign( var_index=irsnowm_ma, var_name="rsnowm_ma", &
             var_description="rsnowm budget: rsnowm vertical mean advection [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rsnowm_sd')
        irsnowm_sd = k
        call stat_assign( var_index=irsnowm_sd, var_name="rsnowm_sd", &
             var_description="rsnowm budget: rsnowm sedimentation [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rsnowm_sd_morr')
        irsnowm_sd_morr = k
        call stat_assign( var_index=irsnowm_sd_morr, var_name="rsnowm_sd_morr", &
             var_description="rsnowm sedimentation when using morrison microphysics &
             &(Not in budget, included in rsnowm_mc) [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('rsnowm_ta')
        irsnowm_ta = k

        call stat_assign( var_index=irsnowm_ta, var_name="rsnowm_ta", &
             var_description="rsnowm budget: rsnowm turbulent advection [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rsnowm_mc')
        irsnowm_mc = k

        call stat_assign( var_index=irsnowm_mc, var_name="rsnowm_mc", &
             var_description="rsnowm budget: Change in rsnowm due to microphysics [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rsnowm_hf')
        irsnowm_hf = k

        call stat_assign( var_index=irsnowm_hf, var_name="rsnowm_hf", &
             var_description="rsnowm budget: rsnowm hole-filling term [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rsnowm_wvhf')
        irsnowm_wvhf = k

        call stat_assign( var_index=irsnowm_wvhf, var_name="rsnowm_wvhf", &
             var_description="rsnowm budget: rsnowm water vapor hole-filling term [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rsnowm_cl')
        irsnowm_cl = k

        call stat_assign( var_index=irsnowm_cl, var_name="rsnowm_cl", &
             var_description="rsnowm budget: rsnowm clipping term [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Nsnowm_bt')
        iNsnowm_bt = k
        call stat_assign( var_index=iNsnowm_bt, var_name="Nsnowm_bt", &
             var_description="Nsnowm budget: [(num/kg)/s]", var_units="(num/kg)/s", &
             l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('Nsnowm_ma')
        iNsnowm_ma = k

        call stat_assign( var_index=iNsnowm_ma, var_name="Nsnowm_ma", &
             var_description="Nsnowm budget: Nsnowm mean advection [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Nsnowm_sd')
        iNsnowm_sd = k

        call stat_assign( var_index=iNsnowm_sd, var_name="Nsnowm_sd", &
             var_description="Nsnowm budget: Nsnowm sedimentation [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('Nsnowm_ta')
        iNsnowm_ta = k
        call stat_assign( var_index=iNsnowm_ta, var_name="Nsnowm_ta", &
             var_description="Nsnowm budget: Nsnowm turbulent advection [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('Nsnowm_mc')
        iNsnowm_mc = k
        call stat_assign( var_index=iNsnowm_mc, var_name="Nsnowm_mc", &
             var_description="Nsnowm budget: Nsnowm microphysics [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('Nsnowm_cl')
        iNsnowm_cl = k

        call stat_assign( var_index=iNsnowm_cl, var_name="Nsnowm_cl", &
             var_description="Nsnowm budget: Nsnowm clipping term [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('ricem_bt')
        iricem_bt = k

        call stat_assign( var_index=iricem_bt, var_name="ricem_bt", &
             var_description="ricem budget: ricem time tendency [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('ricem_ma')
        iricem_ma = k

        call stat_assign( var_index=iricem_ma, var_name="ricem_ma", &
             var_description="ricem budget: ricem vertical mean advection [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('ricem_sd')
        iricem_sd = k

        call stat_assign( var_index=iricem_sd, var_name="ricem_sd", &
             var_description="ricem budget: ricem sedimentation [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('ricem_sd_mg_morr')
        iricem_sd_mg_morr = k

        call stat_assign( var_index=iricem_sd_mg_morr, var_name="ricem_sd_mg_morr", &
             var_description="ricem sedimentation when using morrison or MG microphysics &
             &(not in budget, included in ricem_mc) [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('ricem_ta')
        iricem_ta = k

        call stat_assign( var_index=iricem_ta, var_name="ricem_ta", &
             var_description="ricem budget: ricem turbulent advection [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('ricem_mc')
        iricem_mc = k

        call stat_assign( var_index=iricem_mc, var_name="ricem_mc", &
             var_description="ricem budget: Change in ricem due to microphysics [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('ricem_hf')
        iricem_hf = k

        call stat_assign( var_index=iricem_hf, var_name="ricem_hf", &
             var_description="ricem budget: ricem hole-filling term [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('ricem_wvhf')
        iricem_wvhf = k

        call stat_assign( var_index=iricem_wvhf, var_name="ricem_wvhf", &
             var_description="ricem budget: ricem water vapor hole-filling term [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('ricem_cl')
        iricem_cl = k

        call stat_assign( var_index=iricem_cl, var_name="ricem_cl", &
             var_description="ricem budget: ricem clipping term [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rgraupelm_bt')
        irgraupelm_bt = k

        call stat_assign( var_index=irgraupelm_bt, var_name="rgraupelm_bt", &
             var_description="rgraupelm budget: rgraupelm time tendency [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rgraupelm_ma')
        irgraupelm_ma = k

        call stat_assign( var_index=irgraupelm_ma, var_name="rgraupelm_ma", &
             var_description="rgraupelm budget: rgraupelm vertical mean advection [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rgraupelm_sd')
        irgraupelm_sd = k

        call stat_assign( var_index=irgraupelm_sd, var_name="rgraupelm_sd", &
             var_description="rgraupelm budget: rgraupelm sedimentation [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rgraupelm_sd_morr')
        irgraupelm_sd_morr = k

        call stat_assign( var_index=irgraupelm_sd_morr, var_name="rgraupelm_sd_morr", &
             var_description="rgraupelm sedimentation when using morrison microphysics &
             &(not in budget, included in rgraupelm_mc) [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('rgraupelm_ta')
        irgraupelm_ta = k

        call stat_assign( var_index=irgraupelm_ta, var_name="rgraupelm_ta", &
             var_description="rgraupelm budget: rgraupelm turbulent advection [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rgraupelm_mc')
        irgraupelm_mc = k

        call stat_assign( var_index=irgraupelm_mc, var_name="rgraupelm_mc", &
             var_description="rgraupelm budget: Change in rgraupelm due to microphysics &
             &[(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rgraupelm_hf')
        irgraupelm_hf = k

        call stat_assign( var_index=irgraupelm_hf, var_name="rgraupelm_hf", &
             var_description="rgraupelm budget: rgraupelm hole-filling term [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rgraupelm_wvhf')
        irgraupelm_wvhf = k

        call stat_assign( var_index=irgraupelm_wvhf, var_name="rgraupelm_wvhf", &
             var_description="rgraupelm budget: rgraupelm water vapor hole-filling term &
             &[(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rgraupelm_cl')
        irgraupelm_cl = k

        call stat_assign( var_index=irgraupelm_cl, var_name="rgraupelm_cl", &
             var_description="rgraupelm budget: rgraupelm clipping term [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Ngraupelm_bt')
        iNgraupelm_bt = k
        call stat_assign( var_index=iNgraupelm_bt, var_name="Ngraupelm_bt", &
             var_description="Ngraupelm budget: [(num/kg)/s]", var_units="(num/kg)/s", &
             l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('Ngraupelm_ma')
        iNgraupelm_ma = k

        call stat_assign( var_index=iNgraupelm_ma, var_name="Ngraupelm_ma", &
             var_description="Ngraupelm budget: Ngraupelm mean advection [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Ngraupelm_sd')
        iNgraupelm_sd = k

        call stat_assign( var_index=iNgraupelm_sd, var_name="Ngraupelm_sd", &
             var_description="Ngraupelm budget: Ngraupelm sedimentation [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('Ngraupelm_ta')
        iNgraupelm_ta = k
        call stat_assign( var_index=iNgraupelm_ta, var_name="Ngraupelm_ta", &
             var_description="Ngraupelm budget: Ngraupelm turbulent advection [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('Ngraupelm_mc')
        iNgraupelm_mc = k

        call stat_assign( var_index=iNgraupelm_mc, var_name="Ngraupelm_mc", &
             var_description="Ngraupelm budget: Ngraupelm microphysics term [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Ngraupelm_cl')
        iNgraupelm_cl = k

        call stat_assign( var_index=iNgraupelm_cl, var_name="Ngraupelm_cl", &
             var_description="Ngraupelm budget: Ngraupelm clipping term [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Nim_bt')
        iNim_bt = k
        call stat_assign( var_index=iNim_bt, var_name="Nim_bt", &
             var_description="Nim budget: [(num/kg)/s]", var_units="(num/kg)/s", l_silhs=.false., &
             grid_kind=zt )

        k = k + 1

      case ('Nim_ma')
        iNim_ma = k

        call stat_assign( var_index=iNim_ma, var_name="Nim_ma", &
             var_description="Nim budget: Nim mean advection [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Nim_sd')
        iNim_sd = k

        call stat_assign( var_index=iNim_sd, var_name="Nim_sd", &
             var_description="Nim budget: Nim sedimentation [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('Nim_ta')
        iNim_ta = k
        call stat_assign( var_index=iNim_ta, var_name="Nim_ta", &
             var_description="Nim budget: Nim turbulent advection [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('Nim_mc')
        iNim_mc = k

        call stat_assign( var_index=iNim_mc, var_name="Nim_mc", &
             var_description="Nim budget: Nim microphysics term [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Nim_cl')
        iNim_cl = k

        call stat_assign( var_index=iNim_cl, var_name="Nim_cl", &
             var_description="Nim budget: Nim clipping term [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Ncm_bt')
        iNcm_bt = k
        call stat_assign( var_index=iNcm_bt, var_name="Ncm_bt", &
             var_description="Ncm budget: Cloud droplet number concentration budget [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('Ncm_ma')
        iNcm_ma = k

        call stat_assign( var_index=iNcm_ma, var_name="Ncm_ma", &
             var_description="Ncm budget: Ncm vertical mean advection [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Ncm_act')
        iNcm_act = k

        call stat_assign( var_index=iNcm_act, var_name="Ncm_act", &
             var_description="Ncm budget: Change in Ncm due to activation [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('Ncm_ta')
        iNcm_ta = k
        call stat_assign( var_index=iNcm_ta, var_name="Ncm_ta", &
             var_description="Ncm budget: Ncm turbulent advection [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('Ncm_mc')
        iNcm_mc = k

        call stat_assign( var_index=iNcm_mc, var_name="Ncm_mc", &
             var_description="Ncm budget: Change in Ncm due to microphysics [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('Ncm_cl')
        iNcm_cl = k

        call stat_assign( var_index=iNcm_cl, var_name="Ncm_cl", &
             var_description="Ncm budget: Ncm clipping term [(num/kg)/s]", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('PSMLT')
        iPSMLT = k

        call stat_assign( var_index=iPSMLT, var_name="PSMLT", &
             var_description="Freezing of rain to form snow, +rsnowm, -rrainm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('EVPMS')
        iEVPMS = k

        call stat_assign( var_index=iEVPMS, var_name="EVPMS", &
             var_description="Evaporation of melted snow, +rsnowm, -rvm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PRACS')
        iPRACS = k

        call stat_assign( var_index=iPRACS, var_name="PRACS", &
             var_description="Collection of rain by snow, +rsnowm, -rrainm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('EVPMG')
        iEVPMG = k

        call stat_assign( var_index=iEVPMG, var_name="EVPMG", &
             var_description="Evaporation of melted graupel, +rgraupelm, -rvm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PRACG')
        iPRACG = k

        call stat_assign( var_index=iPRACG, var_name="PRACG", &
             var_description="Negative of collection of rain by graupel, +rrainm, -rgraupelm &
             &[(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PGMLT')
        iPGMLT = k

        call stat_assign( var_index=iPGMLT, var_name="PGMLT", &
             var_description="Negative of melting of graupel, +rgraupelm, -rrainm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('MNUCCC')
        iMNUCCC = k

        call stat_assign( var_index=iMNUCCC, var_name="MNUCCC", &
             var_description="Contact freezing of cloud droplets, +ricem, -rcm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PSACWS')
        iPSACWS = k

        call stat_assign( var_index=iPSACWS, var_name="PSACWS", &
             var_description="Collection of cloud water by snow, +rsnowm, -rcm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PSACWI')
        iPSACWI = k

        call stat_assign( var_index=iPSACWI, var_name="PSACWI", &
             var_description="Collection of cloud water by cloud ice, +ricem, -rcm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('QMULTS')
        iQMULTS = k

        call stat_assign( var_index=iQMULTS, var_name="QMULTS", &
             var_description="Splintering from cloud droplets accreted onto snow, +ricem, -rcm &
             &[(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('QMULTG')
        iQMULTG = k

        call stat_assign( var_index=iQMULTG, var_name="QMULTG", &
             var_description="Splintering from droplets accreted onto graupel, +ricem, -rcm &
             &[(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PSACWG')
        iPSACWG = k

        call stat_assign( var_index=iPSACWG, var_name="PSACWG", &
             var_description="Collection of cloud water by graupel, +rgraupelm, -rcm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PGSACW')
        iPGSACW = k

        call stat_assign( var_index=iPGSACW, var_name="PGSACW", &
             var_description="Reclassification of rimed snow as graupel, +rgraupelm, -rcm &
             &[(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PRD')
        iPRD = k

        call stat_assign( var_index=iPRD, var_name="PRD", &
             var_description="Depositional growth of cloud ice, +ricem, -rvm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PRCI')
        iPRCI = k

        call stat_assign( var_index=iPRCI, var_name="PRCI", &
             var_description="Autoconversion of cloud ice to snow, +rsnowm, -ricem [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PRAI')
        iPRAI = k

        call stat_assign( var_index=iPRAI, var_name="PRAI", &
             var_description="Collection of cloud ice by snow, +rsnowm, -ricem [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('QMULTR')
        iQMULTR = k

        call stat_assign( var_index=iQMULTR, var_name="QMULTR", &
             var_description="Splintering from rain droplets accreted onto snow, +ricem, -rrainm &
             &[(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('QMULTRG')
        iQMULTRG = k

        call stat_assign( var_index=iQMULTRG, var_name="QMULTRG", &
             var_description="Splintering from rain droplets accreted onto graupel, +ricem, -rrainm&
             & [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('MNUCCD')
        iMNUCCD = k

        call stat_assign( var_index=iMNUCCD, var_name="MNUCCD", &
             var_description="Freezing of aerosol, +ricem, -rvm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PRACI')
        iPRACI = k

        call stat_assign( var_index=iPRACI, var_name="PRACI", &
             var_description="Collection of cloud ice by rain, +rgraupelm, -ricem [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PRACIS')
        iPRACIS = k

        call stat_assign( var_index=iPRACIS, var_name="PRACIS", &
             var_description="Collection of cloud ice by rain, +rsnowm, -ricem [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('EPRD')
        iEPRD = k

        call stat_assign( var_index=iEPRD, var_name="EPRD", &
             var_description="Negative of sublimation of cloud ice, +ricem, -rvm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('MNUCCR')
        iMNUCCR = k

        call stat_assign( var_index=iMNUCCR, var_name="MNUCCR", &
             var_description="Contact freezing of rain droplets, +rgraupelm, -rrainm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PIACR')
        iPIACR = k

        call stat_assign( var_index=iPIACR, var_name="PIACR", &
             var_description="Collection of cloud ice by rain, +rgraupelm, -rrainm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PIACRS')
        iPIACRS = k

        call stat_assign( var_index=iPIACRS, var_name="PIACRS", &
             var_description="Collection of cloud ice by rain, +rsnowm, -rrainm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PGRACS')
        iPGRACS = k

        call stat_assign( var_index=iPGRACS, var_name="PGRACS", &
             var_description="Collection of rain by snow, +rgraupelm, -rrainm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PRDS')
        iPRDS = k

        call stat_assign( var_index=iPRDS, var_name="PRDS", &
             var_description="Depositional growth of snow, +rsnowm, -rvm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('EPRDS')
        iEPRDS = k

        call stat_assign( var_index=iEPRDS, var_name="EPRDS", &
             var_description="Negative of sublimation of snow, +rsnowm, -rvm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PSACR')
        iPSACR = k

        call stat_assign( var_index=iPSACR, var_name="PSACR", &
             var_description="Collection of snow by rain, +rgraupelm, -rsnowm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PRDG')
        iPRDG = k

        call stat_assign( var_index=iPRDG, var_name="PRDG", &
             var_description="Depositional growth of graupel, +rgraupelm, -rvm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('EPRDG')
        iEPRDG = k

        call stat_assign( var_index=iEPRDG, var_name="EPRDG", &
             var_description="Negative of sublimation of graupel, +rgraupelm, -rvm [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NGSTEN')
        iNGSTEN = k

        call stat_assign( var_index=iNGSTEN, var_name="NGSTEN", &
             var_description="Graupel sedimentation tendency [(#/kg/s)]", var_units="(#/kg/s)", &
             l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NRSTEN')
        iNRSTEN = k

        call stat_assign( var_index=iNRSTEN, var_name="NRSTEN", &
             var_description="Rain sedimentation tendency [(#/kg/s)]", var_units="(#/kg/s)", &
             l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NISTEN')
        iNISTEN = k

        call stat_assign( var_index=iNISTEN, var_name="NISTEN", &
             var_description="Cloud ice sedimentation tendency [(#/kg/s)]", var_units="(#/kg/s)", &
             l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NSSTEN')
        iNSSTEN = k

        call stat_assign( var_index=iNSSTEN, var_name="NSSTEN", &
             var_description="Snow sedimentation tendency [(#/kg/s)]", var_units="(#/kg/s)", &
             l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NCSTEN')
        iNCSTEN = k

        call stat_assign( var_index=iNCSTEN, var_name="NCSTEN", &
             var_description="Cloud water sedimentation tendency [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1
      
      case ('NPRC1')
        iNPRC1 = k

        call stat_assign( var_index=iNPRC1, var_name="NPRC1", &
             var_description="Change in Nrm due to autoconversion of droplets, +Nrm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1
      
      case ('NRAGG')
        iNRAGG = k

        call stat_assign( var_index=iNRAGG, var_name="NRAGG", &
             var_description="Change in Nrm due to self-collection of raindrops, +Nrm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1
      
      case ('NPRACG')
        iNPRACG = k

        call stat_assign( var_index=iNPRACG, var_name="NPRACG", &
             var_description="Collection of rainwater by graupel, -Nrm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1
      
      case ('NSUBR')
        iNSUBR = k

        call stat_assign( var_index=iNSUBR, var_name="NSUBR", &
             var_description="Loss of Nrm by evaporation, +Nrm [(#/kg/s)]", var_units="(#/kg/s)", &
             l_silhs=.true., grid_kind=zt )
        k = k + 1
      
      case ('NSMLTR')
        iNSMLTR = k

        call stat_assign( var_index=iNSMLTR, var_name="NSMLTR", &
             var_description="Melting of snow to form rain, -Nrm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1
      
      case ('NGMLTR')
        iNGMLTR = k

        call stat_assign( var_index=iNGMLTR, var_name="NGMLTR", &
             var_description="Melting of graupel to form rain, -Nrm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1
      
      case ('NPRACS')
        iNPRACS = k

        call stat_assign( var_index=iNPRACS, var_name="NPRACS", &
             var_description="Collection of rainwater by snow, -Nrm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1
      
      case ('NNUCCR')
        iNNUCCR = k

        call stat_assign( var_index=iNNUCCR, var_name="NNUCCR", &
             var_description="Contact freezing of rain, +Ngraupelm, -Nrm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1
      
      case ('NIACR')
        iNIACR = k

        call stat_assign( var_index=iNIACR, var_name="NIACR", &
             var_description="Collection of cloud ice by rain, +Ngraupelm, -Nrm, -Nim [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1
      
      case ('NIACRS')
        iNIACRS = k

        call stat_assign( var_index=iNIACRS, var_name="NIACRS", &
             var_description="Collection of cloud ice by rain, +Nsnowm, -Nrm, -Nim [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1
      
      case ('NGRACS')
        iNGRACS = k

        call stat_assign( var_index=iNGRACS, var_name="NGRACS", &
             var_description="Collection of rain by snow, +Ngraupelm, -Nrm, -Nsnowm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NSMLTS')
        iNSMLTS= k

        call stat_assign( var_index=iNSMLTS, var_name="NSMLTS", &
             var_description="Melting  of snow, +Nsnowm [(#/kg/s)]", var_units="(#/kg/s)", &
             l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NSAGG')
        iNSAGG= k

        call stat_assign( var_index=iNSAGG, var_name="NSAGG", &
             var_description="Self collection of snow, +Nsnowm [(#/kg/s)]", var_units="(#/kg/s)", &
             l_silhs=.true., grid_kind=zt )

        k = k + 1

      case ('NPRCI')
        iNPRCI= k

        call stat_assign( var_index=iNPRCI, var_name="NPRCI", &
             var_description="Autoconversion of cloud ice to snow, -Nim, +Nsnowm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NSCNG')
        iNSCNG= k

        call stat_assign( var_index=iNSCNG, var_name="NSCNG", &
             var_description="Conversion of snow to graupel, +Ngraupelm, -Nsnowm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NSUBS')
        iNSUBS= k

        call stat_assign( var_index=iNSUBS, var_name="NSUBS", &
             var_description="Loss of snow due to sublimation, +Nsnowm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PRC')
        iPRC= k
             
        call stat_assign( var_index=iPRC, var_name="PRC", &
             var_description="Autoconversion +rrainm -rcm [(kg/kg/s)]", var_units="(kg/kg/s)", &
             l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PRA')
        iPRA= k

        call stat_assign( var_index=iPRA, var_name="PRA", &
             var_description="Accretion +rrainm -rcm [(kg/kg/s)]", var_units="(kg/kg/s)", &
             l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PRE')
        iPRE= k

        call stat_assign( var_index=iPRE, var_name="PRE", &
             var_description="Evaporation of rain -rrainm [(kg/kg/s)]", var_units="(kg/kg/s)", &
             l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('PCC')
        iPCC= k

        call stat_assign( var_index=iPCC, var_name="PCC", &
             var_description="Satuation adjustment -rvm +rcm [(kg/kg/s)]", var_units="(kg/kg/s)", &
             l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NNUCCC')
        iNNUCCC= k

        call stat_assign( var_index=iNNUCCC, var_name="NNUCCC", &
             var_description="Contact freezing of drops, -Ncm + Nim [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NPSACWS')
        iNPSACWS= k

        call stat_assign( var_index=iNPSACWS, var_name="NPSACWS", &
             var_description="Droplet accretion by snow, -Ncm [(#/kg/s)]", var_units="(#/kg/s)", &
             l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NPRA')
        iNPRA= k

        call stat_assign( var_index=iNPRA, var_name="NPRA", &
             var_description="Droplet accretion by rain, -Ncm [(#/kg/s)]", var_units="(#/kg/s)", &
             l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NPRC')
        iNPRC= k

        call stat_assign( var_index=iNPRC, var_name="NPRC", &
             var_description="Autoconversion of cloud drops, -Ncm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NPSACWI')
        iNPSACWI= k

        call stat_assign( var_index=iNPSACWI, var_name="NPSACWI", &
             var_description="Droplet accretion by cloud ice, -Ncm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NPSACWG')
        iNPSACWG= k

        call stat_assign( var_index=iNPSACWG, var_name="NPSACWG", &
             var_description="Collection of cloud droplets by graupel, -Ncm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NPRAI')
        iNPRAI= k

        call stat_assign( var_index=iNPRAI, var_name="NPRAI", &
             var_description="Accretion of cloud ice by snow, -Nim [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NMULTS')
        iNMULTS= k

        call stat_assign( var_index=iNMULTS, var_name="NMULTS", &
             var_description="Ice multiplication due to riming of cloud droplets by snow, +Nim &
             &[(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NMULTG')
        iNMULTG= k

        call stat_assign( var_index=iNMULTG, var_name="NMULTG", &
             var_description="Ice multiplication due to accretion of droplets by graupel, +Nim &
             &[(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NMULTR')
        iNMULTR= k

        call stat_assign( var_index=iNMULTR, var_name="NMULTR", &
             var_description="Ice multiplication due to riming of rain by snow, +Nim [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NMULTRG')
        iNMULTRG= k

        call stat_assign( var_index=iNMULTRG, var_name="NMULTRG", &
             var_description="Ice multiplication due to accretion of rain by graupel, +Nim &
             &[(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NNUCCD')
        iNNUCCD= k

        call stat_assign( var_index=iNNUCCD, var_name="NNUCCD", &
             var_description="Primary ice nucleation, freezing of aerosol, +Nim [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NSUBI')
        iNSUBI= k

        call stat_assign( var_index=iNSUBI, var_name="NSUBI", &
             var_description="Loss of ice due to sublimation, -Nim [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NGMLTG')
        iNGMLTG= k

        call stat_assign( var_index=iNGMLTG, var_name="NGMLTG", &
             var_description="Loss of graupel due to melting, -Ngraupelm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NSUBG')
        iNSUBG= k

        call stat_assign( var_index=iNSUBG, var_name="NSUBG", &
             var_description="Loss of graupel due to sublimation, -Ngraupelm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NACT')
        iNACT= k

        call stat_assign( var_index=iNACT, var_name="NACT", &
             var_description="Cloud drop formation by aerosol activation, +Ncm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('SIZEFIX_NR')
        iSIZEFIX_NR= k

        call stat_assign( var_index=iSIZEFIX_NR, var_name="SIZEFIX_NR", &
             var_description="Adjust rain # conc. for large/small drops, +Nrm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('SIZEFIX_NC')
        iSIZEFIX_NC= k

        call stat_assign( var_index=iSIZEFIX_NC, var_name="SIZEFIX_NC", &
             var_description="Adjust cloud # conc. for large/small drops, +Ncm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('SIZEFIX_NI')
        iSIZEFIX_NI= k

        call stat_assign( var_index=iSIZEFIX_NI, var_name="SIZEFIX_NI", &
             var_description="Adjust ice # conc. for large/small drops, +Nim [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('SIZEFIX_NS')
        iSIZEFIX_NS= k

        call stat_assign( var_index=iSIZEFIX_NS, var_name="SIZEFIX_NS", &
             var_description="Adjust snow # conc. for large/small drops, +Nsnowm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('SIZEFIX_NG')
        iSIZEFIX_NG= k

        call stat_assign( var_index=iSIZEFIX_NG, var_name="SIZEFIX_NG", &
             var_description="Adjust graupel # conc. for large/small drops,+Ngraupelm [(#/kg/s)]",&
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NEGFIX_NR')
        iNEGFIX_NR= k

        call stat_assign( var_index=iNEGFIX_NR, var_name="NEGFIX_NR", &
             var_description="Removal of negative rain drop number conc., -Nrm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NEGFIX_NC')
        iNEGFIX_NC= k

        call stat_assign( var_index=iNEGFIX_NC, var_name="NEGFIX_NC", &
             var_description="Removal of negative cloud drop number conc., -Ncm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NEGFIX_NI')
        iNEGFIX_NI= k

        call stat_assign( var_index=iNEGFIX_NI, var_name="NEGFIX_NI", &
             var_description="Removal of negative ice number conc., -Nim [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NEGFIX_NS')
        iNEGFIX_NS= k

        call stat_assign( var_index=iNEGFIX_NS, var_name="NEGFIX_NS", &
             var_description="Removal of negative snow number conc,, -Nsnowm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NEGFIX_NG')
        iNEGFIX_NG= k

        call stat_assign( var_index=iNEGFIX_NG, var_name="NEGFIX_NG", &
             var_description="Removal of negative graupel number conc., -Ngraupelm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NIM_MORR_CL')
        iNIM_MORR_CL= k

        call stat_assign( var_index=iNIM_MORR_CL, var_name="NIM_MORR_CL", &
             var_description="Clipping of large ice concentrations, -Nim [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('QC_INST')
        iQC_INST= k

        call stat_assign( var_index=iQC_INST, var_name="QC_INST", &
             var_description="Change in mixing ratio due to instantaneous processes," // &
                             " +rcm [(kg/kg/s)]", &
             var_units="(kg/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('QR_INST')
        iQR_INST= k

        call stat_assign( var_index=iQR_INST, var_name="QR_INST", &
             var_description="Change in mixing ratio from instantaneous processes," // &
                             " +rrainm [(kg/kg/s)]", &
             var_units="(kg/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('QI_INST')
        iQI_INST= k

        call stat_assign( var_index=iQI_INST, var_name="QI_INST", &
             var_description="Change in mixing ratio from instantaneous processes," // &
                             " +ricem [(kg/kg/s)]", &
             var_units="(kg/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('QS_INST')
        iQS_INST= k

        call stat_assign( var_index=iQS_INST, var_name="QS_INST", &
             var_description="Change in mixing ratio from instantaneous processes," // &
                             " +rsnowm [(kg/kg/s)]", &
             var_units="(kg/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('QG_INST')
        iQG_INST= k

        call stat_assign( var_index=iQG_INST, var_name="QG_INST", &
             var_description="Change in mixing ratio from instantaneous processes," // &
                             " +rgraupelm [(kg/kg/s)]", &
             var_units="(kg/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NC_INST')
        iNC_INST= k

        call stat_assign( var_index=iNC_INST, var_name="NC_INST", &
             var_description="Change in # conc. from instantaneous processes," // &
                             " +Ncm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NR_INST')
        iNR_INST= k

        call stat_assign( var_index=iNR_INST, var_name="NR_INST", &
             var_description="Change in # conc. from instantaneous processes," // &
                             " +Nrm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NI_INST')
        iNI_INST= k

        call stat_assign( var_index=iNI_INST, var_name="NI_INST", &
             var_description="Change in # conc. from instantaneous processes," // &
                             " +Nim [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NS_INST')
        iNS_INST= k

        call stat_assign( var_index=iNS_INST, var_name="NS_INST", &
             var_description="Change in # conc. from instantaneous processes," // &
                             " +Nsnowm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('NG_INST')
        iNG_INST= k

        call stat_assign( var_index=iNG_INST, var_name="NG_INST", &
             var_description="Change in # conc. from instantaneous processes," // &
                             " +Ngraupelm [(#/kg/s)]", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1


      case ('T_in_K_mc')
        iT_in_K_mc= k

        call stat_assign( var_index=iT_in_K_mc, var_name="T_in_K_mc", &
             var_description="Temperature tendency from Morrison microphysics [(K/s)]", &
             var_units="(K/s)", l_silhs=.true., grid_kind=zt )
        k = k + 1

      case ('w_KK_evap_covar_zt')
       iw_KK_evap_covar_zt = k

        call stat_assign( var_index=iw_KK_evap_covar_zt, var_name="w_KK_evap_covar_zt", &
             var_description="Covariance of w and KK evaporation rate", &
             var_units="m*(kg/kg)/s^2", l_silhs=.false., grid_kind=zt )
       k = k + 1

      case ('rt_KK_evap_covar_zt')
       irt_KK_evap_covar_zt = k

        call stat_assign( var_index=irt_KK_evap_covar_zt, var_name="rt_KK_evap_covar_zt", &
             var_description="Covariance of r_t and KK evaporation rate", &
             var_units="(kg/kg)^2/s", l_silhs=.false., grid_kind=zt )
       k = k + 1

      case ('thl_KK_evap_covar_zt')
       ithl_KK_evap_covar_zt = k

        call stat_assign( var_index=ithl_KK_evap_covar_zt, var_name="thl_KK_evap_covar_zt", &
             var_description="Covariance of theta_l and KK evaporation rate", &
             var_units="K*(kg/kg)/s", l_silhs=.false., grid_kind=zt )
       k = k + 1

      case ('w_KK_auto_covar_zt')
       iw_KK_auto_covar_zt = k

        call stat_assign( var_index=iw_KK_auto_covar_zt, var_name="w_KK_auto_covar_zt", &
             var_description="Covariance of w and KK autoconversion rate", &
             var_units="m*(kg/kg)/s^2", l_silhs=.false., grid_kind=zt )
       k = k + 1

      case ('rt_KK_auto_covar_zt')
       irt_KK_auto_covar_zt = k

        call stat_assign( var_index=irt_KK_auto_covar_zt, var_name="rt_KK_auto_covar_zt", &
             var_description="Covariance of r_t and KK autoconversion rate", &
             var_units="(kg/kg)^2/s", l_silhs=.false., grid_kind=zt )
       k = k + 1

      case ('thl_KK_auto_covar_zt')
       ithl_KK_auto_covar_zt = k

        call stat_assign( var_index=ithl_KK_auto_covar_zt, var_name="thl_KK_auto_covar_zt", &
             var_description="Covariance of theta_l and KK autoconversion rate", &
             var_units="K*(kg/kg)/s", l_silhs=.false., grid_kind=zt )
       k = k + 1

      case ('w_KK_accr_covar_zt')
       iw_KK_accr_covar_zt = k

        call stat_assign( var_index=iw_KK_accr_covar_zt, var_name="w_KK_accr_covar_zt", &
             var_description="Covariance of w and KK accretion rate", var_units="m*(kg/kg)/s^2", &
             l_silhs=.false., grid_kind=zt )
       k = k + 1

      case ('rt_KK_accr_covar_zt')
       irt_KK_accr_covar_zt = k

        call stat_assign( var_index=irt_KK_accr_covar_zt, var_name="rt_KK_accr_covar_zt", &
             var_description="Covariance of r_t and KK accretion rate", var_units="(kg/kg)^2/s", &
             l_silhs=.false., grid_kind=zt )
       k = k + 1

      case ('thl_KK_accr_covar_zt')
       ithl_KK_accr_covar_zt = k

        call stat_assign( var_index=ithl_KK_accr_covar_zt, var_name="thl_KK_accr_covar_zt", &
             var_description="Covariance of theta_l and KK accretion rate", &
             var_units="K*(kg/kg)/s", l_silhs=.false., grid_kind=zt )
       k = k + 1

      case ('rr_KK_mvr_covar_zt')
       irr_KK_mvr_covar_zt = k

        call stat_assign( var_index=irr_KK_mvr_covar_zt, var_name="rr_KK_mvr_covar_zt", &
             var_description="Covariance of r_r and KK rain drop mean volume radius [(kg/kg)m]", &
             var_units="(kg/kg)m", l_silhs=.false., grid_kind=zt )
       k = k + 1

      case ('Nr_KK_mvr_covar_zt')
       iNr_KK_mvr_covar_zt = k

        call stat_assign( var_index=iNr_KK_mvr_covar_zt, var_name="Nr_KK_mvr_covar_zt", &
             var_description="Covariance of N_r and KK rain drop mean volume radius [(num/kg)m]", &
             var_units="(num/kg)m", l_silhs=.false., grid_kind=zt )
       k = k + 1

      case ('KK_mvr_variance_zt')
       iKK_mvr_variance_zt = k

        call stat_assign( var_index=iKK_mvr_variance_zt, var_name="KK_mvr_variance_zt", &
             var_description="Variance of KK rain drop mean volume radius [m^2]", &
             var_units="m^2", l_silhs=.false., grid_kind=zt )
       k = k + 1

      case ('vm_bt')
        ivm_bt = k

        call stat_assign( var_index=ivm_bt, var_name="vm_bt", &
             var_description="vm budget: vm time tendency [m s^{-2}]", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('vm_ma')
        ivm_ma = k
        call stat_assign( var_index=ivm_ma, var_name="vm_ma", &
             var_description="vm budget: vm vertical mean advection [m s^{-2}]", &
             var_units="m s^{-2}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('vm_gf')
        ivm_gf = k

        call stat_assign( var_index=ivm_gf, var_name="vm_gf", &
             var_description="vm budget: vm geostrophic forcing [m s^{-2}]", &
             var_units="m s^{-2}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('vm_cf')
        ivm_cf = k

        call stat_assign( var_index=ivm_cf, var_name="vm_cf", &
             var_description="vm budget: vm coriolis forcing [m s^{-2}]", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('vm_ta')
        ivm_ta = k

        call stat_assign( var_index=ivm_ta, var_name="vm_ta", &
             var_description="vm budget: vm turbulent transport [m s^{-2}]", &
             var_units="m s^{-2}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('vm_f')
        ivm_f = k
        call stat_assign( var_index=ivm_f, var_name="vm_f", &
             var_description="vm budget: vm forcing [m s^{-2}]", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('vm_sdmp')
        ivm_sdmp = k
        call stat_assign( var_index=ivm_sdmp, var_name="vm_sdmp", &
             var_description="vm budget: vm sponge damping [m s^{-2}]", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('vm_ndg')
        ivm_ndg = k
        call stat_assign( var_index=ivm_ndg, var_name="vm_ndg", &
             var_description="vm budget: vm nudging [m s^{-2}]", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('um_bt')
        ium_bt = k

        call stat_assign( var_index=ium_bt, var_name="um_bt", &
             var_description="um budget: um time tendency [m s^{-2}]", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('um_ma')
        ium_ma = k

        call stat_assign( var_index=ium_ma, var_name="um_ma", &
             var_description="um budget: um vertical mean advection [m s^{-2}]", &
             var_units="m s^{-2}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('um_gf')
        ium_gf = k
        call stat_assign( var_index=ium_gf, var_name="um_gf", &
             var_description="um budget: um geostrophic forcing [m s^{-2}]", &
             var_units="m s^{-2}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('um_cf')
        ium_cf = k
        call stat_assign( var_index=ium_cf, var_name="um_cf", &
             var_description="um budget: um coriolis forcing [m s^{-2}]", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('um_ta')
        ium_ta = k
        call stat_assign( var_index=ium_ta, var_name="um_ta", &
             var_description="um budget: um turbulent advection [m s^{-2}]", &
             var_units="m s^{-2}", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('um_f')
        ium_f = k
        call stat_assign( var_index=ium_f, var_name="um_f", &
             var_description="um budget: um forcing [m s^{-2}]", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('um_sdmp')
        ium_sdmp = k
        call stat_assign( var_index=ium_sdmp, var_name="um_sdmp", &
             var_description="um budget: um sponge damping [m s^{-2}]", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('um_ndg')
        ium_ndg = k
        call stat_assign( var_index=ium_ndg, var_name="um_ndg", &
             var_description="um budget: um nudging [m s^{-2}]", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('mixt_frac')
        imixt_frac = k
        call stat_assign( var_index=imixt_frac, var_name="mixt_frac", &
             var_description="pdf parameter: mixture fraction [count]", var_units="count", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('w1')
        iw1 = k
        call stat_assign( var_index=iw1, var_name="w1", &
             var_description="pdf parameter: mean w of component 1 [m/s]", var_units="m/s", &
             l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('w2')
        iw2 = k

        call stat_assign( var_index=iw2, var_name="w2", &
             var_description="pdf paramete: mean w of component 2 [m/s]", var_units="m/s", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('varnce_w1')
        ivarnce_w1 = k
        call stat_assign( var_index=ivarnce_w1, var_name="varnce_w1", &
             var_description="pdf parameter: w variance of component 1 [m^2/s^2]", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('varnce_w2')
        ivarnce_w2 = k

        call stat_assign( var_index=ivarnce_w2, var_name="varnce_w2", &
             var_description="pdf parameter: w variance of component 2 [m^2/s^2]", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('thl1')
        ithl1 = k

        call stat_assign( var_index=ithl1, var_name="thl1", &
             var_description="pdf parameter: mean thl of component 1 [K]", var_units="K", &
             l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('thl2')
        ithl2 = k

        call stat_assign( var_index=ithl2, var_name="thl2", &
             var_description="pdf parameter: mean thl of component 2 [K]", var_units="K", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('varnce_thl1')
        ivarnce_thl1 = k

        call stat_assign( var_index=ivarnce_thl1, var_name="varnce_thl1", &
             var_description="pdf parameter: thl variance of component 1 [K^2]", var_units="K^2", &
             l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('varnce_thl2')
        ivarnce_thl2 = k
        call stat_assign( var_index=ivarnce_thl2, var_name="varnce_thl2", &
             var_description="pdf parameter: thl variance of component 2 [K^2]", var_units="K^2", &
             l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('rt1')
        irt1 = k
        call stat_assign( var_index=irt1, var_name="rt1", &
             var_description="pdf parameter: mean rt of component 1 [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )

        k = k + 1

      case ('rt2')
        irt2 = k

        call stat_assign( var_index=irt2, var_name="rt2", &
             var_description="pdf parameter: mean rt of component 2 [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('varnce_rt1')
        ivarnce_rt1 = k
        call stat_assign( var_index=ivarnce_rt1, var_name="varnce_rt1", &
             var_description="pdf parameter: rt variance of component 1 [(kg^2)/(kg^2)]", &
             var_units="(kg^2)/(kg^2)", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('varnce_rt2')
        ivarnce_rt2 = k

        call stat_assign( var_index=ivarnce_rt2, var_name="varnce_rt2", &
             var_description="pdf parameter: rt variance of component 2 [(kg^2)/(kg^2)]", &
             var_units="(kg^2)/(kg^2)", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rc1')
        irc1 = k

        call stat_assign( var_index=irc1, var_name="rc1", &
             var_description="pdf parameter: mean rc of component 1 [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rc2')
        irc2 = k

        call stat_assign( var_index=irc2, var_name="rc2", &
             var_description="pdf parameter: mean rc of component 2 [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rsl1')
        irsl1 = k

        call stat_assign( var_index=irsl1, var_name="rsl1", &
             var_description="pdf parameter: sat mix rat based on tl1 [kg/kg]", &
             var_units="kg/kg", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rsl2')
        irsl2 = k

        call stat_assign( var_index=irsl2, var_name="rsl2", &
             var_description="pdf parameter: sat mix rat based on tl2 [kg/kg]", &
             var_units="kg/kg", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('cloud_frac1')
        icloud_frac1 = k
        call stat_assign( var_index=icloud_frac1, var_name="cloud_frac1", &
             var_description="pdf parameter cloud_frac1 [count]", var_units="count", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('cloud_frac2')
        icloud_frac2 = k

        call stat_assign( var_index=icloud_frac2, var_name="cloud_frac2", &
             var_description="pdf parameter cloud_frac2 [count]", var_units="count", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('s1')
        is1 = k

        call stat_assign( var_index=is1, var_name="s1", &
             var_description="pdf parameter: Mellor's s (extended liq) for component 1 [kg/kg]", &
             var_units="kg/kg", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('s2')
        is2 = k

        call stat_assign( var_index=is2, var_name="s2", &
             var_description="pdf parameter: Mellor's s (extended liq) for component 2 [kg/kg]", &
             var_units="kg/kg", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('stdev_s1')
        istdev_s1 = k

        call stat_assign( var_index=istdev_s1, var_name="stdev_s1", &
             var_description="pdf parameter: Std dev of s1 [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('stdev_s2')
        istdev_s2 = k

        call stat_assign( var_index=istdev_s2, var_name="stdev_s2", &
             var_description="pdf parameter: Std dev of s2 [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('sp2')
        isp2 = k
        call stat_assign( var_index=isp2, var_name="sp2", &
             var_description="Variance of s (overall) [(kg/kg)^2]", var_units="(kg/kg)^2", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('stdev_t1')
        istdev_t1 = k

        call stat_assign( var_index=istdev_t1, var_name="stdev_t1", &
             var_description="Standard dev. of t (1st PDF component) [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('stdev_t2')
        istdev_t2 = k

        call stat_assign( var_index=istdev_t2, var_name="stdev_t2", &
             var_description="Standard dev. of t (2nd PDF component) [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('covar_st_1')
        icovar_st_1 = k

        call stat_assign( var_index=icovar_st_1, var_name="covar_st_1", &
             var_description="Covariance of s and t (1st PDF component) [kg^2/kg^2]", &
             var_units="kg^2/kg^2", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('covar_st_2')
        icovar_st_2 = k

        call stat_assign( var_index=icovar_st_2, var_name="covar_st_2", &
             var_description="Covariance of s and t (2nd PDF component) [kg^2/kg^2]", &
             var_units="kg^2/kg^2", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('corr_st_1')
        icorr_st_1 = k

        call stat_assign( var_index=icorr_st_1, var_name="corr_st_1", &
             var_description="Correlation btw. s and t (1st PDF component) [-]", var_units="-", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('corr_st_2')
        icorr_st_2 = k

        call stat_assign( var_index=icorr_st_2, var_name="corr_st_2", &
             var_description="Correlation btw. s and t (2nd PDF component) [-]", var_units="-", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rrtthl')
        irrtthl = k

        call stat_assign( var_index=irrtthl, var_name="rrtthl", &
             var_description="Correlation btw. rt and thl (both components) [-]", var_units="-", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('crt1')
        icrt1 = k

        call stat_assign( var_index=icrt1, var_name="crt1", &
             var_description=" Coef. on r_t in s/t eqns. (1st PDF comp.)  [-]", &
             var_units="count", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('crt2')
        icrt2 = k

        call stat_assign( var_index=icrt2, var_name="crt2", &
             var_description=" Coef. on r_t in s/t eqns. (2nd PDF comp.)  [-]", &
             var_units="count", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('cthl1')
        icthl1 = k

        call stat_assign( var_index=icthl1, var_name="cthl1", &
             var_description=" Coef. on theta_l in s/t eqns. (1st PDF comp.)  [kg/kg/K]", &
             var_units="kg/kg/K", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('cthl2')
        icthl2 = k

        call stat_assign( var_index=icthl2, var_name="cthl2", &
             var_description=" Coef. on theta_l in s/t eqns. (2nd PDF comp.)  [kg/kg/K]", &
             var_units="kg/kg/K", l_silhs=.false., grid_kind=zt )
        k = k + 1


      case('wp2_zt')
        iwp2_zt = k

        call stat_assign( var_index=iwp2_zt, var_name="wp2_zt", &
             var_description="w'^2 interpolated to thermodyamic levels [m^2/s^2]", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case('thlp2_zt')
        ithlp2_zt = k

        call stat_assign( var_index=ithlp2_zt, var_name="thlp2_zt", &
             var_description="thl'^2 interpolated to thermodynamic levels [K^2]", &
             var_units="K^2", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case('wpthlp_zt')
        iwpthlp_zt = k

        call stat_assign( var_index=iwpthlp_zt, var_name="wpthlp_zt", &
             var_description="w'thl' interpolated to thermodynamic levels [(m K)/s]", &
             var_units="(m K)/s", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case('wprtp_zt')
        iwprtp_zt = k

        call stat_assign( var_index=iwprtp_zt, var_name="wprtp_zt", &
             var_description="w'rt' interpolated to thermodynamic levels [(m kg)/(s kg)]", &
             var_units="(m kg)/(s kg)", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case('rtp2_zt')
        irtp2_zt = k

        call stat_assign( var_index=irtp2_zt, var_name="rtp2_zt", &
             var_description="rt'^2 interpolated to thermodynamic levels [(kg/kg)^2]", &
             var_units="(kg/kg)^2", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case('rtpthlp_zt')
        irtpthlp_zt = k

        call stat_assign( var_index=irtpthlp_zt, var_name="rtpthlp_zt", &
             var_description="rt'thl' interpolated to thermodynamic levels [(kg K)/kg]", &
             var_units="(kg K)/kg", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('up2_zt')
        iup2_zt = k
        call stat_assign( var_index=iup2_zt, var_name="up2_zt", &
             var_description="u'^2 interpolated to thermodynamic levels [m^2/s^2]", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('vp2_zt')
        ivp2_zt = k
        call stat_assign( var_index=ivp2_zt, var_name="vp2_zt", &
             var_description="v'^2 interpolated to thermodynamic levels [m^2/s^2]", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('upwp_zt')
        iupwp_zt = k
        call stat_assign( var_index=iupwp_zt, var_name="upwp_zt", &
             var_description="u'w' interpolated to thermodynamic levels [m^2/s^2]", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('vpwp_zt')
        ivpwp_zt = k
        call stat_assign( var_index=ivpwp_zt, var_name="vpwp_zt", &
             var_description="v'w' interpolated to thermodynamic levels [m^2/s^2]", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=zt )
        k = k + 1

      ! Hydrometeor overall variances for each hydrometeor type.
      case('hmp2_zt')

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The overall variance of the hydrometeor.
            ihmp2_zt(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=ihmp2_zt(hm_idx), &
                                 var_name=trim( hm_type(1:2) )//"p2_zt", &
                                 var_description="<" &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // "'^2> on thermodyamic levels [(kg/kg)^2]", &
                                 var_units="(kg/kg)^2", &
                                 l_silhs=.false., grid_kind=zt )

            else ! Concentration

               call stat_assign( var_index=ihmp2_zt(hm_idx), &
                                 var_name=trim( hm_type(1:2) )//"p2_zt", &
                                 var_description="<" &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // "'^2> on thermodyamic levels " &
                                 // "[(num/kg)^2]", &
                                 var_units="(num/kg)^2", &
                                 l_silhs=.false., grid_kind=zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ('C11_Skw_fnc')
        iC11_Skw_fnc = k

        call stat_assign( var_index=iC11_Skw_fnc, var_name="C11_Skw_fnc", &
             var_description="C_11 parameter with Sk_w applied [-]", var_units="count", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('s_mellor')
        is_mellor = k

        call stat_assign( var_index=is_mellor, var_name="s_mellor", &
             var_description="Mellor's s (extended liq) [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'a3_coef_zt' )
        ia3_coef_zt = k
        call stat_assign( var_index=ia3_coef_zt, var_name="a3_coef_zt", &
             var_description="The a3 coefficient interpolated the the zt grid [-]", &
             var_units="count", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'wp3_on_wp2_zt' )
        iwp3_on_wp2_zt = k
        call stat_assign( var_index=iwp3_on_wp2_zt, var_name="wp3_on_wp2_zt", &
             var_description="Smoothed version of wp3 / wp2 [m/s]", var_units="m/s", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      ! Hydrometeor component mean values for each PDF component and hydrometeor
      ! type.
      case ( "hmi" )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The mean of the hydrometeor in the 1st PDF component.
            ihm1(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=ihm1(hm_idx), &
                                 var_name=trim( hm_type(1:2) )//"1", &
                                 var_description="Mean of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (1st PDF component) [kg/kg]", &
                                 var_units="kg/kg", &
                                 l_silhs=.false., grid_kind=zt )

            else ! Concentration

               call stat_assign( var_index=ihm1(hm_idx), &
                                 var_name=trim( hm_type(1:2) )//"1", &
                                 var_description="Mean of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (1st PDF component) [num/kg]", &
                                 var_units="num/kg", &
                                 l_silhs=.false., grid_kind=zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

            ! The mean of the hydrometeor in the 2nd PDF component.
            ihm2(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=ihm2(hm_idx), &
                                 var_name=trim( hm_type(1:2) )//"2", &
                                 var_description="Mean of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (2nd PDF component) [kg/kg]", &
                                 var_units="kg/kg", &
                                 l_silhs=.false., grid_kind=zt )

            else ! Concentration

               call stat_assign( var_index=ihm2(hm_idx), &
                                 var_name=trim( hm_type(1:2) )//"2", &
                                 var_description="Mean of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (2nd PDF component) [num/kg]", &
                                 var_units="num/kg", &
                                 l_silhs=.false., grid_kind=zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'LWP1' )
        iLWP1 = k
        call stat_assign( var_index=iLWP1, var_name="LWP1", &
             var_description="Liquid water path (1st PDF component) [kg/m^2]", &
             var_units="kg/m^2", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'LWP2' )
        iLWP2 = k
        call stat_assign( var_index=iLWP2, var_name="LWP2", &
             var_description="Liquid water path (2nd PDF component) [kg/m^2]", &
             var_units="kg/m^2", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'precip_frac' )
        iprecip_frac = k
        call stat_assign( var_index=iprecip_frac, var_name="precip_frac", &
             var_description="Precipitation Fraction [-]", var_units="-", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ( 'precip_frac_1' )
        iprecip_frac_1 = k
        call stat_assign( var_index=iprecip_frac_1, var_name="precip_frac_1", &
             var_description="Precipitation Fraction (1st PDF component) [-]", var_units="-", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'precip_frac_2' )
        iprecip_frac_2 = k
        call stat_assign( var_index=iprecip_frac_2, var_name="precip_frac_2", &
             var_description="Precipitation Fraction (2nd PDF component) [-]", var_units="-", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'Ncnm' )
        iNcnm = k
        call stat_assign( var_index=iNcnm, var_name="Ncnm", &
             var_description="Cloud nuclei concentration (PDF) [num/kg]", &
             var_units="num/kg", l_silhs=.false., grid_kind=zt )
        k = k + 1

      ! Hydrometeor component mean values (in-precip) for each PDF component and
      ! hydrometeor type.
      case ( 'mu_hm_i' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip mean of the hydrometeor in the 1st PDF component.
            imu_hm_1(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=imu_hm_1(hm_idx), &
                                 var_name="mu_"//trim( hm_type(1:2) )//"_1", &
                                 var_description="Mean (in-precip) of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (1st PDF component) [kg/kg]", &
                                 var_units="kg/kg", &
                                 l_silhs=.false., grid_kind=zt )

            else ! Concentration

               call stat_assign( var_index=imu_hm_1(hm_idx), &
                                 var_name="mu_"//trim( hm_type(1:2) )//"_1", &
                                 var_description="Mean (in-precip) of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (1st PDF component) [num/kg]", &
                                 var_units="num/kg", &
                                 l_silhs=.false., grid_kind=zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

            ! The in-precip mean of the hydrometeor in the 2nd PDF component.
            imu_hm_2(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=imu_hm_2(hm_idx), &
                                 var_name="mu_"//trim( hm_type(1:2) )//"_2", &
                                 var_description="Mean (in-precip) of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (2nd PDF component) [kg/kg]", &
                                 var_units="kg/kg", &
                                 l_silhs=.false., grid_kind=zt )

            else ! Concentration

               call stat_assign( var_index=imu_hm_2(hm_idx), &
                                 var_name="mu_"//trim( hm_type(1:2) )//"_2", &
                                 var_description="Mean (in-precip) of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (2nd PDF component) [num/kg]", &
                                 var_units="num/kg", &
                                 l_silhs=.false., grid_kind=zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'mu_Ncn_1' )
        imu_Ncn_1 = k
        call stat_assign( var_index=imu_Ncn_1, var_name="mu_Ncn_1", &
             var_description="Mean of N_cn (1st PDF component) [num/kg]", var_units="num/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'mu_Ncn_2' )
        imu_Ncn_2 = k
        call stat_assign( var_index=imu_Ncn_2, var_name="mu_Ncn_2", &
             var_description="Mean of N_cn (2nd PDF component) [num/kg]", var_units="num/kg", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      ! Hydrometeor component mean values (in-precip) for ln hm for each PDF
      ! component and hydrometeor type.
      case ( 'mu_hm_i_n' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip mean of ln hm in the 1st PDF component.
            imu_hm_1_n(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=imu_hm_1_n(hm_idx), &
                                 var_name="mu_"//trim( hm_type(1:2) )//"_1_n", &
                                 var_description="Mean (in-precip) of ln " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (1st PDF component) [ln(kg/kg)]", &
                                 var_units="ln(kg/kg)", &
                                 l_silhs=.false., grid_kind=zt )

            else ! Concentration

               call stat_assign( var_index=imu_hm_1_n(hm_idx), &
                                 var_name="mu_"//trim( hm_type(1:2) )//"_1_n", &
                                 var_description="Mean (in-precip) of ln " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (1st PDF component) [ln(num/kg)]", &
                                 var_units="ln(num/kg)", &
                                 l_silhs=.false., grid_kind=zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

            ! The in-precip mean of ln hm in the 2nd PDF component.
            imu_hm_2_n(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=imu_hm_2_n(hm_idx), &
                                 var_name="mu_"//trim( hm_type(1:2) )//"_2_n", &
                                 var_description="Mean (in-precip) of ln " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (2nd PDF component) [ln(kg/kg)]", &
                                 var_units="ln(kg/kg)", &
                                 l_silhs=.false., grid_kind=zt )

            else ! Concentration

               call stat_assign( var_index=imu_hm_2_n(hm_idx), &
                                 var_name="mu_"//trim( hm_type(1:2) )//"_2_n", &
                                 var_description="Mean (in-precip) of ln " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (2nd PDF component) [ln(num/kg)]", &
                                 var_units="ln(num/kg)", &
                                 l_silhs=.false., grid_kind=zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'mu_Ncn_1_n' )
        imu_Ncn_1_n = k
        call stat_assign( var_index=imu_Ncn_1_n, var_name="mu_Ncn_1_n", &
             var_description="Mean of ln N_cn (1st PDF component) [ln(num/kg)]", &
             var_units="ln(num/kg)", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'mu_Ncn_2_n' )
        imu_Ncn_2_n = k
        call stat_assign( var_index=imu_Ncn_2_n, var_name="mu_Ncn_2_n", &
             var_description="Mean of ln N_cn (2nd PDF component) [ln(num/kg)]", &
             var_units="ln(num/kg)", l_silhs=.false., grid_kind=zt )
        k = k + 1

      ! Hydrometeor component standard deviations (in-precip) for each PDF
      ! component and hydrometeor type.
      case ( 'sigma_hm_i' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip standard deviation of the hydrometeor in the 1st PDF
            ! component.
            isigma_hm_1(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=isigma_hm_1(hm_idx), &
                                 var_name="sigma_" &
                                 // trim( hm_type(1:2) )//"_1", &
                                 var_description="Standard deviation " &
                                 // "(in-precip) of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (1st PDF component) [kg/kg]", &
                                 var_units="kg/kg", &
                                 l_silhs=.false., grid_kind=zt )

            else ! Concentration

               call stat_assign( var_index=isigma_hm_1(hm_idx), &
                                 var_name="sigma_" &
                                 // trim( hm_type(1:2) )//"_1", &
                                 var_description="Standard deviation " &
                                 // "(in-precip) of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (1st PDF component) [num/kg]", &
                                 var_units="num/kg", &
                                 l_silhs=.false., grid_kind=zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

            ! The in-precip standard deviation of the hydrometeor in the 2nd PDF
            ! component.
            isigma_hm_2(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=isigma_hm_2(hm_idx), &
                                 var_name="sigma_" &
                                 // trim( hm_type(1:2) )//"_2", &
                                 var_description="Standard deviation " &
                                 // "(in-precip) of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (2nd PDF component) [kg/kg]", &
                                 var_units="kg/kg", &
                                 l_silhs=.false., grid_kind=zt )

            else ! Concentration

               call stat_assign( var_index=isigma_hm_2(hm_idx), &
                                 var_name="sigma_" &
                                 // trim( hm_type(1:2) )//"_2", &
                                 var_description="Standard deviation " &
                                 // "(in-precip) of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (2nd PDF component) [num/kg]", &
                                 var_units="num/kg", &
                                 l_silhs=.false., grid_kind=zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'sigma_Ncn_1' )
        isigma_Ncn_1 = k
        call stat_assign( var_index=isigma_Ncn_1, var_name="sigma_Ncn_1", &
             var_description="Standard deviation of N_cn (1st PDF component) [num/kg]", &
             var_units="num/kg", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'sigma_Ncn_2' )
        isigma_Ncn_2 = k
        call stat_assign( var_index=isigma_Ncn_2, var_name="sigma_Ncn_2", &
             var_description="Standard deviation of N_cn (2nd PDF component) [num/kg]", &
             var_units="num/kg", l_silhs=.false., grid_kind=zt )
        k = k + 1

      ! Hydrometeor component standard deviations (in-precip) for ln hm for each
      ! PDF component and hydrometeor type.
      case ( 'sigma_hm_i_n' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip standard deviation of ln hm in the 1st PDF
            ! component.
            isigma_hm_1_n(hm_idx) = k

            call stat_assign( var_index=isigma_hm_1_n(hm_idx), &
                              var_name="sigma_" &
                              // trim( hm_type(1:2) )//"_1_n", &
                              var_description="Standard deviation " &
                              // "(in-precip) of ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (1st PDF component) [-]", &
                              var_units="-", &
                              l_silhs=.false., grid_kind=zt )

            k = k + 1

            ! The in-precip standard deviation of ln hm in the 2nd PDF
            ! component.
            isigma_hm_2_n(hm_idx) = k

            call stat_assign( var_index=isigma_hm_2_n(hm_idx), &
                              var_name="sigma_" &
                              // trim( hm_type(1:2) )//"_2_n", &
                              var_description="Standard deviation " &
                              // "(in-precip) of ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (2nd PDF component) [-]", &
                              var_units="-", &
                              l_silhs=.false., grid_kind=zt )

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'sigma_Ncn_1_n' )
        isigma_Ncn_1_n = k
        call stat_assign( var_index=isigma_Ncn_1_n, var_name="sigma_Ncn_1_n", &
             var_description="Standard deviation of ln N_cn (1st PDF component) [ln(num/kg)]", &
             var_units="ln(num/kg)", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'sigma_Ncn_2_n' )
        isigma_Ncn_2_n = k
        call stat_assign( var_index=isigma_Ncn_2_n, var_name="sigma_Ncn_2_n", &
             var_description="Standard deviation of ln N_cn (2nd PDF component) [ln(num/kg)]", &
             var_units="ln(num/kg)", l_silhs=.false., grid_kind=zt )
        k = k + 1

      ! Correlation between w and a hydrometeor (in-precip) for each PDF
      ! component and hydrometeor type.
      case ( 'corr_whm_i' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip correlation between w and the hydrometeor in the
            ! 1st PDF component.
            icorr_whm_1(hm_idx) = k

            call stat_assign( var_index=icorr_whm_1(hm_idx), &
                              var_name="corr_w"//trim( hm_type(1:2) )//"_1", &
                              var_description="Correlation (in-precip) " &
                              // "between w and " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (1st PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=zt )

            k = k + 1

            ! The in-precip correlation between w and the hydrometeor in the
            ! 2nd PDF component.
            icorr_whm_2(hm_idx) = k

            call stat_assign( var_index=icorr_whm_2(hm_idx), &
                              var_name="corr_w"//trim( hm_type(1:2) )//"_2", &
                              var_description="Correlation (in-precip) " &
                              // "between w and " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (2nd PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=zt )

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'corr_wNcn_1' )
        icorr_wNcn_1 = k
        call stat_assign( var_index=icorr_wNcn_1, var_name="corr_wNcn_1", &
             var_description="Correlation between w and N_cn (1st PDF component) [-]", &
             var_units="-", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'corr_wNcn_2' )
        icorr_wNcn_2 = k
        call stat_assign( var_index=icorr_wNcn_2, var_name="corr_wNcn_2", &
             var_description="Correlation between w and N_cn (2nd PDF component) [-]", &
             var_units="-", l_silhs=.false., grid_kind=zt )
        k = k + 1

      ! Correlation between s and a hydrometeor (in-precip) for each PDF
      ! component and hydrometeor type.
      case ( 'corr_shm_i' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip correlation between s and the hydrometeor in the
            ! 1st PDF component.
            icorr_shm_1(hm_idx) = k

            call stat_assign( var_index=icorr_shm_1(hm_idx), &
                              var_name="corr_s"//trim( hm_type(1:2) )//"_1", &
                              var_description="Correlation (in-precip) " &
                              // "between s and " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (1st PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=zt )

            k = k + 1

            ! The in-precip correlation between s and the hydrometeor in the
            ! 2nd PDF component.
            icorr_shm_2(hm_idx) = k

            call stat_assign( var_index=icorr_shm_2(hm_idx), &
                              var_name="corr_s"//trim( hm_type(1:2) )//"_2", &
                              var_description="Correlation (in-precip) " &
                              // "between s and " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (2nd PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=zt )

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'corr_sNcn_1' )
        icorr_sNcn_1 = k
        call stat_assign( var_index=icorr_sNcn_1, var_name="corr_sNcn_1", &
             var_description="Correlation between s and N_cn (1st PDF component) [-]", &
             var_units="-", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'corr_sNcn_2' )
        icorr_sNcn_2 = k
        call stat_assign( var_index=icorr_sNcn_2, var_name="corr_sNcn_2", &
             var_description="Correlation between s and N_cn (2nd PDF component) [-]", &
             var_units="-", l_silhs=.false., grid_kind=zt )
        k = k + 1

      ! Correlation between t and a hydrometeor (in-precip) for each PDF
      ! component and hydrometeor type.
      case ( 'corr_thm_i' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip correlation between t and the hydrometeor in the
            ! 1st PDF component.
            icorr_thm_1(hm_idx) = k

            call stat_assign( var_index=icorr_thm_1(hm_idx), &
                              var_name="corr_t"//trim( hm_type(1:2) )//"_1", &
                              var_description="Correlation (in-precip) " &
                              // "between t and " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (1st PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=zt )

            k = k + 1

            ! The in-precip correlation between t and the hydrometeor in the
            ! 2nd PDF component.
            icorr_thm_2(hm_idx) = k

            call stat_assign( var_index=icorr_thm_2(hm_idx), &
                              var_name="corr_t"//trim( hm_type(1:2) )//"_2", &
                              var_description="Correlation (in-precip) " &
                              // "between t and " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (2nd PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=zt )

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'corr_tNcn_1' )
        icorr_tNcn_1 = k
        call stat_assign( var_index=icorr_tNcn_1, var_name="corr_tNcn_1", &
             var_description="Correlation between t and N_cn (1st PDF component) [-]", &
             var_units="-", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'corr_tNcn_2' )
        icorr_tNcn_2 = k
        call stat_assign( var_index=icorr_tNcn_2, var_name="corr_tNcn_2", &
             var_description="Correlation between t and N_cn (2nd PDF component) [-]", &
             var_units="-", l_silhs=.false., grid_kind=zt )
        k = k + 1

      ! Correlation between Ncn and a hydrometeor (in-precip) for each PDF
      ! component and hydrometeor type.
      case ( 'corr_Ncnhm_i' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip correlation between Ncn and the hydrometeor in the
            ! 1st PDF component.
            icorr_Ncnhm_1(hm_idx) = k

            call stat_assign( var_index=icorr_Ncnhm_1(hm_idx), &
                              var_name="corr_Ncn"//trim( hm_type(1:2) )//"_1", &
                              var_description="Correlation (in-precip) " &
                              // "between N_cn and " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (1st PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=zt )

            k = k + 1

            ! The in-precip correlation between Ncn and the hydrometeor in the
            ! 2nd PDF component.
            icorr_Ncnhm_2(hm_idx) = k

            call stat_assign( var_index=icorr_Ncnhm_2(hm_idx), &
                              var_name="corr_Ncn"//trim( hm_type(1:2) )//"_2", &
                              var_description="Correlation (in-precip) " &
                              // "between N_cn and " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (2nd PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=zt )

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      ! Correlation (in-precip) between two different hydrometeors (hmx and hmy)
      ! for each PDF component and hydrometeor type.
      case ( 'corr_hmxhmy_i' )

         do hmx_idx = 1, hydromet_dim, 1

            hmx_type = hydromet_list(hmx_idx)

            do hmy_idx = hmx_idx+1, hydromet_dim, 1

               hmy_type = hydromet_list(hmy_idx)

               ! The in-precip correlation between hmx and hmy in the 1st PDF
               ! component.
               icorr_hmxhmy_1(hmy_idx,hmx_idx) = k

               call stat_assign( var_index=icorr_hmxhmy_1(hmy_idx,hmx_idx), &
                                 var_name="corr_"//trim( hmx_type(1:2) ) &
                                 // trim( hmy_type(1:2) )//"_1", &
                                 var_description="Correlation (in-precip) " &
                                 // "between " &
                                 // hmx_type(1:1)//"_"//trim( hmx_type(2:2) ) &
                                 // " and " &
                                 // hmy_type(1:1)//"_"//trim( hmy_type(2:2) ) &
                                 // " (1st PDF component) [-]", &
                                 var_units="-", l_silhs=.false., grid_kind=zt )

               k = k + 1

               ! The in-precip correlation between hmx and hmy in the 2nd PDF
               ! component.
               icorr_hmxhmy_2(hmy_idx,hmx_idx) = k

               call stat_assign( var_index=icorr_hmxhmy_2(hmy_idx,hmx_idx), &
                                 var_name="corr_"//trim( hmx_type(1:2) ) &
                                 // trim( hmy_type(1:2) )//"_2", &
                                 var_description="Correlation (in-precip) " &
                                 // "between " &
                                 // hmx_type(1:1)//"_"//trim( hmx_type(2:2) ) &
                                 // " and " &
                                 // hmy_type(1:1)//"_"//trim( hmy_type(2:2) ) &
                                 // " (2nd PDF component) [-]", &
                                 var_units="-", l_silhs=.false., grid_kind=zt )

               k = k + 1

            enddo ! hmy_idx = hmx_idx+1, hydromet_dim, 1

         enddo ! hmx_idx = 1, hydromet_dim, 1

      ! Correlation (in-precip) between w and ln hm for each PDF component and
      ! hydrometeor type.
      case ( 'corr_whm_i_n' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip correlation between w and ln hm in the 1st PDF
            ! component.
            icorr_whm_1_n(hm_idx) = k

            call stat_assign( var_index=icorr_whm_1_n(hm_idx), &
                              var_name="corr_w"//trim( hm_type(1:2) )//"_1_n", &
                              var_description="Correlation (in-precip) " &
                              // "between w and ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (1st PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=zt )

            k = k + 1

            ! The in-precip correlation between w and ln hm in the 2nd PDF
            ! component.
            icorr_whm_2_n(hm_idx) = k

            call stat_assign( var_index=icorr_whm_2_n(hm_idx), &
                              var_name="corr_w"//trim( hm_type(1:2) )//"_2_n", &
                              var_description="Correlation (in-precip) " &
                              // "between w and ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (2nd PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=zt )

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'corr_wNcn_1_n' )
        icorr_wNcn_1_n = k
        call stat_assign( var_index=icorr_wNcn_1_n, var_name="corr_wNcn_1_n", &
             var_description="Correlation between w and ln N_cn (1st PDF component) [-]", &
             var_units="-", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'corr_wNcn_2_n' )
        icorr_wNcn_2_n = k
        call stat_assign( var_index=icorr_wNcn_2_n, var_name="corr_wNcn_2_n", &
             var_description="Correlation between w and ln N_cn (2nd PDF component) [-]", &
             var_units="-", l_silhs=.false., grid_kind=zt )
        k = k + 1

      ! Correlation (in-precip) between s and ln hm for each PDF component and
      ! hydrometeor type.
      case ( 'corr_shm_i_n' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip correlation between s and ln hm in the 1st PDF
            ! component.
            icorr_shm_1_n(hm_idx) = k

            call stat_assign( var_index=icorr_shm_1_n(hm_idx), &
                              var_name="corr_s"//trim( hm_type(1:2) )//"_1_n", &
                              var_description="Correlation (in-precip) " &
                              // "between s and ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (1st PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=zt )

            k = k + 1

            ! The in-precip correlation between s and ln hm in the 2nd PDF
            ! component.
            icorr_shm_2_n(hm_idx) = k

            call stat_assign( var_index=icorr_shm_2_n(hm_idx), &
                              var_name="corr_s"//trim( hm_type(1:2) )//"_2_n", &
                              var_description="Correlation (in-precip) " &
                              // "between s and ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (2nd PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=zt )

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'corr_sNcn_1_n' )
        icorr_sNcn_1_n = k
        call stat_assign( var_index=icorr_sNcn_1_n, var_name="corr_sNcn_1_n", &
             var_description="Correlation between s and ln N_cn (1st PDF component) [-]", &
             var_units="-", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'corr_sNcn_2_n' )
        icorr_sNcn_2_n = k
        call stat_assign( var_index=icorr_sNcn_2_n, var_name="corr_sNcn_2_n", &
             var_description="Correlation between s and ln N_cn (2nd PDF component) [-]", &
             var_units="-", l_silhs=.false., grid_kind=zt )
        k = k + 1

      ! Correlation (in-precip) between t and ln hm for each PDF component and
      ! hydrometeor type.
      case ( 'corr_thm_i_n' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip correlation between t and ln hm in the 1st PDF
            ! component.
            icorr_thm_1_n(hm_idx) = k

            call stat_assign( var_index=icorr_thm_1_n(hm_idx), &
                              var_name="corr_t"//trim( hm_type(1:2) )//"_1_n", &
                              var_description="Correlation (in-precip) " &
                              // "between t and ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (1st PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=zt )

            k = k + 1

            ! The in-precip correlation between t and ln hm in the 2nd PDF
            ! component.
            icorr_thm_2_n(hm_idx) = k

            call stat_assign( var_index=icorr_thm_2_n(hm_idx), &
                              var_name="corr_t"//trim( hm_type(1:2) )//"_2_n", &
                              var_description="Correlation (in-precip) " &
                              // "between t and ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (2nd PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=zt )

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'corr_tNcn_1_n' )
        icorr_tNcn_1_n = k
        call stat_assign( var_index=icorr_tNcn_1_n, var_name="corr_tNcn_1_n", &
             var_description="Correlation between t and ln N_cn (1st PDF component) [-]", &
             var_units="-", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'corr_tNcn_2_n' )
        icorr_tNcn_2_n = k
        call stat_assign( var_index=icorr_tNcn_2_n, var_name="corr_tNcn_2_n", &
             var_description="Correlation between t and ln N_cn (2nd PDF component) [-]", &
             var_units="-", l_silhs=.false., grid_kind=zt )
        k = k + 1

      ! Correlation (in-precip) between ln Ncn and ln hm for each PDF component
      ! and hydrometeor type.
      case ( 'corr_Ncnhm_i_n' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip correlation between ln Ncn and ln hm in the 1st PDF
            ! component.
            icorr_Ncnhm_1_n(hm_idx) = k

            call stat_assign( var_index=icorr_Ncnhm_1_n(hm_idx), &
                              var_name="corr_Ncn"//trim(hm_type(1:2))//"_1_n", &
                              var_description="Correlation (in-precip) " &
                              // "between ln N_cn and ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (1st PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=zt )

            k = k + 1

            ! The in-precip correlation between ln Ncn and ln hm in the 2nd PDF
            ! component.
            icorr_Ncnhm_2_n(hm_idx) = k

            call stat_assign( var_index=icorr_Ncnhm_2_n(hm_idx), &
                              var_name="corr_Ncn"//trim(hm_type(1:2))//"_2_n", &
                              var_description="Correlation (in-precip) " &
                              // "between ln N_cn and ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (2nd PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=zt )

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      ! Correlation (in-precip) between ln hmx and ln hmy (hmx and hmy are two
      ! different hydrometeors) for each PDF component and hydrometeor type.
      case ( 'corr_hmxhmy_i_n' )

         do hmx_idx = 1, hydromet_dim, 1

            hmx_type = hydromet_list(hmx_idx)

            do hmy_idx = hmx_idx+1, hydromet_dim, 1

               hmy_type = hydromet_list(hmy_idx)

               ! The in-precip correlation between ln hmx and ln hmy in the 1st
               ! PDF component.
               icorr_hmxhmy_1_n(hmy_idx,hmx_idx) = k

               call stat_assign( var_index=icorr_hmxhmy_1_n(hmy_idx,hmx_idx), &
                                 var_name="corr_"//trim( hmx_type(1:2) ) &
                                 // trim( hmy_type(1:2) )//"_1_n", &
                                 var_description="Correlation (in-precip) " &
                                 // "between ln " &
                                 // hmx_type(1:1)//"_"//trim( hmx_type(2:2) ) &
                                 // " and ln " &
                                 // hmy_type(1:1)//"_"//trim( hmy_type(2:2) ) &
                                 // " (1st PDF component) [-]", &
                                 var_units="-", l_silhs=.false., grid_kind=zt )

               k = k + 1

               ! The in-precip correlation between ln hmx and ln hmy in the 2nd
               ! PDF component.
               icorr_hmxhmy_2_n(hmy_idx,hmx_idx) = k

               call stat_assign( var_index=icorr_hmxhmy_2_n(hmy_idx,hmx_idx), &
                                 var_name="corr_"//trim( hmx_type(1:2) ) &
                                 // trim( hmy_type(1:2) )//"_2_n", &
                                 var_description="Correlation (in-precip) " &
                                 // "between ln " &
                                 // hmx_type(1:1)//"_"//trim( hmx_type(2:2) ) &
                                 // " and ln " &
                                 // hmy_type(1:1)//"_"//trim( hmy_type(2:2) ) &
                                 // " (2nd PDF component) [-]", &
                                 var_units="-", l_silhs=.false., grid_kind=zt )

               k = k + 1

            enddo ! hmy_idx = hmx_idx+1, hydromet_dim, 1

         enddo ! hmx_idx = 1, hydromet_dim, 1

      ! changes by janhft 09/25/12
      case ('corr_sw')
        icorr_sw = k
        call stat_assign( var_index=icorr_sw, var_name="corr_sw", &
             var_description="Correlation between s and w [-]", var_units="-", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ( 'corr_srr' )
        icorr_srr = k
        call stat_assign( var_index=icorr_srr, var_name="corr_srr", &
             var_description="Correlation (in-precip) between s and r_r (corr array) [-]", &
             var_units="-", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'corr_sNr' )
        icorr_sNr = k
        call stat_assign( var_index=icorr_sNr, var_name="corr_sNr", &
             var_description="Correlation (in-precip) between s and N_r (corr array) [-]", &
             var_units="-", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'corr_sNcn' )
        icorr_sNcn = k
        call stat_assign( var_index=icorr_sNcn, var_name="corr_sNcn", &
             var_description="Correlation between s and N_cn (corr array) [-]", var_units="-", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ( 'corr_rrNr' )
        icorr_rrNr = k
        call stat_assign( var_index=icorr_rrNr, var_name="corr_rrNr", &
             var_description="Correlation (in-precip) between r_r and N_r (corr array) [-]", &
             var_units="-", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('corr_wrr')
        icorr_wrr = k
        call stat_assign( var_index=icorr_wrr, var_name="corr_wrr", &
             var_description="Correlation between w and rrain [-]", var_units="-", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('corr_wNr')
        icorr_wNr = k
        call stat_assign( var_index=icorr_wNr, var_name="corr_wNr", &
             var_description="Correlation between w and Nr [-]", var_units="-", l_silhs=.false., &
             grid_kind=zt )
        k = k + 1

      case ('corr_wNcn')
        icorr_wNcn = k
        call stat_assign( var_index=icorr_wNcn, var_name="corr_wNcn", &
             var_description="Correlation between w and N_cn [-]", var_units="-", &
             l_silhs=.false., grid_kind=zt )
        k = k + 1
      ! end changes by janhft 09/25/12

      case ('wp2hmp')
        iwp2hmp = k
        call stat_assign( var_index=iwp2hmp, var_name="wp2hmp", &
             var_description="Third moment: <w'^2> * <hydro.'> [(m/s)^2 <hydrometeor units>]", &
             var_units="(m/s)^2 <hydrometeor units>", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('cloud_frac_refined')
        icloud_frac_refined = k
        call stat_assign( var_index=icloud_frac_refined, var_name="cloud_frac_refined", &
                          var_description="Cloud fraction computed on refined grid [-]", &
                          var_units="-", l_silhs=.false., grid_kind=zt )
        k = k + 1

      case ('rcm_refined')
        ircm_refined = k
        call stat_assign( var_index=ircm_refined, var_name="rcm_refined", &
                          var_description="Cloud water mixing ratio computed on refined grid &
                          &[kg/kg]", var_units="kg/kg", l_silhs=.false., grid_kind=zt)
        k = k + 1

      case ('hl_on_Cp_residual')
        ihl_on_Cp_residual = k
        call stat_assign( var_index=ihl_on_Cp_residual, var_name="hl_on_Cp_residual", &
                          var_description="Residual change in HL/Cp from Morrison microphysics &
                          &not due to sedimentation [K]", &
                          var_units="K", l_silhs=.true., grid_kind=zt)
        k = k + 1

      case ('qto_residual')
        iqto_residual = k
        call stat_assign( var_index=iqto_residual, var_name="qto_residual", &
                          var_description="Residual change in total water from Morrison &
                          &microphysics not due to sedimentation [kg/kg]", &
                          var_units="kg/kg", l_silhs=.true., grid_kind=zt)
        k = k + 1


      case default

        l_found = .false.

        j = 1

        do while( j <= sclr_dim .and. .not. l_found)
          write(sclr_idx, * ) j

          sclr_idx = adjustl(sclr_idx)

          if(trim(vars_zt(i)) == "sclr"//trim(sclr_idx)//"m" .and. .not. l_found) then

            isclrm(j) = k

        call stat_assign( var_index=isclrm(j), var_name="sclr"//trim(sclr_idx)//"m", &
             var_description="passive scalar "//trim(sclr_idx), var_units="unknown", &
             l_silhs=.false., grid_kind=zt )

            k = k + 1

            l_found = .true.

          else if(trim(vars_zt(i)) == "sclr"//trim(sclr_idx)//"m_f" .and. .not. l_found) then

            isclrm_f(j) = k

        call stat_assign( var_index=isclrm_f(j), var_name="sclr"//trim(sclr_idx)//"m_f", &
             var_description="passive scalar forcing "//trim(sclr_idx), var_units="unknown", &
             l_silhs=.false., grid_kind=zt )

            k = k + 1

            l_found = .true.

          endif

          j = j + 1
        end do

        j = 1

        do while( j <= edsclr_dim .and. .not. l_found)

          write(sclr_idx, * ) j

          sclr_idx = adjustl(sclr_idx)

          if(trim(vars_zt(i)) == "edsclr"//trim(sclr_idx)//"m" .and. .not. l_found ) then

            iedsclrm(j) = k

        call stat_assign( var_index=iedsclrm(j), var_name="edsclr"//trim(sclr_idx)//"m", &
             var_description="passive scalar "//trim(sclr_idx), var_units="unknown", &
             l_silhs=.false., grid_kind=zt )

            k = k + 1

            l_found = .true.

          else if(trim(vars_zt(i)) == "edsclr"//trim(sclr_idx)//"m_f" .and. .not. l_found) then

            iedsclrm_f(j) = k

        call stat_assign( var_index=iedsclrm_f(j), var_name="edsclr"//trim(sclr_idx)//"m_f", &
             var_description="passive scalar forcing "//trim(sclr_idx), var_units="unknown", &
             l_silhs=.false., grid_kind=zt )

            k = k + 1

            l_found = .true.

          endif

          j = j + 1

        end do

        if (.not. l_found ) then

          write(fstderr,*) 'Error:  unrecognized variable in vars_zt:  ', trim( vars_zt(i) )

          l_error = .true.  ! This will stop the run.

        end if

      end select ! trim( vars_zt )

    end do ! i=1,zt%nn
    return
  end subroutine stats_init_zt

end module stats_zt
