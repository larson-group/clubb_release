!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module parameter_indices

! Description:
!   Since f90/95 lacks enumeration, we're stuck numbering each
!   parameter by hand like this.

!   Adding new parameters is relatively simple.  First, the
!   parameter should be added in the common block of the parameters
!   module so it can be used in other parts of the code. Each
!   variable needs a unique number in this module, and nparams must
!   be incremented for the new variable.  Next, the params_list
!   variable in module parameters should have new variable added to
!   it.  The subroutines pack_parameters and uppack_parameters will
!   need to have the variable added to their list, but the order
!   doesn't actually matter, since the i variables in here determine
!   where in the params vector the number is placed.
!   Finally, the namelists clubb_params_nl and initspread will need to
!   have the parameter added to them.
!-------------------------------------------------------------------------------

  implicit none

  private ! Default Scope

  integer, parameter, public ::  & 
    nparams = 100 ! Total tunable parameters

!***************************************************************
!                    ***** IMPORTANT *****
! If you change the order of these parameters, you will need to
! change the order of params_list as well or the tuner will
! break!
!                    ***** IMPORTANT *****
!***************************************************************

  integer, parameter, public :: & 
    iC1               =  1, & 
    iC1b              =  2, & 
    iC1c              =  3, & 
    iC2               =  4, & 
    iC2b              =  5, & 
    iC2c              =  6, & 
    iC2rt             =  7, & 
    iC2thl            =  8, & 
    iC2rtthl          =  9, & 
    iC4               = 10, & 
    iC_uu_shr         = 11, &
    iC_uu_buoy        = 12, & 
    iC6rt             = 13, & 
    iC6rtb            = 14, & 
    iC6rtc            = 15, & 
    iC6thl            = 16, & 
    iC6thlb           = 17, & 
    iC6thlc           = 18, & 
    iC7               = 19, & 
    iC7b              = 20, & 
    iC7c              = 21, & 
    iC8               = 22, & 
    iC8b              = 23, & 
    iC10              = 24, & 
    iC11              = 25, & 
    iC11b             = 26, & 
    iC11c             = 27, & 
    iC12              = 28, & 
    iC13              = 29, & 
    iC14              = 30, &
    iC_wp2_pr_dfsn    = 31, &
    iC_wp3_pr_tp      = 32, &
    iC_wp3_pr_turb    = 33, &
    iC_wp3_pr_dfsn    = 34, &
    iC_wp2_splat      = 35

  integer, parameter, public :: &
    iC6rt_Lscale0     = 36, &
    iC6thl_Lscale0    = 37, &
    iC7_Lscale0       = 38, &
    iwpxp_L_thresh    = 39

  integer, parameter, public :: & 
    ic_K              = 40, & 
    ic_K1             = 41, & 
    inu1              = 42, & 
    ic_K2             = 43, & 
    inu2              = 44, & 
    ic_K6             = 45, & 
    inu6              = 46, & 
    ic_K8             = 47, & 
    inu8              = 48, & 
    ic_K9             = 49, & 
    inu9              = 50, & 
    inu10             = 51, &
    ic_K_hm           = 52, & 
    ic_K_hmb          = 53, & 
    iK_hm_min_coef    = 54, & 
    inu_hm            = 55 

  integer, parameter, public :: &
    islope_coef_spread_DG_means_w = 56, &
    ipdf_component_stdev_factor_w = 57, &
    icoef_spread_DG_means_rt      = 58, &
    icoef_spread_DG_means_thl     = 59, &
    igamma_coef                   = 60, & 
    igamma_coefb                  = 61, & 
    igamma_coefc                  = 62, & 
    imu                           = 63, & 
    ibeta                         = 64, & 
    ilmin_coef                    = 65, &
    iomicron                      = 66, &
    izeta_vrnce_rat               = 67, &
    iupsilon_precip_frac_rat      = 68, &
    ilambda0_stability_coef       = 69, &
    imult_coef                    = 70, &
    itaumin                       = 71, &
    itaumax                       = 72, &
    iLscale_mu_coef               = 73, &
    iLscale_pert_coef             = 74, &
    ialpha_corr                   = 75, &
    iSkw_denom_coef               = 76, &
    ic_K10                        = 77, &
    ic_K10h                       = 78, &
    ithlp2_rad_coef               = 79, &
    ithlp2_rad_cloud_frac_thresh  = 80, &
    iup2_sfc_coef                 = 81, &
    iSkw_max_mag                  = 82, &
    iC_invrs_tau_bkgnd            = 83, &
    iC_invrs_tau_sfc              = 84, &
    iC_invrs_tau_shear            = 85, &
    iC_invrs_tau_N2               = 86, &
    iC_invrs_tau_N2_wp2           = 87, &
    iC_invrs_tau_N2_xp2           = 88, &
    iC_invrs_tau_N2_wpxp          = 89, &
    iC_invrs_tau_N2_clear_wp3     = 90, &
    iC_invrs_tau_wpxp_Ri          = 91, &
    iC_invrs_tau_wpxp_N2_thresh   = 92, &
    ixp3_coef_base                = 93, &
    ixp3_coef_slope               = 94, &
    ialtitude_threshold           = 95, &
    irtp2_clip_coef               = 96, &
    iCx_min                       = 97, &
    iCx_max                       = 98, &
    iRichardson_num_min           = 99, &
    iRichardson_num_max           = 100

end module parameter_indices
!-----------------------------------------------------------------------
