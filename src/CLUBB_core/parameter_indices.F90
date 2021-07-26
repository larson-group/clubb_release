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
    nparams = 99 ! Total tunable parameters

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
    iC_wp3_pr_turb    = 32, &
    iC_wp3_pr_dfsn    = 33, &
    iC_wp2_splat      = 34

  integer, parameter, public :: &
    iC6rt_Lscale0     = 35, &
    iC6thl_Lscale0    = 36, &
    iC7_Lscale0       = 37, &
    iwpxp_L_thresh    = 38

  integer, parameter, public :: & 
    ic_K              = 39, & 
    ic_K1             = 40, & 
    inu1              = 41, & 
    ic_K2             = 42, & 
    inu2              = 43, & 
    ic_K6             = 44, & 
    inu6              = 45, & 
    ic_K8             = 46, & 
    inu8              = 47, & 
    ic_K9             = 48, & 
    inu9              = 49, & 
    inu10             = 50, &
    ic_K_hm           = 51, & 
    ic_K_hmb          = 52, & 
    iK_hm_min_coef    = 53, & 
    inu_hm            = 54 

  integer, parameter, public :: &
    islope_coef_spread_DG_means_w = 55, &
    ipdf_component_stdev_factor_w = 56, &
    icoef_spread_DG_means_rt      = 57, &
    icoef_spread_DG_means_thl     = 58, &
    igamma_coef                   = 59, & 
    igamma_coefb                  = 60, & 
    igamma_coefc                  = 61, & 
    imu                           = 62, & 
    ibeta                         = 63, & 
    ilmin_coef                    = 64, &
    iomicron                      = 65, &
    izeta_vrnce_rat               = 66, &
    iupsilon_precip_frac_rat      = 67, &
    ilambda0_stability_coef       = 68, &
    imult_coef                    = 69, &
    itaumin                       = 70, &
    itaumax                       = 71, &
    iLscale_mu_coef               = 72, &
    iLscale_pert_coef             = 73, &
    ialpha_corr                   = 74, &
    iSkw_denom_coef               = 75, &
    ic_K10                        = 76, &
    ic_K10h                       = 77, &
    ithlp2_rad_coef               = 78, &
    ithlp2_rad_cloud_frac_thresh  = 79, &
    iup2_sfc_coef                 = 80, &
    iSkw_max_mag                  = 81, &
    iC_invrs_tau_bkgnd            = 82, &
    iC_invrs_tau_sfc              = 83, &
    iC_invrs_tau_shear            = 84, &
    iC_invrs_tau_N2               = 85, &
    iC_invrs_tau_N2_wp2           = 86, &
    iC_invrs_tau_N2_xp2           = 87, &
    iC_invrs_tau_N2_wpxp          = 88, &
    iC_invrs_tau_N2_clear_wp3     = 89, &
    iC_invrs_tau_wpxp_Ri          = 90, &
    iC_invrs_tau_wpxp_N2_thresh   = 91, &
    ixp3_coef_base                = 92, &
    ixp3_coef_slope               = 93, &
    ialtitude_threshold           = 94, &
    irtp2_clip_coef               = 95, &
    iCx_min                       = 96, &
    iCx_max                       = 97, &
    iRichardson_num_min           = 98, &
    iRichardson_num_max           = 99

end module parameter_indices
!-----------------------------------------------------------------------
