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
    nparams = 98 ! Total tunable parameters

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
    iC_wp3_pr_turb    = 31, &
    iC_wp3_pr_dfsn    = 32, &
    iC_wp2_splat      = 33

  integer, parameter, public :: &
    iC6rt_Lscale0     = 34, &
    iC6thl_Lscale0    = 35, &
    iC7_Lscale0       = 36, &
    iwpxp_L_thresh    = 37

  integer, parameter, public :: & 
    ic_K              = 38, & 
    ic_K1             = 39, & 
    inu1              = 40, & 
    ic_K2             = 41, & 
    inu2              = 42, & 
    ic_K6             = 43, & 
    inu6              = 44, & 
    ic_K8             = 45, & 
    inu8              = 46, & 
    ic_K9             = 47, & 
    inu9              = 48, & 
    inu10             = 49, &
    ic_K_hm           = 50, & 
    ic_K_hmb          = 51, & 
    iK_hm_min_coef    = 52, & 
    inu_hm            = 53 

  integer, parameter, public :: &
    islope_coef_spread_DG_means_w = 54, &
    ipdf_component_stdev_factor_w = 55, &
    icoef_spread_DG_means_rt      = 56, &
    icoef_spread_DG_means_thl     = 57, &
    igamma_coef                   = 58, & 
    igamma_coefb                  = 59, & 
    igamma_coefc                  = 60, & 
    imu                           = 61, & 
    ibeta                         = 62, & 
    ilmin_coef                    = 63, &
    iomicron                      = 64, &
    izeta_vrnce_rat               = 65, &
    iupsilon_precip_frac_rat      = 66, &
    ilambda0_stability_coef       = 67, &
    imult_coef                    = 68, &
    itaumin                       = 69, &
    itaumax                       = 70, &
    iLscale_mu_coef               = 71, &
    iLscale_pert_coef             = 72, &
    ialpha_corr                   = 73, &
    iSkw_denom_coef               = 74, &
    ic_K10                        = 75, &
    ic_K10h                       = 76, &
    ithlp2_rad_coef               = 77, &
    ithlp2_rad_cloud_frac_thresh  = 78, &
    iup2_sfc_coef                 = 79, &
    iSkw_max_mag                  = 80, &
    iC_invrs_tau_bkgnd            = 81, &
    iC_invrs_tau_sfc              = 82, &
    iC_invrs_tau_shear            = 83, &
    iC_invrs_tau_N2               = 84, &
    iC_invrs_tau_N2_wp2           = 85, &
    iC_invrs_tau_N2_xp2           = 86, &
    iC_invrs_tau_N2_wpxp          = 87, &
    iC_invrs_tau_N2_clear_wp3     = 88, &
    iC_invrs_tau_wpxp_Ri          = 89, &
    iC_invrs_tau_wpxp_N2_thresh   = 90, &
    ixp3_coef_base                = 91, &
    ixp3_coef_slope               = 92, &
    ialtitude_threshold           = 93, &
    irtp2_clip_coef               = 94, &
    iCx_min                       = 95, &
    iCx_max                       = 96, &
    iRichardson_num_min           = 97, &
    iRichardson_num_max           = 98

end module parameter_indices
!-----------------------------------------------------------------------
