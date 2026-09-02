"""Zero-based JAX mirror of ``src/CLUBB_core/parameter_indices.F90``.

Description:
  Since f90/95 lacks enumeration, we're stuck numbering each
  parameter by hand like this.

  Adding new parameters is relatively simple.  First, the
  parameter should be added in the common block of the parameters
  module so it can be used in other parts of the code. Each
  variable needs a unique number in this module, and nparams must
  be incremented for the new variable.  Next, the params_list
  variable in module parameters should have new variable added to
  it.  The subroutines pack_parameters and uppack_parameters will
  need to have the variable added to their list, but the order
  doesn't actually matter, since the i variables in here determine
  where in the params vector the number is placed.
  Finally, the namelists clubb_params_nl and initspread will need to
  have the parameter added to them.

Porting deviations:
- Fortran numbers the tunable parameters from 1 because it indexes
  ``clubb_params(ngrdcol, nparams)`` with 1-based array subscripts.  The JAX
  ``clubb_params`` array has shape ``(ngrdcol, nparams)`` with no padding
  column, so these constants are the corresponding zero-based Python column
  indices.

IMPORTANT:
  If you change the order of these parameters, you will need to
  change the order of params_list as well or the tuner will
  break!
"""

nparams = 102

iC1 = 0
iC1b = 1
iC1c = 2
iC2rt = 3
iC2thl = 4
iC2rtthl = 5
iC4 = 6
iC_uu_shr = 7
iC_uu_buoy = 8
iC6rt = 9
iC6rtb = 10
iC6rtc = 11
iC6thl = 12
iC6thlb = 13
iC6thlc = 14
iC7 = 15
iC7b = 16
iC7c = 17
iC8 = 18
iC8b = 19
iC10 = 20
iC11 = 21
iC11b = 22
iC11c = 23
iC12 = 24
iC13 = 25
iC14 = 26
iC_wp2_pr_dfsn = 27
iC_wp3_pr_tp = 28
iC_wp3_pr_turb = 29
iC_wp3_pr_dfsn = 30
iC_wp2_splat = 31

iC6rt_Lscale0 = 32
iC6thl_Lscale0 = 33
iC7_Lscale0 = 34
iwpxp_L_thresh = 35

ic_K = 36
ic_K1 = 37
inu1 = 38
ic_K2 = 39
inu2 = 40
ic_K6 = 41
inu6 = 42
ic_K8 = 43
inu8 = 44
ic_K9 = 45
inu9 = 46
inu10 = 47
ic_K_hm = 48
ic_K_hmb = 49
iK_hm_min_coef = 50
inu_hm = 51

islope_coef_spread_DG_means_w = 52
ipdf_component_stdev_factor_w = 53
icoef_spread_DG_means_rt = 54
icoef_spread_DG_means_thl = 55
igamma_coef = 56
igamma_coefb = 57
igamma_coefc = 58
imu = 59
ibeta = 60
ilmin_coef = 61
iomicron = 62
izeta_vrnce_rat = 63
iupsilon_precip_frac_rat = 64
ilambda0_stability_coef = 65
imult_coef = 66
itaumin = 67
itaumax = 68
iLscale_mu_coef = 69
iLscale_pert_coef = 70
ialpha_corr = 71
iSkw_denom_coef = 72
ic_K10 = 73
ic_K10h = 74
ithlp2_rad_coef = 75
ithlp2_rad_cloud_frac_thresh = 76
iup2_sfc_coef = 77
iSkw_max_mag = 78
iC_invrs_tau_bkgnd = 79
iC_invrs_tau_sfc = 80
iC_invrs_tau_shear = 81
iC_invrs_tau_N2 = 82
iC_invrs_tau_N2_wp2 = 83
iC_invrs_tau_N2_xp2 = 84
iC_invrs_tau_N2_wpxp = 85
iC_invrs_tau_N2_clear_wp3 = 86
iC_invrs_tau_wpxp_Ri = 87
iC_invrs_tau_wpxp_N2_thresh = 88
ixp3_coef_base = 89
ixp3_coef_slope = 90
ialtitude_threshold = 91
irtp2_clip_coef = 92
iCx_min = 93
iCx_max = 94
iRichardson_num_min = 95
iRichardson_num_max = 96
ia3_coef_min = 97
ia_const = 98
ibv_efold = 99
iwpxp_Ri_exp = 100
iz_displace = 101
