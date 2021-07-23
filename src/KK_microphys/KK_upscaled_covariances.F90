! $Id$
!===============================================================================
module KK_upscaled_covariances

  implicit none

  private  ! Default scope

  public :: KK_upscaled_covar_driver

  private :: covar_x_KK_evap,         &
             covar_rt_KK_evap,        &
             covar_thl_KK_evap,       &
             covar_x_KK_auto,         &
             covar_rt_KK_auto,        &
             covar_thl_KK_auto,       &
             covar_x_KK_accr,         &
             covar_rt_KK_accr,        &
             covar_thl_KK_accr,       &
             quadrivar_NNLL_covar_eq, &
             trivar_NNL_covar_eq

  contains

  !=============================================================================
  subroutine KK_upscaled_covar_driver( w_mean, rtm, thlm, &
                                       exner, rrainm, Nrm, &
                                       mu_w_1, mu_w_2, mu_chi_1, mu_chi_2, &
                                       mu_eta_1, mu_eta_2, mu_rr_1, mu_rr_2, &
                                       mu_Nr_1, mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &
                                       mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                       mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &
                                       sigma_w_1, sigma_w_2, sigma_chi_1, &
                                       sigma_chi_2, sigma_eta_1, sigma_eta_2, &
                                       sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                       sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &
                                       sigma_rr_1_n, sigma_rr_2_n, &
                                       sigma_Nr_1_n, sigma_Nr_2_n, &
                                       sigma_Ncn_1_n, sigma_Ncn_2_n, &
                                       corr_w_chi_1, corr_w_chi_2, &
                                       corr_w_rr_1_n, corr_w_rr_2_n, &
                                       corr_w_Nr_1_n, corr_w_Nr_2_n, &
                                       corr_w_Ncn_1_n, corr_w_Ncn_2_n, &
                                       corr_chi_eta_1, corr_chi_eta_2, &
                                       corr_chi_rr_1_n, corr_chi_rr_2_n, &
                                       corr_chi_Nr_1_n, corr_chi_Nr_2_n, &
                                       corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, &
                                       corr_eta_rr_1_n, corr_eta_rr_2_n, &
                                       corr_eta_Nr_1_n, corr_eta_Nr_2_n, &
                                       corr_eta_Ncn_1_n, corr_eta_Ncn_2_n, &
                                       corr_rr_Nr_1_n, corr_rr_Nr_2_n, &
                                       mixt_frac, precip_frac_1, &
                                       precip_frac_2, &
                                       KK_evap_coef, KK_auto_coef, &
                                       KK_accr_coef, KK_evap_tndcy, &
                                       KK_auto_tndcy, KK_accr_tndcy, &
                                       mu_rt_1, mu_rt_2, &
                                       mu_thl_1,mu_thl_2, &
                                       crt1, crt2, &
                                       cthl1, cthl2, &
                                       level, l_stats_samp, &
                                       stats_zt, &
                                       wprtp_mc_src_tndcy, &
                                       wpthlp_mc_src_tndcy, &
                                       rtp2_mc_src_tndcy, &
                                       thlp2_mc_src_tndcy, &
                                       rtpthlp_mc_src_tndcy )
    ! Description:

    ! References:
    ! Griffin, B. M., 2016:  Improving the Subgrid-Scale Representation of
    !    Hydrometeors and Microphysical Feedback Effects Using a Multivariate
    !    PDF.  Doctoral dissertation, University of Wisconsin -- Milwaukee,
    !    Milwaukee, WI, Paper 1144, 165 pp., URL
    !    http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Griffin, B. M. and V. E. Larson, 2016:  Parameterizing microphysical
    !    effects on variances and covariances of moisture and heat content using
    !    a multivariate probability density function: a study with CLUBB (tag
    !    MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
    !    doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        two,    & ! Constant(s)
        zero,   &
        Lv,     &
        Cp,     &
        w_tol,  &
        rr_tol, & 
        Nr_tol 

    use constants_clubb, only:  &
        eta_tol  ! Constant

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    use stats_type_utilities, only: & 
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only: &
        iw_KK_evap_covar_zt,   & ! Variable(s)
        irt_KK_evap_covar_zt,  &
        ithl_KK_evap_covar_zt, &
        iw_KK_auto_covar_zt,   &
        irt_KK_auto_covar_zt,  &
        ithl_KK_auto_covar_zt, &
        iw_KK_accr_covar_zt,   &
        irt_KK_accr_covar_zt,  &
        ithl_KK_accr_covar_zt

    use stats_type, only: stats ! Type

    implicit none

    type(stats), target, intent(inout) :: &
      stats_zt

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      w_mean, & ! Mean vertical velocity, w (overall)               [m/s]
      rtm,    & ! Mean total water mixing ratio, rt (overall)       [kg/kg]
      thlm,   & ! Mean liquid water potential temp., thl (overall)  [K]
      exner,  & ! Exner function                                    [-]
      rrainm, & ! Mean rain water mixing ratio (overall)            [kg/kg]
      Nrm       ! Mean rain drop concentration (overall)            [num/kg]

    real( kind = core_rknd ), intent(in) :: &
      mu_w_1,        & ! Mean of w (1st PDF component)                     [m/s]
      mu_w_2,        & ! Mean of w (2nd PDF component)                     [m/s]
      mu_chi_1,      & ! Mean of chi (old s) (1st PDF component)         [kg/kg]
      mu_chi_2,      & ! Mean of chi (old s) (2nd PDF component)         [kg/kg]
      mu_eta_1,      & ! Mean of eta (old t) (1st PDF component)         [kg/kg]
      mu_eta_2,      & ! Mean of eta (old t) (2nd PDF component)         [kg/kg]
      mu_rr_1,       & ! Mean of rr (1st PDF component) in-precip (ip)   [kg/kg]
      mu_rr_2,       & ! Mean of rr (2nd PDF component) ip               [kg/kg]
      mu_Nr_1,       & ! Mean of Nr (1st PDF component) ip              [num/kg]
      mu_Nr_2,       & ! Mean of Nr (2nd PDF component) ip              [num/kg]
      mu_Ncn_1,      & ! Mean of Ncn (1st PDF component)                [num/kg]
      mu_Ncn_2,      & ! Mean of Ncn (2nd PDF component)                [num/kg]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) ip        [ln(kg/kg)]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip        [ln(kg/kg)]
      mu_Nr_1_n,     & ! Mean of ln Nr (1st PDF component) ip       [ln(num/kg)]
      mu_Nr_2_n,     & ! Mean of ln Nr (2nd PDF component) ip       [ln(num/kg)]
      mu_Ncn_1_n,    & ! Mean of ln Ncn (1st PDF component)         [ln(num/kg)]
      mu_Ncn_2_n,    & ! Mean of ln Ncn (2nd PDF component)         [ln(num/kg)]
      sigma_w_1,     & ! Standard deviation of w (1st PDF component)       [m/s]
      sigma_w_2,     & ! Standard deviation of w (2nd PDF component)       [m/s]
      sigma_chi_1,   & ! Standard deviation of chi (1st PDF component)   [kg/kg]
      sigma_chi_2,   & ! Standard deviation of chi (2nd PDF component)   [kg/kg]
      sigma_eta_1,   & ! Standard deviation of eta (1st PDF component)   [kg/kg]
      sigma_eta_2,   & ! Standard deviation of eta (2nd PDF component)   [kg/kg]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_Nr_1,    & ! Standard deviation of Nr (1st PDF comp.) ip    [num/kg]
      sigma_Nr_2,    & ! Standard deviation of Nr (2nd PDF comp.) ip    [num/kg]
      sigma_Ncn_1,   & ! Standard deviation of Ncn (1st PDF component)  [num/kg]
      sigma_Ncn_2,   & ! Standard deviation of Ncn (2nd PDF component)  [num/kg]
      sigma_rr_1_n,  & ! Standard deviation of ln rr (1st PDF component) ip  [-]
      sigma_rr_2_n,  & ! Standard deviation of ln rr (2nd PDF component) ip  [-]
      sigma_Nr_1_n,  & ! Standard deviation of ln Nr (1st PDF component) ip  [-]
      sigma_Nr_2_n,  & ! Standard deviation of ln Nr (2nd PDF component) ip  [-]
      sigma_Ncn_1_n, & ! Standard deviation of ln Ncn (1st PDF component)    [-]
      sigma_Ncn_2_n    ! Standard deviation of ln Ncn (2nd PDF component)    [-]

    real( kind = core_rknd ), intent(in) :: &
      corr_w_chi_1,     & ! Correlation of w and chi (1st PDF component)    [-]
      corr_w_chi_2,     & ! Correlation of w and chi (2nd PDF component)    [-]
      corr_w_rr_1_n,    & ! Correlation of w and ln rr (1st PDF comp.) ip   [-]
      corr_w_rr_2_n,    & ! Correlation of w and ln rr (2nd PDF comp.) ip   [-]
      corr_w_Nr_1_n,    & ! Correlation of w and ln Nr (1st PDF comp.) ip   [-]
      corr_w_Nr_2_n,    & ! Correlation of w and ln Nr (2nd PDF comp.) ip   [-]
      corr_w_Ncn_1_n,   & ! Correlation of w and ln Ncn (1st PDF comp.)     [-]
      corr_w_Ncn_2_n,   & ! Correlation of w and ln Ncn (2nd PDF comp.)     [-]
      corr_chi_eta_1,   & ! Correlation of chi and eta (1st PDF component)  [-]
      corr_chi_eta_2,   & ! Correlation of chi and eta (2nd PDF component)  [-]
      corr_chi_rr_1_n,  & ! Correlation of chi and ln rr (1st PDF comp.) ip [-]
      corr_chi_rr_2_n,  & ! Correlation of chi and ln rr (2nd PDF comp.) ip [-]
      corr_chi_Nr_1_n,  & ! Correlation of chi and ln Nr (1st PDF comp.) ip [-]
      corr_chi_Nr_2_n,  & ! Correlation of chi and ln Nr (2nd PDF comp.) ip [-]
      corr_chi_Ncn_1_n, & ! Correlation of chi and ln Ncn (1st PDF comp.)   [-]
      corr_chi_Ncn_2_n, & ! Correlation of chi and ln Ncn (2nd PDF comp.)   [-]
      corr_eta_rr_1_n,  & ! Correlation of eta and ln rr (1st PDF comp.) ip [-]
      corr_eta_rr_2_n,  & ! Correlation of eta and ln rr (2nd PDF comp.) ip [-]
      corr_eta_Nr_1_n,  & ! Correlation of eta and ln Nr (1st PDF comp.) ip [-]
      corr_eta_Nr_2_n,  & ! Correlation of eta and ln Nr (2nd PDF comp.) ip [-]
      corr_eta_Ncn_1_n, & ! Correlation of eta and ln Ncn (1st PDF comp.)   [-]
      corr_eta_Ncn_2_n, & ! Correlation of eta and ln Ncn (2nd PDF comp.)   [-]
      corr_rr_Nr_1_n,   & ! Correlation of ln rr & ln Nr (1st PDF comp.) ip [-]
      corr_rr_Nr_2_n,   & ! Correlation of ln rr & ln Nr (2nd PDF comp.) ip [-]
      mixt_frac,        & ! Mixture fraction                                [-]
      precip_frac_1,    & ! Precipitation fraction (1st PDF component)      [-]
      precip_frac_2       ! Precipitation fraction (2nd PDF component)      [-]

    real( kind = core_rknd ), intent(in) :: &
      KK_evap_coef, & ! KK evap. coef. [(kg/kg)^(1-s_ex-rr_ex)(num/kg)^-Nr_ex/s]
      KK_auto_coef, & ! KK auto. coef.   [(kg/kg)^(1-s_ex) (num/kg)^-Ncn_ex / s]
      KK_accr_coef    ! KK accr. coef.                [(kg/kg)^(1-s_ex-rr_ex)/s]

    real( kind = core_rknd ), intent(in) :: &
      KK_evap_tndcy, & ! KK evaporation tendency            [(kg/kg)/s]
      KK_auto_tndcy, & ! KK autoconversion tendency         [(kg/kg)/s]
      KK_accr_tndcy    ! KK accretion tendency              [(kg/kg)/s]

    real( kind = core_rknd ), intent(in) :: &
      mu_rt_1,  & ! Mean of rt (PDF component 1)            [kg/kg]
      mu_rt_2,  & ! Mean of rt (PDF component 2)            [kg/kg]
      mu_thl_1, & ! Mean of thl (PDF component 1)           [K]
      mu_thl_2, & ! Mean of thl (PDF component 2)           [K]
      crt1,     & ! Coefficient c_rt (1st PDF component)    [-]
      crt2,     & ! Coefficient c_rt (2nd PDF component)    [-]
      cthl1,    & ! Coefficient c_thl (1st PDF component)   [(kg/kg)/K]
      cthl2       ! Coefficient c_thl (2nd PDF component)   [(kg/kg)/K]

    integer, intent(in) :: &
      level         ! Vertical level index                  [-]

    logical, intent(in) :: &
      l_stats_samp    ! Flag to record statistical output.

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      wprtp_mc_src_tndcy,   & ! Microphysics tendency for w'rt'  [m*(kg/kg)/s^2]
      wpthlp_mc_src_tndcy,  & ! Microphysics tendency for w'thl' [m*K/s^2]
      rtp2_mc_src_tndcy,    & ! Microphysics tendency for rt'^2  [(kg/kg)^2/s]
      thlp2_mc_src_tndcy,   & ! Microphysics tendency for thl'^2 [K^2/s]
      rtpthlp_mc_src_tndcy    ! Microphysics tend. for rt'thl'   [K*(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      w_KK_evap_covar,   & ! Covar. of w and KK evap. tend.    [m*(kg/kg)/s^2]
      rt_KK_evap_covar,  & ! Covar. of rt and KK evap. tend.   [(kg/kg)^2/s]
      thl_KK_evap_covar, & ! Covar. of thl and KK evap. tend.  [K*(kg/kg)/s]
      w_KK_auto_covar,   & ! Covar. of w and KK auto. tend.    [m*(kg/kg)/s^2]
      rt_KK_auto_covar,  & ! Covar. of rt and KK auto. tend.   [(kg/kg)^2/s]
      thl_KK_auto_covar, & ! Covar. of thl and KK auto. tend.  [K*(kg/kg)/s]
      w_KK_accr_covar,   & ! Covar. of w and KK accr. tend.    [m*(kg/kg)/s^2]
      rt_KK_accr_covar,  & ! Covar. of rt and KK accr. tend.   [(kg/kg)^2/s]
      thl_KK_accr_covar    ! Covar. of thl and KK accr. tend.  [K*(kg/kg)/s]

    
    ! Calculate the covariance of vertical velocity and KK evaporation tendency.
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then

       w_KK_evap_covar &
       = covar_x_KK_evap( mu_w_1, mu_w_2, mu_chi_1, mu_chi_2, mu_rr_1, &
                          mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, &
                          mu_Nr_1_n, mu_Nr_2_n, sigma_w_1, sigma_w_2, &
                          sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2, &
                          sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                          sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                          corr_w_chi_1, corr_w_chi_2, corr_w_rr_1_n, &
                          corr_w_rr_2_n, corr_w_Nr_1_n, corr_w_Nr_2_n, &
                          corr_chi_rr_1_n, corr_chi_rr_2_n, corr_chi_Nr_1_n, &
                          corr_chi_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, &
                          w_mean, KK_evap_tndcy, KK_evap_coef, w_tol, &
                          mixt_frac, precip_frac_1, precip_frac_2 )

    else  ! r_r or N_r = 0.

       w_KK_evap_covar = zero

    endif

    ! Calculate the covariance of total water mixing ratio and KK evaporation
    ! tendency.
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then

       rt_KK_evap_covar  &
       = covar_rt_KK_evap( mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_rr_1, &
                           mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, &
                           mu_Nr_1_n, mu_Nr_2_n, sigma_eta_1, sigma_eta_2, &
                           sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2, &
                           sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                           sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                           corr_chi_eta_1, corr_chi_eta_2, corr_eta_rr_1_n, &
                           corr_eta_rr_2_n, corr_eta_Nr_1_n, & 
                           corr_eta_Nr_2_n, corr_chi_rr_1_n, &
                           corr_chi_rr_2_n, corr_chi_Nr_1_n, &
                           corr_chi_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, &
                           mixt_frac, precip_frac_1, precip_frac_2, &
                           rtm, mu_rt_1, mu_rt_2, KK_evap_tndcy, &
                           KK_evap_coef, eta_tol, crt1, crt2 )

    else  ! r_r or N_r = 0.

       rt_KK_evap_covar = zero

    endif

    ! Calculate the covariance of liquid water potential temperature and
    ! KK evaporation tendency.
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then

       thl_KK_evap_covar &
       = covar_thl_KK_evap( mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_rr_1, &
                            mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, &
                            mu_Nr_1_n, mu_Nr_2_n, sigma_eta_1, sigma_eta_2, &
                            sigma_chi_1, sigma_chi_2, sigma_rr_1, &
                            sigma_rr_2, sigma_Nr_1, sigma_Nr_2, &
                            sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                            sigma_Nr_2_n, corr_chi_eta_1, corr_chi_eta_2, &
                            corr_eta_rr_1_n, corr_eta_rr_2_n, &
                            corr_eta_Nr_1_n, corr_eta_Nr_2_n, &
                            corr_chi_rr_1_n, corr_chi_rr_2_n, &
                            corr_chi_Nr_1_n, corr_chi_Nr_2_n, &
                            corr_rr_Nr_1_n, corr_rr_Nr_2_n, &
                            mixt_frac, precip_frac_1, precip_frac_2, &
                            thlm, mu_thl_1, mu_thl_2, KK_evap_tndcy, &
                            KK_evap_coef, eta_tol, cthl1, cthl2 )

    else  ! r_r or N_r = 0.

       thl_KK_evap_covar = zero

    endif

    ! Calculate the covariance of vertical velocity and KK autoconversion
    ! tendency.
    w_KK_auto_covar &
    = covar_x_KK_auto( mu_w_1, mu_w_2, mu_chi_1, mu_chi_2, mu_Ncn_1, &
                       mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n, sigma_w_1, &
                       sigma_w_2, sigma_chi_1, sigma_chi_2, sigma_Ncn_1, &
                       sigma_Ncn_2, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                       corr_w_chi_1, corr_w_chi_2, corr_w_Ncn_1_n, &
                       corr_w_Ncn_2_n, corr_chi_Ncn_1_n, &
                       corr_chi_Ncn_2_n, w_mean, KK_auto_tndcy, &
                       KK_auto_coef, w_tol, mixt_frac )

    ! Calculate the covariance of total water mixing ratio and KK autoconversion
    ! tendency.
    rt_KK_auto_covar &
    = covar_rt_KK_auto( mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_Ncn_1, &
                        mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n, sigma_eta_1, &
                        sigma_eta_2, sigma_chi_1, sigma_chi_2, &
                        sigma_Ncn_1, sigma_Ncn_2, sigma_Ncn_1_n, &
                        sigma_Ncn_2_n, corr_chi_eta_1, corr_chi_eta_2, &
                        corr_eta_Ncn_1_n, corr_eta_Ncn_2_n, &
                        corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, &
                        rtm, mu_rt_1, mu_rt_2, KK_auto_tndcy, &
                        KK_auto_coef, eta_tol, crt1, crt2, mixt_frac )

    ! Calculate the covariance of liquid water potential temperature and
    ! KK autoconversion tendency.
    thl_KK_auto_covar &
    = covar_thl_KK_auto( mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, &
                         mu_Ncn_1, mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n, &
                         sigma_eta_1, sigma_eta_2, sigma_chi_1, &
                         sigma_chi_2, sigma_Ncn_1, sigma_Ncn_2, &
                         sigma_Ncn_1_n, sigma_Ncn_2_n, corr_chi_eta_1, &
                         corr_chi_eta_2, corr_eta_Ncn_1_n, &
                         corr_eta_Ncn_2_n, corr_chi_Ncn_1_n, &
                         corr_chi_Ncn_2_n, thlm, mu_thl_1, mu_thl_2, &
                         KK_auto_tndcy, KK_auto_coef, eta_tol, &
                         cthl1, cthl2, mixt_frac )

    ! Calculate the covariance of vertical velocity and KK accretion tendency.
    if ( rrainm > rr_tol ) then

       w_KK_accr_covar &
       = covar_x_KK_accr( mu_w_1, mu_w_2, mu_chi_1, mu_chi_2, mu_rr_1, &
                          mu_rr_2, mu_rr_1_n, mu_rr_2_n, sigma_w_1, &
                          sigma_w_2, sigma_chi_1, sigma_chi_2, sigma_rr_1, &
                          sigma_rr_2, sigma_rr_1_n, sigma_rr_2_n, &
                          corr_w_chi_1, corr_w_chi_2, corr_w_rr_1_n, &
                          corr_w_rr_2_n, corr_chi_rr_1_n, corr_chi_rr_2_n, &
                          w_mean, KK_accr_tndcy, KK_accr_coef, w_tol, &
                          mixt_frac, precip_frac_1, precip_frac_2 )

    else  ! r_r = 0.

       w_KK_accr_covar = zero

    endif

    ! Calculate the covariance of total water mixing ratio and KK accretion
    ! tendency.
    if ( rrainm > rr_tol ) then

       rt_KK_accr_covar &
       = covar_rt_KK_accr( mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_rr_1, &
                           mu_rr_2, mu_rr_1_n, mu_rr_2_n, sigma_eta_1, &
                           sigma_eta_2, sigma_chi_1, sigma_chi_2, &
                           sigma_rr_1, sigma_rr_2, sigma_rr_1_n, &
                           sigma_rr_2_n, corr_chi_eta_1, corr_chi_eta_2, &
                           corr_eta_rr_1_n, corr_eta_rr_2_n, &
                           corr_chi_rr_1_n, corr_chi_rr_2_n, &
                           rtm, mu_rt_1, mu_rt_2, KK_accr_tndcy, &
                           KK_accr_coef, eta_tol, crt1, crt2, mixt_frac, &
                           precip_frac_1, precip_frac_2 )

    else  ! r_r = 0.

       rt_KK_accr_covar = zero

    endif

    ! Calculate the covariance of liquid water potential temperature and
    ! KK accretion tendency.
    if ( rrainm > rr_tol ) then

       thl_KK_accr_covar &
       = covar_thl_KK_accr( mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_rr_1, &
                            mu_rr_2, mu_rr_1_n, mu_rr_2_n, sigma_eta_1, &
                            sigma_eta_2, sigma_chi_1, sigma_chi_2, &
                            sigma_rr_1, sigma_rr_2, sigma_rr_1_n, &
                            sigma_rr_2_n, corr_chi_eta_1, corr_chi_eta_2, &
                            corr_eta_rr_1_n, corr_eta_rr_2_n, &
                            corr_chi_rr_1_n, corr_chi_rr_2_n, &
                            thlm, mu_thl_1, mu_thl_2, KK_accr_tndcy, &
                            KK_accr_coef, eta_tol, cthl1, cthl2, mixt_frac, &
                            precip_frac_1, precip_frac_2 )

    else  ! r_r = 0.

       thl_KK_accr_covar = zero

    endif


    ! Statistics
    if ( l_stats_samp ) then
       ! All of these covariance variables are being calculated on thermodynamic
       ! grid levels (all inputs are on thermodynamic grid levels, so the output
       ! is also on thermodynamic grid levels).  These covariances will be
       ! combined in various ways to produce the microphysics tendency terms for
       ! various model predictive variances and covariances.  These source
       ! tendency terms will be interpolated to momentum grid levels.

       ! Covariance of w and KK evaporation tendency.
       if ( iw_KK_evap_covar_zt > 0 ) then
          call stat_update_var_pt( iw_KK_evap_covar_zt, level, &
                                   w_KK_evap_covar, stats_zt )
       endif

       ! Covariance of r_t and KK evaporation tendency.
       if ( irt_KK_evap_covar_zt > 0 ) then
          call stat_update_var_pt( irt_KK_evap_covar_zt, level, &
                                   rt_KK_evap_covar, stats_zt )
       endif

       ! Covariance of theta_l and KK evaporation tendency.
       if ( ithl_KK_evap_covar_zt > 0 ) then
          call stat_update_var_pt( ithl_KK_evap_covar_zt, level, &
                                   thl_KK_evap_covar, stats_zt )
       endif

       ! Covariance of w and KK autoconversion tendency.
       if ( iw_KK_auto_covar_zt > 0 ) then
          call stat_update_var_pt( iw_KK_auto_covar_zt, level, &
                                   w_KK_auto_covar, stats_zt )
       endif

       ! Covariance of r_t and KK autoconversion tendency.
       if ( irt_KK_auto_covar_zt > 0 ) then
          call stat_update_var_pt( irt_KK_auto_covar_zt, level, &
                                   rt_KK_auto_covar, stats_zt )
       endif

       ! Covariance of theta_l and KK autoconversion tendency.
       if ( ithl_KK_auto_covar_zt > 0 ) then
          call stat_update_var_pt( ithl_KK_auto_covar_zt, level, &
                                   thl_KK_auto_covar, stats_zt )
       endif

       ! Covariance of w and KK accretion tendency.
       if ( iw_KK_auto_covar_zt > 0 ) then
          call stat_update_var_pt( iw_KK_accr_covar_zt, level, &
                                   w_KK_accr_covar, stats_zt )
       endif

       ! Covariance of r_t and KK accretion tendency.
       if ( irt_KK_auto_covar_zt > 0 ) then
          call stat_update_var_pt( irt_KK_accr_covar_zt, level, &
                                   rt_KK_accr_covar, stats_zt )
       endif

       ! Covariance of theta_l and KK accretion tendency.
       if ( ithl_KK_auto_covar_zt > 0 ) then
          call stat_update_var_pt( ithl_KK_accr_covar_zt, level, &
                                   thl_KK_accr_covar, stats_zt )
       endif

    endif ! l_stats_samp


    ! Calculate the microphysics tendency for <w'r_t'>.
    wprtp_mc_src_tndcy  &
    = - ( w_KK_auto_covar + w_KK_accr_covar + w_KK_evap_covar )

    ! Calculate the microphysics tendency for <w'th_l'>.
    wpthlp_mc_src_tndcy  &
    = ( Lv / ( Cp * exner ) )  &
      * ( w_KK_auto_covar + w_KK_accr_covar + w_KK_evap_covar )

    ! Calculate the microphysics tendency for <r_t'^2>.
    rtp2_mc_src_tndcy  &
    = - two * ( rt_KK_auto_covar + rt_KK_accr_covar + rt_KK_evap_covar )

    ! Calculate the microphysics tendency for <th_l'^2>.
    thlp2_mc_src_tndcy  &
    = two * ( Lv / ( Cp * exner ) )  &
          * ( thl_KK_auto_covar + thl_KK_accr_covar + thl_KK_evap_covar )

    ! Calculate the microphysics tendency for <r_t'th_l'>.
    rtpthlp_mc_src_tndcy  &
    = ( Lv / ( Cp * exner ) )  &
      * ( rt_KK_auto_covar + rt_KK_accr_covar + rt_KK_evap_covar )  &
      - ( thl_KK_auto_covar + thl_KK_accr_covar + thl_KK_evap_covar )


    return

  end subroutine KK_upscaled_covar_driver

  !=============================================================================
  function covar_x_KK_evap( mu_x_1, mu_x_2, mu_chi_1, mu_chi_2, mu_rr_1, &
                            mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, &
                            mu_Nr_1_n, mu_Nr_2_n, sigma_x_1, sigma_x_2, &
                            sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2, &
                            sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                            sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                            corr_x_chi_1, corr_x_chi_2, corr_x_rr_1_n, &
                            corr_x_rr_2_n, corr_x_Nr_1_n, corr_x_Nr_2_n, &
                            corr_chi_rr_1_n, corr_chi_rr_2_n, corr_chi_Nr_1_n, &
                            corr_chi_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, &
                            x_mean, KK_evap_tndcy, KK_evap_coef, x_tol, &
                            mixt_frac, precip_frac_1, precip_frac_2 )

    ! Description:
    ! This function calculates the covariance of x and KK evaporation
    ! tendency, which can be written as < x'((dr_r/dt)_KKevap)' >, or more
    ! simply as < x'KK_evap' >.

    ! References:
    ! Eq. (E21) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (A28) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
    ! microphysical effects on variances and covariances of moisture and heat
    ! content using a multivariate probability density function: a study with
    ! CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one  ! Constant(s)

    use parameters_KK, only: &
        KK_evap_Supersat_exp, & ! Variable(s)
        KK_evap_rr_exp,       &
        KK_evap_Nr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x_1,       & ! Mean of x (1st PDF component)               [units vary]
      mu_x_2,       & ! Mean of x (2nd PDF component)               [units vary]
      mu_chi_1,     & ! Mean of chi (old s) (1st PDF component)          [kg/kg]
      mu_chi_2,     & ! Mean of chi (old s) (2nd PDF component)          [kg/kg]
      mu_rr_1,      & ! Mean of rr (1st PDF component) in-precip (ip)    [kg/kg]
      mu_rr_2,      & ! Mean of rr (2nd PDF component) ip                [kg/kg]
      mu_Nr_1,      & ! Mean of Nr (1st PDF component) ip               [num/kg]
      mu_Nr_2,      & ! Mean of Nr (2nd PDF component) ip               [num/kg]
      mu_rr_1_n,    & ! Mean of ln rr (1st PDF component) ip         [ln(kg/kg)]
      mu_rr_2_n,    & ! Mean of ln rr (2nd PDF component) ip         [ln(kg/kg)]
      mu_Nr_1_n,    & ! Mean of ln Nr (1st PDF component) ip        [ln(num/kg)]
      mu_Nr_2_n,    & ! Mean of ln Nr (2nd PDF component) ip        [ln(num/kg)]
      sigma_x_1,    & ! Standard deviation of x (1st PDF component)   [un. vary]
      sigma_x_2,    & ! Standard deviation of x (2nd PDF component)   [un. vary]
      sigma_chi_1,  & ! Standard deviation of chi (1st PDF component)    [kg/kg]
      sigma_chi_2,  & ! Standard deviation of chi (2nd PDF component)    [kg/kg]
      sigma_rr_1,   & ! Standard deviation of rr (1st PDF component) ip  [kg/kg]
      sigma_rr_2,   & ! Standard deviation of rr (2nd PDF component) ip  [kg/kg]
      sigma_Nr_1,   & ! Standard deviation of Nr (1st PDF component) ip [num/kg]
      sigma_Nr_2,   & ! Standard deviation of Nr (2nd PDF component) ip [num/kg]
      sigma_rr_1_n, & ! Standard deviation of ln rr (1st PDF component) ip   [-]
      sigma_rr_2_n, & ! Standard deviation of ln rr (2nd PDF component) ip   [-]
      sigma_Nr_1_n, & ! Standard deviation of ln Nr (1st PDF component) ip   [-]
      sigma_Nr_2_n    ! Standard deviation of ln Nr (2nd PDF component) ip   [-]

    real( kind = core_rknd ), intent(in) :: &
      corr_x_chi_1,    & ! Correlation of x and chi (1st PDF component)      [-]
      corr_x_chi_2,    & ! Correlation of x and chi (2nd PDF component)      [-]
      corr_x_rr_1_n,   & ! Correlation of x and ln rr (1st PDF comp.) ip     [-]
      corr_x_rr_2_n,   & ! Correlation of x and ln rr (2nd PDF comp.) ip     [-]
      corr_x_Nr_1_n,   & ! Correlation of x and ln Nr (1st PDF comp.) ip     [-]
      corr_x_Nr_2_n,   & ! Correlation of x and ln Nr (2nd PDF comp.) ip     [-]
      corr_chi_rr_1_n, & ! Correlation of chi and ln rr (1st PDF comp.) ip   [-]
      corr_chi_rr_2_n, & ! Correlation of chi and ln rr (2nd PDF comp.) ip   [-]
      corr_chi_Nr_1_n, & ! Correlation of chi and ln Nr (1st PDF comp.) ip   [-]
      corr_chi_Nr_2_n, & ! Correlation of chi and ln Nr (2nd PDF comp.) ip   [-]
      corr_rr_Nr_1_n,  & ! Correlation of ln rr & ln Nr (1st PDF comp.) ip   [-]
      corr_rr_Nr_2_n,  & ! Correlation of ln rr & ln Nr (2nd PDF comp.) ip   [-]
      x_mean,          & ! Mean of x (overall)                      [units vary]
      KK_evap_tndcy,   & ! KK evaporation tendency                   [(kg/kg)/s]
      KK_evap_coef,    & ! KK evap. coef.[(kg/kg)^(1-alpha-beta)(#/kg)^-gamma/s]
      x_tol,           & ! Tolerance value of x                     [units vary]
      mixt_frac,       & ! Mixture fraction                                  [-]
      precip_frac_1,   & ! Precipitation fraction (1st PDF component)        [-]
      precip_frac_2      ! Precipitation fraction (2nd PDF component)        [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_x_KK_evap  ! Covariance of x and KK evap. tendency   [u.v.(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on chi                                         [-]
      beta_exp,  & ! Exponent on r_r                                         [-]
      gamma_exp    ! Exponent on N_r                                         [-]


    ! Values of the KK exponents.
    alpha_exp = KK_evap_Supersat_exp
    beta_exp  = KK_evap_rr_exp
    gamma_exp = KK_evap_Nr_exp

    ! Calculate the covariance of x and KK evaporation tendency.
    covar_x_KK_evap  &
    = mixt_frac &
      * ( KK_evap_coef * precip_frac_1 &
          * quadrivar_NNLL_covar_eq( mu_x_1, mu_chi_1, mu_rr_1, mu_Nr_1, &
                                     mu_rr_1_n, mu_Nr_1_n, sigma_x_1, &
                                     sigma_chi_1, sigma_rr_1, sigma_Nr_1, &
                                     sigma_rr_1_n, sigma_Nr_1_n, &
                                     corr_x_chi_1, corr_x_rr_1_n, &
                                     corr_x_Nr_1_n, corr_chi_rr_1_n, &
                                     corr_chi_Nr_1_n, corr_rr_Nr_1_n, &
                                     x_mean, KK_evap_tndcy, &
                                     KK_evap_coef, x_tol, &
                                     alpha_exp, beta_exp, gamma_exp ) &
          - ( one - precip_frac_1 ) * ( mu_x_1 - x_mean ) * KK_evap_tndcy &
        ) &
      + ( one - mixt_frac ) &
        * ( KK_evap_coef * precip_frac_2 &
            * quadrivar_NNLL_covar_eq( mu_x_2, mu_chi_2, mu_rr_2, mu_Nr_2, &
                                       mu_rr_2_n, mu_Nr_2_n, sigma_x_2, &
                                       sigma_chi_2, sigma_rr_2, sigma_Nr_2, &
                                       sigma_rr_2_n, sigma_Nr_2_n, &
                                       corr_x_chi_2, corr_x_rr_2_n, &
                                       corr_x_Nr_2_n, corr_chi_rr_2_n, &
                                       corr_chi_Nr_2_n, corr_rr_Nr_2_n, &
                                       x_mean, KK_evap_tndcy, &
                                       KK_evap_coef, x_tol, &
                                       alpha_exp, beta_exp, gamma_exp ) &
            - ( one - precip_frac_2 ) * ( mu_x_2 - x_mean ) * KK_evap_tndcy &
          )


    return

  end function covar_x_KK_evap

  !=============================================================================
  function covar_rt_KK_evap( mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_rr_1, &
                             mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, &
                             mu_Nr_1_n, mu_Nr_2_n, sigma_eta_1, sigma_eta_2, &
                             sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2, &
                             sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                             sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                             corr_eta_chi_1, corr_eta_chi_2, corr_eta_rr_1_n, &
                             corr_eta_rr_2_n, corr_eta_Nr_1_n, & 
                             corr_eta_Nr_2_n, corr_chi_rr_1_n, &
                             corr_chi_rr_2_n, corr_chi_Nr_1_n, &
                             corr_chi_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, &
                             mixt_frac, precip_frac_1, precip_frac_2, &
                             rtm, mu_rt_1, mu_rt_2, KK_evap_tndcy, &
                             KK_evap_coef, eta_tol, crt1, crt2 )

    ! Description:

    ! References:
    ! Eq. (E24) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (A31) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
    ! microphysical effects on variances and covariances of moisture and heat
    ! content using a multivariate probability density function: a study with
    ! CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two, & ! Constant(s)
        one

    use KK_upscaled_means, only:  &
        trivar_NLL_mean_eq  ! Procedure

    use parameters_KK, only: &
        KK_evap_Supersat_exp, & ! Variable(s)
        KK_evap_rr_exp,       &
        KK_evap_Nr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_eta_1,        & ! Mean of eta (old t) (1st PDF component)       [kg/kg]
      mu_eta_2,        & ! Mean of eta (old t) (2nd PDF component)       [kg/kg]
      mu_chi_1,        & ! Mean of chi (old s) (1st PDF component)       [kg/kg]
      mu_chi_2,        & ! Mean of chi (old s) (2nd PDF component)       [kg/kg]
      mu_rr_1,         & ! Mean of rr (1st PDF component) in-precip (ip) [kg/kg]
      mu_rr_2,         & ! Mean of rr (2nd PDF component) ip             [kg/kg]
      mu_Nr_1,         & ! Mean of Nr (1st PDF component) ip            [num/kg]
      mu_Nr_2,         & ! Mean of Nr (2nd PDF component) ip            [num/kg]
      mu_rr_1_n,       & ! Mean of ln rr (1st PDF component) ip      [ln(kg/kg)]
      mu_rr_2_n,       & ! Mean of ln rr (2nd PDF component) ip      [ln(kg/kg)]
      mu_Nr_1_n,       & ! Mean of ln Nr (1st PDF component) ip     [ln(num/kg)]
      mu_Nr_2_n,       & ! Mean of ln Nr (2nd PDF component) ip     [ln(num/kg)]
      sigma_eta_1,     & ! Standard deviation of eta (1st PDF component) [kg/kg]
      sigma_eta_2,     & ! Standard deviation of eta (2nd PDF component) [kg/kg]
      sigma_chi_1,     & ! Standard deviation of chi (1st PDF component) [kg/kg]
      sigma_chi_2,     & ! Standard deviation of chi (2nd PDF component) [kg/kg]
      sigma_rr_1,      & ! Standard deviation of rr (1st PDF comp.) ip   [kg/kg]
      sigma_rr_2,      & ! Standard deviation of rr (2nd PDF comp.) ip   [kg/kg]
      sigma_Nr_1,      & ! Standard deviation of Nr (1st PDF comp.) ip    [#/kg]
      sigma_Nr_2,      & ! Standard deviation of Nr (2nd PDF comp.) ip    [#/kg]
      sigma_rr_1_n,    & ! Standard deviation of ln rr (1st PDF comp.) ip    [-]
      sigma_rr_2_n,    & ! Standard deviation of ln rr (2nd PDF comp.) ip    [-]
      sigma_Nr_1_n,    & ! Standard deviation of ln Nr (1st PDF comp.) ip    [-]
      sigma_Nr_2_n,    & ! Standard deviation of ln Nr (2nd PDF comp.) ip    [-]
      corr_eta_chi_1,  & ! Correlation of eta and chi (1st PDF component)    [-]
      corr_eta_chi_2,  & ! Correlation of eta and chi (2nd PDF component)    [-]
      corr_eta_rr_1_n, & ! Correlation of eta and ln rr (1st PDF comp.) ip   [-]
      corr_eta_rr_2_n, & ! Correlation of eta and ln rr (2nd PDF comp.) ip   [-]
      corr_eta_Nr_1_n, & ! Correlation of eta and ln Nr (1st PDF comp.) ip   [-]
      corr_eta_Nr_2_n, & ! Correlation of eta and ln Nr (2nd PDF comp.) ip   [-]
      corr_chi_rr_1_n, & ! Correlation of chi and ln rr (1st PDF comp.) ip   [-]
      corr_chi_rr_2_n, & ! Correlation of chi and ln rr (2nd PDF comp.) ip   [-]
      corr_chi_Nr_1_n, & ! Correlation of chi and ln Nr (1st PDF comp.) ip   [-]
      corr_chi_Nr_2_n, & ! Correlation of chi and ln Nr (2nd PDF comp.) ip   [-]
      corr_rr_Nr_1_n,  & ! Correlation of ln rr & ln Nr (1st PDF comp.) ip   [-]
      corr_rr_Nr_2_n,  & ! Correlation of ln rr & ln Nr (2nd PDF comp.) ip   [-]
      mixt_frac,       & ! Mixture fraction                                  [-]
      precip_frac_1,   & ! Precipitation fraction (1st PDF component)        [-]
      precip_frac_2      ! Precipitation fraction (2nd PDF component)        [-]

    real( kind = core_rknd ), intent(in) :: &
      rtm,           & ! Mean of total water mixing ratio, rt (overall)  [kg/kg]
      mu_rt_1,       & ! Mean of rt (1st PDF component)                  [kg/kg]
      mu_rt_2,       & ! Mean of rt (2nd PDF component)                  [kg/kg]
      KK_evap_tndcy, & ! KK evaporation tendency                     [(kg/kg)/s]
      KK_evap_coef,  & ! KK evap. coef.  [(kg/kg)^(1-alpha-beta)(#/kg)^-gamma/s]
      eta_tol,       & ! Tolerance value of eta                          [kg/kg]
      crt1,          & ! Coefficient c_rt (1st PDF component)                [-]
      crt2             ! Coefficient c_rt (2nd PDF component)                [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_rt_KK_evap  ! Covariance of r_t and KK evap. tendency  [(kg/kg)^2/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp,      & ! Exponent on chi                                    [-]
      beta_exp,       & ! Exponent on r_r                                    [-]
      gamma_exp,      & ! Exponent on N_r                                    [-]
      comp_1_contrib, & ! Contribution to rt'KKevap' (PDF comp. 1) [(kg/kg)^2/s]
      comp_2_contrib    ! Contribution to rt'KKevap' (PDF comp. 2) [(kg/kg)^2/s]


    ! Values of the KK exponents.
    alpha_exp = KK_evap_Supersat_exp
    beta_exp  = KK_evap_rr_exp
    gamma_exp = KK_evap_Nr_exp

    ! Calculate the contribution from PDF component 1 to the covariance of
    ! r_t and KK evaporation tendency.
    comp_1_contrib  &
    = KK_evap_coef * precip_frac_1 &
      * ( ( one / ( two * crt1 ) )  &
          * ( quadrivar_NNLL_covar_eq( mu_eta_1, mu_chi_1, mu_rr_1, mu_Nr_1, &
                                       mu_rr_1_n, mu_Nr_1_n, sigma_eta_1, &
                                       sigma_chi_1, sigma_rr_1, sigma_Nr_1, &
                                       sigma_rr_1_n, sigma_Nr_1_n, &
                                       corr_eta_chi_1, corr_eta_rr_1_n, &
                                       corr_eta_Nr_1_n, corr_chi_rr_1_n, &
                                       corr_chi_Nr_1_n, corr_rr_Nr_1_n, &
                                       mu_eta_1, KK_evap_tndcy, &
                                       KK_evap_coef, eta_tol, &
                                       alpha_exp, beta_exp, gamma_exp )  &
              + trivar_NLL_mean_eq( mu_chi_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, &
                                    mu_Nr_1_n, sigma_chi_1, sigma_rr_1, &
                                    sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, &
                                    corr_chi_rr_1_n, corr_chi_Nr_1_n, &
                                    corr_rr_Nr_1_n, &
                                    alpha_exp + one, beta_exp, gamma_exp )  &
            ) &
          + ( mu_rt_1 - rtm - mu_chi_1 / ( two * crt1 ) )  &
            * trivar_NLL_mean_eq( mu_chi_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, &
                                  mu_Nr_1_n, sigma_chi_1, sigma_rr_1, &
                                  sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, &
                                  corr_chi_rr_1_n, corr_chi_Nr_1_n, &
                                  corr_rr_Nr_1_n, &
                                  alpha_exp, beta_exp, gamma_exp )  &
        )

    ! Calculate the contribution from PDF component 2 to the covariance of
    ! r_t and KK evaporation tendency.
    comp_2_contrib  &
    = KK_evap_coef * precip_frac_2 &
      * ( ( one / ( two * crt2 ) )  &
          * ( quadrivar_NNLL_covar_eq( mu_eta_2, mu_chi_2, mu_rr_2, mu_Nr_2, &
                                       mu_rr_2_n, mu_Nr_2_n, sigma_eta_2, &
                                       sigma_chi_2, sigma_rr_2, sigma_Nr_2, &
                                       sigma_rr_2_n, sigma_Nr_2_n, &
                                       corr_eta_chi_2, corr_eta_rr_2_n, &
                                       corr_eta_Nr_2_n, corr_chi_rr_2_n, &
                                       corr_chi_Nr_2_n, corr_rr_Nr_2_n, &
                                       mu_eta_2, KK_evap_tndcy, &
                                       KK_evap_coef, eta_tol, &
                                       alpha_exp, beta_exp, gamma_exp )  &
              + trivar_NLL_mean_eq( mu_chi_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, &
                                    mu_Nr_2_n, sigma_chi_2, sigma_rr_2, &
                                    sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, &
                                    corr_chi_rr_2_n, corr_chi_Nr_2_n, &
                                    corr_rr_Nr_2_n, &
                                    alpha_exp + one, beta_exp, gamma_exp )  &
            ) &
          + ( mu_rt_2 - rtm - mu_chi_2 / ( two * crt2 ) )  &
            * trivar_NLL_mean_eq( mu_chi_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, &
                                  mu_Nr_2_n, sigma_chi_2, sigma_rr_2, &
                                  sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, &
                                  corr_chi_rr_2_n, corr_chi_Nr_2_n, &
                                  corr_rr_Nr_2_n, &
                                  alpha_exp, beta_exp, gamma_exp )  &
        )

    ! Calculate the covariance of r_t and KK evaporation tendency.
    covar_rt_KK_evap  &
    = mixt_frac * comp_1_contrib + ( one - mixt_frac ) * comp_2_contrib


    return

  end function covar_rt_KK_evap

  !=============================================================================
  function covar_thl_KK_evap( mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_rr_1, &
                              mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, &
                              mu_Nr_1_n, mu_Nr_2_n, sigma_eta_1, sigma_eta_2, &
                              sigma_chi_1, sigma_chi_2, sigma_rr_1, &
                              sigma_rr_2, sigma_Nr_1, sigma_Nr_2, &
                              sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                              sigma_Nr_2_n, corr_eta_chi_1, corr_eta_chi_2, &
                              corr_eta_rr_1_n, corr_eta_rr_2_n, &
                              corr_eta_Nr_1_n, corr_eta_Nr_2_n, &
                              corr_chi_rr_1_n, corr_chi_rr_2_n, &
                              corr_chi_Nr_1_n, corr_chi_Nr_2_n, &
                              corr_rr_Nr_1_n, corr_rr_Nr_2_n, &
                              mixt_frac, precip_frac_1, precip_frac_2, &
                              thlm, mu_thl_1, mu_thl_2, KK_evap_tndcy, &
                              KK_evap_coef, eta_tol, cthl1, cthl2 )

    ! Description:

    ! References:
    ! Eq. (E27) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (A34) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
    ! microphysical effects on variances and covariances of moisture and heat
    ! content using a multivariate probability density function: a study with
    ! CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two, & ! Constant(s)
        one

    use KK_upscaled_means, only:  &
        trivar_NLL_mean_eq  ! Procedure

    use parameters_KK, only: &
        KK_evap_Supersat_exp, & ! Variable(s)
        KK_evap_rr_exp,       &
        KK_evap_Nr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_eta_1,        & ! Mean of eta (old t) (1st PDF component)       [kg/kg]
      mu_eta_2,        & ! Mean of eta (old t) (2nd PDF component)       [kg/kg]
      mu_chi_1,        & ! Mean of chi (old s) (1st PDF component)       [kg/kg]
      mu_chi_2,        & ! Mean of chi (old s) (2nd PDF component)       [kg/kg]
      mu_rr_1,         & ! Mean of rr (1st PDF component) in-precip (ip) [kg/kg]
      mu_rr_2,         & ! Mean of rr (2nd PDF component) ip             [kg/kg]
      mu_Nr_1,         & ! Mean of Nr (1st PDF component) ip            [num/kg]
      mu_Nr_2,         & ! Mean of Nr (2nd PDF component) ip            [num/kg]
      mu_rr_1_n,       & ! Mean of ln rr (1st PDF component) ip      [ln(kg/kg)]
      mu_rr_2_n,       & ! Mean of ln rr (2nd PDF component) ip      [ln(kg/kg)]
      mu_Nr_1_n,       & ! Mean of ln Nr (1st PDF component) ip     [ln(num/kg)]
      mu_Nr_2_n,       & ! Mean of ln Nr (2nd PDF component) ip     [ln(num/kg)]
      sigma_eta_1,     & ! Standard deviation of eta (1st PDF component) [kg/kg]
      sigma_eta_2,     & ! Standard deviation of eta (2nd PDF component) [kg/kg]
      sigma_chi_1,     & ! Standard deviation of chi (1st PDF component) [kg/kg]
      sigma_chi_2,     & ! Standard deviation of chi (2nd PDF component) [kg/kg]
      sigma_rr_1,      & ! Standard deviation of rr (1st PDF comp.) ip   [kg/kg]
      sigma_rr_2,      & ! Standard deviation of rr (2nd PDF comp.) ip   [kg/kg]
      sigma_Nr_1,      & ! Standard deviation of Nr (1st PDF comp.) ip    [#/kg]
      sigma_Nr_2,      & ! Standard deviation of Nr (2nd PDF comp.) ip    [#/kg]
      sigma_rr_1_n,    & ! Standard deviation of ln rr (1st PDF comp.) ip    [-]
      sigma_rr_2_n,    & ! Standard deviation of ln rr (2nd PDF comp.) ip    [-]
      sigma_Nr_1_n,    & ! Standard deviation of ln Nr (1st PDF comp.) ip    [-]
      sigma_Nr_2_n,    & ! Standard deviation of ln Nr (2nd PDF comp.) ip    [-]
      corr_eta_chi_1,  & ! Correlation of eta and chi (1st PDF component)    [-]
      corr_eta_chi_2,  & ! Correlation of eta and chi (2nd PDF component)    [-]
      corr_eta_rr_1_n, & ! Correlation of eta and ln rr (1st PDF comp.) ip   [-]
      corr_eta_rr_2_n, & ! Correlation of eta and ln rr (2nd PDF comp.) ip   [-]
      corr_eta_Nr_1_n, & ! Correlation of eta and ln Nr (1st PDF comp.) ip   [-]
      corr_eta_Nr_2_n, & ! Correlation of eta and ln Nr (2nd PDF comp.) ip   [-]
      corr_chi_rr_1_n, & ! Correlation of chi and ln rr (1st PDF comp.) ip   [-]
      corr_chi_rr_2_n, & ! Correlation of chi and ln rr (2nd PDF comp.) ip   [-]
      corr_chi_Nr_1_n, & ! Correlation of chi and ln Nr (1st PDF comp.) ip   [-]
      corr_chi_Nr_2_n, & ! Correlation of chi and ln Nr (2nd PDF comp.) ip   [-]
      corr_rr_Nr_1_n,  & ! Correlation of ln rr & ln Nr (1st PDF comp.) ip   [-]
      corr_rr_Nr_2_n,  & ! Correlation of ln rr & ln Nr (2nd PDF comp.) ip   [-]
      mixt_frac,       & ! Mixture fraction                                  [-]
      precip_frac_1,   & ! Precipitation fraction (1st PDF component)        [-]
      precip_frac_2      ! Precipitation fraction (2nd PDF component)        [-]

    real( kind = core_rknd ), intent(in) :: &
      thlm,          & ! Mean of liquid water pot. temp., thl (overall)      [K]
      mu_thl_1,      & ! Mean of thl (1st PDF component)                     [K]
      mu_thl_2,      & ! Mean of thl (2nd PDF component)                     [K]
      KK_evap_tndcy, & ! KK evaporation tendency                     [(kg/kg)/s]
      KK_evap_coef,  & ! KK evap. coef.  [(kg/kg)^(1-alpha-beta)(#/kg)^-gamma/s]
      eta_tol,       & ! Tolerance value of eta                          [kg/kg]
      cthl1,         & ! Coefficient c_thl (1st PDF component)       [(kg/kg)/K]
      cthl2            ! Coefficient c_thl (2nd PDF component)       [(kg/kg)/K]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_thl_KK_evap  ! Covariance of th_l and KK evap. tend.   [K*(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp,      & ! Exponent on chi                                    [-]
      beta_exp,       & ! Exponent on r_r                                    [-]
      gamma_exp,      & ! Exponent on N_r                                    [-]
      comp_1_contrib, & ! Contribution to thl'KKevap' (PDF comp. 1) [K(kg/kg)/s]
      comp_2_contrib    ! Contribution to thl'KKevap' (PDF comp. 2) [K(kg/kg)/s]


    ! Values of the KK exponents.
    alpha_exp = KK_evap_Supersat_exp
    beta_exp  = KK_evap_rr_exp
    gamma_exp = KK_evap_Nr_exp

    ! Calculate the contribution from PDF component 1 to the covariance of
    ! th_l and KK evaporation tendency.
    comp_1_contrib  &
    = KK_evap_coef * precip_frac_1 &
      * ( ( one / ( two * cthl1 ) )  &
          * ( quadrivar_NNLL_covar_eq( mu_eta_1, mu_chi_1, mu_rr_1, mu_Nr_1, &
                                       mu_rr_1_n, mu_Nr_1_n, sigma_eta_1, &
                                       sigma_chi_1, sigma_rr_1, sigma_Nr_1, &
                                       sigma_rr_1_n, sigma_Nr_1_n, &
                                       corr_eta_chi_1, corr_eta_rr_1_n, &
                                       corr_eta_Nr_1_n, corr_chi_rr_1_n, &
                                       corr_chi_Nr_1_n, corr_rr_Nr_1_n, &
                                       mu_eta_1, KK_evap_tndcy, &
                                       KK_evap_coef, eta_tol, &
                                       alpha_exp, beta_exp, gamma_exp )  &
              - trivar_NLL_mean_eq( mu_chi_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, &
                                    mu_Nr_1_n, sigma_chi_1, sigma_rr_1, &
                                    sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, &
                                    corr_chi_rr_1_n, corr_chi_Nr_1_n, &
                                    corr_rr_Nr_1_n, &
                                    alpha_exp + one, beta_exp, gamma_exp )  &
            ) &
          + ( mu_thl_1 - thlm + mu_chi_1 / ( two * cthl1 ) ) &
            * trivar_NLL_mean_eq( mu_chi_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, &
                                  mu_Nr_1_n, sigma_chi_1, sigma_rr_1, &
                                  sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, &
                                  corr_chi_rr_1_n, corr_chi_Nr_1_n, &
                                  corr_rr_Nr_1_n, &
                                  alpha_exp, beta_exp, gamma_exp )  &
        )

    ! Calculate the contribution from PDF component 2 to the covariance of
    ! th_l and KK evaporation tendency.
    comp_2_contrib  &
    = KK_evap_coef * precip_frac_2 &
      * ( ( one / ( two * cthl2 ) )  &
          * ( quadrivar_NNLL_covar_eq( mu_eta_2, mu_chi_2, mu_rr_2, mu_Nr_2, &
                                       mu_rr_2_n, mu_Nr_2_n, sigma_eta_2, &
                                       sigma_chi_2, sigma_rr_2, sigma_Nr_2, &
                                       sigma_rr_2_n, sigma_Nr_2_n, &
                                       corr_eta_chi_2, corr_eta_rr_2_n, &
                                       corr_eta_Nr_2_n, corr_chi_rr_2_n, &
                                       corr_chi_Nr_2_n, corr_rr_Nr_2_n, &
                                       mu_eta_2, KK_evap_tndcy, &
                                       KK_evap_coef, eta_tol, &
                                       alpha_exp, beta_exp, gamma_exp )  &
              - trivar_NLL_mean_eq( mu_chi_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, &
                                    mu_Nr_2_n, sigma_chi_2, sigma_rr_2, &
                                    sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, &
                                    corr_chi_rr_2_n, corr_chi_Nr_2_n, &
                                    corr_rr_Nr_2_n, &
                                    alpha_exp + one, beta_exp, gamma_exp )  &
            ) &
          + ( mu_thl_2 - thlm + mu_chi_2 / ( two * cthl2 ) )  &
            * trivar_NLL_mean_eq( mu_chi_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, &
                                  mu_Nr_2_n, sigma_chi_2, sigma_rr_2, &
                                  sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, &
                                  corr_chi_rr_2_n, corr_chi_Nr_2_n, &
                                  corr_rr_Nr_2_n, &
                                  alpha_exp, beta_exp, gamma_exp )  &
        )

    ! Calculate the covariance of th_l and KK evaporation tendency.
    covar_thl_KK_evap  &
    = mixt_frac * comp_1_contrib + ( one - mixt_frac ) * comp_2_contrib


    return

  end function covar_thl_KK_evap

  !=============================================================================
  function covar_x_KK_auto( mu_x_1, mu_x_2, mu_chi_1, mu_chi_2, mu_Ncn_1, &
                            mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n, sigma_x_1, &
                            sigma_x_2, sigma_chi_1, sigma_chi_2, sigma_Ncn_1, &
                            sigma_Ncn_2, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                            corr_x_chi_1, corr_x_chi_2, corr_x_Ncn_1_n, &
                            corr_x_Ncn_2_n, corr_chi_Ncn_1_n, &
                            corr_chi_Ncn_2_n, x_mean, KK_auto_tndcy, &
                            KK_auto_coef, x_tol, mixt_frac )

    ! Description:
    ! This function calculates the covariance of x and KK autoconversion
    ! tendency, which can be written as < x'((dr_r/dt)_KKauto)' >, or more
    ! simply as < x'KK_auto' >.

    ! References:
    ! Eq. (E3) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (A9) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
    ! microphysical effects on variances and covariances of moisture and heat
    ! content using a multivariate probability density function: a study with
    ! CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,     & ! Constant(s)
        Ncn_tol

    use parameters_KK, only: &
        KK_auto_rc_exp, & ! Variable(s)
        KK_auto_Nc_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x_1,           & ! Mean of x (1st PDF component)             [un. vary]
      mu_x_2,           & ! Mean of x (2nd PDF component)             [un. vary]
      mu_chi_1,         & ! Mean of chi (old s) (1st PDF component)      [kg/kg]
      mu_chi_2,         & ! Mean of chi (old s) (2nd PDF component)      [kg/kg]
      mu_Ncn_1,         & ! Mean of Ncn (1st PDF component)             [num/kg]
      mu_Ncn_2,         & ! Mean of Ncn (2nd PDF component)             [num/kg]
      mu_Ncn_1_n,       & ! Mean of ln Ncn (1st PDF component)      [ln(num/kg)]
      mu_Ncn_2_n,       & ! Mean of ln Ncn (2nd PDF component)      [ln(num/kg)]
      sigma_x_1,        & ! Standard deviation of x (1st PDF comp.)   [un. vary]
      sigma_x_2,        & ! Standard deviation of x (2nd PDF comp.)   [un. vary]
      sigma_chi_1,      & ! Standard deviation of chi (1st PDF comp.)    [kg/kg]
      sigma_chi_2,      & ! Standard deviation of chi (2nd PDF comp.)    [kg/kg]
      sigma_Ncn_1,      & ! Standard deviation of Ncn (1st PDF comp.)   [num/kg]
      sigma_Ncn_2,      & ! Standard deviation of Ncn (2nd PDF comp.)   [num/kg]
      sigma_Ncn_1_n,    & ! Standard deviation of ln Ncn (1st PDF component) [-]
      sigma_Ncn_2_n,    & ! Standard deviation of ln Ncn (2nd PDF component) [-]
      corr_x_chi_1,     & ! Correlation of x and chi (1st PDF component)     [-]
      corr_x_chi_2,     & ! Correlation of x and chi (2nd PDF component)     [-]
      corr_x_Ncn_1_n,   & ! Correlation of x and ln Ncn (1st PDF comp.)      [-]
      corr_x_Ncn_2_n,   & ! Correlation of x and ln Ncn (2nd PDF comp.)      [-]
      corr_chi_Ncn_1_n, & ! Correlation of chi and ln Ncn (1st PDF comp.)    [-]
      corr_chi_Ncn_2_n, & ! Correlation of chi and ln Ncn (2nd PDF comp.)    [-]
      x_mean,           & ! Mean of x (overall)                       [un. vary]
      KK_auto_tndcy,    & ! KK autoconversion tendency               [(kg/kg)/s]
      KK_auto_coef,     & ! KK auto. coef.[(kg/kg)^(1-alpha) (num/kg)^-beta / s]
      x_tol,            & ! Tolerance value of x                      [un. vary]
      mixt_frac           ! Mixture fraction                                 [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_x_KK_auto  ! Covariance of x and KK auto. tend. [u.v.(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on chi                                         [-]
      beta_exp     ! Exponent on N_cn                                        [-]


    ! Values of the KK exponents.
    alpha_exp = KK_auto_rc_exp
    beta_exp  = KK_auto_Nc_exp

    ! Calculate the covariance of x and KK autoconversion tendency.
    covar_x_KK_auto  &
    = KK_auto_coef &
      * ( mixt_frac &
          * trivar_NNL_covar_eq( mu_x_1, mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, &
                                 sigma_x_1, sigma_chi_1, sigma_Ncn_1, &
                                 sigma_Ncn_1_n, corr_x_chi_1, &
                                 corr_x_Ncn_1_n, corr_chi_Ncn_1_n, &
                                 x_mean, KK_auto_tndcy, &
                                 KK_auto_coef, x_tol, Ncn_tol, &
                                 alpha_exp, beta_exp ) &
        + ( one - mixt_frac ) &
          * trivar_NNL_covar_eq( mu_x_2, mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, &
                                 sigma_x_2, sigma_chi_2, sigma_Ncn_2, &
                                 sigma_Ncn_2_n, corr_x_chi_2, &
                                 corr_x_Ncn_2_n, corr_chi_Ncn_2_n, &
                                 x_mean, KK_auto_tndcy, &
                                 KK_auto_coef, x_tol, Ncn_tol, &
                                 alpha_exp, beta_exp ) &
        )


    return

  end function covar_x_KK_auto

  !=============================================================================
  function covar_rt_KK_auto( mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_Ncn_1, &
                             mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n, sigma_eta_1, &
                             sigma_eta_2, sigma_chi_1, sigma_chi_2, &
                             sigma_Ncn_1, sigma_Ncn_2, sigma_Ncn_1_n, &
                             sigma_Ncn_2_n, corr_eta_chi_1, corr_eta_chi_2, &
                             corr_eta_Ncn_1_n, corr_eta_Ncn_2_n, &
                             corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, &
                             rtm, mu_rt_1, mu_rt_2, KK_auto_tndcy, &
                             KK_auto_coef, eta_tol, crt1, crt2, mixt_frac )

    ! Description:

    ! References:
    ! Eq. (E6) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (A12) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
    ! microphysical effects on variances and covariances of moisture and heat
    ! content using a multivariate probability density function: a study with
    ! CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two,     & ! Constant(s)
        one,     & 
        Ncn_tol, &
        Nc_tol

    use KK_upscaled_means, only:  &
        bivar_NL_mean_eq    ! Procedure(s)

    use parameters_KK, only: &
        KK_auto_rc_exp, & ! Variable(s)
        KK_auto_Nc_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_eta_1,         & ! Mean of eta (old t) (1st PDF component)      [kg/kg]
      mu_eta_2,         & ! Mean of eta (old t) (2nd PDF component)      [kg/kg]
      mu_chi_1,         & ! Mean of chi (old s) (1st PDF component)      [kg/kg]
      mu_chi_2,         & ! Mean of chi (old s) (2nd PDF component)      [kg/kg]
      mu_Ncn_1,         & ! Mean of Ncn (1st PDF component)             [num/kg]
      mu_Ncn_2,         & ! Mean of Ncn (2nd PDF component)             [num/kg]
      mu_Ncn_1_n,       & ! Mean of ln Ncn (1st PDF component)      [ln(num/kg)]
      mu_Ncn_2_n,       & ! Mean of ln Ncn (2nd PDF component)      [ln(num/kg)]
      sigma_eta_1,      & ! Standard deviation of eta (1st PDF comp.)    [kg/kg]
      sigma_eta_2,      & ! Standard deviation of eta (2nd PDF comp.)    [kg/kg]
      sigma_chi_1,      & ! Standard deviation of chi (1st PDF comp.)    [kg/kg]
      sigma_chi_2,      & ! Standard deviation of chi (2nd PDF comp.)    [kg/kg]
      sigma_Ncn_1,      & ! Standard deviation of Ncn (1st PDF comp.)   [num/kg]
      sigma_Ncn_2,      & ! Standard deviation of Ncn (2nd PDF comp.)   [num/kg]
      sigma_Ncn_1_n,    & ! Standard deviation of ln Ncn (1st PDF component) [-]
      sigma_Ncn_2_n,    & ! Standard deviation of ln Ncn (2nd PDF component) [-]
      corr_eta_chi_1,   & ! Correlation of eta and chi (1st PDF component)   [-]
      corr_eta_chi_2,   & ! Correlation of eta and chi (2nd PDF component)   [-]
      corr_eta_Ncn_1_n, & ! Correlation of eta and ln Ncn (1st PDF comp.)    [-]
      corr_eta_Ncn_2_n, & ! Correlation of eta and ln Ncn (2nd PDF comp.)    [-]
      corr_chi_Ncn_1_n, & ! Correlation of chi and ln Ncn (1st PDF comp.)    [-]
      corr_chi_Ncn_2_n, & ! Correlation of chi and ln Ncn (2nd PDF comp.)    [-]
      rtm,              & ! Mean of total water mix. ratio, rt (overall) [kg/kg]
      mu_rt_1,          & ! Mean of rt (1st PDF component)               [kg/kg]
      mu_rt_2,          & ! Mean of rt (2nd PDF component)               [kg/kg]
      KK_auto_tndcy,    & ! KK autoconversion tendency               [(kg/kg)/s]
      KK_auto_coef,     & ! KK auto. coef.  [(kg/kg)^(1-alpha) (num/kg)^-beta/s]
      eta_tol,          & ! Tolerance value of eta                       [kg/kg]
      crt1,             & ! Coefficient c_rt (1st PDF component)             [-]
      crt2,             & ! Coefficient c_rt (2nd PDF component)             [-]
      mixt_frac           ! Mixture fraction                                 [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_rt_KK_auto  ! Covariance of r_t and KK auto. tendency  [(kg/kg)^2/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp,      & ! Exponent on chi                                    [-]
      beta_exp,       & ! Exponent on N_cn                                   [-]
      comp_1_contrib, & ! Contribution to rt'KKauto' (PDF comp. 1) [(kg/kg)^2/s]
      comp_2_contrib    ! Contribution to rt'KKauto' (PDF comp. 2) [(kg/kg)^2/s]


    ! Values of the KK exponents.
    alpha_exp = KK_auto_rc_exp
    beta_exp  = KK_auto_Nc_exp

    ! Calculate the contribution from PDF component 1 to the covariance of
    ! r_t and KK autoconversion tendency.
    comp_1_contrib  &
    = KK_auto_coef  &
      * ( ( one / ( two * crt1 ) )  &
          * ( trivar_NNL_covar_eq( mu_eta_1, mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, &
                                   sigma_eta_1, sigma_chi_1, sigma_Ncn_1, &
                                   sigma_Ncn_1_n, corr_eta_chi_1, &
                                   corr_eta_Ncn_1_n, corr_chi_Ncn_1_n, &
                                   mu_eta_1, KK_auto_tndcy, &
                                   KK_auto_coef, eta_tol, Ncn_tol, &
                                   alpha_exp, beta_exp )  &
              + bivar_NL_mean_eq( mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, sigma_chi_1, &
                                  sigma_Ncn_1, sigma_Ncn_1_n, corr_chi_Ncn_1_n,&
                                  Nc_tol, alpha_exp + one, beta_exp )  &
            ) &
          + ( mu_rt_1 - rtm - mu_chi_1 / ( two * crt1 ) )  &
            * bivar_NL_mean_eq( mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, sigma_chi_1, &
                                sigma_Ncn_1, sigma_Ncn_1_n, corr_chi_Ncn_1_n, &
                                Nc_tol, alpha_exp, beta_exp )  &
        )

    ! Calculate the contribution from PDF component 2 to the covariance of
    ! r_t and KK autoconversion tendency.
    comp_2_contrib  &
    = KK_auto_coef  &
      * ( ( one / ( two * crt2 ) )  &
          * ( trivar_NNL_covar_eq( mu_eta_2, mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, &
                                   sigma_eta_2, sigma_chi_2, sigma_Ncn_2, &
                                   sigma_Ncn_2_n, corr_eta_chi_2, &
                                   corr_eta_Ncn_2_n, corr_chi_Ncn_2_n, &
                                   mu_eta_2, KK_auto_tndcy, &
                                   KK_auto_coef, eta_tol, Ncn_tol, &
                                   alpha_exp, beta_exp )  &
              + bivar_NL_mean_eq( mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, sigma_chi_2, &
                                  sigma_Ncn_2, sigma_Ncn_2_n, corr_chi_Ncn_2_n,&
                                  Nc_tol, alpha_exp + one, beta_exp )  &
            ) &
          + ( mu_rt_2 - rtm - mu_chi_2 / ( two * crt2 ) )  &
            * bivar_NL_mean_eq( mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, sigma_chi_2, &
                                sigma_Ncn_2, sigma_Ncn_2_n, corr_chi_Ncn_2_n, &
                                Nc_tol, alpha_exp, beta_exp )  &
        )

    ! Calculate the covariance of r_t and KK autoconversion tendency.
    covar_rt_KK_auto  &
    = mixt_frac * comp_1_contrib + ( one - mixt_frac ) * comp_2_contrib


    return

  end function covar_rt_KK_auto

  !=============================================================================
  function covar_thl_KK_auto( mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, &
                              mu_Ncn_1, mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n, &
                              sigma_eta_1, sigma_eta_2, sigma_chi_1, &
                              sigma_chi_2, sigma_Ncn_1, sigma_Ncn_2, &
                              sigma_Ncn_1_n, sigma_Ncn_2_n, corr_eta_chi_1, &
                              corr_eta_chi_2, corr_eta_Ncn_1_n, &
                              corr_eta_Ncn_2_n, corr_chi_Ncn_1_n, &
                              corr_chi_Ncn_2_n, thlm, mu_thl_1, mu_thl_2, &
                              KK_auto_tndcy, KK_auto_coef, eta_tol, &
                              cthl1, cthl2, mixt_frac )

    ! Description:

    ! References:
    ! Eq. (E9) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (A15) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
    ! microphysical effects on variances and covariances of moisture and heat
    ! content using a multivariate probability density function: a study with
    ! CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two,     & ! Constant(s)
        one,     &
        Ncn_tol, &
        Nc_tol

    use KK_upscaled_means, only:  &
        bivar_NL_mean_eq    ! Procedure(s)

    use parameters_KK, only: &
        KK_auto_rc_exp, & ! Variable(s)
        KK_auto_Nc_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_eta_1,         & ! Mean of eta (old t) (1st PDF component)      [kg/kg]
      mu_eta_2,         & ! Mean of eta (old t) (2nd PDF component)      [kg/kg]
      mu_chi_1,         & ! Mean of chi (old s) (1st PDF component)      [kg/kg]
      mu_chi_2,         & ! Mean of chi (old s) (2nd PDF component)      [kg/kg]
      mu_Ncn_1,         & ! Mean of Ncn (1st PDF component)             [num/kg]
      mu_Ncn_2,         & ! Mean of Ncn (2nd PDF component)             [num/kg]
      mu_Ncn_1_n,       & ! Mean of ln Ncn (1st PDF component)      [ln(num/kg)]
      mu_Ncn_2_n,       & ! Mean of ln Ncn (2nd PDF component)      [ln(num/kg)]
      sigma_eta_1,      & ! Standard deviation of eta (1st PDF comp.)    [kg/kg]
      sigma_eta_2,      & ! Standard deviation of eta (2nd PDF comp.)    [kg/kg]
      sigma_chi_1,      & ! Standard deviation of chi (1st PDF comp.)    [kg/kg]
      sigma_chi_2,      & ! Standard deviation of chi (2nd PDF comp.)    [kg/kg]
      sigma_Ncn_1,      & ! Standard deviation of Ncn (1st PDF comp.)   [num/kg]
      sigma_Ncn_2,      & ! Standard deviation of Ncn (2nd PDF comp.)   [num/kg]
      sigma_Ncn_1_n,    & ! Standard deviation of ln Ncn (1st PDF component) [-]
      sigma_Ncn_2_n,    & ! Standard deviation of ln Ncn (2nd PDF component) [-]
      corr_eta_chi_1,   & ! Correlation of eta and chi (1st PDF component)   [-]
      corr_eta_chi_2,   & ! Correlation of eta and chi (2nd PDF component)   [-]
      corr_eta_Ncn_1_n, & ! Correlation of eta and ln Ncn (1st PDF comp.)    [-]
      corr_eta_Ncn_2_n, & ! Correlation of eta and ln Ncn (2nd PDF comp.)    [-]
      corr_chi_Ncn_1_n, & ! Correlation of chi and ln Ncn (1st PDF comp.)    [-]
      corr_chi_Ncn_2_n, & ! Correlation of chi and ln Ncn (2nd PDF comp.)    [-]
      thlm,             & ! Mean of liquid water pot. temp., thl (overall)   [K]
      mu_thl_1,         & ! Mean of thl (1st PDF component)                  [K]
      mu_thl_2,         & ! Mean of thl (2nd PDF component)                  [K]
      KK_auto_tndcy,    & ! KK autoconversion tendency               [(kg/kg)/s]
      KK_auto_coef,     & ! KK auto. coef.  [(kg/kg)^(1-alpha) (num/kg)^-beta/s]
      eta_tol,          & ! Tolerance value of eta                       [kg/kg]
      cthl1,            & ! Coefficient c_thl (1st PDF component)    [(kg/kg)/K]
      cthl2,            & ! Coefficient c_thl (2nd PDF component)    [(kg/kg)/K]
      mixt_frac           ! Mixture fraction                                 [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_thl_KK_auto  ! Covariance of th_l and KK auto. tendency [K(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp,      & ! Exponent on chi                                    [-]
      beta_exp,       & ! Exponent on N_c                                    [-]
      comp_1_contrib, & ! Contribution to thl'KKauto' (PDF comp. 1) [K(kg/kg)/s]
      comp_2_contrib    ! Contribution to thl'KKauto' (PDF comp. 2) [K(kg/kg)/s]


    ! Values of the KK exponents.
    alpha_exp = KK_auto_rc_exp
    beta_exp  = KK_auto_Nc_exp

    ! Calculate the contribution from PDF component 1 to the covariance of
    ! th_l and KK autoconversion tendency.
    comp_1_contrib  &
    = KK_auto_coef  &
      * ( ( one / ( two * cthl1 ) )  &
          * ( trivar_NNL_covar_eq( mu_eta_1, mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, &
                                   sigma_eta_1, sigma_chi_1, sigma_Ncn_1, &
                                   sigma_Ncn_1_n, corr_eta_chi_1, &
                                   corr_eta_Ncn_1_n, corr_chi_Ncn_1_n, &
                                   mu_eta_1, KK_auto_tndcy, &
                                   KK_auto_coef, eta_tol, Ncn_tol, &
                                   alpha_exp, beta_exp )  &
              - bivar_NL_mean_eq( mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, sigma_chi_1, &
                                  sigma_Ncn_1, sigma_Ncn_1_n, corr_chi_Ncn_1_n,&
                                  Nc_tol, alpha_exp + one, beta_exp )  &
            ) &
          + ( mu_thl_1 - thlm + mu_chi_1 / ( two * cthl1 ) )  &
            * bivar_NL_mean_eq( mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, sigma_chi_1, &
                                sigma_Ncn_1, sigma_Ncn_1_n, corr_chi_Ncn_1_n, &
                                Nc_tol, alpha_exp, beta_exp )  &
        )

    ! Calculate the contribution from PDF component 2 to the covariance of
    ! th_l and KK autoconversion tendency.
    comp_2_contrib  &
    = KK_auto_coef  &
      * ( ( one / ( two * cthl2 ) )  &
          * ( trivar_NNL_covar_eq( mu_eta_2, mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, &
                                   sigma_eta_2, sigma_chi_2, sigma_Ncn_2, &
                                   sigma_Ncn_2_n, corr_eta_chi_2, &
                                   corr_eta_Ncn_2_n, corr_chi_Ncn_2_n, &
                                   mu_eta_2, KK_auto_tndcy, &
                                   KK_auto_coef, eta_tol, Ncn_tol, &
                                   alpha_exp, beta_exp )  &
              - bivar_NL_mean_eq( mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, sigma_chi_2, &
                                  sigma_Ncn_2, sigma_Ncn_2_n, corr_chi_Ncn_2_n,&
                                  Nc_tol, alpha_exp + one, beta_exp )  &
            ) &
          + ( mu_thl_2 - thlm + mu_chi_2 / ( two * cthl2 ) )  &
            * bivar_NL_mean_eq( mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, sigma_chi_2, &
                                sigma_Ncn_2, sigma_Ncn_2_n, corr_chi_Ncn_2_n, &
                                Nc_tol, alpha_exp, beta_exp )  &    
        )

    ! Calculate the covariance of th_l and KK autoconversion tendency.
    covar_thl_KK_auto  &
    = mixt_frac * comp_1_contrib + ( one - mixt_frac ) * comp_2_contrib


    return

  end function covar_thl_KK_auto

  !=============================================================================
  function covar_x_KK_accr( mu_x_1, mu_x_2, mu_chi_1, mu_chi_2, mu_rr_1, &
                            mu_rr_2, mu_rr_1_n, mu_rr_2_n, sigma_x_1, &
                            sigma_x_2, sigma_chi_1, sigma_chi_2, sigma_rr_1, &
                            sigma_rr_2, sigma_rr_1_n, sigma_rr_2_n, &
                            corr_x_chi_1, corr_x_chi_2, corr_x_rr_1_n, &
                            corr_x_rr_2_n, corr_chi_rr_1_n, corr_chi_rr_2_n, &
                            x_mean, KK_accr_tndcy, KK_accr_coef, x_tol, &
                            mixt_frac, precip_frac_1, precip_frac_2 )

    ! Description:
    ! This function calculates the covariance of x and KK accretion
    ! tendency, which can be written as < x'((dr_r/dt)_KKaccr)' >, or more
    ! simply as < x'KK_accr' >.

    ! References:
    ! Eq. (E12) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (A18) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
    ! microphysical effects on variances and covariances of moisture and heat
    ! content using a multivariate probability density function: a study with
    ! CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,    & ! Constant(s)
        rr_tol

    use parameters_KK, only: &
        KK_accr_rc_exp, & ! Variable(s)
        KK_accr_rr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x_1,          & ! Mean of x (1st PDF component)   [units vary (un. v.)]
      mu_x_2,          & ! Mean of x (2nd PDF component)   [units vary (un. v.)]
      mu_chi_1,        & ! Mean of chi (old s) (1st PDF component)       [kg/kg]
      mu_chi_2,        & ! Mean of chi (old s) (2nd PDF component)       [kg/kg]
      mu_rr_1,         & ! Mean of rr (1st PDF component) in-precip (ip) [kg/kg]
      mu_rr_2,         & ! Mean of rr (2nd PDF component) ip             [kg/kg]
      mu_rr_1_n,       & ! Mean of ln rr (1st PDF component) ip      [ln(kg/kg)]
      mu_rr_2_n,       & ! Mean of ln rr (2nd PDF component) ip      [ln(kg/kg)]
      sigma_x_1,       & ! Standard deviation of x (1st PDF component)  [un. v.]
      sigma_x_2,       & ! Standard deviation of x (2nd PDF component)  [un. v.]
      sigma_chi_1,     & ! Standard deviation of chi (1st PDF component) [kg/kg]
      sigma_chi_2,     & ! Standard deviation of chi (2nd PDF component) [kg/kg]
      sigma_rr_1,      & ! Standard deviation of rr (1st PDF comp.) ip   [kg/kg]
      sigma_rr_2,      & ! Standard deviation of rr (2nd PDF comp.) ip   [kg/kg]
      sigma_rr_1_n,    & ! Standard deviation of ln rr (1st PDF comp.) ip    [-]
      sigma_rr_2_n,    & ! Standard deviation of ln rr (2nd PDF comp.) ip    [-]
      corr_x_chi_1,    & ! Correlation of x and chi (1st PDF component)      [-]
      corr_x_chi_2,    & ! Correlation of x and chi (2nd PDF component)      [-]
      corr_x_rr_1_n,   & ! Correlation of x and ln rr (1st PDF comp.) ip     [-]
      corr_x_rr_2_n,   & ! Correlation of x and ln rr (2nd PDF comp.) ip     [-]
      corr_chi_rr_1_n, & ! Correlation of chi and ln rr (1st PDF comp.) ip   [-]
      corr_chi_rr_2_n, & ! Correlation of chi and ln rr (2nd PDF comp.) ip   [-]
      x_mean,          & ! Mean of x (overall)                      [units vary]
      KK_accr_tndcy,   & ! KK accretion tendency                     [(kg/kg)/s]
      KK_accr_coef,    & ! KK accretion coefficient                  [(kg/kg)/s]
      x_tol,           & ! Tolerance value of x                     [units vary]
      mixt_frac,       & ! Mixture fraction                                  [-]
      precip_frac_1,   & ! Precipitation fraction (1st PDF component)        [-]
      precip_frac_2      ! Precipitation fraction (2nd PDF component)        [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_x_KK_accr  ! Covariance of x and KK accr. tendency   [u.v.(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on chi                                         [-]
      beta_exp     ! Exponent on r_r                                         [-]


    ! Values of the KK exponents.
    alpha_exp = KK_accr_rc_exp
    beta_exp  = KK_accr_rr_exp

    ! Calculate the covariance of x and KK accretion tendency.
    covar_x_KK_accr  &
    = mixt_frac &
      * ( KK_accr_coef * precip_frac_1 &
          * trivar_NNL_covar_eq( mu_x_1, mu_chi_1, mu_rr_1, mu_rr_1_n, &
                                 sigma_x_1, sigma_chi_1, sigma_rr_1, &
                                 sigma_rr_1_n, corr_x_chi_1, &
                                 corr_x_rr_1_n, corr_chi_rr_1_n, &
                                 x_mean, KK_accr_tndcy, &
                                 KK_accr_coef, x_tol, rr_tol, &
                                 alpha_exp, beta_exp )  &
          - ( one - precip_frac_1 ) * ( mu_x_1 - x_mean ) * KK_accr_tndcy &
        ) &
      + ( one - mixt_frac ) &
        * ( KK_accr_coef * precip_frac_2 &
            * trivar_NNL_covar_eq( mu_x_2, mu_chi_2, mu_rr_2, mu_rr_2_n, &
                                   sigma_x_2, sigma_chi_2, sigma_rr_2, &
                                   sigma_rr_2_n, corr_x_chi_2, &
                                   corr_x_rr_2_n, corr_chi_rr_2_n, &
                                   x_mean, KK_accr_tndcy, &
                                   KK_accr_coef, x_tol, rr_tol, &
                                   alpha_exp, beta_exp )  &
            - ( one - precip_frac_2 ) * ( mu_x_2 - x_mean ) * KK_accr_tndcy &
          )


    return

  end function covar_x_KK_accr

  !=============================================================================
  function covar_rt_KK_accr( mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_rr_1, &
                             mu_rr_2, mu_rr_1_n, mu_rr_2_n, sigma_eta_1, &
                             sigma_eta_2, sigma_chi_1, sigma_chi_2, &
                             sigma_rr_1, sigma_rr_2, sigma_rr_1_n, &
                             sigma_rr_2_n, corr_eta_chi_1, corr_eta_chi_2, &
                             corr_eta_rr_1_n, corr_eta_rr_2_n, &
                             corr_chi_rr_1_n, corr_chi_rr_2_n, &
                             rtm, mu_rt_1, mu_rt_2, KK_accr_tndcy, &
                             KK_accr_coef, eta_tol, crt1, crt2, mixt_frac, &
                             precip_frac_1, precip_frac_2 )

    ! Description:

    ! References:
    ! Eq. (E15) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (A21) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
    ! microphysical effects on variances and covariances of moisture and heat
    ! content using a multivariate probability density function: a study with
    ! CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two,    & ! Constant(s)
        one,    &
        rr_tol

    use KK_upscaled_means, only:  &
        bivar_NL_mean_eq

    use parameters_KK, only: &
        KK_accr_rc_exp, & ! Variable(s)
        KK_accr_rr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_eta_1,        & ! Mean of eta (old t) (1st PDF component)       [kg/kg]
      mu_eta_2,        & ! Mean of eta (old t) (2nd PDF component)       [kg/kg]
      mu_chi_1,        & ! Mean of chi (old s) (1st PDF component)       [kg/kg]
      mu_chi_2,        & ! Mean of chi (old s) (2nd PDF component)       [kg/kg]
      mu_rr_1,         & ! Mean of rr (1st PDF component) in-precip (ip) [kg/kg]
      mu_rr_2,         & ! Mean of rr (2nd PDF component) ip             [kg/kg]
      mu_rr_1_n,       & ! Mean of ln rr (1st PDF component) ip      [ln(kg/kg)]
      mu_rr_2_n,       & ! Mean of ln rr (2nd PDF component) ip      [ln(kg/kg)]
      sigma_eta_1,     & ! Standard deviation of eta (1st PDF component) [kg/kg]
      sigma_eta_2,     & ! Standard deviation of eta (2nd PDF component) [kg/kg]
      sigma_chi_1,     & ! Standard deviation of chi (1st PDF component) [kg/kg]
      sigma_chi_2,     & ! Standard deviation of chi (2nd PDF component) [kg/kg]
      sigma_rr_1,      & ! Standard deviation of rr (1st PDF comp.) ip   [kg/kg]
      sigma_rr_2,      & ! Standard deviation of rr (2nd PDF comp.) ip   [kg/kg]
      sigma_rr_1_n,    & ! Standard deviation of ln rr (1st PDF comp.) ip    [-]
      sigma_rr_2_n,    & ! Standard deviation of ln rr (2nd PDF comp.) ip    [-]
      corr_eta_chi_1,  & ! Correlation of eta and chi (1st PDF component)    [-]
      corr_eta_chi_2,  & ! Correlation of eta and chi (2nd PDF component)    [-]
      corr_eta_rr_1_n, & ! Correlation of eta and ln rr (1st PDF comp.) ip   [-]
      corr_eta_rr_2_n, & ! Correlation of eta and ln rr (2nd PDF comp.) ip   [-]
      corr_chi_rr_1_n, & ! Correlation of chi and ln rr (1st PDF comp.) ip   [-]
      corr_chi_rr_2_n, & ! Correlation of chi and ln rr (2nd PDF comp.) ip   [-]
      rtm,             & ! Mean of total water mix. ratio, rt (overall)  [kg/kg]
      mu_rt_1,         & ! Mean of rt (1st PDF component)                [kg/kg]
      mu_rt_2,         & ! Mean of rt (2nd PDF component)                [kg/kg]
      KK_accr_tndcy,   & ! KK accretion tendency                     [(kg/kg)/s]
      KK_accr_coef,    & ! KK accretion coefficient [(kg/kg)^(1-alpha-beta) / s]
      eta_tol,         & ! Tolerance value of eta                        [kg/kg]
      crt1,            & ! Coefficient c_rt (1st PDF component)              [-]
      crt2,            & ! Coefficient c_rt (2nd PDF component)              [-]
      mixt_frac,       & ! Mixture fraction                                  [-]
      precip_frac_1,   & ! Precipitation fraction (1st PDF component)        [-]
      precip_frac_2      ! Precipitation fraction (2nd PDF component)        [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_rt_KK_accr  ! Covariance of r_t and KK accr. tendency  [(kg/kg)^2/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp,      & ! Exponent on chi                                    [-]
      beta_exp,       & ! Exponent on r_r                                    [-]
      comp_1_contrib, & ! Contribution to rt'KKaccr' (PDF comp. 1) [(kg/kg)^2/s]
      comp_2_contrib    ! Contribution to rt'KKaccr' (PDF comp. 2) [(kg/kg)^2/s]


    ! Values of the KK exponents.
    alpha_exp = KK_accr_rc_exp
    beta_exp  = KK_accr_rr_exp

    ! Calculate the contribution from PDF component 1 to the covariance of
    ! r_t and KK accretion tendency.
    comp_1_contrib  &
    = KK_accr_coef * precip_frac_1 &
      * ( ( one / ( two * crt1 ) )  &
          * ( trivar_NNL_covar_eq( mu_eta_1, mu_chi_1, mu_rr_1, mu_rr_1_n, &
                                   sigma_eta_1, sigma_chi_1, sigma_rr_1, &
                                   sigma_rr_1_n, corr_eta_chi_1, &
                                   corr_eta_rr_1_n, corr_chi_rr_1_n, &
                                   mu_eta_1, KK_accr_tndcy, &
                                   KK_accr_coef, eta_tol, rr_tol, &
                                   alpha_exp, beta_exp )  &
              + bivar_NL_mean_eq( mu_chi_1, mu_rr_1, mu_rr_1_n, sigma_chi_1, &
                                  sigma_rr_1, sigma_rr_1_n, corr_chi_rr_1_n, &
                                  rr_tol, alpha_exp + one, beta_exp )  &
            ) &
          + ( mu_rt_1 - rtm - mu_chi_1 / ( two * crt1 ) )  &
            * bivar_NL_mean_eq( mu_chi_1, mu_rr_1, mu_rr_1_n, sigma_chi_1, &
                                sigma_rr_1, sigma_rr_1_n, corr_chi_rr_1_n, &
                                rr_tol, alpha_exp, beta_exp )  &
        )

    ! Calculate the contribution from PDF component 2 to the covariance of
    ! r_t and KK accretion tendency.
    comp_2_contrib  &
    = KK_accr_coef * precip_frac_2 &
      * ( ( one / ( two * crt2 ) )  &
          * ( trivar_NNL_covar_eq( mu_eta_2, mu_chi_2, mu_rr_2, mu_rr_2_n, &
                                   sigma_eta_2, sigma_chi_2, sigma_rr_2, &
                                   sigma_rr_2_n, corr_eta_chi_2, &
                                   corr_eta_rr_2_n, corr_chi_rr_2_n, &
                                   mu_eta_2, KK_accr_tndcy, &
                                   KK_accr_coef, eta_tol, rr_tol, &
                                   alpha_exp, beta_exp )  &
              + bivar_NL_mean_eq( mu_chi_2, mu_rr_2, mu_rr_2_n, sigma_chi_2, &
                                  sigma_rr_2, sigma_rr_2_n, corr_chi_rr_2_n, &
                                  rr_tol, alpha_exp + one, beta_exp )  &
            ) &
          + ( mu_rt_2 - rtm - mu_chi_2 / ( two * crt2 ) )  &
            * bivar_NL_mean_eq( mu_chi_2, mu_rr_2, mu_rr_2_n, sigma_chi_2, &
                                sigma_rr_2, sigma_rr_2_n, corr_chi_rr_2_n, &
                                rr_tol, alpha_exp, beta_exp )  &
        )

    ! Calculate the covariance of r_t and KK accretion tendency.
    covar_rt_KK_accr  &
    = mixt_frac * comp_1_contrib + ( one - mixt_frac ) * comp_2_contrib


    return

  end function covar_rt_KK_accr

  !=============================================================================
  function covar_thl_KK_accr( mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_rr_1, &
                              mu_rr_2, mu_rr_1_n, mu_rr_2_n, sigma_eta_1, &
                              sigma_eta_2, sigma_chi_1, sigma_chi_2, &
                              sigma_rr_1, sigma_rr_2, sigma_rr_1_n, &
                              sigma_rr_2_n, corr_eta_chi_1, corr_eta_chi_2, &
                              corr_eta_rr_1_n, corr_eta_rr_2_n, &
                              corr_chi_rr_1_n, corr_chi_rr_2_n, &
                              thlm, mu_thl_1, mu_thl_2, KK_accr_tndcy, &
                              KK_accr_coef, eta_tol, cthl1, cthl2, mixt_frac, &
                              precip_frac_1, precip_frac_2 )

    ! Description:

    ! References:
    ! Eq. (E18) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (A24) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
    ! microphysical effects on variances and covariances of moisture and heat
    ! content using a multivariate probability density function: a study with
    ! CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two,    & ! Constant(s)
        one,    &
        rr_tol

    use KK_upscaled_means, only:  &
        bivar_NL_mean_eq

    use parameters_KK, only: &
        KK_accr_rc_exp, & ! Variable(s)
        KK_accr_rr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_eta_1,        & ! Mean of eta (old t) (1st PDF component)       [kg/kg]
      mu_eta_2,        & ! Mean of eta (old t) (2nd PDF component)       [kg/kg]
      mu_chi_1,        & ! Mean of chi (old s) (1st PDF component)       [kg/kg]
      mu_chi_2,        & ! Mean of chi (old s) (2nd PDF component)       [kg/kg]
      mu_rr_1,         & ! Mean of rr (1st PDF component) in-precip (ip) [kg/kg]
      mu_rr_2,         & ! Mean of rr (2nd PDF component) ip             [kg/kg]
      mu_rr_1_n,       & ! Mean of ln rr (1st PDF component) ip      [ln(kg/kg)]
      mu_rr_2_n,       & ! Mean of ln rr (2nd PDF component) ip      [ln(kg/kg)]
      sigma_eta_1,     & ! Standard deviation of eta (1st PDF component) [kg/kg]
      sigma_eta_2,     & ! Standard deviation of eta (2nd PDF component) [kg/kg]
      sigma_chi_1,     & ! Standard deviation of chi (1st PDF component) [kg/kg]
      sigma_chi_2,     & ! Standard deviation of chi (2nd PDF component) [kg/kg]
      sigma_rr_1,      & ! Standard deviation of rr (1st PDF comp.) ip   [kg/kg]
      sigma_rr_2,      & ! Standard deviation of rr (2nd PDF comp.) ip   [kg/kg]
      sigma_rr_1_n,    & ! Standard deviation of ln rr (1st PDF comp.) ip    [-]
      sigma_rr_2_n,    & ! Standard deviation of ln rr (2nd PDF comp.) ip    [-]
      corr_eta_chi_1,  & ! Correlation of eta and chi (1st PDF component)    [-]
      corr_eta_chi_2,  & ! Correlation of eta and chi (2nd PDF component)    [-]
      corr_eta_rr_1_n, & ! Correlation of eta and ln rr (1st PDF component)  [-]
      corr_eta_rr_2_n, & ! Correlation of eta and ln rr (2nd PDF component)  [-]
      corr_chi_rr_1_n, & ! Correlation of chi and ln rr (1st PDF component)  [-]
      corr_chi_rr_2_n, & ! Correlation of chi and ln rr (2nd PDF component)  [-]
      thlm,            & ! Mean of liquid water pot. temp., thl (overall)    [K]
      mu_thl_1,        & ! Mean of thl (1st PDF component)                   [K]
      mu_thl_2,        & ! Mean of thl (2nd PDF component)                   [K]
      KK_accr_tndcy,   & ! KK accretion tendency                     [(kg/kg)/s]
      KK_accr_coef,    & ! KK accretion coefficient [(kg/kg)^(1-alpha-beta) / s]
      eta_tol,         & ! Tolerance value of eta                        [kg/kg]
      cthl1,           & ! Coefficient c_thl (1st PDF component)     [(kg/kg)/K]
      cthl2,           & ! Coefficient c_thl (2nd PDF component)     [(kg/kg)/K]
      mixt_frac,       & ! Mixture fraction                                  [-]
      precip_frac_1,   & ! Precipitation fraction (1st PDF component)        [-]
      precip_frac_2      ! Precipitation fraction (2nd PDF component)        [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_thl_KK_accr  ! Covariance of th_l and KK accr. tendency [K(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp,      & ! Exponent on chi                                    [-]
      beta_exp,       & ! Exponent on r_r                                    [-]
      comp_1_contrib, & ! Contribution to thl'KKaccr' (PDF comp. 1) [K(kg/kg)/s]
      comp_2_contrib    ! Contribution to thl'KKaccr' (PDF comp. 2) [K(kg/kg)/s]


    ! Values of the KK exponents.
    alpha_exp = KK_accr_rc_exp
    beta_exp  = KK_accr_rr_exp

    ! Calculate the contribution from PDF component 1 to the covariance of
    ! th_l and KK evaporation tendency.
    comp_1_contrib  &
    = KK_accr_coef * precip_frac_1 &
      * ( ( one / ( two * cthl1 ) )  &
          * ( trivar_NNL_covar_eq( mu_eta_1, mu_chi_1, mu_rr_1, mu_rr_1_n, &
                                   sigma_eta_1, sigma_chi_1, sigma_rr_1, &
                                   sigma_rr_1_n, corr_eta_chi_1, &
                                   corr_eta_rr_1_n, corr_chi_rr_1_n, &
                                   mu_eta_1, KK_accr_tndcy, &
                                   KK_accr_coef, eta_tol, rr_tol, &
                                   alpha_exp, beta_exp )  &
              - bivar_NL_mean_eq( mu_chi_1, mu_rr_1, mu_rr_1_n, sigma_chi_1, &
                                  sigma_rr_1, sigma_rr_1_n, corr_chi_rr_1_n, &
                                  rr_tol, alpha_exp + one, beta_exp )  &
            ) &
          + ( mu_thl_1 - thlm + mu_chi_1 / ( two * cthl1 ) )  &
            * bivar_NL_mean_eq( mu_chi_1, mu_rr_1, mu_rr_1_n, sigma_chi_1, &
                                sigma_rr_1, sigma_rr_1_n, corr_chi_rr_1_n, &
                                rr_tol, alpha_exp, beta_exp )  &
        )

    ! Calculate the contribution from PDF component 2 to the covariance of
    ! th_l and KK evaporation tendency.
    comp_2_contrib  &
    = KK_accr_coef * precip_frac_2 &
      * ( ( one / ( two * cthl2 ) )  &
          * ( trivar_NNL_covar_eq( mu_eta_2, mu_chi_2, mu_rr_2, mu_rr_2_n, &
                                   sigma_eta_2, sigma_chi_2, sigma_rr_2, &
                                   sigma_rr_2_n, corr_eta_chi_2, &
                                   corr_eta_rr_2_n, corr_chi_rr_2_n, &
                                   mu_eta_2, KK_accr_tndcy, &
                                   KK_accr_coef, eta_tol, rr_tol, &
                                   alpha_exp, beta_exp )  &
              - bivar_NL_mean_eq( mu_chi_2, mu_rr_2, mu_rr_2_n, sigma_chi_2, &
                                  sigma_rr_2, sigma_rr_2_n, corr_chi_rr_2_n, &
                                  rr_tol, alpha_exp + one, beta_exp )  &
            ) &
          + ( mu_thl_2 - thlm + mu_chi_2 / ( two * cthl2 ) )  &
            * bivar_NL_mean_eq( mu_chi_2, mu_rr_2, mu_rr_2_n, sigma_chi_2, &
                                sigma_rr_2, sigma_rr_2_n, corr_chi_rr_2_n, &
                                rr_tol, alpha_exp, beta_exp )  &
        )

    ! Calculate the covariance of th_l and KK accretion tendency.
    covar_thl_KK_accr &
    = mixt_frac * comp_1_contrib + ( one - mixt_frac ) * comp_2_contrib


    return

  end function covar_thl_KK_accr

  !=============================================================================
  function quadrivar_NNLL_covar_eq( mu_x_i, mu_chi_i, mu_rr_i, mu_Nr_i, &
                                    mu_rr_i_n, mu_Nr_i_n, sigma_x_i, &
                                    sigma_chi_i, sigma_rr_i, sigma_Nr_i, &
                                    sigma_rr_i_n, sigma_Nr_i_n, &
                                    corr_x_chi_i, corr_x_rr_i_n, &
                                    corr_x_Nr_i_n, corr_chi_rr_i_n, &
                                    corr_chi_Nr_i_n, corr_rr_Nr_i_n, &
                                    x_mean, mc_tndcy_mean, &
                                    mc_coef, x_tol, &
                                    alpha_exp_in, beta_exp_in, gamma_exp_in )

    ! Description:
    ! This function calculates the contribution by the ith PDF component to the
    ! expression < y1'y2'_(i) >, where y1 = x1 ( = x, which is w or t), and
    ! where y2 = x2^alpha x3^beta x4^gamma ( = chi^alpha r_r^beta N_r^gamma,
    ! which also equals KK_evap_tndcy / KK_evap_coef).  The value of covariance
    ! of x and KK evaporation tendency is:
    !
    ! < x'KK_evap' > = KK_evap_coef
    !                  * ( mixt_frac < y1'y2'_(1) >
    !                      + ( 1 - mixt_frac ) < y1'y2'_(2) > ).
    ! 
    ! One of four functions are called, based on whether x1 and/or x2 (x and/or
    ! chi) vary.  Each one of these four functions is the result of an evaluated
    ! integral based on the specific situation.

    ! References:
    ! Section J.1 of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Section S5 of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use PDF_integrals_covars, only: &
        quadrivar_NNLL_covar,            & ! Procedure(s)
        quadrivar_NNLL_covar_const_x1,   &
        quadrivar_NNLL_covar_const_x2,   &
        quadrivar_NNLL_covar_const_x3,   &
        quadrivar_NNLL_covar_const_x1x2, &
        quadrivar_NNLL_covar_const_x1x3, &
        quadrivar_NNLL_covar_const_x2x3, &
        quadrivar_NNLL_covar_const_x3x4, &
        quadrivar_NNLL_covar_cst_x1x2x3, &
        quadrivar_NNLL_covar_cst_x1x3x4, &
        quadrivar_NNLL_covar_cst_x2x3x4, &
        quadrivar_NNLL_covar_const_all

    use constants_clubb, only: &
        chi_tol,             & ! Constant(s)
        rr_tol,              &
        Nr_tol,              &
        parab_cyl_max_input, &
        zero

    use clubb_precision, only: &
        dp,        & ! double precision
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x_i,          & ! Mean of x (ith PDF component)   [units vary (un. v.)]
      mu_chi_i,        & ! Mean of chi (old s) (ith PDF component)       [kg/kg]
      mu_rr_i,         & ! Mean of rr (ith PDF component) ip             [kg/kg]
      mu_Nr_i,         & ! Mean of Nr (ith PDF component) ip            [num/kg]
      mu_rr_i_n,       & ! Mean of ln rr (ith PDF component) ip      [ln(kg/kg)]
      mu_Nr_i_n,       & ! Mean of ln Nr (ith PDF component) ip     [ln(num/kg)]
      sigma_x_i,       & ! Standard deviation of x (ith PDF component)  [un. v.]
      sigma_chi_i,     & ! Standard deviation of chi (ith PDF component) [kg/kg]
      sigma_rr_i,      & ! Standard deviation of rr (ith PDF comp.) ip   [kg/kg]
      sigma_Nr_i,      & ! Standard deviation of Nr (ith PDF comp.) ip  [num/kg]
      sigma_rr_i_n,    & ! Standard deviation of ln rr (ith PDF comp.) ip    [-]
      sigma_Nr_i_n,    & ! Standard deviation of ln Nr (ith PDF comp.) ip    [-]
      corr_x_chi_i,    & ! Correlation of x and chi (ith PDF component)      [-]
      corr_x_rr_i_n,   & ! Correlation of x and ln rr (ith PDF comp.) ip     [-]
      corr_x_Nr_i_n,   & ! Correlation of x and ln Nr (ith PDF comp.) ip     [-]
      corr_chi_rr_i_n, & ! Correlation of chi and ln rr (ith PDF comp.) ip   [-]
      corr_chi_Nr_i_n, & ! Correlation of chi and ln Nr (ith PDF comp.) ip   [-]
      corr_rr_Nr_i_n     ! Correlation of ln rr & ln Nr (ith PDF comp.) ip   [-]

    real( kind = core_rknd ), intent(in) :: &
      x_mean,        & ! Mean of x (overall)                        [units vary]
      mc_tndcy_mean, & ! Mean of microphysics tendency               [(kg/kg)/s]
      mc_coef,       & ! Coefficient of microphysics tendency        [(kg/kg)/s]
      x_tol            ! Tolerance value of x                       [units vary]

    real( kind = core_rknd ), intent(in) :: &
      alpha_exp_in,  & ! Exponent alpha, corresponding to chi                [-]
      beta_exp_in,   & ! Exponent beta, corresponding to rr                  [-]
      gamma_exp_in     ! Exponent gamma, corresponding to Nr                 [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      quadrivar_NNLL_covar_eq

    ! Local Variables
    real( kind = dp ) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3,      & ! Mean of x3 (ith PDF component)                        [-]
      mu_x4,      & ! Mean of x4 (ith PDF component)                        [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      mu_x4_n,    & ! Mean of ln x4 (ith PDF component)                     [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x2,   & ! Standard deviation of x2 (ith PDF component)          [-]
      sigma_x3,   & ! Standard deviation of x3 (ith PDF component)          [-]
      sigma_x4,   & ! Standard deviation of x4 (ith PDF component)          [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      sigma_x4_n, & ! Standard deviation of ln x4 (ith PDF component)       [-]
      rho_x1x2,   & ! Correlation of x1 and x2 (ith PDF component)          [-]
      rho_x1x3_n, & ! Correlation of x1 and ln x3 (ith PDF component)       [-]
      rho_x1x4_n, & ! Correlation of x1 and ln x4 (ith PDF component)       [-]
      rho_x2x3_n, & ! Correlation of x2 and ln x3 (ith PDF component)       [-]
      rho_x2x4_n, & ! Correlation of x2 and ln x4 (ith PDF component)       [-]
      rho_x3x4_n    ! Correlation of ln x3 & ln x4 (ith PDF component)      [-]

    real( kind = dp ) :: &
      x1_mean,                        & ! Mean of x1 (overall)              [-]
      x2_alpha_x3_beta_x4_gamma_mean    ! Mean of x2^alpha x3^beta x4^gamma [-]
    
    real( kind = dp ) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x3                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x4                   [-]

    real( kind = dp ) :: &
      x1_tol, & ! Tolerance value of x1                                     [-]
      x2_tol, & ! Tolerance value of x2                                     [-]
      x3_tol, & ! Tolerance value of x3                                     [-]
      x4_tol, & ! Tolerance value of x4                                     [-]
      s_cc      ! Parabolic cylinder function input value                   [-]


    ! Means for the ith PDF component. 
    mu_x1 = real( mu_x_i, kind = dp )    ! x is w or eta
    mu_x2 = real( mu_chi_i, kind = dp )
    if ( beta_exp_in >= zero ) then
       mu_x3 = real( mu_rr_i, kind = dp )
    else ! exponent beta < 0
       mu_x3 = real( max( mu_rr_i, rr_tol ), kind = dp )
    endif
    if ( gamma_exp_in >= zero ) then
       mu_x4 = real( mu_Nr_i, kind = dp )
    else ! exponent gamma < 0
       mu_x4 = real( max( mu_Nr_i, Nr_tol ), kind = dp )
    endif
    mu_x3_n = real( mu_rr_i_n, kind = dp )
    mu_x4_n = real( mu_Nr_i_n, kind = dp )

    ! Standard deviations for the ith PDF component.
    sigma_x1   = real( sigma_x_i, kind = dp )    ! x is w or eta
    sigma_x2   = real( sigma_chi_i, kind = dp )
    sigma_x3   = real( sigma_rr_i, kind = dp )
    sigma_x4   = real( sigma_Nr_i, kind = dp )
    sigma_x3_n = real( sigma_rr_i_n, kind = dp )
    sigma_x4_n = real( sigma_Nr_i_n, kind = dp )

    ! Correlations for the ith PDF component.
    rho_x1x2   = real( corr_x_chi_i, kind = dp )    ! x is w or eta
    rho_x1x3_n = real( corr_x_rr_i_n, kind = dp )   ! x is w or eta
    rho_x1x4_n = real( corr_x_Nr_i_n, kind = dp )   ! x is w or eta
    rho_x2x3_n = real( corr_chi_rr_i_n, kind = dp )
    rho_x2x4_n = real( corr_chi_Nr_i_n, kind = dp )
    rho_x3x4_n = real( corr_rr_Nr_i_n, kind = dp )


    ! Overall means.
    x1_mean = real( x_mean, kind = dp )  ! x is w or eta
    x2_alpha_x3_beta_x4_gamma_mean = real( mc_tndcy_mean / mc_coef, kind = dp )

    ! Exponents.
    alpha_exp = real( alpha_exp_in, kind = dp )
    beta_exp  = real( beta_exp_in, kind = dp )
    gamma_exp = real( gamma_exp_in, kind = dp )

    ! Tolerance values.
    ! When the standard deviation of a variable is below the tolerance values,
    ! it is considered to be zero, and the variable is considered to have a
    ! constant value.
    x1_tol = real( x_tol, kind = dp )  ! x is w or eta.
    x2_tol = real( chi_tol, kind = dp )
    x3_tol = real( rr_tol, kind = dp )
    x4_tol = real( Nr_tol, kind = dp )

    ! Determine the value of the parabolic cylinder function input value, s_cc.
    ! The value s_cc is being fed into the parabolic cylinder function.  When
    ! the value of s_cc is too large in magnitude (depending on the order of the
    ! parabolic cylinder function), overflow occurs, and the output of the
    ! parabolic cylinder function is +/-Inf.  This is primarily due to a large
    ! ratio of mu_x2 to sigma_x2.  When the value of s_cc is very large, the
    ! distribution of x2 is basically a spike near the mean, so x2 is treated as
    ! a constant.
    if ( sigma_x2 > x2_tol ) then
       s_cc = ( mu_x2 / sigma_x2 )  &
              + rho_x2x3_n * sigma_x3_n * beta_exp  &
              + rho_x2x4_n * sigma_x4_n * gamma_exp
    else  ! sigma_x2 = 0
       ! Note:  s_cc is +inf when mu_x2 > 0 and sigma_x2 = 0, and s_cc is -inf
       !        when mu_x2 < 0 and sigma_x2 = 0.  Furthermore, s_cc is undefined
       !        when mu_x2 = 0 and sigma_x2 = 0.  However, within the context of
       !        this particular function, only the absolute value of s_cc is
       !        relevant, and furthermore the absolute value of s_cc is only
       !        relevant when sigma_x2 > 0.  Therefore, this statement only
       !        serves as divide-by-zero and compiler warning prevention.
       s_cc = huge( s_cc )
    endif


    ! Based on the values of sigma_x1, sigma_x2 (including the value of s_cc
    ! compared to parab_cyl_max_input), sigma_x3, and sigma_x4, find the correct
    ! form of the quadrivariate equation to use.

    if ( sigma_x1 <= x1_tol .and.  &
         ( sigma_x2 <= x2_tol .or.  &
           abs( s_cc ) > real( parab_cyl_max_input, kind = dp ) ) .and.  &
         sigma_x3 <= x3_tol .and. sigma_x4 <= x4_tol ) then

       ! The ith PDF component variances of x (w or eta), chi, r_r (in-precip),
       ! and N_r (in-precip) are all 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_const_all( mu_x1, mu_x2, mu_x3, mu_x4, &
                                         x1_mean, &
                                         x2_alpha_x3_beta_x4_gamma_mean, &
                                         alpha_exp, beta_exp, gamma_exp ),  &
         kind = core_rknd )


    elseif ( sigma_x1 <= x1_tol .and.  &
             ( sigma_x2 <= x2_tol .or.  &
               abs( s_cc ) > real( parab_cyl_max_input, kind = dp ) ) .and.  &
             sigma_x3 <= x3_tol ) then

       ! The ith PDF component variances of x (w or eta), chi, and r_r
       ! (in-precip) are 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_cst_x1x2x3( mu_x1, mu_x2, mu_x3, mu_x4_n, &
                                          sigma_x4_n, x1_mean, &
                                          x2_alpha_x3_beta_x4_gamma_mean, &
                                          alpha_exp, beta_exp, gamma_exp ),  &
         kind = core_rknd )


    elseif ( sigma_x1 <= x1_tol .and.  &
             ( sigma_x2 <= x2_tol .or.  &
               abs( s_cc ) > real( parab_cyl_max_input, kind = dp ) ) .and.  &
             sigma_x4 <= x4_tol ) then

       ! The ith PDF component variances of x (w or eta), chi, and N_r
       ! (in-precip) are 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_cst_x1x2x3( mu_x1, mu_x2, mu_x4, mu_x3_n, &
                                          sigma_x3_n, x1_mean, &
                                          x2_alpha_x3_beta_x4_gamma_mean, &
                                          alpha_exp, gamma_exp, beta_exp ),  &
         kind = core_rknd )


    elseif ( sigma_x1 <= x1_tol .and. sigma_x3 <= x3_tol .and.  &
             sigma_x4 <= x4_tol ) then

       ! The ith PDF component variances of x (w or eta), r_r (in-precip), and
       ! N_r (in-precip) are 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_cst_x1x3x4( mu_x1, mu_x2, mu_x3, mu_x4, &
                                          sigma_x2, x1_mean, &
                                          x2_alpha_x3_beta_x4_gamma_mean, &
                                          alpha_exp, beta_exp, gamma_exp ),  &
         kind = core_rknd )


    elseif ( ( sigma_x2 <= x2_tol .or.  &
               abs( s_cc ) > real( parab_cyl_max_input, kind = dp ) ) .and.  &
             sigma_x3 <= x3_tol .and. sigma_x4 <= x4_tol ) then

       ! The ith PDF component variances of chi, r_r (in-precip), and N_r
       ! (in-precip) are 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_cst_x2x3x4( mu_x1, mu_x2, mu_x3, mu_x4, &
                                          x1_mean, &
                                          x2_alpha_x3_beta_x4_gamma_mean, &
                                          alpha_exp, beta_exp, gamma_exp ),  &
         kind = core_rknd )


    elseif ( sigma_x1 <= x1_tol .and.  &
             ( sigma_x2 <= x2_tol .or.  &
               abs( s_cc ) > real( parab_cyl_max_input, kind = dp ) ) ) then

       ! The ith PDF component variances of both x (w or eta) and chi are 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_const_x1x2( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                          sigma_x3_n, sigma_x4_n, &
                                          rho_x3x4_n, x1_mean, &
                                          x2_alpha_x3_beta_x4_gamma_mean, &
                                          alpha_exp, beta_exp, gamma_exp ),  &
         kind = core_rknd )


    elseif ( sigma_x1 <= x1_tol .and. sigma_x3 <= x3_tol ) then

       ! The ith PDF component variances of both x (w or eta) and r_r
       ! (in-precip) are 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_const_x1x3( mu_x1, mu_x2, mu_x3, mu_x4_n, &
                                          sigma_x2, sigma_x4_n, rho_x2x4_n, &
                                          x1_mean, &
                                          x2_alpha_x3_beta_x4_gamma_mean, &
                                          alpha_exp, beta_exp, gamma_exp ),  &
         kind = core_rknd )


    elseif ( sigma_x1 <= x1_tol .and. sigma_x4 <= x4_tol ) then

       ! The ith PDF component variances of both x (w or eta) and N_r
       ! (in-precip) are 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_const_x1x3( mu_x1, mu_x2, mu_x4, mu_x3_n, &
                                          sigma_x2, sigma_x3_n, rho_x2x3_n, &
                                          x1_mean, &
                                          x2_alpha_x3_beta_x4_gamma_mean, &
                                          alpha_exp, gamma_exp, beta_exp ),  &
         kind = core_rknd )


    elseif ( ( sigma_x2 <= x2_tol .or.  &
               abs( s_cc ) > real( parab_cyl_max_input, kind = dp ) ) .and.  &
             sigma_x3 <= x3_tol ) then

       ! The ith PDF component variances of both chi and r_r (in-precip) are 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_const_x2x3( mu_x1, mu_x2, mu_x3, mu_x4_n, &
                                          sigma_x1, sigma_x4_n, rho_x1x4_n, &
                                          x1_mean, &
                                          x2_alpha_x3_beta_x4_gamma_mean, &
                                          alpha_exp, beta_exp, gamma_exp ),  &
         kind = core_rknd )


    elseif ( ( sigma_x2 <= x2_tol .or.  &
               abs( s_cc ) > real( parab_cyl_max_input, kind = dp ) ) .and.  &
             sigma_x4 <= x4_tol ) then

       ! The ith PDF component variances of both chi and N_r (in-precip) are 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_const_x2x3( mu_x1, mu_x2, mu_x4, mu_x3_n, &
                                          sigma_x1, sigma_x3_n, rho_x1x3_n, &
                                          x1_mean, &
                                          x2_alpha_x3_beta_x4_gamma_mean, &
                                          alpha_exp, gamma_exp, beta_exp ),  &
         kind = core_rknd )


    elseif ( sigma_x3 <= x3_tol .and. sigma_x4 <= x4_tol ) then

       ! The ith PDF component variances of both r_r (in-precip) and N_r
       ! (in-precip) are 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_const_x3x4( mu_x1, mu_x2, mu_x3, mu_x4, &
                                          sigma_x1, sigma_x2, rho_x1x2, &
                                          x1_mean, &
                                          x2_alpha_x3_beta_x4_gamma_mean, &
                                          alpha_exp, beta_exp, gamma_exp ),  &
         kind = core_rknd )


    elseif ( sigma_x1 <= x1_tol ) then

       ! The ith PDF component variance of x (w or eta) is 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_const_x1( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                        sigma_x2, sigma_x3_n, sigma_x4_n, &
                                        rho_x2x3_n, rho_x2x4_n, rho_x3x4_n, &
                                        x1_mean, &
                                        x2_alpha_x3_beta_x4_gamma_mean, &
                                        alpha_exp, beta_exp, gamma_exp ),  &
         kind = core_rknd )


    elseif ( sigma_x2 <= x2_tol .or.  &
             abs( s_cc ) > real( parab_cyl_max_input, kind = dp ) ) then

       ! The ith PDF component variance of chi is 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_const_x2( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                        sigma_x1, sigma_x3_n, sigma_x4_n, &
                                        rho_x1x3_n, rho_x1x4_n, rho_x3x4_n, &
                                        x1_mean, &
                                        x2_alpha_x3_beta_x4_gamma_mean, &
                                        alpha_exp, beta_exp, gamma_exp ),  &
         kind = core_rknd )


    elseif ( sigma_x3 <= x3_tol ) then

       ! The ith PDF component variance of r_r (in-precip) is 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_const_x3( mu_x1, mu_x2, mu_x3, mu_x4_n, &
                                        sigma_x1, sigma_x2, sigma_x4_n, &
                                        rho_x1x2, rho_x1x4_n, rho_x2x4_n, &
                                        x1_mean, &
                                        x2_alpha_x3_beta_x4_gamma_mean, &
                                        alpha_exp, beta_exp, gamma_exp ),  &
         kind = core_rknd )


    elseif ( sigma_x4 <= x4_tol ) then

       ! The ith PDF component variance of N_r (in-precip) is 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_const_x3( mu_x1, mu_x2, mu_x4, mu_x3_n, &
                                        sigma_x1, sigma_x2, sigma_x3_n, &
                                        rho_x1x2, rho_x1x3_n, rho_x2x3_n, &
                                        x1_mean, &
                                        x2_alpha_x3_beta_x4_gamma_mean, &
                                        alpha_exp, gamma_exp, beta_exp ),  &
         kind = core_rknd )


    else  ! sigma_x1, sigma_x2, sigma_x3, and sigma_x4 > 0.

       ! This is the complete value of the quadrivariate.
       ! All fields vary in the ith PDF component.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                               sigma_x1, sigma_x2, sigma_x3_n, sigma_x4_n, &
                               rho_x1x2, rho_x1x3_n, rho_x1x4_n, &
                               rho_x2x3_n, rho_x2x4_n, rho_x3x4_n, &
                               x1_mean, x2_alpha_x3_beta_x4_gamma_mean, &
                               alpha_exp, beta_exp, gamma_exp ),  &
         kind = core_rknd )


    endif


    return

  end function quadrivar_NNLL_covar_eq

  !=============================================================================
  function trivar_NNL_covar_eq( mu_x_i, mu_chi_i, mu_y_i, mu_y_i_n, &
                                sigma_x_i, sigma_chi_i, sigma_y_i, &
                                sigma_y_i_n, corr_x_chi_i, &
                                corr_x_y_i_n, corr_chi_y_i_n, &
                                x_mean, mc_tndcy_mean, &
                                mc_coef, x_tol, y_tol, &
                                alpha_exp_in, beta_exp_in )

    ! Description:
    ! This function calculates the contribution by the ith PDF component to the
    ! expression < y1'y2'_(i) >, where y1 = x1 ( = x, which is w or t), and
    ! where y2 = x2^alpha x3^beta ( = chi^alpha y^beta, where y is N_cn or r_r
    ! for autoconversion or accretion, respectively, and which also equals
    ! KK_auto_tndcy / KK_auto_coef or KK_accr_tndcy / KK_accr_coef,
    ! respectively).  The value of covariance of x and the KK microphysics
    ! tendency is:
    !
    ! < x'KK_mc' > = KK_mc_coef
    !                * ( mixt_frac < y1'y2'_(1) >
    !                    + ( 1 - mixt_frac ) < y1'y2'_(2) > ).
    ! 
    ! One of four functions are called, based on whether x1 and/or x2 (x and/or
    ! chi) vary.  Each one of these four functions is the result of an evaluated
    ! integral based on the specific situation.

    ! References:
    ! Section J.2 of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Section S6 of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use PDF_integrals_covars, only: &
        trivar_NNL_covar,            & ! Procedure(s)
        trivar_NNL_covar_const_x1,   &
        trivar_NNL_covar_const_x2,   &
        trivar_NNL_covar_const_x3,   &
        trivar_NNL_covar_const_x1x2, &
        trivar_NNL_covar_const_x1x3, &
        trivar_NNL_covar_const_x2x3, &
        trivar_NNL_covar_const_all

    use constants_clubb, only: &
        chi_tol,             & ! Constant(s)
        parab_cyl_max_input, &
        zero

    use clubb_precision, only: &
        dp,        & ! double precision
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x_i,         & ! Mean of x (ith PDF component)    [units vary (un. v.)]
      mu_chi_i,       & ! Mean of chi (old s) (ith PDF component)        [kg/kg]
      mu_y_i,         & ! Mean of y (ith PDF component)             [units vary]
      mu_y_i_n,       & ! Mean of ln y (ith PDF component)     [ln (units vary)]
      sigma_x_i,      & ! Standard deviation of x (ith PDF component)   [un. v.]
      sigma_chi_i,    & ! Standard deviation of chi (ith PDF component)  [kg/kg]
      sigma_y_i,      & ! Standard deviation of y (ith PDF component)   [un. v.]
      sigma_y_i_n,    & ! Standard deviation of ln y (ith PDF component)     [-]
      corr_x_chi_i,   & ! Correlation of x and chi (ith PDF component)       [-]
      corr_x_y_i_n,   & ! Correlation of x and ln y (ith PDF component)      [-]
      corr_chi_y_i_n    ! Correlation of chi and ln y (ith PDF component)    [-]

    real( kind = core_rknd ), intent(in) :: &
      x_mean,        & ! Mean of x (overall)                        [units vary]
      mc_tndcy_mean, & ! Mean of microphysics tendency               [(kg/kg)/s]
      mc_coef,       & ! Coefficient of microphysics                 [(kg/kg)/s]
      x_tol,         & ! Tolerance value of x                       [units vary]
      y_tol            ! Tolerance value of y                       [units vary]

    real( kind = core_rknd ), intent(in) :: &
      alpha_exp_in,  & ! Exponent alpha, corresponding to chi                [-]
      beta_exp_in      ! Exponent beta, corresponding to y                   [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      trivar_NNL_covar_eq

    ! Local Variables
    real( kind = dp ) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3,      & ! Mean of x3 (ith PDF component)                        [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x2,   & ! Standard deviation of x2 (ith PDF component)          [-]
      sigma_x3,   & ! Standard deviation of x3 (ith PDF component)          [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      rho_x1x2,   & ! Correlation of x1 and x2 (ith PDF component)          [-]
      rho_x1x3_n, & ! Correlation of x1 and ln x3 (ith PDF component)       [-]
      rho_x2x3_n    ! Correlation of x2 and ln x3 (ith PDF component)       [-]

    real( kind = dp ) :: &
      x1_mean,               & ! Mean of x1 (overall)                       [-]
      x2_alpha_x3_beta_mean    ! Mean of x2^alpha x3^beta                   [-]
    
    real( kind = dp ) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp      ! Exponent beta, corresponding to x3                    [-]

    real( kind = dp ) :: &
      x1_tol, & ! Tolerance value of x1                                     [-]
      x2_tol, & ! Tolerance value of x2                                     [-]
      x3_tol, & ! Tolerance value of x3                                     [-]
      s_c       ! Parabolic cylinder function input value                   [-]


    ! Means for the ith PDF component. 
    mu_x1 = real( mu_x_i, kind = dp )   ! x is w or eta
    mu_x2 = real( mu_chi_i, kind = dp )
    if ( beta_exp_in >= zero ) then
       mu_x3 = real( mu_y_i, kind = dp )   ! y is N_cn or r_r
    else ! exponent beta < 0
       mu_x3 = real( max( mu_y_i, y_tol ), kind = dp )  ! y is N_cn or r_r
    endif
    mu_x3_n = real( mu_y_i_n, kind = dp ) ! y is N_cn or r_r

    ! Standard deviations for the ith PDF component.
    sigma_x1   = real( sigma_x_i, kind = dp )   ! x is w or eta
    sigma_x2   = real( sigma_chi_i, kind = dp )
    sigma_x3   = real( sigma_y_i, kind = dp )   ! y is N_cn or r_r
    sigma_x3_n = real( sigma_y_i_n, kind = dp ) ! y is N_cn or r_r

    ! Correlations for the ith PDF component.
    rho_x1x2   = real( corr_x_chi_i, kind = dp )   ! x is w or eta
    rho_x1x3_n = real( corr_x_y_i_n, kind = dp )   ! x: w or eta; y: N_cn or r_r
    rho_x2x3_n = real( corr_chi_y_i_n, kind = dp ) ! y is N_cn or r_r

    ! Overall means.
    x1_mean = real( x_mean, kind = dp )  ! x is w or eta
    x2_alpha_x3_beta_mean = real( mc_tndcy_mean / mc_coef, kind = dp )

    ! Exponents.
    alpha_exp = real( alpha_exp_in, kind = dp )
    beta_exp  = real( beta_exp_in, kind = dp )

    ! Tolerance values.
    ! When the standard deviation of a variable is below the tolerance values,
    ! it is considered to be zero, and the variable is considered to have a
    ! constant value.
    x1_tol = real( x_tol, kind = dp )  ! x is w or eta
    x2_tol = real( chi_tol, kind = dp )
    x3_tol = real( y_tol, kind = dp )  ! y is N_cn or r_r

    ! Determine the value of the parabolic cylinder function input value, s_c.
    ! The value s_c is being fed into the parabolic cylinder function.  When
    ! the value of s_c is too large in magnitude (depending on the order of the
    ! parabolic cylinder function), overflow occurs, and the output of the
    ! parabolic cylinder function is +/-Inf.  This is primarily due to a large
    ! ratio of mu_x2 to sigma_x2.  When the value of s_c is very large, the
    ! distribution of x2 is basically a spike near the mean, so x2 is treated as
    ! a constant.
    if ( sigma_x2 > x2_tol ) then
       s_c = ( mu_x2 / sigma_x2 )  + rho_x2x3_n * sigma_x3_n * beta_exp
    else  ! sigma_x2 = 0
       ! Note:  s_c is +inf when mu_x2 > 0 and sigma_x2 = 0, and s_c is -inf
       !        when mu_x2 < 0 and sigma_x2 = 0.  Furthermore, s_c is undefined
       !        when mu_x2 = 0 and sigma_x2 = 0.  However, within the context of
       !        this particular function, only the absolute value of s_c is
       !        relevant, and furthermore the absolute value of s_c is only
       !        relevant when sigma_x2 > 0.  Therefore, this statement only
       !        serves as divide-by-zero and compiler warning prevention.
       s_c = huge( s_c )
    endif


    ! Based on the values of sigma_x1, sigma_x2 (including the value of s_c
    ! compared to parab_cyl_max_input), and sigma_x3, find the correct form of
    ! the trivariate equation to use.

    if ( sigma_x1 <= x1_tol .and.  &
         ( sigma_x2 <= x2_tol .or.  &
           abs( s_c ) > real( parab_cyl_max_input, kind = dp ) ) .and.  &
         sigma_x3 <= x3_tol ) then

       ! The ith PDF component variances of x (w or eta), chi, and y (N_cn or
       ! r_r (in-precip)) are all 0.
       trivar_NNL_covar_eq  &
       = real( trivar_NNL_covar_const_all( mu_x1, mu_x2, mu_x3, &
                                           x1_mean, x2_alpha_x3_beta_mean, &
                                           alpha_exp, beta_exp ),  &
               kind = core_rknd )


    elseif ( sigma_x1 <= x1_tol .and.  &
             ( sigma_x2 <= x2_tol .or.  &
               abs( s_c ) > real( parab_cyl_max_input, kind = dp ) ) ) then

       ! The ith PDF component variances of both x (w or eta) and chi are 0.
       trivar_NNL_covar_eq  &
       = real( trivar_NNL_covar_const_x1x2( mu_x1, mu_x2, mu_x3_n, sigma_x3_n, &
                                            x1_mean, x2_alpha_x3_beta_mean, &
                                            alpha_exp, beta_exp ),  &
               kind = core_rknd )


    elseif ( sigma_x1 <= x1_tol .and. sigma_x3 <= x3_tol ) then

       ! The ith PDF component variances of both x (w or eta) and y (N_cn or
       ! r_r (in-precip)) are 0.
       trivar_NNL_covar_eq  &
       = real( trivar_NNL_covar_const_x1x3( mu_x1, mu_x2, mu_x3, sigma_x2, &
                                            x1_mean, x2_alpha_x3_beta_mean, &
                                            alpha_exp, beta_exp ),  &
               kind = core_rknd )


    elseif ( ( sigma_x2 <= x2_tol .or.  &
               abs( s_c ) > real( parab_cyl_max_input, kind = dp ) ) .and.  &
             sigma_x3 <= x3_tol ) then

       ! The ith PDF component variances of both chi and y (N_cn or
       ! r_r (in-precip)) are 0.
       trivar_NNL_covar_eq  &
       = real( trivar_NNL_covar_const_x2x3( mu_x1, mu_x2, mu_x3, &
                                            x1_mean, x2_alpha_x3_beta_mean, &
                                            alpha_exp, beta_exp ),  &
               kind = core_rknd )


    elseif ( sigma_x1 <= x1_tol ) then

       ! The ith PDF component variance of x (w or eta) is 0.
       trivar_NNL_covar_eq  &
       = real( trivar_NNL_covar_const_x1( mu_x1, mu_x2, mu_x3_n, &
                                          sigma_x2, sigma_x3_n, rho_x2x3_n, &
                                          x1_mean, x2_alpha_x3_beta_mean, &
                                          alpha_exp, beta_exp ),  &
               kind = core_rknd )


    elseif ( sigma_x2 <= x2_tol .or.  &
             abs( s_c ) > real( parab_cyl_max_input, kind = dp ) ) then

       ! The ith PDF component variance of chi is 0.
       trivar_NNL_covar_eq  &
       = real( trivar_NNL_covar_const_x2( mu_x1, mu_x2, mu_x3_n, &
                                          sigma_x1, sigma_x3_n, rho_x1x3_n, &
                                          x1_mean, x2_alpha_x3_beta_mean, &
                                          alpha_exp, beta_exp ),  &
               kind = core_rknd )


    elseif ( sigma_x3 <= x3_tol ) then

       ! The ith PDF component variance of y (N_cn or r_r (in-precip)) is 0.
       trivar_NNL_covar_eq  &
       = real( trivar_NNL_covar_const_x3( mu_x1, mu_x2, mu_x3, &
                                          sigma_x1, sigma_x2, rho_x1x2, &
                                          x1_mean, x2_alpha_x3_beta_mean, &
                                          alpha_exp, beta_exp ),  &
               kind = core_rknd )


    else  ! sigma_x1, sigma_x2, and sigma_x3 > 0.

       ! This is the complete value of the trivariate.
       ! All fields vary in the ith PDF component.
       trivar_NNL_covar_eq  &
       = real( trivar_NNL_covar( mu_x1, mu_x2, mu_x3_n, &
                                 sigma_x1, sigma_x2, sigma_x3_n, &
                                 rho_x1x2, rho_x1x3_n, rho_x2x3_n, &
                                 x1_mean, x2_alpha_x3_beta_mean, &
                                 alpha_exp, beta_exp ),  &
               kind = core_rknd )


    endif


    return

  end function trivar_NNL_covar_eq

!===============================================================================

end module KK_upscaled_covariances
