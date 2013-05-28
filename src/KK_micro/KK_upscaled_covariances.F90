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
             trivar_NNL_covar_eq,     &
             trivar_NNL_covar_eq_Nc0

  contains

  !=============================================================================
  subroutine KK_upscaled_covar_driver( w_mean, exner, rcm, &
                                       rrainm, Nrm, Ncnm, &
                                       mu_s_1, mu_s_2, mu_rr_1, mu_rr_2, &
                                       mu_Nr_1, mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &
                                       mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                       mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &
                                       sigma_s_1, sigma_s_2, sigma_rr_1, &
                                       sigma_rr_2, sigma_Nr_1, sigma_Nr_2, &
                                       sigma_Ncn_1, sigma_Ncn_2, sigma_rr_1_n, &
                                       sigma_rr_2_n, sigma_Nr_1_n, &
                                       sigma_Nr_2_n, sigma_Ncn_1_n, &
                                       sigma_Ncn_2_n, corr_srr_1_n, &
                                       corr_srr_2_n, corr_sNr_1_n, &
                                       corr_sNr_2_n, corr_sNcn_1_n, &
                                       corr_sNcn_2_n, corr_rrNr_1_n, &
                                       corr_rrNr_2_n, mixt_frac, &
                                       precip_frac_1, precip_frac_2, &
                                       Nc0_in_cloud, l_const_Nc_in_cloud, &
                                       KK_evap_coef, KK_auto_coef, &
                                       KK_accr_coef, KK_evap_tndcy, &
                                       KK_auto_tndcy, KK_accr_tndcy, &
                                       pdf_params, level, &
                                       l_stats_samp, &
                                       wprtp_mc_src_tndcy, &
                                       wpthlp_mc_src_tndcy, &
                                       rtp2_mc_src_tndcy, &
                                       thlp2_mc_src_tndcy, &
                                       rtpthlp_mc_src_tndcy )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        two,    & ! Constant(s)
        zero,   &
        Lv,     &
        Cp,     &
        w_tol,  &
        rc_tol, &
        rr_tol, & 
        Nr_tol, & 
        Ncn_tol

    use constants_clubb, only:  &
        t_tol => t_mellor_tol  ! Constant

    use KK_utilities, only: &
        corr_NL2NN  ! Procedure(s)

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s) type

    use KK_fixed_correlations, only: &
        corr_wrr_NL_cloud,  & ! Variable(s)
        corr_wNr_NL_cloud,  &
        corr_wNcn_NL_cloud, &
        corr_wrr_NL_below,  &
        corr_wNr_NL_below,  &
        corr_wNcn_NL_below, &
        corr_trr_NL_cloud,  &
        corr_tNr_NL_cloud,  &
        corr_tNcn_NL_cloud, &
        corr_trr_NL_below,  &
        corr_tNr_NL_below,  &
        corr_tNcn_NL_below, &
        corr_st_NN_cloud,   &
        corr_st_NN_below

    use parameters_microphys, only: &
        l_fix_s_t_correlations ! Variable(s)

    use clubb_precision, only: &
        core_rknd,      & ! Variable(s)
        time_precision

    use stats_type, only: & 
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only: &
        zt,                    & ! Variable(s)
        iw_KK_evap_covar_zt,   &
        irt_KK_evap_covar_zt,  &
        ithl_KK_evap_covar_zt, &
        iw_KK_auto_covar_zt,   &
        irt_KK_auto_covar_zt,  &
        ithl_KK_auto_covar_zt, &
        iw_KK_accr_covar_zt,   &
        irt_KK_accr_covar_zt,  &
        ithl_KK_accr_covar_zt

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      w_mean, & ! Mean vertical velocity                        [m/s]
      exner,  & ! Exner function                                [-]
      rcm,    & ! Mean cloud water mixing ratio (overall)       [kg/kg]
      rrainm, & ! Mean rain water mixing ratio (overall)        [kg/kg]
      Nrm,    & ! Mean rain drop concentration (overall)        [num/kg]
      Ncnm      ! Mean cloud nuclei concentration               [num/kg]

    real( kind = core_rknd ), intent(in) :: &
      mu_s_1,        & ! Mean of s (1st PDF component)                   [kg/kg]
      mu_s_2,        & ! Mean of s (2nd PDF component)                   [kg/kg]
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
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)     [kg/kg]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)     [kg/kg]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_Nr_1,    & ! Standard deviation of Nr (1st PDF component) ip  [#/kg]
      sigma_Nr_2,    & ! Standard deviation of Nr (2nd PDF component) ip  [#/kg]
      sigma_Ncn_1,   & ! Standard deviation of Ncn (1st PDF component)  [num/kg]
      sigma_Ncn_2,   & ! Standard deviation of Ncn (2nd PDF component)  [num/kg]
      sigma_rr_1_n,  & ! Standard dev. of ln rr (1st PDF comp.) ip   [ln(kg/kg)]
      sigma_rr_2_n,  & ! Standard dev. of ln rr (2nd PDF comp.) ip   [ln(kg/kg)]
      sigma_Nr_1_n,  & ! Standard dev. of ln Nr (1st PDF comp.) ip  [ln(num/kg)]
      sigma_Nr_2_n,  & ! Standard dev. of ln Nr (2nd PDF comp.) ip  [ln(num/kg)]
      sigma_Ncn_1_n, & ! Standard dev. of ln Ncn (1st PDF comp.)    [ln(num/kg)]
      sigma_Ncn_2_n, & ! Standard dev. of ln Ncn (2nd PDF comp.)    [ln(num/kg)]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF comp.) ip  [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF comp.) ip  [-]
      corr_sNcn_1_n, & ! Correlation between s and ln Ncn (1st PDF comp.)    [-]
      corr_sNcn_2_n, & ! Correlation between s and ln Ncn (2nd PDF comp.)    [-]
      corr_rrNr_1_n, & ! Correlation btwn. ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rrNr_2_n, & ! Correlation btwn. ln rr & ln Nr (2nd PDF comp.) ip  [-]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2, & ! Precipitation fraction (2nd PDF component)          [-]
      Nc0_in_cloud     ! Constant in-cloud value of cloud droplet conc. [num/kg]

    logical, intent(in) :: &
      l_const_Nc_in_cloud  ! Flag to use a constant value of N_c within cloud

    real( kind = core_rknd ), intent(in) :: &
      KK_evap_coef, & ! KK evaporation coefficient          [(kg/kg)/s]
      KK_auto_coef, & ! KK autoconversion coefficient       [(kg/kg)/s]
      KK_accr_coef    ! KK accretion coefficient            [(kg/kg)/s]

    real( kind = core_rknd ), intent(in) :: &
      KK_evap_tndcy, & ! KK evaporation tendency            [(kg/kg)/s]
      KK_auto_tndcy, & ! KK autoconversion tendency         [(kg/kg)/s]
      KK_accr_tndcy    ! KK accretion tendency              [(kg/kg)/s]

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                        [units vary]

    integer, intent(in) :: &
      level         ! Vertical level index                  [-]

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      wprtp_mc_src_tndcy,   & ! Microphysics tendency for w'rt'  [m*(kg/kg)/s^2]
      wpthlp_mc_src_tndcy,  & ! Microphysics tendency for w'thl' [m*K/s^2]
      rtp2_mc_src_tndcy,    & ! Microphysics tendency for rt'^2  [(kg/kg)^2/s]
      thlp2_mc_src_tndcy,   & ! Microphysics tendency for thl'^2 [K^2/s]
      rtpthlp_mc_src_tndcy    ! Microphysics tend. for rt'thl'   [K*(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      mu_w_1,      & ! Mean of w (1st PDF component)                       [m/s]
      mu_w_2,      & ! Mean of w (2nd PDF component)                       [m/s]
      sigma_w_1,   & ! Standard deviation of w (1st PDF component)         [m/s]
      sigma_w_2,   & ! Standard deviation of w (2nd PDF component)         [m/s]
      sigma_t_1,   & ! Standard deviation of t (1st PDF component)       [kg/kg]
      sigma_t_2,   & ! Standard deviation of t (2nd PDF component)       [kg/kg]
      corr_wrr_1,  & ! Correlation between w and rr (1st PDF component) ip   [-]
      corr_wrr_2,  & ! Correlation between w and rr (2nd PDF component) ip   [-]
      corr_wNr_1,  & ! Correlation between w and Nr (1st PDF component) ip   [-]
      corr_wNr_2,  & ! Correlation between w and Nr (2nd PDF component) ip   [-]
      corr_wNcn_1, & ! Correlation between w and Ncn (1st PDF component)     [-]
      corr_wNcn_2, & ! Correlation between w and Ncn (2nd PDF component)     [-]
      corr_ts_1,   & ! Correlation between t and s (1st PDF component)       [-]
      corr_ts_2,   & ! Correlation between t and s (2nd PDF component)       [-]
      corr_trr_1,  & ! Correlation between t and rr (1st PDF component) ip   [-]
      corr_trr_2,  & ! Correlation between t and rr (2nd PDF component) ip   [-]
      corr_tNr_1,  & ! Correlation between t and Nr (1st PDF component) ip   [-]
      corr_tNr_2,  & ! Correlation between t and Nr (2nd PDF component) ip   [-]
      corr_tNcn_1, & ! Correlation between t and Ncn (1st PDF component)     [-]
      corr_tNcn_2    ! Correlation between t and Ncn (2nd PDF component)     [-]

    real( kind = core_rknd ) :: &
      corr_wrr_1_n,  & ! Correlation between w and ln rr (1st PDF comp.) ip  [-]
      corr_wrr_2_n,  & ! Correlation between w and ln rr (2nd PDF comp.) ip  [-]
      corr_wNr_1_n,  & ! Correlation between w and ln Nr (1st PDF comp.) ip  [-]
      corr_wNr_2_n,  & ! Correlation between w and ln Nr (2nd PDF comp.) ip  [-]
      corr_wNcn_1_n, & ! Correlation between w and ln Ncn (1st PDF comp.)    [-]
      corr_wNcn_2_n, & ! Correlation between w and ln Ncn (2nd PDF comp.)    [-]
      corr_trr_1_n,  & ! Correlation between t and ln rr (1st PDF comp.) ip  [-]
      corr_trr_2_n,  & ! Correlation between t and ln rr (2nd PDF comp.) ip  [-]
      corr_tNr_1_n,  & ! Correlation between t and ln Nr (1st PDF comp.) ip  [-]
      corr_tNr_2_n,  & ! Correlation between t and ln Nr (2nd PDF comp.) ip  [-]
      corr_tNcn_1_n, & ! Correlation between t and ln Ncn (1st PDF comp.)    [-]
      corr_tNcn_2_n    ! Correlation between t and ln Ncn (2nd PDF comp.)    [-]

    real( kind = core_rknd ) :: &
      crt1,  & ! Coefficient c_rt (1st PDF component)    [-]
      crt2,  & ! Coefficient c_rt (2nd PDF component)    [-]
      cthl1, & ! Coefficient c_thl (1st PDF component)   [(kg/kg)/K]
      cthl2    ! Coefficient c_thl (2nd PDF component)   [(kg/kg)/K]

    real( kind = core_rknd ) :: &
      w_KK_evap_covar,   & ! Covar. btw. w and KK evap. tend.    [m*(kg/kg)/s^2]
      rt_KK_evap_covar,  & ! Covar. btw. rt and KK evap. tend.   [(kg/kg)^2/s]
      thl_KK_evap_covar, & ! Covar. btw. thl and KK evap. tend.  [K*(kg/kg)/s]
      w_KK_auto_covar,   & ! Covar. btw. w and KK auto. tend.    [m*(kg/kg)/s^2]
      rt_KK_auto_covar,  & ! Covar. btw. rt and KK auto. tend.   [(kg/kg)^2/s]
      thl_KK_auto_covar, & ! Covar. btw. thl and KK auto. tend.  [K*(kg/kg)/s]
      w_KK_accr_covar,   & ! Covar. btw. w and KK accr. tend.    [m*(kg/kg)/s^2]
      rt_KK_accr_covar,  & ! Covar. btw. rt and KK accr. tend.   [(kg/kg)^2/s]
      thl_KK_accr_covar    ! Covar. btw. thl and KK accr. tend.  [K*(kg/kg)/s]

    ! Constant Parameters
    !
    ! Set the correlations between vertical velocity (w) and extended liquid
    ! water mixing ratio (s_*) to 0.
    ! Note:  The component correlations of w and r_t and the component
    !        correlations of w and theta_l are both set to be 0 within the CLUBB
    !        model code.  In other words, w and r_t (theta_l) have overall
    !        covariance w'r_t' (w'theta_l'), but the single component covariance
    !        and correlation are defined to be 0.  Likewise, the single
    !        component correlation and covariance of w and s, as well as
    !        w and t, are defined to be 0.
    real( kind = core_rknd ), parameter :: &
      corr_ws_1 = zero, & ! Correlation between w and s (1st PDF component) [-]
      corr_ws_2 = zero    ! Correlation between w and s (2nd PDF component) [-]
    !
    ! Set the component mean values of t_* to 0.
    ! Note:  The component mean values of t_* are not important.  They can be
    !        set to anything.  They cancel out in the model code.  However, the
    !        best thing to do is to set them to 0 and avoid any kind of
    !        numerical error.
    real( kind = core_rknd ), parameter :: &
      mu_t_1 = zero, & ! Mean of t (1st PDF component)  [kg/kg]
      mu_t_2 = zero    ! Mean of t (2nd PDF component)  [kg/kg]


    ! Enter the PDF parameters.
    mu_w_1    = pdf_params%w1
    mu_w_2    = pdf_params%w2
    sigma_w_1 = sqrt( pdf_params%varnce_w1 )
    sigma_w_2 = sqrt( pdf_params%varnce_w2 )
    sigma_t_1 = pdf_params%stdev_t1
    sigma_t_2 = pdf_params%stdev_t2

    if ( l_fix_s_t_correlations ) then
       if ( mu_s_1 > zero ) then
          corr_ts_1 = corr_st_NN_cloud
       else
          corr_ts_1 = corr_st_NN_below
       endif
       if ( mu_s_2 > zero ) then
          corr_ts_2 = corr_st_NN_cloud
       else
          corr_ts_2 = corr_st_NN_below
       endif
    else
       corr_ts_1 = pdf_params%corr_st_1
       corr_ts_2 = pdf_params%corr_st_2
    endif

    crt1      = pdf_params%crt1
    crt2      = pdf_params%crt2
    cthl1     = pdf_params%cthl1
    cthl2     = pdf_params%cthl2

    
    ! Set up the values of the statistical correlations and variances.  Since we
    ! currently do not have enough variables to compute the correlations and
    ! variances directly, we have obtained these values by analyzing LES runs of
    ! certain cases.  We have divided those results into an inside-cloud average
    ! and an outside-cloud (or below-cloud) average.  This coding leaves the
    ! software architecture in place in case we ever have the variables in place
    ! to compute these values directly.  It also allows us to use separate
    ! inside-cloud and outside-cloud parameter values.
    !
    ! Set the value of the parameters based on whether the altitude is above or
    ! below cloud base.  Determine whether there is cloud at any given vertical
    ! level.  In order for a vertical level to have cloud, the amount of cloud
    ! water (rcm) must be greater than or equal to the tolerance level (rc_tol).
    ! If there is cloud at a given vertical level, then the ###_cloud value is
    ! used.  Otherwise, the ###_below value is used.
    if ( rcm > rc_tol ) then
       corr_wrr_1  = corr_wrr_NL_cloud
       corr_wrr_2  = corr_wrr_NL_cloud
       corr_wNr_1  = corr_wNr_NL_cloud
       corr_wNr_2  = corr_wNr_NL_cloud
       corr_wNcn_1 = corr_wNcn_NL_cloud
       corr_wNcn_2 = corr_wNcn_NL_cloud
       corr_trr_1  = corr_trr_NL_cloud
       corr_trr_2  = corr_trr_NL_cloud
       corr_tNr_1  = corr_tNr_NL_cloud
       corr_tNr_2  = corr_tNr_NL_cloud
       corr_tNcn_1 = corr_tNcn_NL_cloud
       corr_tNcn_2 = corr_tNcn_NL_cloud
    else
       corr_wrr_1  = corr_wrr_NL_below
       corr_wrr_2  = corr_wrr_NL_below
       corr_wNr_1  = corr_wNr_NL_below
       corr_wNr_2  = corr_wNr_NL_below
       corr_wNcn_1 = corr_wNcn_NL_below
       corr_wNcn_2 = corr_wNcn_NL_below
       corr_trr_1  = corr_trr_NL_below
       corr_trr_2  = corr_trr_NL_below
       corr_tNr_1  = corr_tNr_NL_below
       corr_tNr_2  = corr_tNr_NL_below
       corr_tNcn_1 = corr_tNcn_NL_below
       corr_tNcn_2 = corr_tNcn_NL_below
    endif

    !!! Calculate the normalized correlation between variables that have
    !!! an assumed normal distribution and variables that have an assumed
    !!! (single) lognormal distribution for the ith PDF component, given their
    !!! correlation and the normalized standard deviation of the variable with
    !!! the assumed lognormal distribution.

    if ( rrainm > rr_tol ) then

       ! Normalize the correlation between w and r_r in PDF component 1.
       corr_wrr_1_n = corr_NL2NN( corr_wrr_1, sigma_rr_1_n )

       ! Normalize the correlation between w and r_r in PDF component 2.
       corr_wrr_2_n = corr_NL2NN( corr_wrr_2, sigma_rr_2_n )

       ! Normalize the correlation between t and r_r in PDF component 1.
       corr_trr_1_n = corr_NL2NN( corr_trr_1, sigma_rr_1_n )

       ! Normalize the correlation between t and r_r in PDF component 2.
       corr_trr_2_n = corr_NL2NN( corr_trr_2, sigma_rr_2_n )

    else

       ! Mean rain water mixing ratio is less than the tolerance amount.  It is
       ! considered to have a value of 0.  There is not any rain at this grid
       ! level.  The correlations involving rain water mixing ratio are 0 since
       ! rain water mixing ratio does not vary at this grid level.
       corr_wrr_1_n = zero
       corr_wrr_2_n = zero
       corr_trr_1_n = zero
       corr_trr_2_n = zero

    endif

    if ( Nrm > Nr_tol ) then

       ! Normalize the correlation between w and N_r in PDF component 1.
       corr_wNr_1_n = corr_NL2NN( corr_wNr_1, sigma_Nr_1_n )

       ! Normalize the correlation between w and N_r in PDF component 2.
       corr_wNr_2_n = corr_NL2NN( corr_wNr_2, sigma_Nr_2_n )

       ! Normalize the correlation between t and N_r in PDF component 1.
       corr_tNr_1_n = corr_NL2NN( corr_tNr_1, sigma_Nr_1_n )

       ! Normalize the correlation between t and N_r in PDF component 2.
       corr_tNr_2_n = corr_NL2NN( corr_tNr_2, sigma_Nr_2_n )

    else

       ! Mean rain drop concentration is less than the tolerance amount.  It is
       ! considered to have a value of 0.  There is not any rain at this grid
       ! level.  The correlations involving rain drop concentration are 0 since
       ! rain drop concentration does not vary at this grid level.
       corr_wNr_1_n = zero
       corr_wNr_2_n = zero
       corr_tNr_1_n = zero
       corr_tNr_2_n = zero

    endif

    if ( Ncnm > Ncn_tol ) then

       ! Normalize the correlation between w and N_cn in PDF component 1.
       corr_wNcn_1_n = corr_NL2NN( corr_wNcn_1, sigma_Ncn_1_n )

       ! Normalize the correlation between w and N_cn in PDF component 2.
       corr_wNcn_2_n = corr_NL2NN( corr_wNcn_2, sigma_Ncn_2_n )

       ! Normalize the correlation between t and N_cn in PDF component 1.
       corr_tNcn_1_n = corr_NL2NN( corr_tNcn_1, sigma_Ncn_1_n )

       ! Normalize the correlation between t and N_cn in PDF component 2.
       corr_tNcn_2_n = corr_NL2NN( corr_tNcn_2, sigma_Ncn_2_n )

    else

       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_wNcn_1_n = zero
       corr_wNcn_2_n = zero
       corr_tNcn_1_n = zero
       corr_tNcn_2_n = zero

    endif


    ! Calculate the covariance of vertical velocity and KK evaporation tendency.
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then

       w_KK_evap_covar &
       = covar_x_KK_evap( mu_w_1, mu_w_2, mu_s_1, mu_s_2, mu_rr_1_n, &
                          mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_w_1, &
                          sigma_w_2, sigma_s_1, sigma_s_2, sigma_rr_1_n, &
                          sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                          corr_ws_1, corr_ws_2, corr_wrr_1_n, &
                          corr_wrr_2_n, corr_wNr_1_n, corr_wNr_2_n, &
                          corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                          corr_sNr_2_n, corr_rrNr_1_n, corr_rrNr_2_n, &
                          w_mean, KK_evap_tndcy, KK_evap_coef, w_tol, &
                          mixt_frac, precip_frac_1, precip_frac_2 )

    else  ! r_r or N_r = 0.

       w_KK_evap_covar = zero

    endif

    ! Calculate the covariance of total water mixing ratio and KK evaporation
    ! tendency.
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then

       rt_KK_evap_covar  &
       = covar_rt_KK_evap( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_rr_1, &
                           mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, &
                           mu_Nr_1_n, mu_Nr_2_n, sigma_t_1, sigma_t_2, &
                           sigma_s_1, sigma_s_2, sigma_rr_1, sigma_rr_2, &
                           sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                           sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                           corr_ts_1, corr_ts_2, corr_trr_1_n, &
                           corr_trr_2_n, corr_tNr_1_n, corr_tNr_2_n, &
                           corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                           corr_sNr_2_n, corr_rrNr_1_n, corr_rrNr_2_n, &
                           mixt_frac, precip_frac_1, precip_frac_2, &
                           KK_evap_tndcy, KK_evap_coef, t_tol, crt1, crt2 )

    else  ! r_r or N_r = 0.

       rt_KK_evap_covar = zero

    endif

    ! Calculate the covariance of liquid water potential temperature and
    ! KK evaporation tendency.
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then

       thl_KK_evap_covar &
       = covar_thl_KK_evap( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_rr_1, &
                            mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, &
                            mu_Nr_1_n, mu_Nr_2_n, sigma_t_1, sigma_t_2, &
                            sigma_s_1, sigma_s_2, sigma_rr_1, sigma_rr_2, &
                            sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                            sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                            corr_ts_1, corr_ts_2, corr_trr_1_n, &
                            corr_trr_2_n, corr_tNr_1_n, corr_tNr_2_n, &
                            corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                            corr_sNr_2_n, corr_rrNr_1_n, corr_rrNr_2_n, &
                            mixt_frac, precip_frac_1, precip_frac_2, &
                            KK_evap_tndcy, KK_evap_coef, t_tol, cthl1, cthl2 )

    else  ! r_r or N_r = 0.

       thl_KK_evap_covar = zero

    endif

    ! Calculate the covariance of vertical velocity and KK autoconversion
    ! tendency.
    if ( Ncnm > Ncn_tol ) then

       w_KK_auto_covar &
       = covar_x_KK_auto( mu_w_1, mu_w_2, mu_s_1, mu_s_2, mu_Ncn_1_n, &
                          mu_Ncn_2_n, sigma_w_1, sigma_w_2, sigma_s_1, &
                          sigma_s_2, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                          corr_ws_1, corr_ws_2, corr_wNcn_1_n, &
                          corr_wNcn_2_n, corr_sNcn_1_n, corr_sNcn_2_n, &
                          w_mean, KK_auto_tndcy, KK_auto_coef, w_tol, &
                          mixt_frac, Nc0_in_cloud, l_const_Nc_in_cloud )

    else  ! N_cn = 0.

       w_KK_auto_covar = zero

    endif

    ! Calculate the covariance of total water mixing ratio and KK autoconversion
    ! tendency.
    if ( Ncnm > Ncn_tol ) then

       rt_KK_auto_covar &
       = covar_rt_KK_auto( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_Ncn_1, &
                           mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n, sigma_t_1, &
                           sigma_t_2, sigma_s_1, sigma_s_2, sigma_Ncn_1, &
                           sigma_Ncn_2, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                           corr_ts_1, corr_ts_2, corr_tNcn_1_n, &
                           corr_tNcn_2_n, corr_sNcn_1_n, corr_sNcn_2_n, &
                           KK_auto_tndcy, KK_auto_coef, t_tol, &
                           crt1, crt2, mixt_frac, Nc0_in_cloud, &
                           l_const_Nc_in_cloud )

    else  ! N_cn = 0.

       rt_KK_auto_covar = zero

    endif

    ! Calculate the covariance of liquid water potential temperature and
    ! KK autoconversion tendency.
    if ( Ncnm > Ncn_tol ) then

       thl_KK_auto_covar &
       = covar_thl_KK_auto( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_Ncn_1, &
                            mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n, sigma_t_1, &
                            sigma_t_2, sigma_s_1, sigma_s_2, sigma_Ncn_1, &
                            sigma_Ncn_2, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                            corr_ts_1, corr_ts_2, corr_tNcn_1_n, &
                            corr_tNcn_2_n, corr_sNcn_1_n, corr_sNcn_2_n, &
                            KK_auto_tndcy, KK_auto_coef, t_tol, &
                            cthl1, cthl2, mixt_frac, Nc0_in_cloud, &
                            l_const_Nc_in_cloud )

    else  ! N_cn = 0.

       thl_KK_auto_covar = zero

    endif

    ! Calculate the covariance of vertical velocity and KK accretion tendency.
    if ( rrainm > rr_tol ) then

       w_KK_accr_covar &
       = covar_x_KK_accr( mu_w_1, mu_w_2, mu_s_1, mu_s_2, mu_rr_1_n, &
                          mu_rr_2_n, sigma_w_1, sigma_w_2, sigma_s_1, &
                          sigma_s_2, sigma_rr_1_n, sigma_rr_2_n, &
                          corr_ws_1, corr_ws_2, corr_wrr_1_n, &
                          corr_wrr_2_n, corr_srr_1_n, corr_srr_2_n, &
                          w_mean, KK_accr_tndcy, KK_accr_coef, w_tol, &
                          mixt_frac, precip_frac_1, precip_frac_2 )

    else  ! r_r = 0.

       w_KK_accr_covar = zero

    endif

    ! Calculate the covariance of total water mixing ratio and KK accretion
    ! tendency.
    if ( rrainm > rr_tol ) then

       rt_KK_accr_covar &
       = covar_rt_KK_accr( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_rr_1, &
                           mu_rr_2, mu_rr_1_n, mu_rr_2_n, sigma_t_1, &
                           sigma_t_2, sigma_s_1, sigma_s_2, sigma_rr_1, &
                           sigma_rr_2, sigma_rr_1_n, sigma_rr_2_n, &
                           corr_ts_1, corr_ts_2, corr_trr_1_n, &
                           corr_trr_2_n, corr_srr_1_n, corr_srr_2_n, &
                           KK_accr_tndcy, KK_accr_coef, t_tol, crt1, &
                           crt2, mixt_frac, precip_frac_1, precip_frac_2 )

    else  ! r_r = 0.

       rt_KK_accr_covar = zero

    endif

    ! Calculate the covariance of liquid water potential temperature and
    ! KK accretion tendency.
    if ( rrainm > rr_tol ) then

       thl_KK_accr_covar &
       = covar_thl_KK_accr( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_rr_1, &
                            mu_rr_2, mu_rr_1_n, mu_rr_2_n, sigma_t_1, &
                            sigma_t_2, sigma_s_1, sigma_s_2, sigma_rr_1, &
                            sigma_rr_2, sigma_rr_1_n, sigma_rr_2_n, &
                            corr_ts_1, corr_ts_2, corr_trr_1_n, &
                            corr_trr_2_n, corr_srr_1_n, corr_srr_2_n, &
                            KK_accr_tndcy, KK_accr_coef, t_tol, cthl1, &
                            cthl2, mixt_frac, precip_frac_1, precip_frac_2 )

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
                                   w_KK_evap_covar, zt )
       endif

       ! Covariance of r_t and KK evaporation tendency.
       if ( irt_KK_evap_covar_zt > 0 ) then
          call stat_update_var_pt( irt_KK_evap_covar_zt, level, &
                                   rt_KK_evap_covar, zt )
       endif

       ! Covariance of theta_l and KK evaporation tendency.
       if ( ithl_KK_evap_covar_zt > 0 ) then
          call stat_update_var_pt( ithl_KK_evap_covar_zt, level, &
                                   thl_KK_evap_covar, zt )
       endif

       ! Covariance of w and KK autoconversion tendency.
       if ( iw_KK_auto_covar_zt > 0 ) then
          call stat_update_var_pt( iw_KK_auto_covar_zt, level, &
                                   w_KK_auto_covar, zt )
       endif

       ! Covariance of r_t and KK autoconversion tendency.
       if ( irt_KK_auto_covar_zt > 0 ) then
          call stat_update_var_pt( irt_KK_auto_covar_zt, level, &
                                   rt_KK_auto_covar, zt )
       endif

       ! Covariance of theta_l and KK autoconversion tendency.
       if ( ithl_KK_auto_covar_zt > 0 ) then
          call stat_update_var_pt( ithl_KK_auto_covar_zt, level, &
                                   thl_KK_auto_covar, zt )
       endif

       ! Covariance of w and KK accretion tendency.
       if ( iw_KK_auto_covar_zt > 0 ) then
          call stat_update_var_pt( iw_KK_accr_covar_zt, level, &
                                   w_KK_accr_covar, zt )
       endif

       ! Covariance of r_t and KK accretion tendency.
       if ( irt_KK_auto_covar_zt > 0 ) then
          call stat_update_var_pt( irt_KK_accr_covar_zt, level, &
                                   rt_KK_accr_covar, zt )
       endif

       ! Covariance of theta_l and KK accretion tendency.
       if ( ithl_KK_auto_covar_zt > 0 ) then
          call stat_update_var_pt( ithl_KK_accr_covar_zt, level, &
                                   thl_KK_accr_covar, zt )
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
  function covar_x_KK_evap( mu_x_1, mu_x_2, mu_s_1, mu_s_2, mu_rr_1_n, &
                            mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_x_1, &
                            sigma_x_2, sigma_s_1, sigma_s_2, sigma_rr_1_n, &
                            sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                            corr_xs_1, corr_xs_2, corr_xrr_1_n, &
                            corr_xrr_2_n, corr_xNr_1_n, corr_xNr_2_n, &
                            corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                            corr_sNr_2_n, corr_rrNr_1_n, corr_rrNr_2_n, &
                            x_mean, KK_evap_tndcy, KK_evap_coef, x_tol, &
                            mixt_frac, precip_frac_1, precip_frac_2 )

    ! Description:
    ! This function calculates the covariance between x and KK evaporation
    ! tendency, which can be written as < x'((dr_r/dt)_KKevap)' >, or more
    ! simply as < x'KK_evap' >.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one  ! Constant(s)

    use parameters_microphys, only: &
        KK_evap_Supersat_exp, & ! Variable(s)
        KK_evap_rr_exp,       &
        KK_evap_Nr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x_1,        & ! Mean of x (1st PDF component)                       [-]
      mu_x_2,        & ! Mean of x (2nd PDF component)                       [-]
      mu_s_1,        & ! Mean of s (1st PDF component)                       [-]
      mu_s_2,        & ! Mean of s (2nd PDF component)                       [-]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) in-precip (ip)    [-]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip                [-]
      mu_Nr_1_n,     & ! Mean of ln Nr (1st PDF component) ip                [-]
      mu_Nr_2_n,     & ! Mean of ln Nr (2nd PDF component) ip                [-]
      sigma_x_1,     & ! Standard deviation of x (1st PDF component)         [-]
      sigma_x_2,     & ! Standard deviation of x (2nd PDF component)         [-]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)         [-]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)         [-]
      sigma_rr_1_n,  & ! Standard deviation of ln rr (1st PDF component) ip  [-]
      sigma_rr_2_n,  & ! Standard deviation of ln rr (2nd PDF component) ip  [-]
      sigma_Nr_1_n,  & ! Standard deviation of ln Nr (1st PDF component) ip  [-]
      sigma_Nr_2_n,  & ! Standard deviation of ln Nr (2nd PDF component) ip  [-]
      corr_xs_1,     & ! Correlation between x and s (1st PDF component)     [-]
      corr_xs_2,     & ! Correlation between x and s (2nd PDF component)     [-]
      corr_xrr_1_n,  & ! Correlation between x and ln rr (1st PDF comp.) ip  [-]
      corr_xrr_2_n,  & ! Correlation between x and ln rr (2nd PDF comp.) ip  [-]
      corr_xNr_1_n,  & ! Correlation between x and ln Nr (1st PDF comp.) ip  [-]
      corr_xNr_2_n,  & ! Correlation between x and ln Nr (2nd PDF comp.) ip  [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF comp.) ip  [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF comp.) ip  [-]
      corr_rrNr_1_n, & ! Correlation btwn. ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rrNr_2_n, & ! Correlation btwn. ln rr & ln Nr (2nd PDF comp.) ip  [-]
      x_mean,        & ! Mean of x (overall)                        [units vary]
      KK_evap_tndcy, & ! KK evaporation tendency                     [(kg/kg)/s]
      KK_evap_coef,  & ! KK evaporation coefficient                  [(kg/kg)/s]
      x_tol,         & ! Tolerance value of x                       [units vary]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)          [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_x_KK_evap  ! Covariance between x and KK evaporation tendency    [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on s                                           [-]
      beta_exp,  & ! Exponent on r_r                                         [-]
      gamma_exp    ! Exponent on N_r                                         [-]


    ! Values of the KK exponents.
    alpha_exp = KK_evap_Supersat_exp
    beta_exp  = KK_evap_rr_exp
    gamma_exp = KK_evap_Nr_exp

    ! Calculate the covariance of x and KK evaporation tendency.
    covar_x_KK_evap  &
    = KK_evap_coef  &
      * ( mixt_frac &
          * precip_frac_1 &
          * quadrivar_NNLL_covar_eq( mu_x_1, mu_s_1, mu_rr_1_n, mu_Nr_1_n, &
                                     sigma_x_1, sigma_s_1, sigma_rr_1_n, &
                                     sigma_Nr_1_n, corr_xs_1, corr_xrr_1_n, &
                                     corr_xNr_1_n, corr_srr_1_n, &
                                     corr_sNr_1_n, corr_rrNr_1_n, x_mean, &
                                     KK_evap_tndcy, KK_evap_coef, x_tol, &
                                     alpha_exp, beta_exp, gamma_exp ) &
        + ( one - mixt_frac ) &
          * precip_frac_2 &
          * quadrivar_NNLL_covar_eq( mu_x_2, mu_s_2, mu_rr_2_n, mu_Nr_2_n, &
                                     sigma_x_2, sigma_s_2, sigma_rr_2_n, &
                                     sigma_Nr_2_n, corr_xs_2, corr_xrr_2_n, &
                                     corr_xNr_2_n, corr_srr_2_n, &
                                     corr_sNr_2_n, corr_rrNr_2_n, x_mean, &
                                     KK_evap_tndcy, KK_evap_coef, x_tol, &
                                     alpha_exp, beta_exp, gamma_exp ) &
        )


    return

  end function covar_x_KK_evap

  !=============================================================================
  function covar_rt_KK_evap( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_rr_1, &
                             mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, &
                             mu_Nr_1_n, mu_Nr_2_n, sigma_t_1, sigma_t_2, &
                             sigma_s_1, sigma_s_2, sigma_rr_1, sigma_rr_2, &
                             sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                             sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                             corr_ts_1, corr_ts_2, corr_trr_1_n, &
                             corr_trr_2_n, corr_tNr_1_n, corr_tNr_2_n, &
                             corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                             corr_sNr_2_n, corr_rrNr_1_n, corr_rrNr_2_n, &
                             mixt_frac, precip_frac_1, precip_frac_2, &
                             KK_evap_tndcy, KK_evap_coef, t_tol, crt1, crt2 )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two, & ! Constant(s)
        one

    use KK_upscaled_means, only:  &
        trivar_NLL_mean_eq  ! Procedure

    use parameters_microphys, only: &
        KK_evap_Supersat_exp, & ! Variable(s)
        KK_evap_rr_exp,       &
        KK_evap_Nr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_t_1,        & ! Mean of t (1st PDF component)                   [kg/kg]
      mu_t_2,        & ! Mean of t (2nd PDF component)                   [kg/kg]
      mu_s_1,        & ! Mean of s (1st PDF component)                   [kg/kg]
      mu_s_2,        & ! Mean of s (2nd PDF component)                   [kg/kg]
      mu_rr_1,       & ! Mean of rr (1st PDF component) in-precip (ip)   [kg/kg]
      mu_rr_2,       & ! Mean of rr (2nd PDF component) ip               [kg/kg]
      mu_Nr_1,       & ! Mean of Nr (1st PDF component) ip              [num/kg]
      mu_Nr_2,       & ! Mean of Nr (2nd PDF component) ip              [num/kg]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) ip                [-]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip                [-]
      mu_Nr_1_n,     & ! Mean of ln Nr (1st PDF component) ip                [-]
      mu_Nr_2_n,     & ! Mean of ln Nr (2nd PDF component) ip                [-]
      sigma_t_1,     & ! Standard deviation of t (1st PDF component)     [kg/kg]
      sigma_t_2,     & ! Standard deviation of t (2nd PDF component)     [kg/kg]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)     [kg/kg]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)     [kg/kg]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_Nr_1,    & ! Standard deviation of Nr (1st PDF component) ip  [#/kg]
      sigma_Nr_2,    & ! Standard deviation of Nr (2nd PDF component) ip  [#/kg]
      sigma_rr_1_n,  & ! Standard deviation of ln rr (1st PDF component) ip  [-]
      sigma_rr_2_n,  & ! Standard deviation of ln rr (2nd PDF component) ip  [-]
      sigma_Nr_1_n,  & ! Standard deviation of ln Nr (1st PDF component) ip  [-]
      sigma_Nr_2_n,  & ! Standard deviation of ln Nr (2nd PDF component) ip  [-]
      corr_ts_1,     & ! Correlation between t and s (1st PDF component)     [-]
      corr_ts_2,     & ! Correlation between t and s (2nd PDF component)     [-]
      corr_trr_1_n,  & ! Correlation between t and ln rr (1st PDF comp.) ip  [-]
      corr_trr_2_n,  & ! Correlation between t and ln rr (2nd PDF comp.) ip  [-]
      corr_tNr_1_n,  & ! Correlation between t and ln Nr (1st PDF comp.) ip  [-]
      corr_tNr_2_n,  & ! Correlation between t and ln Nr (2nd PDF comp.) ip  [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF comp.) ip  [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF comp.) ip  [-]
      corr_rrNr_1_n, & ! Correlation btwn. ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rrNr_2_n, & ! Correlation btwn. ln rr & ln Nr (2nd PDF comp.) ip  [-]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)          [-]

    real( kind = core_rknd ), intent(in) :: &
      KK_evap_tndcy, & ! KK evaporation tendency                     [(kg/kg)/s]
      KK_evap_coef,  & ! KK evaporation coefficient                  [(kg/kg)/s]
      t_tol,         & ! Tolerance value of t                                [-]
      crt1,          & ! Coefficient c_rt (1st PDF component)                [-]
      crt2             ! Coefficient c_rt (2nd PDF component)                [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_rt_KK_evap  ! Covariance between r_t and KK evaporation tendency [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on s                                           [-]
      beta_exp,  & ! Exponent on r_r                                         [-]
      gamma_exp    ! Exponent on N_r                                         [-]


    ! Values of the KK exponents.
    alpha_exp = KK_evap_Supersat_exp
    beta_exp  = KK_evap_rr_exp
    gamma_exp = KK_evap_Nr_exp

    ! Calculate the covariance of r_t and KK evaporation tendency.
    covar_rt_KK_evap  &
    = KK_evap_coef  &
      * ( mixt_frac * precip_frac_1 * ( one / ( two * crt1 ) )  &
          * ( quadrivar_NNLL_covar_eq( mu_t_1, mu_s_1, mu_rr_1_n, mu_Nr_1_n, &
                                       sigma_t_1, sigma_s_1, sigma_rr_1_n, &
                                       sigma_Nr_1_n, corr_ts_1, corr_trr_1_n, &
                                       corr_tNr_1_n, corr_srr_1_n, &
                                       corr_sNr_1_n, corr_rrNr_1_n, mu_t_1, &
                                       KK_evap_tndcy, KK_evap_coef, t_tol, &
                                       alpha_exp, beta_exp, gamma_exp )  &
            + trivar_NLL_mean_eq( mu_s_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, &
                                  mu_Nr_1_n, sigma_s_1, sigma_rr_1, &
                                  sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, &
                                  corr_srr_1_n, corr_sNr_1_n, corr_rrNr_1_n, &
                                  alpha_exp + one, beta_exp, gamma_exp )  &
            - mu_s_1  &
              * trivar_NLL_mean_eq( mu_s_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, &
                                    mu_Nr_1_n, sigma_s_1, sigma_rr_1, &
                                    sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, &
                                    corr_srr_1_n, corr_sNr_1_n, corr_rrNr_1_n, &
                                    alpha_exp, beta_exp, gamma_exp )  &
            )  &
        + ( one - mixt_frac ) * precip_frac_2 * ( one / ( two * crt2 ) )  &
          * ( quadrivar_NNLL_covar_eq( mu_t_2, mu_s_2, mu_rr_2_n, mu_Nr_2_n, &
                                       sigma_t_2, sigma_s_2, sigma_rr_2_n, &
                                       sigma_Nr_2_n, corr_ts_2, corr_trr_2_n, &
                                       corr_tNr_2_n, corr_srr_2_n, &
                                       corr_sNr_2_n, corr_rrNr_2_n, mu_t_2, &
                                       KK_evap_tndcy, KK_evap_coef, t_tol, &
                                       alpha_exp, beta_exp, gamma_exp )  &
            + trivar_NLL_mean_eq( mu_s_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, &
                                  mu_Nr_2_n, sigma_s_2, sigma_rr_2, &
                                  sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, &
                                  corr_srr_2_n, corr_sNr_2_n, corr_rrNr_2_n, &
                                  alpha_exp + one, beta_exp, gamma_exp )  &
            - mu_s_2  &
              * trivar_NLL_mean_eq( mu_s_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, &
                                    mu_Nr_2_n, sigma_s_2, sigma_rr_2, &
                                    sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, &
                                    corr_srr_2_n, corr_sNr_2_n, corr_rrNr_2_n, &
                                    alpha_exp, beta_exp, gamma_exp )  &
            )  &
        )


    return

  end function covar_rt_KK_evap

  !=============================================================================
  function covar_thl_KK_evap( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_rr_1, &
                              mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, &
                              mu_Nr_1_n, mu_Nr_2_n, sigma_t_1, sigma_t_2, &
                              sigma_s_1, sigma_s_2, sigma_rr_1, sigma_rr_2, &
                              sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                              sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                              corr_ts_1, corr_ts_2, corr_trr_1_n, &
                              corr_trr_2_n, corr_tNr_1_n, corr_tNr_2_n, &
                              corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                              corr_sNr_2_n, corr_rrNr_1_n, corr_rrNr_2_n, &
                              mixt_frac, precip_frac_1, precip_frac_2, &
                              KK_evap_tndcy, KK_evap_coef, t_tol, cthl1, cthl2 )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two, & ! Constant(s)
        one

    use KK_upscaled_means, only:  &
        trivar_NLL_mean_eq  ! Procedure

    use parameters_microphys, only: &
        KK_evap_Supersat_exp, & ! Variable(s)
        KK_evap_rr_exp,       &
        KK_evap_Nr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_t_1,        & ! Mean of t (1st PDF component)                   [kg/kg]
      mu_t_2,        & ! Mean of t (2nd PDF component)                   [kg/kg]
      mu_s_1,        & ! Mean of s (1st PDF component)                   [kg/kg]
      mu_s_2,        & ! Mean of s (2nd PDF component)                   [kg/kg]
      mu_rr_1,       & ! Mean of rr (1st PDF component) in-precip (ip)   [kg/kg]
      mu_rr_2,       & ! Mean of rr (2nd PDF component) ip               [kg/kg]
      mu_Nr_1,       & ! Mean of Nr (1st PDF component) ip              [num/kg]
      mu_Nr_2,       & ! Mean of Nr (2nd PDF component) ip              [num/kg]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) ip                [-]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip                [-]
      mu_Nr_1_n,     & ! Mean of ln Nr (1st PDF component) ip                [-]
      mu_Nr_2_n,     & ! Mean of ln Nr (2nd PDF component) ip                [-]
      sigma_t_1,     & ! Standard deviation of t (1st PDF component)     [kg/kg]
      sigma_t_2,     & ! Standard deviation of t (2nd PDF component)     [kg/kg]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)     [kg/kg]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)     [kg/kg]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_Nr_1,    & ! Standard deviation of Nr (1st PDF component) ip  [#/kg]
      sigma_Nr_2,    & ! Standard deviation of Nr (2nd PDF component) ip  [#/kg]
      sigma_rr_1_n,  & ! Standard deviation of ln rr (1st PDF component) ip  [-]
      sigma_rr_2_n,  & ! Standard deviation of ln rr (2nd PDF component) ip  [-]
      sigma_Nr_1_n,  & ! Standard deviation of ln Nr (1st PDF component) ip  [-]
      sigma_Nr_2_n,  & ! Standard deviation of ln Nr (2nd PDF component) ip  [-]
      corr_ts_1,     & ! Correlation between t and s (1st PDF component)     [-]
      corr_ts_2,     & ! Correlation between t and s (2nd PDF component)     [-]
      corr_trr_1_n,  & ! Correlation between t and ln rr (1st PDF comp.) ip  [-]
      corr_trr_2_n,  & ! Correlation between t and ln rr (2nd PDF comp.) ip  [-]
      corr_tNr_1_n,  & ! Correlation between t and ln Nr (1st PDF comp.) ip  [-]
      corr_tNr_2_n,  & ! Correlation between t and ln Nr (2nd PDF comp.) ip  [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF comp.) ip  [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF comp.) ip  [-]
      corr_rrNr_1_n, & ! Correlation btwn. ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rrNr_2_n, & ! Correlation btwn. ln rr & ln Nr (2nd PDF comp.) ip  [-]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)          [-]

    real( kind = core_rknd ), intent(in) :: &
      KK_evap_tndcy, & ! KK evaporation tendency                     [(kg/kg)/s]
      KK_evap_coef,  & ! KK evaporation coefficient                  [(kg/kg)/s]
      t_tol,         & ! Tolerance value of t                                [-]
      cthl1,         & ! Coefficient c_thl (1st PDF component)               [-]
      cthl2            ! Coefficient c_thl (2nd PDF component)               [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_thl_KK_evap  ! Covariance between th_l and KK evap. tendency     [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on s                                           [-]
      beta_exp,  & ! Exponent on r_r                                         [-]
      gamma_exp    ! Exponent on N_r                                         [-]


    ! Values of the KK exponents.
    alpha_exp = KK_evap_Supersat_exp
    beta_exp  = KK_evap_rr_exp
    gamma_exp = KK_evap_Nr_exp

    ! Calculate the covariance of th_l and KK evaporation tendency.
    covar_thl_KK_evap  &
    = KK_evap_coef  &
      * ( mixt_frac * precip_frac_1 * ( one / ( two * cthl1 ) )  &
          * ( quadrivar_NNLL_covar_eq( mu_t_1, mu_s_1, mu_rr_1_n, mu_Nr_1_n, &
                                       sigma_t_1, sigma_s_1, sigma_rr_1_n, &
                                       sigma_Nr_1_n, corr_ts_1, corr_trr_1_n, &
                                       corr_tNr_1_n, corr_srr_1_n, &
                                       corr_sNr_1_n, corr_rrNr_1_n, mu_t_1, &
                                       KK_evap_tndcy, KK_evap_coef, t_tol, &
                                       alpha_exp, beta_exp, gamma_exp )  &
            - trivar_NLL_mean_eq( mu_s_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, &
                                  mu_Nr_1_n, sigma_s_1, sigma_rr_1, &
                                  sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, &
                                  corr_srr_1_n, corr_sNr_1_n, corr_rrNr_1_n, &
                                  alpha_exp + one, beta_exp, gamma_exp )  &
            + mu_s_1  &
              * trivar_NLL_mean_eq( mu_s_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, &
                                    mu_Nr_1_n, sigma_s_1, sigma_rr_1, &
                                    sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, &
                                    corr_srr_1_n, corr_sNr_1_n, corr_rrNr_1_n, &
                                    alpha_exp, beta_exp, gamma_exp )  &
            )  &
        + ( one - mixt_frac ) * precip_frac_2 * ( one / ( two * cthl2 ) )  &
          * ( quadrivar_NNLL_covar_eq( mu_t_2, mu_s_2, mu_rr_2_n, mu_Nr_2_n, &
                                       sigma_t_2, sigma_s_2, sigma_rr_2_n, &
                                       sigma_Nr_2_n, corr_ts_2, corr_trr_2_n, &
                                       corr_tNr_2_n, corr_srr_2_n, &
                                       corr_sNr_2_n, corr_rrNr_2_n, mu_t_2, &
                                       KK_evap_tndcy, KK_evap_coef, t_tol, &
                                       alpha_exp, beta_exp, gamma_exp )  &
            - trivar_NLL_mean_eq( mu_s_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, &
                                  mu_Nr_2_n, sigma_s_2, sigma_rr_2, &
                                  sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, &
                                  corr_srr_2_n, corr_sNr_2_n, corr_rrNr_2_n, &
                                  alpha_exp + one, beta_exp, gamma_exp )  &
            + mu_s_2  &
              * trivar_NLL_mean_eq( mu_s_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, &
                                    mu_Nr_2_n, sigma_s_2, sigma_rr_2, &
                                    sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, &
                                    corr_srr_2_n, corr_sNr_2_n, corr_rrNr_2_n, &
                                    alpha_exp, beta_exp, gamma_exp )  &
            )  &
        )


    return

  end function covar_thl_KK_evap

  !=============================================================================
  function covar_x_KK_auto( mu_x_1, mu_x_2, mu_s_1, mu_s_2, mu_Ncn_1_n, &
                            mu_Ncn_2_n, sigma_x_1, sigma_x_2, sigma_s_1, &
                            sigma_s_2, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                            corr_xs_1, corr_xs_2, corr_xNcn_1_n, &
                            corr_xNcn_2_n, corr_sNcn_1_n, corr_sNcn_2_n, &
                            x_mean, KK_auto_tndcy, KK_auto_coef, x_tol, &
                            mixt_frac, Nc0_in_cloud, l_const_Nc_in_cloud )

    ! Description:
    ! This function calculates the correlation between x and KK autoconversion
    ! tendency, which can be written as < x'((dr_r/dt)_KKauto)' >, or more
    ! simply as < x'KK_auto' >.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one  ! Constant(s)

    use parameters_microphys, only: &
        KK_auto_rc_exp, & ! Variable(s)
        KK_auto_Nc_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x_1,        & ! Mean of x (1st PDF component)                [un. vary]
      mu_x_2,        & ! Mean of x (2nd PDF component)                [un. vary]
      mu_s_1,        & ! Mean of s (1st PDF component)                   [kg/kg]
      mu_s_2,        & ! Mean of s (2nd PDF component)                   [kg/kg]
      mu_Ncn_1_n,    & ! Mean of ln Ncn (1st PDF component)         [ln(num/kg)]
      mu_Ncn_2_n,    & ! Mean of ln Ncn (2nd PDF component)         [ln(num/kg)]
      sigma_x_1,     & ! Standard deviation of x (1st PDF component)  [un. vary]
      sigma_x_2,     & ! Standard deviation of x (2nd PDF component)  [un. vary]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)     [kg/kg]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)     [kg/kg]
      sigma_Ncn_1_n, & ! Standard deviation of ln Ncn (1st PDF comp.) [ln(#/kg)]
      sigma_Ncn_2_n, & ! Standard deviation of ln Ncn (2nd PDF comp.) [ln(#/kg)]
      corr_xs_1,     & ! Correlation between x and s (1st PDF component)     [-]
      corr_xs_2,     & ! Correlation between x and s (2nd PDF component)     [-]
      corr_xNcn_1_n, & ! Correlation between x and ln Ncn (1st PDF comp.)    [-]
      corr_xNcn_2_n, & ! Correlation between x and ln Ncn (2nd PDF comp.)    [-]
      corr_sNcn_1_n, & ! Correlation between s and ln Ncn (1st PDF comp.)    [-]
      corr_sNcn_2_n, & ! Correlation between s and ln Ncn (2nd PDF comp.)    [-]
      x_mean,        & ! Mean of x (overall)                          [un. vary]
      KK_auto_tndcy, & ! KK autoconversion tendency                  [(kg/kg)/s]
      KK_auto_coef,  & ! KK autoconversion coefficient               [(kg/kg)/s]
      x_tol,         & ! Tolerance value of x                         [un. vary]
      mixt_frac,     & ! Mixture fraction                                    [-]
      Nc0_in_cloud     ! Constant in-cloud value of cloud droplet conc. [num/kg]

    logical, intent(in) :: &
      l_const_Nc_in_cloud  ! Flag to use a constant value of N_c within cloud

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_x_KK_auto  ! Covariance between x and KK autoconversion tendency [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on s                                           [-]
      beta_exp     ! Exponent on N_c                                         [-]


    ! Values of the KK exponents.
    alpha_exp = KK_auto_rc_exp
    beta_exp  = KK_auto_Nc_exp

    ! Calculate the covariance of x and KK autoconversion tendency.
    if ( l_const_Nc_in_cloud ) then

       covar_x_KK_auto  &
       = KK_auto_coef &
         * ( mixt_frac &
             * trivar_NNL_covar_eq_Nc0( mu_x_1, mu_s_1, Nc0_in_cloud, &
                                        sigma_x_1, sigma_s_1, corr_xs_1, &
                                        x_mean, KK_auto_tndcy, KK_auto_coef, &
                                        x_tol, alpha_exp, beta_exp ) &
           + ( one - mixt_frac ) &
             * trivar_NNL_covar_eq_Nc0( mu_x_2, mu_s_2, Nc0_in_cloud, &
                                        sigma_x_2, sigma_s_2, corr_xs_2, &
                                        x_mean, KK_auto_tndcy, KK_auto_coef, &
                                        x_tol, alpha_exp, beta_exp ) &
           )

    else

       covar_x_KK_auto  &
       = KK_auto_coef &
         * ( mixt_frac &
             * trivar_NNL_covar_eq( mu_x_1, mu_s_1, mu_Ncn_1_n, &
                                    sigma_x_1, sigma_s_1, sigma_Ncn_1_n, &
                                    corr_xs_1, corr_xNcn_1_n, corr_sNcn_1_n, &
                                    x_mean, KK_auto_tndcy, KK_auto_coef, &
                                    x_tol, alpha_exp, beta_exp ) &
           + ( one - mixt_frac ) &
             * trivar_NNL_covar_eq( mu_x_2, mu_s_2, mu_Ncn_2_n, &
                                    sigma_x_2, sigma_s_2, sigma_Ncn_2_n, &
                                    corr_xs_2, corr_xNcn_2_n, corr_sNcn_2_n, &
                                    x_mean, KK_auto_tndcy, KK_auto_coef, &
                                    x_tol, alpha_exp, beta_exp ) &
           )

    endif


    return

  end function covar_x_KK_auto

  !=============================================================================
  function covar_rt_KK_auto( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_Ncn_1, &
                             mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n, sigma_t_1, &
                             sigma_t_2, sigma_s_1, sigma_s_2, sigma_Ncn_1, &
                             sigma_Ncn_2, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                             corr_ts_1, corr_ts_2, corr_tNcn_1_n, &
                             corr_tNcn_2_n, corr_sNcn_1_n, corr_sNcn_2_n, &
                             KK_auto_tndcy, KK_auto_coef, t_tol, &
                             crt1, crt2, mixt_frac, Nc0_in_cloud, &
                             l_const_Nc_in_cloud )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two,    & ! Constant(s)
        one,    & 
        Nc_tol

    use KK_upscaled_means, only:  &
        bivar_NL_mean_eq,     & ! Procedure(s)
        bivar_NL_mean_eq_Nc0

    use parameters_microphys, only: &
        KK_auto_rc_exp, & ! Variable(s)
        KK_auto_Nc_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_t_1,        & ! Mean of t (1st PDF component)                   [kg/kg]
      mu_t_2,        & ! Mean of t (2nd PDF component)                   [kg/kg]
      mu_s_1,        & ! Mean of s (1st PDF component)                   [kg/kg]
      mu_s_2,        & ! Mean of s (2nd PDF component)                   [kg/kg]
      mu_Ncn_1,      & ! Mean of Ncn (1st PDF component)                [num/kg]
      mu_Ncn_2,      & ! Mean of Ncn (2nd PDF component)                [num/kg]
      mu_Ncn_1_n,    & ! Mean of ln Ncn (1st PDF component)         [ln(num/kg)]
      mu_Ncn_2_n,    & ! Mean of ln Ncn (2nd PDF component)         [ln(num/kg)]
      sigma_t_1,     & ! Standard deviation of t (1st PDF component)     [kg/kg]
      sigma_t_2,     & ! Standard deviation of t (2nd PDF component)     [kg/kg]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)     [kg/kg]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)     [kg/kg]
      sigma_Ncn_1,   & ! Standard deviation of Ncn (1st PDF component)  [num/kg]
      sigma_Ncn_2,   & ! Standard deviation of Ncn (2nd PDF component)  [num/kg]
      sigma_Ncn_1_n, & ! Standard deviation of ln Ncn (1st PDF comp.) [ln(#/kg)]
      sigma_Ncn_2_n, & ! Standard deviation of ln Ncn (2nd PDF comp.) [ln(#/kg)]
      corr_ts_1,     & ! Correlation between t and s (1st PDF component)     [-]
      corr_ts_2,     & ! Correlation between t and s (2nd PDF component)     [-]
      corr_tNcn_1_n, & ! Correlation between t and ln Ncn (1st PDF comp.)    [-]
      corr_tNcn_2_n, & ! Correlation between t and ln Ncn (2nd PDF comp.)    [-]
      corr_sNcn_1_n, & ! Correlation between s and ln Ncn (1st PDF comp.)    [-]
      corr_sNcn_2_n, & ! Correlation between s and ln Ncn (2nd PDF comp.)    [-]
      KK_auto_tndcy, & ! KK autoconversion tendency                  [(kg/kg)/s]
      KK_auto_coef,  & ! KK autoconversion coefficient               [(kg/kg)/s]
      t_tol,         & ! Tolerance value of t                            [kg/kg]
      crt1,          & ! Coefficient c_rt (1st PDF component)                [-]
      crt2,          & ! Coefficient c_rt (2nd PDF component)                [-]
      mixt_frac,     & ! Mixture fraction                                    [-]
      Nc0_in_cloud     ! Constant in-cloud value of cloud droplet conc. [num/kg]

    logical, intent(in) :: &
      l_const_Nc_in_cloud  ! Flag to use a constant value of N_c within cloud

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_rt_KK_auto  ! Covariance between r_t and KK autoconv. tendency   [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on s                                           [-]
      beta_exp     ! Exponent on N_c                                         [-]


    ! Values of the KK exponents.
    alpha_exp = KK_auto_rc_exp
    beta_exp  = KK_auto_Nc_exp

    ! Calculate the covariance of r_t and KK autoconversion tendency.
    if ( l_const_Nc_in_cloud ) then

       covar_rt_KK_auto  &
       = KK_auto_coef  &
         * ( mixt_frac * ( one / ( two * crt1 ) )  &
             * ( trivar_NNL_covar_eq_Nc0( mu_t_1, mu_s_1, Nc0_in_cloud, &
                                          sigma_t_1, sigma_s_1, corr_ts_1, &
                                          mu_t_1, KK_auto_tndcy, KK_auto_coef, &
                                          t_tol, alpha_exp, beta_exp )  &
               + bivar_NL_mean_eq_Nc0( mu_s_1, Nc0_in_cloud, sigma_s_1, &
                                       alpha_exp + one, beta_exp )  &
               - mu_s_1  &
                 * bivar_NL_mean_eq_Nc0( mu_s_1, Nc0_in_cloud, sigma_s_1, &
                                         alpha_exp, beta_exp )  &
               )  &
           + ( one - mixt_frac ) * ( one / ( two * crt2 ) )  &
             * ( trivar_NNL_covar_eq_Nc0( mu_t_2, mu_s_2, Nc0_in_cloud, &
                                          sigma_t_2, sigma_s_2, corr_ts_2, &
                                          mu_t_2, KK_auto_tndcy, KK_auto_coef, &
                                          t_tol, alpha_exp, beta_exp )  &
               + bivar_NL_mean_eq_Nc0( mu_s_2, Nc0_in_cloud, sigma_s_2, &
                                       alpha_exp + one, beta_exp )  &
               - mu_s_2  &
                 * bivar_NL_mean_eq_Nc0( mu_s_2, Nc0_in_cloud, sigma_s_2, &
                                         alpha_exp, beta_exp )  &
               )  &
           )

    else

       covar_rt_KK_auto  &
       = KK_auto_coef  &
         * ( mixt_frac * ( one / ( two * crt1 ) )  &
             * ( trivar_NNL_covar_eq( mu_t_1, mu_s_1, mu_Ncn_1_n, &
                                      sigma_t_1, sigma_s_1, sigma_Ncn_1_n, &
                                      corr_ts_1, corr_tNcn_1_n, corr_sNcn_1_n, &
                                      mu_t_1, KK_auto_tndcy, KK_auto_coef, &
                                      t_tol, alpha_exp, beta_exp )  &
               + bivar_NL_mean_eq( mu_s_1, mu_Ncn_1, mu_Ncn_1_n, sigma_s_1, &
                                   sigma_Ncn_1, sigma_Ncn_1_n, corr_sNcn_1_n, &
                                   Nc_tol, alpha_exp + one, beta_exp )  &
               - mu_s_1  &
                 * bivar_NL_mean_eq( mu_s_1, mu_Ncn_1, mu_Ncn_1_n, sigma_s_1, &
                                     sigma_Ncn_1, sigma_Ncn_1_n, corr_sNcn_1_n,&
                                     Nc_tol, alpha_exp, beta_exp )  &
            )  &
           + ( one - mixt_frac ) * ( one / ( two * crt2 ) )  &
             * ( trivar_NNL_covar_eq( mu_t_2, mu_s_2, mu_Ncn_2_n, &
                                      sigma_t_2, sigma_s_2, sigma_Ncn_2_n, &
                                      corr_ts_2, corr_tNcn_2_n, corr_sNcn_2_n, &
                                      mu_t_2, KK_auto_tndcy, KK_auto_coef, &
                                      t_tol, alpha_exp, beta_exp )  &
               + bivar_NL_mean_eq( mu_s_2, mu_Ncn_2, mu_Ncn_2_n, sigma_s_2, &
                                   sigma_Ncn_2, sigma_Ncn_2_n, corr_sNcn_2_n, &
                                   Nc_tol, alpha_exp + one, beta_exp )  &
               - mu_s_2  &
                 * bivar_NL_mean_eq( mu_s_2, mu_Ncn_2, mu_Ncn_2_n, sigma_s_2, &
                                     sigma_Ncn_2, sigma_Ncn_2_n, corr_sNcn_2_n,&
                                     Nc_tol, alpha_exp, beta_exp )  &
               )  &
           )

    endif


    return

  end function covar_rt_KK_auto

  !=============================================================================
  function covar_thl_KK_auto( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_Ncn_1, &
                              mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n, sigma_t_1, &
                              sigma_t_2, sigma_s_1, sigma_s_2, sigma_Ncn_1, &
                              sigma_Ncn_2, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                              corr_ts_1, corr_ts_2, corr_tNcn_1_n, &
                              corr_tNcn_2_n, corr_sNcn_1_n, corr_sNcn_2_n, &
                              KK_auto_tndcy, KK_auto_coef, t_tol, &
                              cthl1, cthl2, mixt_frac, Nc0_in_cloud, &
                              l_const_Nc_in_cloud )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two,    & ! Constant(s)
        one,    &
        Nc_tol

    use KK_upscaled_means, only:  &
        bivar_NL_mean_eq,     & ! Procedure(s)
        bivar_NL_mean_eq_Nc0

    use parameters_microphys, only: &
        KK_auto_rc_exp, & ! Variable(s)
        KK_auto_Nc_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_t_1,        & ! Mean of t (1st PDF component)                   [kg/kg]
      mu_t_2,        & ! Mean of t (2nd PDF component)                   [kg/kg]
      mu_s_1,        & ! Mean of s (1st PDF component)                   [kg/kg]
      mu_s_2,        & ! Mean of s (2nd PDF component)                   [kg/kg]
      mu_Ncn_1,      & ! Mean of Ncn (1st PDF component)                [num/kg]
      mu_Ncn_2,      & ! Mean of Ncn (2nd PDF component)                [num/kg]
      mu_Ncn_1_n,    & ! Mean of ln Ncn (1st PDF component)         [ln(num/kg)]
      mu_Ncn_2_n,    & ! Mean of ln Ncn (2nd PDF component)         [ln(num/kg)]
      sigma_t_1,     & ! Standard deviation of t (1st PDF component)     [kg/kg]
      sigma_t_2,     & ! Standard deviation of t (2nd PDF component)     [kg/kg]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)     [kg/kg]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)     [kg/kg]
      sigma_Ncn_1,   & ! Standard deviation of Ncn (1st PDF component)  [num/kg]
      sigma_Ncn_2,   & ! Standard deviation of Ncn (2nd PDF component)  [num/kg]
      sigma_Ncn_1_n, & ! Standard deviation of ln Ncn (1st PDF comp.) [ln(#/kg)]
      sigma_Ncn_2_n, & ! Standard deviation of ln Ncn (2nd PDF comp.) [ln(#/kg)]
      corr_ts_1,     & ! Correlation between t and s (1st PDF component)     [-]
      corr_ts_2,     & ! Correlation between t and s (2nd PDF component)     [-]
      corr_tNcn_1_n, & ! Correlation between t and ln Ncn (1st PDF comp.)    [-]
      corr_tNcn_2_n, & ! Correlation between t and ln Ncn (2nd PDF comp.)    [-]
      corr_sNcn_1_n, & ! Correlation between s and ln Ncn (1st PDF comp.)    [-]
      corr_sNcn_2_n, & ! Correlation between s and ln Ncn (2nd PDF comp.)    [-]
      KK_auto_tndcy, & ! KK autoconversion tendency                  [(kg/kg)/s]
      KK_auto_coef,  & ! KK autoconversion coefficient               [(kg/kg)/s]
      t_tol,         & ! Tolerance value of t                            [kg/kg]
      cthl1,         & ! Coefficient c_thl (1st PDF component)               [-]
      cthl2,         & ! Coefficient c_thl (2nd PDF component)               [-]
      mixt_frac,     & ! Mixture fraction                                    [-]
      Nc0_in_cloud     ! Constant in-cloud value of cloud droplet conc. [num/kg]

    logical, intent(in) :: &
      l_const_Nc_in_cloud  ! Flag to use a constant value of N_c within cloud

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_thl_KK_auto  ! Covariance between th_l and KK autoconv. tendency [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on s                                           [-]
      beta_exp     ! Exponent on N_c                                         [-]


    ! Values of the KK exponents.
    alpha_exp = KK_auto_rc_exp
    beta_exp  = KK_auto_Nc_exp

    ! Calculate the covariance of th_l and KK autoconversion tendency.
    if ( l_const_Nc_in_cloud ) then

       covar_thl_KK_auto  &
       = KK_auto_coef  &
         * ( mixt_frac * ( one / ( two * cthl1 ) )  &
             * ( trivar_NNL_covar_eq_Nc0( mu_t_1, mu_s_1, Nc0_in_cloud, &
                                          sigma_t_1, sigma_s_1, corr_ts_1, &
                                          mu_t_1, KK_auto_tndcy, KK_auto_coef, &
                                          t_tol, alpha_exp, beta_exp )  &
               - bivar_NL_mean_eq_Nc0( mu_s_1, Nc0_in_cloud, sigma_s_1, &
                                       alpha_exp + one, beta_exp )  &
               + mu_s_1  &
                 * bivar_NL_mean_eq_Nc0( mu_s_1, Nc0_in_cloud, sigma_s_1, &
                                         alpha_exp, beta_exp )  &
               )  &
           + ( one - mixt_frac ) * ( one / ( two * cthl2 ) )  &
             * ( trivar_NNL_covar_eq_Nc0( mu_t_2, mu_s_2, Nc0_in_cloud, &
                                          sigma_t_2, sigma_s_2, corr_ts_2, &
                                          mu_t_2, KK_auto_tndcy, KK_auto_coef, &
                                          t_tol, alpha_exp, beta_exp )  &
               - bivar_NL_mean_eq_Nc0( mu_s_2, Nc0_in_cloud, sigma_s_2, &
                                       alpha_exp + one, beta_exp )  &
               + mu_s_2  &
                 * bivar_NL_mean_eq_Nc0( mu_s_2, Nc0_in_cloud, sigma_s_2, &
                                         alpha_exp, beta_exp )  &    
               )  &
           )

    else

       covar_thl_KK_auto  &
       = KK_auto_coef  &
         * ( mixt_frac * ( one / ( two * cthl1 ) )  &
             * ( trivar_NNL_covar_eq( mu_t_1, mu_s_1, mu_Ncn_1_n, &
                                      sigma_t_1, sigma_s_1, sigma_Ncn_1_n, &
                                      corr_ts_1, corr_tNcn_1_n, corr_sNcn_1_n, &
                                      mu_t_1, KK_auto_tndcy, KK_auto_coef, &
                                      t_tol, alpha_exp, beta_exp )  &
               - bivar_NL_mean_eq( mu_s_1, mu_Ncn_1, mu_Ncn_1_n, sigma_s_1, &
                                   sigma_Ncn_1, sigma_Ncn_1_n, corr_sNcn_1_n, &
                                   Nc_tol, alpha_exp + one, beta_exp )  &
               + mu_s_1  &
                 * bivar_NL_mean_eq( mu_s_1, mu_Ncn_1, mu_Ncn_1_n, sigma_s_1, &
                                     sigma_Ncn_1, sigma_Ncn_1_n, corr_sNcn_1_n,&
                                     Nc_tol, alpha_exp, beta_exp )  &
               )  &
           + ( one - mixt_frac ) * ( one / ( two * cthl2 ) )  &
             * ( trivar_NNL_covar_eq( mu_t_2, mu_s_2, mu_Ncn_2_n, &
                                      sigma_t_2, sigma_s_2, sigma_Ncn_2_n, &
                                      corr_ts_2, corr_tNcn_2_n, corr_sNcn_2_n, &
                                      mu_t_2, KK_auto_tndcy, KK_auto_coef, &
                                      t_tol, alpha_exp, beta_exp )  &
               - bivar_NL_mean_eq( mu_s_2, mu_Ncn_2, mu_Ncn_2_n, sigma_s_2, &
                                   sigma_Ncn_2, sigma_Ncn_2_n, corr_sNcn_2_n, &
                                   Nc_tol, alpha_exp + one, beta_exp )  &
               + mu_s_2  &
                 * bivar_NL_mean_eq( mu_s_2, mu_Ncn_2, mu_Ncn_2_n, sigma_s_2, &
                                     sigma_Ncn_2, sigma_Ncn_2_n, corr_sNcn_2_n,&
                                     Nc_tol, alpha_exp, beta_exp )  &    
               )  &
           )

    endif


    return

  end function covar_thl_KK_auto

  !=============================================================================
  function covar_x_KK_accr( mu_x_1, mu_x_2, mu_s_1, mu_s_2, mu_rr_1_n, &
                            mu_rr_2_n, sigma_x_1, sigma_x_2, sigma_s_1, &
                            sigma_s_2, sigma_rr_1_n, sigma_rr_2_n, &
                            corr_xs_1, corr_xs_2, corr_xrr_1_n, &
                            corr_xrr_2_n, corr_srr_1_n, corr_srr_2_n, &
                            x_mean, KK_accr_tndcy, KK_accr_coef, x_tol, &
                            mixt_frac, precip_frac_1, precip_frac_2 )

    ! Description:
    ! This function calculates the covariance between x and KK accretion
    ! tendency, which can be written as < x'((dr_r/dt)_KKaccr)' >, or more
    ! simply as < x'KK_accr' >.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one  ! Constant(s)

    use parameters_microphys, only: &
        KK_accr_rc_exp, & ! Variable(s)
        KK_accr_rr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x_1,        & ! Mean of x (1st PDF component)                       [-]
      mu_x_2,        & ! Mean of x (2nd PDF component)                       [-]
      mu_s_1,        & ! Mean of s (1st PDF component)                       [-]
      mu_s_2,        & ! Mean of s (2nd PDF component)                       [-]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) in-precip (ip)    [-]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip                [-]
      sigma_x_1,     & ! Standard deviation of x (1st PDF component)         [-]
      sigma_x_2,     & ! Standard deviation of x (2nd PDF component)         [-]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)         [-]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)         [-]
      sigma_rr_1_n,  & ! Standard deviation of ln rr (1st PDF component) ip  [-]
      sigma_rr_2_n,  & ! Standard deviation of ln rr (2nd PDF component) ip  [-]
      corr_xs_1,     & ! Correlation between x and s (1st PDF component)     [-]
      corr_xs_2,     & ! Correlation between x and s (2nd PDF component)     [-]
      corr_xrr_1_n,  & ! Correlation between x and ln rr (1st PDF comp.) ip  [-]
      corr_xrr_2_n,  & ! Correlation between x and ln rr (2nd PDF comp.) ip  [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      x_mean,        & ! Mean of x (overall)                        [units vary]
      KK_accr_tndcy, & ! KK accretion tendency                       [(kg/kg)/s]
      KK_accr_coef,  & ! KK accretion coefficient                    [(kg/kg)/s]
      x_tol,         & ! Tolerance value of x                       [units vary]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)          [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_x_KK_accr  ! Covariance between x and KK accretion tendency      [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on s                                           [-]
      beta_exp     ! Exponent on r_r                                         [-]


    ! Values of the KK exponents.
    alpha_exp = KK_accr_rc_exp
    beta_exp  = KK_accr_rr_exp

    ! Calculate the covariance of x and KK accretion tendency.
    covar_x_KK_accr  &
    = KK_accr_coef &
      * ( mixt_frac &
          * precip_frac_1 &
          * trivar_NNL_covar_eq( mu_x_1, mu_s_1, mu_rr_1_n, &
                                 sigma_x_1, sigma_s_1, sigma_rr_1_n, &
                                 corr_xs_1, corr_xrr_1_n, corr_srr_1_n, &
                                 x_mean, KK_accr_tndcy, KK_accr_coef, &
                                 x_tol, alpha_exp, beta_exp ) &
        + ( one - mixt_frac ) &
          * precip_frac_2 &
          * trivar_NNL_covar_eq( mu_x_2, mu_s_2, mu_rr_2_n, &
                                 sigma_x_2, sigma_s_2, sigma_rr_2_n, &
                                 corr_xs_2, corr_xrr_2_n, corr_srr_2_n, &
                                 x_mean, KK_accr_tndcy, KK_accr_coef, &
                                 x_tol, alpha_exp, beta_exp ) &
        )


    return

  end function covar_x_KK_accr

  !=============================================================================
  function covar_rt_KK_accr( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_rr_1, &
                             mu_rr_2, mu_rr_1_n, mu_rr_2_n, sigma_t_1, &
                             sigma_t_2, sigma_s_1, sigma_s_2, sigma_rr_1, &
                             sigma_rr_2, sigma_rr_1_n, sigma_rr_2_n, &
                             corr_ts_1, corr_ts_2, corr_trr_1_n, &
                             corr_trr_2_n, corr_srr_1_n, corr_srr_2_n, &
                             KK_accr_tndcy, KK_accr_coef, t_tol, crt1, &
                             crt2, mixt_frac, precip_frac_1, precip_frac_2 )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two,    & ! Constant(s)
        one,    &
        rr_tol

    use KK_upscaled_means, only:  &
        bivar_NL_mean_eq

    use parameters_microphys, only: &
        KK_accr_rc_exp, & ! Variable(s)
        KK_accr_rr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_t_1,        & ! Mean of t (1st PDF component)                   [kg/kg]
      mu_t_2,        & ! Mean of t (2nd PDF component)                   [kg/kg]
      mu_s_1,        & ! Mean of s (1st PDF component)                   [kg/kg]
      mu_s_2,        & ! Mean of s (2nd PDF component)                   [kg/kg]
      mu_rr_1,       & ! Mean of rr (1st PDF component) in-precip (ip)   [kg/kg]
      mu_rr_2,       & ! Mean of rr (2nd PDF component) ip               [kg/kg]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) ip                [-]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip                [-]
      sigma_t_1,     & ! Standard deviation of t (1st PDF component)     [kg/kg]
      sigma_t_2,     & ! Standard deviation of t (2nd PDF component)     [kg/kg]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)     [kg/kg]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)     [kg/kg]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_rr_1_n,  & ! Standard deviation of ln rr (1st PDF component) ip  [-]
      sigma_rr_2_n,  & ! Standard deviation of ln rr (2nd PDF component) ip  [-]
      corr_ts_1,     & ! Correlation between t and s (1st PDF component)     [-]
      corr_ts_2,     & ! Correlation between t and s (2nd PDF component)     [-]
      corr_trr_1_n,  & ! Correlation between t and ln rr (1st PDF comp.) ip  [-]
      corr_trr_2_n,  & ! Correlation between t and ln rr (2nd PDF comp.) ip  [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      KK_accr_tndcy, & ! KK accretion tendency                       [(kg/kg)/s]
      KK_accr_coef,  & ! KK accretion coefficient                    [(kg/kg)/s]
      t_tol,         & ! Tolerance value of t                       [units vary]
      crt1,          & ! Coefficient c_rt (1st PDF component)                [-]
      crt2,          & ! Coefficient c_rt (2nd PDF component)                [-]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)          [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_rt_KK_accr  ! Covariance between r_t and KK accretion tendency   [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on s                                           [-]
      beta_exp     ! Exponent on r_r                                         [-]


    ! Values of the KK exponents.
    alpha_exp = KK_accr_rc_exp
    beta_exp  = KK_accr_rr_exp

    ! Calculate the covariance of r_t and KK accretion tendency.
    covar_rt_KK_accr  &
    = KK_accr_coef  &
      * ( mixt_frac * precip_frac_1 * ( one / ( two * crt1 ) )  &
          * ( trivar_NNL_covar_eq( mu_t_1, mu_s_1, mu_rr_1_n, &
                                   sigma_t_1, sigma_s_1, sigma_rr_1_n, &
                                   corr_ts_1, corr_trr_1_n, corr_srr_1_n, &
                                   mu_t_1, KK_accr_tndcy, KK_accr_coef, &
                                   t_tol, alpha_exp, beta_exp )  &
            + bivar_NL_mean_eq( mu_s_1, mu_rr_1, mu_rr_1_n, sigma_s_1, &
                                sigma_rr_1, sigma_rr_1_n, corr_srr_1_n, &
                                rr_tol, alpha_exp + one, beta_exp )  &
            - mu_s_1  &
              * bivar_NL_mean_eq( mu_s_1, mu_rr_1, mu_rr_1_n, sigma_s_1, &
                                  sigma_rr_1, sigma_rr_1_n, corr_srr_1_n, &
                                  rr_tol, alpha_exp, beta_exp )  &
            )  &
        + ( one - mixt_frac ) * precip_frac_2 * ( one / ( two * crt2 ) )  &
          * ( trivar_NNL_covar_eq( mu_t_2, mu_s_2, mu_rr_2_n, &
                                   sigma_t_2, sigma_s_2, sigma_rr_2_n, &
                                   corr_ts_2, corr_trr_2_n, corr_srr_2_n, &
                                   mu_t_2, KK_accr_tndcy, KK_accr_coef, &
                                   t_tol, alpha_exp, beta_exp )  &
            + bivar_NL_mean_eq( mu_s_2, mu_rr_2, mu_rr_2_n, sigma_s_2, &
                                sigma_rr_2, sigma_rr_2_n, corr_srr_2_n, &
                                rr_tol, alpha_exp + one, beta_exp )  &
            - mu_s_2  &
              * bivar_NL_mean_eq( mu_s_2, mu_rr_2, mu_rr_2_n, sigma_s_2, &
                                  sigma_rr_2, sigma_rr_2_n, corr_srr_2_n, &
                                  rr_tol, alpha_exp, beta_exp )  &
            )  &
        )


    return

  end function covar_rt_KK_accr

  !=============================================================================
  function covar_thl_KK_accr( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_rr_1, &
                              mu_rr_2, mu_rr_1_n, mu_rr_2_n, sigma_t_1, &
                              sigma_t_2, sigma_s_1, sigma_s_2, sigma_rr_1, &
                              sigma_rr_2, sigma_rr_1_n, sigma_rr_2_n, &
                              corr_ts_1, corr_ts_2, corr_trr_1_n, &
                              corr_trr_2_n, corr_srr_1_n, corr_srr_2_n, &
                              KK_accr_tndcy, KK_accr_coef, t_tol, cthl1, &
                              cthl2, mixt_frac, precip_frac_1, precip_frac_2 )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two,    & ! Constant(s)
        one,    &
        rr_tol

    use KK_upscaled_means, only:  &
        bivar_NL_mean_eq

    use parameters_microphys, only: &
        KK_accr_rc_exp, & ! Variable(s)
        KK_accr_rr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_t_1,        & ! Mean of t (1st PDF component)                   [kg/kg]
      mu_t_2,        & ! Mean of t (2nd PDF component)                   [kg/kg]
      mu_s_1,        & ! Mean of s (1st PDF component)                   [kg/kg]
      mu_s_2,        & ! Mean of s (2nd PDF component)                   [kg/kg]
      mu_rr_1,       & ! Mean of rr (1st PDF component) in-precip (ip)   [kg/kg]
      mu_rr_2,       & ! Mean of rr (2nd PDF component) ip               [kg/kg]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) ip                [-]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip                [-]
      sigma_t_1,     & ! Standard deviation of t (1st PDF component)     [kg/kg]
      sigma_t_2,     & ! Standard deviation of t (2nd PDF component)     [kg/kg]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)     [kg/kg]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)     [kg/kg]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_rr_1_n,  & ! Standard deviation of ln rr (1st PDF component) ip  [-]
      sigma_rr_2_n,  & ! Standard deviation of ln rr (2nd PDF component) ip  [-]
      corr_ts_1,     & ! Correlation between t and s (1st PDF component)     [-]
      corr_ts_2,     & ! Correlation between t and s (2nd PDF component)     [-]
      corr_trr_1_n,  & ! Correlation between t and ln rr (1st PDF component) [-]
      corr_trr_2_n,  & ! Correlation between t and ln rr (2nd PDF component) [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF component) [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF component) [-]
      KK_accr_tndcy, & ! KK accretion tendency                       [(kg/kg)/s]
      KK_accr_coef,  & ! KK accretion coefficient                    [(kg/kg)/s]
      t_tol,         & ! Tolerance value of t                       [units vary]
      cthl1,         & ! Coefficient c_thl (1st PDF component)               [-]
      cthl2,         & ! Coefficient c_thl (2nd PDF component)               [-]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)          [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_thl_KK_accr  ! Covariance between th_l and KK accretion tendency [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on s                                           [-]
      beta_exp     ! Exponent on r_r                                         [-]


    ! Values of the KK exponents.
    alpha_exp = KK_accr_rc_exp
    beta_exp  = KK_accr_rr_exp

    ! Calculate the covariance of th_l and KK accretion tendency.
    covar_thl_KK_accr  &
    = KK_accr_coef  &
      * ( mixt_frac * precip_frac_1 * ( one / ( two * cthl1 ) )  &
          * ( trivar_NNL_covar_eq( mu_t_1, mu_s_1, mu_rr_1_n, &
                                   sigma_t_1, sigma_s_1, sigma_rr_1_n, &
                                   corr_ts_1, corr_trr_1_n, corr_srr_1_n, &
                                   mu_t_1, KK_accr_tndcy, KK_accr_coef, &
                                   t_tol, alpha_exp, beta_exp )  &
            - bivar_NL_mean_eq( mu_s_1, mu_rr_1, mu_rr_1_n, sigma_s_1, &
                                sigma_rr_1, sigma_rr_1_n, corr_srr_1_n, &
                                rr_tol, alpha_exp + one, beta_exp )  &
            + mu_s_1  &
              * bivar_NL_mean_eq( mu_s_1, mu_rr_1, mu_rr_1_n, sigma_s_1, &
                                  sigma_rr_1, sigma_rr_1_n, corr_srr_1_n, &
                                  rr_tol, alpha_exp, beta_exp )  &
            )  &
        + ( one - mixt_frac ) * precip_frac_2 * ( one / ( two * cthl2 ) )  &
          * ( trivar_NNL_covar_eq( mu_t_2, mu_s_2, mu_rr_2_n, &
                                   sigma_t_2, sigma_s_2, sigma_rr_2_n, &
                                   corr_ts_2, corr_trr_2_n, corr_srr_2_n, &
                                   mu_t_2, KK_accr_tndcy, KK_accr_coef, &
                                   t_tol, alpha_exp, beta_exp )  &
            - bivar_NL_mean_eq( mu_s_2, mu_rr_2, mu_rr_2_n, sigma_s_2, &
                                sigma_rr_2, sigma_rr_2_n, corr_srr_2_n, &
                                rr_tol, alpha_exp + one, beta_exp )  &
            + mu_s_2  &
              * bivar_NL_mean_eq( mu_s_2, mu_rr_2, mu_rr_2_n, sigma_s_2, &
                                  sigma_rr_2, sigma_rr_2_n, corr_srr_2_n, &
                                  rr_tol, alpha_exp, beta_exp )  &
            )  &
        )


    return

  end function covar_thl_KK_accr

  !=============================================================================
  function quadrivar_NNLL_covar_eq( mu_x_i, mu_s_i, mu_rr_i_n, mu_Nr_i_n, &
                                    sigma_x_i, sigma_s_i, sigma_rr_i_n, &
                                    sigma_Nr_i_n, corr_xs_i, corr_xrr_i_n, &
                                    corr_xNr_i_n, corr_srr_i_n, &
                                    corr_sNr_i_n, corr_rrNr_i_n, x_mean, &
                                    mc_tndcy_mean, mc_coef, x_tol, &
                                    alpha_exp_in, beta_exp_in, gamma_exp_in )

    ! Description:
    ! This function calculates the contribution by the ith PDF component to the
    ! expression < y1'y2'_(i) >, where y1 = x1 ( = x, which is w or t), and
    ! where y2 = x2^alpha x3^beta x4^gamma ( = s^alpha r_r^beta N_r^gamma, which
    ! also equals KK_evap_tndcy / KK_evap_coef).  The value of covariance of x
    ! and KK evaporation tendency is:
    !
    ! < x'KK_evap' > = KK_evap_coef
    !                  * ( mixt_frac < y1'y2'_(1) >
    !                      + ( 1 - mixt_frac ) < y1'y2'_(2) > ).
    ! 
    ! One of four functions are called, based on whether x1 and/or x2 (x and/or
    ! s) vary.  Each one of these four functions is the result of an evaluated
    ! integral based on the specific situation.

    ! References:
    !-----------------------------------------------------------------------

    use PDF_integrals_covars, only: &
        quadrivar_NNLL_covar,            & ! Procedure(s)
        quadrivar_NNLL_covar_const_x1,   &
        quadrivar_NNLL_covar_const_x2,   &
        quadrivar_NNLL_covar_const_x1x2

    use constants_clubb, only: &
        s_mellor_tol, & ! Constant(s)
        parab_cyl_max_input

    use clubb_precision, only: &
        dp,        & ! double precision
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x_i,        & ! Mean of x (ith PDF component)                       [-]
      mu_s_i,        & ! Mean of s (ith PDF component)                       [-]
      mu_rr_i_n,     & ! Mean of ln rr (ith PDF component) in-precip (ip)    [-]
      mu_Nr_i_n,     & ! Mean of ln Nr (ith PDF component) ip                [-]
      sigma_x_i,     & ! Standard deviation of x (ith PDF component)         [-]
      sigma_s_i,     & ! Standard deviation of s (ith PDF component)         [-]
      sigma_rr_i_n,  & ! Standard deviation of ln rr (ith PDF component) ip  [-]
      sigma_Nr_i_n,  & ! Standard deviation of ln Nr (ith PDF component) ip  [-]
      corr_xs_i,     & ! Correlation between x and s (ith PDF component)     [-]
      corr_xrr_i_n,  & ! Correlation between x and ln rr (ith PDF comp.) ip  [-]
      corr_xNr_i_n,  & ! Correlation between x and ln Nr (ith PDF comp.) ip  [-]
      corr_srr_i_n,  & ! Correlation between s and ln rr (ith PDF comp.) ip  [-]
      corr_sNr_i_n,  & ! Correlation between s and ln Nr (ith PDF comp.) ip  [-]
      corr_rrNr_i_n    ! Correlation btwn. ln rr & ln Nr (ith PDF comp.) ip  [-]

    real( kind = core_rknd ), intent(in) :: &
      x_mean,        & ! Mean of x (overall)                       [units vary]
      mc_tndcy_mean, & ! Mean of microphysics tendency              [(kg/kg)/s]
      mc_coef,       & ! Coefficient of microphysics tendency       [(kg/kg)/s]
      x_tol            ! Tolerance value of x                      [units vary]

    real( kind = core_rknd ), intent(in) :: &
      alpha_exp_in,  & ! Exponent alpha, corresponding to s                 [-]
      beta_exp_in,   & ! Exponent beta, corresponding to rr                 [-]
      gamma_exp_in     ! Exponent gamma, corresponding to Nr                [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      quadrivar_NNLL_covar_eq

    ! Local Variables
    real( kind = dp ) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      mu_x4_n,    & ! Mean of ln x4 (ith PDF component)                     [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x2,   & ! Standard deviation of x2 (ith PDF component)          [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      sigma_x4_n, & ! Standard deviation of ln x4 (ith PDF component)       [-]
      rho_x1x2,   & ! Correlation between x1 and x2 (ith PDF component)     [-]
      rho_x1x3_n, & ! Correlation between x1 and ln x3 (ith PDF component)  [-]
      rho_x1x4_n, & ! Correlation between x1 and ln x4 (ith PDF component)  [-]
      rho_x2x3_n, & ! Correlation between x2 and ln x3 (ith PDF component)  [-]
      rho_x2x4_n, & ! Correlation between x2 and ln x4 (ith PDF component)  [-]
      rho_x3x4_n    ! Correlation between ln x3 & ln x4 (ith PDF component) [-]

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
      s_cc      ! Parabolic cylinder function input value                   [-]


    ! Means for the ith PDF component. 
    mu_x1   = dble( mu_x_i )  ! x is w or t (ith component).
    mu_x2   = dble( mu_s_i )
    mu_x3_n = dble( mu_rr_i_n )
    mu_x4_n = dble( mu_Nr_i_n )

    ! Standard deviations for the ith PDF component.
    sigma_x1   = dble( sigma_x_i )  ! x is w or t (ith component).
    sigma_x2   = dble( sigma_s_i )
    sigma_x3_n = dble( sigma_rr_i_n )
    sigma_x4_n = dble( sigma_Nr_i_n )

    ! Correlations for the ith PDF component.
    rho_x1x2   = dble( corr_xs_i )    ! x is w or t (ith component).
    rho_x1x3_n = dble( corr_xrr_i_n ) ! x is w or t (ith component).
    rho_x1x4_n = dble( corr_xNr_i_n ) ! x is w or t (ith component).
    rho_x2x3_n = dble( corr_srr_i_n )
    rho_x2x4_n = dble( corr_sNr_i_n )
    rho_x3x4_n = dble( corr_rrNr_i_n )

    ! Overall means.
    x1_mean = dble( x_mean )  ! x is w or t.
    x2_alpha_x3_beta_x4_gamma_mean = dble( mc_tndcy_mean / mc_coef )

    ! Exponents.
    alpha_exp = dble( alpha_exp_in )
    beta_exp  = dble( beta_exp_in )
    gamma_exp = dble( gamma_exp_in )

    ! Tolerance values.
    ! When the standard deviation of a variable is below the tolerance values,
    ! it is considered to be zero, and the variable is considered to have a
    ! constant value.
    x1_tol = dble( x_tol )  ! x is w or t.
    x2_tol = dble( s_mellor_tol )

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


    ! Based on the values of sigma_x1 and sigma_x2 (including the value of s_cc
    ! compared to parab_cyl_max_input), find the correct form of the
    ! quadrivariate equation to use.

    if ( sigma_x1 <= x1_tol .and.  &
         ( sigma_x2 <= x2_tol .or.  &
           abs( s_cc ) > dble( parab_cyl_max_input ) ) ) then

       ! The ith PDF component variance of both x (w or t) and s is 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_const_x1x2( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                          sigma_x3_n, sigma_x4_n, &
                                          rho_x3x4_n, x1_mean, &
                                          x2_alpha_x3_beta_x4_gamma_mean, &
                                          alpha_exp, beta_exp, gamma_exp ),  &
         kind = core_rknd )


    elseif ( sigma_x1 <= x1_tol ) then

       ! The ith PDF component variance of x (w or t) is 0.
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
             abs( s_cc ) > dble( parab_cyl_max_input ) ) then

       ! The ith PDF component variance of s is 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_const_x2( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                        sigma_x1, sigma_x3_n, sigma_x4_n, &
                                        rho_x1x3_n, rho_x1x4_n, rho_x3x4_n, &
                                        x1_mean, &
                                        x2_alpha_x3_beta_x4_gamma_mean, &
                                        alpha_exp, beta_exp, gamma_exp ),  &
         kind = core_rknd )


    else  ! sigma_x1 > 0 and sigma_x2 > 0.

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
  function trivar_NNL_covar_eq( mu_x_i, mu_s_i, mu_y_i_n, &
                                sigma_x_i, sigma_s_i, sigma_y_i_n, &
                                corr_xs_i, corr_xy_i_n, corr_sy_i_n, &
                                x_mean, mc_tndcy_mean, mc_coef, &
                                x_tol, alpha_exp_in, beta_exp_in )

    ! Description:
    ! This function calculates the contribution by the ith PDF component to the
    ! expression < y1'y2'_(i) >, where y1 = x1 ( = x, which is w or t), and
    ! where y2 = x2^alpha x3^beta ( = s^alpha y^beta, where y is N_cn or r_r for
    ! autoconversion or accretion, respectively, and which also equals
    ! KK_auto_tndcy / KK_auto_coef or KK_accr_tndcy / KK_accr_coef,
    ! respectively).  The value of covariance of x and the KK microphysics
    ! tendency is:
    !
    ! < x'KK_mc' > = KK_mc_coef
    !                * ( mixt_frac < y1'y2'_(1) >
    !                    + ( 1 - mixt_frac ) < y1'y2'_(2) > ).
    ! 
    ! One of four functions are called, based on whether x1 and/or x2 (x and/or
    ! s) vary.  Each one of these four functions is the result of an evaluated
    ! integral based on the specific situation.

    ! References:
    !-----------------------------------------------------------------------

    use PDF_integrals_covars, only: &
        trivar_NNL_covar,            & ! Procedure(s)
        trivar_NNL_covar_const_x1,   &
        trivar_NNL_covar_const_x2,   &
        trivar_NNL_covar_const_x1x2

    use constants_clubb, only: &
        s_mellor_tol, & ! Constant(s)
        parab_cyl_max_input

    use clubb_precision, only: &
        dp,        & ! double precision
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x_i,      & ! Mean of x (ith PDF component)                        [-]
      mu_s_i,      & ! Mean of s (ith PDF component)                        [-]
      mu_y_i_n,    & ! Mean of ln y (ith PDF component)                     [-]
      sigma_x_i,   & ! Standard deviation of x (ith PDF component)          [-]
      sigma_s_i,   & ! Standard deviation of s (ith PDF component)          [-]
      sigma_y_i_n, & ! Standard deviation of ln y (ith PDF component)       [-]
      corr_xs_i,   & ! Correlation between x and s (ith PDF component)      [-]
      corr_xy_i_n, & ! Correlation between x and ln y (ith PDF component)   [-]
      corr_sy_i_n    ! Correlation between s and ln y (ith PDF component)   [-]

    real( kind = core_rknd ), intent(in) :: &
      x_mean,        & ! Mean of x (overall)                       [units vary]
      mc_tndcy_mean, & ! Mean of microphysics tendency              [(kg/kg)/s]
      mc_coef,       & ! Coefficient of microphysics                [(kg/kg)/s]
      x_tol            ! Tolerance value of x                      [units vary]

    real( kind = core_rknd ), intent(in) :: &
      alpha_exp_in,  & ! Exponent alpha, corresponding to s                 [-]
      beta_exp_in      ! Exponent beta, corresponding to y                  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      trivar_NNL_covar_eq

    ! Local Variables
    real( kind = dp ) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x2,   & ! Standard deviation of x2 (ith PDF component)          [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      rho_x1x2,   & ! Correlation between x1 and x2 (ith PDF component)     [-]
      rho_x1x3_n, & ! Correlation between x1 and ln x3 (ith PDF component)  [-]
      rho_x2x3_n    ! Correlation between x2 and ln x3 (ith PDF component)  [-]

    real( kind = dp ) :: &
      x1_mean,               & ! Mean of x1 (overall)                       [-]
      x2_alpha_x3_beta_mean    ! Mean of x2^alpha x3^beta                   [-]
    
    real( kind = dp ) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp      ! Exponent beta, corresponding to x3                    [-]

    real( kind = dp ) :: &
      x1_tol, & ! Tolerance value of x1                                     [-]
      x2_tol, & ! Tolerance value of x2                                     [-]
      s_c       ! Parabolic cylinder function input value                   [-]


    ! Means for the ith PDF component. 
    mu_x1   = dble( mu_x_i ) ! x is w or t (ith component).
    mu_x2   = dble( mu_s_i )
    mu_x3_n = dble( mu_y_i_n ) ! y is N_cn (autoconversion) or r_r (accretion).

    ! Standard deviations for the ith PDF component.
    sigma_x1   = dble( sigma_x_i ) ! x is w or t (ith component).
    sigma_x2   = dble( sigma_s_i )
    sigma_x3_n = dble( sigma_y_i_n ) ! y is N_cn (auto.) or r_r (accr.).

    ! Correlations for the ith PDF component.
    rho_x1x2   = dble( corr_xs_i )    ! x is w or t (ith component).
    rho_x1x3_n = dble( corr_xy_i_n )  ! x is w or t (ith component).
                                      ! y is N_cn (auto.) or r_r (accr.).
    rho_x2x3_n = dble( corr_sy_i_n )  ! y is N_cn (auto.) or r_r (accr.).

    ! Overall means.
    x1_mean = dble( x_mean )  ! x is w or t.
    x2_alpha_x3_beta_mean = dble( mc_tndcy_mean / mc_coef )

    ! Exponents.
    alpha_exp = dble( alpha_exp_in )
    beta_exp  = dble( beta_exp_in )

    ! Tolerance values.
    ! When the standard deviation of a variable is below the tolerance values,
    ! it is considered to be zero, and the variable is considered to have a
    ! constant value.
    x1_tol = dble( x_tol )  ! x is w or t.
    x2_tol = dble( s_mellor_tol )

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


    ! Based on the values of sigma_x1 and sigma_x2 (including the value of s_c
    ! compared to parab_cyl_max_input), find the correct form of the trivariate
    ! equation to use.

    if ( sigma_x1 <= x1_tol .and.  &
         ( sigma_x2 <= x2_tol .or.  &
           abs( s_c ) > dble( parab_cyl_max_input ) ) ) then

       ! The ith PDF component variance of both x (w or t) and s is 0.
       trivar_NNL_covar_eq  &
       = real( &
         trivar_NNL_covar_const_x1x2( mu_x1, mu_x2, mu_x3_n, sigma_x3_n, &
                                      x1_mean, x2_alpha_x3_beta_mean, &
                                      alpha_exp, beta_exp ),  &
         kind = core_rknd )


    elseif ( sigma_x1 <= x1_tol ) then

       ! The ith PDF component variance of x (w or t) is 0.
       trivar_NNL_covar_eq  &
       = real( trivar_NNL_covar_const_x1( mu_x1, mu_x2, mu_x3_n, &
                                          sigma_x2, sigma_x3_n, rho_x2x3_n, &
                                          x1_mean, x2_alpha_x3_beta_mean, &
                                          alpha_exp, beta_exp ),  &
               kind = core_rknd )


    elseif ( sigma_x2 <= x2_tol .or.  &
             abs( s_c ) > dble( parab_cyl_max_input ) ) then

       ! The ith PDF component variance of s is 0.
       trivar_NNL_covar_eq  &
       = real( trivar_NNL_covar_const_x2( mu_x1, mu_x2, mu_x3_n, &
                                          sigma_x1, sigma_x3_n, rho_x1x3_n, &
                                          x1_mean, x2_alpha_x3_beta_mean, &
                                          alpha_exp, beta_exp ),  &
               kind = core_rknd )


    else  ! sigma_x1 > 0 and sigma_x2 > 0.

       ! This is the complete value of the trivariate.
       ! All fields vary in the ith PDF component.
       trivar_NNL_covar_eq  &
       = real( trivar_NNL_covar( mu_x1, mu_x2, mu_x3_n, &
                                 sigma_x1, sigma_x2, sigma_x3_n, &
                                 rho_x1x2, rho_x1x3_n, rho_x2x3_n, &
                                 x1_mean, x2_alpha_x3_beta_mean, &
                                 alpha_exp, beta_exp ), kind = core_rknd )


    endif


    return

  end function trivar_NNL_covar_eq

  !=============================================================================
  function trivar_NNL_covar_eq_Nc0( mu_x_i, mu_s_i, Nc0_in_cloud, &
                                    sigma_x_i, sigma_s_i, corr_xs_i, &
                                    x_mean, mc_tndcy_mean, mc_coef, &
                                    x_tol, alpha_exp_in, beta_exp_in )

    ! Description:
    ! This function calculates the contribution by the ith PDF component to the
    ! expression < y1'y2'_(i) >, where y1 = x1 ( = x, which is w or t), and
    ! where y2 = x2^alpha x3^beta ( = s^alpha y^beta, where y is N_c or r_r for
    ! autoconversion or accretion, respectively, and which also equals
    ! KK_auto_tndcy / KK_auto_coef or KK_accr_tndcy / KK_accr_coef,
    ! respectively).  The value of covariance of x and the KK microphysics
    ! tendency is:
    !
    ! < x'KK_mc' > = KK_mc_coef
    !                * ( mixt_frac < y1'y2'_(1) >
    !                    + ( 1 - mixt_frac ) < y1'y2'_(2) > ).
    ! 
    ! One of four functions are called, based on whether x1 and/or x2 (x and/or
    ! s) vary.  Each one of these four functions is the result of an evaluated
    ! integral based on the specific situation.

    ! References:
    !-----------------------------------------------------------------------

    use PDF_integrals_covars, only: &
        trivar_NNL_covar_const_x3,     & ! Procedure(s)
        trivar_NNL_covar_const_x1x3,   &
        trivar_NNL_covar_const_x2x3,   &
        trivar_NNL_covar_const_all

    use constants_clubb, only: &
        s_mellor_tol, & ! Constant(s)
        parab_cyl_max_input

    use clubb_precision, only: &
        dp,        & ! double precision
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x_i,       & ! Mean of x (ith PDF component)                       [-]
      mu_s_i,       & ! Mean of s (ith PDF component)                   [kg/kg]
      Nc0_in_cloud, & ! Constant in-cloud value of cloud droplet conc. [num/kg]
      sigma_x_i,    & ! Standard deviation of x (ith PDF component)         [-]
      sigma_s_i,    & ! Standard deviation of s (ith PDF component)     [kg/kg]
      corr_xs_i       ! Correlation between x and s (ith PDF component)     [-]

    real( kind = core_rknd ), intent(in) :: &
      x_mean,        & ! Mean of x (overall)                       [units vary]
      mc_tndcy_mean, & ! Mean of microphysics tendency              [(kg/kg)/s]
      mc_coef,       & ! Coefficient of microphysics                [(kg/kg)/s]
      x_tol            ! Tolerance value of x                      [units vary]

    real( kind = core_rknd ), intent(in) :: &
      alpha_exp_in,  & ! Exponent alpha, corresponding to s                 [-]
      beta_exp_in      ! Exponent beta, corresponding to y                  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      trivar_NNL_covar_eq_Nc0

    ! Local Variables
    real( kind = dp ) :: &
      mu_x1,    & ! Mean of x1 (ith PDF component)                          [-]
      mu_x2,    & ! Mean of x2 (ith PDF component)                          [-]
      Nc0,      & ! Constant in-cloud value of cloud droplet conc.     [num/kg]
      sigma_x1, & ! Standard deviation of x1 (ith PDF component)            [-]
      sigma_x2, & ! Standard deviation of x2 (ith PDF component)            [-]
      rho_x1x2    ! Correlation between x1 and x2 (ith PDF component)       [-]

    real( kind = dp ) :: &
      x1_mean,               & ! Mean of x1 (overall)                       [-]
      x2_alpha_x3_beta_mean    ! Mean of x2^alpha x3^beta                   [-]
    
    real( kind = dp ) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp      ! Exponent beta, corresponding to x3                    [-]

    real( kind = dp ) :: &
      x1_tol, & ! Tolerance value of x1                                     [-]
      x2_tol, & ! Tolerance value of x2                                     [-]
      s_c       ! Parabolic cylinder function input value                   [-]


    ! Means for the ith PDF component. 
    mu_x1 = dble( mu_x_i ) ! x is w or t (ith component).
    mu_x2 = dble( mu_s_i )
    Nc0   = dble( Nc0_in_cloud )

    ! Standard deviations for the ith PDF component.
    sigma_x1 = dble( sigma_x_i ) ! x is w or t (ith component).
    sigma_x2 = dble( sigma_s_i )

    ! Correlations for the ith PDF component.
    rho_x1x2 = dble( corr_xs_i )    ! x is w or t (ith component).

    ! Overall means.
    x1_mean = dble( x_mean )  ! x is w or t.
    x2_alpha_x3_beta_mean = dble( mc_tndcy_mean / mc_coef )

    ! Exponents.
    alpha_exp = dble( alpha_exp_in )
    beta_exp  = dble( beta_exp_in )

    ! Tolerance values.
    ! When the standard deviation of a variable is below the tolerance values,
    ! it is considered to be zero, and the variable is considered to have a
    ! constant value.
    x1_tol = dble( x_tol )  ! x is w or t.
    x2_tol = dble( s_mellor_tol )

    ! Determine the value of the parabolic cylinder function input value, s_c.
    ! The value s_c is being fed into the parabolic cylinder function.  When
    ! the value of s_c is too large in magnitude (depending on the order of the
    ! parabolic cylinder function), overflow occurs, and the output of the
    ! parabolic cylinder function is +/-Inf.  This is primarily due to a large
    ! ratio of mu_x2 to sigma_x2.  When the value of s_c is very large, the
    ! distribution of x2 is basically a spike near the mean, so x2 is treated as
    ! a constant.
    if ( sigma_x2 > x2_tol ) then
       s_c = mu_x2 / sigma_x2
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


    ! Based on the values of sigma_x1 and sigma_x2 (including the value of s_c
    ! compared to parab_cyl_max_input), find the correct form of the trivariate
    ! equation to use.

    if ( sigma_x1 <= x1_tol .and.  &
         ( sigma_x2 <= x2_tol .or.  &
           abs( s_c ) > dble( parab_cyl_max_input ) ) ) then

       ! The ith PDF component variance of both x (w or t) and s is 0.
       trivar_NNL_covar_eq_Nc0  &
       = real( trivar_NNL_covar_const_all( mu_x1, mu_x2, Nc0, &
                                           x1_mean, x2_alpha_x3_beta_mean, &
                                           alpha_exp, beta_exp ),  &
         kind = core_rknd )


    elseif ( sigma_x1 <= x1_tol ) then

       ! The ith PDF component variance of x (w or t) is 0.
       trivar_NNL_covar_eq_Nc0  &
       = real( trivar_NNL_covar_const_x1x3( mu_x1, mu_x2, Nc0, sigma_x2, &
                                            x1_mean, x2_alpha_x3_beta_mean, &
                                            alpha_exp, beta_exp ),  &
               kind = core_rknd )


    elseif ( sigma_x2 <= x2_tol .or.  &
             abs( s_c ) > dble( parab_cyl_max_input ) ) then

       ! The ith PDF component variance of s is 0.
       trivar_NNL_covar_eq_Nc0  &
       = real( trivar_NNL_covar_const_x2x3( mu_x1, mu_x2, Nc0, &
                                            x1_mean, x2_alpha_x3_beta_mean, &
                                            alpha_exp, beta_exp ),  &
               kind = core_rknd )


    else  ! sigma_x1 > 0 and sigma_x2 > 0.

       ! This is the complete value of the trivariate.
       ! All fields vary in the ith PDF component.
       trivar_NNL_covar_eq_Nc0  &
       = real( trivar_NNL_covar_const_x3( mu_x1, mu_x2, Nc0, &
                                          sigma_x1, sigma_x2, rho_x1x2, &
                                          x1_mean, x2_alpha_x3_beta_mean, &
                                          alpha_exp, beta_exp ),  &
               kind = core_rknd )


    endif


    return

  end function trivar_NNL_covar_eq_Nc0

!===============================================================================

end module KK_upscaled_covariances
