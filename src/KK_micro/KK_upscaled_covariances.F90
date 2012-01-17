! $Id$
!===============================================================================
module KK_upscaled_covariances

  implicit none

  private  ! Default scope

  public :: covar_x_KK_evap, &
            covar_rt_KK_evap, &
            covar_thl_KK_evap, &
            covar_x_KK_auto, &
            covar_rt_KK_auto, &
            covar_thl_KK_auto, &
            covar_x_KK_accr, &
            covar_rt_KK_accr, &
            covar_thl_KK_accr

  private :: quadrivar_NNLL_covar_eq, &
             trivar_NNL_covar_eq

  contains

  !=============================================================================
  function covar_x_KK_evap( mu_x_1, mu_x_2, mu_s_1, mu_s_2, mu_rr_n, &
                            mu_Nr_n, sigma_x_1, sigma_x_2, sigma_s_1, &
                            sigma_s_2, sigma_rr_n, sigma_Nr_n, &
                            corr_xs_1, corr_xs_2, corr_xrr_1_n, &
                            corr_xrr_2_n, corr_xNr_1_n, corr_xNr_2_n, &
                            corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                            corr_sNr_2_n, corr_rrNr_n, x_mean, &
                            KK_evap_tndcy, KK_evap_coef, x_tol, mixt_frac )

    ! Description:
    ! This function calculates the covariance between x and KK evaporation
    ! tendency, which can be written as < x'((dr_r/dt)_KKevap)' >, or more
    ! simply as < x'KK_evap' >.

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) :: &
      mu_x_1,        & ! Mean of x (1st PDF component)                       [-]
      mu_x_2,        & ! Mean of x (2nd PDF component)                       [-]
      mu_s_1,        & ! Mean of s (1st PDF component)                       [-]
      mu_s_2,        & ! Mean of s (2nd PDF component)                       [-]
      mu_rr_n,       & ! Mean of ln rr (both components)                     [-]
      mu_Nr_n,       & ! Mean of ln Nr (both components)                     [-]
      sigma_x_1,     & ! Standard deviation of x (1st PDF component)         [-]
      sigma_x_2,     & ! Standard deviation of x (2nd PDF component)         [-]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)         [-]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)         [-]
      sigma_rr_n,    & ! Standard deviation of ln rr (both components)       [-]
      sigma_Nr_n,    & ! Standard deviation of ln Nr (both components)       [-]
      corr_xs_1,     & ! Correlation between x and s (1st PDF component)     [-]
      corr_xs_2,     & ! Correlation between x and s (2nd PDF component)     [-]
      corr_xrr_1_n,  & ! Correlation between x and ln rr (1st PDF component) [-]
      corr_xrr_2_n,  & ! Correlation between x and ln rr (2nd PDF component) [-]
      corr_xNr_1_n,  & ! Correlation between x and ln Nr (1st PDF component) [-]
      corr_xNr_2_n,  & ! Correlation between x and ln Nr (2nd PDF component) [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF component) [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF component) [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF component) [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF component) [-]
      corr_rrNr_n,   & ! Correlation between ln rr & ln Nr (both components) [-]
      x_mean,        & ! Mean of x (overall)                        [units vary]
      KK_evap_tndcy, & ! KK evaporation tendency                     [(kg/kg)/s]
      KK_evap_coef,  & ! KK evaporation coefficient                  [(kg/kg)/s]
      x_tol,         & ! Tolerance value of x                       [units vary]
      mixt_frac        ! Mixture fraction                                    [-]

    ! Return Variable
    real :: &
      covar_x_KK_evap  ! Covariance between x and KK evaporation tendency    [-]

    ! Constant Parameters
    real, parameter :: &
      alpha_exp = 1.0,       & ! Exponent on s                               [-]
      beta_exp  = (1.0/3.0), & ! Exponent on r_r                             [-]
      gamma_exp = (2.0/3.0)    ! Exponent on N_r                             [-]


    ! Calculate the covariance of x and KK evaporation tendency.
    covar_x_KK_evap  &
    = KK_evap_coef  &
      * ( mixt_frac  &
          * quadrivar_NNLL_covar_eq( mu_x_1, mu_s_1, mu_rr_n, mu_Nr_n, &
                                     sigma_x_1, sigma_s_1, sigma_rr_n, &
                                     sigma_Nr_n, corr_xs_1, corr_xrr_1_n, &
                                     corr_xNr_1_n, corr_srr_1_n, &
                                     corr_sNr_1_n, corr_rrNr_n, x_mean, &
                                     KK_evap_tndcy, KK_evap_coef, x_tol, &
                                     alpha_exp, beta_exp, gamma_exp ) &
        + ( 1.0 - mixt_frac ) &
          * quadrivar_NNLL_covar_eq( mu_x_2, mu_s_2, mu_rr_n, mu_Nr_n, &
                                     sigma_x_2, sigma_s_2, sigma_rr_n, &
                                     sigma_Nr_n, corr_xs_2, corr_xrr_2_n, &
                                     corr_xNr_2_n, corr_srr_2_n, &
                                     corr_sNr_2_n, corr_rrNr_n, x_mean, &
                                     KK_evap_tndcy, KK_evap_coef, x_tol, &
                                     alpha_exp, beta_exp, gamma_exp ) &
        )


    return

  end function covar_x_KK_evap

  !=============================================================================
  function covar_rt_KK_evap( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_rr_n, &
                             mu_Nr_n, sigma_t_1, sigma_t_2, sigma_s_1, &
                             sigma_s_2, sigma_rr_n, sigma_Nr_n, &
                             corr_ts_1, corr_ts_2, corr_trr_1_n, &
                             corr_trr_2_n, corr_tNr_1_n, corr_tNr_2_n, &
                             corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                             corr_sNr_2_n, corr_rrNr_n, KK_evap_tndcy, &
                             KK_evap_coef, t_tol, crt1, crt2, mixt_frac )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use KK_upscaled_means, only:  &
        trivar_NLL_mean_eq  ! Procedure

    implicit none

    ! Input Variables
    real, intent(in) :: &
      mu_t_1,        & ! Mean of t (1st PDF component)                       [-]
      mu_t_2,        & ! Mean of t (2nd PDF component)                       [-]
      mu_s_1,        & ! Mean of s (1st PDF component)                       [-]
      mu_s_2,        & ! Mean of s (2nd PDF component)                       [-]
      mu_rr_n,       & ! Mean of ln rr (both components)                     [-]
      mu_Nr_n,       & ! Mean of ln Nr (both components)                     [-]
      sigma_t_1,     & ! Standard deviation of t (1st PDF component)         [-]
      sigma_t_2,     & ! Standard deviation of t (2nd PDF component)         [-]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)         [-]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)         [-]
      sigma_rr_n,    & ! Standard deviation of ln rr (both components)       [-]
      sigma_Nr_n,    & ! Standard deviation of ln Nr (both components)       [-]
      corr_ts_1,     & ! Correlation between t and s (1st PDF component)     [-]
      corr_ts_2,     & ! Correlation between t and s (2nd PDF component)     [-]
      corr_trr_1_n,  & ! Correlation between t and ln rr (1st PDF component) [-]
      corr_trr_2_n,  & ! Correlation between t and ln rr (2nd PDF component) [-]
      corr_tNr_1_n,  & ! Correlation between t and ln Nr (1st PDF component) [-]
      corr_tNr_2_n,  & ! Correlation between t and ln Nr (2nd PDF component) [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF component) [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF component) [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF component) [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF component) [-]
      corr_rrNr_n,   & ! Correlation between ln rr & ln Nr (both components) [-]
      KK_evap_tndcy, & ! KK evaporation tendency                     [(kg/kg)/s]
      KK_evap_coef,  & ! KK evaporation coefficient                  [(kg/kg)/s]
      t_tol,         & ! Tolerance value of t                                [-]
      crt1,          & ! Coefficient c_rt (1st PDF component)                [-]
      crt2,          & ! Coefficient c_rt (2nd PDF component)                [-]
      mixt_frac        ! Mixture fraction                                    [-]

    ! Return Variable
    real :: &
      covar_rt_KK_evap  ! Covariance between r_t and KK evaporation tendency [-]

    ! Constant Parameters
    real, parameter :: &
      alpha_exp = 1.0,       & ! Exponent on s                               [-]
      beta_exp  = (1.0/3.0), & ! Exponent on r_r                             [-]
      gamma_exp = (2.0/3.0)    ! Exponent on N_r                             [-]


    ! Calculate the covariance of r_t and KK evaporation tendency.
    covar_rt_KK_evap  &
    = KK_evap_coef  &
      * ( mixt_frac * ( 1.0 / ( 2.0 * crt1 ) )  &
          * ( quadrivar_NNLL_covar_eq( mu_t_1, mu_s_1, mu_rr_n, mu_Nr_n, &
                                       sigma_t_1, sigma_s_1, sigma_rr_n, &
                                       sigma_Nr_n, corr_ts_1, corr_trr_1_n, &
                                       corr_tNr_1_n, corr_srr_1_n, &
                                       corr_sNr_1_n, corr_rrNr_n, mu_t_1, &
                                       KK_evap_tndcy, KK_evap_coef, t_tol, &
                                       alpha_exp, beta_exp, gamma_exp )  &
            + trivar_NLL_mean_eq( mu_s_1, mu_rr_n, mu_Nr_n, &
                                  sigma_s_1, sigma_rr_n, sigma_Nr_n, &
                                  corr_srr_1_n, corr_sNr_1_n, corr_rrNr_n, &
                                  alpha_exp + 1.0, beta_exp, gamma_exp )  &
            - mu_s_1  &
              * trivar_NLL_mean_eq( mu_s_1, mu_rr_n, mu_Nr_n, &
                                    sigma_s_1, sigma_rr_n, sigma_Nr_n, &
                                    corr_srr_1_n, corr_sNr_1_n, corr_rrNr_n, &
                                    alpha_exp, beta_exp, gamma_exp )  &
            )  &
        + ( 1.0 - mixt_frac ) * ( 1.0 / ( 2.0 * crt2 ) )  &
          * ( quadrivar_NNLL_covar_eq( mu_t_2, mu_s_2, mu_rr_n, mu_Nr_n, &
                                       sigma_t_2, sigma_s_2, sigma_rr_n, &
                                       sigma_Nr_n, corr_ts_2, corr_trr_2_n, &
                                       corr_tNr_2_n, corr_srr_2_n, &
                                       corr_sNr_2_n, corr_rrNr_n, mu_t_2, &
                                       KK_evap_tndcy, KK_evap_coef, t_tol, &
                                       alpha_exp, beta_exp, gamma_exp )  &
            + trivar_NLL_mean_eq( mu_s_2, mu_rr_n, mu_Nr_n, &
                                  sigma_s_2, sigma_rr_n, sigma_Nr_n, &
                                  corr_srr_2_n, corr_sNr_2_n, corr_rrNr_n, &
                                  alpha_exp + 1.0, beta_exp, gamma_exp )  &
            - mu_s_2  &
              * trivar_NLL_mean_eq( mu_s_2, mu_rr_n, mu_Nr_n, &
                                    sigma_s_2, sigma_rr_n, sigma_Nr_n, &
                                    corr_srr_2_n, corr_sNr_2_n, corr_rrNr_n, &
                                    alpha_exp, beta_exp, gamma_exp )  &
            )  &
        )


    return

  end function covar_rt_KK_evap

  !=============================================================================
  function covar_thl_KK_evap( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_rr_n, &
                              mu_Nr_n, sigma_t_1, sigma_t_2, sigma_s_1, &
                              sigma_s_2, sigma_rr_n, sigma_Nr_n, &
                              corr_ts_1, corr_ts_2, corr_trr_1_n, &
                              corr_trr_2_n, corr_tNr_1_n, corr_tNr_2_n, &
                              corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                              corr_sNr_2_n, corr_rrNr_n, KK_evap_tndcy, &
                              KK_evap_coef, t_tol, cthl1, cthl2, mixt_frac )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use KK_upscaled_means, only:  &
        trivar_NLL_mean_eq  ! Procedure

    implicit none

    ! Input Variables
    real, intent(in) :: &
      mu_t_1,        & ! Mean of t (1st PDF component)                       [-]
      mu_t_2,        & ! Mean of t (2nd PDF component)                       [-]
      mu_s_1,        & ! Mean of s (1st PDF component)                       [-]
      mu_s_2,        & ! Mean of s (2nd PDF component)                       [-]
      mu_rr_n,       & ! Mean of ln rr (both components)                     [-]
      mu_Nr_n,       & ! Mean of ln Nr (both components)                     [-]
      sigma_t_1,     & ! Standard deviation of t (1st PDF component)         [-]
      sigma_t_2,     & ! Standard deviation of t (2nd PDF component)         [-]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)         [-]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)         [-]
      sigma_rr_n,    & ! Standard deviation of ln rr (both components)       [-]
      sigma_Nr_n,    & ! Standard deviation of ln Nr (both components)       [-]
      corr_ts_1,     & ! Correlation between t and s (1st PDF component)     [-]
      corr_ts_2,     & ! Correlation between t and s (2nd PDF component)     [-]
      corr_trr_1_n,  & ! Correlation between t and ln rr (1st PDF component) [-]
      corr_trr_2_n,  & ! Correlation between t and ln rr (2nd PDF component) [-]
      corr_tNr_1_n,  & ! Correlation between t and ln Nr (1st PDF component) [-]
      corr_tNr_2_n,  & ! Correlation between t and ln Nr (2nd PDF component) [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF component) [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF component) [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF component) [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF component) [-]
      corr_rrNr_n,   & ! Correlation between ln rr & ln Nr (both components) [-]
      KK_evap_tndcy, & ! KK evaporation tendency                     [(kg/kg)/s]
      KK_evap_coef,  & ! KK evaporation coefficient                  [(kg/kg)/s]
      t_tol,         & ! Tolerance value of t                                [-]
      cthl1,         & ! Coefficient c_thl (1st PDF component)               [-]
      cthl2,         & ! Coefficient c_thl (2nd PDF component)               [-]
      mixt_frac        ! Mixture fraction                                    [-]

    ! Return Variable
    real :: &
      covar_thl_KK_evap  ! Covariance between th_l and KK evap. tendency     [-]

    ! Constant Parameters
    real, parameter :: &
      alpha_exp = 1.0,       & ! Exponent on s                               [-]
      beta_exp  = (1.0/3.0), & ! Exponent on r_r                             [-]
      gamma_exp = (2.0/3.0)    ! Exponent on N_r                             [-]


    ! Calculate the covariance of th_l and KK evaporation tendency.
    covar_thl_KK_evap  &
    = KK_evap_coef  &
      * ( mixt_frac * ( 1.0 / ( 2.0 * cthl1 ) )  &
          * ( quadrivar_NNLL_covar_eq( mu_t_1, mu_s_1, mu_rr_n, mu_Nr_n, &
                                       sigma_t_1, sigma_s_1, sigma_rr_n, &
                                       sigma_Nr_n, corr_ts_1, corr_trr_1_n, &
                                       corr_tNr_1_n, corr_srr_1_n, &
                                       corr_sNr_1_n, corr_rrNr_n, mu_t_1, &
                                       KK_evap_tndcy, KK_evap_coef, t_tol, &
                                       alpha_exp, beta_exp, gamma_exp )  &
            - trivar_NLL_mean_eq( mu_s_1, mu_rr_n, mu_Nr_n, &
                                  sigma_s_1, sigma_rr_n, sigma_Nr_n, &
                                  corr_srr_1_n, corr_sNr_1_n, corr_rrNr_n, &
                                  alpha_exp + 1.0, beta_exp, gamma_exp )  &
            + mu_s_1  &
              * trivar_NLL_mean_eq( mu_s_1, mu_rr_n, mu_Nr_n, &
                                    sigma_s_1, sigma_rr_n, sigma_Nr_n, &
                                    corr_srr_1_n, corr_sNr_1_n, corr_rrNr_n, &
                                    alpha_exp, beta_exp, gamma_exp )  &
            )  &
        + ( 1.0 - mixt_frac ) * ( 1.0 / ( 2.0 * cthl2 ) )  &
          * ( quadrivar_NNLL_covar_eq( mu_t_2, mu_s_2, mu_rr_n, mu_Nr_n, &
                                       sigma_t_2, sigma_s_2, sigma_rr_n, &
                                       sigma_Nr_n, corr_ts_2, corr_trr_2_n, &
                                       corr_tNr_2_n, corr_srr_2_n, &
                                       corr_sNr_2_n, corr_rrNr_n, mu_t_2, &
                                       KK_evap_tndcy, KK_evap_coef, t_tol, &
                                       alpha_exp, beta_exp, gamma_exp )  &
            - trivar_NLL_mean_eq( mu_s_2, mu_rr_n, mu_Nr_n, &
                                  sigma_s_2, sigma_rr_n, sigma_Nr_n, &
                                  corr_srr_2_n, corr_sNr_2_n, corr_rrNr_n, &
                                  alpha_exp + 1.0, beta_exp, gamma_exp )  &
            + mu_s_2  &
              * trivar_NLL_mean_eq( mu_s_2, mu_rr_n, mu_Nr_n, &
                                    sigma_s_2, sigma_rr_n, sigma_Nr_n, &
                                    corr_srr_2_n, corr_sNr_2_n, corr_rrNr_n, &
                                    alpha_exp, beta_exp, gamma_exp )  &
            )  &
        )


    return

  end function covar_thl_KK_evap

  !=============================================================================
  function covar_x_KK_auto( mu_x_1, mu_x_2, mu_s_1, mu_s_2, mu_Nc_n, &
                            sigma_x_1, sigma_x_2, sigma_s_1, sigma_s_2, &
                            sigma_Nc_n, corr_xs_1, corr_xs_2, &
                            corr_xNc_1_n, corr_xNc_2_n, corr_sNc_1_n, &
                            corr_sNc_2_n, x_mean, KK_auto_tndcy, &
                            KK_auto_coef, x_tol, mixt_frac )

    ! Description:
    ! This function calculates the correlation between x and KK autoconversion
    ! tendency, which can be written as < x'((dr_r/dt)_KKauto)' >, or more
    ! simply as < x'KK_auto' >.

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) :: &
      mu_x_1,        & ! Mean of x (1st PDF component)                       [-]
      mu_x_2,        & ! Mean of x (2nd PDF component)                       [-]
      mu_s_1,        & ! Mean of s (1st PDF component)                       [-]
      mu_s_2,        & ! Mean of s (2nd PDF component)                       [-]
      mu_Nc_n,       & ! Mean of ln Nc (both components)                     [-]
      sigma_x_1,     & ! Standard deviation of x (1st PDF component)         [-]
      sigma_x_2,     & ! Standard deviation of x (2nd PDF component)         [-]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)         [-]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)         [-]
      sigma_Nc_n,    & ! Standard deviation of ln Nc (both components)       [-]
      corr_xs_1,     & ! Correlation between x and s (1st PDF component)     [-]
      corr_xs_2,     & ! Correlation between x and s (2nd PDF component)     [-]
      corr_xNc_1_n,  & ! Correlation between x and ln Nc (1st PDF component) [-]
      corr_xNc_2_n,  & ! Correlation between x and ln Nc (2nd PDF component) [-]
      corr_sNc_1_n,  & ! Correlation between s and ln Nc (1st PDF component) [-]
      corr_sNc_2_n,  & ! Correlation between s and ln Nc (2nd PDF component) [-]
      x_mean,        & ! Mean of x (overall)                        [units vary]
      KK_auto_tndcy, & ! KK autoconversion tendency                  [(kg/kg)/s]
      KK_auto_coef,  & ! KK autoconversion coefficient               [(kg/kg)/s]
      x_tol,         & ! Tolerance value of x                       [units vary]
      mixt_frac        ! Mixture fraction                                    [-]

    ! Return Variable
    real :: &
      covar_x_KK_auto  ! Covariance between x and KK autoconversion tendency [-]

    ! Constant Parameters
    real, parameter :: &
      alpha_exp = 2.47,  & ! Exponent on s                                   [-]
      beta_exp  = -1.79    ! Exponent on r_r                                 [-]


    ! Calculate the covariance of x and KK autoconversion tendency.
    covar_x_KK_auto  &
    = KK_auto_coef &
      * ( mixt_frac &
          * trivar_NNL_covar_eq( mu_x_1, mu_s_1, mu_Nc_n, &
                                 sigma_x_1, sigma_s_1, sigma_Nc_n, &
                                 corr_xs_1, corr_xNc_1_n, corr_sNc_1_n, &
                                 x_mean, KK_auto_tndcy, KK_auto_coef, &
                                 x_tol, alpha_exp, beta_exp ) &
        + ( 1.0 - mixt_frac ) &
          * trivar_NNL_covar_eq( mu_x_2, mu_s_2, mu_Nc_n, &
                                 sigma_x_2, sigma_s_2, sigma_Nc_n, &
                                 corr_xs_2, corr_xNc_2_n, corr_sNc_2_n, &
                                 x_mean, KK_auto_tndcy, KK_auto_coef, &
                                 x_tol, alpha_exp, beta_exp ) &
        )


    return

  end function covar_x_KK_auto

  !=============================================================================
  function covar_rt_KK_auto( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_Nc_n, &
                             sigma_t_1, sigma_t_2, sigma_s_1, sigma_s_2, &
                             sigma_Nc_n, corr_ts_1, corr_ts_2, &
                             corr_tNc_1_n, corr_tNc_2_n, corr_sNc_1_n, &
                             corr_sNc_2_n, KK_auto_tndcy, KK_auto_coef, &
                             t_tol, crt1, crt2, mixt_frac )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use KK_upscaled_means, only:  &
        bivar_NL_mean_eq  ! Procedure

    implicit none

    ! Input Variables
    real, intent(in) :: &
      mu_t_1,        & ! Mean of t (1st PDF component)                       [-]
      mu_t_2,        & ! Mean of t (2nd PDF component)                       [-]
      mu_s_1,        & ! Mean of s (1st PDF component)                       [-]
      mu_s_2,        & ! Mean of s (2nd PDF component)                       [-]
      mu_Nc_n,       & ! Mean of ln Nc (both components)                     [-]
      sigma_t_1,     & ! Standard deviation of t (1st PDF component)         [-]
      sigma_t_2,     & ! Standard deviation of t (2nd PDF component)         [-]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)         [-]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)         [-]
      sigma_Nc_n,    & ! Standard deviation of ln Nc (both components)       [-]
      corr_ts_1,     & ! Correlation between t and s (1st PDF component)     [-]
      corr_ts_2,     & ! Correlation between t and s (2nd PDF component)     [-]
      corr_tNc_1_n,  & ! Correlation between t and ln Nc (1st PDF component) [-]
      corr_tNc_2_n,  & ! Correlation between t and ln Nc (2nd PDF component) [-]
      corr_sNc_1_n,  & ! Correlation between s and ln Nc (1st PDF component) [-]
      corr_sNc_2_n,  & ! Correlation between s and ln Nc (2nd PDF component) [-]
      KK_auto_tndcy, & ! KK autoconversion tendency                  [(kg/kg)/s]
      KK_auto_coef,  & ! KK autoconversion coefficient               [(kg/kg)/s]
      t_tol,         & ! Tolerance value of t                       [units vary]
      crt1,          & ! Coefficient c_rt (1st PDF component)                [-]
      crt2,          & ! Coefficient c_rt (2nd PDF component)                [-]
      mixt_frac        ! Mixture fraction                                    [-]

    ! Return Variable
    real :: &
      covar_rt_KK_auto  ! Covariance between r_t and KK autoconv. tendency   [-]

    ! Constant Parameters
    real, parameter :: &
      alpha_exp = 2.47,  & ! Exponent on s                                   [-]
      beta_exp  = -1.79    ! Exponent on r_r                                 [-]


    ! Calculate the covariance of r_t and KK autoconversion tendency.
    covar_rt_KK_auto  &
    = KK_auto_coef  &
      * ( mixt_frac * ( 1.0 / ( 2.0 * crt1 ) )  &
          * ( trivar_NNL_covar_eq( mu_t_1, mu_s_1, mu_Nc_n, &
                                   sigma_t_1, sigma_s_1, sigma_Nc_n, &
                                   corr_ts_1, corr_tNc_1_n, corr_sNc_1_n, &
                                   mu_t_1, KK_auto_tndcy, KK_auto_coef, &
                                   t_tol, alpha_exp, beta_exp )  &
            + bivar_NL_mean_eq( mu_s_1, mu_Nc_n, sigma_s_1, sigma_Nc_n, &
                                corr_sNc_1_n, alpha_exp + 1.0, beta_exp )  &
            - mu_s_1  &
              * bivar_NL_mean_eq( mu_s_1, mu_Nc_n, sigma_s_1, sigma_Nc_n, &
                                  corr_sNc_1_n, alpha_exp, beta_exp )  &
            )  &
        + ( 1.0 - mixt_frac ) * ( 1.0 / ( 2.0 * crt2 ) )  &
          * ( trivar_NNL_covar_eq( mu_t_2, mu_s_2, mu_Nc_n, &
                                   sigma_t_2, sigma_s_2, sigma_Nc_n, &
                                   corr_ts_2, corr_tNc_2_n, corr_sNc_2_n, &
                                   mu_t_2, KK_auto_tndcy, KK_auto_coef, &
                                   t_tol, alpha_exp, beta_exp )  &
            + bivar_NL_mean_eq( mu_s_2, mu_Nc_n, sigma_s_2, sigma_Nc_n, &
                                corr_sNc_2_n, alpha_exp + 1.0, beta_exp )  &
            - mu_s_2  &
              * bivar_NL_mean_eq( mu_s_2, mu_Nc_n, sigma_s_2, sigma_Nc_n, &
                                  corr_sNc_2_n, alpha_exp, beta_exp )  &
            )  &
        )


    return

  end function covar_rt_KK_auto

  !=============================================================================
  function covar_thl_KK_auto( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_Nc_n, &
                              sigma_t_1, sigma_t_2, sigma_s_1, sigma_s_2, &
                              sigma_Nc_n, corr_ts_1, corr_ts_2, &
                              corr_tNc_1_n, corr_tNc_2_n, corr_sNc_1_n, &
                              corr_sNc_2_n, KK_auto_tndcy, KK_auto_coef, &
                              t_tol, cthl1, cthl2, mixt_frac )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use KK_upscaled_means, only:  &
        bivar_NL_mean_eq  ! Procedure

    implicit none

    ! Input Variables
    real, intent(in) :: &
      mu_t_1,        & ! Mean of t (1st PDF component)                       [-]
      mu_t_2,        & ! Mean of t (2nd PDF component)                       [-]
      mu_s_1,        & ! Mean of s (1st PDF component)                       [-]
      mu_s_2,        & ! Mean of s (2nd PDF component)                       [-]
      mu_Nc_n,       & ! Mean of ln Nc (both components)                     [-]
      sigma_t_1,     & ! Standard deviation of t (1st PDF component)         [-]
      sigma_t_2,     & ! Standard deviation of t (2nd PDF component)         [-]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)         [-]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)         [-]
      sigma_Nc_n,    & ! Standard deviation of ln Nc (both components)       [-]
      corr_ts_1,     & ! Correlation between t and s (1st PDF component)     [-]
      corr_ts_2,     & ! Correlation between t and s (2nd PDF component)     [-]
      corr_tNc_1_n,  & ! Correlation between t and ln Nc (1st PDF component) [-]
      corr_tNc_2_n,  & ! Correlation between t and ln Nc (2nd PDF component) [-]
      corr_sNc_1_n,  & ! Correlation between s and ln Nc (1st PDF component) [-]
      corr_sNc_2_n,  & ! Correlation between s and ln Nc (2nd PDF component) [-]
      KK_auto_tndcy, & ! KK autoconversion tendency                  [(kg/kg)/s]
      KK_auto_coef,  & ! KK autoconversion coefficient               [(kg/kg)/s]
      t_tol,         & ! Tolerance value of t                       [units vary]
      cthl1,         & ! Coefficient c_thl (1st PDF component)               [-]
      cthl2,         & ! Coefficient c_thl (2nd PDF component)               [-]
      mixt_frac        ! Mixture fraction                                    [-]

    ! Return Variable
    real :: &
      covar_thl_KK_auto  ! Covariance between th_l and KK autoconv. tendency [-]

    ! Constant Parameters
    real, parameter :: &
      alpha_exp = 2.47,  & ! Exponent on s                                   [-]
      beta_exp  = -1.79    ! Exponent on r_r                                 [-]


    ! Calculate the covariance of th_l and KK autoconversion tendency.
    covar_thl_KK_auto  &
    = KK_auto_coef  &
      * ( mixt_frac * ( 1.0 / ( 2.0 * cthl1 ) )  &
          * ( trivar_NNL_covar_eq( mu_t_1, mu_s_1, mu_Nc_n, &
                                   sigma_t_1, sigma_s_1, sigma_Nc_n, &
                                   corr_ts_1, corr_tNc_1_n, corr_sNc_1_n, &
                                   mu_t_1, KK_auto_tndcy, KK_auto_coef, &
                                   t_tol, alpha_exp, beta_exp )  &
            - bivar_NL_mean_eq( mu_s_1, mu_Nc_n, sigma_s_1, sigma_Nc_n, &
                                corr_sNc_1_n, alpha_exp + 1.0, beta_exp )  &
            + mu_s_1  &
              * bivar_NL_mean_eq( mu_s_1, mu_Nc_n, sigma_s_1, sigma_Nc_n, &
                                  corr_sNc_1_n, alpha_exp, beta_exp )  &
            )  &
        + ( 1.0 - mixt_frac ) * ( 1.0 / ( 2.0 * cthl2 ) )  &
          * ( trivar_NNL_covar_eq( mu_t_2, mu_s_2, mu_Nc_n, &
                                   sigma_t_2, sigma_s_2, sigma_Nc_n, &
                                   corr_ts_2, corr_tNc_2_n, corr_sNc_2_n, &
                                   mu_t_2, KK_auto_tndcy, KK_auto_coef, &
                                   t_tol, alpha_exp, beta_exp )  &
            - bivar_NL_mean_eq( mu_s_2, mu_Nc_n, sigma_s_2, sigma_Nc_n, &
                                corr_sNc_2_n, alpha_exp + 1.0, beta_exp )  &
            + mu_s_2  &
              * bivar_NL_mean_eq( mu_s_2, mu_Nc_n, sigma_s_2, sigma_Nc_n, &
                                  corr_sNc_2_n, alpha_exp, beta_exp )  &    
            )  &
        )


    return

  end function covar_thl_KK_auto

  !=============================================================================
  function covar_x_KK_accr( mu_x_1, mu_x_2, mu_s_1, mu_s_2, mu_rr_n, &
                            sigma_x_1, sigma_x_2, sigma_s_1, sigma_s_2, &
                            sigma_rr_n, corr_xs_1, corr_xs_2, &
                            corr_xrr_1_n, corr_xrr_2_n, corr_srr_1_n, &
                            corr_srr_2_n, x_mean, KK_accr_tndcy, &
                            KK_accr_coef, x_tol, mixt_frac )

    ! Description:
    ! This function calculates the covariance between x and KK accretion
    ! tendency, which can be written as < x'((dr_r/dt)_KKaccr)' >, or more
    ! simply as < x'KK_accr' >.

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) :: &
      mu_x_1,        & ! Mean of x (1st PDF component)                       [-]
      mu_x_2,        & ! Mean of x (2nd PDF component)                       [-]
      mu_s_1,        & ! Mean of s (1st PDF component)                       [-]
      mu_s_2,        & ! Mean of s (2nd PDF component)                       [-]
      mu_rr_n,       & ! Mean of ln rr (both components)                     [-]
      sigma_x_1,     & ! Standard deviation of x (1st PDF component)         [-]
      sigma_x_2,     & ! Standard deviation of x (2nd PDF component)         [-]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)         [-]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)         [-]
      sigma_rr_n,    & ! Standard deviation of ln rr (both components)       [-]
      corr_xs_1,     & ! Correlation between x and s (1st PDF component)     [-]
      corr_xs_2,     & ! Correlation between x and s (2nd PDF component)     [-]
      corr_xrr_1_n,  & ! Correlation between x and ln rr (1st PDF component) [-]
      corr_xrr_2_n,  & ! Correlation between x and ln rr (2nd PDF component) [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF component) [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF component) [-]
      x_mean,        & ! Mean of x (overall)                        [units vary]
      KK_accr_tndcy, & ! KK accretion tendency                       [(kg/kg)/s]
      KK_accr_coef,  & ! KK accretion coefficient                    [(kg/kg)/s]
      x_tol,         & ! Tolerance value of x                       [units vary]
      mixt_frac        ! Mixture fraction                                    [-]

    ! Return Variable
    real :: &
      covar_x_KK_accr  ! Covariance between x and KK accretion tendency      [-]

    ! Constant Parameters
    real, parameter :: &
      alpha_exp = 1.15, & ! Exponent on s                                    [-]
      beta_exp  = 1.15    ! Exponent on r_r                                  [-]


    ! Calculate the covariance of x and KK accretion tendency.
    covar_x_KK_accr  &
    = KK_accr_coef &
      * ( mixt_frac &
          * trivar_NNL_covar_eq( mu_x_1, mu_s_1, mu_rr_n, &
                                 sigma_x_1, sigma_s_1, sigma_rr_n, &
                                 corr_xs_1, corr_xrr_1_n, corr_srr_1_n, &
                                 x_mean, KK_accr_tndcy, KK_accr_coef, &
                                 x_tol, alpha_exp, beta_exp ) &
        + ( 1.0 - mixt_frac ) &
          * trivar_NNL_covar_eq( mu_x_2, mu_s_2, mu_rr_n, &
                                 sigma_x_2, sigma_s_2, sigma_rr_n, &
                                 corr_xs_2, corr_xrr_2_n, corr_srr_2_n, &
                                 x_mean, KK_accr_tndcy, KK_accr_coef, &
                                 x_tol, alpha_exp, beta_exp ) &
        )


    return

  end function covar_x_KK_accr

  !=============================================================================
  function covar_rt_KK_accr( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_rr_n, &
                             sigma_t_1, sigma_t_2, sigma_s_1, sigma_s_2, &
                             sigma_rr_n, corr_ts_1, corr_ts_2, &
                             corr_trr_1_n, corr_trr_2_n, corr_srr_1_n, &
                             corr_srr_2_n, KK_accr_tndcy, KK_accr_coef, &
                             t_tol, crt1, crt2, mixt_frac )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use KK_upscaled_means, only:  &
        bivar_NL_mean_eq

    implicit none

    ! Input Variables
    real, intent(in) :: &
      mu_t_1,        & ! Mean of t (1st PDF component)                       [-]
      mu_t_2,        & ! Mean of t (2nd PDF component)                       [-]
      mu_s_1,        & ! Mean of s (1st PDF component)                       [-]
      mu_s_2,        & ! Mean of s (2nd PDF component)                       [-]
      mu_rr_n,       & ! Mean of ln rr (both components)                     [-]
      sigma_t_1,     & ! Standard deviation of t (1st PDF component)         [-]
      sigma_t_2,     & ! Standard deviation of t (2nd PDF component)         [-]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)         [-]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)         [-]
      sigma_rr_n,    & ! Standard deviation of ln rr (both components)       [-]
      corr_ts_1,     & ! Correlation between t and s (1st PDF component)     [-]
      corr_ts_2,     & ! Correlation between t and s (2nd PDF component)     [-]
      corr_trr_1_n,  & ! Correlation between t and ln rr (1st PDF component) [-]
      corr_trr_2_n,  & ! Correlation between t and ln rr (2nd PDF component) [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF component) [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF component) [-]
      KK_accr_tndcy, & ! KK accretion tendency                       [(kg/kg)/s]
      KK_accr_coef,  & ! KK accretion coefficient                    [(kg/kg)/s]
      t_tol,         & ! Tolerance value of t                       [units vary]
      crt1,          & ! Coefficient c_rt (1st PDF component)                [-]
      crt2,          & ! Coefficient c_rt (2nd PDF component)                [-]
      mixt_frac        ! Mixture fraction                                    [-]

    ! Return Variable
    real :: &
      covar_rt_KK_accr  ! Covariance between r_t and KK accretion tendency   [-]

    ! Constant Parameters
    real, parameter :: &
      alpha_exp = 1.15, & ! Exponent on s                                    [-]
      beta_exp  = 1.15    ! Exponent on r_r                                  [-]


    ! Calculate the covariance of r_t and KK accretion tendency.
    covar_rt_KK_accr  &
    = KK_accr_coef  &
      * ( mixt_frac * ( 1.0 / ( 2.0 * crt1 ) )  &
          * ( trivar_NNL_covar_eq( mu_t_1, mu_s_1, mu_rr_n, &
                                   sigma_t_1, sigma_s_1, sigma_rr_n, &
                                   corr_ts_1, corr_trr_1_n, corr_srr_1_n, &
                                   mu_t_1, KK_accr_tndcy, KK_accr_coef, &
                                   t_tol, alpha_exp, beta_exp )  &
            + bivar_NL_mean_eq( mu_s_1, mu_rr_n, sigma_s_1, sigma_rr_n, &
                                corr_srr_1_n, alpha_exp + 1.0, beta_exp )  &
            - mu_s_1  &
              * bivar_NL_mean_eq( mu_s_1, mu_rr_n, sigma_s_1, sigma_rr_n, &
                                  corr_srr_1_n, alpha_exp, beta_exp )  &
            )  &
        + ( 1.0 - mixt_frac ) * ( 1.0 / ( 2.0 * crt2 ) )  &
          * ( trivar_NNL_covar_eq( mu_t_2, mu_s_2, mu_rr_n, &
                                   sigma_t_2, sigma_s_2, sigma_rr_n, &
                                   corr_ts_2, corr_trr_2_n, corr_srr_2_n, &
                                   mu_t_2, KK_accr_tndcy, KK_accr_coef, &
                                   t_tol, alpha_exp, beta_exp )  &
            + bivar_NL_mean_eq( mu_s_2, mu_rr_n, sigma_s_2, sigma_rr_n, &
                                corr_srr_2_n, alpha_exp + 1.0, beta_exp )  &
            - mu_s_2  &
              * bivar_NL_mean_eq( mu_s_2, mu_rr_n, sigma_s_2, sigma_rr_n, &
                                  corr_srr_2_n, alpha_exp, beta_exp )  &
            )  &
        )


    return

  end function covar_rt_KK_accr

  !=============================================================================
  function covar_thl_KK_accr( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_rr_n, &
                              sigma_t_1, sigma_t_2, sigma_s_1, sigma_s_2, &
                              sigma_rr_n, corr_ts_1, corr_ts_2, &
                              corr_trr_1_n, corr_trr_2_n, corr_srr_1_n, &
                              corr_srr_2_n, KK_accr_tndcy, KK_accr_coef, &
                              t_tol, cthl1, cthl2, mixt_frac )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use KK_upscaled_means, only:  &
        bivar_NL_mean_eq

    implicit none

    ! Input Variables
    real, intent(in) :: &
      mu_t_1,        & ! Mean of t (1st PDF component)                       [-]
      mu_t_2,        & ! Mean of t (2nd PDF component)                       [-]
      mu_s_1,        & ! Mean of s (1st PDF component)                       [-]
      mu_s_2,        & ! Mean of s (2nd PDF component)                       [-]
      mu_rr_n,       & ! Mean of ln rr (both components)                     [-]
      sigma_t_1,     & ! Standard deviation of t (1st PDF component)         [-]
      sigma_t_2,     & ! Standard deviation of t (2nd PDF component)         [-]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)         [-]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)         [-]
      sigma_rr_n,    & ! Standard deviation of ln rr (both components)       [-]
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
      mixt_frac        ! Mixture fraction                                    [-]

    ! Return Variable
    real :: &
      covar_thl_KK_accr  ! Covariance between th_l and KK accretion tendency [-]

    ! Constant Parameters
    real, parameter :: &
      alpha_exp = 1.15, & ! Exponent on s                                    [-]
      beta_exp  = 1.15    ! Exponent on r_r                                  [-]


    ! Calculate the covariance of th_l and KK accretion tendency.
    covar_thl_KK_accr  &
    = KK_accr_coef  &
      * ( mixt_frac * ( 1.0 / ( 2.0 * cthl1 ) )  &
          * ( trivar_NNL_covar_eq( mu_t_1, mu_s_1, mu_rr_n, &
                                   sigma_t_1, sigma_s_1, sigma_rr_n, &
                                   corr_ts_1, corr_trr_1_n, corr_srr_1_n, &
                                   mu_t_1, KK_accr_tndcy, KK_accr_coef, &
                                   t_tol, alpha_exp, beta_exp )  &
            - bivar_NL_mean_eq( mu_s_1, mu_rr_n, sigma_s_1, sigma_rr_n, &
                                corr_srr_1_n, alpha_exp + 1.0, beta_exp )  &
            + mu_s_1  &
              * bivar_NL_mean_eq( mu_s_1, mu_rr_n, sigma_s_1, sigma_rr_n, &
                                  corr_srr_1_n, alpha_exp, beta_exp )  &
            )  &
        + ( 1.0 - mixt_frac ) * ( 1.0 / ( 2.0 * cthl2 ) )  &
          * ( trivar_NNL_covar_eq( mu_t_2, mu_s_2, mu_rr_n, &
                                   sigma_t_2, sigma_s_2, sigma_rr_n, &
                                   corr_ts_2, corr_trr_2_n, corr_srr_2_n, &
                                   mu_t_2, KK_accr_tndcy, KK_accr_coef, &
                                   t_tol, alpha_exp, beta_exp )  &
            - bivar_NL_mean_eq( mu_s_2, mu_rr_n, sigma_s_2, sigma_rr_n, &
                                corr_srr_2_n, alpha_exp + 1.0, beta_exp )  &
            + mu_s_2  &
              * bivar_NL_mean_eq( mu_s_2, mu_rr_n, sigma_s_2, sigma_rr_n, &
                                  corr_srr_2_n, alpha_exp, beta_exp )  &
            )  &
        )


    return

  end function covar_thl_KK_accr

  !=============================================================================
  function quadrivar_NNLL_covar_eq( mu_x_i, mu_s_i, mu_rr_n, mu_Nr_n, &
                                    sigma_x_i, sigma_s_i, sigma_rr_n, &
                                    sigma_Nr_n, corr_xs_i, corr_xrr_i_n, &
                                    corr_xNr_i_n, corr_srr_i_n, &
                                    corr_sNr_i_n, corr_rrNr_n, x_mean, &
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
        quadrivar_NNLL_covar,              & ! Procedure(s)
        quadrivar_NNLL_covar_const_x1,     &
        quadrivar_NNLL_covar_const_x2,     &
        quadrivar_NNLL_covar_const_x1_x2

    use constants_clubb, only: &
        s_mellor_tol, & ! Constant(s)
        parab_cyl_max_input

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real, intent(in) :: &
      mu_x_i,       & ! Mean of x (ith PDF component)                       [-]
      mu_s_i,       & ! Mean of s (ith PDF component)                       [-]
      mu_rr_n,      & ! Mean of ln rr (both components)                     [-]
      mu_Nr_n,      & ! Mean of ln Nr (both components)                     [-]
      sigma_x_i,    & ! Standard deviation of x (ith PDF component)         [-]
      sigma_s_i,    & ! Standard deviation of s (ith PDF component)         [-]
      sigma_rr_n,   & ! Standard deviation of ln rr (both components)       [-]
      sigma_Nr_n,   & ! Standard deviation of ln Nr (both components)       [-]
      corr_xs_i,    & ! Correlation between x and s (ith PDF component)     [-]
      corr_xrr_i_n, & ! Correlation between x and ln rr (ith PDF component) [-]
      corr_xNr_i_n, & ! Correlation between x and ln Nr (ith PDF component) [-]
      corr_srr_i_n, & ! Correlation between s and ln rr (ith PDF component) [-]
      corr_sNr_i_n, & ! Correlation between s and ln Nr (ith PDF component) [-]
      corr_rrNr_n     ! Correlation between ln rr & ln Nr (both components) [-]

    real, intent(in) :: &
      x_mean,        & ! Mean of x (overall)                       [units vary]
      mc_tndcy_mean, & ! Mean of microphysics tendency              [(kg/kg)/s]
      mc_coef,       & ! Coefficient of microphysics tendency       [(kg/kg)/s]
      x_tol            ! Tolerance value of x                      [units vary]

    real, intent(in) :: &
      alpha_exp_in,  & ! Exponent alpha, corresponding to s                 [-]
      beta_exp_in,   & ! Exponent beta, corresponding to rr                 [-]
      gamma_exp_in     ! Exponent gamma, corresponding to Nr                [-]

    ! Return Variable
    real :: &
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
    mu_x3_n = dble( mu_rr_n ) ! The same for both PDF components.
    mu_x4_n = dble( mu_Nr_n ) ! The same for both PDF components.

    ! Standard deviations for the ith PDF component.
    sigma_x1   = dble( sigma_x_i )  ! x is w or t (ith component).
    sigma_x2   = dble( sigma_s_i )
    sigma_x3_n = dble( sigma_rr_n ) ! The same for both PDF components.
    sigma_x4_n = dble( sigma_Nr_n ) ! The same for both PDF components.

    ! Correlations for the ith PDF component.
    rho_x1x2   = dble( corr_xs_i )    ! x is w or t (ith component).
    rho_x1x3_n = dble( corr_xrr_i_n ) ! x is w or t (ith component).
    rho_x1x4_n = dble( corr_xNr_i_n ) ! x is w or t (ith component).
    rho_x2x3_n = dble( corr_srr_i_n )
    rho_x2x4_n = dble( corr_sNr_i_n )
    rho_x3x4_n = dble( corr_rrNr_n )  ! The same for both PDF components.

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
         ( sigma_x2 <= x2_tol .or. abs( s_cc ) > parab_cyl_max_input ) ) then

       ! The ith PDF component variance of both x (r_t, th_l, or w) and s is 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_const_x1_x2( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                           sigma_x3_n, sigma_x4_n, &
                                           rho_x3x4_n, x1_mean, &
                                           x2_alpha_x3_beta_x4_gamma_mean, &
                                           alpha_exp, beta_exp, gamma_exp ) )


    elseif ( sigma_x1 <= x1_tol ) then

       ! The ith PDF component variance of x (r_t, th_l, or w) is 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_const_x1( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                        sigma_x2, sigma_x3_n, sigma_x4_n, &
                                        rho_x2x3_n, rho_x2x4_n, rho_x3x4_n, &
                                        x1_mean, &
                                        x2_alpha_x3_beta_x4_gamma_mean, &
                                        alpha_exp, beta_exp, gamma_exp ) )


    elseif ( sigma_x2 <= x2_tol .or. abs( s_cc ) > parab_cyl_max_input ) then

       ! The ith PDF component variance of s is 0.
       quadrivar_NNLL_covar_eq  &
       = real( &
         quadrivar_NNLL_covar_const_x2( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                        sigma_x1, sigma_x3_n, sigma_x4_n, &
                                        rho_x1x3_n, rho_x1x4_n, rho_x3x4_n, &
                                        x1_mean, &
                                        x2_alpha_x3_beta_x4_gamma_mean, &
                                        alpha_exp, beta_exp, gamma_exp ) )


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
                               alpha_exp, beta_exp, gamma_exp ) )


    endif


    return

  end function quadrivar_NNLL_covar_eq

  !=============================================================================
  function trivar_NNL_covar_eq( mu_x_i, mu_s_i, mu_y_n, &
                                sigma_x_i, sigma_s_i, sigma_y_n, &
                                corr_xs_i, corr_xy_i_n, corr_sy_i_n, &
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
        trivar_NNL_covar,            & ! Procedure(s)
        trivar_NNL_covar_const_x1,   &
        trivar_NNL_covar_const_x2,   &
        trivar_NNL_covar_const_x1_x2

    use constants_clubb, only: &
        s_mellor_tol, & ! Constant(s)
        parab_cyl_max_input

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real, intent(in) :: &
      mu_x_i,      & ! Mean of x (ith PDF component)                        [-]
      mu_s_i,      & ! Mean of s (ith PDF component)                        [-]
      mu_y_n,      & ! Mean of ln y (both components)                       [-]
      sigma_x_i,   & ! Standard deviation of x (ith PDF component)          [-]
      sigma_s_i,   & ! Standard deviation of s (ith PDF component)          [-]
      sigma_y_n,   & ! Standard deviation of ln y (both components)         [-]
      corr_xs_i,   & ! Correlation between x and s (ith PDF component)      [-]
      corr_xy_i_n, & ! Correlation between x and ln y (ith PDF component)   [-]
      corr_sy_i_n    ! Correlation between s and ln y (ith PDF component)   [-]

    real, intent(in) :: &
      x_mean,        & ! Mean of x (overall)                       [units vary]
      mc_tndcy_mean, & ! Mean of microphysics tendency              [(kg/kg)/s]
      mc_coef,       & ! Coefficient of microphysics                [(kg/kg)/s]
      x_tol            ! Tolerance value of x                      [units vary]

    real, intent(in) :: &
      alpha_exp_in,  & ! Exponent alpha, corresponding to s                 [-]
      beta_exp_in      ! Exponent beta, corresponding to y                  [-]

    ! Return Variable
    real :: &
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
    mu_x3_n = dble( mu_y_n ) ! y is N_c (autoconversion) or r_r (accretion).
                             ! The same for both PDF components.

    ! Standard deviations for the ith PDF component.
    sigma_x1   = dble( sigma_x_i ) ! x is w or t (ith component).
    sigma_x2   = dble( sigma_s_i )
    sigma_x3_n = dble( sigma_y_n ) ! y is N_c (auto.) or r_r (accr.).
                                   ! The same for both PDF components.

    ! Correlations for the ith PDF component.
    rho_x1x2   = dble( corr_xs_i )    ! x is w or t (ith component).
    rho_x1x3_n = dble( corr_xy_i_n )  ! x is w or t (ith component).
                                      ! y is N_c (auto.) or r_r (accr.).
    rho_x2x3_n = dble( corr_sy_i_n )  ! y is N_c (auto.) or r_r (accr.).

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
         ( sigma_x2 <= x2_tol .or. abs( s_c ) > parab_cyl_max_input ) ) then

       ! The ith PDF component variance of both x (r_t, th_l, or w) and s is 0.
       trivar_NNL_covar_eq  &
       = real( &
         trivar_NNL_covar_const_x1_x2( mu_x1, mu_x2, mu_x3_n, sigma_x3_n, &
                                       x1_mean, x2_alpha_x3_beta_mean, &
                                       alpha_exp, beta_exp ) )


    elseif ( sigma_x1 <= x1_tol ) then

       ! The ith PDF component variance of x (r_t, th_l, or w) is 0.
       trivar_NNL_covar_eq  &
       = real( trivar_NNL_covar_const_x1( mu_x1, mu_x2, mu_x3_n, &
                                          sigma_x2, sigma_x3_n, rho_x2x3_n, &
                                          x1_mean, x2_alpha_x3_beta_mean, &
                                          alpha_exp, beta_exp ) )


    elseif ( sigma_x2 <= x2_tol .or. abs( s_c ) > parab_cyl_max_input ) then

       ! The ith PDF component variance of s is 0.
       trivar_NNL_covar_eq  &
       = real( trivar_NNL_covar_const_x2( mu_x1, mu_x2, mu_x3_n, &
                                          sigma_x1, sigma_x3_n, rho_x1x3_n, &
                                          x1_mean, x2_alpha_x3_beta_mean, &
                                          alpha_exp, beta_exp ) )


    else  ! sigma_x1 > 0 and sigma_x2 > 0.

       ! This is the complete value of the trivariate.
       ! All fields vary in the ith PDF component.
       trivar_NNL_covar_eq  &
       = real( trivar_NNL_covar( mu_x1, mu_x2, mu_x3_n, &
                                 sigma_x1, sigma_x2, sigma_x3_n, &
                                 rho_x1x2, rho_x1x3_n, rho_x2x3_n, &
                                 x1_mean, x2_alpha_x3_beta_mean, &
                                 alpha_exp, beta_exp ) )


    endif


    return

  end function trivar_NNL_covar_eq

!===============================================================================

end module KK_upscaled_covariances
