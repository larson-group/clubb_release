! $Id$
!===============================================================================
module KK_upscaled_means

  implicit none

  private

  public :: KK_upscaled_means_driver, &
            KK_evap_upscaled_mean,    &
            KK_auto_upscaled_mean,    &
            KK_accr_upscaled_mean,    &
            KK_mvr_upscaled_mean,     &
            trivar_NLL_mean_eq,       &
            bivar_NL_mean_eq,         &
            bivar_NL_mean_eq_Nc0,     &
            bivar_LL_mean_eq

  contains

  !=============================================================================
  subroutine KK_upscaled_means_driver( mu_s_1, mu_s_2, mu_rr_1, mu_rr_2, &
                                       mu_Nr_1, mu_Nr_2, mu_Nc_1, mu_Nc_2, &
                                       mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                       mu_Nr_2_n, mu_Nc_1_n, mu_Nc_2_n, &
                                       sigma_s_1, sigma_s_2, &
                                       sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                       sigma_Nr_2, sigma_Nc_1, sigma_Nc_2, &
                                       sigma_rr_1_n, sigma_rr_2_n, &
                                       sigma_Nr_1_n, sigma_Nr_2_n, &
                                       sigma_Nc_1_n, sigma_Nc_2_n, &
                                       corr_srr_1_n, corr_srr_2_n, &
                                       corr_sNr_1_n, corr_sNr_2_n, &
                                       corr_sNc_1_n, corr_sNc_2_n, &
                                       corr_rrNr_1_n, corr_rrNr_2_n, &
                                       mixt_frac, precip_frac_1, &
                                       precip_frac_2, Nc0_in_cloud, &
                                       l_const_Nc_in_cloud, &
                                       KK_evap_coef, KK_auto_coef, &
                                       KK_accr_coef, KK_mvr_coef, &
                                       KK_evap_tndcy, KK_auto_tndcy, &
                                       KK_accr_tndcy, KK_mean_vol_rad )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_s_1,        & ! Mean of s (1st PDF component)                   [kg/kg]
      mu_s_2,        & ! Mean of s (2nd PDF component)                   [kg/kg]
      mu_rr_1,       & ! Mean of rr (1st PDF component) in-precip (ip)   [kg/kg]
      mu_rr_2,       & ! Mean of rr (2nd PDF component) ip               [kg/kg]
      mu_Nr_1,       & ! Mean of Nr (1st PDF component) ip              [num/kg]
      mu_Nr_2,       & ! Mean of Nr (2nd PDF component) ip              [num/kg]
      mu_Nc_1,       & ! Mean of Nc (1st PDF component)                 [num/kg]
      mu_Nc_2,       & ! Mean of Nc (2nd PDF component)                 [num/kg]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF comp.) ip            [ln(kg/kg)]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF comp.) ip            [ln(kg/kg)]
      mu_Nr_1_n,     & ! Mean of ln Nr (1st PDF component) ip       [ln(num/kg)]
      mu_Nr_2_n,     & ! Mean of ln Nr (2nd PDF component) ip       [ln(num/kg)]
      mu_Nc_1_n,     & ! Mean of ln Nc (1st PDF component)          [ln(num/kg)]
      mu_Nc_2_n,     & ! Mean of ln Nc (2nd PDF component)          [ln(num/kg)]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)     [kg/kg]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)     [kg/kg]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_Nr_1,    & ! Standard deviation of Nr (1st PDF comp.) ip    [num/kg]
      sigma_Nr_2,    & ! Standard deviation of Nr (2nd PDF comp.) ip    [num/kg]
      sigma_Nc_1,    & ! Standard deviation of Nc (1st PDF component)   [num/kg]
      sigma_Nc_2,    & ! Standard deviation of Nc (2nd PDF component)   [num/kg]
      sigma_rr_1_n,  & ! Standard dev. of ln rr (1st PDF comp.) ip   [ln(kg/kg)]
      sigma_rr_2_n,  & ! Standard dev. of ln rr (2nd PDF comp.) ip   [ln(kg/kg)]
      sigma_Nr_1_n,  & ! Standard dev. of ln Nr (1st PDF comp.) ip  [ln(num/kg)]
      sigma_Nr_2_n,  & ! Standard dev. of ln Nr (2nd PDF comp.) ip  [ln(num/kg)]
      sigma_Nc_1_n,  & ! Standard dev. of ln Nc (1st PDF comp.)     [ln(num/kg)]
      sigma_Nc_2_n,  & ! Standard dev. of ln Nc (2nd PDF comp.)     [ln(num/kg)]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF comp.) ip  [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF comp.) ip  [-]
      corr_sNc_1_n,  & ! Correlation between s and ln Nc (1st PDF comp.) ip  [-]
      corr_sNc_2_n,  & ! Correlation between s and ln Nc (2nd PDF comp.) ip  [-]
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
      KK_accr_coef, & ! KK accretion coefficient            [(kg/kg)/s]
      KK_mvr_coef     ! KK mean volume radius coefficient   [m]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      KK_evap_tndcy,   & ! KK evaporation tendency          [(kg/kg)/s]
      KK_auto_tndcy,   & ! KK autoconversion tendency       [(kg/kg)/s]
      KK_accr_tndcy,   & ! KK accretion tendency            [(kg/kg)/s]
      KK_mean_vol_rad    ! KK rain drop mean volume radius  [m]


    !!! Calculate the upscaled KK evaporation tendency.
    KK_evap_tndcy  &
    = KK_evap_upscaled_mean( mu_s_1, mu_s_2, mu_rr_1, mu_rr_2, mu_Nr_1, &
                             mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                             mu_Nr_2_n, sigma_s_1, sigma_s_2, sigma_rr_1, &
                             sigma_rr_2, sigma_Nr_1, sigma_Nr_2, &
                             sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                             sigma_Nr_2_n, corr_srr_1_n, corr_srr_2_n, &
                             corr_sNr_1_n, corr_sNr_2_n, corr_rrNr_1_n, &
                             corr_rrNr_2_n, KK_evap_coef, mixt_frac, &
                             precip_frac_1, precip_frac_2 )


    !!! Calculate the upscaled KK autoconversion tendency.
    KK_auto_tndcy  &
    = KK_auto_upscaled_mean( mu_s_1, mu_s_2, mu_Nc_1, mu_Nc_2, &
                             mu_Nc_1_n, mu_Nc_2_n, sigma_s_1, &
                             sigma_s_2, sigma_Nc_1, sigma_Nc_2, &
                             sigma_Nc_1_n, sigma_Nc_2_n, corr_sNc_1_n, &
                             corr_sNc_2_n, KK_auto_coef, mixt_frac, &
                             Nc0_in_cloud, l_const_Nc_in_cloud )


    !!! Calculate the upscaled KK accretion tendency.
    KK_accr_tndcy  &
    = KK_accr_upscaled_mean( mu_s_1, mu_s_2, mu_rr_1, mu_rr_2, &
                             mu_rr_1_n, mu_rr_2_n, sigma_s_1, &
                             sigma_s_2, sigma_rr_1, sigma_rr_2, &
                             sigma_rr_1_n, sigma_rr_2_n, corr_srr_1_n, &
                             corr_srr_2_n, KK_accr_coef, mixt_frac, &
                             precip_frac_1, precip_frac_2 )


    !!! Calculate the upscaled KK rain drop mean volume radius.
    KK_mean_vol_rad &
    = KK_mvr_upscaled_mean( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
                            mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                            sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                            sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, &
                            sigma_Nr_1_n, sigma_Nr_2_n, corr_rrNr_1_n, &
                            corr_rrNr_2_n, KK_mvr_coef, mixt_frac, &
                            precip_frac_1, precip_frac_2 )


    return

  end subroutine KK_upscaled_means_driver

  !=============================================================================
  function KK_evap_upscaled_mean( mu_s_1, mu_s_2, mu_rr_1, mu_rr_2, mu_Nr_1, &
                                  mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                  mu_Nr_2_n, sigma_s_1, sigma_s_2, sigma_rr_1, &
                                  sigma_rr_2, sigma_Nr_1, sigma_Nr_2, &
                                  sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                                  sigma_Nr_2_n, corr_srr_1_n, corr_srr_2_n, &
                                  corr_sNr_1_n, corr_sNr_2_n, corr_rrNr_1_n, &
                                  corr_rrNr_2_n, KK_evap_coef, mixt_frac, &
                                  precip_frac_1, precip_frac_2 )

    ! Description:
    ! This function calculates the mean value of the upscaled KK rain water
    ! evaporation tendency.

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
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF comp.) ip  [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF comp.) ip  [-]
      corr_rrNr_1_n, & ! Correlation btwn. ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rrNr_2_n, & ! Correlation btwn. ln rr & ln Nr (2nd PDF comp.) ip  [-]
      KK_evap_coef,  & ! KK evaporation coefficient                  [(kg/kg)/s]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)          [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      KK_evap_upscaled_mean  ! Mean of KK evaporation tendency       [(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on s                                           [-]
      beta_exp,  & ! Exponent on r_r                                         [-]
      gamma_exp    ! Exponent on N_r                                         [-]


    ! Values of the KK exponents.
    alpha_exp = KK_evap_Supersat_exp
    beta_exp  = KK_evap_rr_exp
    gamma_exp = KK_evap_Nr_exp

    ! Calculate the mean KK evaporation tendency.
    KK_evap_upscaled_mean  &
    = KK_evap_coef &
      * ( mixt_frac &
          * precip_frac_1 &
          * trivar_NLL_mean_eq( mu_s_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, &
                                mu_Nr_1_n, sigma_s_1, sigma_rr_1, &
                                sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, &
                                corr_srr_1_n, corr_sNr_1_n, corr_rrNr_1_n, &
                                alpha_exp, beta_exp, gamma_exp ) &
        + ( one - mixt_frac ) &
          * precip_frac_2 &
          * trivar_NLL_mean_eq( mu_s_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, &
                                mu_Nr_2_n, sigma_s_2, sigma_rr_2, &
                                sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, &
                                corr_srr_2_n, corr_sNr_2_n, corr_rrNr_2_n, &
                                alpha_exp, beta_exp, gamma_exp ) &
        )


    return

  end function KK_evap_upscaled_mean

  !=============================================================================
  function KK_auto_upscaled_mean( mu_s_1, mu_s_2, mu_Nc_1, mu_Nc_2, &
                                  mu_Nc_1_n, mu_Nc_2_n, sigma_s_1, &
                                  sigma_s_2, sigma_Nc_1, sigma_Nc_2, &
                                  sigma_Nc_1_n, sigma_Nc_2_n, corr_sNc_1_n, &
                                  corr_sNc_2_n, KK_auto_coef, mixt_frac, &
                                  Nc0_in_cloud, l_const_Nc_in_cloud )

    ! Description:
    ! This function calculates the mean value of the upscaled KK rain water
    ! autoconversion tendency.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,    & ! Constant(s)
        Nc_tol

    use parameters_microphys, only: &
        KK_auto_rc_exp, & ! Variable(s)
        KK_auto_Nc_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_s_1,       & ! Mean of s (1st PDF component)                   [kg/kg]
      mu_s_2,       & ! Mean of s (2nd PDF component)                   [kg/kg]
      mu_Nc_1,      & ! Mean of Nc (1st PDF component)                 [num/kg]
      mu_Nc_2,      & ! Mean of Nc (2nd PDF component)                 [num/kg]
      mu_Nc_1_n,    & ! Mean of ln Nc (1st PDF component)                   [-]
      mu_Nc_2_n,    & ! Mean of ln Nc (2nd PDF component)                   [-]
      sigma_s_1,    & ! Standard deviation of s (1st PDF component)     [kg/kg]
      sigma_s_2,    & ! Standard deviation of s (2nd PDF component)     [kg/kg]
      sigma_Nc_1,   & ! Standard deviation of Nc (1st PDF component)   [num/kg]
      sigma_Nc_2,   & ! Standard deviation of Nc (2nd PDF component)   [num/kg]
      sigma_Nc_1_n, & ! Standard deviation of ln Nc (1st PDF component)     [-]
      sigma_Nc_2_n, & ! Standard deviation of ln Nc (2nd PDF component)     [-]
      corr_sNc_1_n, & ! Correlation between s and ln Nc (1st PDF component) [-]
      corr_sNc_2_n, & ! Correlation between s and ln Nc (2nd PDF component) [-]
      KK_auto_coef, & ! KK autoconversion coefficient               [(kg/kg)/s]
      mixt_frac,    & ! Mixture fraction                                    [-]
      Nc0_in_cloud    ! Constant in-cloud value of cloud droplet conc. [num/kg]

    logical, intent(in) :: &
      l_const_Nc_in_cloud  ! Flag to use a constant value of N_c within cloud

    ! Return Variable
    real( kind = core_rknd ) :: &
      KK_auto_upscaled_mean  ! Mean of KK autoconversion tendency   [(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on s                                          [-]
      beta_exp     ! Exponent on N_c                                        [-]


    ! Values of the KK exponents.
    alpha_exp = KK_auto_rc_exp
    beta_exp  = KK_auto_Nc_exp

    ! Calculate the mean KK autoconversion tendency.
    if ( l_const_Nc_in_cloud ) then

       KK_auto_upscaled_mean  &
       = KK_auto_coef &
         * ( mixt_frac &
             * bivar_NL_mean_eq_Nc0( mu_s_1, Nc0_in_cloud, sigma_s_1, &
                                     alpha_exp, beta_exp ) &
           + ( one - mixt_frac ) &
             * bivar_NL_mean_eq_Nc0( mu_s_2, Nc0_in_cloud, sigma_s_2, &
                                     alpha_exp, beta_exp ) &
           )

    else

       KK_auto_upscaled_mean  &
       = KK_auto_coef &
         * ( mixt_frac &
             * bivar_NL_mean_eq( mu_s_1, mu_Nc_1, mu_Nc_1_n, sigma_s_1, &
                                 sigma_Nc_1, sigma_Nc_1_n, corr_sNc_1_n, &
                                 Nc_tol, alpha_exp, beta_exp ) &
           + ( one - mixt_frac ) &
             * bivar_NL_mean_eq( mu_s_2, mu_Nc_2, mu_Nc_2_n, sigma_s_2, &
                                 sigma_Nc_2, sigma_Nc_2_n, corr_sNc_2_n, &
                                 Nc_tol, alpha_exp, beta_exp ) &
           )

    endif


    return

  end function KK_auto_upscaled_mean

  !=============================================================================
  function KK_accr_upscaled_mean( mu_s_1, mu_s_2, mu_rr_1, mu_rr_2, &
                                  mu_rr_1_n, mu_rr_2_n, sigma_s_1, &
                                  sigma_s_2, sigma_rr_1, sigma_rr_2, &
                                  sigma_rr_1_n, sigma_rr_2_n, corr_srr_1_n, &
                                  corr_srr_2_n, KK_accr_coef, mixt_frac, &
                                  precip_frac_1, precip_frac_2 )

    ! Description:
    ! This function calculates the mean value of the upscaled KK rain water
    ! accretion tendency.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,    & ! Constant(s)
        rr_tol

    use parameters_microphys, only: &
        KK_accr_rc_exp, & ! Variable(s)
        KK_accr_rr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_s_1,        & ! Mean of s (1st PDF component)                   [kg/kg]
      mu_s_2,        & ! Mean of s (2nd PDF component)                   [kg/kg]
      mu_rr_1,       & ! Mean of rr (1st PDF component) in-precip (ip)   [kg/kg]
      mu_rr_2,       & ! Mean of rr (2nd PDF component) ip               [kg/kg]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) ip                [-]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip                [-]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)     [kg/kg]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)     [kg/kg]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_rr_1_n,  & ! Standard deviation of ln rr (1st PDF component) ip  [-]
      sigma_rr_2_n,  & ! Standard deviation of ln rr (2nd PDF component) ip  [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      KK_accr_coef,  & ! KK accretion coefficient                    [(kg/kg)/s]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)          [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      KK_accr_upscaled_mean  ! Mean of KK accretion tendency         [(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on s                                           [-]
      beta_exp     ! Exponent on r_r                                         [-]


    ! Values of the KK exponents.
    alpha_exp = KK_accr_rc_exp
    beta_exp  = KK_accr_rr_exp

    ! Calculate the mean KK accretion tendency.
    KK_accr_upscaled_mean  &
    = KK_accr_coef &
      * ( mixt_frac &
          * precip_frac_1 &
          * bivar_NL_mean_eq( mu_s_1, mu_rr_1, mu_rr_1_n, sigma_s_1, &
                              sigma_rr_1, sigma_rr_1_n, corr_srr_1_n, &
                              rr_tol, alpha_exp, beta_exp ) &
        + ( one - mixt_frac ) &
          * precip_frac_2 &
          * bivar_NL_mean_eq( mu_s_2, mu_rr_2, mu_rr_2_n, sigma_s_2, &
                              sigma_rr_2, sigma_rr_2_n, corr_srr_2_n, &
                              rr_tol, alpha_exp, beta_exp ) &
        )


    return

  end function KK_accr_upscaled_mean

  !=============================================================================
  function KK_mvr_upscaled_mean( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
                                 mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                                 sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                 sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, &
                                 sigma_Nr_1_n, sigma_Nr_2_n, corr_rrNr_1_n, &
                                 corr_rrNr_2_n, KK_mvr_coef, mixt_frac, &
                                 precip_frac_1, precip_frac_2 )

    ! Description:
    ! This function calculates the mean value of the upscaled KK rain drop mean
    ! volume radius.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one  ! Constant(s)

    use parameters_microphys, only: &
        KK_mvr_rr_exp, & ! Variable(s)
        KK_mvr_Nr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_rr_1,       & ! Mean of rr (1st PDF component) in-precip (ip)   [kg/kg]
      mu_rr_2,       & ! Mean of rr (2nd PDF component) ip               [kg/kg]
      mu_Nr_1,       & ! Mean of Nr (1st PDF component) ip              [num/kg]
      mu_Nr_2,       & ! Mean of Nr (2nd PDF component) ip              [num/kg]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) ip                [-]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip                [-]
      mu_Nr_1_n,     & ! Mean of ln Nr (1st PDF component) ip                [-]
      mu_Nr_2_n,     & ! Mean of ln Nr (2nd PDF component) ip                [-]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_Nr_1,    & ! Standard deviation of Nr (1st PDF component) ip  [#/kg]
      sigma_Nr_2,    & ! Standard deviation of Nr (2nd PDF component) ip  [#/kg]
      sigma_rr_1_n,  & ! Standard deviation of ln rr (1st PDF component) ip  [-]
      sigma_rr_2_n,  & ! Standard deviation of ln rr (2nd PDF component) ip  [-]
      sigma_Nr_1_n,  & ! Standard deviation of ln Nr (1st PDF component) ip  [-]
      sigma_Nr_2_n,  & ! Standard deviation of ln Nr (2nd PDF component) ip  [-]
      corr_rrNr_1_n, & ! Correlation btwn. ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rrNr_2_n, & ! Correlation btwn. ln rr & ln Nr (2nd PDF comp.) ip  [-]
      KK_mvr_coef,   & ! KK mean volume radius coefficient                   [m]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)          [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      KK_mvr_upscaled_mean  ! Mean of KK rain drop mean volume radius        [m]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on r_r                                         [-]
      beta_exp     ! Exponent on N_r                                         [-]


    ! Values of the KK exponents.
    alpha_exp = KK_mvr_rr_exp
    beta_exp  = KK_mvr_Nr_exp

    ! Calculate the KK mean volume radius of rain drops.
    KK_mvr_upscaled_mean  &
    = KK_mvr_coef &
      * ( mixt_frac &
          * precip_frac_1 &
          * bivar_LL_mean_eq( mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, &
                              sigma_rr_1, sigma_Nr_1, sigma_rr_1_n, &
                              sigma_Nr_1_n, corr_rrNr_1_n, &
                              alpha_exp, beta_exp ) &
        + ( one - mixt_frac ) &
          * precip_frac_2 &
          * bivar_LL_mean_eq( mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, &
                              sigma_rr_2, sigma_Nr_2, sigma_rr_2_n, &
                              sigma_Nr_2_n, corr_rrNr_2_n, &
                              alpha_exp, beta_exp ) &
        ) 


    return

  end function KK_mvr_upscaled_mean

  !=============================================================================
  function trivar_NLL_mean_eq( mu_s_i, mu_rr_i, mu_Nr_i, mu_rr_i_n, &
                               mu_Nr_i_n, sigma_s_i, sigma_rr_i, &
                               sigma_Nr_i, sigma_rr_i_n, sigma_Nr_i_n, &
                               corr_srr_i_n, corr_sNr_i_n, corr_rrNr_i_n, &
                               alpha_exp_in, beta_exp_in, gamma_exp_in )

    ! Description:
    ! This function calculates the contribution by the ith PDF component to the
    ! expression < x1^alpha x2^beta x3^gamma >, where x1 = s, x2 = r_r, and
    ! x3 = N_r.  The total value of KK mean evaporation tendency is given by:
    !
    ! < KK_evap >
    !    = KK_evap_coef
    !      * ( mixt_frac < s^alpha r_r^beta N_r^gamma (1) >
    !          + ( 1 - mixt_frac ) < s^alpha r_r^beta N_r^gamma (2) > ).
    !
    ! One of two functions is called, based on whether x1 (s) varies.  Each one
    ! of these two functions is the result of an evaluated integral based on the
    ! specific situation.

    ! References:
    !-----------------------------------------------------------------------

    use PDF_integrals_means, only: &
        trivar_NLL_mean,            & ! Procedure(s)
        trivar_NLL_mean_const_x1,   &
        trivar_NLL_mean_const_x2,   &
        trivar_NLL_mean_const_x1x2, &
        trivar_NLL_mean_const_x2x3, &
        trivar_NLL_mean_const_all

    use constants_clubb, only: &
        s_mellor_tol,        & ! Constant(s)
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
      mu_s_i,        & ! Mean of s (ith PDF component)                   [kg/kg]
      mu_rr_i,       & ! Mean of rr (ith PDF component) in-precip (ip)   [kg/kg]
      mu_Nr_i,       & ! Mean of Nr (ith PDF component) ip              [num/kg]
      mu_rr_i_n,     & ! Mean of ln rr (ith PDF component) in-precip (ip)    [-]
      mu_Nr_i_n,     & ! Mean of ln Nr (ith PDF component) ip                [-]
      sigma_s_i,     & ! Standard deviation of s (ith PDF component)     [kg/kg]
      sigma_rr_i,    & ! Standard deviation of rr (ith PDF component) ip [kg/kg]
      sigma_Nr_i,    & ! Standard deviation of Nr (ith PDF component) ip  [#/kg]
      sigma_rr_i_n,  & ! Standard deviation of ln rr (ith PDF component) ip  [-]
      sigma_Nr_i_n,  & ! Standard deviation of ln Nr (ith PDF component) ip  [-]
      corr_srr_i_n,  & ! Correlation between s and ln rr (ith PDF comp.) ip  [-]
      corr_sNr_i_n,  & ! Correlation between s and ln Nr (ith PDF comp.) ip  [-]
      corr_rrNr_i_n    ! Correlation btwn. ln rr & ln Nr (ith PDF comp.) ip  [-]

    real( kind = core_rknd ), intent(in) :: &
      alpha_exp_in,  & ! Exponent alpha, corresponding to s                  [-]
      beta_exp_in,   & ! Exponent beta, corresponding to rr                  [-]
      gamma_exp_in     ! Exponent gamma, corresponding to Nr                 [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      trivar_NLL_mean_eq

    ! Local Variables
    real( kind = dp ) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3,      & ! Mean of x3 (ith PDF component)                        [-]
      mu_x2_n,    & ! Mean of ln x2 (ith PDF component)                     [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x2,   & ! Standard deviation of x2 (ith PDF component)          [-]
      sigma_x3,   & ! Standard deviation of x3 (ith PDF component)          [-]
      sigma_x2_n, & ! Standard deviation of ln x2 (ith PDF component)       [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      rho_x1x2_n, & ! Correlation between x1 and ln x2 (ith PDF component)  [-]
      rho_x1x3_n, & ! Correlation between x1 and ln x3 (ith PDF component)  [-]
      rho_x2x3_n    ! Correlation between ln x2 & ln x3 (ith PDF component) [-]

    real( kind = dp ) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x2                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x3                   [-]

    real( kind = dp ) :: &
      x1_tol, & ! Tolerance value of x1                                     [-]
      x2_tol, & ! Tolerance value of x2                                     [-]
      x3_tol, & ! Tolerance value of x3                                     [-]
      s_cc      ! Parabolic cylinder function input value                   [-]


    ! Means for the ith PDF component. 
    mu_x1 = dble( mu_s_i )
    if ( beta_exp_in >= zero ) then
       mu_x2 = dble( mu_rr_i )
    else ! exponent beta < 0
       mu_x2 = dble( max( mu_rr_i, rr_tol ) )
    endif
    if ( gamma_exp_in >= zero ) then
       mu_x3 = dble( mu_Nr_i )
    else ! exponent gamma < 0
       mu_x3 = dble( max( mu_Nr_i, Nr_tol ) )
    endif
    mu_x2_n = dble( mu_rr_i_n )
    mu_x3_n = dble( mu_Nr_i_n )

    ! Standard deviations for the ith PDF component.
    sigma_x1   = dble( sigma_s_i )
    sigma_x2   = dble( sigma_rr_i )
    sigma_x3   = dble( sigma_Nr_i )
    sigma_x2_n = dble( sigma_rr_i_n )
    sigma_x3_n = dble( sigma_Nr_i_n )

    ! Correlations for the ith PDF component.
    rho_x1x2_n = dble( corr_srr_i_n )
    rho_x1x3_n = dble( corr_sNr_i_n )
    rho_x2x3_n = dble( corr_rrNr_i_n )

    ! Exponents.
    alpha_exp = dble( alpha_exp_in )
    beta_exp  = dble( beta_exp_in )
    gamma_exp = dble( gamma_exp_in )

    ! Tolerance values.
    ! When the standard deviation of a variable is below the tolerance values,
    ! it is considered to be zero, and the variable is considered to have a
    ! constant value.
    x1_tol = dble( s_mellor_tol )
    x2_tol = dble( rr_tol )
    x3_tol = dble( Nr_tol )

    ! Determine the value of the parabolic cylinder function input value, s_cc.
    ! The value s_cc is being fed into the parabolic cylinder function.  When
    ! the value of s_cc is too large in magnitude (depending on the order of the
    ! parabolic cylinder function), overflow occurs, and the output of the
    ! parabolic cylinder function is +/-Inf.  This is primarily due to a large
    ! ratio of mu_x1 to sigma_x1.  When the value of s_cc is very large, the
    ! distribution of x1 is basically a spike near the mean, so x1 is treated as
    ! a constant.
    if ( sigma_x1 > x1_tol ) then
       s_cc = ( mu_x1 / sigma_x1 )  &
              + rho_x1x2_n * sigma_x2_n * beta_exp  &
              + rho_x1x3_n * sigma_x3_n * gamma_exp
    else  ! sigma_x1 = 0
       ! Note:  s_cc is +inf when mu_x1 > 0 and sigma_x1 = 0, and s_cc is -inf
       !        when mu_x1 < 0 and sigma_x1 = 0.  Furthermore, s_cc is undefined
       !        when mu_x1 = 0 and sigma_x1 = 0.  However, within the context of
       !        this particular function, only the absolute value of s_cc is
       !        relevant, and furthermore the absolute value of s_cc is only
       !        relevant when sigma_x1 > 0.  Therefore, this statement only
       !        serves as divide-by-zero and compiler warning prevention.
       s_cc = huge( s_cc )
    endif


    ! Based on the value of sigma_x1 (including the value of s_cc compared to
    ! parab_cyl_max_input), find the correct form of the trivariate equation to
    ! use.

    if ( ( sigma_x1 <= x1_tol .or. abs( s_cc ) > dble( parab_cyl_max_input ) ) &
         .and. sigma_x2 <= x2_tol .and. sigma_x3 <= x3_tol ) then

       ! The ith PDF component variance of each of s, r_r, and N_r is 0.
       trivar_NLL_mean_eq  &
       = real( trivar_NLL_mean_const_all( mu_x1, mu_x2, mu_x3, &
                                          alpha_exp, beta_exp, gamma_exp ), &
               kind = core_rknd )


    elseif ( ( sigma_x1 <= x1_tol &
               .or. abs( s_cc ) > dble( parab_cyl_max_input ) ) &
             .and. sigma_x2 <= x2_tol ) then

       ! The ith PDF component variance of both s and r_r is 0.
       trivar_NLL_mean_eq  &
       = real( trivar_NLL_mean_const_x1x2( mu_x1, mu_x2, mu_x3_n, sigma_x3_n, &
                                           alpha_exp, beta_exp, gamma_exp ), &
               kind = core_rknd )


    elseif ( ( sigma_x1 <= x1_tol &
               .or. abs( s_cc ) > dble( parab_cyl_max_input ) ) &
             .and. sigma_x3 <= x3_tol ) then

       ! The ith PDF component variance of both s and N_r is 0.
       trivar_NLL_mean_eq  &
       = real( trivar_NLL_mean_const_x1x2( mu_x1, mu_x3, mu_x2_n, sigma_x2_n, &
                                           alpha_exp, gamma_exp, beta_exp ), &
               kind = core_rknd )


    elseif ( sigma_x2 <= x2_tol .and. sigma_x3 <= x3_tol ) then

       ! The ith PDF component variance of both r_r and N_r is 0.
       trivar_NLL_mean_eq  &
       = real( trivar_NLL_mean_const_x2x3( mu_x1, mu_x2, mu_x3, sigma_x1, &
                                           alpha_exp, beta_exp, gamma_exp ), &
               kind = core_rknd )


    elseif ( sigma_x1 <= x1_tol &
             .or. abs( s_cc ) > dble( parab_cyl_max_input ) ) then

       ! The ith PDF component variance of s is 0.
       trivar_NLL_mean_eq  &
       = real( trivar_NLL_mean_const_x1( mu_x1, mu_x2_n, mu_x3_n, &
                                         sigma_x2_n, sigma_x3_n, rho_x2x3_n, &
                                         alpha_exp, beta_exp, gamma_exp ),  &
               kind = core_rknd )


    elseif ( sigma_x2 <= x2_tol ) then

       ! The ith PDF component variance of r_r is 0.
       trivar_NLL_mean_eq  &
       = real( trivar_NLL_mean_const_x2( mu_x1, mu_x2, mu_x3_n, &
                                         sigma_x1, sigma_x3_n, rho_x1x3_n, &
                                         alpha_exp, beta_exp, gamma_exp ), &
               kind = core_rknd )


    elseif ( sigma_x3 <= x3_tol ) then

       ! The ith PDF component variance of N_r is 0.
       trivar_NLL_mean_eq  &
       = real( trivar_NLL_mean_const_x2( mu_x1, mu_x3, mu_x2_n, &
                                         sigma_x1, sigma_x2_n, rho_x1x2_n, &
                                         alpha_exp, gamma_exp, beta_exp ), &
               kind = core_rknd )


    else  ! sigma_x1, sigma_x2, and sigma_x3 > 0

       ! All fields vary in the ith PDF component.
       trivar_NLL_mean_eq  &
       = real( trivar_NLL_mean( mu_x1, mu_x2_n, mu_x3_n, &
                                sigma_x1, sigma_x2_n, sigma_x3_n, &
                                rho_x1x2_n, rho_x1x3_n, rho_x2x3_n, &
                                alpha_exp, beta_exp, gamma_exp ),  &
               kind = core_rknd )


    endif


    return

  end function trivar_NLL_mean_eq

  !=============================================================================
  function bivar_NL_mean_eq( mu_s_i, mu_y_i, mu_y_i_n, sigma_s_i, &
                             sigma_y_i, sigma_y_i_n, corr_sy_i_n, &
                             y_tol, alpha_exp_in, beta_exp_in )

    ! Description:
    ! This function calculates the contribution by the ith PDF component to the
    ! expression < x1^alpha x2^beta >, where x1 = s and x2 = N_ccn or r_r,
    ! depending on whether this function is being called for autoconversion or
    ! accretion, respectively.  The total value of KK mean microphysics tendency
    ! is given by:
    !
    ! < KK_mc > = KK_mc_coef
    !             * ( mixt_frac < s^alpha y^beta (1) >
    !                 + ( 1 - mixt_frac ) < s^alpha y^beta (2) > );
    !
    ! where y stands for either N_ccn or r_r.  One of two functions is called,
    ! based on whether x1 (s) and/or x2 (N_ccn or r_r) varies.  Each one of
    ! these two functions is the result of an evaluated integral based on the
    ! specific situation.

    ! References:
    !-----------------------------------------------------------------------

    use PDF_integrals_means, only: &
        bivar_NL_mean,           & ! Procedure(s)
        bivar_NL_mean_const_x1,  &
        bivar_NL_mean_const_x2,  &
        bivar_NL_mean_const_all

    use constants_clubb, only: &
        s_mellor_tol,        & ! Constant(s)
        parab_cyl_max_input, &
        zero

    use clubb_precision, only: &
        dp,        & ! double precision
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_s_i,      & ! Mean of s (ith PDF component)                    [kg/kg]
      mu_y_i,      & ! Mean of y (ith PDF component)                        [-]
      mu_y_i_n,    & ! Mean of ln y (ith PDF component)                     [-]
      sigma_s_i,   & ! Standard deviation of s (ith PDF component)      [kg/kg]
      sigma_y_i,   & ! Standard deviation of y (ith PDF component)          [-]
      sigma_y_i_n, & ! Standard deviation of ln y (ith PDF component)       [-]
      corr_sy_i_n    ! Correlation between s and ln y (ith PDF component)   [-]

    real( kind = core_rknd ), intent(in) :: &
      y_tol          ! Tolerance value of y                        [units vary]

    real( kind = core_rknd ), intent(in) :: &
      alpha_exp_in,  & ! Exponent alpha, corresponding to s                 [-]
      beta_exp_in      ! Exponent beta, corresponding to y                  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      bivar_NL_mean_eq

    ! Local Variables
    real( kind = dp ) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x2_n,    & ! Mean of ln x2 (ith PDF component)                     [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x2,   & ! Standard deviation of x2 (ith PDF component)          [-]
      sigma_x2_n, & ! Standard deviation of ln x2 (ith PDF component)       [-]
      rho_x1x2_n    ! Correlation between x1 and ln x2 (ith PDF component)  [-]
    
    real( kind = dp ) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp      ! Exponent beta, corresponding to x2                    [-]

    real( kind = dp ) :: &
      x1_tol, & ! Tolerance value of x1                                     [-]
      x2_tol, & ! Tolerance value of x2                                     [-]
      s_c       ! Parabolic cylinder function input value                   [-]


    ! Means for the ith PDF component. 
    mu_x1 = dble( mu_s_i )
    if ( beta_exp_in >= zero ) then
       mu_x2 = dble( mu_y_i )   ! y is N_ccn (auto.) or r_r (accr.).
    else ! exponent beta < 0
       mu_x2 = dble( max( mu_y_i, y_tol ) )   ! y is N_ccn or r_r.
    endif
    mu_x2_n = dble( mu_y_i_n ) ! y is N_ccn (autoconversion) or r_r (accretion).

    ! Standard deviations for the ith PDF component.
    sigma_x1   = dble( sigma_s_i )
    sigma_x2   = dble( sigma_y_i )   ! y is N_ccn (auto.) or r_r (accr.).
    sigma_x2_n = dble( sigma_y_i_n ) ! y is N_ccn (auto.) or r_r (accr.).

    ! Correlations for the ith PDF component.
    rho_x1x2_n = dble( corr_sy_i_n ) ! y is N_ccn (auto.) or r_r (accr.).

    ! Exponents.
    alpha_exp = dble( alpha_exp_in )
    beta_exp  = dble( beta_exp_in )

    ! Tolerance values.
    ! When the standard deviation of a variable is below the tolerance values,
    ! it is considered to be zero, and the variable is considered to have a
    ! constant value.
    x1_tol = dble( s_mellor_tol )
    x2_tol = dble( y_tol )  ! y is N_ccn (auto.) or r_r (accr.).

    ! Determine the value of the parabolic cylinder function input value, s_c.
    ! The value s_c is being fed into the parabolic cylinder function.  When
    ! the value of s_c is too large in magnitude (depending on the order of the
    ! parabolic cylinder function), overflow occurs, and the output of the
    ! parabolic cylinder function is +/-Inf.  This is primarily due to a large
    ! ratio of mu_x1 to sigma_x1.  When the value of s_c is very large, the
    ! distribution of x1 is basically a spike near the mean, so x1 is treated as
    ! a constant.
    if ( sigma_x1 > x1_tol ) then
       s_c = ( mu_x1 / sigma_x1 ) + rho_x1x2_n * sigma_x2_n * beta_exp
    else  ! sigma_x1 = 0
       ! Note:  s_c is +inf when mu_x1 > 0 and sigma_x1 = 0, and s_c is -inf
       !        when mu_x1 < 0 and sigma_x1 = 0.  Furthermore, s_c is undefined
       !        when mu_x1 = 0 and sigma_x1 = 0.  However, within the context of
       !        this particular function, only the absolute value of s_c is
       !        relevant, and furthermore the absolute value of s_c is only
       !        relevant when sigma_x1 > 0.  Therefore, this statement only
       !        serves as divide-by-zero and compiler warning prevention.
       s_c = huge( s_c )
    endif


    ! Based on the values of sigma_x1 (including the value of s_c compared to
    ! parab_cyl_max_input) and sigma_x2, find the correct form of the bivariate
    ! equation to use.

    if ( ( sigma_x1 <= x1_tol .or. abs( s_c ) > dble( parab_cyl_max_input ) ) &
         .and. sigma_x2 <= x2_tol ) then

       ! The ith PDF component variance of both s and y (r_r or N_ccn) is 0.
       bivar_NL_mean_eq  &
       = real( bivar_NL_mean_const_all( mu_x1, mu_x2, alpha_exp, beta_exp ), &
               kind = core_rknd )


    elseif ( sigma_x1 <= x1_tol &
             .or. abs( s_c ) > dble( parab_cyl_max_input ) ) then

       ! The ith PDF component variance of s is 0.
       bivar_NL_mean_eq  &
       = real( bivar_NL_mean_const_x1( mu_x1, mu_x2_n, sigma_x2_n, &
                                       alpha_exp, beta_exp ), kind = core_rknd )


    elseif ( sigma_x2 <= x2_tol ) then

       ! The ith PDF component variance of y (r_r or N_ccn) is 0.
       bivar_NL_mean_eq  &
       = real( bivar_NL_mean_const_x2( mu_x1, mu_x2, sigma_x1, &
                                       alpha_exp, beta_exp ), kind = core_rknd )


    else  ! sigma_x1 and sigma_x2 > 0

       ! All fields vary in the ith PDF component.
       bivar_NL_mean_eq  &
       = real( bivar_NL_mean( mu_x1, mu_x2_n, sigma_x1, sigma_x2_n, &
                              rho_x1x2_n, alpha_exp, beta_exp ),  &
               kind = core_rknd )


    endif


    return

  end function bivar_NL_mean_eq

  !=============================================================================
  function bivar_NL_mean_eq_Nc0( mu_s_i, Nc0_in_cloud, sigma_s_i, &
                                 alpha_exp_in, beta_exp_in )

    ! Description:
    ! This function calculates the contribution by the ith PDF component to the
    ! expression < x1^alpha x2^beta >, where x1 = s and x2 = N_c or r_r,
    ! depending on whether this function is being called for autoconversion or
    ! accretion, respectively.  The total value of KK mean microphysics tendency
    ! is given by:
    !
    ! < KK_mc > = KK_mc_coef
    !             * ( mixt_frac < s^alpha y^beta (1) >
    !                 + ( 1 - mixt_frac ) < s^alpha y^beta (2) > );
    !
    ! where y stands for either N_c or r_r.  One of two functions is called,
    ! based on whether x1 (s) varies.  Each one of these two functions is the
    ! result of an evaluated integral based on the specific situation.

    ! References:
    !-----------------------------------------------------------------------

    use PDF_integrals_means, only: &
        bivar_NL_mean_const_all, &
        bivar_NL_mean_const_x2

    use constants_clubb, only: &
        s_mellor_tol,        & ! Constant(s)
        Nc_tol,              &
        parab_cyl_max_input, &
        zero

    use clubb_precision, only: &
        dp,        & ! double precision
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_s_i,       & ! Mean of s (ith PDF component)                   [kg/kg]
      Nc0_in_cloud, & ! Constant in-cloud value of cloud droplet conc. [num/kg]
      sigma_s_i       ! Standard deviation of s (ith PDF component)     [kg/kg]

    real( kind = core_rknd ), intent(in) :: &
      alpha_exp_in,  & ! Exponent alpha, corresponding to s                 [-]
      beta_exp_in      ! Exponent beta, corresponding to y                  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      bivar_NL_mean_eq_Nc0

    ! Local Variables
    real( kind = dp ) :: &
      mu_x1,    & ! Mean of x1 (ith PDF component)                          [-]
      Nc0,      & ! Constant in-cloud value of cloud droplet conc.     [num/kg]
      sigma_x1    ! Standard deviation of x1 (ith PDF component)            [-]
    
    real( kind = dp ) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp      ! Exponent beta, corresponding to x2                    [-]

    real( kind = dp ) :: &
      x1_tol, & ! Tolerance value of x1                                     [-]
      s_c       ! Parabolic cylinder function input value                   [-]


    ! Means for the ith PDF component. 
    mu_x1 = dble( mu_s_i )
    if ( beta_exp_in >= zero ) then
       Nc0 = dble( Nc0_in_cloud )
    else ! exponent beta < 0
       Nc0 = dble( max( Nc0_in_cloud, Nc_tol ) )
    endif

    ! Standard deviations for the ith PDF component.
    sigma_x1   = dble( sigma_s_i )

    ! Exponents.
    alpha_exp = dble( alpha_exp_in )
    beta_exp  = dble( beta_exp_in )

    ! Tolerance values.
    ! When the standard deviation of a variable is below the tolerance values,
    ! it is considered to be zero, and the variable is considered to have a
    ! constant value.
    x1_tol = dble( s_mellor_tol )

    ! Determine the value of the parabolic cylinder function input value, s_c.
    ! The value s_c is being fed into the parabolic cylinder function.  When
    ! the value of s_c is too large in magnitude (depending on the order of the
    ! parabolic cylinder function), overflow occurs, and the output of the
    ! parabolic cylinder function is +/-Inf.  This is primarily due to a large
    ! ratio of mu_x1 to sigma_x1.  When the value of s_c is very large, the
    ! distribution of x1 is basically a spike near the mean, so x1 is treated as
    ! a constant.
    if ( sigma_x1 > x1_tol ) then
       s_c = mu_x1 / sigma_x1
    else  ! sigma_x1 = 0
       ! Note:  s_c is +inf when mu_x1 > 0 and sigma_x1 = 0, and s_c is -inf
       !        when mu_x1 < 0 and sigma_x1 = 0.  Furthermore, s_c is undefined
       !        when mu_x1 = 0 and sigma_x1 = 0.  However, within the context of
       !        this particular function, only the absolute value of s_c is
       !        relevant, and furthermore the absolute value of s_c is only
       !        relevant when sigma_x1 > 0.  Therefore, this statement only
       !        serves as divide-by-zero and compiler warning prevention.
       s_c = huge( s_c )
    endif


    ! Based on the value of sigma_x1 (including the value of s_c compared to
    ! parab_cyl_max_input), find the correct form of the bivariate equation to
    ! use.

    if ( sigma_x1 <= x1_tol .or. abs( s_c ) > dble( parab_cyl_max_input ) ) then

       ! The ith PDF component variance of s is 0.
       bivar_NL_mean_eq_Nc0  &
       = real( bivar_NL_mean_const_all( mu_x1, Nc0, alpha_exp, beta_exp ), &
               kind = core_rknd )


    else  ! sigma_x1 > 0

       ! All fields vary in the ith PDF component.
       bivar_NL_mean_eq_Nc0  &
       = real( bivar_NL_mean_const_x2( mu_x1, Nc0, sigma_x1, &
                                       alpha_exp, beta_exp ), kind = core_rknd )


    endif


    return

  end function bivar_NL_mean_eq_Nc0

  !=============================================================================
  function bivar_LL_mean_eq( mu_rr_i, mu_Nr_i, mu_rr_i_n, mu_Nr_i_n, &
                             sigma_rr_i, sigma_Nr_i, sigma_rr_i_n, &
                             sigma_Nr_i_n, corr_rrNr_i_n, &
                             alpha_exp_in, beta_exp_in )

    ! Description:
    ! This function calculates the expression < x1^alpha x2^beta >, where
    ! x1 = r_r and x2 = N_r.  The value of mean volume radius is given by:
    !
    ! < KK_mvr > = KK_mvr_coef
    !              * ( mixt_frac < r_r^alpha N_r^beta (1) >
    !                  + ( 1 - mixt_frac ) < r_r^alpha N_r^beta (2) > ).
    !
    ! These function is the result of an evaluated integral.

    ! References:
    !-----------------------------------------------------------------------

    use PDF_integrals_means, only: &
        bivar_LL_mean,           & ! Procedure(s)
        bivar_LL_mean_const_x1,  &
        bivar_LL_mean_const_all

    use constants_clubb, only: &
        rr_tol, & ! Constant(s)
        Nr_tol, &
        zero

    use clubb_precision, only: &
        dp,        & ! double precision
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_rr_i,       & ! Mean of rr (ith PDF component) in-precip (ip)   [kg/kg]
      mu_Nr_i,       & ! Mean of Nr (ith PDF component) ip              [num/kg]
      mu_rr_i_n,     & ! Mean of ln rr (ith PDF component) ip                [-]
      mu_Nr_i_n,     & ! Mean of ln Nr (ith PDF component) ip                [-]
      sigma_rr_i,    & ! Standard deviation of rr (ith PDF component) ip [kg/kg]
      sigma_Nr_i,    & ! Standard deviation of Nr (ith PDF component) ip  [#/kg]
      sigma_rr_i_n,  & ! Standard deviation of ln rr (ith PDF component) ip  [-]
      sigma_Nr_i_n,  & ! Standard deviation of ln Nr (ith PDF component) ip  [-]
      corr_rrNr_i_n    ! Correlation btwn. ln rr & ln Nr (ith PDF comp.) ip  [-]

    real( kind = core_rknd ), intent(in) :: &
      alpha_exp_in,  & ! Exponent alpha, corresponding to rr                 [-]
      beta_exp_in      ! Exponent beta, corresponding to Nr                  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      bivar_LL_mean_eq

    ! Local Variables
    real( kind = dp ) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x1_n,    & ! Mean of ln x1 (ith PDF component)                     [-]
      mu_x2_n,    & ! Mean of ln x2 (ith PDF component)                     [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x2,   & ! Standard deviation of x2 (ith PDF component)          [-]
      sigma_x1_n, & ! Standard deviation of ln x1 (ith PDF component)       [-]
      sigma_x2_n, & ! Standard deviation of ln x2 (ith PDF component)       [-]
      rho_x1x2_n    ! Correlation between ln x1 & ln x2 (ith PDF component) [-]

    real( kind = dp ) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp      ! Exponent beta, corresponding to x2                    [-]

    real( kind = dp ) :: &
      x1_tol, & ! Tolerance value of x1                                     [-]
      x2_tol    ! Tolerance value of x2                                     [-]


    ! Means for the ith PDF component.
    if ( alpha_exp_in >= zero ) then
       mu_x1 = dble( mu_rr_i )
    else ! exponent alpha < 0
       mu_x1 = dble( max( mu_rr_i, rr_tol ) )
    endif
    if ( beta_exp_in >= zero ) then
       mu_x2 = dble( mu_Nr_i )
    else ! exponent beta < 0
       mu_x2 = dble( max( mu_Nr_i, Nr_tol ) )
    endif
    mu_x1_n = dble( mu_rr_i_n )
    mu_x2_n = dble( mu_Nr_i_n )

    ! Standard deviations for the ith PDF component.
    sigma_x1   = dble( sigma_rr_i )
    sigma_x2   = dble( sigma_Nr_i )
    sigma_x1_n = dble( sigma_rr_i_n )
    sigma_x2_n = dble( sigma_Nr_i_n )

    ! Correlations for the ith PDF component.
    rho_x1x2_n = dble( corr_rrNr_i_n )

    ! Exponents.
    alpha_exp = dble( alpha_exp_in )
    beta_exp  = dble( beta_exp_in )

    ! Tolerance values.
    ! When the standard deviation of a variable is below the tolerance values,
    ! it is considered to be zero, and the variable is considered to have a
    ! constant value.
    x1_tol = dble( rr_tol )
    x2_tol = dble( Nr_tol )


    ! Calculate the mean of the bivariate lognormal equation.
    if ( sigma_x1 <= x1_tol .and. sigma_x2 <= x2_tol ) then

       ! The ith PDF component variance of both r_r and N_r is 0.
       bivar_LL_mean_eq  &
       = real( bivar_LL_mean_const_all( mu_x1, mu_x2, alpha_exp, beta_exp ), &
               kind = core_rknd )


    elseif ( sigma_x1 <= x1_tol ) then

       ! The ith PDF component variance of r_r is 0.
       bivar_LL_mean_eq  &
       = real( bivar_LL_mean_const_x1( mu_x1, mu_x2_n, sigma_x2_n, &
                                       alpha_exp, beta_exp ), &
               kind = core_rknd )


    elseif ( sigma_x2 <= x2_tol ) then

       ! The ith PDF component variance of N_r is 0.
       bivar_LL_mean_eq  &
       = real( bivar_LL_mean_const_x1( mu_x2, mu_x1_n, sigma_x1_n, &
                                       beta_exp, alpha_exp ), &
               kind = core_rknd )


    else  ! sigma_x1 and sigma_x2 > 0

       ! All fields vary in the ith PDF component.
       bivar_LL_mean_eq  &
       = real( bivar_LL_mean( mu_x1_n, mu_x2_n, sigma_x1_n, sigma_x2_n, &
                              rho_x1x2_n, alpha_exp, beta_exp ),  &
               kind = core_rknd )


    endif


    return

  end function bivar_LL_mean_eq

!===============================================================================

end module KK_upscaled_means
