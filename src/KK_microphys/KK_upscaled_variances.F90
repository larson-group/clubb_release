! $Id$
!===============================================================================
module KK_upscaled_variances

  implicit none

  private

  public :: variance_KK_mvr

  contains

  !=============================================================================
  function variance_KK_mvr( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, &
                            mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1, &
                            sigma_rr_2, sigma_Nr_1, sigma_Nr_2, &
                            sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                            sigma_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, &
                            KK_mean_vol_rad, KK_mvr_coef, mixt_frac, &
                            precip_frac_1, precip_frac_2 )

    ! Description:
    ! This function calculates the variance of KK mean volume radius of rain
    ! drops (R_vr), which can be written as < R_vr'^2 >.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one, & ! Constant(s)
        two

    use KK_upscaled_means, only:  &
        bivar_LL_mean_eq  ! Procedure(s)

    use parameters_KK, only: &
        KK_mvr_rr_exp, & ! Variable(s)
        KK_mvr_Nr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_rr_1,         & ! Mean of rr (1st PDF component) in-precip (ip) [kg/kg]
      mu_rr_2,         & ! Mean of rr (2nd PDF component) ip             [kg/kg]
      mu_Nr_1,         & ! Mean of Nr (1st PDF component) ip            [num/kg]
      mu_Nr_2,         & ! Mean of Nr (2nd PDF component) ip            [num/kg]
      mu_rr_1_n,       & ! Mean of ln rr (1st PDF component) ip              [-]
      mu_rr_2_n,       & ! Mean of ln rr (2nd PDF component) ip              [-]
      mu_Nr_1_n,       & ! Mean of ln Nr (1st PDF component) ip              [-]
      mu_Nr_2_n,       & ! Mean of ln Nr (2nd PDF component) ip              [-]
      sigma_rr_1,      & ! Standard deviation of rr (1st PDF comp.) ip   [kg/kg]
      sigma_rr_2,      & ! Standard deviation of rr (2nd PDF comp.) ip   [kg/kg]
      sigma_Nr_1,      & ! Standard deviation of Nr (1st PDF comp.) ip  [num/kg]
      sigma_Nr_2,      & ! Standard deviation of Nr (2nd PDF comp.) ip  [num/kg]
      sigma_rr_1_n,    & ! Standard deviation of ln rr (1st PDF comp.) ip    [-]
      sigma_rr_2_n,    & ! Standard deviation of ln rr (2nd PDF comp.) ip    [-]
      sigma_Nr_1_n,    & ! Standard deviation of ln Nr (1st PDF comp.) ip    [-]
      sigma_Nr_2_n,    & ! Standard deviation of ln Nr (2nd PDF comp.) ip    [-]
      corr_rr_Nr_1_n,  & ! Correlation of ln rr & ln Nr (1st PDF comp.) ip   [-]
      corr_rr_Nr_2_n,  & ! Correlation of ln rr & ln Nr (2nd PDF comp.) ip   [-]
      KK_mean_vol_rad, & ! KK mean volume radius of rain drops               [m]
      KK_mvr_coef,     & ! KK mean volume radius coefficient                 [m]
      mixt_frac,       & ! Mixture fraction                                  [-]
      precip_frac_1,   & ! Precipitation fraction (1st PDF component)        [-]
      precip_frac_2      ! Precipitation fraction (2nd PDF component)        [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      variance_KK_mvr    ! Variance of KK rain drop mean volume radius     [m^2]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on r_r
      beta_exp     ! Exponent on N_r


    ! Values of the KK exponents.
    alpha_exp = KK_mvr_rr_exp
    beta_exp  = KK_mvr_Nr_exp

    ! Calculate the covariance of r_r and KK mean volume radius of rain drops.
    variance_KK_mvr  &
    = KK_mvr_coef**2 &
      * ( mixt_frac &
          * precip_frac_1 &
          * bivar_LL_mean_eq( mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, &
                              sigma_rr_1, sigma_Nr_1, sigma_rr_1_n, &
                              sigma_Nr_1_n, corr_rr_Nr_1_n, &
                              two * alpha_exp, two * beta_exp ) &
        + ( one - mixt_frac ) &
          * precip_frac_2 &
          * bivar_LL_mean_eq( mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, &
                              sigma_rr_2, sigma_Nr_2, sigma_rr_2_n, &
                              sigma_Nr_2_n, corr_rr_Nr_2_n, &
                              two * alpha_exp, two * beta_exp ) &
        ) &
      - KK_mean_vol_rad**2


    return

    end function variance_KK_mvr

!===============================================================================

end module KK_upscaled_variances
