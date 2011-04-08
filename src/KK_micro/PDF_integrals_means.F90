! $Id$
!===============================================================================
module PDF_integrals_means

  private

  public :: trivar_NLL_mean, &
            trivar_NLL_mean_const_x1, &
            bivar_NL_mean, &
            bivar_NL_mean_const_x1, &
            bivar_LL_mean

  contains

  !=============================================================================
  function trivar_NLL_mean( mu_x1, mu_x2_n, mu_x3_n, &
                            sigma_x1, sigma_x2_n, sigma_x3_n, &
                            rho_x1x2_n, rho_x1x3_n, rho_x2x3_n, &
                            alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    !  Larson, V. E. and B. M. Griffin (2010)
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp  ! Constant(s)

    use KK_utilities, only:  &
        Dv_fnc  ! Procedure(s)

    implicit none

    ! Input Variables
    double precision, intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2_n,    & ! Mean of ln x2 (ith PDF component)                     [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x2_n, & ! Standard deviation of ln x2 (ith PDF component)       [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      rho_x1x2_n, & ! Correlation between x1 and ln x2 (ith PDF component)  [-]
      rho_x1x3_n, & ! Correlation between x1 and ln x3 (ith PDF component)  [-]
      rho_x2x3_n, & ! Correlation between ln x2 & ln x3 (ith PDF component) [-]
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x2                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x3                   [-]

    ! Return Variable
    double precision ::  &
      trivar_NLL_mean

    ! Local Variable
    double precision ::  &
      s_cc    !

    s_cc = ( mu_x1 / sigma_x1 )  &
           + rho_x1x2_n * sigma_x2_n * beta_exp  &
           + rho_x1x3_n * sigma_x3_n * gamma_exp

    trivar_NLL_mean  &
    = ( 1.0 / sqrt( 2.0*pi_dp ) ) * ( - sigma_x1 )**alpha_exp  &
      * exp( mu_x2_n * beta_exp + mu_x3_n * gamma_exp )  &
      * exp( 0.5 *  &
             (   ( 1.0 - rho_x1x2_n**2.0 ) * sigma_x2_n**2.0 * beta_exp**2.0  &
               + ( 1.0 - rho_x1x3_n**2.0 ) * sigma_x3_n**2.0 * gamma_exp**2.0  &
               + 2.0 * ( rho_x2x3_n - rho_x1x2_n * rho_x1x3_n )  &
                     * sigma_x2_n * beta_exp * sigma_x3_n * gamma_exp  &
             )  &
           )  &
      * exp( 0.25 * s_cc**2.0 - ( mu_x1 / sigma_x1 ) * s_cc  &
             + 0.5 * ( mu_x1**2.0 / sigma_x1**2.0 ) )  &
      * gamma( alpha_exp + 1.0 ) * Dv_fnc( -(alpha_exp + 1.0), s_cc )

    return

  end function trivar_NLL_mean
  
  !=============================================================================
  function trivar_NLL_mean_const_x1( mu_x1, mu_x2_n, mu_x3_n, &
                                     sigma_x2_n, sigma_x3_n, rho_x2x3_n, &
                                     alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    !  Larson, V. E. and B. M. Griffin (2010)
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    double precision, intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2_n,    & ! Mean of ln x2 (ith PDF component)                     [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      sigma_x2_n, & ! Standard deviation of ln x2 (ith PDF component)       [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      rho_x2x3_n, & ! Correlation between ln x2 & ln x3 (ith PDF component) [-]
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x2                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x3                   [-]

    ! Return Variable
    double precision ::  &
      trivar_NLL_mean_const_x1

    trivar_NLL_mean_const_x1  &
    = mu_x1**alpha_exp  &
      * exp( mu_x2_n * beta_exp + mu_x3_n * gamma_exp  &
             + 0.5 * sigma_x2_n**2.0 * beta_exp**2.0  &
             + 0.5 * sigma_x3_n**2.0 * gamma_exp**2.0  &
             + rho_x2x3_n * sigma_x2_n * beta_exp * sigma_x3_n * gamma_exp )

    return

  end function trivar_NLL_mean_const_x1

  !=============================================================================
  function bivar_NL_mean( mu_x1, mu_x2_n, sigma_x1, sigma_x2_n, &
                          rho_x1x2_n, alpha_exp, beta_exp )

    ! Description:

    ! References:
    !  Larson, V. E. and B. M. Griffin (2010)
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp  ! Constant(s)

    use KK_utilities, only:  &
        Dv_fnc  ! Procedure(s)

    implicit none

    ! Input Variables
    double precision, intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2_n,    & ! Mean of ln x2 (ith PDF component)                     [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x2_n, & ! Standard deviation of ln x2 (ith PDF component)       [-]
      rho_x1x2_n, & ! Correlation between x1 and ln x2 (ith PDF component)  [-]
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp      ! Exponent beta, corresponding to x2                    [-]

    ! Return Variable
    double precision ::  &
      bivar_NL_mean

    ! Local Variable
    double precision ::  &
      s_c    !

    s_c = ( mu_x1 / sigma_x1 ) + rho_x1x2_n * sigma_x2_n * beta_exp

    bivar_NL_mean  &
    = ( 1.0 / sqrt( 2.0*pi_dp ) ) * sigma_x1**alpha_exp  &
      * exp( mu_x2_n * beta_exp  &
             + 0.5 * sigma_x2_n**2.0 * beta_exp**2.0 - 0.25 * s_c**2.0 )  &
      * gamma( alpha_exp + 1.0 ) * Dv_fnc( -(alpha_exp + 1.0), -s_c )

    return

  end function bivar_NL_mean

  !=============================================================================
  function bivar_NL_mean_const_x1( mu_x1, mu_x2_n, sigma_x2_n, &
                                   alpha_exp, beta_exp )


    ! Description:

    ! References:
    !  Larson, V. E. and B. M. Griffin (2010)
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    double precision, intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2_n,    & ! Mean of ln x2 (ith PDF component)                     [-]
      sigma_x2_n, & ! Standard deviation of ln x2 (ith PDF component)       [-]
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp      ! Exponent beta, corresponding to x2                    [-]

    ! Return Variable
    double precision :: &
      bivar_NL_mean_const_x1

    bivar_NL_mean_const_x1  &
    = mu_x1**alpha_exp  &
      * exp( mu_x2_n * beta_exp  &
             + 0.5 * sigma_x2_n**2.0 * beta_exp**2.0 )

    return

  end function bivar_NL_mean_const_x1

  !=============================================================================
  function bivar_LL_mean( mu_x1_n, mu_x2_n, sigma_x1_n, sigma_x2_n, &
                          rho_x1x2_n, alpha_exp, beta_exp )

    ! Description:

    ! References:
    !  Larson, V. E. and B. M. Griffin (2010)
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    double precision, intent(in) :: &
      mu_x1_n,    & ! Mean of ln x1 (ith PDF component)                     [-]
      mu_x2_n,    & ! Mean of ln x2 (ith PDF component)                     [-]
      sigma_x1_n, & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x2_n, & ! Standard deviation of ln x2 (ith PDF component)       [-]
      rho_x1x2_n, & ! Correlation between ln x1 & ln x2 (ith PDF component) [-]
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp      ! Exponent beta, corresponding to x2                    [-]

    ! Return Variable
    double precision ::  &
      bivar_LL_mean

    bivar_LL_mean  &
    = exp( mu_x1_n * alpha_exp + mu_x2_n * beta_exp  &
           + 0.5 * sigma_x1_n**2.0 * alpha_exp**2.0  &
           + 0.5 * sigma_x2_n**2.0 * beta_exp**2.0  &
           + rho_x1x2_n * sigma_x1_n * alpha_exp * sigma_x2_n * beta_exp )

    return

  end function bivar_LL_mean

!===============================================================================

end module PDF_integrals_means
