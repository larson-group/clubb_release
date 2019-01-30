! $Id$
!===============================================================================
module PDF_integrals_covars

  implicit none

  private

  public :: quadrivar_NNLL_covar, &
            quadrivar_NNLL_covar_const_x1, &
            quadrivar_NNLL_covar_const_x2, &
            quadrivar_NNLL_covar_const_x3, &
            quadrivar_NNLL_covar_const_x1x2, &
            quadrivar_NNLL_covar_const_x1x3, &
            quadrivar_NNLL_covar_const_x2x3, &
            quadrivar_NNLL_covar_const_x3x4, &
            quadrivar_NNLL_covar_cst_x1x2x3, &
            quadrivar_NNLL_covar_cst_x1x3x4, &
            quadrivar_NNLL_covar_cst_x2x3x4, &
            quadrivar_NNLL_covar_const_all, &
            trivar_NNL_covar, &
            trivar_NNL_covar_const_x1, &
            trivar_NNL_covar_const_x2, &
            trivar_NNL_covar_const_x3, &
            trivar_NNL_covar_const_x1x2, &
            trivar_NNL_covar_const_x1x3, &
            trivar_NNL_covar_const_x2x3, &
            trivar_NNL_covar_const_all

  contains

  !=============================================================================
  function quadrivar_NNLL_covar( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                 sigma_x1, sigma_x2, sigma_x3_n, sigma_x4_n, &
                                 rho_x1x2, rho_x1x3_n, rho_x1x4_n, &
                                 rho_x2x3_n, rho_x2x4_n, rho_x3x4_n, &
                                 x1_mean, x2_alpha_x3_beta_x4_gamma_mean, &
                                 alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    ! Eq. (J1) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S9) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp,         &  ! Constant(s)
        two_dp,        &
        one_dp,        &
        one_half_dp,   &
        one_fourth_dp

    use KK_utilities, only:  &
        Dv_fnc  ! Procedure(s)

    use parabolic, only:  &
        gamma  ! Procedure(s)

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
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

    real( kind = dp ), intent(in) :: &
      x1_mean,                        & ! Mean of x1 (overall)              [-]
      x2_alpha_x3_beta_x4_gamma_mean    ! Mean of x2^alpha x3^beta x4^gamma [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x3                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x4                   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      quadrivar_NNLL_covar

    ! Local Variable
    real( kind = dp ) ::  &
      s_cc    !

    s_cc = ( mu_x2 / sigma_x2 )  &
           + rho_x2x3_n * sigma_x3_n * beta_exp  &
           + rho_x2x4_n * sigma_x4_n * gamma_exp

    quadrivar_NNLL_covar  &
    = one_dp / sqrt( two_dp*pi_dp )  &
      * ( - sigma_x2 )**alpha_exp  &
      * exp( mu_x3_n * beta_exp + mu_x4_n * gamma_exp  &
             + one_half_dp * ( one_dp - rho_x2x3_n**2 )  &
                           * sigma_x3_n**2 * beta_exp**2  &
             + one_half_dp * ( one_dp - rho_x2x4_n**2 )  &
                           * sigma_x4_n**2 * gamma_exp**2  &
             + ( rho_x3x4_n - rho_x2x3_n * rho_x2x4_n )  &
                           * sigma_x3_n * beta_exp * sigma_x4_n * gamma_exp  &
           )  &
      * exp( one_fourth_dp * s_cc**2 - ( mu_x2 / sigma_x2 ) * s_cc  &
             + one_half_dp * ( mu_x2**2 / sigma_x2**2 ) )  &
      * ( - rho_x1x2 * sigma_x1 * gamma( alpha_exp + two_dp )  &
          * Dv_fnc( -(alpha_exp + two_dp), s_cc )  &
        + ( mu_x1 - x1_mean  &
            - ( mu_x2 / sigma_x2 ) * rho_x1x2 * sigma_x1  &
            + ( rho_x1x3_n - rho_x1x2 * rho_x2x3_n )  &
                     * sigma_x1 * sigma_x3_n * beta_exp  &
            + ( rho_x1x4_n - rho_x1x2 * rho_x2x4_n )  &
                     * sigma_x1 * sigma_x4_n * gamma_exp )  &
          * gamma( alpha_exp + one_dp )  &
          * Dv_fnc( -(alpha_exp + one_dp), s_cc )  &
        )  &
      - x2_alpha_x3_beta_x4_gamma_mean * ( mu_x1 - x1_mean )

    return

  end function quadrivar_NNLL_covar

  !=============================================================================
  function quadrivar_NNLL_covar_const_x1( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                          sigma_x2, sigma_x3_n, sigma_x4_n, &
                                          rho_x2x3_n, rho_x2x4_n, rho_x3x4_n, &
                                          x1_mean, &
                                          x2_alpha_x3_beta_x4_gamma_mean, &
                                          alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    ! Eq. (J2) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S10) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp,         &  ! Constant(s)
        two_dp,        &
        one_dp,        &
        one_half_dp,   &
        one_fourth_dp

    use KK_utilities, only:  &
        Dv_fnc  ! Procedure(s)

    use parabolic, only:  &
        gamma  ! Procedure(s)

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      mu_x4_n,    & ! Mean of ln x4 (ith PDF component)                     [-]
      sigma_x2,   & ! Standard deviation of x2 (ith PDF component)          [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      sigma_x4_n, & ! Standard deviation of ln x4 (ith PDF component)       [-]
      rho_x2x3_n, & ! Correlation between x2 and ln x3 (ith PDF component)  [-]
      rho_x2x4_n, & ! Correlation between x2 and ln x4 (ith PDF component)  [-]
      rho_x3x4_n    ! Correlation between ln x3 & ln x4 (ith PDF component) [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,                        & ! Mean of x1 (overall)              [-]
      x2_alpha_x3_beta_x4_gamma_mean    ! Mean of x2^alpha x3^beta x4^gamma [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x3                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x4                   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      quadrivar_NNLL_covar_const_x1

    ! Local Variable
    real( kind = dp ) ::  &
      s_cc    !

    s_cc = ( mu_x2 / sigma_x2 )  &
           + rho_x2x3_n * sigma_x3_n * beta_exp  &
           + rho_x2x4_n * sigma_x4_n * gamma_exp

    quadrivar_NNLL_covar_const_x1  &
    = ( one_dp / sqrt( two_dp*pi_dp ) ) * ( mu_x1 - x1_mean )  &
      * ( - sigma_x2 )**alpha_exp  &
      * exp( mu_x3_n * beta_exp + mu_x4_n * gamma_exp  &
             + one_half_dp * ( one_dp - rho_x2x3_n**2 )  &
                           * sigma_x3_n**2 * beta_exp**2  &
             + one_half_dp * ( one_dp - rho_x2x4_n**2 )  &
                           * sigma_x4_n**2 * gamma_exp**2  &
             + ( rho_x3x4_n - rho_x2x3_n * rho_x2x4_n )  &
                           * sigma_x3_n * beta_exp * sigma_x4_n * gamma_exp  &
           )  &
      * exp( one_fourth_dp * s_cc**2 - ( mu_x2 / sigma_x2 ) * s_cc  &
             + one_half_dp * ( mu_x2**2 / sigma_x2**2 ) )  &
      * gamma( alpha_exp + one_dp )  &
      * Dv_fnc( -(alpha_exp + one_dp), s_cc )  &
      - x2_alpha_x3_beta_x4_gamma_mean * ( mu_x1 - x1_mean )

    return

  end function quadrivar_NNLL_covar_const_x1

  !=============================================================================
  function quadrivar_NNLL_covar_const_x2( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                          sigma_x1, sigma_x3_n, sigma_x4_n, &
                                          rho_x1x3_n, rho_x1x4_n, rho_x3x4_n, &
                                          x1_mean, &
                                          x2_alpha_x3_beta_x4_gamma_mean, &
                                          alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    ! Eq. (J3) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S11) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one_half_dp,   &  ! Constant(s)
        zero_dp

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      mu_x4_n,    & ! Mean of ln x4 (ith PDF component)                     [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      sigma_x4_n, & ! Standard deviation of ln x4 (ith PDF component)       [-]
      rho_x1x3_n, & ! Correlation between x1 and ln x3 (ith PDF component)  [-]
      rho_x1x4_n, & ! Correlation between x1 and ln x4 (ith PDF component)  [-]
      rho_x3x4_n    ! Correlation between ln x3 & ln x4 (ith PDF component) [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,                        & ! Mean of x1 (overall)              [-]
      x2_alpha_x3_beta_x4_gamma_mean    ! Mean of x2^alpha x3^beta x4^gamma [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x3                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x4                   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      quadrivar_NNLL_covar_const_x2

    if ( mu_x2 <= zero_dp ) then

       quadrivar_NNLL_covar_const_x2  &
       = mu_x2**alpha_exp  &
         * ( mu_x1 - x1_mean  &
             + rho_x1x3_n * sigma_x1 * sigma_x3_n * beta_exp  &
             + rho_x1x4_n * sigma_x1 * sigma_x4_n * gamma_exp )  &
         * exp( mu_x3_n * beta_exp + mu_x4_n * gamma_exp  &
                + one_half_dp * sigma_x3_n**2 * beta_exp**2  &
                + one_half_dp * sigma_x4_n**2 * gamma_exp**2  &
                + rho_x3x4_n * sigma_x3_n * beta_exp  &
                             * sigma_x4_n * gamma_exp )  &
         - x2_alpha_x3_beta_x4_gamma_mean * ( mu_x1 - x1_mean )

    else ! mu_x2 > 0

       quadrivar_NNLL_covar_const_x2  &
       = - x2_alpha_x3_beta_x4_gamma_mean * ( mu_x1 - x1_mean )

    endif

    return

  end function quadrivar_NNLL_covar_const_x2

  !=============================================================================
  function quadrivar_NNLL_covar_const_x3( mu_x1, mu_x2, mu_x3, mu_x4_n, &
                                          sigma_x1, sigma_x2, sigma_x4_n, &
                                          rho_x1x2, rho_x1x4_n, rho_x2x4_n, &
                                          x1_mean, &
                                          x2_alpha_x3_beta_x4_gamma_mean, &
                                          alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    ! Eq. (J4) and Eq. (J5) of Griffin, B. M., 2016:  Improving the
    ! Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
    ! Effects Using a Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S12) and Eq. (S13) of Griffin, B. M. and V. E. Larson, 2016:
    ! Supplement of Parameterizing microphysical effects on variances and
    ! covariances of moisture and heat content using a multivariate probability
    ! density function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9,
    ! 11, doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp,         &  ! Constant(s)
        two_dp,        &
        one_dp,        &
        one_half_dp,   &
        one_fourth_dp

    use KK_utilities, only:  &
        Dv_fnc  ! Procedure(s)

    use parabolic, only:  &
        gamma  ! Procedure(s)

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3,      & ! Mean of x3 (ith PDF component)                        [-]
      mu_x4_n,    & ! Mean of ln x4 (ith PDF component)                     [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x2,   & ! Standard deviation of x2 (ith PDF component)          [-]
      sigma_x4_n, & ! Standard deviation of ln x4 (ith PDF component)       [-]
      rho_x1x2,   & ! Correlation between x1 and x2 (ith PDF component)     [-]
      rho_x1x4_n, & ! Correlation between x1 and ln x4 (ith PDF component)  [-]
      rho_x2x4_n    ! Correlation between x2 and ln x4 (ith PDF component)  [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,                        & ! Mean of x1 (overall)              [-]
      x2_alpha_x3_beta_x4_gamma_mean    ! Mean of x2^alpha x3^beta x4^gamma [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x3                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x4                   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      quadrivar_NNLL_covar_const_x3

    ! Local Variable
    real( kind = dp ) ::  &
      s_cc    !

    s_cc = ( mu_x2 / sigma_x2 ) + rho_x2x4_n * sigma_x4_n * gamma_exp

    quadrivar_NNLL_covar_const_x3  &
    = one_dp / sqrt( two_dp*pi_dp )  &
      * ( - sigma_x2 )**alpha_exp  &
      * mu_x3**beta_exp  &
      * exp( mu_x4_n * gamma_exp  &
             + one_half_dp * sigma_x4_n**2 * gamma_exp**2  &
             - one_fourth_dp * s_cc**2  &
           )  &
      * ( - rho_x1x2 * sigma_x1 * gamma( alpha_exp + two_dp )  &
          * Dv_fnc( -(alpha_exp + two_dp), s_cc )  &
        + ( mu_x1 - x1_mean  &
            - ( mu_x2 / sigma_x2 ) * rho_x1x2 * sigma_x1  &
            + ( rho_x1x4_n - rho_x1x2 * rho_x2x4_n )  &
              * sigma_x1 * sigma_x4_n * gamma_exp )  &
          * gamma( alpha_exp + one_dp )  &
          * Dv_fnc( -(alpha_exp + one_dp), s_cc )  &
        )  &
      - x2_alpha_x3_beta_x4_gamma_mean * ( mu_x1 - x1_mean )

    return

  end function quadrivar_NNLL_covar_const_x3

  !=============================================================================
  function quadrivar_NNLL_covar_const_x1x2( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                            sigma_x3_n, sigma_x4_n, &
                                            rho_x3x4_n, x1_mean, &
                                            x2_alpha_x3_beta_x4_gamma_mean, &
                                            alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    ! Eq. (J6) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S14) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one_half_dp,   &  ! Constant(s)
        zero_dp

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      mu_x4_n,    & ! Mean of ln x4 (ith PDF component)                     [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      sigma_x4_n, & ! Standard deviation of ln x4 (ith PDF component)       [-]
      rho_x3x4_n    ! Correlation between ln x3 & ln x4 (ith PDF component) [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,                        & ! Mean of x1 (overall)              [-]
      x2_alpha_x3_beta_x4_gamma_mean    ! Mean of x2^alpha x3^beta x4^gamma [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x3                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x4                   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      quadrivar_NNLL_covar_const_x1x2

    if ( mu_x2 <= zero_dp ) then

       quadrivar_NNLL_covar_const_x1x2  &
       = ( mu_x1 - x1_mean )  &
         * mu_x2**alpha_exp  &
         * exp( mu_x3_n * beta_exp + mu_x4_n * gamma_exp  &
                + one_half_dp * sigma_x3_n**2 * beta_exp**2  &
                + one_half_dp * sigma_x4_n**2 * gamma_exp**2  &
                + rho_x3x4_n * sigma_x3_n * beta_exp  &
                             * sigma_x4_n * gamma_exp )  &
         - x2_alpha_x3_beta_x4_gamma_mean * ( mu_x1 - x1_mean )

    else ! mu_x2 > 0

       quadrivar_NNLL_covar_const_x1x2  &
       = - x2_alpha_x3_beta_x4_gamma_mean * ( mu_x1 - x1_mean )

    endif

    return

  end function quadrivar_NNLL_covar_const_x1x2

  !=============================================================================
  function quadrivar_NNLL_covar_const_x1x3( mu_x1, mu_x2, mu_x3, mu_x4_n, &
                                            sigma_x2, sigma_x4_n, rho_x2x4_n, &
                                            x1_mean, &
                                            x2_alpha_x3_beta_x4_gamma_mean, &
                                            alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    ! Eq. (J7) and Eq. (J8) of Griffin, B. M., 2016:  Improving the
    ! Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
    ! Effects Using a Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S15) and Eq. (S16) of Griffin, B. M. and V. E. Larson, 2016:
    ! Supplement of Parameterizing microphysical effects on variances and
    ! covariances of moisture and heat content using a multivariate probability
    ! density function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9,
    ! 11, doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp,         &  ! Constant(s)
        two_dp,        &
        one_dp,        &
        one_half_dp,   &
        one_fourth_dp

    use KK_utilities, only:  &
        Dv_fnc  ! Procedure(s)

    use parabolic, only:  &
        gamma  ! Procedure(s)

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3,      & ! Mean of x3 (ith PDF component)                        [-]
      mu_x4_n,    & ! Mean of ln x4 (ith PDF component)                     [-]
      sigma_x2,   & ! Standard deviation of x2 (ith PDF component)          [-]
      sigma_x4_n, & ! Standard deviation of ln x4 (ith PDF component)       [-]
      rho_x2x4_n    ! Correlation between x2 and ln x4 (ith PDF component)  [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,                        & ! Mean of x1 (overall)              [-]
      x2_alpha_x3_beta_x4_gamma_mean    ! Mean of x2^alpha x3^beta x4^gamma [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x3                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x4                   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      quadrivar_NNLL_covar_const_x1x3

    ! Local Variable
    real( kind = dp ) ::  &
      s_cc    !

    s_cc = ( mu_x2 / sigma_x2 ) + rho_x2x4_n * sigma_x4_n * gamma_exp

    quadrivar_NNLL_covar_const_x1x3  &
    = ( one_dp / sqrt( two_dp*pi_dp ) ) * ( mu_x1 - x1_mean )  &
      * ( - sigma_x2 )**alpha_exp  &
      * mu_x3**beta_exp  &
      * exp( mu_x4_n * gamma_exp  &
             + one_half_dp * sigma_x4_n**2 * gamma_exp**2  &
             - one_fourth_dp * s_cc**2 )  &
      * gamma( alpha_exp + one_dp )  &
      * Dv_fnc( -(alpha_exp + one_dp), s_cc )  &
      - x2_alpha_x3_beta_x4_gamma_mean * ( mu_x1 - x1_mean )

    return

  end function quadrivar_NNLL_covar_const_x1x3

  !=============================================================================
  function quadrivar_NNLL_covar_const_x2x3( mu_x1, mu_x2, mu_x3, mu_x4_n, &
                                            sigma_x1, sigma_x4_n, rho_x1x4_n, &
                                            x1_mean, &
                                            x2_alpha_x3_beta_x4_gamma_mean, &
                                            alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    ! Eq. (J9) and Eq. (J10) of Griffin, B. M., 2016:  Improving the
    ! Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
    ! Effects Using a Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S17) and Eq. (S18) of Griffin, B. M. and V. E. Larson, 2016:
    ! Supplement of Parameterizing microphysical effects on variances and
    ! covariances of moisture and heat content using a multivariate probability
    ! density function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9,
    ! 11, doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one_half_dp,   &  ! Constant(s)
        zero_dp

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3,      & ! Mean of x3 (ith PDF component)                        [-]
      mu_x4_n,    & ! Mean of ln x4 (ith PDF component)                     [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x4_n, & ! Standard deviation of ln x4 (ith PDF component)       [-]
      rho_x1x4_n    ! Correlation between x1 and ln x4 (ith PDF component)  [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,                        & ! Mean of x1 (overall)              [-]
      x2_alpha_x3_beta_x4_gamma_mean    ! Mean of x2^alpha x3^beta x4^gamma [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x3                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x4                   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      quadrivar_NNLL_covar_const_x2x3

    if ( mu_x2 <= zero_dp ) then

       quadrivar_NNLL_covar_const_x2x3  &
       = mu_x2**alpha_exp * mu_x3**beta_exp  &
         * ( mu_x1 - x1_mean  &
             + rho_x1x4_n * sigma_x1 * sigma_x4_n * gamma_exp )  &
         * exp( mu_x4_n * gamma_exp  &
                + one_half_dp * sigma_x4_n**2 * gamma_exp**2 )  &
         - x2_alpha_x3_beta_x4_gamma_mean * ( mu_x1 - x1_mean )

    else ! mu_x2 > 0

       quadrivar_NNLL_covar_const_x2x3  &
       = - x2_alpha_x3_beta_x4_gamma_mean * ( mu_x1 - x1_mean )

    endif

    return

  end function quadrivar_NNLL_covar_const_x2x3

  !=============================================================================
  function quadrivar_NNLL_covar_const_x3x4( mu_x1, mu_x2, mu_x3, mu_x4, &
                                            sigma_x1, sigma_x2, rho_x1x2, &
                                            x1_mean, &
                                            x2_alpha_x3_beta_x4_gamma_mean, &
                                            alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    ! Eq. (J11) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S19) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp,         &  ! Constant(s)
        two_dp,        &
        one_dp,        &
        one_fourth_dp

    use KK_utilities, only:  &
        Dv_fnc  ! Procedure(s)

    use parabolic, only:  &
        gamma  ! Procedure(s)

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,    & ! Mean of x1 (ith PDF component)                          [-]
      mu_x2,    & ! Mean of x2 (ith PDF component)                          [-]
      mu_x3,    & ! Mean of x3 (ith PDF component)                          [-]
      mu_x4,    & ! Mean of x4 (ith PDF component)                          [-]
      sigma_x1, & ! Standard deviation of x1 (ith PDF component)            [-]
      sigma_x2, & ! Standard deviation of x2 (ith PDF component)            [-]
      rho_x1x2    ! Correlation between x1 and x2 (ith PDF component)       [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,                        & ! Mean of x1 (overall)              [-]
      x2_alpha_x3_beta_x4_gamma_mean    ! Mean of x2^alpha x3^beta x4^gamma [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x3                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x4                   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      quadrivar_NNLL_covar_const_x3x4

    quadrivar_NNLL_covar_const_x3x4  &
    = one_dp / sqrt( two_dp*pi_dp )  &
      * ( - sigma_x2 )**alpha_exp  &
      * mu_x3**beta_exp * mu_x4**gamma_exp  &
      * exp( - one_fourth_dp * ( mu_x2**2 / sigma_x2**2 ) )  &
      * ( - rho_x1x2 * sigma_x1 * gamma( alpha_exp + two_dp )  &
          * Dv_fnc( -(alpha_exp + two_dp), ( mu_x2 / sigma_x2 ) )  &
        + ( mu_x1 - x1_mean  &
            - ( mu_x2 / sigma_x2 ) * rho_x1x2 * sigma_x1 )  &
          * gamma( alpha_exp + one_dp )  &
          * Dv_fnc( -(alpha_exp + one_dp), ( mu_x2 / sigma_x2 ) )  &
        )  &
      - x2_alpha_x3_beta_x4_gamma_mean * ( mu_x1 - x1_mean )

    return

  end function quadrivar_NNLL_covar_const_x3x4

  !=============================================================================
  function quadrivar_NNLL_covar_cst_x1x2x3( mu_x1, mu_x2, mu_x3, mu_x4_n, &
                                            sigma_x4_n, x1_mean, &
                                            x2_alpha_x3_beta_x4_gamma_mean, &
                                            alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    ! Eq. (J12) and Eq. (J13) of Griffin, B. M., 2016:  Improving the
    ! Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
    ! Effects Using a Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S20) and Eq. (S21) of Griffin, B. M. and V. E. Larson, 2016:
    ! Supplement of Parameterizing microphysical effects on variances and
    ! covariances of moisture and heat content using a multivariate probability
    ! density function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9,
    ! 11, doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one_half_dp,   &  ! Constant(s)
        zero_dp

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3,      & ! Mean of x3 (ith PDF component)                        [-]
      mu_x4_n,    & ! Mean of ln x4 (ith PDF component)                     [-]
      sigma_x4_n    ! Standard deviation of ln x4 (ith PDF component)       [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,                        & ! Mean of x1 (overall)              [-]
      x2_alpha_x3_beta_x4_gamma_mean    ! Mean of x2^alpha x3^beta x4^gamma [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x3                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x4                   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      quadrivar_NNLL_covar_cst_x1x2x3

    if ( mu_x2 <= zero_dp ) then

       quadrivar_NNLL_covar_cst_x1x2x3  &
       = ( mu_x1 - x1_mean )  &
         * mu_x2**alpha_exp * mu_x3**beta_exp  &
         * exp( mu_x4_n * gamma_exp  &
                + one_half_dp * sigma_x4_n**2 * gamma_exp**2 )  &
         - x2_alpha_x3_beta_x4_gamma_mean * ( mu_x1 - x1_mean )

    else ! mu_x2 > 0

       quadrivar_NNLL_covar_cst_x1x2x3  &
       = - x2_alpha_x3_beta_x4_gamma_mean * ( mu_x1 - x1_mean )

    endif

    return

  end function quadrivar_NNLL_covar_cst_x1x2x3

  !=============================================================================
  function quadrivar_NNLL_covar_cst_x1x3x4( mu_x1, mu_x2, mu_x3, mu_x4, &
                                            sigma_x2, x1_mean, &
                                            x2_alpha_x3_beta_x4_gamma_mean, &
                                            alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    ! Eq. (J14) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S22) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp,         &  ! Constant(s)
        two_dp,        &
        one_dp,        &
        one_fourth_dp

    use KK_utilities, only:  &
        Dv_fnc  ! Procedure(s)

    use parabolic, only:  &
        gamma  ! Procedure(s)

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,    & ! Mean of x1 (ith PDF component)                          [-]
      mu_x2,    & ! Mean of x2 (ith PDF component)                          [-]
      mu_x3,    & ! Mean of x3 (ith PDF component)                          [-]
      mu_x4,    & ! Mean of x4 (ith PDF component)                          [-]
      sigma_x2    ! Standard deviation of x2 (ith PDF component)            [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,                        & ! Mean of x1 (overall)              [-]
      x2_alpha_x3_beta_x4_gamma_mean    ! Mean of x2^alpha x3^beta x4^gamma [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x3                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x4                   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      quadrivar_NNLL_covar_cst_x1x3x4

    quadrivar_NNLL_covar_cst_x1x3x4  &
    = ( one_dp / sqrt( two_dp*pi_dp ) ) * ( mu_x1 - x1_mean )  &
      * ( - sigma_x2 )**alpha_exp  &
      * mu_x3**beta_exp * mu_x4**gamma_exp  &
      * exp( - one_fourth_dp * ( mu_x2**2 / sigma_x2**2 ) )  &
      * gamma( alpha_exp + one_dp )  &
      * Dv_fnc( -(alpha_exp + one_dp), ( mu_x2 / sigma_x2 ) )  &
      - x2_alpha_x3_beta_x4_gamma_mean * ( mu_x1 - x1_mean )

    return

  end function quadrivar_NNLL_covar_cst_x1x3x4

  !=============================================================================
  function quadrivar_NNLL_covar_cst_x2x3x4( mu_x1, mu_x2, mu_x3, mu_x4, &
                                            x1_mean, &
                                            x2_alpha_x3_beta_x4_gamma_mean, &
                                            alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    ! Eq. (J15) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S23) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        zero_dp    ! Constant(s)

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1, & ! Mean of x1 (ith PDF component)                             [-]
      mu_x2, & ! Mean of x2 (ith PDF component)                             [-]
      mu_x3, & ! Mean of x3 (ith PDF component)                             [-]
      mu_x4    ! Mean of x4 (ith PDF component)                             [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,                        & ! Mean of x1 (overall)              [-]
      x2_alpha_x3_beta_x4_gamma_mean    ! Mean of x2^alpha x3^beta x4^gamma [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x3                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x4                   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      quadrivar_NNLL_covar_cst_x2x3x4

    if ( mu_x2 <= zero_dp ) then

       quadrivar_NNLL_covar_cst_x2x3x4  &
       = ( mu_x1 - x1_mean )  &
         * ( mu_x2**alpha_exp * mu_x3**beta_exp * mu_x4**gamma_exp  &
             - x2_alpha_x3_beta_x4_gamma_mean )

    else ! mu_x2 > 0

       quadrivar_NNLL_covar_cst_x2x3x4  &
       = - x2_alpha_x3_beta_x4_gamma_mean * ( mu_x1 - x1_mean )

    endif

    return

  end function quadrivar_NNLL_covar_cst_x2x3x4

  !=============================================================================
  function quadrivar_NNLL_covar_const_all( mu_x1, mu_x2, mu_x3, mu_x4, &
                                           x1_mean, &
                                           x2_alpha_x3_beta_x4_gamma_mean, &
                                           alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    ! Eq. (J16) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S24) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        zero_dp    ! Constant(s) 

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1, & ! Mean of x1 (ith PDF component)                             [-]
      mu_x2, & ! Mean of x2 (ith PDF component)                             [-]
      mu_x3, & ! Mean of x3 (ith PDF component)                             [-]
      mu_x4    ! Mean of x4 (ith PDF component)                             [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,                        & ! Mean of x1 (overall)              [-]
      x2_alpha_x3_beta_x4_gamma_mean    ! Mean of x2^alpha x3^beta x4^gamma [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x3                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x4                   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      quadrivar_NNLL_covar_const_all

    if ( mu_x2 <= zero_dp ) then

       quadrivar_NNLL_covar_const_all  &
       = ( mu_x1 - x1_mean )  &
         * ( mu_x2**alpha_exp * mu_x3**beta_exp * mu_x4**gamma_exp  &
             - x2_alpha_x3_beta_x4_gamma_mean )

    else ! mu_x2 > 0

       quadrivar_NNLL_covar_const_all  &
       = - x2_alpha_x3_beta_x4_gamma_mean * ( mu_x1 - x1_mean )

    endif

    return

  end function quadrivar_NNLL_covar_const_all

  !=============================================================================
  function trivar_NNL_covar( mu_x1, mu_x2, mu_x3_n, &
                             sigma_x1, sigma_x2, sigma_x3_n, &
                             rho_x1x2, rho_x1x3_n, rho_x2x3_n, &
                             x1_mean, x2_alpha_x3_beta_mean, &
                             alpha_exp, beta_exp )

    ! Description:

    ! References:
    ! Eq. (J17) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S25) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp,         &  ! Constant(s)
        two_dp,        &
        one_dp,        &
        one_half_dp,   &
        one_fourth_dp

    use KK_utilities, only:  &
        Dv_fnc  ! Procedure(s)

    use parabolic, only:  &
        gamma  ! Procedure(s)

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x2,   & ! Standard deviation of x2 (ith PDF component)          [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      rho_x1x2,   & ! Correlation between x1 and x2 (ith PDF component)     [-]
      rho_x1x3_n, & ! Correlation between x1 and ln x3 (ith PDF component)  [-]
      rho_x2x3_n    ! Correlation between x2 and ln x3 (ith PDF component)  [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,               & ! Mean of x1 (overall)                       [-]
      x2_alpha_x3_beta_mean    ! Mean of x2^alpha x3^beta                   [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp      ! Exponent beta, corresponding to x3                    [-]

    ! Return Variable
    real( kind = dp ) ::  &
      trivar_NNL_covar

    ! Local Variable
    real( kind = dp ) ::  &
      s_c    !

    s_c = ( mu_x2 / sigma_x2 ) + rho_x2x3_n * sigma_x3_n * beta_exp

    trivar_NNL_covar  &
    = one_dp / sqrt( two_dp*pi_dp )  &
      * sigma_x2**alpha_exp  &
      * exp( mu_x3_n * beta_exp  &
             + one_half_dp * sigma_x3_n**2 * beta_exp**2  &
             - one_fourth_dp * s_c**2 )  &
      * ( rho_x1x2 * sigma_x1 * gamma( alpha_exp + two_dp )  &
          * Dv_fnc( -(alpha_exp + two_dp), -s_c )  &
        + ( mu_x1 - x1_mean  &
            - ( mu_x2 / sigma_x2 ) * rho_x1x2 * sigma_x1  &
            + ( rho_x1x3_n - rho_x1x2 * rho_x2x3_n )  &
                     * sigma_x1 * sigma_x3_n * beta_exp )  &
          * gamma( alpha_exp + one_dp )  &
          * Dv_fnc( -(alpha_exp + one_dp), -s_c )  &
        )  &
      - x2_alpha_x3_beta_mean * ( mu_x1 - x1_mean )

    return

  end function trivar_NNL_covar

  !=============================================================================
  function trivar_NNL_covar_const_x1( mu_x1, mu_x2, mu_x3_n, &
                                      sigma_x2, sigma_x3_n, rho_x2x3_n, &
                                      x1_mean, x2_alpha_x3_beta_mean, &
                                      alpha_exp, beta_exp )

    ! Description:

    ! References:
    ! Eq. (J18) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S26) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp,         &  ! Constant(s)
        two_dp,        &
        one_dp,        &
        one_half_dp,   &
        one_fourth_dp

    use KK_utilities, only:  &
        Dv_fnc  ! Procedure(s)

    use parabolic, only:  &
        gamma  ! Procedure(s)

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      sigma_x2,   & ! Standard deviation of x2 (ith PDF component)          [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      rho_x2x3_n    ! Correlation between x2 and ln x3 (ith PDF component)  [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,               & ! Mean of x1 (overall)                       [-]
      x2_alpha_x3_beta_mean    ! Mean of x2^alpha x3^beta                   [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp      ! Exponent beta, corresponding to x3                    [-]

    ! Return Variable
    real( kind = dp ) ::  &
      trivar_NNL_covar_const_x1

    ! Local Variable
    real( kind = dp ) ::  &
      s_c    !

    s_c = ( mu_x2 / sigma_x2 ) + rho_x2x3_n * sigma_x3_n * beta_exp
   
    trivar_NNL_covar_const_x1  &
    = ( one_dp / sqrt( two_dp*pi_dp ) ) * ( mu_x1 - x1_mean )  &
      * sigma_x2**alpha_exp  &
      * exp( mu_x3_n * beta_exp  & 
             + one_half_dp * sigma_x3_n**2 * beta_exp**2  &
             - one_fourth_dp * s_c**2 )  &
      * gamma( alpha_exp + one_dp ) * Dv_fnc( -(alpha_exp + one_dp), -s_c )  &
      - x2_alpha_x3_beta_mean * ( mu_x1 - x1_mean )

    return

  end function trivar_NNL_covar_const_x1

  !=============================================================================
  function trivar_NNL_covar_const_x2( mu_x1, mu_x2, mu_x3_n, &
                                      sigma_x1, sigma_x3_n, rho_x1x3_n, &
                                      x1_mean, x2_alpha_x3_beta_mean, &
                                      alpha_exp, beta_exp )

    ! Description:

    ! References:
    ! Eq. (J19) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S27) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one_half_dp,   &  ! Constant(s)
        zero_dp

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      rho_x1x3_n    ! Correlation between x1 and ln x3 (ith PDF component)  [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,               & ! Mean of x1 (overall)                       [-]
      x2_alpha_x3_beta_mean    ! Mean of x2^alpha x3^beta                   [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp      ! Exponent beta, corresponding to x3                    [-]

    ! Return Variable
    real( kind = dp ) ::  &
      trivar_NNL_covar_const_x2

    if ( mu_x2 >= zero_dp ) then

       trivar_NNL_covar_const_x2  &
       = mu_x2**alpha_exp  &
         * ( mu_x1 - x1_mean  &
             + rho_x1x3_n * sigma_x1 * sigma_x3_n * beta_exp )  &
         * exp( mu_x3_n * beta_exp  &
                + one_half_dp * sigma_x3_n**2 * beta_exp**2 )  &
         - x2_alpha_x3_beta_mean * ( mu_x1 - x1_mean )

    else ! mu_x2 < 0

       trivar_NNL_covar_const_x2  &
       = - x2_alpha_x3_beta_mean * ( mu_x1 - x1_mean )

    endif
  
    return

  end function trivar_NNL_covar_const_x2

  !=============================================================================
  function trivar_NNL_covar_const_x3( mu_x1, mu_x2, mu_x3, &
                                      sigma_x1, sigma_x2, rho_x1x2, &
                                      x1_mean, x2_alpha_x3_beta_mean, &
                                      alpha_exp, beta_exp )

    ! Description:

    ! References:
    ! Eq. (J20) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S28) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp,         &  ! Constant(s)
        two_dp,        &
        one_dp,        &
        one_fourth_dp

    use KK_utilities, only:  &
        Dv_fnc  ! Procedure(s)

    use parabolic, only:  &
        gamma  ! Procedure(s)

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,    & ! Mean of x1 (ith PDF component)                          [-]
      mu_x2,    & ! Mean of x2 (ith PDF component)                          [-]
      mu_x3,    & ! Mean of x3 (ith PDF component)                          [-]
      sigma_x1, & ! Standard deviation of x1 (ith PDF component)            [-]
      sigma_x2, & ! Standard deviation of x2 (ith PDF component)            [-]
      rho_x1x2    ! Correlation between x1 and x2 (ith PDF component)       [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,               & ! Mean of x1 (overall)                       [-]
      x2_alpha_x3_beta_mean    ! Mean of x2^alpha x3^beta                   [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp      ! Exponent beta, corresponding to x3                    [-]

    ! Return Variable
    real( kind = dp ) ::  &
      trivar_NNL_covar_const_x3

    trivar_NNL_covar_const_x3  &
    = one_dp / sqrt( two_dp*pi_dp )  &
      * sigma_x2**alpha_exp * mu_x3**beta_exp  &
      * exp( - one_fourth_dp * ( mu_x2**2 / sigma_x2**2 ) )  &
      * ( rho_x1x2 * sigma_x1 * gamma( alpha_exp + two_dp )  &
          * Dv_fnc( -(alpha_exp + two_dp), -( mu_x2 / sigma_x2 ) )  &
        + ( mu_x1 - x1_mean  &
            - ( mu_x2 / sigma_x2 ) * rho_x1x2 * sigma_x1 )  &
          * gamma( alpha_exp + one_dp )  &
          * Dv_fnc( -(alpha_exp + one_dp), -( mu_x2 / sigma_x2 ) )  &
        )  &
      - x2_alpha_x3_beta_mean * ( mu_x1 - x1_mean )

    return

  end function trivar_NNL_covar_const_x3

  !=============================================================================
  function trivar_NNL_covar_const_x1x2( mu_x1, mu_x2, mu_x3_n, sigma_x3_n, &
                                        x1_mean, x2_alpha_x3_beta_mean, &
                                        alpha_exp, beta_exp )

    ! Description:

    ! References:
    ! Eq. (J21) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S29) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one_half_dp,   &  ! Constant(s)
        zero_dp

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      sigma_x3_n    ! Standard deviation of ln x3 (ith PDF component)       [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,               & ! Mean of x1 (overall)                       [-]
      x2_alpha_x3_beta_mean    ! Mean of x2^alpha x3^beta                   [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp      ! Exponent beta, corresponding to x3                    [-]

    ! Return Variable
    real( kind = dp ) ::  &
      trivar_NNL_covar_const_x1x2

    if ( mu_x2 >= zero_dp ) then

       trivar_NNL_covar_const_x1x2  &
       = mu_x2**alpha_exp  &
         * ( mu_x1 - x1_mean )  &
         * exp( mu_x3_n * beta_exp  &
                + one_half_dp * sigma_x3_n**2 * beta_exp**2 )  &
         - x2_alpha_x3_beta_mean * ( mu_x1 - x1_mean )

    else ! mu_x2 < 0

       trivar_NNL_covar_const_x1x2  &
       = - x2_alpha_x3_beta_mean * ( mu_x1 - x1_mean )

    endif
  
    return

  end function trivar_NNL_covar_const_x1x2

  !=============================================================================
  function trivar_NNL_covar_const_x1x3( mu_x1, mu_x2, mu_x3, sigma_x2, &
                                        x1_mean, x2_alpha_x3_beta_mean, &
                                        alpha_exp, beta_exp )

    ! Description:

    ! References:
    ! Eq. (J22) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S30) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp,         &  ! Constant(s)
        two_dp,        &
        one_dp,        &
        one_fourth_dp

    use KK_utilities, only:  &
        Dv_fnc  ! Procedure(s)

    use parabolic, only:  &
        gamma  ! Procedure(s)

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,    & ! Mean of x1 (ith PDF component)                          [-]
      mu_x2,    & ! Mean of x2 (ith PDF component)                          [-]
      mu_x3,    & ! Mean of x3 (ith PDF component)                          [-]
      sigma_x2    ! Standard deviation of x2 (ith PDF component)            [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,               & ! Mean of x1 (overall)                       [-]
      x2_alpha_x3_beta_mean    ! Mean of x2^alpha x3^beta                   [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp      ! Exponent beta, corresponding to x3                    [-]

    ! Return Variable
    real( kind = dp ) ::  &
      trivar_NNL_covar_const_x1x3   

    trivar_NNL_covar_const_x1x3  &
    = ( one_dp / sqrt( two_dp*pi_dp ) ) * ( mu_x1 - x1_mean )  &
      * sigma_x2**alpha_exp * mu_x3**beta_exp  &
      * exp( - one_fourth_dp * ( mu_x2**2 / sigma_x2**2 ) )  &
      * gamma( alpha_exp + one_dp )  &
      * Dv_fnc( -(alpha_exp + one_dp), -( mu_x2 / sigma_x2 ) )  &
      - x2_alpha_x3_beta_mean * ( mu_x1 - x1_mean )

    return

  end function trivar_NNL_covar_const_x1x3

  !=============================================================================
  function trivar_NNL_covar_const_x2x3( mu_x1, mu_x2, mu_x3, &
                                        x1_mean, x2_alpha_x3_beta_mean, &
                                        alpha_exp, beta_exp )

    ! Description:

    ! References:
    ! Eq. (J23) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S31) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        zero_dp  ! Constant(s)

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1, & ! Mean of x1 (ith PDF component)                             [-]
      mu_x2, & ! Mean of x2 (ith PDF component)                             [-]
      mu_x3    ! Mean of x3 (ith PDF component)                             [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,               & ! Mean of x1 (overall)                       [-]
      x2_alpha_x3_beta_mean    ! Mean of x2^alpha x3^beta                   [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp      ! Exponent beta, corresponding to x3                    [-]

    ! Return Variable
    real( kind = dp ) ::  &
      trivar_NNL_covar_const_x2x3

    if ( mu_x2 >= zero_dp ) then

       trivar_NNL_covar_const_x2x3  &
       = ( mu_x1 - x1_mean ) &
         * ( mu_x2**alpha_exp * mu_x3**beta_exp - x2_alpha_x3_beta_mean ) 

    else ! mu_x2 < 0

       trivar_NNL_covar_const_x2x3  &
       = - x2_alpha_x3_beta_mean * ( mu_x1 - x1_mean )

    endif
  
    return

  end function trivar_NNL_covar_const_x2x3

!=============================================================================
  function trivar_NNL_covar_const_all( mu_x1, mu_x2, mu_x3, &
                                       x1_mean, x2_alpha_x3_beta_mean, &
                                       alpha_exp, beta_exp )

    ! Description:

    ! References:
    ! Eq. (J24) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S32) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        zero_dp

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1, & ! Mean of x1 (ith PDF component)                             [-]
      mu_x2, & ! Mean of x2 (ith PDF component)                             [-]
      mu_x3    ! Mean of x3 (ith PDF component)                             [-]

    real( kind = dp ), intent(in) :: &
      x1_mean,               & ! Mean of x1 (overall)                       [-]
      x2_alpha_x3_beta_mean    ! Mean of x2^alpha x3^beta                   [-]
    
    real( kind = dp ), intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp      ! Exponent beta, corresponding to x3                    [-]

    ! Return Variable
    real( kind = dp ) ::  &
      trivar_NNL_covar_const_all

    if ( mu_x2 >= zero_dp ) then

       trivar_NNL_covar_const_all  &
       = ( mu_x1 - x1_mean )  &
         * ( mu_x2**alpha_exp * mu_x3**beta_exp - x2_alpha_x3_beta_mean )

    else ! mu_x2 < 0

       trivar_NNL_covar_const_all  &
       = - x2_alpha_x3_beta_mean * ( mu_x1 - x1_mean )

    endif
  
    return

  end function trivar_NNL_covar_const_all

!===============================================================================

end module PDF_integrals_covars
