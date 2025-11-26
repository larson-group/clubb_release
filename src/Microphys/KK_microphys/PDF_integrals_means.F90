! $Id$
!===============================================================================
module PDF_integrals_means

  implicit none

  private

  public :: trivar_NLL_mean,            &
            trivar_NLL_mean_const_x1,   &
            trivar_NLL_mean_const_x2,   &
            trivar_NLL_mean_const_x1x2, &
            trivar_NLL_mean_const_x2x3, &
            trivar_NLL_mean_const_all,  &
            bivar_NL_mean,              &
            bivar_NL_mean_const_x1,     &
            bivar_NL_mean_const_x2,     &
            bivar_NL_mean_const_all,    &
            bivar_LL_mean,              &
            bivar_LL_mean_const_x1,     &
            bivar_LL_mean_const_all

  contains

  !=============================================================================
  function trivar_NLL_mean( mu_x1, mu_x2_n, mu_x3_n, &
                            sigma_x1, sigma_x2_n, sigma_x3_n, &
                            rho_x1x2_n, rho_x1x3_n, rho_x2x3_n, &
                            alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    ! Eq. (60) and Eq. (61) of Larson, V. E. and B. M. Griffin, 2013:  Analytic
    ! upscaling of a local microphysics scheme. Part I: Derivation.
    ! Q. J. Roy. Meteorol. Soc., 139, 670, 46--57,
    ! doi:http://dx.doi.org/10.1002/qj.1967.
    !
    ! Eq. (C24) and Eq. (J25) of Griffin, B. M., 2016:  Improving the
    ! Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
    ! Effects Using a Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S24) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! A new subgrid-scale representation of hydrometeor fields using a
    ! multivariate PDF.  Geosci. Model Dev., 9, 6,
    ! doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
    !
    ! Eq. (S33) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp,         & ! Constant(s)
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
    real( kind = dp ) ::  &
      trivar_NLL_mean

    ! Local Variable
    real( kind = dp ) ::  &
      s_cc    !


    s_cc = ( mu_x1 / sigma_x1 )  &
           + rho_x1x2_n * sigma_x2_n * beta_exp  &
           + rho_x1x3_n * sigma_x3_n * gamma_exp

    trivar_NLL_mean  &
    = ( one_dp / sqrt( two_dp*pi_dp ) ) * ( - sigma_x1 )**alpha_exp  &
      * exp( mu_x2_n * beta_exp + mu_x3_n * gamma_exp )  &
      * exp( one_half_dp *  &
             (   ( one_dp - rho_x1x2_n**2 ) * sigma_x2_n**2 * beta_exp**2  &
               + ( one_dp - rho_x1x3_n**2 ) * sigma_x3_n**2 * gamma_exp**2  &
               + two_dp * ( rho_x2x3_n - rho_x1x2_n * rho_x1x3_n )  &
                        * sigma_x2_n * beta_exp * sigma_x3_n * gamma_exp  &
             )  &
           )  &
      * exp( one_fourth_dp * s_cc**2 - ( mu_x1 / sigma_x1 ) * s_cc  &
             + one_half_dp * ( mu_x1**2 / sigma_x1**2 ) )  &
      * gamma( alpha_exp + one_dp ) * Dv_fnc( -(alpha_exp + one_dp), s_cc )


    return

  end function trivar_NLL_mean
  
  !=============================================================================
  function trivar_NLL_mean_const_x1( mu_x1, mu_x2_n, mu_x3_n, &
                                     sigma_x2_n, sigma_x3_n, rho_x2x3_n, &
                                     alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    ! Eq. (62) and Eq. (63) of Larson, V. E. and B. M. Griffin, 2013:  Analytic
    ! upscaling of a local microphysics scheme. Part I: Derivation.
    ! Q. J. Roy. Meteorol. Soc., 139, 670, 46--57,
    ! doi:http://dx.doi.org/10.1002/qj.1967.
    !
    ! Eq. (C29) and Eq. (J26) of Griffin, B. M., 2016:  Improving the
    ! Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
    ! Effects Using a Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S29) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! A new subgrid-scale representation of hydrometeor fields using a
    ! multivariate PDF.  Geosci. Model Dev., 9, 6,
    ! doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
    !
    ! Eq. (S34) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one_half_dp, & ! Constant(s)
        zero_dp

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
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
    real( kind = dp ) ::  &
      trivar_NLL_mean_const_x1


    if ( mu_x1 <= zero_dp ) then

       trivar_NLL_mean_const_x1  &
       = mu_x1**alpha_exp  &
         * exp( mu_x2_n * beta_exp + mu_x3_n * gamma_exp  &
                + one_half_dp * sigma_x2_n**2 * beta_exp**2  &
                + one_half_dp * sigma_x3_n**2 * gamma_exp**2  &
                + rho_x2x3_n * sigma_x2_n * beta_exp * sigma_x3_n * gamma_exp )

    else ! mu_x1 > 0

       trivar_NLL_mean_const_x1 = zero_dp

    endif


    return

  end function trivar_NLL_mean_const_x1

  !=============================================================================
  function trivar_NLL_mean_const_x2( mu_x1, mu_x2, mu_x3_n, &
                                     sigma_x1, sigma_x3_n, rho_x1x3_n, &
                                     alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    ! Eq. (C26), Eq. (C27), Eq. (J27), and Eq. (J28) of Griffin, B. M., 2016:
    ! Improving the Subgrid-Scale Representation of Hydrometeors and
    ! Microphysical Feedback Effects Using a Multivariate PDF.  Doctoral
    ! dissertation, University of Wisconsin -- Milwaukee, Milwaukee, WI,
    ! Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S26) and Eq. (S27) of Griffin, B. M. and V. E. Larson, 2016:
    ! Supplement of A new subgrid-scale representation of hydrometeor fields
    ! using a multivariate PDF.  Geosci. Model Dev., 9, 6,
    ! doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
    !
    ! Eq. (S35) and Eq. (S36) of Griffin, B. M. and V. E. Larson, 2016:
    ! Supplement of Parameterizing microphysical effects on variances and
    ! covariances of moisture and heat content using a multivariate probability
    ! density function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9,
    ! 11, doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp,         & ! Constant(s)
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
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      rho_x1x3_n, & ! Correlation between x1 and ln x3 (ith PDF component)  [-]
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x2                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x3                   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      trivar_NLL_mean_const_x2

    ! Local Variable
    real( kind = dp ) ::  &
      s_cc    !


    s_cc = ( mu_x1 / sigma_x1 ) + rho_x1x3_n * sigma_x3_n * gamma_exp

    trivar_NLL_mean_const_x2  &
    = ( one_dp / sqrt( two_dp*pi_dp ) ) * ( - sigma_x1 )**alpha_exp  &
      * mu_x2**beta_exp &
      * exp( mu_x3_n * gamma_exp &
             + one_half_dp * sigma_x3_n**2 * gamma_exp**2 &
             - one_fourth_dp * s_cc**2 ) &
      * gamma( alpha_exp + one_dp ) * Dv_fnc( -(alpha_exp + one_dp), s_cc )


    return

  end function trivar_NLL_mean_const_x2
  
  !=============================================================================
  function trivar_NLL_mean_const_x1x2( mu_x1, mu_x2, mu_x3_n, sigma_x3_n, &
                                       alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    ! Eq. (C30), Eq. (C31), Eq. (J29), and Eq. (J30) of Griffin, B. M., 2016:
    ! Improving the Subgrid-Scale Representation of Hydrometeors and
    ! Microphysical Feedback Effects Using a Multivariate PDF.  Doctoral
    ! dissertation, University of Wisconsin -- Milwaukee, Milwaukee, WI,
    ! Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S30) and Eq. (S31) of Griffin, B. M. and V. E. Larson, 2016:
    ! Supplement of A new subgrid-scale representation of hydrometeor fields
    ! using a multivariate PDF.  Geosci. Model Dev., 9, 6,
    ! doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
    !
    ! Eq. (S37) and Eq. (S38) of Griffin, B. M. and V. E. Larson, 2016:
    ! Supplement of Parameterizing microphysical effects on variances and
    ! covariances of moisture and heat content using a multivariate probability
    ! density function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9,
    ! 11, doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one_half_dp, & ! Constant(s)
        zero_dp

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x2                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x3                   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      trivar_NLL_mean_const_x1x2


    if ( mu_x1 <= zero_dp ) then

       trivar_NLL_mean_const_x1x2  &
       = mu_x1**alpha_exp * mu_x2**beta_exp &
         * exp( mu_x3_n * gamma_exp &
                + one_half_dp * sigma_x3_n**2 * gamma_exp**2 )

    else ! mu_x1 > 0

       trivar_NLL_mean_const_x1x2 = zero_dp

    endif


    return

  end function trivar_NLL_mean_const_x1x2

  !=============================================================================
  function trivar_NLL_mean_const_x2x3( mu_x1, mu_x2, mu_x3, sigma_x1, &
                                       alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    ! Eq. (C28) and Eq. (J31) of Griffin, B. M., 2016:  Improving the
    ! Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
    ! Effects Using a Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S28) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! A new subgrid-scale representation of hydrometeor fields using a
    ! multivariate PDF.  Geosci. Model Dev., 9, 6,
    ! doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
    !
    ! Eq. (S39) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp,         & ! Constant(s)
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
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3,      & ! Mean of x3 (ith PDF component)                        [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x2                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x3                   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      trivar_NLL_mean_const_x2x3


    trivar_NLL_mean_const_x2x3  &
    = ( one_dp / sqrt( two_dp*pi_dp ) ) * ( - sigma_x1 )**alpha_exp &
      * mu_x2**beta_exp * mu_x3**gamma_exp &
      * exp( - one_fourth_dp * ( mu_x1**2 / sigma_x1**2 ) ) &
      * gamma( alpha_exp + one_dp ) &
      * Dv_fnc( -(alpha_exp + one_dp), ( mu_x1 / sigma_x1 ) )


    return

  end function trivar_NLL_mean_const_x2x3
  
  !=============================================================================
  function trivar_NLL_mean_const_all( mu_x1, mu_x2, mu_x3, &
                                      alpha_exp, beta_exp, gamma_exp )

    ! Description:

    ! References:
    ! Eq. (C32) and Eq. (J32) of Griffin, B. M., 2016:  Improving the
    ! Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
    ! Effects Using a Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S32) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! A new subgrid-scale representation of hydrometeor fields using a
    ! multivariate PDF.  Geosci. Model Dev., 9, 6,
    ! doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
    !
    ! Eq. (S40) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
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
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3,      & ! Mean of x3 (ith PDF component)                        [-]
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x2                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x3                   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      trivar_NLL_mean_const_all


    if ( mu_x1 <= zero_dp ) then

       trivar_NLL_mean_const_all  &
       = mu_x1**alpha_exp * mu_x2**beta_exp * mu_x3**gamma_exp

    else ! mu_x1 > 0

       trivar_NLL_mean_const_all = zero_dp

    endif


    return

  end function trivar_NLL_mean_const_all
  
  !=============================================================================
  function bivar_NL_mean( mu_x1, mu_x2_n, sigma_x1, sigma_x2_n, &
                          rho_x1x2_n, alpha_exp, beta_exp )

    ! Description:

    ! References:
    ! Eq. (33), Eq. (34), Eq. (41), and Eq. (42) of Larson, V. E. and
    ! B. M. Griffin, 2013:  Analytic upscaling of a local microphysics scheme.
    ! Part I: Derivation.  Q. J. Roy. Meteorol. Soc., 139, 670, 46--57,
    ! doi:http://dx.doi.org/10.1002/qj.1967.
    !
    ! Eq. (C4), Eq. (C14), and Eq. (J33) of Griffin, B. M., 2016:  Improving
    ! the Subgrid-Scale Representation of Hydrometeors and Microphysical
    ! Feedback Effects Using a Multivariate PDF.  Doctoral dissertation,
    ! University of Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp.,
    ! URL http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S4) and Eq. (S14) of Griffin, B. M. and V. E. Larson, 2016:
    ! Supplement of A new subgrid-scale representation of hydrometeor fields
    ! using a multivariate PDF.  Geosci. Model Dev., 9, 6,
    ! doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
    !
    ! Eq. (S41) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp,         & ! Constant(s)
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
      mu_x2_n,    & ! Mean of ln x2 (ith PDF component)                     [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x2_n, & ! Standard deviation of ln x2 (ith PDF component)       [-]
      rho_x1x2_n, & ! Correlation between x1 and ln x2 (ith PDF component)  [-]
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp      ! Exponent beta, corresponding to x2                    [-]

    ! Return Variable
    real( kind = dp ) ::  &
      bivar_NL_mean

    ! Local Variable
    real( kind = dp ) ::  &
      s_c    !


    s_c = ( mu_x1 / sigma_x1 ) + rho_x1x2_n * sigma_x2_n * beta_exp

    bivar_NL_mean  &
    = ( one_dp / sqrt( two_dp*pi_dp ) ) * sigma_x1**alpha_exp  &
      * exp( mu_x2_n * beta_exp  &
             + one_half_dp * sigma_x2_n**2 * beta_exp**2  &
             - one_fourth_dp * s_c**2 )  &
      * gamma( alpha_exp + one_dp ) * Dv_fnc( -(alpha_exp + one_dp), -s_c )


    return

  end function bivar_NL_mean

  !=============================================================================
  function bivar_NL_mean_const_x1( mu_x1, mu_x2_n, sigma_x2_n, &
                                   alpha_exp, beta_exp )


    ! Description:

    ! References:
    ! Eq. (35), Eq. (36), Eq. (43), and Eq. (44) of Larson, V. E. and
    ! B. M. Griffin, 2013:  Analytic upscaling of a local microphysics scheme.
    ! Part I: Derivation.  Q. J. Roy. Meteorol. Soc., 139, 670, 46--57,
    ! doi:http://dx.doi.org/10.1002/qj.1967.
    !
    ! Eq. (C9), Eq. (C16), and Eq. (J34) of Griffin, B. M., 2016:  Improving
    ! the Subgrid-Scale Representation of Hydrometeors and Microphysical
    ! Feedback Effects Using a Multivariate PDF.  Doctoral dissertation,
    ! University of Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp.,
    ! URL http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S9) and Eq. (S16) of Griffin, B. M. and V. E. Larson, 2016:
    ! Supplement of A new subgrid-scale representation of hydrometeor fields
    ! using a multivariate PDF.  Geosci. Model Dev., 9, 6,
    ! doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
    !
    ! Eq. (S42) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! Parameterizing microphysical effects on variances and covariances of
    ! moisture and heat content using a multivariate probability density
    ! function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one_half_dp, & ! Constant(s)
        zero_dp

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2_n,    & ! Mean of ln x2 (ith PDF component)                     [-]
      sigma_x2_n, & ! Standard deviation of ln x2 (ith PDF component)       [-]
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp      ! Exponent beta, corresponding to x2                    [-]

    ! Return Variable
    real( kind = dp ) :: &
      bivar_NL_mean_const_x1


    if ( mu_x1 >= zero_dp ) then

       bivar_NL_mean_const_x1  &
       = mu_x1**alpha_exp  &
         * exp( mu_x2_n * beta_exp  &
                + one_half_dp * sigma_x2_n**2 * beta_exp**2 )

    else ! mu_x1 < 0

       bivar_NL_mean_const_x1 = zero_dp

    endif


    return

  end function bivar_NL_mean_const_x1

  !=============================================================================
  function bivar_NL_mean_const_x2( mu_x1, mu_x2, sigma_x1, &
                                   alpha_exp, beta_exp )

    ! Description:

    ! References:
    ! Eq. (C8), Eq. (C15), and Eq. (J35) of Griffin, B. M., 2016:  Improving
    ! the Subgrid-Scale Representation of Hydrometeors and Microphysical
    ! Feedback Effects Using a Multivariate PDF.  Doctoral dissertation,
    ! University of Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp.,
    ! URL http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S8) and Eq. (S15) of Griffin, B. M. and V. E. Larson, 2016:
    ! Supplement of A new subgrid-scale representation of hydrometeor fields
    ! using a multivariate PDF.  Geosci. Model Dev., 9, 6,
    ! doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
    !
    ! Eq. (S43) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
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
      mu_x1,     & ! Mean of x1 (ith PDF component)                         [-]
      mu_x2,     & ! Mean of x2 (ith PDF component)                         [-]
      sigma_x1,  & ! Standard deviation of x1 (ith PDF component)           [-]
      alpha_exp, & ! Exponent alpha, corresponding to x1                    [-]
      beta_exp     ! Exponent beta, corresponding to x2                     [-]

    ! Return Variable
    real( kind = dp ) ::  &
      bivar_NL_mean_const_x2


    bivar_NL_mean_const_x2  &
    = ( one_dp / sqrt( two_dp*pi_dp ) )  &
      * sigma_x1**alpha_exp * mu_x2**beta_exp  &
      * exp( - one_fourth_dp * ( mu_x1**2 / sigma_x1**2 ) )  &
      * gamma( alpha_exp + one_dp )  &
      * Dv_fnc( -(alpha_exp + one_dp), -( mu_x1 / sigma_x1 ) )


    return

  end function bivar_NL_mean_const_x2

  !=============================================================================
  function bivar_NL_mean_const_all( mu_x1, mu_x2, alpha_exp, beta_exp )


    ! Description:

    ! References:
    ! Eq. (C10), Eq. (C17), and Eq. (J36) of Griffin, B. M., 2016:  Improving
    ! the Subgrid-Scale Representation of Hydrometeors and Microphysical
    ! Feedback Effects Using a Multivariate PDF.  Doctoral dissertation,
    ! University of Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp.,
    ! URL http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S10) and Eq. (S17) of Griffin, B. M. and V. E. Larson, 2016:
    ! Supplement of A new subgrid-scale representation of hydrometeor fields
    ! using a multivariate PDF.  Geosci. Model Dev., 9, 6,
    ! doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
    !
    ! Eq. (S44) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
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
      mu_x1,     & ! Mean of x1 (ith PDF component)                         [-]
      mu_x2,     & ! Mean of x2 (ith PDF component)                         [-]
      alpha_exp, & ! Exponent alpha, corresponding to x1                    [-]
      beta_exp     ! Exponent beta, corresponding to x2                     [-]

    ! Return Variable
    real( kind = dp ) :: &
      bivar_NL_mean_const_all


    if ( mu_x1 >= zero_dp ) then

       bivar_NL_mean_const_all = mu_x1**alpha_exp * mu_x2**beta_exp

    else ! mu_x1 < 0

       bivar_NL_mean_const_all = zero_dp

    endif


    return

  end function bivar_NL_mean_const_all

  !=============================================================================
  function bivar_LL_mean( mu_x1_n, mu_x2_n, sigma_x1_n, sigma_x2_n, &
                          rho_x1x2_n, alpha_exp, beta_exp )

    ! Description:
    ! Calculates the analytic solution to the integral:
    !
    ! INT(0:INF) INT(0:INF) x1^alpha x2^beta P_LL(i)(x1,x2) dx2 dx1;
    !
    ! where x1 and x2 are both variables that have a lognormal marginal
    ! distribution in the ith component of the PDF.
    !
    ! The analytic solution to the integral is:
    !
    ! INT(0:INF) INT(0:INF) x1^alpha x2^beta P_LL(i)(x1,x2) dx2 dx1
    ! = exp{ mu_x1_n * alpha + mu_x2_n * beta
    !        + 1/2 * sigma_x1_n^2 * alpha_exp^2
    !        + 1/2 * sigma_x2_n^2 * beta_exp^2
    !        + rho_x1x2_n * sigma_x1_n * alpha * sigma_x2_n * beta }

    ! References:
    ! Eq. (26) of Larson, V. E. and B. M. Griffin, 2013:  Analytic upscaling of
    ! a local microphysics scheme. Part I: Derivation.  Q. J. Roy. Meteorol.
    ! Soc., 139, 670, 46--57, doi:http://dx.doi.org/10.1002/qj.1967.
    !
    ! Eq. (C39) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S39) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! A new subgrid-scale representation of hydrometeor fields using a
    ! multivariate PDF.  Geosci. Model Dev., 9, 6,
    ! doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one_half_dp  ! Constant(s)

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1_n,    & ! Mean of ln x1 (ith PDF component)                     [-]
      mu_x2_n,    & ! Mean of ln x2 (ith PDF component)                     [-]
      sigma_x1_n, & ! Standard deviation of ln x1 (ith PDF component)       [-]
      sigma_x2_n, & ! Standard deviation of ln x2 (ith PDF component)       [-]
      rho_x1x2_n, & ! Correlation between ln x1 & ln x2 (ith PDF component) [-]
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp      ! Exponent beta, corresponding to x2                    [-]

    ! Return Variable
    real( kind = dp ) ::  &
      bivar_LL_mean


    bivar_LL_mean  &
    = exp( mu_x1_n * alpha_exp + mu_x2_n * beta_exp  &
           + one_half_dp * sigma_x1_n**2 * alpha_exp**2  &
           + one_half_dp * sigma_x2_n**2 * beta_exp**2  &
           + rho_x1x2_n * sigma_x1_n * alpha_exp * sigma_x2_n * beta_exp )


    return

  end function bivar_LL_mean

  !=============================================================================
  function bivar_LL_mean_const_x1( mu_x1, mu_x2_n, sigma_x2_n, &
                                   alpha_exp, beta_exp )

    ! Description:

    ! References:
    ! Eq. (C40) and Eq. (C41) of Griffin, B. M., 2016:  Improving the
    ! Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
    ! Effects Using a Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S40) and Eq. (S41) of Griffin, B. M. and V. E. Larson, 2016:
    ! Supplement of A new subgrid-scale representation of hydrometeor fields
    ! using a multivariate PDF.  Geosci. Model Dev., 9, 6,
    ! doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one_half_dp  ! Constant(s)

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2_n,    & ! Mean of ln x2 (ith PDF component)                     [-]
      sigma_x2_n, & ! Standard deviation of ln x2 (ith PDF component)       [-]
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp      ! Exponent beta, corresponding to x2                    [-]

    ! Return Variable
    real( kind = dp ) ::  &
      bivar_LL_mean_const_x1


    bivar_LL_mean_const_x1  &
    = mu_x1**alpha_exp &
      * exp( mu_x2_n * beta_exp + one_half_dp * sigma_x2_n**2 * beta_exp**2 )


    return

  end function bivar_LL_mean_const_x1

  !=============================================================================
  function bivar_LL_mean_const_all( mu_x1, mu_x2, alpha_exp, beta_exp )

    ! Description:

    ! References:
    ! Eq. (C42) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S42) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! A new subgrid-scale representation of hydrometeor fields using a
    ! multivariate PDF.  Geosci. Model Dev., 9, 6,
    ! doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1,     & ! Mean of ln x1 (ith PDF component)                      [-]
      mu_x2,     & ! Mean of ln x2 (ith PDF component)                      [-]
      alpha_exp, & ! Exponent alpha, corresponding to x1                    [-]
      beta_exp     ! Exponent beta, corresponding to x2                     [-]

    ! Return Variable
    real( kind = dp ) ::  &
      bivar_LL_mean_const_all


    bivar_LL_mean_const_all = mu_x1**alpha_exp * mu_x2**beta_exp


    return

  end function bivar_LL_mean_const_all

!===============================================================================

end module PDF_integrals_means
