! $Id$
!===============================================================================
module PDF_integrals_all_mixed_moments

  implicit none

  private

  public :: quadrivar_NNLL_MM, &
            quadrivar_NNLL_MM_const_x1, &
            quadrivar_NNLL_MM_const_x2, &
            quadrivar_NNLL_MM_const_x1_x2, &
            trivar_NNL_MM, &
            trivar_NNL_MM_const_x1, &
            trivar_NNL_MM_const_x2, &
            trivar_NNL_MM_const_x1_x2

  contains

  !=============================================================================
  function quadrivar_NNLL_MM( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                              sigma_x1, sigma_x2, sigma_x3_n, sigma_x4_n, &
                              rho_x1x2, rho_x1x3_n, rho_x1x4_n, &
                              rho_x2x3_n, rho_x2x4_n, rho_x3x4_n, &
                              x1_mean, x2_alpha_x3_beta_x4_gamma_mean, &
                              alpha_exp, beta_exp, gamma_exp, a_exp, b_exp )

    ! Description:

    ! References:
    !  Griffin, B. M. (2011)
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp  ! Constant(s)

    use KK_utilities, only:  &
        factorial, & ! Procedure(s)
        Dv_fnc

    use parabolic, only:  &
        gamma  ! Procedure(s)

    implicit none

    ! Input Variables
    double precision, intent(in) :: &
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

    double precision, intent(in) :: &
      x1_mean,                        & ! Mean of x1 (overall)              [-]
      x2_alpha_x3_beta_x4_gamma_mean    ! Mean of x2^alpha x3^beta x4^gamma [-]
    
    double precision, intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x3                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x4                   [-]

    integer, intent(in) :: &
      a_exp, & ! Order prime of x1 - < x1 >                                 [-]
      b_exp    ! Order prime of x2^alpha x3^beta - < x2^alpha x3^beta >     [-]

    ! Return Variable
    double precision ::  &
      quadrivar_NNLL_MM

    ! Local Variables
    integer :: &
      p, & !
      q, & !
      r    !

    double precision ::  &
      sigma_sum, & !
      s_cc         !

    ! Initialize sigma_sum
    sigma_sum = 0.0

    do p = 0, a_exp/2, 1
       do r = 0, (a_exp-2*p), 1
          do q = 0, b_exp, 1

             s_cc = ( mu_x2 / sigma_x2 )  &
                    + rho_x2x3_n * sigma_x3_n * beta_exp * q  &
                    + rho_x2x4_n * sigma_x4_n * beta_exp * q

             sigma_sum  &
             = sigma_sum  &
             + 1.0 / sqrt( 2.0*pi_dp )  &
               * factorial( a_exp )  &
                 / ( factorial( a_exp - 2*p ) * factorial( p ) )  &
               * factorial( a_exp - 2*p )  &
                 / ( factorial( a_exp - 2*p - r ) * factorial( r ) )  &
               * factorial( b_exp )  &
                 / ( factorial( b_exp - q ) * factorial( q ) )  &
               * ( - x2_alpha_x3_beta_x4_gamma_mean )**(b_exp-q)  &
               * ( - sigma_x2 )**(alpha_exp*q)  &
               * ( 0.5 * ( 1.0 - rho_x1x2**2 ) * sigma_x1**2 )**p  &
               * ( - rho_x1x2 * sigma_x1 )**r  &
               * ( mu_x1 - x1_mean  &
                   - ( mu_x2 / sigma_x2 ) * rho_x1x2 * sigma_x1  &
                   + ( rho_x1x3_n - rho_x1x2 * rho_x2x3_n )  &
                     * sigma_x1 * sigma_x3_n * beta_exp * q  &
                   + ( rho_x1x4_n - rho_x1x2 * rho_x2x4_n )  &
                     * sigma_x1 * sigma_x4_n * gamma_exp * q )**(a_exp-2*p-r)  &
               * exp( mu_x3_n * beta_exp * q + mu_x4_n * gamma_exp * q  &
                      + 0.5 * ( 1.0 - rho_x2x3_n**2 )  &
                            * sigma_x3_n**2 * beta_exp**2 * q**2  &
                      + 0.5 * ( 1.0 - rho_x2x4_n**2 )  &
                            * sigma_x4_n**2 * gamma_exp**2 * q**2  &
                      + ( rho_x3x4_n - rho_x2x3_n * rho_x2x4_n )  &
                            * sigma_x3_n * beta_exp  &
                            * sigma_x4_n * gamma_exp * q**2  &
                    )  &
               * exp( 0.25 * s_cc**2 - ( mu_x2 / sigma_x2 ) * s_cc  &
                      + 0.5 * ( mu_x2**2 / sigma_x2**2 ) )  &
               * gamma( alpha_exp*q + r + 1 )  &
               * Dv_fnc( -(alpha_exp*q + r + 1), s_cc ) ! Known magic numbers

          enddo
       enddo
    enddo 

    quadrivar_NNLL_MM = sigma_sum

    return

  end function quadrivar_NNLL_MM

  !=============================================================================
  function quadrivar_NNLL_MM_const_x1( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                       sigma_x2, sigma_x3_n, sigma_x4_n, &
                                       rho_x2x3_n, rho_x2x4_n, rho_x3x4_n, &
                                       x1_mean, &
                                       x2_alpha_x3_beta_x4_gamma_mean, &
                                       alpha_exp, beta_exp, gamma_exp, &
                                       a_exp, b_exp )

    ! Description:

    ! References:
    !  Griffin, B. M. (2011)
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp  ! Constant(s)

    use KK_utilities, only:  &
        factorial, & ! Procedure(s)
        Dv_fnc

    use parabolic, only:  &
        gamma  ! Procedure(s)

    implicit none

    ! Input Variables
    double precision, intent(in) :: &
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

    double precision, intent(in) :: &
      x1_mean,                        & ! Mean of x1 (overall)              [-]
      x2_alpha_x3_beta_x4_gamma_mean    ! Mean of x2^alpha x3^beta x4^gamma [-]
    
    double precision, intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x3                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x4                   [-]

    integer, intent(in) :: &
      a_exp, & ! Order prime of x1 - < x1 >                                 [-]
      b_exp    ! Order prime of x2^alpha x3^beta - < x2^alpha x3^beta >     [-]

    ! Return Variable
    double precision ::  &
      quadrivar_NNLL_MM_const_x1

    ! Local Variables
    integer :: &
      q    !

    double precision ::  &
      sigma_sum, & !
      s_cc         !

    ! Initialize sigma_sum
    sigma_sum = 0.0

    do q = 0, b_exp, 1

       s_cc = ( mu_x2 / sigma_x2 )  &
              + rho_x2x3_n * sigma_x3_n * beta_exp * q  &
              + rho_x2x4_n * sigma_x4_n * gamma_exp * q

       sigma_sum  &
       = sigma_sum  &
       + 1.0 / sqrt( 2.0*pi_dp )  &
         * ( mu_x1 - x1_mean )**a_exp  &
         * factorial( b_exp )  &
           / ( factorial( b_exp - q ) * factorial( q ) )  &
         * ( - x2_alpha_x3_beta_x4_gamma_mean )**(b_exp-q)  &
         * ( - sigma_x2 )**(alpha_exp*q)  &
         * exp( mu_x3_n * beta_exp * q + mu_x4_n * gamma_exp * q  &
                + 0.5 * ( 1.0 - rho_x2x3_n**2 )  &
                      * sigma_x3_n**2 * beta_exp**2 * q**2  &
                + 0.5 * ( 1.0 - rho_x2x4_n**2 )  &
                      * sigma_x4_n**2 * gamma_exp**2 * q**2  &
                + ( rho_x3x4_n - rho_x2x3_n * rho_x2x4_n )  &
                      * sigma_x3_n * beta_exp  &
                      * sigma_x4_n * gamma_exp * q**2  &
              )  &
         * exp( 0.25 * s_cc**2 - ( mu_x2 / sigma_x2 ) * s_cc  &
                + 0.5 * ( mu_x2**2 / sigma_x2**2 ) )  &
         * gamma( alpha_exp*q + 1 )  &
         * Dv_fnc( -(alpha_exp*q + 1), s_cc ) ! Known magic number

    enddo

    quadrivar_NNLL_MM_const_x1 = sigma_sum

    return

  end function quadrivar_NNLL_MM_const_x1

  !=============================================================================
  function quadrivar_NNLL_MM_const_x2( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                       sigma_x1, sigma_x3_n, sigma_x4_n, &
                                       rho_x1x3_n, rho_x1x4_n, rho_x3x4_n, &
                                       x1_mean, &
                                       x2_alpha_x3_beta_x4_gamma_mean, &
                                       alpha_exp, beta_exp, gamma_exp, &
                                       a_exp, b_exp )

    ! Description:

    ! References:
    !  Griffin, B. M. (2011)
    !-----------------------------------------------------------------------

    use KK_utilities, only:  &
        factorial    ! Procedure(s)

    implicit none

    ! Input Variables
    double precision, intent(in) :: &
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

    double precision, intent(in) :: &
      x1_mean,                        & ! Mean of x1 (overall)              [-]
      x2_alpha_x3_beta_x4_gamma_mean    ! Mean of x2^alpha x3^beta x4^gamma [-]
    
    double precision, intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x3                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x4                   [-]

    integer, intent(in) :: &
      a_exp, & ! Order prime of x1 - < x1 >                                 [-]
      b_exp    ! Order prime of x2^alpha x3^beta - < x2^alpha x3^beta >     [-]

    ! Return Variable
    double precision ::  &
      quadrivar_NNLL_MM_const_x2

    ! Local Variables
    integer :: &
      p, & !
      q    !

    double precision ::  &
      sigma_sum    !

    ! Initialize sigma_sum
    sigma_sum = 0.0

    do p = 0, a_exp/2, 1
       do q = 0, b_exp, 1

          sigma_sum  &
          = sigma_sum  &
          +   factorial( a_exp )  &
              / ( factorial( a_exp - 2*p ) * factorial( p ) )  &
            * factorial( b_exp )  &
              / ( factorial( b_exp - q ) * factorial( q ) )  &
            * ( - x2_alpha_x3_beta_x4_gamma_mean )**(b_exp-q)  &
            * mu_x2**(alpha_exp*q)  &
            * ( 0.5 * sigma_x1**2 )**p  &
            * ( mu_x1 - x1_mean  &
                + rho_x1x3_n * sigma_x1  &
                  * sigma_x3_n * beta_exp * q  &
                + rho_x1x4_n * sigma_x1  &
                  * sigma_x4_n * gamma_exp * q )**(a_exp-2*p)  &
            * exp( mu_x3_n * beta_exp * q + mu_x4_n * gamma_exp * q  &
                   + 0.5 * sigma_x3_n**2 * beta_exp**2 * q**2  &
                   + 0.5 * sigma_x4_n**2 * gamma_exp**2 * q**2  &
                   + rho_x3x4_n * sigma_x3_n * beta_exp  &
                     * sigma_x4_n * gamma_exp * q**2  &
                 )

       enddo
    enddo

    quadrivar_NNLL_MM_const_x2 = sigma_sum

    return

  end function quadrivar_NNLL_MM_const_x2

  !=============================================================================
  function quadrivar_NNLL_MM_const_x1_x2( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                          sigma_x3_n, sigma_x4_n, &
                                          rho_x3x4_n, x1_mean, &
                                          x2_alpha_x3_beta_x4_gamma_mean, &
                                          alpha_exp, beta_exp, gamma_exp, &
                                          a_exp, b_exp )

    ! Description:

    ! References:
    !  Griffin, B. M. (2011)
    !-----------------------------------------------------------------------

    use KK_utilities, only:  &
        factorial    ! Procedure(s)

    implicit none

    ! Input Variables
    double precision, intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      mu_x4_n,    & ! Mean of ln x4 (ith PDF component)                     [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      sigma_x4_n, & ! Standard deviation of ln x4 (ith PDF component)       [-]
      rho_x3x4_n    ! Correlation between ln x3 & ln x4 (ith PDF component) [-]

    double precision, intent(in) :: &
      x1_mean,                        & ! Mean of x1 (overall)              [-]
      x2_alpha_x3_beta_x4_gamma_mean    ! Mean of x2^alpha x3^beta x4^gamma [-]
    
    double precision, intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x3                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x4                   [-]

    integer, intent(in) :: &
      a_exp, & ! Order prime of x1 - < x1 >                                 [-]
      b_exp    ! Order prime of x2^alpha x3^beta - < x2^alpha x3^beta >     [-]

    ! Return Variable
    double precision ::  &
      quadrivar_NNLL_MM_const_x1_x2

    ! Local Variables
    integer :: &
      q    !

    double precision ::  &
      sigma_sum    !

    ! Initialize sigma_sum
    sigma_sum = 0.0

    do q = 0, b_exp, 1

       sigma_sum  &
       = sigma_sum  &
       + ( mu_x1 - x1_mean )**a_exp  &
         * factorial( b_exp )  &
           / ( factorial( b_exp - q ) * factorial( q ) )  &
         * ( - x2_alpha_x3_beta_x4_gamma_mean )**(b_exp-q)  &
         * mu_x2**(alpha_exp*q)  &
         * exp( mu_x3_n * beta_exp * q + mu_x4_n * gamma_exp * q  &
                + 0.5 * sigma_x3_n**2 * beta_exp**2 * q**2  &
                + 0.5 * sigma_x4_n**2 * gamma_exp**2 * q**2  &
                + rho_x3x4_n * sigma_x3_n * beta_exp  &
                      * sigma_x4_n * gamma_exp * q**2  &
              )

    enddo

    quadrivar_NNLL_MM_const_x1_x2 = sigma_sum

    return

  end function quadrivar_NNLL_MM_const_x1_x2

  !=============================================================================
  function trivar_NNL_MM( mu_x1, mu_x2, mu_x3_n, &
                          sigma_x1, sigma_x2, sigma_x3_n, &
                          rho_x1x2, rho_x1x3_n, rho_x2x3_n, &
                          x1_mean, x2_alpha_x3_beta_mean, &
                          alpha_exp, beta_exp, a_exp, b_exp )

    ! Description:

    ! References:
    !  Griffin, B. M. (2011)
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp  ! Constant(s)

    use KK_utilities, only:  &
        factorial, & ! Procedure(s)
        Dv_fnc

    use parabolic, only:  &
        gamma  ! Procedure(s)

    implicit none

    ! Input Variables
    double precision, intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x2,   & ! Standard deviation of x2 (ith PDF component)          [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      rho_x1x2,   & ! Correlation between x1 and x2 (ith PDF component)     [-]
      rho_x1x3_n, & ! Correlation between x1 and ln x3 (ith PDF component)  [-]
      rho_x2x3_n    ! Correlation between x2 and ln x3 (ith PDF component)  [-]

    double precision, intent(in) :: &
      x1_mean,               & ! Mean of x1 (overall)                       [-]
      x2_alpha_x3_beta_mean    ! Mean of x2^alpha x3^beta                   [-]
    
    double precision, intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp      ! Exponent beta, corresponding to x3                    [-]

    integer, intent(in) :: &
      a_exp, & ! Order prime of x1 - < x1 >                                 [-]
      b_exp    ! Order prime of x2^alpha x3^beta - < x2^alpha x3^beta >     [-]

    ! Return Variable
    double precision ::  &
      trivar_NNL_MM

    ! Local Variables
    integer :: &
      p, & !
      q, & !
      r    !

    double precision ::  &
      sigma_sum, & !
      s_c          !

    ! Initialize sigma_sum
    sigma_sum = 0.0

    do p = 0, a_exp/2, 1
       do r = 0, (a_exp-2*p), 1
          do q = 0, b_exp, 1

             s_c = ( mu_x2 / sigma_x2 )  &
                   + rho_x2x3_n * sigma_x3_n * beta_exp * q

             sigma_sum  &
             = sigma_sum  &
             + 1.0 / sqrt( 2.0*pi_dp )  &
               * factorial( a_exp )  &
                 / ( factorial( a_exp - 2*p ) * factorial( p ) )  &
               * factorial( a_exp - 2*p )  &
                 / ( factorial( a_exp - 2*p - r ) * factorial( r ) )  &
               * factorial( b_exp )  &
                 / ( factorial( b_exp - q ) * factorial( q ) )  &
               * ( - x2_alpha_x3_beta_mean )**(b_exp-q)  &
               * sigma_x2**(alpha_exp*q)  &
               * ( 0.5 * ( 1.0 - rho_x1x2**2 ) * sigma_x1**2 )**p  &
               * ( rho_x1x2 * sigma_x1 )**r  &
               * ( mu_x1 - x1_mean  &
                   - ( mu_x2 / sigma_x2 ) * rho_x1x2 * sigma_x1  &
                   + ( rho_x1x3_n - rho_x1x2 * rho_x2x3_n )  &
                     * sigma_x1 * sigma_x3_n * beta_exp * q )**(a_exp-2*p-r)  &
               * exp( mu_x3_n * beta_exp * q  &
                      + 0.5 * sigma_x3_n**2 * beta_exp**2 * q**2  &
                      - 0.25 * s_c**2 )  &
               * gamma( alpha_exp*q + r + 1 )  &
               * Dv_fnc( -( alpha_exp*q + r + 1 ), -s_c ) ! Known magic number

          enddo
       enddo
    enddo 

    trivar_NNL_MM = sigma_sum

    return

  end function trivar_NNL_MM

  !=============================================================================
  function trivar_NNL_MM_const_x1( mu_x1, mu_x2, mu_x3_n, &
                                   sigma_x2, sigma_x3_n, rho_x2x3_n, &
                                   x1_mean, x2_alpha_x3_beta_mean, &
                                   alpha_exp, beta_exp, a_exp, b_exp )

    ! Description:

    ! References:
    !  Griffin, B. M. (2011)
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        pi_dp  ! Constant(s)

    use KK_utilities, only:  &
        factorial, & ! Procedure(s)
        Dv_fnc

    use parabolic, only:  &
        gamma  ! Procedure(s)

    implicit none

    ! Input Variables
    double precision, intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      sigma_x2,   & ! Standard deviation of x2 (ith PDF component)          [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      rho_x2x3_n    ! Correlation between x2 and ln x3 (ith PDF component)  [-]

    double precision, intent(in) :: &
      x1_mean,               & ! Mean of x1 (overall)                       [-]
      x2_alpha_x3_beta_mean    ! Mean of x2^alpha x3^beta                   [-]
    
    double precision, intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp      ! Exponent beta, corresponding to x3                    [-]

    integer, intent(in) :: &
      a_exp, & ! Order prime of x1 - < x1 >                                 [-]
      b_exp    ! Order prime of x2^alpha x3^beta - < x2^alpha x3^beta >     [-]

    ! Return Variable
    double precision ::  &
      trivar_NNL_MM_const_x1

    ! Local Variables
    integer :: &
      q    !

    double precision ::  &
      sigma_sum, & !
      s_c          !

    ! Initialize sigma_sum
    sigma_sum = 0.0

    do q = 0, b_exp, 1

       s_c = ( mu_x2 / sigma_x2 )  &
             + rho_x2x3_n * sigma_x3_n * beta_exp * q

       sigma_sum  &
       = sigma_sum  &
       + 1.0 / sqrt( 2.0*pi_dp )  &
         * ( mu_x1 - x1_mean )**a_exp  &
         * factorial( b_exp )  &
           / ( factorial( b_exp - q ) * factorial( q ) )  &
         * ( - x2_alpha_x3_beta_mean )**(b_exp-q)  &
         * sigma_x2**(alpha_exp*q)  &
         * exp( mu_x3_n * beta_exp * q  &
                + 0.5 * sigma_x3_n**2 * beta_exp**2 * q**2  &
                - 0.25 * s_c**2 )  &
         * gamma( alpha_exp*q + 1 )  &
         * Dv_fnc( -(alpha_exp*q + 1), -s_c ) ! Known magic number

    enddo

    trivar_NNL_MM_const_x1 = sigma_sum

    return

  end function trivar_NNL_MM_const_x1

  !=============================================================================
  function trivar_NNL_MM_const_x2( mu_x1, mu_x2, mu_x3_n, &
                                   sigma_x1, sigma_x3_n, rho_x1x3_n, &
                                   x1_mean, x2_alpha_x3_beta_mean, &
                                   alpha_exp, beta_exp, a_exp, b_exp )

    ! Description:

    ! References:
    !  Griffin, B. M. (2011)
    !-----------------------------------------------------------------------

    use KK_utilities, only:  &
        factorial    ! Procedure(s)

    implicit none

    ! Input Variables
    double precision, intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      rho_x1x3_n    ! Correlation between x1 and ln x3 (ith PDF component)  [-]

    double precision, intent(in) :: &
      x1_mean,               & ! Mean of x1 (overall)                       [-]
      x2_alpha_x3_beta_mean    ! Mean of x2^alpha x3^beta                   [-]
    
    double precision, intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp      ! Exponent beta, corresponding to x3                    [-]

    integer, intent(in) :: &
      a_exp, & ! Order prime of x1 - < x1 >                                 [-]
      b_exp    ! Order prime of x2^alpha x3^beta - < x2^alpha x3^beta >     [-]

    ! Return Variable
    double precision ::  &
      trivar_NNL_MM_const_x2

    ! Local Variables
    integer :: &
      p, & !
      q    !

    double precision ::  &
      sigma_sum    !

    ! Initialize sigma_sum
    sigma_sum = 0.0

    do p = 0, a_exp/2, 1
       do q = 0, b_exp, 1

          sigma_sum  &
          = sigma_sum  &
          +   factorial( a_exp )  &
              / ( factorial( a_exp - 2*p ) * factorial( p ) )  &
            * factorial( b_exp )  &
              / ( factorial( b_exp - q ) * factorial( q ) )  &
            * ( - x2_alpha_x3_beta_mean )**(b_exp-q)  &
            * mu_x2**(alpha_exp*q)  &
            * ( 0.5 * sigma_x1**2 )**p  &
            * ( mu_x1 - x1_mean  &
                + rho_x1x3_n * sigma_x1  &
                  * sigma_x3_n * beta_exp * q )**(a_exp-2*p)  &
            * exp( mu_x3_n * beta_exp * q  &
                   + 0.5 * sigma_x3_n**2 * beta_exp**2 * q**2 )

       enddo
    enddo 

    trivar_NNL_MM_const_x2 = sigma_sum

    return

  end function trivar_NNL_MM_const_x2

  !=============================================================================
  function trivar_NNL_MM_const_x1_x2( mu_x1, mu_x2, mu_x3_n, sigma_x3_n, &
                                      x1_mean, x2_alpha_x3_beta_mean, &
                                      alpha_exp, beta_exp, a_exp, b_exp )

    ! Description:

    ! References:
    !  Griffin, B. M. (2011)
    !-----------------------------------------------------------------------

    use KK_utilities, only:  &
        factorial    ! Procedure(s)

    implicit none

    ! Input Variables
    double precision, intent(in) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      sigma_x3_n    ! Standard deviation of ln x3 (ith PDF component)       [-]

    double precision, intent(in) :: &
      x1_mean,               & ! Mean of x1 (overall)                       [-]
      x2_alpha_x3_beta_mean    ! Mean of x2^alpha x3^beta                   [-]
    
    double precision, intent(in) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp      ! Exponent beta, corresponding to x3                    [-]

    integer, intent(in) :: &
      a_exp, & ! Order prime of x1 - < x1 >                                 [-]
      b_exp    ! Order prime of x2^alpha x3^beta - < x2^alpha x3^beta >     [-]

    ! Return Variable
    double precision ::  &
      trivar_NNL_MM_const_x1_x2

    ! Local Variables
    integer :: &
      q    !

    double precision ::  &
      sigma_sum    !

    ! Initialize sigma_sum
    sigma_sum = 0.0

    do q = 0, b_exp, 1

       sigma_sum  &
       = sigma_sum  &
       + ( mu_x1 - x1_mean )**a_exp  &
         * factorial( b_exp )  &
           / ( factorial( b_exp - q ) * factorial( q ) )  &
         * ( - x2_alpha_x3_beta_mean )**(b_exp-q)  &
         * mu_x2**(alpha_exp*q)  &
         * exp( mu_x3_n * beta_exp * q  &
                + 0.5 * sigma_x3_n**2 * beta_exp**2 * q**2 )

    enddo

    trivar_NNL_MM_const_x1_x2 = sigma_sum

    return

  end function trivar_NNL_MM_const_x1_x2

!===============================================================================

END MODULE PDF_integrals_all_mixed_moments
