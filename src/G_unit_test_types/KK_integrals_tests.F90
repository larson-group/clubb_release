! $Id$
!===============================================================================
module KK_integrals_tests

  implicit none

  private

  public :: KK_integrals_tests_driver

  private :: KK_covar_tests, &
             KK_mean_tests, &
             quadrivar_NNLL_covar_tests, &
             trivar_NNL_covar_tests, &
             trivar_NLL_mean_tests, &
             bivar_NL_mean_tests, &
             bivar_LL_mean_tests, &
             percent_difference

  contains

  !=============================================================================
  function KK_integrals_tests_driver()

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        fstdout  ! Constant(s)

    use clubb_precision, only: &
        dp ! double precision

    use parabolic, only:  & 
      l_high_accuracy_parab_cyl_fnc ! Variable

    implicit none

    ! Output Vars
    integer :: KK_integrals_tests_driver ! Returns the exit code of the test

    ! Local Variables
    integer :: &
      opt  ! Case option involving differing component mean values of the
           ! saturation variable

    integer :: &
      total_mismatches   ! Total number of mismatches

    real( kind = dp ) :: &
      tol      ! Acceptable percent difference between the results    [-]

  !-----------------------------------------------------------------------
    !----- Begin Code -----

    ! Declare the tolerance value, which is the maximum acceptable percent
    ! difference between the results.
    tol = 1.5e-12_dp

    ! Initialize total number of mismatches.
    total_mismatches = 0

    ! Use high accuracy for tests
    l_high_accuracy_parab_cyl_fnc = .true.

    write(fstdout,'(A)') "Performing KK integrals tests"
    write(fstdout,'(A)') " "
    write(fstdout,'(A,2x,ES12.5)') "Percent Difference Tolerance (maximum " &
                                   // "acceptable percent difference " &
                                   // "between the CLUBB results and " &
                                   // "the MATLAB results using the same " &
                                   // "equations):", tol
    write(fstdout,'(A)') " "

    do opt = 1,3,1

       write(fstdout,'(A,2x,I1)') "Saturation option (x2_opt for covariance " &
                                  // "integrals and x1_opt for mean " &
                                  // "integrals):", opt
       write(fstdout,'(A)') " "

       call KK_covar_tests( opt, tol, total_mismatches )

       write(fstdout,'(A)') "--------------------------------------------------"
       write(fstdout,'(A)') " "

       call KK_mean_tests( opt, tol, total_mismatches )

       write(fstdout,'(A)') "=================================================="
       write(fstdout,'(A)') " "

    enddo

    if ( total_mismatches == 0 ) then

       write(fstdout,'(A)') "Success!"
       KK_integrals_tests_driver = 0 ! Exit Code = 0, Success!

    else  ! total_mismatches > 0

       write(fstdout,'(A,1x,I2,1x,A)') "There were", total_mismatches, &
                                       "total mismatches found."
       KK_integrals_tests_driver = 1 ! Exit Code = 1, Fail

    endif


    return

  end function KK_integrals_tests_driver

  !=============================================================================
  subroutine KK_covar_tests( x2_opt, tol, total_mismatches )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one_dp, &  ! Constant(s)
        two_dp

    use pdf_utilities, only: &
        mean_L2N_dp,   & ! Procedure(s)
        stdev_L2N_dp,  &
        corr_NL2NN_dp, &
        corr_LL2NN_dp

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      x2_opt   ! Case option involving differing component mean values of x2.

    real( kind = dp ), intent(in) :: &
      tol      ! Acceptable percent difference between the results          [-]

    ! Input/Output Variable
    integer, intent(inout) :: &
      total_mismatches   ! Total number of mismatches

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
      rho_x1x2,   & ! Correlation between x1 and x2 (ith PDF component)     [-]
      rho_x1x3,   & ! Correlation between x1 and x3 (ith PDF component)     [-]
      rho_x1x4,   & ! Correlation between x1 and x4 (ith PDF component)     [-]
      rho_x2x3,   & ! Correlation between x2 and x3 (ith PDF component)     [-]
      rho_x2x4,   & ! Correlation between x2 and x4 (ith PDF component)     [-]
      rho_x3x4,   & ! Correlation between x3 and x4 (ith PDF component)     [-]
      rho_x1x3_n, & ! Correlation between x1 and ln x3 (ith PDF component)  [-]
      rho_x1x4_n, & ! Correlation between x1 and ln x4 (ith PDF component)  [-]
      rho_x2x3_n, & ! Correlation between x2 and ln x3 (ith PDF component)  [-]
      rho_x2x4_n, & ! Correlation between x2 and ln x4 (ith PDF component)  [-]
      rho_x3x4_n    ! Correlation between ln x3 & ln x4 (ith PDF component) [-]

    real( kind = dp ) :: &
      x1_mean,                        & ! Mean of x1 (overall)              [-]
      x2_alpha_x3_beta_x4_gamma_mean, & ! Mean of x2^alpha x3^beta x4^gamma [-]
      x2_alpha_x3_beta_mean             ! Mean of x2^alpha x3^beta
    
    real( kind = dp ) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x2                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x3                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x4                   [-]

    real( kind = dp ) :: &
      quadrivar_NNLL_covar_int,      & ! Quadrivar. int. (all fields vary)  [-]
      quadrivar_NNLL_covar_int_c1,   & ! Quadrivar. int. (constant x1)      [-]
      quadrivar_NNLL_covar_int_c2,   & ! Quadrivar. int. (constant x2)      [-]
      quadrivar_NNLL_covar_int_c1c2    ! Quadrivar. int. (constant x1, x2)  [-]

    real( kind = dp ) :: &
      trivar_NNL_covar_int,      & ! Trivariate integral (all fields vary)  [-]
      trivar_NNL_covar_int_c1,   & ! Trivariate integral (constant x1)      [-]
      trivar_NNL_covar_int_c2,   & ! Trivariate integral (constant x2)      [-]
      trivar_NNL_covar_int_c1c2    ! Trivariate integral (constant x1, x2)  [-]

    real( kind = dp ) :: &
      quadrivar_NNLL_cmprsn_res,      & ! Quadrivar. int. comparison result [-]
      quadrivar_NNLL_c1_cmprsn_res,   & ! Quadrivar. (const. x1) comp. res. [-]
      quadrivar_NNLL_c2_cmprsn_res,   & ! Quadrivar. (const. x2) comp. res. [-]
      quadrivar_NNLL_c1c2_cmprsn_res    ! Quadrivar. (const. x1, x2) comp.  [-]

    real( kind = dp ) :: &
      trivar_NNL_cmprsn_res,      & ! Trivariate integral comparison result [-]
      trivar_NNL_c1_cmprsn_res,   & ! Trivariate (const. x1) comp. result   [-]
      trivar_NNL_c2_cmprsn_res,   & ! Trivariate (const. x2) comp. result   [-]
      trivar_NNL_c1c2_cmprsn_res    ! Trivariate (const. x1, x2) comp. res. [-]

    integer :: &
      count_incr    ! Increment total number of mismatches


    !!! Quadrivariate NNLL Covariance Tests
    ! The actual means, variances, and correlations of x1, x2, x3, and x4.
    mu_x1 = -one_dp
    !mu_x2 = -1.0e-3_dp
    mu_x3 = 1.0e-4_dp
    mu_x4 = 100.0_dp
    sigma_x1 = 0.125_dp
    !sigma_x2 = 1.25e-4_dp
    sigma_x3 = 2.0e-5_dp
    sigma_x4 = 25.0_dp
    rho_x1x2 = 0.10_dp
    rho_x1x3 = 0.15_dp
    rho_x1x4 = -0.05_dp
    rho_x2x3 = 0.2_dp
    rho_x2x4 = 0.25_dp
    rho_x3x4 = 0.3_dp
    if ( x2_opt == 1 ) then
       ! The individual marginal PDF of x2 is basically entirely less than 0.
       mu_x2 = -1.0e-3_dp
       sigma_x2 = 1.25e-4_dp
    elseif ( x2_opt == 2 ) then
       ! The individual marginal PDF of x2 is mostly less than 0,
       ! but there is a significant portion that is greater than 0.
       mu_x2 = -1.0e-4_dp
       sigma_x2 = 1.25e-4_dp
    elseif ( x2_opt == 3 ) then
       ! The individual marginal PDF of x2 is slightly more heavily greater
       ! than 0 than it is less than 0.
       mu_x2 = 1.0e-7_dp
       sigma_x2 = 1.25e-4_dp
    else ! default
       ! Same as x2_opt = 1.
       mu_x2 = -1.0e-3_dp
       sigma_x2 = 1.25e-4_dp
    endif

    ! Normalize the means, variances, and correlations involving variables that
    ! have individual marginals that are distributed lognormally.
    mu_x3_n = mean_L2N_dp( mu_x3, (sigma_x3/mu_x3)**2 )
    mu_x4_n = mean_L2N_dp( mu_x4, (sigma_x4/mu_x4)**2 )

    sigma_x3_n = stdev_L2N_dp( (sigma_x3/mu_x3)**2 )
    sigma_x4_n = stdev_L2N_dp( (sigma_x4/mu_x4)**2 )

    rho_x1x3_n = corr_NL2NN_dp( rho_x1x3, sigma_x3_n, (sigma_x3/mu_x3)**2 )
    rho_x1x4_n = corr_NL2NN_dp( rho_x1x4, sigma_x4_n, (sigma_x4/mu_x4)**2 )
    rho_x2x3_n = corr_NL2NN_dp( rho_x2x3, sigma_x3_n, (sigma_x3/mu_x3)**2 )
    rho_x2x4_n = corr_NL2NN_dp( rho_x2x4, sigma_x4_n, (sigma_x4/mu_x4)**2 )
    rho_x3x4_n = corr_LL2NN_dp( rho_x3x4, sigma_x3_n, sigma_x4_n, &
                                (sigma_x3/mu_x3)**2, (sigma_x4/mu_x4)**2  )

    ! Overall mean of x1.
    x1_mean = -0.5_dp
    ! Overall mean of (x2 H(-x2))^alpha x3^beta x4^gamma.
    x2_alpha_x3_beta_x4_gamma_mean = -2.0e-4_dp

    ! Exponent corresponding to x2.
    alpha_exp = one_dp
    ! Exponent corresponding to x3.
    beta_exp  = one_dp/3.0_dp
    ! Exponent corresponding to x4.
    gamma_exp = two_dp/3.0_dp


    ! Obtain the results for the quadrivariate NNLL covariance integral,
    ! including special cases.
    call quadrivar_NNLL_covar_tests( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                     sigma_x1, sigma_x2, sigma_x3_n, &
                                     sigma_x4_n, rho_x1x2, rho_x1x3_n, &
                                     rho_x1x4_n, rho_x2x3_n, rho_x2x4_n, &
                                     rho_x3x4_n, x1_mean, &
                                     x2_alpha_x3_beta_x4_gamma_mean, &
                                     alpha_exp, beta_exp, gamma_exp, &
                                     quadrivar_NNLL_covar_int, &
                                     quadrivar_NNLL_covar_int_c1, &
                                     quadrivar_NNLL_covar_int_c2, &
                                     quadrivar_NNLL_covar_int_c1c2 )


    ! Now that the results from the CLUBB code have been obtained, compare them
    ! to the MATLAB results from the same equations.
    if ( x2_opt == 1 ) then
       quadrivar_NNLL_cmprsn_res      = 3.940566852245361e-04_dp
       quadrivar_NNLL_c1_cmprsn_res   = 3.927115990964175e-04_dp
       quadrivar_NNLL_c2_cmprsn_res   = 3.959157231937519e-04_dp
       quadrivar_NNLL_c1c2_cmprsn_res = 3.961224408774508e-04_dp
    elseif ( x2_opt == 2 ) then
       quadrivar_NNLL_cmprsn_res      = -4.441972023148829e-05_dp
       quadrivar_NNLL_c1_cmprsn_res   = -4.559373570842311e-05_dp
       quadrivar_NNLL_c2_cmprsn_res   = -5.040842768062482e-05_dp
       quadrivar_NNLL_c1c2_cmprsn_res = -5.038775591225493e-05_dp
    elseif ( x2_opt == 3 ) then
       quadrivar_NNLL_cmprsn_res      = -7.622010070476262e-05_dp
       quadrivar_NNLL_c1_cmprsn_res   = -7.695119321547280e-05_dp
       quadrivar_NNLL_c2_cmprsn_res   = -1.0e-04_dp
       quadrivar_NNLL_c1c2_cmprsn_res = -1.0e-04_dp
    else ! default
       quadrivar_NNLL_cmprsn_res      = 3.940566852245361e-04_dp
       quadrivar_NNLL_c1_cmprsn_res   = 3.927115990964175e-04_dp
       quadrivar_NNLL_c2_cmprsn_res   = 3.959157231937519e-04_dp
       quadrivar_NNLL_c1c2_cmprsn_res = 3.961224408774508e-04_dp
    endif
       
    call percent_difference( 'Quadrivariate NNLL', &
                             quadrivar_NNLL_covar_int, &
                             quadrivar_NNLL_cmprsn_res, tol, &
                             count_incr )

    total_mismatches = total_mismatches + count_incr

    call percent_difference( 'Quadrivariate NNLL (const. x1)', &
                             quadrivar_NNLL_covar_int_c1, &
                             quadrivar_NNLL_c1_cmprsn_res, tol, &
                             count_incr )

    total_mismatches = total_mismatches + count_incr

    call percent_difference( 'Quadrivariate NNLL (const. x2)', &
                             quadrivar_NNLL_covar_int_c2, &
                             quadrivar_NNLL_c2_cmprsn_res, tol, &
                             count_incr )

    total_mismatches = total_mismatches + count_incr

    call percent_difference( 'Quadrivariate NNLL (const. x1 and x2)', &
                             quadrivar_NNLL_covar_int_c1c2, &
                             quadrivar_NNLL_c1c2_cmprsn_res, tol, &
                             count_incr )

    total_mismatches = total_mismatches + count_incr


    !!! Trivariate NNL Covariance Tests
    ! The actual means, variances, and correlations of x1, x2, and x3.
    mu_x1 = -one_dp
    !mu_x2 = 1.0e-3_dp
    mu_x3 = 1.0e-4_dp
    sigma_x1 = 0.125_dp
    !sigma_x2 = 1.25e-4_dp
    sigma_x3 = 2.0e-5_dp
    rho_x1x2 = 0.25_dp
    rho_x1x3 = 0.10_dp
    rho_x2x3 = -0.05_dp
    if ( x2_opt == 1 ) then
       ! The individual marginal PDF of x2 is basically entirely greater than 0.
       mu_x2 = 1.0e-3_dp
       sigma_x2 = 1.25e-4_dp
    elseif ( x2_opt == 2 ) then
       ! The individual marginal PDF of x2 is mostly greater than 0,
       ! but there is a significant portion that is less than 0.
       mu_x2 = 1.0e-4_dp
       sigma_x2 = 1.25e-4_dp
    elseif ( x2_opt == 3 ) then
       ! The individual marginal PDF of x2 is slightly more heavily less
       ! than 0 than it is greater than 0.
       mu_x2 = -1.0e-7_dp
       sigma_x2 = 1.25e-4_dp
    else ! default
       ! Same as x2_opt = 1.
       mu_x2 = 1.0e-3_dp
       sigma_x2 = 1.25e-4_dp
    endif

    ! Normalize the means, variances, and correlations involving variables that
    ! have individual marginals that are distributed lognormally.
    mu_x3_n = mean_L2N_dp( mu_x3, (sigma_x3/mu_x3)**2 )

    sigma_x3_n = stdev_L2N_dp( (sigma_x3/mu_x3)**2 )

    rho_x1x3_n = corr_NL2NN_dp( rho_x1x3, sigma_x3_n, (sigma_x3/mu_x3)**2 )
    rho_x2x3_n = corr_NL2NN_dp( rho_x2x3, sigma_x3_n, (sigma_x3/mu_x3)**2 )

    ! Overall mean of x1.
    x1_mean = -0.5_dp
    ! Overall mean of (x2 H(x2))^alpha x3^beta.
    x2_alpha_x3_beta_mean = 1.0e-8_dp

    ! Exponent corresponding to x2.
    alpha_exp = 1.15_dp
    ! Exponent corresponding to x3.
    beta_exp  = 1.15_dp


    ! Obtain the results for the trivariate NNL covariance integral,
    ! including special cases.
    call trivar_NNL_covar_tests( mu_x1, mu_x2, mu_x3_n, &
                                 sigma_x1, sigma_x2, sigma_x3_n, &
                                 rho_x1x2, rho_x1x3_n, rho_x2x3_n, &
                                 x1_mean, x2_alpha_x3_beta_mean, &
                                 alpha_exp, beta_exp, &
                                 trivar_NNL_covar_int, &
                                 trivar_NNL_covar_int_c1, &
                                 trivar_NNL_covar_int_c2, &
                                 trivar_NNL_covar_int_c1c2 )


    ! Now that the results from the CLUBB code have been obtained, compare them
    ! to the MATLAB results from the same equations.
    if ( x2_opt == 1 ) then
       trivar_NNL_cmprsn_res      = 5.957899876478559e-10_dp
       trivar_NNL_c1_cmprsn_res   = 5.299646183560650e-10_dp
       trivar_NNL_c2_cmprsn_res   = 5.543555152546526e-10_dp
       trivar_NNL_c1c2_cmprsn_res = 5.286452253001284e-10_dp
    elseif ( x2_opt == 2 ) then
       trivar_NNL_cmprsn_res      = 4.629614112194527e-09_dp
       trivar_NNL_c1_cmprsn_res   = 4.604609365066509e-09_dp
       trivar_NNL_c2_cmprsn_res   = 4.685272472815394e-09_dp
       trivar_NNL_c1c2_cmprsn_res = 4.683452323676534e-09_dp
    elseif ( x2_opt == 3 ) then
       trivar_NNL_cmprsn_res      = 4.851589269710427e-09_dp
       trivar_NNL_c1_cmprsn_res   = 4.837259568295437e-09_dp
       trivar_NNL_c2_cmprsn_res   = 5.0e-09_dp
       trivar_NNL_c1c2_cmprsn_res = 5.0e-09_dp
    else ! default
       trivar_NNL_cmprsn_res      = 5.957899876478559e-10_dp
       trivar_NNL_c1_cmprsn_res   = 5.299646183560650e-10_dp
       trivar_NNL_c2_cmprsn_res   = 5.543555152546526e-10_dp
       trivar_NNL_c1c2_cmprsn_res = 5.286452253001284e-10_dp
    endif

    call percent_difference( 'Trivariate NNL', &
                             trivar_NNL_covar_int, &
                             trivar_NNL_cmprsn_res, tol, &
                             count_incr )

    total_mismatches = total_mismatches + count_incr

    call percent_difference( 'Trivariate NNL (const. x1)', &
                             trivar_NNL_covar_int_c1, &
                             trivar_NNL_c1_cmprsn_res, tol, &
                             count_incr )

    total_mismatches = total_mismatches + count_incr

    call percent_difference( 'Trivariate NNL (const. x2)', &
                             trivar_NNL_covar_int_c2, &
                             trivar_NNL_c2_cmprsn_res, tol, &
                             count_incr )

    total_mismatches = total_mismatches + count_incr

    call percent_difference( 'Trivariate NNL (const. x1 and x2)', &
                             trivar_NNL_covar_int_c1c2, &
                             trivar_NNL_c1c2_cmprsn_res, tol, &
                             count_incr )

    total_mismatches = total_mismatches + count_incr


    return

  end subroutine KK_covar_tests

  !=============================================================================
  subroutine KK_mean_tests( x1_opt, tol, total_mismatches )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two_dp,  & ! Constant(s)
        one_dp,  &
        zero_dp

    use pdf_utilities, only: &
        mean_L2N_dp,   & ! Procedure(s)
        stdev_L2N_dp,  &
        corr_NL2NN_dp, &
        corr_LL2NN_dp

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      x1_opt   ! Case option involving differing component mean values of x1.

    real( kind = dp ), intent(in) :: &
      tol      ! Acceptable percent difference between the results          [-]

    ! Input/Output Variable
    integer, intent(inout) :: &
      total_mismatches   ! Total number of mismatches

    ! Local Variables
    real( kind = dp ) :: &
      mu_x1,      & ! Mean of x1 (ith PDF component)                        [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      mu_x3,      & ! Mean of x3 (ith PDF component)                        [-]
      mu_x1_n,    & ! Mean of ln x1 (ith PDF component)                     [-]
      mu_x2_n,    & ! Mean of ln x2 (ith PDF component)                     [-]
      mu_x3_n,    & ! Mean of ln x3 (ith PDF component)                     [-]
      sigma_x1,   & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x2,   & ! Standard deviation of x2 (ith PDF component)          [-]
      sigma_x3,   & ! Standard deviation of x3 (ith PDF component)          [-]
      sigma_x1_n, & ! Standard deviation of ln x1 (ith PDF component)       [-]
      sigma_x2_n, & ! Standard deviation of ln x2 (ith PDF component)       [-]
      sigma_x3_n, & ! Standard deviation of ln x3 (ith PDF component)       [-]
      rho_x1x2,   & ! Correlation between x1 and x2 (ith PDF component)     [-]
      rho_x1x3,   & ! Correlation between x1 and x3 (ith PDF component)     [-]
      rho_x2x3,   & ! Correlation between x2 and x3 (ith PDF component)     [-]
      rho_x1x2_n, & ! Correlation between x1 and ln x2 (ith PDF component)  [-]
      rho_x1x3_n, & ! Correlation between x1 and ln x3 (ith PDF component)  [-]
      rho_x2x3_n    ! Correlation between ln x2 & ln x3 (ith PDF component) [-]

    ! For the Bivariate LL mean integral, rho_x1x2_n is the correlation between
    ! ln x1 and ln x2 (ith PDF component).  Additionally, mu_x1_n and sigma_x1_n
    ! are used only for the Bivariate LL mean integral. 
    
    real( kind = dp ) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp,   & ! Exponent beta, corresponding to x2                    [-]
      gamma_exp     ! Exponent gamma, corresponding to x3                   [-]

    real( kind = dp ) :: &
      trivar_NLL_mean_int,    & ! Trivariate integral (all fields vary)     [-]
      trivar_NLL_mean_int_c1    ! Trivariate integral (constant x1)         [-]

    real( kind = dp ) :: &
      bivar_NL_mean_int,    & ! Bivariate NL integral (all fields vary)     [-]
      bivar_NL_mean_int_c1    ! Bivariate NL integral (constant x1)         [-]

    real( kind = dp ) :: &
      bivar_LL_mean_int       ! Bivariate LL integral (all fields vary)     [-]

    real( kind = dp ) :: &
      trivar_NLL_cmprsn_res,    & ! Trivariate integral comparison result   [-]
      trivar_NLL_c1_cmprsn_res    ! Trivariate (const. x1) comp. result     [-]

    real( kind = dp ) :: &
      bivar_NL_cmprsn_res,    & ! Bivariate NL integral comparison result   [-]
      bivar_NL_c1_cmprsn_res    ! Bivariate NL (const. x1) comp. result     [-]

    real( kind = dp ) :: &
      bivar_LL_cmprsn_res       ! Bivariate LL integral comparison result   [-]

    integer :: &
      count_incr    ! Increment total number of mismatches


    !!! Trivariate NLL Mean Tests
    ! The actual means, variances, and correlations of x1, x2, and x3.
    !mu_x1 = -1.0e-3_dp
    mu_x2 = 1.0e-4_dp
    mu_x3 = 100.0_dp
    !sigma_x1 = 1.25e-4_dp
    sigma_x2 = 2.0e-5_dp
    sigma_x3 = 25.0_dp
    rho_x1x2 = 0.2_dp
    rho_x1x3 = 0.25_dp
    rho_x2x3 = 0.3_dp
    if ( x1_opt == 1 ) then
       ! The individual marginal PDF of x1 is basically entirely less than 0.
       mu_x1 = -1.0e-3_dp
       sigma_x1 = 1.25e-4_dp
    elseif ( x1_opt == 2 ) then
       ! The individual marginal PDF of x1 is mostly less than 0,
       ! but there is a significant portion that is greater than 0.
       mu_x1 = -1.0e-4_dp
       sigma_x1 = 1.25e-4_dp
    elseif ( x1_opt == 3 ) then
       ! The individual marginal PDF of x1 is slightly more heavily greater
       ! than 0 than it is less than 0.
       mu_x1 = 1.0e-7_dp
       sigma_x1 = 1.25e-4_dp
    else ! default
       ! Same as x1_opt = 1.
       mu_x1 = -1.0e-3_dp
       sigma_x1 = 1.25e-4_dp
    endif

    ! Normalize the means, variances, and correlations involving variables that
    ! have individual marginals that are distributed lognormally.
    mu_x2_n = mean_L2N_dp( mu_x2, (sigma_x2/mu_x2)**2 )
    mu_x3_n = mean_L2N_dp( mu_x3, (sigma_x3/mu_x3)**2 )

    sigma_x2_n = stdev_L2N_dp( (sigma_x2/mu_x2)**2 )
    sigma_x3_n = stdev_L2N_dp( (sigma_x3/mu_x3)**2 )

    rho_x1x2_n = corr_NL2NN_dp( rho_x1x2, sigma_x2_n, (sigma_x2/mu_x2)**2 )
    rho_x1x3_n = corr_NL2NN_dp( rho_x1x3, sigma_x3_n, (sigma_x3/mu_x3)**2 )
    rho_x2x3_n = corr_LL2NN_dp( rho_x2x3, sigma_x2_n, sigma_x3_n, &
                                (sigma_x2/mu_x2)**2, (sigma_x3/mu_x3)**2 )

    ! Exponent corresponding to x1.
    alpha_exp = one_dp
    ! Exponent corresponding to x2.
    beta_exp  = one_dp/3.0_dp
    ! Exponent corresponding to x3.
    gamma_exp = two_dp/3.0_dp


    ! Obtain the results for the trivariate NLL mean integral,
    ! including special cases.
    call trivar_NLL_mean_tests( mu_x1, mu_x2_n, mu_x3_n, &
                                sigma_x1, sigma_x2_n, sigma_x3_n, &
                                rho_x1x2_n, rho_x1x3_n, rho_x2x3_n, &
                                alpha_exp, beta_exp, gamma_exp, &
                                trivar_NLL_mean_int, &
                                trivar_NLL_mean_int_c1 )


    ! Now that the results from the CLUBB code have been obtained, compare them
    ! to the MATLAB results from the same equations.
    if ( x1_opt == 1 ) then
       trivar_NLL_cmprsn_res    = -9.854231981928350e-04_dp
       trivar_NLL_c1_cmprsn_res = -9.922448817549015e-04_dp
    elseif ( x1_opt == 2 ) then
       trivar_NLL_cmprsn_res    = -1.088125285831538e-04_dp
       trivar_NLL_c1_cmprsn_res = -9.922448817549016e-05_dp
    elseif ( x1_opt == 3 ) then
       trivar_NLL_cmprsn_res    = -4.609761356905441e-05_dp
       trivar_NLL_c1_cmprsn_res = zero_dp
    else ! default
       trivar_NLL_cmprsn_res    = -9.854231981928350e-04_dp
       trivar_NLL_c1_cmprsn_res = -9.922448817549015e-04_dp
    endif

    call percent_difference( 'Trivariate NLL', &
                             trivar_NLL_mean_int, &
                             trivar_NLL_cmprsn_res, tol, count_incr )

    total_mismatches = total_mismatches + count_incr

    call percent_difference( 'Trivariate NLL (const. x1)', &
                             trivar_NLL_mean_int_c1, &
                             trivar_NLL_c1_cmprsn_res, tol, count_incr )

    total_mismatches = total_mismatches + count_incr


    !!! Bivariate NL Mean Tests
    ! The actual means, variances, and correlations of x1 and x2.
    !mu_x1 = 1.0e-3_dp
    mu_x2 = 1.0e-4_dp
    !sigma_x1 = 1.25e-4_dp
    sigma_x2 = 2.0e-5_dp
    rho_x1x2 = -0.05_dp
    if ( x1_opt == 1 ) then
       ! The individual marginal PDF of x1 is basically entirely greater than 0.
       mu_x1 = 1.0e-3_dp
       sigma_x1 = 1.25e-4_dp
    elseif ( x1_opt == 2 ) then
       ! The individual marginal PDF of x1 is mostly greater than 0,
       ! but there is a significant portion that is less than 0.
       mu_x1 = 1.0e-4_dp
       sigma_x1 = 1.25e-4_dp
    elseif ( x1_opt == 3 ) then
       ! The individual marginal PDF of x1 is slightly more heavily less
       ! than 0 than it is greater than 0.
       mu_x1 = -1.0e-7_dp
       sigma_x1 = 1.25e-4_dp
    else ! default
       ! Same as x1_opt = 1.
       mu_x1 = 1.0e-3_dp
       sigma_x1 = 1.25e-4_dp
    endif

    ! Normalize the means, variances, and correlations involving variables that
    ! have individual marginals that are distributed lognormally.
    mu_x2_n = mean_L2N_dp( mu_x2, (sigma_x2/mu_x2)**2 )

    sigma_x2_n = stdev_L2N_dp( (sigma_x2/mu_x2)**2 )

    rho_x1x2_n = corr_NL2NN_dp( rho_x1x2, sigma_x2_n, (sigma_x2/mu_x2)**2 )

    ! Exponent corresponding to x1.
    alpha_exp = 1.15_dp
    ! Exponent corresponding to x2.
    beta_exp  = 1.15_dp


    ! Obtain the results for the bivariate NL mean integral,
    ! including special cases.
    call bivar_NL_mean_tests( mu_x1, mu_x2_n, sigma_x1, sigma_x2_n, &
                              rho_x1x2_n, alpha_exp, beta_exp, &
                              bivar_NL_mean_int, bivar_NL_mean_int_c1 )


    ! Now that the results from the CLUBB code have been obtained, compare them
    ! to the MATLAB results from the same equations.
    if ( x1_opt == 1 ) then
       bivar_NL_cmprsn_res    = 8.940070763287870e-09_dp
       bivar_NL_c1_cmprsn_res = 8.942709549399743e-09_dp
    elseif ( x1_opt == 2 ) then
       bivar_NL_cmprsn_res    = 7.907812698669821e-10_dp
       bivar_NL_c1_cmprsn_res = 6.330953526469324e-10_dp
    elseif ( x1_opt == 3 ) then
       bivar_NL_cmprsn_res    = 3.254808634091251e-10_dp
       bivar_NL_c1_cmprsn_res = zero_dp
    else ! default
       bivar_NL_cmprsn_res    = 8.940070763287870e-09_dp
       bivar_NL_c1_cmprsn_res = 8.942709549399743e-09_dp
    endif

    call percent_difference( 'Bivariate NL', bivar_NL_mean_int, &
                             bivar_NL_cmprsn_res, tol, count_incr )

    total_mismatches = total_mismatches + count_incr

    call percent_difference( 'Bivariate NL (const. x1)', bivar_NL_mean_int_c1, &
                             bivar_NL_c1_cmprsn_res, tol, count_incr )

    total_mismatches = total_mismatches + count_incr

  
    !!! Bivariate LL Mean Tests
    if ( x1_opt == 1 ) then

       ! Only run the Bivariate LL Mean Test once.

       ! The actual means, variances, and correlations of x1 and x2.
       mu_x1 = 1.0e-4_dp
       mu_x2 = 100.0_dp
       sigma_x1 = 2.0e-5_dp
       sigma_x2 = 25.0_dp
       rho_x1x2 = 0.3_dp

       ! Normalize the means, variances, and correlations involving variables
       ! that have individual marginals that are distributed lognormally.
       mu_x1_n = mean_L2N_dp( mu_x1, (sigma_x1/mu_x1)**2 )
       mu_x2_n = mean_L2N_dp( mu_x2, (sigma_x2/mu_x2)**2 )

       sigma_x1_n = stdev_L2N_dp( (sigma_x1/mu_x1)**2 )
       sigma_x2_n = stdev_L2N_dp( (sigma_x2/mu_x2)**2 )

       rho_x1x2_n = corr_LL2NN_dp( rho_x1x2, sigma_x1_n, sigma_x2_n, &
                                   (sigma_x1/mu_x1)**2, (sigma_x2/mu_x2)**2 )

       ! Exponent corresponding to x1.
       alpha_exp = one_dp/3.0_dp
       ! Exponent corresponding to x2.
       beta_exp  = -one_dp/3.0_dp


       ! Obtain the results for the bivariate LL mean integral.
       call bivar_LL_mean_tests( mu_x1_n, mu_x2_n, sigma_x1_n, sigma_x2_n, &
                                 rho_x1x2_n, alpha_exp, beta_exp, &
                                 bivar_LL_mean_int )


       ! Now that the results from the CLUBB code have been obtained, compare
       ! them to the MATLAB results from the same equations.
       bivar_LL_cmprsn_res = 1.007487885941809e-02_dp

       call percent_difference( 'Bivariate LL', bivar_LL_mean_int, &
                                bivar_LL_cmprsn_res, tol, count_incr )

       total_mismatches = total_mismatches + count_incr


    endif


    return

  end subroutine KK_mean_tests

  !=============================================================================
  subroutine quadrivar_NNLL_covar_tests( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                         sigma_x1, sigma_x2, sigma_x3_n, &
                                         sigma_x4_n, rho_x1x2, rho_x1x3_n, &
                                         rho_x1x4_n, rho_x2x3_n, rho_x2x4_n, &
                                         rho_x3x4_n, x1_mean, &
                                         x2_alpha_x3_beta_x4_gamma_mean, &
                                         alpha_exp, beta_exp, gamma_exp, &
                                         quadrivar_NNLL_covar_int, &
                                         quadrivar_NNLL_covar_int_c1, &
                                         quadrivar_NNLL_covar_int_c2, &
                                         quadrivar_NNLL_covar_int_c1c2 )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use PDF_integrals_covars, only: &
        quadrivar_NNLL_covar,            & ! Procedure(s)
        quadrivar_NNLL_covar_const_x1,   &
        quadrivar_NNLL_covar_const_x2,   &
        quadrivar_NNLL_covar_const_x1x2

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

    ! Output Variables
    real( kind = dp ), intent(out) :: &
      quadrivar_NNLL_covar_int,      & ! Quadrivar. int. (all fields vary)  [-]
      quadrivar_NNLL_covar_int_c1,   & ! Quadrivar. int. (constant x1)      [-]
      quadrivar_NNLL_covar_int_c2,   & ! Quadrivar. int. (constant x2)      [-]
      quadrivar_NNLL_covar_int_c1c2    ! Quadrivar. int. (constant x1, x2)  [-]


    ! The quadrivariate NNLL covariance integral where all fields vary in
    ! the ith PDF component.
    quadrivar_NNLL_covar_int  &
    = quadrivar_NNLL_covar( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                            sigma_x1, sigma_x2, sigma_x3_n, sigma_x4_n, &
                            rho_x1x2, rho_x1x3_n, rho_x1x4_n, &
                            rho_x2x3_n, rho_x2x4_n, rho_x3x4_n, &
                            x1_mean, x2_alpha_x3_beta_x4_gamma_mean, &
                            alpha_exp, beta_exp, gamma_exp )

    ! The quadrivariate NNLL covariance integral where (only) x1 is constant in
    ! the ith PDF component.
    quadrivar_NNLL_covar_int_c1  &
    = quadrivar_NNLL_covar_const_x1( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                     sigma_x2, sigma_x3_n, sigma_x4_n, &
                                     rho_x2x3_n, rho_x2x4_n, rho_x3x4_n, &
                                     x1_mean, &
                                     x2_alpha_x3_beta_x4_gamma_mean, &
                                     alpha_exp, beta_exp, gamma_exp )

    ! The quadrivariate NNLL covariance integral where (only) x2 is constant in
    ! the ith PDF component.
    quadrivar_NNLL_covar_int_c2  &
    = quadrivar_NNLL_covar_const_x2( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                     sigma_x1, sigma_x3_n, sigma_x4_n, &
                                     rho_x1x3_n, rho_x1x4_n, rho_x3x4_n, &
                                     x1_mean, &
                                     x2_alpha_x3_beta_x4_gamma_mean, &
                                     alpha_exp, beta_exp, gamma_exp )

    ! The quadrivariate NNLL covariance integral where both x1 and x2 are
    ! constant in the ith PDF component.
    quadrivar_NNLL_covar_int_c1c2  &
    = quadrivar_NNLL_covar_const_x1x2( mu_x1, mu_x2, mu_x3_n, mu_x4_n, &
                                       sigma_x3_n, sigma_x4_n, &
                                       rho_x3x4_n, x1_mean, &
                                       x2_alpha_x3_beta_x4_gamma_mean, &
                                       alpha_exp, beta_exp, gamma_exp )


    return

  end subroutine quadrivar_NNLL_covar_tests

  !=============================================================================
  subroutine trivar_NNL_covar_tests( mu_x1, mu_x2, mu_x3_n, &
                                     sigma_x1, sigma_x2, sigma_x3_n, &
                                     rho_x1x2, rho_x1x3_n, rho_x2x3_n, &
                                     x1_mean, x2_alpha_x3_beta_mean, &
                                     alpha_exp, beta_exp, &
                                     trivar_NNL_covar_int, &
                                     trivar_NNL_covar_int_c1, &
                                     trivar_NNL_covar_int_c2, &
                                     trivar_NNL_covar_int_c1c2 )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use PDF_integrals_covars, only: &
        trivar_NNL_covar,            & ! Procedure(s)
        trivar_NNL_covar_const_x1,   &
        trivar_NNL_covar_const_x2,   &
        trivar_NNL_covar_const_x1x2

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

    ! Output Variables
    real( kind = dp ), intent(out) :: &
      trivar_NNL_covar_int,      & ! Trivariate integral (all fields vary)  [-]
      trivar_NNL_covar_int_c1,   & ! Trivariate integral (constant x1)      [-]
      trivar_NNL_covar_int_c2,   & ! Trivariate integral (constant x2)      [-]
      trivar_NNL_covar_int_c1c2    ! Trivariate integral (constant x1, x2)  [-]


    ! The trivariate NNL covariance integral where all fields vary in
    ! the ith PDF component.
    trivar_NNL_covar_int  &
    = trivar_NNL_covar( mu_x1, mu_x2, mu_x3_n, &
                        sigma_x1, sigma_x2, sigma_x3_n, &
                        rho_x1x2, rho_x1x3_n, rho_x2x3_n, &
                        x1_mean, x2_alpha_x3_beta_mean, &
                        alpha_exp, beta_exp )

    ! The trivariate NNL covariance integral where (only) x1 is constant in
    ! the ith PDF component.
    trivar_NNL_covar_int_c1  &
    = trivar_NNL_covar_const_x1( mu_x1, mu_x2, mu_x3_n, &
                                 sigma_x2, sigma_x3_n, rho_x2x3_n, &
                                 x1_mean, x2_alpha_x3_beta_mean, &
                                 alpha_exp, beta_exp )

    ! The trivariate NNL covariance integral where (only) x2 is constant in
    ! the ith PDF component.
    trivar_NNL_covar_int_c2  &
    = trivar_NNL_covar_const_x2( mu_x1, mu_x2, mu_x3_n, &
                                 sigma_x1, sigma_x3_n, rho_x1x3_n, &
                                 x1_mean, x2_alpha_x3_beta_mean, &
                                 alpha_exp, beta_exp )

    ! The trivariate NNL covariance integral where both x1 and x2 are constant
    ! in the ith PDF component.
    trivar_NNL_covar_int_c1c2  &
    = trivar_NNL_covar_const_x1x2( mu_x1, mu_x2, mu_x3_n, sigma_x3_n, &
                                   x1_mean, x2_alpha_x3_beta_mean, &
                                   alpha_exp, beta_exp )


    return

  end subroutine trivar_NNL_covar_tests

  !=============================================================================
  subroutine trivar_NLL_mean_tests( mu_x1, mu_x2_n, mu_x3_n, &
                                    sigma_x1, sigma_x2_n, sigma_x3_n, &
                                    rho_x1x2_n, rho_x1x3_n, rho_x2x3_n, &
                                    alpha_exp, beta_exp, gamma_exp, &
                                    trivar_NLL_mean_int, &
                                    trivar_NLL_mean_int_c1 )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use PDF_integrals_means, only: &
        trivar_NLL_mean, &  ! Procedure(s)
        trivar_NLL_mean_const_x1

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

    ! Output Variables
    real( kind = dp ), intent(out) :: &
      trivar_NLL_mean_int,    & ! Trivariate integral (all fields vary)     [-]
      trivar_NLL_mean_int_c1    ! Trivariate integral (constant x1)         [-]


    ! The trivariate NLL mean integral where all fields vary in
    ! the ith PDF component.
    trivar_NLL_mean_int  &
    = trivar_NLL_mean( mu_x1, mu_x2_n, mu_x3_n, &
                       sigma_x1, sigma_x2_n, sigma_x3_n, &
                       rho_x1x2_n, rho_x1x3_n, rho_x2x3_n, &
                       alpha_exp, beta_exp, gamma_exp )

    ! The trivariate NLL mean integral where x1 is constant in
    ! the ith PDF component.
    trivar_NLL_mean_int_c1  &
    = trivar_NLL_mean_const_x1( mu_x1, mu_x2_n, mu_x3_n, &
                                sigma_x2_n, sigma_x3_n, rho_x2x3_n, &
                                alpha_exp, beta_exp, gamma_exp )


    return

  end subroutine trivar_NLL_mean_tests

  !=============================================================================
  subroutine bivar_NL_mean_tests( mu_x1, mu_x2_n, sigma_x1, sigma_x2_n, &
                                  rho_x1x2_n, alpha_exp, beta_exp, &
                                  bivar_NL_mean_int, bivar_NL_mean_int_c1 )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use PDF_integrals_means, only: &
        bivar_NL_mean, &  ! Procedure(s)
        bivar_NL_mean_const_x1

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

    ! Output Variables
    real( kind = dp ), intent(out) :: &
      bivar_NL_mean_int,    & ! Bivariate integral (all fields vary)        [-]
      bivar_NL_mean_int_c1    ! Bivariate integral (constant x1)            [-]


    ! The bivariate NL mean integral where all fields vary in
    ! the ith PDF component.
    bivar_NL_mean_int  &
    = bivar_NL_mean( mu_x1, mu_x2_n, sigma_x1, sigma_x2_n, &
                     rho_x1x2_n, alpha_exp, beta_exp )

    ! The bivariate NL mean integral where x1 is constant in
    ! the ith PDF component.
    bivar_NL_mean_int_c1  &
    = bivar_NL_mean_const_x1( mu_x1, mu_x2_n, sigma_x2_n, &
                              alpha_exp, beta_exp )


    return

  end subroutine bivar_NL_mean_tests

  !=============================================================================
  subroutine bivar_LL_mean_tests( mu_x1_n, mu_x2_n, sigma_x1_n, sigma_x2_n, &
                                  rho_x1x2_n, alpha_exp, beta_exp, &
                                  bivar_LL_mean_int )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use PDF_integrals_means, only: &
        bivar_LL_mean  ! Procedure(s)

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1_n,    & ! Mean of ln x1 (ith PDF component)                     [-]
      mu_x2_n,    & ! Mean of ln x2 (ith PDF component)                     [-]
      sigma_x1_n, & ! Standard deviation of x1 (ith PDF component)          [-]
      sigma_x2_n, & ! Standard deviation of ln x2 (ith PDF component)       [-]
      rho_x1x2_n, & ! Correlation between ln x1 & ln x2 (ith PDF component) [-]
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp      ! Exponent beta, corresponding to x2                    [-]

    ! Output Variable
    real( kind = dp ), intent(out) :: &
      bivar_LL_mean_int  ! Bivariate integral (all fields vary)             [-]


    ! The bivariate LL mean integral where all fields vary in
    ! the ith PDF component.
    bivar_LL_mean_int  &
    = bivar_LL_mean( mu_x1_n, mu_x2_n, sigma_x1_n, sigma_x2_n, &
                     rho_x1x2_n, alpha_exp, beta_exp )


    return

  end subroutine bivar_LL_mean_tests

  !=============================================================================
  subroutine percent_difference( integral_type, obtained_result, &
                                 comparison_result, tol, count_incr )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        fstdout, & ! Constant(s)
        eps

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    character(len=*), intent(in) :: &
      integral_type    ! Description of integral being checked.

    real( kind = dp ), intent(in) :: &
      obtained_result,   & ! Result obtained for the integral in CLUBB   [-]
      comparison_result, & ! Result using the same equation in MATLAB    [-]
      tol                  ! Acceptable percent diff. betw. the results  [-]

    ! Output Variable
    integer, intent(out) :: &
      count_incr    ! Increment total number of mismatches.

    ! Local Variable
    real( kind = dp ) :: &
      percent_diff    ! Percent diff. betw. the CLUBB and MATLAB results    [-]


    ! Percent difference between the result obtained for the integral using
    ! the CLUBB code and the result obtained for the same integral using the
    ! same code in MATLAB.
    if ( abs(comparison_result) > eps ) then
       percent_diff = 100.0_dp * abs( ( obtained_result - comparison_result ) &
                                      / comparison_result )
    else  ! Integral should have a value of 0.
       percent_diff = 100.0_dp * abs( ( obtained_result - comparison_result ) &
                                      / epsilon( comparison_result ) )
    endif


    if ( percent_diff <= tol ) then

       ! The percent difference between the obtained result and the MATLAB
       ! result is within an acceptable tolerance.
       write(fstdout,'(A,A)') "Agreement with regards to integral:  ", &
                              trim( integral_type )
       write(fstdout,'(A,ES23.15)') "Obtained result:  ", obtained_result
       write(fstdout,'(A,ES23.15)') "MATLAB result:  ", comparison_result
       write(fstdout,'(A,ES12.5,A)') "Percent difference:  ", percent_diff, "%"
       write(fstdout,'(A)') " "

       count_incr = 0

    else ! percent_diff > tol

       ! The percent difference between the obtained result and the MATLAB
       ! result is beyond an acceptable tolerance.  Print an error message
       ! telling the user to please check for any changes made to the
       ! relevant portion(s) of the CLUBB model code.
       write(fstdout,'(A)') "###################################"
       write(fstdout,'(A,A)') "Mismatch with regards to integral:  ", &
                              trim( integral_type )
       write(fstdout,'(A)') "The percent difference between the obtained "  &
                            // "result and the MATLAB result is beyond an "  &
                            // "acceptable tolerance.  Please check for any "  &
                            // "changes made to the relevant portion(s) of " &
                            // "the CLUBB model code."
       write(fstdout,'(A,ES23.15)') "Obtained result:  ", obtained_result
       write(fstdout,'(A,ES23.15)') "MATLAB result:  ", comparison_result
       write(fstdout,'(A,ES12.5,A)') "Percent difference:  ", percent_diff, "%"
       write(fstdout,'(A)') " "

       count_incr = 1

    endif


    return

    end subroutine percent_difference

!===============================================================================

end module KK_integrals_tests
