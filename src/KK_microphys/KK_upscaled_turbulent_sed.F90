! $Id$
!===============================================================================
module KK_upscaled_turbulent_sed

  implicit none

  private  ! Default scope

  public :: KK_sed_vel_covars

  private :: covar_rr_KK_mvr_coefA, &
             covar_rr_KK_mvr_termB, &
             covar_rr_KK_mvr, &
             covar_Nr_KK_mvr_coefA, &
             covar_Nr_KK_mvr_termB, &
             covar_Nr_KK_mvr, &
             bivar_LL_covar_partial_rr, &
             bivar_LL_covar_partial_Nr, &
             bivar_LL_covar_partial, &
             bivar_LL_covar_const_x2_partial

  contains

  !=============================================================================
  subroutine KK_sed_vel_covars( rrm, rr_1, rr_2, Nrm, &
                                Nr_1, Nr_2, KK_mean_vol_rad, &
                                mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, &
                                mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1, &
                                sigma_rr_2, sigma_Nr_1, sigma_Nr_2, &
                                sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                                sigma_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, &
                                KK_mvr_coef, mixt_frac, precip_frac_1, &
                                precip_frac_2, level, l_stats_samp, &
                                Vrrprrp_impc, Vrrprrp_expc, &
                                VNrpNrp_impc, VNrpNrp_expc )

    ! Description:
    ! The covariance of sedimentation velocity (of rain water mixing ratio) and
    ! rain water mixing ratio, < V_rr'r_r' > is given by the equation:
    !
    ! < V_rr'r_r' > = 0.012 ( 10^6 microm/m ) < r_r'R_vr' >.
    !
    ! Likewise, the covariance of sedimentation velocity (of rain drop
    ! concentration) and rain drop concentration, < V_Nr'N_r' > is given by the
    ! equation:
    !
    ! < V_Nr'N_r' > = 0.007 ( 10^6 microm/m ) < N_r'R_vr' >.
    !
    ! The covariances < r_r'R_vr' > and < N_r'R_vr' > are calculated in this
    ! subroutine.
    !
    ! It is advantageous to write the < V_rr'r_r' > term and the < V_Nr'N_r' >
    ! term in terms of < r_r > and < N_r >, respectively, because < V_rr'r_r' >
    ! and < V_Nr'N_r' > are found in the turbulent sedimentation term as part
    ! of the predictive equations for < r_r > and < N_r >, respectively.  This
    ! allows the turbulent sedimentation term to be handled semi-implicitly in
    ! the < r_r > and < N_r > predictive equations, increasing the numerical
    ! stability of the solutions.  However, the option to handle these terms
    ! in the traditional, completely explicit and less computationally-expensive
    ! manner is included.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        micron_per_m, & ! Constant(s)
        zero

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use stats_type_utilities, only: &
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only: & 
        irr_KK_mvr_covar_zt, & ! Variable(s)
        iNr_KK_mvr_covar_zt, &
        stats_zt

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rrm,             & ! Mean rain water mixing ratio (overall)        [kg/kg]
      rr_1,            & ! Mean rain water mixing ratio (1st PDF comp.)  [kg/kg]
      rr_2,            & ! Mean rain water mixing ratio (2nd PDF comp.)  [kg/kg]
      Nrm,             & ! Mean rain drop concentration (overall)       [num/kg]
      Nr_1,            & ! Mean rain drop concentration (1st PDF comp.) [num/kg]
      Nr_2,            & ! Mean rain drop concentration (2nd PDF comp.) [num/kg]
      KK_mean_vol_rad, & ! KK mean volume radius of rain drops               [m]
      mu_rr_1,         & ! Mean of rr (1st PDF component) in-precip (ip) [kg/kg]
      mu_rr_2,         & ! Mean of rr (2nd PDF component) ip             [kg/kg]
      mu_Nr_1,         & ! Mean of Nr (1st PDF component) ip            [num/kg]
      mu_Nr_2,         & ! Mean of Nr (2nd PDF component) ip            [num/kg]
      mu_rr_1_n,       & ! Mean of ln rr (1st PDF component) ip      [ln(kg/kg)]
      mu_rr_2_n,       & ! Mean of ln rr (2nd PDF component) ip      [ln(kg/kg)]
      mu_Nr_1_n,       & ! Mean of ln Nr (1st PDF component) ip     [ln(num/kg)]
      mu_Nr_2_n,       & ! Mean of ln Nr (2nd PDF component) ip     [ln(num/kg)]
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
      KK_mvr_coef,     & ! KK mean volume radius coefficient                 [m]
      mixt_frac,       & ! Mixture fraction                                  [-]
      precip_frac_1,   & ! Precipitation fraction (1st PDF component)        [-]
      precip_frac_2      ! Precipitation fraction (2nd PDF component)        [-]

    integer, intent(in) :: &
      level   ! Vertical level index 

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      Vrrprrp_impc, & ! Implicit comp. of <V_rr'r_r'>: <r_r> eq. [(m/s)]
      Vrrprrp_expc, & ! Explicit comp. of <V_rr'r_r'>: <r_r> eq. [(m/s)(kg/kg)]
      VNrpNrp_impc, & ! Implicit comp. of <V_Nr'N_r'>: <N_r> eq. [(m/s)]
      VNrpNrp_expc    ! Explicit comp. of <V_Nr'N_r'>: <N_r> eq. [(m/s)(num/kg)]

    ! Local Variables
    real( kind = core_rknd ) :: &
      rr_KK_mvr_covar_coefA, & ! Coefficient of <r_r> in < r_r'R_vr' > eq.   [m]
      rr_KK_mvr_covar_termB, & ! Term in < r_r'R_vr' > equation       [(kg/kg)m]
      Nr_KK_mvr_covar_coefA, & ! Coefficient of <N_r> in < N_r'R_vr' > eq.   [m]
      Nr_KK_mvr_covar_termB    ! Term in < N_r'R_vr' > equation      [(num/kg)m]

    real( kind = core_rknd ) :: &
      rr_KK_mvr_covar, & ! Covariance of r_r and KK mean vol rad     [(kg/kg)m]
      Nr_KK_mvr_covar    ! Covariance of N_r and KK mean vol rad     [(num/kg)m]

    logical, parameter :: &
      l_semi_imp_turbulent_sed = .true. ! Flag to use semi-implicit tur. sed.


    ! Calculate the covariance of rain drop mean volume radius and r_r,
    ! < R_vr'r_r' >.
    if ( l_semi_imp_turbulent_sed ) then

       ! Turbulent sedimentation will be handled semi-implicitly in the
       ! hydrometeor predictive equation set.

       rr_KK_mvr_covar_coefA  &
       = covar_rr_KK_mvr_coefA( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
                                mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                                sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, &
                                sigma_Nr_1_n, sigma_Nr_2_n, corr_rr_Nr_1_n, &
                                corr_rr_Nr_2_n, KK_mean_vol_rad, KK_mvr_coef )

       rr_KK_mvr_covar_termB  &
       = covar_rr_KK_mvr_termB( rr_1, rr_2, mu_rr_1, mu_rr_2, mu_Nr_1, &
                                mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                mu_Nr_2_n, sigma_rr_1, sigma_rr_2, &
                                sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                                sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                                corr_rr_Nr_1_n, corr_rr_Nr_2_n, KK_mvr_coef, &
                                mixt_frac )

    else

       ! Turbulent sedimentation will be handled completely explicitly in the
       ! hydrometeor predictive equation set.

       rr_KK_mvr_covar_coefA = zero

       rr_KK_mvr_covar_termB  &
       = covar_rr_KK_mvr( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, &
                          mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1, &
                          sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                          sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                          corr_rr_Nr_1_n, corr_rr_Nr_2_n, rrm, &
                          KK_mean_vol_rad, KK_mvr_coef, mixt_frac, &
                          precip_frac_1, precip_frac_2 )

    endif


    ! Calculate the covariance of rain drop mean volume radius and N_r,
    ! < R_vr'N_r' >.
    if ( l_semi_imp_turbulent_sed ) then

       ! Turbulent sedimentation will be handled semi-implicitly in the
       ! hydrometeor predictive equation set.

       Nr_KK_mvr_covar_coefA  &
       = covar_Nr_KK_mvr_coefA( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
                                mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                                sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, &
                                sigma_Nr_1_n, sigma_Nr_2_n, corr_rr_Nr_1_n, &
                                corr_rr_Nr_2_n, KK_mean_vol_rad, KK_mvr_coef )

       Nr_KK_mvr_covar_termB  &
       = covar_Nr_KK_mvr_termB( Nr_1, Nr_2, mu_rr_1, mu_rr_2, mu_Nr_1, &
                                mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                mu_Nr_2_n, sigma_rr_1, sigma_rr_2,  &
                                sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                                sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                                corr_rr_Nr_1_n, corr_rr_Nr_2_n, KK_mvr_coef, &
                                mixt_frac )

    else

       ! Turbulent sedimentation will be handled completely explicitly in the
       ! hydrometeor predictive equation set.

       Nr_KK_mvr_covar_coefA = zero

       Nr_KK_mvr_covar_termB  &
       = covar_Nr_KK_mvr( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, &
                          mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1, &
                          sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                          sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                          corr_rr_Nr_1_n, corr_rr_Nr_2_n, Nrm, &
                          KK_mean_vol_rad, KK_mvr_coef, mixt_frac, &
                          precip_frac_1, precip_frac_2 )

    endif


    ! Terms used to calculate the covariance of V_rr and r_r, < V_rr'r_r' >.
    Vrrprrp_impc = - 0.012_core_rknd * micron_per_m * rr_KK_mvr_covar_coefA
    Vrrprrp_expc = - 0.012_core_rknd * micron_per_m * rr_KK_mvr_covar_termB

    ! Terms used to calculate the covariance of V_Nr and N_r, < V_Nr'N_r' >.
    VNrpNrp_impc = - 0.007_core_rknd * micron_per_m * Nr_KK_mvr_covar_coefA
    VNrpNrp_expc = - 0.007_core_rknd * micron_per_m * Nr_KK_mvr_covar_termB


    ! Statistics
    if ( l_stats_samp ) then

       ! Covariance of r_r and KK rain drop mean volume radius.
       if ( irr_KK_mvr_covar_zt > 0 ) then

          rr_KK_mvr_covar &
          = rr_KK_mvr_covar_coefA * rrm + rr_KK_mvr_covar_termB

          call stat_update_var_pt( irr_KK_mvr_covar_zt, level, &
                                   rr_KK_mvr_covar, stats_zt )

       endif

       ! Covariance of N_r and KK rain drop mean volume radius.
       if ( iNr_KK_mvr_covar_zt > 0 ) then

          Nr_KK_mvr_covar &
          = Nr_KK_mvr_covar_coefA * Nrm + Nr_KK_mvr_covar_termB

          call stat_update_var_pt( iNr_KK_mvr_covar_zt, level, &
                                   Nr_KK_mvr_covar, stats_zt )

       endif

    endif ! l_stats_samp


    return

  end subroutine KK_sed_vel_covars

  !=============================================================================
  function covar_rr_KK_mvr_coefA( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
                                  mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                                  sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                  sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, &
                                  sigma_Nr_1_n, sigma_Nr_2_n, corr_rr_Nr_1_n, &
                                  corr_rr_Nr_2_n, KK_mean_vol_rad, KK_mvr_coef )

    ! Description:
    ! This function partially calculates the covariance of r_r and KK mean
    ! volume radius of rain drops (R_vr), which can be written as < r_r'R_vr' >.
    !
    ! The equation for the covariance of r_r and KK rain drop mean volume radius
    ! is:
    !
    ! < r_r'R_vr' >
    ! = KK_mvr_coef
    !   ( a f_p(1)
    !     exp{ mu_rr_1_n ( alpha + 1 ) + mu_Nr_1_n beta
    !          + (1/2) sigma_rr_1_n^2 ( alpha + 1 )^2
    !          + (1/2) sigma_Nr_1_n^2 beta^2
    !          + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
    !   + ( 1 - a ) f_p(2)
    !     exp{ mu_rr_2_n ( alpha + 1 ) + mu_Nr_2_n beta
    !          + (1/2) sigma_rr_2_n^2 ( alpha + 1 )^2
    !          + (1/2) sigma_Nr_2_n^2 beta^2
    !          + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
    !   )
    !   - < R_vr > < r_r >;
    !
    ! where "a" is the mixture fraction (weight of the 1st PDF component), and
    ! f_p(1) and f_p(2) are the precipitation fractions in the 1st and 2nd PDF
    ! components, respectively.  The equation for < r_r'R_vr' > needs to be put
    ! in the form:  < r_r'R_vr' > = coef_A < r_r > + term_B.  Since < r_r > can
    ! be written as:
    !
    ! < r_r > = a f_p(1) exp{ mu_rr_1_n + (1/2) sigma_rr_1_n^2 }
    !           + ( 1 - a ) f_p(2) exp{ mu_rr_2_n + (1/2) sigma_rr_2_n^2 }
    !         = a f_p(1) mu_rr_1 + ( 1 - a ) f_p(2) mu_rr_2
    !         = a rr_1 + ( 1 - a ) rr_2;
    !
    ! the covariance < r_r'R_vr' > can be rewritten in terms of < r_r >:
    !
    ! < r_r'R_vr' >
    ! = ( KK_mvr_coef
    !     ( exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
    !            + (1/2) sigma_rr_1_n^2 ( alpha^2 + 2 alpha )
    !            + (1/2) sigma_Nr_1_n^2 beta^2
    !            + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
    !     + exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
    !            + (1/2) sigma_rr_2_n^2 ( alpha^2 + 2 alpha )
    !            + (1/2) sigma_Nr_2_n^2 beta^2
    !            + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
    !     )
    !     - < R_vr > ) < r_r >
    !   - KK_mvr_coef
    !     ( a rr_1
    !       exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
    !            + (1/2) sigma_rr_2_n^2 ( alpha^2 + 2 alpha )
    !            + (1/2) sigma_Nr_2_n^2 beta^2
    !            + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
    !     + ( 1 - a ) rr_2
    !       exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
    !            + (1/2) sigma_rr_1_n^2 ( alpha^2 + 2 alpha )
    !            + (1/2) sigma_Nr_1_n^2 beta^2
    !            + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
    !     ).
    !
    ! This function calculates the coefficient A, which is given by:
    !
    ! coef_A
    ! = KK_mvr_coef
    !   ( exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
    !          + (1/2) sigma_rr_1_n^2 ( alpha^2 + 2 alpha )
    !          + (1/2) sigma_Nr_1_n^2 beta^2
    !          + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
    !   + exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
    !          + (1/2) sigma_rr_2_n^2 ( alpha^2 + 2 alpha )
    !          + (1/2) sigma_Nr_2_n^2 beta^2
    !          + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
    !   )
    !   - < R_vr >.
    !
    ! The covariance of sedimentation velocity (of rain water mixing ratio) and
    ! rain water mixing ratio, < V_rr'r_r' > is given by the equation:
    !
    ! < V_rr'r_r' > = 0.012 ( 10^6 microm/m ) < r_r'R_vr' >.
    !
    ! It is advantageous to write the < V_rr'r_r' > term in terms of < r_r >
    ! because < V_rr'r_r' > is found in the turbulent sedimentation term as part
    ! of the predictive equation for < r_r >.  This allows the turbulent
    ! sedimentation term to be handled semi-implicitly in the < r_r > predictive
    ! equation, increasing the numerical stability of the solution.

    ! References:
    !-----------------------------------------------------------------------

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
      mu_rr_1_n,       & ! Mean of ln rr (1st PDF component) ip      [ln(kg/kg)]
      mu_rr_2_n,       & ! Mean of ln rr (2nd PDF component) ip      [ln(kg/kg)]
      mu_Nr_1_n,       & ! Mean of ln Nr (1st PDF component) ip     [ln(num/kg)]
      mu_Nr_2_n,       & ! Mean of ln Nr (2nd PDF component) ip     [ln(num/kg)]
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
      KK_mvr_coef        ! KK mean volume radius coefficient                 [m]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_rr_KK_mvr_coefA    ! Coefficient of <r_r> in < r_r'R_vr' > eq.   [m]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on r_r
      beta_exp     ! Exponent on N_r


    ! Values of the KK exponents.
    alpha_exp = KK_mvr_rr_exp
    beta_exp  = KK_mvr_Nr_exp

    ! Calculate the coefficient A in the < r_r'R_vr' > covariance equation,
    ! where < r_r'R_vr' > = coef_A <r_r> + term_B, and:
    !
    ! coef_A
    ! = KK_mvr_coef
    !   ( exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
    !          + (1/2) sigma_rr_1_n^2 ( alpha^2 + 2 alpha )
    !          + (1/2) sigma_Nr_1_n^2 beta^2
    !          + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
    !   + exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
    !          + (1/2) sigma_rr_2_n^2 ( alpha^2 + 2 alpha )
    !          + (1/2) sigma_Nr_2_n^2 beta^2
    !          + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
    !   )
    !   - < R_vr >.
    covar_rr_KK_mvr_coefA  &
    = KK_mvr_coef &
      * ( bivar_LL_covar_partial_rr( mu_rr_1, mu_Nr_1, mu_rr_1_n, &
                                     mu_Nr_1_n, sigma_rr_1, sigma_Nr_1, &
                                     sigma_rr_1_n, sigma_Nr_1_n, &
                                     corr_rr_Nr_1_n, alpha_exp, beta_exp ) &
        + bivar_LL_covar_partial_rr( mu_rr_2, mu_Nr_2, mu_rr_2_n, &
                                     mu_Nr_2_n, sigma_rr_2, sigma_Nr_2, &
                                     sigma_rr_2_n, sigma_Nr_2_n, &
                                     corr_rr_Nr_2_n, alpha_exp, beta_exp ) &
        ) &
      - KK_mean_vol_rad


    return

  end function covar_rr_KK_mvr_coefA

  !=============================================================================
  function covar_rr_KK_mvr_termB( rr_1, rr_2, mu_rr_1, mu_rr_2, mu_Nr_1, &
                                  mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                  mu_Nr_2_n, sigma_rr_1, sigma_rr_2, &
                                  sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                                  sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                                  corr_rr_Nr_1_n, corr_rr_Nr_2_n, KK_mvr_coef, &
                                  mixt_frac )

    ! Description:
    ! This function partially calculates the covariance of r_r and KK mean
    ! volume radius of rain drops (R_vr), which can be written as < r_r'R_vr' >.
    !
    ! The equation for the covariance of r_r and KK rain drop mean volume radius
    ! is:
    !
    ! < r_r'R_vr' >
    ! = KK_mvr_coef
    !   ( a f_p(1)
    !     exp{ mu_rr_1_n ( alpha + 1 ) + mu_Nr_1_n beta
    !          + (1/2) sigma_rr_1_n^2 ( alpha + 1 )^2
    !          + (1/2) sigma_Nr_1_n^2 beta^2
    !          + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
    !   + ( 1 - a ) f_p(2)
    !     exp{ mu_rr_2_n ( alpha + 1 ) + mu_Nr_2_n beta
    !          + (1/2) sigma_rr_2_n^2 ( alpha + 1 )^2
    !          + (1/2) sigma_Nr_2_n^2 beta^2
    !          + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
    !   )
    !   - < R_vr > < r_r >;
    !
    ! where "a" is the mixture fraction (weight of the 1st PDF component), and
    ! f_p(1) and f_p(2) are the precipitation fractions in the 1st and 2nd PDF
    ! components, respectively.  The equation for < r_r'R_vr' > needs to be put
    ! in the form:  < r_r'R_vr' > = coef_A < r_r > + term_B.  Since < r_r > can
    ! be written as:
    !
    ! < r_r > = a f_p(1) exp{ mu_rr_1_n + (1/2) sigma_rr_1_n^2 }
    !           + ( 1 - a ) f_p(2) exp{ mu_rr_2_n + (1/2) sigma_rr_2_n^2 }
    !         = a f_p(1) mu_rr_1 + ( 1 - a ) f_p(2) mu_rr_2
    !         = a rr_1 + ( 1 - a ) rr_2;
    !
    ! the covariance < r_r'R_vr' > can be rewritten in terms of < r_r >:
    !
    ! < r_r'R_vr' >
    ! = ( KK_mvr_coef
    !     ( exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
    !            + (1/2) sigma_rr_1_n^2 ( alpha^2 + 2 alpha )
    !            + (1/2) sigma_Nr_1_n^2 beta^2
    !            + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
    !     + exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
    !            + (1/2) sigma_rr_2_n^2 ( alpha^2 + 2 alpha )
    !            + (1/2) sigma_Nr_2_n^2 beta^2
    !            + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
    !     )
    !     - < R_vr > ) < r_r >
    !   - KK_mvr_coef
    !     ( a rr_1
    !       exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
    !            + (1/2) sigma_rr_2_n^2 ( alpha^2 + 2 alpha )
    !            + (1/2) sigma_Nr_2_n^2 beta^2
    !            + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
    !     + ( 1 - a ) rr_2
    !       exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
    !            + (1/2) sigma_rr_1_n^2 ( alpha^2 + 2 alpha )
    !            + (1/2) sigma_Nr_1_n^2 beta^2
    !            + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
    !     ).
    !
    ! This function calculates term B, which is given by:
    !
    ! term_B
    ! = - KK_mvr_coef
    !     ( a rr_1
    !       exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
    !            + (1/2) sigma_rr_2_n^2 ( alpha^2 + 2 alpha )
    !            + (1/2) sigma_Nr_2_n^2 beta^2
    !            + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
    !     + ( 1 - a ) rr_2
    !       exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
    !            + (1/2) sigma_rr_1_n^2 ( alpha^2 + 2 alpha )
    !            + (1/2) sigma_Nr_1_n^2 beta^2
    !            + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
    !     ).
    !
    ! The covariance of sedimentation velocity (of rain water mixing ratio) and
    ! rain water mixing ratio, < V_rr'r_r' > is given by the equation:
    !
    ! < V_rr'r_r' > = 0.012 ( 10^6 microm/m ) < r_r'R_vr' >.
    !
    ! It is advantageous to write the < V_rr'r_r' > term in terms of < r_r >
    ! because < V_rr'r_r' > is found in the turbulent sedimentation term as part
    ! of the predictive equation for < r_r >.  This allows the turbulent
    ! sedimentation term to be handled semi-implicitly in the < r_r > predictive
    ! equation, increasing the numerical stability of the solution.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one    ! Constant(s)

    use parameters_KK, only: &
        KK_mvr_rr_exp, & ! Variable(s)
        KK_mvr_Nr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rr_1,           & ! Mean rain water mixing ratio (1st PDF comp.)   [kg/kg]
      rr_2,           & ! Mean rain water mixing ratio (2nd PDF comp.)   [kg/kg]
      mu_rr_1,        & ! Mean of rr (1st PDF component) in-precip (ip)  [kg/kg]
      mu_rr_2,        & ! Mean of rr (2nd PDF component) ip              [kg/kg]
      mu_Nr_1,        & ! Mean of Nr (1st PDF component) ip             [num/kg]
      mu_Nr_2,        & ! Mean of Nr (2nd PDF component) ip             [num/kg]
      mu_rr_1_n,      & ! Mean of ln rr (1st PDF component) ip       [ln(kg/kg)]
      mu_rr_2_n,      & ! Mean of ln rr (2nd PDF component) ip       [ln(kg/kg)]
      mu_Nr_1_n,      & ! Mean of ln Nr (1st PDF component) ip      [ln(num/kg)]
      mu_Nr_2_n,      & ! Mean of ln Nr (2nd PDF component) ip      [ln(num/kg)]
      sigma_rr_1,     & ! Standard deviation of rr (1st PDF comp.) ip    [kg/kg]
      sigma_rr_2,     & ! Standard deviation of rr (2nd PDF comp.) ip    [kg/kg]
      sigma_Nr_1,     & ! Standard deviation of Nr (1st PDF comp.) ip   [num/kg]
      sigma_Nr_2,     & ! Standard deviation of Nr (2nd PDF comp.) ip   [num/kg]
      sigma_rr_1_n,   & ! Standard deviation of ln rr (1st PDF component) ip [-]
      sigma_rr_2_n,   & ! Standard deviation of ln rr (2nd PDF component) ip [-]
      sigma_Nr_1_n,   & ! Standard deviation of ln Nr (1st PDF component) ip [-]
      sigma_Nr_2_n,   & ! Standard deviation of ln Nr (2nd PDF component) ip [-]
      corr_rr_Nr_1_n, & ! Correlation of ln rr & ln Nr (1st PDF comp.) ip    [-]
      corr_rr_Nr_2_n, & ! Correlation of ln rr & ln Nr (2nd PDF comp.) ip    [-]
      KK_mvr_coef,    & ! KK mean volume radius coefficient                  [m]
      mixt_frac         ! Mixture fraction                                   [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_rr_KK_mvr_termB    ! Term in < r_r'R_vr' > equation       [(kg/kg)m]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on r_r
      beta_exp     ! Exponent on N_r


    ! Values of the KK exponents.
    alpha_exp = KK_mvr_rr_exp
    beta_exp  = KK_mvr_Nr_exp

    ! Calculate term B in the < r_r'R_vr' > covariance equation, where
    ! < r_r'R_vr' > = coef_A <r_r> + term_B, and:
    !
    ! term_B
    ! = - KK_mvr_coef
    !     ( a rr_1
    !       exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
    !            + (1/2) sigma_rr_2_n^2 ( alpha^2 + 2 alpha )
    !            + (1/2) sigma_Nr_2_n^2 beta^2
    !            + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
    !     + ( 1 - a ) rr_2
    !       exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
    !            + (1/2) sigma_rr_1_n^2 ( alpha^2 + 2 alpha )
    !            + (1/2) sigma_Nr_1_n^2 beta^2
    !            + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
    !     ).
    covar_rr_KK_mvr_termB  &
    = - KK_mvr_coef &
        * ( mixt_frac * rr_1 &
            * bivar_LL_covar_partial_rr( mu_rr_2, mu_Nr_2, mu_rr_2_n, &
                                         mu_Nr_2_n, sigma_rr_2, sigma_Nr_2, &
                                         sigma_rr_2_n, sigma_Nr_2_n, &
                                         corr_rr_Nr_2_n, alpha_exp, beta_exp ) &
          + ( one - mixt_frac ) * rr_2 &
            * bivar_LL_covar_partial_rr( mu_rr_1, mu_Nr_1, mu_rr_1_n, &
                                         mu_Nr_1_n, sigma_rr_1, sigma_Nr_1, &
                                         sigma_rr_1_n, sigma_Nr_1_n, &
                                         corr_rr_Nr_1_n, alpha_exp, beta_exp ) &
          )


    return

  end function covar_rr_KK_mvr_termB

  !=============================================================================
  function covar_rr_KK_mvr( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, &
                            mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1, &
                            sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                            sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                            corr_rr_Nr_1_n, corr_rr_Nr_2_n, rrm, &
                            KK_mean_vol_rad, KK_mvr_coef, mixt_frac, &
                            precip_frac_1, precip_frac_2 )

    ! Description:
    ! This function calculates the covariance of r_r and KK mean volume radius
    ! of rain drops (R_vr), which can be written as < r_r'R_vr' >.
    !
    ! The equation for the covariance of r_r and KK rain drop mean volume radius
    ! is:
    !
    ! < r_r'R_vr' >
    ! = KK_mvr_coef
    !   ( a f_p(1)
    !     exp{ mu_rr_1_n ( alpha + 1 ) + mu_Nr_1_n beta
    !          + (1/2) sigma_rr_1_n^2 ( alpha + 1 )^2
    !          + (1/2) sigma_Nr_1_n^2 beta^2
    !          + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
    !   + ( 1 - a ) f_p(2)
    !     exp{ mu_rr_2_n ( alpha + 1 ) + mu_Nr_2_n beta
    !          + (1/2) sigma_rr_2_n^2 ( alpha + 1 )^2
    !          + (1/2) sigma_Nr_2_n^2 beta^2
    !          + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
    !   )
    !   - < R_vr > < r_r >;
    !
    ! where "a" is the mixture fraction (weight of the 1st PDF component), and
    ! f_p(1) and f_p(2) are the precipitation fractions in the 1st and 2nd PDF
    ! components, respectively.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one  ! Constant(s)

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
      mu_rr_1_n,       & ! Mean of ln rr (1st PDF component) ip      [ln(kg/kg)]
      mu_rr_2_n,       & ! Mean of ln rr (2nd PDF component) ip      [ln(kg/kg)]
      mu_Nr_1_n,       & ! Mean of ln Nr (1st PDF component) ip     [ln(num/kg)]
      mu_Nr_2_n,       & ! Mean of ln Nr (2nd PDF component) ip     [ln(num/kg)]
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
      rrm,             & ! Mean rain water mixing ratio                  [kg/kg]
      KK_mean_vol_rad, & ! KK mean volume radius of rain drops               [m]
      KK_mvr_coef,     & ! KK mean volume radius coefficient                 [m]
      mixt_frac,       & ! Mixture fraction                                  [-]
      precip_frac_1,   & ! Precipitation fraction (1st PDF component)        [-]
      precip_frac_2      ! Precipitation fraction (2nd PDF component)        [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_rr_KK_mvr    ! Covar of rr and KK rain drop mean vol rad  [(kg/kg)m]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on r_r
      beta_exp     ! Exponent on N_r


    ! Values of the KK exponents.
    alpha_exp = KK_mvr_rr_exp
    beta_exp  = KK_mvr_Nr_exp

    ! Calculate the covariance of r_r and KK mean volume radius of rain drops.
    covar_rr_KK_mvr  &
    = KK_mvr_coef &
      * ( mixt_frac &
          * precip_frac_1 &
          * bivar_LL_mean_eq( mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, &
                              sigma_rr_1, sigma_Nr_1, sigma_rr_1_n, &
                              sigma_Nr_1_n, corr_rr_Nr_1_n, &
                              alpha_exp + one, beta_exp ) &
        + ( one - mixt_frac ) &
          * precip_frac_2 &
          * bivar_LL_mean_eq( mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, &
                              sigma_rr_2, sigma_Nr_2, sigma_rr_2_n, &
                              sigma_Nr_2_n, corr_rr_Nr_2_n, &
                              alpha_exp + one, beta_exp ) &
        ) &
      - rrm * KK_mean_vol_rad


    return

  end function covar_rr_KK_mvr

  !=============================================================================
  function covar_Nr_KK_mvr_coefA( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
                                  mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                                  sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                  sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, &
                                  sigma_Nr_1_n, sigma_Nr_2_n, corr_rr_Nr_1_n, &
                                  corr_rr_Nr_2_n, KK_mean_vol_rad, KK_mvr_coef )

    ! Description:
    ! This function partially calculates the covariance of N_r and KK mean
    ! volume radius of rain drops (R_vr), which can be written as < N_r'R_vr' >.
    !
    ! The equation for the covariance of N_r and KK rain drop mean volume radius
    ! is:
    !
    ! < N_r'R_vr' >
    ! = KK_mvr_coef
    !   ( a f_p(1)
    !     exp{ mu_rr_1_n alpha + mu_Nr_1_n ( beta + 1 )
    !          + (1/2) sigma_rr_1_n^2 alpha^2
    !          + (1/2) sigma_Nr_1_n^2 ( beta + 1 )^2
    !          + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
    !   + ( 1 - a ) f_p(2)
    !     exp{ mu_rr_2_n alpha + mu_Nr_2_n ( beta + 1 )
    !          + (1/2) sigma_rr_2_n^2 alpha^2
    !          + (1/2) sigma_Nr_2_n^2 ( beta + 1 )^2
    !          + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
    !   )
    !   - < R_vr > < N_r >;
    !
    ! where "a" is the mixture fraction (weight of the 1st PDF component), and
    ! f_p(1) and f_p(2) are the precipitation fractions in the 1st and 2nd PDF
    ! components, respectively.  The equation for < N_r'R_vr' > needs to be put
    ! in the form:  < N_r'R_vr' > = coef_A < N_r > + term_B.  Since < N_r > can
    ! be written as:
    !
    ! < N_r > = a f_p(1) exp{ mu_Nr_1_n + (1/2) sigma_Nr_1_n^2 }
    !           + ( 1 - a ) f_p(2) exp{ mu_Nr_2_n + (1/2) sigma_Nr_2_n^2 }
    !         = a f_p(1) mu_Nr_1 + ( 1 - a ) f_p(2) mu_Nr_2
    !         = a Nr_1 + ( 1 - a ) Nr_2;
    !
    ! the covariance < N_r'R_vr' > can be rewritten in terms of < N_r >:
    !
    ! < N_r'R_vr' >
    ! = ( KK_mvr_coef
    !     ( exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
    !            + (1/2) sigma_rr_1_n^2 alpha^2
    !            + (1/2) sigma_Nr_1_n^2 ( beta^2 + 2 beta )
    !            + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
    !     + exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
    !            + (1/2) sigma_rr_2_n^2 alpha^2
    !            + (1/2) sigma_Nr_2_n^2 ( beta^2 + 2 beta )
    !            + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
    !     )
    !     - < R_vr > ) < N_r >
    !   - KK_mvr_coef
    !     ( a Nr_1
    !       exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
    !            + (1/2) sigma_rr_2_n^2 alpha^2
    !            + (1/2) sigma_Nr_2_n^2 ( beta^2 + 2 beta )
    !            + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
    !     + ( 1 - a ) Nr_2
    !       exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
    !            + (1/2) sigma_rr_1_n^2 alpha^2
    !            + (1/2) sigma_Nr_1_n^2 ( beta^2 + 2 beta )
    !            + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
    !     ).
    !
    ! This function calculates the coefficient A, which is given by:
    !
    ! coef_A
    ! = KK_mvr_coef
    !   ( exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
    !          + (1/2) sigma_rr_1_n^2 alpha^2
    !          + (1/2) sigma_Nr_1_n^2 ( beta^2 + 2 beta )
    !          + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
    !   + exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
    !          + (1/2) sigma_rr_2_n^2 alpha^2
    !          + (1/2) sigma_Nr_2_n^2 ( beta^2 + 2 beta )
    !          + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
    !   )
    !   - < R_vr >.
    !
    ! The covariance of sedimentation velocity (of rain drop concentration) and
    ! rain drop concentration, < V_Nr'N_r' > is given by the equation:
    !
    ! < V_Nr'N_r' > = 0.007 ( 10^6 microm/m ) < N_r'R_vr' >.
    !
    ! It is advantageous to write the < V_Nr'N_r' > term in terms of < N_r >
    ! because < V_Nr'N_r' > is found in the turbulent sedimentation term as part
    ! of the predictive equation for < N_r >.  This allows the turbulent
    ! sedimentation term to be handled semi-implicitly in the < N_r > predictive
    ! equation, increasing the numerical stability of the solution.

    ! References:
    !-----------------------------------------------------------------------

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
      mu_rr_1_n,       & ! Mean of ln rr (1st PDF component) ip      [ln(kg/kg)]
      mu_rr_2_n,       & ! Mean of ln rr (2nd PDF component) ip      [ln(kg/kg)]
      mu_Nr_1_n,       & ! Mean of ln Nr (1st PDF component) ip     [ln(num/kg)]
      mu_Nr_2_n,       & ! Mean of ln Nr (2nd PDF component) ip     [ln(num/kg)]
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
      KK_mvr_coef        ! KK mean volume radius coefficient                 [m]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_Nr_KK_mvr_coefA    ! Coefficient of <N_r> in < r_r'R_vr' > eq.   [m]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on r_r
      beta_exp     ! Exponent on N_r


    ! Values of the KK exponents.
    alpha_exp = KK_mvr_rr_exp
    beta_exp  = KK_mvr_Nr_exp

    ! Calculate the coefficient A in the < N_r'R_vr' > covariance equation,
    ! where < N_r'R_vr' > = coef_A <N_r> + term_B, and:
    !
    ! coef_A
    ! = KK_mvr_coef
    !   ( exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
    !          + (1/2) sigma_rr_1_n^2 alpha^2
    !          + (1/2) sigma_Nr_1_n^2 ( beta^2 + 2 beta )
    !          + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
    !   + exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
    !          + (1/2) sigma_rr_2_n^2 alpha^2
    !          + (1/2) sigma_Nr_2_n^2 ( beta^2 + 2 beta )
    !          + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
    !   )
    !   - < R_vr >.
    covar_Nr_KK_mvr_coefA  &
    = KK_mvr_coef &
      * ( bivar_LL_covar_partial_Nr( mu_rr_1, mu_Nr_1, mu_rr_1_n, &
                                     mu_Nr_1_n, sigma_rr_1, sigma_Nr_1, &
                                     sigma_rr_1_n, sigma_Nr_1_n, &
                                     corr_rr_Nr_1_n, alpha_exp, beta_exp ) &
        + bivar_LL_covar_partial_Nr( mu_rr_2, mu_Nr_2, mu_rr_2_n, &
                                     mu_Nr_2_n, sigma_rr_2, sigma_Nr_2, &
                                     sigma_rr_2_n, sigma_Nr_2_n, &
                                     corr_rr_Nr_2_n, alpha_exp, beta_exp ) &
        ) &
      - KK_mean_vol_rad


    return

  end function covar_Nr_KK_mvr_coefA

  !=============================================================================
  function covar_Nr_KK_mvr_termB( Nr_1, Nr_2, mu_rr_1, mu_rr_2, mu_Nr_1, &
                                  mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                  mu_Nr_2_n, sigma_rr_1, sigma_rr_2,  &
                                  sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                                  sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                                  corr_rr_Nr_1_n, corr_rr_Nr_2_n, KK_mvr_coef, &
                                  mixt_frac )

    ! Description:
    ! This function partially calculates the covariance of N_r and KK mean
    ! volume radius of rain drops (R_vr), which can be written as < N_r'R_vr' >.
    !
    ! The equation for the covariance of N_r and KK rain drop mean volume radius
    ! is:
    !
    ! < N_r'R_vr' >
    ! = KK_mvr_coef
    !   ( a f_p(1)
    !     exp{ mu_rr_1_n alpha + mu_Nr_1_n ( beta + 1 )
    !          + (1/2) sigma_rr_1_n^2 alpha^2
    !          + (1/2) sigma_Nr_1_n^2 ( beta + 1 )^2
    !          + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
    !   + ( 1 - a ) f_p(2)
    !     exp{ mu_rr_2_n alpha + mu_Nr_2_n ( beta + 1 )
    !          + (1/2) sigma_rr_2_n^2 alpha^2
    !          + (1/2) sigma_Nr_2_n^2 ( beta + 1 )^2
    !          + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
    !   )
    !   - < R_vr > < N_r >;
    !
    ! where "a" is the mixture fraction (weight of the 1st PDF component), and
    ! f_p(1) and f_p(2) are the precipitation fractions in the 1st and 2nd PDF
    ! components, respectively.  The equation for < N_r'R_vr' > needs to be put
    ! in the form:  < N_r'R_vr' > = coef_A < N_r > + term_B.  Since < N_r > can
    ! be written as:
    !
    ! < N_r > = a f_p(1) exp{ mu_Nr_1_n + (1/2) sigma_Nr_1_n^2 }
    !           + ( 1 - a ) f_p(2) exp{ mu_Nr_2_n + (1/2) sigma_Nr_2_n^2 }
    !         = a f_p(1) mu_Nr_1 + ( 1 - a ) f_p(2) mu_Nr_2
    !         = a Nr_1 + ( 1 - a ) Nr_2;
    !
    ! the covariance < N_r'R_vr' > can be rewritten in terms of < N_r >:
    !
    ! < N_r'R_vr' >
    ! = ( KK_mvr_coef
    !     ( exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
    !            + (1/2) sigma_rr_1_n^2 alpha^2
    !            + (1/2) sigma_Nr_1_n^2 ( beta^2 + 2 beta )
    !            + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
    !     + exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
    !            + (1/2) sigma_rr_2_n^2 alpha^2
    !            + (1/2) sigma_Nr_2_n^2 ( beta^2 + 2 beta )
    !            + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
    !     )
    !     - < R_vr > ) < N_r >
    !   - KK_mvr_coef
    !     ( a Nr_1
    !       exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
    !            + (1/2) sigma_rr_2_n^2 alpha^2
    !            + (1/2) sigma_Nr_2_n^2 ( beta^2 + 2 beta )
    !            + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
    !     + ( 1 - a ) Nr_2
    !       exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
    !            + (1/2) sigma_rr_1_n^2 alpha^2
    !            + (1/2) sigma_Nr_1_n^2 ( beta^2 + 2 beta )
    !            + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
    !     ).
    !
    ! This function calculates term B, which is given by:
    !
    ! term_B
    ! = - KK_mvr_coef
    !     ( a Nr_1
    !       exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
    !            + (1/2) sigma_rr_2_n^2 alpha^2
    !            + (1/2) sigma_Nr_2_n^2 ( beta^2 + 2 beta )
    !            + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
    !     + ( 1 - a ) Nr_2
    !       exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
    !            + (1/2) sigma_rr_1_n^2 alpha^2
    !            + (1/2) sigma_Nr_1_n^2 ( beta^2 + 2 beta )
    !            + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
    !     ).
    !
    ! The covariance of sedimentation velocity (of rain drop concentration) and
    ! rain drop concentration, < V_Nr'N_r' > is given by the equation:
    !
    ! < V_Nr'N_r' > = 0.007 ( 10^6 microm/m ) < N_r'R_vr' >.
    !
    ! It is advantageous to write the < V_Nr'N_r' > term in terms of < N_r >
    ! because < V_Nr'N_r' > is found in the turbulent sedimentation term as part
    ! of the predictive equation for < N_r >.  This allows the turbulent
    ! sedimentation term to be handled semi-implicitly in the < N_r > predictive
    ! equation, increasing the numerical stability of the solution.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one    ! Constant(s)

    use parameters_KK, only: &
        KK_mvr_rr_exp, & ! Variable(s)
        KK_mvr_Nr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Nr_1,           & ! Mean rain drop concentration (1st PDF comp.)  [num/kg]
      Nr_2,           & ! Mean rain drop concentration (2nd PDF comp.)  [num/kg]
      mu_rr_1,        & ! Mean of rr (1st PDF component) in-precip (ip)  [kg/kg]
      mu_rr_2,        & ! Mean of rr (2nd PDF component) ip              [kg/kg]
      mu_Nr_1,        & ! Mean of Nr (1st PDF component) ip             [num/kg]
      mu_Nr_2,        & ! Mean of Nr (2nd PDF component) ip             [num/kg]
      mu_rr_1_n,      & ! Mean of ln rr (1st PDF component) ip       [ln(kg/kg)]
      mu_rr_2_n,      & ! Mean of ln rr (2nd PDF component) ip       [ln(kg/kg)]
      mu_Nr_1_n,      & ! Mean of ln Nr (1st PDF component) ip      [ln(num/kg)]
      mu_Nr_2_n,      & ! Mean of ln Nr (2nd PDF component) ip      [ln(num/kg)]
      sigma_rr_1,     & ! Standard deviation of rr (1st PDF comp.) ip    [kg/kg]
      sigma_rr_2,     & ! Standard deviation of rr (2nd PDF comp.) ip    [kg/kg]
      sigma_Nr_1,     & ! Standard deviation of Nr (1st PDF comp.) ip   [num/kg]
      sigma_Nr_2,     & ! Standard deviation of Nr (2nd PDF comp.) ip   [num/kg]
      sigma_rr_1_n,   & ! Standard deviation of ln rr (1st PDF component) ip [-]
      sigma_rr_2_n,   & ! Standard deviation of ln rr (2nd PDF component) ip [-]
      sigma_Nr_1_n,   & ! Standard deviation of ln Nr (1st PDF component) ip [-]
      sigma_Nr_2_n,   & ! Standard deviation of ln Nr (2nd PDF component) ip [-]
      corr_rr_Nr_1_n, & ! Correlation of ln rr & ln Nr (1st PDF comp.) ip    [-]
      corr_rr_Nr_2_n, & ! Correlation of ln rr & ln Nr (2nd PDF comp.) ip    [-]
      KK_mvr_coef,    & ! KK mean volume radius coefficient                  [m]
      mixt_frac         ! Mixture fraction                                   [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_Nr_KK_mvr_termB    ! Term in < N_r'R_vr' > equation      [(num/kg)m]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on r_r
      beta_exp     ! Exponent on N_r


    ! Values of the KK exponents.
    alpha_exp = KK_mvr_rr_exp
    beta_exp  = KK_mvr_Nr_exp

    ! Calculate term B in the < N_r'R_vr' > covariance equation, where
    ! < N_r'R_vr' > = coef_A <N_r> + term_B, and:
    !
    ! term_B
    ! = - KK_mvr_coef
    !     ( a Nr_1
    !       exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
    !            + (1/2) sigma_rr_2_n^2 alpha^2
    !            + (1/2) sigma_Nr_2_n^2 ( beta^2 + 2 beta )
    !            + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
    !     + ( 1 - a ) Nr_2
    !       exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
    !            + (1/2) sigma_rr_1_n^2 alpha^2
    !            + (1/2) sigma_Nr_1_n^2 ( beta^2 + 2 beta )
    !            + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
    !     ).
    covar_Nr_KK_mvr_termB  &
    = - KK_mvr_coef &
        * ( mixt_frac * Nr_1 &
            * bivar_LL_covar_partial_Nr( mu_rr_2, mu_Nr_2, mu_rr_2_n, &
                                         mu_Nr_2_n, sigma_rr_2, sigma_Nr_2, &
                                         sigma_rr_2_n, sigma_Nr_2_n, &
                                         corr_rr_Nr_2_n, alpha_exp, beta_exp ) &
          + ( one - mixt_frac ) * Nr_2 &
            * bivar_LL_covar_partial_Nr( mu_rr_1, mu_Nr_1, mu_rr_1_n, &
                                         mu_Nr_1_n, sigma_rr_1, sigma_Nr_1, &
                                         sigma_rr_1_n, sigma_Nr_1_n, &
                                         corr_rr_Nr_1_n, alpha_exp, beta_exp ) &
          )


    return

  end function covar_Nr_KK_mvr_termB

  !=============================================================================
  function covar_Nr_KK_mvr( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, &
                            mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1, &
                            sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                            sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                            corr_rr_Nr_1_n, corr_rr_Nr_2_n, Nrm, &
                            KK_mean_vol_rad, KK_mvr_coef, mixt_frac, &
                            precip_frac_1, precip_frac_2 )

    ! Description:
    ! This function calculates the covariance of N_r and KK mean volume radius
    ! of rain drops (R_vr), which can be written as < N_r'R_vr' >.
    !
    ! The equation for the covariance of N_r and KK rain drop mean volume radius
    ! is:
    !
    ! < N_r'R_vr' >
    ! = KK_mvr_coef
    !   ( a f_p(1)
    !     exp{ mu_rr_1_n alpha + mu_Nr_1_n ( beta + 1 )
    !          + (1/2) sigma_rr_1_n^2 alpha^2
    !          + (1/2) sigma_Nr_1_n^2 ( beta + 1 )^2
    !          + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
    !   + ( 1 - a ) f_p(2)
    !     exp{ mu_rr_2_n alpha + mu_Nr_2_n ( beta + 1 )
    !          + (1/2) sigma_rr_2_n^2 alpha^2
    !          + (1/2) sigma_Nr_2_n^2 ( beta + 1 )^2
    !          + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
    !   )
    !   - < R_vr > < N_r >;
    !
    ! where "a" is the mixture fraction (weight of the 1st PDF component), and
    ! f_p(1) and f_p(2) are the precipitation fractions in the 1st and 2nd PDF
    ! components, respectively.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one  ! Constant(s)

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
      mu_rr_1_n,       & ! Mean of ln rr (1st PDF component) ip      [ln(kg/kg)]
      mu_rr_2_n,       & ! Mean of ln rr (2nd PDF component) ip      [ln(kg/kg)]
      mu_Nr_1_n,       & ! Mean of ln Nr (1st PDF component) ip     [ln(num/kg)]
      mu_Nr_2_n,       & ! Mean of ln Nr (2nd PDF component) ip     [ln(num/kg)]
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
      Nrm,             & ! Mean rain drop concentration                 [num/kg]
      KK_mean_vol_rad, & ! KK mean volume radius of rain drops               [m]
      KK_mvr_coef,     & ! KK mean volume radius coefficient                 [m]
      mixt_frac,       & ! Mixture fraction                                  [-]
      precip_frac_1,   & ! Precipitation fraction (1st PDF component)        [-]
      precip_frac_2      ! Precipitation fraction (2nd PDF component)        [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_Nr_KK_mvr    ! Covar of Nr and KK rain drop mean vol rad [(num/kg)m]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on r_r
      beta_exp     ! Exponent on N_r


    ! Values of the KK exponents.
    alpha_exp = KK_mvr_rr_exp
    beta_exp  = KK_mvr_Nr_exp

    ! Calculate the covariance of N_r and KK mean volume radius of rain drops.
    covar_Nr_KK_mvr  &
    = KK_mvr_coef &
      * ( mixt_frac &
          * precip_frac_1 &
          * bivar_LL_mean_eq( mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, &
                              sigma_rr_1, sigma_Nr_1, sigma_rr_1_n, &
                              sigma_Nr_1_n, corr_rr_Nr_1_n, &
                              alpha_exp, beta_exp + one ) &
        + ( one - mixt_frac ) &
          * precip_frac_2 &
          * bivar_LL_mean_eq( mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, &
                              sigma_rr_2, sigma_Nr_2, sigma_rr_2_n, &
                              sigma_Nr_2_n, corr_rr_Nr_2_n, &
                              alpha_exp, beta_exp + one ) &
        ) &
      - Nrm * KK_mean_vol_rad


    return

  end function covar_Nr_KK_mvr

  !=============================================================================
  function bivar_LL_covar_partial_rr( mu_rr_i, mu_Nr_i, mu_rr_i_n, &
                                      mu_Nr_i_n, sigma_rr_i, sigma_Nr_i, &
                                      sigma_rr_i_n, sigma_Nr_i_n, &
                                      corr_rr_Nr_i_n, alpha_exp_in, beta_exp_in )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use PDF_integrals_means, only: &
        bivar_LL_mean_const_x1,  & ! Procedure(s)
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
      mu_rr_i,        & ! Mean of rr (ith PDF component) in-precip (ip)  [kg/kg]
      mu_Nr_i,        & ! Mean of Nr (ith PDF component) ip             [num/kg]
      mu_rr_i_n,      & ! Mean of ln rr (ith PDF component) ip       [ln(kg/kg)]
      mu_Nr_i_n,      & ! Mean of ln Nr (ith PDF component) ip      [ln(num/kg)]
      sigma_rr_i,     & ! Standard deviation of rr (ith PDF comp.) ip    [kg/kg]
      sigma_Nr_i,     & ! Standard deviation of Nr (ith PDF comp.) ip   [num/kg]
      sigma_rr_i_n,   & ! Standard deviation of ln rr (ith PDF component) ip [-]
      sigma_Nr_i_n,   & ! Standard deviation of ln Nr (ith PDF component) ip [-]
      corr_rr_Nr_i_n    ! Correlation of ln rr & ln Nr (ith PDF comp.) ip    [-]

    real( kind = core_rknd ), intent(in) :: &
      alpha_exp_in,  & ! Exponent alpha, corresponding to rr                 [-]
      beta_exp_in      ! Exponent beta, corresponding to Nr                  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      bivar_LL_covar_partial_rr

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
      rho_x1x2_n    ! Correlation of ln x1 & ln x2 (ith PDF component)      [-]

    real( kind = dp ) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp      ! Exponent beta, corresponding to x2                    [-]

    real( kind = dp ) :: &
      x1_tol, & ! Tolerance value of x1                                     [-]
      x2_tol    ! Tolerance value of x2                                     [-]


    ! Means for the ith PDF component.
    if ( alpha_exp_in >= zero ) then
       mu_x1 = real(mu_rr_i, kind=dp)
    else ! exponent alpha < 0
       mu_x1 = real(max( mu_rr_i, rr_tol ), kind=dp)
    endif
    if ( beta_exp_in >= zero ) then
       mu_x2 = real(mu_Nr_i, kind=dp)
    else ! exponent beta < 0
       mu_x2 = real(max( mu_Nr_i, Nr_tol ), kind=dp)
    endif
    mu_x1_n = real(mu_rr_i_n, kind=dp)
    mu_x2_n = real(mu_Nr_i_n, kind=dp)

    ! Standard deviations for the ith PDF component.
    sigma_x1   = real(sigma_rr_i, kind=dp)
    sigma_x2   = real(sigma_Nr_i, kind=dp)
    sigma_x1_n = real(sigma_rr_i_n, kind=dp)
    sigma_x2_n = real(sigma_Nr_i_n, kind=dp)

    ! Correlations for the ith PDF component.
    rho_x1x2_n = real(corr_rr_Nr_i_n, kind=dp)

    ! Exponents.
    alpha_exp = real(alpha_exp_in, kind=dp)
    beta_exp  = real(beta_exp_in, kind=dp)

    ! Tolerance values.
    ! When the standard deviation of a variable is below the tolerance values,
    ! it is considered to be zero, and the variable is considered to have a
    ! constant value.
    x1_tol = real(rr_tol, kind=dp)
    x2_tol = real(Nr_tol, kind=dp)


    ! Calculate (partially) the covariance of the bivariate lognormal equation.
    if ( sigma_x1 <= x1_tol .and. sigma_x2 <= x2_tol ) then

       ! The ith PDF component variance of both r_r and N_r is 0.
       bivar_LL_covar_partial_rr  &
       = real( bivar_LL_mean_const_all( mu_x1, mu_x2, alpha_exp, beta_exp ), &
               kind = core_rknd )


    elseif ( sigma_x1 <= x1_tol ) then

       ! The ith PDF component variance of r_r is 0.
       bivar_LL_covar_partial_rr  &
       = real( bivar_LL_mean_const_x1( mu_x1, mu_x2_n, sigma_x2_n, &
                                       alpha_exp, beta_exp ), &
               kind = core_rknd )


    elseif ( sigma_x2 <= x2_tol ) then

       ! The ith PDF component variance of N_r is 0.
       bivar_LL_covar_partial_rr  &
       = real( bivar_LL_covar_const_x2_partial( mu_x1_n, mu_x2, sigma_x1_n, &
                                                alpha_exp, beta_exp ), &
               kind = core_rknd )


    else  ! sigma_x1 and sigma_x2 > 0

       ! All fields vary in the ith PDF component.
       bivar_LL_covar_partial_rr  &
       = real( bivar_LL_covar_partial( mu_x1_n, mu_x2_n, sigma_x1_n, &
                                       sigma_x2_n, rho_x1x2_n, &
                                       alpha_exp, beta_exp ), &
               kind = core_rknd )


    endif


    return

  end function bivar_LL_covar_partial_rr

  !=============================================================================
  function bivar_LL_covar_partial_Nr( mu_rr_i, mu_Nr_i, mu_rr_i_n, &
                                      mu_Nr_i_n, sigma_rr_i, sigma_Nr_i, &
                                      sigma_rr_i_n, sigma_Nr_i_n, &
                                      corr_rr_Nr_i_n, alpha_exp_in, beta_exp_in )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use PDF_integrals_means, only: &
        bivar_LL_mean_const_x1,  & ! Procedure(s)
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
      mu_rr_i,        & ! Mean of rr (ith PDF component) in-precip (ip)  [kg/kg]
      mu_Nr_i,        & ! Mean of Nr (ith PDF component) ip             [num/kg]
      mu_rr_i_n,      & ! Mean of ln rr (ith PDF component) ip       [ln(kg/kg)]
      mu_Nr_i_n,      & ! Mean of ln Nr (ith PDF component) ip      [ln(num/kg)]
      sigma_rr_i,     & ! Standard deviation of rr (ith PDF comp.) ip    [kg/kg]
      sigma_Nr_i,     & ! Standard deviation of Nr (ith PDF comp.) ip   [num/kg]
      sigma_rr_i_n,   & ! Standard deviation of ln rr (ith PDF component) ip [-]
      sigma_Nr_i_n,   & ! Standard deviation of ln Nr (ith PDF component) ip [-]
      corr_rr_Nr_i_n    ! Correlation of ln rr & ln Nr (ith PDF comp.) ip    [-]

    real( kind = core_rknd ), intent(in) :: &
      alpha_exp_in,  & ! Exponent alpha, corresponding to rr                 [-]
      beta_exp_in      ! Exponent beta, corresponding to Nr                  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      bivar_LL_covar_partial_Nr

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
      rho_x1x2_n    ! Correlation of ln x1 & ln x2 (ith PDF component)      [-]

    real( kind = dp ) :: &
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp      ! Exponent beta, corresponding to x2                    [-]

    real( kind = dp ) :: &
      x1_tol, & ! Tolerance value of x1                                     [-]
      x2_tol    ! Tolerance value of x2                                     [-]


    ! Normally, x1 is used for r_r and x2 is used for N_r.  Here, x1 will be
    ! used for N_r and x2 will be used for r_r.  This means that the exponent
    ! on N_r, beta_exp_in, will become alpha_exp, and the exponent on r_r,
    ! alpha_exp_in, will become beta_exp.

    ! Means for the ith PDF component.
    if ( beta_exp_in >= zero ) then
       mu_x1 = real(mu_Nr_i, kind=dp)
    else ! exponent beta < 0
       mu_x1 = real(max( mu_Nr_i, Nr_tol ), kind=dp)
    endif
    if ( alpha_exp_in >= zero ) then
       mu_x2 = real(mu_rr_i, kind=dp)
    else ! exponent alpha < 0
       mu_x2 = real(max( mu_rr_i, rr_tol ), kind=dp)
    endif
    mu_x1_n = real(mu_Nr_i_n, kind=dp)
    mu_x2_n = real(mu_rr_i_n, kind=dp)

    ! Standard deviations for the ith PDF component.
    sigma_x1   =  real(sigma_Nr_i, kind=dp)
    sigma_x2   = real(sigma_rr_i, kind=dp)
    sigma_x1_n = real(sigma_Nr_i_n, kind=dp)
    sigma_x2_n = real(sigma_rr_i_n, kind=dp)

    ! Correlations for the ith PDF component.
    rho_x1x2_n = real(corr_rr_Nr_i_n, kind=dp)

    ! Exponents.
    alpha_exp = real(beta_exp_in, kind=dp)
    beta_exp  = real(alpha_exp_in, kind=dp)

    ! Tolerance values.
    ! When the standard deviation of a variable is below the tolerance values,
    ! it is considered to be zero, and the variable is considered to have a
    ! constant value.
    x1_tol = real(Nr_tol, kind=dp)
    x2_tol = real(rr_tol, kind=dp)


    ! Calculate (partially) the covariance of the bivariate lognormal equation.
    if ( sigma_x1 <= x1_tol .and. sigma_x2 <= x2_tol ) then

       ! The ith PDF component variance of both r_r and N_r is 0.
       bivar_LL_covar_partial_Nr  &
       = real( bivar_LL_mean_const_all( mu_x1, mu_x2, alpha_exp, beta_exp ), &
               kind = core_rknd )


    elseif ( sigma_x1 <= x1_tol ) then

       ! The ith PDF component variance of N_r is 0.
       bivar_LL_covar_partial_Nr  &
       = real( bivar_LL_mean_const_x1( mu_x1, mu_x2_n, sigma_x2_n, &
                                       alpha_exp, beta_exp ), &
               kind = core_rknd )


    elseif ( sigma_x2 <= x2_tol ) then

       ! The ith PDF component variance of r_r is 0.
       bivar_LL_covar_partial_Nr  &
       = real( bivar_LL_covar_const_x2_partial( mu_x1_n, mu_x2, sigma_x1_n, &
                                                alpha_exp, beta_exp ), &
               kind = core_rknd )


    else  ! sigma_x1 and sigma_x2 > 0

       ! All fields vary in the ith PDF component.
       bivar_LL_covar_partial_Nr  &
       = real( bivar_LL_covar_partial( mu_x1_n, mu_x2_n, sigma_x1_n, &
                                       sigma_x2_n, rho_x1x2_n, &
                                       alpha_exp, beta_exp ), &
               kind = core_rknd )


    endif


    return

  end function bivar_LL_covar_partial_Nr

  !=============================================================================
  function bivar_LL_covar_partial( mu_x1_n, mu_x2_n, sigma_x1_n, &
                                   sigma_x2_n, rho_x1x2_n, &
                                   alpha_exp, beta_exp )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one_half_dp, & ! Constant(s)
        one_dp, &
        two_dp

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1_n,    & ! Mean of ln x1 (ith PDF component)                     [-]
      mu_x2_n,    & ! Mean of ln x2 (ith PDF component)                     [-]
      sigma_x1_n, & ! Standard deviation of ln x1 (ith PDF component)       [-]
      sigma_x2_n, & ! Standard deviation of ln x2 (ith PDF component)       [-]
      rho_x1x2_n, & ! Correlation of ln x1 & ln x2 (ith PDF component)      [-]
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp      ! Exponent beta, corresponding to x2                    [-]

    ! Return Variable
    real( kind = dp ) ::  &
      bivar_LL_covar_partial


    bivar_LL_covar_partial  &
    = exp( mu_x1_n * alpha_exp + mu_x2_n * beta_exp  &
           + one_half_dp * sigma_x1_n**2 &
                         * ( alpha_exp**2 + two_dp * alpha_exp ) &
           + one_half_dp * sigma_x2_n**2 * beta_exp**2  &
           + rho_x1x2_n * sigma_x1_n * ( alpha_exp + one_dp ) &
                        * sigma_x2_n * beta_exp )


    return

  end function bivar_LL_covar_partial

  !=============================================================================
  function bivar_LL_covar_const_x2_partial( mu_x1_n, mu_x2, sigma_x1_n, &
                                            alpha_exp, beta_exp )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one_half_dp, & ! Constant(s)
        two_dp

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      mu_x1_n,    & ! Mean of ln x1 (ith PDF component)                     [-]
      mu_x2,      & ! Mean of x2 (ith PDF component)                        [-]
      sigma_x1_n, & ! Standard deviation of ln x1 (ith PDF component)       [-]
      alpha_exp,  & ! Exponent alpha, corresponding to x1                   [-]
      beta_exp      ! Exponent beta, corresponding to x2                    [-]

    ! Return Variable
    real( kind = dp ) ::  &
      bivar_LL_covar_const_x2_partial


    bivar_LL_covar_const_x2_partial &
    = mu_x2**beta_exp &
      * exp( mu_x1_n * alpha_exp &
             + one_half_dp * sigma_x1_n**2 &
                           * ( alpha_exp**2 + two_dp * alpha_exp ) )


    return

  end function bivar_LL_covar_const_x2_partial

!===============================================================================

end module KK_upscaled_turbulent_sed
