! $Id$
!===============================================================================
module KK_upscaled_turbulent_sed

  implicit none

  private  ! Default scope

  public :: KK_sed_vel_covars

  private :: covar_rr_KK_mvr, &
             covar_Nr_KK_mvr

  contains

  !=============================================================================
  subroutine KK_sed_vel_covars( rrainm, Nrm, KK_mean_vol_rad, &
                                mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                                sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                                sigma_Nr_2_n, corr_rrNr_1_n, corr_rrNr_2_n, &
                                KK_mvr_coef, mixt_frac, precip_frac_1, &
                                precip_frac_2, level, l_stats_samp, &
                                Vrrprrp, VNrpNrp )

    ! Description:
    ! 

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        micron_per_m, & ! Constant(s)
        rr_tol, &
        Nr_tol, &
        zero

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use stats_type, only: &
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only: & 
        irr_KK_mvr_covar_zt, & ! Variable(s)
        iNr_KK_mvr_covar_zt, &
        zt

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rrainm,          & ! Mean rain water mixing ratio                  [kg/kg]
      Nrm,             & ! Mean rain drop concentration                 [num/kg]
      KK_mean_vol_rad, & ! KK mean volume radius of rain drops               [m]
      mu_rr_1_n,       & ! Mean of ln rr (1st PDF component) in-precip (ip)  [-]
      mu_rr_2_n,       & ! Mean of ln rr (2nd PDF component) ip              [-]
      mu_Nr_1_n,       & ! Mean of ln Nr (1st PDF component) ip              [-]
      mu_Nr_2_n,       & ! Mean of ln Nr (2nd PDF component) ip              [-]
      sigma_rr_1_n,    & ! Standard deviation of ln rr (1st PDF comp.) ip    [-]
      sigma_rr_2_n,    & ! Standard deviation of ln rr (2nd PDF comp.) ip    [-]
      sigma_Nr_1_n,    & ! Standard deviation of ln Nr (1st PDF comp.) ip    [-]
      sigma_Nr_2_n,    & ! Standard deviation of ln Nr (2nd PDF comp.) ip    [-]
      corr_rrNr_1_n,   & ! Corr. betw. ln rr & ln Nr (1st PDF comp.) ip      [-]
      corr_rrNr_2_n,   & ! Corr. betw. ln rr & ln Nr (2nd PDF comp.) ip      [-]
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
      Vrrprrp, & ! Covariance between V_rr and r_r, <V_rr'r_r'>   [(m/s)(kg/kg)]
      VNrpNrp    ! Covariance between V_Nr and N_r, <V_Nr'N_r'>  [(m/s)(num/kg)]

    ! Local Variables
    real( kind = core_rknd ) :: &
      rr_KK_mvr_covar, & ! Covariance of r_r and KK mean vol rad   [(kg/kg)m]
      Nr_KK_mvr_covar    ! Covariance of N_r and KK mean vol rad   [(num/kg)m]


    ! Calculate the covariance between rain drop mean volume radius and r_r,
    ! < R_vr'r_r' >.
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then

       rr_KK_mvr_covar  &
       = covar_rr_KK_mvr( mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                          sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                          sigma_Nr_2_n, corr_rrNr_1_n, corr_rrNr_2_n, &
                          rrainm, KK_mean_vol_rad, KK_mvr_coef, &
                          mixt_frac, precip_frac_1, precip_frac_2 )

    else  ! r_r or N_r = 0.

       rr_KK_mvr_covar = zero

    endif


    ! Calculate the covariance between rain drop mean volume radius and N_r,
    ! < R_vr'N_r' >.
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then

       Nr_KK_mvr_covar  &
       = covar_Nr_KK_mvr( mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                          sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                          sigma_Nr_2_n, corr_rrNr_1_n, corr_rrNr_2_n, &
                          Nrm, KK_mean_vol_rad, KK_mvr_coef, &
                          mixt_frac, precip_frac_1, precip_frac_2 )

    else  ! r_r or N_r = 0.

       Nr_KK_mvr_covar = zero

    endif


    ! Covariance between V_rr and r_r, < V_rr'r_r' >.
    Vrrprrp = - 0.012_core_rknd * micron_per_m * rr_KK_mvr_covar

    ! Covariance between V_Nr and N_r, < V_Nr'N_r' >.
    VNrpNrp = - 0.007_core_rknd * micron_per_m * Nr_KK_mvr_covar


    ! Statistics
    if ( l_stats_samp ) then

       ! Covariance between r_r and KK rain drop mean volume radius.
       if ( irr_KK_mvr_covar_zt > 0 ) then
          call stat_update_var_pt( irr_KK_mvr_covar_zt, level, &
                                   rr_KK_mvr_covar, zt )
       endif

       ! Covariance between N_r and KK rain drop mean volume radius.
       if ( iNr_KK_mvr_covar_zt > 0 ) then
          call stat_update_var_pt( iNr_KK_mvr_covar_zt, level, &
                                   Nr_KK_mvr_covar, zt )
       endif

    endif ! l_stats_samp


    return

  end subroutine KK_sed_vel_covars

  !=============================================================================
  function covar_rr_KK_mvr( mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                            sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                            sigma_Nr_2_n, corr_rrNr_1_n, corr_rrNr_2_n, &
                            rrm, KK_mean_vol_rad, KK_mvr_coef, &
                            mixt_frac, precip_frac_1, precip_frac_2 )

    ! Description:
    ! This function calculates the covariance between r_r and KK mean
    ! volume
    ! radius of rain drops (R_vr), which can be written as < r_r'R_vr'
    ! >.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one

    use KK_upscaled_means, only:  &
        bivar_LL_mean_eq

    use parameters_microphys, only: &
        KK_mvr_rr_exp, & ! Variable(s)
        KK_mvr_Nr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_rr_1_n,       & ! Mean of ln rr (1st PDF component) in-precip (ip)  [-]
      mu_rr_2_n,       & ! Mean of ln rr (2nd PDF component) ip              [-]
      mu_Nr_1_n,       & ! Mean of ln Nr (1st PDF component) ip              [-]
      mu_Nr_2_n,       & ! Mean of ln Nr (2nd PDF component) ip              [-]
      sigma_rr_1_n,    & ! Standard deviation of ln rr (1st PDF comp.) ip    [-]
      sigma_rr_2_n,    & ! Standard deviation of ln rr (2nd PDF comp.) ip    [-]
      sigma_Nr_1_n,    & ! Standard deviation of ln Nr (1st PDF comp.) ip    [-]
      sigma_Nr_2_n,    & ! Standard deviation of ln Nr (2nd PDF comp.) ip    [-]
      corr_rrNr_1_n,   & ! Corr. betw. ln rr & ln Nr (1st PDF comp.) ip      [-]
      corr_rrNr_2_n,   & ! Corr. betw. ln rr & ln Nr (2nd PDF comp.) ip      [-]
      rrm,             & ! Mean rain water mixing ratio                  [kg/kg]
      KK_mean_vol_rad, & ! KK mean volume radius of rain drops               [m]
      KK_mvr_coef,     & ! KK mean volume radius coefficient                 [m]
      mixt_frac,       & ! Mixture fraction                                  [-]
      precip_frac_1,   & ! Precipitation fraction (1st PDF component)        [-]
      precip_frac_2      ! Precipitation fraction (2nd PDF component)        [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_rr_KK_mvr  ! Covar of rr and KK rain drop mean vol rad    [(kg/kg)m]

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
          * bivar_LL_mean_eq( mu_rr_1_n, mu_Nr_1_n, sigma_rr_1_n, &
                              sigma_Nr_1_n, corr_rrNr_1_n, &
                              alpha_exp + one, beta_exp ) &
        + ( one - mixt_frac ) &
          * precip_frac_2 &
          * bivar_LL_mean_eq( mu_rr_2_n, mu_Nr_2_n, sigma_rr_2_n, &
                              sigma_Nr_2_n, corr_rrNr_2_n, &
                              alpha_exp + one, beta_exp ) &
        ) &
      - rrm * KK_mean_vol_rad


    return

  end function covar_rr_KK_mvr

  !=============================================================================
  function covar_Nr_KK_mvr( mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                            sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                            sigma_Nr_2_n, corr_rrNr_1_n, corr_rrNr_2_n, &
                            Nrm, KK_mean_vol_rad, KK_mvr_coef, &
                            mixt_frac, precip_frac_1, precip_frac_2 )

    ! Description:
    ! This function calculates the covariance between N_r and KK mean volume
    ! radius of rain drops (R_vr), which can be written as < N_r'R_vr' >.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one

    use KK_upscaled_means, only:  &
        bivar_LL_mean_eq

    use parameters_microphys, only: &
        KK_mvr_rr_exp, & ! Variable(s)
        KK_mvr_Nr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_rr_1_n,       & ! Mean of ln rr (1st PDF component) in-precip (ip)  [-]
      mu_rr_2_n,       & ! Mean of ln rr (2nd PDF component) ip              [-]
      mu_Nr_1_n,       & ! Mean of ln Nr (1st PDF component) ip              [-]
      mu_Nr_2_n,       & ! Mean of ln Nr (2nd PDF component) ip              [-]
      sigma_rr_1_n,    & ! Standard deviation of ln rr (1st PDF comp.) ip    [-]
      sigma_rr_2_n,    & ! Standard deviation of ln rr (2nd PDF comp.) ip    [-]
      sigma_Nr_1_n,    & ! Standard deviation of ln Nr (1st PDF comp.) ip    [-]
      sigma_Nr_2_n,    & ! Standard deviation of ln Nr (2nd PDF comp.) ip    [-]
      corr_rrNr_1_n,   & ! Corr. betw. ln rr & ln Nr (1st PDF comp.) ip      [-]
      corr_rrNr_2_n,   & ! Corr. betw. ln rr & ln Nr (2nd PDF comp.) ip      [-]
      Nrm,             & ! Mean rain drop concentration                 [num/kg]
      KK_mean_vol_rad, & ! KK mean volume radius of rain drops               [m]
      KK_mvr_coef,     & ! KK mean volume radius coefficient                 [m]
      mixt_frac,       & ! Mixture fraction                                  [-]
      precip_frac_1,   & ! Precipitation fraction (1st PDF component)        [-]
      precip_frac_2      ! Precipitation fraction (2nd PDF component)        [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      covar_Nr_KK_mvr  ! Covar of Nr and KK rain drop mean vol rad   [(num/kg)m]

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
          * bivar_LL_mean_eq( mu_rr_1_n, mu_Nr_1_n, sigma_rr_1_n, &
                              sigma_Nr_1_n, corr_rrNr_1_n, &
                              alpha_exp, beta_exp + one ) &
        + ( one - mixt_frac ) &
          * precip_frac_2 &
          * bivar_LL_mean_eq( mu_rr_2_n, mu_Nr_2_n, sigma_rr_2_n, &
                              sigma_Nr_2_n, corr_rrNr_2_n, &
                              alpha_exp, beta_exp + one ) &
        ) &
      - Nrm * KK_mean_vol_rad


    return

  end function covar_Nr_KK_mvr

!===============================================================================

end module KK_upscaled_turbulent_sed
