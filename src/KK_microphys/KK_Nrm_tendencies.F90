! $Id$
!===============================================================================
module KK_Nrm_tendencies

  implicit none

  private

  public :: KK_Nrm_evap_upscaled_mean, &
            KK_Nrm_evap_local_mean, &
            KK_Nrm_auto_mean

  contains

  !=============================================================================
  function KK_Nrm_evap_upscaled_mean( mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, &
                                      mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, &
                                      mu_Nr_1_n, mu_Nr_2_n, sigma_chi_1, &
                                      sigma_chi_2, sigma_rr_1, sigma_rr_2, &
                                      sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                                      sigma_rr_2_n, sigma_Nr_1_n, &
                                      sigma_Nr_2_n, corr_chi_rr_1_n, &
                                      corr_chi_rr_2_n, corr_chi_Nr_1_n, &
                                      corr_chi_Nr_2_n, corr_rr_Nr_1_n, &
                                      corr_rr_Nr_2_n, KK_evap_coef, mixt_frac, &
                                      precip_frac_1, precip_frac_2, dt )

    ! Description:
    ! This function calculates the upscaled evaporation rate of rain drop
    ! concentration.  The KK equation for evaporation of rain drop concentration
    ! relates the change of rain drop concentration due to evaporation to the
    ! change of rain water mixing ratio due to evaporation according to:
    !
    ! ( Delta N_rV / N_rV )|_evap = ( ( Delta r_r / r_r )|_evap )^nu;
    !
    ! where N_rV is rain drop concentration per unit volume (units of num/m^3),
    ! r_r is rain water mixing ratio, and nu is a tuned parameter. However, KK
    ! suggest nu = 1 is reasonable choice.  This equation can be written in
    ! terms of rain drop concentration (per unit mass), N_r (in num/kg):
    !
    ! ( Delta N_r / N_r )|_evap = ( ( Delta r_r / r_r )|_evap )^nu.
    !
    ! The change in r_r due to evaporation ( Delta r_r ) is given by the
    ! equation:
    !
    ! ( Delta r_r )|_evap = INT(t_i:t_f) (dr_r/dt)|_evap dt;
    !
    ! where (dr_r/dt)|_evap is the r_r evaporation rate, t_i is the elapsed time
    ! at the start of the time step, and t_f is the elapsed time when the time
    ! step finishes.  Likewise, for N_r:
    !
    ! ( Delta N_r )|_evap = INT(t_i:t_f) (dN_r/dt)|_evap dt;
    !
    ! where (dN_r/dt)|_evap is the N_r evaporation rate.  Both of these integral
    ! equations are approximated as:
    !
    ! ( Delta r_r )|_evap = (dr_r/dt)|_evap Delta t; and
    !
    ! ( Delta N_r )|_evap = (dN_r/dt)|_evap Delta t;
    !
    ! where Delta t = t_f - t_i, which is the duration of one model time step.
    !
    ! The equation for the change in N_r due to evaporation becomes:
    !
    ! (dN_r/dt)|_evap Delta t / N_r = ( (dr_r/dt)|_evap Delta t / r_r )^nu.
    !
    ! This can be rewritten in terms of N_r evaporation rate:
    !
    ! (dN_r/dt)|_evap = Delta_t^(nu-1) ( N_r / r_r^nu ) ( (dr_r/dt)|_evap )^nu.
    !
    ! The equation for r_r evaporation rate is:
    !
    ! (dr_r/dt)|_evap = KK_evap_coef chi^alpha r_r^beta N_r^gamma;
    !
    ! which is substituted into the (dN_r/dt)|_evap equation, which in turn
    ! becomes:
    !
    ! (dN_r/dt)|_evap
    ! = Delta_t^(nu-1) ( N_r / r_r^nu )
    !   ( KK_evap_coef chi^alpha r_r^beta N_r^gamma )^nu;
    !
    ! and can be rewritten as:
    !
    ! (dN_r/dt)|_evap
    ! = Delta_t^(nu-1) KK_evap_coef^nu
    !   chi^(alpha nu) r_r^((beta - 1) nu) N_r^(gamma nu + 1).
    !
    ! The mean value of < N_r > evaporation rate, < (dN_r/dt)|_evap >, is found
    ! by integrating over a trivariate PDF of chi, r_r, and N_r.
    !
    ! When nu is set to 1, as suggested, the equation simply becomes:
    !
    ! (dN_r/dt)|_evap = KK_evap_coef chi^alpha r_r^(beta - 1) N_r^(gamma + 1).

    ! References:
    ! Eq. (C35) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S35) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! A new subgrid-scale representation of hydrometeor fields using a
    ! multivariate PDF.  Geosci. Model Dev., 9, 6,
    ! doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one  ! Constant(s)

    use KK_upscaled_means, only: &
        trivar_NLL_mean_eq  ! Procedure(s)

    use parameters_KK, only: &
        KK_evap_Supersat_exp, & ! Variable(s)
        KK_evap_rr_exp,       &
        KK_evap_Nr_exp,       &
        KK_Nrm_evap_nu

    use clubb_precision, only: &
        core_rknd

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_chi_1,        & ! Mean of chi (old s) (1st PDF component)       [kg/kg]
      mu_chi_2,        & ! Mean of chi (old s) (2nd PDF component)       [kg/kg]
      mu_rr_1,         & ! Mean of rr (1st PDF component) in-precip (ip) [kg/kg]
      mu_rr_2,         & ! Mean of rr (2nd PDF component) ip             [kg/kg]
      mu_Nr_1,         & ! Mean of Nr (1st PDF component) ip            [num/kg]
      mu_Nr_2,         & ! Mean of Nr (2nd PDF component) ip            [num/kg]
      mu_rr_1_n,       & ! Mean of ln rr (1st PDF component) ip      [ln(kg/kg)]
      mu_rr_2_n,       & ! Mean of ln rr (2nd PDF component) ip      [ln(kg/kg)]
      mu_Nr_1_n,       & ! Mean of ln Nr (1st PDF component) ip     [ln(num/kg)]
      mu_Nr_2_n,       & ! Mean of ln Nr (2nd PDF component) ip     [ln(num/kg)]
      sigma_chi_1,     & ! Standard deviation of chi (1st PDF component) [kg/kg]
      sigma_chi_2,     & ! Standard deviation of chi (2nd PDF component) [kg/kg]
      sigma_rr_1,      & ! Standard deviation of rr (1st PDF comp.) ip   [kg/kg]
      sigma_rr_2,      & ! Standard deviation of rr (2nd PDF comp.) ip   [kg/kg]
      sigma_Nr_1,      & ! Standard deviation of Nr (1st PDF comp.) ip  [num/kg]
      sigma_Nr_2,      & ! Standard deviation of Nr (2nd PDF comp.) ip  [num/kg]
      sigma_rr_1_n,    & ! Standard deviation of ln rr (1st PDF comp.) ip    [-]
      sigma_rr_2_n,    & ! Standard deviation of ln rr (2nd PDF comp.) ip    [-]
      sigma_Nr_1_n,    & ! Standard deviation of ln Nr (1st PDF comp.) ip    [-]
      sigma_Nr_2_n,    & ! Standard deviation of ln Nr (2nd PDF comp.) ip    [-]
      corr_chi_rr_1_n, & ! Correlation of chi and ln rr (1st PDF comp.) ip   [-]
      corr_chi_rr_2_n, & ! Correlation of chi and ln rr (2nd PDF comp.) ip   [-]
      corr_chi_Nr_1_n, & ! Correlation of chi and ln Nr (1st PDF comp.) ip   [-]
      corr_chi_Nr_2_n, & ! Correlation of chi and ln Nr (2nd PDF comp.) ip   [-]
      corr_rr_Nr_1_n,  & ! Correlation of ln rr & ln Nr (1st PDF comp.) ip   [-]
      corr_rr_Nr_2_n,  & ! Correlation of ln rr & ln Nr (2nd PDF comp.) ip   [-]
      KK_evap_coef,    & ! KK evaporation coefficient                [(kg/kg)/s]
      mixt_frac,       & ! Mixture fraction                                  [-]
      precip_frac_1,   & ! Precipitation fraction (1st PDF component)        [-]
      precip_frac_2      ! Precipitation fraction (2nd PDF component)        [-]

    real( kind = core_rknd ), intent(in) :: &
      dt                 ! Model time step duration                          [s]

    ! Return Variables
    real( kind = core_rknd ) :: &
      KK_Nrm_evap_upscaled_mean  ! Mean KK <N_r> evaporation rate   [(num/kg)/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on chi                                         [-]
      beta_exp,  & ! Exponent on r_r                                         [-]
      gamma_exp    ! Exponent on N_r                                         [-]


    ! Values of the KK exponents.
    alpha_exp = KK_evap_Supersat_exp * KK_Nrm_evap_nu
    beta_exp  = ( KK_evap_rr_exp - one ) * KK_Nrm_evap_nu
    gamma_exp = KK_evap_Nr_exp * KK_Nrm_evap_nu + one

    ! Calculate the mean KK < N_r > evaporation tendency.
    KK_Nrm_evap_upscaled_mean &
    = dt**( KK_Nrm_evap_nu - one ) &
      * KK_evap_coef**KK_Nrm_evap_nu &
      * ( mixt_frac &
          * precip_frac_1 &
          * trivar_NLL_mean_eq( mu_chi_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, &
                                mu_Nr_1_n, sigma_chi_1, sigma_rr_1, &
                                sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, &
                                corr_chi_rr_1_n, corr_chi_Nr_1_n, &
                                corr_rr_Nr_1_n, &
                                alpha_exp, beta_exp, gamma_exp ) &
        + ( one - mixt_frac ) &
          * precip_frac_2 &
          * trivar_NLL_mean_eq( mu_chi_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, &
                                mu_Nr_2_n, sigma_chi_2, sigma_rr_2, &
                                sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, &
                                corr_chi_rr_2_n, corr_chi_Nr_2_n, &
                                corr_rr_Nr_2_n, &
                                alpha_exp, beta_exp, gamma_exp ) &
        )


    return

  end function KK_Nrm_evap_upscaled_mean

  !=============================================================================
  function KK_Nrm_evap_local_mean( KK_rrm_evap_tndcy, Nrm, rrm, dt )

    ! Description:
    ! This function calculates the local evaporation rate of rain drop
    ! concentration.  The KK equation for evaporation of rain drop concentration
    ! relates the change of rain drop concentration due to evaporation to the
    ! change of rain water mixing ratio due to evaporation according to:
    !
    ! ( Delta N_rV / N_rV )|_evap = ( ( Delta r_r / r_r )|_evap )^nu;
    !
    ! where N_rV is rain drop concentration per unit volume (units of num/m^3),
    ! r_r is rain water mixing ratio, and nu is a tuned parameter. However, KK
    ! suggest nu = 1 is reasonable choice.  This equation can be written in
    ! terms of rain drop concentration (per unit mass), N_r (in num/kg):
    !
    ! ( Delta N_r / N_r )|_evap = ( ( Delta r_r / r_r )|_evap )^nu.
    !
    ! The change in r_r due to evaporation ( Delta r_r ) is given by the
    ! equation:
    !
    ! ( Delta r_r )|_evap = INT(t_i:t_f) (dr_r/dt)|_evap dt;
    !
    ! where (dr_r/dt)|_evap is the r_r evaporation rate, t_i is the elapsed time
    ! at the start of the time step, and t_f is the elapsed time when the time
    ! step finishes.  Likewise, for N_r:
    !
    ! ( Delta N_r )|_evap = INT(t_i:t_f) (dN_r/dt)|_evap dt;
    !
    ! where (dN_r/dt)|_evap is the N_r evaporation rate.  Both of these integral
    ! equations are approximated as:
    !
    ! ( Delta r_r )|_evap = (dr_r/dt)|_evap Delta t; and
    !
    ! ( Delta N_r )|_evap = (dN_r/dt)|_evap Delta t;
    !
    ! where Delta t = t_f - t_i, which is the duration of one model time step.
    !
    ! The equation for the change in N_r due to evaporation becomes:
    !
    ! (dN_r/dt)|_evap Delta t / N_r = ( (dr_r/dt)|_evap Delta t / r_r )^nu.
    !
    ! This can be rewritten in terms of N_r evaporation rate:
    !
    ! (dN_r/dt)|_evap = Delta_t^(nu-1) ( N_r / r_r^nu ) ( (dr_r/dt)|_evap )^nu.
    ! 
    ! When nu is set to 1, as suggested, the equation simply becomes:
    !
    ! (dN_r/dt)|_evap = ( N_r / r_r ) (dr_r/dt)|_evap.

    ! References:
    ! Eq. (23) of Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
    ! Parameterization in a Large-Eddy Simulation Model of Marine Stratocumulus.
    ! Mon. Wea. Rev., 128, 229--243.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one  ! Constant(s)

    use parameters_KK, only: &
        KK_Nrm_evap_nu  ! Variable(s)

    use clubb_precision, only: &
        core_rknd

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      KK_rrm_evap_tndcy, & ! Mean < dr_r/dt > due to evaporation    [(kg/kg)/s]
      Nrm,               & ! Mean rain drop concentration, < N_r >  [num/kg]
      rrm                  ! Mean rain water mixing ratio, < r_r >  [kg/kg]

    real( kind = core_rknd ), intent(in) :: &
      dt                   ! Model time step duration               [s]

    ! Return Variables
    real( kind = core_rknd ) :: &
      KK_Nrm_evap_local_mean  ! Mean of KK <N_r> evaporation rate   [(num/kg)/s]
 

    ! Calculate the local KK < N_r > evaporation tendency.
    KK_Nrm_evap_local_mean &
    = dt**( KK_Nrm_evap_nu - one ) &
      * ( Nrm / rrm**KK_Nrm_evap_nu ) * KK_rrm_evap_tndcy**KK_Nrm_evap_nu


    return

  end function KK_Nrm_evap_local_mean

  !=============================================================================
  function KK_Nrm_auto_mean( KK_rrm_auto_tndcy )

    ! Description:
    ! This function calculates the production rate of rain drop concentration
    ! from the process of autoconversion of cloud droplets into rain drops.
    ! The KK equation for the source of rain drop concentration due to
    ! autoconversion is:
    !
    ! ( dN_rV / dt )_auto
    !    = ( dr_r / dt )_auto / ( ( 4*pi*rho_lw / (3*rho_a) ) * r_0^3 ).
    !
    ! Since Nr (in num/kg) = N_rV (in num/m^3) / rho_a, the equation becomes:
    !
    ! ( dN_r / dt )_auto
    !    = ( dr_r / dt )_auto / ( ( 4*pi*rho_lw / 3 ) * r_0^3 ).
    !
    ! CLUBB uses N_r specified in terms of num/kg.

    ! References:
    ! Eq. (32) of Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
    ! Parameterization in a Large-Eddy Simulation Model of Marine Stratocumulus.
    ! Mon. Wea. Rev., 128, 229--243.
    !-----------------------------------------------------------------------

    use constants_clubb, only: & 
        rho_lw,      & ! Constant(s)
        pi,          &
        four_thirds

    use parameters_KK, only: &
        r_0 ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      KK_rrm_auto_tndcy  ! Mean KK <r_r> autoconversion rate   [(kg/kg)/s]

    ! Return Variable
    real( kind = core_rknd ) :: &
      KK_Nrm_auto_mean  ! Mean KK <N_r> autoconversion rate    [(num/kg)/s]


    ! Production of N_r through autoconversion.
    KK_Nrm_auto_mean &
    = KK_rrm_auto_tndcy / ( four_thirds * pi * rho_lw * r_0**3 )


    return

  end function KK_Nrm_auto_mean

!===============================================================================

end module KK_Nrm_tendencies
