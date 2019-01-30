! $Id$
!===============================================================================
module KK_local_means

  implicit none

  private

  public :: KK_evap_local_mean, &
            KK_auto_local_mean, &
            KK_accr_local_mean, &
            KK_mvr_local_mean

  contains

  !=============================================================================
  function KK_evap_local_mean( mean_s, rrm, Nrm, KK_evap_coef )

    ! Description:

    ! References:
    ! Eq. (22) of Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
    ! Parameterization in a Large-Eddy Simulation Model of Marine Stratocumulus.
    ! Mon. Wea. Rev., 128, 229--243.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero  ! Constant(s)

    use parameters_KK, only: &
        KK_evap_Supersat_exp, & ! Variable(s)
        KK_evap_rr_exp,       &
        KK_evap_Nr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mean_s,       & ! Mean value of extended liquid water mixing ratio     [-]
      rrm,          & ! Mean rain water mixing ratio                         [-]
      Nrm,          & ! Mean rain drop concentration                         [-]
      KK_evap_coef    ! KK evaporation coefficient                   [(kg/kg)/s]

    ! Return Variable
    real( kind = core_rknd ) :: &
      KK_evap_local_mean  ! Mean of KK evaporation tendency          [(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on s                                           [-]
      beta_exp,  & ! Exponent on r_r                                         [-]
      gamma_exp    ! Exponent on N_r                                         [-]


    ! Values of the KK exponents.
    alpha_exp = KK_evap_Supersat_exp
    beta_exp  = KK_evap_rr_exp
    gamma_exp = KK_evap_Nr_exp

    ! Calculate the local KK evaporation tendency.
    if ( mean_s <= zero ) then

       KK_evap_local_mean  &
       = KK_evap_coef * mean_s**alpha_exp * rrm**beta_exp * Nrm**gamma_exp

    else  ! mean_s > 0.

       ! The air is cloudy, and rain does not evaporate.
       KK_evap_local_mean = zero

    endif


    return

  end function KK_evap_local_mean

  !=============================================================================
  function KK_auto_local_mean( mean_s, Ncm, KK_auto_coef )

    ! Description:

    ! References:
    ! Eq. (29) of Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
    ! Parameterization in a Large-Eddy Simulation Model of Marine Stratocumulus.
    ! Mon. Wea. Rev., 128, 229--243.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero  ! Constant(s)

    use parameters_KK, only: &
        KK_auto_rc_exp, & ! Variable(s)
        KK_auto_Nc_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mean_s,       & ! Mean value of extended liquid water mixing ratio    [-]
      Ncm,          & ! Mean cloud droplet concentration                    [-]
      KK_auto_coef    ! KK autoconversion coefficient               [(kg/kg)/s]

    ! Return Variable
    real( kind = core_rknd ) :: &
      KK_auto_local_mean  ! Mean of KK autoconversion tendency      [(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on s                                          [-]
      beta_exp     ! Exponent on N_c                                        [-]


    ! Values of the KK exponents.
    alpha_exp = KK_auto_rc_exp
    beta_exp  = KK_auto_Nc_exp

    ! Calculate the local KK autoconversion tendency.
    if ( mean_s > zero ) then

       KK_auto_local_mean = KK_auto_coef * mean_s**alpha_exp * Ncm**beta_exp

    else  ! mean_s <= 0.

       ! The air is clear, and rain cannot be produced.
       KK_auto_local_mean = zero

    endif


    return

  end function KK_auto_local_mean

  !=============================================================================
  function KK_accr_local_mean( mean_s, rrm, KK_accr_coef )

    ! Description:

    ! References:
    ! Eq. (33) of Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
    ! Parameterization in a Large-Eddy Simulation Model of Marine Stratocumulus.
    ! Mon. Wea. Rev., 128, 229--243.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero  ! Constant(s)

    use parameters_KK, only: &
        KK_accr_rc_exp, & ! Variable(s)
        KK_accr_rr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mean_s,       & ! Mean value of extended liquid water mixing ratio    [-]
      rrm,          & ! Mean rain water mixing ratio                        [-]
      KK_accr_coef    ! KK accretion coefficient                    [(kg/kg)/s]

    ! Return Variable
    real( kind = core_rknd ) :: &
      KK_accr_local_mean  ! Mean of KK accretion tendency           [(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on s                                          [-]
      beta_exp     ! Exponent on r_r                                        [-]


    ! Values of the KK exponents.
    alpha_exp = KK_accr_rc_exp
    beta_exp  = KK_accr_rr_exp

    ! Calculate the local KK accretion tendency.
    if ( mean_s > zero ) then

       KK_accr_local_mean = KK_accr_coef * mean_s**alpha_exp * rrm**beta_exp

    else  ! mean_s <= 0.

       ! The air is clear, and rain cannot be produced.
       KK_accr_local_mean = zero

    endif


    return

  end function KK_accr_local_mean

  !=============================================================================
  function KK_mvr_local_mean( rrm, Nrm, KK_mvr_coef )

    ! Description:

    ! References:
    ! Eq. (3) of Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
    ! Parameterization in a Large-Eddy Simulation Model of Marine Stratocumulus.
    ! Mon. Wea. Rev., 128, 229--243.
    !-----------------------------------------------------------------------

    use parameters_KK, only: &
        KK_mvr_rr_exp, & ! Variable(s)
        KK_mvr_Nr_exp

    use constants_clubb, only: &
        Nr_tol    ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rrm,         & ! Mean rain water mixing ratio        [-]
      Nrm,         & ! Mean rain drop concentration        [-]
      KK_mvr_coef    ! KK mean volume radius coefficient   [m]

    ! Return Variable
    real( kind = core_rknd ) :: &
      KK_mvr_local_mean  ! Rain drop mean volume radius    [m]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on r_r                       [-]
      beta_exp     ! Exponent on N_r                       [-]


    ! Values of the KK exponents.
    alpha_exp = KK_mvr_rr_exp
    beta_exp  = KK_mvr_Nr_exp

    ! Calculate the KK mean volume radius of rain drops
    KK_mvr_local_mean &
    = KK_mvr_coef * rrm**alpha_exp * max( Nrm, Nr_tol )**beta_exp


    return

  end function KK_mvr_local_mean

!===============================================================================

end module KK_local_means
