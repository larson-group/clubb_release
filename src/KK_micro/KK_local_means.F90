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
    !-----------------------------------------------------------------------

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

    ! Constant Parameters
    real( kind = core_rknd ), parameter :: &
      alpha_exp = 1.0,       & ! Exponent on s                               [-]
      beta_exp  = (1.0/3.0), & ! Exponent on r_r                             [-]
      gamma_exp = (2.0/3.0)    ! Exponent on N_r                             [-]


    ! Calculate the local KK evaporation tendency.
    if ( mean_s <= 0.0 ) then

       KK_evap_local_mean  &
       = KK_evap_coef * mean_s**alpha_exp * rrm**beta_exp * Nrm**gamma_exp

    else  ! mean_s > 0.

       ! The air is cloudy, and rain does not evaporate.
       KK_evap_local_mean = 0.0

    endif


    return

  end function KK_evap_local_mean

  !=============================================================================
  function KK_auto_local_mean( mean_s, Ncm, KK_auto_coef )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

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

    ! Constant Parameters
    real( kind = core_rknd ), parameter :: &
      alpha_exp = 2.47,  & ! Exponent on s                                  [-]
      beta_exp  = -1.79    ! Exponent on N_c                                [-]


    ! Calculate the local KK autoconversion tendency.
    if ( mean_s > 0.0 ) then

       KK_auto_local_mean = KK_auto_coef * mean_s**alpha_exp * Ncm**beta_exp

    else  ! mean_s <= 0.

       ! The air is clear, and rain cannot be produced.
       KK_auto_local_mean = 0.0

    endif


    return

  end function KK_auto_local_mean

  !=============================================================================
  function KK_accr_local_mean( mean_s, rrm, KK_accr_coef )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

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

    ! Constant Parameters
    real( kind = core_rknd ), parameter :: &
      alpha_exp = 1.15, & ! Exponent on s                                   [-]
      beta_exp  = 1.15    ! Exponent on r_r                                 [-]


    ! Calculate the local KK accretion tendency.
    if ( mean_s > 0.0 ) then

       KK_accr_local_mean = KK_accr_coef * mean_s**alpha_exp * rrm**beta_exp

    else  ! mean_s <= 0.

       ! The air is clear, and rain cannot be produced.
       KK_accr_local_mean = 0.0

    endif


    return

  end function KK_accr_local_mean

  !=============================================================================
  function KK_mvr_local_mean( rrm, Nrm, KK_mvr_coef )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

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
      KK_mvr_local_mean  ! Mean

    ! Constant Parameters
    real( kind = core_rknd ), parameter :: &
      alpha_exp = (1.0/3.0),  & ! Exponent on r_r          [-]
      beta_exp  = (-1.0/3.0)    ! Exponent on N_r          [-]


    ! Calculate the KK mean volume radius of rain drops
    KK_mvr_local_mean = KK_mvr_coef * rrm**alpha_exp * Nrm**beta_exp


    return

  end function KK_mvr_local_mean

!===============================================================================

end module KK_local_means
