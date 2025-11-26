!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module parameters_KK

! Description:
!   Parameters for KK microphysics

! References:
!   None
!-------------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd

  use constants_clubb, only: &
    one, &
    one_third, &
    two_thirds

  implicit none

  ! Values of exponents in KK microphysics
  real( kind = core_rknd ), parameter, public :: &
    KK_evap_Supersat_exp = one,  & ! Exponent on Supersaturation (S) in KK evap. eq.; 1
    KK_evap_rr_exp = one_third,  & ! Exponent on r_r in KK evaporation eq.; 1/3
    KK_evap_Nr_exp = two_thirds, & ! Exponent on N_r in KK evaporation eq.; 2/3
    KK_auto_rc_exp = 2.47_core_rknd, & ! Exponent on r_c in KK autoconversion eq.; 2.47
    KK_auto_Nc_exp = -1.79_core_rknd, & ! Exponent on N_c in KK autoconversion eq.; -1.79
    KK_accr_rc_exp = 1.15_core_rknd, & ! Exponent on r_c in KK accretion eq.; 1.15
    KK_accr_rr_exp = 1.15_core_rknd, & ! Exponent on r_r in KK accretion eq.; 1.15
    KK_mvr_rr_exp  = one_third, & ! Exponent on r_r in KK mean volume radius
    KK_mvr_Nr_exp  = -one_third, & ! Exponent on N_r in KK mean volume radius
    KK_Nrm_evap_nu = one    ! Exponent (parameter) in <N_r> evaporation eq.; 1

  !--------------

  ! The following variables are not parameters because they can be changed through the
  ! microphysics_setting namelist!

  ! Assumed radius of all new drops; m.
  real( kind = core_rknd ), public :: r_0 = 25.0e-6_core_rknd

  ! Other needed parameters
  ! Khairoutdinov and Kogan (2000) ratio of drizzle drop mean geometric
  ! radius to drizzle drop mean volume radius.
  ! Khairoutdinov and Kogan (2000); p. 233
  real( kind = core_rknd ), public :: C_evap = 0.86_core_rknd

  private ! Default Scope

end module parameters_KK
