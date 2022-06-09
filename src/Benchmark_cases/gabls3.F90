!----------------------------------------------------------------------
! $Id$
module gabls3

  !       Description:
  !       Contains subroutines for the GABLS3.
  !     
  !       References:
  !       http://www.knmi.nl/samenw/gabls/   
  !----------------------------------------------------------------------

  implicit none

  public :: gabls3_sfclyr

  private

  contains

  !-----------------------------------------------------------------------
  subroutine gabls3_sfclyr( ubar, veg_t_in_K, &
                            thlm_sfc, rtm_sfc, lowest_level, exner_sfc, &
                            wpthlp_sfc, wprtp_sfc, ustar )
    !       Description:
    !       This subroutine computes surface fluxes of horizontal momentum,
    !       heat and moisture according to GCSS ATEX specifications

    !       References:
    !       http://www.knmi.nl/samenw/gabls/
    !----------------------------------------------------------------------

    use constants_clubb, only: grav ! Variable(s)

    use diag_ustar_module, only: diag_ustar ! Procedure(s)

    use sfc_flux, only: compute_wpthlp_sfc, compute_wprtp_sfc ! Procedure(s)

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

    ! Constants

    real( kind = core_rknd ), parameter ::  & 
     ! ustar = 0.3_core_rknd,
     ! C_10  = 0.0013_core_rknd, & !ATEX value
     ! C_10  = 0.013_core_rknd, & ! Fudged value
     ! C_10  = 0.0049_core_rknd, & ! Fudged value
     ! C_10  = 0.0039_core_rknd, & ! Fudged value
      C_10 = 0.00195_core_rknd, &
     ! C_10 = 0.001_core_rknd, &
     ! C_10 = 0.003_core_rknd, &
      z0 = 0.15_core_rknd

!    real, parameter, dimension(25) :: T_sfc_given = &
!         (/300._core_rknd, 300.8_core_rknd, 300.9_core_rknd, 301._core_rknd, &
!         300.9_core_rknd, 300.5_core_rknd, 300._core_rknd, 298.5_core_rknd, &
!         297._core_rknd, 296._core_rknd, 295._core_rknd, 294._core_rknd, &
!         293.5_core_rknd, 292.5_core_rknd, 291.5_core_rknd, 291._core_rknd,&
!         290.5_core_rknd, 292.5_core_rknd, 294.5_core_rknd, 296.5_core_rknd, &
!         298._core_rknd, 298.5_core_rknd, 300.5_core_rknd, 301.5_core_rknd, &
!         301._core_rknd/)
!    real, parameter, dimension(25) :: T_sfc_time = (/43200._core_rknd, &
!         46800._core_rknd, 50400._core_rknd, 54000._core_rknd, 57600._core_rknd, &
!         61200._core_rknd, 64800._core_rknd, 68400._core_rknd, 72000._core_rknd, &
!         75600._core_rknd, 79200._core_rknd, 82800._core_rknd, 86400._core_rknd, &
!         90000._core_rknd, 93600._core_rknd, 97200._core_rknd, 100800._core_rknd, &
!         104400._core_rknd, 108000._core_rknd, 111600._core_rknd, 115200._core_rknd, &
!         118800._core_rknd, 122400._core_rknd, 126000._core_rknd, 129600._core_rknd/)


    real( kind = core_rknd ), intent(in) ::  &
      ubar,         & ! mean sfc wind speed    [m/s]
      thlm_sfc,     & ! Theta_l at zt(2)       [K]
      rtm_sfc,      & ! rt at zt(2)            [kg/kg]
      veg_T_in_K,   & ! Vegetation temperature [K]
      lowest_level, & ! gr%zt(1,2)               [m]
      exner_sfc       ! Exner function         [-]

    ! Output variables
    real( kind = core_rknd ), intent(inout):: &
      wpthlp_sfc,  & ! w'theta_l' surface flux   [(m K)/s]
      wprtp_sfc      ! w'rt' surface flux        [(m kg)/(kg s)]

    real( kind = core_rknd ), intent(out) ::  & 
      ustar          ! surface friction velocity [m/s]

    ! Local Variables
    real( kind = core_rknd ) :: veg_theta_in_K, bflx, offset

    ! Compute heat and moisture fluxes

    veg_theta_in_K = veg_T_in_K / exner_sfc

    wpthlp_sfc = compute_wpthlp_sfc( C_10, ubar, thlm_sfc, veg_T_in_K, exner_sfc)
    offset = 9.9e-3_core_rknd
    wprtp_sfc = compute_wprtp_sfc( C_10, ubar, rtm_sfc, offset ) * &
       10._core_rknd ! Known magic number
  
    ! Compute momentum fluxes
    bflx = wpthlp_sfc * grav / veg_theta_in_K

    ustar = diag_ustar( lowest_level, bflx, ubar, z0)

    return
  end subroutine gabls3_sfclyr

end module gabls3
