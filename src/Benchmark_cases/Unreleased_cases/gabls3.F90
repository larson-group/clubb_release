!----------------------------------------------------------------------
! $Id$
module gabls3

  !       Description:
  !       Contains subroutines for the GABLS3.
  !----------------------------------------------------------------------

  implicit none

  public :: gabls3_sfclyr

  private

  contains

  !-----------------------------------------------------------------------
  subroutine gabls3_sfclyr( um_sfc, vm_sfc, veg_t_in_K, &
                            thlm_sfc, rtm_sfc, lowest_level, exner_sfc, & 
                            upwp_sfc, vpwp_sfc, &
                            wpthlp_sfc, wprtp_sfc, ustar )
    !       Description:
    !       This subroutine computes surface fluxes of horizontal momentum,
    !       heat and moisture according to GCSS ATEX specifications

    !       References:

    !----------------------------------------------------------------------

    use constants_clubb, only: kappa, grav, Rd, Cp, p0, Lv ! Variable(s)

    use diag_ustar_module, only: diag_ustar ! Procedure(s)

    use stats_precision, only: time_precision ! Variable(s)

    use surface_flux, only: compute_momentum_flux, &
                            compute_wpthlp_sfc, &
                            compute_wprtp_sfc

    !use surface, only: prognose_soil_T_in_K ! Procedure(s)

    implicit none

    ! Constants

    real, parameter ::  & 
      ubmin = 0.25, & 
     ! ustar = 0.3,
     ! C_10  = 0.0013, & !ATEX value
     ! C_10  = 0.013, & ! Fudged value
     ! C_10  = 0.0049, & ! Fudged value
     ! C_10  = 0.0039, & ! Fudged value
      C_10 = 0.00195, &
     ! C_10 = 0.001, &
     ! C_10 = 0.003, &
      z0 = 0.15

!    real, parameter, dimension(25) :: sst_given = (/300., 300.8, 300.9, 301.,300.9, &
!                                        300.5, 300., 298.5, 297., 296., 295.,&
!                                        294., 293.5, 292.5, 291.5, 291.,&
!                                        290.5, 292.5, 294.5, 296.5, 298.,&
!                                        298.5, 300.5, 301.5, 301./)
!    real, parameter, dimension(25) :: sst_time = (/43200., 46800., 50400., 54000., 57600., &
!                                        61200., 64800., 68400., 72000., 75600.,&
!                                        79200., 82800., 86400., 90000., 93600., &
!                                        97200., 100800., 104400., 108000., 111600.,&
!                                        115200., 118800., 122400., 126000., 129600./)


    ! Input variables


    real, intent(in) ::  & 
      um_sfc,       & ! um at zt(2)            [m/s]
      vm_sfc,       & ! vm at zt(2)            [m/s]
      thlm_sfc,     & ! Theta_l at zt(2)       [K]
      rtm_sfc,      & ! rt at zt(2)            [kg/kg]
      veg_T_in_K,   & ! Vegetation temperature [K]
      lowest_level, & ! gr%zt(2)               [m]
      exner_sfc       ! Exner function         [-]

    ! Output variables
    real, intent(out) ::  & 
      upwp_sfc,    & ! u'w' at surface           [m^2/s^2]
      vpwp_sfc,    & ! v'w' at surface           [m^2/s^2]
      ustar          ! surface friction velocity [m/s]

    real, intent(inout):: &
      wpthlp_sfc,  & ! w'theta_l' surface flux   [(m K)/s]
      wprtp_sfc      ! w'rt' surface flux        [(m kg)/(kg s)]

    ! Local Variables
    real :: ubar, veg_theta_in_K, bflx

    ! Compute heat and moisture fluxes
    ubar = max( ubmin, sqrt( um_sfc**2 + vm_sfc**2 ) )

    veg_theta_in_K = veg_T_in_K / exner_sfc

    wpthlp_sfc = compute_wpthlp_sfc( C_10, ubar, thlm_sfc, veg_T_in_K, exner_sfc)

    wprtp_sfc = compute_wprtp_sfc( C_10, ubar, rtm_sfc, 9.9e-3 ) * 10
  
    ! Compute momentum fluxes
    bflx = wpthlp_sfc * grav / veg_theta_in_K

    ustar = diag_ustar( lowest_level, bflx, ubar, z0)

    call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                                upwp_sfc, vpwp_sfc )

    return
  end subroutine gabls3_sfclyr

end module gabls3
