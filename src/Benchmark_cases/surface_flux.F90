!$Id$
module surface_flux
!
!       This module contains generalized subroutines for determining surface
!       fluxes.
!
  implicit none

  public :: compute_momentum_flux, compute_ubar, compute_wprtp_sfc, &
            compute_wpthlp_sfc

  private

  contains

  real function compute_ubar( um_sfc, vm_sfc )
    !
    !  Description: This function determins the value of ubar
    !               based on the momentum at the surface.
    !
    !-----------------------------------------------------------
    
    implicit none

    ! Input Variable(s)
    real, intent(in) :: um_sfc, & ! u wind at zt(2) [m/s]
                        vm_sfc    ! v wind at zt(2) [m/s]


    real, parameter :: ubmin = 0.25

    compute_ubar = max( ubmin, sqrt( um_sfc**2 + vm_sfc**2 ) )

  end function compute_ubar

  !-----------------------------------------------------------
  subroutine compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                                    upwp_sfc, vpwp_sfc )
    !
    !  Description: This subroutine computes the momentum fluxes
    !               upwp_sfc and vpwp_sfc.
    !
    !-----------------------------------------------------------

    implicit none

    ! Input
    real, intent(in) :: um_sfc, & ! u wind at zt(2) [m/s]
                        vm_sfc, & ! v wind at zt(2) [m/s]
                        ustar, &  ! Surface friction velocity [m/s]
                        ubar      

    ! Output
    real, intent(out) :: upwp_sfc, &  ! Turbulent upward flux of u-momentum [(m/s)^2]
                         vpwp_sfc      ! Turbulent upward flux of v-momentum [(m/s)^2]




    ! Compute momentum fluxes

    upwp_sfc = -um_sfc * ustar**2 / ubar
    vpwp_sfc = -vm_sfc * ustar**2 / ubar

  end subroutine compute_momentum_flux

  real function compute_wpthlp_sfc( Cd, ubar, thlm_sfc, Tsfc, exner_sfc )
  !
  !  Description: This function determins the surface flux of heat.
  !
  !----------------------------------------------------------------------
    implicit none
    ! Intent(in)
    real, intent(in) :: Cd, &     ! Coefficient
                        ubar, &
                        thlm_sfc, & ! Theta_l at zt(2) [K]
                        Tsfc, &     ! Surface temperature [K]
                        exner_sfc   ! Exner function at surface [-]

    compute_wpthlp_sfc = -Cd * ubar * ( thlm_sfc - Tsfc / exner_sfc )


  end function compute_wpthlp_sfc


  !-----------------------------------------------------------------------
  real function compute_wprtp_sfc ( Cd, ubar, rtm_sfc, adjustment )
    !
    !  Description: This function determines the surface flux of moisture.
    !
    !--------------------------------------------------------------------
  use saturation, only:  sat_mixrat_liq

  implicit none
        ! Input(s)
        real, intent(in) :: Cd, &
                            ubar, &
                            rtm_sfc, &  ! Total Water mixing ratio at zt(2) [kg/kg]
                            adjustment

        compute_wprtp_sfc  = -Cd * ubar * ( rtm_sfc - adjustment )


  end function compute_wprtp_sfc
end module surface_flux
