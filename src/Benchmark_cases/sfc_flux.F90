!$Id$
!-------------------------------------------------------------------------------
module sfc_flux
!  Description:
!    This module contains generalized subroutines for determining surface
!    fluxes.
!  References:
!    None
!-------------------------------------------------------------------------------
  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  public :: compute_momentum_flux, compute_ubar, compute_ht_mostr_flux, &
            compute_wprtp_sfc, compute_wpthlp_sfc, set_sclr_sfc_rtm_thlm, &
            convert_sens_ht_to_km_s, convert_latent_ht_to_m_s

  private

  contains

!===============================================================================
  subroutine compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                                    upwp_sfc, vpwp_sfc )
!
! Description:
!   This subroutine computes the momentum fluxes upwp_sfc and vpwp_sfc.
! References:
!   None
!-------------------------------------------------------------------------------

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

    ! Input
    real( kind = core_rknd ), intent(in) :: &
      um_sfc, & ! u wind at zt(2) [m/s]
      vm_sfc, & ! v wind at zt(2) [m/s]
      ustar,  & ! Surface friction velocity [m/s]
      ubar

    ! Output
    real( kind = core_rknd ), intent(out) :: &
      upwp_sfc, &  ! Turbulent upward flux of u-momentum [(m/s)^2]
      vpwp_sfc     ! Turbulent upward flux of v-momentum [(m/s)^2]

    ! ---- Begin Code ----

    ! Compute momentum fluxes

    upwp_sfc = -um_sfc * ustar**2 / ubar
    vpwp_sfc = -vm_sfc * ustar**2 / ubar

    return
  end subroutine compute_momentum_flux
  
!===============================================================================
  real( kind = core_rknd ) function compute_ubar( um_sfc, vm_sfc )
!
!  Description:
!    This function determines the value of ubar based on the momentum at 
!    the surface.
!  References:
!    None
!-------------------------------------------------------------------------------

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: max, sqrt

    ! Constant paramter(s)
    real( kind = core_rknd ), parameter :: ubmin = 0.25_core_rknd ! [m/s]

    ! Input Variable(s)
    real( kind = core_rknd ), intent(in) :: &
      um_sfc, & ! u wind at zt(2) [m/s]
      vm_sfc    ! v wind at zt(2) [m/s]

    ! ---- Begin Code ----

    compute_ubar = max( ubmin, sqrt( um_sfc**2 + vm_sfc**2 ) )

  end function compute_ubar
  
!===============================================================================
  subroutine compute_ht_mostr_flux( time_in, ntimes, &
                                    heat_flx, moisture_flx )

!
! Description:
!   Compute heat and moisture fluxes from ARM data in (W/m2)
! References:
!   None
!-------------------------------------------------------------------------------
    use clubb_precision, only: time_precision, core_rknd ! Variable(s)
    
    use time_dependent_input, only: latent_ht_given, sens_ht_given, time_sfc_given ! Variable(s)
    
    use interpolation, only: linear_interp_factor ! Procedure(s)
    
    implicit none

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      time_in    ! Current time simulation time [s]

    integer, intent(in) :: &
      ntimes

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      heat_flx, &
      moisture_flx

    ! Local variables
    integer :: &
      i1, &
      i2

    real( kind = core_rknd ) ::  & 
      time_frac, time

!---------------------BEGIN CODE-------------------------

    ! Default Initialization
    heat_flx = 0.0_core_rknd
    moisture_flx = 0.0_core_rknd
   
    time = real( time_in, kind = core_rknd ) 

    ! Compute heat and moisture fluxes from ARM data in (W/m2)
    if ( time <= time_sfc_given(1) ) then
      heat_flx     = sens_ht_given(1)
      moisture_flx = latent_ht_given(1)
      
    else if ( time >= time_sfc_given(ntimes) ) then
      heat_flx     = sens_ht_given(ntimes)
      moisture_flx = latent_ht_given(ntimes)
      
    else
      i1 = 1
      do while ( i1 <= ntimes-1 )
        i2 = i1 + 1
        
        if ( time >= time_sfc_given(i1) .and. time < time_sfc_given(i2) ) then
          time_frac = real((time-time_sfc_given(i1))/(time_sfc_given(i2) - &
            time_sfc_given(i1)), kind = core_rknd)
          heat_flx = linear_interp_factor( time_frac, sens_ht_given(i2), sens_ht_given(i1) )
          moisture_flx = linear_interp_factor( time_frac, latent_ht_given(i2), latent_ht_given(i1) )
          i1 = ntimes
        end if
        
        i1 = i2
      end do
    end if ! time <= time_sfc_given(1)

    return
  end subroutine compute_ht_mostr_flux

!===============================================================================
  real( kind = core_rknd ) function compute_wpthlp_sfc( Cd, ubar, thlm_sfc, T_sfc, exner_sfc )
!
! Description:
!   This function determins the surface flux of heat.
! References:
!   None
!-------------------------------------------------------------------------------

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

    ! Intent(in)
    real( kind = core_rknd ), intent(in) :: &
      Cd,       & ! Coefficient
      ubar,     &
      thlm_sfc, & ! Theta_l at zt(2) [K]
      T_sfc,     & ! Surface temperature [K]
      exner_sfc   ! Exner function at surface [-]

    ! ---- Begin Code ----

    compute_wpthlp_sfc = -Cd * ubar * ( thlm_sfc - T_sfc / exner_sfc )

    return
  end function compute_wpthlp_sfc


!===============================================================================
  real( kind = core_rknd ) function compute_wprtp_sfc ( Cd, ubar, rtm_sfc, adjustment )

!
!  Description:
!    This function determines the surface flux of moisture.
!  References:
!    None
!-------------------------------------------------------------------------------

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

    ! Input(s)
    real( kind = core_rknd ), intent(in) :: &
      Cd, &
      ubar, &
      rtm_sfc, &  ! Total Water mixing ratio at zt(2) [kg/kg]
      adjustment

    ! ---- Begin Code ----

    compute_wprtp_sfc  = -Cd * ubar * ( rtm_sfc - adjustment )

    return
  end function compute_wprtp_sfc
  
!===============================================================================
  subroutine set_sclr_sfc_rtm_thlm ( sclr_dim, edsclr_dim, sclr_idx, &
                                     wpthlp_sfc, wprtp_sfc, &
                                     wpsclrp_sfc, wpedsclrp_sfc )

!
!  Description:
!    This function determines the surface flux of moisture.
!  References:
!    None
!-------------------------------------------------------------------------------

    use array_index, only: &
      sclr_idx_type

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      sclr_dim, & 
      edsclr_dim

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    real( kind = core_rknd ), intent(in) ::  & 
      wpthlp_sfc, & ! surface thetal flux        [K m/s]
      wprtp_sfc     ! surface moisture flux      [kg/kg m/s]
      
    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), intent(out), dimension(sclr_dim) ::  & 
      wpsclrp_sfc       ! scalar surface flux            [units m/s]

    real( kind = core_rknd ), intent(out), dimension(edsclr_dim) ::  & 
      wpedsclrp_sfc     ! eddy-scalar surface flux       [units m/s]

    !--------------------- Locally Variables ---------------------


    !--------------------- Begin Code ---------------------
    
    wpsclrp_sfc(:)   = 0.0_core_rknd ! Initialize flux to 0 
    wpedsclrp_sfc(:) = 0.0_core_rknd ! Initialize flux to 0 

    ! Let passive scalars be equal to rt and theta_l for now
    if ( sclr_idx%iisclr_thl > 0 ) wpsclrp_sfc(sclr_idx%iisclr_thl) = wpthlp_sfc
    if ( sclr_idx%iisclr_rt  > 0 ) wpsclrp_sfc(sclr_idx%iisclr_rt)  = wprtp_sfc

    if ( sclr_idx%iiedsclr_thl > 0 ) wpedsclrp_sfc(sclr_idx%iiedsclr_thl) = wpthlp_sfc
    if ( sclr_idx%iiedsclr_rt  > 0 ) wpedsclrp_sfc(sclr_idx%iiedsclr_rt)  = wprtp_sfc

    return
  end subroutine set_sclr_sfc_rtm_thlm

 
!==============================================================================
  real( kind = core_rknd ) function convert_sens_ht_to_km_s ( sens_ht, rho_sfc )

!   This function converts sensible heat flux in W/m^2 to
!   natural units of k m/s for the wpthlp_sfc variable.
!-----------------------------------------------------------------------------  

    use constants_clubb, only: Cp ! Variable(s)

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

    real( kind = core_rknd ), intent(in) :: &
      sens_ht,               & ! Sensible heat flux     [W/m^2]
      rho_sfc                ! Density at the surface [kg/m^3]

    !--------------------BEGIN CODE-----------------------

    convert_sens_ht_to_km_s = sens_ht / ( rho_sfc * Cp)

    return
  
  end function convert_sens_ht_to_km_s


!==============================================================================
  real( kind = core_rknd ) function convert_latent_ht_to_m_s ( latent_ht, rho_sfc )

!   This function converts latent heat flux in W/m^2 to
!   natural units of m/s for the wprtp_sfc variable.
!-----------------------------------------------------------------------------  

    use constants_clubb, only: Lv ! Variable(s)

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none
    
    real( kind = core_rknd ), intent(in) :: &
      latent_ht,               & ! Latent heat flux     [W/m^2]
      rho_sfc                ! Density at the surface [kg/m^3]

    !--------------------BEGIN CODE-----------------------

    convert_latent_ht_to_m_s = latent_ht / ( rho_sfc * Lv)

    return
  
  end function convert_latent_ht_to_m_s

end module sfc_flux
