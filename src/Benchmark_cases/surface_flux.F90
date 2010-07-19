!$Id$
!-------------------------------------------------------------------------------
module surface_flux
!  Description:
!    This module contains generalized subroutines for determining surface
!    fluxes.
!  References:
!    None
!-------------------------------------------------------------------------------
  implicit none

  public :: compute_momentum_flux, compute_ubar, compute_ht_mostr_flux, &
            compute_wprtp_sfc, compute_wpthlp_sfc, set_sclr_sfc_rtm_thlm

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

    implicit none

    ! Input
    real, intent(in) :: &
      um_sfc, & ! u wind at zt(2) [m/s]
      vm_sfc, & ! v wind at zt(2) [m/s]
      ustar,  & ! Surface friction velocity [m/s]
      ubar

    ! Output
    real, intent(out) :: &
      upwp_sfc, &  ! Turbulent upward flux of u-momentum [(m/s)^2]
      vpwp_sfc     ! Turbulent upward flux of v-momentum [(m/s)^2]

    ! ---- Begin Code ----

    ! Compute momentum fluxes

    upwp_sfc = -um_sfc * ustar**2 / ubar
    vpwp_sfc = -vm_sfc * ustar**2 / ubar

    return
  end subroutine compute_momentum_flux
  
!===============================================================================
  real function compute_ubar( um_sfc, vm_sfc )
!
!  Description:
!    This function determines the value of ubar based on the momentum at 
!    the surface.
!  References:
!    None
!-------------------------------------------------------------------------------

    implicit none

    ! External
    intrinsic :: max, sqrt

    ! Constant paramter(s)
    real, parameter :: ubmin = 0.25 ! [m/s]

    ! Input Variable(s)
    real, intent(in) :: &
      um_sfc, & ! u wind at zt(2) [m/s]
      vm_sfc    ! v wind at zt(2) [m/s]

    ! ---- Begin Code ----

    compute_ubar = max( ubmin, sqrt( um_sfc**2 + vm_sfc**2 ) )

  end function compute_ubar
  
!===============================================================================
  subroutine compute_ht_mostr_flux( time, ntimes, &
                                    heat_flx, moisture_flx )

!
! Description:
!   Compute heat and moisture fluxes from ARM data in (W/m2)
! References:
!   None
!-------------------------------------------------------------------------------
    use stats_precision, only: time_precision ! Variable(s)
    
    use time_dependent_input, only: LH_given, SH_given, time_sfc_given ! Variable(s)
    
    use interpolation, only: factor_interp ! Procedure(s)
    
    implicit none

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      time     ! Current time        [s]
    integer, intent(in) :: &
      ntimes

    ! Output Variables
    real, intent(out) :: &
      heat_flx, &
      moisture_flx

    ! Local variables
    integer :: &
      i1, &
      i2
    real ::  & 
      time_frac

!---------------------BEGIN CODE-------------------------

    ! Default Initialization
    heat_flx = 0.0
    moisture_flx = 0.0
    
    ! Compute heat and moisture fluxes from ARM data in (W/m2)
    if ( time <= time_sfc_given(1) ) then
      heat_flx     = SH_given(1)
      moisture_flx = LH_given(1)
      
    else if ( time >= time_sfc_given(ntimes) ) then
      heat_flx     = SH_given(ntimes)
      moisture_flx = LH_given(ntimes)
      
    else
      i1 = 1
      do while ( i1 <= ntimes-1 )
        i2 = i1 + 1
        
        if ( time >= time_sfc_given(i1) .and. time < time_sfc_given(i2) ) then
          time_frac = real((time-time_sfc_given(i1))/(time_sfc_given(i2) - time_sfc_given(i1)))
!          heat_flx     = ( 1. - time_frac ) * SH_given(i1) + time_frac * SH_given(i2)
          heat_flx = factor_interp( time_frac, SH_given(i2), SH_given(i1) )
!          moisture_flx = ( 1. - time_frac ) * LH_given(i1) + time_frac * LH_given(i2)
          moisture_flx = factor_interp( time_frac, LH_given(i2), LH_given(i1) )
          i1 = ntimes
        end if
        
        i1 = i2
      end do
    end if ! time <= time_sfc_given(1)

    return
  end subroutine compute_ht_mostr_flux

!===============================================================================
  real function compute_wpthlp_sfc( Cd, ubar, thlm_sfc, Tsfc, exner_sfc )
!
! Description:
!   This function determins the surface flux of heat.
! References:
!   None
!-------------------------------------------------------------------------------
    implicit none

    ! Intent(in)
    real, intent(in) :: &
      Cd,       & ! Coefficient
      ubar,     &
      thlm_sfc, & ! Theta_l at zt(2) [K]
      Tsfc,     & ! Surface temperature [K]
      exner_sfc   ! Exner function at surface [-]

    ! ---- Begin Code ----

    compute_wpthlp_sfc = -Cd * ubar * ( thlm_sfc - Tsfc / exner_sfc )

    return
  end function compute_wpthlp_sfc


!===============================================================================
  real function compute_wprtp_sfc ( Cd, ubar, rtm_sfc, adjustment )

!
!  Description:
!    This function determines the surface flux of moisture.
!  References:
!    None
!-------------------------------------------------------------------------------

    implicit none

    ! Input(s)
    real, intent(in) :: &
      Cd, &
      ubar, &
      rtm_sfc, &  ! Total Water mixing ratio at zt(2) [kg/kg]
      adjustment

    ! ---- Begin Code ----

    compute_wprtp_sfc  = -Cd * ubar * ( rtm_sfc - adjustment )

    return
  end function compute_wprtp_sfc
  
!===============================================================================
  subroutine set_sclr_sfc_rtm_thlm ( wpthlp_sfc, wprtp_sfc, &
                                       wpsclrp_sfc, wpedsclrp_sfc )

!
!  Description:
!    This function determines the surface flux of moisture.
!  References:
!    None
!-------------------------------------------------------------------------------

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)
    use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl

    implicit none

    ! Input Variables
    real, intent(in) ::  & 
      wpthlp_sfc, & ! surface thetal flux        [K m/s]
      wprtp_sfc     ! surface moisture flux      [kg/kg m/s]
      
    ! Output Variables
    real, intent(out), dimension(sclr_dim) ::  & 
      wpsclrp_sfc       ! scalar surface flux            [units m/s]
    real, intent(out), dimension(edsclr_dim) ::  & 
      wpedsclrp_sfc     ! eddy-scalar surface flux       [units m/s]

    !---------------Begin Code-------------------
    
    wpsclrp_sfc(:)   = 0.0 ! Initialize flux to 0 
    wpedsclrp_sfc(:) = 0.0 ! Initialize flux to 0 

    ! Let passive scalars be equal to rt and theta_l for now
    if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
    if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

    if ( iiedsclr_thl > 0 ) wpedsclrp_sfc(iiedsclr_thl) = wpthlp_sfc
    if ( iiedsclr_rt  > 0 ) wpedsclrp_sfc(iiedsclr_rt)  = wprtp_sfc

    return
  end subroutine set_sclr_sfc_rtm_thlm


end module surface_flux
