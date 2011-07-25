!----------------------------------------------------------------------
! $Id$
module dycoms2_rf02

  ! Description:
  !   Contains subroutines for the DYCOMS II RF02 case.
  !----------------------------------------------------------------------

  implicit none

  public :: dycoms2_rf02_tndcy, dycoms2_rf02_sfclyr

  private ! Default Scope

  contains

  !----------------------------------------------------------------------
  subroutine dycoms2_rf02_tndcy( wm_zt, wm_zm,    &
                                 thlm_forcing, rtm_forcing,  & 
                                 sclrm_forcing, edsclrm_forcing )
    ! Description:
    !   Compute thlm_ls, rtm_ls and adjust subsidence as needed.

    ! References:
    !  "A single-column model intercomparison of a heavily drizzling
    !  stratocumulus-topped boundary layer" Wyant, Matthew C., et al., (2007)
    !  J, Geophys. Res., 112, D24204,
    !  ftp://eos.atmos.washington.edu/pub/breth/papers/2007/GCSS-DYCOMS2-SCM.pdf
    !
    !  ``Dynamics and Chemistry of Marine Stratocumulus -- DYCOMS-II''
    !  Stevens, Bjorn, et al., (2003)
    !  Bull. Amer. Meteorol. Soc., 84, 579-593.
    !  http://www.atmos.ucla.edu/~bstevens/Documents/dycoms.pdf 
    !----------------------------------------------------------------------

    use grid_class, only: gr ! Variable(s)

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use array_index, only:  & 
        iisclr_thl, iisclr_rt, iiedsclr_rt, iiedsclr_thl ! Variable(s)

    implicit none

    ! Input/Output Variables
    real, dimension(gr%nnzp), intent(inout) :: &
      wm_zt, & ! W wind component at thermodynamic levels   [m/s]
      wm_zm    ! W wind component at momentum levels        [m/s]

    ! Output Variables
    real, intent(out), dimension(gr%nnzp) ::  & 
      thlm_forcing, & ! theta_l forcing                [K/s]
      rtm_forcing     ! r_t forcing                    [(kg/kg)/s] 

    real, intent(out), dimension(gr%nnzp,sclr_dim) :: & 
      sclrm_forcing    ! Passive scalar tendency        [units/s]

    real, intent(out), dimension(gr%nnzp,edsclr_dim) :: & 
      edsclrm_forcing  ! Eddy-passive scalar tendency   [units/s]

    ! ---- Begin Code ----

    ! Enter the final thlm and rtm tendency

    thlm_forcing(1:gr%nnzp) = 0.0
    rtm_forcing(1:gr%nnzp) = 0.0

    ! Imposed large-scale subsidence at the uppermost level.
    ! CLUBB used a "one-sided" derivative method to compute mean advection at
    ! the uppermost thermodynamic level.  In order to avoid bringing in large
    ! amounts of various quantities from above the top of the domain, set wm_zt
    ! to 0 at level gr%nnzp.  To stay consistent, set wm_zm to 0 at level
    ! gr%nnzp.
    wm_zt(gr%nnzp) = 0.0
    wm_zm(gr%nnzp) = 0.0

    ! Test scalars with thetal and rt if desired
    if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
    if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

    if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
    if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

    return
  end subroutine dycoms2_rf02_tndcy


!----------------------------------------------------------------------

  subroutine dycoms2_rf02_sfclyr( time, wpthlp_sfc, wprtp_sfc, ustar )
  ! Description:
  !   This subroutine computes surface fluxes of
  !   heat and moisture according to GCSS DYCOMS II RF 02 specifications

  ! References:
  !  ``Dynamics and Chemistry of Marine Stratocumulus -- DYCOMS-II''
  !  Stevens, Bjorn, et al., (2003)
  !  Bull. Amer. Meteorol. Soc., 84, 579-593.
  !  http://www.atmos.ucla.edu/~bstevens/Documents/dycoms.pdf 
  !----------------------------------------------------------------------

    use constants_clubb, only: Cp, Lv ! Variable(s)

    use surface_flux, only: convert_sens_ht_to_km_s, convert_latent_ht_to_m_s ! Procedure(s)

    use time_dependent_input, only: sens_ht_given, latent_ht_given, time_sfc_given,& ! Variable(s)
                                    time_select ! Procedure(s)

    use interpolation, only: factor_interp ! Procedure(s)

    use stats_precision, only: time_precision ! Variable(s)

    implicit none

    real(time_precision), intent(in) :: &
      time ! The current time [s]

    ! Output
    real, intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    ! External
    intrinsic :: sqrt

    ! Constant parameters
    real ::  & 
      sens_ht, &   ! Sensible heat flux
      latent_ht, &   ! Latent heat fluxi
      time_frac ! The time fraction for interpolation

    integer :: &
      before_time, after_time ! The times used for interpolation

    real, parameter :: &
      rho_sfc = 1.21 ! Air density at surface

    !------------------------BEGIN CODE-----------------------------------

    call time_select( time, size(time_sfc_given), time_sfc_given, &
                      before_time, after_time, time_frac )

    sens_ht = factor_interp( time_frac, sens_ht_given(after_time), sens_ht_given(before_time) )
    latent_ht = factor_interp( time_frac, latent_ht_given(after_time), &
                                          latent_ht_given(before_time) )

    ! Declare the value of ustar.
    ustar = 0.25

    wpthlp_sfc = convert_sens_ht_to_km_s( sens_ht, rho_sfc )
    wprtp_sfc  = convert_latent_ht_to_m_s( latent_ht, rho_sfc )

    return
  end subroutine dycoms2_rf02_sfclyr

end module dycoms2_rf02
