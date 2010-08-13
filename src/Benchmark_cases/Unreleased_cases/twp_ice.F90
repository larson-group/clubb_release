!----------------------------------------------------------------------
! $Id: twp_ice.F90 3363 2009-04-02 21:42:22Z dschanen@uwm.edu $
module twp_ice
  !
  !       Description:
  !       Contains subroutines for the Jan. 2006 TWP_ICE case.
  !----------------------------------------------------------------------

  implicit none

  public :: twp_ice_sfclyr

  private ! Default Scope

  contains

  !----------------------------------------------------------------------
  subroutine twp_ice_sfclyr( z, T_sfc, exner_sfc, thlm_sfc, & 
                              ubar, rtm, p_sfc,  & 
                              wpthlp_sfc, wprtp_sfc, ustar )
    !       Description:
    !       This subroutine computes surface fluxes of horizontal momentum,
    !       heat and moisture according to GCSS ARM specifications
    !----------------------------------------------------------------------

    use constants_clubb, only: Cp, Lv, grav ! Variable(s)

    use saturation, only: sat_mixrat_liq ! Procedure(s)

    use stats_precision, only: time_precision ! Variable(s)

    use surface_flux, only: compute_wpthlp_sfc, compute_wprtp_sfc
    
    implicit none

    intrinsic :: max, sqrt

    ! Constants
    real, parameter :: & 
      C_h_20  = 0.001094,  & ! Drag coefficient, defined by RICO 3D specification
      C_q_20  = 0.001133,  & ! Drag coefficient, defined by RICO 3D specification
      z0      = 0.00015      ! Roughness length, defined by ATEX specification

    real, intent(in) ::  & 
      z,             & ! Height at zt=2      [s] 
      T_sfc,          & ! Sea surface temp    [K]
      exner_sfc,     & ! Exner function at (2) 
      ubar,          & ! This is root (u^2 + v^2), per ATEX and RICO spec.
      thlm_sfc,      & ! thlm at (2)         [m/s]
      rtm,           & ! rt at (2)           [kg/kg]
      p_sfc             ! surface pressure    [Pa]

    ! Output variables
    real, intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    ! Internal variables
    real :: & 
      Ch,   & ! This is C_h_20 scaled to the height of the lowest model level.
      Cq      ! This is C_q_20 scaled to the height of the lowest model level.
    !----------------------------------------------------------------------

    ! Declare the value of ustar.
    ustar = 0.3

    ! Modification in case lowest model level isn't at 10 m, from ATEX specification
    !Cm   = C_m_20 * ((log(20/z0))/(log(z/z0))) * &
    !       ((log(20/z0))/(log(z/z0)))
    ! Modification in case lowest model level isn't at 10 m, from ATEX specification
    Ch   = C_h_20 * ((log(20/z0))/(log(z/z0))) * & 
           ((log(20/z0))/(log(z/z0)))
    ! Modification in case lowest model level isn't at 10 m, from ATEX specification
    Cq   = C_q_20 * ((log(20/z0))/(log(z/z0))) * & 
           ((log(20/z0))/(log(z/z0)))

    wpthlp_sfc = compute_wpthlp_sfc( Ch, ubar, thlm_sfc, T_sfc, exner_sfc )
    wprtp_sfc  = compute_wprtp_sfc( Cq, ubar, rtm, sat_mixrat_liq(p_sfc,T_sfc) )
    !upwp_sfc   = -um_sfc * Cm * ubar  ! m^2 s^-2
    !vpwp_sfc   = -vm_sfc * Cm * ubar  ! m^2 s^-2

    return
  end subroutine twp_ice_sfclyr
end module twp_ice
