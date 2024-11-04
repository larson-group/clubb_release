!----------------------------------------------------------------------
! $Id: twp_ice.F90 3363 2009-04-02 21:42:22Z dschanen@uwm.edu $
module twp_ice
  !
  ! Description:
  !   Contains subroutines for the Jan. 2006 TWP_ICE case.
  !
  ! References:
  !   http://users.monash.edu.au/~ladavies/gcss.html
  !----------------------------------------------------------------------

  implicit none

  public :: twp_ice_sfclyr

  private ! Default Scope

  contains

  !----------------------------------------------------------------------
  subroutine twp_ice_sfclyr( ngrdcol, time, z, exner_sfc, & 
                             thlm_sfc, ubar, rtm, p_sfc,  & 
                             saturation_formula, &
                             wpthlp_sfc, wprtp_sfc, ustar, T_sfc )
    ! Description:
    !   This subroutine computes surface fluxes of horizontal momentum,
    !   heat and moisture according to GCSS ARM specifications
    !
    ! References:
    !   http://users.monash.edu.au/~ladavies/gcss.html
    !----------------------------------------------------------------------

    use saturation, only: sat_mixrat_liq ! Procedure(s)

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    use sfc_flux, only: compute_wpthlp_sfc, compute_wprtp_sfc
    
    use time_dependent_input, only: time_sfc_given, T_sfc_given, & ! Variable(s)
                                    time_select                   ! Procedure(s)

    use interpolation, only: linear_interp_factor ! Procedure(s)

    implicit none

    intrinsic :: max, sqrt

    ! Constants
    real( kind = core_rknd ), parameter :: & 
      C_h_20  = 0.001094_core_rknd,  & ! Drag coefficient, defined by RICO 3D specification
      C_q_20  = 0.001133_core_rknd,  & ! Drag coefficient, defined by RICO 3D specification
      z0      = 0.00015_core_rknd      ! Roughness length, defined by ATEX specification

    real( kind = core_rknd ), parameter :: &
      standard_flux_alt = 20._core_rknd ! default height at which the surface flux is computed [m]

    integer, intent(in) :: &
      ngrdcol

    real(time_precision), intent(in) :: &
      time  ! current time [s]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  & 
      z,             & ! Height at zt=2      [s] 
      exner_sfc,     & ! Exner function at (2) 
      ubar,          & ! This is root (u^2 + v^2), per ATEX and RICO spec.
      thlm_sfc,      & ! thlm at (2)         [m/s]
      rtm,           & ! rt at (2)           [kg/kg]
      p_sfc             ! surface pressure    [Pa]

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    ! Output variables
    real( kind = core_rknd ), dimension(ngrdcol), intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar,        & ! surface friction velocity [m/s]
      T_sfc           ! Sea surface temp    [K]

    ! Internal variables
    real( kind = core_rknd ) :: & 
      time_frac, & ! time fraction used for interpolation
      T_sfc_interp

    real( kind = core_rknd ), dimension(ngrdcol) :: & 
        rsat, &
        Ch,   & ! This is C_h_20 scaled to the height of the lowest model level.
        Cq      ! This is C_q_20 scaled to the height of the lowest model level.

    integer :: &
      before_time, after_time, i  ! time indexes used for interpolation

    !----------------------------------------------------------------------

    !$acc enter data create( rsat, Ch, Cq )

    ! interpolate T_sfc from time_dependent_input

    call time_select( time, size(time_sfc_given), time_sfc_given, &
                       before_time, after_time, time_frac )

    T_sfc_interp = linear_interp_factor( time_frac, T_sfc_given(after_time), &
                                         T_sfc_given(before_time) )

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol

      T_sfc(i) = T_sfc_interp

      ! Declare the value of ustar.
      ustar(i) = 0.3_core_rknd

      ! Modification in case lowest model level isn't at 10 m, from ATEX specification
      !Cm   = C_m_20 * ((log(20/z0))/(log(z/z0))) * &
      !       ((log(20/z0))/(log(z/z0)))

      ! (Stevens, et al. 2000, eq 3)
      ! Modification in case lowest model level isn't at 10 m, from ATEX specification
      Ch(i)   = C_h_20 * ((log(standard_flux_alt/z0))/(log(z(i)/z0))) * & 
            ((log(standard_flux_alt/z0))/(log(z(i)/z0)))
      ! Modification in case lowest model level isn't at 10 m, from ATEX specification
      Cq(i)   = C_q_20 * ((log(standard_flux_alt/z0))/(log(z(i)/z0))) * & 
            ((log(standard_flux_alt/z0))/(log(z(i)/z0)))

      rsat(i) = sat_mixrat_liq( p_sfc(i), T_sfc(i), saturation_formula )
      !upwp_sfc   = -um_sfc * Cm * ubar  ! m^2 s^-2
      !vpwp_sfc   = -vm_sfc * Cm * ubar  ! m^2 s^-2
    end do

    ! Compute wpthlp_sfc and wprtp_sfc
    call compute_wpthlp_sfc( ngrdcol, Ch, ubar, thlm_sfc, T_sfc, exner_sfc, &
                             wpthlp_sfc ) 
  
    call compute_wprtp_sfc( ngrdcol, Cq, ubar, rtm, rsat, &
                            wprtp_sfc )

    !$acc exit data delete( rsat, Ch, Cq )

    return

  end subroutine twp_ice_sfclyr
end module twp_ice
