!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module parameters_model

! Description:
!   Contains model parameters that are determined at run time rather than
!   compile time.
!
! References:
!   None
!-------------------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd

  implicit none

  private ! Default scope

  ! Maximum magnitude of PDF parameter 'mixt_frac'. 
  real( kind = core_rknd ), public :: &
    mixt_frac_max_mag

  ! Model parameters and constraints setup in the namelists
  real( kind = core_rknd ), public ::  & 
    T0       = 300._core_rknd, & ! Reference temperature (usually 300)  [K]
    ts_nudge = 0._core_rknd      ! Timescale of u/v nudging             [s]

  real( kind = core_rknd), public :: &
    rtm_min                = epsilon( rtm_min ), & ! Value below which rtm will be nudged [kg/kg]
    rtm_nudge_max_altitude = 10000._core_rknd      ! Highest altitude at which to nudge rtm [m]

  public :: setup_parameters_model 

  contains

!-------------------------------------------------------------------------------
  subroutine setup_parameters_model ( T0_in, ts_nudge_in, Skw_max_mag )

! Description:
!   Sets parameters to their initial values
!
! References:
!   None
!-------------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    !------------------------ Input Variables ------------------------
    real( kind = core_rknd ), intent(in) ::  & 
      T0_in,       & ! Ref. temperature                         [K]
      ts_nudge_in, & ! Timescale for u/v nudging                [s]
      Skw_max_mag    ! Maximum allowable magnitude of Skewness  [-]

    !------------------------ Begin Code ------------------------
     
    ! Formula from subroutine pdf_closure, where sigma_sqd_w = 0.4 and Skw =
    ! Skw_max_mag in this formula.  Note that this is constant, but can't appear
    ! with a Fortran parameter attribute, so we define it here. 
    mixt_frac_max_mag = 1.0_core_rknd &
      - ( 0.5_core_rknd * ( 1.0_core_rknd - Skw_max_mag / &
      sqrt( 4.0_core_rknd * ( 1.0_core_rknd - 0.4_core_rknd )**3 &
      + Skw_max_mag**2 ) ) ) ! Known magic number

    T0       = T0_in
    ts_nudge = ts_nudge_in

    return
  end subroutine setup_parameters_model
!-------------------------------------------------------------------------------

end module parameters_model
