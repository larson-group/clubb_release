! $Id$
module clubb_model_settings  

  use clubb_precision, only: time_precision, core_rknd ! Variable(s)

  implicit none

  private ! Default scope

  public :: initialize_clubb_model_settings ! Initialization subroutine

  ! Model settings

  ! Grid definition
  integer, public ::  & 
    nzmax, &     ! Vertical extent in levels( relevant for grid type 2 and 3 only )  [#]
    grid_type    ! 1 ==> evenly-spaced grid levels
  !                2 ==> stretched (unevenly-spaced) grid entered on
  !                      thermodynamic grid levels; momentum levels
  !                      halfway between thermodynamic levels (style
  !                      of SAM stretched grid).
  !                3 ==> stretched (unevenly-spaced) grid entered on
  !                      momentum grid levels; thermodynamic levels
  !                      halfway between momentum levels (style
  !                      of WRF stretched grid).

! Note: Do not indent these omp directives, they must begin in the 2nd column
!$omp threadprivate(nzmax, grid_type)

  ! Radiation variables
  integer, public :: &
    extended_atmos_bottom_level, & ! Bottom level of the extended atmosphere
    extended_atmos_top_level,    & ! Top level of the extended atmosphere
    extended_atmos_range_size      ! The number of levels in the extended atmosphere

!$omp threadprivate(extended_atmos_bottom_level, extended_atmos_top_level, &
!$omp               extended_atmos_range_size)

  ! The number of interpolated levels between the computational grid
  ! and the extended atmosphere
  integer, public :: &
    lin_int_buffer

!$omp threadprivate(lin_int_buffer)

  real( kind = core_rknd ), public ::  & 
    deltaz,  & ! Change in altitude per grid level     [m]
    zm_init, & ! Initial point on the momentum grid    [m]
    zm_top     ! Maximum point on the momentum grid    [m]

!$omp threadprivate(deltaz, zm_init, zm_top)

  ! For grid_type 2 or 3 (stretched grid cases)
  character(len=100), public :: & 
    zt_grid_fname,  & ! Path and filename of thermodynamic level altitudes
    zm_grid_fname     ! Path and filename of momentum level altitudes

!$omp threadprivate(zt_grid_fname, zm_grid_fname)

  integer, public ::  & 
    day, month, year ! Start time the of simulation

!$omp threadprivate(day, month, year)

  real( kind = core_rknd ), public ::  & 
    lat_vals, & ! Latitude  [Degrees North]
    lon_vals    ! Longitude [Degrees East]

!$omp threadprivate(lat_vals, lon_vals)

  real( kind = core_rknd ), public ::  &
    sfc_elevation ! Elevation of ground level  [m AMSL]

!$omp threadprivate(sfc_elevation)

  character(len=50), public ::  & 
    runtype ! String identifying the model case; e.g. bomex

!$omp threadprivate(runtype)

  integer, public :: &
    sfctype ! 0: fixed sfc sensible and latent heat fluxes as
  !                  given in namelist
  !           1: bulk formula: uses given surface temperature
  !                  and assumes over ocean

!$omp threadprivate(sfctype)

  real(kind = time_precision ), public :: & 
    time_initial, & ! Time of start of simulation     [s]
    time_final,   & ! Time end of simulation          [s]
    time_current   !  Current time of simulation      [s]
!$omp threadprivate(time_initial, time_final, &
!$omp               time_current)

  real(kind = core_rknd ), public ::  & 
    dt_main,  & ! Main model timestep                    [s]
    dt_rad      ! Closure model timestep                 [s]
!$omp threadprivate(dt_main, dt_rad)

  integer, parameter :: &
    sp = selected_real_kind(6)  ! 32-bit floating point number

  real( kind = sp ), public :: &
    PosInf = transfer( 2139095040, 1.0_sp ) ! 2139095040 is a magic number

  contains

!-------------------------------------------------------------------------------
  subroutine initialize_clubb_model_settings( )

! Description:
!   Sets all variables to a default setting.

! References:
!   None
!-------------------------------------------------------------------------------

    implicit none

    nzmax     = 75
    grid_type = 1

    extended_atmos_bottom_level = 0
    extended_atmos_top_level    = 0
    extended_atmos_range_size   = 0

    lin_int_buffer = 0

    deltaz  = 40._core_rknd
    zm_init = 0._core_rknd
    zm_top  = 3500._core_rknd

    zt_grid_fname = ""
    zm_grid_fname = ""

    day = 22; month = 6; year = 1969

    lat_vals = 15._core_rknd
    lon_vals = -56.5_core_rknd

    sfc_elevation = 0._core_rknd

    runtype = "bomex"

    sfctype = 0 

    time_initial = 0._time_precision
    time_final   = 21600._time_precision
    time_current = 0._time_precision

    dt_main = 60._core_rknd
    dt_rad  = 600._core_rknd

    return
  end subroutine initialize_clubb_model_settings

end module clubb_model_settings  
