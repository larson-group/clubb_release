! $Id$
module clubb_model_settings  

  use clubb_precision, only: time_precision, core_rknd ! Variable(s)

  implicit none

  private ! Default scope

  ! Model settings

  ! Grid definition
  integer, public ::  & 
    nzmax = 100,   & ! Vertical extent in levels( relevant for grid type 2 and 3 only )  [#]
    grid_type = 1! 1 ==> evenly-spaced grid levels
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
    extend_atmos_bottom_level = 0, & ! Bottom level of the extended atmosphere
    extend_atmos_top_level    = 0, & ! Top level of the extended atmosphere
    extend_atmos_range_size   = 0    ! The number of levels in the extended atmosphere

!$omp threadprivate(extend_atmos_bottom_level, extend_atmos_top_level, &
!$omp               extend_atmos_range_size)

  ! The number of interpolated levels between the computational grid
  ! and the extended atmosphere
  integer, public :: &
    lin_int_buffer = 0

!$omp threadprivate(lin_int_buffer)

  real( kind = core_rknd ), public ::  & 
    deltaz  = 40._core_rknd, & ! Change in altitude per grid level     [m]
    zm_init = 0._core_rknd,  & ! Initial point on the momentum grid    [m]
    zm_top  = 1000._core_rknd  ! Maximum point on the momentum grid    [m]

!$omp threadprivate(deltaz, zm_init, zm_top)

  ! For grid_type 2 or 3 (stretched grid cases)
  character(len=100), public :: & 
    zt_grid_fname = "", & ! Path and filename of thermodynamic level altitudes
    zm_grid_fname = ""    ! Path and filename of momentum level altitudes

!$omp threadprivate(zt_grid_fname, zm_grid_fname)

  integer, public ::  & 
    day = 1, month = 1, year = 1900 ! Day of start of simulation

!$omp threadprivate(day, month, year)

  real( kind = core_rknd ), public ::  & 
    rlat = 0._core_rknd, & ! Latitude  [Degrees North]
    rlon = 0._core_rknd    ! Longitude [Degrees East]

!$omp threadprivate(rlat, rlon)

  real( kind = core_rknd ), public ::  &
    sfc_elevation = 0._core_rknd ! Elevation of ground level  [m AMSL]

!$omp threadprivate(sfc_elevation)

  character(len=50), public ::  & 
    runtype = "generic" ! String identifying the model case; e.g. bomex

!$omp threadprivate(runtype)

  integer, public :: &
    sfctype = 0 ! 0: fixed sfc sensible and latent heat fluxes as
  !                  given in namelist
  !               1: bulk formula: uses given surface temperature
  !                  and assumes over ocean

!$omp threadprivate(sfctype)

  real(kind=time_precision), public :: & 
    time_initial = 0._time_precision, &    ! Time of start of simulation     [s]
    time_final   = 3600._time_precision, & ! Time end of simulation          [s]
    time_current = 0._time_precision       ! Current time of simulation      [s]
!$omp threadprivate(time_initial, time_final, &
!$omp               time_current)

  real(kind=time_precision), public ::  & 
    dt_main = 60._time_precision, & ! Main model timestep                    [s]
    dt_rad  = 600._time_precision   ! Closure model timestep                 [s]
!$omp threadprivate(dt_main, dt_rad)

end module clubb_model_settings  
