! $Id$
module clubb_model_settings  

  use clubb_precision, only: time_precision, core_rknd ! Variable(s)

  implicit none

  private ! Default scope

  ! Model settings

  ! Grid definition
  integer, public ::  & 
    nzmax,     & ! Vertical extent in levels( relevant for 
  !                                         grid type 2 and 3 only )  [#]
    grid_type    ! 1 ==> evenly-spaced grid levels
  !                2 ==> stretched (unevenly-spaced) grid entered on
  !                      thermodynamic grid levels; momentum levels
  !                      halfway between thermodynamic levels (style
  !                      of SAM stretched grid).
  !                3 ==> stretched (unevenly-spaced) grid entered on
  !                      momentum grid levels; thermodynamic levels
  !                      halfway between momentum levels (style
  !                      of WRF stretched grid).

  ! Radiation variables
  integer, public :: &
    extend_atmos_bottom_level, & ! Bottom level of the extended atmosphere
    extend_atmos_top_level,    & ! Top level of the extended atmosphere
    extend_atmos_range_size      ! The number of levels in the extended atmosphere

  ! The number of interpolated levels between the computational grid
  ! and the extended atmosphere
  integer, public :: &
    lin_int_buffer

! Note: Do not indent these omp directives, they must begin in the 2nd column
!$omp threadprivate(extend_atmos_bottom_level, extend_atmos_top_level)
!$omp threadprivate(extend_atmos_range_size, lin_int_buffer)

  real( kind = core_rknd ), public ::  & 
    deltaz,  & ! Change per grid level                 [m]
    zm_init, & ! Initial point on the momentum grid    [m]
    zm_top     ! Maximum point on the momentum grid    [m]

!$omp threadprivate(nzmax, grid_type, zm_init, deltaz)

  ! For grid_type 2 or 3 (stretched grid cases)
  character(len=100), public :: & 
    zt_grid_fname, & ! Path and filename of thermodynamic level altitudes
    zm_grid_fname    ! Path and filename of momentum level altitudes

!$omp threadprivate(zt_grid_fname, zm_grid_fname)

  integer, public ::  & 
    day, month, year ! Day of start of simulation

!$omp threadprivate(day, month, year)

  real( kind = core_rknd ), public ::  & 
    rlat,  & ! Latitude  [Degrees North]
    rlon     ! Longitude [Degrees East]

!$omp threadprivate(rlat, rlon)

  real( kind = core_rknd ), public ::  &
    sfc_elevation  ! Elevation of ground level  [m AMSL]

!$omp threadprivate(sfc_elevation)

  character(len=50), public ::  & 
    runtype ! String identifying the model case; e.g. bomex

!$omp threadprivate(runtype)

  integer, public :: &
    sfctype ! 0: fixed sfc sensible and latent heat fluxes as
  !              given in namelist
  !           1: bulk formula: uses given surface temperature
  !              and assumes over ocean

!$omp threadprivate(sfctype)

  real(kind=time_precision), public :: & 
    time_initial,  & ! Time of start of simulation     [s]
    time_final,    & ! Time end of simulation          [s]
    time_current     ! Current time of simulation      [s]
!$omp threadprivate(time_initial, time_final, &
!$omp               time_current)

  real(kind=time_precision), public ::  & 
    dt_main,      & ! Main model timestep                      [s]
    dt_rad          ! Closure model timestep                   [s]
!$omp threadprivate(dt_main, dt_rad)

end module clubb_model_settings  
