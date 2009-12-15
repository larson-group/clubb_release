!-----------------------------------------------------------------------
! $Id$

module clubb_driver

! Description:
!   Contains the necessary subroutines to execute individual CLUBB
!   model runs, using one of the driver programs (the simplest case
!   being the clubb_standalone program).
!-----------------------------------------------------------------------

  use stats_precision, only: time_precision ! Variable(s)
  use text_writer, only: write_text, write_date

  implicit none

  ! Setup run_clubb() as the sole external interface
  private ::  & 
    initialize_clubb, & 
    advance_clubb_forcings, & 
    restart_clubb

  public :: &
    run_clubb

  private ! Default to private

  ! Model settings

  ! Grid definition
  integer, private ::  & 
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
  integer, private ::  & 
    extend_atmos_bottom_level, & ! Bottom level of the extended atmosphere
    extend_atmos_top_level,    & ! Top level of the extended atmosphere
    extend_atmos_range_size,   & ! The number of levels in the extended
                                 ! atmosphere
    lin_int_buffer               ! The number of interpolated levels between
                                 ! the computational grid and the extended
                                 ! atmosphere
!$omp threadprivate(extend_atmos_bottom_level, extend_atmos_top_level)
!$omp threadprivate(extend_atmos_range_size, lin_int_buffer)
  
  real, private ::  & 
    deltaz,  & ! Change per grid level                 [m]
    zm_init, & ! Initial point on the momentum grid    [m]
    zm_top     ! Maximum point on the momentum grid    [m]

! Do not indent these omp directives, they must begin in the 2nd column
!$omp threadprivate(nzmax, grid_type, zm_init, deltaz)

  ! For grid_type 2 or 3 (stretched grid cases)
  character(len=100), private :: & 
    zt_grid_fname, & ! Path and filename of thermodynamic level altitudes
    zm_grid_fname    ! Path and filename of momentum level altitudes

!$omp threadprivate(zt_grid_fname, zm_grid_fname)

  integer, private ::  & 
    day, month, year ! Day of start of simulation

!$omp threadprivate(day, month, year)

  real, private ::  & 
    rlat,  & ! Latitude  [Degrees North]
    rlon     ! Longitude [Degrees East]

!$omp threadprivate(rlat, rlon)

  real, private ::  &
    sfc_elevation  ! Elevation of ground level  [m AMSL]

!$omp threadprivate(sfc_elevation)

  character(len=50), private ::  & 
    runtype ! String identifying the model case; e.g. bomex

!$omp threadprivate(runtype)

  integer, private :: &
    sfctype ! 0: fixed sfc sensible and latent heat fluxes as
  !              given in namelist
  !           1: bulk formula: uses given surface temperature
  !              and assumes over ocean

!$omp threadprivate(sfctype)

  real(kind=time_precision), private :: & 
    time_initial,  & ! Time of start of simulation     [s]
    time_final,    & ! Time end of simulation          [s]
    time_spinup,   & ! Time end of spin up period      [s]
    time_current     ! Current time of simulation      [s]
!$omp threadprivate(time_initial, time_final, time_spinup, &
!$omp               time_current)

  real(kind=time_precision), private ::  & 
    dtmain,      & ! Main model timestep                      [s]
    dtclosure,   & ! Closure model timestep                   [s]
    dt             ! Current model timestep (based on spinup) [s]
!$omp threadprivate(dtmain, dtclosure, dt)

  contains

  !-----------------------------------------------------------------------
  subroutine run_clubb & 
             ( params, runfile, err_code, l_stdout )
    ! Description:
    !   Subprogram to integrate the pde equations for pdf closure.

    ! References:
    !   None
    !-----------------------------------------------------------------------

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: read_grid_heights, zt2zm, zm2zt ! Procedure(s)

    use parameter_indices, only: nparams ! Variable(s)

    use variables_diagnostic_module, only: ug, vg, em,  & ! Variable(s)
      tau_zt, thvm, Lscale, Kh_zm, & 
      um_ref, vm_ref, Ncnm, wp2_zt, &
      hydromet, thlm_ref, rtm_ref, &
      Frad, radht, Frad_SW_up, &
      Frad_LW_up, Frad_SW_down, Frad_LW_down

    use variables_prognostic_module, only:  & 
      Tsfc, psfc, SE, LE, thlm, rtm,     & ! Variable(s)
      um, vm, wp2, rcm, wm_zt, wm_zm, exner, & 
      tau_zm, p_in_Pa, rho_zm, upwp, vpwp, wpthlp, wpthvp, & 
      wprcp, Kh_zt, rho, wprtp, wpthlp_sfc, wprtp_sfc, & 
      upwp_sfc, vpwp_sfc, thlm_forcing, & 
      rtm_forcing, um_forcing, vm_forcing, &
      up2, vp2, wp3, rtp2, pdf_params, & 
      thlp2, rtpthlp, sigma_sqd_w, cloud_frac, &
      rcm_in_layer, cloud_cover

    use variables_prognostic_module, only:  & 
      sclrm, sclrp2, sclrprtp, sclrpthlp, sclrm_forcing, & ! Variables
      wpsclrp, wpsclrp_sfc,  & 
      edsclrm, edsclrm_forcing, wpedsclrp_sfc


    use numerical_check, only: invalid_model_arrays ! Procedure(s)

    use inputfields, only: &
      inputfields_init, compute_timestep, stat_fields_reader ! Procedure(s)

    use inputfields, only: stat_file_zt

    use parameters_tunable, only: params_list ! Variable(s)

    use clubb_core, only: & 
      setup_clubb_core,  & ! Procedure(s) 
      cleanup_clubb_core, & 
      advance_clubb_core

    use constants, only: fstdout, fstderr, & ! Variable(s)
      rttol, thltol, wtol, wtol_sqd

    use error_code, only: clubb_var_out_of_bounds,  & ! Variable(s)
      clubb_var_equals_NaN, & 
      clubb_rtm_level_not_found, &
      clubb_at_least_debug_level ! Function

    use error_code, only: fatal_error,  & ! Procedure(s)
                          set_clubb_debug_level

    use stats_precision, only: time_precision ! Variable(s)

    use array_index, only: iisclr_rt, iisclr_thl, iisclr_CO2, & ! Variables
      iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2

    use microphys_driver, only: init_microphys ! Subroutine

    use parameters_radiation, only: init_radiation ! Subroutine

    use model_flags, only: & 
      l_pos_def, l_hole_fill, & ! Constants
      l_single_C2_Skw, l_gamma_Skw, l_byteswap_io

    use stats_variables, only: l_stats_last, l_stats_samp, & ! Variable(s)
      l_output_rad_files
    use stats_variables, only: zt ! Type

    use stats_variables, only: &
      irtm_mc, & ! Variables
      irvm_mc, &
      ircm_mc, &
      ithlm_mc

    use stats_subs, only:  & 
      stats_begin_timestep, stats_end_timestep,  & ! Procedure(s)
      stats_finalize, stats_init

    use stats_type, only: &
      stat_update_var ! Procedure

    use sounding, only: sclr_max ! Variable(s)

    use time_dependent_input, only: &
      l_t_dependent,    & ! Variable(s)
      l_input_xpwp_sfc, &
      finalize_t_dependent_input ! Procedure(s)

    use sponge_layer_damping, only: &
      thlm_sponge_damp_settings, &
      rtm_sponge_damp_settings, &
      uv_sponge_damp_settings, &
      thlm_sponge_damp_profile, &
      rtm_sponge_damp_profile, &
      uv_sponge_damp_profile, &
      finalize_tau_sponge_damp

    use extend_atmosphere_mod, only: &
      total_atmos_dim, &
      complete_alt, &
      complete_momentum, &
      l_use_default_std_atmosphere, &
      finalize_extend_atm

    use parameters_radiation, only: rad_scheme, l_bugsrad ! Variable(s)

    use parameters_microphys, only: &
      l_latin_hypercube_sampling ! Variable

#ifdef UNRELEASED_CODE
    use latin_hypercube_mod, only: &
      latin_hypercube_2D_output, & ! Procedure(s)
      latin_hypercube_2D_close
#endif

    implicit none

    ! Because Fortran I/O is not thread safe, we use this here to
    ! insure that no model uses the same file number simultaneously
    ! when doing a tuning run. -dschanen 31 Jan 2007
#ifdef _OPENMP
    integer :: omp_get_thread_num ! Function
#endif
    ! External
    intrinsic :: mod, real, int, trim, floor, max

    ! Constant Parameters
    logical, parameter :: &
      l_host_applies_sfc_fluxes = .false.

    ! Input Variables
    logical, intent(in) ::  & 
      l_stdout   ! Whether to print output per timestep

    real, intent(in), dimension(nparams) ::  & 
      params  ! Model parameters, C1, nu2, etc.

    ! Subroutine Arguments (Model Setting)
    character(len=*), intent(in) ::  & 
      runfile ! Name of file containing &model_setting and &sounding

    ! Output Variables
    integer, intent(inout) :: &
      err_code ! An error code is returned indicating the status of the run. See error_code.F90

    ! Local Variables
    ! Internal Timing Variables
    integer :: & 
      ifinal, & 
      niterlong

    integer :: & 
      debug_level     ! Amount of debugging information

    real ::  & 
      fcor,            & ! Coriolis parameter            [s^-1]
      T0,              & ! Reference Temperature         [K]
      ts_nudge           ! Timescale for u/v nudging     [s]

    real, dimension(sclr_max) :: & 
      sclr_tol        ! Thresholds on the passive scalars     [units vary]

    real(kind=time_precision) :: & 
      time_restart    ! Time of model restart run     [s]

    logical ::  & 
      l_soil_veg,     & ! Flag for simple surface scheme
      l_uv_nudge,     & ! Whether to adjust the winds within the timestep
      l_restart,      & ! Flag for restarting from GrADS file
      l_input_fields, & ! Whether to set model variables from a file
      l_tke_aniso       ! For anisotropic turbulent kinetic energy,
                        ! i.e. TKE = 1/2 (u'^2 + v'^2 + w'^2)

    character(len=6) :: &
      saturation_formula ! "bolton" approx. or "flatau" approx.

    character(len=128) ::  & 
      restart_path_case,   & ! GrADS file used in case of restart
      forcings_file_path     ! Path to the forcing files

    logical :: & 
      l_stats ! Whether statistics are computed and output to disk

    character(len=10) :: & 
      stats_fmt  ! File format for stats; typically GrADS.

    character(len=100) :: & 
      fname_prefix, &! Prefix of stats filenames, to be followed by, for example "_zt"
      fdir           ! Output directory

    real(kind=time_precision) :: & 
      stats_tsamp,   & ! Stats sampling interval [s]
      stats_tout       ! Stats output interval   [s]

    ! Grid altitude arrays
    real, dimension(:), allocatable ::  & 
      momentum_heights, thermodynamic_heights ! [m]

    ! Dummy dx and dy horizontal grid spacing.
    real :: dummy_dx, dummy_dy  ! [m]

    integer :: i, i1, j ! Internal Loop Variables
    integer :: iinit ! initial iteration

    integer ::  & 
      iunit,           & ! File unit used for I/O
      hydromet_dim,    & ! Number of hydrometeor species        [#]
      sclr_dim,        & ! Number of passive scalars            [#]
      edsclr_dim         ! Number of passive scalars            [#]

    integer :: itime_nearest ! Used for and inputfields run [s]

    logical, parameter :: &
      l_write_to_file = .true. ! If true, will write case information to a file
    character(len=150) :: &
      case_info_file ! The filename for case info

    real, allocatable, dimension(:) :: &
      rcm_mc, & ! Tendency of liquid water due to microphysics     [kg/kg/s]
      rvm_mc, & ! Tendency of vapor water due to microphysics      [kg/kg/s]
      thlm_mc   ! Tendecy of liquid pot. temp. due to microphysics [K/s]

    ! Definition of namelists
    namelist /model_setting/  & 
      runtype, nzmax, grid_type, deltaz, zm_init, zm_top, & 
      zt_grid_fname, zm_grid_fname,  & 
      day, month, year, rlat, rlon, sfc_elevation, & 
      time_initial, time_final, time_spinup, & 
      dtmain, dtclosure, & 
      sfctype, Tsfc, psfc, SE, LE, fcor, T0, ts_nudge, & 
      forcings_file_path, l_t_dependent, l_input_xpwp_sfc, &
      l_use_default_std_atmosphere, &
      saturation_formula, &
      thlm_sponge_damp_settings, rtm_sponge_damp_settings, uv_sponge_damp_settings, &
      l_soil_veg, l_tke_aniso, l_uv_nudge, l_restart, restart_path_case, & 
      time_restart, l_input_fields, debug_level, & 
      sclr_tol, sclr_dim, iisclr_thl, iisclr_rt, iisclr_CO2, &
      edsclr_dim, iiedsclr_thl, iiedsclr_rt, iiedsclr_CO2


    namelist /stats_setting/ & 
      l_stats, fname_prefix, stats_tsamp, stats_tout, stats_fmt

!-----------------------------------------------------------------------

    ! Initialize the model run

    ! Pick some default values for model_setting
    runtype   = "generic"
    nzmax     = 100
    grid_type = 1
    deltaz    = 40.
    zm_init   = 0.
    zm_top    = 1000.
    zt_grid_fname = ''
    zm_grid_fname = ''

    day   = 1
    month = 1
    year  = 1900

    rlat = 0.
    rlon = 0.

    sfc_elevation = 0.

    time_initial = 0.
    time_final   = 3600.
    time_spinup  = 0.

    dtmain    = 30.
    dtclosure = 30.

    sfctype  = 0
    tsfc     = 288.
    psfc     = 1000.e2
    SE       = 0.
    LE       = 0.
    fcor     = 1.e-4
    T0       = 300.
    ts_nudge = 86400.

    forcings_file_path = ''

    l_t_dependent   = .false.
    l_input_xpwp_sfc = .false.

    l_use_default_std_atmosphere = .true.

    thlm_sponge_damp_settings%l_sponge_damping = .false.
    rtm_sponge_damp_settings%l_sponge_damping = .false.
    uv_sponge_damp_settings%l_sponge_damping = .false.

    thlm_sponge_damp_settings%tau_sponge_damp_min = 60.
    thlm_sponge_damp_settings%tau_sponge_damp_max = 1800.
    thlm_sponge_damp_settings%sponge_damp_depth = 0.25

    rtm_sponge_damp_settings%tau_sponge_damp_min = 60.
    rtm_sponge_damp_settings%tau_sponge_damp_max = 1800.
    rtm_sponge_damp_settings%sponge_damp_depth = 0.25

    uv_sponge_damp_settings%tau_sponge_damp_min = 60.
    uv_sponge_damp_settings%tau_sponge_damp_max = 1800.
    uv_sponge_damp_settings%sponge_damp_depth = 0.25

    l_soil_veg     = .false.
    l_tke_aniso    = .false.
    l_uv_nudge     = .false.
    l_restart      = .false.
    l_input_fields  = .false.
    restart_path_case = "none"
    time_restart  = 0.
    debug_level   = 2

    saturation_formula = "flatau" ! Flatau polynomial approx.

    sclr_dim   = 0
    iisclr_thl = -1
    iisclr_rt  = -1
    iisclr_CO2 = -1

    edsclr_dim = 0
    iiedsclr_thl = -1
    iiedsclr_rt  = -1
    iiedsclr_CO2 = -1

    sclr_tol(1:sclr_max) = 1.e-2

    ! Pick some default values for stats_setting
    l_stats       = .false.
    fname_prefix = ''
    stats_fmt    = ''
    stats_tsamp  = 0.
    stats_tout   = 0.

    ! Figure out which I/O unit to use for OpenMP runs
#ifdef _OPENMP
    iunit = omp_get_thread_num( ) + 10
#else
    iunit = 10
#endif

    ! Read namelist file
    open(unit=iunit, file=runfile, status='old')
    read(unit=iunit, nml=model_setting)
    read(unit=iunit, nml=stats_setting)
    close(unit=iunit)

    case_info_file = &
      "../output/" // trim( fname_prefix ) // "_setup.txt" ! The filename for case setup

    ! Sanity check on passive scalars
    ! When adding new 'ii' scalar indices, add them to this list.
    if ( max( iisclr_CO2, iisclr_rt, iisclr_thl ) > sclr_dim ) then
      write(fstderr,*) "Passive scalar index exceeds sclr_dim ", & 
        "iisclr_CO2 = ", iisclr_CO2, "iisclr_rt = ", iisclr_rt,  & 
        "iisclr_thl = ", iisclr_thl, "sclr_dim = ", sclr_dim

      err_code = clubb_var_out_of_bounds
      return

    else if ( max( iiedsclr_CO2, iiedsclr_rt, iiedsclr_thl ) > edsclr_dim ) then
      write(fstderr,*) "Passive scalar index exceeds edsclr_dim ", & 
        "iiedsclr_CO2 = ", iiedsclr_CO2, "iiedsclr_rt = ", iiedsclr_rt,  & 
        "iiedsclr_thl = ", iiedsclr_thl, "edsclr_dim = ", edsclr_dim

      err_code = clubb_var_out_of_bounds
      return
    end if

    ! Set debug level
    call set_clubb_debug_level( debug_level ) ! Intent(in)

    ! Printing Model Inputs
    if ( clubb_at_least_debug_level( 1 ) ) then

      if( l_write_to_file ) open(unit=iunit, file=case_info_file, status='replace', action='write')

      ! Print the date and time
      call write_date( l_write_to_file, iunit )

      ! Print the list of parameters that are being used before the run.
      call write_text( "Parameter          Value", l_write_to_file, iunit, '(4x,A24)')
      call write_text( "---------          -----", l_write_to_file, iunit, '(4x,A24)')
      do j = 1, nparams, 1
        call write_text(params_list(j) // " = ", params(j), & 
          l_write_to_file, iunit, '(A18,F27.20)')
      end do

      call write_text( "--------------------------------------------------", &
        l_write_to_file, iunit )
      call write_text( "Model Settings", l_write_to_file, iunit)
      call write_text( "--------------------------------------------------", &
        l_write_to_file, iunit )
      call write_text( "Preprocessing directives:", l_write_to_file, iunit )
#ifdef NETCDF
      call write_text( "-DNETCDF enabled", l_write_to_file, iunit )
#else
      call write_text( "-DNETCDF disabled", l_write_to_file, iunit )
#endif
#ifdef COAMPS_MICRO
      call write_text( "-DCOAMPS_MICRO enabled", l_write_to_file, iunit )
#else
      call write_text( "-DCOAMPS_MICRO disabled", l_write_to_file, iunit )
#endif
#ifdef TUNER
      call write_text( "-DTUNER enabled", l_write_to_file, iunit )
#else
      call write_text( "-DTUNER disabled", l_write_to_file, iunit )
#endif
#ifdef UNRELEASED_CODE
      call write_text( "-DUNRELEASED_CODE enabled", l_write_to_file, iunit )
#else
      call write_text( "-DUNRELEASED_CODE disabled", l_write_to_file, iunit )
#endif
#ifdef nooverlap
      call write_text( "-Dnooverlap enabled", l_write_to_file, iunit )
#else
      call write_text( "-Dnooverlap disabled", l_write_to_file, iunit )
#endif
#ifdef radoffline
      call write_text( "-Draoffline enabled", l_write_to_file, iunit )
#else
      call write_text( "-Dradoffline disabled", l_write_to_file, iunit )
#endif
#ifdef USE_BUGSrad_ocast_random
      call write_text( "-DUSE_BUGSrad_ocast_random enabled", l_write_to_file, iunit )
#else
      call write_text( "-DUSE_BUGSrad_ocast_random disabled", l_write_to_file, iunit )
#endif

      ! Pick some default values for model_setting
      call write_text( "&model_setting:", l_write_to_file, iunit )
      call write_text( "runtype = " // runtype, l_write_to_file, iunit )
      call write_text( "nzmax = ", nzmax, l_write_to_file, iunit )
      call write_text( "grid_type = ", grid_type, l_write_to_file, iunit )
      call write_text( "deltaz = ", deltaz, l_write_to_file, iunit )
      call write_text( "zm_init = ", zm_init, l_write_to_file, iunit )
      call write_text( "zm_top = ", zm_top, l_write_to_file, iunit )
      call write_text( "zt_grid_fname = " // zt_grid_fname, l_write_to_file, iunit )
      call write_text( "zm_grid_fname = " // zm_grid_fname, l_write_to_file, iunit )

      call write_text( "day = ", day, l_write_to_file, iunit )
      call write_text( "month = ", month, l_write_to_file, iunit )
      call write_text( "year = ", year, l_write_to_file, iunit )

      call write_text( "rlat = ", rlat, l_write_to_file, iunit )
      call write_text( "rlon = ", rlon, l_write_to_file, iunit )

      call write_text( "sfc_elevation = ", sfc_elevation, l_write_to_file, iunit )

      call write_text( "time_initial = ", real( time_initial ), l_write_to_file, iunit )
      call write_text( "time_final = ", real( time_final ), l_write_to_file, iunit )
      call write_text( "time_spinup = ", real( time_spinup ), l_write_to_file, iunit )

      call write_text( "dtmain = ", real( dtmain ), l_write_to_file, iunit )
      call write_text( "dtclosure = ", real( dtclosure ), l_write_to_file, iunit )

      call write_text( "sfctype = ", sfctype, l_write_to_file, iunit )
      call write_text( "tsfc = ", tsfc, l_write_to_file, iunit )
      call write_text( "psfc = ", psfc, l_write_to_file, iunit )
      call write_text( "SE = ", SE, l_write_to_file, iunit )
      call write_text( "LE = ", LE, l_write_to_file, iunit )
      call write_text( "fcor = ", fcor, l_write_to_file, iunit )
      call write_text( "T0 = ", T0, l_write_to_file, iunit )
      call write_text( "ts_nudge = ", ts_nudge, l_write_to_file, iunit )

      call write_text( "forcings_file_path = " // forcings_file_path, l_write_to_file, iunit )

      call write_text( "l_t_dependent = ", l_t_dependent, l_write_to_file, iunit )
      call write_text( "l_input_xpwp_sfc = ", l_input_xpwp_sfc, l_write_to_file, iunit )

      call write_text( "l_use_default_std_atmosphere = ", l_use_default_std_atmosphere, &
        l_write_to_file, iunit )

      call write_text( "saturation_formula = " // saturation_formula, &
        l_write_to_file, iunit )

      call write_text( "thlm_sponge_damp_settings%l_sponge_damping = ", & 
        thlm_sponge_damp_settings%l_sponge_damping, l_write_to_file, iunit )
      call write_text( "rtm_sponge_damp_settings%l_sponge_damping = ", &
        rtm_sponge_damp_settings%l_sponge_damping, l_write_to_file, iunit )

      call write_text( "uv_sponge_damp_settings%l_sponge_damping = ", &
        uv_sponge_damp_settings%l_sponge_damping, l_write_to_file, iunit )

      call write_text( "thlm_sponge_damp_settings%tau_sponge_damp_min = ", &
        thlm_sponge_damp_settings%tau_sponge_damp_min, l_write_to_file, iunit )
                          
      call write_text( "thlm_sponge_damp_settings%tau_sponge_damp_max = ", &
        thlm_sponge_damp_settings%tau_sponge_damp_max, l_write_to_file, iunit )
                          
      call write_text( "thlm_sponge_damp_settings%sponge_damp_depth = ", &
        thlm_sponge_damp_settings%sponge_damp_depth, l_write_to_file, iunit )
                          
      call write_text( "rtm_sponge_damp_settings%tau_sponge_damp_min = ", &
        rtm_sponge_damp_settings%tau_sponge_damp_min, l_write_to_file, iunit )
                          
      call write_text( "rtm_sponge_damp_settings%tau_sponge_damp_max = ", &
        rtm_sponge_damp_settings%tau_sponge_damp_max, l_write_to_file, iunit )

      call write_text( "rtm_sponge_damp_settings%sponge_damp_depth = ", &
        rtm_sponge_damp_settings%sponge_damp_depth, l_write_to_file, iunit )

      call write_text( "uv_sponge_damp_settings%tau_sponge_damp_min = ", &
        uv_sponge_damp_settings%tau_sponge_damp_min, l_write_to_file, iunit )

      call write_text( "uv_sponge_damp_settings%tau_sponge_damp_max = ", &
        uv_sponge_damp_settings%tau_sponge_damp_max, l_write_to_file, iunit )
                          
      call write_text( "uv_sponge_damp_settings%sponge_damp_depth = ", &
        uv_sponge_damp_settings%sponge_damp_depth, l_write_to_file, iunit )

      call write_text( "l_soil_veg = ", l_soil_veg, l_write_to_file, iunit )
      call write_text( "l_tke_aniso = ", l_tke_aniso, l_write_to_file, iunit )
      call write_text( "l_uv_nudge = ", l_uv_nudge, l_write_to_file, iunit )
      call write_text( "l_restart = ", l_restart, l_write_to_file, iunit )
      call write_text( "l_input_fields = ", l_input_fields, l_write_to_file, iunit )
      call write_text( "restart_path_case = " // restart_path_case, l_write_to_file, iunit )
      call write_text( "time_restart = ", real( time_restart ), l_write_to_file, iunit )
      call write_text( "debug_level = ", debug_level, l_write_to_file, iunit )

      call write_text( "sclr_dim = ", sclr_dim, l_write_to_file, iunit )
      call write_text( "edsclr_dim = ", edsclr_dim, l_write_to_file, iunit )
      call write_text( "iisclr_thl = ", iisclr_thl, l_write_to_file, iunit )
      call write_text( "iisclr_rt = ", iisclr_rt, l_write_to_file, iunit )
      call write_text( "iisclr_CO2 = ", iisclr_CO2, l_write_to_file, iunit )

      call write_text( "sclr_tol = ", sclr_tol(1:sclr_dim), l_write_to_file, iunit )

      ! Pick some default values for stats_setting
      call write_text( "&stats_setting:", l_write_to_file, iunit )
      call write_text( "l_stats = ", l_stats, l_write_to_file, iunit )
      call write_text( "fname_prefix = " // fname_prefix, l_write_to_file, iunit )
      call write_text( "stats_fmt = " // stats_fmt, l_write_to_file, iunit)
      call write_text( "stats_tsamp = ", real( stats_tsamp ), l_write_to_file, iunit )
      call write_text( "stats_tout = ", real( stats_tout ), l_write_to_file, iunit )

      call write_text( "Constant flags:", l_write_to_file, iunit )
      call write_text( "l_pos_def = ", l_pos_def, l_write_to_file, iunit )
      call write_text( "l_hole_fill = ", l_hole_fill, l_write_to_file, iunit )
      call write_text( "l_single_C2_Skw = ", l_single_C2_Skw, l_write_to_file, iunit )
      call write_text( "l_gamma_Skw = ", l_gamma_Skw, l_write_to_file, iunit)
      call write_text( "l_byteswap_io = ", l_byteswap_io, l_write_to_file, iunit )

      call write_text( "Constant tolerances [units]", l_write_to_file, iunit )
      call write_text( "rttol [kg/kg] = ", rttol, l_write_to_file, iunit )
      call write_text( "thltol [K] = ", thltol, l_write_to_file, iunit )
      call write_text( "wtol [m/s] = ", wtol, l_write_to_file, iunit )

      call write_text( "--------------------------------------------------", &
        l_write_to_file, iunit )

      if( l_write_to_file) close(unit=iunit);

    end if ! clubb_at_least_debug_level(1)

    !----------------------------------------------------------------------

    ! Allocate stretched grid altitude arrays.
    allocate( momentum_heights(1:nzmax),  & 
              thermodynamic_heights(1:nzmax) )

    ! Handle the reading of grid altitudes for
    ! stretched (unevenly-spaced) grid options.
    ! Do some simple error checking for all grid options.
    call read_grid_heights( nzmax, grid_type, &                 ! Intent(in)
                            zm_grid_fname, zt_grid_fname, &     ! Intent(in)
                            iunit, &                            ! Intent(in) 
                            momentum_heights, &                 ! Intent(out)
                            thermodynamic_heights )             ! Intent(out)

    ! Dummy horizontal grid spacing variables.
    dummy_dx = 0.0
    dummy_dy = 0.0

    ! Setup microphysical fields
    call init_microphys( iunit, runfile, & ! Intent(in)
                         hydromet_dim )    ! Intent(out)

    ! Setup radiation parameters
    call init_radiation( iunit, runfile ) ! Intent(in)

    ! Setup filenames and variables to set for setfields, if enabled
    if ( l_input_fields ) then
      call inputfields_init( iunit, runfile ) ! Intent(in)
    end if

    ! Allocate & initialize variables,
    ! setup grid, setup constants, and setup flags

    call setup_clubb_core                               &
         ( nzmax, T0, ts_nudge,                         & ! Intent(in)
           hydromet_dim, sclr_dim,                      & ! Intent(in)
           sclr_tol(1:sclr_dim), edsclr_dim, params,    & ! Intent(in)
           l_soil_veg, l_host_applies_sfc_fluxes,       & ! Intent(in)
           l_uv_nudge, l_tke_aniso, saturation_formula, & ! Intent(in)
           .false., grid_type, deltaz, zm_init, zm_top, & ! Intent(in)
           momentum_heights, thermodynamic_heights,     & ! Intent(in)
           dummy_dx, dummy_dy, sfc_elevation,           & ! Intent(in)
           err_code )                                     ! Intent(out)


    if ( err_code == clubb_var_out_of_bounds ) return


    ! Deallocate stretched grid altitude arrays
    deallocate( momentum_heights, thermodynamic_heights )

    ! Allocate rvm_mc, rcm_mc, thlm_mc
    allocate( rvm_mc(gr%nnzp), rcm_mc(gr%nnzp), thlm_mc(gr%nnzp) )

    ! Initialize to 0.0
    rvm_mc  = 0.0
    rcm_mc  = 0.0
    thlm_mc = 0.0

    if ( .not. l_restart ) then

      time_current = time_initial
      iinit = 1

      call initialize_clubb &
           ( iunit, trim( forcings_file_path ), psfc, zm_init, & ! Intent(in)
             thlm, rtm, um, vm, ug, vg, wp2, up2, vp2, rcm,    & ! Intent(inout)
             wm_zt, wm_zm, em, exner,                          & ! Intent(inout)
             tau_zt, tau_zm, thvm, p_in_Pa,                    & ! Intent(inout)
             rho, rho_zm, Lscale, rtm_ref, thlm_ref,           & ! Intent(inout) 
             Kh_zt, Kh_zm, um_ref, vm_ref,                     & ! Intent(inout)
             hydromet, Ncnm,                                   & ! Intent(inout)
             sclrm, edsclrm )                                    ! Intent(out)

    else  ! restart


      ! Currently initialize_clubb does more than just read in the initial sounding.
      ! It also includes other important initializations such as um_ref and vm_ref.
      ! There for it should be executed prior to a restart. The restart should overwrite
      ! the initial sounding anyway.
      call initialize_clubb &
           ( iunit, trim( forcings_file_path ), psfc, zm_init, & ! Intent(in)
             thlm, rtm, um, vm, ug, vg, wp2, up2, vp2, rcm,    & ! Intent(inout)
             wm_zt, wm_zm, em, exner,                          & ! Intent(inout)
             tau_zt, tau_zm, thvm, p_in_Pa,                    & ! Intent(inout)
             rho, rho_zm, Lscale, rtm_ref, thlm_ref,           & ! Intent(inout) 
             Kh_zt, Kh_zm, um_ref, vm_ref,                     & ! Intent(inout)
             hydromet, Ncnm,                                   & ! Intent(inout)
             sclrm, edsclrm )                                    ! Intent(out)


      time_current = time_restart

      ! Determining what iteration to restart at.
      ! The value is increased by 1 to sychronize with restart data.
      ! Joshua Fasching February 2008

      ! Ensure that iteration num, iinit, is an integer, so that model time is
      !   incremented correctly by iteration number at end of timestep
      if ( mod( (time_restart-time_initial) , dtmain ) /= 0 ) then

        print*, "Error: (time_restart-time_initial) ",  & 
          "is not a multiple of dtmain."
        print*, "time_restart = ", time_restart
        print*, "time_initial = ", time_initial
        print*, "dtmain = ", dtmain
        stop "Fatal error"

      end if ! mod( (time_restart-time_initial) , dtmain ) /= 0

      iinit = floor( ( time_current - time_initial ) / dtmain ) + 1

      call restart_clubb &
           ( iunit, runfile,                  &            ! Intent(in)
             restart_path_case, time_restart, &            ! Intent(in)
             upwp, vpwp, wm_zt, wm_zm,        &            ! Intent(inout)
             wpthlp, wprtp,   &                            ! Intent(inout)
             wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc )   ! Intent(out)

    end if ! ~l_restart



#ifdef _OPENMP
    iunit = omp_get_thread_num( ) + 50
#else
    iunit = 50
#endif

    fdir = "../output/" ! Output directory

    if ( trim( rad_scheme ) == "bugsrad" ) then
        l_bugsrad = .true.
    else
        l_bugsrad = .false.
    end if

    ! Only output radiation files if using a radiation scheme
    l_output_rad_files = l_bugsrad

    ! This is a kludge added because the grid used by BUGSrad does
    ! not include CLUBB's ghost point. -nielsenb 20 Oct 2009
    if ( l_bugsrad ) then
        ! Initialize statistics output
        call stats_init( iunit, fname_prefix, fdir, l_stats, & ! Intent(in)
                         stats_fmt, stats_tsamp, stats_tout, runfile, & ! Intent(in)
                         gr%nnzp, gr%zt, gr%zm, total_atmos_dim - 1, & ! Intent(in)
                         complete_alt(2:total_atmos_dim), total_atmos_dim, & ! Intent(in)
                         complete_momentum(2:total_atmos_dim + 1), & ! Intent(in)
                         day, month, year, & ! Intent(in)
                         (/rlat/), (/rlon/), time_current, dtmain ) ! Intent(in)
    else  
        ! Initialize statistics output
        call stats_init( iunit, fname_prefix, fdir, l_stats, & ! Intent(in)
                         stats_fmt, stats_tsamp, stats_tout, runfile, & ! Intent(in)
                         gr%nnzp, gr%zt, gr%zm, total_atmos_dim, & ! Intent(in)
                         complete_alt, total_atmos_dim + 1, complete_momentum, & ! Intent(in)
                         day, month, year, & ! Intent(in)
                         (/rlat/), (/rlon/), time_current, dtmain ) ! Intent(in)
    end if

#ifdef UNRELEASED_CODE
    if ( l_latin_hypercube_sampling ) then
      call latin_hypercube_2D_output &
           ( fname_prefix, fdir, stats_tout, gr%nnzp, &
             gr%zt, time_initial  )
    end if
#endif /*UNRELEASED_CODE*/

    ! Time integration
    ! Call advance_clubb_core once per each statistics output time
    ifinal = floor( ( time_final - time_initial ) / dtmain )


!-------------------------------------------------------------------------------
!                         Main Time Stepping Loop
!-------------------------------------------------------------------------------

    do i = iinit, ifinal, 1

      ! When this time step is over, the time will be time + dtmain

      ! We use elapsed time for stats_begin_step
      if ( .not. l_restart ) then
        call stats_begin_timestep( time_current-time_initial+dtmain, dtmain ) ! Intent(in)
      else
        ! Different elapsed time for restart
        ! Joshua Fasching March 2008
        call stats_begin_timestep( time_current-time_restart+dtmain, dtmain ) ! Intent(in)
      end if

      ! If we're doing an inputfields run, get the values for our
      ! model arrays from a netCDF or GrADS file
      if ( l_input_fields ) then
        call compute_timestep( iunit, stat_file_zt, .false., time_current, &    ! Intent(in)
                               itime_nearest )                                  ! Intent(out)

        call stat_fields_reader( max( itime_nearest, 1 ) )                     ! Intent(in)
      end if

      if ( invalid_model_arrays( ) ) then
        err_code = clubb_var_equals_NaN ! Check for NaN values in the model arrays
        exit ! Leave the main loop
      end if

      call advance_clubb_forcings( dtmain, &  ! Intent(in)
                                   err_code ) ! Intent(inout)

      if ( err_code == clubb_rtm_level_not_found ) exit

      if ( l_stats_samp ) then
        ! Total microphysical tendency of vapor and cloud water mixing ratios
        call stat_update_var( irvm_mc, rvm_mc, zt ) ! kg/kg/s
        call stat_update_var( ircm_mc, rcm_mc, zt ) ! kg/kg/s
        call stat_update_var( irtm_mc, rvm_mc+rcm_mc, zt ) ! kg/kg/s
        call stat_update_var( ithlm_mc, thlm_mc, zt ) ! K/s
      end if 

      ! Add microphysical tendencies to the forcings
      rtm_forcing(:) = rtm_forcing(:) + rcm_mc(:) + rvm_mc(:)

      ! This is a kluge added because the _tndcy subroutines will sometimes
      ! compute radht from an analytic formula and add it to thlm_forcing before
      ! we reach this point in the clubb_driver. -dschanen 17 Aug 2009
      if ( l_bugsrad ) then
        thlm_forcing(:) = thlm_forcing(:) + thlm_mc(:) + radht(:)
      else
        thlm_forcing(:) = thlm_forcing(:) + thlm_mc(:)
      end if

      ! Compute number of iterations for closure loop
      if ( time_current > time_spinup ) then
        niterlong = 1
        dt        = dtmain
      else
        niterlong = floor( dtmain / dtclosure )
        dt        = dtclosure
      end if

!-------------------------------------------------------------------------------
!                                Closure loop
!-------------------------------------------------------------------------------

      do i1=1, niterlong
        ! Call the parameterization one timestep
        call advance_clubb_core & 
             ( .false., dt, fcor, sfc_elevation, &                  ! Intent(in)
               thlm_forcing, rtm_forcing, um_forcing, vm_forcing, & ! Intent(in)
               sclrm_forcing, edsclrm_forcing, wm_zm, wm_zt, &      ! Intent(in)
               wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &         ! Intent(in)
               wpsclrp_sfc, wpedsclrp_sfc,  &                       ! Intent(in)
               p_in_Pa, rho_zm, rho, exner, &                       ! Intent(in)
               um, vm, upwp, vpwp, up2, vp2, &                      ! Intent(inout)
               thlm, rtm, wprtp, wpthlp, wpthvp, &                  ! Intent(inout)
               wprcp, Kh_zt, wp2, wp3, &                            ! Intent(inout)
               rtp2, thlp2, rtpthlp, &                              ! Intent(inout)
               sigma_sqd_w, tau_zm, rcm, cloud_frac, &              ! Intent(inout)
               rcm_in_layer, cloud_cover, &                         ! Intent(inout)
               sclrm, sclrp2, sclrprtp, sclrpthlp, &                ! Intent(inout)
               wpsclrp, edsclrm, pdf_params, &                      ! Intent(inout)
               err_code )                                           ! Intent(inout)

        wp2_zt = max( zm2zt( wp2 ), wtol_sqd ) ! Positive definite quantity

        ! Advance a microphysics scheme
        call advance_clubb_microphys &
             ( i, dt, rho, rho_zm, p_in_Pa, exner, cloud_frac, thlm, & ! Intent(in)
               rtm, rcm, wm_zt, wm_zm,                               & ! Intent(in)
               Kh_zm, wp2_zt, pdf_params,                            & ! Intent(in)
               Ncnm, hydromet,                                       & ! Intent(inout)
               rvm_mc, rcm_mc, thlm_mc, err_code )                     ! Intent(out)

         ! Advance a radiation scheme
         ! With this call ordering, snow and ice water mixing ratio will be
         ! updated by the microphysics, but thlm and rtm will not.  This
         ! somewhat inconsistent, but we would need to move the call to
         ! radiation before the call the microphysics to change this.
         ! -dschanen 17 Aug 2009
         call advance_clubb_radiation &
              ( rho_zm, p_in_Pa, exner, cloud_frac, thlm, & ! Intent(in)
                rtm, rcm, hydromet,                       & ! Intent(in)
                radht, Frad, Frad_SW_up, Frad_LW_up,      & ! Intent(out)
                Frad_SW_down, Frad_LW_down )                ! Intent(out)

        call stats_end_timestep( )

        ! Set Time
        ! Advance time here, not in advance_clubb_core,
        ! in order to facilitate use of stats.
        ! A host model, e.g. WRF, would advance time outside
        ! of advance_clubb_core.  Vince Larson 7 Feb 2006
        if ( i1 < niterlong ) then
          time_current = time_initial + (i-1) * dtmain  & 
                       + i1 * dtclosure
        else if ( i1 == niterlong ) then
          time_current = time_initial + i * dtmain
        end if

        ! This was moved from above to be less confusing to the user,
        ! since before it would appear as though the last timestep
        ! was not executed. -dschanen 19 May 08
        if ( l_stats_last .and. l_stdout ) then
          write(unit=fstdout,fmt='(a,i8,a,f10.1)') 'iteration = ',  & 
            i, '; time = ', time_current
        end if

        if ( fatal_error( err_code ) ) exit

      end do ! i1=1..niterlong

!-------------------------------------------------------------------------------
!                             End Closure Loop
!-------------------------------------------------------------------------------

      if ( fatal_error( err_code ) ) exit

    end do ! i=1, ifinal

!-------------------------------------------------------------------------------
!                       End Main Time Stepping Loop
!-------------------------------------------------------------------------------

    ! Free memory
    if( thlm_sponge_damp_settings%l_sponge_damping ) then
      call finalize_tau_sponge_damp( thlm_sponge_damp_profile )
    end if

    if( rtm_sponge_damp_settings%l_sponge_damping ) then
      call finalize_tau_sponge_damp( rtm_sponge_damp_profile )
    end if

    if( uv_sponge_damp_settings%l_sponge_damping ) then
      call finalize_tau_sponge_damp( uv_sponge_damp_profile )
    end if

    if( l_t_dependent ) then
      call finalize_t_dependent_input()
    end if

    call finalize_extend_atm( )

    call cleanup_clubb_core( .false. )

    call stats_finalize( )

#ifdef UNRELEASED_CODE
    if ( l_latin_hypercube_sampling ) then
      call latin_hypercube_2D_close( )
    end if
#endif

    return
  end subroutine run_clubb

  !-----------------------------------------------------------------------
  subroutine initialize_clubb &
             ( iunit, forcings_file_path, psfc, zm_init, &
               thlm, rtm, um, vm, ug, vg, wp2, up2, vp2, rcm, &
               wm_zt, wm_zm, em, exner, &
               tau_zt, tau_zm, thvm, p_in_Pa, & 
               rho, rho_zm, Lscale, rtm_ref, thlm_ref, & 
               Kh_zt, Kh_zm, um_ref, vm_ref, & 
               hydromet, Ncnm, &
               sclrm, edsclrm )
    ! Description:
    !   Execute the necessary steps for the initialization of the
    !   CLUBB model run.
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use constants, only:  & 
      Cp,  &  ! Variable(s)
      Lv,  &
      ep2,  &
      ep1,  &
      emin,  &
      grav, &
      zero_threshold, &
      cm3_per_m3, kappa, p0, Rd

    use parameters_tunable, only:  & 
      taumax,  &  ! Variable(s)
      c_K

    use parameters_model, only:  & 
      T0,  &  ! Variable(s)
      sclr_dim, &
      edsclr_dim, &
      hydromet_dim

    use parameters_microphys, only: &
      Ncm_initial,     & ! Variable(s)
      micro_scheme

    use parameters_radiation, only: radiation_top ! Variable(s)

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: zm2zt, zt2zm ! Procedure(s)

    use sounding, only: n_snd_var, read_sounding ! Procedure(s)

    use model_flags, only: &
        l_uv_nudge, & ! Variable(s)
        l_tke_aniso

#ifdef UNRELEASED_CODE

    use lba, only: lba_init ! Procedure(s)

#endif

    use time_dependent_input, only: &
      initialize_t_dependent_input, & ! Procedure(s)
      l_t_dependent ! Variable(s)

    use extend_atmosphere_mod, only: &
      l_use_default_std_atmosphere, & ! Procedure(s)
      load_extend_std_atm, &
      convert_snd2extend_atm, &
      determine_extend_atmos_bounds

    use mpace_a, only: mpace_a_init ! Procedure(s)

    use mixing_length, only: compute_length ! Procedure(s)

    use error_code, only: clubb_no_error ! Variable(s)

    use constants, only: fstderr ! Variables(s)

    use array_index, only: &
      iisclr_thl, iiedsclr_thl, & ! Variable(s)
      iiNcm

    ! Joshua Fasching
    ! March 2008
    use saturation, only: sat_mixrat_liq, sat_rcm ! Procedure(s)

    use hydrostatic_mod, only: hydrostatic ! Procedure(s)

    use soil_vegetation, only: sfc_soil_T_in_K, deep_soil_T_in_K, veg_T_in_K ! Variable(s)

    use sponge_layer_damping, only: &
      thlm_sponge_damp_settings, &
      rtm_sponge_damp_settings, &
      uv_sponge_damp_settings, &
      thlm_sponge_damp_profile, &
      rtm_sponge_damp_profile, &
      uv_sponge_damp_profile, &
      initialize_tau_sponge_damp ! Procedure(s0

    use input_names, only: &
      z_name, &
      temperature_name, &
      thetal_name, &
      wm_name, &
      omega_name

    use input_reader, only: &
      one_dim_read_var

    implicit none

    intrinsic :: min, max, trim, sqrt, size

    ! Input
    integer, intent(in) :: iunit

    character(len=*), intent(in) :: &
      forcings_file_path ! Path to the .dat files containing the forcings

    real, intent(in) :: psfc, zm_init ! Pressure at the surface [Pa]

    ! Output
    real, dimension(gr%nnzp), intent(inout) ::  & 
      thlm,            & ! Theta l mean                  [K] 
      rtm,             & ! Total water mixing ratio      [kg/kg]
      um,              & ! u wind                        [m/s]
      vm,              & ! v wind                        [m/s]
      ug,              & ! u geostrophic wind            [m/s] 
      vg,              & ! u geostrophic wind            [m/s] 
      wp2,             & ! w'^2                          [m^2/s^2]
      up2,             & ! u'^2                          [m^2/s^2]
      vp2,             & ! v'^2                          [m^2/s^2]
      rcm,             & ! Cloud water mixing ratio      [kg/kg]
      wm_zt, wm_zm,    & ! w wind                        [m/s]
      em,              & ! Turbulence kinetic energy     [m^2/s^2]
      exner,           & ! Exner function                [-] 
      tau_zm, tau_zt,  & ! Dissipation time              [s]
      thvm,            & ! Virtual potential temperature [K]
      p_in_Pa,         & ! Pressure                      [Pa]
      rho, rho_zm,     & ! Density                       [kg/m^3]
      Lscale,          & ! Mixing length                 [m] 
      Kh_zt, Kh_zm,    & ! Eddy diffusivity              [m^2/s]
      um_ref,          & ! Initial profile of u wind     [m/s]
      vm_ref,          & ! Initial profile of v wind     [m/s]
      rtm_ref,         & ! Initial profile of rtm        [kg/kg]
      thlm_ref           ! Initial profile of thlm       [K]

    real, dimension(gr%nnzp,hydromet_dim), intent(inout) :: &
      hydromet ! Hydrometeor species    [kg/kg] or [#/kg]

    real, dimension(gr%nnzp), intent(inout) :: &
      Ncnm ! Cloud nuclei number concentration (COAMPS microphysics)

    ! Output
    real, dimension(gr%nnzp,sclr_dim), intent(out) ::  & 
      sclrm      ! Standard passive scalar [units vary]

    real, dimension(gr%nnzp,edsclr_dim), intent(out) ::  & 
      edsclrm    ! Eddy diffusivity passive scalar [units vary]

    ! Local Variables
    real, dimension(gr%nnzp) :: tmp1

    real :: cloud_top_height ! [m]
    real :: emax

    integer :: k, err_code

    type(one_dim_read_var), dimension(n_snd_var) :: sounding_retVars

    type(one_dim_read_var), dimension(sclr_dim) :: sclr_sounding_retVars

    character(len=50) :: &
      theta_type, & ! Type of temperature sounding 
      alt_type,   & ! Type of altitude sounding
      subs_type     ! Type of large-scale subsidence sounding
    !-----------------------------------------------------------------------

    err_code = clubb_no_error

    ! Read sounding information
    call read_sounding( iunit, runtype, psfc, zm_init, &          ! Intent(in) 
                        thlm, theta_type, rtm, um, vm, ug, vg,  & ! Intent(out)
                        alt_type, p_in_Pa, subs_type, wm_zt, &    ! Intent(out)
                        sclrm, edsclrm, sounding_retVars, sclr_sounding_retVars ) ! Intent(out)

    ! Prepare extended sounding for radiation
    if( l_use_default_std_atmosphere ) then

      call load_extend_std_atm( iunit ) ! Intent (in)

    else

      call convert_snd2extend_atm( n_snd_var, psfc, zm_init, sclr_dim,    & ! Intent(in)
                                sounding_retVars, sclr_sounding_retVars )   ! Intent(in)
    end if

    call determine_extend_atmos_bounds( gr%nnzp, gr%zt,            & ! Intent(in)
                                        gr%zm, gr%dzm,             & ! Intent(in)
                                        radiation_top,             & ! Intent(in)
                                        extend_atmos_bottom_level, & ! Intent(out)
                                        extend_atmos_top_level,    & ! Intent(out)
                                        extend_atmos_range_size,   & ! Intent(out)
                                        lin_int_buffer )             ! Intent(out)

    ! Covert sounding input to CLUBB compatible input
    select case( trim( alt_type ) )
    case ( z_name )

      if (theta_type == temperature_name ) then
        write(fstderr,*) 'Interpetation of sounding files with z as the independent ', &
        'variable and absolute temperature as the temperature variable has not ', &
        'been implemented. Either specify pressure as the independent variable. or ', &
        'thm/thlm as the temperature variable'
        stop
      end if

      ! At this point, thlm actually contains theta (except for DYCOMS).
      ! We need to compute liquid water content, and initilialize thlm properly

      ! First, compute approximate pressure using theta
      call hydrostatic( thlm, psfc, &                         ! Intent(in)
                        p_in_Pa, exner, rho, rho_zm )         ! Intent(out)

      ! Second, use this pressure to compute liquid water
      ! from excess saturation

      do k = 1,gr%nnzp
        rcm(k) = &
           max( rtm(k) - sat_mixrat_liq( p_in_Pa(k), thlm(k) * exner(k) ), &
                zero_threshold )
      enddo

      ! Compute initial theta-l

      select case ( trim( theta_type ) )
        !select case ( trim( runtype ) )
      case ( thetal_name )
        !case ( "dycoms2_rf01", "astex_a209", "nov11_altocu", &
        !      "clex9_nov02", "clex9_oct14", "dycoms2_rf02" )
        ! thlm profile that is initially saturated at points.
        ! thlm profile remains the same as in the input sounding.
        ! use iterative method to find initial rcm.
        do k =1, gr%nnzp, 1
          rcm(k) = sat_rcm( thlm(k), rtm(k), p_in_Pa(k), exner(k) )
        end do

      case default ! ('theta[K]')
        ! Initial profile is non-saturated thlm or any type of theta.
        thlm = thlm - Lv/(Cp*exner) * rcm

        ! Testing of passive scalars
        if ( iisclr_thl > 0 ) then
          sclrm(:,iisclr_thl) = thlm
        end if
        if ( iiedsclr_thl > 0 ) then
          edsclrm(:,iiedsclr_thl) = thlm
        end if

      end select

      ! Now, compute initial thetav

      thvm = thlm + ep1 * T0 * rtm  & 
                  + ( Lv/(Cp*exner) - ep2 * T0 ) * rcm

      ! Recompute more accurate initial exner function and pressure using thvm

      call hydrostatic( thvm, psfc, &                    ! Intent(in)
                        p_in_Pa, exner, rho, rho_zm )    ! Intent(out)

    case default ! ('Press[Pa]')

      ! Set the pressure at the lowest thermodynamic level (k=1), which is below
      ! the model lower boundary, to psfc, which is the pressure at the model
      ! lower boundary (or surface), which is located at momentum level 1.
      ! This is consistent with what is done in subroutine hydrostatic, which is
      ! called when the sounding is given in terms of altitude rather than
      ! pressure.  This is also a good way for the code to keep track of the
      ! surface pressure.
      p_in_Pa(1) = psfc

      ! Set the value of exner.
      exner(1) = ( psfc/p0 )**kappa
      do k=2, gr%nnzp
        exner(k) = (p_in_Pa(k)/p0) ** kappa  ! zt
      end do

      if ( trim( theta_type ) == temperature_name ) then
        thlm = thlm / exner
      end if

      select case ( trim( theta_type ) )
        !select case ( trim( runtype ) )
      case ( thetal_name )
        !case ( "dycoms2_rf01", "astex_a209", "nov11_altocu", &
        !      "clex9_nov02", "clex9_oct14", "dycoms2_rf02" )
        ! thlm profile that is initially saturated at points.
        ! thlm profile remains the same as in the input sounding.
        ! use iterative method to find initial rcm.
        do k =1, gr%nnzp, 1
          rcm(k) = sat_rcm( thlm(k), rtm(k), p_in_Pa(k), exner(k) )
        end do

      case default ! ('theta[K]')
        ! Initial profile is non-saturated thlm or any type of theta.
        thlm = thlm - Lv/(Cp*exner) * rcm
        ! Testing of passive scalars
        if ( iisclr_thl > 0 ) then
          sclrm(:,iisclr_thl) = thlm
        end if
        if ( iiedsclr_thl > 0 ) then
          edsclrm(:,iiedsclr_thl) = thlm
        end if

      end select

      ! Now, compute initial thetav

      thvm = thlm + ep1 * T0 * rtm  & 
                  + ( Lv/(Cp*exner) - ep2 * T0 ) * rcm

      ! Recompute more accurate initial exner function and pressure using thvm

      do k=1,gr%nnzp
        rho(k) = p_in_Pa(k) / ( Rd * thvm(k) * exner(k) )
      end do

      ! Interpolate density back to momentum grid

      rho_zm = max( zt2zm( rho ), zero_threshold )   ! Positive definitequantity
      rho_zm(1) = p_in_Pa(1) / ( Rd * thvm(1) * exner(1) )

    end select ! either 'z[m]' or 'Press[Pa]'

    ! Determine initial value cloud droplet number concentration for the
    ! Morrison microphysics
    select case ( trim( micro_scheme ) )
    case ( "morrison" )
      ! Lower boundary condition
      hydromet(1,iiNcm) = 0.

      hydromet(2:gr%nnzp-1,iiNcm) = cm3_per_m3 * Ncm_initial / rho(2:gr%nnzp-1)

      ! Upper boundary condition
      hydromet(gr%nnzp,iiNcm) = 0.

    case ( "coamps" )
      ! Initialize Ncnm as in COAMPS
      Ncnm(1:gr%nnzp) = 30.0 * (1.0 + exp( -gr%zt(1:gr%nnzp)/2000.0 )) * cm3_per_m3
    end select


    ! Initialize imposed w
    select case ( trim( subs_type ) ) ! Perform different operations based off
    !                                   the sounding file
    case ( wm_name )

      wm_zm = zt2zm( wm_zt )
      wm_zm(1) = 0.0
      wm_zm(gr%nnzp) = 0.0

    case ( omega_name )

      do k=2,gr%nnzp
        wm_zt(k) = -wm_zt(k) / ( grav*rho(k) )
      end do

      wm_zt(1) = 0.0
      wm_zt(gr%nnzp) = 0.0

      wm_zm = zt2zm( wm_zt )
      wm_zm(gr%nnzp) = 0.0

    case default ! This should not happen

      wm_zt = 0.0
      wm_zm = 0.0

    end select

    ! Initialize damping
    if( thlm_sponge_damp_settings%l_sponge_damping ) then
      call initialize_tau_sponge_damp( dt, thlm_sponge_damp_settings, & ! Intent(in)
                                       thlm_sponge_damp_profile )       ! Intent(out)
    end if

    if( rtm_sponge_damp_settings%l_sponge_damping ) then
      call initialize_tau_sponge_damp( dt, rtm_sponge_damp_settings, & ! Intent(in)
                                       rtm_sponge_damp_profile )       ! Intent(out)
    end if

    if(uv_sponge_damp_settings%l_sponge_damping) then
      call initialize_tau_sponge_damp( dt, uv_sponge_damp_settings, &  ! Intent(in)
                                       uv_sponge_damp_profile )        ! Intent(out)
    end if



    ! Initilize Time Dependant Input

    if( l_t_dependent ) then
      call initialize_t_dependent_input &
                   ( iunit, runtype, gr%nnzp, gr%zt, p_in_Pa )
    end if

    ! Initialize TKE and other fields as needed

    select case ( trim( runtype ) )

      ! Generic case
    case ( "generic" )

      em = 1.0

      ! GCSS BOMEX
    case ( "bomex" )

!---> Reduction of initial sounding for stability
!         do k = 1, gr%nnzp
!            em(k) = 1.0 - (gr%zm(k)/3000.0)
!            if ( em(k) < emin ) then
!               em(k) = emin
!            end if
!         end do
!         em(1) = em(2)
!         em(gr%nnzp) = em(gr%nnzp-1)
!<--- End reduction of initial sounding for stability 24 Jan 07

      em(:) = emin

      ! GCSS ARM
    case ( "arm" )

!---> Reduction of initial sounding for stability
!         do k = 1, gr%nnzp
!            if ( gr%zm(k) < 150.0 ) then
!               em(k) = ( 0.15 * (1.0 - gr%zm(k)/150.0) ) / rho_zm(k)
!            else
!               em(k) = emin
!            end if
!         end do
!         em(1) = em(2)
!         em(gr%nnzp) = em(gr%nnzp-1)
!<--- End reduction of initial sounding for stability 24 Jan 07

      em(:) = emin

#ifdef UNRELEASED_CODE
      ! March 2000 ARM case
    case ( "arm_0003" )

      em = 1.0

      ! 3 year ARM case
    case ( "arm_3year" )

      em = 1.0

      ! June 27 1997 ARM case
    case ( "arm_97" )

      em = 1.0
      ! twp_ice
    case ( "twp_ice" )

      em = 1.0
      ! twp_ice case

    case ( "cloud_feedback_s6", "cloud_feedback_s6_p2k",   &
           "cloud_feedback_s11", "cloud_feedback_s11_p2k", &
           "cloud_feedback_s12", "cloud_feedback_s12_p2k" )

      em = 1.0

#endif

      ! GCSS FIRE Sc
    case ( "fire" )

      cloud_top_height = 700. ! 700 m is the top of the cloud in FIRE
      do k=1,gr%nnzp
        if ( gr%zm(k) < cloud_top_height ) then
          em(k) = 1.
        else
          em(k) = emin
        end if
      end do
      em(1) = em(2)
      em(gr%nnzp) = em(gr%nnzp-1)

      ! GCSS ATEX
    case ( "atex" )

      um = max( um, -8. )

!---> Reduction of initial sounding for stability
!         do k = 1, gr%nnzp
!           em(k) = 1.0 - (gr%zm(k)/3000.0)
!           if ( em(k) < emin ) then
!             em(k) = emin
!           end if
!         end do
!         em(1) = em(2)
!         em(gr%nnzp) = em(gr%nnzp-1)
!<--- End reduction of initial sounding for stability 24 Jan 07

      em(:) = emin

      ! GCSS DYCOMS II RF01
    case ( "dycoms2_rf01" )
      cloud_top_height = 800. ! 800 m is the top of the cloud in RF01
      do k=1,gr%nnzp
        if ( gr%zm(k) < cloud_top_height ) then
          em(k) = 0.5
        else
          em(k) = emin
        end if
      end do
      em(1) = em(2)
      em(gr%nnzp) = em(gr%nnzp-1)

      ! GCSS DYCOMS II RF02
    case ( "dycoms2_rf02" )

      em = 1.0

      ! Brian for Nov. 11 altocumulus case.
    case ( "nov11_altocu" )

      ! Vince Larson reduced initial forcing.  4 Nov 2005
!          em = 1.0
!          em = 0.1
      ! 4150 + 2800 m is the top of the cloud in Nov11
      cloud_top_height = 2800. + gr%zm(1)
      do k=1,gr%nnzp
        if ( gr%zm(k) < cloud_top_height ) then

          ! Modification by Adam Smith, 08 April 2008
          ! Reducing the value of em appears to reduce error in the
          ! updated Nov.11 case
          em(k) = 0.01
          ! End of ajsmith4's modification

        else
          em(k) = emin
        end if
      end do
      em(1) = em(2)
      em(gr%nnzp) = em(gr%nnzp-1)
      ! End Vince Larson's change.

      ! Adam Smith addition for June 25 altocumulus case.
    case ( "jun25_altocu" )

      ! Vince Larson reduced initial forcing.  4 Nov 2005
!          em = 1.0
!          em = 0.1
!          do k=1,gr%nnzp
!            if ( gr%zm(k) < 1400. ) then
!               em(k) = 0.1
!            else
!               em(k) = emin
!            end if
!          end do

      ! Note: emin = 1.0e-6, defined in constants.F
      ! Adam Smith, 28 June 2006
      ! Note: now emin = 1.5 * wtol_sqd
      ! Brian Griffin;  Nov. 26, 2008.
      do k = 1, gr%nnzp
        em(k) = 0.01
      end do

      em(1) = em(2)
      em(gr%nnzp) = em(gr%nnzp-1)
      ! End Vince Larson's change.

      ! End of ajsmith4's addition

      ! Adam Smith addition for CLEX-9: Nov. 02 altocumulus case.
    case ( "clex9_nov02" )

      ! Vince Larson reduced initial forcing.  4 Nov 2005
!          em = 1.0
!          em = 0.1
      ! 4150 + 1400 m is the top of the cloud in Nov11
      cloud_top_height = 2200. + gr%zm(1)
      do k=1,gr%nnzp
        if ( gr%zm(k) < cloud_top_height ) then
          em(k) = 0.01
        else
          em(k) = emin
        end if
      end do
      em(1) = em(2)
      em(gr%nnzp) = em(gr%nnzp-1)
      ! End Vince Larson's change.

      ! End of ajsmith4's addition

      ! Adam Smith addition for CLEX-9: Oct. 14 altocumulus case.
    case ( "clex9_oct14" )

      ! Vince Larson reduced initial forcing.  4 Nov 2005
!          em = 1.0
!          em = 0.1
      ! 4150 + 1400 m is the top of the cloud in Nov11
      cloud_top_height = 3500. + gr%zm(1)
      do k=1,gr%nnzp
        if ( gr%zm(k) < cloud_top_height ) then
          em(k) = 0.01
        else
          em(k) = emin
        end if
      end do
      em(1) = em(2)
      em(gr%nnzp) = em(gr%nnzp-1)
      ! End Vince Larson's change.

      ! End of ajsmith4's addition

#ifdef UNRELEASED_CODE
    case ( "lba" )

      em = 0.1
      call lba_init( iunit, forcings_file_path )
#endif

      ! Michael Falk for mpace_a Arctic Stratus case.
    case ( "mpace_a" )

      cloud_top_height = 2000.
      emax = 1.0

      do k=1,gr%nnzp
        if ( gr%zm(k) < cloud_top_height ) then
          em(k) = emax
        else
          em(k) = emin
        end if
      end do
      em(1) = em(2)
      em(gr%nnzp) = em(gr%nnzp-1)

      call mpace_a_init( iunit, forcings_file_path )

      ! Michael Falk for mpace_b Arctic Stratus case.
    case ( "mpace_b" )

      cloud_top_height = 1300. ! 1300 m is the cloud top in mpace_b.  Michael Falk 17 Aug 2006
      emax = 1.0

      do k=1,gr%nnzp

        if ( gr%zm(k) < cloud_top_height ) then
          em(k) = emax
        else
          em(k) = emin
        end if
      enddo
      em(1) = em(2)
      em(gr%nnzp) = em(gr%nnzp-1)

      ! Brian Griffin for COBRA CO2 case.
    case ( "cobra" )

      em = 0.1

      ! Michael Falk for RICO tropical cumulus case, 13 Dec 2006
    case ( "rico" )

      cloud_top_height = 1500.
      emax = 1.0
      do k=1,gr%nnzp
        if ( gr%zm(k) < cloud_top_height ) then
          em(k) = emax
        else
          em(k) = emin
        end if
      enddo

      em(1) = em(2)
      em(gr%nnzp) = em(gr%nnzp-1)

      ! Michael Falk for GABLS2 case, 29 Dec 2006
    case ( "gabls2" )

      cloud_top_height = 800.  ! per GABLS2 specifications
      emax = 0.5
      do k=1,gr%nnzp
        if ( gr%zm(k) < cloud_top_height ) then
          em(k) = emax * (1 - (gr%zm(k)/cloud_top_height))
        else
          em(k) = emin
        end if
      end do

      em(1) = em(2)
      em(gr%nnzp) = em(gr%nnzp-1)


    case ( "gabls3_night" )
      em = 1.0
    case ( "gabls3" )
      em = 1.0

      veg_T_in_K = 300.
      sfc_soil_T_in_K = 300.
      deep_soil_T_in_K = 288.58

    end select

    ! End Initialize TKE and other fields as needed

    !!!! Initialize w'^2 based on initial TKE !!!

    if ( l_tke_aniso ) then

      ! TKE:  em = (1/2) * ( w'^2 + u'^2 + v'^2 )
      ! Evenly divide TKE into its component
      ! contributions (w'^2, u'^2, and v'^2).

      wp2 = (2.0/3.0) * em
      up2 = (2.0/3.0) * em
      vp2 = (2.0/3.0) * em

    else

      ! TKE:  em = (3/2) * w'^2

      wp2 = (2.0/3.0) * em

    endif

    ! Compute mixing length
    call compute_length( thvm, thlm, rtm, rcm, & ! Intent(in)
                         em, p_in_Pa, exner,   & ! Intent(in)    
                         err_code,             & ! Intent(inout)
                         Lscale )                ! Intent(out)

    ! Dissipation time
    tmp1 = sqrt( max( emin, zm2zt( em ) ) )
    tau_zt = min( Lscale / tmp1, taumax )
    tau_zm = min( ( max( zt2zm( Lscale ), zero_threshold ) &
                   / sqrt( max( emin, em ) ) ), taumax )

    ! Modification to damp noise in stable region
! Brian commented this out to match code found in advance_clubb_core.
! Vince Larson commented out because it may prevent turbulence from
!    initiating in unstable regions.  7 Jul 2007
!    do k=1,gr%nnzp
!      if ( wp2(k) <= 0.005 ) then
!        tau_zt(k) = taumin
!        tau_zm(k) = taumin
!      end if
!    end do
! End Vince Larson's commenting.

    ! Eddy diffusivity coefficient
    ! c_K is 0.548 usually (Duynkerke and Driedonks 1987)

    Kh_zt = c_K * Lscale * tmp1
    Kh_zm = c_K * max( zt2zm( Lscale ), zero_threshold )  & 
                * sqrt( max( em, emin ) )

    ! Moved this to be more general -dschanen July 16 2007
    if ( l_uv_nudge .or. uv_sponge_damp_settings%l_sponge_damping ) then
      um_ref = um ! Michael Falk addition for nudging code.  27 Sep/1 Nov 2006
      vm_ref = vm ! ditto
    end if

    if ( thlm_sponge_damp_settings%l_sponge_damping ) then
      thlm_ref = thlm ! Added for nudging code
    end if

    if ( rtm_sponge_damp_settings%l_sponge_damping ) then
      rtm_ref  = rtm
    end if

    return
  end subroutine initialize_clubb
  !-----------------------------------------------------------------------
  subroutine restart_clubb &
             ( iunit, runfile, &
               restart_path_case, time_restart, & 
               upwp, vpwp, wm_zt, wm_zm,  & 
               wpthlp, wprtp, & 
               wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc )
    ! Description:
    !   Execute the necessary steps for the initialization of the
    !   CLUBB model to a designated point in the submitted GrADS file.
    !-----------------------------------------------------------------------
    use inputfields,only:  & 
        input_type, &  ! Variable(s)
        input_um, input_vm, input_rtm, input_thlm, & 
        input_wp2, input_wprtp, input_wpthlp,  & 
        input_wp3, input_rtp2, input_thlp2,  & 
        input_rtpthlp, input_upwp, input_vpwp, & 
        input_ug, input_vg, input_rcm,  & 
        input_wm_zt, input_exner, input_em, & 
        input_p, input_rho, input_rho_zm, & 
        input_Lscale, input_Lscale_up, input_Lscale_down, & 
        input_Kh_zt, input_Kh_zm, input_tau_zm, input_tau_zt, & 
        input_wpthvp, &
        input_thl1, input_thl2, input_a, input_s1, input_s2, &
        input_stdev_s1, input_stdev_s2, input_rc1, input_rc2, &
        input_thvm, input_rrainm,input_Nrm,  & 
        input_rsnowm, input_ricem, input_rgraupelm,  & 
        input_thlm_forcing, input_rtm_forcing, & 
        input_up2, input_vp2, input_sigma_sqd_w, input_Ncm,  & 
        input_Ncnm, input_Nim, input_cloud_frac, input_sigma_sqd_w_zt, &
        input_veg_T_in_K, input_deep_soil_T_in_K, &
        input_sfc_soil_T_in_K, stat_file_zt

    use inputfields, only: &
      compute_timestep,  & ! Procedure(s)
      stat_fields_reader, &
      set_filenames

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: zt2zm ! Procedure(s)

    use constants, only: fstderr ! Variables(s)

    use stats_precision, only: time_precision ! Variable(s)

    use model_flags, only: &
      l_soil_veg ! Variable(s)

    use parameters_microphys, only : &
      micro_scheme, & ! Variable
      l_cloud_sed

    implicit none

    ! Input Variables

    integer, intent(in) :: iunit

    character(len=*), intent(in) ::  & 
      runfile,            & ! Filename for the namelist
      restart_path_case     ! Path to GrADS data for restart

    real(kind=time_precision), intent(in) :: & 
      time_restart

    ! Input/Output Variables
    real, dimension(gr%nnzp), intent(inout) ::  & 
      upwp,            & ! u'w'                         [m^2/s^2]
      vpwp,            & ! v'w'                         [m^2/s^2]
      wm_zt, wm_zm,    & ! w wind                       [m/s]
      wpthlp,          & ! w' th_l'                     [(m K)/s]
      wprtp              ! w' r_t'                      [(kg m)(kg s)]

    ! Output
    real, intent(out) :: & 
      wpthlp_sfc,      & ! w'theta_l' surface flux   [(m K)/s]
      wprtp_sfc,       & ! w'rt' surface flux        [(m kg)/(kg s)]
      upwp_sfc,        & ! u'w' at surface           [m^2/s^2] 
      vpwp_sfc           ! v'w' at surface           [m^2/s^2]

    ! Local variables
    integer :: timestep


    ! --- Begin Code ---

    ! Inform inputfields module
    input_type = "clubb"
    input_um   = .true.
    input_vm   = .true.
    input_rtm  = .true.
    input_thlm = .true.
    input_wp2  = .true.
    input_ug   = .true.
    input_vg   = .true.
    input_rcm  = .true.
    input_wm_zt  = .true.
    input_exner = .true.
    input_em = .true.
    input_p = .true.
    input_rho = .true.
    input_rho_zm = .true.
    input_Lscale = .true.
    input_Lscale_up = .true.
    input_Lscale_down = .true.
    input_Kh_zt = .true.
    input_Kh_zm = .true.
    input_tau_zm = .true.
    input_tau_zt = .true.
    input_thvm = .true.
    input_wpthvp = .true.
    input_thl1 = .true.
    input_thl2 = .true.
    input_a    = .true.
    input_s1   = .true.
    input_s2   = .true.
    input_stdev_s1  = .true.
    input_stdev_s2  = .true.
    input_rc1  = .true.
    input_rc2  = .true.

    select case ( trim( micro_scheme ) )
    case ( "coamps" )
      input_rrainm = .true.
      input_rsnowm = .true.
      input_ricem = .true.
      input_rgraupelm = .true.
      input_Ncnm = .true.
      input_Ncm = .true.
      input_Nrm = .true.
      input_Nim =  .true.

    case ( "morrison" )
      input_rrainm = .true.
      input_rsnowm = .true.
      input_ricem = .true.
      input_rgraupelm = .true.
      input_Ncnm = .false.
      input_Ncm = .true.
      input_Nrm = .true.
      input_Nim =  .true.

    case ( "khairoutdinov_kogan" )
      input_rrainm = .true.
      input_rsnowm = .false.
      input_ricem = .false.
      input_rgraupelm = .false.
      input_Ncnm = .false.
      input_Ncm = .true.
      input_Nrm = .true.
      input_Nim = .false.

    case default
      input_rrainm = .false.
      input_rsnowm = .false.
      input_ricem = .false.
      input_rgraupelm = .false.
      input_Ncnm = .false.
      input_Ncm = .false.
      input_Nrm = .false.
      input_Nim = .false.

    end select

    if ( l_cloud_sed ) then
      input_Ncm = .true.
    end if

    if ( l_soil_veg ) then
      input_veg_T_in_K = .true.
      input_deep_soil_T_in_K = .true.
      input_sfc_soil_T_in_K = .true.
    end if

    input_wprtp = .true.
    input_wpthlp = .true.
    input_wp3 = .true.
    input_rtp2 = .true.
    input_thlp2 = .true.
    input_rtpthlp = .true.
    input_upwp = .true.
    input_vpwp = .true.
    input_thlm_forcing = .true.
    input_rtm_forcing = .true.
    input_up2 = .true.
    input_vp2 = .true.
    input_sigma_sqd_w = .true.
    input_cloud_frac  = .true.
    input_sigma_sqd_w_zt = .true.
    call set_filenames( "../"//trim( restart_path_case ) )
    ! Determine the nearest timestep in the GRADS file to the
    ! restart time.
    call compute_timestep &
      ( iunit, stat_file_zt, .true., time_restart, &! Intent(in)
        timestep )                                                                 ! Intent(out)

    ! Sanity check for input time_restart
    if ( timestep < 0 ) then
      write(fstderr,*) "Invalid time_restart in "// & 
        "file: "//trim( runfile )
      stop
    end if

    ! Read data from stats files
    call stat_fields_reader( timestep )  ! Intent(in)

    wm_zm = zt2zm( wm_zt )

    wpthlp_sfc = wpthlp(1)
    wprtp_sfc  = wprtp(1)
    upwp_sfc   = upwp(1)
    vpwp_sfc   = vpwp(1)

    return
  end subroutine restart_clubb

  !----------------------------------------------------------------------
  subroutine advance_clubb_forcings( dt, err_code )

    ! Description:
    !   Calculate tendency and surface variables
    ! References:
    !   None
    !----------------------------------------------------------------------

    ! Modules to be included
    use model_flags, only:  &
      l_soil_veg

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: zt2zm, zm2zt ! Procedure(s)

    use variables_diagnostic_module, only: &
      radht, um_ref,  & ! Variable(s)
      vm_ref, Frad,  Frad_SW_up,  &
      Frad_SW_down, Frad_LW_down, thvm, ustar, & 
#ifdef UNRELEASED_CODE
    ug, vg, &
#endif
    soil_heat_flux

    use variables_diagnostic_module, only: wpedsclrp ! Passive scalar variables

    use variables_prognostic_module, only: &
      rtm_forcing, thlm_forcing,  & ! Variable(s)
      wm_zt, wm_zm, rho, rtm, thlm, p_in_Pa, & 
      exner, rcm, rho_zm, um, psfc, vm, & 
      upwp_sfc, vpwp_sfc, Tsfc, & 
      wpthlp_sfc, wprtp_sfc, &
#ifdef UNRELEASED_CODE
      um_forcing, vm_forcing, &
#endif
      SE, LE

    use stats_variables, only: &
      ish, & ! Variable(s)
      ilh, &
      iwpthlp_sfc, &
      iwprtp_sfc, &
      iupwp_sfc, &
      ivpwp_sfc, &
      iustar, &
      isoil_heat_flux, &
      l_stats_samp, &
      sfc

    use stats_type, only: stat_update_var_pt ! Procedure(s)

    use constants, only: & 
      Cp, Lv, kappa, p0, & ! Variable(s)
      rc_tol, fstderr, cm3_per_m3

    use variables_prognostic_module, only:  & 
      sclrm_forcing,   & ! Passive scalar variables
      edsclrm_forcing, & ! 
      wpsclrp,  & 
      wpsclrp_sfc,  &
      wpedsclrp_sfc

    use stats_precision, only: time_precision ! Variable(s)

    use time_dependent_input, only: &
      apply_time_dependent_forcings, &
      l_t_dependent
       
    use soil_vegetation, only: advance_soil_veg, veg_T_in_K

    use array_index, only: & 
      iisclr_rt, iisclr_thl, &  ! Variable(s)
      iiedsclr_rt, iiedsclr_thl

    ! Case specific modules
    use arm, only: arm_tndcy, arm_sfclyr ! Procedure(s)

#ifdef UNRELEASED_CODE
    use arm_0003, only: arm_0003_sfclyr ! Procedure(s)

    use arm_3year, only: arm_3year_sfclyr ! Procedure(s)

    use arm_97, only: arm_97_sfclyr ! Procedure(s)

    use astex, only: astex_tndcy, astex_sfclyr ! Procedure(s)
#endif

    use atex, only: atex_tndcy, atex_sfclyr ! Procedure(s)

    use bomex, only: bomex_tndcy, bomex_sfclyr ! Procedure(s)

#ifdef UNRELEASED_CODE
    use clex9_nov02, only: clex9_nov02_tndcy ! Procedure(s)

    use clex9_oct14, only: clex9_oct14_tndcy ! Procedure(s)

    use cobra, only: cobra_tndcy, cobra_sfclyr ! Procedure(s)
#endif

    use dycoms2_rf01, only:  & 
        dycoms2_rf01_tndcy, dycoms2_rf01_sfclyr ! Procedure(s)

    use dycoms2_rf02, only:  & 
        dycoms2_rf02_tndcy, dycoms2_rf02_sfclyr ! Procedure(s)

    use fire, only: fire_tndcy, sfc_momentum_fluxes, & 
                    sfc_thermo_fluxes ! Procedure(s)

    use gabls2, only: gabls2_tndcy, gabls2_sfclyr ! Procedure(s)

#ifdef UNRELEASED_CODE
    use gabls3, only: gabls3_sfclyr ! Procedures(s)

    use gabls3_night, only: gabls3_night_sfclyr

    use rico, only: rico_tndcy, rico_sfclyr ! Procedure(s)

    use lba, only: lba_tndcy, lba_sfclyr ! Procedure(s)

    use cloud_feedback, only: cloud_feedback_sfclyr ! Procedure(s)
#endif

    use mpace_a, only: mpace_a_tndcy, mpace_a_sfclyr ! Procedure(s)

    use mpace_b, only: mpace_b_tndcy, mpace_b_sfclyr ! Procedure(s)

#ifdef UNRELEASED_CODE
    use nov11, only: nov11_altocu_tndcy ! Procedure(s)

    use twp_ice, only: twp_ice_sfclyr ! Procedure(s)

    use jun25, only: jun25_altocu_tndcy ! Procedure(s)
#endif

    use wangara, only: wangara_tndcy, wangara_sfclyr ! Procedure(s)

    use cos_solar_zen_mod, only: cos_solar_zen ! Function

    use parameters_radiation, only: rad_scheme  ! Variable(s)

    implicit none

    ! Input Variables
    real(kind=time_precision), intent(in) :: & 
      dt         ! Model timestep                            [s]

    ! Input/Output Variables
    integer, intent(inout) :: & 
      err_code

    ! Local Variables

    real ::  &
      um_sfc, &  ! um interpolated to momentum level 1 (sfc) [m/s]
      vm_sfc     ! vm interpolated to momentum level 1 (sfc) [m/s]

    real :: &
      wpthep, & ! w'theta_e'                            [m K/s]
      amu0      ! Cosine of the solar zenith angle      [-]

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!                    FIND ALL DIAGNOSTIC VARIABLES
!-----------------------------------------------------------------------

    !----------------------------------------------------------------
    ! Find the cosine of the solar zenith angle if needed
    ! The abs() clipping prevents an error with sunray_sw code.
    !----------------------------------------------------------------
    if ( trim( rad_scheme ) == "simplified" ) then
      amu0 = max( real( cos_solar_zen( day, month, year, time_current, rlat, rlon ) ), 0. )
    else
      amu0 = -999.0
    end if

    !----------------------------------------------------------------
    ! Set vertical velocity, w, and compute large-scale forcings
    !----------------------------------------------------------------


    select case ( runtype )

    case ( "arm" ) ! ARM Cu case
      call arm_tndcy( time_current, &                     ! Intent(in)   
                      thlm_forcing, radht, rtm_forcing, & ! Intent(out)
                      sclrm_forcing, edsclrm_forcing )    ! Intent(out)

#ifdef UNRELEASED_CODE

    case ( "astex_a209" ) ! ASTEX Sc case for K & K
      call astex_tndcy( wm_zt, wm_zm,  &                ! Intent(out) 
                        thlm_forcing, rtm_forcing , &   ! Intent(out)
                        sclrm_forcing, edsclrm_forcing )! Intent(out)
#endif

    case ( "atex" ) ! ATEX case
      call atex_tndcy( time_current, time_initial, &   ! Intent(in)
                       rtm, rho, rcm, exner, &         ! Intent(in)
                       err_code, &                     ! Intent(inout)
                       wm_zt, wm_zm, Frad, radht, &    ! Intent(out)
                       thlm_forcing, rtm_forcing, &    ! Intent(out)
                       sclrm_forcing, edsclrm_forcing )! Intent(out)

    case ( "bomex" ) ! BOMEX Cu case
      call bomex_tndcy( rtm, radht, &                    ! Intent(out)
                        thlm_forcing, rtm_forcing, &     ! Intent(out)
                        sclrm_forcing, edsclrm_forcing ) ! Intent(out)

#ifdef UNRELEASED_CODE
    case ( "clex9_nov02" ) ! CLEX-9: Nov. 02 Altocumulus case.
      call clex9_nov02_tndcy( time_current, time_initial, rlat, rlon, &  ! Intent(in)
                              rcm, exner, rho, &                         ! Intent(in)
                              wm_zt, wm_zm, thlm_forcing, rtm_forcing, & ! Intent(out)
                              Frad, radht, &                             ! Intent(out)
                              sclrm_forcing, edsclrm_forcing )           ! Intent(out)

    case ( "clex9_oct14" ) ! CLEX-9: Oct. 14 Altocumulus case.
      call clex9_oct14_tndcy( time_current, time_initial, rlat, rlon, &    ! Intent(in) 
                              rcm, exner, rho, &                           ! Intent(in)
                              wm_zt, wm_zm, thlm_forcing, rtm_forcing, &   ! Intent(out)
                              Frad, radht, &                               ! Intent(out)
                              sclrm_forcing, edsclrm_forcing )             ! Intent(out)
    case ( "cobra" )
      call cobra_tndcy( thlm_forcing, rtm_forcing, &     ! Intent(out)
                        sclrm_forcing, edsclrm_forcing ) ! Intent(out)
#endif

    case ( "dycoms2_rf01" ) ! DYCOMS2 RF01 case
      call dycoms2_rf01_tndcy( rho, rho_zm, rtm, rcm, exner, &  ! Intent(in)
                               err_code, &                      ! Intent(inout)
                               Frad, radht,  &                  ! Intent(out)
                               thlm_forcing, rtm_forcing,  &    ! Intent(out)
                               sclrm_forcing, edsclrm_forcing ) ! Intent(out)

    case ( "dycoms2_rf02" ) ! DYCOMS2 RF02 case
      call dycoms2_rf02_tndcy( rho,                     &          ! Intent(in)
                               rho_zm, rtm, rcm, exner,  &         ! Intent(in)
                               err_code, wm_zt, wm_zm,   &         ! Intent(inout)
                               thlm_forcing, rtm_forcing, &        ! Intent(out) 
                               Frad, radht, &                      ! Intent(out)
                               sclrm_forcing, edsclrm_forcing )    ! Intent(out)

    case ( "fire" ) ! FIRE Sc case
      call fire_tndcy( rho, rcm, exner,  &                ! Intent(in)
                       Frad, radht, &                     ! Intent(out)
                       thlm_forcing, rtm_forcing, &       ! Intent(out) 
                       sclrm_forcing, edsclrm_forcing )   ! Intent(out)


    case ( "gabls2" ) ! GABLS 2 case
      call gabls2_tndcy( time_current, time_initial,  &   ! Intent(in) 
                         wm_zt, wm_zm, thlm_forcing, &    ! Intent(out)
                         rtm_forcing, radht, &            ! Intent(out)
                         sclrm_forcing, edsclrm_forcing ) ! Intent(out)

#ifdef UNRELEASED_CODE
    case ( "jun25_altocu" ) ! June 25 Altocumulus case.
      call jun25_altocu_tndcy( time_current, time_initial, rlat, rlon,  &  ! Intent(in) 
                               rcm, exner, rho, &                          ! Intent(in)
                               wm_zt, wm_zm, thlm_forcing, rtm_forcing, &  ! Intent(inout)
                               Frad, radht, &                              ! Intent(inout)
                               sclrm_forcing, edsclrm_forcing )            ! Intent(out)

    case ( "lba" )
      call lba_tndcy( time_current, &                  ! Intent(in) 
                      radht,  &                        ! Intent(out)
                      thlm_forcing, rtm_forcing, &     ! Intent(out)
                      sclrm_forcing, edsclrm_forcing ) ! Intent(out)
#endif

    case ( "mpace_a" ) ! mpace_a arctic stratus case

      call mpace_a_tndcy( time_current, amu0, &        ! Intent(in) 
                          rho, p_in_Pa, rcm, &                       ! Intent(in)
                          wm_zt, wm_zm, thlm_forcing, rtm_forcing, & ! Intent(out)
                          Frad, radht, um_ref, vm_ref, &             ! Intent(out)
                          sclrm_forcing, edsclrm_forcing )           ! Intent(out)

    case ( "mpace_b" ) ! mpace_b arctic stratus case

      call mpace_b_tndcy( amu0, &                      ! Intent(in)
                          rho,  p_in_Pa, thvm, rcm, &                ! Intent(in)
                          wm_zt, wm_zm, thlm_forcing, rtm_forcing, & ! Intent(out)
                          Frad, radht,  &                            ! Intent(out)
                          sclrm_forcing, edsclrm_forcing )           ! Intent(out)

#ifdef UNRELEASED_CODE
    case ( "nov11_altocu" ) ! Nov. 11 Altocumulus case.
      call nov11_altocu_tndcy( time_current, time_initial, dt, &           ! Intent(in)
                               day, month, year, rlat, rlon, &             ! Intent(in)
                               rcm, exner, rho, &                          ! Intent(in)
                               rtm, &                                      ! Intent(inout)
                               wm_zt, wm_zm, thlm_forcing, rtm_forcing, &  ! Intent(out)
                               Frad, radht, &                              ! Intent(out)
                               sclrm_forcing, edsclrm_forcing )            ! Intent(out)

    case ( "rico" ) ! RICO case
      call rico_tndcy( exner, &                            ! Intent(in)
                       thlm_forcing, rtm_forcing, radht, & ! Intent(out)   
                       sclrm_forcing, edsclrm_forcing )    ! Intent(out)
#endif

    case ( "wangara" ) ! Wangara dry CBL
      call wangara_tndcy( wm_zt, wm_zm,  &                  ! Intent(out) 
                          thlm_forcing, rtm_forcing, &      ! Intent(out)
                          sclrm_forcing, edsclrm_forcing )  ! Intent(out)

#ifdef UNRELEASED_CODE
    case ( "cloud_feedback_s6", "cloud_feedback_s6_p2k",   &
           "cloud_feedback_s11", "cloud_feedback_s11_p2k", &
           "cloud_feedback_s12", "cloud_feedback_s12_p2k", &
           "gabls3_night", "arm_97", "gabls3", "twp_ice",  &
           "arm_0003", "arm_3year" )
      if ( l_t_dependent ) then 
         call apply_time_dependent_forcings( time_current, gr%nnzp, rtm, rho, exner,&
          thlm_forcing, rtm_forcing, um_ref, vm_ref, um_forcing, vm_forcing, wm_zt, wm_zm, ug, vg, &
          sclrm_forcing, edsclrm_forcing )
      end if
#endif

    case default

      write(unit=fstderr,fmt=*)  & 
         "advance_clubb_forcings: Don't know how to handle " & 
         //"LS forcing for runtype: "//trim( runtype )
      stop

    end select



    !----------------------------------------------------------------
    ! Compute Surface Fluxes
    !----------------------------------------------------------------

    ! Boundary conditions for the second order moments

    ! Find the value of um at the surface (momentum level 1) by interpolating the
    ! values of um found at thermodynamic levels 2 and 1.  This will be helpful in
    ! computing the surface flux u'w'|_sfc.
    um_sfc = zt2zm( um, 1 )
    ! Find the value of vm at the surface (momentum level 1) by interpolating the
    ! values of vm found at thermodynamic levels 2 and 1.  This will be helpful in
    ! computing the surface flux v'w'|_sfc.
    vm_sfc = zt2zm( vm, 1 )

    select case ( trim( runtype ) )

    case ( "arm" )
      call arm_sfclyr( time_current, gr%zt(2), 1.1,  &              ! Intent(in)
                       thlm(2), um(2), vm(2), &                     ! Intent(in)
                       upwp_sfc, vpwp_sfc,  &                       ! Intent(out)
                       wpthlp_sfc, wprtp_sfc, ustar, &              ! Intent(out)
                       wpsclrp_sfc, wpedsclrp_sfc )                 ! Intent(in)

#ifdef UNRELEASED_CODE
    case ( "arm_0003" )
      call arm_0003_sfclyr( time_current, gr%zt(2), rho_zm(1), &   ! Intent(in)
                            thlm(2), um(2), vm(2), &               ! Intent(in)
                            upwp_sfc, vpwp_sfc,  &                 ! Intent(out)
                            wpthlp_sfc, wprtp_sfc, ustar, &        ! Intent(out)
                            wpsclrp_sfc, wpedsclrp_sfc )           ! Intent(out)

    case ( "arm_3year" )
      call arm_3year_sfclyr( time_current, gr%zt(2), rho_zm(1), &  ! Intent(in)
                             thlm(2), um(2), vm(2), &              ! Intent(in)
                             upwp_sfc, vpwp_sfc,  &                ! Intent(out)
                             wpthlp_sfc, wprtp_sfc, ustar, &       ! Intent(out)
                             wpsclrp_sfc, wpedsclrp_sfc )          ! Intent(out)


    case ( "arm_97" )
      call arm_97_sfclyr( time_current, gr%zt(2), rho_zm(1), &     ! Intent(in)
                          thlm(2), um(2), vm(2), &                 ! Intent(in)
                          upwp_sfc, vpwp_sfc,  &                   ! Intent(out)
                          wpthlp_sfc, wprtp_sfc, ustar, &          ! Intent(out)
                          wpsclrp_sfc, wpedsclrp_sfc )             ! Intent(out)

    case ( "astex_a209" )
      call astex_sfclyr( rho_zm(1), &                               ! Intent(in) 
                         upwp_sfc, vpwp_sfc, wpthlp_sfc,  &         ! Intent(out)
                         wprtp_sfc, wpsclrp_sfc, wpedsclrp_sfc )    ! Intent(out)

#endif

    case ( "atex" )
      call atex_sfclyr( um(2), vm(2), thlm(2), rtm(2), &  ! Intent(in)
                        exner(1), Tsfc, &                 ! Intent(in)
                        upwp_sfc, vpwp_sfc, &             ! Intent(out)
                        wpthlp_sfc, wprtp_sfc, ustar, &   ! Intent(out)
                        wpsclrp_sfc, wpedsclrp_sfc )      ! Intent(out)

    case ( "bomex" )
      call bomex_sfclyr( um(2), vm(2), rtm(2),  &                   ! Intent(in) 
                         upwp_sfc, vpwp_sfc, &                      ! Intent(out)
                         wpthlp_sfc, wprtp_sfc, ustar, &            ! Intent(out)
                         wpsclrp_sfc, wpedsclrp_sfc )               ! Intent(out)

#ifdef UNRELEASED_CODE
    case ( "cobra" )
      call cobra_sfclyr( time_current, gr%zt(2), rho_zm(1), &       ! Intent(in)
                         thlm(2), um(2), vm(2), &                   ! Intent(in)
                         upwp_sfc, vpwp_sfc, &                      ! Intent(out)
                         wpthlp_sfc, wprtp_sfc, ustar, &            ! Intent(out)
                         wpsclrp_sfc, wpedsclrp_sfc )               ! Intent(out)

    case ( "clex9_nov02" )
      ! There are no surface momentum or heat fluxes
      ! for the CLEX-9: Nov. 02 Altocumulus case.

      ! Ensure ustar is set
      ustar = 0

    case ( "clex9_oct14" )
      ! There are no surface momentum or heat fluxes
      ! for the CLEX-9: Oct. 14 Altocumulus case.

      ! Ensure ustar is set.
      ustar = 0

    case ( "cloud_feedback_s6", "cloud_feedback_s6_p2k",   &
           "cloud_feedback_s11", "cloud_feedback_s11_p2k", &
           "cloud_feedback_s12", "cloud_feedback_s12_p2k" ) ! Cloud Feedback cases

           call cloud_feedback_sfclyr( runtype, sfctype, & ! Intent(in)
                                       thlm(2), rtm(2), um(2), vm(2), &       ! Intent(in)
                                       exner(1), psfc, Tsfc, &                ! Intent(in)
                                       upwp_sfc, vpwp_sfc, &                  ! Intent(out)
                                       wpthlp_sfc, wprtp_sfc, ustar, &        ! Intent(out)
                                       wpsclrp_sfc, wpedsclrp_sfc )           ! Intent(out)

           ! If the surface type is 0, use fixed fluxes
           if ( sfctype == 0 ) then
             wpthlp_sfc = SE
             wprtp_sfc  = LE
             if ( iisclr_thl > 0 ) wpsclrp(:,iisclr_thl) = SE
             if ( iisclr_rt > 0 ) wpsclrp(:,iisclr_rt)   = LE
             if ( iiedsclr_thl > 0 ) wpedsclrp(:,iiedsclr_thl) = SE
             if ( iiedsclr_rt > 0 ) wpedsclrp(:,iiedsclr_rt)   = LE
           end if
#endif

    case ( "dycoms2_rf01" )
      call dycoms2_rf01_sfclyr( sfctype, Tsfc, psfc,  &             ! Intent(in)
                                exner(1), um(2), vm(2),  &          ! Intent(in)
                                thlm(2), rtm(2), rho_zm(1), &       ! Intent(in) 
                                upwp_sfc, vpwp_sfc,  &              ! Intent(out)
                                wpthlp_sfc, wprtp_sfc, ustar, &     ! Intent(out)
                                wpsclrp_sfc, wpedsclrp_sfc )        ! Intent(out)

    case ( "dycoms2_rf02" )
      call dycoms2_rf02_sfclyr( um(2), vm(2), &                     ! Intent(in)
                                upwp_sfc, vpwp_sfc, &               ! Intent(out)
                                wpthlp_sfc, wprtp_sfc, ustar, &     ! Intent(out)
                                wpsclrp_sfc, wpedsclrp_sfc )        ! Intent(out)

    case( "fire", "generic"  )  ! Generic setup and GCSS FIRE
      call sfc_momentum_fluxes( um_sfc, vm_sfc, &                   ! Intent(in)
                                upwp_sfc, vpwp_sfc, ustar )         ! Intent(out)
      ! sfctype = 0  fixed sfc sensible and latent heat fluxes
      !                   as given in the namelist
      ! sfctype = 1  bulk formula: uses given surface temperature
      !                   and assumes over ocean
      if ( sfctype == 0 ) then

        wpthlp_sfc = SE
        wprtp_sfc  = LE
        if ( iisclr_thl > 0 ) wpsclrp(:,iisclr_thl) = SE
        if ( iisclr_rt > 0 ) wpsclrp(:,iisclr_rt)   = LE
        if ( iiedsclr_thl > 0 ) wpedsclrp(:,iiedsclr_thl) = SE
        if ( iiedsclr_rt > 0 ) wpedsclrp(:,iiedsclr_rt)   = LE

      elseif ( sfctype == 1 ) then

        call sfc_thermo_fluxes( um(2), vm(2), &                     ! Intent(in)
                                Tsfc, psfc,  &                      ! Intent(in)
                                thlm(2), rtm(2), exner(1), &        ! Intent(in)
                                wpthlp_sfc, wprtp_sfc, &            ! Intent(out)
                                wpsclrp_sfc, wpedsclrp_sfc )        ! Intent(out)

      else

        write(unit=fstderr,fmt=*)  & 
          "Invalid value of sfctype = ", sfctype
        stop

      endif

    case ( "gabls2" )
      call gabls2_sfclyr( time_current, time_initial, &             ! Intent(in)
                          gr%zt(2), psfc, &                         ! Intent(in)
                          um(2), vm(2), thlm(2), rtm(2), exner(1), &! Intent(in)     
                          upwp_sfc, vpwp_sfc, &                     ! Intent(out)   
                          wpthlp_sfc, wprtp_sfc, ustar, &           ! Intent(out)
                          wpsclrp_sfc, wpedsclrp_sfc )              ! Intent(out)

#ifdef UNRELEASED_CODE
    case ( "gabls3" )
      call gabls3_sfclyr( um(2), vm(2), veg_T_in_K,             & ! Intent(in)
                          thlm(2), rtm(2), gr%zt(2), exner(1) , & ! Intent(in)
                          upwp_sfc, vpwp_sfc,                   & ! Intent(out)
                          wpthlp_sfc, wprtp_sfc, ustar )          ! Intent(out)
    case ( "gabls3_night" )
      call gabls3_night_sfclyr( time_current, um(2), vm(2),  & ! Intent(in)
                          thlm(2), rtm(2), gr%zt(2),         & ! Intent(in)
                          upwp_sfc, vpwp_sfc,                & ! Intent(out)
                          wpthlp_sfc, wprtp_sfc, ustar )       ! Intent(out)


    case ( "jun25_altocu" )
      ! There are no surface momentum or heat fluxes
      ! for the Jun. 25 Altocumulus case.

      ! Ensure ustar is set
      ustar = 0

    case ( "lba" )
      call lba_sfclyr( time_current, gr%zt(2), rho_zm(1), &        ! Intent(in)
                       thlm(2), um(2), vm(2), &                    ! Intent(in)
                       upwp_sfc, vpwp_sfc,  &                      ! Intent(out)
                       wpthlp_sfc, wprtp_sfc, ustar, &             ! Intent(out)
                       wpsclrp_sfc, wpedsclrp_sfc )                ! Intent(out)

#endif

    case ( "mpace_a" )
      call mpace_a_sfclyr( time_current, rho_zm(1), um(2), vm(2), & ! Intent(in)
                           upwp_sfc, vpwp_sfc, &                    ! Intent(out)
                           wpthlp_sfc, wprtp_sfc, ustar, &          ! Intent(out)
                           wpsclrp_sfc, wpedsclrp_sfc )             ! Intent(out)

    case ( "mpace_b" )
      call mpace_b_sfclyr( rho_zm(1), um(2), vm(2), &               ! Intent(in)
                           upwp_sfc, vpwp_sfc, &                    ! Intent(out)
                           wpthlp_sfc, wprtp_sfc, ustar, &          ! Intent(out)
                           wpsclrp_sfc, wpedsclrp_sfc )             ! Intent(out)

#ifdef UNRELEASED_CODE
    case ( "nov11_altocu" )
      ! There are no surface momentum or heat fluxes
      ! for the Nov. 11 Altocumulus case.

      ! Ensure ustar is set
      ustar = 0

    case ( "rico" )
      call rico_sfclyr( um(2), vm(2), thlm(2), rtm(2), &            ! Intent(in)
                        ! 299.8 K is the RICO SST; 101540 Pa is the sfc pressure.
                        !gr%zt(2), 299.8, 101540.,  &                ! Intent(in)
                        gr%zt(2), Tsfc, psfc, exner(1), &
                        upwp_sfc, vpwp_sfc, wpthlp_sfc, &           ! Intent(out) 
                        wprtp_sfc, ustar, &                         ! Intent(out)
                        wpsclrp_sfc, wpedsclrp_sfc )                ! Intent(out)

    case ( "twp_ice" )
      call twp_ice_sfclyr( gr%zt(2), Tsfc, exner(1), thlm(2), &     ! Intent(in)
                            um(2), vm(2), rtm(2), &                 ! Intent(in)
                            psfc, upwp_sfc, vpwp_sfc,  &            ! Intent(out)
                            wpthlp_sfc, wprtp_sfc, ustar, &         ! Intent(out)
                            wpsclrp_sfc, wpedsclrp_sfc )            ! Intent(out)

#endif

    case ( "wangara" )
      call wangara_sfclyr( time_current, um(2), vm(2), &            ! Intent(in)
                           upwp_sfc, vpwp_sfc, &                    ! Intent(out)
                           wpthlp_sfc, wprtp_sfc, ustar, &          ! Intent(out)
                           wpsclrp_sfc, wpedsclrp_sfc )             ! Intent(out)
    case default

      write(unit=fstderr,fmt=*)  & 
        "Invalid value of runtype = ", runtype
      stop

    end select ! runtype

    !---------------------------------------------------------------
    ! Compute Surface
    !---------------------------------------------------------------
    if ( l_soil_veg ) then
      wpthep = wpthlp_sfc + (Lv/Cp) * ((p0/psfc)**kappa) * wprtp_sfc

      call advance_soil_veg( real( dt ), rho_zm(1), &
                             Frad_SW_down(1) - Frad_SW_up(1), Frad_SW_down(1), &
                             Frad_LW_down(1), wpthep, soil_heat_flux )
    end if

    ! Store values of surface fluxes for statistics
    if ( l_stats_samp ) then

      call stat_update_var_pt( ish, 1, wpthlp_sfc*rho_zm(1)*Cp,& ! intent(in)
                               sfc )                             ! intent(inout)

      call stat_update_var_pt( ilh, 1, wprtp_sfc*rho_zm(1)*Lv, & ! intent(in)
                               sfc )                             ! intent(inout)

      call stat_update_var_pt( iwpthlp_sfc, 1, wpthlp_sfc, & ! intent(in)
                               sfc )                         ! intent(inout)

      call stat_update_var_pt( iwprtp_sfc, 1, wprtp_sfc, & ! intent(in)
                               sfc )                       ! intent(inout)

      call stat_update_var_pt( iupwp_sfc, 1, upwp_sfc, & ! intent(in)
                               sfc )                     ! intent(inout)

      call stat_update_var_pt( ivpwp_sfc, 1, vpwp_sfc, & ! intent(in)
                               sfc )                     ! intent(inout)

      call stat_update_var_pt( iustar, 1, ustar,  & ! intent(in)
                               sfc )                ! intent(inout)

      call stat_update_var_pt( isoil_heat_flux, 1, soil_heat_flux, & ! intent(in)
                               sfc )           ! intent(inout)
    endif


    return
  end subroutine advance_clubb_forcings

!-------------------------------------------------------------------------------
  subroutine advance_clubb_microphys &
             ( iter, dt, rho, rho_zm, p_in_Pa, exner, cloud_frac, thlm, &
               rtm, rcm, wm_zt, wm_zm, &
               Kh_zm, wp2_zt, pdf_params, Ncnm, hydromet, rvm_mc, rcm_mc, &
               thlm_mc, err_code )
! Description:
!   Advance a microphysics scheme
! References:
!   None
!-------------------------------------------------------------------------------

    use parameters_microphys, only: &
      micro_scheme, l_cloud_sed, Ncm_initial  ! Variables

    use constants, only: & 
      rc_tol, fstderr, cm3_per_m3 ! Variable(s)

    use stats_precision, only: time_precision ! Variable(s)

    use microphys_driver, only: advance_microphys ! Procedure(s)

    use cloud_sed_mod, only: cloud_drop_sed ! Procedure(s)

    use parameters_model, only: hydromet_dim ! Variable(s)

    use variables_prognostic_module, only: &
      pdf_parameter ! Type

    use error_code, only: lapack_error  ! Procedure(s)

    use grid_class, only: &
      gr ! Instance of a type

    use array_index, only: & 
      iiNcm

    implicit none

    ! External
    intrinsic :: trim

    ! Input Variables
    integer, intent(in) :: &
      iter ! Model iteration number

    real(kind=time_precision), intent(in) :: & 
      dt ! Model timestep                            [s]

    real, dimension(gr%nnzp), intent(in) :: &
      rho,        & ! Density on thermo. grid                           [kg/m^3] 
      rho_zm,     & ! Density on moment. grid                           [kg/m^3]
      p_in_Pa,    & ! Pressure.                                         [Pa] 
      exner,      & ! Exner function.                                   [-]
      cloud_frac, & ! Cloud fraction (thermodynamic levels)             [-]
      thlm,       & ! Liquid potential temperature                      [K]
      rtm,        & ! Total water mixing ratio, r_t (thermo. levels)    [kg/kg]
      rcm,        & ! Cloud water mixing ratio, r_c (thermo. levels)    [kg/kg]
      wm_zt,      & ! wm on thermo. grid.                               [m/s]
      wm_zm,      & ! wm on moment. grid.                               [m/s]
      Kh_zm,      & ! Eddy-diffusivity on momentum levels               [m^2/s]
      wp2_zt        ! w'^2 interpolated the thermo levels               [m^2/s^2]

    type(pdf_parameter), intent(in) :: & 
      pdf_params      ! PDF parameters   [units vary]

    ! Input/Output Variables
    real, dimension(gr%nnzp), intent(inout) :: &
      Ncnm ! Cloud nuclei number concentration (COAMPS microphyics)     [#/kg]

    real, dimension(gr%nnzp,hydromet_dim), intent(inout) :: &
      hydromet ! Hydrometeor species    [units vary]

    real, dimension(gr%nnzp), intent(inout) :: &
      thlm_mc,   & ! theta_l microphysical tendency [K/s]
      rcm_mc,    & ! r_c microphysical tendency     [(kg/kg)/s]
      rvm_mc       ! r_v microphysical tendency     [(kg/kg)/s]

    integer, intent(inout) :: & 
      err_code

    ! ---- Begin Code ----
    rcm_mc  = 0.0
    rvm_mc  = 0.0
    thlm_mc = 0.0

    !----------------------------------------------------------------
    ! Compute Microphysics
    !----------------------------------------------------------------
    ! Determine Ncm for K&K microphysics or cloud droplet sedimentation
    if ( l_cloud_sed .or. trim( micro_scheme ) == "khairoutdinov_kogan" ) then

      ! The following lines of code specify cloud droplet
      ! concentration (Ncm).  The cloud droplet concentration has
      ! been moved here instead of being stated in KK_microphys
      ! for the following reasons:
      !    a) The effects of cloud droplet sedimentation can be computed
      !       without having to call the precipitation scheme.
      !    b) Ncm tends to be a case-specific parameter.  Therefore, it
      !       is appropriate to declare in the same place as other
      !       case-specific parameters.
      !
      ! Someday, we could move the setting of Ncm to pdf_closure
      ! for the following reasons:
      !    a) The cloud water mixing ratio (rcm) is computed using the
      !       PDF scheme.  Perhaps someday Ncm can also be computed by
      !       the same scheme.
      !    b) It seems more appropriate to declare Ncm in the same place
      !       where rcm is computed.
      !
      ! Since cloud base (zb) is determined by the mixing ratio rc_tol,
      ! so will cloud droplet number concentration (Ncm).
      where ( rcm >= rc_tol )
        hydromet(:,iiNcm) = Ncm_initial * cm3_per_m3 / rho
      else where
        hydromet(:,iiNcm) = 0.0
      end where

    end if ! cloud_sed / K&K

    ! Call Khairoutdinov and Kogan (2000) scheme, or COAMPS for rain microphysics.

    if ( trim( micro_scheme ) /= "none" ) then

      call advance_microphys &
           ( iter, runtype, dt, time_current, &                         ! Intent(in)
             thlm, p_in_Pa, exner, rho, rho_zm, rtm, rcm, cloud_frac, & ! Intent(in)
             wm_zt, wm_zm, Kh_zm, pdf_params, &                         ! Intent(in)
             wp2_zt, &                                                  ! Intent(in)
             Ncnm, hydromet, &                                          ! Intent(inout)
             rvm_mc, rcm_mc, thlm_mc, &                                 ! Intent(inout)
             err_code )                                                 ! Intent(out)

      if ( lapack_error( err_code ) ) return

    end if

    if ( l_cloud_sed ) then

      call cloud_drop_sed( rcm, hydromet(:,iiNcm), rho_zm, rho, exner, & ! Intent(in)
                           rcm_mc, thlm_mc )              ! Intent(inout)

    end if

    return
  end subroutine advance_clubb_microphys

!-------------------------------------------------------------------------------
  subroutine advance_clubb_radiation &
             ( rho_zm, p_in_Pa, exner, cloud_frac, thlm, &
               rtm, rcm, hydromet, radht, Frad, Frad_SW_up, Frad_LW_up, &
               Frad_SW_down, Frad_LW_down )
! Description:
!   Compute a radiation tendency
! References:
!   None
!-------------------------------------------------------------------------------

    use constants, only: fstderr  ! Constant(s)

    use numerical_check, only: isnan2d, rad_check ! Procedure(s)

    use parameters_radiation, only: rad_scheme ! Variable(s)

    use error_code, only: & 
      clubb_at_least_debug_level ! Procedure(s)

    use parameters_model, only: hydromet_dim ! Variable(s)

    use array_index, only: iirsnowm, iiricem ! Variable(s)

    use grid_class, only: gr ! Instance of a type

    use grid_class, only: zt2zm ! Procedure

#ifdef radoffline
    use bugsrad_clubb_mod, only: bugsrad_clubb ! Procedure(s)
#endif
    implicit none

    ! External
    intrinsic :: trim

    ! Input Variables
    real, dimension(gr%nnzp), intent(in) :: &
      rho_zm,     & ! Density on moment. grid                          [kg/m^3]
      p_in_Pa,    & ! Pressure.                                        [Pa] 
      exner,      & ! Exner function.                                  [-]
      cloud_frac, & ! Cloud fraction (thermodynamic levels)            [-]
      thlm,       & ! Liquid potential temperature                     [K]
      rtm,        & ! Total water mixing ratio, r_t (thermo. levels)   [kg/kg]
      rcm           ! Cloud water mixing ratio, r_c (thermo. levels)   [kg/kg]

    real, dimension(gr%nnzp,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    ! Input/Output Variables
    real, dimension(gr%nnzp), intent(inout) :: &
      radht ! Radiative heating rate                   [K/s]

    ! Output Variables
    real, dimension(gr%nnzp), intent(out) :: &
      Frad,         & ! Total radiative flux                   [W/m^2]
      Frad_SW_up,   & ! Short-wave upwelling radiative flux    [W/m^2]
      Frad_LW_up,   & ! Long-wave upwelling radiative flux     [W/m^2]
      Frad_SW_down, & ! Short-wave upwelling radiative flux    [W/m^2]
      Frad_LW_down    ! Long-wave upwelling radiative flux     [W/m^2]

    ! Local Variables
    real, dimension(gr%nnzp) ::  & 
      rsnowm,  & ! Snow mixing ratio                         [kg/kg]
      ricem      ! Prisitine ice water mixing ratio          [kg/kg]


    !----------------------------------------------------------------
    ! BUGSrad Radiation
    !----------------------------------------------------------------

    if ( trim( rad_scheme ) == "bugsrad" ) then

#ifdef radoffline /*This directive is needed for BUGSrad to work with CLUBB.*/

      ! Copy snow and ice
      if ( iirsnowm > 0 ) then
        rsnowm = hydromet(1:gr%nnzp,iirsnowm)
      else
        rsnowm = 0.0
      endif

      if ( iiricem > 0 ) then
        ricem = hydromet(1:gr%nnzp,iiricem)
      else
        ricem = 0.0
      end if

      ! NaN checks added to detect possible errors with BUGSrad
      ! Joshua Fasching November 2007

      if ( clubb_at_least_debug_level( 2 ) ) then

        if ( isnan2d( thlm ) ) then
          write(fstderr,*) "thlm before BUGSrad is NaN"
        endif

        if ( isnan2d( rcm ) ) then
          write(fstderr,*) "rcm before BUGSrad is NaN"
        endif

        if ( isnan2d( rtm ) ) then
          write(fstderr,*) "rtm before BUGSrad is NaN"
        endif

        if ( isnan2d( rsnowm ) ) then
          write(fstderr,*) "rsnowm before BUGSrad is NaN"
        endif

        if ( isnan2d( ricem ) ) then
          write(fstderr,*) "ricem before BUGSrad is NaN"
        endif

        if ( isnan2d( cloud_frac ) ) then
          write(fstderr,*) "cloud_frac before BUGSrad is NaN"
        endif

        if ( isnan2d( p_in_Pa ) ) then
          write(fstderr,*) "p_in_Pa before BUGSrad is NaN"
        endif

        if ( isnan2d( exner ) ) then
          write(fstderr,*) "exner before BUGSrad is NaN"
        endif

        if ( isnan2d( rho_zm ) ) then
          write(fstderr,*) "rho_zm before BUGSrad is NaN"
        endif

        ! Check for impossible negative values
        call rad_check( thlm, rcm, rtm, ricem, &             ! Intent(in)
                        cloud_frac, p_in_Pa, exner, rho_zm ) ! Intent(in)

      endif  ! clubb_at_least_debug_level( 2 )

      call bugsrad_clubb( gr%zm, gr%nnzp, lin_int_buffer,        & ! Intent(in)
                          extend_atmos_range_size,               & ! Intent(in)
                          extend_atmos_bottom_level,             & ! Intent(in)
                          extend_atmos_top_level,                & ! Intent(in)
                          rlat, rlon,                            & ! Intent(in)
                          day, month, year, time_current,        & ! Intent(in)
                          thlm, rcm, rtm, rsnowm, ricem,         & ! Intent(in)
                          cloud_frac, p_in_Pa, zt2zm( p_in_Pa ), & ! Intent(in)
                          exner, rho_zm,                         & ! Intent(in)
                          radht, Frad,                           & ! Intent(out)
                          Frad_SW_up, Frad_LW_up,                & ! Intent(out)
                          Frad_SW_down, Frad_LW_down )             ! Intent(out)

      if ( clubb_at_least_debug_level( 2 ) ) then

        if ( isnan2d( Frad ) ) then
          write(fstderr,*) "Frad after BUGSrad is NaN"
          !write(fstderr,*) Frad
        endif

        if ( isnan2d( radht ) ) then
          write(fstderr,*) "radht after BUGSrad is NaN"
          !write(fstderr,*) radht
        endif

      endif  ! clubb_at_least_debug_level( 2 )

#else

      stop "Cannot call BUGSrad with these compile options."

#endif /*radoffline*/

    end if ! BUGSrad radiation scheme
  end subroutine advance_clubb_radiation

end module clubb_driver
