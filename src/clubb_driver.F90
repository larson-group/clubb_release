!-----------------------------------------------------------------------
! $Id$

module clubb_driver

! Description:
!   Contains the necessary subroutines to execute individual CLUBB
!   model runs, using one of the driver programs (the simplest case
!   being the clubb_standalone program).
!-----------------------------------------------------------------------

  use stats_precision, only: time_precision ! Variable(s)


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
    nzmax,     & ! Vertical extent in levels              [#]
    grid_type    ! 1 ==> evenly-spaced grid levels
  !                2 ==> stretched (unevenly-spaced) grid entered on
  !                      thermodynamic grid levels; momentum levels
  !                      halfway between thermodynamic levels (style
  !                      of SAM stretched grid).
  !                3 ==> stretched (unevenly-spaced) grid entered on
  !                      momentum grid levels; thermodynamic levels
  !                      halfway between momentum levels (style
  !                      of WRF stretched grid).

  real, private ::  & 
    deltaz,  & ! Change per grid level                 [m]
    zm_init    ! Initial point on the momentum grid    [m]

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
             ( params, runfile, err_code, l_stdout, l_input_fields )
    ! Description:
    !   Subprogram to integrate the pde equations for pdf closure.

    ! References:
    !   None
    !-----------------------------------------------------------------------

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: read_grid_heights ! Procedure(s)

    use parameter_indices, only: nparams ! Variable(s)

    use variables_diagnostic_module, only: ug, vg, em,  & ! Variable(s)
      tau_zt, thvm, Lscale, Kh_zm, & 
      um_ref, vm_ref, Ncnm, wp2_zt, &
      hydromet

    use variables_prognostic_module, only:  & 
      Tsfc, psfc, SE, LE, thlm, rtm,     & ! Variable(s)
      um, vm, wp2, rcm, wm_zt, wm_zm, exner, & 
      tau_zm, p_in_Pa, rho_zm, upwp, vpwp, wpthlp, wpthvp, & 
      Kh_zt, rho, wprtp, wpthlp_sfc, wprtp_sfc, & 
      upwp_sfc, vpwp_sfc, thlm_forcing, & 
      rtm_forcing, um_forcing, vm_forcing, &
      up2, vp2, wp3, rtp2, pdf_params, & 
      thlp2, rtpthlp, sigma_sqd_w, cf

    use variables_prognostic_module, only:  & 
      sclrm, sclrp2, sclrprtp, sclrpthlp, sclrm_forcing, & ! Variables
      wpsclrp, wpsclrp_sfc,  & 
      edsclrm, edsclrm_forcing, wpedsclrp_sfc


    use numerical_check, only: invalid_model_arrays ! Procedure(s)

    use inputfields, only: compute_timestep, grads_fields_reader ! Procedure(s)

    use inputfields, only: datafilet ! Variable(s)

    use clubb_core, only: & 
      setup_clubb_core,  & ! Procedure(s) 
      cleanup_clubb_core, & 
      advance_clubb_core

    use constants, only: fstdout, fstderr, & ! Variable(s)
      rttol, thltol, wtol

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

    use stats_variables, only: l_stats_last ! Variable(s

    use stats_subs, only:  & 
      stats_begin_timestep, stats_end_timestep,  & ! Procedure(s)
      stats_finalize, stats_init

    use sounding, only: sclr_max

    use time_dependant_input, only: l_time_dependant, &
      finalize_time_dependant_input

    implicit none

    ! Because Fortran I/O is not thread safe, we use this here to
    ! insure that no model uses the same file number simultaneously
    ! when doing a tuning run. -dschanen 31 Jan 2007
#ifdef _OPENMP
    integer :: omp_get_thread_num ! Function
#endif
    ! External
    intrinsic :: mod, real, int, trim, floor, max

    ! Input Variables
    logical, intent(in) ::  & 
      l_stdout,        & ! Whether to print output per timestep
      l_input_fields    ! Whether to set model variables from a file

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
      l_tke_aniso       ! For anisotropic turbulent kinetic energy,
    ! i.e. TKE = 1/2 (u'^2 + v'^2 + w'^2)

    character(len=6) :: &
      saturation_formula ! "bolton" approx. or "flatau" approx.

    character(len=50) ::  & 
      restart_path_case, & ! GrADS file used in case of restart
      forcings_file_path   ! Path to the forcing files

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

    integer :: i, i1 ! Internal Loop Variables
    integer :: iinit ! initial iteration

    integer ::  & 
      iunit,           & ! File unit used for I/O
      hydromet_dim,    & ! Number of hydrometeor species        [#]
      sclr_dim,        & ! Number of passive scalars            [#]
      edsclr_dim         ! Number of passive scalars            [#]

    integer :: itime_nearest ! Used for and inputfields run [s]

    ! Definition of namelists
    namelist /model_setting/  & 
      runtype, nzmax, grid_type, deltaz, zm_init, & 
      zt_grid_fname, zm_grid_fname,  & 
      day, month, year, rlat, rlon, & 
      time_initial, time_final, time_spinup, & 
      dtmain, dtclosure, & 
      sfctype, Tsfc, psfc, SE, LE, fcor, T0, ts_nudge, & 
      forcings_file_path, l_time_dependant, &
      l_soil_veg, l_tke_aniso, l_uv_nudge, l_restart, restart_path_case, & 
      time_restart, debug_level, & 
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
    zt_grid_fname = ''
    zm_grid_fname = ''

    day   = 1
    month = 1
    year  = 1900

    rlat = 0.
    rlon = 0.

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

    l_time_dependant = .false.

    l_soil_veg     = .false.
    l_tke_aniso    = .false.
    l_uv_nudge     = .false.
    l_restart      = .false.
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
      print *, "--------------------------------------------------"
      print *, "Model Settings"
      print *, "--------------------------------------------------"
      print *, "Preprocessing directives:"
#ifdef NETCDF
      print *, "-DNETCDF enabled"
#else
      print *, "-DNETCDF disabled"
#endif
#ifdef COAMPS_MICRO
      print *, "-DCOAMPS_MICRO enabled"
#else
      print *, "-DCOAMPS_MICRO disabled"
#endif
#ifdef TUNER
      print *, "-DTUNER enabled"
#else
      print *, "-DTUNER disabled"
#endif
#ifdef UNRELEASED_CODE
      print *, "-DUNRELEASED_CODE enabled"
#else
      print *, "-DUNRELEASED_CODE disabled"
#endif

      ! Pick some default values for model_setting
      print *, "&model_setting:"
      print *, "runtype = ", runtype
      print *, "nzmax = ", nzmax
      print *, "grid_type = ", grid_type
      print *, "deltaz = ", deltaz
      print *, "zm_init = ", zm_init
      print *, "zt_grid_fname = ", zt_grid_fname
      print *, "zm_grid_fname = ", zm_grid_fname

      print *, "day = ", day
      print *, "month = ", month
      print *, "year = ", year

      print *, "rlat = ", rlat
      print *, "rlon = ", rlon

      print *, "time_initial = ", time_initial
      print *, "time_final = ", time_final
      print *, "time_spinup = ", time_spinup

      print *, "dtmain = ", dtmain
      print *, "dtclosure = ", dtclosure

      print *, "sfctype = ", sfctype
      print *, "tsfc = ", tsfc
      print *, "psfc = ", psfc
      print *, "SE = ", SE
      print *, "LE = ", LE
      print *, "fcor = ", fcor
      print *, "T0 = ", T0
      print *, "ts_nudge = ", ts_nudge

      print *, "forcings_file_path = ", forcings_file_path

      print *, "l_time_dependant = ", l_time_dependant

      print *, "l_soil_veg = " , l_soil_veg
      print *, "l_tke_aniso = ", l_tke_aniso
      print *, "l_uv_nudge = ", l_uv_nudge
      print *, "l_restart = ", l_restart
      print *, "restart_path_case = ", restart_path_case
      print *, "time_restart = ", time_restart
      print *, "debug_level = ", debug_level

      print *, "sclr_dim = ", sclr_dim
      print *, "edsclr_dim = ", edsclr_dim
      print *, "iisclr_thl = ", iisclr_thl
      print *, "iisclr_rt = ", iisclr_rt
      print *, "iisclr_CO2 = ", iisclr_CO2

      print *, "sclr_tol = ", sclr_tol(1:sclr_dim)

      ! Pick some default values for stats_setting
      print *, "&stats_setting:"
      print *, "l_stats = ", l_stats
      print *, "fname_prefix = ", fname_prefix
      print *, "stats_fmt = ", stats_fmt
      print *, "stats_tsamp = ", stats_tsamp
      print *, "stats_tout = ", stats_tout

      print *, "Constant flags:"
      print *, "l_pos_def = ", l_pos_def
      print *, "l_hole_fill = ", l_hole_fill
      print *, "l_single_C2_Skw = ", l_single_C2_Skw
      print *, "l_gamma_Skw = ", l_gamma_Skw
      print *, "l_byteswap_io = ", l_byteswap_io

      print *, "Constant tolerances", "[units]"
      print *, "rttol = ", rttol, "[kg/kg]"
      print *, "thltol = ", thltol,"[K]"
      print *, "wtol = ", wtol, "[m/s]"

      print *, "--------------------------------------------------"

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

    call init_radiation( iunit, runfile ) ! Intent(in)

    ! Allocate & initialize variables,
    ! setup grid, setup constants, and setup flags

    call setup_clubb_core &                               
         ( nzmax, T0, ts_nudge, &                         ! Intent(in)
           hydromet_dim, sclr_dim, &                      ! Intent(in)
           sclr_tol(1:sclr_dim), edsclr_dim, params, &    ! Intent(in)
           l_soil_veg, &                                  ! Intent(in)
           l_uv_nudge, l_tke_aniso, saturation_formula, & ! Intent(in)
           .false., grid_type, deltaz, zm_init, &         ! Intent(in)
           momentum_heights, thermodynamic_heights, &     ! Intent(in)
           dummy_dx, dummy_dy, &                          ! Intent(in)
           err_code )                                     ! Intent(out)


    if ( err_code == clubb_var_out_of_bounds ) return


    ! Deallocate stretched grid altitude arrays
    deallocate( momentum_heights, thermodynamic_heights )

    if ( .not. l_restart ) then

      time_current = time_initial
      iinit = 1

      call initialize_clubb &
           ( iunit, trim( forcings_file_path ), psfc, zm_init, & ! Intent(in)
             thlm, rtm, um, vm, ug, vg, wp2, wp2_zt, up2, vp2, rcm,  & ! Intent(inout)
             wm_zt, wm_zm, em, exner, &           ! Intent(inout)
             tau_zt, tau_zm, thvm, p_in_Pa, &     ! Intent(inout)
             rho, rho_zm, Lscale, &               ! Intent(inout) 
             Kh_zt, Kh_zm, um_ref, vm_ref, &      ! Intent(inout)
             hydromet, Ncnm, &                    ! Intent(inout)
             sclrm, edsclrm )                     ! Intent(out)

    else  ! restart

      ! Joshua Fasching March 2008
      time_current = time_restart + dtmain
      ! time_current = time_restart

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
           ( iunit, runfile, trim( forcings_file_path ), & ! Intent(in)
             restart_path_case, time_restart, &            ! Intent(in)
             upwp, vpwp, wm_zt, wm_zm,  &                  ! Intent(inout)
             um_ref, vm_ref, wpthlp, wprtp, &        ! Intent(inout)
             wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc )             ! Intent(out)

    end if ! ~l_restart



#ifdef _OPENMP
    iunit = omp_get_thread_num( ) + 50
#else
    iunit = 50
#endif

    fdir = "../output/" ! Output directory

    ! Initialize statistics output
    call stats_init( iunit, fname_prefix, fdir, l_stats, stats_fmt, stats_tsamp, &! Intent(in)
                     stats_tout, runfile, gr%nnzp, gr%zt, gr%zm, &              ! Intent(in)
                     day, month, year, rlat, rlon, time_current, dtmain )       ! Intent(in)


    ! Time integration
    ! Call advance_clubb_core once per each statistics output time
    ifinal = floor( ( time_final - time_initial ) / dtmain )


!-------------------------------------------------------------------------------
!                         Main Time Stepping Loop
!-------------------------------------------------------------------------------

    do i = iinit, ifinal, 1


      ! When this time step is over, the time will be time + dtmain

      ! We use elapsed time for stats_begin_step
      if (.not. l_restart) then
        call stats_begin_timestep( time_current-time_initial+dtmain, dtmain )   ! Intent(in)
      else
        ! Different elapsed time for restart
        ! Joshua Fasching March 2008
        call stats_begin_timestep( time_current-time_restart, dtmain )          ! Intent(in)
      end if



      ! If we're doing an inputfields run, get the values for our
      ! model arrays from a GrADS file
      if ( l_input_fields ) then
        call compute_timestep( iunit, datafilet, .false., time_current, &       ! Intent(in)
                               itime_nearest )                                  ! Intent(out)

        call grads_fields_reader( max( itime_nearest, 1 ) )                     ! Intent(in)
      end if

      if ( invalid_model_arrays( ) ) then
        err_code = clubb_var_equals_NaN ! Check for NaN values in the model arrays
        exit ! Leave the main loop
      end if

      call advance_clubb_forcings( i, dtmain, & ! Intent(in)
                                  err_code )    ! Intent(out)

      if ( err_code == clubb_rtm_level_not_found ) exit

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
        call advance_clubb_core & 
             ( i, .false., dt, fcor, &                             ! Intent(in)
               thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &! Intent(in)
               sclrm_forcing, edsclrm_forcing, wm_zm, wm_zt, &     ! Intent(in)
               wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &        ! Intent(in)
               wpsclrp_sfc, wpedsclrp_sfc,  &                      ! Intent(in)
               p_in_Pa, rho_zm, rho, exner, &                      ! Intent(in)
               um, vm, upwp, vpwp, up2, vp2, &                     ! Intent(inout)
               thlm, rtm, wprtp, wpthlp, wpthvp, &                 ! Intent(inout)
               Kh_zt, wp2, wp3, &                                  ! Intent(inout)
               rtp2, thlp2, rtpthlp, &                             ! Intent(inout)
               sigma_sqd_w, tau_zm, rcm, cf, &                     ! Intent(inout)
               sclrm, sclrp2, sclrprtp, sclrpthlp, &               ! Intent(inout)
               wpsclrp, edsclrm, pdf_params, &                     ! Intent(inout)
               err_code )                                          ! Intent(inout)


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

    if( l_time_dependant ) then
      call finalize_time_dependant_input()
    end if

    call cleanup_clubb_core( .false. )

    call stats_finalize( )

    return
  end subroutine run_clubb

  !-----------------------------------------------------------------------
  subroutine initialize_clubb &
             ( iunit, forcings_file_path, psfc, zm_init, &
               thlm, rtm, um, vm, ug, vg, wp2, wp2_zt, up2, vp2, rcm, &
               wm_zt, wm_zm, em, exner, &
               tau_zt, tau_zm, thvm, p_in_Pa, & 
               rho, rho_zm, Lscale, & 
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

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: zm2zt, zt2zm ! Procedure(s)

    use sounding, only: read_sounding, & ! Procedure(s)
                        z_name, thetal_name, wm_name, omega_name, temp_name ! Variable(s)

    use model_flags, only: &
        l_uv_nudge, & ! Variable(s)
        l_tke_aniso

#ifdef UNRELEASED_CODE
    use arm_0003, only: arm_0003_init ! Procedure(s)

    use arm_3year, only: arm_3year_init ! Procedure(s)

    use time_dependant_input, only: initialize_time_dependant_input,&
                                    l_time_dependant

    use lba, only: lba_init ! Procedure(s)

    use twp_ice, only: twp_ice_init ! Procedure(s)
#endif

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
      wp2, wp2_zt,     & ! w'^2                          [m^2/s^2]
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
      vm_ref             ! Initial profile of v wind     [m/s]

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

    character(len=50) :: &
      theta_type, & ! Type of temperature sounding 
      alt_type,   & ! Type of altitude sounding
      subs_type     ! Type of large-scale subsidence sounding
    !-----------------------------------------------------------------------

    err_code = clubb_no_error

    ! Read sounding information

    call read_sounding( iunit, runtype, psfc, zm_init, &          ! Intent(in) 
                        thlm, theta_type, rtm, um, vm, ug, vg,  & ! Intent(out)
                        alt_type, p_in_Pa, subs_type, wm_zt, &
                        sclrm, edsclrm )                          ! Intent(out)

    select case( trim( alt_type ) )
    case ( z_name )

      if (theta_type == temp_name ) then
        write(fstderr,*) 'Interpetation of sounding files with z as the independant ', &
        'variable and absolute temperature as the temperature variable has not ', &
        'been implemented. Either specify pressure as the independant variable. or ', &
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

      exner(1) = ( psfc/p0 )**kappa
      do k=2, gr%nnzp
        exner(k) = (p_in_Pa(k)/p0) ** kappa  ! zt
      end do

      if ( trim( theta_type ) == temp_name ) then
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


    ! Initilize Time Dependant Input

    if( l_time_dependant ) then
      call initialize_time_dependant_input &
                   ( iunit, runtype, gr%nnzp, gr%zt )
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
      call arm_0003_init( iunit, forcings_file_path )

      ! 3 year ARM case
    case ( "arm_3year" )

      em = 1.0
      call arm_3year_init( iunit, forcings_file_path )

      ! June 27 1997 ARM case
    case ( "arm_97" )

      em = 1.0
      ! twp_ice
    case ( "twp_ice" )

      em = 1.0
      call twp_ice_init( iunit, forcings_file_path )
      ! twp_ice case
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

    wp2_zt = zm2zt( wp2 )

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
    if ( l_uv_nudge ) then
      um_ref = um ! Michael Falk addition for nudging code.  27 Sep/1 Nov 2006
      vm_ref = vm ! ditto
    end if

    return
  end subroutine initialize_clubb
  !-----------------------------------------------------------------------
  subroutine restart_clubb &
             ( iunit, runfile, forcings_file_path, &
               restart_path_case, time_restart, & 
               upwp, vpwp, wm_zt, wm_zm,  & 
               um_ref, vm_ref, wpthlp, wprtp, & 
               wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc )
    ! Description:
    !   Execute the necessary steps for the initialization of the
    !   CLUBB model to a designated point in the submitted GrADS file.
    !-----------------------------------------------------------------------
    use inputfields,only:  & 
        datafile, input_type, &  ! Variable(s)
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
        input_ss1, input_ss2, input_rc1, input_rc2, &
        input_thvm, input_rrainm,input_Nrm,  & 
        input_rsnowm, input_ricem, input_rgraupelm,  & 
        input_thlm_forcing, input_rtm_forcing, & 
        input_up2, input_vp2, input_sigma_sqd_w, input_Ncm,  & 
        input_Ncnm, input_Nim, input_cf, input_sigma_sqd_w_zt, &
        input_veg_T_in_K, input_deep_soil_T_in_K, &
        input_sfc_soil_T_in_K

    use inputfields, only: compute_timestep, grads_fields_reader ! Procedure(s)

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: zt2zm ! Procedure(s)

    use constants, only: fstderr ! Variables(s)

    use stats_precision, only: time_precision ! Variable(s)

    use time_dependant_input, only: l_time_dependant, & ! Variable(s)
                                    initialize_time_dependant_input ! Procedure(s)

    use model_flags, only: &
      l_uv_nudge, & ! Variable(s)
      l_soil_veg

    use parameters_microphys, only : &
      micro_scheme, & ! Variable
      l_cloud_sed

#ifdef UNRELEASED_CODE
    use arm_0003, only: arm_0003_init ! Procedure(s)

    use arm_3year, only: arm_3year_init ! Procedure(s)

    use lba, only: lba_init ! Procedure(s)

    use twp_ice, only: twp_ice_init ! Procedure(s)
#endif

    use mpace_a, only: mpace_a_init ! Procedure(s)

    use input_reader, only: read_one_dim_file, fill_blanks_one_dim_vars, & ! Procedures
      deallocate_one_dim_vars, & 
      one_dim_read_var ! Type

    use sounding, only: read_x_profile ! Procedure(s)

    implicit none
    integer, parameter :: nCol = 7

    ! Input Variables

    integer, intent(in) :: iunit

    character(len=*), intent(in) ::  & 
      runfile,            & ! Filename for the namelist
      forcings_file_path, & ! Path to the forcing files
      restart_path_case     ! Path to GrADS data for restart

    real(kind=time_precision), intent(in) :: & 
      time_restart

    ! Input/Output Variables
    real, dimension(gr%nnzp), intent(inout) ::  & 
      upwp,            & ! u'w'                         [m^2/s^2]
      vpwp,            & ! v'w'                         [m^2/s^2]
      wm_zt, wm_zm,    & ! w wind                       [m/s]
      um_ref,          & ! Initial profile of u wind    [m/s]
      vm_ref,          & ! Initial profile of v wind    [m/s]
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

    type(one_dim_read_var), dimension(nCol) :: retVars

    ! --- Begin Code ---

    ! Inform inputfields module
    datafile = "../"//trim( restart_path_case )
    input_type = "hoc"
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
    input_ss1  = .true.
    input_ss2  = .true.
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
    input_cf  = .true.
    input_sigma_sqd_w_zt = .true.

    ! Determine the nearest timestep in the GRADS file to the
    ! restart time.
    call compute_timestep &
      ( iunit, "../"//trim( restart_path_case )//"_zt.ctl", .true., time_restart, &! Intent(in)
        timestep )                                                                 ! Intent(out)

    ! Sanity check for input time_restart
    if ( timestep < 0 ) then
      write(fstderr,*) "Invalid time_restart in "// & 
        "file: "//trim( runfile )
      stop
    end if

    if ( l_uv_nudge ) then

      call read_one_dim_file( iunit, nCol, &
          '../input/case_setups/'//trim(runtype)//'_sounding.in', retVars )

      call fill_blanks_one_dim_vars( nCol, retVars )

      um_ref = read_x_profile(nCol, nzmax, 'u[m\s]', retVars)
      vm_ref = read_x_profile(nCol, nzmax,'v[m\s]', retVars)

      call deallocate_one_dim_vars( nCol, retVars )
    end if


    ! Read data from GrADS files
    call grads_fields_reader( timestep )                    ! Intent(in)

    ! Initialize forcing files for specific cases
    select case( trim( runtype ) )
#ifdef UNRELEASED_CODE
    case( "arm_3year" )
      call arm_3year_init( iunit, forcings_file_path )

    case( "arm_0003" )
      call arm_0003_init( iunit, forcings_file_path )

    case( "lba" )
      call lba_init( iunit, forcings_file_path )

    case( "twp_ice" )
      call twp_ice_init( iunit, forcings_file_path )
#endif

    case( "mpace_a" )
      call mpace_a_init( iunit, forcings_file_path )

    end select


    if( l_time_dependant ) then
      call initialize_time_dependant_input &
           ( iunit, runtype, gr%nnzp, gr%zt )
    end if

    wm_zm = zt2zm( wm_zt )

    wpthlp_sfc = wpthlp(1)
    wprtp_sfc  = wprtp(1)
    upwp_sfc   = upwp(1)
    vpwp_sfc   = vpwp(1)

    return
  end subroutine restart_clubb

  !----------------------------------------------------------------------
  subroutine advance_clubb_forcings( iter, dt, err_code )

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
      hydromet, radht, um_ref,  & ! Variable(s)
      vm_ref, Frad,  Frad_SW_up,  Frad_LW_up, &
      Frad_SW_down, Frad_LW_down, Ncnm, thvm, ustar, & 
      shf, Kh_zm, ug, vg

    use variables_diagnostic_module, only: wpedsclrp ! Passive scalar variables

    use variables_prognostic_module, only: &
      rtm_forcing, thlm_forcing,  & ! Variable(s)
      wm_zt, wm_zm, rho, rtm, thlm, p_in_Pa, & 
      exner, rcm, rho_zm, um, psfc, vm, & 
      upwp_sfc, vpwp_sfc, Tsfc, & 
      wpthlp_sfc, SE, LE, wprtp_sfc, cf, &
      um_forcing, vm_forcing, pdf_params

    use stats_variables, only: &
      ish, & ! Variable(s)
      ilh, &
      iwpthlp_sfc, &
      iwprtp_sfc, &
      iupwp_sfc, &
      ivpwp_sfc, &
      iustar, &
      ishf, &
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

    use variables_diagnostic_module, only:  & 
      wp2_zt ! w'^2 interpolated the themo levels [m^2/s^2]

    use stats_precision, only: time_precision ! Variable(s)

    use numerical_check, only: isnan2d, rad_check ! Procedure(s)

    use microphys_driver, only: advance_microphys ! Procedure(s)

    use parameters_microphys, only: &
      micro_scheme, l_cloud_sed, Ncm_initial  ! Variables

    use parameters_radiation, only: &
      rad_scheme

    use soil_vegetation, only: advance_soil_veg, veg_T_in_K

    use error_code, only: lapack_error,  & ! Procedure(s)
                          clubb_at_least_debug_level

    use array_index, only: & 
      iirsnowm, iiricem, iiNcm, & 
      iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl

    ! Case specific modules
    use arm, only: arm_tndcy, arm_sfclyr ! Procedure(s)

#ifdef UNRELEASED_CODE
    use arm_0003, only: arm_0003_tndcy, arm_0003_sfclyr ! Procedure(s)

    use arm_3year, only: arm_3year_tndcy, arm_3year_sfclyr ! Procedure(s)

    use arm_97, only: arm_97_tndcy, arm_97_sfclyr ! Procedure(s)

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

    use cloud_sed_mod, only: cloud_drop_sed ! Procedure(s)

    use fire, only: fire_tndcy, sfc_momentum_fluxes, & 
                    sfc_thermo_fluxes ! Procedure(s)

    use gabls2, only: gabls2_tndcy, gabls2_sfclyr ! Procedure(s)

#ifdef UNRELEASED_CODE
    use gabls3, only: gabls3_tndcy, gabls3_sfclyr ! Procedures(s)

    use rico, only: rico_tndcy, rico_sfclyr ! Procedure(s)

    use lba, only: lba_tndcy, lba_sfclyr ! Procedure(s)
#endif

    use mpace_a, only: mpace_a_tndcy, mpace_a_sfclyr ! Procedure(s)

    use mpace_b, only: mpace_b_tndcy, mpace_b_sfclyr ! Procedure(s)

#ifdef UNRELEASED_CODE
    use nov11, only: nov11_altocu_tndcy ! Procedure(s)

    use twp_ice, only: twp_ice_tndcy, twp_ice_sfclyr ! Procedure(s)

    use jun25, only: jun25_altocu_tndcy ! Procedure(s)
#endif

    use wangara, only: wangara_tndcy, wangara_sfclyr ! Procedure(s)

#ifdef radoffline
    use bugsrad_clubb_mod, only: bugsrad_clubb ! Procedure(s)
#endif
    use cos_solar_zen_mod, only: cos_solar_zen ! Function


    implicit none

    ! Input Variables
    integer, intent(in) :: iter ! Model iteration number

    real(kind=time_precision), intent(in) :: & 
      dt         ! Model timestep                            [s]

    ! Input/Output Variables
    integer, intent(inout) :: & 
      err_code

    ! Local Variables
    real, dimension(gr%nnzp) ::  & 
      rsnowm,  & ! Snow mixing ratio                         [kg/kg]
      ricem      ! Prisitine ice water mixing ratio          [kg/kg]

    real ::  &
      um_sfc, &  ! um interpolated to momentum level 1 (sfc) [m/s]
      vm_sfc     ! vm interpolated to momentum level 1 (sfc) [m/s]

    real :: wpthep, amu0

    integer :: lin_int_buffer

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
    case ( "arm_0003" ) ! ARM March 2000 case
      call arm_0003_tndcy( time_current, &                  ! Intent(in)
                           thlm_forcing,  &                 ! Intent(out)
                           rtm_forcing, um_ref, vm_ref, &   ! Intent(out)
                           sclrm_forcing, edsclrm_forcing ) ! Intent(out)

    case ( "arm_3year" ) ! ARM 3 year case
      call arm_3year_tndcy( time_current, &                  ! Intent(in)
                            thlm_forcing,  &                 ! Intent(out)
                            rtm_forcing, um_ref, vm_ref, &   ! Intent(out)
                            sclrm_forcing, edsclrm_forcing ) ! Intent(out)

    case ( "arm_97" ) ! 27 June 1997 ARM case
      call arm_97_tndcy( time_current, &                 ! Intent(in)
                         thlm_forcing,  &                ! Intent(out)
                         rtm_forcing, um_ref, vm_ref, &  ! Intent(out)
                         sclrm_forcing, edsclrm_forcing )! Intent(out)

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
      call bomex_tndcy( radht, &                        ! Intent(out)
                        thlm_forcing, rtm_forcing, &    ! Intent(out)
                        sclrm_forcing, edsclrm_forcing )! Intent(out)

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
      call dycoms2_rf02_tndcy( rho, &          ! Intent(in)
                               rho_zm, rtm, rcm, exner, &                  ! Intent(in)
                               err_code, &                                 ! Intent(inout)
                               thlm_forcing, rtm_forcing, &                ! Intent(out) 
                               Frad, radht, &                              ! Intent(out)
                               sclrm_forcing, edsclrm_forcing )            ! Intent(out)

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
    case ( "gabls3" ) ! GABLS 3 case
      call gabls3_tndcy( time_current, rtm, exner, rho, &                   ! Intent(in)
                         wm_zt, wm_zm, thlm_forcing, rtm_forcing,&          ! Intent(out)
                         um_forcing, vm_forcing, ug, vg )                   ! Intent(out)
#endif

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

    case ( "twp_ice" ) ! TWP_ICE case
      call twp_ice_tndcy( time_current, p_in_Pa, rho, thvm, &    ! Intent(in)
                           wm_zt, wm_zm, thlm_forcing,  &   ! Intent(out)
                           rtm_forcing, um_ref, vm_ref, &   ! Intent(out)
                           sclrm_forcing, edsclrm_forcing ) ! Intent(out)                   
#endif

    case ( "wangara" ) ! Wangara dry CBL
      call wangara_tndcy( wm_zt, wm_zm,  &                  ! Intent(out) 
                          thlm_forcing, rtm_forcing, &      ! Intent(out)
                          sclrm_forcing, edsclrm_forcing )  ! Intent(out)

    case default

      write(unit=fstderr,fmt=*)  & 
         "advance_clubb_forcings: Don't know how to handle " & 
         //"LS forcing for runtype: "//trim( runtype )
      stop

    end select

#ifndef UNRELEASED_CODE
    ! This is just to avoid a compiler warning -dschanen 22 Jan 2008
    um_forcing = 0.
    vm_forcing = 0.
    ug = 0.
    vg = 0.
#endif


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
      call bomex_sfclyr( um(2), vm(2), &                            ! Intent(in) 
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
      call gabls3_sfclyr( um(2), vm(2), veg_T_in_K, &             ! Intent(in)
                          thlm(2), rtm(2), gr%zt(2), exner(1) , & ! Intent(in)
                          upwp_sfc, vpwp_sfc, &                   ! Intent(out)
                          wpthlp_sfc, wprtp_sfc, ustar )     ! Intent(out)

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
                             Frad_LW_down(1), wpthep, shf )
    end if


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
           ( iter, runtype, dt, time_current, &               ! Intent(in)
             thlm, p_in_Pa, exner, rho, rho_zm, rtm, rcm, cf, & ! Intent(in) 
             wm_zt, wm_zm, Kh_zm, pdf_params, & ! Intent(in)
             wp2_zt, &                                        ! Intent(in)
             Ncnm, hydromet, &                                ! Intent(inout)
             rtm_forcing, thlm_forcing, &                     ! Intent(inout)
             err_code )                                       ! Intent(out)

      if ( lapack_error( err_code ) ) return

    end if

    if ( l_cloud_sed ) then

      call cloud_drop_sed( rcm, hydromet(:,iiNcm), rho_zm, rho, exner, & ! Intent(in)
                           rtm_forcing, thlm_forcing )              ! Intent(out)

    end if

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

        if ( isnan2d( cf ) ) then
          write(fstderr,*) "cf before BUGSrad is NaN"
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

        if ( isnan2d( thlm_forcing ) ) then
          write(fstderr,*) "thlm_forcing before BUGSrad is NaN"
        endif

        ! Check for impossible negative values
        call rad_check( thlm, rcm, rtm, ricem, &            ! Intent(in)
                        cf, p_in_Pa, exner, rho_zm )        ! Intent(in)

      endif  ! clubb_at_least_debug_level( 2 )

      ! Initially we will set this to a constant for testing purposes
      ! lin_int_buffer = 20

      ! Use a a new formula that creates and evenly spaced grid
      ! between the model domain top and the standard atmosphere
      ! table.  e.g. if the CLUBB model top is 3200m, and the spacing
      ! between gr%nnzp-1 and gr%nnzp is 40m, then lin_int_buffer is
      ! 19 and each layer of the buffer is 40m deep. -dschanen 14 May 08
      lin_int_buffer =  & 
         max( int( ( 1000.-mod( gr%zm(gr%nnzp), 1000. ) ) & 
                 * gr%dzm(gr%nnzp) ) - 1, 0 )

      !print *, "lin_int_buffer = ", lin_int_buffer !%% debug

      call bugsrad_clubb( gr%zm, gr%nnzp, lin_int_buffer,  & ! In
                          rlat, rlon,                      & ! In
                          day, month, year, time_current,  & ! In
                          thlm, rcm, rtm, rsnowm, ricem,   & ! In
                          cf, p_in_Pa, zt2zm( p_in_Pa ),   & ! In
                          exner, rho_zm,                   & ! In
                          radht, Frad,                     & ! Out
                          Frad_SW_up, Frad_LW_up,          & ! Out
                          Frad_SW_down, Frad_LW_down,      & ! Out
                          thlm_forcing )                     ! In/Out

      if ( clubb_at_least_debug_level( 2 ) ) then

        if ( isnan2d( thlm_forcing ) ) then
          write(fstderr,*) "thlm_forcing after BUGSrad is NaN"
          !write(fstderr,*) thlm_forcing
        endif

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

      call stat_update_var_pt( ishf, 1, shf, & ! intent(in)
                               sfc )           ! intent(inout)
    endif


    return
  end subroutine advance_clubb_forcings

end module clubb_driver
