!-------------------------------------------------------------------------------
! $Id$
! clubb_loss_driver.F90
!
! In-memory CLUBB-vs-benchmark loss evaluation for tuning workflows.
!
! Description:
!   This module owns the reusable loss path used by the tuner and by the
!   standalone loss executable.  It initializes CLUBB once, prepares benchmark
!   truth profiles from a converted LES stats file, runs one or more CLUBB
!   parameter columns through the case, and returns per-window, per-variable,
!   per-column profile-comparison metrics.
!
!
! Namelist format:
!   Loss settings are read from &tuner_loss_nl in the case/run namelist.
!
!   &tuner_loss_nl
!     les_stats_file = "path/to/converted_les_stats.nc"
!     clubb_var_names(1) = "cloud_frac"
!     benchmark_var_name(1) = "cloud_frac"
!     altitude_comparison_range = 0.0, 3000.0
!     time_average_range = 10800, 21600
!     num_time_windows = 1
!   /
!
!   clubb_var_names(:) are variables from the active CLUBB stats registry.
!   benchmark_var_name(:) are paired variables in les_stats_file.  Each active
!   CLUBB variable must have an explicit benchmark variable entry.
!
!
! Usage:
!
!   One-shot loss path:
!     - clubb_get_loss(): initialize CLUBB, score the default parameter matrix,
!       and finalize CLUBB before returning.
!
!     - clubb_standalone_loss calls clubb_get_loss() from an executable and
!       prints one row of metrics for each requested variable and time window.
!
!
!   Reusable tuner path:
!     - init_clubb_loss(): initialize CLUBB, read &tuner_loss_nl, snapshot the
!       stats registry and source grids, and prepare benchmark truth profiles.
!
!     - get_loss_time_window_count(): return the number of configured loss
!       windows so callers can size UI or Python-side outputs.
!
!     - clubb_get_loss_for_params(): run a full CLUBB parameter matrix and fill
!       metric arrays shaped (time_window, variable, parameter_column).
!
!     - finalize_clubb_loss(): clean up CLUBB and release prepared loss state.
!
!
!   Metric helper:
!     - calculate_taylor_metrics(): compute centered correlation, standard
!       deviation ratio, centered RMSE, and bias diagnostics for one pair of
!       profiles.  This helper is public so Python/tuner code can use the same
!       metric definitions as the Fortran loss path.
!
!
! Understanding:
!
!   Request lifecycle:
!     - init_clubb_loss() stores one module-level active_request.  That request
!       contains the parsed loss namelist, the resolved CLUBB stats variable ids,
!       the comparison levels, and benchmark truth profiles already averaged onto
!       the CLUBB source grid.
!
!     - The active request is intentionally reused across repeated
!       clubb_get_loss_for_params() calls.  This keeps benchmark NetCDF reads,
!       vertical interpolation setup, and stats-registry lookups out of the
!       inner tuning loop.
!
!   Stats dependency:
!     - The model profiles being scored come from the in-memory stats snapshot
!       returned by get_stats().  A variable must therefore be present in the
!       stats list used for the run, even if no NetCDF stats file is being
!       inspected by the caller.
!
!     - stats_get_source_grid_api() is used during initialization so benchmark
!       profiles can be interpolated to the same zt or zm grid on which CLUBB
!       accumulates stats.  The loss code compares on the runtime/source grid,
!       not on any remapped NetCDF output grid.
!
!   Benchmark preparation:
!     - les_stats_file is expected to be a stats-like benchmark file, normally
!       produced by the LES benchmark converter.  The driver opens it once at
!       initialization, reads the requested benchmark variables, converts the
!       benchmark time coordinate to seconds, and averages the requested records
!       inside each configured loss window.
!
!     - Vertical interpolation is done from the benchmark altitude coordinate to
!       the selected CLUBB grid for each requested variable.  Only CLUBB zt and
!       zm stats grids are currently handled by this loss path.
!
!   Time windows:
!     - time_average_range is interpreted as an absolute model-time interval in
!       seconds.  Benchmark records use an open/closed convention,
!       (window_start, window_end], which matches the stats sampling convention
!       used elsewhere in CLUBB.
!
!     - num_time_windows controls the loss-window breakdown.  A value of 1
!       scores one average over the full time_average_range; larger values split
!       that range into equal subwindows and return metrics for each subwindow
!       separately.
!
!     - During scoring, CLUBB is advanced only as far as the next requested loss
!       window end.  The next window continues from the previous model state.
!
!   Batching:
!     - The input parameter matrix is sized by total requested parameter columns,
!       but CLUBB may run only runtime_batch_size columns at once.  This module
!       slices the parameter matrix into runtime batches and passes batch_num to
!       the driver so stats output columns and in-memory stats columns stay
!       aligned.
!
!     - Metric arrays are always returned in total-parameter-column space.  The
!       batch loop copies each runtime batch result into the proper column range
!       of the caller's output arrays.
!
!   Metric outputs:
!     - scaled_rmse is the sum of squared model-minus-benchmark profile
!       differences after scaling by the benchmark profile range in the selected
!       altitude window.  The Taylor-style outputs are correlation, std_ratio,
!       centered_rmse_norm, and bias_norm.
!
!     - Non-finite model or benchmark values are converted to finite bad metric
!       values.  This keeps tuner ranking and JSON/NetCDF consumers from
!       inheriting NaN or Inf values.
!
!   Error behavior:
!     - Configuration errors in the loss request are treated as fatal because
!       they indicate the tuner and benchmark setup do not describe a meaningful
!       comparison.  Runtime CLUBB errors are returned through err_info so callers
!       can report the failed parameter set cleanly.
!
!-------------------------------------------------------------------------------
module clubb_loss_driver

  use clubb_driver, only: &
    init_clubb_case, &
    set_case_initial_conditions, &
    advance_clubb_to_end, &
    clean_up_clubb, &
    get_stats, &
    get_runtime_batch_config, &
    get_clubb_params_all

  use stats_netcdf, only: &
    stats_type, &
    stats_get_source_grid_api

  use input_netcdf, only: &
    open_netcdf_read, &
    get_netcdf_var, &
    close_netcdf_read

  use netcdf, only: &
    nf90_get_var, &
    nf90_get_att, &
    nf90_strerror, &
    NF90_NOERR

  use stat_file_module, only: &
    stat_file

  use interpolation, only: &
    lin_interpolate_two_points

  use clubb_precision, only: &
    core_rknd, &
    time_precision

  use constants_clubb, only: &
    fstdout, &
    fstderr

  use error_code, only: &
    clubb_fatal_error, &
    clubb_no_error

  use err_info_type_module, only: &
    err_info_type, &
    init_default_err_info_api

  use clubb_api_module, only: &
    clubb_at_least_debug_level_api

  implicit none

  private

  integer, parameter, public :: &
    loss_name_len = 64, &       ! Max stored length for CLUBB and benchmark variable names.
    max_loss_variables = 32     ! Hard cap on the number of simultaneously scored variables.

  real( kind = core_rknd ), parameter :: &
    invalid_loss_penalty = 1.0e30_core_rknd ! Finite penalty for non-finite loss-driver profiles.

  ! One CLUBB-vs-benchmark profile comparison.
  type :: loss_field_type
    character(len=loss_name_len) :: clubb_var_name = ""            ! Requested CLUBB stats variable name.
    character(len=loss_name_len) :: benchmark_var_name = ""        ! Paired benchmark variable name.
    integer :: stats_var_id = 0                                    ! Matching variable id in the CLUBB stats registry.
    integer :: k_min = 0                                           ! Lower resolved model level for this field.
    integer :: k_max = 0                                           ! Upper resolved model level for this field.
    real( kind = core_rknd ), allocatable, dimension(:,:) :: &
      truth_profile                                                ! Benchmark profiles by time window on the CLUBB grid.
  end type loss_field_type

  ! Full configuration for one loss evaluation.
  type :: loss_request_type
    integer :: num_variables = 0                                   ! Number of requested loss variables.
    integer :: time_average_range(2) = 0                           ! Absolute model-time window [s] used for averaging.
    integer :: num_time_windows = 1                                ! Number of averaged comparison subwindows.
    integer, allocatable, dimension(:,:) :: &
      time_window_ranges                                           ! Absolute model-time windows shaped (2,num_time_windows).
    integer :: total_param_sets = 0                                ! Total requested parameter columns for this prepared request.
    real( kind = core_rknd ) :: &
      dt_main_seconds = 0.0_core_rknd, &                           ! Runtime main timestep [s].
      time_initial_seconds = 0.0_core_rknd                         ! Runtime model initial time [s].
    real( kind = core_rknd ) :: &
      altitude_comparison_range(2) = 0.0_core_rknd                ! Height window [m] used in the loss.
    character(len=256) :: les_stats_file = ""                      ! Benchmark stats file to compare against.
    logical :: l_initialized = .false.                             ! True once CLUBB and this request are prepared for repeated runs.
    type(loss_field_type), allocatable, dimension(:) :: &
      fields                                                       ! Per-variable request metadata and prepared truth data.
  end type loss_request_type

  type(loss_request_type), save :: &
    active_request                                                 ! Module-owned prepared loss request reused across manual loss calls.

  public :: &
    init_clubb_loss, &
    get_loss_time_window_count, &
    clubb_get_loss_for_params, &
    finalize_clubb_loss, &
    clubb_get_loss, &
    calculate_taylor_metrics

contains

  logical function is_finite_core_value( value )

    ! Description:
    !   Return true only for finite CLUBB-core real values.  This avoids
    !   propagating NaN/Inf into tuner ranking and JSON output.

    implicit none

    real( kind = core_rknd ), intent(in) :: &
      value

    is_finite_core_value = ( value == value ) .and. ( abs( value ) <= huge( value ) )

  end function is_finite_core_value

  subroutine set_invalid_field_metric_outputs( scaled_rmse_value, correlation_value, std_ratio_value, &
                                               centered_rmse_norm_value, bias_norm_value )

    ! Description:
    !   Fill one field/column result with finite values that rank as very bad.

    implicit none

    real( kind = core_rknd ), intent(out) :: &
      scaled_rmse_value, &
      correlation_value, &
      std_ratio_value, &
      centered_rmse_norm_value, &
      bias_norm_value

    scaled_rmse_value = invalid_loss_penalty
    correlation_value = 0.0_core_rknd
    std_ratio_value = 0.0_core_rknd
    centered_rmse_norm_value = invalid_loss_penalty
    bias_norm_value = invalid_loss_penalty

  end subroutine set_invalid_field_metric_outputs

  subroutine calculate_taylor_metrics( model_profile, benchmark_profile, correlation, std_ratio, &
                                       centered_rmse_norm, bias_norm )

    ! Description:
    !   Compute Taylor-diagram profile diagnostics from one model profile and
    !   its paired benchmark profile.

    implicit none

    ! ------------------------ Inputs ------------------------
    real( kind = core_rknd ), dimension(:), intent(in) :: &
      model_profile, &                                             ! CLUBB profile over compared levels.
      benchmark_profile                                            ! Benchmark profile over compared levels.

    ! ------------------------ Outputs ------------------------
    real( kind = core_rknd ), intent(out) :: &
      correlation, &                                                ! Centered profile correlation.
      std_ratio, &                                                  ! Model profile standard deviation / benchmark standard deviation.
      centered_rmse_norm, &                                         ! Centered RMS difference normalized by benchmark stddev.
      bias_norm                                                     ! Mean bias normalized by benchmark stddev.

    ! ------------------------- Locals -------------------------
    integer :: &
      k, &                                                          ! Profile level index.
      num_levels                                                    ! Number of compared levels.

    real( kind = core_rknd ) :: &
      model_mean, &                                                 ! Model profile mean.
      benchmark_mean, &                                             ! Benchmark profile mean.
      model_centered, &                                             ! Model value with mean removed.
      benchmark_centered, &                                         ! Benchmark value with mean removed.
      model_centered_sumsq, &                                       ! Sum of squared centered model values.
      benchmark_centered_sumsq, &                                   ! Sum of squared centered benchmark values.
      covariance_sum, &                                             ! Sum of centered model times centered benchmark.
      centered_diff_sumsq, &                                        ! Sum of squared centered-profile differences.
      model_stddev, &                                               ! Model profile standard deviation.
      benchmark_stddev, &                                           ! Benchmark profile standard deviation.
      centered_rmse, &                                              ! Centered profile RMS difference.
      bias                                                          ! Model mean minus benchmark mean.

    if ( size( model_profile ) /= size( benchmark_profile ) ) then
      call stop_with_error( "Taylor metric profiles must have matching sizes" )
    end if

    num_levels = size( model_profile )
    if ( num_levels <= 0 ) then
      call stop_with_error( "Taylor metric profiles must contain at least one level" )
    end if

    model_mean = sum( model_profile ) / real( num_levels, kind = core_rknd )
    benchmark_mean = sum( benchmark_profile ) / real( num_levels, kind = core_rknd )

    model_centered_sumsq = 0.0_core_rknd
    benchmark_centered_sumsq = 0.0_core_rknd
    covariance_sum = 0.0_core_rknd
    centered_diff_sumsq = 0.0_core_rknd

    do k = 1, num_levels
      model_centered = model_profile(k) - model_mean
      benchmark_centered = benchmark_profile(k) - benchmark_mean
      model_centered_sumsq = model_centered_sumsq + model_centered**2
      benchmark_centered_sumsq = benchmark_centered_sumsq + benchmark_centered**2
      covariance_sum = covariance_sum + model_centered * benchmark_centered
      centered_diff_sumsq = centered_diff_sumsq + ( model_centered - benchmark_centered )**2
    end do

    model_stddev = sqrt( model_centered_sumsq / real( num_levels, kind = core_rknd ) )
    benchmark_stddev = sqrt( benchmark_centered_sumsq / real( num_levels, kind = core_rknd ) )
    centered_rmse = sqrt( centered_diff_sumsq / real( num_levels, kind = core_rknd ) )
    bias = model_mean - benchmark_mean

    if ( model_stddev > 0.0_core_rknd .and. benchmark_stddev > 0.0_core_rknd ) then
      correlation = covariance_sum / sqrt( model_centered_sumsq * benchmark_centered_sumsq )
      correlation = max( -1.0_core_rknd, min( 1.0_core_rknd, correlation ) )
    else
      correlation = 0.0_core_rknd
    end if

    if ( benchmark_stddev > 0.0_core_rknd ) then
      std_ratio = model_stddev / benchmark_stddev
      centered_rmse_norm = centered_rmse / benchmark_stddev
      bias_norm = bias / benchmark_stddev
    else
      std_ratio = 0.0_core_rknd
      centered_rmse_norm = centered_rmse
      bias_norm = bias
    end if

  end subroutine calculate_taylor_metrics

  subroutine stop_with_error( message )

    ! Description:
    !   Write a fatal message and stop.

    implicit none

    character(len=*), intent(in) :: &
      message

    write( fstderr, '(a)' ) trim( message )
    stop 1

  end subroutine stop_with_error

  subroutine init_loss_request( runfile, total_param_sets, request, requested_clubb_var_names )

    ! Description:
    !   Read tuner_loss_nl, validate it, prepare benchmark truth profiles, and
    !   return the reusable loss-request state.

    implicit none

    ! ------------------------ Inputs ------------------------
    character(len=*), intent(in) :: &
      runfile                                                      ! Aggregate runfile that contains tuner_loss_nl.

    integer, intent(in) :: &
      total_param_sets                                             ! Total number of requested parameter-set columns.

    ! ------------------------ Outputs ------------------------
    type(loss_request_type), intent(out) :: &
      request                                                      ! Normalized loss request state.

    character(len=loss_name_len), allocatable, dimension(:), intent(out) :: &
      requested_clubb_var_names                                    ! Requested CLUBB variable names in printed row order.

    ! ------------------------- Locals -------------------------
    character(len=256) :: &
      les_stats_file                                               ! Benchmark stats file to compare against.

    character(len=loss_name_len), dimension(max_loss_variables) :: &
      clubb_var_names, &                                           ! Raw CLUBB variable names from the namelist.
      benchmark_var_name                                           ! Raw paired benchmark variable names.

    real( kind = core_rknd ) :: &
      altitude_comparison_range(2)                                ! Requested height window [m].

    integer :: &
      time_average_range(2), &                                     ! Requested absolute averaging window [s].
      num_time_windows, &                                          ! Number of averaged comparison subwindows.
      window_width, &                                              ! Equal split-window width [s].
      iunit_local, &                                               ! Local unit used to scan for tuner_loss_nl.
      ios, &                                                       ! I/O status while scanning or reading the namelist.
      i, &                                                         ! Loop index over requested variables.
      j, &                                                         ! Loop index for duplicate-request checks.
      nvars                                                        ! Number of active nonblank variable entries.

    logical :: &
      l_found, &                                                   ! True once tuner_loss_nl is found and read.
      l_seen_blank                                                 ! True after the first blank CLUBB variable entry.

    character(len=512) :: &
      line, &                                                      ! Raw runfile line used while scanning for tuner_loss_nl.
      iomsg                                                        ! Namelist read failure text.

    type(stats_type) :: &
      stats_snapshot                                               ! Snapshot of the configured CLUBB stats registry.

    namelist /tuner_loss_nl/ les_stats_file, clubb_var_names, benchmark_var_name, &
                             altitude_comparison_range, time_average_range, &
                             num_time_windows

    les_stats_file = ""
    clubb_var_names = ""
    benchmark_var_name = ""
    altitude_comparison_range = 0.0_core_rknd
    time_average_range = 0
    num_time_windows = 1

    iunit_local = 907
    l_found = .false.

    ! Find the tuner-specific namelist inside the aggregate runfile.
    open( unit = iunit_local, file = trim( runfile ), status = "old", action = "read" )
    do
      read( iunit_local, '(A)', iostat = ios ) line
      if ( ios /= 0 ) exit
      if ( index( adjustl( line ), "&tuner_loss_nl" ) == 1 .or. &
           index( adjustl( line ), "&TUNER_LOSS_NL" ) == 1 ) then
        backspace( iunit_local )
        iomsg = ""
        ! Read the tuner request once the right namelist block is found.
        read( iunit_local, nml = tuner_loss_nl, iostat = ios, iomsg = iomsg )
        if ( ios /= 0 ) then
          close( iunit_local )
          call stop_with_error( "Failed to read &tuner_loss_nl from "//trim( runfile )//"; iomsg="//trim( iomsg ) )
        end if
        l_found = .true.
        exit
      end if
    end do
    close( iunit_local )

    if ( .not. l_found ) then
      call stop_with_error( "Missing &tuner_loss_nl in "//trim( runfile ) )
    end if

    ! Request validation
    ! Reject malformed high-level loss input before touching CLUBB state.
    if ( len_trim( les_stats_file ) == 0 ) then
      call stop_with_error( "tuner_loss_nl requires les_stats_file" )
    end if

    if ( altitude_comparison_range(1) < 0.0_core_rknd ) then
      call stop_with_error( "tuner_loss_nl requires altitude_comparison_range(1) >= 0" )
    end if

    if ( altitude_comparison_range(2) < altitude_comparison_range(1) ) then
      call stop_with_error( "tuner_loss_nl requires altitude_comparison_range(2) >= altitude_comparison_range(1)" )
    end if

    if ( time_average_range(1) < 0 .or. time_average_range(2) <= time_average_range(1) ) then
      call stop_with_error( "tuner_loss_nl requires time_average_range(1) >= 0 and time_average_range(2) > time_average_range(1)" )
    end if

    if ( num_time_windows < 1 ) then
      call stop_with_error( "tuner_loss_nl requires num_time_windows >= 1" )
    end if
    if ( mod( time_average_range(2) - time_average_range(1), num_time_windows ) /= 0 ) then
      call stop_with_error( "time_average_range must divide evenly by num_time_windows" )
    end if

    ! Collapse the fixed-size namelist arrays down to the active requests.
    nvars = 0
    l_seen_blank = .false.
    do i = 1, size( clubb_var_names )
      if ( len_trim( clubb_var_names(i) ) == 0 ) then
        l_seen_blank = .true.
        if ( len_trim( benchmark_var_name(i) ) /= 0 ) then
          call stop_with_error( "benchmark_var_name must not contain nonblank entries after the first blank entry" )
        end if
        cycle
      end if

      if ( l_seen_blank ) then
        call stop_with_error( "clubb_var_names must not contain nonblank entries after the first blank entry" )
      end if

      if ( len_trim( benchmark_var_name(i) ) == 0 ) then
        call stop_with_error( "Each CLUBB variable in clubb_var_names must have a paired benchmark_var_name entry" )
      end if

      do j = 1, i - 1
        if ( trim( clubb_var_names(i) ) == trim( clubb_var_names(j) ) ) then
          call stop_with_error( "Duplicate CLUBB variable requested: "//trim( clubb_var_names(i) ) )
        end if
      end do

      nvars = nvars + 1
    end do

    if ( nvars == 0 ) then
      call stop_with_error( "clubb_var_names must contain at least one CLUBB variable" )
    end if

    ! Store only the active request entries.
    request%les_stats_file = trim( les_stats_file )
    request%altitude_comparison_range = altitude_comparison_range
    request%time_average_range = time_average_range
    request%num_time_windows = num_time_windows
    request%num_variables = nvars
    request%total_param_sets = total_param_sets
    request%l_initialized = .false.

    allocate( request%fields(nvars) )
    allocate( request%time_window_ranges(2,num_time_windows) )
    window_width = ( time_average_range(2) - time_average_range(1) ) / num_time_windows
    do i = 1, num_time_windows
      request%time_window_ranges(1,i) = time_average_range(1) + ( i - 1 ) * window_width
      request%time_window_ranges(2,i) = time_average_range(1) + i * window_width
    end do

    ! Store each requested comparison in one field descriptor.
    do i = 1, nvars
      request%fields(i)%clubb_var_name = trim( adjustl( clubb_var_names(i) ) )
      request%fields(i)%benchmark_var_name = trim( adjustl( benchmark_var_name(i) ) )

    end do

    ! CLUBB-dependent preparation
    ! Snapshot the configured stats registry used for name lookup and scoring.
    call get_stats( stats_snapshot )
    request%dt_main_seconds = real( stats_snapshot%dt_main, kind = core_rknd )
    request%time_initial_seconds = real( stats_snapshot%time_initial, kind = core_rknd )

    ! Finish the CLUBB-dependent request preparation before batching starts.
    call prepare_loss_request_for_scoring( request, stats_snapshot )

    ! Return row labels in request order.
    allocate( requested_clubb_var_names( request%num_variables ) )
    requested_clubb_var_names(:) = request%fields(:)%clubb_var_name
    request%l_initialized = .true.

  end subroutine init_loss_request

  subroutine init_clubb_loss( runfile, clubb_var_names, err_info, clubb_params_all )

    ! Description:
    !   Initialize CLUBB once, prepare the reusable loss request, and
    !   optionally return the full default parameter matrix.

    implicit none

    character(len=*), intent(in) :: &
      runfile                                                      ! Aggregated runfile used to initialize CLUBB.

    character(len=loss_name_len), allocatable, dimension(:), intent(out) :: &
      clubb_var_names                                              ! Requested CLUBB variable names, in printed loss-row order.

    type(err_info_type), intent(out) :: &
      err_info                                                     ! CLUBB/runtime error state for the caller.

    real( kind = core_rknd ), allocatable, dimension(:,:), optional, intent(out) :: &
      clubb_params_all                                             ! Full configured parameter matrix returned for manual reruns.

    integer :: &
      total_param_sets, &                                          ! Total number of requested parameter sets.
      runtime_batch_size                                           ! Number of columns advanced in one CLUBB batch.

    ! Phase 1: initialize CLUBB and prepare the loss request.
    if ( active_request%l_initialized ) then
      call stop_with_error( "init_clubb_loss was called while a prepared loss request is still active" )
    end if

    call init_default_err_info_api( 1, err_info )
    call init_clubb_case( trim( runfile ), err_info )

    if ( clubb_at_least_debug_level_api( 0 ) ) then
      if ( any( err_info%err_code == clubb_fatal_error ) ) then
        call stop_with_error( "init_clubb_case failed during init_clubb_loss" )
      end if
    end if

    ! Size the prepared request from the initialized CLUBB runtime.
    call get_runtime_batch_config( total_param_sets, runtime_batch_size )
    if ( runtime_batch_size <= 0 ) then
      call stop_with_error( "init_clubb_loss requires a positive runtime batch size" )
    end if

    ! Build the prepared request after CLUBB init.
    call init_loss_request( trim( runfile ), total_param_sets, active_request, clubb_var_names )

    if ( present( clubb_params_all ) ) then
      call get_clubb_params_all( clubb_params_all )
    end if

  end subroutine init_clubb_loss

  subroutine get_loss_time_window_count( num_time_windows )

    ! Description:
    !   Return the time-window count for the active prepared loss request.

    implicit none

    integer, intent(out) :: &
      num_time_windows

    if ( active_request%l_initialized ) then
      num_time_windows = active_request%num_time_windows
    else
      num_time_windows = 1
    end if

  end subroutine get_loss_time_window_count

  subroutine clubb_get_loss_for_params( clubb_params_all, scaled_rmse, correlation, std_ratio, &
                                        centered_rmse_norm, bias_norm, err_info )

    ! Description:
    !   Reuse an initialized CLUBB/loss setup to score one full parameter
    !   matrix, potentially advancing it in runtime-sized sub-batches.

    implicit none

    real( kind = core_rknd ), dimension(:,:), intent(in) :: &
      clubb_params_all                                             ! Full parameter matrix to score against the prepared request.

    real( kind = core_rknd ), allocatable, dimension(:,:,:), intent(out) :: &
      scaled_rmse, &                                               ! Scaled RMSE diagnostic shaped (window, variable, parameter set).
      correlation, &                                               ! Taylor correlation shaped (window, variable, parameter set).
      std_ratio, &                                                 ! Taylor standard-deviation ratio shaped (window, variable, parameter set).
      centered_rmse_norm, &                                        ! Taylor centered RMSE shaped (window, variable, parameter set).
      bias_norm                                                    ! Taylor normalized bias shaped (window, variable, parameter set).

    type(err_info_type), intent(out) :: &
      err_info                                                     ! CLUBB/runtime error state for the caller.

    type(stats_type) :: &
      stats_snapshot                                               ! Snapshot of the completed final stats window.

    integer :: &
      total_param_sets, &                                          ! Total number of requested parameter sets.
      runtime_batch_size, &                                        ! Number of columns advanced in one CLUBB batch.
      num_batches, &                                               ! Number of reruns needed to cover all parameter sets.
      batch_idx, &                                                 ! Current batch loop index.
      window_idx, &                                                ! Current time-window index.
      itime_start, &                                               ! First timestep for this advance chunk.
      itime_end, &                                                 ! Last timestep for this advance chunk.
      batch_start, &                                               ! First full-matrix column covered by the active batch.
      batch_end                                                    ! Last full-matrix column covered by the active batch.

    if ( .not. active_request%l_initialized ) then
      call stop_with_error( "clubb_get_loss_for_params requires init_clubb_loss to be called first" )
    end if

    call get_runtime_batch_config( total_param_sets, runtime_batch_size )
    call init_default_err_info_api( runtime_batch_size, err_info )

    if ( total_param_sets /= active_request%total_param_sets ) then
      call stop_with_error( "Initialized CLUBB column count does not match the prepared loss request" )
    end if

    if ( size( clubb_params_all, 1 ) /= active_request%total_param_sets ) then
      call stop_with_error( "clubb_get_loss_for_params requires the full parameter matrix to match total_param_sets" )
    end if

    allocate( scaled_rmse( active_request%num_time_windows, active_request%num_variables, &
                           active_request%total_param_sets ) )
    allocate( correlation( active_request%num_time_windows, active_request%num_variables, &
                           active_request%total_param_sets ) )
    allocate( std_ratio( active_request%num_time_windows, active_request%num_variables, &
                         active_request%total_param_sets ) )
    allocate( centered_rmse_norm( active_request%num_time_windows, active_request%num_variables, &
                                  active_request%total_param_sets ) )
    allocate( bias_norm( active_request%num_time_windows, active_request%num_variables, &
                         active_request%total_param_sets ) )

    scaled_rmse = 0.0_core_rknd
    correlation = 0.0_core_rknd
    std_ratio = 0.0_core_rknd
    centered_rmse_norm = 0.0_core_rknd
    bias_norm = 0.0_core_rknd

    num_batches = total_param_sets / runtime_batch_size

    ! Run each parameter batch through the full case.
    do batch_idx = 1, num_batches
      batch_start = ( batch_idx - 1 ) * runtime_batch_size + 1
      batch_end = batch_start + runtime_batch_size - 1

      err_info%err_code = clubb_no_error

      ! Reset the state and activate the current parameter batch.
      call set_case_initial_conditions( err_info, clubb_params_in = clubb_params_all(batch_start:batch_end,:), &
                                        batch_num = batch_idx )

      if ( clubb_at_least_debug_level_api( 0 ) ) then
        if ( any( err_info%err_code == clubb_fatal_error ) ) then
          call stop_with_error( "set_case_initial_conditions failed during clubb_get_loss_for_params" )
        end if
      end if

      itime_start = 1
      do window_idx = 1, active_request%num_time_windows
        itime_end = nint( ( real( active_request%time_window_ranges(2,window_idx), kind = core_rknd ) - &
                            active_request%time_initial_seconds ) / active_request%dt_main_seconds )

        ! Advance through this subwindow for the current batch.
        call advance_clubb_to_end( .false., err_info, itime_start = itime_start, itime_end = itime_end )

        ! Snapshot and score this completed stats subwindow.
        call get_stats( stats_snapshot )
        call calculate_field_loss( active_request, stats_snapshot, window_idx, &
                                   scaled_rmse(window_idx,:,batch_start:batch_end), &
                                   correlation(window_idx,:,batch_start:batch_end), &
                                   std_ratio(window_idx,:,batch_start:batch_end), &
                                   centered_rmse_norm(window_idx,:,batch_start:batch_end), &
                                   bias_norm(window_idx,:,batch_start:batch_end) )

        itime_start = itime_end + 1
      end do
    end do

  end subroutine clubb_get_loss_for_params

  subroutine finalize_clubb_loss( err_info )

    ! Description:
    !   Release CLUBB state and clear the prepared module-local loss request.

    implicit none

    type(err_info_type), intent(out) :: &
      err_info                                                     ! CLUBB/runtime error state for the caller.

    integer :: &
      field_idx                                                    ! Loop index used to release prepared truth profiles.

    call init_default_err_info_api( 1, err_info )
    if ( active_request%l_initialized ) then
      call clean_up_clubb( err_info )
    end if

    if ( allocated( active_request%fields ) ) then
      do field_idx = 1, size( active_request%fields )
        if ( allocated( active_request%fields(field_idx)%truth_profile ) ) then
          deallocate( active_request%fields(field_idx)%truth_profile )
        end if
      end do
      deallocate( active_request%fields )
    end if
    if ( allocated( active_request%time_window_ranges ) ) then
      deallocate( active_request%time_window_ranges )
    end if

    active_request%num_variables = 0
    active_request%time_average_range = 0
    active_request%num_time_windows = 1
    active_request%total_param_sets = 0
    active_request%altitude_comparison_range = 0.0_core_rknd
    active_request%les_stats_file = ""
    active_request%l_initialized = .false.

  end subroutine finalize_clubb_loss

  subroutine clubb_get_loss( runfile, clubb_var_names, scaled_rmse, correlation, std_ratio, &
                             centered_rmse_norm, bias_norm, err_info )

    ! Description:
    !   Initialize CLUBB, run all batches for the default parameter matrix, and
    !   clean up after the one-shot loss evaluation.

    implicit none

    character(len=*), intent(in) :: &
      runfile                                                      ! Aggregated runfile used to initialize CLUBB.

    character(len=loss_name_len), allocatable, dimension(:), intent(out) :: &
      clubb_var_names                                              ! Requested CLUBB variable names, in printed loss-row order.

    real( kind = core_rknd ), allocatable, dimension(:,:,:), intent(out) :: &
      scaled_rmse, &                                               ! Scaled RMSE diagnostic shaped (window, variable, parameter set).
      correlation, &                                               ! Taylor correlation shaped (window, variable, parameter set).
      std_ratio, &                                                 ! Taylor standard-deviation ratio shaped (window, variable, parameter set).
      centered_rmse_norm, &                                        ! Taylor centered RMSE shaped (window, variable, parameter set).
      bias_norm                                                    ! Taylor normalized bias shaped (window, variable, parameter set).

    type(err_info_type), intent(out) :: &
      err_info                                                     ! CLUBB/runtime error state for the caller.

    real( kind = core_rknd ), allocatable, dimension(:,:) :: &
      clubb_params_all                                             ! Full default parameter matrix used by the one-shot wrapper.

    ! Phase 1: initialize CLUBB and prepare the reusable loss request.
    call init_clubb_loss( trim( runfile ), clubb_var_names, err_info, clubb_params_all )
    ! Phase 2: run the default parameter matrix through the reusable loss path.
    call clubb_get_loss_for_params( clubb_params_all, scaled_rmse, correlation, std_ratio, &
                                    centered_rmse_norm, bias_norm, err_info )
    ! Phase 3: release CLUBB state after the one-shot evaluation.
    call finalize_clubb_loss( err_info )

  end subroutine clubb_get_loss

  subroutine prepare_loss_request_for_scoring( request, stats_snapshot )

    ! Description:
    !   Resolve stats bindings, model levels, and benchmark truth profiles for
    !   all requested fields.

    implicit none

    type(loss_request_type), intent(inout) :: &
      request                                                      ! Request enriched in-place with stats ids, bounds, and truth profiles.

    type(stats_type), intent(in) :: &
      stats_snapshot                                               ! Snapshot of the configured CLUBB stats registry.

    type(stat_file) :: &
      benchmark_file                                               ! Open benchmark file handle reused for all requested fields.

    integer :: &
      field_idx, &                                                 ! Loop index over requested loss variables.
      stats_var_id, &                                              ! Matching stats registry entry for the active variable.
      nz, &                                                        ! Active vertical size of the requested field.
      benchmark_time_range(2), &                                   ! Benchmark record interval matching the requested time window.
      window_idx, &                                                ! Time-window loop index.
      k, &                                                         ! Vertical loop index used for CLUBB bounds and interpolation.
      j, &                                                         ! Benchmark vertical loop index used to build interpolation stencils.
      itime, &                                                     ! Benchmark time-coordinate / profile loop index.
      stats_idx, &                                                 ! Stats-registry loop index used for name lookup.
      ierr                                                         ! NetCDF status code while reading benchmark metadata.

    character(len=16) :: &
      grid_name                                                    ! Native grid name declared by the stats registry.

    character(len=80) :: &
      time_units                                                   ! Raw time-units attribute from the benchmark time coordinate.

    real( kind = core_rknd ), pointer :: &
      active_grid(:)                                               ! Native CLUBB grid used for this field's interpolation and bounds.

    real( kind = core_rknd ), allocatable, target, dimension(:,:) :: &
      runtime_zt, &                                                ! Initialized CLUBB thermodynamic grid heights.
      runtime_zm                                                   ! Initialized CLUBB momentum grid heights.

    real( kind = core_rknd ), allocatable, dimension(:) :: &
      raw_time_values, &                                           ! Native benchmark time-coordinate values.
      benchmark_time_values_seconds, &                             ! Benchmark record times converted to seconds.
      benchmark_profile                                            ! One benchmark profile on the native benchmark grid.

    integer, allocatable, dimension(:) :: &
      lower_idx, &                                                 ! Lower benchmark interpolation levels.
      upper_idx                                                    ! Upper benchmark interpolation levels.

    logical :: &
      l_error                                                      ! Benchmark file open failure flag.

    real( kind = core_rknd ) :: &
      multiplier, &                                                ! Unit conversion factor from the benchmark time coordinate to seconds.
      requested_time_range_seconds(2), &                           ! Requested absolute benchmark time window [s].
      z_min, &                                                     ! Requested lower comparison height [m].
      z_max                                                        ! Requested upper comparison height [m].

    ! Use the native CLUBB grids for height-range resolution and interpolation.
    call stats_get_source_grid_api( stats_snapshot, runtime_zt, runtime_zm )

    ! Benchmark time setup
    ! Open the benchmark dataset once for the whole request.
    call open_netcdf_read( trim( request%fields(1)%benchmark_var_name ), trim( request%les_stats_file ), &
                           benchmark_file, l_error )
    if ( l_error ) then
      call stop_with_error( "Failed to open benchmark file: " // trim( request%les_stats_file ) )
    end if

    multiplier = 0.0_core_rknd
    benchmark_time_range = 0

    ! Convert the benchmark time coordinate to seconds.
    ierr = nf90_get_att( ncid = benchmark_file%iounit, varid = benchmark_file%TimeVarId, &
                         name = "units", values = time_units )
    if ( ierr /= NF90_NOERR ) then
      call stop_with_error( "Failed to read benchmark time units: "//trim( nf90_strerror( ierr ) ) )
    end if

    select case ( time_units( 1:index( time_units, ' ' ) ) )
    case ( "hours" )
      multiplier = 3600.0_core_rknd
    case ( "minutes" )
      multiplier = 60.0_core_rknd
    case ( "seconds" )
      multiplier = 1.0_core_rknd
    case default
      call stop_with_error( "Benchmark time units are unsupported: "//trim( time_units ) )
    end select

    allocate( raw_time_values( benchmark_file%ntimes ) )
    allocate( benchmark_time_values_seconds( benchmark_file%ntimes ) )

    ierr = nf90_get_var( ncid = benchmark_file%iounit, varid = benchmark_file%TimeVarId, &
                         start = (/ 1 /), count = (/ benchmark_file%ntimes /), values = raw_time_values )
    if ( ierr /= NF90_NOERR ) then
      call stop_with_error( "Failed to read benchmark time coordinate: "//trim( nf90_strerror( ierr ) ) )
    end if

    benchmark_time_values_seconds = raw_time_values * multiplier

    z_min = request%altitude_comparison_range(1)
    z_max = request%altitude_comparison_range(2)

    ! Per-field setup
    ! Prepare each requested field once before batching starts.
    do field_idx = 1, request%num_variables
      associate( field => request%fields(field_idx) )

        ! Resolve the requested CLUBB name to one stats field and native grid.
        stats_var_id = 0
        do stats_idx = 1, stats_snapshot%nvars
          if ( trim( stats_snapshot%vars(stats_idx)%name ) == trim( field%clubb_var_name ) ) then
            stats_var_id = stats_idx
            exit
          end if
        end do

        if ( stats_var_id == 0 ) then
          call stop_with_error( "Requested stats variable was not found for "//trim( field%clubb_var_name ) )
        end if

        field%stats_var_id = stats_var_id
        grid_name = trim( adjustl( stats_snapshot%vars(field%stats_var_id)%grid ) )
        nz = stats_snapshot%vars(field%stats_var_id)%nz

        ! Pick the native CLUBB grid used by this field.
        select case ( trim( grid_name ) )
        case ( "zt" )
          active_grid => runtime_zt(1,1:nz)
        case ( "zm" )
          active_grid => runtime_zm(1,1:nz)
        case default
          call stop_with_error( "Requested CLUBB variable must live on the zt or zm grid: "// &
                                trim( field%clubb_var_name ) )
        end select

        field%k_min = 0
        field%k_max = 0

        ! Convert the requested height range into model levels.
        do k = 1, size( active_grid )
          if ( active_grid(k) >= z_min ) then
            field%k_min = k
            exit
          end if
        end do

        if ( field%k_min == 0 ) then
          call stop_with_error( "Requested lower height bound is above the model top" )
        end if

        do k = size( active_grid ), 1, -1
          if ( active_grid(k) <= z_max ) then
            field%k_max = k
            exit
          end if
        end do

        if ( field%k_max == 0 ) then
          call stop_with_error( "Requested upper height bound is below the model bottom" )
        end if

        if ( field%k_max < field%k_min ) then
          call stop_with_error( "Requested height window does not include any model levels" )
        end if

        allocate( field%truth_profile( request%num_time_windows, nz ) )
        field%truth_profile = 0.0_core_rknd

        if ( active_grid( field%k_min ) < benchmark_file%z( benchmark_file%ia ) ) then
          call stop_with_error( "Requested lower height bound is below the benchmark domain for "// &
                                trim( field%benchmark_var_name ) )
        end if

        if ( active_grid( field%k_max ) > benchmark_file%z( benchmark_file%iz ) ) then
          call stop_with_error( "Requested upper height bound is above the benchmark domain for "// &
                                trim( field%benchmark_var_name ) )
        end if

        allocate( benchmark_profile( benchmark_file%iz ) )
        allocate( lower_idx( size( active_grid ) ), upper_idx( size( active_grid ) ) )

        lower_idx = 0
        upper_idx = 0

        ! Build the vertical interpolation stencil once for this field.
        do k = field%k_min, field%k_max
          do j = 1, size( benchmark_file%z )
            if ( abs( benchmark_file%z(j) - active_grid(k) ) < 1.e-6_core_rknd ) then
              lower_idx(k) = j
              upper_idx(k) = j
              exit
            else if ( benchmark_file%z(j) < active_grid(k) ) then
              lower_idx(k) = j
            else
              upper_idx(k) = j
              exit
            end if
          end do
        end do

        do window_idx = 1, request%num_time_windows
          benchmark_time_range = 0
          requested_time_range_seconds(1) = real( request%time_window_ranges(1,window_idx), kind = core_rknd )
          requested_time_range_seconds(2) = real( request%time_window_ranges(2,window_idx), kind = core_rknd )

          do itime = 1, size( benchmark_time_values_seconds )
            if ( benchmark_time_values_seconds(itime) > requested_time_range_seconds(1) + 1.0e-6_core_rknd &
                 .and. benchmark_time_range(1) == 0 ) then
              benchmark_time_range(1) = itime
            end if
            if ( benchmark_time_values_seconds(itime) <= requested_time_range_seconds(2) + 1.0e-6_core_rknd ) then
              benchmark_time_range(2) = itime
            end if
          end do

          if ( any( benchmark_time_range <= 0 ) .or. benchmark_time_range(2) < benchmark_time_range(1) ) then
            call stop_with_error( "Requested time subwindow does not align with the benchmark time coordinate" )
          end if

          if ( field_idx == 1 ) then
            write( fstdout, * ) "benchmark timing:", "window", window_idx, &
                                "records", benchmark_time_range(1), benchmark_time_range(2), &
                                "absolute_window", request%time_window_ranges(1,window_idx), &
                                request%time_window_ranges(2,window_idx), &
                                "benchmark_window", requested_time_range_seconds(1), requested_time_range_seconds(2), &
                                "matched_window", benchmark_time_values_seconds( benchmark_time_range(1) ), &
                                benchmark_time_values_seconds( benchmark_time_range(2) )
          end if

          ! Accumulate this benchmark time window directly on the CLUBB grid.
          do itime = benchmark_time_range(1), benchmark_time_range(2)
            call get_netcdf_var( benchmark_file, trim( field%benchmark_var_name ), itime, .true., &
                                 benchmark_profile, l_error )
            if ( l_error ) then
              call stop_with_error( "Failed to read benchmark variable "//trim( field%benchmark_var_name ) )
            end if

            do k = field%k_min, field%k_max
              if ( lower_idx(k) == upper_idx(k) ) then
                field%truth_profile(window_idx,k) = field%truth_profile(window_idx,k) + benchmark_profile( lower_idx(k) )
              else
                field%truth_profile(window_idx,k) = field%truth_profile(window_idx,k) + &
                  lin_interpolate_two_points( active_grid(k), benchmark_file%z( upper_idx(k) ), &
                                              benchmark_file%z( lower_idx(k) ), benchmark_profile( upper_idx(k) ), &
                                              benchmark_profile( lower_idx(k) ) )
              end if
            end do
          end do

          field%truth_profile( window_idx, field%k_min:field%k_max ) = &
            field%truth_profile( window_idx, field%k_min:field%k_max ) / &
            real( benchmark_time_range(2) - benchmark_time_range(1) + 1, kind = core_rknd )
        end do

        deallocate( benchmark_profile, lower_idx, upper_idx )

      end associate
    end do

    ! Close the shared benchmark file.
    call close_netcdf_read( benchmark_file )

  end subroutine prepare_loss_request_for_scoring

  subroutine calculate_field_loss( request, stats_snapshot, window_idx, scaled_rmse, correlation, std_ratio, &
                                   centered_rmse_norm, bias_norm )

    ! Description:
    !   Compute loss and Taylor diagnostics for one completed CLUBB stats window.

    implicit none

    ! ------------------------ Inputs ------------------------
    type(loss_request_type), intent(in) :: &
      request                                                      ! Fully validated loss request metadata.

    type(stats_type), intent(in) :: &
      stats_snapshot                                               ! Completed final CLUBB stats window.

    integer, intent(in) :: &
      window_idx                                                   ! Prepared comparison time-window index.

    ! ------------------------ Outputs ------------------------
    real( kind = core_rknd ), dimension(:,:), intent(out) :: &
      scaled_rmse, &                                               ! Scaled RMSE shaped (variable, active batch column).
      correlation, &                                               ! Taylor correlation shaped (variable, active batch column).
      std_ratio, &                                                 ! Taylor standard-deviation ratio shaped (variable, active batch column).
      centered_rmse_norm, &                                        ! Taylor centered RMSE shaped (variable, active batch column).
      bias_norm                                                    ! Taylor normalized bias shaped (variable, active batch column).

    ! ------------------------- Locals -------------------------
    integer :: &
      row_idx, &                                                   ! Loss-matrix row / requested-variable index.
      i, &                                                         ! Parameter-column index.
      k, &                                                         ! Vertical level index.
      level_idx, &                                                 ! Compact profile index over the compared height window.
      num_levels                                                   ! Number of vertical levels compared for one field.

    real( kind = core_rknd ) :: &
      diff, &                                                      ! Columnwise CLUBB minus benchmark difference.
      norm, &                                                      ! Benchmark variability over the compared height range.
      model_value, &                                               ! CLUBB profile value for one compared level.
      truth_value, &                                               ! Benchmark profile value for one compared level.
      field_correlation, &                                         ! Centered model/benchmark profile correlation.
      field_std_ratio, &                                           ! Model stddev normalized by benchmark stddev.
      field_centered_rmse_norm, &                                  ! Centered RMS difference normalized by benchmark stddev.
      field_bias_norm                                              ! Mean bias normalized by benchmark stddev.

    logical :: &
      l_invalid_profile                                            ! True when one compared value is NaN/Inf.

    real( kind = core_rknd ), allocatable, dimension(:) :: &
      model_profile, &                                             ! Compact CLUBB profile over compared levels.
      benchmark_profile                                            ! Compact benchmark profile over compared levels.

    scaled_rmse = 0.0_core_rknd
    correlation = 0.0_core_rknd
    std_ratio = 0.0_core_rknd
    centered_rmse_norm = 0.0_core_rknd
    bias_norm = 0.0_core_rknd

    ! Score each prepared field independently.
    do row_idx = 1, request%num_variables
      associate( field => request%fields(row_idx) )

        num_levels = field%k_max - field%k_min + 1
        norm = maxval( field%truth_profile(window_idx,field%k_min:field%k_max) ) - &
               minval( field%truth_profile(window_idx,field%k_min:field%k_max) )
        if ( .not. is_finite_core_value( norm ) .or. norm <= 0.0_core_rknd ) norm = 1.0_core_rknd

        allocate( model_profile(num_levels), benchmark_profile(num_levels) )
        benchmark_profile = field%truth_profile(window_idx,field%k_min:field%k_max)

        ! Sum the profile mismatch for each parameter column.
        do i = 1, size( scaled_rmse, 2 )
          l_invalid_profile = .false.

          ! Collapse the requested height window for this one column.
          do k = field%k_min, field%k_max
            level_idx = k - field%k_min + 1
            model_value = real( stats_snapshot%vars(field%stats_var_id)%buffer(i,k), kind = core_rknd )
            truth_value = field%truth_profile(window_idx,k)
            model_profile(level_idx) = model_value

            if ( .not. is_finite_core_value( model_value ) .or. &
                 .not. is_finite_core_value( truth_value ) ) then
              l_invalid_profile = .true.
            else
              diff = model_value - truth_value
              scaled_rmse(row_idx,i) = scaled_rmse(row_idx,i) + ( diff / norm )**2
              if ( .not. is_finite_core_value( scaled_rmse(row_idx,i) ) ) l_invalid_profile = .true.
            end if
          end do

          if ( l_invalid_profile ) then
            call set_invalid_field_metric_outputs( scaled_rmse(row_idx,i), correlation(row_idx,i), &
                                                   std_ratio(row_idx,i), centered_rmse_norm(row_idx,i), &
                                                   bias_norm(row_idx,i) )
            cycle
          end if

          call calculate_taylor_metrics( model_profile, benchmark_profile, field_correlation, field_std_ratio, &
                                         field_centered_rmse_norm, field_bias_norm )
          if ( .not. is_finite_core_value( field_correlation ) .or. &
               .not. is_finite_core_value( field_std_ratio ) .or. &
               .not. is_finite_core_value( field_centered_rmse_norm ) .or. &
               .not. is_finite_core_value( field_bias_norm ) ) then
            call set_invalid_field_metric_outputs( scaled_rmse(row_idx,i), correlation(row_idx,i), &
                                                   std_ratio(row_idx,i), centered_rmse_norm(row_idx,i), &
                                                   bias_norm(row_idx,i) )
            cycle
          end if

          correlation(row_idx,i) = field_correlation
          std_ratio(row_idx,i) = field_std_ratio
          centered_rmse_norm(row_idx,i) = field_centered_rmse_norm
          bias_norm(row_idx,i) = field_bias_norm
        end do

        deallocate( model_profile, benchmark_profile )

      end associate
    end do

  end subroutine calculate_field_loss

end module clubb_loss_driver
