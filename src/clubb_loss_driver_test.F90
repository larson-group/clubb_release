!----------------------------------------------------------------------
! $Id$
! clubb_loss_driver_test.F90
!
! Native Fortran regression harness for the in-memory CLUBB loss driver.
!
! Usage:
!
!   Normally run through run_scm_loss.py:
!
!     python run_scripts/run_scm_loss.py -driver_test -cases bomex
!
!   run_scm_loss.py builds the loss-driver namelist and passes that generated
!   input file to this executable.  A direct run can also pass the namelist path
!   as the first command-line argument.  If no argument is provided, this test
!   looks for clubb.in.
!
!   Like the normal standalone drivers, exit code 6 means success.
!
! Understanding:
!
!   This test is about the loss-driver contract used by the tuner, not about
!   validating every physical result in a case.
!
!   The first check is deliberately small: calculate_profile_loss() is tested
!   against simple profiles where the expected correlation, standard-deviation ratio,
!   weighted normalized RMSE, and bias are obvious.
!
!   After that, the test runs the real loss driver in two forms:
!
!     - clubb_get_loss() is the one-shot path.  It initializes CLUBB, runs the
!       requested case, calculates the loss, and finalizes everything.
!
!     - init_clubb_loss(), clubb_get_loss_for_params(), and
!       finalize_clubb_loss() are the reusable path used by tuning.  The tuner
!       initializes once, repeatedly changes the parameter matrix, and asks for
!       new loss values without rebuilding all state from scratch.
!
!   The reusable path must reproduce the one-shot result BFB when the parameters
!   are unchanged.  Repeated reusable calls must also stay BFB, changing C8 must
!   change the loss, and finalizing/reinitializing must return to the same
!   baseline.  Those checks catch stale state, ignored parameter updates, and
!   leaked state between tuner evaluations.
!
!   The Python harness covers the same reusable-driver contract through the
!   f2py layer.  This file tests the native Fortran API directly.
!-----------------------------------------------------------------------

program clubb_loss_driver_test

  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan, ieee_is_finite

  use clubb_loss_driver, only: &
    loss_name_len, &
    clubb_get_loss, &
    init_clubb_loss, &
    clubb_get_loss_for_params, &
    finalize_clubb_loss, &
    calculate_profile_loss

  use err_info_type_module, only: &
    err_info_type, &
    cleanup_err_info_api

  use clubb_precision, only: &
    core_rknd

  use parameter_indices, only: &
    iC8

  implicit none

  integer, parameter :: &
    success_code = 6

  integer :: &
    arg_len, &
    arg_status

  character(len=256) :: &
    namelist_filename

  real( kind = core_rknd ), allocatable, dimension(:,:,:) :: &
      rmse, &
      correlation, &
      std_ratio, &
      bias_norm
  real( kind = core_rknd ), allocatable, dimension(:) :: &
      case_loss

  character(len=loss_name_len), allocatable, dimension(:) :: &
    clubb_var_names

  type(err_info_type) :: &
    err_info

  character(len=256) :: arg

  namelist_filename = "clubb.in"
  if ( command_argument_count() >= 1 ) then
    call get_command_argument( 1, arg, length = arg_len, status = arg_status )
    if ( arg_status == 0 .and. arg_len > 0 ) then
      namelist_filename = trim( arg(1:arg_len) )
    end if
  end if

  call run_taylor_metric_known_profile_test()

  ! Run the normal one-shot loss path first.
  call clubb_get_loss( trim( namelist_filename ), clubb_var_names, rmse, correlation, std_ratio, &
                       bias_norm, case_loss, err_info )

  ! Reuse the same request manually and verify that repeated calls stay BFB.
  call run_manual_loss_bfb_test( trim( namelist_filename ), clubb_var_names, rmse, correlation, &
                                 std_ratio, bias_norm, case_loss, err_info )
  call run_tweaked_parameter_test( trim( namelist_filename ), clubb_var_names, rmse, err_info )
  call run_reinitialize_loss_bfb_test( trim( namelist_filename ), clubb_var_names, rmse, correlation, &
                                      std_ratio, bias_norm, case_loss, err_info )

  call print_loss_matrix( clubb_var_names, rmse, correlation, std_ratio, bias_norm, case_loss )

  call cleanup_err_info_api( err_info )
  call exit( success_code )

contains

  subroutine assert_close( test_name, actual, expected )

    ! Description:
    !   Require one scalar diagnostic to match an expected value to roundoff.

    implicit none

    character(len=*), intent(in) :: &
      test_name

    real( kind = core_rknd ), intent(in) :: &
      actual, &
      expected

    real( kind = core_rknd ), parameter :: &
      tolerance = 1.0e-12_core_rknd

    if ( ieee_is_nan( expected ) ) then
      if ( .not. ieee_is_nan( actual ) ) stop "Expected an undefined diagnostic"
      return
    end if

    if ( .not. ieee_is_finite( actual ) ) stop "Expected a finite diagnostic"

    if ( abs( actual - expected ) > tolerance ) then
      write( *, '(a,1x,es16.8,1x,a,1x,es16.8)' ) &
        trim( test_name )//" failed:", actual, "expected", expected
      stop 1
    end if

  end subroutine assert_close

  subroutine run_taylor_metric_known_profile_test()

    ! Description:
    !   Check the four profile metrics, normalized total RMSE and fixed-weight
    !   behavior on exact, rescaled, flat and opposite-pattern profiles.

    implicit none

    !--------------------------- Local Variables ---------------------------
    real( kind = core_rknd ), dimension(3) :: &
      model_profile, &                                             ! Synthetic compared values.
      benchmark_profile, &                                         ! Paired reference values.
      weights                                                      ! Fixed mean/std/correlation weights.

    real( kind = core_rknd ) :: &
      rmse_value, &                                                ! Weighted normalized profile error.
      correlation_value, &                                         ! Centered profile correlation.
      std_ratio_value, &                                           ! Standard-deviation ratio.
      bias_value                                                   ! Normalized mean bias.

    !--------------------------- Begin Code ---------------------------
    weights = 1.0_core_rknd
    model_profile = [1.0_core_rknd, 2.0_core_rknd, 3.0_core_rknd]
    benchmark_profile = model_profile

    call calculate_profile_loss( model_profile, benchmark_profile, weights, & ! In
                                 rmse_value, correlation_value, std_ratio_value, bias_value ) ! Out
    call assert_close( "Identical RMSE", rmse_value, 0.0_core_rknd )
    call assert_close( "Identical correlation", correlation_value, 1.0_core_rknd )
    call assert_close( "Identical std ratio", std_ratio_value, 1.0_core_rknd )
    call assert_close( "Identical bias", bias_value, 0.0_core_rknd )

    model_profile = 2.0_core_rknd * benchmark_profile

    call calculate_profile_loss( model_profile, benchmark_profile, weights, & ! In
                                 rmse_value, correlation_value, std_ratio_value, bias_value ) ! Out
    call assert_close( "Rescaled total RMSE", rmse_value, sqrt( 7.0_core_rknd ) )
    call assert_close( "Rescaled correlation", correlation_value, 1.0_core_rknd )
    call assert_close( "Rescaled std ratio", std_ratio_value, 2.0_core_rknd )

    model_profile = 2.0_core_rknd

    call calculate_profile_loss( model_profile, benchmark_profile, weights, & ! In
                                 rmse_value, correlation_value, std_ratio_value, bias_value ) ! Out
    call assert_close( "Flat model RMSE", rmse_value, 1.0_core_rknd )
    call assert_close( "Flat model correlation", correlation_value, ieee_value( 0.0_core_rknd, ieee_quiet_nan ) )
    call assert_close( "Flat model std ratio", std_ratio_value, 0.0_core_rknd )

    benchmark_profile = model_profile

    call calculate_profile_loss( model_profile, benchmark_profile, weights, & ! In
                                 rmse_value, correlation_value, std_ratio_value, bias_value ) ! Out
    call assert_close( "Flat exact RMSE", rmse_value, 0.0_core_rknd )
    call assert_close( "Flat correlation", correlation_value, ieee_value( 0.0_core_rknd, ieee_quiet_nan ) )
    call assert_close( "Flat std ratio", std_ratio_value, ieee_value( 0.0_core_rknd, ieee_quiet_nan ) )
    call assert_close( "Flat bias", bias_value, 0.0_core_rknd )

    model_profile = [1.0_core_rknd, 2.0_core_rknd, 3.0_core_rknd]

    call calculate_profile_loss( model_profile, benchmark_profile, weights, & ! In
                                 rmse_value, correlation_value, std_ratio_value, bias_value ) ! Out
    call assert_close( "Flat reference variability error", rmse_value, sqrt( 2.0_core_rknd / 3.0_core_rknd ) )

    model_profile = 0.2_core_rknd
    benchmark_profile = 0.1_core_rknd

    call calculate_profile_loss( model_profile, benchmark_profile, weights, & ! In
                                 rmse_value, correlation_value, std_ratio_value, bias_value ) ! Out
    call assert_close( "Rounded constant RMSE", rmse_value, 0.1_core_rknd )
    call assert_close( "Rounded constant correlation", correlation_value, ieee_value( 0.0_core_rknd, ieee_quiet_nan ) )
    call assert_close( "Rounded constant std ratio", std_ratio_value, ieee_value( 0.0_core_rknd, ieee_quiet_nan ) )
    call assert_close( "Rounded constant bias", bias_value, 0.1_core_rknd )

    model_profile = [1.0_core_rknd, 2.0_core_rknd, 3.0_core_rknd]
    benchmark_profile = 2.0_core_rknd

    weights = [2.0_core_rknd, 3.0_core_rknd, 4.0_core_rknd]

    call calculate_profile_loss( model_profile, benchmark_profile, weights, & ! In
                                 rmse_value, correlation_value, std_ratio_value, bias_value ) ! Out
    call assert_close( "Flat reference weighted variability error", rmse_value, sqrt( 2.0_core_rknd ) )

    benchmark_profile = [3.0_core_rknd, 2.0_core_rknd, 1.0_core_rknd]
    weights = [2.0_core_rknd, 0.25_core_rknd, 3.0_core_rknd]

    call calculate_profile_loss( model_profile, benchmark_profile, weights, & ! In
                                 rmse_value, correlation_value, std_ratio_value, bias_value ) ! Out
    call assert_close( "Opposite pattern weighted RMSE", rmse_value, sqrt( 12.0_core_rknd ) )
    call assert_close( "Opposite pattern correlation", correlation_value, -1.0_core_rknd )

  end subroutine run_taylor_metric_known_profile_test

  subroutine assert_loss_outputs_match( test_name, test_var_names, test_rmse, test_correlation, &
                                        test_std_ratio, test_bias_norm, test_case_loss, &
                                        reference_var_names, reference_rmse, reference_correlation, &
                                        reference_std_ratio, reference_bias_norm, reference_case_loss )

    ! Description:
    !   Require two loss-driver outputs to match bit-for-bit.

    implicit none

    character(len=*), intent(in) :: &
      test_name

    character(len=*), dimension(:), intent(in) :: &
      test_var_names, &
      reference_var_names

    real( kind = core_rknd ), dimension(:,:,:), intent(in) :: &
      test_rmse, &
      test_correlation, &
      test_std_ratio, &
      test_bias_norm
    real( kind = core_rknd ), dimension(:), intent(in) :: &
      test_case_loss
    real( kind = core_rknd ), dimension(:,:,:), intent(in) :: &
      reference_rmse, &
      reference_correlation, &
      reference_std_ratio, &
      reference_bias_norm
    real( kind = core_rknd ), dimension(:), intent(in) :: &
      reference_case_loss

    if ( size( test_var_names ) /= size( reference_var_names ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: variable-name count changed"
      stop 1
    end if

    if ( any( test_var_names /= reference_var_names ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: variable-name ordering changed"
      stop 1
    end if

    if ( any( shape( test_rmse ) /= shape( reference_rmse ) ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: rmse shapes changed"
      stop 1
    end if

    if ( any( shape( test_correlation ) /= shape( reference_correlation ) ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: correlation shapes changed"
      stop 1
    end if

    if ( any( shape( test_std_ratio ) /= shape( reference_std_ratio ) ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: std_ratio shapes changed"
      stop 1
    end if

    if ( any( shape( test_case_loss ) /= shape( reference_case_loss ) ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: case_loss shapes changed"
      stop 1
    end if

    if ( any( shape( test_bias_norm ) /= shape( reference_bias_norm ) ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: bias_norm shapes changed"
      stop 1
    end if

    if ( any( test_rmse /= reference_rmse ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: rmse changed"
      stop 1
    end if

    if ( any( ieee_is_nan( test_correlation ) .neqv. ieee_is_nan( reference_correlation ) ) .or. &
         any( merge( 0.0_core_rknd, test_correlation, ieee_is_nan( test_correlation ) ) /= &
              merge( 0.0_core_rknd, reference_correlation, ieee_is_nan( reference_correlation ) ) ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: correlation changed"
      stop 1
    end if

    if ( any( ieee_is_nan( test_std_ratio ) .neqv. ieee_is_nan( reference_std_ratio ) ) .or. &
         any( merge( 0.0_core_rknd, test_std_ratio, ieee_is_nan( test_std_ratio ) ) /= &
              merge( 0.0_core_rknd, reference_std_ratio, ieee_is_nan( reference_std_ratio ) ) ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: std_ratio changed"
      stop 1
    end if

    if ( any( test_case_loss /= reference_case_loss ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: case_loss changed"
      stop 1
    end if

    if ( any( test_bias_norm /= reference_bias_norm ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: bias_norm changed"
      stop 1
    end if

  end subroutine assert_loss_outputs_match

  subroutine run_manual_loss_bfb_test( namelist_filename, reference_var_names, reference_rmse, &
                                       reference_correlation, reference_std_ratio, &
                                       reference_bias_norm, reference_case_loss, err_info )

    ! Description:
    !   Reuse the explicit init/run/finalize loss path twice with the same
    !   parameter matrix and require BFB agreement with the one-shot path.

    implicit none

    character(len=*), intent(in) :: &
      namelist_filename

    character(len=*), dimension(:), intent(in) :: &
      reference_var_names

    real( kind = core_rknd ), dimension(:,:,:), intent(in) :: &
      reference_rmse, &
      reference_correlation, &
      reference_std_ratio, &
      reference_bias_norm
    real( kind = core_rknd ), dimension(:), intent(in) :: &
      reference_case_loss

    type(err_info_type), intent(out) :: &
      err_info

    character(len=loss_name_len), allocatable, dimension(:) :: &
      manual_var_names

    real( kind = core_rknd ), allocatable, dimension(:,:) :: &
      clubb_params_all

    real( kind = core_rknd ), allocatable, dimension(:,:,:) :: &
      manual_rmse_1, &
      manual_correlation_1, &
      manual_std_ratio_1, &
      manual_bias_norm_1
    real( kind = core_rknd ), allocatable, dimension(:) :: &
      manual_case_loss_1
    real( kind = core_rknd ), allocatable, dimension(:,:,:) :: &
      manual_rmse_2, &
      manual_correlation_2, &
      manual_std_ratio_2, &
      manual_bias_norm_2
    real( kind = core_rknd ), allocatable, dimension(:) :: &
      manual_case_loss_2

    ! Initialize once and capture the full default parameter matrix.
    call init_clubb_loss( trim( namelist_filename ), manual_var_names, err_info, clubb_params_all )

    ! Score the same parameter matrix twice without reinitializing CLUBB.
    call clubb_get_loss_for_params( clubb_params_all, manual_rmse_1, manual_correlation_1, &
                                    manual_std_ratio_1, manual_bias_norm_1, manual_case_loss_1, err_info )
    call clubb_get_loss_for_params( clubb_params_all, manual_rmse_2, manual_correlation_2, &
                                    manual_std_ratio_2, manual_bias_norm_2, manual_case_loss_2, err_info )

    ! Release the manual loss session after both reruns are done.
    call finalize_clubb_loss( err_info )

    call assert_loss_outputs_match( "Manual loss BFB check (first rerun)", manual_var_names, &
                                    manual_rmse_1, manual_correlation_1, manual_std_ratio_1, &
                                    manual_bias_norm_1, manual_case_loss_1, reference_var_names, &
                                    reference_rmse, reference_correlation, reference_std_ratio, &
                                    reference_bias_norm, reference_case_loss )
    call assert_loss_outputs_match( "Manual loss BFB check (second rerun)", manual_var_names, &
                                    manual_rmse_2, manual_correlation_2, manual_std_ratio_2, &
                                    manual_bias_norm_2, manual_case_loss_2, reference_var_names, &
                                    reference_rmse, reference_correlation, reference_std_ratio, &
                                    reference_bias_norm, reference_case_loss )

    call assert_loss_outputs_match( "Manual repeated metrics", manual_var_names, & ! In
                                    manual_rmse_2, manual_correlation_2, manual_std_ratio_2, & ! In
                                    manual_bias_norm_2, manual_case_loss_2, manual_var_names, & ! In
                                    manual_rmse_1, manual_correlation_1, manual_std_ratio_1, & ! In
                                    manual_bias_norm_1, manual_case_loss_1 ) ! In

    write( *, '(a)' ) "Manual loss BFB check passed"

  end subroutine run_manual_loss_bfb_test

  subroutine run_tweaked_parameter_test( namelist_filename, reference_var_names, reference_rmse, err_info )

    ! Description:
    !   Verify that a deliberate parameter perturbation changes the loss.

    implicit none

    character(len=*), intent(in) :: &
      namelist_filename

    character(len=*), dimension(:), intent(in) :: &
      reference_var_names

    real( kind = core_rknd ), dimension(:,:,:), intent(in) :: &
      reference_rmse

    type(err_info_type), intent(out) :: &
      err_info

    character(len=loss_name_len), allocatable, dimension(:) :: &
      tweaked_var_names

    real( kind = core_rknd ), allocatable, dimension(:,:) :: &
      clubb_params_all

    real( kind = core_rknd ), allocatable, dimension(:,:,:) :: &
      tweaked_rmse, &
      tweaked_correlation, &
      tweaked_std_ratio, &
      tweaked_bias_norm
    real( kind = core_rknd ), allocatable, dimension(:) :: &
      tweaked_case_loss

    call init_clubb_loss( trim( namelist_filename ), tweaked_var_names, err_info, clubb_params_all )

    if ( size( tweaked_var_names ) /= size( reference_var_names ) ) then
      write( *, '(a)' ) "Tweaked-parameter test failed: variable-name count changed"
      stop 1
    end if

    if ( any( tweaked_var_names /= reference_var_names ) ) then
      write( *, '(a)' ) "Tweaked-parameter test failed: variable-name ordering changed"
      stop 1
    end if

    ! Perturb one tunable parameter enough that the resulting loss should move.
    clubb_params_all(1,iC8) = clubb_params_all(1,iC8) + 0.5_core_rknd

    call clubb_get_loss_for_params( clubb_params_all, tweaked_rmse, tweaked_correlation, &
                                    tweaked_std_ratio, tweaked_bias_norm, tweaked_case_loss, err_info )
    call finalize_clubb_loss( err_info )

    if ( all( tweaked_rmse == reference_rmse ) ) then
      write( *, '(a)' ) "Tweaked-parameter test failed: changing C8 did not change rmse"
      stop 1
    end if

    write( *, '(a)' ) "Tweaked-parameter test passed"

  end subroutine run_tweaked_parameter_test

  subroutine run_reinitialize_loss_bfb_test( namelist_filename, reference_var_names, reference_rmse, &
                                             reference_correlation, reference_std_ratio, &
                                             reference_bias_norm, reference_case_loss, err_info )

    ! Description:
    !   Reinitialize the reusable loss path after finalization and require that
    !   it reproduces the one-shot loss bit-for-bit.

    implicit none

    character(len=*), intent(in) :: &
      namelist_filename

    character(len=*), dimension(:), intent(in) :: &
      reference_var_names

    real( kind = core_rknd ), dimension(:,:,:), intent(in) :: &
      reference_rmse, &
      reference_correlation, &
      reference_std_ratio, &
      reference_bias_norm
    real( kind = core_rknd ), dimension(:), intent(in) :: &
      reference_case_loss

    type(err_info_type), intent(out) :: &
      err_info

    character(len=loss_name_len), allocatable, dimension(:) :: &
      reinit_var_names_1, &
      reinit_var_names_2

    real( kind = core_rknd ), allocatable, dimension(:,:) :: &
      clubb_params_all_1, &
      clubb_params_all_2

    real( kind = core_rknd ), allocatable, dimension(:,:,:) :: &
      reinit_rmse_1, &
      reinit_correlation_1, &
      reinit_std_ratio_1, &
      reinit_bias_norm_1
    real( kind = core_rknd ), allocatable, dimension(:) :: &
      reinit_case_loss_1
    real( kind = core_rknd ), allocatable, dimension(:,:,:) :: &
      reinit_rmse_2, &
      reinit_correlation_2, &
      reinit_std_ratio_2, &
      reinit_bias_norm_2
    real( kind = core_rknd ), allocatable, dimension(:) :: &
      reinit_case_loss_2

    call init_clubb_loss( trim( namelist_filename ), reinit_var_names_1, err_info, clubb_params_all_1 )
    call clubb_get_loss_for_params( clubb_params_all_1, reinit_rmse_1, reinit_correlation_1, &
                                    reinit_std_ratio_1, reinit_bias_norm_1, reinit_case_loss_1, err_info )
    call finalize_clubb_loss( err_info )

    call init_clubb_loss( trim( namelist_filename ), reinit_var_names_2, err_info, clubb_params_all_2 )
    call clubb_get_loss_for_params( clubb_params_all_2, reinit_rmse_2, reinit_correlation_2, &
                                    reinit_std_ratio_2, reinit_bias_norm_2, reinit_case_loss_2, err_info )
    call finalize_clubb_loss( err_info )

    call assert_loss_outputs_match( "Reinitialize loss BFB check (first reinit)", reinit_var_names_1, &
                                    reinit_rmse_1, reinit_correlation_1, reinit_std_ratio_1, &
                                    reinit_bias_norm_1, reinit_case_loss_1, reference_var_names, &
                                    reference_rmse, reference_correlation, reference_std_ratio, &
                                    reference_bias_norm, reference_case_loss )
    call assert_loss_outputs_match( "Reinitialize loss BFB check (second reinit)", reinit_var_names_2, &
                                    reinit_rmse_2, reinit_correlation_2, reinit_std_ratio_2, &
                                    reinit_bias_norm_2, reinit_case_loss_2, reference_var_names, &
                                    reference_rmse, reference_correlation, reference_std_ratio, &
                                    reference_bias_norm, reference_case_loss )

    if ( any( reinit_var_names_1 /= reinit_var_names_2 ) ) then
      write( *, '(a)' ) "Reinitialize loss BFB check failed: variable-name ordering changed across reinit"
      stop 1
    end if

    if ( any( clubb_params_all_1 /= clubb_params_all_2 ) ) then
      write( *, '(a)' ) "Reinitialize loss BFB check failed: default parameter matrix changed across reinit"
      stop 1
    end if

    call assert_loss_outputs_match( "Reinitialized repeated metrics", reinit_var_names_1, & ! In
                                    reinit_rmse_1, reinit_correlation_1, reinit_std_ratio_1, & ! In
                                    reinit_bias_norm_1, reinit_case_loss_1, reinit_var_names_2, & ! In
                                    reinit_rmse_2, reinit_correlation_2, reinit_std_ratio_2, & ! In
                                    reinit_bias_norm_2, reinit_case_loss_2 ) ! In

    write( *, '(a)' ) "Reinitialize loss BFB check passed"

  end subroutine run_reinitialize_loss_bfb_test

  subroutine print_loss_matrix( clubb_var_names, rmse, correlation, std_ratio, bias_norm, case_loss )

    ! Description:
  !   Prints one row per variable/window together with the explicit
  !   loss-driver metrics and the best-performing parameter column.

    implicit none

    character(len=*), dimension(:), intent(in) :: &
      clubb_var_names ! CLUBB variable names, one per loss row.

    real( kind = core_rknd ), dimension(:,:,:), intent(in) :: &
      rmse, & ! RMSE shaped (window, variable, parameter column).
      correlation, & ! Taylor correlation shaped (window, variable, parameter column).
      std_ratio ! Taylor standard-deviation ratio shaped (window, variable, parameter column).
    real( kind = core_rknd ), dimension(:), intent(in) :: &
      case_loss ! Taylor centered RMSE shaped (window, variable, parameter column).
    real( kind = core_rknd ), dimension(:,:,:), intent(in) :: &
      bias_norm ! Taylor normalized bias shaped (window, variable, parameter column).

    integer :: &
      row_idx, &      ! Current loss-variable row being printed.
      col_idx, &      ! Current parameter column being printed.
      window_idx, &   ! Current time-window row being printed.
      best_col_idx    ! Parameter column with the minimum loss for this row.

    real( kind = core_rknd ) :: &
      best_col_loss   ! Lowest loss found across all parameter columns.

    write( *, '(a)', advance = "no" ) "variable"
    if ( size( rmse, 1 ) > 1 ) write( *, '(1x,a)', advance = "no" ) "window"
    do col_idx = 1, size( rmse, 3 )
      write( *, '(1x,a,i0)', advance = "no" ) "rmse_col", col_idx
    end do
    do col_idx = 1, size( rmse, 3 )
      write( *, '(1x,a,i0)' , advance = "no" ) "corr_col", col_idx
    end do
    do col_idx = 1, size( rmse, 3 )
      write( *, '(1x,a,i0)' , advance = "no" ) "std_ratio_col", col_idx
    end do
    do col_idx = 1, size( rmse, 3 )
      write( *, '(1x,a,i0)' , advance = "no" ) "bias_norm_col", col_idx
    end do
    write( *, '(1x,a)' ) "best_col"

    do window_idx = 1, size( rmse, 1 )
      do row_idx = 1, size( rmse, 2 )
        ! Summarize the best-performing parameter column directly in the printed
        ! table so manual tuning runs are easier to scan.
        best_col_idx = 1
        best_col_loss = rmse(window_idx,row_idx,1)
        do col_idx = 2, size( rmse, 3 )
          if ( rmse(window_idx,row_idx,col_idx) < best_col_loss ) then
            best_col_idx = col_idx
            best_col_loss = rmse(window_idx,row_idx,col_idx)
          end if
        end do

        write( *, '(a)', advance = "no" ) trim( clubb_var_names(row_idx) )
        if ( size( rmse, 1 ) > 1 ) write( *, '(1x,a,i0)', advance = "no" ) "window", window_idx
        do col_idx = 1, size( rmse, 3 )
          write( *, '(1x,es16.8)', advance = "no" ) rmse(window_idx,row_idx,col_idx)
        end do
        do col_idx = 1, size( rmse, 3 )
          write( *, '(1x,es16.8)', advance = "no" ) correlation(window_idx,row_idx,col_idx)
        end do
        do col_idx = 1, size( rmse, 3 )
          write( *, '(1x,es16.8)', advance = "no" ) std_ratio(window_idx,row_idx,col_idx)
        end do
        do col_idx = 1, size( rmse, 3 )
          write( *, '(1x,es16.8)', advance = "no" ) bias_norm(window_idx,row_idx,col_idx)
        end do
        write( *, '(1x,a,i0)' ) "col", best_col_idx
      end do
    end do

    write( *, '(a)', advance = "no" ) "case_loss"
    do col_idx = 1, size( case_loss )
      write( *, '(1x,es16.8)', advance = "no" ) case_loss(col_idx)
    end do
    write( *, * )

  end subroutine print_loss_matrix
end program clubb_loss_driver_test
