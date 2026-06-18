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
!   The first check is deliberately small: calculate_taylor_metrics() is tested
!   against simple profiles where the expected correlation, variance ratio,
!   centered RMSE, and bias are obvious.
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

  use clubb_loss_driver, only: &
    loss_name_len, &
    clubb_get_loss, &
    init_clubb_loss, &
    clubb_get_loss_for_params, &
    finalize_clubb_loss, &
    calculate_taylor_metrics

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
    scaled_rmse, &
    correlation, &
    std_ratio, &
    centered_rmse_norm, &
    bias_norm

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
  call clubb_get_loss( trim( namelist_filename ), clubb_var_names, scaled_rmse, correlation, std_ratio, &
                       centered_rmse_norm, bias_norm, err_info )

  ! Reuse the same request manually and verify that repeated calls stay BFB.
  call run_manual_loss_bfb_test( trim( namelist_filename ), clubb_var_names, scaled_rmse, correlation, &
                                 std_ratio, centered_rmse_norm, bias_norm, err_info )
  call run_tweaked_parameter_test( trim( namelist_filename ), clubb_var_names, scaled_rmse, err_info )
  call run_reinitialize_loss_bfb_test( trim( namelist_filename ), clubb_var_names, scaled_rmse, correlation, &
                                      std_ratio, centered_rmse_norm, bias_norm, err_info )

  call print_loss_matrix( clubb_var_names, scaled_rmse, correlation, std_ratio, centered_rmse_norm, bias_norm )

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

    if ( abs( actual - expected ) > tolerance ) then
      write( *, '(a,1x,es16.8,1x,a,1x,es16.8)' ) &
        trim( test_name )//" failed:", actual, "expected", expected
      stop 1
    end if

  end subroutine assert_close

  subroutine run_taylor_metric_known_profile_test()

    ! Description:
    !   Check Taylor diagnostic arithmetic against simple profile cases with
    !   obvious expected values.

    implicit none

    real( kind = core_rknd ), dimension(3) :: &
      model_profile, &
      benchmark_profile

    real( kind = core_rknd ) :: &
      correlation, &
      std_ratio, &
      centered_rmse_norm, &
      bias_norm

    model_profile = [ 1.0_core_rknd, 2.0_core_rknd, 3.0_core_rknd ]
    benchmark_profile = model_profile

    call calculate_taylor_metrics( model_profile, benchmark_profile, correlation, std_ratio, &
                                   centered_rmse_norm, bias_norm )
    call assert_close( "Taylor identical correlation", correlation, 1.0_core_rknd )
    call assert_close( "Taylor identical std_ratio", std_ratio, 1.0_core_rknd )
    call assert_close( "Taylor identical centered_rmse_norm", centered_rmse_norm, 0.0_core_rknd )
    call assert_close( "Taylor identical bias_norm", bias_norm, 0.0_core_rknd )

    model_profile = [ 2.0_core_rknd, 4.0_core_rknd, 6.0_core_rknd ]
    benchmark_profile = [ 1.0_core_rknd, 2.0_core_rknd, 3.0_core_rknd ]

    call calculate_taylor_metrics( model_profile, benchmark_profile, correlation, std_ratio, &
                                   centered_rmse_norm, bias_norm )
    call assert_close( "Taylor scaled correlation", correlation, 1.0_core_rknd )
    call assert_close( "Taylor scaled std_ratio", std_ratio, 2.0_core_rknd )
    call assert_close( "Taylor scaled centered_rmse_norm", centered_rmse_norm, 1.0_core_rknd )

    model_profile = [ 2.0_core_rknd, 2.0_core_rknd, 2.0_core_rknd ]
    benchmark_profile = [ 2.0_core_rknd, 2.0_core_rknd, 2.0_core_rknd ]

    call calculate_taylor_metrics( model_profile, benchmark_profile, correlation, std_ratio, &
                                   centered_rmse_norm, bias_norm )
    call assert_close( "Taylor flat correlation", correlation, 0.0_core_rknd )
    call assert_close( "Taylor flat std_ratio", std_ratio, 0.0_core_rknd )
    call assert_close( "Taylor flat centered_rmse_norm", centered_rmse_norm, 0.0_core_rknd )
    call assert_close( "Taylor flat bias_norm", bias_norm, 0.0_core_rknd )

  end subroutine run_taylor_metric_known_profile_test

  subroutine assert_loss_outputs_match( test_name, test_var_names, test_scaled_rmse, test_correlation, &
                                        test_std_ratio, test_centered_rmse_norm, test_bias_norm, &
                                        reference_var_names, reference_scaled_rmse, reference_correlation, &
                                        reference_std_ratio, reference_centered_rmse_norm, reference_bias_norm )

    ! Description:
    !   Require two loss-driver outputs to match bit-for-bit.

    implicit none

    character(len=*), intent(in) :: &
      test_name

    character(len=*), dimension(:), intent(in) :: &
      test_var_names, &
      reference_var_names

    real( kind = core_rknd ), dimension(:,:,:), intent(in) :: &
      test_scaled_rmse, &
      test_correlation, &
      test_std_ratio, &
      test_centered_rmse_norm, &
      test_bias_norm, &
      reference_scaled_rmse, &
      reference_correlation, &
      reference_std_ratio, &
      reference_centered_rmse_norm, &
      reference_bias_norm

    if ( size( test_var_names ) /= size( reference_var_names ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: variable-name count changed"
      stop 1
    end if

    if ( any( test_var_names /= reference_var_names ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: variable-name ordering changed"
      stop 1
    end if

    if ( any( shape( test_scaled_rmse ) /= shape( reference_scaled_rmse ) ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: scaled_rmse shapes changed"
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

    if ( any( shape( test_centered_rmse_norm ) /= shape( reference_centered_rmse_norm ) ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: centered_rmse_norm shapes changed"
      stop 1
    end if

    if ( any( shape( test_bias_norm ) /= shape( reference_bias_norm ) ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: bias_norm shapes changed"
      stop 1
    end if

    if ( any( test_scaled_rmse /= reference_scaled_rmse ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: scaled_rmse changed"
      stop 1
    end if

    if ( any( test_correlation /= reference_correlation ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: correlation changed"
      stop 1
    end if

    if ( any( test_std_ratio /= reference_std_ratio ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: std_ratio changed"
      stop 1
    end if

    if ( any( test_centered_rmse_norm /= reference_centered_rmse_norm ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: centered_rmse_norm changed"
      stop 1
    end if

    if ( any( test_bias_norm /= reference_bias_norm ) ) then
      write( *, '(a)' ) trim( test_name )//" failed: bias_norm changed"
      stop 1
    end if

  end subroutine assert_loss_outputs_match

  subroutine run_manual_loss_bfb_test( namelist_filename, reference_var_names, reference_scaled_rmse, &
                                       reference_correlation, reference_std_ratio, &
                                       reference_centered_rmse_norm, reference_bias_norm, err_info )

    ! Description:
    !   Reuse the explicit init/run/finalize loss path twice with the same
    !   parameter matrix and require BFB agreement with the one-shot path.

    implicit none

    character(len=*), intent(in) :: &
      namelist_filename

    character(len=*), dimension(:), intent(in) :: &
      reference_var_names

    real( kind = core_rknd ), dimension(:,:,:), intent(in) :: &
      reference_scaled_rmse, &
      reference_correlation, &
      reference_std_ratio, &
      reference_centered_rmse_norm, &
      reference_bias_norm

    type(err_info_type), intent(out) :: &
      err_info

    character(len=loss_name_len), allocatable, dimension(:) :: &
      manual_var_names

    real( kind = core_rknd ), allocatable, dimension(:,:) :: &
      clubb_params_all

    real( kind = core_rknd ), allocatable, dimension(:,:,:) :: &
      manual_scaled_rmse_1, &
      manual_correlation_1, &
      manual_std_ratio_1, &
      manual_centered_rmse_norm_1, &
      manual_bias_norm_1, &
      manual_scaled_rmse_2, &
      manual_correlation_2, &
      manual_std_ratio_2, &
      manual_centered_rmse_norm_2, &
      manual_bias_norm_2

    ! Initialize once and capture the full default parameter matrix.
    call init_clubb_loss( trim( namelist_filename ), manual_var_names, err_info, clubb_params_all )

    ! Score the same parameter matrix twice without reinitializing CLUBB.
    call clubb_get_loss_for_params( clubb_params_all, manual_scaled_rmse_1, manual_correlation_1, &
                                    manual_std_ratio_1, manual_centered_rmse_norm_1, manual_bias_norm_1, err_info )
    call clubb_get_loss_for_params( clubb_params_all, manual_scaled_rmse_2, manual_correlation_2, &
                                    manual_std_ratio_2, manual_centered_rmse_norm_2, manual_bias_norm_2, err_info )

    ! Release the manual loss session after both reruns are done.
    call finalize_clubb_loss( err_info )

    call assert_loss_outputs_match( "Manual loss BFB check (first rerun)", manual_var_names, &
                                    manual_scaled_rmse_1, manual_correlation_1, manual_std_ratio_1, &
                                    manual_centered_rmse_norm_1, manual_bias_norm_1, reference_var_names, &
                                    reference_scaled_rmse, reference_correlation, reference_std_ratio, &
                                    reference_centered_rmse_norm, reference_bias_norm )
    call assert_loss_outputs_match( "Manual loss BFB check (second rerun)", manual_var_names, &
                                    manual_scaled_rmse_2, manual_correlation_2, manual_std_ratio_2, &
                                    manual_centered_rmse_norm_2, manual_bias_norm_2, reference_var_names, &
                                    reference_scaled_rmse, reference_correlation, reference_std_ratio, &
                                    reference_centered_rmse_norm, reference_bias_norm )

    if ( any( manual_scaled_rmse_2 /= manual_scaled_rmse_1 ) .or. &
         any( manual_correlation_2 /= manual_correlation_1 ) .or. &
         any( manual_std_ratio_2 /= manual_std_ratio_1 ) .or. &
         any( manual_centered_rmse_norm_2 /= manual_centered_rmse_norm_1 ) .or. &
         any( manual_bias_norm_2 /= manual_bias_norm_1 ) ) then
      write( *, '(a)' ) "Manual loss BFB check failed: repeated manual metrics are not identical"
      stop 1
    end if

    write( *, '(a)' ) "Manual loss BFB check passed"

  end subroutine run_manual_loss_bfb_test

  subroutine run_tweaked_parameter_test( namelist_filename, reference_var_names, reference_scaled_rmse, err_info )

    ! Description:
    !   Verify that a deliberate parameter perturbation changes the loss.

    implicit none

    character(len=*), intent(in) :: &
      namelist_filename

    character(len=*), dimension(:), intent(in) :: &
      reference_var_names

    real( kind = core_rknd ), dimension(:,:,:), intent(in) :: &
      reference_scaled_rmse

    type(err_info_type), intent(out) :: &
      err_info

    character(len=loss_name_len), allocatable, dimension(:) :: &
      tweaked_var_names

    real( kind = core_rknd ), allocatable, dimension(:,:) :: &
      clubb_params_all

    real( kind = core_rknd ), allocatable, dimension(:,:,:) :: &
      tweaked_scaled_rmse, &
      tweaked_correlation, &
      tweaked_std_ratio, &
      tweaked_centered_rmse_norm, &
      tweaked_bias_norm

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

    call clubb_get_loss_for_params( clubb_params_all, tweaked_scaled_rmse, tweaked_correlation, &
                                    tweaked_std_ratio, tweaked_centered_rmse_norm, tweaked_bias_norm, err_info )
    call finalize_clubb_loss( err_info )

    if ( all( tweaked_scaled_rmse == reference_scaled_rmse ) ) then
      write( *, '(a)' ) "Tweaked-parameter test failed: changing C8 did not change scaled_rmse"
      stop 1
    end if

    write( *, '(a)' ) "Tweaked-parameter test passed"

  end subroutine run_tweaked_parameter_test

  subroutine run_reinitialize_loss_bfb_test( namelist_filename, reference_var_names, reference_scaled_rmse, &
                                             reference_correlation, reference_std_ratio, &
                                             reference_centered_rmse_norm, reference_bias_norm, err_info )

    ! Description:
    !   Reinitialize the reusable loss path after finalization and require that
    !   it reproduces the one-shot loss bit-for-bit.

    implicit none

    character(len=*), intent(in) :: &
      namelist_filename

    character(len=*), dimension(:), intent(in) :: &
      reference_var_names

    real( kind = core_rknd ), dimension(:,:,:), intent(in) :: &
      reference_scaled_rmse, &
      reference_correlation, &
      reference_std_ratio, &
      reference_centered_rmse_norm, &
      reference_bias_norm

    type(err_info_type), intent(out) :: &
      err_info

    character(len=loss_name_len), allocatable, dimension(:) :: &
      reinit_var_names_1, &
      reinit_var_names_2

    real( kind = core_rknd ), allocatable, dimension(:,:) :: &
      clubb_params_all_1, &
      clubb_params_all_2

    real( kind = core_rknd ), allocatable, dimension(:,:,:) :: &
      reinit_scaled_rmse_1, &
      reinit_correlation_1, &
      reinit_std_ratio_1, &
      reinit_centered_rmse_norm_1, &
      reinit_bias_norm_1, &
      reinit_scaled_rmse_2, &
      reinit_correlation_2, &
      reinit_std_ratio_2, &
      reinit_centered_rmse_norm_2, &
      reinit_bias_norm_2

    call init_clubb_loss( trim( namelist_filename ), reinit_var_names_1, err_info, clubb_params_all_1 )
    call clubb_get_loss_for_params( clubb_params_all_1, reinit_scaled_rmse_1, reinit_correlation_1, &
                                    reinit_std_ratio_1, reinit_centered_rmse_norm_1, reinit_bias_norm_1, err_info )
    call finalize_clubb_loss( err_info )

    call init_clubb_loss( trim( namelist_filename ), reinit_var_names_2, err_info, clubb_params_all_2 )
    call clubb_get_loss_for_params( clubb_params_all_2, reinit_scaled_rmse_2, reinit_correlation_2, &
                                    reinit_std_ratio_2, reinit_centered_rmse_norm_2, reinit_bias_norm_2, err_info )
    call finalize_clubb_loss( err_info )

    call assert_loss_outputs_match( "Reinitialize loss BFB check (first reinit)", reinit_var_names_1, &
                                    reinit_scaled_rmse_1, reinit_correlation_1, reinit_std_ratio_1, &
                                    reinit_centered_rmse_norm_1, reinit_bias_norm_1, reference_var_names, &
                                    reference_scaled_rmse, reference_correlation, reference_std_ratio, &
                                    reference_centered_rmse_norm, reference_bias_norm )
    call assert_loss_outputs_match( "Reinitialize loss BFB check (second reinit)", reinit_var_names_2, &
                                    reinit_scaled_rmse_2, reinit_correlation_2, reinit_std_ratio_2, &
                                    reinit_centered_rmse_norm_2, reinit_bias_norm_2, reference_var_names, &
                                    reference_scaled_rmse, reference_correlation, reference_std_ratio, &
                                    reference_centered_rmse_norm, reference_bias_norm )

    if ( any( reinit_var_names_1 /= reinit_var_names_2 ) ) then
      write( *, '(a)' ) "Reinitialize loss BFB check failed: variable-name ordering changed across reinit"
      stop 1
    end if

    if ( any( clubb_params_all_1 /= clubb_params_all_2 ) ) then
      write( *, '(a)' ) "Reinitialize loss BFB check failed: default parameter matrix changed across reinit"
      stop 1
    end if

    if ( any( reinit_scaled_rmse_1 /= reinit_scaled_rmse_2 ) .or. &
         any( reinit_correlation_1 /= reinit_correlation_2 ) .or. &
         any( reinit_std_ratio_1 /= reinit_std_ratio_2 ) .or. &
         any( reinit_centered_rmse_norm_1 /= reinit_centered_rmse_norm_2 ) .or. &
         any( reinit_bias_norm_1 /= reinit_bias_norm_2 ) ) then
      write( *, '(a)' ) "Reinitialize loss BFB check failed: repeated reinitialization changed metrics"
      stop 1
    end if

    write( *, '(a)' ) "Reinitialize loss BFB check passed"

  end subroutine run_reinitialize_loss_bfb_test

  subroutine print_loss_matrix( clubb_var_names, scaled_rmse, correlation, std_ratio, centered_rmse_norm, bias_norm )

    ! Description:
  !   Prints one row per variable/window together with the explicit
  !   loss-driver metrics and the best-performing parameter column.

    implicit none

    character(len=*), dimension(:), intent(in) :: &
      clubb_var_names ! CLUBB variable names, one per loss row.

    real( kind = core_rknd ), dimension(:,:,:), intent(in) :: &
      scaled_rmse, &        ! Scaled RMSE shaped (window, variable, parameter column).
      correlation, &        ! Taylor correlation shaped (window, variable, parameter column).
      std_ratio, &          ! Taylor standard-deviation ratio shaped (window, variable, parameter column).
      centered_rmse_norm, & ! Taylor centered RMSE shaped (window, variable, parameter column).
      bias_norm             ! Taylor normalized bias shaped (window, variable, parameter column).

    integer :: &
      row_idx, &      ! Current loss-variable row being printed.
      col_idx, &      ! Current parameter column being printed.
      window_idx, &   ! Current time-window row being printed.
      best_col_idx    ! Parameter column with the minimum loss for this row.

    real( kind = core_rknd ) :: &
      best_col_loss   ! Lowest loss found across all parameter columns.

    write( *, '(a)', advance = "no" ) "variable"
    if ( size( scaled_rmse, 1 ) > 1 ) write( *, '(1x,a)', advance = "no" ) "window"
    do col_idx = 1, size( scaled_rmse, 3 )
      write( *, '(1x,a,i0)', advance = "no" ) "scaled_rmse_col", col_idx
    end do
    do col_idx = 1, size( scaled_rmse, 3 )
      write( *, '(1x,a,i0)' , advance = "no" ) "corr_col", col_idx
    end do
    do col_idx = 1, size( scaled_rmse, 3 )
      write( *, '(1x,a,i0)' , advance = "no" ) "std_ratio_col", col_idx
    end do
    do col_idx = 1, size( scaled_rmse, 3 )
      write( *, '(1x,a,i0)' , advance = "no" ) "crmse_norm_col", col_idx
    end do
    do col_idx = 1, size( scaled_rmse, 3 )
      write( *, '(1x,a,i0)' , advance = "no" ) "bias_norm_col", col_idx
    end do
    write( *, '(1x,a)' ) "best_col"

    do window_idx = 1, size( scaled_rmse, 1 )
      do row_idx = 1, size( scaled_rmse, 2 )
        ! Summarize the best-performing parameter column directly in the printed
        ! table so manual tuning runs are easier to scan.
        best_col_idx = 1
        best_col_loss = scaled_rmse(window_idx,row_idx,1)
        do col_idx = 2, size( scaled_rmse, 3 )
          if ( scaled_rmse(window_idx,row_idx,col_idx) < best_col_loss ) then
            best_col_idx = col_idx
            best_col_loss = scaled_rmse(window_idx,row_idx,col_idx)
          end if
        end do

        write( *, '(a)', advance = "no" ) trim( clubb_var_names(row_idx) )
        if ( size( scaled_rmse, 1 ) > 1 ) write( *, '(1x,a,i0)', advance = "no" ) "window", window_idx
        do col_idx = 1, size( scaled_rmse, 3 )
          write( *, '(1x,es16.8)', advance = "no" ) scaled_rmse(window_idx,row_idx,col_idx)
        end do
        do col_idx = 1, size( scaled_rmse, 3 )
          write( *, '(1x,es16.8)', advance = "no" ) correlation(window_idx,row_idx,col_idx)
        end do
        do col_idx = 1, size( scaled_rmse, 3 )
          write( *, '(1x,es16.8)', advance = "no" ) std_ratio(window_idx,row_idx,col_idx)
        end do
        do col_idx = 1, size( scaled_rmse, 3 )
          write( *, '(1x,es16.8)', advance = "no" ) centered_rmse_norm(window_idx,row_idx,col_idx)
        end do
        do col_idx = 1, size( scaled_rmse, 3 )
          write( *, '(1x,es16.8)', advance = "no" ) bias_norm(window_idx,row_idx,col_idx)
        end do
        write( *, '(1x,a,i0)' ) "col", best_col_idx
      end do
    end do

  end subroutine print_loss_matrix
end program clubb_loss_driver_test
