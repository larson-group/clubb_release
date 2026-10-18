!-------------------------------------------------------------------------------
! Run complete CLUBB cases from aggregate namelists.
!
! Usage:
!   clubb_standalone [namelist_filename] runs once (default: clubb.in).
!   clubb_standalone -worker reads one namelist filename per stdin line until EOF.
!   run_case() calls the driver's full lifecycle and checks its error container;
!   worker_main() repeats that operation without restarting the executable.
!
! Understanding:
!   run_clubb owns initialization, all runtime batches, advancement and cleanup.
!   Model allocations are released between requests; only the process persists.
!   Workers flush CLUBB_WORKER_READY once, then CLUBB_WORKER_DONE after each
!   successful case and its output. Fatal errors exit without a completion marker.
!   EOF ends an idle worker. Native successful shutdown uses exit status 6;
!   utilities/clubb_case_runner.py recognizes it alongside Python's status 0.
!-------------------------------------------------------------------------------
program clubb_standalone

  ! Description:
  !   Run one namelist through initialization, advancement and cleanup.
  !   With -worker, read one namelist filename per line until EOF. Each case
  !   completes the same lifecycle; only the process remains alive between runs.

  use clubb_driver, only: &
    run_clubb

  use iso_fortran_env, only: &
    input_unit, output_unit

  use error_code, only: &
    clubb_generalized_grd_test_err, clubb_fatal_error

  use model_flags, only: &
    l_test_grid_generalization

  use constants_clubb, only: &
    fstderr

  use err_info_type_module, only: &
    err_info_type, cleanup_err_info_api

  implicit none

  !--------------------------- Local Variables ---------------------------
  integer, parameter :: &
    success_code = 6                       ! Native successful exit status [-]

  integer :: &
    arg_len, &                            ! Command argument length [characters]
    arg_status                            ! Command argument status [-]

  character( len = 1024 ) :: &
    arg                                   ! Namelist filename or -worker [-]

  !--------------------------- Begin Code ---------------------------

  ! Select one namelist to run, or -worker to read namelist paths from stdin.
  arg = 'clubb.in'

  if ( command_argument_count() >= 1 ) then

    call get_command_argument( 1, arg, length = arg_len, status = arg_status )

    if ( arg_status /= 0 .or. arg_len == 0 ) arg = 'clubb.in'

  end if

  if ( trim( arg ) == '-worker' ) then

    ! Run each incoming namelist in the same process until stdin reaches EOF.
    call worker_main()

  else

    ! Run one case and exit after its cleanup.
    call run_case( trim( arg ) ) ! In

  end if

  call exit( success_code )

contains

!-------------------------------------------------------------------------------
  subroutine worker_main()

    ! Description:
    !   Read aggregate namelist paths from input_unit until EOF, running each
    !   case before reading the next. Flush readiness/completion on output_unit;
    !   completion follows closed stats files and any enabled timing output.

    implicit none

    !--------------------------- Local Variables ---------------------------
    character( len = 1024 ) :: &
      namelist_filename                   ! Aggregate input path [-]

    integer :: &
      read_status                         ! Input status; negative at EOF [-]

    !--------------------------- Begin Code ---------------------------
    ! Let the parent finish starting all peers before it sends the first case.
    write( output_unit, '(a)' ) 'CLUBB_WORKER_READY'
    flush( output_unit )

    do

      read( input_unit, '(a)', iostat = read_status ) namelist_filename

      if ( read_status < 0 ) exit

      if ( read_status > 0 ) error stop 'Unable to read worker namelist filename'

      call run_case( trim( namelist_filename ) ) ! In

      write( output_unit, '(a)' ) 'CLUBB_WORKER_DONE'
      flush( output_unit )

    end do

  end subroutine worker_main

!-------------------------------------------------------------------------------
  subroutine run_case( namelist_filename )

    ! Description:
    !   Run namelist_filename through run_clubb, which writes configured output
    !   and releases model allocations. Stop on fatal/grid errors; on success,
    !   release the caller-owned error container before another case starts.

    implicit none

    !--------------------------- Input Variables ---------------------------
    character( len = * ), intent(in) :: &
      namelist_filename                   ! Aggregate input path [-]

    !--------------------------- Local Variables ---------------------------
    type(err_info_type) :: &
      err_info                            ! Driver-initialized error container [-]

    !--------------------------- Begin Code ---------------------------
    ! Keep allocation, runtime batching and output lifetime inside the driver.
    call run_clubb( namelist_filename, .true., & ! In
                    err_info )                 ! Out

    ! Preserve the standalone's existing fatal and generalized-grid diagnostics.
    if ( .not. l_test_grid_generalization ) then

      if ( any( err_info%err_code == clubb_fatal_error ) ) then
        error stop 'Fatal error in clubb, check your parameter values and timestep'
      end if

      write( fstderr, * ) 'Program exited normally'

    else

      if ( any( err_info%err_code == clubb_generalized_grd_test_err ) ) then
        error stop 'Error in generalized grid test; check the error messages'
      end if

      write( fstderr, * ) 'Generalized grid test passed'

    end if

    call cleanup_err_info_api( err_info ) ! InOut

  end subroutine run_case

end program clubb_standalone
