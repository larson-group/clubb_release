!----------------------------------------------------------------------
! $Id$

program clubb_driver_test

! Description:
! This is essentially a minimalist frontend for CLUBB.
!
!----------------------------------------------------------------------

  use clubb_driver, only: &
    init_clubb_case, &
    set_case_initial_conditions, &
    advance_clubb_to_end, &
    clean_up_clubb

  use error_code, only: &
    clubb_generalized_grd_test_err, & ! Constants
    clubb_fatal_error

  use model_flags, only: &
    l_test_grid_generalization

  use constants_clubb, only: &
    fstderr ! Constant

  use err_info_type_module, only: &
    err_info_type,                  & ! Type
    cleanup_err_info_api              ! Procedure(s)
    
  implicit none

  ! External
  intrinsic :: trim

  ! Constant parameters
  integer, parameter :: &
    success_code = 6

  integer :: &
    arg_len, &
    arg_status

  type(err_info_type) :: &
    err_info        ! err_info struct containing err_code and err_header
                    ! Initialization is done within run_clubb

  character(len=256) :: &
    namelist_filename  ! Text file containing namelists
                       ! concatenated from various input files such as
                       ! model.in, tunable_parameters.in, error....in.
  logical, parameter :: &
    l_stdout = .true.

  character(len=256) :: arg

  !--------------------------------- Begin Code ---------------------------------

  ! Set default namelist values.
  namelist_filename = "clubb.in"
  if ( command_argument_count() >= 1 ) then
    call get_command_argument(1, arg, length=arg_len, status=arg_status)
    if ( arg_status == 0 .and. arg_len > 0 ) then
      namelist_filename = trim(arg(1:arg_len))
    end if
  end if

  !---------------------------------------- Test Sections ----------------------------------------

  write(fstderr, *) "======================== REINITIALIZATION TEST ========================"
  write(fstderr, *) "This section ensures that everything allocated in init_clubb_case will be deallocated"
  write(fstderr, *) "in clean_up_clubb. This may cause a runtime error if there is a mismatch between"
  write(fstderr, *) "allocate/deallocate statements, but could "

  ! Read in model parameter values by initializing, cleaning up, and initializing again.
  call init_clubb_case(trim(namelist_filename), err_info)
  call clean_up_clubb(err_info)
  call init_clubb_case(trim(namelist_filename), err_info)


  write(fstderr, *) "======================== DOUBLE TIMESTEP RUN ========================"
  write(fstderr, *) "Calling advance_clubb_to_end, then set_case_initial_conditions, "
  write(fstderr, *) "then advance_clubb_to_end again. This could result in a runtime"
  write(fstderr, *) "error if an allocation statement appears in these routines, since"
  write(fstderr, *) "we don't call clean_up_clubb in between the advance_ calls."

  ! Run once with stats off. Turning stats on writes output to disk, and we only want
  ! that for the second run used for the standalone comparison.
  call advance_clubb_to_end( l_stdout, err_info, l_suppress_stats=.true. )

  ! Reset model back to initial conditions
  call set_case_initial_conditions(err_info)

  ! Turn stats on and run. If this produces different output than what calling run_clubb
  ! does, then set_case_initial_conditions is not correctly reseting to initial conditions
  call advance_clubb_to_end( l_stdout, err_info )

  write(fstderr, *) "======================== RUN OVER ========================"
  write(fstderr, *) "WARNING: The double timestep test is not complete until the "
  write(fstderr, *) "output is compared with the output from calling clubb_standalone."
  write(fstderr, *) "This driver should produce BFB output with clubb_standalone."
  write(fstderr, *) "Save the output, rerun with clubb_standalone, then compare."
  write(fstderr, *) "If there are differences, set_case_initial_conditions is likely"
  write(fstderr, *) "not resetting everything it needs to."

  call clean_up_clubb(err_info)

  ! Clean up err_info

  if ( any(err_info%err_code == clubb_fatal_error) ) then
    error stop ""
  else
    write(fstderr,*) "Program exited normally"
    call exit(success_code)
  end if
  call cleanup_err_info_api(err_info)

end program clubb_driver_test
