!----------------------------------------------------------------------
! $Id$

program clubb_standalone

! Description:
! This is essentially a minimalist frontend for CLUBB.
!
!----------------------------------------------------------------------

  use clubb_driver, only: run_clubb ! Procedure(s)

  use error_code, only: &
      clubb_generalized_grd_test_err, & ! Constants
      clubb_fatal_error

  use model_flags, only: &
      l_test_grid_generalization

  use constants_clubb, only: fstderr ! Constant

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

  ! Read in model parameter values and run the model.
  call run_clubb( trim(namelist_filename), l_stdout, err_info )

  if ( .not. l_test_grid_generalization ) then
    ! Normal CLUBB run
    if ( any(err_info%err_code == clubb_fatal_error) ) then
      error stop "Fatal error in clubb, check your parameter values and timestep"
    else
      write(fstderr,*) "Program exited normally"
      call exit(success_code)
    end if
  else ! l_test_grid_generalization
    ! CLUBB grid generalization test
    ! A different success or error return is required
    if ( any(err_info%err_code == clubb_generalized_grd_test_err) ) then
      error stop "Error in generalized grid test; check the error messages"
    else
      write(fstderr,*) "Generalized grid test passed"
      call exit(success_code)
    end if
  endif

  ! Clean up err_info
  call cleanup_err_info_api(err_info)

end program clubb_standalone
