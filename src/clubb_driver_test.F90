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

  use parameter_indices, only: &
    nparams ! Variable(s)

  use parameters_tunable, only: &
    set_default_parameters, & ! Procedure(s)
    init_clubb_params_api

  use clubb_precision, only: &
    core_rknd ! Constant

  use constants_clubb, only: &
    fstderr, & ! Constant
    fstdout

  use err_info_type_module, only: &
    err_info_type,                  & ! Type
    cleanup_err_info_api              ! Procedure(s)
    
  implicit none

  ! External
  intrinsic :: trim

  ! Constant parameters
  integer, parameter :: &
    iunit = 10

  integer, parameter :: &
    success_code = 6

  integer :: &
    ngrdcol, &
    calls_per_out, &
    iostat

  type(err_info_type) :: &
    err_info        ! err_info struct containing err_code and err_header
                    ! Initialization is done within run_clubb

  character(len=13), parameter :: &
    namelist_filename = "clubb.in"  ! Text file containing namelists
                                    ! concatenated from various input files such as
                                    ! model.in, tunable_parameters.in, error....in.
  logical, parameter :: &
    l_stdout = .true.

  logical :: &
    l_output_multi_col, &
    l_output_double_prec    ! Flag to enable double precision

  character(len=10) :: arg

  real( kind = core_rknd ), dimension(:,:), allocatable :: & 
    clubb_params  ! Array of the model constants

  namelist /multicol_def/  & 
    ngrdcol, &
    l_output_multi_col, &
    l_output_double_prec, &
    calls_per_out

  !--------------------------------- Begin Code ---------------------------------

  ! Set default namelist values
  ngrdcol = 1
  l_output_multi_col = .false.
  l_output_double_prec = .true.
  calls_per_out = 1

  ! Read the namelist for ngrdcol only
  open(unit=iunit, file=namelist_filename, status='old', action='read')
  read(unit=iunit, iostat=iostat, nml=multicol_def)
  close(unit=iunit)

  if ( iostat /= 0 ) then
    write(fstderr,*) "multicol_def namelist not found in clubb.in -- defaulting to ngrdcol = 1"
  end if

  allocate( clubb_params(ngrdcol,nparams) )

  !Read in model parameter values
  call init_clubb_params_api( ngrdcol, iunit, namelist_filename, &
                              clubb_params )

  !---------------------------------------- Test Sections ----------------------------------------

  write(fstderr, *) "======================== REINITIALIZATION TEST ========================"
  write(fstderr, *) "This section ensures that everything allocated in init_clubb_case will be deallocated"
  write(fstderr, *) "in clean_up_clubb. This may cause a runtime error if there is a mismatch between"
  write(fstderr, *) "allocate/deallocate statements, but could "

  ! initialize, cleaup, then initialize again
  call init_clubb_case(namelist_filename, ngrdcol, clubb_params, err_info)
  call clean_up_clubb(ngrdcol, clubb_params, err_info)
  call init_clubb_case(namelist_filename, ngrdcol, clubb_params, err_info)


  write(fstderr, *) "======================== DOUBLE TIMESTEP RUN ========================"
  write(fstderr, *) "Calling advance_clubb_to_end, then set_case_initial_conditions, "
  write(fstderr, *) "then advance_clubb_to_end again. This could result in a runtime"
  write(fstderr, *) "error if an allocation statement appears in these routines, since"
  write(fstderr, *) "we don't call clean_up_clubb in between the advance_ calls."

  ! Run once with stats off. Turning stats on will output to disk, which we only want to do
  ! for the second run. This sets l_output_multi_col=.false. for the same reason
  call advance_clubb_to_end( ngrdcol, calls_per_out, clubb_params, &
                       l_stdout, .false., l_output_double_prec, &
                       err_info, &
                       l_suppress_stats=.true. )

  ! Reset model back to initial conditions
  call set_case_initial_conditions(ngrdcol, clubb_params, err_info)

  ! Turn stats on and run. If this produces different output than what calling run_clubb
  ! does, then set_case_initial_conditions is not correctly reseting to initial conditions
  call advance_clubb_to_end( ngrdcol, calls_per_out, clubb_params, &
                       l_stdout, l_output_multi_col, l_output_double_prec, &
                       err_info )

  write(fstderr, *) "======================== RUN OVER ========================"
  write(fstderr, *) "WARNING: The double timestep test is not complete until the "
  write(fstderr, *) "output is compared with the output from calling clubb_standalone."
  write(fstderr, *) "This driver should produce BFB output with clubb_standalone."
  write(fstderr, *) "Save the output, rerun with clubb_standalone, then compare."
  write(fstderr, *) "If there are differences, set_case_initial_conditions is likely"
  write(fstderr, *) "not resetting everything it needs to."

  call clean_up_clubb(ngrdcol, clubb_params, err_info)

  ! Clean up err_info

  if ( any(err_info%err_code == clubb_fatal_error) ) then
    error stop ""
  else
    write(fstderr,*) "Program exited normally"
    call exit(success_code)
  end if
  call cleanup_err_info_api(err_info)

end program clubb_driver_test
