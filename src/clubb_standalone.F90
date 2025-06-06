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

  use parameter_indices, only: nparams ! Variable(s)

  use parameters_tunable, only: &
      set_default_parameters, & ! Procedure(s)
      init_clubb_params_api

  use clubb_precision, only: core_rknd ! Constant

  use constants_clubb, only: fstderr ! Constant

  use err_info_type_module, only: &
    err_info_type,                  & ! Type
    init_default_err_info_api,      & ! Procedure(s)
    cleanup_err_info_api

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

  ! Initialize err_info with default values
  call init_default_err_info_api(ngrdcol, err_info)

  allocate( clubb_params(ngrdcol,nparams) )

  !Read in model parameter values
  call init_clubb_params_api( ngrdcol, iunit, namelist_filename, &
                              clubb_params )

  ! Run the model
  call run_clubb( ngrdcol, calls_per_out, l_output_multi_col, l_output_double_prec, &
                  clubb_params, namelist_filename, l_stdout, err_info )

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
