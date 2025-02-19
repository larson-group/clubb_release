!----------------------------------------------------------------------
! $Id$

program clubb_standalone

! Description:
! This is essentially a minimalist frontend for CLUBB.
!
!----------------------------------------------------------------------

  use clubb_driver, only: run_clubb ! Procedure(s)

  use error_code, only: &
        clubb_no_error, &               ! Constants
        clubb_fatal_error

  use parameter_indices, only: nparams ! Variable(s)

  use parameters_tunable, only: &
      set_default_parameters, & ! Procedure(s)
      init_clubb_params

  use clubb_precision, only: core_rknd ! Constant

  use constants_clubb, only: fstderr ! Constant

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
    iostat, &
    err_code

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
  call init_clubb_params( ngrdcol, iunit, namelist_filename, &
                          clubb_params )

  ! Initialize status of run 
  err_code = clubb_no_error

  ! Run the model
  call run_clubb( ngrdcol, calls_per_out, l_output_multi_col, l_output_double_prec, &
                  clubb_params, namelist_filename, l_stdout, err_code )

  if ( err_code == clubb_fatal_error ) then
    error stop "Fatal error in clubb, check your parameter values and timestep"
  else
    write(fstderr,*) "Program exited normally"
    call exit(success_code)
  end if 

end program clubb_standalone
