!----------------------------------------------------------------------
! $Id$

program clubb_standalone

! Description:
! This is essentially a minimalist frontend for CLUBB.
!
!----------------------------------------------------------------------

  use clubb_driver, only: run_clubb ! Procedure(s)
  
  use error_code, only: clubb_no_error ! Variable(s)

  use error_code, only: fatal_error ! Procedure(s)

  use parameter_indices, only: nparams ! Variable(s)

  use parameters_tunable, only: read_parameters ! Procedure(s)

  use clubb_precision, only: core_rknd ! Constant

  use constants_clubb, only: fstderr ! Constant

  implicit none

  ! External
  intrinsic :: trim

  ! Constant parameters
  integer, parameter :: iunit = 10

  integer, parameter :: success_code = 6

  character(len=13), parameter :: &
    namelist_filename = "clubb.in"  ! Text file containing namelists
                                    ! concatenated from various input files such as
                                    ! model.in, tunable_parameters.in, error....in.

  logical, parameter :: &
    l_stdout = .true.

  ! Run information
  real( kind = core_rknd ), dimension(nparams) :: & 
    params  ! Array of the model constants

  ! Internal variables
  integer :: err_code 

!-----------------------------------------------------------------------

  ! --- Begin Code ---

  ! Read in model parameter values
  call read_parameters( iunit, namelist_filename, params )

  ! Initialize status of run 
  err_code = clubb_no_error

  ! Run the model
  call run_clubb( params, namelist_filename, l_stdout, err_code )

  if ( fatal_error( err_code ) ) then
    stop "Model wasn't valid, check your parameters and timestep"

  else
    write(fstderr,*) "Program exited normally"
    call exit(success_code)
  end if 

end program clubb_standalone