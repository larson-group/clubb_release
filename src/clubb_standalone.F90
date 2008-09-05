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

  use parameters_tunable, only: params_list ! Variable(s)

  implicit none

  ! External
  intrinsic :: trim

  ! Constant parameters
  integer, parameter :: iunit = 10

  character(len=13), parameter :: &
    namelist_filename = "clubb.in"  ! Text file with the namelists

  logical, parameter :: &
    l_stdout       = .true., &
    l_input_fields = .false.

  ! Run information
  real, dimension(nparams) :: & 
    params  ! Array of the model constants

  ! Internal variables
  integer :: err_code 

  integer :: i ! Loop iterator

!-----------------------------------------------------------------------

  ! --- Begin Code ---

  ! Read in model parameter values
  call read_parameters( iunit, namelist_filename, params )

  ! If standard output (stdout) is selected, print the list of
  ! parameters that are being used to the screen before the run.
  if ( l_stdout ) then
    write(unit=*,fmt='(4x,A9,5x,11x,A5)') "Parameter", "Value"
    write(unit=*,fmt='(4x,A9,5x,11x,A5)') "---------", "-----"
    do i = 1, nparams, 1
       write(unit=*,fmt='(A18,F27.20)') params_list(i) // " = ", params(i)
    end do
  end if

  ! Initialize status of run 
  err_code = clubb_no_error

  ! Run the model
  call run_clubb( params, namelist_filename, err_code, l_stdout, l_input_fields )

  if ( fatal_error( err_code ) ) then
    stop "Model wasn't valid, check your parameters and timestep"

  else
    stop "Program exited normally"

  end if 

end program clubb_standalone
