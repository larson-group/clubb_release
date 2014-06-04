!-------------------------------------------------------------------------------
! $Id$

program clubb_thread_test

! Description:
!   This is a unit test to CLUBB to determine whether clubb is threadsafe
!
! References:
!   None
!
! Notes:
!   This program doesn't work consistently with the version of PGI or GNU
!   Fortran at UWM.  Use Intel Fortran to test it.
!-------------------------------------------------------------------------------

  use clubb_driver, only: run_clubb ! Procedure(s)
  
  use error_code, only: clubb_no_error ! Variable(s)

  use error_code, only: fatal_error ! Procedure(s)

  use parameter_indices, only: nparams ! Variable(s)

  use parameters_tunable, only: read_parameters ! Procedure(s)

  use clubb_precision, only: core_rknd ! Variable(s)

  use constants_clubb, only: fstderr ! Constant(s)

  implicit none

#ifdef _OPENMP
  integer :: omp_get_thread_num
#endif

  ! Constant Parameters

  integer, parameter :: ncases = 4 ! Number of cases to run

  ! Text files containing namelists
  ! concatenated from various input files such as
  ! model.in, tunable_parameters.in, error....in.
  character(len=10), dimension(ncases), parameter :: &
    namelist_filename = (/"clubb_1.in", "clubb_2.in", "clubb_3.in", "clubb_4.in" /)

  logical, parameter :: &
    l_stdout = .false.
  
  ! Local Variables
  ! Run information
  real( kind = core_rknd ), dimension(nparams) :: & 
    params  ! Array of the model constants

  ! Internal variables
  integer, dimension(ncases) :: err_code

  integer :: iter, iunit

  !-----------------------------------------------------------------------------

  ! --- Begin Code ---

#ifndef _OPENMP
  stop "This program needs to be compiled with OpenMP enabled to test if CLUBB is threadsafe"
#endif

  ! Initialize status of run 
  err_code = clubb_no_error

  ! Run the model in parallel
!$omp parallel do default(shared), private(iter, params, iunit), &
!$omp   shared(err_code)
  do iter = 1, ncases
#ifdef _OPENMP
    iunit = omp_get_thread_num() + 10
#else
    iunit = 10
#endif
    ! Read in model parameter values
    call read_parameters( iunit, namelist_filename(iter), params )
    ! Run the model
    call run_clubb( params, namelist_filename(iter), l_stdout, err_code(iter) )
  end do ! 1 .. ncases
!$omp end parallel do

  do iter = 1, ncases
    if ( fatal_error( err_code(iter) ) ) then
      write(fstderr,*) "Simulation ", iter, " failed (multi-threaded)"
    end if
  end do

end program clubb_thread_test
!-------------------------------------------------------------------------------
