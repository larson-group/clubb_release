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
  
  use error_code, only: &
        err_code,                    & ! Error Indicator
        clubb_no_error,              & ! Constant
        clubb_fatal_error

  use parameter_indices, only: nparams ! Variable(s)

  use parameters_tunable, only: &
      init_clubb_params ! Procedure(s)

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
    l_stdout = .false., &
    l_output_multi_col = .false.
  
  ! Local Variables
  ! Run information
  real( kind = core_rknd ), dimension(1,nparams) :: & 
    clubb_params  ! Array of the model constants

  ! Internal variables
  integer, dimension(ncases) :: err_code_saves

  integer :: iter, iunit

  !-----------------------------------------------------------------------------

  ! --- Begin Code ---

#ifndef _OPENMP
  error stop "This program needs to be compiled with OpenMP enabled to test if CLUBB is threadsafe"
#endif

  ! Initialize status of run 
  err_code_saves = clubb_no_error

  ! Run the model in parallel
!$omp parallel do default(shared), private(iter, clubb_params, iunit), &
!$omp   shared(err_code_saves)
  do iter = 1, ncases
#ifdef _OPENMP
    iunit = omp_get_thread_num() + 10
#else
    iunit = 10
#endif
    ! Read in model parameter values
    call init_clubb_params( 1, iunit, namelist_filename(iter), &
                          clubb_params )

    ! Run the model
    call run_clubb( 1, 1, clubb_params, namelist_filename(iter), l_stdout, l_output_multi_col )

    err_code_saves(iter) = err_code

  end do ! 1 .. ncases
!$omp end parallel do

  do iter = 1, ncases
    if ( any( err_code_saves(:) == clubb_fatal_error ) ) then
      write(fstderr,*) "Simulation ", iter, " failed (multi-threaded)"
    end if
  end do

end program clubb_thread_test
!-------------------------------------------------------------------------------
