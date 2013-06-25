!-------------------------------------------------------------------------------
! $Id$

program clubb_thread_test

! Description:
!   This is a unit test to CLUBB to determine whether clubb is threadsafe
!
!-------------------------------------------------------------------------------

  use clubb_driver, only: run_clubb ! Procedure(s)
  
  use error_code, only: clubb_no_error ! Variable(s)

  use error_code, only: fatal_error ! Procedure(s)

  use parameter_indices, only: nparams ! Variable(s)

  use parameters_tunable, only: read_parameters ! Procedure(s)

  use clubb_precision, only: core_rknd ! Variable(s)

  implicit none

  ! External
  intrinsic :: size

#ifdef _OPENMP
  integer :: omp_get_thread_num
#endif

  ! Constant Parameters

  ! Text files containing namelists
  ! concatenated from various input files such as
  ! model.in, tunable_parameters.in, error....in.
  character(len=13), dimension(2), parameter :: &
    namelist_filename = (/"clubb_1.in", "clubb_2.in" /)

  logical, parameter :: &
    l_stdout = .false.
  
  ! Local Variables
  ! Run information
  real( kind = core_rknd ), dimension(nparams) :: & 
    params  ! Array of the model constants

  ! Internal variables
  integer, dimension(2) :: err_code

  integer :: iter, iunit

  !-----------------------------------------------------------------------------

  ! --- Begin Code ---

  ! Initialize status of run 
! err_code = clubb_no_error

! do iter = 1, size( err_code )
!   iunit = 10
!   ! Read in model parameter values
!   call read_parameters( iunit, namelist_filename(iter), params )
!   ! Run the model
!   call run_clubb( params, namelist_filename(iter), err_code(iter), l_stdout )
! end do

! if ( fatal_error( err_code(1) ) ) then
!   stop "The first simulation failed (single-threaded)"
! else if ( fatal_error( err_code(2) ) ) then
!   stop "The 2nd simulation failed (single-threaded)"
! end if

  ! Initialize status of run 
  err_code = clubb_no_error
  ! Run the model in parallel

!$omp parallel do default(shared), private(iter, params, iunit), &
!$omp   shared(err_code)
  do iter = 1, size( err_code )
#ifdef _OPENMP
    iunit = omp_get_thread_num() + 10
#else
    iunit = 10
#endif
    ! Read in model parameter values
    call read_parameters( iunit, namelist_filename(iter), params )
    ! Run the model
    call run_clubb( params, namelist_filename(iter), err_code(iter), l_stdout )
  end do
!$omp end parallel do

  if ( fatal_error( err_code(1) ) ) then
    stop "The first simulation failed (multi-threaded)"
  else if ( fatal_error( err_code(2) ) ) then
    stop "The 2nd simulation failed (multi-threaded)"
  end if

end program clubb_thread_test
