!---------------------------------------------------------------------------
! $Id$
!===============================================================================
program G_unit_tests

  ! Description:
  ! The unit testing framework for various pieces of CLUBB code.

  ! References:

  ! Musical References:
  !
  !                               USIN' CLUBB
  !                               ===========
  !         (Parody to be rapped to the beat of "In Da Club" by 50 Cent)
  !
  ! Go, go, go, go, go, go, go, go CLUBB!
  ! Commit 7000!
  ! We gonna party like its a birthday!
  ! We gonna sip Bicardi like its a birthday!
  ! And we don't care cause it's commit 7000!
  !
  ! (Chorus 2x)
  ! You will find me usin' CLUBB, chewin' on a sub,
  ! Mumma, I got what ya need if you need beneath the grid.
  ! I'm into PDFs, I ain't into mass flux junk,
  ! So come checkout CLUBB, if ya like the stuff we did.
  !
  ! (Verse)
  ! My code's as pretty as a Benz on Dubs,
  ! But the Sun Compiler's squaking -- always drama usin' CLUBB.
  ! When we match the LES everybody shows us love,
  ! And when we don't have low clouds nobody shows us love.
  ! But homie ain't nothing change, ...
  !
  ! Need more verses ... place yours here!
  ! 
  !
  ! (Chorus 2x)
  !
  ! (Bridge)
  ! The micro handles rain, ice, and snow,
  ! plus we got all our fancy schemes,
  ! An-al-y-tic micro and SILHS,
  ! We're best on the subgrid and that won't change.
  !
  ! (Verse)
  ! And you're gonna love it, way more than you hate it,
  ! Oh, you mad?  I thought you'd be happy we made it.
  ! We're the best darn param. for subgrid clouds,
  ! And when you use us for micro ya know you'll be wowwed.
  !
  ! Need more verses ... place yours here!
  !
  !
  ! (Chorus 2x)
  !
  ! Don't try to act like you don't know what we be usin' either,
  ! We be usin' CLUBB all the time, it's about the blast off!
  ! G-UNIT (tests)
  !-------------------------------------------------------------------------

  use constants_clubb, only: &
      fstdout  ! Constant(s)

  use KK_integrals_tests, only: &
      KK_integrals_tests_driver  ! Procedure(s)

  use corr_cholesky_mtx_tests, only: &
      corr_cholesky_mtx_tests_driver  ! Procedure(s)

  use hole_filling_tests, only: &
      hole_filling_tests_driver ! Procedure(s)

  use Nc_Ncn_test, only: &
      Nc_Ncn_unit_test ! Procedure(s)

  use read_corr_mtx_test, only: &
      read_corr_mtx_unit_test ! Procedure(s)

  implicit none

  ! Local Constants
  integer, parameter :: iunit = 25

  ! Local Variables
  integer :: exit_code = 0 ! Assume all of the tests worked.

  ! Initialize flags
  logical :: &
    l_KK_unit_tests = .true.,           & ! Flag for KK integrals tests
    l_corr_cholesky_mtx_tests = .true., & ! Flag for corr_cholesky_mtx_tests
    l_hole_filling_tests = .true.,      & ! Flag for hole filling tests
    l_Nc_Ncn_test = .true.,             & ! Flag for Nc-Ncn Equations tests
    l_read_corr_mtx_test = .true.,      & ! Flag for corr matrix read test
      show_read_test_arrays = .true.      ! If true, the arrays used in
                                          !   read_corr_mtx_test will be shown

  ! Definition of namelist
  namelist /G_unit_namelist/ &
    l_KK_unit_tests, l_corr_cholesky_mtx_tests, l_hole_filling_tests, &
    l_Nc_Ncn_test, l_read_corr_mtx_test


  ! Read namelist file
  open(unit=iunit, file="G_unit_tests.in", status='old')
  read(unit=iunit, nml=G_unit_namelist)
  close(unit=iunit)


  write(fstdout,'(A)') "Running G_unit_tests"
  write(fstdout,'(A)') " "


  if ( l_KK_unit_tests ) then
     if (KK_integrals_tests_driver() /= 0) then
       exit_code = 1
     end if
  end if

  if ( l_corr_cholesky_mtx_tests ) then
     if (corr_cholesky_mtx_tests_driver() /= 0) then
       exit_code = 1
     end if
  end if

  if ( l_hole_filling_tests ) then
     if (hole_filling_tests_driver() /= 0) then
       exit_code = 1
     end if
  end if

  if ( l_Nc_Ncn_test ) then
     if (Nc_Ncn_unit_test() /= 0) then
       exit_code = 1
     end if
  end if

  if ( l_read_corr_mtx_test ) then
     if (read_corr_mtx_unit_test(show_read_test_arrays) /= 0) then
       exit_code = 1
     end if
  end if

  ! Stop with exit code if error found
  if (exit_code /= 0) then
    call exit(1)
  end if

!===============================================================================

end program G_unit_tests

