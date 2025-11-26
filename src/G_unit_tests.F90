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
  !
  ! -----
  !
  ! (Tune switches to "Handlebars" by Flobots)
  !
  ! I can ride my CLUBB going standalone.
  !                     going standalone.
  !                     going standalone.
  ! I can ride my CLUBB with no host model.
  !                     with no host model.
  !                     with no host model.
  !
  ! Hey, hey, look at me. Advancing prognostic moments like I'm supposed to be.
  ! On top... and the bottom of the domain; staggered grid for conveniency.
  ! I can write code to conserve water; my residual has no negativity.
  ! I discretize implicitly to prevent numerical instability!
  !
  ! I can advance a point-based microphysics scheme; no I don't ignore subgrid variability!
  ! I can upscale it by numerical integration, or derive the integral analytically.
  ! Just give me some sounding terms... once I find the right PDF family,
  ! I can close those higher order moments, so making assumptions there's no need cause
  !
  ! I can resolve a grid box with a PDF!
  !                          with a PDF!
  !                          with a PDF!
  !                          with a PDF!
  !                          with a PDFFFFFFFFFFFFFFFFFFF!!!!!!!!!!!!!!!!
  ! (interlude)
  !
  ! -----
  !
  !  (Tune switches to "Singing in the Rain" by Ignacio Herbert Brown)
  !
  !  I'm singing in the rain,
  !  just singing (and parameterizing subgrid variability) in the rain
  !  What a glorious feeling (when the Bitten tests run and then)
  !  I'm happy again.
  !  I'm laughing at clouds
  !  so dark up above.
  !  The sun's in my heart
  !  And I'm ready for love (or at least another simulation of RICO).
  !  Let the stormy clouds (from ARM97) chase
  !  everyone (especially parameterizations lacking a rigorous theoretical basis) from the place.
  !  (They have too much convective inhibition
  !  to produce satisfactory results in the GCSS intercomparison.)
  !  I walk down the lane
  !  with a happy refrain,
  !  just singing,
  !  singing (and also parameterizing) in the rain.
  !
  !  Dancing in the rain (because my co-workers were tired of my singing and kicked me out of W404)
  !  Dee-ah dee-ah dee-ah
  !  Dee-ah dee-ah dee-ah
  !  I'm happy again (as long as the 2 lines above weren't a new and unfamiliar compiler warning)!
  !  I'm singing and dancing in the rain (are any of the doors to EMS unlocked on a weekend?). 
  !-------------------------------------------------------------------------

  use constants_clubb, only: &
      fstdout  ! Constant(s)

  use model_flags, only: &
      iiPDF_new,        & ! Variable(s)
      iiPDF_ADG1,       &
      iiPDF_TSDADG,     &
      iiPDF_LY93,       &
      iiPDF_new_hybrid

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

  use silhs_category_test, only: &
      silhs_category_test_driver ! Procedure

  use mu_sigma_hm_tests, only: &
      mu_sigma_hm_unit_tests  ! Procedure(s)

  use pdf_parameter_tests, only: &
      pdf_parameter_unit_tests  ! Procedure(s)

  use spurious_source_test, only: &
      spurious_source_unit_test  ! Procedure(s)

  use tuner_tests, only: &
      tuner_tests_driver        ! Procedure
      
  use w_up_in_cloud_tests, only: &
      w_up_in_cloud_tests_driver
      
  use smooth_heaviside_tests, only: &
      smooth_heaviside_tests_driver
      
  use smooth_min_max_tests, only: &
      smooth_min_max_tests_driver

  use rev_direction_grid_test, only: &
      rev_direction_grid_unit_test    ! Procedure
      
  use fill_holes_tests, only: &
      fill_holes_tests_driver    ! Procedure

  use grid_class, only: grid ! Type

  use stats_type, only: stats ! Type

  implicit none

  type(grid), target :: gr

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
    l_silhs_category_test = .true.,     & ! Flag for silhs category test
    show_read_test_arrays = .true.,     & ! If true, the arrays used in
                                          !   read_corr_mtx_test will be shown
    l_mu_sigma_hm_tests = .true.,       & ! Flag for the hydromet mu/sigma tests
    l_pdf_parameter_tests = .true.,     & ! Flag for the PDF parameter tests
    l_spurious_source_test = .true.,    & ! Flag for the spurious source test
    l_tuner_tests = .true.,             & ! Flag for the tuner tests
    l_w_up_in_cloud_test = .true.,      & ! Flag for the calc_w_up_in_cloud test
    l_smooth_heaviside_test = .true.,   & ! Flag for the smooth_heaviside test
    l_smooth_min_max_test = .true.,     & ! Flag for the smooth_min_max test
    l_rev_direction_grid_test = .true., & ! Flag for reverse direction grid test
    l_fill_holes_test = .true.            ! Flag for hole filling test

  ! Definition of namelist
  namelist /G_unit_namelist/ &
    l_KK_unit_tests, l_corr_cholesky_mtx_tests, l_hole_filling_tests, &
    l_Nc_Ncn_test, l_read_corr_mtx_test, l_silhs_category_test, &
    l_mu_sigma_hm_tests, l_pdf_parameter_tests, l_spurious_source_test, &
    l_tuner_tests, l_w_up_in_cloud_test, l_smooth_heaviside_test, &
    l_smooth_min_max_test, l_rev_direction_grid_test, l_fill_holes_test


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

  if ( l_silhs_category_test ) then
     if (silhs_category_test_driver() /= 0) then
       exit_code = 1
     end if
  end if

  if ( l_mu_sigma_hm_tests ) then
     if ( mu_sigma_hm_unit_tests( ) /= 0 ) then
        exit_code = 1
     endif
  endif

  if ( l_pdf_parameter_tests ) then
     if ( pdf_parameter_unit_tests( gr, iiPDF_ADG1 ) /= 0 ) then
        exit_code = 1
     endif
     if ( pdf_parameter_unit_tests( gr, iiPDF_LY93 ) /= 0 ) then
        exit_code = 1
     endif
     if ( pdf_parameter_unit_tests( gr, iiPDF_TSDADG ) /= 0 ) then
        exit_code = 1
     endif
     if ( pdf_parameter_unit_tests( gr, iiPDF_new ) /= 0 ) then
        exit_code = 1
     endif
     if ( pdf_parameter_unit_tests( gr, iiPDF_new_hybrid ) /= 0 ) then
        exit_code = 1
     endif
  endif

  if ( l_spurious_source_test ) then
     if ( spurious_source_unit_test( ) /= 0 ) then
        exit_code = 1
     endif
  endif

  if ( l_tuner_tests ) then
     if ( tuner_tests_driver( ) /= 0 ) then
        exit_code = 1
     endif
  endif
  
  if ( l_w_up_in_cloud_test ) then
     if ( w_up_in_cloud_tests_driver(gr) /= 0 ) then
        exit_code = 1
     endif
  endif
  
  if ( l_smooth_heaviside_test ) then
     if ( smooth_heaviside_tests_driver() /= 0 ) then
        exit_code = 1
     endif
  endif
  
  if ( l_smooth_min_max_test ) then
     if ( smooth_min_max_tests_driver() /= 0 ) then
        exit_code = 1
     endif
  endif

  if ( l_rev_direction_grid_test ) then
     if ( rev_direction_grid_unit_test() /= 0 ) then
        exit_code = 1
     endif
  endif

  if ( l_fill_holes_test ) then
     if ( fill_holes_tests_driver() /= 0 ) then
        exit_code = 1
     endif
  endif

  ! Stop with exit code if error found
  if (exit_code /= 0) then
    call exit(1)
  end if

!===============================================================================

end program G_unit_tests
