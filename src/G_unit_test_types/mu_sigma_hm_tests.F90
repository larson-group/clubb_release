! $Id$
!===============================================================================
module mu_sigma_hm_tests

  implicit none

  public :: mu_sigma_hm_unit_tests

  private :: produce_seed

  private  ! default scope

  contains

  !=============================================================================
  function mu_sigma_hm_unit_tests( )

    ! Description:
    ! Unit testing framework for the code that calculates the PDF component
    ! mean (in-precip) and standard deviation (in-precip) values for
    ! hydrometeors.
    !
    ! There are 10 different PDF parameter sets specified.  For each PDF
    ! parameter set, 11 different values of omicron are used (0, 0.1, 0.2, 0.3,
    ! 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0).  Furthermore, for every value of
    ! omicron, 11 different values of zeta are used (0, -1/2, 1, -1/3, 1/2,
    ! -1/5, 1/4, -9/10, 9, -3/4, and 3).  So, for every PDF parameter set, 121
    ! combinations of omicron and zeta are used.  Overall, 1210 parameter sets
    ! (combinations of values of the PDF parameters, omicron, and zeta) are used
    ! to call the subroutine.
    !
    ! Each subroutine call produces values of mu_hm_1, mu_hm_2, sigma_hm_1,
    ! sigma_hm_2, hm_1, hm_2, sigma_hm_1_sqd_on_mu_hm_1_sqd, and
    ! sigma_hm_2_sqd_on_mu_hm_2_sqd.  These values are evaluated by five tests.
    !
    ! Test 1
    ! Recalculate the mean (overall) value of the hydrometeor using the PDF
    ! parameters calculated by the subroutine.  Compare that result to the value
    ! that was used to start with.  They should match numerically within a very
    ! small tolerance.
    !
    ! Test 2
    ! Recalculate the variance (overall) of the hydrometeor using the PDF
    ! parameters calculated by the subroutine.  Compare that result to the value
    ! that was used to start with.  They should match numerically within a very
    ! small tolerance.
    !
    ! Test 3
    ! Calculate the "ratio of ratios", which is:
    ! sigma_hm_1_sqd_on_mu_hm_1_sqd / sigma_hm_2_sqd_on_mu_hm_2_sqd.
    ! Compare that result to the value of 1 + zeta.  They should match
    ! numerically within a very small tolerance.  This is not called when
    ! omicron = 0.
    !
    ! Test 4
    ! Check that mu_hm_1 >= hm_tol / precip_frac_1,
    ! mu_hm_2 >= hm_tol / precip_frac_2, hm_1 >= hm_tol, and hm_2 >= hm_tol.
    !
    ! Test 5
    ! Check that sigma_hm_1 >= 0, sigma_hm_2 >= 0,
    ! sigma_hm_1_sqd_on_mu_hm_1_sqd >= 0,
    ! and sigma_hm_2_sqd_on_mu_hm_2_sqd >= 0.
    !
    ! If each of the 1210 parameter sets pass all of the tests, then the unit
    ! tests are passed!

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        three,         & ! Constant(s)
        one,           &
        three_fourths, &
        one_half,      &
        one_third,     &
        one_fourth,    &
        zero,          &
        fstdout,       &
        fstderr

    use setup_clubb_pdf_params, only: &
        calc_comp_mu_sigma_hm  ! Procedure(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    real( kind = core_rknd ) :: &
      hmm,                & ! Hydrometeor mean (overall), <hm>           [hm un]
      hmm_ip,             & ! Hydrometeor mean (in-precip.), <hm|_ip>    [hm un]
      hmp2,               & ! Hydrometeor variance (overall), <hm'^2>  [hm un^2]
      hmp2_ip_on_hmm2_ip, & ! Ratio <hm|_ip'^2> / <hm|_ip>^2                 [-]
      mixt_frac,          & ! Mixture fraction                               [-]
      precip_frac,        & ! Precipitation fraction (overall)               [-]
      precip_frac_1,      & ! Precipitation fraction (1st PDF component)     [-]
      precip_frac_2,      & ! Precipitation fraction (2nd PDF component)     [-]
      hm_tol,             & ! Tolerance value of hydrometeor             [hm un]
      precip_frac_tol       ! Min. precip. frac. when hydromet. are present  [-]

    real( kind = core_rknd ) :: &
      omicron,        & ! Relative width parameter, omicron = R / Rmax       [-]
      zeta_vrnce_rat    ! Width parameter for sigma_hm_1^2 / mu_hm_1^2       [-]

    real( kind = core_rknd ) :: &
      mu_hm_1,    & ! Mean of hm (1st PDF component) in-precip (ip)      [hm un]
      mu_hm_2,    & ! Mean of hm (2nd PDF component) ip                  [hm un]
      sigma_hm_1, & ! Standard deviation of hm (1st PDF component) ip    [hm un]
      sigma_hm_2, & ! Standard deviation of hm (2nd PDF component) ip    [hm un]
      hm_1,       & ! Mean of hm (1st PDF component)                     [hm un]
      hm_2          ! Mean of hm (2nd PDF component)                     [hm un]

    real( kind = core_rknd ) :: &
      sigma_hm_1_sqd_on_mu_hm_1_sqd, & ! Ratio sigma_hm_1**2 / mu_hm_1**2    [-]
      sigma_hm_2_sqd_on_mu_hm_2_sqd    ! Ratio sigma_hm_2**2 / mu_hm_2**2    [-]

    real( kind = core_rknd ) :: &
      recalc_hmm,   & ! Recalculation of <hm> using PDF parameters       [hm un]
      recalc_hmp2,  & ! Recalculation of <hm'^2> using PDF parameters  [hm un^2]
      calc_rat_rat    ! Calculation of the "ratio of ratios"                 [-]

    ! Tiny tolerance for acceptable numerical difference between two results.
    real( kind = core_rknd ), parameter :: &
      tol = 1.0e-11_core_rknd

    logical :: &
      l_pass_test_1, & ! Flag for passing test 1
      l_pass_test_2, & ! Flag for passing test 2
      l_pass_test_3, & ! Flag for passing test 3
      l_pass_test_4, & ! Flag for passing test 4
      l_pass_test_5    ! Flag for passing test 4

    integer :: &
      mu_sigma_hm_unit_tests  ! Returns pass or fail

    integer :: &
      num_failed_sets    ! Records the number of failed parameter sets

    integer :: &
      seed_size    ! The size of the random seed array expected by the system

    integer, dimension(:), allocatable :: &
      seed_vals    ! Values used to seed the random number generator

    real( kind = core_rknd ) :: &
      precip_frac_1_min, & ! Minimum value for precip_frac_1 in parameter set 10
      rand1,             & ! Random number 1 used for PDF parameter set 10
      rand2,             & ! Random number 2 used for PDF parameter set 10
      rand3,             & ! Random number 3 used for PDF parameter set 10
      rand4,             & ! Random number 4 used for PDF parameter set 10
      rand5                ! Random number 5 used for PDF parameter set 10

    integer, parameter :: &
      num_param_sets = 10,     & ! Number of different PDF parameter sets used
      num_omicron = 11,        & ! Number of different values of omicron used
      num_zeta_vrnce_rat = 11    ! Number of different values of zeta used

    integer :: &
      iter_param_sets,     & ! Loop index for PDF parameter set
      iter_omicron,        & ! Loop index for value of omicron
      iter_zeta_vrnce_rat    ! Loop index for value of zeta


    write(fstdout,*) ""
    write(fstdout,*) "Performing hydrometeor mu and sigma values unit test"
    write(fstdout,*) "===================================================="
    write(fstdout,*) ""

    ! Initialize number of failed parameter sets.
    num_failed_sets = 0

    ! Perform unit tests for situations where there is precipitation found in
    ! both PDF components (precip_frac_1 > 0 and precip_frac_2 > 0).
    do iter_param_sets = 1, num_param_sets, 1

       if ( iter_param_sets == 1 ) then
          write(fstdout,*) "PDF parameter set 1:"
          hmm = 2.0e-5_core_rknd
          mixt_frac = 0.027_core_rknd
          precip_frac = 0.08_core_rknd
          precip_frac_1 = 1.0_core_rknd
          hmp2_ip_on_hmm2_ip = 3.6_core_rknd
          hm_tol = 1.0e-10_core_rknd
       elseif ( iter_param_sets == 2 ) then
          write(fstdout,*) "PDF parameter set 2:"
          hmm = 5.0e-6_core_rknd
          mixt_frac = 0.5_core_rknd
          precip_frac = 0.5_core_rknd
          precip_frac_1 = 0.5_core_rknd
          hmp2_ip_on_hmm2_ip = 1.0_core_rknd
          hm_tol = 1.0e-10_core_rknd
       elseif ( iter_param_sets == 3 ) then
          write(fstdout,*) "PDF parameter set 3:"
          hmm = 1.0e-10_core_rknd
          mixt_frac = 0.25_core_rknd
          precip_frac = 0.4_core_rknd
          precip_frac_1 = 0.75_core_rknd
          hmp2_ip_on_hmm2_ip = 0.5_core_rknd
          hm_tol = 1.0e-10_core_rknd
       elseif ( iter_param_sets == 4 ) then
          write(fstdout,*) "PDF parameter set 4:"
          hmm = 4.0e-5_core_rknd
          mixt_frac = 0.01_core_rknd
          precip_frac = 0.1_core_rknd
          precip_frac_1 = 0.01_core_rknd
          hmp2_ip_on_hmm2_ip = 2.5_core_rknd
          hm_tol = 1.0e-10_core_rknd
       elseif ( iter_param_sets == 5 ) then
          write(fstdout,*) "PDF parameter set 5:"
          hmm = 4.0e-5_core_rknd
          mixt_frac = 0.99_core_rknd
          precip_frac = 0.1_core_rknd
          precip_frac_1 = ( precip_frac - 0.0001_core_rknd ) / mixt_frac
          hmp2_ip_on_hmm2_ip = 2.5_core_rknd
          hm_tol = 1.0e-10_core_rknd
       elseif ( iter_param_sets == 6 ) then
          write(fstdout,*) "PDF parameter set 6:"
          hmm = 1.0e-5_core_rknd
          mixt_frac = 0.125_core_rknd
          precip_frac = 1.0_core_rknd
          precip_frac_1 = 1.0_core_rknd
          hmp2_ip_on_hmm2_ip = 1.5_core_rknd
          hm_tol = 1.0e-10_core_rknd
       elseif ( iter_param_sets == 7 ) then
          write(fstdout,*) "PDF parameter set 7:"
          hmm = 6.0e-5_core_rknd
          mixt_frac = 0.05_core_rknd
          precip_frac = 0.11_core_rknd
          precip_frac_1 = 0.9_core_rknd
          hmp2_ip_on_hmm2_ip = 5.0_core_rknd
          hm_tol = 1.0e-10_core_rknd
       elseif ( iter_param_sets == 8 ) then
          write(fstdout,*) "PDF parameter set 8:"
          hmm = 8.0e-7_core_rknd
          mixt_frac = 0.8_core_rknd
          precip_frac = 0.25_core_rknd
          precip_frac_1 = 0.1_core_rknd
          hmp2_ip_on_hmm2_ip = 0.75_core_rknd
          hm_tol = 1.0e-10_core_rknd
       elseif ( iter_param_sets == 9 ) then
          write(fstdout,*) "PDF parameter set 9:"
          hmm = 1.35e-5_core_rknd
          mixt_frac = 0.1_core_rknd
          precip_frac = 0.10_core_rknd
          precip_frac_1 = 0.5_core_rknd
          hmp2_ip_on_hmm2_ip = 2.0_core_rknd
          hm_tol = 1.0e-10_core_rknd
       elseif ( iter_param_sets == 10 ) then
          write(fstdout,*) "PDF parameter set 10 (randomly generated):"
          call random_seed( size=seed_size )
          allocate( seed_vals(1:seed_size) )
          seed_vals = produce_seed( seed_size )
          write(fstdout,*) "Random seed values = ", seed_vals
          call random_seed( put=seed_vals )
          deallocate( seed_vals )
          call random_number( rand1 )
          call random_number( rand2 )
          call random_number( rand3 )
          call random_number( rand4 )
          call random_number( rand5 )
          ! The value of hmm can range from 1.0 x 10^-6 to 5.1 x 10^-5.
          hmm = 5.0e-5_core_rknd * rand1 + 1.0e-6_core_rknd
          ! The value of mixt_frac can range from 0.01 to 0.99.
          mixt_frac = 0.98_core_rknd * rand2 + 0.01_core_rknd
          ! The value of precip_frac must be greater than mixt_frac (so that
          ! precipitation must be found in both PDF components).  The upper
          ! limit of precip_frac is 1.
          precip_frac = mixt_frac + epsilon( mixt_frac ) &
                        + rand3 * ( one - ( mixt_frac + epsilon( mixt_frac ) ) )
          ! The minimum value of precip_frac_1 is set such that the resulting
          ! precip_frac_2 is not greater than 1.  Alternatively, it must have a
          ! value of at least 0.01.
          precip_frac_1_min &
          = max( ( precip_frac - ( one - mixt_frac ) ) / mixt_frac, &
                 0.01_core_rknd )
          ! The value of precip_frac_1 can range from precip_frac_1_min to 1.
          precip_frac_1 &
          = precip_frac_1_min + rand4 * ( one - precip_frac_1_min )
          ! The value of hmp2_ip_on_hmm2_ip can range from 0.2 to 5.0.
          hmp2_ip_on_hmm2_ip = 0.2_core_rknd + rand5 * 4.8_core_rknd
          ! The value of hm_tol remains constant.
          hm_tol = 1.0e-10_core_rknd
       endif ! iter_param_sets == index

       ! Calculate precip_frac_2.
       precip_frac_2 = ( precip_frac - mixt_frac * precip_frac_1 ) &
                       / ( one - mixt_frac )

       ! Set the minimum precipitation fraction allowable when hydrometeors are
       ! present, precip_frac_tol.
       precip_frac_tol = 0.01_core_rknd

       ! Calculate <hm|_ip>.
       hmm_ip = hmm / precip_frac

       ! Calculate <hm'^2>
       hmp2 = precip_frac * ( hmp2_ip_on_hmm2_ip + one ) * hmm_ip**2 - hmm**2

       ! Print PDF parameters
       write(fstdout,*) "hmm = ", hmm
       write(fstdout,*) "mixt_frac = ", mixt_frac
       write(fstdout,*) "precip_frac = ", precip_frac
       write(fstdout,*) "precip_frac_1 = ", precip_frac_1
       write(fstdout,*) "precip_frac_2 = ", precip_frac_2
       write(fstdout,*) "hmp2_ip_on_hmm2_ip = ", hmp2_ip_on_hmm2_ip
       write(fstdout,*) "hmm_ip = ", hmm_ip
       write(fstdout,*) "hmp2 = ", hmp2
       write(fstdout,*) "hm_tol = ", hm_tol
       write(fstdout,*) "precip_frac_tol = ", precip_frac_tol
       write(fstdout,*) ""

       write(fstdout,*) "Running tests for the above parameter set for all " &
                        // "combinations of omicron and zeta.  Omicron " &
                        // "values are 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, " &
                        // "0.8, 0.9, and 1.  Zeta values are 0, -1/2, 1, " &
                        // "-1/3, 1/2, -1/5, 1/4, -9/10, 9, -3/4, and 3."
       write(fstdout,*) ""

       do iter_omicron = 1, num_omicron, 1

         ! Set the value of omicron.
         omicron = 0.1_core_rknd * real( iter_omicron - 1, kind = core_rknd )

         do iter_zeta_vrnce_rat = 1, num_zeta_vrnce_rat, 1

           ! Set the value of zeta_vrnce_rat.
           if ( iter_zeta_vrnce_rat == 1 ) then
              zeta_vrnce_rat = zero
           elseif ( iter_zeta_vrnce_rat == 2 ) then
              zeta_vrnce_rat = -one_half
           elseif ( iter_zeta_vrnce_rat == 3 ) then
              zeta_vrnce_rat = one
           elseif ( iter_zeta_vrnce_rat == 4 ) then
              zeta_vrnce_rat = -one_third
           elseif ( iter_zeta_vrnce_rat == 5 ) then
              zeta_vrnce_rat = one_half
           elseif ( iter_zeta_vrnce_rat == 6 ) then
              zeta_vrnce_rat = -0.2_core_rknd
           elseif ( iter_zeta_vrnce_rat == 7 ) then
              zeta_vrnce_rat = one_fourth
           elseif ( iter_zeta_vrnce_rat == 8 ) then
              zeta_vrnce_rat = -0.9_core_rknd
           elseif ( iter_zeta_vrnce_rat == 9 ) then
              zeta_vrnce_rat = 9.0_core_rknd
           elseif ( iter_zeta_vrnce_rat == 10 ) then
              zeta_vrnce_rat = -three_fourths
           elseif ( iter_zeta_vrnce_rat == 11 ) then
              zeta_vrnce_rat = three
           endif

           ! Call the subroutine for calculating mu_hm_1, mu_hm_2, sigma_hm_1,
           ! sigma_hm_2, hm_1, hm_2, sigma_hm_1_sqd_on_mu_hm_1_sqd, and
           ! sigma_hm_2_sqd_on_mu_hm_2_sqd.
           ! This subroutine is called only when precip_frac_1 > 0 and
           ! precip_frac_2 > 0 (there is precipitation in both PDF components).
           call calc_comp_mu_sigma_hm( hmm, hmp2_ip_on_hmm2_ip, &
                                       mixt_frac, precip_frac, precip_frac_1, &
                                       precip_frac_2, hm_tol, &
                                       precip_frac_tol, &
                                       omicron, zeta_vrnce_rat, &
                                       mu_hm_1, mu_hm_2, sigma_hm_1, &
                                       sigma_hm_2, hm_1, hm_2, &
                                       sigma_hm_1_sqd_on_mu_hm_1_sqd, &
                                       sigma_hm_2_sqd_on_mu_hm_2_sqd )

           ! Recalculate the value of <hm> using the hydrometeor PDF parameters
           ! output from the subroutine.
           recalc_hmm = mixt_frac * precip_frac_1 * mu_hm_1 &
                        + ( one - mixt_frac ) * precip_frac_2 * mu_hm_2

           ! Recalculate the value of <hm'^2> using the hydrometeor PDF
           ! parameters output from the subroutine.
           recalc_hmp2 = mixt_frac * precip_frac_1 &
                         * ( mu_hm_1**2 + sigma_hm_1**2 ) &
                         + ( one - mixt_frac ) * precip_frac_2 &
                           * ( mu_hm_2**2 + sigma_hm_2**2 ) &
                         - hmm**2

           ! Calculate the "ratio of ratios".
           ! When omicron > 0, calculate the ratio of ratios:
           ! sigma_hm_1_sqd_on_mu_hm_1_sqd / sigma_hm_2_sqd_on_mu_hm_2_sqd.
           if ( omicron > zero ) then
              calc_rat_rat &
              = sigma_hm_1_sqd_on_mu_hm_1_sqd / sigma_hm_2_sqd_on_mu_hm_2_sqd
           endif ! omicron > 0

           ! Test 1
           ! Compare the original hmm to its recalculated value.
           !    | ( <hm>|_recalc - <hm> ) / <hm> |  <=  tol;
           ! which can be rewritten as:
           !    | <hm>|_recalc - <hm> |  <=  | <hm> | * tol;
           ! and since <hm> is always positive:
           !    | <hm>|_recalc - <hm> |  <=  <hm> * tol.
           if ( abs( recalc_hmm - hmm ) <= hmm * tol ) then
              l_pass_test_1 = .true.
           else
              l_pass_test_1 = .false.
              write(fstderr,*) "Test 1 failed"
              write(fstderr,*) "hmm = ", hmm
              write(fstderr,*) "recalc_hmm = ", recalc_hmm
              write(fstderr,*) ""
           endif

           ! Test 2
           ! Compare the original hmp2 to its recalculated value.
           !    | ( <hm'^2>|_recalc - <hm'^2> ) / <hm'^2> |  <=  tol;
           ! which can be rewritten as:
           !    | <hm'^2>|_recalc - <hm'^2> |  <=  | <hm'^2> | * tol;
           ! and since <hm'^2> is always positive:
           !    | <hm'^2>|_recalc - <hm'^2> |  <=  <hm'^2> * tol.
           if ( abs( recalc_hmp2 - hmp2 ) <= hmp2 * tol ) then
              l_pass_test_2 = .true.
           else
              l_pass_test_2 = .false.
              write(fstderr,*) "Test 2 failed"
              write(fstderr,*) "hmp2 = ", hmp2
              write(fstderr,*) "recalc_hp2 = ", recalc_hmp2
              write(fstderr,*) ""
           endif

           ! Test 3
           ! Compare the "ratio of ratios" to 1 + zeta.
           !    | ( calc_rat_rat - ( 1 + zeta ) ) / ( 1 + zeta ) |  <=  tol;
           ! which can be rewritten as:
           !    | calc_rat_rat - ( 1 + zeta ) |  <=  | ( 1 + zeta ) | * tol;
           ! and since 1 + zeta is always positive:
           !    | calc_rat_rat - ( 1 + zeta ) |  <=  ( 1 + zeta ) * tol.
           if ( omicron > zero ) then
              if ( abs( calc_rat_rat - ( one + zeta_vrnce_rat ) ) &
                   <= ( one + zeta_vrnce_rat ) * tol ) then
                 l_pass_test_3 = .true.
              else
                 l_pass_test_3 = .false.
                 write(fstderr,*) "Test 3 failed"
                 write(fstderr,*) "sigma_hm_1_sqd_on_mu_hm_1_sqd " &
                                  // "/ sigma_hm_2_sqd_on_mu_hm_2_sqd = ", &
                                  calc_rat_rat
                 write(fstderr,*) "1 + zeta = ", one + zeta_vrnce_rat
                 write(fstderr,*) ""
              endif
           else ! omicron = 0
              ! When omicron = 0, both ratios, sigma_hm_1_sqd_on_mu_hm_1_sqd and
              ! sigma_hm_2_sqd_on_mu_hm_2_sqd, have a value of 0.  Their ratio,
              ! the "ratio of ratios" is undefined.  This test cannot be
              ! evaluated, so automatically pass it.
              l_pass_test_3 = .true.
           endif ! omicron > 0

           ! Test 4
           ! Check that mu_hm_1, mu_hm_2, hm_1, and hm_2 have at least threshold
           ! values.
           if ( ( mu_hm_1 >= ( hm_tol / precip_frac_1 ) * (one-tol) ) &
                .and. ( mu_hm_2 >= ( hm_tol / precip_frac_2 ) * (one-tol) ) &
                .and. ( hm_1 >= hm_tol ) .and. ( hm_2 >= hm_tol ) ) then
              l_pass_test_4 = .true.
           else
              l_pass_test_4 = .false.
              write(fstderr,*) "Test 4 failed"
              write(fstderr,*) "mu_hm_1 = ", mu_hm_1
              write(fstderr,*) "hm_tol/precip_frac_1 = ", hm_tol / precip_frac_1
              write(fstderr,*) "mu_hm_2 = ", mu_hm_2
              write(fstderr,*) "hm_tol/precip_frac_2 = ", hm_tol / precip_frac_2
              write(fstderr,*) "hm_1 = ", hm_1
              write(fstderr,*) "hm_2 = ", hm_2
              write(fstderr,*) "hm_tol = ", hm_tol
              write(fstderr,*) ""
           endif

           ! Test 5
           ! Check that sigma_hm_1, sigma_hm_2, sigma_hm_1_sqd_on_mu_hm_1_sqd,
           ! and sigma_hm_2_sqd_on_mu_hm_2_sqd have a value of at least zero.
           if ( ( sigma_hm_1 >= zero ) .and. ( sigma_hm_2 >= zero ) &
                .and. ( sigma_hm_1_sqd_on_mu_hm_1_sqd >= zero ) &
                .and. ( sigma_hm_2_sqd_on_mu_hm_2_sqd >= zero ) ) then
              l_pass_test_5 = .true.
           else
              l_pass_test_5 = .false.
              write(fstderr,*) "Test 5 failed"
              write(fstderr,*) "sigma_hm_1 = ", sigma_hm_1
              write(fstderr,*) "sigma_hm_2 = ", sigma_hm_2
              write(fstderr,*) "sigma_hm_1^2 / mu_hm_1^2 = ", &
                               sigma_hm_1_sqd_on_mu_hm_1_sqd
              write(fstderr,*) "sigma_hm_2^2 / mu_hm_2^2 = ", &
                               sigma_hm_2_sqd_on_mu_hm_2_sqd
              write(fstderr,*) ""
           endif


           if ( l_pass_test_1 .and. l_pass_test_2 .and. l_pass_test_3 &
                .and. l_pass_test_4 .and. l_pass_test_5 ) then
              ! All tests pass
              num_failed_sets = num_failed_sets
           else
              ! At least one test failed
              num_failed_sets = num_failed_sets + 1
              write(fstderr,*) "At least one test or check failed for " &
                               // "the following parameter set:  "
              write(fstderr,*) "PDF parameter set index = ", iter_param_sets
              write(fstderr,*) "omicron = ", omicron
              write(fstderr,*) "zeta = ", zeta_vrnce_rat
              write(fstderr,*) ""
           endif


         enddo ! iter_zeta_vrnce_rat = 1, num_zeta_vrnce_rat, 1

       enddo ! iter_omicron = 1, num_omicron, 1

    enddo ! iter_param_sets = 1, num_param_sets, 1


    ! Perform unit tests for situations where there is precipitation found in
    ! only one PDF component.
    !
    ! Precipitation found in 1st PDF component only.
    write(fstdout,*) "Checking values with precipitation in 1st PDF " &
                     // "component only."
    write(fstdout,*) ""
    hmm = 1.0e-5_core_rknd
    mixt_frac = 0.25_core_rknd
    precip_frac = 0.10_core_rknd
    precip_frac_1 = precip_frac / mixt_frac
    precip_frac_2 = zero
    hmp2_ip_on_hmm2_ip = 3.2_core_rknd
    hm_tol = 1.0e-10_core_rknd
    precip_frac_tol = 0.01_core_rknd

    ! Calculate <hm|_ip>.
    hmm_ip = hmm / precip_frac

    ! Calculate <hm'^2>
    hmp2 = precip_frac * ( hmp2_ip_on_hmm2_ip + one ) * hmm_ip**2 - hmm**2

    ! Call the subroutine for calculating mu_hm_1, mu_hm_2, sigma_hm_1,
    ! sigma_hm_2, hm_1, hm_2, sigma_hm_1_sqd_on_mu_hm_1_sqd, and
    ! sigma_hm_2_sqd_on_mu_hm_2_sqd.
    call calc_comp_mu_sigma_hm( hmm, hmp2_ip_on_hmm2_ip, &
                                mixt_frac, precip_frac, precip_frac_1, &
                                precip_frac_2, hm_tol, &
                                precip_frac_tol, &
                                omicron, zeta_vrnce_rat, &
                                mu_hm_1, mu_hm_2, sigma_hm_1, &
                                sigma_hm_2, hm_1, hm_2, &
                                sigma_hm_1_sqd_on_mu_hm_1_sqd, &
                                sigma_hm_2_sqd_on_mu_hm_2_sqd )

    ! Recalculate the value of <hm> using the hydrometeor PDF parameters output
    ! from the subroutine.
    recalc_hmm = mixt_frac * precip_frac_1 * mu_hm_1 &
                 + ( one - mixt_frac ) * precip_frac_2 * mu_hm_2

    ! Recalculate the value of <hm'^2> using the hydrometeor PDF parameters
    ! output from the subroutine.
    recalc_hmp2 = mixt_frac * precip_frac_1 &
                  * ( mu_hm_1**2 + sigma_hm_1**2 ) &
                  + ( one - mixt_frac ) * precip_frac_2 &
                    * ( mu_hm_2**2 + sigma_hm_2**2 ) &
                  - hmm**2

    ! Test 1
    ! Compare the original hmm to its recalculated value.
    !    | ( <hm>|_recalc - <hm> ) / <hm> |  <=  tol;
    ! which can be rewritten as:
    !    | <hm>|_recalc - <hm> |  <=  | <hm> | * tol;
    ! and since <hm> is always positive:
    !    | <hm>|_recalc - <hm> |  <=  <hm> * tol.
    if ( abs( recalc_hmm - hmm ) <= hmm * tol ) then
       l_pass_test_1 = .true.
    else
       l_pass_test_1 = .false.
       write(fstderr,*) "Test 1 failed"
       write(fstderr,*) "hmm = ", hmm
       write(fstderr,*) "recalc_hmm = ", recalc_hmm
       write(fstderr,*) ""
    endif

    ! Test 2
    ! Compare the original hmp2 to its recalculated value.
    !    | ( <hm'^2>|_recalc - <hm'^2> ) / <hm'^2> |  <=  tol;
    ! which can be rewritten as:
    !    | <hm'^2>|_recalc - <hm'^2> |  <=  | <hm'^2> | * tol;
    ! and since <hm'^2> is always positive:
    !    | <hm'^2>|_recalc - <hm'^2> |  <=  <hm'^2> * tol.
    if ( abs( recalc_hmp2 - hmp2 ) <= hmp2 * tol ) then
       l_pass_test_2 = .true.
    else
       l_pass_test_2 = .false.
       write(fstderr,*) "Test 2 failed"
       write(fstderr,*) "hmp2 = ", hmp2
       write(fstderr,*) "recalc_hp2 = ", recalc_hmp2
       write(fstderr,*) ""
    endif

    if ( l_pass_test_1 .and. l_pass_test_2 ) then
       ! Both tests pass
       num_failed_sets = num_failed_sets
    else
       ! At least one test failed
       num_failed_sets = num_failed_sets + 1
       write(fstderr,*) "At least one test or check failed for " &
                         // "the precip. in 1st PDF component only " &
                         // "parameter set."
       write(fstderr,*) ""
    endif

    ! Precipitation found in 2nd PDF component only.
    write(fstdout,*) "Checking values with precipitation in 2nd PDF " &
                     // "component only."
    write(fstdout,*) ""
    hmm = 1.0e-5_core_rknd
    mixt_frac = 0.25_core_rknd
    precip_frac = 0.10_core_rknd
    precip_frac_1 = zero
    precip_frac_2 = precip_frac / ( one - mixt_frac )
    hmp2_ip_on_hmm2_ip = 3.2_core_rknd
    hm_tol = 1.0e-10_core_rknd
    precip_frac_tol = 0.01_core_rknd

    ! Calculate <hm|_ip>.
    hmm_ip = hmm / precip_frac

    ! Calculate <hm'^2>
    hmp2 = precip_frac * ( hmp2_ip_on_hmm2_ip + one ) * hmm_ip**2 - hmm**2

    ! Call the subroutine for calculating mu_hm_1, mu_hm_2, sigma_hm_1,
    ! sigma_hm_2, hm_1, hm_2, sigma_hm_1_sqd_on_mu_hm_1_sqd, and
    ! sigma_hm_2_sqd_on_mu_hm_2_sqd.
    call calc_comp_mu_sigma_hm( hmm, hmp2_ip_on_hmm2_ip, &
                                mixt_frac, precip_frac, precip_frac_1, &
                                precip_frac_2, hm_tol, &
                                precip_frac_tol, &
                                omicron, zeta_vrnce_rat, &
                                mu_hm_1, mu_hm_2, sigma_hm_1, &
                                sigma_hm_2, hm_1, hm_2, &
                                sigma_hm_1_sqd_on_mu_hm_1_sqd, &
                                sigma_hm_2_sqd_on_mu_hm_2_sqd )

    ! Recalculate the value of <hm> using the hydrometeor PDF parameters output
    ! from the subroutine.
    recalc_hmm = mixt_frac * precip_frac_1 * mu_hm_1 &
                 + ( one - mixt_frac ) * precip_frac_2 * mu_hm_2

    ! Recalculate the value of <hm'^2> using the hydrometeor PDF parameters
    ! output from the subroutine.
    recalc_hmp2 = mixt_frac * precip_frac_1 &
                  * ( mu_hm_1**2 + sigma_hm_1**2 ) &
                  + ( one - mixt_frac ) * precip_frac_2 &
                    * ( mu_hm_2**2 + sigma_hm_2**2 ) &
                  - hmm**2

    ! Test 1
    ! Compare the original hmm to its recalculated value.
    !    | ( <hm>|_recalc - <hm> ) / <hm> |  <=  tol;
    ! which can be rewritten as:
    !    | <hm>|_recalc - <hm> |  <=  | <hm> | * tol;
    ! and since <hm> is always positive:
    !    | <hm>|_recalc - <hm> |  <=  <hm> * tol.
    if ( abs( recalc_hmm - hmm ) <= hmm * tol ) then
       l_pass_test_1 = .true.
    else
       l_pass_test_1 = .false.
       write(fstderr,*) "Test 1 failed"
       write(fstderr,*) "hmm = ", hmm
       write(fstderr,*) "recalc_hmm = ", recalc_hmm
       write(fstderr,*) ""
    endif

    ! Test 2
    ! Compare the original hmp2 to its recalculated value.
    !    | ( <hm'^2>|_recalc - <hm'^2> ) / <hm'^2> |  <=  tol;
    ! which can be rewritten as:
    !    | <hm'^2>|_recalc - <hm'^2> |  <=  | <hm'^2> | * tol;
    ! and since <hm'^2> is always positive:
    !    | <hm'^2>|_recalc - <hm'^2> |  <=  <hm'^2> * tol.
    if ( abs( recalc_hmp2 - hmp2 ) <= hmp2 * tol ) then
       l_pass_test_2 = .true.
    else
       l_pass_test_2 = .false.
       write(fstderr,*) "Test 2 failed"
       write(fstderr,*) "hmp2 = ", hmp2
       write(fstderr,*) "recalc_hp2 = ", recalc_hmp2
       write(fstderr,*) ""
    endif

    if ( l_pass_test_1 .and. l_pass_test_2 ) then
       ! Both tests pass
       num_failed_sets = num_failed_sets
    else
       ! At least one test failed
       num_failed_sets = num_failed_sets + 1
       write(fstderr,*) "At least one test or check failed for " &
                         // "the precip. in 2nd PDF component only " &
                         // "parameter set."
       write(fstderr,*) ""
    endif


    ! Print results and return exit code.
    if ( num_failed_sets == 0 ) then
       write(fstdout,'(1x,A)') "Success!"
       write(fstdout,*) ""
       mu_sigma_hm_unit_tests = 0 ! Exit Code = 0, Success!
    else ! num_failed_sets > 0
       write(fstdout,'(1x,A,I4,A)') "There were ", num_failed_sets, &
                                    " failed parameter sets."
       write(fstdout,*) ""
       mu_sigma_hm_unit_tests = 1 ! Exit Code = 1, Fail
    endif ! num_failed_sets = 0


    return

  end function mu_sigma_hm_unit_tests

  !=============================================================================
  function produce_seed( seed_size ) result( seed_vals )

    ! Description:
    ! Produces a random seed based on the date and time for a variety of seed
    ! sizes.

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: &
      seed_size    ! The size of the random seed array expected by the system

    ! Return Variable(s)
    integer, dimension(seed_size) :: &
      seed_vals    ! The values of the random seed array

    ! Local Variable(s)
    integer, dimension(8) :: &
      date_vals    ! Date values used to seed the random number generator

    integer :: &
      indx,   & ! Loop index
      offset    ! Array index for scenario that seed_size > 8


    ! Get the current date and time
    call date_and_time( values=date_vals )

    ! The intrisic subroutine date_and_time returns an array that has size 8.
    ! date_vals(1):  year
    ! date_vals(2):  month
    ! date_vals(3):  day
    ! date_vals(4):  time difference in minutes vs. UTC
    ! date_vals(5):  hour
    ! date_vals(6):  minute
    ! date_vals(7):  second
    ! date_vals(8):  millisecond
    if ( seed_size == 8 ) then

       ! The date_vals can be directly used as the seed_vals.
       seed_vals = date_vals

    elseif ( seed_size < 8 ) then

       ! The size of the seed_vals array is smaller than the size of the
       ! date_vals array.  Use the most highly variant date_vals in the
       ! seed_vals array.
       do indx = 1, seed_size, 1
          seed_vals(indx) = date_vals(8-seed_size+indx)
       enddo ! indx = 1, seed_size, 1

    else ! seed_size > 8

       ! The size of the seed_vals array is larger than the size of the
       ! date_vals array.

       ! Set the first 8 seed_vals to the date_vals.
       seed_vals(1:8) = date_vals

       ! Set the remaining seed_vals (array indices > 8) by adding the remainder
       ! of the loop index when divided by 8 to date_vals(8).  When the loop
       ! index is evenly divisible by 8, use 8 as the remainder (rather than 0).
       do indx = 9, seed_size, 1
          ! The value of offset is the integer remainder of indx when divided
          ! by 8.
          offset = mod( indx, 8 )
          ! When indx is evenly divisible by 8, the remainder is 0.  The array
          ! date_vals doesn't have an index 0, so set offset to 8, which
          ! otherwise wouldn't be produced by the mod command.
          if ( offset == 0 ) then
             offset = 8
          endif ! offset = 0
          ! Add date_vals(offset) to date_vals(8) to produce the seed value at
          ! indx.
          seed_vals(indx) = date_vals(8) + date_vals(offset)
       enddo ! indx = 1, seed_size, 1

    endif ! seed_size


    return

  end function produce_seed

!===============================================================================

end module mu_sigma_hm_tests
