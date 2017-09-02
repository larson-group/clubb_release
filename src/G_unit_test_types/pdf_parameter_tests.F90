! $Id$
!===============================================================================
module pdf_parameter_tests

  implicit none

  public :: pdf_parameter_unit_tests

  private :: recalc_single_var

  private  ! default scope

  contains

  !=============================================================================
  function pdf_parameter_unit_tests( )

    ! Description:
    ! Unit testing framework for the code that calculates the mixture fraction,
    ! the PDF component means, and the PDF component standard deviations for the
    ! trivariate PDF of w, rt, and theta-l.
    !
    ! There are 10 different PDF parameter sets specified.  For each PDF
    ! parameter set, 11 different values of F_w are used (0, 0.1, 0.2, 0.3,
    ! 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0).  Furthermore, for every value of
    ! F_w, 11 different values of zeta_w are used (0, -1/2, 1, -1/3, 1/2,
    ! -1/5, 1/4, -9/10, 9, -3/4, and 3).  So, for every PDF parameter set, 121
    ! combinations of F_w and zeta_w are used.  Overall, 1210 parameter sets
    ! (combinations of values of the PDF parameters, F_w, and zeta_w) are used
    ! to call the subroutine.
    !
    ! Each subroutine call produces values of mu_w_1, mu_w_2, sigma_w_1,
    ! sigma_w_2, and mixt_frac.  These values are evaluated by five tests.
    !
    ! Test 1
    ! Recalculate the mean (overall) value of the w using the PDF parameters
    ! calculated by the subroutine.  Compare that result to the value that was
    ! used to start with.  They should match numerically within a very small
    ! tolerance.
    !
    ! Test 2
    ! Recalculate the variance (overall) of the w using the PDF parameters
    ! calculated by the subroutine.  Compare that result to the value that was
    ! used to start with.  They should match numerically within a very small
    ! tolerance.
    !
    ! Test 3
    ! Recalculate the value of <w'^3> using the PDF parameters calculated by the
    ! subroutine.  Compare that result to the value that was used to start with.
    ! They should match numerically within a very small tolerance.  This test
    ! can only be run when F_w > 0 or F_w = 0 with Skw = 0.
    !
    ! Test 4
    ! Recalculate the skewness of the w using the PDF parameters calculated by
    ! the subroutine.  Compare that result to the value that was used to start
    ! with.  They should match numerically within a very small tolerance.  This
    ! test can only be run when F_w > 0 or F_w = 0 with Skw = 0.
    !
    ! Test 5
    ! Check that sigma_w_1 >= 0, sigma_w_2 >= 0, and 0 < mixt_frac < 1.
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
        w_tol,         &
        w_tol_sqd,     &
        fstdout,       &
        fstderr

    use new_pdf, only: &
        calc_setter_var_params    ! Procedure(s)

    use mu_sigma_hm_tests, only: &
        produce_seed    ! Procedure(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    real( kind = core_rknd ) :: &
      wm,     & ! Mean of w (overall)                [m/s]
      rtm,    & ! Mean of rt (overall)               [kg/kg]
      thlm,   & ! Mean of thl (overall)              [K]
      wp2,    & ! Variance of w (overall)            [m^2/s^2]
      rtp2,   & ! Variance of rt (overall)           [kg^2/kg^2]
      thlp2,  & ! Variance of thl (overall)          [K^2]
      wp3,    & ! <w'^3>                             [m^3/s^3]
      rtp3,   & ! <rt'^3>                            [kg^3/kg^3]
      thlp3,  & ! <thl'^3>                           [K^3]
      Skw,    & ! Skewness of w (overall)            [-]
      Skrt,   & ! Skewness of rt (overall)           [-]
      Skthl,  & ! Skewness of thl (overall)          [-]
      wprtp,  & ! Covariance of w and rt (overall)   [(m/s)kg/kg]
      wpthlp    ! Covariance of w and thl (overall)  [(m/s)K]

    real( kind = core_rknd ) :: &
      F_w,    & !
      zeta_w, & !
      F_rt,   & !
      F_thl     !

    real( kind = core_rknd ) :: &
      mu_w_1,          & ! Mean of w (1st PDF component)             [m/s]
      mu_w_2,          & ! Mean of w (2nd PDF component)             [m/s]
      mu_rt_1,         & ! Mean of rt (1st PDF component)            [kg/kg]
      mu_rt_2,         & ! Mean of rt (2nd PDF component)            [kg/kg]
      mu_thl_1,        & ! Mean of thl (1st PDF component)           [K]
      mu_thl_2,        & ! Mean of thl (2nd PDF component)           [K]
      sigma_w_1,       & ! Standard deviation of w (1st PDF comp.)   [m/s]
      sigma_w_2,       & ! Standard deviation of w (2nd PDF comp.)   [m/s]
      sigma_w_1_sqd,   & ! Variance of w (1st PDF component)         [m^2/s^2]
      sigma_w_2_sqd,   & ! Variance of w (2nd PDF component)         [m^2/s^2]
      sigma_rt_1_sqd,  & ! Variance of rt (1st PDF component)        [kg^2/kg^2]
      sigma_rt_2_sqd,  & ! Variance of rt (2nd PDF component)        [kg^2/kg^2]
      sigma_thl_1_sqd, & ! Variance of thl (1st PDF component)       [K^2]
      sigma_thl_2_sqd, & ! Variance of thl (2nd PDF component)       [K^2]
      mixt_frac          ! Mixture fraction                          [-]

    real( kind = core_rknd ) :: &
      recalc_wm,    & ! Recalculation of <w> using PDF parameters          [m/s]
      recalc_wp2,   & ! Recalculation of <w'^2> using PDF parameters   [m^2/s^2]
      recalc_wp3,   & ! Recalculation of <w'^3> using PDF parameters   [m^3/s^3]
      recalc_Skw,   & ! Recalculation of Skw using PDF parameters            [-]
      recalc_rtm,   & ! Recalculation of <rt> using PDF parameters       [kg/kg]
      recalc_rtp2,  & ! Recalculation of <rt'^2> using PDF params.   [kg^2/kg^2]
      recalc_rtp3,  & ! Recalculation of <rt'^3> using PDF params.   [kg^3/kg^3]
      recalc_Skrt,  & ! Recalculation of Skrt using PDF parameters           [-]
      recalc_thlm,  & ! Recalculation of <thl> using PDF parameters          [K]
      recalc_thlp2, & ! Recalculation of <thl'^2> using PDF parameters     [K^2]
      recalc_thlp3, & ! Recalculation of <thl'^3> using PDF parameters     [K^3]
      recalc_Skthl    ! Recalculation of Skthl using PDF parameters          [-]

    ! Tiny tolerance for acceptable numerical difference between two results.
    real( kind = core_rknd ), parameter :: &
      tol = 1.0e-8_core_rknd

    logical :: &
      l_pass_test_1, & ! Flag for passing test 1
      l_pass_test_2, & ! Flag for passing test 2
      l_pass_test_3, & ! Flag for passing test 3
      l_pass_test_4, & ! Flag for passing test 4
      l_pass_test_5    ! Flag for passing test 5

    integer :: &
      pdf_parameter_unit_tests  ! Returns pass or fail

    integer :: &
      num_failed_sets    ! Records the number of failed parameter sets

    integer :: &
      seed_size    ! The size of the random seed array expected by the system

    integer, dimension(:), allocatable :: &
      seed_vals    ! Values used to seed the random number generator

    real( kind = core_rknd ) :: &
      rand1,  & ! Random number 1 used for PDF parameter set 10
      rand2,  & ! Random number 2 used for PDF parameter set 10
      rand3,  & ! Random number 3 used for PDF parameter set 10
      rand4,  & ! Random number 4 used for PDF parameter set 10
      rand5,  & ! Random number 5 used for PDF parameter set 10
      rand6,  & ! Random number 6 used for PDF parameter set 10
      rand7,  & ! Random number 7 used for PDF parameter set 10
      rand8,  & ! Random number 8 used for PDF parameter set 10
      rand9,  & ! Random number 9 used for PDF parameter set 10
      rand10, & ! Random number 10 used for PDF parameter set 10
      rand11    ! Random number 11 used for PDF parameter set 10

    integer, parameter :: &
      num_param_sets = 10, & ! Number of different PDF parameter sets used
      num_F_w = 11,        & ! Number of different values of F_w used
      num_zeta_w = 11        ! Number of different values of zeta_w used

    integer :: &
      iter_param_sets, & ! Loop index for PDF parameter set
      iter_F_w,        & ! Loop index for value of F_w
      iter_zeta_w        ! Loop index for value of zeta_w


    write(fstdout,*) ""
    write(fstdout,*) "Performing PDF parameter values unit test"
    write(fstdout,*) "========================================="
    write(fstdout,*) ""

    ! Initialize number of failed parameter sets.
    num_failed_sets = 0

    ! Perform unit tests.
    do iter_param_sets = 1, num_param_sets, 1

       if ( iter_param_sets == 1 ) then
          write(fstdout,*) "PDF parameter set 1:"
          ! SAM LES of RICO: t = 4200 minutes and z = 2780 meters.
          wm = -0.005_core_rknd
          wp2 = 0.0337419_core_rknd
          wp3 = 0.0538501_core_rknd
          Skw = wp3 / wp2**1.5
          rtm = 4.28684e-3_core_rknd
          rtp2 = 1.17865e-6_core_rknd
          Skrt = 1.48809_core_rknd
          rtp3 = Skrt * rtp2**1.5
          thlm = 309.131_core_rknd
          thlp2 = 0.05992_core_rknd
          Skthl = -1.45809_core_rknd
          thlp3 = Skthl * thlp2**1.5
          wprtp = 2.37254e-5_core_rknd
          wpthlp = -5.55127e-3_core_rknd
       elseif ( iter_param_sets == 2 ) then
          write(fstdout,*) "PDF parameter set 2:"
          wm = -0.001_core_rknd
          wp2 = 0.2_core_rknd
          Skw = 4.5_core_rknd
          wp3 = Skw * wp2**1.5
          rtm = 1.0e-2_core_rknd
          rtp2 = 5.0e-6_core_rknd
          Skrt = 2.5_core_rknd
          rtp3 = Skrt * rtp2**1.5
          thlm = 305.0_core_rknd
          thlp2 = 0.1_core_rknd
          Skthl = -2.5_core_rknd
          thlp3 = Skthl * thlp2**1.5
          wprtp = 5.0e-5_core_rknd
          wpthlp = -5.0e-3_core_rknd
       elseif ( iter_param_sets == 3 ) then
          write(fstdout,*) "PDF parameter set 3:"
          wm = -0.001_core_rknd
          wp2 = 0.2_core_rknd
          Skw = 4.5_core_rknd
          wp3 = Skw * wp2**1.5
          rtm = 1.0e-2_core_rknd
          rtp2 = 5.0e-6_core_rknd
          Skrt = 2.5_core_rknd
          rtp3 = Skrt * rtp2**1.5
          thlm = 305.0_core_rknd
          thlp2 = 0.1_core_rknd
          Skthl = -2.5_core_rknd
          thlp3 = Skthl * thlp2**1.5
          wprtp = 5.0e-5_core_rknd
          wpthlp = -5.0e-3_core_rknd
       elseif ( iter_param_sets == 4 ) then
          write(fstdout,*) "PDF parameter set 4:"
          wm = -0.001_core_rknd
          wp2 = 0.2_core_rknd
          Skw = 4.5_core_rknd
          wp3 = Skw * wp2**1.5
          rtm = 1.0e-2_core_rknd
          rtp2 = 5.0e-6_core_rknd
          Skrt = 2.5_core_rknd
          rtp3 = Skrt * rtp2**1.5
          thlm = 305.0_core_rknd
          thlp2 = 0.1_core_rknd
          Skthl = -2.5_core_rknd
          thlp3 = Skthl * thlp2**1.5
          wprtp = 5.0e-5_core_rknd
          wpthlp = -5.0e-3_core_rknd
       elseif ( iter_param_sets == 5 ) then
          write(fstdout,*) "PDF parameter set 5:"
          wm = -0.001_core_rknd
          wp2 = 0.2_core_rknd
          Skw = 4.5_core_rknd
          wp3 = Skw * wp2**1.5
          rtm = 1.0e-2_core_rknd
          rtp2 = 5.0e-6_core_rknd
          Skrt = 2.5_core_rknd
          rtp3 = Skrt * rtp2**1.5
          thlm = 305.0_core_rknd
          thlp2 = 0.1_core_rknd
          Skthl = -2.5_core_rknd
          thlp3 = Skthl * thlp2**1.5
          wprtp = 5.0e-5_core_rknd
          wpthlp = -5.0e-3_core_rknd
       elseif ( iter_param_sets == 6 ) then
          write(fstdout,*) "PDF parameter set 6:"
          wm = -0.001_core_rknd
          wp2 = 0.2_core_rknd
          Skw = 4.5_core_rknd
          wp3 = Skw * wp2**1.5
          rtm = 1.0e-2_core_rknd
          rtp2 = 5.0e-6_core_rknd
          Skrt = 2.5_core_rknd
          rtp3 = Skrt * rtp2**1.5
          thlm = 305.0_core_rknd
          thlp2 = 0.1_core_rknd
          Skthl = -2.5_core_rknd
          thlp3 = Skthl * thlp2**1.5
          wprtp = 5.0e-5_core_rknd
          wpthlp = -5.0e-3_core_rknd
       elseif ( iter_param_sets == 7 ) then
          write(fstdout,*) "PDF parameter set 7:"
          wm = -0.001_core_rknd
          wp2 = 0.2_core_rknd
          Skw = 4.5_core_rknd
          wp3 = Skw * wp2**1.5
          rtm = 1.0e-2_core_rknd
          rtp2 = 5.0e-6_core_rknd
          Skrt = 2.5_core_rknd
          rtp3 = Skrt * rtp2**1.5
          thlm = 305.0_core_rknd
          thlp2 = 0.1_core_rknd
          Skthl = -2.5_core_rknd
          thlp3 = Skthl * thlp2**1.5
          wprtp = 5.0e-5_core_rknd
          wpthlp = -5.0e-3_core_rknd
       elseif ( iter_param_sets == 8 ) then
          write(fstdout,*) "PDF parameter set 8:"
          wm = -0.001_core_rknd
          wp2 = 0.2_core_rknd
          Skw = 4.5_core_rknd
          wp3 = Skw * wp2**1.5
          rtm = 1.0e-2_core_rknd
          rtp2 = 5.0e-6_core_rknd
          Skrt = 2.5_core_rknd
          rtp3 = Skrt * rtp2**1.5
          thlm = 305.0_core_rknd
          thlp2 = 0.1_core_rknd
          Skthl = -2.5_core_rknd
          thlp3 = Skthl * thlp2**1.5
          wprtp = 5.0e-5_core_rknd
          wpthlp = -5.0e-3_core_rknd
       elseif ( iter_param_sets == 9 ) then
          write(fstdout,*) "PDF parameter set 9:"
          wm = -0.001_core_rknd
          wp2 = 0.2_core_rknd
          Skw = 4.5_core_rknd
          wp3 = Skw * wp2**1.5
          rtm = 1.0e-2_core_rknd
          rtp2 = 5.0e-6_core_rknd
          Skrt = 2.5_core_rknd
          rtp3 = Skrt * rtp2**1.5
          thlm = 305.0_core_rknd
          thlp2 = 0.1_core_rknd
          Skthl = -2.5_core_rknd
          thlp3 = Skthl * thlp2**1.5
          wprtp = 5.0e-5_core_rknd
          wpthlp = -5.0e-3_core_rknd
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
          call random_number( rand6 )
          call random_number( rand7 )
          call random_number( rand8 )
          call random_number( rand9 )
          call random_number( rand10 )
          call random_number( rand11 )
          ! The value of wm can range from -0.5 m/s to 0.5 m/s.
          wm = 1.0_core_rknd * rand1 - 0.5_core_rknd
          ! The value of rtm can range from 1.0e-3 kg/kg to 2.1e-2 kg/kg.
          rtm = 2.0e-2_core_rknd * rand2 + 1.0e-3_core_rknd
          ! The value of thlm can range from 290 K to 310 K.
          thlm = 20.0_core_rknd * rand3 + 290.0_core_rknd
          ! The value of wp2 can range from 0 m^2/s^2 to 1.0 m^2/s^2.
          wp2 = 1.0_core_rknd * rand4
          ! The value of rtp2 can range from 0 kg^2/kg^2 to 0.25 * rtm^2
          rtp2 = 0.25_core_rknd * rtm**2 * rand5
          ! The value of thlp2 can range from 0 K^2 to 3.0 K^2.
          wp2 = 3.0_core_rknd * rand6
          ! The value of Skw can range from -4.5 to 4.5.
          Skw = 9.0_core_rknd * rand7 - 4.5_core_rknd
          ! Calculate wp3.
          wp3 = Skw * wp2**1.5
          ! The value of Skrt can range from -4.5 to 4.5.
          Skrt = 9.0_core_rknd * rand8 - 4.5_core_rknd
          ! Calculate rtp3.
          rtp3 = Skrt * rtp2**1.5
          ! The value of Skthl can range from -4.5 to 4.5.
          Skthl = 9.0_core_rknd * rand9 - 4.5_core_rknd
          ! Calculate thlp3.
          thlp3 = Skthl * thlp2**1.5
          ! Use a random number to calculate the value of wprtp.
          wprtp = rand10 * sqrt( wp2 ) * sqrt( rtp2 ) * sign( one, Skrt )
          ! Use a random number to calculate the value of wpthlp.
          wpthlp = rand11 * sqrt( wp2 ) * sqrt( thlp2 ) * sign( one, Skthl )
       endif ! iter_param_sets == index

       ! Print PDF parameters
       write(fstdout,*) "wm = ", wm
       write(fstdout,*) "rtm = ", rtm
       write(fstdout,*) "thlm = ", thlm
       write(fstdout,*) "wp2 = ", wp2
       write(fstdout,*) "rtp2 = ", rtp2
       write(fstdout,*) "thlp2 = ", thlp2
       write(fstdout,*) "Skw = ", Skw
       write(fstdout,*) "wp3 = ", wp3
       write(fstdout,*) "Skrt = ", Skrt
       write(fstdout,*) "rtp3 = ", rtp3
       write(fstdout,*) "Skthl = ", Skthl
       write(fstdout,*) "thlp3 = ", thlp3
       write(fstdout,*) "wprtp = ", wprtp
       write(fstdout,*) "wpthlp = ", wpthlp
       write(fstdout,*) ""

       write(fstdout,*) "Running tests for the above parameter set for all " &
                        // "combinations of F_w and zeta_w.  F_w values are " &
                        // "0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, " &
                        // "and 1.  Zeta_w values are 0, -1/2, 1, -1/3, 1/2, " &
                        // "-1/5, 1/4, -9/10, 9, -3/4, and 3."
       write(fstdout,*) ""

       do iter_F_w = 1, num_F_w, 1

          ! Set the value of F_w.
          F_w = max( 0.1_core_rknd * real( iter_F_w - 1, kind = core_rknd ), &
                     1.0e-5_core_rknd )

          do iter_zeta_w = 1, num_zeta_w, 1

             ! Set the value of zeta_w.
             if ( iter_zeta_w == 1 ) then
                zeta_w = zero
             elseif ( iter_zeta_w == 2 ) then
                zeta_w = -one_half
             elseif ( iter_zeta_w == 3 ) then
                zeta_w = one
             elseif ( iter_zeta_w == 4 ) then
                zeta_w = -one_third
             elseif ( iter_zeta_w == 5 ) then
                zeta_w = one_half
             elseif ( iter_zeta_w == 6 ) then
                zeta_w = -0.2_core_rknd
             elseif ( iter_zeta_w == 7 ) then
                zeta_w = one_fourth
             elseif ( iter_zeta_w == 8 ) then
                zeta_w = -0.9_core_rknd
             elseif ( iter_zeta_w == 9 ) then
                zeta_w = 9.0_core_rknd
             elseif ( iter_zeta_w == 10 ) then
                zeta_w = -three_fourths
             elseif ( iter_zeta_w == 11 ) then
                zeta_w = three
             endif

             ! Call the subroutine for calculating mu_w_1, mu_w_2, sigma_w_1,
             ! sigma_w_2, and mixt_frac.
             call calc_setter_var_params( wm, wp2, Skw, F_w, zeta_w, & ! In
                                          mu_w_1, mu_w_2, sigma_w_1, & ! Out
                                          sigma_w_2, mixt_frac )       ! Out

             ! Recalculate the values of <wm>, <w'^2>, <w'^3>, and Skw using the
             ! PDF parameters that are output from the subroutine.
             call recalc_single_var( wm, mu_w_1, mu_w_2,              & ! In
                                     sigma_w_1, sigma_w_2, mixt_frac, & ! In
                                     recalc_wm, recalc_wp2,           & ! Out
                                     recalc_wp3, recalc_Skw           ) ! Out

             ! Test 1
             ! Compare the original wm to its recalculated value.
             !    | ( <w>|_recalc - <w> ) / <w> |  <=  tol;
             ! which can be rewritten as:
             !    | <w>|_recalc - <w> |  <=  | <w> | * tol.
             if ( abs( recalc_wm - wm ) <= max( abs( wm ), w_tol ) * tol ) then
                l_pass_test_1 = .true.
             else
                l_pass_test_1 = .false.
                write(fstderr,*) "Test 1 failed"
                write(fstderr,*) "wm = ", wm
                write(fstderr,*) "recalc_wm = ", recalc_wm
                write(fstderr,*) ""
             endif

             ! Test 2
             ! Compare the original wp2 to its recalculated value.
             !    | ( <w'^2>|_recalc - <w'^2> ) / <w'^2> |  <=  tol;
             ! which can be rewritten as:
             !    | <w'^2>|_recalc - <w'^2> |  <=  | <w'^2> | * tol;
             ! and since <w'^2> is always positive:
             !    | <w'^2>|_recalc - <w'^2> |  <=  <w'^2> * tol.
             if ( abs( recalc_wp2 - wp2 ) <= max( wp2, w_tol_sqd ) * tol ) then
                l_pass_test_2 = .true.
             else
                l_pass_test_2 = .false.
                write(fstderr,*) "Test 2 failed"
                write(fstderr,*) "wp2 = ", wp2
                write(fstderr,*) "recalc_wp2 = ", recalc_wp2
                write(fstderr,*) ""
             endif

             ! Test 3
             ! Compare the original wp3 to its recalculated value.
             !    | ( <w'^3>|_recalc - <w'^3> ) / <w'^3> |  <=  tol;
             ! which can be rewritten as:
             !    | <w'^3>|_recalc - <w'^3> |  <=  | <w'^3> | * tol.
             if ( F_w > zero ) then
                if ( abs( recalc_wp3 - wp3 ) &
                     <= max( abs( wp3 ), w_tol**3 ) * tol ) then
                   l_pass_test_3 = .true.
                else
                   l_pass_test_3 = .false.
                   write(fstderr,*) "Test 3 failed"
                   write(fstderr,*) "wp3 = ", wp3
                   write(fstderr,*) "recalc_wp3 = ", recalc_wp3
                   write(fstderr,*) ""
                endif
             else ! F_w = 0
                l_pass_test_3 = .true.
             endif ! F_w > 0

             ! Test 4
             ! Compare the original Skw to its recalculated value.
             !    | ( Skw|_recalc - Skw ) / Skw |  <=  tol;
             ! which can be rewritten as:
             !    | Skw|_recalc - Skw |  <=  | Skw | * tol.
             if ( F_w > zero ) then
                if ( abs( recalc_Skw - Skw ) &
                     <= max( abs( Skw ), 1.0e-3_core_rknd ) * tol ) then
                   l_pass_test_4 = .true.
                else
                   l_pass_test_4 = .false.
                   write(fstderr,*) "Test 4 failed"
                   write(fstderr,*) "Skw = ", Skw
                   write(fstderr,*) "recalc_Skw = ", recalc_Skw
                   write(fstderr,*) ""
                endif
             else ! F_w = 0
                l_pass_test_4 = .true.
             endif ! F_w > 0

             ! Test 5
             ! Check that sigma_w_1 and sigma_w_2 have a value of at least zero,
             ! and that mixture fraction has a value between 0 and 1.
             if ( ( sigma_w_1 >= zero ) .and. ( sigma_w_2 >= zero ) &
                  .and. ( mixt_frac > zero ) .and. ( mixt_frac < one ) ) then
                l_pass_test_5 = .true.
             else
                l_pass_test_5 = .false.
                write(fstderr,*) "Test 5 failed"
                write(fstderr,*) "sigma_w_1 = ", sigma_w_1
                write(fstderr,*) "sigma_w_2 = ", sigma_w_2
                write(fstderr,*) "mixt_frac = ", mixt_frac
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
                write(fstderr,*) "F_w = ", F_w
                write(fstderr,*) "zeta_w = ", zeta_w
                write(fstderr,*) ""
             endif


          enddo ! iter_zeta_w = 1, num_zeta_w, 1

       enddo ! iter_F_w = 1, num_F_w, 1

    enddo ! iter_param_sets = 1, num_param_sets, 1


    ! Print results and return exit code.
    if ( num_failed_sets == 0 ) then
       write(fstdout,'(1x,A)') "Success!"
       write(fstdout,*) ""
       pdf_parameter_unit_tests = 0 ! Exit Code = 0, Success!
    else ! num_failed_sets > 0
       write(fstdout,'(1x,A,I4,A)') "There were ", num_failed_sets, &
                                    " failed parameter sets."
       write(fstdout,*) ""
       pdf_parameter_unit_tests = 1 ! Exit Code = 1, Fail
    endif ! num_failed_sets = 0


    return

  end function pdf_parameter_unit_tests

  !=============================================================================
  subroutine recalc_single_var( xm, mu_x_1, mu_x_2,              & ! In
                                sigma_x_1, sigma_x_2, mixt_frac, & ! In
                                recalc_xm, recalc_xp2,           & ! Out
                                recalc_xp3, recalc_Skx           ) ! Out

    use constants_clubb, only: &
        three, & ! Variable(s)
        one

    use clubb_precision, only: &
        core_rknd

    implicit none

    real( kind = core_rknd ), intent(in) :: &
      xm,        & ! Mean of x (overall)                  [units vary]
      mu_x_1,    & ! Mean of x (1st PDF component)        [units vary]
      mu_x_2,    & ! Mean of x (2nd PDF component)        [units vary]
      sigma_x_1, & ! Variance of x (1st PDF component)    [units vary]
      sigma_x_2, & ! Variance of x (2nd PDF component)    [units vary]
      mixt_frac    ! Mixture fraction                     [-]

    real( kind = core_rknd ), intent(out) :: &
      recalc_xm,  & ! Recalculation of <x> using PDF parameters    [units vary]
      recalc_xp2, & ! Recalculation of <x'^2> using PDF parameters [(un vary)^2]
      recalc_xp3, & ! Recalculation of <x'^3> using PDF parameters [(un vary)^3]
      recalc_Skx    ! Recalculation of Skx using PDF parameters    [-]


    recalc_xm = mixt_frac * mu_x_1 + ( one - mixt_frac ) * mu_x_2

    recalc_xp2 &
    = mixt_frac * ( ( mu_x_1 - xm )**2 + sigma_x_1**2 ) &
      + ( one - mixt_frac ) * ( ( mu_x_2 - xm )**2 + sigma_x_2**2 )

    recalc_xp3 &
    = mixt_frac * ( mu_x_1 - xm ) &
      * ( ( mu_x_1 - xm )**2 + three * sigma_x_1**2 ) &
      + ( one - mixt_frac ) * ( mu_x_2 - xm ) &
        * ( ( mu_x_2 - xm )**2 + three * sigma_x_2**2 )

    recalc_Skx = recalc_xp3 / recalc_xp2**1.5


    return

  end subroutine recalc_single_var

  !=============================================================================

end module pdf_parameter_tests
