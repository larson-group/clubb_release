! $Id$
!===============================================================================
module pdf_parameter_tests

  implicit none

  public :: pdf_parameter_unit_tests    ! Procedure(s)

  private :: setter_var_tests,  & ! Procedure(s)
             recalc_single_var

  private  ! default scope

  contains

  !=============================================================================
  function pdf_parameter_unit_tests( gr, test_PDF_type )

    ! Description:
    ! Unit testing framework for the code that calculates the mixture fraction,
    ! the PDF component means, and the PDF component standard deviations for the
    ! trivariate PDF of w, rt, and theta-l.
    !
    !
    ! Description of tests for the new PDF:
    !
    ! There are 10 different input parameter sets specified.  The 10th input
    ! parameter set is a randomized set.  There are two different sets of tests.
    ! The first set of tests is used for the "setting" variable (the variable
    ! that is used to set the mixture fraction).  The second set of tests is
    ! for the full PDF.
    !
    ! In the first set of tests (for only the setting variable), for each PDF
    ! parameter set, 11 different values of F_w are used (0, 0.1, 0.2, 0.3,
    ! 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0).  Furthermore, for every value of
    ! F_w, 11 different values of zeta_w are used (0, -1/2, 1, -1/3, 1/2,
    ! -1/5, 1/4, -9/10, 9, -3/4, and 3).  So, for every PDF parameter set, 121
    ! combinations of F_w and zeta_w are used.
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
    ! Check that sigma_w_1^2 >= 0, sigma_w_2^2 >= 0, and 0 < mixt_frac < 1.
    !
    ! Test 6
    ! Check that mu_w_1 >= mu_w_2.
    !
    ! Test 7
    ! Calculate <w'^4> as found by integrating over the PDF.  Then, calculate
    ! <w'^4> according to the equation used to handle its calculation implicitly
    ! as part of the <w'^3> turbulent advection term, which is of the form:
    ! <w'^4> = coef_wp4_implicit * <w'^2>^2.  The two calculations of <w'^4>
    ! should match numerically within a very small tolerance.
    !
    ! In the second set of tests, the full PDF is called for each of the 10
    ! aforementioned input parameter sets.  The values of parameters F and zeta
    ! are set internally.  Each subroutine call produces values of mu_w_1,
    ! mu_w_2, mu_rt_1, mu_rt_2, mu_thl_1, mu_thl_2, sigma_w_1^2, sigma_w_2^2,
    ! sigma_rt_1^2, sigma_rt_2^2, sigma_thl_1^2, sigma_thl_2^2, and mixt_frac.
    ! These values are evaluated by 13 tests.
    !
    ! Test 8, Test 12, and Test 16 are analogous to Test 1 for w, rt, and
    ! theta-l respectively.  Test 9, Test 13, and Test 17 are analogous to
    ! Test 2 for w, rt, and theta-l respectively.  Test 10, Test 14, and Test 18
    ! are analogous to Test 3 for w, rt, and theta-l, respectively.  Test 11,
    ! Test 15, and Test 19 are analogous to Test 4 for w, rt, and theta-l,
    ! respectively.  Test 20 is analogous to Test 5, but also checks that
    ! sigma_rt_1 >= 0, sigma_rt_2 >= 0, sigma_thl_1 >= 0, and sigma_thl_2 >= 0.
    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: grid ! Type

    use constants_clubb, only: &
        three,         & ! Constant(s)
        two,           &
        one,           &
        three_fourths, &
        two_thirds,    &
        one_half,      &
        one_third,     &
        one_fourth,    &
        zero,          &
        w_tol,         &
        w_tol_sqd,     &
        rt_tol,        &
        thl_tol,       &
        fstdout,       &
        fstderr

    use new_pdf, only: &
        calc_setter_var_params, & ! Procedure(s)
        calc_coef_wp4_implicit

    use new_pdf_main, only: &
        new_pdf_driver    ! Procedure(s)

    use new_hybrid_pdf, only: &
        calculate_w_params,          & ! Procedure(s)
        calculate_coef_wp4_implicit

    use new_hybrid_pdf_main, only: &
        new_hybrid_pdf_driver    ! Procedure(s)

    use adg1_adg2_3d_luhar_pdf, only: &
        ADG1_w_closure,  & ! Procedure(s)
        ADG1_pdf_driver

    use sigma_sqd_w_module, only: &
        compute_sigma_sqd_w    ! Procedure(s)

    use new_tsdadg_pdf, only: &
        calc_setter_parameters, & ! Procedure(s) 
        calc_L_x_Skx_fnc

    use LY93_pdf, only: &
        calc_mixt_frac_LY93, & ! Procedure(s)
        calc_params_LY93,    &
        LY93_driver

    use pdf_closure_module, only: &
        calc_wp4_pdf        ! Procedure(s)

    use pdf_parameter_module, only: &
        implicit_coefs_terms    ! Variable Type

    use model_flags, only: &
        set_default_clubb_config_flags, & ! Procedure(s)
        iiPDF_new,        & ! Variable(s)
        iiPDF_ADG1,       &
        iiPDF_TSDADG,     &
        iiPDF_LY93,       &
        iiPDF_new_hybrid, &
        l_gamma_Skw

    use parameters_model, only: &
        sclr_dim    ! Variable(s)

    use parameters_tunable, only: &
        gamma_coef,  & ! Variable(s)
        gamma_coefb, &
        gamma_coefc

    use mu_sigma_hm_tests, only: &
        produce_seed    ! Procedure(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(inout) :: gr

    ! Input Variable
    integer, intent(in) :: &
      test_pdf_type    ! The PDF type being tested

    ! Return Variable
    integer :: &
      pdf_parameter_unit_tests  ! Returns pass or fail

    ! Local Variables
    integer, parameter :: &
      nz = 1

    real( kind = core_rknd ), dimension(nz) :: &
      wm,      & ! Mean of w (overall)                 [m/s]
      rtm,     & ! Mean of rt (overall)                [kg/kg]
      thlm,    & ! Mean of thl (overall)               [K]
      um,      & ! Mean of eastward wind (overall)     [m/s]
      vm,      & ! Mean of northward wind (overall)    [m/s]
      wp2,     & ! Variance of w (overall)             [m^2/s^2]
      rtp2,    & ! Variance of rt (overall)            [kg^2/kg^2]
      thlp2,   & ! Variance of thl (overall)           [K^2]
      up2,     & ! Variance of u (overall)             [m^2/s^2]
      vp2,     & ! Variance of v (overall)             [m^2/s^2]
      wp3,     & ! <w'^3>                              [m^3/s^3]
      rtp3,    & ! <rt'^3>                             [kg^3/kg^3]
      thlp3,   & ! <thl'^3>                            [K^3]
      Skw,     & ! Skewness of w (overall)             [-]
      Skrt,    & ! Skewness of rt (overall)            [-]
      Skthl,   & ! Skewness of thl (overall)           [-]
      Sku,     & ! Skewness of u (overall)             [-]
      Skv,     & ! Skewness of v (overall)             [-]
      wprtp,   & ! Covariance of w and rt (overall)    [(m/s)kg/kg]
      wpthlp,  & ! Covariance of w and thl (overall)   [(m/s)K]
      upwp,    & ! Covariance of u and w (overall)     [(m/s)^2]
      vpwp,    & ! Covariance of v and w (overall)     [(m/s)^2]
      rtpthlp    ! Covariance of rt and thl (overall)  [(kg/kg)K]

    real( kind = core_rknd ), dimension(nz) :: &
      F_w,    & ! Parameter for the spread of the PDF component means of w   [-]
      zeta_w, & ! Parameter for the PDF component variances of w             [-]
      F_rt,   & ! Parameter for the spread of the PDF component means of rt  [-]
      F_thl     ! Parameter for the spread of the PDF component means of thl [-]

    real( kind = core_rknd ), dimension(nz) :: &
      min_F_w,   & ! Minimum allowable value of parameter F_w      [-]
      max_F_w,   & ! Maximum allowable value of parameter F_w      [-]
      min_F_rt,  & ! Minimum allowable value of parameter F_rt     [-]
      max_F_rt,  & ! Maximum allowable value of parameter F_rt     [-]
      min_F_thl, & ! Minimum allowable value of parameter F_thl    [-]
      max_F_thl    ! Maximum allowable value of parameter F_thl    [-]

    real( kind = core_rknd ), dimension(nz) :: &
      sgn_wp2     ! Sign of the variance of w (overall); always positive [-]

    real( kind = core_rknd ), dimension(nz) :: &
      mu_w_1,          & ! Mean of w (1st PDF component)             [m/s]
      mu_w_2,          & ! Mean of w (2nd PDF component)             [m/s]
      mu_rt_1,         & ! Mean of rt (1st PDF component)            [kg/kg]
      mu_rt_2,         & ! Mean of rt (2nd PDF component)            [kg/kg]
      mu_thl_1,        & ! Mean of thl (1st PDF component)           [K]
      mu_thl_2,        & ! Mean of thl (2nd PDF component)           [K]
      mu_u_1,          & ! Mean of u (1st PDF component)             [m/s]
      mu_u_2,          & ! Mean of u (2nd PDF component)             [m/s]
      mu_v_1,          & ! Mean of v (1st PDF component)             [m/s]
      mu_v_2,          & ! Mean of v (2nd PDF component)             [m/s]
      sigma_w_1,       & ! Standard deviation of w (1st PDF comp.)   [m/s]
      sigma_w_2,       & ! Standard deviation of w (2nd PDF comp.)   [m/s]
      sigma_rt_1,      & ! Standard deviation of rt (1st PDF comp.)  [kg/kg]
      sigma_rt_2,      & ! Standard deviation of rt (2nd PDF comp.)  [kg/kg]
      sigma_thl_1,     & ! Standard deviation of thl (1st PDF comp.) [K]
      sigma_thl_2,     & ! Standard deviation of thl (2nd PDF comp.) [K]
      sigma_w_1_sqd,   & ! Variance of w (1st PDF component)         [m^2/s^2]
      sigma_w_2_sqd,   & ! Variance of w (2nd PDF component)         [m^2/s^2]
      sigma_rt_1_sqd,  & ! Variance of rt (1st PDF component)        [kg^2/kg^2]
      sigma_rt_2_sqd,  & ! Variance of rt (2nd PDF component)        [kg^2/kg^2]
      sigma_thl_1_sqd, & ! Variance of thl (1st PDF component)       [K^2]
      sigma_thl_2_sqd, & ! Variance of thl (2nd PDF component)       [K^2]
      sigma_u_1_sqd,   & ! Variance of u (1st PDF component)         [m^2/s^2]
      sigma_u_2_sqd,   & ! Variance of u (2nd PDF component)         [m^2/s^2]
      sigma_v_1_sqd,   & ! Variance of v (1st PDF component)         [m^2/s^2]
      sigma_v_2_sqd,   & ! Variance of v (2nd PDF component)         [m^2/s^2]
      mixt_frac          ! Mixture fraction                          [-]

    real( kind = core_rknd ), dimension(nz) :: &
      coef_sigma_w_1_sqd, & ! sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>      [-]
      coef_sigma_w_2_sqd    ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>      [-]

    real( kind = core_rknd ), dimension(nz) :: &
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

    real( kind = core_rknd ), dimension(nz) :: &
      coef_wp4_implicit, & ! <w'^4> = coef_wp4_implicit * <w'^2>^2           [-]
      wp4_implicit_calc, & ! <w'^4> calculated by coef_wp4_implicit eq [m^4/s^4]
      wp4_pdf_calc         ! <w'^4> calculated by PDF                  [m^4/s^4]

    type(implicit_coefs_terms) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    ! Tiny tolerance for acceptable numerical difference between two results.
    real( kind = core_rknd ), parameter :: &
      tol = 1.0e-8_core_rknd

    logical, dimension(nz) :: &
      l_pass_test_1,  & ! Flag for passing test 1
      l_pass_test_2,  & ! Flag for passing test 2
      l_pass_test_3,  & ! Flag for passing test 3
      l_pass_test_4,  & ! Flag for passing test 4
      l_pass_test_5,  & ! Flag for passing test 5
      l_pass_test_6,  & ! Flag for passing test 6
      l_pass_test_7,  & ! Flag for passing test 7
      l_pass_test_8,  & ! Flag for passing test 8
      l_pass_test_9,  & ! Flag for passing test 9
      l_pass_test_10, & ! Flag for passing test 10
      l_pass_test_11, & ! Flag for passing test 11
      l_pass_test_12, & ! Flag for passing test 12
      l_pass_test_13, & ! Flag for passing test 13
      l_pass_test_14, & ! Flag for passing test 14
      l_pass_test_15, & ! Flag for passing test 15
      l_pass_test_16, & ! Flag for passing test 16
      l_pass_test_17, & ! Flag for passing test 17
      l_pass_test_18, & ! Flag for passing test 18
      l_pass_test_19, & ! Flag for passing test 19
      l_pass_test_20    ! Flag for passing test 20

    integer :: &
      total_num_failed_sets    ! Records total number of failed parameter sets

    integer, dimension(nz) :: &
      num_failed_sets    ! Records the number of failed parameter sets

    logical, dimension(nz) :: &
      l_failed_sets    ! Flag recording when there are failed sets

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

    logical :: &
      l_check_mu_w_1_gte_mu_w_2    ! Flag to check whether mu_w_1 >= mu_w_2

    integer, parameter :: &
      num_param_sets = 10, & ! Number of different PDF parameter sets used
      num_F_w = 11,        & ! Number of different values of F_w used
      num_zeta_w = 11        ! Number of different values of zeta_w used

    integer :: &
      iter_param_sets, & ! Loop index for PDF parameter set
      iter_F_w,        & ! Loop index for value of F_w
      iter_zeta_w        ! Loop index for value of zeta_w

    ! Variables for ADG1
    real( kind = core_rknd ), dimension(nz) :: &
      sqrt_wp2,      & ! Square root of w (ADG1)                         [m/s]
      w_1_n,         & ! Normalized mean of w (1st PDF comp.) (ADG1)     [-]
      w_2_n,         & ! Normalized mean of w (2nd PDF comp.) (ADG1)     [-]
      alpha_thl,     & ! Factor relating to normalized variance for th_l [-]
      alpha_rt,      & ! Factor relating to normalized variance for r_t  [-]
      alpha_u,       & ! Factor relating to normalized variance for u    [-]
      alpha_v,       & ! Factor relating to normalized variance for v    [-]
      gamma_Skw_fnc, & ! Skewness function for tunable parameter gamma   [-]
      sigma_sqd_w      ! Width of individual w plumes (ADG1)             [-]

    real( kind = core_rknd ) :: &
      mixt_frac_max_mag    ! Maximum magnitude of mixture fraction (ADG1)    [-]

    integer, parameter :: &
      num_sigma_sqd_w = 11  ! Number of diff. values of sigma_sqd_w used (ADG1)

    integer :: &
      iter_sigma_sqd_w    ! Loop index for value of sigma_sqd_w (ADG1)

    ! Variables for TSDADG
    real( kind = core_rknd ), dimension(nz) :: &
      big_L_w_1, & ! Parameter for the spread of the 1st PDF comp. mean of w [-]
      big_L_w_2    ! Parameter for the spread of the 2nd PDF comp. mean of w [-]

    real( kind = core_rknd ) :: &
      small_l_w_1, & ! Param. for the spread of the 1st PDF comp. mean of w  [-]
      small_l_w_2    ! Param. for the spread of the 2nd PDF comp. mean of w  [-]

    integer, parameter :: &
      num_small_l_w_1 = 5,  & ! Number of different values of small_l_w_1 used
      num_small_l_w_2 = 11    ! Number of different values of small_l_w_2 used

    integer :: &
      iter_small_l_w_1, & ! Loop index for value of small_l_w_1
      iter_small_l_w_2    ! Loop index for value of small_l_w_2

    ! Scalar variables
    real( kind = core_rknd ), dimension(nz,sclr_dim) :: &
      sclrm,            & ! Mean of passive scalar (overall)        [units vary]
      sclrp2,           & ! Variance of pass. scalar (overall)  [(units vary)^2]
      wpsclrp,          & ! Covariance of w and pass. scalar  [m/s (units vary)]
      Sksclr,           & ! Skewness of sclr (overall)                       [-]
      mu_sclr_1,        & ! Mean of passive scalar (1st PDF comp.)  [units vary]
      mu_sclr_2,        & ! Mean of passive scalar (2nd PDF comp.)  [units vary]
      sigma_sclr_1_sqd, & ! Variance pass. sclr (1st PDF comp.) [(units vary)^2]
      sigma_sclr_2_sqd, & ! Variance pass. sclr (2nd PDF comp.) [(units vary)^2]
      alpha_sclr          ! Factor relating to normalized variance for sclr  [-]

    logical, parameter :: &
      l_scalar_calc = .false. ! Flag to perform calculations for passive scalars

    integer :: idx    ! Loop index
  
    integer :: &
      iiPDF_type,          & ! Selected option for the two-component normal
                             ! (double Gaussian) PDF type to use for the w, rt,
                             ! and theta-l (or w, chi, and eta) portion of
                             ! CLUBB's multivariate, two-component PDF.
      ipdf_call_placement    ! Selected option for the placement of the call to
                             ! CLUBB's PDF.

    logical :: &
      l_use_precip_frac,            & ! Flag to use precipitation fraction in KK microphysics. The
                                      ! precipitation fraction is automatically set to 1 when this
                                      ! flag is turned off.
      l_predict_upwp_vpwp,          & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v>
                                      ! alongside the advancement of <rt>, <w'rt'>, <thl>,
                                      ! <wpthlp>, <sclr>, and <w'sclr'> in subroutine
                                      ! advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still
                                      ! approximated by eddy diffusivity when <u> and <v> are
                                      ! advanced in subroutine advance_windm_edsclrm.
      l_min_wp2_from_corr_wx,       & ! Flag to base the threshold minimum value of wp2 on keeping
                                      ! the overall correlation of w and x (w and rt, as well as w
                                      ! and theta-l) within the limits of -max_mag_correlation_flux
                                      ! to max_mag_correlation_flux.
      l_min_xp2_from_corr_wx,       & ! Flag to base the threshold minimum value of xp2 (rtp2 and
                                      ! thlp2) on keeping the overall correlation of w and x within
                                      ! the limits of -max_mag_correlation_flux to
                                      ! max_mag_correlation_flux.
      l_C2_cloud_frac,              & ! Flag to use cloud fraction to adjust the value of the
                                      ! turbulent dissipation coefficient, C2.
      l_diffuse_rtm_and_thlm,       & ! Diffuses rtm and thlm
      l_stability_correct_Kh_N2_zm, & ! Divides Kh_N2_zm by a stability factor
      l_calc_thlp2_rad,             & ! Include the contribution of radiation to thlp2
      l_upwind_wpxp_ta,             & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered
                                      ! differencing for turbulent or mean advection terms. It
                                      ! affects wprtp, wpthlp, & wpsclrp.
      l_upwind_xpyp_ta,             & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered
                                      ! differencing for turbulent or mean advection terms. It
                                      ! affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp, &
                                      ! sclrpthlp.
      l_upwind_xm_ma,               & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered
                                      ! differencing for turbulent or mean advection terms. It
                                      ! affects rtm, thlm, sclrm, um and vm.
      l_uv_nudge,                   & ! For wind speed nudging.
      l_rtm_nudge,                  & ! For rtm nudging
      l_tke_aniso,                  & ! For anisotropic turbulent kinetic energy, i.e.
                                      ! TKE = 1/2 (u'^2 + v'^2 + w'^2)
      l_vert_avg_closure,           & ! Use 2 calls to pdf_closure and the trapezoidal rule to
                                      ! compute the varibles that are output from high order
                                      ! closure
      l_trapezoidal_rule_zt,        & ! If true, the trapezoidal rule is called for the
                                      ! thermodynamic-level variables output from pdf_closure.
      l_trapezoidal_rule_zm,        & ! If true, the trapezoidal rule is called for three
                                      ! momentum-level variables - wpthvp, thlpthvp, and rtpthvp -
                                      ! output from pdf_closure.
      l_call_pdf_closure_twice,     & ! This logical flag determines whether or not to call
                                      ! subroutine pdf_closure twice.  If true, pdf_closure is
                                      ! called first on thermodynamic levels and then on momentum
                                      ! levels so that each variable is computed on its native
                                      ! level.  If false, pdf_closure is only called on
                                      ! thermodynamic levels, and variables which belong on
                                      ! momentum levels are interpolated.
      l_standard_term_ta,           & ! Use the standard discretization for the turbulent advection
                                      ! terms.  Setting to .false. means that a_1 and a_3 are
                                      ! pulled outside of the derivative in
                                      ! advance_wp2_wp3_module.F90 and in
                                      ! advance_xp2_xpyp_module.F90.
      l_partial_upwind_wp3,         & ! Flag to use an "upwind" discretization rather
                                      ! than a centered discretization for the portion
                                      ! of the wp3 turbulent advection term for ADG1
                                      ! that is linearized in terms of wp3<t+1>.
                                      ! (Requires ADG1 PDF and l_standard_term_ta).
      l_godunov_upwind_wpxp_ta,     & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered
                                      ! differencing for turbulent advection terms.
                                      ! It affects  wpxp only.
      l_godunov_upwind_xpyp_ta,     & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered 
                                      ! differencing for turbulent advection terms. It affects
                                      ! xpyp only.
      l_bc_at_constant_height,      & ! Flag for having CLUBB calculate boundary conditions at 
                                      ! a constant height level
      l_linear_Kh_dp_term,          & ! This flag detrmines whether we ignore the part of dp 
                                      ! term that is related to dKh/dz 
      l_mono_cubic_sounding,        & ! This flag determines whether we want to use the mono cubic
                                      ! spline interpolartion instead of linear interpolation to 
                                      ! map the sounding profile to the clubb grid 
      l_use_cloud_cover,            & ! Use cloud_cover and rcm_in_layer to help boost cloud_frac
                                      ! and rcm to help increase cloudiness at coarser grid
                                      ! resolutions.
      l_diagnose_correlations,      & ! Diagnose correlations instead of using fixed ones
      l_calc_w_corr,                & ! Calculate the correlations between w and the hydrometeors
      l_const_Nc_in_cloud,          & ! Use a constant cloud droplet conc. within cloud (K&K)
      l_fix_w_chi_eta_correlations, & ! Use a fixed correlation for s and t Mellor(chi/eta)
      l_stability_correct_tau_zm,   & ! Use tau_N2_zm instead of tau_zm in wpxp_pr1 stability
                                      ! correction
      l_damp_wp2_using_em,          & ! In wp2 equation, use a dissipation formula of
                                      ! -(2/3)*em/tau_zm, as in Bougeault (1981)
      l_do_expldiff_rtm_thlm,       & ! Diffuse rtm and thlm explicitly
      l_Lscale_plume_centered,      & ! Alternate that uses the PDF to compute the perturbed values
      l_diag_Lscale_from_tau,       & ! First diagnose dissipation time tau, and then diagnose the
                                      ! mixing length scale as Lscale = tau * tke
      l_use_C7_Richardson,          & ! Parameterize C7 based on Richardson number
      l_use_C11_Richardson,         & ! Parameterize C11 and C16 based on Richardson number
      l_use_shear_Richardson,       & ! Use shear in the calculation of Richardson number
      l_brunt_vaisala_freq_moist,   & ! Use a different formula for the Brunt-Vaisala frequency in
                                      ! saturated atmospheres (from Durran and Klemp, 1982)
      l_use_thvm_in_bv_freq,        & ! Use thvm in the calculation of Brunt-Vaisala frequency
      l_rcm_supersat_adj,           & ! Add excess supersaturated vapor to cloud water
      l_single_C2_Skw,              & ! Use a single Skewness dependent C2 for rtp2, thlp2, and
                                      ! rtpthlp
      l_damp_wp3_Skw_squared,       & ! Set damping on wp3 to use Skw^2 rather than Skw^4
      l_prescribed_avg_deltaz,      & ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz
      l_update_pressure,            & ! Flag for having CLUBB update pressure and exner
      l_lmm_stepping,               & ! Apply Linear Multistep Method (LMM) Stepping
      l_e3sm_config,                & ! Run model with E3SM settings
      l_use_tke_in_wp3_pr_turb_term   ! Use TKE formulation for wp3 pr_turb term

    call set_default_clubb_config_flags( iiPDF_type, &
                                         ipdf_call_placement, &
                                         l_use_precip_frac, &
                                         l_predict_upwp_vpwp, &
                                         l_min_wp2_from_corr_wx, &
                                         l_min_xp2_from_corr_wx, &
                                         l_C2_cloud_frac, &
                                         l_diffuse_rtm_and_thlm, &
                                         l_stability_correct_Kh_N2_zm, &
                                         l_calc_thlp2_rad, &
                                         l_upwind_wpxp_ta, &
                                         l_upwind_xpyp_ta, &
                                         l_upwind_xm_ma, &
                                         l_uv_nudge, &
                                         l_rtm_nudge, &
                                         l_tke_aniso, &
                                         l_vert_avg_closure, &
                                         l_trapezoidal_rule_zt, &
                                         l_trapezoidal_rule_zm, &
                                         l_call_pdf_closure_twice, &
                                         l_standard_term_ta, &
                                         l_partial_upwind_wp3, &
                                         l_godunov_upwind_wpxp_ta, &
                                         l_godunov_upwind_xpyp_ta, &
                                         l_bc_at_constant_height, & 
                                         l_linear_Kh_dp_term, & 
                                         l_mono_cubic_sounding, & 
                                         l_use_cloud_cover, &
                                         l_diagnose_correlations, &
                                         l_calc_w_corr, &
                                         l_const_Nc_in_cloud, &
                                         l_fix_w_chi_eta_correlations, &
                                         l_stability_correct_tau_zm, &
                                         l_damp_wp2_using_em, &
                                         l_do_expldiff_rtm_thlm, &
                                         l_Lscale_plume_centered, &
                                         l_diag_Lscale_from_tau, &
                                         l_use_C7_Richardson, &
                                         l_use_C11_Richardson, &
                                         l_use_shear_Richardson, &
                                         l_brunt_vaisala_freq_moist, &
                                         l_use_thvm_in_bv_freq, &
                                         l_rcm_supersat_adj, &
                                         l_single_C2_Skw, &
                                         l_damp_wp3_Skw_squared, &
                                         l_prescribed_avg_deltaz, &
                                         l_update_pressure, &
                                         l_lmm_stepping, &
                                         l_use_tke_in_wp3_pr_turb_term, &
                                         l_e3sm_config )

    iiPDF_type = test_pdf_type

    write(fstdout,*) ""
    write(fstdout,*) "Performing PDF parameter values unit test"
    write(fstdout,*) "========================================="
    write(fstdout,*) ""


    if ( test_PDF_type == iiPDF_new ) then
       write(fstdout,*) "Performing PDF parameter unit tests for the new PDF"
       write(fstdout,*) ""
    elseif ( test_PDF_type == iiPDF_new_hybrid ) then
       write(fstdout,*) "Performing PDF parameter unit tests for the new hybrid PDF"
       write(fstdout,*) ""
    elseif ( test_PDF_type == iiPDF_ADG1 ) then
       write(fstdout,*) "Performing PDF parameter unit tests for ADG1"
       write(fstdout,*) ""
    elseif ( test_PDF_type == iiPDF_TSDADG ) then
       write(fstdout,*) "Performing PDF parameter unit tests for TSDADG"
       write(fstdout,*) ""
    elseif ( test_PDF_type == iiPDF_LY93 ) then
       write(fstdout,*) "Performing PDF parameter unit tests for LY93"
       write(fstdout,*) ""
    else ! All other values of test_PDF_type.
       write(fstderr,*) "The PDF parameter unit tests cannot be run for the " &
                        // "selected type of PDF."
       write(fstderr,*) "Selected type of PDF: ", test_PDF_type
       write(fstderr,*) "See src/CLUBB_core/pdf_parameter_module.F90 for an" &
                        // "index of PDF types."
       write(fstderr,*) ""
       pdf_parameter_unit_tests = 1 ! Exit Code = 1, Fail
       return
    endif

    ! Initialize number of failed parameter sets.
    num_failed_sets = 0
    l_failed_sets = .false.
    total_num_failed_sets = 0

    ! The sign of the variance of vertical velocity is always positive.
    sgn_wp2 = one

    ! Set gr%nz to nz.  This is for use within PDF subroutines that input and
    ! output variable arrays of size gr%nz.
    gr%nz = nz

    ! Set dummy values of horz wind variables for ADG1 call
    um =   0.00_core_rknd
    vm =   0.00_core_rknd
    up2 =  1.00_core_rknd
    vp2 =  1.00_core_rknd
    upwp = 0.00_core_rknd
    vpwp = 0.00_core_rknd

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
          ! Large, negative skewness of w.
          wm = 0.01_core_rknd
          wp2 = 0.1_core_rknd
          Skw = -4.5_core_rknd
          wp3 = Skw * wp2**1.5
          rtm = 1.0e-2_core_rknd
          rtp2 = 5.0e-6_core_rknd
          Skrt = 1.0_core_rknd
          rtp3 = Skrt * rtp2**1.5
          thlm = 305.0_core_rknd
          thlp2 = 0.1_core_rknd
          Skthl = -2.5_core_rknd
          thlp3 = Skthl * thlp2**1.5
          wprtp = -5.0e-5_core_rknd
          wpthlp = 5.0e-3_core_rknd
       elseif ( iter_param_sets == 3 ) then
          write(fstdout,*) "PDF parameter set 3:"
          ! Large, positive skewness of w.
          wm = -0.01_core_rknd
          wp2 = 0.5_core_rknd
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
          ! Large variance of w; moderate, positive skewness of w.
          wm = 0.1_core_rknd
          wp2 = 5.0_core_rknd
          Skw = 2.0_core_rknd
          wp3 = Skw * wp2**1.5
          rtm = 1.0e-2_core_rknd
          rtp2 = 7.0e-6_core_rknd
          Skrt = 1.5_core_rknd
          rtp3 = Skrt * rtp2**1.5
          thlm = 305.0_core_rknd
          thlp2 = 1.0_core_rknd
          Skthl = -1.5_core_rknd
          thlp3 = Skthl * thlp2**1.5
          wprtp = 5.0e-5_core_rknd
          wpthlp = -1.0e-2_core_rknd
       elseif ( iter_param_sets == 5 ) then
          write(fstdout,*) "PDF parameter set 5:"
          ! Moderate, negative skewness of w; theta-l with greatest magnitude
          ! of skewness.
          wm = -0.001_core_rknd
          wp2 = 0.3_core_rknd
          Skw = -2.0_core_rknd
          wp3 = Skw * wp2**1.5
          rtm = 1.0e-2_core_rknd
          rtp2 = 2.0e-6_core_rknd
          Skrt = -0.25_core_rknd
          rtp3 = Skrt * rtp2**1.5
          thlm = 305.0_core_rknd
          thlp2 = 0.2_core_rknd
          Skthl = -3.5_core_rknd
          thlp3 = Skthl * thlp2**1.5
          wprtp = 5.0e-5_core_rknd
          wpthlp = 5.0e-3_core_rknd
       elseif ( iter_param_sets == 6 ) then
          write(fstdout,*) "PDF parameter set 6:"
          ! w, rt, and theta-l are all unskewed.
          wm = -0.01_core_rknd
          wp2 = 1.0_core_rknd
          Skw = 0.0_core_rknd
          wp3 = Skw * wp2**1.5
          rtm = 5.0e-3_core_rknd
          rtp2 = 5.0e-7_core_rknd
          Skrt = 0.0_core_rknd
          rtp3 = Skrt * rtp2**1.5
          thlm = 305.0_core_rknd
          thlp2 = 0.1_core_rknd
          Skthl = 0.0_core_rknd
          thlp3 = Skthl * thlp2**1.5
          wprtp = 5.0e-5_core_rknd
          wpthlp = -5.0e-3_core_rknd
       elseif ( iter_param_sets == 7 ) then
          write(fstdout,*) "PDF parameter set 7:"
          ! Small, positive skewness of w; rt with greatest magnitude of
          ! skewness.
          wm = 0.2_core_rknd
          wp2 = 0.25_core_rknd
          Skw = 0.75_core_rknd
          wp3 = Skw * wp2**1.5
          rtm = 1.0e-2_core_rknd
          rtp2 = 5.0e-6_core_rknd
          Skrt = 3.0_core_rknd
          rtp3 = Skrt * rtp2**1.5
          thlm = 305.0_core_rknd
          thlp2 = 0.15_core_rknd
          Skthl = -0.5_core_rknd
          thlp3 = Skthl * thlp2**1.5
          wprtp = 5.0e-5_core_rknd
          wpthlp = -5.0e-3_core_rknd
       elseif ( iter_param_sets == 8 ) then
          write(fstdout,*) "PDF parameter set 8:"
          ! w, rt, and theta-l are all constant.
          wm = -0.002_core_rknd
          wp2 = 0.0_core_rknd
          Skw = 0.0_core_rknd
          wp3 = Skw * wp2**1.5
          rtm = 1.0e-2_core_rknd
          rtp2 = 0.0_core_rknd
          Skrt = 0.0_core_rknd
          rtp3 = Skrt * rtp2**1.5
          thlm = 305.0_core_rknd
          thlp2 = 0.0_core_rknd
          Skthl = 0.0_core_rknd
          thlp3 = Skthl * thlp2**1.5
          wprtp = 0.0e-5_core_rknd
          wpthlp = 0.0e-3_core_rknd
       elseif ( iter_param_sets == 9 ) then
          write(fstdout,*) "PDF parameter set 9:"
          ! w is constant.
          wm = 0.001_core_rknd
          wp2 = 0.0_core_rknd
          Skw = 0.0_core_rknd
          wp3 = Skw * wp2**1.5
          rtm = 1.0e-2_core_rknd
          rtp2 = 2.0e-6_core_rknd
          Skrt = -1.0_core_rknd
          rtp3 = Skrt * rtp2**1.5
          thlm = 305.0_core_rknd
          thlp2 = 0.7_core_rknd
          Skthl = -2.0_core_rknd
          thlp3 = Skthl * thlp2**1.5
          wprtp = 0.0e-5_core_rknd
          wpthlp = 0.0e-3_core_rknd
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
          thlp2 = 3.0_core_rknd * rand6
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
          ! The random number term is the overall correlation of w and rt, and
          ! its value can range from -1 to 1.
          wprtp = ( two * rand10 - one ) * sqrt( wp2 ) * sqrt( rtp2 )
          ! Use a random number to calculate the value of wpthlp.
          ! The random number term is the overall correlation of w and thl, and
          ! its value can range from -1 to 1.
          wpthlp = ( two * rand11 - one ) * sqrt( wp2 ) * sqrt( thlp2 )
       endif ! iter_param_sets == index

       where ( sign( one, wprtp ) * sign( one, wpthlp ) < zero )
          rtpthlp = -0.9_core_rknd * sqrt( rtp2 ) * sqrt( thlp2 )
       elsewhere
          rtpthlp = 0.9_core_rknd * sqrt( rtp2 ) * sqrt( thlp2 )
       endwhere

       where ( sign( one, upwp ) * sign( one, wp3 ) < zero )
          Sku = -0.5_core_rknd * Skw 
       elsewhere
          Sku = 0.5_core_rknd * Skw
       endwhere

       where ( sign( one, vpwp ) * sign( one, wp3 ) < zero )
          Skv = -0.5_core_rknd * Skw 
       elsewhere
          Skv = 0.5_core_rknd * Skw
       endwhere

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
       write(fstdout,*) "rtpthlp = ", rtpthlp
       write(fstdout,*) ""

       !====================================================================
       ! Perform the first set of tests for only the setting variable.
       ! Here, w will be used as the setting variable.
       !====================================================================

       if ( test_PDF_type == iiPDF_new ) then

          write(fstdout,*) "Running tests for the above parameter set for " &
                           // "all combinations of F_w and zeta_w for the " &
                           // "setting variable.  F_w values are 0 (or " &
                           // "1.0 x 10^-5), 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, " &
                           // "0.7, 0.8, 0.9, and 1.  Zeta_w values are 0, " &
                           // "-1/2, 1, -1/3, 1/2, -1/5, 1/4, -9/10, 9, " &
                           // "-3/4, and 3."
          write(fstdout,*) ""

          do iter_F_w = 1, num_F_w, 1

             ! Set the value of F_w.
             if ( iter_F_w == 1 ) then
                where ( abs( Skw ) > zero )
                   ! The value of F_w needs to be greater than 0 when
                   ! | Skw | > 0 in order for the PDF to be valid.
                   F_w = 1.0e-5_core_rknd
                elsewhere ! Skw = 0
                   ! F_w can have a value of 0 when Skw = 0.
                   F_w = zero
                endwhere ! | Skw | > 0
             else ! iter_F_w > 1
                ! F_w ranges from 0 to 1.
                F_w = 0.1_core_rknd * real( iter_F_w - 1, kind = core_rknd )
             endif ! iter_F_w

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
                call calc_setter_var_params( gr, wm, wp2, Skw, sgn_wp2,     & ! In
                                             F_w, zeta_w,               & ! In
                                             mu_w_1, mu_w_2, sigma_w_1, & ! Out
                                             sigma_w_2, mixt_frac,      & ! Out
                                             coef_sigma_w_1_sqd,        & ! Out
                                             coef_sigma_w_2_sqd         ) ! Out

                sigma_w_1_sqd = sigma_w_1**2
                sigma_w_2_sqd = sigma_w_2**2

                l_check_mu_w_1_gte_mu_w_2 = .true.

                ! Perform the tests for the "setter" variable, which is the
                ! variable that is used to set the mixture fraction.
                call setter_var_tests( nz, wm, wp2, wp3, Skw,        & ! In
                                       mu_w_1, mu_w_2, sigma_w_1,    & ! In
                                       sigma_w_2, mixt_frac, tol,    & ! In
                                       sigma_w_1_sqd, sigma_w_2_sqd, & ! In
                                       l_pass_test_1, l_pass_test_2, & ! Out
                                       l_pass_test_3, l_pass_test_4, & ! Out
                                       l_pass_test_5, l_pass_test_6, & ! Out
                                       l_check_mu_w_1_gte_mu_w_2     ) ! Out

                ! Calculate <w'^4> by integrating over the PDF.
                wp4_pdf_calc = calc_wp4_pdf( gr, wm, mu_w_1, mu_w_2, sigma_w_1**2, &
                                             sigma_w_2**2, mixt_frac )

                ! Calculate <w'^4> by <w'^4> = coef_wp4_implicit * <w'^2>^2.
                coef_wp4_implicit &
                = calc_coef_wp4_implicit( gr, mixt_frac, F_w, &
                                          coef_sigma_w_1_sqd, &
                                          coef_sigma_w_2_sqd )

                wp4_implicit_calc = coef_wp4_implicit * wp2**2

                ! Test 7
                ! Compare the <w'^4> calculated by the PDF to its value
                ! calculated by <w'^4> = wp4_implicit_calc * <w'^2>^2 (which was
                ! derived from the PDF).
                !    | ( <w'^4>|_pdf - <w'^4>|_impl ) / <w'^4>|_pdf |  <=  tol;
                ! which can be rewritten as:
                !    | <w'^4>|_pdf - <w'^4>|_impl |  <=  <w'^4>|_pdf * tol.
                where ( abs( wp4_pdf_calc - wp4_implicit_calc ) &
                        <= max( wp4_pdf_calc, w_tol**4 ) * tol )
                   l_pass_test_7 = .true.
                elsewhere
                   l_pass_test_7 = .false.
                endwhere

                if ( any( .not. l_pass_test_7 ) ) then
                   do idx = 1, nz, 1
                      if ( .not. l_pass_test_7(idx) ) then
                         write(fstderr,*) "Test 7 failed"
                         !write(fstderr,*) "index = ", idx
                         write(fstderr,*) "wp4 pdf = ", wp4_pdf_calc(idx)
                         write(fstderr,*) "wp4 implicit calc = ", &
                                          wp4_implicit_calc(idx)
                         write(fstderr,*) ""
                      endif ! .not. l_pass_test_7(idx)
                   enddo ! idx = 1, nz, 1
                endif

                where ( l_pass_test_1 .and. l_pass_test_2 .and. l_pass_test_3 &
                        .and. l_pass_test_4 .and. l_pass_test_5 &
                        .and. l_pass_test_6 .and. l_pass_test_7 )
                   ! All tests pass
                   num_failed_sets = num_failed_sets
                   l_failed_sets = .false.
                elsewhere
                   ! At least one test failed
                   num_failed_sets = num_failed_sets + 1
                   l_failed_sets = .true.
                endwhere

                if ( any( l_failed_sets ) ) then
                   do idx = 1, nz, 1
                      if ( l_failed_sets(idx) ) then
                         write(fstderr,*) "At least one test or check for " &
                                          // "the setting variable PDF " &
                                          // "failed for the following " &
                                          // "parameter set:  "
                         write(fstderr,*) "PDF parameter set index = ", &
                                          iter_param_sets
                         write(fstderr,*) "F_w = ", F_w(idx)
                         write(fstderr,*) "zeta_w = ", zeta_w(idx)
                         write(fstderr,*) ""
                      endif ! l_failed_sets(idx)
                   enddo ! idx = 1, nz, 1
                endif

             enddo ! iter_zeta_w = 1, num_zeta_w, 1

          enddo ! iter_F_w = 1, num_F_w, 1

       elseif ( test_PDF_type == iiPDF_new_hybrid ) then

          write(fstdout,*) "Running tests for the above parameter set for " &
                           // "all combinations of F_w and zeta_w for the " &
                           // "setting variable.  F_w values are 0 (or " &
                           // "1.0 x 10^-5), 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, " &
                           // "0.7, 0.8, 0.9, and 1.  Zeta_w values are 0, " &
                           // "-1/2, 1, -1/3, 1/2, -1/5, 1/4, -9/10, 9, " &
                           // "-3/4, and 3."
          write(fstdout,*) ""

          do iter_F_w = 1, num_F_w, 1

             ! Set the value of F_w.
             if ( iter_F_w == 1 ) then
                where ( abs( Skw ) > zero )
                   ! The value of F_w needs to be greater than 0 when
                   ! | Skw | > 0 in order for the PDF to be valid.
                   F_w = 1.0e-5_core_rknd
                elsewhere ! Skw = 0
                   ! F_w can have a value of 0 when Skw = 0.
                   F_w = zero
                endwhere ! | Skw | > 0
             else ! iter_F_w > 1
                ! F_w ranges from 0 to 1.
                F_w = 0.1_core_rknd * real( iter_F_w - 1, kind = core_rknd )
             endif ! iter_F_w

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
                call calculate_w_params( wm, wp2, Skw, F_w, zeta_w, & ! In
                                         mu_w_1, mu_w_2, sigma_w_1, & ! Out
                                         sigma_w_2, mixt_frac,      & ! Out
                                         coef_sigma_w_1_sqd,        & ! Out
                                         coef_sigma_w_2_sqd         ) ! Out

                sigma_w_1_sqd = sigma_w_1**2
                sigma_w_2_sqd = sigma_w_2**2

                l_check_mu_w_1_gte_mu_w_2 = .true.

                ! Perform the tests for the "setter" variable, which is the
                ! variable that is used to set the mixture fraction.
                call setter_var_tests( nz, wm, wp2, wp3, Skw,        & ! In
                                       mu_w_1, mu_w_2, sigma_w_1,    & ! In
                                       sigma_w_2, mixt_frac, tol,    & ! In
                                       sigma_w_1_sqd, sigma_w_2_sqd, & ! In
                                       l_pass_test_1, l_pass_test_2, & ! Out
                                       l_pass_test_3, l_pass_test_4, & ! Out
                                       l_pass_test_5, l_pass_test_6, & ! Out
                                       l_check_mu_w_1_gte_mu_w_2     ) ! Out

                ! Calculate <w'^4> by integrating over the PDF.
                wp4_pdf_calc = calc_wp4_pdf( gr, wm, mu_w_1, mu_w_2, sigma_w_1**2, &
                                             sigma_w_2**2, mixt_frac )

                ! Calculate <w'^4> by <w'^4> = coef_wp4_implicit * <w'^2>^2.
                coef_wp4_implicit &
                = calculate_coef_wp4_implicit( mixt_frac, F_w, &
                                               coef_sigma_w_1_sqd, &
                                               coef_sigma_w_2_sqd )

                wp4_implicit_calc = coef_wp4_implicit * wp2**2

                ! Test 7
                ! Compare the <w'^4> calculated by the PDF to its value
                ! calculated by <w'^4> = wp4_implicit_calc * <w'^2>^2 (which was
                ! derived from the PDF).
                !    | ( <w'^4>|_pdf - <w'^4>|_impl ) / <w'^4>|_pdf |  <=  tol;
                ! which can be rewritten as:
                !    | <w'^4>|_pdf - <w'^4>|_impl |  <=  <w'^4>|_pdf * tol.
                where ( abs( wp4_pdf_calc - wp4_implicit_calc ) &
                        <= max( wp4_pdf_calc, w_tol**4 ) * tol )
                   l_pass_test_7 = .true.
                elsewhere
                   l_pass_test_7 = .false.
                endwhere

                if ( any( .not. l_pass_test_7 ) ) then
                   do idx = 1, nz, 1
                      if ( .not. l_pass_test_7(idx) ) then
                         write(fstderr,*) "Test 7 failed"
                         !write(fstderr,*) "index = ", idx
                         write(fstderr,*) "wp4 pdf = ", wp4_pdf_calc(idx)
                         write(fstderr,*) "wp4 implicit calc = ", &
                                          wp4_implicit_calc(idx)
                         write(fstderr,*) ""
                      endif ! .not. l_pass_test_7(idx)
                   enddo ! idx = 1, nz, 1
                endif

                where ( l_pass_test_1 .and. l_pass_test_2 .and. l_pass_test_3 &
                        .and. l_pass_test_4 .and. l_pass_test_5 &
                        .and. l_pass_test_6 .and. l_pass_test_7 )
                   ! All tests pass
                   num_failed_sets = num_failed_sets
                   l_failed_sets = .false.
                elsewhere
                   ! At least one test failed
                   num_failed_sets = num_failed_sets + 1
                   l_failed_sets = .true.
                endwhere

                if ( any( l_failed_sets ) ) then
                   do idx = 1, nz, 1
                      if ( l_failed_sets(idx) ) then
                         write(fstderr,*) "At least one test or check for " &
                                          // "the setting variable PDF " &
                                          // "failed for the following " &
                                          // "parameter set:  "
                         write(fstderr,*) "PDF parameter set index = ", &
                                          iter_param_sets
                         write(fstderr,*) "F_w = ", F_w(idx)
                         write(fstderr,*) "zeta_w = ", zeta_w(idx)
                         write(fstderr,*) ""
                      endif ! l_failed_sets(idx)
                   enddo ! idx = 1, nz, 1
                endif

             enddo ! iter_zeta_w = 1, num_zeta_w, 1

          enddo ! iter_F_w = 1, num_F_w, 1

       elseif ( test_PDF_type == iiPDF_ADG1 ) then

          write(fstdout,*) "Running tests for the above parameter set for " &
                           // "various values of sigma_sqd_w for the setting " &
                           // "variable (w).  Values of sigma_sqd_w are 0, " &
                           // "0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, " &
                           // "and 1 (or 0.99)."
          write(fstdout,*) ""

          sqrt_wp2 = sqrt( wp2 )
          mixt_frac_max_mag = 1.0_core_rknd

          do iter_sigma_sqd_w = 1, num_sigma_sqd_w, 1

             ! Set the value of sigma_sqd_w.
             if ( iter_sigma_sqd_w == num_sigma_sqd_w ) then
                where ( abs( Skw ) > zero )
                   ! The value of sigma_sqd_w needs to be less than than 1 when
                   ! | Skw | > 0 in order for the PDF to be valid.
                   sigma_sqd_w = 0.99_core_rknd
                elsewhere ! Skw = 0
                   ! sigma_sqd_w can have a value of 0 when Skw = 0.
                   sigma_sqd_w = one
                endwhere ! | Skw | > 0
             else ! iter_sigma_sqd_w > 1
                ! sigma_sqd_w ranges from 0 to 1.
                sigma_sqd_w = 0.1_core_rknd * real( iter_sigma_sqd_w - 1, &
                                                    kind = core_rknd )
             endif ! iter_sigma_sqd_w
 
             call ADG1_w_closure( wm, wp2, Skw, sigma_sqd_w,              &! In
                                  sqrt_wp2, mixt_frac_max_mag,            &! In
                                  mu_w_1, mu_w_2, w_1_n, w_2_n,           &! Out
                                  sigma_w_1_sqd, sigma_w_2_sqd, mixt_frac )! Out

             sigma_w_1 = sqrt( max( sigma_w_1_sqd, zero ) )
             sigma_w_2 = sqrt( max( sigma_w_2_sqd, zero ) )

             l_check_mu_w_1_gte_mu_w_2 = .true.

             ! Perform the tests for the "setter" variable, which is the
             ! variable that is used to set the mixture fraction.  This is
             ! always w for ADG1.
             call setter_var_tests( nz, wm, wp2, wp3, Skw,        & ! In
                                    mu_w_1, mu_w_2, sigma_w_1,    & ! In
                                    sigma_w_2, mixt_frac, tol,    & ! In
                                    sigma_w_1_sqd, sigma_w_2_sqd, & ! In
                                    l_pass_test_1, l_pass_test_2, & ! Out
                                    l_pass_test_3, l_pass_test_4, & ! Out
                                    l_pass_test_5, l_pass_test_6, & ! Out
                                    l_check_mu_w_1_gte_mu_w_2     ) ! Out


             where ( l_pass_test_1 .and. l_pass_test_2 .and. l_pass_test_3 &
                     .and. l_pass_test_4 .and. l_pass_test_5 &
                     .and. l_pass_test_6 )
                ! All tests pass
                num_failed_sets = num_failed_sets
                l_failed_sets = .false.
             elsewhere
                ! At least one test failed
                num_failed_sets = num_failed_sets + 1
                l_failed_sets = .true.
             endwhere

             if ( any( l_failed_sets ) ) then
                do idx = 1, nz, 1
                   if ( l_failed_sets(idx) ) then
                      write(fstderr,*) "At least one test or check for the " &
                                       // "setting variable PDF failed for " &
                                       // "the following parameter set:  "
                      write(fstderr,*) "PDF parameter set index = ", &
                                       iter_param_sets
                      write(fstderr,*) "sigma_sqd_w = ", sigma_sqd_w(idx)
                      write(fstderr,*) ""
                   endif ! l_failed_sets(idx)
                enddo ! idx = 1, nz, 1
             endif

          enddo ! iter_sigma_sqd_w = 1, num_sigma_sqd_w, 1

       elseif ( test_PDF_type == iiPDF_TSDADG ) then

          do iter_small_l_w_1 = 1, num_small_l_w_1, 1

            small_l_w_1 &
            = two_thirds + one_third &
                           * real( iter_small_l_w_1, kind = core_rknd ) &
                             / real( num_small_l_w_1, kind = core_rknd )

            do iter_small_l_w_2 = 1, num_small_l_w_2, 1

              small_l_w_2 = max( 0.1_core_rknd * real( iter_small_l_w_2 - 1, &
                                                       kind = core_rknd ), &
                                 5.0e-3_core_rknd )

              do idx = 1, nz, 1

                 call calc_L_x_Skx_fnc( Skw(idx), sgn_wp2(idx),        & ! In
                                        small_l_w_1, small_l_w_2,      & ! In
                                        big_L_w_1(idx), big_L_w_2(idx) ) ! Out

                 call calc_setter_parameters( wm(idx), wp2(idx),        & ! In
                                              Skw(idx), sgn_wp2(idx),   & ! In
                                              big_L_w_1(idx),           & ! In
                                              big_L_w_2(idx),           & ! In
                                              mu_w_1(idx), mu_w_2(idx), & ! Out
                                              sigma_w_1_sqd(idx),       & ! Out
                                              sigma_w_2_sqd(idx),       & ! Out
                                              mixt_frac(idx),           & ! Out
                                              coef_sigma_w_1_sqd(idx),  & ! Out
                                              coef_sigma_w_2_sqd(idx)   ) ! Out

              enddo ! idx = 1, nz, 1

              sigma_w_1 = sqrt( max( sigma_w_1_sqd, zero ) )
              sigma_w_2 = sqrt( max( sigma_w_2_sqd, zero ) )

              l_check_mu_w_1_gte_mu_w_2 = .true.

              ! Perform the tests for the "setter" variable, which is the
              ! variable that is used to set the mixture fraction.
              call setter_var_tests( nz, wm, wp2, wp3, Skw,        & ! In
                                     mu_w_1, mu_w_2, sigma_w_1,    & ! In
                                     sigma_w_2, mixt_frac, tol,    & ! In
                                     sigma_w_1_sqd, sigma_w_2_sqd, & ! In
                                     l_pass_test_1, l_pass_test_2, & ! Out
                                     l_pass_test_3, l_pass_test_4, & ! Out
                                     l_pass_test_5, l_pass_test_6, & ! Out
                                     l_check_mu_w_1_gte_mu_w_2     ) ! Out


              where ( l_pass_test_1 .and. l_pass_test_2 .and. l_pass_test_3 &
                   .and. l_pass_test_4 .and. l_pass_test_5 &
                   .and. l_pass_test_6 )
                 ! All tests pass
                 num_failed_sets = num_failed_sets
                 l_failed_sets = .false.
              elsewhere
                 ! At least one test failed
                 num_failed_sets = num_failed_sets + 1
                 l_failed_sets = .true.
              endwhere

              if ( any( l_failed_sets ) ) then
                 do idx = 1, nz, 1
                    if ( l_failed_sets(idx) ) then
                       write(fstderr,*) "At least one test or check for the " &
                                        // "setting variable PDF failed for " &
                                        // "the following parameter set:  "
                         write(fstderr,*) "PDF parameter set index = ", &
                                          iter_param_sets
                         write(fstderr,*) "l_w_1 = ", small_l_w_1
                         write(fstderr,*) "l_w_2 = ", small_l_w_2
                         write(fstderr,*) "L_w_1 = ", big_L_w_1(idx)
                         write(fstderr,*) "L_w_2 = ", big_L_w_2(idx)
                         write(fstderr,*) ""
                    endif ! l_failed_sets(idx)
                 enddo ! idx = 1, nz, 1
              endif

            enddo ! iter_small_l_w_2 = 1, num_small_l_w_2, 1

          enddo ! iter_small_l_w_1 = 1, num_small_l_w_1, 1

       elseif ( test_PDF_type == iiPDF_LY93 ) then

          write(fstdout,*) "Running tests for the above parameter set for " &
                           // "the setting variable (w is used here)."
          write(fstdout,*) ""

          mixt_frac = calc_mixt_frac_LY93( gr, abs( Skw ) )

          call calc_params_LY93( gr, wm, wp2, Skw, mixt_frac,     & ! In
                                 mu_w_1, mu_w_2,              & ! Out
                                 sigma_w_1_sqd, sigma_w_2_sqd ) ! Out

          sigma_w_1 = sqrt( max( sigma_w_1_sqd, zero ) )
          sigma_w_2 = sqrt( max( sigma_w_2_sqd, zero ) )

          l_check_mu_w_1_gte_mu_w_2 = .false.

          ! Perform the tests for the "setter" variable, which is the variable
          ! that is used to set the mixture fraction.
          call setter_var_tests( nz, wm, wp2, wp3, Skw,        & ! In
                                 mu_w_1, mu_w_2, sigma_w_1,    & ! In
                                 sigma_w_2, mixt_frac, tol,    & ! In
                                 sigma_w_1_sqd, sigma_w_2_sqd, & ! In
                                 l_pass_test_1, l_pass_test_2, & ! Out
                                 l_pass_test_3, l_pass_test_4, & ! Out
                                 l_pass_test_5, l_pass_test_6, & ! Out
                                 l_check_mu_w_1_gte_mu_w_2     ) ! Out


          ! Note:  l_pass_test_6 -- LY93 doesn't necessarily make
          !        mu_w_1 >= mu_w_2, so that test is thrown out.
          where ( l_pass_test_1 .and. l_pass_test_2 .and. l_pass_test_3 &
                  .and. l_pass_test_4 .and. l_pass_test_5 )
             ! All tests pass
             num_failed_sets = num_failed_sets
             l_failed_sets = .false.
          elsewhere
             ! At least one test failed
             num_failed_sets = num_failed_sets + 1
             l_failed_sets = .true.
          endwhere

          if ( any( l_failed_sets ) ) then
             do idx = 1, nz, 1
                if ( l_failed_sets(idx) ) then
                   write(fstderr,*) "At least one test or check for the " &
                                    // "setting variable PDF failed for the " &
                                    // "following parameter set:  "
                   write(fstderr,*) "PDF parameter set index = ", &
                                    iter_param_sets
                   write(fstderr,*) ""
                endif ! l_failed_sets(idx)
             enddo ! idx = 1, nz, 1
          endif

       endif ! test_PDF_type

       !====================================================================
       ! Perform the second set of tests for the full PDF.
       !====================================================================

       if ( test_PDF_type == iiPDF_new ) then

          write(fstdout,*) "Running tests for the above parameter set for " &
                           // "the full PDF (the setting of F and zeta " &
                           // "values is handled internally)."
          write(fstdout,*) ""

          call new_pdf_driver( gr, wm, rtm, thlm, wp2, rtp2, thlp2, Skw,    & ! In
                               wprtp, wpthlp, rtpthlp,                  & ! In
                               Skrt, Skthl,                             & ! I/O
                               mu_w_1, mu_w_2, mu_rt_1, mu_rt_2,        & ! Out
                               mu_thl_1, mu_thl_2, sigma_w_1_sqd,       & ! Out
                               sigma_w_2_sqd, sigma_rt_1_sqd,           & ! Out
                               sigma_rt_2_sqd, sigma_thl_1_sqd,         & ! Out
                               sigma_thl_2_sqd, mixt_frac,              & ! Out
                               pdf_implicit_coefs_terms,                & ! Out
                               F_w, F_rt, F_thl, min_F_w, max_F_w,      & ! Out
                               min_F_rt, max_F_rt, min_F_thl, max_F_thl ) ! Out

          ! Recalculate <rt'^3> and <thl'^3> just in case Skrt and Skthl needed
          ! to be clipped in new_pdf_driver.
          rtp3 = Skrt * rtp2**1.5
          thlp3 = Skthl * thlp2**1.5

       elseif ( test_PDF_type == iiPDF_new_hybrid ) then

          write(fstdout,*) "Running tests for the above parameter set for " &
                           // "the full PDF (the setting of F and zeta " &
                           // "values is handled internally)."
          write(fstdout,*) ""

          call new_hybrid_pdf_driver( gr, wm, rtm, thlm, um, vm,              &! In
                                      wp2, rtp2, thlp2, up2, vp2,         &! In
                                      Skw, wprtp, wpthlp, upwp, vpwp,     &! In
                                      sclrm, sclrp2, wpsclrp,             &! In
                                      Skrt, Skthl, Sku, Skv, Sksclr,      &! I/O
                                      mu_w_1, mu_w_2,                     &! Out
                                      mu_rt_1, mu_rt_2,                   &! Out
                                      mu_thl_1, mu_thl_2,                 &! Out
                                      mu_u_1, mu_u_2, mu_v_1, mu_v_2,     &! Out
                                      sigma_w_1_sqd, sigma_w_2_sqd,       &! Out
                                      sigma_rt_1_sqd, sigma_rt_2_sqd,     &! Out
                                      sigma_thl_1_sqd, sigma_thl_2_sqd,   &! Out
                                      sigma_u_1_sqd, sigma_u_2_sqd,       &! Out
                                      sigma_v_1_sqd, sigma_v_2_sqd,       &! Out
                                      mu_sclr_1, mu_sclr_2,               &! Out
                                      sigma_sclr_1_sqd, sigma_sclr_2_sqd, &! Out
                                      mixt_frac,                          &! Out
                                      pdf_implicit_coefs_terms,           &! Out
                                      F_w, min_F_w, max_F_w               )! Out

          ! Recalculate <rt'^3> and <thl'^3> just in case Skrt and Skthl needed
          ! to be clipped in new_pdf_driver.
          rtp3 = Skrt * rtp2**1.5
          thlp3 = Skthl * thlp2**1.5

       elseif ( test_PDF_type == iiPDF_ADG1 ) then

          write(fstdout,*) "Running tests for the above parameter set for " &
                           // "the full PDF (the setting of sigma_sqd_w is " &
                           // "handled the same way as in the model code)."
          write(fstdout,*) ""

          if ( l_gamma_Skw ) then
             gamma_Skw_fnc = gamma_coefb &
                             + ( gamma_coef - gamma_coefb ) &
                               * exp( -one_half * ( Skw / gamma_coefc )**2 )
          else
             gamma_Skw_fnc = gamma_coef
          endif

          sigma_sqd_w &
          = compute_sigma_sqd_w( gamma_Skw_fnc, wp2, thlp2, rtp2, &
                                 up2, vp2, wpthlp, wprtp, upwp, vpwp, &
                                 l_predict_upwp_vpwp )

          call ADG1_pdf_driver( gr, wm, rtm, thlm, um, vm,                  & ! In 
                               wp2, rtp2, thlp2, up2, vp2,              & ! In 
                               Skw, wprtp, wpthlp, upwp, vpwp, sqrt_wp2,& ! In 
                               sigma_sqd_w, mixt_frac_max_mag,          & ! In 
                               sclrm, sclrp2, wpsclrp, l_scalar_calc,   & ! In 
                               mu_w_1, mu_w_2, mu_rt_1, mu_rt_2, mu_thl_1, mu_thl_2,& ! Out
                               mu_u_1, mu_u_2, mu_v_1, mu_v_2,          & ! Out
                               sigma_w_1_sqd, sigma_w_2_sqd, sigma_rt_1_sqd, & ! Out
                               sigma_rt_2_sqd, sigma_thl_1_sqd, sigma_thl_2_sqd, & ! Out
                               sigma_u_1_sqd, sigma_u_2_sqd,            & ! Out
                               sigma_v_1_sqd, sigma_v_2_sqd,            & ! Out
                               mixt_frac, alpha_rt, alpha_thl,          & ! Out
                               alpha_u, alpha_v,                        & ! Out
                               mu_sclr_1, mu_sclr_2, sigma_sclr_1_sqd,  & ! Out
                               sigma_sclr_2_sqd, alpha_sclr )             ! Out

       elseif ( test_PDF_type == iiPDF_TSDADG ) then

          ! Temporarily skip this step.
          cycle

       elseif ( test_PDF_type == iiPDF_LY93 ) then

          write(fstdout,*) "Running tests for the above parameter set for " &
                           // "the full PDF."
          write(fstdout,*) ""

          call LY93_driver( gr, wm, rtm, thlm, wp2, rtp2,          & ! In
                            thlp2, Skw, Skrt, Skthl,           & ! In
                            mu_w_1, mu_w_2, mu_rt_1, mu_rt_2,  & ! Out
                            mu_thl_1, mu_thl_2, sigma_w_1_sqd, & ! Out
                            sigma_w_2_sqd, sigma_rt_1_sqd,     & ! Out
                            sigma_rt_2_sqd, sigma_thl_1_sqd,   & ! Out
                            sigma_thl_2_sqd, mixt_frac         ) ! Out

       endif ! test_PDF_type

       sigma_w_1 = sqrt( sigma_w_1_sqd )
       sigma_w_2 = sqrt( sigma_w_2_sqd )
       sigma_rt_1 = sqrt( sigma_rt_1_sqd )
       sigma_rt_2 = sqrt( sigma_rt_2_sqd )
       sigma_thl_1 = sqrt( sigma_thl_1_sqd )
       sigma_thl_2 = sqrt( sigma_thl_2_sqd )

       ! Recalculate the values of <wm>, <w'^2>, <w'^3>, and Skw using the
       ! PDF parameters that are output from the subroutine.
       call recalc_single_var( nz, wm, mu_w_1, mu_w_2, sigma_w_1, & ! In
                               sigma_w_2, mixt_frac,              & ! In
                               recalc_wm, recalc_wp2,             & ! Out
                               recalc_wp3, recalc_Skw             ) ! Out

       ! Recalculate the values of <rtm>, <rt'^2>, <rt'^3>, and Skrt using the
       ! PDF parameters that are output from the subroutine.
       call recalc_single_var( nz, rtm, mu_rt_1, mu_rt_2, sigma_rt_1, & ! In
                               sigma_rt_2, mixt_frac,                 & ! In
                               recalc_rtm, recalc_rtp2,               & ! Out
                               recalc_rtp3, recalc_Skrt               ) ! Out

       ! Recalculate the values of <thlm>, <thl'^2>, <thl'^3>, and Skthl using
       ! the PDF parameters that are output from the subroutine.
       call recalc_single_var( nz, thlm, mu_thl_1, mu_thl_2, sigma_thl_1, & ! In
                               sigma_thl_2, mixt_frac,                    & ! In
                               recalc_thlm, recalc_thlp2,                & ! Out
                               recalc_thlp3, recalc_Skthl                ) ! Out

       ! Test 8
       ! Compare the original wm to its recalculated value.
       !    | ( <w>|_recalc - <w> ) / <w> |  <=  tol;
       ! which can be rewritten as:
       !    | <w>|_recalc - <w> |  <=  | <w> | * tol.
       where ( abs( recalc_wm - wm ) <= max( abs( wm ), w_tol ) * tol )
          l_pass_test_8 = .true.
       elsewhere
          l_pass_test_8 = .false.
       endwhere

       if ( any( .not. l_pass_test_8 ) ) then
          do idx = 1, nz, 1
             if ( .not. l_pass_test_8(idx) ) then
                write(fstderr,*) "Test 8 failed"
                !write(fstderr,*) "index = ", idx
                write(fstderr,*) "wm = ", wm(idx)
                write(fstderr,*) "recalc_wm = ", recalc_wm(idx)
                write(fstderr,*) ""
             endif ! .not. l_pass_test_8(idx)
          enddo ! idx = 1, nz, 1
       endif

       ! Test 9
       ! Compare the original wp2 to its recalculated value.
       !    | ( <w'^2>|_recalc - <w'^2> ) / <w'^2> |  <=  tol;
       ! which can be rewritten as:
       !    | <w'^2>|_recalc - <w'^2> |  <=  | <w'^2> | * tol;
       ! and since <w'^2> is always positive:
       !    | <w'^2>|_recalc - <w'^2> |  <=  <w'^2> * tol.
       where ( abs( recalc_wp2 - wp2 ) <= max( wp2, w_tol_sqd ) * tol )
          l_pass_test_9 = .true.
       elsewhere
          l_pass_test_9 = .false.
       endwhere

       if ( any( .not. l_pass_test_9 ) ) then
          do idx = 1, nz, 1
             if ( .not. l_pass_test_9(idx) ) then
                write(fstderr,*) "Test 9 failed"
                !write(fstderr,*) "index = ", idx
                write(fstderr,*) "wp2 = ", wp2(idx)
                write(fstderr,*) "recalc_wp2 = ", recalc_wp2(idx)
                write(fstderr,*) ""
             endif ! .not. l_pass_test_9(idx)
          enddo ! idx = 1, nz, 1
       endif

       ! Test 10
       ! Compare the original wp3 to its recalculated value.
       !    | ( <w'^3>|_recalc - <w'^3> ) / <w'^3> |  <=  tol;
       ! which can be rewritten as:
       !    | <w'^3>|_recalc - <w'^3> |  <=  | <w'^3> | * tol.
       where ( abs( recalc_wp3 - wp3 ) <= max( abs( wp3 ), w_tol**3 ) * tol )
          l_pass_test_10 = .true.
       elsewhere
          l_pass_test_10 = .false.
       endwhere

       if ( any( .not. l_pass_test_10 ) ) then
          do idx = 1, nz, 1
             if ( .not. l_pass_test_10(idx) ) then
                write(fstderr,*) "Test 10 failed"
                !write(fstderr,*) "index = ", idx
                write(fstderr,*) "wp3 = ", wp3(idx)
                write(fstderr,*) "recalc_wp3 = ", recalc_wp3(idx)
                write(fstderr,*) ""
             endif ! .not. l_pass_test_10(idx)
          enddo ! idx = 1, nz, 1
       endif

       ! Test 11
       ! Compare the original Skw to its recalculated value.
       !    | ( Skw|_recalc - Skw ) / Skw |  <=  tol;
       ! which can be rewritten as:
       !    | Skw|_recalc - Skw |  <=  | Skw | * tol.
       where ( abs( recalc_Skw - Skw ) &
               <= max( abs( Skw ), 1.0e-3_core_rknd ) * tol )
          l_pass_test_11 = .true.
       elsewhere
          l_pass_test_11 = .false.
       endwhere

       if ( any( .not. l_pass_test_11 ) ) then
          do idx = 1, nz, 1
             if ( .not. l_pass_test_11(idx) ) then
                write(fstderr,*) "Test 11 failed"
                !write(fstderr,*) "index = ", idx
                write(fstderr,*) "Skw = ", Skw(idx)
                write(fstderr,*) "recalc_Skw = ", recalc_Skw(idx)
                write(fstderr,*) ""
             endif ! .not. l_pass_test_11(idx)
          enddo ! idx = 1, nz, 1
       endif

       ! Test 12
       ! Compare the original rtm to its recalculated value.
       !    | ( <rt>|_recalc - <rt> ) / <rt> |  <=  tol;
       ! which can be rewritten as:
       !    | <rt>|_recalc - <rt> |  <=  | <rt> | * tol.
       where ( abs( recalc_rtm - rtm ) <= max( abs( rtm ), rt_tol ) * tol )
          l_pass_test_12 = .true.
       elsewhere
          l_pass_test_12 = .false.
       endwhere

       if ( any( .not. l_pass_test_12 ) ) then
          do idx = 1, nz, 1
             if ( .not. l_pass_test_12(idx) ) then
                write(fstderr,*) "Test 12 failed"
                !write(fstderr,*) "index = ", idx
                write(fstderr,*) "rtm = ", rtm(idx)
                write(fstderr,*) "recalc_rtm = ", recalc_rtm(idx)
                write(fstderr,*) ""
             endif ! .not. l_pass_test_12(idx)
          enddo ! idx = 1, nz, 1
       endif

       ! Test 13
       ! Compare the original rtp2 to its recalculated value.
       !    | ( <rt'^2>|_recalc - <rt'^2> ) / <rt'^2> |  <=  tol;
       ! which can be rewritten as:
       !    | <rt'^2>|_recalc - <rt'^2> |  <=  | <rt'^2> | * tol;
       ! and since <rt'^2> is always positive:
       !    | <rt'^2>|_recalc - <rt'^2> |  <=  <rt'^2> * tol.
       where ( abs( recalc_rtp2 - rtp2 ) <= max( rtp2, rt_tol**2 ) * tol )
          l_pass_test_13 = .true.
       elsewhere
          l_pass_test_13 = .false.
       endwhere

       if ( any( .not. l_pass_test_13 ) ) then
          do idx = 1, nz, 1
             if ( .not. l_pass_test_13(idx) ) then
                write(fstderr,*) "Test 13 failed"
                !write(fstderr,*) "index = ", idx
                write(fstderr,*) "rtp2 = ", rtp2(idx)
                write(fstderr,*) "recalc_rtp2 = ", recalc_rtp2(idx)
                write(fstderr,*) ""
             endif ! .not. l_pass_test_13(idx)
          enddo ! idx = 1, nz, 1
       endif

       if ( ( test_PDF_type == iiPDF_new ) &
            .or. ( test_PDF_type == iiPDF_new_hybrid ) &
            .or. ( test_PDF_type == iiPDF_LY93 ) ) then

          ! Test 14
          ! Compare the original rtp3 to its recalculated value.
          !    | ( <rt'^3>|_recalc - <rt'^3> ) / <rt'^3> |  <=  tol;
          ! which can be rewritten as:
          !    | <rt'^3>|_recalc - <rt'^3> |  <=  | <rt'^3> | * tol.
          where ( abs( recalc_rtp3 - rtp3 ) &
                  <= max( abs( rtp3 ), 1.0e-15_core_rknd ) * tol )
             l_pass_test_14 = .true.
          elsewhere
             l_pass_test_14 = .false.
          endwhere

          if ( any( .not. l_pass_test_14 ) ) then
             do idx = 1, nz, 1
                if ( .not. l_pass_test_14(idx) ) then
                   write(fstderr,*) "Test 14 failed"
                   !write(fstderr,*) "index = ", idx
                   write(fstderr,*) "rtp3 = ", rtp3(idx)
                   write(fstderr,*) "recalc_rtp3 = ", recalc_rtp3(idx)
                   write(fstderr,*) ""
                endif ! .not. l_pass_test_14(idx)
             enddo ! idx = 1, nz, 1
          endif

          ! Test 15
          ! Compare the original Skrt to its recalculated value.
          !    | ( Skrt|_recalc - Skrt ) / Skrt |  <=  tol;
          ! which can be rewritten as:
          !    | Skrt|_recalc - Skrt |  <=  | Skrt | * tol.
          where ( abs( recalc_Skrt - Skrt ) &
                  <= max( abs( Skrt ), 1.0e-3_core_rknd ) * tol )
             l_pass_test_15 = .true.
          elsewhere
             l_pass_test_15 = .false.
          endwhere

          if ( any( .not. l_pass_test_15 ) ) then
             do idx = 1, nz, 1
                if ( .not. l_pass_test_15(idx) ) then
                   write(fstderr,*) "Test 15 failed"
                   !write(fstderr,*) "index = ", idx
                   write(fstderr,*) "Skrt = ", Skrt(idx)
                   write(fstderr,*) "recalc_Skrt = ", recalc_Skrt(idx)
                   write(fstderr,*) ""
                   endif ! .not. l_pass_test_15(idx)
             enddo ! idx = 1, nz, 1
          endif

       else

          ! Automatically set the pass test flags to true for PDF types that
          ! don't necessarily preserve <rt'^3> and Skrt.
          l_pass_test_14 = .true.
          l_pass_test_15 = .true.

       endif ! test_PDF_type

       ! Test 16
       ! Compare the original thlm to its recalculated value.
       !    | ( <thl>|_recalc - <thl> ) / <thl> |  <=  tol;
       ! which can be rewritten as:
       !    | <thl>|_recalc - <thl> |  <=  | <thl> | * tol.
       where ( abs( recalc_thlm - thlm ) &
               <= max( abs( thlm ), thl_tol ) * tol )
          l_pass_test_16 = .true.
       elsewhere
          l_pass_test_16 = .false.
       endwhere

       if ( any( .not. l_pass_test_16 ) ) then
          do idx = 1, nz, 1
             if ( .not. l_pass_test_16(idx) ) then
                write(fstderr,*) "Test 16 failed"
                !write(fstderr,*) "index = ", idx
                write(fstderr,*) "thlm = ", thlm(idx)
                write(fstderr,*) "recalc_thlm = ", recalc_thlm(idx)
                write(fstderr,*) ""
             endif ! .not. l_pass_test_16(idx)
          enddo ! idx = 1, nz, 1
       endif

       ! Test 17
       ! Compare the original thlp2 to its recalculated value.
       !    | ( <thl'^2>|_recalc - <thl'^2> ) / <thl'^2> |  <=  tol;
       ! which can be rewritten as:
       !    | <thl'^2>|_recalc - <thl'^2> |  <=  | <thl'^2> | * tol;
       ! and since <thl'^2> is always positive:
       !    | <thl'^2>|_recalc - <thl'^2> |  <=  <thl'^2> * tol.
       where ( abs( recalc_thlp2 - thlp2 ) <= max( thlp2, thl_tol**2 ) * tol )
          l_pass_test_17 = .true.
       elsewhere
          l_pass_test_17 = .false.
       endwhere

       if ( any( .not. l_pass_test_17 ) ) then
          do idx = 1, nz, 1
             if ( .not. l_pass_test_17(idx) ) then
                write(fstderr,*) "Test 17 failed"
                !write(fstderr,*) "index = ", idx
                write(fstderr,*) "thlp2 = ", thlp2(idx)
                write(fstderr,*) "recalc_thlp2 = ", recalc_thlp2(idx)
                write(fstderr,*) ""
             endif ! .not. l_pass_test_17(idx)
          enddo ! idx = 1, nz, 1
       endif

       if ( ( test_PDF_type == iiPDF_new ) &
            .or. ( test_PDF_type == iiPDF_new_hybrid ) &
            .or. ( test_PDF_type == iiPDF_LY93 ) ) then

          ! Test 18
          ! Compare the original thlp3 to its recalculated value.
          !    | ( <thl'^3>|_recalc - <thl'^3> ) / <thl'^3> |  <=  tol;
          ! which can be rewritten as:
          !    | <thl'^3>|_recalc - <thl'^3> |  <=  | <thl'^3> | * tol.
          where ( abs( recalc_thlp3 - thlp3 ) &
                  <= max( abs( thlp3 ), thl_tol**3 ) * tol )
             l_pass_test_18 = .true.
          elsewhere
             l_pass_test_18 = .false.
          endwhere

          if ( any( .not. l_pass_test_18 ) ) then
             do idx = 1, nz, 1
                if ( .not. l_pass_test_18(idx) ) then
                   write(fstderr,*) "Test 18 failed"
                   !write(fstderr,*) "index = ", idx
                   write(fstderr,*) "thlp3 = ", thlp3(idx)
                   write(fstderr,*) "recalc_thlp3 = ", recalc_thlp3(idx)
                   write(fstderr,*) ""
                endif ! .not. l_pass_test_18(idx)
             enddo ! idx = 1, nz, 1
          endif

          ! Test 19
          ! Compare the original Skthl to its recalculated value.
          !    | ( Skthl|_recalc - Skthl ) / Skthl |  <=  tol;
          ! which can be rewritten as:
          !    | Skthl|_recalc - Skthl |  <=  | Skthl | * tol.
          where ( abs( recalc_Skthl - Skthl ) &
                  <= max( abs( Skthl ), 1.0e-3_core_rknd ) * tol )
             l_pass_test_19 = .true.
          elsewhere
             l_pass_test_19 = .false.
          endwhere

          if ( any( .not. l_pass_test_19 ) ) then
             do idx = 1, nz, 1
                if ( .not. l_pass_test_19(idx) ) then
                   write(fstderr,*) "Test 19 failed"
                   !write(fstderr,*) "index = ", idx
                   write(fstderr,*) "Skthl = ", Skthl(idx)
                   write(fstderr,*) "recalc_Skthl = ", recalc_Skthl(idx)
                   write(fstderr,*) ""
                endif ! .not. l_pass_test_19(idx)
             enddo ! idx = 1, nz, 1
          endif

       else

          ! Automatically set the pass test flags to true for PDF types that
          ! don't necessarily preserve <thl'^3> and Skthl.
          l_pass_test_18 = .true.
          l_pass_test_19 = .true.

       endif ! test_PDF_type

       ! Test 20
       ! Check that sigma_w_1, sigma_w_2, sigma_rt_1, sigma_rt_2, sigma_thl_1,
       ! and sigma_thl_2 have a value of at least zero, and that mixture
       ! fraction has a value between 0 and 1.
       where ( ( sigma_w_1 >= zero ) .and. ( sigma_w_2 >= zero ) &
               .and. ( sigma_rt_1 >= zero ) .and. ( sigma_rt_2 >= zero ) &
               .and. ( sigma_thl_1 >= zero ) .and. ( sigma_thl_2 >= zero ) &
               .and. ( mixt_frac > zero ) .and. ( mixt_frac < one ) )
          l_pass_test_20 = .true.
       elsewhere
          l_pass_test_20 = .false.
       endwhere

       if ( any( .not. l_pass_test_20 ) ) then
          do idx = 1, nz, 1
             if ( .not. l_pass_test_20(idx) ) then
                write(fstderr,*) "Test 20 failed"
                !write(fstderr,*) "index = ", idx
                write(fstderr,*) "sigma_w_1 = ", sigma_w_1(idx)
                write(fstderr,*) "sigma_w_2 = ", sigma_w_2(idx)
                write(fstderr,*) "sigma_rt_1 = ", sigma_rt_1(idx)
                write(fstderr,*) "sigma_rt_2 = ", sigma_rt_2(idx)
                write(fstderr,*) "sigma_thl_1 = ", sigma_thl_1(idx)
                write(fstderr,*) "sigma_thl_2 = ", sigma_thl_2(idx)
                write(fstderr,*) "mixt_frac = ", mixt_frac(idx)
                write(fstderr,*) ""
             endif ! .not. l_pass_test_20(idx)
          enddo ! idx = 1, nz, 1
       endif

       where ( l_pass_test_8 .and. l_pass_test_9 &
               .and. l_pass_test_10 .and. l_pass_test_11 &
               .and. l_pass_test_12 .and. l_pass_test_13 &
               .and. l_pass_test_14 .and. l_pass_test_15 &
               .and. l_pass_test_16 .and. l_pass_test_17 &
               .and. l_pass_test_18 .and. l_pass_test_19 &
               .and. l_pass_test_20 )

          ! All tests pass
          num_failed_sets = num_failed_sets
          l_failed_sets = .false.

       elsewhere

          ! At least one test failed
          num_failed_sets = num_failed_sets + 1
          l_failed_sets = .true.

       endwhere

       if ( any( l_failed_sets ) ) then
          do idx = 1, nz, 1
             if ( l_failed_sets(idx) ) then
                write(fstderr,*) "At least one test or check for the full PDF "&
                                 // "failed for the following parameter set:  "
                write(fstderr,*) "PDF parameter set index = ", iter_param_sets
                write(fstderr,*) ""
             endif ! l_failed_sets(idx)
          enddo ! idx = 1, nz, 1
       endif

    enddo ! iter_param_sets = 1, num_param_sets, 1


    total_num_failed_sets = sum( num_failed_sets )

    ! Print results and return exit code.
    if ( total_num_failed_sets == 0 ) then
       write(fstdout,'(1x,A)') "Success!"
       write(fstdout,*) ""
       pdf_parameter_unit_tests = 0 ! Exit Code = 0, Success!
    else ! total_num_failed_sets > 0
       write(fstdout,'(1x,A,I4,A)') "There were ", total_num_failed_sets, &
                                    " failed parameter sets."
       write(fstdout,*) ""
       pdf_parameter_unit_tests = 1 ! Exit Code = 1, Fail
    endif ! total_num_failed_sets = 0


    return

  end function pdf_parameter_unit_tests

  !=============================================================================
  subroutine setter_var_tests( nz, wm, wp2, wp3, Skw,        & ! In
                               mu_w_1, mu_w_2, sigma_w_1,    & ! In
                               sigma_w_2, mixt_frac, tol,    & ! In
                               sigma_w_1_sqd, sigma_w_2_sqd, & ! In
                               l_pass_test_1, l_pass_test_2, & ! Out
                               l_pass_test_3, l_pass_test_4, & ! Out
                               l_pass_test_5, l_pass_test_6, & ! Out
                               l_check_mu_w_1_gte_mu_w_2     ) ! Out

    ! Description:
    ! Tests 1 through 5 as described above.  These tests are all applied to the
    ! "setter" variable, which is the variable that is used to set the mixture
    ! fraction.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,       & ! Variable(s)
        zero,      &
        w_tol,     &
        w_tol_sqd, &
        fstderr

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      wm,            & ! Mean of w (overall)                           [m/s]
      wp2,           & ! Variance of w (overall)                       [m^2/s^2]
      wp3,           & ! <w'^3>                                        [m^3/s^3]
      Skw,           & ! Skewness of w (overall)                       [-]
      mu_w_1,        & ! Mean of w (1st PDF component)                 [m/s]
      mu_w_2,        & ! Mean of w (2nd PDF component)                 [m/s]
      sigma_w_1,     & ! Standard deviation of w (1st PDF component)   [m/s]
      sigma_w_2,     & ! Standard deviation of w (2nd PDF component)   [m/s]
      mixt_frac,     & ! Mixture fraction                              [-]
      sigma_w_1_sqd, & ! Variance of w (1st PDF component)             [m^2/s^2]
      sigma_w_2_sqd    ! Variance of w (2nd PDF component)             [m^2/s^2]

    real( kind = core_rknd ), intent(in) :: &
      tol    ! Tolerance for acceptable numerical difference    [-]

    logical, intent(in) :: &
      l_check_mu_w_1_gte_mu_w_2    ! Flag to check whether mu_w_1 >= mu_w_2

    ! Output Variables
    logical, dimension(nz), intent(out) :: &
      l_pass_test_1, & ! Flag for passing test 1
      l_pass_test_2, & ! Flag for passing test 2
      l_pass_test_3, & ! Flag for passing test 3
      l_pass_test_4, & ! Flag for passing test 4
      l_pass_test_5, & ! Flag for passing test 5
      l_pass_test_6    ! Flag for passing test 6

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) :: &
      recalc_wm,  & ! Recalculation of <w> using PDF parameters          [m/s]
      recalc_wp2, & ! Recalculation of <w'^2> using PDF parameters   [m^2/s^2]
      recalc_wp3, & ! Recalculation of <w'^3> using PDF parameters   [m^3/s^3]
      recalc_Skw    ! Recalculation of Skw using PDF parameters            [-]

    integer :: idx    ! Loop index


    ! Recalculate the values of <wm>, <w'^2>, <w'^3>, and Skw using the PDF
    ! parameters that are output from the subroutine.
    call recalc_single_var( nz, wm, mu_w_1, mu_w_2, sigma_w_1, & ! In
                            sigma_w_2, mixt_frac,              & ! In
                            recalc_wm, recalc_wp2,             & ! Out
                            recalc_wp3, recalc_Skw             ) ! Out

    ! Test 1
    ! Compare the original wm to its recalculated value.
    !    | ( <w>|_recalc - <w> ) / <w> |  <=  tol;
    ! which can be rewritten as:
    !    | <w>|_recalc - <w> |  <=  | <w> | * tol.
    where ( abs( recalc_wm - wm ) <= max( abs( wm ), w_tol ) * tol )
       l_pass_test_1 = .true.
    elsewhere
       l_pass_test_1 = .false.
    endwhere

    if ( any( .not. l_pass_test_1 ) ) then
       do idx = 1, nz, 1
          if ( .not. l_pass_test_1(idx) ) then
             write(fstderr,*) "Test 1 failed"
             !write(fstderr,*) "index = ", idx
             write(fstderr,*) "wm = ", wm(idx)
             write(fstderr,*) "recalc_wm = ", recalc_wm(idx)
             write(fstderr,*) ""
          endif ! .not. l_pass_test_1(idx)
       enddo ! idx = 1, nz, 1
    endif

    ! Test 2
    ! Compare the original wp2 to its recalculated value.
    !    | ( <w'^2>|_recalc - <w'^2> ) / <w'^2> |  <=  tol;
    ! which can be rewritten as:
    !    | <w'^2>|_recalc - <w'^2> |  <=  | <w'^2> | * tol;
    ! and since <w'^2> is always positive:
    !    | <w'^2>|_recalc - <w'^2> |  <=  <w'^2> * tol.
    where ( abs( recalc_wp2 - wp2 ) <= max( wp2, w_tol_sqd ) * tol )
       l_pass_test_2 = .true.
    elsewhere
       l_pass_test_2 = .false.
    endwhere

    if ( any( .not. l_pass_test_2 ) ) then
       do idx = 1, nz, 1
          if ( .not. l_pass_test_2(idx) ) then
             write(fstderr,*) "Test 2 failed"
             !write(fstderr,*) "index = ", idx
             write(fstderr,*) "wp2 = ", wp2(idx)
             write(fstderr,*) "recalc_wp2 = ", recalc_wp2(idx)
             write(fstderr,*) ""
          endif ! .not. l_pass_test_2(idx)
       enddo ! idx = 1, nz, 1
    endif

    ! Test 3
    ! Compare the original wp3 to its recalculated value.
    !    | ( <w'^3>|_recalc - <w'^3> ) / <w'^3> |  <=  tol;
    ! which can be rewritten as:
    !    | <w'^3>|_recalc - <w'^3> |  <=  | <w'^3> | * tol.
    where ( abs( recalc_wp3 - wp3 ) <= max( abs( wp3 ), w_tol**3 ) * tol )
       l_pass_test_3 = .true.
    elsewhere
       l_pass_test_3 = .false.
    endwhere

    if ( any( .not. l_pass_test_3 ) ) then
       do idx = 1, nz, 1
          if ( .not. l_pass_test_3(idx) ) then
             write(fstderr,*) "Test 3 failed"
             !write(fstderr,*) "index = ", idx
             write(fstderr,*) "wp3 = ", wp3(idx)
             write(fstderr,*) "recalc_wp3 = ", recalc_wp3(idx)
             write(fstderr,*) ""
          endif ! .not. l_pass_test_3(idx)
       enddo ! idx = 1, nz, 1
    endif

    ! Test 4
    ! Compare the original Skw to its recalculated value.
    !    | ( Skw|_recalc - Skw ) / Skw |  <=  tol;
    ! which can be rewritten as:
    !    | Skw|_recalc - Skw |  <=  | Skw | * tol.
    where ( abs( recalc_Skw - Skw ) &
            <= max( abs( Skw ), 1.0e-3_core_rknd ) * tol )
       l_pass_test_4 = .true.
    elsewhere
       l_pass_test_4 = .false.
    endwhere

    if ( any( .not. l_pass_test_4 ) ) then
       do idx = 1, nz, 1
          if ( .not. l_pass_test_4(idx) ) then
             write(fstderr,*) "Test 4 failed"
             !write(fstderr,*) "index = ", idx
             write(fstderr,*) "Skw = ", Skw(idx)
             write(fstderr,*) "recalc_Skw = ", recalc_Skw(idx)
             write(fstderr,*) ""
          endif ! .not. l_pass_test_4(idx)
       enddo ! idx = 1, nz, 1
    endif

    ! Test 5
    ! Check that sigma_w_1^2 and sigma_w_2^2 have a value of at least zero, and
    ! that mixture fraction has a value between 0 and 1.
    where ( ( sigma_w_1_sqd >= zero ) .and. ( sigma_w_2_sqd >= zero ) &
            .and. ( mixt_frac > zero ) .and. ( mixt_frac < one ) )
       l_pass_test_5 = .true.
    elsewhere
       l_pass_test_5 = .false.
    endwhere

    if ( any( .not. l_pass_test_5 ) ) then
       do idx = 1, nz, 1
          if ( .not. l_pass_test_5(idx) ) then
             write(fstderr,*) "Test 5 failed"
             !write(fstderr,*) "index = ", idx
             write(fstderr,*) "sigma_w_1^2 = ", sigma_w_1_sqd(idx)
             write(fstderr,*) "sigma_w_2^2 = ", sigma_w_2_sqd(idx)
             write(fstderr,*) "mixt_frac = ", mixt_frac(idx)
             write(fstderr,*) ""
          endif ! .not. l_pass_test_5(idx)
       enddo ! idx = 1, nz, 1
    endif

    if ( l_check_mu_w_1_gte_mu_w_2 ) then

       ! Test 6
       ! Check that mu_w_1 >= mu_w_2
       where ( mu_w_1 >= mu_w_2 )
          l_pass_test_6 = .true.
       elsewhere
          l_pass_test_6 = .false.
       endwhere

       if ( any( .not. l_pass_test_6 ) ) then
          do idx = 1, nz, 1
             if ( .not. l_pass_test_6(idx) ) then
                write(fstderr,*) "Test 6 failed"
                !write(fstderr,*) "index = ", idx
                write(fstderr,*) "mu_w_1 = ", mu_w_1(idx)
                write(fstderr,*) "mu_w_2 = ", mu_w_2(idx)
                write(fstderr,*) ""
             endif ! .not. l_pass_test_6(idx)
          enddo ! idx = 1, nz, 1
       endif

    else

       l_pass_test_6 = .true.

    endif


    return

  end subroutine setter_var_tests

  !=============================================================================
  subroutine recalc_single_var( nz, xm, mu_x_1, mu_x_2, sigma_x_1, & ! In
                                sigma_x_2, mixt_frac,              & ! In
                                recalc_xm, recalc_xp2,             & ! Out
                                recalc_xp3, recalc_Skx             ) ! Out

    use constants_clubb, only: &
        three, & ! Variable(s)
        one

    use clubb_precision, only: &
        core_rknd

    implicit none

    integer, intent(in) :: &
      nz

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      xm,        & ! Mean of x (overall)                  [units vary]
      mu_x_1,    & ! Mean of x (1st PDF component)        [units vary]
      mu_x_2,    & ! Mean of x (2nd PDF component)        [units vary]
      sigma_x_1, & ! Variance of x (1st PDF component)    [units vary]
      sigma_x_2, & ! Variance of x (2nd PDF component)    [units vary]
      mixt_frac    ! Mixture fraction                     [-]

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
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

    recalc_Skx = recalc_xp3 / max( recalc_xp2**1.5, epsilon( recalc_xp2 ) )


    return

  end subroutine recalc_single_var

  !=============================================================================

end module pdf_parameter_tests
