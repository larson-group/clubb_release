!-----------------------------------------------------------------------
! $Id$

module error

! Description:

!   subroutine tuner_init: reads in namelists /stats/, /cases/,
!   /clubb_params_nl/, & /variance/ from 'error.in'
!   It then uses them to setup the initial param_vals_matrix of independent
!   variables, i.e. the CLUBB constants, and allocate the runtime arrays
!   for each of the model runs and each of the variables

!   function min_les_clubb_diff:  A driver for the CLUBB program/module.
!   Calls run_clubb, reads in les & CLUBB results from GRADS files, and
!   calculates the average difference between the two over all z-levels.

!   subroutine write_results :
!   Prints the results of tuning to the terminal or to a file.

!   subroutine output_nml_tuner :
!   Generates the error.in file using the current constants.

!   subroutine output_nml_standalone :
!   Generates the standalone.in file using the current constants.
!   Standalone CLUBB is only configured to run a single model, and so only
!   the first model is used to make the namelist.

!-----------------------------------------------------------------------
  use parameter_indices, only: nparams ! Variable(s)

  use parameters_tunable, only:  &
    set_default_parameters, & ! Procedure(s)
    read_parameters, & ! Procedure(s)
    read_param_minmax, &
    read_param_constraints

  use parameters_tunable, only: &
    params_list  ! Variable(s)

  use mt95, only: genrand_real ! Constant

  use constants_clubb, only: eps

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  ! Constant Parameters
  integer, parameter, private ::  & 
    max_run = 12,  & ! Maximum model runs the tuner can handle at a time
    max_variables = 32 ! This number / 2 is maximum variables to tune for

  integer, parameter, public :: &
    iamoeba = 0, & ! Numerical Recipes downhill simplex
    iamebsa = 1, & ! Numerical Recipes simulated annealing
    iesa    = 2, & ! ESA tuning type
    iflags  = 3, & ! Model flags "tuner"
    iploops = 4    ! Parameter loop tuning (loop over range of parameter vals.)

  ! Variables
  integer, public :: &
    ndim ! Number of variables, e.g. rcm, to be tuned. Dimension of the init simplex

  ! inv_count is a modular counter [1-3] used to determine
  ! which file to output to if l_stdout_on_invalid is true.
  integer, private ::  & 
    inv_count = 0

  !-----------------------------------------------------------------------


  real( kind = core_rknd ), public ::  & 
    f_tol = 1e-5_core_rknd,                 & ! The precision to tune for
    anneal_temp = 100._core_rknd,           & ! Initial temperature for the SA algorithm
    max_final_temp = 1._core_rknd,          & ! Maxmimum final temperature for the SA algorithm
    stp_adjst_shift_in = 0.5_core_rknd,     & ! Shift parameter for step size calculation
    stp_adjst_factor_in = 1._core_rknd        ! Linear coefficient for step size calculation

  integer, public :: & 
    anneal_iter = 0, &    ! Number of annealing iterations to perform
    tune_type = iesa, &   ! Toggle for downhill simplex of simulated annealing
    c_total = 0, &        ! Total number of simulation cases to tune over
    v_total = 0, &        ! Total number of variables to tune over
    prescribed_rand_seed = 1, & ! Default value for prescribed random seed for esa methods
                                ! Value read from stats nml in input_misc/tuner/error*.in
    max_iters_in = 2000         ! Max number of iteration steps for tuner. Read from error*.in

  logical, public :: & 
    l_results_stdout = .false.,    & ! Whether to print tuning results to the terminal
    l_results_file = .false.,      & ! Whether to generate a new error.in based on
                                     ! the new tuning constants
    l_stdout_on_invalid = .false., & ! Generate a new error.in when the simulation crashes
    l_keep_params_equal = .false., & ! Whether to keep parameters like C1, C1b,
                                     ! etc. equal throughout tuning run
    l_use_prescribed_rand_seed = .false.    ! Whether to use a fixed random seed for
                                            ! esa methods. Value in prescribed_rand_seed

  logical, parameter, public :: &
    l_save_tuning_run = .true.  ! If true, writes the results of the tuning run to a file

  character(len=50), public :: &
    tuning_filename = '' ! File where results of tuning run are
                         ! written if l_save_tuning_run = .true.

  integer, parameter, public :: &
    file_unit = 15  ! File unit number connected with tuning_filename


  character(len=10), dimension(:), allocatable, private ::  &
    hoc_v,  & ! Variables in CLUBB GrADS files
    les_v  ! Variables in LES GrADS files

  integer, dimension(:,:), allocatable, private ::  & 
    timestep_intvls ! Time intervals

  ! Additions for using imposed weights as scaling factors
  logical :: l_initialize_sigma = .true.

  real( kind = core_rknd ), dimension(:,:), allocatable, private ::  & 
    err_terms, & 
    invsigma2, & 
    min_err_terms, & 
    init_err_terms

  real( kind = core_rknd ), dimension(:), allocatable, private ::  & 
    weight_case, weight_var

  ! End additions for using imposed weights

  integer, dimension(:), allocatable, private :: & 
    z_i, & ! Initial z level for tuning purposes
    z_f    ! Final z level for tuning purposes

  character(len=150), dimension(:), allocatable, private ::  & 
    run_file,        & ! Model run files
    hoc_stats_file,  & ! Model GrADS files
    les_stats_file  ! Model GrADS files

  ! Various Variables for returning results
  integer, public :: &
    iter = 0 ! Total number of iterations amoeba spent calculating optimal values

  real( kind = core_rknd ), public :: & 
    init_err    = -999._core_rknd,  & ! Error for the initial constants
    min_err     = -999._core_rknd,  & ! The lowest the minimization algorithm could go
    min_err_old = -999._core_rknd     ! Same as above, used to find min_err_terms

  real( kind = core_rknd ), dimension(nparams), private :: & 
    params = -999._core_rknd  ! Vector of all CLUBB tunable parameter values

  integer, dimension(nparams), private :: & 
    params_index = 0  ! Index of the params elements that are used in the simplex

  character(len=28), dimension(nparams), private :: &
    param_constraints  ! Array of parameters meant to be kept equal to each other
                       ! throughout tuning run

  real( kind = core_rknd ), allocatable, dimension(:,:), public ::  & 
    param_vals_matrix ! Holds 2D simplex the CLUBB constant parameters
    ! The first row contains the initial values of the tunable parameters.
    ! The remaining rows contain random perturbations of those parameter values.

  real( kind = core_rknd ), allocatable, dimension(:,:), public :: &
    param_vals_minmax  ! Amount to vary each respec. constant by

  real( kind = core_rknd ), allocatable, dimension(:), public :: &
    cost_fnc_vector    ! cache of differences between the LES and CLUBB

  real(kind=genrand_real), allocatable, dimension(:), private ::  & 
    rand_vect ! A vector of random reals for initializing the x array

  logical, public, dimension(:,:), allocatable :: &
    model_flags_array ! Flags we're checking for


  ! Procedures
  public :: tuner_init, min_les_clubb_diff, write_results, & 
    output_nml_standalone, output_nml_tuner

  private ! Default Scope

  contains

  !-----------------------------------------------------------------------
  subroutine tuner_init( l_read_files )

! Description:
!   Initializes param_vals_matrix with constants from error.in
!   Allocates arrays for cases and tuning variables.
!   Initializes grads file names to read in.

! References:
!   None
!-----------------------------------------------------------------------
    use constants_clubb, only: fstderr ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constant

    use text_writer, only: write_text ! Subroutine(s)

    use mt95, only: genrand_real1 ! Procedure

    ! These variables are only used for the parameter loop tuner.
    ! They need to be changed to use different tunable parameters.
    use parameter_indices, only: &
        islope_coef_spread_DG_means_w, & ! Variable(s)
        ipdf_component_stdev_factor_w, &
        icoef_spread_DG_means_rt, &
        icoef_spread_DG_means_thl

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Constant Variables

    integer, parameter ::  & 
      max_times = 50 ! max number of timesteps to compare (arbitrary)

    character(len=8), parameter :: & 
      filename = "error.in"

    ! Input  Variables

    ! Determines whether to read in the namelists and do the
    ! initial allocation of arrays
    logical, intent(in) ::  & 
      l_read_files

    ! Local Variables

    real( kind = core_rknd ), dimension(2,nparams) :: params_minmax_all_vals

    integer :: i, j, & ! looping variables
               iunit ! file unit

    real( kind = core_rknd ) :: & 
      C1, C1b, C1c, C2rt, C2thl, C2rtthl, & 
      C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, & 
      C7, C7b, C7c, C8, C8b, C10, & 
      C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, &
      C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, & 
      C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
      c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8,  & 
      c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, &
      slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, &
      coef_spread_DG_means_rt, coef_spread_DG_means_thl, &
      gamma_coef, gamma_coefb, gamma_coefc, mu, beta, lmin_coef, &
      omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
      lambda0_stability_coef, mult_coef, taumin, taumax, Lscale_mu_coef, &
      Lscale_pert_coef, alpha_corr, Skw_denom_coef, c_K10, c_K10h, &
      thlp2_rad_coef, thlp2_rad_cloud_frac_thresh, up2_sfc_coef, &
      Skw_max_mag, xp3_coef_base, xp3_coef_slope, altitude_threshold, &
      rtp2_clip_coef, C_invrs_tau_bkgnd, C_invrs_tau_sfc, &
      C_invrs_tau_shear, C_invrs_tau_N2, C_invrs_tau_N2_wp2, &
      C_invrs_tau_N2_xp2, C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3, &
      C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh, &
      Cx_min, Cx_max, Richardson_num_min, Richardson_num_max, a3_coef_min, a_const, bv_efold

    !-----------------------------------------------------------------------

    ! Namelist vars for determining which variable to tune for:
    ! F95 does not allow allocatable arrays in namelists. These variables are used
    ! for reading the error file, then assigned to their non-namelist equivalents below.

    !  timestep_intvls_nl:   Order pairs of time intervals to analyze
    !  z_i_nl:    initial z-level to begin reading in for tuning
    !  z_f_nl:    final z-level to end reading in for tuning

    integer, dimension(max_run, max_times):: timestep_intvls_nl

    integer, dimension(max_run) :: z_i_nl, z_f_nl

    ! Addition to use imposed weights as scaling factors
    real( kind = core_rknd ), dimension(max_run) :: weight_case_nl

    real( kind = core_rknd ), dimension(max_variables) :: weight_var_nl

    character(len=150), dimension(max_run) ::  & 
      run_file_nl, & 
      hoc_stats_file_nl, & 
      les_stats_file_nl

    character(len=10), dimension(max_variables) :: &
      t_variables ! List of variables to be read from the GrADS output

    ! Variables for parameter loops (ploops) tuner
    integer :: &
      num_var1, & ! Number of parameter values for variable 1
      num_var2, & ! Number of parameter values for variable 2
      num_var3, & ! Number of parameter values for variable 3
      num_var4, & ! Number of parameter values for variable 4
      num_runs    ! Total number of tuning runs

    integer :: &
      idx_var1, & ! Index of parameter value for variable 1
      idx_var2, & ! Index of parameter value for variable 2
      idx_var3, & ! Index of parameter value for variable 3
      idx_var4, & ! Index of parameter value for variable 4
      mod1,     & ! Remainder used to calculate idx_var1
      mod2,     & ! Remainder used to calculate idx_var2
      mod3,     & ! Remainder used to calculate idx_var3
      mod4        ! Remainder used to calculate idx_var4

    real( kind = core_rknd ), dimension(:), allocatable :: &
      tune_var1, & ! Values of variable 1
      tune_var2, & ! Values of variable 2
      tune_var3, & ! Values of variable 3
      tune_var4    ! Values of variable 4

    ! Namelists read from error.in
    namelist /stats/  & 
      f_tol, tune_type, anneal_temp, max_final_temp, anneal_iter, & 
      l_results_stdout, l_results_file, l_stdout_on_invalid, &
      l_keep_params_equal, &
      t_variables, weight_var_nl, stp_adjst_shift_in, stp_adjst_factor_in, &
      l_use_prescribed_rand_seed, prescribed_rand_seed, max_iters_in

    namelist /cases/  & 
      les_stats_file_nl, hoc_stats_file_nl, & 
      run_file_nl, z_i_nl, z_f_nl, timestep_intvls_nl, weight_case_nl

    ! Reset iteration counter (set by amoeba)
    iter = 0

    ! Reset invalid run counter
    inv_count = 0

    ! Set the default tunable parameter values
    call set_default_parameters( &
               C1, C1b, C1c, C2rt, C2thl, C2rtthl, &
               C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, &
               C6thl, C6thlb, C6thlc, C7, C7b, C7c, C8, C8b, C10, &
               C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, &
               C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, &
               C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
               c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8, &
               c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, &
               slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, &
               coef_spread_DG_means_rt, coef_spread_DG_means_thl, &
               gamma_coef, gamma_coefb, gamma_coefc, mu, beta, lmin_coef, &
               omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
               lambda0_stability_coef, mult_coef, taumin, taumax, &
               Lscale_mu_coef, Lscale_pert_coef, alpha_corr, &
               Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, &
               thlp2_rad_cloud_frac_thresh, up2_sfc_coef, &
               Skw_max_mag, xp3_coef_base, xp3_coef_slope, &
               altitude_threshold, rtp2_clip_coef, C_invrs_tau_bkgnd, &
               C_invrs_tau_sfc, C_invrs_tau_shear, C_invrs_tau_N2, & 
               C_invrs_tau_N2_wp2, C_invrs_tau_N2_xp2, &
               C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3, &
               C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh, &
               Cx_min, Cx_max, Richardson_num_min, &
               Richardson_num_max, a3_coef_min, a_const, bv_efold )

    ! Re-read namelists if requested
    if ( l_read_files ) then

      ! Initialize all compile time arrays to zero
      timestep_intvls_nl = 0

      ! Imposed weights as scaling factors
      weight_case_nl = 0.0_core_rknd
      weight_var_nl  = 0.0_core_rknd

      z_i_nl = 0
      z_f_nl = 0

      ! Initialize variable names to spaces
      t_variables(1:max_variables)  = "          "

      iunit = 10

      ! Open our namelist input file
      open(unit=iunit, file=filename, status='old')

      ! Determine which files to read data from based on namelist
      read(unit=iunit, nml=stats)

      ! Read in the models to be run
      read(unit=iunit, nml=cases)

      ! Close our input namelist file
      close(unit=iunit)

      ! Read in initial constant values
      call read_parameters( iunit, filename, &
                            C1, C1b, C1c, C2rt, C2thl, C2rtthl, &
                            C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, &
                            C6thl, C6thlb, C6thlc, C7, C7b, C7c, C8, C8b, C10, &
                            C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, &
                            C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, &
                            C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
                            c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8, &
                            c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, &
                            slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, &
                            coef_spread_DG_means_rt, coef_spread_DG_means_thl, &
                            gamma_coef, gamma_coefb, gamma_coefc, mu, beta, lmin_coef, &
                            omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
                            lambda0_stability_coef, mult_coef, taumin, taumax, &
                            Lscale_mu_coef, Lscale_pert_coef, alpha_corr, &
                            Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, &
                            thlp2_rad_cloud_frac_thresh, up2_sfc_coef, &
                            Skw_max_mag, xp3_coef_base, xp3_coef_slope, &
                            altitude_threshold, rtp2_clip_coef, C_invrs_tau_bkgnd, &
                            C_invrs_tau_sfc, C_invrs_tau_shear, C_invrs_tau_N2, & 
                            C_invrs_tau_N2_wp2, C_invrs_tau_N2_xp2, &
                            C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3, &
                            C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh, &
                            Cx_min, Cx_max, Richardson_num_min, &
                            Richardson_num_max, a3_coef_min, a_const, bv_efold, &
                            params )

      ! Allocate the arrays for the tuning variables

      do i = 1, max_variables, 2 ! 1, 3, 5, 7
        if (t_variables(i) == "          ") exit
        v_total = (i + 1) / 2
      end do

      allocate( hoc_v(v_total), les_v(v_total) )

      ! Allocate the arrays for the run cases
      do i=1, max_run
        if (z_f_nl(i) == 0 ) exit
        c_total = i
      end do

      allocate( & 
        z_i(c_total), z_f(c_total), timestep_intvls(c_total, max_times),  &
        run_file(c_total),  & 
        les_stats_file(c_total), hoc_stats_file(c_total) )

      allocate( & 
        err_terms(c_total, v_total), & 
        min_err_terms(c_total, v_total), & 
        init_err_terms(c_total, v_total), & 
        invsigma2(c_total, v_total), & 
        weight_case(c_total), weight_var(v_total) )

      ! Transfer the variable numbers to hoc_v and les_v
      do i=1, v_total
        hoc_v(i) = t_variables(i*2 - 1)
        les_v(i) = t_variables(i*2)
      end do

      ! Transfer the case information to run-time arrays
      do i = 1, c_total
        z_i(i)              = z_i_nl(i)
        z_f(i)              = z_f_nl(i)
        les_stats_file(i)   = les_stats_file_nl(i)
        hoc_stats_file(i)   = hoc_stats_file_nl(i)
        run_file(i)         = run_file_nl(i)
        timestep_intvls(i,1:max_times) = timestep_intvls_nl(i,1:max_times)
      end do

      ! Use imposed weights as scaling factors
      weight_case(1:c_total) = weight_case_nl(1:c_total)
      weight_var(1:v_total)  = weight_var_nl(1:v_total)

      ! Setup the simplex

      if ( tune_type /= iploops ) then
         call read_param_minmax( iunit, filename, params_index,  &
                                 params_minmax_all_vals, ndim )
      else ! tune_type == iploops

         ! Number of tunable parameters that are looped over.
         ndim = 4
         ! Initialize to 0.
         params_index(1:nparams) = 0
         ! Set the indices of the variables to be tuned.
         params_index(1) = islope_coef_spread_DG_means_w
         params_index(2) = ipdf_component_stdev_factor_w
         params_index(3) = icoef_spread_DG_means_rt
         params_index(4) = icoef_spread_DG_means_thl

         ! Set to 0 (params_minmax_all_vals doesn't matter for parameter loops tuning).
         params_minmax_all_vals(1:2,1:nparams) = 0.0_core_rknd

         ! Number of parameter values to be looped over for each variable.
         ! Must be at least 1.
         ! var1 corresponds to the variable in params_index(1).
         ! var2 corresponds to the variable in params_index(2).
         ! var3 corresponds to the variable in params_index(3).
         ! var4 corresponds to the variable in params_index(4).
         num_var1 = 10
         num_var2 = 10
         num_var3 = 10
         num_var4 = 10

         ! Total number of tuning runs.
         num_runs = num_var1 * num_var2 * num_var3 * num_var4 + 1

         ! Allocate tunable variables.
         allocate( tune_var1(num_var1), tune_var2(num_var2), &
                   tune_var3(num_var3), tune_var4(num_var4) )

         ! Initialize variables
         tune_var1 = 0.0_core_rknd
         tune_var2 = 0.0_core_rknd
         tune_var3 = 0.0_core_rknd
         tune_var4 = 0.0_core_rknd

         ! Set parameter values for each variable.
         ! The number of values for each variable must equal the num_var value
         ! above.
         ! var1 corresponds to the variable in params_index(1).
         ! var2 corresponds to the variable in params_index(2).
         ! var3 corresponds to the variable in params_index(3).
         ! var4 corresponds to the variable in params_index(4).
         tune_var1 = (/ 0.5, 1.0, 1.5, 2.0, 3.0, &
                        4.0, 5.0, 7.5, 20.0, 45.0 /)
         tune_var2 = (/ 0.25, 0.5, 1.0, 2.0, 3.0, &
                        4.0, 5.0, 6.5, 8.0, 10.0 /)
         tune_var3 = (/ 0.05, 0.15, 0.25, 0.35, 0.45, &
                        0.55, 0.65, 0.75, 0.85, 0.95 /)
         tune_var4 = (/ 0.05, 0.15, 0.25, 0.35, 0.45, &
                        0.55, 0.65, 0.75, 0.85, 0.95 /)

      endif

      if ( ndim == 0 .and. tune_type /= iflags ) then
        write(fstderr,*) "You must vary at least one parameter"
        error stop
      end if

      ! Read in the parameter constraint namelist and do some simple error
      ! checks.  param_constraints is an array whose indices designate the
      ! tunable parameters, in the same order as params_list.  Its values are
      ! strings which should be the names of other tunable parameters. Each
      ! index-string pair represents equality between two parameters, 
      ! such as C1b (index) = C1 (string).
      if ( l_keep_params_equal ) then
        call read_param_constraints( iunit, filename, param_constraints )

        ! Error check to ensure that if l_keep_params_equal = .true. then at
        ! least one parameter is set equal to another.
        if ( all( param_constraints == "" ) ) then
          write(fstderr,*) "Check parameter_constraints namelist: l_keep_params_equal = " // &
                           ".true. but no parameters are set equal to each other."
          error stop
        end if

        ! Two more error checks:
        do i = 1, nparams, 1
          if ( len(trim(param_constraints(i))) > 0 ) then
            ! First check to make sure that parameters are only set equal to others
            ! that are actually being tuned.  This makes the output process simpler at
            ! the end of the tuning run, avoids user confusion, and should also avoid
            ! some unnecessary processing.
            if ( .not. any( param_constraints(i) == params_list(params_index(1:ndim)) ) ) then
              write(fstderr,*) "Check parameter_constraints namelist: a parameter " // &
                "has been set equal to " // trim(param_constraints(i)) // ", but " // &
                trim(param_constraints(i)) // " is not currently being tuned."
              error stop
            end if
            ! Second, check to make sure that if a parameter is set equal to another,
            ! the first parameter is not also being actively tuned. E.g. if
            ! parameter A is set equal to parameter B, parameter A should not be tuned.
            if ( any( params_index == i ) ) then
              write(*,*) "Check parameter_constraints namelist: " // trim(params_list(i)) // &
                " is currently being tuned but is also set equal to another parameter."
              error stop
            end if
          end if
        end do
      end if

      if ( tune_type == iamoeba .or. tune_type == iamebsa ) then
        ! Numerical recipes simulated annealing or downhill simplex
        allocate( rand_vect(ndim), param_vals_matrix(ndim+1,ndim), & 
                  param_vals_minmax(1,ndim), cost_fnc_vector(ndim+1) )
        param_vals_minmax(1,1:ndim) = params_minmax_all_vals(2,params_index(1:ndim))

      elseif ( tune_type == iploops ) then

        allocate( param_vals_matrix(num_runs,ndim), & 
                  param_vals_minmax(1,ndim), cost_fnc_vector(num_runs) )
        param_vals_minmax(1,1:ndim) = params_minmax_all_vals(2,params_index(1:ndim))

      else ! ESA algorithm or model flags
        allocate( param_vals_matrix(1,ndim),param_vals_minmax(2,ndim), cost_fnc_vector(1) )
        param_vals_minmax(1:2,1:ndim) = params_minmax_all_vals(1:2,params_index(1:ndim))
      end if

      ! Copy tunable parameter values into the first row of the simplex
      param_vals_matrix(1,1:ndim) = params(params_index(1:ndim))

      ! Attempt to generate a pseudo-random seed using a file
      ! generated from /dev/random.  File is an ASCII text file
      ! and can be edited manually.
      call read_random_seed( "rand_seed.dat" )

    end if  ! l_read_files
    !-----------------------------------------------------------------------

    if ( tune_type == iamoeba .or. tune_type == iamebsa ) then
      ! Fill in the remaining values of the array by varying the initial
      ! vector (i.e. the first row of the array) by a small random perturbation
      do j = 1, ndim

        call genrand_real1( rand_vect(1:ndim) )

        do i = 2, ndim+1, 1
            param_vals_matrix(i,j) = param_vals_matrix(1,j)* &
               (1.0_core_rknd+((real(i, kind = core_rknd)-1._core_rknd)/ &
               real(ndim, kind = core_rknd)*0.5_core_rknd))
        end do ! i..ndim+1

      end do ! j..ndim

    elseif ( tune_type == iploops ) then

       do i = 2, num_runs, 1

          ! Calculate the parameter value index for variable 1.
          if ( mod( i - 1, num_var2 * num_var3 * num_var4 ) > 0 ) then
             mod1 = mod( ( ( i - 1 ) / ( num_var2 * num_var3 * num_var4 ) ) &
                         + 1, num_var1 )
          else
             mod1 = mod( ( ( i - 1 ) / ( num_var2 * num_var3 * num_var4 ) ), &
                         num_var1 )
          endif

          if ( mod1 > 0 ) then
             idx_var1 = mod1
          else
             idx_var1 = num_var1
          endif

          ! Calculate the parameter value index for variable 2.
          if ( mod( i - 1, num_var3 * num_var4 ) > 0 ) then
             mod2 = mod( ( ( i - 1 ) / ( num_var3 * num_var4 ) ) + 1, num_var2 )
          else
             mod2 = mod( ( ( i - 1 ) / ( num_var3 * num_var4 ) ), num_var2 )
          endif

          if ( mod2 > 0 ) then
             idx_var2 = mod2
          else
             idx_var2 = num_var2
          endif

          ! Calculate the parameter value index for variable 3.
          if ( mod( i - 1, num_var4 ) > 0 ) then
             mod3 = mod( ( ( i - 1 ) / num_var4 ) + 1, num_var3 )
          else
             mod3 = mod( ( ( i - 1 ) / num_var4 ), num_var3 )
          endif
          
          if ( mod3 > 0 ) then
             idx_var3 = mod3
          else
             idx_var3 = num_var3
          endif

          ! Calculate the parameter value index for variable 4.
          mod4 = mod( i - 1, num_var4 )
          if ( mod4 > 0 ) then
             idx_var4 = mod4
          else
             idx_var4 = num_var4
          endif

          param_vals_matrix(i,1) = tune_var1(idx_var1)
          param_vals_matrix(i,2) = tune_var2(idx_var2)
          param_vals_matrix(i,3) = tune_var3(idx_var3)
          param_vals_matrix(i,4) = tune_var4(idx_var4)

       enddo ! i = 2, num_runs, 1

    end if


    ! First call is used to initialize weights

    l_initialize_sigma = .true.
    cost_fnc_vector(1) =  real(min_les_clubb_diff( &
        real(param_vals_matrix(1,1:ndim)) ), kind = core_rknd)
    l_initialize_sigma = .false.

    ! Note: min_les_clubb_diff is written to deal with undefined and
    ! invalid values for variations on the initial vector, but that
    ! algorithm relies on the initial vector being valid.

    if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then
          write(fstderr,*) "Initial variable values must be valid."
          error stop
        end if
    end if

    ! Save initial error

    init_err = cost_fnc_vector(1)
    min_err_old = init_err
    init_err_terms = err_terms

    ! Other initialization runs

    if ( tune_type == iamoeba .or. tune_type == iamebsa ) then
      ! Initialize the 'y' vector for amoeba
      ! This is done by calling min_les_clubb_diff with the initial vector
      do i = 2, ndim+1, 1
        cost_fnc_vector(i) =  & 
             real(min_les_clubb_diff( real(param_vals_matrix(i,1:ndim))), kind = core_rknd )
      end do
    elseif ( tune_type == iploops ) then
       do i = 2, num_runs, 1
           cost_fnc_vector(i) &
           = real( min_les_clubb_diff( real( param_vals_matrix(i,1:ndim) ) ), &
                   kind = core_rknd )
       enddo
    end if

    ! Save tuning results in file if specified
    if( l_save_tuning_run ) open(unit=file_unit, file=tuning_filename, &
      action="write", position='append')
    call write_text( "cost_fnc_vector:", l_save_tuning_run, file_unit )
    call write_text( '', cost_fnc_vector, l_save_tuning_run, file_unit )
    if( l_save_tuning_run ) close(unit=file_unit)

    return
  end subroutine tuner_init

  !-----------------------------------------------------------------------
  real function min_les_clubb_diff( param_vals_vector_r4 )

! Description:
!   Function that returns the sum of the error between the dependent
!   variable (i.e. the variable we want to match) in each of the models

! References:
!   _Numerical Recipes in Fortran 77_ P.402-406 (Description)
!   _Numerical Recipes in Fortran 90_ source code (Routine)
!-----------------------------------------------------------------------

    use clubb_driver, only: run_clubb ! Procedure(s)

    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_no_error,              & ! Constant
        clubb_fatal_error

    use stat_file_utils, only: &
      stat_file_num_vertical_levels, & ! Procedure(s)
      stat_file_vertical_levels, &
      stat_file_average_interval

    use numerical_check, only: is_nan_2d ! Procedure(s)

    use text_writer, only: write_text ! Subroutine(s)

    use constants_clubb, only: fstderr ! Constant(s)

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use stat_file_module, only: &
      stat_file ! Type(s)

#ifdef NETCDF
    use input_netcdf, only: &
      open_netcdf_read, & ! Procedure(s)
      close_netcdf_read
#endif /* NETCDF */

    implicit none

    ! External
    intrinsic :: achar, modulo, minval, maxval, trim, any, real

    ! Constant Parameters
    logical, parameter :: l_stdout = .false.

    ! Input Variables

    ! The interface declaration in the nr module prevents us from making
    ! this a fixed length array declaration
    real, dimension(:), intent(in) ::  & 
      param_vals_vector_r4 ! Tuning vector(ndim dimension) contains
                        ! parameterization constants (C1, C2 etc.)

    ! Local Variables

    real( kind = core_rknd ), dimension(size(param_vals_vector_r4)) :: &
      param_vals_vector ! core_rknd version of param_vals_vector_r4

    real( kind = core_rknd ), dimension(nparams) :: & 
      params_local ! Local copy of the CLUBB parameters fed into run_clubb

    ! These are read after each run from the GrADS control files
    integer ::  & 
      clubb_nz   ! Extent of the CLUBB domain in the z dimension

    character(50) ::  & 
      errorfile ! nml filename for invalid runs

    integer ::  & 
      AllocateStatus ! For hoc_zl, les_zl

    real( kind = core_rknd ) ::  & 
      err_sum  ! scalar sum of all z-levels

    integer ::  & 
      c_terms ! num of terms in err_sum (for normalization)

    ! LES and CLUBB values over nz z-levels
    real( kind = core_rknd ), dimension(:), allocatable ::  & 
      clubb_zl, & 
      clubb2_zl, & 
      les_zl, &
      clubb_grid_heights

    real( kind = core_rknd ), dimension(:,:), allocatable ::  & 
      err_sums  ! To save breakdown of cost function

    real( kind = core_rknd ) ::  & 
      les_minmax ! The largest LES value subtracted from the
    ! smallest value of all zlvl's (for normalization)

    logical ::  & 
      l_error, & ! Used to det. if reading the variable failed
      l_file_error ! Determines if reading the .netCDF file fails when checking that
                   ! The models start at the same time

    integer, dimension(c_total) ::  & 
      run_stat ! isValid over each model case

    type (stat_file) ::  &
      les_netcdf_file,   & ! Data file derived type
      clubb_netcdf_file

    integer :: &
#ifdef NETCDF
      len_file, & ! The length of the currently read in Grads or NetCDF file.
                  ! Used in a NetCDF assertion check below.
#endif /* NETCDF */
      i, j, c_run ! looping variables

    !-----------------------------------------------------------------------

    param_vals_vector = real(param_vals_vector_r4, kind = core_rknd)

    l_error = .false.

    ! Output information every 10 iterations if stdout is enabled;
    ! Amoeba's unusual calling convention makes this happen less
    ! often than might be expected.
    if ( l_results_stdout & 
        .and. ( modulo( iter, 10 ) == 0 ) & 
        .and. iter /= 0  ) then

      write(unit=*,fmt='(A12,I10)') "Iteration: ", iter
      write(unit=*,fmt='(A12)') "Parameters: "

      do i = 1, ndim, 1
        j = params_index(i)
        write(unit=*,fmt='(A30,F27.20)') params_list(j)//" = ", & 
          param_vals_vector(i)
      end do

    end if

    allocate( err_sums(c_total, v_total) )

    ! Initialize
    err_sum  = 0.0_core_rknd
    err_sums = 0.0_core_rknd
    c_terms  = 0
    err_code = clubb_no_error
    run_stat(1:c_total) = clubb_no_error

    ! Copy simplex into a vector of all possible CLUBB parameters
    do i=1, nparams, 1
      ! If the variable isn't in the simplex, leave it as is
      params_local(i) = params(i)
      do j=1, ndim, 1
        if ( i == params_index(j) ) then
          ! Copy variable from param_vals_vector argument
          params_local(i) = param_vals_vector(j)
          exit
        else
          ! Continue searching the list
          cycle
        end if
      end do
    end do

    ! Copy the values of parameters that are being held equal during tuning.
    if ( l_keep_params_equal ) then
      do i = 1, nparams, 1
        do j = 1, nparams, 1
          if ( param_constraints(i) == params_list(j) ) then
            params_local(i) = params_local(j)
          end if
        end do
      end do
    end if

!-----------------------------------------------------------------------

    ! Cycle through all the model cases specified for CLUBB

    ! OpenMP directives should work as expected now, assuming new
    ! model variables are declared threadprivate -dschanen 31 Jan 2007

!$omp parallel do default(none), private(c_run), &
!$omp   shared(params_local, run_file, run_stat, c_total, model_flags_array, iter)
    do c_run=1, c_total, 1

#ifndef _OPENMP 
      ! Write a message about which case we're calling if OpenMP is not enabled
      ! Save tuning results in file if specified
      if( l_save_tuning_run ) open(unit=file_unit, file=tuning_filename, &
        action='write', position='append')
      call write_text( "Calling CLUBB with case "//trim( run_file(c_run) ), &
        l_save_tuning_run, file_unit )
      if( l_save_tuning_run ) close(unit=file_unit)
      
#endif
      ! Run the CLUBB model with parameters as input

      if ( allocated( model_flags_array ) ) then
        call run_clubb &
             ( params_local, run_file(c_run), l_stdout, &
               model_flags_array(iter,:) )
      else
        call run_clubb & 
             ( params_local, run_file(c_run), l_stdout )
      end if

      run_stat(c_run) = err_code

      ! Reset error code for next iteration
      err_code = clubb_no_error
    end do ! 1..c_run
!$omp end parallel do

    !-----------------------------------------------------------------------

    ! Now check if CLUBB has blown up, i.e. if CLUBB has set a variable to NaN,
    ! or encountered a failure in the matrix solver routines

    ! If it has, it returns higher value than those previous to
    ! Amoeba (the downhill simplex)
    if ( any( run_stat(:) == clubb_fatal_error ) ) then
      write(fstderr,*) "Warning: the parameter set has caused CLUBB to crash"
      min_les_clubb_diff = real(2._core_rknd * maxval( cost_fnc_vector )  & 
                       - minval( cost_fnc_vector ))

      if ( l_stdout_on_invalid ) then
        inv_count = modulo( inv_count, 3 ) + 1 ! 1,2,3,1,2,3...
        errorfile = "error_crash_"// achar( inv_count+48 ) // ".in" ! Known magic number
        call output_nml_tuner( errorfile,  & 
                               param_vals_vector(1:ndim) )
      end if

      return ! return value was set, return after fatal error
    end if

    !----------------------------------------------------------------------- 

    do c_run=1, c_total, 1

      ! Determine how large the GrADS input is
      clubb_nz = stat_file_num_vertical_levels( hoc_v(1), hoc_stats_file(c_run) )

      ! Allocate the arrays for reading in the GrADS plot data
      allocate( clubb_zl(clubb_nz), clubb2_zl(clubb_nz),  & 
                les_zl(clubb_nz), clubb_grid_heights(clubb_nz), stat=AllocateStatus )

      if ( AllocateStatus /= 0 ) then
        error stop "Allocation of arrays in minimization function failed"
      end if

      ! Determine the height of GrADS input
      clubb_grid_heights = stat_file_vertical_levels( hoc_v(1), hoc_stats_file(c_run), clubb_nz )

      ! Start with first CLUBB & LES variables, then loop through and
      ! calculate the mean squared difference for all the variables
      do i=1, v_total, 1  
        ! Read in LES grads data for one variable, averaged
        ! over specified time intervals
        les_zl =  & 
        stat_file_average_interval &
        ( les_stats_file(c_run), clubb_nz,  & 
          timestep_intvls(c_run,:), les_v(i), clubb_grid_heights, 1, l_error )

#ifdef NETCDF
        ! Verify that the CLUBB and LES runs start at the same time and
        ! have the same timestep length

        ! First, be sure we are dealing with a netCDF file
        len_file = LEN_TRIM(les_stats_file(c_run))
        if (les_stats_file(c_run)(len_file-2: len_file) == ".nc") then
          call open_netcdf_read( les_v(i), les_stats_file(c_run), les_netcdf_file, l_file_error);
          call close_netcdf_read(les_netcdf_file);
          if (l_file_error) then
            l_error = .true. ! This may or may not be necessary
          end if
        else
          l_file_error = .true. ! This will cause the following assertion check to be skipped
        end if

        ! Call open_netcdf_read for CLUBB files too, to get date/time, etc. for next step.
        call open_netcdf_read( hoc_v(i), hoc_stats_file(c_run), clubb_netcdf_file, l_file_error)
        call close_netcdf_read(clubb_netcdf_file)

        if ( .not. l_file_error) then
          ! If the file could not be read, then it is most likely that the file is a GrADS file.
          ! This assertion check only supports NetCDF files. The case that the file could not
          ! be found is handled elsewhere in the tuner.
          if ( &
            clubb_netcdf_file%day /= les_netcdf_file%day &
            .or. clubb_netcdf_file%month /= les_netcdf_file%month &
            .or. clubb_netcdf_file%year /= les_netcdf_file%year &
            .or. abs(clubb_netcdf_file%time - les_netcdf_file%time) > &
                abs(clubb_netcdf_file%time + les_netcdf_file%time) / 2 * eps &
            .or. abs(clubb_netcdf_file%dtwrite - les_netcdf_file%dtwrite) > &
                abs(clubb_netcdf_file%dtwrite + les_netcdf_file%dtwrite) / 2 * eps ) then
              write(*,*) "Error: The CLUBB run and LES run do not start at the same time &
                  &or have different stat output intervals. Here are the currently set &
                  &start times and stat output intervals."

              write(*,*) "CLUBB start day: ", clubb_netcdf_file%day
              write(*,*) "CLUBB start month: ", clubb_netcdf_file%month
              write(*,*) "CLUBB start year: ", clubb_netcdf_file%year
              write(*,*) "CLUBB start time: ", clubb_netcdf_file%time
              write(*,*) "CLUBB dtwrite: ", clubb_netcdf_file%dtwrite

              write(*,*) "LES start day: ", les_netcdf_file%day
              write(*,*) "LES start month: ", les_netcdf_file%month
              write(*,*) "LES start year: ", les_netcdf_file%year
              write(*,*) "LES start time: ", les_netcdf_file%time
              write(*,*) "LES dtwrite: ", les_netcdf_file%dtwrite

              error stop "Please modify the models so that they start at the same time and &
                &have the same stat output interval." 
          end if
        end if
#endif /* NETCDF */

        if ( l_error ) then
          if( l_save_tuning_run ) then
            open(unit=file_unit, file=tuning_filename, &
            action='write', position='append')
            call write_text( "LES variable was invalid or the GrADS file does not exist", &
              l_save_tuning_run, file_unit )
            close(unit=file_unit)
          end if
          write(fstderr,*) "The specified LES variable "//trim( les_v(i) )//" was invalid, "// &
            "or the GrADS file did not exist."
          error stop "Missing or improperly formatted file"
        end if

        ! Verify that the domain that we're tuning CLUBB over is fully defined in
        ! the LES data.  If not, some points will be NaN
        if ( is_nan_2d( les_zl(z_i(c_run):z_f(c_run)) ) ) then
          write(*,*)
          write(fstderr,*) "The tuning domain exceeds the size of the LES data, "// &
            "or the LES data is NaN"
          write(*,*)
          write(fstderr,*) trim( les_v(i) )//" = ", les_zl(z_i(c_run):z_f(c_run))
          error stop "Fatal error"
        end if

        ! Read in CLUBB grads data for one variable, averaged
        ! over specified time intervals
        clubb_zl =  & 
        stat_file_average_interval & 
        ( hoc_stats_file(c_run), clubb_nz,  & 
          timestep_intvls(c_run,:), hoc_v(i), clubb_grid_heights, 1, l_error )

        if ( l_error ) then
          error stop "The specified CLUBB variable was invalid 1"
        end if

        ! The same variable, with npower = 2
        clubb2_zl =  & 
        stat_file_average_interval & 
        ( hoc_stats_file(c_run), clubb_nz, & 
          timestep_intvls(c_run,:), hoc_v(i), clubb_grid_heights, 2, l_error )

        if ( l_error ) then
          error stop "The specified CLUBB variable was invalid 2"
        end if

        !-----------------------------------------------------------------------

        ! Calculate the mean squared difference between the CLUBB
        ! and the LES variables

        ! In order to deal with differences in order of magnitude
        ! between the variables, the err_sum equation has been
        ! modified to normalize the values with respect to the
        ! the minimum and maximum in the LES. -Dave Schanen

        les_minmax = maxval( les_zl(z_i(c_run):z_f(c_run)) ) &
          - minval( les_zl(z_i(c_run):z_f(c_run)) )


        ! The following three lines are commented out and the five lines
        ! afterward are added, to avoid an error when the LES experiment has no
        ! cloud or cloud water.  This would give an error about the LES variable
        ! being zero and shut down the tuning process before, but with the new
        ! lines it will keep going and just use a minimum threshold instead of
        ! zero.

        if (hoc_v(i) == 'cloud_frac') then
          les_minmax = max( 0.01_core_rknd, les_minmax )
        elseif (hoc_v(i) == 'rcm' ) then
          les_minmax = max( 1.e-6_core_rknd, les_minmax )
        end if

        if ( abs(les_minmax) < eps) then
          error stop "An LES variable was 0 from z_i to z_f."
        end if

        ! Old code
!     err_sum = err_sum &
!      + mean_sqr_diff_zt( clubb_nz, clubb_zl, les_zl, les_minmax )

        ! Chris Golaz modification: mean_sqr_diff_2 was designed to try
        ! to limit time noise in tuning simulations.
        ! New code with Modification for weighting
        err_sums(c_run,i) = mean_sqr_diff_2( clubb_nz, z_i(c_run), z_f(c_run), clubb_zl,  & 
                              clubb2_zl, les_zl, les_minmax )

        c_terms = c_terms + 1

      end do ! i=1..v_total

      ! De-allocate the arrays for reading in the GrADS plot data
      deallocate( clubb_zl, clubb2_zl, les_zl, clubb_grid_heights )

    end do     ! end of do c_run=1, c_total

!----------------------------------------------------------------------

    ! Return error averaged over all cases, variables,
    ! and vertical levels
    ! Old Code
!       min_les_clubb_diff = err_sum / real( c_terms )

    !---------------------------------------------------------------
    ! Compute normalization factors to satisfy imposed weights
    ! This non-dimensionlizes each term in the cost function so 
    ! that the units in the variables(eg. rcm, cloud_frac) don't matter
    ! Invsigma2 is computed once and then reused in all later 
    ! calculations because invsigma2 is declared in the module
    !---------------------------------------------------------------
    if ( l_initialize_sigma ) then
      do c_run=1,c_total
        do i=1,v_total
          invsigma2(c_run,i)  & 
          = weight_case(c_run)*weight_var(i) / err_sums(c_run,i)
        end do
      end do
    end if

    !---------------------------------------------------------------
    ! Compute normalized error
    !---------------------------------------------------------------
    err_sums = invsigma2 * err_sums
    err_sum  = sum( err_sums )

    !---------------------------------------------------------------
    ! Save new minimums (used in write_results)
    !---------------------------------------------------------------    
    if ( err_sum < min_err_old ) then
      min_err_old = err_sum
      min_err_terms = err_sums
    end if

    !---------------------------------------------------------------
    ! Save total error and error contributions breakdown
    !---------------------------------------------------------------
    err_terms = err_sums
    min_les_clubb_diff = real(err_sum)

    deallocate( err_sums )

    ! Save tuning results in file if specified
    if( l_save_tuning_run ) then

        open(unit=file_unit, file=tuning_filename, action='write', position='append')
        write(file_unit,*) "Iter ", int(iter), "  Cost = ", min_les_clubb_diff, &
                                                 "  Params = ", param_vals_vector
        close(unit=file_unit)

    end if

    return
  end function min_les_clubb_diff

  !----------------------------------------------------------------------
  subroutine write_results( iunit )

! Description:
!   Outputs the results of a tuning run to the terminal if iunit = fstdout
!   or to a file which has already been opened and connected with
!   integer iunit.

! References:
!   None
!----------------------------------------------------------------------

    implicit none

    ! Input variables
    integer, intent(in) :: iunit ! = fstdout for printing to the terminal
                                 ! or = a number connected with a file

    ! Local variables
    integer :: i, j, k, c_run ! Loop iterators
    character(50) :: case_name, variable
    integer :: underscore_idx

    if ( tune_type == 0 ) then
      write(unit=iunit,fmt=*) "Number of iterations past initialization:",  iter
    end if ! tune_type == 0

    write(unit=iunit,fmt='(4x,A9,25x,A7,10x,A7)') & 
        "Parameter", "Initial", "Optimal"

    do i = 1, ndim, 1
      write(unit=iunit,fmt='(A31,2F17.10)')  & 
        params_list(params_index(i))//" = ",  & 
        params(params_index(i)), param_vals_matrix(1,i)
    end do

    ! Also print the initial and final values of the parameters that were held
    ! equal to each other.
    if ( l_keep_params_equal ) then
      do i = 1, nparams, 1
        do j = 1, ndim, 1
          k = params_index(j)
          if ( param_constraints(i) == params_list(k) ) then
            write(unit=iunit,fmt='(A31,2F17.10,5x,A)') &
              params_list(i)//" = ", params(i), param_vals_matrix(1,j), &
                        "( held equal to " // trim(params_list(k)) // " )"
          end if
        end do
      end do
    end if


    write(unit=iunit,fmt='(A20)') "Initial cost: "
    write(unit=iunit,fmt='(F15.6)') init_err
    write(unit=iunit,fmt='(A20)') "Optimal cost: "
    ! The $$ is here to make it easy to find with grep
    write(unit=iunit,fmt='(A3,F15.6)') "$$ ", min_err

    write(unit=iunit,fmt='(A,F6.3,A2)') "Approx. percent increase in accuracy: ",  &
      ((init_err - min_err) / init_err*100.0_core_rknd), "%" ! Known magic number

    write(unit=iunit,fmt=*)
    write(unit=iunit,fmt='(A)') "Approx. percent increase in accuracy by variable:"
    write(unit=iunit,fmt='(A9,5x,6x,A8,12x,A11)') &
        "Case name", "Variable", "Improvement"

    do c_run = 1, c_total
      underscore_idx = index( trim(run_file(c_run)) , "_" , .true. )
      read( run_file(c_run)(1:underscore_idx-1), * ) case_name
      do i = 1 , v_total
         variable = hoc_v(i)
          write(unit=iunit,fmt="(A20,A20,F6.3,A2)") adjustl(case_name), adjustl(variable) , &
            (init_err_terms(c_run,i)-min_err_terms(c_run,i))/init_err &
             * 100.0_core_rknd, "%"
      end do
    end do

    return
  end subroutine write_results
  !-----------------------------------------------------------------------

  subroutine output_nml_tuner( results_f, param_vals_vector )

! Description:
!   Output namelists to a formatted text file

! References:
!   None

! Notes:
!   You can do the same thing with
!   write(unit=<UNIT>,nml=<NAMELIST>), but at the time I wrote this
!   I was more ambitious and didn't like that fact that it came out
!   in all caps on pgf90. -dschanen 28 July 2006
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: achar, trim

    ! Parameter Constants
    character(len=1), parameter :: dbqt = '"'

    ! Input Variables
    character(len=*), intent(in) ::  & 
    results_f ! Name of the results file to write to

    real( kind = core_rknd ), intent(in), dimension(:) ::  & 
    param_vals_vector ! The current tuning parameters

    ! Local Variables
    real( kind = core_rknd ), dimension(nparams) :: params_local

    integer :: i, j, & ! loop variable
               iunit   ! file unit

    character :: i_c   ! loop variable in ASCII
    
    iunit = 20

    ! Open a new file
    open(unit=iunit, file=results_f,  & 
         action="write", access="sequential")

    ! Write variables to results file
    ! All this is based on the previous namelists in error.in,
    ! except for the constants parameters for CLUBB

    write(unit=iunit,fmt=*) "! Parameter file " // results_f
    write(unit=iunit,fmt=*) "&stats"
    write(unit=iunit,fmt=*) "f_tol = ", f_tol
    if (l_results_stdout) then
      write(unit=iunit,fmt=*) "l_results_stdout = " // ".true."
    else
      write(unit=iunit,fmt=*) "l_results_stdout = " // ".false."
    end if
    if (l_results_file) then
      write(unit=iunit,fmt=*) "l_results_file = " // ".true."
    else
      write(unit=iunit,fmt=*) "l_results_file = " // ".false."
    end if
    if (l_stdout_on_invalid) then
      write(unit=iunit,fmt=*) "l_stdout_on_invalid = " // ".true."
    else
      write(unit=iunit,fmt=*) "l_stdout_on_invalid = " // ".false."
    end if
    write(unit=iunit,fmt=*) "t_variables = "
    do i=1, v_total
      write(unit=iunit,fmt=*) dbqt, hoc_v(i), dbqt, ",",  & 
        dbqt, les_v(i), dbqt, ","
    end do
    write(unit=iunit,fmt=*) "weight_var_nl = ", weight_var
    write(unit=iunit,fmt=*) "anneal_temp = " , anneal_temp
    write(unit=iunit,fmt=*) "max_final_temp = " , max_final_temp
    write(unit=iunit,fmt=*) "anneal_iter = " , anneal_iter
    write(unit=iunit,fmt=*) "tune_type   = " , tune_type
    write(unit=iunit,fmt=*) "/"

    write(unit=iunit,fmt=*) "&cases"
    do i = 1, c_total
      i_c = achar( i + 48 )
      write(unit=iunit,fmt=*) "les_stats_file_nl("// i_c //") = ",  & 
        dbqt, trim( les_stats_file(i) ), dbqt
      write(unit=iunit,fmt=*) "hoc_stats_file_nl("// i_c //") = ",  & 
        dbqt, trim( hoc_stats_file(i) ), dbqt
      write(unit=iunit,fmt=*) "run_file_nl("// i_c //") = ",  & 
        dbqt, trim( run_file(i) ), dbqt
      write(unit=iunit,fmt=*) "z_i_nl("// i_c //")  = " , z_i(i)
      write(unit=iunit,fmt=*) "z_f_nl("// i_c //")  = " ,  z_f(i)
      write(unit=iunit,fmt=*) "weight_case_nl("// i_c //") = " ,  & 
        weight_case(i)
      write(unit=iunit,fmt=*) "timestep_intvls_nl("// i_c //",:)  = " , timestep_intvls(i,:)
    end do
    write(unit=iunit,fmt=*) "/"

    write(unit=iunit,fmt=*) "&clubb_params_nl"

    ! Copy simplex into a vector of all possible CLUBB parameters
    do i=1, nparams, 1
      ! If the variable isn't in the simplex, leave it as is
      params_local(i) = params(i)
      do j=1, ndim, 1
        if ( i == params_index(j) ) then
          ! Copy variable from param_vals_vector argument
          params_local(i) = param_vals_vector(j)
          exit
        else
          cycle
        end if
      end do
    end do

    ! Output optimal values and all possible CLUBB parameters
    do i=1, nparams, 1
      write(unit=iunit,fmt='(A30,F27.20)') & 
        trim( params_list(i) )//" = ", params_local(i)
    end do
    write(unit=iunit,fmt=*) "/"

    write(unit=iunit,fmt=*) "&initmax"
    ! Copy the max values into a vector of all possible CLUBB parameters
    do i=1, nparams, 1
      ! If the variable isn't being changed, set it to zero
      params_local(i) = 0.0_core_rknd
      do j=1, ndim, 1
        if ( i == params_index(j) ) then
          ! Copy variable from param_vals_vector argument
          params_local(i) = param_vals_minmax(2,j)
          exit
        else
          cycle
        end if
      end do
    end do

    ! Output the amount each variable was changed for the simplex
    do i=1, nparams, 1
      write(unit=iunit,fmt='(a30,f12.5)') & 
        trim( params_list(i) )//" = ", params_local(i)
    end do
    write(unit=iunit,fmt=*) "/"

    ! Close new namelist file
    close(unit=iunit)

    return
  end subroutine output_nml_tuner

!-----------------------------------------------------------------------
  subroutine output_nml_standalone ( results_f,  & 
                                     param_vals_vector )

! Description:
!   Output namelists to a formatted file

! References:
!   None

! Notes:
!   See note for output_nml_tuner, above. dschanen 28 July 2006
!-----------------------------------------------------------------------

    implicit none

    ! External
    intrinsic :: achar, trim

    ! Input variables
    character(len=*), intent(in) :: & 
      results_f ! Results file to write to

    real( kind = core_rknd ), intent(in), dimension(ndim) :: & 
      param_vals_vector ! the current constants

    ! Local variables
    real( kind = core_rknd ), dimension(nparams) :: params_local

    integer   :: i, j  ! loop variables

    integer :: iunit ! file unit

    ! ---- Begin Code ---
    write(6,*) "Writing parameters for namelist to "//results_f

    iunit = 20
    ! Open new file
    open(unit=iunit, file=results_f,  & 
         action="write", access="sequential")

    ! Write variables to namelist for standalone CLUBB.
    ! All this is based on the previous error.in, except the constants

    write(unit=iunit,fmt=*) "! Parameter file " // results_f

    write(unit=iunit,fmt=*) "&clubb_params_nl"
    ! Copy simplex into a vector of all possible CLUBB parameters
    do i=1, nparams, 1
      ! If the variable isn't in the simplex, leave it as is
      params_local(i) = params(i)
      do j=1, ndim, 1
        if ( i == params_index(j) ) then
          ! Copy variable from param_vals_vector argument
          params_local(i) = param_vals_vector(j)
          exit
        else
          cycle
        end if
      end do
    end do

    ! Copy the values of parameters that were held equal during tuning.
    if ( l_keep_params_equal ) then
      do i = 1, nparams, 1
        do j = 1, nparams, 1
          if ( param_constraints(i) == params_list(j) ) then
            params_local(i) = params_local(j)
          end if
        end do
      end do
    end if

    ! Output optimal values and all possible CLUBB parameters
    do i=1, nparams, 1
      write(unit=iunit,fmt='(A30,F27.20)') & 
        trim( params_list(i) )//" = ", params_local(i)
    end do
    write(unit=iunit,fmt=*) "/"

    ! Close new file
    close(unit=iunit)

    return
  end subroutine output_nml_standalone

!-----------------------------------------------------------------------
  real( kind = core_rknd ) function mean_sqr_diff & 
                ( nz, z_init, z_final, scm_zl, les_zl, norm_term )
! Description
!   Calculate the mean squared difference between two input vectors,
!   then normalize.

! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: sum

    ! Input Variables
    integer, intent(in) ::  & 
      nz,      & ! Vertical extent for CLUBB
      z_init,  & ! Initial point for the purposes of computing the sum
      z_final    ! Final point for the purpose of computign the sum

    real( kind = core_rknd ), intent(in), dimension(nz) ::  & 
      scm_zl, &! CLUBB GrADS variable   [units vary]
      les_zl   ! The LES GrADS variable [units vary]

    real( kind = core_rknd ), intent(in) ::  & 
      norm_term ! normalization term; typically maxval(les) - minval(les)

    real( kind = core_rknd ), dimension(nz) ::  & 
      tmp_zl

    integer :: k

!----------------------------------------------------------------------

    ! ---- Begin Code ----

    tmp_zl = 0.0_core_rknd

    do k = z_init, z_final, 1
      tmp_zl(k) = ( ( scm_zl(k) - les_zl(k) ) / norm_term )**2
    end do

    mean_sqr_diff = sum( tmp_zl )

    return
  end function mean_sqr_diff
!-----------------------------------------------------------------------
  real( kind = core_rknd ) function mean_sqr_diff_2 & 
                ( nz, z_init, z_final, scm_zl,  & 
                  scm2_zl, les_zl, norm_term )
! Description:
!   Alternate function to compute mean difference between input
!   fields.
!   It computes:
!     scm2_zl - 2 * scm_zl * les_zl + les_zl**2
!     where scm2_zl = avg( scm_zl**2 )
!      scm_zl  = avg( scm_zl )
!      les_zl  = avg( les_zl )
!   This alternate formulation adds a penalty to the cost function
!   from the time varying noise that might be present in a simulation.
!   It allows the tuner to avoid very noisy simulations, although some
!   noise might still be present.
!
!   Configured to do interpolation on LES / CLUBB comparisons on the
!   CLUBB grid

! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: sum

    ! Input Variables
    integer, intent(in) ::  & 
      nz,      & ! Vertical extent for CLUBB
      z_init,  & ! Initial point for the purposes of computing the sum
      z_final    ! Final point for the purpose of computign the sum

    real( kind = core_rknd ), intent(in), dimension(nz) ::  & 
      scm_zl,  & ! CLUBB GrADS variable [units vary]
      scm2_zl, & ! CLUBB GrADS variable [units vary]
      les_zl     ! The LES GrADS variable [units vary]

    real( kind = core_rknd ), intent(in) ::  & 
      norm_term ! normalization term. Typically maxval(les) - minval(les)

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) ::  & 
      tmp_zl

    integer :: k

!----------------------------------------------------------------------

    ! ---- Begin Code ----

    tmp_zl = 0.0_core_rknd

    do k = z_init, z_final, 1
      tmp_zl(k) = ( scm2_zl(k) - 2.0_core_rknd * scm_zl(k)*les_zl(k) + les_zl(k)**2 ) &
        / norm_term**2
    end do

    mean_sqr_diff_2 = sum( tmp_zl )

    return
  end function mean_sqr_diff_2

!-----------------------------------------------------------------------

  subroutine read_random_seed( seed_file )

! Description:
!   Reads an ASCII flat file passed as an argument

! References:
!   None
!-----------------------------------------------------------------------
    use constants_clubb, only: fstderr ! Constant

    use mt95, only: genrand_init ! Procedure

    use mt95, only: genrand_intg ! Constant

    implicit none

    ! Parameter constants
    integer, parameter :: rand_size = 34

    ! Input
    character(len=*), intent(in) ::  & 
      seed_file ! This should usually contain >= 34 integers

    ! Local Variables
    integer(kind=genrand_intg), dimension(rand_size) ::  & 
      rand_seed  ! Set of 32 bit integers for seeding the generator

    integer :: & 
      InputStatus, iunit

!-----------------------------------------------------------------------

    iunit = 30

    ! ASCII formatted file, usually generated by int2txt
    open(unit=iunit, file=seed_file, action='read')

    read(unit=iunit, fmt=*, iostat=InputStatus) rand_seed(1:rand_size)
    if ( InputStatus /= 0 ) then
      write(fstderr,*) "Error reading "//seed_file
      error stop
    end if

    close(unit=iunit)

    call genrand_init( put=rand_seed )

    return
  end subroutine read_random_seed

!-----------------------------------------------------------------------

end module error
!-----------------------------------------------------------------------
