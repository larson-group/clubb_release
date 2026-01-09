!----------------------------------------------------------------------
! $Id$

module clubb_driver

! Description:
!   Contains the necessary subroutines to execute individual CLUBB
!   model runs, using one of the driver programs (the simplest case
!   being the clubb_standalone program).
!-----------------------------------------------------------------------

  !--------------------------------- Procedures ---------------------------------

  use text_writer, only: &
    write_text, &
    write_date

#ifdef _OPENMP
  ! Because Fortran I/O is not thread safe, we use this here to
  ! ensure that no model uses the same file number simultaneously
  ! when doing a tuning run. -dschanen 31 Jan 2007
  use omp_lib, only: &
    omp_get_thread_num ! Function
#endif

  !--------------------------------- Types ---------------------------------

  use grid_class, only: &
    grid

  use stats_type, only: &
    stats

  use parameters_tunable, only: &
    nu_vertical_res_dep

  use model_flags, only: &
    clubb_config_flags_type

  use corr_varnce_module, only: &
    hm_metadata_type 
    
  use hydromet_pdf_parameter_module, only: &
    hydromet_pdf_parameter

  use parameters_silhs, only: &
    silhs_config_flags_type

  use array_index, only: &
    sclr_idx_type

  use err_info_type_module, only: &
    err_info_type 

  use clubb_api_module, only: &
      precipitation_fractions

  use pdf_parameter_module, only: &
      pdf_parameter, &
      implicit_coefs_terms

  use stats_variables, only: &
      stats_metadata_type

  !--------------------------------- Constants ---------------------------------

  use clubb_precision, only: &
    time_precision, &
    core_rknd

  use constants_clubb, only: &
    fstdout, & 
    fstderr, &
    three_halves, &
    one, &
    one_half, &
    zero, &
    two, &
    radians_per_deg, &
    rt_tol, &
    thl_tol, &
    w_tol, &
    w_tol_sqd, &
    em_min, &
    eps, &
    Lv, &
    kappa, &
    Cp, &
    p0, &
    omega_planet

  use error_code, only: &
    clubb_fatal_error, & 
    clubb_generalized_grd_test_err, &
    clubb_no_error

  use parameters_tunable, only: &
    params_list

  use model_flags, only: &
    l_silhs_rad, &
    l_pos_def, &
    l_gamma_Skw, &
    l_byteswap_io, &
    saturation_flatau, &
    l_test_grid_generalization, &
    no_grid_adaptation

  use parameter_indices, only: &
    nparams, &
    ic_K, &
    iSkw_max_mag, &
    ibv_efold

  use sounding, only: &
    sclr_max

  !--------------------------------- Setting Variables ---------------------------------
  ! These are variables meant to be set once at the beginning and should be left unchanged.
  ! Module variables are frowned upon, but these are the more acceptable type.
  ! Each of these should be set within the call to "init_clubb_case"
  
  use inputfields, only: &
    l_input_wp3, &
    stat_files

  use sponge_layer_damping, only: &
    thlm_sponge_damp_settings,    &
    rtm_sponge_damp_settings,     &
    uv_sponge_damp_settings,      &
    wp2_sponge_damp_settings,     &
    wp3_sponge_damp_settings,     &
    up2_vp2_sponge_damp_settings

  use time_dependent_input, only: &
    l_t_dependent,     &
    l_input_xpwp_sfc,  &
    l_ignore_forcings, &
    l_sfc_already_initialized

  use parameters_radiation, only: &
    rad_scheme

  use soil_vegetation, only: &
    l_soil_veg 

  use parameters_microphys, only: &
    microphys_scheme, & 
    sigma_g,              & ! Parameter used in the cloud droplet sedimentation code
    l_cloud_sed             ! Cloud water sedimentation (K&K or no microphysics)

  use extended_atmosphere_module, only: &
    total_atmos_dim, &
    complete_alt, &
    complete_momentum

#ifdef SILHS
  use parameters_microphys, only: &
    lh_microphys_type,     &
    lh_microphys_disabled, &
    microphys_scheme,      &
    lh_num_samples
#endif

  implicit none

  public :: &
    run_clubb, &
    init_clubb_case, &
    set_case_initial_conditions, &
    advance_clubb_to_end, &
    clean_up_clubb

  private ! Default to private

  !--------------------------------- Constants ---------------------------------

  logical, parameter :: &
    l_implemented     = .false.,  & ! CLUBB is not being run inside of a host model
    l_ascending_grid  = .true.,   & ! CLUBB runs with an ascending grid by default
    l_write_to_file   = .true.      ! If true, will write case information to a file

  integer, parameter :: &
    nlon = 1, & ! Number of points in the X/Y [-]
    nlat = 1 

  character(len=100), parameter :: &
    output_dir = "../output/"   ! Output directory

  real( kind = core_rknd ) , parameter ::  &
    timing_tol = 0.01_core_rknd   ! allowed tolerance for the timing budget check

  
  !--------------------------------- Setup Variables ---------------------------------
  ! Things that should be set in "init_clubb_case" and never updated

  ! 'public' only to fix compiler error for nvhpc versions <24.9, please remove in future
  logical, public :: &
    l_restart,                  & ! Flag for restarting from GrADS file
    l_input_fields,             & ! Whether to set model variables from a file
    l_modify_ic_with_cubic_int, & ! Flag for interpolating the sounding profile with Steffen's monotone cubic 
                                  ! method to obtain smoother initial condition profile, which is found to be 
                                  ! beneficial to achive a better numerical solution convergence. If this flag 
                                  ! is turned off, the initial conditions will be generated with linear interpolation.
                                  ! This is done on a case-by-case basis, since using the monotone cubic method
                                  ! requires a special sounding.in file with many additional sounding levels.
    l_modify_bc_for_cnvg_test, &  ! Flag to activate modifications on boundary condition for convergence test
                                  ! (surface fluxes computed at fixed 25 m height).
    l_stats,                    & ! Whether statistics are computed and output to disk
    l_silhs_out,                & ! Whether to output SILHS files
    l_rad_itime,                & ! .true. if we calculate radiation at the current time step
    l_first_write_to_file,      & ! whether this is the first time data gets written
                                  ! to the netcdf file or not
    l_allow_small_stats_tout,   &
    l_restart_input,            &
    l_last_timestep

  ! 'public' only to fix compiler error for nvhpc versions <24.9, please remove in future
  integer, public :: &
    ifinal,                 & ! Total number of iterations (final iter value)
    stats_nsamp,            & ! Stats sampling interval [timestep]
    stats_nout,             & ! Stats output interval   [timestep]
    iunit,                  & ! File unit used for I/O
    iunit_grid_adaptation,  & ! Different file unit for grid_adaptation file, which stores
                              ! how grid adapts
    hydromet_dim,           & ! Number of hydrometeor species        [#]
    sclr_dim,               & ! Number of passive scalars            [#]
    edsclr_dim,             & ! Number of passive scalars            [#]
    pdf_dim

  ! 'public' only to fix compiler error for nvhpc versions <24.9, please remove in future
  integer, public ::  & 
    nzmax, &      ! Vertical extent in levels( relevant for grid type 2 and 3 only )  [#]
    grid_type     ! 1 ==> evenly-spaced grid levels
                  ! 2 ==> stretched (unevenly-spaced) grid entered on
                  !       thermodynamic grid levels; momentum levels
                  !       halfway between thermodynamic levels (style
                  !       of SAM stretched grid).
                  ! 3 ==> stretched (unevenly-spaced) grid entered on
                  !       momentum grid levels; thermodynamic levels
                  !       halfway between momentum levels (style
                  !       of WRF stretched grid).

  ! Radiation variables
  integer :: &
    extended_atmos_bottom_level, & ! Bottom level of the extended atmosphere
    extended_atmos_top_level,    & ! Top level of the extended atmosphere
    extended_atmos_range_size      ! The number of levels in the extended atmosphere

  ! 'public' only to fix compiler error for nvhpc versions <24.9, please remove in future
  integer, public ::  &
    day, month, year ! Start time the of simulation

  ! The number of interpolated levels between the computational grid
  ! and the extended atmosphere
  integer :: &
    lin_int_buffer

  ! 'public' only to fix compiler error for nvhpc versions <24.9, please remove in future
  integer, public :: &
    sfctype ! 0: fixed sfc sensible and latent heat fluxes as
            !    given in namelist
            ! 1: bulk formula: uses given surface temperature
            !    and assumes over ocean

  type (stats_metadata_type) :: &
    stats_metadata

  type (sclr_idx_type) :: &
    sclr_idx

  type (hm_metadata_type) :: &
    hm_metadata

  type(clubb_config_flags_type) :: &
    clubb_config_flags ! Derived type holding all configurable CLUBB flags

  type(silhs_config_flags_type) :: &
    silhs_config_flags ! Flags for the SILHS sampling code

  ! 'public' only to fix compiler error for nvhpc versions <24.9, please remove in future
  real( kind = core_rknd ), public :: &
    stats_tsamp,   & ! Stats sampling interval [s]
    stats_tout       ! Stats output interval   [s]

  ! 'public' only to fix compiler error for nvhpc versions <24.9, please remove in future
  real( kind = time_precision ), public :: &
    time_restart    ! Time of model restart run     [s]
  
  ! 'public' only to fix compiler error for nvhpc versions <24.9, please remove in future
  character(len=128), public :: &
    restart_path_case,  & ! GrADS file used in case of restart
    forcings_file_path    ! Path to the forcing files

  ! 'public' only to fix compiler error for nvhpc versions <24.9, please remove in future
  character(len=10), public :: &
    stats_fmt             ! File format for stats; typically GrADS.

  ! 'public' only to fix compiler error for nvhpc versions <24.9, please remove in future
  character(len=100), public :: &
    fname_prefix           ! Prefix of stats filenames, to be followed by, for example "_zt"

  character(len=150) :: &
    case_info_file, &     ! The filename for case info
    fname_grid_adaptation ! The filename for the grid adaptation file which stores
                          ! how the grid is adapted

  character(len=100) :: &
    output_file_prefix

  ! 'public' only to fix compiler error for nvhpc versions <24.9, please remove in future
  character(len=50), public ::  & 
    runtype ! String identifying the model case; e.g. bomex

  ! 'public' only to fix compiler error for nvhpc versions <24.9, please remove in future
  ! For grid_type 2 or 3 (stretched grid cases)
  character(len=100), public :: & 
    zt_grid_fname,  & ! Path and filename of thermodynamic level altitudes
    zm_grid_fname     ! Path and filename of momentum level altitudes

  ! 'public' only to fix compiler error for nvhpc versions <24.9, please remove in future
  real( kind = core_rknd ), public ::  & 
    lat_vals, & ! Latitude  [Degrees North]
    lon_vals    ! Longitude [Degrees East]

  ! 'public' only to fix compiler error for nvhpc versions <24.9, please remove in future
  real(kind = time_precision ), public :: & 
    time_initial, & ! Time of start of simulation     [s]
    time_final,   & ! Time end of simulation          [s]
    time_current   !  Current time of simulation      [s]

  ! 'public' only to fix compiler error for nvhpc versions <24.9, please remove in future
  real(kind = core_rknd ), public ::  & 
    dt_main,  & ! Main model timestep                    [s]
    dt_rad      ! Closure model timestep                 [s]

  ! 'public' only to fix compiler error for nvhpc versions <24.9, please remove in future
  real( kind = core_rknd ), public ::  & 
    deltaz_nl,  & ! Change in altitude per grid level     [m]
    zm_init_nl, & ! Initial point on the momentum grid    [m]
    zm_top_nl     ! Maximum point on the momentum grid    [m]

  ! 'public' only to fix compiler error for nvhpc versions <24.9, please remove in future
  real( kind = core_rknd ), public :: &
    mixt_frac_max_mag, &    ! Maximum allowable mag. of mixt_frac   [-]
    T0, &                   ! Reference temperature (usually 300)  [K]
    ts_nudge, &             ! Timescale of u/v nudging             [s]
    rtm_min, &              ! Value below which rtm will be nudged [kg/kg]
    rtm_nudge_max_altitude  ! Highest altitude at which to nudge rtm [m]

  real( kind = core_rknd ), dimension(1) ::&
    rad_dummy ! Dummy variable for radiation levels

  real( kind = core_rknd ), dimension(1) :: &
    deep_soil_T_in_K_init, &
    sfc_soil_T_in_K_init, &
    veg_T_in_K_init
    
  real( kind = core_rknd ), dimension(:), allocatable :: &
    dummy_dx, &
    dummy_dy  ! [m]
    
  real( kind = core_rknd ), dimension(:), allocatable :: &
    sclr_tol        ! Thresholds on the passive scalars     [units vary]

  ! Variables used to initialize multi_col version of variables
  real( kind = core_rknd ), dimension(:), allocatable :: &
    thlm_init, &
    rtm_init, &
    um_init, &
    vm_init, &
    ug_init, &
    vg_init, &
    wp2_init, &
    up2_init, &
    vp2_init, &
    upwp_init, &
    rcm_init, &
    wm_zt_init, &
    wm_zm_init, &
    em_init, &
    exner_init, &
    thvm_init, &
    p_in_Pa_init, &
    rho_init, &
    rho_zm_init, &
    rho_ds_zm_init, &
    rho_ds_zt_init, &
    invrs_rho_ds_zm_init, &
    invrs_rho_ds_zt_init, &
    thv_ds_zm_init, &
    thv_ds_zt_init, &
    rtm_ref_init, &
    thlm_ref_init, &
    um_ref_init, &
    vm_ref_init, &
    Ncm_init, &
    Nc_in_cloud_init, &
    Nccnm_init, &
    rho_ds_zm_dycore_init
    
  real( kind = core_rknd ), dimension(:), allocatable :: &
    p_sfc,          & ! surface pressure              [Pa]
    T_sfc,          &
    fcor,           & ! Traditional Coriolis parameter    [s^-1]
                      ! Vertical planetary vorticity.   Proportional to sin(latitude) 
    fcor_y,         & ! Nontraditional Coriolis parameter [s^-1]
                      ! Meridional planetary vorticity. Proportional to cos(latitude)
    sfc_elevation,  &
    zm_init,        &
    zm_top,         &
    deltaz

  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    sclrm_init, &
    edsclrm_init
    
  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    corr_array_n_cloud, &
    corr_array_n_below

  !----------------------------------- Dynamic Variables -----------------------------------
  ! These variables are used to store field values that can update each timestep.
  ! They should be allocated in "init_clubb_case", "set_case_initial_conditions"
  ! should reset them, and "advance_clubb_to_end" can modify them.

  integer :: &
    iinit
  
  integer, dimension(:,:,:), allocatable :: &
    X_mixt_comp_all_levs ! Which mixture component a sample is in

  type(grid) :: &
    gr, &
    gr_desc, &
    gr_dycore, & ! only set and used, if l_add_dycore_grid is .true.
    gr_output, &
    gr_fixed_min ! the read in grid for the minimum grid density profile
                 ! for the grid adaptation; only set and used,
                 ! if grid_adapt_in_time_method > no_grid_adaptation
    
  type(pdf_parameter) :: &
    pdf_params,     & ! PDF parameters (thermodynamic levels)    [units vary]
    pdf_params_zm     ! PDF parameters on momentum levels        [units vary]

  type(implicit_coefs_terms) :: &
    pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

  type(precipitation_fractions) :: &
    precip_fracs           ! Precipitation fractions      [-]

  type(nu_vertical_res_dep) :: &
    nu_vert_res_dep    ! Vertical resolution dependent nu values

  type (stats), dimension(:), allocatable :: &
    stats_zt,      & ! stats_zt grid
    stats_zm,      & ! stats_zm grid
    stats_lh_zt,   & ! stats_lh_zt grid
    stats_lh_sfc,  & ! stats_lh_sfc grid
    stats_rad_zt,  & ! stats_rad_zt grid
    stats_rad_zm,  & ! stats_rad_zm grid
    stats_sfc        ! stats_sfc

  type(hydromet_pdf_parameter), dimension(:,:), allocatable :: &
    hydromet_pdf_params    ! Hydrometeor PDF parameters      [units vary]

  ! 'public' only to fix compiler error for nvhpc versions <24.9, please remove in future
  real( kind = core_rknd ), public :: &
    sens_ht,          & ! sensible heat flux                          [K m/s]
    latent_ht,        & ! latent heat flux                            [m/s]
    lmin,             & ! Min. value for the length scale             [m]
    vert_decorr_coef    ! Empirically defined de-correlation constant [-]

  real( kind = core_rknd ), dimension(:), allocatable :: &
    deep_soil_T_in_K, &
    sfc_soil_T_in_K, &
    veg_T_in_K

  real( kind = core_rknd ), dimension(:), allocatable :: &
    wpthlp_sfc,     & ! w' theta_l' at surface   [(m K)/s]
    wprtp_sfc,      & ! w' r_t' at surface       [(kg m)/( kg s)]
    upwp_sfc,       & ! u'w' at surface          [m^2/s^2]
    vpwp_sfc,       & ! v'w' at surface          [m^2/s^2]
    upwp_sfc_pert,  & ! pertubed u'w' at surface    [m^2/s^2]
    vpwp_sfc_pert     ! pertubed v'w' at surface    [m^2/s^2]

  real( kind = core_rknd ), dimension(:), allocatable :: &
    gr_dens_z, &     ! levels at which the grid density is given      [m]
    gr_dens          ! grid density at the given levels in gr_dens_z  [1/m]

  ! Save the normalized grid density from the last grid adaptation, to check if new grid would
  ! be significantly different (grid adaptation trigger)
  real( kind = core_rknd ), dimension(:,:), allocatable ::  &
    ! The two arrays gr_dens_old_z and gr_dens_old build the piecewise linear
    ! normalized grid density, from the last time the grid was adapted
    gr_dens_old_z,  & ! altitudes [m]
    gr_dens_old       ! densities [1/m]

  ! Save the calculated values for the refinement criterion terms as a cumulative sum starting
  ! from the last grid adaptation, to take the time average over that region, to have a
  ! smoother evolution in time (x_counter + cumulative_x)
  integer :: &
    Lscale_counter, &          ! the number of values cumulated in cumulative_Lscale
    richardson_num_counter, &  ! the number of values cumulated in cumulative_richardson_num
    chi_counter                ! the number of values cumulated in cumulative_chi

  real( kind = core_rknd ), dimension(:), allocatable :: &
    cumulative_Lscale, &          ! cumulative sum of the Lscale term in the
                                  ! refinement cirterion                                  [1/m]
    cumulative_richardson_num, &  ! cumulative sum of the richardson number term in the
                                  ! refinement cirterion                                  [1/m]
    cumulative_chi                ! cumulative sum of the chi term in the
                                  ! refinement cirterion                                  [1/m]
    
  real( kind = core_rknd ), dimension(:), allocatable :: &
    alt_term, &
    lscale_term, &
    lscale_term_time_avg, &
    chi_term, &
    chi_term_time_avg, &
    richardson_num_term, &
    richardson_num_term_time_avg, &
    norm_min_grid_dens, &
    norm_grid_dens

  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    um,      & ! eastward grid-mean wind component (thermo. levs.)  [m/s]
    upwp,    & ! u'w' (momentum levels)                             [m^2/s^2]
    vm,      & ! northward grid-mean wind component (thermo. levs.) [m/s]
    vpwp,    & ! v'w' (momentum levels)                             [m^2/s^2]
    up2,     & ! u'^2 (momentum levels)                             [m^2/s^2]
    vp2,     & ! v'^2 (momentum levels)                             [m^2/s^2]
    up3,     & ! u'^3 (thermodynamic levels)                        [m^3/s^3]
    vp3,     & ! v'^3 (thermodynamic levels)                        [m^3/s^3]
    wprtp,   & ! w' r_t' (momentum levels)                          [(kg/kg) m/s]
    wpthlp,  & ! w'th_l' (momentum levels)                          [(m/s) K]
    rtp2,    & ! r_t'^2 (momentum levels)                           [(kg/kg)^2]
    rtp3,    & ! r_t'^3 (thermodynamic levels)                      [(kg/kg)^3]
    thlp2,   & ! th_l'^2 (momentum levels)                          [K^2]
    thlp3,   & ! th_l'^3 (thermodynamic levels)                     [K^3]
    rtpthlp, & ! r_t'th_l' (momentum levels)                        [(kg/kg) K]
    wp2,     & ! w'^2 (momentum levels)                             [m^2/s^2]
    wp2_zt     ! w'^2 on thermo. grid                                  [m^2/s^2]
  
  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    thvm,     & ! Virtual potential temperature                  [K]
    exner,    & ! Exner function (thermodynamic levels)          [-]
    rtm,      & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
    thlm,     & ! liq. water pot. temp., th_l (thermo. levels)   [K]
    rcm,      & ! cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
    wp3,      & ! w'^3 (thermodynamic levels)                    [m^3/s^3]
    wp3_zm,   & ! w'^3 (momentum levels)                         [m^3/s^3]
    delta_zm
  
  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    p_in_Pa,    & ! Air pressure (thermodynamic levels)         [Pa]
    cloud_frac, & ! cloud fraction (thermodynamic levels)       [-]
    wpthvp,     & ! < w' th_v' > (momentum levels)              [kg/kg K]
    wp2thvp,    & ! < w'^2 th_v' > (thermodynamic levels)       [m^2/s^2 K]
    wp2up,      & ! < w'^2 u' > (thermodynamic levels)          [m^3/s^3]
    rtpthvp,    & ! < r_t' th_v' > (momentum levels)            [kg/kg K]
    thlpthvp,   & ! < th_l' th_v' > (momentum levels)           [K^2]
    wp2rtp,     & ! w'^2 rt' (thermodynamic levels)             [m^2/s^2 kg/kg]
    wp2thlp,    & ! w'^2 thl' (thermodynamic levels)            [m^2/s^2 K]
    uprcp,      & ! < u' r_c' > (momentum levels)               [(m/s)(kg/kg)]
    vprcp,      & ! < v' r_c' > (momentum levels)               [(m/s)(kg/kg)]
    rc_coef_zm, & ! Coefficient of X'r_c' in Eq. (34) (m-levs.) [K/(kg/kg)]
    wp4,        & ! w'^4 (momentum levels)                      [m^4/s^4]
    wpup2,      & ! w'u'^2 (thermodynamic levels)               [m^3/s^3]
    wpvp2,      & ! w'v'^2 (thermodynamic levels)               [m^3/s^3]
    wp2up2,     & ! w'^2 u'^2 (momentum levels)                 [m^4/s^4]
    wp2vp2        ! w'^2 v'^2 (momentum levels)                 [m^4/s^4]

  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    rho_ds_zt,  & ! Dry, static density on thermo. levels      [kg/m^3]
    thv_ds_zt     ! Dry, base-state theta_v on thermo levs.    [K]
  
  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    wm_zm,           & ! vertical mean wind comp. on momentum levs  [m/s]
    wm_zt,           & ! vertical mean wind comp. on thermo. levs   [m/s]
    rho,             & ! Air density on thermodynamic levels        [kg/m^3]
    rho_zm,          & ! Air density on momentum levels             [kg/m^3]
    rho_ds_zm,       & ! Dry, static density on momentum levels     [kg/m^3]
    invrs_rho_ds_zm, & ! Inverse dry, static density on m-levs.     [m^3/kg]
    invrs_rho_ds_zt, & ! Inverse dry, static density on thermo levs.[m^3/kg]
    thv_ds_zm          ! Dry, base-state theta_v on momentum levs.  [K]

  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    thlm_forcing,    & ! liq. wat. pot. temp. forcing (thermo. levs)[K/s]
    rtm_forcing,     & ! total water forcing (thermo. levels)       [(kg/kg)/s]
    um_forcing,      & ! eastward wind forcing (thermo. levels)     [m/s/s]
    vm_forcing,      & ! northward wind forcing (thermo. levels)    [m/s/s]
    wprtp_forcing,   & ! total water turbulent flux forcing (m-levs)[m*K/s^2]
    wpthlp_forcing,  & ! liq pot temp turb flux forcing (m-levs)    [m(kg/kg)/s^2]
    rtp2_forcing,    & ! total water variance forcing (m-levs)      [(kg/kg)^2/s]
    thlp2_forcing,   & ! liq pot temp variance forcing (m-levs)     [K^2/s]
    rtpthlp_forcing    ! <r_t'th_l'> covariance forcing (m-levs)    [K(kg/kg)/s]

  ! Variables used to track perturbed version of winds.
  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    um_pert,        & ! perturbed <u>       [m/s]
    vm_pert,        & ! perturbed <v>       [m/s]
    upwp_pert,      & ! perturbed <u'w'>    [m^2/s^2]
    vpwp_pert         ! perturbed <v'w'>    [m^2/s^2]

  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    Skw_zm,         & ! Skewness of w on momentum levels  [-]
    Skw_zm_smooth,  & ! Smoothed version of Skw_zm
    Lscale            ! Length scale                      [m]
    
  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    thlprcp,        & ! thl'rc'                                      [K kg/kg]
    sigma_sqd_w,    & ! PDF width parameter (momentum levels)        [-]
    sigma_sqd_w_zt    ! PDF width parameter interpolated to t-levs.  [-]

  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    wprcp,                 & ! w'r_c' (momentum levels)              [(kg/kg) m/s]
    w_up_in_cloud,         & ! Average cloudy updraft velocity       [m/s]
    w_down_in_cloud,       & ! Average cloudy downdraft velocity     [m/s]
    cloudy_updraft_frac,   & ! cloudy updraft fraction               [-]
    cloudy_downdraft_frac, & ! cloudy downdraft fraction             [-]
    ice_supersat_frac,     & ! ice cloud fraction (thermo. levels)   [-]
    rcm_in_layer,          & ! rcm within cloud layer                [kg/kg]
    cloud_cover,           & ! cloud cover                           [-]
    invrs_tau_zm             ! One divided by tau on zm levels       [1/s]

  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    ug,       & ! u geostrophic wind                           [m/s]
    vg,       & ! v geostrophic wind                           [m/s]
    rtm_ref,  & ! Initial total water mixing ratio             [kg/kg]
    thlm_ref, & ! Initial liquid water potential temperature   [K]
    um_ref,   & ! Initial u wind                               [m/s]
    vm_ref      ! Initial v wind                               [m/s]

  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    rcm_mc,     & ! Tendency of liquid water due to microphysics      [kg/kg/s]
    rvm_mc,     & ! Tendency of vapor water due to microphysics       [kg/kg/s]
    thlm_mc,    & ! Tendency of liquid pot. temp. due to microphysics [K/s]
    wprtp_mc,   & ! Microphysics tendency for <w'rt'>                 [m*(kg/kg)/s^2]
    wpthlp_mc,  & ! Microphysics tendency for <w'thl'>                [m*K/s^2]
    rtp2_mc,    & ! Microphysics tendency for <rt'^2>                 [(kg/kg)^2/s]
    thlp2_mc,   & ! Microphysics tendency for <thl'^2>                [K^2/s]
    rtpthlp_mc, & ! Microphysics tendency for <rt'thl'>               [K*(kg/kg)/s]
    Ncm_mc        ! Microphysics tendency for Ncm                     [num/kg/s]

  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    Ncm,          & ! Mean cloud droplet concentration, <N_c> (t-levs.)     [num/kg]
    Nccnm,        & ! Cloud condensation nuclei concentration (COAMPS/MG)   [num/kg]
    em,           & ! Turbulent Kinetic Energy (TKE)                        [m^2/s^2]
    tau_zm,       & ! Eddy dissipation time scale on momentum levels        [s]
    tau_zt,       & ! Eddy dissipation time scale on thermodynamic levels   [s]
    Kh_zt,        & ! Eddy diffusivity coefficient on thermodynamic levels  [m^2/s]
    Kh_zm,        & ! Eddy diffusivity coefficient on momentum levels       [m^2/s]
    rfrzm,        & ! Total ice-phase water mixing ratio                    [kg/kg]
    rrm,          & ! Overall mean rain water mixing ratio                  [kg/kg]
    Nrm,          & ! Overall mean rain drop concentration                  [num/kg]
    Nc_in_cloud,  & ! Mean (in-cloud) cloud droplet concentration           [num/kg]
    wpNcp           ! Covariance of w and N_c, <w'N_c'> (momentum levels)   [(m/s)(num/kg)]

  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    radht,        & ! SW + LW heating rate               [K/s]
    Frad,         & ! Radiative flux (momentum levels)   [W/m^2]
    Frad_SW_up,   & ! SW radiative upwelling flux        [W/m^2]
    Frad_LW_up,   & ! LW radiative upwelling flux        [W/m^2]
    Frad_SW_down, & ! SW radiative downwelling flux      [W/m^2]
    Frad_LW_down    ! LW radiative downwelling flux      [W/m^2]

  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    wpsclrp_sfc,    & ! Passive scalar flux at surface                [{units vary} m/s]
    wpedsclrp_sfc     ! Eddy-diffusion passive scalar flux at surface [{units vary} m/s]

  real( kind = core_rknd ), dimension(:,:,:), allocatable :: &
    sclrm,          & ! Passive scalar mean (thermo. levels) [{units vary}]
    wpsclrp,        & ! w'sclr' (momentum levels)            [{units vary} m/s]
    sclrp2,         & ! sclr'^2 (momentum levels)            [{units vary}^2]
    sclrp3,         & ! sclr'^3 (thermodynamic levels)       [{units vary}^3]
    sclrprtp,       & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
    sclrpthlp,      & ! sclr'thl' (momentum levels)          [{units vary} K]
    sclrpthvp,      & ! < sclr' th_v' > (momentum levels)   [units vary]
    sclrm_forcing,  & ! Passive scalar forcing          [{units vary}/s]
    edsclrm,        & ! Eddy passive scalar grid-mean (thermo. levels)   [units vary]
    edsclrm_forcing   ! Eddy-diffusion passive scalar forcing    [{units vary}/s]
    
  real( kind = core_rknd ), dimension(:,:,:), allocatable :: &
    lh_rt_clipped,          & ! rt generated from silhs sample points
    lh_thl_clipped,         & ! thl generated from silhs sample points
    lh_rc_clipped,          & ! rc generated from silhs sample points
    lh_rv_clipped,          & ! rv generated from silhs sample points
    lh_Nc_clipped,          & ! Nc generated from silhs sample points
    lh_sample_point_weights   ! Weights for cloud weighted sampling

  real( kind = core_rknd ), dimension(:,:,:), allocatable :: &
    K_hm,                       & ! Eddy diffusivity coef. for hydrometeors on mom. levs. [m^2 s^-1]
    wp2hmp,                     & ! Third moment:  <w'^2> * <hm'>       [(m/s)^2 <hm units>]
    rtphmp_zt,                  & ! Covariance of rt and a hydrometeor  [(kg/kg) <hm units>]
    thlphmp_zt,                 & ! Covariance of thl and a hydrometeor [K <hm units>]
    hydromet,                   & ! Array of hydrometeors                [hm units]
    hydrometp2,                 & ! Variance of a hydrometeor (m-levs.)  [<hm units>^2]
    wphydrometp,                & ! Covariance of w and a hydrometeor    [(m/s) <hm units>]
    hydromet_mc,                & ! Microphysics tendency for mean hydrometeors  [units/s]
    hydromet_vel_zt,            & ! Mean hydrometeor sed. velocity on thermo. levs. [m/s]
    hydromet_vel_covar_zt_impc, & ! Imp. comp. <V_xx'x_x'> t-levs [m/s]
    hydromet_vel_covar_zt_expc    ! Exp. comp. <V_xx'x_x'> t-levs [units(m/s)]

  real( kind = core_rknd ), dimension(:,:,:), allocatable :: &
    mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
    mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
    sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
    sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

  real( kind = core_rknd ), dimension(:,:,:,:), allocatable :: &
    corr_array_1_n,       & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
    corr_array_2_n,       & ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]
    corr_cholesky_mtx_1,  & ! Transposed corr. cholesky matrix, 1st comp. [-]
    corr_cholesky_mtx_2     ! Transposed corr. cholesky matrix, 2nd comp. [-]

  real( kind = core_rknd ), dimension(:,:,:,:), allocatable :: &
    X_nl_all_levs    ! Lognormally distributed hydrometeors

  !$omp threadprivate( l_restart, l_input_fields, l_modify_ic_with_cubic_int, &
  !$omp  l_modify_bc_for_cnvg_test, l_stats, l_silhs_out, l_rad_itime, l_first_write_to_file, &
  !$omp  l_allow_small_stats_tout, &
  !$omp  l_restart_input, l_last_timestep, ifinal, stats_nsamp, stats_nout, iunit, &
  !$omp  iunit_grid_adaptation, hydromet_dim, sclr_dim, edsclr_dim, pdf_dim, nzmax, &
  !$omp  grid_type, extended_atmos_bottom_level, extended_atmos_top_level, &
  !$omp  extended_atmos_range_size, day, month, year, lin_int_buffer, sfctype, stats_metadata, &
  !$omp  sclr_idx, hm_metadata, clubb_config_flags, silhs_config_flags, stats_tsamp, stats_tout, &
  !$omp  time_restart, restart_path_case, forcings_file_path, stats_fmt, fname_prefix, &
  !$omp  case_info_file, fname_grid_adaptation, output_file_prefix, runtype, zt_grid_fname, &
  !$omp  zm_grid_fname, lat_vals, lon_vals, time_initial, time_final, time_current, dt_main, &
  !$omp  dt_rad, deltaz_nl, zm_init_nl, zm_top_nl, mixt_frac_max_mag, T0, ts_nudge, rtm_min, &
  !$omp  rtm_nudge_max_altitude, rad_dummy, deep_soil_T_in_K_init, sfc_soil_T_in_K_init, &
  !$omp  veg_T_in_K_init, dummy_dx, dummy_dy, sclr_tol, thlm_init, rtm_init, um_init, &
  !$omp  vm_init, ug_init, vg_init, wp2_init, up2_init, vp2_init, upwp_init, rcm_init, wm_zt_init, &
  !$omp  wm_zm_init, em_init, exner_init, thvm_init, p_in_Pa_init, rho_init, rho_zm_init, &
  !$omp  rho_ds_zm_init, rho_ds_zt_init, invrs_rho_ds_zm_init, invrs_rho_ds_zt_init, &
  !$omp  thv_ds_zm_init, thv_ds_zt_init, rtm_ref_init, thlm_ref_init, um_ref_init, &
  !$omp  vm_ref_init, Ncm_init, Nc_in_cloud_init, Nccnm_init, rho_ds_zm_dycore_init, &
  !$omp  p_sfc, T_sfc, fcor, fcor_y, sfc_elevation, zm_init, zm_top, deltaz, sclrm_init, &
  !$omp  edsclrm_init, corr_array_n_cloud, corr_array_n_below, iinit, X_mixt_comp_all_levs, &
  !$omp  gr, gr_desc, gr_dycore, gr_output, gr_fixed_min, pdf_params, pdf_params_zm, &
  !$omp  pdf_implicit_coefs_terms, precip_fracs, nu_vert_res_dep, stats_zt, stats_zm, &
  !$omp  stats_lh_zt, stats_lh_sfc, stats_rad_zt, stats_rad_zm, stats_sfc, hydromet_pdf_params, &
  !$omp  sens_ht, latent_ht, lmin, vert_decorr_coef, deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K, &
  !$omp  wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, upwp_sfc_pert, vpwp_sfc_pert, gr_dens_z, &
  !$omp  gr_dens, gr_dens_old_z, gr_dens_old, Lscale_counter, richardson_num_counter, chi_counter, &
  !$omp  cumulative_Lscale, cumulative_richardson_num, cumulative_chi, alt_term, lscale_term, &
  !$omp  lscale_term_time_avg, chi_term, chi_term_time_avg, richardson_num_term, &
  !$omp  richardson_num_term_time_avg, norm_min_grid_dens, norm_grid_dens, um, upwp, vm, vpwp, &
  !$omp  up2, vp2, up3, vp3, wprtp, wpthlp, rtp2, rtp3, thlp2, thlp3, rtpthlp, wp2, wp2_zt, thvm, &
  !$omp  exner, rtm, thlm, rcm, wp3, wp3_zm, delta_zm, p_in_Pa, cloud_frac, wpthvp, wp2thvp, &
  !$omp  rtpthvp, thlpthvp, wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, wpup2, wpvp2, wp2up2, &
  !$omp  wp2vp2, rho_ds_zt, thv_ds_zt, wm_zm, wm_zt, rho, rho_zm, rho_ds_zm, invrs_rho_ds_zm, &
  !$omp  invrs_rho_ds_zt, thv_ds_zm, thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &
  !$omp  wprtp_forcing, wpthlp_forcing, rtp2_forcing, thlp2_forcing, rtpthlp_forcing, um_pert, &
  !$omp  vm_pert, upwp_pert, vpwp_pert, Skw_zm, Skw_zm_smooth, Lscale, thlprcp, sigma_sqd_w, &
  !$omp  sigma_sqd_w_zt, wprcp, w_up_in_cloud, w_down_in_cloud, cloudy_updraft_frac, &
  !$omp  cloudy_downdraft_frac, ice_supersat_frac, rcm_in_layer, cloud_cover, invrs_tau_zm, &
  !$omp  ug, vg, rtm_ref, thlm_ref, um_ref, vm_ref, rcm_mc, rvm_mc, thlm_mc, wprtp_mc, wpthlp_mc, &
  !$omp  rtp2_mc, thlp2_mc, rtpthlp_mc, Ncm_mc, Ncm, Nccnm, em, tau_zm, tau_zt, Kh_zt, Kh_zm, rfrzm, &
  !$omp  rrm, Nrm, Nc_in_cloud, wpNcp, radht, Frad, Frad_SW_up, Frad_LW_up, Frad_SW_down, Frad_LW_down, &
  !$omp  wpsclrp_sfc, wpedsclrp_sfc, sclrm, wpsclrp, sclrp2, sclrp3, sclrprtp, sclrpthlp, sclrpthvp, &
  !$omp  sclrm_forcing, edsclrm, edsclrm_forcing, lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, &
  !$omp  lh_rv_clipped, lh_Nc_clipped, lh_sample_point_weights, K_hm, wp2hmp, rtphmp_zt, &
  !$omp  thlphmp_zt, hydromet, hydrometp2, wphydrometp, hydromet_mc, hydromet_vel_zt, &
  !$omp  hydromet_vel_covar_zt_impc, hydromet_vel_covar_zt_expc, mu_x_1_n, mu_x_2_n, sigma_x_1_n, &
  !$omp  sigma_x_2_n, corr_array_1_n, corr_array_2_n, corr_cholesky_mtx_1, corr_cholesky_mtx_2, &
  !$omp  X_nl_all_levs )

  contains

  !=============================================================================================
  !                                       run_clubb
  !=============================================================================================
  subroutine run_clubb ( ngrdcol, calls_per_out, l_output_multi_col, l_output_double_prec, &
                          clubb_params, runfile, l_stdout, err_info, &
                          model_flags_array )
  !
  ! Description:
  !   Perform a standard run of CLUBB, defined mainly by the contents of runfile
  !
  !---------------------------------------------------------------------------------------------
            
    implicit none

    !----------------------------------- Input Variables -----------------------------------
    integer, intent(in) :: &
      ngrdcol, &
      calls_per_out

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params  ! Model parameters, C1, nu2, etc.

    logical, intent(in) ::  & 
      l_stdout,             & ! Whether to print output per timestep
      l_output_multi_col,   & ! Determines whether mutlicolumn data is saved
      l_output_double_prec    ! Flag to enable double precision

    ! Subroutine Arguments (Model Setting)
    character(len=*), intent(in) :: &
      runfile ! Name of file containing &model_setting and &sounding

    !----------------------------------- Output Variables -----------------------------------
    type(err_info_type), intent(out) :: &
      err_info        ! err_info struct containing err_code and err_header

    !----------------------------------- Optional Input Variables -----------------------------------
    logical, optional, dimension(:), intent(in) :: &
      model_flags_array ! Array containing model flags (for the clubb_tuner only)
      
    !----------------------------------- Begin code -----------------------------------

    ! Read namelist, initialize all settings, and allocate all variables
    call init_clubb_case(runfile, ngrdcol, clubb_params, err_info, model_flags_array)

    ! Run the case to completion, advancing all fields to the time set in init_clubb_case
    call advance_clubb_to_end( ngrdcol, calls_per_out, clubb_params, &
                          l_stdout, l_output_multi_col, l_output_double_prec, &
                          err_info )

    ! Deallcoate and cleanup all memory
    call clean_up_clubb(ngrdcol, clubb_params, err_info)

    return

  end subroutine run_clubb

  !=============================================================================================
  !                                       init_clubb_case
  !=============================================================================================
  subroutine init_clubb_case(runfile, ngrdcol, clubb_params, err_info, model_flags_array)
  !
  ! Description:
  !   This subroutine sets reads the runfile namelist, sets up all model settings, allocates
  !   all variables used on the model at the top level, then calls 'set_case_initial_conditions'
  !   to apply the initialization values to the fields variables that are advanced throughout. 
  !
  !---------------------------------------------------------------------------------------------
    
    use bugsrad_driver, only: &
      init_radiation !------------------------------- Subroutine

    use stats_clubb_utilities, only: & 
      stats_init_api, &
      stats_init_w_diff_output_gr
    
    use silhs_api_module, only: &
      latin_hypercube_2D_output_api

    use variables_radiation_module, only: &
      setup_radiation_variables    !---------------------------------------- Procedure(s)

    use soil_vegetation, only: &
      initialize_soil_veg      !--------------------------------------------- Procedure(s)

    use microphys_init_cleanup, only: &
      init_microphys     !------------------------------------------------- Procedure

#ifdef SILHS
    use simple_rad_module, only: &
      simple_rad_lba_init !---------------------- Procedure(s)
#endif

    use err_info_type_module, only: &
      set_err_info_values_api, &
      init_default_err_info_api

    use clubb_api_module, only: &
      init_pdf_implicit_coefs_terms_api, & ! Procedure(s)
      init_precip_fracs_api, &
      setup_grid_api, &
      set_clubb_debug_level_api, &
      print_clubb_config_flags_api, &
      check_clubb_settings_api, &
      clubb_at_least_debug_level_api, &
      set_default_clubb_config_flags_api, &  
      initialize_clubb_config_flags_type_api, &
      init_pdf_params_api, &
      check_parameters_api, &
      calc_derrived_params_api, &
      iiPDF_new, &        ! Constants
      iiPDF_new_hybrid

    use inputfields, only: &
      inputfields_init !------------ Procedure(s)

    use grid_class, only: &
      read_grid_heights, &   !---------------------------------------------- Procedure(s)
      flip

    use grid_adaptation_module, only: &
      setup_gr_fixed_min, &
      setup_gr_dycore

    implicit none

    character(len=*), intent(in) :: &
      runfile ! Name of file containing &model_setting and &sounding

    integer, intent(in) :: &
      ngrdcol

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params  ! Model parameters, C1, nu2, etc.

    logical, optional, dimension(:), intent(in) :: &
      model_flags_array ! Array containing model flags (for the clubb_tuner only)

    !----------------------------------- Output Variables -----------------------------------
    type(err_info_type), intent(out) :: &
      err_info        ! err_info struct containing err_code and err_header

    !----------------------------------- Local Variables -----------------------------------

    integer :: &
      iiPDF_type,                     & ! Selected option for the two-component normal
                                        ! (double Gaussian) PDF type to use for the w, rt,
                                        ! and theta-l (or w, chi, and eta) portion of
                                        ! CLUBB's multivariate, two-component PDF.
      ipdf_call_placement,            & ! Selected option for the placement of the call to
                                        ! CLUBB's PDF.
      penta_solve_method,             & ! Option to set the penta-diagonal matrix solving method
      tridiag_solve_method,           & ! Option to set the tri-diagonal matrix solving method
      saturation_formula,             & ! Integer that stores the saturation formula to be used
      grid_remap_method,              & ! Integer that stores what remapping technique should
                                        ! be used to remap the values from one grid to another
                                        ! (starts with 1, so 0 is invalid value for flag)
      grid_adapt_in_time_method,      & ! Integer that stores how the grid should be adapted every
                                        ! timestep or if the grid should not be adapted at all
      fill_holes_type                   ! Option for which type of hole filler to use in the 
                                        ! fill_holes_vertical procedure

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
      l_ho_nontrad_coriolis,        & ! Flag to implement the nontraditional Coriolis terms in the
                                      ! prognostic equations of <w'w'>, <u'w'>, and <u'u'>.
      l_ho_trad_coriolis,           & ! Flag to implement the traditional Coriolis terms in the
                                      ! prognostic equations of <v'w'> and <u'w'>.
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
      l_damp_wp3_Skw_squared,       & ! Set damping on wp3 to use Skw^2 rather than Skw^4
      l_prescribed_avg_deltaz,      & ! used in calc_derrived_params_api. If .true., avg_deltaz = deltaz
      l_lmm_stepping,               & ! Apply Linear Multistep Method (LMM) Stepping
      l_e3sm_config,                & ! Run model with E3SM settings
      l_vary_convect_depth,         & ! Flag used to calculate convective velocity using
                                      ! a variable estimate of layer depth based on the depth
                                      ! over which wpthlp is positive near the ground when true
                                      ! More information can be found by
                                      ! Looking at issue #905 on the clubb repo
      l_use_tke_in_wp3_pr_turb_term,& ! Use TKE formulation for wp3 pr_turb term
      l_use_tke_in_wp2_wp3_K_dfsn,  & ! Use TKE in eddy diffusion for wp2 and wp3
      l_use_wp3_lim_with_smth_Heaviside, & ! Flag to activate mods on wp3 limiters for conv test
      l_smooth_Heaviside_tau_wpxp,  & ! Use smoothed Heaviside 'Preskin' function
                                      ! in the calculation of H_invrs_tau_wpxp_N2
                                      ! in src/CLUBB_core/mixing_length.F90
      l_modify_limiters_for_cnvg_test, & ! Flag to activate mods on limiters for conv test
      l_enable_relaxed_clipping,    & ! Flag to relax clipping on wpxp in
                                      ! xm_wpxp_clipping_and_stats
      l_linearize_pbl_winds,        & ! Code to linearize PBL winds
      l_mono_flux_lim_thlm,         & ! Flag to turn on monotonic flux limiter for thlm
      l_mono_flux_lim_rtm,          & ! Flag to turn on monotonic flux limiter for rtm
      l_mono_flux_lim_um,           & ! Flag to turn on monotonic flux limiter for um
      l_mono_flux_lim_vm,           & ! Flag to turn on monotonic flux limiter for vm
      l_mono_flux_lim_spikefix,     & ! Flag to implement monotonic flux limiter code that
                                      ! eliminates spurious drying tendencies at model top
      l_host_applies_sfc_fluxes,    & ! Use to determine whether a host model has already applied 
                                      ! the surface flux, to avoid double counting.
      l_wp2_fill_holes_tke,         & ! Turn on additional hole-filling for wp2
                                      ! that takes TKE from up2 and vp2, if necessary
      l_add_dycore_grid               ! Turn on remapping values from a dycore grid

    integer :: &
      debug_level,        & ! Amount of debugging information
      iisclr_rt    = -1,  & ! Scalar indices
      iisclr_thl   = -1,  &
      iisclr_CO2   = -1,  &
      iiedsclr_rt  = -1,  &
      iiedsclr_thl = -1,  &
      iiedsclr_CO2 = -1,  &
      i, j, k               ! Loop variables
      
    real( kind = core_rknd ), dimension(sclr_max) :: &
      sclr_tol_nl        ! Thresholds on the passive scalars     [units vary]

    real( kind = core_rknd ) :: &
      sfc_elevation_nl, &
      p_sfc_nl, &
      T_sfc_nl, &
      fcor_nl
    
    ! Grid altitude arrays
    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      momentum_heights, thermodynamic_heights ! [m]
      
    ! Brian
    real( kind = core_rknd ), dimension(:,:), allocatable ::  & 
      momentum_heights_flip, thermodynamic_heights_flip ! [m]
    ! remove above

    ! Definition of namelists
    namelist /model_setting/ &
      runtype, nzmax, grid_type, deltaz_nl, zm_init_nl, zm_top_nl, &
      zt_grid_fname, zm_grid_fname, &
      day, month, year, lat_vals, lon_vals, sfc_elevation_nl, &
      time_initial, time_final, &
      dt_main, dt_rad, &
      sfctype, T_sfc_nl, p_sfc_nl, sens_ht, latent_ht, fcor_nl, T0, ts_nudge, &
      forcings_file_path, l_t_dependent, l_input_xpwp_sfc, &
      l_ignore_forcings, l_modify_ic_with_cubic_int, &
      l_modify_bc_for_cnvg_test, &
      thlm_sponge_damp_settings, rtm_sponge_damp_settings, &
      uv_sponge_damp_settings, wp2_sponge_damp_settings, &
      wp3_sponge_damp_settings, up2_vp2_sponge_damp_settings, &
      l_soil_veg, l_uv_nudge, l_restart, restart_path_case, &
      time_restart, l_input_fields, debug_level, &
      sclr_tol_nl, sclr_dim, iisclr_thl, iisclr_rt, iisclr_CO2, &
      edsclr_dim, iiedsclr_thl, iiedsclr_rt, iiedsclr_CO2, &
      l_rtm_nudge, rtm_min, rtm_nudge_max_altitude, &
      l_diagnose_correlations, l_calc_w_corr

  namelist /stats_setting/ &
    l_stats, fname_prefix, stats_tsamp, stats_tout, stats_fmt, &
    l_allow_small_stats_tout

  namelist /configurable_clubb_flags_nl/ &
    iiPDF_type, ipdf_call_placement, penta_solve_method, tridiag_solve_method, &
    saturation_formula, grid_remap_method, &
    grid_adapt_in_time_method, fill_holes_type, &
    l_upwind_xpyp_ta, l_upwind_xm_ma, &
    l_tke_aniso, l_vert_avg_closure, l_standard_term_ta, &
    l_partial_upwind_wp3, l_godunov_upwind_wpxp_ta, l_godunov_upwind_xpyp_ta, &
    l_use_cloud_cover, l_rcm_supersat_adj, &
    l_damp_wp3_Skw_squared, l_min_wp2_from_corr_wx, l_min_xp2_from_corr_wx, &
    l_C2_cloud_frac, l_predict_upwp_vpwp, l_ho_nontrad_coriolis, l_ho_trad_coriolis, &
    l_diag_Lscale_from_tau, &
    l_stability_correct_tau_zm, l_damp_wp2_using_em, l_use_C7_Richardson, &
    l_use_precip_frac, l_do_expldiff_rtm_thlm, l_use_C11_Richardson, &
    l_use_shear_Richardson, l_prescribed_avg_deltaz, l_diffuse_rtm_and_thlm, &
    l_stability_correct_Kh_N2_zm, l_trapezoidal_rule_zt, l_trapezoidal_rule_zm, &
    l_call_pdf_closure_twice, l_Lscale_plume_centered, &
    l_brunt_vaisala_freq_moist, l_use_thvm_in_bv_freq, &
    l_lmm_stepping, l_e3sm_config, l_vary_convect_depth, l_use_tke_in_wp3_pr_turb_term, &
    l_use_tke_in_wp2_wp3_K_dfsn, l_use_wp3_lim_with_smth_Heaviside, l_smooth_Heaviside_tau_wpxp, &
    l_modify_limiters_for_cnvg_test, l_enable_relaxed_clipping, &
    l_linearize_pbl_winds, l_mono_flux_lim_thlm, &
    l_mono_flux_lim_rtm, l_mono_flux_lim_um, l_mono_flux_lim_vm, l_mono_flux_lim_spikefix, &
    l_host_applies_sfc_fluxes, l_wp2_fill_holes_tke, &
    l_add_dycore_grid

    !----------------------------------- Begin Code -----------------------------------
    
    ! Pick some default values for model_setting.  Some variables are initialized
    ! at declaration in module clubb_model_settings.

    T_sfc_nl  = 288._core_rknd
    p_sfc_nl  = 1000.e2_core_rknd
    sens_ht   = 0._core_rknd
    latent_ht = 0._core_rknd
    fcor_nl   = 1.e-4_core_rknd
    
    T0        = 300._core_rknd
    ts_nudge  = 86400._core_rknd
    rtm_min   = epsilon( rtm_min )
    rtm_nudge_max_altitude = 10000._core_rknd

    nzmax     = 75
    grid_type = 1

    extended_atmos_bottom_level = 0
    extended_atmos_top_level    = 0
    extended_atmos_range_size   = 0

    lin_int_buffer = 0

    deltaz_nl  = 40._core_rknd
    zm_init_nl = 0._core_rknd
    zm_top_nl  = 3500._core_rknd

    zt_grid_fname = ""
    zm_grid_fname = ""

    day = 22; month = 6; year = 1969

    lat_vals = 15._core_rknd
    lon_vals = -56.5_core_rknd

    sfctype = 0 

    time_initial = 0._time_precision
    time_final   = 21600._time_precision
    time_current = 0._time_precision

    dt_main = 60._core_rknd
    dt_rad  = 600._core_rknd

    forcings_file_path = ''
    l_t_dependent   = .false. 
    l_input_xpwp_sfc = .false. 
    l_ignore_forcings = .false. 
    l_modify_ic_with_cubic_int = .false.
    l_modify_bc_for_cnvg_test = .false.

    l_sfc_already_initialized = .false.

    thlm_sponge_damp_settings%l_sponge_damping = .false.
    rtm_sponge_damp_settings%l_sponge_damping = .false.
    uv_sponge_damp_settings%l_sponge_damping = .false.
    wp2_sponge_damp_settings%l_sponge_damping = .false.
    wp3_sponge_damp_settings%l_sponge_damping = .false.
    up2_vp2_sponge_damp_settings%l_sponge_damping = .false.

    thlm_sponge_damp_settings%tau_sponge_damp_min = 60._core_rknd
    thlm_sponge_damp_settings%tau_sponge_damp_max = 1800._core_rknd
    thlm_sponge_damp_settings%sponge_damp_depth = 0.25_core_rknd

    rtm_sponge_damp_settings%tau_sponge_damp_min = 60._core_rknd
    rtm_sponge_damp_settings%tau_sponge_damp_max = 1800._core_rknd
    rtm_sponge_damp_settings%sponge_damp_depth = 0.25_core_rknd

    uv_sponge_damp_settings%tau_sponge_damp_min = 60._core_rknd
    uv_sponge_damp_settings%tau_sponge_damp_max = 1800._core_rknd
    uv_sponge_damp_settings%sponge_damp_depth = 0.25_core_rknd

    wp2_sponge_damp_settings%tau_sponge_damp_min = 60._core_rknd
    wp2_sponge_damp_settings%tau_sponge_damp_max = 1800._core_rknd
    wp2_sponge_damp_settings%sponge_damp_depth = 0.25_core_rknd

    wp3_sponge_damp_settings%tau_sponge_damp_min = 60._core_rknd
    wp3_sponge_damp_settings%tau_sponge_damp_max = 1800._core_rknd
    wp3_sponge_damp_settings%sponge_damp_depth = 0.25_core_rknd

    up2_vp2_sponge_damp_settings%tau_sponge_damp_min = 60._core_rknd
    up2_vp2_sponge_damp_settings%tau_sponge_damp_max = 1800._core_rknd
    up2_vp2_sponge_damp_settings%sponge_damp_depth = 0.25_core_rknd

    l_restart      = .false.
    l_input_fields  = .false.
    restart_path_case = "none"
    time_restart  = 0._time_precision
    debug_level   = 2

    ! Use the Flatau polynomial approximation for computing saturation in advance_clubb_core_module
    saturation_formula = saturation_flatau

    sclr_dim   = 0
    sclr_idx%iisclr_thl = -1
    sclr_idx%iisclr_rt  = -1
    sclr_idx%iisclr_CO2 = -1

    edsclr_dim = 0
    sclr_idx%iiedsclr_thl = -1
    sclr_idx%iiedsclr_rt  = -1
    sclr_idx%iiedsclr_CO2 = -1

    l_first_write_to_file = .true.

    ! Pick some default values for stats_setting; other variables are set in
    ! module stats_variables
    fname_prefix = ''
    stats_fmt    = ''

    ! Initialize err_info with default values for one column
    call init_default_err_info_api(ngrdcol, err_info)

    ! Fill err_info with dummy values
    ! The reshapes create arrays of dimension(ngrdcol) filled with values
    ! lat_vals and lon_vals, respectively
    ! Populating err_info with the available info
    ! We don't want to populate the parallelization info
    ! since that is only available outside of run_clubb
#ifdef _OPENMP
    iunit = omp_get_thread_num( ) + 10 ! Known magic number

    call set_err_info_values_api( ngrdcol, err_info, &
                                  chunk_idx_in=omp_get_thread_num(), &
                                  lat_in = reshape((/lat_vals/), (/ngrdcol/), (/lat_vals/)), &
                                  lon_in = reshape((/lon_vals/), (/ngrdcol/), (/lon_vals/)) )
#else
    iunit = 10

    call set_err_info_values_api( ngrdcol, err_info, &
                                  lat_in = reshape((/lat_vals/), (/ngrdcol/), (/lat_vals/)), &
                                  lon_in = reshape((/lon_vals/), (/ngrdcol/), (/lon_vals/)) )
#endif

    ! Default values for the soil scheme
    call initialize_soil_veg( 1, deep_soil_T_in_K_init, sfc_soil_T_in_K_init, veg_T_in_K_init )

    call set_default_clubb_config_flags_api( iiPDF_type, & ! Intent(out)
                                            ipdf_call_placement, & ! Intent(out)
                                            penta_solve_method, & ! Intent(out)
                                            tridiag_solve_method, & ! Intent(out)
                                            saturation_formula, &  ! Intent(out)
                                            grid_remap_method, & ! Intent(out)
                                            grid_adapt_in_time_method, & ! Intent(out)
                                            fill_holes_type, & ! Intent(out)
                                            l_use_precip_frac, & ! Intent(out)
                                            l_predict_upwp_vpwp, & ! Intent(out)
                                            l_ho_nontrad_coriolis, & ! Intent(out)
                                            l_ho_trad_coriolis, & ! Intent(out)
                                            l_min_wp2_from_corr_wx, & ! Intent(out)
                                            l_min_xp2_from_corr_wx, & ! Intent(out)
                                            l_C2_cloud_frac, & ! Intent(out)
                                            l_diffuse_rtm_and_thlm, & ! Intent(out)
                                            l_stability_correct_Kh_N2_zm, & ! Intent(out)
                                            l_calc_thlp2_rad, & ! Intent(out)
                                            l_upwind_xpyp_ta, & ! Intent(out)
                                            l_upwind_xm_ma, & ! Intent(out)
                                            l_uv_nudge, & ! Intent(out)
                                            l_rtm_nudge, & ! Intent(out)
                                            l_tke_aniso, & ! Intent(out)
                                            l_vert_avg_closure, & ! Intent(out)
                                            l_trapezoidal_rule_zt, & ! Intent(out)
                                            l_trapezoidal_rule_zm, & ! Intent(out)
                                            l_call_pdf_closure_twice, & ! Intent(out)
                                            l_standard_term_ta, & ! Intent(out)
                                            l_partial_upwind_wp3, & ! Intent(out)
                                            l_godunov_upwind_wpxp_ta, & ! Intent(out)
                                            l_godunov_upwind_xpyp_ta, & ! Intent(out)
                                            l_use_cloud_cover, & ! Intent(out)
                                            l_diagnose_correlations, & ! Intent(out)
                                            l_calc_w_corr, & ! Intent(out)
                                            l_const_Nc_in_cloud, & ! Intent(out)
                                            l_fix_w_chi_eta_correlations, & ! Intent(out)
                                            l_stability_correct_tau_zm, & ! Intent(out)
                                            l_damp_wp2_using_em, & ! Intent(out)
                                            l_do_expldiff_rtm_thlm, & ! Intent(out)
                                            l_Lscale_plume_centered, & ! Intent(out)
                                            l_diag_Lscale_from_tau, & ! Intent(out)
                                            l_use_C7_Richardson, & ! Intent(out)
                                            l_use_C11_Richardson, & ! Intent(out)
                                            l_use_shear_Richardson, & ! Intent(out)
                                            l_brunt_vaisala_freq_moist, & ! Intent(out)
                                            l_use_thvm_in_bv_freq, & ! Intent(out)
                                            l_rcm_supersat_adj, & ! Intent(out)
                                            l_damp_wp3_Skw_squared, & ! Intent(out)
                                            l_prescribed_avg_deltaz, & ! Intent(out)
                                            l_lmm_stepping, & ! Intent(out)
                                            l_e3sm_config, & ! Intent(out)
                                            l_vary_convect_depth, & ! Intent(out)
                                            l_use_tke_in_wp3_pr_turb_term, & ! Intent(out)
                                            l_use_tke_in_wp2_wp3_K_dfsn, & ! Intent(out)
                                            l_use_wp3_lim_with_smth_Heaviside, & ! Intent(out)
                                            l_smooth_Heaviside_tau_wpxp, & ! Intent(out)
                                            l_modify_limiters_for_cnvg_test, & ! Intent(out)
                                            l_enable_relaxed_clipping, & ! Intent(out)
                                            l_linearize_pbl_winds, & ! Intent(out)
                                            l_mono_flux_lim_thlm, & ! Intent(out)
                                            l_mono_flux_lim_rtm, & ! Intent(out)
                                            l_mono_flux_lim_um, & ! Intent(out)
                                            l_mono_flux_lim_vm, & ! Intent(out)
                                            l_mono_flux_lim_spikefix, & ! Intent(out)
                                            l_host_applies_sfc_fluxes, & ! Intent(out)
                                            l_wp2_fill_holes_tke, & ! Intent(out)
                                            l_add_dycore_grid ) ! Intent(out)

    ! Read namelist file
    open(unit=iunit, file=trim( runfile ), status='old')
    read(unit=iunit, nml=model_setting)
    read(unit=iunit, nml=stats_setting)
    close(unit=iunit)

    sclr_idx%iisclr_thl = iisclr_thl
    sclr_idx%iisclr_rt  = iisclr_rt
    sclr_idx%iisclr_CO2 = iisclr_CO2

    sclr_idx%iiedsclr_thl = iiedsclr_thl
    sclr_idx%iiedsclr_rt  = iiedsclr_rt
    sclr_idx%iiedsclr_CO2 = iiedsclr_CO2

    allocate( sclr_tol(sclr_dim) )
    sclr_tol(1:sclr_dim) = sclr_tol_nl(1:sclr_dim)

    stats_metadata%l_allow_small_stats_tout = l_allow_small_stats_tout

    allocate( &
      stats_zt(ngrdcol),      &
      stats_zm(ngrdcol),      &
      stats_lh_zt(ngrdcol),   &
      stats_lh_sfc(ngrdcol),  &
      stats_rad_zt(ngrdcol),  &
      stats_rad_zm(ngrdcol),  &
      stats_sfc(ngrdcol),     &
      p_sfc(ngrdcol),         &
      T_sfc(ngrdcol),         &
      fcor(ngrdcol),          &
      fcor_y(ngrdcol),        &
      dummy_dx(ngrdcol),      &
      dummy_dy(ngrdcol),      &
      deep_soil_T_in_K(ngrdcol), &
      sfc_soil_T_in_K(ngrdcol),  &
      veg_T_in_K(ngrdcol),    &
      upwp_sfc_pert(ngrdcol), &
      vpwp_sfc_pert(ngrdcol), &
      sfc_elevation(ngrdcol), &
      zm_init(ngrdcol),       &
      zm_top(ngrdcol),        &
      deltaz(ngrdcol)         )

    sfc_elevation = sfc_elevation_nl
    deltaz = deltaz_nl
    zm_init = zm_init_nl
    zm_top = zm_top_nl
    p_sfc = p_sfc_nl
    T_sfc = T_sfc_nl
    fcor = fcor_nl

    ! Suppose the vector Coriolis term is written as 2*Omega_vector x v_vector, where 
    ! Omega_vector is the Earth's angular velocity vector (omega_planet in code), which 
    ! points along the North Pole, and v_vector is the full air velocity vector in the 
    ! rotating frame of reference. Now suppose we break Omega_vector into a local upward
    ! component, Omega_z*z_hat, and a local northward component, Omega_y*y_hat. 
    ! Then we define fcor = 2*Omega_z and fcor_y = 2*Omega_y. 
    ! Then fcor_y = 2 * Omega_vector \dot y_hat = 2 * |Omega_vector| * cos( angle between north 
    ! pole and local northward direction ) 
    fcor_y = two * omega_planet * cos ( lat_vals * radians_per_deg )

    ! fcor could reasonably be defined as below, but currently it is input via namelist
    ! Hing Ong, 25 November 2025
    !fcor   = two * omega_planet * sin ( lat_vals * radians_per_deg )

    open(unit=iunit, file=runfile, status='old', action='read')
    read(unit=iunit, nml=configurable_clubb_flags_nl)
    close(unit=iunit)

    if ( l_vert_avg_closure ) then
      l_trapezoidal_rule_zt    = .true.
      l_trapezoidal_rule_zm    = .true.
      l_call_pdf_closure_twice = .true.
    end if

    ! The case prefix including output di
    output_file_prefix  = trim( output_dir ) // trim( fname_prefix )        

    ! The filename for case setup
    case_info_file      = trim( output_dir ) // trim( fname_prefix ) // "_setup.txt"

    ! Sanity check on passive scalars
    ! When adding new 'ii' scalar indices, add them to this list.
    if ( max( iisclr_CO2, iisclr_rt, iisclr_thl ) > sclr_dim ) then
      write(fstderr, *) err_info%err_header_global
      write(fstderr,*) "Passive scalar index exceeds sclr_dim ", &
        "iisclr_CO2 = ", iisclr_CO2, "iisclr_rt = ", iisclr_rt,  &
        "iisclr_thl = ", iisclr_thl, "sclr_dim = ", sclr_dim

      ! General error -> set all entries to clubb_fatal_error
      err_info%err_code = clubb_fatal_error
      return

    else if ( max( iiedsclr_CO2, iiedsclr_rt, iiedsclr_thl ) > edsclr_dim ) then
      write(fstderr, *) err_info%err_header_global
      write(fstderr,*) "Passive scalar index exceeds edsclr_dim ", &
        "iiedsclr_CO2 = ", iiedsclr_CO2, "iiedsclr_rt = ", iiedsclr_rt, &
        "iiedsclr_thl = ", iiedsclr_thl, "edsclr_dim = ", edsclr_dim

      ! General error -> set all entries to clubb_fatal_error
      err_info%err_code = clubb_fatal_error
      return

    end if

#ifdef _OPENACC
    if ( l_restart ) stop "l_restart with OpenACC is not enabled"
#endif

    ! Set debug level
    call set_clubb_debug_level_api( debug_level ) ! Intent(in)

    ! Printing Model Inputs
    if ( clubb_at_least_debug_level_api( 1 ) ) then

      if ( l_write_to_file ) then
        open(unit=iunit, file=case_info_file, status='replace', action='write')
      end if

      ! Print the date and time
      call write_date( l_write_to_file, iunit )

      call write_text( "", l_write_to_file, iunit )
      call write_text( "--------------------------------------------------", &
                        l_write_to_file, iunit )
      call write_text( "Latest git log entry", l_write_to_file, iunit )
      call write_text( "--------------------------------------------------", &
                        l_write_to_file, iunit )
      call write_text( "", l_write_to_file, iunit )
      if ( l_write_to_file ) close( unit=iunit )

      if ( l_write_to_file ) then
        call execute_command_line( 'echo "Branch name:" >> '//case_info_file )
        call execute_command_line( 'git rev-parse --abbrev-ref HEAD >> '//case_info_file )
        call execute_command_line( 'echo "(If no branch name is shown HEAD may be detached)"'// &
                                    ' >> '//case_info_file )
        call execute_command_line( 'echo "" >> '//case_info_file )
        call execute_command_line( 'git log -n 1 >> '//case_info_file )
      end if

      if ( l_write_to_file ) then
        open(unit=iunit, file=case_info_file, status='OLD', action='write', position='APPEND')
      end if

      call write_text( NEW_LINE('A')//"A detailed git diff can be found "// &
                        "at the end of this file"//NEW_LINE('A'), &
                        l_write_to_file, iunit )

      call write_text( "--------------------------------------------------", &
        l_write_to_file, iunit )
      call write_text( "Tunable Parameters:", l_write_to_file, iunit)
      call write_text( "--------------------------------------------------", &
        l_write_to_file, iunit )


      call write_text( "--------------------------------------------------", &
        l_write_to_file, iunit )
      call write_text( "Preprocessing Directives:", l_write_to_file, iunit)
      call write_text( "--------------------------------------------------", &
        l_write_to_file, iunit )

#ifdef NETCDF
      call write_text( "-DNETCDF enabled", l_write_to_file, iunit )
#else
      call write_text( "-DNETCDF disabled", l_write_to_file, iunit )
#endif
#ifdef COAMPS_MICRO
      call write_text( "-DCOAMPS_MICRO enabled", l_write_to_file, iunit )
#else
      call write_text( "-DCOAMPS_MICRO disabled", l_write_to_file, iunit )
#endif
#ifdef TUNER
      call write_text( "-DTUNER enabled", l_write_to_file, iunit )
#else
      call write_text( "-DTUNER disabled", l_write_to_file, iunit )
#endif
#ifdef SILHS
      call write_text( "-DSILHS enabled", l_write_to_file, iunit )
#else
      call write_text( "-DSILHS disabled", l_write_to_file, iunit )
#endif
#ifdef nooverlap
      call write_text( "-Dnooverlap enabled", l_write_to_file, iunit )
#else
      call write_text( "-Dnooverlap disabled", l_write_to_file, iunit )
#endif
#ifdef radoffline
      call write_text( "-Draoffline enabled", l_write_to_file, iunit )
#else
      call write_text( "-Dradoffline disabled", l_write_to_file, iunit )
#endif
#ifdef USE_BUGSrad_ocast_random
      call write_text( "-DUSE_BUGSrad_ocast_random enabled", l_write_to_file, iunit )
#else
      call write_text( "-DUSE_BUGSrad_ocast_random disabled", l_write_to_file, iunit )
#endif
#ifdef BYTESWAP_IO
      call write_text( "-DBYTESWAP_IO enabled", l_write_to_file, iunit )
#else
      call write_text( "-DBYTESWAP_IO disabled", l_write_to_file, iunit )
#endif

      ! Pick some default values for model_setting
      call write_text( "--------------------------------------------------", &
        l_write_to_file, iunit )
      call write_text( "&model_setting", l_write_to_file, iunit)
      call write_text( "--------------------------------------------------", &
        l_write_to_file, iunit )
      call write_text( "runtype = " // runtype, l_write_to_file, iunit )
      call write_text( "nzmax = ", nzmax, l_write_to_file, iunit )
      call write_text( "grid_type = ", grid_type, l_write_to_file, iunit )
      call write_text( "deltaz = ", deltaz, l_write_to_file, iunit )
      call write_text( "zm_init = ", zm_init, l_write_to_file, iunit )
      call write_text( "zm_top = ", zm_top, l_write_to_file, iunit )
      call write_text( "zt_grid_fname = " // zt_grid_fname, l_write_to_file, iunit )
      call write_text( "zm_grid_fname = " // zm_grid_fname, l_write_to_file, iunit )

      call write_text( "day = ", day, l_write_to_file, iunit )
      call write_text( "month = ", month, l_write_to_file, iunit )
      call write_text( "year = ", year, l_write_to_file, iunit )

      call write_text( "lat_vals = ", lat_vals, l_write_to_file, iunit )
      call write_text( "lon_vals = ", lon_vals, l_write_to_file, iunit )

      call write_text( "sfc_elevation = ", sfc_elevation_nl, l_write_to_file, iunit )

      call write_text( "time_initial = ", real( time_initial, kind = core_rknd ), &
            l_write_to_file, iunit )
      call write_text( "time_final = ", real( time_final, kind = core_rknd ), &
            l_write_to_file, iunit )

      call write_text( "dt_main = ", dt_main, l_write_to_file, iunit )
      call write_text( "dt_rad = ",  dt_rad, l_write_to_file, iunit )

      call write_text( "sfctype = ", sfctype, l_write_to_file, iunit )
      call write_text( "T_sfc = ", T_sfc, l_write_to_file, iunit )
      call write_text( "p_sfc = ", p_sfc, l_write_to_file, iunit )
      call write_text( "sens_ht = ", sens_ht, l_write_to_file, iunit )
      call write_text( "latent_ht = ", latent_ht, l_write_to_file, iunit )
      call write_text( "fcor = ", fcor, l_write_to_file, iunit )
      call write_text( "fcor_y = ", fcor_y, l_write_to_file, iunit )
      call write_text( "T0 = ", T0, l_write_to_file, iunit )
      call write_text( "ts_nudge = ", ts_nudge, l_write_to_file, iunit )

      call write_text( "forcings_file_path = " // forcings_file_path, l_write_to_file, iunit )

      call write_text( "l_t_dependent = ", l_t_dependent, l_write_to_file, iunit )
      call write_text( "l_ignore_forcings = ", l_ignore_forcings, l_write_to_file, iunit )
      call write_text( "l_modify_ic_with_cubic_int = ", l_modify_ic_with_cubic_int, &
                        l_write_to_file, iunit )
      call write_text( "l_modify_bc_for_cnvg_test = ", l_modify_bc_for_cnvg_test, &
                        l_write_to_file, iunit )
      call write_text( "l_input_xpwp_sfc = ", l_input_xpwp_sfc, l_write_to_file, iunit )

      call write_text( "saturation_formula = ", saturation_formula, &
        l_write_to_file, iunit )

      call write_text( "thlm_sponge_damp_settings%l_sponge_damping = ", &
        thlm_sponge_damp_settings%l_sponge_damping, l_write_to_file, iunit )

      call write_text( "rtm_sponge_damp_settings%l_sponge_damping = ", &
        rtm_sponge_damp_settings%l_sponge_damping, l_write_to_file, iunit )

      call write_text( "uv_sponge_damp_settings%l_sponge_damping = ", &
        uv_sponge_damp_settings%l_sponge_damping, l_write_to_file, iunit )

      call write_text( "wp2_sponge_damp_settings%l_sponge_damping = ", &
        wp2_sponge_damp_settings%l_sponge_damping, l_write_to_file, iunit )

      call write_text( "wp3_sponge_damp_settings%l_sponge_damping = ", &
        wp3_sponge_damp_settings%l_sponge_damping, l_write_to_file, iunit )

      call write_text( "up2_vp2_sponge_damp_settings%l_sponge_damping = ", &
        up2_vp2_sponge_damp_settings%l_sponge_damping, l_write_to_file, iunit )

      call write_text( "thlm_sponge_damp_settings%tau_sponge_damp_min = ", &
        thlm_sponge_damp_settings%tau_sponge_damp_min, l_write_to_file, iunit )

      call write_text( "thlm_sponge_damp_settings%tau_sponge_damp_max = ", &
        thlm_sponge_damp_settings%tau_sponge_damp_max, l_write_to_file, iunit )

      call write_text( "thlm_sponge_damp_settings%sponge_damp_depth = ", &
        thlm_sponge_damp_settings%sponge_damp_depth, l_write_to_file, iunit )

      call write_text( "rtm_sponge_damp_settings%tau_sponge_damp_min = ", &
        rtm_sponge_damp_settings%tau_sponge_damp_min, l_write_to_file, iunit )

      call write_text( "rtm_sponge_damp_settings%tau_sponge_damp_max = ", &
        rtm_sponge_damp_settings%tau_sponge_damp_max, l_write_to_file, iunit )

      call write_text( "rtm_sponge_damp_settings%sponge_damp_depth = ", &
        rtm_sponge_damp_settings%sponge_damp_depth, l_write_to_file, iunit )

      call write_text( "uv_sponge_damp_settings%tau_sponge_damp_min = ", &
        uv_sponge_damp_settings%tau_sponge_damp_min, l_write_to_file, iunit )

      call write_text( "uv_sponge_damp_settings%tau_sponge_damp_max = ", &
        uv_sponge_damp_settings%tau_sponge_damp_max, l_write_to_file, iunit )

      call write_text( "uv_sponge_damp_settings%sponge_damp_depth = ", &
        uv_sponge_damp_settings%sponge_damp_depth, l_write_to_file, iunit )

      call write_text( "wp2_sponge_damp_settings%tau_sponge_damp_min = ", &
        wp2_sponge_damp_settings%tau_sponge_damp_min, l_write_to_file, iunit )

      call write_text( "wp2_sponge_damp_settings%tau_sponge_damp_max = ", &
        wp2_sponge_damp_settings%tau_sponge_damp_max, l_write_to_file, iunit )

      call write_text( "wp2_sponge_damp_settings%sponge_damp_depth = ", &
        wp2_sponge_damp_settings%sponge_damp_depth, l_write_to_file, iunit )

      call write_text( "wp3_sponge_damp_settings%tau_sponge_damp_min = ", &
        wp3_sponge_damp_settings%tau_sponge_damp_min, l_write_to_file, iunit )

      call write_text( "wp3_sponge_damp_settings%tau_sponge_damp_max = ", &
        wp3_sponge_damp_settings%tau_sponge_damp_max, l_write_to_file, iunit )

      call write_text( "wp3_sponge_damp_settings%sponge_damp_depth = ", &
        wp3_sponge_damp_settings%sponge_damp_depth, l_write_to_file, iunit )

      call write_text( "up2_vp2_sponge_damp_settings%tau_sponge_damp_min = ", &
        up2_vp2_sponge_damp_settings%tau_sponge_damp_min, l_write_to_file, &
        iunit )

      call write_text( "up2_vp2_sponge_damp_settings%tau_sponge_damp_max = ", &
        up2_vp2_sponge_damp_settings%tau_sponge_damp_max, l_write_to_file, &
        iunit )

      call write_text( "up2_vp2_sponge_damp_settings%sponge_damp_depth = ", &
        up2_vp2_sponge_damp_settings%sponge_damp_depth, l_write_to_file, iunit )

      call write_text( "l_soil_veg = ", l_soil_veg, l_write_to_file, iunit )
      call write_text( "l_restart = ", l_restart, l_write_to_file, iunit )
      call write_text( "l_input_fields = ", l_input_fields, l_write_to_file, iunit )
      call write_text( "restart_path_case = " // restart_path_case, l_write_to_file, iunit )
      call write_text( "time_restart = ", real( time_restart, kind = core_rknd ), &
            l_write_to_file, iunit )
      call write_text( "debug_level = ", debug_level, l_write_to_file, iunit )

      call write_text( "sclr_dim = ", sclr_dim, l_write_to_file, iunit )
      call write_text( "edsclr_dim = ", edsclr_dim, l_write_to_file, iunit )
      call write_text( "iisclr_thl = ", iisclr_thl, l_write_to_file, iunit )
      call write_text( "iisclr_rt = ", iisclr_rt, l_write_to_file, iunit )
      call write_text( "iisclr_CO2 = ", iisclr_CO2, l_write_to_file, iunit )

      call write_text( "sclr_tol = ", sclr_tol(1:sclr_dim), l_write_to_file, iunit )

      ! Pick some default values for stats_setting
      call write_text( "--------------------------------------------------", &
        l_write_to_file, iunit )
      call write_text( "&stats_setting", l_write_to_file, iunit)
      call write_text( "--------------------------------------------------", &
        l_write_to_file, iunit )
      call write_text( "l_stats = ", l_stats, l_write_to_file, iunit )
      call write_text( "fname_prefix = " // fname_prefix, l_write_to_file, iunit )
      call write_text( "stats_fmt = " // stats_fmt, l_write_to_file, iunit)
      call write_text( "stats_tsamp = ", real( stats_tsamp, kind = core_rknd ),&
              l_write_to_file, iunit )
      call write_text( "stats_tout = ", real( stats_tout, kind = core_rknd ), &
              l_write_to_file, iunit )
      call write_text( "l_allow_small_stats_tout = ", l_allow_small_stats_tout, &
              l_write_to_file, iunit)
      call write_text( "Constant flags:", l_write_to_file, iunit )
      call write_text( "l_pos_def = ", l_pos_def, l_write_to_file, iunit )
      call write_text( "l_gamma_Skw = ", l_gamma_Skw, l_write_to_file, iunit)
      call write_text( "l_byteswap_io = ", l_byteswap_io, l_write_to_file, iunit )

      call write_text( "Constant tolerances [units]", l_write_to_file, iunit )
      call write_text( "rt_tol [kg/kg] = ", rt_tol, l_write_to_file, iunit )
      call write_text( "thl_tol [K] = ", thl_tol, l_write_to_file, iunit )
      call write_text( "w_tol [m/s] = ", w_tol, l_write_to_file, iunit )

      if ( l_write_to_file ) close(unit=iunit)

      iunit_grid_adaptation = iunit + 10
      fname_grid_adaptation = '../output/'//trim( runtype )//'_grid_adapt.txt'
      if ( clubb_config_flags%grid_adapt_in_time_method > no_grid_adaptation ) then
        open(unit=iunit_grid_adaptation, file=fname_grid_adaptation, status='replace', action='write')
      end if

    end if ! clubb_at_least_debug_level_api( 1 )

    if ( clubb_at_least_debug_level_api( 1 ) ) then

      ! Print the list of parameters that are being used before the run.
      call write_text( "Parameter          Value", l_write_to_file, iunit, '(4x,A24)')
      call write_text( "---------          -----", l_write_to_file, iunit, '(4x,A24)')
      do j = 1, nparams, 1
        call write_text(params_list(j) // " = ", clubb_params(1,j), &
          l_write_to_file, iunit, '(A31,F27.20)')
      end do

    end if

    ! Allocate stretched grid altitude arrays.
    allocate( momentum_heights(ngrdcol,nzmax), &
              thermodynamic_heights(ngrdcol,nzmax-1) )

    ! Generalized grid test
    ! This block of code is used when a generalized grid test
    ! (ascending vs. descending grid) is being performed.
    if ( l_test_grid_generalization ) then
      allocate( momentum_heights_flip(ngrdcol,nzmax),  & 
                thermodynamic_heights_flip(ngrdcol,nzmax-1) )
    endif

    ! Handle the reading of grid altitudes for
    ! stretched (unevenly-spaced) grid options.
    ! Do some simple error checking for all grid options.
    do i = 1, ngrdcol
      call read_grid_heights( nzmax, grid_type, &                 ! Intent(in)
                              zm_grid_fname, zt_grid_fname, &     ! Intent(in)
                              iunit, &                            ! Intent(in)
                              momentum_heights(i,:), &            ! Intent(out)
                              thermodynamic_heights(i,:), &       ! Intent(out)
                              err_info )                          ! Intent(inout)
      if ( err_info%err_code(i) == clubb_fatal_error ) then
        write(fstderr, *) err_info%err_header(i)
      end if

      ! Generalized grid test
      ! This block of code is used when a generalized grid test
      ! (ascending vs. descending grid) is being performed.
      if ( l_test_grid_generalization ) then
          momentum_heights_flip(i,:) = flip( momentum_heights(i,:), nzmax )
          thermodynamic_heights_flip(i,:) = flip( thermodynamic_heights(i,:), nzmax-1 )
      endif
    end do

    if ( any(err_info%err_code == clubb_fatal_error) ) then
      write(fstderr, *) "Error in read_grid_heights"
      error stop
    end if

    ! These numbers represent host model horizontal grid spacing
    ! which for a single column simulation is effectively infinite
    dummy_dx = 1.0e6_core_rknd ! known magic number
    dummy_dy = 1.0e6_core_rknd ! known magic number

    ! Setup microphysical fields
    call init_microphys( iunit, trim( runtype ), runfile, case_info_file, & ! Intent(in)
                          dummy_dx(1), dummy_dy(1), &                        ! Intent(in)
                          clubb_params(1,:), &                               ! Intent(in)
                          l_diagnose_correlations, &                         ! Intent(in)
                          l_const_Nc_in_cloud, &                             ! Intent(inout)
                          l_fix_w_chi_eta_correlations, &                    ! Intent(inout)
                          hydromet_dim, pdf_dim, hm_metadata, &              ! Intent(out)
                          silhs_config_flags, &                              ! Intent(out)
                          vert_decorr_coef, &                                ! Intent(out)
                          corr_array_n_cloud, corr_array_n_below )           ! Intent(out)

    ! Setup radiation parameters
    call init_radiation( iunit, runfile, case_info_file, & ! Intent(in)
                          l_calc_thlp2_rad )                ! Intent(inout)

    if ( trim( rad_scheme ) == "lba" ) then
      call simple_rad_lba_init( iunit, trim( forcings_file_path ) )
    end if

    ! Initialize CLUBB configurable flags type
    call initialize_clubb_config_flags_type_api( iiPDF_type, & ! Intent(in)
                                                  ipdf_call_placement, & ! Intent(in)
                                                  penta_solve_method, & ! Intent(in)
                                                  tridiag_solve_method, & ! Intent(in)
                                                  saturation_formula, & ! Intent(in)
                                                  grid_remap_method, & ! Intent(in)
                                                  grid_adapt_in_time_method, & ! Intent(in)
                                                  fill_holes_type, & ! Intent(in)
                                                  l_use_precip_frac, & ! Intent(in)
                                                  l_predict_upwp_vpwp, & ! Intent(in)
                                                  l_ho_nontrad_coriolis, & ! Intent(in)
                                                  l_ho_trad_coriolis, & ! Intent(in)
                                                  l_min_wp2_from_corr_wx, & ! Intent(in)
                                                  l_min_xp2_from_corr_wx, & ! Intent(in)
                                                  l_C2_cloud_frac, & ! Intent(in)
                                                  l_diffuse_rtm_and_thlm, & ! Intent(in)
                                                  l_stability_correct_Kh_N2_zm, & ! Intent(in)
                                                  l_calc_thlp2_rad, & ! Intent(in)
                                                  l_upwind_xpyp_ta, & ! Intent(in)
                                                  l_upwind_xm_ma, & ! Intent(in)
                                                  l_uv_nudge, & ! Intent(in)
                                                  l_rtm_nudge, & ! Intent(in)
                                                  l_tke_aniso, & ! Intent(in)
                                                  l_vert_avg_closure, & ! Intent(in)
                                                  l_trapezoidal_rule_zt, & ! Intent(in)
                                                  l_trapezoidal_rule_zm, & ! Intent(in)
                                                  l_call_pdf_closure_twice, & ! Intent(in)
                                                  l_standard_term_ta, & ! Intent(in)
                                                  l_partial_upwind_wp3, & ! Intent(in)
                                                  l_godunov_upwind_wpxp_ta, & ! Intent(in)
                                                  l_godunov_upwind_xpyp_ta, & ! Intent(in)
                                                  l_use_cloud_cover, & ! Intent(in)
                                                  l_diagnose_correlations, & ! Intent(in)
                                                  l_calc_w_corr, & ! Intent(in)
                                                  l_const_Nc_in_cloud, & ! Intent(in)
                                                  l_fix_w_chi_eta_correlations, & ! Intent(in)
                                                  l_stability_correct_tau_zm, & ! Intent(in)
                                                  l_damp_wp2_using_em, & ! Intent(in)
                                                  l_do_expldiff_rtm_thlm, & ! Intent(in)
                                                  l_Lscale_plume_centered, & ! Intent(in)
                                                  l_diag_Lscale_from_tau, & ! Intent(in)
                                                  l_use_C7_Richardson, & ! Intent(in)
                                                  l_use_C11_Richardson, & ! Intent(in)
                                                  l_use_shear_Richardson, & ! Intent(in)
                                                  l_brunt_vaisala_freq_moist, & ! Intent(in)
                                                  l_use_thvm_in_bv_freq, & ! Intent(in)
                                                  l_rcm_supersat_adj, & ! Intent(in)
                                                  l_damp_wp3_Skw_squared, & ! Intent(in)
                                                  l_prescribed_avg_deltaz, & ! Intent(in)
                                                  l_lmm_stepping, & ! Intent(in)
                                                  l_e3sm_config, & ! Intent(in)
                                                  l_vary_convect_depth, & ! Intent(in)
                                                  l_use_tke_in_wp3_pr_turb_term, & ! Intent(in)
                                                  l_use_tke_in_wp2_wp3_K_dfsn, & ! Intent(in)
                                                  l_use_wp3_lim_with_smth_Heaviside, & ! Intent(in)
                                                  l_smooth_Heaviside_tau_wpxp, & ! Intent(in)
                                                  l_modify_limiters_for_cnvg_test, & ! Intent(in)
                                                  l_enable_relaxed_clipping, & ! Intent(in)
                                                  l_linearize_pbl_winds, & ! Intent(in)
                                                  l_mono_flux_lim_thlm, & ! Intent(in)
                                                  l_mono_flux_lim_rtm, & ! Intent(in)
                                                  l_mono_flux_lim_um, & ! Intent(in)
                                                  l_mono_flux_lim_vm, & ! Intent(in)
                                                  l_mono_flux_lim_spikefix, & ! Intent(in)
                                                  l_host_applies_sfc_fluxes, & ! Intent(in)
                                                  l_wp2_fill_holes_tke, & ! Intent(in)
                                                  l_add_dycore_grid, & ! Intent(out)
                                                  clubb_config_flags ) ! Intent(out)

    ! Printing configurable CLUBB flags Inputs
    if ( clubb_at_least_debug_level_api( 1 ) ) then

      if ( l_write_to_file ) then
        open(unit=iunit, file=case_info_file, status='old', action='write', position='append')
      end if

      call write_text( "--------------------------------------------------", &
                        l_write_to_file, iunit )
      call write_text( "&configurable_clubb_flags_nl", l_write_to_file, iunit )
      call write_text( "--------------------------------------------------", &
                        l_write_to_file, iunit )

      call print_clubb_config_flags_api( fstdout, clubb_config_flags ) ! Intent(in)
      call print_clubb_config_flags_api( iunit, clubb_config_flags ) ! Intent(in)

      call write_text( "--------------------------------------------------", &
                        l_write_to_file, iunit )
      call write_text( "git diff src/", l_write_to_file, iunit )
      call write_text( "--------------------------------------------------", &
                        l_write_to_file, iunit )

      if ( l_write_to_file ) then
        write(fstdout, *) "See *setup.txt file in output folder"//NEW_LINE('A')
        close(unit=iunit)
      end if

      if ( l_write_to_file ) then
        call execute_command_line( 'git --no-pager diff >> '//case_info_file )
      end if

    end if ! clubb_at_least_debug_level_api( 1 )

    ! Allocate & initialize variables,
    ! setup grid, setup constants, and setup flags

    ! Setup grid
    call setup_grid_api( nzmax, ngrdcol, sfc_elevation, l_implemented,         & ! intent(in)
                          l_ascending_grid, grid_type, deltaz, zm_init, zm_top, & ! intent(in)
                          momentum_heights, thermodynamic_heights,              & ! intent(in)
                          gr,                                                   & ! intent(out)
                          err_info )                                              ! intent(inout)

    if ( clubb_at_least_debug_level_api( 0 ) ) then
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr, *) "Fatal error calling setup_grid_api in clubb_driver"
        return
      end if
    end if

    ! Generalized grid test
    ! This block of code is used when a generalized grid test
    ! (ascending vs. descending grid) is being performed.
    if ( l_test_grid_generalization ) then
      call setup_grid_api( nzmax, ngrdcol, sfc_elevation, l_implemented,                 & ! in
                            (.not. l_ascending_grid), grid_type, deltaz, zm_init, zm_top, & ! in
                            momentum_heights_flip, thermodynamic_heights_flip,            & ! in
                            gr_desc,                                                      & ! out
                            err_info )                                                      ! inout

      if ( clubb_at_least_debug_level_api( 0 ) ) then
        if ( any(err_info%err_code == clubb_fatal_error) ) then
          write(fstderr, *) err_info%err_header_global
          write(fstderr, *) "Fatal error calling setup_grid_api for gr_desc in clubb_driver"
          return
        end if
      end if

    endif

    ! Allocate and initialize variables

    allocate( um(ngrdcol, gr%nzt) )        ! u wind
    allocate( vm(ngrdcol, gr%nzt) )        ! v wind

    allocate( upwp(ngrdcol, gr%nzm) )      ! vertical u momentum flux
    allocate( vpwp(ngrdcol, gr%nzm) )      ! vertical v momentum flux

    allocate( up2(ngrdcol, gr%nzm) )
    allocate( up3(ngrdcol, gr%nzt) )
    allocate( vp2(ngrdcol, gr%nzm) )
    allocate( vp3(ngrdcol, gr%nzt) )

    allocate( thlm(ngrdcol, gr%nzt) )      ! liquid potential temperature
    allocate( rtm(ngrdcol, gr%nzt) )       ! total water mixing ratio
    allocate( wprtp(ngrdcol, gr%nzm) )     ! w'rt'
    allocate( wpthlp(ngrdcol, gr%nzm) )    ! w'thl'
    allocate( wprcp(ngrdcol, gr%nzm) )     ! w'rc'
    allocate( w_up_in_cloud(ngrdcol, gr%nzt) )
    allocate( w_down_in_cloud(ngrdcol, gr%nzt) )
    allocate( cloudy_updraft_frac(ngrdcol, gr%nzt) )
    allocate( cloudy_downdraft_frac(ngrdcol, gr%nzt) )
    allocate( wp2(ngrdcol, gr%nzm) )       ! w'^2
    allocate( wp3(ngrdcol, gr%nzt) )       ! w'^3
    allocate( wp3_zm(ngrdcol, gr%nzm) )       ! w'^3
    allocate( rtp2(ngrdcol, gr%nzm) )      ! rt'^2
    allocate( thlp2(ngrdcol, gr%nzm) )     ! thl'^2
    allocate( rtpthlp(ngrdcol, gr%nzm) )   ! rt'thlp'

    allocate( p_in_Pa(ngrdcol, gr%nzt) )         ! pressure (pascals)
    allocate( exner(ngrdcol, gr%nzt) )           ! exner function
    allocate( rho(ngrdcol, gr%nzt) )             ! density: t points
    allocate( rho_zm(ngrdcol, gr%nzm) )          ! density: m points
    allocate( rho_ds_zm(ngrdcol, gr%nzm) )       ! dry, static density: m-levs
    allocate( rho_ds_zt(ngrdcol, gr%nzt) )     ! dry, static density: t-levs
    allocate( invrs_rho_ds_zm(ngrdcol, gr%nzm) ) ! inv. dry, static density: m-levs
    allocate( invrs_rho_ds_zt(ngrdcol, gr%nzt) ) ! inv. dry, static density: t-levs
    allocate( thv_ds_zm(ngrdcol, gr%nzm) )       ! dry, base-state theta_v: m-levs
    allocate( thv_ds_zt(ngrdcol, gr%nzt) )       ! dry, base-state theta_v: t-levs

    allocate( thlm_forcing(ngrdcol, gr%nzt) )    ! thlm ls forcing
    allocate( rtm_forcing(ngrdcol, gr%nzt) )     ! rtm ls forcing
    allocate( um_forcing(ngrdcol, gr%nzt) )      ! u forcing
    allocate( vm_forcing(ngrdcol, gr%nzt) )      ! v forcing
    allocate( wprtp_forcing(ngrdcol, gr%nzm) )   ! <w'r_t'> forcing (microphysics)
    allocate( wpthlp_forcing(ngrdcol, gr%nzm) )  ! <w'th_l'> forcing (microphysics)
    allocate( rtp2_forcing(ngrdcol, gr%nzm) )    ! <r_t'^2> forcing (microphysics)
    allocate( thlp2_forcing(ngrdcol, gr%nzm) )   ! <th_l'^2> forcing (microphysics)
    allocate( rtpthlp_forcing(ngrdcol, gr%nzm) ) ! <r_t'th_l'> forcing (microphysics)

    ! Variables used to track perturbed version of winds.
    allocate( um_pert(ngrdcol, gr%nzt) )
    allocate( vm_pert(ngrdcol, gr%nzt) )
    allocate( upwp_pert(ngrdcol, gr%nzm) )
    allocate( vpwp_pert(ngrdcol, gr%nzm) )

    ! Imposed large scale w
    allocate( wm_zm(ngrdcol, gr%nzm) )       ! momentum levels
    allocate( wm_zt(ngrdcol, gr%nzt) )       ! thermodynamic levels

    ! Cloud water variables
    allocate( rcm(ngrdcol, gr%nzt) )
    allocate( delta_zm(ngrdcol, gr%nzm) )
    allocate( cloud_frac(ngrdcol, gr%nzt) )
    allocate( ice_supersat_frac(ngrdcol, gr%nzt) )
    allocate( rcm_in_layer(ngrdcol, gr%nzt) )
    allocate( cloud_cover(ngrdcol, gr%nzt) )
    allocate( invrs_tau_zm(ngrdcol, gr%nzm) )

    ! Passive scalar variables
    ! Note that sclr_dim can be 0
    allocate( wpsclrp_sfc(ngrdcol, sclr_dim) )
    allocate( sclrm(ngrdcol, gr%nzt, sclr_dim) )
    allocate( sclrp2(ngrdcol, gr%nzm, sclr_dim) )
    allocate( sclrp3(ngrdcol, gr%nzt, sclr_dim) )
    allocate( sclrm_forcing(ngrdcol, gr%nzt, sclr_dim) )
    allocate( sclrprtp(ngrdcol, gr%nzm, sclr_dim) )
    allocate( sclrpthlp(ngrdcol, gr%nzm, sclr_dim) )

    allocate( wpedsclrp_sfc(ngrdcol, edsclr_dim) )
    allocate( edsclrm_forcing(ngrdcol, gr%nzt, edsclr_dim) )

    allocate( edsclrm(ngrdcol, gr%nzt, edsclr_dim) )
    allocate( wpsclrp(ngrdcol, gr%nzm, sclr_dim) )

    allocate( sigma_sqd_w(ngrdcol, gr%nzm) )    ! PDF width parameter (momentum levels)
    allocate( sigma_sqd_w_zt(ngrdcol, gr%nzt) ) ! PDF width parameter interp. to t-levs.
    allocate( Skw_zm(ngrdcol, gr%nzm) )         ! Skewness of w on momentum levels
    allocate( Skw_zm_smooth(ngrdcol, gr%nzm) )  ! Skewness of w on momentum levels
    allocate( wp2_zt(ngrdcol, gr%nzt) )         ! wp2 interpolated to thermo. levels
    allocate( ug(ngrdcol, gr%nzt) )             ! u geostrophic wind
    allocate( vg(ngrdcol, gr%nzt) )             ! v geostrophic wind
    allocate( um_ref(ngrdcol, gr%nzt) )         ! Reference u wind for nudging; Michael Falk, 17 Oct 2007
    allocate( vm_ref(ngrdcol, gr%nzt) )         ! Reference v wind for nudging; Michael Falk, 17 Oct 2007
    allocate( thlm_ref(ngrdcol, gr%nzt) )       ! Reference liquid water potential for nudging
    allocate( rtm_ref(ngrdcol, gr%nzt) )        ! Reference total water mixing ratio for nudging
    allocate( thvm(ngrdcol, gr%nzt) )           ! Virtual potential temperature
    allocate( radht(ngrdcol, gr%nzt) )          ! SW + LW heating rate
    allocate( Frad(ngrdcol, gr%nzm) )           ! radiative flux (momentum point)
    allocate( Frad_SW_up(ngrdcol, gr%nzm) )
    allocate( Frad_LW_up(ngrdcol, gr%nzm) )
    allocate( Frad_SW_down(ngrdcol, gr%nzm) )
    allocate( Frad_LW_down(ngrdcol, gr%nzm) )
    allocate( thlprcp(ngrdcol, gr%nzm) )   ! thl'rc'
    allocate( thlp3(ngrdcol, gr%nzt) )     ! thl'^3
    allocate( rtp3(ngrdcol, gr%nzt) )      ! rt'^3

    ! Buoyancy related moments
    allocate( rtpthvp(ngrdcol, gr%nzm) )  ! rt'thv'
    allocate( thlpthvp(ngrdcol, gr%nzm) ) ! thl'thv'
    allocate( wpthvp(ngrdcol, gr%nzm) )   ! w'thv'
    allocate( wp2thvp(ngrdcol, gr%nzt) )  ! w'^2thv'

    allocate( wp2up(ngrdcol, gr%nzt) )   ! w'^2u'
    allocate( wp2rtp(ngrdcol, gr%nzt) )  ! w'^2 rt'
    allocate( wp2thlp(ngrdcol, gr%nzt) ) ! w'^2 thl'
    allocate( uprcp(ngrdcol, gr%nzm) )   ! u'rc'
    allocate( vprcp(ngrdcol, gr%nzm) )   ! v'rc'
    allocate( rc_coef_zm(ngrdcol, gr%nzm) ) ! Coefficient of X'r_c' in Eq. (34)
    allocate( wp4(ngrdcol, gr%nzm) )     ! w'^4
    allocate( wpup2(ngrdcol, gr%nzt) )   ! w'u'^2
    allocate( wpvp2(ngrdcol, gr%nzt) )   ! w'v'^2
    allocate( wp2up2(ngrdcol, gr%nzm) )  ! w'^2 u'^2
    allocate( wp2vp2(ngrdcol, gr%nzm) )  ! w'^2 v'^2

    allocate( Kh_zt(ngrdcol, gr%nzt) )  ! Eddy diffusivity coefficient: thermo. levels
    allocate( Kh_zm(ngrdcol, gr%nzm) )  ! Eddy diffusivity coefficient: momentum levels
    allocate( K_hm(ngrdcol, gr%nzm, hydromet_dim) ) ! Eddy diff. coef. for hydromets.: mom. levs.

    allocate( em(ngrdcol, gr%nzm) )
    allocate( Lscale(ngrdcol, gr%nzt) )

    allocate( tau_zm(ngrdcol, gr%nzm) ) ! Eddy dissipation time scale: momentum levels
    allocate( tau_zt(ngrdcol, gr%nzt) ) ! Eddy dissipation time scale: thermo. levels

    allocate( Nccnm(ngrdcol, gr%nzt) )
    allocate( hydromet(ngrdcol, gr%nzt, hydromet_dim) )    ! All hydrometeor mean fields
    allocate( hydrometp2(ngrdcol, gr%nzm, hydromet_dim) )  ! All < h_m'^2 > fields
    allocate( wphydrometp(ngrdcol, gr%nzm, hydromet_dim) ) ! All < w'h_m' > fields
    allocate( Ncm(ngrdcol, gr%nzt) )   ! Mean cloud droplet concentration, < N_c >
    allocate( wpNcp(ngrdcol, gr%nzm) ) ! < w'N_c' >

    ! High-order passive scalars
    allocate( sclrpthvp(ngrdcol, gr%nzm, sclr_dim) )

    allocate( wp2hmp(ngrdcol, gr%nzt,hydromet_dim), rtphmp_zt(ngrdcol, gr%nzt,hydromet_dim), &
              thlphmp_zt(ngrdcol, gr%nzt,hydromet_dim) )

    ! Allocate rvm_mc, rcm_mc, thlm_mc
    allocate( rvm_mc(ngrdcol, gr%nzt), rcm_mc(ngrdcol, gr%nzt), thlm_mc(ngrdcol, gr%nzt) )

    ! Allocate hydrometeor variables.
    allocate( rrm(ngrdcol, gr%nzt) )
    allocate( Nrm(ngrdcol, gr%nzt) )

    ! Allocate hydromet_pdf_params
    allocate( hydromet_pdf_params(ngrdcol, gr%nzt) )
    

    ! Allocate the correlation arrays
    allocate(corr_array_1_n(ngrdcol, gr%nzt,pdf_dim,pdf_dim))
    allocate(corr_array_2_n(ngrdcol, gr%nzt,pdf_dim,pdf_dim))
    allocate(corr_cholesky_mtx_1(ngrdcol, gr%nzt,pdf_dim,pdf_dim))
    allocate(corr_cholesky_mtx_2(ngrdcol, gr%nzt,pdf_dim,pdf_dim))

    ! Allocate the mean and stddev arrays
    allocate(mu_x_1_n(ngrdcol, gr%nzt,pdf_dim))
    allocate(mu_x_2_n(ngrdcol, gr%nzt,pdf_dim))
    allocate(sigma_x_1_n(ngrdcol, gr%nzt,pdf_dim))
    allocate(sigma_x_2_n(ngrdcol, gr%nzt,pdf_dim))


    ! Allocate microphysics tendencies for <w'rt'>, <w'thl'>,
    ! <rt'^2>, <thl'^2>, and <rt'thl'>.
    allocate( wprtp_mc(ngrdcol, gr%nzm) )
    allocate( wpthlp_mc(ngrdcol, gr%nzm) )
    allocate( rtp2_mc(ngrdcol, gr%nzm) )
    allocate( thlp2_mc(ngrdcol, gr%nzm) )
    allocate( rtpthlp_mc(ngrdcol, gr%nzm) )

    ! Allocate microphysics tendencies for mean hydrometeors and mean cloud
    ! droplet concentration.
    allocate( hydromet_mc(ngrdcol, gr%nzt,hydromet_dim) )
    allocate( Ncm_mc(ngrdcol, gr%nzt) )

    ! Allocate mean sedimentation velocity for hydrometeors
    allocate( hydromet_vel_zt(ngrdcol, gr%nzt,hydromet_dim) )

    ! Allocate covariance sedimentation velocities for <V_hm'hm'>, implicit and
    ! explicit components, such that:
    ! <V_hm'hm'> = <V_hm'hm'>|_impc * <hm> + <V_hm'hm'>|_expc.
    allocate( hydromet_vel_covar_zt_impc(ngrdcol, gr%nzt,hydromet_dim) )
    allocate( hydromet_vel_covar_zt_expc(ngrdcol, gr%nzt,hydromet_dim) )

    allocate( rfrzm(ngrdcol, gr%nzt) )
              
    allocate( X_nl_all_levs(ngrdcol, lh_num_samples, gr%nzt,pdf_dim), &
              X_mixt_comp_all_levs(ngrdcol, lh_num_samples, gr%nzt), &
              lh_rt_clipped(ngrdcol, lh_num_samples, gr%nzt), &
              lh_thl_clipped(ngrdcol, lh_num_samples, gr%nzt), &
              lh_rc_clipped(ngrdcol, lh_num_samples, gr%nzt), &
              lh_rv_clipped(ngrdcol, lh_num_samples, gr%nzt), &
              lh_Nc_clipped(ngrdcol, lh_num_samples, gr%nzt), &
              lh_sample_point_weights(ngrdcol, lh_num_samples, gr%nzt), &
              Nc_in_cloud(ngrdcol,gr%nzt) )

    allocate( wpthlp_sfc(ngrdcol), wprtp_sfc(ngrdcol), upwp_sfc(ngrdcol), vpwp_sfc(ngrdcol) )

    allocate( thlm_init(gr%nzt), &
              rtm_init(gr%nzt), &
              um_init(gr%nzt), &
              vm_init(gr%nzt), &
              ug_init(gr%nzt), &
              vg_init(gr%nzt), &
              wp2_init(gr%nzm), &
              up2_init(gr%nzm), &
              vp2_init(gr%nzm), &
              upwp_init(gr%nzm), &
              rcm_init(gr%nzt), &
              wm_zt_init(gr%nzt), &
              wm_zm_init(gr%nzm), &
              em_init(gr%nzm), &
              exner_init(gr%nzt), &
              thvm_init(gr%nzt), &
              p_in_Pa_init(gr%nzt), &
              rho_init(gr%nzt), &
              rho_zm_init(gr%nzm), &
              rho_ds_zm_init(gr%nzm), &
              rho_ds_zt_init(gr%nzt), &
              invrs_rho_ds_zm_init(gr%nzm), &
              invrs_rho_ds_zt_init(gr%nzt), &
              thv_ds_zm_init(gr%nzm), &
              thv_ds_zt_init(gr%nzt), &
              rtm_ref_init(gr%nzt), &
              thlm_ref_init(gr%nzt), &
              um_ref_init(gr%nzt), &
              vm_ref_init(gr%nzt), &
              Ncm_init(gr%nzt), &
              Nc_in_cloud_init(gr%nzt), &
              Nccnm_init(gr%nzt), &
              sclrm_init(gr%nzt,sclr_dim), &
              edsclrm_init(gr%nzt,edsclr_dim) )

    if ( clubb_config_flags%l_add_dycore_grid ) then
      ! Initialize the dycore grid
      call setup_gr_dycore( iunit, ngrdcol, &                 ! Intent(in)
                            gr%zm(:,1), gr%zm(:,gr%nzm), &    ! Intent(in)
                            l_ascending_grid, &               ! Intent(in)
                            gr_dycore, &                      ! Intent(out)
                            err_info )                        ! Intent(inout)
      gr_output = gr_dycore

    else

      gr_output = gr

    end if

    if (clubb_config_flags%grid_adapt_in_time_method > no_grid_adaptation ) then
      allocate( gr_dens_z(gr%nzm) )
      allocate( gr_dens(gr%nzm) )

      allocate( alt_term(gr%nzm) )
      allocate( lscale_term(gr%nzm) )
      allocate( lscale_term_time_avg(gr%nzm) )
      allocate( chi_term_time_avg(gr%nzm) )
      allocate( richardson_num_term_time_avg(gr%nzm) )
      allocate( norm_min_grid_dens(gr%nzm) )
      allocate( norm_grid_dens(gr%nzm) )
      allocate( chi_term(gr%nzm) )
      allocate( richardson_num_term(gr%nzm) )

      allocate( gr_dens_old_z(ngrdcol,gr%nzm) )
      allocate( gr_dens_old(ngrdcol,gr%nzm) )

      gr_dens_old_z = gr%zm
      gr_dens_old = gr%invrs_dzm

      allocate( cumulative_Lscale(gr%nzm) )
      allocate( cumulative_chi(gr%nzm) )
      allocate( cumulative_richardson_num(gr%nzm) )

      Lscale_counter = 0
      chi_counter = 0
      richardson_num_counter = 0
      do k = 1, gr%nzm
        cumulative_Lscale(k) = zero
        cumulative_chi(k) = zero
        cumulative_richardson_num(k) = zero
      end do

      call setup_gr_fixed_min( iunit, ngrdcol, &               ! Intent(in)
                               gr%zm(:,1), gr%zm(:,gr%nzm), &  ! Intent(in)
                               l_ascending_grid, &             ! Intent(in)
                               gr_fixed_min, err_info )        ! Intent(out)

    end if

    ! This special purpose code only applies to tuner runs where the tune_type
    ! is setup to try all permutations of our model flags
    if ( present( model_flags_array ) ) then
      clubb_config_flags%l_godunov_upwind_wpxp_ta = model_flags_array(1)
      clubb_config_flags%l_godunov_upwind_xpyp_ta = model_flags_array(2)
      clubb_config_flags%l_upwind_xpyp_ta         = model_flags_array(3)
      clubb_config_flags%l_upwind_xm_ma           = model_flags_array(4)
      clubb_config_flags%l_vert_avg_closure       = model_flags_array(5)
      clubb_config_flags%l_standard_term_ta       = model_flags_array(6)
      clubb_config_flags%l_tke_aniso              = model_flags_array(7)
      clubb_config_flags%l_use_cloud_cover        = model_flags_array(8)
      clubb_config_flags%l_rcm_supersat_adj       = model_flags_array(9)

      if ( clubb_config_flags%l_vert_avg_closure ) then
        clubb_config_flags%l_trapezoidal_rule_zt    = .true.
        clubb_config_flags%l_trapezoidal_rule_zm    = .true.
        clubb_config_flags%l_call_pdf_closure_twice = .true.
      end if
    end if

#ifdef _OPENMP
    iunit = omp_get_thread_num( ) + 50 ! Known magic number
#else
    iunit = 50
#endif

    ! Only output radiation files if using a radiation scheme
    if ( trim( rad_scheme ) == "bugsrad" ) then
      stats_metadata%l_output_rad_files = .true.
    end if

#ifdef SILHS
    if ( lh_microphys_type /= lh_microphys_disabled ) then
      l_silhs_out = .true.
    else
      l_silhs_out = .false.
    end if
#else
    l_silhs_out = .false.
#endif

    if ( .not. l_restart ) then
      time_current = time_initial
    else  ! restart
      time_current = time_restart
    end if

    ! Time integration
    ! Call advance_clubb_core once per each statistics output time
    ifinal = floor( ( time_final - time_initial ) / real(dt_main, kind=time_precision) )

    ! Setup filenames and variables to set for setfields, if enabled
    if ( l_input_fields .or. l_restart ) then
      call inputfields_init( iunit, runfile, day, month, year ) ! Intent(in)
    end if

    ! check to make sure dt_rad is a mutliple of dt_main
    if ( abs(dt_rad/dt_main - real(floor(dt_rad/dt_main), kind=core_rknd)) &
          > 1.e-8_core_rknd) then
      error stop "dt_rad must be a multiple of dt_main"
    end if

    if ( l_silhs_rad .and. clubb_config_flags%l_calc_thlp2_rad ) then

      write(fstderr, *) err_info%err_header_global
      write(fstderr, *) "The options l_silhs_rad and l_calc_thlp2_rad are incompatible."
      ! General error -> set all entries to clubb_fatal_error
      err_info%err_code = clubb_fatal_error
      return

    end if

    if ( clubb_config_flags%l_calc_thlp2_rad .and. rad_scheme == "none" ) then
      write(fstderr, *) err_info%err_header_global
      write(fstderr, *) "The options rad_scheme == none and l_calc_thlp2_rad are incompatible."
      ! General error -> set all entries to clubb_fatal_error
      err_info%err_code = clubb_fatal_error
      return
    end if

    ! Currently initialize_clubb does more than just read in the initial sounding.
    ! It also includes other important initializations such as um_ref and vm_ref.
    ! Therefore it should be executed prior to a restart. The restart should overwrite
    ! the initial sounding anyway.
    ! NOTE: initialize_clubb has the capability to be run on the whole 2D arrays we have
    !       but there is no difference in the initializations between columns at the 
    !       moment, so we pass in a value of "1" for ngrdcol and use these _init
    !       arrays that are then copied to the 2D versions. This reduces the time 
    !       spent in this routine and allows us to avoid large CPU to GPU data 
    !       transfers since we need only copy these 1D _init variables in.
    call initialize_clubb( &
          gr, 1, iunit, trim( forcings_file_path ), p_sfc(1), zm_init(1), & ! Intent(in)
          sclr_dim, edsclr_dim, sclr_idx,                                 & ! Intent(in)
          clubb_config_flags,                                             & ! Intent(in)
          l_modify_ic_with_cubic_int,                                     & ! Intent(in)
          clubb_config_flags%l_add_dycore_grid,                           & ! Intent(in)
          clubb_config_flags%grid_adapt_in_time_method,                   & ! Intent(in)
          l_ascending_grid, fcor_y,                                       & ! Intent(in)
          thlm_init, rtm_init, um_init, vm_init, ug_init, vg_init,        & ! Intent(out)
          wp2_init, up2_init, vp2_init, upwp_init, rcm_init,              & ! Intent(out)
          wm_zt_init, wm_zm_init, em_init, exner_init,                    & ! Intent(out)
          thvm_init, p_in_Pa_init,                                        & ! Intent(out)
          rho_init, rho_zm_init, rho_ds_zm_init, rho_ds_zt_init,          & ! Intent(out)
          invrs_rho_ds_zm_init, invrs_rho_ds_zt_init,                     & ! Intent(out)
          thv_ds_zm_init, thv_ds_zt_init,                                 & ! Intent(out)
          rtm_ref_init, thlm_ref_init,                                    & ! Intent(out) 
          um_ref_init, vm_ref_init,                                       & ! Intent(out)
          Ncm_init, Nc_in_cloud_init, Nccnm_init, err_info,               & ! Intent(out)
          deep_soil_T_in_K_init, sfc_soil_T_in_K_init, veg_T_in_K_init,   & ! Intent(out)
          sclrm_init, edsclrm_init )                                        ! Intent(out)

    if ( clubb_at_least_debug_level_api( 0 ) ) then
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr,*) "FATAL ERROR calling initialize_clubb in run_clubb"
        return
      end if
    end if

!$OMP CRITICAL
    if ( stats_metadata%l_output_rad_files ) then
      ! Initialize statistics output, note that this will allocate/initialize stats
      ! variables for all columns, but only create the stats files for the first columns
      if ( clubb_config_flags%grid_adapt_in_time_method == no_grid_adaptation &
            .and. .not. clubb_config_flags%l_add_dycore_grid ) then

        call stats_init_api( iunit, fname_prefix, output_dir, l_stats, & ! In
                              stats_fmt, stats_tsamp, stats_tout, runfile, & ! In
                              hydromet_dim, sclr_dim, edsclr_dim, sclr_tol, & ! In
                              hm_metadata%hydromet_list, hm_metadata%l_mix_rat_hm, & ! In
                              gr%nzm, ngrdcol, nlon, nlat, gr%zt, gr%zm, total_atmos_dim - 1, & ! In
                              complete_alt(1:total_atmos_dim), total_atmos_dim, & ! In
                              complete_momentum(1:total_atmos_dim + 1), day, month, year, & ! In
                              (/lon_vals/), (/lat_vals/), time_current, dt_main, l_silhs_out,&!In
                              T0, ts_nudge, & ! In
                              clubb_params, & ! In
                              clubb_config_flags%l_uv_nudge, & ! In
                              clubb_config_flags%l_tke_aniso, & ! In
                              clubb_config_flags%l_standard_term_ta, & ! In
                              stats_metadata, & ! In/Out
                              stats_zt, stats_zm, stats_sfc, & ! In/Out
                              stats_lh_zt, stats_lh_sfc, & ! In/Out
                              stats_rad_zt, stats_rad_zm, & ! In/Out
                              err_info ) ! In/Out

      else

        call stats_init_w_diff_output_gr( iunit, fname_prefix, output_dir, l_stats, & ! In
                                          stats_fmt, stats_tsamp, stats_tout, runfile, & ! In
                                          hydromet_dim, sclr_dim, edsclr_dim, sclr_tol, & ! In
                                          hm_metadata%hydromet_list, & ! In
                                          hm_metadata%l_mix_rat_hm, & ! In
                                          gr%nzm, ngrdcol, nlon, nlat, gr%zt, gr%zm, & ! In
                                          total_atmos_dim - 1, & ! In
                                          complete_alt(1:total_atmos_dim), total_atmos_dim, & ! In
                                          complete_momentum(1:total_atmos_dim + 1), & ! In
                                          day, month, year, & ! In
                                          (/lon_vals/), (/lat_vals/), time_current, & ! In
                                          dt_main, l_silhs_out,& !In
                                          T0, ts_nudge, &
                                          clubb_params, & ! In
                                          clubb_config_flags%l_uv_nudge, & ! In
                                          clubb_config_flags%l_tke_aniso, & ! In
                                          clubb_config_flags%l_standard_term_ta, & ! In
                                          gr_output%nzm, &
                                          gr_output%zt, &
                                          gr_output%zm, &
                                          stats_metadata, & ! In/Out
                                          stats_zt, stats_zm, stats_sfc, & ! In/Out
                                          stats_lh_zt, stats_lh_sfc, & ! In/Out
                                          stats_rad_zt, stats_rad_zm, & ! In/Out
                                          err_info ) ! In/Out
      endif

    else
      ! Initialize statistics output, note that this will allocate/initialize stats
      ! variables for all columns, but only create the stats files for the first columns
      if ( clubb_config_flags%grid_adapt_in_time_method == no_grid_adaptation &
          .and. .not. clubb_config_flags%l_add_dycore_grid ) then

        call stats_init_api( iunit, fname_prefix, output_dir, l_stats, & ! In
                              stats_fmt, stats_tsamp, stats_tout, runfile, & ! In
                              hydromet_dim, sclr_dim, edsclr_dim, sclr_tol, & ! In
                              hm_metadata%hydromet_list, hm_metadata%l_mix_rat_hm, & ! In
                              gr%nzm, ngrdcol, nlon, nlat, gr%zt, gr%zm, 0, & ! In
                              rad_dummy, 0, rad_dummy, day, month, year, & ! In
                              (/lon_vals/), (/lat_vals/), time_current, dt_main, l_silhs_out,&!In
                              T0, ts_nudge, & ! In
                              clubb_params, & ! In
                              clubb_config_flags%l_uv_nudge, & ! In
                              clubb_config_flags%l_tke_aniso, & ! In
                              clubb_config_flags%l_standard_term_ta, & ! In
                              stats_metadata, & ! In/Out
                              stats_zt, stats_zm, stats_sfc, & ! In/Out
                              stats_lh_zt, stats_lh_sfc, & ! In/Out
                              stats_rad_zt, stats_rad_zm, & ! In/Out
                              err_info ) ! In/Out
      
      else

        call stats_init_w_diff_output_gr( iunit, fname_prefix, output_dir, l_stats, & ! In
                                          stats_fmt, stats_tsamp, stats_tout, runfile, & ! In
                                          hydromet_dim, sclr_dim, edsclr_dim, sclr_tol, & ! In
                                          hm_metadata%hydromet_list, & ! In
                                          hm_metadata%l_mix_rat_hm, & ! In
                                          gr%nzm, ngrdcol, nlon, nlat, gr%zt, gr%zm, 0, & ! In
                                          rad_dummy, 0, rad_dummy, day, month, year, & ! In
                                          (/lon_vals/), (/lat_vals/), time_current, & ! In
                                          dt_main, l_silhs_out,& ! In
                                          T0, ts_nudge, & ! In
                                          clubb_params, & ! In
                                          clubb_config_flags%l_uv_nudge, & ! In
                                          clubb_config_flags%l_tke_aniso, & ! In
                                          clubb_config_flags%l_standard_term_ta, & ! In
                                          gr_output%nzm, &
                                          gr_output%zt, &
                                          gr_output%zm, &
                                          stats_metadata, & ! In/Out
                                          stats_zt, stats_zm, stats_sfc, & ! In/Out
                                          stats_lh_zt, stats_lh_sfc, & ! In/Out
                                          stats_rad_zt, stats_rad_zm, & ! In/Out
                                          err_info ) ! In/Out

      endif

    end if
!$OMP END CRITICAL

    if ( clubb_at_least_debug_level_api( 0 ) ) then
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr,*) "FATAL ERROR in stats_init"
        return
      end if
    end if

#ifdef SILHS
    if ( lh_microphys_type /= lh_microphys_disabled ) then

      ! Setup 2D output of all subcolumns (if enabled)
      call latin_hypercube_2D_output_api( &
              fname_prefix, output_dir, stats_metadata%stats_tout, &
              gr%nzt, pdf_dim, & ! Intent(in)
              gr%zt, time_initial, lh_num_samples, & ! Intent(in)
              nlon, nlat, (/lon_vals/), (/lat_vals/), &
              hm_metadata, &
              T0, ts_nudge, &
              clubb_params(1,:), &
              sclr_dim, sclr_tol, &
              clubb_config_flags%l_uv_nudge, &
              clubb_config_flags%l_tke_aniso, &
              clubb_config_flags%l_standard_term_ta, &  ! Intent(in)
              err_info )    ! Intent(inout)

    end if

    if ( any(err_info%err_code == clubb_fatal_error) ) then
      write(fstderr, *) err_info%err_header_global
      write(fstderr, *) "Fatal error calling latin_hypercube_2D_output_api in run_clubb"
      return
    end if

#endif /* SILHS */

    call setup_radiation_variables( gr%nzm, lin_int_buffer, &
                                    extended_atmos_range_size )

    if( stats_metadata%l_stats ) then

      if ( .not. (( abs(dt_rad/stats_metadata%stats_tout &
                        - real(floor(dt_rad/stats_metadata%stats_tout), kind=core_rknd)) &
            < 1.e-8_core_rknd) .or. &
          ( abs(stats_metadata%stats_tout/dt_rad &
                - real(floor(stats_metadata%stats_tout/dt_rad), kind=core_rknd)) &
            < 1.e-8_core_rknd)) ) then
        error stop &
              "dt_rad must be a multiple of stats_tout or stats_tout must be a mulitple of dt_rad"
      end if

    end if
    
    stats_nsamp = nint( stats_metadata%stats_tsamp / dt_main )
    stats_nout = nint( stats_metadata%stats_tout / dt_main )

    call init_pdf_params_api( gr%nzt, ngrdcol, pdf_params )
    call init_pdf_params_api( gr%nzm, ngrdcol, pdf_params_zm )

    if ( clubb_config_flags%iiPDF_type == iiPDF_new .or. clubb_config_flags%iiPDF_type == iiPDF_new_hybrid ) then
      call init_pdf_implicit_coefs_terms_api( gr%nzt, ngrdcol, sclr_dim, &   ! Intent(in)
                                              pdf_implicit_coefs_terms )     ! Intent(out)
    end if

    if ( .not. trim( microphys_scheme ) == "none" ) then

      call init_precip_fracs_api( gr%nzt, ngrdcol, &
                                  precip_fracs )

    end if

    ! Calculate clubb parameters that should be derrived from other quantities
    call calc_derrived_params_api( gr, ngrdcol, grid_type, deltaz,              & ! Intent(in)
                                   clubb_params,                                & ! Intent(in)
                                   clubb_config_flags%l_prescribed_avg_deltaz,  & ! Intent(in)
                                   nu_vert_res_dep, lmin,                       & ! intent(inout)
                                   mixt_frac_max_mag )                            ! intent(inout)

    ! This checks many things to ensure that all flags and parameters are valid and 
    ! any dependencies between them are handled correctly.
    ! clubb_params might not remain constant throughout the model call, and in theory 
    ! this should be called whenever they change, but in the spirit of performance
    ! we only check these during the initial case setup
    call check_clubb_settings_api( ngrdcol,             & ! Intent(in)
                                   clubb_params,        & ! Intent(in)
                                   l_implemented,       & ! Intent(in)
                                   l_input_fields,      & ! Intent(in)
                                   clubb_config_flags,  & ! intent(in)
                                   err_info )             ! Intent(inout)

    ! Similar call to above, checking parameter values, but because we've placed this in "init_clubb_case"
    ! it will only be checked when initializing the case
    call check_parameters_api( ngrdcol, clubb_params, lmin, & ! Intent(in)
                               err_info )                     ! Intent(inout)

    if ( clubb_at_least_debug_level_api( 0 ) ) then
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr, *) "Fatal error calling check_clubb_settings_api and/or check_parameters_api in clubb_driver"
        return
      end if
    end if

    !$acc enter data copyin( gr, gr%zm, gr%zt, gr%dzm, gr%dzt, gr%invrs_dzt, gr%invrs_dzm, &
    !$acc              gr%weights_zt2zm, gr%weights_zm2zt, &
    !$acc              nu_vert_res_dep, nu_vert_res_dep%nu2, nu_vert_res_dep%nu9, &
    !$acc              nu_vert_res_dep%nu1, nu_vert_res_dep%nu8, nu_vert_res_dep%nu10, &
    !$acc              nu_vert_res_dep%nu6, &
    !$acc              pdf_params, pdf_params_zm, &
    !$acc              sclr_idx, clubb_params, hm_metadata, &
    !$acc              T_sfc, p_sfc, fcor, fcor_y, sfc_elevation, &
    !$acc              dummy_dx, dummy_dy, deltaz, &
    !$acc              pdf_params%w_1, pdf_params%w_2, pdf_params%varnce_w_1, &
    !$acc              pdf_params%varnce_w_2, pdf_params%mixt_frac, &
    !$acc              pdf_params_zm%w_1, pdf_params_zm%w_2, pdf_params_zm%varnce_w_1, &
    !$acc              pdf_params_zm%varnce_w_2, pdf_params_zm%mixt_frac, &
    !$acc              thlm_init, rtm_init, um_init, vm_init, ug_init, vg_init, &
    !$acc              wp2_init, up2_init, vp2_init, upwp_init, rcm_init, wm_zt_init, &
    !$acc              wm_zm_init, em_init, exner_init, thvm_init, p_in_Pa_init, &
    !$acc              rho_init, rho_zm_init, rho_ds_zm_init, rho_ds_zt_init, &
    !$acc              invrs_rho_ds_zm_init, invrs_rho_ds_zt_init, thv_ds_zm_init, &
    !$acc              thv_ds_zt_init, rtm_ref_init, thlm_ref_init, um_ref_init, &
    !$acc              deep_soil_T_in_K_init, sfc_soil_T_in_K_init, veg_T_in_K_init, &
    !$acc              vm_ref_init, rho_ds_zm_dycore_init, err_info, err_info%err_header ) &
    !$acc      create( rtm, wm_zt, ug, vg, um_ref, vm_ref, thlm_forcing, rtm_forcing, &
    !$acc              um_forcing, vm_forcing, wm_zm, wprtp_forcing, wpthlp_forcing, &
    !$acc              rtp2_forcing, thlp2_forcing, rtpthlp_forcing, wpthlp_sfc, &
    !$acc              wprtp_sfc, upwp_sfc, vpwp_sfc, &
    !$acc              deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K, rcm_mc, rvm_mc, &
    !$acc              thlm_mc, radht, wprtp_mc, wpthlp_mc, rtp2_mc, thlp2_mc, rtpthlp_mc, &
    !$acc              thlprcp, wprtp, wpthlp, rtp2, thlp2, rtpthlp, &
    !$acc              wp2, wp3, wp2thvp, rtpthvp, thlpthvp, &
    !$acc              rho_zm, rho, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
    !$acc              invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &
    !$acc              upwp_sfc_pert, vpwp_sfc_pert, rtm_ref, &
    !$acc              thlm_ref, um, upwp, &
    !$acc              vm, vpwp, up2, vp2, up3, vp3, thlm, &
    !$acc              rtp3, thlp3, p_in_Pa, exner, &
    !$acc              rcm, cloud_frac, wpthvp, wp2rtp, wp2up, &
    !$acc              wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, wpup2, wpvp2, &
    !$acc              wp2up2, wp2vp2, ice_supersat_frac, um_pert, vm_pert, upwp_pert, vpwp_pert, &
    !$acc              rcm_in_layer, cloud_cover, wprcp, w_up_in_cloud, w_down_in_cloud, cloudy_updraft_frac, &
    !$acc              cloudy_downdraft_frac, invrs_tau_zm, Kh_zt, Kh_zm, Lscale, &
    !$acc              X_nl_all_levs, x_mixt_comp_all_levs, lh_sample_point_weights, &
    !$acc              lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, lh_rv_clipped, lh_Nc_clipped, &
    !$acc              mu_x_1_n, mu_x_2_n, sigma_x_1_n, sigma_x_2_n, &
    !$acc              corr_array_1_n, corr_array_2_n, corr_cholesky_mtx_1, corr_cholesky_mtx_2, &
    !$acc              rfrzm, thvm, &
    !$acc              err_info%err_code, &
    !$acc              pdf_params%rt_1, pdf_params%rt_2, &
    !$acc              pdf_params%varnce_rt_1, pdf_params%varnce_rt_2, pdf_params%thl_1, &
    !$acc              pdf_params%thl_2, pdf_params%varnce_thl_1, pdf_params%varnce_thl_2, &
    !$acc              pdf_params%corr_w_rt_1, pdf_params%corr_w_rt_2, pdf_params%corr_w_thl_1, &
    !$acc              pdf_params%corr_w_thl_2, pdf_params%corr_rt_thl_1, pdf_params%corr_rt_thl_2, &
    !$acc              pdf_params%alpha_thl, pdf_params%alpha_rt, pdf_params%crt_1, pdf_params%crt_2, &
    !$acc              pdf_params%cthl_1, pdf_params%cthl_2, pdf_params%chi_1, pdf_params%chi_2, &
    !$acc              pdf_params%stdev_chi_1, pdf_params%stdev_chi_2, pdf_params%stdev_eta_1, &
    !$acc              pdf_params%stdev_eta_2, pdf_params%covar_chi_eta_1, pdf_params%covar_chi_eta_2, &
    !$acc              pdf_params%corr_w_chi_1, pdf_params%corr_w_chi_2, pdf_params%corr_w_eta_1, &
    !$acc              pdf_params%corr_w_eta_2, pdf_params%corr_chi_eta_1, pdf_params%corr_chi_eta_2, &
    !$acc              pdf_params%rsatl_1, pdf_params%rsatl_2, pdf_params%rc_1, pdf_params%rc_2, &
    !$acc              pdf_params%cloud_frac_1, pdf_params%cloud_frac_2, &
    !$acc              pdf_params%ice_supersat_frac_1, pdf_params%ice_supersat_frac_2, &
    !$acc              pdf_params_zm%rt_1, pdf_params_zm%rt_2, &
    !$acc              pdf_params_zm%varnce_rt_1, pdf_params_zm%varnce_rt_2, pdf_params_zm%thl_1, &
    !$acc              pdf_params_zm%thl_2, pdf_params_zm%varnce_thl_1, pdf_params_zm%varnce_thl_2, &
    !$acc              pdf_params_zm%corr_w_rt_1, pdf_params_zm%corr_w_rt_2, pdf_params_zm%corr_w_thl_1, &
    !$acc              pdf_params_zm%corr_w_thl_2, pdf_params_zm%corr_rt_thl_1, pdf_params_zm%corr_rt_thl_2, &
    !$acc              pdf_params_zm%alpha_thl, pdf_params_zm%alpha_rt, pdf_params_zm%crt_1, pdf_params_zm%crt_2, &
    !$acc              pdf_params_zm%cthl_1, pdf_params_zm%cthl_2, pdf_params_zm%chi_1, pdf_params_zm%chi_2, &
    !$acc              pdf_params_zm%stdev_chi_1, pdf_params_zm%stdev_chi_2, pdf_params_zm%stdev_eta_1, &
    !$acc              pdf_params_zm%stdev_eta_2, pdf_params_zm%covar_chi_eta_1, pdf_params_zm%covar_chi_eta_2, &
    !$acc              pdf_params_zm%corr_w_chi_1, pdf_params_zm%corr_w_chi_2, pdf_params_zm%corr_w_eta_1, &
    !$acc              pdf_params_zm%corr_w_eta_2, pdf_params_zm%corr_chi_eta_1, pdf_params_zm%corr_chi_eta_2, &
    !$acc              pdf_params_zm%rsatl_1, pdf_params_zm%rsatl_2, pdf_params_zm%rc_1, pdf_params_zm%rc_2, &
    !$acc              pdf_params_zm%cloud_frac_1, pdf_params_zm%cloud_frac_2, &
    !$acc              pdf_params_zm%ice_supersat_frac_1, pdf_params_zm%ice_supersat_frac_2 )

    !$acc enter data if( sclr_dim > 0 ) &
    !$acc      copyin( sclr_tol, sclrm_init ) &
    !$acc      create( sclrm_forcing, wpsclrp_sfc, sclrm, &
    !$acc              wpsclrp, sclrp2, sclrp3, sclrprtp, sclrpthlp, sclrpthvp )

    !$acc enter data if( edsclr_dim > 0 ) &
    !$acc      copyin( edsclrm_init ) &
    !$acc      create( wpedsclrp_sfc, edsclrm_forcing, edsclrm )

    !$acc enter data if( hydromet_dim > 0 ) &
    !$acc      copyin( hm_metadata%l_mix_rat_hm ) &
    !$acc      create( wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt )
    
    call set_case_initial_conditions(ngrdcol, clubb_params, err_info)

  end subroutine init_clubb_case

  !=============================================================================================
  !                                  set_case_initial_conditions
  !=============================================================================================
  subroutine set_case_initial_conditions(ngrdcol, clubb_params, err_info)
  !
  ! Description:
  !   Calling 'init_clubb_case' defines up the settings and initial field values.
  !   This subroutines sets the field variables that get advanced by 'advance_clubb_to_end' to
  !   to those initial values. 
  !
  !---------------------------------------------------------------------------------------------

    use clubb_api_module, only: &
      zero_pdf_implicit_coefs_terms_api, & ! Procedure(s)
      zero_pdf_params_api, &
      calc_derrived_params_api, &
      zero_precip_fracs_api, &
      init_hydromet_pdf_params, & 
      iiPDF_new, &      ! Constants
      iiPDF_new_hybrid

    implicit none

    integer, intent(in) :: &
      ngrdcol

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params  ! Model parameters, C1, nu2, etc.

    !----------------------------------- Input/Output Variables -----------------------------------

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    !----------------------------------- Local Variables -----------------------------------
    integer :: &
      i, b, k,      & ! Loop variables
      sclr, edsclr, &
      sample, hm

    !----------------------------------- Begin Code -----------------------------------

    ! These variables are not used on the GPU and the commented out ones
    ! are not used at all beyond being initialized
    !sigma_sqd_w_zt         = zero        ! PDF width parameter interp. to t-levs.
    !sigma_sqd_w            = zero          ! PDF width parameter (momentum levels)
    !tau_zt                 = zero        ! Eddy dissipation time scale: thermo. levels
    !tau_zm                 = zero          ! Eddy dissipation time scale: momentum levels
    Ncm_mc                 = zero
    Frad                   = zero          ! Radiative flux
    Frad_SW_up             = zero
    Frad_LW_up             = zero
    Frad_SW_down           = zero
    Frad_LW_down           = zero
    wpNcp                  = zero
    do k = 1, gr%nzt
      do i = 1, ngrdcol
        Nccnm(i,k)           = Nccnm_init(k)        ! CCN concentration (COAMPS/MG)
        Ncm(i,k)             = Ncm_init(k)
        Nc_in_cloud(i,k)     = Nc_in_cloud_init(k)
        !em                  = em_init(k)
      end do
    end do
    
    ! Initialize silhs samples to indicate unused status, these are overwritten if silhs is used
    !$acc parallel loop gang vector collapse(3) default(present)
    do sample = 1, lh_num_samples
      do k = 1, gr%nzt
        do i = 1, ngrdcol
          X_nl_all_levs(i,sample,k,:)         = -999._core_rknd
          X_mixt_comp_all_levs(i,sample,k)    = -999
          lh_sample_point_weights(i,sample,k) = -999._core_rknd
        end do
      end do
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, gr%nzt
      do i = 1, ngrdcol
        um(i,k)                     = um_init(k)        ! u wind
        vm(i,k)                     = vm_init(k)        ! v wind
        up3(i,k)                    = zero        ! u'^3
        vp3(i,k)                    = zero        ! v'^3
        thlm(i,k)                   = thlm_init(k)         ! liquid potential temperature
        rtm(i,k)                    = rtm_init(k)        ! total water mixing ratio
        w_up_in_cloud(i,k)          = zero
        w_down_in_cloud(i,k)        = zero
        cloudy_updraft_frac(i,k)    = zero
        cloudy_downdraft_frac(i,k)  = zero
        wp3(i,k)                    = zero        ! w'^3
        p_in_Pa(i,k)                = p_in_Pa_init(k)        ! pressure 
        exner(i,k)                  = exner_init(k)        ! exner
        rho(i,k)                    = rho_init(k)        ! density on thermo. levels
        rho_ds_zt(i,k)              = rho_ds_zt_init(k)         ! dry, static density: t-levs
        invrs_rho_ds_zt(i,k)        = invrs_rho_ds_zt_init(k)        ! inv. dry, static density: t-levs
        thv_ds_zt(i,k)              = thv_ds_zt_init(k)        ! dry, base-state theta_v: t-levs
        thlm_forcing(i,k)           = zero        ! thlm large-scale forcing
        rtm_forcing(i,k)            = zero        ! rtm large-scale forcing
        um_forcing(i,k)             = zero        ! u forcing
        vm_forcing(i,k)             = zero        ! v forcing
        um_pert(i,k)                = zero        ! Variables used to track perturbed version of winds.
        vm_pert(i,k)                = zero
        wm_zt(i,k)                  = wm_zt_init(k)        ! Imposed large scale w - Thermodynamic levels
        rcm(i,k)                    = rcm_init(k)
        cloud_frac(i,k)             = zero
        ice_supersat_frac(i,k)      = zero
        rcm_in_layer(i,k)           = zero
        cloud_cover(i,k)            = zero
        ug(i,k)                     = ug_init(k)        ! u geostrophic wind
        vg(i,k)                     = vg_init(k)        ! v geostrophic wind
        um_ref(i,k)                 = um_ref_init(k)
        vm_ref(i,k)                 = vm_ref_init(k)
        thlm_ref(i,k)               = thlm_ref_init(k)
        rtm_ref(i,k)                = rtm_ref_init(k)
        thvm(i,k)                   = thvm_init(k)         ! Virtual potential temperature
        radht(i,k)                  = zero        ! Heating rate
        thlp3(i,k)                  = zero
        rtp3(i,k)                   = zero
        wp2thvp(i,k)                = zero        ! w'^2thv'
        wp2up(i,k)                  = zero        ! w'^2u'
        wp2rtp(i,k)                 = zero        ! w'^2 rt'
        wp2thlp(i,k)                = zero        ! w'^2 thl'
        wpup2(i,k)                  = zero        ! w'u'^2
        wpvp2(i,k)                  = zero        ! w'v'^2
        Kh_zt(i,k)                  = zero        ! Eddy diffusivity coefficient: thermo. levels
        Lscale(i,k)                 = zero
        rvm_mc(i,k)                 = zero
        rcm_mc(i,k)                 = zero
        thlm_mc(i,k)                = zero
        rfrzm(i,k)                  = zero
      end do
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, gr%nzm
      do i = 1, ngrdcol
        upwp(i,k)             = upwp_init(k)  ! vertical u momentum flux
        vpwp(i,k)             = zero          ! vertical v momentum flux
        up2(i,k)              = up2_init(k)     ! u'^2
        vp2(i,k)              = vp2_init(k)      ! v'^2
        wprtp(i,k)            = zero          ! w'rt'
        wpthlp(i,k)           = zero          ! w'thl'
        wprcp(i,k)            = zero          ! w'rc'
        wp2(i,k)              = wp2_init(k)     ! w'^2
        rtp2(i,k)             = rt_tol**2     ! rt'^2
        thlp2(i,k)            = thl_tol**2    ! thl'^2
        rtpthlp(i,k)          = zero          ! rt'thl'
        rho_zm(i,k)           = rho_zm_init(k)          ! density on moment. levels
        rho_ds_zm(i,k)        = rho_ds_zm_init(k)          ! dry, static density: m-levs
        invrs_rho_ds_zm(i,k)  = invrs_rho_ds_zm_init(k)          ! inv. dry, static density: m-levs
        thv_ds_zm(i,k)        = thv_ds_zm_init(k)          ! dry, base-state theta_v: m-levs
        wprtp_forcing(i,k)    = zero          ! <w'r_t'> forcing 
        wpthlp_forcing(i,k)   = zero          ! <w'th_l'> forcing 
        rtp2_forcing(i,k)     = zero          ! <r_t'^2> forcing 
        thlp2_forcing(i,k)    = zero          ! <th_l'^2> forcing 
        rtpthlp_forcing(i,k)  = zero          ! <r_t'th_l'> forcing 
        upwp_pert(i,k)        = zero          ! Variables used to track perturbed version of winds.
        vpwp_pert(i,k)        = zero
        wm_zm(i,k)            = wm_zm_init(k)          ! Imposed large scale w - Momentum levels
        invrs_tau_zm(i,k)     = zero
        thlprcp(i,k)          = zero
        rtpthvp(i,k)          = zero          ! rt'thv'
        thlpthvp(i,k)         = zero          ! thl'thv'
        wpthvp(i,k)           = zero          ! w'thv'
        uprcp(i,k)            = zero          ! u'rc'
        vprcp(i,k)            = zero          ! v'rc'
        rc_coef_zm(i,k)       = zero          ! Coefficient of X'r_c' in Eq. (34)
        wp4(i,k)              = zero          ! w'^4
        wp2up2(i,k)           = zero          ! w'^2 u'^2
        wp2vp2(i,k)           = zero          ! w'^2 v'^2
        Kh_zm(i,k)            = zero          ! Eddy diffusivity coefficient: momentum levels
        wprtp_mc(i,k)         = zero
        wpthlp_mc(i,k)        = zero
        rtp2_mc(i,k)          = zero
        thlp2_mc(i,k)         = zero
        rtpthlp_mc(i,k)       = zero
      end do
    end do

    ! Surface fluxes
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      wpthlp_sfc(i)       = zero
      wprtp_sfc(i)        = zero
      upwp_sfc(i)         = zero
      vpwp_sfc(i)         = zero
      deep_soil_T_in_K(i) = deep_soil_T_in_K_init(1)
      sfc_soil_T_in_K(i)  = sfc_soil_T_in_K_init(1)
      veg_T_in_K(i)       = veg_T_in_K_init(1)
    end do

    ! Passive scalars
    if ( sclr_dim > 0 ) then

      !$acc parallel loop gang vector collapse(3) default(present)
      do sclr = 1, sclr_dim
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            sclrm(i,k,sclr)         = sclrm_init(k,sclr)
            sclrm_forcing(i,k,sclr) = zero
            sclrp3(i,k,sclr)        = zero
          end do
        end do
      end do

      !$acc parallel loop gang vector collapse(3) default(present)
      do sclr = 1, sclr_dim
        do k = 1, gr%nzm
          do i = 1, ngrdcol
            sclrpthvp(i,k,sclr) = zero
            sclrp2(i,k,sclr)    = sclr_tol(sclr)**2
            sclrprtp(i,k,sclr)  = zero
            sclrpthlp(i,k,sclr) = zero
            wpsclrp(i,k,sclr)   = zero
          end do
        end do
      end do

      !$acc parallel loop gang vector collapse(2) default(present)
      do sclr = 1, sclr_dim
        do i = 1, ngrdcol
          wpsclrp_sfc(i,sclr) = zero
        end do
      end do

    end if

    if ( edsclr_dim > 0 ) then

      !$acc parallel loop gang vector collapse(3) default(present)
      do edsclr = 1, edsclr_dim
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            edsclrm(i,k,edsclr)         = edsclrm_init(k,edsclr)
            edsclrm_forcing(i,k,edsclr) = zero
          end do
        end do
      end do

      !$acc parallel loop gang vector collapse(2) default(present)
      do edsclr = 1, edsclr_dim
        do i = 1, ngrdcol
          wpedsclrp_sfc(i,edsclr)   = zero
        end do
      end do

    end if

    if ( hydromet_dim > 0 ) then

      !$acc parallel loop gang vector collapse(3) default(present)
      do hm = 1, hydromet_dim
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            thlphmp_zt(i,k,hm)  = zero
            rtphmp_zt(i,k,hm)   = zero
            wp2hmp(i,k,hm)      = zero
          end do
        end do
      end do

      !$acc parallel loop gang vector collapse(3) default(present)
      do hm = 1, hydromet_dim
        do k = 1, gr%nzm
          do i = 1, ngrdcol
            wphydrometp(i,k,hm) = zero
          end do
        end do
      end do

      ! These are not used in any GPU code yet
      K_hm            = zero ! Eddy diff. coef. for hydromets.: mom. levs.
      hydromet        = zero
      hydromet_mc     = zero
      hydrometp2      = zero
      hydromet_vel_zt = zero
      hydromet_vel_covar_zt_impc = zero
      hydromet_vel_covar_zt_expc = zero

    end if

    if ( clubb_config_flags%grid_adapt_in_time_method > no_grid_adaptation .and. l_stats ) then
      write(iunit_grid_adaptation, *) 'g', 0, gr%zm

      open( unit=iunit_grid_adaptation+10, file='../output/'//trim( runtype )//'_grid.txt', &
            status='replace', action='write' )
      do b = 1, size(gr%zm(1,:))
          write(iunit_grid_adaptation+10, *) gr%zm(1,b)
      end do
      close(unit=iunit_grid_adaptation+10)
    end if

    ! Zero out our data structures
    call zero_pdf_params_api( pdf_params )
    call zero_pdf_params_api( pdf_params_zm )

    if ( clubb_config_flags%iiPDF_type == iiPDF_new .or. clubb_config_flags%iiPDF_type == iiPDF_new_hybrid ) then
      call zero_pdf_implicit_coefs_terms_api( pdf_implicit_coefs_terms )
    end if

    if ( .not. trim( microphys_scheme ) == "none" ) then

      call zero_precip_fracs_api( precip_fracs )

      ! Currently there is no specific "zero" setting subroutine for hydromet_pdf_params.
      ! It contains only scalars, so no need to allocate anything, but if this ever changes
      ! we want the init_ routine in init_clubb_case and the zero_ call should go here
      do k = 1, gr%nzt
        do i = 1, ngrdcol
          call init_hydromet_pdf_params( hydromet_pdf_params(i,k) )
        end do
      end do

    end if

    if ( .not. l_restart ) then

      time_current = time_initial
      iinit = 1

    else

      time_current = time_restart

      ! Determining what iteration to restart at.
      ! The value is increased by 1 to sychronize with restart data.
      ! Joshua Fasching February 2008

      ! Ensure that iteration num, iinit, is an integer, so that model time is
      !   incremented correctly by iteration number at end of timestep
      if ( abs(mod((time_restart-time_initial),real(dt_main, kind=time_precision))) > &
            real(eps, kind=time_precision) ) then

        write(fstderr,*) "Error: (time_restart-time_initial) ", &
          "is not a multiple of dt_main."
        write(fstderr,*) "time_restart = ", time_restart
        write(fstderr,*) "time_initial = ", time_initial
        write(fstderr,*) "dt_main = ", dt_main
        error stop "Fatal error"

      end if ! mod( (time_restart-time_initial) , dt_main ) /= 0

      iinit = floor( ( time_current - time_initial ) / real(dt_main,kind=time_precision) ) + 1

      do i = 1, ngrdcol
        call restart_clubb &
            ( gr, iunit, hydromet_dim, hm_metadata,                                & ! Intent(in)
              restart_path_case, time_restart,                                     & ! Intent(in)
              um(i,:), upwp(i,:), vm(i,:), vpwp(i,:), up2(i,:), vp2(i,:), rtm(i,:),& ! Intent(inout)
              wprtp(i,:), thlm(i,:), wpthlp(i,:), rtp2(i,:), rtp3(i,:),            & ! Intent(inout)
              thlp2, thlp3, rtpthlp(i,:), wp2, wp3(i,:),                           & ! Intent(inout)
              p_in_Pa(i,:), exner(i,:), rcm(i,:), cloud_frac(i,:),                 & ! Intent(inout)
              wpthvp(i,:), wp2thvp(i,:), wp2up(i,:), rtpthvp(i,:), thlpthvp(i,:),  & ! Intent(inout)
              wp2rtp(i,:), wp2thlp(i,:), uprcp(i,:), vprcp(i,:),                   & ! Intent(inout)
              rc_coef_zm(i,:), wp4(i,:), wpup2(i,:), wpvp2(i,:), wp2up2(i,:),      & ! Intent(inout)
              wp2vp2(i,:), ice_supersat_frac(i,:),                                 & ! Intent(inout)
              wm_zt(i,:), rho(i,:), rho_zm(i,:), rho_ds_zm(i,:),                   & ! Intent(inout)
              rho_ds_zt(i,:), thv_ds_zm(i,:), thv_ds_zt(i,:),                      & ! Intent(inout)
              thlm_forcing(i,:), rtm_forcing(i,:), wprtp_forcing(i,:),             & ! Intent(inout)
              wpthlp_forcing(i,:), rtp2_forcing(i,:),                              & ! Intent(inout)
              thlp2_forcing(i,:), rtpthlp_forcing(i,:),                            & ! Intent(inout)
              hydromet(i,:,:), hydrometp2(i,:,:), wphydrometp(i,:,:),              & ! Intent(inout)
              Ncm(i,:), Nccnm(i,:), thvm(i,:), em(i,:), tau_zm(i,:), tau_zt(i,:),  & ! Intent(inout)
              Kh_zt(i,:), Kh_zm(i,:), ug(i,:), vg(i,:),                            & ! Intent(inout)
              thlprcp(i,:),                                                        & ! Intent(inout)
              sigma_sqd_w(i,:), sigma_sqd_w_zt(i,:), radht(i,:),                   & ! Intent(inout)
              deep_soil_T_in_K(i), sfc_soil_T_in_K(i), veg_T_in_K(i),              & ! Intent(inout)
              pdf_params, pdf_params_zm,                                           & ! Intent(inout)
              rcm_mc(i,:), rvm_mc(i,:), thlm_mc(i,:),                              & ! Intent(out)
              wprtp_mc(i,:), wpthlp_mc(i,:), rtp2_mc(i,:),                         & ! Intent(out)
              thlp2_mc(i,:), rtpthlp_mc(i,:),                                      & ! Intent(out)
              wpthlp_sfc(i), wprtp_sfc(i), upwp_sfc(i), vpwp_sfc(i) )                ! Intent(out)
      end do
  
      ! Calculate invrs_rho_ds_zm and invrs_rho_ds_zt from the values of
      ! rho_ds_zm and rho_ds_zt, respectively, which were read in from the input
      ! file during the call to subroutine restart_clubb.
      invrs_rho_ds_zm = 1.0_core_rknd/rho_ds_zm
      invrs_rho_ds_zt = 1.0_core_rknd/rho_ds_zt

    end if ! ~l_restart

    ! This should be GPUized
    !$acc update host(  deltaz, clubb_params, gr%zt, gr%zm )
    call calc_derrived_params_api( gr, ngrdcol, grid_type, deltaz,              & ! Intent(in)
                                   clubb_params,                                & ! Intent(in)
                                   clubb_config_flags%l_prescribed_avg_deltaz,  & ! Intent(in)
                                   nu_vert_res_dep, lmin,                       & ! intent(inout)
                                   mixt_frac_max_mag )                            ! intent(inout)
    !$acc update device( nu_vert_res_dep%nu1,nu_vert_res_dep%nu2, nu_vert_res_dep%nu6, &
    !$acc                nu_vert_res_dep%nu8, nu_vert_res_dep%nu9, nu_vert_res_dep%nu10 )

  end subroutine set_case_initial_conditions

  !=============================================================================================
  !                                       advance_clubb_to_end
  !=============================================================================================
  subroutine advance_clubb_to_end( ngrdcol, calls_per_out, clubb_params, &
                              l_stdout, l_output_multi_col, l_output_double_prec, &
                              err_info, &
                              l_suppress_stats )
  !
  ! Description:
  !   Calling 'init_clubb_case', sets up the case settings and initial state of all fields.
  !   This subroutine advances those fields from iteration 'iinit' to 'ifinal'. 
  !
  !---------------------------------------------------------------------------------------------

    use radiation_module, only: &
      advance_clubb_radiation, &
      update_radiation_variables, &
      silhs_radiation_driver

    use microphys_driver, only: &
      calc_microphys_scheme_tendcies  !------------------------------------ Procedure(s)

    use stats_type_utilities, only: &
      stat_update_var    !---------------------------------------------------- Procedure

    use stats_clubb_utilities, only: &
      stats_end_timestep_w_diff_output_gr
      
    use prescribe_forcings_module, only: &
      prescribe_forcings

    use Skx_module, only: &
      Skx_func !----------------------------------------- Procedure(s)
  
    use calc_pressure, only: &
      calculate_thvm !-------------------------------- Procedure(s)
  
    use clip_explicit, only: &
      clip_skewness_core !---------------------------- Procedure(s)
    
    use numerical_check, only: &
      invalid_model_arrays !------------------------ Procedure(s)

    use clubb_api_module, only: &
      advance_clubb_core_api, &
      calculate_thlp2_rad_api, &
      stats_begin_timestep_api, &
      stats_end_timestep_api, &
      clubb_at_least_debug_level_api

    use advance_microphys_module, only: &
      advance_microphys  !------------------------------------------------- Procedure(s)

    use soil_vegetation, only: &
      advance_soil_veg    !--------------------------------------------- Procedure(s)
      
    use generalized_grid_test, only: &
      clubb_generalized_grid_testing, & ! Procedure(s)
      silhs_generalized_grid_testing

#ifdef NETCDF
  use output_netcdf, only: &
    output_multi_col_fields
#endif

    use pdf_hydromet_microphys_wrapper, only: &
      pdf_hydromet_microphys_prep    ! Procedure(s)

    use cloud_sed_module, only: &
      cloud_drop_sed  ! Procedure(s)

    use time_dependent_input, only: &
      initialize_t_dependent_input ! Procedure(s)

    use inputfields, only: &
      compute_timestep, stat_fields_reader !------------ Procedure(s)

    use grid_class, only: &
      zt2zm_api, &
      zm2zt_api, &
      zm2zt2zm

    use grid_adaptation_module, only: &
      normalize_grid_density, &
      adapt_grid, &
      calc_grid_dens, &                       ! Procedure(s)
      alt_term_weight, &
      chi_term_weight, &
      richardson_num_term_weight              ! Parameters
      
#ifdef GPTL
    use gptl
#endif

    implicit none

    integer, intent(in) :: &
      ngrdcol, &
      calls_per_out

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params  ! Model parameters, C1, nu2, etc.

    logical, intent(in) ::  & 
      l_stdout,             & ! Whether to print output per timestep
      l_output_multi_col,   & ! Determines whether mutlicolumn data is saved
      l_output_double_prec    ! Flag to enable double precision

    !----------------------------------- Input/Output Variables -----------------------------------
    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header
    
    !----------------------------------- Optional Input Variables -----------------------------------
    logical, optional, intent(in) :: &
      l_suppress_stats

    !----------------------------------- Local Variables -----------------------------------
    logical :: &
      l_stats

    ! coarse-grained timing budget of main time stepping loop
    real( kind = core_rknd ) :: &
      time_loop_init,  &           ! time spent in the beginning part of the main loop [s]
      time_loop_end, &             ! time spent in the end part of the main loop [s]
      time_output_multi_col, &     ! Time spent outputting multi_col data
      time_adapt_grid, &           ! Time spent adapting the grid and remapping all the values
      time_clubb_advance, &        ! time spent in advance_clubb_core [s]
      time_clubb_pdf, &            ! time spent in setup_pdf_parameters
                                    !    and hydrometeor_mixed_moments [s]
      time_SILHS,    &             ! time needed to compute subcolumns [s]q
      time_microphys_scheme, &     ! time needed for calc_microphys_scheme_tendcies [s]
      time_microphys_advance, &    ! time needed for advance_microphys [s]
      time_stop, time_start,  &    ! help variables to measure the time [s]
      time_total                   ! control timer for the overall time spent in the main loop [s]
          
    real( kind = core_rknd ) :: &
      lambda

    integer :: &
      ret_code, &
      itime,          & ! Iteration counters
      itime_nearest,  & ! Used for and inputfields run [s]       
      i, k              ! Local Loop Variables
    !----------------------------------- Begin Code -----------------------------------

    ! If l_suppress_stats is true, turn off stats, otherwise just go by stats_metadata%l_stats
    if ( present(l_suppress_stats) ) then
      l_stats = stats_metadata%l_stats .and. .not. l_suppress_stats
    else
      l_stats = stats_metadata%l_stats
    end if

    !initialize timers    
    time_loop_init = 0.0_core_rknd
    time_clubb_advance = 0.0_core_rknd
    time_clubb_pdf = 0.0_core_rknd
    time_SILHS = 0.0_core_rknd
    time_microphys_advance = 0.0_core_rknd
    time_microphys_scheme = 0.0_core_rknd
    time_loop_end = 0.0_core_rknd
    time_output_multi_col = 0.0_core_rknd
    time_adapt_grid = 0.0_core_rknd
    time_total = 0.0_core_rknd
    time_stop = 0.0_core_rknd
    time_start = 0.0_core_rknd
    
#ifdef GPTL
    ret_code = GPTLsetoption(GPTLprint_method, GPTLfull_tree)
    ret_code = GPTLsetoption(GPTLabort_on_error, 1) ! Abort on GPTL error
    ret_code = GPTLsetoption(GPTLoverhead, 0)       ! Turn off overhead estimate
    ret_code = GPTLinitialize()                     ! Initialize GPTL
#endif

    ! Save time before main loop starts
    call cpu_time( time_start )
    time_total = time_start

    if ( l_test_grid_generalization ) then
      write(fstdout,*) "Performing grid generalization test"
    endif

    !-------------------------------------------------------------------------------
    !                         Main Time Stepping Loop
    !-------------------------------------------------------------------------------

    mainloop: do itime = iinit, ifinal
      
      call cpu_time( time_start ) ! start timer for initial part of main loop
      
      if ( l_stats ) then
        ! When this time step is over, the time will be time + dt_main
        ! We use integer timestep for stats_begin_step
        call stats_begin_timestep_api( itime, stats_nsamp, stats_nout, & ! Intent(in)
                                        stats_metadata )                  ! Intent(inout)
      end if

      if ( l_input_fields ) then

        ! If we're doing an inputfields run, get the values for our
        ! model arrays from a netCDF or GrADS file.
        ! Note:  the time of the 1st LES statistical output is time_initial_LES
        !        plus the time of one LES statistical time step.  For example, a
        !        LES run that starts at 00:00Z and outputs stats every minute has
        !        its first statistical output at 00:01Z.  In order to match the
        !        LES stat output times to CLUBB timesteps, CLUBB's time_current
        !        + dt_main needs to be passed into subroutine compute_timestep.

        l_restart_input = .false.
        call compute_timestep( &
              iunit, stat_files(1), l_restart_input,            & ! Intent(in)
              time_current + real(dt_main,kind=time_precision), & ! Intent(in)
              itime_nearest )                                     ! Intent(out)

        do i = 1, ngrdcol
          call stat_fields_reader( gr, max( itime_nearest, 1 ), hydromet_dim, hm_metadata, & ! In
                                  um(i,:), upwp(i,:), vm(i,:), vpwp(i,:), & ! Inout
                                  up2(i,:), vp2(i,:), rtm(i,:), & ! Inout
                                  wprtp(i,:), thlm(i,:), wpthlp(i,:), & ! Inout
                                  rtp2(i,:), rtp3(i,:), & ! Inout
                                  thlp2(i,:), thlp3(i,:), rtpthlp(i,:), & ! Inout
                                  wp2(i,:), wp3(i,:), & ! Inout
                                  p_in_Pa(i,:), exner(i,:), rcm(i,:), cloud_frac(i,:), & ! Inout
                                  wpthvp(i,:), wp2thvp(i,:), wp2up(i,:), rtpthvp(i,:), & ! Inout
                                  thlpthvp(i,:), & ! Inout
                                  wp2rtp(i,:), wp2thlp(i,:), uprcp(i,:), vprcp(i,:), & ! Inout
                                  rc_coef_zm(i,:), wp4(i,:), wpup2(i,:), & ! Inout
                                  wpvp2(i,:), wp2up2(i,:), & ! Inout
                                  wp2vp2(i,:), ice_supersat_frac(i,:), & ! Inout
                                  wm_zt(i,:), rho(i,:), rho_zm(i,:), rho_ds_zm(i,:), & ! Inout
                                  rho_ds_zt(i,:), thv_ds_zm(i,:), thv_ds_zt(i,:), & ! Inout
                                  thlm_forcing(i,:), rtm_forcing(i,:), wprtp_forcing(i,:), & !""
                                  wpthlp_forcing(i,:), rtp2_forcing(i,:), & ! Inout
                                  thlp2_forcing(i,:), rtpthlp_forcing(i,:), & ! Inout
                                  hydromet(i,:,:), hydrometp2(i,:,:), wphydrometp(i,:,:), & ! Inout
                                  Ncm(i,:), Nccnm(i,:), thvm(i,:), em(i,:), & ! Inout
                                  tau_zm(i,:), tau_zt(i,:), & ! Inout
                                  Kh_zt(i,:), Kh_zm(i,:), ug(i,:), vg(i,:), & ! Inout
                                  thlprcp(i,:), & ! Inout
                                  sigma_sqd_w(i,:), sigma_sqd_w_zt(i,:), radht(i,:), & ! Inout
                                  deep_soil_T_in_K(i), sfc_soil_T_in_K(i), veg_T_in_K(i), & ! Inout
                                  pdf_params, pdf_params_zm ) ! Inout
        end do

        ! clip wp3 if it is input from inputfields
        ! this helps restrict the skewness of wp3_on_wp2
        if( l_input_wp3 ) then
            ! Positive definite quantity
            wp2_zt = max( zm2zt_api( gr%nzm, gr%nzt, ngrdcol, gr, wp2 ), w_tol_sqd )

            call clip_skewness_core( gr%nzt, ngrdcol, gr, sfc_elevation(:), &
                                    clubb_params(:,iSkw_max_mag), wp2_zt, &
                                    clubb_config_flags%l_use_wp3_lim_with_smth_Heaviside, &
                                    wp3 )
        end if
      end if

      if ( clubb_at_least_debug_level_api( 2 ) ) then

        !$acc update host( um, vm, rtm, wprtp, thlm, wpthlp, rtp2, thlp2, rtpthlp, wp2, wp3, &
        !$acc              wp2thvp, rtpthvp, thlpthvp )

        !$acc if( sclr_dim > 0   ) update host( sclrm )
        !$acc if( edsclr_dim > 0 ) update host( edsclrm )

        do i = 1, ngrdcol

          ! Check for NaN values in the model arrays
          if ( invalid_model_arrays( gr%nzm, gr%nzt, hydromet_dim, hm_metadata%hydromet_list, &
                                      sclr_dim, edsclr_dim, &
                                      um(i,:), vm(i,:), rtm(i,:), wprtp(i,:), thlm(i,:), wpthlp(i,:), &
                                      rtp2(i,:), thlp2(i,:), rtpthlp(i,:), wp2(i,:), wp3(i,:), &
                                      wp2thvp(i,:), wp2up(i,:), rtpthvp(i,:), thlpthvp(i,:), &
                                      hydromet(i,:,:), sclrm(i,:,:), edsclrm(i,:,:) ) ) then

            err_info%err_code(i) = clubb_fatal_error
            write(fstderr, *) err_info%err_header(i)
            write(fstderr,*) "Fatal error: a CLUBB variable is NaN in main time stepping loop."
            exit mainloop

          end if

        end do
      end if

      ! Calculate radiation only once in a while
      l_rad_itime = (mod( itime, floor(dt_rad/dt_main) ) == 0 .or. itime == 1)

      ! Calculate thvm for use in prescribe_forcings.
      call calculate_thvm( gr%nzt, ngrdcol, &
                           thlm, rtm, rcm, exner, thv_ds_zt, &
                           thvm )
      
      ! Set large-scale tendencies and subsidence profiles
      call prescribe_forcings( gr, gr%nzm, gr%nzt, ngrdcol, &
                                sclr_dim, edsclr_dim, sclr_idx, & ! In
                                runtype, sfctype, &
                                time_current, time_initial, dt_main, &
                                um, vm, thlm, & ! In
                                p_in_Pa, exner, rho, rho_zm, thvm, & ! In
                                veg_T_in_K, & ! In
                                l_modify_bc_for_cnvg_test, & ! In
                                clubb_config_flags%saturation_formula, & ! In
                                stats_metadata, stats_sfc, & ! In
                                clubb_config_flags%l_add_dycore_grid, & ! In
                                clubb_config_flags%grid_remap_method, & ! In
                                gr%nzm, rho_ds_zm, &
                                gr%zm, &
                                gr_dycore, & ! In
                                rtm, wm_zm, wm_zt, ug, vg, um_ref, vm_ref, & ! Inout
                                thlm_forcing, rtm_forcing, um_forcing, & ! Inout
                                vm_forcing, wprtp_forcing, wpthlp_forcing, & ! Inout
                                rtp2_forcing, thlp2_forcing, rtpthlp_forcing, & ! Inout
                                wpsclrp, sclrm_forcing, edsclrm_forcing, & ! Inout
                                wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, & ! Inout
                                T_sfc, p_sfc, sens_ht, latent_ht, & ! Inout
                                wpsclrp_sfc, wpedsclrp_sfc, err_info ) ! Inout

      if ( clubb_at_least_debug_level_api( 0 ) ) then
        if ( any(err_info%err_code == clubb_fatal_error) ) then
          write(fstderr, *) err_info%err_header_global
          write(fstderr,*) "Fatal error in prescribe_forcings:"
          exit mainloop
        end if
      end if

      !---------------------------------------------------------------
      ! Compute Surface
      !---------------------------------------------------------------
      if ( l_soil_veg ) then

        !$acc update host( rho_zm, wpthlp_sfc, wprtp_sfc, p_sfc, &
        !$acc              deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K )

        call advance_soil_veg( ngrdcol, dt_main, rho_zm(:,1), &
                                Frad_SW_up(:,1), Frad_SW_down(:,1), &
                                Frad_LW_down(:,1), &
                                wpthlp_sfc, wprtp_sfc, p_sfc, &
                                stats_metadata, &
                                stats_sfc, &
                                deep_soil_T_in_K, sfc_soil_T_in_K, &
                                veg_T_in_K )

        !$acc update device( deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K )

      end if

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, gr%nzt
        do i = 1, ngrdcol

          ! Compute total water in ice phase mixing ratio
          rfrzm(i,k) = zero

          ! Add microphysical tendencies to rtm_forcing
          rtm_forcing(i,k) = rtm_forcing(i,k) + rcm_mc(i,k) + rvm_mc(i,k)

          ! Add radiation and microphysical tendencies to thlm_forcing
          thlm_forcing(i,k) = thlm_forcing(i,k) + thlm_mc(i,k) + radht(i,k)

        end do
      end do

      ! Add microphysical tendencies to the forcings for the predictive
      ! variances and covariances.
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, gr%nzm
        do i = 1, ngrdcol
          wprtp_forcing(i,k)   = wprtp_forcing(i,k)   + wprtp_mc(i,k)
          wpthlp_forcing(i,k)  = wpthlp_forcing(i,k)  + wpthlp_mc(i,k)
          rtp2_forcing(i,k)    = rtp2_forcing(i,k)    + rtp2_mc(i,k)
          thlp2_forcing(i,k)   = thlp2_forcing(i,k)   + thlp2_mc(i,k)
          rtpthlp_forcing(i,k) = rtpthlp_forcing(i,k) + rtpthlp_mc(i,k)
        end do
      end do

      if ( hydromet_dim > 0 ) then

        !$acc update host( rfrzm )

        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( hm_metadata%iiri > 0 ) rfrzm(i,k) = rfrzm(i,k) + hydromet(i,k,hm_metadata%iiri)
            if ( hm_metadata%iirs > 0 ) rfrzm(i,k) = rfrzm(i,k) + hydromet(i,k,hm_metadata%iirs)
            if ( hm_metadata%iirg > 0 ) rfrzm(i,k) = rfrzm(i,k) + hydromet(i,k,hm_metadata%iirg)
          end do
        end do

        !$acc update device( rfrzm )

      end if

      ! Add effects of radiation on thlp2
      if ( clubb_config_flags%l_calc_thlp2_rad ) then

        call calculate_thlp2_rad_api( ngrdcol, gr%nzm, gr%nzt, gr,        & ! intent(in)
                                      rcm, thlprcp, radht, clubb_params,  & ! intent(in)
                                      thlp2_forcing )                       ! intent(inout)

      end if
      
      ! Measure time in the beginning part of the main loop
      call cpu_time(time_stop)      
      time_loop_init = time_loop_init + time_stop - time_start      
      call cpu_time(time_start) ! initialize timer for advance_clubb_core

      if ( .not. l_test_grid_generalization ) then

        ! Normal CLUBB run (with ascending grid)
        ! Call the clubb core api for all columns
        call advance_clubb_core_api( &
                gr, gr%nzm, gr%nzt, ngrdcol, &
                l_implemented, dt_main, fcor, fcor_y, sfc_elevation, &!Intent(in)
                hydromet_dim, &                                      ! intent(in)
                sclr_dim, sclr_tol, edsclr_dim, sclr_idx, &          ! intent(in)
                thlm_forcing, rtm_forcing, um_forcing, vm_forcing, & ! Intent(in)
                sclrm_forcing, edsclrm_forcing, wprtp_forcing, &     ! Intent(in)
                wpthlp_forcing, rtp2_forcing, thlp2_forcing, &       ! Intent(in)
                rtpthlp_forcing, wm_zm, wm_zt, &                     ! Intent(in)
                wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, &  ! Intent(in)
                wpsclrp_sfc, wpedsclrp_sfc, &                        ! Intent(in)
                upwp_sfc_pert, vpwp_sfc_pert, &                      ! intent(in)
                rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &         ! Intent(in)
                p_in_Pa, rho_zm, rho, exner, &                       ! Intent(in)
                rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &             ! Intent(in)
                invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &             ! Intent(in)
                hm_metadata%l_mix_rat_hm, &                          ! Intent(in)
                rfrzm, wphydrometp, &                                ! Intent(in)
                wp2hmp, rtphmp_zt, thlphmp_zt, &                     ! Intent(in)
                dummy_dx, dummy_dy, &                                ! Intent(in)
                clubb_params, nu_vert_res_dep, lmin, &               ! Intent(in)
                mixt_frac_max_mag, T0, ts_nudge, &                   ! Intent(in)
                rtm_min, rtm_nudge_max_altitude, &                   ! Intent(in)
                clubb_config_flags, &                                ! Intent(in)
                stats_metadata, &                                    ! Intent(in)
                stats_zt, stats_zm, stats_sfc, &                     ! intent(inout)
                um, vm, upwp, vpwp, up2, vp2, up3, vp3, &            ! Intent(inout)
                thlm, rtm, wprtp, wpthlp, &                          ! Intent(inout)
                wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &       ! Intent(inout)
                sclrm, sclrp2, sclrp3, sclrprtp, sclrpthlp, &        ! Intent(inout)
                wpsclrp, edsclrm, err_info, &                        ! Intent(inout)
                rcm, cloud_frac, &                                   ! Intent(inout)
                wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, &         ! Intent(inout)
                sclrpthvp, &                                         ! Intent(inout)
                wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, &    ! intent(inout)
                wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &   ! intent(inout)
                um_pert, vm_pert, upwp_pert, vpwp_pert, &            ! intent(inout)
                pdf_params, pdf_params_zm, &                         ! Intent(inout)
                pdf_implicit_coefs_terms, &                          ! intent(inout)
                Kh_zm, Kh_zt, &                                      ! intent(out)
                thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, &    ! Intent(out)
                cloudy_updraft_frac, cloudy_downdraft_frac, &        ! Intent(out)
                rcm_in_layer, cloud_cover, invrs_tau_zm, &           ! Intent(out)
                Lscale )                                             ! Intent(out)

      else ! l_test_grid_generalization

        ! Generalized grid test
        ! This block of code is used when a generalized grid test
        ! (ascending vs. descending grid) is being performed.
        call clubb_generalized_grid_testing &
            ( gr, gr_desc, gr%nzm, gr%nzt, ngrdcol, &              ! Intent(in)
              l_implemented, dt_main, fcor, fcor_y, sfc_elevation, &!Intent(in)
              hydromet_dim, &                                      ! intent(in)
              sclr_dim, sclr_tol, edsclr_dim, sclr_idx, &          ! intent(in)
              thlm_forcing, rtm_forcing, um_forcing, vm_forcing, & ! Intent(in)
              sclrm_forcing, edsclrm_forcing, wprtp_forcing, &     ! Intent(in)
              wpthlp_forcing, rtp2_forcing, thlp2_forcing, &       ! Intent(in)
              rtpthlp_forcing, wm_zm, wm_zt, &                     ! Intent(in)
              wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, &  ! Intent(in)
              wpsclrp_sfc, wpedsclrp_sfc,  &                       ! Intent(in)
              upwp_sfc_pert, vpwp_sfc_pert, &                      ! intent(in)
              rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &         ! Intent(in)
              p_in_Pa, rho_zm, rho, exner, &                       ! Intent(in)
              rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &             ! Intent(in)
              invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &             ! Intent(in)
              hm_metadata%l_mix_rat_hm, &                          ! Intent(in)
              rfrzm, wphydrometp, &                                ! Intent(in)
              wp2hmp, rtphmp_zt, thlphmp_zt, &                     ! Intent(in)
              dummy_dx, dummy_dy, &                                ! Intent(in)
              clubb_params, nu_vert_res_dep, lmin, &               ! Intent(in)
              mixt_frac_max_mag, T0, ts_nudge, &                   ! Intent(in)
              rtm_min, rtm_nudge_max_altitude, &                   ! Intent(in)
              clubb_config_flags, &                                ! Intent(in)
              stats_metadata, &                                    ! Intent(in)
              stats_zt, stats_zm, stats_sfc, &                     ! intent(inout)
              um, vm, upwp, vpwp, up2, vp2, up3, vp3, &            ! Intent(inout)
              thlm, rtm, wprtp, wpthlp, &                          ! Intent(inout)
              wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &       ! Intent(inout)
              sclrm, sclrp2, sclrp3, sclrprtp, sclrpthlp, &        ! Intent(inout)
              wpsclrp, edsclrm, err_info, &                        ! Intent(inout)
              rcm, cloud_frac, &                                   ! Intent(inout)
              wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, &         ! Intent(inout)
              sclrpthvp, &                                         ! Intent(inout)
              wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, &    ! intent(inout)
              wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &   ! intent(inout)
              um_pert, vm_pert, upwp_pert, vpwp_pert, &            ! intent(inout)
              pdf_params, pdf_params_zm, &                         ! Intent(inout)
              pdf_implicit_coefs_terms, &                          ! intent(inout)
              Kh_zm, Kh_zt, &                                      ! intent(out)
              thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, &    ! Intent(out)
              cloudy_updraft_frac, cloudy_downdraft_frac, &        ! Intent(out)
              rcm_in_layer, cloud_cover, invrs_tau_zm, &           ! Intent(out)
              Lscale )                                             ! Intent(out)

        if ( clubb_at_least_debug_level_api( 0 ) ) then
          if ( any(err_info%err_code == clubb_generalized_grd_test_err) ) then
            write(fstderr,*) "Mismatch in generalized grid test; Stopping run"
            exit mainloop
          endif
        endif

      endif

      if ( clubb_at_least_debug_level_api( 0 ) ) then
        if ( any(err_info%err_code == clubb_fatal_error) ) then
          write(fstderr, *) err_info%err_header_global
          write(fstderr, *) "Fatal error in advance_clubb_core, ", &
                            "check your parameter values and timestep"
          exit mainloop
        end if
      end if

      ! Measure time in advance_clubb_core
      call cpu_time(time_stop)
      time_clubb_advance = time_clubb_advance + time_stop - time_start
      call cpu_time(time_start) ! initialize timer for setup_pdf_parameters

      if ( .not. l_test_grid_generalization ) then

          ! Call CLUBB's hydrometeor PDF subroutines for cases that use the PDF
          ! for microphysics (either upscaled microphysics or SILHS).
          call pdf_hydromet_microphys_prep &
              ( gr, ngrdcol, pdf_dim, hydromet_dim,              & ! In
                itime, dt_main, vert_decorr_coef,                & ! In
                Nc_in_cloud, cloud_frac, ice_supersat_frac,      & ! In
                rho_ds_zt, Lscale, Kh_zm, hydromet, wphydrometp, & ! In
                corr_array_n_cloud, corr_array_n_below,          & ! In
                hm_metadata, pdf_params, clubb_params,           & ! In
                clubb_config_flags, silhs_config_flags,          & ! In
                l_rad_itime, stats_metadata,                     & ! In
                stats_zt, stats_zm, stats_sfc,                   & ! In/Out
                stats_lh_zt, stats_lh_sfc, err_info,             & ! In/Out
                time_clubb_pdf, time_stop, time_start,           & ! In/Out
                hydrometp2,                                      & ! Out
                mu_x_1_n, mu_x_2_n,                              & ! Out
                sigma_x_1_n, sigma_x_2_n,                        & ! Out
                corr_array_1_n, corr_array_2_n,                  & ! Out
                corr_cholesky_mtx_1, corr_cholesky_mtx_2,        & ! Out
                precip_fracs,                                    & ! In/Out
                rtphmp_zt, thlphmp_zt, wp2hmp,                   & ! Out
                X_nl_all_levs, X_mixt_comp_all_levs,             & ! Out
                lh_sample_point_weights,                         & ! Out
                lh_rt_clipped, lh_thl_clipped, lh_rc_clipped,    & ! Out
                lh_rv_clipped, lh_Nc_clipped,                    & ! Out
                hydromet_pdf_params )                              ! Optional(out)

      else ! l_test_grid_generalization

          call silhs_generalized_grid_testing &
              ( gr, gr_desc, ngrdcol, pdf_dim, hydromet_dim,     & ! In
                itime, dt_main, vert_decorr_coef,                & ! In
                Nc_in_cloud, cloud_frac, ice_supersat_frac,      & ! In
                rho_ds_zt, Lscale, Kh_zm, hydromet, wphydrometp, & ! In
                corr_array_n_cloud, corr_array_n_below,          & ! In
                hm_metadata, pdf_params, clubb_params,           & ! In
                clubb_config_flags, silhs_config_flags,          & ! In
                l_rad_itime, stats_metadata,                     & ! In
                stats_zt, stats_zm, stats_sfc,                   & ! In/Out
                stats_lh_zt, stats_lh_sfc, err_info,             & ! In/Out
                time_clubb_pdf, time_stop, time_start,           & ! In/Out
                hydrometp2,                                      & ! Out
                mu_x_1_n, mu_x_2_n,                              & ! Out
                sigma_x_1_n, sigma_x_2_n,                        & ! Out
                corr_array_1_n, corr_array_2_n,                  & ! Out
                corr_cholesky_mtx_1, corr_cholesky_mtx_2,        & ! Out
                precip_fracs,                                    & ! In/Out
                rtphmp_zt, thlphmp_zt, wp2hmp,                   & ! Out
                X_nl_all_levs, X_mixt_comp_all_levs,             & ! Out
                lh_sample_point_weights,                         & ! Out
                lh_rt_clipped, lh_thl_clipped, lh_rc_clipped,    & ! Out
                lh_rv_clipped, lh_Nc_clipped,                    & ! Out
                hydromet_pdf_params )                              ! Optional(out)

      endif

      ! Error check after pdf_hydromet_microphys_prep
      if ( clubb_at_least_debug_level_api( 0 ) ) then
        if ( any(err_info%err_code == clubb_fatal_error) ) then
          write(fstderr, *) err_info%err_header_global
          write(fstderr,*) "Fatal error after pdf_hydromet_microphys_prep"
          exit mainloop
        end if
      end if

      ! Measure time in SILHS
      call cpu_time(time_stop)
      time_SILHS = time_SILHS + time_stop - time_start
      call cpu_time(time_start) ! initialize timer for calc_microphys_scheme_tendcies

      !----------------------------------------------------------------
      ! Compute Microphysics
      !----------------------------------------------------------------

      if ( .not. trim( microphys_scheme ) == "none" ) then

        !$acc update host( thlm, p_in_Pa, exner, rho, rho_zm, rtm, rcm, cloud_frac, wm_zt, wm_zm, &
        !$acc              wp2, wp3, Kh_zm, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
        !$acc              X_nl_all_levs, X_mixt_comp_all_levs,  lh_sample_point_weights, mu_x_1_n, mu_x_2_n, &
        !$acc              sigma_x_1_n, sigma_x_2_n, corr_array_1_n, corr_array_2_n, &
        !$acc              lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, lh_rv_clipped, lh_Nc_clipped, &
        !$acc              pdf_params%mixt_frac, pdf_params%chi_1, pdf_params%chi_2, &
        !$acc              pdf_params%rt_1, pdf_params%rt_2, pdf_params%rc_1, pdf_params%rc_2, &
        !$acc              pdf_params%thl_1, pdf_params%thl_2, pdf_params%w_1, pdf_params%w_2, &
        !$acc              pdf_params%cloud_frac_1, pdf_params%cloud_frac_2, pdf_params%cthl_1, pdf_params%cthl_2, &
        !$acc              pdf_params%crt_1, pdf_params%crt_2, pdf_params%stdev_chi_1, pdf_params%stdev_chi_2, &
        !$acc              pdf_params%varnce_w_1, pdf_params%varnce_w_2 ) 

        wp3_zm = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, wp3 )

        ! Calculate Skw_zm for use in advance_microphys.
        !$acc data copyin( wp2, wp3_zm, clubb_params ) copyout( Skw_zm )
        call Skx_func( gr%nzm, ngrdcol, wp2, wp3_zm, &
                        w_tol, clubb_params, &
                        Skw_zm )
        !$acc end data
        
        ! This field is smoothed by interpolating to thermodynamic levels and then
        ! interpolating back to momentum levels.
        Skw_zm_smooth = zm2zt2zm( gr%nzm, gr%nzt, ngrdcol, gr, Skw_zm )

        wp2_zt = zm2zt_api( gr%nzm, gr%nzt, ngrdcol, gr, wp2, w_tol_sqd ) ! Positive definite quantity

        ! Call microphysics scheme and produce microphysics tendencies.
        do i = 1, ngrdcol
          call calc_microphys_scheme_tendcies( gr, dt_main, time_current, &               ! In
                                  pdf_dim, hydromet_dim, runtype, &                       ! In
                                  thlm(i,:), p_in_Pa(i,:), exner(i,:), rho(i,:), rho_zm(i,:), rtm(i,:), &! In
                                  rcm(i,:), cloud_frac(i,:), wm_zt(i,:), wm_zm(i,:), wp2_zt(i,:), &      ! In
                                  hydromet(i,:,:), Nc_in_cloud(i,:), &                                ! In
                                  hm_metadata, &                                       ! In
                                  pdf_params, hydromet_pdf_params(i,:), &                 ! In
                                  precip_fracs, &                                         ! In
                                  X_nl_all_levs(i,:,:,:), X_mixt_comp_all_levs(i,:,:), &  ! In
                                  lh_sample_point_weights(i,:,:), &                       ! In
                                  mu_x_1_n(i,:,:), mu_x_2_n(i,:,:), &                     ! In
                                  sigma_x_1_n(i,:,:), sigma_x_2_n(i,:,:), &               ! In
                                  corr_array_1_n(i,:,:,:), corr_array_2_n(i,:,:,:), &     ! In
                                  lh_rt_clipped(i,:,:), lh_thl_clipped(i,:,:), &          ! In
                                  lh_rc_clipped(i,:,:), lh_rv_clipped(i,:,:), &           ! In
                                  lh_Nc_clipped(i,:,:), &                                 ! In
                                  silhs_config_flags%l_lh_importance_sampling, &          ! In
                                  silhs_config_flags%l_lh_instant_var_covar_src, &        ! In
                                  clubb_config_flags%saturation_formula, &                ! In
                                  stats_metadata, &                                       ! In
                                  stats_zt(i), stats_zm(i), &                                   ! Inout
                                  stats_sfc(i), stats_lh_zt(i), &                               ! Inout
                                  Nccnm(i,:), &                                                ! Inout
                                  hydromet_mc(i,:,:), Ncm_mc(i,:), rcm_mc(i,:), rvm_mc(i,:), &                  ! Out
                                  thlm_mc(i,:), hydromet_vel_zt(i,:,:), &                             ! Out
                                  hydromet_vel_covar_zt_impc(i,:,:), &                           ! Out
                                  hydromet_vel_covar_zt_expc(i,:,:), &                           ! Out
                                  wprtp_mc(i,:), wpthlp_mc(i,:), rtp2_mc(i,:), &                         ! Out
                                  thlp2_mc(i,:), rtpthlp_mc(i,:) )                                  ! Out
        end do


        ! Measure time in calc_microphys_scheme_tendcies
        call cpu_time(time_stop)
        time_microphys_scheme = time_microphys_scheme + time_stop - time_start
        call cpu_time(time_start) ! initialize timer for advance_microphys

        ! Advance predictive microphysics fields one model timestep.
        do i = 1, ngrdcol
          call advance_microphys( gr, dt_main, time_current,                                  & ! In
                                  hydromet_dim, hm_metadata,                                  & ! In
                                  wm_zt(i,:), wp2(i,:),                                       & ! In
                                  exner(i,:), rho(i,:), rho_zm(i,:), rcm(i,:),                & ! In
                                  cloud_frac(i,:), Kh_zm(i,:), Skw_zm_smooth(i,:),            & ! In
                                  rho_ds_zm(i,:), rho_ds_zt(i,:), invrs_rho_ds_zt(i,:),       & ! In
                                  hydromet_mc(i,:,:), Ncm_mc(i,:), Lscale(i,:),               & ! In
                                  hydromet_vel_covar_zt_impc(i,:,:),                          & ! In
                                  hydromet_vel_covar_zt_expc(i,:,:),                          & ! In
                                  clubb_params(i,:), nu_vert_res_dep,                         & ! In
                                  clubb_config_flags%tridiag_solve_method,                    & ! In
                                  clubb_config_flags%fill_holes_type,                         & ! In
                                  clubb_config_flags%l_upwind_xm_ma,                          & ! In
                                  stats_metadata,                                             & ! In
                                  stats_zt(i), stats_zm(i), stats_sfc(i),                     & ! Inout
                                  hydromet(i,:,:), hydromet_vel_zt(i,:,:), hydrometp2(i,:,:), & ! Inout
                                  K_hm(i,:,:), Ncm(i,:), Nc_in_cloud(i,:), rvm_mc(i,:),       & ! Inout
                                  thlm_mc(i,:), err_info,                                     & ! Inout
                                  wphydrometp(i,:,:), wpNcp(i,:) )                              ! Out

          if ( clubb_at_least_debug_level_api( 0 ) ) then
            if ( any(err_info%err_code == clubb_fatal_error) ) then
              write(fstderr, *) err_info%err_header(i)
            endif
          endif
        end do

        !$acc update device( rcm_mc, rvm_mc, thlm_mc, wprtp_mc, wpthlp_mc, rtp2_mc, thlp2_mc, rtpthlp_mc )

        !$acc if( hydromet_dim > 0 ) update device( wphydrometp )

        ! Measure time in calc_microphys_scheme_tendcies
        call cpu_time(time_stop)
        time_microphys_advance = time_microphys_advance + time_stop - time_start
        call cpu_time(time_start) ! initialize timer for the end part of the main loop

        if ( clubb_at_least_debug_level_api( 0 ) ) then
          if ( any(err_info%err_code == clubb_fatal_error) ) then
            write(fstderr,*) "Fatal error in advance_microphys:"
            exit mainloop
          endif
        end if

      else

        ! Cloud water sedimentation.
        if ( l_cloud_sed ) then

          !$acc update host( cloud_frac, rcm, rho_zm, rho, exner )
        
          Ncm = Nc_in_cloud * cloud_frac
          rcm_mc = zero
          thlm_mc = zero

          ! Note:  it would be very easy to upscale the cloud water sedimentation
          !        flux, so we should look into adding an upscaled option.
          do i = 1, ngrdcol 
            call cloud_drop_sed( gr, rcm(i,:), Ncm(i,:),                      & ! Intent(in)
                                  rho_zm(i,:), rho(i,:), exner(i,:), sigma_g,  & ! Intent(in)
                                  stats_metadata,                              & ! Intent(in)
                                  stats_zt(i), stats_zm(i),                    & ! Intent(inout)
                                  rcm_mc(i,:), thlm_mc(i,:) )                    ! Intent(inout)
          end do

          !$acc update device( rcm_mc, thlm_mc )

        endif ! l_cloud_sed
      
        if ( stats_metadata%l_stats_samp ) then

          !$acc update host( cloud_frac )

          Ncm = Nc_in_cloud * cloud_frac

          do i = 1, ngrdcol
            call stat_update_var( stats_metadata%iNcm, Ncm(i,:), stats_zt(i) )
            call stat_update_var( stats_metadata%iNc_in_cloud, Nc_in_cloud(i,:), stats_zt(i) )
          end do
          
        endif ! stats_metadata%l_stats_samp

      end if

      if ( stats_metadata%l_stats_samp ) then

        !$acc update host( rvm_mc, rcm_mc, thlm_mc, wprtp_mc, wpthlp_mc, rtp2_mc, thlp2_mc, rtpthlp_mc )

        ! Total microphysical tendency of vapor and cloud water mixing ratios
        do i = 1, ngrdcol
          call stat_update_var( stats_metadata%irvm_mc, rvm_mc(i,:), stats_zt(i) )         ! kg/kg/s
          call stat_update_var( stats_metadata%ircm_mc, rcm_mc(i,:), stats_zt(i) )         ! kg/kg/s
          call stat_update_var( stats_metadata%irtm_mc, rvm_mc(i,:)+rcm_mc(i,:), stats_zt(i) )  ! kg/kg/s
          call stat_update_var( stats_metadata%ithlm_mc, thlm_mc(i,:), stats_zt(i) )       ! K/s
          call stat_update_var( stats_metadata%iwprtp_mc, wprtp_mc(i,:), stats_zm(i) )     ! m*(kg/kg)/s^2
          call stat_update_var( stats_metadata%iwpthlp_mc, wpthlp_mc(i,:), stats_zm(i) )   ! K*m/s^2
          call stat_update_var( stats_metadata%irtp2_mc, rtp2_mc(i,:), stats_zm(i) )       ! (kg/kg)^2/s
          call stat_update_var( stats_metadata%ithlp2_mc, thlp2_mc(i,:), stats_zm(i) )     ! K^2/s
          call stat_update_var( stats_metadata%irtpthlp_mc, rtpthlp_mc(i,:), stats_zm(i) ) ! K*(kg/kg)/s
        end do
      endif

      ! Radiation is always called on the first timestep in order to ensure
      ! that the simulation is subject to radiative heating and cooling from
      ! the first timestep.
      if ( l_rad_itime .and. trim( rad_scheme ) /= "none" ) then

        ! Advance a radiation scheme
        ! With this call ordering, snow and ice water mixing ratio will be
        ! updated by the microphysics, but thlm and rtm will not.  This
        ! somewhat inconsistent, but we would need to move the call to
        ! radiation before the call the microphysics to change this.
        ! -dschanen 17 Aug 2009
        if ( l_silhs_rad ) then

          !$acc update host( rho, rho_zm, p_in_Pa, exner, cloud_frac, ice_supersat_frac, X_nl_all_levs, &
          !$acc              lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, lh_sample_point_weights )

          do i = 1, ngrdcol
            call silhs_radiation_driver( &
                  gr, gr%nzm, gr%nzt, lh_num_samples, pdf_dim, hydromet_dim, hm_metadata, & ! In
                  day, month, year,                                                       & ! In
                  lin_int_buffer,                                                         & ! In
                  extended_atmos_bottom_level,                                            & ! In
                  extended_atmos_top_level,                                               & ! In
                  extended_atmos_range_size,                                              & ! In
                  lat_vals, lon_vals, &
                  time_current, time_initial, rho(i,:), rho_zm(i,:),                      & ! In
                  p_in_Pa(i,:), exner(i,:), cloud_frac(i,:), ice_supersat_frac(i,:),      & ! In
                  X_nl_all_levs(i,:,:,:),lh_rt_clipped(i,:,:), lh_thl_clipped(i,:,:),     & ! In
                  lh_rc_clipped(i,:,:), lh_sample_point_weights(i,:,:), hydromet(i,:,:),  & ! In
                  stats_metadata,                                                         & ! In
                  stats_sfc(i), err_info,                                                 & ! InOut
                  radht(i,:), Frad(i,:), Frad_SW_up(i,:), Frad_LW_up(i,:),                & ! Out
                  Frad_SW_down(i,:), Frad_LW_down(i,:) )                                    ! Out
          end do

          !$acc update device( radht )

        else

          !$acc update host( rho, rho_zm, p_in_Pa, exner, cloud_frac, ice_supersat_frac, thlm, rtm, rcm )

          do i = 1, ngrdcol
            call advance_clubb_radiation( &
                  gr, time_current, time_initial, hydromet_dim,            & ! In
                  day, month, year,                                        & ! In
                  lin_int_buffer,                                          & ! In
                  extended_atmos_bottom_level,                             & ! In
                  extended_atmos_top_level,                                & ! In
                  extended_atmos_range_size,                               & ! In
                  lat_vals, lon_vals, &
                  rho(i,:), rho_zm(i,:), p_in_Pa(i,:),                     & ! In
                  exner(i,:), cloud_frac(i,:), ice_supersat_frac(i,:),     & ! In
                  thlm(i,:), rtm(i,:), rcm(i,:), hydromet(i,:,:),          & ! In
                  hm_metadata, stats_metadata,                             & ! In
                  stats_sfc(i), err_info,                                  & ! InOut
                  radht(i,:), Frad(i,:), Frad_SW_up(i,:), Frad_LW_up(i,:), & ! Out
                  Frad_SW_down(i,:), Frad_LW_down(i,:) )                     ! Out

          end do

          !$acc update device( radht )

        end if ! l_silhs_rad

        if ( clubb_at_least_debug_level_api( 0 ) ) then
          if ( any(err_info%err_code == clubb_fatal_error) ) then
            write(fstderr, *) err_info%err_header_global
            write(fstderr, *) "Fatal error in radiation, " &
                              // "check your parameter values and timestep"
            exit mainloop
          end if
        end if

      end if ! mod( itime, floor(dt_rad/dt_main) ) == 0 .or. itime == 1

      if ( stats_metadata%l_stats_samp ) then

        !$acc update host( radht )

        do i = 1, ngrdcol
          call stat_update_var( stats_metadata%iFrad, Frad(i,:), stats_zm(i) )
        end do

        ! Update the radiation variables here so they are updated every timestep
        do i = 1, ngrdcol
          call update_radiation_variables( gr%nzm, gr%nzt, radht(i,:), Frad_SW_up(i,:), Frad_LW_up(i,:), &
                                            Frad_SW_down(i,:), Frad_LW_down(i,:), &
                                            extended_atmos_range_size, lin_int_buffer, &
                                            stats_metadata, stats_zt(i), stats_zm(i), &
                                            stats_rad_zt(i), stats_rad_zm(i) )
        end do
      end if

      if ( clubb_config_flags%grid_adapt_in_time_method > no_grid_adaptation ) then

        ! interrupt stopping time for time_loop_end if grid gets adapted to save time specific
        ! for calculating the normalized grid density and updating the stats vars
        call cpu_time(time_stop)
        time_loop_end = time_loop_end + time_stop - time_start
        call cpu_time(time_start)

        ! calculate the non-normalized grid density using the refinement criterion
        call calc_grid_dens( ngrdcol, gr, &
                             um, vm, &
                             Lscale, wp2, &
                             pdf_params, &
                             thlm, exner, rtm, &
                             rcm, p_in_Pa, thvm, &
                             ice_supersat_frac, &
                             clubb_params(:,ibv_efold), T0, &
                             clubb_config_flags, &
                             gr_dens_z, gr_dens, &
                             alt_term, lscale_term, &
                             lscale_term_time_avg, &
                             chi_term, &
                             chi_term_time_avg, &
                             richardson_num_term, &
                             richardson_num_term_time_avg, &
                             Lscale_counter, &
                             chi_counter, &
                             richardson_num_counter, &
                             cumulative_Lscale, &
                             cumulative_chi, &
                             cumulative_richardson_num )

        lambda = 0.5

        ! normalize the grid density
        call normalize_grid_density( ngrdcol, &
                                     gr%nzm, &
                                     gr_dens_z, gr_dens, &
                                     lambda, &
                                     gr%nzm, &
                                     gr_fixed_min, &
                                     norm_min_grid_dens, &
                                     norm_grid_dens )

        ! update the stats variables
        if ( stats_metadata%l_stats_samp ) then
          do i = 1, ngrdcol
            call stat_update_var( stats_metadata%igrid_density, gr_dens, stats_zm(i) )
            call stat_update_var( stats_metadata%ialt_term, alt_term_weight*alt_term, stats_zm(i) )
            call stat_update_var( stats_metadata%ilscale_term, lscale_term, stats_zm(i) )
            call stat_update_var( stats_metadata%ilscale_term_time_avg, lscale_term_time_avg, &
                                  stats_zm(i) )
            call stat_update_var( stats_metadata%ichi_term_time_avg, &
                                  chi_term_weight*chi_term_time_avg, stats_zm(i) )
            call stat_update_var( stats_metadata%ibrunt_term_time_avg, &
                                  richardson_num_term_weight*richardson_num_term_time_avg, &
                                  stats_zm(i) )
            call stat_update_var( stats_metadata%inorm_min_grid_dens, norm_min_grid_dens, &
                                  stats_zm(i) )
            call stat_update_var( stats_metadata%inorm_grid_dens, norm_grid_dens, &
                                  stats_zm(i) )
            call stat_update_var( stats_metadata%ichi_term, chi_term, stats_zm(i) )
            call stat_update_var( stats_metadata%ibrunt_term, richardson_num_term, stats_zm(i) )
          end do
        endif

        call cpu_time(time_stop)
        time_adapt_grid = time_adapt_grid + time_stop - time_start
        call cpu_time(time_start)

      endif

      ! End statistics timestep
      ! Only end stats for the first column of values, this writes the values to file,
      ! but since the stats isn't setup to use multiple columns, it will just attempt
      ! write to the same file
      if ( stats_metadata%l_stats_last ) then

        if ( clubb_config_flags%grid_adapt_in_time_method == no_grid_adaptation &
              .and. .not. clubb_config_flags%l_add_dycore_grid ) then
          call stats_end_timestep_api( stats_metadata,                            & ! intent(in)
                                        stats_zt(1), stats_zm(1), stats_sfc(1),    & ! intent(inout)
                                        stats_lh_zt(1), stats_lh_sfc(1),           & ! intent(inout)
                                        stats_rad_zt(1), stats_rad_zm(1), err_info ) ! intent(inout)
        else
          call stats_end_timestep_w_diff_output_gr( stats_metadata,    & ! intent(in)
                                                    gr, gr_output,     & ! intent(in)
                                                    gr%nzm,            & ! intent(in)
                                                    rho_ds_zm(1,:),    & ! intent(in)
                                                    gr%zm(1,:),        & ! intent(in)
                                                    p_sfc(1),          & ! intent(in)
                                                    clubb_config_flags%grid_remap_method, & ! intent(in)
                                                    stats_zt(1),       & ! intent(inout)
                                                    stats_zm(1),       & ! intent(inout)
                                                    stats_sfc(1),      & ! intent(inout)
                                                    stats_lh_zt(1),    & ! intent(inout)
                                                    stats_lh_sfc(1),   & ! intent(inout)
                                                    stats_rad_zt(1),   & ! intent(inout)
                                                    stats_rad_zm(1),   & ! intent(inout)
                                                    err_info           & ! intent(inout)
                                                  )
        endif

        if ( any(err_info%err_code == clubb_fatal_error) ) then
          write(fstderr, *) err_info%err_header_global
          write(fstderr, *) "Fatal error calling stats_end_timestep_api in run_clubb"
          exit mainloop
        end if

      end if

      ! Set Time
      ! Advance time here, not in advance_clubb_core,
      ! in order to facilitate use of stats.
      ! A host model, e.g. WRF, would advance time outside
      ! of advance_clubb_core.  Vince Larson 7 Feb 2006
      time_current = time_initial + real( itime, kind=time_precision ) * &
                                    real(dt_main, kind=time_precision)

      ! This was moved from above to be less confusing to the user,
      ! since before it would appear as though the last timestep
      ! was not executed. -dschanen 19 May 08
      if ( l_stdout ) then
        write(unit=fstdout,fmt='(a,i8,a,f10.1)') 'iteration = ',  &
          itime, '; time = ', time_current
      end if

      ! Measure time in the end part
      call cpu_time(time_stop)
      time_loop_end = time_loop_end + time_stop - time_start
      call cpu_time(time_start)

#ifdef NETCDF
      if ( ngrdcol > 1 .and. l_output_multi_col ) then

        l_last_timestep = itime == ifinal

        call output_multi_col_fields( gr%nzm, gr%nzt, ngrdcol, sclr_dim, edsclr_dim, &
                                      calls_per_out, l_output_double_prec, l_last_timestep, &
                                      gr, dt_main, output_file_prefix, &
                                      day, month, year, time_initial, &
                                      um, vm, up3, vp3, rtm, thlm, rtp3, thlp3, wp3, upwp, vpwp, &
                                      up2, vp2, wprtp, wpthlp, rtp2, thlp2, rtpthlp, wp2, &
                                      sclrm, sclrp3, wpsclrp, sclrp2, sclrprtp, sclrpthlp, &
                                      p_in_Pa, exner, rcm, cloud_frac, wp2thvp, wp2up, wpthvp, &
                                      rtpthvp, &
                                      thlpthvp, sclrpthvp, wp2rtp, wp2thlp, wpup2, wpvp2, &
                                      ice_supersat_frac, uprcp, vprcp, rc_coef_zm, wp4, wp2up2, &
                                      wp2vp2, um_pert, vm_pert, upwp_pert, vpwp_pert, edsclrm, &
                                      rcm_in_layer, cloud_cover, w_up_in_cloud, w_down_in_cloud, &
                                      cloudy_updraft_frac, cloudy_downdraft_frac, wprcp, &
                                      invrs_tau_zm, Kh_zt, Kh_zm, thlprcp )
      end if                      
#endif

      call cpu_time(time_stop)
      time_output_multi_col = time_output_multi_col + time_stop - time_start

      if ( ( clubb_config_flags%grid_adapt_in_time_method > no_grid_adaptation ) &
            .and. (modulo(itime, 120) == 0 .or. l_first_write_to_file) & 
            .and. (stats_metadata%l_stats_last) ) then ! only adapt grid when the average of the last
                                                      ! iterations was just written to file,
                                                      ! in order not to change the grid in between
                                                      ! iterations that get averaged before written
                                                      ! to file

        call cpu_time(time_start)

        l_first_write_to_file = .false.

        call adapt_grid( ngrdcol, gr%nzm, &                                      ! Intent(in)
                         gr%zm, norm_grid_dens, &                                ! Intent(in)
                         gr%zm, norm_min_grid_dens, &                            ! Intent(in)
                         sfc_elevation, l_implemented, &                         ! Intent(in)
                         hydromet_dim, sclr_dim, edsclr_dim, &                   ! Intent(in)
                         thvm, gr%nzt, p_sfc, &                                  ! Intent(in)
                         clubb_config_flags%grid_remap_method, &                 ! Intent(in)
                         gr, &                                                   ! Intent(inout)
                         gr_dens_old_z, gr_dens_old, &                           ! Intent(inout)
                         Lscale_counter, &                                       ! Intent(inout)
                         chi_counter, &                                          ! Intent(inout)
                         richardson_num_counter, &                               ! Intent(inout)
                         cumulative_Lscale, &                                    ! Intent(inout)
                         cumulative_chi, &                                       ! Intent(inout)
                         cumulative_richardson_num, &                            ! Intent(inout)
                         thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &    ! Intent(inout)
                         sclrm_forcing, edsclrm_forcing, wprtp_forcing, &        ! Intent(inout)
                         wpthlp_forcing, rtp2_forcing, thlp2_forcing, &          ! Intent(inout)
                         rtpthlp_forcing, wm_zm, wm_zt, &                        ! Intent(inout)
                         rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &            ! Intent(inout)
                         p_in_Pa, rho_zm, rho, exner, &                          ! Intent(inout)
                         rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &                ! Intent(inout)
                         invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &                ! Intent(inout)
                         rfrzm, wphydrometp, &                                   ! Intent(inout)
                         wp2hmp, rtphmp_zt, thlphmp_zt, &                        ! Intent(inout)
                         um, vm, upwp, vpwp, up2, vp2, up3, vp3, &               ! Intent(inout)
                         thlm, rtm, wprtp, wpthlp, &                             ! Intent(inout)
                         wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &          ! Intent(inout)
                         sclrm, sclrp2, sclrp3, sclrprtp, sclrpthlp, &           ! Intent(inout)
                         wpsclrp, edsclrm, &                                     ! Intent(inout)
                         rcm, cloud_frac, &                                      ! Intent(inout)
                         wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, &            ! Intent(inout)
                         sclrpthvp, &                                            ! Intent(inout)
                         wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, &       ! Intent(inout)
                         wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &      ! Intent(inout)
                         um_pert, vm_pert, upwp_pert, vpwp_pert, &               ! Intent(inout)
                         Kh_zm, Kh_zt, &                                         ! Intent(inout)
                         thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, &       ! Intent(inout)
                         cloudy_updraft_frac, cloudy_downdraft_frac, &           ! Intent(inout)
                         rcm_in_layer, cloud_cover, invrs_tau_zm, &              ! Intent(inout)
                         Lscale, err_info )                                      ! Intent(inout)

        if ( any(err_info%err_code == clubb_fatal_error) ) then
          write(fstderr, *) err_info%err_header_global
          write(fstderr, *) "Fatal error calling adapt_grid in run_clubb"
          exit mainloop
        end if

        if ( .not. clubb_config_flags%l_add_dycore_grid .and. l_t_dependent ) then
          ! if we want to adapt the grid over time, but don't want to simulate the host model by
          ! getting forcings on the dycore grid, we need to adjust the forcings everytime for the
          ! newly created grid from the raw data in the forcings file
          ! TODO does this only work with 1D input?
          call initialize_t_dependent_input &
                        ( iunit, runtype, gr%nzt, gr%zt(1,:), p_in_Pa, &
                          clubb_config_flags%l_add_dycore_grid, &
                          clubb_config_flags%grid_adapt_in_time_method )
        end if

        call cpu_time(time_stop)
        time_adapt_grid = time_adapt_grid + time_stop - time_start
      end if

      if ( clubb_config_flags%grid_adapt_in_time_method > no_grid_adaptation .and. l_stats ) then
        call cpu_time(time_start)
        write(iunit_grid_adaptation, *) 'g', itime, gr%zm
        call cpu_time(time_stop)
        time_adapt_grid = time_adapt_grid + time_stop - time_start
      end if ! grid_adapt_in_time_method > 0 .and. modulo(itime, 1) == 0

      ! Checking error status for stats_end_timestep_api
      if ( clubb_at_least_debug_level_api( 0 ) ) then
        if ( any(err_info%err_code == clubb_fatal_error) ) then
          write(fstderr, *) err_info%err_header_global
          write(fstderr, *) "Fatal error calling stats_end_timestep_api in run_clubb"
          exit mainloop
        endif
      end if

    end do mainloop ! itime=1, ifinal
    

    !-------------------------------------------------------------------------------
    !                       End Main Time Stepping Loop
    !-------------------------------------------------------------------------------

    ! Measure overall time in the main loop
    call cpu_time(time_stop)
    time_total = time_stop - time_total ! subtract previously saved start time 

    ! Print timers 
    write(unit=fstdout, fmt='(a,f10.4)') 'CLUBB-TIMER time_loop_init =         ', time_loop_init
    write(unit=fstdout, fmt='(a,f10.4)') 'CLUBB-TIMER time_clubb_advance =     ', time_clubb_advance
    write(unit=fstdout, fmt='(a,f10.4)') 'CLUBB-TIMER time_clubb_pdf =         ', time_clubb_pdf
    write(unit=fstdout, fmt='(a,f10.4)') 'CLUBB-TIMER time_SILHS =             ', time_SILHS
    write(unit=fstdout, fmt='(a,f10.4)') 'CLUBB-TIMER time_microphys_scheme =  ', time_microphys_scheme
    write(unit=fstdout, fmt='(a,f10.4)') 'CLUBB-TIMER time_microphys_advance = ', time_microphys_advance
    write(unit=fstdout, fmt='(a,f10.4)') 'CLUBB-TIMER time_loop_end =          ', time_loop_end
    write(unit=fstdout, fmt='(a,f10.4)') 'CLUBB-TIMER time_output_multi_col =  ', time_output_multi_col
    write(unit=fstdout, fmt='(a,f10.4)') 'CLUBB-TIMER time_adapt_grid =        ', time_adapt_grid
    write(unit=fstdout, fmt='(a,f10.4)') 'CLUBB-TIMER time_total =             ', time_total

#ifdef GPTL
    ret_code = GPTLpr(1)
    ret_code = GPTLfinalize()
#endif

  end subroutine advance_clubb_to_end

  !=============================================================================================
  !                                       clean_up_clubb
  !=============================================================================================
  subroutine clean_up_clubb(ngrdcol, clubb_params, err_info)
  !
  ! Description:
  !   Deallocate everything that was allocated in 'init_clubb_case'. After running this
  !   it should be possible to rerun 'init_clubb_case' without errors.
  !
  !---------------------------------------------------------------------------------------------
  
    
    use grid_class, only: &
      grid

    use inputfields, only: &
      cleanup_input_fields

    use sponge_layer_damping, only: &
      thlm_sponge_damp_settings,    & !---------------- Variable(s)
      rtm_sponge_damp_settings,     &
      uv_sponge_damp_settings,      &
      wp2_sponge_damp_settings,     &
      wp3_sponge_damp_settings,     &
      up2_vp2_sponge_damp_settings, &
      thlm_sponge_damp_profile,     &
      rtm_sponge_damp_profile,      &
      uv_sponge_damp_profile,       &
      wp2_sponge_damp_profile,      &
      wp3_sponge_damp_profile,      &
      up2_vp2_sponge_damp_profile

    use clubb_api_module, only: &
      finalize_tau_sponge_damp_api, &    !---------------- Procedure(s)
      cleanup_grid_api

    use time_dependent_input, only: &
      l_t_dependent,    & !---------------------------- Variable(s)
      finalize_t_dependent_input !--------------------- Procedure(s)

    use extended_atmosphere_module, only: &
      finalize_extended_atm

    use variables_radiation_module, only: &
      cleanup_radiation_variables

    use microphys_init_cleanup, only: &
      cleanup_microphys

    use stats_clubb_utilities, only: &
      stats_finalize_api

#ifdef SILHS
    use parameters_microphys, only: &
        lh_microphys_type,     & !------------------------ Variable(s)
        lh_microphys_disabled

    use silhs_api_module, only: &
        latin_hypercube_2D_close_api

    use latin_hypercube_arrays, only: &
        cleanup_latin_hypercube_arrays !------------------ Procedure(s)
#endif

    implicit none

    integer, intent(in) :: &
      ngrdcol

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params  ! Model parameters, C1, nu2, etc.

    !----------------------------------- Input/Output Variables -----------------------------------
    type(err_info_type), intent(in) :: &
      err_info        ! err_info struct containing err_code and err_header

    !--------------------------------- Begin Code ---------------------------------
      
    !$acc exit data delete( gr, gr%zm, gr%zt, gr%dzm, gr%dzt, gr%invrs_dzt, gr%invrs_dzm, &
    !$acc              gr%weights_zt2zm, gr%weights_zm2zt, &
    !$acc              nu_vert_res_dep, nu_vert_res_dep%nu2, nu_vert_res_dep%nu9, &
    !$acc              nu_vert_res_dep%nu1, nu_vert_res_dep%nu8, nu_vert_res_dep%nu10, &
    !$acc              nu_vert_res_dep%nu6, &
    !$acc              pdf_params, pdf_params_zm, &
    !$acc              sclr_idx, clubb_params, hm_metadata, &
    !$acc              T_sfc, p_sfc, fcor, fcor_y, sfc_elevation, &
    !$acc              dummy_dx, dummy_dy, &
    !$acc              pdf_params%w_1, pdf_params%w_2, pdf_params%varnce_w_1, &
    !$acc              pdf_params%varnce_w_2, pdf_params%mixt_frac, &
    !$acc              pdf_params_zm%w_1, pdf_params_zm%w_2, pdf_params_zm%varnce_w_1, &
    !$acc              pdf_params_zm%varnce_w_2, pdf_params_zm%mixt_frac, &
    !$acc              thlm_init, rtm_init, um_init, vm_init, ug_init, vg_init, &
    !$acc              wp2_init, up2_init, vp2_init, upwp_init, rcm_init, wm_zt_init, &
    !$acc              wm_zm_init, em_init, exner_init, thvm_init, p_in_Pa_init, &
    !$acc              rho_init, rho_zm_init, rho_ds_zm_init, rho_ds_zt_init, &
    !$acc              invrs_rho_ds_zm_init, invrs_rho_ds_zt_init, thv_ds_zm_init, &
    !$acc              thv_ds_zt_init, rtm_ref_init, thlm_ref_init, um_ref_init, &
    !$acc              deep_soil_T_in_K_init, sfc_soil_T_in_K_init, veg_T_in_K_init, &
    !$acc              vm_ref_init, rho_ds_zm_dycore_init, err_info, err_info%err_header, &
    !$acc              rtm, wm_zt, ug, vg, um_ref, vm_ref, thlm_forcing, rtm_forcing, &
    !$acc              um_forcing, vm_forcing, wm_zm, wprtp_forcing, wpthlp_forcing, &
    !$acc              rtp2_forcing, thlp2_forcing, rtpthlp_forcing, wpthlp_sfc, &
    !$acc              wprtp_sfc, upwp_sfc, vpwp_sfc, &
    !$acc              deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K, rcm_mc, rvm_mc, &
    !$acc              thlm_mc, radht, wprtp_mc, wpthlp_mc, rtp2_mc, thlp2_mc, rtpthlp_mc, &
    !$acc              thlprcp, wprtp, wpthlp, rtp2, thlp2, rtpthlp, &
    !$acc              wp2, wp3, wp2thvp, rtpthvp, thlpthvp, &
    !$acc              rho_zm, rho, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
    !$acc              invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &
    !$acc              upwp_sfc_pert, vpwp_sfc_pert, rtm_ref, &
    !$acc              thlm_ref, um, upwp, &
    !$acc              vm, vpwp, up2, vp2, up3, vp3, thlm, &
    !$acc              rtp3, thlp3, p_in_Pa, exner, &
    !$acc              rcm, cloud_frac, wpthvp, wp2rtp, wp2up, &
    !$acc              wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, wpup2, wpvp2, &
    !$acc              wp2up2, wp2vp2, ice_supersat_frac, um_pert, vm_pert, upwp_pert, vpwp_pert, &
    !$acc              rcm_in_layer, cloud_cover, wprcp, w_up_in_cloud, w_down_in_cloud, cloudy_updraft_frac, &
    !$acc              cloudy_downdraft_frac, invrs_tau_zm, Kh_zt, Kh_zm, Lscale, &
    !$acc              X_nl_all_levs, x_mixt_comp_all_levs, lh_sample_point_weights, &
    !$acc              lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, lh_rv_clipped, lh_Nc_clipped, &
    !$acc              mu_x_1_n, mu_x_2_n, sigma_x_1_n, sigma_x_2_n, &
    !$acc              corr_array_1_n, corr_array_2_n, corr_cholesky_mtx_1, corr_cholesky_mtx_2, &
    !$acc              rfrzm, thvm, deltaz, &
    !$acc              err_info%err_code, &
    !$acc              pdf_params%rt_1, pdf_params%rt_2, &
    !$acc              pdf_params%varnce_rt_1, pdf_params%varnce_rt_2, pdf_params%thl_1, &
    !$acc              pdf_params%thl_2, pdf_params%varnce_thl_1, pdf_params%varnce_thl_2, &
    !$acc              pdf_params%corr_w_rt_1, pdf_params%corr_w_rt_2, pdf_params%corr_w_thl_1, &
    !$acc              pdf_params%corr_w_thl_2, pdf_params%corr_rt_thl_1, pdf_params%corr_rt_thl_2, &
    !$acc              pdf_params%alpha_thl, pdf_params%alpha_rt, pdf_params%crt_1, pdf_params%crt_2, &
    !$acc              pdf_params%cthl_1, pdf_params%cthl_2, pdf_params%chi_1, pdf_params%chi_2, &
    !$acc              pdf_params%stdev_chi_1, pdf_params%stdev_chi_2, pdf_params%stdev_eta_1, &
    !$acc              pdf_params%stdev_eta_2, pdf_params%covar_chi_eta_1, pdf_params%covar_chi_eta_2, &
    !$acc              pdf_params%corr_w_chi_1, pdf_params%corr_w_chi_2, pdf_params%corr_w_eta_1, &
    !$acc              pdf_params%corr_w_eta_2, pdf_params%corr_chi_eta_1, pdf_params%corr_chi_eta_2, &
    !$acc              pdf_params%rsatl_1, pdf_params%rsatl_2, pdf_params%rc_1, pdf_params%rc_2, &
    !$acc              pdf_params%cloud_frac_1, pdf_params%cloud_frac_2, &
    !$acc              pdf_params%ice_supersat_frac_1, pdf_params%ice_supersat_frac_2, &
    !$acc              pdf_params_zm%rt_1, pdf_params_zm%rt_2, &
    !$acc              pdf_params_zm%varnce_rt_1, pdf_params_zm%varnce_rt_2, pdf_params_zm%thl_1, &
    !$acc              pdf_params_zm%thl_2, pdf_params_zm%varnce_thl_1, pdf_params_zm%varnce_thl_2, &
    !$acc              pdf_params_zm%corr_w_rt_1, pdf_params_zm%corr_w_rt_2, pdf_params_zm%corr_w_thl_1, &
    !$acc              pdf_params_zm%corr_w_thl_2, pdf_params_zm%corr_rt_thl_1, pdf_params_zm%corr_rt_thl_2, &
    !$acc              pdf_params_zm%alpha_thl, pdf_params_zm%alpha_rt, pdf_params_zm%crt_1, pdf_params_zm%crt_2, &
    !$acc              pdf_params_zm%cthl_1, pdf_params_zm%cthl_2, pdf_params_zm%chi_1, pdf_params_zm%chi_2, &
    !$acc              pdf_params_zm%stdev_chi_1, pdf_params_zm%stdev_chi_2, pdf_params_zm%stdev_eta_1, &
    !$acc              pdf_params_zm%stdev_eta_2, pdf_params_zm%covar_chi_eta_1, pdf_params_zm%covar_chi_eta_2, &
    !$acc              pdf_params_zm%corr_w_chi_1, pdf_params_zm%corr_w_chi_2, pdf_params_zm%corr_w_eta_1, &
    !$acc              pdf_params_zm%corr_w_eta_2, pdf_params_zm%corr_chi_eta_1, pdf_params_zm%corr_chi_eta_2, &
    !$acc              pdf_params_zm%rsatl_1, pdf_params_zm%rsatl_2, pdf_params_zm%rc_1, pdf_params_zm%rc_2, &
    !$acc              pdf_params_zm%cloud_frac_1, pdf_params_zm%cloud_frac_2, &
    !$acc              pdf_params_zm%ice_supersat_frac_1, pdf_params_zm%ice_supersat_frac_2 )


    !$acc exit data if( sclr_dim > 0 ) &
    !$acc      delete( sclr_tol, sclrm_init, &
    !$acc              sclrm_forcing, wpsclrp_sfc, sclrm, &
    !$acc              wpsclrp, sclrp2, sclrp3, sclrprtp, sclrpthlp, sclrpthvp )

    !$acc exit data if( edsclr_dim > 0 ) &
    !$acc      delete( edsclrm_init, &
    !$acc              wpedsclrp_sfc, edsclrm_forcing, edsclrm )

    !$acc exit data if( hydromet_dim > 0 ) &
    !$acc      delete( hm_metadata%l_mix_rat_hm, &
    !$acc              wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt )

    ! Close stats files
    call stats_finalize_api( ngrdcol, stats_metadata, &
                              stats_zt, stats_zm, stats_sfc, &
                              stats_lh_zt, stats_lh_sfc, &
                              stats_rad_zt, stats_rad_zm )
                              
    deallocate( stats_zt, stats_zm, stats_sfc, &
                stats_lh_zt, stats_lh_sfc, &
                stats_rad_zt, stats_rad_zm )

    ! Free memory
    if ( thlm_sponge_damp_settings%l_sponge_damping )     call finalize_tau_sponge_damp_api( thlm_sponge_damp_profile )
    if ( rtm_sponge_damp_settings%l_sponge_damping )      call finalize_tau_sponge_damp_api( rtm_sponge_damp_profile )
    if ( uv_sponge_damp_settings%l_sponge_damping )       call finalize_tau_sponge_damp_api( uv_sponge_damp_profile )
    if ( wp2_sponge_damp_settings%l_sponge_damping )      call finalize_tau_sponge_damp_api( wp2_sponge_damp_profile )
    if ( wp3_sponge_damp_settings%l_sponge_damping )      call finalize_tau_sponge_damp_api( wp3_sponge_damp_profile )
    if ( up2_vp2_sponge_damp_settings%l_sponge_damping )  call finalize_tau_sponge_damp_api( up2_vp2_sponge_damp_profile )

    if( l_t_dependent ) call finalize_t_dependent_input()

    call finalize_extended_atm( )

    ! De-allocate the arrays for the grid
    call cleanup_grid_api( gr )

    call cleanup_radiation_variables( )

    call cleanup_microphys( )
    
    if( l_input_fields .or. l_restart ) then
      call cleanup_input_fields()
    end if

#ifdef SILHS
    if ( lh_microphys_type /= lh_microphys_disabled ) then
      call latin_hypercube_2D_close_api( )
      call cleanup_latin_hypercube_arrays( )
    end if
#endif

    deallocate( um )        ! u wind
    deallocate( vm )        ! v wind

    deallocate( upwp )      ! vertical u momentum flux
    deallocate( vpwp )      ! vertical v momentum flux

    deallocate( up2 )
    deallocate( up3 )
    deallocate( vp2 )
    deallocate( vp3 )

    deallocate( thlm )      ! liquid potential temperature
    deallocate( rtm )       ! total water mixing ratio
    deallocate( wprtp )     ! w'rt'
    deallocate( wpthlp )    ! w'thl'
    deallocate( wprcp )     ! w'rc'
    deallocate( w_up_in_cloud )
    deallocate( w_down_in_cloud )
    deallocate( cloudy_updraft_frac )
    deallocate( cloudy_downdraft_frac )
    deallocate( wp2 )       ! w'^2
    deallocate( wp3 )       ! w'^3
    deallocate( wp3_zm )    ! w'^3
    deallocate( rtp2 )      ! rt'^2
    deallocate( thlp2 )     ! thl'^2
    deallocate( rtpthlp )   ! rt'thlp'

    deallocate( p_in_Pa )         ! pressure (pascals)
    deallocate( exner )           ! exner function
    deallocate( rho )             ! density: t points
    deallocate( rho_zm )          ! density: m points
    deallocate( rho_ds_zm )       ! dry, static density: m-levs
    deallocate( rho_ds_zt )       ! dry, static density: t-levs
    deallocate( invrs_rho_ds_zm ) ! inv. dry, static density: m-levs
    deallocate( invrs_rho_ds_zt ) ! inv. dry, static density: t-levs
    deallocate( thv_ds_zm )       ! dry, base-state theta_v: m-levs
    deallocate( thv_ds_zt )       ! dry, base-state theta_v: t-levs

    deallocate( thlm_forcing )    ! thlm ls forcing
    deallocate( rtm_forcing )     ! rtm ls forcing
    deallocate( um_forcing )      ! u forcing
    deallocate( vm_forcing )      ! v forcing
    deallocate( wprtp_forcing )   ! <w'r_t'> forcing (microphysics)
    deallocate( wpthlp_forcing )  ! <w'th_l'> forcing (microphysics)
    deallocate( rtp2_forcing )    ! <r_t'^2> forcing (microphysics)
    deallocate( thlp2_forcing )   ! <th_l'^2> forcing (microphysics)
    deallocate( rtpthlp_forcing ) ! <r_t'th_l'> forcing (microphysics)

    ! Variables used to track perturbed version of winds.
    deallocate( um_pert )
    deallocate( vm_pert )
    deallocate( upwp_pert )
    deallocate( vpwp_pert )

    ! Imposed large scale w
    deallocate( wm_zm )       ! momentum levels
    deallocate( wm_zt )       ! thermodynamic levels

    ! Cloud water variables
    deallocate( rcm )
    deallocate( delta_zm )
    deallocate( cloud_frac )
    deallocate( ice_supersat_frac )
    deallocate( rcm_in_layer )
    deallocate( cloud_cover )
    deallocate( invrs_tau_zm )

    ! Passive scalar variables
    ! Note that sclr_dim can be 0
    deallocate( wpsclrp_sfc )
    deallocate( sclrm )
    deallocate( sclrp2 )
    deallocate( sclrp3 )
    deallocate( sclrm_forcing )
    deallocate( sclrprtp )
    deallocate( sclrpthlp )

    deallocate( wpedsclrp_sfc )
    deallocate( edsclrm_forcing )

    deallocate( edsclrm )
    deallocate( wpsclrp )

    deallocate( sigma_sqd_w )    ! PDF width parameter (momentum levels)
    deallocate( sigma_sqd_w_zt ) ! PDF width parameter interp. to t-levs.
    deallocate( Skw_zm )         ! Skewness of w on momentum levels
    deallocate( wp2_zt )         ! wp2 interpolated to thermo. levels
    deallocate( Skw_zm_smooth )         ! Skewness of w on momentum levels
    deallocate( ug )             ! u geostrophic wind
    deallocate( vg )             ! v geostrophic wind
    deallocate( um_ref )         ! Reference u wind for nudging; Michael Falk, 17 Oct 2007
    deallocate( vm_ref )         ! Reference v wind for nudging; Michael Falk, 17 Oct 2007
    deallocate( thlm_ref )       ! Reference liquid water potential for nudging
    deallocate( rtm_ref )        ! Reference total water mixing ratio for nudging
    deallocate( thvm )           ! Virtual potential temperature
    deallocate( radht )          ! SW + LW heating rate
    deallocate( Frad )           ! radiative flux (momentum point)
    deallocate( Frad_SW_up )
    deallocate( Frad_LW_up )
    deallocate( Frad_SW_down )
    deallocate( Frad_LW_down )
    deallocate( thlprcp )   ! thl'rc'
    deallocate( thlp3 )     ! thl'^3
    deallocate( rtp3 )      ! rt'^3

    ! Buoyancy related moments
    deallocate( rtpthvp )  ! rt'thv'
    deallocate( thlpthvp ) ! thl'thv'
    deallocate( wpthvp )   ! w'thv'
    deallocate( wp2thvp )  ! w'^2thv'

    deallocate( wp2up )   ! w'^2u'
    deallocate( wp2rtp )  ! w'^2 rt'
    deallocate( wp2thlp ) ! w'^2 thl'
    deallocate( uprcp )   ! u'rc'
    deallocate( vprcp )   ! v'rc'
    deallocate( rc_coef_zm ) ! Coefficient of X'r_c' in Eq. (34)
    deallocate( wp4 )     ! w'^4
    deallocate( wpup2 )   ! w'u'^2
    deallocate( wpvp2 )   ! w'v'^2
    deallocate( wp2up2 )  ! w'^2 u'^2
    deallocate( wp2vp2 )  ! w'^2 v'^2

    deallocate( Kh_zt )  ! Eddy diffusivity coefficient: thermo. levels
    deallocate( Kh_zm )  ! Eddy diffusivity coefficient: momentum levels
    deallocate( K_hm ) ! Eddy diff. coef. for hydromets.: mom. levs.

    deallocate( em )
    deallocate( Lscale )

    deallocate( tau_zm ) ! Eddy dissipation time scale: momentum levels
    deallocate( tau_zt ) ! Eddy dissipation time scale: thermo. levels

    deallocate( Nccnm )
    deallocate( hydromet )    ! All hydrometeor mean fields
    deallocate( hydrometp2 )  ! All < h_m'^2 > fields
    deallocate( wphydrometp ) ! All < w'h_m' > fields
    deallocate( Ncm )   ! Mean cloud droplet concentration, < N_c >
    deallocate( wpNcp ) ! < w'N_c' >

    ! High-order passive scalars
    deallocate( sclrpthvp )

    deallocate( thlm_mc, rvm_mc, rcm_mc, wprtp_mc, wpthlp_mc, rtp2_mc, &
                thlp2_mc, rtpthlp_mc, hydromet_mc, Ncm_mc, hydromet_vel_zt, &
                hydromet_vel_covar_zt_impc, hydromet_vel_covar_zt_expc )

    deallocate( X_nl_all_levs, &
                lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, &
                lh_rv_clipped, lh_Nc_clipped, &
                X_mixt_comp_all_levs, lh_sample_point_weights, Nc_in_cloud )

    if ( clubb_config_flags%grid_adapt_in_time_method > no_grid_adaptation ) then

      close(iunit_grid_adaptation)

      deallocate( gr_dens_z )
      deallocate( gr_dens )

      deallocate( gr_dens_old_z )
      deallocate( gr_dens_old )

      deallocate( cumulative_Lscale )
      deallocate( cumulative_chi )
      deallocate( cumulative_richardson_num )
      deallocate( alt_term )
      deallocate( lscale_term )
      deallocate( lscale_term_time_avg )
      deallocate( norm_min_grid_dens )
      deallocate( norm_grid_dens )
      deallocate( chi_term )
      deallocate( chi_term_time_avg )
      deallocate( richardson_num_term )
      deallocate( richardson_num_term_time_avg )

    end if

    deallocate( sclr_tol )

    deallocate( p_sfc,            &
                T_sfc,            &
                fcor,             &
                fcor_y,           &
                dummy_dx,         &
                dummy_dy,         &
                deep_soil_T_in_K, &
                sfc_soil_T_in_K,  &
                veg_T_in_K,       &
                upwp_sfc_pert,    &
                vpwp_sfc_pert,    &
                sfc_elevation,    &
                zm_init,          &
                zm_top,           &
                deltaz )
                
    deallocate( rrm )
    deallocate( Nrm )

    deallocate( wp2hmp, rtphmp_zt, thlphmp_zt )
    deallocate( hydromet_pdf_params )
    deallocate( corr_array_1_n )
    deallocate( corr_array_2_n )
    deallocate( corr_cholesky_mtx_1 )
    deallocate( corr_cholesky_mtx_2 )
    deallocate( mu_x_1_n )
    deallocate( mu_x_2_n )
    deallocate( sigma_x_1_n )
    deallocate( sigma_x_2_n )
    deallocate( rfrzm )

    deallocate( wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc )

    deallocate( thlm_init, &
              rtm_init, &
              um_init, &
              vm_init, &
              ug_init, &
              vg_init, &
              wp2_init, &
              up2_init, &
              vp2_init, &
              upwp_init, &
              rcm_init, &
              wm_zt_init, &
              wm_zm_init, &
              em_init, &
              exner_init, &
              thvm_init, &
              p_in_Pa_init, &
              rho_init, &
              rho_zm_init, &
              rho_ds_zm_init, &
              rho_ds_zt_init, &
              invrs_rho_ds_zm_init, &
              invrs_rho_ds_zt_init, &
              thv_ds_zm_init, &
              thv_ds_zt_init, &
              rtm_ref_init, &
              thlm_ref_init, &
              um_ref_init, &
              vm_ref_init, &
              Ncm_init, &
              Nc_in_cloud_init, &
              Nccnm_init, &
              sclrm_init, &
              edsclrm_init )

  end subroutine clean_up_clubb

  !-----------------------------------------------------------------------
  subroutine initialize_clubb( &
                gr, ngrdcol, iunit, forcings_file_path, p_sfc, zm_init, &
                sclr_dim, edsclr_dim, sclr_idx, &
                clubb_config_flags, &
                l_modify_ic_with_cubic_int, &
                l_add_dycore_grid, &
                grid_adapt_in_time_method, &
                l_ascending_grid, fcor_y, &
                thlm, rtm, um, vm, ug, vg, wp2, up2, vp2, upwp, rcm, &
                wm_zt, wm_zm, em, exner, &
                thvm, p_in_Pa, &
                rho, rho_zm, rho_ds_zm, rho_ds_zt, &
                invrs_rho_ds_zm, invrs_rho_ds_zt, &
                thv_ds_zm, thv_ds_zt, &
                rtm_ref, thlm_ref, &
                um_ref, vm_ref, &
                Ncm, Nc_in_cloud, Nccnm, err_info, &
                deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K, &
                sclrm, edsclrm )
    ! Description:
    !   Execute the necessary steps for the initialization of the
    !   CLUBB model run.
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,            & !--------------------------------------------- Constant(s)
        zero,           &
        pi,             &
        em_min,         &
        w_tol_sqd,      &
        grav,           &
        cm3_per_m3,     &
        cloud_frac_min, &
        fstderr         !----------------------------------------------- Variable(s)

    use grid_class, only: grid ! Type

    use parameters_microphys, only: &
        Nc0_in_cloud, & !----------------------------------------------- Variable(s)
        microphys_scheme

    use parameters_radiation, only: radiation_top, rad_scheme !--------- Variable(s)

    use grid_class, only: &
      zm2zt_api, &   !----------------------------------------------------- Procedure(s)
      zt2zm_api

    use sounding, only: read_sounding !--------------------------------- Procedure(s)

    use time_dependent_input, only: &
      initialize_t_dependent_input, & !-------------------------------- Procedure(s)
      l_t_dependent !-------------------------------------------------- Variable(s)

    use extended_atmosphere_module, only: &
      determine_extended_atmos_bounds !-------------------------------- Procedure(s)

    use mpace_a, only: mpace_a_init !---------------------------------- Procedure(s)

    use sponge_layer_damping, only: &
        thlm_sponge_damp_settings,    & !------------------------------ Variable(s)
        rtm_sponge_damp_settings,     &
        uv_sponge_damp_settings,      &
        wp2_sponge_damp_settings,     &
        wp3_sponge_damp_settings,     &
        up2_vp2_sponge_damp_settings, &
        thlm_sponge_damp_profile,     &
        rtm_sponge_damp_profile,      &
        uv_sponge_damp_profile,       &
        wp2_sponge_damp_profile,      &
        wp3_sponge_damp_profile,      &
        up2_vp2_sponge_damp_profile

    use clubb_api_module, only: &
        clubb_fatal_error,              & !------------------------------ Constant(s)
        initialize_tau_sponge_damp_api    !------------------------------ Procedure(s)

    use input_names, only: &
      wm_name, &
      omega_name

    use clubb_precision, only: &
      core_rknd !------------------------------------------------------ Variable(s)

    use model_flags, only: &
      clubb_config_flags_type

    use array_index, only: &
      sclr_idx_type

    use interpolation, only: &
        lin_interp_between_grids

    use grid_adaptation_module, only: &
        setup_gr_dycore

    use err_info_type_module, only: &
      err_info_type        ! Type

    implicit none

    intrinsic :: min, max, trim, sqrt, size

    ! Input
    
    type (grid), intent(in) :: &
      gr

    integer, intent(in) :: &
      ngrdcol

    integer, intent(in) :: &
      iunit

    character(len=*), intent(in) :: &
      forcings_file_path ! Path to the .dat files containing the forcings

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      p_sfc,   & ! Pressure at the surface        [Pa]
      zm_init, & ! Initial moment. level altitude [m]
      fcor_y     ! Nontraditional Coriolis parameter [s^-1]
                 ! Meridional planetary vorticity. Proportional to cos(latitude)

    integer, intent(in) :: &
      sclr_dim, &
      edsclr_dim

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    type (clubb_config_flags_type), intent(in) :: &
      clubb_config_flags

    ! Flag for interpolating the sounding profile with Steffen's monotone cubic 
    ! method to obtain smoother initial condition profile, which is found to be 
    ! beneficial to achive a better numerical solution convergence. If this flag 
    ! is turned off, the initial conditions will be generated with linear interpolation.
    ! This is done on a case-by-case basis, since using the monotone cubic method
    ! requires a special sounding.in file with many additional sounding levels.
    logical, intent(in) :: &
      l_modify_ic_with_cubic_int

    logical, intent(in) :: &
      l_add_dycore_grid ! flag to set remapping from dycore to on or off

    integer, intent(in) :: &
      grid_adapt_in_time_method   ! Integer flag to see if grid should be adapted over time and if
                                  ! so what parameters should be used to setup the grid
                                  ! density function

    logical, intent(in) :: &
      l_ascending_grid

    ! Output
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) :: &
      exner,           & ! Exner function (thermodynamic levels)     [-] 
      thvm,            & ! Virtual potential temp. (thermo. levs.)   [K]
      rcm                ! Cloud water mixing ratio (thermo. levs.)  [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) :: &
      thlm,            & ! Grid mean of liquid water pot. temp               [K]
      rtm,             & ! Grid mean of total water mixing ratio             [kg/kg]
      um,              & ! Grid mean of eastward wind                        [m/s]
      vm,              & ! Grid mean of northbound wind                      [m/s]
      ug,              & ! Eastward geostrophic wind                         [m/s]
      vg,              & ! Northbound geostrophic wind                       [m/s]
      wm_zt,           & ! Vertical wind                                     [m/s]
      p_in_Pa,         & ! Pressure                                          [Pa]
      rho,             & ! Density (thermodynamic levels)                    [kg/m^3]
      rho_ds_zt,       & ! Dry, static density (thermo. levs.)               [kg/m^3]
      invrs_rho_ds_zt, & ! Inv. dry, static density (t-levs.)                [m^3/kg]
      thv_ds_zt,       & ! Dry, base-state theta_v (t-levs.)                 [K]
      um_ref,          & ! Initial profile of u wind                         [m/s]
      vm_ref,          & ! Initial profile of v wind                         [m/s]
      rtm_ref,         & ! Initial profile of rtm                            [kg/kg]
      thlm_ref           ! Initial profile of thlm                           [K]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(inout) :: &
      wp2,             & ! Vertical velocity variance (w'^2)                 [m^2/s^2]
      up2,             & ! East-west velocity variance (u'^2)                [m^2/s^2]
      vp2,             & ! North-south velocity variance (v'^2)              [m^2/s^2]
      upwp,            & ! Vertical east-west velocity covariance (u'w')     [m^2/s^2]
      wm_zm,           & ! Vertical wind                                     [m/s]
      em,              & ! Turbulence kinetic energy                         [m^2/s^2]
      rho_zm,          & ! Density on momentum levels                        [kg/m^3]
      rho_ds_zm,       & ! Dry, static density (moment. levs.)               [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density (m-levs.)                [m^3/kg]
      thv_ds_zm          ! Dry, base-state theta_v (m-levs.)                 [K]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) :: &
      Ncm,         & ! Mean cloud droplet conc., <N_c> (thermo. levs.)  [num/kg]
      Nc_in_cloud    ! Mean (in-cloud) cloud droplet concentration      [num/kg]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) :: &
      Nccnm    ! Cloud condensation nuclei concentration (COAMPS/MG)  [num/kg]

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    real( kind = core_rknd ), dimension(ngrdcol), intent(inout) :: &
      deep_soil_T_in_K, &
      sfc_soil_T_in_K, &
      veg_T_in_K

    ! Output
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,sclr_dim), intent(out) :: &
      sclrm      ! Standard passive scalar [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,edsclr_dim), intent(out) :: &
      edsclrm    ! Eddy diffusivity passive scalar [units vary]

    ! Local Variables

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm) :: &
      p_in_Pa_zm ! Pressure on momentum levels  [Pa]

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      rtm_sfc,          & ! Surface total water mixing ratio       [kg/kg]
      thlm_sfc            ! Surface liq. water potential temp.     [K]

    real( kind = core_rknd ) :: &
      cloud_top_height, & ! Cloud top altitude in initial profile  [m]
      domain_depth,     & ! Domain depth                           [m]
      em_max              ! Maximum value of initial subgrid TKE   [m^2/s^2]

    character(len=50) :: &
      theta_type, & ! Type of temperature sounding
      alt_type,   & ! Type of altitude sounding
      subs_type     ! Type of large-scale subsidence sounding

    integer :: i, k ! Loop index

    type( grid ) :: gr_dycore ! only used if l_add_dycore_grid = .true.

    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      p_in_Pa_dycore  ! p_in_Pa interpolated to the dycore grid

    !---- Begin code ----

    ! Read sounding information
    ! Only use first value of p_sfc and zm_init because sounding files are not configured
    ! to use multiple columns yet
    call read_sounding( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, &        ! Intent(in) 
                        gr, iunit, runtype, p_sfc(1), zm_init(1), &    ! Intent(in) 
                        clubb_config_flags%saturation_formula, & ! Intent(in) 
                        l_modify_ic_with_cubic_int, &            ! Intent(in)
                        thlm, theta_type, rtm, um, vm, ug, vg, & ! Intent(out)
                        alt_type, p_in_Pa, subs_type, wm_zt, &   ! Intent(out)
                        rtm_sfc, thlm_sfc, sclrm, edsclrm )      ! Intent(out)

    ! Covert sounding input to CLUBB compatible input
    call initialize_clubb_variables( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, & ! Intent(in)
                                    gr, alt_type, theta_type,                 & ! Intent(in)
                                    p_sfc, rtm_sfc, rtm,                      & ! Intent(in)
                                    clubb_config_flags%saturation_formula,    & ! Intent(in)
                                    thlm, p_in_Pa, p_in_Pa_zm,                & ! Intent(inout)
                                    exner, rho, rho_zm,                       & ! Intent(out)
                                    rcm, thvm, rho_ds_zm,                     & ! Intent(out)
                                    rho_ds_zt, invrs_rho_ds_zm,               & ! Intent(out)
                                    invrs_rho_ds_zt, thv_ds_zm,               & ! Intent(out)
                                    thv_ds_zt, sclrm, edsclrm )                 ! Intent(out)

    if ( trim( rad_scheme ) == "bugsrad" ) then
      ! Currently clubb does not support different grid heights, use only the first column to
      ! determine the size of the extended atmosphere
      call determine_extended_atmos_bounds( gr%nzm, gr%zt(1,:),              & ! Intent(in)
                                          gr%zm(1,:), gr%dzt(1,:), p_in_Pa_zm(1,:), &   ! Intent(in)
                                          radiation_top,             &   ! Intent(in)
                                          extended_atmos_bottom_level, & ! Intent(out)
                                          extended_atmos_top_level,    & ! Intent(out)
                                          extended_atmos_range_size,   & ! Intent(out)
                                          lin_int_buffer )               ! Intent(out)

    else
      ! lin_int_buffer et al. are set to zero in clubb_model_settings.
    end if
    
    ! Determine initial value cloud droplet number concentration when Nc
    ! is predicted.
    Nc_in_cloud = Nc0_in_cloud / rho

    do i = 1, ngrdcol
      do k = 1, gr%nzt, 1

        if ( rcm(i,k) > zero ) then

            ! The initial profile at this level is entirely saturated (due to
            ! constant moisture and temperature over the level).  The level is
            ! entirely cloudy.
            Ncm(i,k) = Nc_in_cloud(i,k)

        else ! rcm = 0

            ! The initial profile at this level is entirely unsaturated (due to
            ! constant moisture and temperature over the level).  The level is
            ! entirely clear.
            Ncm(i,k) = Nc_in_cloud(i,k) * cloud_frac_min

        endif

      end do ! k = 1, gr%nzt, 1
    end do

    select case ( trim( microphys_scheme ) )

    case ( "coamps" )
        ! Initialize Nccnm as in COAMPS-LES
        Nccnm(:,1:gr%nzt) &
        = 30.0_core_rknd &
          * ( one + exp( -gr%zt(:,1:gr%nzt) / 2000.0_core_rknd ) ) &
          * cm3_per_m3 / rho    ! Known magic number

    case default

      Nccnm = zero

    end select


    ! Initialize imposed w
    select case ( trim( subs_type ) ) ! Perform different operations based off
      !                                   the sounding file
    case ( wm_name )

      wm_zm = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, wm_zt )
      wm_zm(:,1) = 0.0_core_rknd
      wm_zm(:,gr%nzm) = 0.0_core_rknd

    case ( omega_name )

      do k=1,gr%nzt
        wm_zt(:,k) = -wm_zt(:,k) / ( grav*rho(:,k) )
      end do

      wm_zt(:,gr%nzt) = 0.0_core_rknd

      wm_zm = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, wm_zt )
      wm_zm(:,gr%nzm) = 0.0_core_rknd

    case default ! This should not happen

      wm_zt = 0.0_core_rknd
      wm_zm = 0.0_core_rknd

    end select

    ! Initialize damping
    if ( thlm_sponge_damp_settings%l_sponge_damping ) then
      call initialize_tau_sponge_damp_api( gr, gr%nzt, dt_main, gr%zt(1,:), & ! Intent(in)
                                            thlm_sponge_damp_settings, & ! Intent(in)
                                            thlm_sponge_damp_profile )   ! Intent(out)
    endif

    if ( rtm_sponge_damp_settings%l_sponge_damping ) then
      call initialize_tau_sponge_damp_api( gr, gr%nzt, dt_main, gr%zt(1,:), & ! Intent(in)
                                            rtm_sponge_damp_settings, & ! Intent(in)
                                            rtm_sponge_damp_profile )   ! Intent(out)
    endif

    if ( uv_sponge_damp_settings%l_sponge_damping ) then
      call initialize_tau_sponge_damp_api( gr, gr%nzt, dt_main, gr%zt(1,:), & ! Intent(in)
                                            uv_sponge_damp_settings, & ! Intent(in)
                                            uv_sponge_damp_profile )   ! Intent(out)
    endif

    if ( wp2_sponge_damp_settings%l_sponge_damping ) then
      call initialize_tau_sponge_damp_api( gr, gr%nzm, dt_main, gr%zm(1,:), & ! Intent(in)
                                            wp2_sponge_damp_settings, & ! Intent(in)
                                            wp2_sponge_damp_profile )   ! Intent(out)
    endif

    if ( wp3_sponge_damp_settings%l_sponge_damping ) then
      call initialize_tau_sponge_damp_api( gr, gr%nzt, dt_main, gr%zt(1,:), & ! Intent(in)
                                            wp3_sponge_damp_settings, & ! Intent(in)
                                            wp3_sponge_damp_profile )   ! Intent(out)
    endif

    if ( up2_vp2_sponge_damp_settings%l_sponge_damping ) then
      call initialize_tau_sponge_damp_api( gr, gr%nzm, dt_main, gr%zm(1,:), & ! Intent(in)
                                            up2_vp2_sponge_damp_settings, & ! Intent(in)
                                            up2_vp2_sponge_damp_profile )   ! Intent(out)
    endif

    ! Initialize Time Dependent Input
    
    if( l_t_dependent ) then

      if ( clubb_config_flags%l_add_dycore_grid ) then
        ! get dycore grid
        call setup_gr_dycore( iunit, ngrdcol, &                 ! Intent(in)
                              gr%zm(:,1), gr%zm(:,gr%nzm), &    ! Intent(in)
                              l_ascending_grid, &               ! Intent(in)
                              gr_dycore, &                      ! Intent(out)
                              err_info )                        ! Intent(inout)

        allocate( p_in_Pa_dycore(ngrdcol,gr_dycore%nzt) )

        ! get pressure for the dycore grid
        do i = 1, ngrdcol
          p_in_Pa_dycore(i,:) = lin_interp_between_grids( gr_dycore%nzt, gr%nzm, &
                                                          gr_dycore%zt, gr%zm, &
                                                          p_in_Pa_zm(i,:) )
        end do
        ! if we want to simulate the forcings from the host model on the dycore grid,
        ! we first read in all the forcings from the file to the dycore grid and every
        ! timestep remap the results to the current physics grid, so we need to pay
        ! attention on how the forcings are read in from object t_dependent_forcing_data
        ! in time_dependent_input
        call initialize_t_dependent_input &
                    ( iunit, runtype, gr_dycore%nzt, gr_dycore%zt(1,:), p_in_Pa_dycore(1,:), &
                    clubb_config_flags%l_add_dycore_grid, clubb_config_flags%grid_adapt_in_time_method )

        deallocate( p_in_Pa_dycore )
      else
        ! no simulating forcings input from host model on dycore grid
        ! l_add_dycore_grid=.false.
        call initialize_t_dependent_input &
                      ( iunit, runtype, gr%nzt, gr%zt(1,:), p_in_Pa(1,:), &
                      clubb_config_flags%l_add_dycore_grid, clubb_config_flags%grid_adapt_in_time_method )
      end if

    else

      if ( clubb_config_flags%l_add_dycore_grid .and. trim( runtype ) /= 'atex' ) then
        if ( trim( runtype ) /= 'gabls2' ) then
          ! works for gabls2 since forcings are just set to 0, so it doesn't make any difference
          ! if we first set the forcings on dycore grid to zero and then remap to the physics
          ! grid or if we just directly set the forcings on the physics grid to 0,
          ! at least if remapping scheme is consistent
          ! (which is the case for ullrich remapping scheme (method 1) and ppm (method 2))
          write(fstderr,*) 'WARNING! The option l_add_dycore_grid=.true. or ', &
                            'grid_adapt_in_time_method>0 can currently only ', &
                            'be used for cases with forcings from an input file and for the atex ', &
                            'case, so for this case l_add_dycore_grid must be set to .false. ', &
                            'and grid_adapt_in_time_method=0'
          err_info%err_code = clubb_fatal_error
          write(fstderr,*) 'Flags have to be changed like mentioned.'
          return
        end if

      end if

    end if
    
    ! Initialize TKE and other fields as needed

    select case ( trim( runtype ) )

      ! Generic case
    case ( "generic" )

      em = 1.0_core_rknd
      em(:,gr%nzm) = em_min

      ! GCSS BOMEX
    case ( "bomex", "ekman", "atex_long" )

!---> Reduction of initial sounding for stability
!         do k = 1, gr%nz
!            em(k) = 1.0_core_rknd - (gr%zm(1,k)/3000.0_core_rknd)
!            if ( em(k) < em_min ) then
!               em(k) = em_min
!            end if
!         end do
!         em(1) = em(2)
!         em(gr%nz) = em(gr%nz-1)
!<--- End reduction of initial sounding for stability 24 Jan 07

      em(:,:) = em_min

      ! GCSS ARM
    case ( "arm" )

!---> Reduction of initial sounding for stability
!         do k = 1, gr%nz
!            if ( gr%zm(1,k) < 150.0_core_rknd ) then
!               em(k) = ( 0.15_core_rknd * (1.0_core_rknd - gr%zm(1,k)/150.0_core_rknd) ) 
!                       / rho_zm(k)
!            else
!               em(k) = em_min
!            end if
!         end do
!         em(1) = em(2)
!         em(gr%nz) = em(gr%nz-1)
!<--- End reduction of initial sounding for stability 24 Jan 07

      em(:,:) = em_min

      ! June 27 1997 ARM case
    case ( "arm_97" )

      em = 1.0_core_rknd
      em(:,gr%nzm) = em_min

      ! twp_ice
    case ( "twp_ice" )

      em = 1.0_core_rknd
      em(:,gr%nzm) = em_min

      ! March 2000 ARM case
    case ( "arm_0003" )

      em = 1.0_core_rknd
      em(:,gr%nzm) = em_min

      ! 3 year ARM case
    case ( "arm_3year" )

      em = 1.0_core_rknd
      em(:,gr%nzm) = em_min

      ! cloud feedback cases
    case ( "cloud_feedback_s6", "cloud_feedback_s6_p2k",   &
            "cloud_feedback_s11", "cloud_feedback_s11_p2k", &
            "cloud_feedback_s12", "cloud_feedback_s12_p2k" )

      em = 1.0_core_rknd
      em(:,gr%nzm) = em_min

      ! ASTEX_A209 case 16 Jul, 2010 kcwhite
    case ( "astex_a209" )

      cloud_top_height = 700._core_rknd
      em_max = 1.0_core_rknd

      do i = 1, ngrdcol
        do k=1,gr%nzm
          if ( gr%zm(i,k) < cloud_top_height ) then
            em(i,k) = em_max
          else
            em(i,k) = em_min
          end if
        end do
      end do
      em(:,1) = em(:,2)
      em(:,gr%nzm) = em_min

      ! GCSS FIRE Sc
    case ( "fire" )

      cloud_top_height = 700._core_rknd ! 700 m is the top of the cloud in FIRE
      do i = 1, ngrdcol
        do k=1,gr%nzm
          if ( gr%zm(i,k) < cloud_top_height ) then
            !em(k) = 1._core_rknd
            em(i,k) = 4.5_core_rknd
          else
            em(i,k) = em_min
          end if
        end do
      end do
      em(:,1) = em(:,2)
      em(:,gr%nzm) = em_min

      ! GCSS ATEX
    case ( "atex" )

      um = max( um, -8._core_rknd ) ! Known magic number (Stevens, et al. 2001, eq. 1)

!---> Reduction of initial sounding for stability
!         do k = 1, gr%nz
!           em(k) = 1.0_core_rknd - (gr%zm(1,k)/3000.0_core_rknd)
!           if ( em(k) < em_min ) then
!             em(k) = em_min
!           end if
!         end do
!         em(1) = em(2)
!         em(gr%nz) = em(gr%nz-1)
!<--- End reduction of initial sounding for stability 24 Jan 07

      em(:,:) = em_min

      ! GCSS DYCOMS II RF01
    case ( "dycoms2_rf01" )
      cloud_top_height = 800._core_rknd ! 800 m is the top of the cloud in RF01
      do i = 1, ngrdcol
        do k=1,gr%nzm
          if ( gr%zm(i,k) < cloud_top_height ) then
            !em(k) = 0.5_core_rknd
            em(i,k) = 1.1_core_rknd
          else
            em(i,k) = em_min
          end if
        end do
      end do
      em(:,1) = em(:,2)
      em(:,gr%nzm) = em_min

      ! GCSS DYCOMS II RF02
    case ( "dycoms2_rf02" )

      em = 1.0_core_rknd
      em(:,gr%nzm) = em_min

      ! Brian for Nov. 11 altocumulus case.
    case ( "nov11_altocu" )

      ! Vince Larson reduced initial forcing.  4 Nov 2005
!          em = 1.0_core_rknd
!          em = 0.1_core_rknd
      ! 4150 + 2800 m is the top of the cloud in Nov11
      cloud_top_height = 2800._core_rknd + gr%zm(1,1) ! Known magic number
      do i = 1, ngrdcol
        do k=1,gr%nzm
          if ( gr%zm(i,k) < cloud_top_height ) then

            ! Modification by Adam Smith, 08 April 2008
            ! Reducing the value of em appears to reduce error in the
            ! updated Nov.11 case
            em(i,k) = 0.01_core_rknd
            ! End of ajsmith4's modification

          else
            em(i,k) = em_min
          end if
        end do
      end do
      em(:,1) = em(:,2)
      em(:,gr%nzm) = em_min
      ! End Vince Larson's change.

      ! Adam Smith addition for June 25 altocumulus case.
    case ( "jun25_altocu" )

      ! Vince Larson reduced initial forcing.  4 Nov 2005
!          em = 1.0_core_rknd
!          em = 0.1_core_rknd
!          do k=1,gr%nz
!            if ( gr%zm(1,k) < 1400._core_rknd ) then
!               em(k) = 0.1_core_rknd
!            else
!               em(k) = em_min
!            end if
!          end do

      ! Note: em_min = 1.0e-6, defined in constants.F
      ! Adam Smith, 28 June 2006
      ! Note: now em_min = 1.5 * w_tol_sqd
      ! Brian Griffin;  Nov. 26, 2008.
      do k = 1, gr%nzm
        em(:,k) = 0.01_core_rknd
      end do

      em(:,1) = em(:,2)
      ! End Vince Larson's change.

      em(:,gr%nzm) = em_min

      ! End of ajsmith4's addition

      ! Adam Smith addition for CLEX-9: Nov. 02 altocumulus case.
    case ( "clex9_nov02" )

      ! Vince Larson reduced initial forcing.  4 Nov 2005
!          em = 1.0_core_rknd
!          em = 0.1_core_rknd
      ! 4150 + 1400 m is the top of the cloud in Nov11
      cloud_top_height = 2200._core_rknd + gr%zm(1,1) ! Known magic number
      do i = 1, ngrdcol
        do k=1,gr%nzm
          if ( gr%zm(i,k) < cloud_top_height ) then
            em(i,k) = 0.01_core_rknd
          else
            em(i,k) = em_min
          end if
        end do
      end do
      em(:,1) = em(:,2)
      ! End Vince Larson's change.
      em(:,gr%nzm) = em_min

      ! End of ajsmith4's addition

      ! Adam Smith addition for CLEX-9: Oct. 14 altocumulus case.
    case ( "clex9_oct14" )

      ! Vince Larson reduced initial forcing.  4 Nov 2005
!          em = 1.0_core_rknd
!          em = 0.1_core_rknd
      ! 4150 + 1400 m is the top of the cloud in Nov11
      cloud_top_height = 3500._core_rknd + gr%zm(1,1) ! Known magic number
      do i = 1, ngrdcol
        do k=1,gr%nzm
          if ( gr%zm(i,k) < cloud_top_height ) then
            em(:,k) = 0.01_core_rknd
          else
            em(:,k) = em_min
          end if
        end do
      end do
      em(:,1) = em(:,2)
      ! End Vince Larson's change.
      em(:,gr%nzm) = em_min

      ! End of ajsmith4's addition

    case ( "lba" )

      em = 0.1_core_rknd
      em(:,gr%nzm) = em_min

      ! Michael Falk for mpace_a Arctic Stratus case.
    case ( "mpace_a" )

      cloud_top_height = 2000._core_rknd
      em_max = 1.0_core_rknd

      do i = 1, ngrdcol
        do k=1,gr%nzm
          if ( gr%zm(i,k) < cloud_top_height ) then
            em(:,k) = em_max
          else
            em(:,k) = em_min
          end if
        end do
      end do
      em(:,1) = em(:,2)
      em(:,gr%nzm) = em_min

      call mpace_a_init( iunit, forcings_file_path )

      ! Michael Falk for mpace_b Arctic Stratus case.
    case ( "mpace_b" )

      cloud_top_height = 1300._core_rknd ! 1300 m is the cloud top in mpace_b.
                                          ! Michael Falk 17 Aug 2006
      em_max = 1.0_core_rknd

      do i = 1, ngrdcol
        do k=1,gr%nzm

          if ( gr%zm(i,k) < cloud_top_height ) then
            em(:,k) = em_max
          else
            em(:,k) = em_min
          end if
        enddo
      enddo
      em(:,1) = em(:,2)
      em(:,gr%nzm) = em_min

      ! Brian Griffin for COBRA CO2 case.
    case ( "cobra" )

      em = 0.1_core_rknd
      em(:,gr%nzm) = em_min

      ! Michael Falk for RICO tropical cumulus case, 13 Dec 2006
    case ( "rico" )

      cloud_top_height = 1500._core_rknd
      em_max = 1.0_core_rknd
      do i = 1, ngrdcol
        do k=1,gr%nzm
          if ( gr%zm(i,k) < cloud_top_height ) then
            em(:,k) = em_max
          else
            em(:,k) = em_min
          end if
        enddo
      enddo

      em(:,1) = em(:,2)
      em(:,gr%nzm) = em_min

      ! Michael Falk for GABLS2 case, 29 Dec 2006
    case ( "gabls2" )

      cloud_top_height = 800._core_rknd  ! per GABLS2 specifications
      em_max = 0.5_core_rknd
      do i = 1, ngrdcol
        do k=1,gr%nzm
          if ( gr%zm(i,k) < cloud_top_height ) then
            em(:,k) = em_max * (1._core_rknd - (gr%zm(i,k)/cloud_top_height))
          else
            em(:,k) = em_min
          end if
        end do
      end do

      em(:,1) = em(:,2)
      em(:,gr%nzm) = em_min

    case ( "gabls3_night" )
      em = 1.0_core_rknd

    case ( "gabls3" )
      em = 1.0_core_rknd
      em(:,gr%nzm) = em_min

      veg_T_in_K = 300._core_rknd
      sfc_soil_T_in_K = 300._core_rknd
      deep_soil_T_in_K = 288.58_core_rknd

    ! Set the amplitude of the harmonic oscillator in the "coriolis_test".
    ! The maximum amplitude is 2 * w_tol_sqd at the middle vertical level.
    ! Hing Ong, 8 August 2025
    case ( "coriolis_test" )

      do i = 1, ngrdcol
        domain_depth = gr%zm(i,gr%nzm) - gr%zm(i,1)
        do k=1,gr%nzm
        em(i,k) = sin( pi * gr%zm(i,k) / domain_depth ) * w_tol_sqd * 6.0_core_rknd
        end do
      end do

    case default

      em = em_min

    end select

    ! End Initialize TKE and other fields as needed

    !!!! Initialize w'^2 based on initial TKE !!!

    if ( clubb_config_flags%l_tke_aniso ) then

      ! TKE:  em = (1/2) * ( w'^2 + u'^2 + v'^2 )
      ! Evenly divide TKE into its component
      ! contributions (w'^2, u'^2, and v'^2).

      wp2 = (2.0_core_rknd/3.0_core_rknd) * em
      up2 = (2.0_core_rknd/3.0_core_rknd) * em
      vp2 = (2.0_core_rknd/3.0_core_rknd) * em
      upwp = zero

      ! The "coriolis_test" is a benchmark against analytic solutions where up2, wp2, and upwp
      ! interact as a harmonic oscillator analogous to the Foucault pendulum.
      ! Here are the solutions to test the implementation of the nontraditional Coriolis terms:
      !              upwp = amplitude(grdcol,zm) * sin( 2 * fcor_y * time )
      ! ( up2 - wp2 ) / 2 = amplitude(grdcol,zm) * cos( 2 * fcor_y * time )
      ! Here, set up2, wp2, upwp at time = 0.
      ! The setting can be modified to test the implementation of the traditional Coriolis terms,
      ! where the solutions are as follows:
      ! upwp = amplitude(grdcol,zm) * cos(   fcor * time )
      ! vpwp = amplitude(grdcol,zm) * sin( - fcor * time )
      ! Hing Ong, 8 August 2025
      if ( trim( runtype ) == "coriolis_test" ) then

        wp2  = (1.0_core_rknd/3.0_core_rknd) * em + w_tol_sqd
        up2  = (3.0_core_rknd/3.0_core_rknd) * em + w_tol_sqd
        vp2  = (2.0_core_rknd/3.0_core_rknd) * em + w_tol_sqd
        em   = em + 1.5_core_rknd * w_tol_sqd

        ! Advance upwp by half a time step like the leapfrog time integration.
        ! This improves the benchmarking when initializing at time = 0, where the time tendency is
        ! large for upwp but small for up2 and wp2.
        ! Hing Ong, 8 August 2025
        do i = 1, ngrdcol
          do k=1,gr%nzm
            upwp(i,k) = 0.5_core_rknd * dt_main * fcor_y(i) * ( up2(i,k) - wp2(i,k) )
          end do
        end do

      end if

    else

      ! TKE:  em = (3/2) * w'^2

      wp2  = (2.0_core_rknd/3.0_core_rknd) * em
      up2  = zero
      vp2  = zero
      upwp = zero

    end if ! l_tke_aniso

    ! Moved this to be more general -dschanen July 16 2007
    if ( clubb_config_flags%l_uv_nudge .or. uv_sponge_damp_settings%l_sponge_damping ) then
      um_ref = um ! Michael Falk addition for nudging code.  27 Sep/1 Nov 2006
      vm_ref = vm ! ditto
    else
      um_ref = zero
      vm_ref = zero
    end if

    if ( thlm_sponge_damp_settings%l_sponge_damping ) then
      thlm_ref = thlm ! Added for nudging code
    else
      thlm_ref = zero
    end if

    if ( rtm_sponge_damp_settings%l_sponge_damping ) then
      rtm_ref = rtm
    else
      rtm_ref = zero
    end if

    return
  end subroutine initialize_clubb

  !-----------------------------------------------------------------------------
  subroutine initialize_clubb_variables( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, &
                                          gr, alt_type, theta_type, &
                                          p_sfc, rtm_sfc, rtm, & !thlm_sfc, &
                                          saturation_formula, &
                                          thlm, p_in_Pa, p_in_Pa_zm, &
                                          exner, rho, rho_zm, &
                                          rcm, thvm, rho_ds_zm, &
                                          rho_ds_zt, invrs_rho_ds_zm, &
                                          invrs_rho_ds_zt, thv_ds_zm, &
                                          thv_ds_zt, sclrm, edsclrm )

    ! Description:
    !   Given inital sounding data (already interpolated onto model thermodynamic
    !   levels) for rtm, thlm (which can be temperature, theta, or theta_l at this
    !   point), and pressure (in the case that the sounding is given in pressure
    !   coordinates), as well as surface data on surface pressure, rtm at the
    !   surface, and thlm (temp., theta, or theta_l) at the surface, calculate
    !   many initial profiles of variables used in CLUBB.  Pressure is calculated
    !   (in the case that the sounding is given in altitude coordinates), as well
    !   as exner and density.  Initial rcm, theta, and theta_l are calculated.
    !   Additionally, the dry profiles (dry densities and dry, base-state theta_v)
    !   for the anelastic equation set are calculated.

    ! References:
    !   None
    !---------------------------------------------------------------------------

    use grid_class, only: grid ! Type

    use grid_class, only: &
        zt2zm_api ! Procedure(s)

    use constants_clubb, only: & ! Constant(s)
        one,   & ! 1
        Rd,    & ! Gas constant for dry air          [J/(kg K)]
        Cp,    & ! Specific heat of dry air          [J/(kg K)]
        Lv,    & ! Latent heat of vaporization       [J/kg]
        ep2,   & ! R_v/R_d                           [-]
        ep1,   & ! [ 1 - (R_d/R_v) ] / (R_d/R_v)     [-]
        kappa, & ! R_d/C_p                           [-]
        p0,    & ! Reference pressure of 10^5 Pa     [Pa]
        zero_threshold, & ! A threshold value of 0   [units vary]
        fstderr    ! Output to error output stream

    use calc_pressure, only: &
        calculate_thvm    ! Procedure(s)

    use model_flags, only: &
        l_use_boussinesq

    use input_names, only: &
        z_name, & !----------------------------------- Variable(s)
        pressure_name, &
        temperature_name, &
        theta_name, &
        thetal_name

    use hydrostatic_module, only: &
        hydrostatic !--------------------------------- Procedure(s)

    use saturation, only: &
        sat_mixrat_liq_api, & !--------------------------- Procedure(s)
        rcm_sat_adj

    use array_index, only: &
        sclr_idx_type

    use clubb_precision, only: &
        core_rknd !----------------------------------- Variable(s)

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      sclr_dim, &
      edsclr_dim

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    type (grid), intent(in) :: gr

    character(len=*), intent(in) :: &
      alt_type,   & ! Type of altitude sounding (altitude or pressure)
      theta_type    ! Type of temperature sounding (temp., theta, or theta_l)

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      p_sfc,     & ! Surface pressure                              [Pa]
      rtm_sfc !,& ! Surface total water mixing ratio               [kg/kg]
      !thlm_sfc    ! Surface liquid water potential temperature     [K]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(in) :: &
      rtm    ! Total water mixing ratio (thermodynamic levels)    [kg/kg]

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    !--------------------- InOut Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) :: &
      thlm       ! Liquid water potential temperature (thermo. levs.)  [K]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) :: &
      p_in_Pa    ! Pressure (thermodynamic levels)                     [Pa]

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(out) :: &
      exner,           & ! Exner function (thermodynamic levels)     [-]
      thvm,            & ! Virtual potential temp. (thermo. levs.)   [K]
      rcm                ! Cloud water mixing ratio (thermo. levs.)  [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(out) :: &
      rho,             & ! Density (thermodynamic levels)            [kg/m^3]
      rho_ds_zt,       & ! Dry, static density (thermodynamic levs.) [kg/m^3]
      invrs_rho_ds_zt, & ! Inverse dry, static density (t-levs.)     [m^3/kg]
      thv_ds_zt          ! Dry, base-state theta_v (thermo. levels)  [K]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(out) :: &
      p_in_Pa_zm,      & ! Pressure (momentum levels)                [Pa]
      rho_zm,          & ! Density on momentum levels                [kg/m^3]
      rho_ds_zm,       & ! Dry, static density (momentum levels)     [kg/m^3]
      invrs_rho_ds_zm, & ! Inverse dry, static density (m-levs.)     [m^3/kg]
      thv_ds_zm          ! Dry, base-state theta_v (momentum levels) [K]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,sclr_dim), intent(inout) :: &
      sclrm  ! Standard passive scalar           [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,edsclr_dim), intent(inout) :: &
      edsclrm ! Eddy-diffusivity passive scalar   [units vary]

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol) :: &
      pd_sfc, & ! Dry surface pressure                [Pa]
      rv_sfc    ! Surface water vapor mixing ratio    [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt) :: &
      thm             ! Potential temperature (thermodynamic levels)   [K]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt) :: &
      th_dry,       & ! Dry potential temperature (thermo. levels)     [K]
      p_dry,        & ! Dry air pressure (thermodynamic levels)        [Pa]
      exner_dry,    & ! Exner of dry air (thermodynamic levels)        [-]
      rho_dry         ! Dry air density (thermodynamic levels)         [kg/m^3]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm) :: &
      tmp, &
      thm_zm, &
      exner_zm,     & ! Exner on momentum levels                       [-]
      th_dry_zm,    & ! Dry potential temperature on momentum levels   [K]
      p_dry_zm,     & ! Dry air pressure on momentum levels            [Pa]
      exner_dry_zm, & ! Exner of dry air on momentum levels            [-]
      rho_dry_zm      ! Dry air density on momentum levels             [kg/m^3]

    integer :: i, k   ! Array index

    !--------------------- Begin Code ---------------------

    ! The value of rtm at the surface is output from the sounding, as long as
    ! the initial sounding extendeds to the model surface (at gr%zm(1,1)).
    do i = 1, ngrdcol

      if ( rtm_sfc(i) < 0.0_core_rknd ) then
        ! The sounding doesn't extended to the surface, so rtm_sfc is set to a
        ! negative number.  Use rtm(1) as rv_sfc.
        rv_sfc(i) = rtm(i,1)
      else ! rtm_sfc >= 0.0_core_rknd
        ! The sounding does extended to the surface, so rtm_sfc is the initial value
        ! of total water mixing ratio at the surface.
        rv_sfc(i) = rtm_sfc(i)
      end if

      ! Calculate dry surface pressure from surface pressure and surface water
      ! vapor mixing ratio, such that p_d = p / [ 1 + (R_v/R_d)*r_v ].
      pd_sfc(i) = p_sfc(i) / ( one + ep2 * rv_sfc(i) )

    end do

    if ( theta_type == temperature_name ) then

      if ( trim( alt_type ) == z_name ) then

        write(fstderr,*) 'Interpetation of sounding files with z as the ', &
                          'independent variable and absolute temperature ', &
                          'as the temperature variable has not been ', &
                          'implemented.  Either specify pressure as the ', &
                          'independent variable or thm/thlm as the ', &
                          'temperature variable.'
        error stop "Fatal error."

      elseif ( trim( alt_type ) == pressure_name ) then

        ! The variable "thlm" actually contains temperature (in Kelvin) at
        ! this point.

        ! Calculate initial potential temperature from temperature and exner.
        ! Again, the variable "thlm" actually contains temperature before this
        ! calculation.  After this calculation, the variable "thlm" will
        ! actually contain potential temperature.
        do k = 1, gr%nzt, 1
          do i = 1, ngrdcol
            thlm(i,k) = thlm(i,k) / ( p_in_Pa(i,k) / p0 )**kappa
          end do
        end do

      else

        error stop "Invalid sounding vertical-coordinate variable"

      endif

    endif

    ! At this point, thlm may actually contain either theta or theta-l.

    ! Calculate approximate thvm, given initial thm/thlm and rtm.
    !
    ! The exact form of the equaton for theta_v (with only water vapor
    ! included) is:
    !
    ! theta_v = theta * [ ( 1 + (R_v/R_d)*r_v ) / ( 1 + r_v ) ];
    !
    ! which can be rearranged as:
    !
    ! theta_v = theta * [ 1 + { (R_v/R_d) - 1 } * { r_v / ( 1 + r_v ) } ].
    !
    ! This can be approximated by using r_t instead of r_v.  The value of
    ! initial thvm (including water vapor and cloud water) will be
    ! recalculated more accurately once the value of initial r_c has been
    ! computed.
    !
    ! As stated above, thlm may actually contain either theta or theta-l at
    ! this point.  The equation above is based on theta.  If the variable
    ! 'thlm' actually does contain theta-l at this point, then that is
    ! another source of inaccuracy in the thvm approximation.  However, theta
    ! cannot be found from theta-l until r_c has been computed.  That being
    ! said, calling subroutine hydrostatic with an approximate thvm, rather
    ! than thm/thlm, will allow for a better calculation of exner, pressure,
    ! and density in subroutine hydrostatic.  While exner, pressure, and
    ! density are all recalculated more accurately later, in the second call
    ! to subroutine hydrostatic -- after a more accurate thvm has been
    ! calculated -- it is still important to obtain the best values for
    ! pressure and exner from the first call to subroutine hydrostatic.  This
    ! is important to allow the ensuing computation of initial r_c is done as
    ! accurately as possible.
    do k = 1, gr%nzt
      do i = 1, ngrdcol
        thvm(i,k) = thlm(i,k) * ( one + ep1 * ( rtm(i,k) / ( one + rtm(i,k) ) ) )
      end do
    end do

    ! Compute approximate pressure, exner, and density using an approximate
    ! value of theta_v.
    call hydrostatic( ngrdcol, gr, thvm, p_sfc, & ! Intent(in)
                      p_in_Pa, p_in_Pa_zm,      & ! Intent(out)
                      exner, exner_zm,          & ! Intent(out)
                      rho, rho_zm               ) ! Intent(out)

    select case( trim( theta_type ) )

    case ( theta_name, temperature_name )

      ! The variable "thlm" actually contains potential temperature (theta)
      ! at this point.
      do k = 1, gr%nzt
        do i = 1, ngrdcol
          thm(i,k) = thlm(i,k)
        end do
      end do

      ! Calculate cloud water mixing ratio based on total water mixing ratio
      ! and saturation mixing ratio, which based total pressure and
      ! temperature, which is equal to theta * exner.
      do k = 1, gr%nzt
        do i = 1, ngrdcol
          rcm(i,k) &
          = max( rtm(i,k) &
                  - sat_mixrat_liq_api( p_in_Pa(i,k), thm(i,k) * exner(i,k), saturation_formula ), &
                  zero_threshold )
        end do
      end do

      ! Compute initial theta_l based on the theta profile (currently stored
      ! in variable thlm) and cloud water mixing ratio (rcm), such that:
      !  theta_l = theta - [Lv/(Cp*exner)]*rcm.
      do k = 1, gr%nzt
        do i = 1, ngrdcol
          thlm(i,k) = thlm(i,k) - Lv/(Cp*exner(i,k)) * rcm(i,k)
        end do
      end do

      ! Testing of passive scalars
      if ( sclr_idx%iisclr_thl > 0 ) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            sclrm(i,k,sclr_idx%iisclr_thl) = thlm(i,k)
          end do
        end do
      endif
      if ( sclr_idx%iiedsclr_thl > 0 ) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            edsclrm(i,k,sclr_idx%iiedsclr_thl) = thlm(i,k)
          end do
        end do
      endif

    case ( thetal_name )

        ! The value of variable thlm that was just used to call subroutine
        ! hydrostatic is indeed thlm.

        ! Find theta based on the given profile of theta_l.  If the profile
        ! is unsaturated, then theta = theta_l.  If this initial profile is
        ! saturated at any level, then initial r_c must be determined using an
        ! iterative method involving theta_l, r_t, pressure, and exner.  Once
        ! initial r_c is found, initial theta can be found, such that:
        ! theta = theta_l + [Lv/(Cp*exner)]*rcm.
        
        ! Find mean cloud water mixing ratio.
        do k = 1, gr%nzt, 1
          do i = 1, ngrdcol
            ! Compute cloud water mixing ratio using an iterative method.
            rcm(i,k) = rcm_sat_adj( thlm(i,k), rtm(i,k), p_in_Pa(i,k), exner(i,k), &
                                    saturation_formula )
          end do
        end do

        ! Compute initial theta.
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            thm(i,k) = thlm(i,k) + Lv/(Cp*exner(i,k)) * rcm(i,k)
        end do
        end do

        ! Testing of passive scalars
        if ( sclr_idx%iisclr_thl > 0 ) then
          do k = 1, gr%nzt
            do i = 1, ngrdcol
              sclrm(i,k,sclr_idx%iisclr_thl) = thlm(i,k)
            end do
          end do
        endif

        if ( sclr_idx%iiedsclr_thl > 0 ) then
          do k = 1, gr%nzt
            do i = 1, ngrdcol
              edsclrm(i,k,sclr_idx%iiedsclr_thl) = thlm(i,k)
            end do
          end do
        endif


    case default

        write(fstderr,*) "Invalid theta_type: ", theta_type
        error stop


    end select


    ! Now, compute initial thvm, given initial thm, rtm, and rcm.
    !
    ! The exact form of the equaton for theta_v (with only water vapor
    ! included) is:
    !
    ! theta_v = theta * [ ( 1 + (R_v/R_d)*r_v ) / ( 1 + r_v ) ];
    !
    ! which can be rearranged as:
    !
    ! theta_v = theta * [ 1 + { (R_v/R_d) - 1 } * { r_v / ( 1 + r_v ) } ].
    !
    ! The exact form of the equation for theta_v (including water vapor and
    ! cloud water) is:
    !
    ! theta_v = theta * [ ( 1 + (R_v/R_d)*r_v ) / ( 1 + r_v + r_c ) ];
    !
    ! which can be rearranged as:
    !
    ! theta_v = theta * [ 1 + { (R_v/R_d) - 1 } * { r_v / ( 1 + r_v + r_c ) }
    !                     - { r_c / ( 1 + r_v + r_c ) } ].
    !
    ! This version is written with r_t = r_v + r_c, such that:
    !
    ! theta_v = theta * [ 1 + { (R_v/R_d) - 1 }
    !                         * { ( r_t - r_c ) / ( 1 + r_t ) }
    !                     - { r_c / ( 1 + r_t ) } ].
    !
    ! To use theta_l instead of theta, simply substitute the following for
    ! theta in the above expression:  theta_l + {L_v/(C_p*exner)} * r_c.
    !
    ! The CLUBB code uses a linearized version of the above equation (in order
    ! to calculate thv' terms -- such as w'thv', etc.) for theta_v throughout
    ! the model code, such that:
    !
    ! theta_v = theta_l + { (R_v/R_d) - 1 } * thv_ds * r_t
    !                   + [ {L_v/(C_p*exner)} - (R_v/R_d) * thv_ds ] * r_c;
    !
    ! where thv_ds is used as a reference value to approximate theta_l.
    call calculate_thvm( gr%nzt, ngrdcol, &
                          thlm, rtm, rcm, exner, &
                          thm * ( one + ep2 * ( rtm - rcm ) )**kappa, &
                          thvm )

    ! Recompute more accurate initial exner function, pressure, and density
    ! using thvm, which includes the effects of water vapor and cloud water.
    call hydrostatic( ngrdcol, gr, thvm, p_sfc, & ! Intent(in)
                      p_in_Pa, p_in_Pa_zm,      &  ! Intent(out)
                      exner, exner_zm,          &  ! Intent(out)
                      rho, rho_zm               )  ! Intent(out)

    !#### Calculate dry, static base-state density for the anelastic ####
    !#### equation set.  Calculate dry pressure from total pressure, ####
    !#### rtm, and rcm, and calculate dry exner from dry pressure.   ####

    !!! Calculate dry density on thermodynamic levels

    do k = 1, gr%nzt
      do i = 1, ngrdcol

        ! Calculate dry pressure from total pressure and water vapor mixing
        ! ratio, such that:  p_d = p / [ 1 + (R_v/R_d)*r_v ].
        p_dry(i,k) = p_in_Pa(i,k) / ( one + ep2 * ( rtm(i,k) - rcm(i,k) ) )

        ! Calculate dry exner from dry pressure.
        exner_dry(i,k) = ( p_dry(i,k) / p0 )**kappa

      end do
    end do

    ! Calculate dry potential temperature, theta_d, which is defined as:
    !
    ! theta_d = T / exner_d;
    !
    ! where exner_d = ( p_d / p0 )^(R_d/C_p).
    !
    ! Also note that standard potential temperature, theta:
    !
    ! theta = T / exner;
    !
    ! where exner = ( p / p0 )^(R_d/C_p).
    !
    ! Therefore, since both equations can be written in terms of temperature,
    ! T:  theta_d * exner_d = theta * exner.  Thus, theta_d can be written as:
    ! theta_d = theta * ( exner / exner_d ).  Furthermore, exner can be
    ! written in terms of exner_d, such that:
    !
    ! exner = exner_d * [ 1 + (R_v/R_d)*r_v ]^(R_d/C_p).
    !
    ! Thus, the equation for theta_d becomes:
    !
    ! theta_d = theta * [ 1 + (R_v/R_d)*r_v ]^(R_d/C_p).
    !
    ! In other words, there is a given mass of air that has temperature, T,
    ! pressure, p, a dry component of pressure, p_d, density rho, and a dry
    ! component of density, rho_d.  Pressure (or total pressure, which is
    ! dry air pressure plus water vapor pressure) is used to determine total
    ! exner.  Dividing temperature by exner yields theta.  Likewise, dry
    ! pressure is used to determine dry exner.  Dividing temperature by dry
    ! exner yields dry theta, which differs by actual theta by a small
    ! amount, which is given by the equations above.
    do k = 1, gr%nzt
      do i = 1, ngrdcol
        th_dry(i,k) = thm(i,k) * ( one + ep2 * ( rtm(i,k) - rcm(i,k) ) )**kappa
      end do
    end do

    ! Compute dry density using dry pressure, dry exner, and theta_d.
    do k = 1, gr%nzt
      do i = 1, ngrdcol
        rho_dry(i,k) = p_dry(i,k) / ( Rd * th_dry(i,k) * exner_dry(i,k) )
      end do
    end do

    !!! Calculate dry density on momentum levels

    ! Dry pressure at momentum level k = 1 is the dry pressure at the surface.
    do i = 1, ngrdcol
      p_dry_zm(i,1) = pd_sfc(i)
    end do

    tmp(:,:) = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, rtm(:,:) - rcm(:,:) )

    do k = 2, gr%nzm
      do i = 1, ngrdcol

        ! Calculate dry pressure on momentum levels from total pressure (on
        ! momentum levels) and water vapor mixing ratio (interpolated to
        ! momentum levels), such that:  p_d = p / [ 1 + (R_v/R_d)*r_v ].
        p_dry_zm(i,k) = p_in_Pa_zm(i,k) / ( one + ep2 * max( tmp(i,k), zero_threshold ) )
      end do
    end do

    do k = 1, gr%nzm
      do i = 1, ngrdcol
        ! Calculate dry exner on momentum levels from dry pressure on momentum
        ! levels.
        exner_dry_zm(i,k) = ( p_dry_zm(i,k) / p0 )**kappa
      end do
    end do

    thm_zm(:,:) = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, thm(:,:) )

    ! Calculate theta_d on momentum levels by interpolating theta and water
    ! vapor mixing ratio to momentum levels.
    do k = 1, gr%nzm
      do i = 1, ngrdcol
        th_dry_zm(i,k) = thm_zm(i,k) * ( one + ep2 * max( tmp(i,k), zero_threshold ) )**kappa
      end do
    end do

    ! Compute dry density on momentum levels using dry pressure on momentum
    ! levels, dry exner on momentum levels, and theta_d interpolated to
    ! momentum levels.
    do k = 1, gr%nzm
      do i = 1, ngrdcol
        rho_dry_zm(i,k) = p_dry_zm(i,k) / ( Rd * th_dry_zm(i,k) * exner_dry_zm(i,k) )
      end do
    end do

    ! The values of rho_dry and rho_dry_zm that were just calculated are dry
    ! (they do not take into account water vapor or cloud water).  These are
    ! the values of dry, static, base-state density that are needed for the
    ! anelastic equation set.
    rho_ds_zt = rho_dry
    rho_ds_zm = rho_dry_zm

    ! Since theta_d does not include water in any phase, the value of dry,
    ! static, base-state theta is the same as the values of both dry, static,
    ! base-state theta_l and dry, static, base-state theta_v.  Thus, for use
    ! with the anelastic equation set, thv_ds = thl_ds = theta_ds.
    thv_ds_zt = th_dry
    thv_ds_zm = th_dry_zm

    ! At this point, the values of the dry, static, base-state variables
    ! rho_ds_zt, rho_ds_zm, thv_ds_zt, and thv_ds_zm have been calculated.

    ! The CLUBB code is set up to be an anelastic code.  If use of the
    ! Boussinesq approximation is desired, rather than the anelastic
    ! approximation, reset the values of rho_ds_zt and rho_ds_zm to 1.  Also,
    ! reset the values of thv_ds_zt and thv_ds_zm to reference temperature T0.
    if ( l_use_boussinesq ) then
      rho_ds_zt = one
      rho_ds_zm = one
      thv_ds_zt = T0
      thv_ds_zm = T0
    endif ! otherwise, the code remains anelastic.

    ! Set the values of inverse, dry, static, base-state density.
    invrs_rho_ds_zm = one / rho_ds_zm
    invrs_rho_ds_zt = one / rho_ds_zt

    return

  end subroutine initialize_clubb_variables

  !-----------------------------------------------------------------------
  subroutine restart_clubb &
              ( gr, iunit, hydromet_dim, hm_metadata, & ! In
                restart_path_case, time_restart, & ! In
                um, upwp, vm, vpwp, up2, vp2, rtm, & ! Inout
                wprtp, thlm, wpthlp, rtp2, rtp3, & ! Inout
                thlp2, thlp3, rtpthlp, wp2, wp3, & ! Inout
                p_in_Pa, exner, rcm, cloud_frac, & ! Inout
                wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, & ! Inout
                wp2rtp, wp2thlp, uprcp, vprcp, & ! Inout
                rc_coef_zm, wp4, wpup2, wpvp2, wp2up2, & ! Inout
                wp2vp2, ice_supersat_frac, & ! Inout
                wm_zt, rho, rho_zm, rho_ds_zm, & ! Inout
                rho_ds_zt, thv_ds_zm, thv_ds_zt, & ! Inout
                thlm_forcing, rtm_forcing, wprtp_forcing, & ! Inout
                wpthlp_forcing, rtp2_forcing, & ! Inout
                thlp2_forcing, rtpthlp_forcing, & ! Inout
                hydromet, hydrometp2, wphydrometp, & ! Inout
                Ncm, Nccnm, thvm, em, tau_zm, tau_zt, & ! Inout
                Kh_zt, Kh_zm, ug, vg, & ! Inout
                thlprcp, & ! Inout
                sigma_sqd_w, sigma_sqd_w_zt, radht, & ! Inout
                deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K, & ! Inout
                pdf_params, pdf_params_zm, & ! Inout
                rcm_mc, rvm_mc, thlm_mc, & ! Out
                wprtp_mc, wpthlp_mc, rtp2_mc, & ! Out
                thlp2_mc, rtpthlp_mc, & ! Out
                wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc ) ! Out

    ! Description:
    !   Execute the necessary steps for the initialization of the
    !   CLUBB model to a designated point in the submitted GrADS file.
    !-----------------------------------------------------------------------
    use inputfields,only: &
        input_type, &  ! Variable(s)
        l_input_um, l_input_vm, l_input_rtm, l_input_thlm, &
        l_input_wp2, l_input_wprtp, l_input_wpthlp, &
        l_input_wp3, l_input_rtp2, l_input_rtp3, l_input_thlp2, &
        l_input_thlp3, l_input_rtpthlp, l_input_upwp, l_input_vpwp, &
        l_input_ug, l_input_vg, l_input_rcm, &
        l_input_wm_zt, l_input_exner, l_input_em, &
        l_input_p, l_input_rho, l_input_rho_zm, &
        l_input_rho_ds_zm, l_input_rho_ds_zt, &
        l_input_thv_ds_zm, l_input_thv_ds_zt, &
        l_input_Lscale, l_input_Lscale_up, l_input_Lscale_down, &
        l_input_Kh_zt, l_input_Kh_zm, l_input_tau_zm, l_input_tau_zt, &
        l_input_wpthvp, l_input_wp2thvp, l_input_wp2up, l_input_rtpthvp, l_input_thlpthvp, &
        l_input_wp2rtp, l_input_wp2thlp, l_input_uprcp, l_input_vprcp, &
        l_input_rc_coef_zm, l_input_wp4, l_input_wpup2, l_input_wpvp2, &
        l_input_wp2up2, l_input_wp2vp2, l_input_iss_frac, &
        l_input_radht, &
        l_input_w_1, l_input_w_2, l_input_varnce_w_1, l_input_varnce_w_2, &
        l_input_rt_1, l_input_rt_2, l_input_varnce_rt_1, l_input_varnce_rt_2, &
        l_input_thl_1, l_input_thl_2, l_input_varnce_thl_1, &
        l_input_varnce_thl_2, l_input_mixt_frac, l_input_chi_1, l_input_chi_2, &
        l_input_stdev_chi_1, l_input_stdev_chi_2, l_input_rc_1, l_input_rc_2, &
        l_input_w_1_zm, l_input_w_2_zm, l_input_varnce_w_1_zm, &
        l_input_varnce_w_2_zm, l_input_mixt_frac_zm, &
        l_input_thvm, l_input_rrm, l_input_Nrm, &
        l_input_rsm, l_input_rim, l_input_rgm, &
        l_input_thlm_forcing, l_input_rtm_forcing, &
        l_input_up2, l_input_vp2, l_input_sigma_sqd_w, l_input_Ncm, &
        l_input_Nccnm, l_input_Nim, l_input_cloud_frac, &
        l_input_sigma_sqd_w_zt, l_input_veg_T_in_K, l_input_deep_soil_T_in_K, &
        l_input_sfc_soil_T_in_K, l_input_wprtp_forcing, &
        l_input_wpthlp_forcing, l_input_rtp2_forcing, l_input_thlp2_forcing, &
        l_input_rtpthlp_forcing, l_input_thlprcp, l_input_rcm_mc, &
        l_input_rvm_mc, l_input_thlm_mc, l_input_wprtp_mc, l_input_wpthlp_mc, &
        l_input_rtp2_mc, l_input_thlp2_mc, l_input_rtpthlp_mc, &
        stat_files

    use inputfields, only: &
        compute_timestep, & !------------------------------------ Procedure(s)
        stat_fields_reader, &
        set_filenames, &
        get_clubb_variable_interpolated


    use grid_class, only: &
        grid,   & !--------------------------------------------------- Type
        zt2zm_api !--------------------------------------------------- Procedure(s)

    use constants_clubb, only: fstderr !-------------------------- Variables(s)

    use pdf_parameter_module, only: pdf_parameter !--------------- Type(s)

    use clubb_precision, only: time_precision, core_rknd !-------- Variable(s)

    use soil_vegetation, only: &
        l_soil_veg !---------------------------------------------- Variable(s)

    use parameters_microphys, only : &
        microphys_scheme, &  !------------------------------------ Variable(s)
        l_predict_Nc, &
        l_ice_microphys, &
        l_graupel

    use corr_varnce_module, only: &
        hm_metadata_type

    implicit none

    ! Input Variables
    type(grid), intent(in) :: gr

    integer, intent(in) :: iunit

    character(len=*), intent(in) :: &
      restart_path_case     ! Path to GrADS data for restart

    integer, intent(in) :: &
      hydromet_dim

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    real(kind=time_precision), intent(in) :: &
      time_restart

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(gr%nzt), intent(inout) :: &
      um,                & ! eastward grid-mean wind component (thermo. levs.)  [m/s]
      vm,                & ! northward grid-mean wind component (thermo. levs.) [m/s]
      rtm,               & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      thlm,              & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      rtp3,              & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
      thlp3,             & ! th_l'^3 (thermodynamic levels)                 [K^3]
      wp3,               & ! w'^3 (thermodynamic levels)                    [m^3/s^3]
      p_in_Pa,           & ! Air pressure (thermodynamic levels)            [Pa]
      exner,             & ! Exner function (thermodynamic levels)          [-]
      rcm,               & ! cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
      cloud_frac,        & ! cloud fraction (thermodynamic levels)          [-]
      wp2thvp,           & ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]
      wp2up,             & ! < w'^2 u' > (thermodynamic levels)             [m^3/s^3]
      wp2rtp,            & ! w'^2 rt' (thermodynamic levels)      [m^2/s^2 kg/kg]
      wp2thlp,           & ! w'^2 thl' (thermodynamic levels)     [m^2/s^2 K]
      wpup2,             & ! w'u'^2 (thermodynamic levels)        [m^3/s^3]
      wpvp2,             & ! w'v'^2 (thermodynamic levels)        [m^3/s^3]
      ice_supersat_frac    ! ice cloud fraction (thermo. levels)  [-]

    real( kind = core_rknd ), dimension(gr%nzm), intent(inout) :: &
      upwp,       & ! u'w' (momentum levels)                         [m^2/s^2]
      vpwp,       & ! v'w' (momentum levels)                         [m^2/s^2]
      up2,        & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2,        & ! v'^2 (momentum levels)                         [m^2/s^2]
      wprtp,      & ! w' r_t' (momentum levels)                      [kg/kg m/s]
      wpthlp,     & ! w'th_l' (momentum levels)                      [(m/s) K]
      rtp2,       & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      thlp2,      & ! th_l'^2 (momentum levels)                      [K^2]
      rtpthlp,    & ! r_t'th_l' (momentum levels)                    [(kg/kg) K]
      wp2,        & ! w'^2 (momentum levels)                         [m^2/s^2]
      wpthvp,     & ! < w' th_v' > (momentum levels)                 [kg/kg K]
      rtpthvp,    & ! < r_t' th_v' > (momentum levels)               [kg/kg K]
      thlpthvp,   & ! < th_l' th_v' > (momentum levels)              [K^2]
      uprcp,      & ! < u' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      vprcp,      & ! < v' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      rc_coef_zm, & ! Coef of X'r_c' in Eq. (34) (m-levs.) [K/(kg/kg)]
      wp4,        & ! w'^4 (momentum levels)               [m^4/s^4]
      wp2up2,     & ! w'^2 u'^2 (momentum levels)          [m^4/s^4]
      wp2vp2        ! w'^2 v'^2 (momentum levels)          [m^4/s^4]

    real( kind = core_rknd ), dimension(gr%nzt), intent(inout) :: &
      wm_zt,        & ! vertical mean wind component on thermo. levels  [m/s]
      rho,          & ! Air density on thermodynamic levels             [kg/m^3]
      rho_ds_zt,    & ! Dry, static density on thermo. levels           [kg/m^3]
      thv_ds_zt,    & ! Dry, base-state theta_v on thermo levels        [K]
      thlm_forcing, & ! liquid potential temp. forcing (thermo. levels) [K/s]
      rtm_forcing     ! total water forcing (thermo. levels)      [(kg/kg)/s]

    real( kind = core_rknd ), dimension(gr%nzm), intent(inout) :: &
      rho_zm,          & ! Air density on momentum levels               [kg/m^3]
      rho_ds_zm,       & ! Dry, static density on momentum levels       [kg/m^3]
      thv_ds_zm,       & ! Dry, base-state theta_v on momentum levels   [K]
      wprtp_forcing,   & ! total water turbulent flux forcing (m-levs) [m*K/s^2]
      wpthlp_forcing,  & ! liq pot temp turb flux forcing (m-levs)[m(kg/kg)/s^2]
      rtp2_forcing,    & ! total water variance forcing (m-levs)   [(kg/kg)^2/s]
      thlp2_forcing,   & ! liq pot temp variance forcing (m-levs)  [K^2/s]
      rtpthlp_forcing    ! <r_t'th_l'> covariance forcing (m-levs) [K*(kg/kg)/s]

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(inout) :: &
      hydromet       ! Array of hydrometeors                [hm units]

    real( kind = core_rknd ), dimension(gr%nzm,hydromet_dim), intent(inout) :: &
      hydrometp2,  & ! Variance of a hydrometeor (m-levs.)  [<hm units>^2]
      wphydrometp    ! Covariance of w and a hydrometeor    [(m/s) <hm units>]

    real( kind = core_rknd ), dimension(gr%nzt), intent(inout) :: &
      Ncm,            & ! Mean cloud droplet concentration, <N_c> (t-levs.)    [num/kg]
      Nccnm,          & ! Cloud condensation nuclei concentration (COAMPS/MG)  [num/kg]
      thvm,           & ! Virtual potential temperature                        [K]
      tau_zt,         & ! Eddy dissipation time scale on thermodynamic levels  [s]
      Kh_zt,          & ! Eddy diffusivity coefficient on thermodynamic levels [m^2/s]
      ug,             & ! u geostrophic wind                                   [m/s]
      vg,             & ! v geostrophic wind                                   [m/s]
      sigma_sqd_w_zt, & ! PDF width parameter interpolated to t-levs.          [-]
      radht             ! SW + LW heating rate                                 [K/s]

    real( kind = core_rknd ), dimension(gr%nzm), intent(inout) :: &
      em,          & ! Turbulent Kinetic Energy (TKE)                       [m^2/s^2]
      tau_zm,      & ! Eddy dissipation time scale on momentum levels       [s]
      Kh_zm,       & ! Eddy diffusivity coefficient on momentum levels      [m^2/s]
      thlprcp,     & ! thl'rc'                                              [K kg/kg]
      sigma_sqd_w    ! PDF width parameter (momentum levels)                [-]

    real( kind = core_rknd ), intent(inout) :: &
      deep_soil_T_in_K, &
      sfc_soil_T_in_K, &
      veg_T_in_K

    type(pdf_parameter), intent(inout) :: &
      pdf_params,    & ! PDF parameters (thermodynamic levels)    [units vary]
      pdf_params_zm    ! PDF parameters on momentum levels        [units vary]

    ! Output
    real( kind = core_rknd ), intent(out) :: &
      wpthlp_sfc,      & ! w'theta_l' surface flux   [(m K)/s]
      wprtp_sfc,       & ! w'rt' surface flux        [(m kg)/(kg s)]
      upwp_sfc,        & ! u'w' at surface           [m^2/s^2]
      vpwp_sfc           ! v'w' at surface           [m^2/s^2]

    real( kind = core_rknd ), dimension(gr%nzt), intent(out) :: &
      rcm_mc, &    ! Tendency of liquid water due to microphysics      [kg/kg/s]
      rvm_mc, &    ! Tendency of vapor water due to microphysics       [kg/kg/s]
      thlm_mc      ! Tendency of liquid pot. temp. due to microphysics [K/s]

    real( kind = core_rknd ), dimension(gr%nzm), intent(out) :: &
      wprtp_mc, &  ! Microphysics tendency for <w'rt'>   [m*(kg/kg)/s^2]
      wpthlp_mc, & ! Microphysics tendency for <w'thl'>  [m*K/s^2]
      rtp2_mc, &   ! Microphysics tendency for <rt'^2>   [(kg/kg)^2/s]
      thlp2_mc, &  ! Microphysics tendency for <thl'^2>  [K^2/s]
      rtpthlp_mc   ! Microphysics tendency for <rt'thl'> [K*(kg/kg)/s]

    ! Local variables
    integer :: timestep

    logical :: l_restart, l_read_error, l_fatal_error


    ! --- Begin Code ---

    l_fatal_error = .false.

    ! Inform inputfields module
    input_type = "clubb"
    l_input_um   = .true.
    l_input_vm   = .true.
    l_input_rtm  = .true.
    l_input_thlm = .true.
    l_input_wp2  = .true.
    l_input_ug   = .true.
    l_input_vg   = .true.
    l_input_rcm  = .true.
    l_input_wm_zt  = .true.
    l_input_exner = .true.
    l_input_em = .true.
    l_input_p = .true.
    l_input_rho = .true.
    l_input_rho_zm = .true.
    l_input_rho_ds_zm = .true.
    l_input_rho_ds_zt = .true.
    l_input_thv_ds_zm = .true.
    l_input_thv_ds_zt = .true.
    l_input_Lscale = .true.
    l_input_Lscale_up = .true.
    l_input_Lscale_down = .true.
    l_input_Kh_zt = .true.
    l_input_Kh_zm = .true.
    l_input_tau_zm = .true.
    l_input_tau_zt = .true.
    l_input_thvm = .true.
    l_input_wpthvp = .true.
    l_input_wp2thvp = .true.
    l_input_wp2up = .true.
    l_input_rtpthvp = .true.
    l_input_thlpthvp = .true.
    l_input_wp2rtp = .true.
    l_input_wp2thlp = .true.
    l_input_uprcp = .true.
    l_input_vprcp = .true.
    l_input_rc_coef_zm = .true.
    l_input_wp4 = .true.
    l_input_wpup2 = .true.
    l_input_wpvp2 = .true.
    l_input_wp2up2 = .true.
    l_input_wp2vp2 = .true.
    l_input_iss_frac = .true.
    l_input_w_1 = .true.
    l_input_w_2 = .true.
    l_input_varnce_w_1 = .true.
    l_input_varnce_w_2 = .true.
    l_input_rt_1 = .true.
    l_input_rt_2 = .true.
    l_input_varnce_rt_1 = .true.
    l_input_varnce_rt_2 = .true.
    l_input_thl_1 = .true.
    l_input_thl_2 = .true.
    l_input_varnce_thl_1 = .true.
    l_input_varnce_thl_2 = .true.
    l_input_mixt_frac = .true.
    l_input_chi_1   = .true.
    l_input_chi_2   = .true.
    l_input_stdev_chi_1  = .true.
    l_input_stdev_chi_2  = .true.
    l_input_rc_1  = .true.
    l_input_rc_2  = .true.
    l_input_w_1_zm = .true.
    l_input_w_2_zm = .true.
    l_input_varnce_w_1_zm = .true.
    l_input_varnce_w_2_zm = .true.
    l_input_mixt_frac_zm = .true.
    l_input_radht = .true.

    select case ( trim( microphys_scheme ) )
    case ( "coamps" )
      l_input_rrm = .true.
      l_input_rsm = .true.
      l_input_rim = .true.
      l_input_rgm = .true.
      l_input_Nccnm = .true.
      l_input_Ncm = .true.
      l_input_Nrm = .true.
      l_input_Nim =  .true.

    case ( "morrison" )
      l_input_rrm = .true.
      l_input_Nrm = .true.
      if ( l_ice_microphys ) then
        l_input_rsm = .true.
        l_input_rim = .true.
        l_input_Nim = .true.
        if ( l_graupel ) then
          l_input_rgm = .true.
        else
          l_input_rgm = .false.
        end if
      else
        l_input_rsm = .false.
        l_input_rim = .false.
        l_input_Nim = .false.
        l_input_rgm = .false.
      end if
      l_input_Nccnm = .false.
      if ( l_predict_Nc ) then
        l_input_Ncm = .true.
      else
        l_input_Ncm = .false.
      end if

    case ( "khairoutdinov_kogan" )
      l_input_rrm = .true.
      l_input_rsm = .false.
      l_input_rim = .false.
      l_input_rgm = .false.
      l_input_Nccnm = .false.
      l_input_Ncm = .false.
      l_input_Nrm = .true.
      l_input_Nim = .false.

    case default
      l_input_rrm = .false.
      l_input_rsm = .false.
      l_input_rim = .false.
      l_input_rgm = .false.
      l_input_Nccnm = .false.
      l_input_Ncm = .false.
      l_input_Nrm = .false.
      l_input_Nim = .false.

    end select

    if ( l_soil_veg ) then
      l_input_veg_T_in_K = .true.
      l_input_deep_soil_T_in_K = .true.
      l_input_sfc_soil_T_in_K = .true.
    end if

    l_input_wprtp = .true.
    l_input_wpthlp = .true.
    l_input_wp3 = .true.
    l_input_rtp2 = .true.
    l_input_rtp3 = .true.
    l_input_thlp2 = .true.
    l_input_thlp3 = .true.
    l_input_rtpthlp = .true.
    l_input_upwp = .true.
    l_input_vpwp = .true.
    l_input_thlm_forcing = .true.
    l_input_rtm_forcing = .true.
    l_input_up2 = .true.
    l_input_vp2 = .true.
    l_input_sigma_sqd_w = .true.
    l_input_cloud_frac  = .true.
    l_input_sigma_sqd_w_zt = .true.
    l_input_wprtp_forcing = .true.
    l_input_wpthlp_forcing = .true.
    l_input_rtp2_forcing = .true.
    l_input_thlp2_forcing = .true.
    l_input_rtpthlp_forcing = .true.
    l_input_thlprcp = .true.
    l_input_rcm_mc = .true.
    l_input_rvm_mc = .true.
    l_input_thlm_mc = .true.
    l_input_wprtp_mc = .true.
    l_input_wpthlp_mc = .true.
    l_input_rtp2_mc = .true.
    l_input_thlp2_mc = .true.
    l_input_rtpthlp_mc = .true.

    call set_filenames( "../"//trim( restart_path_case ) )
    ! Determine the nearest timestep in the GRADS file to the
    ! restart time.
    l_restart = .true.

    call compute_timestep &
      ( iunit, stat_files(1), l_restart, time_restart, &    ! Intent(in)
        timestep )                                          ! Intent(out)

    ! Sanity check for input time_restart
    if ( timestep < 0 ) then
      write(fstderr,*) "Invalid time_restart"
      error stop
    end if

    ! Read data from stats files
    call stat_fields_reader( gr, timestep, hydromet_dim, hm_metadata, & ! In
                              um, upwp, vm, vpwp, & ! Inout
                              up2, vp2, rtm, & ! Inout
                              wprtp, thlm, wpthlp, & ! Inout
                              rtp2, rtp3, & ! Inout
                              thlp2, thlp3, rtpthlp, & ! Inout
                              wp2, wp3, & ! Inout
                              p_in_Pa, exner, rcm, cloud_frac, & ! Inout
                              wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, & ! Inout
                              wp2rtp, wp2thlp, uprcp, vprcp, & ! Inout
                              rc_coef_zm, wp4, wpup2, & ! Inout
                              wpvp2, wp2up2, & ! Inout
                              wp2vp2, ice_supersat_frac, & ! Inout
                              wm_zt, rho, rho_zm, rho_ds_zm, & ! Inout
                              rho_ds_zt, thv_ds_zm, thv_ds_zt, & ! Inout
                              thlm_forcing, rtm_forcing, wprtp_forcing, & ! Inout
                              wpthlp_forcing, rtp2_forcing, & ! Inout
                              thlp2_forcing, rtpthlp_forcing, & ! Inout
                              hydromet, hydrometp2, wphydrometp, & ! Inout
                              Ncm, Nccnm, thvm, em, & ! Inout
                              tau_zm, tau_zt, & ! Inout
                              Kh_zt, Kh_zm, ug, vg, & ! Inout
                              thlprcp, & ! Inout
                              sigma_sqd_w, sigma_sqd_w_zt, radht, & ! Inout
                              deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K, & ! Inout
                              pdf_params, pdf_params_zm ) ! Inout

      call get_clubb_variable_interpolated &
            ( l_input_rcm_mc, stat_files(1), "rcm_mc", gr%nzt, timestep, &
              gr%zt(1,:), rcm_mc, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
            ( l_input_rvm_mc, stat_files(1), "rvm_mc", gr%nzt, timestep, &
              gr%zt(1,:), rvm_mc, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
            ( l_input_thlm_mc, stat_files(1), "thlm_mc", gr%nzt, timestep, &
              gr%zt(1,:), thlm_mc, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
            ( l_input_wprtp_mc, stat_files(2), "wprtp_mc", gr%nzm, timestep, &
              gr%zm(1,:), wprtp_mc, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
            ( l_input_wpthlp_mc, stat_files(2), "wpthlp_mc", gr%nzm, timestep, &
              gr%zm(1,:), wpthlp_mc, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
            ( l_input_rtp2_mc, stat_files(2), "rtp2_mc", gr%nzm, timestep, &
              gr%zm(1,:), rtp2_mc, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
            ( l_input_thlp2_mc, stat_files(2), "thlp2_mc", gr%nzm, timestep, &
              gr%zm(1,:), thlp2_mc, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
            ( l_input_rtpthlp_mc, stat_files(2), "rtpthlp_mc", gr%nzm, timestep, &
              gr%zm(1,:), rtpthlp_mc, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error


    wpthlp_sfc = wpthlp(1)
    wprtp_sfc  = wprtp(1)
    upwp_sfc   = upwp(1)
    vpwp_sfc   = vpwp(1)

    return
  end subroutine restart_clubb

end module clubb_driver
