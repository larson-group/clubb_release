!-----------------------------------------------------------------------
! $Id$

module clubb_driver

! Description:
!   Contains the necessary subroutines to execute individual CLUBB
!   model runs, using one of the driver programs (the simplest case
!   being the clubb_standalone program).
!-----------------------------------------------------------------------

  use clubb_precision, only: time_precision ! Variable(s)
  use text_writer, only: write_text

  use grid_class, only: grid ! Type

  use stats_type, only: stats ! Type

  use mt95, only: &
    genrand_intg

  implicit none

  ! Variables that contains all the statistics

  type (stats), target, public, save :: stats_zt,      & ! stats_zt grid
                                        stats_zm,      & ! stats_zm grid
                                        stats_lh_zt,   & ! stats_lh_zt grid
                                        stats_lh_sfc,  & ! stats_lh_sfc grid
                                        stats_rad_zt,  & ! stats_rad_zt grid
                                        stats_rad_zm,  & ! stats_rad_zm grid
                                        stats_sfc        ! stats_sfc

!$omp threadprivate(stats_zt, stats_zm, stats_lh_zt, stats_lh_sfc)
!$omp threadprivate(stats_rad_zt, stats_rad_zm, stats_sfc)

  type(grid), target :: gr
!$omp threadprivate(gr)

  ! Setup run_clubb() as the sole external interface
  private ::  &
    initialize_clubb, &
    initialize_clubb_variables, &
    prescribe_forcings, &
    restart_clubb

  public :: &
    run_clubb


  private ! Default to private

  contains

  !-----------------------------------------------------------------------
  subroutine run_clubb & 
             ( params, runfile, l_stdout, &
               model_flags_array )
    ! Description:
    !   Subprogram to integrate the partial differential equations for pdf
    !   closure.

    ! References:
    !   None
    !---------------------------------------------------------------------

    use grid_class, only: read_grid_heights, zt2zm, zm2zt !------------------ Procedure(s)

    use parameter_indices, only: &
      nparams, ic_K, iSkw_denom_coef, iSkw_max_mag !------------------------- Variable(s)

    use numerical_check, only: invalid_model_arrays !------------------------ Procedure(s)

    use inputfields, only: &
      inputfields_init, compute_timestep, stat_fields_reader, & !------------ Procedure(s)
      l_input_wp3                                                 !---------- Variable(s)

    use inputfields, only: stat_files

    use parameters_tunable, only: &
        params_list, & !--------------------------------- Variable(s)
        nu_vertical_res_dep !---------------------------- Type(s)

    use advance_clubb_core_module, only: &
      setup_clubb_core,  & !------------------------------------------------- Procedure(s)
      advance_clubb_core, &
      calculate_thlp2_rad

    use constants_clubb, only: &
        fstdout, & !--------------------------------------------------------- Constant(s)
        fstderr, &
        three_halves, &
        one, &
        one_half, &
        zero, &
        rt_tol, &
        thl_tol, &
        w_tol, &
        w_tol_sqd, &
        em_min, &
        eps
        
    use clubb_api_module, only: &
      setup_pdf_parameters_api, &
      precipitation_fractions, &
      init_precip_fracs_api

    use pdf_parameter_module, only: &
        pdf_parameter,                 & !----------------------------------- Variable Type(s)
        implicit_coefs_terms,          &
        init_pdf_params,               & !----------------------------------- Procedure(s)
        init_pdf_implicit_coefs_terms

    use error_code, only: &
        clubb_at_least_debug_level,  & ! ------------------------------------ Procedures
        set_clubb_debug_level,       &
        err_code,                    & ! ------------------------------------ Error Indicator
        clubb_fatal_error              ! ------------------------------------ Constant

    use clubb_precision, only: time_precision, core_rknd !------------------- Constants

    use array_index, only: iisclr_rt, iisclr_thl, iisclr_CO2, & !------------ Variables
      iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2, &
      iiri, iirs, iirg

    use microphys_driver, only: &
        calc_microphys_scheme_tendcies  !------------------------------------ Procedure(s)

    use advance_microphys_module, only: &
        advance_microphys  !------------------------------------------------- Procedure(s)

    use microphys_init_cleanup, only: &
        init_microphys     !------------------------------------------------- Procedure

    use bugsrad_driver, only: init_radiation !------------------------------- Subroutine

    use model_flags, only: &
      set_default_clubb_config_flags, & !-------------------------------- Procedure(s)
      initialize_clubb_config_flags_type, &
      print_clubb_config_flags, &
      clubb_config_flags_type, & !--------------------------------------- Type(s)
      l_silhs_rad, & !--------------------------------------------------- Constants
      l_pos_def, l_hole_fill, &
      l_gamma_Skw, l_byteswap_io, &
      l_quintic_poly_interp !-------------------------------------------- Variable(s)

    use stats_variables, only: l_stats_last, l_stats_samp, & !--------------- Variable(s)
      l_output_rad_files

    use stats_variables, only: &
        irtm_mc,     & !----------------------------------------------------- Variables
        irvm_mc,     &
        ircm_mc,     &
        ithlm_mc,    &
        iwprtp_mc,   &
        iwpthlp_mc,  &
        irtp2_mc,    &
        ithlp2_mc,   &
        irtpthlp_mc, &
        iFrad

    use stats_variables, only: &
        l_allow_small_stats_tout

    use stats_clubb_utilities, only:  & 
        stats_begin_timestep, stats_end_timestep,  & !----------------------- Procedure(s)
        stats_init

    use stats_type_utilities, only: &
        stat_update_var !---------------------------------------------------- Procedure

    use sounding, only: sclr_max !------------------------------------------- Variable(s)

    use time_dependent_input, only: &
        l_t_dependent,    & !------------------------------------------------ Variable(s)
        l_input_xpwp_sfc, &
        l_ignore_forcings

    use sponge_layer_damping, only: &
        thlm_sponge_damp_settings,    & !------------------------------------ Variable(s)
        rtm_sponge_damp_settings,     &
        uv_sponge_damp_settings,      &
        wp2_sponge_damp_settings,     &
        wp3_sponge_damp_settings,     &
        up2_vp2_sponge_damp_settings

    use extended_atmosphere_module, only: &
        total_atmos_dim, & !------------------------------------------------- Variable(s)
        complete_alt, &
        complete_momentum

    use parameters_radiation, only: rad_scheme !----------------------------- Variable(s)

    use Skx_module, only: Skx_func !----------------------------------------- Procedure(s)

    use calc_pressure, only: calculate_thvm !-------------------------------- Procedure(s)

    use clip_explicit, only: clip_skewness_core !---------------------------- Procedure(s)

#ifdef SILHS
    use parameters_microphys, only: &
      lh_microphys_type,     & !--------------------------------------------- Variable(s)
      lh_microphys_disabled, &
      lh_seed,               &
      microphys_scheme,      &
      lh_num_samples,    &
      lh_sequence_length
      
    use silhs_api_module, only: &
      generate_silhs_sample_api, & !----------------------------------------- Procedure(s)
      clip_transform_silhs_output_api, &
      latin_hypercube_2D_output_api

    use latin_hypercube_driver_module, only: &
      stats_accumulate_lh

    use latin_hypercube_arrays, only: &
      cleanup_latin_hypercube_arrays !-------------------------------------- Procedure(s)

    use simple_rad_module, only: simple_rad_lba_init !---------------------- Procedure(s)

#endif

    use variables_radiation_module, only: &
      setup_radiation_variables    !---------------------------------------- Procedure(s)

    use text_writer, only: write_text, write_date !------------------------- Procedure(s)

    use clubb_model_settings, only: &
      initialize_clubb_model_settings !------------------------------------- Procedure(s)

    use clubb_model_settings, only: &
      time_initial, & !----------------------------------------------------- Variable(s)
      time_final, &
      time_current, &
      nzmax, &
      grid_type, &
      extended_atmos_range_size, &
      lin_int_buffer, &
      zt_grid_fname, &
      zm_grid_fname, &
      zm_top, &
      zm_init, &
      deltaz, &
      day, month, year, &
      lat_vals, &
      lon_vals, &
      sfc_elevation, &
      runtype, &
      sfctype, &
      dt_rad, &
      dt_main

    use soil_vegetation, only: &
        l_soil_veg !------------------------------------------------------ Variable(s)

    use soil_vegetation, only: &
        initialize_soil_veg !--------------------------------------------- Procedure(s)

    use parameters_model, only: &
        rtm_min, &
        rtm_nudge_max_altitude

    use corr_varnce_module, only: &
        corr_array_n_cloud, & !------------------------------------------- Variable(s)
        corr_array_n_below, &
        pdf_dim

    use setup_clubb_pdf_params, only: &
        setup_pdf_parameters    !----------------------------------------- Procedure(s)

    use mixed_moment_PDF_integrals, only: &
        hydrometeor_mixed_moments    !------------------------------------ Procedure(s)

    use hydromet_pdf_parameter_module, only: &
        hydromet_pdf_parameter,   & !------------------------------------- Type(s)
        precipitation_fractions,  &
        init_hydromet_pdf_params, & !------------------------------------- Procedure(s)
        init_precip_fracs

    use fill_holes, only: &
        vertical_avg  !--------------------------------------------------- Procedure(s)

    use parameters_silhs, only: &
        silhs_config_flags_type !----------------------------------------- Type(s)

#ifdef _OPENMP
    ! Because Fortran I/O is not thread safe, we use this here to
    ! ensure that no model uses the same file number simultaneously
    ! when doing a tuning run. -dschanen 31 Jan 2007
    use omp_lib, only: &
      omp_get_thread_num ! Function
#endif

    implicit none

    ! External
    intrinsic :: mod, real, int, trim, floor, max, sqrt

    ! Constant Parameters
    logical, parameter :: &
      l_host_applies_sfc_fluxes = .false., &
      l_implemented = .false.

    logical, parameter :: &
      l_write_to_file = .true. ! If true, will write case information to a file

    integer, parameter :: nlon = 1, nlat = 1 ! Number of points in the X/Y [-]

    ! Input Variables
    logical, intent(in) ::  & 
      l_stdout   ! Whether to print output per timestep

    ! Subroutine Arguments (Model Setting)
    character(len=*), intent(in) ::  & 
      runfile ! Name of file containing &model_setting and &sounding

    logical, optional, dimension(:), intent(in) :: &
      model_flags_array ! Array containing model flags (for the clubb_tuner only)

    ! Local Variables
    ! Internal Timing Variables
    integer :: & 
      ifinal

    integer :: & 
      debug_level     ! Amount of debugging information

    real( kind = core_rknd ) ::  & 
      fcor,            & ! Coriolis parameter            [s^-1]
      T0,              & ! Reference Temperature         [K]
      ts_nudge           ! Timescale for u/v nudging     [s]

    real( kind = core_rknd ), dimension(sclr_max) :: & 
      sclr_tol        ! Thresholds on the passive scalars     [units vary]

    real(kind=time_precision) :: & 
      time_restart    ! Time of model restart run     [s]

    logical :: &
      l_restart,      & ! Flag for restarting from GrADS file
      l_input_fields    ! Whether to set model variables from a file

    logical :: l_use_Ncn_to_Nc ! Whether to call Ncn_to_Nc (.true.) or not (.false.);
                               ! Ncn_to_Nc might cause problems with the MG microphysics 
                               ! since the changes made here (Nc-tendency) are not fed into 
                               ! the microphysics

    character(len=6) :: &
      saturation_formula ! "bolton" approx. or "flatau" approx.

    character(len=128) ::  & 
      restart_path_case,   & ! GrADS file used in case of restart
      forcings_file_path     ! Path to the forcing files

    logical :: & 
      l_stats ! Whether statistics are computed and output to disk

    character(len=10) :: & 
      stats_fmt  ! File format for stats; typically GrADS.

    character(len=100) :: & 
      fname_prefix, & ! Prefix of stats filenames, to be followed by, for example "_zt"
      fdir            ! Output directory

    real( kind = core_rknd ) :: & 
      stats_tsamp,   & ! Stats sampling interval [s]
      stats_tout       ! Stats output interval   [s]

    integer :: &
      stats_nsamp,   & ! Stats sampling interval [timestep]
      stats_nout       ! Stats output interval   [timestep]

    ! Grid altitude arrays
    real( kind = core_rknd ), dimension(:), allocatable ::  & 
      momentum_heights, thermodynamic_heights ! [m]

    ! Dummy dx and dy horizontal grid spacing.
    real( kind = core_rknd ) :: dummy_dx, dummy_dy  ! [m]

    integer :: &
      itime, i, j, & ! Local Loop Variables
      iinit          ! initial iteration

    integer ::  & 
      iunit,           & ! File unit used for I/O
      hydromet_dim,    & ! Number of hydrometeor species        [#]
      sclr_dim,        & ! Number of passive scalars            [#]
      edsclr_dim         ! Number of passive scalars            [#]

    integer :: itime_nearest ! Used for and inputfields run [s]

    character(len=150) :: &
      case_info_file ! The filename for case info

    real( kind = core_rknd ), dimension(0) :: rad_dummy ! Dummy variable for radiation levels

    real( kind = core_rknd ), dimension(:), allocatable ::  &
      um,      & ! eastward grid-mean wind component (thermo. levs.)  [m/s]
      upwp,    & ! u'w' (momentum levels)                         [m^2/s^2]
      vm,      & ! northward grid-mean wind component (thermo. levs.) [m/s]
      vpwp,    & ! v'w' (momentum levels)                         [m^2/s^2]
      up2,     & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2,     & ! v'^2 (momentum levels)                         [m^2/s^2]
      up3,     & ! u'^3 (thermodynamic levels)                    [m^3/s^3]
      vp3,     & ! v'^3 (thermodynamic levels)                    [m^3/s^3]
      rtm,     & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      wprtp,   & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
      thlm,    & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      wpthlp,  & ! w'th_l' (momentum levels)                      [(m/s) K]
      rtp2,    & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      rtp3,    & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
      thlp2,   & ! th_l'^2 (momentum levels)                      [K^2]
      thlp3,   & ! th_l'^3 (thermodynamic levels)                 [K^3]
      rtpthlp, & ! r_t'th_l' (momentum levels)                    [(kg/kg) K]
      wp2,     & ! w'^2 (momentum levels)                         [m^2/s^2]
      wp3        ! w'^3 (thermodynamic levels)                    [m^3/s^3]

    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      rcm,      & ! cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
      delta_zm
    
    real( kind = core_rknd ), dimension(:), allocatable :: &
      p_in_Pa,    & ! Air pressure (thermodynamic levels)       [Pa]
      exner,      & ! Exner function (thermodynamic levels)     [-]
      cloud_frac, & ! cloud fraction (thermodynamic levels)     [-]
      wpthvp,     & ! < w' th_v' > (momentum levels)            [kg/kg K]
      wp2thvp,    & ! < w'^2 th_v' > (thermodynamic levels)     [m^2/s^2 K]
      rtpthvp,    & ! < r_t' th_v' > (momentum levels)          [kg/kg K]
      thlpthvp,   & ! < th_l' th_v' > (momentum levels)         [K^2]
      uprcp,      & ! < u' r_c' >                               [(m kg)/(s kg)]
      vprcp         ! < v' r_c' >                               [(m kg)/(s kg)]

    real( kind = core_rknd ), dimension(:,:), allocatable ::  &
      rho_ds_zt  ! Dry, static density on thermo. levels      [kg/m^3]
      
    
    real( kind = core_rknd ), dimension(:), allocatable ::  &
      wm_zm,           & ! vertical mean wind comp. on momentum levs  [m/s]
      wm_zt,           & ! vertical mean wind comp. on thermo. levs   [m/s]
      rho,             & ! Air density on thermodynamic levels        [kg/m^3]
      rho_zm,          & ! Air density on momentum levels             [kg/m^3]
      rho_ds_zm,       & ! Dry, static density on momentum levels     [kg/m^3]
      invrs_rho_ds_zm, & ! Inverse dry, static density on m-levs.     [m^3/kg]
      invrs_rho_ds_zt, & ! Inverse dry, static density on thermo levs.[m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v on momentum levs.  [K]
      thv_ds_zt          ! Dry, base-state theta_v on thermo levs.    [K]

    real( kind = core_rknd ), dimension(:), allocatable ::  &
      thlm_forcing,    & ! liq. wat. pot. temp. forcing (thermo. levs)[K/s]
      rtm_forcing,     & ! total water forcing (thermo. levels)       [(kg/kg)/s]
      um_forcing,      & ! eastward wind forcing (thermo. levels)     [m/s/s]
      vm_forcing,      & ! northward wind forcing (thermo. levels)    [m/s/s]
      wprtp_forcing,   & ! total water turbulent flux forcing (m-levs)[m*K/s^2]
      wpthlp_forcing,  & ! liq pot temp turb flux forcing (m-levs)    [m(kg/kg)/s^2]
      rtp2_forcing,    & ! total water variance forcing (m-levs)      [(kg/kg)^2/s]
      thlp2_forcing,   & ! liq pot temp variance forcing (m-levs)     [K^2/s]
      rtpthlp_forcing    ! <r_t'th_l'> covariance forcing (m-levs)    [K(kg/kg)/s]

    type(pdf_parameter), allocatable :: &
      pdf_params ! PDF parameters (thermodynamic levels)    [units vary]
      
    type(pdf_parameter), allocatable :: &
      pdf_params_zm    ! PDF parameters on momentum levels        [units vary]

    type(implicit_coefs_terms), allocatable :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      hydromet,    & ! Array of hydrometeors                [hm units]
      hydrometp2,  & ! Variance of a hydrometeor (m-levs.)  [<hm units>^2]
      wphydrometp    ! Covariance of w and a hydrometeor    [(m/s) <hm units>]

    real( kind = core_rknd ), dimension(:), allocatable :: &
      Ncm,    & ! Mean cloud droplet concentration, <N_c> (t-levs.)    [num/kg]
      Nccnm,  & ! Cloud condensation nuclei concentration (COAMPS/MG)  [num/kg]
      thvm,   & ! Virtual potential temperature                        [K]
      em,     & ! Turbulent Kinetic Energy (TKE)                       [m^2/s^2]
      tau_zm, & ! Eddy dissipation time scale on momentum levels       [s]
      tau_zt, & ! Eddy dissipation time scale on thermodynamic levels  [s]
      Kh_zt,  & ! Eddy diffusivity coefficient on thermodynamic levels [m^2/s]
      Kh_zm     ! Eddy diffusivity coefficient on momentum levels      [m^2/s]

    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      Lscale   ! Length scale                                 [m]
      
    real( kind = core_rknd ), dimension(:), allocatable :: &
      Lscale_up,      & ! Length scale (upwards component)             [m]
      Lscale_down,    & ! Length scale (downwards component)           [m]
      thlprcp,        & ! thl'rc'                                      [K kg/kg]
      sigma_sqd_w,    & ! PDF width parameter (momentum levels)        [-]
      sigma_sqd_w_zt    ! PDF width parameter interpolated to t-levs.  [-]

    real( kind = core_rknd ), dimension(:), allocatable ::  &
      wprcp,             & ! w'r_c' (momentum levels)              [(kg/kg) m/s]
      w_up_in_cloud,     & ! Average upward velocity within liquid cloud   [m/s]
      ice_supersat_frac, & ! ice cloud fraction (thermo. levels)   [-]
      rcm_in_layer,      & ! rcm within cloud layer                [kg/kg]
      cloud_cover,       & ! cloud cover                           [-]
      invrs_tau_zm         ! One divided by tau on zm levels       [1/s]

    real( kind = core_rknd ), dimension(:), allocatable :: &
      ug,       & ! u geostrophic wind                           [m/s]
      vg,       & ! v geostrophic wind                           [m/s]
      rtm_ref,  & ! Initial total water mixing ratio             [kg/kg]
      thlm_ref, & ! Initial liquid water potential temperature   [K]
      um_ref,   & ! Initial u wind                               [m/s]
      vm_ref      ! Initial v wind                               [m/s]

    real( kind = core_rknd ) ::  &
      wpthlp_sfc, & ! w' theta_l' at surface   [(m K)/s]
      wprtp_sfc,  & ! w' r_t' at surface       [(kg m)/( kg s)]
      upwp_sfc,   & ! u'w' at surface          [m^2/s^2]
      vpwp_sfc      ! v'w' at surface          [m^2/s^2]

    real( kind = core_rknd ), allocatable, dimension(:) :: &
      rcm_mc, & ! Tendency of liquid water due to microphysics      [kg/kg/s]
      rvm_mc, & ! Tendency of vapor water due to microphysics       [kg/kg/s]
      thlm_mc   ! Tendency of liquid pot. temp. due to microphysics [K/s]

    real( kind = core_rknd ), allocatable, dimension(:) :: &
      wprtp_mc,   & ! Microphysics tendency for <w'rt'>   [m*(kg/kg)/s^2]
      wpthlp_mc,  & ! Microphysics tendency for <w'thl'>  [m*K/s^2]
      rtp2_mc,    & ! Microphysics tendency for <rt'^2>   [(kg/kg)^2/s]
      thlp2_mc,   & ! Microphysics tendency for <thl'^2>  [K^2/s]
      rtpthlp_mc    ! Microphysics tendency for <rt'thl'> [K*(kg/kg)/s]

    real( kind = core_rknd ), allocatable, dimension(:,:) :: &
      hydromet_mc     ! Microphysics tendency for mean hydrometeors  [units/s]

    real( kind = core_rknd ), allocatable, dimension(:) :: &
      Ncm_mc     ! Microphysics tendency for Ncm                     [num/kg/s]

    real( kind = core_rknd ), allocatable, dimension(:,:) :: &
      hydromet_vel_zt   ! Mean hydrometeor sed. velocity on thermo. levs. [m/s]

    real( kind = core_rknd ), allocatable, dimension(:,:) :: &
      hydromet_vel_covar_zt_impc, & ! Imp. comp. <V_xx'x_x'> t-levs [m/s]
      hydromet_vel_covar_zt_expc    ! Exp. comp. <V_xx'x_x'> t-levs [units(m/s)]

    real( kind = core_rknd ), allocatable, dimension(:) :: &
      rfrzm    ! Total ice-phase water mixing ratio        [kg/kg]

    real( kind = core_rknd ), allocatable, dimension(:) :: &
      radf     ! Buoyancy production at CL top due to LW radiative cooling [m^2/s^3]
               ! This is currently set to zero for CLUBB standalone

    real( kind = core_rknd ), dimension(:), allocatable :: &
      Skw_zm, & ! Skewness of w on momentum levels                      [-]
      wp2_zt, & ! w'^2 on thermo. grid                                  [m^2/s^2]
      wpNcp     ! Covariance of w and N_c, <w'N_c'> (momentum levels)   [(m/s)(#/kg)]

    real( kind = core_rknd ), allocatable, dimension(:,:) :: &
      K_hm    ! Eddy diffusivity coef. for hydrometeors on mom. levs. [m^2 s^-1]

    real( kind = core_rknd ) :: &
      T_sfc,     & ! surface temperature     [K]
      p_sfc,     & ! surface pressure        [Pa]
      sens_ht,   & ! sensible heat flux      [K m/s]
      latent_ht    ! latent heat flux        [m/s]

    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      sclrm,     & ! Passive scalar mean (thermo. levels) [units vary]
      wpsclrp,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
      sclrp2,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
      sclrp3,    & ! sclr'^3 (thermodynamic levels)       [{units vary}^3]
      sclrprtp,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp    ! sclr'thl' (momentum levels)          [{units vary} K]

    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      sclrpthvp    ! < sclr' th_v' > (momentum levels)   [units vary]

    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      sclrm_forcing    ! Passive scalar forcing          [{units vary}/s]

    real( kind = core_rknd ), dimension(:), allocatable ::  &
      wpsclrp_sfc      ! Passive scalar flux at surface         [{units vary} m/s]

    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      edsclrm   ! Eddy passive scalar grid-mean (thermo. levels)   [units vary]

    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      edsclrm_forcing  ! Eddy-diffusion passive scalar forcing    [{units vary}/s]

    real( kind = core_rknd ), dimension(:), allocatable ::  &
      wpedsclrp_sfc    ! Eddy-diffusion passive scalar flux at surface [{un vary}m/s]

    real( kind = core_rknd ), dimension(:), allocatable :: &
      radht,        & ! SW + LW heating rate               [K/s]
      Frad,         & ! Radiative flux (momentum levels)   [W/m^2]
      Frad_SW_up,   & ! SW radiative upwelling flux        [W/m^2]
      Frad_LW_up,   & ! LW radiative upwelling flux        [W/m^2]
      Frad_SW_down, & ! SW radiative downwelling flux      [W/m^2]
      Frad_LW_down    ! LW radiative downwelling flux      [W/m^2]

    logical :: l_restart_input

    integer :: k ! Loop iterator(s)

    real( kind = core_rknd ), intent(in), dimension(nparams) ::  & 
      params  ! Model parameters, C1, nu2, etc.

    real( kind = core_rknd ) :: &
      lmin    ! Min. value for the length scale    [m]

    real( kind = core_rknd ), dimension(:), allocatable :: &
      rrm, & ! Overall mean rain water mixing ratio                  [kg/kg]
      Nrm       ! Overall mean rain drop concentration               [num/kg]

    real( kind = core_rknd ), dimension(:,:,:), allocatable :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(:,:,:,:), allocatable :: &
      corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
      corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

    real( kind = core_rknd ), dimension(:,:,:,:), allocatable :: &
      corr_cholesky_mtx_1, & ! Transposed corr. cholesky matrix, 1st comp. [-]
      corr_cholesky_mtx_2    ! Transposed corr. cholesky matrix, 2nd comp. [-]

    real( kind = core_rknd ), dimension(:), allocatable :: &
      radht_zm, &
      rcm_zm

    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      wp2hmp,     & ! Third moment:  <w'^2> * <hm'>       [(m/s)^2 <hm units>]
      rtphmp_zt,  & ! Covariance of rt and a hydrometeor  [(kg/kg) <hm units>]
      thlphmp_zt    ! Covariance of thl and a hydrometeor [K <hm units>]

    real( kind = core_rknd ), dimension(:,:,:,:), allocatable :: &
      X_nl_all_levs    ! Lognormally distributed hydrometeors

    integer, dimension(:,:,:), allocatable :: &
      X_mixt_comp_all_levs ! Which mixture component a sample is in
      
    real( kind = core_rknd ), dimension(:,:,:), allocatable :: &
      lh_rt_clipped,  & ! rt generated from silhs sample points
      lh_thl_clipped, & ! thl generated from silhs sample points
      lh_rc_clipped,  & ! rc generated from silhs sample points
      lh_rv_clipped,  & ! rv generated from silhs sample points
      lh_Nc_clipped     ! Nc generated from silhs sample points

    real( kind = core_rknd ), dimension(:), allocatable :: &
      Nc_in_cloud        ! Mean (in-cloud) cloud droplet concentration  [num/kg]

    real( kind = core_rknd ), dimension(:,:,:), allocatable :: &
      lh_sample_point_weights ! Weights for cloud weighted sampling

    logical :: l_silhs_out    ! Whether to output SILHS files

    type(hydromet_pdf_parameter), dimension(:,:), allocatable :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters      [units vary]
      
    type(precipitation_fractions), allocatable :: &
      precip_fracs           ! Precipitation fractions      [-]
      
    ! coarse-grained timing budget of main time stepping loop
    real( kind = core_rknd ) :: & 
      time_loop_init,  &	   ! time spent in the beginning part of the main loop [s]
      time_loop_end, &             ! time spent in the end part of the main loop [s]
      time_clubb_advance, &        ! time spent in advance_clubb_core [s]
      time_clubb_pdf, &      	   ! time spent in setup_pdf_parameters 
				   !	and hydrometeor_mixed_moments [s]
      time_SILHS,    &             ! time needed to compute subcolumns [s]
      time_microphys_scheme, &     ! time needed for calc_microphys_scheme_tendcies [s]
      time_microphys_advance, &    ! time needed for advance_microphys [s]
      time_stop, time_start,  &    ! help variables to measure the time [s]
      time_total                   ! control timer for the overall time spent in the main loop [s]
      
    ! allowed tolerance for the timing budget check
    real( kind = core_rknd ) , parameter ::  timing_tol = 0.01_core_rknd

    logical, parameter :: &
      l_calc_weights_all_levs = .false. ! .false. if all time steps use the same weights at all grid
                                        ! levels
    
    logical :: &
      l_calc_weights_all_levs_itime, & ! .true. if we calculate sample weights separately at all 
                                       ! grid levels at the current time step
      l_rad_itime                      ! .true. if we calculate radiation at the current time step

    type(silhs_config_flags_type) :: &
      silhs_config_flags ! Flags for the SILHS sampling code

    real( kind = core_rknd ) :: &
      vert_decorr_coef    ! Empirically defined de-correlation constant [-]

    type(nu_vertical_res_dep) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values
  
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
      l_prescribed_avg_deltaz,      & ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz
      l_lmm_stepping,               & ! Apply Linear Multistep Method (LMM) Stepping
      l_e3sm_config,                & ! Run model with E3SM settings
      l_vary_convect_depth,         & ! Flag used to calculate convective velocity using
                                      ! a variable estimate of layer depth based on the depth
                                      ! over which wpthlp is positive near the ground when true
                                      ! More information can be found by
                                      ! Looking at issue #905 on the clubb repo
      l_use_tke_in_wp3_pr_turb_term,& ! Use TKE formulation for wp3 pr_turb term
      l_use_tke_in_wp2_wp3_K_dfsn     ! Use TKE in eddy diffusion for wp2 and wp3

    type(clubb_config_flags_type) :: &
      clubb_config_flags ! Derived type holding all configurable CLUBB flags
    
    ! Definition of namelists
    namelist /model_setting/  &
      runtype, nzmax, grid_type, deltaz, zm_init, zm_top, &
      zt_grid_fname, zm_grid_fname,  &
      day, month, year, lat_vals, lon_vals, sfc_elevation, &
      time_initial, time_final, &
      dt_main, dt_rad, &
      sfctype, T_sfc, p_sfc, sens_ht, latent_ht, fcor, T0, ts_nudge, &
      forcings_file_path, l_t_dependent, l_input_xpwp_sfc, &
      l_ignore_forcings, saturation_formula, &
      thlm_sponge_damp_settings, rtm_sponge_damp_settings, &
      uv_sponge_damp_settings, wp2_sponge_damp_settings, &
      wp3_sponge_damp_settings, up2_vp2_sponge_damp_settings, &
      l_soil_veg, l_uv_nudge, l_restart, restart_path_case, &
      time_restart, l_input_fields, debug_level, &
      sclr_tol, sclr_dim, iisclr_thl, iisclr_rt, iisclr_CO2, &
      edsclr_dim, iiedsclr_thl, iiedsclr_rt, iiedsclr_CO2, &
      l_rtm_nudge, rtm_min, rtm_nudge_max_altitude, &
      l_diagnose_correlations, l_calc_w_corr


    namelist /stats_setting/ &
      l_stats, fname_prefix, stats_tsamp, stats_tout, stats_fmt, &
      l_allow_small_stats_tout

    namelist /configurable_clubb_flags_nl/ &
      iiPDF_type, ipdf_call_placement, &
      l_upwind_xpyp_ta, l_upwind_xm_ma, l_quintic_poly_interp, &
      l_tke_aniso, l_vert_avg_closure, l_standard_term_ta, &
      l_partial_upwind_wp3, l_godunov_upwind_wpxp_ta, l_godunov_upwind_xpyp_ta, &
      l_use_cloud_cover, l_rcm_supersat_adj, &
      l_damp_wp3_Skw_squared, l_min_wp2_from_corr_wx, l_min_xp2_from_corr_wx, &
      l_C2_cloud_frac, l_predict_upwp_vpwp, l_diag_Lscale_from_tau, &
      l_stability_correct_tau_zm, l_damp_wp2_using_em, l_use_C7_Richardson, &
      l_use_precip_frac, l_do_expldiff_rtm_thlm, l_use_C11_Richardson, &
      l_use_shear_Richardson, l_prescribed_avg_deltaz, l_diffuse_rtm_and_thlm, &
      l_stability_correct_Kh_N2_zm, l_trapezoidal_rule_zt, l_trapezoidal_rule_zm, &
      l_call_pdf_closure_twice, l_Lscale_plume_centered, &
      l_brunt_vaisala_freq_moist, l_use_thvm_in_bv_freq, &
      l_lmm_stepping, l_e3sm_config, l_vary_convect_depth, l_use_tke_in_wp3_pr_turb_term, &
      l_use_tke_in_wp2_wp3_K_dfsn
      
    integer :: &
      err_code_dummy ! Host models use an error code that comes out of some API routines, but
                     ! here we have access to the global version

!-----------------------------------------------------------------------
    ! Begin code

    ! Initialize the model run

    ! Pick some default values for model_setting.  Some variables are initialized
    ! at declaration in module clubb_model_settings.

    T_sfc     = 288._core_rknd
    p_sfc     = 1000.e2_core_rknd
    sens_ht   = 0._core_rknd
    latent_ht = 0._core_rknd
    fcor      = 1.e-4_core_rknd
    T0        = 300._core_rknd
    ts_nudge  = 86400._core_rknd

    forcings_file_path = ''
    l_t_dependent   = .false. 
    l_input_xpwp_sfc = .false. 
    l_ignore_forcings = .false. 

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

    rtm_min = 0.0_core_rknd
    rtm_nudge_max_altitude = 0.0_core_rknd

    ! Use the Flatau polynomial approximation for computing saturation in advance_clubb_core_module
    saturation_formula = "flatau"

    sclr_dim   = 0
    iisclr_thl = -1
    iisclr_rt  = -1
    iisclr_CO2 = -1

    edsclr_dim = 0
    iiedsclr_thl = -1
    iiedsclr_rt  = -1
    iiedsclr_CO2 = -1

    sclr_tol(1:sclr_max) = 1.e-2_core_rknd

    l_use_Ncn_to_Nc = .true.

    ! Pick some default values for stats_setting; other variables are set in
    ! module stats_variables
    fname_prefix = ''
    stats_fmt    = ''

    ! Default values for the soil scheme
    call initialize_soil_veg()

    ! Default values for generic model settings
    call initialize_clubb_model_settings()

    ! Figure out which I/O unit to use for OpenMP runs
#ifdef _OPENMP
    iunit = omp_get_thread_num( ) + 10 ! Known magic number
#else
    iunit = 10
#endif

    call set_default_clubb_config_flags( iiPDF_type, & ! Intent(out)
                                         ipdf_call_placement, & ! Intent(out)
                                         l_use_precip_frac, & ! Intent(out)
                                         l_predict_upwp_vpwp, & ! Intent(out)
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
                                         l_use_tke_in_wp2_wp3_K_dfsn ) ! Intent(out)

    ! Read namelist file
    open(unit=iunit, file=trim( runfile ), status='old')
    read(unit=iunit, nml=model_setting)
    read(unit=iunit, nml=stats_setting)
    close(unit=iunit)

    open(unit=iunit, file=runfile, status='old', action='read')
    read(unit=iunit, nml=configurable_clubb_flags_nl)
    close(unit=iunit)

    if ( l_vert_avg_closure ) then
      l_trapezoidal_rule_zt    = .true.
      l_trapezoidal_rule_zm    = .true.
      l_call_pdf_closure_twice = .true.
    else
      l_trapezoidal_rule_zt    = .false.
      l_trapezoidal_rule_zm    = .false.
      l_call_pdf_closure_twice = .false.
    end if

    case_info_file = &
      "../output/" // trim( fname_prefix ) // "_setup.txt" ! The filename for case setup

    ! Sanity check on passive scalars
    ! When adding new 'ii' scalar indices, add them to this list.
    if ( max( iisclr_CO2, iisclr_rt, iisclr_thl ) > sclr_dim ) then
      write(fstderr,*) "Passive scalar index exceeds sclr_dim ", & 
        "iisclr_CO2 = ", iisclr_CO2, "iisclr_rt = ", iisclr_rt,  & 
        "iisclr_thl = ", iisclr_thl, "sclr_dim = ", sclr_dim

      err_code = clubb_fatal_error
      return

    else if ( max( iiedsclr_CO2, iiedsclr_rt, iiedsclr_thl ) > edsclr_dim ) then
      write(fstderr,*) "Passive scalar index exceeds edsclr_dim ", & 
        "iiedsclr_CO2 = ", iiedsclr_CO2, "iiedsclr_rt = ", iiedsclr_rt,  & 
        "iiedsclr_thl = ", iiedsclr_thl, "edsclr_dim = ", edsclr_dim

      err_code = clubb_fatal_error
      return

    end if

    ! Set debug level
    call set_clubb_debug_level( debug_level ) ! Intent(in)

    ! Printing Model Inputs
    if ( clubb_at_least_debug_level( 1 ) ) then

      if ( l_write_to_file ) then
        open(unit=iunit, file=case_info_file, status='replace', action='write')
      end if

      ! Print the date and time
      call write_date( l_write_to_file, iunit )

      ! Print the list of parameters that are being used before the run.
      call write_text( "Parameter          Value", l_write_to_file, iunit, '(4x,A24)')
      call write_text( "---------          -----", l_write_to_file, iunit, '(4x,A24)')
      do j = 1, nparams, 1
        call write_text(params_list(j) // " = ", params(j), & 
          l_write_to_file, iunit, '(A31,F27.20)')
      end do

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

      call write_text( "sfc_elevation = ", sfc_elevation, l_write_to_file, iunit )

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
      call write_text( "T0 = ", T0, l_write_to_file, iunit )
      call write_text( "ts_nudge = ", ts_nudge, l_write_to_file, iunit )

      call write_text( "forcings_file_path = " // forcings_file_path, l_write_to_file, iunit )

      call write_text( "l_t_dependent = ", l_t_dependent, l_write_to_file, iunit )
      call write_text( "l_ignore_forcings = ", l_ignore_forcings, l_write_to_file, iunit )
      call write_text( "l_input_xpwp_sfc = ", l_input_xpwp_sfc, l_write_to_file, iunit )

      call write_text( "saturation_formula = " // saturation_formula, &
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
      call write_text( "l_hole_fill = ", l_hole_fill, l_write_to_file, iunit )
      call write_text( "l_gamma_Skw = ", l_gamma_Skw, l_write_to_file, iunit)
      call write_text( "l_byteswap_io = ", l_byteswap_io, l_write_to_file, iunit )

      call write_text( "Constant tolerances [units]", l_write_to_file, iunit )
      call write_text( "rt_tol [kg/kg] = ", rt_tol, l_write_to_file, iunit )
      call write_text( "thl_tol [K] = ", thl_tol, l_write_to_file, iunit )
      call write_text( "w_tol [m/s] = ", w_tol, l_write_to_file, iunit )

      if ( l_write_to_file ) close(unit=iunit)

    end if ! clubb_at_least_debug_level( 1 )

    !----------------------------------------------------------------------

    ! Allocate stretched grid altitude arrays.
    allocate( momentum_heights(nzmax),  & 
              thermodynamic_heights(nzmax) )

    ! Handle the reading of grid altitudes for
    ! stretched (unevenly-spaced) grid options.
    ! Do some simple error checking for all grid options.
    call read_grid_heights( nzmax, grid_type, &                 ! Intent(in)
                            zm_grid_fname, zt_grid_fname, &     ! Intent(in)
                            iunit, &                            ! Intent(in) 
                            momentum_heights, &                 ! Intent(out)
                            thermodynamic_heights )             ! Intent(out)

    if ( err_code == clubb_fatal_error ) error stop 

    ! These numbers represent host model horizontal grid spacing
    ! which for a single column simulation is effectively infinite
    dummy_dx = 1.0e6_core_rknd ! known magic number
    dummy_dy = 1.0e6_core_rknd ! known magic number

    ! Setup microphysical fields
    call init_microphys( iunit, trim( runtype ), runfile, case_info_file, & ! Intent(in)
                         dummy_dx, dummy_dy, &                              ! Intent(in)
                         params, &                                          ! Intent(in)
                         l_diagnose_correlations, &                         ! Intent(in)
                         l_const_Nc_in_cloud, &                             ! Intent(inout)
                         l_fix_w_chi_eta_correlations, &                    ! Intent(inout)
                         hydromet_dim, silhs_config_flags, &                ! Intent(out)
                         vert_decorr_coef )                                 ! Intent(out)

    ! Setup radiation parameters
    call init_radiation( iunit, runfile, case_info_file, & ! Intent(in)
                         l_calc_thlp2_rad )                ! Intent(inout)

    if ( trim( rad_scheme ) == "lba" ) then
      call simple_rad_lba_init( iunit, trim( forcings_file_path ) )
    end if

    ! Initialize CLUBB configurable flags type
    call initialize_clubb_config_flags_type( iiPDF_type, & ! Intent(in)
                                             ipdf_call_placement, & ! Intent(in)
                                             l_use_precip_frac, & ! Intent(in)
                                             l_predict_upwp_vpwp, & ! Intent(in)
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
                                             clubb_config_flags ) ! Intent(out)

    ! Printing configurable CLUBB flags Inputs
    if ( clubb_at_least_debug_level( 1 ) ) then

      if ( l_write_to_file ) then
        open(unit=iunit, file=case_info_file, status='old', action='write', position='append')
      end if

      call write_text( "--------------------------------------------------", &
                       l_write_to_file, iunit )
      call write_text( "&configurable_clubb_flags_nl", l_write_to_file, iunit )
      call write_text( "--------------------------------------------------", &
                       l_write_to_file, iunit )

       call print_clubb_config_flags( iunit, clubb_config_flags ) ! Intent(in)

      if ( l_write_to_file ) close(unit=iunit)

    end if ! clubb_at_least_debug_level( 1 )

    ! Allocate & initialize variables,
    ! setup grid, setup constants, and setup flags

    call setup_clubb_core                                     & ! Intent(in)
         ( nzmax, T0, ts_nudge,                               & ! Intent(in)
           hydromet_dim, sclr_dim,                            & ! Intent(in)
           sclr_tol(1:sclr_dim), edsclr_dim, params,          & ! Intent(in)
           l_host_applies_sfc_fluxes,                         & ! Intent(in)
           saturation_formula,                                & ! Intent(in)
           l_input_fields,                                    & ! Intent(in)
           l_implemented, grid_type, deltaz, zm_init, zm_top, & ! Intent(in)
           momentum_heights, thermodynamic_heights,           & ! Intent(in)
           sfc_elevation,                                     & ! Intent(in)
           iiPDF_type,                                        & ! intent(in)
           ipdf_call_placement,                               & ! intent(in)
           l_predict_upwp_vpwp,                               & ! intent(in)
           l_min_xp2_from_corr_wx,                            & ! intent(in)
           l_prescribed_avg_deltaz,                           & ! intent(in)
           l_damp_wp2_using_em,                               & ! intent(in)
           l_stability_correct_tau_zm,                        & ! intent(in)
           gr, lmin, nu_vert_res_dep, err_code_dummy )          ! Intent(out)

    ! Allocate and initialize variables

    allocate( um(1:gr%nz) )        ! u wind
    allocate( vm(1:gr%nz) )        ! v wind

    allocate( upwp(1:gr%nz) )      ! vertical u momentum flux
    allocate( vpwp(1:gr%nz) )      ! vertical v momentum flux

    allocate( up2(1:gr%nz) )
    allocate( up3(1:gr%nz) )
    allocate( vp2(1:gr%nz) )
    allocate( vp3(1:gr%nz) )

    allocate( thlm(1:gr%nz) )      ! liquid potential temperature
    allocate( rtm(1:gr%nz) )       ! total water mixing ratio
    allocate( wprtp(1:gr%nz) )     ! w'rt'
    allocate( wpthlp(1:gr%nz) )    ! w'thl'
    allocate( wprcp(1:gr%nz) )     ! w'rc'
    allocate( w_up_in_cloud(1:gr%nz) )
    allocate( wp2(1:gr%nz) )       ! w'^2
    allocate( wp3(1:gr%nz) )       ! w'^3
    allocate( rtp2(1:gr%nz) )      ! rt'^2
    allocate( thlp2(1:gr%nz) )     ! thl'^2
    allocate( rtpthlp(1:gr%nz) )   ! rt'thlp'

    allocate( p_in_Pa(1:gr%nz) )         ! pressure (pascals)
    allocate( exner(1:gr%nz) )           ! exner function
    allocate( rho(1:gr%nz) )             ! density: t points
    allocate( rho_zm(1:gr%nz) )          ! density: m points
    allocate( rho_ds_zm(1:gr%nz) )       ! dry, static density: m-levs
    allocate( rho_ds_zt(1,1:gr%nz) )     ! dry, static density: t-levs
    allocate( invrs_rho_ds_zm(1:gr%nz) ) ! inv. dry, static density: m-levs
    allocate( invrs_rho_ds_zt(1:gr%nz) ) ! inv. dry, static density: t-levs
    allocate( thv_ds_zm(1:gr%nz) )       ! dry, base-state theta_v: m-levs
    allocate( thv_ds_zt(1:gr%nz) )       ! dry, base-state theta_v: t-levs

    allocate( thlm_forcing(1:gr%nz) )    ! thlm ls forcing
    allocate( rtm_forcing(1:gr%nz) )     ! rtm ls forcing
    allocate( um_forcing(1:gr%nz) )      ! u forcing
    allocate( vm_forcing(1:gr%nz) )      ! v forcing
    allocate( wprtp_forcing(1:gr%nz) )   ! <w'r_t'> forcing (microphysics)
    allocate( wpthlp_forcing(1:gr%nz) )  ! <w'th_l'> forcing (microphysics)
    allocate( rtp2_forcing(1:gr%nz) )    ! <r_t'^2> forcing (microphysics)
    allocate( thlp2_forcing(1:gr%nz) )   ! <th_l'^2> forcing (microphysics)
    allocate( rtpthlp_forcing(1:gr%nz) ) ! <r_t'th_l'> forcing (microphysics)

    ! Imposed large scale w
    allocate( wm_zm(1:gr%nz) )       ! momentum levels
    allocate( wm_zt(1:gr%nz) )       ! thermodynamic levels

    ! Cloud water variables
    allocate( rcm(1,1:gr%nz) )
    allocate( delta_zm(1,1:gr%nz) )
    allocate( cloud_frac(1:gr%nz) )
    allocate( ice_supersat_frac(1:gr%nz) )
    allocate( rcm_in_layer(1:gr%nz) )
    allocate( cloud_cover(1:gr%nz) )
    allocate( invrs_tau_zm(1:gr%nz) )

    ! Passive scalar variables
    ! Note that sclr_dim can be 0
    allocate( wpsclrp_sfc(1:sclr_dim) )
    allocate( sclrm(1:gr%nz, 1:sclr_dim) )
    allocate( sclrp2(1:gr%nz, 1:sclr_dim) )
    allocate( sclrp3(1:gr%nz, 1:sclr_dim) )
    allocate( sclrm_forcing(1:gr%nz, 1:sclr_dim) )
    allocate( sclrprtp(1:gr%nz, 1:sclr_dim) )
    allocate( sclrpthlp(1:gr%nz, 1:sclr_dim) )

    allocate( wpedsclrp_sfc(1:edsclr_dim) )
    allocate( edsclrm_forcing(1:gr%nz, 1:edsclr_dim) )

    allocate( edsclrm(1:gr%nz, 1:edsclr_dim) )
    allocate( wpsclrp(1:gr%nz, 1:sclr_dim) )

    allocate( sigma_sqd_w(1:gr%nz) )    ! PDF width parameter (momentum levels)
    allocate( sigma_sqd_w_zt(1:gr%nz) ) ! PDF width parameter interp. to t-levs.
    allocate( Skw_zm(1:gr%nz) )         ! Skewness of w on momentum levels
    allocate( wp2_zt(1:gr%nz) )         ! wp2 interpolated to thermo. levels
    allocate( ug(1:gr%nz) )             ! u geostrophic wind
    allocate( vg(1:gr%nz) )             ! v geostrophic wind
    allocate( um_ref(1:gr%nz) )         ! Reference u wind for nudging; Michael Falk, 17 Oct 2007
    allocate( vm_ref(1:gr%nz) )         ! Reference v wind for nudging; Michael Falk, 17 Oct 2007
    allocate( thlm_ref(1:gr%nz) )       ! Reference liquid water potential for nudging
    allocate( rtm_ref(1:gr%nz) )        ! Reference total water mixing ratio for nudging
    allocate( thvm(1:gr%nz) )           ! Virtual potential temperature
    allocate( radht(1:gr%nz) )          ! SW + LW heating rate
    allocate( Frad(1:gr%nz) )           ! radiative flux (momentum point)
    allocate( Frad_SW_up(1:gr%nz) )
    allocate( Frad_LW_up(1:gr%nz) )
    allocate( Frad_SW_down(1:gr%nz) )
    allocate( Frad_LW_down(1:gr%nz) )
    allocate( thlprcp(1:gr%nz) )   ! thl'rc'
    allocate( thlp3(1:gr%nz) )     ! thl'^3
    allocate( rtp3(1:gr%nz) )      ! rt'^3

    ! Buoyancy related moments
    allocate( rtpthvp(1:gr%nz) )  ! rt'thv'
    allocate( thlpthvp(1:gr%nz) ) ! thl'thv'
    allocate( wpthvp(1:gr%nz) )   ! w'thv'
    allocate( wp2thvp(1:gr%nz) )  ! w'^2thv'

    allocate( uprcp(1:gr%nz) )  ! u'rc'
    allocate( vprcp(1:gr%nz) )  ! v'rc'

    allocate( Kh_zt(1:gr%nz) )  ! Eddy diffusivity coefficient: thermo. levels
    allocate( Kh_zm(1:gr%nz) )  ! Eddy diffusivity coefficient: momentum levels
    allocate( K_hm(1:gr%nz,1:hydromet_dim) ) ! Eddy diff. coef. for hydromets.: mom. levs.

    allocate( em(1:gr%nz) )
    allocate( Lscale(1,1:gr%nz) )
    allocate( Lscale_up(1:gr%nz) )
    allocate( Lscale_down(1:gr%nz) )

    allocate( tau_zm(1:gr%nz) ) ! Eddy dissipation time scale: momentum levels
    allocate( tau_zt(1:gr%nz) ) ! Eddy dissipation time scale: thermo. levels

    allocate( Nccnm(1:gr%nz) )
    allocate( hydromet(1:gr%nz,1:hydromet_dim) )    ! All hydrometeor mean fields
    allocate( hydrometp2(1:gr%nz,1:hydromet_dim) )  ! All < h_m'^2 > fields
    allocate( wphydrometp(1:gr%nz,1:hydromet_dim) ) ! All < w'h_m' > fields
    allocate( Ncm(1:gr%nz) )   ! Mean cloud droplet concentration, < N_c >
    allocate( wpNcp(1:gr%nz) ) ! < w'N_c' >

    ! High-order passive scalars
    allocate( sclrpthvp(1:gr%nz, 1:sclr_dim) )

    ! Variables for PDF closure scheme
    allocate( pdf_params )
    call init_pdf_params( gr%nz, 1, pdf_params )
    allocate( pdf_params_zm )
    call init_pdf_params( gr%nz, 1, pdf_params_zm )

    allocate( pdf_implicit_coefs_terms )
    call init_pdf_implicit_coefs_terms( gr%nz, sclr_dim, &         ! Intent(in)
                                        pdf_implicit_coefs_terms ) ! Intent(out)

    um(1:gr%nz)      = zero          ! u wind
    vm (1:gr%nz)     = zero          ! v wind
    upwp(1:gr%nz)    = zero          ! vertical u momentum flux
    vpwp(1:gr%nz)    = zero          ! vertical v momentum flux
    up2(1:gr%nz)     = w_tol_sqd     ! u'^2
    up3(1:gr%nz)     = zero          ! u'^3
    vp2(1:gr%nz)     = w_tol_sqd     ! v'^2
    vp3(1:gr%nz)     = zero          ! v'^3

    thlm(1:gr%nz)    = zero          ! liquid potential temperature
    rtm(1:gr%nz)     = zero          ! total water mixing ratio
    wprtp(1:gr%nz)   = zero          ! w'rt'
    wpthlp(1:gr%nz)  = zero          ! w'thl'
    wp2(1:gr%nz)     = w_tol_sqd     ! w'^2
    wp3(1:gr%nz)     = zero          ! w'^3
    rtp2(1:gr%nz)    = rt_tol**2     ! rt'^2
    thlp2(1:gr%nz)   = thl_tol**2    ! thl'^2
    rtpthlp(1:gr%nz) = zero          ! rt'thl'
    wprcp(1:gr%nz)   = zero          ! w'rc'
    w_up_in_cloud(1:gr%nz) = zero

    p_in_Pa(1:gr%nz)= zero           ! pressure (Pa)
    exner(1:gr%nz) = zero            ! exner
    rho(1:gr%nz)  = zero             ! density on thermo. levels
    rho_zm(1:gr%nz)  = zero          ! density on moment. levels
    rho_ds_zm(1:gr%nz) = zero        ! dry, static density: m-levs
    rho_ds_zt(1,1:gr%nz) = zero      ! dry, static density: t-levs
    invrs_rho_ds_zm(1:gr%nz) = zero  ! inv. dry, static density: m-levs
    invrs_rho_ds_zt(1:gr%nz) = zero  ! inv. dry, static density: t-levs
    thv_ds_zm(1:gr%nz) = zero        ! dry, base-state theta_v: m-levs
    thv_ds_zt(1:gr%nz) = zero        ! dry, base-state theta_v: t-levs

    thlm_forcing(1:gr%nz)    = zero  ! thlm large-scale forcing
    rtm_forcing(1:gr%nz)     = zero  ! rtm large-scale forcing
    um_forcing(1:gr%nz)      = zero  ! u forcing
    vm_forcing(1:gr%nz)      = zero  ! v forcing
    wprtp_forcing(1:gr%nz)   = zero  ! <w'r_t'> forcing (microphysics)
    wpthlp_forcing(1:gr%nz)  = zero  ! <w'th_l'> forcing (microphysics)
    rtp2_forcing(1:gr%nz)    = zero  ! <r_t'^2> forcing (microphysics)
    thlp2_forcing(1:gr%nz)   = zero  ! <th_l'^2> forcing (microphysics)
    rtpthlp_forcing(1:gr%nz) = zero  ! <r_t'th_l'> forcing (microphysics)

    ! Imposed large scale w
    wm_zm(1:gr%nz) = zero      ! Momentum levels
    wm_zt(1:gr%nz) = zero      ! Thermodynamic levels

    ! Cloud water variables
    rcm(1,1:gr%nz)             = zero
    cloud_frac(1:gr%nz)        = zero
    ice_supersat_frac(1:gr%nz) = zero
    rcm_in_layer(1:gr%nz)      = zero
    cloud_cover(1:gr%nz)       = zero
    invrs_tau_zm(1:gr%nz)      = zero

    sigma_sqd_w    = zero ! PDF width parameter (momentum levels)
    sigma_sqd_w_zt = zero ! PDF width parameter interp. to t-levs.
    Skw_zm         = zero ! Skewness of w on momentum levels
    wp2_zt         = w_tol_sqd ! wp2 interpolated to thermo. levels
    ug             = zero ! u geostrophic wind
    vg             = zero ! v geostrophic wind
    um_ref         = zero
    vm_ref         = zero
    thlm_ref       = zero
    rtm_ref        = zero
    thvm           = zero ! Virtual potential temperature

    radht = zero ! Heating rate
    Frad  = zero ! Radiative flux
    Frad_SW_up = zero
    Frad_LW_up = zero
    Frad_SW_down = zero
    Frad_LW_down = zero
    thlprcp = zero
    thlp3   = zero
    rtp3    = zero

    ! Buoyancy related moments
    rtpthvp  = zero ! rt'thv'
    thlpthvp = zero ! thl'thv'
    wpthvp   = zero ! w'thv'
    wp2thvp  = zero ! w'^2thv'

    uprcp  = zero ! u'rc'
    vprcp  = zero ! v'rc'

    ! Eddy diffusivity
    Kh_zt = zero  ! Eddy diffusivity coefficient: thermo. levels
    Kh_zm = zero  ! Eddy diffusivity coefficient: momentum levels

    do i = 1, hydromet_dim, 1
      K_hm(1:gr%nz,i) = zero ! Eddy diff. coef. for hydromets.: mom. levs.
    end do

    ! TKE
    em = em_min

    ! Length scale
    Lscale      = zero
    Lscale_up   = zero
    Lscale_down = zero

    ! Dissipation time
    tau_zm = zero ! Eddy dissipation time scale: momentum levels
    tau_zt = zero ! Eddy dissipation time scale: thermo. levels

    ! Hydrometer types
    Nccnm(1:gr%nz) = zero ! CCN concentration (COAMPS/MG)

    do i = 1, hydromet_dim, 1
       hydromet(1:gr%nz,i)    = zero
       hydrometp2(1:gr%nz,i)  = zero
       wphydrometp(1:gr%nz,i) = zero
    enddo

    ! Cloud droplet concentration
    Ncm(1:gr%nz)   = zero
    wpNcp(1:gr%nz) = zero

    ! Passive scalars
    if ( sclr_dim > 0 ) then
      sclrpthvp(:,:) = zero
    end if

    ! Surface fluxes
    wpthlp_sfc = zero
    wprtp_sfc  = zero
    upwp_sfc   = zero
    vpwp_sfc   = zero

    ! Passive scalars
    do i = 1, sclr_dim, 1
       wpsclrp_sfc(i)           = zero
       sclrm(1:gr%nz,i)         = zero
       sclrp2(1:gr%nz,i)        = sclr_tol(i)**2
       sclrp3(1:gr%nz,i)        = zero
       sclrprtp(1:gr%nz,i)      = zero
       sclrpthlp(1:gr%nz,i)     = zero
       sclrm_forcing(1:gr%nz,i) = zero
       wpsclrp(1:gr%nz,i)       = zero
    enddo

    do i = 1, edsclr_dim, 1
      wpedsclrp_sfc(i)           = zero
      edsclrm(1:gr%nz,i)         = zero
      edsclrm_forcing(1:gr%nz,i) = zero
    end do

    ! Allocate a correctly-sized array for radf and zero it
    allocate( radf(gr%nz) )

    ! Zero all elements of radf
    radf(1:gr%nz) = 0.0_core_rknd

    allocate( wp2hmp(gr%nz,hydromet_dim), rtphmp_zt(gr%nz,hydromet_dim), &
              thlphmp_zt(gr%nz,hydromet_dim) )

    wp2hmp(:,:)     = 0._core_rknd
    rtphmp_zt(:,:)  = 0._core_rknd
    thlphmp_zt(:,:) = 0._core_rknd

    if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then
          ! At this point, input fields haven't been set up, so don't clean them up.
          call cleanup_clubb( l_input_fields=.false. )
          return
        end if
    end if

    ! This special purpose code only applies to tuner runs where the tune_type
    ! is setup to try all permutations of our model flags
    if ( present( model_flags_array ) ) then
      clubb_config_flags%l_godunov_upwind_wpxp_ta = model_flags_array(1)
      clubb_config_flags%l_godunov_upwind_xpyp_ta = model_flags_array(2)
      clubb_config_flags%l_upwind_xpyp_ta = model_flags_array(3)
      clubb_config_flags%l_upwind_xm_ma = model_flags_array(4)
      l_quintic_poly_interp = model_flags_array(5)
      clubb_config_flags%l_vert_avg_closure = model_flags_array(6)
      clubb_config_flags%l_standard_term_ta = model_flags_array(7)
      clubb_config_flags%l_tke_aniso = model_flags_array(8)
      clubb_config_flags%l_use_cloud_cover = model_flags_array(9)
      clubb_config_flags%l_rcm_supersat_adj = model_flags_array(10)

      if ( clubb_config_flags%l_vert_avg_closure ) then
        clubb_config_flags%l_trapezoidal_rule_zt    = .true.
        clubb_config_flags%l_trapezoidal_rule_zm    = .true.
        clubb_config_flags%l_call_pdf_closure_twice = .true.
      else
        clubb_config_flags%l_trapezoidal_rule_zt    = .false.
        clubb_config_flags%l_trapezoidal_rule_zm    = .false.
        clubb_config_flags%l_call_pdf_closure_twice = .false.
      end if
    end if

    ! Deallocate stretched grid altitude arrays
    deallocate( momentum_heights, thermodynamic_heights )

    ! Allocate rvm_mc, rcm_mc, thlm_mc
    allocate( rvm_mc(gr%nz), rcm_mc(gr%nz), thlm_mc(gr%nz) )

    ! Allocate hydrometeor variables.
    allocate( rrm(gr%nz) )
    allocate( Nrm(gr%nz) )

    ! Allocate hydromet_pdf_params
    allocate( hydromet_pdf_params(1,gr%nz) )
    
    allocate( precip_fracs )

    ! Initialize to 0.
    do k=1,gr%nz
      call init_hydromet_pdf_params( hydromet_pdf_params(1,k) )
    end do
    
    call init_precip_fracs_api( gr%nz, 1, &
                                precip_fracs )

    ! Allocate the correlation arrays
    allocate(corr_array_1_n(1,gr%nz,pdf_dim,pdf_dim))
    allocate(corr_array_2_n(1,gr%nz,pdf_dim,pdf_dim))
    allocate(corr_cholesky_mtx_1(1,gr%nz,pdf_dim,pdf_dim))
    allocate(corr_cholesky_mtx_2(1,gr%nz,pdf_dim,pdf_dim))

    ! Allocate the mean and stddev arrays
    allocate(mu_x_1_n(1,gr%nz,pdf_dim))
    allocate(mu_x_2_n(1,gr%nz,pdf_dim))
    allocate(sigma_x_1_n(1,gr%nz,pdf_dim))
    allocate(sigma_x_2_n(1,gr%nz,pdf_dim))

    ! Initialize to 0.
    rvm_mc  = zero
    rcm_mc  = zero
    thlm_mc = zero

    ! Allocate microphysics tendencies for <w'rt'>, <w'thl'>,
    ! <rt'^2>, <thl'^2>, and <rt'thl'>.
    allocate( wprtp_mc(gr%nz) )
    allocate( wpthlp_mc(gr%nz) )
    allocate( rtp2_mc(gr%nz) )
    allocate( thlp2_mc(gr%nz) )
    allocate( rtpthlp_mc(gr%nz) )

    ! Initialize to 0.
    wprtp_mc   = zero
    wpthlp_mc  = zero
    rtp2_mc    = zero
    thlp2_mc   = zero
    rtpthlp_mc = zero

    ! Allocate microphysics tendencies for mean hydrometeors and mean cloud
    ! droplet concentration.
    allocate( hydromet_mc(gr%nz,hydromet_dim) )
    allocate( Ncm_mc(gr%nz) )

    ! Allocate mean sedimentation velocity for hydrometeors
    allocate( hydromet_vel_zt(gr%nz,hydromet_dim) )

    ! Allocate covariance sedimentation velocities for <V_hm'hm'>, implicit and
    ! explicit components, such that:
    ! <V_hm'hm'> = <V_hm'hm'>|_impc * <hm> + <V_hm'hm'>|_expc.
    allocate( hydromet_vel_covar_zt_impc(gr%nz,hydromet_dim) )
    allocate( hydromet_vel_covar_zt_expc(gr%nz,hydromet_dim) )

    ! Initialize to 0.
    hydromet_mc = zero
    Ncm_mc      = zero

    hydromet_vel_zt = zero

    hydromet_vel_covar_zt_impc = zero
    hydromet_vel_covar_zt_expc = zero

    allocate( rfrzm(gr%nz) )
    rfrzm = zero

    allocate( rcm_zm(gr%nz), radht_zm(gr%nz) )
              
    allocate( X_nl_all_levs(1,lh_num_samples,gr%nz,pdf_dim), &
              X_mixt_comp_all_levs(1,lh_num_samples,gr%nz), &
              lh_rt_clipped(1,lh_num_samples,gr%nz), &
              lh_thl_clipped(1,lh_num_samples,gr%nz), &
              lh_rc_clipped(1,lh_num_samples,gr%nz), &
              lh_rv_clipped(1,lh_num_samples,gr%nz), &
              lh_Nc_clipped(1,lh_num_samples,gr%nz), &
              lh_sample_point_weights(1,lh_num_samples,gr%nz), &
              Nc_in_cloud(gr%nz) )

    if ( .not. l_restart ) then

      time_current = time_initial
      iinit = 1

      call initialize_clubb &
           ( gr, iunit, trim( forcings_file_path ), p_sfc, zm_init,  & ! Intent(in)
             clubb_config_flags%l_uv_nudge,                      & ! Intent(in)
             clubb_config_flags%l_tke_aniso,                     & ! Intent(in)
             thlm, rtm, um, vm, ug, vg, wp2, up2, vp2, rcm(1,:), & ! Intent(inout)
             wm_zt, wm_zm, em, exner,                            & ! Intent(inout)
             thvm, p_in_Pa,                                      & ! Intent(inout)
             rho, rho_zm, rho_ds_zm, rho_ds_zt(1,:),             & ! Intent(inout)
             invrs_rho_ds_zm, invrs_rho_ds_zt,                   & ! Intent(inout)
             thv_ds_zm, thv_ds_zt,                               & ! Intent(inout)
             rtm_ref, thlm_ref,                                  & ! Intent(inout) 
             um_ref, vm_ref,                                     & ! Intent(inout)
             Ncm, Nc_in_cloud, Nccnm,                            & ! Intent(inout)
             sclrm, edsclrm )                                      ! Intent(out)

      if ( clubb_at_least_debug_level( 0 ) ) then
          if ( err_code == clubb_fatal_error ) then
              ! At this point, input fields haven't been set up, so don't clean them up.
              call cleanup_clubb( l_input_fields=.false. )
              return
          end if
      end if

    else  ! restart


      ! Currently initialize_clubb does more than just read in the initial sounding.
      ! It also includes other important initializations such as um_ref and vm_ref.
      ! Therefore it should be executed prior to a restart. The restart should overwrite
      ! the initial sounding anyway.
      call initialize_clubb &
           ( gr, iunit, trim( forcings_file_path ), p_sfc, zm_init,  & ! Intent(in)
             clubb_config_flags%l_uv_nudge,                      & ! Intent(in)
             clubb_config_flags%l_tke_aniso,                     & ! Intent(in)
             thlm, rtm, um, vm, ug, vg, wp2, up2, vp2, rcm(1,:), & ! Intent(inout)
             wm_zt, wm_zm, em, exner,                            & ! Intent(inout)
             thvm, p_in_Pa,                                      & ! Intent(inout)
             rho, rho_zm, rho_ds_zm, rho_ds_zt(1,:),             & ! Intent(inout)
             invrs_rho_ds_zm, invrs_rho_ds_zt,                   & ! Intent(inout)
             thv_ds_zm, thv_ds_zt,                               & ! Intent(inout)
             rtm_ref, thlm_ref,                                  & ! Intent(inout) 
             um_ref, vm_ref,                                     & ! Intent(inout)
             Ncm, Nc_in_cloud, Nccnm,                            & ! Intent(inout)
             sclrm, edsclrm )                                      ! Intent(out)

      if ( clubb_at_least_debug_level( 0 ) ) then
          if ( err_code == clubb_fatal_error ) then
              ! At this point, input fields haven't been set up, so don't clean them up.
              call cleanup_clubb( l_input_fields=.false. )
              return
          end if
      end if

      time_current = time_restart

      ! Determining what iteration to restart at.
      ! The value is increased by 1 to sychronize with restart data.
      ! Joshua Fasching February 2008

      ! Ensure that iteration num, iinit, is an integer, so that model time is
      !   incremented correctly by iteration number at end of timestep
      if ( abs(mod((time_restart-time_initial),real(dt_main, kind=time_precision))) > &
            real(eps, kind=time_precision) ) then

        write(fstderr,*) "Error: (time_restart-time_initial) ",  & 
          "is not a multiple of dt_main."
        write(fstderr,*) "time_restart = ", time_restart
        write(fstderr,*) "time_initial = ", time_initial
        write(fstderr,*) "dt_main = ", dt_main
        error stop "Fatal error"

      end if ! mod( (time_restart-time_initial) , dt_main ) /= 0

      iinit = floor( ( time_current - time_initial ) / real(dt_main,kind=time_precision) ) + 1

      call restart_clubb &
           ( iunit, runfile,                           & ! Intent(in)
             restart_path_case, time_restart,          & ! Intent(in)
             um, upwp, vm, vpwp, up2, vp2, rtm,        & ! Intent(inout)
             wprtp, thlm, wpthlp, rtp2, rtp3,          & ! Intent(inout)
             thlp2, thlp3, rtpthlp, wp2, wp3,          & ! Intent(inout)
             p_in_Pa, exner, rcm(1,:), cloud_frac,     & ! Intent(inout)
             wpthvp, wp2thvp, rtpthvp, thlpthvp,       & ! Intent(inout)
             wm_zt, rho, rho_zm, rho_ds_zm,            & ! Intent(inout)
             rho_ds_zt(1,:), thv_ds_zm, thv_ds_zt,     & ! Intent(inout)
             thlm_forcing, rtm_forcing, wprtp_forcing, & ! Intent(inout)
             wpthlp_forcing, rtp2_forcing,             & ! Intent(inout)
             thlp2_forcing, rtpthlp_forcing,           & ! Intent(inout)
             hydromet, hydrometp2, wphydrometp,        & ! Intent(inout)
             Ncm, Nccnm, thvm, em, tau_zm, tau_zt,     & ! Intent(inout)
             Kh_zt, Kh_zm, ug, vg, Lscale(1,:),        & ! Intent(inout)
             Lscale_up, Lscale_down, thlprcp,          & ! Intent(inout)
             sigma_sqd_w, sigma_sqd_w_zt, radht,       & ! Intent(inout)
             pdf_params, pdf_params_zm,                & ! Intent(inout)
             rcm_mc, rvm_mc, thlm_mc,                  & ! Intent(out)
             wprtp_mc, wpthlp_mc, rtp2_mc,             & ! Intent(out)
             thlp2_mc, rtpthlp_mc,                     & ! Intent(out)
             wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc ) ! Intent(out)
 
      ! Calculate invrs_rho_ds_zm and invrs_rho_ds_zt from the values of
      ! rho_ds_zm and rho_ds_zt, respectively, which were read in from the input
      ! file during the call to subroutine restart_clubb.
      invrs_rho_ds_zm = 1.0_core_rknd/rho_ds_zm
      invrs_rho_ds_zt = 1.0_core_rknd/rho_ds_zt(1,:)

    end if ! ~l_restart

    call setup_radiation_variables( gr%nz, lin_int_buffer, &
                                    extended_atmos_range_size )

#ifdef _OPENMP
    iunit = omp_get_thread_num( ) + 50 ! Known magic number
#else
    iunit = 50
#endif

    fdir = "../output/" ! Output directory

    ! Only output radiation files if using a radiation scheme
    if ( trim( rad_scheme ) == "bugsrad" ) then
      l_output_rad_files = .true.
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

    ! This is a kludge added because the grid used by BUGSrad does
    ! not include CLUBB's ghost point. -nielsenb 20 Oct 2009
    if ( l_output_rad_files ) then
      ! Initialize statistics output
      call stats_init( iunit, fname_prefix, fdir, l_stats, & ! Intent(in)
                       stats_fmt, stats_tsamp, stats_tout, runfile, & ! Intent(in)
                       gr%nz, nlon, nlat, gr%zt, gr%zm, total_atmos_dim - 1, & ! Intent(in)
                       complete_alt(2:total_atmos_dim), total_atmos_dim, & ! Intent(in)
                       complete_momentum(2:total_atmos_dim + 1), day, month, year, & ! Intent(in)
                       (/lon_vals/), (/lat_vals/), time_current, dt_main, l_silhs_out,&!intent(in)
                       stats_zt, stats_zm, stats_sfc, & ! intent(inout)
                       stats_lh_zt, stats_lh_sfc, & ! intent(inout)
                       stats_rad_zt, stats_rad_zm ) ! intent(inout)
    else
      ! Initialize statistics output
      call stats_init( iunit, fname_prefix, fdir, l_stats, & ! Intent(in)
                       stats_fmt, stats_tsamp, stats_tout, runfile, & ! Intent(in)
                       gr%nz, nlon, nlat, gr%zt, gr%zm, 0, & ! Intent(in)
                       rad_dummy, 0, rad_dummy, day, month, year, & ! Intent(in)
                       (/lon_vals/), (/lat_vals/), time_current, dt_main, l_silhs_out,&!intent(in)
                       stats_zt, stats_zm, stats_sfc, & ! intent(inout)
                       stats_lh_zt, stats_lh_sfc, & ! intent(inout)
                       stats_rad_zt, stats_rad_zm ) ! intent(inout)
    end if
 

#ifdef SILHS
    if ( lh_microphys_type /= lh_microphys_disabled ) then

      ! Setup 2D output of all subcolumns (if enabled)
      call latin_hypercube_2D_output_api &
           ( fname_prefix, fdir, stats_tout, gr%nz, & ! Intent(in)
             gr%zt, time_initial, lh_num_samples, & ! Intent(in)
             nlon, nlat, (/lon_vals/), (/lat_vals/) )    ! Intent(in)

    end if
#endif /* SILHS */

    ! Time integration
    ! Call advance_clubb_core once per each statistics output time
    ifinal = floor( ( time_final - time_initial ) / real(dt_main, kind=time_precision) )

    ! Setup filenames and variables to set for setfields, if enabled
    if ( l_input_fields ) then
      call inputfields_init( iunit, runfile ) ! Intent(in)
    end if

    ! check to make sure dt_rad is a mutliple of dt_main
    if ( abs(dt_rad/dt_main - real(floor(dt_rad/dt_main), kind=core_rknd)) &
          > 1.e-8_core_rknd) then
      error stop "dt_rad must be a multiple of dt_main"
    end if

    if( l_stats ) then

      if ( .not. (( abs(dt_rad/stats_tout - real(floor(dt_rad/stats_tout), kind=core_rknd)) &
            < 1.e-8_core_rknd) .or. &
         ( abs(stats_tout/dt_rad - real(floor(stats_tout/dt_rad), kind=core_rknd)) &
            < 1.e-8_core_rknd)) ) then
        error stop &
              "dt_rad must be a multiple of stats_tout or stats_tout must be a mulitple of dt_rad"
      end if

    end if

    if ( l_silhs_rad .and. clubb_config_flags%l_calc_thlp2_rad ) then

      write(fstderr,*) "The options l_silhs_rad and l_calc_thlp2_rad are incompatible."
      err_code = clubb_fatal_error
      call cleanup_clubb( l_input_fields )
      return

    end if

    if ( clubb_config_flags%l_calc_thlp2_rad .and. rad_scheme == "none" ) then

      error stop "The options rad_scheme == none and l_calc_thlp2_rad are incompatible."

    end if
    stats_nsamp = nint( stats_tsamp / dt_main )
    stats_nout = nint( stats_tout / dt_main )

    !initialize timers    
    time_loop_init = 0.0_core_rknd
    time_clubb_advance = 0.0_core_rknd
    time_clubb_pdf = 0.0_core_rknd
    time_SILHS = 0.0_core_rknd
    time_microphys_advance = 0.0_core_rknd
    time_microphys_scheme = 0.0_core_rknd
    time_loop_end = 0.0_core_rknd
    time_total = 0.0_core_rknd
    time_stop = 0.0_core_rknd
    time_start = 0.0_core_rknd
    
    ! Save time before main loop starts
    call cpu_time( time_start )
    time_total = time_start
!-------------------------------------------------------------------------------
!                         Main Time Stepping Loop
!-------------------------------------------------------------------------------

    do itime = iinit, ifinal, 1
      
      call cpu_time( time_start ) ! start timer for initial part of main loop
      
      ! When this time step is over, the time will be time + dt_main
      ! We use integer timestep for stats_begin_step
        call stats_begin_timestep( itime, stats_nsamp, stats_nout ) ! Intent(in)

      ! If we're doing an inputfields run, get the values for our
      ! model arrays from a netCDF or GrADS file.
      ! Note:  the time of the 1st LES statistical output is time_initial_LES
      !        plus the time of one LES statistical time step.  For example, a
      !        LES run that starts at 00:00Z and outputs stats every minute has
      !        its first statistical output at 00:01Z.  In order to match the
      !        LES stat output times to CLUBB timesteps, CLUBB's time_current
      !        + dt_main needs to be passed into subroutine compute_timestep.
      if ( l_input_fields ) then
        l_restart_input = .false.
        call compute_timestep( &
             iunit, stat_files(1), l_restart_input, &            ! Intent(in)
             time_current + real(dt_main,kind=time_precision), & ! Intent(in)
             itime_nearest )                                     ! Intent(out)

        call stat_fields_reader( gr, max( itime_nearest, 1 ), & ! In
                                 um, upwp, vm, vpwp, up2, vp2, rtm, & ! Inout
                                 wprtp, thlm, wpthlp, rtp2, rtp3, & ! Inout
                                 thlp2, thlp3, rtpthlp, wp2, wp3, & ! Inout
                                 p_in_Pa, exner, rcm(1,:), cloud_frac, & ! Inout
                                 wpthvp, wp2thvp, rtpthvp, thlpthvp, & ! Inout
                                 wm_zt, rho, rho_zm, rho_ds_zm, & ! Inout
                                 rho_ds_zt(1,:), thv_ds_zm, thv_ds_zt, & ! Inout
                                 thlm_forcing, rtm_forcing, wprtp_forcing, & !""
                                 wpthlp_forcing, rtp2_forcing, & ! Inout
                                 thlp2_forcing, rtpthlp_forcing, & ! Inout
                                 hydromet, hydrometp2, wphydrometp, & ! Inout
                                 Ncm, Nccnm, thvm, em, tau_zm, tau_zt, & ! Inout
                                 Kh_zt, Kh_zm, ug, vg, Lscale(1,:), & ! Inout
                                 Lscale_up, Lscale_down, thlprcp, & ! Inout
                                 sigma_sqd_w, sigma_sqd_w_zt, radht, & ! Inout
                                 pdf_params, pdf_params_zm ) ! Inout

        ! clip wp3 if it is input from inputfields
        ! this helps restrict the skewness of wp3_on_wp2
        if( l_input_wp3 ) then
          wp2_zt = max( zm2zt( gr, wp2 ), w_tol_sqd ) ! Positive definite quantity
          call clip_skewness_core( gr, sfc_elevation, params(iSkw_max_mag), &
                                   wp2_zt, wp3 )
        end if
      end if

      ! Check for NaN values in the model arrays
      if ( invalid_model_arrays( gr, um, vm, rtm, wprtp, thlm, wpthlp, &
                                 rtp2, thlp2, rtpthlp, wp2, wp3, &
                                 wp2thvp, rtpthvp, thlpthvp, &
                                 hydromet, sclrm, edsclrm ) ) then
        err_code = clubb_fatal_error
        write(fstderr,*) "Fatal error: a CLUBB variable is NaN in main time stepping loop."
        error stop
      end if

      ! Calculate radiation only once in a while
      l_rad_itime = (mod( itime, floor(dt_rad/dt_main) ) == 0 .or. itime == 1)
      
      ! Calculate sample weights separately at all grid levels when radiation is not called
      l_calc_weights_all_levs_itime = l_calc_weights_all_levs .and. .not. l_rad_itime

      ! Calculate thvm for use in prescribe_forcings.
      thvm = calculate_thvm( thlm, rtm, rcm(1,:), exner, thv_ds_zt )

      ! Set large-scale tendencies and subsidence profiles
      call prescribe_forcings( gr, dt_main, um, vm, thlm, & ! In
                               p_in_Pa, exner, rho, rho_zm, thvm, & ! In
                               Frad_SW_up, Frad_SW_down, Frad_LW_down, & ! In
                               rtm, wm_zm, wm_zt, ug, vg, um_ref, vm_ref, & ! Inout
                               thlm_forcing, rtm_forcing, um_forcing, & ! Inout
                               vm_forcing, wprtp_forcing, wpthlp_forcing, & ! Inout
                               rtp2_forcing, thlp2_forcing, rtpthlp_forcing, & ! Inout
                               wpsclrp, sclrm_forcing, edsclrm_forcing, & ! Inout
                               wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, & ! Inout
                               T_sfc, p_sfc, sens_ht, latent_ht, & ! Inout
                               wpsclrp_sfc, wpedsclrp_sfc ) ! Inout

      if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "Fatal error in prescribe_forcings:"
            error stop
        end if
      end if

      ! Add microphysical tendencies to rtm_forcing
      rtm_forcing(:) = rtm_forcing(:) + rcm_mc(:) + rvm_mc(:)

      ! Add radiation and microphysical tendencies to thlm_forcing
      thlm_forcing(:) = thlm_forcing(:) + thlm_mc(:) + radht(:)

      ! Add microphysical tendencies to the forcings for the predictive
      ! variances and covariances.
      wprtp_forcing(:)   = wprtp_forcing(:) + wprtp_mc(:)
      wpthlp_forcing(:)  = wpthlp_forcing(:) + wpthlp_mc(:)
      rtp2_forcing(:)    = rtp2_forcing(:) + rtp2_mc(:)
      thlp2_forcing(:)   = thlp2_forcing(:) + thlp2_mc(:)
      rtpthlp_forcing(:) = rtpthlp_forcing(:) + rtpthlp_mc(:)

      ! Compute total water in ice phase mixing ratio
      rfrzm = zero
      if ( iiri > 0 ) then
        rfrzm = rfrzm + hydromet(:,iiri)
      end if
      if ( iirs > 0 ) then
        rfrzm = rfrzm + hydromet(:,iirs)
      end if
      if ( iirg > 0 ) then
        rfrzm = rfrzm + hydromet(:,iirg)
      end if

      rcm_zm = zt2zm( gr, rcm(1,:) )
      radht_zm = zt2zm( gr, radht )

      ! Add effects of radiation on thlp2
      if ( clubb_config_flags%l_calc_thlp2_rad ) then

        call calculate_thlp2_rad( gr%nz, rcm_zm, thlprcp, radht_zm, & ! intent(in)
                                  params,                           & ! intent(in)
                                  thlp2_forcing )                     ! intent(inout)

      end if

      
      ! Measure time in the beginning part of the main loop
      call cpu_time(time_stop)      
      time_loop_init = time_loop_init + time_stop - time_start      
      call cpu_time(time_start) ! initialize timer for advance_clubb_core
      
      ! Call the parameterization one timestep
      call advance_clubb_core &
           ( gr, l_implemented, dt_main, fcor, sfc_elevation, hydromet_dim, & ! Intent(in)
             thlm_forcing, rtm_forcing, um_forcing, vm_forcing, & ! Intent(in)
             sclrm_forcing, edsclrm_forcing, wprtp_forcing, &     ! Intent(in)
             wpthlp_forcing, rtp2_forcing, thlp2_forcing, &       ! Intent(in)
             rtpthlp_forcing, wm_zm, wm_zt, &                     ! Intent(in)
             wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &         ! Intent(in)
             wpsclrp_sfc, wpedsclrp_sfc,  &                       ! Intent(in)
             rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &         ! Intent(in)
             p_in_Pa, rho_zm, rho, exner, &                       ! Intent(in)
             rho_ds_zm, rho_ds_zt(1,:), invrs_rho_ds_zm, &        ! Intent(in)
             invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, hydromet, &   ! Intent(in)
             rfrzm, radf, wphydrometp, &                          ! Intent(in)
             wp2hmp, rtphmp_zt, thlphmp_zt, &                     ! Intent(in)
             dummy_dx, dummy_dy, &                                ! Intent(in)
             params, nu_vert_res_dep, lmin, &                     ! Intent(in)
             clubb_config_flags, &                                ! Intent(in)
             stats_zt, stats_zm, stats_sfc, &                     ! intent(inout)
             um, vm, upwp, vpwp, up2, vp2, up3, vp3, &            ! Intent(inout)
             thlm, rtm, wprtp, wpthlp, &                          ! Intent(inout)
             wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &       ! Intent(inout)
             sclrm, sclrp2, sclrp3, sclrprtp, sclrpthlp, &        ! Intent(inout)
             wpsclrp, edsclrm, &                                  ! Intent(inout)
             rcm(1,:), cloud_frac, &                              ! Intent(inout)
             wpthvp, wp2thvp, rtpthvp, thlpthvp, &                ! Intent(inout)
             sclrpthvp, &                                         ! Intent(inout)
             uprcp, vprcp, &                                      ! intent(inout)
             pdf_params, pdf_params_zm, &                         ! Intent(inout)
             pdf_implicit_coefs_terms, &                          ! intent(inout)
             Kh_zm, Kh_zt, &                                      ! intent(out)
             thlprcp, wprcp, w_up_in_cloud, ice_supersat_frac, &  ! Intent(out)
             rcm_in_layer, cloud_cover, invrs_tau_zm, &           ! Intent(out)
             err_code_dummy )                                     ! Intent(out)

      if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then
            error stop "Fatal error in clubb, check your parameter values and timestep"
        end if
      end if

      ! Measure time in advance_clubb_core
      call cpu_time(time_stop)
      time_clubb_advance = time_clubb_advance + time_stop - time_start
      call cpu_time(time_start) ! initialize timer for setup_pdf_parameters
      
      wp2_zt = max( zm2zt( gr, wp2 ), w_tol_sqd ) ! Positive definite quantity

      
      if ( .not. trim( microphys_scheme ) == "none" ) then

         !!! Setup the PDF parameters.
         call setup_pdf_parameters_api( gr, gr%nz, pdf_dim, dt_main,                & ! Intent(in)
                        Nc_in_cloud, rcm(1,:), cloud_frac, Kh_zm,                   & ! Intent(in)
                        ice_supersat_frac, hydromet, wphydrometp,                   & ! Intent(in)
                        corr_array_n_cloud, corr_array_n_below,                     & ! Intent(in)
                        pdf_params, l_stats_samp,                                   & ! Intent(in)
                        params,                                                     & ! Intent(in)
                        clubb_config_flags%iiPDF_type,                              & ! Intent(in)
                        l_use_precip_frac,                                          & ! Intent(in)
                        clubb_config_flags%l_predict_upwp_vpwp,                     & ! Intent(in)
                        clubb_config_flags%l_diagnose_correlations,                 & ! Intent(in)
                        clubb_config_flags%l_calc_w_corr,                           & ! Intent(in)
                        clubb_config_flags%l_const_Nc_in_cloud,                     & ! Intent(in)
                        clubb_config_flags%l_fix_w_chi_eta_correlations,            & ! Intent(in)
                        stats_zt, stats_zm, stats_sfc,                           & ! intent(inout)
                        hydrometp2,                                                 & ! Intent(out)
                        mu_x_1_n(1,:,:), mu_x_2_n(1,:,:),                           & ! Intent(out)
                        sigma_x_1_n(1,:,:), sigma_x_2_n(1,:,:),                     & ! Intent(out)
                        corr_array_1_n(1,:,:,:), corr_array_2_n(1,:,:,:),           & ! Intent(out)
                        corr_cholesky_mtx_1(1,:,:,:), corr_cholesky_mtx_2(1,:,:,:), & ! Intent(out)
                        precip_fracs,                                             & ! Intent(inout)
                        hydromet_pdf_params(1,:) )                                  ! Optional(out)

         if ( err_code == clubb_fatal_error ) error stop

         ! Calculate < rt'hm' >, < thl'hm' >, and < w'^2 hm' >.
         call hydrometeor_mixed_moments( gr, gr%nz, pdf_dim, hydromet,               & ! Intent(in)
                                   mu_x_1_n(1,:,:), mu_x_2_n(1,:,:),                 & ! Intent(in)
                                   sigma_x_1_n(1,:,:), sigma_x_2_n(1,:,:),           & ! Intent(in)
                                   corr_array_1_n(1,:,:,:), corr_array_2_n(1,:,:,:), & ! Intent(in)
                                   pdf_params, hydromet_pdf_params(1,:),             & ! Intent(in)
                                   precip_fracs,                                     & ! Intent(in)
                                   stats_zt, stats_zm,                        & ! intent(inout)
                                   rtphmp_zt, thlphmp_zt, wp2hmp )                     ! Intent(out)

      endif ! not microphys_scheme == "none"
      
      ! Measure time in setup_pdf_parameters and hydrometeor_mixed_moments
      call cpu_time(time_stop)
      time_clubb_pdf = time_clubb_pdf + time_stop - time_start
      call cpu_time(time_start) ! initialize timer for SILHS
      
#ifdef SILHS
      !----------------------------------------------------------------
      ! Compute subcolumns if enabled
      !----------------------------------------------------------------

      if ( lh_microphys_type /= lh_microphys_disabled .or. l_silhs_rad ) then

        ! The profile of CLUBB's mixing length, Lscale, is passed into
        ! subroutine generate_silhs_sample for use in calculating the vertical
        ! correlation coefficient.  Rather than output Lscale directly, its
        ! value can be calculated from other fields that are already output from
        ! subroutine advance_clubb_core.  The equation relating Lscale to eddy
        ! diffusivity is:
        !
        ! Kh = c_K * Lscale * sqrt( TKE ).
        !
        ! Kh is available, TKE can be calculated from wp2, up2, and vp2, all of
        ! which are available, and c_K is easily extracted from CLUBB's tunable
        ! parameters.  The equation for Lscale is:
        !
        ! Lscale = Kh / ( c_K * sqrt( TKE ) ).
        !
        ! Since Kh and TKE are output on momentum grid levels, the resulting
        ! calculation of Lscale is also found on momentum levels.  It needs to
        ! be interpolated back to thermodynamic (midpoint) grid levels for
        ! further use.

        ! Calculate TKE
        if ( .not. clubb_config_flags%l_tke_aniso ) then
           ! tke is assumed to be 3/2 of wp2
           em = three_halves * wp2
        else
           em = one_half * ( wp2 + vp2 + up2 )
        endif

        ! Calculate Lscale on momentum levels and then interpolate back to
        ! thermodynamic levels.
        Lscale(1,:) &
        = max( zm2zt( gr, Kh_zm / ( params(ic_K) * sqrt( max( em, em_min ) ) ) ), &
               0.01_core_rknd )
               
        ! Copy grid dzt to variable with column index as 1
        delta_zm(1,:) = gr%dzt
                 
        !$acc data copyout( X_mixt_comp_all_levs, X_nl_all_levs, lh_sample_point_weights, &
        !$acc&              lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, lh_rv_clipped, &
        !$acc&              lh_Nc_clipped ) &
        !$acc& async(1)

        call generate_silhs_sample_api( &
               itime, pdf_dim, lh_num_samples, lh_sequence_length, gr%nz, 1, & ! In
               l_calc_weights_all_levs_itime,                                & ! In
               pdf_params, delta_zm, rcm, Lscale,                            & ! In
               lh_seed,                                                      & ! In
               rho_ds_zt,                                                    & ! In
               mu_x_1_n, mu_x_2_n, sigma_x_1_n, sigma_x_2_n,                 & ! In
               corr_cholesky_mtx_1, corr_cholesky_mtx_2,                     & ! In
               precip_fracs, silhs_config_flags,                             & ! In
               params,                                                       & ! In
               clubb_config_flags%l_uv_nudge,                                & ! In
               clubb_config_flags%l_tke_aniso,                               & ! In
               clubb_config_flags%l_standard_term_ta,                        & ! In
               vert_decorr_coef,                                             & ! In
               stats_lh_zt, stats_lh_sfc,                                    & ! intent(inout)
               X_nl_all_levs, X_mixt_comp_all_levs,                          & ! Out
               lh_sample_point_weights ) ! Out
       
       
        call clip_transform_silhs_output_api( gr, gr%nz, 1, lh_num_samples,       & ! In
                                              pdf_dim, hydromet_dim,          & ! In
                                              X_mixt_comp_all_levs,           & ! In
                                              X_nl_all_levs,                  & ! Inout
                                              pdf_params, l_use_Ncn_to_Nc,    & ! In
                                              lh_rt_clipped, lh_thl_clipped,  & ! Out
                                              lh_rc_clipped, lh_rv_clipped,   & ! Out
                                              lh_Nc_clipped                   ) ! Out
        !$acc end data
        !$acc wait
                                          
        call stats_accumulate_lh &
             ( gr, gr%nz, lh_num_samples, pdf_dim, rho_ds_zt(1,:),      & ! In
               lh_sample_point_weights(1,:,:),  X_nl_all_levs(1,:,:,:), & ! In
               lh_rt_clipped(1,:,:), lh_thl_clipped(1,:,:),             & ! In
               lh_rc_clipped(1,:,:), lh_rv_clipped(1,:,:),              & ! In
               lh_Nc_clipped(1,:,:),                                    & ! In
               stats_lh_zt, stats_lh_sfc )                                ! intent(inout)
          
               

      end if ! lh_microphys_enabled
      
      

#else
      ! Alleviate compiler warnings
      X_nl_all_levs = -999._core_rknd
      X_mixt_comp_all_levs = -999
      lh_sample_point_weights = -999._core_rknd
      if ( .false. .or. Lscale(1,1) < 0._core_rknd ) print *, ""
#endif /* SILHS */
      
      ! Measure time in SILHS
      call cpu_time(time_stop)
      time_SILHS = time_SILHS + time_stop - time_start
      call cpu_time(time_start) ! initialize timer for calc_microphys_scheme_tendcies
      
      !----------------------------------------------------------------
      ! Compute Microphysics
      !----------------------------------------------------------------

      ! Call microphysics scheme and produce microphysics tendencies.
      call calc_microphys_scheme_tendcies( gr, dt_main, time_current, pdf_dim, runtype, & ! In
                              thlm, p_in_Pa, exner, rho, rho_zm, rtm, &               ! In
                              rcm(1,:), cloud_frac, wm_zt, wm_zm, wp2_zt, &           ! In
                              hydromet, Nc_in_cloud, &                                ! In
                              pdf_params, hydromet_pdf_params(1,:), &                 ! In
                              precip_fracs, &                                         ! In
                              X_nl_all_levs(1,:,:,:), X_mixt_comp_all_levs(1,:,:), &  ! In
                              lh_sample_point_weights(1,:,:), &                       ! In
                              mu_x_1_n(1,:,:), mu_x_2_n(1,:,:), &                     ! In
                              sigma_x_1_n(1,:,:), sigma_x_2_n(1,:,:), &               ! In
                              corr_array_1_n(1,:,:,:), corr_array_2_n(1,:,:,:), &     ! In
                              lh_rt_clipped(1,:,:), lh_thl_clipped(1,:,:), &          ! In
                              lh_rc_clipped(1,:,:), lh_rv_clipped(1,:,:), &           ! In
                              lh_Nc_clipped(1,:,:), &                                 ! In
                              silhs_config_flags%l_lh_importance_sampling, &          ! In
                              silhs_config_flags%l_lh_instant_var_covar_src, &        ! In
                              stats_zt, stats_zm, stats_sfc, stats_lh_zt, &  ! intent(inout)
                              Nccnm, &                                                ! Inout
                              hydromet_mc, Ncm_mc, rcm_mc, rvm_mc, &                  ! Out
                              thlm_mc, hydromet_vel_zt, &                             ! Out
                              hydromet_vel_covar_zt_impc, &                           ! Out
                              hydromet_vel_covar_zt_expc, &                           ! Out
                              wprtp_mc, wpthlp_mc, rtp2_mc, &                         ! Out
                              thlp2_mc, rtpthlp_mc )                                  ! Out

      ! Measure time in calc_microphys_scheme_tendcies
      call cpu_time(time_stop)
      time_microphys_scheme = time_microphys_scheme + time_stop - time_start
      call cpu_time(time_start) ! initialize timer for advance_microphys

      ! Calculate Skw_zm for use in advance_microphys.
      ! This field is smoothed by interpolating to thermodynamic levels and then
      ! interpolating back to momentum levels.
      Skw_zm &
      = zt2zm( gr, zm2zt( gr, Skx_func( gr, wp2, zt2zm( gr, wp3 ), w_tol, &
                                        params(iSkw_denom_coef), &
                                        params(iSkw_max_mag) ) ) )

      ! Advance predictive microphysics fields one model timestep.
      call advance_microphys( gr, dt_main, time_current, wm_zt, wp2,      & ! In
                              exner, rho, rho_zm, rcm(1,:),               & ! In
                              cloud_frac, Kh_zm, Skw_zm,                  & ! In
                              rho_ds_zm, rho_ds_zt(1,:), invrs_rho_ds_zt, & ! In
                              hydromet_mc, Ncm_mc, Lscale(1,:),           & ! In
                              hydromet_vel_covar_zt_impc,                 & ! In
                              hydromet_vel_covar_zt_expc,                 & ! In
                              params, nu_vert_res_dep,                    & ! In
                              clubb_config_flags%l_upwind_xm_ma,          & ! In
                              stats_zt, stats_zm, stats_sfc,              & ! intent(inout)
                              hydromet, hydromet_vel_zt, hydrometp2,      & ! Inout
                              K_hm, Ncm, Nc_in_cloud, rvm_mc, thlm_mc,    & ! Inout
                              wphydrometp, wpNcp )                          ! Out

      ! Measure time in calc_microphys_scheme_tendcies
      call cpu_time(time_stop)
      time_microphys_advance = time_microphys_advance + time_stop - time_start
      call cpu_time(time_start) ! initialize timer for the end part of the main loop
      
      if ( clubb_at_least_debug_level( 0 ) ) then
          if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "Fatal error in advance_microphys:"
            error stop   
          endif
      end if

      if ( l_stats_samp ) then
        ! Total microphysical tendency of vapor and cloud water mixing ratios
        call stat_update_var( irvm_mc, rvm_mc, stats_zt )         ! kg/kg/s
        call stat_update_var( ircm_mc, rcm_mc, stats_zt )         ! kg/kg/s
        call stat_update_var( irtm_mc, rvm_mc+rcm_mc, stats_zt )  ! kg/kg/s
        call stat_update_var( ithlm_mc, thlm_mc, stats_zt )       ! K/s
        call stat_update_var( iwprtp_mc, wprtp_mc, stats_zm )     ! m*(kg/kg)/s^2
        call stat_update_var( iwpthlp_mc, wpthlp_mc, stats_zm )   ! K*m/s^2
        call stat_update_var( irtp2_mc, rtp2_mc, stats_zm )       ! (kg/kg)^2/s
        call stat_update_var( ithlp2_mc, thlp2_mc, stats_zm )     ! K^2/s
        call stat_update_var( irtpthlp_mc, rtpthlp_mc, stats_zm ) ! K*(kg/kg)/s
      endif

      ! Radiation is always called on the first timestep in order to ensure
      ! that the simulation is subject to radiative heating and cooling from
      ! the first timestep.
      if ( l_rad_itime ) then

        ! Advance a radiation scheme
        ! With this call ordering, snow and ice water mixing ratio will be
        ! updated by the microphysics, but thlm and rtm will not.  This
        ! somewhat inconsistent, but we would need to move the call to
        ! radiation before the call the microphysics to change this.
        ! -dschanen 17 Aug 2009
        if ( l_silhs_rad ) then

          call silhs_radiation_driver &
               ( gr%nz, lh_num_samples, pdf_dim, hydromet_dim,                          & !In
                 time_current, time_initial, rho, rho_zm,                               & !In
                 p_in_Pa, exner, cloud_frac, ice_supersat_frac, X_nl_all_levs(1,:,:,:), & !In
                 lh_rt_clipped(1,:,:), lh_thl_clipped(1,:,:), lh_rc_clipped(1,:,:),     & !In
                 lh_sample_point_weights(1,:,:), hydromet,                              & !In
                 radht, Frad, Frad_SW_up, Frad_LW_up, Frad_SW_down, Frad_LW_down )        !out

        else

          call advance_clubb_radiation &
               ( gr, time_current, time_initial, rho, rho_zm, p_in_Pa,                    & ! In
                 exner, cloud_frac, ice_supersat_frac, thlm, rtm, rcm(1,:), hydromet, & ! In
                 radht, Frad, Frad_SW_up, Frad_LW_up,                                 & ! Out
                 Frad_SW_down, Frad_LW_down )                                           ! Out

        end if ! l_silhs_rad

      end if ! mod( itime, floor(dt_rad/dt_main) ) == 0 .or. itime == 1

      if ( l_stats_samp ) then
         call stat_update_var( iFrad, Frad, stats_zm )
      endif

      ! Update the radiation variables here so they are updated every timestep
      call update_radiation_variables( gr%nz, radht, Frad_SW_up, Frad_LW_up, &
                                       Frad_SW_down, Frad_LW_down )

      ! End statistics timestep
      call stats_end_timestep( params, &                    ! intent(in)
                               stats_zt, stats_zm, stats_sfc, & ! intent(inout)
                               stats_lh_zt, stats_lh_sfc, & ! intent(inout)
                               stats_rad_zt, stats_rad_zm & ! intent(inout)
#ifdef NETCDF
                               , clubb_config_flags%l_uv_nudge, &
                               clubb_config_flags%l_tke_aniso, &
                               clubb_config_flags%l_standard_term_ta &
#endif
                                )

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
      if ( ( l_stats_last .or. l_stats ) .and. l_stdout ) then
        write(unit=fstdout,fmt='(a,i8,a,f10.1)') 'iteration = ',  & 
          itime, '; time = ', time_current
      end if

      ! Measure time in the end part
      call cpu_time(time_stop)
      time_loop_end = time_loop_end + time_stop - time_start
    
      if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) exit
      end if

    end do ! itime=1, ifinal
    
    ! Measure overall time in the main loop
    call cpu_time(time_stop)
    time_total = time_stop - time_total ! subtract previously saved start time 
    
    ! Check if time budgets sum up
    if ((time_total * (one + timing_tol) < (time_loop_init + time_clubb_advance +  &
      time_clubb_pdf + time_SILHS + time_microphys_scheme +    &
      time_microphys_advance + time_loop_end)) .or.                  &
      (time_total * (one - timing_tol) > (time_loop_init + time_clubb_advance +    &
      time_clubb_pdf + time_SILHS + time_microphys_scheme +    &
      time_microphys_advance + time_loop_end))) then
      
      ! print warning to stderror and stdout 
      write(unit=fstderr, fmt='(a)') 'WARNING! Main loop timing budget has errors &
        & in excess of 1 per cent!'
      write(unit=fstdout, fmt='(a)') 'WARNING! Main loop timing budget has errors &
        & in excess of 1 per cent!'
      
    endif
    
    ! Print timers
    write(unit=fstdout, fmt='(a,f10.4)') 'CLUBB-TIMER time_loop_init =         ', &
      time_loop_init
    write(unit=fstdout, fmt='(a,f10.4)') 'CLUBB-TIMER time_clubb_advance =     ', &
      time_clubb_advance
    write(unit=fstdout, fmt='(a,f10.4)') 'CLUBB-TIMER time_clubb_pdf =         ', &
      time_clubb_pdf
    write(unit=fstdout, fmt='(a,f10.4)') 'CLUBB-TIMER time_SILHS =             ', & 
      time_SILHS
    write(unit=fstdout, fmt='(a,f10.4)') 'CLUBB-TIMER time_microphys_scheme =  ', &
      time_microphys_scheme
    write(unit=fstdout, fmt='(a,f10.4)') 'CLUBB-TIMER time_microphys_advance = ', & 
      time_microphys_advance
    write(unit=fstdout, fmt='(a,f10.4)') 'CLUBB-TIMER time_loop_end =          ', &
      time_loop_end
    write(unit=fstdout, fmt='(a,f10.4)') 'CLUBB-TIMER time_total =             ', &
      time_total
    

!-------------------------------------------------------------------------------
!                       End Main Time Stepping Loop
!-------------------------------------------------------------------------------

    ! Free memory
    call cleanup_clubb( l_input_fields )

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
    deallocate( wp2 )       ! w'^2
    deallocate( wp3 )       ! w'^3
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

    ! Imposed large scale w
    deallocate( wm_zm )       ! momentum levels
    deallocate( wm_zt )       ! thermodynamic levels

    ! Cloud water variables
    deallocate( rcm )
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

    deallocate( uprcp )  ! u'rc'
    deallocate( vprcp )  ! v'rc'

    deallocate( Kh_zt )  ! Eddy diffusivity coefficient: thermo. levels
    deallocate( Kh_zm )  ! Eddy diffusivity coefficient: momentum levels
    deallocate( K_hm ) ! Eddy diff. coef. for hydromets.: mom. levs.

    deallocate( em )
    deallocate( Lscale )
    deallocate( Lscale_up )
    deallocate( Lscale_down )

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

    ! Variables for PDF closure scheme
    deallocate( pdf_params )
    deallocate( pdf_params_zm )
    deallocate( pdf_implicit_coefs_terms )

    deallocate( thlm_mc, rvm_mc, rcm_mc, wprtp_mc, wpthlp_mc, rtp2_mc, &
                thlp2_mc, rtpthlp_mc, hydromet_mc, Ncm_mc, hydromet_vel_zt, &
                hydromet_vel_covar_zt_impc, hydromet_vel_covar_zt_expc )

    deallocate( radf, rcm_zm, radht_zm, X_nl_all_levs, &
                lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, &
                lh_rv_clipped, lh_Nc_clipped, &
                X_mixt_comp_all_levs, lh_sample_point_weights, Nc_in_cloud )

    return

  end subroutine run_clubb

  !-----------------------------------------------------------------------
  subroutine initialize_clubb &
             ( gr, iunit, forcings_file_path, p_sfc, zm_init, &
               l_uv_nudge, &
               l_tke_aniso, &
               thlm, rtm, um, vm, ug, vg, wp2, up2, vp2, rcm, &
               wm_zt, wm_zm, em, exner, &
               thvm, p_in_Pa, &
               rho, rho_zm, rho_ds_zm, rho_ds_zt, &
               invrs_rho_ds_zm, invrs_rho_ds_zt, &
               thv_ds_zm, thv_ds_zt, &
               rtm_ref, thlm_ref, &
               um_ref, vm_ref, &
               Ncm, Nc_in_cloud, Nccnm, &
               sclrm, edsclrm )
    ! Description:
    !   Execute the necessary steps for the initialization of the
    !   CLUBB model run.
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,            & !--------------------------------------------- Constant(s)
        zero,           &
        em_min,         &
        grav,           &
        cm3_per_m3,     &
        cloud_frac_min


    use grid_class, only: grid ! Type

    use parameters_model, only:  & 
        sclr_dim, &
        edsclr_dim

    use parameters_microphys, only: &
        Nc0_in_cloud, & !----------------------------------------------- Variable(s)
        microphys_scheme

    use parameters_radiation, only: radiation_top, rad_scheme !--------- Variable(s)


    use grid_class, only: zm2zt, zt2zm !-------------------------------- Procedure(s)

    use sounding, only: read_sounding !--------------------------------- Procedure(s)

    use time_dependent_input, only: &
      initialize_t_dependent_input, & !-------------------------------- Procedure(s)
      l_t_dependent !-------------------------------------------------- Variable(s)

    use extended_atmosphere_module, only: &
      determine_extended_atmos_bounds !-------------------------------- Procedure(s)

    use mpace_a, only: mpace_a_init !---------------------------------- Procedure(s)

    ! Joshua Fasching
    ! March 2008
    use soil_vegetation, only: & 
      sfc_soil_T_in_K, & !--------------------------------------------- Variable(s)
      deep_soil_T_in_K, &
      veg_T_in_K

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
        up2_vp2_sponge_damp_profile,  &
        initialize_tau_sponge_damp      !------------------------------ Procedure(s)

    use input_names, only: &
      wm_name, &
      omega_name

    use clubb_model_settings, only: &
      extended_atmos_bottom_level, &
      extended_atmos_top_level, &
      extended_atmos_range_size, &
      lin_int_buffer, &
      runtype, &
      dt_main

    use clubb_precision, only: &
      core_rknd !------------------------------------------------------ Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    intrinsic :: min, max, trim, sqrt, size

    ! Input
    integer, intent(in) :: iunit

    character(len=*), intent(in) :: &
      forcings_file_path ! Path to the .dat files containing the forcings

    real( kind = core_rknd ), intent(in) :: &
      p_sfc,  & ! Pressure at the surface        [Pa]
      zm_init   ! Initial moment. level altitude [m]

    logical, intent(in) :: &
      l_uv_nudge, & ! For wind speed nudging
      l_tke_aniso   ! For anisotropic turbulent kinetic energy, i.e. TKE = 1/2 (u'^2 + v'^2 + w'^2)

    ! Output
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) ::  & 
      thlm,            & ! Grid mean of liquid water pot. temp               [K] 
      rtm,             & ! Grid mean of total water mixing ratio             [kg/kg]
      um,              & ! Grid mean of eastward wind                        [m/s]
      vm,              & ! Grid mean of northbound wind                      [m/s]
      ug,              & ! Eastward geostrophic wind                         [m/s] 
      vg,              & ! Northbound geostrophic wind                       [m/s] 
      wp2,             & ! Vertical velocity variance (w'^2)                 [m^2/s^2]
      up2,             & ! East-west velocity variance (u'^2)                [m^2/s^2]
      vp2,             & ! North-south velocity variance (v'^2)              [m^2/s^2]
      rcm,             & ! Cloud water mixing ratio                          [kg/kg]
      wm_zt, wm_zm,    & ! Vertical wind                                     [m/s]
      em,              & ! Turbulence kinetic energy                         [m^2/s^2]
      exner,           & ! Exner function                                    [-] 
      thvm,            & ! Virtual potential temperature                     [K]
      p_in_Pa,         & ! Pressure                                          [Pa]
      rho,             & ! Density (thermodynamic levels)                    [kg/m^3]
      rho_zm,          & ! Density on momentum levels                        [kg/m^3]
      rho_ds_zm,       & ! Dry, static density (moment. levs.)               [kg/m^3]
      rho_ds_zt,       & ! Dry, static density (thermo. levs.)               [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density (m-levs.)                [m^3/kg]
      invrs_rho_ds_zt, & ! Inv. dry, static density (t-levs.)                [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v (m-levs.)                 [K]
      thv_ds_zt,       & ! Dry, base-state theta_v (t-levs.)                 [K]
      um_ref,          & ! Initial profile of u wind                         [m/s]
      vm_ref,          & ! Initial profile of v wind                         [m/s]
      rtm_ref,         & ! Initial profile of rtm                            [kg/kg]
      thlm_ref           ! Initial profile of thlm                           [K]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      Ncm,         & ! Mean cloud droplet conc., <N_c> (thermo. levs.)  [num/kg]
      Nc_in_cloud    ! Mean (in-cloud) cloud droplet concentration      [num/kg]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      Nccnm    ! Cloud condensation nuclei concentration (COAMPS/MG)  [num/kg]

    ! Output
    real( kind = core_rknd ), dimension(gr%nz,sclr_dim), intent(out) ::  & 
      sclrm      ! Standard passive scalar [units vary]

    real( kind = core_rknd ), dimension(gr%nz,edsclr_dim), intent(out) ::  & 
      edsclrm    ! Eddy diffusivity passive scalar [units vary]

    ! Local Variables

    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      p_in_Pa_zm ! Pressure on momentum levels  [Pa]

    real( kind = core_rknd ) ::  &
      rtm_sfc,          & ! Surface total water mixing ratio       [kg/kg]
      thlm_sfc,         & ! Surface liq. water potential temp.     [K]
      cloud_top_height, & ! Cloud top altitude in initial profile  [m]
      em_max                ! Maximum value of initial subgrid TKE   [m^2/s^2]

    character(len=50) :: &
      theta_type, & ! Type of temperature sounding
      alt_type,   & ! Type of altitude sounding
      subs_type     ! Type of large-scale subsidence sounding

    integer :: k ! Loop index

    !---- Begin code ----

    ! Read sounding information
    call read_sounding( gr, iunit, runtype, p_sfc, zm_init, &        ! Intent(in) 
                        thlm, theta_type, rtm, um, vm, ug, vg, & ! Intent(out)
                        alt_type, p_in_Pa, subs_type, wm_zt, &   ! Intent(out)
                        rtm_sfc, thlm_sfc, sclrm, edsclrm )      ! Intent(out)

    ! Covert sounding input to CLUBB compatible input
    call initialize_clubb_variables( gr, alt_type, theta_type,         & ! Intent(in)
                                     p_sfc, rtm_sfc, rtm,          & ! Intent(in)
                                     thlm, p_in_Pa, p_in_Pa_zm,    & ! Intent(inout)
                                     exner, rho, rho_zm,           & ! Intent(out)
                                     rcm, thvm, rho_ds_zm,         & ! Intent(out)
                                     rho_ds_zt, invrs_rho_ds_zm,   & ! Intent(out)
                                     invrs_rho_ds_zt, thv_ds_zm,   & ! Intent(out)
                                     thv_ds_zt, sclrm, edsclrm )     ! Intent(out)

    if ( trim( rad_scheme ) == "bugsrad" ) then
      call determine_extended_atmos_bounds( gr%nz, gr%zt,              & ! Intent(in)
                                          gr%zm, gr%dzt, p_in_Pa_zm, &   ! Intent(in)
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
    do k = 1, gr%nz, 1

       if ( rcm(k) > zero ) then

          ! The initial profile at this level is entirely saturated (due to
          ! constant moisture and temperature over the level).  The level is
          ! entirely cloudy.
          Ncm(k) = Nc_in_cloud(k)

       else ! rcm = 0

          ! The initial profile at this level is entirely unsaturated (due to
          ! constant moisture and temperature over the level).  The level is
          ! entirely clear.
          Ncm(k) = Nc_in_cloud(k) * cloud_frac_min

       endif

    enddo ! k = 1, gr%nz, 1

    select case ( trim( microphys_scheme ) )

    case ( "coamps" )
       ! Initialize Nccnm as in COAMPS-LES
       Nccnm(1:gr%nz) &
       = 30.0_core_rknd &
         * ( one + exp( -gr%zt(1:gr%nz) / 2000.0_core_rknd ) )  &
         * cm3_per_m3 / rho    ! Known magic number
    end select


    ! Initialize imposed w
    select case ( trim( subs_type ) ) ! Perform different operations based off
      !                                   the sounding file
    case ( wm_name )

      wm_zm = zt2zm( gr, wm_zt )
      wm_zm(1) = 0.0_core_rknd
      wm_zm(gr%nz) = 0.0_core_rknd

    case ( omega_name )

      do k=2,gr%nz
        wm_zt(k) = -wm_zt(k) / ( grav*rho(k) )
      end do

      wm_zt(1) = 0.0_core_rknd
      wm_zt(gr%nz) = 0.0_core_rknd

      wm_zm = zt2zm( gr, wm_zt )
      wm_zm(gr%nz) = 0.0_core_rknd

    case default ! This should not happen

      wm_zt = 0.0_core_rknd
      wm_zm = 0.0_core_rknd

    end select

    ! Initialize damping
    if ( thlm_sponge_damp_settings%l_sponge_damping ) then
      call initialize_tau_sponge_damp( gr, dt_main, gr%zt,            & ! Intent(in)
                                       thlm_sponge_damp_settings, & ! Intent(in)
                                       thlm_sponge_damp_profile )   ! Intent(out)
    endif

    if ( rtm_sponge_damp_settings%l_sponge_damping ) then
      call initialize_tau_sponge_damp( gr, dt_main, gr%zt,           & ! Intent(in)
                                       rtm_sponge_damp_settings, & ! Intent(in)
                                       rtm_sponge_damp_profile )   ! Intent(out)
    endif

    if ( uv_sponge_damp_settings%l_sponge_damping ) then
      call initialize_tau_sponge_damp( gr, dt_main, gr%zt,          & ! Intent(in)
                                       uv_sponge_damp_settings, & ! Intent(in)
                                       uv_sponge_damp_profile )   ! Intent(out)
    endif

    if ( wp2_sponge_damp_settings%l_sponge_damping ) then
      call initialize_tau_sponge_damp( gr, dt_main, gr%zm,           & ! Intent(in)
                                       wp2_sponge_damp_settings, & ! Intent(in)
                                       wp2_sponge_damp_profile )   ! Intent(out)
    endif

    if ( wp3_sponge_damp_settings%l_sponge_damping ) then
      call initialize_tau_sponge_damp( gr, dt_main, gr%zt,           & ! Intent(in)
                                       wp3_sponge_damp_settings, & ! Intent(in)
                                       wp3_sponge_damp_profile )   ! Intent(out)
    endif

    if ( up2_vp2_sponge_damp_settings%l_sponge_damping ) then
      call initialize_tau_sponge_damp( gr, dt_main, gr%zm,           & ! Intent(in)
                                   up2_vp2_sponge_damp_settings, & ! Intent(in)
                                   up2_vp2_sponge_damp_profile )   ! Intent(out)
    endif



    ! Initilize Time Dependant Input

    if( l_t_dependent ) then
      call initialize_t_dependent_input &
                   ( iunit, runtype, gr%nz, gr%zt, p_in_Pa )
    end if

    ! Initialize TKE and other fields as needed

    select case ( trim( runtype ) )

      ! Generic case
    case ( "generic" )

      em = 1.0_core_rknd
      em(gr%nz) = em_min

      ! GCSS BOMEX
    case ( "bomex" )

!---> Reduction of initial sounding for stability
!         do k = 1, gr%nz
!            em(k) = 1.0_core_rknd - (gr%zm(k)/3000.0_core_rknd)
!            if ( em(k) < em_min ) then
!               em(k) = em_min
!            end if
!         end do
!         em(1) = em(2)
!         em(gr%nz) = em(gr%nz-1)
!<--- End reduction of initial sounding for stability 24 Jan 07

      em(:) = em_min

      ! GCSS ARM
    case ( "arm" )

!---> Reduction of initial sounding for stability
!         do k = 1, gr%nz
!            if ( gr%zm(k) < 150.0_core_rknd ) then
!               em(k) = ( 0.15_core_rknd * (1.0_core_rknd - gr%zm(k)/150.0_core_rknd) ) / rho_zm(k)
!            else
!               em(k) = em_min
!            end if
!         end do
!         em(1) = em(2)
!         em(gr%nz) = em(gr%nz-1)
!<--- End reduction of initial sounding for stability 24 Jan 07

      em(:) = em_min

      ! June 27 1997 ARM case
    case ( "arm_97" )

      em = 1.0_core_rknd
      em(gr%nz) = em_min

      ! twp_ice
    case ( "twp_ice" )

      em = 1.0_core_rknd
      em(gr%nz) = em_min

      ! March 2000 ARM case
    case ( "arm_0003" )

      em = 1.0_core_rknd
      em(gr%nz) = em_min

      ! 3 year ARM case
    case ( "arm_3year" )

      em = 1.0_core_rknd
      em(gr%nz) = em_min

      ! cloud feedback cases
    case ( "cloud_feedback_s6", "cloud_feedback_s6_p2k",   &
           "cloud_feedback_s11", "cloud_feedback_s11_p2k", &
           "cloud_feedback_s12", "cloud_feedback_s12_p2k" )

      em = 1.0_core_rknd
      em(gr%nz) = em_min

      ! ASTEX_A209 case 16 Jul, 2010 kcwhite
    case ( "astex_a209" )

      cloud_top_height = 700._core_rknd
      em_max = 1.0_core_rknd

      do k=1,gr%nz
        if ( gr%zm(k) < cloud_top_height ) then
          em(k) = em_max
        else
          em(k) = em_min
        end if
      end do
      em(1) = em(2)
      em(gr%nz) = em_min

      ! GCSS FIRE Sc
    case ( "fire" )

      cloud_top_height = 700._core_rknd ! 700 m is the top of the cloud in FIRE
      do k=1,gr%nz
        if ( gr%zm(k) < cloud_top_height ) then
          !em(k) = 1._core_rknd
          em(k) = 4.5_core_rknd
        else
          em(k) = em_min
        end if
      end do
      em(1) = em(2)
      em(gr%nz) = em_min

      ! GCSS ATEX
    case ( "atex" )

      um = max( um, -8._core_rknd ) ! Known magic number (Stevens, et al. 2001, eq. 1)

!---> Reduction of initial sounding for stability
!         do k = 1, gr%nz
!           em(k) = 1.0_core_rknd - (gr%zm(k)/3000.0_core_rknd)
!           if ( em(k) < em_min ) then
!             em(k) = em_min
!           end if
!         end do
!         em(1) = em(2)
!         em(gr%nz) = em(gr%nz-1)
!<--- End reduction of initial sounding for stability 24 Jan 07

      em(:) = em_min

      ! GCSS DYCOMS II RF01
    case ( "dycoms2_rf01" )
      cloud_top_height = 800._core_rknd ! 800 m is the top of the cloud in RF01
      do k=1,gr%nz
        if ( gr%zm(k) < cloud_top_height ) then
          !em(k) = 0.5_core_rknd
          em(k) = 1.1_core_rknd
        else
          em(k) = em_min
        end if
      end do
      em(1) = em(2)
      em(gr%nz) = em_min

      ! GCSS DYCOMS II RF02
    case ( "dycoms2_rf02" )

      em = 1.0_core_rknd
      em(gr%nz) = em_min

      ! Brian for Nov. 11 altocumulus case.
    case ( "nov11_altocu" )

      ! Vince Larson reduced initial forcing.  4 Nov 2005
!          em = 1.0_core_rknd
!          em = 0.1_core_rknd
      ! 4150 + 2800 m is the top of the cloud in Nov11
      cloud_top_height = 2800._core_rknd + gr%zm(1) ! Known magic number
      do k=1,gr%nz
        if ( gr%zm(k) < cloud_top_height ) then

          ! Modification by Adam Smith, 08 April 2008
          ! Reducing the value of em appears to reduce error in the
          ! updated Nov.11 case
          em(k) = 0.01_core_rknd
          ! End of ajsmith4's modification

        else
          em(k) = em_min
        end if
      end do
      em(1) = em(2)
      em(gr%nz) = em_min
      ! End Vince Larson's change.

      ! Adam Smith addition for June 25 altocumulus case.
    case ( "jun25_altocu" )

      ! Vince Larson reduced initial forcing.  4 Nov 2005
!          em = 1.0_core_rknd
!          em = 0.1_core_rknd
!          do k=1,gr%nz
!            if ( gr%zm(k) < 1400._core_rknd ) then
!               em(k) = 0.1_core_rknd
!            else
!               em(k) = em_min
!            end if
!          end do

      ! Note: em_min = 1.0e-6, defined in constants.F
      ! Adam Smith, 28 June 2006
      ! Note: now em_min = 1.5 * w_tol_sqd
      ! Brian Griffin;  Nov. 26, 2008.
      do k = 1, gr%nz
        em(k) = 0.01_core_rknd
      end do

      em(1) = em(2)
      ! End Vince Larson's change.

      em(gr%nz) = em_min

      ! End of ajsmith4's addition

      ! Adam Smith addition for CLEX-9: Nov. 02 altocumulus case.
    case ( "clex9_nov02" )

      ! Vince Larson reduced initial forcing.  4 Nov 2005
!          em = 1.0_core_rknd
!          em = 0.1_core_rknd
      ! 4150 + 1400 m is the top of the cloud in Nov11
      cloud_top_height = 2200._core_rknd + gr%zm(1) ! Known magic number
      do k=1,gr%nz
        if ( gr%zm(k) < cloud_top_height ) then
          em(k) = 0.01_core_rknd
        else
          em(k) = em_min
        end if
      end do
      em(1) = em(2)
      ! End Vince Larson's change.
      em(gr%nz) = em_min

      ! End of ajsmith4's addition

      ! Adam Smith addition for CLEX-9: Oct. 14 altocumulus case.
    case ( "clex9_oct14" )

      ! Vince Larson reduced initial forcing.  4 Nov 2005
!          em = 1.0_core_rknd
!          em = 0.1_core_rknd
      ! 4150 + 1400 m is the top of the cloud in Nov11
      cloud_top_height = 3500._core_rknd + gr%zm(1) ! Known magic number
      do k=1,gr%nz
        if ( gr%zm(k) < cloud_top_height ) then
          em(k) = 0.01_core_rknd
        else
          em(k) = em_min
        end if
      end do
      em(1) = em(2)
      ! End Vince Larson's change.
      em(gr%nz) = em_min

      ! End of ajsmith4's addition

    case ( "lba" )

      em = 0.1_core_rknd
      em(gr%nz) = em_min

      ! Michael Falk for mpace_a Arctic Stratus case.
    case ( "mpace_a" )

      cloud_top_height = 2000._core_rknd
      em_max = 1.0_core_rknd

      do k=1,gr%nz
        if ( gr%zm(k) < cloud_top_height ) then
          em(k) = em_max
        else
          em(k) = em_min
        end if
      end do
      em(1) = em(2)
      em(gr%nz) = em_min

      call mpace_a_init( iunit, forcings_file_path )

      ! Michael Falk for mpace_b Arctic Stratus case.
    case ( "mpace_b" )

      cloud_top_height = 1300._core_rknd ! 1300 m is the cloud top in mpace_b.
                                         ! Michael Falk 17 Aug 2006
      em_max = 1.0_core_rknd

      do k=1,gr%nz

        if ( gr%zm(k) < cloud_top_height ) then
          em(k) = em_max
        else
          em(k) = em_min
        end if
      enddo
      em(1) = em(2)
      em(gr%nz) = em_min

      ! Brian Griffin for COBRA CO2 case.
    case ( "cobra" )

      em = 0.1_core_rknd
      em(gr%nz) = em_min

      ! Michael Falk for RICO tropical cumulus case, 13 Dec 2006
    case ( "rico" )

      cloud_top_height = 1500._core_rknd
      em_max = 1.0_core_rknd
      do k=1,gr%nz
        if ( gr%zm(k) < cloud_top_height ) then
          em(k) = em_max
        else
          em(k) = em_min
        end if
      enddo

      em(1) = em(2)
      em(gr%nz) = em_min

      ! Michael Falk for GABLS2 case, 29 Dec 2006
    case ( "gabls2" )

      cloud_top_height = 800._core_rknd  ! per GABLS2 specifications
      em_max = 0.5_core_rknd
      do k=1,gr%nz
        if ( gr%zm(k) < cloud_top_height ) then
          em(k) = em_max * (1._core_rknd - (gr%zm(k)/cloud_top_height))
        else
          em(k) = em_min
        end if
      end do

      em(1) = em(2)
      em(gr%nz) = em_min

    case ( "gabls3_night" )
      em = 1.0_core_rknd

    case ( "gabls3" )
      em = 1.0_core_rknd
      em(gr%nz) = em_min

      veg_T_in_K = 300._core_rknd
      sfc_soil_T_in_K = 300._core_rknd
      deep_soil_T_in_K = 288.58_core_rknd

    end select

    ! End Initialize TKE and other fields as needed

    !!!! Initialize w'^2 based on initial TKE !!!

    if ( l_tke_aniso ) then

      ! TKE:  em = (1/2) * ( w'^2 + u'^2 + v'^2 )
      ! Evenly divide TKE into its component
      ! contributions (w'^2, u'^2, and v'^2).

      wp2 = (2.0_core_rknd/3.0_core_rknd) * em
      up2 = (2.0_core_rknd/3.0_core_rknd) * em
      vp2 = (2.0_core_rknd/3.0_core_rknd) * em

    else

      ! TKE:  em = (3/2) * w'^2

      wp2 = (2.0_core_rknd/3.0_core_rknd) * em

    end if ! l_tke_aniso

    ! Moved this to be more general -dschanen July 16 2007
    if ( l_uv_nudge .or. uv_sponge_damp_settings%l_sponge_damping ) then
      um_ref = um ! Michael Falk addition for nudging code.  27 Sep/1 Nov 2006
      vm_ref = vm ! ditto
    end if

    if ( thlm_sponge_damp_settings%l_sponge_damping ) then
      thlm_ref = thlm ! Added for nudging code
    end if

    if ( rtm_sponge_damp_settings%l_sponge_damping ) then
      rtm_ref  = rtm
    end if

    return
  end subroutine initialize_clubb
  !-----------------------------------------------------------------------------
  subroutine initialize_clubb_variables( gr, alt_type, theta_type, &
                                         p_sfc, rtm_sfc, rtm, & !thlm_sfc, &
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


    use grid_class, only: grid ! Type

    ! References:
    !   None
    !---------------------------------------------------------------------------

    use grid_class, only: &
        zt2zm ! Procedure(s)


    use constants_clubb, only:  & ! Constant(s)
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

    use parameters_model, only:  & 
        T0,  &  !------------------------------------- Variable(s)
        sclr_dim, &
        edsclr_dim

    use model_flags, only:  &
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
        sat_mixrat_liq, & !--------------------------- Procedure(s)
        rcm_sat_adj

    use array_index, only: &
        iisclr_thl, & !------------------------------- Variable(s)
        iiedsclr_thl

    use clubb_precision, only: &
        core_rknd !----------------------------------- Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    character(len=*), intent(in) :: &
      alt_type,   & ! Type of altitude sounding (altitude or pressure)
      theta_type    ! Type of temperature sounding (temp., theta, or theta_l)

    real( kind = core_rknd ), intent(in) ::  &
      p_sfc,     & ! Surface pressure                              [Pa]
      rtm_sfc !,& ! Surface total water mixing ratio               [kg/kg]
!     thlm_sfc    ! Surface liquid water potential temperature     [K]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      rtm    ! Total water mixing ratio (thermodynamic levels)    [kg/kg]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) ::  &
      thlm,    & ! Liquid water potential temperature (thermo. levs.)  [K] 
      p_in_Pa    ! Pressure (thermodynamic levels)                     [Pa]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) ::  &
      p_in_Pa_zm,      & ! Pressure (momentum levels)                [Pa]
      exner,           & ! Exner function (thermodynamic levels)     [-] 
      rho,             & ! Density (thermodynamic levels)            [kg/m^3]
      rho_zm,          & ! Density on momentum levels                [kg/m^3]
      rcm,             & ! Cloud water mixing ratio (thermo. levs.)  [kg/kg]
      thvm,            & ! Virtual potential temp. (thermo. levs.)   [K]
      rho_ds_zm,       & ! Dry, static density (momentum levels)     [kg/m^3]
      rho_ds_zt,       & ! Dry, static density (thermodynamic levs.) [kg/m^3]
      invrs_rho_ds_zm, & ! Inverse dry, static density (m-levs.)     [m^3/kg]
      invrs_rho_ds_zt, & ! Inverse dry, static density (t-levs.)     [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v (momentum levels) [K]
      thv_ds_zt          ! Dry, base-state theta_v (thermo. levels)  [K]

    real( kind = core_rknd ), dimension(gr%nz,sclr_dim), intent(inout) ::  & 
      sclrm  ! Standard passive scalar           [units vary]

    real( kind = core_rknd ), dimension(gr%nz,edsclr_dim), intent(inout) ::  & 
      edsclrm ! Eddy-diffusivity passive scalar   [units vary]

    ! Local Variables
    real( kind = core_rknd ) ::  &
      pd_sfc, & ! Dry surface pressure                [Pa]
      rv_sfc    ! Surface water vapor mixing ratio    [kg/kg]

    real( kind = core_rknd ), dimension(gr%nz) ::  &
      thm,          & ! Potential temperature (thermodynamic levels)   [K]
      exner_zm,     & ! Exner on momentum levels                       [-]
      th_dry,       & ! Dry potential temperature (thermo. levels)     [K]
      p_dry,        & ! Dry air pressure (thermodynamic levels)        [Pa]
      exner_dry,    & ! Exner of dry air (thermodynamic levels)        [-]
      rho_dry,      & ! Dry air density (thermodynamic levels)         [kg/m^3]
      th_dry_zm,    & ! Dry potential temperature on momentum levels   [K]
      p_dry_zm,     & ! Dry air pressure on momentum levels            [Pa]
      exner_dry_zm, & ! Exner of dry air on momentum levels            [-]
      rho_dry_zm      ! Dry air density on momentum levels             [kg/m^3]

    integer :: k   ! Array index

    ! ---- Begin Code ----

    ! The value of rtm at the surface is output from the sounding, as long as
    ! the initial sounding extendeds to the model surface (at gr%zm(1)).
    if ( rtm_sfc < 0.0_core_rknd ) then
      ! The sounding doesn't extended to the surface, so rtm_sfc is set to a
      ! negative number.  Use rtm(1) as rv_sfc.
      rv_sfc = rtm(1)
    else ! rtm_sfc >= 0.0_core_rknd
      ! The sounding does extended to the surface, so rtm_sfc is the initial value
      ! of total water mixing ratio at the surface.
      rv_sfc = rtm_sfc
    end if

    ! Calculate dry surface pressure from surface pressure and surface water
    ! vapor mixing ratio, such that p_d = p / [ 1 + (R_v/R_d)*r_v ].
    pd_sfc = p_sfc / ( one + ep2 * rv_sfc )

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
          do k = 1, gr%nz, 1
             thlm(k) = thlm(k) / ( p_in_Pa(k) / p0 )**kappa
          enddo

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
    do k = 1, gr%nz, 1
       thvm(k) = thlm(k) * ( one + ep1 * ( rtm(k) / ( one + rtm(k) ) ) )
    enddo

    ! Compute approximate pressure, exner, and density using an approximate
    ! value of theta_v.
    call hydrostatic( gr, thvm, p_sfc,         & ! Intent(in)
                      p_in_Pa, p_in_Pa_zm, & ! Intent(out)
                      exner, exner_zm,     & ! Intent(out)
                      rho, rho_zm          ) ! Intent(out)


    select case( trim( theta_type ) )

    case ( theta_name, temperature_name )

       ! The variable "thlm" actually contains potential temperature (theta)
       ! at this point.
       thm = thlm

       ! Calculate cloud water mixing ratio based on total water mixing ratio
       ! and saturation mixing ratio, which based total pressure and
       ! temperature, which is equal to theta * exner.
       do k = 1, gr%nz
          rcm(k) &
          = max( rtm(k) - sat_mixrat_liq( p_in_Pa(k), thm(k) * exner(k) ), &
                 zero_threshold )
       enddo

       ! Compute initial theta_l based on the theta profile (currently stored
       ! in variable thlm) and cloud water mixing ratio (rcm), such that:
       !  theta_l = theta - [Lv/(Cp*exner)]*rcm.
       do k = 1, gr%nz
          thlm(k) = thlm(k) - Lv/(Cp*exner(k)) * rcm(k)
       enddo

       ! Testing of passive scalars
       if ( iisclr_thl > 0 ) then
          sclrm(:,iisclr_thl) = thlm
       endif
       if ( iiedsclr_thl > 0 ) then
          edsclrm(:,iiedsclr_thl) = thlm
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
       do k = 1, gr%nz, 1
          ! Compute cloud water mixing ratio using an iterative method.
          rcm(k) = rcm_sat_adj( thlm(k), rtm(k), p_in_Pa(k), exner(k) )
       enddo

       ! Compute initial theta.
       do k = 1, gr%nz, 1
          thm(k) = thlm(k) + Lv/(Cp*exner(k)) * rcm(k)
       enddo

       ! Testing of passive scalars
       if ( iisclr_thl > 0 ) then
          sclrm(:,iisclr_thl) = thlm
       endif
       if ( iiedsclr_thl > 0 ) then
          edsclrm(:,iiedsclr_thl) = thlm
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
    thvm = calculate_thvm( thlm, rtm, rcm, exner, &
                           thm * ( one + ep2 * ( rtm - rcm ) )**kappa )

    ! Recompute more accurate initial exner function, pressure, and density
    ! using thvm, which includes the effects of water vapor and cloud water.
    call hydrostatic( gr, thvm, p_sfc,          & ! Intent(in)
                      p_in_Pa, p_in_Pa_zm, &  ! Intent(out)
                      exner, exner_zm,     &  ! Intent(out)
                      rho, rho_zm          )  ! Intent(out)


    !#### Calculate dry, static base-state density for the anelastic ####
    !#### equation set.  Calculate dry pressure from total pressure, ####
    !#### rtm, and rcm, and calculate dry exner from dry pressure.   ####

    !!! Calculate dry density on thermodynamic levels

    do k = 1, gr%nz, 1

       ! Calculate dry pressure from total pressure and water vapor mixing
       ! ratio, such that:  p_d = p / [ 1 + (R_v/R_d)*r_v ].
       p_dry(k) = p_in_Pa(k) / ( one + ep2 * ( rtm(k) - rcm(k) ) )

       ! Calculate dry exner from dry pressure.
       exner_dry(k) = ( p_dry(k) / p0 )**kappa

    enddo

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
    do k = 1, gr%nz, 1
       th_dry(k) = thm(k) * ( one + ep2 * ( rtm(k) - rcm(k) ) )**kappa
    enddo

    ! Compute dry density using dry pressure, dry exner, and theta_d.
    do k = 1, gr%nz, 1
       rho_dry(k) = p_dry(k) / ( Rd * th_dry(k) * exner_dry(k) )
    enddo

    !!! Calculate dry density on momentum levels

    ! Dry pressure at momentum level k = 1 is the dry pressure at the surface.
    p_dry_zm(1) = pd_sfc

    do k = 2, gr%nz, 1
       ! Calculate dry pressure on momentum levels from total pressure (on
       ! momentum levels) and water vapor mixing ratio (interpolated to
       ! momentum levels), such that:  p_d = p / [ 1 + (R_v/R_d)*r_v ].
       p_dry_zm(k) = p_in_Pa_zm(k) &
                     / ( one + ep2 * max( zt2zm( gr, rtm - rcm, k ), &
                                          zero_threshold ) )
    enddo

    do k = 1, gr%nz, 1
       ! Calculate dry exner on momentum levels from dry pressure on momentum
       ! levels.
       exner_dry_zm(k) = ( p_dry_zm(k) / p0 )**kappa
    enddo

    ! Calculate theta_d on momentum levels by interpolating theta and water
    ! vapor mixing ratio to momentum levels.
    do k = 1, gr%nz, 1
       th_dry_zm(k) = zt2zm( gr, thm, k ) &
                      * ( one + ep2 * max( zt2zm( gr, rtm - rcm, k ), &
                                           zero_threshold ) )**kappa
    enddo

    ! Compute dry density on momentum levels using dry pressure on momentum
    ! levels, dry exner on momentum levels, and theta_d interpolated to
    ! momentum levels.
    do k = 1, gr%nz, 1
       rho_dry_zm(k) = p_dry_zm(k) / ( Rd * th_dry_zm(k) * exner_dry_zm(k) )
    enddo

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

  subroutine cleanup_clubb( l_input_fields )

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
        up2_vp2_sponge_damp_profile,  &
        finalize_tau_sponge_damp        !---------------- Procedure(s)

    use time_dependent_input, only: &
        l_t_dependent,    & !---------------------------- Variable(s)
        finalize_t_dependent_input !--------------------- Procedure(s)

    use extended_atmosphere_module, only: &
        finalize_extended_atm

    use advance_clubb_core_module, only: &
        cleanup_clubb_core

    use variables_radiation_module, only: &
        cleanup_radiation_variables

    use microphys_init_cleanup, only: &
        cleanup_microphys

    use corr_varnce_module, only: &
        cleanup_corr_matrix_arrays

    use stats_clubb_utilities, only:  &
        stats_finalize


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

    ! Input Variables
    logical, intent(in) :: &
      l_input_fields

    !----- Begin Code -----

    ! Free memory
    if ( thlm_sponge_damp_settings%l_sponge_damping ) then
       call finalize_tau_sponge_damp( thlm_sponge_damp_profile )
    endif

    if ( rtm_sponge_damp_settings%l_sponge_damping ) then
       call finalize_tau_sponge_damp( rtm_sponge_damp_profile )
    endif

    if ( uv_sponge_damp_settings%l_sponge_damping ) then
       call finalize_tau_sponge_damp( uv_sponge_damp_profile )
    endif

    if ( wp2_sponge_damp_settings%l_sponge_damping ) then
       call finalize_tau_sponge_damp( wp2_sponge_damp_profile )
    endif

    if ( wp3_sponge_damp_settings%l_sponge_damping ) then
       call finalize_tau_sponge_damp( wp3_sponge_damp_profile )
    endif

    if ( up2_vp2_sponge_damp_settings%l_sponge_damping ) then
       call finalize_tau_sponge_damp( up2_vp2_sponge_damp_profile )
    endif

    if( l_t_dependent ) then
       call finalize_t_dependent_input()
    end if

    call finalize_extended_atm( )

    call cleanup_clubb_core( gr )

    call cleanup_radiation_variables( )

    call cleanup_microphys( )

    call cleanup_corr_matrix_arrays( )

    if( l_input_fields ) then
      call cleanup_input_fields()
    end if

    call stats_finalize( stats_zt, stats_zm, stats_sfc, &
                         stats_lh_zt, stats_lh_sfc, &
                         stats_rad_zt, stats_rad_zm )

#ifdef SILHS
    if ( lh_microphys_type /= lh_microphys_disabled ) then
      call latin_hypercube_2D_close_api( )
      call cleanup_latin_hypercube_arrays( )
    end if
#endif


  end subroutine cleanup_clubb

  !-----------------------------------------------------------------------
  subroutine restart_clubb &
             ( iunit, runfile, & ! In
               restart_path_case, time_restart, & ! In
               um, upwp, vm, vpwp, up2, vp2, rtm, & ! Inout
               wprtp, thlm, wpthlp, rtp2, rtp3, & ! Inout
               thlp2, thlp3, rtpthlp, wp2, wp3, & ! Inout
               p_in_Pa, exner, rcm, cloud_frac, & ! Inout
               wpthvp, wp2thvp, rtpthvp, thlpthvp, & ! Inout
               wm_zt, rho, rho_zm, rho_ds_zm, & ! Inout
               rho_ds_zt, thv_ds_zm, thv_ds_zt, & ! Inout
               thlm_forcing, rtm_forcing, wprtp_forcing, & ! Inout
               wpthlp_forcing, rtp2_forcing, & ! Inout
               thlp2_forcing, rtpthlp_forcing, & ! Inout
               hydromet, hydrometp2, wphydrometp, & ! Inout
               Ncm, Nccnm, thvm, em, tau_zm, tau_zt, & ! Inout
               Kh_zt, Kh_zm, ug, vg, Lscale, & ! Inout
               Lscale_up, Lscale_down, thlprcp, & ! Inout
               sigma_sqd_w, sigma_sqd_w_zt, radht, & ! Inout
               pdf_params, pdf_params_zm, & ! Inout
               rcm_mc, rvm_mc, thlm_mc, & ! Out
               wprtp_mc, wpthlp_mc, rtp2_mc, & ! Out
               thlp2_mc, rtpthlp_mc, & ! Out
               wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc ) ! Out

    ! Description:
    !   Execute the necessary steps for the initialization of the
    !   CLUBB model to a designated point in the submitted GrADS file.
    !-----------------------------------------------------------------------
    use inputfields,only:  & 
        input_type, &  ! Variable(s)
        l_input_um, l_input_vm, l_input_rtm, l_input_thlm, & 
        l_input_wp2, l_input_wprtp, l_input_wpthlp,  & 
        l_input_wp3, l_input_rtp2, l_input_rtp3, l_input_thlp2, &
        l_input_thlp3, l_input_rtpthlp, l_input_upwp, l_input_vpwp, & 
        l_input_ug, l_input_vg, l_input_rcm,  & 
        l_input_wm_zt, l_input_exner, l_input_em, & 
        l_input_p, l_input_rho, l_input_rho_zm, &
        l_input_rho_ds_zm, l_input_rho_ds_zt, &
        l_input_thv_ds_zm, l_input_thv_ds_zt, &
        l_input_Lscale, l_input_Lscale_up, l_input_Lscale_down, & 
        l_input_Kh_zt, l_input_Kh_zm, l_input_tau_zm, l_input_tau_zt, & 
        l_input_wpthvp, l_input_wp2thvp, l_input_rtpthvp, l_input_thlpthvp, &
        l_input_radht, &
        l_input_w_1, l_input_w_2, l_input_varnce_w_1, l_input_varnce_w_2, &
        l_input_rt_1, l_input_rt_2, l_input_varnce_rt_1, l_input_varnce_rt_2, &
        l_input_thl_1, l_input_thl_2, l_input_varnce_thl_1, &
        l_input_varnce_thl_2, l_input_mixt_frac, l_input_chi_1, l_input_chi_2, &
        l_input_stdev_chi_1, l_input_stdev_chi_2, l_input_rc_1, l_input_rc_2, &
        l_input_w_1_zm, l_input_w_2_zm, l_input_varnce_w_1_zm, &
        l_input_varnce_w_2_zm, l_input_mixt_frac_zm, &
        l_input_thvm, l_input_rrm, l_input_Nrm,  & 
        l_input_rsm, l_input_rim, l_input_rgm,  & 
        l_input_thlm_forcing, l_input_rtm_forcing, & 
        l_input_up2, l_input_vp2, l_input_sigma_sqd_w, l_input_Ncm,  & 
        l_input_Nccnm, l_input_Nim, l_input_cloud_frac, &
        l_input_sigma_sqd_w_zt, l_input_veg_T_in_K, l_input_deep_soil_T_in_K, &
        l_input_sfc_soil_T_in_K, l_input_wprtp_forcing, &
        l_input_wpthlp_forcing, l_input_rtp2_forcing, l_input_thlp2_forcing, &
        l_input_rtpthlp_forcing, l_input_thlprcp, l_input_rcm_mc, &
        l_input_rvm_mc, l_input_thlm_mc, l_input_wprtp_mc, l_input_wpthlp_mc, &
        l_input_rtp2_mc, l_input_thlp2_mc, l_input_rtpthlp_mc, &
        stat_files

    use inputfields, only: &
        compute_timestep,  & !------------------------------------ Procedure(s)
        stat_fields_reader, &
        set_filenames, &
        get_clubb_variable_interpolated


    use grid_class, only: zt2zm !--------------------------------- Procedure(s)

    use constants_clubb, only: fstderr !-------------------------- Variables(s)

    use pdf_parameter_module, only: pdf_parameter !--------------- Type(s)

    use parameters_model, only: hydromet_dim ! ------------------- Integer

    use clubb_precision, only: time_precision, core_rknd !-------- Variable(s)

    use soil_vegetation, only: &
        l_soil_veg !---------------------------------------------- Variable(s)

    use parameters_microphys, only : &
        microphys_scheme, &  !------------------------------------ Variable(s)
        l_predict_Nc, &
        l_ice_microphys, &
        l_graupel


    implicit none

    ! Input Variables

    integer, intent(in) :: iunit

    character(len=*), intent(in) ::  & 
      runfile,            & ! Filename for the namelist
      restart_path_case     ! Path to GrADS data for restart

    real(kind=time_precision), intent(in) :: & 
      time_restart

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      um,         & ! eastward grid-mean wind component (thermo. levs.)  [m/s]
      upwp,       & ! u'w' (momentum levels)                         [m^2/s^2]
      vm,         & ! northward grid-mean wind component (thermo. levs.) [m/s]
      vpwp,       & ! v'w' (momentum levels)                         [m^2/s^2]
      up2,        & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2,        & ! v'^2 (momentum levels)                         [m^2/s^2]
      rtm,        & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      wprtp,      & ! w' r_t' (momentum levels)                      [kg/kg m/s]
      thlm,       & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      wpthlp,     & ! w'th_l' (momentum levels)                      [(m/s) K]
      rtp2,       & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      rtp3,       & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
      thlp2,      & ! th_l'^2 (momentum levels)                      [K^2]
      thlp3,      & ! th_l'^3 (thermodynamic levels)                 [K^3]
      rtpthlp,    & ! r_t'th_l' (momentum levels)                    [(kg/kg) K]
      wp2,        & ! w'^2 (momentum levels)                         [m^2/s^2]
      wp3,        & ! w'^3 (thermodynamic levels)                    [m^3/s^3]
      p_in_Pa,    & ! Air pressure (thermodynamic levels)            [Pa]
      exner,      & ! Exner function (thermodynamic levels)          [-]
      rcm,        & ! cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
      cloud_frac, & ! cloud fraction (thermodynamic levels)          [-]
      wpthvp,     & ! < w' th_v' > (momentum levels)                 [kg/kg K]
      wp2thvp,    & ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]
      rtpthvp,    & ! < r_t' th_v' > (momentum levels)               [kg/kg K]
      thlpthvp      ! < th_l' th_v' > (momentum levels)              [K^2]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      wm_zt,     & ! vertical mean wind component on thermo. levels  [m/s]
      rho,       & ! Air density on thermodynamic levels             [kg/m^3]
      rho_zm,    & ! Air density on momentum levels                  [kg/m^3]
      rho_ds_zm, & ! Dry, static density on momentum levels          [kg/m^3]
      rho_ds_zt, & ! Dry, static density on thermo. levels           [kg/m^3]
      thv_ds_zm, & ! Dry, base-state theta_v on momentum levels      [K]
      thv_ds_zt    ! Dry, base-state theta_v on thermo levels        [K]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      thlm_forcing,    & ! liquid potential temp. forcing (thermo. levels) [K/s]
      rtm_forcing,     & ! total water forcing (thermo. levels)      [(kg/kg)/s]
      wprtp_forcing,   & ! total water turbulent flux forcing (m-levs) [m*K/s^2]
      wpthlp_forcing,  & ! liq pot temp turb flux forcing (m-levs)[m(kg/kg)/s^2]
      rtp2_forcing,    & ! total water variance forcing (m-levs)   [(kg/kg)^2/s]
      thlp2_forcing,   & ! liq pot temp variance forcing (m-levs)  [K^2/s]
      rtpthlp_forcing    ! <r_t'th_l'> covariance forcing (m-levs) [K*(kg/kg)/s]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(inout) :: &
      hydromet,    & ! Array of hydrometeors                [hm units]
      hydrometp2,  & ! Variance of a hydrometeor (m-levs.)  [<hm units>^2]
      wphydrometp    ! Covariance of w and a hydrometeor    [(m/s) <hm units>]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      Ncm,    & ! Mean cloud droplet concentration, <N_c> (t-levs.)    [num/kg]
      Nccnm,  & ! Cloud condensation nuclei concentration (COAMPS/MG)  [num/kg]
      thvm,   & ! Virtual potential temperature                        [K]
      em,     & ! Turbulent Kinetic Energy (TKE)                       [m^2/s^2]
      tau_zm, & ! Eddy dissipation time scale on momentum levels       [s]
      tau_zt, & ! Eddy dissipation time scale on thermodynamic levels  [s]
      Kh_zt,  & ! Eddy diffusivity coefficient on thermodynamic levels [m^2/s]
      Kh_zm,  & ! Eddy diffusivity coefficient on momentum levels      [m^2/s]
      ug,     & ! u geostrophic wind                                   [m/s]
      vg        ! v geostrophic wind                                   [m/s]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      Lscale,         & ! Length scale                                 [m]
      Lscale_up,      & ! Length scale (upwards component)             [m]
      Lscale_down,    & ! Length scale (downwards component)           [m]
      thlprcp,        & ! thl'rc'                                      [K kg/kg]
      sigma_sqd_w,    & ! PDF width parameter (momentum levels)        [-]
      sigma_sqd_w_zt, & ! PDF width parameter interpolated to t-levs.  [-]
      radht             ! SW + LW heating rate                         [K/s]

    type(pdf_parameter), intent(inout) :: &
      pdf_params,    & ! PDF parameters (thermodynamic levels)    [units vary]
      pdf_params_zm    ! PDF parameters on momentum levels        [units vary]

    ! Output
    real( kind = core_rknd ), intent(out) :: & 
      wpthlp_sfc,      & ! w'theta_l' surface flux   [(m K)/s]
      wprtp_sfc,       & ! w'rt' surface flux        [(m kg)/(kg s)]
      upwp_sfc,        & ! u'w' at surface           [m^2/s^2] 
      vpwp_sfc           ! v'w' at surface           [m^2/s^2]

    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      rcm_mc, &    ! Tendency of liquid water due to microphysics      [kg/kg/s]
      rvm_mc, &    ! Tendency of vapor water due to microphysics       [kg/kg/s]
      thlm_mc, &   ! Tendency of liquid pot. temp. due to microphysics [K/s]
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
    l_input_rtpthvp = .true.
    l_input_thlpthvp = .true.
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
      write(fstderr,*) "Invalid time_restart in "// & 
        "file: "//trim( runfile )
      error stop
    end if

    ! Read data from stats files
    call stat_fields_reader( gr, timestep, & ! In
                             um, upwp, vm, vpwp, up2, vp2, rtm, & ! Inout
                             wprtp, thlm, wpthlp, rtp2, rtp3, & ! Inout
                             thlp2, thlp3, rtpthlp, wp2, wp3, & ! Inout
                             p_in_Pa, exner, rcm, cloud_frac, & ! Inout
                             wpthvp, wp2thvp, rtpthvp, thlpthvp, & ! Inout
                             wm_zt, rho, rho_zm, rho_ds_zm, & ! Inout
                             rho_ds_zt, thv_ds_zm, thv_ds_zt, & ! Inout
                             thlm_forcing, rtm_forcing, wprtp_forcing, & ! Inout
                             wpthlp_forcing, rtp2_forcing, & ! Inout
                             thlp2_forcing, rtpthlp_forcing, & ! Inout
                             hydromet, hydrometp2, wphydrometp, & ! Inout
                             Ncm, Nccnm, thvm, em, tau_zm, tau_zt, & ! Inout
                             Kh_zt, Kh_zm, ug, vg, Lscale, & ! Inout
                             Lscale_up, Lscale_down, thlprcp, & ! Inout
                             sigma_sqd_w, sigma_sqd_w_zt, radht, & ! Inout
                             pdf_params, pdf_params_zm ) ! Inout

      call get_clubb_variable_interpolated &
           ( l_input_rcm_mc, stat_files(1), "rcm_mc", gr%nz, timestep, &
             gr%zt, rcm_mc, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rvm_mc, stat_files(1), "rvm_mc", gr%nz, timestep, &
             gr%zt, rvm_mc, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thlm_mc, stat_files(1), "thlm_mc", gr%nz, timestep, &
             gr%zt, thlm_mc, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wprtp_mc, stat_files(2), "wprtp_mc", gr%nz, timestep, &
             gr%zm, wprtp_mc, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wpthlp_mc, stat_files(2), "wpthlp_mc", gr%nz, timestep, &
             gr%zm, wpthlp_mc, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rtp2_mc, stat_files(2), "rtp2_mc", gr%nz, timestep, &
             gr%zm, rtp2_mc, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thlp2_mc, stat_files(2), "thlp2_mc", gr%nz, timestep, &
             gr%zm, thlp2_mc, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rtpthlp_mc, stat_files(2), "rtpthlp_mc", gr%nz, timestep, &
             gr%zm, rtpthlp_mc, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error


    wpthlp_sfc = wpthlp(1)
    wprtp_sfc  = wprtp(1)
    upwp_sfc   = upwp(1)
    vpwp_sfc   = vpwp(1)

    return
  end subroutine restart_clubb

  !----------------------------------------------------------------------
  subroutine prescribe_forcings( gr, dt, um, vm, thlm, &
                                 p_in_Pa, exner, rho, rho_zm, thvm, &
                                 Frad_SW_up, Frad_SW_down, Frad_LW_down, &
                                 rtm, wm_zm, wm_zt, ug, vg, um_ref, vm_ref, &
                                 thlm_forcing, rtm_forcing, um_forcing, &
                                 vm_forcing, wprtp_forcing, wpthlp_forcing, &
                                 rtp2_forcing, thlp2_forcing, rtpthlp_forcing, &
                                 wpsclrp, sclrm_forcing, edsclrm_forcing, &
                                 wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &
                                 T_sfc, p_sfc, sens_ht, latent_ht, &
                                 wpsclrp_sfc, wpedsclrp_sfc )

    ! Description:
    !   Calculate tendency and surface variables
    ! References:
    !   None
    !----------------------------------------------------------------------

    ! Modules to be included
    use soil_vegetation, only:  &
      l_soil_veg

    use grid_class, only: grid ! Type


    use grid_class, only: zt2zm, zm2zt !---------------------- Procedure(s)

    use stats_variables, only: &
      ish, & !------------------------------------------------ Variable(s)
      ilh, &
      iwpthlp_sfc, &
      iwprtp_sfc, &
      iupwp_sfc, &
      ivpwp_sfc, &
      iustar, &
      isoil_heat_flux, &
      l_stats_samp, &
      iT_sfc

    use stats_type_utilities, only: stat_update_var_pt !------ Procedure(s)

    use constants_clubb, only: & 
      Cp, Lv, kappa, p0, & !---------------------------------- Variable(s)
      zero, fstderr

    use parameters_model, only: &
        sclr_dim,   & !--------------------------------------- Variable(s)
        edsclr_dim

    use clubb_precision, only: core_rknd !-------------------- Variable(s)

    use time_dependent_input, only: &
      apply_time_dependent_forcings, &
      l_t_dependent, &
      l_ignore_forcings

    use soil_vegetation, only: advance_soil_veg, veg_T_in_K

    use array_index, only: & 
        iisclr_rt, iisclr_thl !------------------------------ Variable(s)

    ! Case specific modules
    use arm, only: arm_sfclyr !------------------------------ Procedure(s)

    use arm_0003, only: arm_0003_sfclyr !-------------------- Procedure(s)

    use arm_3year, only: arm_3year_sfclyr !------------------ Procedure(s)

    use astex_a209, only: astex_a209_sfclyr !---------------- astex_a209_tndcy ! Procedure(s)

    use atex, only: atex_tndcy, atex_sfclyr !---------------- Procedure(s)

    use arm_97, only: arm_97_sfclyr !------------------------ Procedure(s)

    use bomex, only: bomex_tndcy, bomex_sfclyr !------------- Procedure(s)

    use clex9_nov02, only: clex9_nov02_read_t_dependent !---- Procedure(s)

    use clex9_oct14, only: clex9_oct14_read_t_dependent !---- Procedure(s)

    use cloud_feedback, only: cloud_feedback_sfclyr !-------- Procedure(s)

    use cobra, only: cobra_sfclyr !-------------------------- Procedure(s)

    use dycoms2_rf01, only:     &           !---------------- Procedure(s)
        dycoms2_rf01_tndcy, dycoms2_rf01_sfclyr

    use dycoms2_rf02, only:  & 
        dycoms2_rf02_tndcy, dycoms2_rf02_sfclyr !------------ Procedure(s)

    use fire, only: &
      fire_sfclyr !------------------------------------------ Procedure(s)

    use gabls2, only: gabls2_tndcy, gabls2_sfclyr !---------- Procedure(s)

    use gabls3, only: gabls3_sfclyr !------------------------ Procedures(s)

    use gabls3_night, only: gabls3_night_sfclyr

    use jun25, only: jun25_altocu_read_t_dependent !--------- Procedure(s)
    
    use lba, only: lba_tndcy, lba_sfclyr !------------------- Procedure(s)

    use mpace_a, only: mpace_a_tndcy, mpace_a_sfclyr !------- Procedure(s)

    use mpace_b, only: mpace_b_tndcy, mpace_b_sfclyr !------- Procedure(s)

    use nov11, only: nov11_altocu_rtm_adjust, nov11_altocu_read_t_dependent ! Procedure(s)

    use rico, only: rico_tndcy, rico_sfclyr !---------------- Procedure(s)

    use twp_ice, only: twp_ice_sfclyr !---------------------- Procedure(s)

    use wangara, only: wangara_tndcy, wangara_sfclyr !------- Procedure(s)

    use sfc_flux, only:  &   !--------------------------- Procedure(s)
      compute_momentum_flux, &
      compute_ubar,          &
      set_sclr_sfc_rtm_thlm

    use clubb_model_settings, only: &
      runtype, & ! Variable(s)
      sfctype, &
      time_current, &
      time_initial

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real(kind=core_rknd), intent(in) :: & 
      dt         ! Model timestep         [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      um,           & ! eastward grid-mean wind component (thermo. levs.)  [m/s]
      vm,           & ! northward grid-mean wind component (thermo. levs.) [m/s]
      thlm,         & ! liq. water pot. temp., th_l (thermo. levels)       [K]
      p_in_Pa,      & ! Air pressure (thermodynamic levels)                [Pa]
      exner,        & ! Exner function (thermodynamic levels)              [-]
      rho,          & ! Air density on thermodynamic levels                [kg/m^3]
      rho_zm,       & ! Air density on momentum levels                     [kg/m^3]
      thvm,         & ! Virtual potential temperature                      [K]
      Frad_SW_up,   & ! SW radiative upwelling flux                        [W/m^2]
      Frad_SW_down, & ! SW radiative downwelling flux                      [W/m^2]
      Frad_LW_down    ! LW radiative downwelling flux                      [W/m^2]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      rtm,             & ! total water mixing ratio, r_t (thermo. levs.) [kg/kg]
      wm_zm,           & ! vertical mean wind comp. on momentum levs     [m/s]
      wm_zt,           & ! vertical mean wind comp. on thermo. levs      [m/s]
      ug,              & ! u geostrophic wind                            [m/s]
      vg,              & ! v geostrophic wind                            [m/s]
      um_ref,          & ! Initial u wind                                [m/s]
      vm_ref             ! Initial v wind                                [m/s]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      thlm_forcing,    & ! liquid potential temp. forcing (thermodynamic levels)[K/s]
      rtm_forcing,     & ! total water forcing (thermodynamic levels)           [(kg/kg)/s]
      um_forcing,      & ! eastward wind forcing (thermodynamic levels)         [m/s/s]
      vm_forcing,      & ! northward wind forcing (thermodynamic levels)        [m/s/s]
      wprtp_forcing,   & ! total water turbulent flux forcing (momentum levels) [m*K/s^2]
      wpthlp_forcing,  & ! liq pot temp turb flux forcing (momentum levels)     [m(kg/kg)/s^2]
      rtp2_forcing,    & ! total water variance forcing (momentum levels)       [(kg/kg)^2/s]
      thlp2_forcing,   & ! liq pot temp variance forcing (momentum levels)      [K^2/s]
      rtpthlp_forcing    ! <r_t'th_l'> covariance forcing (momentum levels)     [K(kg/kg)/s]

    real( kind = core_rknd ), dimension(gr%nz,sclr_dim), intent(inout) :: &
      wpsclrp,       & ! w'sclr' (momentum levels)       [{units vary} m/s]
      sclrm_forcing    ! Passive scalar forcing          [{units vary}/s]

    real( kind = core_rknd ), dimension(gr%nz,edsclr_dim), intent(inout) :: &
      edsclrm_forcing  ! Eddy-diffusion passive scalar forcing    [{units vary}/s]

    real( kind = core_rknd ), intent(inout) :: &
      wpthlp_sfc, & ! w' theta_l' at surface   [(m K)/s]
      wprtp_sfc,  & ! w' r_t' at surface       [(kg m)/( kg s)]
      upwp_sfc,   & ! u'w' at surface          [m^2/s^2]
      vpwp_sfc,   & ! v'w' at surface          [m^2/s^2]
      T_sfc,      & ! surface temperature      [K]
      p_sfc,      & ! surface pressure         [Pa]
      sens_ht,    & ! sensible heat flux       [K m/s]
      latent_ht     ! latent heat flux         [m/s]

    real( kind = core_rknd ), dimension(sclr_dim), intent(inout) :: &
      wpsclrp_sfc      ! Passive scalar flux at surface         [{units vary} m/s]

    real( kind = core_rknd ), dimension(edsclr_dim), intent(inout) :: &
      wpedsclrp_sfc    ! Eddy-diffusion passive scalar flux at surface [{un vary}m/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      wpthep  ! w'theta_e'                [m K/s]

    real( kind = core_rknd ) :: &
      ubar, &     ! mean sfc wind speed   [m/s]
      rho_sfc     ! Density at zm(1)      [kg/m^3]

    real( kind = core_rknd ) :: &
      ustar,          & ! Average value of friction velocity [m/s]
      soil_heat_flux    ! Soil Heat Flux [W/m^2]

    ! Flags to help avoid code duplication
    logical :: &
      l_compute_momentum_flux, &
      l_set_sclr_sfc_rtm_thlm, &
      l_fixed_flux            

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!                    FIND ALL DIAGNOSTIC VARIABLES
!-----------------------------------------------------------------------
    l_compute_momentum_flux = .false.
    l_set_sclr_sfc_rtm_thlm = .false.
    l_fixed_flux            = .false.

    !----------------------------------------------------------------
    ! Set vertical velocity, w, and compute large-scale forcings
    !----------------------------------------------------------------

    ! These lines were added to reset the forcing arrays to 0 each iteration.
    ! This was previously done in the <case>_tndcy subroutine.
    rtm_forcing = zero
    thlm_forcing = zero
    wprtp_forcing = zero
    wpthlp_forcing = zero
    rtp2_forcing = zero
    thlp2_forcing = zero
    rtpthlp_forcing = zero

    if ( l_t_dependent .and. .not. l_ignore_forcings ) then
      ! This should include the following:
      ! "cloud_feedback_s6", "cloud_feedback_s6_p2k",
      !   "cloud_feedback_s11", "cloud_feedback_s11_p2k",
      !   "cloud_feedback_s12", "cloud_feedback_s12_p2k",
      !   "gabls3_night", "arm_97", "gabls3", "twp_ice",
      !   "arm", "arm_0003", "arm_3year", "astex_a209", & "cobra".

      call apply_time_dependent_forcings &
          ( gr, time_current, gr%nz, rtm, rho, exner,& ! In
            thlm_forcing, rtm_forcing, um_ref, vm_ref, um_forcing, vm_forcing, & ! In/Out
            wm_zt, wm_zm, ug, vg, & ! In/Out
            sclrm_forcing, edsclrm_forcing ) ! In/Out

      ! Vince Larson set forcing to zero at the top point so that we don't need
      ! so much sponge damping, which is associated with sawtooth noise
      ! in the cloud_feedback cases.  I don't know how it will affect
      ! the other cases.
      rtm_forcing(gr%nz) = zero
      thlm_forcing(gr%nz) = zero
      ! End Vince Larson's addition

    else ! Legacy method of setting the forcings

      select case ( runtype )

!    case ( "astex_a209" ) ! ASTEX Sc case for K & K
!      call astex_a209_tndcy( wm_zt, wm_zm,  &           ! Intent(out)
!                       thlm_forcing, rtm_forcing , &    ! Intent(out)
!                       sclrm_forcing, edsclrm_forcing ) ! Intent(out)

      case ( "atex" ) ! ATEX case
        call atex_tndcy( gr, time_current, time_initial, &   ! Intent(in)
                         rtm, &                          ! Intent(in)
                         wm_zt, wm_zm, &                 ! Intent(out)
                         thlm_forcing, rtm_forcing, &    ! Intent(out)
                         sclrm_forcing, edsclrm_forcing )! Intent(out)

      case ( "bomex" ) ! BOMEX Cu case
        call bomex_tndcy( gr, rtm, &                           ! Intent(in)
                          thlm_forcing, rtm_forcing, &     ! Intent(out)
                          sclrm_forcing, edsclrm_forcing ) ! Intent(out)

      case ( "dycoms2_rf01" ) ! DYCOMS2 RF01 case
        call dycoms2_rf01_tndcy( gr, thlm_forcing, rtm_forcing,  &    ! Intent(out)
                                 sclrm_forcing, edsclrm_forcing ) ! Intent(out)

      case ( "dycoms2_rf02" ) ! DYCOMS2 RF02 case
        call dycoms2_rf02_tndcy( gr, wm_zt, wm_zm,   &                ! Intent(inout)
                                 thlm_forcing, rtm_forcing, &     ! Intent(out) 
                                 sclrm_forcing, edsclrm_forcing ) ! Intent(out)

      case ( "fire", "generic" ) ! FIRE Sc case
        ! Analytic radiation is computed elsewhere
        thlm_forcing = 0._core_rknd
        rtm_forcing = 0._core_rknd

      case ( "gabls2" ) ! GABLS 2 case
        call gabls2_tndcy( gr, time_current, time_initial,  &   ! Intent(in) 
                           wm_zt, wm_zm, thlm_forcing, &    ! Intent(out)
                           rtm_forcing, &                   ! Intent(out)
                           sclrm_forcing, edsclrm_forcing ) ! Intent(out)

      case ( "lba" )
        call lba_tndcy( gr, thlm_forcing, rtm_forcing, &     ! Intent(out)
                        sclrm_forcing, edsclrm_forcing ) ! Intent(out)

      case ( "mpace_a" ) ! mpace_a arctic stratus case

        call mpace_a_tndcy( gr, time_current, p_in_Pa, &                   ! Intent(in) 
                            wm_zt, wm_zm, thlm_forcing, rtm_forcing, & ! Intent(out)
                            um_ref, vm_ref, &                          ! Intent(out)
                            sclrm_forcing, edsclrm_forcing )           ! Intent(out)

      case ( "mpace_b" ) ! mpace_b arctic stratus case

        call mpace_b_tndcy( gr, p_in_Pa, thvm, &                           ! Intent(in)
                            wm_zt, wm_zm, thlm_forcing, rtm_forcing, & ! Intent(out)
                            sclrm_forcing, edsclrm_forcing )           ! Intent(out)
      case ( "rico" ) ! RICO case
        call rico_tndcy( gr, rtm, exner, &                    ! Intent(in)
                         thlm_forcing, rtm_forcing, &     ! Intent(out)   
                         sclrm_forcing, edsclrm_forcing ) ! Intent(out)

      case ( "wangara" ) ! Wangara dry CBL
        call wangara_tndcy( gr, wm_zt, wm_zm,  &                  ! Intent(out) compute_momentum
                            thlm_forcing, rtm_forcing, &      ! Intent(out)
                            sclrm_forcing, edsclrm_forcing )  ! Intent(out)

      case default

        write(unit=fstderr,fmt=*)  & 
           "prescribe_forcings: Don't know how to handle " &
           //"LS forcing for runtype: "//trim( runtype )
        error stop

      end select

    end if ! l_t_dependent



    !----------------------------------------------------------------
    ! Compute Surface Fluxes
    !----------------------------------------------------------------

    ! Boundary conditions for the second order moments

    ubar = compute_ubar( um(2), vm(2) )

    select case ( trim( runtype ) )

    case ( "rico" )
      l_set_sclr_sfc_rtm_thlm = .true.
      call rico_sfclyr( time_current, um(2), vm(2), thlm(2), rtm(2),  & ! Intent(in)
                          ! 299.8_core_rknd K is the RICO T_sfc;
                          ! 101540 Pa is the sfc pressure.
                          !gr%zt(2), 299.8_core_rknd, 101540._core_rknd,  &           ! Intent(in)
                        gr%zt(2), p_sfc, exner(1), &                    ! Intent(in)
                        upwp_sfc, vpwp_sfc, wpthlp_sfc, &               ! Intent(out)
                        wprtp_sfc, ustar, T_sfc )                       ! Intent(out)

    case ( "gabls3" )
      l_compute_momentum_flux = .true.
      call gabls3_sfclyr( ubar, veg_T_in_K,      &               ! Intent(in)
                          thlm(2), rtm(2), gr%zt(2), exner(1), & ! Intent(in)
                          wpthlp_sfc, wprtp_sfc, ustar )         ! Intent(out)

    case ( "gabls3_night" )
      call gabls3_night_sfclyr( time_current, um(2), vm(2),    & ! Intent(in)
                          thlm(2), rtm(2), gr%zt(2),           & ! Intent(in)
                          upwp_sfc, vpwp_sfc,                  & ! Intent(out)
                          wpthlp_sfc, wprtp_sfc, ustar )         ! Intent(out)
    case ( "jun25_altocu" )
      ! There are no surface momentum or heat fluxes
      ! for the Jun. 25 Altocumulus case.

      ! Ensure ustar is set
      ustar = 0._core_rknd

      ! Read in time dependent inputs
      call jun25_altocu_read_t_dependent( time_current, &       ! Intent(in)
                                          sens_ht, latent_ht )  ! Intent(out)

    case ( "cobra" )
      l_compute_momentum_flux = .true.
      call cobra_sfclyr( time_current, gr%zt(2), rho_zm(1), &      ! Intent(in)
                         thlm(2), ubar,                     &      ! Intent(in)
                         wpthlp_sfc, wprtp_sfc, ustar,      &      ! Intent(out)
                         wpsclrp_sfc, wpedsclrp_sfc, T_sfc )       ! Intent(out)
    case ( "clex9_nov02" )
      ! There are no surface momentum or heat fluxes
      ! for the CLEX-9: Nov. 02 Altocumulus case.

      ! Ensure ustar is set
      ustar = 0._core_rknd

      ! Read in time dependent inputs
      call clex9_nov02_read_t_dependent( time_current, &      ! Intent(in)
                                         sens_ht, latent_ht ) ! Intent(out)

    case ( "clex9_oct14" )
      ! There are no surface momentum or heat fluxes
      ! for the CLEX-9: Oct. 14 Altocumulus case.

      ! Ensure ustar is set.
      ustar = 0._core_rknd

      ! Read in time dependent inputs
      call clex9_oct14_read_t_dependent( time_current, &      ! Intent(in)
                                         sens_ht, latent_ht ) ! Intent(out)

    case ( "astex_a209" )
      l_compute_momentum_flux = .true.
      call astex_a209_sfclyr( time_current, ubar, rtm(2), &     ! Intent(in)
                          thlm(2), gr%zt(2), exner(1), p_sfc, & ! Intent(in)
                          wpthlp_sfc, wprtp_sfc, ustar, T_sfc ) ! Intent(out)

    case ( "nov11_altocu" )
      ! There are no surface momentum or heat fluxes
      ! for the Nov. 11 Altocumulus case.

      ! Ensure ustar is set
      ustar = 0._core_rknd

      ! However, the Nov. 11 Altocumulus case has a one-time adjustment
      ! of rtm at t=3600s after the start of the simulation.
      ! As the nov11_altocu_tndcy subroutine is now obsolete, this was
      ! moved to a separate subroutine, nov11_altocu_rtm_adjust.
      ! This subroutine is called here, as the surface momentum/heat fluxes
      ! are called every timestep.
      ! ~EIHoppe/20110104
      call nov11_altocu_rtm_adjust( gr, time_current, time_initial, dt, & ! (in)
                                    rtm )                             ! (inout)

      ! Read in time dependent inputs
      call nov11_altocu_read_t_dependent( time_current, &             ! Intent(in)
                                          sens_ht, latent_ht )        ! Intent(out)

    case ( "fire", "generic" )  ! Generic setup, and GCSS FIRE
      l_compute_momentum_flux = .true.
      l_set_sclr_sfc_rtm_thlm = .true.
      l_fixed_flux            = .true.
      call fire_sfclyr( time_current, ubar, p_sfc,            &       ! Intent(in)
                        thlm(2), rtm(2), exner(1),   &                ! Intent(in)
                        wpthlp_sfc, wprtp_sfc, ustar, T_sfc )         ! Intent(out)

    case ( "cloud_feedback_s6", "cloud_feedback_s6_p2k",   &
           "cloud_feedback_s11", "cloud_feedback_s11_p2k", &
           "cloud_feedback_s12", "cloud_feedback_s12_p2k", &
           "cgils_s6", "cgils_s6_p2k", "cgils_s11",        &
           "cgils_s11_p2k", "cgils_s12", "cgils_s12_p2k"  ) ! Cloud Feedback cases
      l_compute_momentum_flux = .true.
      l_set_sclr_sfc_rtm_thlm = .true.
      l_fixed_flux            = .true.
      call cloud_feedback_sfclyr( time_current, sfctype, &  ! Intent(in)
                                  thlm(2), rtm(2), gr%zt(2),   &     ! Intent(in)
                                  ubar, p_sfc, T_sfc,            &   ! Intent(in)
                                  wpthlp_sfc, wprtp_sfc, ustar)      ! Intent(out)

    case ( "arm" )
      l_compute_momentum_flux = .true.
      l_set_sclr_sfc_rtm_thlm = .true.
      rho_sfc = 1.1_core_rknd
      call arm_sfclyr( time_current, gr%zt(2), rho_sfc,  &   ! Intent(in)
                        thlm(2), ubar,               &       ! Intent(in)
                        wpthlp_sfc, wprtp_sfc, ustar )       ! Intent(out)

    case ( "arm_0003" )
      l_compute_momentum_flux = .true.
      l_set_sclr_sfc_rtm_thlm = .true.
      call arm_0003_sfclyr( time_current, gr%zt(2), rho_zm(1), &  ! Intent(in)
                            thlm(2), ubar,                     &  ! Intent(in)
                            wpthlp_sfc, wprtp_sfc, ustar )        ! Intent(out)
    case ( "arm_3year" )
      l_compute_momentum_flux = .true.
      l_set_sclr_sfc_rtm_thlm = .true.
      call arm_3year_sfclyr( time_current, gr%zt(2), rho_zm(1), & ! Intent(in)
                              thlm(2), ubar,                    & ! Intent(in)
                              wpthlp_sfc, wprtp_sfc, ustar) ! Intent(out)
    case ( "mc3e" )
      l_compute_momentum_flux = .true.
      l_set_sclr_sfc_rtm_thlm = .true.
      call arm_97_sfclyr( time_current, gr%zt(2), rho_zm(1), &   ! Intent(in)
                          thlm(2), ubar,                     &   ! Intent(in)
                          wpthlp_sfc, wprtp_sfc, ustar )         ! Intent(out)

    case ( "arm_97" )
      l_compute_momentum_flux = .true.
      l_set_sclr_sfc_rtm_thlm = .true.
      call arm_97_sfclyr( time_current, gr%zt(2), rho_zm(1), &   ! Intent(in)
                          thlm(2), ubar,                     &   ! Intent(in)
                          wpthlp_sfc, wprtp_sfc, ustar )         ! Intent(out)


    case ( "atex" )
      l_compute_momentum_flux = .true.
      l_set_sclr_sfc_rtm_thlm = .true.
      call atex_sfclyr( time_current, ubar,          &         ! Intent(in)
                        thlm(2), rtm(2), exner(1),   &         ! Intent(in)
                        wpthlp_sfc, wprtp_sfc, ustar, T_sfc )  ! Intent(out)

    case ( "bomex" )
      l_compute_momentum_flux = .true.
      l_set_sclr_sfc_rtm_thlm = .true.
      call bomex_sfclyr( time_current, rtm(2),                      &  ! Intent(in) 
                         wpthlp_sfc, wprtp_sfc, ustar )                ! Intent(out)

    case ( "dycoms2_rf01" )
      l_compute_momentum_flux = .true.
      l_set_sclr_sfc_rtm_thlm = .true.
      call dycoms2_rf01_sfclyr( time_current, sfctype, p_sfc, &       ! Intent(in)
                                exner(1), ubar,              &        ! Intent(in)
                                thlm(2), rtm(2), rho_zm(1),  &        ! Intent(in)
                                wpthlp_sfc, wprtp_sfc, ustar, T_sfc ) ! Intent(out)
    case ( "dycoms2_rf02" )
      l_compute_momentum_flux = .true.
      l_set_sclr_sfc_rtm_thlm = .true.
      call dycoms2_rf02_sfclyr( time_current, & ! Intent(in)
                                wpthlp_sfc, wprtp_sfc, ustar ) ! Intent(out)

    case ( "gabls2" )
      l_compute_momentum_flux = .true.
      l_set_sclr_sfc_rtm_thlm = .true.
      call gabls2_sfclyr( time_current, time_initial, &              ! Intent(in)
                          gr%zt(2), p_sfc, &                         ! Intent(in)
                          ubar, thlm(2), rtm(2), exner(1), &         ! Intent(in)
                          wpthlp_sfc, wprtp_sfc, ustar, T_sfc )      ! Intent(out)

    case ( "lba" )
      l_compute_momentum_flux = .true.
      l_set_sclr_sfc_rtm_thlm = .true.
      call lba_sfclyr( time_current, time_initial, gr%zt(2), &  ! Intent(in)
                       rho_zm(1), thlm(2), ubar, &              ! Intent(in)
                       wpthlp_sfc, wprtp_sfc, ustar )           ! Intent(out)

    case ( "mpace_a" )
      l_compute_momentum_flux = .true.
      l_set_sclr_sfc_rtm_thlm = .true.
      call mpace_a_sfclyr( time_current, rho_zm(1),      & ! Intent(in)
                            wpthlp_sfc, wprtp_sfc, ustar ) ! Intent(out)
    case ( "mpace_b" )
      l_compute_momentum_flux = .true.
      l_set_sclr_sfc_rtm_thlm = .true.
      call mpace_b_sfclyr( time_current, rho_zm(1), &      ! Intent(in)
                            wpthlp_sfc, wprtp_sfc, ustar ) ! Intent(out)

    case ( "twp_ice" )
      l_compute_momentum_flux = .true.
      l_set_sclr_sfc_rtm_thlm = .true.
      call twp_ice_sfclyr( time_current, gr%zt(2), exner(1), thlm(2), & ! Intent(in)
                            ubar, rtm(2), p_sfc,               &        ! Intent(in)
                            wpthlp_sfc, wprtp_sfc, ustar, T_sfc )       ! Intent(out)

    case ( "wangara" )
      l_compute_momentum_flux = .true.
      l_set_sclr_sfc_rtm_thlm = .true.
      call wangara_sfclyr( time_current, &                 ! Intent(in)
                            wpthlp_sfc, wprtp_sfc, ustar ) ! Intent(out)
    case default
      write(unit=fstderr,fmt=*)  & 
        "Invalid value of runtype = ", runtype
      error stop

    end select ! runtype

    ! These have been placed here to help avoid repetition in the cases
    if( l_compute_momentum_flux ) then
      call compute_momentum_flux( um(2), vm(2), ubar, ustar, &   ! Intent(in)
                                  upwp_sfc, vpwp_sfc )           ! Intent(out)
    end if

    if( l_set_sclr_sfc_rtm_thlm ) then
      call set_sclr_sfc_rtm_thlm( wpthlp_sfc, wprtp_sfc, &      ! Intent(in)
                                  wpsclrp_sfc, wpedsclrp_sfc )  ! Intent(out)
    end if

    ! If the surface type is 0, use fixed fluxes
    if ( sfctype == 0 .and. l_fixed_flux ) then
      wpthlp_sfc = sens_ht
      wprtp_sfc  = latent_ht
      if ( iisclr_thl > 0 ) wpsclrp(:,iisclr_thl) = sens_ht
      if ( iisclr_rt > 0 ) wpsclrp(:,iisclr_rt)   = latent_ht
    end if

    !---------------------------------------------------------------
    ! Compute Surface
    !---------------------------------------------------------------
    if ( l_soil_veg ) then
      wpthep = wpthlp_sfc + (Lv/Cp) * ((p0/p_sfc)**kappa) * wprtp_sfc

      call advance_soil_veg( dt, rho_zm(1), &
                             Frad_SW_down(1) - Frad_SW_up(1), Frad_SW_down(1), &
                             Frad_LW_down(1), wpthep, &
                             stats_sfc, &
                             soil_heat_flux )
    else
      ! Here the value is undefined
      soil_heat_flux = -999._core_rknd
    end if

    ! Store values of surface fluxes for statistics
    if ( l_stats_samp ) then

      call stat_update_var_pt( ish, 1, wpthlp_sfc*rho_zm(1)*Cp,&       ! intent(in)
                               stats_sfc )                             ! intent(inout)

      call stat_update_var_pt( ilh, 1, wprtp_sfc*rho_zm(1)*Lv, &       ! intent(in)
                               stats_sfc )                             ! intent(inout)

      call stat_update_var_pt( iwpthlp_sfc, 1, wpthlp_sfc, &           ! intent(in)
                               stats_sfc )                             ! intent(inout)

      call stat_update_var_pt( iwprtp_sfc, 1, wprtp_sfc, &             ! intent(in)
                               stats_sfc )                             ! intent(inout)

      call stat_update_var_pt( iupwp_sfc, 1, upwp_sfc, &               ! intent(in)
                               stats_sfc )                             ! intent(inout)

      call stat_update_var_pt( ivpwp_sfc, 1, vpwp_sfc, &               ! intent(in)
                               stats_sfc )                             ! intent(inout)

      call stat_update_var_pt( iustar, 1, ustar,  &                    ! intent(in)
                               stats_sfc )                             ! intent(inout)

      call stat_update_var_pt( isoil_heat_flux, 1, soil_heat_flux, &   ! intent(in)
                               stats_sfc )                             ! intent(inout)
      call stat_update_var_pt( iT_sfc, 1, T_sfc, &                     ! intent(in)
                               stats_sfc )                             ! intent(inout)

    endif


    return
  end subroutine prescribe_forcings

!-------------------------------------------------------------------------------
  subroutine advance_clubb_radiation &
             ( gr, time_current, time_initial, rho, rho_zm, p_in_Pa, &
               exner, cloud_frac, ice_supersat_frac, thlm, rtm, rcm, hydromet, &
               radht, Frad, Frad_SW_up, Frad_LW_up, &
               Frad_SW_down, Frad_LW_down )
! Description:
!   Compute a radiation tendency.

! References:
!   None
!-------------------------------------------------------------------------------

    use constants_clubb, only: fstderr  !-------------------------------- Constant(s)


    use grid_class, only: grid ! Type

    use numerical_check, only: is_nan_2d, rad_check !-------------------- Procedure(s)

    use parameters_radiation, only: &
      rad_scheme, & !---------------------------------------------------- Variable(s)
      nparam, &
      l_fix_cos_solar_zen, &
      l_sw_radiation, &
      Fs_values, &
      cos_solar_zen_values, &
      cos_solar_zen_times

    use cos_solar_zen_module, only: cos_solar_zen !--------------------- Procedure(s)

    use simple_rad_module, only: &
      simple_rad, simple_rad_bomex, simple_rad_lba, sunray_sw_wrap

    use error_code, only: &
        clubb_at_least_debug_level,  & !-------------------------------- Procedure
        err_code,                    & !-------------------------------- Error Indicator
        clubb_fatal_error              !-------------------------------- Constant

    use parameters_model, only: hydromet_dim !-------------------------- Variable(s)

    use array_index, only: iirs, iiri !--------------------------------- Variable(s)


    use grid_class, only: zt2zm !--------------------------------------- Procedure

    use interpolation, only: binary_search, lin_interpolate_on_grid !--- Procdure(s)

#ifdef radoffline
    use bugsrad_driver, only: compute_bugsrad_radiation !--------------- Procedure(s)
#endif

    use variables_radiation_module, only: &
      radht_LW, radht_SW, Frad_SW, Frad_LW

    use clubb_precision, only: &
      dp, & !----------------------------------------------------------- double precision
      time_precision, & !----------------------------------------------- Variable(s)
      core_rknd

    use clubb_model_settings, only: &
      day, month, year, & !--------------------------------------------- Variable(s)
      extended_atmos_bottom_level, &
      extended_atmos_top_level, &
      extended_atmos_range_size, &
      lin_int_buffer, &
      lat_vals, &
      lon_vals


    implicit none

    type (grid), target, intent(in) :: gr

    ! External
    intrinsic :: trim

    ! Input Variables
    real(kind=time_precision), intent(in) :: &
      time_current, & ! Current time (UTC)               [s]
      time_initial    ! Start time of model run (UTC)    [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      rho,              & ! Density on thermo. grid                          [kg/m^3]
      rho_zm,           & ! Density on moment. grid                          [kg/m^3]
      p_in_Pa,          & ! Pressure.                                        [Pa] 
      exner,            & ! Exner function.                                  [-]
      cloud_frac,       & ! Cloud fraction (thermodynamic levels)            [-]
      ice_supersat_frac,& ! Ice cloud fraction (thermodynamic levels)        [-]
      thlm,             & ! Liquid potential temperature                     [K]
      rtm,              & ! Total water mixing ratio, r_t (thermo. levels)   [kg/kg]
      rcm                 ! Cloud water mixing ratio, r_c (thermo. levels)   [kg/kg]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species                                         [units vary]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      radht ! Radiative heating rate                                         [K/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      Frad,         & ! Total radiative flux                   [W/m^2]
      Frad_SW_up,   & ! Short-wave upwelling radiative flux    [W/m^2]
      Frad_LW_up,   & ! Long-wave upwelling radiative flux     [W/m^2]
      Frad_SW_down, & ! Short-wave upwelling radiative flux    [W/m^2]
      Frad_LW_down    ! Long-wave upwelling radiative flux     [W/m^2]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      rsm,   & ! Snow mixing ratio                             [kg/kg]
      rim       ! Prisitine ice water mixing ratio             [kg/kg]

    real( kind = core_rknd ) :: Fs0, amu0_core_rknd

    real( kind = dp ) :: amu0 ! Cosine of the solar zenith angle [-]

    integer :: i

    ! ---- Begin Code ----

    ! Initialize all outputs to 0.
    Frad = 0._core_rknd   ! The addition is to prevent an Intel compiler warning of an unused
    Frad_SW_up = 0._core_rknd       ! variable.  May be removed if rho is used below.  -meyern
    Frad_LW_up = 0._core_rknd
    Frad_SW_down = 0._core_rknd
    Frad_LW_down = 0._core_rknd

    radht = 0._core_rknd ! Initialize the radiative heating rate to 0.

    ! If l_fix_cos_solar_zen is not set in the model.in, calculate amu0
    ! Otherwise, it was defined in cos_solar_zen_list file
    if ( l_sw_radiation ) then
      if ( l_fix_cos_solar_zen ) then
        if ( nparam > 1 ) then
          ! Find the closest time value greater than or equal to time_current
          i = binary_search( nparam, cos_solar_zen_times(1:nparam), &
                real( time_current, kind = core_rknd ) )
        else
          i = 1
        end if
        if ( i /= -1 ) then
          amu0 = dble( cos_solar_zen_values(i) )
        else
          write(fstderr,*) "Time not found in cos_solar_zen_times"
          error stop "Critical error."
        end if

      else ! Compute using the formula
        amu0 = cos_solar_zen( day, month, year, time_current, lat_vals, lon_vals )

      end if ! l_fix_cos_solar_zen
    else
      amu0 = 0._dp ! This should disable shortwave radiation

    end if ! l_sw_radiation

    select case ( trim( rad_scheme ) )

    case ( "bugsrad" )
      !----------------------------------------------------------------
      ! BUGSrad Radiation
      !----------------------------------------------------------------
#ifdef radoffline /*This directive is needed for BUGSrad to work with CLUBB.*/

      ! Copy snow and ice
      if ( iirs > 0 ) then
        rsm = hydromet(1:gr%nz,iirs)
      else
        rsm = 0.0_core_rknd
      endif

      if ( iiri > 0 ) then
        rim = hydromet(1:gr%nz,iiri)
      else
        rim = 0.0_core_rknd
      end if

      ! NaN checks added to detect possible errors with BUGSrad
      ! Joshua Fasching November 2007

      if ( clubb_at_least_debug_level( 0 ) ) then

        if ( is_nan_2d( thlm ) ) then
          write(fstderr,*) "thlm before BUGSrad is NaN"
          err_code = clubb_fatal_error
        end if

        if ( is_nan_2d( rcm ) ) then
          write(fstderr,*) "rcm before BUGSrad is NaN"
          err_code = clubb_fatal_error
        end if

        if ( is_nan_2d( rtm ) ) then
          write(fstderr,*) "rtm before BUGSrad is NaN"
          err_code = clubb_fatal_error
        end if

        if ( is_nan_2d( rsm ) ) then
          write(fstderr,*) "rsm before BUGSrad is NaN"
          err_code = clubb_fatal_error
        end if

        if ( is_nan_2d( rim ) ) then
          write(fstderr,*) "rim before BUGSrad is NaN"
          err_code = clubb_fatal_error
        end if

        if ( is_nan_2d( cloud_frac ) ) then
          write(fstderr,*) "cloud_frac before BUGSrad is NaN"
          err_code = clubb_fatal_error
        end if

        if ( is_nan_2d( p_in_Pa ) ) then
          write(fstderr,*) "p_in_Pa before BUGSrad is NaN"
          err_code = clubb_fatal_error
        end if

        if ( is_nan_2d( exner ) ) then
          write(fstderr,*) "exner before BUGSrad is NaN"
          err_code = clubb_fatal_error
        end if

        if ( is_nan_2d( rho_zm ) ) then
          write(fstderr,*) "rho_zm before BUGSrad is NaN"
          err_code = clubb_fatal_error
        end if

        ! Check for impossible negative values
        call rad_check( gr, thlm, rcm, rtm, rim, &               ! Intent(in)
                        cloud_frac, p_in_Pa, exner, rho_zm ) ! Intent(in)

      end if  ! clubb_at_least_debug_level( 0 )

      call compute_bugsrad_radiation &
           ( gr%zm, gr%nz, lin_int_buffer,            &   ! Intent(in)
             extended_atmos_range_size,                 & ! Intent(in)
             extended_atmos_bottom_level,               & ! Intent(in)
             extended_atmos_top_level,                  & ! Intent(in)
             amu0,                                    &   ! Intent(in)
             thlm, rcm, rtm, rsm, rim,           &        ! Intent(in)
             cloud_frac, ice_supersat_frac,           &   ! Intent(in)
             p_in_Pa, zt2zm( gr, p_in_Pa ), exner, rho_zm,&   ! Intent(in)
             radht, Frad,                             &   ! Intent(out)
             Frad_SW_up, Frad_LW_up,                  &   ! Intent(out)
             Frad_SW_down, Frad_LW_down )                 ! Intent(out)

      if ( clubb_at_least_debug_level( 0 ) ) then

        if ( is_nan_2d( Frad ) ) then
          write(fstderr,*) "Frad after BUGSrad is NaN"
          err_code = clubb_fatal_error
        end if

        if ( is_nan_2d( radht ) ) then
          write(fstderr,*) "radht after BUGSrad is NaN"
          err_code = clubb_fatal_error
        end if

      end if  ! clubb_at_least_debug_level( 2 )

#else

      error stop "Cannot call BUGSrad with these compile options."

#endif /*radoffline*/

    case ( "simplified" )
      !----------------------------------------------------------------
      ! Simplified radiation
      !----------------------------------------------------------------

      ! The sunray_sw code cannot handle negative values of cosine
      ! so we check that the value of amu0 is positive here.
      if ( l_sw_radiation .and. amu0 > 0._dp ) then
        amu0_core_rknd = real( amu0, kind = core_rknd )
        if ( nparam > 1 ) then
          call lin_interpolate_on_grid( nparam, cos_solar_zen_values(1:nparam), &
                                    Fs_values(1:nparam), amu0_core_rknd, Fs0 )
        else
          Fs0 = Fs_values(1)
        end if
        call sunray_sw_wrap( gr, Fs0, amu0_core_rknd, rho, rcm, & ! In
                             Frad_SW, radht_SW )              ! Out
      else
        radht_SW = 0._core_rknd
        Frad_SW  = 0._core_rknd
      end if

      call simple_rad( gr, rho, rho_zm, rtm, rcm, exner, & ! In
                       stats_sfc,                        & ! intent(inout)
                       Frad_LW, radht_LW )                 ! Out


      Frad = Frad_SW + Frad_LW
      radht = radht_SW + radht_LW

    case ( "simplified_bomex" )
      !----------------------------------------------------------------
      ! GCSS BOMEX specifiction radiation
      !----------------------------------------------------------------

      call simple_rad_bomex( gr, radht ) ! Out

    case ( "lba"  )
      call simple_rad_lba( gr, time_current, time_initial, & ! In
                           radht )   ! Out

    case ( "none" )
      radht_SW = 0._core_rknd
      Frad_SW  = 0._core_rknd
      radht    = 0._core_rknd

    case default
      write(fstderr,*) "Undefined value for namelist variable rad_scheme: "//trim( rad_scheme )
      error stop "Fatal error encountered in advance_clubb_radiation."

    end select ! Radiation scheme

    if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "Fatal error in advance_clubb_radiation:"
            error stop
        end if
    end if

    return
  end subroutine advance_clubb_radiation


  !-----------------------------------------------------------------------------
  subroutine update_radiation_variables( nz, radht, Frad_SW_up, Frad_LW_up, &
                                         Frad_SW_down, Frad_LW_down )

    ! Description:
    !   Updates the radiation variables using the stat_var_update() subroutine.
    !
    ! References:
    !   None
    !---------------------------------------------------------------------------

    use stats_variables, only: &
      iradht_LW, iradht_SW, iFrad_SW, iFrad_LW, iFrad_SW_up, & !------------------------ Variables
      iFrad_LW_up, iFrad_SW_down, iFrad_LW_down, iT_in_k_rad, ircil_rad, &
      io3l_rad, irsm_rad, ircm_in_cloud_rad, icloud_frac_rad, iice_supersat_frac_rad, &
      iradht_rad, iradht_LW_rad, iFrad_SW_rad, &
      iFrad_LW_rad, iFrad_SW_up_rad, iFrad_LW_up_rad, iFrad_SW_down_rad, &
      iFrad_LW_down_rad, ifdswcl, ifuswcl, ifdlwcl, ifulwcl, iradht, &
      ip_in_mb_rad, isp_humidity_rad

    use variables_radiation_module, only: &
      radht_LW, radht_SW, Frad_SW, Frad_LW, T_in_k, rcil, o3l, & !---------------------- Variables
      rsm_2d, rcm_in_cloud_2d, cloud_frac_2d, ice_supersat_frac_2d, radht_LW_2d, &
      radht_SW_2d, p_in_mb, sp_humidity, Frad_uLW, Frad_dLW, Frad_uSW, Frad_dSW, &
      fdswcl, fuswcl, fdlwcl, fulwcl

    use grid_class, only: &
      flip !------------------------------------------------------------------------- Prodecure(s)

    use stats_variables, only: l_stats_samp, l_output_rad_files !--------------------- Variable(s)

    use stats_type_utilities, only: &
      stat_update_var !----------------------------------------------------------------- Procedure

    use clubb_model_settings, only: &
      extended_atmos_range_size, &
      lin_int_buffer

    use clubb_precision, only: &
      core_rknd

    implicit none

    ! Input Variables

    integer, intent(in) :: nz ! Model domain / # of vertical levels     [-]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      radht,        & ! SW + LW heating rate               [K/s]
      Frad_SW_up,   & ! SW radiative upwelling flux        [W/m^2]
      Frad_LW_up,   & ! LW radiative upwelling flux        [W/m^2]
      Frad_SW_down, & ! SW radiative downwelling flux      [W/m^2]
      Frad_LW_down    ! LW radiative downwelling flux      [W/m^2]

    ! Local Variables

    integer :: rad_zt_dim, rad_zm_dim ! Dimensions of the radiation grid

    ! ---- Begin Code ----

    if ( l_stats_samp ) then

      call stat_update_var( iradht, radht, stats_zt )

      call stat_update_var( iradht_LW, radht_LW, stats_zt )

      call stat_update_var( iradht_SW, radht_SW, stats_zt )

      call stat_update_var( iFrad_SW, Frad_SW, stats_zm )

      call stat_update_var( iFrad_LW, Frad_LW, stats_zm )

      call stat_update_var( iFrad_SW_up, Frad_SW_up, stats_zm )

      call stat_update_var( iFrad_LW_up, Frad_LW_up, stats_zm )

      call stat_update_var( iFrad_SW_down, Frad_SW_down, stats_zm )

      call stat_update_var( iFrad_LW_down, Frad_LW_down, stats_zm )

      if ( l_output_rad_files ) then

        rad_zt_dim = (nz-1)+lin_int_buffer+extended_atmos_range_size
        rad_zm_dim = (nz-1)+lin_int_buffer+extended_atmos_range_size+1

        call stat_update_var( iT_in_K_rad, real( flip(T_in_K(1,:), rad_zt_dim),&
                 kind = core_rknd ), stats_rad_zt )

        call stat_update_var( ircil_rad, real( flip(rcil(1,:), rad_zt_dim), &
                 kind = core_rknd  ), stats_rad_zt )

        call stat_update_var( io3l_rad, real( flip(o3l(1,:), rad_zt_dim), &
                 kind = core_rknd  ), stats_rad_zt )

        call stat_update_var( irsm_rad, real( flip(rsm_2d(1,:), rad_zt_dim), &
                 kind = core_rknd  ), stats_rad_zt )

        call stat_update_var( ircm_in_cloud_rad, &
          real( flip(rcm_in_cloud_2d(1,:), rad_zt_dim), kind = core_rknd  ), stats_rad_zt )

        call stat_update_var( icloud_frac_rad, &
          real( flip(cloud_frac_2d(1,:), rad_zt_dim), kind = core_rknd  ), stats_rad_zt )
        
        call stat_update_var( iice_supersat_frac_rad, &
          real( flip(ice_supersat_frac_2d(1,:), rad_zt_dim), kind = core_rknd  ), stats_rad_zt )

        call stat_update_var( iradht_rad, real(flip((radht_SW_2d(1,:) + &
               radht_LW_2d(1,:)), rad_zt_dim), kind = core_rknd  ), stats_rad_zt )

        call stat_update_var( iradht_LW_rad, &
          real( flip(radht_LW_2d(1,:), rad_zt_dim), kind = core_rknd  ),stats_rad_zt )

        call stat_update_var( ip_in_mb_rad, &
          real( flip(p_in_mb(1,:), rad_zt_dim), kind = core_rknd  ), stats_rad_zt )

        call stat_update_var( isp_humidity_rad, &
          real( flip(sp_humidity(1,:), rad_zt_dim), kind = core_rknd  ), stats_rad_zt )

        call stat_update_var( iFrad_SW_rad, real( flip((Frad_uSW(1,:) - &
             Frad_dSW(1,:)), rad_zm_dim), kind = core_rknd  ), stats_rad_zm )

        call stat_update_var( iFrad_LW_rad, real( flip((Frad_uLW(1,:) - &
             Frad_dLW(1,:)), rad_zm_dim), kind = core_rknd ), stats_rad_zm )

        call stat_update_var( iFrad_SW_up_rad, &
             real( flip(Frad_uSW(1,:), rad_zm_dim), kind = core_rknd ), stats_rad_zm )

        call stat_update_var( iFrad_LW_up_rad, real( flip(Frad_uLW(1,:), &
             rad_zm_dim), kind = core_rknd ), stats_rad_zm )

        call stat_update_var( iFrad_SW_down_rad, real( flip(Frad_dSW(1,:), &
             rad_zm_dim), kind = core_rknd ), stats_rad_zm )

        call stat_update_var( iFrad_LW_down_rad, real( flip(Frad_dLW(1,:), &
             rad_zm_dim), kind = core_rknd ), stats_rad_zm )

        call stat_update_var( ifdswcl, real( flip(fdswcl(1,:), rad_zm_dim), &
             kind = core_rknd ), stats_rad_zm )
        call stat_update_var( ifuswcl, real( flip(fuswcl(1,:), rad_zm_dim), &
             kind = core_rknd ), stats_rad_zm )
        call stat_update_var( ifdlwcl, real( flip(fdlwcl(1,:), rad_zm_dim), &
             kind = core_rknd ), stats_rad_zm )
        call stat_update_var( ifulwcl, real( flip(fulwcl(1,:), rad_zm_dim), &
             kind = core_rknd ), stats_rad_zm )

      end if ! l_output_rad_files

    end if ! lstats_samp

  end subroutine update_radiation_variables

  !-----------------------------------------------------------------------
  subroutine silhs_radiation_driver &
             ( nz, lh_num_samples, pdf_dim, hydromet_dim, time_current, &
               time_initial, rho, rho_zm, p_in_Pa, exner, &
               cloud_frac, ice_supersat_frac, X_nl_all_levs, &
               lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, &
               lh_sample_point_weights, hydromet, &
               radht, Frad, Frad_SW_up, Frad_LW_up, Frad_SW_down, Frad_LW_down )

  ! Description:
  !   Computes radiation over a set of sample points and averages the
  !   results

  ! References
  !   clubb:ticket:663
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd

    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constant

    use latin_hypercube_driver_module, only: &
      copy_X_nl_into_hydromet_all_pts   !--------------------- Procedure

    use constants_clubb, only: &
      fstderr        !------------------------------------------- Constant

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz, &                 ! Number of vertical levels
      lh_num_samples, &     ! Number of SILHS sample points
      pdf_dim, &            ! Number of lognormal variates
      hydromet_dim          ! Number of hydrometeor species

    real( kind = time_precision ), intent(in) :: &
      time_current, & ! Current time of simulation               [s]
      time_initial    ! Start time of simulation                 [s]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rho,               & ! Density on thermo. grid                   [kg/m^3]
      rho_zm,            & ! Density on moment. grid                   [kg/m^3]
      p_in_Pa,           & ! Pressure.                                 [Pa] 
      exner,             & ! Exner function.                           [-]
      cloud_frac,        & ! Cloud fraction (thermodynamic levels)     [-]
      ice_supersat_frac    ! Ice cloud fraction (thermodynamic levels) [-]

    real( kind = core_rknd ), dimension(lh_num_samples,nz,pdf_dim), intent(in) :: &
      X_nl_all_levs        ! Normal-lognormal samples                  [units vary]

    real( kind = core_rknd ), dimension(lh_num_samples,nz), intent(in) :: &
      lh_rt_clipped,  & ! rt generated from silhs sample points
      lh_thl_clipped, & ! thl generated from silhs sample points
      lh_rc_clipped     ! rc generated from silhs sample points

    real( kind = core_rknd ), dimension(lh_num_samples,nz), intent(in) :: &
      lh_sample_point_weights ! Weight of each SILHS sample point      [-]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet             ! Hydrometeor mean fields

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      radht,        & ! Radiative heating rate                         [K/s]
      Frad,         & ! Total radiative flux                           [W/m^2]
      Frad_SW_up,   & ! Short-wave upwelling radiative flux            [W/m^2]
      Frad_LW_up,   & ! Long-wave upwelling radiative flux             [W/m^2]
      Frad_SW_down, & ! Short-wave downwelling radiative flux          [W/m^2]
      Frad_LW_down    ! Long-wave downwelling radiative flux           [W/m^2]

    ! Local Variables
    real( kind = core_rknd ), dimension(lh_num_samples,nz,hydromet_dim) :: &
      hydromet_all_pts ! SILHS sample of hydrometeors for each column  [units vary]

    real( kind = core_rknd ), dimension(lh_num_samples,nz) :: &
      Ncn_all_points   ! SILHS sample of Ncn for each column           [#/kg]
                       ! (not used)

    real( kind = core_rknd ), dimension(lh_num_samples,nz) :: &
      radht_samples,        &    ! radht evaluated at each sample point
      Frad_samples,         &    ! Frad evaluated at each sample point
      Frad_SW_up_samples,   &    ! Frad_SW_up evaluated at each sample point
      Frad_LW_up_samples,   &    ! Frad_LW_up evaluated at each sample point
      Frad_SW_down_samples, &    ! Frad_SW_down evaluated at each sample point
      Frad_LW_down_samples       ! Frad_LW_down evaluated at each sample point

    integer :: isample, k ! Looping variates

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    call copy_X_nl_into_hydromet_all_pts &
         ( nz, pdf_dim, lh_num_samples, &                 ! Intent(in)
           X_nl_all_levs, &                               ! Intent(in)
           hydromet, &                                    ! Intent(in)
           hydromet_all_pts, &                            ! Intent(out)
           Ncn_all_points )                               ! Intent(out)

    do isample=1, lh_num_samples
      ! Call a radiation scheme
      call advance_clubb_radiation &
           ( gr, time_current, time_initial, rho, rho_zm, p_in_Pa, &           ! Intent(in)
             exner, cloud_frac, ice_supersat_frac, lh_thl_clipped(isample,:), & ! Intent(in)
             lh_rt_clipped(isample,:), lh_rc_clipped(isample,:), &         ! Intent(in)
             hydromet_all_pts(isample,:,:), &                                        ! Intent(in)
             radht_samples(isample,:), Frad_samples(isample,:), &                    ! Intent(out)
             Frad_SW_up_samples(isample,:), Frad_LW_up_samples(isample,:), &         ! Intent(out)
             Frad_SW_down_samples(isample,:), Frad_LW_down_samples(isample,:) )      ! Intent(out)
    end do

    ! Average results
    forall ( k = 1:nz )

      radht(k) = sum( radht_samples(:,k) * lh_sample_point_weights(:,k) ) / &
                  real( lh_num_samples, kind=core_rknd )
      Frad(k)  = sum( Frad_samples(:,k) * lh_sample_point_weights(:,k) ) / &
                  real( lh_num_samples, kind=core_rknd )
      Frad_SW_up(k) = sum( Frad_SW_up_samples(:,k) * lh_sample_point_weights(:,k) ) / &
                       real( lh_num_samples, kind=core_rknd )
      Frad_LW_up(k) = sum( Frad_LW_up_samples(:,k) * lh_sample_point_weights(:,k) ) / &
                       real( lh_num_samples, kind=core_rknd )
      Frad_SW_down(k)  = sum( Frad_SW_down_samples(:,k) * lh_sample_point_weights(:,k) ) / &
                          real( lh_num_samples, kind=core_rknd )
      Frad_LW_down(k)  = sum( Frad_LW_down_samples(:,k) * lh_sample_point_weights(:,k) ) / &
                          real( lh_num_samples, kind=core_rknd )

    end forall

    if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then
          write(fstderr,*) "Fatal error in silhs_radiation_driver:"
          error stop
        end if
    end if

    return
  end subroutine silhs_radiation_driver
  !-----------------------------------------------------------------------

end module clubb_driver



