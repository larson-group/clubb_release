// Standalone CLUBB lifecycle and case-owned model state.
//
// This file follows src/clubb_driver.F90. run_clubb performs a complete case;
// init_clubb_case, set_case_initial_conditions, advance_clubb_to_end and
// clean_up_clubb expose the same stages for reruns and bounded advances.
// get_stats returns a read-only statistics snapshot, get_runtime_batch_config
// returns total/runtime column counts, and get_clubb_params_all copies the
// full parameter matrix. The three private scientific helpers retain their
// native names and argument order.
//
// Initialize Kokkos before constructing a driver. Arrays are column-first
// LayoutLeft Views. Initialization uses host fields; timestepping retains
// case-owned execution fields across calls. The host microphysics engine and
// sampled statistics own explicit transfers. Only optional BFB PDF matrices
// retain src/f2k calls.
// Named timing regions use the C++ code_timer_module service.
// Mutable physics/input state belongs to this driver; no Native model module
// or saved bridge copy holds the C++ case state.
// Storage indices are zero-based; timestep and batch numbers retain Native ordinals.
// Shared scalar/species metadata keeps its Native contract at explicit access boundaries.
//
#include "scratch_workspace.hpp"
#include "CLUBB_core/error_code.hpp"
#include "CLUBB_core/code_timer_module.hpp"
#include "CLUBB_core/numerical_check.hpp"
#include "CLUBB_core/calc_pressure.hpp"
#include "clubb_driver.hpp"
#include "CLUBB_core/grid_adaptation_module.hpp"
#include "Microphys/pdf_hydromet_microphys_wrapper.hpp"
#include "Microphys/microphys_init_cleanup.hpp"
#include "Microphys/microphys_driver.hpp"
#include "Microphys/advance_microphys_module.hpp"
#include "Microphys/cloud_sed_module.hpp"
#include "SILHS/latin_hypercube_driver_module.hpp"
#include "Input_fields/hydrostatic_module.hpp"
#include "Radiation/radiation_module.hpp"
#include "Input_fields/sounding.hpp"
#include "Input_fields/input_names.hpp"
#include "Input_fields/input_fields.hpp"
#include "Benchmark_cases/prescribe_forcings.hpp"
#include "Benchmark_cases/time_dependent_input.hpp"
#include "Benchmark_cases/mpace_a.hpp"
#include "readiopdata_module.hpp"

#include "CLUBB_core/saturation.hpp"
#include "CLUBB_core/clubb_api_module.hpp"
#include "CLUBB_core/interpolation.hpp"
#include "driver_arrays.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace clubb_kokkos {
ClubbDriver::ClubbDriver() {
  // Allocate the three fixed soil profiles; package settings belong to this case.
  deep_soil_T_in_K_init = DriverArray1("deep soil initial", 1);
  sfc_soil_T_in_K_init  = DriverArray1("surface soil initial", 1);
  veg_T_in_K_init       = DriverArray1("vegetation initial", 1);
}

void ClubbDriver::run_clubb(const std::string& namelist_filename,   // In
                            bool l_stdout,                          // In
                            err_info_type& err_info,                // Out
                            const DriverArray2* clubb_params_in,    // Optional In
                            const DriverBool1* model_flags_array) { // Optional In

  // Run initialization, each parameter batch, and cleanup in native order.
  // The namelist and optional parameter/flag overrides define the case. err_info
  // returns status from all stages; stats writes each batch into its output slice.
  int num_batches{}; // Number of reruns needed to cover all parameter columns [-]
  int batch_idx{};   // Current runtime batch index [-]

  // Build the full model state from the namelist.  init_clubb_case owns the
  // expensive one-time setup: namelist reads, grid construction, allocation,
  // tunable-parameter loading, stats file creation, radiation/microphysics
  // initialization. Runtime execution allocations survive every batch reset.
  // model_flags_array is optional because normal standalone runs use the flags
  // from the namelist/defaults, while tuner-style callers can inject a flag set
  // without editing input files.
  if (model_flags_array) {
    init_clubb_case(namelist_filename, err_info, clubb_params_in, model_flags_array);
  } else {
    init_clubb_case(namelist_filename, err_info, clubb_params_in, nullptr);
  }

  if (any_equal(err_info.err_code, clubb_fatal_error)) {
    return;
  }

  // init_clubb_case validates that total_param_sets is evenly divisible by
  // batch_size.  A no-batching run is therefore just one batch whose width is
  // total_param_sets.  A batched run reuses the same allocated runtime arrays
  // for each slice of clubb_params_all.
  num_batches = ((total_param_sets / batch_size));

  for (batch_idx = 1; batch_idx < num_batches + 1; batch_idx += 1) {
    // Reset the model to the case initial condition for this logical batch.
    // This call also selects the requested columns from clubb_params_all into
    // clubb_params and retargets stats output to the matching column slice.
    set_case_initial_conditions(err_info, nullptr, batch_idx);

    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      break;
    }

    // Run this batch to completion, advancing all fields from time_initial to
    // time_final.  advance_clubb_to_end contains the timestep loop and is where
    // the core model, diagnostics, physics packages, and stats output are called.
    advance_clubb_to_end(l_stdout, err_info, false, std::nullopt, std::nullopt);

    if (clubb_at_least_debug_level_api(0)) {

      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        break;
      }
    }
  }

  // Release arrays, stats handles, physics package state, and any other
  // persistent module-level resources.  This runs after successful completion
  // and after fatal batch exits.
  clean_up_clubb(err_info);
  return;
}

void ClubbDriver::get_stats(stats_type& stats_out) const { // Out

  // Return native values without transferring ownership of the live output file.
  copy_stats_snapshot(stats, stats_out);
}

void ClubbDriver::get_runtime_batch_config(int& total_param_sets_out,   // Out
                                           int& batch_size_out) const { // Out

  // Return total requested columns and the number advanced together.
  total_param_sets_out = total_param_sets;
  batch_size_out       = batch_size;
}

void ClubbDriver::get_clubb_params_all(DriverArray2& clubb_params_all_out) const { // Out

  // Copy the full (parameter-set column, parameter) matrix into independent storage.
  clubb_params_all_out =
      DriverArray2("clubb_params_all_out", clubb_params_all.extent(0), clubb_params_all.extent(1));
  Kokkos::deep_copy(clubb_params_all_out, clubb_params_all);
}

void ClubbDriver::init_clubb_case(const std::string& runfile, // In
                                  err_info_type& err_info,    // Out
                                  const DriverArray2* clubb_params_in,
                                  const DriverBool1* model_flags_array) { // Optional In

  // Read case settings, construct grids, allocate state, initialize the sounding and stats.
  // Native INTENT(OUT) releases the caller's previous allocatable components.
  Kokkos::fence("clubb_driver::reinitialize");
  workspace.clear();
  auto workspace_scope         = workspace.scope();
  err_info                     = err_info_type{};
  execution                    = {};
  execution_fields_initialized = false;

  auto& l_t_dependent             = time_dependent_input.l_t_dependent;
  auto& l_input_xpwp_sfc          = time_dependent_input.l_input_xpwp_sfc;
  auto& l_ignore_forcings         = time_dependent_input.l_ignore_forcings;
  auto& l_sfc_already_initialized = time_dependent_input.l_sfc_already_initialized;
  auto& rad_scheme                = parameters_radiation.rad_scheme;
  auto& l_soil_veg                = soil_vegetation.l_soil_veg;
  auto& microphys_scheme          = microphys.microphys_scheme;
  auto& lh_microphys_type         = microphys.lh_microphys_type;
  auto& lh_num_samples            = microphys.lh_num_samples;

  // Local variables.
  // Selected option for the two-component normal (double Gaussian) PDF type to use for the w, rt,
  // and theta-l (or w, chi, and eta) portion of CLUBB's multivariate, two-component PDF.
  int iiPDF_type{};
  int ipdf_call_placement{};  // Selected option for the placement of the call to CLUBB's PDF.
  int penta_solve_method{};   // Option to set the penta-diagonal matrix solving method
  int tridiag_solve_method{}; // Option to set the tri-diagonal matrix solving method
  int saturation_formula{};   // Integer that stores the saturation formula to be used
  int Lscale_method{};        // Mixing-length formulation

  // Integer that stores what remapping technique should be used to remap the values from one grid
  // to another (starts with 1, so 0 is invalid value for flag)
  int grid_remap_method{};

  // Integer that stores how the grid should be adapted every timestep or if the grid should not be
  // adapted at all
  int grid_adapt_in_time_method{};
  int fill_holes_type{}; // Hole-filling option for the fill_holes_vertical procedure.
  std::string default_stats_output_filename{"__CLUBB_DEFAULT_STATS_OUTPUT_FILENAME__"};
  std::string stats_registry{};
  std::string stats_output_filename{};
  std::string stats_output_path{};
  int noutdir{};

  // Flag to use precipitation fraction in KK microphysics. The precipitation fraction is
  // automatically set to 1 when this flag is turned off.
  bool l_use_precip_frac{};

  // Flag to predict <u'w'> and <v'w'> along with <u> and <v> alongside the advancement of <rt>,
  // <w'rt'>, <thl>, <wpthlp>, <sclr>, and <w'sclr'> in subroutine advance_xm_wpxp.  Otherwise,
  // <u'w'> and <v'w'> are still approximated by eddy diffusivity when <u> and <v> are advanced in
  // subroutine advance_windm_edsclrm.
  bool l_predict_upwp_vpwp{};

  // Flag to implement the nontraditional Coriolis terms in the prognostic equations of <w'w'>,
  // <u'w'>, and <u'u'>.
  bool l_ho_nontrad_coriolis{};

  // Flag to implement the traditional Coriolis terms in the prognostic equations of <v'w'> and
  // <u'w'>.
  bool l_ho_trad_coriolis{};

  // Flag to base the threshold minimum value of wp2 on keeping the overall correlation of w and x
  // (w and rt, as well as w and theta-l) within the limits of -max_mag_correlation_flux to
  // max_mag_correlation_flux.
  bool l_min_wp2_from_corr_wx{};

  // Flag to base the threshold minimum value of xp2 (rtp2 and thlp2) on keeping the overall
  // correlation of w and x within the limits of -max_mag_correlation_flux to
  // max_mag_correlation_flux.
  bool l_min_xp2_from_corr_wx{};

  // Flag to use cloud fraction to adjust the value of the turbulent dissipation coefficient, C2.
  bool l_C2_cloud_frac{};
  bool l_diffuse_rtm_and_thlm{};       // Diffuses rtm and thlm
  bool l_stability_correct_Kh_N2_zm{}; // Divides Kh_N2_zm by a stability factor
  bool l_calc_thlp2_rad{};             // Include the contribution of radiation to thlp2

  // This flag determines whether we want to use an upwind differencing approximation rather than a
  // centered differencing for turbulent or mean advection terms. It affects rtp2, thlp2, up2, vp2,
  // sclrp2, rtpthlp, sclrprtp, & sclrpthlp.
  bool l_upwind_xpyp_ta{};

  // This flag determines whether we want to use an upwind differencing approximation rather than a
  // centered differencing for turbulent or mean advection terms. It affects rtm, thlm, sclrm, um
  // and vm.
  bool l_upwind_xm_ma{};
  bool l_uv_nudge{};  // For wind speed nudging.
  bool l_rtm_nudge{}; // For rtm nudging

  // Use 2 calls to pdf_closure and the trapezoidal rule to compute the varibles that are output
  // from high order closure
  bool l_vert_avg_closure{};

  // If true, the trapezoidal rule is called for the thermodynamic-level variables output from
  // pdf_closure.
  bool l_trapezoidal_rule_zt{};

  // If true, the trapezoidal rule is called for three momentum-level variables - wpthvp, thlpthvp,
  // and rtpthvp - output from pdf_closure.
  bool l_trapezoidal_rule_zm{};

  // This logical flag determines whether or not to call subroutine pdf_closure twice.  If true,
  // pdf_closure is called first on thermodynamic levels and then on momentum levels so that each
  // variable is computed on its native level.  If false, pdf_closure is only called on
  // thermodynamic levels, and variables which belong on momentum levels are interpolated.
  bool l_call_pdf_closure_twice{};

  // Use the standard discretization for the turbulent advection terms.  Setting to .false. means
  // that a_1 and a_3 are pulled outside of the derivative in advance_wp2_wp3_module.F90 and in
  // advance_xp2_xpyp_module.F90.
  bool l_standard_term_ta{};

  // Flag to use an "upwind" discretization rather than a centered discretization for the portion of
  // the wp3 turbulent advection term for ADG1 that is linearized in terms of wp3<t+1>. (Requires
  // ADG1 PDF and l_standard_term_ta).
  bool l_partial_upwind_wp3{};

  // This flag determines whether we want to use an upwind differencing approximation rather than a
  // centered differencing for turbulent advection terms. It affects  wpxp only.
  bool l_godunov_upwind_wpxp_ta{};

  // This flag determines whether we want to use an upwind differencing approximation rather than a
  // centered differencing for turbulent advection terms. It affects xpyp only.
  bool l_godunov_upwind_xpyp_ta{};

  // Use cloud_cover and rcm_in_layer to help boost cloud_frac and rcm to help increase cloudiness
  // at coarser grid resolutions.
  bool l_use_cloud_cover{};
  bool l_diagnose_correlations{};      // Diagnose correlations instead of using fixed ones
  bool l_calc_w_corr{};                // Calculate the correlations between w and the hydrometeors
  bool l_const_Nc_in_cloud{};          // Use a constant cloud droplet conc. within cloud (K&K)
  bool l_fix_w_chi_eta_correlations{}; // Use a fixed correlation for s and t Mellor(chi/eta)

  // Use tau_N2_zm instead of tau_zm in wpxp_pr1 stability correction
  bool l_stability_correct_tau_zm{};

  // In wp2 equation, use a dissipation formula of -(2/3)*em/tau_zm, as in Bougeault (1981)
  bool l_damp_wp2_using_em{};
  bool l_do_expldiff_rtm_thlm{};  // Diffuse rtm and thlm explicitly
  bool l_Lscale_plume_centered{}; // Alternate that uses the PDF to compute the perturbed values
  bool l_use_C7_Richardson{};     // Parameterize C7 based on Richardson number
  bool l_use_C11_Richardson{};    // Parameterize C11 and C16 based on Richardson number
  bool l_use_shear_Richardson{};  // Use shear in the calculation of Richardson number

  // Use a different formula for the Brunt-Vaisala frequency in saturated atmospheres (from Durran
  // and Klemp, 1982)
  bool l_brunt_vaisala_freq_moist{};
  bool l_use_thvm_in_bv_freq{};   // Use thvm in the calculation of Brunt-Vaisala frequency
  bool l_rcm_supersat_adj{};      // Add excess supersaturated vapor to cloud water
  bool l_damp_wp3_Skw_squared{};  // Set damping on wp3 to use Skw^2 rather than Skw^4
  bool l_prescribed_avg_deltaz{}; // If .true., avg_deltaz = deltaz
  bool l_lmm_stepping{};          // Apply Linear Multistep Method (LMM) Stepping
  bool l_e3sm_config{};           // Run model with E3SM settings

  // Flag used to calculate convective velocity using a variable estimate of layer depth based on
  // the depth over which wpthlp is positive near the ground when true More information can be found
  // by Looking at issue #905 on the clubb repo
  bool l_vary_convect_depth{};
  bool l_use_tke_in_wp3_pr_turb_term{};     // Use TKE formulation for wp3 pr_turb term
  bool l_use_tke_in_wp2_wp3_K_dfsn{};       // Use TKE in eddy diffusion for wp2 and wp3
  bool l_use_wp3_lim_with_smth_Heaviside{}; // Flag to activate mods on wp3 limiters for conv test

  // Use smoothed Heaviside 'Preskin' function in the calculation of H_invrs_tau_wpxp_N2 in
  // src/CLUBB_core/mixing_length.F90
  bool l_smooth_Heaviside_tau_wpxp{};
  bool l_modify_limiters_for_cnvg_test{}; // Flag to activate mods on limiters for conv test
  bool l_enable_relaxed_clipping{}; // Flag to relax clipping on wpxp in xm_wpxp_clipping_and_stats
  bool l_linearize_pbl_winds{};     // Code to linearize PBL winds
  bool l_mono_flux_lim_thlm{};      // Flag to turn on monotonic flux limiter for thlm
  bool l_mono_flux_lim_rtm{};       // Flag to turn on monotonic flux limiter for rtm
  bool l_mono_flux_lim_um{};        // Flag to turn on monotonic flux limiter for um
  bool l_mono_flux_lim_vm{};        // Flag to turn on monotonic flux limiter for vm

  // Flag to implement monotonic flux limiter code that eliminates spurious drying tendencies at
  // model top
  bool l_mono_flux_lim_spikefix{};

  // Use to determine whether a host model has already applied the surface flux, to avoid double
  // counting.
  bool l_host_applies_sfc_fluxes{};

  // Turn on additional hole-filling for wp2 that takes TKE from up2 and vp2, if necessary
  bool l_wp2_fill_holes_tke{};
  bool l_add_dycore_grid{}; // Turn on remapping values from a dycore grid
  int debug_level{};        // Amount of debugging information
  int iostat{};             // iostat for optional namelist reads
  int iisclr_rt{-1};        // Native scalar IDs copied into shared metadata
  int iisclr_thl{-1};
  int iisclr_CO2{-1};
  int iiedsclr_rt{-1};
  int iiedsclr_thl{-1};
  int iiedsclr_CO2{-1};
  int i{};
  int j{};
  DriverArray1 sclr_tol_nl = workspace_scope.view<DriverArray1>(
      "sclr_tol_nl", sclr_max); // Thresholds on the passive scalars     [units vary]
  real_t sfc_elevation_nl{};
  real_t p_sfc_nl{};
  real_t T_sfc_nl{};
  real_t fcor_nl{};
  DriverArray2 momentum_heights{};
  DriverArray2 thermodynamic_heights{};
  DriverArray2 momentum_heights_flip{};
  DriverArray2 thermodynamic_heights_flip{};

  // ----------------------------------- Begin Code -----------------------------------
  timer_initialize();
  timer_start("init_clubb_case");

  //  Pick some default values for model_setting.  Some variables are initialized
  //  at declaration in module clubb_model_settings.
  T_sfc_nl                   = real_t(288.);
  p_sfc_nl                   = real_t(1000.e2);
  sens_ht                    = real_t(0.);
  latent_ht                  = real_t(0.);
  fcor_nl                    = real_t(1.e-4);
  T0                         = real_t(300.);
  ts_nudge                   = real_t(86400.);
  rtm_min                    = std::numeric_limits<real_t>::epsilon();
  rtm_nudge_max_altitude     = real_t(10000.);
  nzmax                      = 75;
  grid_type                  = 1;
  deltaz_nl                  = real_t(40.);
  zm_init_nl                 = real_t(0.);
  zm_top_nl                  = real_t(3500.);
  zt_grid_fname              = "";
  zm_grid_fname              = "";
  day                        = 22;
  month                      = 6;
  year                       = 1969;
  lat_vals                   = real_t(15.);
  lon_vals                   = (-real_t(56.5));
  sfctype                    = 0;
  time_initial               = real_t(0.);
  time_final                 = real_t(21600.);
  time_current               = real_t(0.);
  dt_main                    = real_t(60.);
  dt_rad                     = real_t(600.);
  forcings_file_path         = "";
  l_t_dependent              = false;
  l_input_xpwp_sfc           = false;
  l_ignore_forcings          = false;
  l_modify_ic_with_cubic_int = false;
  l_modify_bc_for_cnvg_test  = false;
  l_sfc_already_initialized  = false;

  // Sponge settings are namelist-only native types in the supported port scope.
  thlm_sponge_damp_settings    = {real_t(60), real_t(1800), real_t(0.25), false};
  rtm_sponge_damp_settings     = {real_t(60), real_t(1800), real_t(0.25), false};
  uv_sponge_damp_settings      = {real_t(60), real_t(1800), real_t(0.25), false};
  wp2_sponge_damp_settings     = {real_t(60), real_t(1800), real_t(0.25), false};
  wp3_sponge_damp_settings     = {real_t(60), real_t(1800), real_t(0.25), false};
  up2_vp2_sponge_damp_settings = {real_t(60), real_t(1800), real_t(0.25), false};
  l_restart                    = false;
  l_input_fields               = false;
  l_readiopdata                = false;
  restart_path_case            = "none";
  time_restart                 = real_t(0.);
  debug_level                  = 2;

  //  Use the Flatau polynomial approximation for computing saturation in advance_clubb_core_module
  saturation_formula    = saturation_flatau;
  sclr_dim              = 0;
  sclr_idx.iisclr_thl   = -2;
  sclr_idx.iisclr_rt    = -2;
  sclr_idx.iisclr_CO2   = -2;
  edsclr_dim            = 0;
  sclr_idx.iiedsclr_thl = -2;
  sclr_idx.iiedsclr_rt  = -2;
  sclr_idx.iiedsclr_CO2 = -2;
  l_first_write_to_file = true;

  // Pick some default values for stats_setting; other variables are set in
  // module stats_variables
  fname_prefix          = "";
  stats_output_filename = default_stats_output_filename;
  stats_tstart          = real_t(0.);
  stats_tend            = real_t(21600.);
  output_dir            = "../output/";
  ngrdcol               = 1;
  total_param_sets      = 1;
  batch_size            = (-1);
  iunit                 = 10;

  // Resolve total parameter count, runtime batch size, and the CLUBB parameter
  // matrices in one place:
  // - If clubb_params_in is provided, it defines both the total and runtime
  // column count, so batching is inactive.
  // - Otherwise, read total ngrdcol and optional batch_size from
  // &multicol_def, read the full parameter matrix, and select the first
  // runtime batch for initialization.
  if (clubb_params_in) {

    if ((*clubb_params_in).extent_int(0) < 1) {
      init_default_err_info_api(workspace, 1, err_info);
      std::cerr << "clubb_params_in must have at least one column" << "\n";

      //  General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      timer_stop("init_clubb_case");
      return;
    }

    if ((*clubb_params_in).extent_int(1) != nparams) {
      init_default_err_info_api(workspace, std::max({(*clubb_params_in).extent_int(0), 1}),
                                err_info);
      std::cerr << "clubb_params_in has invalid shape: size(2) must equal nparams"
                << "\n";

      //  General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      timer_stop("init_clubb_case");
      return;
    }

    total_param_sets = (*clubb_params_in).extent_int(0);
    batch_size       = total_param_sets;
    ngrdcol          = batch_size;

    if (clubb_params_all.is_allocated()) {
      clubb_params_all = {};
    }

    clubb_params_all = DriverArray2("clubb_params_all", total_param_sets, nparams);

    for (int iparam = 0; iparam < clubb_params_all.extent_int(1); ++iparam) {

      for (int icol = 0; icol < clubb_params_all.extent_int(0); ++icol) {
        clubb_params_all(icol, iparam) = (*clubb_params_in)(icol, iparam);
      }
    }

    if (clubb_params.is_allocated()) {
      clubb_params = {};
    }

    clubb_params = DriverArray2("clubb_params", ngrdcol, nparams);

    for (int iparam = 0; iparam < clubb_params.extent_int(1); ++iparam) {

      for (int icol = 0; icol < clubb_params.extent_int(0); ++icol) {
        clubb_params(icol, iparam) = clubb_params_all(icol, iparam);
      }
    }
  } else {
    //  Optional backwards-compat read: if &multicol_def is absent, default to single column.
    read_multicol_def(runfile, ngrdcol, batch_size, iostat);

    if (iostat != 0) {
      total_param_sets = 1;
      batch_size       = 1;
    } else {

      if (batch_size == -1) {
        batch_size = ngrdcol;
      }

      total_param_sets = ngrdcol;
    }

    if (total_param_sets < 1) {
      init_default_err_info_api(workspace, 1, err_info);
      std::cerr << "ngrdcol in &multicol_def must be >= 1" << "\n";

      //  General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      timer_stop("init_clubb_case");
      return;
    }

    if (batch_size < 1) {
      init_default_err_info_api(workspace, 1, err_info);
      std::cerr << "batch_size in &multicol_def must be >= 1" << "\n";

      //  General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      timer_stop("init_clubb_case");
      return;
    }

    if (batch_size > total_param_sets) {
      init_default_err_info_api(workspace, 1, err_info);
      std::cerr << "batch_size in &multicol_def must be <= ngrdcol" << "\n";

      //  General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      timer_stop("init_clubb_case");
      return;
    }

    if (total_param_sets % batch_size != 0) {
      init_default_err_info_api(workspace, 1, err_info);
      std::cerr << "ngrdcol in &multicol_def must be evenly divisible by batch_size"
                << "\n";

      //  General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      timer_stop("init_clubb_case");
      return;
    }

    if (clubb_params_all.is_allocated()) {
      clubb_params_all = {};
    }

    clubb_params_all = DriverArray2("clubb_params_all", total_param_sets, nparams);

    // This call reads the full parameter matrix from &clubb_params_nl.
    init_clubb_params_api(total_param_sets, iunit, runfile, clubb_params_all);
    ngrdcol = batch_size;

    if (clubb_params.is_allocated()) {
      clubb_params = {};
    }

    clubb_params = DriverArray2("clubb_params", ngrdcol, nparams);

    for (int iparam = 0; iparam < clubb_params.extent_int(1); ++iparam) {

      for (int icol = 0; icol < clubb_params.extent_int(0); ++icol) {
        clubb_params(icol, iparam) = clubb_params_all(icol, iparam);
      }
    }
  }

  // Initialize err_info with default values for one column or more,
  // depending on the resolved ngrdcol.
  init_default_err_info_api(workspace, ngrdcol, err_info);

  // Fill err_info with dummy values
  // The reshapes create arrays of dimension(ngrdcol) filled with values
  // lat_vals and lon_vals, respectively
  // Populating err_info with the available info
  // We don't want to populate the parallelization info
  // since that is only available outside of run_clubb
  DriverArray1 host_lat = workspace_scope.view<DriverArray1>("host_lat", ngrdcol),
               host_lon = workspace_scope.view<DriverArray1>("host_lon", ngrdcol);
  Kokkos::deep_copy(host_lat, lat_vals);
  Kokkos::deep_copy(host_lon, lon_vals);
  set_err_info_values_api(workspace, ngrdcol,         // In
                          err_info,                   // InOut
                          std::nullopt, std::nullopt, // Optional In
                          host_lat, host_lon);        // Optional In

  // Default values for the soil scheme
  initialize_soil_veg(1,                                                            // In
                      deep_soil_T_in_K_init, sfc_soil_T_in_K_init, veg_T_in_K_init, // Out
                      soil_vegetation);                                             // Out

  set_default_clubb_config_flags_api(iiPDF_type,                        // Out
                                     ipdf_call_placement,               // Out
                                     penta_solve_method,                // Out
                                     tridiag_solve_method,              // Out
                                     saturation_formula,                // Out
                                     Lscale_method,                     // Out
                                     grid_remap_method,                 // Out
                                     grid_adapt_in_time_method,         // Out
                                     fill_holes_type,                   // Out
                                     l_use_precip_frac,                 // Out
                                     l_predict_upwp_vpwp,               // Out
                                     l_ho_nontrad_coriolis,             // Out
                                     l_ho_trad_coriolis,                // Out
                                     l_min_wp2_from_corr_wx,            // Out
                                     l_min_xp2_from_corr_wx,            // Out
                                     l_C2_cloud_frac,                   // Out
                                     l_diffuse_rtm_and_thlm,            // Out
                                     l_stability_correct_Kh_N2_zm,      // Out
                                     l_calc_thlp2_rad,                  // Out
                                     l_upwind_xpyp_ta,                  // Out
                                     l_upwind_xm_ma,                    // Out
                                     l_uv_nudge,                        // Out
                                     l_rtm_nudge,                       // Out
                                     l_vert_avg_closure,                // Out
                                     l_trapezoidal_rule_zt,             // Out
                                     l_trapezoidal_rule_zm,             // Out
                                     l_call_pdf_closure_twice,          // Out
                                     l_standard_term_ta,                // Out
                                     l_partial_upwind_wp3,              // Out
                                     l_godunov_upwind_wpxp_ta,          // Out
                                     l_godunov_upwind_xpyp_ta,          // Out
                                     l_use_cloud_cover,                 // Out
                                     l_diagnose_correlations,           // Out
                                     l_calc_w_corr,                     // Out
                                     l_const_Nc_in_cloud,               // Out
                                     l_fix_w_chi_eta_correlations,      // Out
                                     l_stability_correct_tau_zm,        // Out
                                     l_damp_wp2_using_em,               // Out
                                     l_do_expldiff_rtm_thlm,            // Out
                                     l_Lscale_plume_centered,           // Out
                                     l_use_C7_Richardson,               // Out
                                     l_use_C11_Richardson,              // Out
                                     l_use_shear_Richardson,            // Out
                                     l_brunt_vaisala_freq_moist,        // Out
                                     l_use_thvm_in_bv_freq,             // Out
                                     l_rcm_supersat_adj,                // Out
                                     l_damp_wp3_Skw_squared,            // Out
                                     l_prescribed_avg_deltaz,           // Out
                                     l_lmm_stepping,                    // Out
                                     l_e3sm_config,                     // Out
                                     l_vary_convect_depth,              // Out
                                     l_use_tke_in_wp3_pr_turb_term,     // Out
                                     l_use_tke_in_wp2_wp3_K_dfsn,       // Out
                                     l_use_wp3_lim_with_smth_Heaviside, // Out
                                     l_smooth_Heaviside_tau_wpxp,       // Out
                                     l_modify_limiters_for_cnvg_test,   // Out
                                     l_enable_relaxed_clipping,         // Out
                                     l_linearize_pbl_winds,             // Out
                                     l_mono_flux_lim_thlm,              // Out
                                     l_mono_flux_lim_rtm,               // Out
                                     l_mono_flux_lim_um,                // Out
                                     l_mono_flux_lim_vm,                // Out
                                     l_mono_flux_lim_spikefix,          // Out
                                     l_host_applies_sfc_fluxes,         // Out
                                     l_wp2_fill_holes_tke,              // Out
                                     l_add_dycore_grid);                // Out

  // Read namelist file
  timer_start("read_namelists");

  read_model_setting(runfile,                                                // In
                     runtype,                                                // InOut
                     nzmax,                                                  // InOut
                     grid_type,                                              // InOut
                     deltaz_nl,                                              // InOut
                     zm_init_nl,                                             // InOut
                     zm_top_nl,                                              // InOut
                     zt_grid_fname,                                          // InOut
                     zm_grid_fname,                                          // InOut
                     day,                                                    // InOut
                     month,                                                  // InOut
                     year,                                                   // InOut
                     lat_vals,                                               // InOut
                     lon_vals,                                               // InOut
                     sfc_elevation_nl,                                       // InOut
                     time_initial,                                           // InOut
                     time_final,                                             // InOut
                     dt_main,                                                // InOut
                     dt_rad,                                                 // InOut
                     sfctype,                                                // InOut
                     T_sfc_nl,                                               // InOut
                     p_sfc_nl,                                               // InOut
                     sens_ht,                                                // InOut
                     latent_ht,                                              // InOut
                     fcor_nl,                                                // InOut
                     T0,                                                     // InOut
                     ts_nudge,                                               // InOut
                     forcings_file_path,                                     // InOut
                     l_t_dependent,                                          // InOut
                     l_input_xpwp_sfc,                                       // InOut
                     l_ignore_forcings,                                      // InOut
                     l_modify_ic_with_cubic_int,                             // InOut
                     l_modify_bc_for_cnvg_test,                              // InOut
                     l_soil_veg,                                             // InOut
                     l_uv_nudge,                                             // InOut
                     l_restart,                                              // InOut
                     restart_path_case,                                      // InOut
                     time_restart,                                           // InOut
                     l_input_fields,                                         // InOut
                     l_readiopdata,                                          // InOut
                     debug_level,                                            // InOut
                     sclr_tol_nl,                                            // InOut
                     sclr_dim,                                               // InOut
                     iisclr_thl,                                             // InOut
                     iisclr_rt,                                              // InOut
                     iisclr_CO2,                                             // InOut
                     edsclr_dim,                                             // InOut
                     iiedsclr_thl,                                           // InOut
                     iiedsclr_rt,                                            // InOut
                     iiedsclr_CO2,                                           // InOut
                     l_rtm_nudge,                                            // InOut
                     rtm_min,                                                // InOut
                     rtm_nudge_max_altitude,                                 // InOut
                     l_diagnose_correlations,                                // InOut
                     l_calc_w_corr,                                          // InOut
                     thlm_sponge_damp_settings, rtm_sponge_damp_settings,    // InOut
                     uv_sponge_damp_settings, wp2_sponge_damp_settings,      // InOut
                     wp3_sponge_damp_settings, up2_vp2_sponge_damp_settings, // InOut
                     iostat);                                                // Out

  if (iostat != 0)
    throw std::runtime_error("Unable to read &model_setting namelist");

  // The JAX/Kokkos scope excludes these native branches; reject them before allocation.
  if (thlm_sponge_damp_settings.l_sponge_damping || rtm_sponge_damp_settings.l_sponge_damping ||
      uv_sponge_damp_settings.l_sponge_damping || wp2_sponge_damp_settings.l_sponge_damping ||
      wp3_sponge_damp_settings.l_sponge_damping || up2_vp2_sponge_damp_settings.l_sponge_damping ||
      l_input_fields || l_test_grid_generalization)
    throw std::runtime_error("Requested model option is outside the JAX/Kokkos port scope");
  stats_tstart = time_initial;
  stats_tend   = time_final;

  read_stats_setting(runfile,                  // In
                     l_stats,                  // InOut
                     fname_prefix,             // InOut
                     stats_output_filename,    // InOut
                     stats_tsamp,              // InOut
                     stats_tout,               // InOut
                     stats_tstart,             // InOut
                     stats_tend,               // InOut
                     l_allow_small_stats_tout, // InOut
                     output_dir,               // InOut
                     iostat);                  // Out

  if (iostat != 0)
    throw std::runtime_error("Unable to read &stats_setting namelist");

  if (stats_output_filename == default_stats_output_filename) {
    stats_output_filename = (fname_prefix + "_stats.nc");
  }

  read_configurable_clubb_flags_nl(runfile,                           // In
                                   iiPDF_type,                        // InOut
                                   ipdf_call_placement,               // InOut
                                   penta_solve_method,                // InOut
                                   tridiag_solve_method,              // InOut
                                   saturation_formula,                // InOut
                                   Lscale_method,                     // InOut
                                   grid_remap_method,                 // InOut
                                   grid_adapt_in_time_method,         // InOut
                                   fill_holes_type,                   // InOut
                                   l_upwind_xpyp_ta,                  // InOut
                                   l_upwind_xm_ma,                    // InOut
                                   l_vert_avg_closure,                // InOut
                                   l_standard_term_ta,                // InOut
                                   l_partial_upwind_wp3,              // InOut
                                   l_godunov_upwind_wpxp_ta,          // InOut
                                   l_godunov_upwind_xpyp_ta,          // InOut
                                   l_use_cloud_cover,                 // InOut
                                   l_rcm_supersat_adj,                // InOut
                                   l_damp_wp3_Skw_squared,            // InOut
                                   l_min_wp2_from_corr_wx,            // InOut
                                   l_min_xp2_from_corr_wx,            // InOut
                                   l_C2_cloud_frac,                   // InOut
                                   l_predict_upwp_vpwp,               // InOut
                                   l_ho_nontrad_coriolis,             // InOut
                                   l_ho_trad_coriolis,                // InOut
                                   l_stability_correct_tau_zm,        // InOut
                                   l_damp_wp2_using_em,               // InOut
                                   l_use_C7_Richardson,               // InOut
                                   l_use_precip_frac,                 // InOut
                                   l_do_expldiff_rtm_thlm,            // InOut
                                   l_use_C11_Richardson,              // InOut
                                   l_use_shear_Richardson,            // InOut
                                   l_prescribed_avg_deltaz,           // InOut
                                   l_diffuse_rtm_and_thlm,            // InOut
                                   l_stability_correct_Kh_N2_zm,      // InOut
                                   l_trapezoidal_rule_zt,             // InOut
                                   l_trapezoidal_rule_zm,             // InOut
                                   l_call_pdf_closure_twice,          // InOut
                                   l_Lscale_plume_centered,           // InOut
                                   l_brunt_vaisala_freq_moist,        // InOut
                                   l_use_thvm_in_bv_freq,             // InOut
                                   l_lmm_stepping,                    // InOut
                                   l_e3sm_config,                     // InOut
                                   l_vary_convect_depth,              // InOut
                                   l_use_tke_in_wp3_pr_turb_term,     // InOut
                                   l_use_tke_in_wp2_wp3_K_dfsn,       // InOut
                                   l_use_wp3_lim_with_smth_Heaviside, // InOut
                                   l_smooth_Heaviside_tau_wpxp,       // InOut
                                   l_modify_limiters_for_cnvg_test,   // InOut
                                   l_enable_relaxed_clipping,         // InOut
                                   l_linearize_pbl_winds,             // InOut
                                   l_mono_flux_lim_thlm,              // InOut
                                   l_mono_flux_lim_rtm,               // InOut
                                   l_mono_flux_lim_um,                // InOut
                                   l_mono_flux_lim_vm,                // InOut
                                   l_mono_flux_lim_spikefix,          // InOut
                                   l_host_applies_sfc_fluxes,         // InOut
                                   l_wp2_fill_holes_tke,              // InOut
                                   l_add_dycore_grid,                 // InOut
                                   iostat);                           // Out

  if (iostat != 0)
    throw std::runtime_error("Unable to read &configurable_clubb_flags_nl namelist");
  timer_stop("read_namelists");

  // Ensure output directory has a trailing slash.
  if (static_cast<int>(output_dir.size()) <= 0) {
    output_dir = "../output/";
  }

  noutdir = static_cast<int>(output_dir.size());

  if (output_dir.substr((noutdir - 1), 1) != "/") {
    output_dir = (output_dir + "/");
  }

  sclr_idx.iisclr_thl   = (iisclr_thl - 1);
  sclr_idx.iisclr_rt    = (iisclr_rt - 1);
  sclr_idx.iisclr_CO2   = (iisclr_CO2 - 1);
  sclr_idx.iiedsclr_thl = (iiedsclr_thl - 1);
  sclr_idx.iiedsclr_rt  = (iiedsclr_rt - 1);
  sclr_idx.iiedsclr_CO2 = (iiedsclr_CO2 - 1);
  sclr_tol              = DriverArray1("sclr_tol", sclr_dim);

  for (int sclr = 0; sclr < sclr_dim; ++sclr) {
    sclr_tol(sclr) = sclr_tol_nl(sclr);
  }

  p_sfc            = DriverArray1("p_sfc", ngrdcol);
  T_sfc            = DriverArray1("T_sfc", ngrdcol);
  fcor             = DriverArray1("fcor", ngrdcol);
  fcor_y           = DriverArray1("fcor_y", ngrdcol);
  dummy_dx         = DriverArray1("dummy_dx", ngrdcol);
  dummy_dy         = DriverArray1("dummy_dy", ngrdcol);
  deep_soil_T_in_K = DriverArray1("deep_soil_T_in_K", ngrdcol);
  sfc_soil_T_in_K  = DriverArray1("sfc_soil_T_in_K", ngrdcol);
  veg_T_in_K       = DriverArray1("veg_T_in_K", ngrdcol);
  // TODO: Native standalone allocates these two inputs without initializing
  // them. Views start at zero; resolve native initialization before requiring
  // standalone BFB for l_linearize_pbl_winds. Equal-input core calls agree.
  upwp_sfc_pert = DriverArray1("upwp_sfc_pert", ngrdcol);
  vpwp_sfc_pert = DriverArray1("vpwp_sfc_pert", ngrdcol);
  sfc_elevation = DriverArray1("sfc_elevation", ngrdcol);
  zm_init       = DriverArray1("zm_init", ngrdcol);
  zm_top        = DriverArray1("zm_top", ngrdcol);
  deltaz        = DriverArray1("deltaz", ngrdcol);

  for (int icol = 0; icol < sfc_elevation.extent_int(0); ++icol) {
    sfc_elevation(icol) = sfc_elevation_nl;
  }

  for (int icol = 0; icol < deltaz.extent_int(0); ++icol) {
    deltaz(icol) = deltaz_nl;
  }

  for (int icol = 0; icol < zm_init.extent_int(0); ++icol) {
    zm_init(icol) = zm_init_nl;
  }

  for (int icol = 0; icol < zm_top.extent_int(0); ++icol) {
    zm_top(icol) = zm_top_nl;
  }

  for (int icol = 0; icol < p_sfc.extent_int(0); ++icol) {
    p_sfc(icol) = p_sfc_nl;
  }

  for (int icol = 0; icol < T_sfc.extent_int(0); ++icol) {
    T_sfc(icol) = T_sfc_nl;
  }

  for (int icol = 0; icol < fcor.extent_int(0); ++icol) {
    fcor(icol) = fcor_nl;
  }

  // Suppose the vector Coriolis term is written as 2*Omega_vector x v_vector, where
  // Omega_vector is the Earth's angular velocity vector (omega_planet in code), which
  // points along the North Pole, and v_vector is the full air velocity vector in the
  // rotating frame of reference. Now suppose we break Omega_vector into a local upward
  // component, Omega_z*z_hat, and a local northward component, Omega_y*y_hat.
  // Then we define fcor = 2*Omega_z and fcor_y = 2*Omega_y.
  // Then fcor_y = 2 * Omega_vector \dot y_hat = 2 * |Omega_vector| * cos( angle between north
  // pole and local northward direction )
  for (int icol = 0; icol < fcor_y.extent_int(0); ++icol) {
    fcor_y(icol) = ((two * omega_planet) * std::cos((lat_vals * radians_per_deg)));
  }

  // fcor could reasonably be defined as below, but currently it is input via namelist
  // Hing Ong, 25 November 2025
  // fcor   = two * omega_planet * sin ( lat_vals * radians_per_deg )
  if (l_vert_avg_closure) {
    l_trapezoidal_rule_zt    = true;
    l_trapezoidal_rule_zm    = true;
    l_call_pdf_closure_twice = true;
  }

  // The case prefix including output directory.
  output_file_prefix = (output_dir + fname_prefix);

  // Resolve stats output path from output_dir plus the namelist-provided filename.
  if (static_cast<int>(stats_output_filename.size()) > 0) {

    if (stats_output_filename.substr(0, 1) == "/") {
      stats_output_path = stats_output_filename;
    } else {
      stats_output_path = (output_dir + stats_output_filename);
    }
  } else {
    stats_output_path = "";
  }

  // The filename for case setup
  case_info_file = ((output_dir + fname_prefix) + "_setup.txt");

  // Sanity check on passive scalars
  // When adding new 'ii' scalar indices, add them to this list.
  if (std::max({iisclr_CO2, iisclr_rt, iisclr_thl}) > sclr_dim) {
    std::cerr << err_info.err_header_global << "\n";
    std::cerr << "Passive scalar index exceeds sclr_dim " << " "
              << "iisclr_CO2 = " << " " << iisclr_CO2 << " "
              << "iisclr_rt = " << " " << iisclr_rt << " "
              << "iisclr_thl = " << " " << iisclr_thl << " "
              << "sclr_dim = " << " " << sclr_dim << "\n";

    //  General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    timer_stop("init_clubb_case");
    return;
  } else if (std::max({iiedsclr_CO2, iiedsclr_rt, iiedsclr_thl}) > edsclr_dim) {
    std::cerr << err_info.err_header_global << "\n";
    std::cerr << "Passive scalar index exceeds edsclr_dim " << " "
              << "iiedsclr_CO2 = " << " " << iiedsclr_CO2 << " "
              << "iiedsclr_rt = " << " " << iiedsclr_rt << " "
              << "iiedsclr_thl = " << " " << iiedsclr_thl << " "
              << "edsclr_dim = " << " " << edsclr_dim << "\n";

    //  General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    timer_stop("init_clubb_case");
    return;
  }

  // Set debug level
  set_clubb_debug_level_api(debug_level);

  if (debug_level < 0) {
    timer_stop("init_clubb_case");
    timer_disable();
  }

  if (clubb_at_least_debug_level_api(1)) {
    timer_start("write_case_info");

    if (l_write_to_file) {
      case_log.open(case_info_file, std::ios::out);
    }

    // Print the date and time
    write_date(case_log);
    write_text("", case_log);
    write_text("--------------------------------------------------", case_log);
    write_text("Latest git log entry", case_log);
    write_text("--------------------------------------------------", case_log);
    write_text("", case_log);

    if (l_write_to_file) {

      if (case_log.is_open())
        case_log.close();
    }

    if (l_write_to_file) {
      // Native command-line redirection is expressed with fixed commands and a file stream.
      std::ofstream metadata(case_info_file, std::ios::app);
      metadata << "Branch name:\n";
      append_command_output("git rev-parse --abbrev-ref HEAD", metadata);
      metadata << "(If no branch name is shown HEAD may be detached)\n\n";
      append_command_output("git log -n 1", metadata);
    }

    if (l_write_to_file) {
      case_log.open(case_info_file, std::ios::app);
    }

    write_text("\nA detailed git diff can be found at the end of this file\n", case_log);
    write_text("--------------------------------------------------", case_log);
    write_text("Tunable Parameters:", case_log);
    write_text("--------------------------------------------------", case_log);

    // Print the list of parameters that are being used before the run.
    write_text("Parameter          Value", l_write_to_file, case_log);
    write_text("---------          -----", l_write_to_file, case_log);

    for (j = 0; j < nparams; ++j) {
      write_text((text_at(params_list, j) + " = "), clubb_params(0, j), case_log);
    }

    write_text("--------------------------------------------------", case_log);
    write_text("Preprocessing Directives:", case_log);
    write_text("--------------------------------------------------", case_log);
#ifdef COAMPS_MICRO
    write_text("-DCOAMPS_MICRO enabled", case_log);
#else
    write_text("-DCOAMPS_MICRO disabled", case_log);
#endif
#ifdef TUNER
    write_text("-DTUNER enabled", case_log);
#else
    write_text("-DTUNER disabled", case_log);
#endif
#ifdef nooverlap
    write_text("-Dnooverlap enabled", case_log);
#else
    write_text("-Dnooverlap disabled", case_log);
#endif
#ifdef radoffline
    write_text("-Dradoffline enabled", case_log);
#else
    write_text("-Dradoffline disabled", case_log);
#endif
    write_text("-DUSE_BUGSrad_ocast_random disabled", case_log);
    write_text("--------------------------------------------------", case_log);
    write_text("&model_setting", case_log);
    write_text("--------------------------------------------------", case_log);
    write_text(("runtype = " + runtype), case_log);
    write_text("nzmax = ", nzmax, case_log);
    write_text("grid_type = ", grid_type, case_log);
    write_text("deltaz = ", deltaz, case_log);
    write_text("zm_init = ", zm_init, case_log);
    write_text("zm_top = ", zm_top, case_log);
    write_text(("zt_grid_fname = " + zt_grid_fname), case_log);
    write_text(("zm_grid_fname = " + zm_grid_fname), case_log);
    write_text("day = ", day, case_log);
    write_text("month = ", month, case_log);
    write_text("year = ", year, case_log);
    write_text("lat_vals = ", lat_vals, case_log);
    write_text("lon_vals = ", lon_vals, case_log);
    write_text("sfc_elevation = ", sfc_elevation_nl, case_log);
    write_text("time_initial = ", static_cast<real_t>(time_initial), case_log);
    write_text("time_final = ", static_cast<real_t>(time_final), case_log);
    write_text("dt_main = ", dt_main, case_log);
    write_text("dt_rad = ", dt_rad, case_log);
    write_text("sfctype = ", sfctype, case_log);
    write_text("T_sfc = ", T_sfc, case_log);
    write_text("p_sfc = ", p_sfc, case_log);
    write_text("sens_ht = ", sens_ht, case_log);
    write_text("latent_ht = ", latent_ht, case_log);
    write_text("fcor = ", fcor, case_log);
    write_text("fcor_y = ", fcor_y, case_log);
    write_text("T0 = ", T0, case_log);
    write_text("ts_nudge = ", ts_nudge, case_log);
    write_text(("forcings_file_path = " + forcings_file_path), case_log);
    write_text("l_t_dependent = ", l_t_dependent, case_log);
    write_text("l_ignore_forcings = ", l_ignore_forcings, case_log);
    write_text("l_modify_ic_with_cubic_int = ", l_modify_ic_with_cubic_int, case_log);
    write_text("l_modify_bc_for_cnvg_test = ", l_modify_bc_for_cnvg_test, case_log);
    write_text("l_input_xpwp_sfc = ", l_input_xpwp_sfc, case_log);
    write_text("saturation_formula = ", saturation_formula, case_log);
    write_text("thlm_sponge_damp_settings%l_sponge_damping = ",
               thlm_sponge_damp_settings.l_sponge_damping, case_log);
    write_text("rtm_sponge_damp_settings%l_sponge_damping = ",
               rtm_sponge_damp_settings.l_sponge_damping, case_log);
    write_text("uv_sponge_damp_settings%l_sponge_damping = ",
               uv_sponge_damp_settings.l_sponge_damping, case_log);
    write_text("wp2_sponge_damp_settings%l_sponge_damping = ",
               wp2_sponge_damp_settings.l_sponge_damping, case_log);
    write_text("wp3_sponge_damp_settings%l_sponge_damping = ",
               wp3_sponge_damp_settings.l_sponge_damping, case_log);
    write_text("up2_vp2_sponge_damp_settings%l_sponge_damping = ",
               up2_vp2_sponge_damp_settings.l_sponge_damping, case_log);
    write_text("thlm_sponge_damp_settings%tau_sponge_damp_min = ",
               thlm_sponge_damp_settings.tau_sponge_damp_min, case_log);
    write_text("thlm_sponge_damp_settings%tau_sponge_damp_max = ",
               thlm_sponge_damp_settings.tau_sponge_damp_max, case_log);
    write_text("thlm_sponge_damp_settings%sponge_damp_depth = ",
               thlm_sponge_damp_settings.sponge_damp_depth, case_log);
    write_text("rtm_sponge_damp_settings%tau_sponge_damp_min = ",
               rtm_sponge_damp_settings.tau_sponge_damp_min, case_log);
    write_text("rtm_sponge_damp_settings%tau_sponge_damp_max = ",
               rtm_sponge_damp_settings.tau_sponge_damp_max, case_log);
    write_text("rtm_sponge_damp_settings%sponge_damp_depth = ",
               rtm_sponge_damp_settings.sponge_damp_depth, case_log);
    write_text("uv_sponge_damp_settings%tau_sponge_damp_min = ",
               uv_sponge_damp_settings.tau_sponge_damp_min, case_log);
    write_text("uv_sponge_damp_settings%tau_sponge_damp_max = ",
               uv_sponge_damp_settings.tau_sponge_damp_max, case_log);
    write_text("uv_sponge_damp_settings%sponge_damp_depth = ",
               uv_sponge_damp_settings.sponge_damp_depth, case_log);
    write_text("wp2_sponge_damp_settings%tau_sponge_damp_min = ",
               wp2_sponge_damp_settings.tau_sponge_damp_min, case_log);
    write_text("wp2_sponge_damp_settings%tau_sponge_damp_max = ",
               wp2_sponge_damp_settings.tau_sponge_damp_max, case_log);
    write_text("wp2_sponge_damp_settings%sponge_damp_depth = ",
               wp2_sponge_damp_settings.sponge_damp_depth, case_log);
    write_text("wp3_sponge_damp_settings%tau_sponge_damp_min = ",
               wp3_sponge_damp_settings.tau_sponge_damp_min, case_log);
    write_text("wp3_sponge_damp_settings%tau_sponge_damp_max = ",
               wp3_sponge_damp_settings.tau_sponge_damp_max, case_log);
    write_text("wp3_sponge_damp_settings%sponge_damp_depth = ",
               wp3_sponge_damp_settings.sponge_damp_depth, case_log);
    write_text("up2_vp2_sponge_damp_settings%tau_sponge_damp_min = ",
               up2_vp2_sponge_damp_settings.tau_sponge_damp_min, case_log);
    write_text("up2_vp2_sponge_damp_settings%tau_sponge_damp_max = ",
               up2_vp2_sponge_damp_settings.tau_sponge_damp_max, case_log);
    write_text("up2_vp2_sponge_damp_settings%sponge_damp_depth = ",
               up2_vp2_sponge_damp_settings.sponge_damp_depth, case_log);
    write_text("l_soil_veg = ", l_soil_veg, case_log);
    write_text("l_restart = ", l_restart, case_log);
    write_text("l_input_fields = ", l_input_fields, case_log);
    write_text(("restart_path_case = " + restart_path_case), case_log);
    write_text("time_restart = ", static_cast<real_t>(time_restart), case_log);
    write_text("debug_level = ", debug_level, case_log);
    write_text("sclr_dim = ", sclr_dim, case_log);
    write_text("edsclr_dim = ", edsclr_dim, case_log);
    write_text("iisclr_thl = ", iisclr_thl, case_log);
    write_text("iisclr_rt = ", iisclr_rt, case_log);
    write_text("iisclr_CO2 = ", iisclr_CO2, case_log);
    write_text("sclr_tol = ", Kokkos::subview(sclr_tol, std::make_pair(0, sclr_dim)), case_log);
    write_text("--------------------------------------------------", case_log);
    write_text("&stats_setting", case_log);
    write_text("--------------------------------------------------", case_log);
    write_text("l_stats = ", l_stats, case_log);
    write_text(("output_dir = " + output_dir), case_log);
    write_text(("fname_prefix = " + fname_prefix), case_log);
    write_text(("stats_output_filename = " + stats_output_filename), case_log);
    write_text("stats_tsamp = ", static_cast<real_t>(stats_tsamp), case_log);
    write_text("stats_tout = ", static_cast<real_t>(stats_tout), case_log);
    write_text("stats_tstart = ", static_cast<real_t>(stats_tstart), case_log);
    write_text("stats_tend = ", static_cast<real_t>(stats_tend), case_log);
    write_text("l_allow_small_stats_tout = ", l_allow_small_stats_tout, case_log);
    write_text("Constant flags:", case_log);
    write_text("l_pos_def = ", l_pos_def, case_log);
    write_text("l_gamma_Skw = ", l_gamma_Skw, case_log);
    write_text("Constant tolerances [units]", case_log);
    write_text("rt_tol [kg/kg] = ", rt_tol, case_log);
    write_text("thl_tol [K] = ", thl_tol, case_log);
    write_text("w_tol [m/s] = ", w_tol, case_log);

    if (l_write_to_file) {

      if (case_log.is_open())
        case_log.close();
    }

    iunit_grid_adaptation = ((iunit + 10));
    fname_grid_adaptation = ((output_dir + runtype) + "_grid_adapt.txt");
    timer_stop("write_case_info");
  }

  //  Allocate stretched grid altitude arrays.
  momentum_heights = workspace_scope.view<DriverArray2>("momentum_heights", ngrdcol, nzmax);
  thermodynamic_heights =
      workspace_scope.view<DriverArray2>("thermodynamic_heights", ngrdcol, (nzmax - 1));

  for (i = 0; i < ngrdcol; ++i) {

    read_grid_heights(workspace, nzmax,                                       // In
                      grid_type,                                              // InOut
                      zm_grid_fname, zt_grid_fname,                           // InOut
                      iunit,                                                  // InOut
                      Kokkos::subview(momentum_heights, i, Kokkos::ALL),      // InOut
                      Kokkos::subview(thermodynamic_heights, i, Kokkos::ALL), // InOut
                      err_info);                                              // Out

    if (host_element(err_info.err_code, i) == clubb_fatal_error) {
      std::cerr << text_at(err_info.err_header, i) << "\n";
    }
  }

  if (any_equal(err_info.err_code, clubb_fatal_error)) {
    std::cerr << "Error in read_grid_heights" << "\n";
    throw std::runtime_error("Fatal error");
  }

  // These numbers represent host model horizontal grid spacing
  // which for a single column simulation is effectively infinite
  for (int icol = 0; icol < dummy_dx.extent_int(0); ++icol) {
    dummy_dx(icol) = real_t(1.0e6);
  }

  for (int icol = 0; icol < dummy_dy.extent_int(0); ++icol) {
    dummy_dy(icol) = real_t(1.0e6);
  }

  // Setup microphysical fields
  init_microphys(workspace, iunit, runtype, runfile, case_info_file, // In
                 dummy_dx(0), dummy_dy(0),                           // In
                 Kokkos::subview(clubb_params, 0, Kokkos::ALL),      // In
                 l_diagnose_correlations,                            // In
                 l_const_Nc_in_cloud,                                // InOut
                 l_fix_w_chi_eta_correlations,                       // InOut
                 hydromet_dim, pdf_dim, hm_metadata,                 // Out
                 silhs_config_flags,                                 // Out
                 vert_decorr_coef,                                   // Out
                 corr_array_n_cloud, corr_array_n_below,             // Out
                 microphys, kk, morrison, silhs);                    // Case state

  if (microphys_scheme != "none" && microphys_scheme != "khairoutdinov_kogan" &&
      microphys_scheme != "morrison")
    throw std::runtime_error("Unsupported microphysics scheme in the Kokkos driver");

  if (microphys.l_gfdl_activation)
    throw std::runtime_error("GFDL activation is outside the port scope");

  // Setup radiation parameters
  init_radiation(iunit, runfile, case_info_file, // In
                 l_calc_thlp2_rad,               // InOut
                 parameters_radiation);          // InOut

  if (rad_scheme != "none" && rad_scheme != "simplified" && rad_scheme != "simplified_bomex" &&
      rad_scheme != "lba")
    throw std::runtime_error("Unsupported radiation scheme in the Kokkos driver");

  if (l_silhs_rad)
    throw std::runtime_error("SILHS radiation is outside the port scope");

  if (rad_scheme == "lba") {
    simple_rad_lba_init(iunit, forcings_file_path, // In
                        simple_radiation);         // Out
  }

  // Initialize CLUBB configurable flags type
  initialize_clubb_config_flags_type_api(iiPDF_type,                        // In
                                         ipdf_call_placement,               // In
                                         penta_solve_method,                // In
                                         tridiag_solve_method,              // In
                                         saturation_formula,                // In
                                         Lscale_method,                     // In
                                         grid_remap_method,                 // In
                                         grid_adapt_in_time_method,         // In
                                         fill_holes_type,                   // In
                                         l_use_precip_frac,                 // In
                                         l_predict_upwp_vpwp,               // In
                                         l_ho_nontrad_coriolis,             // In
                                         l_ho_trad_coriolis,                // In
                                         l_min_wp2_from_corr_wx,            // In
                                         l_min_xp2_from_corr_wx,            // In
                                         l_C2_cloud_frac,                   // In
                                         l_diffuse_rtm_and_thlm,            // In
                                         l_stability_correct_Kh_N2_zm,      // In
                                         l_calc_thlp2_rad,                  // In
                                         l_upwind_xpyp_ta,                  // In
                                         l_upwind_xm_ma,                    // In
                                         l_uv_nudge,                        // In
                                         l_rtm_nudge,                       // In
                                         l_vert_avg_closure,                // In
                                         l_trapezoidal_rule_zt,             // In
                                         l_trapezoidal_rule_zm,             // In
                                         l_call_pdf_closure_twice,          // In
                                         l_standard_term_ta,                // In
                                         l_partial_upwind_wp3,              // In
                                         l_godunov_upwind_wpxp_ta,          // In
                                         l_godunov_upwind_xpyp_ta,          // In
                                         l_use_cloud_cover,                 // In
                                         l_diagnose_correlations,           // In
                                         l_calc_w_corr,                     // In
                                         l_const_Nc_in_cloud,               // In
                                         l_fix_w_chi_eta_correlations,      // In
                                         l_stability_correct_tau_zm,        // In
                                         l_damp_wp2_using_em,               // In
                                         l_do_expldiff_rtm_thlm,            // In
                                         l_Lscale_plume_centered,           // In
                                         l_use_C7_Richardson,               // In
                                         l_use_C11_Richardson,              // In
                                         l_use_shear_Richardson,            // In
                                         l_brunt_vaisala_freq_moist,        // In
                                         l_use_thvm_in_bv_freq,             // In
                                         l_rcm_supersat_adj,                // In
                                         l_damp_wp3_Skw_squared,            // In
                                         l_prescribed_avg_deltaz,           // In
                                         l_lmm_stepping,                    // In
                                         l_e3sm_config,                     // In
                                         l_vary_convect_depth,              // In
                                         l_use_tke_in_wp3_pr_turb_term,     // In
                                         l_use_tke_in_wp2_wp3_K_dfsn,       // In
                                         l_use_wp3_lim_with_smth_Heaviside, // In
                                         l_smooth_Heaviside_tau_wpxp,       // In
                                         l_modify_limiters_for_cnvg_test,   // In
                                         l_enable_relaxed_clipping,         // In
                                         l_linearize_pbl_winds,             // In
                                         l_mono_flux_lim_thlm,              // In
                                         l_mono_flux_lim_rtm,               // In
                                         l_mono_flux_lim_um,                // In
                                         l_mono_flux_lim_vm,                // In
                                         l_mono_flux_lim_spikefix,          // In
                                         l_host_applies_sfc_fluxes,         // In
                                         l_wp2_fill_holes_tke,              // In
                                         l_add_dycore_grid,                 // In
                                         clubb_config_flags);               // Out

  // Adaptive-grid bodies are omitted under the same feature scope as JAX.
  if (clubb_config_flags.grid_adapt_in_time_method > no_grid_adaptation)
    throw std::runtime_error("Adaptive grids are outside the port scope");

  if (clubb_at_least_debug_level_api(1)) {
    timer_start("write_case_info");

    if (l_write_to_file) {
      case_log.open(case_info_file, std::ios::app);
    }

    write_text("--------------------------------------------------", case_log);
    write_text("&configurable_clubb_flags_nl", case_log);
    write_text("--------------------------------------------------", case_log);

    // Native prints the same flag rows to stdout and the optional case log.
    print_clubb_config_flags_api(std::cout, clubb_config_flags);
    if (case_log.is_open())
      print_clubb_config_flags_api(case_log, clubb_config_flags);
    write_text("--------------------------------------------------", case_log);
    write_text("git diff src/", case_log);
    write_text("--------------------------------------------------", case_log);

    if (l_write_to_file) {
      std::cout << "See *setup.txt file in output folder\n" << "\n";

      if (case_log.is_open())
        case_log.close();
    }

    if (l_write_to_file) {
      std::ofstream metadata(case_info_file, std::ios::app);
      append_command_output("git --no-pager diff", metadata);
    }

    timer_stop("write_case_info"); // In
  }

  // Allocate & initialize variables,
  // setup grid, setup constants, and setup flags
  // Setup grid
  setup_grid_api(nzmax, ngrdcol, sfc_elevation, l_implemented,         // In
                 l_ascending_grid, grid_type, deltaz, zm_init, zm_top, // In
                 momentum_heights, thermodynamic_heights,              // In
                 gr,                                                   // InOut
                 err_info);                                            // InOut

  if (clubb_at_least_debug_level_api(0)) {

    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      std::cerr << err_info.err_header_global << "\n";
      std::cerr << "Fatal error calling setup_grid_api in clubb_driver" << "\n";
      timer_stop("init_clubb_case");
      return;
    }
  }

  //  Allocate and initialize variables
  lmin                  = DriverArray1("lmin", ngrdcol);
  mixt_frac_max_mag     = DriverArray1("mixt_frac_max_mag", ngrdcol);
  um                    = DriverArray2("um", ngrdcol, gr.nzt);
  vm                    = DriverArray2("vm", ngrdcol, gr.nzt);
  upwp                  = DriverArray2("upwp", ngrdcol, gr.nzm);
  vpwp                  = DriverArray2("vpwp", ngrdcol, gr.nzm);
  up2                   = DriverArray2("up2", ngrdcol, gr.nzm);
  up3                   = DriverArray2("up3", ngrdcol, gr.nzt);
  vp2                   = DriverArray2("vp2", ngrdcol, gr.nzm);
  vp3                   = DriverArray2("vp3", ngrdcol, gr.nzt);
  thlm                  = DriverArray2("thlm", ngrdcol, gr.nzt);
  rtm                   = DriverArray2("rtm", ngrdcol, gr.nzt);
  wprtp                 = DriverArray2("wprtp", ngrdcol, gr.nzm);
  wpthlp                = DriverArray2("wpthlp", ngrdcol, gr.nzm);
  wprcp                 = DriverArray2("wprcp", ngrdcol, gr.nzm);
  w_up_in_cloud         = DriverArray2("w_up_in_cloud", ngrdcol, gr.nzt);
  w_down_in_cloud       = DriverArray2("w_down_in_cloud", ngrdcol, gr.nzt);
  cloudy_updraft_frac   = DriverArray2("cloudy_updraft_frac", ngrdcol, gr.nzt);
  cloudy_downdraft_frac = DriverArray2("cloudy_downdraft_frac", ngrdcol, gr.nzt);
  wp2                   = DriverArray2("wp2", ngrdcol, gr.nzm);
  wp3                   = DriverArray2("wp3", ngrdcol, gr.nzt);
  rtp2                  = DriverArray2("rtp2", ngrdcol, gr.nzm);
  thlp2                 = DriverArray2("thlp2", ngrdcol, gr.nzm);
  rtpthlp               = DriverArray2("rtpthlp", ngrdcol, gr.nzm);
  p_in_Pa               = DriverArray2("p_in_Pa", ngrdcol, gr.nzt);
  exner                 = DriverArray2("exner", ngrdcol, gr.nzt);
  rho                   = DriverArray2("rho", ngrdcol, gr.nzt);
  rho_zm                = DriverArray2("rho_zm", ngrdcol, gr.nzm);
  rho_ds_zm             = DriverArray2("rho_ds_zm", ngrdcol, gr.nzm);
  rho_ds_zt             = DriverArray2("rho_ds_zt", ngrdcol, gr.nzt);
  invrs_rho_ds_zm       = DriverArray2("invrs_rho_ds_zm", ngrdcol, gr.nzm);
  invrs_rho_ds_zt       = DriverArray2("invrs_rho_ds_zt", ngrdcol, gr.nzt);
  thv_ds_zm             = DriverArray2("thv_ds_zm", ngrdcol, gr.nzm);
  thv_ds_zt             = DriverArray2("thv_ds_zt", ngrdcol, gr.nzt);
  thlm_forcing          = DriverArray2("thlm_forcing", ngrdcol, gr.nzt);
  rtm_forcing           = DriverArray2("rtm_forcing", ngrdcol, gr.nzt);
  um_forcing            = DriverArray2("um_forcing", ngrdcol, gr.nzt);
  vm_forcing            = DriverArray2("vm_forcing", ngrdcol, gr.nzt);
  wprtp_forcing         = DriverArray2("wprtp_forcing", ngrdcol, gr.nzm);
  wpthlp_forcing        = DriverArray2("wpthlp_forcing", ngrdcol, gr.nzm);
  rtp2_forcing          = DriverArray2("rtp2_forcing", ngrdcol, gr.nzm);
  thlp2_forcing         = DriverArray2("thlp2_forcing", ngrdcol, gr.nzm);
  rtpthlp_forcing       = DriverArray2("rtpthlp_forcing", ngrdcol, gr.nzm);

  //  Variables used to track perturbed version of winds.
  um_pert   = DriverArray2("um_pert", ngrdcol, gr.nzt);
  vm_pert   = DriverArray2("vm_pert", ngrdcol, gr.nzt);
  upwp_pert = DriverArray2("upwp_pert", ngrdcol, gr.nzm);
  vpwp_pert = DriverArray2("vpwp_pert", ngrdcol, gr.nzm);

  //  Imposed large scale w
  wm_zm = DriverArray2("wm_zm", ngrdcol, gr.nzm);
  wm_zt = DriverArray2("wm_zt", ngrdcol, gr.nzt);

  //  Cloud water variables
  rcm               = DriverArray2("rcm", ngrdcol, gr.nzt);
  delta_zm          = DriverArray2("delta_zm", ngrdcol, gr.nzm);
  cloud_frac        = DriverArray2("cloud_frac", ngrdcol, gr.nzt);
  ice_supersat_frac = DriverArray2("ice_supersat_frac", ngrdcol, gr.nzt);
  rcm_in_layer      = DriverArray2("rcm_in_layer", ngrdcol, gr.nzt);
  cloud_cover       = DriverArray2("cloud_cover", ngrdcol, gr.nzt);
  invrs_tau_zm      = DriverArray2("invrs_tau_zm", ngrdcol, gr.nzm);

  //  Passive scalar variables
  //  Note that sclr_dim can be 0
  wpsclrp_sfc     = DriverArray2("wpsclrp_sfc", ngrdcol, sclr_dim);
  sclrm           = DriverArray3("sclrm", ngrdcol, gr.nzt, sclr_dim);
  sclrp2          = DriverArray3("sclrp2", ngrdcol, gr.nzm, sclr_dim);
  sclrp3          = DriverArray3("sclrp3", ngrdcol, gr.nzt, sclr_dim);
  sclrm_forcing   = DriverArray3("sclrm_forcing", ngrdcol, gr.nzt, sclr_dim);
  sclrprtp        = DriverArray3("sclrprtp", ngrdcol, gr.nzm, sclr_dim);
  sclrpthlp       = DriverArray3("sclrpthlp", ngrdcol, gr.nzm, sclr_dim);
  wpedsclrp_sfc   = DriverArray2("wpedsclrp_sfc", ngrdcol, edsclr_dim);
  edsclrm_forcing = DriverArray3("edsclrm_forcing", ngrdcol, gr.nzt, edsclr_dim);
  edsclrm         = DriverArray3("edsclrm", ngrdcol, gr.nzt, edsclr_dim);
  wpsclrp         = DriverArray3("wpsclrp", ngrdcol, gr.nzm, sclr_dim);
  sigma_sqd_w     = DriverArray2("sigma_sqd_w", ngrdcol, gr.nzm);
  sigma_sqd_w_zt  = DriverArray2("sigma_sqd_w_zt", ngrdcol, gr.nzt);
  Skw_zm_smooth   = DriverArray2("Skw_zm_smooth", ngrdcol, gr.nzm);
  wp2_zt          = DriverArray2("wp2_zt", ngrdcol, gr.nzt);
  ug              = DriverArray2("ug", ngrdcol, gr.nzt);
  vg              = DriverArray2("vg", ngrdcol, gr.nzt);
  um_ref          = DriverArray2("um_ref", ngrdcol, gr.nzt);
  vm_ref          = DriverArray2("vm_ref", ngrdcol, gr.nzt);
  thlm_ref        = DriverArray2("thlm_ref", ngrdcol, gr.nzt);
  rtm_ref         = DriverArray2("rtm_ref", ngrdcol, gr.nzt);
  thvm            = DriverArray2("thvm", ngrdcol, gr.nzt);
  radht           = DriverArray2("radht", ngrdcol, gr.nzt);
  thlprcp         = DriverArray2("thlprcp", ngrdcol, gr.nzm);
  thlp3           = DriverArray2("thlp3", ngrdcol, gr.nzt);
  rtp3            = DriverArray2("rtp3", ngrdcol, gr.nzt);

  //  Buoyancy related moments
  rtpthvp     = DriverArray2("rtpthvp", ngrdcol, gr.nzm);
  thlpthvp    = DriverArray2("thlpthvp", ngrdcol, gr.nzm);
  wpthvp      = DriverArray2("wpthvp", ngrdcol, gr.nzm);
  wp2thvp     = DriverArray2("wp2thvp", ngrdcol, gr.nzt);
  wp2up       = DriverArray2("wp2up", ngrdcol, gr.nzt);
  wp2rtp      = DriverArray2("wp2rtp", ngrdcol, gr.nzt);
  wp2thlp     = DriverArray2("wp2thlp", ngrdcol, gr.nzt);
  uprcp       = DriverArray2("uprcp", ngrdcol, gr.nzm);
  vprcp       = DriverArray2("vprcp", ngrdcol, gr.nzm);
  rc_coef_zm  = DriverArray2("rc_coef_zm", ngrdcol, gr.nzm);
  wp4         = DriverArray2("wp4", ngrdcol, gr.nzm);
  wpup2       = DriverArray2("wpup2", ngrdcol, gr.nzt);
  wpvp2       = DriverArray2("wpvp2", ngrdcol, gr.nzt);
  wp2up2      = DriverArray2("wp2up2", ngrdcol, gr.nzm);
  wp2vp2      = DriverArray2("wp2vp2", ngrdcol, gr.nzm);
  Kh_zt       = DriverArray2("Kh_zt", ngrdcol, gr.nzt);
  Kh_zm       = DriverArray2("Kh_zm", ngrdcol, gr.nzm);
  K_hm        = DriverArray3("K_hm", ngrdcol, gr.nzm, hydromet_dim);
  em          = DriverArray2("em", ngrdcol, gr.nzm);
  Lscale      = DriverArray2("Lscale", ngrdcol, gr.nzt);
  tau_zm      = DriverArray2("tau_zm", ngrdcol, gr.nzm);
  tau_zt      = DriverArray2("tau_zt", ngrdcol, gr.nzt);
  Nccnm       = DriverArray2("Nccnm", ngrdcol, gr.nzt);
  hydromet    = DriverArray3("hydromet", ngrdcol, gr.nzt, hydromet_dim);
  hydrometp2  = DriverArray3("hydrometp2", ngrdcol, gr.nzm, hydromet_dim);
  wphydrometp = DriverArray3("wphydrometp", ngrdcol, gr.nzm, hydromet_dim);
  Ncm         = DriverArray2("Ncm", ngrdcol, gr.nzt);
  wpNcp       = DriverArray2("wpNcp", ngrdcol, gr.nzm);

  //  High-order passive scalars
  sclrpthvp  = DriverArray3("sclrpthvp", ngrdcol, gr.nzm, sclr_dim);
  wp2hmp     = DriverArray3("wp2hmp", ngrdcol, gr.nzt, hydromet_dim);
  rtphmp_zt  = DriverArray3("rtphmp_zt", ngrdcol, gr.nzt, hydromet_dim);
  thlphmp_zt = DriverArray3("thlphmp_zt", ngrdcol, gr.nzt, hydromet_dim);

  //  Allocate rvm_mc, rcm_mc, thlm_mc
  rvm_mc  = DriverArray2("rvm_mc", ngrdcol, gr.nzt);
  rcm_mc  = DriverArray2("rcm_mc", ngrdcol, gr.nzt);
  thlm_mc = DriverArray2("thlm_mc", ngrdcol, gr.nzt);

  //  Allocate hydrometeor variables.
  rrm = DriverArray2("rrm", ngrdcol, gr.nzt);
  Nrm = DriverArray2("Nrm", ngrdcol, gr.nzt);

  //  Allocate hydromet_pdf_params
  hydromet_pdf_params.allocate(ngrdcol, gr.nzt);

  //  Allocate the correlation arrays
  corr_array_1_n      = DriverArray4("corr_array_1_n", ngrdcol, gr.nzt, pdf_dim, pdf_dim);
  corr_array_2_n      = DriverArray4("corr_array_2_n", ngrdcol, gr.nzt, pdf_dim, pdf_dim);
  corr_cholesky_mtx_1 = DriverArray4("corr_cholesky_mtx_1", ngrdcol, gr.nzt, pdf_dim, pdf_dim);
  corr_cholesky_mtx_2 = DriverArray4("corr_cholesky_mtx_2", ngrdcol, gr.nzt, pdf_dim, pdf_dim);

  //  Allocate the mean and stddev arrays
  mu_x_1_n    = DriverArray3("mu_x_1_n", ngrdcol, gr.nzt, pdf_dim);
  mu_x_2_n    = DriverArray3("mu_x_2_n", ngrdcol, gr.nzt, pdf_dim);
  sigma_x_1_n = DriverArray3("sigma_x_1_n", ngrdcol, gr.nzt, pdf_dim);
  sigma_x_2_n = DriverArray3("sigma_x_2_n", ngrdcol, gr.nzt, pdf_dim);

  //  Allocate microphysics tendencies for <w'rt'>, <w'thl'>,
  //  <rt'^2>, <thl'^2>, and <rt'thl'>.
  wprtp_mc   = DriverArray2("wprtp_mc", ngrdcol, gr.nzm);
  wpthlp_mc  = DriverArray2("wpthlp_mc", ngrdcol, gr.nzm);
  rtp2_mc    = DriverArray2("rtp2_mc", ngrdcol, gr.nzm);
  thlp2_mc   = DriverArray2("thlp2_mc", ngrdcol, gr.nzm);
  rtpthlp_mc = DriverArray2("rtpthlp_mc", ngrdcol, gr.nzm);

  //  Allocate microphysics tendencies for mean hydrometeors and mean cloud
  //  droplet concentration.
  hydromet_mc = DriverArray3("hydromet_mc", ngrdcol, gr.nzt, hydromet_dim);
  Ncm_mc      = DriverArray2("Ncm_mc", ngrdcol, gr.nzt);

  //  Allocate mean sedimentation velocity for hydrometeors
  hydromet_vel_zt = DriverArray3("hydromet_vel_zt", ngrdcol, gr.nzt, hydromet_dim);

  //  Allocate covariance sedimentation velocities for <V_hm'hm'>, implicit and
  //  explicit components, such that:
  //  <V_hm'hm'> = <V_hm'hm'>|_impc * <hm> + <V_hm'hm'>|_expc.
  hydromet_vel_covar_zt_impc =
      DriverArray3("hydromet_vel_covar_zt_impc", ngrdcol, gr.nzt, hydromet_dim);
  hydromet_vel_covar_zt_expc =
      DriverArray3("hydromet_vel_covar_zt_expc", ngrdcol, gr.nzt, hydromet_dim);
  rfrzm                = DriverArray2("rfrzm", ngrdcol, gr.nzt);
  X_nl_all_levs        = DriverArray4("X_nl_all_levs", ngrdcol, lh_num_samples, gr.nzt, pdf_dim);
  X_mixt_comp_all_levs = DriverInt3("X_mixt_comp_all_levs", ngrdcol, lh_num_samples, gr.nzt);
  lh_rt_clipped        = DriverArray3("lh_rt_clipped", ngrdcol, lh_num_samples, gr.nzt);
  lh_thl_clipped       = DriverArray3("lh_thl_clipped", ngrdcol, lh_num_samples, gr.nzt);
  lh_rc_clipped        = DriverArray3("lh_rc_clipped", ngrdcol, lh_num_samples, gr.nzt);
  lh_rv_clipped        = DriverArray3("lh_rv_clipped", ngrdcol, lh_num_samples, gr.nzt);
  lh_Nc_clipped        = DriverArray3("lh_Nc_clipped", ngrdcol, lh_num_samples, gr.nzt);
  lh_sample_point_weights =
      DriverArray3("lh_sample_point_weights", ngrdcol, lh_num_samples, gr.nzt);
  Nc_in_cloud          = DriverArray2("Nc_in_cloud", ngrdcol, gr.nzt);
  wpthlp_sfc           = DriverArray1("wpthlp_sfc", ngrdcol);
  wprtp_sfc            = DriverArray1("wprtp_sfc", ngrdcol);
  upwp_sfc             = DriverArray1("upwp_sfc", ngrdcol);
  vpwp_sfc             = DriverArray1("vpwp_sfc", ngrdcol);
  thlm_init            = DriverArray1("thlm_init", gr.nzt);
  rtm_init             = DriverArray1("rtm_init", gr.nzt);
  um_init              = DriverArray1("um_init", gr.nzt);
  vm_init              = DriverArray1("vm_init", gr.nzt);
  ug_init              = DriverArray1("ug_init", gr.nzt);
  vg_init              = DriverArray1("vg_init", gr.nzt);
  wp2_init             = DriverArray1("wp2_init", gr.nzm);
  up2_init             = DriverArray1("up2_init", gr.nzm);
  vp2_init             = DriverArray1("vp2_init", gr.nzm);
  upwp_init            = DriverArray1("upwp_init", gr.nzm);
  rcm_init             = DriverArray1("rcm_init", gr.nzt);
  wm_zt_init           = DriverArray1("wm_zt_init", gr.nzt);
  wm_zm_init           = DriverArray1("wm_zm_init", gr.nzm);
  em_init              = DriverArray1("em_init", gr.nzm);
  exner_init           = DriverArray1("exner_init", gr.nzt);
  thvm_init            = DriverArray1("thvm_init", gr.nzt);
  p_in_Pa_init         = DriverArray1("p_in_Pa_init", gr.nzt);
  p_in_Pa_zm_init      = DriverArray1("p_in_Pa_zm_init", gr.nzm);
  rho_init             = DriverArray1("rho_init", gr.nzt);
  rho_zm_init          = DriverArray1("rho_zm_init", gr.nzm);
  rho_ds_zm_init       = DriverArray1("rho_ds_zm_init", gr.nzm);
  rho_ds_zt_init       = DriverArray1("rho_ds_zt_init", gr.nzt);
  invrs_rho_ds_zm_init = DriverArray1("invrs_rho_ds_zm_init", gr.nzm);
  invrs_rho_ds_zt_init = DriverArray1("invrs_rho_ds_zt_init", gr.nzt);
  thv_ds_zm_init       = DriverArray1("thv_ds_zm_init", gr.nzm);
  thv_ds_zt_init       = DriverArray1("thv_ds_zt_init", gr.nzt);
  rtm_ref_init         = DriverArray1("rtm_ref_init", gr.nzt);
  thlm_ref_init        = DriverArray1("thlm_ref_init", gr.nzt);
  um_ref_init          = DriverArray1("um_ref_init", gr.nzt);
  vm_ref_init          = DriverArray1("vm_ref_init", gr.nzt);
  Ncm_init             = DriverArray1("Ncm_init", gr.nzt);
  Nc_in_cloud_init     = DriverArray1("Nc_in_cloud_init", gr.nzt);
  Nccnm_init           = DriverArray1("Nccnm_init", gr.nzt);
  sclrm_init           = DriverArray2("sclrm_init", gr.nzt, sclr_dim);
  edsclrm_init         = DriverArray2("edsclrm_init", gr.nzt, edsclr_dim);

  if (clubb_config_flags.l_add_dycore_grid) {
    // Initialize the dycore grid
    setup_gr_dycore(workspace, iunit, ngrdcol,                         // In
                    Kokkos::subview(gr.zm, Kokkos::ALL, 0),            // In
                    Kokkos::subview(gr.zm, Kokkos::ALL, (gr.nzm - 1)), // In
                    l_ascending_grid,                                  // In
                    gr_dycore,                                         // Out
                    err_info);                                         // InOut

    gr_output = clone_grid(gr_dycore);
  } else {
    gr_output = clone_grid(gr);
  }

  // This special purpose code only applies to tuner runs where the tune_type
  // is setup to try all permutations of our model flags
  if (model_flags_array) {
    clubb_config_flags.l_godunov_upwind_wpxp_ta = (*model_flags_array)(0);
    clubb_config_flags.l_godunov_upwind_xpyp_ta = (*model_flags_array)(1);
    clubb_config_flags.l_upwind_xpyp_ta         = (*model_flags_array)(2);
    clubb_config_flags.l_upwind_xm_ma           = (*model_flags_array)(3);
    clubb_config_flags.l_vert_avg_closure       = (*model_flags_array)(4);
    clubb_config_flags.l_standard_term_ta       = (*model_flags_array)(5);
    clubb_config_flags.l_use_cloud_cover        = (*model_flags_array)(6);
    clubb_config_flags.l_rcm_supersat_adj       = (*model_flags_array)(7);

    if (clubb_config_flags.l_vert_avg_closure) {
      clubb_config_flags.l_trapezoidal_rule_zt    = true;
      clubb_config_flags.l_trapezoidal_rule_zm    = true;
      clubb_config_flags.l_call_pdf_closure_twice = true;
    }
  }

  iunit = 50;

  //  Only output radiation files if using a radiation scheme
  if (lh_microphys_type != lh_microphys_disabled) {
    l_silhs_out = true;
  } else {
    l_silhs_out = false;
  }

  if (!l_restart) {
    time_current = time_initial;
  } else {
    time_current = time_restart;
  }

  // Time integration
  // Call advance_clubb_core once per each statistics output time
  ifinal =
      static_cast<int>(std::floor(((time_final - time_initial) / static_cast<double>(dt_main))));

  // Setup filenames and variables to set for setfields, if enabled
  if (l_input_fields || l_restart) {
    inputfields_init(iunit, runfile, day, month, year, // In
                     inputfields);                     // InOut
  }

  if (std::abs(((dt_rad / dt_main) - static_cast<real_t>(static_cast<int>(
                                         std::floor((dt_rad / dt_main)))))) > real_t(1.e-8)) {
    throw std::runtime_error("dt_rad must be a multiple of dt_main");
  }

  if (l_silhs_rad && clubb_config_flags.l_calc_thlp2_rad) {
    std::cerr << err_info.err_header_global << "\n";
    std::cerr << "The options l_silhs_rad and l_calc_thlp2_rad are incompatible."
              << "\n";

    //  General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    timer_stop("init_clubb_case");
    return;
  }

  if (clubb_config_flags.l_calc_thlp2_rad && (rad_scheme == "none")) {
    std::cerr << err_info.err_header_global << "\n";
    std::cerr << "The options rad_scheme == none and l_calc_thlp2_rad are incompatible."
              << "\n";

    //  General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    timer_stop("init_clubb_case");
    return;
  }

  if (l_soil_veg && (rad_scheme == "none")) {
    std::cerr << err_info.err_header_global << "\n";
    std::cerr << "The options l_soil_veg and rad_scheme == none are incompatible."
              << "\n";

    //  General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    timer_stop("init_clubb_case");
    return;
  }

  // Currently initialize_clubb does more than just read in the initial sounding.
  // It also includes other important initializations such as um_ref and vm_ref.
  // Therefore it should be executed prior to a restart. The restart should overwrite
  // the initial sounding anyway.
  // NOTE: initialize_clubb has the capability to be run on the whole 2D arrays we have
  // but there is no difference in the initializations between columns at the
  // moment, so we pass in a value of "1" for ngrdcol and use these _init
  // arrays that are then copied to the 2D versions. This reduces the time
  // spent in this routine and allows us to avoid large CPU to GPU data
  // transfers since we need only copy these 1D _init variables in.
  if (l_readiopdata) {
    readiopdata_dims(workspace, lat_vals, lon_vals, // Out
                     iop);                          // InOut
  }

  initialize_clubb(
      gr, 1, iunit, forcings_file_path,                                           // In
      p_sfc(0),                                                                   // InOut
      host_values({zm_init(0)}),                                                  // In
      sclr_dim, edsclr_dim, sclr_idx,                                             // In
      clubb_config_flags,                                                         // In
      l_modify_ic_with_cubic_int,                                                 // In
      l_ascending_grid, fcor_y, l_readiopdata,                                    // In
      as_columns(thlm_init), as_columns(rtm_init), as_columns(um_init),           // InOut
      as_columns(vm_init), as_columns(ug_init), as_columns(vg_init),              // InOut
      as_columns(wp2_init), as_columns(up2_init), as_columns(vp2_init),           // InOut
      as_columns(upwp_init), as_columns(rcm_init),                                // InOut
      as_columns(wm_zt_init), as_columns(wm_zm_init), as_columns(em_init),        // InOut
      as_columns(exner_init),                                                     // InOut
      as_columns(thvm_init), as_columns(p_in_Pa_init),                            // InOut
      p_in_Pa_zm_init,                                                            // Out
      as_columns(rho_init), as_columns(rho_zm_init), as_columns(rho_ds_zm_init),  // InOut
      as_columns(rho_ds_zt_init),                                                 // InOut
      as_columns(invrs_rho_ds_zm_init), as_columns(invrs_rho_ds_zt_init),         // InOut
      as_columns(thv_ds_zm_init), as_columns(thv_ds_zt_init),                     // InOut
      as_columns(rtm_ref_init), as_columns(thlm_ref_init),                        // InOut
      as_columns(um_ref_init), as_columns(vm_ref_init),                           // InOut
      as_columns(Ncm_init), as_columns(Nc_in_cloud_init), as_columns(Nccnm_init), // InOut
      err_info,                                                                   // InOut
      deep_soil_T_in_K_init, sfc_soil_T_in_K_init, veg_T_in_K_init,               // InOut
      as_columns(sclrm_init), as_columns(edsclrm_init));                          // Out

  if (clubb_at_least_debug_level_api(0)) {

    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      std::cerr << err_info.err_header_global << "\n";
      std::cerr << "FATAL ERROR calling initialize_clubb in run_clubb" << "\n";
      timer_stop("init_clubb_case");
      return;
    }
  }

  if (rad_scheme != "none") {
    setup_radiation_variables(ngrdcol, gr.nzm, gr.nzt, // In
                              radiation_variables);    // Out
  }

  if (l_stats) {
    //  Initialize statistics output; this allocates/initializes stats metadata
    //  for all columns but writes to one shared stats file.
    stats_registry = runfile;
    timer_start("stats_init_api");

    if ((clubb_config_flags.grid_adapt_in_time_method == no_grid_adaptation) &&
        (!clubb_config_flags.l_add_dycore_grid)) {
      StatsInitOptions stats_options;
      stats_options.clubb_params  = to_host_space(clubb_params_all);
      stats_options.param_names   = text_values(params_list);
      stats_options.ncol_total    = total_param_sets;
      stats_options.stats_tstart  = stats_tstart;
      stats_options.stats_tend    = stats_tend;
      stats_options.hydromet_list = text_values(hm_metadata.hydromet_list);
      stats_options.sclr_dim      = sclr_dim;
      stats_options.edsclr_dim    = edsclr_dim;
      stats_init_api(
          stats_registry, stats_output_path, ngrdcol, stats_tsamp, stats_tout, dt_main, day, month,
          year, time_initial, to_host_space(Kokkos::subview(gr.zt, 0, Kokkos::ALL)),
          to_host_space(Kokkos::subview(gr.zm, 0, Kokkos::ALL)), stats, err_info, stats_options);
    } else {
      StatsInitOptions stats_options;
      stats_options.clubb_params  = to_host_space(clubb_params_all);
      stats_options.param_names   = text_values(params_list);
      stats_options.ncol_total    = total_param_sets;
      stats_options.stats_tstart  = stats_tstart;
      stats_options.stats_tend    = stats_tend;
      stats_options.hydromet_list = text_values(hm_metadata.hydromet_list);
      stats_options.sclr_dim      = sclr_dim;
      stats_options.edsclr_dim    = edsclr_dim;
      stats_options.output_zt     = to_host_space(Kokkos::subview(gr_output.zt, 0, Kokkos::ALL));
      stats_options.output_zm     = to_host_space(Kokkos::subview(gr_output.zm, 0, Kokkos::ALL));
      stats_options.grid_remap_method = clubb_config_flags.grid_remap_method;
      stats_init_api(
          stats_registry, stats_output_path, ngrdcol, stats_tsamp, stats_tout, dt_main, day, month,
          year, time_initial, to_host_space(Kokkos::subview(gr.zt, 0, Kokkos::ALL)),
          to_host_space(Kokkos::subview(gr.zm, 0, Kokkos::ALL)), stats, err_info, stats_options);
    }

    timer_stop("stats_init_api"); // In

    if (clubb_at_least_debug_level_api(0)) {

      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        std::cerr << err_info.err_header_global << "\n";
        std::cerr << "FATAL ERROR in stats_init_api" << "\n";
        timer_stop("init_clubb_case");
        return;
      }
    }
  }

  if (((total_param_sets > batch_size) && stats.l_netcdf_output) && l_silhs_out) {
    std::cerr << err_info.err_header_global << "\n";
    std::cerr << "Batch-mode stats NetCDF output does not yet support SILHS sample output."
              << "\n";

    //  General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    timer_stop("init_clubb_case");
    return;
  }

  if (lh_microphys_type != lh_microphys_disabled) {
    // Setup 2D output of all subcolumns (if enabled)
    latin_hypercube_2D_output_api(workspace, gr.nzt, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
                                  pdf_dim, lh_num_samples, // In
                                  hm_metadata,             // In
                                  &stats,                  // Optional InOut
                                  err_info);               // InOut
  }

  if (any_equal(err_info.err_code, clubb_fatal_error)) {
    std::cerr << err_info.err_header_global << "\n";
    std::cerr << "Fatal error calling latin_hypercube_2D_output_api in run_clubb"
              << "\n";
    timer_stop("init_clubb_case");
    return;
  }

  if (l_stats) {

    if (!((std::abs(((dt_rad / stats_tout) - static_cast<real_t>(static_cast<int>(std::floor(
                                                 (dt_rad / stats_tout)))))) < real_t(1.e-8)) ||
          (std::abs(((stats_tout / dt_rad) - static_cast<real_t>(static_cast<int>(std::floor(
                                                 (stats_tout / dt_rad)))))) < real_t(1.e-8)))) {
      throw std::runtime_error(
          "dt_rad must be a multiple of stats_tout or stats_tout must be a mulitple of dt_rad");
    }
  }

  stats_nsamp = static_cast<int>(std::round((stats_tsamp / dt_main)));
  stats_nout  = static_cast<int>(std::round((stats_tout / dt_main)));
  init_pdf_params_api(gr.nzt, ngrdcol, pdf_params);
  init_pdf_params_api(gr.nzm, ngrdcol, pdf_params_zm);

  if ((clubb_config_flags.iiPDF_type == iiPDF_new) ||
      (clubb_config_flags.iiPDF_type == iiPDF_new_hybrid)) {
    init_pdf_implicit_coefs_terms_api(gr.nzt, ngrdcol, sclr_dim, pdf_implicit_coefs_terms);
  }

  if (!(microphys_scheme == "none")) {
    init_precip_fracs_api(gr.nzt, ngrdcol, precip_fracs);
  }

  // Calculate CLUBB parameters that should be derived from other quantities
  calc_derived_params_api(gr, ngrdcol, grid_type, deltaz,             // In
                          clubb_params,                               // In
                          clubb_config_flags.l_prescribed_avg_deltaz, // In
                          nu_vert_res_dep, lmin,                      // InOut
                          mixt_frac_max_mag);                         // InOut

  set_case_initial_conditions(err_info, nullptr, std::nullopt);
  timer_stop("init_clubb_case");
}

void ClubbDriver::set_case_initial_conditions(err_info_type& err_info,             // InOut
                                              const DriverArray2* clubb_params_in, // Optional In
                                              std::optional<int> batch_num) {      // Optional In

  // Restore initial profiles, physics caches and sampling state for a rerun or batch.
  auto& rad_scheme       = parameters_radiation.rad_scheme;
  auto& microphys_scheme = microphys.microphys_scheme;
  auto& lh_num_samples   = microphys.lh_num_samples;

  // Local variables.
  int i{};
  int k{};
  int sclr{};
  int edsclr{};
  int sample{};
  int hm{};
  int batch_start{}; // Zero-based first column in the requested runtime batch.
  int batch_end{};   // Exclusive upper column bound in the requested runtime batch.
  int num_batches{}; // Total number of runtime batches available.

  // ----------------------------------- Begin Code -----------------------------------
  timer_start("set_case_initial_conditions");
  stats_reset_api(stats);

  if (execution_fields_initialized) {
    // Preserve entries retained by inactive reset branches before host writes.
    for_each_execution_field([](auto& field, const auto&) { field.sync_host(); });
  }

  if (batch_num.has_value()) {
    //  init_clubb_case validates that total_param_sets is evenly divisible by
    //  batch_size.  A no-batching run is therefore just one batch whose width is
    //  total_param_sets.  A batched run reuses the same allocated runtime arrays
    //  for each slice of clubb_params_all.
    num_batches = ((total_param_sets / batch_size));

    if (((*batch_num) < 1) || ((*batch_num) > num_batches)) {
      //  General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      std::cerr << "set_case_initial_conditions batch_num is outside the valid runtime batch range"
                << "\n";
      timer_stop("set_case_initial_conditions");
      return;
    }

    if ((*batch_num) > 1) {
      // Set stats to output the next batch
      start_next_stats_batch_api(stats, err_info);

      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        timer_stop("set_case_initial_conditions");
        return;
      }
    }
  }

  // This should be GPUized. Preserve any newly supplied host-side parameter
  // matrix instead of overwriting it from the device before recomputing the
  // parameter-dependent derived state for the next run.
  if (clubb_params_in) {

    if (((*clubb_params_in).extent_int(0) != ngrdcol) ||
        ((*clubb_params_in).extent_int(1) != nparams)) {
      //  General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      std::cerr << "set_case_initial_conditions requires clubb_params_in to match the initialized"
                << " " << " (ngrdcol,nparams) shape" << "\n";
      timer_stop("set_case_initial_conditions");
      return;
    }

    for (int iparam = 0; iparam < clubb_params.extent_int(1); ++iparam) {

      for (int icol = 0; icol < clubb_params.extent_int(0); ++icol) {
        clubb_params(icol, iparam) = (*clubb_params_in)(icol, iparam);
      }
    }
  } else if (batch_num.has_value()) {
    batch_start = ((*batch_num) - 1) * batch_size;
    batch_end   = (batch_size + batch_start);
    Kokkos::deep_copy(
        clubb_params,
        Kokkos::subview(clubb_params_all, std::make_pair(batch_start, batch_end), Kokkos::ALL));
  }

  if (rad_scheme != "none") {
    reset_radiation_variables(radiation_variables); // InOut
  }

  // These variables are not used on the GPU and the commented out ones
  // are not used at all beyond being initialized
  // sigma_sqd_w_zt         = zero        ! PDF width parameter interp. to t-levs.
  // sigma_sqd_w            = zero          ! PDF width parameter (momentum levels)
  // tau_zt                 = zero        ! Eddy dissipation time scale: thermo. levels
  // tau_zm                 = zero          ! Eddy dissipation time scale: momentum levels
  for (int k = 0; k < Ncm_mc.extent_int(1); ++k) {

    for (int icol = 0; icol < Ncm_mc.extent_int(0); ++icol) {
      Ncm_mc(icol, k) = zero;
    }
  }

  for (int k = 0; k < wpNcp.extent_int(1); ++k) {

    for (int icol = 0; icol < wpNcp.extent_int(0); ++icol) {
      wpNcp(icol, k) = zero;
    }
  }

  for (k = 0; k < gr.nzt; ++k) {

    for (i = 0; i < ngrdcol; ++i) {
      Nccnm(i, k)       = Nccnm_init(k);
      Ncm(i, k)         = Ncm_init(k);
      Nc_in_cloud(i, k) = Nc_in_cloud_init(k);
    }
  }

  for (sample = 0; sample < lh_num_samples; ++sample) {

    for (k = 0; k < gr.nzt; ++k) {

      for (i = 0; i < ngrdcol; ++i) {

        for (int pdf_var = 0; pdf_var < X_nl_all_levs.extent_int(3); ++pdf_var) {
          X_nl_all_levs(i, sample, k, pdf_var) = (-real_t(999.));
        }

        X_mixt_comp_all_levs(i, sample, k)    = (-999);
        lh_sample_point_weights(i, sample, k) = (-real_t(999.));
      }
    }
  }

  for (k = 0; k < gr.nzt; ++k) {

    for (i = 0; i < ngrdcol; ++i) {
      um(i, k)                    = um_init(k);
      vm(i, k)                    = vm_init(k);
      up3(i, k)                   = zero;
      vp3(i, k)                   = zero;
      thlm(i, k)                  = thlm_init(k);
      rtm(i, k)                   = rtm_init(k);
      w_up_in_cloud(i, k)         = zero;
      w_down_in_cloud(i, k)       = zero;
      cloudy_updraft_frac(i, k)   = zero;
      cloudy_downdraft_frac(i, k) = zero;
      wp3(i, k)                   = zero;
      p_in_Pa(i, k)               = p_in_Pa_init(k);
      exner(i, k)                 = exner_init(k);
      rho(i, k)                   = rho_init(k);
      rho_ds_zt(i, k)             = rho_ds_zt_init(k);
      invrs_rho_ds_zt(i, k)       = invrs_rho_ds_zt_init(k);
      thv_ds_zt(i, k)             = thv_ds_zt_init(k);
      thlm_forcing(i, k)          = zero;
      rtm_forcing(i, k)           = zero;
      um_forcing(i, k)            = zero;
      vm_forcing(i, k)            = zero;
      um_pert(i, k)               = zero;
      vm_pert(i, k)               = zero;
      wm_zt(i, k)                 = wm_zt_init(k);
      rcm(i, k)                   = rcm_init(k);
      cloud_frac(i, k)            = zero;
      ice_supersat_frac(i, k)     = zero;
      rcm_in_layer(i, k)          = zero;
      cloud_cover(i, k)           = zero;
      ug(i, k)                    = ug_init(k);
      vg(i, k)                    = vg_init(k);
      um_ref(i, k)                = um_ref_init(k);
      vm_ref(i, k)                = vm_ref_init(k);
      thlm_ref(i, k)              = thlm_ref_init(k);
      rtm_ref(i, k)               = rtm_ref_init(k);
      thvm(i, k)                  = thvm_init(k);
      radht(i, k)                 = zero;
      thlp3(i, k)                 = zero;
      rtp3(i, k)                  = zero;
      wp2thvp(i, k)               = zero;
      wp2up(i, k)                 = zero;
      wp2rtp(i, k)                = zero;
      wp2thlp(i, k)               = zero;
      wpup2(i, k)                 = zero;
      wpvp2(i, k)                 = zero;
      Kh_zt(i, k)                 = zero;
      Lscale(i, k)                = zero;
      rvm_mc(i, k)                = zero;
      rcm_mc(i, k)                = zero;
      thlm_mc(i, k)               = zero;
      rfrzm(i, k)                 = zero;
    }
  }

  for (k = 0; k < gr.nzm; ++k) {

    for (i = 0; i < ngrdcol; ++i) {
      upwp(i, k)            = upwp_init(k);
      vpwp(i, k)            = zero;
      up2(i, k)             = up2_init(k);
      vp2(i, k)             = vp2_init(k);
      wprtp(i, k)           = zero;
      wpthlp(i, k)          = zero;
      wprcp(i, k)           = zero;
      wp2(i, k)             = wp2_init(k);
      rtp2(i, k)            = std::pow(rt_tol, 2);
      thlp2(i, k)           = std::pow(thl_tol, 2);
      rtpthlp(i, k)         = zero;
      rho_zm(i, k)          = rho_zm_init(k);
      rho_ds_zm(i, k)       = rho_ds_zm_init(k);
      invrs_rho_ds_zm(i, k) = invrs_rho_ds_zm_init(k);
      thv_ds_zm(i, k)       = thv_ds_zm_init(k);
      wprtp_forcing(i, k)   = zero;
      wpthlp_forcing(i, k)  = zero;
      rtp2_forcing(i, k)    = zero;
      thlp2_forcing(i, k)   = zero;
      rtpthlp_forcing(i, k) = zero;
      upwp_pert(i, k)       = zero;
      vpwp_pert(i, k)       = zero;
      wm_zm(i, k)           = wm_zm_init(k);
      invrs_tau_zm(i, k)    = zero;
      thlprcp(i, k)         = zero;
      rtpthvp(i, k)         = zero;
      thlpthvp(i, k)        = zero;
      wpthvp(i, k)          = zero;
      uprcp(i, k)           = zero;
      vprcp(i, k)           = zero;
      rc_coef_zm(i, k)      = zero;
      wp4(i, k)             = zero;
      wp2up2(i, k)          = zero;
      wp2vp2(i, k)          = zero;
      Kh_zm(i, k)           = zero;
      wprtp_mc(i, k)        = zero;
      wpthlp_mc(i, k)       = zero;
      rtp2_mc(i, k)         = zero;
      thlp2_mc(i, k)        = zero;
      rtpthlp_mc(i, k)      = zero;
    }
  }

  for (i = 0; i < ngrdcol; ++i) {
    wpthlp_sfc(i)       = zero;
    wprtp_sfc(i)        = zero;
    upwp_sfc(i)         = zero;
    vpwp_sfc(i)         = zero;
    deep_soil_T_in_K(i) = deep_soil_T_in_K_init(0);
    sfc_soil_T_in_K(i)  = sfc_soil_T_in_K_init(0);
    veg_T_in_K(i)       = veg_T_in_K_init(0);
  }

  // Passive scalars
  if (sclr_dim > 0) {

    for (sclr = 0; sclr < sclr_dim; ++sclr) {

      for (k = 0; k < gr.nzt; ++k) {

        for (i = 0; i < ngrdcol; ++i) {
          sclrm(i, k, sclr)         = sclrm_init(k, sclr);
          sclrm_forcing(i, k, sclr) = zero;
          sclrp3(i, k, sclr)        = zero;
        }
      }
    }

    for (sclr = 0; sclr < sclr_dim; ++sclr) {

      for (k = 0; k < gr.nzm; ++k) {

        for (i = 0; i < ngrdcol; ++i) {
          sclrpthvp(i, k, sclr) = zero;
          sclrp2(i, k, sclr)    = std::pow(sclr_tol(sclr), 2);
          sclrprtp(i, k, sclr)  = zero;
          sclrpthlp(i, k, sclr) = zero;
          wpsclrp(i, k, sclr)   = zero;
        }
      }
    }

    for (sclr = 0; sclr < sclr_dim; ++sclr) {

      for (i = 0; i < ngrdcol; ++i) {
        wpsclrp_sfc(i, sclr) = zero;
      }
    }
  }

  if (edsclr_dim > 0) {

    for (edsclr = 0; edsclr < edsclr_dim; ++edsclr) {

      for (k = 0; k < gr.nzt; ++k) {

        for (i = 0; i < ngrdcol; ++i) {
          edsclrm(i, k, edsclr)         = edsclrm_init(k, edsclr);
          edsclrm_forcing(i, k, edsclr) = zero;
        }
      }
    }

    for (edsclr = 0; edsclr < edsclr_dim; ++edsclr) {

      for (i = 0; i < ngrdcol; ++i) {
        wpedsclrp_sfc(i, edsclr) = zero;
      }
    }
  }

  if (hydromet_dim > 0) {

    for (hm = 0; hm < hydromet_dim; ++hm) {

      for (k = 0; k < gr.nzt; ++k) {

        for (i = 0; i < ngrdcol; ++i) {
          thlphmp_zt(i, k, hm) = zero;
          rtphmp_zt(i, k, hm)  = zero;
          wp2hmp(i, k, hm)     = zero;
        }
      }
    }

    for (hm = 0; hm < hydromet_dim; ++hm) {

      for (k = 0; k < gr.nzm; ++k) {

        for (i = 0; i < ngrdcol; ++i) {
          wphydrometp(i, k, hm) = zero;
        }
      }
    }

    // These are not used in any GPU code yet
    for (int hydromet_index = 0; hydromet_index < K_hm.extent_int(2); ++hydromet_index) {

      for (int k = 0; k < K_hm.extent_int(1); ++k) {

        for (int icol = 0; icol < K_hm.extent_int(0); ++icol) {
          K_hm(icol, k, hydromet_index) = zero;
        }
      }
    }

    for (int hydromet_index = 0; hydromet_index < hydromet.extent_int(2); ++hydromet_index) {

      for (int k = 0; k < hydromet.extent_int(1); ++k) {

        for (int icol = 0; icol < hydromet.extent_int(0); ++icol) {
          hydromet(icol, k, hydromet_index) = zero;
        }
      }
    }

    for (int hydromet_index = 0; hydromet_index < hydromet_mc.extent_int(2); ++hydromet_index) {

      for (int k = 0; k < hydromet_mc.extent_int(1); ++k) {

        for (int icol = 0; icol < hydromet_mc.extent_int(0); ++icol) {
          hydromet_mc(icol, k, hydromet_index) = zero;
        }
      }
    }

    for (int hydromet_index = 0; hydromet_index < hydrometp2.extent_int(2); ++hydromet_index) {

      for (int k = 0; k < hydrometp2.extent_int(1); ++k) {

        for (int icol = 0; icol < hydrometp2.extent_int(0); ++icol) {
          hydrometp2(icol, k, hydromet_index) = zero;
        }
      }
    }

    for (int hydromet_index = 0; hydromet_index < hydromet_vel_zt.extent_int(2); ++hydromet_index) {

      for (int k = 0; k < hydromet_vel_zt.extent_int(1); ++k) {

        for (int icol = 0; icol < hydromet_vel_zt.extent_int(0); ++icol) {
          hydromet_vel_zt(icol, k, hydromet_index) = zero;
        }
      }
    }

    for (int hydromet_index = 0; hydromet_index < hydromet_vel_covar_zt_impc.extent_int(2);
         ++hydromet_index) {

      for (int k = 0; k < hydromet_vel_covar_zt_impc.extent_int(1); ++k) {

        for (int icol = 0; icol < hydromet_vel_covar_zt_impc.extent_int(0); ++icol) {
          hydromet_vel_covar_zt_impc(icol, k, hydromet_index) = zero;
        }
      }
    }

    for (int hydromet_index = 0; hydromet_index < hydromet_vel_covar_zt_expc.extent_int(2);
         ++hydromet_index) {

      for (int k = 0; k < hydromet_vel_covar_zt_expc.extent_int(1); ++k) {

        for (int icol = 0; icol < hydromet_vel_covar_zt_expc.extent_int(0); ++icol) {
          hydromet_vel_covar_zt_expc(icol, k, hydromet_index) = zero;
        }
      }
    }
  }

  // Zero out our data structures
  zero_pdf_params_api(pdf_params);
  zero_pdf_params_api(pdf_params_zm);

  if ((clubb_config_flags.iiPDF_type == iiPDF_new) ||
      (clubb_config_flags.iiPDF_type == iiPDF_new_hybrid)) {
    zero_pdf_implicit_coefs_terms_api(pdf_implicit_coefs_terms);
  }

  if (!(microphys_scheme == "none")) {
    zero_precip_fracs_api(precip_fracs);
    init_hydromet_pdf_params(hydromet_pdf_params);
  }

  if (!l_restart) {
    time_current = time_initial;
    iinit        = 1;
  } else {
    time_current = time_restart;

    // Determining what iteration to restart at.
    // The value is increased by 1 to sychronize with restart data.
    // Joshua Fasching February 2008
    // Ensure that iteration num, iinit, is an integer, so that model time is
    // incremented correctly by iteration number at end of timestep
    if (std::abs(std::fmod((time_restart - time_initial), static_cast<double>(dt_main))) >
        static_cast<double>(eps)) {
      std::cerr << "Error: (time_restart-time_initial) " << " "
                << "is not a multiple of dt_main." << "\n";
      std::cerr << "time_restart = " << " " << time_restart << "\n";
      std::cerr << "time_initial = " << " " << time_initial << "\n";
      std::cerr << "dt_main = " << " " << dt_main << "\n";
      throw std::runtime_error("Fatal error");
    }

    iinit = (static_cast<int>(
                 std::floor(((time_current - time_initial) / static_cast<double>(dt_main)))) +
             1);

    for (i = 0; i < ngrdcol; ++i) {

      restart_clubb(gr,                                                        // In
                    hydromet_dim,                                              // In
                    hm_metadata,                                               // In
                    restart_path_case,                                         // In
                    time_restart,                                              // In
                    Kokkos::subview(um, i, Kokkos::ALL),                       // InOut
                    Kokkos::subview(upwp, i, Kokkos::ALL),                     // InOut
                    Kokkos::subview(vm, i, Kokkos::ALL),                       // InOut
                    Kokkos::subview(vpwp, i, Kokkos::ALL),                     // InOut
                    Kokkos::subview(up2, i, Kokkos::ALL),                      // InOut
                    Kokkos::subview(vp2, i, Kokkos::ALL),                      // InOut
                    Kokkos::subview(rtm, i, Kokkos::ALL),                      // InOut
                    Kokkos::subview(wprtp, i, Kokkos::ALL),                    // InOut
                    Kokkos::subview(thlm, i, Kokkos::ALL),                     // InOut
                    Kokkos::subview(wpthlp, i, Kokkos::ALL),                   // InOut
                    Kokkos::subview(rtp2, i, Kokkos::ALL),                     // InOut
                    Kokkos::subview(rtp3, i, Kokkos::ALL),                     // InOut
                    Kokkos::subview(thlp2, i, Kokkos::ALL),                    // InOut
                    Kokkos::subview(thlp3, i, Kokkos::ALL),                    // InOut
                    Kokkos::subview(rtpthlp, i, Kokkos::ALL),                  // InOut
                    Kokkos::subview(wp2, i, Kokkos::ALL),                      // InOut
                    Kokkos::subview(wp3, i, Kokkos::ALL),                      // InOut
                    Kokkos::subview(p_in_Pa, i, Kokkos::ALL),                  // InOut
                    Kokkos::subview(exner, i, Kokkos::ALL),                    // InOut
                    Kokkos::subview(rcm, i, Kokkos::ALL),                      // InOut
                    Kokkos::subview(cloud_frac, i, Kokkos::ALL),               // InOut
                    Kokkos::subview(wpthvp, i, Kokkos::ALL),                   // InOut
                    Kokkos::subview(wp2thvp, i, Kokkos::ALL),                  // InOut
                    Kokkos::subview(wp2up, i, Kokkos::ALL),                    // InOut
                    Kokkos::subview(rtpthvp, i, Kokkos::ALL),                  // InOut
                    Kokkos::subview(thlpthvp, i, Kokkos::ALL),                 // InOut
                    Kokkos::subview(wp2rtp, i, Kokkos::ALL),                   // InOut
                    Kokkos::subview(wp2thlp, i, Kokkos::ALL),                  // InOut
                    Kokkos::subview(uprcp, i, Kokkos::ALL),                    // InOut
                    Kokkos::subview(vprcp, i, Kokkos::ALL),                    // InOut
                    Kokkos::subview(rc_coef_zm, i, Kokkos::ALL),               // InOut
                    Kokkos::subview(wp4, i, Kokkos::ALL),                      // InOut
                    Kokkos::subview(wpup2, i, Kokkos::ALL),                    // InOut
                    Kokkos::subview(wpvp2, i, Kokkos::ALL),                    // InOut
                    Kokkos::subview(wp2up2, i, Kokkos::ALL),                   // InOut
                    Kokkos::subview(wp2vp2, i, Kokkos::ALL),                   // InOut
                    Kokkos::subview(ice_supersat_frac, i, Kokkos::ALL),        // InOut
                    Kokkos::subview(wm_zt, i, Kokkos::ALL),                    // InOut
                    Kokkos::subview(rho, i, Kokkos::ALL),                      // InOut
                    Kokkos::subview(rho_zm, i, Kokkos::ALL),                   // InOut
                    Kokkos::subview(rho_ds_zm, i, Kokkos::ALL),                // InOut
                    Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),                // InOut
                    Kokkos::subview(thv_ds_zm, i, Kokkos::ALL),                // InOut
                    Kokkos::subview(thv_ds_zt, i, Kokkos::ALL),                // InOut
                    Kokkos::subview(thlm_forcing, i, Kokkos::ALL),             // InOut
                    Kokkos::subview(rtm_forcing, i, Kokkos::ALL),              // InOut
                    Kokkos::subview(wprtp_forcing, i, Kokkos::ALL),            // InOut
                    Kokkos::subview(wpthlp_forcing, i, Kokkos::ALL),           // InOut
                    Kokkos::subview(rtp2_forcing, i, Kokkos::ALL),             // InOut
                    Kokkos::subview(thlp2_forcing, i, Kokkos::ALL),            // InOut
                    Kokkos::subview(rtpthlp_forcing, i, Kokkos::ALL),          // InOut
                    Kokkos::subview(hydromet, i, Kokkos::ALL, Kokkos::ALL),    // InOut
                    Kokkos::subview(hydrometp2, i, Kokkos::ALL, Kokkos::ALL),  // InOut
                    Kokkos::subview(wphydrometp, i, Kokkos::ALL, Kokkos::ALL), // InOut
                    Kokkos::subview(Ncm, i, Kokkos::ALL),                      // InOut
                    Kokkos::subview(Nccnm, i, Kokkos::ALL),                    // InOut
                    Kokkos::subview(thvm, i, Kokkos::ALL),                     // InOut
                    Kokkos::subview(em, i, Kokkos::ALL),                       // InOut
                    Kokkos::subview(tau_zm, i, Kokkos::ALL),                   // InOut
                    Kokkos::subview(tau_zt, i, Kokkos::ALL),                   // InOut
                    Kokkos::subview(Kh_zt, i, Kokkos::ALL),                    // InOut
                    Kokkos::subview(Kh_zm, i, Kokkos::ALL),                    // InOut
                    Kokkos::subview(ug, i, Kokkos::ALL),                       // InOut
                    Kokkos::subview(vg, i, Kokkos::ALL),                       // InOut
                    Kokkos::subview(thlprcp, i, Kokkos::ALL),                  // InOut
                    Kokkos::subview(sigma_sqd_w, i, Kokkos::ALL),              // InOut
                    Kokkos::subview(sigma_sqd_w_zt, i, Kokkos::ALL),           // InOut
                    Kokkos::subview(radht, i, Kokkos::ALL),                    // InOut
                    deep_soil_T_in_K(i),                                       // InOut
                    sfc_soil_T_in_K(i),                                        // InOut
                    veg_T_in_K(i),                                             // InOut
                    pdf_params,                                                // InOut
                    pdf_params_zm,                                             // InOut
                    Kokkos::subview(rcm_mc, i, Kokkos::ALL),                   // Out
                    Kokkos::subview(rvm_mc, i, Kokkos::ALL),                   // Out
                    Kokkos::subview(thlm_mc, i, Kokkos::ALL),                  // Out
                    Kokkos::subview(wprtp_mc, i, Kokkos::ALL),                 // Out
                    Kokkos::subview(wpthlp_mc, i, Kokkos::ALL),                // Out
                    Kokkos::subview(rtp2_mc, i, Kokkos::ALL),                  // Out
                    Kokkos::subview(thlp2_mc, i, Kokkos::ALL),                 // Out
                    Kokkos::subview(rtpthlp_mc, i, Kokkos::ALL),               // Out
                    wpthlp_sfc(i),                                             // Out
                    wprtp_sfc(i),                                              // Out
                    upwp_sfc(i),                                               // Out
                    vpwp_sfc(i));                                              // Out
    }

    // Calculate invrs_rho_ds_zm and invrs_rho_ds_zt from the values of
    // rho_ds_zm and rho_ds_zt, respectively, which were read in from the input
    // file during the call to subroutine restart_clubb.
    for (int k = 0; k < invrs_rho_ds_zm.extent_int(1); ++k) {

      for (int icol = 0; icol < invrs_rho_ds_zm.extent_int(0); ++icol) {
        invrs_rho_ds_zm(icol, k) = (real_t(1.0) / rho_ds_zm(icol, k));
      }
    }

    for (int k = 0; k < invrs_rho_ds_zt.extent_int(1); ++k) {

      for (int icol = 0; icol < invrs_rho_ds_zt.extent_int(0); ++icol) {
        invrs_rho_ds_zt(icol, k) = (real_t(1.0) / rho_ds_zt(icol, k));
      }
    }
  }

  // Host-resident Views need no native OpenACC host refresh before this call.
  //  Calculate CLUBB parameters that should be derived from other quantities
  calc_derived_params_api(gr, ngrdcol, grid_type, deltaz,             // In
                          clubb_params,                               // In
                          clubb_config_flags.l_prescribed_avg_deltaz, // In
                          nu_vert_res_dep, lmin,                      // InOut
                          mixt_frac_max_mag);                         // InOut

  // Re-run the parameter sanity checks when debugging so reruns with updated
  // parameter sets are validated the same way as the initial setup.
  if (clubb_at_least_debug_level_api(0)) {

    check_clubb_settings_api(ngrdcol,                                                 // In
                             edsclr_dim,                                              // In
                             sclr_idx,                                                // In
                             clubb_params,                                            // In
                             l_implemented,                                           // In
                             l_input_fields,                                          // In
                             clubb_config_flags,                                      // In
                             err_info,                                                // InOut
                             thlm_sponge_damp_settings, rtm_sponge_damp_settings,     // In
                             uv_sponge_damp_settings, wp2_sponge_damp_settings,       // In
                             wp3_sponge_damp_settings, up2_vp2_sponge_damp_settings); // In

    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      std::cerr << err_info.err_header_global << "\n";
      std::cerr << "Fatal error calling check_clubb_settings_api" << "\n";
      std::cerr << "in set_case_initial_conditions" << "\n";
      timer_stop("set_case_initial_conditions");
      return;
    }
  }

  // Match Native case data lifetime: allocate once and refresh initial fields
  // for this reset/batch. Existing execution allocations survive later resets.
  for_each_execution_field([](auto& field, const auto& host) {
    field.bind(host);
    field.modify_host();
    field.sync_device();
  });
  execution_fields_initialized = true;

  timer_stop("set_case_initial_conditions"); // In
}

void ClubbDriver::advance_clubb_to_end(bool l_stdout,                  // In
                                       err_info_type& err_info,        // InOut
                                       bool l_suppress_stats,          // Optional In
                                       std::optional<int> itime_start, // Optional In
                                       std::optional<int> itime_end) { // Optional In

  // Advance forcings, CLUBB core, PDF, microphysics, radiation and stats in native order.
  auto& rad_scheme        = parameters_radiation.rad_scheme;
  auto& microphys_scheme  = microphys.microphys_scheme;
  auto& sigma_g           = microphys.sigma_g;
  auto& l_cloud_sed       = microphys.l_cloud_sed;
  auto& lh_microphys_type = microphys.lh_microphys_type;
  auto& lh_num_samples    = microphys.lh_num_samples;

  // Runtime calls write execution fields; host refreshes occur only on reset.
  for_each_execution_field([](auto& field, const auto&) { field.modify_device(); });

  // Bind model names to the persistent execution Views used by this advance.
  // These handle copies retain allocations without transferring field data.
  const auto clubb_params               = execution.clubb_params.execution_view();
  const auto dummy_dx                   = execution.dummy_dx.execution_view();
  const auto dummy_dy                   = execution.dummy_dy.execution_view();
  const auto sclr_tol                   = execution.sclr_tol.execution_view();
  const auto p_sfc                      = execution.p_sfc.execution_view();
  const auto T_sfc                      = execution.T_sfc.execution_view();
  const auto fcor                       = execution.fcor.execution_view();
  const auto fcor_y                     = execution.fcor_y.execution_view();
  const auto sfc_elevation              = execution.sfc_elevation.execution_view();
  const auto corr_array_n_cloud         = execution.corr_array_n_cloud.execution_view();
  const auto corr_array_n_below         = execution.corr_array_n_below.execution_view();
  const auto X_mixt_comp_all_levs       = execution.X_mixt_comp_all_levs.execution_view();
  const auto lmin                       = execution.lmin.execution_view();
  const auto mixt_frac_max_mag          = execution.mixt_frac_max_mag.execution_view();
  const auto deep_soil_T_in_K           = execution.deep_soil_T_in_K.execution_view();
  const auto sfc_soil_T_in_K            = execution.sfc_soil_T_in_K.execution_view();
  const auto veg_T_in_K                 = execution.veg_T_in_K.execution_view();
  const auto wpthlp_sfc                 = execution.wpthlp_sfc.execution_view();
  const auto wprtp_sfc                  = execution.wprtp_sfc.execution_view();
  const auto upwp_sfc                   = execution.upwp_sfc.execution_view();
  const auto vpwp_sfc                   = execution.vpwp_sfc.execution_view();
  const auto upwp_sfc_pert              = execution.upwp_sfc_pert.execution_view();
  const auto vpwp_sfc_pert              = execution.vpwp_sfc_pert.execution_view();
  const auto um                         = execution.um.execution_view();
  const auto upwp                       = execution.upwp.execution_view();
  const auto vm                         = execution.vm.execution_view();
  const auto vpwp                       = execution.vpwp.execution_view();
  const auto up2                        = execution.up2.execution_view();
  const auto vp2                        = execution.vp2.execution_view();
  const auto up3                        = execution.up3.execution_view();
  const auto vp3                        = execution.vp3.execution_view();
  const auto wprtp                      = execution.wprtp.execution_view();
  const auto wpthlp                     = execution.wpthlp.execution_view();
  const auto rtp2                       = execution.rtp2.execution_view();
  const auto rtp3                       = execution.rtp3.execution_view();
  const auto thlp2                      = execution.thlp2.execution_view();
  const auto thlp3                      = execution.thlp3.execution_view();
  const auto rtpthlp                    = execution.rtpthlp.execution_view();
  const auto wp2                        = execution.wp2.execution_view();
  const auto thvm                       = execution.thvm.execution_view();
  const auto exner                      = execution.exner.execution_view();
  const auto rtm                        = execution.rtm.execution_view();
  const auto thlm                       = execution.thlm.execution_view();
  const auto rcm                        = execution.rcm.execution_view();
  const auto wp3                        = execution.wp3.execution_view();
  const auto p_in_Pa                    = execution.p_in_Pa.execution_view();
  const auto cloud_frac                 = execution.cloud_frac.execution_view();
  const auto wpthvp                     = execution.wpthvp.execution_view();
  const auto wp2thvp                    = execution.wp2thvp.execution_view();
  const auto wp2up                      = execution.wp2up.execution_view();
  const auto rtpthvp                    = execution.rtpthvp.execution_view();
  const auto thlpthvp                   = execution.thlpthvp.execution_view();
  const auto wp2rtp                     = execution.wp2rtp.execution_view();
  const auto wp2thlp                    = execution.wp2thlp.execution_view();
  const auto uprcp                      = execution.uprcp.execution_view();
  const auto vprcp                      = execution.vprcp.execution_view();
  const auto rc_coef_zm                 = execution.rc_coef_zm.execution_view();
  const auto wp4                        = execution.wp4.execution_view();
  const auto wpup2                      = execution.wpup2.execution_view();
  const auto wpvp2                      = execution.wpvp2.execution_view();
  const auto wp2up2                     = execution.wp2up2.execution_view();
  const auto wp2vp2                     = execution.wp2vp2.execution_view();
  const auto rho_ds_zt                  = execution.rho_ds_zt.execution_view();
  const auto thv_ds_zt                  = execution.thv_ds_zt.execution_view();
  const auto wm_zm                      = execution.wm_zm.execution_view();
  const auto wm_zt                      = execution.wm_zt.execution_view();
  const auto rho                        = execution.rho.execution_view();
  const auto rho_zm                     = execution.rho_zm.execution_view();
  const auto rho_ds_zm                  = execution.rho_ds_zm.execution_view();
  const auto invrs_rho_ds_zm            = execution.invrs_rho_ds_zm.execution_view();
  const auto invrs_rho_ds_zt            = execution.invrs_rho_ds_zt.execution_view();
  const auto thv_ds_zm                  = execution.thv_ds_zm.execution_view();
  const auto thlm_forcing               = execution.thlm_forcing.execution_view();
  const auto rtm_forcing                = execution.rtm_forcing.execution_view();
  const auto um_forcing                 = execution.um_forcing.execution_view();
  const auto vm_forcing                 = execution.vm_forcing.execution_view();
  const auto wprtp_forcing              = execution.wprtp_forcing.execution_view();
  const auto wpthlp_forcing             = execution.wpthlp_forcing.execution_view();
  const auto rtp2_forcing               = execution.rtp2_forcing.execution_view();
  const auto thlp2_forcing              = execution.thlp2_forcing.execution_view();
  const auto rtpthlp_forcing            = execution.rtpthlp_forcing.execution_view();
  const auto um_pert                    = execution.um_pert.execution_view();
  const auto vm_pert                    = execution.vm_pert.execution_view();
  const auto upwp_pert                  = execution.upwp_pert.execution_view();
  const auto vpwp_pert                  = execution.vpwp_pert.execution_view();
  const auto Skw_zm_smooth              = execution.Skw_zm_smooth.execution_view();
  const auto Lscale                     = execution.Lscale.execution_view();
  const auto thlprcp                    = execution.thlprcp.execution_view();
  const auto wprcp                      = execution.wprcp.execution_view();
  const auto w_up_in_cloud              = execution.w_up_in_cloud.execution_view();
  const auto w_down_in_cloud            = execution.w_down_in_cloud.execution_view();
  const auto cloudy_updraft_frac        = execution.cloudy_updraft_frac.execution_view();
  const auto cloudy_downdraft_frac      = execution.cloudy_downdraft_frac.execution_view();
  const auto ice_supersat_frac          = execution.ice_supersat_frac.execution_view();
  const auto rcm_in_layer               = execution.rcm_in_layer.execution_view();
  const auto cloud_cover                = execution.cloud_cover.execution_view();
  const auto invrs_tau_zm               = execution.invrs_tau_zm.execution_view();
  const auto ug                         = execution.ug.execution_view();
  const auto vg                         = execution.vg.execution_view();
  const auto rtm_ref                    = execution.rtm_ref.execution_view();
  const auto thlm_ref                   = execution.thlm_ref.execution_view();
  const auto um_ref                     = execution.um_ref.execution_view();
  const auto vm_ref                     = execution.vm_ref.execution_view();
  const auto rcm_mc                     = execution.rcm_mc.execution_view();
  const auto rvm_mc                     = execution.rvm_mc.execution_view();
  const auto thlm_mc                    = execution.thlm_mc.execution_view();
  const auto wprtp_mc                   = execution.wprtp_mc.execution_view();
  const auto wpthlp_mc                  = execution.wpthlp_mc.execution_view();
  const auto rtp2_mc                    = execution.rtp2_mc.execution_view();
  const auto thlp2_mc                   = execution.thlp2_mc.execution_view();
  const auto rtpthlp_mc                 = execution.rtpthlp_mc.execution_view();
  const auto Ncm_mc                     = execution.Ncm_mc.execution_view();
  const auto Ncm                        = execution.Ncm.execution_view();
  const auto Nccnm                      = execution.Nccnm.execution_view();
  const auto Kh_zt                      = execution.Kh_zt.execution_view();
  const auto Kh_zm                      = execution.Kh_zm.execution_view();
  const auto rfrzm                      = execution.rfrzm.execution_view();
  const auto Nc_in_cloud                = execution.Nc_in_cloud.execution_view();
  const auto wpNcp                      = execution.wpNcp.execution_view();
  const auto radht                      = execution.radht.execution_view();
  const auto wpsclrp_sfc                = execution.wpsclrp_sfc.execution_view();
  const auto wpedsclrp_sfc              = execution.wpedsclrp_sfc.execution_view();
  const auto sclrm                      = execution.sclrm.execution_view();
  const auto wpsclrp                    = execution.wpsclrp.execution_view();
  const auto sclrp2                     = execution.sclrp2.execution_view();
  const auto sclrp3                     = execution.sclrp3.execution_view();
  const auto sclrprtp                   = execution.sclrprtp.execution_view();
  const auto sclrpthlp                  = execution.sclrpthlp.execution_view();
  const auto sclrpthvp                  = execution.sclrpthvp.execution_view();
  const auto sclrm_forcing              = execution.sclrm_forcing.execution_view();
  const auto edsclrm                    = execution.edsclrm.execution_view();
  const auto edsclrm_forcing            = execution.edsclrm_forcing.execution_view();
  const auto lh_rt_clipped              = execution.lh_rt_clipped.execution_view();
  const auto lh_thl_clipped             = execution.lh_thl_clipped.execution_view();
  const auto lh_rc_clipped              = execution.lh_rc_clipped.execution_view();
  const auto lh_rv_clipped              = execution.lh_rv_clipped.execution_view();
  const auto lh_Nc_clipped              = execution.lh_Nc_clipped.execution_view();
  const auto lh_sample_point_weights    = execution.lh_sample_point_weights.execution_view();
  const auto K_hm                       = execution.K_hm.execution_view();
  const auto wp2hmp                     = execution.wp2hmp.execution_view();
  const auto rtphmp_zt                  = execution.rtphmp_zt.execution_view();
  const auto thlphmp_zt                 = execution.thlphmp_zt.execution_view();
  const auto hydromet                   = execution.hydromet.execution_view();
  const auto hydrometp2                 = execution.hydrometp2.execution_view();
  const auto wphydrometp                = execution.wphydrometp.execution_view();
  const auto hydromet_mc                = execution.hydromet_mc.execution_view();
  const auto hydromet_vel_zt            = execution.hydromet_vel_zt.execution_view();
  const auto hydromet_vel_covar_zt_impc = execution.hydromet_vel_covar_zt_impc.execution_view();
  const auto hydromet_vel_covar_zt_expc = execution.hydromet_vel_covar_zt_expc.execution_view();
  const auto mu_x_1_n                   = execution.mu_x_1_n.execution_view();
  const auto mu_x_2_n                   = execution.mu_x_2_n.execution_view();
  const auto sigma_x_1_n                = execution.sigma_x_1_n.execution_view();
  const auto sigma_x_2_n                = execution.sigma_x_2_n.execution_view();
  const auto corr_array_1_n             = execution.corr_array_1_n.execution_view();
  const auto corr_array_2_n             = execution.corr_array_2_n.execution_view();
  const auto corr_cholesky_mtx_1        = execution.corr_cholesky_mtx_1.execution_view();
  const auto corr_cholesky_mtx_2        = execution.corr_cholesky_mtx_2.execution_view();
  const auto X_nl_all_levs              = execution.X_nl_all_levs.execution_view();

  // Local variables.
  bool l_stats{};
  int itime{}; // Iteration counters
  int i{};
  int loop_iinit{};
  int loop_ifinal{};
  real_t stats_time{};

  // An absent native l_suppress_stats is equivalent to the C++ default false.
  l_stats = stats.enabled && !l_suppress_stats;
  timer_start("advance_clubb_to_end");

  // -------------------------------------------------------------------------------
  //                          Main Time Stepping Loop
  // -------------------------------------------------------------------------------
  loop_iinit  = iinit;
  loop_ifinal = ifinal;

  if (itime_start.has_value()) {
    loop_iinit = (*itime_start);
  }

  if (itime_end.has_value()) {
    loop_ifinal = (*itime_end);
  }

  for (itime = loop_iinit; itime < loop_ifinal + 1; itime += 1) {
    // Returned temporary fields belong to this timestep; nested routines keep
    // their own reservations until their queued users precede the next reuse.
    auto workspace_scope = workspace.scope();

    if (l_stats) {
      timer_start("stats_begin_timestep_api");
      stats_begin_timestep_api(itime, stats);
      timer_stop("stats_begin_timestep_api");
    }

    if (clubb_at_least_debug_level_api(2)) {
      timer_start("invalid_model_arrays");

      for (i = 0; i < ngrdcol; ++i) {
        // Check for NaN values in the model arrays
        if (invalid_model_arrays(
                workspace, gr.nzm, gr.nzt, hydromet_dim, hm_metadata.hydromet_list, sclr_dim,
                edsclr_dim, Kokkos::subview(um, i, Kokkos::ALL),
                Kokkos::subview(vm, i, Kokkos::ALL), Kokkos::subview(rtm, i, Kokkos::ALL),
                Kokkos::subview(wprtp, i, Kokkos::ALL), Kokkos::subview(thlm, i, Kokkos::ALL),
                Kokkos::subview(wpthlp, i, Kokkos::ALL), Kokkos::subview(rtp2, i, Kokkos::ALL),
                Kokkos::subview(thlp2, i, Kokkos::ALL), Kokkos::subview(rtpthlp, i, Kokkos::ALL),
                Kokkos::subview(wp2, i, Kokkos::ALL), Kokkos::subview(wp3, i, Kokkos::ALL),
                Kokkos::subview(wp2thvp, i, Kokkos::ALL), Kokkos::subview(wp2up, i, Kokkos::ALL),
                Kokkos::subview(rtpthvp, i, Kokkos::ALL), Kokkos::subview(thlpthvp, i, Kokkos::ALL),
                Kokkos::subview(hydromet, i, Kokkos::ALL, Kokkos::ALL),
                Kokkos::subview(sclrm, i, Kokkos::ALL, Kokkos::ALL),
                Kokkos::subview(edsclrm, i, Kokkos::ALL, Kokkos::ALL))) {
          set_host_element(err_info.err_code, i, clubb_fatal_error);
          std::cerr << text_at(err_info.err_header, i) << "\n";
          std::cerr << "Fatal error: a CLUBB variable is NaN in main time stepping loop."
                    << "\n";
          timer_stop("invalid_model_arrays");
          goto advance_finished;
        }
      }

      timer_stop("invalid_model_arrays"); // In
    }

    // Calculate thvm for use in prescribe_forcings.
    calculate_thvm(gr.nzt, ngrdcol, thlm, rtm, rcm, exner, thv_ds_zt, thvm);

    // Set large-scale tendencies and subsidence profiles
    timer_start("prescribe_forcings");

    prescribe_forcings(workspace, gr, gr.nzm, gr.nzt, ngrdcol,       // In
                       sclr_dim, edsclr_dim, sclr_idx,               // In
                       clubb_config_flags.l_do_expldiff_rtm_thlm,    // In
                       runtype, sfctype,                             // In
                       time_current, time_initial, dt_main,          // In
                       l_readiopdata,                                // In
                       um, vm, thlm,                                 // In
                       p_in_Pa, exner, rho, rho_zm, thvm,            // In
                       veg_T_in_K,                                   // In
                       l_modify_bc_for_cnvg_test,                    // In
                       clubb_config_flags.saturation_formula,        // In
                       stats,                                        // InOut
                       clubb_config_flags.l_add_dycore_grid,         // In
                       clubb_config_flags.grid_remap_method,         // In
                       gr.nzm, rho_ds_zm,                            // In
                       gr.zm,                                        // In
                       gr_dycore,                                    // In
                       rtm, wm_zm, wm_zt, ug, vg, um_ref, vm_ref,    // InOut
                       thlm_forcing, rtm_forcing, um_forcing,        // InOut
                       vm_forcing, wprtp_forcing, wpthlp_forcing,    // InOut
                       rtp2_forcing, thlp2_forcing, rtpthlp_forcing, // InOut
                       wpsclrp, sclrm_forcing, edsclrm_forcing,      // InOut
                       wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,    // InOut
                       T_sfc, p_sfc, sens_ht, latent_ht,             // InOut
                       wpsclrp_sfc, wpedsclrp_sfc, err_info,         // InOut
                       time_dependent_input, mpace);                 // In

    timer_stop("prescribe_forcings");

    if (clubb_at_least_debug_level_api(0)) {

      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        std::cerr << err_info.err_header_global << "\n";
        std::cerr << "Fatal error in prescribe_forcings:" << "\n";
        goto advance_finished;
      }
    }

    Kokkos::parallel_for(
        "clubb_driver::thermodynamic_forcing", (gr.nzt * ngrdcol), KOKKOS_LAMBDA(const int index) {
          const int i = (index % rfrzm.extent_int(0));
          const int k = (index / rfrzm.extent_int(0));
          // Compute total water in ice phase mixing ratio
          rfrzm(i, k) = zero;

          // Add microphysical tendencies to rtm_forcing
          rtm_forcing(i, k) = ((rtm_forcing(i, k) + rcm_mc(i, k)) + rvm_mc(i, k));

          // Add radiation and microphysical tendencies to thlm_forcing
          thlm_forcing(i, k) = ((thlm_forcing(i, k) + thlm_mc(i, k)) + radht(i, k));
        });

    Kokkos::parallel_for(
        "clubb_driver::momentum_forcing", (gr.nzm * ngrdcol), KOKKOS_LAMBDA(const int index) {
          const int i           = (index % wprtp_forcing.extent_int(0));
          const int k           = (index / wprtp_forcing.extent_int(0));
          wprtp_forcing(i, k)   = (wprtp_forcing(i, k) + wprtp_mc(i, k));
          wpthlp_forcing(i, k)  = (wpthlp_forcing(i, k) + wpthlp_mc(i, k));
          rtp2_forcing(i, k)    = (rtp2_forcing(i, k) + rtp2_mc(i, k));
          thlp2_forcing(i, k)   = (thlp2_forcing(i, k) + thlp2_mc(i, k));
          rtpthlp_forcing(i, k) = (rtpthlp_forcing(i, k) + rtpthlp_mc(i, k));
        });

    if (hydromet_dim > 0) {
      // Normalize local species slots; absent Native IDs remain negative.
      const int iiri = hm_metadata.iiri;
      const int iirs = hm_metadata.iirs;
      const int iirg = hm_metadata.iirg;
      Kokkos::parallel_for(
          "clubb_driver::ice_water", (gr.nzt * ngrdcol), KOKKOS_LAMBDA(const int index) {
            const int i = (index % rfrzm.extent_int(0));
            const int k = (index / rfrzm.extent_int(0));
            if (iiri >= 0) {
              rfrzm(i, k) = (rfrzm(i, k) + hydromet(i, k, iiri));
            }
            if (iirs >= 0) {
              rfrzm(i, k) = (rfrzm(i, k) + hydromet(i, k, iirs));
            }
            if (iirg >= 0) {
              rfrzm(i, k) = (rfrzm(i, k) + hydromet(i, k, iirg));
            }
          });
    }

    // Add effects of radiation on thlp2
    if (clubb_config_flags.l_calc_thlp2_rad) {
      timer_start("calculate_thlp2_rad_api");

      calculate_thlp2_rad_api(workspace, ngrdcol, gr.nzm, gr.nzt, gr, // In
                              rcm, thlprcp, radht, clubb_params,      // In
                              thlp2_forcing);                         // InOut

      timer_stop("calculate_thlp2_rad_api");
    }

    // Normal CLUBB run (with ascending grid)
    // Call the clubb core api for all columns
    timer_start("advance_clubb_core_api");

    advance_clubb_core_api(workspace, gr, gr.nzm, gr.nzt, ngrdcol,              // In
                           l_implemented, dt_main, fcor, fcor_y, sfc_elevation, // In
                           hydromet_dim,                                        // In
                           sclr_dim, sclr_tol, edsclr_dim, sclr_idx,            // In
                           thlm_forcing, rtm_forcing, um_forcing, vm_forcing,   // In
                           sclrm_forcing, edsclrm_forcing, wprtp_forcing,       // In
                           wpthlp_forcing, rtp2_forcing, thlp2_forcing,         // In
                           rtpthlp_forcing, wm_zm, wm_zt,                       // In
                           wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc,    // In
                           wpsclrp_sfc, wpedsclrp_sfc,                          // In
                           upwp_sfc_pert, vpwp_sfc_pert,                        // In
                           rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg,           // In
                           p_in_Pa,                                             // InOut
                           rho_zm, rho,                                         // In
                           exner,                                               // InOut
                           rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,               // In
                           invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt,               // In
                           hm_metadata.l_mix_rat_hm,                            // In
                           rfrzm, wphydrometp,                                  // In
                           wp2hmp, rtphmp_zt, thlphmp_zt,                       // In
                           dummy_dx, dummy_dy,                                  // In
                           clubb_params, nu_vert_res_dep, lmin,                 // In
                           mixt_frac_max_mag, T0, ts_nudge,                     // In
                           rtm_min, rtm_nudge_max_altitude,                     // In
                           clubb_config_flags,                                  // In
                           stats,                                               // InOut
                           um, vm, upwp, vpwp, up2, vp2, up3, vp3,              // InOut
                           thlm, rtm, wprtp, wpthlp,                            // InOut
                           wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp,         // InOut
                           sclrm, sclrp2, sclrp3, sclrprtp, sclrpthlp,          // InOut
                           wpsclrp, edsclrm, err_info,                          // InOut
                           rcm, cloud_frac,                                     // InOut
                           wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp,           // InOut
                           sclrpthvp,                                           // InOut
                           wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4,      // InOut
                           wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac,     // InOut
                           um_pert, vm_pert, upwp_pert, vpwp_pert,              // InOut
                           pdf_params, pdf_params_zm,                           // InOut
                           pdf_implicit_coefs_terms,                            // InOut
                           Kh_zm, Kh_zt,                                        // Out
                           thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud,      // Out
                           cloudy_updraft_frac, cloudy_downdraft_frac,          // Out
                           rcm_in_layer, cloud_cover, invrs_tau_zm,             // Out
                           Lscale);                                             // Out

    timer_stop("advance_clubb_core_api");

    if (clubb_at_least_debug_level_api(0)) {

      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        std::cerr << err_info.err_header_global << "\n";
        std::cerr << "Fatal error in advance_clubb_core, " << " "
                  << "check your parameter values and timestep" << "\n";
        goto advance_finished;
      }
    }

    // Radiation is always called on the first timestep in order to ensure
    // that the simulation is subject to radiative heating and cooling from
    // the first timestep.
    l_rad_itime = ((itime % static_cast<int>(std::floor((dt_rad / dt_main))) == 0) || (itime == 1));

    if (((microphys_scheme != "none") || (lh_microphys_type != lh_microphys_disabled)) ||
        l_silhs_rad) {
      // Call CLUBB's hydrometeor PDF subroutines for cases that use the PDF
      // for microphysics (either upscaled microphysics or 1).
      timer_start("pdf_hydromet_microphys_prep");

      pdf_hydromet_microphys_prep(workspace, gr,                          // In
                                  ngrdcol,                                // In
                                  pdf_dim,                                // In
                                  hydromet_dim,                           // In
                                  itime,                                  // In
                                  vert_decorr_coef,                       // In
                                  Nc_in_cloud,                            // In
                                  cloud_frac,                             // In
                                  ice_supersat_frac,                      // In
                                  rho_ds_zt,                              // In
                                  Lscale,                                 // In
                                  Kh_zm,                                  // In
                                  hydromet,                               // In
                                  wphydrometp,                            // In
                                  corr_array_n_cloud,                     // In
                                  corr_array_n_below,                     // In
                                  hm_metadata,                            // In
                                  pdf_params,                             // In
                                  clubb_params,                           // In
                                  clubb_config_flags,                     // In
                                  silhs_config_flags,                     // In
                                  l_rad_itime,                            // In
                                  stats,                                  // InOut
                                  err_info,                               // InOut
                                  hydrometp2,                             // Out
                                  mu_x_1_n,                               // Out
                                  mu_x_2_n,                               // Out
                                  sigma_x_1_n,                            // Out
                                  sigma_x_2_n,                            // Out
                                  corr_array_1_n,                         // Out
                                  corr_array_2_n,                         // Out
                                  corr_cholesky_mtx_1,                    // Out
                                  corr_cholesky_mtx_2,                    // Out
                                  precip_fracs,                           // InOut
                                  rtphmp_zt,                              // Out
                                  thlphmp_zt,                             // Out
                                  wp2hmp,                                 // Out
                                  X_nl_all_levs,                          // Out
                                  X_mixt_comp_all_levs,                   // Out
                                  lh_sample_point_weights,                // Out
                                  lh_rt_clipped,                          // Out
                                  lh_thl_clipped,                         // Out
                                  lh_rc_clipped,                          // Out
                                  lh_rv_clipped,                          // Out
                                  lh_Nc_clipped,                          // Out
                                  &hydromet_pdf_params,                   // Optional Out
                                  microphys, silhs, sampling_state, rng); // Case state

      timer_stop("pdf_hydromet_microphys_prep");
    }

    if (clubb_at_least_debug_level_api(0)) {

      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        std::cerr << err_info.err_header_global << "\n";
        std::cerr << "Fatal error after pdf_hydromet_microphys_prep" << "\n";
        goto advance_finished;
      }
    }

    // Compute Microphysics
    if (!(microphys_scheme == "none")) {
      timer_start("calc_microphys_scheme_tendcies");

      // Call microphysics scheme and produce microphysics tendencies.
      calc_microphys_scheme_tendcies(workspace, gr, ngrdcol, dt_main, time_current, // In
                                     pdf_dim, hydromet_dim, runtype,                // In
                                     thlm, p_in_Pa, exner, rho, rho_zm, rtm,        // In
                                     rcm, cloud_frac, wm_zt, wm_zm, wp2, wp3,       // In
                                     clubb_params,                                  // In
                                     hydromet, Nc_in_cloud,                         // In
                                     hm_metadata,                                   // In
                                     pdf_params, hydromet_pdf_params,               // In
                                     precip_fracs,                                  // In
                                     X_nl_all_levs, X_mixt_comp_all_levs,           // In
                                     lh_sample_point_weights,                       // In
                                     mu_x_1_n, mu_x_2_n,                            // In
                                     sigma_x_1_n, sigma_x_2_n,                      // In
                                     corr_array_1_n, corr_array_2_n,                // In
                                     lh_rt_clipped, lh_thl_clipped,                 // In
                                     lh_rc_clipped, lh_rv_clipped,                  // In
                                     lh_Nc_clipped,                                 // In
                                     silhs_config_flags.l_lh_importance_sampling,   // In
                                     silhs_config_flags.l_lh_instant_var_covar_src, // In
                                     clubb_config_flags.saturation_formula,         // In
                                     stats,                                         // InOut
                                     Nccnm,                                         // InOut
                                     hydromet_mc, Ncm_mc, rcm_mc, rvm_mc,           // Out
                                     thlm_mc, hydromet_vel_zt,                      // Out
                                     hydromet_vel_covar_zt_impc,                    // Out
                                     hydromet_vel_covar_zt_expc,                    // Out
                                     wprtp_mc, wpthlp_mc, rtp2_mc,                  // Out
                                     thlp2_mc, rtpthlp_mc,                          // Out
                                     Skw_zm_smooth,                                 // Out
                                     microphys, kk, morrison);                      // Case state

      timer_stop("calc_microphys_scheme_tendcies");
      timer_start("advance_microphys");

      // Advance predictive microphysics fields one model timestep.
      advance_microphys(workspace, gr, ngrdcol, dt_main, time_current, // In
                        hydromet_dim, hm_metadata,                     // In
                        wm_zt, wp2, exner, rho, rho_zm, rcm,           // In
                        cloud_frac, Kh_zm, Skw_zm_smooth,              // In
                        rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,         // In
                        hydromet_mc, Ncm_mc, Lscale,                   // In
                        hydromet_vel_covar_zt_impc,                    // In
                        hydromet_vel_covar_zt_expc,                    // In
                        clubb_params, nu_vert_res_dep,                 // In
                        clubb_config_flags.tridiag_solve_method,       // In
                        clubb_config_flags.fill_holes_type,            // In
                        clubb_config_flags.l_upwind_xm_ma,             // In
                        stats, hydromet, hydromet_vel_zt, hydrometp2,  // InOut
                        K_hm, Ncm, Nc_in_cloud, rvm_mc, thlm_mc,       // InOut
                        err_info,                                      // InOut
                        wphydrometp, wpNcp,                            // Out
                        microphys); // Case state                          // Out

      timer_stop("advance_microphys");

      if (clubb_at_least_debug_level_api(0)) {

        if (any_equal(err_info.err_code, clubb_fatal_error)) {

          for (i = 0; i < ngrdcol; ++i) {

            if (host_element(err_info.err_code, i) == clubb_fatal_error) {
              std::cerr << text_at(err_info.err_header, i) << "\n";
            }
          }

          std::cerr << "Fatal error in advance_microphys:" << "\n";
          goto advance_finished;
        }
      }
    } else {
      // Cloud water sedimentation.
      if (l_cloud_sed) {
        timer_start("cloud_drop_sed");

        Kokkos::parallel_for(
            "clubb_driver::cloud_number", Ncm.size(), KOKKOS_LAMBDA(const int index) {
              const int i = (index % Ncm.extent_int(0));
              const int k = (index / Ncm.extent_int(0));
              Ncm(i, k)   = (Nc_in_cloud(i, k) * cloud_frac(i, k));
            });

        Kokkos::deep_copy(rcm_mc, zero);

        Kokkos::deep_copy(thlm_mc, zero);

        // Note:  it would be very easy to upscale the cloud water sedimentation
        // flux, so we should look into adding an upscaled option.
        cloud_drop_sed(workspace, gr, ngrdcol, rcm, Ncm, // In
                       rho_zm, rho, exner, sigma_g,      // In
                       stats, rcm_mc, thlm_mc);          // InOut

        timer_stop("cloud_drop_sed");
      }

      if (stats.l_sample) {

        Kokkos::parallel_for(
            "clubb_driver::cloud_number", Ncm.size(), KOKKOS_LAMBDA(const int index) {
              const int i = (index % Ncm.extent_int(0));
              const int k = (index / Ncm.extent_int(0));
              Ncm(i, k)   = (Nc_in_cloud(i, k) * cloud_frac(i, k));
            });

        stats_update("Ncm", Ncm, stats);
        stats_update("Nc_in_cloud", Nc_in_cloud, stats);
      }
    }

    if (stats.l_sample) {
      // Update stats variables for total microphysical tendencies and
      // associated second-order moments.
      stats_update("rvm_mc", rvm_mc, stats);
      stats_update("rcm_mc", rcm_mc, stats);
      Real2D stats_rtm_mc_1 =
          workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "stats_rtm_mc_1",
                                       rvm_mc.extent_int(0), rvm_mc.extent_int(1));

      Kokkos::parallel_for(
          "clubb_driver::stats_rtm_mc", rvm_mc.size(), KOKKOS_LAMBDA(const int index) {
            const int i          = (index % rvm_mc.extent_int(0));
            const int k          = (index / rvm_mc.extent_int(0));
            stats_rtm_mc_1(i, k) = (rvm_mc(i, k) + rcm_mc(i, k));
          });

      stats_update("rtm_mc", stats_rtm_mc_1, stats);
      stats_update("thlm_mc", thlm_mc, stats);
      stats_update("wprtp_mc", wprtp_mc, stats);
      stats_update("wpthlp_mc", wpthlp_mc, stats);
      stats_update("rtp2_mc", rtp2_mc, stats);
      stats_update("thlp2_mc", thlp2_mc, stats);
      stats_update("rtpthlp_mc", rtpthlp_mc, stats);
    }

    if (rad_scheme != "none") {
      timer_start("advance_clubb_radiation");

      advance_clubb_radiation(workspace, gr, ngrdcol, hydromet_dim, pdf_dim, lh_num_samples, // In
                              l_rad_itime, dt_main, day, month, year,                        // In
                              lat_vals, lon_vals, time_current, time_initial,                // In
                              rho, rho_zm, p_in_Pa, exner,                                   // In
                              wpthlp_sfc, wprtp_sfc, p_sfc,                                  // In
                              cloud_frac, ice_supersat_frac,                                 // In
                              thlm, rtm, rcm,                                                // In
                              X_nl_all_levs, lh_rt_clipped, lh_thl_clipped,
                              lh_rc_clipped,                                  // In
                              lh_sample_point_weights, hydromet, hm_metadata, // In
                              stats, err_info,                                // InOut
                              deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K,  // InOut
                              radht,                                          // InOut
                              parameters_radiation,                           // In
                              radiation_variables,                            // InOut
                              simple_radiation, soil_vegetation);             // In

      timer_stop("advance_clubb_radiation");

      if (clubb_at_least_debug_level_api(0)) {

        if (any_equal(err_info.err_code, clubb_fatal_error)) {
          goto advance_finished;
        }
      }
    }

    if ((clubb_config_flags.grid_adapt_in_time_method > no_grid_adaptation) ||
        (clubb_config_flags.l_add_dycore_grid && stats.l_sample)) {
      // Update stats with the current grid information needed to correctly
      // map the output to the netcdf file
      timer_start("stats_update_grid");
      stats_update_grid(gr.zt, gr.zm, rho_ds_zm, gr.zm, p_sfc, stats);
      timer_stop("stats_update_grid");
    }

    if (l_stats) {
      stats_time = static_cast<real_t>(
          (time_initial + (static_cast<double>(itime) * static_cast<double>(dt_main))));

      // End statistics timestep and flush sampled buffers to file.
      // Flush sub-timestep averages every timestep and write at output boundaries.
      timer_start("stats_end_timestep_api");
      stats_end_timestep_api(stats_time, stats, err_info);
      timer_stop("stats_end_timestep_api");

      if (clubb_at_least_debug_level_api(0)) {

        if (any_equal(err_info.err_code, clubb_fatal_error)) {
          std::cerr << err_info.err_header_global << "\n";
          std::cerr << "Fatal error calling stats_end_timestep_api in run_clubb"
                    << "\n";
          goto advance_finished;
        }
      }
    }

    // Set Time
    // Advance time here, not in advance_clubb_core,
    // in order to facilitate use of stats.
    // A host model, e.g. WRF, would advance time outside
    // of advance_clubb_core.  Vince Larson 7 Feb 2006
    time_current = (time_initial + (static_cast<double>(itime) * static_cast<double>(dt_main)));

    // This was moved from above to be less confusing to the user,
    // since before it would appear as though the last timestep
    // was not executed. -dschanen 19 May 08
    if (l_stdout && (clubb_at_least_debug_level_api(0))) {
      std::cout << "iteration: " << " " << itime << " " << " / " << " " << ifinal << " "
                << " -- time = " << " " << time_current << " "
                << "s / " << " " << time_final << " " << "s" << "\n";
    }
  }

// -------------------------------------------------------------------------------
//                        End Main Time Stepping Loop
// -------------------------------------------------------------------------------
advance_finished:;

  // End Main Time Stepping Loop
  timer_stop("advance_clubb_to_end");
}

void ClubbDriver::clean_up_clubb(err_info_type& err_info) { // InOut

  // Close case resources and release the host and execution View allocations.
  auto& l_t_dependent     = time_dependent_input.l_t_dependent;
  auto& rad_scheme        = parameters_radiation.rad_scheme;
  auto& lh_microphys_type = microphys.lh_microphys_type;

  // Local variables.

  // --------------------------------- Begin Code ---------------------------------
  timer_start("clean_up_clubb");

  // Complete case kernels before releasing device and host allocations.
  Kokkos::fence("clubb_driver::cleanup");
  workspace.clear();
  execution                    = {};
  execution_fields_initialized = false;

  // Close stats files
  timer_start("stats_finalize_api");
  stats_finalize_api(stats, err_info);
  timer_stop("stats_finalize_api");

  if (l_t_dependent) {
    finalize_t_dependent_input(time_dependent_input); // InOut
  }

  finalize_extended_atm(extended_atmosphere); // InOut
  // De-allocate the arrays for the grid
  cleanup_grid_api(gr);

  if (rad_scheme != "none") {
    cleanup_radiation_variables(radiation_variables); // InOut
  }

  cleanup_microphys(microphys);
  if (l_input_fields || l_restart) {
    cleanup_input_fields(inputfields); // InOut
  }

  if (lh_microphys_type != lh_microphys_disabled) {
    cleanup_latin_hypercube_arrays(sampling_state);
  }

  lmin                  = {};
  mixt_frac_max_mag     = {};
  um                    = {};
  vm                    = {};
  upwp                  = {};
  vpwp                  = {};
  up2                   = {};
  up3                   = {};
  vp2                   = {};
  vp3                   = {};
  thlm                  = {};
  rtm                   = {};
  wprtp                 = {};
  wpthlp                = {};
  wprcp                 = {};
  w_up_in_cloud         = {};
  w_down_in_cloud       = {};
  cloudy_updraft_frac   = {};
  cloudy_downdraft_frac = {};
  wp2                   = {};
  wp3                   = {};
  rtp2                  = {};
  thlp2                 = {};
  rtpthlp               = {};
  p_in_Pa               = {};
  exner                 = {};
  rho                   = {};
  rho_zm                = {};
  rho_ds_zm             = {};
  rho_ds_zt             = {};
  invrs_rho_ds_zm       = {};
  invrs_rho_ds_zt       = {};
  thv_ds_zm             = {};
  thv_ds_zt             = {};
  thlm_forcing          = {};
  rtm_forcing           = {};
  um_forcing            = {};
  vm_forcing            = {};
  wprtp_forcing         = {};
  wpthlp_forcing        = {};
  rtp2_forcing          = {};
  thlp2_forcing         = {};
  rtpthlp_forcing       = {};

  //  Variables used to track perturbed version of winds.
  um_pert   = {};
  vm_pert   = {};
  upwp_pert = {};
  vpwp_pert = {};

  //  Imposed large scale w
  wm_zm = {};
  wm_zt = {};

  //  Cloud water variables
  rcm               = {};
  delta_zm          = {};
  cloud_frac        = {};
  ice_supersat_frac = {};
  rcm_in_layer      = {};
  cloud_cover       = {};
  invrs_tau_zm      = {};

  //  Passive scalar variables
  //  Note that sclr_dim can be 0
  wpsclrp_sfc     = {};
  sclrm           = {};
  sclrp2          = {};
  sclrp3          = {};
  sclrm_forcing   = {};
  sclrprtp        = {};
  sclrpthlp       = {};
  wpedsclrp_sfc   = {};
  edsclrm_forcing = {};
  edsclrm         = {};
  wpsclrp         = {};
  sigma_sqd_w     = {};
  sigma_sqd_w_zt  = {};
  wp2_zt          = {};
  Skw_zm_smooth   = {};
  ug              = {};
  vg              = {};
  um_ref          = {};
  vm_ref          = {};
  thlm_ref        = {};
  rtm_ref         = {};
  thvm            = {};
  radht           = {};
  thlprcp         = {};
  thlp3           = {};
  rtp3            = {};

  //  Buoyancy related moments
  rtpthvp     = {};
  thlpthvp    = {};
  wpthvp      = {};
  wp2thvp     = {};
  wp2up       = {};
  wp2rtp      = {};
  wp2thlp     = {};
  uprcp       = {};
  vprcp       = {};
  rc_coef_zm  = {};
  wp4         = {};
  wpup2       = {};
  wpvp2       = {};
  wp2up2      = {};
  wp2vp2      = {};
  Kh_zt       = {};
  Kh_zm       = {};
  K_hm        = {};
  em          = {};
  Lscale      = {};
  tau_zm      = {};
  tau_zt      = {};
  Nccnm       = {};
  hydromet    = {};
  hydrometp2  = {};
  wphydrometp = {};
  Ncm         = {};
  wpNcp       = {};

  //  High-order passive scalars
  sclrpthvp                  = {};
  thlm_mc                    = {};
  rvm_mc                     = {};
  rcm_mc                     = {};
  wprtp_mc                   = {};
  wpthlp_mc                  = {};
  rtp2_mc                    = {};
  thlp2_mc                   = {};
  rtpthlp_mc                 = {};
  hydromet_mc                = {};
  Ncm_mc                     = {};
  hydromet_vel_zt            = {};
  hydromet_vel_covar_zt_impc = {};
  hydromet_vel_covar_zt_expc = {};
  X_nl_all_levs              = {};
  lh_rt_clipped              = {};
  lh_thl_clipped             = {};
  lh_rc_clipped              = {};
  lh_rv_clipped              = {};
  lh_Nc_clipped              = {};
  X_mixt_comp_all_levs       = {};
  lh_sample_point_weights    = {};
  Nc_in_cloud                = {};
  sclr_tol                   = {};
  p_sfc                      = {};
  T_sfc                      = {};
  fcor                       = {};
  fcor_y                     = {};
  dummy_dx                   = {};
  dummy_dy                   = {};
  deep_soil_T_in_K           = {};
  sfc_soil_T_in_K            = {};
  veg_T_in_K                 = {};
  upwp_sfc_pert              = {};
  vpwp_sfc_pert              = {};
  sfc_elevation              = {};
  zm_init                    = {};
  zm_top                     = {};
  deltaz                     = {};
  rrm                        = {};
  Nrm                        = {};
  wp2hmp                     = {};
  rtphmp_zt                  = {};
  thlphmp_zt                 = {};
  hydromet_pdf_params        = {};
  corr_array_1_n             = {};
  corr_array_2_n             = {};
  corr_cholesky_mtx_1        = {};
  corr_cholesky_mtx_2        = {};
  mu_x_1_n                   = {};
  mu_x_2_n                   = {};
  sigma_x_1_n                = {};
  sigma_x_2_n                = {};
  rfrzm                      = {};
  wpthlp_sfc                 = {};
  wprtp_sfc                  = {};
  upwp_sfc                   = {};
  vpwp_sfc                   = {};
  thlm_init                  = {};
  rtm_init                   = {};
  um_init                    = {};
  vm_init                    = {};
  ug_init                    = {};
  vg_init                    = {};
  wp2_init                   = {};
  up2_init                   = {};
  vp2_init                   = {};
  upwp_init                  = {};
  rcm_init                   = {};
  wm_zt_init                 = {};
  wm_zm_init                 = {};
  em_init                    = {};
  exner_init                 = {};
  thvm_init                  = {};
  p_in_Pa_init               = {};
  p_in_Pa_zm_init            = {};
  rho_init                   = {};
  rho_zm_init                = {};
  rho_ds_zm_init             = {};
  rho_ds_zt_init             = {};
  invrs_rho_ds_zm_init       = {};
  invrs_rho_ds_zt_init       = {};
  thv_ds_zm_init             = {};
  thv_ds_zt_init             = {};
  rtm_ref_init               = {};
  thlm_ref_init              = {};
  um_ref_init                = {};
  vm_ref_init                = {};
  Ncm_init                   = {};
  Nc_in_cloud_init           = {};
  Nccnm_init                 = {};
  sclrm_init                 = {};
  edsclrm_init               = {};

  if (clubb_params.is_allocated()) {
    clubb_params = {};
  }

  if (clubb_params_all.is_allocated()) {
    clubb_params_all = {};
  }

  ngrdcol          = 1;
  total_param_sets = 1;
  batch_size       = 1;
  timer_stop("clean_up_clubb");
  timer_finalize((output_file_prefix + ".timing"));
}

void ClubbDriver::initialize_clubb(
    const grid& gr, int ngrdcol, int iunit, const std::string& forcings_file_path, // In
    real_t& p_sfc,                                                                 // InOut
    DriverSlice1 zm_init,                                                          // In
    int sclr_dim, int edsclr_dim, const sclr_idx_type& sclr_idx,                   // In
    const clubb_config_flags_type& clubb_config_flags,                             // In
    bool l_modify_ic_with_cubic_int,                                               // In
    bool l_ascending_grid, DriverSlice1 fcor_y, bool l_readiopdata,                // In
    DriverSlice2 thlm, DriverSlice2 rtm, DriverSlice2 um, DriverSlice2 vm,         // InOut
    DriverSlice2 ug, DriverSlice2 vg,                                              // InOut
    DriverSlice2 wp2, DriverSlice2 up2, DriverSlice2 vp2, DriverSlice2 upwp,       // InOut
    DriverSlice2 rcm,                                                              // InOut
    DriverSlice2 wm_zt, DriverSlice2 wm_zm, DriverSlice2 em, DriverSlice2 exner,   // InOut
    DriverSlice2 thvm, DriverSlice2 p_in_Pa,                                       // InOut
    DriverSlice1 p_in_Pa_zm_init,                                                  // Out
    DriverSlice2 rho, DriverSlice2 rho_zm, DriverSlice2 rho_ds_zm,                 // InOut
    DriverSlice2 rho_ds_zt,                                                        // InOut
    DriverSlice2 invrs_rho_ds_zm, DriverSlice2 invrs_rho_ds_zt,                    // InOut
    DriverSlice2 thv_ds_zm, DriverSlice2 thv_ds_zt,                                // InOut
    DriverSlice2 rtm_ref, DriverSlice2 thlm_ref,                                   // InOut
    DriverSlice2 um_ref, DriverSlice2 vm_ref,                                      // InOut
    DriverSlice2 Ncm, DriverSlice2 Nc_in_cloud, DriverSlice2 Nccnm,                // InOut
    err_info_type& err_info,                                                       // InOut
    DriverSlice1 deep_soil_T_in_K, DriverSlice1 sfc_soil_T_in_K,                   // InOut
    DriverSlice1 veg_T_in_K,                                                       // InOut
    DriverSlice3 sclrm, DriverSlice3 edsclrm) {                                    // Out

  // Read sounding data and prepare the pressure, cloud and turbulence initial profiles.
  auto workspace_scope   = workspace.scope();
  auto& l_t_dependent    = time_dependent_input.l_t_dependent;
  auto& microphys_scheme = microphys.microphys_scheme;
  auto& Nc0_in_cloud     = microphys.Nc0_in_cloud;

  // Local variables.
  DriverArray2 p_in_Pa_zm = workspace_scope.view<DriverArray2>(
      "p_in_Pa_zm", ngrdcol, gr.nzm); // Pressure on momentum levels  [Pa]
  DriverArray1 rtm_sfc = workspace_scope.view<DriverArray1>(
      "rtm_sfc", ngrdcol); // Surface total water mixing ratio       [kg/kg]
  DriverArray1 thlm_sfc = workspace_scope.view<DriverArray1>(
      "thlm_sfc", ngrdcol);       // Surface liq. water potential temp.     [K]
  real_t cloud_top_height{};      // Cloud top altitude in initial profile  [m]
  real_t domain_depth{};          // Domain depth                           [m]
  real_t em_max{};                // Maximum value of initial subgrid TKE   [m^2/s^2]
  std::string temperature_type{}; // Type of temperature sounding
  std::string altitude_type{};    // Type of altitude sounding
  std::string subs_type{};        // Type of large-scale subsidence sounding
  int i{};
  int k{};
  grid gr_dycore{};
  DriverArray2 p_in_Pa_dycore{}; // p_in_Pa interpolated to the dycore grid
  auto gr_zm_host        = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, gr.zm);
  auto gr_zt_host        = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, gr.zt);
  auto gr_dycore_zt_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, gr_dycore.zt);

  // Read sounding information
  // Only use first value of p_sfc and zm_init because sounding files are not configured
  // to use multiple columns yet
  read_sounding(workspace, ngrdcol,                          // In
                sclr_dim, edsclr_dim, sclr_idx,              // InOut
                gr, iunit, runtype, zm_init(0),              // InOut
                clubb_config_flags.saturation_formula,       // InOut
                l_readiopdata, l_modify_ic_with_cubic_int,   // InOut
                p_sfc,                                       // InOut
                thlm, temperature_type, rtm, um, vm, ug, vg, // InOut
                altitude_type, p_in_Pa, subs_type, wm_zt,    // InOut
                rtm_sfc, thlm_sfc, sclrm,                    // InOut
                edsclrm,                                     // Out
                parameters_radiation,                        // In
                extended_atmosphere, iop);                   // InOut

  // Covert sounding input to CLUBB compatible input
  initialize_clubb_variables(ngrdcol, sclr_dim, edsclr_dim, sclr_idx,           // In
                             gr, altitude_type, temperature_type,               // In
                             host_values({p_sfc}), rtm_sfc, rtm, l_readiopdata, // In
                             clubb_config_flags.saturation_formula,             // In
                             thlm, p_in_Pa,                                     // InOut
                             p_in_Pa_zm,                                        // Out
                             exner, rho, rho_zm,                                // Out
                             rcm, thvm, rho_ds_zm,                              // Out
                             rho_ds_zt, invrs_rho_ds_zm,                        // Out
                             invrs_rho_ds_zt, thv_ds_zm,                        // Out
                             thv_ds_zt,                                         // Out
                             sclrm, edsclrm);                                   // InOut

  for (int k = 0; k < p_in_Pa_zm_init.extent_int(0); ++k) {
    p_in_Pa_zm_init(k) = p_in_Pa_zm(0, k);
  }

  // Determine initial value cloud droplet number concentration when Nc
  // is predicted.
  for (int k = 0; k < Nc_in_cloud.extent_int(1); ++k) {

    for (int icol = 0; icol < Nc_in_cloud.extent_int(0); ++icol) {
      Nc_in_cloud(icol, k) = (Nc0_in_cloud / rho(icol, k));
    }
  }

  for (i = 0; i < ngrdcol; ++i) {

    for (k = 0; k < gr.nzt; ++k) {

      if (rcm(i, k) > zero) {
        // The initial profile at this level is entirely saturated (due to
        // constant moisture and temperature over the level).  The level is
        // entirely cloudy.
        Ncm(i, k) = Nc_in_cloud(i, k);
      } else {
        // The initial profile at this level is entirely unsaturated (due to
        // constant moisture and temperature over the level).  The level is
        // entirely clear.
        Ncm(i, k) = (Nc_in_cloud(i, k) * cloud_frac_min);
      }
    }
  }

  if (microphys_scheme == "coamps") {
    // Initialize Nccnm as in COAMPS-LES
    for (int k = 0; k < gr.nzt; ++k) {

      for (int icol = 0; icol < Nccnm.extent_int(0); ++icol) {
        Nccnm(icol, k) =
            (((real_t(30.0) * (one + std::exp(((-gr_zt_host(icol, k)) / real_t(2000.0))))) *
              cm3_per_m3) /
             rho(icol, k));
      }
    }
  } else {

    for (int k = 0; k < Nccnm.extent_int(1); ++k) {

      for (int icol = 0; icol < Nccnm.extent_int(0); ++icol) {
        Nccnm(icol, k) = zero;
      }
    }
  }

  if (subs_type == wm_name) {
    zt2zm_api(gr.nzm, gr.nzt, ngrdcol, gr, wm_zt, wm_zm);

    for (int icol = 0; icol < wm_zm.extent_int(0); ++icol) {
      wm_zm(icol, 0) = real_t(0.0);
    }

    for (int icol = 0; icol < wm_zm.extent_int(0); ++icol) {
      wm_zm(icol, (gr.nzm - 1)) = real_t(0.0);
    }
  } else if (subs_type == omega_name) {

    for (k = 0; k < gr.nzt; ++k) {

      for (int icol = 0; icol < wm_zt.extent_int(0); ++icol) {
        wm_zt(icol, k) = ((-wm_zt(icol, k)) / (grav * rho(icol, k)));
      }
    }

    for (int icol = 0; icol < wm_zt.extent_int(0); ++icol) {
      wm_zt(icol, (gr.nzt - 1)) = real_t(0.0);
    }

    zt2zm_api(gr.nzm, gr.nzt, ngrdcol, gr, wm_zt, wm_zm);

    for (int icol = 0; icol < wm_zm.extent_int(0); ++icol) {
      wm_zm(icol, (gr.nzm - 1)) = real_t(0.0);
    }
  } else {

    for (int k = 0; k < wm_zt.extent_int(1); ++k) {

      for (int icol = 0; icol < wm_zt.extent_int(0); ++icol) {
        wm_zt(icol, k) = real_t(0.0);
      }
    }

    for (int k = 0; k < wm_zm.extent_int(1); ++k) {

      for (int icol = 0; icol < wm_zm.extent_int(0); ++icol) {
        wm_zm(icol, k) = real_t(0.0);
      }
    }
  }

  // Initialize Time Dependent Input
  if (l_t_dependent) {

    if (clubb_config_flags.l_add_dycore_grid) {
      //  Initialize the dycore grid
      setup_gr_dycore(workspace, iunit, ngrdcol,                              // In
                      Kokkos::subview(gr_zm_host, Kokkos::ALL, 0),            // In
                      Kokkos::subview(gr_zm_host, Kokkos::ALL, (gr.nzm - 1)), // In
                      l_ascending_grid,                                       // In
                      gr_dycore,                                              // Out
                      err_info);                                              // InOut

      gr_dycore_zt_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, gr_dycore.zt);
      p_in_Pa_dycore = workspace_scope.view<DriverArray2>("p_in_Pa_dycore", ngrdcol, gr_dycore.nzt);

      for (i = 0; i < ngrdcol; ++i) {
        lin_interp_between_grids(workspace, gr_dycore.nzt, gr.nzm,
                                 Kokkos::subview(gr_dycore_zt_host, i, Kokkos::ALL),
                                 Kokkos::subview(gr_zm_host, i, Kokkos::ALL),
                                 Kokkos::subview(p_in_Pa_zm, i, Kokkos::ALL),
                                 Kokkos::subview(p_in_Pa_dycore, i, Kokkos::ALL));
      }

      // if we want to simulate the forcings from the host model on the dycore grid,
      // we first read in all the forcings from the file to the dycore grid and every
      // timestep remap the results to the current physics grid, so we need to pay
      // attention on how the forcings are read in from object t_dependent_forcing_data
      // in time_dependent_input
      initialize_t_dependent_input(workspace,
                                   iunit,                                              // In
                                   runtype,                                            // In
                                   l_readiopdata,                                      // In
                                   rho_zm,                                             // In
                                   gr_dycore.nzt,                                      // In
                                   Kokkos::subview(gr_dycore_zt_host, 0, Kokkos::ALL), // In
                                   Kokkos::subview(p_in_Pa_dycore, 0, Kokkos::ALL),    // In
                                   clubb_config_flags.l_add_dycore_grid,               // In
                                   clubb_config_flags.grid_adapt_in_time_method,       // In
                                   time_dependent_input, iop);                         // InOut

      p_in_Pa_dycore = {};
    } else {
      // no simulating forcings input from host model on dycore grid
      // l_add_dycore_grid=.false.
      initialize_t_dependent_input(workspace, iunit,                             // In
                                   runtype,                                      // In
                                   l_readiopdata,                                // In
                                   rho,                                          // In
                                   gr.nzt,                                       // In
                                   Kokkos::subview(gr_zt_host, 0, Kokkos::ALL),  // In
                                   Kokkos::subview(p_in_Pa, 0, Kokkos::ALL),     // In
                                   clubb_config_flags.l_add_dycore_grid,         // In
                                   clubb_config_flags.grid_adapt_in_time_method, // In
                                   time_dependent_input, iop);                   // InOut
    }
  } else {

    if (clubb_config_flags.l_add_dycore_grid && (runtype != "atex")) {

      if (runtype != "gabls2") {
        std::cerr << "WARNING! The option l_add_dycore_grid=.true. or " << " "
                  << "grid_adapt_in_time_method>0 can currently only " << " "
                  << "be used for cases with forcings from an input file and for the atex "
                  << " "
                  << "case, so for this case l_add_dycore_grid must be set to .false. "
                  << " " << "and grid_adapt_in_time_method=0" << "\n";

        //  General error -> set all entries to clubb_fatal_error
        Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
        std::cerr << "Flags have to be changed like mentioned." << "\n";
        return;
      }
    }
  }

  // Initialize TKE and other fields as needed.
  if (runtype == "generic") {
    // Generic case

    for (int k = 0; k < em.extent_int(1); ++k) {

      for (int icol = 0; icol < em.extent_int(0); ++icol) {
        em(icol, k) = real_t(1.0);
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }
  } else if (runtype == "bomex" || runtype == "ekman" || runtype == "atex_long") {
    // GCSS BOMEX
    //          do k = 1, gr%nz
    //             em(k) = 1.0_core_rknd - (gr%zm(1,k)/3000.0_core_rknd)
    //             if ( em(k) < em_min ) then
    //                em(k) = em_min
    //             end if
    //          end do
    //          em(1) = em(2)
    //          em(gr%nz) = em(gr%nz-1)
    for (int k = 0; k < em.extent_int(1); ++k) {

      for (int icol = 0; icol < em.extent_int(0); ++icol) {
        em(icol, k) = em_min;
      }
    }
  } else if (runtype == "arm") {
    // GCSS ARM
    // Reduction of initial sounding for stability (inactive native alternative).
    // do k = 1, gr%nz
    //   if ( gr%zm(1,k) < 150.0_core_rknd ) then
    //     em(k) = ( 0.15_core_rknd * (1.0_core_rknd - gr%zm(1,k)/150.0_core_rknd) )
    //             / rho_zm(k)
    //   else
    //     em(k) = em_min
    //   end if
    // end do
    // em(1) = em(2)
    // em(gr%nz) = em(gr%nz-1)
    // End reduction of initial sounding for stability, 24 Jan 2007.
    for (int k = 0; k < em.extent_int(1); ++k) {

      for (int icol = 0; icol < em.extent_int(0); ++icol) {
        em(icol, k) = em_min;
      }
    }
  } else if (runtype == "arm_97") {
    // June 27 1997 ARM case

    for (int k = 0; k < em.extent_int(1); ++k) {

      for (int icol = 0; icol < em.extent_int(0); ++icol) {
        em(icol, k) = real_t(1.0);
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }
  } else if (runtype == "twp_ice") {
    // twp_ice

    for (int k = 0; k < em.extent_int(1); ++k) {

      for (int icol = 0; icol < em.extent_int(0); ++icol) {
        em(icol, k) = real_t(1.0);
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }
  } else if (runtype == "arm_0003") {
    // March 2000 ARM case

    for (int k = 0; k < em.extent_int(1); ++k) {

      for (int icol = 0; icol < em.extent_int(0); ++icol) {
        em(icol, k) = real_t(1.0);
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }
  } else if (runtype == "arm_3year") {
    // 3 year ARM case

    for (int k = 0; k < em.extent_int(1); ++k) {

      for (int icol = 0; icol < em.extent_int(0); ++icol) {
        em(icol, k) = real_t(1.0);
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }
  } else if (runtype == "cloud_feedback_s6" || runtype == "cloud_feedback_s6_p2k" ||
             runtype == "cloud_feedback_s11" || runtype == "cloud_feedback_s11_p2k" ||
             runtype == "cloud_feedback_s12" || runtype == "cloud_feedback_s12_p2k") {

    for (int k = 0; k < em.extent_int(1); ++k) {

      for (int icol = 0; icol < em.extent_int(0); ++icol) {
        em(icol, k) = real_t(1.0);
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }
  } else if (runtype == "astex_a209") {
    // ASTEX_A209 case 16 Jul, 2010 kcwhite
    cloud_top_height = real_t(700.);

    em_max = real_t(1.0);

    for (i = 0; i < ngrdcol; ++i) {

      for (k = 0; k < gr.nzm; ++k) {

        if (gr_zm_host(i, k) < cloud_top_height) {
          em(i, k) = em_max;
        } else {
          em(i, k) = em_min;
        }
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, 0) = em(icol, 1);
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }
  } else if (runtype == "fire") {
    // GCSS FIRE Sc
    cloud_top_height = real_t(700.);

    for (i = 0; i < ngrdcol; ++i) {

      for (k = 0; k < gr.nzm; ++k) {

        if (gr_zm_host(i, k) < cloud_top_height) {
          // em(k) = 1._core_rknd
          em(i, k) = real_t(4.5);
        } else {
          em(i, k) = em_min;
        }
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, 0) = em(icol, 1);
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }
  } else if (runtype == "atex") {
    // GCSS ATEX
    // Known magic number (Stevens, et al. 2001, eq. 1).
    for (int k = 0; k < um.extent_int(1); ++k) {

      for (int icol = 0; icol < um.extent_int(0); ++icol) {
        um(icol, k) = std::max({um(icol, k), (-real_t(8.))});
      }
    }

    //          do k = 1, gr%nz
    //             em(k) = 1.0_core_rknd - (gr%zm(1,k)/3000.0_core_rknd)
    //             if ( em(k) < em_min ) then
    //                em(k) = em_min
    //             end if
    //          end do
    //          em(1) = em(2)
    //          em(gr%nz) = em(gr%nz-1)
    for (int k = 0; k < em.extent_int(1); ++k) {

      for (int icol = 0; icol < em.extent_int(0); ++icol) {
        em(icol, k) = em_min;
      }
    }
  } else if (runtype == "dycoms2_rf01") {
    // GCSS DYCOMS II RF01
    cloud_top_height = real_t(800.);

    for (i = 0; i < ngrdcol; ++i) {

      for (k = 0; k < gr.nzm; ++k) {

        if (gr_zm_host(i, k) < cloud_top_height) {
          // em(k) = 0.5_core_rknd
          em(i, k) = real_t(1.1);
        } else {
          em(i, k) = em_min;
        }
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, 0) = em(icol, 1);
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }
  } else if (runtype == "dycoms2_rf02") {
    // GCSS DYCOMS II RF02

    for (int k = 0; k < em.extent_int(1); ++k) {

      for (int icol = 0; icol < em.extent_int(0); ++icol) {
        em(icol, k) = real_t(1.0);
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }
  } else if (runtype == "nov11_altocu") {
    // Brian for Nov. 11 altocumulus case.
    // Vince Larson reduced initial forcing.  4 Nov 2005
    // em = 1.0_core_rknd
    // em = 0.1_core_rknd
    // 4150 + 2800 m is the top of the cloud in Nov11
    cloud_top_height = (real_t(2800.) + gr_zm_host(0, 0));

    for (i = 0; i < ngrdcol; ++i) {

      for (k = 0; k < gr.nzm; ++k) {

        if (gr_zm_host(i, k) < cloud_top_height) {
          //  Modification by Adam Smith, 08 April 2008
          //  Reducing the value of em appears to reduce error in the
          //  updated Nov.11 case
          em(i, k) = real_t(0.01);
        } else {
          em(i, k) = em_min;
        }
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, 0) = em(icol, 1);
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }

    // End Vince Larson's change.

  } else if (runtype == "jun25_altocu") {
    // Adam Smith addition for June 25 altocumulus case.

    // Note: em_min was 1.0e-6, defined in constants.F (Adam Smith, 28 June 2006).
    // Now em_min = 1.5 * w_tol_sqd (Brian Griffin, 26 Nov. 2008).
    for (k = 0; k < gr.nzm; ++k) {

      for (int icol = 0; icol < em.extent_int(0); ++icol) {
        em(icol, k) = real_t(0.01);
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, 0) = em(icol, 1);
    }

    // End Vince Larson's change.
    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }
  } else if (runtype == "clex9_nov02") {
    // Adam Smith addition for CLEX-9: Nov. 02 altocumulus case.
    // Vince Larson reduced initial forcing.  4 Nov 2005
    // em = 1.0_core_rknd
    // em = 0.1_core_rknd
    // 4150 + 1400 m is the top of the cloud in Nov11
    cloud_top_height = (real_t(2200.) + gr_zm_host(0, 0));

    for (i = 0; i < ngrdcol; ++i) {

      for (k = 0; k < gr.nzm; ++k) {

        if (gr_zm_host(i, k) < cloud_top_height) {
          em(i, k) = real_t(0.01);
        } else {
          em(i, k) = em_min;
        }
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, 0) = em(icol, 1);
    }

    // End Vince Larson's change.
    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }
  } else if (runtype == "clex9_oct14") {
    // Adam Smith addition for CLEX-9: Oct. 14 altocumulus case.
    // Vince Larson reduced initial forcing.  4 Nov 2005
    // em = 1.0_core_rknd
    // em = 0.1_core_rknd
    // 4150 + 1400 m is the top of the cloud in Nov11
    cloud_top_height = (real_t(3500.) + gr_zm_host(0, 0));

    for (i = 0; i < ngrdcol; ++i) {

      for (k = 0; k < gr.nzm; ++k) {

        if (gr_zm_host(i, k) < cloud_top_height) {

          for (int icol = 0; icol < em.extent_int(0); ++icol) {
            em(icol, k) = real_t(0.01);
          }
        } else {

          for (int icol = 0; icol < em.extent_int(0); ++icol) {
            em(icol, k) = em_min;
          }
        }
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, 0) = em(icol, 1);
    }

    // End Vince Larson's change.
    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }
  } else if (runtype == "lba") {

    for (int k = 0; k < em.extent_int(1); ++k) {

      for (int icol = 0; icol < em.extent_int(0); ++icol) {
        em(icol, k) = real_t(0.1);
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }
  } else if (runtype == "mpace_a") {
    // Michael Falk for mpace_a Arctic Stratus case.
    cloud_top_height = real_t(2000.);

    em_max = real_t(1.0);

    for (i = 0; i < ngrdcol; ++i) {

      for (k = 0; k < gr.nzm; ++k) {

        if (gr_zm_host(i, k) < cloud_top_height) {

          for (int icol = 0; icol < em.extent_int(0); ++icol) {
            em(icol, k) = em_max;
          }
        } else {

          for (int icol = 0; icol < em.extent_int(0); ++icol) {
            em(icol, k) = em_min;
          }
        }
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, 0) = em(icol, 1);
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }

    mpace_a_init(iunit, forcings_file_path, // In
                 mpace);                    // Out

  } else if (runtype == "mpace_b") {
    // Michael Falk for mpace_b Arctic Stratus case.
    cloud_top_height = real_t(1300.);

    //  Michael Falk 17 Aug 2006
    em_max = real_t(1.0);

    for (i = 0; i < ngrdcol; ++i) {

      for (k = 0; k < gr.nzm; ++k) {

        if (gr_zm_host(i, k) < cloud_top_height) {

          for (int icol = 0; icol < em.extent_int(0); ++icol) {
            em(icol, k) = em_max;
          }
        } else {

          for (int icol = 0; icol < em.extent_int(0); ++icol) {
            em(icol, k) = em_min;
          }
        }
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, 0) = em(icol, 1);
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }
  } else if (runtype == "cobra") {
    // Brian Griffin for COBRA CO2 case.

    for (int k = 0; k < em.extent_int(1); ++k) {

      for (int icol = 0; icol < em.extent_int(0); ++icol) {
        em(icol, k) = real_t(0.1);
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }
  } else if (runtype == "rico") {
    // Michael Falk for RICO tropical cumulus case, 13 Dec 2006
    cloud_top_height = real_t(1500.);

    em_max = real_t(1.0);

    for (i = 0; i < ngrdcol; ++i) {

      for (k = 0; k < gr.nzm; ++k) {

        if (gr_zm_host(i, k) < cloud_top_height) {

          for (int icol = 0; icol < em.extent_int(0); ++icol) {
            em(icol, k) = em_max;
          }
        } else {

          for (int icol = 0; icol < em.extent_int(0); ++icol) {
            em(icol, k) = em_min;
          }
        }
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, 0) = em(icol, 1);
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }
  } else if (runtype == "gabls2") {
    // Michael Falk for GABLS2 case, 29 Dec 2006
    cloud_top_height = real_t(800.);
    em_max           = real_t(0.5);

    for (i = 0; i < ngrdcol; ++i) {

      for (k = 0; k < gr.nzm; ++k) {

        if (gr_zm_host(i, k) < cloud_top_height) {

          for (int icol = 0; icol < em.extent_int(0); ++icol) {
            em(icol, k) = (em_max * (real_t(1.) - (gr_zm_host(i, k) / cloud_top_height)));
          }
        } else {

          for (int icol = 0; icol < em.extent_int(0); ++icol) {
            em(icol, k) = em_min;
          }
        }
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, 0) = em(icol, 1);
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }
  } else if (runtype == "gabls3_night") {

    for (int k = 0; k < em.extent_int(1); ++k) {

      for (int icol = 0; icol < em.extent_int(0); ++icol) {
        em(icol, k) = real_t(1.0);
      }
    }
  } else if (runtype == "gabls3") {

    for (int k = 0; k < em.extent_int(1); ++k) {

      for (int icol = 0; icol < em.extent_int(0); ++icol) {
        em(icol, k) = real_t(1.0);
      }
    }

    for (int icol = 0; icol < em.extent_int(0); ++icol) {
      em(icol, (gr.nzm - 1)) = em_min;
    }

    for (int icol = 0; icol < veg_T_in_K.extent_int(0); ++icol) {
      veg_T_in_K(icol) = real_t(300.);
    }

    for (int icol = 0; icol < sfc_soil_T_in_K.extent_int(0); ++icol) {
      sfc_soil_T_in_K(icol) = real_t(300.);
    }

    for (int icol = 0; icol < deep_soil_T_in_K.extent_int(0); ++icol) {
      deep_soil_T_in_K(icol) = real_t(288.58);
    }
  } else if (runtype == "coriolis_test") {
    // Set the amplitude of the harmonic oscillator in the "coriolis_test".
    // The maximum amplitude is 2 * w_tol_sqd at the middle vertical level.
    // Hing Ong, 8 August 2025

    for (i = 0; i < ngrdcol; ++i) {
      domain_depth = (gr_zm_host(i, (gr.nzm - 1)) - gr_zm_host(i, 0));

      for (k = 0; k < gr.nzm; ++k) {
        em(i, k) = ((std::sin(((pi * gr_zm_host(i, k)) / domain_depth)) * w_tol_sqd) * real_t(6.0));
      }
    }
  } else {

    for (int k = 0; k < em.extent_int(1); ++k) {

      for (int icol = 0; icol < em.extent_int(0); ++icol) {
        em(icol, k) = em_min;
      }
    }
  }

  // End Initialize TKE and other fields as needed.
  // Initialize w'^2 based on initial TKE:
  // TKE: em = (1/2) * (w'^2 + u'^2 + v'^2).
  // Evenly divide TKE into its component contributions (w'^2, u'^2, and v'^2).
  for (int k = 0; k < wp2.extent_int(1); ++k) {

    for (int icol = 0; icol < wp2.extent_int(0); ++icol) {
      wp2(icol, k) = ((real_t(2.0) / real_t(3.0)) * em(icol, k));
    }
  }

  for (int k = 0; k < up2.extent_int(1); ++k) {

    for (int icol = 0; icol < up2.extent_int(0); ++icol) {
      up2(icol, k) = ((real_t(2.0) / real_t(3.0)) * em(icol, k));
    }
  }

  for (int k = 0; k < vp2.extent_int(1); ++k) {

    for (int icol = 0; icol < vp2.extent_int(0); ++icol) {
      vp2(icol, k) = ((real_t(2.0) / real_t(3.0)) * em(icol, k));
    }
  }

  for (int k = 0; k < upwp.extent_int(1); ++k) {

    for (int icol = 0; icol < upwp.extent_int(0); ++icol) {
      upwp(icol, k) = zero;
    }
  }

  // The "coriolis_test" is a benchmark against analytic solutions where up2, wp2, and upwp
  // interact as a harmonic oscillator analogous to the Foucault pendulum.
  // Here are the solutions to test the implementation of the nontraditional Coriolis terms:
  //              upwp = amplitude(grdcol,zm) * sin(2 * fcor_y * time)
  // (up2 - wp2) / 2 = amplitude(grdcol,zm) * cos(2 * fcor_y * time)
  // Here, set up2, wp2, upwp at time = 0.
  // The setting can be modified to test the implementation of the traditional Coriolis terms,
  // where the solutions are as follows:
  // upwp = amplitude(grdcol,zm) * cos( fcor * time)
  // vpwp = amplitude(grdcol,zm) * sin(-fcor * time)
  // Hing Ong, 8 August 2025
  if (runtype == "coriolis_test") {

    for (int k = 0; k < wp2.extent_int(1); ++k) {

      for (int icol = 0; icol < wp2.extent_int(0); ++icol) {
        wp2(icol, k) = (((real_t(1.0) / real_t(3.0)) * em(icol, k)) + w_tol_sqd);
      }
    }

    for (int k = 0; k < up2.extent_int(1); ++k) {

      for (int icol = 0; icol < up2.extent_int(0); ++icol) {
        up2(icol, k) = (((real_t(3.0) / real_t(3.0)) * em(icol, k)) + w_tol_sqd);
      }
    }

    for (int k = 0; k < vp2.extent_int(1); ++k) {

      for (int icol = 0; icol < vp2.extent_int(0); ++icol) {
        vp2(icol, k) = (((real_t(2.0) / real_t(3.0)) * em(icol, k)) + w_tol_sqd);
      }
    }

    for (int k = 0; k < em.extent_int(1); ++k) {

      for (int icol = 0; icol < em.extent_int(0); ++icol) {
        em(icol, k) = (em(icol, k) + (real_t(1.5) * w_tol_sqd));
      }
    }

    // Advance upwp by half a time step like the leapfrog time integration.
    // This improves the benchmark at time = 0, where the tendency is large for
    // upwp but small for up2 and wp2 (Hing Ong, 8 August 2025).
    for (i = 0; i < ngrdcol; ++i) {

      for (k = 0; k < gr.nzm; ++k) {
        upwp(i, k) = (((real_t(0.5) * dt_main) * fcor_y(i)) * (up2(i, k) - wp2(i, k)));
      }
    }
  }

  // Moved this to be more general -dschanen July 16 2007.
  if (clubb_config_flags.l_uv_nudge || uv_sponge_damp_settings.l_sponge_damping) {

    for (int k = 0; k < um_ref.extent_int(1); ++k) {

      for (int icol = 0; icol < um_ref.extent_int(0); ++icol) {
        um_ref(icol, k) = um(icol, k);
      }
    }

    for (int k = 0; k < vm_ref.extent_int(1); ++k) {

      for (int icol = 0; icol < vm_ref.extent_int(0); ++icol) {
        vm_ref(icol, k) = vm(icol, k);
      }
    }
  } else {

    for (int k = 0; k < um_ref.extent_int(1); ++k) {

      for (int icol = 0; icol < um_ref.extent_int(0); ++icol) {
        um_ref(icol, k) = zero;
      }
    }

    for (int k = 0; k < vm_ref.extent_int(1); ++k) {

      for (int icol = 0; icol < vm_ref.extent_int(0); ++icol) {
        vm_ref(icol, k) = zero;
      }
    }
  }

  for (int k = 0; k < thlm_ref.extent_int(1); ++k) {

    for (int icol = 0; icol < thlm_ref.extent_int(0); ++icol) {
      thlm_ref(icol, k) = zero;
    }
  }

  for (int k = 0; k < rtm_ref.extent_int(1); ++k) {

    for (int icol = 0; icol < rtm_ref.extent_int(0); ++icol) {
      rtm_ref(icol, k) = zero;
    }
  }

  return;
}

void ClubbDriver::initialize_clubb_variables(
    int ngrdcol, int sclr_dim, int edsclr_dim, const sclr_idx_type& sclr_idx,       // In
    const grid& gr, const std::string& altitude_type,                               // In
    const std::string& temperature_type,                                            // In
    DriverSlice1 p_sfc, DriverSlice1 rtm_sfc, DriverSlice2 rtm, bool l_readiopdata, // In
    int saturation_formula,                                                         // In
    DriverSlice2 thlm, DriverSlice2 p_in_Pa,                                        // InOut
    DriverSlice2 p_in_Pa_zm,                                                        // Out
    DriverSlice2 exner, DriverSlice2 rho, DriverSlice2 rho_zm,                      // Out
    DriverSlice2 rcm, DriverSlice2 thvm, DriverSlice2 rho_ds_zm,                    // Out
    DriverSlice2 rho_ds_zt, DriverSlice2 invrs_rho_ds_zm,                           // Out
    DriverSlice2 invrs_rho_ds_zt, DriverSlice2 thv_ds_zm,                           // Out
    DriverSlice2 thv_ds_zt,                                                         // Out
    DriverSlice3 sclrm, DriverSlice3 edsclrm) {                                     // InOut

  // Convert sounding thermodynamics and compute the dry anelastic reference state.
  auto workspace_scope = workspace.scope();

  // Local variables.
  DriverArray1 pd_sfc = workspace_scope.view<DriverArray1>(
      "pd_sfc", ngrdcol); // Dry surface pressure                [Pa]
  DriverArray1 rv_sfc = workspace_scope.view<DriverArray1>(
      "rv_sfc", ngrdcol); // Surface water vapor mixing ratio    [kg/kg]
  DriverArray2 thm = workspace_scope.view<DriverArray2>(
      "thm", ngrdcol, gr.nzt); // Potential temperature (thermodynamic levels)   [K]
  DriverArray2 th_dry = workspace_scope.view<DriverArray2>(
      "th_dry", ngrdcol, gr.nzt); // Dry potential temperature (thermo. levels)     [K]
  DriverArray2 p_dry = workspace_scope.view<DriverArray2>(
      "p_dry", ngrdcol, gr.nzt); // Dry air pressure (thermodynamic levels)        [Pa]
  DriverArray2 exner_dry = workspace_scope.view<DriverArray2>(
      "exner_dry", ngrdcol, gr.nzt); // Exner of dry air (thermodynamic levels)        [-]
  DriverArray2 rho_dry = workspace_scope.view<DriverArray2>(
      "rho_dry", ngrdcol, gr.nzt); // Dry air density (thermodynamic levels)         [kg/m^3]
  DriverArray2 tmp      = workspace_scope.view<DriverArray2>("tmp", ngrdcol, gr.nzm);
  DriverArray2 thm_zm   = workspace_scope.view<DriverArray2>("thm_zm", ngrdcol, gr.nzm);
  DriverArray2 exner_zm = workspace_scope.view<DriverArray2>(
      "exner_zm", ngrdcol, gr.nzm); // Exner on momentum levels                       [-]
  DriverArray2 th_dry_zm = workspace_scope.view<DriverArray2>(
      "th_dry_zm", ngrdcol, gr.nzm); // Dry potential temperature on momentum levels   [K]
  DriverArray2 p_dry_zm = workspace_scope.view<DriverArray2>(
      "p_dry_zm", ngrdcol, gr.nzm); // Dry air pressure on momentum levels            [Pa]
  DriverArray2 exner_dry_zm = workspace_scope.view<DriverArray2>(
      "exner_dry_zm", ngrdcol, gr.nzm); // Exner of dry air on momentum levels            [-]
  DriverArray2 rho_dry_zm = workspace_scope.view<DriverArray2>(
      "rho_dry_zm", ngrdcol, gr.nzm); // Dry air density on momentum levels             [kg/m^3]
  int i{};
  int k{};

  for (i = 0; i < ngrdcol; ++i) {

    if (rtm_sfc(i) < real_t(0.0)) {
      // The sounding doesn't extended to the surface, so rtm_sfc is set to a
      // negative number.  Use rtm(1) as rv_sfc.
      rv_sfc(i) = rtm(i, 0);
    } else {
      // The sounding does extended to the surface, so rtm_sfc is the initial value
      // of total water mixing ratio at the surface.
      rv_sfc(i) = rtm_sfc(i);
    }

    // Calculate dry surface pressure from surface pressure and surface water
    // vapor mixing ratio, such that p_d = p / [ 1 + (R_v/R_d)*r_v ].
    pd_sfc(i) = (p_sfc(i) / (one + (ep2 * rv_sfc(i))));
  }

  if (temperature_type == temperature_name) {

    if (altitude_type == z_name) {
      std::cerr << "Interpetation of sounding files with z as the " << " "
                << "independent variable and absolute temperature " << " "
                << "as the temperature variable has not been " << " "
                << "implemented.  Either specify pressure as the " << " "
                << "independent variable or thm/thlm as the " << " "
                << "temperature variable." << "\n";
      throw std::runtime_error("Fatal error.");
    } else if (altitude_type == pressure_name) {

      for (k = 0; k < gr.nzt; ++k) {

        for (i = 0; i < ngrdcol; ++i) {
          thlm(i, k) = (thlm(i, k) / std::pow((p_in_Pa(i, k) / p0), kappa));
        }
      }
    } else {
      throw std::runtime_error("Invalid sounding vertical-coordinate variable");
    }
  }

  for (k = 0; k < gr.nzt; ++k) {

    for (i = 0; i < ngrdcol; ++i) {
      thvm(i, k) = (thlm(i, k) * (one + (ep1 * (rtm(i, k) / (one + rtm(i, k))))));
    }
  }

  //  Compute approximate pressure, exner, and density using an approximate
  //  value of theta_v.
  hydrostatic(workspace, ngrdcol, gr, thvm, p_sfc, p_in_Pa, p_in_Pa_zm, exner, exner_zm, rho,
              rho_zm);

  if (temperature_type == theta_name || temperature_type == temperature_name) {

    for (k = 0; k < gr.nzt; ++k) {

      for (i = 0; i < ngrdcol; ++i) {
        thm(i, k) = thlm(i, k);
      }
    }

    for (k = 0; k < gr.nzt; ++k) {

      for (i = 0; i < ngrdcol; ++i) {
        rcm(i, k) =
            std::max({(rtm(i, k) - sat_mixrat_liq_api(p_in_Pa(i, k), (thm(i, k) * exner(i, k)),
                                                      saturation_formula)),
                      zero_threshold});
      }
    }

    for (k = 0; k < gr.nzt; ++k) {

      for (i = 0; i < ngrdcol; ++i) {
        thlm(i, k) = (thlm(i, k) - ((Lv / (Cp * exner(i, k))) * rcm(i, k)));
      }
    }

    if (sclr_idx.iisclr_thl > -1) {

      for (k = 0; k < gr.nzt; ++k) {

        for (i = 0; i < ngrdcol; ++i) {
          sclrm(i, k, sclr_idx.iisclr_thl) = thlm(i, k);
        }
      }
    }

    if (sclr_idx.iiedsclr_thl > -1) {

      for (k = 0; k < gr.nzt; ++k) {

        for (i = 0; i < ngrdcol; ++i) {
          edsclrm(i, k, sclr_idx.iiedsclr_thl) = thlm(i, k);
        }
      }
    }
  } else if (temperature_type == thetal_name) {

    for (k = 0; k < gr.nzt; ++k) {

      for (i = 0; i < ngrdcol; ++i) {
        // Compute cloud water mixing ratio using an iterative method.
        rcm(i, k) =
            rcm_sat_adj(thlm(i, k), rtm(i, k), p_in_Pa(i, k), exner(i, k), saturation_formula);
      }
    }

    for (k = 0; k < gr.nzt; ++k) {

      for (i = 0; i < ngrdcol; ++i) {
        thm(i, k) = (thlm(i, k) + ((Lv / (Cp * exner(i, k))) * rcm(i, k)));
      }
    }

    if (sclr_idx.iisclr_thl > -1) {

      for (k = 0; k < gr.nzt; ++k) {

        for (i = 0; i < ngrdcol; ++i) {
          sclrm(i, k, sclr_idx.iisclr_thl) = thlm(i, k);
        }
      }
    }

    if (sclr_idx.iiedsclr_thl > -1) {

      for (k = 0; k < gr.nzt; ++k) {

        for (i = 0; i < ngrdcol; ++i) {
          edsclrm(i, k, sclr_idx.iiedsclr_thl) = thlm(i, k);
        }
      }
    }
  } else {
    std::cerr << "Invalid temperature_type: " << " " << temperature_type << "\n";
    throw std::runtime_error("Fatal error");
  }

  //  Now, compute initial thvm, given initial thm, rtm, and rcm.
  //
  //  The exact form of the equaton for theta_v (with only water vapor
  //  included) is:
  //
  //  theta_v = theta * [ ( 1 + (R_v/R_d)*r_v ) / ( 1 + r_v ) ];
  //
  //  which can be rearranged as:
  //
  //  theta_v = theta * [ 1 + { (R_v/R_d) - 1 } * { r_v / ( 1 + r_v ) } ].
  //
  //  The exact form of the equation for theta_v (including water vapor and
  //  cloud water) is:
  //
  //  theta_v = theta * [ ( 1 + (R_v/R_d)*r_v ) / ( 1 + r_v + r_c ) ];
  //
  //  which can be rearranged as:
  //
  //  theta_v = theta * [ 1 + { (R_v/R_d) - 1 } * { r_v / ( 1 + r_v + r_c ) }
  //                      - { r_c / ( 1 + r_v + r_c ) } ].
  //
  //  This version is written with r_t = r_v + r_c, such that:
  //
  //  theta_v = theta * [ 1 + { (R_v/R_d) - 1 }
  //                          * { ( r_t - r_c ) / ( 1 + r_t ) }
  //                      - { r_c / ( 1 + r_t ) } ].
  //
  //  To use theta_l instead of theta, simply substitute the following for
  //  theta in the above expression:  theta_l + {L_v/(C_p*exner)} * r_c.
  //
  //  The CLUBB code uses a linearized version of the above equation (in order
  //  to calculate thv' terms -- such as w'thv', etc.) for theta_v throughout
  //  the model code, such that:
  //
  //  theta_v = theta_l + { (R_v/R_d) - 1 } * thv_ds * r_t
  //                    + [ {L_v/(C_p*exner)} - (R_v/R_d) * thv_ds ] * r_c;
  //
  //  where thv_ds is used as a reference value to approximate theta_l.
  DriverArray2 thv_ds_zt_1 =
      workspace_scope.view<DriverArray2>("thv_ds_zt_1", thm.extent_int(0), thm.extent_int(1));

  for (int k = 0; k < thm.extent_int(1); ++k) {

    for (int icol = 0; icol < thm.extent_int(0); ++icol) {
      thv_ds_zt_1(icol, k) =
          (thm(icol, k) * std::pow((one + (ep2 * (rtm(icol, k) - rcm(icol, k)))), kappa));
    }
  }

  // Now, compute initial thvm, given initial thm, rtm, and rcm.
  // The exact form of the equaton for theta_v (with only water vapor
  // included) is:
  // theta_v = theta * [ ( 1 + (R_v/R_d)*r_v ) / ( 1 + r_v ) ];
  // which can be rearranged as:
  // theta_v = theta * [ 1 + { (R_v/R_d) - 1 } * { r_v / ( 1 + r_v ) } ].
  // The exact form of the equation for theta_v (including water vapor and
  // cloud water) is:
  // theta_v = theta * [ ( 1 + (R_v/R_d)*r_v ) / ( 1 + r_v + r_c ) ];
  // which can be rearranged as:
  // theta_v = theta * [ 1 + { (R_v/R_d) - 1 } * { r_v / ( 1 + r_v + r_c ) }
  // - { r_c / ( 1 + r_v + r_c ) } ].
  // This version is written with r_t = r_v + r_c, such that:
  // theta_v = theta * [ 1 + { (R_v/R_d) - 1 }
  // * { ( r_t - r_c ) / ( 1 + r_t ) }
  // - { r_c / ( 1 + r_t ) } ].
  // To use theta_l instead of theta, simply substitute the following for
  // theta in the above expression:  theta_l + {L_v/(C_p*exner)} * r_c.
  // The CLUBB code uses a linearized version of the above equation (in order
  // to calculate thv' terms -- such as w'thv', etc.) for theta_v throughout
  // the model code, such that:
  // theta_v = theta_l + { (R_v/R_d) - 1 } * thv_ds * r_t
  // + [ {L_v/(C_p*exner)} - (R_v/R_d) * thv_ds ] * r_c;
  // where thv_ds is used as a reference value to approximate theta_l.
  calculate_thvm(gr.nzt, ngrdcol, thlm, rtm, rcm, exner, thv_ds_zt_1, thvm);

  //  Compute approximate pressure, exner, and density using an approximate
  //  value of theta_v.
  hydrostatic(workspace, ngrdcol, gr, thvm, p_sfc, p_in_Pa, p_in_Pa_zm, exner, exner_zm, rho,
              rho_zm);

  for (k = 0; k < gr.nzt; ++k) {

    for (i = 0; i < ngrdcol; ++i) {
      // Calculate dry pressure from total pressure and water vapor mixing
      // ratio, such that:  p_d = p / [ 1 + (R_v/R_d)*r_v ].
      p_dry(i, k) = (p_in_Pa(i, k) / (one + (ep2 * (rtm(i, k) - rcm(i, k)))));

      // Calculate dry exner from dry pressure.
      exner_dry(i, k) = std::pow((p_dry(i, k) / p0), kappa);
    }
  }

  for (k = 0; k < gr.nzt; ++k) {

    for (i = 0; i < ngrdcol; ++i) {
      th_dry(i, k) = (thm(i, k) * std::pow((one + (ep2 * (rtm(i, k) - rcm(i, k)))), kappa));
    }
  }

  for (k = 0; k < gr.nzt; ++k) {

    for (i = 0; i < ngrdcol; ++i) {
      rho_dry(i, k) = (p_dry(i, k) / ((Rd * th_dry(i, k)) * exner_dry(i, k)));
    }
  }

  for (i = 0; i < ngrdcol; ++i) {
    p_dry_zm(i, 0) = pd_sfc(i);
  }

  DriverArray2 azt_2 =
      workspace_scope.view<DriverArray2>("azt_2", rtm.extent_int(0), rtm.extent_int(1));

  for (int k = 0; k < rtm.extent_int(1); ++k) {

    for (int icol = 0; icol < rtm.extent_int(0); ++icol) {
      azt_2(icol, k) = (rtm(icol, k) - rcm(icol, k));
    }
  }

  zt2zm_api(gr.nzm, gr.nzt, ngrdcol, gr, azt_2, tmp);

  for (k = 1; k < gr.nzm; ++k) {

    for (i = 0; i < ngrdcol; ++i) {
      // Calculate dry pressure on momentum levels from total pressure (on
      // momentum levels) and water vapor mixing ratio (interpolated to
      // momentum levels), such that:  p_d = p / [ 1 + (R_v/R_d)*r_v ].
      p_dry_zm(i, k) = (p_in_Pa_zm(i, k) / (one + (ep2 * std::max({tmp(i, k), zero_threshold}))));
    }
  }

  for (k = 0; k < gr.nzm; ++k) {

    for (i = 0; i < ngrdcol; ++i) {
      // Calculate dry exner on momentum levels from dry pressure on momentum
      // levels.
      exner_dry_zm(i, k) = std::pow((p_dry_zm(i, k) / p0), kappa);
    }
  }

  zt2zm_api(gr.nzm, gr.nzt, ngrdcol, gr, thm, thm_zm);

  for (k = 0; k < gr.nzm; ++k) {

    for (i = 0; i < ngrdcol; ++i) {
      th_dry_zm(i, k) =
          (thm_zm(i, k) * std::pow((one + (ep2 * std::max({tmp(i, k), zero_threshold}))), kappa));
    }
  }

  for (k = 0; k < gr.nzm; ++k) {

    for (i = 0; i < ngrdcol; ++i) {
      rho_dry_zm(i, k) = (p_dry_zm(i, k) / ((Rd * th_dry_zm(i, k)) * exner_dry_zm(i, k)));
    }
  }

  // The values of rho_dry and rho_dry_zm that were just calculated are dry
  // (they do not take into account water vapor or cloud water).  These are
  // the values of dry, static, base-state density that are needed for the
  // anelastic equation set.
  for (int k = 0; k < rho_ds_zt.extent_int(1); ++k) {

    for (int icol = 0; icol < rho_ds_zt.extent_int(0); ++icol) {
      rho_ds_zt(icol, k) = rho_dry(icol, k);
    }
  }

  for (int k = 0; k < rho_ds_zm.extent_int(1); ++k) {

    for (int icol = 0; icol < rho_ds_zm.extent_int(0); ++icol) {
      rho_ds_zm(icol, k) = rho_dry_zm(icol, k);
    }
  }

  // Since theta_d does not include water in any phase, the value of dry,
  // static, base-state theta is the same as the values of both dry, static,
  // base-state theta_l and dry, static, base-state theta_v.  Thus, for use
  // with the anelastic equation set, thv_ds = thl_ds = theta_ds.
  for (int k = 0; k < thv_ds_zt.extent_int(1); ++k) {

    for (int icol = 0; icol < thv_ds_zt.extent_int(0); ++icol) {
      thv_ds_zt(icol, k) = th_dry(icol, k);
    }
  }

  for (int k = 0; k < thv_ds_zm.extent_int(1); ++k) {

    for (int icol = 0; icol < thv_ds_zm.extent_int(0); ++icol) {
      thv_ds_zm(icol, k) = th_dry_zm(icol, k);
    }
  }

  // At this point, the values of the dry, static, base-state variables
  // rho_ds_zt, rho_ds_zm, thv_ds_zt, and thv_ds_zm have been calculated.
  // The CLUBB code is set up to be an anelastic code.  If use of the
  // Boussinesq approximation is desired, rather than the anelastic
  // approximation, reset the values of rho_ds_zt and rho_ds_zm to 1.  Also,
  // reset the values of thv_ds_zt and thv_ds_zm to reference temperature T0.
  if (l_use_boussinesq) {

    for (int k = 0; k < rho_ds_zt.extent_int(1); ++k) {

      for (int icol = 0; icol < rho_ds_zt.extent_int(0); ++icol) {
        rho_ds_zt(icol, k) = one;
      }
    }

    for (int k = 0; k < rho_ds_zm.extent_int(1); ++k) {

      for (int icol = 0; icol < rho_ds_zm.extent_int(0); ++icol) {
        rho_ds_zm(icol, k) = one;
      }
    }

    for (int k = 0; k < thv_ds_zt.extent_int(1); ++k) {

      for (int icol = 0; icol < thv_ds_zt.extent_int(0); ++icol) {
        thv_ds_zt(icol, k) = T0;
      }
    }

    for (int k = 0; k < thv_ds_zm.extent_int(1); ++k) {

      for (int icol = 0; icol < thv_ds_zm.extent_int(0); ++icol) {
        thv_ds_zm(icol, k) = T0;
      }
    }
  }

  // Set the values of inverse, dry, static, base-state density.
  for (int k = 0; k < invrs_rho_ds_zm.extent_int(1); ++k) {

    for (int icol = 0; icol < invrs_rho_ds_zm.extent_int(0); ++icol) {
      invrs_rho_ds_zm(icol, k) = (one / rho_ds_zm(icol, k));
    }
  }

  for (int k = 0; k < invrs_rho_ds_zt.extent_int(1); ++k) {

    for (int icol = 0; icol < invrs_rho_ds_zt.extent_int(0); ++icol) {
      invrs_rho_ds_zt(icol, k) = (one / rho_ds_zt(icol, k));
    }
  }

  return;
}

void ClubbDriver::restart_clubb(const grid& gr,                       // In
                                int hydromet_dim,                     // In
                                const hm_metadata_type& hm_metadata,  // In
                                const std::string& restart_path_case, // In
                                double time_restart,                  // In
                                DriverSlice1 um,                      // InOut
                                DriverSlice1 upwp,                    // InOut
                                DriverSlice1 vm,                      // InOut
                                DriverSlice1 vpwp,                    // InOut
                                DriverSlice1 up2,                     // InOut
                                DriverSlice1 vp2,                     // InOut
                                DriverSlice1 rtm,                     // InOut
                                DriverSlice1 wprtp,                   // InOut
                                DriverSlice1 thlm,                    // InOut
                                DriverSlice1 wpthlp,                  // InOut
                                DriverSlice1 rtp2,                    // InOut
                                DriverSlice1 rtp3,                    // InOut
                                DriverSlice1 thlp2,                   // InOut
                                DriverSlice1 thlp3,                   // InOut
                                DriverSlice1 rtpthlp,                 // InOut
                                DriverSlice1 wp2,                     // InOut
                                DriverSlice1 wp3,                     // InOut
                                DriverSlice1 p_in_Pa,                 // InOut
                                DriverSlice1 exner,                   // InOut
                                DriverSlice1 rcm,                     // InOut
                                DriverSlice1 cloud_frac,              // InOut
                                DriverSlice1 wpthvp,                  // InOut
                                DriverSlice1 wp2thvp,                 // InOut
                                DriverSlice1 wp2up,                   // InOut
                                DriverSlice1 rtpthvp,                 // InOut
                                DriverSlice1 thlpthvp,                // InOut
                                DriverSlice1 wp2rtp,                  // InOut
                                DriverSlice1 wp2thlp,                 // InOut
                                DriverSlice1 uprcp,                   // InOut
                                DriverSlice1 vprcp,                   // InOut
                                DriverSlice1 rc_coef_zm,              // InOut
                                DriverSlice1 wp4,                     // InOut
                                DriverSlice1 wpup2,                   // InOut
                                DriverSlice1 wpvp2,                   // InOut
                                DriverSlice1 wp2up2,                  // InOut
                                DriverSlice1 wp2vp2,                  // InOut
                                DriverSlice1 ice_supersat_frac,       // InOut
                                DriverSlice1 wm_zt,                   // InOut
                                DriverSlice1 rho,                     // InOut
                                DriverSlice1 rho_zm,                  // InOut
                                DriverSlice1 rho_ds_zm,               // InOut
                                DriverSlice1 rho_ds_zt,               // InOut
                                DriverSlice1 thv_ds_zm,               // InOut
                                DriverSlice1 thv_ds_zt,               // InOut
                                DriverSlice1 thlm_forcing,            // InOut
                                DriverSlice1 rtm_forcing,             // InOut
                                DriverSlice1 wprtp_forcing,           // InOut
                                DriverSlice1 wpthlp_forcing,          // InOut
                                DriverSlice1 rtp2_forcing,            // InOut
                                DriverSlice1 thlp2_forcing,           // InOut
                                DriverSlice1 rtpthlp_forcing,         // InOut
                                DriverSlice2 hydromet,                // InOut
                                DriverSlice2 hydrometp2,              // InOut
                                DriverSlice2 wphydrometp,             // InOut
                                DriverSlice1 Ncm,                     // InOut
                                DriverSlice1 Nccnm,                   // InOut
                                DriverSlice1 thvm,                    // InOut
                                DriverSlice1 em,                      // InOut
                                DriverSlice1 tau_zm,                  // InOut
                                DriverSlice1 tau_zt,                  // InOut
                                DriverSlice1 Kh_zt,                   // InOut
                                DriverSlice1 Kh_zm,                   // InOut
                                DriverSlice1 ug,                      // InOut
                                DriverSlice1 vg,                      // InOut
                                DriverSlice1 thlprcp,                 // InOut
                                DriverSlice1 sigma_sqd_w,             // InOut
                                DriverSlice1 sigma_sqd_w_zt,          // InOut
                                DriverSlice1 radht,                   // InOut
                                real_t& deep_soil_T_in_K,             // InOut
                                real_t& sfc_soil_T_in_K,              // InOut
                                real_t& veg_T_in_K,                   // InOut
                                pdf_parameter& pdf_params,            // InOut
                                pdf_parameter& pdf_params_zm,         // InOut
                                DriverSlice1 rcm_mc,                  // Out
                                DriverSlice1 rvm_mc,                  // Out
                                DriverSlice1 thlm_mc,                 // Out
                                DriverSlice1 wprtp_mc,                // Out
                                DriverSlice1 wpthlp_mc,               // Out
                                DriverSlice1 rtp2_mc,                 // Out
                                DriverSlice1 thlp2_mc,                // Out
                                DriverSlice1 rtpthlp_mc,              // Out
                                real_t& wpthlp_sfc,                   // Out
                                real_t& wprtp_sfc,                    // Out
                                real_t& upwp_sfc,                     // Out
                                real_t& vpwp_sfc) {                   // Out

  // Configure case-owned field readers and restore a saved time on the current grid.
  auto& l_soil_veg       = soil_vegetation.l_soil_veg;
  auto& microphys_scheme = microphys.microphys_scheme;
  auto& l_predict_Nc     = microphys.l_predict_Nc;
  auto& l_ice_microphys  = microphys.l_ice_microphys;
  auto& l_graupel        = microphys.l_graupel;
  auto& stat_files       = inputfields.stat_files;

  // Local variables.
  int timestep{};
  bool l_restart{};
  bool l_read_error{};
  bool l_fatal_error{};
  auto gr_zm_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, gr.zm);
  auto gr_zt_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, gr.zt);

  //  --- Begin Code ---
  l_fatal_error = false;

  // Inform inputfields module
  inputfields.input_type            = "clubb";
  inputfields.l_input_um            = true;
  inputfields.l_input_vm            = true;
  inputfields.l_input_rtm           = true;
  inputfields.l_input_thlm          = true;
  inputfields.l_input_wp2           = true;
  inputfields.l_input_ug            = true;
  inputfields.l_input_vg            = true;
  inputfields.l_input_rcm           = true;
  inputfields.l_input_wm_zt         = true;
  inputfields.l_input_exner         = true;
  inputfields.l_input_em            = true;
  inputfields.l_input_p             = true;
  inputfields.l_input_rho           = true;
  inputfields.l_input_rho_zm        = true;
  inputfields.l_input_rho_ds_zm     = true;
  inputfields.l_input_rho_ds_zt     = true;
  inputfields.l_input_thv_ds_zm     = true;
  inputfields.l_input_thv_ds_zt     = true;
  inputfields.l_input_Lscale        = true;
  inputfields.l_input_Lscale_up     = true;
  inputfields.l_input_Lscale_down   = true;
  inputfields.l_input_Kh_zt         = true;
  inputfields.l_input_Kh_zm         = true;
  inputfields.l_input_tau_zm        = true;
  inputfields.l_input_tau_zt        = true;
  inputfields.l_input_thvm          = true;
  inputfields.l_input_wpthvp        = true;
  inputfields.l_input_wp2thvp       = true;
  inputfields.l_input_wp2up         = true;
  inputfields.l_input_rtpthvp       = true;
  inputfields.l_input_thlpthvp      = true;
  inputfields.l_input_wp2rtp        = true;
  inputfields.l_input_wp2thlp       = true;
  inputfields.l_input_uprcp         = true;
  inputfields.l_input_vprcp         = true;
  inputfields.l_input_rc_coef_zm    = true;
  inputfields.l_input_wp4           = true;
  inputfields.l_input_wpup2         = true;
  inputfields.l_input_wpvp2         = true;
  inputfields.l_input_wp2up2        = true;
  inputfields.l_input_wp2vp2        = true;
  inputfields.l_input_iss_frac      = true;
  inputfields.l_input_w_1           = true;
  inputfields.l_input_w_2           = true;
  inputfields.l_input_varnce_w_1    = true;
  inputfields.l_input_varnce_w_2    = true;
  inputfields.l_input_rt_1          = true;
  inputfields.l_input_rt_2          = true;
  inputfields.l_input_varnce_rt_1   = true;
  inputfields.l_input_varnce_rt_2   = true;
  inputfields.l_input_thl_1         = true;
  inputfields.l_input_thl_2         = true;
  inputfields.l_input_varnce_thl_1  = true;
  inputfields.l_input_varnce_thl_2  = true;
  inputfields.l_input_mixt_frac     = true;
  inputfields.l_input_chi_1         = true;
  inputfields.l_input_chi_2         = true;
  inputfields.l_input_stdev_chi_1   = true;
  inputfields.l_input_stdev_chi_2   = true;
  inputfields.l_input_rc_1          = true;
  inputfields.l_input_rc_2          = true;
  inputfields.l_input_w_1_zm        = true;
  inputfields.l_input_w_2_zm        = true;
  inputfields.l_input_varnce_w_1_zm = true;
  inputfields.l_input_varnce_w_2_zm = true;
  inputfields.l_input_mixt_frac_zm  = true;
  inputfields.l_input_radht         = true;

  if (microphys_scheme == "coamps") {
    inputfields.l_input_rrm   = true;
    inputfields.l_input_rsm   = true;
    inputfields.l_input_rim   = true;
    inputfields.l_input_rgm   = true;
    inputfields.l_input_Nccnm = true;
    inputfields.l_input_Ncm   = true;
    inputfields.l_input_Nrm   = true;
    inputfields.l_input_Nim   = true;
  } else if (microphys_scheme == "morrison") {
    inputfields.l_input_rrm = true;
    inputfields.l_input_Nrm = true;

    if (l_ice_microphys) {
      inputfields.l_input_rsm = true;
      inputfields.l_input_rim = true;
      inputfields.l_input_Nim = true;

      if (l_graupel) {
        inputfields.l_input_rgm = true;
      } else {
        inputfields.l_input_rgm = false;
      }
    } else {
      inputfields.l_input_rsm = false;
      inputfields.l_input_rim = false;
      inputfields.l_input_Nim = false;
      inputfields.l_input_rgm = false;
    }

    inputfields.l_input_Nccnm = false;

    if (l_predict_Nc) {
      inputfields.l_input_Ncm = true;
    } else {
      inputfields.l_input_Ncm = false;
    }
  } else if (microphys_scheme == "khairoutdinov_kogan") {
    inputfields.l_input_rrm   = true;
    inputfields.l_input_rsm   = false;
    inputfields.l_input_rim   = false;
    inputfields.l_input_rgm   = false;
    inputfields.l_input_Nccnm = false;
    inputfields.l_input_Ncm   = false;
    inputfields.l_input_Nrm   = true;
    inputfields.l_input_Nim   = false;
  } else {
    inputfields.l_input_rrm   = false;
    inputfields.l_input_rsm   = false;
    inputfields.l_input_rim   = false;
    inputfields.l_input_rgm   = false;
    inputfields.l_input_Nccnm = false;
    inputfields.l_input_Ncm   = false;
    inputfields.l_input_Nrm   = false;
    inputfields.l_input_Nim   = false;
  }

  if (l_soil_veg) {
    inputfields.l_input_veg_T_in_K       = true;
    inputfields.l_input_deep_soil_T_in_K = true;
    inputfields.l_input_sfc_soil_T_in_K  = true;
  }

  inputfields.l_input_wprtp           = true;
  inputfields.l_input_wpthlp          = true;
  inputfields.l_input_wp3             = true;
  inputfields.l_input_rtp2            = true;
  inputfields.l_input_rtp3            = true;
  inputfields.l_input_thlp2           = true;
  inputfields.l_input_thlp3           = true;
  inputfields.l_input_rtpthlp         = true;
  inputfields.l_input_upwp            = true;
  inputfields.l_input_vpwp            = true;
  inputfields.l_input_thlm_forcing    = true;
  inputfields.l_input_rtm_forcing     = true;
  inputfields.l_input_up2             = true;
  inputfields.l_input_vp2             = true;
  inputfields.l_input_sigma_sqd_w     = true;
  inputfields.l_input_cloud_frac      = true;
  inputfields.l_input_sigma_sqd_w_zt  = true;
  inputfields.l_input_wprtp_forcing   = true;
  inputfields.l_input_wpthlp_forcing  = true;
  inputfields.l_input_rtp2_forcing    = true;
  inputfields.l_input_thlp2_forcing   = true;
  inputfields.l_input_rtpthlp_forcing = true;
  inputfields.l_input_thlprcp         = true;
  inputfields.l_input_rcm_mc          = true;
  inputfields.l_input_rvm_mc          = true;
  inputfields.l_input_thlm_mc         = true;
  inputfields.l_input_wprtp_mc        = true;
  inputfields.l_input_wpthlp_mc       = true;
  inputfields.l_input_rtp2_mc         = true;
  inputfields.l_input_thlp2_mc        = true;
  inputfields.l_input_rtpthlp_mc      = true;
  set_filenames(("../" + restart_path_case), // In
                inputfields);                // InOut

  //  Determine the nearest timestep in the netCDF file to the
  //  restart time.
  l_restart = true;
  compute_timestep(workspace, text_at(stat_files, 0), l_restart, time_restart, // In
                   timestep,                                                   // Out
                   inputfields);                                               // In

  // Sanity check for input time_restart
  if (timestep < 0) {
    std::cerr << "Invalid time_restart" << "\n";
    throw std::runtime_error("Fatal error");
  }

  // Read data from stats files
  stat_fields_reader(workspace, gr, timestep, hydromet_dim, hm_metadata, // In
                     microphys_scheme, l_predict_Nc,                     // In
                     um, upwp, vm, vpwp,                                 // InOut
                     up2, vp2, rtm,                                      // InOut
                     wprtp, thlm, wpthlp,                                // InOut
                     rtp2, rtp3,                                         // InOut
                     thlp2, thlp3, rtpthlp,                              // InOut
                     wp2, wp3,                                           // InOut
                     p_in_Pa, exner, rcm, cloud_frac,                    // InOut
                     wpthvp, wp2thvp, wp2up, rtpthvp,                    // InOut
                     thlpthvp,                                           // InOut
                     wp2rtp, wp2thlp, uprcp, vprcp,                      // InOut
                     rc_coef_zm, wp4, wpup2,                             // InOut
                     wpvp2, wp2up2,                                      // InOut
                     wp2vp2, ice_supersat_frac,                          // InOut
                     wm_zt, rho, rho_zm, rho_ds_zm,                      // InOut
                     rho_ds_zt, thv_ds_zm, thv_ds_zt,                    // InOut
                     thlm_forcing, rtm_forcing, wprtp_forcing,           // InOut
                     wpthlp_forcing, rtp2_forcing,                       // InOut
                     thlp2_forcing, rtpthlp_forcing,                     // InOut
                     hydromet, hydrometp2, wphydrometp,                  // InOut
                     Ncm, Nccnm, thvm, em,                               // InOut
                     tau_zm, tau_zt,                                     // InOut
                     Kh_zt, Kh_zm, ug, vg,                               // InOut
                     thlprcp,                                            // InOut
                     sigma_sqd_w, sigma_sqd_w_zt, radht,                 // InOut
                     deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K,      // InOut
                     pdf_params, pdf_params_zm,                          // InOut
                     inputfields);                                       // In

  get_clubb_variable_interpolated(workspace, true,                             // In
                                  text_at(stat_files, 0),                      // In
                                  "rcm_mc",                                    // In
                                  gr.nzt,                                      // In
                                  timestep,                                    // In
                                  Kokkos::subview(gr_zt_host, 0, Kokkos::ALL), // In
                                  rcm_mc,                                      // InOut
                                  l_read_error,                                // Out
                                  inputfields);                                // In

  l_fatal_error = (l_fatal_error || l_read_error);

  get_clubb_variable_interpolated(workspace, true,                             // In
                                  text_at(stat_files, 0),                      // In
                                  "rvm_mc",                                    // In
                                  gr.nzt,                                      // In
                                  timestep,                                    // In
                                  Kokkos::subview(gr_zt_host, 0, Kokkos::ALL), // In
                                  rvm_mc,                                      // InOut
                                  l_read_error,                                // Out
                                  inputfields);                                // In

  l_fatal_error = (l_fatal_error || l_read_error);

  get_clubb_variable_interpolated(workspace, true,                             // In
                                  text_at(stat_files, 0),                      // In
                                  "thlm_mc",                                   // In
                                  gr.nzt,                                      // In
                                  timestep,                                    // In
                                  Kokkos::subview(gr_zt_host, 0, Kokkos::ALL), // In
                                  thlm_mc,                                     // InOut
                                  l_read_error,                                // Out
                                  inputfields);                                // In

  l_fatal_error = (l_fatal_error || l_read_error);

  get_clubb_variable_interpolated(workspace, true,                             // In
                                  text_at(stat_files, 1),                      // In
                                  "wprtp_mc",                                  // In
                                  gr.nzm,                                      // In
                                  timestep,                                    // In
                                  Kokkos::subview(gr_zm_host, 0, Kokkos::ALL), // In
                                  wprtp_mc,                                    // InOut
                                  l_read_error,                                // Out
                                  inputfields);                                // In

  l_fatal_error = (l_fatal_error || l_read_error);

  get_clubb_variable_interpolated(workspace, true,                             // In
                                  text_at(stat_files, 1),                      // In
                                  "wpthlp_mc",                                 // In
                                  gr.nzm,                                      // In
                                  timestep,                                    // In
                                  Kokkos::subview(gr_zm_host, 0, Kokkos::ALL), // In
                                  wpthlp_mc,                                   // InOut
                                  l_read_error,                                // Out
                                  inputfields);                                // In

  l_fatal_error = (l_fatal_error || l_read_error);

  get_clubb_variable_interpolated(workspace, true,                             // In
                                  text_at(stat_files, 1),                      // In
                                  "rtp2_mc",                                   // In
                                  gr.nzm,                                      // In
                                  timestep,                                    // In
                                  Kokkos::subview(gr_zm_host, 0, Kokkos::ALL), // In
                                  rtp2_mc,                                     // InOut
                                  l_read_error,                                // Out
                                  inputfields);                                // In

  l_fatal_error = (l_fatal_error || l_read_error);

  get_clubb_variable_interpolated(workspace, true,                             // In
                                  text_at(stat_files, 1),                      // In
                                  "thlp2_mc",                                  // In
                                  gr.nzm,                                      // In
                                  timestep,                                    // In
                                  Kokkos::subview(gr_zm_host, 0, Kokkos::ALL), // In
                                  thlp2_mc,                                    // InOut
                                  l_read_error,                                // Out
                                  inputfields);                                // In

  l_fatal_error = (l_fatal_error || l_read_error);

  get_clubb_variable_interpolated(workspace, true,                             // In
                                  text_at(stat_files, 1),                      // In
                                  "rtpthlp_mc",                                // In
                                  gr.nzm,                                      // In
                                  timestep,                                    // In
                                  Kokkos::subview(gr_zm_host, 0, Kokkos::ALL), // In
                                  rtpthlp_mc,                                  // InOut
                                  l_read_error,                                // Out
                                  inputfields);                                // In

  l_fatal_error = (l_fatal_error || l_read_error);
  wpthlp_sfc    = wpthlp(0);
  wprtp_sfc     = wprtp(0);
  upwp_sfc      = upwp(0);
  vpwp_sfc      = vpwp(0);
  return;
}

} // namespace clubb_kokkos
