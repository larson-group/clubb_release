// Initialize and release one case's supported microphysics configuration.
//
// init_microphys restores Native defaults, reads the four Native namelist groups,
// configures KK/Morrison and SILHS, initializes hydrometeor/PDF metadata and
// correlation matrices, then prints the same setup diagnostics. cleanup_microphys
// releases the species sedimentation View. Mutable package state is appended to
// the unchanged Native argument prefix and belongs to the calling ClubbDriver.
// Namelist assignments and setup logic use call-local C++ readers.
// Supported schemes/options match JAX: none, KK and Morrison; COAMPS, simplified
// ice, GFDL activation and the named JAX-excluded Morrison selectors are rejected.
// References: Morrison et al.(2005), JAS62,1665–1677; KK(2000), MWR128,229–243.
#pragma once
#include "scratch_workspace.hpp"
#include "Microphys/parameters_microphys.hpp"
#include "Microphys/KK_microphys/parameters_KK.hpp"
#include "Microphys/Morrison_microphys/module_mp_graupel.hpp"
#include "SILHS/parameters_silhs.hpp"
#include "CLUBB_core/corr_varnce_module.hpp"
#include "CLUBB_core/setup_clubb_pdf_params.hpp"
#include "CLUBB_core/matrix_operations.hpp"
#include "CLUBB_core/parameter_indices.hpp"
#include "CLUBB_core/pdf_utilities.hpp"
#include "CLUBB_core/index_mapping.hpp"
#include "CLUBB_core/error_code.hpp"
#include "CLUBB_core/precipitation_fraction.hpp"
#include "Input_fields/corr_varnce_input_reader.hpp"
#include "Microphys/microphys_io.hpp"
#include "driver_arrays.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
namespace clubb_kokkos {
inline void init_microphys(ScratchWorkspace& workspace, int iunit, const std::string& runtype,
                           const std::string& namelist_file, const std::string& case_info_file,
                           real_t host_dx, real_t host_dy, Real1D clubb_params,
                           bool l_diagnose_correlations, bool& l_const_Nc_in_cloud,
                           bool& l_fix_w_chi_eta_correlations, int& hydromet_dim, int& pdf_dim,
                           hm_metadata_type& hm_metadata,
                           silhs_config_flags_type& silhs_config_flags,
                           real_t& vert_decorr_coef_out, Real2D& corr_array_n_cloud,
                           Real2D& corr_array_n_below, parameters_microphys& microphys,
                           parameters_KK& kk, morrison_microphysics_state& morrison,
                           parameters_silhs_state& silhs) {
  // Set hydrometeor species indices and dimensions for allocation. Restore all
  // namelist defaults on every initialization so repeated cases cannot inherit
  // settings. Appended references replace the implicit Native package modules.
  auto workspace_scope                  = workspace.scope();
  constexpr bool l_write_to_file        = true;
  constexpr const char* corr_input_path = "../input/case_data/";
  constexpr const char* cloud_file_ext  = "_corr_array_cloud.in";
  constexpr const char* below_file_ext  = "_corr_array_below.in";
  std::string lh_microphys_type         = "disabled";
  hmp2_ip_on_hmm2_ip_slope_type hmp2_ip_on_hmm2_ip_slope;
  hmp2_ip_on_hmm2_ip_intrcpt_type hmp2_ip_on_hmm2_ip_intrcpt;
  real_t Ncnp2_on_Ncnm2 = one;
  int iirr, iirs, iiri, iirg, iiNr, iiNs, iiNi, iiNg;
  constexpr int one_lev = 1;
  int cluster_allocation_strategy;
  bool l_lh_importance_sampling, l_Lscale_vert_avg, l_lh_straight_mc, l_lh_clustered_sampling;
  bool l_rcm_in_cloud_k_lh_start, l_random_k_lh_start, l_max_overlap_in_cloud;
  bool l_lh_instant_var_covar_src, l_lh_limit_weights, l_lh_var_frac, l_lh_normalize_weights;
  bool l_lh_deterministic_test;
  int read_status = 0;
  std::ofstream case_log;
  // Restore declaration defaults before reading each case's namelists.
  // C++ owns explicit case records instead of Native saved module variables.
  microphys                  = parameters_microphys{};
  kk                         = parameters_KK{};
  silhs                      = parameters_silhs_state{};
  hmp2_ip_on_hmm2_ip_slope   = hmp2_ip_on_hmm2_ip_slope_type{};
  hmp2_ip_on_hmm2_ip_intrcpt = hmp2_ip_on_hmm2_ip_intrcpt_type{};
  // GFDL settings are parsed even when disabled, as in Native. The activation
  // engine itself remains excluded with the agreed JAX option surface.
  bool nooc         = false;
  real_t sul_concen = real_t(.1), low_concen = real_t(.1), high_concen = one;
  real_t lowup = real_t(.3), highup = real_t(10.), lowup2 = real_t(.01), highup2 = real_t(.3);
  real_t lowmass2 = real_t(.01), highmass2 = real_t(100.), lowmass3 = real_t(.01),
         highmass3 = real_t(100.);
  real_t lowmass4 = real_t(.01), highmass4 = real_t(100.), lowmass5 = real_t(.01),
         highmass5 = real_t(100.);
  real_t lowT2 = real_t(243.15), highT2 = real_t(308.15), aeromass_value = real_t(2.25e-12);
  // Parameters for NNUCCD/NNUCCC coefficients in clex9_oct14 (Marc Pilon,2011).
  morrison.NNUCCD_REDUCE_COEF = 1.0f;
  morrison.NNUCCC_REDUCE_COEF = 1.0f;
  if (microphys_detail::trim(runtype) == "clex9_oct14") {
    morrison.NNUCCD_REDUCE_COEF = .01f;
    morrison.NNUCCC_REDUCE_COEF = .01f;
  }
  // Aerosol for RF02 from Mikhail Ovtchinnikov: preserve default-real rounding.
  morrison.aer_rm1    = .011e-6f;
  morrison.aer_rm2    = .06e-6f;
  morrison.aer_sig1   = 1.2f;
  morrison.aer_sig2   = 1.7f;
  morrison.aer_n1     = 125.e6f;
  morrison.aer_n2     = 65.e6f;
  morrison.ccnconst   = 120.f;
  morrison.ccnexpnt   = .4f;
  morrison.pgam_fixed = 5.f;
  // Replace defaults with the case's Native microphysics namelist entries.
  read_microphysics_setting(
      namelist_file, iunit, read_status, microphys.microphys_scheme, microphys.l_cloud_sed,
      microphys.sigma_g, microphys.l_ice_microphys, microphys.l_graupel, microphys.l_hail,
      microphys.l_var_covar_src, microphys.l_upwind_diff_sed, microphys.l_seifert_beheng,
      microphys.l_predict_Nc, l_const_Nc_in_cloud, microphys.specify_aerosol, microphys.l_subgrid_w,
      microphys.l_arctic_nucl, microphys.l_cloud_edge_activation, microphys.l_fix_pgam,
      microphys.l_in_cloud_Nc_diff, lh_microphys_type, microphys.l_local_kk,
      microphys.lh_num_samples, microphys.lh_sequence_length, microphys.lh_seed,
      l_fix_w_chi_eta_correlations, microphys.l_silhs_KK_convergence_adj_mean,
      hmp2_ip_on_hmm2_ip_slope.rr, hmp2_ip_on_hmm2_ip_slope.Nr, hmp2_ip_on_hmm2_ip_slope.ri,
      hmp2_ip_on_hmm2_ip_slope.Ni, hmp2_ip_on_hmm2_ip_slope.rs, hmp2_ip_on_hmm2_ip_slope.Ns,
      hmp2_ip_on_hmm2_ip_slope.rg, hmp2_ip_on_hmm2_ip_slope.Ng, hmp2_ip_on_hmm2_ip_intrcpt.rr,
      hmp2_ip_on_hmm2_ip_intrcpt.Nr, hmp2_ip_on_hmm2_ip_intrcpt.ri, hmp2_ip_on_hmm2_ip_intrcpt.Ni,
      hmp2_ip_on_hmm2_ip_intrcpt.rs, hmp2_ip_on_hmm2_ip_intrcpt.Ns, hmp2_ip_on_hmm2_ip_intrcpt.rg,
      hmp2_ip_on_hmm2_ip_intrcpt.Ng, Ncnp2_on_Ncnm2, kk.C_evap, kk.r_0,
      microphys.microphys_start_time, microphys.Nc0_in_cloud, morrison.ccnconst, morrison.ccnexpnt,
      morrison.aer_rm1, morrison.aer_rm2, morrison.aer_n1, morrison.aer_n2, morrison.aer_sig1,
      morrison.aer_sig2, morrison.pgam_fixed);
  if (read_status != 0)
    throw std::runtime_error("Unable to read &microphysics_setting namelist");
  set_default_silhs_config_flags_api(
      cluster_allocation_strategy, l_lh_importance_sampling, l_Lscale_vert_avg, l_lh_straight_mc,
      l_lh_clustered_sampling, l_rcm_in_cloud_k_lh_start, l_random_k_lh_start,
      l_max_overlap_in_cloud, l_lh_instant_var_covar_src, l_lh_limit_weights, l_lh_var_frac,
      l_lh_normalize_weights, l_lh_deterministic_test);
  // Read SILHS parameters only when sampled microphysics is enabled.
  if (microphys_detail::trim(lh_microphys_type) != "disabled") {
    read_silhs_params_nl(namelist_file, iunit, read_status,
                              silhs.eight_cluster_presc_probs.cloud_precip_comp1,
                              silhs.eight_cluster_presc_probs.cloud_precip_comp2,
                              silhs.eight_cluster_presc_probs.nocloud_precip_comp1,
                              silhs.eight_cluster_presc_probs.nocloud_precip_comp2,
                              silhs.eight_cluster_presc_probs.cloud_noprecip_comp1,
                              silhs.eight_cluster_presc_probs.cloud_noprecip_comp2,
                              silhs.eight_cluster_presc_probs.nocloud_noprecip_comp1,
                              silhs.eight_cluster_presc_probs.nocloud_noprecip_comp2,
                              silhs.vert_decorr_coef, silhs.importance_prob_thresh);
    if (read_status != 0)
      throw std::runtime_error("Unable to read &silhs_params_nl namelist");
    read_configurable_silhs_flags_nl(
        namelist_file, iunit, read_status, cluster_allocation_strategy, l_lh_importance_sampling,
        l_Lscale_vert_avg, l_lh_straight_mc, l_lh_clustered_sampling, l_rcm_in_cloud_k_lh_start,
        l_random_k_lh_start, l_max_overlap_in_cloud, l_lh_instant_var_covar_src, l_lh_limit_weights,
        l_lh_var_frac, l_lh_normalize_weights, l_lh_deterministic_test);
    if (read_status != 0)
      throw std::runtime_error("Unable to read &configurable_silhs_flags_nl namelist");
  }
  // Deterministic inputs do not replace independent importance/start-level draws.
  if (microphys_detail::trim(lh_microphys_type) != "disabled" && l_lh_deterministic_test &&
      (l_lh_importance_sampling || l_random_k_lh_start))
    throw std::runtime_error(
        "Deterministic SILHS testing requires importance sampling and random starts disabled");
  initialize_silhs_config_flags_type_api(
      cluster_allocation_strategy, l_lh_importance_sampling, l_Lscale_vert_avg, l_lh_straight_mc,
      l_lh_clustered_sampling, l_rcm_in_cloud_k_lh_start, l_random_k_lh_start,
      l_max_overlap_in_cloud, l_lh_instant_var_covar_src, l_lh_limit_weights, l_lh_var_frac,
      l_lh_normalize_weights, l_lh_deterministic_test, silhs_config_flags);
  vert_decorr_coef_out = silhs.vert_decorr_coef;
  // Printing microphysics inputs, in Native namelist field order.
  if (clubb_at_least_debug_level_api(1)) {
    if (l_write_to_file)
      case_log.open(case_info_file, std::ios::app);
    write_text("--------------------------------------------------", case_log);
    write_text("&microphysics_setting", case_log);
    write_text("--------------------------------------------------", case_log);
    write_text("microphys_scheme = " + microphys.microphys_scheme, case_log);
    write_text("l_cloud_sed = ", microphys.l_cloud_sed, case_log);
    write_text("sigma_g = ", microphys.sigma_g, case_log);
    write_text("l_graupel = ", microphys.l_graupel, case_log);
    write_text("l_hail = ", microphys.l_hail, case_log);
    write_text("l_seifert_beheng = ", microphys.l_seifert_beheng, case_log);
    write_text("l_predict_Nc = ", microphys.l_predict_Nc, case_log);
    write_text("l_const_Nc_in_cloud = ", l_const_Nc_in_cloud, case_log);
    write_text("specify_aerosol = " + microphys.specify_aerosol, case_log);
    write_text("l_subgrid_w = ", microphys.l_subgrid_w, case_log);
    write_text("l_arctic_nucl = ", microphys.l_arctic_nucl, case_log);
    write_text("l_cloud_edge_activation = ", microphys.l_cloud_edge_activation, case_log);
    write_text("l_fix_pgam = ", microphys.l_fix_pgam, case_log);
    write_text("l_in_cloud_Nc_diff = ", microphys.l_in_cloud_Nc_diff, case_log);
    write_text("l_var_covar_src = ", microphys.l_var_covar_src, case_log);
    write_text("l_upwind_diff_sed = ", microphys.l_upwind_diff_sed, case_log);
    write_text("lh_microphys_type = " + lh_microphys_type, case_log);
    write_text("lh_num_samples = ", microphys.lh_num_samples, case_log);
    write_text("lh_sequence_length = ", microphys.lh_sequence_length, case_log);
    write_text("lh_seed = ", microphys.lh_seed, case_log);
    write_text("l_fix_w_chi_eta_correlations = ", l_fix_w_chi_eta_correlations, case_log);
    write_text("l_silhs_KK_convergence_adj_mean = ", microphys.l_silhs_KK_convergence_adj_mean,
               case_log);
    write_text("importance_prob_thresh = ", silhs.importance_prob_thresh, case_log);
    write_text("host_dx = ", host_dx, case_log);
    write_text("host_dy = ", host_dy, case_log);
    write_text("hmp2_ip_on_hmm2_ip_slope%rr = ", hmp2_ip_on_hmm2_ip_slope.rr, case_log);
    write_text("hmp2_ip_on_hmm2_ip_slope%ri = ", hmp2_ip_on_hmm2_ip_slope.ri, case_log);
    write_text("hmp2_ip_on_hmm2_ip_slope%rs = ", hmp2_ip_on_hmm2_ip_slope.rs, case_log);
    write_text("hmp2_ip_on_hmm2_ip_slope%rg = ", hmp2_ip_on_hmm2_ip_slope.rg, case_log);
    write_text("hmp2_ip_on_hmm2_ip_slope%Nr = ", hmp2_ip_on_hmm2_ip_slope.Nr, case_log);
    write_text("hmp2_ip_on_hmm2_ip_slope%Ni = ", hmp2_ip_on_hmm2_ip_slope.Ni, case_log);
    write_text("hmp2_ip_on_hmm2_ip_slope%Ni = ", hmp2_ip_on_hmm2_ip_slope.Ni, case_log);
    write_text("hmp2_ip_on_hmm2_ip_slope%Ng = ", hmp2_ip_on_hmm2_ip_slope.Ng, case_log);
    write_text("hmp2_ip_on_hmm2_ip_intrcpt%rr = ", hmp2_ip_on_hmm2_ip_intrcpt.rr, case_log);
    write_text("hmp2_ip_on_hmm2_ip_intrcpt%ri = ", hmp2_ip_on_hmm2_ip_intrcpt.ri, case_log);
    write_text("hmp2_ip_on_hmm2_ip_intrcpt%rs = ", hmp2_ip_on_hmm2_ip_intrcpt.rs, case_log);
    write_text("hmp2_ip_on_hmm2_ip_intrcpt%rg = ", hmp2_ip_on_hmm2_ip_intrcpt.rg, case_log);
    write_text("hmp2_ip_on_hmm2_ip_intrcpt%Nr = ", hmp2_ip_on_hmm2_ip_intrcpt.Nr, case_log);
    write_text("hmp2_ip_on_hmm2_ip_intrcpt%Ni = ", hmp2_ip_on_hmm2_ip_intrcpt.Ni, case_log);
    write_text("hmp2_ip_on_hmm2_ip_intrcpt%Ni = ", hmp2_ip_on_hmm2_ip_intrcpt.Ni, case_log);
    write_text("hmp2_ip_on_hmm2_ip_intrcpt%Ng = ", hmp2_ip_on_hmm2_ip_intrcpt.Ng, case_log);
    write_text("Ncnp2_on_Ncnm2 = ", Ncnp2_on_Ncnm2, case_log);
    write_text("C_evap = ", kk.C_evap, case_log);
    write_text("r_0 = ", kk.r_0, case_log);
    write_text("microphys_start_time = ", static_cast<real_t>(microphys.microphys_start_time),
               case_log);
    write_text("Nc0_in_cloud = ", microphys.Nc0_in_cloud, case_log);
    write_text("ccnconst = ", static_cast<real_t>(morrison.ccnconst), case_log);
    write_text("ccnexpnt = ", static_cast<real_t>(morrison.ccnexpnt), case_log);
    write_text("aer_rm1 = ", static_cast<real_t>(morrison.aer_rm1), case_log);
    write_text("aer_rm2 = ", static_cast<real_t>(morrison.aer_rm2), case_log);
    write_text("aer_n1 = ", static_cast<real_t>(morrison.aer_n1), case_log);
    write_text("aer_n2 = ", static_cast<real_t>(morrison.aer_n2), case_log);
    write_text("aer_sig1 = ", static_cast<real_t>(morrison.aer_sig1), case_log);
    write_text("aer_sig2 = ", static_cast<real_t>(morrison.aer_sig2), case_log);
    write_text("pgam_fixed = ", static_cast<real_t>(morrison.pgam_fixed), case_log);
    write_text("precip_frac_calc_type = ", precip_frac_calc_type, case_log);
    write_text("--------------------------------------------------", case_log);
    write_text("&SILHS_setting", case_log);
    write_text("--------------------------------------------------", case_log);
    print_silhs_config_flags_api(std::cout, silhs_config_flags);
    if (l_write_to_file)
      print_silhs_config_flags_api(case_log, silhs_config_flags);
    if (l_write_to_file)
      case_log.close();
  }
  read_gfdl_activation_setting(
      namelist_file, iunit, read_status, nooc, sul_concen, low_concen, high_concen, lowup, highup,
      lowup2, highup2, lowmass2, highmass2, lowmass3, highmass3, lowmass4, highmass4, lowmass5,
      highmass5, lowT2, highT2, aeromass_value, microphys.l_gfdl_activation);
  if (read_status != 0)
    throw std::runtime_error("Unable to read &gfdl_activation_setting namelist");
  // TODO: GFDL activation (Loading/aer_ccn_act_k_init) is excluded with JAX.
  // Reject its namelist selector before allocating lookup tables or timestepping.
  if (microphys.l_gfdl_activation)
    throw std::runtime_error("GFDL activation core is outside the supported JAX/Kokkos options");
  // Hydrometeor positions need not match between schemes if indices agree.
  if (microphys_detail::trim(microphys.microphys_scheme) == "morrison") {
    iirr = 0;
    iiNr = 1;
    if (microphys.l_ice_microphys) {
      iiri                = 2;
      iiNi                = 3;
      iirs                = 4;
      iiNs                = 5;
      morrison.doicemicro = true;
      if (microphys.l_graupel) {
        iirg               = 6;
        iiNg               = 7;
        hydromet_dim       = 8;
        morrison.dograupel = true;
      } else {
        iirg               = -2;
        iiNg               = -2;
        hydromet_dim       = 6;
        morrison.dograupel = false;
      }
    } else {
      iirs                = -2;
      iiri                = -2;
      iiNs                = -2;
      iiNi                = -2;
      iirg                = -2;
      iiNg                = -2;
      hydromet_dim        = 2;
      morrison.doicemicro = false;
      morrison.dograupel  = false;
    }
    // Nc0 is in cm^-3 inside Morrison; retain Native real4 conversion here.
    morrison.Nc0         = static_cast<float>(microphys.Nc0_in_cloud / cm3_per_m3);
    morrison.dopredictNc = microphys.l_predict_Nc;
    if (microphys_detail::trim(microphys.specify_aerosol) == "morrison_no_aerosol")
      morrison.aerosol_mode = morrison_no_aerosol;
    else if (microphys_detail::trim(microphys.specify_aerosol) == "morrison_power_law")
      morrison.aerosol_mode = morrison_power_law;
    else if (microphys_detail::trim(microphys.specify_aerosol) == "morrison_lognormal")
      morrison.aerosol_mode = morrison_lognormal;
    else
      throw std::runtime_error("Unknown Morrison aerosol mode.");
    morrison.docloudedgeactivation = microphys.l_cloud_edge_activation;
    morrison.doarcticicenucl       = microphys.l_arctic_nucl;
    morrison.dohail                = microphys.l_hail;
    morrison.dosb_warm_rain        = microphys.l_seifert_beheng;
    morrison.dofix_pgam            = microphys.l_fix_pgam;
    morrison.dosubgridw            = microphys.l_subgrid_w;
    if (microphys.l_cloud_sed) {
      std::cerr << "Morrison microphysics has seperate code for cloud water sedimentation, "
                   "therefore l_cloud_sed should be set to .false.\n";
      throw std::runtime_error("Fatal error.");
    }
    if (!l_fix_w_chi_eta_correlations && microphys.l_ice_microphys &&
        microphys_detail::trim(lh_microphys_type) != "disabled") {
      std::cerr << "The flag l_fix_w_chi_eta_correlations must be true in order to enable latin "
                   "hypercube sampling and ice microphysics.\n";
      std::cerr << "The flag l_ice_microphys must be set to false to use this option.\n";
      throw std::runtime_error("Fatal error.");
    }
    // Sedimentation is handled within Morrison's numerical core.
    microphys.l_hydromet_sed = Bool1D("init_microphys::l_hydromet_sed", hydromet_dim);
    Kokkos::deep_copy(microphys.l_hydromet_sed, std::uint8_t(0));
    // The same options are excluded by the actual JAX initialization gate.
    if (microphys.l_hail || microphys.l_seifert_beheng || microphys.l_arctic_nucl ||
        microphys.l_cloud_edge_activation || microphys.l_fix_pgam)
      throw std::runtime_error(
          "Requested Morrison selector is outside the supported JAX/Kokkos options");
    GRAUPEL_INIT(workspace, morrison);
  } else if (microphys_detail::trim(microphys.microphys_scheme) == "coamps") {
    // TODO: COAMPS is outside the agreed JAX scheme surface.
    throw std::runtime_error("COAMPS microphysics is outside the supported JAX/Kokkos options");
  } else if (microphys_detail::trim(microphys.microphys_scheme) == "khairoutdinov_kogan") {
    if (microphys.l_predict_Nc) {
      std::cerr << "Khairoutdinov-Kogan microphysics does not support l_predict_Nc = T\n";
      throw std::runtime_error("Fatal Error");
    }
    iirr                     = 0;
    iirs                     = -2;
    iiri                     = -2;
    iirg                     = -2;
    iiNr                     = 1;
    iiNs                     = -2;
    iiNi                     = -2;
    iiNg                     = -2;
    hydromet_dim             = 2;
    microphys.l_hydromet_sed = Bool1D("init_microphys::l_hydromet_sed", hydromet_dim);
    Kokkos::deep_copy(microphys.l_hydromet_sed, std::uint8_t(1));
  } else if (microphys_detail::trim(microphys.microphys_scheme) == "simplified_ice") {
    // TODO: simplified_ice is excluded by actual JAX initialization.
    throw std::runtime_error("simplified_ice is outside the supported JAX/Kokkos options");
  } else if (microphys_detail::trim(microphys.microphys_scheme) == "none") {
    iirr                   = -2;
    iirs                   = -2;
    iiri                   = -2;
    iirg                   = -2;
    iiNr                   = -2;
    iiNs                   = -2;
    iiNi                   = -2;
    iiNg                   = -2;
    hydromet_dim           = 0;
    microphys.l_predict_Nc = false;
  } else {
    std::cerr << "Unknown microphys_scheme: " << microphys.microphys_scheme << '\n';
    throw std::runtime_error("Fatal error");
  }
  // Set PDF/species indices, hydrometeor metadata and prescribed variance ratios.
  init_pdf_hydromet_arrays_api(host_dx, host_dy, hydromet_dim, iirr, iiNr, iiri, iiNi, iirs, iiNs,
                               iirg, iiNg, Ncnp2_on_Ncnm2, hm_metadata, pdf_dim,
                               hmp2_ip_on_hmm2_ip_slope, hmp2_ip_on_hmm2_ip_intrcpt);
  if (microphys_detail::trim(lh_microphys_type) == "interactive")
    microphys.lh_microphys_type = lh_microphys_interactive;
  else if (microphys_detail::trim(lh_microphys_type) == "non-interactive")
    microphys.lh_microphys_type = lh_microphys_non_interactive;
  else if (microphys_detail::trim(lh_microphys_type) == "disabled")
    microphys.lh_microphys_type = lh_microphys_disabled;
  else
    throw std::runtime_error("Error determining lh_microphys_type");
  // Native incompatible combinations are rejected before allocation/timestepping.
  if (microphys.lh_microphys_type != lh_microphys_disabled &&
      (microphys_detail::trim(microphys.microphys_scheme) == "coamps" ||
       microphys_detail::trim(microphys.microphys_scheme) == "simplified_ice"))
    throw std::runtime_error(
        "LH sampling can not be enabled when using coamps, or simplified_ice microphysics types");
  if (microphys.l_silhs_KK_convergence_adj_mean &&
      microphys_detail::trim(microphys.microphys_scheme) != "khairoutdinov_kogan")
    throw std::runtime_error(
        "l_silhs_KK_convergence_adj_mean requires khairoutdinov_kogan microphysics");
  if (l_diagnose_correlations &&
      microphys_detail::trim(microphys.microphys_scheme) != "khairoutdinov_kogan" &&
      microphys.lh_microphys_type == lh_microphys_disabled) {
    std::cerr << "Error: The diagnose_corr algorithm only works for KK microphysics by now.\n";
    throw std::runtime_error("Fatal error");
  }
  if (!microphys.l_local_kk &&
      microphys_detail::trim(microphys.microphys_scheme) == "khairoutdinov_kogan" &&
      microphys.lh_microphys_type == lh_microphys_interactive) {
    std::cerr << "Error: KK upscaled microphysics (l_local_kk = .false.) and interactive Latin "
                 "Hypercube (lh_microphys_type = interactive) are incompatible.\n";
    throw std::runtime_error("Fatal error");
  }
  // Preserve the Native compile-time checks even though l_morr_xp2_mc is
  // disabled in this standalone configuration and excluded by JAX.
  if (l_morr_xp2_mc && microphys.lh_microphys_type != lh_microphys_disabled) {
    std::cerr
        << "Error: The code to include the effects of rain evaporation on rtp2 and thlp2 in "
           "Morrison microphysics (l_morr_xp2_mc = .true.) and Latin Hypercube are incompatible.\n";
    throw std::runtime_error("Fatal error");
  }
  if (l_morr_xp2_mc && microphys.l_var_covar_src) {
    std::cerr << "Error: The code l_morr_xp2_mc and l_var_covar_src are incompatible, since they "
                 "both are used to determine the effect of microphysics on variances.\n";
    throw std::runtime_error("Fatal error");
  }
  if (l_morr_xp2_mc && l_evaporate_cold_rcm) {
    std::cerr << "Error: l_morr_xp2_mc and l_evaporate_cold_rcm are currently incompatible.\n";
    throw std::runtime_error("Fatal error");
  }
  corr_array_n_cloud = Real2D("init_microphys::corr_array_n_cloud", pdf_dim, pdf_dim);
  corr_array_n_below = Real2D("init_microphys::corr_array_n_below", pdf_dim, pdf_dim);
  const std::string corr_file_path_cloud =
      std::string(corr_input_path) + microphys_detail::trim(runtype) + cloud_file_ext;
  const std::string corr_file_path_below =
      std::string(corr_input_path) + microphys_detail::trim(runtype) + below_file_ext;
  const bool l_corr_file_1_exist = std::filesystem::exists(corr_file_path_cloud);
  const bool l_corr_file_2_exist = std::filesystem::exists(corr_file_path_below);
  if (l_corr_file_1_exist && l_corr_file_2_exist) {
    Real2D corr_array_n_cloud_in =
        workspace_scope.view<Real2D>("init_microphys::corr_array_n_cloud_in", pdf_dim, pdf_dim);
    Real2D corr_array_n_below_in =
        workspace_scope.view<Real2D>("init_microphys::corr_array_n_below_in", pdf_dim, pdf_dim);
    read_correlation_matrix(workspace, iunit, corr_file_path_cloud, pdf_dim, hm_metadata,
                            corr_array_n_cloud_in);
    read_correlation_matrix(workspace, iunit, corr_file_path_below, pdf_dim, hm_metadata,
                            corr_array_n_below_in);
    setup_corr_varnce_array_api(workspace, pdf_dim, hm_metadata, l_fix_w_chi_eta_correlations,
                                corr_array_n_cloud, corr_array_n_below, corr_array_n_cloud_in,
                                corr_array_n_below_in);
  } else {
    if (clubb_at_least_debug_level_api(1)) {
      std::cerr << "Warning: missing correlation input file(s): " << corr_file_path_cloud
                << " and/or " << corr_file_path_below << '\n';
      std::cerr << "The default correlation arrays will be used.\n";
    }
    setup_corr_varnce_array_api(workspace, pdf_dim, hm_metadata, l_fix_w_chi_eta_correlations,
                                corr_array_n_cloud, corr_array_n_below);
  }
  // Print approximate real-space correlations only when zeta_vrnce_rat=0.
  // Emergency component-mean adjustment can still change the variance ratio;
  // use saved NetCDF fields for actual diagnosed correlations.
  const auto clubb_params_host = to_host_space(workspace, clubb_params);
  if (clubb_at_least_debug_level_api(1) &&
      std::abs(clubb_params_host(izeta_vrnce_rat)) < eps &&
      microphys_detail::trim(microphys.microphys_scheme) != "none") {
    Real1D sigma2_on_mu2_ip_cloud =
        workspace_scope.view<Real1D>("init_microphys::sigma2_on_mu2_ip_cloud", pdf_dim);
    Real1D sigma2_on_mu2_ip_below =
        workspace_scope.view<Real1D>("init_microphys::sigma2_on_mu2_ip_below", pdf_dim);
    Real1D sigma_x_n_cloud =
        workspace_scope.view<Real1D>("init_microphys::sigma_x_n_cloud", pdf_dim);
    Real1D sigma_x_n_below =
        workspace_scope.view<Real1D>("init_microphys::sigma_x_n_below", pdf_dim);
    Real2D corr_array_cloud =
        workspace_scope.view<Real2D>("init_microphys::corr_array_cloud", pdf_dim, pdf_dim);
    Real2D corr_array_below =
        workspace_scope.view<Real2D>("init_microphys::corr_array_below", pdf_dim, pdf_dim);
    auto sigma_cloud         = Kokkos::create_mirror_view(sigma2_on_mu2_ip_cloud);
    auto sigma_below         = Kokkos::create_mirror_view(sigma2_on_mu2_ip_below);
    auto width_cloud         = Kokkos::create_mirror_view(sigma_x_n_cloud);
    auto width_below         = Kokkos::create_mirror_view(sigma_x_n_below);
    sigma_cloud((pdf_dim - 1)) = zero;
    sigma_below((pdf_dim - 1)) = zero;
    Kokkos::deep_copy(width_cloud, zero);
    Kokkos::deep_copy(width_below, zero);
    Kokkos::deep_copy(corr_array_cloud, zero);
    Kokkos::deep_copy(corr_array_below, zero);
    if (!l_const_Nc_in_cloud) {
      sigma_cloud(hm_metadata.iiPDF_Ncn) = hm_metadata.Ncnp2_on_Ncnm2;
      sigma_below(hm_metadata.iiPDF_Ncn) = hm_metadata.Ncnp2_on_Ncnm2;
    } else {
      sigma_cloud(hm_metadata.iiPDF_Ncn) = zero;
      sigma_below(hm_metadata.iiPDF_Ncn) = zero;
    }
    width_cloud(hm_metadata.iiPDF_Ncn) = stdev_L2N(sigma_cloud(hm_metadata.iiPDF_Ncn));
    width_below(hm_metadata.iiPDF_Ncn) = width_cloud(hm_metadata.iiPDF_Ncn);
    const auto hmp2_ip_on_hmm2_ip      = to_host_space(workspace, hm_metadata.hmp2_ip_on_hmm2_ip);
    for (int ivar = (hm_metadata.iiPDF_Ncn + 1); ivar < pdf_dim; ++ivar) {
      sigma_cloud(ivar) = clubb_params_host(iomicron) *
                          hmp2_ip_on_hmm2_ip((pdf2hydromet_idx(ivar, hm_metadata) + 1) - 1);
      sigma_below(ivar) = sigma_cloud(ivar);
      width_cloud(ivar) = stdev_L2N(sigma_cloud(ivar));
      width_below(ivar) = width_cloud(ivar);
    }
    Kokkos::deep_copy(sigma2_on_mu2_ip_cloud, sigma_cloud);
    Kokkos::deep_copy(sigma2_on_mu2_ip_below, sigma_below);
    Kokkos::deep_copy(sigma_x_n_cloud, width_cloud);
    Kokkos::deep_copy(sigma_x_n_below, width_below);
    // Native sequence association supplies one column and one vertical level.
    // The explicit C++ reshaping copies bits only and preserves that same shape.
    Real3D width_cloud_3d =
               workspace_scope.view<Real3D>("init_microphys::width_cloud_3d", 1, one_lev, pdf_dim),
           width_below_3d =
               workspace_scope.view<Real3D>("init_microphys::width_below_3d", 1, one_lev, pdf_dim);
    Real3D sigma_cloud_3d =
               workspace_scope.view<Real3D>("init_microphys::sigma_cloud_3d", 1, one_lev, pdf_dim),
           sigma_below_3d =
               workspace_scope.view<Real3D>("init_microphys::sigma_below_3d", 1, one_lev, pdf_dim);
    Real4D normal_cloud_4d = workspace_scope.view<Real4D>("init_microphys::normal_cloud_4d", 1,
                                                          one_lev, pdf_dim, pdf_dim),
           normal_below_4d = workspace_scope.view<Real4D>("init_microphys::normal_below_4d", 1,
                                                          one_lev, pdf_dim, pdf_dim);
    Real4D corr_cloud_4d = workspace_scope.view<Real4D>("init_microphys::corr_cloud_4d", 1, one_lev,
                                                        pdf_dim, pdf_dim),
           corr_below_4d = workspace_scope.view<Real4D>("init_microphys::corr_below_4d", 1, one_lev,
                                                        pdf_dim, pdf_dim);
    Kokkos::parallel_for(
        "init_microphys::sequence_association", pdf_dim, KOKKOS_LAMBDA(int j) {
          width_cloud_3d(0, 0, j) = sigma_x_n_cloud(j);
          width_below_3d(0, 0, j) = sigma_x_n_below(j);
          sigma_cloud_3d(0, 0, j) = sigma2_on_mu2_ip_cloud(j);
          sigma_below_3d(0, 0, j) = sigma2_on_mu2_ip_below(j);
          for (int k = 0; k < pdf_dim; ++k) {
            normal_cloud_4d(0, 0, j, k) = corr_array_n_cloud(j, k);
            normal_below_4d(0, 0, j, k) = corr_array_n_below(j, k);
          }
        });
    denorm_transform_corr(1, one_lev, pdf_dim, hm_metadata, width_cloud_3d, width_below_3d,
                          sigma_cloud_3d, sigma_below_3d, normal_cloud_4d, normal_below_4d,
                          corr_cloud_4d, corr_below_4d);
    Kokkos::parallel_for(
        "init_microphys::restore_correlations", (pdf_dim * pdf_dim), KOKKOS_LAMBDA(int cell) {
          const int j = (cell % pdf_dim), k = (cell / pdf_dim);
          corr_array_cloud(j, k) = corr_cloud_4d(0, 0, j, k);
          corr_array_below(j, k) = corr_below_4d(0, 0, j, k);
        });
    mirror_lower_triangular_matrix(pdf_dim, corr_array_cloud);
    mirror_lower_triangular_matrix(pdf_dim, corr_array_below);
    const auto cloud          = to_host_space(workspace, corr_array_cloud),
               below          = to_host_space(workspace, corr_array_below);
    const auto cout_flags     = std::cout.flags();
    const auto cout_precision = std::cout.precision();
    std::cout << " Correlation array (approximate); in cloud:\n";
    for (int ivar = 0; ivar < pdf_dim; ++ivar) {
      for (int j = 0; j < pdf_dim; ++j)
        std::cout << std::setw(7) << std::fixed << std::setprecision(3) << cloud(ivar, j);
      std::cout << '\n';
    }
    std::cout << " Correlation array (approximate); below cloud:\n";
    for (int ivar = 0; ivar < pdf_dim; ++ivar) {
      for (int j = 0; j < pdf_dim; ++j)
        std::cout << std::setw(7) << std::fixed << std::setprecision(3) << below(ivar, j);
      std::cout << '\n';
    }
    std::cout.flags(cout_flags);
    std::cout.precision(cout_precision);
    if (l_write_to_file) {
      case_log.open(case_info_file, std::ios::app);
      case_log << " Correlation array (approximate); in cloud:\n";
      for (int ivar = 0; ivar < pdf_dim; ++ivar) {
        for (int j = 0; j < pdf_dim; ++j)
          case_log << std::setw(7) << std::fixed << std::setprecision(3) << cloud(ivar, j);
        case_log << '\n';
      }
      case_log << " Correlation array (approximate); below cloud:\n";
      for (int ivar = 0; ivar < pdf_dim; ++ivar) {
        for (int j = 0; j < pdf_dim; ++j)
          case_log << std::setw(7) << std::fixed << std::setprecision(3) << below(ivar, j);
        case_log << '\n';
      }
      case_log.close();
    }
  }
}
// Preserve the Native allocation outputs when the caller owns host/strided Views.
// The appended records are the same case state used by the concrete API above.
template <class Params, class Cloud, class Below>
inline void
init_microphys(ScratchWorkspace& workspace, int iunit, const std::string& runtype,
               const std::string& namelist_file, const std::string& case_info_file, real_t host_dx,
               real_t host_dy, Params clubb_params, bool l_diagnose_correlations,
               bool& l_const_Nc_in_cloud, bool& l_fix_w_chi_eta_correlations, int& hydromet_dim,
               int& pdf_dim, hm_metadata_type& hm_metadata,
               silhs_config_flags_type& silhs_config_flags, real_t& vert_decorr_coef_out,
               Cloud& corr_array_n_cloud, Below& corr_array_n_below,
               parameters_microphys& microphys, parameters_KK& kk,
               morrison_microphysics_state& morrison, parameters_silhs_state& silhs) {
  auto workspace_scope = workspace.scope();
  auto params          = to_execution_space(workspace, clubb_params);
  Real2D cloud, below;
  init_microphys(workspace, iunit, runtype, namelist_file, case_info_file, host_dx, host_dy, params,
                 l_diagnose_correlations, l_const_Nc_in_cloud, l_fix_w_chi_eta_correlations,
                 hydromet_dim, pdf_dim, hm_metadata, silhs_config_flags, vert_decorr_coef_out,
                 cloud, below, microphys, kk, morrison, silhs);
  corr_array_n_cloud = Cloud("init_microphys::corr_array_n_cloud", pdf_dim, pdf_dim);
  corr_array_n_below = Below("init_microphys::corr_array_n_below", pdf_dim, pdf_dim);
  from_execution_space(workspace, corr_array_n_cloud, cloud);
  from_execution_space(workspace, corr_array_n_below, below);
}
inline void cleanup_microphys(parameters_microphys& microphys) {
  // Deallocate arrays used by microphysics. The appended record replaces
  // Native parameters_microphys module storage and leaves unrelated state intact.
  if (microphys.l_hydromet_sed.is_allocated())
    microphys.l_hydromet_sed = Bool1D{};
}
} // namespace clubb_kokkos
