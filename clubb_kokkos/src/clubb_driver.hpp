// Declare the case-owned standalone driver ported from src/clubb_driver.F90.
//
// Construct a driver after Kokkos initialization. run_clubb executes a complete
// case. For controlled advances, call init_clubb_case, set_case_initial_conditions,
// advance_clubb_to_end and clean_up_clubb in lifecycle order. Initialization
// prepares the first batch; resets also support reruns and later batches.
// get_stats returns a read-only snapshot, get_runtime_batch_config returns the
// requested/runtime column counts, and get_clubb_params_all copies all parameters.
//
// Initial conditions use column-first host Views. Case-owned execution Views
// retain runtime fields across timesteps and bounded advances. Statistics and
// output run in C++; only optional BFB PDF matrix operations retain a Fortran
// boundary. Mutable state belongs to a case.
#pragma once
#include "CLUBB_core/error_code.hpp"
#include "CLUBB_core/corr_varnce_module.hpp"
#include "CLUBB_core/hydromet_pdf_parameter_module.hpp"
#include "CLUBB_core/model_flags.hpp"
#include "CLUBB_core/parameters_tunable.hpp"
#include "CLUBB_core/pdf_parameter_module.hpp"
#include "CLUBB_core/stats_netcdf.hpp"
#include "driver_io.hpp"
#include "driver_data.hpp"
#include "scratch_workspace.hpp"
#include "CLUBB_core/sponge_layer_damping.hpp"
#include "Benchmark_cases/time_dependent_input.hpp"
#include "Benchmark_cases/mpace_a.hpp"
#include "Input_fields/inputfields_state.hpp"
#include "readiopdata_module.hpp"
#include "Radiation/parameters_radiation.hpp"
#include "Radiation/extended_atmosphere_module.hpp"
#include "Radiation/radiation_variables_module.hpp"
#include "Radiation/simple_rad_module.hpp"
#include "Radiation/soil_vegetation.hpp"

#include "Microphys/parameters_microphys.hpp"
#include "Microphys/KK_microphys/parameters_KK.hpp"
#include "Microphys/Morrison_microphys/module_mp_graupel_state.hpp"
#include "SILHS/parameters_silhs.hpp"
#include "SILHS/latin_hypercube_arrays.hpp"
#include "SILHS/mt95.hpp"
#include <fstream>
#include <optional>
#include <string>

namespace clubb_kokkos {
using DriverArray1 = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using DriverSlice1 = Kokkos::View<real_t*, Kokkos::LayoutStride, Kokkos::HostSpace>;
using DriverBool1  = Kokkos::View<std::uint8_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using DriverArray2 = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using DriverSlice2 = Kokkos::View<real_t**, Kokkos::LayoutStride, Kokkos::HostSpace>;
using DriverArray3 = Kokkos::View<real_t***, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using DriverSlice3 = Kokkos::View<real_t***, Kokkos::LayoutStride, Kokkos::HostSpace>;
using DriverInt3   = Kokkos::View<int***, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using DriverArray4 = Kokkos::View<real_t****, Kokkos::LayoutLeft, Kokkos::HostSpace>;

// Host initialization and execution-space timestepping share one case lifetime.
// Host physics and sampled statistics cross explicit storage boundaries.
class ClubbDriver {
public:
  ClubbDriver();
  ClubbDriver(const ClubbDriver&)            = delete;
  ClubbDriver& operator=(const ClubbDriver&) = delete;
  void run_clubb(const std::string& namelist_filename,                 // In
                 bool l_stdout,                                        // In
                 err_info_type& err_info,                              // Out
                 const DriverArray2* clubb_params_in  = nullptr,       // Optional In
                 const DriverBool1* model_flags_array = nullptr);      // Optional In
  void get_stats(stats_type& stats_out) const;                         // Out
  void get_runtime_batch_config(int& total_param_sets_out,             // Out
                                int& batch_size_out) const;            // Out
  void get_clubb_params_all(DriverArray2& clubb_params_all_out) const; // Out

  void init_clubb_case(const std::string& runfile,                                // In
                       err_info_type& err_info,                                   // Out
                       const DriverArray2* clubb_params_in  = nullptr,            // Optional In
                       const DriverBool1* model_flags_array = nullptr);           // Optional In
  void set_case_initial_conditions(err_info_type& err_info,                       // InOut
                                   const DriverArray2* clubb_params_in = nullptr, // Optional In
                                   std::optional<int> batch_num = std::nullopt);  // Optional In
  void advance_clubb_to_end(bool l_stdout,                                        // In
                            err_info_type& err_info,                              // InOut
                            bool l_suppress_stats          = false,               // Optional In
                            std::optional<int> itime_start = std::nullopt,        // Optional In
                            std::optional<int> itime_end   = std::nullopt);       // Optional In
  void clean_up_clubb(err_info_type& err_info);                                   // InOut

private:
  ScratchWorkspace workspace; // Case-owned temporary storage, reused by timestep scopes.

  void initialize_clubb(
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
      DriverSlice3 sclrm, DriverSlice3 edsclrm);                                     // Out
  void initialize_clubb_variables(int ngrdcol, int sclr_dim, int edsclr_dim,
                                  const sclr_idx_type& sclr_idx,                    // In
                                  const grid& gr, const std::string& altitude_type, // In
                                  const std::string& temperature_type,              // In
                                  DriverSlice1 p_sfc, DriverSlice1 rtm_sfc, DriverSlice2 rtm,
                                  bool l_readiopdata,                      // In
                                  int saturation_formula,                  // In
                                  DriverSlice2 thlm, DriverSlice2 p_in_Pa, // InOut
                                  DriverSlice2 p_in_Pa_zm,                 // Out
                                  DriverSlice2 exner, DriverSlice2 rho, DriverSlice2 rho_zm, // Out
                                  DriverSlice2 rcm, DriverSlice2 thvm,
                                  DriverSlice2 rho_ds_zm,                               // Out
                                  DriverSlice2 rho_ds_zt, DriverSlice2 invrs_rho_ds_zm, // Out
                                  DriverSlice2 invrs_rho_ds_zt, DriverSlice2 thv_ds_zm, // Out
                                  DriverSlice2 thv_ds_zt,                               // Out
                                  DriverSlice3 sclrm, DriverSlice3 edsclrm);            // InOut
  void restart_clubb(const grid& gr,                                                    // In
                     int hydromet_dim,                                                  // In
                     const hm_metadata_type& hm_metadata,                               // In
                     const std::string& restart_path_case,                              // In
                     double time_restart,                                               // In
                     DriverSlice1 um,                                                   // InOut
                     DriverSlice1 upwp,                                                 // InOut
                     DriverSlice1 vm,                                                   // InOut
                     DriverSlice1 vpwp,                                                 // InOut
                     DriverSlice1 up2,                                                  // InOut
                     DriverSlice1 vp2,                                                  // InOut
                     DriverSlice1 rtm,                                                  // InOut
                     DriverSlice1 wprtp,                                                // InOut
                     DriverSlice1 thlm,                                                 // InOut
                     DriverSlice1 wpthlp,                                               // InOut
                     DriverSlice1 rtp2,                                                 // InOut
                     DriverSlice1 rtp3,                                                 // InOut
                     DriverSlice1 thlp2,                                                // InOut
                     DriverSlice1 thlp3,                                                // InOut
                     DriverSlice1 rtpthlp,                                              // InOut
                     DriverSlice1 wp2,                                                  // InOut
                     DriverSlice1 wp3,                                                  // InOut
                     DriverSlice1 p_in_Pa,                                              // InOut
                     DriverSlice1 exner,                                                // InOut
                     DriverSlice1 rcm,                                                  // InOut
                     DriverSlice1 cloud_frac,                                           // InOut
                     DriverSlice1 wpthvp,                                               // InOut
                     DriverSlice1 wp2thvp,                                              // InOut
                     DriverSlice1 wp2up,                                                // InOut
                     DriverSlice1 rtpthvp,                                              // InOut
                     DriverSlice1 thlpthvp,                                             // InOut
                     DriverSlice1 wp2rtp,                                               // InOut
                     DriverSlice1 wp2thlp,                                              // InOut
                     DriverSlice1 uprcp,                                                // InOut
                     DriverSlice1 vprcp,                                                // InOut
                     DriverSlice1 rc_coef_zm,                                           // InOut
                     DriverSlice1 wp4,                                                  // InOut
                     DriverSlice1 wpup2,                                                // InOut
                     DriverSlice1 wpvp2,                                                // InOut
                     DriverSlice1 wp2up2,                                               // InOut
                     DriverSlice1 wp2vp2,                                               // InOut
                     DriverSlice1 ice_supersat_frac,                                    // InOut
                     DriverSlice1 wm_zt,                                                // InOut
                     DriverSlice1 rho,                                                  // InOut
                     DriverSlice1 rho_zm,                                               // InOut
                     DriverSlice1 rho_ds_zm,                                            // InOut
                     DriverSlice1 rho_ds_zt,                                            // InOut
                     DriverSlice1 thv_ds_zm,                                            // InOut
                     DriverSlice1 thv_ds_zt,                                            // InOut
                     DriverSlice1 thlm_forcing,                                         // InOut
                     DriverSlice1 rtm_forcing,                                          // InOut
                     DriverSlice1 wprtp_forcing,                                        // InOut
                     DriverSlice1 wpthlp_forcing,                                       // InOut
                     DriverSlice1 rtp2_forcing,                                         // InOut
                     DriverSlice1 thlp2_forcing,                                        // InOut
                     DriverSlice1 rtpthlp_forcing,                                      // InOut
                     DriverSlice2 hydromet,                                             // InOut
                     DriverSlice2 hydrometp2,                                           // InOut
                     DriverSlice2 wphydrometp,                                          // InOut
                     DriverSlice1 Ncm,                                                  // InOut
                     DriverSlice1 Nccnm,                                                // InOut
                     DriverSlice1 thvm,                                                 // InOut
                     DriverSlice1 em,                                                   // InOut
                     DriverSlice1 tau_zm,                                               // InOut
                     DriverSlice1 tau_zt,                                               // InOut
                     DriverSlice1 Kh_zt,                                                // InOut
                     DriverSlice1 Kh_zm,                                                // InOut
                     DriverSlice1 ug,                                                   // InOut
                     DriverSlice1 vg,                                                   // InOut
                     DriverSlice1 thlprcp,                                              // InOut
                     DriverSlice1 sigma_sqd_w,                                          // InOut
                     DriverSlice1 sigma_sqd_w_zt,                                       // InOut
                     DriverSlice1 radht,                                                // InOut
                     real_t& deep_soil_T_in_K,                                          // InOut
                     real_t& sfc_soil_T_in_K,                                           // InOut
                     real_t& veg_T_in_K,                                                // InOut
                     pdf_parameter& pdf_params,                                         // InOut
                     pdf_parameter& pdf_params_zm,                                      // InOut
                     DriverSlice1 rcm_mc,                                               // Out
                     DriverSlice1 rvm_mc,                                               // Out
                     DriverSlice1 thlm_mc,                                              // Out
                     DriverSlice1 wprtp_mc,                                             // Out
                     DriverSlice1 wpthlp_mc,                                            // Out
                     DriverSlice1 rtp2_mc,                                              // Out
                     DriverSlice1 thlp2_mc,                                             // Out
                     DriverSlice1 rtpthlp_mc,                                           // Out
                     real_t& wpthlp_sfc,                                                // Out
                     real_t& wprtp_sfc,                                                 // Out
                     real_t& upwp_sfc,                                                  // Out
                     real_t& vpwp_sfc);                                                 // Out
  // TODO: Native enter/update/exit data regions need explicit typed storage in
  // C++. Keep these host/execution pairs until initialization also runs on device.
  struct ExecutionFields {
    HostDeviceField<DriverArray2> clubb_params;
    HostDeviceField<DriverArray1> dummy_dx;
    HostDeviceField<DriverArray1> dummy_dy;
    HostDeviceField<DriverArray1> sclr_tol;
    HostDeviceField<DriverArray1> p_sfc;
    HostDeviceField<DriverArray1> T_sfc;
    HostDeviceField<DriverArray1> fcor;
    HostDeviceField<DriverArray1> fcor_y;
    HostDeviceField<DriverArray1> sfc_elevation;
    HostDeviceField<DriverArray2> corr_array_n_cloud;
    HostDeviceField<DriverArray2> corr_array_n_below;
    HostDeviceField<DriverInt3> X_mixt_comp_all_levs;
    HostDeviceField<DriverArray1> lmin;
    HostDeviceField<DriverArray1> mixt_frac_max_mag;
    HostDeviceField<DriverArray1> deep_soil_T_in_K;
    HostDeviceField<DriverArray1> sfc_soil_T_in_K;
    HostDeviceField<DriverArray1> veg_T_in_K;
    HostDeviceField<DriverArray1> wpthlp_sfc;
    HostDeviceField<DriverArray1> wprtp_sfc;
    HostDeviceField<DriverArray1> upwp_sfc;
    HostDeviceField<DriverArray1> vpwp_sfc;
    HostDeviceField<DriverArray1> upwp_sfc_pert;
    HostDeviceField<DriverArray1> vpwp_sfc_pert;
    HostDeviceField<DriverArray2> um;
    HostDeviceField<DriverArray2> upwp;
    HostDeviceField<DriverArray2> vm;
    HostDeviceField<DriverArray2> vpwp;
    HostDeviceField<DriverArray2> up2;
    HostDeviceField<DriverArray2> vp2;
    HostDeviceField<DriverArray2> up3;
    HostDeviceField<DriverArray2> vp3;
    HostDeviceField<DriverArray2> wprtp;
    HostDeviceField<DriverArray2> wpthlp;
    HostDeviceField<DriverArray2> rtp2;
    HostDeviceField<DriverArray2> rtp3;
    HostDeviceField<DriverArray2> thlp2;
    HostDeviceField<DriverArray2> thlp3;
    HostDeviceField<DriverArray2> rtpthlp;
    HostDeviceField<DriverArray2> wp2;
    HostDeviceField<DriverArray2> thvm;
    HostDeviceField<DriverArray2> exner;
    HostDeviceField<DriverArray2> rtm;
    HostDeviceField<DriverArray2> thlm;
    HostDeviceField<DriverArray2> rcm;
    HostDeviceField<DriverArray2> wp3;
    HostDeviceField<DriverArray2> p_in_Pa;
    HostDeviceField<DriverArray2> cloud_frac;
    HostDeviceField<DriverArray2> wpthvp;
    HostDeviceField<DriverArray2> wp2thvp;
    HostDeviceField<DriverArray2> wp2up;
    HostDeviceField<DriverArray2> rtpthvp;
    HostDeviceField<DriverArray2> thlpthvp;
    HostDeviceField<DriverArray2> wp2rtp;
    HostDeviceField<DriverArray2> wp2thlp;
    HostDeviceField<DriverArray2> uprcp;
    HostDeviceField<DriverArray2> vprcp;
    HostDeviceField<DriverArray2> rc_coef_zm;
    HostDeviceField<DriverArray2> wp4;
    HostDeviceField<DriverArray2> wpup2;
    HostDeviceField<DriverArray2> wpvp2;
    HostDeviceField<DriverArray2> wp2up2;
    HostDeviceField<DriverArray2> wp2vp2;
    HostDeviceField<DriverArray2> rho_ds_zt;
    HostDeviceField<DriverArray2> thv_ds_zt;
    HostDeviceField<DriverArray2> wm_zm;
    HostDeviceField<DriverArray2> wm_zt;
    HostDeviceField<DriverArray2> rho;
    HostDeviceField<DriverArray2> rho_zm;
    HostDeviceField<DriverArray2> rho_ds_zm;
    HostDeviceField<DriverArray2> invrs_rho_ds_zm;
    HostDeviceField<DriverArray2> invrs_rho_ds_zt;
    HostDeviceField<DriverArray2> thv_ds_zm;
    HostDeviceField<DriverArray2> thlm_forcing;
    HostDeviceField<DriverArray2> rtm_forcing;
    HostDeviceField<DriverArray2> um_forcing;
    HostDeviceField<DriverArray2> vm_forcing;
    HostDeviceField<DriverArray2> wprtp_forcing;
    HostDeviceField<DriverArray2> wpthlp_forcing;
    HostDeviceField<DriverArray2> rtp2_forcing;
    HostDeviceField<DriverArray2> thlp2_forcing;
    HostDeviceField<DriverArray2> rtpthlp_forcing;
    HostDeviceField<DriverArray2> um_pert;
    HostDeviceField<DriverArray2> vm_pert;
    HostDeviceField<DriverArray2> upwp_pert;
    HostDeviceField<DriverArray2> vpwp_pert;
    HostDeviceField<DriverArray2> Skw_zm_smooth;
    HostDeviceField<DriverArray2> Lscale;
    HostDeviceField<DriverArray2> thlprcp;
    HostDeviceField<DriverArray2> wprcp;
    HostDeviceField<DriverArray2> w_up_in_cloud;
    HostDeviceField<DriverArray2> w_down_in_cloud;
    HostDeviceField<DriverArray2> cloudy_updraft_frac;
    HostDeviceField<DriverArray2> cloudy_downdraft_frac;
    HostDeviceField<DriverArray2> ice_supersat_frac;
    HostDeviceField<DriverArray2> rcm_in_layer;
    HostDeviceField<DriverArray2> cloud_cover;
    HostDeviceField<DriverArray2> invrs_tau_zm;
    HostDeviceField<DriverArray2> ug;
    HostDeviceField<DriverArray2> vg;
    HostDeviceField<DriverArray2> rtm_ref;
    HostDeviceField<DriverArray2> thlm_ref;
    HostDeviceField<DriverArray2> um_ref;
    HostDeviceField<DriverArray2> vm_ref;
    HostDeviceField<DriverArray2> rcm_mc;
    HostDeviceField<DriverArray2> rvm_mc;
    HostDeviceField<DriverArray2> thlm_mc;
    HostDeviceField<DriverArray2> wprtp_mc;
    HostDeviceField<DriverArray2> wpthlp_mc;
    HostDeviceField<DriverArray2> rtp2_mc;
    HostDeviceField<DriverArray2> thlp2_mc;
    HostDeviceField<DriverArray2> rtpthlp_mc;
    HostDeviceField<DriverArray2> Ncm_mc;
    HostDeviceField<DriverArray2> Ncm;
    HostDeviceField<DriverArray2> Nccnm;
    HostDeviceField<DriverArray2> Kh_zt;
    HostDeviceField<DriverArray2> Kh_zm;
    HostDeviceField<DriverArray2> rfrzm;
    HostDeviceField<DriverArray2> Nc_in_cloud;
    HostDeviceField<DriverArray2> wpNcp;
    HostDeviceField<DriverArray2> radht;
    HostDeviceField<DriverArray2> wpsclrp_sfc;
    HostDeviceField<DriverArray2> wpedsclrp_sfc;
    HostDeviceField<DriverArray3> sclrm;
    HostDeviceField<DriverArray3> wpsclrp;
    HostDeviceField<DriverArray3> sclrp2;
    HostDeviceField<DriverArray3> sclrp3;
    HostDeviceField<DriverArray3> sclrprtp;
    HostDeviceField<DriverArray3> sclrpthlp;
    HostDeviceField<DriverArray3> sclrpthvp;
    HostDeviceField<DriverArray3> sclrm_forcing;
    HostDeviceField<DriverArray3> edsclrm;
    HostDeviceField<DriverArray3> edsclrm_forcing;
    HostDeviceField<DriverArray3> lh_rt_clipped;
    HostDeviceField<DriverArray3> lh_thl_clipped;
    HostDeviceField<DriverArray3> lh_rc_clipped;
    HostDeviceField<DriverArray3> lh_rv_clipped;
    HostDeviceField<DriverArray3> lh_Nc_clipped;
    HostDeviceField<DriverArray3> lh_sample_point_weights;
    HostDeviceField<DriverArray3> K_hm;
    HostDeviceField<DriverArray3> wp2hmp;
    HostDeviceField<DriverArray3> rtphmp_zt;
    HostDeviceField<DriverArray3> thlphmp_zt;
    HostDeviceField<DriverArray3> hydromet;
    HostDeviceField<DriverArray3> hydrometp2;
    HostDeviceField<DriverArray3> wphydrometp;
    HostDeviceField<DriverArray3> hydromet_mc;
    HostDeviceField<DriverArray3> hydromet_vel_zt;
    HostDeviceField<DriverArray3> hydromet_vel_covar_zt_impc;
    HostDeviceField<DriverArray3> hydromet_vel_covar_zt_expc;
    HostDeviceField<DriverArray3> mu_x_1_n;
    HostDeviceField<DriverArray3> mu_x_2_n;
    HostDeviceField<DriverArray3> sigma_x_1_n;
    HostDeviceField<DriverArray3> sigma_x_2_n;
    HostDeviceField<DriverArray4> corr_array_1_n;
    HostDeviceField<DriverArray4> corr_array_2_n;
    HostDeviceField<DriverArray4> corr_cholesky_mtx_1;
    HostDeviceField<DriverArray4> corr_cholesky_mtx_2;
    HostDeviceField<DriverArray4> X_nl_all_levs;
  } execution;
  bool execution_fields_initialized{};

  template <class Function> void for_each_execution_field(Function&& function) {
    // Apply lifecycle transfers to the matching host/execution allocations.
    function(execution.clubb_params, clubb_params);
    function(execution.dummy_dx, dummy_dx);
    function(execution.dummy_dy, dummy_dy);
    function(execution.sclr_tol, sclr_tol);
    function(execution.p_sfc, p_sfc);
    function(execution.T_sfc, T_sfc);
    function(execution.fcor, fcor);
    function(execution.fcor_y, fcor_y);
    function(execution.sfc_elevation, sfc_elevation);
    function(execution.corr_array_n_cloud, corr_array_n_cloud);
    function(execution.corr_array_n_below, corr_array_n_below);
    function(execution.X_mixt_comp_all_levs, X_mixt_comp_all_levs);
    function(execution.lmin, lmin);
    function(execution.mixt_frac_max_mag, mixt_frac_max_mag);
    function(execution.deep_soil_T_in_K, deep_soil_T_in_K);
    function(execution.sfc_soil_T_in_K, sfc_soil_T_in_K);
    function(execution.veg_T_in_K, veg_T_in_K);
    function(execution.wpthlp_sfc, wpthlp_sfc);
    function(execution.wprtp_sfc, wprtp_sfc);
    function(execution.upwp_sfc, upwp_sfc);
    function(execution.vpwp_sfc, vpwp_sfc);
    function(execution.upwp_sfc_pert, upwp_sfc_pert);
    function(execution.vpwp_sfc_pert, vpwp_sfc_pert);
    function(execution.um, um);
    function(execution.upwp, upwp);
    function(execution.vm, vm);
    function(execution.vpwp, vpwp);
    function(execution.up2, up2);
    function(execution.vp2, vp2);
    function(execution.up3, up3);
    function(execution.vp3, vp3);
    function(execution.wprtp, wprtp);
    function(execution.wpthlp, wpthlp);
    function(execution.rtp2, rtp2);
    function(execution.rtp3, rtp3);
    function(execution.thlp2, thlp2);
    function(execution.thlp3, thlp3);
    function(execution.rtpthlp, rtpthlp);
    function(execution.wp2, wp2);
    function(execution.thvm, thvm);
    function(execution.exner, exner);
    function(execution.rtm, rtm);
    function(execution.thlm, thlm);
    function(execution.rcm, rcm);
    function(execution.wp3, wp3);
    function(execution.p_in_Pa, p_in_Pa);
    function(execution.cloud_frac, cloud_frac);
    function(execution.wpthvp, wpthvp);
    function(execution.wp2thvp, wp2thvp);
    function(execution.wp2up, wp2up);
    function(execution.rtpthvp, rtpthvp);
    function(execution.thlpthvp, thlpthvp);
    function(execution.wp2rtp, wp2rtp);
    function(execution.wp2thlp, wp2thlp);
    function(execution.uprcp, uprcp);
    function(execution.vprcp, vprcp);
    function(execution.rc_coef_zm, rc_coef_zm);
    function(execution.wp4, wp4);
    function(execution.wpup2, wpup2);
    function(execution.wpvp2, wpvp2);
    function(execution.wp2up2, wp2up2);
    function(execution.wp2vp2, wp2vp2);
    function(execution.rho_ds_zt, rho_ds_zt);
    function(execution.thv_ds_zt, thv_ds_zt);
    function(execution.wm_zm, wm_zm);
    function(execution.wm_zt, wm_zt);
    function(execution.rho, rho);
    function(execution.rho_zm, rho_zm);
    function(execution.rho_ds_zm, rho_ds_zm);
    function(execution.invrs_rho_ds_zm, invrs_rho_ds_zm);
    function(execution.invrs_rho_ds_zt, invrs_rho_ds_zt);
    function(execution.thv_ds_zm, thv_ds_zm);
    function(execution.thlm_forcing, thlm_forcing);
    function(execution.rtm_forcing, rtm_forcing);
    function(execution.um_forcing, um_forcing);
    function(execution.vm_forcing, vm_forcing);
    function(execution.wprtp_forcing, wprtp_forcing);
    function(execution.wpthlp_forcing, wpthlp_forcing);
    function(execution.rtp2_forcing, rtp2_forcing);
    function(execution.thlp2_forcing, thlp2_forcing);
    function(execution.rtpthlp_forcing, rtpthlp_forcing);
    function(execution.um_pert, um_pert);
    function(execution.vm_pert, vm_pert);
    function(execution.upwp_pert, upwp_pert);
    function(execution.vpwp_pert, vpwp_pert);
    function(execution.Skw_zm_smooth, Skw_zm_smooth);
    function(execution.Lscale, Lscale);
    function(execution.thlprcp, thlprcp);
    function(execution.wprcp, wprcp);
    function(execution.w_up_in_cloud, w_up_in_cloud);
    function(execution.w_down_in_cloud, w_down_in_cloud);
    function(execution.cloudy_updraft_frac, cloudy_updraft_frac);
    function(execution.cloudy_downdraft_frac, cloudy_downdraft_frac);
    function(execution.ice_supersat_frac, ice_supersat_frac);
    function(execution.rcm_in_layer, rcm_in_layer);
    function(execution.cloud_cover, cloud_cover);
    function(execution.invrs_tau_zm, invrs_tau_zm);
    function(execution.ug, ug);
    function(execution.vg, vg);
    function(execution.rtm_ref, rtm_ref);
    function(execution.thlm_ref, thlm_ref);
    function(execution.um_ref, um_ref);
    function(execution.vm_ref, vm_ref);
    function(execution.rcm_mc, rcm_mc);
    function(execution.rvm_mc, rvm_mc);
    function(execution.thlm_mc, thlm_mc);
    function(execution.wprtp_mc, wprtp_mc);
    function(execution.wpthlp_mc, wpthlp_mc);
    function(execution.rtp2_mc, rtp2_mc);
    function(execution.thlp2_mc, thlp2_mc);
    function(execution.rtpthlp_mc, rtpthlp_mc);
    function(execution.Ncm_mc, Ncm_mc);
    function(execution.Ncm, Ncm);
    function(execution.Nccnm, Nccnm);
    function(execution.Kh_zt, Kh_zt);
    function(execution.Kh_zm, Kh_zm);
    function(execution.rfrzm, rfrzm);
    function(execution.Nc_in_cloud, Nc_in_cloud);
    function(execution.wpNcp, wpNcp);
    function(execution.radht, radht);
    function(execution.wpsclrp_sfc, wpsclrp_sfc);
    function(execution.wpedsclrp_sfc, wpedsclrp_sfc);
    function(execution.sclrm, sclrm);
    function(execution.wpsclrp, wpsclrp);
    function(execution.sclrp2, sclrp2);
    function(execution.sclrp3, sclrp3);
    function(execution.sclrprtp, sclrprtp);
    function(execution.sclrpthlp, sclrpthlp);
    function(execution.sclrpthvp, sclrpthvp);
    function(execution.sclrm_forcing, sclrm_forcing);
    function(execution.edsclrm, edsclrm);
    function(execution.edsclrm_forcing, edsclrm_forcing);
    function(execution.lh_rt_clipped, lh_rt_clipped);
    function(execution.lh_thl_clipped, lh_thl_clipped);
    function(execution.lh_rc_clipped, lh_rc_clipped);
    function(execution.lh_rv_clipped, lh_rv_clipped);
    function(execution.lh_Nc_clipped, lh_Nc_clipped);
    function(execution.lh_sample_point_weights, lh_sample_point_weights);
    function(execution.K_hm, K_hm);
    function(execution.wp2hmp, wp2hmp);
    function(execution.rtphmp_zt, rtphmp_zt);
    function(execution.thlphmp_zt, thlphmp_zt);
    function(execution.hydromet, hydromet);
    function(execution.hydrometp2, hydrometp2);
    function(execution.wphydrometp, wphydrometp);
    function(execution.hydromet_mc, hydromet_mc);
    function(execution.hydromet_vel_zt, hydromet_vel_zt);
    function(execution.hydromet_vel_covar_zt_impc, hydromet_vel_covar_zt_impc);
    function(execution.hydromet_vel_covar_zt_expc, hydromet_vel_covar_zt_expc);
    function(execution.mu_x_1_n, mu_x_1_n);
    function(execution.mu_x_2_n, mu_x_2_n);
    function(execution.sigma_x_1_n, sigma_x_1_n);
    function(execution.sigma_x_2_n, sigma_x_2_n);
    function(execution.corr_array_1_n, corr_array_1_n);
    function(execution.corr_array_2_n, corr_array_2_n);
    function(execution.corr_cholesky_mtx_1, corr_cholesky_mtx_1);
    function(execution.corr_cholesky_mtx_2, corr_cholesky_mtx_2);
    function(execution.X_nl_all_levs, X_nl_all_levs);
  }

  // Each former Native module stores mutable state in this driver case.
  time_dependent_input_state time_dependent_input;
  mpace_a_state mpace;
  inputfields_state inputfields;
  readiopdata_state iop;
  parameters_radiation_state parameters_radiation;
  extended_atmosphere_state extended_atmosphere;
  radiation_variables_state radiation_variables;
  simple_rad_state simple_radiation;
  soil_vegetation_state soil_vegetation;
  sponge_damp_settings thlm_sponge_damp_settings;
  sponge_damp_settings rtm_sponge_damp_settings;
  sponge_damp_settings uv_sponge_damp_settings;
  sponge_damp_settings wp2_sponge_damp_settings;
  sponge_damp_settings wp3_sponge_damp_settings;
  sponge_damp_settings up2_vp2_sponge_damp_settings;
  parameters_microphys microphys;
  parameters_KK kk;
  morrison_microphysics_state morrison;
  parameters_silhs_state silhs;
  latin_hypercube_arrays_state sampling_state;
  mt95 rng;
  std::ofstream case_log;
  static constexpr bool l_implemented{false};   // CLUBB is not being run inside of a host model
  static constexpr bool l_ascending_grid{true}; // CLUBB runs with an ascending grid by default
  static constexpr bool l_write_to_file{true};  // If true, will write case information to a file
  static constexpr int nlon{1};                 // Number of points in the X/Y [-]
  static constexpr int nlat{1};
  bool l_restart{};                  // Flag for restarting from netCDF file
  bool l_input_fields{};             // Whether to set model variables from a file
  bool l_readiopdata{};              // Whether to read in IOP data from a .nc file
  bool l_modify_ic_with_cubic_int{}; // Use Steffen monotone cubic interpolation for sounding
                                     // profiles.
  bool l_modify_bc_for_cnvg_test{}; // Compute surface fluxes at a fixed 25 m height for convergence
                                    // tests.
  bool l_stats{};                   // Whether statistics are computed and output to disk
  bool l_silhs_out{};               // Whether to output SILHS files
  bool l_rad_itime{};               // .true. if we calculate radiation at the current time step
  bool l_first_write_to_file{};     // True until the first grid-adaptation output write.
  bool l_allow_small_stats_tout{};
  bool l_restart_input{};
  bool l_last_timestep{};
  int ifinal{};                // Total number of iterations (final iter value)
  int stats_nsamp{};           // Stats sampling interval [timestep]
  int stats_nout{};            // Stats output interval   [timestep]
  int ngrdcol{};               // Number of columns
  int iunit{};                 // File unit used for I/O
  int iunit_grid_adaptation{}; // Native file unit for grid-adaptation output.
  int hydromet_dim{};          // Number of hydrometeor species        [#]
  int sclr_dim{};              // Number of passive scalars            [#]
  int edsclr_dim{};            // Number of passive scalars            [#]
  int pdf_dim{};
  int nzmax{};     // Vertical extent in levels( relevant for grid type 2 and 3 only )  [#]
  int grid_type{}; // 1: uniform; 2: prescribed zt heights; 3: prescribed zm heights.
  int day{};
  int month{};
  int year{};
  int batch_size{1}; // Runtime column count used for one batch advance.
  int sfctype{};     // 0: prescribed heat fluxes; 1: bulk ocean flux formula.
  stats_type stats{};
  DriverArray2 clubb_params{};     // Active runtime batch of CLUBB parameters.
  DriverArray2 clubb_params_all{}; // Full parameter matrix across all requested sets.
  int total_param_sets{1};         // Total number of requested parameter columns.
  sclr_idx_type sclr_idx{};
  hm_metadata_type hm_metadata{};
  clubb_config_flags_type clubb_config_flags{}; // Derived type holding all configurable CLUBB flags
  silhs_config_flags_type silhs_config_flags{}; // Flags for the SILHS sampling code
  real_t stats_tsamp{};                         // Stats sampling interval [s]
  real_t stats_tout{};                          // Stats output interval   [s]
  double stats_tstart{};                        // Stats sampling window start time [s]
  double stats_tend{};                          // Stats sampling window end time   [s]
  double time_restart{};                        // Time of model restart run     [s]
  std::string restart_path_case{};              // netCDF file used in case of restart
  std::string forcings_file_path{};             // Path to the forcing files
  std::string fname_prefix{};   // Prefix of stats filenames, to be followed by, for example "_zt"
  std::string output_dir{};     // Output directory for setup/log/stats files
  std::string case_info_file{}; // The filename for case info
  std::string fname_grid_adaptation{}; // Text output of changing grid heights.
  std::string output_file_prefix{};
  std::string runtype{};           // String identifying the model case; e.g. bomex
  std::string zt_grid_fname{};     // Path and filename of thermodynamic level altitudes
  std::string zm_grid_fname{};     // Path and filename of momentum level altitudes
  real_t lat_vals{};               // Latitude  [Degrees North]
  real_t lon_vals{};               // Longitude [Degrees East]
  double time_initial{};           // Time of start of simulation     [s]
  double time_final{};             // Time end of simulation          [s]
  double time_current{};           // Current time of simulation      [s]
  real_t dt_main{};                // Main model timestep                    [s]
  real_t dt_rad{};                 // Closure model timestep                 [s]
  real_t deltaz_nl{};              // Change in altitude per grid level     [m]
  real_t zm_init_nl{};             // Initial point on the momentum grid    [m]
  real_t zm_top_nl{};              // Maximum point on the momentum grid    [m]
  real_t T0{};                     // Reference temperature (usually 300)  [K]
  real_t ts_nudge{};               // Timescale of u/v nudging             [s]
  real_t rtm_min{};                // Value below which rtm will be nudged [kg/kg]
  real_t rtm_nudge_max_altitude{}; // Highest altitude at which to nudge rtm [m]
  DriverArray1 deep_soil_T_in_K_init{};
  DriverArray1 sfc_soil_T_in_K_init{};
  DriverArray1 veg_T_in_K_init{};
  DriverArray1 dummy_dx{};
  DriverArray1 dummy_dy{}; // [m]
  DriverArray1 sclr_tol{}; // Thresholds on the passive scalars     [units vary]
  DriverArray1 thlm_init{};
  DriverArray1 rtm_init{};
  DriverArray1 um_init{};
  DriverArray1 vm_init{};
  DriverArray1 ug_init{};
  DriverArray1 vg_init{};
  DriverArray1 wp2_init{};
  DriverArray1 up2_init{};
  DriverArray1 vp2_init{};
  DriverArray1 upwp_init{};
  DriverArray1 rcm_init{};
  DriverArray1 wm_zt_init{};
  DriverArray1 wm_zm_init{};
  DriverArray1 em_init{};
  DriverArray1 exner_init{};
  DriverArray1 thvm_init{};
  DriverArray1 p_in_Pa_init{};
  DriverArray1 p_in_Pa_zm_init{};
  DriverArray1 rho_init{};
  DriverArray1 rho_zm_init{};
  DriverArray1 rho_ds_zm_init{};
  DriverArray1 rho_ds_zt_init{};
  DriverArray1 invrs_rho_ds_zm_init{};
  DriverArray1 invrs_rho_ds_zt_init{};
  DriverArray1 thv_ds_zm_init{};
  DriverArray1 thv_ds_zt_init{};
  DriverArray1 rtm_ref_init{};
  DriverArray1 thlm_ref_init{};
  DriverArray1 um_ref_init{};
  DriverArray1 vm_ref_init{};
  DriverArray1 Ncm_init{};
  DriverArray1 Nc_in_cloud_init{};
  DriverArray1 Nccnm_init{};
  DriverArray1 p_sfc{}; // surface pressure              [Pa]
  DriverArray1 T_sfc{};
  DriverArray1 fcor{};   // Traditional Coriolis parameter    [s^-1]
  DriverArray1 fcor_y{}; // Nontraditional Coriolis parameter [s^-1]
  DriverArray1 sfc_elevation{};
  DriverArray1 zm_init{};
  DriverArray1 zm_top{};
  DriverArray1 deltaz{};
  DriverArray2 sclrm_init{};
  DriverArray2 edsclrm_init{};
  DriverArray2 corr_array_n_cloud{};
  DriverArray2 corr_array_n_below{};
  int iinit{};
  DriverInt3 X_mixt_comp_all_levs{}; // Which mixture component a sample is in
  grid gr{};
  grid gr_desc{};
  grid gr_dycore{}; // only set and used, if l_add_dycore_grid is .true.
  grid gr_output{};
  grid gr_fixed_min{};              // the read in grid for the minimum grid density profile
  pdf_parameter pdf_params{};       // PDF parameters (thermodynamic levels)    [units vary]
  pdf_parameter pdf_params_zm{};    // PDF parameters on momentum levels        [units vary]
  DriverArray1 lmin{};              // Min. value for the length scale [m]
  DriverArray1 mixt_frac_max_mag{}; // Maximum allowable magnitude of mixture fraction [-]
  implicit_coefs_terms pdf_implicit_coefs_terms{}; // Implicit coefs / explicit terms [units vary]
  precipitation_fractions precip_fracs{};          // Precipitation fractions      [-]
  nu_vertical_res_dep nu_vert_res_dep{};           // Vertical resolution dependent nu values
  hydromet_pdf_parameter_collection
      hydromet_pdf_params{}; // Hydrometeor PDF parameters      [units vary]
  real_t sens_ht{};          // sensible heat flux                          [K m/s]
  real_t latent_ht{};        // latent heat flux                            [m/s]
  real_t vert_decorr_coef{}; // Empirically defined de-correlation constant [-]
  DriverArray1 deep_soil_T_in_K{};
  DriverArray1 sfc_soil_T_in_K{};
  DriverArray1 veg_T_in_K{};
  DriverArray1 wpthlp_sfc{};        // w' theta_l' at surface   [(m K)/s]
  DriverArray1 wprtp_sfc{};         // w' r_t' at surface       [(kg m)/( kg s)]
  DriverArray1 upwp_sfc{};          // u'w' at surface          [m^2/s^2]
  DriverArray1 vpwp_sfc{};          // v'w' at surface          [m^2/s^2]
  DriverArray1 upwp_sfc_pert{};     // pertubed u'w' at surface    [m^2/s^2]
  DriverArray1 vpwp_sfc_pert{};     // pertubed v'w' at surface    [m^2/s^2]
  DriverArray1 gr_dens_z{};         // levels at which the grid density is given      [m]
  DriverArray1 gr_dens{};           // grid density at the given levels in gr_dens_z  [1/m]
  DriverArray2 gr_dens_old_z{};     // altitudes [m]
  DriverArray2 gr_dens_old{};       // densities [1/m]
  int Lscale_counter{};             // the number of values cumulated in cumulative_Lscale
  int richardson_num_counter{};     // the number of values cumulated in cumulative_richardson_num
  int chi_counter{};                // the number of values cumulated in cumulative_chi
  DriverArray1 cumulative_Lscale{}; // Sum of mixing-length refinement terms [1/m].
  DriverArray1 cumulative_richardson_num{}; // Sum of Richardson-number refinement terms [1/m].
  DriverArray1 cumulative_chi{};            // Sum of chi refinement terms [1/m].
  DriverArray1 alt_term{};
  DriverArray1 lscale_term{};
  DriverArray1 lscale_term_time_avg{};
  DriverArray1 chi_term{};
  DriverArray1 chi_term_time_avg{};
  DriverArray1 richardson_num_term{};
  DriverArray1 richardson_num_term_time_avg{};
  DriverArray1 norm_min_grid_dens{};
  DriverArray1 norm_grid_dens{};
  DriverArray2 um{};      // eastward grid-mean wind component (thermo. levs.)  [m/s]
  DriverArray2 upwp{};    // u'w' (momentum levels)                             [m^2/s^2]
  DriverArray2 vm{};      // northward grid-mean wind component (thermo. levs.) [m/s]
  DriverArray2 vpwp{};    // v'w' (momentum levels)                             [m^2/s^2]
  DriverArray2 up2{};     // u'^2 (momentum levels)                             [m^2/s^2]
  DriverArray2 vp2{};     // v'^2 (momentum levels)                             [m^2/s^2]
  DriverArray2 up3{};     // u'^3 (thermodynamic levels)                        [m^3/s^3]
  DriverArray2 vp3{};     // v'^3 (thermodynamic levels)                        [m^3/s^3]
  DriverArray2 wprtp{};   // w' r_t' (momentum levels)                          [(kg/kg) m/s]
  DriverArray2 wpthlp{};  // w'th_l' (momentum levels)                          [(m/s) K]
  DriverArray2 rtp2{};    // r_t'^2 (momentum levels)                           [(kg/kg)^2]
  DriverArray2 rtp3{};    // r_t'^3 (thermodynamic levels)                      [(kg/kg)^3]
  DriverArray2 thlp2{};   // th_l'^2 (momentum levels)                          [K^2]
  DriverArray2 thlp3{};   // th_l'^3 (thermodynamic levels)                     [K^3]
  DriverArray2 rtpthlp{}; // r_t'th_l' (momentum levels)                        [(kg/kg) K]
  DriverArray2 wp2{};     // w'^2 (momentum levels)                             [m^2/s^2]
  DriverArray2 wp2_zt{};  // w'^2 on thermo. grid                                  [m^2/s^2]
  DriverArray2 thvm{};    // Virtual potential temperature                  [K]
  DriverArray2 exner{};   // Exner function (thermodynamic levels)          [-]
  DriverArray2 rtm{};     // total water mixing ratio, r_t (thermo. levels) [kg/kg]
  DriverArray2 thlm{};    // liq. water pot. temp., th_l (thermo. levels)   [K]
  DriverArray2 rcm{};     // cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
  DriverArray2 wp3{};     // w'^3 (thermodynamic levels)                    [m^3/s^3]
  DriverArray2 delta_zm{};
  DriverArray2 p_in_Pa{};             // Air pressure (thermodynamic levels)         [Pa]
  DriverArray2 cloud_frac{};          // cloud fraction (thermodynamic levels)       [-]
  DriverArray2 wpthvp{};              // < w' th_v' > (momentum levels)              [kg/kg K]
  DriverArray2 wp2thvp{};             // < w'^2 th_v' > (thermodynamic levels)       [m^2/s^2 K]
  DriverArray2 wp2up{};               // < w'^2 u' > (thermodynamic levels)          [m^3/s^3]
  DriverArray2 rtpthvp{};             // < r_t' th_v' > (momentum levels)            [kg/kg K]
  DriverArray2 thlpthvp{};            // < th_l' th_v' > (momentum levels)           [K^2]
  DriverArray2 wp2rtp{};              // w'^2 rt' (thermodynamic levels)             [m^2/s^2 kg/kg]
  DriverArray2 wp2thlp{};             // w'^2 thl' (thermodynamic levels)            [m^2/s^2 K]
  DriverArray2 uprcp{};               // < u' r_c' > (momentum levels)               [(m/s)(kg/kg)]
  DriverArray2 vprcp{};               // < v' r_c' > (momentum levels)               [(m/s)(kg/kg)]
  DriverArray2 rc_coef_zm{};          // Coefficient of X'r_c' in Eq. (34) (m-levs.) [K/(kg/kg)]
  DriverArray2 wp4{};                 // w'^4 (momentum levels)                      [m^4/s^4]
  DriverArray2 wpup2{};               // w'u'^2 (thermodynamic levels)               [m^3/s^3]
  DriverArray2 wpvp2{};               // w'v'^2 (thermodynamic levels)               [m^3/s^3]
  DriverArray2 wp2up2{};              // w'^2 u'^2 (momentum levels)                 [m^4/s^4]
  DriverArray2 wp2vp2{};              // w'^2 v'^2 (momentum levels)                 [m^4/s^4]
  DriverArray2 rho_ds_zt{};           // Dry, static density on thermo. levels      [kg/m^3]
  DriverArray2 thv_ds_zt{};           // Dry, base-state theta_v on thermo levs.    [K]
  DriverArray2 wm_zm{};               // vertical mean wind comp. on momentum levs  [m/s]
  DriverArray2 wm_zt{};               // vertical mean wind comp. on thermo. levs   [m/s]
  DriverArray2 rho{};                 // Air density on thermodynamic levels        [kg/m^3]
  DriverArray2 rho_zm{};              // Air density on momentum levels             [kg/m^3]
  DriverArray2 rho_ds_zm{};           // Dry, static density on momentum levels     [kg/m^3]
  DriverArray2 invrs_rho_ds_zm{};     // Inverse dry, static density on m-levs.     [m^3/kg]
  DriverArray2 invrs_rho_ds_zt{};     // Inverse dry, static density on thermo levs.[m^3/kg]
  DriverArray2 thv_ds_zm{};           // Dry, base-state theta_v on momentum levs.  [K]
  DriverArray2 thlm_forcing{};        // liq. wat. pot. temp. forcing (thermo. levs)[K/s]
  DriverArray2 rtm_forcing{};         // total water forcing (thermo. levels)       [(kg/kg)/s]
  DriverArray2 um_forcing{};          // eastward wind forcing (thermo. levels)     [m/s/s]
  DriverArray2 vm_forcing{};          // northward wind forcing (thermo. levels)    [m/s/s]
  DriverArray2 wprtp_forcing{};       // total water turbulent flux forcing (m-levs)[m*K/s^2]
  DriverArray2 wpthlp_forcing{};      // liq pot temp turb flux forcing (m-levs)    [m(kg/kg)/s^2]
  DriverArray2 rtp2_forcing{};        // total water variance forcing (m-levs)      [(kg/kg)^2/s]
  DriverArray2 thlp2_forcing{};       // liq pot temp variance forcing (m-levs)     [K^2/s]
  DriverArray2 rtpthlp_forcing{};     // <r_t'th_l'> covariance forcing (m-levs)    [K(kg/kg)/s]
  DriverArray2 um_pert{};             // perturbed <u>       [m/s]
  DriverArray2 vm_pert{};             // perturbed <v>       [m/s]
  DriverArray2 upwp_pert{};           // perturbed <u'w'>    [m^2/s^2]
  DriverArray2 vpwp_pert{};           // perturbed <v'w'>    [m^2/s^2]
  DriverArray2 Skw_zm_smooth{};       // Smoothed version of Skw_zm
  DriverArray2 Lscale{};              // Length scale                      [m]
  DriverArray2 thlprcp{};             // thl'rc'                                      [K kg/kg]
  DriverArray2 sigma_sqd_w{};         // PDF width parameter (momentum levels)        [-]
  DriverArray2 sigma_sqd_w_zt{};      // PDF width parameter interpolated to t-levs.  [-]
  DriverArray2 wprcp{};               // w'r_c' (momentum levels)              [(kg/kg) m/s]
  DriverArray2 w_up_in_cloud{};       // Average cloudy updraft velocity       [m/s]
  DriverArray2 w_down_in_cloud{};     // Average cloudy downdraft velocity     [m/s]
  DriverArray2 cloudy_updraft_frac{}; // cloudy updraft fraction               [-]
  DriverArray2 cloudy_downdraft_frac{}; // cloudy downdraft fraction             [-]
  DriverArray2 ice_supersat_frac{};     // ice cloud fraction (thermo. levels)   [-]
  DriverArray2 rcm_in_layer{};          // rcm within cloud layer                [kg/kg]
  DriverArray2 cloud_cover{};           // cloud cover                           [-]
  DriverArray2 invrs_tau_zm{};          // One divided by tau on zm levels       [1/s]
  DriverArray2 ug{};                    // u geostrophic wind                           [m/s]
  DriverArray2 vg{};                    // v geostrophic wind                           [m/s]
  DriverArray2 rtm_ref{};               // Initial total water mixing ratio             [kg/kg]
  DriverArray2 thlm_ref{};              // Initial liquid water potential temperature   [K]
  DriverArray2 um_ref{};                // Initial u wind                               [m/s]
  DriverArray2 vm_ref{};                // Initial v wind                               [m/s]
  DriverArray2 rcm_mc{};      // Tendency of liquid water due to microphysics      [kg/kg/s]
  DriverArray2 rvm_mc{};      // Tendency of vapor water due to microphysics       [kg/kg/s]
  DriverArray2 thlm_mc{};     // Tendency of liquid pot. temp. due to microphysics [K/s]
  DriverArray2 wprtp_mc{};    // Microphysics tendency for <w'rt'>                 [m*(kg/kg)/s^2]
  DriverArray2 wpthlp_mc{};   // Microphysics tendency for <w'thl'>                [m*K/s^2]
  DriverArray2 rtp2_mc{};     // Microphysics tendency for <rt'^2>                 [(kg/kg)^2/s]
  DriverArray2 thlp2_mc{};    // Microphysics tendency for <thl'^2>                [K^2/s]
  DriverArray2 rtpthlp_mc{};  // Microphysics tendency for <rt'thl'>               [K*(kg/kg)/s]
  DriverArray2 Ncm_mc{};      // Microphysics tendency for Ncm                     [num/kg/s]
  DriverArray2 Ncm{};         // Mean cloud droplet concentration, <N_c> (t-levs.)     [num/kg]
  DriverArray2 Nccnm{};       // Cloud condensation nuclei concentration (COAMPS/MG)   [num/kg]
  DriverArray2 em{};          // Turbulent Kinetic Energy (TKE)                        [m^2/s^2]
  DriverArray2 tau_zm{};      // Eddy dissipation time scale on momentum levels        [s]
  DriverArray2 tau_zt{};      // Eddy dissipation time scale on thermodynamic levels   [s]
  DriverArray2 Kh_zt{};       // Eddy diffusivity coefficient on thermodynamic levels  [m^2/s]
  DriverArray2 Kh_zm{};       // Eddy diffusivity coefficient on momentum levels       [m^2/s]
  DriverArray2 rfrzm{};       // Total ice-phase water mixing ratio                    [kg/kg]
  DriverArray2 rrm{};         // Overall mean rain water mixing ratio                  [kg/kg]
  DriverArray2 Nrm{};         // Overall mean rain drop concentration                  [num/kg]
  DriverArray2 Nc_in_cloud{}; // Mean (in-cloud) cloud droplet concentration           [num/kg]
  DriverArray2 wpNcp{}; // Covariance of w and N_c, <w'N_c'> (momentum levels)   [(m/s)(num/kg)]
  DriverArray2 radht{}; // SW + LW heating rate               [K/s]
  DriverArray2 wpsclrp_sfc{};   // Passive scalar flux at surface                [{units vary} m/s]
  DriverArray2 wpedsclrp_sfc{}; // Eddy-diffusion passive scalar flux at surface [{units vary} m/s]
  DriverArray3 sclrm{};         // Passive scalar mean (thermo. levels) [{units vary}]
  DriverArray3 wpsclrp{};       // w'sclr' (momentum levels)            [{units vary} m/s]
  DriverArray3 sclrp2{};        // sclr'^2 (momentum levels)            [{units vary}^2]
  DriverArray3 sclrp3{};        // sclr'^3 (thermodynamic levels)       [{units vary}^3]
  DriverArray3 sclrprtp{};      // sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
  DriverArray3 sclrpthlp{};     // sclr'thl' (momentum levels)          [{units vary} K]
  DriverArray3 sclrpthvp{};     // < sclr' th_v' > (momentum levels)   [units vary]
  DriverArray3 sclrm_forcing{}; // Passive scalar forcing          [{units vary}/s]
  DriverArray3 edsclrm{};       // Eddy passive scalar grid-mean (thermo. levels)   [units vary]
  DriverArray3 edsclrm_forcing{}; // Eddy-diffusion passive scalar forcing    [{units vary}/s]
  DriverArray3 lh_rt_clipped{};   // rt generated from silhs sample points
  DriverArray3 lh_thl_clipped{};  // thl generated from silhs sample points
  DriverArray3 lh_rc_clipped{};   // rc generated from silhs sample points
  DriverArray3 lh_rv_clipped{};   // rv generated from silhs sample points
  DriverArray3 lh_Nc_clipped{};   // Nc generated from silhs sample points
  DriverArray3 lh_sample_point_weights{}; // Weights for cloud weighted sampling
  DriverArray3 K_hm{};        // Eddy diffusivity coef. for hydrometeors on mom. levs. [m^2 s^-1]
  DriverArray3 wp2hmp{};      // Third moment:  <w'^2> * <hm'>       [(m/s)^2 <hm units>]
  DriverArray3 rtphmp_zt{};   // Covariance of rt and a hydrometeor  [(kg/kg) <hm units>]
  DriverArray3 thlphmp_zt{};  // Covariance of thl and a hydrometeor [K <hm units>]
  DriverArray3 hydromet{};    // Array of hydrometeors                [hm units]
  DriverArray3 hydrometp2{};  // Variance of a hydrometeor (m-levs.)  [<hm units>^2]
  DriverArray3 wphydrometp{}; // Covariance of w and a hydrometeor    [(m/s) <hm units>]
  DriverArray3 hydromet_mc{}; // Microphysics tendency for mean hydrometeors  [units/s]
  DriverArray3 hydromet_vel_zt{}; // Mean hydrometeor sed. velocity on thermo. levs. [m/s]
  DriverArray3 hydromet_vel_covar_zt_impc{}; // Imp. comp. <V_xx'x_x'> t-levs [m/s]
  DriverArray3 hydromet_vel_covar_zt_expc{}; // Exp. comp. <V_xx'x_x'> t-levs [units(m/s)]
  DriverArray3 mu_x_1_n{};            // Mean array (normal space): PDF vars. (comp. 1) [un. vary]
  DriverArray3 mu_x_2_n{};            // Mean array (normal space): PDF vars. (comp. 2) [un. vary]
  DriverArray3 sigma_x_1_n{};         // Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
  DriverArray3 sigma_x_2_n{};         // Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]
  DriverArray4 corr_array_1_n{};      // Corr. array (normal space) of PDF vars. (comp. 1)  [-]
  DriverArray4 corr_array_2_n{};      // Corr. array (normal space) of PDF vars. (comp. 2)  [-]
  DriverArray4 corr_cholesky_mtx_1{}; // Transposed corr. cholesky matrix, 1st comp. [-]
  DriverArray4 corr_cholesky_mtx_2{}; // Transposed corr. cholesky matrix, 2nd comp. [-]
  DriverArray4 X_nl_all_levs{};       // Lognormally distributed hydrometeors
};
} // namespace clubb_kokkos
