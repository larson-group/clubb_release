// Own the persistent Input_fields settings for one CLUBB case.
//
// Native module defaults are retained. inputfields_init performs only its source
// assignments; it does not infer a new reset of flags omitted by Native. The
// three simulation-date fields replace input_netcdf module globals. ClubbDriver
// owns this state; no C++ process-global cache exists.
#pragma once
#include <string>
#include <vector>
namespace clubb_kokkos {
struct inputfields_state {
  std::vector<std::string> stat_files; // Native fixed120 file names.
  std::string input_type;              // Native fixed10 input convention.
  bool l_input_um = false;
  bool l_input_vm = false;
  bool l_input_rtm = false;
  bool l_input_thlm = false;
  bool l_input_wp2 = false;
  bool l_input_wprtp = false;
  bool l_input_wpthlp = false;
  bool l_input_wp3 = false;
  bool l_input_rtp2 = false;
  bool l_input_rtp3 = false;
  bool l_input_thlp2 = false;
  bool l_input_thlp3 = false;
  bool l_input_rtpthlp = false;
  bool l_input_upwp = false;
  bool l_input_vpwp = false;
  bool l_input_ug = false;
  bool l_input_vg = false;
  bool l_input_rcm = false;
  bool l_input_wm_zt = false;
  bool l_input_exner = false;
  bool l_input_em = false;
  bool l_input_p = false;
  bool l_input_rho = false;
  bool l_input_rho_zm = false;
  bool l_input_rho_ds_zm = false;
  bool l_input_rho_ds_zt = false;
  bool l_input_thv_ds_zm = false;
  bool l_input_thv_ds_zt = false;
  bool l_input_Lscale = false;
  bool l_input_Lscale_up = false;
  bool l_input_Lscale_down = false;
  bool l_input_Kh_zt = false;
  bool l_input_Kh_zm = false;
  bool l_input_tau_zm = false;
  bool l_input_tau_zt = false;
  bool l_input_wpthvp = false;
  bool l_input_wp2thvp = false;
  bool l_input_wp2up = false;
  bool l_input_rtpthvp = false;
  bool l_input_thlpthvp = false;
  bool l_input_wp2rtp = false;
  bool l_input_wp2thlp = false;
  bool l_input_uprcp = false;
  bool l_input_vprcp = false;
  bool l_input_rc_coef_zm = false;
  bool l_input_wp4 = false;
  bool l_input_wpup2 = false;
  bool l_input_wpvp2 = false;
  bool l_input_wp2up2 = false;
  bool l_input_wp2vp2 = false;
  bool l_input_iss_frac = false;
  bool l_input_radht = false;
  bool l_input_w_1 = false;
  bool l_input_w_2 = false;
  bool l_input_varnce_w_1 = false;
  bool l_input_varnce_w_2 = false;
  bool l_input_rt_1 = false;
  bool l_input_rt_2 = false;
  bool l_input_varnce_rt_1 = false;
  bool l_input_varnce_rt_2 = false;
  bool l_input_thl_1 = false;
  bool l_input_thl_2 = false;
  bool l_input_varnce_thl_1 = false;
  bool l_input_varnce_thl_2 = false;
  bool l_input_mixt_frac = false;
  bool l_input_chi_1 = false;
  bool l_input_chi_2 = false;
  bool l_input_stdev_chi_1 = false;
  bool l_input_stdev_chi_2 = false;
  bool l_input_rc_1 = false;
  bool l_input_rc_2 = false;
  bool l_input_w_1_zm = false;
  bool l_input_w_2_zm = false;
  bool l_input_varnce_w_1_zm = false;
  bool l_input_varnce_w_2_zm = false;
  bool l_input_mixt_frac_zm = false;
  bool l_input_thvm = false;
  bool l_input_rrm = false;
  bool l_input_Nrm = false;
  bool l_input_Ncm = false;
  bool l_input_rsm = false;
  bool l_input_rim = false;
  bool l_input_Nsm = false;
  bool l_input_Ngm = false;
  bool l_input_rgm = false;
  bool l_input_Nccnm = false;
  bool l_input_Nim = false;
  bool l_input_rrp2 = false;
  bool l_input_Nrp2 = false;
  bool l_input_wprrp = false;
  bool l_input_wpNrp = false;
  bool l_input_thlm_forcing = false;
  bool l_input_rtm_forcing = false;
  bool l_input_up2 = false;
  bool l_input_vp2 = false;
  bool l_input_sigma_sqd_w = false;
  bool l_input_cloud_frac = false;
  bool l_input_sigma_sqd_w_zt = false;
  bool l_input_veg_T_in_K = false;
  bool l_input_deep_soil_T_in_K = false;
  bool l_input_sfc_soil_T_in_K = false;
  bool l_input_wprtp_forcing = false;
  bool l_input_wpthlp_forcing = false;
  bool l_input_rtp2_forcing = false;
  bool l_input_thlp2_forcing = false;
  bool l_input_rtpthlp_forcing = false;
  bool l_input_thlprcp = false;
  bool l_input_rcm_mc = false;
  bool l_input_rvm_mc = false;
  bool l_input_thlm_mc = false;
  bool l_input_wprtp_mc = false;
  bool l_input_wpthlp_mc = false;
  bool l_input_rtp2_mc = false;
  bool l_input_thlp2_mc = false;
  bool l_input_rtpthlp_mc = false;
  int stats_input_type = -999;            // Native private format selector, used by package bodies.
  int clubb_day, clubb_month, clubb_year; // Caller sets the simulation start date.
};
} // namespace clubb_kokkos
