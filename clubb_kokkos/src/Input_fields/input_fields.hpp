// Restore model fields from Native-compatible NetCDF profiles.
//
// inputfields_init reads /setfields/; set_filenames selects files and their
// convention. compute_timestep locates a Native record, stat_fields_reader
// restores one column, and get_clubb_variable_interpolated reads one optional
// profile. cleanup_input_fields releases file names at case shutdown.
// General inputfields execution remains driver-gated as in JAX; restart is
// supported. All Native source branches are retained in the implementation.
// Fields use the Native one-column(level[,species]) boundary; PDF fields retain
// Native row1 writes. The case context owns settings instead of C++ globals.
#pragma once
#include "scratch_workspace.hpp"
#include "inputfields_state.hpp"
#include "stat_file_utils.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "CLUBB_core/corr_varnce_module.hpp"
#include "CLUBB_core/pdf_parameter_module.hpp"
namespace clubb_kokkos {
void set_filenames(const std::string& file_prefix,  // In
                   inputfields_state& inputfields); // Native module InOut
void stat_fields_reader(ScratchWorkspace& workspace, const grid& gr, int timestep, int hydromet_dim,
                        const hm_metadata_type& hm_metadata,                         // In
                        const std::string& microphys_scheme, bool l_predict_Nc,      // In
                        Real1D um, Real1D upwp, Real1D vm, Real1D vpwp,              // InOut
                        Real1D up2, Real1D vp2, Real1D rtm,                          // InOut
                        Real1D wprtp, Real1D thlm, Real1D wpthlp,                    // InOut
                        Real1D rtp2, Real1D rtp3,                                    // InOut
                        Real1D thlp2, Real1D thlp3, Real1D rtpthlp,                  // InOut
                        Real1D wp2, Real1D wp3,                                      // InOut
                        Real1D p_in_Pa, Real1D exner, Real1D rcm, Real1D cloud_frac, // InOut
                        Real1D wpthvp, Real1D wp2thvp, Real1D wp2up, Real1D rtpthvp,
                        Real1D thlpthvp,                                               // InOut
                        Real1D wp2rtp, Real1D wp2thlp, Real1D uprcp, Real1D vprcp,     // InOut
                        Real1D rc_coef_zm, Real1D wp4, Real1D wpup2,                   // InOut
                        Real1D wpvp2, Real1D wp2up2,                                   // InOut
                        Real1D wp2vp2, Real1D ice_supersat_frac,                       // InOut
                        Real1D wm_zt, Real1D rho, Real1D rho_zm, Real1D rho_ds_zm,     // InOut
                        Real1D rho_ds_zt, Real1D thv_ds_zm, Real1D thv_ds_zt,          // InOut
                        Real1D thlm_forcing, Real1D rtm_forcing, Real1D wprtp_forcing, // InOut
                        Real1D wpthlp_forcing, Real1D rtp2_forcing,                    // InOut
                        Real1D thlp2_forcing, Real1D rtpthlp_forcing,                  // InOut
                        Real2D hydromet, Real2D hydrometp2, Real2D wphydrometp,        // InOut
                        Real1D Ncm, Real1D Nccnm, Real1D thvm, Real1D em,              // InOut
                        Real1D tau_zm, Real1D tau_zt,                                  // InOut
                        Real1D Kh_zt, Real1D Kh_zm, Real1D ug, Real1D vg,              // InOut
                        Real1D thlprcp,                                                // InOut
                        Real1D sigma_sqd_w, Real1D sigma_sqd_w_zt, Real1D radht,       // InOut
                        real_t& deep_soil_T_in_K, real_t& sfc_soil_T_in_K,
                        real_t& veg_T_in_K,                                      // InOut
                        pdf_parameter& pdf_params, pdf_parameter& pdf_params_zm, // InOut
                        const inputfields_state& inputfields);                   // Native module In
void compute_timestep(ScratchWorkspace& workspace, const std::string& filename, bool l_restart,
                      double time,                           // In
                      int& nearest_timestep,                 // Out
                      const inputfields_state& inputfields); // Native module In
void get_clubb_variable_interpolated(ScratchWorkspace& workspace, bool l_input_var,
                                     const std::string& filename,
                                     const std::string& varname, // In
                                     int vardim, int timestep,
                                     ConstInputProfile clubb_heights,       // In
                                     InputProfile variable_interpolated,    // InOut
                                     bool& l_read_error,                    // Out
                                     const inputfields_state& inputfields); // Native module In
void inputfields_init(int iunit, const std::string& namelist_filename, int day, int month,
                      int year,                            // In
                      inputfields_state& inputfields);     // Native module InOut
void cleanup_input_fields(inputfields_state& inputfields); // Native module InOut

template <class Profile, class Hydromet>
void stat_fields_reader(ScratchWorkspace& workspace, // InOut
                        const grid& gr, int timestep, int hydromet_dim,
                        const hm_metadata_type& hm_metadata,                             // In
                        const std::string& microphys_scheme, bool l_predict_Nc,          // In
                        Profile um, Profile upwp, Profile vm, Profile vpwp,              // InOut
                        Profile up2, Profile vp2, Profile rtm,                           // InOut
                        Profile wprtp, Profile thlm, Profile wpthlp,                     // InOut
                        Profile rtp2, Profile rtp3,                                      // InOut
                        Profile thlp2, Profile thlp3, Profile rtpthlp,                   // InOut
                        Profile wp2, Profile wp3,                                        // InOut
                        Profile p_in_Pa, Profile exner, Profile rcm, Profile cloud_frac, // InOut
                        Profile wpthvp, Profile wp2thvp, Profile wp2up, Profile rtpthvp,
                        Profile thlpthvp,                                                 // InOut
                        Profile wp2rtp, Profile wp2thlp, Profile uprcp, Profile vprcp,    // InOut
                        Profile rc_coef_zm, Profile wp4, Profile wpup2,                   // InOut
                        Profile wpvp2, Profile wp2up2,                                    // InOut
                        Profile wp2vp2, Profile ice_supersat_frac,                        // InOut
                        Profile wm_zt, Profile rho, Profile rho_zm, Profile rho_ds_zm,    // InOut
                        Profile rho_ds_zt, Profile thv_ds_zm, Profile thv_ds_zt,          // InOut
                        Profile thlm_forcing, Profile rtm_forcing, Profile wprtp_forcing, // InOut
                        Profile wpthlp_forcing, Profile rtp2_forcing,                     // InOut
                        Profile thlp2_forcing, Profile rtpthlp_forcing,                   // InOut
                        Hydromet hydromet, Hydromet hydrometp2, Hydromet wphydrometp,     // InOut
                        Profile Ncm, Profile Nccnm, Profile thvm, Profile em,             // InOut
                        Profile tau_zm, Profile tau_zt,                                   // InOut
                        Profile Kh_zt, Profile Kh_zm, Profile ug, Profile vg,             // InOut
                        Profile thlprcp,                                                  // InOut
                        Profile sigma_sqd_w, Profile sigma_sqd_w_zt, Profile radht,       // InOut
                        real_t& deep_soil_T_in_K, real_t& sfc_soil_T_in_K,
                        real_t& veg_T_in_K,                                      // InOut
                        pdf_parameter& pdf_params, pdf_parameter& pdf_params_zm, // InOut
                        const inputfields_state& inputfields)                    // Native module In
{
  auto workspace_scope        = workspace.scope();
  auto um_exec                = to_execution_space(workspace, um);
  auto upwp_exec              = to_execution_space(workspace, upwp);
  auto vm_exec                = to_execution_space(workspace, vm);
  auto vpwp_exec              = to_execution_space(workspace, vpwp);
  auto up2_exec               = to_execution_space(workspace, up2);
  auto vp2_exec               = to_execution_space(workspace, vp2);
  auto rtm_exec               = to_execution_space(workspace, rtm);
  auto wprtp_exec             = to_execution_space(workspace, wprtp);
  auto thlm_exec              = to_execution_space(workspace, thlm);
  auto wpthlp_exec            = to_execution_space(workspace, wpthlp);
  auto rtp2_exec              = to_execution_space(workspace, rtp2);
  auto rtp3_exec              = to_execution_space(workspace, rtp3);
  auto thlp2_exec             = to_execution_space(workspace, thlp2);
  auto thlp3_exec             = to_execution_space(workspace, thlp3);
  auto rtpthlp_exec           = to_execution_space(workspace, rtpthlp);
  auto wp2_exec               = to_execution_space(workspace, wp2);
  auto wp3_exec               = to_execution_space(workspace, wp3);
  auto p_in_Pa_exec           = to_execution_space(workspace, p_in_Pa);
  auto exner_exec             = to_execution_space(workspace, exner);
  auto rcm_exec               = to_execution_space(workspace, rcm);
  auto cloud_frac_exec        = to_execution_space(workspace, cloud_frac);
  auto wpthvp_exec            = to_execution_space(workspace, wpthvp);
  auto wp2thvp_exec           = to_execution_space(workspace, wp2thvp);
  auto wp2up_exec             = to_execution_space(workspace, wp2up);
  auto rtpthvp_exec           = to_execution_space(workspace, rtpthvp);
  auto thlpthvp_exec          = to_execution_space(workspace, thlpthvp);
  auto wp2rtp_exec            = to_execution_space(workspace, wp2rtp);
  auto wp2thlp_exec           = to_execution_space(workspace, wp2thlp);
  auto uprcp_exec             = to_execution_space(workspace, uprcp);
  auto vprcp_exec             = to_execution_space(workspace, vprcp);
  auto rc_coef_zm_exec        = to_execution_space(workspace, rc_coef_zm);
  auto wp4_exec               = to_execution_space(workspace, wp4);
  auto wpup2_exec             = to_execution_space(workspace, wpup2);
  auto wpvp2_exec             = to_execution_space(workspace, wpvp2);
  auto wp2up2_exec            = to_execution_space(workspace, wp2up2);
  auto wp2vp2_exec            = to_execution_space(workspace, wp2vp2);
  auto ice_supersat_frac_exec = to_execution_space(workspace, ice_supersat_frac);
  auto wm_zt_exec             = to_execution_space(workspace, wm_zt);
  auto rho_exec               = to_execution_space(workspace, rho);
  auto rho_zm_exec            = to_execution_space(workspace, rho_zm);
  auto rho_ds_zm_exec         = to_execution_space(workspace, rho_ds_zm);
  auto rho_ds_zt_exec         = to_execution_space(workspace, rho_ds_zt);
  auto thv_ds_zm_exec         = to_execution_space(workspace, thv_ds_zm);
  auto thv_ds_zt_exec         = to_execution_space(workspace, thv_ds_zt);
  auto thlm_forcing_exec      = to_execution_space(workspace, thlm_forcing);
  auto rtm_forcing_exec       = to_execution_space(workspace, rtm_forcing);
  auto wprtp_forcing_exec     = to_execution_space(workspace, wprtp_forcing);
  auto wpthlp_forcing_exec    = to_execution_space(workspace, wpthlp_forcing);
  auto rtp2_forcing_exec      = to_execution_space(workspace, rtp2_forcing);
  auto thlp2_forcing_exec     = to_execution_space(workspace, thlp2_forcing);
  auto rtpthlp_forcing_exec   = to_execution_space(workspace, rtpthlp_forcing);
  auto hydromet_exec          = to_execution_space(workspace, hydromet);
  auto hydrometp2_exec        = to_execution_space(workspace, hydrometp2);
  auto wphydrometp_exec       = to_execution_space(workspace, wphydrometp);
  auto Ncm_exec               = to_execution_space(workspace, Ncm);
  auto Nccnm_exec             = to_execution_space(workspace, Nccnm);
  auto thvm_exec              = to_execution_space(workspace, thvm);
  auto em_exec                = to_execution_space(workspace, em);
  auto tau_zm_exec            = to_execution_space(workspace, tau_zm);
  auto tau_zt_exec            = to_execution_space(workspace, tau_zt);
  auto Kh_zt_exec             = to_execution_space(workspace, Kh_zt);
  auto Kh_zm_exec             = to_execution_space(workspace, Kh_zm);
  auto ug_exec                = to_execution_space(workspace, ug);
  auto vg_exec                = to_execution_space(workspace, vg);
  auto thlprcp_exec           = to_execution_space(workspace, thlprcp);
  auto sigma_sqd_w_exec       = to_execution_space(workspace, sigma_sqd_w);
  auto sigma_sqd_w_zt_exec    = to_execution_space(workspace, sigma_sqd_w_zt);
  auto radht_exec             = to_execution_space(workspace, radht);
  stat_fields_reader(
      workspace, gr, timestep, hydromet_dim, hm_metadata, microphys_scheme, l_predict_Nc, um_exec,
      upwp_exec, vm_exec, vpwp_exec, up2_exec, vp2_exec, rtm_exec, wprtp_exec, thlm_exec,
      wpthlp_exec, rtp2_exec, rtp3_exec, thlp2_exec, thlp3_exec, rtpthlp_exec, wp2_exec, wp3_exec,
      p_in_Pa_exec, exner_exec, rcm_exec, cloud_frac_exec, wpthvp_exec, wp2thvp_exec, wp2up_exec,
      rtpthvp_exec, thlpthvp_exec, wp2rtp_exec, wp2thlp_exec, uprcp_exec, vprcp_exec,
      rc_coef_zm_exec, wp4_exec, wpup2_exec, wpvp2_exec, wp2up2_exec, wp2vp2_exec,
      ice_supersat_frac_exec, wm_zt_exec, rho_exec, rho_zm_exec, rho_ds_zm_exec, rho_ds_zt_exec,
      thv_ds_zm_exec, thv_ds_zt_exec, thlm_forcing_exec, rtm_forcing_exec, wprtp_forcing_exec,
      wpthlp_forcing_exec, rtp2_forcing_exec, thlp2_forcing_exec, rtpthlp_forcing_exec,
      hydromet_exec, hydrometp2_exec, wphydrometp_exec, Ncm_exec, Nccnm_exec, thvm_exec, em_exec,
      tau_zm_exec, tau_zt_exec, Kh_zt_exec, Kh_zm_exec, ug_exec, vg_exec, thlprcp_exec,
      sigma_sqd_w_exec, sigma_sqd_w_zt_exec, radht_exec, deep_soil_T_in_K, sfc_soil_T_in_K,
      veg_T_in_K, pdf_params, pdf_params_zm, inputfields);
  from_execution_space(workspace, um, um_exec);
  from_execution_space(workspace, upwp, upwp_exec);
  from_execution_space(workspace, vm, vm_exec);
  from_execution_space(workspace, vpwp, vpwp_exec);
  from_execution_space(workspace, up2, up2_exec);
  from_execution_space(workspace, vp2, vp2_exec);
  from_execution_space(workspace, rtm, rtm_exec);
  from_execution_space(workspace, wprtp, wprtp_exec);
  from_execution_space(workspace, thlm, thlm_exec);
  from_execution_space(workspace, wpthlp, wpthlp_exec);
  from_execution_space(workspace, rtp2, rtp2_exec);
  from_execution_space(workspace, rtp3, rtp3_exec);
  from_execution_space(workspace, thlp2, thlp2_exec);
  from_execution_space(workspace, thlp3, thlp3_exec);
  from_execution_space(workspace, rtpthlp, rtpthlp_exec);
  from_execution_space(workspace, wp2, wp2_exec);
  from_execution_space(workspace, wp3, wp3_exec);
  from_execution_space(workspace, p_in_Pa, p_in_Pa_exec);
  from_execution_space(workspace, exner, exner_exec);
  from_execution_space(workspace, rcm, rcm_exec);
  from_execution_space(workspace, cloud_frac, cloud_frac_exec);
  from_execution_space(workspace, wpthvp, wpthvp_exec);
  from_execution_space(workspace, wp2thvp, wp2thvp_exec);
  from_execution_space(workspace, wp2up, wp2up_exec);
  from_execution_space(workspace, rtpthvp, rtpthvp_exec);
  from_execution_space(workspace, thlpthvp, thlpthvp_exec);
  from_execution_space(workspace, wp2rtp, wp2rtp_exec);
  from_execution_space(workspace, wp2thlp, wp2thlp_exec);
  from_execution_space(workspace, uprcp, uprcp_exec);
  from_execution_space(workspace, vprcp, vprcp_exec);
  from_execution_space(workspace, rc_coef_zm, rc_coef_zm_exec);
  from_execution_space(workspace, wp4, wp4_exec);
  from_execution_space(workspace, wpup2, wpup2_exec);
  from_execution_space(workspace, wpvp2, wpvp2_exec);
  from_execution_space(workspace, wp2up2, wp2up2_exec);
  from_execution_space(workspace, wp2vp2, wp2vp2_exec);
  from_execution_space(workspace, ice_supersat_frac, ice_supersat_frac_exec);
  from_execution_space(workspace, wm_zt, wm_zt_exec);
  from_execution_space(workspace, rho, rho_exec);
  from_execution_space(workspace, rho_zm, rho_zm_exec);
  from_execution_space(workspace, rho_ds_zm, rho_ds_zm_exec);
  from_execution_space(workspace, rho_ds_zt, rho_ds_zt_exec);
  from_execution_space(workspace, thv_ds_zm, thv_ds_zm_exec);
  from_execution_space(workspace, thv_ds_zt, thv_ds_zt_exec);
  from_execution_space(workspace, thlm_forcing, thlm_forcing_exec);
  from_execution_space(workspace, rtm_forcing, rtm_forcing_exec);
  from_execution_space(workspace, wprtp_forcing, wprtp_forcing_exec);
  from_execution_space(workspace, wpthlp_forcing, wpthlp_forcing_exec);
  from_execution_space(workspace, rtp2_forcing, rtp2_forcing_exec);
  from_execution_space(workspace, thlp2_forcing, thlp2_forcing_exec);
  from_execution_space(workspace, rtpthlp_forcing, rtpthlp_forcing_exec);
  from_execution_space(workspace, hydromet, hydromet_exec);
  from_execution_space(workspace, hydrometp2, hydrometp2_exec);
  from_execution_space(workspace, wphydrometp, wphydrometp_exec);
  from_execution_space(workspace, Ncm, Ncm_exec);
  from_execution_space(workspace, Nccnm, Nccnm_exec);
  from_execution_space(workspace, thvm, thvm_exec);
  from_execution_space(workspace, em, em_exec);
  from_execution_space(workspace, tau_zm, tau_zm_exec);
  from_execution_space(workspace, tau_zt, tau_zt_exec);
  from_execution_space(workspace, Kh_zt, Kh_zt_exec);
  from_execution_space(workspace, Kh_zm, Kh_zm_exec);
  from_execution_space(workspace, ug, ug_exec);
  from_execution_space(workspace, vg, vg_exec);
  from_execution_space(workspace, thlprcp, thlprcp_exec);
  from_execution_space(workspace, sigma_sqd_w, sigma_sqd_w_exec);
  from_execution_space(workspace, sigma_sqd_w_zt, sigma_sqd_w_zt_exec);
  from_execution_space(workspace, radht, radht_exec);
}

template <class Heights, class Variable>
void get_clubb_variable_interpolated(ScratchWorkspace& workspace, bool l_input_var,
                                     const std::string& filename,
                                     const std::string& varname,                      // In
                                     int vardim, int timestep, Heights clubb_heights, // In
                                     Variable variable_interpolated,                  // InOut
                                     bool& l_read_error,                              // Out
                                     const inputfields_state& inputfields) // Native module In
{
  // Preserve InOut fields when disabled; otherwise stage only inaccessible memory.
  if (!l_input_var) {
    l_read_error = false;
    return;
  }
  auto scope         = workspace.scope();
  const auto heights = [&]() {
    if constexpr (Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                             typename Heights::memory_space>::accessible)
      return ConstInputProfile(clubb_heights);
    else
      return ConstInputProfile(to_execution_space(workspace, clubb_heights));
  }();
  const auto output = [&]() {
    if constexpr (Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                             typename Variable::memory_space>::accessible)
      return InputProfile(variable_interpolated);
    else
      return InputProfile(to_execution_space(workspace, variable_interpolated));
  }();
  get_clubb_variable_interpolated(workspace, l_input_var, filename, varname, vardim, timestep,
                                  heights, output, l_read_error, inputfields);
  from_execution_space(workspace, variable_interpolated, output);
}
} // namespace clubb_kokkos
