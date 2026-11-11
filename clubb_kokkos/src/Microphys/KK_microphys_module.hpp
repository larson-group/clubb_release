// Coordinate local and analytically upscaled Khairoutdinov-Kogan microphysics.
//
// KK_local_microphys_driver calculates local tendencies and weighted statistics;
// KK_upscaled_microphys integrates means/covariances over the CLUBB PDF;
// KK_microphys_adjust limits source and evaporation tendencies to available water.
// KK_microphys_adj_terms_type records those adjustments [(hydrometeor units)/s].
// Column-first LayoutLeft arrays use (column,level[,species]); metadata indices
// retain Native one-based values. Appended parameter records replace Native
// module settings. Call-local host Views isolate the initial C++ engine boundary;
// the public adapters restore InOut/Out fields to their existing allocations.
// References: Khairoutdinov and Kogan(2000), MWR128,229–243; Griffin(2016), UWM1144.
#pragma once

#include "scratch_workspace.hpp"
#include "Microphys/microphys_host_types.hpp"
#include "Microphys/KK_microphys/KK_microphys.hpp"
#include "Microphys/advance_microphys_module.hpp"
#include "Microphys/Morrison_microphys/module_mp_graupel.hpp"
namespace clubb_kokkos::microphys_host {
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24>
inline void KK_local_microphys_driver(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int nzt, int hydromet_dim,
    const hm_metadata_type& hm_metadata, bool l_latin_hypercube, Array0 thlm, Array1 wm_zt,
    Array2 p_in_Pa, Array3 exner, Array4 rho, Array5 cloud_frac, Array6 w_std_dev, Array7 dzq,
    Array8 rcm, Array9 Ncm, Array10 chi, Array11 rvm, Array12 hydromet, int saturation_formula,
    Array13 sample_weight, stats_type& stats, Array14 hydromet_mc, Array15 hydromet_vel,
    Array16 Ncm_mc, Array17 rcm_mc, Array18 rvm_mc, Array19 thlm_mc, Array20 rrm_auto_diag,
    Array21 rrm_accr_diag, Array22 rrm_evap_diag, Array23 Nrm_auto_diag, Array24 Nrm_evap_diag,
    const parameters_microphys& microphys, const parameters_KK& kk);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24>
inline void KK_local_microphys_core(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int nzt, int hydromet_dim,
    const hm_metadata_type& hm_metadata, bool l_latin_hypercube, Array0 thlm, Array1 p_in_Pa,
    Array2 exner, Array3 rho, Array4 rcm, Array5 Ncm, Array6 chi, Array7 hydromet,
    int saturation_formula, Array8 hydromet_mc, Array9 hydromet_vel, Array10 Ncm_mc, Array11 rcm_mc,
    Array12 rvm_mc, Array13 thlm_mc, Array14 rrm_auto_diag, Array15 rrm_accr_diag,
    Array16 rrm_evap_diag, Array17 Nrm_auto_diag, Array18 Nrm_evap_diag, Array19 mvrr_diag,
    Array20 rrm_src_adj_diag, Array21 Nrm_src_adj_diag, Array22 rrm_evap_adj_diag,
    Array23 Nrm_evap_adj_diag, Array24 rrm_mc_nonadj_diag, const parameters_microphys& microphys,
    const parameters_KK& kk);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25>
inline void KK_upscaled_microphys(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int nzt, int nzm,
    int pdf_dim, int hydromet_dim, const hm_metadata_type& hm_metadata, Array0 wm_zt, Array1 rtm,
    Array2 thlm, Array3 p_in_Pa, Array4 exner, Array5 rho, Array6 rcm,
    const pdf_parameter& pdf_params, const hydromet_pdf_parameter_collection& hydromet_pdf_params,
    const precipitation_fractions& precip_fracs, Array7 hydromet, Array8 mu_x_1_n, Array9 mu_x_2_n,
    Array10 sigma_x_1_n, Array11 sigma_x_2_n, Array12 corr_array_1_n, Array13 corr_array_2_n,
    int saturation_formula, stats_type& stats, Array14 hydromet_mc, Array15 hydromet_vel,
    Array16 rcm_mc, Array17 rvm_mc, Array18 thlm_mc, Array19 hydromet_vel_covar_zt_impc,
    Array20 hydromet_vel_covar_zt_expc, Array21 wprtp_mc, Array22 wpthlp_mc, Array23 rtp2_mc,
    Array24 thlp2_mc, Array25 rtpthlp_mc, const parameters_microphys& microphys,
    const parameters_KK& kk);
inline void KK_tendency_coefs(real_t thlm, real_t exner, real_t p_in_Pa, real_t rho,
                              real_t& KK_evap_coef, real_t& KK_auto_coef, real_t& KK_accr_coef,
                              real_t& KK_mvr_coef, int saturation_formula, const parameters_KK& kk);
inline void KK_microphys_adjust(real_t dt, real_t exner, real_t rcm, real_t rrm, real_t Nrm,
                                real_t KK_evap_tndcy, real_t KK_auto_tndcy, real_t KK_accr_tndcy,
                                real_t KK_Nrm_evap_tndcy, real_t KK_Nrm_auto_tndcy,
                                bool l_src_adj_enabled, bool l_evap_adj_enabled, real_t& rrm_mc,
                                real_t& Nrm_mc, real_t& rvm_mc, real_t& rcm_mc, real_t& thlm_mc,
                                KK_microphys_adj_terms_type& adj_terms, const parameters_KK& kk);
template <class Array0, class Array1, class Array2, class Array3>
inline void KK_upscaled_sedimentation(int ngrdcol, int nzt, Array0 cloud_top_level,
                                      Array1 KK_mean_vol_rad, bool l_clip_positive_sed, Array2 Vrr,
                                      Array3 VNr);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5>
inline void KK_microphys_output(int ngrdcol, int nzt, int hydromet_dim,
                                const hm_metadata_type& hm_metadata, Array0 Vrr, Array1 VNr,
                                Array2 rrm_mc, Array3 Nrm_mc, Array4 hydromet_mc,
                                Array5 hydromet_vel);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24>
inline void KK_local_microphys_driver(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int nzt, int hydromet_dim,
    const hm_metadata_type& hm_metadata, bool l_latin_hypercube, Array0 thlm, Array1 wm_zt,
    Array2 p_in_Pa, Array3 exner, Array4 rho, Array5 cloud_frac, Array6 w_std_dev, Array7 dzq,
    Array8 rcm, Array9 Ncm, Array10 chi, Array11 rvm, Array12 hydromet, int saturation_formula,
    Array13 sample_weight, stats_type& stats, Array14 hydromet_mc, Array15 hydromet_vel,
    Array16 Ncm_mc, Array17 rcm_mc, Array18 rvm_mc, Array19 thlm_mc, Array20 rrm_auto_diag,
    Array21 rrm_accr_diag, Array22 rrm_evap_diag, Array23 Nrm_auto_diag, Array24 Nrm_evap_diag,
    const parameters_microphys& microphys, const parameters_KK& kk) {
  // Description:
  //   Run local K-K microphysics for all model columns and update statistics
  //   with full-column fields after the microphysics calculations finish.
  // In: gr
  // In: ngrdcol — Number of model columns
  // In: dt — Model time step duration                 [s]
  // In: nzt — Number of model thermodynamic vertical grid levels
  // In: hydromet_dim
  // In: hm_metadata
  // In: l_latin_hypercube — Flag to use Latin Hypercube interface
  // In: thlm (ngrdcol,nzt) — Mean liquid water potential temperature         [K]
  // In: wm_zt (ngrdcol,nzt) — Mean vertical velocity on thermodynamic levels  [m/s]
  // In: p_in_Pa (ngrdcol,nzt) — Pressure                                        [Pa]
  // In: exner (ngrdcol,nzt) — Exner function                           [-]
  // In: rho (ngrdcol,nzt) — Density                                         [kg/m^3]
  // In: cloud_frac (ngrdcol,nzt) — Cloud fraction                                     [-]
  // In: w_std_dev (ngrdcol,nzt) — Standard deviation of w (for LH interface)          [m/s]
  // In: dzq (ngrdcol,nzt) — Thickness between thermo. levels (for LH interface) [m]
  // In: rcm (ngrdcol,nzt) — Mean cloud water mixing ratio            [kg/kg]
  // In: Ncm (ngrdcol,nzt) — Mean cloud droplet conc., < N_c >               [num/kg]
  // In: chi (ngrdcol,nzt) — Mean extended liquid water mixing ratio         [kg/kg]
  // In: rvm (ngrdcol,nzt) — Mean water vapor mixing ratio (for LH interface)    [kg/kg]
  // In: hydromet (ngrdcol,nzt,hydromet_dim) — Hydrometeor mean, < h_m > (thermodynamic levels)  [units]
  // In: saturation_formula — Integer that stores the saturation formula to be used
  // In: sample_weight (ngrdcol,nzt) — SILHS sample weight; unity for an ordinary call     [-]
  // Inout: stats
  // Out: hydromet_mc (ngrdcol,nzt,hydromet_dim) — Hydrometeor time tendency          [(units vary)/s]
  // Out: hydromet_vel (ngrdcol,nzt,hydromet_dim) — Hydrometeor sedimentation velocity [m/s]
  // Out: Ncm_mc (ngrdcol,nzt) — Time tendency of cloud droplet concentration  [num/kg/s]
  // Out: rcm_mc (ngrdcol,nzt) — Time tendency of liquid water mixing ratio       [kg/kg/s]
  // Out: rvm_mc (ngrdcol,nzt) — Time tendency of vapor water mixing ratio        [kg/kg/s]
  // Out: thlm_mc (ngrdcol,nzt) — Time tendency of liquid potential temperature    [K/s]
  // Out: rrm_auto_diag (ngrdcol,nzt)
  // Out: rrm_accr_diag (ngrdcol,nzt)
  // Out: rrm_evap_diag (ngrdcol,nzt)
  // Out: Nrm_auto_diag (ngrdcol,nzt)
  // Out: Nrm_evap_diag (ngrdcol,nzt)
  // Case state: microphys replaces the Native parameters_microphys module.
  // Case state: kk replaces the Native parameters_KK module.

  auto workspace_scope = workspace.scope();
  HostReal2D mvrr_diag =
      workspace_scope.view<HostReal2D>("KK_local_microphys_driver::mvrr_diag", ngrdcol, nzt);
  HostReal2D rrm_src_adj_diag =
      workspace_scope.view<HostReal2D>("KK_local_microphys_driver::rrm_src_adj_diag", ngrdcol, nzt);
  HostReal2D Nrm_src_adj_diag =
      workspace_scope.view<HostReal2D>("KK_local_microphys_driver::Nrm_src_adj_diag", ngrdcol, nzt);
  HostReal2D rrm_evap_adj_diag = workspace_scope.view<HostReal2D>(
      "KK_local_microphys_driver::rrm_evap_adj_diag", ngrdcol, nzt);
  HostReal2D Nrm_evap_adj_diag = workspace_scope.view<HostReal2D>(
      "KK_local_microphys_driver::Nrm_evap_adj_diag", ngrdcol, nzt);
  HostReal2D rrm_mc_nonadj_diag = workspace_scope.view<HostReal2D>(
      "KK_local_microphys_driver::rrm_mc_nonadj_diag", ngrdcol, nzt);
  // Update weighted SILHS diagnostics under their source names.
  auto update_microphys_stat = [&](std::string name, auto value) {
    auto stats_scope = workspace.scope();
    std::string output_name;

    output_name = name;
    // When SILHS is non-interactive, output sample-mean diagnostics under
    // their SILHS names (e.g., lh_rrm_auto), leaving ordinary stats unchanged.
    if ((l_latin_hypercube && (microphys.lh_microphys_type == lh_microphys_non_interactive))) {
      if (microphys_detail::trim(name) == std::string("rrm_auto") ||
          microphys_detail::trim(name) == std::string("rrm_accr") ||
          microphys_detail::trim(name) == std::string("rrm_evap") ||
          microphys_detail::trim(name) == std::string("Nrm_auto") ||
          microphys_detail::trim(name) == std::string("Nrm_evap")) {
        output_name = (std::string("lh_") + microphys_detail::trim(name));
      } else if (microphys_detail::trim(name) == std::string("mvrr")) {
        output_name = std::string("lh_m_vol_rad_rain");
      } else if (microphys_detail::trim(name) == std::string("rrm_mc_nonadj")) {
        output_name = std::string("lh_rrm_mc_nonadj");
      } else {
        return;
      }
    }

    if (l_latin_hypercube) {
      HostReal2D array_actual_1 = stats_scope.view<HostReal2D>(
          Kokkos::WithoutInitializing,
          "kk_local_microphys_driver/update_microphys_stat::array_actual_1", value.extent_int(0),
          value.extent_int(1));
      for (int index_1 = 0; index_1 < value.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < value.extent_int(0); ++index_0) {
          array_actual_1(index_0, index_1) =
              (value(index_0, index_1) * sample_weight(index_0, index_1));
        }
      }
      stats_update(microphys_detail::trim(output_name), array_actual_1, stats, true);
    } else {
      stats_update(microphys_detail::trim(output_name), value, stats);
    }
  };

  microphys_host::KK_local_microphys_core(
      workspace, gr, ngrdcol, dt, nzt, hydromet_dim, hm_metadata, l_latin_hypercube, thlm, p_in_Pa,
      exner, rho, rcm, Ncm, chi, hydromet, saturation_formula, hydromet_mc, hydromet_vel, Ncm_mc,
      rcm_mc, rvm_mc, thlm_mc, rrm_auto_diag, rrm_accr_diag, rrm_evap_diag, Nrm_auto_diag,
      Nrm_evap_diag, mvrr_diag, rrm_src_adj_diag, Nrm_src_adj_diag, rrm_evap_adj_diag,
      Nrm_evap_adj_diag, rrm_mc_nonadj_diag, microphys, kk);

  //!! Output values for statistics
  update_microphys_stat(std::string("rrm_evap"), rrm_evap_diag);
  update_microphys_stat(std::string("rrm_auto"), rrm_auto_diag);
  update_microphys_stat(std::string("rrm_accr"), rrm_accr_diag);
  update_microphys_stat(std::string("mvrr"), mvrr_diag);
  update_microphys_stat(std::string("Nrm_evap"), Nrm_evap_diag);
  update_microphys_stat(std::string("Nrm_auto"), Nrm_auto_diag);
  update_microphys_stat(std::string("rrm_src_adj"), rrm_src_adj_diag);
  update_microphys_stat(std::string("Nrm_src_adj"), Nrm_src_adj_diag);
  update_microphys_stat(std::string("rrm_evap_adj"), rrm_evap_adj_diag);
  update_microphys_stat(std::string("Nrm_evap_adj"), Nrm_evap_adj_diag);
  update_microphys_stat(std::string("rrm_mc_nonadj"), rrm_mc_nonadj_diag);
}
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24>
inline void KK_local_microphys_core(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int nzt, int hydromet_dim,
    const hm_metadata_type& hm_metadata, bool l_latin_hypercube, Array0 thlm, Array1 p_in_Pa,
    Array2 exner, Array3 rho, Array4 rcm, Array5 Ncm, Array6 chi, Array7 hydromet,
    int saturation_formula, Array8 hydromet_mc, Array9 hydromet_vel, Array10 Ncm_mc, Array11 rcm_mc,
    Array12 rvm_mc, Array13 thlm_mc, Array14 rrm_auto_diag, Array15 rrm_accr_diag,
    Array16 rrm_evap_diag, Array17 Nrm_auto_diag, Array18 Nrm_evap_diag, Array19 mvrr_diag,
    Array20 rrm_src_adj_diag, Array21 Nrm_src_adj_diag, Array22 rrm_evap_adj_diag,
    Array23 Nrm_evap_adj_diag, Array24 rrm_mc_nonadj_diag, const parameters_microphys& microphys,
    const parameters_KK& kk) {
  // Description:
  // References:
  // Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
  //    Parameterization in a Large-Eddy Simulation Model of Marine
  //    Stratocumulus.  Mon. Wea. Rev., 128, 229--243.
  // Input Variables
  // Output Variables
  // Output Variables
  // Local Variables
  // In: gr
  // In: ngrdcol — Number of model columns
  // In: dt — Model time step duration                 [s]
  // In: nzt — Number of model vertical grid levels
  // In: hydromet_dim
  // In: hm_metadata
  // In: l_latin_hypercube — Flag to use Latin Hypercube interface
  // In: thlm (ngrdcol,nzt) — Mean liquid water potential temperature         [K]
  // In: p_in_Pa (ngrdcol,nzt) — Pressure                                        [Pa]
  // In: exner (ngrdcol,nzt) — Exner function                                  [-]
  // In: rho (ngrdcol,nzt) — Density                                         [kg/m^3]
  // In: rcm (ngrdcol,nzt) — Mean cloud water mixing ratio                   [kg/kg]
  // In: Ncm (ngrdcol,nzt) — Mean cloud droplet conc., < N_c >               [num/kg]
  // In: chi (ngrdcol,nzt) — Mean extended liquid water mixing ratio         [kg/kg]
  // In: hydromet (ngrdcol,nzt,hydromet_dim) — Hydrometeor species                      [units vary]
  // In: saturation_formula — Integer that stores the saturation formula to be used
  // Out: hydromet_mc (ngrdcol,nzt,hydromet_dim) — Hydrometeor time tendency          [(units vary)/s]
  // Out: hydromet_vel (ngrdcol,nzt,hydromet_dim) — Hydrometeor sedimentation velocity [m/s]
  // Out: Ncm_mc (ngrdcol,nzt) — Time tendency of cloud droplet concentration  [num/kg/s]
  // Out: rcm_mc (ngrdcol,nzt) — Time tendency of liquid water mixing ratio    [kg/kg/s]
  // Out: rvm_mc (ngrdcol,nzt) — Time tendency of vapor water mixing ratio     [kg/kg/s]
  // Out: thlm_mc (ngrdcol,nzt) — Time tendency of liquid potential temperature [K/s]
  // Out: rrm_auto_diag (ngrdcol,nzt)
  // Out: rrm_accr_diag (ngrdcol,nzt)
  // Out: rrm_evap_diag (ngrdcol,nzt)
  // Out: Nrm_auto_diag (ngrdcol,nzt)
  // Out: Nrm_evap_diag (ngrdcol,nzt)
  // Out: mvrr_diag (ngrdcol,nzt)
  // Out: rrm_src_adj_diag (ngrdcol,nzt)
  // Out: Nrm_src_adj_diag (ngrdcol,nzt)
  // Out: rrm_evap_adj_diag (ngrdcol,nzt)
  // Out: Nrm_evap_adj_diag (ngrdcol,nzt)
  // Out: rrm_mc_nonadj_diag (ngrdcol,nzt)
  // Case state: microphys replaces the Native parameters_microphys module.
  // Case state: kk replaces the Native parameters_KK module.

  auto workspace_scope        = workspace.scope();
  const auto hm_metadata_iirr = hm_metadata.iirr;
  const auto hm_metadata_iiri = hm_metadata.iiri;
  const auto hm_metadata_iiNr = hm_metadata.iiNr;
  HostReal2D KK_auto_tndcy =
      workspace_scope.view<HostReal2D>("KK_local_microphys_core::KK_auto_tndcy", ngrdcol,
                                       nzt); // Mean KK (dr_r/dt) due to autoconversion  [(kg/kg)/s]
  HostReal2D KK_accr_tndcy =
      workspace_scope.view<HostReal2D>("KK_local_microphys_core::KK_accr_tndcy", ngrdcol,
                                       nzt); // Mean KK (dr_r/dt) due to accretion       [(kg/kg)/s]
  HostReal2D KK_evap_tndcy =
      workspace_scope.view<HostReal2D>("KK_local_microphys_core::KK_evap_tndcy", ngrdcol,
                                       nzt); // Mean KK (dr_r/dt) due to evaporation     [(kg/kg)/s]
  HostReal2D KK_Nrm_evap_tndcy =
      workspace_scope.view<HostReal2D>("KK_local_microphys_core::KK_Nrm_evap_tndcy", ngrdcol,
                                       nzt); // Mean KK (dN_r/dt) due to evaporation  [(num/kg)/s]
  HostReal2D KK_Nrm_auto_tndcy =
      workspace_scope.view<HostReal2D>("KK_local_microphys_core::KK_Nrm_auto_tndcy", ngrdcol,
                                       nzt); // Mean KK (dN_r/dt) due to autoconv.    [(num/kg)/s]
  HostReal2D rrm =
      workspace_scope.view<HostReal2D>("KK_local_microphys_core::rrm", ngrdcol,
                                       nzt); // Mean rain water mixing ratio, < r_r > [kg/kg]
  HostReal2D Nrm =
      workspace_scope.view<HostReal2D>("KK_local_microphys_core::Nrm", ngrdcol,
                                       nzt); // Mean rain drop concentration, < N_r > [num/kg]
  HostReal2D Vrr = workspace_scope.view<HostReal2D>(
      "KK_local_microphys_core::Vrr", ngrdcol, nzt); // Mean sedimentation velocity of r_r    [m/s]
  HostReal2D VNr = workspace_scope.view<HostReal2D>(
      "KK_local_microphys_core::VNr", ngrdcol, nzt); // Mean sedimentation velocity of N_r    [m/s]
  HostReal2D rrm_mc =
      workspace_scope.view<HostReal2D>("KK_local_microphys_core::rrm_mc", ngrdcol,
                                       nzt); // Mean (dr_r/dt) due to microphysics    [(kg/kg)/s]
  HostReal2D Nrm_mc =
      workspace_scope.view<HostReal2D>("KK_local_microphys_core::Nrm_mc", ngrdcol,
                                       nzt); // Mean (dN_r/dt) due to microphysics    [(num/kg)/s]
  HostReal2D KK_mean_vol_rad =
      workspace_scope.view<HostReal2D>("KK_local_microphys_core::KK_mean_vol_rad", ngrdcol,
                                       nzt); // Mean KK rain drop mean volume radius     [m]
  real_t KK_evap_coef; // KK evaporation coefficient                  [(kg/kg)/s]
  real_t KK_auto_coef; // KK autoconversion coefficient               [(kg/kg)/s]
  real_t KK_accr_coef; // KK accretion coefficient                    [(kg/kg)/s]
  real_t KK_mvr_coef;  // KK mean volume radius coefficient           [m]
  HostAdj2D adj_terms("KK_local_microphys_core::adj_terms", ngrdcol,
                      nzt); // Adjustments to microphysics terms
  bool l_src_adj_enabled;   // Flag to enable rrm/Nrm source adjustment
  bool l_evap_adj_enabled;  // Flag to enable rrm/Nrm evaporation adjustment
  bool l_clip_positive_sed; // Clip positive Vrr and VNr terms to zero
  int i;                    // Column loop iterator
  int k;                    // Vertical loop iterator
  HostInt1D cloud_top_level = workspace_scope.view<HostInt1D>(
      "KK_local_microphys_core::cloud_top_level", ngrdcol); // Vertical level index of cloud top
  //----- Begin Code -----
  // Set up mean field variables for <r_r> and <N_r>.
  for (int index_1 = 0; index_1 < rrm.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm.extent_int(0); ++index_0) {
      rrm(index_0, index_1) = hydromet(index_0, index_1, hm_metadata_iirr);
    }
  }
  for (int index_1 = 0; index_1 < Nrm.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Nrm.extent_int(0); ++index_0) {
      Nrm(index_0, index_1) = hydromet(index_0, index_1, hm_metadata_iiNr);
    }
  }
  // Set KK microphysics tendency adjustment flags
  l_src_adj_enabled   = true;
  l_evap_adj_enabled  = true;
  l_clip_positive_sed = true;

  // The time tendency of cloud droplet concentration is only present because
  // of Latin Hypercube interface.  Ncm_mc is needed for other microphysics
  // schemes, but not KK.  Simply set Ncm_mc to 0.
  for (int index_1 = 0; index_1 < Ncm_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Ncm_mc.extent_int(0); ++index_0) {
      Ncm_mc(index_0, index_1) = zero;
    }
  }

  // Do not use these adjustments when l_silhs_KK_convergence_adj_mean is
  // true. We need to do this to get Latin hypercube to converge to KK
  // analytic. The adjustment code will be called by Latin hypercube for the
  // means only.
  if ((microphys.l_silhs_KK_convergence_adj_mean && l_latin_hypercube)) {
    l_src_adj_enabled   = false;
    l_evap_adj_enabled  = false;
    l_clip_positive_sed = false;
  }

  //!! Microphysics tendency loop.
  // Loop over all model thermodynamic level above the model lower boundary.
  for (k = 0; k < nzt; k += 1) {
    for (i = 0; i < ngrdcol; i += 1) {

      //!! Calculate the coefficients for the KK microphysics tendencies.
      microphys_host::KK_tendency_coefs(thlm(i, k), exner(i, k), p_in_Pa(i, k), rho(i, k),
                                        KK_evap_coef, KK_auto_coef, KK_accr_coef, KK_mvr_coef,
                                        saturation_formula, kk);

      //!! Calculate the local KK rain drop mean volume radius.
      if ((rrm(i, k) > rr_tol)) {
        KK_mean_vol_rad(i, k) = kk_engine::KK_mvr_local_mean(rrm(i, k), Nrm(i, k), KK_mvr_coef, kk);
      } else {
        KK_mean_vol_rad(i, k) = zero;
      }

      //!! Calculate the local KK evaporation tendency.
      if (((rrm(i, k) > rr_tol) && (Nrm(i, k) > Nr_tol))) {
        KK_evap_tndcy(i, k) =
            kk_engine::KK_evap_local_mean(chi(i, k), rrm(i, k), Nrm(i, k), KK_evap_coef, kk);
      } else {
        KK_evap_tndcy(i, k) = zero;
      }

      //!! Calculate the local KK autoconversion tendency.
      if ((Ncm(i, k) > Nc_tol)) {
        KK_auto_tndcy(i, k) = kk_engine::KK_auto_local_mean(chi(i, k), Ncm(i, k), KK_auto_coef, kk);
      } else {
        KK_auto_tndcy(i, k) = zero;
      }

      //!! Calculate the local KK accretion tendency.
      if ((rrm(i, k) > rr_tol)) {
        KK_accr_tndcy(i, k) = kk_engine::KK_accr_local_mean(chi(i, k), rrm(i, k), KK_accr_coef, kk);
      } else {
        KK_accr_tndcy(i, k) = zero;
      }

      //!! Calculate the KK rain-number tendencies.
      if (((rrm(i, k) > rr_tol) && (Nrm(i, k) > Nr_tol))) {
        KK_Nrm_evap_tndcy(i, k) =
            kk_engine::KK_Nrm_evap_local_mean(KK_evap_tndcy(i, k), Nrm(i, k), rrm(i, k), dt, kk);
      } else {
        KK_Nrm_evap_tndcy(i, k) = zero;
      }
      KK_Nrm_auto_tndcy(i, k) = kk_engine::KK_Nrm_auto_mean(KK_auto_tndcy(i, k), kk);

      //!! Calculate any necessary adjustments to KK microphysics tendencies.
      microphys_host::KK_microphys_adjust(
          dt, exner(i, k), rcm(i, k), rrm(i, k), Nrm(i, k), KK_evap_tndcy(i, k),
          KK_auto_tndcy(i, k), KK_accr_tndcy(i, k), KK_Nrm_evap_tndcy(i, k),
          KK_Nrm_auto_tndcy(i, k), l_src_adj_enabled, l_evap_adj_enabled, rrm_mc(i, k),
          Nrm_mc(i, k), rvm_mc(i, k), rcm_mc(i, k), thlm_mc(i, k), adj_terms(i, k), kk);
    }

    // Microphysics tendency loop: k = 1, nzt, 1
  }

  //!! Boundary conditions for microphysics tendencies.
  for (i = 0; i < ngrdcol; i += 1) {
    rrm_mc(i, ((nzt - 1)))          = zero;
    Nrm_mc(i, ((nzt - 1)))          = zero;
    KK_mean_vol_rad(i, ((nzt - 1))) = zero;
    rvm_mc(i, ((nzt - 1)))          = zero;
    rcm_mc(i, ((nzt - 1)))          = zero;
    thlm_mc(i, ((nzt - 1)))         = zero;
  }

  // Find the vertical level index of cloud top.
  Kokkos::deep_copy(cloud_top_level,
                    to_host_space(workspace, clubb_kokkos::get_cloud_top_level(
                                                 nzt, ngrdcol, to_execution_space(workspace, rcm),
                                                 to_execution_space(workspace, hydromet),
                                                 hydromet_dim, hm_metadata_iiri)));

  //!! Microphysics sedimentation velocities.
  microphys_host::KK_upscaled_sedimentation(ngrdcol, nzt, cloud_top_level, KK_mean_vol_rad,
                                            l_clip_positive_sed, Vrr, VNr);

  //!! Output hydrometeor mean tendencies and mean sedimentation velocities
  //!! in output arrays.
  microphys_host::KK_microphys_output(ngrdcol, nzt, hydromet_dim, hm_metadata, Vrr, VNr, rrm_mc,
                                      Nrm_mc, hydromet_mc, hydromet_vel);

  for (int index_1 = 0; index_1 < rrm_auto_diag.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm_auto_diag.extent_int(0); ++index_0) {
      rrm_auto_diag(index_0, index_1) = KK_auto_tndcy(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < rrm_accr_diag.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm_accr_diag.extent_int(0); ++index_0) {
      rrm_accr_diag(index_0, index_1) = KK_accr_tndcy(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < rrm_evap_diag.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm_evap_diag.extent_int(0); ++index_0) {
      rrm_evap_diag(index_0, index_1) = KK_evap_tndcy(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < Nrm_auto_diag.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Nrm_auto_diag.extent_int(0); ++index_0) {
      Nrm_auto_diag(index_0, index_1) = KK_Nrm_auto_tndcy(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < Nrm_evap_diag.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Nrm_evap_diag.extent_int(0); ++index_0) {
      Nrm_evap_diag(index_0, index_1) = KK_Nrm_evap_tndcy(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < mvrr_diag.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < mvrr_diag.extent_int(0); ++index_0) {
      mvrr_diag(index_0, index_1) = KK_mean_vol_rad(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < rrm_src_adj_diag.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm_src_adj_diag.extent_int(0); ++index_0) {
      rrm_src_adj_diag(index_0, index_1) = adj_terms(index_0, index_1).rrm_src_adj;
    }
  }
  for (int index_1 = 0; index_1 < Nrm_src_adj_diag.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Nrm_src_adj_diag.extent_int(0); ++index_0) {
      Nrm_src_adj_diag(index_0, index_1) = adj_terms(index_0, index_1).Nrm_src_adj;
    }
  }
  for (int index_1 = 0; index_1 < rrm_evap_adj_diag.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm_evap_adj_diag.extent_int(0); ++index_0) {
      rrm_evap_adj_diag(index_0, index_1) = adj_terms(index_0, index_1).rrm_evap_adj;
    }
  }
  for (int index_1 = 0; index_1 < Nrm_evap_adj_diag.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Nrm_evap_adj_diag.extent_int(0); ++index_0) {
      Nrm_evap_adj_diag(index_0, index_1) = adj_terms(index_0, index_1).Nrm_evap_adj;
    }
  }
  for (int index_1 = 0; index_1 < rrm_mc_nonadj_diag.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm_mc_nonadj_diag.extent_int(0); ++index_0) {
      rrm_mc_nonadj_diag(index_0, index_1) =
          ((KK_auto_tndcy(index_0, index_1) + KK_accr_tndcy(index_0, index_1)) +
           KK_evap_tndcy(index_0, index_1));
    }
  }

  return;
}
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25>
inline void KK_upscaled_microphys(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int nzt, int nzm,
    int pdf_dim, int hydromet_dim, const hm_metadata_type& hm_metadata, Array0 wm_zt, Array1 rtm,
    Array2 thlm, Array3 p_in_Pa, Array4 exner, Array5 rho, Array6 rcm,
    const pdf_parameter& pdf_params, const hydromet_pdf_parameter_collection& hydromet_pdf_params,
    const precipitation_fractions& precip_fracs, Array7 hydromet, Array8 mu_x_1_n, Array9 mu_x_2_n,
    Array10 sigma_x_1_n, Array11 sigma_x_2_n, Array12 corr_array_1_n, Array13 corr_array_2_n,
    int saturation_formula, stats_type& stats, Array14 hydromet_mc, Array15 hydromet_vel,
    Array16 rcm_mc, Array17 rvm_mc, Array18 thlm_mc, Array19 hydromet_vel_covar_zt_impc,
    Array20 hydromet_vel_covar_zt_expc, Array21 wprtp_mc, Array22 wpthlp_mc, Array23 rtp2_mc,
    Array24 thlp2_mc, Array25 rtpthlp_mc, const parameters_microphys& microphys,
    const parameters_KK& kk) {
  // Description:
  // Version of KK microphysics scheme that is analytically upscaled by
  // integrating over the product of the microphysics tendency and the
  // functional form of the PDF.
  // References:
  // Larson, V. E. and B. M. Griffin, 2013:  Analytic upscaling of a local
  //    microphysics scheme. Part I: Derivation.  Q. J. Roy. Meteorol. Soc.,
  //    139, 670, 46--57, doi:http://dx.doi.org/10.1002/qj.1967.
  //
  // Griffin, B. M. and V. E. Larson, 2013:  Analytic upscaling of a local
  //    microphysics scheme. Part II: Simulations.  Q. J. Roy. Meteorol. Soc.,
  //    139, 670, 58--69, doi:http://dx.doi.org/10.1002/qj.1966.
  //
  // Griffin, B. M., 2016:  Improving the Subgrid-Scale Representation of
  //    Hydrometeors and Microphysical Feedback Effects Using a Multivariate
  //    PDF.  Doctoral dissertation, University of Wisconsin -- Milwaukee,
  //    Milwaukee, WI, Paper 1144, 165 pp., URL
  //    http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Griffin, B. M. and V. E. Larson, 2016:  Supplement of A new subgrid-scale
  //    representation of hydrometeor fields using a multivariate PDF.
  //    Geosci. Model Dev., 9, 6,
  //    doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //
  // Griffin, B. M. and V. E. Larson, 2016:  Parameterizing microphysical
  //    effects on variances and covariances of moisture and heat content using
  //    a multivariate probability density function: a study with CLUBB (tag
  //    MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
  //    doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
  // Local Constants
  //-------------------------- Input Variables --------------------------
  //-------------------------- Output Variables --------------------------
  //-------------------------- Local Variables --------------------------
  // In: gr
  // In: ngrdcol — Number of model columns
  // In: dt — Model time step duration                 [s]
  // In: nzt — Number of model thermodynamic vertical grid levels
  // In: nzm — Number of model momentum vertical grid levels
  // In: pdf_dim — Number of variables in the correlation arrays
  // In: hydromet_dim
  // In: hm_metadata
  // In: wm_zt (ngrdcol,nzt) — Mean vertical velocity on thermodynamic levels  [m/s]
  // In: rtm (ngrdcol,nzt) — Mean total water mixing ratio                   [kg/kg]
  // In: thlm (ngrdcol,nzt) — Mean liquid water potential temperature         [K]
  // In: p_in_Pa (ngrdcol,nzt) — Pressure                                        [Pa]
  // In: exner (ngrdcol,nzt) — Exner function                                  [-]
  // In: rho (ngrdcol,nzt) — Density                                         [kg/m^3]
  // In: rcm (ngrdcol,nzt) — Mean cloud water mixing ratio                   [kg/kg]
  // In: pdf_params — PDF parameters                         [units vary]
  // In: hydromet_pdf_params (ngrdcol,nzt)
  // In: precip_fracs — Precipitation fractions      [-]
  // In: hydromet (ngrdcol,nzt,hydromet_dim) — Hydrometeor mean, < h_m > (thermodynamic levels)  [units]
  // In: mu_x_1_n (ngrdcol,nzt,pdf_dim) — Mean array (normal space): PDF vars. (comp. 1) [un. vary]
  // In: mu_x_2_n (ngrdcol,nzt,pdf_dim) — Mean array (normal space): PDF vars. (comp. 2) [un. vary]
  // In: sigma_x_1_n (ngrdcol,nzt,pdf_dim) — Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
  // In: sigma_x_2_n (ngrdcol,nzt,pdf_dim) — Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]
  // In: corr_array_1_n (ngrdcol,nzt,pdf_dim,pdf_dim) — Corr. array (normal space) of PDF vars. (comp. 1)  [-]
  // In: corr_array_2_n (ngrdcol,nzt,pdf_dim,pdf_dim) — Corr. array (normal space) of PDF vars. (comp. 2)  [-]
  // In: saturation_formula — Integer that stores the saturation formula to be used
  // Inout: stats
  // Out: hydromet_mc (ngrdcol,nzt,hydromet_dim) — Hydrometeor time tendency          [(units vary)/s]
  // Out: hydromet_vel (ngrdcol,nzt,hydromet_dim) — Hydrometeor sedimentation velocity [m/s]
  // Out: rcm_mc (ngrdcol,nzt) — Time tendency of liquid water mixing ratio    [kg/kg/s]
  // Out: rvm_mc (ngrdcol,nzt) — Time tendency of vapor water mixing ratio     [kg/kg/s]
  // Out: thlm_mc (ngrdcol,nzt) — Time tendency of liquid potential temperature [K/s]
  // Out: hydromet_vel_covar_zt_impc (ngrdcol,nzt,hydromet_dim) — Imp. comp. <V_hm'h_m'> t-levs [m/s]
  // Out: hydromet_vel_covar_zt_expc (ngrdcol,nzt,hydromet_dim) — Exp. comp. <V_hm'h_m'> t-levs [units(m/s)]
  // Out: wprtp_mc (ngrdcol,nzm) — Microphysics tendency for <w'rt'>   [m*(kg/kg)/s^2]
  // Out: wpthlp_mc (ngrdcol,nzm) — Microphysics tendency for <w'thl'>  [m*K/s^2]
  // Out: rtp2_mc (ngrdcol,nzm) — Microphysics tendency for <rt'^2>   [(kg/kg)^2/s]
  // Out: thlp2_mc (ngrdcol,nzm) — Microphysics tendency for <thl'^2>  [K^2/s]
  // Out: rtpthlp_mc (ngrdcol,nzm) — Microphysics tendency for <rt'thl'> [K*(kg/kg)/s]
  // Case state: microphys replaces the Native parameters_microphys module.
  // Case state: kk replaces the Native parameters_KK module.

  auto workspace_scope                   = workspace.scope();
  const auto hm_metadata_iirr            = hm_metadata.iirr;
  const auto hm_metadata_iiri            = hm_metadata.iiri;
  const auto hm_metadata_iiNr            = hm_metadata.iiNr;
  const auto hm_metadata_iiPDF_chi       = hm_metadata.iiPDF_chi;
  const auto hm_metadata_iiPDF_eta       = hm_metadata.iiPDF_eta;
  const auto hm_metadata_iiPDF_w         = hm_metadata.iiPDF_w;
  const auto hm_metadata_iiPDF_rr        = hm_metadata.iiPDF_rr;
  const auto hm_metadata_iiPDF_Nr        = hm_metadata.iiPDF_Nr;
  const auto hm_metadata_iiPDF_Ncn       = hm_metadata.iiPDF_Ncn;
  const auto pdf_params_rt_1             = to_host_space(workspace, pdf_params.rt_1);
  const auto pdf_params_rt_2             = to_host_space(workspace, pdf_params.rt_2);
  const auto pdf_params_thl_1            = to_host_space(workspace, pdf_params.thl_1);
  const auto pdf_params_thl_2            = to_host_space(workspace, pdf_params.thl_2);
  const auto pdf_params_crt_1            = to_host_space(workspace, pdf_params.crt_1);
  const auto pdf_params_crt_2            = to_host_space(workspace, pdf_params.crt_2);
  const auto pdf_params_cthl_1           = to_host_space(workspace, pdf_params.cthl_1);
  const auto pdf_params_cthl_2           = to_host_space(workspace, pdf_params.cthl_2);
  const auto pdf_params_mixt_frac        = to_host_space(workspace, pdf_params.mixt_frac);
  const auto hydromet_pdf_params_hm_1    = to_host_space(workspace, hydromet_pdf_params.hm_1);
  const auto hydromet_pdf_params_hm_2    = to_host_space(workspace, hydromet_pdf_params.hm_2);
  const auto hydromet_pdf_params_mu_hm_1 = to_host_space(workspace, hydromet_pdf_params.mu_hm_1);
  const auto hydromet_pdf_params_mu_hm_2 = to_host_space(workspace, hydromet_pdf_params.mu_hm_2);
  const auto hydromet_pdf_params_sigma_hm_1 =
      to_host_space(workspace, hydromet_pdf_params.sigma_hm_1);
  const auto hydromet_pdf_params_sigma_hm_2 =
      to_host_space(workspace, hydromet_pdf_params.sigma_hm_2);
  const auto hydromet_pdf_params_mu_Ncn_1 = to_host_space(workspace, hydromet_pdf_params.mu_Ncn_1);
  const auto hydromet_pdf_params_mu_Ncn_2 = to_host_space(workspace, hydromet_pdf_params.mu_Ncn_2);
  const auto hydromet_pdf_params_sigma_Ncn_1 =
      to_host_space(workspace, hydromet_pdf_params.sigma_Ncn_1);
  const auto hydromet_pdf_params_sigma_Ncn_2 =
      to_host_space(workspace, hydromet_pdf_params.sigma_Ncn_2);
  const auto precip_fracs_precip_frac_1 = to_host_space(workspace, precip_fracs.precip_frac_1);
  const auto precip_fracs_precip_frac_2 = to_host_space(workspace, precip_fracs.precip_frac_2);
  constexpr bool l_clip_positive_sed    = true;
  HostReal2D rrm                        = workspace_scope.view<HostReal2D>(
      "KK_upscaled_microphys::rrm", ngrdcol, nzt); // Mean rain water mixing ratio, < r_r > [kg/kg]
  HostReal2D Nrm = workspace_scope.view<HostReal2D>(
      "KK_upscaled_microphys::Nrm", ngrdcol, nzt); // Mean rain drop concentration, < N_r > [num/kg]
  HostReal2D Vrr = workspace_scope.view<HostReal2D>(
      "KK_upscaled_microphys::Vrr", ngrdcol, nzt); // Mean sedimentation velocity of r_r    [m/s]
  HostReal2D VNr = workspace_scope.view<HostReal2D>(
      "KK_upscaled_microphys::VNr", ngrdcol, nzt); // Mean sedimentation velocity of N_r    [m/s]
  HostReal2D rrm_mc =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::rrm_mc", ngrdcol,
                                       nzt); // Mean (dr_r/dt) due to microphysics    [(kg/kg)/s]
  HostReal2D Nrm_mc =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::Nrm_mc", ngrdcol,
                                       nzt); // Mean (dN_r/dt) due to microphysics    [(num/kg)/s]
  HostReal2D KK_evap_tndcy =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::KK_evap_tndcy", ngrdcol,
                                       nzt); // Mean KK (dr_r/dt) due to evaporation     [(kg/kg)/s]
  HostReal2D KK_auto_tndcy =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::KK_auto_tndcy", ngrdcol,
                                       nzt); // Mean KK (dr_r/dt) due to autoconversion  [(kg/kg)/s]
  HostReal2D KK_accr_tndcy =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::KK_accr_tndcy", ngrdcol,
                                       nzt); // Mean KK (dr_r/dt) due to accretion       [(kg/kg)/s]
  HostReal2D KK_mean_vol_rad =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::KK_mean_vol_rad", ngrdcol,
                                       nzt); // Mean KK rain drop mean volume radius     [m]
  HostReal2D KK_Nrm_evap_tndcy =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::KK_Nrm_evap_tndcy", ngrdcol,
                                       nzt); // Mean KK (dN_r/dt) due to evaporation  [(num/kg)/s]
  HostReal2D KK_Nrm_auto_tndcy =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::KK_Nrm_auto_tndcy", ngrdcol,
                                       nzt); // Mean KK (dN_r/dt) due to autoconv.    [(num/kg)/s]
  real_t KK_evap_coef; // KK evaporation coefficient                  [(kg/kg)/s]
  real_t KK_auto_coef; // KK autoconversion coefficient               [(kg/kg)/s]
  real_t KK_accr_coef; // KK accretion coefficient                    [(kg/kg)/s]
  real_t KK_mvr_coef;  // KK mean volume radius coefficient           [m]
  HostReal2D mixt_frac =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::mixt_frac", ngrdcol,
                                       nzt); // Mixture fraction                                [-]
  real_t mu_w_1;        // Mean of w (1st PDF component)                     [m/s]
  real_t mu_w_2;        // Mean of w (2nd PDF component)                     [m/s]
  real_t mu_chi_1;      // Mean of chi (1st PDF component)                 [kg/kg]
  real_t mu_chi_2;      // Mean of chi (2nd PDF component)                 [kg/kg]
  real_t mu_eta_1;      // Mean of eta (1st PDF component)                 [kg/kg]
  real_t mu_eta_2;      // Mean of eta (2nd PDF component)                 [kg/kg]
  real_t mu_rr_1;       // Mean of rr (1st PDF component) in-precip (ip)   [kg/kg]
  real_t mu_rr_2;       // Mean of rr (2nd PDF component) ip               [kg/kg]
  real_t mu_Nr_1;       // Mean of Nr (1st PDF component) ip              [num/kg]
  real_t mu_Nr_2;       // Mean of Nr (2nd PDF component) ip              [num/kg]
  real_t mu_Ncn_1;      // Mean of Ncn (1st PDF component)                [num/kg]
  real_t mu_Ncn_2;      // Mean of Ncn (2nd PDF component)                [num/kg]
  real_t mu_rr_1_n;     // Mean of ln rr (1st PDF component) ip        [ln(kg/kg)]
  real_t mu_rr_2_n;     // Mean of ln rr (2nd PDF component) ip        [ln(kg/kg)]
  real_t mu_Nr_1_n;     // Mean of ln Nr (1st PDF component) ip       [ln(num/kg)]
  real_t mu_Nr_2_n;     // Mean of ln Nr (2nd PDF component) ip       [ln(num/kg)]
  real_t mu_Ncn_1_n;    // Mean of ln Ncn (1st PDF component)         [ln(num/kg)]
  real_t mu_Ncn_2_n;    // Mean of ln Ncn (2nd PDF component)         [ln(num/kg)]
  real_t sigma_w_1;     // Standard deviation of w (1st PDF component)       [m/s]
  real_t sigma_w_2;     // Standard deviation of w (2nd PDF component)       [m/s]
  real_t sigma_chi_1;   // Standard deviation of chi (1st PDF component)   [kg/kg]
  real_t sigma_chi_2;   // Standard deviation of chi (2nd PDF component)   [kg/kg]
  real_t sigma_eta_1;   // Standard deviation of eta (1st PDF component)   [kg/kg]
  real_t sigma_eta_2;   // Standard deviation of eta (2nd PDF component)   [kg/kg]
  real_t sigma_rr_1;    // Standard deviation of rr (1st PDF component) ip [kg/kg]
  real_t sigma_rr_2;    // Standard deviation of rr (2nd PDF component) ip [kg/kg]
  real_t sigma_Nr_1;    // Standard deviation of Nr (1st PDF comp.) ip    [num/kg]
  real_t sigma_Nr_2;    // Standard deviation of Nr (2nd PDF comp.) ip    [num/kg]
  real_t sigma_Ncn_1;   // Standard deviation of Ncn (1st PDF component)  [num/kg]
  real_t sigma_Ncn_2;   // Standard deviation of Ncn (2nd PDF component)  [num/kg]
  real_t sigma_rr_1_n;  // Standard deviation of ln rr (1st PDF component) ip  [-]
  real_t sigma_rr_2_n;  // Standard deviation of ln rr (2nd PDF component) ip  [-]
  real_t sigma_Nr_1_n;  // Standard deviation of ln Nr (1st PDF component) ip  [-]
  real_t sigma_Nr_2_n;  // Standard deviation of ln Nr (2nd PDF component) ip  [-]
  real_t sigma_Ncn_1_n; // Standard deviation of ln Ncn (1st PDF component)    [-]
  real_t sigma_Ncn_2_n; // Standard deviation of ln Ncn (2nd PDF component)    [-]
  real_t corr_w_chi_1;  // Correlation of w and chi (1st PDF component)     [-]
  real_t
      corr_w_chi_2; // Correlation of w and chi (2nd PDF component)     [-] corr_w_rr_1,    & ! Correlation of w and rr (1st PDF component) ip   [-] corr_w_rr_2,    & ! Correlation of w and rr (2nd PDF component) ip   [-] corr_w_Nr_1,    & ! Correlation of w and Nr (1st PDF component) ip   [-] corr_w_Nr_2,    & ! Correlation of w and Nr (2nd PDF component) ip   [-] corr_w_Ncn_1,   & ! Correlation of w and Ncn (1st PDF component)     [-] corr_w_Ncn_2,   & ! Correlation of w and Ncn (2nd PDF component)     [-]
  real_t corr_chi_eta_1; // Correlation of chi and eta (1st PDF component)   [-]
  real_t
      corr_chi_eta_2; // ,& ! Correlation of chi and eta (2nd PDF component)   [-] corr_chi_rr_1,  & ! Correlation of chi and rr (1st PDF component) ip [-] corr_chi_rr_2,  & ! Correlation of chi and rr (2nd PDF component) ip [-] corr_chi_Nr_1,  & ! Correlation of chi and Nr (1st PDF component) ip [-] corr_chi_Nr_2,  & ! Correlation of chi and Nr (2nd PDF component) ip [-] corr_chi_Ncn_1, & ! Correlation of chi and Ncn (1st PDF component)   [-] corr_chi_Ncn_2, & ! Correlation of chi and Ncn (2nd PDF component)   [-] corr_eta_rr_1,  & ! Correlation of eta and rr (1st PDF component) ip [-] corr_eta_rr_2,  & ! Correlation of eta and rr (2nd PDF component) ip [-] corr_eta_Nr_1,  & ! Correlation of eta and Nr (1st PDF component) ip [-] corr_eta_Nr_2,  & ! Correlation of eta and Nr (2nd PDF component) ip [-] corr_eta_Ncn_1, & ! Correlation of eta and Ncn (1st PDF component)   [-] corr_eta_Ncn_2, & ! Correlation of eta and Ncn (2nd PDF component)   [-] corr_rr_Nr_1,   & ! Correlation of rr and Nr (1st PDF component) ip  [-] corr_rr_Nr_2      ! Correlation of rr and Nr (2nd PDF component) ip  [-]
  real_t corr_w_rr_1_n;    // Correlation of w and ln rr (1st PDF comp.) ip    [-]
  real_t corr_w_rr_2_n;    // Correlation of w and ln rr (2nd PDF comp.) ip    [-]
  real_t corr_w_Nr_1_n;    // Correlation of w and ln Nr (1st PDF comp.) ip    [-]
  real_t corr_w_Nr_2_n;    // Correlation of w and ln Nr (2nd PDF comp.) ip    [-]
  real_t corr_w_Ncn_1_n;   // Correlation of w and ln Ncn (1st PDF component)  [-]
  real_t corr_w_Ncn_2_n;   // Correlation of w and ln Ncn (2nd PDF component)  [-]
  real_t corr_chi_rr_1_n;  // Correlation of chi and ln rr (1st PDF comp.) ip  [-]
  real_t corr_chi_rr_2_n;  // Correlation of chi and ln rr (2nd PDF comp.) ip  [-]
  real_t corr_chi_Nr_1_n;  // Correlation of chi and ln Nr (1st PDF comp.) ip  [-]
  real_t corr_chi_Nr_2_n;  // Correlation of chi and ln Nr (2nd PDF comp.) ip  [-]
  real_t corr_chi_Ncn_1_n; // Correlation of chi and ln Ncn (1st PDF comp.)    [-]
  real_t corr_chi_Ncn_2_n; // Correlation of chi and ln Ncn (2nd PDF comp.)    [-]
  real_t corr_eta_rr_1_n;  // Correlation of eta and ln rr (1st PDF comp.) ip  [-]
  real_t corr_eta_rr_2_n;  // Correlation of eta and ln rr (2nd PDF comp.) ip  [-]
  real_t corr_eta_Nr_1_n;  // Correlation of eta and ln Nr (1st PDF comp.) ip  [-]
  real_t corr_eta_Nr_2_n;  // Correlation of eta and ln Nr (2nd PDF comp.) ip  [-]
  real_t corr_eta_Ncn_1_n; // Correlation of eta and ln Ncn (1st PDF comp.)    [-]
  real_t corr_eta_Ncn_2_n; // Correlation of eta and ln Ncn (2nd PDF comp.)    [-]
  real_t corr_rr_Nr_1_n;   // Correlation of ln rr & ln Nr (1st PDF comp.) ip  [-]
  real_t corr_rr_Nr_2_n;   // Correlation of ln rr & ln Nr (2nd PDF comp.) ip  [-]
  real_t rr_1;             // Mean rain water mixing ratio (1st PDF component)      [kg/kg]
  real_t rr_2;             // Mean rain water mixing ratio (2nd PDF component)      [kg/kg]
  real_t Nr_1;             // Mean rain drop concentration (1st PDF component)      [num/kg]
  real_t Nr_2;             // Mean rain drop concentration (2nd PDF component)      [num/kg]
  real_t precip_frac_1;    // Precipitation fraction (1st PDF component) [-]
  real_t precip_frac_2;    // Precipitation fraction (2nd PDF component) [-]
  HostReal2D Vrrprrp_zt_impc =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::Vrrprrp_zt_impc", ngrdcol,
                                       nzt); // Imp. comp. of <V_rr'r_r'>: <r_r> eq.  [(m/s)]
  HostReal2D Vrrprrp_zt_expc =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::Vrrprrp_zt_expc", ngrdcol,
                                       nzt); // Exp. comp. of <V_rr'r_r'>: <r_r> eq.  [(m/s)(kg/kg)]
  HostReal2D VNrpNrp_zt_impc =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::VNrpNrp_zt_impc", ngrdcol,
                                       nzt); // Imp. comp. of <V_Nr'N_r'>: <N_r> eq.  [(m/s)]
  HostReal2D VNrpNrp_zt_expc = workspace_scope.view<HostReal2D>(
      "KK_upscaled_microphys::VNrpNrp_zt_expc", ngrdcol,
      nzt); // Exp. comp. of <V_Nr'N_r'>: <N_r> eq.  [(m/s)(num/kg)]
  HostReal2D rr_KK_mvr_covar =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::rr_KK_mvr_covar", ngrdcol,
                                       nzt); // Covariance of r_r and KK mean vol. radius [(kg/kg)m]
  HostReal2D Nr_KK_mvr_covar = workspace_scope.view<HostReal2D>(
      "KK_upscaled_microphys::Nr_KK_mvr_covar", ngrdcol,
      nzt); // Covariance of N_r and KK mean vol. radius [(num/kg)m]
  HostReal2D KK_mvr_variance =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::KK_mvr_variance", ngrdcol,
                                       nzt); // Variance of KK rain drop mean vol. radius [m^2]
  HostReal2D w_KK_evap_covar =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::w_KK_evap_covar", ngrdcol,
                                       nzt); // Covar. of w and KK evap. tend.    [m*(kg/kg)/s^2]
  HostReal2D rt_KK_evap_covar =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::rt_KK_evap_covar", ngrdcol,
                                       nzt); // Covar. of rt and KK evap. tend.   [(kg/kg)^2/s]
  HostReal2D thl_KK_evap_covar =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::thl_KK_evap_covar", ngrdcol,
                                       nzt); // Covar. of thl and KK evap. tend.  [K*(kg/kg)/s]
  HostReal2D w_KK_auto_covar =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::w_KK_auto_covar", ngrdcol,
                                       nzt); // Covar. of w and KK auto. tend.    [m*(kg/kg)/s^2]
  HostReal2D rt_KK_auto_covar =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::rt_KK_auto_covar", ngrdcol,
                                       nzt); // Covar. of rt and KK auto. tend.   [(kg/kg)^2/s]
  HostReal2D thl_KK_auto_covar =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::thl_KK_auto_covar", ngrdcol,
                                       nzt); // Covar. of thl and KK auto. tend.  [K*(kg/kg)/s]
  HostReal2D w_KK_accr_covar =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::w_KK_accr_covar", ngrdcol,
                                       nzt); // Covar. of w and KK accr. tend.    [m*(kg/kg)/s^2]
  HostReal2D rt_KK_accr_covar =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::rt_KK_accr_covar", ngrdcol,
                                       nzt); // Covar. of rt and KK accr. tend.   [(kg/kg)^2/s]
  HostReal2D thl_KK_accr_covar =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::thl_KK_accr_covar", ngrdcol,
                                       nzt); // Covar. of thl and KK accr. tend.  [K*(kg/kg)/s]
  HostReal2D wprtp_mc_zt =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::wprtp_mc_zt", ngrdcol,
                                       nzt); // Micro. tend. for <w'rt'>; t-lev   [m*(kg/kg)/s^2]
  HostReal2D wpthlp_mc_zt =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::wpthlp_mc_zt", ngrdcol,
                                       nzt); // Micro. tend. for <w'thl'>; t-lev  [m*K/s^2]
  HostReal2D rtp2_mc_zt =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::rtp2_mc_zt", ngrdcol,
                                       nzt); // Micro. tend. for <rt'^2>; t-lev   [(kg/kg)^2/s]
  HostReal2D thlp2_mc_zt =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::thlp2_mc_zt", ngrdcol,
                                       nzt); // Micro. tend. for <thl'^2>; t-lev  [K^2/s]
  HostReal2D rtpthlp_mc_zt =
      workspace_scope.view<HostReal2D>("KK_upscaled_microphys::rtpthlp_mc_zt", ngrdcol,
                                       nzt); // Micro. tend. for <rt'thl'>; t-lev [K*(kg/kg)/s]
  HostAdj2D adj_terms("KK_upscaled_microphys::adj_terms", ngrdcol,
                      nzt); // Adjustments to microphysics terms
  bool l_src_adj_enabled;   // Flag to enable rrm/Nrm source adjustment
  bool l_evap_adj_enabled;  // Flag to enable rrm/Nrm evaporation adjustment
  HostInt1D cloud_top_level = workspace_scope.view<HostInt1D>(
      "KK_upscaled_microphys::cloud_top_level", ngrdcol); // Vertical level index of cloud top
  int i;                                                  // Column loop index
  int k;

  //-------------------------- Begin Code --------------------------

  for (int index_1 = 0; index_1 < rrm.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm.extent_int(0); ++index_0) {
      rrm(index_0, index_1) = hydromet(index_0, index_1, hm_metadata_iirr);
    }
  }
  for (int index_1 = 0; index_1 < Nrm.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Nrm.extent_int(0); ++index_0) {
      Nrm(index_0, index_1) = hydromet(index_0, index_1, hm_metadata_iiNr);
    }
  }
  l_src_adj_enabled  = true;
  l_evap_adj_enabled = true;

  // Setup mixture fraction.
  for (int index_1 = 0; index_1 < mixt_frac.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < mixt_frac.extent_int(0); ++index_0) {
      mixt_frac(index_0, index_1) = pdf_params_mixt_frac(index_0, index_1);
    }
  }

  //!! Microphysics tendency loop.
  // Loop over all model thermodynamic level above the model lower boundary.
  for (k = 0; k < nzt; ++k) {
    for (i = 0; i < ngrdcol; ++i) {

      //!! Calculate the coefficients for the KK microphysics tendencies.
      microphys_host::KK_tendency_coefs(thlm(i, k), exner(i, k), p_in_Pa(i, k), rho(i, k),
                                        KK_evap_coef, KK_auto_coef, KK_accr_coef, KK_mvr_coef,
                                        saturation_formula, kk);

      //!! Unpack the PDF parameters directly from the multicolumn arrays.
      // Unpack mu_x_i and sigma_x_i into Means and Standard Deviations.
      mu_w_1        = mu_x_1_n(i, k, hm_metadata_iiPDF_w);
      mu_w_2        = mu_x_2_n(i, k, hm_metadata_iiPDF_w);
      mu_chi_1      = mu_x_1_n(i, k, hm_metadata_iiPDF_chi);
      mu_chi_2      = mu_x_2_n(i, k, hm_metadata_iiPDF_chi);
      mu_eta_1      = mu_x_1_n(i, k, hm_metadata_iiPDF_eta);
      mu_eta_2      = mu_x_2_n(i, k, hm_metadata_iiPDF_eta);
      mu_rr_1_n     = mu_x_1_n(i, k, hm_metadata_iiPDF_rr);
      mu_rr_2_n     = mu_x_2_n(i, k, hm_metadata_iiPDF_rr);
      mu_Nr_1_n     = mu_x_1_n(i, k, hm_metadata_iiPDF_Nr);
      mu_Nr_2_n     = mu_x_2_n(i, k, hm_metadata_iiPDF_Nr);
      mu_Ncn_1_n    = mu_x_1_n(i, k, hm_metadata_iiPDF_Ncn);
      mu_Ncn_2_n    = mu_x_2_n(i, k, hm_metadata_iiPDF_Ncn);
      sigma_w_1     = sigma_x_1_n(i, k, hm_metadata_iiPDF_w);
      sigma_w_2     = sigma_x_2_n(i, k, hm_metadata_iiPDF_w);
      sigma_chi_1   = sigma_x_1_n(i, k, hm_metadata_iiPDF_chi);
      sigma_chi_2   = sigma_x_2_n(i, k, hm_metadata_iiPDF_chi);
      sigma_eta_1   = sigma_x_1_n(i, k, hm_metadata_iiPDF_eta);
      sigma_eta_2   = sigma_x_2_n(i, k, hm_metadata_iiPDF_eta);
      sigma_rr_1_n  = sigma_x_1_n(i, k, hm_metadata_iiPDF_rr);
      sigma_rr_2_n  = sigma_x_2_n(i, k, hm_metadata_iiPDF_rr);
      sigma_Nr_1_n  = sigma_x_1_n(i, k, hm_metadata_iiPDF_Nr);
      sigma_Nr_2_n  = sigma_x_2_n(i, k, hm_metadata_iiPDF_Nr);
      sigma_Ncn_1_n = sigma_x_1_n(i, k, hm_metadata_iiPDF_Ncn);
      sigma_Ncn_2_n = sigma_x_2_n(i, k, hm_metadata_iiPDF_Ncn);

      // Unpack variables from hydromet_pdf_params
      mu_rr_1     = hydromet_pdf_params_mu_hm_1(i, k, hm_metadata_iirr);
      mu_rr_2     = hydromet_pdf_params_mu_hm_2(i, k, hm_metadata_iirr);
      mu_Nr_1     = hydromet_pdf_params_mu_hm_1(i, k, hm_metadata_iiNr);
      mu_Nr_2     = hydromet_pdf_params_mu_hm_2(i, k, hm_metadata_iiNr);
      mu_Ncn_1    = hydromet_pdf_params_mu_Ncn_1(i, k);
      mu_Ncn_2    = hydromet_pdf_params_mu_Ncn_2(i, k);
      sigma_rr_1  = hydromet_pdf_params_sigma_hm_1(i, k, hm_metadata_iirr);
      sigma_rr_2  = hydromet_pdf_params_sigma_hm_2(i, k, hm_metadata_iirr);
      sigma_Nr_1  = hydromet_pdf_params_sigma_hm_1(i, k, hm_metadata_iiNr);
      sigma_Nr_2  = hydromet_pdf_params_sigma_hm_2(i, k, hm_metadata_iiNr);
      sigma_Ncn_1 = hydromet_pdf_params_sigma_Ncn_1(i, k);
      sigma_Ncn_2 = hydromet_pdf_params_sigma_Ncn_2(i, k);
      rr_1        = hydromet_pdf_params_hm_1(i, k, hm_metadata_iirr);
      rr_2        = hydromet_pdf_params_hm_2(i, k, hm_metadata_iirr);
      Nr_1        = hydromet_pdf_params_hm_1(i, k, hm_metadata_iiNr);
      Nr_2        = hydromet_pdf_params_hm_2(i, k, hm_metadata_iiNr);

      // Unpack corr_array_1_n into correlations (1st PDF component).
      corr_chi_eta_1   = corr_array_1_n(i, k, hm_metadata_iiPDF_eta, hm_metadata_iiPDF_chi);
      corr_w_chi_1     = corr_array_1_n(i, k, hm_metadata_iiPDF_w, hm_metadata_iiPDF_chi);
      corr_chi_rr_1_n  = corr_array_1_n(i, k, hm_metadata_iiPDF_rr, hm_metadata_iiPDF_chi);
      corr_chi_Nr_1_n  = corr_array_1_n(i, k, hm_metadata_iiPDF_Nr, hm_metadata_iiPDF_chi);
      corr_chi_Ncn_1_n = corr_array_1_n(i, k, hm_metadata_iiPDF_Ncn, hm_metadata_iiPDF_chi);
      corr_eta_rr_1_n  = corr_array_1_n(i, k, hm_metadata_iiPDF_rr, hm_metadata_iiPDF_eta);
      corr_eta_Nr_1_n  = corr_array_1_n(i, k, hm_metadata_iiPDF_Nr, hm_metadata_iiPDF_eta);
      corr_eta_Ncn_1_n = corr_array_1_n(i, k, hm_metadata_iiPDF_Ncn, hm_metadata_iiPDF_eta);
      corr_w_rr_1_n    = corr_array_1_n(i, k, hm_metadata_iiPDF_rr, hm_metadata_iiPDF_w);
      corr_w_Nr_1_n    = corr_array_1_n(i, k, hm_metadata_iiPDF_Nr, hm_metadata_iiPDF_w);
      corr_w_Ncn_1_n   = corr_array_1_n(i, k, hm_metadata_iiPDF_Ncn, hm_metadata_iiPDF_w);
      corr_rr_Nr_1_n   = corr_array_1_n(i, k, hm_metadata_iiPDF_Nr, hm_metadata_iiPDF_rr);

      // Unpack corr_array_2_n into correlations (2nd PDF component).
      corr_chi_eta_2   = corr_array_2_n(i, k, hm_metadata_iiPDF_eta, hm_metadata_iiPDF_chi);
      corr_w_chi_2     = corr_array_2_n(i, k, hm_metadata_iiPDF_w, hm_metadata_iiPDF_chi);
      corr_chi_rr_2_n  = corr_array_2_n(i, k, hm_metadata_iiPDF_rr, hm_metadata_iiPDF_chi);
      corr_chi_Nr_2_n  = corr_array_2_n(i, k, hm_metadata_iiPDF_Nr, hm_metadata_iiPDF_chi);
      corr_chi_Ncn_2_n = corr_array_2_n(i, k, hm_metadata_iiPDF_Ncn, hm_metadata_iiPDF_chi);
      corr_eta_rr_2_n  = corr_array_2_n(i, k, hm_metadata_iiPDF_rr, hm_metadata_iiPDF_eta);
      corr_eta_Nr_2_n  = corr_array_2_n(i, k, hm_metadata_iiPDF_Nr, hm_metadata_iiPDF_eta);
      corr_eta_Ncn_2_n = corr_array_2_n(i, k, hm_metadata_iiPDF_Ncn, hm_metadata_iiPDF_eta);
      corr_w_rr_2_n    = corr_array_2_n(i, k, hm_metadata_iiPDF_rr, hm_metadata_iiPDF_w);
      corr_w_Nr_2_n    = corr_array_2_n(i, k, hm_metadata_iiPDF_Nr, hm_metadata_iiPDF_w);
      corr_w_Ncn_2_n   = corr_array_2_n(i, k, hm_metadata_iiPDF_Ncn, hm_metadata_iiPDF_w);
      corr_rr_Nr_2_n   = corr_array_2_n(i, k, hm_metadata_iiPDF_Nr, hm_metadata_iiPDF_rr);

      precip_frac_1 = precip_fracs_precip_frac_1(i, k);
      precip_frac_2 = precip_fracs_precip_frac_2(i, k);

      //!! Calculate the values of the upscaled KK microphysics tendencies.
      kk_engine::KK_upscaled_means_driver(
          mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_Ncn_1, mu_Ncn_2, mu_rr_1_n,
          mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, sigma_chi_1, sigma_chi_2,
          sigma_rr_1, sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, sigma_rr_1_n,
          sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, corr_chi_rr_1_n,
          corr_chi_rr_2_n, corr_chi_Nr_1_n, corr_chi_Nr_2_n, corr_chi_Ncn_1_n, corr_chi_Ncn_2_n,
          corr_rr_Nr_1_n, corr_rr_Nr_2_n, mixt_frac(i, k), precip_frac_1, precip_frac_2,
          KK_evap_coef, KK_auto_coef, KK_accr_coef, KK_mvr_coef, KK_evap_tndcy(i, k),
          KK_auto_tndcy(i, k), KK_accr_tndcy(i, k), KK_mean_vol_rad(i, k), kk);

      //!! Calculate the values of the turbulent sedimentation terms,
      //!! <V_rr'rr'> and <V_Nr'Nr'>.  Each turbulent sedimentation term has an
      //!! implicit component and an explicit component to be used in the
      //!! predictive equation for the hydrometeor.
      kk_engine::KK_sed_vel_covars(
          rrm(i, k), rr_1, rr_2, Nrm(i, k), Nr_1, Nr_2, KK_mean_vol_rad(i, k), mu_rr_1, mu_rr_2,
          mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1, sigma_rr_2,
          sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n,
          corr_rr_Nr_1_n, corr_rr_Nr_2_n, KK_mvr_coef, mixt_frac(i, k), precip_frac_1,
          precip_frac_2, Vrrprrp_zt_impc(i, k), Vrrprrp_zt_expc(i, k), VNrpNrp_zt_impc(i, k),
          VNrpNrp_zt_expc(i, k), rr_KK_mvr_covar(i, k), Nr_KK_mvr_covar(i, k), kk);

      if (microphys.l_var_covar_src) {

        kk_engine::KK_upscaled_covar_driver(
            wm_zt(i, k), rtm(i, k), thlm(i, k), exner(i, k), rrm(i, k), Nrm(i, k), mu_w_1, mu_w_2,
            mu_chi_1, mu_chi_2, mu_eta_1, mu_eta_2, mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_Ncn_1,
            mu_Ncn_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, sigma_w_1,
            sigma_w_2, sigma_chi_1, sigma_chi_2, sigma_eta_1, sigma_eta_2, sigma_rr_1, sigma_rr_2,
            sigma_Nr_1, sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, sigma_rr_1_n, sigma_rr_2_n,
            sigma_Nr_1_n, sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, corr_w_chi_1, corr_w_chi_2,
            corr_w_rr_1_n, corr_w_rr_2_n, corr_w_Nr_1_n, corr_w_Nr_2_n, corr_w_Ncn_1_n,
            corr_w_Ncn_2_n, corr_chi_eta_1, corr_chi_eta_2, corr_chi_rr_1_n, corr_chi_rr_2_n,
            corr_chi_Nr_1_n, corr_chi_Nr_2_n, corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, corr_eta_rr_1_n,
            corr_eta_rr_2_n, corr_eta_Nr_1_n, corr_eta_Nr_2_n, corr_eta_Ncn_1_n, corr_eta_Ncn_2_n,
            corr_rr_Nr_1_n, corr_rr_Nr_2_n, mixt_frac(i, k), precip_frac_1, precip_frac_2,
            KK_evap_coef, KK_auto_coef, KK_accr_coef, KK_evap_tndcy(i, k), KK_auto_tndcy(i, k),
            KK_accr_tndcy(i, k), pdf_params_rt_1(i, k), pdf_params_rt_2(i, k),
            pdf_params_thl_1(i, k), pdf_params_thl_2(i, k), pdf_params_crt_1(i, k),
            pdf_params_crt_2(i, k), pdf_params_cthl_1(i, k), pdf_params_cthl_2(i, k),
            wprtp_mc_zt(i, k), wpthlp_mc_zt(i, k), rtp2_mc_zt(i, k), thlp2_mc_zt(i, k),
            rtpthlp_mc_zt(i, k), w_KK_evap_covar(i, k), rt_KK_evap_covar(i, k),
            thl_KK_evap_covar(i, k), w_KK_auto_covar(i, k), rt_KK_auto_covar(i, k),
            thl_KK_auto_covar(i, k), w_KK_accr_covar(i, k), rt_KK_accr_covar(i, k),
            thl_KK_accr_covar(i, k), kk);
      }

      //!! KK rain drop concentration microphysics tendencies.

      //!! Calculate the KK N_r evaporation tendency.
      KK_Nrm_evap_tndcy(i, k) = kk_engine::KK_Nrm_evap_upscaled_mean(
          mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n,
          mu_Nr_2_n, sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2, sigma_Nr_1, sigma_Nr_2,
          sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, corr_chi_rr_1_n, corr_chi_rr_2_n,
          corr_chi_Nr_1_n, corr_chi_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, KK_evap_coef,
          mixt_frac(i, k), precip_frac_1, precip_frac_2, dt, kk);

      //!! Calculate the KK N_r autoconversion tendency.
      KK_Nrm_auto_tndcy(i, k) = kk_engine::KK_Nrm_auto_mean(KK_auto_tndcy(i, k), kk);

      //!! Calculate any necessary adjustments to KK microphysics tendencies.
      microphys_host::KK_microphys_adjust(
          dt, exner(i, k), rcm(i, k), rrm(i, k), Nrm(i, k), KK_evap_tndcy(i, k),
          KK_auto_tndcy(i, k), KK_accr_tndcy(i, k), KK_Nrm_evap_tndcy(i, k),
          KK_Nrm_auto_tndcy(i, k), l_src_adj_enabled, l_evap_adj_enabled, rrm_mc(i, k),
          Nrm_mc(i, k), rvm_mc(i, k), rcm_mc(i, k), thlm_mc(i, k), adj_terms(i, k), kk);

      if ((stats.l_sample && var_on_stats_list(stats, std::string("KK_mvr_variance_zt")))) {
        // Calculate the variance of KK rain drop mean volume radius,
        // < R_vr'^2 >.
        KK_mvr_variance(i, k) = kk_engine::variance_KK_mvr(
            mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n,
            sigma_rr_1, sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n,
            sigma_Nr_1_n, sigma_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, KK_mean_vol_rad(i, k),
            KK_mvr_coef, mixt_frac(i, k), precip_frac_1, precip_frac_2, kk);
      }
    }
    // Microphysics tendency loops
  }

  // Statistics
  if (stats.l_sample) {
    // Covariance of r_r and KK rain drop mean volume radius.
    stats_update(std::string("rr_KK_mvr_covar_zt"), rr_KK_mvr_covar, stats);
    // Covariance of N_r and KK rain drop mean volume radius.
    stats_update(std::string("Nr_KK_mvr_covar_zt"), Nr_KK_mvr_covar, stats);
    if (microphys.l_var_covar_src) {
      // All of these covariance variables are being calculated on thermodynamic
      // grid levels (all inputs are on thermodynamic grid levels, so the output
      // is also on thermodynamic grid levels).  These covariances will be
      // combined in various ways to produce the microphysics tendency terms for
      // various model predictive variances and covariances.  These source
      // tendency terms will be interpolated to momentum grid levels.
      // Covariance of w and KK evaporation tendency.
      stats_update(std::string("w_KK_evap_covar_zt"), w_KK_evap_covar, stats);
      // Covariance of r_t and KK evaporation tendency.
      stats_update(std::string("rt_KK_evap_covar_zt"), rt_KK_evap_covar, stats);
      // Covariance of theta_l and KK evaporation tendency.
      stats_update(std::string("thl_KK_evap_covar_zt"), thl_KK_evap_covar, stats);
      // Covariance of w and KK autoconversion tendency.
      stats_update(std::string("w_KK_auto_covar_zt"), w_KK_auto_covar, stats);
      // Covariance of r_t and KK autoconversion tendency.
      stats_update(std::string("rt_KK_auto_covar_zt"), rt_KK_auto_covar, stats);
      // Covariance of theta_l and KK autoconversion tendency.
      stats_update(std::string("thl_KK_auto_covar_zt"), thl_KK_auto_covar, stats);
      // Covariance of w and KK accretion tendency.
      stats_update(std::string("w_KK_accr_covar_zt"), w_KK_accr_covar, stats);
      // Covariance of r_t and KK accretion tendency.
      stats_update(std::string("rt_KK_accr_covar_zt"), rt_KK_accr_covar, stats);
      // Covariance of theta_l and KK accretion tendency.
      stats_update(std::string("thl_KK_accr_covar_zt"), thl_KK_accr_covar, stats);
    }
    stats_update(std::string("rrm_evap"), KK_evap_tndcy, stats);
    stats_update(std::string("rrm_auto"), KK_auto_tndcy, stats);
    stats_update(std::string("rrm_accr"), KK_accr_tndcy, stats);
    stats_update(std::string("mvrr"), KK_mean_vol_rad, stats);
    stats_update(std::string("Nrm_evap"), KK_Nrm_evap_tndcy, stats);
    stats_update(std::string("Nrm_auto"), KK_Nrm_auto_tndcy, stats);
    if (var_on_stats_list(stats, std::string("KK_mvr_variance_zt"))) {
      stats_update(std::string("KK_mvr_variance_zt"), KK_mvr_variance, stats);
    }
  }

  if (microphys.l_var_covar_src) {

    // Output microphysics tendency terms for
    // model variances and covariances on momentum levels.
    zt2zm_api(nzm, nzt, ngrdcol, gr, wprtp_mc_zt, wprtp_mc);
    zt2zm_api(nzm, nzt, ngrdcol, gr, wpthlp_mc_zt, wpthlp_mc);
    zt2zm_api(nzm, nzt, ngrdcol, gr, rtp2_mc_zt, rtp2_mc);
    zt2zm_api(nzm, nzt, ngrdcol, gr, thlp2_mc_zt, thlp2_mc);
    zt2zm_api(nzm, nzt, ngrdcol, gr, rtpthlp_mc_zt, rtpthlp_mc);

    // This section has not been fully GPUized yet, but these variables have GPU copies.
    // Copy the GPU interpolation results back before setting boundaries on the CPU.

    // Set values of microphysics tendency terms to 0 at model lower boundary.
    for (int index_0 = 0; index_0 < wprtp_mc.extent_int(0); ++index_0) {
      wprtp_mc(index_0, 0) = zero;
    }
    for (int index_0 = 0; index_0 < wpthlp_mc.extent_int(0); ++index_0) {
      wpthlp_mc(index_0, 0) = zero;
    }
    for (int index_0 = 0; index_0 < rtp2_mc.extent_int(0); ++index_0) {
      rtp2_mc(index_0, 0) = zero;
    }
    for (int index_0 = 0; index_0 < thlp2_mc.extent_int(0); ++index_0) {
      thlp2_mc(index_0, 0) = zero;
    }
    for (int index_0 = 0; index_0 < rtpthlp_mc.extent_int(0); ++index_0) {
      rtpthlp_mc(index_0, 0) = zero;
    }
    // Set values of microphysics tendency terms to 0 at model upper boundary.
    for (int index_0 = 0; index_0 < wprtp_mc.extent_int(0); ++index_0) {
      wprtp_mc(index_0, ((nzm - 1))) = zero;
    }
    for (int index_0 = 0; index_0 < wpthlp_mc.extent_int(0); ++index_0) {
      wpthlp_mc(index_0, ((nzm - 1))) = zero;
    }
    for (int index_0 = 0; index_0 < rtp2_mc.extent_int(0); ++index_0) {
      rtp2_mc(index_0, ((nzm - 1))) = zero;
    }
    for (int index_0 = 0; index_0 < thlp2_mc.extent_int(0); ++index_0) {
      thlp2_mc(index_0, ((nzm - 1))) = zero;
    }
    for (int index_0 = 0; index_0 < rtpthlp_mc.extent_int(0); ++index_0) {
      rtpthlp_mc(index_0, ((nzm - 1))) = zero;
    }

  } else {

    // Microphysics tendency terms for model variances and covariances
    // are set to 0.
    for (int index_1 = 0; index_1 < wprtp_mc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < wprtp_mc.extent_int(0); ++index_0) {
        wprtp_mc(index_0, index_1) = zero;
      }
    }
    for (int index_1 = 0; index_1 < wpthlp_mc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < wpthlp_mc.extent_int(0); ++index_0) {
        wpthlp_mc(index_0, index_1) = zero;
      }
    }
    for (int index_1 = 0; index_1 < rtp2_mc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < rtp2_mc.extent_int(0); ++index_0) {
        rtp2_mc(index_0, index_1) = zero;
      }
    }
    for (int index_1 = 0; index_1 < thlp2_mc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < thlp2_mc.extent_int(0); ++index_0) {
        thlp2_mc(index_0, index_1) = zero;
      }
    }
    for (int index_1 = 0; index_1 < rtpthlp_mc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < rtpthlp_mc.extent_int(0); ++index_0) {
        rtpthlp_mc(index_0, index_1) = zero;
      }
    }
  }

  //!! Boundary conditions for microphysics tendencies.
  for (int index_0 = 0; index_0 < rrm_mc.extent_int(0); ++index_0) {
    rrm_mc(index_0, (nzt - 1)) = zero;
  }
  for (int index_0 = 0; index_0 < Nrm_mc.extent_int(0); ++index_0) {
    Nrm_mc(index_0, (nzt - 1)) = zero;
  }

  // Boundary conditions
  for (int index_0 = 0; index_0 < KK_mean_vol_rad.extent_int(0); ++index_0) {
    KK_mean_vol_rad(index_0, (nzt - 1)) = zero;
  }
  for (int index_0 = 0; index_0 < rvm_mc.extent_int(0); ++index_0) {
    rvm_mc(index_0, (nzt - 1)) = zero;
  }
  for (int index_0 = 0; index_0 < rcm_mc.extent_int(0); ++index_0) {
    rcm_mc(index_0, (nzt - 1)) = zero;
  }
  for (int index_0 = 0; index_0 < thlm_mc.extent_int(0); ++index_0) {
    thlm_mc(index_0, (nzt - 1)) = zero;
  }

  // Find the vertical level index of cloud top.
  Kokkos::deep_copy(cloud_top_level,
                    to_host_space(workspace, clubb_kokkos::get_cloud_top_level(
                                                 nzt, ngrdcol, to_execution_space(workspace, rcm),
                                                 to_execution_space(workspace, hydromet),
                                                 hydromet_dim, hm_metadata_iiri)));

  //!! Microphysics sedimentation velocities.
  microphys_host::KK_upscaled_sedimentation(ngrdcol, nzt, cloud_top_level, KK_mean_vol_rad,
                                            l_clip_positive_sed, Vrr, VNr);

  //!! Output hydrometeor mean tendencies and mean sedimentation velocities
  //!! in output arrays.
  for (int index_2 = 0; index_2 < hydromet_mc.extent_int(2); ++index_2) {
    for (int index_1 = 0; index_1 < hydromet_mc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_mc.extent_int(0); ++index_0) {
        hydromet_mc(index_0, index_1, index_2) = zero;
      }
    }
  }
  for (int index_2 = 0; index_2 < hydromet_vel.extent_int(2); ++index_2) {
    for (int index_1 = 0; index_1 < hydromet_vel.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_vel.extent_int(0); ++index_0) {
        hydromet_vel(index_0, index_1, index_2) = zero;
      }
    }
  }
  for (int index_1 = 0; index_1 < hydromet_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < hydromet_mc.extent_int(0); ++index_0) {
      hydromet_mc(index_0, index_1, hm_metadata_iirr) = rrm_mc(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < hydromet_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < hydromet_mc.extent_int(0); ++index_0) {
      hydromet_mc(index_0, index_1, hm_metadata_iiNr) = Nrm_mc(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < hydromet_vel.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < hydromet_vel.extent_int(0); ++index_0) {
      hydromet_vel(index_0, index_1, hm_metadata_iirr) = Vrr(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < hydromet_vel.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < hydromet_vel.extent_int(0); ++index_0) {
      hydromet_vel(index_0, index_1, hm_metadata_iiNr) = VNr(index_0, index_1);
    }
  }

  //!! Turbulent sedimentation above cloud top should have a value of 0.
  for (k = 0; k < (nzt - 1); ++k) {
    for (i = 0; i < ngrdcol; ++i) {
      if (((cloud_top_level(i) > 0) && (k > cloud_top_level(i)))) {
        Vrrprrp_zt_impc(i, k) = zero;
        Vrrprrp_zt_expc(i, k) = zero;
        VNrpNrp_zt_impc(i, k) = zero;
        VNrpNrp_zt_expc(i, k) = zero;
      }
    }
  }

  //!! Boundary conditions (upper) for the covariances of hydrometeor
  //!! sedimentation velocities and their associated hydrometeors
  //!! (<V_rr'r_r'> and <V_Nr'N_r'>).
  for (int index_0 = 0; index_0 < Vrrprrp_zt_impc.extent_int(0); ++index_0) {
    Vrrprrp_zt_impc(index_0, (nzt - 1)) = zero;
  }
  for (int index_0 = 0; index_0 < Vrrprrp_zt_expc.extent_int(0); ++index_0) {
    Vrrprrp_zt_expc(index_0, (nzt - 1)) = zero;
  }
  for (int index_0 = 0; index_0 < VNrpNrp_zt_impc.extent_int(0); ++index_0) {
    VNrpNrp_zt_impc(index_0, (nzt - 1)) = zero;
  }
  for (int index_0 = 0; index_0 < VNrpNrp_zt_expc.extent_int(0); ++index_0) {
    VNrpNrp_zt_expc(index_0, (nzt - 1)) = zero;
  }

  // The implicit and explicit components used to calculate the covariances of
  // hydrometeor sedimentation velocities and their associated hydrometeors
  // (<V_rr'r_r'> and <V_Nr'N_r'>) are fed into the output arrays.
  for (int index_2 = 0; index_2 < hydromet_vel_covar_zt_impc.extent_int(2); ++index_2) {
    for (int index_1 = 0; index_1 < hydromet_vel_covar_zt_impc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_vel_covar_zt_impc.extent_int(0); ++index_0) {
        hydromet_vel_covar_zt_impc(index_0, index_1, index_2) = zero;
      }
    }
  }
  for (int index_2 = 0; index_2 < hydromet_vel_covar_zt_expc.extent_int(2); ++index_2) {
    for (int index_1 = 0; index_1 < hydromet_vel_covar_zt_expc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_vel_covar_zt_expc.extent_int(0); ++index_0) {
        hydromet_vel_covar_zt_expc(index_0, index_1, index_2) = zero;
      }
    }
  }
  for (int index_1 = 0; index_1 < hydromet_vel_covar_zt_impc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < hydromet_vel_covar_zt_impc.extent_int(0); ++index_0) {
      hydromet_vel_covar_zt_impc(index_0, index_1, hm_metadata_iirr) =
          Vrrprrp_zt_impc(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < hydromet_vel_covar_zt_expc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < hydromet_vel_covar_zt_expc.extent_int(0); ++index_0) {
      hydromet_vel_covar_zt_expc(index_0, index_1, hm_metadata_iirr) =
          Vrrprrp_zt_expc(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < hydromet_vel_covar_zt_impc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < hydromet_vel_covar_zt_impc.extent_int(0); ++index_0) {
      hydromet_vel_covar_zt_impc(index_0, index_1, hm_metadata_iiNr) =
          VNrpNrp_zt_impc(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < hydromet_vel_covar_zt_expc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < hydromet_vel_covar_zt_expc.extent_int(0); ++index_0) {
      hydromet_vel_covar_zt_expc(index_0, index_1, hm_metadata_iiNr) =
          VNrpNrp_zt_expc(index_0, index_1);
    }
  }

  if (stats.l_sample) {
    //!! Statistical output for upscaled KK.
    // The covariance < V_rr'r_r' > calculated completely explicitly.
    // When semi-implicit turbulent advection is used, this result can be
    // compared to the < V_rr'r_r' > results used in the code, which are
    // calculated semi-implicitly.

    HostReal2D array_actual_2 = workspace_scope.view<HostReal2D>(
        "kk_upscaled_microphys::array_actual_2", Vrrprrp_zt_impc.extent_int(0),
        Vrrprrp_zt_impc.extent_int(1));
    for (int index_1 = 0; index_1 < Vrrprrp_zt_impc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < Vrrprrp_zt_impc.extent_int(0); ++index_0) {
        array_actual_2(index_0, index_1) =
            ((Vrrprrp_zt_impc(index_0, index_1) * rrm(index_0, index_1)) +
             Vrrprrp_zt_expc(index_0, index_1));
      }
    }
    const auto array_actual_1 = zt2zm_api(nzm, nzt, ngrdcol, gr, array_actual_2);
    stats_update(std::string("Vrrprrp_expcalc"), array_actual_1, stats);
    // The covariance < V_Nr'N_r' > calculated completely explicitly.
    // When semi-implicit turbulent advection is used, this result can be
    // compared to the < V_Nr'N_r' > results used in the code, which are
    // calculated semi-implicitly.

    HostReal2D array_actual_4 = workspace_scope.view<HostReal2D>(
        "kk_upscaled_microphys::array_actual_4", VNrpNrp_zt_impc.extent_int(0),
        VNrpNrp_zt_impc.extent_int(1));
    for (int index_1 = 0; index_1 < VNrpNrp_zt_impc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < VNrpNrp_zt_impc.extent_int(0); ++index_0) {
        array_actual_4(index_0, index_1) =
            ((VNrpNrp_zt_impc(index_0, index_1) * Nrm(index_0, index_1)) +
             VNrpNrp_zt_expc(index_0, index_1));
      }
    }
    const auto array_actual_3 = zt2zm_api(nzm, nzt, ngrdcol, gr, array_actual_4);
    stats_update(std::string("VNrpNrp_expcalc"), array_actual_3, stats);
    //!! Statistical output for mean microphysics tendencies.
    // Stats output for microphysics adjustment terms.
    HostReal2D array_actual_5 = workspace_scope.view<HostReal2D>(
        "kk_upscaled_microphys::array_actual_5", adj_terms.extent_int(0), adj_terms.extent_int(1));
    for (int index_1 = 0; index_1 < adj_terms.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < adj_terms.extent_int(0); ++index_0) {
        array_actual_5(index_0, index_1) = adj_terms(index_0, index_1).rrm_src_adj;
      }
    }
    stats_update(std::string("rrm_src_adj"), array_actual_5, stats);
    HostReal2D array_actual_6 = workspace_scope.view<HostReal2D>(
        "kk_upscaled_microphys::array_actual_6", adj_terms.extent_int(0), adj_terms.extent_int(1));
    for (int index_1 = 0; index_1 < adj_terms.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < adj_terms.extent_int(0); ++index_0) {
        array_actual_6(index_0, index_1) = adj_terms(index_0, index_1).Nrm_src_adj;
      }
    }
    stats_update(std::string("Nrm_src_adj"), array_actual_6, stats);
    HostReal2D array_actual_7 = workspace_scope.view<HostReal2D>(
        "kk_upscaled_microphys::array_actual_7", adj_terms.extent_int(0), adj_terms.extent_int(1));
    for (int index_1 = 0; index_1 < adj_terms.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < adj_terms.extent_int(0); ++index_0) {
        array_actual_7(index_0, index_1) = adj_terms(index_0, index_1).rrm_evap_adj;
      }
    }
    stats_update(std::string("rrm_evap_adj"), array_actual_7, stats);
    HostReal2D array_actual_8 = workspace_scope.view<HostReal2D>(
        "kk_upscaled_microphys::array_actual_8", adj_terms.extent_int(0), adj_terms.extent_int(1));
    for (int index_1 = 0; index_1 < adj_terms.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < adj_terms.extent_int(0); ++index_0) {
        array_actual_8(index_0, index_1) = adj_terms(index_0, index_1).Nrm_evap_adj;
      }
    }
    stats_update(std::string("Nrm_evap_adj"), array_actual_8, stats);
    HostReal2D array_actual_9 =
        workspace_scope.view<HostReal2D>("kk_upscaled_microphys::array_actual_9",
                                         KK_auto_tndcy.extent_int(0), KK_auto_tndcy.extent_int(1));
    for (int index_1 = 0; index_1 < KK_auto_tndcy.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < KK_auto_tndcy.extent_int(0); ++index_0) {
        array_actual_9(index_0, index_1) =
            ((KK_auto_tndcy(index_0, index_1) + KK_accr_tndcy(index_0, index_1)) +
             KK_evap_tndcy(index_0, index_1));
      }
    }
    stats_update(std::string("rrm_mc_nonadj"), array_actual_9, stats);
  }

  return;
}
inline void KK_tendency_coefs(real_t thlm, real_t exner, real_t p_in_Pa, real_t rho,
                              real_t& KK_evap_coef, real_t& KK_auto_coef, real_t& KK_accr_coef,
                              real_t& KK_mvr_coef, int saturation_formula,
                              const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (3), Eq. (22), Eq. (29), and Eq. (33) of Khairoutdinov, M. and
  // Y. Kogan, 2000:  A New Cloud Physics Parameterization in a Large-Eddy
  // Simulation Model of Marine Stratocumulus.  Mon. Wea. Rev., 128, 229--243.
  //
  // Eq. (22), Eq. (28), Eq. (38), and Eq. (51) of Larson, V. E. and
  // B. M. Griffin, 2013:  Analytic upscaling of a local microphysics scheme.
  // Part I: Derivation.  Q. J. Roy. Meteorol. Soc., 139, 670, 46--57,
  // doi:http://dx.doi.org/10.1002/qj.1967.
  //
  // Eq. (C21) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S21) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // A new subgrid-scale representation of hydrometeor fields using a
  // multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //
  // Eq. (A27) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
  // microphysical effects on variances and covariances of moisture and heat
  // content using a multivariate probability density function: a study with
  // CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
  // Input Variables
  // Output Variables
  // Local Variables

  // In: thlm — Mean liquid water potential temperature         [K]
  // In: exner — Exner function                                  [-]
  // In: p_in_Pa — Pressure                                        [Pa]
  // In: rho — Density                                         [kg/m^3]
  // Out: KK_evap_coef — KK evaporation coefficient                 [(kg/kg)/s]
  // Out: KK_auto_coef — KK autoconversion coefficient              [(kg/kg)/s]
  // Out: KK_accr_coef — KK accretion coefficient                   [(kg/kg)/s]
  // Out: KK_mvr_coef — KK mean volume radius coefficient          [m]
  // In: saturation_formula — Integer that stores the saturation formula to be used
  // Case state: kk replaces the Native parameters_KK module.

  real_t T_liq_in_K; // Mean liquid water temperature, T_l           [K]
  real_t r_sl;       // Liquid water sat. mixing ratio, r_s(T_l,p)   [kg/kg]
  real_t Beta_Tl;    // Parameter Beta, Beta(T_l)                    [1/(kg/kg)]

  // Liquid water temperature.
  T_liq_in_K = (thlm * exner);

  // Saturation mixing ratio (based on liquid water temperature and
  // pressure), r_sl = r_s(T_l,p).
  r_sl = sat_mixrat_liq_api(p_in_Pa, T_liq_in_K, saturation_formula);

  // Beta(T_l).
  Beta_Tl = (((Rd / Rv) * (Lv / (Rd * T_liq_in_K))) * (Lv / (Cp * T_liq_in_K)));

  // Coefficient for KK evaporation.
  KK_evap_coef =
      ((((three * kk.C_evap) * kk_engine::G_T_p(T_liq_in_K, p_in_Pa, saturation_formula, kk)) *
        Kokkos::pow(((four_thirds * pi) * rho_lw), two_thirds)) *
       ((one + (Beta_Tl * r_sl)) / r_sl));

  // Coefficient for KK autoconversion.
  KK_auto_coef = (real_t(1350.0) * Kokkos::pow((rho / cm3_per_m3), KK_auto_Nc_exp));

  // Coefficient for KK accretion.
  KK_accr_coef = real_t(67.0);

  // Coefficient for KK rain drop mean volume radius.
  KK_mvr_coef = Kokkos::pow(((four_thirds * pi) * rho_lw), (-one_third));

  return;
}
inline void KK_microphys_adjust(real_t dt, real_t exner, real_t rcm, real_t rrm, real_t Nrm,
                                real_t KK_evap_tndcy, real_t KK_auto_tndcy, real_t KK_accr_tndcy,
                                real_t KK_Nrm_evap_tndcy, real_t KK_Nrm_auto_tndcy,
                                bool l_src_adj_enabled, bool l_evap_adj_enabled, real_t& rrm_mc,
                                real_t& Nrm_mc, real_t& rvm_mc, real_t& rcm_mc, real_t& thlm_mc,
                                KK_microphys_adj_terms_type& adj_terms, const parameters_KK& kk) {
  // Description:
  // References:
  // Input Variables
  // Output Variables
  // Local Variables

  // In: dt — Model time step duration                 [s]
  // In: exner — Exner function                           [-]
  // In: rcm — Mean cloud water mixing ratio            [kg/kg]
  // In: rrm — Mean rain water mixing ratio, < r_r >    [kg/kg]
  // In: Nrm — Mean rain drop concentration, < N_r >    [num/kg]
  // In: KK_evap_tndcy — Mean KK <r_r> evaporation tendency     [(kg/kg)/s]
  // In: KK_auto_tndcy — Mean KK <r_r> autoconversion tendency  [(kg/kg)/s]
  // In: KK_accr_tndcy — Mean KK <r_r> accretion tendency       [(kg/kg)/s]
  // In: KK_Nrm_evap_tndcy — Mean KK <N_r> evaporation tendency     [(num/kg)/s]
  // In: KK_Nrm_auto_tndcy — Mean KK <N_r> autoconversion tendency  [(num/kg)/s]
  // In: l_src_adj_enabled — Flag to enable rrm/Nrm source adjustment
  // In: l_evap_adj_enabled — Flag to enable rrm/Nrm evaporation adjustment
  // Out: rrm_mc — Mean (dr_r/dt) due to microphysics       [(kg/kg)/s]
  // Out: Nrm_mc — Mean (dN_r/dt) due to microphysics       [(num/kg)/s]
  // Out: rvm_mc — Time tendency of vapor water mixing ratio        [kg/kg/s]
  // Out: rcm_mc — Time tendency of liquid water mixing ratio       [kg/kg/s]
  // Out: thlm_mc — Time tendency of liquid potential temperature    [K/s]
  // Out: adj_terms — Structure containing information about adjustment of microphysics terms
  // Case state: kk replaces the Native parameters_KK module.

  real_t rrm_source;     // Total source term rate for rrm           [(kg/kg)/s]
  real_t Nrm_source;     // Total source term rate for Nrm           [(num/kg)/s]
  real_t rrm_src_adj;    // Total adjustment to rrm source terms     [(kg/kg)/s]
  real_t Nrm_src_adj;    // Total adjustment to Nrm source terms     [(num/kg)/s]
  real_t rrm_evap_net;   // Net evaporation rate of <r_r>            [(kg/kg)/s]
  real_t Nrm_evap_net;   // Net evaporation rate of <N_r>            [(num/kg)/s]
  real_t rrm_src_max;    // Maximum allowable rrm source rate    [(kg/kg)/s]
  real_t rrm_auto_ratio; // Ratio of rrm autoconv to overall source term [-]
  real_t
      total_rc_needed; // Amount of r_c needed to over the timestep for rain source terms                    [kg/kg] -----------------------------------------------------------------------
  //---- Begin Code -----

  //!! Source-adjustment code for rrm and Nrm.
  rrm_source = (KK_auto_tndcy + KK_accr_tndcy);
  Nrm_source = KK_Nrm_auto_tndcy;

  if (l_src_adj_enabled) {

    // The increase of rain due to autoconversion and accretion both draw
    // water from the available cloud water.  Over a long time step, these
    // rates may over-deplete cloud water.  In other words, these processes
    // may draw more cloud water than there is available.  Thus, the total
    // source rate multiplied by the duration of the time step cannot exceed
    // the total amount of cloud water available.  If it does, then the rate
    // must be adjusted.
    total_rc_needed = (rrm_source * dt);

    if ((total_rc_needed > rcm)) {

      // The maximum allowable rate of the source terms is rcm/dt.
      rrm_src_max = (rcm / dt);

      // The amount of adjustment to the source terms.
      // This value should always be negative.
      rrm_src_adj = (rrm_src_max - rrm_source);

      // Reset the value of the source terms to the maximum allowable value
      // of the source terms.
      rrm_source = rrm_src_max;

      // The rrm source terms are made up of autoconversion and accretion.
      // Only the sum of those two terms is corrected.  However, Nrm has only
      // an autoconversion term for a source term.  Assume that change in the
      // rrm autoconversion term is proportional to the total rrm adjustment
      // rate by the ratio of rrm autoconversion to the overall source term.
      // Then, plug the rrm autoconversion adjustment into the equation for
      // Nrm autoconversion to determine the effect on the Nrm source term.
      rrm_auto_ratio = (KK_auto_tndcy / (KK_auto_tndcy + KK_accr_tndcy));
      Nrm_src_adj    = kk_engine::KK_Nrm_auto_mean((rrm_auto_ratio * rrm_src_adj), kk);

      // Change Nrm by Nrm_src_adj.  Nrm_src_adj will always be negative.
      Nrm_source = (Nrm_source + Nrm_src_adj);

    } else {
      // total_rc_needed <= rcm: enough available cloud water.

      rrm_src_adj = zero;
      Nrm_src_adj = zero;
    }

  } else {
    // l_src_adj_enabled is false: source adjustment is disabled.

    rrm_src_adj = zero;
    Nrm_src_adj = zero;
  }

  //!! Evaporation-adjustment code for rrm and Nrm.
  if (l_evap_adj_enabled) {

    // Prevent over-evaporation of rain over a long model time step.
    // Limit the evaporation rate of both rain water mixing ratio and rain
    // drop concentration.  The total amount of rain lost due to evaporation
    // cannot be so great as to result in negative rain.

    // Calculate net evaporation rate of <r_r>.
    rrm_evap_net = ((KK_evap_tndcy) > (((-rrm) / dt)) ? (KK_evap_tndcy) : (((-rrm) / dt)));

    // Recalcuate the net evaporation rate of <N_r> based on the net
    // evaporation rate of <r_r>.
    if ((((Kokkos::abs((KK_evap_tndcy - rrm_evap_net)) >
           ((Kokkos::abs((KK_evap_tndcy + rrm_evap_net)) * eps) / 2)) &&
          (rrm > rr_tol)) &&
         (Nrm > Nr_tol))) {
      Nrm_evap_net = kk_engine::KK_Nrm_evap_local_mean(rrm_evap_net, Nrm, rrm, dt, kk);
    } else {
      Nrm_evap_net = KK_Nrm_evap_tndcy;
    }

    Nrm_evap_net = ((Nrm_evap_net) > (((-Nrm) / dt)) ? (Nrm_evap_net) : (((-Nrm) / dt)));

  } else {
    // l_evap_adj_enabled is false: evaporation adjustment is disabled.

    rrm_evap_net = KK_evap_tndcy;
    Nrm_evap_net = KK_Nrm_evap_tndcy;
  }

  //!! Calculate overall KK microphysics tendencies.
  rrm_mc = (rrm_evap_net + rrm_source);
  Nrm_mc = (Nrm_evap_net + Nrm_source);

  //!! Explicit contributions to thlm and rtm from the microphysics
  rvm_mc = (-rrm_evap_net);
  // Accretion + Autoconversion
  rcm_mc  = (-rrm_source);
  thlm_mc = ((Lv / (Cp * exner)) * rrm_mc);

  // Return adjustment terms
  adj_terms.rrm_src_adj  = rrm_src_adj;
  adj_terms.Nrm_src_adj  = Nrm_src_adj;
  adj_terms.rrm_evap_adj = (rrm_evap_net - KK_evap_tndcy);
  adj_terms.Nrm_evap_adj = (Nrm_evap_net - KK_Nrm_evap_tndcy);

  return;
}
template <class Array0, class Array1, class Array2, class Array3>
inline void KK_upscaled_sedimentation(int ngrdcol, int nzt, Array0 cloud_top_level,
                                      Array1 KK_mean_vol_rad, bool l_clip_positive_sed, Array2 Vrr,
                                      Array3 VNr) {
  // Description:
  // Calculate sedimentation velocities for every upscaled K-K model column.
  // References:
  // Eq. (37) of Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
  // Parameterization in a Large-Eddy Simulation Model of Marine Stratocumulus.
  // Mon. Wea. Rev., 128, 229--243.
  // Input Variables
  // Output Variables
  // Local Variables

  // In: ngrdcol — Number of model columns
  // In: nzt — Number of model thermodynamic vertical grid levels
  // In: cloud_top_level (ngrdcol) — Vertical level index of cloud top
  // In: KK_mean_vol_rad (ngrdcol,nzt) — KK rain drop mean volume radius       [m]
  // In: l_clip_positive_sed — Clip positive values of Vrr and VNr   [T/F]
  // Out: Vrr (ngrdcol,nzt) — Mean sedimentation velocity of < r_r >   [m/s]
  // Out: VNr (ngrdcol,nzt) — Mean sedimentation velocity of < N_r >   [m/s]

  int i; // Column loop iterator
  int k;

  //!! Mean sedimentation velocities
  for (k = 0; k < (nzt - 1); ++k) {
    for (i = 0; i < ngrdcol; ++i) {

      // Mean sedimentation velocity of rain water mixing ratio.
      Vrr(i, k) = (-((real_t(0.012) * (micron_per_m * KK_mean_vol_rad(i, k))) - real_t(0.2)));

      // Mean sedimentation velocity of rain drop concentration.
      VNr(i, k) = (-((real_t(0.007) * (micron_per_m * KK_mean_vol_rad(i, k))) - real_t(0.1)));

      if (l_clip_positive_sed) {

        // Mean sedimentation velocity of rain water mixing ratio cannot
        // have a
        // positive value.
        if ((Vrr(i, k) > zero)) {
          Vrr(i, k) = zero;
        }

        // Mean sedimentation velocity of rain drop concentration cannot have a
        // positive value.
        if ((VNr(i, k) > zero)) {
          VNr(i, k) = zero;
        }

        if (((cloud_top_level(i) > 0) && (k > cloud_top_level(i)))) {
          Vrr(i, k) = zero;
          VNr(i, k) = zero;
        }
      }
    }
  }

  //!! Boundary conditions for sedimentation velocities.

  // The flux of rain water through the model top is 0.
  // Vrr and VNr are set to 0 at the highest model level.
  for (int index_0 = 0; index_0 < Vrr.extent_int(0); ++index_0) {
    Vrr(index_0, (nzt - 1)) = zero;
  }
  for (int index_0 = 0; index_0 < VNr.extent_int(0); ++index_0) {
    VNr(index_0, (nzt - 1)) = zero;
  }

  return;
}
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5>
inline void KK_microphys_output(int ngrdcol, int nzt, int hydromet_dim,
                                const hm_metadata_type& hm_metadata, Array0 Vrr, Array1 VNr,
                                Array2 rrm_mc, Array3 Nrm_mc, Array4 hydromet_mc,
                                Array5 hydromet_vel) {
  // Description:
  // References:
  // Input Variables
  // Output Variables
  // In: ngrdcol — Number of model columns
  // In: nzt — Number of model thermodynamic vertical grid levels
  // In: hydromet_dim
  // In: hm_metadata
  // In: Vrr (ngrdcol,nzt) — Mean sedimentation velocity of < r_r >   [m/s]
  // In: VNr (ngrdcol,nzt) — Mean sedimentation velocity of < N_r >   [m/s]
  // In: rrm_mc (ngrdcol,nzt) — Mean (dr_r/dt) due to microphysics       [(kg/kg)/s]
  // In: Nrm_mc (ngrdcol,nzt) — Mean (dN_r/dt) due to microphysics       [(num/kg)/s]
  // Out: hydromet_mc (ngrdcol,nzt,hydromet_dim) — Hydrometeor time tendency          [(units vary)/s]
  // Out: hydromet_vel (ngrdcol,nzt,hydromet_dim) — Hydrometeor sedimentation velocity [m/s]

  const auto hm_metadata_iirr = hm_metadata.iirr;
  const auto hm_metadata_iiNr = hm_metadata.iiNr;

  //!! Output mean hydrometeor tendencies and mean sedimentation velocities.

  for (int index_2 = 0; index_2 < hydromet_mc.extent_int(2); ++index_2) {
    for (int index_1 = 0; index_1 < hydromet_mc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_mc.extent_int(0); ++index_0) {
        hydromet_mc(index_0, index_1, index_2) = zero;
      }
    }
  }
  for (int index_2 = 0; index_2 < hydromet_vel.extent_int(2); ++index_2) {
    for (int index_1 = 0; index_1 < hydromet_vel.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_vel.extent_int(0); ++index_0) {
        hydromet_vel(index_0, index_1, index_2) = zero;
      }
    }
  }

  // Mean field tendencies.
  for (int index_1 = 0; index_1 < hydromet_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < hydromet_mc.extent_int(0); ++index_0) {
      hydromet_mc(index_0, index_1, hm_metadata_iirr) = rrm_mc(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < hydromet_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < hydromet_mc.extent_int(0); ++index_0) {
      hydromet_mc(index_0, index_1, hm_metadata_iiNr) = Nrm_mc(index_0, index_1);
    }
  }

  // Sedimentation Velocities.
  for (int index_1 = 0; index_1 < hydromet_vel.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < hydromet_vel.extent_int(0); ++index_0) {
      hydromet_vel(index_0, index_1, hm_metadata_iirr) = Vrr(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < hydromet_vel.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < hydromet_vel.extent_int(0); ++index_0) {
      hydromet_vel(index_0, index_1, hm_metadata_iiNr) = VNr(index_0, index_1);
    }
  }

  return;
}
} // namespace clubb_kokkos::microphys_host

namespace clubb_kokkos {
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24>
inline void KK_local_microphys_driver(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int nzt, int hydromet_dim,
    const hm_metadata_type& hm_metadata, bool l_latin_hypercube, Array0 thlm, Array1 wm_zt,
    Array2 p_in_Pa, Array3 exner, Array4 rho, Array5 cloud_frac, Array6 w_std_dev, Array7 dzq,
    Array8 rcm, Array9 Ncm, Array10 chi, Array11 rvm, Array12 hydromet, int saturation_formula,
    Array13 sample_weight, stats_type& stats, Array14 hydromet_mc, Array15 hydromet_vel,
    Array16 Ncm_mc, Array17 rcm_mc, Array18 rvm_mc, Array19 thlm_mc, Array20 rrm_auto_diag,
    Array21 rrm_accr_diag, Array22 rrm_evap_diag, Array23 Nrm_auto_diag, Array24 Nrm_evap_diag,
    const parameters_microphys& microphys, const parameters_KK& kk) {
  // Stage the Native engine call boundary into call-local host Views.
  // The C++ science implementation above retains Native float casts/order;
  // model arrays and every InOut/Out allocation remain caller owned.
  auto workspace_scope    = workspace.scope();
  auto thlm_host          = to_host_space(workspace, thlm);
  auto wm_zt_host         = to_host_space(workspace, wm_zt);
  auto p_in_Pa_host       = to_host_space(workspace, p_in_Pa);
  auto exner_host         = to_host_space(workspace, exner);
  auto rho_host           = to_host_space(workspace, rho);
  auto cloud_frac_host    = to_host_space(workspace, cloud_frac);
  auto w_std_dev_host     = to_host_space(workspace, w_std_dev);
  auto dzq_host           = to_host_space(workspace, dzq);
  auto rcm_host           = to_host_space(workspace, rcm);
  auto Ncm_host           = to_host_space(workspace, Ncm);
  auto chi_host           = to_host_space(workspace, chi);
  auto rvm_host           = to_host_space(workspace, rvm);
  auto hydromet_host      = to_host_space(workspace, hydromet);
  auto sample_weight_host = to_host_space(workspace, sample_weight);
  auto hydromet_mc_host   = to_host_space(workspace, hydromet_mc);
  auto hydromet_vel_host  = to_host_space(workspace, hydromet_vel);
  auto Ncm_mc_host        = to_host_space(workspace, Ncm_mc);
  auto rcm_mc_host        = to_host_space(workspace, rcm_mc);
  auto rvm_mc_host        = to_host_space(workspace, rvm_mc);
  auto thlm_mc_host       = to_host_space(workspace, thlm_mc);
  auto rrm_auto_diag_host = to_host_space(workspace, rrm_auto_diag);
  auto rrm_accr_diag_host = to_host_space(workspace, rrm_accr_diag);
  auto rrm_evap_diag_host = to_host_space(workspace, rrm_evap_diag);
  auto Nrm_auto_diag_host = to_host_space(workspace, Nrm_auto_diag);
  auto Nrm_evap_diag_host = to_host_space(workspace, Nrm_evap_diag);
  microphys_host::KK_local_microphys_driver(
      workspace, gr, ngrdcol, dt, nzt, hydromet_dim, hm_metadata, l_latin_hypercube, thlm_host,
      wm_zt_host, p_in_Pa_host, exner_host, rho_host, cloud_frac_host, w_std_dev_host, dzq_host,
      rcm_host, Ncm_host, chi_host, rvm_host, hydromet_host, saturation_formula, sample_weight_host,
      stats, hydromet_mc_host, hydromet_vel_host, Ncm_mc_host, rcm_mc_host, rvm_mc_host,
      thlm_mc_host, rrm_auto_diag_host, rrm_accr_diag_host, rrm_evap_diag_host, Nrm_auto_diag_host,
      Nrm_evap_diag_host, microphys, kk);
  from_execution_space(workspace, hydromet_mc, hydromet_mc_host);
  from_execution_space(workspace, hydromet_vel, hydromet_vel_host);
  from_execution_space(workspace, Ncm_mc, Ncm_mc_host);
  from_execution_space(workspace, rcm_mc, rcm_mc_host);
  from_execution_space(workspace, rvm_mc, rvm_mc_host);
  from_execution_space(workspace, thlm_mc, thlm_mc_host);
  from_execution_space(workspace, rrm_auto_diag, rrm_auto_diag_host);
  from_execution_space(workspace, rrm_accr_diag, rrm_accr_diag_host);
  from_execution_space(workspace, rrm_evap_diag, rrm_evap_diag_host);
  from_execution_space(workspace, Nrm_auto_diag, Nrm_auto_diag_host);
  from_execution_space(workspace, Nrm_evap_diag, Nrm_evap_diag_host);
}
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25>
inline void KK_upscaled_microphys(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int nzt, int nzm,
    int pdf_dim, int hydromet_dim, const hm_metadata_type& hm_metadata, Array0 wm_zt, Array1 rtm,
    Array2 thlm, Array3 p_in_Pa, Array4 exner, Array5 rho, Array6 rcm,
    const pdf_parameter& pdf_params, const hydromet_pdf_parameter_collection& hydromet_pdf_params,
    const precipitation_fractions& precip_fracs, Array7 hydromet, Array8 mu_x_1_n, Array9 mu_x_2_n,
    Array10 sigma_x_1_n, Array11 sigma_x_2_n, Array12 corr_array_1_n, Array13 corr_array_2_n,
    int saturation_formula, stats_type& stats, Array14 hydromet_mc, Array15 hydromet_vel,
    Array16 rcm_mc, Array17 rvm_mc, Array18 thlm_mc, Array19 hydromet_vel_covar_zt_impc,
    Array20 hydromet_vel_covar_zt_expc, Array21 wprtp_mc, Array22 wpthlp_mc, Array23 rtp2_mc,
    Array24 thlp2_mc, Array25 rtpthlp_mc, const parameters_microphys& microphys,
    const parameters_KK& kk) {
  // Stage the Native engine call boundary into call-local host Views.
  // The C++ science implementation above retains Native float casts/order;
  // model arrays and every InOut/Out allocation remain caller owned.
  auto workspace_scope                 = workspace.scope();
  auto wm_zt_host                      = to_host_space(workspace, wm_zt);
  auto rtm_host                        = to_host_space(workspace, rtm);
  auto thlm_host                       = to_host_space(workspace, thlm);
  auto p_in_Pa_host                    = to_host_space(workspace, p_in_Pa);
  auto exner_host                      = to_host_space(workspace, exner);
  auto rho_host                        = to_host_space(workspace, rho);
  auto rcm_host                        = to_host_space(workspace, rcm);
  auto hydromet_host                   = to_host_space(workspace, hydromet);
  auto mu_x_1_n_host                   = to_host_space(workspace, mu_x_1_n);
  auto mu_x_2_n_host                   = to_host_space(workspace, mu_x_2_n);
  auto sigma_x_1_n_host                = to_host_space(workspace, sigma_x_1_n);
  auto sigma_x_2_n_host                = to_host_space(workspace, sigma_x_2_n);
  auto corr_array_1_n_host             = to_host_space(workspace, corr_array_1_n);
  auto corr_array_2_n_host             = to_host_space(workspace, corr_array_2_n);
  auto hydromet_mc_host                = to_host_space(workspace, hydromet_mc);
  auto hydromet_vel_host               = to_host_space(workspace, hydromet_vel);
  auto rcm_mc_host                     = to_host_space(workspace, rcm_mc);
  auto rvm_mc_host                     = to_host_space(workspace, rvm_mc);
  auto thlm_mc_host                    = to_host_space(workspace, thlm_mc);
  auto hydromet_vel_covar_zt_impc_host = to_host_space(workspace, hydromet_vel_covar_zt_impc);
  auto hydromet_vel_covar_zt_expc_host = to_host_space(workspace, hydromet_vel_covar_zt_expc);
  auto wprtp_mc_host                   = to_host_space(workspace, wprtp_mc);
  auto wpthlp_mc_host                  = to_host_space(workspace, wpthlp_mc);
  auto rtp2_mc_host                    = to_host_space(workspace, rtp2_mc);
  auto thlp2_mc_host                   = to_host_space(workspace, thlp2_mc);
  auto rtpthlp_mc_host                 = to_host_space(workspace, rtpthlp_mc);
  microphys_host::KK_upscaled_microphys(
      workspace, gr, ngrdcol, dt, nzt, nzm, pdf_dim, hydromet_dim, hm_metadata, wm_zt_host,
      rtm_host, thlm_host, p_in_Pa_host, exner_host, rho_host, rcm_host, pdf_params,
      hydromet_pdf_params, precip_fracs, hydromet_host, mu_x_1_n_host, mu_x_2_n_host,
      sigma_x_1_n_host, sigma_x_2_n_host, corr_array_1_n_host, corr_array_2_n_host,
      saturation_formula, stats, hydromet_mc_host, hydromet_vel_host, rcm_mc_host, rvm_mc_host,
      thlm_mc_host, hydromet_vel_covar_zt_impc_host, hydromet_vel_covar_zt_expc_host, wprtp_mc_host,
      wpthlp_mc_host, rtp2_mc_host, thlp2_mc_host, rtpthlp_mc_host, microphys, kk);
  from_execution_space(workspace, hydromet_mc, hydromet_mc_host);
  from_execution_space(workspace, hydromet_vel, hydromet_vel_host);
  from_execution_space(workspace, rcm_mc, rcm_mc_host);
  from_execution_space(workspace, rvm_mc, rvm_mc_host);
  from_execution_space(workspace, thlm_mc, thlm_mc_host);
  from_execution_space(workspace, hydromet_vel_covar_zt_impc, hydromet_vel_covar_zt_impc_host);
  from_execution_space(workspace, hydromet_vel_covar_zt_expc, hydromet_vel_covar_zt_expc_host);
  from_execution_space(workspace, wprtp_mc, wprtp_mc_host);
  from_execution_space(workspace, wpthlp_mc, wpthlp_mc_host);
  from_execution_space(workspace, rtp2_mc, rtp2_mc_host);
  from_execution_space(workspace, thlp2_mc, thlp2_mc_host);
  from_execution_space(workspace, rtpthlp_mc, rtpthlp_mc_host);
}
using microphys_host::KK_microphys_adjust;
} // namespace clubb_kokkos
