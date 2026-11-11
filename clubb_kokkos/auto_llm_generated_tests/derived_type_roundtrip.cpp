#include "f2k/derived_type_bridge.hpp"

#include <Kokkos_Core.hpp>

#include <cstring>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace {

template <class View>
void fill(View& view, int seed) {
  using T = typename View::non_const_value_type;
  auto host = Kokkos::create_mirror_view(view);
  int p = 0;
  const auto assign = [&](T& item) {
    if constexpr (std::is_same_v<T, char>) item = static_cast<char>('A' + (seed + p++) % 26);
    else if constexpr (std::is_same_v<T, std::uint8_t>) item = static_cast<T>((seed + p++) % 2);
    else item = static_cast<T>(seed + p++);
  };
  if constexpr (View::rank == 1) {
    for (int i = 0; i < view.extent(0); ++i) assign(host(i));
  } else if constexpr (View::rank == 2) {
    for (int j = 0; j < view.extent(1); ++j)
      for (int i = 0; i < view.extent(0); ++i) assign(host(i, j));
  } else if constexpr (View::rank == 3) {
    for (int k = 0; k < view.extent(2); ++k)
      for (int j = 0; j < view.extent(1); ++j)
        for (int i = 0; i < view.extent(0); ++i) assign(host(i, j, k));
  }
  Kokkos::deep_copy(view, host);
}

template <class View>
void same_view(const View& before, const View& after, const char* name) {
  if (before.is_allocated() != after.is_allocated()) throw std::runtime_error(name);
  if (!before.is_allocated()) return;
  for (int d = 0; d < View::rank; ++d)
    if (before.extent(d) != after.extent(d)) throw std::runtime_error(name);
  auto a = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, before);
  auto b = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, after);
  if constexpr (View::rank == 1) {
    for (int i = 0; i < before.extent(0); ++i)
      if (std::memcmp(&a(i), &b(i), sizeof(a(i)))) throw std::runtime_error(name);
  } else if constexpr (View::rank == 2) {
    for (int j = 0; j < before.extent(1); ++j)
      for (int i = 0; i < before.extent(0); ++i)
        if (std::memcmp(&a(i, j), &b(i, j), sizeof(a(i, j)))) throw std::runtime_error(name);
  } else if constexpr (View::rank == 3) {
    for (int k = 0; k < before.extent(2); ++k)
      for (int j = 0; j < before.extent(1); ++j)
        for (int i = 0; i < before.extent(0); ++i)
          if (std::memcmp(&a(i, j, k), &b(i, j, k), sizeof(a(i, j, k))))
            throw std::runtime_error(name);
  }
}

template <class T>
void same_scalar(const T& before, const T& after, const char* name) {
  if constexpr (std::is_floating_point_v<T>) {
    if (std::memcmp(&before, &after, sizeof(T))) throw std::runtime_error(name);
  } else if (before != after) {
    throw std::runtime_error(name);
  }
}

}  // namespace

int main() {
  Kokkos::initialize();
  {
// Generated fieldwise checks follow.
    {
      clubb_kokkos::grid value;
      value.nzm = 101;
      value.nzt = 102;
      value.zm = decltype(value.zm)("zm", 2, 3);
      fill(value.zm, 3);
      value.zt = decltype(value.zt)("zt", 2, 3);
      fill(value.zt, 4);
      value.invrs_dzm = decltype(value.invrs_dzm)("invrs_dzm", 2, 3);
      fill(value.invrs_dzm, 5);
      value.invrs_dzt = decltype(value.invrs_dzt)("invrs_dzt", 2, 3);
      fill(value.invrs_dzt, 6);
      value.dzm = decltype(value.dzm)("dzm", 2, 3);
      fill(value.dzm, 7);
      value.dzt = decltype(value.dzt)("dzt", 2, 3);
      fill(value.dzt, 8);
      value.weights_zm2zt = decltype(value.weights_zm2zt)("weights_zm2zt", 2, 3, 4);
      fill(value.weights_zm2zt, 9);
      value.weights_zt2zm = decltype(value.weights_zt2zm)("weights_zt2zm", 2, 3, 4);
      fill(value.weights_zt2zm, 10);
      value.k_lb_zm = 111;
      value.k_lb_zt = 112;
      value.k_ub_zm = 113;
      value.k_ub_zt = 114;
      value.grid_dir_indx = 115;
      value.grid_dir = static_cast<clubb_kokkos::real_t>(16 + 0.125);
      value.grid_dir = static_cast<clubb_kokkos::real_t>(-0.0);
      auto host_zm = Kokkos::create_mirror_view(value.zm);
      host_zm(0, 0) = static_cast<clubb_kokkos::real_t>(-0.0);
      Kokkos::deep_copy(value.zm, host_zm);
      auto before = value;
      f2k::roundtrip(value);
      same_scalar(before.nzm, value.nzm, "grid%nzm");
      same_scalar(before.nzt, value.nzt, "grid%nzt");
      same_view(before.zm, value.zm, "grid%zm");
      same_view(before.zt, value.zt, "grid%zt");
      same_view(before.invrs_dzm, value.invrs_dzm, "grid%invrs_dzm");
      same_view(before.invrs_dzt, value.invrs_dzt, "grid%invrs_dzt");
      same_view(before.dzm, value.dzm, "grid%dzm");
      same_view(before.dzt, value.dzt, "grid%dzt");
      same_view(before.weights_zm2zt, value.weights_zm2zt, "grid%weights_zm2zt");
      same_view(before.weights_zt2zm, value.weights_zt2zm, "grid%weights_zt2zm");
      same_scalar(before.k_lb_zm, value.k_lb_zm, "grid%k_lb_zm");
      same_scalar(before.k_lb_zt, value.k_lb_zt, "grid%k_lb_zt");
      same_scalar(before.k_ub_zm, value.k_ub_zm, "grid%k_ub_zm");
      same_scalar(before.k_ub_zt, value.k_ub_zt, "grid%k_ub_zt");
      same_scalar(before.grid_dir_indx, value.grid_dir_indx, "grid%grid_dir_indx");
      same_scalar(before.grid_dir, value.grid_dir, "grid%grid_dir");
      value.zm = {};
      f2k::roundtrip(value);
      if (value.zm.is_allocated())
        throw std::runtime_error("grid%zm allocation");
      value.zm = decltype(value.zm)("empty", 0, 3);
      f2k::roundtrip(value);
      if (!value.zm.is_allocated() || value.zm.extent(0) != 0)
        throw std::runtime_error("grid%zm empty shape");
    }
    {
      clubb_kokkos::sclr_idx_type defaults;
      auto packed = f2k::to_fortran(defaults);
      same_scalar(-1, defaults.iisclr_rt, "sclr_idx_type%iisclr_rt default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[0].data, sizeof(native_value));
        same_scalar(0, native_value, "sclr_idx_type%iisclr_rt default wire ordinal");
      }
      same_scalar(-1, defaults.iisclr_thl, "sclr_idx_type%iisclr_thl default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[1].data, sizeof(native_value));
        same_scalar(0, native_value, "sclr_idx_type%iisclr_thl default wire ordinal");
      }
      same_scalar(-1, defaults.iisclr_CO2, "sclr_idx_type%iisclr_CO2 default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[2].data, sizeof(native_value));
        same_scalar(0, native_value, "sclr_idx_type%iisclr_CO2 default wire ordinal");
      }
      same_scalar(-1, defaults.iiedsclr_rt, "sclr_idx_type%iiedsclr_rt default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[3].data, sizeof(native_value));
        same_scalar(0, native_value, "sclr_idx_type%iiedsclr_rt default wire ordinal");
      }
      same_scalar(-1, defaults.iiedsclr_thl, "sclr_idx_type%iiedsclr_thl default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[4].data, sizeof(native_value));
        same_scalar(0, native_value, "sclr_idx_type%iiedsclr_thl default wire ordinal");
      }
      same_scalar(-1, defaults.iiedsclr_CO2, "sclr_idx_type%iiedsclr_CO2 default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[5].data, sizeof(native_value));
        same_scalar(0, native_value, "sclr_idx_type%iiedsclr_CO2 default wire ordinal");
      }
    }
    for (int native_index : {-1, 0, 1, 2, 19}) {
      clubb_kokkos::sclr_idx_type value;
      value.iisclr_rt    = native_index - 1;
      value.iisclr_thl   = native_index - 1;
      value.iisclr_CO2   = native_index - 1;
      value.iiedsclr_rt  = native_index - 1;
      value.iiedsclr_thl = native_index - 1;
      value.iiedsclr_CO2 = native_index - 1;
      auto packed        = f2k::to_fortran(value);
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[0].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "sclr_idx_type%iisclr_rt wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[1].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "sclr_idx_type%iisclr_thl wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[2].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "sclr_idx_type%iisclr_CO2 wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[3].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "sclr_idx_type%iiedsclr_rt wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[4].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "sclr_idx_type%iiedsclr_thl wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[5].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "sclr_idx_type%iiedsclr_CO2 wire ordinal");
      }
      f2k::roundtrip(value);
      same_scalar(native_index - 1, value.iisclr_rt, "sclr_idx_type%iisclr_rt sentinel roundtrip");
      same_scalar(native_index - 1, value.iisclr_thl,
                  "sclr_idx_type%iisclr_thl sentinel roundtrip");
      same_scalar(native_index - 1, value.iisclr_CO2,
                  "sclr_idx_type%iisclr_CO2 sentinel roundtrip");
      same_scalar(native_index - 1, value.iiedsclr_rt,
                  "sclr_idx_type%iiedsclr_rt sentinel roundtrip");
      same_scalar(native_index - 1, value.iiedsclr_thl,
                  "sclr_idx_type%iiedsclr_thl sentinel roundtrip");
      same_scalar(native_index - 1, value.iiedsclr_CO2,
                  "sclr_idx_type%iiedsclr_CO2 sentinel roundtrip");
    }
    {
      clubb_kokkos::sclr_idx_type value;
      value.iisclr_rt = 101;
      value.iisclr_thl = 102;
      value.iisclr_CO2 = 103;
      value.iiedsclr_rt = 104;
      value.iiedsclr_thl = 105;
      value.iiedsclr_CO2 = 106;
      auto before = value;
      f2k::roundtrip(value);
      same_scalar(before.iisclr_rt, value.iisclr_rt, "sclr_idx_type%iisclr_rt");
      same_scalar(before.iisclr_thl, value.iisclr_thl, "sclr_idx_type%iisclr_thl");
      same_scalar(before.iisclr_CO2, value.iisclr_CO2, "sclr_idx_type%iisclr_CO2");
      same_scalar(before.iiedsclr_rt, value.iiedsclr_rt, "sclr_idx_type%iiedsclr_rt");
      same_scalar(before.iiedsclr_thl, value.iiedsclr_thl, "sclr_idx_type%iiedsclr_thl");
      same_scalar(before.iiedsclr_CO2, value.iiedsclr_CO2, "sclr_idx_type%iiedsclr_CO2");
    }
    {
      clubb_kokkos::clubb_config_flags_type value;
      value.iiPDF_type = 101;
      value.ipdf_call_placement = 102;
      value.penta_solve_method = 103;
      value.tridiag_solve_method = 104;
      value.saturation_formula = 105;
      value.Lscale_method = 106;
      value.grid_remap_method = 107;
      value.grid_adapt_in_time_method = 108;
      value.fill_holes_type = 109;
      value.l_use_precip_frac = 0;
      value.l_predict_upwp_vpwp = 1;
      value.l_ho_nontrad_coriolis = 0;
      value.l_ho_trad_coriolis = 1;
      value.l_min_wp2_from_corr_wx = 0;
      value.l_min_xp2_from_corr_wx = 1;
      value.l_C2_cloud_frac = 0;
      value.l_diffuse_rtm_and_thlm = 1;
      value.l_stability_correct_Kh_N2_zm = 0;
      value.l_calc_thlp2_rad = 1;
      value.l_upwind_xpyp_ta = 0;
      value.l_upwind_xm_ma = 1;
      value.l_uv_nudge = 0;
      value.l_rtm_nudge = 1;
      value.l_vert_avg_closure = 0;
      value.l_trapezoidal_rule_zt = 1;
      value.l_trapezoidal_rule_zm = 0;
      value.l_call_pdf_closure_twice = 1;
      value.l_standard_term_ta = 0;
      value.l_partial_upwind_wp3 = 1;
      value.l_godunov_upwind_wpxp_ta = 0;
      value.l_godunov_upwind_xpyp_ta = 1;
      value.l_use_cloud_cover = 0;
      value.l_diagnose_correlations = 1;
      value.l_calc_w_corr = 0;
      value.l_const_Nc_in_cloud = 1;
      value.l_fix_w_chi_eta_correlations = 0;
      value.l_stability_correct_tau_zm = 1;
      value.l_damp_wp2_using_em = 0;
      value.l_do_expldiff_rtm_thlm = 1;
      value.l_Lscale_plume_centered = 0;
      value.l_use_C7_Richardson = 1;
      value.l_use_C11_Richardson = 0;
      value.l_use_shear_Richardson = 1;
      value.l_brunt_vaisala_freq_moist = 0;
      value.l_use_thvm_in_bv_freq = 1;
      value.l_rcm_supersat_adj = 0;
      value.l_damp_wp3_Skw_squared = 1;
      value.l_prescribed_avg_deltaz = 0;
      value.l_lmm_stepping = 1;
      value.l_e3sm_config = 0;
      value.l_vary_convect_depth = 1;
      value.l_use_tke_in_wp3_pr_turb_term = 0;
      value.l_use_tke_in_wp2_wp3_K_dfsn = 1;
      value.l_use_wp3_lim_with_smth_Heaviside = 0;
      value.l_smooth_Heaviside_tau_wpxp = 1;
      value.l_modify_limiters_for_cnvg_test = 0;
      value.l_enable_relaxed_clipping = 1;
      value.l_linearize_pbl_winds = 0;
      value.l_mono_flux_lim_thlm = 1;
      value.l_mono_flux_lim_rtm = 0;
      value.l_mono_flux_lim_um = 1;
      value.l_mono_flux_lim_vm = 0;
      value.l_mono_flux_lim_spikefix = 1;
      value.l_host_applies_sfc_fluxes = 0;
      value.l_wp2_fill_holes_tke = 1;
      value.l_add_dycore_grid = 0;
      auto before = value;
      f2k::roundtrip(value);
      same_scalar(before.iiPDF_type, value.iiPDF_type, "clubb_config_flags_type%iiPDF_type");
      same_scalar(before.ipdf_call_placement, value.ipdf_call_placement, "clubb_config_flags_type%ipdf_call_placement");
      same_scalar(before.penta_solve_method, value.penta_solve_method, "clubb_config_flags_type%penta_solve_method");
      same_scalar(before.tridiag_solve_method, value.tridiag_solve_method, "clubb_config_flags_type%tridiag_solve_method");
      same_scalar(before.saturation_formula, value.saturation_formula, "clubb_config_flags_type%saturation_formula");
      same_scalar(before.Lscale_method, value.Lscale_method, "clubb_config_flags_type%Lscale_method");
      same_scalar(before.grid_remap_method, value.grid_remap_method, "clubb_config_flags_type%grid_remap_method");
      same_scalar(before.grid_adapt_in_time_method, value.grid_adapt_in_time_method, "clubb_config_flags_type%grid_adapt_in_time_method");
      same_scalar(before.fill_holes_type, value.fill_holes_type, "clubb_config_flags_type%fill_holes_type");
      same_scalar(before.l_use_precip_frac, value.l_use_precip_frac, "clubb_config_flags_type%l_use_precip_frac");
      same_scalar(before.l_predict_upwp_vpwp, value.l_predict_upwp_vpwp, "clubb_config_flags_type%l_predict_upwp_vpwp");
      same_scalar(before.l_ho_nontrad_coriolis, value.l_ho_nontrad_coriolis, "clubb_config_flags_type%l_ho_nontrad_coriolis");
      same_scalar(before.l_ho_trad_coriolis, value.l_ho_trad_coriolis, "clubb_config_flags_type%l_ho_trad_coriolis");
      same_scalar(before.l_min_wp2_from_corr_wx, value.l_min_wp2_from_corr_wx, "clubb_config_flags_type%l_min_wp2_from_corr_wx");
      same_scalar(before.l_min_xp2_from_corr_wx, value.l_min_xp2_from_corr_wx, "clubb_config_flags_type%l_min_xp2_from_corr_wx");
      same_scalar(before.l_C2_cloud_frac, value.l_C2_cloud_frac, "clubb_config_flags_type%l_C2_cloud_frac");
      same_scalar(before.l_diffuse_rtm_and_thlm, value.l_diffuse_rtm_and_thlm, "clubb_config_flags_type%l_diffuse_rtm_and_thlm");
      same_scalar(before.l_stability_correct_Kh_N2_zm, value.l_stability_correct_Kh_N2_zm, "clubb_config_flags_type%l_stability_correct_Kh_N2_zm");
      same_scalar(before.l_calc_thlp2_rad, value.l_calc_thlp2_rad, "clubb_config_flags_type%l_calc_thlp2_rad");
      same_scalar(before.l_upwind_xpyp_ta, value.l_upwind_xpyp_ta, "clubb_config_flags_type%l_upwind_xpyp_ta");
      same_scalar(before.l_upwind_xm_ma, value.l_upwind_xm_ma, "clubb_config_flags_type%l_upwind_xm_ma");
      same_scalar(before.l_uv_nudge, value.l_uv_nudge, "clubb_config_flags_type%l_uv_nudge");
      same_scalar(before.l_rtm_nudge, value.l_rtm_nudge, "clubb_config_flags_type%l_rtm_nudge");
      same_scalar(before.l_vert_avg_closure, value.l_vert_avg_closure, "clubb_config_flags_type%l_vert_avg_closure");
      same_scalar(before.l_trapezoidal_rule_zt, value.l_trapezoidal_rule_zt, "clubb_config_flags_type%l_trapezoidal_rule_zt");
      same_scalar(before.l_trapezoidal_rule_zm, value.l_trapezoidal_rule_zm, "clubb_config_flags_type%l_trapezoidal_rule_zm");
      same_scalar(before.l_call_pdf_closure_twice, value.l_call_pdf_closure_twice, "clubb_config_flags_type%l_call_pdf_closure_twice");
      same_scalar(before.l_standard_term_ta, value.l_standard_term_ta, "clubb_config_flags_type%l_standard_term_ta");
      same_scalar(before.l_partial_upwind_wp3, value.l_partial_upwind_wp3, "clubb_config_flags_type%l_partial_upwind_wp3");
      same_scalar(before.l_godunov_upwind_wpxp_ta, value.l_godunov_upwind_wpxp_ta, "clubb_config_flags_type%l_godunov_upwind_wpxp_ta");
      same_scalar(before.l_godunov_upwind_xpyp_ta, value.l_godunov_upwind_xpyp_ta, "clubb_config_flags_type%l_godunov_upwind_xpyp_ta");
      same_scalar(before.l_use_cloud_cover, value.l_use_cloud_cover, "clubb_config_flags_type%l_use_cloud_cover");
      same_scalar(before.l_diagnose_correlations, value.l_diagnose_correlations, "clubb_config_flags_type%l_diagnose_correlations");
      same_scalar(before.l_calc_w_corr, value.l_calc_w_corr, "clubb_config_flags_type%l_calc_w_corr");
      same_scalar(before.l_const_Nc_in_cloud, value.l_const_Nc_in_cloud, "clubb_config_flags_type%l_const_Nc_in_cloud");
      same_scalar(before.l_fix_w_chi_eta_correlations, value.l_fix_w_chi_eta_correlations, "clubb_config_flags_type%l_fix_w_chi_eta_correlations");
      same_scalar(before.l_stability_correct_tau_zm, value.l_stability_correct_tau_zm, "clubb_config_flags_type%l_stability_correct_tau_zm");
      same_scalar(before.l_damp_wp2_using_em, value.l_damp_wp2_using_em, "clubb_config_flags_type%l_damp_wp2_using_em");
      same_scalar(before.l_do_expldiff_rtm_thlm, value.l_do_expldiff_rtm_thlm, "clubb_config_flags_type%l_do_expldiff_rtm_thlm");
      same_scalar(before.l_Lscale_plume_centered, value.l_Lscale_plume_centered, "clubb_config_flags_type%l_Lscale_plume_centered");
      same_scalar(before.l_use_C7_Richardson, value.l_use_C7_Richardson, "clubb_config_flags_type%l_use_C7_Richardson");
      same_scalar(before.l_use_C11_Richardson, value.l_use_C11_Richardson, "clubb_config_flags_type%l_use_C11_Richardson");
      same_scalar(before.l_use_shear_Richardson, value.l_use_shear_Richardson, "clubb_config_flags_type%l_use_shear_Richardson");
      same_scalar(before.l_brunt_vaisala_freq_moist, value.l_brunt_vaisala_freq_moist, "clubb_config_flags_type%l_brunt_vaisala_freq_moist");
      same_scalar(before.l_use_thvm_in_bv_freq, value.l_use_thvm_in_bv_freq, "clubb_config_flags_type%l_use_thvm_in_bv_freq");
      same_scalar(before.l_rcm_supersat_adj, value.l_rcm_supersat_adj, "clubb_config_flags_type%l_rcm_supersat_adj");
      same_scalar(before.l_damp_wp3_Skw_squared, value.l_damp_wp3_Skw_squared, "clubb_config_flags_type%l_damp_wp3_Skw_squared");
      same_scalar(before.l_prescribed_avg_deltaz, value.l_prescribed_avg_deltaz, "clubb_config_flags_type%l_prescribed_avg_deltaz");
      same_scalar(before.l_lmm_stepping, value.l_lmm_stepping, "clubb_config_flags_type%l_lmm_stepping");
      same_scalar(before.l_e3sm_config, value.l_e3sm_config, "clubb_config_flags_type%l_e3sm_config");
      same_scalar(before.l_vary_convect_depth, value.l_vary_convect_depth, "clubb_config_flags_type%l_vary_convect_depth");
      same_scalar(before.l_use_tke_in_wp3_pr_turb_term, value.l_use_tke_in_wp3_pr_turb_term, "clubb_config_flags_type%l_use_tke_in_wp3_pr_turb_term");
      same_scalar(before.l_use_tke_in_wp2_wp3_K_dfsn, value.l_use_tke_in_wp2_wp3_K_dfsn, "clubb_config_flags_type%l_use_tke_in_wp2_wp3_K_dfsn");
      same_scalar(before.l_use_wp3_lim_with_smth_Heaviside, value.l_use_wp3_lim_with_smth_Heaviside, "clubb_config_flags_type%l_use_wp3_lim_with_smth_Heaviside");
      same_scalar(before.l_smooth_Heaviside_tau_wpxp, value.l_smooth_Heaviside_tau_wpxp, "clubb_config_flags_type%l_smooth_Heaviside_tau_wpxp");
      same_scalar(before.l_modify_limiters_for_cnvg_test, value.l_modify_limiters_for_cnvg_test, "clubb_config_flags_type%l_modify_limiters_for_cnvg_test");
      same_scalar(before.l_enable_relaxed_clipping, value.l_enable_relaxed_clipping, "clubb_config_flags_type%l_enable_relaxed_clipping");
      same_scalar(before.l_linearize_pbl_winds, value.l_linearize_pbl_winds, "clubb_config_flags_type%l_linearize_pbl_winds");
      same_scalar(before.l_mono_flux_lim_thlm, value.l_mono_flux_lim_thlm, "clubb_config_flags_type%l_mono_flux_lim_thlm");
      same_scalar(before.l_mono_flux_lim_rtm, value.l_mono_flux_lim_rtm, "clubb_config_flags_type%l_mono_flux_lim_rtm");
      same_scalar(before.l_mono_flux_lim_um, value.l_mono_flux_lim_um, "clubb_config_flags_type%l_mono_flux_lim_um");
      same_scalar(before.l_mono_flux_lim_vm, value.l_mono_flux_lim_vm, "clubb_config_flags_type%l_mono_flux_lim_vm");
      same_scalar(before.l_mono_flux_lim_spikefix, value.l_mono_flux_lim_spikefix, "clubb_config_flags_type%l_mono_flux_lim_spikefix");
      same_scalar(before.l_host_applies_sfc_fluxes, value.l_host_applies_sfc_fluxes, "clubb_config_flags_type%l_host_applies_sfc_fluxes");
      same_scalar(before.l_wp2_fill_holes_tke, value.l_wp2_fill_holes_tke, "clubb_config_flags_type%l_wp2_fill_holes_tke");
      same_scalar(before.l_add_dycore_grid, value.l_add_dycore_grid, "clubb_config_flags_type%l_add_dycore_grid");
    }
    {
      clubb_kokkos::nu_vertical_res_dep value;
      value.nu1 = decltype(value.nu1)("nu1", 2);
      fill(value.nu1, 1);
      value.nu2 = decltype(value.nu2)("nu2", 2);
      fill(value.nu2, 2);
      value.nu6 = decltype(value.nu6)("nu6", 2);
      fill(value.nu6, 3);
      value.nu8 = decltype(value.nu8)("nu8", 2);
      fill(value.nu8, 4);
      value.nu9 = decltype(value.nu9)("nu9", 2);
      fill(value.nu9, 5);
      value.nu10 = decltype(value.nu10)("nu10", 2);
      fill(value.nu10, 6);
      value.nu_hm = decltype(value.nu_hm)("nu_hm", 2);
      fill(value.nu_hm, 7);
      auto before = value;
      f2k::roundtrip(value);
      same_view(before.nu1, value.nu1, "nu_vertical_res_dep%nu1");
      same_view(before.nu2, value.nu2, "nu_vertical_res_dep%nu2");
      same_view(before.nu6, value.nu6, "nu_vertical_res_dep%nu6");
      same_view(before.nu8, value.nu8, "nu_vertical_res_dep%nu8");
      same_view(before.nu9, value.nu9, "nu_vertical_res_dep%nu9");
      same_view(before.nu10, value.nu10, "nu_vertical_res_dep%nu10");
      same_view(before.nu_hm, value.nu_hm, "nu_vertical_res_dep%nu_hm");
      value.nu1 = {};
      f2k::roundtrip(value);
      if (value.nu1.is_allocated())
        throw std::runtime_error("nu_vertical_res_dep%nu1 allocation");
      value.nu1 = decltype(value.nu1)("empty", 0);
      f2k::roundtrip(value);
      if (!value.nu1.is_allocated() || value.nu1.extent(0) != 0)
        throw std::runtime_error("nu_vertical_res_dep%nu1 empty shape");
    }
    {
      clubb_kokkos::pdf_parameter value;
      value.ngrdcol = 101;
      value.nz = 102;
      value.w_1 = decltype(value.w_1)("w_1", 2, 3);
      fill(value.w_1, 3);
      value.w_2 = decltype(value.w_2)("w_2", 2, 3);
      fill(value.w_2, 4);
      value.varnce_w_1 = decltype(value.varnce_w_1)("varnce_w_1", 2, 3);
      fill(value.varnce_w_1, 5);
      value.varnce_w_2 = decltype(value.varnce_w_2)("varnce_w_2", 2, 3);
      fill(value.varnce_w_2, 6);
      value.rt_1 = decltype(value.rt_1)("rt_1", 2, 3);
      fill(value.rt_1, 7);
      value.rt_2 = decltype(value.rt_2)("rt_2", 2, 3);
      fill(value.rt_2, 8);
      value.varnce_rt_1 = decltype(value.varnce_rt_1)("varnce_rt_1", 2, 3);
      fill(value.varnce_rt_1, 9);
      value.varnce_rt_2 = decltype(value.varnce_rt_2)("varnce_rt_2", 2, 3);
      fill(value.varnce_rt_2, 10);
      value.thl_1 = decltype(value.thl_1)("thl_1", 2, 3);
      fill(value.thl_1, 11);
      value.thl_2 = decltype(value.thl_2)("thl_2", 2, 3);
      fill(value.thl_2, 12);
      value.varnce_thl_1 = decltype(value.varnce_thl_1)("varnce_thl_1", 2, 3);
      fill(value.varnce_thl_1, 13);
      value.varnce_thl_2 = decltype(value.varnce_thl_2)("varnce_thl_2", 2, 3);
      fill(value.varnce_thl_2, 14);
      value.corr_w_rt_1 = decltype(value.corr_w_rt_1)("corr_w_rt_1", 2, 3);
      fill(value.corr_w_rt_1, 15);
      value.corr_w_rt_2 = decltype(value.corr_w_rt_2)("corr_w_rt_2", 2, 3);
      fill(value.corr_w_rt_2, 16);
      value.corr_w_thl_1 = decltype(value.corr_w_thl_1)("corr_w_thl_1", 2, 3);
      fill(value.corr_w_thl_1, 17);
      value.corr_w_thl_2 = decltype(value.corr_w_thl_2)("corr_w_thl_2", 2, 3);
      fill(value.corr_w_thl_2, 18);
      value.corr_rt_thl_1 = decltype(value.corr_rt_thl_1)("corr_rt_thl_1", 2, 3);
      fill(value.corr_rt_thl_1, 19);
      value.corr_rt_thl_2 = decltype(value.corr_rt_thl_2)("corr_rt_thl_2", 2, 3);
      fill(value.corr_rt_thl_2, 20);
      value.alpha_thl = decltype(value.alpha_thl)("alpha_thl", 2, 3);
      fill(value.alpha_thl, 21);
      value.alpha_rt = decltype(value.alpha_rt)("alpha_rt", 2, 3);
      fill(value.alpha_rt, 22);
      value.crt_1 = decltype(value.crt_1)("crt_1", 2, 3);
      fill(value.crt_1, 23);
      value.crt_2 = decltype(value.crt_2)("crt_2", 2, 3);
      fill(value.crt_2, 24);
      value.cthl_1 = decltype(value.cthl_1)("cthl_1", 2, 3);
      fill(value.cthl_1, 25);
      value.cthl_2 = decltype(value.cthl_2)("cthl_2", 2, 3);
      fill(value.cthl_2, 26);
      value.chi_1 = decltype(value.chi_1)("chi_1", 2, 3);
      fill(value.chi_1, 27);
      value.chi_2 = decltype(value.chi_2)("chi_2", 2, 3);
      fill(value.chi_2, 28);
      value.stdev_chi_1 = decltype(value.stdev_chi_1)("stdev_chi_1", 2, 3);
      fill(value.stdev_chi_1, 29);
      value.stdev_chi_2 = decltype(value.stdev_chi_2)("stdev_chi_2", 2, 3);
      fill(value.stdev_chi_2, 30);
      value.stdev_eta_1 = decltype(value.stdev_eta_1)("stdev_eta_1", 2, 3);
      fill(value.stdev_eta_1, 31);
      value.stdev_eta_2 = decltype(value.stdev_eta_2)("stdev_eta_2", 2, 3);
      fill(value.stdev_eta_2, 32);
      value.covar_chi_eta_1 = decltype(value.covar_chi_eta_1)("covar_chi_eta_1", 2, 3);
      fill(value.covar_chi_eta_1, 33);
      value.covar_chi_eta_2 = decltype(value.covar_chi_eta_2)("covar_chi_eta_2", 2, 3);
      fill(value.covar_chi_eta_2, 34);
      value.corr_w_chi_1 = decltype(value.corr_w_chi_1)("corr_w_chi_1", 2, 3);
      fill(value.corr_w_chi_1, 35);
      value.corr_w_chi_2 = decltype(value.corr_w_chi_2)("corr_w_chi_2", 2, 3);
      fill(value.corr_w_chi_2, 36);
      value.corr_w_eta_1 = decltype(value.corr_w_eta_1)("corr_w_eta_1", 2, 3);
      fill(value.corr_w_eta_1, 37);
      value.corr_w_eta_2 = decltype(value.corr_w_eta_2)("corr_w_eta_2", 2, 3);
      fill(value.corr_w_eta_2, 38);
      value.corr_chi_eta_1 = decltype(value.corr_chi_eta_1)("corr_chi_eta_1", 2, 3);
      fill(value.corr_chi_eta_1, 39);
      value.corr_chi_eta_2 = decltype(value.corr_chi_eta_2)("corr_chi_eta_2", 2, 3);
      fill(value.corr_chi_eta_2, 40);
      value.rsatl_1 = decltype(value.rsatl_1)("rsatl_1", 2, 3);
      fill(value.rsatl_1, 41);
      value.rsatl_2 = decltype(value.rsatl_2)("rsatl_2", 2, 3);
      fill(value.rsatl_2, 42);
      value.rc_1 = decltype(value.rc_1)("rc_1", 2, 3);
      fill(value.rc_1, 43);
      value.rc_2 = decltype(value.rc_2)("rc_2", 2, 3);
      fill(value.rc_2, 44);
      value.cloud_frac_1 = decltype(value.cloud_frac_1)("cloud_frac_1", 2, 3);
      fill(value.cloud_frac_1, 45);
      value.cloud_frac_2 = decltype(value.cloud_frac_2)("cloud_frac_2", 2, 3);
      fill(value.cloud_frac_2, 46);
      value.mixt_frac = decltype(value.mixt_frac)("mixt_frac", 2, 3);
      fill(value.mixt_frac, 47);
      value.ice_supersat_frac_1 = decltype(value.ice_supersat_frac_1)("ice_supersat_frac_1", 2, 3);
      fill(value.ice_supersat_frac_1, 48);
      value.ice_supersat_frac_2 = decltype(value.ice_supersat_frac_2)("ice_supersat_frac_2", 2, 3);
      fill(value.ice_supersat_frac_2, 49);
      auto before = value;
      f2k::roundtrip(value);
      same_scalar(before.ngrdcol, value.ngrdcol, "pdf_parameter%ngrdcol");
      same_scalar(before.nz, value.nz, "pdf_parameter%nz");
      same_view(before.w_1, value.w_1, "pdf_parameter%w_1");
      same_view(before.w_2, value.w_2, "pdf_parameter%w_2");
      same_view(before.varnce_w_1, value.varnce_w_1, "pdf_parameter%varnce_w_1");
      same_view(before.varnce_w_2, value.varnce_w_2, "pdf_parameter%varnce_w_2");
      same_view(before.rt_1, value.rt_1, "pdf_parameter%rt_1");
      same_view(before.rt_2, value.rt_2, "pdf_parameter%rt_2");
      same_view(before.varnce_rt_1, value.varnce_rt_1, "pdf_parameter%varnce_rt_1");
      same_view(before.varnce_rt_2, value.varnce_rt_2, "pdf_parameter%varnce_rt_2");
      same_view(before.thl_1, value.thl_1, "pdf_parameter%thl_1");
      same_view(before.thl_2, value.thl_2, "pdf_parameter%thl_2");
      same_view(before.varnce_thl_1, value.varnce_thl_1, "pdf_parameter%varnce_thl_1");
      same_view(before.varnce_thl_2, value.varnce_thl_2, "pdf_parameter%varnce_thl_2");
      same_view(before.corr_w_rt_1, value.corr_w_rt_1, "pdf_parameter%corr_w_rt_1");
      same_view(before.corr_w_rt_2, value.corr_w_rt_2, "pdf_parameter%corr_w_rt_2");
      same_view(before.corr_w_thl_1, value.corr_w_thl_1, "pdf_parameter%corr_w_thl_1");
      same_view(before.corr_w_thl_2, value.corr_w_thl_2, "pdf_parameter%corr_w_thl_2");
      same_view(before.corr_rt_thl_1, value.corr_rt_thl_1, "pdf_parameter%corr_rt_thl_1");
      same_view(before.corr_rt_thl_2, value.corr_rt_thl_2, "pdf_parameter%corr_rt_thl_2");
      same_view(before.alpha_thl, value.alpha_thl, "pdf_parameter%alpha_thl");
      same_view(before.alpha_rt, value.alpha_rt, "pdf_parameter%alpha_rt");
      same_view(before.crt_1, value.crt_1, "pdf_parameter%crt_1");
      same_view(before.crt_2, value.crt_2, "pdf_parameter%crt_2");
      same_view(before.cthl_1, value.cthl_1, "pdf_parameter%cthl_1");
      same_view(before.cthl_2, value.cthl_2, "pdf_parameter%cthl_2");
      same_view(before.chi_1, value.chi_1, "pdf_parameter%chi_1");
      same_view(before.chi_2, value.chi_2, "pdf_parameter%chi_2");
      same_view(before.stdev_chi_1, value.stdev_chi_1, "pdf_parameter%stdev_chi_1");
      same_view(before.stdev_chi_2, value.stdev_chi_2, "pdf_parameter%stdev_chi_2");
      same_view(before.stdev_eta_1, value.stdev_eta_1, "pdf_parameter%stdev_eta_1");
      same_view(before.stdev_eta_2, value.stdev_eta_2, "pdf_parameter%stdev_eta_2");
      same_view(before.covar_chi_eta_1, value.covar_chi_eta_1, "pdf_parameter%covar_chi_eta_1");
      same_view(before.covar_chi_eta_2, value.covar_chi_eta_2, "pdf_parameter%covar_chi_eta_2");
      same_view(before.corr_w_chi_1, value.corr_w_chi_1, "pdf_parameter%corr_w_chi_1");
      same_view(before.corr_w_chi_2, value.corr_w_chi_2, "pdf_parameter%corr_w_chi_2");
      same_view(before.corr_w_eta_1, value.corr_w_eta_1, "pdf_parameter%corr_w_eta_1");
      same_view(before.corr_w_eta_2, value.corr_w_eta_2, "pdf_parameter%corr_w_eta_2");
      same_view(before.corr_chi_eta_1, value.corr_chi_eta_1, "pdf_parameter%corr_chi_eta_1");
      same_view(before.corr_chi_eta_2, value.corr_chi_eta_2, "pdf_parameter%corr_chi_eta_2");
      same_view(before.rsatl_1, value.rsatl_1, "pdf_parameter%rsatl_1");
      same_view(before.rsatl_2, value.rsatl_2, "pdf_parameter%rsatl_2");
      same_view(before.rc_1, value.rc_1, "pdf_parameter%rc_1");
      same_view(before.rc_2, value.rc_2, "pdf_parameter%rc_2");
      same_view(before.cloud_frac_1, value.cloud_frac_1, "pdf_parameter%cloud_frac_1");
      same_view(before.cloud_frac_2, value.cloud_frac_2, "pdf_parameter%cloud_frac_2");
      same_view(before.mixt_frac, value.mixt_frac, "pdf_parameter%mixt_frac");
      same_view(before.ice_supersat_frac_1, value.ice_supersat_frac_1, "pdf_parameter%ice_supersat_frac_1");
      same_view(before.ice_supersat_frac_2, value.ice_supersat_frac_2, "pdf_parameter%ice_supersat_frac_2");
      value.w_1 = {};
      f2k::roundtrip(value);
      if (value.w_1.is_allocated())
        throw std::runtime_error("pdf_parameter%w_1 allocation");
      value.w_1 = decltype(value.w_1)("empty", 0, 3);
      f2k::roundtrip(value);
      if (!value.w_1.is_allocated() || value.w_1.extent(0) != 0)
        throw std::runtime_error("pdf_parameter%w_1 empty shape");
    }
    {
      clubb_kokkos::implicit_coefs_terms value;
      value.ngrdcol = 101;
      value.nz = 102;
      value.sclr_dim = 103;
      value.coef_wp4_implicit = decltype(value.coef_wp4_implicit)("coef_wp4_implicit", 2, 3);
      fill(value.coef_wp4_implicit, 4);
      value.coef_wp2rtp_implicit = decltype(value.coef_wp2rtp_implicit)("coef_wp2rtp_implicit", 2, 3);
      fill(value.coef_wp2rtp_implicit, 5);
      value.term_wp2rtp_explicit = decltype(value.term_wp2rtp_explicit)("term_wp2rtp_explicit", 2, 3);
      fill(value.term_wp2rtp_explicit, 6);
      value.coef_wp2thlp_implicit = decltype(value.coef_wp2thlp_implicit)("coef_wp2thlp_implicit", 2, 3);
      fill(value.coef_wp2thlp_implicit, 7);
      value.term_wp2thlp_explicit = decltype(value.term_wp2thlp_explicit)("term_wp2thlp_explicit", 2, 3);
      fill(value.term_wp2thlp_explicit, 8);
      value.coef_wp2up_implicit = decltype(value.coef_wp2up_implicit)("coef_wp2up_implicit", 2, 3);
      fill(value.coef_wp2up_implicit, 9);
      value.term_wp2up_explicit = decltype(value.term_wp2up_explicit)("term_wp2up_explicit", 2, 3);
      fill(value.term_wp2up_explicit, 10);
      value.coef_wp2vp_implicit = decltype(value.coef_wp2vp_implicit)("coef_wp2vp_implicit", 2, 3);
      fill(value.coef_wp2vp_implicit, 11);
      value.term_wp2vp_explicit = decltype(value.term_wp2vp_explicit)("term_wp2vp_explicit", 2, 3);
      fill(value.term_wp2vp_explicit, 12);
      value.coef_wprtp2_implicit = decltype(value.coef_wprtp2_implicit)("coef_wprtp2_implicit", 2, 3);
      fill(value.coef_wprtp2_implicit, 13);
      value.term_wprtp2_explicit = decltype(value.term_wprtp2_explicit)("term_wprtp2_explicit", 2, 3);
      fill(value.term_wprtp2_explicit, 14);
      value.coef_wpthlp2_implicit = decltype(value.coef_wpthlp2_implicit)("coef_wpthlp2_implicit", 2, 3);
      fill(value.coef_wpthlp2_implicit, 15);
      value.term_wpthlp2_explicit = decltype(value.term_wpthlp2_explicit)("term_wpthlp2_explicit", 2, 3);
      fill(value.term_wpthlp2_explicit, 16);
      value.coef_wprtpthlp_implicit = decltype(value.coef_wprtpthlp_implicit)("coef_wprtpthlp_implicit", 2, 3);
      fill(value.coef_wprtpthlp_implicit, 17);
      value.term_wprtpthlp_explicit = decltype(value.term_wprtpthlp_explicit)("term_wprtpthlp_explicit", 2, 3);
      fill(value.term_wprtpthlp_explicit, 18);
      value.coef_wpup2_implicit = decltype(value.coef_wpup2_implicit)("coef_wpup2_implicit", 2, 3);
      fill(value.coef_wpup2_implicit, 19);
      value.term_wpup2_explicit = decltype(value.term_wpup2_explicit)("term_wpup2_explicit", 2, 3);
      fill(value.term_wpup2_explicit, 20);
      value.coef_wpvp2_implicit = decltype(value.coef_wpvp2_implicit)("coef_wpvp2_implicit", 2, 3);
      fill(value.coef_wpvp2_implicit, 21);
      value.term_wpvp2_explicit = decltype(value.term_wpvp2_explicit)("term_wpvp2_explicit", 2, 3);
      fill(value.term_wpvp2_explicit, 22);
      value.coef_wp2sclrp_implicit = decltype(value.coef_wp2sclrp_implicit)("coef_wp2sclrp_implicit", 2, 3, 4);
      fill(value.coef_wp2sclrp_implicit, 23);
      value.term_wp2sclrp_explicit = decltype(value.term_wp2sclrp_explicit)("term_wp2sclrp_explicit", 2, 3, 4);
      fill(value.term_wp2sclrp_explicit, 24);
      value.coef_wpsclrp2_implicit = decltype(value.coef_wpsclrp2_implicit)("coef_wpsclrp2_implicit", 2, 3, 4);
      fill(value.coef_wpsclrp2_implicit, 25);
      value.term_wpsclrp2_explicit = decltype(value.term_wpsclrp2_explicit)("term_wpsclrp2_explicit", 2, 3, 4);
      fill(value.term_wpsclrp2_explicit, 26);
      value.coef_wprtpsclrp_implicit = decltype(value.coef_wprtpsclrp_implicit)("coef_wprtpsclrp_implicit", 2, 3, 4);
      fill(value.coef_wprtpsclrp_implicit, 27);
      value.term_wprtpsclrp_explicit = decltype(value.term_wprtpsclrp_explicit)("term_wprtpsclrp_explicit", 2, 3, 4);
      fill(value.term_wprtpsclrp_explicit, 28);
      value.coef_wpthlpsclrp_implicit = decltype(value.coef_wpthlpsclrp_implicit)("coef_wpthlpsclrp_implicit", 2, 3, 4);
      fill(value.coef_wpthlpsclrp_implicit, 29);
      value.term_wpthlpsclrp_explicit = decltype(value.term_wpthlpsclrp_explicit)("term_wpthlpsclrp_explicit", 2, 3, 4);
      fill(value.term_wpthlpsclrp_explicit, 30);
      auto before = value;
      f2k::roundtrip(value);
      same_scalar(before.ngrdcol, value.ngrdcol, "implicit_coefs_terms%ngrdcol");
      same_scalar(before.nz, value.nz, "implicit_coefs_terms%nz");
      same_scalar(before.sclr_dim, value.sclr_dim, "implicit_coefs_terms%sclr_dim");
      same_view(before.coef_wp4_implicit, value.coef_wp4_implicit, "implicit_coefs_terms%coef_wp4_implicit");
      same_view(before.coef_wp2rtp_implicit, value.coef_wp2rtp_implicit, "implicit_coefs_terms%coef_wp2rtp_implicit");
      same_view(before.term_wp2rtp_explicit, value.term_wp2rtp_explicit, "implicit_coefs_terms%term_wp2rtp_explicit");
      same_view(before.coef_wp2thlp_implicit, value.coef_wp2thlp_implicit, "implicit_coefs_terms%coef_wp2thlp_implicit");
      same_view(before.term_wp2thlp_explicit, value.term_wp2thlp_explicit, "implicit_coefs_terms%term_wp2thlp_explicit");
      same_view(before.coef_wp2up_implicit, value.coef_wp2up_implicit, "implicit_coefs_terms%coef_wp2up_implicit");
      same_view(before.term_wp2up_explicit, value.term_wp2up_explicit, "implicit_coefs_terms%term_wp2up_explicit");
      same_view(before.coef_wp2vp_implicit, value.coef_wp2vp_implicit, "implicit_coefs_terms%coef_wp2vp_implicit");
      same_view(before.term_wp2vp_explicit, value.term_wp2vp_explicit, "implicit_coefs_terms%term_wp2vp_explicit");
      same_view(before.coef_wprtp2_implicit, value.coef_wprtp2_implicit, "implicit_coefs_terms%coef_wprtp2_implicit");
      same_view(before.term_wprtp2_explicit, value.term_wprtp2_explicit, "implicit_coefs_terms%term_wprtp2_explicit");
      same_view(before.coef_wpthlp2_implicit, value.coef_wpthlp2_implicit, "implicit_coefs_terms%coef_wpthlp2_implicit");
      same_view(before.term_wpthlp2_explicit, value.term_wpthlp2_explicit, "implicit_coefs_terms%term_wpthlp2_explicit");
      same_view(before.coef_wprtpthlp_implicit, value.coef_wprtpthlp_implicit, "implicit_coefs_terms%coef_wprtpthlp_implicit");
      same_view(before.term_wprtpthlp_explicit, value.term_wprtpthlp_explicit, "implicit_coefs_terms%term_wprtpthlp_explicit");
      same_view(before.coef_wpup2_implicit, value.coef_wpup2_implicit, "implicit_coefs_terms%coef_wpup2_implicit");
      same_view(before.term_wpup2_explicit, value.term_wpup2_explicit, "implicit_coefs_terms%term_wpup2_explicit");
      same_view(before.coef_wpvp2_implicit, value.coef_wpvp2_implicit, "implicit_coefs_terms%coef_wpvp2_implicit");
      same_view(before.term_wpvp2_explicit, value.term_wpvp2_explicit, "implicit_coefs_terms%term_wpvp2_explicit");
      same_view(before.coef_wp2sclrp_implicit, value.coef_wp2sclrp_implicit, "implicit_coefs_terms%coef_wp2sclrp_implicit");
      same_view(before.term_wp2sclrp_explicit, value.term_wp2sclrp_explicit, "implicit_coefs_terms%term_wp2sclrp_explicit");
      same_view(before.coef_wpsclrp2_implicit, value.coef_wpsclrp2_implicit, "implicit_coefs_terms%coef_wpsclrp2_implicit");
      same_view(before.term_wpsclrp2_explicit, value.term_wpsclrp2_explicit, "implicit_coefs_terms%term_wpsclrp2_explicit");
      same_view(before.coef_wprtpsclrp_implicit, value.coef_wprtpsclrp_implicit, "implicit_coefs_terms%coef_wprtpsclrp_implicit");
      same_view(before.term_wprtpsclrp_explicit, value.term_wprtpsclrp_explicit, "implicit_coefs_terms%term_wprtpsclrp_explicit");
      same_view(before.coef_wpthlpsclrp_implicit, value.coef_wpthlpsclrp_implicit, "implicit_coefs_terms%coef_wpthlpsclrp_implicit");
      same_view(before.term_wpthlpsclrp_explicit, value.term_wpthlpsclrp_explicit, "implicit_coefs_terms%term_wpthlpsclrp_explicit");
      value.coef_wp4_implicit = {};
      f2k::roundtrip(value);
      if (value.coef_wp4_implicit.is_allocated())
        throw std::runtime_error("implicit_coefs_terms%coef_wp4_implicit allocation");
      value.coef_wp4_implicit = decltype(value.coef_wp4_implicit)("empty", 0, 3);
      f2k::roundtrip(value);
      if (!value.coef_wp4_implicit.is_allocated() || value.coef_wp4_implicit.extent(0) != 0)
        throw std::runtime_error("implicit_coefs_terms%coef_wp4_implicit empty shape");
    }
    {
      clubb_kokkos::err_info_type value;
      value.err_code = decltype(value.err_code)("err_code", 2);
      fill(value.err_code, 1);
      value.host_lat = decltype(value.host_lat)("host_lat", 2);
      fill(value.host_lat, 2);
      value.host_lon = decltype(value.host_lon)("host_lon", 2);
      fill(value.host_lon, 3);
      value.chunk_idx = 104;
      value.mpi_rank = 105;
      value.err_header = decltype(value.err_header)("err_header", 200, 2);
      fill(value.err_header, 6);
      value.err_header_global = std::string(214, 'H');
      auto before = value;
      f2k::roundtrip(value);
      same_view(before.err_code, value.err_code, "err_info_type%err_code");
      same_view(before.host_lat, value.host_lat, "err_info_type%host_lat");
      same_view(before.host_lon, value.host_lon, "err_info_type%host_lon");
      same_scalar(before.chunk_idx, value.chunk_idx, "err_info_type%chunk_idx");
      same_scalar(before.mpi_rank, value.mpi_rank, "err_info_type%mpi_rank");
      same_view(before.err_header, value.err_header, "err_info_type%err_header");
      same_scalar(before.err_header_global, value.err_header_global, "err_info_type%err_header_global");
      value.err_code = {};
      f2k::roundtrip(value);
      if (value.err_code.is_allocated())
        throw std::runtime_error("err_info_type%err_code allocation");
      value.err_code = decltype(value.err_code)("empty", 0);
      f2k::roundtrip(value);
      if (!value.err_code.is_allocated() || value.err_code.extent(0) != 0)
        throw std::runtime_error("err_info_type%err_code empty shape");
      value.err_header = {};
      f2k::roundtrip(value);
      if (value.err_header.is_allocated())
        throw std::runtime_error("err_info_type%err_header allocation");
      value.err_header = decltype(value.err_header)("empty_headers", 200, 0);
      f2k::roundtrip(value);
      if (!value.err_header.is_allocated() || value.err_header.extent(1) != 0)
        throw std::runtime_error("err_info_type%err_header empty shape");
    }
    {
      clubb_kokkos::hm_metadata_type defaults;
      auto packed = f2k::to_fortran(defaults);
      same_scalar(-1, defaults.iirr, "hm_metadata_type%iirr default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[0].data, sizeof(native_value));
        same_scalar(0, native_value, "hm_metadata_type%iirr default wire ordinal");
      }
      same_scalar(-1, defaults.iirs, "hm_metadata_type%iirs default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[1].data, sizeof(native_value));
        same_scalar(0, native_value, "hm_metadata_type%iirs default wire ordinal");
      }
      same_scalar(-1, defaults.iiri, "hm_metadata_type%iiri default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[2].data, sizeof(native_value));
        same_scalar(0, native_value, "hm_metadata_type%iiri default wire ordinal");
      }
      same_scalar(-1, defaults.iirg, "hm_metadata_type%iirg default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[3].data, sizeof(native_value));
        same_scalar(0, native_value, "hm_metadata_type%iirg default wire ordinal");
      }
      same_scalar(-1, defaults.iiNr, "hm_metadata_type%iiNr default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[4].data, sizeof(native_value));
        same_scalar(0, native_value, "hm_metadata_type%iiNr default wire ordinal");
      }
      same_scalar(-1, defaults.iiNs, "hm_metadata_type%iiNs default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[5].data, sizeof(native_value));
        same_scalar(0, native_value, "hm_metadata_type%iiNs default wire ordinal");
      }
      same_scalar(-1, defaults.iiNi, "hm_metadata_type%iiNi default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[6].data, sizeof(native_value));
        same_scalar(0, native_value, "hm_metadata_type%iiNi default wire ordinal");
      }
      same_scalar(-1, defaults.iiNg, "hm_metadata_type%iiNg default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[7].data, sizeof(native_value));
        same_scalar(0, native_value, "hm_metadata_type%iiNg default wire ordinal");
      }
      same_scalar(-2, defaults.iiPDF_chi, "hm_metadata_type%iiPDF_chi default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[12].data, sizeof(native_value));
        same_scalar(-1, native_value, "hm_metadata_type%iiPDF_chi default wire ordinal");
      }
      same_scalar(-2, defaults.iiPDF_eta, "hm_metadata_type%iiPDF_eta default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[13].data, sizeof(native_value));
        same_scalar(-1, native_value, "hm_metadata_type%iiPDF_eta default wire ordinal");
      }
      same_scalar(-2, defaults.iiPDF_w, "hm_metadata_type%iiPDF_w default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[14].data, sizeof(native_value));
        same_scalar(-1, native_value, "hm_metadata_type%iiPDF_w default wire ordinal");
      }
      same_scalar(-2, defaults.iiPDF_rr, "hm_metadata_type%iiPDF_rr default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[15].data, sizeof(native_value));
        same_scalar(-1, native_value, "hm_metadata_type%iiPDF_rr default wire ordinal");
      }
      same_scalar(-2, defaults.iiPDF_rs, "hm_metadata_type%iiPDF_rs default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[16].data, sizeof(native_value));
        same_scalar(-1, native_value, "hm_metadata_type%iiPDF_rs default wire ordinal");
      }
      same_scalar(-2, defaults.iiPDF_ri, "hm_metadata_type%iiPDF_ri default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[17].data, sizeof(native_value));
        same_scalar(-1, native_value, "hm_metadata_type%iiPDF_ri default wire ordinal");
      }
      same_scalar(-2, defaults.iiPDF_rg, "hm_metadata_type%iiPDF_rg default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[18].data, sizeof(native_value));
        same_scalar(-1, native_value, "hm_metadata_type%iiPDF_rg default wire ordinal");
      }
      same_scalar(-2, defaults.iiPDF_Nr, "hm_metadata_type%iiPDF_Nr default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[19].data, sizeof(native_value));
        same_scalar(-1, native_value, "hm_metadata_type%iiPDF_Nr default wire ordinal");
      }
      same_scalar(-2, defaults.iiPDF_Ns, "hm_metadata_type%iiPDF_Ns default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[20].data, sizeof(native_value));
        same_scalar(-1, native_value, "hm_metadata_type%iiPDF_Ns default wire ordinal");
      }
      same_scalar(-2, defaults.iiPDF_Ni, "hm_metadata_type%iiPDF_Ni default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[21].data, sizeof(native_value));
        same_scalar(-1, native_value, "hm_metadata_type%iiPDF_Ni default wire ordinal");
      }
      same_scalar(-2, defaults.iiPDF_Ng, "hm_metadata_type%iiPDF_Ng default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[22].data, sizeof(native_value));
        same_scalar(-1, native_value, "hm_metadata_type%iiPDF_Ng default wire ordinal");
      }
      same_scalar(-2, defaults.iiPDF_Ncn, "hm_metadata_type%iiPDF_Ncn default slot");
      {
        int native_value;
        std::memcpy(&native_value, packed.fields[23].data, sizeof(native_value));
        same_scalar(-1, native_value, "hm_metadata_type%iiPDF_Ncn default wire ordinal");
      }
    }
    for (int native_index : {-1, 0, 1, 2, 19}) {
      clubb_kokkos::hm_metadata_type value;
      value.iirr      = native_index - 1;
      value.iirs      = native_index - 1;
      value.iiri      = native_index - 1;
      value.iirg      = native_index - 1;
      value.iiNr      = native_index - 1;
      value.iiNs      = native_index - 1;
      value.iiNi      = native_index - 1;
      value.iiNg      = native_index - 1;
      value.iiPDF_chi = native_index - 1;
      value.iiPDF_eta = native_index - 1;
      value.iiPDF_w   = native_index - 1;
      value.iiPDF_rr  = native_index - 1;
      value.iiPDF_rs  = native_index - 1;
      value.iiPDF_ri  = native_index - 1;
      value.iiPDF_rg  = native_index - 1;
      value.iiPDF_Nr  = native_index - 1;
      value.iiPDF_Ns  = native_index - 1;
      value.iiPDF_Ni  = native_index - 1;
      value.iiPDF_Ng  = native_index - 1;
      value.iiPDF_Ncn = native_index - 1;
      auto packed     = f2k::to_fortran(value);
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[0].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iirr wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[1].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iirs wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[2].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iiri wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[3].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iirg wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[4].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iiNr wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[5].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iiNs wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[6].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iiNi wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[7].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iiNg wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[12].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iiPDF_chi wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[13].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iiPDF_eta wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[14].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iiPDF_w wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[15].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iiPDF_rr wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[16].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iiPDF_rs wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[17].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iiPDF_ri wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[18].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iiPDF_rg wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[19].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iiPDF_Nr wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[20].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iiPDF_Ns wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[21].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iiPDF_Ni wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[22].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iiPDF_Ng wire ordinal");
      }
      {
        int wire_index;
        std::memcpy(&wire_index, packed.fields[23].data, sizeof(wire_index));
        same_scalar(native_index, wire_index, "hm_metadata_type%iiPDF_Ncn wire ordinal");
      }
      f2k::roundtrip(value);
      same_scalar(native_index - 1, value.iirr, "hm_metadata_type%iirr sentinel roundtrip");
      same_scalar(native_index - 1, value.iirs, "hm_metadata_type%iirs sentinel roundtrip");
      same_scalar(native_index - 1, value.iiri, "hm_metadata_type%iiri sentinel roundtrip");
      same_scalar(native_index - 1, value.iirg, "hm_metadata_type%iirg sentinel roundtrip");
      same_scalar(native_index - 1, value.iiNr, "hm_metadata_type%iiNr sentinel roundtrip");
      same_scalar(native_index - 1, value.iiNs, "hm_metadata_type%iiNs sentinel roundtrip");
      same_scalar(native_index - 1, value.iiNi, "hm_metadata_type%iiNi sentinel roundtrip");
      same_scalar(native_index - 1, value.iiNg, "hm_metadata_type%iiNg sentinel roundtrip");
      same_scalar(native_index - 1, value.iiPDF_chi,
                  "hm_metadata_type%iiPDF_chi sentinel roundtrip");
      same_scalar(native_index - 1, value.iiPDF_eta,
                  "hm_metadata_type%iiPDF_eta sentinel roundtrip");
      same_scalar(native_index - 1, value.iiPDF_w, "hm_metadata_type%iiPDF_w sentinel roundtrip");
      same_scalar(native_index - 1, value.iiPDF_rr, "hm_metadata_type%iiPDF_rr sentinel roundtrip");
      same_scalar(native_index - 1, value.iiPDF_rs, "hm_metadata_type%iiPDF_rs sentinel roundtrip");
      same_scalar(native_index - 1, value.iiPDF_ri, "hm_metadata_type%iiPDF_ri sentinel roundtrip");
      same_scalar(native_index - 1, value.iiPDF_rg, "hm_metadata_type%iiPDF_rg sentinel roundtrip");
      same_scalar(native_index - 1, value.iiPDF_Nr, "hm_metadata_type%iiPDF_Nr sentinel roundtrip");
      same_scalar(native_index - 1, value.iiPDF_Ns, "hm_metadata_type%iiPDF_Ns sentinel roundtrip");
      same_scalar(native_index - 1, value.iiPDF_Ni, "hm_metadata_type%iiPDF_Ni sentinel roundtrip");
      same_scalar(native_index - 1, value.iiPDF_Ng, "hm_metadata_type%iiPDF_Ng sentinel roundtrip");
      same_scalar(native_index - 1, value.iiPDF_Ncn,
                  "hm_metadata_type%iiPDF_Ncn sentinel roundtrip");
    }
    {
      clubb_kokkos::hm_metadata_type value;
      value.iirr = 101;
      value.iirs = 102;
      value.iiri = 103;
      value.iirg = 104;
      value.iiNr = 105;
      value.iiNs = 106;
      value.iiNi = 107;
      value.iiNg = 108;
      value.l_frozen_hm = decltype(value.l_frozen_hm)("l_frozen_hm", 2);
      fill(value.l_frozen_hm, 9);
      value.l_mix_rat_hm = decltype(value.l_mix_rat_hm)("l_mix_rat_hm", 2);
      fill(value.l_mix_rat_hm, 10);
      value.hydromet_list = decltype(value.hydromet_list)("hydromet_list", 10, 2);
      fill(value.hydromet_list, 11);
      value.hydromet_tol = decltype(value.hydromet_tol)("hydromet_tol", 2);
      fill(value.hydromet_tol, 12);
      value.iiPDF_chi = 113;
      value.iiPDF_eta = 114;
      value.iiPDF_w = 115;
      value.iiPDF_rr = 116;
      value.iiPDF_rs = 117;
      value.iiPDF_ri = 118;
      value.iiPDF_rg = 119;
      value.iiPDF_Nr = 120;
      value.iiPDF_Ns = 121;
      value.iiPDF_Ni = 122;
      value.iiPDF_Ng = 123;
      value.iiPDF_Ncn = 124;
      value.hmp2_ip_on_hmm2_ip = decltype(value.hmp2_ip_on_hmm2_ip)("hmp2_ip_on_hmm2_ip", 2);
      fill(value.hmp2_ip_on_hmm2_ip, 25);
      value.Ncnp2_on_Ncnm2 = static_cast<clubb_kokkos::real_t>(26 + 0.125);
      auto before = value;
      f2k::roundtrip(value);
      same_scalar(before.iirr, value.iirr, "hm_metadata_type%iirr");
      same_scalar(before.iirs, value.iirs, "hm_metadata_type%iirs");
      same_scalar(before.iiri, value.iiri, "hm_metadata_type%iiri");
      same_scalar(before.iirg, value.iirg, "hm_metadata_type%iirg");
      same_scalar(before.iiNr, value.iiNr, "hm_metadata_type%iiNr");
      same_scalar(before.iiNs, value.iiNs, "hm_metadata_type%iiNs");
      same_scalar(before.iiNi, value.iiNi, "hm_metadata_type%iiNi");
      same_scalar(before.iiNg, value.iiNg, "hm_metadata_type%iiNg");
      same_view(before.l_frozen_hm, value.l_frozen_hm, "hm_metadata_type%l_frozen_hm");
      same_view(before.l_mix_rat_hm, value.l_mix_rat_hm, "hm_metadata_type%l_mix_rat_hm");
      same_view(before.hydromet_list, value.hydromet_list, "hm_metadata_type%hydromet_list");
      same_view(before.hydromet_tol, value.hydromet_tol, "hm_metadata_type%hydromet_tol");
      same_scalar(before.iiPDF_chi, value.iiPDF_chi, "hm_metadata_type%iiPDF_chi");
      same_scalar(before.iiPDF_eta, value.iiPDF_eta, "hm_metadata_type%iiPDF_eta");
      same_scalar(before.iiPDF_w, value.iiPDF_w, "hm_metadata_type%iiPDF_w");
      same_scalar(before.iiPDF_rr, value.iiPDF_rr, "hm_metadata_type%iiPDF_rr");
      same_scalar(before.iiPDF_rs, value.iiPDF_rs, "hm_metadata_type%iiPDF_rs");
      same_scalar(before.iiPDF_ri, value.iiPDF_ri, "hm_metadata_type%iiPDF_ri");
      same_scalar(before.iiPDF_rg, value.iiPDF_rg, "hm_metadata_type%iiPDF_rg");
      same_scalar(before.iiPDF_Nr, value.iiPDF_Nr, "hm_metadata_type%iiPDF_Nr");
      same_scalar(before.iiPDF_Ns, value.iiPDF_Ns, "hm_metadata_type%iiPDF_Ns");
      same_scalar(before.iiPDF_Ni, value.iiPDF_Ni, "hm_metadata_type%iiPDF_Ni");
      same_scalar(before.iiPDF_Ng, value.iiPDF_Ng, "hm_metadata_type%iiPDF_Ng");
      same_scalar(before.iiPDF_Ncn, value.iiPDF_Ncn, "hm_metadata_type%iiPDF_Ncn");
      same_view(before.hmp2_ip_on_hmm2_ip, value.hmp2_ip_on_hmm2_ip, "hm_metadata_type%hmp2_ip_on_hmm2_ip");
      same_scalar(before.Ncnp2_on_Ncnm2, value.Ncnp2_on_Ncnm2, "hm_metadata_type%Ncnp2_on_Ncnm2");
      value.l_frozen_hm = {};
      f2k::roundtrip(value);
      if (value.l_frozen_hm.is_allocated())
        throw std::runtime_error("hm_metadata_type%l_frozen_hm allocation");
      value.l_frozen_hm = decltype(value.l_frozen_hm)("empty", 0);
      f2k::roundtrip(value);
      if (!value.l_frozen_hm.is_allocated() || value.l_frozen_hm.extent(0) != 0)
        throw std::runtime_error("hm_metadata_type%l_frozen_hm empty shape");
    }
    {
      clubb_kokkos::silhs_config_flags_type value;
      value.cluster_allocation_strategy = 101;
      value.l_lh_importance_sampling = 0;
      value.l_Lscale_vert_avg = 1;
      value.l_lh_straight_mc = 0;
      value.l_lh_clustered_sampling = 1;
      value.l_rcm_in_cloud_k_lh_start = 0;
      value.l_random_k_lh_start = 1;
      value.l_max_overlap_in_cloud = 0;
      value.l_lh_instant_var_covar_src = 1;
      value.l_lh_limit_weights = 0;
      value.l_lh_var_frac = 1;
      value.l_lh_normalize_weights = 0;
      value.l_lh_deterministic_test = 1;
      auto before = value;
      f2k::roundtrip(value);
      same_scalar(before.cluster_allocation_strategy, value.cluster_allocation_strategy, "silhs_config_flags_type%cluster_allocation_strategy");
      same_scalar(before.l_lh_importance_sampling, value.l_lh_importance_sampling, "silhs_config_flags_type%l_lh_importance_sampling");
      same_scalar(before.l_Lscale_vert_avg, value.l_Lscale_vert_avg, "silhs_config_flags_type%l_Lscale_vert_avg");
      same_scalar(before.l_lh_straight_mc, value.l_lh_straight_mc, "silhs_config_flags_type%l_lh_straight_mc");
      same_scalar(before.l_lh_clustered_sampling, value.l_lh_clustered_sampling, "silhs_config_flags_type%l_lh_clustered_sampling");
      same_scalar(before.l_rcm_in_cloud_k_lh_start, value.l_rcm_in_cloud_k_lh_start, "silhs_config_flags_type%l_rcm_in_cloud_k_lh_start");
      same_scalar(before.l_random_k_lh_start, value.l_random_k_lh_start, "silhs_config_flags_type%l_random_k_lh_start");
      same_scalar(before.l_max_overlap_in_cloud, value.l_max_overlap_in_cloud, "silhs_config_flags_type%l_max_overlap_in_cloud");
      same_scalar(before.l_lh_instant_var_covar_src, value.l_lh_instant_var_covar_src, "silhs_config_flags_type%l_lh_instant_var_covar_src");
      same_scalar(before.l_lh_limit_weights, value.l_lh_limit_weights, "silhs_config_flags_type%l_lh_limit_weights");
      same_scalar(before.l_lh_var_frac, value.l_lh_var_frac, "silhs_config_flags_type%l_lh_var_frac");
      same_scalar(before.l_lh_normalize_weights, value.l_lh_normalize_weights, "silhs_config_flags_type%l_lh_normalize_weights");
      same_scalar(before.l_lh_deterministic_test, value.l_lh_deterministic_test, "silhs_config_flags_type%l_lh_deterministic_test");
    }
    {
      clubb_kokkos::hydromet_pdf_parameter value;
      value.hm_1 = decltype(value.hm_1)("hm_1", 8);
      fill(value.hm_1, 1);
      value.hm_2 = decltype(value.hm_2)("hm_2", 8);
      fill(value.hm_2, 2);
      value.mu_hm_1 = decltype(value.mu_hm_1)("mu_hm_1", 8);
      fill(value.mu_hm_1, 3);
      value.mu_hm_2 = decltype(value.mu_hm_2)("mu_hm_2", 8);
      fill(value.mu_hm_2, 4);
      value.sigma_hm_1 = decltype(value.sigma_hm_1)("sigma_hm_1", 8);
      fill(value.sigma_hm_1, 5);
      value.sigma_hm_2 = decltype(value.sigma_hm_2)("sigma_hm_2", 8);
      fill(value.sigma_hm_2, 6);
      value.corr_w_hm_1 = decltype(value.corr_w_hm_1)("corr_w_hm_1", 8);
      fill(value.corr_w_hm_1, 7);
      value.corr_w_hm_2 = decltype(value.corr_w_hm_2)("corr_w_hm_2", 8);
      fill(value.corr_w_hm_2, 8);
      value.corr_chi_hm_1 = decltype(value.corr_chi_hm_1)("corr_chi_hm_1", 8);
      fill(value.corr_chi_hm_1, 9);
      value.corr_chi_hm_2 = decltype(value.corr_chi_hm_2)("corr_chi_hm_2", 8);
      fill(value.corr_chi_hm_2, 10);
      value.corr_eta_hm_1 = decltype(value.corr_eta_hm_1)("corr_eta_hm_1", 8);
      fill(value.corr_eta_hm_1, 11);
      value.corr_eta_hm_2 = decltype(value.corr_eta_hm_2)("corr_eta_hm_2", 8);
      fill(value.corr_eta_hm_2, 12);
      value.corr_hmx_hmy_1 = decltype(value.corr_hmx_hmy_1)("corr_hmx_hmy_1", 8, 8);
      fill(value.corr_hmx_hmy_1, 13);
      value.corr_hmx_hmy_2 = decltype(value.corr_hmx_hmy_2)("corr_hmx_hmy_2", 8, 8);
      fill(value.corr_hmx_hmy_2, 14);
      value.mu_Ncn_1 = static_cast<clubb_kokkos::real_t>(15 + 0.125);
      value.mu_Ncn_2 = static_cast<clubb_kokkos::real_t>(16 + 0.125);
      value.sigma_Ncn_1 = static_cast<clubb_kokkos::real_t>(17 + 0.125);
      value.sigma_Ncn_2 = static_cast<clubb_kokkos::real_t>(18 + 0.125);
      auto before = value;
      f2k::roundtrip(value);
      same_view(before.hm_1, value.hm_1, "hydromet_pdf_parameter%hm_1");
      same_view(before.hm_2, value.hm_2, "hydromet_pdf_parameter%hm_2");
      same_view(before.mu_hm_1, value.mu_hm_1, "hydromet_pdf_parameter%mu_hm_1");
      same_view(before.mu_hm_2, value.mu_hm_2, "hydromet_pdf_parameter%mu_hm_2");
      same_view(before.sigma_hm_1, value.sigma_hm_1, "hydromet_pdf_parameter%sigma_hm_1");
      same_view(before.sigma_hm_2, value.sigma_hm_2, "hydromet_pdf_parameter%sigma_hm_2");
      same_view(before.corr_w_hm_1, value.corr_w_hm_1, "hydromet_pdf_parameter%corr_w_hm_1");
      same_view(before.corr_w_hm_2, value.corr_w_hm_2, "hydromet_pdf_parameter%corr_w_hm_2");
      same_view(before.corr_chi_hm_1, value.corr_chi_hm_1, "hydromet_pdf_parameter%corr_chi_hm_1");
      same_view(before.corr_chi_hm_2, value.corr_chi_hm_2, "hydromet_pdf_parameter%corr_chi_hm_2");
      same_view(before.corr_eta_hm_1, value.corr_eta_hm_1, "hydromet_pdf_parameter%corr_eta_hm_1");
      same_view(before.corr_eta_hm_2, value.corr_eta_hm_2, "hydromet_pdf_parameter%corr_eta_hm_2");
      same_view(before.corr_hmx_hmy_1, value.corr_hmx_hmy_1, "hydromet_pdf_parameter%corr_hmx_hmy_1");
      same_view(before.corr_hmx_hmy_2, value.corr_hmx_hmy_2, "hydromet_pdf_parameter%corr_hmx_hmy_2");
      same_scalar(before.mu_Ncn_1, value.mu_Ncn_1, "hydromet_pdf_parameter%mu_Ncn_1");
      same_scalar(before.mu_Ncn_2, value.mu_Ncn_2, "hydromet_pdf_parameter%mu_Ncn_2");
      same_scalar(before.sigma_Ncn_1, value.sigma_Ncn_1, "hydromet_pdf_parameter%sigma_Ncn_1");
      same_scalar(before.sigma_Ncn_2, value.sigma_Ncn_2, "hydromet_pdf_parameter%sigma_Ncn_2");
    }
    {
      clubb_kokkos::precipitation_fractions value;
      value.ngrdcol = 101;
      value.nzt = 102;
      value.precip_frac = decltype(value.precip_frac)("precip_frac", 2, 3);
      fill(value.precip_frac, 3);
      value.precip_frac_1 = decltype(value.precip_frac_1)("precip_frac_1", 2, 3);
      fill(value.precip_frac_1, 4);
      value.precip_frac_2 = decltype(value.precip_frac_2)("precip_frac_2", 2, 3);
      fill(value.precip_frac_2, 5);
      auto before = value;
      f2k::roundtrip(value);
      same_scalar(before.ngrdcol, value.ngrdcol, "precipitation_fractions%ngrdcol");
      same_scalar(before.nzt, value.nzt, "precipitation_fractions%nzt");
      same_view(before.precip_frac, value.precip_frac, "precipitation_fractions%precip_frac");
      same_view(before.precip_frac_1, value.precip_frac_1, "precipitation_fractions%precip_frac_1");
      same_view(before.precip_frac_2, value.precip_frac_2, "precipitation_fractions%precip_frac_2");
      value.precip_frac = {};
      f2k::roundtrip(value);
      if (value.precip_frac.is_allocated())
        throw std::runtime_error("precipitation_fractions%precip_frac allocation");
      value.precip_frac = decltype(value.precip_frac)("empty", 0, 3);
      f2k::roundtrip(value);
      if (!value.precip_frac.is_allocated() || value.precip_frac.extent(0) != 0)
        throw std::runtime_error("precipitation_fractions%precip_frac empty shape");
    }
  }
  Kokkos::finalize();
  return 0;
}
