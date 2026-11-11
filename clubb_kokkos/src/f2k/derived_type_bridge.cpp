#include "f2k/derived_type_bridge.hpp"
#include "f2k/array_bridge.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace f2k {

Transfer::Transfer(std::size_t count, bool output) : fields(count), owns_output_(output) {}

Transfer::Transfer(Transfer&& other) noexcept
    : fields(std::move(other.fields)), buffers(std::move(other.buffers)),
      owns_output_(std::exchange(other.owns_output_, false)) {}

Transfer::~Transfer() {
  if (owns_output_) {
    for (const auto& field : fields) std::free(field.data);
  }
}

namespace {

using Shape = std::array<std::int64_t, 4>;

WireField descriptor(int rank, Shape shape, int width = 0) {
  WireField field{};
  field.rank = rank;
  field.allocated = 1;
  field.char_len = width;
  for (int i = 0; i < 4; ++i) field.extent[i] = shape[i];
  return field;
}

void check(const WireField& field, int rank, Shape expected, int width, const char* name) {
  if (field.rank != rank || field.char_len != width)
    throw std::runtime_error(std::string("f2k field metadata mismatch: ") + name);
  if (!field.allocated) return;
  for (int d = 0; d < 4; ++d)
    if (field.extent[d] < 0 || (expected[d] >= 0 && field.extent[d] != expected[d]))
      throw std::runtime_error(std::string("f2k field shape mismatch: ") + name);
}

int count(const WireField& field) {
  std::int64_t n = 1;
  for (int d = 0; d < field.rank; ++d) {
    if (n && field.extent[d] > std::numeric_limits<int>::max() / n)
      throw std::length_error("f2k field is too large");
    n *= field.extent[d];
  }
  return static_cast<int>(n);
}

template <class T>
WireField pack_scalar(const T& value, Transfer& transfer) {
  auto field = descriptor(0, {1, 1, 1, 1});
  auto& bytes = transfer.buffers.emplace_back(sizeof(T));
  std::memcpy(bytes.data(), &value, sizeof(T));
  field.data = bytes.data();
  return field;
}

template <class T>
void unpack_scalar(T& value, const WireField& field, const char* name) {
  check(field, 0, {1, 1, 1, 1}, 0, name);
  if (!field.allocated || !field.data) throw std::runtime_error("f2k missing scalar");
  std::memcpy(&value, field.data, sizeof(T));
}

WireField pack_string(const std::string& value, int width, Transfer& transfer) {
  auto field = descriptor(0, {1, 1, 1, 1}, width);
  auto& bytes = transfer.buffers.emplace_back(width, ' ');
  std::memcpy(bytes.data(), value.data(), std::min<std::size_t>(width, value.size()));
  field.data = bytes.data();
  return field;
}

void unpack_string(std::string& value, const WireField& field, int width, const char* name) {
  check(field, 0, {1, 1, 1, 1}, width, name);
  if (!field.data) throw std::runtime_error("f2k missing string");
  value.assign(static_cast<const char*>(field.data), width);
}

template <class View>
WireField pack_view(const View& view, Transfer& transfer, Shape expected, const char* name) {
  auto field = pack_array(view, transfer);
  check(field, View::rank, expected, 0, name);
  if (!field.allocated && expected[0] >= 0)
    throw std::runtime_error(std::string("f2k missing fixed array: ") + name);
  return field;
}

template <class View>
void unpack_view(View& view, const WireField& field, Shape expected, const char* name) {
  check(field, View::rank, expected, 0, name);
  if (!field.allocated && expected[0] >= 0)
    throw std::runtime_error(std::string("f2k missing fixed array: ") + name);
  unpack_array(view, field, name);
}

WireField pack_char_view(const Kokkos::View<char**, Kokkos::LayoutLeft>& view,
                         int width, Transfer& transfer, const char* name) {
  auto field = descriptor(1, {0, 1, 1, 1}, width);
  if (!view.is_allocated()) {
    field.allocated = 0;
    return field;
  }
  if (view.extent(0) != static_cast<std::size_t>(width))
    throw std::runtime_error(std::string("f2k character width: ") + name);
  field.extent[0] = view.extent(1);
  const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, view);
  auto& bytes = transfer.buffers.emplace_back(static_cast<std::size_t>(count(field)) * width);
  for (int j = 0; j < field.extent[0]; ++j)
    for (int k = 0; k < width; ++k) bytes[static_cast<std::size_t>(j) * width + k] = host(k, j);
  field.data = bytes.empty() ? nullptr : bytes.data();
  return field;
}

void unpack_char_view(Kokkos::View<char**, Kokkos::LayoutLeft>& view,
                      const WireField& field, int width, const char* name) {
  check(field, 1, {-1, 1, 1, 1}, width, name);
  if (!field.allocated) {
    view = {};
    return;
  }
  if (count(field) && !field.data) throw std::runtime_error("f2k missing character data");
  view = Kokkos::View<char**, Kokkos::LayoutLeft>(name, width, field.extent[0]);
  auto host = Kokkos::create_mirror_view(view);
  const auto* bytes = static_cast<const char*>(field.data);
  for (int j = 0; j < field.extent[0]; ++j)
    for (int k = 0; k < width; ++k) host(k, j) = bytes[static_cast<std::size_t>(j) * width + k];
  Kokkos::deep_copy(view, host);
}

}  // namespace

WireField pack_real_array(const Real1D& value, Transfer& transfer) {
  // Preserve allocation, shape and first-index order at the native boundary.
  return pack_view(value, transfer, {-1, 1, 1, 1}, "real array");
}

WireField pack_real_array(const Real2D& value, Transfer& transfer) {
  // Preserve allocation, shape and first-index order at the native boundary.
  return pack_view(value, transfer, {-1, -1, 1, 1}, "real array");
}

WireField pack_strings(const std::vector<std::string>& value, Transfer& transfer) {
  // Blank-pad host metadata strings to one Fortran character width.
  std::size_t width = 1;
  for (const auto& item : value) width = std::max(width, item.size());
  if (width > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("f2k string is too long");
  }
  auto field = descriptor(1, {static_cast<std::int64_t>(value.size()), 1, 1, 1},
                          static_cast<int>(width));
  auto& bytes = transfer.buffers.emplace_back(value.size() * width, ' ');
  for (std::size_t i = 0; i < value.size(); ++i) {
    std::memcpy(bytes.data() + i * width, value[i].data(), value[i].size());
  }
  field.data = bytes.data();
  return field;
}

WireField pack_char_array(const Kokkos::View<char**, Kokkos::LayoutLeft>& value,
                          Transfer& transfer) {
  // Character arrays use the first View index for the Fortran string width.
  return pack_char_view(value, static_cast<int>(value.extent(0)), transfer, "character array");
}

// Generated per-type field lists follow.

extern "C" void f2k_roundtrip_grid_c(const WireField*, WireField*);

Transfer to_fortran(const clubb_kokkos::grid& value) {
  Transfer transfer(16);
  transfer.fields[0] = pack_scalar(value.nzm, transfer);
  transfer.fields[1] = pack_scalar(value.nzt, transfer);
  transfer.fields[2] = pack_view(value.zm, transfer, {-1, -1, 1, 1}, "zm");
  transfer.fields[3] = pack_view(value.zt, transfer, {-1, -1, 1, 1}, "zt");
  transfer.fields[4] = pack_view(value.invrs_dzm, transfer, {-1, -1, 1, 1}, "invrs_dzm");
  transfer.fields[5] = pack_view(value.invrs_dzt, transfer, {-1, -1, 1, 1}, "invrs_dzt");
  transfer.fields[6] = pack_view(value.dzm, transfer, {-1, -1, 1, 1}, "dzm");
  transfer.fields[7] = pack_view(value.dzt, transfer, {-1, -1, 1, 1}, "dzt");
  transfer.fields[8] = pack_view(value.weights_zm2zt, transfer, {-1, -1, -1, 1}, "weights_zm2zt");
  transfer.fields[9] = pack_view(value.weights_zt2zm, transfer, {-1, -1, -1, 1}, "weights_zt2zm");
  transfer.fields[10] = pack_scalar(value.k_lb_zm + 1, transfer);
  transfer.fields[11] = pack_scalar(value.k_lb_zt + 1, transfer);
  transfer.fields[12] = pack_scalar(value.k_ub_zm + 1, transfer);
  transfer.fields[13] = pack_scalar(value.k_ub_zt + 1, transfer);
  transfer.fields[14] = pack_scalar(value.grid_dir_indx, transfer);
  transfer.fields[15] = pack_scalar(value.grid_dir, transfer);
  return transfer;
}

Transfer make_output(const clubb_kokkos::grid&) {
  return Transfer(16, true);
}

void to_kokkos(clubb_kokkos::grid& value, const Transfer& output) {
  if (output.fields.size() != 16) throw std::runtime_error("f2k field count");
  unpack_scalar(value.nzm, output.fields[0], "nzm");
  unpack_scalar(value.nzt, output.fields[1], "nzt");
  unpack_view(value.zm, output.fields[2], {-1, -1, 1, 1}, "zm");
  unpack_view(value.zt, output.fields[3], {-1, -1, 1, 1}, "zt");
  unpack_view(value.invrs_dzm, output.fields[4], {-1, -1, 1, 1}, "invrs_dzm");
  unpack_view(value.invrs_dzt, output.fields[5], {-1, -1, 1, 1}, "invrs_dzt");
  unpack_view(value.dzm, output.fields[6], {-1, -1, 1, 1}, "dzm");
  unpack_view(value.dzt, output.fields[7], {-1, -1, 1, 1}, "dzt");
  unpack_view(value.weights_zm2zt, output.fields[8], {-1, -1, -1, 1}, "weights_zm2zt");
  unpack_view(value.weights_zt2zm, output.fields[9], {-1, -1, -1, 1}, "weights_zt2zm");
  unpack_scalar(value.k_lb_zm, output.fields[10], "k_lb_zm");
  value.k_lb_zm -= 1;
  unpack_scalar(value.k_lb_zt, output.fields[11], "k_lb_zt");
  value.k_lb_zt -= 1;
  unpack_scalar(value.k_ub_zm, output.fields[12], "k_ub_zm");
  value.k_ub_zm -= 1;
  unpack_scalar(value.k_ub_zt, output.fields[13], "k_ub_zt");
  value.k_ub_zt -= 1;
  unpack_scalar(value.grid_dir_indx, output.fields[14], "grid_dir_indx");
  unpack_scalar(value.grid_dir, output.fields[15], "grid_dir");
}

void roundtrip(clubb_kokkos::grid& value) {
  auto input = to_fortran(value);
  auto output = make_output(value);
  f2k_roundtrip_grid_c(input.data(), output.data());
  to_kokkos(value, output);
}

extern "C" void f2k_roundtrip_sclr_idx_type_c(const WireField*, WireField*);

Transfer to_fortran(const clubb_kokkos::sclr_idx_type& value) {
  Transfer transfer(6);
  transfer.fields[0] = pack_scalar(value.iisclr_rt + 1, transfer);
  transfer.fields[1] = pack_scalar(value.iisclr_thl + 1, transfer);
  transfer.fields[2] = pack_scalar(value.iisclr_CO2 + 1, transfer);
  transfer.fields[3] = pack_scalar(value.iiedsclr_rt + 1, transfer);
  transfer.fields[4] = pack_scalar(value.iiedsclr_thl + 1, transfer);
  transfer.fields[5] = pack_scalar(value.iiedsclr_CO2 + 1, transfer);
  return transfer;
}

Transfer make_output(const clubb_kokkos::sclr_idx_type&) {
  return Transfer(6, true);
}

void to_kokkos(clubb_kokkos::sclr_idx_type& value, const Transfer& output) {
  if (output.fields.size() != 6) throw std::runtime_error("f2k field count");
  unpack_scalar(value.iisclr_rt, output.fields[0], "iisclr_rt");
  value.iisclr_rt -= 1;
  unpack_scalar(value.iisclr_thl, output.fields[1], "iisclr_thl");
  value.iisclr_thl -= 1;
  unpack_scalar(value.iisclr_CO2, output.fields[2], "iisclr_CO2");
  value.iisclr_CO2 -= 1;
  unpack_scalar(value.iiedsclr_rt, output.fields[3], "iiedsclr_rt");
  value.iiedsclr_rt -= 1;
  unpack_scalar(value.iiedsclr_thl, output.fields[4], "iiedsclr_thl");
  value.iiedsclr_thl -= 1;
  unpack_scalar(value.iiedsclr_CO2, output.fields[5], "iiedsclr_CO2");
  value.iiedsclr_CO2 -= 1;
}

void roundtrip(clubb_kokkos::sclr_idx_type& value) {
  auto input = to_fortran(value);
  auto output = make_output(value);
  f2k_roundtrip_sclr_idx_type_c(input.data(), output.data());
  to_kokkos(value, output);
}

extern "C" void f2k_roundtrip_clubb_config_flags_type_c(const WireField*, WireField*);

Transfer to_fortran(const clubb_kokkos::clubb_config_flags_type& value) {
  Transfer transfer(66);
  transfer.fields[0] = pack_scalar(value.iiPDF_type, transfer);
  transfer.fields[1] = pack_scalar(value.ipdf_call_placement, transfer);
  transfer.fields[2] = pack_scalar(value.penta_solve_method, transfer);
  transfer.fields[3] = pack_scalar(value.tridiag_solve_method, transfer);
  transfer.fields[4] = pack_scalar(value.saturation_formula, transfer);
  transfer.fields[5] = pack_scalar(value.Lscale_method, transfer);
  transfer.fields[6] = pack_scalar(value.grid_remap_method, transfer);
  transfer.fields[7] = pack_scalar(value.grid_adapt_in_time_method, transfer);
  transfer.fields[8] = pack_scalar(value.fill_holes_type, transfer);
  transfer.fields[9] = pack_scalar(value.l_use_precip_frac, transfer);
  transfer.fields[10] = pack_scalar(value.l_predict_upwp_vpwp, transfer);
  transfer.fields[11] = pack_scalar(value.l_ho_nontrad_coriolis, transfer);
  transfer.fields[12] = pack_scalar(value.l_ho_trad_coriolis, transfer);
  transfer.fields[13] = pack_scalar(value.l_min_wp2_from_corr_wx, transfer);
  transfer.fields[14] = pack_scalar(value.l_min_xp2_from_corr_wx, transfer);
  transfer.fields[15] = pack_scalar(value.l_C2_cloud_frac, transfer);
  transfer.fields[16] = pack_scalar(value.l_diffuse_rtm_and_thlm, transfer);
  transfer.fields[17] = pack_scalar(value.l_stability_correct_Kh_N2_zm, transfer);
  transfer.fields[18] = pack_scalar(value.l_calc_thlp2_rad, transfer);
  transfer.fields[19] = pack_scalar(value.l_upwind_xpyp_ta, transfer);
  transfer.fields[20] = pack_scalar(value.l_upwind_xm_ma, transfer);
  transfer.fields[21] = pack_scalar(value.l_uv_nudge, transfer);
  transfer.fields[22] = pack_scalar(value.l_rtm_nudge, transfer);
  transfer.fields[23] = pack_scalar(value.l_vert_avg_closure, transfer);
  transfer.fields[24] = pack_scalar(value.l_trapezoidal_rule_zt, transfer);
  transfer.fields[25] = pack_scalar(value.l_trapezoidal_rule_zm, transfer);
  transfer.fields[26] = pack_scalar(value.l_call_pdf_closure_twice, transfer);
  transfer.fields[27] = pack_scalar(value.l_standard_term_ta, transfer);
  transfer.fields[28] = pack_scalar(value.l_partial_upwind_wp3, transfer);
  transfer.fields[29] = pack_scalar(value.l_godunov_upwind_wpxp_ta, transfer);
  transfer.fields[30] = pack_scalar(value.l_godunov_upwind_xpyp_ta, transfer);
  transfer.fields[31] = pack_scalar(value.l_use_cloud_cover, transfer);
  transfer.fields[32] = pack_scalar(value.l_diagnose_correlations, transfer);
  transfer.fields[33] = pack_scalar(value.l_calc_w_corr, transfer);
  transfer.fields[34] = pack_scalar(value.l_const_Nc_in_cloud, transfer);
  transfer.fields[35] = pack_scalar(value.l_fix_w_chi_eta_correlations, transfer);
  transfer.fields[36] = pack_scalar(value.l_stability_correct_tau_zm, transfer);
  transfer.fields[37] = pack_scalar(value.l_damp_wp2_using_em, transfer);
  transfer.fields[38] = pack_scalar(value.l_do_expldiff_rtm_thlm, transfer);
  transfer.fields[39] = pack_scalar(value.l_Lscale_plume_centered, transfer);
  transfer.fields[40] = pack_scalar(value.l_use_C7_Richardson, transfer);
  transfer.fields[41] = pack_scalar(value.l_use_C11_Richardson, transfer);
  transfer.fields[42] = pack_scalar(value.l_use_shear_Richardson, transfer);
  transfer.fields[43] = pack_scalar(value.l_brunt_vaisala_freq_moist, transfer);
  transfer.fields[44] = pack_scalar(value.l_use_thvm_in_bv_freq, transfer);
  transfer.fields[45] = pack_scalar(value.l_rcm_supersat_adj, transfer);
  transfer.fields[46] = pack_scalar(value.l_damp_wp3_Skw_squared, transfer);
  transfer.fields[47] = pack_scalar(value.l_prescribed_avg_deltaz, transfer);
  transfer.fields[48] = pack_scalar(value.l_lmm_stepping, transfer);
  transfer.fields[49] = pack_scalar(value.l_e3sm_config, transfer);
  transfer.fields[50] = pack_scalar(value.l_vary_convect_depth, transfer);
  transfer.fields[51] = pack_scalar(value.l_use_tke_in_wp3_pr_turb_term, transfer);
  transfer.fields[52] = pack_scalar(value.l_use_tke_in_wp2_wp3_K_dfsn, transfer);
  transfer.fields[53] = pack_scalar(value.l_use_wp3_lim_with_smth_Heaviside, transfer);
  transfer.fields[54] = pack_scalar(value.l_smooth_Heaviside_tau_wpxp, transfer);
  transfer.fields[55] = pack_scalar(value.l_modify_limiters_for_cnvg_test, transfer);
  transfer.fields[56] = pack_scalar(value.l_enable_relaxed_clipping, transfer);
  transfer.fields[57] = pack_scalar(value.l_linearize_pbl_winds, transfer);
  transfer.fields[58] = pack_scalar(value.l_mono_flux_lim_thlm, transfer);
  transfer.fields[59] = pack_scalar(value.l_mono_flux_lim_rtm, transfer);
  transfer.fields[60] = pack_scalar(value.l_mono_flux_lim_um, transfer);
  transfer.fields[61] = pack_scalar(value.l_mono_flux_lim_vm, transfer);
  transfer.fields[62] = pack_scalar(value.l_mono_flux_lim_spikefix, transfer);
  transfer.fields[63] = pack_scalar(value.l_host_applies_sfc_fluxes, transfer);
  transfer.fields[64] = pack_scalar(value.l_wp2_fill_holes_tke, transfer);
  transfer.fields[65] = pack_scalar(value.l_add_dycore_grid, transfer);
  return transfer;
}

Transfer make_output(const clubb_kokkos::clubb_config_flags_type&) {
  return Transfer(66, true);
}

void to_kokkos(clubb_kokkos::clubb_config_flags_type& value, const Transfer& output) {
  if (output.fields.size() != 66) throw std::runtime_error("f2k field count");
  unpack_scalar(value.iiPDF_type, output.fields[0], "iiPDF_type");
  unpack_scalar(value.ipdf_call_placement, output.fields[1], "ipdf_call_placement");
  unpack_scalar(value.penta_solve_method, output.fields[2], "penta_solve_method");
  unpack_scalar(value.tridiag_solve_method, output.fields[3], "tridiag_solve_method");
  unpack_scalar(value.saturation_formula, output.fields[4], "saturation_formula");
  unpack_scalar(value.Lscale_method, output.fields[5], "Lscale_method");
  unpack_scalar(value.grid_remap_method, output.fields[6], "grid_remap_method");
  unpack_scalar(value.grid_adapt_in_time_method, output.fields[7], "grid_adapt_in_time_method");
  unpack_scalar(value.fill_holes_type, output.fields[8], "fill_holes_type");
  unpack_scalar(value.l_use_precip_frac, output.fields[9], "l_use_precip_frac");
  unpack_scalar(value.l_predict_upwp_vpwp, output.fields[10], "l_predict_upwp_vpwp");
  unpack_scalar(value.l_ho_nontrad_coriolis, output.fields[11], "l_ho_nontrad_coriolis");
  unpack_scalar(value.l_ho_trad_coriolis, output.fields[12], "l_ho_trad_coriolis");
  unpack_scalar(value.l_min_wp2_from_corr_wx, output.fields[13], "l_min_wp2_from_corr_wx");
  unpack_scalar(value.l_min_xp2_from_corr_wx, output.fields[14], "l_min_xp2_from_corr_wx");
  unpack_scalar(value.l_C2_cloud_frac, output.fields[15], "l_C2_cloud_frac");
  unpack_scalar(value.l_diffuse_rtm_and_thlm, output.fields[16], "l_diffuse_rtm_and_thlm");
  unpack_scalar(value.l_stability_correct_Kh_N2_zm, output.fields[17], "l_stability_correct_Kh_N2_zm");
  unpack_scalar(value.l_calc_thlp2_rad, output.fields[18], "l_calc_thlp2_rad");
  unpack_scalar(value.l_upwind_xpyp_ta, output.fields[19], "l_upwind_xpyp_ta");
  unpack_scalar(value.l_upwind_xm_ma, output.fields[20], "l_upwind_xm_ma");
  unpack_scalar(value.l_uv_nudge, output.fields[21], "l_uv_nudge");
  unpack_scalar(value.l_rtm_nudge, output.fields[22], "l_rtm_nudge");
  unpack_scalar(value.l_vert_avg_closure, output.fields[23], "l_vert_avg_closure");
  unpack_scalar(value.l_trapezoidal_rule_zt, output.fields[24], "l_trapezoidal_rule_zt");
  unpack_scalar(value.l_trapezoidal_rule_zm, output.fields[25], "l_trapezoidal_rule_zm");
  unpack_scalar(value.l_call_pdf_closure_twice, output.fields[26], "l_call_pdf_closure_twice");
  unpack_scalar(value.l_standard_term_ta, output.fields[27], "l_standard_term_ta");
  unpack_scalar(value.l_partial_upwind_wp3, output.fields[28], "l_partial_upwind_wp3");
  unpack_scalar(value.l_godunov_upwind_wpxp_ta, output.fields[29], "l_godunov_upwind_wpxp_ta");
  unpack_scalar(value.l_godunov_upwind_xpyp_ta, output.fields[30], "l_godunov_upwind_xpyp_ta");
  unpack_scalar(value.l_use_cloud_cover, output.fields[31], "l_use_cloud_cover");
  unpack_scalar(value.l_diagnose_correlations, output.fields[32], "l_diagnose_correlations");
  unpack_scalar(value.l_calc_w_corr, output.fields[33], "l_calc_w_corr");
  unpack_scalar(value.l_const_Nc_in_cloud, output.fields[34], "l_const_Nc_in_cloud");
  unpack_scalar(value.l_fix_w_chi_eta_correlations, output.fields[35], "l_fix_w_chi_eta_correlations");
  unpack_scalar(value.l_stability_correct_tau_zm, output.fields[36], "l_stability_correct_tau_zm");
  unpack_scalar(value.l_damp_wp2_using_em, output.fields[37], "l_damp_wp2_using_em");
  unpack_scalar(value.l_do_expldiff_rtm_thlm, output.fields[38], "l_do_expldiff_rtm_thlm");
  unpack_scalar(value.l_Lscale_plume_centered, output.fields[39], "l_Lscale_plume_centered");
  unpack_scalar(value.l_use_C7_Richardson, output.fields[40], "l_use_C7_Richardson");
  unpack_scalar(value.l_use_C11_Richardson, output.fields[41], "l_use_C11_Richardson");
  unpack_scalar(value.l_use_shear_Richardson, output.fields[42], "l_use_shear_Richardson");
  unpack_scalar(value.l_brunt_vaisala_freq_moist, output.fields[43], "l_brunt_vaisala_freq_moist");
  unpack_scalar(value.l_use_thvm_in_bv_freq, output.fields[44], "l_use_thvm_in_bv_freq");
  unpack_scalar(value.l_rcm_supersat_adj, output.fields[45], "l_rcm_supersat_adj");
  unpack_scalar(value.l_damp_wp3_Skw_squared, output.fields[46], "l_damp_wp3_Skw_squared");
  unpack_scalar(value.l_prescribed_avg_deltaz, output.fields[47], "l_prescribed_avg_deltaz");
  unpack_scalar(value.l_lmm_stepping, output.fields[48], "l_lmm_stepping");
  unpack_scalar(value.l_e3sm_config, output.fields[49], "l_e3sm_config");
  unpack_scalar(value.l_vary_convect_depth, output.fields[50], "l_vary_convect_depth");
  unpack_scalar(value.l_use_tke_in_wp3_pr_turb_term, output.fields[51], "l_use_tke_in_wp3_pr_turb_term");
  unpack_scalar(value.l_use_tke_in_wp2_wp3_K_dfsn, output.fields[52], "l_use_tke_in_wp2_wp3_K_dfsn");
  unpack_scalar(value.l_use_wp3_lim_with_smth_Heaviside, output.fields[53], "l_use_wp3_lim_with_smth_Heaviside");
  unpack_scalar(value.l_smooth_Heaviside_tau_wpxp, output.fields[54], "l_smooth_Heaviside_tau_wpxp");
  unpack_scalar(value.l_modify_limiters_for_cnvg_test, output.fields[55], "l_modify_limiters_for_cnvg_test");
  unpack_scalar(value.l_enable_relaxed_clipping, output.fields[56], "l_enable_relaxed_clipping");
  unpack_scalar(value.l_linearize_pbl_winds, output.fields[57], "l_linearize_pbl_winds");
  unpack_scalar(value.l_mono_flux_lim_thlm, output.fields[58], "l_mono_flux_lim_thlm");
  unpack_scalar(value.l_mono_flux_lim_rtm, output.fields[59], "l_mono_flux_lim_rtm");
  unpack_scalar(value.l_mono_flux_lim_um, output.fields[60], "l_mono_flux_lim_um");
  unpack_scalar(value.l_mono_flux_lim_vm, output.fields[61], "l_mono_flux_lim_vm");
  unpack_scalar(value.l_mono_flux_lim_spikefix, output.fields[62], "l_mono_flux_lim_spikefix");
  unpack_scalar(value.l_host_applies_sfc_fluxes, output.fields[63], "l_host_applies_sfc_fluxes");
  unpack_scalar(value.l_wp2_fill_holes_tke, output.fields[64], "l_wp2_fill_holes_tke");
  unpack_scalar(value.l_add_dycore_grid, output.fields[65], "l_add_dycore_grid");
}

void roundtrip(clubb_kokkos::clubb_config_flags_type& value) {
  auto input = to_fortran(value);
  auto output = make_output(value);
  f2k_roundtrip_clubb_config_flags_type_c(input.data(), output.data());
  to_kokkos(value, output);
}

extern "C" void f2k_roundtrip_nu_vertical_res_dep_c(const WireField*, WireField*);

Transfer to_fortran(const clubb_kokkos::nu_vertical_res_dep& value) {
  Transfer transfer(7);
  transfer.fields[0] = pack_view(value.nu1, transfer, {-1, 1, 1, 1}, "nu1");
  transfer.fields[1] = pack_view(value.nu2, transfer, {-1, 1, 1, 1}, "nu2");
  transfer.fields[2] = pack_view(value.nu6, transfer, {-1, 1, 1, 1}, "nu6");
  transfer.fields[3] = pack_view(value.nu8, transfer, {-1, 1, 1, 1}, "nu8");
  transfer.fields[4] = pack_view(value.nu9, transfer, {-1, 1, 1, 1}, "nu9");
  transfer.fields[5] = pack_view(value.nu10, transfer, {-1, 1, 1, 1}, "nu10");
  transfer.fields[6] = pack_view(value.nu_hm, transfer, {-1, 1, 1, 1}, "nu_hm");
  return transfer;
}

Transfer make_output(const clubb_kokkos::nu_vertical_res_dep&) {
  return Transfer(7, true);
}

void to_kokkos(clubb_kokkos::nu_vertical_res_dep& value, const Transfer& output) {
  if (output.fields.size() != 7) throw std::runtime_error("f2k field count");
  unpack_view(value.nu1, output.fields[0], {-1, 1, 1, 1}, "nu1");
  unpack_view(value.nu2, output.fields[1], {-1, 1, 1, 1}, "nu2");
  unpack_view(value.nu6, output.fields[2], {-1, 1, 1, 1}, "nu6");
  unpack_view(value.nu8, output.fields[3], {-1, 1, 1, 1}, "nu8");
  unpack_view(value.nu9, output.fields[4], {-1, 1, 1, 1}, "nu9");
  unpack_view(value.nu10, output.fields[5], {-1, 1, 1, 1}, "nu10");
  unpack_view(value.nu_hm, output.fields[6], {-1, 1, 1, 1}, "nu_hm");
}

void roundtrip(clubb_kokkos::nu_vertical_res_dep& value) {
  auto input = to_fortran(value);
  auto output = make_output(value);
  f2k_roundtrip_nu_vertical_res_dep_c(input.data(), output.data());
  to_kokkos(value, output);
}

extern "C" void f2k_roundtrip_pdf_parameter_c(const WireField*, WireField*);

Transfer to_fortran(const clubb_kokkos::pdf_parameter& value) {
  Transfer transfer(49);
  transfer.fields[0] = pack_scalar(value.ngrdcol, transfer);
  transfer.fields[1] = pack_scalar(value.nz, transfer);
  transfer.fields[2] = pack_view(value.w_1, transfer, {-1, -1, 1, 1}, "w_1");
  transfer.fields[3] = pack_view(value.w_2, transfer, {-1, -1, 1, 1}, "w_2");
  transfer.fields[4] = pack_view(value.varnce_w_1, transfer, {-1, -1, 1, 1}, "varnce_w_1");
  transfer.fields[5] = pack_view(value.varnce_w_2, transfer, {-1, -1, 1, 1}, "varnce_w_2");
  transfer.fields[6] = pack_view(value.rt_1, transfer, {-1, -1, 1, 1}, "rt_1");
  transfer.fields[7] = pack_view(value.rt_2, transfer, {-1, -1, 1, 1}, "rt_2");
  transfer.fields[8] = pack_view(value.varnce_rt_1, transfer, {-1, -1, 1, 1}, "varnce_rt_1");
  transfer.fields[9] = pack_view(value.varnce_rt_2, transfer, {-1, -1, 1, 1}, "varnce_rt_2");
  transfer.fields[10] = pack_view(value.thl_1, transfer, {-1, -1, 1, 1}, "thl_1");
  transfer.fields[11] = pack_view(value.thl_2, transfer, {-1, -1, 1, 1}, "thl_2");
  transfer.fields[12] = pack_view(value.varnce_thl_1, transfer, {-1, -1, 1, 1}, "varnce_thl_1");
  transfer.fields[13] = pack_view(value.varnce_thl_2, transfer, {-1, -1, 1, 1}, "varnce_thl_2");
  transfer.fields[14] = pack_view(value.corr_w_rt_1, transfer, {-1, -1, 1, 1}, "corr_w_rt_1");
  transfer.fields[15] = pack_view(value.corr_w_rt_2, transfer, {-1, -1, 1, 1}, "corr_w_rt_2");
  transfer.fields[16] = pack_view(value.corr_w_thl_1, transfer, {-1, -1, 1, 1}, "corr_w_thl_1");
  transfer.fields[17] = pack_view(value.corr_w_thl_2, transfer, {-1, -1, 1, 1}, "corr_w_thl_2");
  transfer.fields[18] = pack_view(value.corr_rt_thl_1, transfer, {-1, -1, 1, 1}, "corr_rt_thl_1");
  transfer.fields[19] = pack_view(value.corr_rt_thl_2, transfer, {-1, -1, 1, 1}, "corr_rt_thl_2");
  transfer.fields[20] = pack_view(value.alpha_thl, transfer, {-1, -1, 1, 1}, "alpha_thl");
  transfer.fields[21] = pack_view(value.alpha_rt, transfer, {-1, -1, 1, 1}, "alpha_rt");
  transfer.fields[22] = pack_view(value.crt_1, transfer, {-1, -1, 1, 1}, "crt_1");
  transfer.fields[23] = pack_view(value.crt_2, transfer, {-1, -1, 1, 1}, "crt_2");
  transfer.fields[24] = pack_view(value.cthl_1, transfer, {-1, -1, 1, 1}, "cthl_1");
  transfer.fields[25] = pack_view(value.cthl_2, transfer, {-1, -1, 1, 1}, "cthl_2");
  transfer.fields[26] = pack_view(value.chi_1, transfer, {-1, -1, 1, 1}, "chi_1");
  transfer.fields[27] = pack_view(value.chi_2, transfer, {-1, -1, 1, 1}, "chi_2");
  transfer.fields[28] = pack_view(value.stdev_chi_1, transfer, {-1, -1, 1, 1}, "stdev_chi_1");
  transfer.fields[29] = pack_view(value.stdev_chi_2, transfer, {-1, -1, 1, 1}, "stdev_chi_2");
  transfer.fields[30] = pack_view(value.stdev_eta_1, transfer, {-1, -1, 1, 1}, "stdev_eta_1");
  transfer.fields[31] = pack_view(value.stdev_eta_2, transfer, {-1, -1, 1, 1}, "stdev_eta_2");
  transfer.fields[32] = pack_view(value.covar_chi_eta_1, transfer, {-1, -1, 1, 1}, "covar_chi_eta_1");
  transfer.fields[33] = pack_view(value.covar_chi_eta_2, transfer, {-1, -1, 1, 1}, "covar_chi_eta_2");
  transfer.fields[34] = pack_view(value.corr_w_chi_1, transfer, {-1, -1, 1, 1}, "corr_w_chi_1");
  transfer.fields[35] = pack_view(value.corr_w_chi_2, transfer, {-1, -1, 1, 1}, "corr_w_chi_2");
  transfer.fields[36] = pack_view(value.corr_w_eta_1, transfer, {-1, -1, 1, 1}, "corr_w_eta_1");
  transfer.fields[37] = pack_view(value.corr_w_eta_2, transfer, {-1, -1, 1, 1}, "corr_w_eta_2");
  transfer.fields[38] = pack_view(value.corr_chi_eta_1, transfer, {-1, -1, 1, 1}, "corr_chi_eta_1");
  transfer.fields[39] = pack_view(value.corr_chi_eta_2, transfer, {-1, -1, 1, 1}, "corr_chi_eta_2");
  transfer.fields[40] = pack_view(value.rsatl_1, transfer, {-1, -1, 1, 1}, "rsatl_1");
  transfer.fields[41] = pack_view(value.rsatl_2, transfer, {-1, -1, 1, 1}, "rsatl_2");
  transfer.fields[42] = pack_view(value.rc_1, transfer, {-1, -1, 1, 1}, "rc_1");
  transfer.fields[43] = pack_view(value.rc_2, transfer, {-1, -1, 1, 1}, "rc_2");
  transfer.fields[44] = pack_view(value.cloud_frac_1, transfer, {-1, -1, 1, 1}, "cloud_frac_1");
  transfer.fields[45] = pack_view(value.cloud_frac_2, transfer, {-1, -1, 1, 1}, "cloud_frac_2");
  transfer.fields[46] = pack_view(value.mixt_frac, transfer, {-1, -1, 1, 1}, "mixt_frac");
  transfer.fields[47] = pack_view(value.ice_supersat_frac_1, transfer, {-1, -1, 1, 1}, "ice_supersat_frac_1");
  transfer.fields[48] = pack_view(value.ice_supersat_frac_2, transfer, {-1, -1, 1, 1}, "ice_supersat_frac_2");
  return transfer;
}

Transfer make_output(const clubb_kokkos::pdf_parameter&) {
  return Transfer(49, true);
}

void to_kokkos(clubb_kokkos::pdf_parameter& value, const Transfer& output) {
  if (output.fields.size() != 49) throw std::runtime_error("f2k field count");
  unpack_scalar(value.ngrdcol, output.fields[0], "ngrdcol");
  unpack_scalar(value.nz, output.fields[1], "nz");
  unpack_view(value.w_1, output.fields[2], {-1, -1, 1, 1}, "w_1");
  unpack_view(value.w_2, output.fields[3], {-1, -1, 1, 1}, "w_2");
  unpack_view(value.varnce_w_1, output.fields[4], {-1, -1, 1, 1}, "varnce_w_1");
  unpack_view(value.varnce_w_2, output.fields[5], {-1, -1, 1, 1}, "varnce_w_2");
  unpack_view(value.rt_1, output.fields[6], {-1, -1, 1, 1}, "rt_1");
  unpack_view(value.rt_2, output.fields[7], {-1, -1, 1, 1}, "rt_2");
  unpack_view(value.varnce_rt_1, output.fields[8], {-1, -1, 1, 1}, "varnce_rt_1");
  unpack_view(value.varnce_rt_2, output.fields[9], {-1, -1, 1, 1}, "varnce_rt_2");
  unpack_view(value.thl_1, output.fields[10], {-1, -1, 1, 1}, "thl_1");
  unpack_view(value.thl_2, output.fields[11], {-1, -1, 1, 1}, "thl_2");
  unpack_view(value.varnce_thl_1, output.fields[12], {-1, -1, 1, 1}, "varnce_thl_1");
  unpack_view(value.varnce_thl_2, output.fields[13], {-1, -1, 1, 1}, "varnce_thl_2");
  unpack_view(value.corr_w_rt_1, output.fields[14], {-1, -1, 1, 1}, "corr_w_rt_1");
  unpack_view(value.corr_w_rt_2, output.fields[15], {-1, -1, 1, 1}, "corr_w_rt_2");
  unpack_view(value.corr_w_thl_1, output.fields[16], {-1, -1, 1, 1}, "corr_w_thl_1");
  unpack_view(value.corr_w_thl_2, output.fields[17], {-1, -1, 1, 1}, "corr_w_thl_2");
  unpack_view(value.corr_rt_thl_1, output.fields[18], {-1, -1, 1, 1}, "corr_rt_thl_1");
  unpack_view(value.corr_rt_thl_2, output.fields[19], {-1, -1, 1, 1}, "corr_rt_thl_2");
  unpack_view(value.alpha_thl, output.fields[20], {-1, -1, 1, 1}, "alpha_thl");
  unpack_view(value.alpha_rt, output.fields[21], {-1, -1, 1, 1}, "alpha_rt");
  unpack_view(value.crt_1, output.fields[22], {-1, -1, 1, 1}, "crt_1");
  unpack_view(value.crt_2, output.fields[23], {-1, -1, 1, 1}, "crt_2");
  unpack_view(value.cthl_1, output.fields[24], {-1, -1, 1, 1}, "cthl_1");
  unpack_view(value.cthl_2, output.fields[25], {-1, -1, 1, 1}, "cthl_2");
  unpack_view(value.chi_1, output.fields[26], {-1, -1, 1, 1}, "chi_1");
  unpack_view(value.chi_2, output.fields[27], {-1, -1, 1, 1}, "chi_2");
  unpack_view(value.stdev_chi_1, output.fields[28], {-1, -1, 1, 1}, "stdev_chi_1");
  unpack_view(value.stdev_chi_2, output.fields[29], {-1, -1, 1, 1}, "stdev_chi_2");
  unpack_view(value.stdev_eta_1, output.fields[30], {-1, -1, 1, 1}, "stdev_eta_1");
  unpack_view(value.stdev_eta_2, output.fields[31], {-1, -1, 1, 1}, "stdev_eta_2");
  unpack_view(value.covar_chi_eta_1, output.fields[32], {-1, -1, 1, 1}, "covar_chi_eta_1");
  unpack_view(value.covar_chi_eta_2, output.fields[33], {-1, -1, 1, 1}, "covar_chi_eta_2");
  unpack_view(value.corr_w_chi_1, output.fields[34], {-1, -1, 1, 1}, "corr_w_chi_1");
  unpack_view(value.corr_w_chi_2, output.fields[35], {-1, -1, 1, 1}, "corr_w_chi_2");
  unpack_view(value.corr_w_eta_1, output.fields[36], {-1, -1, 1, 1}, "corr_w_eta_1");
  unpack_view(value.corr_w_eta_2, output.fields[37], {-1, -1, 1, 1}, "corr_w_eta_2");
  unpack_view(value.corr_chi_eta_1, output.fields[38], {-1, -1, 1, 1}, "corr_chi_eta_1");
  unpack_view(value.corr_chi_eta_2, output.fields[39], {-1, -1, 1, 1}, "corr_chi_eta_2");
  unpack_view(value.rsatl_1, output.fields[40], {-1, -1, 1, 1}, "rsatl_1");
  unpack_view(value.rsatl_2, output.fields[41], {-1, -1, 1, 1}, "rsatl_2");
  unpack_view(value.rc_1, output.fields[42], {-1, -1, 1, 1}, "rc_1");
  unpack_view(value.rc_2, output.fields[43], {-1, -1, 1, 1}, "rc_2");
  unpack_view(value.cloud_frac_1, output.fields[44], {-1, -1, 1, 1}, "cloud_frac_1");
  unpack_view(value.cloud_frac_2, output.fields[45], {-1, -1, 1, 1}, "cloud_frac_2");
  unpack_view(value.mixt_frac, output.fields[46], {-1, -1, 1, 1}, "mixt_frac");
  unpack_view(value.ice_supersat_frac_1, output.fields[47], {-1, -1, 1, 1}, "ice_supersat_frac_1");
  unpack_view(value.ice_supersat_frac_2, output.fields[48], {-1, -1, 1, 1}, "ice_supersat_frac_2");
}

void roundtrip(clubb_kokkos::pdf_parameter& value) {
  auto input = to_fortran(value);
  auto output = make_output(value);
  f2k_roundtrip_pdf_parameter_c(input.data(), output.data());
  to_kokkos(value, output);
}

extern "C" void f2k_roundtrip_implicit_coefs_terms_c(const WireField*, WireField*);

Transfer to_fortran(const clubb_kokkos::implicit_coefs_terms& value) {
  Transfer transfer(30);
  transfer.fields[0] = pack_scalar(value.ngrdcol, transfer);
  transfer.fields[1] = pack_scalar(value.nz, transfer);
  transfer.fields[2] = pack_scalar(value.sclr_dim, transfer);
  transfer.fields[3] = pack_view(value.coef_wp4_implicit, transfer, {-1, -1, 1, 1}, "coef_wp4_implicit");
  transfer.fields[4] = pack_view(value.coef_wp2rtp_implicit, transfer, {-1, -1, 1, 1}, "coef_wp2rtp_implicit");
  transfer.fields[5] = pack_view(value.term_wp2rtp_explicit, transfer, {-1, -1, 1, 1}, "term_wp2rtp_explicit");
  transfer.fields[6] = pack_view(value.coef_wp2thlp_implicit, transfer, {-1, -1, 1, 1}, "coef_wp2thlp_implicit");
  transfer.fields[7] = pack_view(value.term_wp2thlp_explicit, transfer, {-1, -1, 1, 1}, "term_wp2thlp_explicit");
  transfer.fields[8] = pack_view(value.coef_wp2up_implicit, transfer, {-1, -1, 1, 1}, "coef_wp2up_implicit");
  transfer.fields[9] = pack_view(value.term_wp2up_explicit, transfer, {-1, -1, 1, 1}, "term_wp2up_explicit");
  transfer.fields[10] = pack_view(value.coef_wp2vp_implicit, transfer, {-1, -1, 1, 1}, "coef_wp2vp_implicit");
  transfer.fields[11] = pack_view(value.term_wp2vp_explicit, transfer, {-1, -1, 1, 1}, "term_wp2vp_explicit");
  transfer.fields[12] = pack_view(value.coef_wprtp2_implicit, transfer, {-1, -1, 1, 1}, "coef_wprtp2_implicit");
  transfer.fields[13] = pack_view(value.term_wprtp2_explicit, transfer, {-1, -1, 1, 1}, "term_wprtp2_explicit");
  transfer.fields[14] = pack_view(value.coef_wpthlp2_implicit, transfer, {-1, -1, 1, 1}, "coef_wpthlp2_implicit");
  transfer.fields[15] = pack_view(value.term_wpthlp2_explicit, transfer, {-1, -1, 1, 1}, "term_wpthlp2_explicit");
  transfer.fields[16] = pack_view(value.coef_wprtpthlp_implicit, transfer, {-1, -1, 1, 1}, "coef_wprtpthlp_implicit");
  transfer.fields[17] = pack_view(value.term_wprtpthlp_explicit, transfer, {-1, -1, 1, 1}, "term_wprtpthlp_explicit");
  transfer.fields[18] = pack_view(value.coef_wpup2_implicit, transfer, {-1, -1, 1, 1}, "coef_wpup2_implicit");
  transfer.fields[19] = pack_view(value.term_wpup2_explicit, transfer, {-1, -1, 1, 1}, "term_wpup2_explicit");
  transfer.fields[20] = pack_view(value.coef_wpvp2_implicit, transfer, {-1, -1, 1, 1}, "coef_wpvp2_implicit");
  transfer.fields[21] = pack_view(value.term_wpvp2_explicit, transfer, {-1, -1, 1, 1}, "term_wpvp2_explicit");
  transfer.fields[22] = pack_view(value.coef_wp2sclrp_implicit, transfer, {-1, -1, -1, 1}, "coef_wp2sclrp_implicit");
  transfer.fields[23] = pack_view(value.term_wp2sclrp_explicit, transfer, {-1, -1, -1, 1}, "term_wp2sclrp_explicit");
  transfer.fields[24] = pack_view(value.coef_wpsclrp2_implicit, transfer, {-1, -1, -1, 1}, "coef_wpsclrp2_implicit");
  transfer.fields[25] = pack_view(value.term_wpsclrp2_explicit, transfer, {-1, -1, -1, 1}, "term_wpsclrp2_explicit");
  transfer.fields[26] = pack_view(value.coef_wprtpsclrp_implicit, transfer, {-1, -1, -1, 1}, "coef_wprtpsclrp_implicit");
  transfer.fields[27] = pack_view(value.term_wprtpsclrp_explicit, transfer, {-1, -1, -1, 1}, "term_wprtpsclrp_explicit");
  transfer.fields[28] = pack_view(value.coef_wpthlpsclrp_implicit, transfer, {-1, -1, -1, 1}, "coef_wpthlpsclrp_implicit");
  transfer.fields[29] = pack_view(value.term_wpthlpsclrp_explicit, transfer, {-1, -1, -1, 1}, "term_wpthlpsclrp_explicit");
  return transfer;
}

Transfer make_output(const clubb_kokkos::implicit_coefs_terms&) {
  return Transfer(30, true);
}

void to_kokkos(clubb_kokkos::implicit_coefs_terms& value, const Transfer& output) {
  if (output.fields.size() != 30) throw std::runtime_error("f2k field count");
  unpack_scalar(value.ngrdcol, output.fields[0], "ngrdcol");
  unpack_scalar(value.nz, output.fields[1], "nz");
  unpack_scalar(value.sclr_dim, output.fields[2], "sclr_dim");
  unpack_view(value.coef_wp4_implicit, output.fields[3], {-1, -1, 1, 1}, "coef_wp4_implicit");
  unpack_view(value.coef_wp2rtp_implicit, output.fields[4], {-1, -1, 1, 1}, "coef_wp2rtp_implicit");
  unpack_view(value.term_wp2rtp_explicit, output.fields[5], {-1, -1, 1, 1}, "term_wp2rtp_explicit");
  unpack_view(value.coef_wp2thlp_implicit, output.fields[6], {-1, -1, 1, 1}, "coef_wp2thlp_implicit");
  unpack_view(value.term_wp2thlp_explicit, output.fields[7], {-1, -1, 1, 1}, "term_wp2thlp_explicit");
  unpack_view(value.coef_wp2up_implicit, output.fields[8], {-1, -1, 1, 1}, "coef_wp2up_implicit");
  unpack_view(value.term_wp2up_explicit, output.fields[9], {-1, -1, 1, 1}, "term_wp2up_explicit");
  unpack_view(value.coef_wp2vp_implicit, output.fields[10], {-1, -1, 1, 1}, "coef_wp2vp_implicit");
  unpack_view(value.term_wp2vp_explicit, output.fields[11], {-1, -1, 1, 1}, "term_wp2vp_explicit");
  unpack_view(value.coef_wprtp2_implicit, output.fields[12], {-1, -1, 1, 1}, "coef_wprtp2_implicit");
  unpack_view(value.term_wprtp2_explicit, output.fields[13], {-1, -1, 1, 1}, "term_wprtp2_explicit");
  unpack_view(value.coef_wpthlp2_implicit, output.fields[14], {-1, -1, 1, 1}, "coef_wpthlp2_implicit");
  unpack_view(value.term_wpthlp2_explicit, output.fields[15], {-1, -1, 1, 1}, "term_wpthlp2_explicit");
  unpack_view(value.coef_wprtpthlp_implicit, output.fields[16], {-1, -1, 1, 1}, "coef_wprtpthlp_implicit");
  unpack_view(value.term_wprtpthlp_explicit, output.fields[17], {-1, -1, 1, 1}, "term_wprtpthlp_explicit");
  unpack_view(value.coef_wpup2_implicit, output.fields[18], {-1, -1, 1, 1}, "coef_wpup2_implicit");
  unpack_view(value.term_wpup2_explicit, output.fields[19], {-1, -1, 1, 1}, "term_wpup2_explicit");
  unpack_view(value.coef_wpvp2_implicit, output.fields[20], {-1, -1, 1, 1}, "coef_wpvp2_implicit");
  unpack_view(value.term_wpvp2_explicit, output.fields[21], {-1, -1, 1, 1}, "term_wpvp2_explicit");
  unpack_view(value.coef_wp2sclrp_implicit, output.fields[22], {-1, -1, -1, 1}, "coef_wp2sclrp_implicit");
  unpack_view(value.term_wp2sclrp_explicit, output.fields[23], {-1, -1, -1, 1}, "term_wp2sclrp_explicit");
  unpack_view(value.coef_wpsclrp2_implicit, output.fields[24], {-1, -1, -1, 1}, "coef_wpsclrp2_implicit");
  unpack_view(value.term_wpsclrp2_explicit, output.fields[25], {-1, -1, -1, 1}, "term_wpsclrp2_explicit");
  unpack_view(value.coef_wprtpsclrp_implicit, output.fields[26], {-1, -1, -1, 1}, "coef_wprtpsclrp_implicit");
  unpack_view(value.term_wprtpsclrp_explicit, output.fields[27], {-1, -1, -1, 1}, "term_wprtpsclrp_explicit");
  unpack_view(value.coef_wpthlpsclrp_implicit, output.fields[28], {-1, -1, -1, 1}, "coef_wpthlpsclrp_implicit");
  unpack_view(value.term_wpthlpsclrp_explicit, output.fields[29], {-1, -1, -1, 1}, "term_wpthlpsclrp_explicit");
}

void roundtrip(clubb_kokkos::implicit_coefs_terms& value) {
  auto input = to_fortran(value);
  auto output = make_output(value);
  f2k_roundtrip_implicit_coefs_terms_c(input.data(), output.data());
  to_kokkos(value, output);
}

extern "C" void f2k_roundtrip_err_info_type_c(const WireField*, WireField*);

Transfer to_fortran(const clubb_kokkos::err_info_type& value) {
  Transfer transfer(7);
  transfer.fields[0] = pack_view(value.err_code, transfer, {-1, 1, 1, 1}, "err_code");
  transfer.fields[1] = pack_view(value.host_lat, transfer, {-1, 1, 1, 1}, "host_lat");
  transfer.fields[2] = pack_view(value.host_lon, transfer, {-1, 1, 1, 1}, "host_lon");
  transfer.fields[3] = pack_scalar(value.chunk_idx, transfer);
  transfer.fields[4] = pack_scalar(value.mpi_rank, transfer);
  transfer.fields[5] = pack_char_view(value.err_header, 200, transfer, "err_header");
  transfer.fields[6] = pack_string(value.err_header_global, 214, transfer);
  return transfer;
}

Transfer make_output(const clubb_kokkos::err_info_type&) {
  return Transfer(7, true);
}

void to_kokkos(clubb_kokkos::err_info_type& value, const Transfer& output) {
  if (output.fields.size() != 7) throw std::runtime_error("f2k field count");
  unpack_view(value.err_code, output.fields[0], {-1, 1, 1, 1}, "err_code");
  unpack_view(value.host_lat, output.fields[1], {-1, 1, 1, 1}, "host_lat");
  unpack_view(value.host_lon, output.fields[2], {-1, 1, 1, 1}, "host_lon");
  unpack_scalar(value.chunk_idx, output.fields[3], "chunk_idx");
  unpack_scalar(value.mpi_rank, output.fields[4], "mpi_rank");
  unpack_char_view(value.err_header, output.fields[5], 200, "err_header");
  unpack_string(value.err_header_global, output.fields[6], 214, "err_header_global");
}

void roundtrip(clubb_kokkos::err_info_type& value) {
  auto input = to_fortran(value);
  auto output = make_output(value);
  f2k_roundtrip_err_info_type_c(input.data(), output.data());
  to_kokkos(value, output);
}

extern "C" void f2k_roundtrip_hm_metadata_type_c(const WireField*, WireField*);

Transfer to_fortran(const clubb_kokkos::hm_metadata_type& value) {
  Transfer transfer(26);
  transfer.fields[0] = pack_scalar(value.iirr + 1, transfer);
  transfer.fields[1] = pack_scalar(value.iirs + 1, transfer);
  transfer.fields[2] = pack_scalar(value.iiri + 1, transfer);
  transfer.fields[3] = pack_scalar(value.iirg + 1, transfer);
  transfer.fields[4] = pack_scalar(value.iiNr + 1, transfer);
  transfer.fields[5] = pack_scalar(value.iiNs + 1, transfer);
  transfer.fields[6] = pack_scalar(value.iiNi + 1, transfer);
  transfer.fields[7] = pack_scalar(value.iiNg + 1, transfer);
  transfer.fields[8] = pack_view(value.l_frozen_hm, transfer, {-1, 1, 1, 1}, "l_frozen_hm");
  transfer.fields[9] = pack_view(value.l_mix_rat_hm, transfer, {-1, 1, 1, 1}, "l_mix_rat_hm");
  transfer.fields[10] = pack_char_view(value.hydromet_list, 10, transfer, "hydromet_list");
  transfer.fields[11] = pack_view(value.hydromet_tol, transfer, {-1, 1, 1, 1}, "hydromet_tol");
  transfer.fields[12] = pack_scalar(value.iiPDF_chi + 1, transfer);
  transfer.fields[13] = pack_scalar(value.iiPDF_eta + 1, transfer);
  transfer.fields[14] = pack_scalar(value.iiPDF_w + 1, transfer);
  transfer.fields[15] = pack_scalar(value.iiPDF_rr + 1, transfer);
  transfer.fields[16] = pack_scalar(value.iiPDF_rs + 1, transfer);
  transfer.fields[17] = pack_scalar(value.iiPDF_ri + 1, transfer);
  transfer.fields[18] = pack_scalar(value.iiPDF_rg + 1, transfer);
  transfer.fields[19] = pack_scalar(value.iiPDF_Nr + 1, transfer);
  transfer.fields[20] = pack_scalar(value.iiPDF_Ns + 1, transfer);
  transfer.fields[21] = pack_scalar(value.iiPDF_Ni + 1, transfer);
  transfer.fields[22] = pack_scalar(value.iiPDF_Ng + 1, transfer);
  transfer.fields[23] = pack_scalar(value.iiPDF_Ncn + 1, transfer);
  transfer.fields[24] = pack_view(value.hmp2_ip_on_hmm2_ip, transfer, {-1, 1, 1, 1}, "hmp2_ip_on_hmm2_ip");
  transfer.fields[25] = pack_scalar(value.Ncnp2_on_Ncnm2, transfer);
  return transfer;
}

Transfer make_output(const clubb_kokkos::hm_metadata_type&) {
  return Transfer(26, true);
}

void to_kokkos(clubb_kokkos::hm_metadata_type& value, const Transfer& output) {
  if (output.fields.size() != 26) throw std::runtime_error("f2k field count");
  unpack_scalar(value.iirr, output.fields[0], "iirr");
  value.iirr -= 1;
  unpack_scalar(value.iirs, output.fields[1], "iirs");
  value.iirs -= 1;
  unpack_scalar(value.iiri, output.fields[2], "iiri");
  value.iiri -= 1;
  unpack_scalar(value.iirg, output.fields[3], "iirg");
  value.iirg -= 1;
  unpack_scalar(value.iiNr, output.fields[4], "iiNr");
  value.iiNr -= 1;
  unpack_scalar(value.iiNs, output.fields[5], "iiNs");
  value.iiNs -= 1;
  unpack_scalar(value.iiNi, output.fields[6], "iiNi");
  value.iiNi -= 1;
  unpack_scalar(value.iiNg, output.fields[7], "iiNg");
  value.iiNg -= 1;
  unpack_view(value.l_frozen_hm, output.fields[8], {-1, 1, 1, 1}, "l_frozen_hm");
  unpack_view(value.l_mix_rat_hm, output.fields[9], {-1, 1, 1, 1}, "l_mix_rat_hm");
  unpack_char_view(value.hydromet_list, output.fields[10], 10, "hydromet_list");
  unpack_view(value.hydromet_tol, output.fields[11], {-1, 1, 1, 1}, "hydromet_tol");
  unpack_scalar(value.iiPDF_chi, output.fields[12], "iiPDF_chi");
  value.iiPDF_chi -= 1;
  unpack_scalar(value.iiPDF_eta, output.fields[13], "iiPDF_eta");
  value.iiPDF_eta -= 1;
  unpack_scalar(value.iiPDF_w, output.fields[14], "iiPDF_w");
  value.iiPDF_w -= 1;
  unpack_scalar(value.iiPDF_rr, output.fields[15], "iiPDF_rr");
  value.iiPDF_rr -= 1;
  unpack_scalar(value.iiPDF_rs, output.fields[16], "iiPDF_rs");
  value.iiPDF_rs -= 1;
  unpack_scalar(value.iiPDF_ri, output.fields[17], "iiPDF_ri");
  value.iiPDF_ri -= 1;
  unpack_scalar(value.iiPDF_rg, output.fields[18], "iiPDF_rg");
  value.iiPDF_rg -= 1;
  unpack_scalar(value.iiPDF_Nr, output.fields[19], "iiPDF_Nr");
  value.iiPDF_Nr -= 1;
  unpack_scalar(value.iiPDF_Ns, output.fields[20], "iiPDF_Ns");
  value.iiPDF_Ns -= 1;
  unpack_scalar(value.iiPDF_Ni, output.fields[21], "iiPDF_Ni");
  value.iiPDF_Ni -= 1;
  unpack_scalar(value.iiPDF_Ng, output.fields[22], "iiPDF_Ng");
  value.iiPDF_Ng -= 1;
  unpack_scalar(value.iiPDF_Ncn, output.fields[23], "iiPDF_Ncn");
  value.iiPDF_Ncn -= 1;
  unpack_view(value.hmp2_ip_on_hmm2_ip, output.fields[24], {-1, 1, 1, 1}, "hmp2_ip_on_hmm2_ip");
  unpack_scalar(value.Ncnp2_on_Ncnm2, output.fields[25], "Ncnp2_on_Ncnm2");
}

void roundtrip(clubb_kokkos::hm_metadata_type& value) {
  auto input = to_fortran(value);
  auto output = make_output(value);
  f2k_roundtrip_hm_metadata_type_c(input.data(), output.data());
  to_kokkos(value, output);
}

extern "C" void f2k_roundtrip_silhs_config_flags_type_c(const WireField*, WireField*);

Transfer to_fortran(const clubb_kokkos::silhs_config_flags_type& value) {
  Transfer transfer(13);
  transfer.fields[0] = pack_scalar(value.cluster_allocation_strategy, transfer);
  transfer.fields[1] = pack_scalar(value.l_lh_importance_sampling, transfer);
  transfer.fields[2] = pack_scalar(value.l_Lscale_vert_avg, transfer);
  transfer.fields[3] = pack_scalar(value.l_lh_straight_mc, transfer);
  transfer.fields[4] = pack_scalar(value.l_lh_clustered_sampling, transfer);
  transfer.fields[5] = pack_scalar(value.l_rcm_in_cloud_k_lh_start, transfer);
  transfer.fields[6] = pack_scalar(value.l_random_k_lh_start, transfer);
  transfer.fields[7] = pack_scalar(value.l_max_overlap_in_cloud, transfer);
  transfer.fields[8] = pack_scalar(value.l_lh_instant_var_covar_src, transfer);
  transfer.fields[9] = pack_scalar(value.l_lh_limit_weights, transfer);
  transfer.fields[10] = pack_scalar(value.l_lh_var_frac, transfer);
  transfer.fields[11] = pack_scalar(value.l_lh_normalize_weights, transfer);
  transfer.fields[12] = pack_scalar(value.l_lh_deterministic_test, transfer);
  return transfer;
}

Transfer make_output(const clubb_kokkos::silhs_config_flags_type&) {
  return Transfer(13, true);
}

void to_kokkos(clubb_kokkos::silhs_config_flags_type& value, const Transfer& output) {
  if (output.fields.size() != 13) throw std::runtime_error("f2k field count");
  unpack_scalar(value.cluster_allocation_strategy, output.fields[0], "cluster_allocation_strategy");
  unpack_scalar(value.l_lh_importance_sampling, output.fields[1], "l_lh_importance_sampling");
  unpack_scalar(value.l_Lscale_vert_avg, output.fields[2], "l_Lscale_vert_avg");
  unpack_scalar(value.l_lh_straight_mc, output.fields[3], "l_lh_straight_mc");
  unpack_scalar(value.l_lh_clustered_sampling, output.fields[4], "l_lh_clustered_sampling");
  unpack_scalar(value.l_rcm_in_cloud_k_lh_start, output.fields[5], "l_rcm_in_cloud_k_lh_start");
  unpack_scalar(value.l_random_k_lh_start, output.fields[6], "l_random_k_lh_start");
  unpack_scalar(value.l_max_overlap_in_cloud, output.fields[7], "l_max_overlap_in_cloud");
  unpack_scalar(value.l_lh_instant_var_covar_src, output.fields[8], "l_lh_instant_var_covar_src");
  unpack_scalar(value.l_lh_limit_weights, output.fields[9], "l_lh_limit_weights");
  unpack_scalar(value.l_lh_var_frac, output.fields[10], "l_lh_var_frac");
  unpack_scalar(value.l_lh_normalize_weights, output.fields[11], "l_lh_normalize_weights");
  unpack_scalar(value.l_lh_deterministic_test, output.fields[12], "l_lh_deterministic_test");
}

void roundtrip(clubb_kokkos::silhs_config_flags_type& value) {
  auto input = to_fortran(value);
  auto output = make_output(value);
  f2k_roundtrip_silhs_config_flags_type_c(input.data(), output.data());
  to_kokkos(value, output);
}

extern "C" void f2k_roundtrip_hydromet_pdf_parameter_c(const WireField*, WireField*);

Transfer to_fortran(const clubb_kokkos::hydromet_pdf_parameter& value) {
  Transfer transfer(18);
  transfer.fields[0] = pack_view(value.hm_1, transfer, {8, 1, 1, 1}, "hm_1");
  transfer.fields[1] = pack_view(value.hm_2, transfer, {8, 1, 1, 1}, "hm_2");
  transfer.fields[2] = pack_view(value.mu_hm_1, transfer, {8, 1, 1, 1}, "mu_hm_1");
  transfer.fields[3] = pack_view(value.mu_hm_2, transfer, {8, 1, 1, 1}, "mu_hm_2");
  transfer.fields[4] = pack_view(value.sigma_hm_1, transfer, {8, 1, 1, 1}, "sigma_hm_1");
  transfer.fields[5] = pack_view(value.sigma_hm_2, transfer, {8, 1, 1, 1}, "sigma_hm_2");
  transfer.fields[6] = pack_view(value.corr_w_hm_1, transfer, {8, 1, 1, 1}, "corr_w_hm_1");
  transfer.fields[7] = pack_view(value.corr_w_hm_2, transfer, {8, 1, 1, 1}, "corr_w_hm_2");
  transfer.fields[8] = pack_view(value.corr_chi_hm_1, transfer, {8, 1, 1, 1}, "corr_chi_hm_1");
  transfer.fields[9] = pack_view(value.corr_chi_hm_2, transfer, {8, 1, 1, 1}, "corr_chi_hm_2");
  transfer.fields[10] = pack_view(value.corr_eta_hm_1, transfer, {8, 1, 1, 1}, "corr_eta_hm_1");
  transfer.fields[11] = pack_view(value.corr_eta_hm_2, transfer, {8, 1, 1, 1}, "corr_eta_hm_2");
  transfer.fields[12] = pack_view(value.corr_hmx_hmy_1, transfer, {8, 8, 1, 1}, "corr_hmx_hmy_1");
  transfer.fields[13] = pack_view(value.corr_hmx_hmy_2, transfer, {8, 8, 1, 1}, "corr_hmx_hmy_2");
  transfer.fields[14] = pack_scalar(value.mu_Ncn_1, transfer);
  transfer.fields[15] = pack_scalar(value.mu_Ncn_2, transfer);
  transfer.fields[16] = pack_scalar(value.sigma_Ncn_1, transfer);
  transfer.fields[17] = pack_scalar(value.sigma_Ncn_2, transfer);
  return transfer;
}

Transfer make_output(const clubb_kokkos::hydromet_pdf_parameter&) {
  return Transfer(18, true);
}

void to_kokkos(clubb_kokkos::hydromet_pdf_parameter& value, const Transfer& output) {
  if (output.fields.size() != 18) throw std::runtime_error("f2k field count");
  unpack_view(value.hm_1, output.fields[0], {8, 1, 1, 1}, "hm_1");
  unpack_view(value.hm_2, output.fields[1], {8, 1, 1, 1}, "hm_2");
  unpack_view(value.mu_hm_1, output.fields[2], {8, 1, 1, 1}, "mu_hm_1");
  unpack_view(value.mu_hm_2, output.fields[3], {8, 1, 1, 1}, "mu_hm_2");
  unpack_view(value.sigma_hm_1, output.fields[4], {8, 1, 1, 1}, "sigma_hm_1");
  unpack_view(value.sigma_hm_2, output.fields[5], {8, 1, 1, 1}, "sigma_hm_2");
  unpack_view(value.corr_w_hm_1, output.fields[6], {8, 1, 1, 1}, "corr_w_hm_1");
  unpack_view(value.corr_w_hm_2, output.fields[7], {8, 1, 1, 1}, "corr_w_hm_2");
  unpack_view(value.corr_chi_hm_1, output.fields[8], {8, 1, 1, 1}, "corr_chi_hm_1");
  unpack_view(value.corr_chi_hm_2, output.fields[9], {8, 1, 1, 1}, "corr_chi_hm_2");
  unpack_view(value.corr_eta_hm_1, output.fields[10], {8, 1, 1, 1}, "corr_eta_hm_1");
  unpack_view(value.corr_eta_hm_2, output.fields[11], {8, 1, 1, 1}, "corr_eta_hm_2");
  unpack_view(value.corr_hmx_hmy_1, output.fields[12], {8, 8, 1, 1}, "corr_hmx_hmy_1");
  unpack_view(value.corr_hmx_hmy_2, output.fields[13], {8, 8, 1, 1}, "corr_hmx_hmy_2");
  unpack_scalar(value.mu_Ncn_1, output.fields[14], "mu_Ncn_1");
  unpack_scalar(value.mu_Ncn_2, output.fields[15], "mu_Ncn_2");
  unpack_scalar(value.sigma_Ncn_1, output.fields[16], "sigma_Ncn_1");
  unpack_scalar(value.sigma_Ncn_2, output.fields[17], "sigma_Ncn_2");
}

void roundtrip(clubb_kokkos::hydromet_pdf_parameter& value) {
  auto input = to_fortran(value);
  auto output = make_output(value);
  f2k_roundtrip_hydromet_pdf_parameter_c(input.data(), output.data());
  to_kokkos(value, output);
}

extern "C" void f2k_roundtrip_precipitation_fractions_c(const WireField*, WireField*);

Transfer to_fortran(const clubb_kokkos::precipitation_fractions& value) {
  Transfer transfer(5);
  transfer.fields[0] = pack_scalar(value.ngrdcol, transfer);
  transfer.fields[1] = pack_scalar(value.nzt, transfer);
  transfer.fields[2] = pack_view(value.precip_frac, transfer, {-1, -1, 1, 1}, "precip_frac");
  transfer.fields[3] = pack_view(value.precip_frac_1, transfer, {-1, -1, 1, 1}, "precip_frac_1");
  transfer.fields[4] = pack_view(value.precip_frac_2, transfer, {-1, -1, 1, 1}, "precip_frac_2");
  return transfer;
}

Transfer make_output(const clubb_kokkos::precipitation_fractions&) {
  return Transfer(5, true);
}

void to_kokkos(clubb_kokkos::precipitation_fractions& value, const Transfer& output) {
  if (output.fields.size() != 5) throw std::runtime_error("f2k field count");
  unpack_scalar(value.ngrdcol, output.fields[0], "ngrdcol");
  unpack_scalar(value.nzt, output.fields[1], "nzt");
  unpack_view(value.precip_frac, output.fields[2], {-1, -1, 1, 1}, "precip_frac");
  unpack_view(value.precip_frac_1, output.fields[3], {-1, -1, 1, 1}, "precip_frac_1");
  unpack_view(value.precip_frac_2, output.fields[4], {-1, -1, 1, 1}, "precip_frac_2");
}

void roundtrip(clubb_kokkos::precipitation_fractions& value) {
  auto input = to_fortran(value);
  auto output = make_output(value);
  f2k_roundtrip_precipitation_fractions_c(input.data(), output.data());
  to_kokkos(value, output);
}

extern "C" void f2k_roundtrip_hydromet_pdf_parameter_collection_c(const WireField*, WireField*);

Transfer to_fortran(const clubb_kokkos::hydromet_pdf_parameter_collection& value) {
  Transfer transfer(18);
  transfer.fields[0] = pack_view(value.hm_1, transfer, {-1, -1, 8, 1}, "hm_1");
  transfer.fields[1] = pack_view(value.hm_2, transfer, {-1, -1, 8, 1}, "hm_2");
  transfer.fields[2] = pack_view(value.mu_hm_1, transfer, {-1, -1, 8, 1}, "mu_hm_1");
  transfer.fields[3] = pack_view(value.mu_hm_2, transfer, {-1, -1, 8, 1}, "mu_hm_2");
  transfer.fields[4] = pack_view(value.sigma_hm_1, transfer, {-1, -1, 8, 1}, "sigma_hm_1");
  transfer.fields[5] = pack_view(value.sigma_hm_2, transfer, {-1, -1, 8, 1}, "sigma_hm_2");
  transfer.fields[6] = pack_view(value.corr_w_hm_1, transfer, {-1, -1, 8, 1}, "corr_w_hm_1");
  transfer.fields[7] = pack_view(value.corr_w_hm_2, transfer, {-1, -1, 8, 1}, "corr_w_hm_2");
  transfer.fields[8] = pack_view(value.corr_chi_hm_1, transfer, {-1, -1, 8, 1}, "corr_chi_hm_1");
  transfer.fields[9] = pack_view(value.corr_chi_hm_2, transfer, {-1, -1, 8, 1}, "corr_chi_hm_2");
  transfer.fields[10] = pack_view(value.corr_eta_hm_1, transfer, {-1, -1, 8, 1}, "corr_eta_hm_1");
  transfer.fields[11] = pack_view(value.corr_eta_hm_2, transfer, {-1, -1, 8, 1}, "corr_eta_hm_2");
  transfer.fields[12] = pack_view(value.corr_hmx_hmy_1, transfer, {-1, -1, 8, 8}, "corr_hmx_hmy_1");
  transfer.fields[13] = pack_view(value.corr_hmx_hmy_2, transfer, {-1, -1, 8, 8}, "corr_hmx_hmy_2");
  transfer.fields[14] = pack_view(value.mu_Ncn_1, transfer, {-1, -1, 1, 1}, "mu_Ncn_1");
  transfer.fields[15] = pack_view(value.mu_Ncn_2, transfer, {-1, -1, 1, 1}, "mu_Ncn_2");
  transfer.fields[16] = pack_view(value.sigma_Ncn_1, transfer, {-1, -1, 1, 1}, "sigma_Ncn_1");
  transfer.fields[17] = pack_view(value.sigma_Ncn_2, transfer, {-1, -1, 1, 1}, "sigma_Ncn_2");
  return transfer;
}

Transfer make_output(const clubb_kokkos::hydromet_pdf_parameter_collection&) {
  return Transfer(18, true);
}

void to_kokkos(clubb_kokkos::hydromet_pdf_parameter_collection& value, const Transfer& output) {
  if (output.fields.size() != 18) throw std::runtime_error("f2k collection field count");
  unpack_view(value.hm_1, output.fields[0], {-1, -1, 8, 1}, "hm_1");
  unpack_view(value.hm_2, output.fields[1], {-1, -1, 8, 1}, "hm_2");
  unpack_view(value.mu_hm_1, output.fields[2], {-1, -1, 8, 1}, "mu_hm_1");
  unpack_view(value.mu_hm_2, output.fields[3], {-1, -1, 8, 1}, "mu_hm_2");
  unpack_view(value.sigma_hm_1, output.fields[4], {-1, -1, 8, 1}, "sigma_hm_1");
  unpack_view(value.sigma_hm_2, output.fields[5], {-1, -1, 8, 1}, "sigma_hm_2");
  unpack_view(value.corr_w_hm_1, output.fields[6], {-1, -1, 8, 1}, "corr_w_hm_1");
  unpack_view(value.corr_w_hm_2, output.fields[7], {-1, -1, 8, 1}, "corr_w_hm_2");
  unpack_view(value.corr_chi_hm_1, output.fields[8], {-1, -1, 8, 1}, "corr_chi_hm_1");
  unpack_view(value.corr_chi_hm_2, output.fields[9], {-1, -1, 8, 1}, "corr_chi_hm_2");
  unpack_view(value.corr_eta_hm_1, output.fields[10], {-1, -1, 8, 1}, "corr_eta_hm_1");
  unpack_view(value.corr_eta_hm_2, output.fields[11], {-1, -1, 8, 1}, "corr_eta_hm_2");
  unpack_view(value.corr_hmx_hmy_1, output.fields[12], {-1, -1, 8, 8}, "corr_hmx_hmy_1");
  unpack_view(value.corr_hmx_hmy_2, output.fields[13], {-1, -1, 8, 8}, "corr_hmx_hmy_2");
  unpack_view(value.mu_Ncn_1, output.fields[14], {-1, -1, 1, 1}, "mu_Ncn_1");
  unpack_view(value.mu_Ncn_2, output.fields[15], {-1, -1, 1, 1}, "mu_Ncn_2");
  unpack_view(value.sigma_Ncn_1, output.fields[16], {-1, -1, 1, 1}, "sigma_Ncn_1");
  unpack_view(value.sigma_Ncn_2, output.fields[17], {-1, -1, 1, 1}, "sigma_Ncn_2");
  value.ngrdcol = static_cast<int>(value.hm_1.extent(0));
  value.nzt = static_cast<int>(value.hm_1.extent(1));
}

void roundtrip(clubb_kokkos::hydromet_pdf_parameter_collection& value) {
  auto input = to_fortran(value);
  auto output = make_output(value);
  f2k_roundtrip_hydromet_pdf_parameter_collection_c(input.data(), output.data());
  to_kokkos(value, output);
}

}  // namespace f2k
