// Compare ported temperature, skewness and PDF width kernels to compiled Fortran.
// Covers scalar/profile/column interfaces, zero and small variances, unequal and
// nearly equal gamma coefficients, positive/negative moments and the one-level
// PDF branch with/without momentum fluxes. Every floating-point check is exact.
#include "scratch_workspace_owner.hpp"
#include "CLUBB_core/T_in_K_module.hpp"
#include "CLUBB_core/Skx_module.hpp"
#include "CLUBB_core/sigma_sqd_w_module.hpp"
#include "CLUBB_core/parameter_indices.hpp"

#include <array>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using clubb_kokkos::real_t;
extern "C" {
void f2k_all_core_constants_c(real_t*, double*, int*, int*);
void f2k_temperature_reference_c(int, int, const real_t*, const real_t*, const real_t*, real_t*);
void f2k_skewness_reference_c(int, int, const real_t*, const real_t*, const real_t*, const real_t*,
                              const real_t*, const real_t*, const real_t*, real_t, real_t*);
void f2k_sigma_one_level_reference_c(int, const real_t*, const real_t*, const real_t*,
                                     const real_t*, const real_t*, const real_t*, const real_t*,
                                     const real_t*, const real_t*, const real_t*, const real_t*,
                                     int, real_t*);
}

template <class View, class Reference>
void compare(const char* name, View result, Reference reference) {
  // Transfer device output and identify the first differing column/level bit pattern.
  auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, result);
  if (host.extent(0) != reference.extent(0) || host.extent(1) != reference.extent(1))
    throw std::runtime_error("native oracle shape mismatch");

  for (std::size_t k = 0; k < host.extent(1); ++k)
    for (std::size_t i = 0; i < host.extent(0); ++i)
      if (std::memcmp(&host(i, k), &reference(i, k), sizeof(real_t))) {
        std::cerr << name << " column=" << i << " level=" << k << std::hexfloat
                  << " C++=" << host(i, k) << " Fortran=" << reference(i, k) << '\n';
        throw std::runtime_error("kernel differs from compiled Fortran");
      }
}

int main(int argc, char* argv[]) {
  // Exercise real native routines with column-varying inputs and preserve strict BFB.
  try {
    Kokkos::ScopeGuard kokkos(argc, argv);
    clubb_kokkos::test::ScratchWorkspaceOwner scratch_owner;
    auto& workspace = scratch_owner.workspace;
    using namespace clubb_kokkos;
    constexpr int nz = 17, nc = 4;
    using Host2D = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
    using Host3D = Kokkos::View<real_t***, Kokkos::LayoutLeft, Kokkos::HostSpace>;

    const std::array<real_t, 96> constants{parab_cyl_max_input,
                                           gamma_over_implicit_ts,
                                           pi,
                                           invrs_pi,
                                           radians_per_deg,
                                           sqrt_2pi,
                                           sqrt_2,
                                           one_hundred,
                                           fifty,
                                           thirty_six,
                                           twenty,
                                           eighteen,
                                           twelve,
                                           ten,
                                           nine,
                                           eight,
                                           seven,
                                           six,
                                           five,
                                           four,
                                           three,
                                           two,
                                           three_halves,
                                           four_thirds,
                                           one,
                                           three_fourths,
                                           two_thirds,
                                           one_half,
                                           one_third,
                                           one_fourth,
                                           three_sixteenth,
                                           one_eighth,
                                           zero,
                                           Cp,
                                           Lv,
                                           Ls,
                                           Lf,
                                           Rd,
                                           Rv,
                                           stefan_boltzmann,
                                           T_freeze_K,
                                           ep,
                                           ep1,
                                           ep2,
                                           kappa,
                                           omega_planet,
                                           grav,
                                           p0,
                                           vonk,
                                           rho_lw,
                                           rho_ice,
                                           w_tol,
                                           thl_tol,
                                           rt_tol,
                                           chi_tol,
                                           eta_tol,
                                           thl_tol_mfl,
                                           rt_tol_mfl,
                                           w_tol_sqd,
                                           rc_tol,
                                           Nc_tol,
                                           Ncn_tol,
                                           mvr_cloud_max,
                                           Nc_in_cloud_min,
                                           rr_tol,
                                           ri_tol,
                                           rs_tol,
                                           rg_tol,
                                           mvr_rain_max,
                                           mvr_ice_max,
                                           mvr_snow_max,
                                           mvr_graupel_max,
                                           Nr_tol,
                                           Ni_tol,
                                           Ns_tol,
                                           Ng_tol,
                                           em_min,
                                           eps,
                                           max_num_stdevs,
                                           zero_threshold,
                                           max_mag_correlation,
                                           max_mag_correlation_flux,
                                           min_max_smth_mag,
                                           cloud_frac_min,
                                           wp2_max,
                                           sec_per_day,
                                           sec_per_hr,
                                           sec_per_min,
                                           min_per_hr,
                                           g_per_kg,
                                           pascal_per_mb,
                                           cm3_per_m3,
                                           micron_per_m,
                                           cm_per_m,
                                           mm_per_m,
                                           unused_var};
    std::array<real_t, 96> native_constants;
    const std::array<double, 9> double_constants{
        pi_dp,  radians_per_deg_dp, sqrt_2pi_dp,   sqrt_2_dp, two_dp,
        one_dp, one_half_dp,        one_fourth_dp, zero_dp};
    std::array<double, 9> native_double_constants;
    const std::array<int, 5> integer_constants{fstderr, fstdin, fstdout, var_length,
                                               num_hf_draw_points};
    std::array<int, 5> native_integer_constants;
    const std::array<int, 106> all_indices{nparams,
                                           iC1 + 1,
                                           iC1b + 1,
                                           iC1c + 1,
                                           iC2rt + 1,
                                           iC2thl + 1,
                                           iC2rtthl + 1,
                                           iC4 + 1,
                                           iC_uu_shr + 1,
                                           iC_uu_buoy + 1,
                                           iC6rt + 1,
                                           iC6rtb + 1,
                                           iC6rtc + 1,
                                           iC6thl + 1,
                                           iC6thlb + 1,
                                           iC6thlc + 1,
                                           iC7 + 1,
                                           iC7b + 1,
                                           iC7c + 1,
                                           iC8 + 1,
                                           iC8b + 1,
                                           iC10 + 1,
                                           iC11 + 1,
                                           iC11b + 1,
                                           iC11c + 1,
                                           iC12 + 1,
                                           iC13 + 1,
                                           iC14 + 1,
                                           iC_wp2_pr_dfsn + 1,
                                           iC_wp3_pr_tp + 1,
                                           iC_wp3_pr_turb + 1,
                                           iC_wp3_pr_dfsn + 1,
                                           iC_wp2_splat + 1,
                                           iC6rt_Lscale0 + 1,
                                           iC6thl_Lscale0 + 1,
                                           iC7_Lscale0 + 1,
                                           iwpxp_L_thresh + 1,
                                           ic_K + 1,
                                           ic_K1 + 1,
                                           inu1 + 1,
                                           ic_K2 + 1,
                                           inu2 + 1,
                                           ic_K6 + 1,
                                           inu6 + 1,
                                           ic_K8 + 1,
                                           inu8 + 1,
                                           ic_K9 + 1,
                                           inu9 + 1,
                                           inu10 + 1,
                                           ic_K_hm + 1,
                                           ic_K_hmb + 1,
                                           iK_hm_min_coef + 1,
                                           inu_hm + 1,
                                           islope_coef_spread_DG_means_w + 1,
                                           ipdf_component_stdev_factor_w + 1,
                                           icoef_spread_DG_means_rt + 1,
                                           icoef_spread_DG_means_thl + 1,
                                           igamma_coef + 1,
                                           igamma_coefb + 1,
                                           igamma_coefc + 1,
                                           iLscale_mu_entrain + 1,
                                           ibeta + 1,
                                           ilmin_coef + 1,
                                           iomicron + 1,
                                           izeta_vrnce_rat + 1,
                                           iupsilon_precip_frac_rat + 1,
                                           ilambda0_stability_coef + 1,
                                           imult_coef + 1,
                                           itaumin + 1,
                                           itaumax + 1,
                                           iLscale_mu_coef + 1,
                                           iLscale_pert_coef + 1,
                                           ialpha_corr + 1,
                                           iSkw_denom_coef + 1,
                                           ic_K10 + 1,
                                           ic_K10h + 1,
                                           ithlp2_rad_coef + 1,
                                           ithlp2_rad_cloud_frac_thresh + 1,
                                           iup2_sfc_coef + 1,
                                           iSkw_max_mag + 1,
                                           iC_invrs_tau_bkgnd + 1,
                                           iC_invrs_tau_sfc + 1,
                                           iC_invrs_tau_shear + 1,
                                           iC_invrs_tau_N2 + 1,
                                           iC_invrs_tau_N2_wp2 + 1,
                                           iC_invrs_tau_N2_xp2 + 1,
                                           iC_invrs_tau_N2_wpxp + 1,
                                           iC_invrs_tau_N2_clear_wp3 + 1,
                                           iC_invrs_tau_wpxp_Ri + 1,
                                           iC_invrs_tau_wpxp_N2_thresh + 1,
                                           ixp3_coef_base + 1,
                                           ixp3_coef_slope + 1,
                                           ialtitude_threshold + 1,
                                           irtp2_clip_coef + 1,
                                           iCx_min + 1,
                                           iCx_max + 1,
                                           iRichardson_num_min + 1,
                                           iRichardson_num_max + 1,
                                           ia3_coef_min + 1,
                                           ia_const + 1,
                                           ibv_efold + 1,
                                           iwpxp_Ri_exp + 1,
                                           iz_displace + 1,
                                           iLscale_tke_multiplier + 1,
                                           iLscale_pdf_launch_frac + 1,
                                           iLscale_up_weight + 1};
    std::array<int, 106> native_all_indices;
    f2k_all_core_constants_c(native_constants.data(), native_double_constants.data(),
                             native_integer_constants.data(), native_all_indices.data());
    if (std::memcmp(constants.data(), native_constants.data(), sizeof(constants)))
      throw std::runtime_error("complete native constants mismatch");
    if (std::memcmp(double_constants.data(), native_double_constants.data(),
                    sizeof(double_constants)))
      throw std::runtime_error("complete native double_constants mismatch");
    if (std::memcmp(integer_constants.data(), native_integer_constants.data(),
                    sizeof(integer_constants)))
      throw std::runtime_error("complete native integer_constants mismatch");
    if (std::memcmp(all_indices.data(), native_all_indices.data(), sizeof(all_indices)))
      throw std::runtime_error("complete native all_indices mismatch");

    Host2D thlm("thlm", nc, nz), exner("exner", nc, nz), rcm("rcm", nc, nz);
    for (int k = 0; k < nz; ++k)
      for (int i = 0; i < nc; ++i) {
        thlm(i, k)  = real_t(280) + real_t(0.73) * real_t(k + i);
        exner(i, k) = real_t(0.57) + real_t(0.019) * real_t(k + i);
        rcm(i, k)   = (k % 3 == 0) ? zero : real_t(0.000013) * real_t(k + i + 1);
      }

    Host3D temperature("native_temperature", nc, nz, 4);
    f2k_temperature_reference_c(nz, nc, thlm.data(), exner.data(), rcm.data(), temperature.data());
    auto thlm_d  = to_execution_space(thlm);
    auto exner_d = to_execution_space(exner);
    auto rcm_d   = to_execution_space(rcm);
    Real2D scalar("scalar_temperature", nc, nz);
    Real2D scalar_inverse("scalar_inverse", nc, nz);
    Kokkos::parallel_for(
        "temperature::scalar", Kokkos::RangePolicy<>(0, nc * nz), KOKKOS_LAMBDA(const int cell) {
          const int i = cell % nc, k = cell / nc;
          scalar(i, k)         = thlm2T_in_K_api(thlm_d(i, k), exner_d(i, k), rcm_d(i, k));
          scalar_inverse(i, k) = T_in_K2thlm_api(scalar(i, k), exner_d(i, k), rcm_d(i, k));
        });
    compare("scalar temperature", scalar,
            Kokkos::subview(temperature, Kokkos::ALL, Kokkos::ALL, 0));
    compare("scalar inverse", scalar_inverse,
            Kokkos::subview(temperature, Kokkos::ALL, Kokkos::ALL, 3));

    auto columns = thlm2T_in_K_api(nz, nc, thlm_d, exner_d, rcm_d);
    compare("column temperature", columns,
            Kokkos::subview(temperature, Kokkos::ALL, Kokkos::ALL, 2));
    compare("column inverse", T_in_K2thlm_api(columns, exner_d, rcm_d),
            Kokkos::subview(temperature, Kokkos::ALL, Kokkos::ALL, 3));

    // Pack a vertical profile on device; a selected column is strided in LayoutLeft.
    Real1D thlm_1d("thlm_profile", nz), exner_1d("exner_profile", nz), rcm_1d("rcm_profile", nz);
    Kokkos::deep_copy(thlm_1d, Kokkos::subview(thlm_d, 2, Kokkos::ALL));
    Kokkos::deep_copy(exner_1d, Kokkos::subview(exner_d, 2, Kokkos::ALL));
    Kokkos::deep_copy(rcm_1d, Kokkos::subview(rcm_d, 2, Kokkos::ALL));
    auto profile   = thlm2T_in_K_api(nz, thlm_1d, exner_1d, rcm_1d);
    auto inverse   = T_in_K2thlm_api(profile, exner_1d, rcm_1d);
    auto profile_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, profile);
    auto inverse_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, inverse);
    for (int k = 0; k < nz; ++k)
      if (profile_h(k) != temperature(2, k, 1) || inverse_h(k) != temperature(2, k, 3))
        throw std::runtime_error("profile temperature or inverse differs");

    Host2D xp2("xp2", nc, nz), xp3("xp3", nc, nz), wp2("wp2", nc, nz), wpxp("wpxp", nc, nz);
    Host2D width("width", nc, nz), Skw("Skw", nc, nz), params("params", nc, nparams);
    for (int i = 0; i < nc; ++i) {
      params(i, iSkw_denom_coef) = real_t(0.25) + real_t(i);
      params(i, iSkw_max_mag)    = real_t(12);
      params(i, ibeta)           = real_t(0.2) * real_t(i + 1);
      params(i, igamma_coef)     = real_t(0.2) + real_t(0.1) * real_t(i);
      params(i, igamma_coefb)    = (i == 0 || i == 3) ? params(i, igamma_coef) : real_t(0.15);
      if (i == 3)
        params(i, igamma_coefb) += real_t(1.e-12);
      params(i, igamma_coefc) = (i == 0) ? zero : real_t(0.75) * real_t(i);
      for (int k = 0; k < nz; ++k) {
        xp2(i, k)   = (k % 4 == 0) ? zero : real_t(0.0013) * real_t(i + k + 1);
        xp3(i, k)   = real_t(0.0000017) * real_t((k % 2 ? -1 : 1) * (i + k));
        wp2(i, k)   = (k % 3 == 0) ? real_t(1.e-5) : real_t(0.019) * real_t(i + k + 1);
        wpxp(i, k)  = real_t(0.00013) * real_t((k % 2 ? -1 : 1) * (i + k));
        width(i, k) = real_t(0.025) * real_t(i + k);
        Skw(i, k)   = real_t(0.173) * real_t(k - 8 + i);
      }
    }
    auto xp2_d = to_execution_space(xp2), xp3_d = to_execution_space(xp3);
    auto wp2_d = to_execution_space(wp2), wpxp_d = to_execution_space(wpxp);
    auto width_d = to_execution_space(width), Skw_d = to_execution_space(Skw);
    auto params_d = to_execution_space(params);
    Host3D native_skew("native_skew", nc, nz, 4);
    Real2D result("skew_result", nc, nz);
    Real1D beta("beta", nc);
    Kokkos::deep_copy(beta, Kokkos::subview(params_d, Kokkos::ALL, ibeta));
    for (real_t x_tol : {real_t(0.02), real_t(0.01), real_t(1.e-8)}) {
      f2k_skewness_reference_c(nz, nc, xp2.data(), xp3.data(), wp2.data(), wpxp.data(),
                               width.data(), Skw.data(), params.data(), x_tol, native_skew.data());
      Skx_func(nz, nc, xp2_d, xp3_d, x_tol, params_d, result);
      compare("Skx_func", result, Kokkos::subview(native_skew, Kokkos::ALL, Kokkos::ALL, 0));
      compute_gamma_Skw(nz, nc, Skw_d, params_d, result);
      compare("compute_gamma_Skw", result,
              Kokkos::subview(native_skew, Kokkos::ALL, Kokkos::ALL, 1));
      LG_2005_ansatz(nz, nc, Skw_d, wpxp_d, wp2_d, xp2_d, beta, width_d, x_tol, result);
      compare("LG_2005_ansatz", result, Kokkos::subview(native_skew, Kokkos::ALL, Kokkos::ALL, 2));
      xp3_LG_2005_ansatz(workspace, nz, nc, Skw_d, wpxp_d, wp2_d, xp2_d, width_d, params_d, x_tol,
                         result);
      compare("xp3_LG_2005_ansatz", result,
              Kokkos::subview(native_skew, Kokkos::ALL, Kokkos::ALL, 3));
    }

    // A one-level grid avoids interpolation and isolates the full width formula.
    Host2D wp3_s("wp3", nc, 1), wp2_s("wp2", nc, 1), thlp2("thlp2", nc, 1), rtp2("rtp2", nc, 1);
    Host2D up2("up2", nc, 1), vp2("vp2", nc, 1), wpthlp("wpthlp", nc, 1), wprtp("wprtp", nc, 1);
    Host2D upwp("upwp", nc, 1), vpwp("vpwp", nc, 1), native_sigma("native_sigma", nc, 1);
    for (int i = 0; i < nc; ++i) {
      wp3_s(i, 0)  = real_t(0.007) * real_t(i - 1);
      wp2_s(i, 0)  = real_t(0.19) * real_t(i + 1);
      thlp2(i, 0)  = real_t(0.13) * real_t(i + 1);
      rtp2(i, 0)   = real_t(1.e-8) * real_t(i + 1);
      up2(i, 0)    = real_t(0.2);
      vp2(i, 0)    = real_t(0.3);
      wpthlp(i, 0) = real_t(0.01) * real_t(i);
      wprtp(i, 0)  = real_t(0.00001) * real_t(i);
      upwp(i, 0)   = (i == 3) ? real_t(2) : real_t(0.01) * real_t(i);
      vpwp(i, 0)   = (i == 2) ? real_t(0.2) : real_t(-0.01) * real_t(i);
    }
    Real2D sigma("sigma", nc, 1);
    grid gr;
    for (int predict : {0, 1}) {
      f2k_sigma_one_level_reference_c(nc, wp3_s.data(), wp2_s.data(), thlp2.data(), rtp2.data(),
                                      up2.data(), vp2.data(), wpthlp.data(), wprtp.data(),
                                      upwp.data(), vpwp.data(), params.data(), predict,
                                      native_sigma.data());
      compute_sigma_sqd_w(
          workspace, 1, 1, nc, gr, to_execution_space(wp3_s), to_execution_space(wp2_s),
          to_execution_space(thlp2), to_execution_space(rtp2), to_execution_space(up2),
          to_execution_space(vp2), to_execution_space(wpthlp), to_execution_space(wprtp),
          to_execution_space(upwp), to_execution_space(vpwp), params_d, predict != 0, sigma);
      compare("compute_sigma_sqd_w one level", sigma, native_sigma);
    }
    std::cout << "Temperature, skewness and one-level PDF kernels match compiled Fortran exactly\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
