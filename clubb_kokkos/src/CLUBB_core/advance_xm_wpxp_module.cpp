// Advance CLUBB means and turbulent fluxes with the native coupled equations.
//
// advance_xm_wpxp updates total water, liquid-water potential temperature,
// passive scalars, and optionally horizontal winds and their turbulent fluxes.
// Inputs/outputs use column-first execution-space Views. Private matrix arrays
// retain Native bands-first storage; the existing solver boundary transposes it.
// Sampled budgets use the C++ statistics context. Sponge damping follows JAX scope
// and is rejected by standalone initialization before this routine.
// References: Golaz et al. (2002), JAS 59, 3540–3551, Eqs.16–17;
// https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:wpxp_eqns
#include "scratch_workspace.hpp"
#include "fill_holes.hpp"
#include "pos_definite_module.hpp"
#include "mono_flux_limiter.hpp"
#include "turbulent_adv_pdf.hpp"
#include "mean_adv.hpp"
#include "diffusion.hpp"
#include "advance_xm_wpxp_module.hpp"
#include "constants_clubb.hpp"
#include "parameter_indices.hpp"
#include "model_flags.hpp"
#include "clip_explicit.hpp"
#include "Skx_module.hpp"
#include "advance_helper_module.hpp"
#include "matrix_solver_wrapper.hpp"
#include "error_code.hpp"
#include "stats_netcdf.hpp"
#include <Kokkos_MathematicalFunctions.hpp>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include <algorithm>
namespace clubb_kokkos {
namespace {
// Parameter constants
constexpr int nsub                     = 2; // Number of subdiagonals in the LHS matrix
constexpr int nsup                     = 2; // Number of superdiagonals in the LHS matrix
constexpr int xm_wpxp_thlm             = 1; // Named constant for thlm and wpthlp solving
constexpr int xm_wpxp_rtm              = 2; // Named constant for rtm and wprtp solving
constexpr int xm_wpxp_scalar           = 3; // Named constant for sclrm and wpsclrp solving
constexpr int xm_wpxp_um               = 4; // Named constant for optional um and upwp solving
constexpr int xm_wpxp_vm               = 5; // Named constant for optional vm and vpwp solving
constexpr int ndiags2                  = 2;
constexpr int ndiags3                  = 3;
[[maybe_unused]] constexpr int ndiags5 = 5;
template <class View> bool any_less(View values, real_t threshold) {
  // Express native ANY(values<threshold) across device cells for a host branch.
  int count = 0;
  Kokkos::parallel_reduce(
      "advance_xm_wpxp::any_less", values.extent(0) * values.extent(1),
      KOKKOS_LAMBDA(int cell, int& hits) {
        if (values(cell % values.extent(0), cell / values.extent(0)) < threshold)
          ++hits;
      },
      count);
  return count != 0;
}
template <class View> std::string dump_values(View values) {
  // Format Native diagnostics in first-index-fastest order. Host snapshots
  // already supplied by the caller can be read directly, including strided windows.
  const auto format_host = [](const auto& host) {
    using Host = std::decay_t<decltype(host)>;
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    for (std::size_t n = 0; n < host.size(); ++n) {
      if constexpr (Host::rank == 1)
        out << host(n) << ' ';
      else if constexpr (Host::rank == 2)
        out << host(n % host.extent(0), n / host.extent(0)) << ' ';
      else if constexpr (Host::rank == 3)
        out << host(n % host.extent(0), (n / host.extent(0)) % host.extent(1),
                    n / (host.extent(0) * host.extent(1)))
            << ' ';
      else {
        static_assert(Host::rank == 4);
        out << host(n % host.extent(0), (n / host.extent(0)) % host.extent(1),
                    (n / (host.extent(0) * host.extent(1))) % host.extent(2),
                    n / (host.extent(0) * host.extent(1) * host.extent(2)))
            << ' ';
      }
    }
    return out.str();
  };
  if constexpr (std::is_same_v<typename View::memory_space, Kokkos::HostSpace>)
    return format_host(values);
  else
    return format_host(to_host_space(values));
}
template <class View> std::string dump_header(View values, int column) {
  // Convert one fixed-width Native character diagnostic to a host string.
  auto host = to_host_space(values);
  std::string result;
  for (std::size_t j = 0; j < host.extent(0); ++j)
    result += host(j, column);
  return result;
}
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9>
void xm_wpxp_lhs(int nzm, int nzt, int ngrdcol, bool l_iter, real_t dt,                    // In
                 bool l_implemented, Array0 lhs_diff_zm, Array1 lhs_diff_zt,               // In
                 Array2 lhs_ma_zm, Array3 lhs_ma_zt, Array4 lhs_ta_wpxp, Array5 lhs_ta_xm, // In
                 Array6 lhs_tp, Array7 lhs_pr1, Array8 lhs_ac_pr2,                         // In
                 bool l_diffuse_rtm_and_thlm,                                              // In
                 Array9 lhs                                                                // Out
);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17>
void calc_xm_wpxp_lhs_terms(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol,
                            const grid& gr, Array0 wm_zm,
                            Array1 wm_zt, // In
                            Array2 wp2,   // In
                            Array3 Kh_zm, Array4 stability_correction, Array5 Kw6,
                            Array6 C7_Skw_fnc,                                                // In
                            Array7 invrs_rho_ds_zt, Array8 invrs_rho_ds_zm, Array9 rho_ds_zt, // In
                            Array10 rho_ds_zm, bool l_implemented,                            // In
                            const nu_vertical_res_dep& nu_vert_res_dep,                       // In
                            bool l_diffuse_rtm_and_thlm,                                      // In
                            bool l_stability_correct_Kh_N2_zm,                                // In
                            bool l_upwind_xm_ma,                                              // In
                            Array11 lhs_diff_zm, Array12 lhs_diff_zt, Array13 lhs_ma_zt,
                            Array14 lhs_ma_zm,                                    // Out
                            Array15 lhs_tp, Array16 lhs_ta_xm, Array17 lhs_ac_pr2 // Out
);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10>
void xm_wpxp_rhs(ScratchWorkspace& workspace, // InOut
                 int nzm, int nzt, int ngrdcol, const grid& gr, int solve_type, bool l_iter, // In
                 real_t dt,                                                                  // In
                 Array0 xm, Array1 wpxp, Array2 xm_forcing, Array3 wpxp_forcing,             // In
                 Array4 C7_Skw_fnc,                                                          // In
                 Array5 xpthvp, Array6 rhs_ta, Array7 thv_ds_zm,                             // In
                 Array8 lhs_pr1, Array9 lhs_ta_wpxp,                                         // In
                 stats_type& stats, // InOut
                 Array10 rhs        // Out
);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17>
void calc_xm_wpxp_ta_terms(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, int sclr_dim,
                           const grid& gr,
                           Array0 wp2rtp,                                                  // In
                           Array1 wp2thlp, Array2 wp2sclrp,                                // In
                           Array3 rho_ds_zt, Array4 invrs_rho_ds_zm, Array5 rho_ds_zm,     // In
                           Array6 sigma_sqd_w, Array7 wp3_on_wp2_zt,                       // In
                           const implicit_coefs_terms& pdf_implicit_coefs_terms,           // In
                           int iiPDF_type,                                                 // In
                           bool l_explicit_turbulent_adv_wpxp, bool l_predict_upwp_vpwp,   // In
                           bool l_scalar_calc,                                             // In
                           bool l_godunov_upwind_wpxp_ta,                                  // In
                           stats_type& stats,                                              // InOut
                           Array8 lhs_ta_wprtp, Array9 lhs_ta_wpthlp, Array10 lhs_ta_wpup, // Out
                           Array11 lhs_ta_wpvp, Array12 lhs_ta_wpsclrp,                    // Out
                           Array13 rhs_ta_wprtp, Array14 rhs_ta_wpthlp, Array15 rhs_ta_wpup, // Out
                           Array16 rhs_ta_wpvp, Array17 rhs_ta_wpsclrp                       // Out
);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26, class Array27, class Array28, class Array29,
          class Array30, class Array31, class Array32, class Array33, class Array34, class Array35,
          class Array36, class Array37, class Array38, class Array39, class Array40, class Array41,
          class Array42, class Array43, class Array44, class Array45, class Array46, class Array47,
          class Array48, class Array49, class Array50, class Array51, class Array52, class Array53,
          class Array54, class Array55, class Array56, class Array57, class Array58, class Array59,
          class Array60, class Array61, class Array62, class Array63, class Array64, class Array65,
          class Array66, class Array67>
void solve_xm_wpxp_with_single_lhs(
    ScratchWorkspace& workspace,                                                    // InOut
    int nzm, int nzt, int ngrdcol, int sclr_dim, Array0 sclr_tol, const grid& gr,   // In
    real_t dt,                                                                      // In
    bool l_iter, int nrhs,                                                          // In
    Array1 wm_zt, Array2 wp2, Array3 invrs_tau_C6_zm, Array4 tau_max_zm,            // In
    Array5 rtpthvp, Array6 rtm_forcing, Array7 wprtp_forcing, Array8 thlpthvp,      // In
    Array9 thlm_forcing, Array10 wpthlp_forcing, Array11 rho_ds_zm,                 // In
    Array12 rho_ds_zt, Array13 invrs_rho_ds_zm, Array14 invrs_rho_ds_zt,            // In
    Array15 thv_ds_zm, Array16 rtp2, Array17 thlp2, bool l_implemented,             // In
    Array18 sclrpthvp, Array19 sclrm_forcing, Array20 sclrp2, Array21 um_forcing,   // In
    Array22 vm_forcing, Array23 ug, Array24 vg, Array25 uprcp, Array26 vprcp,       // In
    Array27 rc_coef_zm, Array28 fcor,                                               // In
    Array29 fcor_y, Array30 up2, Array31 vp2,                                       // In
    Array32 low_lev_effect, Array33 high_lev_effect,                                // In
    Array34 C6rt_Skw_fnc, Array35 C6thl_Skw_fnc, Array36 C7_Skw_fnc,                // In
    Array37 lhs_diff_zm, Array38 lhs_diff_zt, Array39 lhs_ma_zt, Array40 lhs_ma_zm, // In
    Array41 lhs_ta_wpxp,                                                            // In
    Array42 rhs_ta_wprtp, Array43 rhs_ta_wpthlp, Array44 rhs_ta_wpup,               // In
    Array45 rhs_ta_wpvp, Array46 rhs_ta_wpsclrp,                                    // In
    Array47 lhs_tp, Array48 lhs_ta_xm, Array49 lhs_ac_pr2, Array50 lhs_pr1_wprtp,   // In
    Array51 lhs_pr1_wpthlp, Array52 lhs_pr1_wpsclrp,                                // In
    Array53 C_uu_shr,                                                               // In
    int penta_solve_method,                                                         // In
    int tridiag_solve_method,                                                       // In
    int fill_holes_type,                                                            // In
    bool l_predict_upwp_vpwp,                                                       // In
    bool l_ho_nontrad_coriolis,                                                     // In
    bool l_ho_trad_coriolis,                                                        // In
    bool l_diffuse_rtm_and_thlm,                                                    // In
    bool l_upwind_xm_ma,                                                            // In
    bool l_enable_relaxed_clipping,                                                 // In
    bool l_perturbed_wind,                                                          // In
    bool l_mono_flux_lim_thlm,                                                      // In
    bool l_mono_flux_lim_rtm,                                                       // In
    bool l_mono_flux_lim_um,                                                        // In
    bool l_mono_flux_lim_vm,                                                        // In
    bool l_mono_flux_lim_spikefix,                                                  // In
    int& wprtp_cl_num,                                                              // InOut
    int& wpthlp_cl_num, int& upwp_cl_num, int& vpwp_cl_num,                         // InOut
    stats_type& stats,                                                              // InOut
    Array54 rtm, Array55 wprtp, Array56 thlm, Array57 wpthlp,                       // InOut
    Array58 sclrm, Array59 wpsclrp, Array60 um, Array61 upwp, Array62 vm,           // InOut
    Array63 vpwp,                                                                   // InOut
    Array64 um_pert, Array65 vm_pert, Array66 upwp_pert, Array67 vpwp_pert,         // InOut
    err_info_type& err_info                                                         // InOut
);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26, class Array27, class Array28, class Array29,
          class Array30, class Array31, class Array32, class Array33, class Array34, class Array35,
          class Array36, class Array37, class Array38, class Array39, class Array40, class Array41,
          class Array42, class Array43>
void solve_xm_wpxp_with_multiple_lhs(
    ScratchWorkspace& workspace,                                                    // InOut
    int nzm, int nzt, int ngrdcol, int sclr_dim, Array0 sclr_tol, const grid& gr,   // In
    real_t dt,                                                                      // In
    bool l_iter, int nrhs, Array1 wm_zt, Array2 wp2,                                // In
    Array3 rtpthvp, Array4 rtm_forcing, Array5 wprtp_forcing, Array6 thlpthvp,      // In
    Array7 thlm_forcing, Array8 wpthlp_forcing, Array9 rho_ds_zm,                   // In
    Array10 rho_ds_zt, Array11 invrs_rho_ds_zm, Array12 invrs_rho_ds_zt,            // In
    Array13 thv_ds_zm, Array14 rtp2, Array15 thlp2, bool l_implemented,             // In
    Array16 sclrpthvp, Array17 sclrm_forcing, Array18 sclrp2,                       // In
    Array19 low_lev_effect, Array20 high_lev_effect, Array21 C7_Skw_fnc,            // In
    Array22 lhs_diff_zm, Array23 lhs_diff_zt, Array24 lhs_ma_zt, Array25 lhs_ma_zm, // In
    Array26 lhs_ta_wprtp, Array27 lhs_ta_wpthlp, Array28 lhs_ta_wpsclrp,            // In
    Array29 rhs_ta_wprtp, Array30 rhs_ta_wpthlp, Array31 rhs_ta_wpsclrp,            // In
    Array32 lhs_tp, Array33 lhs_ta_xm, Array34 lhs_ac_pr2, Array35 lhs_pr1_wprtp,   // In
    Array36 lhs_pr1_wpthlp, Array37 lhs_pr1_wpsclrp,                                // In
    int penta_solve_method,                                                         // In
    int tridiag_solve_method,                                                       // In
    int fill_holes_type,                                                            // In
    bool l_diffuse_rtm_and_thlm,                                                    // In
    bool l_upwind_xm_ma,                                                            // In
    bool l_enable_relaxed_clipping,                                                 // In
    bool l_mono_flux_lim_thlm,                                                      // In
    bool l_mono_flux_lim_rtm,                                                       // In
    bool l_mono_flux_lim_um,                                                        // In
    bool l_mono_flux_lim_vm,                                                        // In
    bool l_mono_flux_lim_spikefix,                                                  // In
    int& wprtp_cl_num, int& wpthlp_cl_num,                                          // InOut
    stats_type& stats,                                                              // InOut
    Array38 rtm, Array39 wprtp, Array40 thlm, Array41 wpthlp, Array42 sclrm,        // InOut
    Array43 wpsclrp, err_info_type& err_info                                        // InOut
);
template <class Array0, class Array1, class Array2, class Array3>
void xm_wpxp_solve(ScratchWorkspace& workspace,                     // InOut
                   int nzm, int ngrdcol, int nrhs, const grid& gr,  // In
                   int penta_solve_method,                          // In
                   bool l_implemented,                              // In
                   Array0 old_solution,                             // InOut
                   Array1 lhs, Array2 rhs, err_info_type& err_info, // InOut
                   Array3 solution,                                 // Out
                   Real1D* rcond = nullptr                          // Optional Out
);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21>
void xm_wpxp_clipping_and_stats(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol,
                                const grid& gr, int solve_type,
                                real_t dt,                                                 // In
                                Array0 wp2, Array1 xp2, Array2 wm_zt,                      // In
                                Array3 xm_forcing, Array4 rho_ds_zm, Array5 rho_ds_zt,     // In
                                Array6 invrs_rho_ds_zm, Array7 invrs_rho_ds_zt,            // In
                                real_t xp2_threshold, real_t xm_threshold, Array8 rcond,   // In
                                Array9 low_lev_effect, Array10 high_lev_effect,            // In
                                Array11 lhs_ma_zt, Array12 lhs_ma_zm, Array13 lhs_ta_wpxp, // In
                                Array14 lhs_diff_zm, Array15 C7_Skw_fnc,                   // In
                                Array16 lhs_tp, Array17 lhs_ta_xm, Array18 lhs_pr1,        // In
                                bool l_implemented, Array19 solution,                      // In
                                int tridiag_solve_method,                                  // In
                                int fill_holes_type,                                       // In
                                bool l_upwind_xm_ma,                                       // In
                                bool l_enable_relaxed_clipping,                            // In
                                bool l_mono_flux_lim_thlm,                                 // In
                                bool l_mono_flux_lim_rtm,                                  // In
                                bool l_mono_flux_lim_um,                                   // In
                                bool l_mono_flux_lim_vm,                                   // In
                                bool l_mono_flux_lim_spikefix,                             // In
                                stats_type& stats,                                         // InOut
                                Array20 xm,                                                // InOut
                                real_t xm_tol,                                             // In
                                Array21 wpxp, err_info_type& err_info,                     // InOut
                                int* wpxp_cl_num = nullptr // Optional InOut
);
template <class Array0, class Array1, class Array2>
void xm_term_ta_lhs(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                    Array0 rho_ds_zm, Array1 invrs_rho_ds_zt,      // In
                    Array2 lhs_ta_xm                               // Out
);
template <class Array0, class Array1>
void wpxp_term_tp_lhs(int nzm, int ngrdcol, const grid& gr, Array0 wp2, // In
                      Array1 lhs_tp                                     // Out
);
template <class Array0, class Array1, class Array2, class Array3>
void wpxp_terms_ac_pr2_lhs(int nzm, int nzt, int ngrdcol, const grid& gr, Array0 C7_Skw_fnc, // In
                           Array1 wm_zt, Array2 invrs_dzm,                                   // In
                           Array3 lhs_ac_pr2                                                 // Out
);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6>
void wpxp_term_pr1_lhs(int nzm, int ngrdcol, const grid& gr, Array0 C6rt_Skw_fnc, // In
                       Array1 C6thl_Skw_fnc, Array2 C7_Skw_fnc,                   // In
                       Array3 invrs_tau_C6_zm, bool l_scalar_calc,                // In
                       Array4 lhs_pr1_wprtp, Array5 lhs_pr1_wpthlp,               // Out
                       Array6 lhs_pr1_wpsclrp                                     // Out
);
template <class Array0, class Array1, class Array2, class Array3>
void wpxp_terms_bp_pr3_rhs(int nzm, int ngrdcol, const grid& gr, Array0 C7_Skw_fnc,
                           Array1 thv_ds_zm, // In
                           Array2 xpthvp,    // In
                           Array3 rhs_bp_pr3 // Out
);
template <class Array0, class Array1, class Array2>
void xm_correction_wpxp_cl(ScratchWorkspace& workspace,                              // InOut
                           int nzm, int nzt, int ngrdcol, int solve_type, real_t dt, // In
                           Array0 wpxp_chnge, Array1 invrs_dzt,                      // In
                           stats_type& stats,                                        // InOut
                           Array2 xm                                                 // InOut
);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6>
void damp_coefficient(int nzm, int ngrdcol, const grid& gr, Array0 coefficient,
                      Array1 Cx_Skw_fnc,                                 // In
                      Array2 max_coeff_value, Array3 altitude_threshold, // In
                      Array4 threshold, Array5 Lscale_zm,                // In
                      Array6 damped_value                                // Out
);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7>
void diagnose_upxp(ScratchWorkspace& workspace,                                           // InOut
                   int nzm, int nzt, int ngrdcol, const grid& gr, Array0 ypwp, Array1 xm, // In
                   Array2 wpxp, Array3 ym,                                                // In
                   Array4 C6x_Skw_fnc, Array5 tau_C6_zm, Array6 C7_Skw_fnc,               // In
                   Array7 ypxp                                                            // Out
);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26, class Array27, class Array28, class Array29,
          class Array30, class Array31, class Array32, class Array33, class Array34, class Array35,
          class Array36, class Array37, class Array38, class Array39, class Array40, class Array41,
          class Array42, class Array43, class Array44, class Array45, class Array46, class Array47,
          class Array48, class Array49, class Array50, class Array51, class Array52, class Array53,
          class Array54, class Array55, class Array56, class Array57, class Array58, class Array59,
          class Array60, class Array61, class Array62, class Array63, class Array64, class Array65,
          class Array66, class Array67, class Array68, class Array69, class Array70, class Array71,
          class Array72>
void error_prints_xm_wpxp(int nzm, int nzt, int sclr_dim, Array0 zm, Array1 zt,      // In
                          real_t dt, Array2 sigma_sqd_w, Array3 wm_zm, Array4 wm_zt, // In
                          Array5 wp2,                                                // In
                          Array6 Lscale_zm, Array7 wp3_on_wp2,                       // In
                          Array8 wp3_on_wp2_zt,                                      // In
                          Array9 Kh_zt, Array10 Kh_zm, Array11 invrs_tau_C6_zm,      // In
                          Array12 Skw_zm,                                            // In
                          Array13 wp2rtp, Array14 rtpthvp, Array15 rtm_forcing,      // In
                          Array16 wprtp_forcing, Array17 rtm_ref, Array18 wp2thlp,   // In
                          Array19 thlpthvp, Array20 thlm_forcing,                    // In
                          Array21 wpthlp_forcing, Array22 thlm_ref,                  // In
                          Array23 rho_ds_zm,                                         // In
                          Array24 rho_ds_zt, Array25 invrs_rho_ds_zm,                // In
                          Array26 invrs_rho_ds_zt, Array27 thv_ds_zm, Array28 rtp2,  // In
                          Array29 thlp2, Array30 w_1_zm, Array31 w_2_zm,             // In
                          Array32 varnce_w_1_zm, Array33 varnce_w_2_zm,              // In
                          Array34 mixt_frac_zm, bool l_implemented, Array35 em,      // In
                          Array36 wp2sclrp, Array37 sclrpthvp,                       // In
                          Array38 sclrm_forcing,                                     // In
                          Array39 sclrp2, Array40 Cx_fnc_Richardson,                 // In
                          const implicit_coefs_terms& pdf_implicit_coefs_terms,      // In
                          Array41 um_forcing, Array42 vm_forcing,                    // In
                          Array43 ug, Array44 vg,                                    // In
                          Array45 wpthvp, real_t fcor, Array46 um_ref,               // In
                          Array47 vm_ref, Array48 up2,                               // In
                          Array49 vp2, Array50 uprcp, Array51 vprcp,                 // In
                          Array52 rc_coef_zm, Array53 rtm,                           // In
                          Array54 wprtp, Array55 thlm, Array56 wpthlp,               // In
                          Array57 sclrm, Array58 wpsclrp,                            // In
                          Array59 um, Array60 upwp, Array61 vm,                      // In
                          Array62 vpwp, Array63 rtm_old,                             // In
                          Array64 wprtp_old, Array65 thlm_old,                       // In
                          Array66 wpthlp_old,                                        // In
                          Array67 sclrm_old, Array68 wpsclrp_old,                    // In
                          Array69 um_old,                                            // In
                          Array70 upwp_old, Array71 vm_old, Array72 vpwp_old,        // In
                          bool l_predict_upwp_vpwp, bool l_lmm_stepping              // In
);
} // namespace

// nzm: Number of momentum vertical levels
// nzt: Number of thermodynamic vertical levels
// ngrdcol: Number of grid columns
// sclr_dim: Number of passive scalars
// sclr_tol: Threshold(s) on the passive scalars  [units vary]
// dt: Timestep                                 [s]
// sigma_sqd_w: sigma_sqd_w on momentum levels           [-]
// wm_zm: w wind component on momentum levels      [m/s]
// wm_zt: w wind component on thermodynamic levels [m/s]
// wp2: w'^2 (momentum levels)                   [m^2/s^2]
// Lscale_zm: Turbulent mixing length interp to m levs [m]
// wp3: <w'^3> on thermodynamic levels           [m^3/s^3]
// Kh_zt: Eddy diffusivity on thermodynamic levels [m^2/s]
// Kh_zm: Eddy diffusivity on momentum levels
// stability_correction: Stability correction factor
// invrs_tau_C6_zm: Inverse time-scale on mom. levels applied to C6 term [1/s]
// tau_max_zm: Max. allowable eddy dissipation time scale on m-levs  [s]
// wp2rtp: <w'^2 r_t'> (thermodynamic levels)    [m^2/s^2 kg/kg]
// rtpthvp: r_t'th_v' (momentum levels)              [(kg/kg) K]
// rtm_forcing: r_t forcing (thermodynamic levels)       [(kg/kg)/s]
// wprtp_forcing: <w'r_t'> forcing (momentum levels)       [(kg/kg)/s^2]
// rtm_ref: rtm for nudging                          [kg/kg]
// wp2thlp: <w'^2 th_l'> (thermodynamic levels)      [m^2/s^2 K]
// thlpthvp: th_l'th_v' (momentum levels)             [K^2]
// thlm_forcing: th_l forcing (thermodynamic levels)      [K/s]
// wpthlp_forcing: <w'th_l'> forcing (momentum levels)      [K/s^2]
// thlm_ref: thlm for nudging                         [K]
// rho_ds_zm: Dry, static density on momentum levels   [kg/m^3]
// rho_ds_zt: Dry, static density on thermo. levels    [kg/m^3]
// invrs_rho_ds_zm: Inv. dry, static density @ moment. levs. [m^3/kg]
// invrs_rho_ds_zt: Inv. dry, static density @ thermo. levs. [m^3/kg]
// thv_ds_zm: Dry, base-state theta_v on moment. levs. [K]
// rtp2: r_t'^2 (momentum levels)                 [(kg/kg)^2]
// thlp2: th_l'^2 (momentum levels)                [K^2]
// w_1_zm: Mean w (1st PDF component)              [m/s]
// w_2_zm: Mean w (2nd PDF component)              [m/s]
// varnce_w_1_zm: Variance of w (1st PDF component)       [m^2/s^2]
// varnce_w_2_zm: Variance of w (2nd PDF component)       [m^2/s^2]
// mixt_frac_zm: Weight of 1st PDF component (Sk_w dependent) [-]
// l_implemented: True if CLUBB is being implemented and run in a host model
// em: Turbulent Kinetic Energy (TKE)           [m^2/s^2]
// wp2sclrp: <w'^2 sclr'> (thermodynamic levels)   [Units vary]
// sclrpthvp: <sclr' th_v'> (momentum levels)       [Units vary]
// sclrm_forcing: sclrm forcing (thermodynamic levels)  [Units vary]
// sclrp2: For clipping Vince Larson             [Units vary]
// Cx_fnc_Richardson: Cx_fnc computed from Richardson_num       [-]
// pdf_implicit_coefs_terms: Implicit coefs / explicit terms [units vary]
// um_forcing: <u> forcing term (thermodynamic levels)      [m/s^2]
// vm_forcing: <v> forcing term (thermodynamic levels)      [m/s^2]
// ug: <u> geostrophic wind (thermodynamic levels)  [m/s]
// vg: <v> geostrophic wind (thermodynamic levels)  [m/s]
// wpthvp: <w'thv'> (momentum levels)                   [m/s K]
// fcor: Traditional Coriolis parameter               [s^-1] Vertical planetary vorticity.
// Proportional to sin(latitude)
// fcor_y: Nontraditional Coriolis parameter            [s^-1] Meridional planetary vorticity.
// Proportional to cos(latitude)
// um_ref: Reference u wind component for nudging       [m/s]
// vm_ref: Reference v wind component for nudging       [m/s]
// up2: Variance of the u wind component             [m^2/s^2]
// vp2: Variance of the v wind component             [m^2/s^2]
// uprcp: < u' r_c' >                               [(m kg)/(s kg)]
// vprcp: < v' r_c' >                               [(m kg)/(s kg)]
// rc_coef_zm: Coefficient on X'r_c' in X'th_v' equation [K/(kg/kg)]
// clubb_params: Array of CLUBB's tunable parameters    [units vary]
// nu_vert_res_dep: Vertical resolution dependent nu values
// ts_nudge: Timescale of u/v nudging             [s]
// iiPDF_type: Selected option for the two-component normal (double Gaussian) PDF type to use
// for the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate,
// two-component PDF.
// penta_solve_method: Method to solve then penta-diagonal system
// tridiag_solve_method: Specifier for method to solve tridiagonal systems
// fill_holes_type: Specifier for which hole filling method to use
// l_predict_upwp_vpwp: Flag to predict <u'w'> and <v'w'> along with <u> and <v> alongside the
// advancement of <rt>, <w'rt'>, <thl>, <wpthlp>, <sclr>, and <w'sclr'> in subroutine
// advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still approximated by eddy diffusivity
// when <u> and <v> are advanced in subroutine advance_windm_edsclrm.
// l_ho_nontrad_coriolis: Flag to implement the nontraditional Coriolis terms in the prognostic
// equations of <w'w'>, <u'w'>, and <u'u'>.
// l_ho_trad_coriolis: Flag to implement the traditional Coriolis terms in the prognostic
// equations of <v'w'> and <u'w'>.
// l_diffuse_rtm_and_thlm: This flag determines whether or not we want CLUBB to do diffusion on
// rtm and thlm
// l_stability_correct_Kh_N2_zm: This flag determines whether or not we want CLUBB to apply a
// stability correction
// l_godunov_upwind_wpxp_ta: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent advection terms. It affects
// wpxp only.
// l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent or mean advection terms. It
// affects rtm, thlm, sclrm, um and vm.
// l_uv_nudge: For wind speed nudging
// l_diag_Lscale_from_tau: First diagnose dissipation time tau, and then diagnose the mixing
// length scale as Lscale = tau * tke
// l_use_C7_Richardson: Parameterize C7 based on Richardson number
// l_lmm_stepping: Apply Linear Multistep Method (LMM) Stepping
// l_enable_relaxed_clipping: Flag to relax clipping on wpxp in xm_wpxp_clipping_and_stats
// l_linearize_pbl_winds: Flag (used by E3SM) to linearize PBL winds
// l_mono_flux_lim_thlm: Flag to turn on monotonic flux limiter for thlm
// l_mono_flux_lim_rtm: Flag to turn on monotonic flux limiter for rtm
// l_mono_flux_lim_um: Flag to turn on monotonic flux limiter for um
// l_mono_flux_lim_vm: Flag to turn on monotonic flux limiter for vm
// l_mono_flux_lim_spikefix: Flag to implement monotonic flux limiter code that eliminates
// spurious drying tendencies at model top
// rtm: r_t  (total water mixing ratio)           [kg/kg]
// wprtp: w'r_t'                                    [(kg/kg) m/s]
// thlm: th_l (liquid water potential temperature) [K]
// wpthlp: w'th_l'                                   [K m/s]
// sclrm: [Units vary]
// wpsclrp: [Units vary]
// um: <u>:  mean west-east horiz. velocity (thermo. levs.)   [m/s]
// upwp: <u'w'>:  momentum flux (momentum levels)               [m^2/s^2]
// vm: <v>:  mean south-north horiz. velocity (thermo. levs.) [m/s]
// vpwp: <v'w'>:  momentum flux (momentum levels)               [m^2/s^2]
// um_pert: perturbed <u>    [m/s]
// vm_pert: perturbed <v>    [m/s]
// upwp_pert: perturbed <u'w'> [m^2/s^2]
// vpwp_pert: perturbed <v'w'> [m^2/s^2]
// err_info: err_info struct containing err_code and err_header
void advance_xm_wpxp(
    ScratchWorkspace& workspace,                                                  // InOut
    int nzm, int nzt, int ngrdcol, int sclr_dim, Real1D sclr_tol, const grid& gr, // In
    real_t dt,                                                                    // In
    Real2D sigma_sqd_w, Real2D wm_zm, Real2D wm_zt, Real2D wp2, Real2D Lscale_zm, // In
    Real2D wp3, Real2D Kh_zt, Real2D Kh_zm,                                       // In
    Real2D stability_correction,                                                  // In
    Real2D invrs_tau_C6_zm, Real2D tau_max_zm, Real2D wp2rtp, Real2D rtpthvp,     // In
    Real2D rtm_forcing, Real2D wprtp_forcing, Real2D rtm_ref, Real2D wp2thlp,     // In
    Real2D thlpthvp, Real2D thlm_forcing, Real2D wpthlp_forcing, Real2D thlm_ref, // In
    Real2D rho_ds_zm, Real2D rho_ds_zt, Real2D invrs_rho_ds_zm,                   // In
    Real2D invrs_rho_ds_zt, Real2D thv_ds_zm, Real2D rtp2, Real2D thlp2,          // In
    Real2D w_1_zm, Real2D w_2_zm, Real2D varnce_w_1_zm, Real2D varnce_w_2_zm,     // In
    Real2D mixt_frac_zm, bool l_implemented, Real2D em, Real3D wp2sclrp,          // In
    Real3D sclrpthvp, Real3D sclrm_forcing, Real3D sclrp2,                        // In
    Real2D Cx_fnc_Richardson,                                                     // In
    const implicit_coefs_terms& pdf_implicit_coefs_terms,                         // In
    Real2D um_forcing, Real2D vm_forcing, Real2D ug, Real2D vg, Real2D wpthvp,    // In
    Real1D fcor, Real1D fcor_y, Real2D um_ref, Real2D vm_ref, Real2D up2,         // In
    Real2D vp2,                                                                   // In
    Real2D uprcp, Real2D vprcp, Real2D rc_coef_zm,                                // In
    Real2D clubb_params, const nu_vertical_res_dep& nu_vert_res_dep,              // In
    real_t ts_nudge,                                                              // In
    int iiPDF_type,                                                               // In
    int penta_solve_method,                                                       // In
    int tridiag_solve_method,                                                     // In
    int fill_holes_type,                                                          // In
    bool l_predict_upwp_vpwp,                                                     // In
    bool l_ho_nontrad_coriolis,                                                   // In
    bool l_ho_trad_coriolis,                                                      // In
    bool l_diffuse_rtm_and_thlm,                                                  // In
    bool l_stability_correct_Kh_N2_zm,                                            // In
    bool l_godunov_upwind_wpxp_ta,                                                // In
    bool l_upwind_xm_ma,                                                          // In
    bool l_uv_nudge,                                                              // In
    bool l_diag_Lscale_from_tau,                                                  // In
    bool l_use_C7_Richardson,                                                     // In
    bool l_lmm_stepping,                                                          // In
    bool l_enable_relaxed_clipping,                                               // In
    bool l_linearize_pbl_winds,                                                   // In
    bool l_mono_flux_lim_thlm,                                                    // In
    bool l_mono_flux_lim_rtm,                                                     // In
    bool l_mono_flux_lim_um,                                                      // In
    bool l_mono_flux_lim_vm,                                                      // In
    bool l_mono_flux_lim_spikefix,                                                // In
    int& wprtp_cl_num, int& wpthlp_cl_num, int& upwp_cl_num, int& vpwp_cl_num,    // InOut
    stats_type& stats,                                                            // InOut
    Real2D rtm, Real2D wprtp, Real2D thlm, Real2D wpthlp,                         // InOut
    Real3D sclrm, Real3D wpsclrp, Real2D um, Real2D upwp, Real2D vm, Real2D vpwp, // InOut
    Real2D um_pert, Real2D vm_pert, Real2D upwp_pert, Real2D vpwp_pert,           // InOut
    err_info_type& err_info                                                       // InOut
) {
  // Description:
  // Advance the mean and flux terms by one timestep.
  // References:
  // https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:wpxp_eqns
  //
  // Eqn. 16 & 17 on p. 3546 of
  // ``A PDF-Based Model for Boundary Layer Clouds. Part I:
  //   Method and Model Description'' Golaz, et al. (2002)
  //   JAS, Vol. 59, pp. 3540--3551.
  // See Also
  // ``Equations for CLUBB'' Section 5:
  //   /Implicit solutions for the means and fluxes/
  //-----------------------------------------------------------------------
  // -------------------- Input Variables --------------------
  // Added for clipping by Vince Larson 29 Sep 2007
  // End of Vince Larson's addition.
  // Additional variables for passive scalars
  // Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
  // Vertical planetary vorticity.   Proportional to sin(latitude)
  // Meridional planetary vorticity. Proportional to cos(latitude)
  // Gaussian) PDF type to use for the w, rt, and theta-l (or
  // w, chi, and eta) portion of CLUBB's multivariate,
  // two-component PDF.
  // alongside the advancement of <rt>, <w'rt'>, <thl>,
  // <wpthlp>, <sclr>, and <w'sclr'> in subroutine
  // advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still
  // approximated by eddy diffusivity when <u> and <v> are
  // advanced in subroutine advance_windm_edsclrm.
  // prognostic equations of <w'w'>, <u'w'>, and <u'u'>.
  // prognostic equations of <v'w'> and <u'w'>.
  // diffusion on rtm and thlm
  // a stability correction
  // differencing approximation rather than a centered
  // differencing for turbulent advection terms.
  // It affects  wpxp only.
  // differencing approximation rather than a centered
  // differencing for turbulent or mean advection terms.
  // It affects rtm, thlm, sclrm, um and vm.
  // mixing length scale as Lscale = tau * tke
  // eliminates spurious drying tendencies at model top
  //------------------- InOut Variables -------------------
  // Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
  // Variables used to track perturbed version of winds.
  // -------------------- Local Variables --------------------
  // Parameter Constants
  // Eddy Diffusion for wpthlp and wprtp.
  // Variables used as part of the monotonic turbulent advection scheme.
  // Find the lowermost and uppermost grid levels that can have an effect
  // on the central thermodynamic level during the course of a time step,
  // due to the effects of turbulent advection only.
  // Constant parameters as a function of Skw.
  // Saved values of predictive fields, prior to being advanced, for use in
  // print statements in case of fatal error.
  // -------------------- Input/Output Variables --------------------
  // Input/Output Variables
  // Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
  // LHS/RHS terms
  // Whether preturbed winds are being solved.
  auto workspace_scope  = workspace.scope();
  constexpr bool l_iter = true;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> C6rt_Skw_fnc =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::C6rt_Skw_fnc", ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> C6thl_Skw_fnc =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::C6thl_Skw_fnc", ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> C7_Skw_fnc =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::C7_Skw_fnc", ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> C6_term =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::C6_term", ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wp3_zm =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::wp3_zm", ngrdcol,
          nzm); // w'^3 on momentum levels [m^3/s^3]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> Skw_zm =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::Skw_zm", ngrdcol,
          nzm); // Skewness of w on momentum levels [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wp3_on_wp2 =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::wp3_on_wp2", ngrdcol,
          nzm); // Smoothed wp3 / wp2 on momentum levels [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> Kw6 =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::Kw6", ngrdcol,
          nzt); // wpxp eddy diff. [m^2/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wp3_on_wp2_zt =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::wp3_on_wp2_zt", ngrdcol,
          nzt); // Smoothed wp3 / wp2 on thermo. levels [m/s]
  Kokkos::View<int**, Kokkos::LayoutLeft> low_lev_effect =
      workspace_scope.view<Kokkos::View<int**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::low_lev_effect", ngrdcol,
          nzt); // Index of the lowest level that has an effect.
  Kokkos::View<int**, Kokkos::LayoutLeft> high_lev_effect =
      workspace_scope.view<Kokkos::View<int**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::high_lev_effect", ngrdcol,
          nzt); // Index of the highest level that has an effect.
  int nrhs;     // Number of RHS vectors
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rtm_old =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::rtm_old", ngrdcol,
          nzt); // Saved value of r_t        [kg/kg]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> thlm_old =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::thlm_old", ngrdcol,
          nzt); // Saved value of th_l       [K]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wprtp_old =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::wprtp_old", ngrdcol,
          nzm); // Saved value of w'r_t'     [(kg/kg) m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wpthlp_old =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::wpthlp_old", ngrdcol,
          nzm); // Saved value of w'th_l'    [K m/s]
  Kokkos::View<real_t***, Kokkos::LayoutLeft> sclrm_old =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::sclrm_old", ngrdcol, nzt,
          sclr_dim); // Saved value of sclr      [units vary]
  Kokkos::View<real_t***, Kokkos::LayoutLeft> wpsclrp_old =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::wpsclrp_old", ngrdcol, nzm,
          sclr_dim); // Saved value of wpsclrp   [units vary]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> um_old =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::um_old", ngrdcol,
          nzt); // Saved value of <u>       [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vm_old =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::vm_old", ngrdcol,
          nzt); // Saved value of <v>       [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> upwp_old =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::upwp_old", ngrdcol,
          nzm); // Saved value of <u'w'>    [m^2/s^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vpwp_old =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::vpwp_old", ngrdcol,
          nzm); // Saved value of <v'w'>    [m^2/s^2]
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_diff_zt =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::lhs_diff_zt", ndiags3, ngrdcol,
          nzt); // Diffusion term for w'x'
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_ma_zt =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::lhs_ma_zt", ndiags3, ngrdcol,
          nzt); // Mean advection contributions to lhs
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_diff_zm =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::lhs_diff_zm", ndiags3, ngrdcol,
          nzm); // Diffusion term for w'x'
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_ma_zm =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::lhs_ma_zm", ndiags3, ngrdcol,
          nzm); // Mean advection contributions to lhs
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_ta_wprtp =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::lhs_ta_wprtp", ndiags3, ngrdcol,
          nzm); // w'r_t' turbulent advection contributions to lhs
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_ta_wpthlp =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::lhs_ta_wpthlp", ndiags3, ngrdcol,
          nzm); // w'thl' turbulent advection contributions to lhs
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_ta_wpup =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::lhs_ta_wpup", ndiags3, ngrdcol,
          nzm); // w'u' turbulent advection contributions to lhs
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_ta_wpvp =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::lhs_ta_wpvp", ndiags3, ngrdcol,
          nzm); // w'v' turbulent advection contributions to lhs
  Kokkos::View<real_t****, Kokkos::LayoutLeft> lhs_ta_wpsclrp =
      workspace_scope.view<Kokkos::View<real_t****, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::lhs_ta_wpsclrp", ndiags3, ngrdcol, nzm,
          sclr_dim); // w'sclr' turbulent advection contributions to lhs
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rhs_ta_wprtp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::rhs_ta_wprtp", ngrdcol,
          nzm); // w'r_t' turbulent advection contributions to rhs
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rhs_ta_wpthlp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::rhs_ta_wpthlp", ngrdcol,
          nzm); // w'thl' turbulent advection contributions to rhs
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rhs_ta_wpup =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::rhs_ta_wpup", ngrdcol,
          nzm); // w'u' turbulent advection contributions to rhs
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rhs_ta_wpvp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::rhs_ta_wpvp", ngrdcol,
          nzm); // w'v' turbulent advection contributions to rhs
  Kokkos::View<real_t***, Kokkos::LayoutLeft> rhs_ta_wpsclrp =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::rhs_ta_wpsclrp", ngrdcol, nzm,
          sclr_dim); // w'sclr' turbulent advection contributions to rhs
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_tp =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::lhs_tp", ndiags2, ngrdcol,
          nzm); // Turbulent production terms of w'x'
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_ta_xm =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::lhs_ta_xm", ndiags2, ngrdcol,
          nzt); // Turbulent advection terms of xm
  Kokkos::View<real_t**, Kokkos::LayoutLeft> lhs_ac_pr2 =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::lhs_ac_pr2", ngrdcol,
          nzm); // Accumulation of w'x' and w'x' pressure term 2
  Kokkos::View<real_t**, Kokkos::LayoutLeft> lhs_pr1_wprtp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::lhs_pr1_wprtp", ngrdcol,
          nzm); // Pressure term 1 for w'r_t' for all grid levels
  Kokkos::View<real_t**, Kokkos::LayoutLeft> lhs_pr1_wpthlp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::lhs_pr1_wpthlp", ngrdcol,
          nzm); // Pressure term 1 for w'thl' for all grid levels
  Kokkos::View<real_t**, Kokkos::LayoutLeft> lhs_pr1_wpsclrp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_xm_wpxp::lhs_pr1_wpsclrp", ngrdcol,
          nzm);       // Pressure term 1 for w'sclr' for all grid levels
  bool l_scalar_calc; // True if sclr_dim > 0
  int i;
  int k;
  bool l_perturbed_wind;
  const auto err_info_err_code = err_info.err_code;
  // -------------------- Begin Code --------------------
  l_perturbed_wind = l_predict_upwp_vpwp && l_linearize_pbl_winds;
  zt2zm_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(wp3, Kokkos::ALL, Kokkos::ALL),
            Kokkos::subview(wp3_zm, Kokkos::ALL, Kokkos::ALL));
  Skx_func(nzm, ngrdcol, wp2, wp3_zm, // In
           w_tol, clubb_params,       // In
           Skw_zm                     // Out
  );
  calc_wp3_on_wp2(workspace, nzm, nzt, ngrdcol, gr, wp2, wp3, // In
                  wp3_on_wp2, wp3_on_wp2_zt                   // Out
  );
  // Check whether monotonic flux limiter flags are set appropriately
  if (clubb_at_least_debug_level_api(0)) {
    if (l_mono_flux_lim_rtm && !l_mono_flux_lim_spikefix) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << err_info.err_header_global << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "l_mono_flux_lim_rtm=T with l_mono_flux_lim_spikefix=F can lead to spikes aloft."
                << "\n";
        std::cerr << message.str();
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
  }
  // Check whether the passive scalars are present.
  if (sclr_dim > 0) {
    l_scalar_calc = true;
  } else {
    l_scalar_calc = false;
  }
  if (iiPDF_type == iiPDF_new && (!l_explicit_turbulent_adv_wpxp)) {
    nrhs = 1;
  } else {
    nrhs = sclr_dim + 2;
    if (l_predict_upwp_vpwp) {
      nrhs = nrhs + 2;
      if (l_perturbed_wind) {
        nrhs = nrhs + 2;
      }
    }
  }
  // Save values of predictive fields to be printed in case of crash.
  if (l_lmm_stepping) {
    Kokkos::parallel_for(
        "advance_xm_wpxp::loop1", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i    = (cell % ngrdcol);
          const int k    = (cell / ngrdcol);
          rtm_old(i, k)  = rtm(i, k);
          thlm_old(i, k) = thlm(i, k);
        });
    Kokkos::parallel_for(
        "advance_xm_wpxp::loop2", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i      = (cell % ngrdcol);
          const int k      = (cell / ngrdcol);
          wprtp_old(i, k)  = wprtp(i, k);
          wpthlp_old(i, k) = wpthlp(i, k);
        });
    if (sclr_dim > 0) {
      Kokkos::parallel_for(
          "advance_xm_wpxp::loop3", std::max(0, sclr_dim) * std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i           = (cell % ngrdcol);
            const int k           = ((cell / ngrdcol) % nzt);
            const int sclr        = (cell / (ngrdcol * nzt));
            sclrm_old(i, k, sclr) = sclrm(i, k, sclr);
          });
      Kokkos::parallel_for(
          "advance_xm_wpxp::loop4", std::max(0, sclr_dim) * std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i             = (cell % ngrdcol);
            const int k             = ((cell / ngrdcol) % nzm);
            const int sclr          = (cell / (ngrdcol * nzm));
            wpsclrp_old(i, k, sclr) = wpsclrp(i, k, sclr);
          });
    }
    if (l_predict_upwp_vpwp) {
      Kokkos::parallel_for(
          "advance_xm_wpxp::loop5", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i  = (cell % ngrdcol);
            const int k  = (cell / ngrdcol);
            um_old(i, k) = um(i, k);
            vm_old(i, k) = vm(i, k);
          });
      Kokkos::parallel_for(
          "advance_xm_wpxp::loop6", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i    = (cell % ngrdcol);
            const int k    = (cell / ngrdcol);
            upwp_old(i, k) = upwp(i, k);
            vpwp_old(i, k) = vpwp(i, k);
          });
    }
  }
  if (!l_diag_Lscale_from_tau) {
    Kokkos::parallel_for(
        "advance_xm_wpxp::loop7", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          real_t C6rt;  // CLUBB tunable parameter C6rt
          real_t C6rtb; // CLUBB tunable parameter C6rtb
          real_t C6rtc; // CLUBB tunable parameter C6rtc
          C6rt  = clubb_params(i, iC6rt);
          C6rtb = clubb_params(i, iC6rtb);
          C6rtc = clubb_params(i, iC6rtc);
          // Compute C6 as a function of Skw
          // The if...then is just here to save compute time
          if (Kokkos::abs(C6rt - C6rtb) > Kokkos::abs(C6rt + C6rtb) * eps / 2) {
            C6rt_Skw_fnc(i, k) =
                C6rtb + (C6rt - C6rtb) * Kokkos::exp(-one_half * ((Skw_zm(i, k) / C6rtc) *
                                                                  (Skw_zm(i, k) / C6rtc)));
          } else {
            C6rt_Skw_fnc(i, k) = C6rtb;
          }
        });
    Kokkos::parallel_for(
        "advance_xm_wpxp::loop8", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          real_t C6thl;  // CLUBB tunable parameter C6thl
          real_t C6thlb; // CLUBB tunable parameter C6thlb
          real_t C6thlc; // CLUBB tunable parameter C6thlc
          C6thl  = clubb_params(i, iC6thl);
          C6thlb = clubb_params(i, iC6thlb);
          C6thlc = clubb_params(i, iC6thlc);
          if (Kokkos::abs(C6thl - C6thlb) > Kokkos::abs(C6thl + C6thlb) * eps / 2) {
            C6thl_Skw_fnc(i, k) =
                C6thlb + (C6thl - C6thlb) * Kokkos::exp(-one_half * ((Skw_zm(i, k) / C6thlc) *
                                                                     (Skw_zm(i, k) / C6thlc)));
          } else {
            C6thl_Skw_fnc(i, k) = C6thlb;
          }
        });
    // Damp C6 as a function of Lscale in stably stratified regions
    damp_coefficient(nzm, ngrdcol, gr, Kokkos::subview(clubb_params, Kokkos::ALL, iC6rt),   // In
                     C6rt_Skw_fnc,                                                          // In
                     Kokkos::subview(clubb_params, Kokkos::ALL, iC6rt_Lscale0),             // In
                     Kokkos::subview(clubb_params, Kokkos::ALL, ialtitude_threshold),       // In
                     Kokkos::subview(clubb_params, Kokkos::ALL, iwpxp_L_thresh), Lscale_zm, // In
                     C6rt_Skw_fnc                                                           // Out
    );
    damp_coefficient(nzm, ngrdcol, gr, Kokkos::subview(clubb_params, Kokkos::ALL, iC6thl),  // In
                     C6thl_Skw_fnc,                                                         // In
                     Kokkos::subview(clubb_params, Kokkos::ALL, iC6thl_Lscale0),            // In
                     Kokkos::subview(clubb_params, Kokkos::ALL, ialtitude_threshold),       // In
                     Kokkos::subview(clubb_params, Kokkos::ALL, iwpxp_L_thresh), Lscale_zm, // In
                     C6thl_Skw_fnc                                                          // Out
    );
  } else {
    Kokkos::parallel_for(
        "advance_xm_wpxp::loop9", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i         = (cell % ngrdcol);
          const int k         = (cell / ngrdcol);
          C6rt_Skw_fnc(i, k)  = clubb_params(i, iC6rt);
          C6thl_Skw_fnc(i, k) = clubb_params(i, iC6thl);
        });
  }
  // Compute C7_Skw_fnc
  if (l_use_C7_Richardson) {
    // New formulation based on Richardson number
    Kokkos::parallel_for(
        "advance_xm_wpxp::loop10", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i      = (cell % ngrdcol);
          const int k      = (cell / ngrdcol);
          C7_Skw_fnc(i, k) = Cx_fnc_Richardson(i, k);
        });
  } else {
    Kokkos::parallel_for(
        "advance_xm_wpxp::loop11", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          real_t C7;  // CLUBB tunable parameter C7
          real_t C7b; // CLUBB tunable parameter C7b
          real_t C7c; // CLUBB tunable parameter C7c
          C7  = clubb_params(i, iC7);
          C7b = clubb_params(i, iC7b);
          C7c = clubb_params(i, iC7c);
          // Compute C7 as a function of Skw
          if (Kokkos::abs(C7 - C7b) > Kokkos::abs(C7 + C7b) * eps / 2) {
            C7_Skw_fnc(i, k) =
                C7b +
                (C7 - C7b) * Kokkos::exp(-one_half * ((Skw_zm(i, k) / C7c) * (Skw_zm(i, k) / C7c)));
          } else {
            C7_Skw_fnc(i, k) = C7b;
          }
        });
    // Damp C7 as a function of Lscale in stably stratified regions
    damp_coefficient(nzm, ngrdcol, gr, Kokkos::subview(clubb_params, Kokkos::ALL, iC7),     // In
                     C7_Skw_fnc,                                                            // In
                     Kokkos::subview(clubb_params, Kokkos::ALL, iC7_Lscale0),               // In
                     Kokkos::subview(clubb_params, Kokkos::ALL, ialtitude_threshold),       // In
                     Kokkos::subview(clubb_params, Kokkos::ALL, iwpxp_L_thresh), Lscale_zm, // In
                     C7_Skw_fnc                                                             // Out
    );
  }
  if (stats.l_sample) {
    stats_update("C7_Skw_fnc", C7_Skw_fnc, // In
                 stats                     // InOut
    );
    stats_update("C6rt_Skw_fnc", C6rt_Skw_fnc, // In
                 stats                         // InOut
    );
    stats_update("C6thl_Skw_fnc", C6thl_Skw_fnc, // In
                 stats                           // InOut
    );
  }
  if (clubb_at_least_debug_level_api(0)) {
    // Assertion check for C7_Skw_fnc
    // One work item checks a column and writes its status once, avoiding races
    // between bad levels. Transfer the field only when diagnostics must print.
    int invalid_C7 = 0;
    Kokkos::parallel_reduce(
        "advance_xm_wpxp::check_C7", Kokkos::RangePolicy<>(0, ngrdcol),
        KOKKOS_LAMBDA(const int column, int& found) {
          bool invalid = false;
          for (int level = 0; level < nzm; ++level) {
            const real_t coefficient = C7_Skw_fnc(column, level);
            invalid |= coefficient > one || coefficient < zero;
          }
          if (invalid) {
            err_info_err_code(column) = clubb_fatal_error;
            found |= 1;
          }
        },
        Kokkos::BOr<int>(invalid_C7));
    if (invalid_C7) {
      auto C7_Skw_fnc_host = to_host_space(workspace, C7_Skw_fnc);
      for (k = 0; k < nzm; ++k) {
        for (i = 0; i < ngrdcol; ++i) {
          if (C7_Skw_fnc_host(i, k) > one || C7_Skw_fnc_host(i, k) < zero) {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << dump_header(err_info.err_header, i) << "\n";
            std::cerr << message.str();
          }
        }
      }
    }
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "The C7_Skw_fnc variable is outside the valid range" << "\n";
        std::cerr << message.str();
      }
      return;
    }
  }
  // Define the Coefficent of Eddy Diffusivity for the wpthlp and wprtp.
  // Kw6 is used for wpthlp and wprtp, which are located on momentum levels.
  // Kw6 is located on thermodynamic levels.
  // Kw6 = c_K6 * Kh_zt
  Kokkos::parallel_for(
      "advance_xm_wpxp::loop12", std::max(0, nzt) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i = (cell % ngrdcol);
        const int k = (cell / ngrdcol);
        Kw6(i, k)   = clubb_params(i, ic_K6) * Kh_zt(i, k);
      });
  // Find the number of grid levels, both upwards and downwards, that can
  // have an effect on the central thermodynamic level during the course of
  // one time step due to turbulent advection.  This is used as part of the
  // monotonic turbulent advection scheme.
  calc_turb_adv_range(workspace, nzm, nzt, ngrdcol, gr, dt,         // In
                      w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, // In
                      mixt_frac_zm,                                 // In
                      stats,                                        // InOut
                      low_lev_effect, high_lev_effect               // Out
  );
  // Calculate 1st pressure terms for w'r_t', w'thl', and w'sclr'.
  wpxp_term_pr1_lhs(nzm, ngrdcol, gr, C6rt_Skw_fnc, // In
                    C6thl_Skw_fnc, C7_Skw_fnc,      // In
                    invrs_tau_C6_zm, l_scalar_calc, // In
                    lhs_pr1_wprtp, lhs_pr1_wpthlp,  // Out
                    lhs_pr1_wpsclrp                 // Out
  );
  Kokkos::parallel_for(
      "advance_xm_wpxp::loop13", std::max(0, nzm) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i   = (cell % ngrdcol);
        const int k   = (cell / ngrdcol);
        C6_term(i, k) = C6rt_Skw_fnc(i, k) * invrs_tau_C6_zm(i, k);
      });
  if (stats.l_sample) {
    stats_update("C6_term", C6_term, // In
                 stats               // InOut
    );
  }
  calc_xm_wpxp_ta_terms(workspace, nzm, nzt, ngrdcol, sclr_dim, gr, wp2rtp, // In
                        wp2thlp, wp2sclrp,                                  // In
                        rho_ds_zt, invrs_rho_ds_zm, rho_ds_zm,              // In
                        sigma_sqd_w, wp3_on_wp2_zt,                         // In
                        pdf_implicit_coefs_terms,                           // In
                        iiPDF_type,                                         // In
                        l_explicit_turbulent_adv_wpxp, l_predict_upwp_vpwp, // In
                        l_scalar_calc,                                      // In
                        l_godunov_upwind_wpxp_ta,                           // In
                        stats,                                              // InOut
                        lhs_ta_wprtp, lhs_ta_wpthlp, lhs_ta_wpup,           // Out
                        lhs_ta_wpvp, lhs_ta_wpsclrp,                        // Out
                        rhs_ta_wprtp, rhs_ta_wpthlp, rhs_ta_wpup,           // Out
                        rhs_ta_wpvp, rhs_ta_wpsclrp                         // Out
  );
  // Calculate various terms that are the same between all LHS matricies
  calc_xm_wpxp_lhs_terms(workspace, nzm, nzt, ngrdcol, gr, wm_zm, wm_zt, wp2, // In
                         Kh_zm, stability_correction, Kw6, C7_Skw_fnc,        // In
                         invrs_rho_ds_zt, invrs_rho_ds_zm, rho_ds_zt,         // In
                         rho_ds_zm, l_implemented, nu_vert_res_dep,           // In
                         l_diffuse_rtm_and_thlm,                              // In
                         l_stability_correct_Kh_N2_zm,                        // In
                         l_upwind_xm_ma,                                      // In
                         lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm,      // Out
                         lhs_tp, lhs_ta_xm, lhs_ac_pr2                        // Out
  );
  // Setup and decompose matrix for each variable.
  if ((iiPDF_type == iiPDF_new) && (!l_explicit_turbulent_adv_wpxp)) {
    // LHS matrices are unique, multiple band solves required
    solve_xm_wpxp_with_multiple_lhs(workspace, nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, dt, // In
                                    l_iter, nrhs, wm_zt, wp2,                                 // In
                                    rtpthvp, rtm_forcing, wprtp_forcing, thlpthvp,            // In
                                    thlm_forcing, wpthlp_forcing, rho_ds_zm,                  // In
                                    rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,              // In
                                    thv_ds_zm, rtp2, thlp2, l_implemented,                    // In
                                    sclrpthvp, sclrm_forcing, sclrp2,                         // In
                                    low_lev_effect, high_lev_effect, C7_Skw_fnc,              // In
                                    lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm,           // In
                                    lhs_ta_wprtp, lhs_ta_wpthlp, lhs_ta_wpsclrp,              // In
                                    rhs_ta_wprtp, rhs_ta_wpthlp, rhs_ta_wpsclrp,              // In
                                    lhs_tp, lhs_ta_xm, lhs_ac_pr2, lhs_pr1_wprtp,             // In
                                    lhs_pr1_wpthlp, lhs_pr1_wpsclrp,                          // In
                                    penta_solve_method,                                       // In
                                    tridiag_solve_method,                                     // In
                                    fill_holes_type,                                          // In
                                    l_diffuse_rtm_and_thlm,                                   // In
                                    l_upwind_xm_ma,                                           // In
                                    l_enable_relaxed_clipping,                                // In
                                    l_mono_flux_lim_thlm,                                     // In
                                    l_mono_flux_lim_rtm,                                      // In
                                    l_mono_flux_lim_um,                                       // In
                                    l_mono_flux_lim_vm,                                       // In
                                    l_mono_flux_lim_spikefix,                                 // In
                                    wprtp_cl_num, wpthlp_cl_num,              // InOut
                                    stats,                                    // InOut
                                    rtm, wprtp, thlm, wpthlp, sclrm, wpsclrp, // InOut
                                    err_info                                  // InOut
    );
  } else {
    // LHS matrices are equivalent, only one solve required
    solve_xm_wpxp_with_single_lhs(workspace, nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, dt, // In
                                  l_iter, nrhs,                                             // In
                                  wm_zt, wp2, invrs_tau_C6_zm, tau_max_zm,                  // In
                                  rtpthvp, rtm_forcing, wprtp_forcing, thlpthvp,            // In
                                  thlm_forcing, wpthlp_forcing, rho_ds_zm,                  // In
                                  rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,              // In
                                  thv_ds_zm, rtp2, thlp2, l_implemented,                    // In
                                  sclrpthvp, sclrm_forcing, sclrp2, um_forcing,             // In
                                  vm_forcing, ug, vg, uprcp, vprcp, rc_coef_zm, fcor,       // In
                                  fcor_y, up2, vp2,                                         // In
                                  low_lev_effect, high_lev_effect,                          // In
                                  C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc,                  // In
                                  lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm,           // In
                                  lhs_ta_wprtp,                                             // In
                                  rhs_ta_wprtp, rhs_ta_wpthlp, rhs_ta_wpup,                 // In
                                  rhs_ta_wpvp, rhs_ta_wpsclrp,                              // In
                                  lhs_tp, lhs_ta_xm, lhs_ac_pr2, lhs_pr1_wprtp,             // In
                                  lhs_pr1_wpthlp, lhs_pr1_wpsclrp,                          // In
                                  Kokkos::subview(clubb_params, Kokkos::ALL, iC_uu_shr),    // In
                                  penta_solve_method,                                       // In
                                  tridiag_solve_method,                                     // In
                                  fill_holes_type,                                          // In
                                  l_predict_upwp_vpwp,                                      // In
                                  l_ho_nontrad_coriolis,                                    // In
                                  l_ho_trad_coriolis,                                       // In
                                  l_diffuse_rtm_and_thlm,                                   // In
                                  l_upwind_xm_ma,                                           // In
                                  l_enable_relaxed_clipping,                                // In
                                  l_perturbed_wind,                                         // In
                                  l_mono_flux_lim_thlm,                                     // In
                                  l_mono_flux_lim_rtm,                                      // In
                                  l_mono_flux_lim_um,                                       // In
                                  l_mono_flux_lim_vm,                                       // In
                                  l_mono_flux_lim_spikefix,                                 // In
                                  wprtp_cl_num, wpthlp_cl_num,                              // InOut
                                  upwp_cl_num, vpwp_cl_num,                                 // InOut
                                  stats,                                                    // InOut
                                  rtm, wprtp, thlm, wpthlp,                                 // InOut
                                  sclrm, wpsclrp, um, upwp, vm, vpwp,                       // InOut
                                  um_pert, vm_pert, upwp_pert, vpwp_pert, err_info          // InOut
    );
  }
  if (l_lmm_stepping) {
    Kokkos::parallel_for(
        "advance_xm_wpxp::loop14", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          thlm(i, k)  = one_half * (thlm_old(i, k) + thlm(i, k));
          rtm(i, k)   = one_half * (rtm_old(i, k) + rtm(i, k));
        });
    Kokkos::parallel_for(
        "advance_xm_wpxp::loop15", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i  = (cell % ngrdcol);
          const int k  = (cell / ngrdcol);
          wpthlp(i, k) = one_half * (wpthlp_old(i, k) + wpthlp(i, k));
          wprtp(i, k)  = one_half * (wprtp_old(i, k) + wprtp(i, k));
        });
    if (sclr_dim > 0) {
      Kokkos::parallel_for(
          "advance_xm_wpxp::loop16",
          std::max(0, sclr_dim) * std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i       = (cell % ngrdcol);
            const int k       = ((cell / ngrdcol) % nzt);
            const int sclr    = (cell / (ngrdcol * nzt));
            sclrm(i, k, sclr) = one_half * (sclrm_old(i, k, sclr) + sclrm(i, k, sclr));
          });
      Kokkos::parallel_for(
          "advance_xm_wpxp::loop17",
          std::max(0, sclr_dim) * std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i         = (cell % ngrdcol);
            const int k         = ((cell / ngrdcol) % nzm);
            const int sclr      = (cell / (ngrdcol * nzm));
            wpsclrp(i, k, sclr) = one_half * (wpsclrp_old(i, k, sclr) + wpsclrp(i, k, sclr));
          });
    }
    if (l_predict_upwp_vpwp) {
      Kokkos::parallel_for(
          "advance_xm_wpxp::loop18", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            um(i, k)    = one_half * (um_old(i, k) + um(i, k));
            vm(i, k)    = one_half * (vm_old(i, k) + vm(i, k));
          });
      Kokkos::parallel_for(
          "advance_xm_wpxp::loop19", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            upwp(i, k)  = one_half * (upwp_old(i, k) + upwp(i, k));
            vpwp(i, k)  = one_half * (vpwp_old(i, k) + vpwp(i, k));
          });
    }
  }
  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << err_info.err_header_global << "\n";
        std::cerr << message.str();
      }
      for (i = 0; i < ngrdcol; i += 1) {
        error_prints_xm_wpxp(nzm, nzt, sclr_dim, Kokkos::subview(gr.zm, i, Kokkos::ALL),   // In
                             Kokkos::subview(gr.zt, i, Kokkos::ALL),                       // In
                             dt, Kokkos::subview(sigma_sqd_w, i, Kokkos::ALL),             // In
                             Kokkos::subview(wm_zm, i, Kokkos::ALL),                       // In
                             Kokkos::subview(wm_zt, i, Kokkos::ALL),                       // In
                             Kokkos::subview(wp2, i, Kokkos::ALL),                         // In
                             Kokkos::subview(Lscale_zm, i, Kokkos::ALL),                   // In
                             Kokkos::subview(wp3_on_wp2, i, Kokkos::ALL),                  // In
                             Kokkos::subview(wp3_on_wp2_zt, i, Kokkos::ALL),               // In
                             Kokkos::subview(Kh_zt, i, Kokkos::ALL),                       // In
                             Kokkos::subview(Kh_zm, i, Kokkos::ALL),                       // In
                             Kokkos::subview(invrs_tau_C6_zm, i, Kokkos::ALL),             // In
                             Kokkos::subview(Skw_zm, i, Kokkos::ALL),                      // In
                             Kokkos::subview(wp2rtp, i, Kokkos::ALL),                      // In
                             Kokkos::subview(rtpthvp, i, Kokkos::ALL),                     // In
                             Kokkos::subview(rtm_forcing, i, Kokkos::ALL),                 // In
                             Kokkos::subview(wprtp_forcing, i, Kokkos::ALL),               // In
                             Kokkos::subview(rtm_ref, i, Kokkos::ALL),                     // In
                             Kokkos::subview(wp2thlp, i, Kokkos::ALL),                     // In
                             Kokkos::subview(thlpthvp, i, Kokkos::ALL),                    // In
                             Kokkos::subview(thlm_forcing, i, Kokkos::ALL),                // In
                             Kokkos::subview(wpthlp_forcing, i, Kokkos::ALL),              // In
                             Kokkos::subview(thlm_ref, i, Kokkos::ALL),                    // In
                             Kokkos::subview(rho_ds_zm, i, Kokkos::ALL),                   // In
                             Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),                   // In
                             Kokkos::subview(invrs_rho_ds_zm, i, Kokkos::ALL),             // In
                             Kokkos::subview(invrs_rho_ds_zt, i, Kokkos::ALL),             // In
                             Kokkos::subview(thv_ds_zm, i, Kokkos::ALL),                   // In
                             Kokkos::subview(rtp2, i, Kokkos::ALL),                        // In
                             Kokkos::subview(thlp2, i, Kokkos::ALL),                       // In
                             Kokkos::subview(w_1_zm, i, Kokkos::ALL),                      // In
                             Kokkos::subview(w_2_zm, i, Kokkos::ALL),                      // In
                             Kokkos::subview(varnce_w_1_zm, i, Kokkos::ALL),               // In
                             Kokkos::subview(varnce_w_2_zm, i, Kokkos::ALL),               // In
                             Kokkos::subview(mixt_frac_zm, i, Kokkos::ALL), l_implemented, // In
                             Kokkos::subview(em, i, Kokkos::ALL),                          // In
                             Kokkos::subview(wp2sclrp, i, Kokkos::ALL, Kokkos::ALL),       // In
                             Kokkos::subview(sclrpthvp, i, Kokkos::ALL, Kokkos::ALL),      // In
                             Kokkos::subview(sclrm_forcing, i, Kokkos::ALL, Kokkos::ALL),  // In
                             Kokkos::subview(sclrp2, i, Kokkos::ALL, Kokkos::ALL),         // In
                             Kokkos::subview(Cx_fnc_Richardson, i, Kokkos::ALL),           // In
                             pdf_implicit_coefs_terms,                                     // In
                             Kokkos::subview(um_forcing, i, Kokkos::ALL),                  // In
                             Kokkos::subview(vm_forcing, i, Kokkos::ALL),                  // In
                             Kokkos::subview(ug, i, Kokkos::ALL),                          // In
                             Kokkos::subview(vg, i, Kokkos::ALL),                          // In
                             Kokkos::subview(wpthvp, i, Kokkos::ALL),                      // In
                             to_host_space(workspace, fcor)(i),                            // In
                             Kokkos::subview(um_ref, i, Kokkos::ALL),                      // In
                             Kokkos::subview(vm_ref, i, Kokkos::ALL),                      // In
                             Kokkos::subview(up2, i, Kokkos::ALL),                         // In
                             Kokkos::subview(vp2, i, Kokkos::ALL),                         // In
                             Kokkos::subview(uprcp, i, Kokkos::ALL),                       // In
                             Kokkos::subview(vprcp, i, Kokkos::ALL),                       // In
                             Kokkos::subview(rc_coef_zm, i, Kokkos::ALL),                  // In
                             Kokkos::subview(rtm, i, Kokkos::ALL),                         // In
                             Kokkos::subview(wprtp, i, Kokkos::ALL),                       // In
                             Kokkos::subview(thlm, i, Kokkos::ALL),                        // In
                             Kokkos::subview(wpthlp, i, Kokkos::ALL),                      // In
                             Kokkos::subview(sclrm, i, Kokkos::ALL, Kokkos::ALL),          // In
                             Kokkos::subview(wpsclrp, i, Kokkos::ALL, Kokkos::ALL),        // In
                             Kokkos::subview(um, i, Kokkos::ALL),                          // In
                             Kokkos::subview(upwp, i, Kokkos::ALL),                        // In
                             Kokkos::subview(vm, i, Kokkos::ALL),                          // In
                             Kokkos::subview(vpwp, i, Kokkos::ALL),                        // In
                             Kokkos::subview(rtm_old, i, Kokkos::ALL),                     // In
                             Kokkos::subview(wprtp_old, i, Kokkos::ALL),                   // In
                             Kokkos::subview(thlm_old, i, Kokkos::ALL),                    // In
                             Kokkos::subview(wpthlp_old, i, Kokkos::ALL),                  // In
                             Kokkos::subview(sclrm_old, i, Kokkos::ALL, Kokkos::ALL),      // In
                             Kokkos::subview(wpsclrp_old, i, Kokkos::ALL, Kokkos::ALL),    // In
                             Kokkos::subview(um_old, i, Kokkos::ALL),                      // In
                             Kokkos::subview(upwp_old, i, Kokkos::ALL),                    // In
                             Kokkos::subview(vm_old, i, Kokkos::ALL),                      // In
                             Kokkos::subview(vpwp_old, i, Kokkos::ALL),                    // In
                             l_predict_upwp_vpwp, l_lmm_stepping                           // In
        );
      }
    }
  }
  // Native rtm_sponge_damp_settings%l_sponge_damping:
  // Sponge damping is excluded with the JAX scope. Standalone initialization
  // rejects sponge-enabled cases before model allocation.
  // Native thlm_sponge_damp_settings%l_sponge_damping:
  // Sponge damping is excluded with the JAX scope. Standalone initialization
  // rejects sponge-enabled cases before model allocation.
  if (l_predict_upwp_vpwp) {
    // Native uv_sponge_damp_settings%l_sponge_damping:
    // Sponge damping is excluded with the JAX scope. Standalone initialization
    // rejects sponge-enabled cases before model allocation.
    // Adjust um and vm if nudging is turned on.
    if (l_uv_nudge) {
      // Reflect nudging in budget
      if (stats.l_sample) {
        {
          Real2D stats_values_20 = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                                "advance_xm_wpxp::stats_values_20",
                                                                um.extent(0), um.extent(1));
          Kokkos::parallel_for(
              "advance_xm_wpxp::stats_value", um.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = (cell % um.extent(0));
                const int k           = (cell / um.extent(0));
                stats_values_20(i, k) = um(i, k) / dt;
              });
          stats_begin_budget("um_ndg", stats_values_20, // In
                             stats                      // InOut
          );
        }
        {
          Real2D stats_values_21 = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                                "advance_xm_wpxp::stats_values_21",
                                                                vm.extent(0), vm.extent(1));
          Kokkos::parallel_for(
              "advance_xm_wpxp::stats_value", vm.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = (cell % vm.extent(0));
                const int k           = (cell / vm.extent(0));
                stats_values_21(i, k) = vm(i, k) / dt;
              });
          stats_begin_budget("vm_ndg", stats_values_21, // In
                             stats                      // InOut
          );
        }
      }
      Kokkos::parallel_for(
          "advance_xm_wpxp::loop22", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            um(i, k)    = um(i, k) - ((um(i, k) - um_ref(i, k)) * (dt / ts_nudge));
            vm(i, k)    = vm(i, k) - ((vm(i, k) - vm_ref(i, k)) * (dt / ts_nudge));
          });
      // Reflect nudging in budget
      if (stats.l_sample) {
        {
          Real2D stats_values_23 = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                                "advance_xm_wpxp::stats_values_23",
                                                                um.extent(0), um.extent(1));
          Kokkos::parallel_for(
              "advance_xm_wpxp::stats_value", um.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = (cell % um.extent(0));
                const int k           = (cell / um.extent(0));
                stats_values_23(i, k) = um(i, k) / dt;
              });
          stats_finalize_budget("um_ndg", stats_values_23, // In
                                stats                      // InOut
          );
        }
        {
          Real2D stats_values_24 = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                                "advance_xm_wpxp::stats_values_24",
                                                                vm.extent(0), vm.extent(1));
          Kokkos::parallel_for(
              "advance_xm_wpxp::stats_value", vm.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = (cell % vm.extent(0));
                const int k           = (cell / vm.extent(0));
                stats_values_24(i, k) = vm(i, k) / dt;
              });
          stats_finalize_budget("vm_ndg", stats_values_24, // In
                                stats                      // InOut
          );
        }
      }
    }
    if (stats.l_sample) {
      stats_update("um_ref", um_ref, // In
                   stats             // InOut
      );
      stats_update("vm_ref", vm_ref, // In
                   stats             // InOut
      );
    }
  }
  return;
}

namespace {
// dt: Timestep                                  [s]
// l_implemented: True if CLUBB is being implemented and run in a host model
// lhs_diff_zm: Diffusion term for w'x'
// lhs_diff_zt: Diffusion term for xm
// lhs_ma_zm: Mean advection contributions to lhs
// lhs_ma_zt: Mean advection contributions to lhs
// lhs_ta_wpxp: Turbulent advection contributions to lhs
// lhs_ta_xm: Turbulent advection terms of xm
// lhs_tp: Turbulent production terms of w'x'
// lhs_pr1: Pressure term 1 for w'x'
// lhs_ac_pr2: Accumulation of w'x' and w'x' pressure term 2
// l_diffuse_rtm_and_thlm: This flag determines whether or not we want CLUBB to do diffusion on
// rtm and thlm
// lhs: Implicit contributions to wpxp/xm (band diag. matrix) (LAPACK)
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9>
void xm_wpxp_lhs(int nzm, int nzt, int ngrdcol, bool l_iter, real_t dt,                    // In
                 bool l_implemented, Array0 lhs_diff_zm, Array1 lhs_diff_zt,               // In
                 Array2 lhs_ma_zm, Array3 lhs_ma_zt, Array4 lhs_ta_wpxp, Array5 lhs_ta_xm, // In
                 Array6 lhs_tp, Array7 lhs_pr1, Array8 lhs_ac_pr2,                         // In
                 bool l_diffuse_rtm_and_thlm,                                              // In
                 Array9 lhs                                                                // Out
) {
  // Description:
  //   Compute LHS band diagonal matrix for xm and w'x'.
  //   This subroutine computes the implicit portion of
  //   the xm and w'x' equations.
  //
  //
  // Notes:
  //
  //   Boundary conditions:
  //       The turbulent flux (wpxp) use fixed-point boundary conditions at both the
  //       upper and lower boundaries.  Therefore, anything set in the wpxp loop
  //       at both the upper and lower boundaries would be overwritten here.
  //       However, the wpxp loop does not extend to the boundary levels.  An array
  //       with a value of 1 at the main diagonal on the left-hand side and with
  //       values of 0 at all other diagonals on the left-hand side will preserve the
  //       right-hand side value at that level.
  //
  //           wpxp(1) ... wpxp(nzm)
  //         [   0.0         0.0    ]
  //         [   0.0         0.0    ]
  //         [   1.0   ...   1.0    ]
  //         [   0.0         0.0    ]
  //         [   0.0         0.0    ]
  //
  //
  //   LHS turbulent advection (ta) term:
  //        An "over-implicit" weighted time step is applied to this term.
  //        The weight of the implicit portion of this term is controlled by
  //        the factor gamma_over_implicit_ts (abbreviated "gamma" in the
  //        equation in order to balance a weight that is not equal to 1,
  //        such that:
  //             -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
  //        where X is the variable that is being solved for in a predictive
  //        equation (<w'x'> in this case), y(t) is the linearized portion of
  //        the term that gets treated implicitly, and RHS is the portion of
  //        the term that is always treated explicitly.  A weight of greater
  //        than 1 can be applied to make the term more numerically stable.
  //
  //
  //   xm: Left-hand side (implicit xm portion of the code).
  //
  //   Thermodynamic subdiagonal (lhs index: t_km1_tdiag)
  //         [ x xm(k-1,<t+1>) ]
  //   Momentum subdiagonal (lhs index: t_k_mdiag)
  //         [ x wpxp(k,<t+1>) ]
  //   Thermodynamic main diagonal (lhs index: t_k_tdiag)
  //         [ x xm(k,<t+1>) ]
  //   Momentum superdiagonal (lhs index: t_kp1_mdiag)
  //         [ x wpxp(k+1,<t+1>) ]
  //   Thermodynamic superdiagonal (lhs index: t_kp1_tdiag)
  //         [ x xm(k+1,<t+1>) ]
  //
  //
  //   w'x': Left-hand side (implicit w'x' portion of the code).
  //
  //   Momentum subdiagonal (lhs index: m_km1_mdiag)
  //         [ x wpxp(k-1,<t+1>) ]
  //   Thermodynamic subdiagonal (lhs index: m_km1_tdiag)
  //         [ x xm(k-1,<t+1>) ]
  //   Momentum main diagonal (lhs index: m_k_mdiag)
  //         [ x wpxp(k,<t+1>) ]
  //   Thermodynamic superdiagonal (lhs index: m_k_tdiag)
  //         [ x xm(k,<t+1>) ]
  //   Momentum superdiagonal (lhs index: m_kp1_mdiag)
  //         [ x wpxp(k+1,<t+1>) ]
  //
  //----------------------------------------------------------------------------------
  //------------------- Input Variables -------------------
  //------------------- Output Variable -------------------
  //------------------- Local Variables -------------------
  // Indices
  // on rtm and thlm
  real_t invrs_dt;
  //------------------- Begin Code -------------------
  // Initializations/precalculations
  invrs_dt = real_t(1.0) / dt;
  // Lower (upper) boundary for w'x' on ascending (descending) grid.
  Kokkos::parallel_for(
      "xm_wpxp_lhs::loop1", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        const int i  = cell;
        lhs(0, i, 0) = real_t(0.0);
        lhs(1, i, 0) = real_t(0.0);
        lhs(2, i, 0) = real_t(1.0);
        lhs(3, i, 0) = real_t(0.0);
        lhs(4, i, 0) = real_t(0.0);
      });
  // Combine xm and w'x' terms into LHS
  Kokkos::parallel_for(
      "xm_wpxp_lhs::loop2", std::max(0, nzt) * std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        const int i = (cell % ngrdcol);
        const int k = (cell / ngrdcol);
        int k_xm;
        k_xm = 2 * k + 1;
        // ---- sum xm terms ----
        lhs(0, i, k_xm) = zero;
        lhs(1, i, k_xm) = lhs_ta_xm(0, i, k);
        lhs(2, i, k_xm) = invrs_dt;
        lhs(3, i, k_xm) = lhs_ta_xm(1, i, k);
        lhs(4, i, k_xm) = zero;
      });
  Kokkos::parallel_for(
      "xm_wpxp_lhs::loop3", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i = (cell % ngrdcol);
        const int k = (cell / ngrdcol) + 1;
        int k_wpxp;
        k_wpxp = 2 * k;
        // ---- sum w'x' terms ----
        lhs(0, i, k_wpxp) = lhs_ma_zm(0, i, k) + lhs_diff_zm(0, i, k) +
                            gamma_over_implicit_ts * lhs_ta_wpxp(0, i, k);
        lhs(1, i, k_wpxp) = lhs_tp(0, i, k);
        // Preserve the Native CPU contraction before adding MA and diffusion (BFB-017).
        const real_t lhs_ma_diff = lhs_ma_zm(1, i, k) + lhs_diff_zm(1, i, k);
        const real_t lhs_ta_pr1  = lhs_ta_wpxp(1, i, k) + lhs_pr1(i, k);
        lhs(2, i, k_wpxp) =
            Kokkos::fma(gamma_over_implicit_ts, lhs_ta_pr1, lhs_ac_pr2(i, k)) + lhs_ma_diff;
        lhs(3, i, k_wpxp) = lhs_tp(1, i, k);
        lhs(4, i, k_wpxp) = lhs_ma_zm(2, i, k) + lhs_diff_zm(2, i, k) +
                            gamma_over_implicit_ts * lhs_ta_wpxp(2, i, k);
      });
  // Upper (lower) boundary for w'x' on ascending (descending) grid.
  // These were set in the loop above for simplicity, so they must be set properly here
  Kokkos::parallel_for(
      "xm_wpxp_lhs::loop4", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        const int i            = cell;
        lhs(0, i, 2 * nzm - 2) = real_t(0.0);
        lhs(1, i, 2 * nzm - 2) = real_t(0.0);
        lhs(2, i, 2 * nzm - 2) = real_t(1.0);
        lhs(3, i, 2 * nzm - 2) = real_t(0.0);
        lhs(4, i, 2 * nzm - 2) = real_t(0.0);
      });
  // LHS time tendency
  if (l_iter) {
    Kokkos::parallel_for(
        "xm_wpxp_lhs::loop5", std::max(0, nzm - 2) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol) + 1;
          int k_wpxp;
          k_wpxp            = 2 * k;
          lhs(2, i, k_wpxp) = lhs(2, i, k_wpxp) + invrs_dt;
        });
  }
  // Calculate diffusion terms for all thermodynamic grid level
  if (l_diffuse_rtm_and_thlm) {
    Kokkos::parallel_for(
        "xm_wpxp_lhs::loop6", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          int k_xm;
          k_xm            = 2 * k + 1;
          lhs(0, i, k_xm) = lhs(0, i, k_xm) + lhs_diff_zt(0, i, k);
          lhs(2, i, k_xm) = lhs(2, i, k_xm) + lhs_diff_zt(1, i, k);
          lhs(4, i, k_xm) = lhs(4, i, k_xm) + lhs_diff_zt(2, i, k);
        });
  }
  // Calculate mean advection terms for all momentum grid level
  if (!l_implemented) {
    Kokkos::parallel_for(
        "xm_wpxp_lhs::loop7", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          int k_xm;
          k_xm            = 2 * k + 1;
          lhs(0, i, k_xm) = lhs(0, i, k_xm) + lhs_ma_zt(0, i, k);
          lhs(2, i, k_xm) = lhs(2, i, k_xm) + lhs_ma_zt(1, i, k);
          lhs(4, i, k_xm) = lhs(4, i, k_xm) + lhs_ma_zt(2, i, k);
        });
  }
  return;
}

// wm_zm: w wind component on momentum levels      [m/s]
// wm_zt: w wind component on thermo. levels       [m/s]
// wp2: w'^2 (momentum levels)               [m^2/s^2]
// Kh_zm: Eddy diffusivity on momentum levels    [m^2/s]
// stability_correction: Stability correction factor
// Kw6: Coef. of eddy diffusivity for w'x'     [m^2/s]
// C7_Skw_fnc: C_7 parameter with Sk_w applied            [-]
// invrs_rho_ds_zt: Inv. dry, static density at t-levs.   [m^3/kg]
// rho_ds_zm: Dry, static density on momentum levs. [kg/m^3]
// l_implemented: True if CLUBB is being implemented and run in a host model
// nu_vert_res_dep: Vertical resolution dependent nu values
// l_diffuse_rtm_and_thlm: This flag determines whether or not we want CLUBB to do diffusion on
// rtm and thlm
// l_stability_correct_Kh_N2_zm: This flag determines whether or not we want CLUBB to apply a
// stability correction
// l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent or mean advection terms. It
// affects rtm, thlm, sclrm, um and vm.
// lhs_diff_zm: Diffusion term for w'x'
// lhs_diff_zt: Diffusion term for xm
// lhs_ma_zt: Mean advection contributions to lhs
// lhs_ma_zm: Mean advection contributions to lhs
// lhs_tp: Turbulent production terms of w'x'
// lhs_ta_xm: Turbulent advection terms of xm
// lhs_ac_pr2: Accumulation of w'x' and w'x' pressure term 2
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17>
void calc_xm_wpxp_lhs_terms(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol,
                            const grid& gr, Array0 wm_zm,
                            Array1 wm_zt, // In
                            Array2 wp2,   // In
                            Array3 Kh_zm, Array4 stability_correction, Array5 Kw6,
                            Array6 C7_Skw_fnc,                                                // In
                            Array7 invrs_rho_ds_zt, Array8 invrs_rho_ds_zm, Array9 rho_ds_zt, // In
                            Array10 rho_ds_zm, bool l_implemented,                            // In
                            const nu_vertical_res_dep& nu_vert_res_dep,                       // In
                            bool l_diffuse_rtm_and_thlm,                                      // In
                            bool l_stability_correct_Kh_N2_zm,                                // In
                            bool l_upwind_xm_ma,                                              // In
                            Array11 lhs_diff_zm, Array12 lhs_diff_zt, Array13 lhs_ma_zt,
                            Array14 lhs_ma_zm,                                    // Out
                            Array15 lhs_tp, Array16 lhs_ta_xm, Array17 lhs_ac_pr2 // Out
) {
  // Description:
  //   Calculate various xm and w'x' terms. These are general terms that are the same
  //   for multiple LHS matrices, so we save computations by calculating them once
  //   here, then reusing them where needed.
  //
  //-------------------------------------------------------------------------------------------
  //------------------- Input Variables -------------------
  // diffusion on rtm and thlm
  // a stability correction
  // differencing approximation rather than a centered
  // differencing for turbulent or mean advection terms.
  // It affects rtm, thlm, sclrm, um and vm.
  //------------------- Output Variables -------------------
  //------------------- Local Variables -------------------
  auto workspace_scope = workspace.scope();
  Kokkos::View<real_t**, Kokkos::LayoutLeft> Kh_N2_zm =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_lhs_terms::Kh_N2_zm", ngrdcol,
          nzm); // Stability correction factor
  Kokkos::View<real_t**, Kokkos::LayoutLeft> K_zm =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_lhs_terms::K_zm", ngrdcol,
          nzm); // Coef. of eddy diffusivity at momentum level (k)   [m^2/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> Kw6_zm =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_lhs_terms::Kw6_zm", ngrdcol,
          nzm); // Eddy diffusivity coefficient, momentum levels [m2/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> K_zt =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_lhs_terms::K_zt", ngrdcol,
          nzt);       // Eddy diffusivity coefficient, thermo. levels [m2/s]
  real_t constant_nu; // controls the magnitude of diffusion
  Kokkos::View<real_t*, Kokkos::LayoutLeft> zeros_array =
      workspace_scope.view<Kokkos::View<real_t*, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_lhs_terms::zeros_array", ngrdcol);
  //------------------- Begin Code -------------------
  // Initializations/precalculations
  constant_nu = real_t(0.1);
  zt2zm_api(nzm, nzt, ngrdcol, gr, Kw6, Kw6_zm, zero_threshold);
  // Calculate turbulent advection terms of xm for all grid levels
  xm_term_ta_lhs(nzm, nzt, ngrdcol, gr,      // In
                 rho_ds_zm, invrs_rho_ds_zt, // In
                 lhs_ta_xm                   // Out
  );
  // Calculate turbulent production terms of w'x' for all grid level
  wpxp_term_tp_lhs(nzm, ngrdcol, gr, wp2, // In
                   lhs_tp                 // Out
  );
  // Calculate accumulation of w'x' and w'x' pressure term 2 of w'x' for all grid level
  // https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:wpxp_pr
  wpxp_terms_ac_pr2_lhs(nzm, nzt, ngrdcol, gr, C7_Skw_fnc, // In
                        wm_zt, gr.invrs_dzm,               // In
                        lhs_ac_pr2                         // Out
  );
  // Calculate diffusion terms for all momentum grid level
  diffusion_zm_lhs(workspace, nzm, nzt, ngrdcol, gr, Kw6, Kw6_zm, nu_vert_res_dep.nu6, // In
                   invrs_rho_ds_zm, rho_ds_zt,                                         // In
                   lhs_diff_zm                                                         // Out
  );
  // Calculate mean advection terms for all momentum grid level
  term_ma_zm_lhs(nzm, nzt, ngrdcol, wm_zm,       // In
                 gr.invrs_dzm, gr.weights_zm2zt, // In
                 lhs_ma_zm                       // Out
  );
  // Calculate diffusion terms for all thermodynamic grid level
  if (l_diffuse_rtm_and_thlm) {
    if (l_stability_correct_Kh_N2_zm) {
      Kokkos::parallel_for(
          "calc_xm_wpxp_lhs_terms::loop1", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i    = (cell % ngrdcol);
            const int k    = (cell / ngrdcol);
            Kh_N2_zm(i, k) = Kh_zm(i, k) / stability_correction(i, k);
          });
    } else {
      Kokkos::parallel_for(
          "calc_xm_wpxp_lhs_terms::loop2", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i    = (cell % ngrdcol);
            const int k    = (cell / ngrdcol);
            Kh_N2_zm(i, k) = Kh_zm(i, k);
          });
    }
    Kokkos::parallel_for(
        "calc_xm_wpxp_lhs_terms::loop3", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          K_zm(i, k)  = Kh_N2_zm(i, k) + constant_nu;
        });
    zm2zt_api(nzm, nzt, ngrdcol, gr, K_zm, K_zt, zero_threshold);
    Kokkos::parallel_for(
        "calc_xm_wpxp_lhs_terms::loop4", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i    = cell;
          zeros_array(i) = zero;
        });
    diffusion_zt_lhs(workspace, nzm, nzt, ngrdcol, gr, K_zm, K_zt, zeros_array, // In
                     invrs_rho_ds_zt, rho_ds_zm,                                // In
                     lhs_diff_zt                                                // Out
    );
  }
  // Calculate mean advection terms for all thermodynamic grid level
  if (!l_implemented) {
    term_ma_zt_lhs(nzm, nzt, ngrdcol, wm_zt, gr.weights_zt2zm, // In
                   gr.invrs_dzt, gr.invrs_dzm,                 // In
                   l_upwind_xm_ma, gr.grid_dir,                // In
                   lhs_ma_zt                                   // Out
    );
  }
  return;
}

// solve_type: Variables being solved for.
// dt: Timestep                                  [s]
// xm: xm (thermodynamic levels)               [x un]
// wpxp: <w'x'> (momentum levels)          [{x un} m/s]
// xm_forcing: xm forcings (thermodynamic levels)  [{x un}/s]
// wpxp_forcing: <w'x'> forcing (momentum levs)  [{x un} m/s^2]
// C7_Skw_fnc: C_7 parameter with Sk_w applied            [-]
// xpthvp: x'th_v' (momentum levels)           [{x un} K]
// thv_ds_zm: Dry, base-state theta_v on mom. levs.      [K]
// lhs_pr1: Pressure term 1 for w'x'
// lhs_ta_wpxp: Turbulent advection terms of w'x'
// rhs: Right-hand side of band diag. matrix. (LAPACK)
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10>
void xm_wpxp_rhs(ScratchWorkspace& workspace, // InOut
                 int nzm, int nzt, int ngrdcol, const grid& gr, int solve_type, bool l_iter, // In
                 real_t dt,                                                                  // In
                 Array0 xm, Array1 wpxp, Array2 xm_forcing, Array3 wpxp_forcing,             // In
                 Array4 C7_Skw_fnc,                                                          // In
                 Array5 xpthvp, Array6 rhs_ta, Array7 thv_ds_zm,                             // In
                 Array8 lhs_pr1, Array9 lhs_ta_wpxp,                                         // In
                 stats_type& stats, // InOut
                 Array10 rhs        // Out
) {
  // Description:
  // Compute RHS vector for xm and w'x'.
  // This subroutine computes the explicit portion of
  // the xm and w'x' equations.
  //
  // Notes:
  //   For LHS turbulent advection (ta) term.
  //       An "over-implicit" weighted time step is applied to this term.
  //       The weight of the implicit portion of this term is controlled by
  //       the factor gamma_over_implicit_ts (abbreviated "gamma" in the
  //       expression below).  A factor is added to the right-hand side of
  //       the equation in order to balance a weight that is not equal to 1,
  //       such that:
  //            -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
  //       where X is the variable that is being solved for in a predictive
  //       equation (<w'x'> in this case), y(t) is the linearized portion of
  //       the term that gets treated implicitly, and RHS is the portion of
  //       the term that is always treated explicitly.  A weight of greater
  //       than 1 can be applied to make the term more numerically stable.
  //
  //   --- THIS SUBROUTINE HAS BEEN OPTIMIZED ---
  //   Significant changes to this routine may adversely affect computational speed
  //       - Gunther Huebler, Aug. 2018, clubb:ticket:834
  //----------------------------------------------------------------------------------
  //------------------- Input Variables -------------------
  // For "over-implicit" weighted time step.
  // This vector holds output from the LHS (implicit) portion of a term at a
  // given vertical level.  This output is weighted and applied to the RHS.
  // This is used if the implicit portion of the term is "over-implicit", which
  // means that the LHS contribution is given extra weight (>1) in order to
  // increase numerical stability.  A weighted factor must then be applied to
  // the RHS in order to balance the weight.
  //------------------- Output Variable -------------------
  //------------------- Local Variables -------------------
  // Indices
  auto workspace_scope = workspace.scope();
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rhs_bp_pr3 =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_rhs::rhs_bp_pr3", ngrdcol,
          nzm); // Buoyancy production of w'x' and w'x' pressure term 3
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rhs_bp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_rhs::rhs_bp", ngrdcol,
          nzm); // Buoyancy production of w'x' (stats only)
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rhs_pr3 =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_rhs::rhs_pr3", ngrdcol,
          nzm); // w'x' pressure term 3 (stats only)
  real_t invrs_dt;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> C7_Skw_fnc_plus_one =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_rhs::C7_Skw_fnc_plus_one", ngrdcol,
          nzm); // stats only
  Kokkos::View<real_t**, Kokkos::LayoutLeft> C7_Skw_fnc_zeros =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_rhs::C7_Skw_fnc_zeros", ngrdcol, nzm); // stats only
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wpxp_forcing_stats =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_rhs::wpxp_forcing_stats", ngrdcol,
          nzm); // stats only
  Kokkos::View<real_t**, Kokkos::LayoutLeft> ta_over =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_rhs::ta_over", ngrdcol, nzm); // stats only
  Kokkos::View<real_t**, Kokkos::LayoutLeft> pr1_over =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_rhs::pr1_over", ngrdcol, nzm); // stats only
  int rhs_lb_idx_zm; // Index of (momentum) lower boundary level in rhs array
  int rhs_ub_idx_zm; // Index of (momentum) upper boundary level in rhs array
  std::string name_xm_f;
  std::string name_wpxp_bp;
  std::string name_wpxp_pr3;
  std::string name_wpxp_f;
  std::string name_wpxp_ta;
  std::string name_wpxp_pr1;
  const auto gr_grid_dir_indx = gr.grid_dir_indx;
  const auto gr_k_lb_zm       = gr.k_lb_zm;
  //------------------- Begin Code -------------------
  // Initialize output array and precalculate the reciprocal of dt
  invrs_dt = real_t(1.0) / dt;
  // Calculate buoyancy production of w'x' and w'x' pressure term 3
  wpxp_terms_bp_pr3_rhs(nzm, ngrdcol, gr, C7_Skw_fnc, thv_ds_zm, xpthvp, // In
                        rhs_bp_pr3                                       // Out
  );
  // Index of the momentum boundary levels
  if (gr.grid_dir_indx > 0) {
    // Ascending Grid
    rhs_lb_idx_zm = 0;
    rhs_ub_idx_zm = 2 * nzm - 2;
  } else {
    // Descending Grid
    rhs_lb_idx_zm = 2 * nzm - 2;
    rhs_ub_idx_zm = 0;
  }
  Kokkos::parallel_for(
      "xm_wpxp_rhs::loop1", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        const int i = cell;
        // Set lower boundary for w'x'
        rhs(i, rhs_lb_idx_zm) = wpxp(i, gr_k_lb_zm);
      });
  // Combine terms to calculate other values, rhs(3) to rhs(gr%nz-2)
  Kokkos::parallel_for(
      "xm_wpxp_rhs::loop2", std::max(0, nzt) * std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        const int i = (cell % ngrdcol);
        const int k = (cell / ngrdcol);
        int k_xm;
        k_xm = 2 * k + 1;
        // RHS time tendency and forcings for xm
        // Note: xm forcings include the effects of microphysics,
        //       cloud water sedimentation, radiation, and any
        //       imposed forcings on xm.
        rhs(i, k_xm) = xm(i, k) * invrs_dt + xm_forcing(i, k);
      });
  Kokkos::parallel_for(
      "xm_wpxp_rhs::loop3", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i = (cell % ngrdcol);
        const int k = (cell / ngrdcol) + 1;
        int k_wpxp;
        k_wpxp = 2 * k;
        // Calculate rhs values for w'x' using precalculated terms
        rhs(i, k_wpxp) =
            rhs_bp_pr3(i, k) + wpxp_forcing(i, k) + rhs_ta(i, k) +
            (one - gamma_over_implicit_ts) *
                (-lhs_ta_wpxp(-gr_grid_dir_indx + 1, i, k) * wpxp(i, gr_grid_dir_indx + k) -
                 lhs_ta_wpxp(1, i, k) * wpxp(i, k) -
                 lhs_ta_wpxp(gr_grid_dir_indx + 1, i, k) * wpxp(i, -gr_grid_dir_indx + k) -
                 lhs_pr1(i, k) * wpxp(i, k));
      });
  Kokkos::parallel_for(
      "xm_wpxp_rhs::loop4", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        const int i = cell;
        // Upper boundary for w'x'
        rhs(i, rhs_ub_idx_zm) = real_t(0.0);
      });
  // RHS time tendency.
  if (l_iter) {
    Kokkos::parallel_for(
        "xm_wpxp_rhs::loop5", std::max(0, nzm - 2) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol) + 1;
          int k_wpxp;
          k_wpxp         = 2 * k;
          rhs(i, k_wpxp) = rhs(i, k_wpxp) + wpxp(i, k) * invrs_dt;
        });
  }
  if (stats.l_sample) {
    switch (solve_type) {
    case xm_wpxp_rtm: {
      name_xm_f     = "rtm_forcing";
      name_wpxp_bp  = "wprtp_bp";
      name_wpxp_pr3 = "wprtp_pr3";
      name_wpxp_f   = "wprtp_forcing";
      name_wpxp_ta  = "wprtp_ta";
      name_wpxp_pr1 = "wprtp_pr1";
      break;
    }
    case xm_wpxp_thlm: {
      name_xm_f     = "thlm_forcing";
      name_wpxp_bp  = "wpthlp_bp";
      name_wpxp_pr3 = "wpthlp_pr3";
      name_wpxp_f   = "wpthlp_forcing";
      name_wpxp_ta  = "wpthlp_ta";
      name_wpxp_pr1 = "wpthlp_pr1";
      break;
    }
    case xm_wpxp_um: {
      name_xm_f     = "";
      name_wpxp_bp  = "upwp_bp";
      name_wpxp_pr3 = "upwp_pr3";
      name_wpxp_f   = "";
      name_wpxp_ta  = "upwp_ta";
      name_wpxp_pr1 = "upwp_pr1";
      break;
    }
    case xm_wpxp_vm: {
      name_xm_f     = "";
      name_wpxp_bp  = "vpwp_bp";
      name_wpxp_pr3 = "vpwp_pr3";
      name_wpxp_f   = "";
      name_wpxp_ta  = "vpwp_ta";
      name_wpxp_pr1 = "vpwp_pr1";
      break;
    }
    default: {
      name_xm_f     = "";
      name_wpxp_bp  = "";
      name_wpxp_pr3 = "";
      name_wpxp_f   = "";
      name_wpxp_ta  = "";
      name_wpxp_pr1 = "";
      break;
    }
    }
    Kokkos::parallel_for(
        "xm_wpxp_rhs::array6", C7_Skw_fnc_zeros.extent(0) * C7_Skw_fnc_zeros.extent(1),
        KOKKOS_LAMBDA(int cell) {
          const int slice_0                  = (cell % C7_Skw_fnc_zeros.extent(0));
          const int slice_1                  = (cell / C7_Skw_fnc_zeros.extent(0));
          C7_Skw_fnc_zeros(slice_0, slice_1) = zero;
        });
    Kokkos::parallel_for(
        "xm_wpxp_rhs::array7", C7_Skw_fnc_plus_one.extent(0) * C7_Skw_fnc_plus_one.extent(1),
        KOKKOS_LAMBDA(int cell) {
          const int slice_0                     = (cell % C7_Skw_fnc_plus_one.extent(0));
          const int slice_1                     = (cell / C7_Skw_fnc_plus_one.extent(0));
          C7_Skw_fnc_plus_one(slice_0, slice_1) = C7_Skw_fnc(slice_0, slice_1) + one;
        });
    // Statistics: explicit contributions for wpxp.
    // w'x' term bp is completely explicit; call stat_update_var.
    // Note:  To find the contribution of w'x' term bp, substitute 0 for the
    //        C_7 skewness function input to function wpxp_terms_bp_pr3_rhs.
    wpxp_terms_bp_pr3_rhs(nzm, ngrdcol, gr, C7_Skw_fnc_zeros, // In
                          thv_ds_zm, xpthvp,                  // In
                          rhs_bp                              // Out
    );
    // w'x' term pr3 is completely explicit; call stat_update_var.
    // Note:  To find the contribution of w'x' term pr3, add 1 to the
    //        C_7 skewness function input to function wpxp_terms_bp_pr2_rhs.
    wpxp_terms_bp_pr3_rhs(nzm, ngrdcol, gr, C7_Skw_fnc_plus_one, // In
                          thv_ds_zm, xpthvp,                     // In
                          rhs_pr3                                // Out
    );
    // Keep the bp/pr3 split consistent with the legacy stats decomposition.
    // rhs_bp uses C7_Skw_fnc=0 and rhs_pr3 uses C7_Skw_fnc+1.
    stats_update(name_wpxp_bp, rhs_bp, // In
                 stats                 // InOut
    );
    stats_update(name_wpxp_pr3, rhs_pr3, // In
                 stats                   // InOut
    );
    if (name_wpxp_f.size() > 0) {
      Kokkos::parallel_for(
          "xm_wpxp_rhs::array8", wpxp_forcing_stats.extent(0), KOKKOS_LAMBDA(int cell) {
            const int slice_0              = cell;
            wpxp_forcing_stats(slice_0, 0) = real_t(0.0);
          });
      Kokkos::parallel_for(
          "xm_wpxp_rhs::array9", wpxp_forcing_stats.extent(0), KOKKOS_LAMBDA(int cell) {
            const int slice_0                      = cell;
            wpxp_forcing_stats(slice_0, (nzm - 1)) = real_t(0.0);
          });
      Kokkos::parallel_for(
          "xm_wpxp_rhs::array10",
          wpxp_forcing_stats.extent(0) * nzm - 2 * wpxp_forcing_stats.extent(0),
          KOKKOS_LAMBDA(int cell) {
            const int slice_0                          = (cell % wpxp_forcing_stats.extent(0));
            const int slice_1                          = (cell / wpxp_forcing_stats.extent(0));
            wpxp_forcing_stats(slice_0, (slice_1 + 1)) = wpxp_forcing(slice_0, (slice_1 + 1));
          });
      // w'x' forcing term is completely explicit; call stat_update_var_pt.
      stats_update(name_wpxp_f, wpxp_forcing_stats, // In
                   stats                            // InOut
      );
    }
    if (name_wpxp_ta.size() > 0) {
      // <w'x'> term ta has both implicit and explicit components; call
      // stat_begin_update_pt.  Since stat_begin_update_pt automatically
      // subtracts the value sent in, reverse the sign on
      // xpyp_term_ta_pdf_rhs.
      // Start with the explicit rhs term via stats_begin_budget.
      stats_begin_budget(name_wpxp_ta, rhs_ta, // In
                         stats                 // InOut
      );
      Kokkos::parallel_for(
          "xm_wpxp_rhs::array11", ta_over.extent(0), KOKKOS_LAMBDA(int cell) {
            const int slice_0   = cell;
            ta_over(slice_0, 0) = real_t(0.0);
          });
      Kokkos::parallel_for(
          "xm_wpxp_rhs::loop12", std::max(0, nzm - 2) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol) + 1;
            ta_over(i, k) =
                (one - gamma_over_implicit_ts) *
                (-lhs_ta_wpxp(-gr_grid_dir_indx + 1, i, k) * wpxp(i, gr_grid_dir_indx + k) -
                 lhs_ta_wpxp(1, i, k) * wpxp(i, k) -
                 lhs_ta_wpxp(gr_grid_dir_indx + 1, i, k) * wpxp(i, -gr_grid_dir_indx + k));
          });
      Kokkos::parallel_for(
          "xm_wpxp_rhs::array13", ta_over.extent(0), KOKKOS_LAMBDA(int cell) {
            const int slice_0           = cell;
            ta_over(slice_0, (nzm - 1)) = real_t(0.0);
          });
      // Note:  An "over-implicit" weighted time step is applied to this term.
      //        A weighting factor of greater than 1 may be used to make the
      //        term more numerically stable (see note above for RHS
      //        contribution from "over-implicit" weighted time step for LHS
      //        turbulent advection (ta) term).
      stats_update_budget(name_wpxp_ta, ta_over, // In
                          stats                  // InOut
      );
    }
    if (name_wpxp_pr1.size() > 0) {
      // w'x' term pr1 is normally completely implicit.  However, there is a
      // RHS contribution from the "over-implicit" weighted time step.
      // RHS contribution from the "over-implicit" weighted time step.  A
      // weighting factor of greater than 1 may be used to make the term more
      // numerically stable (see note above for RHS contribution from
      // "over-implicit" weighted time step for LHS turbulent advection (ta)
      // term).  Therefore, w'x' term pr1 has both implicit and explicit
      // components; call stat_begin_update_pt.  Since stat_begin_update_pt
      // Therefore, w'x' term pr1 has both implicit and explicit components;
      // call stat_begin_update_pt.  Since stat_begin_update_pt
      // automatically subtracts the value sent in, reverse the sign on the
      // input value.
      Kokkos::parallel_for(
          "xm_wpxp_rhs::array14", pr1_over.extent(0), KOKKOS_LAMBDA(int cell) {
            const int slice_0    = cell;
            pr1_over(slice_0, 0) = real_t(0.0);
          });
      Kokkos::parallel_for(
          "xm_wpxp_rhs::loop15", std::max(0, nzm - 2) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i    = (cell % ngrdcol);
            const int k    = (cell / ngrdcol) + 1;
            pr1_over(i, k) = (one - gamma_over_implicit_ts) * lhs_pr1(i, k) * wpxp(i, k);
          });
      Kokkos::parallel_for(
          "xm_wpxp_rhs::array16", pr1_over.extent(0), KOKKOS_LAMBDA(int cell) {
            const int slice_0            = cell;
            pr1_over(slice_0, (nzm - 1)) = real_t(0.0);
          });
      stats_begin_budget(name_wpxp_pr1, pr1_over, // In
                         stats                    // InOut
      );
    }
    // Statistics: explicit contributions for xm
    //             (including microphysics/radiation).
    // xm forcings term is completely explicit; call stat_update_var_pt.
    // xm forcing (e.g., large-scale forcing and source terms) is explicit.
    stats_update(name_xm_f, xm_forcing, // In
                 stats                  // InOut
    );
  }
  return;
}

// pdf_implicit_coefs_terms: Implicit coefs / explicit terms [units vary]
// iiPDF_type: Selected option for the two-component normal (double Gaussian) PDF type to use
// for the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate,
// two-component PDF.
// l_godunov_upwind_wpxp_ta: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent advection terms. It affects
// wpxp only.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17>
void calc_xm_wpxp_ta_terms(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, int sclr_dim,
                           const grid& gr,
                           Array0 wp2rtp,                                                  // In
                           Array1 wp2thlp, Array2 wp2sclrp,                                // In
                           Array3 rho_ds_zt, Array4 invrs_rho_ds_zm, Array5 rho_ds_zm,     // In
                           Array6 sigma_sqd_w, Array7 wp3_on_wp2_zt,                       // In
                           const implicit_coefs_terms& pdf_implicit_coefs_terms,           // In
                           int iiPDF_type,                                                 // In
                           bool l_explicit_turbulent_adv_wpxp, bool l_predict_upwp_vpwp,   // In
                           bool l_scalar_calc,                                             // In
                           bool l_godunov_upwind_wpxp_ta,                                  // In
                           stats_type& stats,                                              // InOut
                           Array8 lhs_ta_wprtp, Array9 lhs_ta_wpthlp, Array10 lhs_ta_wpup, // Out
                           Array11 lhs_ta_wpvp, Array12 lhs_ta_wpsclrp,                    // Out
                           Array13 rhs_ta_wprtp, Array14 rhs_ta_wpthlp, Array15 rhs_ta_wpup, // Out
                           Array16 rhs_ta_wpvp, Array17 rhs_ta_wpsclrp                       // Out
) {
  //
  // Description: This subroutine calculates the turbulent advection terms for
  //              the left and right hand side matrices. Solutions may be entirely
  //              explicit, entirely implicit, or mixed between, depending on
  //              various flags and the PDF type.
  //---------------------------------------------------------------------------------------------
  //------------------- Input Variables -------------------
  // Gaussian) PDF type to use for the w, rt, and theta-l (or
  // w, chi, and eta) portion of CLUBB's multivariate,
  // two-component PDF.
  // differencing approximation rather than a centered
  // differencing for turbulent advection terms.
  // It affects  wpxp only.
  // passed into the xpyp_term_ta_pdf_rhs subroutine.
  // This stems from removing the l_upwind_wpxp_ta flag.
  // More information on this can be found on issue #926
  // on the clubb repository.
  //------------------- Output Variables -------------------
  //------------------- Local Variables -------------------
  // Variables for turbulent advection of predictive variances and covariances.
  // <w'^2 rt'> = coef_wp2rtp_implicit * <w'rt'> + term_wp2rtp_explicit
  // <w'^2 thl'> = coef_wp2thlp_implicit * <w'thl'> + term_wp2thlp_explicit
  // <w'^2 sclr'> = coef_wp2sclrp_implicit * <w'sclr'> + term_wp2sclrp_explicit
  // Sign of turbulent velocity (used for "upwind" turbulent advection)
  auto workspace_scope         = workspace.scope();
  constexpr bool l_dummy_false = false;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wp2rtp_implicit =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_ta_terms::coef_wp2rtp_implicit", ngrdcol,
          nzt); // Coefficient that is multiplied by <w'rt'>  [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wp2rtp_explicit =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_ta_terms::term_wp2rtp_explicit", ngrdcol,
          nzt); // Term that is on the RHS          [m^2/s^2 kg/kg]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wp2rtp_implicit_zm =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_ta_terms::coef_wp2rtp_implicit_zm", ngrdcol,
          nzm); // coef_wp2rtp_implicit interp. to m-levs. [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wp2rtp_explicit_zm =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_ta_terms::term_wp2rtp_explicit_zm", ngrdcol,
          nzm); // term_wp2rtp_expl intrp m-levs [m^2/s^2 kg/kg]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wp2thlp_implicit =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_ta_terms::coef_wp2thlp_implicit", ngrdcol,
          nzt); // Coef. that is multiplied by <w'thl'>      [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wp2thlp_explicit =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_ta_terms::term_wp2thlp_explicit", ngrdcol,
          nzt); // Term that is on the RHS             [m^2/s^2 K]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wp2thlp_implicit_zm =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_ta_terms::coef_wp2thlp_implicit_zm", ngrdcol,
          nzm); // coef_wp2thlp_implicit interp. m-levs.  [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wp2thlp_explicit_zm =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_ta_terms::term_wp2thlp_explicit_zm", ngrdcol,
          nzm); // term_wp2thlp_expl interp. m-levs [m^2/s^2 K]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wp2sclrp_explicit =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_ta_terms::term_wp2sclrp_explicit", ngrdcol,
          nzt); // Term that is on the RHS    [m^2/s^2(un. vary)]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wp2sclrp_explicit_zm =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_ta_terms::term_wp2sclrp_explicit_zm", ngrdcol,
          nzm); // term_wp2sclrp_expl intrp zm [m^2/s^2(un v)]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sgn_t_vel_wprtp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_ta_terms::sgn_t_vel_wprtp", ngrdcol,
          nzm); // Sign of the turbulent velocity for <w'rt'>       [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sgn_t_vel_wpthlp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_ta_terms::sgn_t_vel_wpthlp", ngrdcol,
          nzm); // Sign of the turbulent velocity for <w'thl'>      [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sgn_t_vel_wpsclrp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_ta_terms::sgn_t_vel_wpsclrp", ngrdcol,
          nzm); // Sign of the turbulent velocity for <w'sclr'>    [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> a1_coef =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_ta_terms::a1_coef", ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> a1_coef_zt =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "calc_xm_wpxp_ta_terms::a1_coef_zt", ngrdcol, nzt);
  int i;
  int sclr;
  const auto pdf_implicit_coefs_terms_coef_wp2rtp_implicit =
      pdf_implicit_coefs_terms.coef_wp2rtp_implicit;
  const auto pdf_implicit_coefs_terms_coef_wp2thlp_implicit =
      pdf_implicit_coefs_terms.coef_wp2thlp_implicit;
  const auto pdf_implicit_coefs_terms_term_wp2rtp_explicit =
      pdf_implicit_coefs_terms.term_wp2rtp_explicit;
  const auto pdf_implicit_coefs_terms_term_wp2thlp_explicit =
      pdf_implicit_coefs_terms.term_wp2thlp_explicit;
  //------------------- Begin Code -------------------
  // Set up the implicit coefficients and explicit terms for turbulent
  // advection of <w'rt'>, <w'thl'>, and <w'sclr'>.
  if (l_explicit_turbulent_adv_wpxp) {
    // The turbulent advection of <w'x'> is handled explicitly
    // The turbulent advection of <w'x'> is handled explicitly, the
    // terms are calculated only for the RHS matrices. The
    // term_wp2xp_explicit terms are equal to <w'x'> as calculated using PDF
    // parameters, which are general for any PDF type. The values of
    // <w'x'> are calculated on thermodynamic levels.
    // These coefficients only need to be set if stats output is on
    if (stats.l_sample) {
      Kokkos::parallel_for(
          "calc_xm_wpxp_ta_terms::array1",
          coef_wp2rtp_implicit.extent(0) * coef_wp2rtp_implicit.extent(1), KOKKOS_LAMBDA(int cell) {
            const int slice_0                      = (cell % coef_wp2rtp_implicit.extent(0));
            const int slice_1                      = (cell / coef_wp2rtp_implicit.extent(0));
            coef_wp2rtp_implicit(slice_0, slice_1) = zero;
          });
      Kokkos::parallel_for(
          "calc_xm_wpxp_ta_terms::array2",
          coef_wp2thlp_implicit.extent(0) * coef_wp2thlp_implicit.extent(1),
          KOKKOS_LAMBDA(int cell) {
            const int slice_0                       = (cell % coef_wp2thlp_implicit.extent(0));
            const int slice_1                       = (cell / coef_wp2thlp_implicit.extent(0));
            coef_wp2thlp_implicit(slice_0, slice_1) = zero;
          });
    }
    // The turbulent advection terms are handled entirely explicitly. Thus the LHS
    // terms can be set to zero.
    Kokkos::parallel_for(
        "calc_xm_wpxp_ta_terms::loop3",
        std::max(0, nzm) * std::max(0, ngrdcol) * std::max(0, ndiags3),
        KOKKOS_LAMBDA(const int cell) {
          const int b            = (cell % ndiags3);
          const int i            = ((cell / ndiags3) % ngrdcol);
          const int k            = (cell / (ndiags3 * ngrdcol));
          lhs_ta_wprtp(b, i, k)  = zero;
          lhs_ta_wpthlp(b, i, k) = zero;
        });
    if (l_scalar_calc) {
      Kokkos::parallel_for(
          "calc_xm_wpxp_ta_terms::loop4",
          std::max(0, sclr_dim) * std::max(0, nzm) * std::max(0, ngrdcol) * std::max(0, ndiags3),
          KOKKOS_LAMBDA(const int cell) {
            const int b                           = (cell % ndiags3);
            const int i                           = ((cell / ndiags3) % ngrdcol);
            const int k                           = ((cell / (ndiags3 * ngrdcol)) % nzm);
            [[maybe_unused]] const int sclr       = (cell / (ndiags3 * ngrdcol * nzm));
            lhs_ta_wpsclrp(b, i, k, sclr_dim - 1) = zero;
          });
    }
    Kokkos::parallel_for(
        "calc_xm_wpxp_ta_terms::loop5", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i                 = (cell % ngrdcol);
          const int k                 = (cell / ngrdcol);
          term_wp2rtp_explicit(i, k)  = wp2rtp(i, k);
          term_wp2thlp_explicit(i, k) = wp2thlp(i, k);
        });
    // Calculate the RHS turbulent advection term for <w'r_t'>
    xpyp_term_ta_pdf_rhs(nzm, nzt, ngrdcol, gr, term_wp2rtp_explicit, // In
                         rho_ds_zt, rho_ds_zm,                        // In
                         invrs_rho_ds_zm,                             // In
                         l_dummy_false,                               // In
                         sgn_t_vel_wprtp,                             // In
                         term_wp2rtp_explicit_zm,                     // In
                         rhs_ta_wprtp                                 // Out
    );
    // Calculate the RHS turbulent advection term for <w'thl'>
    xpyp_term_ta_pdf_rhs(nzm, nzt, ngrdcol, gr, term_wp2thlp_explicit, // In
                         rho_ds_zt, rho_ds_zm,                         // In
                         invrs_rho_ds_zm,                              // In
                         l_dummy_false,                                // In
                         sgn_t_vel_wpthlp,                             // In
                         term_wp2thlp_explicit_zm,                     // In
                         rhs_ta_wpthlp                                 // Out
    );
    for (sclr = 0; sclr < sclr_dim; sclr += 1) {
      Kokkos::parallel_for(
          "calc_xm_wpxp_ta_terms::loop6", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i                  = (cell % ngrdcol);
            const int k                  = (cell / ngrdcol);
            term_wp2sclrp_explicit(i, k) = wp2sclrp(i, k, sclr);
          });
      // Calculate the RHS turbulent advection term for <w'thl'>
      { // Native sequence association: this source call writes the first
        // ngrdcol*nzm elements of the rank-3 LHS scalar section, not RHS.
        // Preserve that source mapping; this disabled branch leaves scalar RHS undefined.
        auto lhs_scalar =
            Kokkos::subview(lhs_ta_wpsclrp, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, sclr);
        Kokkos::View<real_t**, Kokkos::LayoutLeft, typename decltype(lhs_scalar)::memory_space,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            lhs_sequence(lhs_scalar.data(), ngrdcol, nzm);
        xpyp_term_ta_pdf_rhs(nzm, nzt, ngrdcol, gr, term_wp2sclrp_explicit, // In
                             rho_ds_zt, rho_ds_zm,                          // In
                             invrs_rho_ds_zm,                               // In
                             l_dummy_false,                                 // In
                             sgn_t_vel_wpsclrp,                             // In
                             term_wp2sclrp_explicit_zm,                     // In
                             lhs_sequence                                   // Out
        );
      }
    }
  } else {
    // The turbulent advection of <w'x'> is handled implicitly or
    // semi-implicitly.
    if (iiPDF_type == iiPDF_ADG1) {
      // The ADG1 PDF is used.
      // Calculate the implicit coefficients and explicit terms on
      // thermodynamic grid levels.
      // Calculate a_1.
      // It is a variable that is a function of sigma_sqd_w (where
      // sigma_sqd_w is located on momentum levels).
      Kokkos::parallel_for(
          "calc_xm_wpxp_ta_terms::loop7", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i   = (cell % ngrdcol);
            const int k   = (cell / ngrdcol);
            a1_coef(i, k) = one / (one - sigma_sqd_w(i, k));
          });
      // Interpolate a_1 from momentum levels to thermodynamic levels.  This
      // will be used for the <w'x'> turbulent advection (ta) term.
      zm2zt_api(nzm, nzt, ngrdcol, gr, a1_coef,
                Kokkos::subview(a1_coef_zt, Kokkos::ALL, Kokkos::ALL), zero_threshold);
      Kokkos::parallel_for(
          "calc_xm_wpxp_ta_terms::loop8", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i                 = (cell % ngrdcol);
            const int k                 = (cell / ngrdcol);
            coef_wp2rtp_implicit(i, k)  = a1_coef_zt(i, k) * wp3_on_wp2_zt(i, k);
            coef_wp2thlp_implicit(i, k) = coef_wp2rtp_implicit(i, k);
          });
      if (!l_godunov_upwind_wpxp_ta) {
        // Calculate the LHS turbulent advection term for <w'r_t'>
        xpyp_term_ta_pdf_lhs(nzm, nzt, ngrdcol, gr, coef_wp2rtp_implicit, // In
                             rho_ds_zt, rho_ds_zm,                        // In
                             invrs_rho_ds_zm,                             // In
                             l_dummy_false,                               // In
                             sgn_t_vel_wprtp,                             // In
                             coef_wp2rtp_implicit_zm,                     // In
                             lhs_ta_wprtp                                 // Out
        );
      } else {
        // Godunov-like method for the vertical discretization of ta term
        Kokkos::parallel_for(
            "calc_xm_wpxp_ta_terms::loop9", std::max(0, nzt) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              const int i                 = (cell % ngrdcol);
              const int k                 = (cell / ngrdcol);
              coef_wp2rtp_implicit(i, k)  = a1_coef_zt(i, k) * wp3_on_wp2_zt(i, k);
              coef_wp2thlp_implicit(i, k) = coef_wp2rtp_implicit(i, k);
            });
        xpyp_term_ta_pdf_lhs_godunov(nzm, nzt, ngrdcol, gr,      // In
                                     coef_wp2rtp_implicit,       // In
                                     invrs_rho_ds_zm, rho_ds_zm, // In
                                     lhs_ta_wprtp                // Out
        );
      }
      // For ADG1, the LHS turbulent advection terms for
      // <w'r_t'>, <w'thl'>, <w'sclr'> are all equal
      Kokkos::parallel_for(
          "calc_xm_wpxp_ta_terms::loop10",
          std::max(0, nzm) * std::max(0, ngrdcol) * std::max(0, ndiags3),
          KOKKOS_LAMBDA(const int cell) {
            const int b            = (cell % ndiags3);
            const int i            = ((cell / ndiags3) % ngrdcol);
            const int k            = (cell / (ndiags3 * ngrdcol));
            lhs_ta_wpthlp(b, i, k) = lhs_ta_wprtp(b, i, k);
          });
      if (l_scalar_calc) {
        Kokkos::parallel_for(
            "calc_xm_wpxp_ta_terms::loop11",
            std::max(0, sclr_dim) * std::max(0, nzm) * std::max(0, ngrdcol) * std::max(0, ndiags3),
            KOKKOS_LAMBDA(const int cell) {
              const int b                   = (cell % ndiags3);
              const int i                   = ((cell / ndiags3) % ngrdcol);
              const int k                   = ((cell / (ndiags3 * ngrdcol)) % nzm);
              const int sclr                = (cell / (ndiags3 * ngrdcol * nzm));
              lhs_ta_wpsclrp(b, i, k, sclr) = lhs_ta_wprtp(b, i, k);
            });
      }
      if (stats.l_sample) {
        Kokkos::parallel_for(
            "calc_xm_wpxp_ta_terms::loop12", std::max(0, nzt) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              const int i                 = (cell % ngrdcol);
              const int k                 = (cell / ngrdcol);
              term_wp2rtp_explicit(i, k)  = zero;
              term_wp2thlp_explicit(i, k) = zero;
            });
      }
      // The <w'r_t'>, <w'thl'>, <w'sclr'> turbulent advection terms are entirely implicit.
      // Set the RHS turbulent advection terms to 0
      Kokkos::parallel_for(
          "calc_xm_wpxp_ta_terms::loop13", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i         = (cell % ngrdcol);
            const int k         = (cell / ngrdcol);
            rhs_ta_wprtp(i, k)  = zero;
            rhs_ta_wpthlp(i, k) = zero;
          });
      if (l_scalar_calc) {
        Kokkos::parallel_for(
            "calc_xm_wpxp_ta_terms::loop14",
            std::max(0, sclr_dim) * std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              const int i                = (cell % ngrdcol);
              const int k                = ((cell / ngrdcol) % nzm);
              const int sclr             = (cell / (ngrdcol * nzm));
              rhs_ta_wpsclrp(i, k, sclr) = zero;
            });
      }
      if (l_predict_upwp_vpwp) {
        // Predict <u> and <u'w'>, as well as <v> and <v'w'>.
        // These terms are equal to the <w'r_t'> terms as well in this case
        Kokkos::parallel_for(
            "calc_xm_wpxp_ta_terms::loop15",
            std::max(0, nzm) * std::max(0, ngrdcol) * std::max(0, ndiags3),
            KOKKOS_LAMBDA(const int cell) {
              const int b          = (cell % ndiags3);
              const int i          = ((cell / ndiags3) % ngrdcol);
              const int k          = (cell / (ndiags3 * ngrdcol));
              lhs_ta_wpup(b, i, k) = lhs_ta_wprtp(b, i, k);
              lhs_ta_wpvp(b, i, k) = lhs_ta_wprtp(b, i, k);
            });
        // The <w'u'> and <w'v'> turbulent advection terms are entirely implicit.
        // Set the RHS turbulent advection terms to 0
        Kokkos::parallel_for(
            "calc_xm_wpxp_ta_terms::loop16", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              const int i       = (cell % ngrdcol);
              const int k       = (cell / ngrdcol);
              rhs_ta_wpup(i, k) = zero;
              rhs_ta_wpvp(i, k) = zero;
            });
      }
    } else if (iiPDF_type == iiPDF_new) {
      // The new PDF is used.
      // Unpack the variables coef_wp2rtp_implicit, term_wp2rtp_explicit,
      // coef_wp2thlp_implicit, and term_wp2thlp_explicit from
      // pdf_implicit_coefs_terms.  The PDF parameters and the resulting
      // implicit coefficients and explicit terms are calculated on
      // thermodynamic levels.
      for (i = 0; i < ngrdcol; i += 1) {
        Kokkos::parallel_for(
            "calc_xm_wpxp_ta_terms::array17", coef_wp2rtp_implicit.extent(1),
            KOKKOS_LAMBDA(int cell) {
              const int slice_0 = cell;
              coef_wp2rtp_implicit(i, slice_0) =
                  pdf_implicit_coefs_terms_coef_wp2rtp_implicit(i, slice_0);
            });
        Kokkos::parallel_for(
            "calc_xm_wpxp_ta_terms::array18", coef_wp2thlp_implicit.extent(1),
            KOKKOS_LAMBDA(int cell) {
              const int slice_0 = cell;
              coef_wp2thlp_implicit(i, slice_0) =
                  pdf_implicit_coefs_terms_coef_wp2thlp_implicit(i, slice_0);
            });
        Kokkos::parallel_for(
            "calc_xm_wpxp_ta_terms::array19", term_wp2rtp_explicit.extent(1),
            KOKKOS_LAMBDA(int cell) {
              const int slice_0 = cell;
              term_wp2rtp_explicit(i, slice_0) =
                  pdf_implicit_coefs_terms_term_wp2rtp_explicit(i, slice_0);
            });
        Kokkos::parallel_for(
            "calc_xm_wpxp_ta_terms::array20", term_wp2thlp_explicit.extent(1),
            KOKKOS_LAMBDA(int cell) {
              const int slice_0 = cell;
              term_wp2thlp_explicit(i, slice_0) =
                  pdf_implicit_coefs_terms_term_wp2thlp_explicit(i, slice_0);
            });
      }
      // Calculate the LHS turbulent advection term for <w'rt'>
      xpyp_term_ta_pdf_lhs(nzm, nzt, ngrdcol, gr, coef_wp2rtp_implicit, // In
                           rho_ds_zt, rho_ds_zm,                        // In
                           invrs_rho_ds_zm,                             // In
                           l_dummy_false,                               // In
                           sgn_t_vel_wprtp,                             // In
                           coef_wp2rtp_implicit_zm,                     // In
                           lhs_ta_wprtp                                 // Out
      );
      // Calculate the RHS turbulent advection term for <w'rt'>
      xpyp_term_ta_pdf_rhs(nzm, nzt, ngrdcol, gr, term_wp2rtp_explicit, // In
                           rho_ds_zt, rho_ds_zm,                        // In
                           invrs_rho_ds_zm,                             // In
                           l_dummy_false,                               // In
                           sgn_t_vel_wprtp,                             // In
                           term_wp2rtp_explicit_zm,                     // In
                           rhs_ta_wprtp                                 // Out
      );
      // Calculate the LHS turbulent advection term for <w'thl'>
      xpyp_term_ta_pdf_lhs(nzm, nzt, ngrdcol, gr, coef_wp2thlp_implicit, // In
                           rho_ds_zt, rho_ds_zm,                         // In
                           invrs_rho_ds_zm,                              // In
                           l_dummy_false,                                // In
                           sgn_t_vel_wpthlp,                             // In
                           coef_wp2thlp_implicit_zm,                     // In
                           lhs_ta_wpthlp                                 // Out
      );
      // Calculate the RHS turbulent advection term for <w'thl'>
      xpyp_term_ta_pdf_rhs(nzm, nzt, ngrdcol, gr, term_wp2thlp_explicit, // In
                           rho_ds_zt, rho_ds_zm,                         // In
                           invrs_rho_ds_zm,                              // In
                           l_dummy_false,                                // In
                           sgn_t_vel_wpthlp,                             // In
                           term_wp2thlp_explicit_zm,                     // In
                           rhs_ta_wpthlp                                 // Out
      );
      // The code for the scalar variables will be set up later.
      Kokkos::parallel_for(
          "calc_xm_wpxp_ta_terms::array21",
          ((lhs_ta_wpsclrp.extent(0) - 1) / 1 + 1) * ((lhs_ta_wpsclrp.extent(1) - 1) / 1 + 1) *
              ((lhs_ta_wpsclrp.extent(2) - 1) / 1 + 1) * ((lhs_ta_wpsclrp.extent(3) - 1) / 1 + 1),
          KOKKOS_LAMBDA(int cell) {
            const int slice_0 = cell % ((lhs_ta_wpsclrp.extent(0) - 1) / 1 + 1);
            const int slice_1 = (cell / (((lhs_ta_wpsclrp.extent(0) - 1) / 1 + 1))) %
                                ((lhs_ta_wpsclrp.extent(1) - 1) / 1 + 1);
            const int slice_2 = (cell / (((lhs_ta_wpsclrp.extent(0) - 1) / 1 + 1) *
                                         ((lhs_ta_wpsclrp.extent(1) - 1) / 1 + 1))) %
                                ((lhs_ta_wpsclrp.extent(2) - 1) / 1 + 1);
            const int slice_3 = cell / (((lhs_ta_wpsclrp.extent(0) - 1) / 1 + 1) *
                                        ((lhs_ta_wpsclrp.extent(1) - 1) / 1 + 1) *
                                        ((lhs_ta_wpsclrp.extent(2) - 1) / 1 + 1));
            lhs_ta_wpsclrp(slice_0, slice_1, slice_2, slice_3) = zero;
          });
      Kokkos::parallel_for(
          "calc_xm_wpxp_ta_terms::array22",
          ((rhs_ta_wpsclrp.extent(0) - 1) / 1 + 1) * ((rhs_ta_wpsclrp.extent(1) - 1) / 1 + 1) *
              ((rhs_ta_wpsclrp.extent(2) - 1) / 1 + 1),
          KOKKOS_LAMBDA(int cell) {
            const int slice_0 = cell % ((rhs_ta_wpsclrp.extent(0) - 1) / 1 + 1);
            const int slice_1 = (cell / (((rhs_ta_wpsclrp.extent(0) - 1) / 1 + 1))) %
                                ((rhs_ta_wpsclrp.extent(1) - 1) / 1 + 1);
            const int slice_2 = cell / (((rhs_ta_wpsclrp.extent(0) - 1) / 1 + 1) *
                                        ((rhs_ta_wpsclrp.extent(1) - 1) / 1 + 1));
            rhs_ta_wpsclrp(slice_0, slice_1, slice_2) = zero;
          });
    } else if (iiPDF_type == iiPDF_new_hybrid) {
      // The new hybrid PDF is used.
      // Unpack the variable coef_wp2rtp_implicit from the structure
      // pdf_implicit_coefs_terms.  The values of coef_wp2thlp_implicit,
      // coef_wp2up_implicit, coef_wp2vp_implict, and coef_wp2sclrp_implicit
      // are all equal to coef_wp2rtp_implicit.  The PDF parameters and the
      // resulting implicit coefficients are calculated on thermodynamic
      // levels.
      for (i = 0; i < ngrdcol; i += 1) {
        Kokkos::parallel_for(
            "calc_xm_wpxp_ta_terms::array23", coef_wp2rtp_implicit.extent(1),
            KOKKOS_LAMBDA(int cell) {
              const int slice_0 = cell;
              coef_wp2rtp_implicit(i, slice_0) =
                  pdf_implicit_coefs_terms_coef_wp2rtp_implicit(i, slice_0);
            });
        Kokkos::parallel_for(
            "calc_xm_wpxp_ta_terms::array24", coef_wp2thlp_implicit.extent(1),
            KOKKOS_LAMBDA(int cell) {
              const int slice_0                 = cell;
              coef_wp2thlp_implicit(i, slice_0) = coef_wp2rtp_implicit(i, slice_0);
            });
      }
      // Calculate the LHS turbulent advection term for <w'rt'>
      xpyp_term_ta_pdf_lhs(nzm, nzt, ngrdcol, gr, coef_wp2rtp_implicit, // In
                           rho_ds_zt, rho_ds_zm,                        // In
                           invrs_rho_ds_zm,                             // In
                           l_dummy_false,                               // In
                           sgn_t_vel_wprtp,                             // In
                           coef_wp2rtp_implicit_zm,                     // In
                           lhs_ta_wprtp                                 // Out
      );
      // For the new hybrid PDF, the LHS turbulent advection terms for
      // <w'r_t'>, <w'thl'>, and <w'sclr'> are all the same.
      Kokkos::parallel_for(
          "calc_xm_wpxp_ta_terms::array25",
          ((lhs_ta_wpthlp.extent(0) - 1) / 1 + 1) * ((lhs_ta_wpthlp.extent(1) - 1) / 1 + 1) *
              ((lhs_ta_wpthlp.extent(2) - 1) / 1 + 1),
          KOKKOS_LAMBDA(int cell) {
            const int slice_0 = cell % ((lhs_ta_wpthlp.extent(0) - 1) / 1 + 1);
            const int slice_1 = (cell / (((lhs_ta_wpthlp.extent(0) - 1) / 1 + 1))) %
                                ((lhs_ta_wpthlp.extent(1) - 1) / 1 + 1);
            const int slice_2 = cell / (((lhs_ta_wpthlp.extent(0) - 1) / 1 + 1) *
                                        ((lhs_ta_wpthlp.extent(1) - 1) / 1 + 1));
            lhs_ta_wpthlp(slice_0, slice_1, slice_2) = lhs_ta_wprtp(slice_0, slice_1, slice_2);
          });
      if (l_scalar_calc) {
        for (sclr = 0; sclr < sclr_dim; sclr += 1) {
          Kokkos::parallel_for(
              "calc_xm_wpxp_ta_terms::array26",
              ((lhs_ta_wpsclrp.extent(0) - 1) / 1 + 1) * ((lhs_ta_wpsclrp.extent(1) - 1) / 1 + 1) *
                  ((lhs_ta_wpsclrp.extent(2) - 1) / 1 + 1),
              KOKKOS_LAMBDA(int cell) {
                const int slice_0 = cell % ((lhs_ta_wpsclrp.extent(0) - 1) / 1 + 1);
                const int slice_1 = (cell / (((lhs_ta_wpsclrp.extent(0) - 1) / 1 + 1))) %
                                    ((lhs_ta_wpsclrp.extent(1) - 1) / 1 + 1);
                const int slice_2 = cell / (((lhs_ta_wpsclrp.extent(0) - 1) / 1 + 1) *
                                            ((lhs_ta_wpsclrp.extent(1) - 1) / 1 + 1));
                lhs_ta_wpsclrp(slice_0, slice_1, slice_2, sclr) =
                    lhs_ta_wprtp(slice_0, slice_1, slice_2);
              });
        }
      }
      if (stats.l_sample) {
        Kokkos::parallel_for(
            "calc_xm_wpxp_ta_terms::array27",
            term_wp2rtp_explicit.extent(0) * term_wp2rtp_explicit.extent(1),
            KOKKOS_LAMBDA(int cell) {
              const int slice_0                      = (cell % term_wp2rtp_explicit.extent(0));
              const int slice_1                      = (cell / term_wp2rtp_explicit.extent(0));
              term_wp2rtp_explicit(slice_0, slice_1) = zero;
            });
        Kokkos::parallel_for(
            "calc_xm_wpxp_ta_terms::array28",
            term_wp2thlp_explicit.extent(0) * term_wp2thlp_explicit.extent(1),
            KOKKOS_LAMBDA(int cell) {
              const int slice_0                       = (cell % term_wp2thlp_explicit.extent(0));
              const int slice_1                       = (cell / term_wp2thlp_explicit.extent(0));
              term_wp2thlp_explicit(slice_0, slice_1) = zero;
            });
      }
      // The <w'r_t'>, <w'thl'>, <w'sclr'> turbulent advection terms are
      // entirely implicit.  Set the RHS turbulent advection terms to 0
      Kokkos::parallel_for(
          "calc_xm_wpxp_ta_terms::array29",
          ((rhs_ta_wprtp.extent(0) - 1) / 1 + 1) * ((rhs_ta_wprtp.extent(1) - 1) / 1 + 1),
          KOKKOS_LAMBDA(int cell) {
            const int slice_0              = cell % ((rhs_ta_wprtp.extent(0) - 1) / 1 + 1);
            const int slice_1              = cell / (((rhs_ta_wprtp.extent(0) - 1) / 1 + 1));
            rhs_ta_wprtp(slice_0, slice_1) = zero;
          });
      Kokkos::parallel_for(
          "calc_xm_wpxp_ta_terms::array30",
          ((rhs_ta_wpthlp.extent(0) - 1) / 1 + 1) * ((rhs_ta_wpthlp.extent(1) - 1) / 1 + 1),
          KOKKOS_LAMBDA(int cell) {
            const int slice_0               = cell % ((rhs_ta_wpthlp.extent(0) - 1) / 1 + 1);
            const int slice_1               = cell / (((rhs_ta_wpthlp.extent(0) - 1) / 1 + 1));
            rhs_ta_wpthlp(slice_0, slice_1) = zero;
          });
      Kokkos::parallel_for(
          "calc_xm_wpxp_ta_terms::array31",
          ((rhs_ta_wpsclrp.extent(0) - 1) / 1 + 1) * ((rhs_ta_wpsclrp.extent(1) - 1) / 1 + 1) *
              ((rhs_ta_wpsclrp.extent(2) - 1) / 1 + 1),
          KOKKOS_LAMBDA(int cell) {
            const int slice_0 = cell % ((rhs_ta_wpsclrp.extent(0) - 1) / 1 + 1);
            const int slice_1 = (cell / (((rhs_ta_wpsclrp.extent(0) - 1) / 1 + 1))) %
                                ((rhs_ta_wpsclrp.extent(1) - 1) / 1 + 1);
            const int slice_2 = cell / (((rhs_ta_wpsclrp.extent(0) - 1) / 1 + 1) *
                                        ((rhs_ta_wpsclrp.extent(1) - 1) / 1 + 1));
            rhs_ta_wpsclrp(slice_0, slice_1, slice_2) = zero;
          });
      if (l_predict_upwp_vpwp) {
        // Predict <u> and <u'w'>, as well as <v> and <v'w'>.
        // These terms are equal to the <w'r_t'> terms as well in this case
        Kokkos::parallel_for(
            "calc_xm_wpxp_ta_terms::array32",
            ((lhs_ta_wpup.extent(0) - 1) / 1 + 1) * ((lhs_ta_wpup.extent(1) - 1) / 1 + 1) *
                ((lhs_ta_wpup.extent(2) - 1) / 1 + 1),
            KOKKOS_LAMBDA(int cell) {
              const int slice_0 = cell % ((lhs_ta_wpup.extent(0) - 1) / 1 + 1);
              const int slice_1 = (cell / (((lhs_ta_wpup.extent(0) - 1) / 1 + 1))) %
                                  ((lhs_ta_wpup.extent(1) - 1) / 1 + 1);
              const int slice_2 = cell / (((lhs_ta_wpup.extent(0) - 1) / 1 + 1) *
                                          ((lhs_ta_wpup.extent(1) - 1) / 1 + 1));
              lhs_ta_wpup(slice_0, slice_1, slice_2) = lhs_ta_wprtp(slice_0, slice_1, slice_2);
            });
        Kokkos::parallel_for(
            "calc_xm_wpxp_ta_terms::array33",
            ((lhs_ta_wpvp.extent(0) - 1) / 1 + 1) * ((lhs_ta_wpvp.extent(1) - 1) / 1 + 1) *
                ((lhs_ta_wpvp.extent(2) - 1) / 1 + 1),
            KOKKOS_LAMBDA(int cell) {
              const int slice_0 = cell % ((lhs_ta_wpvp.extent(0) - 1) / 1 + 1);
              const int slice_1 = (cell / (((lhs_ta_wpvp.extent(0) - 1) / 1 + 1))) %
                                  ((lhs_ta_wpvp.extent(1) - 1) / 1 + 1);
              const int slice_2 = cell / (((lhs_ta_wpvp.extent(0) - 1) / 1 + 1) *
                                          ((lhs_ta_wpvp.extent(1) - 1) / 1 + 1));
              lhs_ta_wpvp(slice_0, slice_1, slice_2) = lhs_ta_wprtp(slice_0, slice_1, slice_2);
            });
        // The <w'u'> and <w'v'> turbulent advection terms are entirely
        // implicit.  Set the RHS turbulent advection terms to 0
        Kokkos::parallel_for(
            "calc_xm_wpxp_ta_terms::array34",
            ((rhs_ta_wpup.extent(0) - 1) / 1 + 1) * ((rhs_ta_wpup.extent(1) - 1) / 1 + 1),
            KOKKOS_LAMBDA(int cell) {
              const int slice_0             = cell % ((rhs_ta_wpup.extent(0) - 1) / 1 + 1);
              const int slice_1             = cell / (((rhs_ta_wpup.extent(0) - 1) / 1 + 1));
              rhs_ta_wpup(slice_0, slice_1) = zero;
            });
        Kokkos::parallel_for(
            "calc_xm_wpxp_ta_terms::array35",
            ((rhs_ta_wpvp.extent(0) - 1) / 1 + 1) * ((rhs_ta_wpvp.extent(1) - 1) / 1 + 1),
            KOKKOS_LAMBDA(int cell) {
              const int slice_0             = cell % ((rhs_ta_wpvp.extent(0) - 1) / 1 + 1);
              const int slice_1             = cell / (((rhs_ta_wpvp.extent(0) - 1) / 1 + 1));
              rhs_ta_wpvp(slice_0, slice_1) = zero;
            });
      }
    }
  }
  if (stats.l_sample) {
    stats_update("coef_wp2rtp_implicit", coef_wp2rtp_implicit, // In
                 stats                                         // InOut
    );
    stats_update("term_wp2rtp_explicit", term_wp2rtp_explicit, // In
                 stats                                         // InOut
    );
    stats_update("coef_wp2thlp_implicit", coef_wp2thlp_implicit, // In
                 stats                                           // InOut
    );
    stats_update("term_wp2thlp_explicit", term_wp2thlp_explicit, // In
                 stats                                           // InOut
    );
  }
}

// nzm: Number of momentum vertical levels
// nzt: Number of thermodynamic vertical levels
// ngrdcol: Number of grid columns
// sclr_dim: Number of passive scalars
// sclr_tol: Threshold(s) on the passive scalars  [units vary]
// dt: Timestep                                 [s]
// nrhs: Number of RHS vectors
// wm_zt: w wind component on thermodynamic levels [m/s]
// wp2: w'^2 (momentum levels)                   [m^2/s^2]
// invrs_tau_C6_zm: Inverse tau on momentum levels applied to C6 term [1/s]
// tau_max_zm: Max. allowable eddy dissipation time scale on m-levs [s]
// rtpthvp: r_t'th_v' (momentum levels)              [(kg/kg) K]
// rtm_forcing: r_t forcing (thermodynamic levels)       [(kg/kg)/s]
// wprtp_forcing: <w'r_t'> forcing (momentum levels)       [(kg/kg)/s^2]
// thlpthvp: th_l'th_v' (momentum levels)             [K^2]
// thlm_forcing: th_l forcing (thermodynamic levels)      [K/s]
// wpthlp_forcing: <w'th_l'> forcing (momentum levels)      [K/s^2]
// rho_ds_zm: Dry, static density on momentum levels   [kg/m^3]
// rho_ds_zt: Dry, static density on thermo. levels    [kg/m^3]
// invrs_rho_ds_zm: Inv. dry, static density @ moment. levs. [m^3/kg]
// invrs_rho_ds_zt: Inv. dry, static density @ thermo. levs. [m^3/kg]
// thv_ds_zm: Dry, base-state theta_v on moment. levs. [K]
// rtp2: r_t'^2 (momentum levels)                 [(kg/kg)^2]
// thlp2: th_l'^2 (momentum levels)                [K^2]
// l_implemented: True if CLUBB is being implemented and run in a host model
// sclrpthvp: <sclr' th_v'> (momentum levels)       [Units vary]
// sclrm_forcing: sclrm forcing (thermodynamic levels)  [Units vary]
// sclrp2: For clipping Vince Larson             [Units vary]
// um_forcing: <u> forcing term (thermodynamic levels)      [m/s^2]
// vm_forcing: <v> forcing term (thermodynamic levels)      [m/s^2]
// ug: <u> geostrophic wind (thermodynamic levels)  [m/s]
// vg: <v> geostrophic wind (thermodynamic levels)  [m/s]
// uprcp: < u' r_c' >                               [(m kg)/(s kg)]
// vprcp: < v' r_c' >                               [(m kg)/(s kg)]
// rc_coef_zm: Coefficient on X'r_c' in X'th_v' equation [K/(kg/kg)]
// fcor: Traditional Coriolis parameter               [s^-1] Vertical planetary vorticity.
// Proportional to sin(latitude)
// fcor_y: Nontraditional Coriolis parameter            [s^-1] Meridional planetary vorticity.
// Proportional to cos(latitude)
// up2: Variance of the u wind component             [m^2/s^2]
// vp2: Variance of the v wind component             [m^2/s^2]
// low_lev_effect: Index of the lowest level that has an effect.
// high_lev_effect: Index of the highest level that has an effect.
// lhs_diff_zm: Diffusion term for w'x'
// lhs_diff_zt: Diffusion term for w'x'
// lhs_ma_zt: Mean advection contributions to lhs
// lhs_ma_zm: Mean advection contributions to lhs
// lhs_ta_wpxp: w'r_t' turbulent advection contributions to lhs
// rhs_ta_wprtp: w'r_t' turbulent advection contributions to rhs
// rhs_ta_wpthlp: w'thl' turbulent advection contributions to rhs
// rhs_ta_wpup: w'u' turbulent advection contributions to rhs
// rhs_ta_wpvp: w'v' turbulent advection contributions to rhs
// rhs_ta_wpsclrp: w'sclr' turbulent advection contributions to rhs
// lhs_tp: Turbulent production terms of w'x'
// lhs_ta_xm: Turbulent advection terms of xm
// lhs_ac_pr2: Accumulation of w'x' and w'x' pressure term 2
// lhs_pr1_wprtp: Pressure term 1 for w'r_t' for all grid levels
// lhs_pr1_wpthlp: Pressure term 1 for w'thl' for all grid levels
// lhs_pr1_wpsclrp: Pressure term 1 for w'sclr' for all grid levels
// C_uu_shr: CLUBB tunable parameter C_uu_shr
// penta_solve_method: Method to solve then penta-diagonal system
// tridiag_solve_method: Specifier for method to solve tridiagonal systems, used for monotonic
// flux limiter
// fill_holes_type: Specifier for which hole filling method to use
// l_predict_upwp_vpwp: Flag to predict <u'w'> and <v'w'> along with <u> and <v> alongside the
// advancement of <rt>, <w'rt'>, <thl>, <wpthlp>, <sclr>, and <w'sclr'> in subroutine
// advance_xm_wpxp. Otherwise, <u'w'> and <v'w'> are still approximated by eddy diffusivity when
// <u> and <v> are advanced in subroutine advance_windm_edsclrm.
// l_ho_nontrad_coriolis: Flag to implement the nontraditional Coriolis terms in the prognostic
// equations of <w'w'>, <u'w'>, and <u'u'>.
// l_ho_trad_coriolis: Flag to implement the traditional Coriolis terms in the prognostic
// equations of <v'w'> and <u'w'>.
// l_diffuse_rtm_and_thlm: This flag determines whether or not we want CLUBB to do diffusion on
// rtm and thlm
// l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent or mean advection terms. It
// affects rtm, thlm, sclrm, um and vm.
// l_enable_relaxed_clipping: Flag to relax clipping on wpxp in xm_wpxp_clipping_and_stats
// l_perturbed_wind: Whether perturbed winds are being solved
// l_mono_flux_lim_thlm: Flag to turn on monotonic flux limiter for thlm
// l_mono_flux_lim_rtm: Flag to turn on monotonic flux limiter for rtm
// l_mono_flux_lim_um: Flag to turn on monotonic flux limiter for um
// l_mono_flux_lim_vm: Flag to turn on monotonic flux limiter for vm
// l_mono_flux_lim_spikefix: Flag to implement monotonic flux limiter code that eliminates
// spurious drying tendencies at model top
// rtm: r_t  (total water mixing ratio)           [kg/kg]
// wprtp: w'r_t'                                    [(kg/kg) m/s]
// thlm: th_l (liquid water potential temperature) [K]
// wpthlp: w'th_l'                                   [K m/s]
// sclrm: [Units vary]
// wpsclrp: [Units vary]
// um: <u>:  mean west-east horiz. velocity (thermo. levs.)   [m/s]
// upwp: <u'w'>:  momentum flux (momentum levels)               [m^2/s^2]
// vm: <v>:  mean south-north horiz. velocity (thermo. levs.) [m/s]
// vpwp: <v'w'>:  momentum flux (momentum levels)               [m^2/s^2]
// um_pert: perturbed <u>    [m/s]
// vm_pert: perturbed <v>    [m/s]
// upwp_pert: perturbed <u'w'> [m^2/s^2]
// vpwp_pert: perturbed <v'w'> [m^2/s^2]
// err_info: err_info struct containing err_code and err_header
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26, class Array27, class Array28, class Array29,
          class Array30, class Array31, class Array32, class Array33, class Array34, class Array35,
          class Array36, class Array37, class Array38, class Array39, class Array40, class Array41,
          class Array42, class Array43, class Array44, class Array45, class Array46, class Array47,
          class Array48, class Array49, class Array50, class Array51, class Array52, class Array53,
          class Array54, class Array55, class Array56, class Array57, class Array58, class Array59,
          class Array60, class Array61, class Array62, class Array63, class Array64, class Array65,
          class Array66, class Array67>
void solve_xm_wpxp_with_single_lhs(
    ScratchWorkspace& workspace,                                                    // InOut
    int nzm, int nzt, int ngrdcol, int sclr_dim, Array0 sclr_tol, const grid& gr,   // In
    real_t dt,                                                                      // In
    bool l_iter, int nrhs,                                                          // In
    Array1 wm_zt, Array2 wp2, Array3 invrs_tau_C6_zm, Array4 tau_max_zm,            // In
    Array5 rtpthvp, Array6 rtm_forcing, Array7 wprtp_forcing, Array8 thlpthvp,      // In
    Array9 thlm_forcing, Array10 wpthlp_forcing, Array11 rho_ds_zm,                 // In
    Array12 rho_ds_zt, Array13 invrs_rho_ds_zm, Array14 invrs_rho_ds_zt,            // In
    Array15 thv_ds_zm, Array16 rtp2, Array17 thlp2, bool l_implemented,             // In
    Array18 sclrpthvp, Array19 sclrm_forcing, Array20 sclrp2, Array21 um_forcing,   // In
    Array22 vm_forcing, Array23 ug, Array24 vg, Array25 uprcp, Array26 vprcp,       // In
    Array27 rc_coef_zm, Array28 fcor,                                               // In
    Array29 fcor_y, Array30 up2, Array31 vp2,                                       // In
    Array32 low_lev_effect, Array33 high_lev_effect,                                // In
    Array34 C6rt_Skw_fnc, Array35 C6thl_Skw_fnc, Array36 C7_Skw_fnc,                // In
    Array37 lhs_diff_zm, Array38 lhs_diff_zt, Array39 lhs_ma_zt, Array40 lhs_ma_zm, // In
    Array41 lhs_ta_wpxp,                                                            // In
    Array42 rhs_ta_wprtp, Array43 rhs_ta_wpthlp, Array44 rhs_ta_wpup,               // In
    Array45 rhs_ta_wpvp, Array46 rhs_ta_wpsclrp,                                    // In
    Array47 lhs_tp, Array48 lhs_ta_xm, Array49 lhs_ac_pr2, Array50 lhs_pr1_wprtp,   // In
    Array51 lhs_pr1_wpthlp, Array52 lhs_pr1_wpsclrp,                                // In
    Array53 C_uu_shr,                                                               // In
    int penta_solve_method,                                                         // In
    int tridiag_solve_method,                                                       // In
    int fill_holes_type,                                                            // In
    bool l_predict_upwp_vpwp,                                                       // In
    bool l_ho_nontrad_coriolis,                                                     // In
    bool l_ho_trad_coriolis,                                                        // In
    bool l_diffuse_rtm_and_thlm,                                                    // In
    bool l_upwind_xm_ma,                                                            // In
    bool l_enable_relaxed_clipping,                                                 // In
    bool l_perturbed_wind,                                                          // In
    bool l_mono_flux_lim_thlm,                                                      // In
    bool l_mono_flux_lim_rtm,                                                       // In
    bool l_mono_flux_lim_um,                                                        // In
    bool l_mono_flux_lim_vm,                                                        // In
    bool l_mono_flux_lim_spikefix,                                                  // In
    int& wprtp_cl_num,                                                              // InOut
    int& wpthlp_cl_num, int& upwp_cl_num, int& vpwp_cl_num,                         // InOut
    stats_type& stats,                                                              // InOut
    Array54 rtm, Array55 wprtp, Array56 thlm, Array57 wpthlp,                       // InOut
    Array58 sclrm, Array59 wpsclrp, Array60 um, Array61 upwp, Array62 vm,           // InOut
    Array63 vpwp,                                                                   // InOut
    Array64 um_pert, Array65 vm_pert, Array66 upwp_pert, Array67 vpwp_pert,         // InOut
    err_info_type& err_info                                                         // InOut
) {
  //
  // Description: This subroutine solves all xm_wpxp when all the LHS matrices are equal.
  //              The LHS matrices being equivalent allows for only a single solve, rather
  //              than a seperate solve for each field.
  //----------------------------------------------------------------------------------------
  // ------------------- Input Variables -------------------
  // Additional variables for passive scalars
  // Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
  // Vertical planetary vorticity.   Proportional to sin(latitude)
  // Meridional planetary vorticity. Proportional to cos(latitude)
  // LHS/RHS terms
  // Variables used as part of the monotonic turbulent advection scheme.
  // Find the lowermost and uppermost grid levels that can have an effect
  // on the central thermodynamic level during the course of a time step,
  // due to the effects of turbulent advection only.
  // used for monotonic flux limiter
  // with <u> and <v> alongside the advancement
  // of <rt>, <w'rt'>, <thl>, <wpthlp>, <sclr>,
  // and <w'sclr'> in subroutine advance_xm_wpxp.
  // Otherwise, <u'w'> and <v'w'> are still
  // approximated by eddy diffusivity when <u>
  // and <v> are advanced in subroutine
  // advance_windm_edsclrm.
  // prognostic equations of <w'w'>, <u'w'>, and <u'u'>.
  // prognostic equations of <v'w'> and <u'w'>.
  // CLUBB to do diffusion on rtm and thlm
  // an upwind differencing approximation rather
  // than a centered differencing for turbulent
  // or mean advection terms. It affects rtm,
  // thlm, sclrm, um and vm.
  // i.e. TKE = 1/2 (u'^2 + v'^2 + w'^2)
  // xm_wpxp_clipping_and_stats
  // eliminates spurious drying tendencies at model top
  // ------------------- Input/Output Variables -------------------
  // Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
  // Variables used to track perturbed version of winds.
  // ------------------- Local Variables -------------------
  // Additional variables for passive scalars
  // Constant parameters as a function of Skw.
  auto workspace_scope     = workspace.scope();
  const auto sclr_tol_host = to_host_space(workspace, sclr_tol);
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::lhs", nsub + nsup + 1,
          ngrdcol,
          2 * nzm - 1); // Implicit contributions to wpxp/xm (band diag. matrix) (LAPACK)
  Kokkos::View<real_t***, Kokkos::LayoutLeft> wpsclrp_forcing =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::wpsclrp_forcing", ngrdcol,
          nzm,
          sclr_dim); // <w'sclr'> forcing (momentum levels)  [m/s{un vary}]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> um_tndcy =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::um_tndcy", ngrdcol,
          nzt); // <u> forcing term + coriolis (thermo levs)        [m/s^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vm_tndcy =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::vm_tndcy", ngrdcol,
          nzt); // <v> forcing term + coriolis (thermo levs)        [m/s^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> upwp_forcing =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::upwp_forcing", ngrdcol,
          nzm); // <u'w'> extra RHS pressure term (mom levs)        [m^2/s^3]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vpwp_forcing =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::vpwp_forcing", ngrdcol,
          nzm); // <v'w'> extra RHS pressure term (mom levs)        [m^2/s^3]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> upthvp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::upthvp", ngrdcol,
          nzm); // <u'thv'> (momentum levels)                       [m/s K]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vpthvp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::vpthvp", ngrdcol,
          nzm); // <v'thv'> (momentum levels)                       [m/s K]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> upthlp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::upthlp", ngrdcol,
          nzm); // eastward horz turb flux of theta_l (mom levs)    [m/s K]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vpthlp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::vpthlp", ngrdcol,
          nzm); // northward horz turb flux of theta_l (mom levs)   [m/s K]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> uprtp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::uprtp", ngrdcol,
          nzm); // eastward horz turb flux of tot water (mom levs)  [m/s kg/kg]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vprtp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::vprtp", ngrdcol,
          nzm); // northward horz turb flux of tot water (mom levs) [m/s kg/kg]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> tau_C6_zm =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::tau_C6_zm", ngrdcol,
          nzm); // Time-scale tau on momentum levels applied to C6 term [s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> upwp_forcing_pert =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::upwp_forcing_pert", ngrdcol,
          nzm); // perturbed extra RHS term (mom levs)         [m^2/s^3]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vpwp_forcing_pert =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::vpwp_forcing_pert", ngrdcol,
          nzm); // perturbed extra RHS term (mom levs)         [m^2/s^3]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> upthvp_pert =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::upthvp_pert", ngrdcol,
          nzm); // perturbed <u'thv'> (momentum levels)        [m/s K]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vpthvp_pert =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::vpthvp_pert", ngrdcol,
          nzm); // perturbed <v'thv'> (momentum levels)        [m/s K]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> upthlp_pert =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::upthlp_pert", ngrdcol,
          nzm); // perturbed horz flux of theta_l (mom levs)   [m/s K]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vpthlp_pert =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::vpthlp_pert", ngrdcol,
          nzm); // perturbed horz flux of theta_l (mom levs)   [m/s K]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> uprtp_pert =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::uprtp_pert", ngrdcol,
          nzm); // perturbed horz flux of tot water (mom levs) [m/s kg/kg]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vprtp_pert =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::vprtp_pert", ngrdcol,
          nzm); // perturbed horz flux of tot water (mom levs) [m/s kg/kg]
  Kokkos::View<real_t***, Kokkos::LayoutLeft> rhs =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::rhs", ngrdcol, 2 * nzm - 1,
          nrhs); // Right-hand sides of band diag. matrix. (LAPACK)
  Kokkos::View<real_t***, Kokkos::LayoutLeft> rhs_save =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::rhs_save", ngrdcol,
          2 * nzm - 1,
          nrhs); // Saved Right-hand sides of band diag. matrix. (LAPACK)
  Kokkos::View<real_t***, Kokkos::LayoutLeft> solution =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::solution", ngrdcol,
          2 * nzm - 1,
          nrhs); // solution vectors of band diag. matrix. (LAPACK)
  Kokkos::View<real_t***, Kokkos::LayoutLeft> old_solution =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::old_solution", ngrdcol,
          2 * nzm - 1, nrhs); // previous solutions
  Kokkos::View<real_t*, Kokkos::LayoutLeft> rcond =
      workspace_scope.view<Kokkos::View<real_t*, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::rcond", ngrdcol);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> zeros_vector =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::zeros_vector", ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> ddzt_um =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::ddzt_um", ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> ddzt_vm =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::ddzt_vm", ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> ddzt_um_pert =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::ddzt_um_pert", ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> ddzt_vm_pert =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::ddzt_vm_pert", ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> upthvp_tmp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::upthvp_tmp", ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vpthvp_tmp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::vpthvp_tmp", ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> um_smth =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::um_smth", ngrdcol, nzt);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vm_smth =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::vm_smth", ngrdcol, nzt);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> tmp_zt =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::tmp_zt", ngrdcol, nzt);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> tmp_zm =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_single_lhs::tmp_zm", ngrdcol, nzm);
  int i;
  int k;
  int j;
  int sclr;
  // ------------------- Begin Code -------------------
  // This is initialized solely for the purpose of avoiding a compiler
  // warning about uninitialized variables.
  Kokkos::parallel_for(
      "solve_xm_wpxp_with_single_lhs::loop1", std::max(0, nzm) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i        = (cell % ngrdcol);
        const int k        = (cell / ngrdcol);
        zeros_vector(i, k) = zero;
      });
  // Simple case, where the new PDF is
  // used, l_explicit_turbulent_adv_wpxp is enabled.
  // Create the lhs once
  xm_wpxp_lhs(nzm, nzt, ngrdcol, l_iter, dt,                // In
              l_implemented, lhs_diff_zm, lhs_diff_zt,      // In
              lhs_ma_zm, lhs_ma_zt, lhs_ta_wpxp, lhs_ta_xm, // In
              lhs_tp, lhs_pr1_wprtp, lhs_ac_pr2,            // In
              l_diffuse_rtm_and_thlm,                       // In
              lhs                                           // Out
  );
  // Compute the explicit portion of the r_t and w'r_t' equations.
  // Build the right-hand side vector.
  xm_wpxp_rhs(workspace, nzm, nzt, ngrdcol, gr, xm_wpxp_rtm, l_iter, dt, // In
              rtm, wprtp, rtm_forcing, wprtp_forcing, C7_Skw_fnc,        // In
              rtpthvp, rhs_ta_wprtp, thv_ds_zm,                          // In
              lhs_pr1_wprtp, lhs_ta_wpxp,                                // In
              stats,                                                     // InOut
              Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, 0)          // Out
  );
  // Compute the explicit portion of the th_l and w'th_l' equations.
  // Build the right-hand side vector.
  xm_wpxp_rhs(workspace, nzm, nzt, ngrdcol, gr, xm_wpxp_thlm, l_iter, dt, // In
              thlm, wpthlp, thlm_forcing, wpthlp_forcing, C7_Skw_fnc,     // In
              thlpthvp, rhs_ta_wpthlp, thv_ds_zm,                         // In
              lhs_pr1_wpthlp, lhs_ta_wpxp,                                // In
              stats,                                                      // InOut
              Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, 1)           // Out
  );
  // ---> h1g, 2010-06-15
  // scalar transport, e.g, droplet and ice number concentration
  // are handled in  " advance_sclrm_Nd_module.F90 "
  for (sclr = 0; sclr < sclr_dim; sclr += 1) {
    // <--- h1g, 2010-06-15
    // Set <w'sclr'> forcing to 0 unless unless testing the wpsclrp code
    // using wprtp or wpthlp (then use wprtp_forcing or wpthlp_forcing).
    Kokkos::parallel_for(
        "solve_xm_wpxp_with_single_lhs::loop2", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i                 = (cell % ngrdcol);
          const int k                 = (cell / ngrdcol);
          wpsclrp_forcing(i, k, sclr) = zero;
        });
    xm_wpxp_rhs(workspace, nzm, nzt, ngrdcol, gr, xm_wpxp_scalar, l_iter, dt,     // In
                Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, sclr),           // In
                Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, sclr),         // In
                Kokkos::subview(sclrm_forcing, Kokkos::ALL, Kokkos::ALL, sclr),   // In
                Kokkos::subview(wpsclrp_forcing, Kokkos::ALL, Kokkos::ALL, sclr), // In
                C7_Skw_fnc,                                                       // In
                Kokkos::subview(sclrpthvp, Kokkos::ALL, Kokkos::ALL, sclr),       // In
                Kokkos::subview(rhs_ta_wpsclrp, Kokkos::ALL, Kokkos::ALL, sclr),  // In
                thv_ds_zm,                                                        // In
                lhs_pr1_wpsclrp, lhs_ta_wpxp,                                     // In
                stats,                                                            // InOut
                Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, (sclr + 2))        // Out
    );
  }
  if (l_predict_upwp_vpwp) {
    // Predict <u> and <u'w'>, as well as <v> and <v'w'>.
    // Currently, this requires the ADG1 PDF with implicit turbulent advection.
    // l_explicit_turbulent_adv_wpxp = false
    // and ( iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_new_hybrid )
    // Coriolis term for <u> and <v>
    if (!l_implemented) {
      // Only compute the Coriolis term if the model is running on its own,
      // and is not part of a larger, host model.
      Kokkos::parallel_for(
          "solve_xm_wpxp_with_single_lhs::loop3", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i    = (cell % ngrdcol);
            const int k    = (cell / ngrdcol);
            um_tndcy(i, k) = um_forcing(i, k) - fcor(i) * (vg(i, k) - vm(i, k));
            vm_tndcy(i, k) = vm_forcing(i, k) + fcor(i) * (ug(i, k) - um(i, k));
          });
      if (stats.l_sample) {
        // um or vm term gf is completely explicit; call stat_update_var.
        // Geostrophic-forcing contribution to um tendency.
        Kokkos::parallel_for(
            "solve_xm_wpxp_with_single_lhs::loop4", std::max(0, nzt) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              const int i  = (cell % ngrdcol);
              const int k  = (cell / ngrdcol);
              tmp_zt(i, k) = -fcor(i) * vg(i, k);
            });
        stats_update("um_gf", tmp_zt, // In
                     stats            // InOut
        );
        // Geostrophic-forcing contribution to vm tendency.
        Kokkos::parallel_for(
            "solve_xm_wpxp_with_single_lhs::loop5", std::max(0, nzt) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              const int i  = (cell % ngrdcol);
              const int k  = (cell / ngrdcol);
              tmp_zt(i, k) = fcor(i) * ug(i, k);
            });
        stats_update("vm_gf", tmp_zt, // In
                     stats            // InOut
        );
        // um or vm term cf is completely explicit; call stat_update_var.
        // Coriolis-coupling contribution to um tendency.
        Kokkos::parallel_for(
            "solve_xm_wpxp_with_single_lhs::loop6", std::max(0, nzt) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              const int i  = (cell % ngrdcol);
              const int k  = (cell / ngrdcol);
              tmp_zt(i, k) = fcor(i) * vm(i, k);
            });
        stats_update("um_cf", tmp_zt, // In
                     stats            // InOut
        );
        // Coriolis-coupling contribution to vm tendency.
        Kokkos::parallel_for(
            "solve_xm_wpxp_with_single_lhs::loop7", std::max(0, nzt) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              const int i  = (cell % ngrdcol);
              const int k  = (cell / ngrdcol);
              tmp_zt(i, k) = -fcor(i) * um(i, k);
            });
        stats_update("vm_cf", tmp_zt, // In
                     stats            // InOut
        );
        // um or vm forcing term
        // External forcing terms are explicit.
        stats_update("um_f", um_forcing, // In
                     stats               // InOut
        );
        stats_update("vm_f", vm_forcing, // In
                     stats               // InOut
        );
      }
    } else {
      Kokkos::parallel_for(
          "solve_xm_wpxp_with_single_lhs::loop8", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i    = (cell % ngrdcol);
            const int k    = (cell / ngrdcol);
            um_tndcy(i, k) = zero;
            vm_tndcy(i, k) = zero;
          });
    }
    ddzt(nzm, nzt, ngrdcol, gr, to_execution_space(workspace, um), ddzt_um);
    ddzt(nzm, nzt, ngrdcol, gr, to_execution_space(workspace, vm), ddzt_vm);
    // Add "extra term" and optional Coriolis term for <u'w'> and <v'w'>.
    Kokkos::parallel_for(
        "solve_xm_wpxp_with_single_lhs::loop9", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i        = (cell % ngrdcol);
          const int k        = (cell / ngrdcol);
          upwp_forcing(i, k) = C_uu_shr(i) * wp2(i, k) * ddzt_um(i, k);
          vpwp_forcing(i, k) = C_uu_shr(i) * wp2(i, k) * ddzt_vm(i, k);
        });
    // Add optional traditional Coriolis terms for <u'w'> and <v'w'>
    // Hing Ong, 19 July 2025
    if (l_ho_trad_coriolis) {
      Kokkos::parallel_for(
          "solve_xm_wpxp_with_single_lhs::loop10", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i        = (cell % ngrdcol);
            const int k        = (cell / ngrdcol);
            upwp_forcing(i, k) = upwp_forcing(i, k) + fcor(i) * vpwp(i, k);
            vpwp_forcing(i, k) = vpwp_forcing(i, k) - fcor(i) * upwp(i, k);
          });
    }
    // Add optional nontraditional Coriolis term for <u'w'>
    // Hing Ong, 19 July 2025
    if (l_ho_nontrad_coriolis) {
      Kokkos::parallel_for(
          "solve_xm_wpxp_with_single_lhs::loop11", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i        = (cell % ngrdcol);
            const int k        = (cell / ngrdcol);
            upwp_forcing(i, k) = upwp_forcing(i, k) + fcor_y(i) * (up2(i, k) - wp2(i, k));
          });
    }
    if (l_perturbed_wind) {
      ddzt(nzm, nzt, ngrdcol, gr, to_execution_space(workspace, um_pert), ddzt_um_pert);
      ddzt(nzm, nzt, ngrdcol, gr, to_execution_space(workspace, vm_pert), ddzt_vm_pert);
      Kokkos::parallel_for(
          "solve_xm_wpxp_with_single_lhs::loop12", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i             = (cell % ngrdcol);
            const int k             = (cell / ngrdcol);
            upwp_forcing_pert(i, k) = C_uu_shr(i) * wp2(i, k) * ddzt_um_pert(i, k);
            vpwp_forcing_pert(i, k) = C_uu_shr(i) * wp2(i, k) * ddzt_vm_pert(i, k);
          });
    }
    if (stats.l_sample) {
      Kokkos::parallel_for(
          "solve_xm_wpxp_with_single_lhs::loop13", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i  = (cell % ngrdcol);
            const int k  = (cell / ngrdcol);
            tmp_zm(i, k) = C_uu_shr(i) * wp2(i, k) * ddzt_um(i, k);
          });
      stats_update("upwp_pr4", tmp_zm, // In
                   stats               // InOut
      );
      Kokkos::parallel_for(
          "solve_xm_wpxp_with_single_lhs::loop14", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i  = (cell % ngrdcol);
            const int k  = (cell / ngrdcol);
            tmp_zm(i, k) = C_uu_shr(i) * wp2(i, k) * ddzt_vm(i, k);
          });
      stats_update("vpwp_pr4", tmp_zm, // In
                   stats               // InOut
      );
      if (l_ho_trad_coriolis) {
        Kokkos::parallel_for(
            "solve_xm_wpxp_with_single_lhs::loop15", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              const int i  = (cell % ngrdcol);
              const int k  = (cell / ngrdcol);
              tmp_zm(i, k) = fcor(i) * vpwp(i, k);
            });
        stats_update("upwp_tct", tmp_zm, // In
                     stats               // InOut
        );
        Kokkos::parallel_for(
            "solve_xm_wpxp_with_single_lhs::loop16", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              const int i  = (cell % ngrdcol);
              const int k  = (cell / ngrdcol);
              tmp_zm(i, k) = -fcor(i) * upwp(i, k);
            });
        stats_update("vpwp_tct", tmp_zm, // In
                     stats               // InOut
        );
      }
      if (l_ho_nontrad_coriolis) {
        Kokkos::parallel_for(
            "solve_xm_wpxp_with_single_lhs::loop17", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              const int i  = (cell % ngrdcol);
              const int k  = (cell / ngrdcol);
              tmp_zm(i, k) = fcor_y(i) * (up2(i, k) - wp2(i, k));
            });
        stats_update("upwp_nct", tmp_zm, // In
                     stats               // InOut
        );
      }
    }
    // need tau_C6_zm for these calls
    Kokkos::parallel_for(
        "solve_xm_wpxp_with_single_lhs::loop18", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          tau_C6_zm(i, k) =
              ((one / invrs_tau_C6_zm(i, k)) < (tau_max_zm(i, k)) ? (one / invrs_tau_C6_zm(i, k))
                                                                  : (tau_max_zm(i, k)));
        });
    zt2zm2zt(workspace, nzm, nzt, ngrdcol, gr, um, um_smth);
    zt2zm2zt(workspace, nzm, nzt, ngrdcol, gr, vm, vm_smth);
    diagnose_upxp(workspace, nzm, nzt, ngrdcol, gr, upwp, thlm, wpthlp, um_smth, // In
                  C6thl_Skw_fnc, tau_C6_zm, C7_Skw_fnc,                          // In
                  upthlp                                                         // Out
    );
    diagnose_upxp(workspace, nzm, nzt, ngrdcol, gr, upwp, rtm, wprtp, um_smth, // In
                  C6rt_Skw_fnc, tau_C6_zm, C7_Skw_fnc,                         // In
                  uprtp                                                        // Out
    );
    diagnose_upxp(workspace, nzm, nzt, ngrdcol, gr, vpwp, thlm, wpthlp, vm_smth, // In
                  C6thl_Skw_fnc, tau_C6_zm, C7_Skw_fnc,                          // In
                  vpthlp                                                         // Out
    );
    diagnose_upxp(workspace, nzm, nzt, ngrdcol, gr, vpwp, rtm, wprtp, vm_smth, // In
                  C6rt_Skw_fnc, tau_C6_zm, C7_Skw_fnc,                         // In
                  vprtp                                                        // Out
    );
    if (l_perturbed_wind) {
      diagnose_upxp(workspace, nzm, nzt, ngrdcol, gr, upwp_pert, thlm, wpthlp, um_pert, // In
                    C6thl_Skw_fnc, tau_C6_zm, C7_Skw_fnc,                               // In
                    upthlp_pert                                                         // Out
      );
      diagnose_upxp(workspace, nzm, nzt, ngrdcol, gr, upwp_pert, rtm, wprtp, um_pert, // In
                    C6rt_Skw_fnc, tau_C6_zm, C7_Skw_fnc,                              // In
                    uprtp_pert                                                        // Out
      );
      diagnose_upxp(workspace, nzm, nzt, ngrdcol, gr, vpwp_pert, thlm, wpthlp, vm_pert, // In
                    C6thl_Skw_fnc, tau_C6_zm, C7_Skw_fnc,                               // In
                    vpthlp_pert                                                         // Out
      );
      diagnose_upxp(workspace, nzm, nzt, ngrdcol, gr, vpwp_pert, rtm, wprtp, vm_pert, // In
                    C6rt_Skw_fnc, tau_C6_zm, C7_Skw_fnc,                              // In
                    vprtp_pert                                                        // Out
      );
    }
    // Use a crude approximation for buoyancy terms <u'thv'> and <v'thv'>.
    //upthvp = upwp * wpthvp / max( wp2, w_tol_sqd )
    //vpthvp = vpwp * wpthvp / max( wp2, w_tol_sqd )
    //upthvp = 0.3_core_rknd * ( upthlp + 200.0_core_rknd * uprtp ) &
    //         + 200._core_rknd * sign( one, upwp) * sqrt( up2 * rcm**2 )
    //vpthvp = 0.3_core_rknd * ( vpthlp + 200.0_core_rknd * vprtp ) &
    //         + 200._core_rknd * sign( one, vpwp ) * sqrt( vp2 * rcm**2 )
    Kokkos::parallel_for(
        "solve_xm_wpxp_with_single_lhs::loop19", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          upthvp_tmp(i, k) =
              upthlp(i, k) + ep1 * thv_ds_zm(i, k) * uprtp(i, k) + rc_coef_zm(i, k) * uprcp(i, k);
          vpthvp_tmp(i, k) =
              vpthlp(i, k) + ep1 * thv_ds_zm(i, k) * vprtp(i, k) + rc_coef_zm(i, k) * vprcp(i, k);
        });
    zm2zt2zm(workspace, nzm, nzt, ngrdcol, gr, upthvp_tmp, upthvp);
    zm2zt2zm(workspace, nzm, nzt, ngrdcol, gr, vpthvp_tmp, vpthvp);
    if (l_perturbed_wind) {
      Kokkos::parallel_for(
          "solve_xm_wpxp_with_single_lhs::loop20", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i       = (cell % ngrdcol);
            const int k       = (cell / ngrdcol);
            upthvp_pert(i, k) = upthlp_pert(i, k) + ep1 * thv_ds_zm(i, k) * uprtp_pert(i, k) +
                                rc_coef_zm(i, k) * uprcp(i, k);
            vpthvp_pert(i, k) = vpthlp_pert(i, k) + ep1 * thv_ds_zm(i, k) * vprtp_pert(i, k) +
                                rc_coef_zm(i, k) * vprcp(i, k);
          });
    }
    if (stats.l_sample) {
      stats_update("upthlp", upthlp, // In
                   stats             // InOut
      );
      stats_update("uprtp", uprtp, // In
                   stats           // InOut
      );
      stats_update("vpthlp", vpthlp, // In
                   stats             // InOut
      );
      stats_update("vprtp", vprtp, // In
                   stats           // InOut
      );
      stats_update("upthvp", upthvp, // In
                   stats             // InOut
      );
      stats_update("vpthvp", vpthvp, // In
                   stats             // InOut
      );
    }
    xm_wpxp_rhs(workspace, nzm, nzt, ngrdcol, gr, xm_wpxp_um, l_iter, dt,      // In
                um, upwp, um_tndcy, upwp_forcing, C7_Skw_fnc,                  // In
                upthvp, rhs_ta_wpup, thv_ds_zm,                                // In
                lhs_pr1_wprtp, lhs_ta_wpxp,                                    // In
                stats,                                                         // InOut
                Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, (sclr_dim + 2)) // Out
    );
    xm_wpxp_rhs(workspace, nzm, nzt, ngrdcol, gr, xm_wpxp_vm, l_iter, dt,      // In
                vm, vpwp, vm_tndcy, vpwp_forcing, C7_Skw_fnc,                  // In
                vpthvp, rhs_ta_wpvp, thv_ds_zm,                                // In
                lhs_pr1_wprtp, lhs_ta_wpxp,                                    // In
                stats,                                                         // InOut
                Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, (sclr_dim + 3)) // Out
    );
    if (l_perturbed_wind) {
      xm_wpxp_rhs(workspace, nzm, nzt, ngrdcol, gr, xm_wpxp_um, l_iter, dt, um_pert, // In
                  upwp_pert, um_tndcy, upwp_forcing_pert, C7_Skw_fnc,                // In
                  upthvp_pert, rhs_ta_wpup, thv_ds_zm,                               // In
                  lhs_pr1_wprtp, lhs_ta_wpxp,                                        // In
                  stats,                                                             // InOut
                  Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, (sclr_dim + 4))     // Out
      );
      xm_wpxp_rhs(workspace, nzm, nzt, ngrdcol, gr, xm_wpxp_vm, l_iter, dt, vm_pert, // In
                  vpwp_pert, vm_tndcy, vpwp_forcing_pert, C7_Skw_fnc,                // In
                  vpthvp_pert, rhs_ta_wpvp, thv_ds_zm,                               // In
                  lhs_pr1_wprtp, lhs_ta_wpxp,                                        // In
                  stats,                                                             // InOut
                  Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, (sclr_dim + 5))     // Out
      );
    }
  }
  // Save the value of rhs, which will be overwritten with the solution as
  // part of the solving routine.
  Kokkos::parallel_for(
      "solve_xm_wpxp_with_single_lhs::loop21",
      std::max(0, nrhs) * std::max(0, 2 * nzm - 1) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i       = (cell % ngrdcol);
        const int k       = ((cell / ngrdcol) % (2 * nzm - 1));
        const int n       = (cell / (ngrdcol * (2 * nzm - 1)));
        rhs_save(i, k, n) = rhs(i, k, n);
      });
  // Use the previous solution as an initial guess for the bicgstab method
  if (penta_solve_method == penta_bicgstab) {
    Kokkos::parallel_for(
        "solve_xm_wpxp_with_single_lhs::loop22", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i                   = (cell % ngrdcol);
          const int k                   = (cell / ngrdcol);
          old_solution(i, 2 * k + 1, 0) = rtm(i, k);
          old_solution(i, 2 * k + 1, 1) = thlm(i, k);
        });
    Kokkos::parallel_for(
        "solve_xm_wpxp_with_single_lhs::loop23", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i               = (cell % ngrdcol);
          const int k               = (cell / ngrdcol);
          old_solution(i, 2 * k, 0) = wprtp(i, k);
          old_solution(i, 2 * k, 1) = wpthlp(i, k);
        });
    Kokkos::parallel_for(
        "solve_xm_wpxp_with_single_lhs::loop24",
        std::max(0, sclr_dim) * std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i                          = (cell % ngrdcol);
          const int k                          = ((cell / ngrdcol) % nzt);
          const int sclr                       = (cell / (ngrdcol * nzt));
          old_solution(i, 2 * k + 1, sclr + 2) = sclrm(i, k, sclr);
        });
    Kokkos::parallel_for(
        "solve_xm_wpxp_with_single_lhs::loop25",
        std::max(0, sclr_dim) * std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i                      = (cell % ngrdcol);
          const int k                      = ((cell / ngrdcol) % nzm);
          const int sclr                   = (cell / (ngrdcol * nzm));
          old_solution(i, 2 * k, sclr + 2) = wpsclrp(i, k, sclr);
        });
    if (l_predict_upwp_vpwp) {
      Kokkos::parallel_for(
          "solve_xm_wpxp_with_single_lhs::loop26", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i                              = (cell % ngrdcol);
            const int k                              = (cell / ngrdcol);
            old_solution(i, 2 * k + 1, sclr_dim + 2) = um(i, k);
            old_solution(i, 2 * k + 1, sclr_dim + 3) = vm(i, k);
          });
      Kokkos::parallel_for(
          "solve_xm_wpxp_with_single_lhs::loop27", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i                          = (cell % ngrdcol);
            const int k                          = (cell / ngrdcol);
            old_solution(i, 2 * k, sclr_dim + 2) = upwp(i, k);
            old_solution(i, 2 * k, sclr_dim + 3) = vpwp(i, k);
          });
    }
  }
  // Solve for all fields
  if (stats.l_sample && (var_on_stats_list(stats, "thlm_matrix_condt_num") ||
                                 var_on_stats_list(stats, "rtm_matrix_condt_num"))) {
    xm_wpxp_solve(workspace, nzm, ngrdcol, nrhs, gr, // In
                  penta_solve_method,                // In
                  l_implemented,                     // In
                  old_solution,                      // InOut
                  lhs, rhs, err_info,                // InOut
                  solution,                          // Out
                  &rcond                             // Optional Out
    );
  } else {
    xm_wpxp_solve(workspace, nzm, ngrdcol, nrhs, gr, // In
                  penta_solve_method,                // In
                  l_implemented,                     // In
                  old_solution,                      // InOut
                  lhs, rhs, err_info,                // InOut
                  solution,                          // Out
                  nullptr                            // Optional
    );
  }
  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      // Snapshot the failed system once before printing individual entries.
      const auto zm_host  = to_host_space(workspace, gr.zm);
      const auto zt_host  = to_host_space(workspace, gr.zt);
      const auto lhs_host = to_host_space(workspace, lhs);
      const auto rhs_host = to_host_space(workspace, rhs_save);
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << err_info.err_header_global << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "xm & wpxp LU decomp. failed" << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "General xm and wpxp LHS" << "\n";
        std::cerr << message.str();
      }
      for (k = 0; k < nzm - 1; k += 1) {
        for (i = 0; i < ngrdcol; i += 1) {
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "grid col = " << " " << (i + 1) << " " << "zm level = " << " " << (k + 1)
                    << " "
                    << "height [m] = " << " " << zm_host(i, k) << " "
                    << "LHS = " << " "
                    << dump_values(Kokkos::subview(lhs_host, std::make_pair(0, nsub + nsup + 1), i,
                                                   (2 * k)))
                    << "\n";
            std::cerr << message.str();
          }
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "grid col = " << " " << (i + 1) << " " << "zt level = " << " " << (k + 1)
                    << " "
                    << "height [m] = " << " " << zt_host(i, k) << " "
                    << "LHS = " << " "
                    << dump_values(Kokkos::subview(lhs_host, std::make_pair(0, nsub + nsup + 1), i,
                                                   (2 * k + 1)))
                    << "\n";
            std::cerr << message.str();
          }
        }
      }
      k = nzm - 1;
      for (i = 0; i < ngrdcol; i += 1) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "grid col = " << " " << (i + 1) << " " << "zm level = " << " " << (k + 1)
                  << " "
                  << "height [m] = " << " " << zm_host(i, k) << " "
                  << "LHS = " << " "
                  << dump_values(
                         Kokkos::subview(lhs_host, std::make_pair(0, nsub + nsup + 1), i, (2 * k)))
                  << "\n";
          std::cerr << message.str();
        }
      }
      for (j = 0; j < nrhs; j += 1) {
        if (j == 0) {
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "rtm and wprtp RHS" << "\n";
            std::cerr << message.str();
          }
        } else if (j == 1) {
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "thlm and wpthlp RHS" << "\n";
            std::cerr << message.str();
          }
        } else {
          if (sclr_dim > 0) {
            if (j < sclr_dim + 2) {
              {
                std::ostringstream message;
                message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
                message << "sclrm and wpsclrp RHS for sclr" << " " << j - 1 << "\n";
                std::cerr << message.str();
              }
            }
          }
          if (l_predict_upwp_vpwp) {
            if (j == sclr_dim + 2) {
              {
                std::ostringstream message;
                message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
                message << "um and upwp RHS" << "\n";
                std::cerr << message.str();
              }
            } else if (j == sclr_dim + 3) {
              {
                std::ostringstream message;
                message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
                message << "vm and vpwp RHS" << "\n";
                std::cerr << message.str();
              }
            }
          }
        }
        for (k = 0; k < nzm - 1; k += 1) {
          for (i = 0; i < ngrdcol; i += 1) {
            {
              std::ostringstream message;
              message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
              message << "grid col = " << " " << (i + 1) << " " << "zm level = " << " " << (k + 1)
                      << " "
                      << "height [m] = " << " " << zm_host(i, k) << " "
                      << "RHS = " << " " << rhs_host(i, 2 * k, j) << "\n";
              std::cerr << message.str();
            }
            {
              std::ostringstream message;
              message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
              message << "grid col = " << " " << (i + 1) << " " << "zt level = " << " " << (k + 1)
                      << " "
                      << "height [m] = " << " " << zt_host(i, k) << " "
                      << "RHS = " << " " << rhs_host(i, 2 * k + 1, j) << "\n";
              std::cerr << message.str();
            }
          }
        }
        k = nzm - 1;
        for (i = 0; i < ngrdcol; i += 1) {
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "grid col = " << " " << (i + 1) << " " << "zm level = " << " " << (k + 1)
                    << " "
                    << "height [m] = " << " " << zm_host(i, k) << " "
                    << "RHS = " << " " << rhs_host(i, 2 * k, j) << "\n";
            std::cerr << message.str();
          }
        }
      }
      return;
    }
  }
  xm_wpxp_clipping_and_stats(workspace, nzm, nzt, ngrdcol,          // In
                             gr, xm_wpxp_rtm, dt, wp2, rtp2, wm_zt, // In
                             rtm_forcing, rho_ds_zm, rho_ds_zt,     // In
                             invrs_rho_ds_zm, invrs_rho_ds_zt,      // In
                             (rt_tol * rt_tol), rt_tol, rcond,      // In
                             low_lev_effect, high_lev_effect,       // In
                             lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,     // In
                             lhs_diff_zm, C7_Skw_fnc,               // In
                             lhs_tp, lhs_ta_xm, lhs_pr1_wprtp,      // In
                             l_implemented,
                             Kokkos::subview(solution, Kokkos::ALL, Kokkos::ALL, 0), // In
                             tridiag_solve_method,                                   // In
                             fill_holes_type,                                        // In
                             l_upwind_xm_ma,                                         // In
                             l_enable_relaxed_clipping,                              // In
                             l_mono_flux_lim_thlm,                                   // In
                             l_mono_flux_lim_rtm,                                    // In
                             l_mono_flux_lim_um,                                     // In
                             l_mono_flux_lim_vm,                                     // In
                             l_mono_flux_lim_spikefix,                               // In
                             stats,                                                  // InOut
                             rtm,                                                    // InOut
                             rt_tol_mfl,                                             // In
                             wprtp, err_info,                                        // InOut
                             &wprtp_cl_num // Optional InOut
  );
  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << err_info.err_header_global << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "rtm monotonic flux limiter:  tridiag failed" << "\n";
        std::cerr << message.str();
      }
      return;
    }
  }
  xm_wpxp_clipping_and_stats(workspace, nzm, nzt, ngrdcol,            // In
                             gr, xm_wpxp_thlm, dt, wp2, thlp2, wm_zt, // In
                             thlm_forcing, rho_ds_zm, rho_ds_zt,      // In
                             invrs_rho_ds_zm, invrs_rho_ds_zt,        // In
                             (thl_tol * thl_tol), thl_tol, rcond,     // In
                             low_lev_effect, high_lev_effect,         // In
                             lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,       // In
                             lhs_diff_zm, C7_Skw_fnc,                 // In
                             lhs_tp, lhs_ta_xm, lhs_pr1_wprtp,        // In
                             l_implemented,
                             Kokkos::subview(solution, Kokkos::ALL, Kokkos::ALL, 1), // In
                             tridiag_solve_method,                                   // In
                             fill_holes_type,                                        // In
                             l_upwind_xm_ma,                                         // In
                             l_enable_relaxed_clipping,                              // In
                             l_mono_flux_lim_thlm,                                   // In
                             l_mono_flux_lim_rtm,                                    // In
                             l_mono_flux_lim_um,                                     // In
                             l_mono_flux_lim_vm,                                     // In
                             l_mono_flux_lim_spikefix,                               // In
                             stats,                                                  // InOut
                             thlm,                                                   // InOut
                             thl_tol_mfl,                                            // In
                             wpthlp, err_info,                                       // InOut
                             &wpthlp_cl_num // Optional InOut
  );
  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << err_info.err_header_global << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "thlm monotonic flux limiter:  tridiag failed" << "\n";
        std::cerr << message.str();
      }
      return;
    }
  }
  // ---> h1g, 2010-06-15
  // scalar transport, e.g, droplet and ice number concentration
  // are handled in  " advance_sclrm_Nd_module.F90 "
  for (sclr = 0; sclr < sclr_dim; sclr += 1) {
    // <--- h1g, 2010-06-15
    xm_wpxp_clipping_and_stats(
        workspace, nzm, nzt, ngrdcol,                                       // In
        gr, xm_wpxp_scalar, dt, wp2,                                        // In
        Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, sclr), wm_zt,     // In
        Kokkos::subview(sclrm_forcing, Kokkos::ALL, Kokkos::ALL, sclr),     // In
        rho_ds_zm, rho_ds_zt,                                               // In
        invrs_rho_ds_zm, invrs_rho_ds_zt,                                   // In
        (sclr_tol_host(sclr) * sclr_tol_host(sclr)),                        // In
        sclr_tol_host(sclr), rcond,                                         // In
        low_lev_effect, high_lev_effect,                                    // In
        lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,                                  // In
        lhs_diff_zm, C7_Skw_fnc,                                            // In
        lhs_tp, lhs_ta_xm, lhs_pr1_wprtp,                                   // In
        l_implemented,                                                      // In
        Kokkos::subview(solution, Kokkos::ALL, Kokkos::ALL, (sclr + 2)),    // In
        tridiag_solve_method,                                               // In
        fill_holes_type,                                                    // In
        l_upwind_xm_ma,                                                     // In
        l_enable_relaxed_clipping,                                          // In
        l_mono_flux_lim_thlm,                                               // In
        l_mono_flux_lim_rtm,                                                // In
        l_mono_flux_lim_um,                                                 // In
        l_mono_flux_lim_vm,                                                 // In
        l_mono_flux_lim_spikefix,                                           // In
        stats,                                                              // InOut
        Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, sclr),             // InOut
        sclr_tol_host(sclr),                                                // In
        Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, sclr), err_info, // InOut
        nullptr                                                             // Optional
    );
    if (clubb_at_least_debug_level_api(0)) {
      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << err_info.err_header_global << "\n";
          std::cerr << message.str();
        }
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "sclrm # " << " " << (sclr + 1) << " "
                  << "monotonic flux limiter: tridiag failed"
                  << "\n";
          std::cerr << message.str();
        }
        return;
      }
    }
  }
  if (l_predict_upwp_vpwp) {
    // Predict <u> and <u'w'>, as well as <v> and <v'w'>.
    xm_wpxp_clipping_and_stats(
        workspace, nzm, nzt, ngrdcol,                                        // In
        gr, xm_wpxp_um, dt, wp2, up2, wm_zt,                                 // In
        um_tndcy, rho_ds_zm, rho_ds_zt,                                      // In
        invrs_rho_ds_zm, invrs_rho_ds_zt,                                    // In
        w_tol_sqd, w_tol, rcond,                                             // In
        low_lev_effect, high_lev_effect,                                     // In
        lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,                                   // In
        lhs_diff_zm, C7_Skw_fnc,                                             // In
        lhs_tp, lhs_ta_xm, lhs_pr1_wprtp,                                    // In
        l_implemented,                                                       // In
        Kokkos::subview(solution, Kokkos::ALL, Kokkos::ALL, (sclr_dim + 2)), // In
        tridiag_solve_method,                                                // In
        fill_holes_type,                                                     // In
        l_upwind_xm_ma,                                                      // In
        l_enable_relaxed_clipping,                                           // In
        l_mono_flux_lim_thlm,                                                // In
        l_mono_flux_lim_rtm,                                                 // In
        l_mono_flux_lim_um,                                                  // In
        l_mono_flux_lim_vm,                                                  // In
        l_mono_flux_lim_spikefix,                                            // In
        stats,                                                               // InOut
        um,                                                                  // InOut
        w_tol,                                                               // In
        upwp, err_info,                                                      // InOut
        &upwp_cl_num                                                         // Optional InOut
    );
    if (clubb_at_least_debug_level_api(0)) {
      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << err_info.err_header_global << "\n";
          std::cerr << message.str();
        }
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "um monotonic flux limiter:  tridiag failed" << "\n";
          std::cerr << message.str();
        }
        return;
      }
    }
    xm_wpxp_clipping_and_stats(
        workspace, nzm, nzt, ngrdcol,                                        // In
        gr, xm_wpxp_vm, dt, wp2, vp2, wm_zt,                                 // In
        vm_tndcy, rho_ds_zm, rho_ds_zt,                                      // In
        invrs_rho_ds_zm, invrs_rho_ds_zt,                                    // In
        w_tol_sqd, w_tol, rcond,                                             // In
        low_lev_effect, high_lev_effect,                                     // In
        lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,                                   // In
        lhs_diff_zm, C7_Skw_fnc,                                             // In
        lhs_tp, lhs_ta_xm, lhs_pr1_wprtp,                                    // In
        l_implemented,                                                       // In
        Kokkos::subview(solution, Kokkos::ALL, Kokkos::ALL, (sclr_dim + 3)), // In
        tridiag_solve_method,                                                // In
        fill_holes_type,                                                     // In
        l_upwind_xm_ma,                                                      // In
        l_enable_relaxed_clipping,                                           // In
        l_mono_flux_lim_thlm,                                                // In
        l_mono_flux_lim_rtm,                                                 // In
        l_mono_flux_lim_um,                                                  // In
        l_mono_flux_lim_vm,                                                  // In
        l_mono_flux_lim_spikefix,                                            // In
        stats,                                                               // InOut
        vm,                                                                  // InOut
        w_tol,                                                               // In
        vpwp, err_info,                                                      // InOut
        &vpwp_cl_num                                                         // Optional InOut
    );
    if (clubb_at_least_debug_level_api(0)) {
      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << err_info.err_header_global << "\n";
          std::cerr << message.str();
        }
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "vm monotonic flux limiter:  tridiag failed" << "\n";
          std::cerr << message.str();
        }
        return;
      }
    }
    if (l_perturbed_wind) {
      xm_wpxp_clipping_and_stats(
          workspace, nzm, nzt, ngrdcol,                                        // In
          gr, xm_wpxp_um, dt, wp2, up2, wm_zt,                                 // In
          um_tndcy, rho_ds_zm, rho_ds_zt,                                      // In
          invrs_rho_ds_zm, invrs_rho_ds_zt,                                    // In
          w_tol_sqd, w_tol, rcond,                                             // In
          low_lev_effect, high_lev_effect,                                     // In
          lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,                                   // In
          lhs_diff_zm, C7_Skw_fnc,                                             // In
          lhs_tp, lhs_ta_xm, lhs_pr1_wprtp,                                    // In
          l_implemented,                                                       // In
          Kokkos::subview(solution, Kokkos::ALL, Kokkos::ALL, (sclr_dim + 4)), // In
          tridiag_solve_method,                                                // In
          fill_holes_type,                                                     // In
          l_upwind_xm_ma,                                                      // In
          l_enable_relaxed_clipping,                                           // In
          l_mono_flux_lim_thlm,                                                // In
          l_mono_flux_lim_rtm,                                                 // In
          l_mono_flux_lim_um,                                                  // In
          l_mono_flux_lim_vm,                                                  // In
          l_mono_flux_lim_spikefix,                                            // In
          stats,                                                               // InOut
          um_pert,                                                             // InOut
          w_tol,                                                               // In
          upwp_pert, err_info,                                                 // InOut
          nullptr                                                              // Optional
      );
      if (clubb_at_least_debug_level_api(0)) {
        if (any_equal(err_info.err_code, clubb_fatal_error)) {
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << err_info.err_header_global << "\n";
            std::cerr << message.str();
          }
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "um_pert monotonic flux limiter:  tridiag failed" << "\n";
            std::cerr << message.str();
          }
          return;
        }
      }
      xm_wpxp_clipping_and_stats(
          workspace, nzm, nzt, ngrdcol,                                        // In
          gr, xm_wpxp_vm, dt, wp2, vp2, wm_zt,                                 // In
          vm_tndcy, rho_ds_zm, rho_ds_zt,                                      // In
          invrs_rho_ds_zm, invrs_rho_ds_zt,                                    // In
          w_tol_sqd, w_tol, rcond,                                             // In
          low_lev_effect, high_lev_effect,                                     // In
          lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,                                   // In
          lhs_diff_zm, C7_Skw_fnc,                                             // In
          lhs_tp, lhs_ta_xm, lhs_pr1_wprtp,                                    // In
          l_implemented,                                                       // In
          Kokkos::subview(solution, Kokkos::ALL, Kokkos::ALL, (sclr_dim + 5)), // In
          tridiag_solve_method,                                                // In
          fill_holes_type,                                                     // In
          l_upwind_xm_ma,                                                      // In
          l_enable_relaxed_clipping,                                           // In
          l_mono_flux_lim_thlm,                                                // In
          l_mono_flux_lim_rtm,                                                 // In
          l_mono_flux_lim_um,                                                  // In
          l_mono_flux_lim_vm,                                                  // In
          l_mono_flux_lim_spikefix,                                            // In
          stats,                                                               // InOut
          vm_pert,                                                             // InOut
          w_tol,                                                               // In
          vpwp_pert, err_info,                                                 // InOut
          nullptr                                                              // Optional
      );
      if (clubb_at_least_debug_level_api(0)) {
        if (any_equal(err_info.err_code, clubb_fatal_error)) {
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << err_info.err_header_global << "\n";
            std::cerr << message.str();
          }
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "vm_pert monotonic flux limiter:  tridiag failed" << "\n";
            std::cerr << message.str();
          }
          return;
        }
      }
    }
  }
}

// nzm: Number of momentum vertical levels
// nzt: Number of thermodynamic vertical levels
// ngrdcol: Number of grid columns
// sclr_dim: Number of passive scalars
// sclr_tol: Threshold(s) on the passive scalars  [units vary]
// dt: Timestep                                 [s]
// nrhs: Number of RHS vectors
// wm_zt: w wind component on thermodynamic levels [m/s]
// wp2: w'^2 (momentum levels)                   [m^2/s^2]
// rtpthvp: r_t'th_v' (momentum levels)              [(kg/kg) K]
// rtm_forcing: r_t forcing (thermodynamic levels)       [(kg/kg)/s]
// wprtp_forcing: <w'r_t'> forcing (momentum levels)       [(kg/kg)/s^2]
// thlpthvp: th_l'th_v' (momentum levels)             [K^2]
// thlm_forcing: th_l forcing (thermodynamic levels)      [K/s]
// wpthlp_forcing: <w'th_l'> forcing (momentum levels)      [K/s^2]
// rho_ds_zm: Dry, static density on momentum levels   [kg/m^3]
// rho_ds_zt: Dry, static density on thermo. levels    [kg/m^3]
// invrs_rho_ds_zm: Inv. dry, static density @ moment. levs. [m^3/kg]
// invrs_rho_ds_zt: Inv. dry, static density @ thermo. levs. [m^3/kg]
// thv_ds_zm: Dry, base-state theta_v on moment. levs. [K]
// rtp2: r_t'^2 (momentum levels)                 [(kg/kg)^2]
// thlp2: th_l'^2 (momentum levels)                [K^2]
// l_implemented: True if CLUBB is being implemented and run in a host model
// sclrpthvp: <sclr' th_v'> (momentum levels)       [Units vary]
// sclrm_forcing: sclrm forcing (thermodynamic levels)  [Units vary]
// sclrp2: For clipping Vince Larson             [Units vary]
// low_lev_effect: Index of the lowest level that has an effect.
// high_lev_effect: Index of the highest level that has an effect.
// lhs_diff_zm: Diffusion term for w'x'
// lhs_diff_zt: Diffusion term for w'x'
// lhs_ma_zt: Mean advection contributions to lhs
// lhs_ma_zm: Mean advection contributions to lhs
// lhs_ta_wprtp: w'r_t' turbulent advection contributions to lhs
// lhs_ta_wpthlp: w'thl' turbulent advection contributions to lhs
// lhs_ta_wpsclrp: w'sclr' turbulent advection contributions to lhs
// rhs_ta_wprtp: w'r_t' turbulent advection contributions to rhs
// rhs_ta_wpthlp: w'thl' turbulent advection contributions to rhs
// rhs_ta_wpsclrp: w'sclr' turbulent advection contributions to rhs
// lhs_tp: Turbulent production terms of w'x'
// lhs_ta_xm: Turbulent advection terms of xm
// lhs_ac_pr2: Accumulation of w'x' and w'x' pressure term 2
// lhs_pr1_wprtp: Pressure term 1 for w'r_t' for all grid levels
// lhs_pr1_wpthlp: Pressure term 1 for w'thl' for all grid levels
// lhs_pr1_wpsclrp: Pressure term 1 for w'sclr' for all grid levels
// penta_solve_method: Method to solve then penta-diagonal system
// tridiag_solve_method: Specifier for method to solve tridiagonal systems
// fill_holes_type: Specifier for which hole filling method to use
// l_diffuse_rtm_and_thlm: This flag determines whether or not we want CLUBB to do diffusion on
// rtm and thlm
// l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent or mean advection terms. It
// affects rtm, thlm, sclrm, um and vm.
// l_enable_relaxed_clipping: Flag to relax clipping on wpxp in xm_wpxp_clipping_and_stats
// l_mono_flux_lim_thlm: Flag to turn on monotonic flux limiter for thlm
// l_mono_flux_lim_rtm: Flag to turn on monotonic flux limiter for rtm
// l_mono_flux_lim_um: Flag to turn on monotonic flux limiter for um
// l_mono_flux_lim_vm: Flag to turn on monotonic flux limiter for vm
// l_mono_flux_lim_spikefix: Flag to implement monotonic flux limiter code that eliminates
// spurious drying tendencies at model top
// rtm: r_t  (total water mixing ratio)           [kg/kg]
// wprtp: w'r_t'                                    [(kg/kg) m/s]
// thlm: th_l (liquid water potential temperature) [K]
// wpthlp: w'th_l'                                   [K m/s]
// sclrm: [Units vary]
// wpsclrp: [Units vary]
// err_info: err_info struct containing err_code and err_header
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26, class Array27, class Array28, class Array29,
          class Array30, class Array31, class Array32, class Array33, class Array34, class Array35,
          class Array36, class Array37, class Array38, class Array39, class Array40, class Array41,
          class Array42, class Array43>
void solve_xm_wpxp_with_multiple_lhs(
    ScratchWorkspace& workspace,                                                    // InOut
    int nzm, int nzt, int ngrdcol, int sclr_dim, Array0 sclr_tol, const grid& gr,   // In
    real_t dt,                                                                      // In
    bool l_iter, int nrhs, Array1 wm_zt, Array2 wp2,                                // In
    Array3 rtpthvp, Array4 rtm_forcing, Array5 wprtp_forcing, Array6 thlpthvp,      // In
    Array7 thlm_forcing, Array8 wpthlp_forcing, Array9 rho_ds_zm,                   // In
    Array10 rho_ds_zt, Array11 invrs_rho_ds_zm, Array12 invrs_rho_ds_zt,            // In
    Array13 thv_ds_zm, Array14 rtp2, Array15 thlp2, bool l_implemented,             // In
    Array16 sclrpthvp, Array17 sclrm_forcing, Array18 sclrp2,                       // In
    Array19 low_lev_effect, Array20 high_lev_effect, Array21 C7_Skw_fnc,            // In
    Array22 lhs_diff_zm, Array23 lhs_diff_zt, Array24 lhs_ma_zt, Array25 lhs_ma_zm, // In
    Array26 lhs_ta_wprtp, Array27 lhs_ta_wpthlp, Array28 lhs_ta_wpsclrp,            // In
    Array29 rhs_ta_wprtp, Array30 rhs_ta_wpthlp, Array31 rhs_ta_wpsclrp,            // In
    Array32 lhs_tp, Array33 lhs_ta_xm, Array34 lhs_ac_pr2, Array35 lhs_pr1_wprtp,   // In
    Array36 lhs_pr1_wpthlp, Array37 lhs_pr1_wpsclrp,                                // In
    int penta_solve_method,                                                         // In
    int tridiag_solve_method,                                                       // In
    int fill_holes_type,                                                            // In
    bool l_diffuse_rtm_and_thlm,                                                    // In
    bool l_upwind_xm_ma,                                                            // In
    bool l_enable_relaxed_clipping,                                                 // In
    bool l_mono_flux_lim_thlm,                                                      // In
    bool l_mono_flux_lim_rtm,                                                       // In
    bool l_mono_flux_lim_um,                                                        // In
    bool l_mono_flux_lim_vm,                                                        // In
    bool l_mono_flux_lim_spikefix,                                                  // In
    int& wprtp_cl_num, int& wpthlp_cl_num,                                          // InOut
    stats_type& stats,                                                              // InOut
    Array38 rtm, Array39 wprtp, Array40 thlm, Array41 wpthlp, Array42 sclrm,        // InOut
    Array43 wpsclrp, err_info_type& err_info                                        // InOut
) {
  //
  // Description: This subroutine solves all xm_wpxp when all the LHS matrices are NOT equal.
  //              This means multiple solves are required, one for each unique LHS.
  //
  //----------------------------------------------------------------------------------------
  // ------------------- Input Variables -------------------
  // Additional variables for passive scalars
  // LHS/RHS terms
  // Variables used as part of the monotonic turbulent advection scheme.
  // Find the lowermost and uppermost grid levels that can have an effect
  // on the central thermodynamic level during the course of a time step,
  // due to the effects of turbulent advection only.
  // CLUBB to do diffusion on rtm and thlm
  // an upwind differencing approximation rather
  // than a centered differencing for turbulent
  // or mean advection terms. It affects rtm,
  // thlm, sclrm, um and vm.
  // i.e. TKE = 1/2 (u'^2 + v'^2 + w'^2)
  // xm_wpxp_clipping_and_stats
  // eliminates spurious drying tendencies at model top
  // ------------------- Input/Output Variables -------------------
  // ------------------- Local Variables -------------------
  // Additional variables for passive scalars
  // Constant parameters as a function of Skw.
  auto workspace_scope = workspace.scope();
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_multiple_lhs::lhs", nsub + nsup + 1,
          ngrdcol, 2 * nzm - 1); // Implicit contributions to wpxp/xm (band diag. matrix) (LAPACK)
  Kokkos::View<real_t***, Kokkos::LayoutLeft> rhs =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_multiple_lhs::rhs", ngrdcol, 2 * nzm - 1,
          nrhs); // Right-hand sides of band diag. matrix. (LAPACK)
  Kokkos::View<real_t***, Kokkos::LayoutLeft> rhs_save =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_multiple_lhs::rhs_save", ngrdcol,
          2 * nzm - 1, nrhs); // Saved Right-hand sides of band diag. matrix. (LAPACK)
  Kokkos::View<real_t***, Kokkos::LayoutLeft> solution =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_multiple_lhs::solution", ngrdcol,
          2 * nzm - 1, nrhs); // solution vectors of band diag. matrix. (LAPACK)
  Kokkos::View<real_t***, Kokkos::LayoutLeft> old_solution =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_multiple_lhs::old_solution", ngrdcol,
          2 * nzm - 1, nrhs); // previous solutions
  Kokkos::View<real_t***, Kokkos::LayoutLeft> wpsclrp_forcing =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_multiple_lhs::wpsclrp_forcing", ngrdcol,
          nzm, sclr_dim); // <w'sclr'> forcing (momentum levels)  [m/s{un vary}]
  Kokkos::View<real_t*, Kokkos::LayoutLeft> rcond =
      workspace_scope.view<Kokkos::View<real_t*, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "solve_xm_wpxp_with_multiple_lhs::rcond", ngrdcol);
  int i;
  int sclr;
  int k;
  // ------------------- Begin Code -------------------
  // Compute the implicit portion of the r_t and w'r_t' equations.
  // Build the left-hand side matrix.
  xm_wpxp_lhs(nzm, nzt, ngrdcol, l_iter, dt,                 // In
              l_implemented, lhs_diff_zm, lhs_diff_zt,       // In
              lhs_ma_zm, lhs_ma_zt, lhs_ta_wprtp, lhs_ta_xm, // In
              lhs_tp, lhs_pr1_wprtp, lhs_ac_pr2,             // In
              l_diffuse_rtm_and_thlm,                        // In
              lhs                                            // Out
  );
  // Compute the explicit portion of the r_t and w'r_t' equations.
  // Build the right-hand side vector.
  xm_wpxp_rhs(workspace, nzm, nzt, ngrdcol, gr, xm_wpxp_rtm, l_iter, dt, // In
              rtm, wprtp, rtm_forcing, wprtp_forcing, C7_Skw_fnc,        // In
              rtpthvp, rhs_ta_wprtp, thv_ds_zm,                          // In
              lhs_pr1_wprtp, lhs_ta_wprtp,                               // In
              stats,                                                     // InOut
              Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, 0)          // Out
  );
  // Save the value of rhs, which will be overwritten with the solution as
  // part of the solving routine.
  Kokkos::parallel_for(
      "solve_xm_wpxp_with_multiple_lhs::array1",
      rhs_save.extent(0) * rhs_save.extent(1) * rhs_save.extent(2), KOKKOS_LAMBDA(int cell) {
        const int slice_0                   = (cell % rhs_save.extent(0));
        const int slice_1                   = ((cell / rhs_save.extent(0)) % rhs_save.extent(1));
        const int slice_2                   = (cell / (rhs_save.extent(0) * rhs_save.extent(1)));
        rhs_save(slice_0, slice_1, slice_2) = rhs(slice_0, slice_1, slice_2);
      });
  // Use the previous solution as an initial guess for the bicgstab method
  if (penta_solve_method == penta_bicgstab) {
    for (k = 0; k < nzt; k += 1) {
      Kokkos::parallel_for(
          "solve_xm_wpxp_with_multiple_lhs::array2", old_solution.extent(0),
          KOKKOS_LAMBDA(int cell) {
            const int slice_0                     = cell;
            old_solution(slice_0, (2 * k + 1), 0) = rtm(slice_0, k);
          });
    }
    for (k = 0; k < nzm; k += 1) {
      Kokkos::parallel_for(
          "solve_xm_wpxp_with_multiple_lhs::array3", old_solution.extent(0),
          KOKKOS_LAMBDA(int cell) {
            const int slice_0                 = cell;
            old_solution(slice_0, (2 * k), 0) = wprtp(slice_0, k);
          });
    }
  }
  // Solve r_t / w'r_t'
  if (stats.l_sample && var_on_stats_list(stats, "rtm_matrix_condt_num")) {
    xm_wpxp_solve(workspace, nzm, ngrdcol, nrhs, gr, // In
                  penta_solve_method,                // In
                  l_implemented,                     // In
                  old_solution,                      // InOut
                  lhs, rhs, err_info,                // InOut
                  solution,                          // Out
                  &rcond                             // Optional Out
    );
  } else {
    xm_wpxp_solve(workspace, nzm, ngrdcol, nrhs, gr, // In
                  penta_solve_method,                // In
                  l_implemented,                     // In
                  old_solution,                      // InOut
                  lhs, rhs, err_info,                // InOut
                  solution,                          // Out
                  nullptr                            // Optional
    );
  }
  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << err_info.err_header_global << "\n";
        std::cerr << message.str();
      }
      for (i = 0; i < ngrdcol; i += 1) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "Mean total water & total water flux LU decomp. failed" << "\n";
          std::cerr << message.str();
        }
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "rtm and wprtp LHS" << "\n";
          std::cerr << message.str();
        }
        for (k = 0; k < nzm - 1; k += 1) {
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "grid col = " << " " << (i + 1) << " " << "zm level = " << " " << (k + 1)
                    << " "
                    << "height [m] = " << " " << to_host_space(workspace, gr.zm)(i, k) << " "
                    << "LHS = " << " "
                    << dump_values(
                           Kokkos::subview(lhs, std::make_pair(0, nsub + nsup + 1), i, (2 * k)))
                    << "\n";
            std::cerr << message.str();
          }
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "grid col = " << " " << (i + 1) << " " << "zt level = " << " " << (k + 1)
                    << " "
                    << "height [m] = " << " " << to_host_space(workspace, gr.zt)(i, k) << " "
                    << "LHS = " << " "
                    << dump_values(
                           Kokkos::subview(lhs, std::make_pair(0, nsub + nsup + 1), i, (2 * k + 1)))
                    << "\n";
            std::cerr << message.str();
          }
        }
        k = nzm - 1;
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "grid col = " << " " << (i + 1) << " " << "zm level = " << " " << (k + 1)
                  << " "
                  << "height [m] = " << " " << to_host_space(workspace, gr.zm)(i, k) << " "
                  << "LHS = " << " "
                  << dump_values(
                         Kokkos::subview(lhs, std::make_pair(0, nsub + nsup + 1), i, (2 * k)))
                  << "\n";
          std::cerr << message.str();
        }
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "rtm and wprtp RHS" << "\n";
          std::cerr << message.str();
        }
        for (k = 0; k < nzm - 1; k += 1) {
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "grid col = " << " " << (i + 1) << " " << "zm level = " << " " << (k + 1)
                    << " "
                    << "height [m] = " << " " << to_host_space(workspace, gr.zm)(i, k) << " "
                    << "RHS = " << " " << to_host_space(workspace, rhs_save)(i, 2 * k, 0) << "\n";
            std::cerr << message.str();
          }
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "grid col = " << " " << (i + 1) << " " << "zt level = " << " " << (k + 1)
                    << " "
                    << "height [m] = " << " " << to_host_space(workspace, gr.zt)(i, k) << " "
                    << "RHS = " << " " << to_host_space(workspace, rhs_save)(i, 2 * k + 1, 0)
                    << "\n";
            std::cerr << message.str();
          }
        }
        k = nzm - 1;
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "grid col = " << " " << (i + 1) << " " << "zm level = " << " " << (k + 1)
                  << " "
                  << "height [m] = " << " " << to_host_space(workspace, gr.zm)(i, k) << " "
                  << "RHS = " << " " << to_host_space(workspace, rhs_save)(i, 2 * k, 0) << "\n";
          std::cerr << message.str();
        }
      }
      return;
    }
  }
  xm_wpxp_clipping_and_stats(workspace, nzm, nzt, ngrdcol,          // In
                             gr, xm_wpxp_rtm, dt, wp2, rtp2, wm_zt, // In
                             rtm_forcing, rho_ds_zm, rho_ds_zt,     // In
                             invrs_rho_ds_zm, invrs_rho_ds_zt,      // In
                             (rt_tol * rt_tol), rt_tol, rcond,      // In
                             low_lev_effect, high_lev_effect,       // In
                             lhs_ma_zt, lhs_ma_zm, lhs_ta_wprtp,    // In
                             lhs_diff_zm, C7_Skw_fnc,               // In
                             lhs_tp, lhs_ta_xm, lhs_pr1_wprtp,      // In
                             l_implemented,
                             Kokkos::subview(solution, Kokkos::ALL, Kokkos::ALL, 0), // In
                             tridiag_solve_method,                                   // In
                             fill_holes_type,                                        // In
                             l_upwind_xm_ma,                                         // In
                             l_enable_relaxed_clipping,                              // In
                             l_mono_flux_lim_thlm,                                   // In
                             l_mono_flux_lim_rtm,                                    // In
                             l_mono_flux_lim_um,                                     // In
                             l_mono_flux_lim_vm,                                     // In
                             l_mono_flux_lim_spikefix,                               // In
                             stats,                                                  // InOut
                             rtm,                                                    // InOut
                             rt_tol_mfl,                                             // In
                             wprtp, err_info,                                        // InOut
                             &wprtp_cl_num // Optional InOut
  );
  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << err_info.err_header_global << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "rtm monotonic flux limiter:  tridiag failed" << "\n";
        std::cerr << message.str();
      }
      return;
    }
  }
  // Compute the implicit portion of the th_l and w'th_l' equations.
  // Build the left-hand side matrix.
  xm_wpxp_lhs(nzm, nzt, ngrdcol, l_iter, dt,                  // In
              l_implemented, lhs_diff_zm, lhs_diff_zt,        // In
              lhs_ma_zm, lhs_ma_zt, lhs_ta_wpthlp, lhs_ta_xm, // In
              lhs_tp, lhs_pr1_wpthlp, lhs_ac_pr2,             // In
              l_diffuse_rtm_and_thlm,                         // In
              lhs                                             // Out
  );
  // Compute the explicit portion of the th_l and w'th_l' equations.
  // Build the right-hand side vector.
  xm_wpxp_rhs(workspace, nzm, nzt, ngrdcol, gr, xm_wpxp_thlm, l_iter, dt, // In
              thlm, wpthlp, thlm_forcing, wpthlp_forcing, C7_Skw_fnc,     // In
              thlpthvp, rhs_ta_wpthlp, thv_ds_zm,                         // In
              lhs_pr1_wpthlp, lhs_ta_wpthlp,                              // In
              stats,                                                      // InOut
              Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, 0)           // Out
  );
  // Save the value of rhs, which will be overwritten with the solution as
  // part of the solving routine.
  Kokkos::parallel_for(
      "solve_xm_wpxp_with_multiple_lhs::array4",
      rhs_save.extent(0) * rhs_save.extent(1) * rhs_save.extent(2), KOKKOS_LAMBDA(int cell) {
        const int slice_0                   = (cell % rhs_save.extent(0));
        const int slice_1                   = ((cell / rhs_save.extent(0)) % rhs_save.extent(1));
        const int slice_2                   = (cell / (rhs_save.extent(0) * rhs_save.extent(1)));
        rhs_save(slice_0, slice_1, slice_2) = rhs(slice_0, slice_1, slice_2);
      });
  // Use the previous solution as an initial guess for the bicgstab method
  if (penta_solve_method == penta_bicgstab) {
    for (k = 0; k < nzt; k += 1) {
      Kokkos::parallel_for(
          "solve_xm_wpxp_with_multiple_lhs::array5", old_solution.extent(0),
          KOKKOS_LAMBDA(int cell) {
            const int slice_0                     = cell;
            old_solution(slice_0, (2 * k + 1), 0) = thlm(slice_0, k);
          });
    }
    for (k = 0; k < nzm; k += 1) {
      Kokkos::parallel_for(
          "solve_xm_wpxp_with_multiple_lhs::array6", old_solution.extent(0),
          KOKKOS_LAMBDA(int cell) {
            const int slice_0                 = cell;
            old_solution(slice_0, (2 * k), 0) = wpthlp(slice_0, k);
          });
    }
  }
  // Solve for th_l / w'th_l'
  if (stats.l_sample && var_on_stats_list(stats, "thlm_matrix_condt_num")) {
    xm_wpxp_solve(workspace, nzm, ngrdcol, nrhs, gr, // In
                  penta_solve_method,                // In
                  l_implemented,                     // In
                  old_solution,                      // InOut
                  lhs, rhs, err_info,                // InOut
                  solution,                          // Out
                  &rcond                             // Optional Out
    );
  } else {
    xm_wpxp_solve(workspace, nzm, ngrdcol, nrhs, gr, // In
                  penta_solve_method,                // In
                  l_implemented,                     // In
                  old_solution,                      // InOut
                  lhs, rhs, err_info,                // InOut
                  solution,                          // Out
                  nullptr                            // Optional
    );
  }
  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << err_info.err_header_global << "\n";
        std::cerr << message.str();
      }
      for (i = 0; i < ngrdcol; i += 1) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "Liquid pot. temp & thetal flux LU decomp. failed" << "\n";
          std::cerr << message.str();
        }
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "thlm and wpthlp LHS" << "\n";
          std::cerr << message.str();
        }
        for (k = 0; k < nzm - 1; k += 1) {
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "grid col = " << " " << (i + 1) << " " << "zm level = " << " " << (k + 1)
                    << " "
                    << "height [m] = " << " " << to_host_space(workspace, gr.zm)(i, k) << " "
                    << "LHS = " << " "
                    << dump_values(
                           Kokkos::subview(lhs, std::make_pair(0, nsub + nsup + 1), i, (2 * k)))
                    << "\n";
            std::cerr << message.str();
          }
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "grid col = " << " " << (i + 1) << " " << "zt level = " << " " << (k + 1)
                    << " "
                    << "height [m] = " << " " << to_host_space(workspace, gr.zt)(i, k) << " "
                    << "LHS = " << " "
                    << dump_values(
                           Kokkos::subview(lhs, std::make_pair(0, nsub + nsup + 1), i, (2 * k + 1)))
                    << "\n";
            std::cerr << message.str();
          }
        }
        k = nzm - 1;
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "grid col = " << " " << (i + 1) << " " << "zm level = " << " " << (k + 1)
                  << " "
                  << "height [m] = " << " " << to_host_space(workspace, gr.zm)(i, k) << " "
                  << "LHS = " << " "
                  << dump_values(
                         Kokkos::subview(lhs, std::make_pair(0, nsub + nsup + 1), i, (2 * k)))
                  << "\n";
          std::cerr << message.str();
        }
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "thlm and wpthlp RHS" << "\n";
          std::cerr << message.str();
        }
        for (k = 0; k < nzm - 1; k += 1) {
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "grid col = " << " " << (i + 1) << " " << "zm level = " << " " << (k + 1)
                    << " "
                    << "height [m] = " << " " << to_host_space(workspace, gr.zm)(i, k) << " "
                    << "RHS = " << " " << to_host_space(workspace, rhs_save)(i, 2 * k, 0) << "\n";
            std::cerr << message.str();
          }
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "grid col = " << " " << (i + 1) << " " << "zt level = " << " " << (k + 1)
                    << " "
                    << "height [m] = " << " " << to_host_space(workspace, gr.zt)(i, k) << " "
                    << "RHS = " << " " << to_host_space(workspace, rhs_save)(i, 2 * k + 1, 0)
                    << "\n";
            std::cerr << message.str();
          }
        }
        k = nzm - 1;
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "grid col = " << " " << (i + 1) << " " << "zm level = " << " " << (k + 1)
                  << " "
                  << "height [m] = " << " " << to_host_space(workspace, gr.zm)(i, k) << " "
                  << "RHS = " << " " << to_host_space(workspace, rhs_save)(i, 2 * k, 0) << "\n";
          std::cerr << message.str();
        }
      }
      return;
    }
  }
  xm_wpxp_clipping_and_stats(workspace, nzm, nzt, ngrdcol,            // In
                             gr, xm_wpxp_thlm, dt, wp2, thlp2, wm_zt, // In
                             thlm_forcing, rho_ds_zm, rho_ds_zt,      // In
                             invrs_rho_ds_zm, invrs_rho_ds_zt,        // In
                             (thl_tol * thl_tol), thl_tol, rcond,     // In
                             low_lev_effect, high_lev_effect,         // In
                             lhs_ma_zt, lhs_ma_zm, lhs_ta_wpthlp,     // In
                             lhs_diff_zm, C7_Skw_fnc,                 // In
                             lhs_tp, lhs_ta_xm, lhs_pr1_wpthlp,       // In
                             l_implemented,
                             Kokkos::subview(solution, Kokkos::ALL, Kokkos::ALL, 0), // In
                             tridiag_solve_method,                                   // In
                             fill_holes_type,                                        // In
                             l_upwind_xm_ma,                                         // In
                             l_enable_relaxed_clipping,                              // In
                             l_mono_flux_lim_thlm,                                   // In
                             l_mono_flux_lim_rtm,                                    // In
                             l_mono_flux_lim_um,                                     // In
                             l_mono_flux_lim_vm,                                     // In
                             l_mono_flux_lim_spikefix,                               // In
                             stats,                                                  // InOut
                             thlm,                                                   // InOut
                             thl_tol_mfl,                                            // In
                             wpthlp, err_info,                                       // InOut
                             &wpthlp_cl_num // Optional InOut
  );
  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << err_info.err_header_global << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "thlm monotonic flux limiter:  tridiag failed" << "\n";
        std::cerr << message.str();
      }
      return;
    }
  }
  // Solve sclrm / wpsclrp
  // If sclr_dim is 0, then this loop will execute 0 times.
  // ---> h1g, 2010-06-15
  // scalar transport, e.g, droplet and ice number concentration
  // are handled in  " advance_sclrm_Nd_module.F90 "
  for (sclr = 0; sclr < sclr_dim; sclr += 1) {
    // <--- h1g, 2010-06-15
    // Set <w'sclr'> forcing to 0 unless unless testing the wpsclrp code
    // using wprtp or wpthlp (then use wprtp_forcing or wpthlp_forcing).
    Kokkos::parallel_for(
        "solve_xm_wpxp_with_multiple_lhs::array7",
        wpsclrp_forcing.extent(0) * wpsclrp_forcing.extent(1), KOKKOS_LAMBDA(int cell) {
          const int slice_0                       = (cell % wpsclrp_forcing.extent(0));
          const int slice_1                       = (cell / wpsclrp_forcing.extent(0));
          wpsclrp_forcing(slice_0, slice_1, sclr) = zero;
        });
    // Compute the implicit portion of the sclr and w'sclr' equations.
    // Build the left-hand side matrix.
    xm_wpxp_lhs(nzm, nzt, ngrdcol, l_iter, dt,                                                // In
                l_implemented, lhs_diff_zm, lhs_diff_zt,                                      // In
                lhs_ma_zm, lhs_ma_zt,                                                         // In
                Kokkos::subview(lhs_ta_wpsclrp, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, sclr), // In
                lhs_ta_xm,                                                                    // In
                lhs_tp, lhs_pr1_wpsclrp, lhs_ac_pr2,                                          // In
                l_diffuse_rtm_and_thlm,                                                       // In
                lhs                                                                           // Out
    );
    // Compute the explicit portion of the sclrm and w'sclr' equations.
    // Build the right-hand side vector.
    xm_wpxp_rhs(workspace, nzm, nzt, ngrdcol, gr, xm_wpxp_scalar, l_iter, dt,                 // In
                Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, sclr),                       // In
                Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, sclr),                     // In
                Kokkos::subview(sclrm_forcing, Kokkos::ALL, Kokkos::ALL, sclr),               // In
                Kokkos::subview(wpsclrp_forcing, Kokkos::ALL, Kokkos::ALL, sclr),             // In
                C7_Skw_fnc,                                                                   // In
                Kokkos::subview(sclrpthvp, Kokkos::ALL, Kokkos::ALL, sclr),                   // In
                Kokkos::subview(rhs_ta_wpsclrp, Kokkos::ALL, Kokkos::ALL, sclr),              // In
                thv_ds_zm,                                                                    // In
                lhs_pr1_wpsclrp,                                                              // In
                Kokkos::subview(lhs_ta_wpsclrp, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, sclr), // In
                stats,                                            // InOut
                Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, 0) // Out
    );
    // Save the value of rhs, which will be overwritten with the solution as
    // part of the solving routine.
    Kokkos::parallel_for(
        "solve_xm_wpxp_with_multiple_lhs::array8",
        rhs_save.extent(0) * rhs_save.extent(1) * rhs_save.extent(2), KOKKOS_LAMBDA(int cell) {
          const int slice_0                   = (cell % rhs_save.extent(0));
          const int slice_1                   = ((cell / rhs_save.extent(0)) % rhs_save.extent(1));
          const int slice_2                   = (cell / (rhs_save.extent(0) * rhs_save.extent(1)));
          rhs_save(slice_0, slice_1, slice_2) = rhs(slice_0, slice_1, slice_2);
        });
    // Use the previous solution as an initial guess for the bicgstab method
    if (penta_solve_method == penta_bicgstab) {
      for (k = 0; k < nzt; k += 1) {
        Kokkos::parallel_for(
            "solve_xm_wpxp_with_multiple_lhs::array9", old_solution.extent(0),
            KOKKOS_LAMBDA(int cell) {
              const int slice_0                     = cell;
              old_solution(slice_0, (2 * k + 1), 0) = sclrm(slice_0, k, sclr);
            });
      }
      for (k = 0; k < nzm; k += 1) {
        Kokkos::parallel_for(
            "solve_xm_wpxp_with_multiple_lhs::array10", old_solution.extent(0),
            KOKKOS_LAMBDA(int cell) {
              const int slice_0                 = cell;
              old_solution(slice_0, (2 * k), 0) = wpsclrp(slice_0, k, sclr);
            });
      }
    }
    // Solve for sclrm / w'sclr'
    xm_wpxp_solve(workspace, nzm, ngrdcol, nrhs, gr, // In
                  penta_solve_method,                // In
                  l_implemented,                     // In
                  old_solution,                      // InOut
                  lhs, rhs, err_info,                // InOut
                  solution,                          // Out
                  nullptr                            // Optional
    );
    if (clubb_at_least_debug_level_api(0)) {
      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << err_info.err_header_global << "\n";
          std::cerr << message.str();
        }
        for (i = 0; i < ngrdcol; i += 1) {
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "Passive scalar # " << " " << (sclr + 1) << " " << " LU decomp. failed."
                    << "\n";
            std::cerr << message.str();
          }
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "sclrm and wpsclrp LHS" << "\n";
            std::cerr << message.str();
          }
          for (k = 0; k < nzm - 1; k += 1) {
            {
              std::ostringstream message;
              message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
              message << "grid col = " << " " << (i + 1) << " " << "zm level = " << " " << (k + 1)
                      << " "
                      << "height [m] = " << " " << to_host_space(workspace, gr.zm)(i, k) << " "
                      << "LHS = " << " "
                      << dump_values(
                             Kokkos::subview(lhs, std::make_pair(0, nsub + nsup + 1), i, (2 * k)))
                      << "\n";
              std::cerr << message.str();
            }
            {
              std::ostringstream message;
              message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
              message << "grid col = " << " " << (i + 1) << " " << "zt level = " << " " << (k + 1)
                      << " "
                      << "height [m] = " << " " << to_host_space(workspace, gr.zt)(i, k) << " "
                      << "LHS = " << " "
                      << dump_values(Kokkos::subview(lhs, std::make_pair(0, nsub + nsup + 1), i,
                                                     (2 * k + 1)))
                      << "\n";
              std::cerr << message.str();
            }
          }
          k = nzm - 1;
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "grid col = " << " " << (i + 1) << " " << "zm level = " << " " << (k + 1)
                    << " "
                    << "height [m] = " << " " << to_host_space(workspace, gr.zm)(i, k) << " "
                    << "LHS = " << " "
                    << dump_values(
                           Kokkos::subview(lhs, std::make_pair(0, nsub + nsup + 1), i, (2 * k)))
                    << "\n";
            std::cerr << message.str();
          }
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "sclrm and wpsclrp RHS" << "\n";
            std::cerr << message.str();
          }
          for (k = 0; k < nzm - 1; k += 1) {
            {
              std::ostringstream message;
              message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
              message << "grid col = " << " " << (i + 1) << " " << "zm level = " << " " << (k + 1)
                      << " "
                      << "height [m] = " << " " << to_host_space(workspace, gr.zm)(i, k) << " "
                      << "RHS = " << " " << to_host_space(workspace, rhs_save)(i, 2 * k, 0) << "\n";
              std::cerr << message.str();
            }
            {
              std::ostringstream message;
              message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
              message << "grid col = " << " " << (i + 1) << " " << "zt level = " << " " << (k + 1)
                      << " "
                      << "height [m] = " << " " << to_host_space(workspace, gr.zt)(i, k) << " "
                      << "RHS = " << " " << to_host_space(workspace, rhs_save)(i, 2 * k + 1, 0)
                      << "\n";
              std::cerr << message.str();
            }
          }
          k = nzm - 1;
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "grid col = " << " " << (i + 1) << " " << "zm level = " << " " << (k + 1)
                    << " "
                    << "height [m] = " << " " << to_host_space(workspace, gr.zm)(i, k) << " "
                    << "RHS = " << " " << to_host_space(workspace, rhs_save)(i, 2 * k, 0) << "\n";
            std::cerr << message.str();
          }
        }
        return;
      }
    }
    xm_wpxp_clipping_and_stats(
        workspace, nzm, nzt, ngrdcol,                                                          // In
        gr, xm_wpxp_scalar, dt, wp2,                                                           // In
        Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, sclr), wm_zt,                        // In
        Kokkos::subview(sclrm_forcing, Kokkos::ALL, Kokkos::ALL, sclr),                        // In
        rho_ds_zm, rho_ds_zt,                                                                  // In
        invrs_rho_ds_zm, invrs_rho_ds_zt,                                                      // In
        (to_host_space(workspace, sclr_tol)(sclr) * to_host_space(workspace, sclr_tol)(sclr)), // In
        to_host_space(workspace, sclr_tol)(sclr), rcond,                                       // In
        low_lev_effect, high_lev_effect,                                                       // In
        lhs_ma_zt, lhs_ma_zm,                                                                  // In
        Kokkos::subview(lhs_ta_wpsclrp, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, sclr),          // In
        lhs_diff_zm, C7_Skw_fnc,                                                               // In
        lhs_tp, lhs_ta_xm, lhs_pr1_wpsclrp,                                                    // In
        l_implemented,                                                                         // In
        Kokkos::subview(solution, Kokkos::ALL, Kokkos::ALL, 0),                                // In
        tridiag_solve_method,                                                                  // In
        fill_holes_type,                                                                       // In
        l_upwind_xm_ma,                                                                        // In
        l_enable_relaxed_clipping,                                                             // In
        l_mono_flux_lim_thlm,                                                                  // In
        l_mono_flux_lim_rtm,                                                                   // In
        l_mono_flux_lim_um,                                                                    // In
        l_mono_flux_lim_vm,                                                                    // In
        l_mono_flux_lim_spikefix,                                                              // In
        stats,                                                              // InOut
        Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, sclr),             // InOut
        to_host_space(workspace, sclr_tol)(sclr),                           // In
        Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, sclr), err_info, // InOut
        nullptr                                                             // Optional
    );
    if (clubb_at_least_debug_level_api(0)) {
      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << err_info.err_header_global << "\n";
          std::cerr << message.str();
        }
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "sclrm # " << " " << (sclr + 1) << " "
                  << "monotonic flux limiter: tridiag failed"
                  << "\n";
          std::cerr << message.str();
        }
        return;
      }
    }
  }
}

// nrhs: Number of rhs vectors
// penta_solve_method: Method to solve then penta-diagonal system
// l_implemented: True if CLUBB is being implemented and run in a host model
// lhs: Implicit contributions to wpxp/xm (band diag. matrix in LAPACK storage)
// rhs: Right-hand side of band diag. matrix. (LAPACK storage)
// err_info: err_info struct containing err_code and err_header
// solution: Solution to band diagonal system (LAPACK storage)
// rcond: Est. of the reciprocal of the condition #
template <class Array0, class Array1, class Array2, class Array3>
void xm_wpxp_solve(ScratchWorkspace& workspace,                     // InOut
                   int nzm, int ngrdcol, int nrhs, const grid& gr,  // In
                   int penta_solve_method,                          // In
                   bool l_implemented,                              // In
                   Array0 old_solution,                             // InOut
                   Array1 lhs, Array2 rhs, err_info_type& err_info, // InOut
                   Array3 solution,                                 // Out
                   Real1D* rcond                                    // Optional Out
) {
  // Description:
  //   Solve for xm / w'x' using the band diagonal solver.
  // References:
  //   None
  //------------------------------------------------------------------------
  //------------------------- Input Variables -------------------------
  //------------------------- Input/Output Variables -------------------------
  //------------------------- Output Variables -------------------------
  // Solve for xm / w'x' using the native interleaved matrix equations.
  // lhs is Native (band,column,row); matrix_solver_wrapper owns
  // (column,band,row), so transpose the first two axes only at this boundary.
  auto workspace_scope = workspace.scope();
  const bool reverse   = l_force_descending_solves && gr.grid_dir_indx > 0;
  Real3D lhs_solver =
      workspace_scope.view<Real3D>(Kokkos::WithoutInitializing, "xm_wpxp_solve::lhs_solver",
                                   ngrdcol, nsup + nsub + 1, 2 * nzm - 1);
  Real3D rhs_input;
  if (reverse) {
    rhs_input = workspace_scope.view<Real3D>(
        Kokkos::WithoutInitializing, "xm_wpxp_solve::rhs_input", ngrdcol, 2 * nzm - 1, nrhs);
    Kokkos::deep_copy(rhs_input, rhs);
  }
  Kokkos::parallel_for(
      "xm_wpxp_solve::matrix_order", ngrdcol * (2 * nzm - 1), KOKKOS_LAMBDA(int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        const int source_k = reverse ? 2 * nzm - 2 - k : k;
        for (int b = 0; b < nsup + nsub + 1; ++b)
          lhs_solver(i, b, k) = lhs(reverse ? nsup + nsub - b : b, i, source_k);
        if (reverse)
          for (int n = 0; n < nrhs; ++n)
            rhs(i, k, n) = rhs_input(i, source_k, n);
      });
  if (reverse && penta_solve_method == penta_bicgstab) {
    // Keep the native optional initial-guess permutation. The current lower
    // solver rejects BiCGStab with the same supported scope as JAX.
    Real3D old_input = workspace_scope.view<Real3D>(
        Kokkos::WithoutInitializing, "xm_wpxp_solve::old_input", ngrdcol, 2 * nzm - 1, nrhs);
    Kokkos::deep_copy(old_input, old_solution);
    Kokkos::parallel_for(
        "xm_wpxp_solve::old_order", ngrdcol * (2 * nzm - 1), KOKKOS_LAMBDA(int cell) {
          int i = cell % ngrdcol, k = cell / ngrdcol;
          for (int n = 0; n < nrhs; ++n)
            old_solution(i, k, n) = old_input(i, 2 * nzm - 2 - k, n);
        });
  }
  Real3D rhs_solver      = to_execution_space(workspace, rhs),
         solution_solver = workspace_scope.view<Real3D>(
             Kokkos::WithoutInitializing, "xm_wpxp_solve::solution", ngrdcol, 2 * nzm - 1, nrhs);
  Real3D old_solver = to_execution_space(workspace, old_solution);
  Real1D condition;
  if (rcond)
    condition =
        workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "xm_wpxp_solve::rcond", ngrdcol);
  band_solve(workspace, "xm_wpxp", penta_solve_method, ngrdcol, nsup, nsub, 2 * nzm - 1, nrhs,
             l_implemented, lhs_solver, rhs_solver, err_info, solution_solver, &old_solver,
             rcond ? &condition : nullptr);
  from_execution_space(workspace, rhs, rhs_solver);
  if (rcond)
    from_execution_space(workspace, *rcond, condition);
  Kokkos::parallel_for(
      "xm_wpxp_solve::restore", ngrdcol * (2 * nzm - 1), KOKKOS_LAMBDA(int cell) {
        int i = cell % ngrdcol, k = cell / ngrdcol;
        for (int b = 0; b < nsup + nsub + 1; ++b)
          lhs(b, i, k) = lhs_solver(i, b, k);
        for (int n = 0; n < nrhs; ++n)
          solution(i, k, n) = solution_solver(i, reverse ? 2 * nzm - 2 - k : k, n);
      });
  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      std::cerr << err_info.err_header_global << "\nError in xm_wpxp_solve\n";
      return;
    }
  }
}

// solve_type: Variables being solved for.
// dt: Timestep   [s]
// wp2: w'^2 (momentum levels)                   [m^2/s^2]
// xp2: x'^2 (momentum levels)                   [{xm units}^2]
// wm_zt: w wind component on thermodynamic levels [m/s]
// xm_forcing: xm forcings (thermodynamic levels)       [units vary]
// rho_ds_zm: Dry, static density on momentum levels   [kg/m^3]
// rho_ds_zt: Dry, static density on thermo. levels    [kg/m^3]
// invrs_rho_ds_zm: Inv. dry, static density @ moment. levs. [m^3/kg]
// invrs_rho_ds_zt: Inv. dry, static density @ thermo. levs. [m^3/kg]
// xp2_threshold: Minimum allowable value of x'^2   [units vary]
// xm_threshold: Minimum allowable value of xm     [units vary]
// rcond: Reciprocal of the estimated condition number (from computing A^-1)
// low_lev_effect: Index of the lowest level that has an effect.
// high_lev_effect: Index of the highest level that has an effect.
// lhs_ma_zt: Mean advection contributions to lhs
// lhs_ma_zm: Mean advection contributions to lhs
// lhs_ta_wpxp: Turbulent advection contributions to lhs
// lhs_diff_zm: Diffusion term for w'x'
// lhs_tp: Turbulent production terms of w'x'
// lhs_ta_xm: Turbulent advection terms of xm
// lhs_pr1: Pressure term 1 for w'x'
// l_implemented: True if CLUBB is being implemented and run in a host model
// solution: The <t+1> value of xm and wpxp   [units vary]
// tridiag_solve_method: Specifier for method to solve tridiagonal systems
// fill_holes_type: Specifier for which hole filling method to use
// l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent or mean advection terms. It
// affects rtm, thlm, sclrm, um and vm.
// l_enable_relaxed_clipping: Flag to relax clipping on wpxp in xm_wpxp_clipping_and_stats
// l_mono_flux_lim_thlm: Flag to turn on monotonic flux limiter for thlm
// l_mono_flux_lim_rtm: Flag to turn on monotonic flux limiter for rtm
// l_mono_flux_lim_um: Flag to turn on monotonic flux limiter for um
// l_mono_flux_lim_vm: Flag to turn on monotonic flux limiter for vm
// l_mono_flux_lim_spikefix: Flag to implement monotonic flux limiter code that eliminates
// spurious drying tendencies at model top
// xm: The mean x field  [units vary]
// xm_tol: Minimum allowable deviation of xm [units vary]
// wpxp: The flux of x     [units vary m/s]
// err_info: err_info struct containing err_code and err_header
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21>
void xm_wpxp_clipping_and_stats(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol,
                                const grid& gr, int solve_type,
                                real_t dt,                                                 // In
                                Array0 wp2, Array1 xp2, Array2 wm_zt,                      // In
                                Array3 xm_forcing, Array4 rho_ds_zm, Array5 rho_ds_zt,     // In
                                Array6 invrs_rho_ds_zm, Array7 invrs_rho_ds_zt,            // In
                                real_t xp2_threshold, real_t xm_threshold, Array8 rcond,   // In
                                Array9 low_lev_effect, Array10 high_lev_effect,            // In
                                Array11 lhs_ma_zt, Array12 lhs_ma_zm, Array13 lhs_ta_wpxp, // In
                                Array14 lhs_diff_zm, Array15 C7_Skw_fnc,                   // In
                                Array16 lhs_tp, Array17 lhs_ta_xm, Array18 lhs_pr1,        // In
                                bool l_implemented, Array19 solution,                      // In
                                int tridiag_solve_method,                                  // In
                                int fill_holes_type,                                       // In
                                bool l_upwind_xm_ma,                                       // In
                                bool l_enable_relaxed_clipping,                            // In
                                bool l_mono_flux_lim_thlm,                                 // In
                                bool l_mono_flux_lim_rtm,                                  // In
                                bool l_mono_flux_lim_um,                                   // In
                                bool l_mono_flux_lim_vm,                                   // In
                                bool l_mono_flux_lim_spikefix,                             // In
                                stats_type& stats,                                         // InOut
                                Array20 xm,                                                // InOut
                                real_t xm_tol,                                             // In
                                Array21 wpxp, err_info_type& err_info,                     // InOut
                                int* wpxp_cl_num // Optional InOut
) {
  // Description:
  // Clips and computes implicit stats for an artitrary xm and wpxp
  //
  // References:
  //   None
  //-----------------------------------------------------------------------
  //--------------------------- Input Variables ---------------------------
  // Variables used as part of the monotonic turbulent advection scheme.
  // Find the lowermost and uppermost grid levels that can have an effect
  // on the central thermodynamic level during the course of a time step,
  // due to the effects of turbulent advection only.
  // an upwind differencing approximation rather
  // than a centered differencing for turbulent
  // or mean advection terms. It affects rtm,
  // thlm, sclrm, um and vm.
  // i.e. TKE = 1/2 (u'^2 + v'^2 + w'^2)
  // xm_wpxp_clipping_and_stats
  // eliminates spurious drying tendencies at model top
  //--------------------------- Input/Output Variables ---------------------------
  //--------------------------- Local Variables ---------------------------
  // Indices
  auto workspace_scope = workspace.scope();
  int solve_type_cl;          // solve_type used for clipping statistics.
  int wpxp_cl_max;            // Number of clipping events in the timestep
  std::string solve_type_str; // solve_type as a string for debug output purposes
  Kokkos::View<real_t**, Kokkos::LayoutLeft> xm_old =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::xm_old", ngrdcol,
          nzt); // Old value of xm for positive definite scheme     [units vary]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> xm_pd =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::xm_pd", ngrdcol,
          nzt); // Change in xm due to the pos. def. scheme
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wpxp_pd =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::wpxp_pd", ngrdcol,
          nzm); // Change in wpxp due to the pos. def. scheme
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wpxp_chnge =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::wpxp_chnge", ngrdcol,
          nzm); // Net change in w'x' due to clipping       [units vary]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> xp2_relaxed =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::xp2_relaxed", ngrdcol,
          nzm); // Value of x'^2 * clip_factor               [units vary]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> C7_Skw_fnc_plus_one =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::C7_Skw_fnc_plus_one", ngrdcol,
          nzm); // stats only
  Kokkos::View<real_t**, Kokkos::LayoutLeft> C7_Skw_fnc_zeros =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::C7_Skw_fnc_zeros", ngrdcol,
          nzm); // stats only
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wpxp_ac =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::wpxp_ac", ngrdcol,
          nzm); // stats only
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wpxp_pr2 =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::wpxp_pr2", ngrdcol,
          nzm); // stats only
  int k;
  int i;
  std::string name_xm_ta;
  std::string name_xm_ma;
  std::string name_xm_matrix_condt_num;
  std::string name_xm_pd;
  std::string name_xm_cl;
  std::string name_xm_hf_before;
  std::string name_xm_hf_after;
  std::string name_wpxp_ma;
  std::string name_wpxp_ta;
  std::string name_wpxp_tp;
  std::string name_wpxp_ac;
  std::string name_wpxp_pr1;
  std::string name_wpxp_pr2;
  std::string name_wpxp_dp1;
  std::string name_wpxp_pd;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> tmp_zt_stats =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::tmp_zt_stats", ngrdcol, nzt);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> tmp_zm_stats =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::tmp_zm_stats", ngrdcol, nzm);
  const auto gr_grid_dir_indx = gr.grid_dir_indx;
  // --------------------------- Begin code ---------------------------
  switch (solve_type) {
  case xm_wpxp_rtm: {
    // This is a diagnostic from inverting the matrix, not a budget
    name_xm_ta               = "rtm_ta";
    name_xm_ma               = "rtm_ma";
    name_xm_pd               = "rtm_pd";
    name_xm_cl               = "rtm_cl";
    name_xm_hf_before        = "rtm_before_hf";
    name_xm_hf_after         = "rtm_after_hf";
    name_xm_matrix_condt_num = "rtm_matrix_condt_num";
    name_wpxp_ma             = "wprtp_ma";
    name_wpxp_ta             = "wprtp_ta";
    name_wpxp_tp             = "wprtp_tp";
    name_wpxp_ac             = "wprtp_ac";
    name_wpxp_pr1            = "wprtp_pr1";
    name_wpxp_pr2            = "wprtp_pr2";
    name_wpxp_dp1            = "wprtp_dp1";
    name_wpxp_pd             = "wprtp_pd";
    break;
  }
  case xm_wpxp_thlm: {
    // This is a diagnostic from inverting the matrix, not a budget
    name_xm_ta               = "thlm_ta";
    name_xm_ma               = "thlm_ma";
    name_xm_pd               = "";
    name_xm_cl               = "thlm_cl";
    name_xm_hf_before        = "thlm_before_hf";
    name_xm_hf_after         = "thlm_after_hf";
    name_xm_matrix_condt_num = "thlm_matrix_condt_num";
    name_wpxp_ma             = "wpthlp_ma";
    name_wpxp_ta             = "wpthlp_ta";
    name_wpxp_tp             = "wpthlp_tp";
    name_wpxp_ac             = "wpthlp_ac";
    name_wpxp_pr1            = "wpthlp_pr1";
    name_wpxp_pr2            = "wpthlp_pr2";
    name_wpxp_dp1            = "wpthlp_dp1";
    name_wpxp_pd             = "";
    break;
  }
  case xm_wpxp_um: {
    // This is a diagnostic from inverting the matrix, not a budget
    name_xm_ta               = "um_ta";
    name_xm_ma               = "um_ma";
    name_xm_pd               = "";
    name_xm_cl               = "";
    name_xm_hf_before        = "";
    name_xm_hf_after         = "";
    name_xm_matrix_condt_num = "";
    name_wpxp_ma             = "upwp_ma";
    name_wpxp_ta             = "upwp_ta";
    name_wpxp_tp             = "upwp_tp";
    name_wpxp_ac             = "upwp_ac";
    name_wpxp_pr1            = "upwp_pr1";
    name_wpxp_pr2            = "upwp_pr2";
    name_wpxp_dp1            = "upwp_dp1";
    name_wpxp_pd             = "";
    break;
  }
  case xm_wpxp_vm: {
    // This is a diagnostic from inverting the matrix, not a budget
    name_xm_ta               = "vm_ta";
    name_xm_ma               = "vm_ma";
    name_xm_pd               = "";
    name_xm_cl               = "";
    name_xm_hf_before        = "";
    name_xm_hf_after         = "";
    name_xm_matrix_condt_num = "";
    name_wpxp_ma             = "vpwp_ma";
    name_wpxp_ta             = "vpwp_ta";
    name_wpxp_tp             = "vpwp_tp";
    name_wpxp_ac             = "vpwp_ac";
    name_wpxp_pr1            = "vpwp_pr1";
    name_wpxp_pr2            = "vpwp_pr2";
    name_wpxp_dp1            = "vpwp_dp1";
    name_wpxp_pd             = "";
    break;
  }
  default: {
    name_xm_ta               = "";
    name_xm_ma               = "";
    name_xm_pd               = "";
    name_xm_cl               = "";
    name_xm_hf_before        = "";
    name_xm_hf_after         = "";
    name_xm_matrix_condt_num = "";
    name_wpxp_ma             = "";
    name_wpxp_ta             = "";
    name_wpxp_tp             = "";
    name_wpxp_ac             = "";
    name_wpxp_pr1            = "";
    name_wpxp_pr2            = "";
    name_wpxp_dp1            = "";
    name_wpxp_pd             = "";
    break;
  }
  }
  Kokkos::parallel_for(
      "xm_wpxp_clipping_and_stats::loop1", std::max(0, nzt) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i = (cell % ngrdcol);
        const int k = (cell / ngrdcol);
        int k_xm;
        k_xm         = 2 * k + 1;
        xm_old(i, k) = xm(i, k);
        xm(i, k)     = solution(i, k_xm);
      });
  // Copy result into output arrays
  Kokkos::parallel_for(
      "xm_wpxp_clipping_and_stats::loop2", std::max(0, nzm) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i = (cell % ngrdcol);
        const int k = (cell / ngrdcol);
        int k_wpxp;
        k_wpxp     = 2 * k;
        wpxp(i, k) = solution(i, k_wpxp);
      });
  if (stats.l_sample) {
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::array3",
        C7_Skw_fnc_zeros.extent(0) * C7_Skw_fnc_zeros.extent(1), KOKKOS_LAMBDA(int cell) {
          const int slice_0                  = (cell % C7_Skw_fnc_zeros.extent(0));
          const int slice_1                  = (cell / C7_Skw_fnc_zeros.extent(0));
          C7_Skw_fnc_zeros(slice_0, slice_1) = zero;
        });
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::array4",
        C7_Skw_fnc_plus_one.extent(0) * C7_Skw_fnc_plus_one.extent(1), KOKKOS_LAMBDA(int cell) {
          const int slice_0                     = (cell % C7_Skw_fnc_plus_one.extent(0));
          const int slice_1                     = (cell / C7_Skw_fnc_plus_one.extent(0));
          C7_Skw_fnc_plus_one(slice_0, slice_1) = C7_Skw_fnc(slice_0, slice_1) + one;
        });
    // Note:  To find the contribution of w'x' term ac,
    //        substitute 0 for the C_7 skewness function input
    //        to function wpxp_terms_ac_pr2_lhs.
    wpxp_terms_ac_pr2_lhs(nzm, nzt, ngrdcol, gr, C7_Skw_fnc_zeros, // In
                          wm_zt, gr.invrs_dzm,                     // In
                          wpxp_ac                                  // Out
    );
    // Note:  To find the contribution of w'x' term pr2,
    //        add 1 to the C_7 skewness function input
    //        to function wpxp_terms_ac_pr2_lhs.
    wpxp_terms_ac_pr2_lhs(nzm, nzt, ngrdcol, gr, C7_Skw_fnc_plus_one, // In
                          wm_zt, gr.invrs_dzm,                        // In
                          wpxp_pr2                                    // Out
    );
  }
  if (stats.l_sample) {
    // Est. of the condition number of the mean/flux LHS matrix.
    if (var_on_stats_list(stats, name_xm_matrix_condt_num)) {
      {
        auto stats_value_scope = workspace.scope();
        Real1D stats_values_5  = stats_value_scope.view<Real1D>(
            Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_5",
            rcond.extent(0));
        Kokkos::parallel_for(
            "xm_wpxp_clipping_and_stats::stats_value", rcond.size(), KOKKOS_LAMBDA(int cell) {
              const int i       = cell;
              stats_values_5(i) = one / rcond(i);
            });
        stats_update(name_xm_matrix_condt_num, stats_values_5, // In
                     stats                                     // InOut
        );
      }
    }
    // Finalize implicit contributions for xm.
    if (!l_implemented) {
      Kokkos::parallel_for(
          "xm_wpxp_clipping_and_stats::array6", tmp_zt_stats.extent(0) * tmp_zt_stats.extent(1),
          KOKKOS_LAMBDA(int cell) {
            const int slice_0              = (cell % tmp_zt_stats.extent(0));
            const int slice_1              = (cell / tmp_zt_stats.extent(0));
            tmp_zt_stats(slice_0, slice_1) = zero;
          });
      Kokkos::parallel_for(
          "xm_wpxp_clipping_and_stats::loop7", std::max(0, nzt), KOKKOS_LAMBDA(const int cell) {
            const int k = cell;
            int i;
            int km1;
            int kp1;
            km1 = (-gr_grid_dir_indx + k > 0 ? -gr_grid_dir_indx + k : 0);
            kp1 = (gr_grid_dir_indx + k < nzt - 1 ? gr_grid_dir_indx + k : nzt - 1);
            for (i = 0; i < ngrdcol; i += 1) {
              tmp_zt_stats(i, k) = (-lhs_ma_zt(gr_grid_dir_indx + 1, i, k)) * xm(i, km1) +
                                   (-lhs_ma_zt(1, i, k)) * xm(i, k) +
                                   (-lhs_ma_zt(-gr_grid_dir_indx + 1, i, k)) * xm(i, kp1);
            }
          });
      // xm term ma is completely implicit; call stat_update_var_pt.
      stats_update(name_xm_ma, tmp_zt_stats, // In
                   stats                     // InOut
      );
    }
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::array8", tmp_zt_stats.extent(0) * tmp_zt_stats.extent(1),
        KOKKOS_LAMBDA(int cell) {
          const int slice_0              = (cell % tmp_zt_stats.extent(0));
          const int slice_1              = (cell / tmp_zt_stats.extent(0));
          tmp_zt_stats(slice_0, slice_1) = zero;
        });
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::loop9", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          tmp_zt_stats(i, k) =
              (-lhs_ta_xm(1, i, k)) * wpxp(i, k) + (-lhs_ta_xm(0, i, k)) * wpxp(i, k + 1);
        });
    // xm term ta is completely implicit; call stat_update_var_pt.
    stats_update(name_xm_ta, tmp_zt_stats, // In
                 stats                     // InOut
    );
    // The wpxp loop runs between k = 2 and k = nz-1.  The value of wpxp
    // is set to specified values at both the lowest level, k = 1, and the
    // highest level, k = nz.  Thus, the statistical code will run from
    // levels 2 through nz-1.
    // Finalize implicit contributions for wpxp
    // Finalize implicit contributions for wpxp.
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::array10", tmp_zm_stats.extent(0) * tmp_zm_stats.extent(1),
        KOKKOS_LAMBDA(int cell) {
          const int slice_0              = (cell % tmp_zm_stats.extent(0));
          const int slice_1              = (cell / tmp_zm_stats.extent(0));
          tmp_zm_stats(slice_0, slice_1) = zero;
        });
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::loop11", std::max(0, nzm - 2), KOKKOS_LAMBDA(const int cell) {
          const int k = cell + 1;
          int i;
          int km1;
          int kp1;
          km1 = (-gr_grid_dir_indx + k > 0 ? -gr_grid_dir_indx + k : 0);
          kp1 = (gr_grid_dir_indx + k < nzm - 1 ? gr_grid_dir_indx + k : nzm - 1);
          for (i = 0; i < ngrdcol; i += 1) {
            tmp_zm_stats(i, k) = (-lhs_ma_zm(gr_grid_dir_indx + 1, i, k)) * wpxp(i, km1) +
                                 (-lhs_ma_zm(1, i, k)) * wpxp(i, k) +
                                 (-lhs_ma_zm(-gr_grid_dir_indx + 1, i, k)) * wpxp(i, kp1);
          }
        });
    // w'x' term ma is completely implicit; call stat_update_var_pt.
    stats_update(name_wpxp_ma, tmp_zm_stats, // In
                 stats                       // InOut
    );
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::array12", tmp_zm_stats.extent(0) * tmp_zm_stats.extent(1),
        KOKKOS_LAMBDA(int cell) {
          const int slice_0              = (cell % tmp_zm_stats.extent(0));
          const int slice_1              = (cell / tmp_zm_stats.extent(0));
          tmp_zm_stats(slice_0, slice_1) = zero;
        });
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::loop13", std::max(0, nzm - 2), KOKKOS_LAMBDA(const int cell) {
          const int k = cell + 1;
          int i;
          int km1;
          int kp1;
          km1 = (-gr_grid_dir_indx + k > 0 ? -gr_grid_dir_indx + k : 0);
          kp1 = (gr_grid_dir_indx + k < nzm - 1 ? gr_grid_dir_indx + k : nzm - 1);
          for (i = 0; i < ngrdcol; i += 1) {
            tmp_zm_stats(i, k) =
                (-gamma_over_implicit_ts * lhs_ta_wpxp(gr_grid_dir_indx + 1, i, k)) * wpxp(i, km1) +
                (-gamma_over_implicit_ts * lhs_ta_wpxp(1, i, k)) * wpxp(i, k) +
                (-gamma_over_implicit_ts * lhs_ta_wpxp(-gr_grid_dir_indx + 1, i, k)) * wpxp(i, kp1);
          }
        });
    // w'x' term ta has implicit and explicit components.
    stats_finalize_budget(name_wpxp_ta, tmp_zm_stats, // In
                          stats                       // InOut
    );
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::array14", tmp_zm_stats.extent(0) * tmp_zm_stats.extent(1),
        KOKKOS_LAMBDA(int cell) {
          const int slice_0              = (cell % tmp_zm_stats.extent(0));
          const int slice_1              = (cell / tmp_zm_stats.extent(0));
          tmp_zm_stats(slice_0, slice_1) = zero;
        });
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::loop15", std::max(0, nzm - 2) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i        = (cell % ngrdcol);
          const int k        = (cell / ngrdcol) + 1;
          tmp_zm_stats(i, k) = (-lhs_tp(1, i, k)) * xm(i, k - 1) + (-lhs_tp(0, i, k)) * xm(i, k);
        });
    // w'x' term tp is completely implicit; call stat_update_var_pt.
    stats_update(name_wpxp_tp, tmp_zm_stats, // In
                 stats                       // InOut
    );
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::array16", tmp_zm_stats.extent(0) * tmp_zm_stats.extent(1),
        KOKKOS_LAMBDA(int cell) {
          const int slice_0              = (cell % tmp_zm_stats.extent(0));
          const int slice_1              = (cell / tmp_zm_stats.extent(0));
          tmp_zm_stats(slice_0, slice_1) = zero;
        });
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::loop17", std::max(0, nzm - 2) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i        = (cell % ngrdcol);
          const int k        = (cell / ngrdcol) + 1;
          tmp_zm_stats(i, k) = -wpxp_ac(i, k) * wpxp(i, k);
        });
    // w'x' term ac is completely implicit; call stat_update_var_pt.
    stats_update(name_wpxp_ac, tmp_zm_stats, // In
                 stats                       // InOut
    );
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::array18", tmp_zm_stats.extent(0) * tmp_zm_stats.extent(1),
        KOKKOS_LAMBDA(int cell) {
          const int slice_0              = (cell % tmp_zm_stats.extent(0));
          const int slice_1              = (cell / tmp_zm_stats.extent(0));
          tmp_zm_stats(slice_0, slice_1) = zero;
        });
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::loop19", std::max(0, nzm - 2) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i        = (cell % ngrdcol);
          const int k        = (cell / ngrdcol) + 1;
          tmp_zm_stats(i, k) = (-gamma_over_implicit_ts * lhs_pr1(i, k)) * wpxp(i, k);
        });
    // w'x' term pr1 is normally completely implicit.  However, due to the
    // RHS contribution from the "over-implicit" weighted time step,
    // w'x' term pr1 has both implicit and explicit components.
    // w'x' term pr1 has both implicit and explicit components;
    // call stat_end_update_pt.
    // Note:  An "over-implicit" weighted time step is applied to this term.
    //        A weighting factor of greater than 1 may be used to make the
    //        term more numerically stable (see note above for LHS turbulent
    //        advection (ta) term).
    stats_finalize_budget(name_wpxp_pr1, tmp_zm_stats, // In
                          stats                        // InOut
    );
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::array20", tmp_zm_stats.extent(0) * tmp_zm_stats.extent(1),
        KOKKOS_LAMBDA(int cell) {
          const int slice_0              = (cell % tmp_zm_stats.extent(0));
          const int slice_1              = (cell / tmp_zm_stats.extent(0));
          tmp_zm_stats(slice_0, slice_1) = zero;
        });
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::loop21", std::max(0, nzm - 2) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i        = (cell % ngrdcol);
          const int k        = (cell / ngrdcol) + 1;
          tmp_zm_stats(i, k) = -wpxp_pr2(i, k) * wpxp(i, k);
        });
    // w'x' pr2 is fully implicit in this solve path.
    stats_update(name_wpxp_pr2, tmp_zm_stats, // In
                 stats                        // InOut
    );
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::array22", tmp_zm_stats.extent(0) * tmp_zm_stats.extent(1),
        KOKKOS_LAMBDA(int cell) {
          const int slice_0              = (cell % tmp_zm_stats.extent(0));
          const int slice_1              = (cell / tmp_zm_stats.extent(0));
          tmp_zm_stats(slice_0, slice_1) = zero;
        });
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::loop23", std::max(0, nzm - 2), KOKKOS_LAMBDA(const int cell) {
          const int k = cell + 1;
          int i;
          int km1;
          int kp1;
          km1 = (-gr_grid_dir_indx + k > 0 ? -gr_grid_dir_indx + k : 0);
          kp1 = (gr_grid_dir_indx + k < nzm - 1 ? gr_grid_dir_indx + k : nzm - 1);
          for (i = 0; i < ngrdcol; i += 1) {
            tmp_zm_stats(i, k) = (-lhs_diff_zm(gr_grid_dir_indx + 1, i, k)) * wpxp(i, km1) +
                                 (-lhs_diff_zm(1, i, k)) * wpxp(i, k) +
                                 (-lhs_diff_zm(-gr_grid_dir_indx + 1, i, k)) * wpxp(i, kp1);
          }
        });
    // w'x' term dp1 is completely implicit; call stat_update_var_pt.
    stats_update(name_wpxp_dp1, tmp_zm_stats, // In
                 stats                        // InOut
    );
  }
  // Apply a monotonic turbulent flux limiter to xm/w'x'.
  if ((l_mono_flux_lim_thlm && solve_type == xm_wpxp_thlm) ||
      (l_mono_flux_lim_rtm && solve_type == xm_wpxp_rtm) ||
      (l_mono_flux_lim_um && solve_type == xm_wpxp_um) ||
      (l_mono_flux_lim_vm && solve_type == xm_wpxp_vm)) {
    monotonic_turbulent_flux_limit(workspace, nzm, nzt, ngrdcol, gr, solve_type, dt, // In
                                   xm_old, xp2, wm_zt, xm_forcing,                   // In
                                   rho_ds_zm, rho_ds_zt,                             // In
                                   invrs_rho_ds_zm, invrs_rho_ds_zt,                 // In
                                   xp2_threshold, xm_tol, l_implemented,             // In
                                   low_lev_effect, high_lev_effect,                  // In
                                   tridiag_solve_method,                             // In
                                   l_upwind_xm_ma,                                   // In
                                   l_mono_flux_lim_spikefix,                         // In
                                   stats,                                            // InOut
                                   xm, wpxp, err_info                                // InOut
    );
    if (clubb_at_least_debug_level_api(0)) {
      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << err_info.err_header_global << "\n";
          std::cerr << message.str();
        }
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "calling monotonic_turbulent_flux_limit in xm_wpxp_clipping_and_stats" << "\n";
          std::cerr << message.str();
        }
        return;
      }
    }
  }
  // Apply a flux limiting positive definite scheme if the solution
  // for the mean field is negative and we're determining total water
  if (solve_type == xm_wpxp_rtm && l_pos_def) {
    // If any xm values are negative and the values at the previous
    // timestep were all non-negative, then call pos_definite_adj
    if (any_less(xm, zero) && !any_less(xm_old, zero)) {
      pos_definite_adj(workspace, nzm, nzt, ngrdcol, gr, dt, // In
                       xm, wpxp,                             // InOut
                       xm_old,                               // In
                       xm_pd, wpxp_pd                        // Out
      );
    }
    if (stats.l_sample) {
      stats_update(name_wpxp_pd, wpxp_pd, // In
                   stats                  // InOut
      );
      stats_update(name_xm_pd, xm_pd, // In
                   stats              // InOut
      );
    }
  }
  // For stats purposes
  if (stats.l_sample) {
    // Computed value before clipping.
    {
      auto stats_value_scope = workspace.scope();
      Real2D stats_values_24 = stats_value_scope.view<Real2D>(
          Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_24", xm.extent(0),
          xm.extent(1));
      Kokkos::parallel_for(
          "xm_wpxp_clipping_and_stats::stats_value", xm.size(), KOKKOS_LAMBDA(int cell) {
            const int i           = cell % xm.extent(0);
            const int k           = cell / xm.extent(0);
            stats_values_24(i, k) = xm(i, k) / dt;
          });
      stats_begin_budget(name_xm_cl, stats_values_24, // In
                         stats                        // InOut
      );
    }
  }
  if (fill_holes_type != 0 && solve_type != xm_wpxp_um && solve_type != xm_wpxp_vm) {
    if (stats.l_sample && name_xm_hf_before.size() > 0) {
      stats_update(name_xm_hf_before, xm, // In
                   stats                  // InOut
      );
    }
    if (clubb_at_least_debug_level_api(3)) {
      if (any_less(xm, xm_threshold)) {
        // Snapshot once for the Native host diagnostic loop.
        auto xm_host = to_host_space(workspace, xm);
        switch (solve_type) {
        case xm_wpxp_rtm: {
          solve_type_str = "rtm";
          break;
        }
        case xm_wpxp_thlm: {
          solve_type_str = "thlm";
          break;
        }
        default: {
          solve_type_str = "scalars";
          break;
        }
        }
        for (i = 0; i < ngrdcol; i += 1) {
          for (k = 0; k < nzt; k += 1) {
            if (xm_host(i, k) < xm_threshold) {
              {
                std::ostringstream message;
                message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
                message << solve_type_str + " < " << " " << xm_threshold << " "
                        << " in advance_xm_wpxp_module at k= " << " " << (k + 1) << " "
                        << "i=" << " " << (i + 1) << "\n";
                std::cerr << message.str();
              }
            }
          }
        }
      }
    }
    // Include all thermodynamic levels in hole filling.
    fill_holes_vertical_api(workspace, nzt, ngrdcol, xm_threshold, // In
                            gr.k_lb_zt, gr.k_ub_zt,                // In
                            gr.dzt, rho_ds_zt, gr.grid_dir_indx,   // In
                            fill_holes_type,                       // In
                            xm                                     // InOut
    );
    if (stats.l_sample && name_xm_hf_after.size() > 0) {
      stats_update(name_xm_hf_after, xm, // In
                   stats                 // InOut
      );
    }
  }
  if (stats.l_sample) {
    {
      auto stats_value_scope = workspace.scope();
      Real2D stats_values_25 = stats_value_scope.view<Real2D>(
          Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_25", xm.extent(0),
          xm.extent(1));
      Kokkos::parallel_for(
          "xm_wpxp_clipping_and_stats::stats_value", xm.size(), KOKKOS_LAMBDA(int cell) {
            const int i           = cell % xm.extent(0);
            const int k           = cell / xm.extent(0);
            stats_values_25(i, k) = xm(i, k) / dt;
          });
      stats_finalize_budget(name_xm_cl, stats_values_25, // In
                            stats                        // InOut
      );
    }
  }
  // Clipping for w'x'
  // Clipping w'x' at each vertical level, based on the
  // correlation of w and x at each vertical level, such that:
  // corr_(w,x) = w'x' / [ sqrt(w'^2) * sqrt(x'^2) ];
  // -1 <= corr_(w,x) <= 1.
  // Since w'^2, x'^2, and w'x' are updated in different places
  // from each other, clipping for w'x' has to be done three times
  // (three times each for w'r_t', w'th_l', and w'sclr').  This is
  // the second instance of w'x' clipping.
  // Compute a slightly larger value of rt'^2 for clipping purposes.  This was
  // added to prevent a situation in which both the variance and flux are small
  // and the simulation gets "stuck" at the rt_tol^2 value.
  // See ticket #389 on the CLUBB TRAC for further details.
  // -dschanen 10 Jan 2011
  if (l_enable_relaxed_clipping) {
    if (solve_type == xm_wpxp_rtm) {
      Kokkos::parallel_for(
          "xm_wpxp_clipping_and_stats::loop26", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i       = (cell % ngrdcol);
            const int k       = (cell / ngrdcol);
            xp2_relaxed(i, k) = ((real_t(1e-7)) > (xp2(i, k)) ? (real_t(1e-7)) : (xp2(i, k)));
          });
    } else if (solve_type == xm_wpxp_thlm) {
      Kokkos::parallel_for(
          "xm_wpxp_clipping_and_stats::loop27", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i       = (cell % ngrdcol);
            const int k       = (cell / ngrdcol);
            xp2_relaxed(i, k) = ((real_t(0.01)) > (xp2(i, k)) ? (real_t(0.01)) : (xp2(i, k)));
          });
    } else {
      Kokkos::parallel_for(
          "xm_wpxp_clipping_and_stats::loop28", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i       = (cell % ngrdcol);
            const int k       = (cell / ngrdcol);
            xp2_relaxed(i, k) = ((real_t(1e-7)) > (xp2(i, k)) ? (real_t(1e-7)) : (xp2(i, k)));
          });
    }
  } else {
    Kokkos::parallel_for(
        "xm_wpxp_clipping_and_stats::loop29", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i       = (cell % ngrdcol);
          const int k       = (cell / ngrdcol);
          xp2_relaxed(i, k) = xp2(i, k);
        });
  }
  // Use solve_type to find solve_type_cl, which is used
  // in subroutine clip_covar.
  switch (solve_type) {
  case xm_wpxp_rtm: {
    solve_type_cl = clip_wprtp;
    wpxp_cl_max   = wprtp_cl_max;
    break;
  }
  case xm_wpxp_thlm: {
    solve_type_cl = clip_wpthlp;
    wpxp_cl_max   = wpthlp_cl_max;
    break;
  }
  case xm_wpxp_um: {
    solve_type_cl = clip_upwp;
    wpxp_cl_max   = upwp_cl_max;
    break;
  }
  case xm_wpxp_vm: {
    solve_type_cl = clip_vpwp;
    wpxp_cl_max   = vpwp_cl_max;
    break;
  }
  default: {
    solve_type_cl = clip_wpsclrp;
    wpxp_cl_max   = 0;
    break;
  }
  }
  if ((wpxp_cl_num != nullptr)) {
    if (stats.l_sample) {
      if (solve_type_cl == clip_wprtp) {
        if ((*wpxp_cl_num) == 0) {
          {
            auto stats_value_scope = workspace.scope();
            Real2D stats_values_30 = stats_value_scope.view<Real2D>(
                Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_30",
                wpxp.extent(0), wpxp.extent(1));
            Kokkos::parallel_for(
                "xm_wpxp_clipping_and_stats::stats_value", wpxp.size(), KOKKOS_LAMBDA(int cell) {
                  const int i           = cell % wpxp.extent(0);
                  const int k           = cell / wpxp.extent(0);
                  stats_values_30(i, k) = wpxp(i, k) / dt;
                });
            stats_begin_budget("wprtp_cl", stats_values_30, // In
                               stats                        // InOut
            );
          }
        } else {
          {
            auto stats_value_scope = workspace.scope();
            Real2D stats_values_31 = stats_value_scope.view<Real2D>(
                Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_31",
                wpxp.extent(0), wpxp.extent(1));
            Kokkos::parallel_for(
                "xm_wpxp_clipping_and_stats::stats_value", wpxp.size(), KOKKOS_LAMBDA(int cell) {
                  const int i           = cell % wpxp.extent(0);
                  const int k           = cell / wpxp.extent(0);
                  stats_values_31(i, k) = -wpxp(i, k) / dt;
                });
            stats_update_budget("wprtp_cl", stats_values_31, // In
                                stats                        // InOut
            );
          }
        }
      } else if (solve_type_cl == clip_wpthlp) {
        if ((*wpxp_cl_num) == 0) {
          {
            auto stats_value_scope = workspace.scope();
            Real2D stats_values_32 = stats_value_scope.view<Real2D>(
                Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_32",
                wpxp.extent(0), wpxp.extent(1));
            Kokkos::parallel_for(
                "xm_wpxp_clipping_and_stats::stats_value", wpxp.size(), KOKKOS_LAMBDA(int cell) {
                  const int i           = cell % wpxp.extent(0);
                  const int k           = cell / wpxp.extent(0);
                  stats_values_32(i, k) = wpxp(i, k) / dt;
                });
            stats_begin_budget("wpthlp_cl", stats_values_32, // In
                               stats                         // InOut
            );
          }
        } else {
          {
            auto stats_value_scope = workspace.scope();
            Real2D stats_values_33 = stats_value_scope.view<Real2D>(
                Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_33",
                wpxp.extent(0), wpxp.extent(1));
            Kokkos::parallel_for(
                "xm_wpxp_clipping_and_stats::stats_value", wpxp.size(), KOKKOS_LAMBDA(int cell) {
                  const int i           = cell % wpxp.extent(0);
                  const int k           = cell / wpxp.extent(0);
                  stats_values_33(i, k) = -wpxp(i, k) / dt;
                });
            stats_update_budget("wpthlp_cl", stats_values_33, // In
                                stats                         // InOut
            );
          }
        }
      } else if (solve_type_cl == clip_upwp) {
        if ((*wpxp_cl_num) == 0) {
          {
            auto stats_value_scope = workspace.scope();
            Real2D stats_values_34 = stats_value_scope.view<Real2D>(
                Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_34",
                wpxp.extent(0), wpxp.extent(1));
            Kokkos::parallel_for(
                "xm_wpxp_clipping_and_stats::stats_value", wpxp.size(), KOKKOS_LAMBDA(int cell) {
                  const int i           = cell % wpxp.extent(0);
                  const int k           = cell / wpxp.extent(0);
                  stats_values_34(i, k) = wpxp(i, k) / dt;
                });
            stats_begin_budget("upwp_cl", stats_values_34, // In
                               stats                       // InOut
            );
          }
        } else {
          {
            auto stats_value_scope = workspace.scope();
            Real2D stats_values_35 = stats_value_scope.view<Real2D>(
                Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_35",
                wpxp.extent(0), wpxp.extent(1));
            Kokkos::parallel_for(
                "xm_wpxp_clipping_and_stats::stats_value", wpxp.size(), KOKKOS_LAMBDA(int cell) {
                  const int i           = cell % wpxp.extent(0);
                  const int k           = cell / wpxp.extent(0);
                  stats_values_35(i, k) = -wpxp(i, k) / dt;
                });
            stats_update_budget("upwp_cl", stats_values_35, // In
                                stats                       // InOut
            );
          }
        }
      } else if (solve_type_cl == clip_vpwp) {
        if ((*wpxp_cl_num) == 0) {
          {
            auto stats_value_scope = workspace.scope();
            Real2D stats_values_36 = stats_value_scope.view<Real2D>(
                Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_36",
                wpxp.extent(0), wpxp.extent(1));
            Kokkos::parallel_for(
                "xm_wpxp_clipping_and_stats::stats_value", wpxp.size(), KOKKOS_LAMBDA(int cell) {
                  const int i           = cell % wpxp.extent(0);
                  const int k           = cell / wpxp.extent(0);
                  stats_values_36(i, k) = wpxp(i, k) / dt;
                });
            stats_begin_budget("vpwp_cl", stats_values_36, // In
                               stats                       // InOut
            );
          }
        } else {
          {
            auto stats_value_scope = workspace.scope();
            Real2D stats_values_37 = stats_value_scope.view<Real2D>(
                Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_37",
                wpxp.extent(0), wpxp.extent(1));
            Kokkos::parallel_for(
                "xm_wpxp_clipping_and_stats::stats_value", wpxp.size(), KOKKOS_LAMBDA(int cell) {
                  const int i           = cell % wpxp.extent(0);
                  const int k           = cell / wpxp.extent(0);
                  stats_values_37(i, k) = -wpxp(i, k) / dt;
                });
            stats_update_budget("vpwp_cl", stats_values_37, // In
                                stats                       // InOut
            );
          }
        }
      }
    }
  }
  if (solve_type != xm_wpxp_um && solve_type != xm_wpxp_vm) {
    if ((wpxp_cl_num != nullptr)) {
      (*wpxp_cl_num) = (*wpxp_cl_num) + 1;
    }
    clip_covar(nzm, ngrdcol, solve_type_cl, wp2, xp2_relaxed, // In
               wpxp,                                          // InOut
               wpxp_chnge                                     // Out
    );
    if ((wpxp_cl_num != nullptr) && stats.l_sample) {
      if (solve_type_cl == clip_wprtp) {
        if ((*wpxp_cl_num) == wpxp_cl_max) {
          {
            auto stats_value_scope = workspace.scope();
            Real2D stats_values_38 = stats_value_scope.view<Real2D>(
                Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_38",
                wpxp.extent(0), wpxp.extent(1));
            Kokkos::parallel_for(
                "xm_wpxp_clipping_and_stats::stats_value", wpxp.size(), KOKKOS_LAMBDA(int cell) {
                  const int i           = cell % wpxp.extent(0);
                  const int k           = cell / wpxp.extent(0);
                  stats_values_38(i, k) = wpxp(i, k) / dt;
                });
            stats_finalize_budget("wprtp_cl", stats_values_38, // In
                                  stats                        // InOut
            );
          }
        } else {
          {
            auto stats_value_scope = workspace.scope();
            Real2D stats_values_39 = stats_value_scope.view<Real2D>(
                Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_39",
                wpxp.extent(0), wpxp.extent(1));
            Kokkos::parallel_for(
                "xm_wpxp_clipping_and_stats::stats_value", wpxp.size(), KOKKOS_LAMBDA(int cell) {
                  const int i           = cell % wpxp.extent(0);
                  const int k           = cell / wpxp.extent(0);
                  stats_values_39(i, k) = wpxp(i, k) / dt;
                });
            stats_update_budget("wprtp_cl", stats_values_39, // In
                                stats                        // InOut
            );
          }
        }
      } else if (solve_type_cl == clip_wpthlp) {
        if ((*wpxp_cl_num) == wpxp_cl_max) {
          {
            auto stats_value_scope = workspace.scope();
            Real2D stats_values_40 = stats_value_scope.view<Real2D>(
                Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_40",
                wpxp.extent(0), wpxp.extent(1));
            Kokkos::parallel_for(
                "xm_wpxp_clipping_and_stats::stats_value", wpxp.size(), KOKKOS_LAMBDA(int cell) {
                  const int i           = cell % wpxp.extent(0);
                  const int k           = cell / wpxp.extent(0);
                  stats_values_40(i, k) = wpxp(i, k) / dt;
                });
            stats_finalize_budget("wpthlp_cl", stats_values_40, // In
                                  stats                         // InOut
            );
          }
        } else {
          {
            auto stats_value_scope = workspace.scope();
            Real2D stats_values_41 = stats_value_scope.view<Real2D>(
                Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_41",
                wpxp.extent(0), wpxp.extent(1));
            Kokkos::parallel_for(
                "xm_wpxp_clipping_and_stats::stats_value", wpxp.size(), KOKKOS_LAMBDA(int cell) {
                  const int i           = cell % wpxp.extent(0);
                  const int k           = cell / wpxp.extent(0);
                  stats_values_41(i, k) = wpxp(i, k) / dt;
                });
            stats_update_budget("wpthlp_cl", stats_values_41, // In
                                stats                         // InOut
            );
          }
        }
      }
    }
  } else {
    if ((wpxp_cl_num != nullptr)) {
      (*wpxp_cl_num) = (*wpxp_cl_num) + 1;
    }
    clip_covar(nzm, ngrdcol, solve_type_cl, wp2, xp2, // In
               wpxp,                                  // InOut
               wpxp_chnge                             // Out
    );
    if ((wpxp_cl_num != nullptr) && stats.l_sample) {
      if (solve_type_cl == clip_upwp) {
        if ((*wpxp_cl_num) == wpxp_cl_max) {
          {
            auto stats_value_scope = workspace.scope();
            Real2D stats_values_42 = stats_value_scope.view<Real2D>(
                Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_42",
                wpxp.extent(0), wpxp.extent(1));
            Kokkos::parallel_for(
                "xm_wpxp_clipping_and_stats::stats_value", wpxp.size(), KOKKOS_LAMBDA(int cell) {
                  const int i           = cell % wpxp.extent(0);
                  const int k           = cell / wpxp.extent(0);
                  stats_values_42(i, k) = wpxp(i, k) / dt;
                });
            stats_finalize_budget("upwp_cl", stats_values_42, // In
                                  stats                       // InOut
            );
          }
        } else {
          {
            auto stats_value_scope = workspace.scope();
            Real2D stats_values_43 = stats_value_scope.view<Real2D>(
                Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_43",
                wpxp.extent(0), wpxp.extent(1));
            Kokkos::parallel_for(
                "xm_wpxp_clipping_and_stats::stats_value", wpxp.size(), KOKKOS_LAMBDA(int cell) {
                  const int i           = cell % wpxp.extent(0);
                  const int k           = cell / wpxp.extent(0);
                  stats_values_43(i, k) = wpxp(i, k) / dt;
                });
            stats_update_budget("upwp_cl", stats_values_43, // In
                                stats                       // InOut
            );
          }
        }
      } else if (solve_type_cl == clip_vpwp) {
        if ((*wpxp_cl_num) == wpxp_cl_max) {
          {
            auto stats_value_scope = workspace.scope();
            Real2D stats_values_44 = stats_value_scope.view<Real2D>(
                Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_44",
                wpxp.extent(0), wpxp.extent(1));
            Kokkos::parallel_for(
                "xm_wpxp_clipping_and_stats::stats_value", wpxp.size(), KOKKOS_LAMBDA(int cell) {
                  const int i           = cell % wpxp.extent(0);
                  const int k           = cell / wpxp.extent(0);
                  stats_values_44(i, k) = wpxp(i, k) / dt;
                });
            stats_finalize_budget("vpwp_cl", stats_values_44, // In
                                  stats                       // InOut
            );
          }
        } else {
          {
            auto stats_value_scope = workspace.scope();
            Real2D stats_values_45 = stats_value_scope.view<Real2D>(
                Kokkos::WithoutInitializing, "xm_wpxp_clipping_and_stats::stats_values_45",
                wpxp.extent(0), wpxp.extent(1));
            Kokkos::parallel_for(
                "xm_wpxp_clipping_and_stats::stats_value", wpxp.size(), KOKKOS_LAMBDA(int cell) {
                  const int i           = cell % wpxp.extent(0);
                  const int k           = cell / wpxp.extent(0);
                  stats_values_45(i, k) = wpxp(i, k) / dt;
                });
            stats_update_budget("vpwp_cl", stats_values_45, // In
                                stats                       // InOut
            );
          }
        }
      }
    }
  }
  // Adjusting xm based on clipping for w'x'.
  if (l_clip_turb_adv) {
    xm_correction_wpxp_cl(workspace, nzm, nzt, ngrdcol, solve_type, dt, // In
                          wpxp_chnge, gr.invrs_dzt,                     // In
                          stats,                                        // InOut
                          xm                                            // InOut
    );
  }
  return;
}

// rho_ds_zm: Dry, static density at momentum levels     [kg/m^3]
// invrs_rho_ds_zt: Inverse dry, static density at thermo levs [m^3/kg]
// lhs_ta_xm: LHS coefficient of xm turbulent advection  [1/m]
template <class Array0, class Array1, class Array2>
void xm_term_ta_lhs(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                    Array0 rho_ds_zm, Array1 invrs_rho_ds_zt,      // In
                    Array2 lhs_ta_xm                               // Out
) {
  // Description:
  // Turbulent advection of xm:  implicit portion of the code.
  //
  // The d(xm)/dt equation contains a turbulent advection term:
  //
  // - (1/rho_ds) * d( rho_ds * w'x' )/dz.
  //
  // This term is solved for completely implicitly, such that:
  //
  // - (1/rho_ds) * d( rho_ds * w'x'(t+1) )/dz.
  //
  // Note:  When the term is brought over to the left-hand side, the sign
  //        is reversed and the leading "-" in front of the term is changed
  //        to a "+".
  //
  // The timestep index (t+1) means that the value of w'x' being used is from
  // the next timestep, which is being advanced to in solving the d(xm)/dt and
  // d(w'x')/dt equations.
  //
  // This term is discretized as follows:
  //
  // While the values of xm are found on the thermodynamic levels, the values
  // of w'x' are found on the momentum levels.  Additionally, the values of
  // rho_ds_zm are found on the momentum levels, and the values of
  // invrs_rho_ds_zt are found on the thermodynamic levels.  On the momentum
  // levels, the values of rho_ds_zm are multiplied by the values of w'x'.  The
  // derivative of (rho_ds_zm * w'x') is taken over the intermediate (central)
  // thermodynamic level, where it is multiplied by invrs_rho_ds_zt, yielding
  // the desired results.
  //
  // =====rho_ds_zm=====wpxp================================== m(k+1)
  //
  // ------invrs_rho_ds_zt--------d(rho_ds*wpxp)/dz----------- t(k)
  //
  // =====rho_ds_zm=====wpxp================================== m(k)
  //
  // The vertical indices m(k+1), t(k), and m(k) correspond with altitudes
  // zm(k+1), zt(k), and zm(k), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) )
  // References:
  //-----------------------------------------------------------------------
  // Constant parameters
  // Input Variables
  // Return Variable
  // Local Variable
  constexpr int kp1_mdiag = 0;
  constexpr int k_mdiag   = 1;
  // k: Vertical level index
  const auto gr_invrs_dzt = gr.invrs_dzt;
  // Calculate term at all grid levels.
  Kokkos::parallel_for(
      "xm_term_ta_lhs::loop1", std::max(0, nzt) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i = (cell % ngrdcol);
        const int k = (cell / ngrdcol);
        // Momentum superdiagonal [ x wpxp(k+1,<t+1>) ]
        lhs_ta_xm(kp1_mdiag, i, k) =
            +invrs_rho_ds_zt(i, k) * gr_invrs_dzt(i, k) * rho_ds_zm(i, k + 1);
        // Momentum subdiagonal [ x wpxp(k,<t+1>) ]
        lhs_ta_xm(k_mdiag, i, k) = -invrs_rho_ds_zt(i, k) * gr_invrs_dzt(i, k) * rho_ds_zm(i, k);
      });
  return;
}

// wp2: w'^2                       [m^2/s^2]
// lhs_tp: LHS coefficient of xm for turbulent production  [m/s^2]
template <class Array0, class Array1>
void wpxp_term_tp_lhs(int nzm, int ngrdcol, const grid& gr, Array0 wp2, // In
                      Array1 lhs_tp                                     // Out
) {
  // Description:
  // Turbulent production of w'x':  implicit portion of the code.
  //
  // The d(w'x')/dt equation contains a turbulent production term:
  //
  // - w'^2 d(xm)/dz.
  //
  // This term is solved for completely implicitly, such that:
  //
  // - w'^2 * d( xm(t+1) )/dz.
  //
  // Note:  When the term is brought over to the left-hand side, the sign
  //        is reversed and the leading "-" in front of the term is changed
  //        to a "+".
  //
  // The timestep index (t+1) means that the value of xm being used is from the
  // next timestep, which is being advanced to in solving the d(w'x')/dt and
  // d(xm)/dt equations.
  //
  // This term is discretized as follows:
  //
  // The values of xm are found on thermodynamic levels, while the values of
  // w'^2 are found on momentum levels.  The derivative of xm is taken over the
  // intermediate (central) momentum level, where it is multiplied by w'^2,
  // yielding the desired result.
  //
  // ---------------------------xm---------------------------- t(k)
  //
  // ==========wp2=====================d(xm)/dz=============== m(k)
  //
  // ---------------------------xm---------------------------- t(k-1)
  //
  // The vertical indices t(k), m(k), and t(k-1) correspond with altitudes
  // zt(k), zm(k), and zt(k-1), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // invrs_dzm(k) = 1 / ( zt(k) - zt(k-1) )
  // References:
  //-----------------------------------------------------------------------
  // Constant parameters
  // Input Variables
  // Return Variable
  // Local Variable
  constexpr int k_tdiag   = 0;
  constexpr int km1_tdiag = 1;
  // k: Vertical level index
  const auto gr_k_lb_zm   = gr.k_lb_zm;
  const auto gr_invrs_dzm = gr.invrs_dzm;
  const auto gr_k_ub_zm   = gr.k_ub_zm;
  // Set lower boundary to 0
  Kokkos::parallel_for(
      "wpxp_term_tp_lhs::loop1", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        const int i              = cell;
        lhs_tp(0, i, gr_k_lb_zm) = zero;
        lhs_tp(1, i, gr_k_lb_zm) = zero;
      });
  // Calculate term at all interior grid levels.
  Kokkos::parallel_for(
      "wpxp_term_tp_lhs::loop2", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i = (cell % ngrdcol);
        const int k = (cell / ngrdcol) + 1;
        // Thermodynamic superdiagonal [ x xm(k,<t+1>) ]
        lhs_tp(k_tdiag, i, k) = +wp2(i, k) * gr_invrs_dzm(i, k);
        // Thermodynamic subdiagonal [ x xm(k-1,<t+1>) ]
        lhs_tp(km1_tdiag, i, k) = -wp2(i, k) * gr_invrs_dzm(i, k);
      });
  // Set upper boundary to 0
  Kokkos::parallel_for(
      "wpxp_term_tp_lhs::loop3", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        const int i              = cell;
        lhs_tp(0, i, gr_k_ub_zm) = real_t(0.0);
        lhs_tp(1, i, gr_k_ub_zm) = real_t(0.0);
      });
  return;
}

// C7_Skw_fnc: C_7 parameter with Sk_w applied              [-]
// wm_zt: w wind component on thermodynamic levels     [m/s]
// invrs_dzm: Inverse of grid spacing                      [1/m]
// lhs_ac_pr2: LHS coefficient of accumulation and pressure term 2 [1/s]
template <class Array0, class Array1, class Array2, class Array3>
void wpxp_terms_ac_pr2_lhs(int nzm, int nzt, int ngrdcol, const grid& gr, Array0 C7_Skw_fnc, // In
                           Array1 wm_zt, Array2 invrs_dzm,                                   // In
                           Array3 lhs_ac_pr2                                                 // Out
) {
  // Description:
  // Accumulation of w'x' and w'x' pressure term 2:  implicit portion of the
  // code.
  //
  // The d(w'x')/dt equation contains an accumulation term:
  //
  // - w'x' dw/dz;
  //
  // and pressure term 2:
  //
  // + C_7 w'x' dw/dz.
  //
  // Both the w'x' accumulation term and pressure term 2 are completely
  // implicit.  The accumulation term and pressure term 2 are combined and
  // solved together as:
  //
  // - ( 1 - C_7 ) * w'x'(t+1) * dw/dz.
  //
  // Note:  When the term is brought over to the left-hand side, the sign
  //        is reversed and the leading "-" in front of the term is changed
  //        to a "+".
  //
  // The timestep index (t+1) means that the value of w'x' being used is from
  // the next timestep, which is being advanced to in solving the d(w'x')/dt
  // equation.
  //
  // The terms are discretized as follows:
  //
  // The values of w'x' are found on momentum levels, while the values of wm_zt
  // (mean vertical velocity on thermodynamic levels) are found on
  // thermodynamic levels.  The vertical derivative of wm_zt is taken over the
  // intermediate (central) momentum level.  It is then multiplied by w'x'
  // (implicitly calculated at timestep (t+1)) and the coefficients to yield
  // the desired results.
  //
  // -------wm_zt--------------------------------------------- t(k)
  //
  // ===============d(wm_zt)/dz============wpxp=============== m(k)
  //
  // -------wm_zt--------------------------------------------- t(k-1)
  //
  // The vertical indices t(k), m(k), and t(k-1) correspond with altitudes
  // zt(k), zm(k), and zt(k-1), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // invrs_dzm(k) = 1 / ( zt(k) - zt(k-1) )
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variable
  // k: Vertical level index
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_k_ub_zm = gr.k_ub_zm;
  // Set lower boundary to 0
  Kokkos::parallel_for(
      "wpxp_terms_ac_pr2_lhs::loop1", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        const int i               = cell;
        lhs_ac_pr2(i, gr_k_lb_zm) = zero;
      });
  // Calculate term at all interior grid levels.
  Kokkos::parallel_for(
      "wpxp_terms_ac_pr2_lhs::loop2", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i = (cell % ngrdcol);
        const int k = (cell / ngrdcol) + 1;
        // Momentum main diagonal: [ x wpxp(k,<t+1>) ]
        lhs_ac_pr2(i, k) =
            (one - C7_Skw_fnc(i, k)) * invrs_dzm(i, k) * (wm_zt(i, k) - wm_zt(i, k - 1));
      });
  // Set upper boundary to 0
  Kokkos::parallel_for(
      "wpxp_terms_ac_pr2_lhs::loop3", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        const int i               = cell;
        lhs_ac_pr2(i, gr_k_ub_zm) = zero;
      });
  return;
}

// C6rt_Skw_fnc: C_6rt parameter with Sk_w applied        [-]
// C6thl_Skw_fnc: C_6thl parameter with Sk_w applied       [-]
// C7_Skw_fnc: C_7 parameter with Sk_w applied          [-]
// invrs_tau_C6_zm: Inverse time-scale tau at momentum levels   [1/s]
// l_scalar_calc: True if sclr_dim > 0
// lhs_pr1_wprtp: LHS coefficient for w'r_t' pressure term 1   [1/s]
// lhs_pr1_wpthlp: LHS coefficient for w'thl' pressure term 1   [1/s]
// lhs_pr1_wpsclrp: LHS coefficient for w'sclr' pressure term 1  [1/s]
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6>
void wpxp_term_pr1_lhs(int nzm, int ngrdcol, const grid& gr, Array0 C6rt_Skw_fnc, // In
                       Array1 C6thl_Skw_fnc, Array2 C7_Skw_fnc,                   // In
                       Array3 invrs_tau_C6_zm, bool l_scalar_calc,                // In
                       Array4 lhs_pr1_wprtp, Array5 lhs_pr1_wpthlp,               // Out
                       Array6 lhs_pr1_wpsclrp                                     // Out
) {
  // Description
  // Pressure term 1 for w'x':  implicit portion of the code.
  //
  // The d(w'x')/dt equation contains pressure term 1:
  //
  // - ( C_6 / tau_m ) w'x'.
  //
  // This term is solved for completely implicitly, such that:
  //
  // - ( C_6 / tau_m ) w'x'(t+1)
  //
  // Note:  When the term is brought over to the left-hand side, the sign
  //        is reversed and the leading "-" in front of the term is changed
  //        to a "+".
  //
  // The timestep index (t+1) means that the value of w'x' being used is from
  // the next timestep, which is being advanced to in solving the d(w'x')/dt
  // equation.
  //
  // The values of w'x' are found on the momentum levels.  The values of the
  // C_6 skewness function and time-scale tau_m are also found on the momentum
  // levels.
  //
  //-----------------------------------------------------------------------
  //--------------------------- Input Variables ---------------------------
  //--------------------------- Output Variables ---------------------------
  //--------------------------- Local Variables ---------------------------
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_k_ub_zm = gr.k_ub_zm;
  //--------------------------- Begin Code ---------------------------
  Kokkos::parallel_for(
      "wpxp_term_pr1_lhs::loop1", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i = (cell % ngrdcol);
        const int k = (cell / ngrdcol) + 1;
        // Momentum main diagonals: [ x wpxp(k,<t+1>) ]
        lhs_pr1_wprtp(i, k) = C6rt_Skw_fnc(i, k) * invrs_tau_C6_zm(i, k);
        // Momentum main diagonals: [ x wpxp(k,<t+1>) ]
        lhs_pr1_wpthlp(i, k) = C6thl_Skw_fnc(i, k) * invrs_tau_C6_zm(i, k);
      });
  Kokkos::parallel_for(
      "wpxp_term_pr1_lhs::loop2", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        const int i = cell;
        // Set lower boundary to 0
        lhs_pr1_wprtp(i, gr_k_lb_zm) = zero;
        // Set upper boundary to 0
        lhs_pr1_wprtp(i, gr_k_ub_zm) = zero;
        // Set lower boundary to 0
        lhs_pr1_wpthlp(i, gr_k_lb_zm) = zero;
        // Set upper boundary to 0
        lhs_pr1_wpthlp(i, gr_k_ub_zm) = zero;
      });
  if (l_scalar_calc) {
    Kokkos::parallel_for(
        "wpxp_term_pr1_lhs::loop3", std::max(0, nzm - 2) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol) + 1;
          // Momentum main diagonals: [ x wpxp(k,<t+1>) ]
          lhs_pr1_wpsclrp(i, k) = C7_Skw_fnc(i, k) * invrs_tau_C6_zm(i, k);
        });
    Kokkos::parallel_for(
        "wpxp_term_pr1_lhs::loop4", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i = cell;
          // Set lower boundary to 0
          lhs_pr1_wpsclrp(i, gr_k_lb_zm) = zero;
          // Set upper boundary to 0
          lhs_pr1_wpsclrp(i, gr_k_ub_zm) = zero;
        });
  }
  return;
}

// C7_Skw_fnc: C_7 parameter with Sk_w applied       [-]
// thv_ds_zm: Dry, base-state theta_v on mom. levs. [K]
// xpthvp: x'th_v'                               [K {xm units}]
// rhs_bp_pr3: RHS portion of bouyancy prod and pressure term 3
template <class Array0, class Array1, class Array2, class Array3>
void wpxp_terms_bp_pr3_rhs(int nzm, int ngrdcol, const grid& gr, Array0 C7_Skw_fnc,
                           Array1 thv_ds_zm, // In
                           Array2 xpthvp,    // In
                           Array3 rhs_bp_pr3 // Out
) {
  // Description:
  // Buoyancy production of w'x' and w'x' pressure term 3:  explicit portion of
  // the code.
  //
  // The d(w'x')/dt equation contains a buoyancy production term:
  //
  // + (g/thv_ds) x'th_v';
  //
  // and pressure term 3:
  //
  // - C_7 (g/thv_ds) x'th_v'.
  //
  // Both the w'x' buoyancy production term and pressure term 3 are completely
  // explicit.  The buoyancy production term and pressure term 3 are combined
  // and solved together as:
  //
  // + ( 1 - C_7 ) * (g/thv_ds) * x'th_v'.
  // References:
  //-----------------------------------------------------------------------
  //---------------------------- Input Variables ----------------------------
  //---------------------------- Output Variables ----------------------------
  //---------------------------- Local Variables ----------------------------
  // k: Vertical level index
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_k_ub_zm = gr.k_ub_zm;
  //---------------------------- Begin Code ----------------------------
  // Set lower boundary to 0
  Kokkos::parallel_for(
      "wpxp_terms_bp_pr3_rhs::loop1", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        const int i               = cell;
        rhs_bp_pr3(i, gr_k_lb_zm) = zero;
      });
  // Calculate term at all interior grid levels.
  Kokkos::parallel_for(
      "wpxp_terms_bp_pr3_rhs::loop2", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i      = (cell % ngrdcol);
        const int k      = (cell / ngrdcol) + 1;
        rhs_bp_pr3(i, k) = (grav / thv_ds_zm(i, k)) * (one - C7_Skw_fnc(i, k)) * xpthvp(i, k);
      });
  // Set upper boundary to 0
  Kokkos::parallel_for(
      "wpxp_terms_bp_pr3_rhs::loop3", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        const int i               = cell;
        rhs_bp_pr3(i, gr_k_ub_zm) = zero;
      });
  return;
}

// solve_type: Variable that is being solved for.
// dt: Model timestep                            [s]
// wpxp_chnge: Amount of change in w'x' due to clipping  [m/s {xm units}]
// invrs_dzt: Inverse of grid spacing                   [1/m]
// xm: xm (thermodynamic levels)                 [{xm units}]
template <class Array0, class Array1, class Array2>
void xm_correction_wpxp_cl(ScratchWorkspace& workspace,                              // InOut
                           int nzm, int nzt, int ngrdcol, int solve_type, real_t dt, // In
                           Array0 wpxp_chnge, Array1 invrs_dzt,                      // In
                           stats_type& stats,                                        // InOut
                           Array2 xm                                                 // InOut
) {
  // Description:
  // Corrects the value of xm if w'x' needed to be clipped, for xm is partially
  // based on the derivative of w'x' with respect to altitude.
  //
  // The time-tendency equation for xm is:
  //
  // d(xm)/dt = -w d(xm)/dz - d(w'x')/dz + d(xm)/dt|_ls;
  //
  // where d(xm)/dt|_ls is the rate of change of xm over time due to radiation,
  // microphysics, and/or any other large-scale forcing(s).
  //
  // The time-tendency equation for xm is solved in conjunction with the
  // time-tendency equation for w'x'.  Both equations are solved together in a
  // semi-implicit manner.  However, after both equations have been solved (and
  // thus both xm and w'x' have been advanced to the next timestep with
  // timestep index {t+1}), the value of covariance w'x' may be clipped at any
  // level in order to prevent the correlation of w and x from becoming greater
  // than 1 or less than -1.
  //
  // The correlation between w and x is:
  //
  // corr_(w,x) = w'x' / [ sqrt(w'^2) * sqrt(x'^2) ].
  //
  // The correlation must always have a value between -1 and 1, such that:
  //
  // -1 <= corr_(w,x) <= 1.
  //
  // Therefore, there is an upper limit on w'x', such that:
  //
  // w'x' <=  [ sqrt(w'^2) * sqrt(x'^2) ];
  //
  // and a lower limit on w'x', such that:
  //
  // w'x' >= -[ sqrt(w'^2) * sqrt(x'^2) ].
  //
  // The aforementioned time-tendency equation for xm is based on the value of
  // w'x' without being clipped (w'x'{t+1}_unclipped), such that:
  //
  // d(xm)/dt = -w d(xm{t+1})/dz - d(w'x'{t+1}_unclipped)/dz + d(xm{t})/dt|_ls;
  //
  // where the both the mean advection term, -w d(xm{t+1})/dz, and the
  // turbulent advection term, -d(w'x'{t+1}_unclipped)/dz, are solved
  // completely implicitly.  The xm forcing term, +d(xm{t})/dt|_ls, is solved
  // completely explicitly.
  //
  // However, if w'x' needs to be clipped after being advanced one timestep,
  // then xm needs to be altered to reflect the fact that w'x' has a different
  // value than the value used while both were being solved together.  Ideally,
  // the xm time-tendency equation that should be used is:
  //
  // d(xm)/dt = -w d(xm{t+1})/dz - d(w'x'{t+1}_clipped)/dz + d(xm{t})/dt|_ls.
  //
  // However, w'x'{t+1}_clipped isn't known until after the w'x' and xm
  // equations have been solved together.  However, a proper adjuster can be
  // applied to xm through the use of the following relationship:
  //
  // w'x'{t+1}_clipped = w'x'{t+1}_unclipped + w'x'{t+1}_amount_clipped;
  //
  // at any given vertical level.
  //
  // When the expression above is substituted into the preceeding xm
  // time-tendency equation, the resulting equation for xm time-tendency is:
  //
  // d(xm)/dt = -w d(xm{t+1})/dz - d(w'x'{t+1}_unclipped)/dz
  //               - d(w'x'{t+1}_amount_clipped)/dz + d(xm{t})/dt|_ls.
  //
  // Thus, the resulting xm time-tendency equation is the same as the original
  // xm time-tendency equation, but with added adjuster term:
  //
  // -d(w'x'{t+1}_amount_clipped)/dz.
  //
  // Since the adjuster term needs to be applied after xm has already been
  // solved, it needs to be multiplied by the timestep length and added on to
  // xm{t+1}, such that:
  //
  // xm{t+1}_after_adjustment =
  //    xm{t+1}_before_adjustment + ( -d(w'x'{t+1}_amount_clipped)/dz ) * dt.
  //
  // The adjuster term is discretized as follows:
  //
  // The values of w'x' are located on the momentum levels.  Thus, the values
  // of w'x'_amount_clipped are also located on the momentum levels.  The
  // values of xm are located on the thermodynamic levels.  The derivatives
  // (d/dz) of w'x'_amount_clipped are taken over the intermediate
  // thermodynamic levels, where they are applied to xm.
  //
  // =======wpxp_amount_clipped=============================== m(k+1)
  //
  // -----------------------------d(wpxp_amount_clipped)/dz--- t(k)
  //
  // =======wpxp_amount_clipped=============================== m(k)
  //
  // The vertical indices m(k+1), t(k), and m(k) correspond with altitudes
  // zm(k+1), zt(k), and zm(k), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) )
  // Note:  The results of this xm adjustment are highly dependent on the
  //        numerical stability and the smoothness of the w'^2 and x'^2 fields.
  //        An unstable "sawtooth" profile for w'^2 and/or x'^2 causes an
  //        unstable "sawtooth" profile for the upper and lower limits on w'x'.
  //        In turn, this causes an unstable "sawtooth" profile for
  //        w'x'_amount_clipped.  Taking the derivative of that such a "noisy"
  //        field and applying the results to xm causes the xm field to become
  //        more "noisy" and unstable.
  // References:
  //-----------------------------------------------------------------------
  //---------------------------- Input Variables ----------------------------
  //---------------------------- Input/Output Variable ----------------------------
  //---------------------------- Local Variables ----------------------------
  auto workspace_scope = workspace.scope();
  Kokkos::View<real_t**, Kokkos::LayoutLeft> xm_tndcy_wpxp_cl =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_correction_wpxp_cl::xm_tndcy_wpxp_cl", ngrdcol,
          nzt); // d(xm)/dt due to clipping of w'x'       [{xm units}/s]
  // k: Array index
  std::string name_xm_tacl;
  Kokkos::View<std::uint8_t*, Kokkos::LayoutLeft> l_clipping_needed =
      workspace_scope.view<Kokkos::View<std::uint8_t*, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "xm_correction_wpxp_cl::l_clipping_needed", ngrdcol);
  bool l_any_clipping_needed;
  //---------------------------- Begin Code ----------------------------
  l_any_clipping_needed = false;
  // Native leaves untouched per-column mask entries undefined. Define them
  // false as in JAX; no Native BFB claim applies to those undefined entries.
  // Reduce integer column hits and keep each column's vertical scan ordered.
  int clipping_columns = 0;
  Kokkos::parallel_reduce(
      "xm_correction_wpxp_cl::clipping_columns", ngrdcol,
      KOKKOS_LAMBDA(int column, int& count) {
        bool needed = false;
        for (int level = 0; level < nzm; ++level) {
          if (Kokkos::abs(wpxp_chnge(column, level)) > eps)
            needed = true;
        }
        l_clipping_needed(column) = needed;
        if (needed)
          ++count;
      },
      clipping_columns);
  l_any_clipping_needed = clipping_columns != 0;
  if (!l_any_clipping_needed) {
    return;
  }
  switch (solve_type) {
  case xm_wpxp_rtm: {
    name_xm_tacl = "rtm_tacl";
    break;
  }
  case xm_wpxp_thlm: {
    name_xm_tacl = "thlm_tacl";
    break;
  }
  default: {
    name_xm_tacl = "";
    break;
  }
  }
  // Adjusting xm based on clipping for w'x'.
  // Loop over all thermodynamic levels between the second-lowest and the
  // highest.
  // The source also leaves nonclipping-column tacl statistics unwritten.
  // Define their contribution as zero, matching JAX; physical updates below
  // retain the Native masked formula and early-return placement.
  Kokkos::deep_copy(xm_tndcy_wpxp_cl, zero);
  Kokkos::parallel_for(
      "xm_correction_wpxp_cl::loop1", std::max(0, nzt) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i = (cell % ngrdcol);
        const int k = (cell / ngrdcol);
        if (l_clipping_needed(i)) {
          xm_tndcy_wpxp_cl(i, k) = -invrs_dzt(i, k) * (wpxp_chnge(i, k + 1) - wpxp_chnge(i, k));
          xm(i, k)               = xm(i, k) + xm_tndcy_wpxp_cl(i, k) * dt;
        }
      });
  if (stats.l_sample) {
    // The adjustment to xm due to turbulent advection term clipping
    // (xm term tacl) is completely explicit; call stat_update_var.
    // (xm term tacl) is completely explicit.
    stats_update(name_xm_tacl, xm_tndcy_wpxp_cl, // In
                 stats                           // InOut
    );
  }
  return;
}

// coefficient: The coefficient to be damped
// Cx_Skw_fnc: Initial skewness function before damping
// max_coeff_value: Maximum value the damped coefficient should have
// altitude_threshold: Minimum altitude where damping should occur
// threshold: Value of Lscale below which the damping should occur
// Lscale_zm: Current value of Lscale (interp to momentum levels)
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6>
void damp_coefficient(int nzm, int ngrdcol, const grid& gr, Array0 coefficient,
                      Array1 Cx_Skw_fnc,                                 // In
                      Array2 max_coeff_value, Array3 altitude_threshold, // In
                      Array4 threshold, Array5 Lscale_zm,                // In
                      Array6 damped_value                                // Out
) {
  // Description:
  // Damps a given coefficient linearly based on the value of Lscale.
  // For additional information see CLUBB ticket #431.
  // Input variables
  // Return Variable
  // Local Variables
  const auto gr_zm = gr.zm;
  Kokkos::parallel_for(
      "damp_coefficient::loop1", std::max(0, nzm) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i = (cell % ngrdcol);
        const int k = (cell / ngrdcol);
        if (Lscale_zm(i, k) < threshold(i) && gr_zm(i, k) > altitude_threshold(i)) {
          damped_value(i, k) =
              max_coeff_value(i) +
              ((coefficient(i) - max_coeff_value(i)) / threshold(i)) * Lscale_zm(i, k);
        } else {
          damped_value(i, k) = Cx_Skw_fnc(i, k);
        }
      });
  return;
}

// ypwp: momentum flux component, either upwp or vpwp  [m^2/s^2]
// xm: grid-mean conserved thermodynamic variable, either thlm or rtm [varies]
// wpxp: vertical scalar flux, either wpthlp or wprtp [varies]
// ym: grid-mean velocity component, either um or vm [m/s]
// C6x_Skw_fnc: C_6 pressure parameter with effects of Sk_w incorporated (k)  [-]
// tau_C6_zm: Time-scale tau on momentum levels applied to C6 term [s]
// C7_Skw_fnc: C_7 pressure parameter with effects of Sk_w incorporated (k)  [-]
// ypxp: horizontal flux of a conserved scalar, either upthlp, uprtp, vpthlp, or vprtp
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7>
void diagnose_upxp(ScratchWorkspace& workspace,                                           // InOut
                   int nzm, int nzt, int ngrdcol, const grid& gr, Array0 ypwp, Array1 xm, // In
                   Array2 wpxp, Array3 ym,                                                // In
                   Array4 C6x_Skw_fnc, Array5 tau_C6_zm, Array6 C7_Skw_fnc,               // In
                   Array7 ypxp                                                            // Out
) {
  // Description:
  //   Diagnose turbulent horizontal flux of a conserved scalar.
  //
  // References:
  //   Eqn. 7 of Andre et al. (1978)
  //   Eqn. 4 of Bougeault et al. (1981)
  //   github issue #841
  //
  //-------------------------------------------------------------------------------------
  //------------------------------ Input Variables ------------------------------
  //------------------------------ Return Variables ------------------------------
  //------------------------------ Local Variables ------------------------------
  auto workspace_scope = workspace.scope();
  Kokkos::View<real_t**, Kokkos::LayoutLeft> ddzt_xm =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "diagnose_upxp::ddzt_xm", ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> ddzt_ym =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "diagnose_upxp::ddzt_ym", ngrdcol, nzm);
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_k_ub_zm = gr.k_ub_zm;
  //----------------------------- Begin Code ------------------------------
  ddzt(nzm, nzt, ngrdcol, gr, to_execution_space(workspace, xm), ddzt_xm);
  ddzt(nzm, nzt, ngrdcol, gr, to_execution_space(workspace, ym), ddzt_ym);
  Kokkos::parallel_for(
      "diagnose_upxp::loop1", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i = (cell % ngrdcol);
        const int k = (cell / ngrdcol) + 1;
        ypxp(i, k) =
            (tau_C6_zm(i, k) / C6x_Skw_fnc(i, k)) *
            (-ypwp(i, k) * ddzt_xm(i, k) - (one - C7_Skw_fnc(i, k)) * (wpxp(i, k) * ddzt_ym(i, k)));
      });
  // The value of ypxp is irrelevant to the calculations at the upper and
  // lower boundaries
  Kokkos::parallel_for(
      "diagnose_upxp::loop2", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        const int i         = cell;
        ypxp(i, gr_k_lb_zm) = real_t(0.0);
        ypxp(i, gr_k_ub_zm) = real_t(0.0);
      });
  return;
}

// zm: Momentum grid
// zt: Thermo grid
// dt: Timestep                                 [s]
// sigma_sqd_w: sigma_sqd_w on momentum levels           [-]
// wm_zm: w wind component on momentum levels      [m/s]
// wm_zt: w wind component on thermodynamic levels [m/s]
// wp2: w'^2 (momentum levels)                   [m^2/s^2]
// Lscale_zm: Turbulent mixing length interp to m levs [m]
// wp3_on_wp2: Smoothed wp3 / wp2 on momentum levels    [m/s]
// wp3_on_wp2_zt: Smoothed wp3 / wp2 on thermo. levels     [m/s]
// Kh_zt: Eddy diffusivity on thermodynamic levels [m^2/s]
// Kh_zm: Eddy diffusivity on momentum levels
// invrs_tau_C6_zm: Inverse time-scale tau on momentum levels applied to C6 term [1/s]
// Skw_zm: Skewness of w on momentum levels         [-]
// wp2rtp: <w'^2 r_t'> (thermodynamic levels)    [m^2/s^2 kg/kg]
// rtpthvp: r_t'th_v' (momentum levels)              [(kg/kg) K]
// rtm_forcing: r_t forcing (thermodynamic levels)       [(kg/kg)/s]
// wprtp_forcing: <w'r_t'> forcing (momentum levels)       [(kg/kg)/s^2]
// rtm_ref: rtm for nudging                          [kg/kg]
// wp2thlp: <w'^2 th_l'> (thermodynamic levels)      [m^2/s^2 K]
// thlpthvp: th_l'th_v' (momentum levels)             [K^2]
// thlm_forcing: th_l forcing (thermodynamic levels)      [K/s]
// wpthlp_forcing: <w'th_l'> forcing (momentum levels)      [K/s^2]
// thlm_ref: thlm for nudging                         [K]
// rho_ds_zm: Dry, static density on momentum levels   [kg/m^3]
// rho_ds_zt: Dry, static density on thermo. levels    [kg/m^3]
// invrs_rho_ds_zm: Inv. dry, static density @ moment. levs. [m^3/kg]
// invrs_rho_ds_zt: Inv. dry, static density @ thermo. levs. [m^3/kg]
// thv_ds_zm: Dry, base-state theta_v on moment. levs. [K]
// rtp2: r_t'^2 (momentum levels)                 [(kg/kg)^2]
// thlp2: th_l'^2 (momentum levels)                [K^2]
// w_1_zm: Mean w (1st PDF component)              [m/s]
// w_2_zm: Mean w (2nd PDF component)              [m/s]
// varnce_w_1_zm: Variance of w (1st PDF component)       [m^2/s^2]
// varnce_w_2_zm: Variance of w (2nd PDF component)       [m^2/s^2]
// mixt_frac_zm: Weight of 1st PDF component (Sk_w dependent) [-]
// l_implemented: True if CLUBB is being implemented and run in a host model
// em: Turbulent Kinetic Energy (TKE)           [m^2/s^2]
// wp2sclrp: <w'^2 sclr'> (thermodynamic levels)   [Units vary]
// sclrpthvp: <sclr' th_v'> (momentum levels)       [Units vary]
// sclrm_forcing: sclrm forcing (thermodynamic levels)  [Units vary]
// sclrp2: For clipping Vince Larson             [Units vary]
// Cx_fnc_Richardson: Cx_fnc computed from Richardson_num       [-]
// pdf_implicit_coefs_terms: Implicit coefs / explicit terms [units vary]
// um_forcing: <u> forcing term (thermodynamic levels)      [m/s^2]
// vm_forcing: <v> forcing term (thermodynamic levels)      [m/s^2]
// ug: <u> geostrophic wind (thermodynamic levels)  [m/s]
// vg: <v> geostrophic wind (thermodynamic levels)  [m/s]
// wpthvp: <w'thv'> (momentum levels)                   [m/s K]
// fcor: Coriolis parameter                           [s^-1]
// um_ref: Reference u wind component for nudging       [m/s]
// vm_ref: Reference v wind component for nudging       [m/s]
// up2: Variance of the u wind component             [m^2/s^2]
// vp2: Variance of the v wind component             [m^2/s^2]
// uprcp: < u' r_c' >                                  [(m kg)/(s kg)]
// vprcp: < v' r_c' >                                  [(m kg)/(s kg)]
// rc_coef_zm: Coefficient on X'r_c' in X'th_v' equation    [K/(kg/kg)]
// rtm: r_t  (total water mixing ratio)           [kg/kg]
// wprtp: w'r_t'                                    [(kg/kg) m/s]
// thlm: th_l (liquid water potential temperature) [K]
// wpthlp: w'th_l'                                   [K m/s]
// sclrm: [Units vary]
// wpsclrp: [Units vary]
// um: <u>:  mean west-east horiz. velocity (thermo. levs.)   [m/s]
// upwp: <u'w'>:  momentum flux (momentum levels)               [m^2/s^2]
// vm: <v>:  mean south-north horiz. velocity (thermo. levs.) [m/s]
// vpwp: <v'w'>:  momentum flux (momentum levels)               [m^2/s^2]
// rtm_old: Saved value of r_t        [kg/kg]
// wprtp_old: Saved value of w'r_t'     [(kg/kg) m/s]
// thlm_old: Saved value of th_l       [K]
// wpthlp_old: Saved value of w'th_l'    [K m/s]
// sclrm_old: Saved value of sclrm     [units vary]
// wpsclrp_old: Saved value of wpsclrp   [units vary]
// um_old: Saved value of <u>       [m/s]
// upwp_old: Saved value of <u'w'>    [m^2/s^2]
// vm_old: Saved value of <v>       [m/s]
// vpwp_old: Saved value of <v'w'>    [m^2/s^2]
// l_predict_upwp_vpwp: Flag to predict <u'w'> and <v'w'> along with <u> and <v> alongside the
// advancement of <rt>, <w'rt'>, <thl>, <wpthlp>, <sclr>, and <w'sclr'> in subroutine
// advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still approximated by eddy diffusivity
// when <u> and <v> are advanced in subroutine advance_windm_edsclrm.
// l_lmm_stepping: Apply Linear Multistep Method (LMM) Stepping
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26, class Array27, class Array28, class Array29,
          class Array30, class Array31, class Array32, class Array33, class Array34, class Array35,
          class Array36, class Array37, class Array38, class Array39, class Array40, class Array41,
          class Array42, class Array43, class Array44, class Array45, class Array46, class Array47,
          class Array48, class Array49, class Array50, class Array51, class Array52, class Array53,
          class Array54, class Array55, class Array56, class Array57, class Array58, class Array59,
          class Array60, class Array61, class Array62, class Array63, class Array64, class Array65,
          class Array66, class Array67, class Array68, class Array69, class Array70, class Array71,
          class Array72>
void error_prints_xm_wpxp(int nzm, int nzt, int sclr_dim, Array0 zm, Array1 zt,      // In
                          real_t dt, Array2 sigma_sqd_w, Array3 wm_zm, Array4 wm_zt, // In
                          Array5 wp2,                                                // In
                          Array6 Lscale_zm, Array7 wp3_on_wp2,                       // In
                          Array8 wp3_on_wp2_zt,                                      // In
                          Array9 Kh_zt, Array10 Kh_zm, Array11 invrs_tau_C6_zm,      // In
                          Array12 Skw_zm,                                            // In
                          Array13 wp2rtp, Array14 rtpthvp, Array15 rtm_forcing,      // In
                          Array16 wprtp_forcing, Array17 rtm_ref, Array18 wp2thlp,   // In
                          Array19 thlpthvp, Array20 thlm_forcing,                    // In
                          Array21 wpthlp_forcing, Array22 thlm_ref,                  // In
                          Array23 rho_ds_zm,                                         // In
                          Array24 rho_ds_zt, Array25 invrs_rho_ds_zm,                // In
                          Array26 invrs_rho_ds_zt, Array27 thv_ds_zm, Array28 rtp2,  // In
                          Array29 thlp2, Array30 w_1_zm, Array31 w_2_zm,             // In
                          Array32 varnce_w_1_zm, Array33 varnce_w_2_zm,              // In
                          Array34 mixt_frac_zm, bool l_implemented, Array35 em,      // In
                          Array36 wp2sclrp, Array37 sclrpthvp,                       // In
                          Array38 sclrm_forcing,                                     // In
                          Array39 sclrp2, Array40 Cx_fnc_Richardson,                 // In
                          const implicit_coefs_terms& pdf_implicit_coefs_terms,      // In
                          Array41 um_forcing, Array42 vm_forcing,                    // In
                          Array43 ug, Array44 vg,                                    // In
                          Array45 wpthvp, real_t fcor, Array46 um_ref,               // In
                          Array47 vm_ref, Array48 up2,                               // In
                          Array49 vp2, Array50 uprcp, Array51 vprcp,                 // In
                          Array52 rc_coef_zm, Array53 rtm,                           // In
                          Array54 wprtp, Array55 thlm, Array56 wpthlp,               // In
                          Array57 sclrm, Array58 wpsclrp,                            // In
                          Array59 um, Array60 upwp, Array61 vm,                      // In
                          Array62 vpwp, Array63 rtm_old,                             // In
                          Array64 wprtp_old, Array65 thlm_old,                       // In
                          Array66 wpthlp_old,                                        // In
                          Array67 sclrm_old, Array68 wpsclrp_old,                    // In
                          Array69 um_old,                                            // In
                          Array70 upwp_old, Array71 vm_old, Array72 vpwp_old,        // In
                          bool l_predict_upwp_vpwp, bool l_lmm_stepping              // In
) {
  // Description:
  // Prints values of model fields when fatal errors (LU decomp.) occur.
  // All field that are passed into and out of subroutine advance_xm_wpxp are
  // printed.  If additional fields are added to the call to subroutine
  // advance_xm_wpxp, they should also be added here.
  // Input Variables
  // Added for clipping by Vince Larson 29 Sep 2007
  // End of Vince Larson's addition.
  // Additional variables for passive scalars
  // Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
  // Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
  // Saved values of predictive fields, prior to being advanced, for use in
  // print statements in case of fatal error.
  // Input/Output Variables
  // Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
  // and <v> alongside the advancement of <rt>,
  // <w'rt'>, <thl>, <wpthlp>, <sclr>, and <w'sclr'>
  // in subroutine advance_xm_wpxp.  Otherwise, <u'w'>
  // and <v'w'> are still approximated by eddy
  // diffusivity when <u> and <v> are advanced in
  // subroutine advance_windm_edsclrm.
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "Error in advance_xm_wpxp" << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "Intent(in)" << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "zt = " << " " << dump_values(zt) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "zm = " << " " << dump_values(zm) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "dt = " << " " << dt << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "sigma_sqd_w = " << " " << dump_values(sigma_sqd_w) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "wm_zm = " << " " << dump_values(wm_zm) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "wm_zt = " << " " << dump_values(wm_zt) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "wp2 = " << " " << dump_values(wp2) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "Lscale_zm = " << " " << dump_values(Lscale_zm) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "wp3_on_wp2 = " << " " << dump_values(wp3_on_wp2) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "wp3_on_wp2_zt = " << " " << dump_values(wp3_on_wp2_zt) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "Kh_zt = " << " " << dump_values(Kh_zt) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "Kh_zm = " << " " << dump_values(Kh_zm) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "invrs_tau_C6_zm = " << " " << dump_values(invrs_tau_C6_zm) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "Skw_zm = " << " " << dump_values(Skw_zm) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "wp2rtp = " << " " << dump_values(wp2rtp) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "rtpthvp = " << " " << dump_values(rtpthvp) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "rtm_forcing = " << " " << dump_values(rtm_forcing) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "wprtp_forcing = " << " " << dump_values(wprtp_forcing) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "rtm_ref = " << " " << dump_values(rtm_ref) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "wp2thlp = " << " " << dump_values(wp2thlp) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "thlpthvp = " << " " << dump_values(thlpthvp) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "thlm_forcing = " << " " << dump_values(thlm_forcing) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "wpthlp_forcing = " << " " << dump_values(wpthlp_forcing) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "thlm_ref = " << " " << dump_values(thlm_ref) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "rho_ds_zm = " << " " << dump_values(rho_ds_zm) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "rho_ds_zt = " << " " << dump_values(rho_ds_zt) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "invrs_rho_ds_zm = " << " " << dump_values(invrs_rho_ds_zm) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "invrs_rho_ds_zt = " << " " << dump_values(invrs_rho_ds_zt) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "thv_ds_zm = " << " " << dump_values(thv_ds_zm) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "rtp2 = " << " " << dump_values(rtp2) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "thlp2 = " << " " << dump_values(thlp2) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "w_1_zm = " << " " << dump_values(w_1_zm) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "w_2_zm = " << " " << dump_values(w_2_zm) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "varnce_w_1_zm = " << " " << dump_values(varnce_w_1_zm) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "varnce_w_2_zm = " << " " << dump_values(varnce_w_2_zm) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "mixt_frac_zm = " << " " << dump_values(mixt_frac_zm) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "l_implemented = " << " " << l_implemented << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "em = " << " " << dump_values(em) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "Cx_fnc_Richardson = " << " " << dump_values(Cx_fnc_Richardson) << " " << "\n"
            << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "pdf_implicit_coefs_terms%coef_wp2rtp_implicit = " << " "
            << dump_values(pdf_implicit_coefs_terms.coef_wp2rtp_implicit) << " " << "\n"
            << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "pdf_implicit_coefs_terms%term_wp2rtp_explicit = " << " "
            << dump_values(pdf_implicit_coefs_terms.term_wp2rtp_explicit) << " " << "\n"
            << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "pdf_implicit_coefs_terms%coef_wp2thlp_implicit = " << " "
            << dump_values(pdf_implicit_coefs_terms.coef_wp2thlp_implicit) << " " << "\n"
            << "\n";
    std::cerr << message.str();
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "pdf_implicit_coefs_terms%term_wp2thlp_explicit = " << " "
            << dump_values(pdf_implicit_coefs_terms.term_wp2thlp_explicit) << " " << "\n"
            << "\n";
    std::cerr << message.str();
  }
  if (sclr_dim > 0) {
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "sclrp2 = " << " " << dump_values(sclrp2) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "wp2sclrp = " << " " << dump_values(wp2sclrp) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "sclrpthvp = " << " " << dump_values(sclrpthvp) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "sclrm_forcing = " << " " << dump_values(sclrm_forcing) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "pdf_implicit_coefs_terms%coef_wp2sclrp_implicit = " << " "
              << dump_values(pdf_implicit_coefs_terms.coef_wp2sclrp_implicit) << " " << "\n"
              << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "pdf_implicit_coefs_terms%term_wp2sclrp_explicit = " << " "
              << dump_values(pdf_implicit_coefs_terms.term_wp2sclrp_explicit) << " " << "\n"
              << "\n";
      std::cerr << message.str();
    }
  }
  if (l_predict_upwp_vpwp) {
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "um_forcing = " << " " << dump_values(um_forcing) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "vm_forcing = " << " " << dump_values(vm_forcing) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "ug = " << " " << dump_values(ug) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "vg = " << " " << dump_values(vg) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "wpthvp = " << " " << dump_values(wpthvp) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "fcor = " << " " << fcor << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "um_ref = " << " " << dump_values(um_ref) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "vm_ref = " << " " << dump_values(vm_ref) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "up2 = " << " " << dump_values(up2) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "vp2 = " << " " << dump_values(vp2) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "uprcp = " << " " << dump_values(uprcp) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "vprcp = " << " " << dump_values(vprcp) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "rc_coef_zm = " << " " << dump_values(rc_coef_zm) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "pdf_implicit_coefs_terms%coef_wp2up_implicit = " << " "
              << dump_values(pdf_implicit_coefs_terms.coef_wp2up_implicit) << " " << "\n"
              << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "pdf_implicit_coefs_terms%term_wp2up_explicit = " << " "
              << dump_values(pdf_implicit_coefs_terms.term_wp2up_explicit) << " " << "\n"
              << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "pdf_implicit_coefs_terms%coef_wp2vp_implicit = " << " "
              << dump_values(pdf_implicit_coefs_terms.coef_wp2vp_implicit) << " " << "\n"
              << "\n";
      std::cerr << message.str();
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "pdf_implicit_coefs_terms%term_wp2vp_explicit = " << " "
              << dump_values(pdf_implicit_coefs_terms.term_wp2vp_explicit) << " " << "\n"
              << "\n";
      std::cerr << message.str();
    }
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "Intent(inout)" << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  if (l_lmm_stepping) {
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "rtm (pre-solve) = " << " " << dump_values(rtm_old) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "rtm = " << " " << dump_values(rtm) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  if (l_lmm_stepping) {
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "wprtp (pre-solve) = " << " " << dump_values(wprtp_old) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "wprtp = " << " " << dump_values(wprtp) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  if (l_lmm_stepping) {
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "thlm (pre-solve) = " << " " << dump_values(thlm_old) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "thlm = " << " " << dump_values(thlm) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  if (l_lmm_stepping) {
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "wpthlp (pre-solve) =" << " " << dump_values(wpthlp_old) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
  }
  {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    message << "wpthlp =" << " " << dump_values(wpthlp) << " " << "\n" << "\n";
    std::cerr << message.str();
  }
  if (sclr_dim > 0) {
    if (l_lmm_stepping) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "sclrm (pre-solve) = " << " " << dump_values(sclrm_old) << " " << "\n" << "\n";
        std::cerr << message.str();
      }
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "sclrm = " << " " << dump_values(sclrm) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    if (l_lmm_stepping) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "wpsclrp (pre-solve) = " << " " << dump_values(wpsclrp_old) << " " << "\n"
                << "\n";
        std::cerr << message.str();
      }
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "wpsclrp = " << " " << dump_values(wpsclrp) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
  }
  if (l_predict_upwp_vpwp) {
    if (l_lmm_stepping) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "um (pre-solve) = " << " " << dump_values(um_old) << " " << "\n" << "\n";
        std::cerr << message.str();
      }
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "um = " << " " << dump_values(um) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    if (l_lmm_stepping) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "upwp (pre-solve) = " << " " << dump_values(upwp_old) << " " << "\n" << "\n";
        std::cerr << message.str();
      }
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "upwp = " << " " << dump_values(upwp) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    if (l_lmm_stepping) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "vm (pre-solve) = " << " " << dump_values(vm_old) << " " << "\n" << "\n";
        std::cerr << message.str();
      }
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "vm = " << " " << dump_values(vm) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
    if (l_lmm_stepping) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "vpwp (pre-solve) = " << " " << dump_values(vpwp_old) << " " << "\n" << "\n";
        std::cerr << message.str();
      }
    }
    {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      message << "vpwp = " << " " << dump_values(vpwp) << " " << "\n" << "\n";
      std::cerr << message.str();
    }
  }
  return;
}

} // namespace
} // namespace clubb_kokkos
