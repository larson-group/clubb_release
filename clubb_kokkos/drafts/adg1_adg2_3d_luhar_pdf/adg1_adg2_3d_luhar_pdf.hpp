// ADG1, ADG2 and 3-D Luhar component PDF closures.
//
// ADG1_pdf_driver and ADG2_pdf_driver close two-dimensional column/level fields;
// Luhar_3D_pdf_driver closes one profile. ADG1_w_closure supplies the vertical-
// velocity mixture, while calc_Luhar_params and close_Luhar_pdf supply Luhar
// mixture/response parameters. Means retain field units, variances squared units,
// and normalized means/mixture/width factors are dimensionless.
// ADG1 kernels use column-first LayoutLeft execution-space Views. The three
// Luhar profile APIs currently stage arbitrary Views to host C++ calculations
// while complex calc_roots remains a transient native dependency. No stats or
// persistent/global model storage is owned here. See DRAFT_STATUS.md for scope,
// missing lower bridges, definedness deviation and absent runtime validation.
#pragma once
#include "array_types.hpp"
#include "err_info_type_module.hpp"
namespace clubb_kokkos {
using HostReal1D = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
void ADG1_pdf_driver(int nz, int ngrdcol, int sclr_dim, Real1D sclr_tol,            // In
                     Real2D wm, Real2D rtm, Real2D thlm, Real2D um, Real2D vm,      // In
                     Real2D wp2, Real2D rtp2, Real2D thlp2, Real2D up2, Real2D vp2, // In
                     Real2D Skw, Real2D wprtp, Real2D wpthlp, Real2D upwp, Real2D vpwp,
                     Real2D sqrt_wp2,                                                 // In
                     Real2D sigma_sqd_w, Real1D beta, Real1D mixt_frac_max_mag,       // In
                     Real3D sclrm, Real3D sclrp2, Real3D wpsclrp, bool l_scalar_calc, // In
                     err_info_type& err_info,                                         // InOut
                     Real2D w_1, Real2D w_2,                                          // Out
                     Real2D rt_1, Real2D rt_2,                                        // Out
                     Real2D thl_1, Real2D thl_2,                                      // Out
                     Real2D u_1, Real2D u_2, Real2D v_1, Real2D v_2,                  // Out
                     Real2D varnce_w_1, Real2D varnce_w_2,                            // Out
                     Real2D varnce_rt_1, Real2D varnce_rt_2,                          // Out
                     Real2D varnce_thl_1, Real2D varnce_thl_2,                        // Out
                     Real2D varnce_u_1, Real2D varnce_u_2,                            // Out
                     Real2D varnce_v_1, Real2D varnce_v_2,                            // Out
                     Real2D mixt_frac,                                                // Out
                     Real2D alpha_rt, Real2D alpha_thl,                               // Out
                     Real2D alpha_u, Real2D alpha_v,                                  // Out
                     Real3D sclr_1, Real3D sclr_2, Real3D varnce_sclr_1,              // Out
                     Real3D varnce_sclr_2, Real3D alpha_sclr);

void ADG2_pdf_driver(int nz, int ngrdcol, int sclr_dim, Real1D sclr_tol, // In
                     Real2D wm, Real2D rtm, Real2D thlm, Real2D wp2, Real2D rtp2,
                     Real2D thlp2,                                                          // In
                     Real2D Skw, Real2D wprtp, Real2D wpthlp, Real2D sqrt_wp2, Real1D beta, // In
                     Real3D sclrm, Real3D sclrp2, Real3D wpsclrp, bool l_scalar_calc,       // In
                     err_info_type& err_info,                                               // InOut
                     Real2D w_1, Real2D w_2,                                                // Out
                     Real2D rt_1, Real2D rt_2,                                              // Out
                     Real2D thl_1, Real2D thl_2,                                            // Out
                     Real2D varnce_w_1, Real2D varnce_w_2,                                  // Out
                     Real2D varnce_rt_1, Real2D varnce_rt_2,                                // Out
                     Real2D varnce_thl_1, Real2D varnce_thl_2,                              // Out
                     Real2D mixt_frac,                                                      // Out
                     Real2D alpha_rt, Real2D alpha_thl,                                     // Out
                     Real2D sigma_sqd_w, Real3D sclr_1, Real3D sclr_2,                      // Out
                     Real3D varnce_sclr_1, Real3D varnce_sclr_2, Real3D alpha_sclr);

void Luhar_3D_pdf_driver(int nz, HostReal1D wm, HostReal1D rtm, HostReal1D thlm, HostReal1D wp2,
                         HostReal1D rtp2, HostReal1D thlp2, // In
                         HostReal1D Skw, HostReal1D Skrt, HostReal1D Skthl, HostReal1D wprtp,
                         HostReal1D wpthlp,                                // In
                         HostReal1D w_1, HostReal1D w_2,                   // Out
                         HostReal1D rt_1, HostReal1D rt_2,                 // Out
                         HostReal1D thl_1, HostReal1D thl_2,               // Out
                         HostReal1D varnce_w_1, HostReal1D varnce_w_2,     // Out
                         HostReal1D varnce_rt_1, HostReal1D varnce_rt_2,   // Out
                         HostReal1D varnce_thl_1, HostReal1D varnce_thl_2, // Out
                         HostReal1D mixt_frac);

void ADG1_w_closure(int nz, int ngrdcol, Real2D wm, Real2D wp2, Real2D Skw,
                    Real2D sigma_sqd_w,                                 // In
                    Real2D sqrt_wp2, Real1D mixt_frac_max_mag,          // In
                    err_info_type& err_info,                            // InOut
                    Real2D w_1, Real2D w_2, Real2D w_1_n, Real2D w_2_n, // Out
                    Real2D varnce_w_1, Real2D varnce_w_2, Real2D mixt_frac);

void calc_Luhar_params(int nz, HostReal1D Skx, HostReal1D wpxp, // In
                       HostReal1D xp2, real_t x_tol_sqd,        // In
                       HostReal1D mixt_frac, HostReal1D big_m, HostReal1D small_m);

void close_Luhar_pdf(int nz, HostReal1D xm, HostReal1D xp2, HostReal1D mixt_frac, // In
                     HostReal1D small_m, HostReal1D wpxp, real_t x_tol_sqd,       // In
                     HostReal1D sigma_sqd_x_1, HostReal1D sigma_sqd_x_2,          // Out
                     HostReal1D varnce_x_1, HostReal1D varnce_x_2,                // Out
                     HostReal1D x_1_n, HostReal1D x_2_n, HostReal1D x_1, HostReal1D x_2);

template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23>
void Luhar_3D_pdf_driver(int nz, Array0 wm, Array1 rtm, Array2 thlm, Array3 wp2, Array4 rtp2,
                         Array5 thlp2, Array6 Skw, Array7 Skrt, Array8 Skthl, Array9 wprtp,
                         Array10 wpthlp, Array11 w_1, Array12 w_2, Array13 rt_1, Array14 rt_2,
                         Array15 thl_1, Array16 thl_2, Array17 varnce_w_1, Array18 varnce_w_2,
                         Array19 varnce_rt_1, Array20 varnce_rt_2, Array21 varnce_thl_1,
                         Array22 varnce_thl_2, Array23 mixt_frac) {
  // Stage profiles at the temporary host/complex-root boundary.
  auto wm_host     = to_host_space(wm);
  auto rtm_host    = to_host_space(rtm);
  auto thlm_host   = to_host_space(thlm);
  auto wp2_host    = to_host_space(wp2);
  auto rtp2_host   = to_host_space(rtp2);
  auto thlp2_host  = to_host_space(thlp2);
  auto Skw_host    = to_host_space(Skw);
  auto Skrt_host   = to_host_space(Skrt);
  auto Skthl_host  = to_host_space(Skthl);
  auto wprtp_host  = to_host_space(wprtp);
  auto wpthlp_host = to_host_space(wpthlp);
  HostReal1D w_1_host("Luhar_3D_pdf_driver::w_1", nz);
  HostReal1D w_2_host("Luhar_3D_pdf_driver::w_2", nz);
  HostReal1D rt_1_host("Luhar_3D_pdf_driver::rt_1", nz);
  HostReal1D rt_2_host("Luhar_3D_pdf_driver::rt_2", nz);
  HostReal1D thl_1_host("Luhar_3D_pdf_driver::thl_1", nz);
  HostReal1D thl_2_host("Luhar_3D_pdf_driver::thl_2", nz);
  HostReal1D varnce_w_1_host("Luhar_3D_pdf_driver::varnce_w_1", nz);
  HostReal1D varnce_w_2_host("Luhar_3D_pdf_driver::varnce_w_2", nz);
  HostReal1D varnce_rt_1_host("Luhar_3D_pdf_driver::varnce_rt_1", nz);
  HostReal1D varnce_rt_2_host("Luhar_3D_pdf_driver::varnce_rt_2", nz);
  HostReal1D varnce_thl_1_host("Luhar_3D_pdf_driver::varnce_thl_1", nz);
  HostReal1D varnce_thl_2_host("Luhar_3D_pdf_driver::varnce_thl_2", nz);
  HostReal1D mixt_frac_host("Luhar_3D_pdf_driver::mixt_frac", nz);
  Luhar_3D_pdf_driver(nz, wm_host, rtm_host, thlm_host, wp2_host, rtp2_host, thlp2_host, Skw_host,
                      Skrt_host, Skthl_host, wprtp_host, wpthlp_host, w_1_host, w_2_host, rt_1_host,
                      rt_2_host, thl_1_host, thl_2_host, varnce_w_1_host, varnce_w_2_host,
                      varnce_rt_1_host, varnce_rt_2_host, varnce_thl_1_host, varnce_thl_2_host,
                      mixt_frac_host);
  from_execution_space(w_1, w_1_host);
  from_execution_space(w_2, w_2_host);
  from_execution_space(rt_1, rt_1_host);
  from_execution_space(rt_2, rt_2_host);
  from_execution_space(thl_1, thl_1_host);
  from_execution_space(thl_2, thl_2_host);
  from_execution_space(varnce_w_1, varnce_w_1_host);
  from_execution_space(varnce_w_2, varnce_w_2_host);
  from_execution_space(varnce_rt_1, varnce_rt_1_host);
  from_execution_space(varnce_rt_2, varnce_rt_2_host);
  from_execution_space(varnce_thl_1, varnce_thl_1_host);
  from_execution_space(varnce_thl_2, varnce_thl_2_host);
  from_execution_space(mixt_frac, mixt_frac_host);
}

template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5>
void calc_Luhar_params(int nz, Array0 Skx, Array1 wpxp, Array2 xp2, real_t x_tol_sqd,
                       Array3 mixt_frac, Array4 big_m, Array5 small_m) {
  // Stage profiles at the temporary host/complex-root boundary.
  auto Skx_host  = to_host_space(Skx);
  auto wpxp_host = to_host_space(wpxp);
  auto xp2_host  = to_host_space(xp2);
  HostReal1D mixt_frac_host("calc_Luhar_params::mixt_frac", nz);
  HostReal1D big_m_host("calc_Luhar_params::big_m", nz);
  HostReal1D small_m_host("calc_Luhar_params::small_m", nz);
  calc_Luhar_params(nz, Skx_host, wpxp_host, xp2_host, x_tol_sqd, mixt_frac_host, big_m_host,
                    small_m_host);
  from_execution_space(mixt_frac, mixt_frac_host);
  from_execution_space(big_m, big_m_host);
  from_execution_space(small_m, small_m_host);
}

template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12>
void close_Luhar_pdf(int nz, Array0 xm, Array1 xp2, Array2 mixt_frac, Array3 small_m, Array4 wpxp,
                     real_t x_tol_sqd, Array5 sigma_sqd_x_1, Array6 sigma_sqd_x_2,
                     Array7 varnce_x_1, Array8 varnce_x_2, Array9 x_1_n, Array10 x_2_n, Array11 x_1,
                     Array12 x_2) {
  // Stage profiles at the temporary host/complex-root boundary.
  auto xm_host        = to_host_space(xm);
  auto xp2_host       = to_host_space(xp2);
  auto mixt_frac_host = to_host_space(mixt_frac);
  auto small_m_host   = to_host_space(small_m);
  auto wpxp_host      = to_host_space(wpxp);
  HostReal1D sigma_sqd_x_1_host("close_Luhar_pdf::sigma_sqd_x_1", nz);
  HostReal1D sigma_sqd_x_2_host("close_Luhar_pdf::sigma_sqd_x_2", nz);
  HostReal1D varnce_x_1_host("close_Luhar_pdf::varnce_x_1", nz);
  HostReal1D varnce_x_2_host("close_Luhar_pdf::varnce_x_2", nz);
  HostReal1D x_1_n_host("close_Luhar_pdf::x_1_n", nz);
  HostReal1D x_2_n_host("close_Luhar_pdf::x_2_n", nz);
  HostReal1D x_1_host("close_Luhar_pdf::x_1", nz);
  HostReal1D x_2_host("close_Luhar_pdf::x_2", nz);
  close_Luhar_pdf(nz, xm_host, xp2_host, mixt_frac_host, small_m_host, wpxp_host, x_tol_sqd,
                  sigma_sqd_x_1_host, sigma_sqd_x_2_host, varnce_x_1_host, varnce_x_2_host,
                  x_1_n_host, x_2_n_host, x_1_host, x_2_host);
  from_execution_space(sigma_sqd_x_1, sigma_sqd_x_1_host);
  from_execution_space(sigma_sqd_x_2, sigma_sqd_x_2_host);
  from_execution_space(varnce_x_1, varnce_x_1_host);
  from_execution_space(varnce_x_2, varnce_x_2_host);
  from_execution_space(x_1_n, x_1_n_host);
  from_execution_space(x_2_n, x_2_n_host);
  from_execution_space(x_1, x_1_host);
  from_execution_space(x_2, x_2_host);
}

} // namespace clubb_kokkos
