// ADG1, ADG2 and 3-D Luhar component PDF closures.
//
// ADG1_pdf_driver and ADG2_pdf_driver close two-dimensional column/level fields;
// Luhar_3D_pdf_driver closes one profile. ADG1_w_closure supplies the vertical-
// velocity mixture, while calc_Luhar_params and close_Luhar_pdf supply Luhar
// mixture/response parameters. Means retain field units, variances squared units,
// and normalized means/mixture/width factors are dimensionless.
// Kernels use column-first LayoutLeft execution-space Views. Luhar profile
// adapters stage arbitrary Views through the shared array owner; independent
// levels and polynomial roots are evaluated on the execution space. No stats
// or persistent/global model storage is owned here. close_Luhar_pdf documents
// the defined covariance sign used for Native's unwritten nonvarying sign.
#pragma once
#include "scratch_workspace.hpp"
#include "array_types.hpp"
#include "err_info_type_module.hpp"
namespace clubb_kokkos {
void ADG1_pdf_driver(ScratchWorkspace& workspace,                                   // InOut
                     int nz, int ngrdcol, int sclr_dim, Real1D sclr_tol,            // In
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

void ADG2_pdf_driver(ScratchWorkspace& workspace,                        // InOut
                     int nz, int ngrdcol, int sclr_dim, Real1D sclr_tol, // In
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

void Luhar_3D_pdf_driver(int nz, Real1D wm, Real1D rtm, Real1D thlm, Real1D wp2, Real1D rtp2,
                         Real1D thlp2, // In
                         Real1D Skw, Real1D Skrt, Real1D Skthl, Real1D wprtp,
                         Real1D wpthlp,                            // In
                         Real1D w_1, Real1D w_2,                   // Out
                         Real1D rt_1, Real1D rt_2,                 // Out
                         Real1D thl_1, Real1D thl_2,               // Out
                         Real1D varnce_w_1, Real1D varnce_w_2,     // Out
                         Real1D varnce_rt_1, Real1D varnce_rt_2,   // Out
                         Real1D varnce_thl_1, Real1D varnce_thl_2, // Out
                         Real1D mixt_frac);

void ADG1_w_closure(int nz, int ngrdcol, Real2D wm, Real2D wp2, Real2D Skw,
                    Real2D sigma_sqd_w,                                 // In
                    Real2D sqrt_wp2, Real1D mixt_frac_max_mag,          // In
                    err_info_type& err_info,                            // InOut
                    Real2D w_1, Real2D w_2, Real2D w_1_n, Real2D w_2_n, // Out
                    Real2D varnce_w_1, Real2D varnce_w_2, Real2D mixt_frac);

void calc_Luhar_params(int nz, Real1D Skx, Real1D wpxp, // In
                       Real1D xp2, real_t x_tol_sqd,    // In
                       Real1D mixt_frac, Real1D big_m, Real1D small_m);

void close_Luhar_pdf(int nz, Real1D xm, Real1D xp2, Real1D mixt_frac, // In
                     Real1D small_m, Real1D wpxp, real_t x_tol_sqd,   // In
                     Real1D sigma_sqd_x_1, Real1D sigma_sqd_x_2,      // Out
                     Real1D varnce_x_1, Real1D varnce_x_2,            // Out
                     Real1D x_1_n, Real1D x_2_n, Real1D x_1, Real1D x_2);

template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23>
void Luhar_3D_pdf_driver(ScratchWorkspace& workspace, int nz, Array0 wm, Array1 rtm, Array2 thlm,
                         Array3 wp2, Array4 rtp2, Array5 thlp2, Array6 Skw, Array7 Skrt,
                         Array8 Skthl, Array9 wprtp, Array10 wpthlp, Array11 w_1, Array12 w_2,
                         Array13 rt_1, Array14 rt_2, Array15 thl_1, Array16 thl_2,
                         Array17 varnce_w_1, Array18 varnce_w_2, Array19 varnce_rt_1,
                         Array20 varnce_rt_2, Array21 varnce_thl_1, Array22 varnce_thl_2,
                         Array23 mixt_frac) {
  // Preserve Native profile arguments while reusing caller-owned staging.
  // The scalar closure writes every output level in every setter/tolerance branch.
  auto staging_scope      = workspace.scope();
  auto wm_device          = to_execution_space(workspace, wm);
  auto rtm_device         = to_execution_space(workspace, rtm);
  auto thlm_device        = to_execution_space(workspace, thlm);
  auto wp2_device         = to_execution_space(workspace, wp2);
  auto rtp2_device        = to_execution_space(workspace, rtp2);
  auto thlp2_device       = to_execution_space(workspace, thlp2);
  auto Skw_device         = to_execution_space(workspace, Skw);
  auto Skrt_device        = to_execution_space(workspace, Skrt);
  auto Skthl_device       = to_execution_space(workspace, Skthl);
  auto wprtp_device       = to_execution_space(workspace, wprtp);
  auto wpthlp_device      = to_execution_space(workspace, wpthlp);
  auto w_1_device         = to_execution_space(workspace, Kokkos::WithoutInitializing, w_1);
  auto w_2_device         = to_execution_space(workspace, Kokkos::WithoutInitializing, w_2);
  auto rt_1_device        = to_execution_space(workspace, Kokkos::WithoutInitializing, rt_1);
  auto rt_2_device        = to_execution_space(workspace, Kokkos::WithoutInitializing, rt_2);
  auto thl_1_device       = to_execution_space(workspace, Kokkos::WithoutInitializing, thl_1);
  auto thl_2_device       = to_execution_space(workspace, Kokkos::WithoutInitializing, thl_2);
  auto varnce_w_1_device  = to_execution_space(workspace, Kokkos::WithoutInitializing, varnce_w_1);
  auto varnce_w_2_device  = to_execution_space(workspace, Kokkos::WithoutInitializing, varnce_w_2);
  auto varnce_rt_1_device = to_execution_space(workspace, Kokkos::WithoutInitializing, varnce_rt_1);
  auto varnce_rt_2_device = to_execution_space(workspace, Kokkos::WithoutInitializing, varnce_rt_2);
  auto varnce_thl_1_device =
      to_execution_space(workspace, Kokkos::WithoutInitializing, varnce_thl_1);
  auto varnce_thl_2_device =
      to_execution_space(workspace, Kokkos::WithoutInitializing, varnce_thl_2);
  auto mixt_frac_device = to_execution_space(workspace, Kokkos::WithoutInitializing, mixt_frac);
  Luhar_3D_pdf_driver(nz, wm_device, rtm_device, thlm_device, wp2_device, rtp2_device, thlp2_device,
                      Skw_device, Skrt_device, Skthl_device, wprtp_device, wpthlp_device,
                      w_1_device, w_2_device, rt_1_device, rt_2_device, thl_1_device, thl_2_device,
                      varnce_w_1_device, varnce_w_2_device, varnce_rt_1_device, varnce_rt_2_device,
                      varnce_thl_1_device, varnce_thl_2_device, mixt_frac_device);
  from_execution_space(workspace, w_1, w_1_device);
  from_execution_space(workspace, w_2, w_2_device);
  from_execution_space(workspace, rt_1, rt_1_device);
  from_execution_space(workspace, rt_2, rt_2_device);
  from_execution_space(workspace, thl_1, thl_1_device);
  from_execution_space(workspace, thl_2, thl_2_device);
  from_execution_space(workspace, varnce_w_1, varnce_w_1_device);
  from_execution_space(workspace, varnce_w_2, varnce_w_2_device);
  from_execution_space(workspace, varnce_rt_1, varnce_rt_1_device);
  from_execution_space(workspace, varnce_rt_2, varnce_rt_2_device);
  from_execution_space(workspace, varnce_thl_1, varnce_thl_1_device);
  from_execution_space(workspace, varnce_thl_2, varnce_thl_2_device);
  from_execution_space(workspace, mixt_frac, mixt_frac_device);
}

template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5>
void calc_Luhar_params(int nz, Array0 Skx, Array1 wpxp, Array2 xp2, real_t x_tol_sqd,
                       Array3 mixt_frac, Array4 big_m, Array5 small_m) {
  // Preserve Native profile arguments at the shared memory-space boundary.
  auto Skx_device       = to_execution_space(Skx);
  auto wpxp_device      = to_execution_space(wpxp);
  auto xp2_device       = to_execution_space(xp2);
  auto mixt_frac_device = to_execution_space(mixt_frac, false);
  auto big_m_device     = to_execution_space(big_m, false);
  auto small_m_device   = to_execution_space(small_m, false);
  calc_Luhar_params(nz, Skx_device, wpxp_device, xp2_device, x_tol_sqd, mixt_frac_device,
                    big_m_device, small_m_device);
  from_execution_space(mixt_frac, mixt_frac_device);
  from_execution_space(big_m, big_m_device);
  from_execution_space(small_m, small_m_device);
}

template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12>
void close_Luhar_pdf(int nz, Array0 xm, Array1 xp2, Array2 mixt_frac, Array3 small_m, Array4 wpxp,
                     real_t x_tol_sqd, Array5 sigma_sqd_x_1, Array6 sigma_sqd_x_2,
                     Array7 varnce_x_1, Array8 varnce_x_2, Array9 x_1_n, Array10 x_2_n, Array11 x_1,
                     Array12 x_2) {
  // Preserve Native profile arguments at the shared memory-space boundary.
  auto xm_device            = to_execution_space(xm);
  auto xp2_device           = to_execution_space(xp2);
  auto mixt_frac_device     = to_execution_space(mixt_frac);
  auto small_m_device       = to_execution_space(small_m);
  auto wpxp_device          = to_execution_space(wpxp);
  auto sigma_sqd_x_1_device = to_execution_space(sigma_sqd_x_1, false);
  auto sigma_sqd_x_2_device = to_execution_space(sigma_sqd_x_2, false);
  auto varnce_x_1_device    = to_execution_space(varnce_x_1, false);
  auto varnce_x_2_device    = to_execution_space(varnce_x_2, false);
  auto x_1_n_device         = to_execution_space(x_1_n, false);
  auto x_2_n_device         = to_execution_space(x_2_n, false);
  auto x_1_device           = to_execution_space(x_1, false);
  auto x_2_device           = to_execution_space(x_2, false);
  close_Luhar_pdf(nz, xm_device, xp2_device, mixt_frac_device, small_m_device, wpxp_device,
                  x_tol_sqd, sigma_sqd_x_1_device, sigma_sqd_x_2_device, varnce_x_1_device,
                  varnce_x_2_device, x_1_n_device, x_2_n_device, x_1_device, x_2_device);
  from_execution_space(sigma_sqd_x_1, sigma_sqd_x_1_device);
  from_execution_space(sigma_sqd_x_2, sigma_sqd_x_2_device);
  from_execution_space(varnce_x_1, varnce_x_1_device);
  from_execution_space(varnce_x_2, varnce_x_2_device);
  from_execution_space(x_1_n, x_1_n_device);
  from_execution_space(x_2_n, x_2_n_device);
  from_execution_space(x_1, x_1_device);
  from_execution_space(x_2, x_2_device);
}

} // namespace clubb_kokkos
