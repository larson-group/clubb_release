// Clip moments and cloud water while retaining Native budgets and diagnostics.
//
// Covariance clipping preserves boundary fluxes; variance clipping retains the
// optional scalar maximum and Native termination guard. Skewness clipping keeps
// both height-dependent limiters and the absolute wp3 cap. All physical fields
// use column-first LayoutLeft Views and caller-owned outputs; scratch is local.
// Statistics remains Native. Host warning/termination boundaries are explicit.
#include "clip_explicit.hpp"
#include "constants_clubb.hpp"
#include "advance_helper_module.hpp" // Reuse the stable Native helper owner at promotion.
#include "f2k/procedure_bridge.hpp"
#include "f2k_dependencies.hpp"
#include <Kokkos_MathematicalFunctions.hpp>
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
namespace clubb_kokkos {
//=============================================================================
// rtp2: r_t'^2         [(kg/kg)^2]
// thlp2: theta_l'^2     [K^2]
// up2: u'^2           [m^2/s^2]
// vp2: v'^2           [m^2/s^2]
// wp2: w'^2           [m^2/s^2]
// sclrp2: sclr'^2  [{units vary}^2]
// l_linearize_pbl_winds: Flag (used by E3SM) to linearize PBL winds
// l_predict_upwp_vpwp: Prognose u'w' and v'w'
// wprtp: w'r_t'        [(kg/kg) m/s]
// wpthlp: w'theta_l'    [K m/s]
// upwp: u'w'          [m^2/s^2]
// vpwp: v'w'          [m^2/s^2]
// wpsclrp: w'sclr'         [units m/s]
// upwp_pert: perturbed <u'w'> [m^2/s^2]
// vpwp_pert: perturbed <v'w'> [m^2/s^2]
void clip_covars_denom(int nzm, int ngrdcol, int sclr_dim,                            // In
                       real_t dt,                                                     // In
                       Real2D rtp2, Real2D thlp2, Real2D up2, Real2D vp2, Real2D wp2, // In
                       Real3D sclrp2,                                                 // In
                       bool l_linearize_pbl_winds,                                    // In
                       bool l_predict_upwp_vpwp,                                      // In
                       f2k::stats_type& stats,                                        // InOut
                       int& wprtp_cl_num, int& wpthlp_cl_num, int& upwp_cl_num,
                       int& vpwp_cl_num, // InOut
                       Real2D wprtp, Real2D wpthlp, Real2D upwp, Real2D vpwp,
                       Real3D wpsclrp,                     // InOut
                       Real2D upwp_pert, Real2D vpwp_pert) // InOut
{
  // Description:
  // Some of the covariances found in the CLUBB model code need to be clipped
  // multiple times during each timestep to ensure that the correlation between
  // the two relevant variables stays between -1 and 1 at all times during the
  // model run.  The covariances that need to be clipped multiple times are
  // w'r_t', w'th_l', w'sclr', u'w', and v'w'.  One of the times that each one
  // of these covariances is clipped is immediately after each one is set.
  // However, each covariance still needs to be clipped two more times during
  // each timestep (once after advance_xp2_xpyp is called and once after
  // advance_wp2_wp3 is called).  This subroutine handles the times that the
  // covariances are clipped away from the time that they are set.  In other
  // words, this subroutine clips the covariances after the denominator terms
  // in the relevant correlation equation have been altered, ensuring that
  // all correlations will remain between -1 and 1 at all times.
  // References:
  // None
  //-----------------------------------------------------------------------
  // --------------------- Input Variables ---------------------
  // --------------------- Input/Output Variables ---------------------
  // Variables used to track perturbed version of winds.
  // --------------------- Local Variables ---------------------
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wprtp_chnge(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::wprtp_chnge"), ngrdcol,
      nzm); // Net change in w'r_t' due to clipping  [(kg/kg) m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wpthlp_chnge(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::wpthlp_chnge"), ngrdcol,
      nzm); // Net change in w'th_l' due to clipping [K m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> upwp_chnge(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::upwp_chnge"), ngrdcol,
      nzm); // Net change in u'w' due to clipping    [m^2/s^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vpwp_chnge(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::vpwp_chnge"), ngrdcol,
      nzm); // Net change in v'w' due to clipping    [m^2/s^2]
  Kokkos::View<real_t***, Kokkos::LayoutLeft> wpsclrp_chnge(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::wpsclrp_chnge"), ngrdcol,
      nzm, sclr_dim); // Net change in w'sclr' due to clipping [{units vary}]
  int sclr;           // scalar array index.

  // --------------------- Begin Code ---------------------

  if (stats.flags().l_sample) {

    if (wprtp_cl_num == 0) {
      {
        Real2D stats_values_1(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::stats_values_1"),
            wprtp.extent(0), wprtp.extent(1));
        Kokkos::parallel_for(
            "clip_covars_denom::stats_value", wprtp.size(), KOKKOS_LAMBDA(int cell) {
              const int i          = cell % wprtp.extent(0);
              const int k          = cell / wprtp.extent(0);
              stats_values_1(i, k) = wprtp(i, k) / dt;
            });
        f2k::stats_begin_budget("wprtp_cl", stats_values_1, // In
                                stats);                     // InOut
      }
    } else {
      {
        Real2D stats_values_2(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::stats_values_2"),
            wprtp.extent(0), wprtp.extent(1));
        Kokkos::parallel_for(
            "clip_covars_denom::stats_value", wprtp.size(), KOKKOS_LAMBDA(int cell) {
              const int i          = cell % wprtp.extent(0);
              const int k          = cell / wprtp.extent(0);
              stats_values_2(i, k) = -wprtp(i, k) / dt;
            });
        f2k::stats_update_budget("wprtp_cl", stats_values_2, // In
                                 stats);                     // InOut
      }
    }

    if (wpthlp_cl_num == 0) {
      {
        Real2D stats_values_3(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::stats_values_3"),
            wpthlp.extent(0), wpthlp.extent(1));
        Kokkos::parallel_for(
            "clip_covars_denom::stats_value", wpthlp.size(), KOKKOS_LAMBDA(int cell) {
              const int i          = cell % wpthlp.extent(0);
              const int k          = cell / wpthlp.extent(0);
              stats_values_3(i, k) = wpthlp(i, k) / dt;
            });
        f2k::stats_begin_budget("wpthlp_cl", stats_values_3, // In
                                stats);                      // InOut
      }
    } else {
      {
        Real2D stats_values_4(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::stats_values_4"),
            wpthlp.extent(0), wpthlp.extent(1));
        Kokkos::parallel_for(
            "clip_covars_denom::stats_value", wpthlp.size(), KOKKOS_LAMBDA(int cell) {
              const int i          = cell % wpthlp.extent(0);
              const int k          = cell / wpthlp.extent(0);
              stats_values_4(i, k) = -wpthlp(i, k) / dt;
            });
        f2k::stats_update_budget("wpthlp_cl", stats_values_4, // In
                                 stats);                      // InOut
      }
    }

    if (l_predict_upwp_vpwp) {
      if (upwp_cl_num == 0) {
        {
          Real2D stats_values_5(
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::stats_values_5"),
              upwp.extent(0), upwp.extent(1));
          Kokkos::parallel_for(
              "clip_covars_denom::stats_value", upwp.size(), KOKKOS_LAMBDA(int cell) {
                const int i          = cell % upwp.extent(0);
                const int k          = cell / upwp.extent(0);
                stats_values_5(i, k) = upwp(i, k) / dt;
              });
          f2k::stats_begin_budget("upwp_cl", stats_values_5, // In
                                  stats);                    // InOut
        }
      } else {
        {
          Real2D stats_values_6(
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::stats_values_6"),
              upwp.extent(0), upwp.extent(1));
          Kokkos::parallel_for(
              "clip_covars_denom::stats_value", upwp.size(), KOKKOS_LAMBDA(int cell) {
                const int i          = cell % upwp.extent(0);
                const int k          = cell / upwp.extent(0);
                stats_values_6(i, k) = -upwp(i, k) / dt;
              });
          f2k::stats_update_budget("upwp_cl", stats_values_6, // In
                                   stats);                    // InOut
        }
      }

      if (vpwp_cl_num == 0) {
        {
          Real2D stats_values_7(
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::stats_values_7"),
              vpwp.extent(0), vpwp.extent(1));
          Kokkos::parallel_for(
              "clip_covars_denom::stats_value", vpwp.size(), KOKKOS_LAMBDA(int cell) {
                const int i          = cell % vpwp.extent(0);
                const int k          = cell / vpwp.extent(0);
                stats_values_7(i, k) = vpwp(i, k) / dt;
              });
          f2k::stats_begin_budget("vpwp_cl", stats_values_7, // In
                                  stats);                    // InOut
        }
      } else {
        {
          Real2D stats_values_8(
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::stats_values_8"),
              vpwp.extent(0), vpwp.extent(1));
          Kokkos::parallel_for(
              "clip_covars_denom::stats_value", vpwp.size(), KOKKOS_LAMBDA(int cell) {
                const int i          = cell % vpwp.extent(0);
                const int k          = cell / vpwp.extent(0);
                stats_values_8(i, k) = -vpwp(i, k) / dt;
              });
          f2k::stats_update_budget("vpwp_cl", stats_values_8, // In
                                   stats);                    // InOut
        }
      }
    }
  }

  wprtp_cl_num  = wprtp_cl_num + 1;
  wpthlp_cl_num = wpthlp_cl_num + 1;
  upwp_cl_num   = upwp_cl_num + 1;
  vpwp_cl_num   = vpwp_cl_num + 1;

  //!! Clipping for w'r_t'
  //
  // Clipping w'r_t' at each vertical level, based on the
  // correlation of w and r_t at each vertical level, such that:
  // corr_(w,r_t) = w'r_t' / [ sqrt(w'^2) * sqrt(r_t'^2) ];
  // -1 <= corr_(w,r_t) <= 1.
  //
  // Since w'^2, r_t'^2, and w'r_t' are each advanced in different
  // subroutines from each other in advance_clubb_core, clipping for w'r_t'
  // is done three times during each timestep (once after each variable has
  // been updated).
  //
  // This subroutine handles the first and third instances of
  // w'r_t' clipping.
  // The first instance of w'r_t' clipping takes place after
  // r_t'^2 is updated in advance_xp2_xpyp.
  // The third instance of w'r_t' clipping takes place after
  // w'^2 is updated in advance_wp2_wp3.

  // Clip w'r_t'
  clip_covar(nzm, ngrdcol, clip_wprtp, wp2, rtp2, // In
             wprtp,                               // InOut
             wprtp_chnge);                        // Out

  // Clip w'th_l'
  clip_covar(nzm, ngrdcol, clip_wpthlp, wp2, thlp2, // In
             wpthlp,                                // InOut
             wpthlp_chnge);                         // Out

  //!! Clipping for w'sclr'
  //
  // Clipping w'sclr' at each vertical level, based on the
  // correlation of w and sclr at each vertical level, such that:
  // corr_(w,sclr) = w'sclr' / [ sqrt(w'^2) * sqrt(sclr'^2) ];
  // -1 <= corr_(w,sclr) <= 1.
  //
  // Since w'^2, sclr'^2, and w'sclr' are each advanced in different
  // subroutines from each other in advance_clubb_core, clipping for w'sclr'
  // is done three times during each timestep (once after each variable has
  // been updated).
  //
  // This subroutine handles the first and third instances of
  // w'sclr' clipping.
  // The first instance of w'sclr' clipping takes place after
  // sclr'^2 is updated in advance_xp2_xpyp.
  // The third instance of w'sclr' clipping takes place after
  // w'^2 is updated in advance_wp2_wp3.

  // Clip w'sclr'
  for (sclr = 1; sclr <= sclr_dim; sclr += 1) {
    clip_covar(nzm, ngrdcol, clip_wpsclrp, Kokkos::subview(wp2, Kokkos::ALL, Kokkos::ALL), // In
               Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),              // In
               Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),             // InOut
               Kokkos::subview(wpsclrp_chnge, Kokkos::ALL, Kokkos::ALL, (sclr - 1)));      // Out
  }

  //!! Clipping for u'w'
  //
  // Clipping u'w' at each vertical level, based on the
  // correlation of u and w at each vertical level, such that:
  // corr_(u,w) = u'w' / [ sqrt(u'^2) * sqrt(w'^2) ];
  // -1 <= corr_(u,w) <= 1.
  //
  // Since w'^2, u'^2, and u'w' are each advanced in different
  // subroutines from each other in advance_clubb_core, clipping for u'w'
  // is done three times during each timestep (once after each variable has
  // been updated).
  //
  // This subroutine handles the first and second instances of
  // u'w' clipping.
  // The first instance of u'w' clipping takes place after
  // u'^2 is updated in advance_xp2_xpyp.
  // The second instance of u'w' clipping takes place after
  // w'^2 is updated in advance_wp2_wp3.

  // Clip u'w'
  clip_covar(nzm, ngrdcol, clip_upwp, wp2, up2, // In
             upwp,                              // InOut
             upwp_chnge);                       // Out

  if (l_linearize_pbl_winds) {
    clip_covar(nzm, ngrdcol, clip_upwp, wp2, up2, // In
               upwp_pert,                         // InOut
               upwp_chnge);                       // Out
    // l_linearize_pbl_winds
  }

  //!! Clipping for v'w'
  //
  // Clipping v'w' at each vertical level, based on the
  // correlation of v and w at each vertical level, such that:
  // corr_(v,w) = v'w' / [ sqrt(v'^2) * sqrt(w'^2) ];
  // -1 <= corr_(v,w) <= 1.
  //
  // Since w'^2, v'^2, and v'w' are each advanced in different
  // subroutines from each other in advance_clubb_core, clipping for v'w'
  // is done three times during each timestep (once after each variable has
  // been updated).
  //
  // This subroutine handles the first and second instances of
  // v'w' clipping.
  // The first instance of v'w' clipping takes place after
  // v'^2 is updated in advance_xp2_xpyp.
  // The second instance of v'w' clipping takes place after
  // w'^2 is updated in advance_wp2_wp3.

  clip_covar(nzm, ngrdcol, clip_vpwp, wp2, vp2, // In
             vpwp,                              // InOut
             vpwp_chnge);                       // Out

  if (l_linearize_pbl_winds) {
    clip_covar(nzm, ngrdcol, clip_vpwp, wp2, vp2, // In
               vpwp_pert,                         // InOut
               vpwp_chnge);                       // Out
    // l_linearize_pbl_winds
  }

  if (stats.flags().l_sample) {

    if (wprtp_cl_num == wprtp_cl_max) {
      {
        Real2D stats_values_9(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::stats_values_9"),
            wprtp.extent(0), wprtp.extent(1));
        Kokkos::parallel_for(
            "clip_covars_denom::stats_value", wprtp.size(), KOKKOS_LAMBDA(int cell) {
              const int i          = cell % wprtp.extent(0);
              const int k          = cell / wprtp.extent(0);
              stats_values_9(i, k) = wprtp(i, k) / dt;
            });
        f2k::stats_finalize_budget("wprtp_cl", stats_values_9, // In
                                   stats);                     // InOut
      }
    } else {
      {
        Real2D stats_values_10(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::stats_values_10"),
            wprtp.extent(0), wprtp.extent(1));
        Kokkos::parallel_for(
            "clip_covars_denom::stats_value", wprtp.size(), KOKKOS_LAMBDA(int cell) {
              const int i           = cell % wprtp.extent(0);
              const int k           = cell / wprtp.extent(0);
              stats_values_10(i, k) = wprtp(i, k) / dt;
            });
        f2k::stats_update_budget("wprtp_cl", stats_values_10, // In
                                 stats);                      // InOut
      }
    }

    if (wpthlp_cl_num == wpthlp_cl_max) {
      {
        Real2D stats_values_11(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::stats_values_11"),
            wpthlp.extent(0), wpthlp.extent(1));
        Kokkos::parallel_for(
            "clip_covars_denom::stats_value", wpthlp.size(), KOKKOS_LAMBDA(int cell) {
              const int i           = cell % wpthlp.extent(0);
              const int k           = cell / wpthlp.extent(0);
              stats_values_11(i, k) = wpthlp(i, k) / dt;
            });
        f2k::stats_finalize_budget("wpthlp_cl", stats_values_11, // In
                                   stats);                       // InOut
      }
    } else {
      {
        Real2D stats_values_12(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::stats_values_12"),
            wpthlp.extent(0), wpthlp.extent(1));
        Kokkos::parallel_for(
            "clip_covars_denom::stats_value", wpthlp.size(), KOKKOS_LAMBDA(int cell) {
              const int i           = cell % wpthlp.extent(0);
              const int k           = cell / wpthlp.extent(0);
              stats_values_12(i, k) = wpthlp(i, k) / dt;
            });
        f2k::stats_update_budget("wpthlp_cl", stats_values_12, // In
                                 stats);                       // InOut
      }
    }

    if (l_predict_upwp_vpwp) {
      if (upwp_cl_num == upwp_cl_max) {
        {
          Real2D stats_values_13(
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::stats_values_13"),
              upwp.extent(0), upwp.extent(1));
          Kokkos::parallel_for(
              "clip_covars_denom::stats_value", upwp.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = cell % upwp.extent(0);
                const int k           = cell / upwp.extent(0);
                stats_values_13(i, k) = upwp(i, k) / dt;
              });
          f2k::stats_finalize_budget("upwp_cl", stats_values_13, // In
                                     stats);                     // InOut
        }
      } else {
        {
          Real2D stats_values_14(
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::stats_values_14"),
              upwp.extent(0), upwp.extent(1));
          Kokkos::parallel_for(
              "clip_covars_denom::stats_value", upwp.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = cell % upwp.extent(0);
                const int k           = cell / upwp.extent(0);
                stats_values_14(i, k) = upwp(i, k) / dt;
              });
          f2k::stats_update_budget("upwp_cl", stats_values_14, // In
                                   stats);                     // InOut
        }
      }

      if (vpwp_cl_num == vpwp_cl_max) {
        {
          Real2D stats_values_15(
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::stats_values_15"),
              vpwp.extent(0), vpwp.extent(1));
          Kokkos::parallel_for(
              "clip_covars_denom::stats_value", vpwp.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = cell % vpwp.extent(0);
                const int k           = cell / vpwp.extent(0);
                stats_values_15(i, k) = vpwp(i, k) / dt;
              });
          f2k::stats_finalize_budget("vpwp_cl", stats_values_15, // In
                                     stats);                     // InOut
        }
      } else {
        {
          Real2D stats_values_16(
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_covars_denom::stats_values_16"),
              vpwp.extent(0), vpwp.extent(1));
          Kokkos::parallel_for(
              "clip_covars_denom::stats_value", vpwp.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = cell % vpwp.extent(0);
                const int k           = cell / vpwp.extent(0);
                stats_values_16(i, k) = vpwp(i, k) / dt;
              });
          f2k::stats_update_budget("vpwp_cl", stats_values_16, // In
                                   stats);                     // InOut
        }
      }
    }
  }

  return;
}

//=============================================================================
// solve_type: Variable being solved; used for STATS.
// xp2: Variance of x, x'^2 (momentum levels)         [{x units}^2]
// yp2: Variance of y, y'^2 (momentum levels)         [{y units}^2]
// xpyp: Covariance of x and y, x'y' (momentum levels) [{x units}*{y units}]
// xpyp_chnge: Net change in x'y' due to clipping [{x units}*{y units}]
void clip_covar(int nzm, int ngrdcol, int solve_type, Real2D xp2, Real2D yp2, // In
                Real2D xpyp,                                                  // InOut
                Real2D xpyp_chnge)                                            // Out
{
  // Description:
  // Clipping the value of covariance x'y' based on the correlation between x
  // and y.
  // and y.  This is all done at a single vertical level.
  //
  // The correlation between variables x and y is:
  //
  // corr_(x,y) = x'y' / [ sqrt(x'^2) * sqrt(y'^2) ];
  //
  // where x'^2 is the variance of x, y'^2 is the variance of y, and x'y' is
  // the covariance of x and y.
  //
  // The correlation of two variables must always have a value between -1
  // and 1, such that:
  //
  // -1 <= corr_(x,y) <= 1.
  //
  // Therefore, there is an upper limit on x'y', such that:
  //
  // x'y' <=  [ sqrt(x'^2) * sqrt(y'^2) ];
  //
  // and a lower limit on x'y', such that:
  //
  // x'y' >= -[ sqrt(x'^2) * sqrt(y'^2) ].
  //
  // The values of x'y', x'^2, and y'^2 are all found on momentum levels.
  //
  // The value of x'y' may need to be clipped whenever x'y', x'^2, or y'^2 is
  // updated.
  //
  // The following covariances are found in the code:
  //
  // w'r_t', w'th_l', w'sclr', (computed in advance_xm_wpxp);
  // r_t'th_l', sclr'r_t', sclr'th_l', (computed in advance_xp2_xpyp);
  // u'w', v'w', w'edsclr' (computed in advance_windm_edsclrm);
  // and w'hm' (computed in setup_pdf_parameters).
  // References:
  // None
  //-----------------------------------------------------------------------
  // -------------------------- Input Variables --------------------------
  // ------------------------- InOut Variable -------------------------
  //-------------------------- Output Variable --------------------------
  // -------------------------- Local Variables --------------------------
  real_t max_mag_corr; // Maximum magnitude of a correlation allowed
  real_t xpyp_bound;
  int i;
  int k; // Array index

  // -------------------------- Begin Code --------------------------

  // When clipping for wprtp or wpthlp, use the special value for
  // max_mag_correlation_flux.  For all other correlations, use
  // max_mag_correlation.
  if ((solve_type == clip_wprtp) || (solve_type == clip_wpthlp)) {
    max_mag_corr = max_mag_correlation_flux;
  } else {
    // All other covariances
    max_mag_corr = max_mag_correlation;
    // solve_type
  }

  // The value of x'y' at the surface (or lower boundary) is a set value that
  // is either specified or determined elsewhere in a surface subroutine.  It
  // is ensured elsewhere that the correlation between x and y at the surface
  // (or lower boundary) is between -1 and 1.  Thus, the covariance clipping
  // code does not need to be invoked at the lower boundary.  Likewise, the
  // value of x'y' is set at the upper boundary, so the covariance clipping
  // code does not need to be invoked at the upper boundary.
  // Note that if clipping were applied at the lower boundary, momentum will
  // not be conserved, therefore it should never be added.
  Kokkos::parallel_for(
      "clip_covar::loop1", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        real_t xpyp_bound;
        xpyp_bound = max_mag_corr * Kokkos::sqrt(xp2(i - 1, k - 1) * yp2(i - 1, k - 1));

        // Clipping for xpyp at an upper limit corresponding with a correlation
        // between x and y of max_mag_corr.
        if (xpyp(i - 1, k - 1) > xpyp_bound) {

          xpyp_chnge(i - 1, k - 1) = xpyp_bound - xpyp(i - 1, k - 1);
          xpyp(i - 1, k - 1)       = xpyp_bound;

          // Clipping for xpyp at a lower limit corresponding with a correlation
          // between x and y of -max_mag_corr.
        } else if (xpyp(i - 1, k - 1) < -xpyp_bound) {

          xpyp_chnge(i - 1, k - 1) = -xpyp_bound - xpyp(i - 1, k - 1);
          xpyp(i - 1, k - 1)       = -xpyp_bound;

        } else {

          xpyp_chnge(i - 1, k - 1) = real_t(0.0);
        }
      });

  // Since there is no covariance clipping at the upper or lower boundaries,
  // the change in x'y' due to covariance clipping at those levels is 0.
  Kokkos::parallel_for(
      "clip_covar::loop2", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        xpyp_chnge(i - 1, 1 - 1)   = real_t(0.0);
        xpyp_chnge(i - 1, nzm - 1) = real_t(0.0);
      });

  return;
}

//=============================================================================
// solve_type: Variable being solved; used for STATS
// dt: Model timestep; used here for STATS     [s]
// threshold_lo: Minimum value of x'^2                   [{x units}^2]
// xp2: Variance of x, x'^2 (momentum levels)   [{x units}^2]
// threshold_hi: Maximum value of x'^2                   [{x units}^2]
void clip_variance(int nzm, int ngrdcol, const grid& gr, int solve_type, real_t dt, // In
                   Real2D threshold_lo,                                             // In
                   f2k::stats_type& stats,                                          // InOut
                   Real2D xp2,                                                      // InOut
                   const real_t* threshold_hi)                                      // Optional In
{
  // Description:
  // Clipping the value of variance x'^2 based on a minimum and (optional) max threshold value.
  // The threshold value must be greater than or equal to 0.
  //
  // The values of x'^2 are found on the momentum levels.
  //
  // The following variances are found in the code:
  //
  // r_t'^2, th_l'^2, u'^2, v'^2, sclr'^2, (computed in advance_xp2_xpyp);
  // w'^2 (computed in advance_wp2_wp3).
  // References:
  // None
  //-----------------------------------------------------------------------
  // -------------------- Input Variables --------------------
  // -------------------- Output Variable --------------------
  // -------------------- Optional Input Variables --------------------
  // -------------------- Local Variables --------------------
  int i;
  int k; // Array index
  std::string name_xp2_cl;
  bool l_clubb_at_least_debug_level_3;
  const auto gr_k_lb_zm       = gr.k_lb_zm;
  const auto gr_k_ub_zm       = gr.k_ub_zm;
  const auto gr_grid_dir_indx = gr.grid_dir_indx;

  // -------------------- Begin Code --------------------

  // If a threshold_hi was passed, check if any element of threshold_lo
  // is larger than threshold_hi.
  // If so, terminate the program
  if ((threshold_hi != nullptr)) {
    // Copy the optional host scalar once before device use.
    const real_t threshold_hi_value = *threshold_hi;
    int invalid_thresholds          = 0;
    Kokkos::parallel_reduce(
        "clip_variance::invalid_thresholds", ngrdcol * std::max(0, nzm - 1),
        KOKKOS_LAMBDA(const int cell, int& hits) {
          const int i = cell % ngrdcol, k = cell / ngrdcol;
          if (threshold_lo(i, k) > threshold_hi_value)
            ++hits;
        },
        invalid_thresholds);
    if (invalid_thresholds != 0) {
      auto threshold_lo_host = to_host_space(threshold_lo);
      for (int k = 1; k <= nzm - 1; ++k) {
        for (int i = 1; i <= ngrdcol; ++i) {
          if (threshold_lo_host(i - 1, k - 1) > threshold_hi_value) {
            std::cerr << std::setprecision(std::numeric_limits<real_t>::max_digits10)
                      << "Error in clip_variance! threshold_lo(" << i << "," << k
                      << ") =" << threshold_lo_host(i - 1, k - 1) << ">" << threshold_hi_value
                      << "= threshold_hi. Terminating.." << "\n";
            // Preserve Native unqualified STOP, including normal process termination.
            std::exit(EXIT_SUCCESS);
          }
        }
      }
    }
  }

  name_xp2_cl = "";
  switch (solve_type) {
  case clip_wp2: {
    // wp2 clipping budget term
    name_xp2_cl = "wp2_cl";
    break;
  }
  case clip_rtp2: {
    // rtp2 clipping budget term
    name_xp2_cl = "rtp2_cl";
    break;
  }
  case clip_thlp2: {
    // thlp2 clipping budget term
    name_xp2_cl = "thlp2_cl";
    break;
  }
  case clip_up2: {
    // up2 clipping budget term
    name_xp2_cl = "up2_cl";
    break;
  }
  case clip_vp2: {
    // vp2 clipping budget term
    name_xp2_cl = "vp2_cl";
    break;
  }
  default: {
    // scalars are involved
    name_xp2_cl = "";
    break;
  }
  }
  if (stats.flags().l_sample) {
    {
      Real2D stats_values_3(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_variance::stats_values_3"),
          xp2.extent(0), xp2.extent(1));
      Kokkos::parallel_for(
          "clip_variance::stats_value", xp2.size(), KOKKOS_LAMBDA(int cell) {
            const int i          = cell % xp2.extent(0);
            const int k          = cell / xp2.extent(0);
            stats_values_3(i, k) = xp2(i, k) / dt;
          });
      f2k::stats_begin_budget(name_xp2_cl, stats_values_3, // In
                              stats);                      // InOut
    }
  }

  l_clubb_at_least_debug_level_3 = f2k::clubb_at_least_debug_level_api(3);

  // Limit the value of x'^2 at threshold.
  // The value of x'^2 at the surface (or lower boundary) is a set value that
  // is determined elsewhere in a surface subroutine.  Thus, the variance
  // clipping code does not need to be invoked at the lower boundary.
  // Likewise, the value of x'^2 is set at the upper boundary, so the variance
  // clipping code does not need to be invoked at the upper boundary.
  //
  // charlass on 09/11/2013: I changed the clipping so that also the surface
  // level is clipped. I did this because we discovered that there are slightly
  // negative values in thlp2(1) and rtp2(1) when running quarter_ss case with
  // WRF-CLUBB (see wrf:ticket:51#comment:33)
  if (l_clubb_at_least_debug_level_3) {
    auto xp2_host          = to_host_space(xp2);
    auto threshold_lo_host = to_host_space(threshold_lo);
    for (k = gr.k_lb_zm; (gr.grid_dir_indx > 0 ? k <= gr.k_ub_zm - gr.grid_dir_indx
                                               : k >= gr.k_ub_zm - gr.grid_dir_indx);
         k += gr.grid_dir_indx) {
      for (i = 1; i <= ngrdcol; i += 1) {
        if (xp2_host(i - 1, k - 1) < threshold_lo_host(i - 1, k - 1)) {
          if (l_clubb_at_least_debug_level_3) {
            {
              std::ostringstream diagnostic_message;
              diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
              diagnostic_message << "Warning: (solve_type==" << " " << solve_type << " "
                                 << ") xp2 =" << " " << xp2_host(i - 1, k - 1) << " " << "<" << " "
                                 << threshold_lo_host(i - 1, k - 1) << " " << " @ k =" << " " << k
                                 << " " << ". Small values are clipped." << "\n";
              std::cerr << diagnostic_message.str();
            }
          }
          xp2_host(i - 1, k - 1) = threshold_lo_host(i - 1, k - 1);
        }
      }
    }
    Kokkos::deep_copy(xp2, xp2_host);
  } else {
    Kokkos::parallel_for(
        "clip_variance::loop1",
        std::max(0, ((gr_k_ub_zm - gr_grid_dir_indx - (gr_k_lb_zm)) / (gr_grid_dir_indx) + 1)) *
            std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k =
              gr_k_lb_zm +
              remaining %
                  (((gr_k_ub_zm - gr_grid_dir_indx - (gr_k_lb_zm)) / (gr_grid_dir_indx) + 1)) *
                  (gr_grid_dir_indx);
          remaining /= (((gr_k_ub_zm - gr_grid_dir_indx - (gr_k_lb_zm)) / (gr_grid_dir_indx) + 1));
          if (xp2(i - 1, k - 1) < threshold_lo(i - 1, k - 1)) {
            xp2(i - 1, k - 1) = threshold_lo(i - 1, k - 1);
          }
        });
  }

  // Optional clipping of large values
  // This was added so the <var>_cl budget term also contains the upper clipping contributions
  if ((threshold_hi != nullptr)) {
    const real_t threshold_hi_value = *threshold_hi;
    if (l_clubb_at_least_debug_level_3) {
      auto xp2_host = to_host_space(xp2);
      for (k = 1; k <= nzm - 1; k += 1) {
        for (i = 1; i <= ngrdcol; i += 1) {
          if (xp2_host(i - 1, k - 1) > threshold_hi_value) {
            xp2_host(i - 1, k - 1) = threshold_hi_value;
            if (l_clubb_at_least_debug_level_3) {
              {
                std::ostringstream diagnostic_message;
                diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
                diagnostic_message << "Warning: (solve_type==" << " " << solve_type << " "
                                   << ") xp2 =" << " " << xp2_host(i - 1, k - 1) << " " << ">"
                                   << " " << threshold_hi_value << " " << " @ k =" << " " << k
                                   << " " << ". Large values are clipped." << "\n";
                std::cerr << diagnostic_message.str();
              }
            }
          }
        }
      }
      Kokkos::deep_copy(xp2, xp2_host);
    } else {
      Kokkos::parallel_for(
          "clip_variance::loop2", std::max(0, nzm - 1) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzm - 1);
            remaining /= (nzm - 1);
            if (xp2(i - 1, k - 1) > threshold_hi_value) {
              xp2(i - 1, k - 1) = threshold_hi_value;
            }
          });
    }
  }

  if (stats.flags().l_sample) {
    {
      Real2D stats_values_4(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_variance::stats_values_4"),
          xp2.extent(0), xp2.extent(1));
      Kokkos::parallel_for(
          "clip_variance::stats_value", xp2.size(), KOKKOS_LAMBDA(int cell) {
            const int i          = cell % xp2.extent(0);
            const int k          = cell / xp2.extent(0);
            stats_values_4(i, k) = xp2(i, k) / dt;
          });
      f2k::stats_finalize_budget(name_xp2_cl, stats_values_4, // In
                                 stats);                      // InOut
    }
  }

  return;
}

//=============================================================================
// dt: Model timestep; used here for STATS        [s]
// sfc_elevation: Elevation of ground level                  [m AMSL]
// Skw_max_mag: Maximum allowable magnitude of Skewness    [-]
// wp2_zt: w'^2 interpolated to thermodyamic levels   [m^2/s^2]
// wp3: w'^3 (thermodynamic levels)                [m^3/s^3]
void clip_skewness(int nzt, int ngrdcol, const grid& gr, real_t dt, Real1D sfc_elevation, // In
                   Real1D Skw_max_mag, Real2D wp2_zt,                                     // In
                   bool l_use_wp3_lim_with_smth_Heaviside,                                // In
                   f2k::stats_type& stats,                                                // InOut
                   Real2D wp3)                                                            // InOut
{
  // Description:
  // Clipping the value of w'^3 based on the skewness of w, Sk_w.
  //
  // Aditionally, to prevent possible crashes due to wp3 growing too large,
  // abs(wp3) will be clipped to 100.
  //
  // The skewness of w is:
  //
  // Sk_w = w'^3 / (w'^2)^(3/2).
  //
  // The value of Sk_w is limited to a range between an upper limit and a lower
  // limit.  The values of the limits depend on whether the level altitude is
  // within 100 meters of the surface.
  //
  // For altitudes less than or equal to 100 meters above ground level (AGL):
  //
  // -0.2_core_rknd*sqrt(2) <= Sk_w <= 0.2_core_rknd*sqrt(2);
  //
  // while for all altitudes greater than 100 meters AGL:
  //
  // -4.5_core_rknd <= Sk_w <= 4.5_core_rknd.
  //
  // Therefore, there is an upper limit on w'^3, such that:
  //
  // w'^3  <=  threshold_magnitude * (w'^2)^(3/2);
  //
  // and a lower limit on w'^3, such that:
  //
  // w'^3  >= -threshold_magnitude * (w'^2)^(3/2).
  //
  // The values of w'^3 are found on the thermodynamic levels, while the values
  // of w'^2 are found on the momentum levels.  Therefore, the values of w'^2
  // are interpolated to the thermodynamic levels before being used to
  // calculate the upper and lower limits for w'^3.
  // References:
  // None
  //-----------------------------------------------------------------------
  // ----------------------- Input Variables -----------------------
  // Flag to activate modifications on wp3 limiters for convergence test
  // (use smooth Heaviside 'Preskin' function in the calculation of
  // clip_skewness for wp3)
  // ----------------------- Input/Output Variables -----------------------

  // ----------------------- Local Variables -----------------------
  // ----------------------- Begin Code -----------------------

  if (stats.flags().l_sample) {
    {
      Real2D stats_values_1(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_skewness::stats_values_1"),
          wp3.extent(0), wp3.extent(1));
      Kokkos::parallel_for(
          "clip_skewness::stats_value", wp3.size(), KOKKOS_LAMBDA(int cell) {
            const int i          = cell % wp3.extent(0);
            const int k          = cell / wp3.extent(0);
            stats_values_1(i, k) = wp3(i, k) / dt;
          });
      f2k::stats_begin_budget("wp3_cl", stats_values_1, // In
                              stats);                   // InOut
    }
  }

  clip_skewness_core(nzt, ngrdcol, gr, sfc_elevation,   // In
                     Skw_max_mag, wp2_zt,               // In
                     l_use_wp3_lim_with_smth_Heaviside, // In
                     wp3);                              // InOut

  if (stats.flags().l_sample) {
    {
      Real2D stats_values_2(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_skewness::stats_values_2"),
          wp3.extent(0), wp3.extent(1));
      Kokkos::parallel_for(
          "clip_skewness::stats_value", wp3.size(), KOKKOS_LAMBDA(int cell) {
            const int i          = cell % wp3.extent(0);
            const int k          = cell / wp3.extent(0);
            stats_values_2(i, k) = wp3(i, k) / dt;
          });
      f2k::stats_finalize_budget("wp3_cl", stats_values_2, // In
                                 stats);                   // InOut
    }
  }

  return;
}

//=============================================================================
// sfc_elevation: Elevation of ground level                  [m AMSL]
// Skw_max_mag: Maximum allowable magnitude of Skewness    [-]
// wp2_zt: w'^2 interpolated to thermodyamic levels   [m^2/s^2]
// wp3: w'^3 (thermodynamic levels)                [m^3/s^3]
void clip_skewness_core(int nzt, int ngrdcol, const grid& gr, Real1D sfc_elevation, // In
                        Real1D Skw_max_mag, Real2D wp2_zt,                          // In
                        bool l_use_wp3_lim_with_smth_Heaviside,                     // In
                        Real2D wp3)                                                 // InOut
{
  // Description:
  // Apply the Native height-dependent skewness limit and absolute wp3 cap.
  //----------------------- Input Variables -----------------------
  // Flag to activate modifications on wp3 limiters for convergence test
  // (use smooth Heaviside 'Preskin' function in the calculation of clip_skewness for wp3)
  //----------------------- Input/Output Variables -----------------------
  //----------------------- Local Variables -----------------------
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wp2_zt_cubed(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_skewness_core::wp2_zt_cubed"), ngrdcol,
      nzt); // Variance of vertical velocity cubed (w^2_{zt}^3)   [m^6/s^6]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wp3_lim_sqd(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_skewness_core::wp3_lim_sqd"), ngrdcol,
      nzt); // Keeps absolute value of Sk_w from becoming > limit [m^6/s^6]
  int i;
  int k; // Vertical array index.
  constexpr real_t wp3_max = real_t(100.);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> zagl_thresh(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_skewness_core::zagl_thresh"), ngrdcol,
      nzt); // temporatory array
  Kokkos::View<real_t**, Kokkos::LayoutLeft> H_zagl(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "clip_skewness_core::H_zagl"), ngrdcol,
      nzt); // Heaviside function for clippings of wp3_lim_sqd
  const auto gr_zt = gr.zt;

  //----------------------- Begin Code-----------------------

  // Compute the upper and lower limits of w'^3 at every level,
  // based on the skewness of w, Sk_w, such that:
  // Sk_w = w'^3 / (w'^2)^(3/2);
  // -4.5 <= Sk_w <= 4.5;
  // or, if the level altitude is within 100 meters of the surface,
  // -0.2*sqrt(2) <= Sk_w <= 0.2*sqrt(2).

  // The normal magnitude limit of skewness of w in the CLUBB code is 4.5.
  // However, according to Andre et al. (1976b & 1978), wp3 should not exceed
  // [2*(wp2^3)]^(1/2) at any level.  However, this term should be multiplied
  // by 0.2 close to the surface to include surface effects.  There already is
  // a wp3 clipping term in place for all other altitudes, but this term will
  // be included for the surface layer only.  Therefore, the lowest level wp3
  // should not exceed 0.2 * sqrt(2) * wp2^(3/2).  Brian Griffin.  12/18/05.

  // To lower compute time, we squared both sides of the equation and compute
  // wp2^3 only once. -dschanen 9 Oct 2008
  Kokkos::parallel_for(
      "clip_skewness_core::loop1", std::max(0, nzt) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzt);
        remaining /= (nzt);
        wp2_zt_cubed(i - 1, k - 1) =
            ((wp2_zt(i - 1, k - 1) * wp2_zt(i - 1, k - 1)) * wp2_zt(i - 1, k - 1));
      });

  if (l_use_wp3_lim_with_smth_Heaviside) {

    // implement a smoothed Heaviside function to avoid discontinuities
    Kokkos::parallel_for(
        "clip_skewness_core::loop2", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzt);
          remaining /= (nzt);
          zagl_thresh(i - 1, k - 1) = (gr_zt(i - 1, k - 1) - sfc_elevation(i - 1)) / real_t(100.0);
          zagl_thresh(i - 1, k - 1) = zagl_thresh(i - 1, k - 1) - real_t(1.0);
        });

    Kokkos::deep_copy(
        Kokkos::subview(H_zagl, Kokkos::ALL, Kokkos::ALL),
        smooth_heaviside_peskin(
            nzt, ngrdcol, Kokkos::subview(zagl_thresh, Kokkos::ALL, Kokkos::ALL), real_t(0.6)));

    Kokkos::parallel_for(
        "clip_skewness_core::loop3", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzt);
          remaining /= (nzt);
          wp3_lim_sqd(i - 1, k - 1) =
              wp2_zt_cubed(i - 1, k - 1) *
              (H_zagl(i - 1, k - 1) * (Skw_max_mag(i - 1) * Skw_max_mag(i - 1)) +
               (real_t(1.0) - H_zagl(i - 1, k - 1)) * real_t(0.0021) *
                   (Skw_max_mag(i - 1) * Skw_max_mag(i - 1)));
        });

  } else {
    // default method

    Kokkos::parallel_for(
        "clip_skewness_core::loop4", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzt);
          remaining /= (nzt);
          // Clip for 100 m. AGL.
          if (gr_zt(i - 1, k - 1) - sfc_elevation(i - 1) <= real_t(100.0)) {
            // wp3_upper_lim(k) =  0.2_core_rknd * sqrt_2 * wp2_zt(k)**(3.0_core_rknd/2.0_core_rknd)
            // wp3_lower_lim(k) = -0.2_core_rknd * sqrt_2 * wp2_zt(k)**(3.0_core_rknd/2.0_core_rknd)
            wp3_lim_sqd(i - 1, k - 1) = real_t(0.0021) * (Skw_max_mag(i - 1) * Skw_max_mag(i - 1)) *
                                        wp2_zt_cubed(i - 1, k - 1);
          } else {
            // Clip skewness consistently with a.
            // wp3_upper_lim(k) =  4.5_core_rknd * wp2_zt(k)**(3.0_core_rknd/2.0_core_rknd)
            // wp3_lower_lim(k) = -4.5_core_rknd * wp2_zt(k)**(3.0_core_rknd/2.0_core_rknd)
            // Skw_max_mag = 4.5_core_rknd^2
            wp3_lim_sqd(i - 1, k - 1) =
                (Skw_max_mag(i - 1) * Skw_max_mag(i - 1)) * wp2_zt_cubed(i - 1, k - 1);
          }
        });
  }

  // Clipping for w'^3 at an upper and lower limit corresponding with
  // the appropriate value of Sk_w.
  Kokkos::parallel_for(
      "clip_skewness_core::loop5", std::max(0, nzt) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzt);
        remaining /= (nzt);
        // Set the magnitude to the wp3 limit and apply the sign of the current wp3
        if ((wp3(i - 1, k - 1) * wp3(i - 1, k - 1)) > wp3_lim_sqd(i - 1, k - 1)) {
          wp3(i - 1, k - 1) =
              Kokkos::copysign(Kokkos::sqrt(wp3_lim_sqd(i - 1, k - 1)), wp3(i - 1, k - 1));
        }
      });

  // Clipping abs(wp3) to 100. This keeps wp3 from growing too large in some
  // deep convective cases, which helps prevent these cases from blowing up.
  Kokkos::parallel_for(
      "clip_skewness_core::loop6", std::max(0, nzt) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzt);
        remaining /= (nzt);
        if (Kokkos::abs(wp3(i - 1, k - 1)) > wp3_max) {
          // Known magic number
          wp3(i - 1, k - 1) = Kokkos::copysign(wp3_max, wp3(i - 1, k - 1));
        }
      });
}

//=============================================================================
// rtm: Total water mixing ratio             [kg/kg]
// rcm: Cloud water mixing ratio  [kg/kg]
void clip_rcm(int nzt, int ngrdcol, Real2D rtm, // In
              const std::string& message,       // In
              Real2D rcm)                       // InOut
{
  //
  // Description:
  //   Subroutine that reduces cloud water (rcm) whenever
  //   it exceeds total water (rtm = vapor + liquid).
  //   This avoids negative values of rvm = water vapor mixing ratio.
  //   However, it will not ensure that rcm <= rtm if rtm <= 0.
  //
  // References:
  //   None
  //---------------------------------------------------------------------
  // Input variables
  int i;
  int k;

  // ------------ Begin code ---------------

  if (f2k::clubb_at_least_debug_level_api(3)) {
    auto rtm_host = to_host_space(rtm);
    auto rcm_host = to_host_space(rcm);

    for (k = 1; k <= nzt; k += 1) {
      for (i = 1; i <= ngrdcol; i += 1) {

        if (rtm_host(i - 1, k - 1) < rcm_host(i - 1, k - 1)) {

          {
            std::ostringstream diagnostic_message;
            diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            diagnostic_message << message << " " << " at k=" << " " << k << " " << " at i=" << " "
                               << i << " " << "rcm(k) = " << " " << rcm_host(i - 1, k - 1) << " "
                               << "rtm(k) = " << " " << rtm_host(i - 1, k - 1) << " " << "." << " "
                               << " Clipping rcm." << "\n";
            std::cerr << diagnostic_message.str();
          }

          // rtm(k) < rcm(k)
        }
      }
    }
    // clubb_at_least_debug_level_api( 3 )
  }

  // Vince Larson clipped rcm in order to prevent rvm < 0.  5 Apr 2008.
  // This code won't work unless rtm >= 0 !!!
  // We do not clip rcm_in_layer because rcm_in_layer only influences
  // radiation, and we do not want to bother recomputing it.  6 Aug 2009
  Kokkos::parallel_for(
      "clip_rcm::loop1", std::max(0, nzt) * std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzt);
        remaining /= (nzt);
        if (rtm(i - 1, k - 1) < rcm(i - 1, k - 1)) {
          rcm(i - 1, k - 1) =
              ((zero_threshold) > (rtm(i - 1, k - 1) - std::numeric_limits<real_t>::epsilon())
                   ? (zero_threshold)
                   : (rtm(i - 1, k - 1) - std::numeric_limits<real_t>::epsilon()));
          // rtm(k) < rcm(k)
        }
      });

  return;
}

} // namespace clubb_kokkos
