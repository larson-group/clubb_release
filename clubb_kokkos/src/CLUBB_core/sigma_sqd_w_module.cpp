// PDF width calculation; see sigma_sqd_w_module.hpp for its public interface.
#include "scratch_workspace.hpp"
#include "sigma_sqd_w_module.hpp"
#include "Skx_module.hpp"
#include "constants_clubb.hpp"

namespace clubb_kokkos {

void compute_sigma_sqd_w(ScratchWorkspace& workspace,                       // InOut
                         int nzm, int nzt, int ngrdcol, const grid& gr,     // In
                         Real2D wp3, Real2D wp2, Real2D thlp2, Real2D rtp2, // In
                         Real2D up2, Real2D vp2, Real2D wpthlp, Real2D wprtp, Real2D upwp,
                         Real2D vpwp,              // In
                         Real2D clubb_params,      // In
                         bool l_predict_upwp_vpwp, // In
                         Real2D sigma_sqd_w) {     // Out
  // Compute the variable sigma_sqd_w (PDF width parameter).
  //
  // The value of sigma_sqd_w is restricted in the ADG1 PDF in order to keep
  // the marginal PDFs of all responder variables (variables that do not set
  // the mixture fraction) valid. The limits on sigma_sqd_w in order to keep
  // the PDF of a responder variable, x, valid are:
  //
  // 0 <= sigma_sqd_w <= 1 - <w'x'>^2 / ( <w'^2> * <x'^2> ).
  //
  // The overall limits on sigma_sqd_w must be applied based on the most
  // restrictive case so that all Double Gaussian PDF responder variables, x,
  // have realizable PDFs. The overall limits on sigma_sqd_w are:
  //
  // 0 <= sigma_sqd_w <= 1 - max( <w'x'>^2 / ( <w'^2> * <x'^2> ), for all x).
  //
  // The equation used for sigma_sqd_w is:
  //
  // sigma_sqd_w = gamma_Skw_fnc
  //               * ( 1 - max( <w'x'>^2 / ( <w'^2> * <x'^2> ), for all x) );
  //
  // where 0 <= gamma_Skw_fnc <= 1.
  // References: Eqn 22 in Equations for CLUBB.

  // Local Variable
  auto workspace_scope = workspace.scope();
  Real2D wp3_zm =
      workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "compute_sigma_sqd_w::wp3_zm",
                                   ngrdcol, nzm); // Third moment of w [m^3/s^3]
  Real2D Skw_zm =
      workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "compute_sigma_sqd_w::Skw_zm",
                                   ngrdcol, nzm); // Skewness of w [-]
  Real2D gamma_Skw_fnc = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "compute_sigma_sqd_w::gamma_Skw_fnc", ngrdcol, nzm); // Gamma [-]
  Real2D sigma_sqd_w_tmp = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                        "compute_sigma_sqd_w::sigma_sqd_w_tmp",
                                                        ngrdcol, nzm); // Unsmoothed width [-]

  // ---- Begin Code ----
  if (nzm > 1) {
    zt2zm_api(nzm, nzt, ngrdcol, gr, wp3, wp3_zm);
  } else {
    // Skip interpolation if nzm == 1, this only occurs for testing purposes.
    Kokkos::deep_copy(wp3_zm, wp3);
  }

  Skx_func(nzm, ngrdcol, wp2, wp3_zm, // In
           w_tol, clubb_params,       // In
           Skw_zm);                   // Out

  compute_gamma_Skw(nzm, ngrdcol, Skw_zm, clubb_params, // In
                    gamma_Skw_fnc);                     // Out

  //----------------------------------------------------------------
  // Compute sigma_sqd_w with new formula from Vince
  //----------------------------------------------------------------
  // Find the maximum value of <w'x'>^2 / ( <w'^2> * <x'^2> ) for all
  // variables x that are Double Gaussian PDF responder variables. This
  // includes rt and theta-l. When l_predict_upwp_vpwp is enabled, u and v are
  // also calculated as part of the PDF, and they are included as well.
  // Additionally, when sclr_dim > 0, passive scalars (sclr) are also included.
  Kokkos::parallel_for(
      "compute_sigma_sqd_w::width", Kokkos::RangePolicy<>(0, ngrdcol * nzm),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;

        // Named factors retain the native integer-square operation without pow.
        const real_t corr_w_thl =
            wpthlp(i, k) / (Kokkos::sqrt(wp2(i, k) * thlp2(i, k)) + one_hundred * w_tol * thl_tol);
        const real_t corr_w_rt =
            wprtp(i, k) / (Kokkos::sqrt(wp2(i, k) * rtp2(i, k)) + one_hundred * w_tol * rt_tol);
        real_t max_corr_w_x_sqd = Kokkos::max(corr_w_thl * corr_w_thl, corr_w_rt * corr_w_rt);

        if (l_predict_upwp_vpwp) {
          const real_t corr_w_u =
              upwp(i, k) / (Kokkos::sqrt(up2(i, k) * wp2(i, k)) + one_hundred * w_tol_sqd);
          const real_t corr_w_v =
              vpwp(i, k) / (Kokkos::sqrt(vp2(i, k) * wp2(i, k)) + one_hundred * w_tol_sqd);
          max_corr_w_x_sqd =
              Kokkos::max(Kokkos::max(max_corr_w_x_sqd, corr_w_u * corr_w_u), corr_w_v * corr_w_v);
        }

        sigma_sqd_w_tmp(i, k) =
            gamma_Skw_fnc(i, k) * (one - Kokkos::min(max_corr_w_x_sqd, real_t(one)));
      });

  if (nzm > 1) {
    // Smooth in the vertical using interpolation.
    zm2zt2zm(workspace, nzm, nzt, ngrdcol, gr, sigma_sqd_w_tmp, sigma_sqd_w, zero_threshold);
  } else {
    // Skip interpolation if nzm == 1, this only occurs for testing purposes.
    Kokkos::deep_copy(sigma_sqd_w, sigma_sqd_w_tmp);
  }
}

} // namespace clubb_kokkos
