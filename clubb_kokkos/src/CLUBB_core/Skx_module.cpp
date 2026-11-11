// Skewness and third-moment closures; see Skx_module.hpp for public routines.
// Cell kernels retain the Native formulas; third-moment scratch reuses case storage.
#include "scratch_workspace.hpp"
#include "Skx_module.hpp"
#include "constants_clubb.hpp"
#include "model_flags.hpp"
#include "parameter_indices.hpp"

namespace clubb_kokkos {

void Skx_func(int nz, int ngrdcol, Real2D xp2, Real2D xp3, // In
              real_t x_tol, Real2D clubb_params,           // In
              Real2D Skx) {                                // Out
  // Calculate the skewness of x from its second and third moments. x_tol has
  // x's units; the output is dimensionless.

  // Parameter Constants: whether to apply clipping to the final result.
  constexpr bool l_clipping_kluge = false;

  // ---- Begin Code ----

  // Skx = xp3 / ( max( xp2, x_tol**two ) )**three_halves
  // Calculation of skewness to help reduce the sensitivity of this value to
  // small values of xp2.
  Kokkos::parallel_for(
      "Skx_func::skewness", Kokkos::RangePolicy<>(0, ngrdcol * nz), KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        // Native integer power **(-3) squares, multiplies and then reciprocates.
        // A real-exponent pow changes rounding even at -O0.
        const real_t Skx_denom_tol =
            clubb_params(i, iSkw_denom_coef) * (x_tol * x_tol); // Variance tolerance
        const real_t sqrt_xp2 = Kokkos::sqrt(xp2(i, k) + Skx_denom_tol);
        Skx(i, k)             = xp3(i, k) * (real_t(1) / (sqrt_xp2 * (sqrt_xp2 * sqrt_xp2)));
      });

  // This is no longer needed since clipping is already
  // imposed on wp2 and wp3 elsewhere in the code.
  // I turned clipping on in this local copy since thlp3 and rtp3 are not clipped.
  if (l_clipping_kluge) {
    Kokkos::parallel_for(
        "Skx_func::clip", Kokkos::RangePolicy<>(0, ngrdcol * nz), KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = cell / ngrdcol;
          Skx(i, k)   = Kokkos::min(Kokkos::max(Skx(i, k), -clubb_params(i, iSkw_max_mag)),
                                    clubb_params(i, iSkw_max_mag));
        });
  }
}

void compute_gamma_Skw(int nz, int ngrdcol, Real2D Skw, Real2D clubb_params, // In
                       Real2D gamma_Skw_fnc) {                               // Out
  // Compute gamma as a function of Skw on a single vertical grid.
  // Both skewness and the returned PDF width coefficient are dimensionless.

  // ----------------------------- Begin Code ------------------------------
  if (l_gamma_Skw) {
    Kokkos::parallel_for(
        "compute_gamma_Skw::varying", Kokkos::RangePolicy<>(0, ngrdcol * nz),
        KOKKOS_LAMBDA(const int cell) {
          const int i              = cell % ngrdcol;
          const int k              = cell / ngrdcol;
          const real_t gamma_coef  = clubb_params(i, igamma_coef);
          const real_t gamma_coefb = clubb_params(i, igamma_coefb);
          const real_t gamma_coefc = clubb_params(i, igamma_coefc);

          if (Kokkos::abs(gamma_coef - gamma_coefb) >
              Kokkos::abs(gamma_coef + gamma_coefb) * eps / 2) {
            gamma_Skw_fnc(i, k) =
                gamma_coefb +
                (gamma_coef - gamma_coefb) * Kokkos::exp(-one_half * ((Skw(i, k) / gamma_coefc) *
                                                                      (Skw(i, k) / gamma_coefc)));
          } else {
            gamma_Skw_fnc(i, k) = gamma_coef;
          }
        });
  } else {
    Kokkos::parallel_for(
        "compute_gamma_Skw::constant", Kokkos::RangePolicy<>(0, ngrdcol * nz),
        KOKKOS_LAMBDA(const int cell) {
          const int i         = cell % ngrdcol;
          const int k         = cell / ngrdcol;
          gamma_Skw_fnc(i, k) = clubb_params(i, igamma_coef);
        });
  }
}

void LG_2005_ansatz(int nz, int ngrdcol, Real2D Skw, Real2D wpxp, Real2D wp2,  // In
                    Real2D xp2, Real1D beta, Real2D sigma_sqd_w, real_t x_tol, // In
                    Real2D Skx) {                                              // Out
  // Calculate skewness of x using Larson and Golaz (2005), Mon. Wea. Rev.,
  // 133, 1023–1042. wpxp is the turbulent flux [m/s (x units)], wp2 and xp2
  // are variances, beta and sigma_sqd_w are dimensionless, x_tol has x's units.

  // -------------------------- Begin Code --------------------------
  // weberjk, 8-July 2015. Commented this out for now. cgils was failing during some tests.

  // Larson and Golaz (2005) eq. 16
  Kokkos::parallel_for(
      "LG_2005_ansatz", Kokkos::RangePolicy<>(0, ngrdcol * nz), KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        const real_t nrmlzd_corr_wx =
            wpxp(i, k) /
            Kokkos::sqrt(Kokkos::max(wp2(i, k), real_t(w_tol_sqd)) *
                         Kokkos::max(xp2(i, k), x_tol * x_tol) * (one - sigma_sqd_w(i, k)));

        // Larson and Golaz (2005) eq. 11
        const real_t nrmlzd_Skw =
            Skw(i, k) / ((one - sigma_sqd_w(i, k)) * Kokkos::sqrt(one - sigma_sqd_w(i, k)));

        // Larson and Golaz (2005) eq. 33
        Skx(i, k) = nrmlzd_Skw * nrmlzd_corr_wx *
                    (beta(i) + (one - beta(i)) * (nrmlzd_corr_wx * nrmlzd_corr_wx));
      });
}

void xp3_LG_2005_ansatz(ScratchWorkspace& workspace,                                        // InOut
                        int nzt, int ngrdcol, Real2D Skw_zt, Real2D wpxp_zt, Real2D wp2_zt, // In
                        Real2D xp2_zt, Real2D sigma_sqd_w_zt,                               // In
                        Real2D clubb_params, real_t x_tol,                                  // In
                        Real2D xp3) {                                                       // Out
  // Calculate <x'^3> [(x units)^3] after diagnosing skewness using the ansatz
  // of Larson and Golaz (2005). Input moments are on thermodynamic levels.

  // -------------------------- Local Variable --------------------------
  auto workspace_scope = workspace.scope();
  Real2D Skx_zt        = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "xp3_LG_2005_ansatz::Skx_zt", ngrdcol, nzt); // Skewness of x [-]

  // -------------------------- Begin Code --------------------------
  // Calculate skewness of x using the ansatz of LG05.
  LG_2005_ansatz(nzt, ngrdcol, Skw_zt, wpxp_zt, wp2_zt, // In
                 xp2_zt, Kokkos::subview(clubb_params, Kokkos::ALL, ibeta), sigma_sqd_w_zt,
                 x_tol,   // In
                 Skx_zt); // Out

  // Calculate <x'^3> using the reverse of the special sensitivity reduction
  // formula in function Skx_func above.
  Kokkos::parallel_for(
      "xp3_LG_2005_ansatz::third_moment", Kokkos::RangePolicy<>(0, ngrdcol * nzt),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        const real_t Skx_denom_tol =
            clubb_params(i, iSkw_denom_coef) * (x_tol * x_tol); // Variance tolerance
        xp3(i, k) = Skx_zt(i, k) * (xp2_zt(i, k) + Skx_denom_tol) *
                    Kokkos::sqrt(xp2_zt(i, k) + Skx_denom_tol);
      });
}

} // namespace clubb_kokkos
