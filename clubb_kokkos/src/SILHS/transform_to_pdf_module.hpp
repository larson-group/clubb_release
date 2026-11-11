// Convert SILHS uniform variates to the CLUBB normal/lognormal PDF.
// Public functions: transform_uniform_samples_to_pdf, ltqnorm,
// multiply_Cholesky and chi_eta_2_rtthl. Scratch Cholesky/standard-normal arrays
// retain Native (variate,column,level,...) order; physical samples are column first.
#pragma once
#include "scratch_workspace.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/corr_varnce_module.hpp"
#include "CLUBB_core/constants_clubb.hpp"
namespace clubb_kokkos {
void transform_uniform_samples_to_pdf(
    ScratchWorkspace& workspace,                                             // InOut
    int nzt, int ngrdcol, int num_samples, int pdf_dim, int d_uniform_extra, // In
    const hm_metadata_type& hm_metadata,                                     // In
    Real4D Sigma_Cholesky1, Real4D Sigma_Cholesky2,                          // In
    Real3D mu1, Real3D mu2,                                                  // In
    Kokkos::View<int***, Kokkos::LayoutLeft> X_mixt_comp_all_levs,           // In
    Real4D X_u_all_levs, Real3D cloud_frac,                                  // In
    Kokkos::View<std::uint8_t***, Kokkos::LayoutLeft> l_in_precip_all_levs,  // In
    Real4D X_nl_all_levs                                                     // Out
);                                                                           // Out

KOKKOS_INLINE_FUNCTION real_t ltqnorm(real_t p_core_rknd // In
) {
  // Peter J. Acklam (2003-04-23), http://home.online.no/~pjacklam.
  // Hongli Jiang ported the original MATLAB quantile (17Feb2004).
  // Acklam's lower-tail standard-normal quantile, converted to double by
  // V. Larson (22Feb2004), with the Native Halley refinement (23Aug2014).
  // The rational minimax approximation has relative error below 1.15e-9;
  // Halley's rational method improves this to full machine precision.
  // Native invalid/end-point branches are retained, including subsequent
  // Halley evaluation; do not silently substitute JAX's guarded endpoint result.
  constexpr bool l_apply_halley_method = true;
  constexpr double a1 = -3.969683028665376E+01, a2 = 2.209460984245205E+02,
                   a3 = -2.759285104469687E+02, a4 = 1.383577518672690E+02,
                   a5 = -3.066479806614716E+01, a6 = 2.506628277459239E+00;
  constexpr double b1 = -5.447609879822406E+01, b2 = 1.615858368580409E+02,
                   b3 = -1.556989798598866E+02, b4 = 6.680131188771972E+01,
                   b5 = -1.328068155288572E+01;
  constexpr double c1 = -7.784894002430293E-03, c2 = -3.223964580411365E-01,
                   c3 = -2.400758277161838E+00, c4 = -2.549732539343734E+00,
                   c5 = 4.374664141464968E+00, c6 = 2.938163982698783E+00;
  constexpr double d1 = 7.784695709041462E-03, d2 = 3.224671290700398E-01,
                   d3 = 2.445134137142996E+00, d4 = 3.754408661907416E+00;
  const double p     = double(p_core_rknd);
  double z           = 0.0;
  const double plow  = .02425;
  const double phigh = 1.0 - plow;
  if (p > 0.0 && p < plow) {
    const double q = Kokkos::sqrt(-2.0 * Kokkos::log(p));
    z              = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
                     ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
  } else if (p >= plow && p <= phigh) {
    const double q = p - .5;
    const double r = q * q;
    z              = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
                     (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0);
  } else if (p > phigh && p < 1.0) {
    const double q = Kokkos::sqrt(-2.0 * Kokkos::log(1.0 - p));
    z              = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
                     ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
  }
  // Ordinary standalone Native uses arithmetic infinities rather than CAM's
  // external constants; preserve its comparisons and operation sequence.
  const double z1 = 0.0;
  if (Kokkos::abs(p) < eps)
    z = (-1.0) / z1;
  if (Kokkos::abs(p - 1.0) < Kokkos::abs(p + 1.0) / 2 * eps)
    z = 1.0 / z1;
  if (p < 0.0 || p > 1.0)
    z = (1.0 / z1) * (1.0 / z1);
  if (l_apply_halley_method) {
    const double e = one_half_dp * Kokkos::erfc(-z / sqrt_2_dp) - p;
    const double u = e * sqrt_2pi_dp * Kokkos::exp((z * z) / two_dp);
    z              = z - u / (one_dp + z * u / two_dp);
  }
  return real_t(z);
}

void multiply_Cholesky(int nzt, int ngrdcol, int num_samples, int pdf_dim, Real4D std_normal, // In
                       Real4D Sigma_Cholesky1, Real4D Sigma_Cholesky2,                        // In
                       Real3D mu1, Real3D mu2,                                                // In
                       Kokkos::View<int***, Kokkos::LayoutLeft> X_mixt_comp_all_levs,         // In
                       Real4D X_nl_all_levs                                                   // Out
);                                                                                            // Out
void chi_eta_2_rtthl(int nzt, int ngrdcol, int num_samples,                                   // In
                     Real2D rt_1, Real2D thl_1,                                               // In
                     Real2D rt_2, Real2D thl_2,                                               // In
                     Real2D crt_1, Real2D cthl_1,                                             // In
                     Real2D crt_2, Real2D cthl_2,                                             // In
                     Real2D mu_chi_1, Real2D mu_chi_2,                                        // In
                     Real3D chi, Real3D eta,                                                  // In
                     Kokkos::View<int***, Kokkos::LayoutLeft> X_mixt_comp_all_levs,           // In
                     Real3D lh_rt, Real3D lh_thl                                              // Out
);                                                                                            // Out
} // namespace clubb_kokkos
