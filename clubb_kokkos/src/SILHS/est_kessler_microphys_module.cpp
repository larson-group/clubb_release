// SILHS and analytic gridbox Kessler autoconversion estimates.
// Public function: est_kessler_microphys_api; private helper: calc_estimate.
#include "scratch_workspace.hpp"
#include "SILHS/est_kessler_microphys_module.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace clubb_kokkos {
namespace {
// The private estimate is declared here and defined in Native source order
// after its public parent. Host validation retains source-ordered diagnostics;
// the device helper retains Native guards for the reached calculation.
template <class Rc, class Component, class Weight>
KOKKOS_INLINE_FUNCTION void calc_estimate(int num_samples, real_t mixt_frac,               // In
                                          real_t cloud_frac_1, real_t cloud_frac_2, Rc rc, // In
                                          Component X_mixt_comp_one_lev,
                                          Weight lh_sample_point_weights, // In
                                          bool l_lh_importance_sampling,  // In
                                          real_t coeff, real_t r_crit,    // In
                                          real_t& est_m                   // Out
);                                                                        // Out
} // namespace

//-----------------------------------------------------------------------
void est_kessler_microphys_api(ScratchWorkspace& workspace,                        // InOut
                               int nzt, int num_samples, int pdf_dim, int ngrdcol, // In
                               Real4D X_nl_all_levs, const pdf_parameter& pdf_params, Real2D rcm,
                               Real2D cloud_frac,                                             // In
                               Kokkos::View<int***, Kokkos::LayoutLeft> X_mixt_comp_all_levs, // In
                               Real3D lh_sample_point_weights,                                // In
                               bool l_lh_importance_sampling,                                 // In
                               Real2D lh_AKm, Real2D AKm, Real2D AKstd, Real2D AKstd_cld,     // Out
                               Real2D AKm_rcm, Real2D AKm_rcc, Real2D lh_rcm_avg              // Out
) {
  // Compute Kessler gridbox means and standard deviations by Latin-hypercube
  // sampling and analytic integration. Liquid mixing ratios are kg/kg;
  // autoconversion means and deviations are kg/kg/s.
  // nzt: Number of vertical levels
  // num_samples: Number of sample points
  // pdf_dim: Number of variates
  // ngrdcol: Number of model columns
  // X_nl_all_levs: Sample that is transformed ultimately to normal-lognormal
  // pdf_params: PDF parameters       [units vary]
  // rcm: Liquid water mixing ratio                [kg/kg]
  // cloud_frac: Cloud fraction           [-]
  // X_mixt_comp_all_levs: Whether we're in mixture component 1 or 2
  // lh_sample_point_weights: Weight for cloud weighted sampling
  // l_lh_importance_sampling: Do importance sampling (SILHS) [-]
  // lh_AKm: Monte Carlo estimate of Kessler autoconversion [kg/kg/s]
  // AKm: Exact Kessler autoconversion, AKm,             [kg/kg/s]
  // AKstd: Exact standard deviation of gba Kessler        [kg/kg/s]
  // AKstd_cld: Exact w/in cloud std of gba Kessler            [kg/kg/s]
  // AKm_rcm: Exact local gba Kessler auto based on rcm      [kg/kg/s]
  // AKm_rcc: Exact local gba Kessler based on w/in cloud rc [kg/kg/s]
  // lh_rcm_avg: LH estimate of grid box avg liquid water [kg/kg]
  // Kokkos::max takes references; capture addressable values for device code.
  auto workspace_scope        = workspace.scope();
  const real_t zero           = clubb_kokkos::zero;
  const real_t chi_tol        = clubb_kokkos::chi_tol;
  const real_t zero_threshold = clubb_kokkos::zero_threshold;
  const auto pdf_mixt_frac    = pdf_params.mixt_frac;
  const auto pdf_chi_1        = pdf_params.chi_1;
  const auto pdf_chi_2        = pdf_params.chi_2;
  const auto pdf_stdev_chi_1  = pdf_params.stdev_chi_1;
  const auto pdf_stdev_chi_2  = pdf_params.stdev_chi_2;

  // Native calc_estimate rejects invalid fractions and an empty sample at the
  // first reached (level,column). Component cloud fractions are forced to one
  // below. Keep error I/O ordered on the host, before launching device science.
  const auto mixt_frac_host = to_host_space(workspace, pdf_mixt_frac);
  for (int level = 0; level < nzt; ++level) {
    for (int i = 0; i < ngrdcol; ++i) {
      if (mixt_frac_host(i, level) > one || mixt_frac_host(i, level) < zero) {
        std::cerr
            << "Error in calc_estimate:  mixture fraction, mixt_frac, does not lie in [0,1].\n"
            << "mixt_frac = " << std::setprecision(std::numeric_limits<real_t>::max_digits10)
            << mixt_frac_host(i, level) << '\n';
        throw std::runtime_error("Error stop in calc_estimate");
      }
      if (num_samples == 0)
        throw std::runtime_error("Error:  no sample points in calc_estimate");
    }
  }
  Real3D rcm_sample = workspace_scope.view<Real3D>(Kokkos::WithoutInitializing, "rcm_sample",
                                                   ngrdcol, num_samples, nzt);
  Kokkos::parallel_for(
      "est_kessler_microphys::samples", ngrdcol * num_samples * nzt, KOKKOS_LAMBDA(int cell) {
        const int i                  = cell % ngrdcol;
        const int sample             = (cell / ngrdcol) % num_samples;
        const int level              = cell / (ngrdcol * num_samples);
        rcm_sample(i, sample, level) = Kokkos::max(X_nl_all_levs(i, sample, level, 0), zero);
      });
  Kokkos::parallel_for(
      "est_kessler_microphys_api", ngrdcol * nzt, KOKKOS_LAMBDA(int cell) {
        const int i     = cell % ngrdcol;
        const int level = cell / ngrdcol;
        // Extract PDF parameters. Unit component cloud fractions implement Native
        // in-and-out-of-cloud sampling (D. Schanen, 30 July 2009).
        const real_t mixt_frac    = pdf_mixt_frac(i, level);
        const real_t cloud_frac_1 = real_t(1);
        const real_t cloud_frac_2 = real_t(1);
        const real_t chi_1        = pdf_chi_1(i, level);
        const real_t chi_2        = pdf_chi_2(i, level);
        const real_t stdev_chi_1  = pdf_stdev_chi_1(i, level);
        const real_t stdev_chi_2  = pdf_stdev_chi_2(i, level);

        // Call Kessler autoconversion with the Latin-hypercube sample, then compute
        // its analytic result. A_K=coeff*(rc-r_crit)*H(rc-r_crit).
        real_t coeff         = real_t(1.e-3);
        real_t r_crit        = real_t(.2) / g_per_kg;
        const auto rc        = Kokkos::subview(rcm_sample, i, Kokkos::ALL, level);
        const auto component = Kokkos::subview(X_mixt_comp_all_levs, i, Kokkos::ALL, level);
        const auto weights   = Kokkos::subview(lh_sample_point_weights, i, Kokkos::ALL, level);
        calc_estimate(num_samples, mixt_frac,         // In
                      cloud_frac_1, cloud_frac_2, rc, // In
                      component, weights,             // In
                      l_lh_importance_sampling,       // In
                      coeff, r_crit,                  // In
                      lh_AKm(i, level)                // Out
        );                                            // Out
        // Monte Carlo liquid estimate for comparison with the analytic mean.
        coeff  = one;
        r_crit = zero;
        calc_estimate(num_samples, mixt_frac,         // In
                      cloud_frac_1, cloud_frac_2, rc, // In
                      component, weights,             // In
                      l_lh_importance_sampling,       // In
                      coeff, r_crit,                  // In
                      lh_rcm_avg(i, level)            // Out
        );                                            // Out

        // Exact Kessler autoconversion [kg/kg/s].
        r_crit                    = real_t(.2e-3);
        const real_t K_one        = real_t(1.e-3);
        const real_t chi_n_1_crit = (chi_1 - r_crit) / Kokkos::max(stdev_chi_1, chi_tol);
        const real_t cloud_frac_1_crit =
            real_t(.5) * (real_t(1) + Kokkos::erf(chi_n_1_crit / Kokkos::sqrt(real_t(2))));
        const real_t AK1 =
            K_one * ((chi_1 - r_crit) * cloud_frac_1_crit +
                     stdev_chi_1 * Kokkos::exp(-real_t(.5) * (chi_n_1_crit * chi_n_1_crit)) /
                         Kokkos::sqrt(real_t(2) * pi));
        const real_t chi_n_2_crit = (chi_2 - r_crit) / Kokkos::max(stdev_chi_2, chi_tol);
        const real_t cloud_frac_2_crit =
            real_t(.5) * (real_t(1) + Kokkos::erf(chi_n_2_crit / Kokkos::sqrt(real_t(2))));
        const real_t AK2 =
            K_one * ((chi_2 - r_crit) * cloud_frac_2_crit +
                     stdev_chi_2 * Kokkos::exp(-real_t(.5) * (chi_n_2_crit * chi_n_2_crit)) /
                         Kokkos::sqrt(real_t(2) * pi));
        AKm(i, level) = mixt_frac * AK1 + (real_t(1) - mixt_frac) * AK2;

        // Exact standard deviations [kg/kg/s]. Native floors component variances
        // at zero_threshold because roundoff can otherwise make them negative.
        const real_t AK1var = Kokkos::max(
            zero_threshold, K_one * (chi_1 - r_crit) * AK1 +
                                K_one * K_one * (stdev_chi_1 * stdev_chi_1) * cloud_frac_1_crit -
                                (AK1 * AK1));
        const real_t AK2var = Kokkos::max(
            zero_threshold, K_one * (chi_2 - r_crit) * AK2 +
                                K_one * K_one * (stdev_chi_2 * stdev_chi_2) * cloud_frac_2_crit -
                                (AK2 * AK2));
        // Gridbox average.
        AKstd(i, level) = Kokkos::sqrt(
            mixt_frac * ((AK1 - AKm(i, level)) * (AK1 - AKm(i, level)) + AK1var) +
            (real_t(1) - mixt_frac) * ((AK2 - AKm(i, level)) * (AK2 - AKm(i, level)) + AK2var));
        // Within-cloud average; keep Native conditional evaluation, not extra guards.
        if (cloud_frac(i, level) > real_t(0)) {
          // BFB-029: square the rounded Native quotient; duplicate divisions reassociate.
          const real_t AKm_over_cloud_frac = AKm(i, level) / cloud_frac(i, level);
          AKstd_cld(i, level)              = Kokkos::sqrt(Kokkos::max(
              real_t(zero_threshold), (real_t(1) / cloud_frac(i, level)) *
                                              (mixt_frac * (AK1 * AK1 + AK1var) +
                                               (real_t(1) - mixt_frac) * (AK2 * AK2 + AK2var)) -
                                          (AKm_over_cloud_frac * AKm_over_cloud_frac)));
        } else
          AKstd_cld(i, level) = zero_threshold;
        AKm_rcm(i, level) = K_one * Kokkos::max(zero_threshold, rcm(i, level) - r_crit);
        // The Native .001 cloud threshold avoids NaNs for tiny cloud fractions
        // (D. Schanen, 3 June 2009).
        if (cloud_frac(i, level) > real_t(.001))
          AKm_rcc(i, level) =
              cloud_frac(i, level) * K_one *
              Kokkos::max(zero_threshold, rcm(i, level) / cloud_frac(i, level) - r_crit);
        else
          AKm_rcc(i, level) = zero_threshold;
      });
}

namespace {
//-----------------------------------------------------------------------
template <class Rc, class Component, class Weight>
KOKKOS_INLINE_FUNCTION void calc_estimate(int num_samples, real_t mixt_frac,               // In
                                          real_t cloud_frac_1, real_t cloud_frac_2, Rc rc, // In
                                          Component X_mixt_comp_one_lev,
                                          Weight lh_sample_point_weights, // In
                                          bool l_lh_importance_sampling,  // In
                                          real_t coeff, real_t r_crit,    // In
                                          real_t& est_m                   // Out
) {                                                                       // Out
  // Compute a Monte Carlo gridbox estimate; divide by cloud fraction for a
  // within-cloud estimate. The false local option retains the Native alternative
  // plume-weighted averaging branch without introducing a runtime option.
  // num_samples: Number of calls to microphysics (normally=2)
  // mixt_frac: Mixture fraction of Gaussians
  // cloud_frac_1: Cloud fraction associated w/ 1st, 2nd mixture component
  // cloud_frac_2: Cloud fraction associated w/ 1st, 2nd mixture component
  // rc: , & ! n in-cloud values of spec liq water content (when positive) [kg/kg]. w,  & ! n
  //   in-cloud values of vertical velocity (m/s) Nc, & ! n in-cloud values of droplet number (#/mg
  //   air) rr    ! n in-cloud values of specific rain content (g/kg)
  // X_mixt_comp_one_lev: Whether we're in the first or second mixture component
  // lh_sample_point_weights: Weight for cloud weighted sampling
  // l_lh_importance_sampling: Do importance sampling (SILHS)
  // Keep the reference argument to Kokkos::max in the calling execution space.
  const real_t zero                         = clubb_kokkos::zero;
  constexpr bool l_cloud_weighted_averaging = false;
  if (mixt_frac > one || mixt_frac < zero)
    Kokkos::abort("Error in calc_estimate: mixture fraction, mixt_frac, does not lie in [0,1].");
  if (cloud_frac_1 > one || cloud_frac_1 < zero)
    Kokkos::abort("Error in calc_estimate: cloud fraction 1, cloud_frac_1, does not lie in [0,1].");
  if (cloud_frac_2 > one || cloud_frac_2 < zero)
    Kokkos::abort("Error in calc_estimate: cloud fraction 2, cloud_frac_2, does not lie in [0,1].");
  real_t est_m1 = zero;
  real_t est_m2 = zero;
  int n1        = 0;
  int n2        = 0;
  for (int sample = 0; sample < num_samples; ++sample) {
    // M. E. Johnson (1987), p.56: accumulate each PDF mixture component in
    // sample order. Importance weights are applied only when enabled.
    if (X_mixt_comp_one_lev(sample) == 1) {
      if (l_lh_importance_sampling)
        est_m1 = est_m1 +
                 coeff * Kokkos::max(zero, rc(sample) - r_crit) * lh_sample_point_weights(sample);
      else
        est_m1 = est_m1 + coeff * Kokkos::max(zero, rc(sample) - r_crit);
      ++n1;
    } else {
      if (l_lh_importance_sampling)
        est_m2 = est_m2 +
                 coeff * Kokkos::max(zero, rc(sample) - r_crit) * lh_sample_point_weights(sample);
      else
        est_m2 = est_m2 + coeff * Kokkos::max(zero, rc(sample) - r_crit);
      ++n2;
    }
  }
  if (n1 == 0 && n2 == 0)
    Kokkos::abort("Error:  no sample points in calc_estimate");
  if (l_cloud_weighted_averaging) {
    // If a plume has no points, estimate its liquid from the other plume.
    if (!(n1 == 0))
      est_m1 = est_m1 / real_t(n1);
    if (!(n2 == 0))
      est_m2 = est_m2 / real_t(n2);
    if (n1 == 0)
      est_m1 = est_m2;
    if (n2 == 0)
      est_m2 = est_m1;
    est_m = mixt_frac * cloud_frac_1 * est_m1 + (one - mixt_frac) * cloud_frac_2 * est_m2;
  } else
    est_m = (est_m1 + est_m2) / real_t(num_samples);
}
} // namespace
} // namespace clubb_kokkos
