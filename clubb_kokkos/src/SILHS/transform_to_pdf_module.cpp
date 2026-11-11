// SILHS uniform-to-PDF transformations in Native routine/call order.
#include "scratch_workspace.hpp"
#include "SILHS/transform_to_pdf_module.hpp"
#include <limits>
#include <algorithm>
namespace clubb_kokkos {
namespace {
void cdfnorminv(int pdf_dim, int nzt, int ngrdcol, int num_samples, Real4D X_u_all_levs, // In
                Real4D std_normal                                                        // Out
);
}

//-------------------------------------------------------------------------------
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
) {                                                                          // Out
  // Transform a uniform sample to CLUBB's PDF. Larson et al.(2005),JAS62,
  // 4010–4026; https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:uniform2pdf.
  // chi,eta,w are normal; later variates are lognormal hydrometeors/nuclei.
  // nzt: Number of vertical grid levels
  // ngrdcol: Number of grid columns
  // num_samples: Number of subcolumn samples
  // pdf_dim: `d' Number of variates (normally 3 + microphysics specific variables)
  // d_uniform_extra: Number of variates included in uniform sample only (often 2)
  // Sigma_Cholesky1: Correlations Cholesky matrix, 1st component [-]
  // Sigma_Cholesky2: Correlations Cholesky matrix, 2nd component [-]
  // mu1: Means of the hydrometeors,(chi, eta, w, <hydrometeors>), 1st component [units vary]
  // mu2: Means of the hydrometeors,(chi, eta, w, <hydrometeors>), 2nd component [units vary]
  // X_mixt_comp_all_levs: Whether we're in the first or second mixture component
  // X_u_all_levs: Sample drawn from uniform distribution from a particular grid level
  // cloud_frac: Cloud fraction [-]
  // l_in_precip_all_levs: Whether we are in precipitation (T/F)
  // X_nl_all_levs: Sample that is transformed ultimately to normal-lognormal
  // Capture addressable values; epsilon retains Native core-kind precision.
  auto workspace_scope = workspace.scope();
  const real_t zero    = clubb_kokkos::zero;
  const real_t epsilon = std::numeric_limits<real_t>::epsilon();
  const int iiPDF_chi  = hm_metadata.iiPDF_chi;
  const int iiPDF_Ncn  = hm_metadata.iiPDF_Ncn;
  const int p          = std::max({(((hm_metadata.iiPDF_chi + 1))), (((hm_metadata.iiPDF_eta + 1))),
                                   (((hm_metadata.iiPDF_w + 1)))});
  constexpr bool l_clip_extreme_chi_sample_pts = true;
  // Native's l_d_variable_lognormal metadata is never read. The exact same
  // p-based range below controls the actual exponential conversion.
  Real4D std_normal = workspace_scope.view<Real4D>(Kokkos::WithoutInitializing, "std_normal",
                                                   pdf_dim, ngrdcol, nzt, num_samples);
  cdfnorminv(pdf_dim, nzt, ngrdcol, num_samples, X_u_all_levs,      // In
             std_normal                                             // Out
  );                                                                // Out
  multiply_Cholesky(nzt, ngrdcol, num_samples, pdf_dim, std_normal, // In
                    Sigma_Cholesky1, Sigma_Cholesky2,               // In
                    mu1, mu2, X_mixt_comp_all_levs,                 // In
                    X_nl_all_levs                                   // Out
  );                                                                // Out
  Kokkos::parallel_for(
      "transform_to_pdf::lognormal",
      ((-ngrdcol * num_samples * nzt * p + ngrdcol * num_samples * nzt * pdf_dim)),
      KOKKOS_LAMBDA(int cell) {
        const int i                           = ((cell % ngrdcol));
        const int sample                      = (((cell / ngrdcol) % num_samples));
        const int k                           = (((cell / (ngrdcol * num_samples)) % nzt));
        const int variable                    = (((cell / (ngrdcol * num_samples * nzt)) + p));
        X_nl_all_levs(i, sample, k, variable) = Kokkos::exp(X_nl_all_levs(i, sample, k, variable));
      });
  Kokkos::parallel_for(
      "transform_to_pdf::precipitation",
      ((-iiPDF_Ncn * ngrdcol * num_samples * nzt + ngrdcol * num_samples * nzt * pdf_dim -
        ngrdcol * num_samples * nzt)),
      KOKKOS_LAMBDA(int cell) {
        const int i        = ((cell % ngrdcol));
        const int sample   = (((cell / ngrdcol) % num_samples));
        const int k        = (((cell / (ngrdcol * num_samples)) % nzt));
        const int variable = ((iiPDF_Ncn + (cell / (ngrdcol * num_samples * nzt)) + 1));
        if (!l_in_precip_all_levs(i, sample, k))
          X_nl_all_levs(i, sample, k, variable) = zero;
      });
  // Closure clips component cloud fraction under extreme conditions. Force
  // chi to the same saturated/unsaturated side, retaining both strict ties.
  if (l_clip_extreme_chi_sample_pts) {
    Kokkos::parallel_for(
        "transform_to_pdf::extreme_chi", ((ngrdcol * num_samples * nzt)), KOKKOS_LAMBDA(int cell) {
          const int i      = ((cell % ngrdcol));
          const int sample = (((cell / ngrdcol) % num_samples));
          const int k      = ((cell / (ngrdcol * num_samples)));
          if (cloud_frac(i, sample, k) < epsilon)
            X_nl_all_levs(i, sample, k, iiPDF_chi) =
                Kokkos::min(X_nl_all_levs(i, sample, k, iiPDF_chi), zero);
          else if (cloud_frac(i, sample, k) > one - epsilon)
            X_nl_all_levs(i, sample, k, iiPDF_chi) =
                Kokkos::max(X_nl_all_levs(i, sample, k, iiPDF_chi), epsilon);
        });
  }
}

namespace {
//-----------------------------------------------------------------------
void cdfnorminv(int pdf_dim, int nzt, int ngrdcol, int num_samples, Real4D X_u_all_levs, // In
                Real4D std_normal                                                        // Out
) {                                                                                      // Out
  // Inverse normal CDF via Giles(2010), Approximating the erfinv function.
  // It has single-precision accuracy even when accepting double arguments.
  // Native coefficients are unsuffixed default REAL literals: round as float
  // before promotion to core precision. This follows the source kind contract.
  // pdf_dim: `d' Number of variates (normally 3 + microphysics specific variables)
  // nzt: Number of vertical grid levels
  // ngrdcol: Number of grid columns
  // num_samples: Number of subcolumn samples
  const Kokkos::Array<real_t, 9> a = {
      real_t(2.81022636e-8f),  real_t(3.43273939e-7f), real_t(-3.5233877e-6f),
      real_t(-4.39150654e-6f), real_t(.00021858087f),  real_t(-.00125372503f),
      real_t(-.00417768164f),  real_t(.246640727f),    real_t(1.50140941f)};
  const Kokkos::Array<real_t, 9> b = {
      real_t(-.000200214257f), real_t(.000100950558f), real_t(.00134934322f),
      real_t(-.00367342844f),  real_t(.00573950773f),  real_t(-.0076224613f),
      real_t(.00943887047f),   real_t(1.00167406f),    real_t(2.83297682f)};
  Kokkos::parallel_for(
      "cdfnorminv", ((ngrdcol * num_samples * nzt * pdf_dim)), KOKKOS_LAMBDA(int cell) {
        const int i      = ((cell % ngrdcol));
        const int p      = (((cell / ngrdcol) % pdf_dim));
        const int k      = (((cell / (ngrdcol * pdf_dim)) % nzt));
        const int sample = ((cell / (ngrdcol * nzt * pdf_dim)));
        const real_t x   = two * X_u_all_levs(i, sample, k, p) - one;
        real_t w         = -Kokkos::log((one - x) * (one + x));
        if (w < real_t(5)) {
          w = w - real_t(2.5);
          std_normal(p, i, k, sample) =
              sqrt_2 * x *
              ((((((((a[0] * w + a[1]) * w + a[2]) * w + a[3]) * w + a[4]) * w + a[5]) * w + a[6]) *
                    w +
                a[7]) *
                   w +
               a[8]);
        } else {
          w = Kokkos::sqrt(w) - real_t(3);
          std_normal(p, i, k, sample) =
              sqrt_2 * x *
              ((((((((b[0] * w + b[1]) * w + b[2]) * w + b[3]) * w + b[4]) * w + b[5]) * w + b[6]) *
                    w +
                b[7]) *
                   w +
               b[8]);
        }
      });
}
} // namespace

// ltqnorm is inline in the public header for device-callable scalar use.

//-------------------------------------------------------------------------------
void multiply_Cholesky(int nzt, int ngrdcol, int num_samples, int pdf_dim, Real4D std_normal, // In
                       Real4D Sigma_Cholesky1, Real4D Sigma_Cholesky2,                        // In
                       Real3D mu1, Real3D mu2,                                                // In
                       Kokkos::View<int***, Kokkos::LayoutLeft> X_mixt_comp_all_levs,         // In
                       Real4D X_nl_all_levs                                                   // Out
) {                                                                                           // Out
  // X_nl=Sigma_Cholesky*std_normal+mu. Each triangular product accumulates j
  // sequentially from 1 through p, including the component mean first.
  // M.E.Johnson(1987), Multivariate Normal and Related Distributions,pp.50–55.
  // nzt: Number of vertical grid levels
  // ngrdcol: Number of grid columns
  // num_samples: Number of samples
  // pdf_dim: Number of variates (normally=5)
  // std_normal: vector of d-variate standard normal distribution [-]
  // Sigma_Cholesky1: Cholesky factorization of the Sigma matrix, 1st component [units vary]
  // Sigma_Cholesky2: Cholesky factorization of the Sigma matrix, 2nd component [units vary]
  // mu1: d-dimensional column vector of means of Gaussian, 1st component [units vary]
  // mu2: d-dimensional column vector of means of Gaussian, 2nd component [units vary]
  // X_mixt_comp_all_levs: Whether we're in the first or second mixture component
  Kokkos::parallel_for(
      "multiply_Cholesky", ((ngrdcol * num_samples * nzt * pdf_dim)), KOKKOS_LAMBDA(int cell) {
        const int i             = ((cell % ngrdcol));
        const int k             = (((cell / ngrdcol) % nzt));
        const int sample        = (((cell / (ngrdcol * nzt)) % num_samples));
        const int p             = ((cell / (ngrdcol * num_samples * nzt)));
        const bool l_first_comp = X_mixt_comp_all_levs(i, sample, k) == 1;
        real_t X_nl_k_sample_i_tmp;
        if (l_first_comp)
          X_nl_k_sample_i_tmp = mu1(i, k, p);
        else
          X_nl_k_sample_i_tmp = mu2(i, k, p);
        for (int j = 0; j < (p + 1); ++j) {
          if (l_first_comp)
            X_nl_k_sample_i_tmp =
                X_nl_k_sample_i_tmp + Sigma_Cholesky1(j, i, k, p) * std_normal(j, i, k, sample);
          else
            X_nl_k_sample_i_tmp =
                X_nl_k_sample_i_tmp + Sigma_Cholesky2(j, i, k, p) * std_normal(j, i, k, sample);
        }
        X_nl_all_levs(i, sample, k, p) = X_nl_k_sample_i_tmp;
      });
}

//-----------------------------------------------------------------------
void chi_eta_2_rtthl(int nzt, int ngrdcol, int num_samples,                         // In
                     Real2D rt_1, Real2D thl_1,                                     // In
                     Real2D rt_2, Real2D thl_2,                                     // In
                     Real2D crt_1, Real2D cthl_1,                                   // In
                     Real2D crt_2, Real2D cthl_2,                                   // In
                     Real2D mu_chi_1, Real2D mu_chi_2,                              // In
                     Real3D chi, Real3D eta,                                        // In
                     Kokkos::View<int***, Kokkos::LayoutLeft> X_mixt_comp_all_levs, // In
                     Real3D lh_rt, Real3D lh_thl                                    // Out
) {                                                                                 // Out
  // Convert Mellor chi/eta to rt[kg/kg] and thl[K], limiting each temperature
  // deviation to 5K. Native limits the deviation, not the cthl coefficient.
  // nzt: Vertical grid levels
  // ngrdcol: Columns
  // num_samples: Number of subcolumn samples
  // rt_1: n dimensional column vector of rt         [kg/kg]
  // thl_1: n dimensional column vector of thetal     [K]
  // rt_2: n dimensional column vector of rt         [kg/kg]
  // thl_2: n dimensional column vector of thetal     [K]
  // crt_1: Constants from plumes 1 & 2 of rt
  // cthl_1: Constants from plumes 1 & 2 of thetal
  // crt_2: Constants from plumes 1 & 2 of rt
  // cthl_2: Constants from plumes 1 & 2 of thetal
  // mu_chi_1: Mean for chi_1 and chi_2         [kg/kg]
  // mu_chi_2: Mean for chi_1 and chi_2         [kg/kg]
  // chi: [kg/kg]
  // eta: [-]
  // X_mixt_comp_all_levs: Whether we're in the first or second mixture component
  // lh_rt: n-dimensional column vectors of rt and thl, including mean and perturbation
  // lh_thl: n-dimensional column vectors of rt and thl, including mean and perturbation
  constexpr real_t thl_dev_lim = real_t(5);
  Kokkos::parallel_for(
      "chi_eta_2_rtthl", ((ngrdcol * num_samples * nzt)), KOKKOS_LAMBDA(int cell) {
        const int i      = ((cell % ngrdcol));
        const int k      = (((cell / ngrdcol) % nzt));
        const int sample = ((cell / (ngrdcol * nzt)));
        real_t lh_dev_thl_lim;
        if (X_mixt_comp_all_levs(i, sample, k) == 1) {
          lh_rt(i, sample, k) = rt_1(i, k) +
                                (real_t(.5) / crt_1(i, k)) * (chi(i, sample, k) - mu_chi_1(i, k)) +
                                (real_t(.5) / crt_1(i, k)) * eta(i, sample, k);
          lh_dev_thl_lim = (-real_t(.5) / cthl_1(i, k)) * (chi(i, sample, k) - mu_chi_1(i, k)) +
                           (real_t(.5) / cthl_1(i, k)) * eta(i, sample, k);
          lh_dev_thl_lim = Kokkos::max(Kokkos::min(lh_dev_thl_lim, thl_dev_lim), -thl_dev_lim);
          lh_thl(i, sample, k) = thl_1(i, k) + lh_dev_thl_lim;
        } else {
          lh_rt(i, sample, k) = rt_2(i, k) +
                                (real_t(.5) / crt_2(i, k)) * (chi(i, sample, k) - mu_chi_2(i, k)) +
                                (real_t(.5) / crt_2(i, k)) * eta(i, sample, k);
          lh_dev_thl_lim = (-real_t(.5) / cthl_2(i, k)) * (chi(i, sample, k) - mu_chi_2(i, k)) +
                           (real_t(.5) / cthl_2(i, k)) * eta(i, sample, k);
          lh_dev_thl_lim = Kokkos::max(Kokkos::min(lh_dev_thl_lim, thl_dev_lim), -thl_dev_lim);
          lh_thl(i, sample, k) = thl_2(i, k) + lh_dev_thl_lim;
        }
      });
}
} // namespace clubb_kokkos
