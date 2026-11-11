// Microphysics calculations translated from the same-named Native source.
//
// Public coordinators take caller-owned column-first LayoutLeft Views and retain
// Native call, statistics and recurrence order. Scalar helpers execute on the
// Kokkos backend; scratch allocations live only for the current call.
#pragma once
#include "scratch_workspace.hpp"
#include "Microphys/KK_microphys/KK_utilities.hpp"
#include "Microphys/microphys_diagnostics.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/pdf_parameter_module.hpp"
#include "CLUBB_core/hydromet_pdf_parameter_module.hpp"
#include "CLUBB_core/corr_varnce_module.hpp"
#include "CLUBB_core/pdf_utilities.hpp"
#include "CLUBB_core/index_mapping.hpp"
#include "CLUBB_core/T_in_K_module.hpp"
#include "CLUBB_core/saturation.hpp"
#include "CLUBB_core/stats_netcdf.hpp"
#include <Kokkos_MathematicalFunctions.hpp>
#include <stdexcept>
#include <limits>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "CLUBB_core/error_code.hpp"
#include "CLUBB_core/parameter_indices.hpp"
#include "CLUBB_core/parameters_tunable.hpp"
#include "CLUBB_core/fill_holes.hpp"
#include "CLUBB_core/advance_helper_module.hpp"
#include "CLUBB_core/matrix_solver_wrapper.hpp"
#include "CLUBB_core/diffusion.hpp"
#include "CLUBB_core/mean_adv.hpp"
#include "Microphys/parameters_microphys.hpp"
namespace clubb_kokkos {
using Int1D = Kokkos::View<int*, Kokkos::LayoutLeft>;
using Int3D = Kokkos::View<int***, Kokkos::LayoutLeft>;

inline void hydrometeor_mixed_moments(ScratchWorkspace& workspace, const grid& gr, int ngrdcol,
                                      int nzt, int pdf_dim, int hydromet_dim, Real3D hydromet,
                                      const hm_metadata_type& hm_metadata, Real3D mu_x_1_n,
                                      Real3D mu_x_2_n, Real3D sigma_x_1_n, Real3D sigma_x_2_n,
                                      Real4D corr_array_1_n, Real4D corr_array_2_n,
                                      const pdf_parameter& pdf_params,
                                      const hydromet_pdf_parameter_collection& hydromet_pdf_params,
                                      const precipitation_fractions& precip_fracs,
                                      stats_type& stats, Real3D rtphmp_zt, Real3D thlphmp_zt,
                                      Real3D wp2hmp);
KOKKOS_INLINE_FUNCTION real_t xphmp_integral_covar(real_t mu_x_1, real_t mu_x_2, real_t mu_hm_1,
                                                   real_t mu_hm_2, real_t sigma_x_1,
                                                   real_t sigma_x_2, real_t sigma_hm_1,
                                                   real_t sigma_hm_2, real_t corr_x_hm_1,
                                                   real_t corr_x_hm_2, real_t mixt_frac,
                                                   real_t precip_frac_1, real_t precip_frac_2,
                                                   real_t x_mean, real_t x_tol, real_t hm_tol);
KOKKOS_INLINE_FUNCTION real_t hmxphmyp_integral_covar(
    real_t mu_hmx_1, real_t mu_hmx_2, real_t mu_hmy_1, real_t mu_hmy_2, real_t sigma_hmx_1,
    real_t sigma_hmx_2, real_t sigma_hmy_1, real_t sigma_hmy_2, real_t corr_hmx_hmy_1,
    real_t corr_hmx_hmy_2, real_t mixt_frac, real_t precip_frac_1, real_t precip_frac_2,
    real_t hmx_mean, real_t hmy_mean, real_t hmx_tol, real_t hmy_tol);
KOKKOS_INLINE_FUNCTION real_t xp_a_hmpb_integrals_all_MM(
    real_t mu_x_1, real_t mu_x_2, real_t mu_hm_1, real_t mu_hm_2, real_t mu_hm_1_n,
    real_t mu_hm_2_n, real_t sigma_x_1, real_t sigma_x_2, real_t sigma_hm_1, real_t sigma_hm_2,
    real_t sigma_hm_1_n, real_t sigma_hm_2_n, real_t corr_x_hm_1_n, real_t corr_x_hm_2_n,
    real_t mixt_frac, real_t precip_frac_1, real_t precip_frac_2, real_t x_mean, real_t hm_mean,
    real_t x_tol, real_t hm_tol, int a_exp, int b_exp);
KOKKOS_INLINE_FUNCTION real_t bivar_NL_x_hm_all_MM_comp_eq(
    real_t mu_x_i, real_t mu_hm_i, real_t mu_hm_i_n, real_t sigma_x_i, real_t sigma_hm_i,
    real_t sigma_hm_i_n, real_t corr_x_hm_i_n, real_t precip_frac_i, real_t x_mean, real_t hm_mean,
    real_t x_tol, real_t hm_tol, int a_exp, int b_exp);
KOKKOS_INLINE_FUNCTION real_t bivar_NL_int_PDF_comp_all_MM(real_t mu_x1_i, real_t mu_x2_i_n,
                                                           real_t sigma_x1_i, real_t sigma_x2_i_n,
                                                           real_t corr_x1_x2_i_n, real_t x1_mean,
                                                           real_t x2_mean, int a_exp, int b_exp);
KOKKOS_INLINE_FUNCTION real_t univar_N_int_PDF_comp_all_MM(real_t mu_x_i, real_t sigma_x_i,
                                                           real_t x_mean, int a_exp);
KOKKOS_INLINE_FUNCTION real_t univar_L_int_PDF_comp_all_MM(real_t mu_x_i_n, real_t sigma_x_i_n,
                                                           real_t x_mean, int b_exp);
inline void hydrometeor_mixed_moments(ScratchWorkspace& workspace, const grid& gr, int ngrdcol,
                                      int nzt, int pdf_dim, int hydromet_dim, Real3D hydromet,
                                      const hm_metadata_type& hm_metadata, Real3D mu_x_1_n,
                                      Real3D mu_x_2_n, Real3D sigma_x_1_n, Real3D sigma_x_2_n,
                                      Real4D corr_array_1_n, Real4D corr_array_2_n,
                                      const pdf_parameter& pdf_params,
                                      const hydromet_pdf_parameter_collection& hydromet_pdf_params,
                                      const precipitation_fractions& precip_fracs,
                                      stats_type& stats, Real3D rtphmp_zt, Real3D thlphmp_zt,
                                      Real3D wp2hmp) {
  // Description:
  // Calculates <rt'hm'>, <thl'hm'>, and <w'^2 hm'> for every model column
  // and hydrometeor. The column index is innermost in each calculation loop.
  // These terms are used in the liquid/ice water loading term as part of the
  // buoyancy term in some of the CLUBB predictive equations.
  // References:
  //-----------------------------------------------------------------------
  //------------------------ Input Variables ------------------------
  //------------------------ Input/Output Variables ------------------------
  //------------------------ Output Variables ------------------------
  //------------------------ Local Variables ------------------------
  // Unpacked parameters.
  // Calculated or recalculated values.
  auto workspace_scope = workspace.scope();
  const bool l_sample = stats.l_sample;

  const auto hm_metadata_hydromet_list_host = to_host_space(workspace, hm_metadata.hydromet_list);
  real_t mu_hm_1;       // Mean of hm (1st PDF component) in-precip (ip) [hm units]
  real_t mu_hm_2;       // Mean of hm (2nd PDF component) ip             [hm units]
  real_t sigma_hm_1;    // Standard deviation of hm (1st PDF comp.) ip   [hm units]
  real_t sigma_hm_2;    // Standard deviation of hm (2nd PDF comp.) ip   [hm units]
  real_t corr_chi_hm_1; // Correlation of chi and hm (1st PDF component) ip    [-]
  real_t corr_chi_hm_2; // Correlation of chi and hm (2nd PDF component) ip    [-]
  real_t corr_eta_hm_1; // Correlation of eta and hm (1st PDF component) ip    [-]
  real_t corr_eta_hm_2; // Correlation of eta and hm (2nd PDF component) ip    [-]
  real_t mixt_frac;     // Mixture fraction                                    [-]
  real_t precip_frac_1; // Precipitation fraction (1st PDF component)          [-]
  real_t precip_frac_2; // Precipitation fraction (2nd PDF component)          [-]
  real_t crt_1;         // Coef. of rt in chi/eta eqns. (1st PDF component)    [-]
  real_t crt_2;         // Coef. of rt in chi/eta eqns. (2nd PDF component)    [-]
  real_t cthl_1;        // Coef. of thl: chi/eta eqns. (1st PDF comp.) [(kg/kg)/K]
  real_t cthl_2;        // Coef. of thl: chi/eta eqns. (2nd PDF comp.) [(kg/kg)/K]
  Real4D hmxphmyp_zt = workspace_scope.view<Real4D>(
      "hydrometeor_mixed_moments::hmxphmyp_zt", ngrdcol, nzt, hydromet_dim,
      hydromet_dim); // Covariance (overall) of two hydrometeors  [hmx*hmy units]
  Real1D hydromet_tols = workspace_scope.view<Real1D>(
      "hydrometeor_mixed_moments::hydromet_tols",
      hydromet_dim); // Hydrometeor tolerance by species                 [units vary]
  int hm_idx;        // Index of a hydrometeor in the hydrometeor set of indices
  int hmy_idx;       // Index of second hydrometeor in the hydrometeor set of indices
  Int1D pdf_idx_by_hm = workspace_scope.view<Int1D>(
      "hydrometeor_mixed_moments::pdf_idx_by_hm", hydromet_dim); // PDF index by hydrometeor species
  std::string hm_type;
  std::string hmx_type;
  std::string hmy_type;
  std::string var_name;
  const auto hm_metadata_iiPDF_chi              = hm_metadata.iiPDF_chi;
  const auto hm_metadata_iiPDF_eta              = hm_metadata.iiPDF_eta;
  const auto hm_metadata_iiPDF_w                = hm_metadata.iiPDF_w;
  const auto pdf_params_rt_1                    = pdf_params.rt_1;
  const auto pdf_params_rt_2                    = pdf_params.rt_2;
  const auto pdf_params_varnce_rt_1             = pdf_params.varnce_rt_1;
  const auto pdf_params_varnce_rt_2             = pdf_params.varnce_rt_2;
  const auto pdf_params_thl_1                   = pdf_params.thl_1;
  const auto pdf_params_thl_2                   = pdf_params.thl_2;
  const auto pdf_params_varnce_thl_1            = pdf_params.varnce_thl_1;
  const auto pdf_params_varnce_thl_2            = pdf_params.varnce_thl_2;
  const auto pdf_params_crt_1                   = pdf_params.crt_1;
  const auto pdf_params_crt_2                   = pdf_params.crt_2;
  const auto pdf_params_cthl_1                  = pdf_params.cthl_1;
  const auto pdf_params_cthl_2                  = pdf_params.cthl_2;
  const auto pdf_params_mixt_frac               = pdf_params.mixt_frac;
  const auto hydromet_pdf_params_mu_hm_1        = hydromet_pdf_params.mu_hm_1;
  const auto hydromet_pdf_params_mu_hm_2        = hydromet_pdf_params.mu_hm_2;
  const auto hydromet_pdf_params_sigma_hm_1     = hydromet_pdf_params.sigma_hm_1;
  const auto hydromet_pdf_params_sigma_hm_2     = hydromet_pdf_params.sigma_hm_2;
  const auto hydromet_pdf_params_corr_chi_hm_1  = hydromet_pdf_params.corr_chi_hm_1;
  const auto hydromet_pdf_params_corr_chi_hm_2  = hydromet_pdf_params.corr_chi_hm_2;
  const auto hydromet_pdf_params_corr_eta_hm_1  = hydromet_pdf_params.corr_eta_hm_1;
  const auto hydromet_pdf_params_corr_eta_hm_2  = hydromet_pdf_params.corr_eta_hm_2;
  const auto hydromet_pdf_params_corr_hmx_hmy_1 = hydromet_pdf_params.corr_hmx_hmy_1;
  const auto hydromet_pdf_params_corr_hmx_hmy_2 = hydromet_pdf_params.corr_hmx_hmy_2;
  const auto precip_fracs_precip_frac_1         = precip_fracs.precip_frac_1;
  const auto precip_fracs_precip_frac_2         = precip_fracs.precip_frac_2;

  //--------------------------- Begin Code ---------------------------

  // Species-index mapping is host metadata; stage the resulting indices.
  auto pdf_idx_by_hm_host = Kokkos::create_mirror_view(pdf_idx_by_hm);
  auto hydromet_tols_host = Kokkos::create_mirror_view(hydromet_tols);
  auto tolerance_values   = to_host_space(workspace, hm_metadata.hydromet_tol);
  for (int hm_idx = 0; hm_idx < hydromet_dim; ++hm_idx) {
    pdf_idx_by_hm_host(hm_idx) = ((hydromet2pdf_idx(hm_idx, hm_metadata) + 1)) - 1;
    hydromet_tols_host(hm_idx) = tolerance_values(hm_idx);
  }
  Kokkos::deep_copy(pdf_idx_by_hm, pdf_idx_by_hm_host);
  Kokkos::deep_copy(hydromet_tols, hydromet_tols_host);

  // Calculate <rt'hm'>, <thl'hm'>, and <w'^2 hm'> for each hydrometeor
  // species.
  // Full fields are stored as (column,level,species), so keeping i innermost
  // makes each hot-loop traversal contiguous and vectorization-friendly.
  Kokkos::parallel_for(
      "hydrometeor_mixed_moments::loop2", ngrdcol, KOKKOS_LAMBDA(int column) {
        real_t mu_w_1;        // Mean of w (1st PDF component)                      [m/s]
        real_t mu_w_2;        // Mean of w (2nd PDF component)                      [m/s]
        real_t mu_rt_1;       // Mean of rt (1st PDF component)                   [kg/kg]
        real_t mu_rt_2;       // Mean of rt (2nd PDF component)                   [kg/kg]
        real_t mu_thl_1;      // Mean of thl (1st PDF component)                      [K]
        real_t mu_thl_2;      // Mean of thl (2nd PDF component)                      [K]
        real_t mu_hm_1;       // Mean of hm (1st PDF component) in-precip (ip) [hm units]
        real_t mu_hm_2;       // Mean of hm (2nd PDF component) ip             [hm units]
        real_t mu_hm_1_n;     // Mean of ln hm (1st PDF component) ip      [ln(hm units)]
        real_t mu_hm_2_n;     // Mean of ln hm (2nd PDF component) ip      [ln(hm units)]
        real_t sigma_w_1;     // Standard deviation of w (1st PDF component)        [m/s]
        real_t sigma_w_2;     // Standard deviation of w (2nd PDF component)        [m/s]
        real_t sigma_rt_1;    // Standard deviation of rt (1st PDF component)     [kg/kg]
        real_t sigma_rt_2;    // Standard deviation of rt (2nd PDF component)     [kg/kg]
        real_t sigma_thl_1;   // Standard deviation of thl (1st PDF component)        [K]
        real_t sigma_thl_2;   // Standard deviation of thl (2nd PDF component)        [K]
        real_t sigma_chi_1;   // Standard deviation of chi (1st PDF component)    [kg/kg]
        real_t sigma_chi_2;   // Standard deviation of chi (2nd PDF component)    [kg/kg]
        real_t sigma_eta_1;   // Standard deviation of eta (1st PDF component)    [kg/kg]
        real_t sigma_eta_2;   // Standard deviation of eta (2nd PDF component)    [kg/kg]
        real_t sigma_hm_1;    // Standard deviation of hm (1st PDF comp.) ip   [hm units]
        real_t sigma_hm_2;    // Standard deviation of hm (2nd PDF comp.) ip   [hm units]
        real_t sigma_hm_1_n;  // Standard deviation of ln hm (1st PDF component) ip   [-]
        real_t sigma_hm_2_n;  // Standard deviation of ln hm (2nd PDF component) ip   [-]
        real_t corr_chi_hm_1; // Correlation of chi and hm (1st PDF component) ip    [-]
        real_t corr_chi_hm_2; // Correlation of chi and hm (2nd PDF component) ip    [-]
        real_t corr_eta_hm_1; // Correlation of eta and hm (1st PDF component) ip    [-]
        real_t corr_eta_hm_2; // Correlation of eta and hm (2nd PDF component) ip    [-]
        real_t corr_w_hm_1_n; // Correlation of w and ln hm (1st PDF component) ip   [-]
        real_t corr_w_hm_2_n; // Correlation of w and ln hm (2nd PDF component) ip   [-]
        real_t mixt_frac;     // Mixture fraction                                    [-]
        real_t precip_frac_1; // Precipitation fraction (1st PDF component)          [-]
        real_t precip_frac_2; // Precipitation fraction (2nd PDF component)          [-]
        real_t crt_1;         // Coef. of rt in chi/eta eqns. (1st PDF component)    [-]
        real_t crt_2;         // Coef. of rt in chi/eta eqns. (2nd PDF component)    [-]
        real_t cthl_1;        // Coef. of thl: chi/eta eqns. (1st PDF comp.) [(kg/kg)/K]
        real_t cthl_2;        // Coef. of thl: chi/eta eqns. (2nd PDF comp.) [(kg/kg)/K]
        real_t hm_mean;       // Mean of hydrometeor, hm (overall)            [hm units]
        real_t hm_tol;        // Tolerance value of hydrometeor, hm           [hm units]
        real_t corr_rt_hm_1;  // Correlation of rt and hm (1st PDF component)        [-]
        real_t corr_rt_hm_2;  // Correlation of rt and hm (2nd PDF component)        [-]
        real_t corr_thl_hm_1; // Correlation of thl and hm (1st PDF component)       [-]
        real_t corr_thl_hm_2; // Correlation of thl and hm (2nd PDF component)       [-]
        real_t wm;            // Mean of w (overall)                               [m/s]
        real_t rtm;           // Mean of rt (overall)                            [kg/kg]
        real_t thlm;          // Mean of thl (overall)                               [K]
        int hm_idx;           // Index of a hydrometeor in the hydrometeor set of indices
        int a_exp;            // Exponent on w' in < w'^a hm'^b >
        int b_exp;            // Exponent on hm' in < w'^a hm'^b >
        int k;                // Index of a vertical level
        const int i = column;
        for (hm_idx = 0; hm_idx < hydromet_dim; hm_idx += 1) {
          hm_tol = hydromet_tols(hm_idx);
          for (k = 0; k < nzt; k += 1) {
            // Unpack the PDF and precipitation parameters for this column.
            mu_rt_1  = pdf_params_rt_1((((i + 1) - 1)), k);
            mu_rt_2  = pdf_params_rt_2((((i + 1) - 1)), k);
            mu_thl_1 = pdf_params_thl_1((((i + 1) - 1)), k);
            mu_thl_2 = pdf_params_thl_2((((i + 1) - 1)), k);
            // Unpack the standard deviations of rt and thl in each PDF component.
            sigma_rt_1  = Kokkos::sqrt(pdf_params_varnce_rt_1((((i + 1) - 1)), k));
            sigma_rt_2  = Kokkos::sqrt(pdf_params_varnce_rt_2((((i + 1) - 1)), k));
            sigma_thl_1 = Kokkos::sqrt(pdf_params_varnce_thl_1((((i + 1) - 1)), k));
            sigma_thl_2 = Kokkos::sqrt(pdf_params_varnce_thl_2((((i + 1) - 1)), k));
            // Unpack the mixture fraction.
            mixt_frac = pdf_params_mixt_frac((((i + 1) - 1)), k);
            // Unpack the precipitation fraction in each PDF component.
            precip_frac_1 = precip_fracs_precip_frac_1((((i + 1) - 1)), k);
            precip_frac_2 = precip_fracs_precip_frac_2((((i + 1) - 1)), k);
            // Unpack the coefficients of rt and thl in the chi/eta PDF transformation
            // equations for each PDF component.
            crt_1  = pdf_params_crt_1((((i + 1) - 1)), k);
            crt_2  = pdf_params_crt_2((((i + 1) - 1)), k);
            cthl_1 = pdf_params_cthl_1((((i + 1) - 1)), k);
            cthl_2 = pdf_params_cthl_2((((i + 1) - 1)), k);
            // Unpack the means and standard deviations of w in each PDF component.
            mu_w_1    = mu_x_1_n((((i + 1) - 1)), k, hm_metadata_iiPDF_w);
            mu_w_2    = mu_x_2_n((((i + 1) - 1)), k, hm_metadata_iiPDF_w);
            sigma_w_1 = sigma_x_1_n((((i + 1) - 1)), k, hm_metadata_iiPDF_w);
            sigma_w_2 = sigma_x_2_n((((i + 1) - 1)), k, hm_metadata_iiPDF_w);
            // Unpack the standard deviations of chi and eta in each PDF component.
            sigma_chi_1 = sigma_x_1_n((((i + 1) - 1)), k, hm_metadata_iiPDF_chi);
            sigma_chi_2 = sigma_x_2_n((((i + 1) - 1)), k, hm_metadata_iiPDF_chi);
            sigma_eta_1 = sigma_x_1_n((((i + 1) - 1)), k, hm_metadata_iiPDF_eta);
            sigma_eta_2 = sigma_x_2_n((((i + 1) - 1)), k, hm_metadata_iiPDF_eta);

            // Re-calculate rtm, thlm, and wm from PDF parameters.
            // This needs to be done because rtm and thlm have been advanced since
            // the PDF parameters have been calculated.  It is necessary to use values
            // of the mean fields (rtm, thlm, and wm) that are consistent with the
            // PDF parameters.  This does not need to be done for hydromet because
            // hydrometeors have not been advanced since the hydrometeor PDF
            // parameters were set up.
            wm   = compute_mean_binormal(mu_w_1, mu_w_2, mixt_frac);
            rtm  = compute_mean_binormal(mu_rt_1, mu_rt_2, mixt_frac);
            thlm = compute_mean_binormal(mu_thl_1, mu_thl_2, mixt_frac);

            // Unpack the mean (in-precip) of hm in each PDF component.
            mu_hm_1 = hydromet_pdf_params_mu_hm_1((((i + 1) - 1)), k, hm_idx);
            mu_hm_2 = hydromet_pdf_params_mu_hm_2((((i + 1) - 1)), k, hm_idx);
            // Unpack the standard deviation (in-precip) of hm in each PDF
            // component.
            sigma_hm_1 = hydromet_pdf_params_sigma_hm_1((((i + 1) - 1)), k, hm_idx);
            sigma_hm_2 = hydromet_pdf_params_sigma_hm_2((((i + 1) - 1)), k, hm_idx);
            // Calculate the correlation (in-precip) of rt/thl and hm for each PDF
            // component.  Since CLUBB uses a PDF transformation from rt and
            // theta-l coordinates to chi and eta coordinates for each PDF
            // component, the correlation arrays are written in terms of chi and
            // eta correlations.  This makes a calculation necessary for these
            // correlations.
            corr_chi_hm_1 = hydromet_pdf_params_corr_chi_hm_1((((i + 1) - 1)), k, hm_idx);
            corr_chi_hm_2 = hydromet_pdf_params_corr_chi_hm_2((((i + 1) - 1)), k, hm_idx);
            corr_eta_hm_1 = hydromet_pdf_params_corr_eta_hm_1((((i + 1) - 1)), k, hm_idx);
            corr_eta_hm_2 = hydromet_pdf_params_corr_eta_hm_2((((i + 1) - 1)), k, hm_idx);
            corr_rt_hm_1  = calc_corr_rt_x(crt_1, sigma_rt_1, sigma_chi_1, sigma_eta_1,
                                           corr_chi_hm_1, corr_eta_hm_1);
            corr_rt_hm_2  = calc_corr_rt_x(crt_2, sigma_rt_2, sigma_chi_2, sigma_eta_2,
                                           corr_chi_hm_2, corr_eta_hm_2);
            corr_thl_hm_1 = calc_corr_thl_x(cthl_1, sigma_thl_1, sigma_chi_1, sigma_eta_1,
                                            corr_chi_hm_1, corr_eta_hm_1);
            corr_thl_hm_2 = calc_corr_thl_x(cthl_2, sigma_thl_2, sigma_chi_2, sigma_eta_2,
                                            corr_chi_hm_2, corr_eta_hm_2);

            // Calculate <rt'hm'>.
            rtphmp_zt((((i + 1) - 1)), k, hm_idx) =
                xphmp_integral_covar(mu_rt_1, mu_rt_2, mu_hm_1, mu_hm_2, sigma_rt_1, sigma_rt_2,
                                     sigma_hm_1, sigma_hm_2, corr_rt_hm_1, corr_rt_hm_2, mixt_frac,
                                     precip_frac_1, precip_frac_2, rtm, rt_tol, hm_tol);
            // Calculate <thl'hm'>.
            thlphmp_zt((((i + 1) - 1)), k, hm_idx) = xphmp_integral_covar(
                mu_thl_1, mu_thl_2, mu_hm_1, mu_hm_2, sigma_thl_1, sigma_thl_2, sigma_hm_1,
                sigma_hm_2, corr_thl_hm_1, corr_thl_hm_2, mixt_frac, precip_frac_1, precip_frac_2,
                thlm, thl_tol, hm_tol);

            // Unpack the mean (in-precip) of ln hm in each PDF component.
            mu_hm_1_n = mu_x_1_n((((i + 1) - 1)), k, ((pdf_idx_by_hm(hm_idx) + 1) - 1));
            mu_hm_2_n = mu_x_2_n((((i + 1) - 1)), k, ((pdf_idx_by_hm(hm_idx) + 1) - 1));
            // Unpack the standard deviation (in-precip) of ln hm in each PDF
            // component.
            sigma_hm_1_n = sigma_x_1_n((((i + 1) - 1)), k, ((pdf_idx_by_hm(hm_idx) + 1) - 1));
            sigma_hm_2_n = sigma_x_2_n((((i + 1) - 1)), k, ((pdf_idx_by_hm(hm_idx) + 1) - 1));
            // Unpack the correlation (in-precip) of w and ln hm in each PDF
            // component.
            corr_w_hm_1_n = corr_array_1_n((((i + 1) - 1)), k, ((pdf_idx_by_hm(hm_idx) + 1) - 1),
                                           hm_metadata_iiPDF_w);
            corr_w_hm_2_n = corr_array_2_n((((i + 1) - 1)), k, ((pdf_idx_by_hm(hm_idx) + 1) - 1),
                                           hm_metadata_iiPDF_w);
            // Unpack the mean (overall) value of the hydrometeor.
            hm_mean = hydromet((((i + 1) - 1)), k, hm_idx);
            // The general form of the mixed moment equation is <w'^a hm'^b>.
            // For <w'^2 hm'>, a = 2 and b = 1.
            a_exp = 2;
            b_exp = 1;
            // Calculate <w'^2 hm'>.
            wp2hmp((((i + 1) - 1)), k, hm_idx) = xp_a_hmpb_integrals_all_MM(
                mu_w_1, mu_w_2, mu_hm_1, mu_hm_2, mu_hm_1_n, mu_hm_2_n, sigma_w_1, sigma_w_2,
                sigma_hm_1, sigma_hm_2, sigma_hm_1_n, sigma_hm_2_n, corr_w_hm_1_n, corr_w_hm_2_n,
                mixt_frac, precip_frac_1, precip_frac_2, wm, hm_mean, w_tol, hm_tol, a_exp, b_exp);
          }
        }
      });

  Kokkos::parallel_for(
      "hydrometeor_mixed_moments::loop3", ngrdcol, KOKKOS_LAMBDA(int column) {
        real_t mu_hm_1;       // Mean of hm (1st PDF component) in-precip (ip) [hm units]
        real_t mu_hm_2;       // Mean of hm (2nd PDF component) ip             [hm units]
        real_t mu_hmy_1;      // Mean of hmy (1st PDF component) ip           [hmy units]
        real_t mu_hmy_2;      // Mean of hmy (2nd PDF component) ip           [hmy units]
        real_t sigma_hm_1;    // Standard deviation of hm (1st PDF comp.) ip   [hm units]
        real_t sigma_hm_2;    // Standard deviation of hm (2nd PDF comp.) ip   [hm units]
        real_t sigma_hmy_1;   // Standard deviation of hmy (1st PDF comp.) ip [hmy units]
        real_t sigma_hmy_2;   // Standard deviation of hmy (2nd PDF comp.) ip [hmy units]
        real_t corr_hm_hmy_1; // Correlation of hm and hmy (1st PDF component) ip    [-]
        real_t corr_hm_hmy_2; // Correlation of hm and hmy (2nd PDF component) ip    [-]
        real_t mixt_frac;     // Mixture fraction                                    [-]
        real_t precip_frac_1; // Precipitation fraction (1st PDF component)          [-]
        real_t precip_frac_2; // Precipitation fraction (2nd PDF component)          [-]
        real_t hm_mean;       // Mean of hydrometeor, hm (overall)            [hm units]
        real_t hmy_mean;      // Mean of second hydrometeor, hmy (overall)   [hmy units]
        real_t hm_tol;        // Tolerance value of hydrometeor, hm           [hm units]
        real_t hmy_tol;       // Tolerance value of second hydrometeor, hmy  [hmy units]
        int hm_idx;           // Index of a hydrometeor in the hydrometeor set of indices
        int hmy_idx;          // Index of second hydrometeor in the hydrometeor set of indices
        int k;                // Index of a vertical level
        const int i = column;
        for (hm_idx = 0; hm_idx < hydromet_dim; hm_idx += 1) {
          hm_tol = hydromet_tols(hm_idx);
          // Calculate the covariance (overall) of two hydrometeors, <hmx'hmy'>,
          // for each unique set of two different hydrometeors.
          for (hmy_idx = (hm_idx + 1); hmy_idx < hydromet_dim; hmy_idx += 1) {
            // Unpack the tolerance value for the second hydrometeor, hmy.
            hmy_tol = hydromet_tols(hmy_idx);
            for (k = 0; k < nzt; k += 1) {
              mixt_frac     = pdf_params_mixt_frac((((i + 1) - 1)), k);
              precip_frac_1 = precip_fracs_precip_frac_1((((i + 1) - 1)), k);
              precip_frac_2 = precip_fracs_precip_frac_2((((i + 1) - 1)), k);
              mu_hm_1       = hydromet_pdf_params_mu_hm_1((((i + 1) - 1)), k, hm_idx);
              mu_hm_2       = hydromet_pdf_params_mu_hm_2((((i + 1) - 1)), k, hm_idx);
              // Unpack the mean (in-precip) of the second hydrometeor, hmy, in
              // each PDF component.
              mu_hmy_1   = hydromet_pdf_params_mu_hm_1((((i + 1) - 1)), k, hmy_idx);
              mu_hmy_2   = hydromet_pdf_params_mu_hm_2((((i + 1) - 1)), k, hmy_idx);
              sigma_hm_1 = hydromet_pdf_params_sigma_hm_1((((i + 1) - 1)), k, hm_idx);
              sigma_hm_2 = hydromet_pdf_params_sigma_hm_2((((i + 1) - 1)), k, hm_idx);
              // Unpack the standard deviation (in-precip) of hmy in each PDF
              // component.
              sigma_hmy_1 = hydromet_pdf_params_sigma_hm_1((((i + 1) - 1)), k, hmy_idx);
              sigma_hmy_2 = hydromet_pdf_params_sigma_hm_2((((i + 1) - 1)), k, hmy_idx);
              // Unpack the correlation (in-precip) of hm and hmy in each PDF
              // component.
              corr_hm_hmy_1 =
                  hydromet_pdf_params_corr_hmx_hmy_1((((i + 1) - 1)), k, hm_idx, hmy_idx);
              corr_hm_hmy_2 =
                  hydromet_pdf_params_corr_hmx_hmy_2((((i + 1) - 1)), k, hm_idx, hmy_idx);
              hm_mean = hydromet((((i + 1) - 1)), k, hm_idx);
              // Unpack the mean (overall) value of hmy.
              hmy_mean = hydromet((((i + 1) - 1)), k, hmy_idx);
              // Calculate the covariance <hmx'hmy'>.
              hmxphmyp_zt((((i + 1) - 1)), k, hmy_idx, hm_idx) = hmxphmyp_integral_covar(
                  mu_hm_1, mu_hm_2, mu_hmy_1, mu_hmy_2, sigma_hm_1, sigma_hm_2, sigma_hmy_1,
                  sigma_hmy_2, corr_hm_hmy_1, corr_hm_hmy_2, mixt_frac, precip_frac_1,
                  precip_frac_2, hm_mean, hmy_mean, hm_tol, hmy_tol);
            }
          }
        }
      });

  // Statistics
  if (stats.l_sample) {
    for (hm_idx = 0; hm_idx < hydromet_dim; hm_idx += 1) {
      hm_type = std::string(hm_metadata_hydromet_list_host.data() +
                                (hm_metadata_hydromet_list_host.extent(0) * hm_idx),
                            hm_metadata_hydromet_list_host.extent(0));

      var_name =
          ((std::string("wp2") + microphys_detail::trim(hm_type.substr(0, 2))) + std::string("p"));
      stats_update(var_name, Kokkos::subview(wp2hmp, Kokkos::ALL, Kokkos::ALL, hm_idx), stats);

      var_name =
          ((std::string("rtp") + microphys_detail::trim(hm_type.substr(0, 2))) + std::string("p"));
      stats_update(var_name,
                   zt2zm_api(gr.nzm, nzt, ngrdcol, gr,
                             Kokkos::subview(rtphmp_zt, Kokkos::ALL, Kokkos::ALL, hm_idx)),
                   stats);

      var_name =
          ((std::string("thlp") + microphys_detail::trim(hm_type.substr(0, 2))) + std::string("p"));
      stats_update(var_name,
                   zt2zm_api(gr.nzm, nzt, ngrdcol, gr,
                             Kokkos::subview(thlphmp_zt, Kokkos::ALL, Kokkos::ALL, hm_idx)),
                   stats);
    }

    for (hm_idx = 0; hm_idx < hydromet_dim; hm_idx += 1) {
      hmx_type = std::string(hm_metadata_hydromet_list_host.data() +
                                 (hm_metadata_hydromet_list_host.extent(0) * hm_idx),
                             hm_metadata_hydromet_list_host.extent(0));
      for (hmy_idx = (hm_idx + 1); hmy_idx < hydromet_dim; hmy_idx += 1) {
        hmy_type = std::string(hm_metadata_hydromet_list_host.data() +
                                   (hm_metadata_hydromet_list_host.extent(0) * hmy_idx),
                               hm_metadata_hydromet_list_host.extent(0));
        var_name = (((microphys_detail::trim(hmx_type.substr(0, 2)) + std::string("p")) +
                     microphys_detail::trim(hmy_type.substr(0, 2))) +
                    std::string("p"));

        const auto array_actual_6 =
            zt2zm_api(gr.nzm, nzt, ngrdcol, gr,
                      Kokkos::subview(hmxphmyp_zt, Kokkos::ALL, Kokkos::ALL, hmy_idx, hm_idx));
        stats_update(var_name, array_actual_6, stats);
      }
    }
  }

  return;
}
KOKKOS_INLINE_FUNCTION real_t xphmp_integral_covar(real_t mu_x_1, real_t mu_x_2, real_t mu_hm_1,
                                                   real_t mu_hm_2, real_t sigma_x_1,
                                                   real_t sigma_x_2, real_t sigma_hm_1,
                                                   real_t sigma_hm_2, real_t corr_x_hm_1,
                                                   real_t corr_x_hm_2, real_t mixt_frac,
                                                   real_t precip_frac_1, real_t precip_frac_2,
                                                   real_t x_mean, real_t x_tol, real_t hm_tol) {
  // Description:
  // Solves for the covariance < x'hm' >.  The variable "x" is a variable that
  // is distributed binormally (meaning that it has a normally-distributed
  // individual marginal in each PDF component).  This applies to such
  // variables as w, rt, thl, sclr, etc.  The variable "hm" stands for any
  // precipitating hydrometeor.  This applies to such variables as rr, Nr, ri,
  // Ni, etc.
  //
  // The covariance < x'hm' > can also be found by passing a = 1 and b = 1 into
  // function xp_a_hmpb_integrals_all_MM, which is found below.  However,
  // the code found here has been streamlined for the special case of
  // covariances.
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable

  real_t xphmp_integral_covar;

  if ((((sigma_x_1 <= x_tol) || (sigma_hm_1 <= hm_tol)) &&
       ((sigma_x_2 <= x_tol) || (sigma_hm_2 <= hm_tol)))) {

    xphmp_integral_covar = ((((mixt_frac * precip_frac_1) * (mu_x_1 - x_mean)) * mu_hm_1) +
                            ((((one - mixt_frac) * precip_frac_2) * (mu_x_2 - x_mean)) * mu_hm_2));

  } else {
    if (((sigma_x_1 <= x_tol) || (sigma_hm_1 <= hm_tol))) {

      xphmp_integral_covar =
          ((((mixt_frac * precip_frac_1) * (mu_x_1 - x_mean)) * mu_hm_1) +
           (((one - mixt_frac) * precip_frac_2) *
            (((mu_x_2 - x_mean) * mu_hm_2) + ((corr_x_hm_2 * sigma_x_2) * sigma_hm_2))));

    } else {
      if (((sigma_x_2 <= x_tol) || (sigma_hm_2 <= hm_tol))) {

        xphmp_integral_covar =
            (((mixt_frac * precip_frac_1) *
              (((mu_x_1 - x_mean) * mu_hm_1) + ((corr_x_hm_1 * sigma_x_1) * sigma_hm_1))) +
             ((((one - mixt_frac) * precip_frac_2) * (mu_x_2 - x_mean)) * mu_hm_2));

      } else {
        // and sigma_x_2 > x_tol and sigma_hm_2 > hm_tol

        xphmp_integral_covar =
            (((mixt_frac * precip_frac_1) *
              (((mu_x_1 - x_mean) * mu_hm_1) + ((corr_x_hm_1 * sigma_x_1) * sigma_hm_1))) +
             (((one - mixt_frac) * precip_frac_2) *
              (((mu_x_2 - x_mean) * mu_hm_2) + ((corr_x_hm_2 * sigma_x_2) * sigma_hm_2))));
      }
    }
  }

  return xphmp_integral_covar;
}
KOKKOS_INLINE_FUNCTION real_t hmxphmyp_integral_covar(
    real_t mu_hmx_1, real_t mu_hmx_2, real_t mu_hmy_1, real_t mu_hmy_2, real_t sigma_hmx_1,
    real_t sigma_hmx_2, real_t sigma_hmy_1, real_t sigma_hmy_2, real_t corr_hmx_hmy_1,
    real_t corr_hmx_hmy_2, real_t mixt_frac, real_t precip_frac_1, real_t precip_frac_2,
    real_t hmx_mean, real_t hmy_mean, real_t hmx_tol, real_t hmy_tol) {
  // Description:
  // Solves for the covariance of two precipitating hydrometoers < hmx'hmy' >.
  // The variable "hmx" stands for any precipitating hydrometeor.  This applies
  // to such variables as rr, Nr, ri, Ni, etc.  The variable "hmy" stands for
  // any different precipitating hydrometeor.
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable

  real_t hmxphmyp_integral_covar; // Covariance of hmx and hmy     [hmx*hmy units]

  if ((((sigma_hmx_1 <= hmx_tol) || (sigma_hmy_1 <= hmy_tol)) &&
       ((sigma_hmx_2 <= hmx_tol) || (sigma_hmy_2 <= hmy_tol)))) {

    hmxphmyp_integral_covar = (((((mixt_frac * precip_frac_1) * mu_hmx_1) * mu_hmy_1) +
                                ((((one - mixt_frac) * precip_frac_2) * mu_hmx_2) * mu_hmy_2)) -
                               (hmx_mean * hmy_mean));

  } else {
    if (((sigma_hmx_1 <= hmx_tol) || (sigma_hmy_1 <= hmy_tol))) {

      hmxphmyp_integral_covar =
          (((((mixt_frac * precip_frac_1) * mu_hmx_1) * mu_hmy_1) +
            (((one - mixt_frac) * precip_frac_2) *
             ((mu_hmx_2 * mu_hmy_2) + ((corr_hmx_hmy_2 * sigma_hmx_2) * sigma_hmy_2)))) -
           (hmx_mean * hmy_mean));

    } else {
      if (((sigma_hmx_2 <= hmx_tol) || (sigma_hmy_2 <= hmy_tol))) {

        hmxphmyp_integral_covar =
            ((((mixt_frac * precip_frac_1) *
               ((mu_hmx_1 * mu_hmy_1) + ((corr_hmx_hmy_1 * sigma_hmx_1) * sigma_hmy_1))) +
              ((((one - mixt_frac) * precip_frac_2) * mu_hmx_2) * mu_hmy_2)) -
             (hmx_mean * hmy_mean));

      } else {
        // and sigma_hmx_2 > hmx_tol and sigma_hmy_2 > hmy_tol

        hmxphmyp_integral_covar =
            ((((mixt_frac * precip_frac_1) *
               ((mu_hmx_1 * mu_hmy_1) + ((corr_hmx_hmy_1 * sigma_hmx_1) * sigma_hmy_1))) +
              (((one - mixt_frac) * precip_frac_2) *
               ((mu_hmx_2 * mu_hmy_2) + ((corr_hmx_hmy_2 * sigma_hmx_2) * sigma_hmy_2)))) -
             (hmx_mean * hmy_mean));
      }
    }
  }

  return hmxphmyp_integral_covar;
}
KOKKOS_INLINE_FUNCTION real_t xp_a_hmpb_integrals_all_MM(
    real_t mu_x_1, real_t mu_x_2, real_t mu_hm_1, real_t mu_hm_2, real_t mu_hm_1_n,
    real_t mu_hm_2_n, real_t sigma_x_1, real_t sigma_x_2, real_t sigma_hm_1, real_t sigma_hm_2,
    real_t sigma_hm_1_n, real_t sigma_hm_2_n, real_t corr_x_hm_1_n, real_t corr_x_hm_2_n,
    real_t mixt_frac, real_t precip_frac_1, real_t precip_frac_2, real_t x_mean, real_t hm_mean,
    real_t x_tol, real_t hm_tol, int a_exp, int b_exp) {
  // Description:
  // Solves for any mixed moment < x'^a hm'^b >, where a and b are both
  // integers with values greater than or equal to 0.  The variable "x" is a
  // variable that is distributed binormally (meaning that it has a
  // normally-distributed individual marginal in each PDF component).  This
  // applies to such variables as w, rt, thl, sclr, etc.  The variable "hm"
  // stands for any precipitating hydrometeor.  This applies to such variables
  // as rr, Nr, ri, Ni, etc.
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable

  real_t xp_a_hmpb_integrals_all_MM;

  xp_a_hmpb_integrals_all_MM =
      ((mixt_frac * bivar_NL_x_hm_all_MM_comp_eq(mu_x_1, mu_hm_1, mu_hm_1_n, sigma_x_1, sigma_hm_1,
                                                 sigma_hm_1_n, corr_x_hm_1_n, precip_frac_1, x_mean,
                                                 hm_mean, x_tol, hm_tol, a_exp, b_exp)) +
       ((one - mixt_frac) * bivar_NL_x_hm_all_MM_comp_eq(mu_x_2, mu_hm_2, mu_hm_2_n, sigma_x_2,
                                                         sigma_hm_2, sigma_hm_2_n, corr_x_hm_2_n,
                                                         precip_frac_2, x_mean, hm_mean, x_tol,
                                                         hm_tol, a_exp, b_exp)));

  return xp_a_hmpb_integrals_all_MM;
}
KOKKOS_INLINE_FUNCTION real_t bivar_NL_x_hm_all_MM_comp_eq(
    real_t mu_x_i, real_t mu_hm_i, real_t mu_hm_i_n, real_t sigma_x_i, real_t sigma_hm_i,
    real_t sigma_hm_i_n, real_t corr_x_hm_i_n, real_t precip_frac_i, real_t x_mean, real_t hm_mean,
    real_t x_tol, real_t hm_tol, int a_exp, int b_exp) {
  // Description:
  // Solves the portion of the integral for < x'^a hm'^b > that relates to the
  // ith PDF component.  This takes into account the portions of the component
  // that are within-precipitation and outside-precipitation.  The equation
  // and function that are used depends on whether or not x and/or hm vary
  // inside of a PDF component.
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable

  real_t bivar_NL_x_hm_all_MM_comp_eq;

  if (((sigma_x_i <= x_tol) && (sigma_hm_i <= hm_tol))) {

    // Both x and hm do not vary in the ith PDF component.

    bivar_NL_x_hm_all_MM_comp_eq = (Kokkos::pow((mu_x_i - x_mean), a_exp) *
                                    ((precip_frac_i * Kokkos::pow((mu_hm_i - hm_mean), b_exp)) +
                                     ((one - precip_frac_i) * Kokkos::pow((-hm_mean), b_exp))));

  } else {
    if ((sigma_x_i <= x_tol)) {

      // Only x does not vary in the ith PDF component.

      bivar_NL_x_hm_all_MM_comp_eq =
          (Kokkos::pow((mu_x_i - x_mean), a_exp) *
           ((precip_frac_i *
             univar_L_int_PDF_comp_all_MM(mu_hm_i_n, sigma_hm_i_n, hm_mean, b_exp)) +
            ((one - precip_frac_i) * Kokkos::pow((-hm_mean), b_exp))));

    } else {
      if ((sigma_hm_i <= hm_tol)) {

        // Only hm does not vary in the ith PDF component.

        bivar_NL_x_hm_all_MM_comp_eq =
            (((precip_frac_i * Kokkos::pow((mu_hm_i - hm_mean), b_exp)) +
              ((one - precip_frac_i) * Kokkos::pow((-hm_mean), b_exp))) *
             univar_N_int_PDF_comp_all_MM(mu_x_i, sigma_x_i, x_mean, a_exp));

      } else {

        // Both x and hm vary in the ith PDF component.

        bivar_NL_x_hm_all_MM_comp_eq =
            ((precip_frac_i * bivar_NL_int_PDF_comp_all_MM(mu_x_i, mu_hm_i_n, sigma_x_i,
                                                           sigma_hm_i_n, corr_x_hm_i_n, x_mean,
                                                           hm_mean, a_exp, b_exp)) +
             (((one - precip_frac_i) * Kokkos::pow((-hm_mean), b_exp)) *
              univar_N_int_PDF_comp_all_MM(mu_x_i, sigma_x_i, x_mean, a_exp)));
      }
    }
  }

  return bivar_NL_x_hm_all_MM_comp_eq;
}
KOKKOS_INLINE_FUNCTION real_t bivar_NL_int_PDF_comp_all_MM(real_t mu_x1_i, real_t mu_x2_i_n,
                                                           real_t sigma_x1_i, real_t sigma_x2_i_n,
                                                           real_t corr_x1_x2_i_n, real_t x1_mean,
                                                           real_t x2_mean, int a_exp, int b_exp) {
  // Description:
  // This function is the evaluated form of the following integral:
  //
  // INT(-inf:inf) INT(0:inf) ( x1 - <x1> )^a ( x2 - <x2> )^b
  //                          * P_NL_i( x1, x2 ) dx2 dx1;
  //
  // where P_NL_i( x1, x2 ) is the functional form of a bivariate
  // normal-lognormal PDF (in the PDF component), x1 is a variable that has an
  // individual marginal that is distributed normally (in the PDF component),
  // and x2 is a variable that has an individual marginal that is distributed
  // lognormally (in the PDF component).  Additionally, <x1> is the overall
  // mean of x1, <x2> is the overall mean of x2, "a" is the integer (>= 0)
  // order of the mixed moment with respect to x1, and b is the integer (>= 0)
  // order of the mixed moment with respect to x2.
  //
  // When the integral is evaluated, the equation is:
  //
  // INT(-inf:inf) INT(0:inf) ( x1 - <x1> )^a ( x2 - <x2> )^b
  //                          * P_NL_i( x1, x2 ) dx2 dx1
  // = SUM( p = 0:floor(a/2) ) SUM ( q = 0:b )
  //   ( a! / ( ( a - 2p )! p! ) ) * ( b! / ( ( b - q )! q! ) )
  //   * [ (1/2) * sigma_x1_i**2 ]**p
  //   * ( mu_x1_i - <x1>
  //       + corr_x1_x2_i_n * sigma_x1_i * sigma_x2_i_n * q )**(a-2p)
  //   * ( -<x2> )**(b-q)
  //   * exp{ mu_x2_i_n * q + (1/2) * sigma_x2_i_n**2 * q**2 }.
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables

  real_t bivar_NL_int_PDF_comp_all_MM;
  int p;                 // Index relating to exponent a (a_exp)
  int q;                 // Index relating to exponent b (b_exp)
  real_t sigma_sum;      // Running sum total value
  real_t q_fp;           // Floating point version of index q
  real_t factorial_term; // Floating point version of factorial term

  // Initialize sigma_sum
  sigma_sum = zero;

  for (p = 0; p < (a_exp / 2) + 1; p += 1) {
    for (q = 0; q < b_exp + 1; q += 1) {

      q_fp = static_cast<real_t>(q);

      factorial_term =
          static_cast<real_t>(((factorial(a_exp) / (factorial(((a_exp - 2 * p))) * factorial(p))) *
                               (factorial(b_exp) / (factorial(((b_exp - q))) * factorial(q)))));

      sigma_sum = (sigma_sum +
                   ((((factorial_term * Kokkos::pow((one_half * (sigma_x1_i * sigma_x1_i)), p)) *
                      Kokkos::pow(((mu_x1_i - x1_mean) +
                                   (((corr_x1_x2_i_n * sigma_x1_i) * sigma_x2_i_n) * q_fp)),
                                  ((a_exp - 2 * p)))) *
                     Kokkos::pow((-x2_mean), ((b_exp - q)))) *
                    Kokkos::exp(((mu_x2_i_n * q_fp) +
                                 ((one_half * (sigma_x2_i_n * sigma_x2_i_n)) * (q_fp * q_fp))))));
    }
  }

  bivar_NL_int_PDF_comp_all_MM = sigma_sum;

  return bivar_NL_int_PDF_comp_all_MM;
}
KOKKOS_INLINE_FUNCTION real_t univar_N_int_PDF_comp_all_MM(real_t mu_x_i, real_t sigma_x_i,
                                                           real_t x_mean, int a_exp) {
  // Description:
  // This function is the evaluated form of the following integral:
  //
  // INT(-inf:inf) ( x - <x> )^a * P_N_i( x ) dx;
  //
  // where P_N_i( x ) is the functional form of a single-variable normal PDF
  // (in the PDF component), and x is a variable that has an individual
  // marginal that is distributed normally (in the PDF component).
  // Additionally, <x> is the overall mean of x, and "a" is the integer (>= 0)
  // order of the mixed moment with respect to x.
  //
  // When the integral is evaluated, the equation is:
  //
  // INT(-inf:inf) ( x - <x> )^a * P_N_i( x ) dx
  // = SUM( p = 0:floor(a/2) ) ( a! / ( ( a - 2p )! p! ) )
  //   * [ (1/2) * sigma_x_i**2 ]**p * ( mu_x_i - <x> )**(a-2p).
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables

  real_t univar_N_int_PDF_comp_all_MM;
  int p;                 // Index relating to exponent a (a_exp)
  real_t sigma_sum;      // Running sum total value
  real_t factorial_term; // Floating point version of factorial term

  // Initialize sigma_sum
  sigma_sum = zero;

  for (p = 0; p < (a_exp / 2) + 1; p += 1) {

    factorial_term =
        static_cast<real_t>((factorial(a_exp) / (factorial(((a_exp - 2 * p))) * factorial(p))));

    sigma_sum =
        (sigma_sum + ((factorial_term * Kokkos::pow((one_half * (sigma_x_i * sigma_x_i)), p)) *
                      Kokkos::pow((mu_x_i - x_mean), ((a_exp - 2 * p)))));
  }

  // Total integral within a component.
  univar_N_int_PDF_comp_all_MM = sigma_sum;

  return univar_N_int_PDF_comp_all_MM;
}
KOKKOS_INLINE_FUNCTION real_t univar_L_int_PDF_comp_all_MM(real_t mu_x_i_n, real_t sigma_x_i_n,
                                                           real_t x_mean, int b_exp) {
  // Description:
  // This function is the evaluated form of the following integral:
  //
  // INT(0:inf) ( x - <x> )^b * P_L_i( x ) dx;
  //
  // where P_L_i( x ) is the functional form of a single-variable lognormal PDF
  // (in the PDF component), and x is a variable that has an individual
  // marginal that is distributed lognormally (in the PDF component).
  // Additionally, <x> is the overall mean of x, and b is the integer (>= 0)
  // order of the mixed moment with respect to x.
  //
  // When the integral is evaluated, the equation is:
  //
  // INT(0:inf) ( x - <x> )^b * P_L_i( x ) dx;
  // = SUM ( q = 0:b ) ( b! / ( ( b - q )! q! ) )
  //   * ( -<x> )**(b-q)
  //   * exp{ mu_x_i_n * q + (1/2) * sigma_x_i_n**2 * q**2 }.
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables

  real_t univar_L_int_PDF_comp_all_MM;
  int q;                 // Index relating to exponent b (b_exp)
  real_t sigma_sum;      // Running sum total value
  real_t q_fp;           // Floating point version of index q
  real_t factorial_term; // Floating point version of factorial term

  // Initialize sigma_sum
  sigma_sum = zero;

  for (q = 0; q < b_exp + 1; q += 1) {

    q_fp = static_cast<real_t>(q);

    factorial_term =
        static_cast<real_t>((factorial(b_exp) / (factorial(((b_exp - q))) * factorial(q))));

    sigma_sum =
        (sigma_sum + ((factorial_term * Kokkos::pow((-x_mean), ((b_exp - q)))) *
                      Kokkos::exp(((mu_x_i_n * q_fp) +
                                   ((one_half * (sigma_x_i_n * sigma_x_i_n)) * (q_fp * q_fp))))));
  }

  univar_L_int_PDF_comp_all_MM = sigma_sum;

  return univar_L_int_PDF_comp_all_MM;
}
// Host/strided callers stage through the existing array owner; device Views alias.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9>
inline void
hydrometeor_mixed_moments(ScratchWorkspace& workspace, const grid& gr, int ngrdcol, int nzt,
                          int pdf_dim, int hydromet_dim, Array0 hydromet,
                          const hm_metadata_type& hm_metadata, Array1 mu_x_1_n, Array2 mu_x_2_n,
                          Array3 sigma_x_1_n, Array4 sigma_x_2_n, Array5 corr_array_1_n,
                          Array6 corr_array_2_n, const pdf_parameter& pdf_params,
                          const hydromet_pdf_parameter_collection& hydromet_pdf_params,
                          const precipitation_fractions& precip_fracs, stats_type& stats,
                          Array7 rtphmp_zt, Array8 thlphmp_zt, Array9 wp2hmp) {
  // Native argument order is unchanged; only memory-space staging occurs here.
  auto workspace_scope          = workspace.scope();
  auto hydromet_execution       = to_execution_space(workspace, hydromet, true);
  auto mu_x_1_n_execution       = to_execution_space(workspace, mu_x_1_n, true);
  auto mu_x_2_n_execution       = to_execution_space(workspace, mu_x_2_n, true);
  auto sigma_x_1_n_execution    = to_execution_space(workspace, sigma_x_1_n, true);
  auto sigma_x_2_n_execution    = to_execution_space(workspace, sigma_x_2_n, true);
  auto corr_array_1_n_execution = to_execution_space(workspace, corr_array_1_n, true);
  auto corr_array_2_n_execution = to_execution_space(workspace, corr_array_2_n, true);
  auto rtphmp_zt_execution      = to_execution_space(workspace, rtphmp_zt, false);
  auto thlphmp_zt_execution     = to_execution_space(workspace, thlphmp_zt, false);
  auto wp2hmp_execution         = to_execution_space(workspace, wp2hmp, false);
  hydrometeor_mixed_moments(workspace, gr, ngrdcol, nzt, pdf_dim, hydromet_dim, hydromet_execution,
                            hm_metadata, mu_x_1_n_execution, mu_x_2_n_execution,
                            sigma_x_1_n_execution, sigma_x_2_n_execution, corr_array_1_n_execution,
                            corr_array_2_n_execution, pdf_params, hydromet_pdf_params, precip_fracs,
                            stats, rtphmp_zt_execution, thlphmp_zt_execution, wp2hmp_execution);
  from_execution_space(workspace, rtphmp_zt, rtphmp_zt_execution);
  from_execution_space(workspace, thlphmp_zt, thlphmp_zt_execution);
  from_execution_space(workspace, wp2hmp, wp2hmp_execution);
}
} // namespace clubb_kokkos
