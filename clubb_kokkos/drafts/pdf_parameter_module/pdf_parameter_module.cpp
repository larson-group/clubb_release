// Allocate, reset, inspect and serialize the parameters of CLUBB's PDF.
//
// Public routines:
//   init_pdf_params: allocate and zero all 47 PDF parameter fields.
//   zero_pdf_params_api: reset PDF fields while retaining their allocations.
//   init_pdf_implicit_coefs_terms_api: allocate and zero closure coefficients and terms.
//   zero_pdf_implicit_coefs_terms_api: reset existing closure coefficients and terms.
//   print_pdf_params: print every parameter profile for the requested columns.
//   pack_pdf_params_api: serialize the first column into 47 ordered field slots.
//   unpack_pdf_params_api: restore those slots into the first column only.
//
// Existing generated types own (column,level[,scalar]) LayoutLeft Views.
// Initializers replace the caller's prior storage; resets preserve it. Packing
// uses a (level,field) buffer. A level range is honored only when both optional
// Native one-based bounds are supplied, and its length must equal nz. Unpacking
// leaves other columns and metadata unchanged. Scalar coefficients are allocated
// only when sclr_dim is positive. Diagnostic text spacing is not guaranteed to
// reproduce Native list-directed output.
#include "pdf_parameter_module.hpp"
#include "array_types.hpp"
#include "constants_clubb.hpp"
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace clubb_kokkos {

// =============================================================================
void init_pdf_params(
    int nz, int ngrdcol, // In
    pdf_parameter& pdf_params // Out
) {
  // Description:
  // Initializes all PDF parameters in the variable type pdf_parameter.
  // References:
  //--------------------------------------------------------------------
  // nz: Number of vertical grid levels [-].
  // ngrdcol: Number of grid columns [-].
  // pdf_params: PDF parameters [units vary].

  // Native IntentOut releases allocated components of this caller object.
  pdf_params = {};

  // Allocate pdf_params in Native field order.
  pdf_params.w_1 = decltype(pdf_params.w_1)("pdf_params::w_1", ngrdcol, nz);
  pdf_params.w_2 = decltype(pdf_params.w_2)("pdf_params::w_2", ngrdcol, nz);
  pdf_params.varnce_w_1 = decltype(pdf_params.varnce_w_1)("pdf_params::varnce_w_1", ngrdcol, nz);
  pdf_params.varnce_w_2 = decltype(pdf_params.varnce_w_2)("pdf_params::varnce_w_2", ngrdcol, nz);
  pdf_params.rt_1 = decltype(pdf_params.rt_1)("pdf_params::rt_1", ngrdcol, nz);
  pdf_params.rt_2 = decltype(pdf_params.rt_2)("pdf_params::rt_2", ngrdcol, nz);
  pdf_params.varnce_rt_1 = decltype(pdf_params.varnce_rt_1)("pdf_params::varnce_rt_1", ngrdcol, nz);
  pdf_params.varnce_rt_2 = decltype(pdf_params.varnce_rt_2)("pdf_params::varnce_rt_2", ngrdcol, nz);
  pdf_params.thl_1 = decltype(pdf_params.thl_1)("pdf_params::thl_1", ngrdcol, nz);
  pdf_params.thl_2 = decltype(pdf_params.thl_2)("pdf_params::thl_2", ngrdcol, nz);
  pdf_params.varnce_thl_1 =
      decltype(pdf_params.varnce_thl_1)("pdf_params::varnce_thl_1", ngrdcol, nz);
  pdf_params.varnce_thl_2 =
      decltype(pdf_params.varnce_thl_2)("pdf_params::varnce_thl_2", ngrdcol, nz);
  pdf_params.corr_w_rt_1 = decltype(pdf_params.corr_w_rt_1)("pdf_params::corr_w_rt_1", ngrdcol, nz);
  pdf_params.corr_w_rt_2 = decltype(pdf_params.corr_w_rt_2)("pdf_params::corr_w_rt_2", ngrdcol, nz);
  pdf_params.corr_w_thl_1 =
      decltype(pdf_params.corr_w_thl_1)("pdf_params::corr_w_thl_1", ngrdcol, nz);
  pdf_params.corr_w_thl_2 =
      decltype(pdf_params.corr_w_thl_2)("pdf_params::corr_w_thl_2", ngrdcol, nz);
  pdf_params.corr_rt_thl_1 =
      decltype(pdf_params.corr_rt_thl_1)("pdf_params::corr_rt_thl_1", ngrdcol, nz);
  pdf_params.corr_rt_thl_2 =
      decltype(pdf_params.corr_rt_thl_2)("pdf_params::corr_rt_thl_2", ngrdcol, nz);
  pdf_params.alpha_thl = decltype(pdf_params.alpha_thl)("pdf_params::alpha_thl", ngrdcol, nz);
  pdf_params.alpha_rt = decltype(pdf_params.alpha_rt)("pdf_params::alpha_rt", ngrdcol, nz);
  pdf_params.crt_1 = decltype(pdf_params.crt_1)("pdf_params::crt_1", ngrdcol, nz);
  pdf_params.crt_2 = decltype(pdf_params.crt_2)("pdf_params::crt_2", ngrdcol, nz);
  pdf_params.cthl_1 = decltype(pdf_params.cthl_1)("pdf_params::cthl_1", ngrdcol, nz);
  pdf_params.cthl_2 = decltype(pdf_params.cthl_2)("pdf_params::cthl_2", ngrdcol, nz);
  pdf_params.chi_1 = decltype(pdf_params.chi_1)("pdf_params::chi_1", ngrdcol, nz);
  pdf_params.chi_2 = decltype(pdf_params.chi_2)("pdf_params::chi_2", ngrdcol, nz);
  pdf_params.stdev_chi_1 = decltype(pdf_params.stdev_chi_1)("pdf_params::stdev_chi_1", ngrdcol, nz);
  pdf_params.stdev_chi_2 = decltype(pdf_params.stdev_chi_2)("pdf_params::stdev_chi_2", ngrdcol, nz);
  pdf_params.stdev_eta_1 = decltype(pdf_params.stdev_eta_1)("pdf_params::stdev_eta_1", ngrdcol, nz);
  pdf_params.stdev_eta_2 = decltype(pdf_params.stdev_eta_2)("pdf_params::stdev_eta_2", ngrdcol, nz);
  pdf_params.covar_chi_eta_1 =
      decltype(pdf_params.covar_chi_eta_1)("pdf_params::covar_chi_eta_1", ngrdcol, nz);
  pdf_params.covar_chi_eta_2 =
      decltype(pdf_params.covar_chi_eta_2)("pdf_params::covar_chi_eta_2", ngrdcol, nz);
  pdf_params.corr_w_chi_1 =
      decltype(pdf_params.corr_w_chi_1)("pdf_params::corr_w_chi_1", ngrdcol, nz);
  pdf_params.corr_w_chi_2 =
      decltype(pdf_params.corr_w_chi_2)("pdf_params::corr_w_chi_2", ngrdcol, nz);
  pdf_params.corr_w_eta_1 =
      decltype(pdf_params.corr_w_eta_1)("pdf_params::corr_w_eta_1", ngrdcol, nz);
  pdf_params.corr_w_eta_2 =
      decltype(pdf_params.corr_w_eta_2)("pdf_params::corr_w_eta_2", ngrdcol, nz);
  pdf_params.corr_chi_eta_1 =
      decltype(pdf_params.corr_chi_eta_1)("pdf_params::corr_chi_eta_1", ngrdcol, nz);
  pdf_params.corr_chi_eta_2 =
      decltype(pdf_params.corr_chi_eta_2)("pdf_params::corr_chi_eta_2", ngrdcol, nz);
  pdf_params.rsatl_1 = decltype(pdf_params.rsatl_1)("pdf_params::rsatl_1", ngrdcol, nz);
  pdf_params.rsatl_2 = decltype(pdf_params.rsatl_2)("pdf_params::rsatl_2", ngrdcol, nz);
  pdf_params.rc_1 = decltype(pdf_params.rc_1)("pdf_params::rc_1", ngrdcol, nz);
  pdf_params.rc_2 = decltype(pdf_params.rc_2)("pdf_params::rc_2", ngrdcol, nz);
  pdf_params.cloud_frac_1 =
      decltype(pdf_params.cloud_frac_1)("pdf_params::cloud_frac_1", ngrdcol, nz);
  pdf_params.cloud_frac_2 =
      decltype(pdf_params.cloud_frac_2)("pdf_params::cloud_frac_2", ngrdcol, nz);
  pdf_params.mixt_frac = decltype(pdf_params.mixt_frac)("pdf_params::mixt_frac", ngrdcol, nz);
  pdf_params.ice_supersat_frac_1 =
      decltype(pdf_params.ice_supersat_frac_1)("pdf_params::ice_supersat_frac_1", ngrdcol, nz);
  pdf_params.ice_supersat_frac_2 =
      decltype(pdf_params.ice_supersat_frac_2)("pdf_params::ice_supersat_frac_2", ngrdcol, nz);

  // Set metadata values.
  pdf_params.ngrdcol = ngrdcol;
  pdf_params.nz = nz;

  // Zero pdf_params.
  zero_pdf_params_api(
      pdf_params // InOut
  );
}

// =============================================================================
void zero_pdf_params_api(
    pdf_parameter& pdf_params // InOut
) {
  // Description:
  // Zeros all PDF parameters in the variable type pdf_parameter.
  // References:
  //--------------------------------------------------------------------
  // pdf_params: PDF parameters [units vary].

  // Capture component Views; keep the owning structure on the host.
  const auto pdf_params_w_1 = pdf_params.w_1;
  const auto pdf_params_w_2 = pdf_params.w_2;
  const auto pdf_params_varnce_w_1 = pdf_params.varnce_w_1;
  const auto pdf_params_varnce_w_2 = pdf_params.varnce_w_2;
  const auto pdf_params_rt_1 = pdf_params.rt_1;
  const auto pdf_params_rt_2 = pdf_params.rt_2;
  const auto pdf_params_varnce_rt_1 = pdf_params.varnce_rt_1;
  const auto pdf_params_varnce_rt_2 = pdf_params.varnce_rt_2;
  const auto pdf_params_thl_1 = pdf_params.thl_1;
  const auto pdf_params_thl_2 = pdf_params.thl_2;
  const auto pdf_params_varnce_thl_1 = pdf_params.varnce_thl_1;
  const auto pdf_params_varnce_thl_2 = pdf_params.varnce_thl_2;
  const auto pdf_params_corr_w_rt_1 = pdf_params.corr_w_rt_1;
  const auto pdf_params_corr_w_rt_2 = pdf_params.corr_w_rt_2;
  const auto pdf_params_corr_w_thl_1 = pdf_params.corr_w_thl_1;
  const auto pdf_params_corr_w_thl_2 = pdf_params.corr_w_thl_2;
  const auto pdf_params_corr_rt_thl_1 = pdf_params.corr_rt_thl_1;
  const auto pdf_params_corr_rt_thl_2 = pdf_params.corr_rt_thl_2;
  const auto pdf_params_alpha_thl = pdf_params.alpha_thl;
  const auto pdf_params_alpha_rt = pdf_params.alpha_rt;
  const auto pdf_params_crt_1 = pdf_params.crt_1;
  const auto pdf_params_crt_2 = pdf_params.crt_2;
  const auto pdf_params_cthl_1 = pdf_params.cthl_1;
  const auto pdf_params_cthl_2 = pdf_params.cthl_2;
  const auto pdf_params_chi_1 = pdf_params.chi_1;
  const auto pdf_params_chi_2 = pdf_params.chi_2;
  const auto pdf_params_stdev_chi_1 = pdf_params.stdev_chi_1;
  const auto pdf_params_stdev_chi_2 = pdf_params.stdev_chi_2;
  const auto pdf_params_stdev_eta_1 = pdf_params.stdev_eta_1;
  const auto pdf_params_stdev_eta_2 = pdf_params.stdev_eta_2;
  const auto pdf_params_covar_chi_eta_1 = pdf_params.covar_chi_eta_1;
  const auto pdf_params_covar_chi_eta_2 = pdf_params.covar_chi_eta_2;
  const auto pdf_params_corr_w_chi_1 = pdf_params.corr_w_chi_1;
  const auto pdf_params_corr_w_chi_2 = pdf_params.corr_w_chi_2;
  const auto pdf_params_corr_w_eta_1 = pdf_params.corr_w_eta_1;
  const auto pdf_params_corr_w_eta_2 = pdf_params.corr_w_eta_2;
  const auto pdf_params_corr_chi_eta_1 = pdf_params.corr_chi_eta_1;
  const auto pdf_params_corr_chi_eta_2 = pdf_params.corr_chi_eta_2;
  const auto pdf_params_rsatl_1 = pdf_params.rsatl_1;
  const auto pdf_params_rsatl_2 = pdf_params.rsatl_2;
  const auto pdf_params_rc_1 = pdf_params.rc_1;
  const auto pdf_params_rc_2 = pdf_params.rc_2;
  const auto pdf_params_cloud_frac_1 = pdf_params.cloud_frac_1;
  const auto pdf_params_cloud_frac_2 = pdf_params.cloud_frac_2;
  const auto pdf_params_mixt_frac = pdf_params.mixt_frac;
  const auto pdf_params_ice_supersat_frac_1 = pdf_params.ice_supersat_frac_1;
  const auto pdf_params_ice_supersat_frac_2 = pdf_params.ice_supersat_frac_2;
  const int ngrdcol = pdf_params.ngrdcol;
  const int nz = pdf_params.nz;

  // Initialize pdf_params in Native per-cell field order.
  Kokkos::parallel_for(
      "zero_pdf_params_api", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        pdf_params_w_1(i, k) = zero;
        pdf_params_w_2(i, k) = zero;
        pdf_params_varnce_w_1(i, k) = zero;
        pdf_params_varnce_w_2(i, k) = zero;
        pdf_params_rt_1(i, k) = zero;
        pdf_params_rt_2(i, k) = zero;
        pdf_params_varnce_rt_1(i, k) = zero;
        pdf_params_varnce_rt_2(i, k) = zero;
        pdf_params_thl_1(i, k) = zero;
        pdf_params_thl_2(i, k) = zero;
        pdf_params_varnce_thl_1(i, k) = zero;
        pdf_params_varnce_thl_2(i, k) = zero;
        pdf_params_corr_w_rt_1(i, k) = zero;
        pdf_params_corr_w_rt_2(i, k) = zero;
        pdf_params_corr_w_thl_1(i, k) = zero;
        pdf_params_corr_w_thl_2(i, k) = zero;
        pdf_params_corr_rt_thl_1(i, k) = zero;
        pdf_params_corr_rt_thl_2(i, k) = zero;
        pdf_params_alpha_thl(i, k) = zero;
        pdf_params_alpha_rt(i, k) = zero;
        pdf_params_crt_1(i, k) = zero;
        pdf_params_crt_2(i, k) = zero;
        pdf_params_cthl_1(i, k) = zero;
        pdf_params_cthl_2(i, k) = zero;
        pdf_params_chi_1(i, k) = zero;
        pdf_params_chi_2(i, k) = zero;
        pdf_params_stdev_chi_1(i, k) = zero;
        pdf_params_stdev_chi_2(i, k) = zero;
        pdf_params_stdev_eta_1(i, k) = zero;
        pdf_params_stdev_eta_2(i, k) = zero;
        pdf_params_covar_chi_eta_1(i, k) = zero;
        pdf_params_covar_chi_eta_2(i, k) = zero;
        pdf_params_corr_w_chi_1(i, k) = zero;
        pdf_params_corr_w_chi_2(i, k) = zero;
        pdf_params_corr_w_eta_1(i, k) = zero;
        pdf_params_corr_w_eta_2(i, k) = zero;
        pdf_params_corr_chi_eta_1(i, k) = zero;
        pdf_params_corr_chi_eta_2(i, k) = zero;
        pdf_params_rsatl_1(i, k) = zero;
        pdf_params_rsatl_2(i, k) = zero;
        pdf_params_rc_1(i, k) = zero;
        pdf_params_rc_2(i, k) = zero;
        pdf_params_cloud_frac_1(i, k) = zero;
        pdf_params_cloud_frac_2(i, k) = zero;
        pdf_params_mixt_frac(i, k) = zero;
        pdf_params_ice_supersat_frac_1(i, k) = zero;
        pdf_params_ice_supersat_frac_2(i, k) = zero;
      });
}

// =============================================================================
void init_pdf_implicit_coefs_terms_api(
    int nz, int ngrdcol, int sclr_dim, // In
    implicit_coefs_terms& pdf_implicit_coefs_terms // Out
) {
  // Description:
  // Initializes all PDF implicit coefficients and explicit terms in the
  // variable type implicit_coefs_terms.
  // References:
  //--------------------------------------------------------------------
  // nz: Number of vertical grid levels [-].
  // ngrdcol: Number of grid columns [-].
  // sclr_dim: Number of scalar variables [-].
  // pdf_implicit_coefs_terms: Implicit coefficients and explicit terms [units vary].

  // Native IntentOut releases allocated components of this caller object.
  pdf_implicit_coefs_terms = {};

  // Allocate pdf_implicit_coefs_terms in Native field order.
  pdf_implicit_coefs_terms.coef_wp4_implicit = decltype(pdf_implicit_coefs_terms.coef_wp4_implicit)(
      "pdf_implicit_coefs_terms::coef_wp4_implicit", ngrdcol, nz);
  pdf_implicit_coefs_terms.coef_wp2rtp_implicit =
      decltype(pdf_implicit_coefs_terms.coef_wp2rtp_implicit)(
          "pdf_implicit_coefs_terms::coef_wp2rtp_implicit", ngrdcol, nz);
  pdf_implicit_coefs_terms.term_wp2rtp_explicit =
      decltype(pdf_implicit_coefs_terms.term_wp2rtp_explicit)(
          "pdf_implicit_coefs_terms::term_wp2rtp_explicit", ngrdcol, nz);
  pdf_implicit_coefs_terms.coef_wp2thlp_implicit =
      decltype(pdf_implicit_coefs_terms.coef_wp2thlp_implicit)(
          "pdf_implicit_coefs_terms::coef_wp2thlp_implicit", ngrdcol, nz);
  pdf_implicit_coefs_terms.term_wp2thlp_explicit =
      decltype(pdf_implicit_coefs_terms.term_wp2thlp_explicit)(
          "pdf_implicit_coefs_terms::term_wp2thlp_explicit", ngrdcol, nz);
  pdf_implicit_coefs_terms.coef_wp2up_implicit =
      decltype(pdf_implicit_coefs_terms.coef_wp2up_implicit)(
          "pdf_implicit_coefs_terms::coef_wp2up_implicit", ngrdcol, nz);
  pdf_implicit_coefs_terms.term_wp2up_explicit =
      decltype(pdf_implicit_coefs_terms.term_wp2up_explicit)(
          "pdf_implicit_coefs_terms::term_wp2up_explicit", ngrdcol, nz);
  pdf_implicit_coefs_terms.coef_wp2vp_implicit =
      decltype(pdf_implicit_coefs_terms.coef_wp2vp_implicit)(
          "pdf_implicit_coefs_terms::coef_wp2vp_implicit", ngrdcol, nz);
  pdf_implicit_coefs_terms.term_wp2vp_explicit =
      decltype(pdf_implicit_coefs_terms.term_wp2vp_explicit)(
          "pdf_implicit_coefs_terms::term_wp2vp_explicit", ngrdcol, nz);
  pdf_implicit_coefs_terms.coef_wprtp2_implicit =
      decltype(pdf_implicit_coefs_terms.coef_wprtp2_implicit)(
          "pdf_implicit_coefs_terms::coef_wprtp2_implicit", ngrdcol, nz);
  pdf_implicit_coefs_terms.term_wprtp2_explicit =
      decltype(pdf_implicit_coefs_terms.term_wprtp2_explicit)(
          "pdf_implicit_coefs_terms::term_wprtp2_explicit", ngrdcol, nz);
  pdf_implicit_coefs_terms.coef_wpthlp2_implicit =
      decltype(pdf_implicit_coefs_terms.coef_wpthlp2_implicit)(
          "pdf_implicit_coefs_terms::coef_wpthlp2_implicit", ngrdcol, nz);
  pdf_implicit_coefs_terms.term_wpthlp2_explicit =
      decltype(pdf_implicit_coefs_terms.term_wpthlp2_explicit)(
          "pdf_implicit_coefs_terms::term_wpthlp2_explicit", ngrdcol, nz);
  pdf_implicit_coefs_terms.coef_wprtpthlp_implicit =
      decltype(pdf_implicit_coefs_terms.coef_wprtpthlp_implicit)(
          "pdf_implicit_coefs_terms::coef_wprtpthlp_implicit", ngrdcol, nz);
  pdf_implicit_coefs_terms.term_wprtpthlp_explicit =
      decltype(pdf_implicit_coefs_terms.term_wprtpthlp_explicit)(
          "pdf_implicit_coefs_terms::term_wprtpthlp_explicit", ngrdcol, nz);
  pdf_implicit_coefs_terms.coef_wpup2_implicit =
      decltype(pdf_implicit_coefs_terms.coef_wpup2_implicit)(
          "pdf_implicit_coefs_terms::coef_wpup2_implicit", ngrdcol, nz);
  pdf_implicit_coefs_terms.term_wpup2_explicit =
      decltype(pdf_implicit_coefs_terms.term_wpup2_explicit)(
          "pdf_implicit_coefs_terms::term_wpup2_explicit", ngrdcol, nz);
  pdf_implicit_coefs_terms.coef_wpvp2_implicit =
      decltype(pdf_implicit_coefs_terms.coef_wpvp2_implicit)(
          "pdf_implicit_coefs_terms::coef_wpvp2_implicit", ngrdcol, nz);
  pdf_implicit_coefs_terms.term_wpvp2_explicit =
      decltype(pdf_implicit_coefs_terms.term_wpvp2_explicit)(
          "pdf_implicit_coefs_terms::term_wpvp2_explicit", ngrdcol, nz);

  if (sclr_dim > 0) {
    pdf_implicit_coefs_terms.coef_wp2sclrp_implicit =
        decltype(pdf_implicit_coefs_terms.coef_wp2sclrp_implicit)(
            "pdf_implicit_coefs_terms::coef_wp2sclrp_implicit", ngrdcol, nz, sclr_dim);
    pdf_implicit_coefs_terms.term_wp2sclrp_explicit =
        decltype(pdf_implicit_coefs_terms.term_wp2sclrp_explicit)(
            "pdf_implicit_coefs_terms::term_wp2sclrp_explicit", ngrdcol, nz, sclr_dim);
    pdf_implicit_coefs_terms.coef_wpsclrp2_implicit =
        decltype(pdf_implicit_coefs_terms.coef_wpsclrp2_implicit)(
            "pdf_implicit_coefs_terms::coef_wpsclrp2_implicit", ngrdcol, nz, sclr_dim);
    pdf_implicit_coefs_terms.term_wpsclrp2_explicit =
        decltype(pdf_implicit_coefs_terms.term_wpsclrp2_explicit)(
            "pdf_implicit_coefs_terms::term_wpsclrp2_explicit", ngrdcol, nz, sclr_dim);
    pdf_implicit_coefs_terms.coef_wprtpsclrp_implicit =
        decltype(pdf_implicit_coefs_terms.coef_wprtpsclrp_implicit)(
            "pdf_implicit_coefs_terms::coef_wprtpsclrp_implicit", ngrdcol, nz, sclr_dim);
    pdf_implicit_coefs_terms.term_wprtpsclrp_explicit =
        decltype(pdf_implicit_coefs_terms.term_wprtpsclrp_explicit)(
            "pdf_implicit_coefs_terms::term_wprtpsclrp_explicit", ngrdcol, nz, sclr_dim);
    pdf_implicit_coefs_terms.coef_wpthlpsclrp_implicit =
        decltype(pdf_implicit_coefs_terms.coef_wpthlpsclrp_implicit)(
            "pdf_implicit_coefs_terms::coef_wpthlpsclrp_implicit", ngrdcol, nz, sclr_dim);
    pdf_implicit_coefs_terms.term_wpthlpsclrp_explicit =
        decltype(pdf_implicit_coefs_terms.term_wpthlpsclrp_explicit)(
            "pdf_implicit_coefs_terms::term_wpthlpsclrp_explicit", ngrdcol, nz, sclr_dim);
  }

  // Set metadata values.
  pdf_implicit_coefs_terms.ngrdcol = ngrdcol;
  pdf_implicit_coefs_terms.nz = nz;
  pdf_implicit_coefs_terms.sclr_dim = sclr_dim;

  // Zero pdf_implicit_coefs_terms.
  zero_pdf_implicit_coefs_terms_api(
      pdf_implicit_coefs_terms // InOut
  );
}

// =============================================================================
void zero_pdf_implicit_coefs_terms_api(
    implicit_coefs_terms& pdf_implicit_coefs_terms // InOut
) {
  // Description:
  // Zeros all PDF implicit coefficients and explicit terms in the
  // variable type implicit_coefs_terms.
  // References:
  //--------------------------------------------------------------------
  // pdf_implicit_coefs_terms: Implicit coefficients and explicit terms [units vary].

  // Capture component Views; keep the owning structure on the host.
  const auto pdf_implicit_coefs_terms_coef_wp4_implicit =
      pdf_implicit_coefs_terms.coef_wp4_implicit;
  const auto pdf_implicit_coefs_terms_coef_wp2rtp_implicit =
      pdf_implicit_coefs_terms.coef_wp2rtp_implicit;
  const auto pdf_implicit_coefs_terms_term_wp2rtp_explicit =
      pdf_implicit_coefs_terms.term_wp2rtp_explicit;
  const auto pdf_implicit_coefs_terms_coef_wp2thlp_implicit =
      pdf_implicit_coefs_terms.coef_wp2thlp_implicit;
  const auto pdf_implicit_coefs_terms_term_wp2thlp_explicit =
      pdf_implicit_coefs_terms.term_wp2thlp_explicit;
  const auto pdf_implicit_coefs_terms_coef_wp2up_implicit =
      pdf_implicit_coefs_terms.coef_wp2up_implicit;
  const auto pdf_implicit_coefs_terms_term_wp2up_explicit =
      pdf_implicit_coefs_terms.term_wp2up_explicit;
  const auto pdf_implicit_coefs_terms_coef_wp2vp_implicit =
      pdf_implicit_coefs_terms.coef_wp2vp_implicit;
  const auto pdf_implicit_coefs_terms_term_wp2vp_explicit =
      pdf_implicit_coefs_terms.term_wp2vp_explicit;
  const auto pdf_implicit_coefs_terms_coef_wprtp2_implicit =
      pdf_implicit_coefs_terms.coef_wprtp2_implicit;
  const auto pdf_implicit_coefs_terms_term_wprtp2_explicit =
      pdf_implicit_coefs_terms.term_wprtp2_explicit;
  const auto pdf_implicit_coefs_terms_coef_wpthlp2_implicit =
      pdf_implicit_coefs_terms.coef_wpthlp2_implicit;
  const auto pdf_implicit_coefs_terms_term_wpthlp2_explicit =
      pdf_implicit_coefs_terms.term_wpthlp2_explicit;
  const auto pdf_implicit_coefs_terms_coef_wprtpthlp_implicit =
      pdf_implicit_coefs_terms.coef_wprtpthlp_implicit;
  const auto pdf_implicit_coefs_terms_term_wprtpthlp_explicit =
      pdf_implicit_coefs_terms.term_wprtpthlp_explicit;
  const auto pdf_implicit_coefs_terms_coef_wpup2_implicit =
      pdf_implicit_coefs_terms.coef_wpup2_implicit;
  const auto pdf_implicit_coefs_terms_term_wpup2_explicit =
      pdf_implicit_coefs_terms.term_wpup2_explicit;
  const auto pdf_implicit_coefs_terms_coef_wpvp2_implicit =
      pdf_implicit_coefs_terms.coef_wpvp2_implicit;
  const auto pdf_implicit_coefs_terms_term_wpvp2_explicit =
      pdf_implicit_coefs_terms.term_wpvp2_explicit;
  const auto pdf_implicit_coefs_terms_coef_wp2sclrp_implicit =
      pdf_implicit_coefs_terms.coef_wp2sclrp_implicit;
  const auto pdf_implicit_coefs_terms_term_wp2sclrp_explicit =
      pdf_implicit_coefs_terms.term_wp2sclrp_explicit;
  const auto pdf_implicit_coefs_terms_coef_wpsclrp2_implicit =
      pdf_implicit_coefs_terms.coef_wpsclrp2_implicit;
  const auto pdf_implicit_coefs_terms_term_wpsclrp2_explicit =
      pdf_implicit_coefs_terms.term_wpsclrp2_explicit;
  const auto pdf_implicit_coefs_terms_coef_wprtpsclrp_implicit =
      pdf_implicit_coefs_terms.coef_wprtpsclrp_implicit;
  const auto pdf_implicit_coefs_terms_term_wprtpsclrp_explicit =
      pdf_implicit_coefs_terms.term_wprtpsclrp_explicit;
  const auto pdf_implicit_coefs_terms_coef_wpthlpsclrp_implicit =
      pdf_implicit_coefs_terms.coef_wpthlpsclrp_implicit;
  const auto pdf_implicit_coefs_terms_term_wpthlpsclrp_explicit =
      pdf_implicit_coefs_terms.term_wpthlpsclrp_explicit;
  const int ngrdcol = pdf_implicit_coefs_terms.ngrdcol;
  const int nz = pdf_implicit_coefs_terms.nz;

  // Initialize pdf_implicit_coefs_terms in Native per-cell field order.
  Kokkos::parallel_for(
      "zero_pdf_implicit_coefs_terms_api", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        pdf_implicit_coefs_terms_coef_wp4_implicit(i, k) = zero;
        pdf_implicit_coefs_terms_coef_wp2rtp_implicit(i, k) = zero;
        pdf_implicit_coefs_terms_term_wp2rtp_explicit(i, k) = zero;
        pdf_implicit_coefs_terms_coef_wp2thlp_implicit(i, k) = zero;
        pdf_implicit_coefs_terms_term_wp2thlp_explicit(i, k) = zero;
        pdf_implicit_coefs_terms_coef_wp2up_implicit(i, k) = zero;
        pdf_implicit_coefs_terms_term_wp2up_explicit(i, k) = zero;
        pdf_implicit_coefs_terms_coef_wp2vp_implicit(i, k) = zero;
        pdf_implicit_coefs_terms_term_wp2vp_explicit(i, k) = zero;
        pdf_implicit_coefs_terms_coef_wprtp2_implicit(i, k) = zero;
        pdf_implicit_coefs_terms_term_wprtp2_explicit(i, k) = zero;
        pdf_implicit_coefs_terms_coef_wpthlp2_implicit(i, k) = zero;
        pdf_implicit_coefs_terms_term_wpthlp2_explicit(i, k) = zero;
        pdf_implicit_coefs_terms_coef_wprtpthlp_implicit(i, k) = zero;
        pdf_implicit_coefs_terms_term_wprtpthlp_explicit(i, k) = zero;
        pdf_implicit_coefs_terms_coef_wpup2_implicit(i, k) = zero;
        pdf_implicit_coefs_terms_term_wpup2_explicit(i, k) = zero;
        pdf_implicit_coefs_terms_coef_wpvp2_implicit(i, k) = zero;
        pdf_implicit_coefs_terms_term_wpvp2_explicit(i, k) = zero;
      });

  if (pdf_implicit_coefs_terms.sclr_dim > 0) {
    const int sclr_dim = pdf_implicit_coefs_terms.sclr_dim;
    Kokkos::parallel_for(
        "zero_pdf_implicit_coefs_terms_api::scalars", ngrdcol * nz * sclr_dim,
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = (cell / ngrdcol) % nz;
          const int sclr = cell / (ngrdcol * nz);
          pdf_implicit_coefs_terms_coef_wp2sclrp_implicit(i, k, sclr) = zero;
          pdf_implicit_coefs_terms_term_wp2sclrp_explicit(i, k, sclr) = zero;
          pdf_implicit_coefs_terms_coef_wpsclrp2_implicit(i, k, sclr) = zero;
          pdf_implicit_coefs_terms_term_wpsclrp2_explicit(i, k, sclr) = zero;
          pdf_implicit_coefs_terms_coef_wprtpsclrp_implicit(i, k, sclr) = zero;
          pdf_implicit_coefs_terms_term_wprtpsclrp_explicit(i, k, sclr) = zero;
          pdf_implicit_coefs_terms_coef_wpthlpsclrp_implicit(i, k, sclr) = zero;
          pdf_implicit_coefs_terms_term_wpthlpsclrp_explicit(i, k, sclr) = zero;
        });
  }
}

// =============================================================================
void print_pdf_params(
    const pdf_parameter& pdf_params, int ngrdcol // In
) {
  // Description:
  // Prints every value found in type pdf_parameter
  //-------------------------------------------------------------------
  // pdf_params: PDF parameters [units vary].
  // ngrdcol: Number of grid columns [-].

  // Host I/O snapshots retain all Native labels and field order.
  const auto w_1_host = to_host_space(pdf_params.w_1);
  const auto w_2_host = to_host_space(pdf_params.w_2);
  const auto varnce_w_1_host = to_host_space(pdf_params.varnce_w_1);
  const auto varnce_w_2_host = to_host_space(pdf_params.varnce_w_2);
  const auto rt_1_host = to_host_space(pdf_params.rt_1);
  const auto rt_2_host = to_host_space(pdf_params.rt_2);
  const auto varnce_rt_1_host = to_host_space(pdf_params.varnce_rt_1);
  const auto varnce_rt_2_host = to_host_space(pdf_params.varnce_rt_2);
  const auto thl_1_host = to_host_space(pdf_params.thl_1);
  const auto thl_2_host = to_host_space(pdf_params.thl_2);
  const auto varnce_thl_1_host = to_host_space(pdf_params.varnce_thl_1);
  const auto varnce_thl_2_host = to_host_space(pdf_params.varnce_thl_2);
  const auto corr_w_rt_1_host = to_host_space(pdf_params.corr_w_rt_1);
  const auto corr_w_rt_2_host = to_host_space(pdf_params.corr_w_rt_2);
  const auto corr_w_thl_1_host = to_host_space(pdf_params.corr_w_thl_1);
  const auto corr_w_thl_2_host = to_host_space(pdf_params.corr_w_thl_2);
  const auto corr_rt_thl_1_host = to_host_space(pdf_params.corr_rt_thl_1);
  const auto corr_rt_thl_2_host = to_host_space(pdf_params.corr_rt_thl_2);
  const auto alpha_thl_host = to_host_space(pdf_params.alpha_thl);
  const auto alpha_rt_host = to_host_space(pdf_params.alpha_rt);
  const auto crt_1_host = to_host_space(pdf_params.crt_1);
  const auto crt_2_host = to_host_space(pdf_params.crt_2);
  const auto cthl_1_host = to_host_space(pdf_params.cthl_1);
  const auto cthl_2_host = to_host_space(pdf_params.cthl_2);
  const auto chi_1_host = to_host_space(pdf_params.chi_1);
  const auto chi_2_host = to_host_space(pdf_params.chi_2);
  const auto stdev_chi_1_host = to_host_space(pdf_params.stdev_chi_1);
  const auto stdev_chi_2_host = to_host_space(pdf_params.stdev_chi_2);
  const auto stdev_eta_1_host = to_host_space(pdf_params.stdev_eta_1);
  const auto stdev_eta_2_host = to_host_space(pdf_params.stdev_eta_2);
  const auto covar_chi_eta_1_host = to_host_space(pdf_params.covar_chi_eta_1);
  const auto covar_chi_eta_2_host = to_host_space(pdf_params.covar_chi_eta_2);
  const auto corr_w_chi_1_host = to_host_space(pdf_params.corr_w_chi_1);
  const auto corr_w_chi_2_host = to_host_space(pdf_params.corr_w_chi_2);
  const auto corr_w_eta_1_host = to_host_space(pdf_params.corr_w_eta_1);
  const auto corr_w_eta_2_host = to_host_space(pdf_params.corr_w_eta_2);
  const auto corr_chi_eta_1_host = to_host_space(pdf_params.corr_chi_eta_1);
  const auto corr_chi_eta_2_host = to_host_space(pdf_params.corr_chi_eta_2);
  const auto rsatl_1_host = to_host_space(pdf_params.rsatl_1);
  const auto rsatl_2_host = to_host_space(pdf_params.rsatl_2);
  const auto rc_1_host = to_host_space(pdf_params.rc_1);
  const auto rc_2_host = to_host_space(pdf_params.rc_2);
  const auto cloud_frac_1_host = to_host_space(pdf_params.cloud_frac_1);
  const auto cloud_frac_2_host = to_host_space(pdf_params.cloud_frac_2);
  const auto mixt_frac_host = to_host_space(pdf_params.mixt_frac);
  const auto ice_supersat_frac_1_host = to_host_space(pdf_params.ice_supersat_frac_1);
  const auto ice_supersat_frac_2_host = to_host_space(pdf_params.ice_supersat_frac_2);
  std::ostringstream output;
  output << std::setprecision(std::numeric_limits<real_t>::max_digits10);
  for (int i = 1; i <= ngrdcol; ++i) {
    output << "PDF parameters for column " << i << "\n";
    output << "---------------------------------------------------\n";
    output << "w_1(i,:) = ";
    for (std::size_t k = 0; k < w_1_host.extent(1); ++k)
      output << w_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "w_2(i,:) = ";
    for (std::size_t k = 0; k < w_2_host.extent(1); ++k)
      output << w_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "varnce_w_1(i,:) = ";
    for (std::size_t k = 0; k < varnce_w_1_host.extent(1); ++k)
      output << varnce_w_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "varnce_w_2(i,:) = ";
    for (std::size_t k = 0; k < varnce_w_2_host.extent(1); ++k)
      output << varnce_w_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "rt_1(i,:) = ";
    for (std::size_t k = 0; k < rt_1_host.extent(1); ++k)
      output << rt_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "rt_2(i,:) = ";
    for (std::size_t k = 0; k < rt_2_host.extent(1); ++k)
      output << rt_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "varnce_rt_1(i,:) = ";
    for (std::size_t k = 0; k < varnce_rt_1_host.extent(1); ++k)
      output << varnce_rt_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "varnce_rt_2(i,:) = ";
    for (std::size_t k = 0; k < varnce_rt_2_host.extent(1); ++k)
      output << varnce_rt_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "thl_1(i,:) = ";
    for (std::size_t k = 0; k < thl_1_host.extent(1); ++k)
      output << thl_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "thl_2(i,:) = ";
    for (std::size_t k = 0; k < thl_2_host.extent(1); ++k)
      output << thl_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "varnce_thl_1(i,:) = ";
    for (std::size_t k = 0; k < varnce_thl_1_host.extent(1); ++k)
      output << varnce_thl_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "varnce_thl_2(i,:) = ";
    for (std::size_t k = 0; k < varnce_thl_2_host.extent(1); ++k)
      output << varnce_thl_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "corr_w_rt_1(i,:) = ";
    for (std::size_t k = 0; k < corr_w_rt_1_host.extent(1); ++k)
      output << corr_w_rt_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "corr_w_rt_2(i,:) = ";
    for (std::size_t k = 0; k < corr_w_rt_2_host.extent(1); ++k)
      output << corr_w_rt_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "corr_w_thl_1(i,:) = ";
    for (std::size_t k = 0; k < corr_w_thl_1_host.extent(1); ++k)
      output << corr_w_thl_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "corr_w_thl_2(i,:) = ";
    for (std::size_t k = 0; k < corr_w_thl_2_host.extent(1); ++k)
      output << corr_w_thl_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "corr_rt_thl_1(i,:) = ";
    for (std::size_t k = 0; k < corr_rt_thl_1_host.extent(1); ++k)
      output << corr_rt_thl_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "corr_rt_thl_2(i,:) = ";
    for (std::size_t k = 0; k < corr_rt_thl_2_host.extent(1); ++k)
      output << corr_rt_thl_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "alpha_thl(i,:) = ";
    for (std::size_t k = 0; k < alpha_thl_host.extent(1); ++k)
      output << alpha_thl_host(i - 1, k) << ' ';
    output << "\n";
    output << "alpha_rt(i,:) = ";
    for (std::size_t k = 0; k < alpha_rt_host.extent(1); ++k)
      output << alpha_rt_host(i - 1, k) << ' ';
    output << "\n";
    output << "crt_1(i,:) = ";
    for (std::size_t k = 0; k < crt_1_host.extent(1); ++k)
      output << crt_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "crt_2(i,:) = ";
    for (std::size_t k = 0; k < crt_2_host.extent(1); ++k)
      output << crt_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "cthl_1(i,:) = ";
    for (std::size_t k = 0; k < cthl_1_host.extent(1); ++k)
      output << cthl_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "cthl_2(i,:) = ";
    for (std::size_t k = 0; k < cthl_2_host.extent(1); ++k)
      output << cthl_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "chi_1(i,:) = ";
    for (std::size_t k = 0; k < chi_1_host.extent(1); ++k)
      output << chi_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "chi_2(i,:) = ";
    for (std::size_t k = 0; k < chi_2_host.extent(1); ++k)
      output << chi_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "stdev_chi_1(i,:) = ";
    for (std::size_t k = 0; k < stdev_chi_1_host.extent(1); ++k)
      output << stdev_chi_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "stdev_chi_2(i,:) = ";
    for (std::size_t k = 0; k < stdev_chi_2_host.extent(1); ++k)
      output << stdev_chi_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "stdev_eta_1(i,:) = ";
    for (std::size_t k = 0; k < stdev_eta_1_host.extent(1); ++k)
      output << stdev_eta_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "stdev_eta_2(i,:) = ";
    for (std::size_t k = 0; k < stdev_eta_2_host.extent(1); ++k)
      output << stdev_eta_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "covar_chi_eta_1(i,:) = ";
    for (std::size_t k = 0; k < covar_chi_eta_1_host.extent(1); ++k)
      output << covar_chi_eta_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "covar_chi_eta_2(i,:) = ";
    for (std::size_t k = 0; k < covar_chi_eta_2_host.extent(1); ++k)
      output << covar_chi_eta_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "corr_w_chi_1(i,:) = ";
    for (std::size_t k = 0; k < corr_w_chi_1_host.extent(1); ++k)
      output << corr_w_chi_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "corr_w_chi_2(i,:) = ";
    for (std::size_t k = 0; k < corr_w_chi_2_host.extent(1); ++k)
      output << corr_w_chi_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "corr_w_eta_1(i,:) = ";
    for (std::size_t k = 0; k < corr_w_eta_1_host.extent(1); ++k)
      output << corr_w_eta_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "corr_w_eta_2(i,:) = ";
    for (std::size_t k = 0; k < corr_w_eta_2_host.extent(1); ++k)
      output << corr_w_eta_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "corr_chi_eta_1(i,:) = ";
    for (std::size_t k = 0; k < corr_chi_eta_1_host.extent(1); ++k)
      output << corr_chi_eta_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "corr_chi_eta_2(i,:) = ";
    for (std::size_t k = 0; k < corr_chi_eta_2_host.extent(1); ++k)
      output << corr_chi_eta_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "rsatl_1(i,:) = ";
    for (std::size_t k = 0; k < rsatl_1_host.extent(1); ++k)
      output << rsatl_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "rsatl_2(i,:) = ";
    for (std::size_t k = 0; k < rsatl_2_host.extent(1); ++k)
      output << rsatl_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "rc_1(i,:) = ";
    for (std::size_t k = 0; k < rc_1_host.extent(1); ++k)
      output << rc_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "rc_2(i,:) = ";
    for (std::size_t k = 0; k < rc_2_host.extent(1); ++k)
      output << rc_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "cloud_frac_1(i,:) = ";
    for (std::size_t k = 0; k < cloud_frac_1_host.extent(1); ++k)
      output << cloud_frac_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "cloud_frac_2(i,:) = ";
    for (std::size_t k = 0; k < cloud_frac_2_host.extent(1); ++k)
      output << cloud_frac_2_host(i - 1, k) << ' ';
    output << "\n";
    output << "mixt_frac(i,:) = ";
    for (std::size_t k = 0; k < mixt_frac_host.extent(1); ++k)
      output << mixt_frac_host(i - 1, k) << ' ';
    output << "\n";
    output << "ice_supersat_frac_1(i,:) = ";
    for (std::size_t k = 0; k < ice_supersat_frac_1_host.extent(1); ++k)
      output << ice_supersat_frac_1_host(i - 1, k) << ' ';
    output << "\n";
    output << "ice_supersat_frac_2(i,:) = ";
    for (std::size_t k = 0; k < ice_supersat_frac_2_host.extent(1); ++k)
      output << ice_supersat_frac_2_host(i - 1, k) << ' ';
    output << "\n";
  }
  std::cerr << output.str();
}

// =============================================================================
void pack_pdf_params_api(
    const pdf_parameter& pdf_params, int nz, // In
    Real2D r_param_array, // Out
    std::optional<int> k_start_in, std::optional<int> k_end_in // Optional In
) {
  // Description: Pack the first PDF column into its Native 47-slot field array.
  // pdf_params: PDF parameters [units vary].
  // nz: Number of vertical grid levels [-].
  // r_param_array: Packed first-column PDF values, shape (nz,num_pdf_params) [units vary].
  // k_start_in: Optional Native one-based first level; used only with k_end_in [-].
  // k_end_in: Optional Native one-based last level; used only with k_start_in [-].

  int k_start;
  int k_end;

  if (k_start_in.has_value() && k_end_in.has_value()) {
    k_start = *k_start_in;
    k_end = *k_end_in;
  } else {
    k_start = 1;
    k_end = nz;
  }

  // Native array assignment requires a slice conformable with nz. Reject an
  // invalid shape here rather than substituting a different number of levels.
  if (k_end - k_start + 1 != nz)
    throw std::invalid_argument("pack_pdf_params_api: Native slice must contain nz levels");
  const auto pdf_params_w_1 = pdf_params.w_1;
  const auto pdf_params_w_2 = pdf_params.w_2;
  const auto pdf_params_varnce_w_1 = pdf_params.varnce_w_1;
  const auto pdf_params_varnce_w_2 = pdf_params.varnce_w_2;
  const auto pdf_params_rt_1 = pdf_params.rt_1;
  const auto pdf_params_rt_2 = pdf_params.rt_2;
  const auto pdf_params_varnce_rt_1 = pdf_params.varnce_rt_1;
  const auto pdf_params_varnce_rt_2 = pdf_params.varnce_rt_2;
  const auto pdf_params_thl_1 = pdf_params.thl_1;
  const auto pdf_params_thl_2 = pdf_params.thl_2;
  const auto pdf_params_varnce_thl_1 = pdf_params.varnce_thl_1;
  const auto pdf_params_varnce_thl_2 = pdf_params.varnce_thl_2;
  const auto pdf_params_corr_w_rt_1 = pdf_params.corr_w_rt_1;
  const auto pdf_params_corr_w_rt_2 = pdf_params.corr_w_rt_2;
  const auto pdf_params_corr_w_thl_1 = pdf_params.corr_w_thl_1;
  const auto pdf_params_corr_w_thl_2 = pdf_params.corr_w_thl_2;
  const auto pdf_params_corr_rt_thl_1 = pdf_params.corr_rt_thl_1;
  const auto pdf_params_corr_rt_thl_2 = pdf_params.corr_rt_thl_2;
  const auto pdf_params_alpha_thl = pdf_params.alpha_thl;
  const auto pdf_params_alpha_rt = pdf_params.alpha_rt;
  const auto pdf_params_crt_1 = pdf_params.crt_1;
  const auto pdf_params_crt_2 = pdf_params.crt_2;
  const auto pdf_params_cthl_1 = pdf_params.cthl_1;
  const auto pdf_params_cthl_2 = pdf_params.cthl_2;
  const auto pdf_params_chi_1 = pdf_params.chi_1;
  const auto pdf_params_chi_2 = pdf_params.chi_2;
  const auto pdf_params_stdev_chi_1 = pdf_params.stdev_chi_1;
  const auto pdf_params_stdev_chi_2 = pdf_params.stdev_chi_2;
  const auto pdf_params_stdev_eta_1 = pdf_params.stdev_eta_1;
  const auto pdf_params_stdev_eta_2 = pdf_params.stdev_eta_2;
  const auto pdf_params_covar_chi_eta_1 = pdf_params.covar_chi_eta_1;
  const auto pdf_params_covar_chi_eta_2 = pdf_params.covar_chi_eta_2;
  const auto pdf_params_corr_w_chi_1 = pdf_params.corr_w_chi_1;
  const auto pdf_params_corr_w_chi_2 = pdf_params.corr_w_chi_2;
  const auto pdf_params_corr_w_eta_1 = pdf_params.corr_w_eta_1;
  const auto pdf_params_corr_w_eta_2 = pdf_params.corr_w_eta_2;
  const auto pdf_params_corr_chi_eta_1 = pdf_params.corr_chi_eta_1;
  const auto pdf_params_corr_chi_eta_2 = pdf_params.corr_chi_eta_2;
  const auto pdf_params_rsatl_1 = pdf_params.rsatl_1;
  const auto pdf_params_rsatl_2 = pdf_params.rsatl_2;
  const auto pdf_params_rc_1 = pdf_params.rc_1;
  const auto pdf_params_rc_2 = pdf_params.rc_2;
  const auto pdf_params_cloud_frac_1 = pdf_params.cloud_frac_1;
  const auto pdf_params_cloud_frac_2 = pdf_params.cloud_frac_2;
  const auto pdf_params_mixt_frac = pdf_params.mixt_frac;
  const auto pdf_params_ice_supersat_frac_1 = pdf_params.ice_supersat_frac_1;
  const auto pdf_params_ice_supersat_frac_2 = pdf_params.ice_supersat_frac_2;

  // Only Native first-column levels change; other columns/metadata remain.
  Kokkos::parallel_for(
      "pack_pdf_params_api", nz, KOKKOS_LAMBDA(const int level) {
        const int k = k_start - 1 + level;
        r_param_array(level, 0) = pdf_params_w_1(0, k);
        r_param_array(level, 1) = pdf_params_w_2(0, k);
        r_param_array(level, 2) = pdf_params_varnce_w_1(0, k);
        r_param_array(level, 3) = pdf_params_varnce_w_2(0, k);
        r_param_array(level, 4) = pdf_params_rt_1(0, k);
        r_param_array(level, 5) = pdf_params_rt_2(0, k);
        r_param_array(level, 6) = pdf_params_varnce_rt_1(0, k);
        r_param_array(level, 7) = pdf_params_varnce_rt_2(0, k);
        r_param_array(level, 8) = pdf_params_thl_1(0, k);
        r_param_array(level, 9) = pdf_params_thl_2(0, k);
        r_param_array(level, 10) = pdf_params_varnce_thl_1(0, k);
        r_param_array(level, 11) = pdf_params_varnce_thl_2(0, k);
        r_param_array(level, 12) = pdf_params_corr_w_rt_1(0, k);
        r_param_array(level, 13) = pdf_params_corr_w_rt_2(0, k);
        r_param_array(level, 14) = pdf_params_corr_w_thl_1(0, k);
        r_param_array(level, 15) = pdf_params_corr_w_thl_2(0, k);
        r_param_array(level, 16) = pdf_params_corr_rt_thl_1(0, k);
        r_param_array(level, 17) = pdf_params_corr_rt_thl_2(0, k);
        r_param_array(level, 18) = pdf_params_alpha_thl(0, k);
        r_param_array(level, 19) = pdf_params_alpha_rt(0, k);
        r_param_array(level, 20) = pdf_params_crt_1(0, k);
        r_param_array(level, 21) = pdf_params_crt_2(0, k);
        r_param_array(level, 22) = pdf_params_cthl_1(0, k);
        r_param_array(level, 23) = pdf_params_cthl_2(0, k);
        r_param_array(level, 24) = pdf_params_chi_1(0, k);
        r_param_array(level, 25) = pdf_params_chi_2(0, k);
        r_param_array(level, 26) = pdf_params_stdev_chi_1(0, k);
        r_param_array(level, 27) = pdf_params_stdev_chi_2(0, k);
        r_param_array(level, 28) = pdf_params_stdev_eta_1(0, k);
        r_param_array(level, 29) = pdf_params_stdev_eta_2(0, k);
        r_param_array(level, 30) = pdf_params_covar_chi_eta_1(0, k);
        r_param_array(level, 31) = pdf_params_covar_chi_eta_2(0, k);
        r_param_array(level, 32) = pdf_params_corr_w_chi_1(0, k);
        r_param_array(level, 33) = pdf_params_corr_w_chi_2(0, k);
        r_param_array(level, 34) = pdf_params_corr_w_eta_1(0, k);
        r_param_array(level, 35) = pdf_params_corr_w_eta_2(0, k);
        r_param_array(level, 36) = pdf_params_corr_chi_eta_1(0, k);
        r_param_array(level, 37) = pdf_params_corr_chi_eta_2(0, k);
        r_param_array(level, 38) = pdf_params_rsatl_1(0, k);
        r_param_array(level, 39) = pdf_params_rsatl_2(0, k);
        r_param_array(level, 40) = pdf_params_rc_1(0, k);
        r_param_array(level, 41) = pdf_params_rc_2(0, k);
        r_param_array(level, 42) = pdf_params_cloud_frac_1(0, k);
        r_param_array(level, 43) = pdf_params_cloud_frac_2(0, k);
        r_param_array(level, 44) = pdf_params_mixt_frac(0, k);
        r_param_array(level, 45) = pdf_params_ice_supersat_frac_1(0, k);
        r_param_array(level, 46) = pdf_params_ice_supersat_frac_2(0, k);
      });
}

// =============================================================================
void unpack_pdf_params_api(
    Real2D r_param_array, int nz, // In
    pdf_parameter& pdf_params, // InOut
    std::optional<int> k_start_in, std::optional<int> k_end_in // Optional In
) {
  // Description: Restore the selected first-column levels from the Native 47-slot field array.
  // r_param_array: Packed first-column PDF values, shape (nz,num_pdf_params) [units vary].
  // nz: Number of vertical grid levels [-].
  // pdf_params: PDF parameters [units vary].
  // k_start_in: Optional Native one-based first level; used only with k_end_in [-].
  // k_end_in: Optional Native one-based last level; used only with k_start_in [-].

  int k_start;
  int k_end;

  if (k_start_in.has_value() && k_end_in.has_value()) {
    k_start = *k_start_in;
    k_end = *k_end_in;
  } else {
    k_start = 1;
    k_end = nz;
  }

  // Native array assignment requires a slice conformable with nz. Reject an
  // invalid shape here rather than substituting a different number of levels.
  if (k_end - k_start + 1 != nz)
    throw std::invalid_argument("unpack_pdf_params_api: Native slice must contain nz levels");
  const auto pdf_params_w_1 = pdf_params.w_1;
  const auto pdf_params_w_2 = pdf_params.w_2;
  const auto pdf_params_varnce_w_1 = pdf_params.varnce_w_1;
  const auto pdf_params_varnce_w_2 = pdf_params.varnce_w_2;
  const auto pdf_params_rt_1 = pdf_params.rt_1;
  const auto pdf_params_rt_2 = pdf_params.rt_2;
  const auto pdf_params_varnce_rt_1 = pdf_params.varnce_rt_1;
  const auto pdf_params_varnce_rt_2 = pdf_params.varnce_rt_2;
  const auto pdf_params_thl_1 = pdf_params.thl_1;
  const auto pdf_params_thl_2 = pdf_params.thl_2;
  const auto pdf_params_varnce_thl_1 = pdf_params.varnce_thl_1;
  const auto pdf_params_varnce_thl_2 = pdf_params.varnce_thl_2;
  const auto pdf_params_corr_w_rt_1 = pdf_params.corr_w_rt_1;
  const auto pdf_params_corr_w_rt_2 = pdf_params.corr_w_rt_2;
  const auto pdf_params_corr_w_thl_1 = pdf_params.corr_w_thl_1;
  const auto pdf_params_corr_w_thl_2 = pdf_params.corr_w_thl_2;
  const auto pdf_params_corr_rt_thl_1 = pdf_params.corr_rt_thl_1;
  const auto pdf_params_corr_rt_thl_2 = pdf_params.corr_rt_thl_2;
  const auto pdf_params_alpha_thl = pdf_params.alpha_thl;
  const auto pdf_params_alpha_rt = pdf_params.alpha_rt;
  const auto pdf_params_crt_1 = pdf_params.crt_1;
  const auto pdf_params_crt_2 = pdf_params.crt_2;
  const auto pdf_params_cthl_1 = pdf_params.cthl_1;
  const auto pdf_params_cthl_2 = pdf_params.cthl_2;
  const auto pdf_params_chi_1 = pdf_params.chi_1;
  const auto pdf_params_chi_2 = pdf_params.chi_2;
  const auto pdf_params_stdev_chi_1 = pdf_params.stdev_chi_1;
  const auto pdf_params_stdev_chi_2 = pdf_params.stdev_chi_2;
  const auto pdf_params_stdev_eta_1 = pdf_params.stdev_eta_1;
  const auto pdf_params_stdev_eta_2 = pdf_params.stdev_eta_2;
  const auto pdf_params_covar_chi_eta_1 = pdf_params.covar_chi_eta_1;
  const auto pdf_params_covar_chi_eta_2 = pdf_params.covar_chi_eta_2;
  const auto pdf_params_corr_w_chi_1 = pdf_params.corr_w_chi_1;
  const auto pdf_params_corr_w_chi_2 = pdf_params.corr_w_chi_2;
  const auto pdf_params_corr_w_eta_1 = pdf_params.corr_w_eta_1;
  const auto pdf_params_corr_w_eta_2 = pdf_params.corr_w_eta_2;
  const auto pdf_params_corr_chi_eta_1 = pdf_params.corr_chi_eta_1;
  const auto pdf_params_corr_chi_eta_2 = pdf_params.corr_chi_eta_2;
  const auto pdf_params_rsatl_1 = pdf_params.rsatl_1;
  const auto pdf_params_rsatl_2 = pdf_params.rsatl_2;
  const auto pdf_params_rc_1 = pdf_params.rc_1;
  const auto pdf_params_rc_2 = pdf_params.rc_2;
  const auto pdf_params_cloud_frac_1 = pdf_params.cloud_frac_1;
  const auto pdf_params_cloud_frac_2 = pdf_params.cloud_frac_2;
  const auto pdf_params_mixt_frac = pdf_params.mixt_frac;
  const auto pdf_params_ice_supersat_frac_1 = pdf_params.ice_supersat_frac_1;
  const auto pdf_params_ice_supersat_frac_2 = pdf_params.ice_supersat_frac_2;

  // Only Native first-column levels change; other columns/metadata remain.
  Kokkos::parallel_for(
      "unpack_pdf_params_api", nz, KOKKOS_LAMBDA(const int level) {
        const int k = k_start - 1 + level;
        pdf_params_w_1(0, k) = r_param_array(level, 0);
        pdf_params_w_2(0, k) = r_param_array(level, 1);
        pdf_params_varnce_w_1(0, k) = r_param_array(level, 2);
        pdf_params_varnce_w_2(0, k) = r_param_array(level, 3);
        pdf_params_rt_1(0, k) = r_param_array(level, 4);
        pdf_params_rt_2(0, k) = r_param_array(level, 5);
        pdf_params_varnce_rt_1(0, k) = r_param_array(level, 6);
        pdf_params_varnce_rt_2(0, k) = r_param_array(level, 7);
        pdf_params_thl_1(0, k) = r_param_array(level, 8);
        pdf_params_thl_2(0, k) = r_param_array(level, 9);
        pdf_params_varnce_thl_1(0, k) = r_param_array(level, 10);
        pdf_params_varnce_thl_2(0, k) = r_param_array(level, 11);
        pdf_params_corr_w_rt_1(0, k) = r_param_array(level, 12);
        pdf_params_corr_w_rt_2(0, k) = r_param_array(level, 13);
        pdf_params_corr_w_thl_1(0, k) = r_param_array(level, 14);
        pdf_params_corr_w_thl_2(0, k) = r_param_array(level, 15);
        pdf_params_corr_rt_thl_1(0, k) = r_param_array(level, 16);
        pdf_params_corr_rt_thl_2(0, k) = r_param_array(level, 17);
        pdf_params_alpha_thl(0, k) = r_param_array(level, 18);
        pdf_params_alpha_rt(0, k) = r_param_array(level, 19);
        pdf_params_crt_1(0, k) = r_param_array(level, 20);
        pdf_params_crt_2(0, k) = r_param_array(level, 21);
        pdf_params_cthl_1(0, k) = r_param_array(level, 22);
        pdf_params_cthl_2(0, k) = r_param_array(level, 23);
        pdf_params_chi_1(0, k) = r_param_array(level, 24);
        pdf_params_chi_2(0, k) = r_param_array(level, 25);
        pdf_params_stdev_chi_1(0, k) = r_param_array(level, 26);
        pdf_params_stdev_chi_2(0, k) = r_param_array(level, 27);
        pdf_params_stdev_eta_1(0, k) = r_param_array(level, 28);
        pdf_params_stdev_eta_2(0, k) = r_param_array(level, 29);
        pdf_params_covar_chi_eta_1(0, k) = r_param_array(level, 30);
        pdf_params_covar_chi_eta_2(0, k) = r_param_array(level, 31);
        pdf_params_corr_w_chi_1(0, k) = r_param_array(level, 32);
        pdf_params_corr_w_chi_2(0, k) = r_param_array(level, 33);
        pdf_params_corr_w_eta_1(0, k) = r_param_array(level, 34);
        pdf_params_corr_w_eta_2(0, k) = r_param_array(level, 35);
        pdf_params_corr_chi_eta_1(0, k) = r_param_array(level, 36);
        pdf_params_corr_chi_eta_2(0, k) = r_param_array(level, 37);
        pdf_params_rsatl_1(0, k) = r_param_array(level, 38);
        pdf_params_rsatl_2(0, k) = r_param_array(level, 39);
        pdf_params_rc_1(0, k) = r_param_array(level, 40);
        pdf_params_rc_2(0, k) = r_param_array(level, 41);
        pdf_params_cloud_frac_1(0, k) = r_param_array(level, 42);
        pdf_params_cloud_frac_2(0, k) = r_param_array(level, 43);
        pdf_params_mixt_frac(0, k) = r_param_array(level, 44);
        pdf_params_ice_supersat_frac_1(0, k) = r_param_array(level, 45);
        pdf_params_ice_supersat_frac_2(0, k) = r_param_array(level, 46);
      });
}

} // namespace clubb_kokkos
